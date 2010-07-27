#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <getopt.h>
#include <misc.h>
#include <gff.h>
#include <lists.h>
#include <stringsplus.h>
#include <dmotif_phmm.h>
#include "dmcompare.help"

#define DEFAULT_OFFSET 0

void print_gff_sets(GFF_Set *matches, GFF_Set *mismatch,
		    GFF_Set *unique_to_query, GFF_Set *unique_to_target,
		    String *fname_root);

int main(int argc, char *argv[]) {
  char c;
  int i, opt_idx, len1, len2, nwins, *stats;
  double fnr, fpr, sens, spec;
  GFF_Set *target_gff, *query_gff, *matches, *mismatch, *unique_to_query,
    *unique_to_target;
  String *fname_root;
  List *tmpl;
  TreeModel *tm;
  FILE *f;

#ifdef RPHAST
  GetRNGstate(); //seed R's random number generator
#endif

  matches = mismatch = unique_to_query = unique_to_target = NULL;

  struct option long_opts[] = {
    {"gff-out", 0, 0, 'g'},
    {"compute-fpr", 1, 0, 'f'},
    {"offset", 1, 0, 'o'},
    {"inexact", 1, 0, 'i'},
    {"present-in-posteriors", 0, 0, 'p'},
    {"thresh", 1, 0, 't'},
    {"help", 0, 0, 'h'},
    {0,0,0,0}
  };

  /* arguments and defaults for options */
  int gff_out = FALSE, seqlen = 0, w = 0, offset = DEFAULT_OFFSET,
    present_in_posteriors = FALSE;
  double thresh = 0;
  tm = NULL;

  while ((c = getopt_long(argc, argv, "g, f, h, o", long_opts, &opt_idx)) 
	 != -1) {
    switch (c) {
    case 'g':
      gff_out = TRUE;
      break;
    case 'f':
      tmpl = get_arg_list_int(optarg);
      seqlen = lst_get_int(tmpl, 0);
      w = lst_get_int(tmpl, 1);
      if (seqlen <= 0 || w <= 0)
	die("Error: bad argument to --compute-fpr.\n");
      lst_free(tmpl);
      break;
    case 'o':
      offset = atoi(optarg);
      break;
    case 'i':
      f = fopen_fname(optarg, "r");
      tm = tm_new_from_file(f);
      tr_name_ancestors(tm->tree);
      fclose(f);
      break;
    case 'p':
      present_in_posteriors = TRUE;
      break;
    case 't':
      thresh = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'h':
      printf(HELP);
      exit(0);
    case '?':
      die("Bad argument. Try 'dmcompare -h'.\n");
    }
  }
  
  if (optind != argc - 2)
    die("Two arguments required. Try 'dmcompare -h'.\n");
  if (thresh > 0 && present_in_posteriors == FALSE)
    die("ERROR: --thresh requires --present-in-posteriors");
  if (present_in_posteriors && tm != NULL)
    die("ERROR: --inexact and --present-in-posteriors are not compatible!\n");

  fprintf(stderr, "Reading target gff from %s...\n", argv[optind]);
  target_gff = gff_read_set(fopen_fname(argv[optind], "r"));
  
  fprintf(stderr, "Reading query gff from %s...\n", argv[optind+1]);
  query_gff = gff_read_set(fopen_fname(argv[optind+1], "r"));

  if (gff_out == TRUE) {
    matches = gff_new_set();
    mismatch = gff_new_set();
    unique_to_query = gff_new_set();
    unique_to_target = gff_new_set();
  }

  len1 = lst_size(target_gff->features);
  len2 = lst_size(query_gff->features);
  stats = smalloc(4 * sizeof(int));

  for (i = 0; i < 4; i++)
    stats[i] = 0;

  dms_compare_gffs(target_gff, query_gff, stats, offset, matches, mismatch,
		   unique_to_query, unique_to_target, tm,
		   present_in_posteriors, thresh);

  /* Sensitivity is (# true pos) / (# true pos + # false neg) -- computed here
     as (# true pos) / (# features in target set) */
  sens = (double)stats[0] / (double)len1;
  /* False negative rate is (# false neg) / (# features in target set) */
  fnr = (double)stats[2] / (double)len1;
  
  printf("Features in target set: %d\n", len1);
  printf("Features in query set: %d\n", len2);
  printf("True positives: %d of %d  (Sensitivity = %f)\n", 
	 stats[0], len1, sens);
  printf("PPV: %f\n", ((double)stats[0]/(double)len2));
  printf("False negatives: %d (FNR = %f)\n", stats[2], fnr);
  
  if (seqlen > 0) {
    nwins = seqlen - (w - 1);
    /* False positive rate is (# false pos) / (total # neg) */
    fpr = (double)(stats[3] + stats[1]) / (double)(nwins - len1);
    /* Specificity is (# true neg) / (# true neg + # false pos) */
    spec = (double)(nwins - len2) / (double)((nwins - len2) + stats[3] 
					     + stats[1]);
    printf("False positives: %d (FPR = %f) \n\t%d misidentified features\n\t%d in locations with no annotated feature)\n"
	   , (stats[3] + stats[1]), fpr, stats[1], stats[3]);
    printf("Specificity = %f\n", spec);
  } else {
    printf("Misidentified features: %d\n", stats[1]);
  }

  /* Print features to gff files, if specified */
  if (gff_out == TRUE) {
    fname_root = str_new(STR_SHORT_LEN);
    fname_root = str_new_charstr(argv[optind+1]);
    str_root(fname_root, '.');
    
    print_gff_sets(matches, mismatch, unique_to_query, unique_to_target, 
		   fname_root);
  }
  return(0);
}

/* Print GFF_Sets to four separate outfiles, based on the query filename */
void print_gff_sets(GFF_Set *matches, GFF_Set *mismatch,
		    GFF_Set *unique_to_query, GFF_Set *unique_to_target,
		    String *fname_root) {
  String *fname_out = str_new(STR_SHORT_LEN);
  FILE *out;
    
  str_append(fname_out, fname_root);
  str_append_charstr(fname_out, ".match.gff");
  out = fopen_fname(fname_out->chars, "w");
  gff_print_set(out, matches);
  
  str_clear(fname_out);
  str_append(fname_out, fname_root);
  str_append_charstr(fname_out, ".mismatch.gff");
  out = fopen_fname(fname_out->chars, "w");
  gff_print_set(out, mismatch);
  
  str_clear(fname_out);
  str_append(fname_out, fname_root);
  str_append_charstr(fname_out, ".unique_to_query.gff");
  out = fopen_fname(fname_out->chars, "w");
  gff_print_set(out, unique_to_query);
  
  str_clear(fname_out);
  str_append(fname_out, fname_root);
  str_append_charstr(fname_out, ".unique_to_target.gff");
  out = fopen_fname(fname_out->chars, "w");
  gff_print_set(out, unique_to_target);  
}
