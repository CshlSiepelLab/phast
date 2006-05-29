/* $Id: pasteCds.c,v 1.3 2006-05-29 02:11:48 acs Exp $ */

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <misc.h>
#include <msa.h>
#include <maf.h>
#include <gff.h>
#include <sufficient_stats.h>
#include "pasteCds.help"

/* for use with --mask-frame-shifts; if gaps are farther apart than
   DIST we won't consider compensatory frame shifts */
#define DIST 30

/* for use in has_stops below */
inline int is_stop_codon(char *str) {
 return (strncmp(str, "TAA", 3) == 0 || strncmp(str, "TAG", 3) == 0 ||
         strncmp(str, "TGA", 3) == 0);
}

/* check alignment for in-frame stops */
int has_stops(MSA *msa) {
  int i, j;
  for (i = 0; i < msa->nseqs; i++)
    for (j = 0; j <= (msa->length - 3); j += 3)
      if (is_stop_codon(&msa->seqs[i][j]))
        return TRUE;
  return FALSE;
}

/* mask frame-shifted regions */
void do_mask_frame_shifts(MSA *msa, int *err) {

  int i, j, k, l;
  *err = FALSE;

  /* sequence by sequence, identify regions that have at least one
     non-multiple-of-three-length gap and are flanked by gapless
     regions, then mask them with 'N's */

  for (i = 0; i < msa->nseqs; i++) {
    int count = -1;
    int start[256], end[256], ngaps[256], has_fshift[256];
    for (j = 0; j < msa->length; j++) {
      if (msa->seqs[i][j] == GAP_CHAR) {
        for (k = j+1; k < msa->length && msa->seqs[i][k] == GAP_CHAR; k++);

        if (count < 0 || j > end[count] + DIST) { /* start new interval */
          count++;
          assert(count < 256);
          start[count] = j;
          ngaps[count] = 0;
          has_fshift[count] = FALSE;
        }

        /* update current interval */
        end[count] = k;
        ngaps[count] += k - j;
        if ((k - j) % 3 != 0) has_fshift[count] = TRUE;
        j = k;
      }
    }

    for (l = 0; l <= count; l++) {
      if (has_fshift[l]) {
        if (ngaps[count] % 3 != 0) {
          *err = TRUE;          /* whole alignment will be shifted out
                                   of frame; could repair by inserting
                                   or deleting columns but for now
                                   just punt */
          return;
        }
        for (j = start[count]; j < end[count]; j++)
          msa->seqs[i][j] = 'N';
      }
    }
  }
}

int main(int argc, char *argv[]) {
  FILE *MSAF, *GFFF, *OUTF;
  MSA *msa;
  GFF_Set *gff;
  char c;
  int opt_idx;
  int i, j;
  char outname[50];
  List *l = lst_new_ptr(1);

  /* options and defaults */
  char *outroot = NULL;
  int check_orfs = FALSE, mask_frame_shifts = FALSE;
  msa_format_type out_format = FASTA, in_format = SS;

  struct option long_opts[] = {
    {"in-format", 1, 0, 'i'},
    {"out-format", 1, 0, 'o'},
    {"check-orfs", 0, 0, 'c'},
    {"mask-frame-shifts", 0, 0, 'm'},
    {"out-root", 1, 0, 'r'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "o:i:r:cmh", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'i':
      in_format = msa_str_to_format(optarg);
      if (in_format == -1)
        die("ERROR: bad input format.  Try 'pasteCds -h'.\n");
      break;
    case 'o':
      out_format = msa_str_to_format(optarg);
      if (out_format == -1)
        die("ERROR: bad output format.  Try 'pasteCds -h'.\n");
      break;
    case 'c':
      check_orfs = TRUE;
      break;
    case 'm':
      mask_frame_shifts = TRUE;
      break;
    case 'r':
      outroot = optarg;
      break;
    case 'h':
      printf(HELP);
      exit(0);
    case '?':
      die("Bad argument.  Try '%s -h'.\n", argv[0]);
    }
  }

  if (optind != argc - 2) 
    die("ERROR: two arguments expected.  Try '%s -h'.\n", argv[0]);
    
  MSAF = fopen_fname(argv[optind], "r");
  GFFF = fopen_fname(argv[optind+1], "r");

  fprintf(stderr, "Reading alignment...\n");
  if (in_format == MAF) 
    msa = maf_read(MSAF, NULL, 1, NULL, NULL, NULL, -1, TRUE, NULL, 
                   NO_STRIP, FALSE);
  else 
    msa = msa_new_from_file(MSAF, in_format, NULL);

  fprintf(stderr, "Reading GFF...\n");
  gff = gff_read_set(GFFF);

  /* make sure stops outside of CDS */
  gff_group(gff, "transcript_id");
  gff_create_signals(gff);
  gff_fix_start_stop(gff);
  gff_ungroup(gff);             /* will redo after coord mapping, filtering */

  fprintf(stderr, "Mapping coordinates...\n");
  msa_map_gff_coords(msa, gff, 1, 0, 0, NULL);

  lst_push_ptr(l, str_new_charstr("CDS"));
  gff_filter_by_type(gff, l, FALSE, NULL);
  gff_group(gff, "transcript_id");

  /* iterate through the gff groups */ 
  for (i = 0; i < lst_size(gff->groups); i++) {
    MSA *gene_msa = NULL;
    GFF_FeatureGroup *group = lst_get_ptr(gff->groups, i);
    char strand = '.';
    fprintf(stderr, "Processing gene '%s'...\n", group->name->chars);
    for (j = 0; j < lst_size(group->features); j++) {
      GFF_Feature *f = lst_get_ptr(group->features, j); 
      MSA *subaln = msa_sub_alignment(msa, NULL, FALSE, f->start-1, f->end);
      ss_to_msa(subaln);
      ss_free(subaln->ss);
      subaln->ss = NULL;
     
      if (j == 0)               /* initialize gene_msa */
        gene_msa = subaln;
      else {
        msa_concatenate(gene_msa, subaln);
        msa_free(subaln); 
      }
      if (strand == '.') strand = f->strand;
      else if (strand != f->strand) 
        die("ERROR: features in group don't have same strand.\n");
    }

    if (strand == '-')
      msa_reverse_compl(gene_msa);

    /* this is a hack to address the problem of stop codons that span
       splice sites, which are not dealt with properly by
       gff_fix_start_stop */
    if (gene_msa->length > 3 && 
        is_stop_codon(&gene_msa->seqs[0][gene_msa->length-3]))
          gene_msa->length -= 3;

    if (mask_frame_shifts) {
      int err;
      do_mask_frame_shifts(gene_msa, &err);
      if (err) {
        fprintf(stderr, "WARNING: gene '%s' has frame-shift that can't be repaired by masking; skipping.\n", 
                group->name->chars);
        msa_free(gene_msa); 
        continue; 
      }
    }

    if (check_orfs && has_stops(gene_msa)) {
      fprintf(stderr, "WARNING: gene '%s' has stop codons; skipping.\n", 
              group->name->chars);
      msa_free(gene_msa); 
      continue; 
    }
    
    if (outroot != NULL)
      sprintf(outname, "%s.%s.%s", outroot, group->name->chars,
              msa_suffix_for_format(out_format));
    else 
      sprintf(outname, "%s.%s", group->name->chars,
              msa_suffix_for_format(out_format));

    OUTF = fopen_fname(outname, "w+");
    msa_print(OUTF, gene_msa, out_format, FALSE);
    msa_free(gene_msa);
  }

  fprintf(stderr, "Done.\n");
  return 0;
}

