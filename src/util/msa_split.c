/* $Id: msa_split.c,v 1.2 2004-06-14 03:06:21 acs Exp $
   Written by Adam Siepel, 2002
   Copyright 2002, Adam Siepel, University of California */

/* to do: possibly combine this program with 'msa_view' */

#include <stdlib.h>
#include <stdio.h>
#include "sufficient_stats.h"
#include "msa.h"
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "gff.h"
#include "maf.h"

#define DOWNSTREAM_OTHER "other"
#define NSITES_BETWEEN_BLOCKS 30

void print_usage() {
    printf("\n\
USAGE: msa_split [OPTIONS] <msa_fname> \n\
\n\
Partitions a multiple sequence alignment either at designated columns,\n\
or according to specified category labels, and outputs sub-alignments\n\
for the partitions.  Optionally splits an associated annotations file.\n\
\n\
Options:\n\
\n\
    -i PHYLIP|FASTA|PSU|SS|MAF\n\
        Input alignment file format.  Default is PHYLIP.\n\
\n\
    -M <rseq_fname>\n\
        (for use with -i MAF) Name of file containing reference sequence, \n\
        in FASTA format.\n\
\n\
    -g <gff_fname>\n\
        Name of GFF file.  Frame of reference of feature indices is\n\
        determined feature-by-feature according to 'seqname'\n\
        attribute.  WARNING: features that span partition indices will\n\
        be discarded.\n\
\n\
    -c <feature_frame>  \n\
        Index of frame of reference for feature coordinates, as\n\
        defined in the GFF file.  Use an integer 1-N (if N seqs) or 0\n\
        to indicate the coordinate system of the alignment as a whole.\n\
        Default behavior is to match features with alignment sequences\n\
        by name (feature by feature).  This option is only for use\n\
        with -g.  Features in output files will be in the same frame\n\
        of reference as corresponding features in the input file.\n\
\n\
    -p <partition_indices>\n\
        List of explicit indices at which to split alignment\n\
        (comma-separated).  If the list of indices is \"10,20\",\n\
        then sub-alignments will be output for columns 1-9, 10-19, and\n\
        20-<msa_len>.  Default is the empty list (single partition).\n\
\n\
    -P\n\
        Split by \"groups\" in GFF file (the last column, or\n\
        \"attribute\" column, is taken to indicate group).  Every time\n\
        two consecutive features belong to different groups (no exact\n\
        string match), split midway between the end of the first and\n\
        the start of the second.  The order of the features is\n\
        important (no sorting occurs).  For use with -g.\n\
\n\
    -d <partition_frame>\n\
        Index of frame of reference for split indices.  Default is\n\
        0 (entire MSA).  For use with -p or -w.\n\
\n\
    -n <npartitions>\n\
        Split the alignment equally into specified number of partitions.\n\
\n\
    -w <win_size,win_overlap>\n\
        Split the alignment into \"windows\" of size <win_size> bases,\n\
        overlapping by <win_overlap>\n\
\n\
    -L <cat_map_fname>\n\
        Split according to category labels, as defined by category map\n\
        and GFF file.  For use with -g.\n\
\n\
    -C <cat_list>\n\
        (For use with -g and -L) Output sub-alignments for only the\n\
        specified categories (column-delimited list).\n\
\n\
    -s\n\
        Strand-sensitive mode.  Reverse complement all segments having\n\
        at least one feature on the reverse strand and none on the\n\
        positive strand.  For use with -g.\n\
\n\
    -f\n\
        Retain the coordinates of the original GFF file, rather than\n\
        resetting them for each partition (such that the first column\n\
        of each new alignment has coordinate 1).  For use with -g.\n\
\n\
    -G ALL|ANY|<seqno>\n\
        Strip columns in output alignments containing all gaps, any\n\
        gaps, or gaps in the specified sequence (<seqno>; indexing\n\
        begins with one).  Default is not to strip any columns.\n\
\n\
    -r <output_fname_root>\n\
        Root of filename for output files.  Suffixes will be .i.msa and\n\
        .i.gff, for i between 1 and the total number of partitions.\n\
        Default filename root is \"msa_split\".\n\
\n\
    -o PHYLIP|FASTA|PSU|SS\n\
        Output alignment file format.  Default is PHYLIP.\n\
\n\
    -O <name_list>\n\
        Change order of rows in alignment to match sequence names\n\
        specified in name_list.  If a name appears in name_list but\n\
        not in the alignment, a row of gaps will be inserted.  Order\n\
        is changed *before* -c and -d are applied.\n\
\n\
    -S \n\
        Output summary of each partition to \"<output_fname_root>.sum\"\n\
        (includes base frequencies and numbers of gapped columns).\n\
\n\
    -T <tuple_size>\n\
        (for use with -L and -g or -o SS)  With -L and -g, insert \n\
        tuple_size-1 columns of missing data ('N' characters) between \n\
        sites that were not adjacent in the original alignment, to avoid\n\
        the creation of artificial context.  With -o SS, express sufficient\n\
        statistics in terms of specified tuple size.\n\
\n\
    -z  \n\
        (For use with -o SS)  Suppress the portion of the sufficient \n\
        statistics concerned with the order in which columns appear.  Useful \n\
        for analyses for which order is unimportant.\n\
\n\
    -I <ninf_sites> \n\
        Only output partitions having at least <ninf_sites> informative sites\n\
        (sites at which at least two non-gap and non-N gaps are present).\n\
\n\
    -B <radius>\n\
        (for use with -p, -P, -n, or -w) Try to partition at sites\n\
        between alignment blocks.  Assumes a reference sequence\n\
        alignment, with the first sequence as the reference (as\n\
        created by multiz).  Blocks of %d sites with gaps in all\n\
        sequences but the reference are assumed to indicate boundaries\n\
        between alignment blocks.  Partition indices will not be moved \n\
        more than <radius> sites.\n\
\n\
    -q\n\
        Proceed quietly.\n\
\n\
    -h\n\
        Print this help message.\n\n", NSITES_BETWEEN_BLOCKS);


    /* defunct option */
/*     -D <old_cm_fname>,<new_cm_fname>\n\ */
/*         Distinguish between upstream and downstream \"other\"\n\ */
/*         (non-annotated) regions in each partition.  Adds a new\n\ */
/*         category to the category map, called \"other2\", and adds a\n\ */
/*         feature of this type at the end of each partition (extends\n\ */
/*         from last feature to end of partition).  Useful for training\n\ */
/*         single-gene HMMs.  Revises the category map in <old_cm_fname>\n\ */
/*         and writes the new version to <new_cm_fname>.  For use with -g\n\ */
/*         and -p or -P.\n\ */

}

void write_sub_msa(MSA *submsa, char *fname_root, int idx, 
                   msa_format_type output_format, int tuple_size, 
                   int ordered_stats) {
  String *outfname = str_new(STR_MED_LEN);
  FILE *F;

  str_cpy_charstr(outfname, fname_root);
  if (idx >= 0) {
    str_append_charstr(outfname, ".");
    str_append_int(outfname, idx);
  }
  str_append_charstr(outfname, output_format == SS ? ".ss" : ".msa");
  if ((F = fopen(outfname->chars, "w+")) == NULL) {
    fprintf(stderr, "ERROR: unable to open %s for writing.\n", 
            outfname->chars);
    exit(1);
  }

  /* create sufficient stats, if necessary */
  if (output_format == SS) {
    if (submsa->ss == NULL)
      ss_from_msas(submsa, tuple_size, ordered_stats, NULL, NULL, NULL, -1);
    else {
      if (submsa->ss->tuple_size != tuple_size) {
        fprintf(stderr, "ERROR: tuple size in SS file does not match desired tuple size for output.\nConversion not supported.\n");
        exit(1);
      }
      if (submsa->ss->tuple_idx != NULL && ordered_stats == 0) {
        free(submsa->ss->tuple_idx);
        submsa->ss->tuple_idx = NULL;
      }
    }
  }

  /* write sub_msa */
  msa_print(F, submsa, output_format, 0);
  fclose(F);

  str_free(outfname);
}

/* print header for summary file */
void write_summary_header(FILE* F, char *alphabet, int gap_strip_mode) {
  int i;
  fprintf(F, "%-20s", "name");
  for (i = 0; i < strlen(alphabet); i++) {
    if (alphabet[i] != GAP_CHAR) {
      if (gap_strip_mode != NO_STRIP) 
        fprintf(F, "         ");
      fprintf(F, "%10c", alphabet[i]);
    }    
  }
  if (gap_strip_mode != NO_STRIP) 
    fprintf(F, "          ");
  fprintf(F, "%10s", "length");
  if (gap_strip_mode != NO_STRIP) 
      fprintf(F, "          "); 
  fprintf(F, "%10s", "any_gaps");
  if (gap_strip_mode != NO_STRIP) 
    fprintf(F, "          ");
  fprintf(F, "%10s", "all_gaps");
  fprintf(F, "\n");    
}

/* print line to summary file; allow for pre- and post-gap-strip statistics */
void write_summary_line(FILE *SUM_F, char *label, char *alphabet, 
                        gsl_vector *freqs, gsl_vector *freqs_strip, int length, 
                        int length_strip, int nallgaps, int nallgaps_strip, 
                        int nanygaps, int nanygaps_strip) {
  int j;
  char tmp[STR_MED_LEN];
  fprintf(SUM_F, "%-20s", label);
  for (j = 0; j < strlen(alphabet); j++) {
    if (alphabet[j] != GAP_CHAR) {
      fprintf(SUM_F, "%10.4f", gsl_vector_get(freqs, j));
      if (freqs_strip != NULL)
        fprintf(SUM_F, " (%6.4f)", gsl_vector_get(freqs_strip, j));
    }    
  }
  fprintf(SUM_F, "%10d", length);
  if (length_strip != -1) {
    sprintf(tmp, "(%d)", length_strip);
    fprintf(SUM_F, " %9s", tmp);
  }

  fprintf(SUM_F, "%10d", nallgaps);
  if (nallgaps_strip != -1) {
    sprintf(tmp, "(%d)", nallgaps_strip);
    fprintf(SUM_F, " %9s", tmp);
  }

  fprintf(SUM_F, "%10d", nanygaps);
  if (nanygaps_strip != -1) {
    sprintf(tmp, "(%d)", nanygaps_strip);
    fprintf(SUM_F, " %9s", tmp);
  }

  fprintf(SUM_F, "\n");
}

/* return 1 if alignment has gaps in all seqs but the reference seq at
   specified position (refseq assumed to be first one); otherwise
   returns 0 */
int all_gaps_but_ref(MSA *msa, int pos) {
  int i;
  if (msa_get_char(msa, 0, pos) == GAP_CHAR) return 0;
  for (i = 1; i < msa->nseqs; i++) 
    if (msa_get_char(msa, i, pos) != GAP_CHAR)
      return 0;
  return 1;
}

/* try to adjust split indices to fall between alignment blocks,
   assuming a reference alignment with respect to the first sequence.
   A block of NSITES_BETWEEN_BLOCKS sites with no gaps in the
   reference sequence and gaps in all other sequences is assumed to
   indicate a region between alignment blocks.  */
void adjust_split_indices_for_blocks(MSA *msa, List *split_indices_list, 
                                     int radius) {  
  int i, j, k = -1, idx, last_idx;
  for (i = 0; i < lst_size(split_indices_list); i++) {
    int new_idx, okay, count, range_beg, range_end;
    idx = lst_get_int(split_indices_list, i) - 1; /* convert to 0-based idx */

    /* first see if we're already in a good place */
    okay = 0;
    j = k = idx;
    for (; j >= 0; j--) { /* look to left */
      if (!all_gaps_but_ref(msa, j)) break;
      if (idx - j + 1 >= NSITES_BETWEEN_BLOCKS) { okay = 1; break; } 
                                /* we're done -- we know we're okay */
    }
    if (j != idx && !okay) {    /* only look right if it's possible
                                   that we're between blocks but we
                                   haven't yet seen enough sites */
      for (; k < msa->length; k++) {
        if (!all_gaps_but_ref(msa, k)) break;
        if (k - j >= NSITES_BETWEEN_BLOCKS) { okay = 1; break; }
      }
    }
    if (okay) continue;         /* no change to index necessary */
        

    /* scan left for a sequence of NSITES_BETWEEN_BLOCKS gaps in all
       seqs but the reference seq */
    range_beg = max(0, idx - radius); 
    range_end = min(msa->length-1, idx + radius);

    /* j currently points to first site equal to or to the left of idx
       s.t. all_gaps_but_ref is false */
    count = 0; new_idx = -1;
    for (; new_idx < 0 && j >= range_beg; j--) {
      if (all_gaps_but_ref(msa, j)) {
        count++;
        if (count == NSITES_BETWEEN_BLOCKS) 
          new_idx = j + NSITES_BETWEEN_BLOCKS / 2 ;
      }
      else count = 0;
    }

    /* k currently points to first site equal to or to the right of
       idx s.t. all_gaps_but_ref is false */
    count = 0;
    for (; new_idx < 0 && k <= range_end; k++) {
      if (all_gaps_but_ref(msa, k)) {
        count++;
        if (count == NSITES_BETWEEN_BLOCKS) 
          new_idx = k - NSITES_BETWEEN_BLOCKS / 2 ;
      }
      else count = 0;
    }

    if (new_idx >= 0)
      lst_set_int(split_indices_list, i, new_idx);
  }

  /* scan for redundant indices */
  last_idx = -1;
  for (i = 0; i < lst_size(split_indices_list); i++) {
    idx = lst_get_int(split_indices_list, i);
    if (idx <= last_idx) lst_delete_idx(split_indices_list, i);
    last_idx = idx;
  }  
}

/* returns List of Strings derived from an argument that may either be
   a literal comma-separated list or a reference to a file (using the
   "*" convention).  Note: List and all Strings are newly allocated
   (should be freed externally) */
/* List *get_arg_list(char *arg) { */
/*   String *argstr = str_new_charstr(arg); */
/*   List *l = lst_new_ptr(10); */
/*   if (str_starts_with_charstr(argstr, "*")) { */
/*     FILE *F; */
/*     String *fname_str; */
/*     if ((F = fopen(&argstr->chars[1], "r")) == NULL) { */
/*       fprintf(stderr, "ERROR: Cannot open file %s.\n", &argstr->chars[1]); */
/*       exit(1); */
/*     } */
/*     fname_str = str_new(STR_MED_LEN); */
/*     str_slurp(fname_str, F); */
/*     str_split(fname_str, NULL, l); */
/*     fclose(F); */
/*     str_free(fname_str); */
/*   } */
/*   else  */
/*     str_split(argstr, ",", l); */

/*   str_free(argstr); */
/*   return l; */
/* } */

int main(int argc, char* argv[]) {
  FILE* F;
  MSA *msa;
  msa_format_type input_format = PHYLIP, output_format = PHYLIP;
  char *msa_fname = NULL, *gff_fname = NULL, *split_indices_str = NULL, 
    *out_fname_root = "msa_split", *cat_map_fname = NULL, *rseq_fname = NULL;
  GFF_Set *gff = NULL;
  int npartitions = -1, feature_frame = -1, strand_sensitive = 0, 
    faithful = 0, partition_frame = 0, quiet_mode = 0, gap_strip_mode = NO_STRIP,
    group_mode = 0, output_summary = 0, tuple_size = 1, win_size = -1, 
    win_overlap = -1, ordered_stats = 1, min_ninf_sites = -1, 
    adjust_radius = -1;
  List *split_indices_list, *cats_to_do = NULL, *order_list = NULL;  
  String *outfname, *sum_fname = NULL;
  FILE *SUM_F = NULL;
  char c;
  int nallgaps, nallgaps_strip, nanygaps, nanygaps_strip, length, 
    length_strip, i;
  gsl_vector *freqs, *freqs_strip;
  msa_coord_map *map = NULL;
  

  while ((c = getopt(argc, argv, "i:M:g:c:p:d:n:sfG:r:o:L:C:T:w:I:O:B:SPzqh")) != -1) {
    switch(c) {
    case 'i':
      if (!strcmp(optarg, "PSU")) input_format = PSU;
      else if (!strcmp(optarg, "FASTA")) input_format = FASTA;
      else if (!strcmp(optarg, "SS")) input_format = SS;
      else if (!strcmp(optarg, "MAF")) input_format = MAF;
      else if (strcmp(optarg, "PHYLIP") != 0) { 
        fprintf(stderr, "ERROR: -i must specify PHYLIP, FASTA, or PSU.  Try msa_split -h for help.\n");
        exit(1); 
      }
      break;
    case 'M':
      rseq_fname = optarg;
      break;
    case 'g':
      gff_fname = optarg;
      break;
    case 'c':
      feature_frame = atoi(optarg);
      break;
    case 'p':
      split_indices_str = optarg;
      break;
    case 'P':
      group_mode = 1;
      break;
    case 'd':
      partition_frame = atoi(optarg);
      break;
    case 'n':
      npartitions = atoi(optarg);
      break;
    case 'w':
      {
        List *l = get_arg_list(optarg);
        if (lst_size(l) != 2 || 
            str_as_int(lst_get_ptr(l, 0), &win_size) != 0 ||
            str_as_int(lst_get_ptr(l, 1), &win_overlap) != 0) {
          fprintf(stderr, "ERROR: illegal arguments to -w.  Try \"msa_split -h\" for help.\n");
          exit(1);
        }
        str_free(lst_get_ptr(l, 0)); str_free(lst_get_ptr(l, 1));
        lst_free(l);
      }
      break;
    case 'L':
      cat_map_fname = optarg;
      break;
    case 's':
      strand_sensitive = 1;
      break;
    case 'f':
      faithful = 1;
      break;
    case 'G':
      if (!strcmp(optarg, "ALL")) gap_strip_mode = STRIP_ALL_GAPS;
      else if (!strcmp(optarg, "ANY")) gap_strip_mode = STRIP_ANY_GAPS;
      else gap_strip_mode = atoi(optarg);
      break;
    case 'r':
      out_fname_root = optarg;
      break;
    case 'T':
      tuple_size = atoi(optarg);
      break;
    case 'z':
      ordered_stats = 0;
      break;
    case 'o':
      if (!strcmp(optarg, "FASTA")) output_format = FASTA;
      else if (!strcmp(optarg, "PSU")) output_format = PSU;
      else if (!strcmp(optarg, "SS")) output_format = SS;
      else if (strcmp(optarg, "PHYLIP") != 0) { 
        fprintf(stderr, "ERROR: -o must specify PHYLIP, FASTA, or PSU.  Try msa_split -h for help.\n");
        exit(1); 
      }
      break;
    case 'O': 
      order_list = get_arg_list(optarg);
      break;
    case 'C':
      {
        List *l = get_arg_list(optarg);
        cats_to_do = lst_new_int(lst_size(l));
        for (i = 0; i < lst_size(l); i++) {
          int cat;
          if (str_as_int(lst_get_ptr(l, i), &cat) != 0 || cat < 0) {
            fprintf(stderr, "ERROR: illegal value in argument to -C.  Try msa_split -h for help.\n");
            exit(1);
          }
          lst_push_int(cats_to_do, cat);
          str_free((String*)lst_get_ptr(l, i));
        }
        lst_free(l);
      }
      break;
    case 'S':
      output_summary = 1;
      break;
    case 'I':
      min_ninf_sites = atoi(optarg);
      break;
    case 'B':
      adjust_radius = atoi(optarg);
      break;
    case 'q':
      quiet_mode = 1;
      break;
    case 'h':
      print_usage();
      exit(0);
      break;
    case '?':
      fprintf(stderr, "ERROR: unrecognized option.  Try \"msa_split -h\" for help.\n");
      exit(1);
    }
  }

  if (optind >= argc) {
    fprintf(stderr, "Missing alignment filename.  Try \"msa_split -h\" for help.\n");
    exit(1);
  }
  else msa_fname = argv[optind];

  if ((strand_sensitive || faithful || feature_frame != -1 || 
       group_mode == 1 || cat_map_fname != NULL) && 
      gff_fname == NULL) {
    fprintf(stderr, "ERROR: -P, -L, -s, -f, and -c require -g.  Try \"msa_split -h\" for help.\n");
    exit(1);
  }

  if ((split_indices_str == NULL && npartitions == -1 
       && group_mode == 0 && cat_map_fname == NULL && win_size == -1) ||
      (split_indices_str != NULL && (npartitions != -1 || group_mode != 0 || 
                                     cat_map_fname != NULL || win_size != -1)) ||
      (npartitions != -1 && (split_indices_str != NULL || group_mode != 0 || 
                             cat_map_fname != NULL || win_size != -1)) ||
      (group_mode != 0 && (split_indices_str != NULL || npartitions != -1 || 
                           cat_map_fname != NULL || win_size != -1)) ||
      (cat_map_fname != NULL && (split_indices_str != NULL || npartitions != -1 || 
                                 group_mode != 0 || win_size != -1))) {
    fprintf(stderr, "ERROR: must specify exactly one of -p, -P, -n, -L, and -w.\nTry \"msa_split -h\" for help.\n");
    exit(1);
  }

  if (npartitions != -1 && npartitions <= 0) {
    fprintf(stderr, "ERROR: number of partitions must be greater than 0.\nTry \"msa_split -h\" for help.\n");
    exit(1);
  }

  if (gff_fname != NULL) {
    if (!quiet_mode)
      fprintf(stderr, "Reading features from %s ...\n", gff_fname);
    if ((F = fopen(gff_fname, "r")) == NULL) {
      fprintf(stderr, "ERROR: cannot open %s.\n", gff_fname);
      exit(1);
    }
    if ((gff = gff_read_set(F)) == NULL) { 
      fprintf(stderr, "ERROR: error reading %s.\n", gff_fname);
      exit(1);
    }
    fclose(F);
  }

  if (!quiet_mode)
    fprintf(stderr, "Reading alignment from %s ...\n", 
            !strcmp(msa_fname, "-") ? "stdin" : msa_fname);

  if (input_format == MAF) {
    assert(gff == NULL);        /* use with GFF not yet supported */
    msa = maf_read(fopen_fname(msa_fname, "r"), 
                   rseq_fname == NULL ? NULL : fopen_fname(rseq_fname, "r"), 
                   NULL, tuple_size, NULL, NULL, -1, 1, 0, NO_STRIP); 
  }
  else 
    msa = msa_new_from_file(fopen_fname(msa_fname, "r"),
                            input_format, NULL); 

  assert(msa->length > 0);

  if (order_list != NULL)
    msa_reorder_rows(msa, order_list);

  if (input_format == SS && msa->ss->tuple_idx == NULL) {
    fprintf(stderr, "ERROR: ordered representation of alignment required.\n");
    exit(1);
  }

  if (cat_map_fname != NULL) {
    CategoryMap *cm;

    if (!quiet_mode)
      fprintf(stderr, "Reading category map from %s ...\n", cat_map_fname);

    cm_read(fopen_fname(cat_map_fname, "r"));

    /* convert GFF to frame of ref of entire alignment */
    msa_map_gff_coords(msa, gff, -1, 0, 0, NULL);

    if (!quiet_mode)
      fprintf(stderr, "Reverse complementing groups on reverse strand ...\n");

    /* FIXME: should this really be done by default? */
    gff_group(gff, "transcript_id"); /* FIXME: parameterize */
    msa_reverse_compl_feats(msa, gff, NULL);

    if (!quiet_mode)
      fprintf(stderr, "Labeling columns of alignment by category ...\n");
    msa_label_categories(msa, gff, cm);

    cm_free(cm);
  }

  split_indices_list = lst_new_int(10);

  if (partition_frame < 0 || partition_frame >= msa->nseqs) {
    fprintf(stderr, "ERROR: illegal argument to -d.  Try \"msa_split -h\" for help.\n");
    exit(1);
  }

  if (partition_frame != 0) 
    map = msa_build_coord_map(msa, partition_frame);

  lst_push_int(split_indices_list, map != NULL ? 
               msa_map_seq_to_msa(map, 1) : 1); 
                                /* always include 1st index */

  if (split_indices_str != NULL) {
    List *l = lst_new_ptr(10);
    String *tmpstr = str_new_charstr(split_indices_str);

    str_split(tmpstr, ",", l);
    
    for (i = 0; i < lst_size(l); i++) {
      int idx;
      String *idxstr = (String*)lst_get_ptr(l, i);
      if (str_as_int(idxstr, &idx) != 0 || idx <= 0) {
        fprintf(stderr, "ERROR: illegal value in argument to -p.  Try \"msa_split -h\" for help.\n");
        exit(1);
      }
      if (map != NULL)          /* convert to frame of entire alignment */
        idx = msa_map_seq_to_msa(map, idx);
      lst_push_int(split_indices_list, idx);
      str_free(idxstr);
    }

    lst_qsort_int(split_indices_list, ASCENDING);

    str_free(tmpstr);
    lst_free(l);
    if (map != NULL) msa_map_free(map);
  }

  if (npartitions > 1) {
    /* NOTE: currently ignores partition frame */
    double split_size = (double)msa->length/npartitions;
    for (i = 1; i < npartitions; i++) { 
                                /* note: will loop npartitions-1
                                   times */
      int idx = (int)ceil(i * split_size);
      lst_push_int(split_indices_list, idx);
    }
  }

  if (win_size != -1) {
    int startidx, mapped_startidx;
    for (startidx = 1 + win_size-win_overlap; startidx < msa->length; 
         startidx += win_size-win_overlap) {
      if (map != NULL) {
        mapped_startidx = msa_map_seq_to_msa(map, startidx);
        if (mapped_startidx == -1) break;
      }
      else mapped_startidx = startidx;
      lst_push_int(split_indices_list, mapped_startidx);
    }
  }

  /* map coords in GFF to frame of ref of alignment */
  if (gff != NULL && cat_map_fname == NULL) {
    if (!quiet_mode)
      fprintf(stderr, "Mapping GFF coordinates to frame of alignment  ...\n");
    msa_map_gff_coords(msa, gff, feature_frame, 0, 0, NULL);
  }

  if (group_mode == 1) {
    GFF_Feature *prevfeat = NULL;
    for (i = 0; i < lst_size(gff->features); i++) {
      GFF_Feature *thisfeat = (GFF_Feature*)lst_get_ptr(gff->features, i);
      if (prevfeat != NULL && 
          !str_equals(thisfeat->attribute, prevfeat->attribute))
        lst_push_int(split_indices_list, prevfeat->end + 
                     (thisfeat->start - prevfeat->end)/2);
      prevfeat = thisfeat;
    }
  }

  if (adjust_radius >= 0)       /* try to adjust split indices to fall
                                   between alignment blocks */
    adjust_split_indices_for_blocks(msa, split_indices_list, adjust_radius);

  /* open summary file for writing, if necessary */
  if (output_summary) {
    sum_fname = str_new_charstr(out_fname_root);
    str_append_charstr(sum_fname, ".sum");
    if ((SUM_F = fopen(sum_fname->chars, "w+")) == NULL) {
      fprintf(stderr, "ERROR: cannot open %s for writing.\n", 
              sum_fname->chars);
      exit(1);
    }

    /* print header */
    write_summary_header(SUM_F, msa->alphabet, gap_strip_mode);
  }

  if (cat_map_fname == NULL) {  /* not using features; splitting
                                   by position (split_indices_list) */
    outfname = str_new(STR_MED_LEN);
    for (i = 0; i < lst_size(split_indices_list); i++) {
      MSA *sub_msa;
      GFF_Set *sub_gff;
      int start = lst_get_int(split_indices_list, i);
      int end;

      end = (i == lst_size(split_indices_list)-1 ? msa->length :
             lst_get_int(split_indices_list, i+1) - 1);

      if (win_size != -1 && end != msa->length) {
        end = (map == NULL ? 
               end + win_overlap :
               msa_map_seq_to_msa(map, msa_map_msa_to_seq(map, end) + 
                                  win_overlap));
        if (end == -1 || end >= msa->length) end = msa->length;
      }

      if (!quiet_mode)
        fprintf(stderr, "Creating partition %d (column %d to column %d) ...\n",
                i+1, start, end);
      
      sub_msa = msa_sub_alignment(msa, NULL, 0, start - 1, end);

      if (map != NULL)
        sub_msa->idx_offset = msa->idx_offset + 
          msa_map_msa_to_seq(map, start) - 1;
                                /* in this case, we'll let the offset
                                   to be wrt the specified reference
                                   sequence */        

      /* collect summary information; do this *before* stripping gaps */
      if (SUM_F != NULL) {
        freqs = msa_get_base_freqs(sub_msa, -1, -1);
        nallgaps = msa_num_gapped_cols(sub_msa, STRIP_ALL_GAPS, -1, -1);
        nanygaps = msa_num_gapped_cols(sub_msa, STRIP_ANY_GAPS, -1, -1);      
        freqs_strip = NULL; nallgaps_strip = -1; nanygaps_strip = -1;
      }

      if (gap_strip_mode != NO_STRIP) {
        strip_gaps(sub_msa, gap_strip_mode);

        /* collect new summary information (post gap strip) */
        if (SUM_F != NULL) {
          freqs_strip = msa_get_base_freqs(sub_msa, -1, -1);
          nallgaps_strip = msa_num_gapped_cols(sub_msa, STRIP_ALL_GAPS, -1, -1);
          nanygaps_strip = msa_num_gapped_cols(sub_msa, STRIP_ANY_GAPS, -1, -1);      
        }
      }

      /* check number of informative sites, if necessary; do after
         stripping gaps */
      if (min_ninf_sites != -1 && 
          msa_ninformative_sites(sub_msa, -1) < min_ninf_sites) {
        fprintf(stderr, "WARNING: skipping partition %d; insufficient informative sites.\n", i+1);
        msa_free(sub_msa);
        continue;
      }
      
      if (gff != NULL) {
        if (gap_strip_mode != NO_STRIP) {
          fprintf(stderr, "ERROR: generation of GFF files for partitions not supported in gap-stripping mode.\n");
          exit(1);
        }

        /* create fname for gff subset */
        sub_gff = gff_subset_range(gff, start, end, !faithful);

        /* reverse complement, if necessary */
        if (strand_sensitive && gff_reverse_strand_only(sub_gff->features)) {
          msa_reverse_compl(sub_msa);
          gff_reverse_compl(sub_gff->features, 
                            faithful ? start : 1, 
                            faithful ? end : end - start + 1);
        }

        /* map coords back to original frame(s) of ref */
        msa_map_gff_coords(faithful ? msa : sub_msa, sub_gff, 0, 
                           feature_frame, 0, NULL);

        /* write gff file for subset */
        str_cpy_charstr(outfname, out_fname_root);
        str_append_charstr(outfname, ".");
        str_append_int(outfname, i+1);
        str_append_charstr(outfname, ".gff");

        if ((F = fopen(outfname->chars, "w+")) == NULL) {
          fprintf(stderr, "ERROR: unable to open %s for writing.\n", 
                  outfname->chars);
          exit(1);
        }
        if (!quiet_mode)
          fprintf(stderr, "Writing GFF subset %d to %s...\n", i+1, 
                  outfname->chars);
        gff_print_set(F, sub_gff);
        fclose(F);
        gff_free_set(sub_gff);
      }

      if (!quiet_mode)
        fprintf(stderr, "Writing partition %d to %s.%d.%s ...\n",
                i+1, out_fname_root, i+1, output_format == SS ? "ss" : "msa");
      write_sub_msa(sub_msa, out_fname_root, i+1, output_format,
                    tuple_size, ordered_stats);

      /* if necessary, write line to summary file */
      if (SUM_F != NULL) {
        char tmpstr[STR_MED_LEN];
        sprintf(tmpstr, "%s.%d.msa", out_fname_root, i+1);

        write_summary_line(SUM_F, tmpstr, sub_msa->alphabet, freqs, 
                           freqs_strip, length, length_strip, nallgaps, 
                           nallgaps_strip, nanygaps, nanygaps_strip);
      }

      msa_free(sub_msa);
    }
    str_free(outfname);
  }
  else {                        /* cat_map_fname != NULL (splitting by
                                   category) */
    List *submsas = lst_new_ptr(10);
    if (!quiet_mode)
      fprintf(stderr, "Partitioning alignment by category ...\n");
    msa_partition_by_category(msa, submsas, cats_to_do, tuple_size);
    for (i = 0; i < lst_size(submsas); i++) {
      MSA *sub = (MSA*)lst_get_ptr(submsas, i);
      int cat = (cats_to_do == NULL ? i : lst_get_int(cats_to_do, i));
      if (lst_size(submsas) == 1) cat = -1;

      /* collect summary information; do this *before* stripping gaps */
      if (SUM_F != NULL) {
        freqs = msa_get_base_freqs(sub, -1, -1);
        length = sub->length;
        nallgaps = msa_num_gapped_cols(sub, STRIP_ALL_GAPS, -1, -1);
        nanygaps = msa_num_gapped_cols(sub, STRIP_ANY_GAPS, -1, -1);      
        freqs_strip = NULL; length_strip = -1; nallgaps_strip = -1; 
        nanygaps_strip = -1; 
      }

      if (!quiet_mode) {
        if (cat >= 0)
          fprintf(stderr, "Writing partition %d to %s.%d.msa ...\n",
                  cat, out_fname_root, cat);
        else
          fprintf(stderr, "Writing partition to %s.msa ...\n",
                  out_fname_root);
      }

      if (gap_strip_mode != NO_STRIP) {
        strip_gaps(sub, gap_strip_mode);

        /* collect new summary information (post gap strip) */
        if (SUM_F != NULL) {
          freqs_strip = msa_get_base_freqs(sub, -1, -1);
          length_strip = sub->length;
          nallgaps_strip = msa_num_gapped_cols(sub, STRIP_ALL_GAPS, -1, -1);
          nanygaps_strip = msa_num_gapped_cols(sub, STRIP_ANY_GAPS, -1, -1);      
        }
      }

      write_sub_msa(sub, out_fname_root, cat, output_format,
                    tuple_size, ordered_stats);

      /* if necessary, write line to summary file */
      if (SUM_F != NULL) {
        char tmpstr[STR_MED_LEN];
        if (cat >= 0)
          sprintf(tmpstr, "%s.%d.msa", out_fname_root, cat);
        else
          sprintf(tmpstr, "%s.msa", out_fname_root);

        write_summary_line(SUM_F, tmpstr, sub->alphabet, freqs, freqs_strip, 
                           length, length_strip, nallgaps, nallgaps_strip, 
                           nanygaps, nanygaps_strip);
      }
      

      msa_free(sub);
    }
    lst_free(submsas);
  }

  if (SUM_F != NULL) {
    if (!quiet_mode) 
      fprintf(stderr, "Writing summary to %s ...\n", sum_fname->chars);
    fclose(SUM_F);
    str_free(sum_fname);
  }

  if (!quiet_mode)
    fprintf(stderr, "Done.\n");

  msa_free(msa);  
  lst_free(split_indices_list);
  if (gff != NULL) gff_free_set(gff);
  if (cats_to_do != NULL) lst_free(cats_to_do);

  return 0;
}

