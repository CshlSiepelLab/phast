/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: msa_split.c,v 1.29 2009-01-09 22:01:00 mt269 Exp $ */

#include <stdlib.h>
#include <stdio.h>
#include "sufficient_stats.h"
#include "msa.h"
#include <getopt.h>
#include <ctype.h>
#include <math.h>
#include "gff.h"
#include "maf.h"

#define DOWNSTREAM_OTHER "other"
#define NSITES_BETWEEN_BLOCKS 30

void print_usage() {
    printf("\n\
USAGE: msa_split [OPTIONS] <fname> \n\
\n\
DESCRIPTION:\n\
\n\
    Partitions a multiple sequence alignment either at designated\n\
    columns, or according to specified category labels, and outputs\n\
    sub-alignments for the partitions.  Optionally splits an\n\
    associated annotations file.\n\
\n\
EXAMPLES:\n\
\n\
    (See below for details on options)\n\
\n\
    1. Read an alignment for a whole human chromosome from a MAF file\n\
    and extract sub-alignments in 1Mb windows overlapping by 1kb.  Use\n\
    sufficient statistics (SS) format for output (can be used by\n\
    phyloFit, phastCons, or exoniphy).  Set window boundaries between\n\
    alignment blocks, if possible.\n\
\n\
        msa_split chr1.maf --refseq chr1.fa --in-format MAF \\\n\
            --windows 1000000,1000 --out-format SS \\\n\
            --between-blocks 5000 --out-root chr1\n\
\n\
    (Windows will be defined using the coordinate system of the first\n\
    sequence in the alignment, assumed to be the reference sequence;\n\
    output will be to chr1.1-1000000.ss, chr1.999001-1999000.ss, ...)\n\
\n\
    2. As in (1), but report unordered sufficient statistics (much\n\
    more compact and adequate for use with phyloFit).\n\
\n\
        msa_split chr1.maf --refseq chr1.fa --in-format MAF \\\n\
            --windows 1000000,1000 --out-format SS \\\n\
            --between-blocks 5000 --out-root chr1 --unordered-ss\n\
\n\
    3. Extract sub-alignments of sites in conserved elements and not\n\
    in conserved elements, as defined by a BED file (coordinates\n\
    assumed to be for 1st sequence).  Read multiple alignment in FASTA\n\
    format.\n\
\n\
        msa_split mydata.fa --features conserved.bed --by-category \\\n\
            --out-root mydata\n\
\n\
    (Output will be to mydata.background-0.fa and mydata.bed_feature-1.fa\n\
    [latter has sites of category number 1, defined by bed file]\n\
\n\
    3. Extract sub-alignments of sites in each of the three codon\n\
    positions, as defined by a GFF file (coordinates assumed to be for\n\
    1st sequence).  Reverse complement genes on minus strand.\n\
\n\
        msa_split chr22.maf --in-format MAF --features chr22.gff \\\n\
            --by-category --catmap \"NCATS 3 ; CDS 1-3\" --do-cats CDS \\\n\
            --reverse-compl --out-root chr22 --out-format SS\n\
\n\
    (Output will be to chr22.cds-1.ss, chr22.cds-2.ss, chr22.cds-3.ss)\n\
\n\
    4. Split an alignment into pieces corresponding to the genes in a\n\
    GFF file.  Assume genes are defined by the tag \"transcript_id\".\n\
\n\
        msa_split cftr.fa --features cftr.gff --by-group transcript_id\n\
\n\
    5. Obtain a sub-alignment for each of a set of regulatory regions,\n\
    as defined in a BED file.\n\
\n\
        msa_split chr22.maf --in-format MAF --refseq chr22.fa \\\n\
	    --features chr22.reg.bed --for-features \\\n\
	    --out-root chr22.reg\n\
\n\
OPTIONS:\n\
\n\
 (Splitting options)\n\
    --windows, -w <win_size,win_overlap>\n\
        Split the alignment into \"windows\" of size <win_size> bases,\n\
        overlapping by <win_overlap>.\n\
\n\
    --by-category, -L\n\
        (Requires --features) Split by category, as defined by\n\
        annotations file and (optionally) category map (see\n\
        --catmap)\n\
\n\
    --by-group, -P <tag>\n\
        (Requires --features) Split by groups in annotation file,\n\
        as defined by specified tag.  Splits midway between every\n\
        pair of consecutive groups.  Features will be sorted by group.\n\
        There should be no overlapping features (see 'refeature\n\
        --unique').\n\
\n\
    --for-features, -F\n\
        (Requires --features) Extract section of alignment\n\
        corresponding to every feature.  There will be no output for\n\
        regions not covered by features.\n\
\n\
    --by-index, -p <indices>\n\
        List of explicit indices at which to split alignment\n\
        (comma-separated).  If the list of indices is \"10,20\",\n\
        then sub-alignments will be output for sites 1-9, 10-19, and\n\
        20-<msa_len>.  Note that the indices are relative to the \n\
        input alignment, and not necessarily in genomic coordinates.\n\
\n\
    --npartitions, -n <number>\n\
        Split alignment equally into specified number of partitions.\n\
\n\
    --between-blocks, -B <radius>\n\
        (Not for use with --by-category or --for-features) Try to\n\
        partition at sites between alignment blocks.  Assumes a\n\
        reference sequence alignment, with the first sequence as the\n\
        reference seq (as created by multiz).  Blocks of %d sites with\n\
        gaps in all sequences but the reference seq are assumed to\n\
        indicate boundaries between alignment blocks.  Partition\n\
        indices will not be moved more than <radius> sites.\n\
\n\
    --features, -g <fname>\n\
        (For use with --by-category, --by-group, --for-features, or \n\
	--windows) Annotations file.  May be GFF, BED, or genepred\n\
        format.  Coordinates are assumed to be in the coordinate frame of \n\
        the first sequence in the alignment (assumed to be the reference\n\
        sequence).\n\
\n\
    --catmap, -c <fname>|<string>\n\
        (Optionally use with --by-category) Mapping of feature types\n\
        to category numbers.  Can either give a filename or an\n\
        \"inline\" description of a simple category map, e.g.,\n\
        --catmap \"NCATS = 3 ; CDS 1-3\" or --catmap \"NCATS = 1 ; UTR\n\
        1\".\n\
\n\
    --refidx, -d <frame_index>\n\
        (For use with --windows or --by-index) Index of frame of\n\
        reference for split indices.  Default is 1 (1st sequence\n\
        assumed reference).\n\
\n\
 (File names & formats, type of output, etc.)\n\
    --in-format, -i FASTA|PHYLIP|MPM|MAF|SS\n\
        Input alignment file format.  Default is to guess format from \n\
        file contents.\n\
\n\
    --refseq, -M <fname>\n\
        (For use with --in-format MAF) Name of file containing\n\
        reference sequence, in FASTA format.\n\
\n\
    --out-format, -o FASTA|PHYLIP|MPM|SS\n\
        Output alignment file format.  Default is FASTA.\n\
\n\
    --out-root, -r <name>\n\
        Filename root for output files (default \"msa_split\").\n\
\n\
    --sub-features, -f\n\
	(For use with --features)  Output subsets of features\n\
	corresponding to subalignments.  Features overlapping\n\
	partition boundaries will be discarded.  Not permitted with\n\
	--by-category.\n\
\n\
    --reverse-compl, -s\n\
        Reverse complement all segments having at least one feature on\n\
        the reverse strand and none on the positive strand.  For use\n\
        with --by-group.  Can also be used with --by-category to ensure\n\
        all sites in a category are represented in the same strand\n\
        orientation.\n\
\n\
    --gap-strip, -G ALL|ANY|<seqno>\n\
        Strip columns in output alignments containing all gaps, any\n\
        gaps, or gaps in the specified sequence (<seqno>; indexing\n\
        begins with one).  Default is not to strip any columns.\n\
\n\
    --seqs, -l <seq_list>\n\
        Include only specified sequences in output.  Indicate by \n\
        sequence number or name (numbering starts with 1 and is\n\
        evaluated *after* --order is applied).\n\
\n\
    --exclude, -x\n\
        Exclude rather than include specified sequences.\n\
\n\
    --order, -O <name_list>\n\
        Change order of rows in alignment to match sequence names\n\
        specified in name_list.  If a name appears in name_list but\n\
        not in the alignment, a row of gaps will be inserted.\n\
\n\
    --min-informative, -I <n> \n\
        Only output alignments having at least <n> informative sites\n\
        (sites at which at least two non-gap and non-N gaps are present).\n\
\n\
    --do-cats, -C <cat_list>\n\
        (For use with --by-category) Output sub-alignments for only the\n\
        specified categories (column-delimited list).\n\
\n\
    --tuple-size, -T <tuple_size>\n\
        (for use with --by-category or --out-format SS) Size of tuples\n\
        of columns to consider in downstream analysis (e.g., with\n\
        context-dependent phylogenetic models; see 'phyloFit').  With\n\
        --by-category, insert tuple_size-1 columns of missing data\n\
        between sites that were not adjacent in the original alignment,\n\
        to avoid creating artificial context.  With --out-format SS,\n\
        express sufficient statistics in terms of tuples of specified size.\n\
\n\
    --unordered-ss, -z  \n\
        (For use with --out-format SS)  Suppress the portion of the\n\
        sufficient statistics concerned with the order in which columns\n\
        appear.\n\
\n\
    --summary, -S \n\
        Output summary of each output alignment to a file with suffix\n\
        \".sum\" (includes base frequencies and numbers of gapped columns).\n\
\n\
 (Other)\n\
    --quiet, -q\n\
        Proceed quietly.\n\
\n\
    --help, -h\n\
        Print this help message.\n\n", NSITES_BETWEEN_BLOCKS);
}

void write_sub_msa(MSA *submsa, char *fname, msa_format_type output_format, 
                   int tuple_size, int ordered_stats) {
  FILE *F = phast_fopen(fname, "w+");

  /* create sufficient stats, if necessary */
  if (output_format == SS) {
    if (submsa->ss == NULL)
      ss_from_msas(submsa, tuple_size, ordered_stats, NULL, NULL, NULL, -1, 0);
    else if (submsa->ss->tuple_size != tuple_size) 
      die("ERROR: tuple size in SS file does not match desired tuple size for output.\nConversion not supported.\n");
    ss_write(submsa, F, ordered_stats);
  }
  else 
    msa_print(F, submsa, output_format, 0);

  phast_fclose(F);
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
                        Vector *freqs, Vector *freqs_strip, int length, 
                        int length_strip, int nallgaps, int nallgaps_strip, 
                        int nanygaps, int nanygaps_strip) {
  int j;
  char tmp[STR_MED_LEN];
  fprintf(SUM_F, "%-20s", label);
  for (j = 0; j < strlen(alphabet); j++) {
    if (alphabet[j] != GAP_CHAR) {
      fprintf(SUM_F, "%10.4f", vec_get(freqs, j));
      if (freqs_strip != NULL)
        fprintf(SUM_F, " (%6.4f)", vec_get(freqs_strip, j));
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

/* try to adjust split indices to fall between alignment blocks,
   assuming a reference alignment with respect to the first sequence.
   A block of NSITES_BETWEEN_BLOCKS sites with no gaps in the
   reference sequence and missing data in all other sequences is assumed to
   indicate a region between alignment blocks.  */
void adjust_split_indices_for_blocks(MSA *msa, List *split_indices_list, 
                                     int radius) {  
  int i, j, k = -1, idx, last_idx;
  for (i = 0; i < lst_size(split_indices_list); i++) {
    int new_idx, okay, count, range_beg, range_end;
    idx = lst_get_int(split_indices_list, i) - 1; /* convert to 0-based idx */

    if (idx == 0) continue;     /* don't do this at the beginning of
                                   an alignment */

    /* first see if we're already in a good place */
    /* NOTE: it's best to avoid breaking where there's an alignment
       gap in the reference sequence (can happen when there are Ns in
       other seqs); these places can be a problem in certain cases,
       e.g., with win-overlap 0, because multiple indices in the full
       alignment map to the same index in the reference sequence */
    okay = 0;
    j = k = idx;
    for (; j >= 0; j--) { /* look to left */
      if (msa_get_char(msa, 0, j) == GAP_CHAR ||
          !msa_missing_col(msa, 1, j)) break;
      if (idx - j + 1 >= NSITES_BETWEEN_BLOCKS) { okay = 1; break; } 
                                /* we're done -- we know we're okay */
    }
    if (j != idx && !okay) {    /* only look right if it's possible
                                   that we're between blocks but we
                                   haven't yet seen enough sites */
      for (; k < msa->length; k++) {
        if (msa_get_char(msa, 0, k) == GAP_CHAR || 
            !msa_missing_col(msa, 1, k)) break;
        if (k - j >= NSITES_BETWEEN_BLOCKS) { okay = 1; break; }
      }
    }
    if (okay) continue;         /* no change to index necessary */
        

    /* scan left for a sequence of NSITES_BETWEEN_BLOCKS missing-data
       columns */
    range_beg = max(0, idx - radius); 
    range_end = min(msa->length-1, idx + radius);

    /* j currently points to first non-missing-data col equal to or to the
       left of idx */
    count = 0; new_idx = -1;
    for (; new_idx < 0 && j >= range_beg; j--) {
      if (msa_missing_col(msa, 1, j) && msa_get_char(msa, 0, j) != GAP_CHAR) {
        count++;
        if (count == NSITES_BETWEEN_BLOCKS) 
          new_idx = j + NSITES_BETWEEN_BLOCKS / 2 ;
      }
      else count = 0;
    }

    /* k currently points to first non-missing-data col equal to or to the
       right of idx */
    count = 0;
    for (; new_idx < 0 && k <= range_end; k++) {
      if (msa_missing_col(msa, 1, k) && msa_get_char(msa, 0, k) != GAP_CHAR) {
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

int main(int argc, char* argv[]) {
  FILE* F;
  MSA *msa;
  msa_format_type input_format = UNKNOWN_FORMAT, output_format = FASTA;
  char *msa_fname = NULL, *split_indices_str = NULL, 
    *out_fname_root = "msa_split", *rseq_fname = NULL, *group_tag = NULL;
  GFF_Set *gff = NULL;
  int npartitions = -1, strand_sensitive = 0, 
    partition_frame = 1, quiet_mode = 0, gap_strip_mode = NO_STRIP,
    output_summary = 0, tuple_size = 1, win_size = -1, 
    win_overlap = -1, ordered_stats = 1, min_ninf_sites = -1, 
    adjust_radius = -1, opt_idx, by_category = FALSE, for_features = FALSE,
    exclude_seqs = FALSE, sub_features = FALSE;
  List *split_indices_list, *cats_to_do_str = NULL, *order_list = NULL, 
    *segment_ends_list = NULL, *seqlist_str = NULL, *seqlist = NULL, 
    *cats_to_do = NULL;  
  String *sum_fname = NULL;
  FILE *SUM_F = NULL;
  char c;
  int nallgaps, nallgaps_strip, nanygaps, nanygaps_strip, 
    length_strip, i;
  Vector *freqs, *freqs_strip;
  msa_coord_map *map = NULL;
  CategoryMap *cm = NULL;
  char subfname[STR_MED_LEN];

  struct option long_opts[] = {
    {"windows", 1, 0, 'w'},
    {"by-category", 0, 0, 'L'},
    {"by-group", 1, 0, 'P'},
    {"for-features", 0, 0, 'F'},
    {"by-index", 1, 0, 'p'},
    {"npartitions", 1, 0, 'n'},
    {"between-blocks", 1, 0, 'B'},
    {"features", 1, 0, 'g'},
    {"catmap", 1, 0, 'c'},
    {"refidx", 1, 0, 'd'},
    {"in-format", 1, 0, 'i'},
    {"refseq", 1, 0, 'M'},
    {"out-format", 1, 0, 'o'},
    {"out-root", 1, 0, 'r'},
    {"sub-features", 0, 0, 'f'},
    {"reverse-compl", 0, 0, 's'},
    {"gap-strip", 1, 0, 'G'},
    {"seqs", 1, 0, 'l'},
    {"exclude", 0, 0, 'x'},
    {"order", 1, 0, 'O'},
    {"min-informative", 1, 0, 'I'},
    {"do-cats", 1, 0, 'C'},
    {"tuple-size", 1, 0, 'T'},
    {"unordered-ss", 0, 0, 'z'},
    {"summary", 0, 0, 'S'},
    {"quiet", 0, 0, 'q'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "i:M:g:c:p:d:n:sfG:r:o:L:C:T:w:I:O:B:P:F:l:xSzqh", long_opts, &opt_idx)) != -1) {
    switch(c) {
    case 'i':
      input_format = msa_str_to_format(optarg);
      if (input_format == UNKNOWN_FORMAT) die("ERROR: bad input alignment format.\n");
      break;
    case 'M':
      rseq_fname = optarg;
      break;
    case 'g':
      gff = gff_read_set(phast_fopen(optarg, "r"));
      break;
    case 'c':
      cm = cm_new_string_or_file(optarg);
      break;
    case 'p':
      split_indices_str = optarg;
      break;
    case 'P':
      group_tag = optarg;
      break;
    case 'd':
      partition_frame = get_arg_int(optarg);
      break;
    case 'n':
      npartitions = get_arg_int(optarg);
      break;
    case 'w':
      {
        List *l = get_arg_list(optarg);
        if (lst_size(l) != 2 || 
            str_as_int(lst_get_ptr(l, 0), &win_size) != 0 ||
            str_as_int(lst_get_ptr(l, 1), &win_overlap) != 0) 
          die("ERROR: illegal arguments to -w.  Try \"msa_split -h\" for help.\n");
        str_free(lst_get_ptr(l, 0)); str_free(lst_get_ptr(l, 1));
        lst_free(l);
      }
      break;
    case 'L':
      by_category = TRUE;
      break;
    case 's':
      strand_sensitive = 1;
      break;
    case 'G':
      if (!strcmp(optarg, "ALL")) gap_strip_mode = STRIP_ALL_GAPS;
      else if (!strcmp(optarg, "ANY")) gap_strip_mode = STRIP_ANY_GAPS;
      else gap_strip_mode = get_arg_int(optarg);
      break;
    case 'l':
      seqlist_str = get_arg_list(optarg);
      break;
    case 'x':
      exclude_seqs = TRUE;
      break;
    case 'r':
      out_fname_root = optarg;
      break;
    case 'T':
      tuple_size = get_arg_int(optarg);
      break;
    case 'z':
      ordered_stats = 0;
      break;
    case 'o':
      output_format = msa_str_to_format(optarg);
      if (output_format == UNKNOWN_FORMAT) die("ERROR: bad output alignment format.\n");
      break;
    case 'O': 
      order_list = get_arg_list(optarg);
      break;
    case 'C':
      cats_to_do_str = get_arg_list(optarg);
      break;
    case 'F':
      for_features = TRUE;
      break;
    case 'f':
      sub_features = TRUE;
      break;
    case 'S':
      output_summary = 1;
      break;
    case 'I':
      min_ninf_sites = get_arg_int(optarg);
      break;
    case 'B':
      adjust_radius = get_arg_int(optarg);
      break;
    case 'q':
      quiet_mode = 1;
      break;
    case 'h':
      print_usage();
      exit(0);
      break;
    case '?':
      die("ERROR: unrecognized option.  Try \"msa_split -h\" for help.\n");
    }
  }

  if (optind >= argc) 
    die("Missing alignment filename.  Try \"msa_split -h\" for help.\n");

  set_seed(-1);

  msa_fname = argv[optind];

  if (split_indices_str == NULL && npartitions == -1 && !by_category && 
      win_size == -1 && !for_features && group_tag == NULL)
    die ("ERROR: must specify one of --windows, --by-category, --by-group,\n--for-features, --by-index, and --npartitions.  Try \"msa_split -h\" for help.\n");

  if ((split_indices_str != NULL && (npartitions != -1 || by_category || 
                                     win_size != -1 || for_features || 
                                     group_tag != NULL)) ||
      (npartitions != -1 && (by_category || win_size != -1 || for_features ||
                             group_tag != NULL)) ||
      (by_category && (win_size != -1 || for_features || group_tag != NULL)) ||
      (win_size != -1 && (for_features || group_tag != NULL)) ||
      (for_features && group_tag != NULL))
    die("ERROR: cannot select more than one of --windows, --by-category, --by-group,\n--for-features, --by-index, and --npartitions.  Try \"msa_split -h\" for help.\n");

  if ((group_tag != NULL || by_category || for_features) && gff == NULL) 
    die("ERROR: --by_category, --by-group, and --for-features require --features.  Try \"msa_split -h\" for help.\n");

  if (strand_sensitive && group_tag == NULL && !by_category)
    die("ERROR: --reverse-compl requires --by-group or --by-category.\n");

  if (npartitions != -1 && npartitions <= 0) 
    die("ERROR: number of partitions must be greater than 0.\nTry \"msa_split -h\" for help.\n");

  if (adjust_radius >= 0 && (for_features || by_category))
    die("ERROR: can't use --between-blocks with --by-category or --for-features.\nTry \"msa_split -h\" for help.\n");

  if (!quiet_mode)
    fprintf(stderr, "Reading alignment from %s...\n", 
            !strcmp(msa_fname, "-") ? "stdin" : msa_fname);

  if (gff != NULL && cm == NULL) cm = cm_new_from_features(gff);

  if (cats_to_do_str != NULL) cats_to_do = cm_get_category_list(cm, cats_to_do_str, FALSE);

  FILE *infile = phast_fopen(msa_fname, "r");
  if (input_format == UNKNOWN_FORMAT)
    input_format = msa_format_for_content(infile, 1);
  if (input_format == MAF) {
    if (gff != NULL) fprintf(stderr, "WARNING: use of --features with a MAF file currently forces a projection onto the reference sequence.\n");

    msa = maf_read_cats(infile, 
                   rseq_fname == NULL ? NULL : phast_fopen(rseq_fname, "r"), 
                   tuple_size, NULL, gff, cm, -1, TRUE, NULL, NO_STRIP, FALSE,
		   cats_to_do); 
    /* NOTE: no support yet for reverse complementing groups on
       reverse strand in MAF case */
  }
  else {
    msa = msa_new_from_file_define_format(infile,
                            input_format, NULL); 

    if (gff != NULL) {
      if (!quiet_mode)
	fprintf(stderr, "Mapping feature coordinates to frame of alignment...\n");
      msa_map_gff_coords(msa, gff, 1, 0, 0);

      if (strand_sensitive) {
        if (!quiet_mode)
          fprintf(stderr, "Reverse complementing groups on reverse strand...\n");
        msa_reverse_compl_feats(msa, gff, NULL);
      }
    }
  }
  if (msa->length <= 0) 
    die("ERROR: msa->length is %i\n", msa->length);

  if (gff != NULL) {
    if (!quiet_mode)
      fprintf(stderr, "Labeling columns of alignment by category...\n");
    msa_label_categories(msa, gff, cm);
  }

  if (order_list != NULL)
    msa_reorder_rows(msa, order_list);

  if (seqlist_str != NULL) 
    seqlist = msa_seq_indices(msa, seqlist_str);

  if (input_format == SS && msa->ss->tuple_idx == NULL) 
    die("ERROR: ordered representation of alignment required.\n");

  split_indices_list = lst_new_int(10);

  if (partition_frame < 0 || partition_frame >= msa->nseqs) 
    die("ERROR: illegal argument to -d.  Try \"msa_split -h\" for help.\n");

  if (partition_frame != 0) 
    map = msa_build_coord_map(msa, partition_frame);

  if (!for_features)            /* always add position 1 to list
                                   *unless* --for-features */
    //    lst_push_int(split_indices_list, map != NULL ? 
    //                 msa_map_seq_to_msa(map, 1) : 1); 
    //bug fix: don't use map here, since this will cut off gaps in reference seq at
    //beginning of alignment
    lst_push_int(split_indices_list, 1);

  if (split_indices_str != NULL) { /* --by-index */
    List *l = lst_new_ptr(10);
    String *tmpstr = str_new_charstr(split_indices_str);

    str_split(tmpstr, ",", l);
    
    for (i = 0; i < lst_size(l); i++) {
      int idx;
      String *idxstr = (String*)lst_get_ptr(l, i);
      if (str_as_int(idxstr, &idx) != 0 || idx <= 0) 
        die("ERROR: illegal value in argument to -p.  Try \"msa_split -h\" for help.\n");

      if (map != NULL)          /* convert to frame of entire alignment */
        idx = msa_map_seq_to_msa(map, idx);

      lst_push_int(split_indices_list, idx);
      str_free(idxstr);
    }

    lst_qsort_int(split_indices_list, ASCENDING);

    str_free(tmpstr);
    lst_free(l);
  }

  else if (npartitions > 1) {        /* --npartitions */
    /* NOTE: currently ignores partition frame */
    double split_size = (double)msa->length/npartitions;
    for (i = 1; i < npartitions; i++) {  
                                /* note: will loop npartitions-1
                                   times */
      int idx = (int)ceil(i * split_size);
      lst_push_int(split_indices_list, idx);
    }
  }

  else if (win_size != -1) {         /* --windows */
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

  else if (group_tag != NULL) {    /* --by-group */
    GFF_FeatureGroup *prevg = NULL;
    gff_group(gff, group_tag);
    gff_sort(gff);

    for (i = 0; i < lst_size(gff->groups); i++) {
      GFF_FeatureGroup *g = lst_get_ptr(gff->groups, i);
      if (g->start < 0 || g->end < 0) {
	fprintf(stderr, "Ignoring group '%s' (bad coordinates).\n", g->name->chars);
	continue;
      }
      if (prevg != NULL) {
	if (prevg->end >= g->start)
	  die("ERROR: feature groups overlap.\n");
	lst_push_int(split_indices_list, prevg->end + 
		     (g->start - prevg->end)/2);
      }
      prevg = g;
    }
  }

  else if (for_features) {    /* --for-features */
    segment_ends_list = lst_new_int(lst_size(gff->features));
      
    for (i = 0; i < lst_size(gff->features); i++) {
      GFF_Feature *f = lst_get_ptr(gff->features, i);
      lst_push_int(split_indices_list, f->start);
      lst_push_int(segment_ends_list, f->end);
    }
  }


  if (adjust_radius >= 0)       /* try to adjust split indices to fall
                                   between alignment blocks */
    adjust_split_indices_for_blocks(msa, split_indices_list, adjust_radius);

  /* open summary file for writing, if necessary */
  if (output_summary) {
    sum_fname = str_new_charstr(out_fname_root);
    str_append_charstr(sum_fname, ".sum");
    phast_fopen(sum_fname->chars, "w+");

    /* print header */
    write_summary_header(SUM_F, msa->alphabet, gap_strip_mode);
  }

  if (!by_category) {           /* splitting by position
                                   (split_indices_list) */
    msa_free_categories(msa);
    for (i = 0; i < lst_size(split_indices_list); i++) {
      MSA *sub_msa;
      GFF_Set *sub_gff;
      int start = lst_get_int(split_indices_list, i);
      int orig_start = map == NULL ? start : msa_map_msa_to_seq(map, start);
                                /* keep track of orig. coords also --
                                   report these to user */
      int end, orig_end;

      /* patch for rare bug that occurs when start coincides with gap in
         reference sequence: orig_start will be off by one in this
         case because of the way the coord mapping is done */
      if (map != NULL && 
          msa_get_char(msa, partition_frame-1, start - 1) == GAP_CHAR) {
	//if orig_start == -1 this is presumably because we are at
	//very beginning of alignment and coordinate map didn't include this
	//position
	if (orig_start == -1) orig_start = 1;
        else orig_start++;
      }
      
      if (segment_ends_list == NULL) {
        end = (i == lst_size(split_indices_list)-1 ? msa->length :
               lst_get_int(split_indices_list, i+1) - 1);
        if (win_size != -1 && end != msa->length) 
          end = (map == NULL ? 
                 end + win_overlap :
                 msa_map_seq_to_msa(map, msa_map_msa_to_seq(map, end) + 
                                    win_overlap));
      }
      else 
        end = lst_get_int(segment_ends_list, i);

      if (end == -1 || end > msa->length) end = msa->length;
      orig_end = map == NULL ? end : msa_map_msa_to_seq(map, end);

      if (!quiet_mode)
        fprintf(stderr, "Creating partition %d (column %d to column %d)...\n",
                i+1, orig_start, orig_end);
      
      sub_msa = msa_sub_alignment(msa, seqlist, !exclude_seqs, start - 1, end);

      if (map != NULL)
        sub_msa->idx_offset = msa->idx_offset + orig_start - 1;
                                /* in this case, we'll let the offset
                                   be wrt the specified reference
                                   sequence */        

      /* collect summary information; do this *before* stripping gaps */
      if (SUM_F != NULL) {
        freqs = msa_get_base_freqs(sub_msa, -1, -1);
        nallgaps = msa_num_gapped_cols(sub_msa, STRIP_ALL_GAPS, -1, -1);
        nanygaps = msa_num_gapped_cols(sub_msa, STRIP_ANY_GAPS, -1, -1);      
        freqs_strip = NULL; nallgaps_strip = -1; nanygaps_strip = -1;
      }

      if (gap_strip_mode != NO_STRIP) {
        msa_strip_gaps(sub_msa, gap_strip_mode);

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
      
      sprintf(subfname, "%s.%d-%d.%s", out_fname_root, orig_start, 
              orig_end, msa_suffix_for_format(output_format));

      /* Avoid complaints about msa->seqs and msa->categories being non-null
	 when using maf input sequences */
      //      sub_msa->seqs = NULL;
      sub_msa->categories = NULL;
      sub_msa->ncats = -1;

      if (!quiet_mode)
        fprintf(stderr, "Writing partition %d to %s...\n", i+1, subfname);
      write_sub_msa(sub_msa, subfname, output_format, tuple_size, ordered_stats);

      if (SUM_F != NULL) 
        write_summary_line(SUM_F, subfname, sub_msa->alphabet, freqs, 
                           freqs_strip, sub_msa->length, -1, nallgaps, 
                           nallgaps_strip, nanygaps, nanygaps_strip);

      if (sub_features && gff != NULL) {
        if (gap_strip_mode != NO_STRIP) 
          die("ERROR: generation of GFF files for partitions not supported in gap-stripping mode.\n");

        /* create fname for gff subset */
        sub_gff = gff_subset_range(gff, start, end, TRUE);

	if (lst_size(sub_gff->features) == 0) {
	  if (!quiet_mode)
	    fprintf(stderr, "(No features for subset %d)\n", i+1);
	}
	else {  /* write gff file for subset */
	  /* map coords back to original frame(s) of ref */
	  msa_map_gff_coords(sub_msa, sub_gff, 0, 1, 
			     output_format == SS ? sub_msa->idx_offset : 0);
			     /* if output SS, add offset */

	  sprintf(subfname, "%s.%d-%d.gff", out_fname_root, orig_start, orig_end);
	  F = phast_fopen(subfname, "w+");
	  if (!quiet_mode)
	    fprintf(stderr, "Writing GFF subset %d to %s...\n", i+1, subfname);
	  gff_print_set(F, sub_gff);
	  phast_fclose(F);
	}
	gff_free_set(sub_gff);
      }

      msa_free(sub_msa);
    }
  }
  else {                        /* by_category == TRUE */
    List *submsas = lst_new_ptr(10);
    if (!quiet_mode)
      fprintf(stderr, "Partitioning alignment by category...\n");
    msa_partition_by_category(msa, submsas, cats_to_do, tuple_size);
    for (i = 0; i < lst_size(submsas); i++) {
      MSA *sub = (MSA*)lst_get_ptr(submsas, i);
      int cat = (cats_to_do == NULL ? i : lst_get_int(cats_to_do, i));
      
      /* collect summary information; do this *before* stripping gaps */
      if (SUM_F != NULL) {
        freqs = msa_get_base_freqs(sub, -1, -1);
        nallgaps = msa_num_gapped_cols(sub, STRIP_ALL_GAPS, -1, -1);
        nanygaps = msa_num_gapped_cols(sub, STRIP_ANY_GAPS, -1, -1);      
        freqs_strip = NULL; length_strip = -1; nallgaps_strip = -1; 
        nanygaps_strip = -1; 
      }

      if (gap_strip_mode != NO_STRIP) {
        msa_strip_gaps(sub, gap_strip_mode);

        /* collect new summary information (post gap strip) */
        if (SUM_F != NULL) {
          freqs_strip = msa_get_base_freqs(sub, -1, -1);
          length_strip = sub->length;
          nallgaps_strip = msa_num_gapped_cols(sub, STRIP_ALL_GAPS, -1, -1);
          nanygaps_strip = msa_num_gapped_cols(sub, STRIP_ANY_GAPS, -1, -1);      
        }
      }

      sprintf(subfname, "%s.%s-%d.%s", out_fname_root, 
              ((String*)cm_get_feature(cm, cat))->chars,
              cat, msa_suffix_for_format(output_format));

      if (!quiet_mode) 
        fprintf(stderr, "Writing partition (category) %d to %s...\n", 
                cat, subfname);

      write_sub_msa(sub, subfname, output_format, tuple_size, ordered_stats);

      /* if necessary, write line to summary file */
      if (SUM_F != NULL) 
        write_summary_line(SUM_F, subfname, sub->alphabet, freqs, freqs_strip, 
                           sub->length, length_strip, nallgaps, nallgaps_strip, 
                           nanygaps, nanygaps_strip);

      msa_free(sub);
    }
    lst_free(submsas);
  }

  if (SUM_F != NULL) {
    if (!quiet_mode) 
      fprintf(stderr, "Writing summary to %s...\n", sum_fname->chars);
    phast_fclose(SUM_F);
  }

  if (!quiet_mode)
    fprintf(stderr, "Done.\n");

  return 0;
}

