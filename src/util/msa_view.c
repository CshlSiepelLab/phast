/* $Id: msa_view.c,v 1.15 2004-07-27 20:39:41 acs Exp $
   Written by Adam Siepel, 2002
   Copyright 2002, Adam Siepel, University of California */

#include <stdlib.h>
#include <stdio.h>
#include <msa.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include <string.h>
#include <ctype.h>
#include <sufficient_stats.h>
#include <local_alignment.h>
#include <maf.h>

/* minimum number of codons required for -L */
#define MIN_NCODONS 10

/* number of bases considered to border each indel and minimum number
   of gapless bases required for -I */
#define INDEL_BORDER 3
#define MIN_NBASES 15

void print_usage() {
    printf("\n\
USAGE: msa_view [OPTIONS] <infile>\n\
\n\
DESCRIPTION:\n\
\n\
    Provides various kinds of \"views\" of one or more multiple\n\
    alignments.  Can extract a sub-alignment from an alignment (by row\n\
    or by column) or combine several alignments into one.  Also can\n\
    extract the sufficient statistics for phylogenetic analysis from\n\
    an alignment, optionally accounting for site categories that are\n\
    defined by an auxiliary annotations file.  Supports various other\n\
    functions, including gap stripping, column randomization, and\n\
    reordering of sequences.  Capable of reading and writing in a few\n\
    common formats.  Can be used for file conversion (by default,\n\
    output is the entire input alignment).\n\
\n\
EXAMPLES:\n\
\n\
    (See below for more details on options)\n\
\n\
    1. Convert alignment formats (default input and output is FASTA)\n\
\n\
        msa_view myfile.fa --out-format PHYLIP > myfile.ph\n\
        msa_view myfile2.raw --in-format MPM > myfile2.fa\n\
\n\
    2. Obtain a sub-alignment by position, using the coordinate frame\n\
    of the first sequence in the alignment.\n\
\n\
        msa_view myfile.fa --start 1234 --end 5678 --refidx 1 > mysub.fa\n\
\n\
    3. Obtain a sub-alignment by sequence.\n\
\n\
        msa_view myfile.fa --seqs 1,4,5 > seqs145.fa\n\
        msa_view myfile.fa --seqs 1,4,5 --exclude > seqs236.fa\n\
\n\
    (can also specify sequences by name, e.g., --seqs cow,rat,pig)\n\
\n\
    4. Concatenate alignments. \n\
\n\
        msa_view --aggregate human,mouse,rat myf1.fa myf2.fa myf3.fa \n\
            > concat.fa\n\
\n\
    (source alignments may have different subsets of sequences and may\n\
    use different sequence orders; here, human,mouse,rat defines full\n\
    set and order in output alignment)\n\
\n\
    5. Extract sufficient statistics from a FASTA file.\n\
\n\
        msa_view myfile.fa --out-format SS > myfile.ss\n\
\n\
    6. Extract sufficient statistics from a MAF file for a complete\n\
    human chromosome.  (Can be used by phyloFit.)\n\
\n\
        msa_view chr1.maf --in-format MAF --out-format SS > chr1.ss\n\
\n\
    7. As in (6), but include information about regions of the\n\
    reference sequence not present in the MAF file, and include a\n\
    representation of the order in which alignment columns occur\n\
    (needed by programs such as phastCons or exoniphy).  \n\
\n\
        msa_view chr1.maf --in-format MAF --refseq chr1.fa\n\
            --out-format SS > chr1.ordered.ss\n\
\n\
    8. As in (6), but collect statistics for pairs of adjacent sites\n\
    (can be used by phyloFit to estimate a dinucleotide model).\n\
\n\
        msa_view chr1.maf --in-format MAF --out-format SS \n\
            --tuple-size 2 > chr1.pairs.ss\n\
\n\
    9. Pool sufficient statistics from several human chromosomes.\n\
\n\
        msa_view --in-format SS --aggregate human,mouse,rat \n\
            --out-format SS chr1.ss chr2.ss chr3.ss > chr123.ss\n\
\n\
    10. Extract separate sufficient statistics for the three codon\n\
    positions, as defined by annotations in a GFF file.\n\
\n\
        msa_view chr1.maf --in-format MAF --features chr22.gff \n\
            --catmap \"NCATS = 3; CDS 1-3\" --out-format SS \n\
            > chr22.pos.ss\n\
\n\
    11. As in (10), but re-orient genes on - strand so that stats\n\
    reflect + strand.  Assume genes are defined by tag \"transcript_id\".\n\
\n\
        msa_view chr1.maf --in-format MAF --features chr22.gff \n\
            --catmap \"NCATS = 3; CDS 1-3\" --reverse-groups transcript_id\n\
            --out-format SS > chr22.pos.ss\n\
\n\
OPTIONS:\n\
\n\
 (Obtaining sub-alignments and combining alignments)\n\
    --start, -s <start_col>\n\
        Starting column of sub-alignment (indexing starts with 1).\n\
        Default is 1.\n\
\n\
    --end, -e <end_col>\n\
        Ending column of sub-alignment.  Default is length of\n\
        alignment.\n\
\n\
    --seqs, -l <seq_list>\n\
        Comma-separated list of sequences to include (default)\n\
        exclude (if --exclude).  Indicate by sequence number or name\n\
        (numbering starts with 1 and is evaluated *after* --order is\n\
        applied).\n\
\n\
    --exclude, -x\n\
        Exclude rather than include specified sequences.\n\
\n\
    --refidx, -r <ref_seq>\n\
        Index of reference sequence for coordinates.  Use 0 to\n\
        indicate the coordinate system of the alignment as a whole\n\
        (this is the default).\n\
\n\
    --aggregate, -A <name_list>\n\
        (Incompatible with -i MAF) Create an aggregate alignment from\n\
        a set of alignment files, by concatenating individual\n\
        alignments.  If used with --out-format SS and --unordered-ss (and\n\
        without --start, --end, or --seqs), the aggregate\n\
        alignment will never be created explicitly (recommended for\n\
        large data sets).  The argument <name_list> must be a list of\n\
        sequence names, including all names in all specified\n\
        alignments (missing sequences will be replaced by rows of\n\
        gaps).  The standard <msa_fname> argument should be replaced\n\
        with a list of file names.\n\
\n\
 (File formats, gap stripping, reordering, etc.)\n\
    --in-format, -i PHYLIP|FASTA|MPM|MAF|SS\n\
        (Default FASTA) Input file format.  FASTA is as usual.  PHYLIP\n\
        is compatible with the formats used in the PHYLIP and PAML\n\
        packages.  MPM is the format used by the MultiPipMaker aligner\n\
        and some other of Webb Miller's older tools.  MAF (\"Multiple\n\
        Alignment Format\") is used by MULTIZ/TBA and the UCSC Genome\n\
        Browser.  SS is a simple format describing the sufficient\n\
        statistics for phylogenetic inference (distinct columns or\n\
        tuple of columns and their counts).  Use --out-format SS with\n\
        --in-format MAF for best efficiency (explicit alignment is\n\
        never created).  Also, use --unordered-ss if possible.\n\
\n\
    --out-format, -o PHYLIP|FASTA|MPM|SS\n\
        (Default FASTA)  Output file format.\n\
\n\
    --pretty, -P\n\
        Pretty-print alignment (use '.' when character matches\n\
        corresponding character in first sequence).  Ignored if\n\
        --out-format SS is selected.\n\
\n\
    --gap-strip, -G ALL|ANY|<s>\n\
        Strip columns containing all gaps, any gaps, or a gap in the\n\
        specified sequence (<s>).  Indexing starts at one and refers\n\
        to the list *after* any sequences have been added or\n\
        subtracted (via --seqs and --exclude or --order).\n\
\n\
    --collapse-missing, -p\n\
        (For use with -o SS) Convert all missing-data characters and\n\
        gaps to \"*\" characters.  Can be used to make sufficient\n\
        statistics more compact, which can improve the performance of\n\
        phyloFit (all missing data and gap characters are typically\n\
        treated the same by phyloFit anyway).\n\
\n\
    --order, -O <name_list>\n\
        Change order of rows in alignment to match sequence names\n\
        specified in name_list.  If a name appears in name_list but\n\
        not in the alignment, a row of gaps will be inserted.  This\n\
        option is applied to the alignment *before* --seqs,\n\
        --refidx, and --gap-strip are applied.\n\
\n\
    --reverse-complement, -V\n\
        Reverse complement output alignment.\n\
\n\
    --randomize, -R\n\
        Randomly permute the columns of the source alignment (done\n\
        *before* taking sub-alignment).  Requires an ordered\n\
        representation of the alignment (careful using with\n\
        --in-format SS|MAF -- will create full alignment from\n\
        sufficient statistics).\n\
\n\
    --fill-Ns, -N <s:b-e>\n\
        Fill sequence no. <s> with Ns, from <b> to <e>. Applied before\n\
        --start, --end, --seqs, --gap-strip, but after --order.\n\
        Coordinate frame depends on --refidx.  Can be used\n\
        multiple times.\n\
\n\
    --summary-only -S\n\
        Report only summary statistics, rather than complete\n\
        alignment.  Statistics are for alignment that would otherwise\n\
        be output (i.e., after other options have been applied).\n\
\n\
    --window-summary, -w <win_size>\n\
        Like -S, but output summary statistics for non-overlapping\n\
        windows of the specified size.\n\
\n\
 (Sufficient statistics)\n\
    --tuple-size, -T <tup_size>\n\
        (For use with --out-format SS).  Represent an alignment in\n\
        terms of tuples of columns of the designated size.  Useful\n\
        with context-dependent phylogenetic models.\n\
\n\
    --unordered-ss, -z\n\
        (For use with --out-format SS).  Suppress the portion of the\n\
        sufficient statistics concerned with the order in which\n\
        columns appear.  Useful for analyses for which order is\n\
        unimportant.\n\
\n\
 (MAF input)\n\
    --refseq, -M <fname>\n\
        Read the complete text of the reference sequence from\n\
        <fname> (FASTA format) and combine it with the contents of\n\
        the MAF file to produce a complete, ordered representation of\n\
        the alignment (unaligned regions will be represented by gaps).\n\
        Best used with --out-format SS.  The reference sequence of the\n\
        MAF file is assumed to be the one that appears first in each\n\
        block.\n\
\n\
    --keep-overlapping, -k\n\
        Keep blocks in MAF that have overlapping coordinates in the\n\
        reference (1st) sequence (by default, only the first one is\n\
        kept.  Useful in extracting unordered stats from a jumbled\n\
        collection of MAF blocks (e.g., output of Jim Kent's mafFrags\n\
         program).  Cannot be used with --refseq, --features, or\n\
        --cats-cycle.\n\
\n\
 (Site categories: all options require --out-format SS)\n\
    --features, -g <gff_fname>\n\
        (Requires --catmap) Read sequence annotations from the\n\
        specified file (GFF) and label the columns of the alignment\n\
        accordingly.  Note: UCSC BED and genepred formats are now\n\
        recognized as well.\n\
\n\
    --catmap, -c <fname>|<string>\n\
        (optionally use with --features) Mapping of feature types to\n\
        category numbers.  Can either give a filename or an \"inline\"\n\
        description of a simple category map, e.g., --catmap \"NCATS =\n\
        3 ; CDS 1-3\" or --catmap \"NCATS = 1 ; UTR 1\".\n\
\n\
    --cats-cycle, -Y <cycle_size>\n\
        (alternative to --features and --catmap) Assign site categories in\n\
        cycles of the specified size, e.g., as 1,2,3,...,1,2,3 (for\n\
        cycle_size == 3).  Useful for in-frame coding sequence, or to\n\
        partition a data set into nonoverlapping tuples of columns\n\
        (use with --do-cats).\n\
\n\
    --do-cats, -C <cat_list>\n\
        (For use with --features or --cats-cycle)  Obtain\n\
        sufficient statistics only for the specified categories\n\
        (comma-delimited list, by number).  Currently has no effect with\n\
        --in-format MAF (stats for all categories are automatically\n\
        collected as file is read).\n\
\n\
    --codons, -D\n\
        (For use with --features/--catmap or --cats-cycle.) Implies\n\
        --tuple-size 3 --cats-cycle 3 --do-cats 3.  Extract sufficient\n\
        statistics for in-frame codons.\n\
\n\
    --reverse-groups, -W <tag>\n\
        (For use with --features) Group features by <tag> (e.g.,\n\
        \"transcript_id\" or \"exon_id\") and reverse complement\n\
        segments of the alignment corresponding to groups on the\n\
        reverse strand.  Groups must be non-overlapping (see refeature\n\
        --unique).  Useful when extracting sufficient statistics for\n\
        strand-specific site categories (e.g., codon positions).\n\
\n\
 (Alignment cleaning)\n\
    --clean-coding, -L <seqname>\n\
        Clean an alignment of coding sequences with respect to a named\n\
        reference sequence.  Removes sites with gaps and blocks of\n\
        gapless sites smaller than 10 codons in length, ensures\n\
        everything is in-frame wrt reference sequence, prohibits\n\
        in-frame stop codons.  Reference sequence must begin with a\n\
        start codon and end with a stop codon.\n\
\n\
    --clean-indels, -I <nseqs>\n\
        Clean an alignment with special attention to indels.  Sites\n\
        with fewer than <nseqs> bases are removed; bases adjacent to\n\
        indels, and short gapless subsequences, are replaced with Ns.\n\
        If used with --tuple-size, then <tup_size>-1 columns of Ns\n\
        will be retained between columns not adjacent in the original\n\
        alignment.  Frame is not considered.\n\
\n\
 (Other)\n\
    --help, -h\n\
        Print this help message.\n\n");
}

void fill_with_Ns(MSA *msa, List *fill_N_list, msa_coord_map *map) {
  int i, j, nseq, nstart, nend;
  Regex* fill_N_re = str_re_new("([[:digit:]]+):([[:digit:]]+)-([[:digit:]]+)");
  List *word_list = lst_new_ptr(4);
  for (i = 0; i < lst_size(fill_N_list); i++) {
    String *s = lst_get_ptr(fill_N_list, i);
    if (str_re_match(s, fill_N_re, word_list, 3) <= 0 ||  
        str_as_int(lst_get_ptr(word_list, 1), &nseq) != 0 ||
        str_as_int(lst_get_ptr(word_list, 2), &nstart) != 0 ||
        str_as_int(lst_get_ptr(word_list, 3), &nend) != 0) {
      fprintf(stderr, "ERROR: cannot parse option to -N ('%s').\n", s->chars);
      exit(1);
    }
    nstart -= msa->idx_offset;
    nend -= msa->idx_offset;
    if (map != NULL) {
      nstart = msa_map_seq_to_msa(map, nstart);
      nend = msa_map_seq_to_msa(map, nend);
    }
    if (nstart < 1 || nstart > msa->length || nend < nstart || 
        nend > msa->length || nseq < 1 || nseq > msa->nseqs) {
      fprintf(stderr, "ERROR: bad sequence no. or coords. (-N option).\n");
      exit(1);
    }

    /* for now, require explicit alignment */
    if (msa->seqs == NULL) {
      fprintf(stderr, "WARNING: converting to explicit alignment (for -N).\n");
      ss_to_msa(msa);
      ss_free(msa->ss);
      msa->ss = NULL;
    }

    for (j = nstart-1; j < nend; j++) 
      msa->seqs[nseq-1][j] = 'N';

    lst_free_strings(word_list);
    str_free(s);
  }
  lst_free(word_list);
  str_re_free(fill_N_re);
}


int main(int argc, char* argv[]) {
  MSA *msa = NULL, *sub_msa = NULL;
  msa_format_type input_format = FASTA, output_format = FASTA;
  List *seqlist_str = NULL, *l = NULL;
  char *infname = NULL, *clean_seqname = NULL, *rseq_fname = NULL,
    *reverse_groups_tag = NULL;
  int i, opt_idx, startcol = 1, endcol = -1, include = 1, gap_strip_mode = NO_STRIP,
    pretty_print = FALSE, refseq = 0, tuple_size = 1, 
    ordered_stats = TRUE, cds_mode = FALSE, indel_clean_nseqs = -1, cats_done = FALSE,
    rand_perm = FALSE, reverse_compl = FALSE, stats_only = FALSE, win_size = -1, 
    cycle_size = -1, maf_keep_overlapping = FALSE, collapse_missing = FALSE;
  char c;
  List *cats_to_do = NULL, *aggregate_list = NULL, *msa_fname_list = NULL, 
    *order_list = NULL, *fill_N_list = NULL;
  msa_coord_map *map = NULL;
  GFF_Set *gff = NULL;
  CategoryMap *cm = NULL;

  struct option long_opts[] = {
    {"start", 1, 0, 's'},
    {"end", 1, 0, 'e'},
    {"refidx", 1, 0, 'r'},
    {"seqs", 1, 0, 'l'},
    {"exclude", 0, 0, 'x'},
    {"gap-strip", 1, 0, 'G'},
    {"collapse-missing", 0, 0, 'p'},
    {"in-format", 1, 0, 'i'},
    {"out-format", 1, 0, 'o'},
    {"pretty", 0, 0, 'P'},
    {"tuple-size", 1, 0, 'T'},
    {"unordered-ss", 0, 0, 'z'},
    {"features", 1, 0, 'g'},
    {"catmap", 1, 0, 'c'},
    {"cats-cycle", 1, 0, 'Y'},
    {"do-cats", 1, 0, 'C'},
    {"codons", 0, 0, 'D'},
    {"aggregate", 1, 0, 'A'},
    {"refseq", 1, 0, 'M'},
    {"order", 1, 0, 'O'},
    {"summary-only", 0, 0, 'S'},
    {"window-summary", 1, 0, 'w'},
    {"fill-Ns", 1, 0, 'N'},
    {"reverse-complement", 0, 0, 'V'},
    {"reverse-groups", 1, 0, 'W'},
    {"clean-coding", 1, 0, 'L'},
    {"clean-indels", 1, 0, 'I'},
    {"randomize", 0, 0, 'R'},
    {"fill-Ns", 1, 0, 'N'},
    {"keep-overlapping", 0, 0, 'k'},
    {"help", 0, 0, 'h'}
  };

  while ((c = getopt_long(argc, argv, "m:i:o:s:e:l:G:r:T:a:g:c:C:L:I:A:M:O:w:N:Y:DVxPzRSkh", long_opts, &opt_idx)) != -1) {
    switch(c) {
    case 'm':
      infname = optarg;
      break;
    case 'i':
      input_format = msa_str_to_format(optarg);
      if (input_format == -1) die("ERROR: bad input format.  Try 'msa_view -h' for help.\n");
      break;
    case 's':
      startcol = atoi(optarg);
      break;
    case 'e':
      endcol = atoi(optarg);
      break;
    case 'l':
      seqlist_str = get_arg_list(optarg);
      break;
    case 'x':
      include = FALSE;
      break;
    case 'G':
      if (!strcmp(optarg, "ALL")) gap_strip_mode = STRIP_ALL_GAPS;
      else if (!strcmp(optarg, "ANY")) gap_strip_mode = STRIP_ANY_GAPS;
      else gap_strip_mode = atoi(optarg);
      break;
    case 'p':
      collapse_missing = TRUE;
      break;
    case 'o':
      output_format = msa_str_to_format(optarg);
      if (output_format == -1) die("ERROR: bad output format.  Try 'msa_view -h' for help.\n");
      break;
    case 'r':
      refseq = atoi(optarg);
      break;
    case 'T':
      tuple_size = atoi(optarg);
      break;
    case 'A':
      aggregate_list = get_arg_list(optarg);
      break;
    case 'g':
      gff = gff_read_set(fopen_fname(optarg, "r"));
      break;
    case 'c':
      cm = cm_new_string_or_file(optarg);
      break;
    case 'P':
      pretty_print = TRUE;
      break;
    case 'Y':
      cycle_size = atoi(optarg);
      break;
    case 'D':
      cycle_size = tuple_size = 3;
      if (cats_to_do == NULL) cats_to_do = lst_new_int(1);
      lst_push_int(cats_to_do, 3);
      break;
    case 'z':
      ordered_stats = FALSE;
      break;
    case 'L':
      clean_seqname = optarg;
      break;
    case 'O': 
      order_list = get_arg_list(optarg);
      break;
    case 'I':
      indel_clean_nseqs = atoi(optarg);
      break;
    case 'R':
      rand_perm = TRUE;
      break;
    case 'M':
      rseq_fname = optarg;
      break;
    case 'V':
      reverse_compl = TRUE;
      break;
    case 'W':
      reverse_groups_tag = optarg;
      break;
    case 'S':
      stats_only = TRUE;
      break;
    case 'w':
      stats_only = TRUE;
      win_size = atoi(optarg);
      break;
    case 'N':
      if (fill_N_list == NULL) fill_N_list = lst_new_ptr(10);
      lst_push_ptr(fill_N_list, str_new_charstr(optarg));
      break;
    case 'h':
      print_usage();
      exit(0);
    case 'C':
      cats_to_do = get_arg_list_int(optarg);
      break;
    case 'k':
      maf_keep_overlapping = TRUE;
      break;
    case '?':
      fprintf(stderr, "Bad argument.  Try 'msa_view -h' for help.\n");
      exit(1); 
    }
  }

  if (optind >= argc) {
    fprintf(stderr, "Missing alignment filename.  Try 'msa_view -h' for help.\n");
    exit(1);
  }
  else infname = argv[optind];

  if (aggregate_list != NULL) {
    msa_fname_list = get_arg_list(infname);

    if (input_format == SS && output_format == SS && ordered_stats == 1) {
      fprintf(stderr, "WARNING: ignoring request for ordered sufficient statistics (not supported with --aggregate).\n");
      ordered_stats = 0;
    }

    if (output_format == SS && ordered_stats == 0 && 
        startcol == 1 && endcol == -1 && seqlist_str == NULL) {
      msa = ss_aggregate_from_files(msa_fname_list, input_format, 
                                    aggregate_list, NULL, tuple_size, 
                                    cats_to_do, cycle_size);
                                /* avoid creating aggregate alignment
                                   explicitly, if possible */
      cats_done = 1;            /* in this case, cats are taken care of */
    }
    else 
      msa = msa_concat_from_files(msa_fname_list, input_format, 
                                  aggregate_list, NULL);
  }

  else if (input_format == MAF) {
    FILE *RSEQF = NULL;

    if (rseq_fname != NULL) RSEQF = fopen_fname(rseq_fname, "r");

    if (output_format == SS && RSEQF == NULL && ordered_stats == 1 && 
        gff == NULL && startcol == 1 && endcol == -1)
      ordered_stats = 0;        /* in this case, assume unordered
                                   stats are desired; can't think of
                                   any value in collecting ordered
                                   stats, and it's a common mistake to
                                   forget -z */

    msa = maf_read(fopen_fname(infname, "r"), RSEQF, tuple_size, 
                   gff, cm, cycle_size, 
                   output_format != SS || ordered_stats, 
                   reverse_groups_tag, gap_strip_mode, maf_keep_overlapping);
                                /* store order unless output is SS and
                                   no ordered stats */
  }

  else {
    msa = msa_new_from_file(fopen_fname(infname, "r"), input_format, NULL);
    if (msa == NULL) die ("ERROR reading %s.\n", infname);
  }

  if (order_list != NULL)
    msa_reorder_rows(msa, order_list);

  if (seqlist_str != NULL)
    l = msa_seq_indices(msa, seqlist_str);

  if (rand_perm) msa_permute(msa);

  if (startcol < 1 || (endcol != -1 && endcol > msa->length) || 
      (endcol != -1 && endcol < startcol)) 
                                /* careful: msa->length is unsigned */
    die("ERROR: must have 1 <= start <= end <= [msa_length]\n");


  if (refseq < 0 || refseq > msa->nseqs) { 
    fprintf(stderr, "ERROR: reference sequence out of range.\n");
    exit(1); 
  }
  if (refseq != 0) {
    if ((input_format == SS || aggregate_list != NULL) && 
        msa->ss->tuple_idx == NULL) {
      fprintf(stderr, "ERROR: an ordered representation of the alignment columns is required.\n");
      exit(1);
    }
    map = msa_build_coord_map(msa, refseq);
    startcol = msa_map_seq_to_msa(map, startcol);
    if (endcol != -1) endcol = msa_map_seq_to_msa(map, endcol);
  }

  if (endcol == -1) endcol = msa->length;

  /* fill with Ns, if necessary */
  if (fill_N_list != NULL) fill_with_Ns(msa, fill_N_list, map);

  if (gff != NULL && cm == NULL) 
    cm = cm_new_from_features(gff);

  /* read annotations and label columns, if necessary */
  if (gff != NULL && input_format != MAF) {
    if (input_format == SS || input_format == MAF || aggregate_list != NULL) {
      if (msa->ss->tuple_idx == NULL) {
        fprintf(stderr, "ERROR: ordered representation of alignment required with --features.\n");
        exit(1);
      }
    }

    /* convert GFF to coordinate frame of alignment */
    msa_map_gff_coords(msa, gff, -1, 0, 0, NULL);

    if (reverse_groups_tag != NULL) { /* reverse complement by group */
      if (input_format == SS) {
        fprintf(stderr, "ERROR: need an explicit representation of the alignment to reverse complement.\n");
        exit(1);
      }

      gff_group(gff, reverse_groups_tag);
      msa_reverse_compl_feats(msa, gff, NULL);
    }

    /* label categories */
    msa_label_categories(msa, gff, cm);
  }
  else if (cycle_size > 0 && clean_seqname == NULL && !cats_done) {
    msa->categories = (int*)smalloc(msa->length * sizeof(int));
    msa->ncats = cycle_size;
    for (i = 0; i < msa->length; i++)
      msa->categories[i] = (i % cycle_size) + 1;
  }  

  if (startcol > 1 || endcol < msa->length || l != NULL)
    sub_msa = msa_sub_alignment(msa, l, include, startcol - 1, endcol);
  else
    sub_msa = msa;

  /* if we're working with transformed coords, we'll mark the offset
     using the same coord frame that was used to specify the
     startcol (usually produces the desired result) */
  if (map != NULL) 
    sub_msa->idx_offset = msa_map_msa_to_seq(map, startcol) - 1 + 
      msa->idx_offset;

  if (clean_seqname != NULL) {
    String *errstr = str_new(STR_MED_LEN);
    int i;

    assert(startcol == 1);

    /* find seq number corresponding to name */
    for (i = 0; i < sub_msa->nseqs && 
           strcmp(clean_seqname, sub_msa->names[i]) != 0; i++);
    if (i == msa->nseqs) { 
      fprintf(stderr, "ERROR: no match for argument of --clean-coding.\n");
      exit(1);
    }
    if (msa_coding_clean(sub_msa, i, MIN_NCODONS, errstr) != 0) {
      fprintf(stderr, "ERROR: unable to clean coding alignment (%s).\n%s\n", 
              infname, errstr->chars);
      exit(1);
    }

    /* in this case, assign categories *after* cleaning */
    if (cds_mode && gff == NULL) {
      sub_msa->categories = (int*)smalloc(sub_msa->length * sizeof(int));
      sub_msa->ncats = 3;
      for (i = 0; i < sub_msa->length; i++)
        sub_msa->categories[i] = (i % 3) + 1;
    }

    str_free(errstr);
  }
  else if (indel_clean_nseqs >= 0) {
    msa_indel_clean(sub_msa, INDEL_BORDER, MIN_NBASES, indel_clean_nseqs,
                    tuple_size, 'N');
  }
  else if (gap_strip_mode != NO_STRIP && input_format != MAF)
    msa_strip_gaps(sub_msa, gap_strip_mode);

  /* create sufficient stats, if necessary */
  if (output_format == SS) {
    if (sub_msa->ss == NULL)
      ss_from_msas(sub_msa, tuple_size, ordered_stats, cats_to_do, NULL, NULL, -1);
    else {
      if (sub_msa->ss->tuple_size != tuple_size) {
        fprintf(stderr, "ERROR: tuple size in SS file does not match desired tuple size for output.\nConversion not supported.\n");
        exit(1);
      }
      if (sub_msa->ss->tuple_idx != NULL && ordered_stats == 0) {
        free(sub_msa->ss->tuple_idx);
        sub_msa->ss->tuple_idx = NULL;
      }
    }

    if (collapse_missing) ss_collapse_missing(sub_msa, TRUE);
  }

  if (gff == NULL && reverse_compl) {
    if (msa->seqs == NULL && msa->ss->tuple_idx == NULL) {
      fprintf(stderr, "ERROR: an ordered representation of the alignment is required to\nreverse complement.\n");
      exit(1);
    }
    msa_reverse_compl(sub_msa);
  }

  if (stats_only) {             /* only print summary stats */
    msa_print_stats(sub_msa, stdout, NULL, 1, -1, -1);
    if (win_size > 0) {         /* print summary for each window */
      for (i = 0; i < sub_msa->length; i += win_size) {
        char name[STR_SHORT_LEN];
        sprintf(name, "%d-%d", i, i+win_size-1);
        msa_print_stats(sub_msa, stdout, name, 0, i, 
                        min(i+win_size, msa->length));
      }
    }
    else 
      msa_print_stats(sub_msa, stdout, infname, 0, -1, -1);
  }

  else                          /* print alignment */
    msa_print(stdout, sub_msa, output_format, pretty_print);

  if (sub_msa != msa) msa_free(sub_msa);
  msa_free(msa);
  if (map != NULL) msa_map_free(map);
  if (l != NULL) lst_free(l);

  return 0;
}
