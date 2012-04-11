/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <misc.h>
#include <gff.h>
#include <bed.h>
#include <genepred.h>
#include <hashtable.h>
#include <wig.h>

/* to do: add an option to insert features for splice sites or
   start/stop coords at exon boundaries ('addsignals'); */

void usage(char *prog) {
  printf("\n\
PROGRAM:        %s\n\
\n\
DESCRIPTION:    Read a file representing a set of features, optionally\n\
                alter the set in one or more of several possible ways, then\n\
                output it in the desired format.  Input and output formats\n\
                may be GFF, BED, or genepred.\n\
\n\
                The input format is recognized automatically, but auto-\n\
                recognition requires a 'seekable' input stream (e.g., an\n\
                actual file rather than a pipe from stdin).\n\
\n\
USAGE:          %s [OPTIONS] <infile>\n\
\n\
OPTIONS:\n\
    --include-only, -i <types>\n\
        Include only features of the specified types (comma-delimited list);\n\
        filter out everything else.\n\
\n\
    --include-groups, -l <file>\n\
        Include only groups whose names are listed in the specified file.\n\
        Group names in file must be delimited by white-space (can be on\n\
        any number of lines).\n\
\n\
    --sort, -s\n\
        Sort features primarily by start position and secondarily\n\
        by end position (usually has desired effect in case of short\n\
        overlapping features, e.g., start & stop codons).  Features\n\
        will be sorted both across groups and within groups, but members\n\
        of a group will be kept together.\n\
\n\
    --unique, -u\n\
        Ensures that output contains no overlapping groups (or\n\
        subgroups, if -e).  If groups overlap, the one with the highest\n\
        score (if available) or longest length (if no score) is kept and\n\
        others are discarded.  Warning: long UTRs can have undesirable\n\
        results; filter out UTR exons to avoid.\n\
\n\
    --groupby, -g <tag>\n\
        Group features according to specified tag (default \"transcript_id\")\n\
\n\
    --exongroup, -e <tag>\n\
        Sub-group features into contiguous sets, and define\n\
        sub-groups using specified tag (e.g., \"exon_id\").  Can be\n\
        used to group the features describing individual exons, e.g.,\n\
        each CDS and its flanking splice sites.  Only features in the\n\
        same major group will be included in the same minor group\n\
        (e.g., exons of the same transcript).\n\
\n\
    --fix-start-stop, -f\n\
        Ensure that CDS features include start codons and exclude stop\n\
        codons, as required by the GTF2 standard.  Assumes at most one\n\
        start_codon and at most one stop_codon per group.\n\
\n\
    --add-utrs, -U\n\
        Create UTR features for portions of exons outside CDS (only\n\
        useful with GFF output; features must be grouped at level\n\
        of transcript).\n\
\n\
    --add-introns, -I\n\
        Create intron features between exons (only useful with GFF output;\n\
        features must be grouped at level of transcript).\n\
\n\
    --add-signals, -S\n\
        Adds features for start and stop codons and 3' and 5' splice\n\
        sites (only useful with GFF output; features must be grouped\n\
        at level of transcript).  Start and stop codons will be added\n\
        as required by the GTF2 standard (--fix-start-stop is not\n\
        necessary).  Warning: does not correctly handle case of splice\n\
        site in middle of start or stop codon.\n\
\n\
    --output, -o gff|bed|genepred|wig\n\
        Output format (default gff).  Note that wig output is fixedStep\n\
        can only be used if all elements have a score and are of equal\n\
        length.\n\
\n\
    --simplebed, -b\n\
        (for use with --output bed) Create a separate line for each\n\
        feature in bed output (by default, all features of a group are\n\
        described by a single line).\n\
\n\
    --discards, -d <fname>\n\
        Write any discarded features to specified file.\n\
\n\
    --help, -h\n\
        Print this help message.\n\n", prog, prog);
  exit(0);
}

int main(int argc, char *argv[]) {
  char c;
  int opt_idx;
  GFF_Set *gff;
  List *include = NULL;
  char *groupby = "transcript_id", *exongroup_tag = NULL;
  int unique = FALSE, sort = FALSE, simplebed = FALSE, fix_start_stop = FALSE,
    add_utrs = FALSE, add_introns = FALSE, add_signals = FALSE;
  enum {GFF, BED, GENEPRED, WIG} output_format = GFF;
  FILE *discards_f = NULL, *groups_f = NULL;

  struct option long_opts[] = {
    {"output", 1, 0, 'o'},
    {"include-only", 1, 0, 'i'},
    {"include-groups", 1, 0, 'l'},
    {"groupby", 1, 0, 'g'},
    {"exongroup", 1, 0, 'e'},
    {"add-utrs", 0, 0, 'U'},
    {"add-introns", 0, 0, 'I'},
    {"add-signals", 0, 0, 'S'},
    {"fix-start-stop", 0, 0, 'f'},
    {"unique", 0, 0, 'u'},
    {"sort", 0, 0, 's'},
    {"simplebed", 0, 0, 'b'},
    {"discards", 1, 0, 'd'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "o:i:l:g:e:d:UISfusbh", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'o':
      if (!strcmp("bed", optarg)) output_format = BED;
      else if (!strcmp("genepred", optarg)) output_format = GENEPRED;
      else if (!strcmp("wig", optarg)) output_format = WIG;
      else if (strcmp("gff", optarg)) die("ERROR: bad output format.\n");
      break;
    case 'i':
      include = get_arg_list(optarg);
      break;
    case 'l':
      groups_f = phast_fopen(optarg, "r");
      break;
    case 'g':
      groupby = optarg;
      break;
    case 'e':
      exongroup_tag = optarg;
      break;
    case 'U':
      add_utrs = TRUE;
      break;
    case 'I':
      add_introns = TRUE;
      break;
    case 'S':
      add_signals = TRUE;
      break;
    case 'f':
      fix_start_stop = TRUE;
      break;
    case 'u':
      unique = TRUE;
      break;
    case 'b':
      simplebed = TRUE;
      output_format = BED;
      break;
    case 'd':
      discards_f = phast_fopen(optarg, "w+");
      break;
    case 's':
      sort = TRUE;
      break;
    case 'h':
      usage(argv[0]);
    case '?':
      die("Bad argument.  Try '%s -h'.\n", argv[0]);
    }
  }

  if (optind != argc - 1) 
    die("Input filename required.  Try '%s -h'.\n", argv[0]);

  set_seed(-1);

  gff = gff_read_set(phast_fopen(argv[optind], "r"));

  if (lst_size(gff->features) == 0) exit(0); /* helps avoid unexpected
                                                behavior below */

  /* filter by type */
  if (include != NULL) gff_filter_by_type(gff, include, FALSE, discards_f);

  /* group */
  gff_group(gff, groupby);

  /* utrs, introns, & signals */
  if (add_utrs) gff_create_utrs(gff);
  if (add_introns) gff_create_introns(gff);
  if (add_signals) gff_create_signals(gff);

  /* subgroup */
  if (exongroup_tag != NULL) gff_exon_group(gff, exongroup_tag);

  /* filter by group */
  if (groups_f != NULL) {
    String *s = str_new(STR_LONG_LEN);
    List *groups = lst_new_ptr(10000);
    str_slurp(s, groups_f);
    str_split(s, NULL, groups);
    gff_filter_by_group(gff, groups);
    lst_free_strings(groups); lst_free(groups);
    str_free(s);
  }

  /* sort */
  if (sort) gff_sort(gff);

  /* make unique */
  if (unique) gff_remove_overlaps(gff, discards_f);
  
  if (fix_start_stop) gff_fix_start_stop(gff);

  if (output_format == BED)
    gff_print_bed(stdout, gff, !simplebed);
  else if (output_format == GENEPRED)
    gff_print_genepred(stdout, gff);
  else if (output_format == WIG)
    wig_print(stdout, gff);
  else 
    gff_print_set(stdout, gff);
  gff_free_set(gff);
  
  return 0;
}
