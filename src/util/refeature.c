#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <misc.h>
#include <gff.h>
#include <bed.h>
#include <genepred.h>
#include <hashtable.h>

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
    --output, -o gff|bed|genepred\n\
        Output format (default gff).\n\
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
  int unique = 0, sort = 0, simplebed = 0;
  enum {GFF, BED, GENEPRED} output_format = GFF;
  FILE *discards_f = NULL, *groups_f = NULL;

  struct option long_opts[] = {
    {"output", 1, 0, 'o'},
    {"include-only", 1, 0, 'i'},
    {"include-groups", 1, 0, 'l'},
    {"groupby", 1, 0, 'g'},
    {"exongroup", 1, 0, 'e'},
    {"unique", 0, 0, 'u'},
    {"sort", 0, 0, 's'},
    {"simplebed", 0, 0, 'b'},
    {"discards", 1, 0, 'd'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "o:i:l:g:e:d:usbh", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'o':
      if (!strcmp("bed", optarg)) output_format = BED;
      else if (!strcmp("genepred", optarg)) output_format = GENEPRED;
      else if (strcmp("gff", optarg)) die("ERROR: bad output format.\n");
      break;
    case 'i':
      include = get_arg_list(optarg);
      break;
    case 'l':
      groups_f = fopen_fname(optarg, "r");
      break;
    case 'g':
      groupby = optarg;
      break;
    case 'e':
      exongroup_tag = optarg;
      break;
    case 'u':
      unique = 1;
      break;
    case 'b':
      simplebed = 1;
      output_format = BED;
      break;
    case 'd':
      discards_f = fopen_fname(optarg, "w+");
      break;
    case 's':
      sort = 1;
      break;
    case 'h':
      usage(argv[0]);
    case '?':
      die("Bad argument.  Try '%s -h'.\n", argv[0]);
    }
  }

  if (optind != argc - 1) 
    die("Input filename required.  Try '%s -h'.\n", argv[0]);

  gff = gff_read_set(fopen_fname(argv[optind], "r"));

  if (lst_size(gff->features) == 0) exit(0); /* helps avoid unexpected
                                                behavior below */

  /* filter by type */
  if (include != NULL) gff_filter_by_type(gff, include, FALSE, discards_f);

  /* group */
  gff_group(gff, groupby);

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

  if (output_format == BED)
    gff_print_bed(stdout, gff, !simplebed);
  else if (output_format == GENEPRED)
    gff_print_genepred(stdout, gff);
  else 
    gff_print_set(stdout, gff);
  
  return 0;
}
