#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <misc.h>
#include <gff.h>
#include <bed.h>
#include <hashtable.h>

/* to do: support reading and writing of genepred format (should be
   easy -- very similar to bed 12); add an option to insert features
   for splice sites or start/stop coords at exon boundaries
   ('addsignals'); make 'unique' option give priority to highest
   scoring group (or longest, if no score is available) */

void usage(char *prog) {
  printf("\n\
PROGRAM:        %s\n\
\n\
DESCRIPTION:    Read a file representing a set of features, optionally\n\
                alter the set in one or more of several possible ways, then\n\
                output it in the desired format.  Input and output formats\n\
                may be GFF, BED, or genepred.  Input format is automatically\n\
                recognized.  The full, 12-column BED format is used for\n\
                output, but abbreviated versions (e.g., 3-column, 4-column,\n\
                or 6-column) are accepted as input.\n\
\n\
USAGE:          %s [OPTIONS] <infile>\n\
\n\
OPTIONS:\n\
    --include, -i <types>\n\
        Include only features of the specified types (comma-delimited list).\n\
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
    --output, -o GFF|BED|genepred\n\
        Output format (default GFF).\n\
\n\
    --simplebed, -b\n\
        (for use with --output BED) Create a separate line for each\n\
        feature in BED output (by default, all features of a group are\n\
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
  enum {GFF, BED, genepred} output_format = GFF;
  FILE *discards_f = NULL;

  struct option long_opts[] = {
    {"output", 1, 0, 'o'},
    {"include", 1, 0, 'i'},
    {"groupby", 1, 0, 'g'},
    {"exongroup", 1, 0, 'e'},
    {"unique", 0, 0, 'u'},
    {"sort", 0, 0, 's'},
    {"simplebed", 0, 0, 'b'},
    {"discards", 1, 0, 'd'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "o:i:g:e:d:usbh", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'o':
      if (!strcmp("BED", optarg)) output_format = BED;
      else if (strcmp("GFF", optarg)) die("ERROR: bad output format.\n");
      break;
    case 'i':
      include = get_arg_list(optarg);
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

  /* filter */
  if (include != NULL) gff_filter_by_type(gff, include, discards_f);

  /* group */
  gff_group(gff, groupby);

  /* subgroup */
  if (exongroup_tag != NULL) gff_exon_group(gff, exongroup_tag);

  /* sort */
  if (sort) gff_sort(gff);

  /* make unique */
  if (unique) gff_remove_overlaps(gff, discards_f);

  if (output_format == BED)
    gff_print_bed(stdout, gff, simplebed ? NULL : groupby, NULL);
  else if (output_format == genepred)
    die("Sorry, no genepred output yet. :(\n");
  else 
    gff_print_set(stdout, gff);
  
  return 0;
}
