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
                output it in the desired format.  Currently, input and\n\
                output formats may be GFF or BED.  Input format is\n\
                automatically recognized.  The full, 12-column BED format \n\
                is used for output, but abbreviated versions are accepted \n\
                as input.\n\
\n\
USAGE:          %s [OPTIONS] <infile>\n\
\n\
OPTIONS:\n\
    --output, -o GFF|BED\n\
        Output format (default GFF).\n\
\n\
    --include, -i <types>\n\
        Include only features of the specified types (comma-delimited list).\n\
\n\
    --groupby, -g <tag>\n\
        Group features according to specified tag (e.g., \"transcript_id\")\n\
\n\
    --exongroup, -e <tag>\n\
        Sub-group features into contiguous sets, and define\n\
        sub-groups with specified tag (e.g., \"exon_id\").  Can be\n\
        used to group the features describing individual exons, for\n\
        example, each CDS and its flanking splice sites.  Only\n\
        features in the same major group will be included in the same\n\
        minor group (e.g., exons of the same transcript).\n\
\n\
    --unique, -u\n\
        Ensures that output contains no overlapping groups or\n\
        subgroups (if -e).  If groups overlap, the first one encountered\n\
        is retained and subsequent ones are discarded.\n\
\n\
    --sort, -s\n\
        Sort features by start position.  If features are grouped,\n\
        they will be sorted within groups, and groups will be sorted\n\
        by first start position.\n\
\n\
    --simplebed, -b\n\
        (for use with --output BED) Create a separate line for each\n\
        feature in BED output (by default, all features of a group are\n\
        described by a single line).\n\
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
  String *groupby = NULL, *exongroup_tag = NULL;
  int unique = 0, sort = 0, simplebed = 0;
  enum {GFF, BED} output_format = GFF;

  struct option long_opts[] = {
    {"output", 1, 0, 'o'},
    {"include", 1, 0, 'i'},
    {"groupby", 1, 0, 'g'},
    {"exongroup", 1, 0, 'e'},
    {"unique", 0, 0, 'u'},
    {"sort", 0, 0, 's'},
    {"simplebed", 0, 0, 'b'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "o:i:g:e:usbh", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'o':
      if (!strcmp("BED", optarg)) output_format = BED;
      else if (strcmp("GFF", optarg)) die("ERROR: bad output format.\n");
      break;
    case 'i':
      include = get_arg_list(optarg);
      break;
    case 'g':
      groupby = str_new_charstr(optarg);
      break;
    case 'e':
      exongroup_tag = str_new_charstr(optarg);
      break;
    case 'u':
      unique = 1;
      break;
    case 'b':
      simplebed = 1;
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
  if (include != NULL) gff_filter_by_type(gff, include);

  /* group */
  if (groupby != NULL) gff_group(gff, groupby);

  /* subgroup */
  if (exongroup_tag != NULL) gff_exon_group(gff, exongroup_tag);

  /* make unique */
  if (unique) {
    die("Sorry, --unique not yet implemented...\n");
  }

  /* sort */
  if (sort) gff_sort(gff);

  if (output_format == BED)
    gff_print_bed(stdout, gff, simplebed ? NULL : groupby, NULL);
  else 
    gff_print_set(stdout, gff);
  
  return 0;
}
