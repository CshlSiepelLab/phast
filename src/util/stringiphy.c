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

typedef enum {INITIAL, INTERNAL, TERMINAL, SINGLETON} ExonType;

/* meta-data for feature group representing single exon prediction */
typedef struct {
  int start;
  int end;
  char strand;
  String *seqname;
  int start_phase;
  int end_phase;
  ExonType type;
} GroupData;

/* threshold for joining exons, in base pairs (don't join exons if
   farther apart than JOIN_THRESHOLD bases, even if strand, phase,
   etc. are okay) */
#define JOIN_THRESHOLD 20000

void usage(char *prog) {
  printf("\n\
PROGRAM:      %s\n\
\n\
DESCRIPTION:  Attempt to string together predicted exons into full transcripts\n\
              or gene fragments.\n\
\n\
USAGE:        %s [OPTIONS] <exons.gff>\n\
              (Predictions should be sorted and non-overlapping; use\n\
              'refeature --sort --unique')\n\
\n\
OPTIONS:\n\
    --help, -h\n\
        Print this help message.\n\n", prog, prog);
  exit(0);
}

/* redefine group for set of features, by replacing 'attribute' field
   of all features with new tag and value */
void regroup(GFF_FeatureGroup *g, char *gene_tag, char *gene_val) {
  int i;
  char tmpstr[STR_LONG_LEN];
  for (i = 0; i < lst_size(g->features); i++) {
    GFF_Feature *f = lst_get_ptr(g->features, i);
    sprintf(tmpstr, "%s \"%s\"", gene_tag, gene_val);
    str_cpy_charstr(f->attribute, tmpstr);
  }
}

int main(int argc, char *argv[]) {
  GFF_Set *exons;
  char c;
  int i, j, opt_idx;
  List *group_data;

  struct option long_opts[] = {
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "h", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'h':
      usage(argv[0]);
    case '?':
      die("Bad argument.  Try '%s -h'.\n", argv[0]);
    }
  }

  if (optind != argc - 1) 
    die("Input filename required.  Try '%s -h'.\n", argv[0]);

  set_seed(-1);
    
  exons = gff_read_set(phast_fopen(argv[optind], "r"));

  gff_group(exons, "transcript_id");
  gff_sort(exons);
  group_data = lst_new_ptr(lst_size(exons->groups));
  /* first get end coord, strand, and start/end phase of each feature
     group (exon) and classify as initial, internal, or terminal */
  for (i = 0; i < lst_size(exons->groups); i++) {
    GFF_FeatureGroup *g = lst_get_ptr(exons->groups, i);
    GroupData *d = smalloc(sizeof(GroupData));
    d->start = g->start;
    d->end = -1;
    d->strand = '.';
    d->seqname = NULL;
    d->start_phase = d->end_phase = -1;
    d->type = INTERNAL;
    for (j = 0; j < lst_size(g->features); j++) {
      GFF_Feature *f = lst_get_ptr(g->features, j);

      if (f->strand == '.') 
        die("ERROR: missing strand in group '%s'\n", g->name->chars);
      if (j == 0) d->strand = f->strand;
      else if (d->strand != f->strand)
        die("ERROR: inconsistent strand in group '%s'\n", g->name->chars);

      if (j == 0) d->seqname = f->seqname;
      else if (!str_equals(d->seqname, f->seqname))
        die("ERROR: inconsistent seqname in group '%s'\n", g->name->chars);      

      if (f->end > d->end) d->end = f->end;

      if (str_equals_charstr(f->feature, GFF_CDS_TYPE)) { 
        if (d->start_phase != -1 || d->end_phase != -1) /* should only be one */
          die("ERROR: multiple CDS features in group '%s'\n", g->name->chars);
        d->start_phase = f->frame;
        d->end_phase = (f->frame + f->end - f->start + 1) % 3;
      }
      if (str_equals_charstr(f->feature, GFF_START_TYPE)) {
        if (d->type == TERMINAL) d->type = SINGLETON;
        else d->type = INITIAL; 
      }
      else if (str_equals_charstr(f->feature, GFF_STOP_TYPE)) {
        if (d->type == INITIAL) d->type = SINGLETON;
        else d->type = TERMINAL;
      }
    }
    lst_push_ptr(group_data, d);
  }

  /* now join groups */
  for (i = 0; i < lst_size(exons->groups); ) {
    GFF_FeatureGroup *g1 = lst_get_ptr(exons->groups, i);
    GroupData *d1 = lst_get_ptr(group_data, i);
    String *groupname = g1->name;
    regroup(g1, "transcript_id", groupname->chars);

    for (j = i+1; j < lst_size(exons->groups); j++) {
      GFF_FeatureGroup *g2 = lst_get_ptr(exons->groups, j);
      GroupData *d2 = lst_get_ptr(group_data, j);

      if (str_equals(d1->seqname, d2->seqname) &&
          d1->strand == d2->strand && d2->start < d1->end + JOIN_THRESHOLD &&
          ((d1->strand == '+' && d1->end_phase == d2->start_phase &&
            (d1->type == INITIAL || d1->type == INTERNAL) && 
            (d2->type == INTERNAL || d2->type == TERMINAL)) || 
           (d1->strand == '-' && d2->end_phase == d1->start_phase &&
            (d2->type == INITIAL || d2->type == INTERNAL) && 
            (d1->type == INTERNAL || d1->type == TERMINAL)))) {

        regroup(g2, "transcript_id", groupname->chars);
        g1 = g2;
        d1 = d2;
      }

      else break;
    }
    i = j;
  }

  gff_group(exons, "transcript_id");
  gff_add_gene_id(exons);
  gff_print_set(stdout, exons);
  return 0;
}


