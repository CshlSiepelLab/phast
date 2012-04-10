/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: hmm_view.c,v 1.7 2008-11-12 02:07:58 acs Exp $ */

#include <stdlib.h>
#include <stdio.h>
#include "hmm.h"
#include <getopt.h>
#include "category_map.h"
#include "gap_patterns.h"

void print_usage() {
    printf("\n\
PROGRAM: hmm_view\n\
\n\
DESCRIPTION: produces a graphical description of the state-transition\n\
structure of a phylo-HMM, which can be converted to a viewable image\n\
using the 'dot' program. (http://www.graphviz.org)\n\
\n\
USAGE: hmm_view [OPTIONS] <hmm_fname> <cat_map_fname>|<cat_map_str>\n\
\n\
OPTIONS:\n\
    -k <nrcats>   Assume a separate version of each state for each of \n\
                  <nrcats> rate categories. \n\
    -i <icats>    Assume use of indel model for specified category names.\n\
    -t <tree>     (Required with -i) Tree topology to assume for indel\n\
                  model (.nh file).\n\
    -C <cats>     Show only the states corresponding to the specified\n\
                  category names.\n\
    -R <piv>      Reflect the HMM about the specified 'pivot' categories.\n\
                  (Not yet implemented.)\n\
    -x            Don't show unconnected states.\n\n");
}

int main(int argc, char *argv[]) {
  HMM *hmm;
  CategoryMap *cm;
  List *indel_cats = NULL, *cats_to_show = NULL, *pivots = NULL;
  GapPatternMap *gpm = NULL;
  int i, j, nratecats = 1, gapped_cat, basecat, suppress_unconnected = 0;
  int *new_to_old=NULL, *show_cat;
  TreeNode *tree = NULL;
  double t;
  char c;
  String *source, *sink;

  while ((c = getopt(argc, argv, "k:i:t:C:xh")) != -1) {
    switch(c) {
    case 'k':
      nratecats = get_arg_int(optarg);
      break;
    case 'i':
      indel_cats = get_arg_list(optarg);
      break;
    case 't':
      tree = tr_new_from_file(phast_fopen(optarg, "r"));
      break;
    case 'C':
      cats_to_show = get_arg_list(optarg);
      break;
    case 'R':
      pivots = get_arg_list(optarg);
      break;
    case 'x':
      suppress_unconnected = 1;
      break;
    case 'h':
      print_usage();
      exit(0);
    case '?':
      die("Bad argument.  Try 'hmm_view -h' for help.\n");
    }
  }

  if (optind != argc - 2) 
    die("Bad arguments.  Try 'hmm_view -h' for help.\n");

  if (indel_cats != NULL && tree == NULL) 
    die("Must specify -t with -i.  Try 'hmm_view -h' for help.\n");

  set_seed(-1);

  hmm = hmm_new_from_file(phast_fopen(argv[optind], "r"));
  cm = cm_new_string_or_file(argv[optind+1]);

  show_cat = smalloc((cm->ncats+1) * sizeof(int));
  if (cats_to_show != NULL) {
    List *l = cm_get_category_list(cm, cats_to_show, 0);
    for (i = 0; i <= cm->ncats; i++) show_cat[i] = 0;
    for (i = 0; i < lst_size(l); i++)
      show_cat[lst_get_int(l, i)] = 1;
    lst_free(l);
  }
  else 
    for (i = 0; i <= cm->ncats; i++) show_cat[i] = 1;

  if (indel_cats != NULL)
    gpm = gp_create_gapcats(cm, indel_cats, tree, FALSE);  

  if (hmm->nstates != (cm->unspooler == NULL ? cm->ncats + 1 : 
                       cm->unspooler->nstates_unspooled) * nratecats) {
    die("ERROR: number of states in HMM must equal number of site categories (unspooled).\n");
  }

  if (pivots != NULL) 
    hmm = hmm_reverse_compl(hmm, pivots, new_to_old);

  printf("\n\
digraph hmm {\n\
        rankdir=LR;\n\
        size=\"10,7.5\";\n\
        ratio=\"compress\";\n\
        orientation=land;\n\
        node [shape = box];\n\
");

  for (i = 0; i < hmm->begin_transitions->size; i++) {
    
    gapped_cat = cm_unspooled_to_spooled_cat(cm, i/nratecats);
    basecat = gpm != NULL ? gpm->gapcat_to_cat[gapped_cat] : gapped_cat;

    if (!show_cat[basecat]) continue;

    t = vec_get(hmm->begin_transitions, i);
    sink = cm_get_feature(cm, gapped_cat);
    if (t != 0) {
      if (nratecats > 1)
        printf("        begin -> \"%s-%d(%d)\" [ label = \"%.6f\" ];\n", 
               sink->chars, (i % nratecats) + 1, i, t);
      else
        printf("        begin -> \"%s(%d)\" [ label = \"%.6f\" ];\n", 
               sink->chars, i, t);
    }
  }

  for (i = 0; i < hmm->transition_matrix->size; i++) {
    gapped_cat = cm_unspooled_to_spooled_cat(cm, i/nratecats);
    basecat = gpm != NULL ? gpm->gapcat_to_cat[gapped_cat] : gapped_cat;

    if (!show_cat[basecat]) continue;
    source = cm_get_feature(cm, gapped_cat);  

    for (j = 0; j < hmm->transition_matrix->size; j++) {
      gapped_cat = cm_unspooled_to_spooled_cat(cm, j/nratecats);
      basecat = gpm != NULL ? gpm->gapcat_to_cat[gapped_cat] : gapped_cat;

      if (!show_cat[basecat]) continue;
      t = mm_get(hmm->transition_matrix, i, j);

      if (suppress_unconnected && i == j && t == 1) continue;

      sink = cm_get_feature(cm, gapped_cat);
      if (t != 0) {
        if (nratecats > 1)
          printf("        \"%s-%d(%d)\" -> \"%s-%d(%d)\" [ label = \"%.6f\" ];\n", 
                 source->chars, (i % nratecats) + 1, i, sink->chars, 
                 (j % nratecats) + 1, j, t);      
        else
          printf("        \"%s(%d)\" -> \"%s(%d)\" [ label = \"%.6f\" ];\n", 
                 source->chars, i, sink->chars, j, t);      
      }
    }
  }

  printf ("}\n");

  return 0;
}

