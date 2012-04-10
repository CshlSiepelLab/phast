/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: draw_tree.c,v 1.5 2008-11-12 02:07:58 acs Exp $ */

#include <stdlib.h>
#include <stdio.h>
#include <trees.h>
#include <tree_model.h>
#include <getopt.h>

void print_usage() {
  fprintf(stderr, "\n\
PROGRAM:        draw_tree\n\
DESCRIPTION:    Produces a simple postscript rendering of a tree.\n\
USAGE:          draw_tree [-dbvs] <tree.nh>|<tree.mod>\n\
\n\
OPTIONS:\n\
    <tree_fname>    (Required) File name of tree (expected to be in \n\
                    Newick format, unless filename ends with '.mod', in\n\
                    which case expected to be a tree model file).\n\
    -d              Print \"diagonal\" branches, instead of \"right-angle\" or \n\
                    \"square\" ones (produces a \"cladogram\", as opposed to a \n\
                    \"phenogram\").  This option implies -s.\n\
    -b              Suppress branch lengths.\n\
    -v              Vertical layout.\n\
    -s              Don't draw branches to scale.\n\n");
}

int main(int argc, char *argv[]) {
  int horizontal = 1;
  int square_branches = 1;
  int show_branch_lens = 1;
  int draw_to_scale = 1;

  TreeNode *tree;
  char c;
  String *suffix;

  while ((c = getopt(argc, argv, "dbvsh")) != -1) {
    switch(c) {
    case 'd':
      square_branches = 0;
      draw_to_scale = 0;
      break;
    case 'b':
      show_branch_lens = 0;
      break;
    case 'v':
      horizontal = 0;
      break;
    case 's':
      draw_to_scale = 0;
      break;
    case 'h':
      print_usage();
      exit(1);
    case '?':
      print_usage();
      exit(1);
    }
  }

  if (optind >= argc) {
    print_usage();
    exit(1);
  }
  
  set_seed(-1);

  suffix = str_new_charstr(argv[optind]);
  str_suffix(suffix, '.');
  if (str_equals_charstr(suffix, "mod")) {
    TreeModel *tmp = tm_new_from_file(phast_fopen(argv[optind], "r"), 1);
    tree = tmp->tree;
  }
  else 
    tree = tr_new_from_file(phast_fopen(argv[optind], "r"));

  tr_print_ps(stdout, tree, show_branch_lens, square_branches, draw_to_scale, 
              horizontal);

  return 0;
}
