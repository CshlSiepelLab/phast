/* $Id: draw_tree.c,v 1.3 2004-07-29 23:31:38 acs Exp $
   Written by Adam Siepel, 2002
   Copyright 2002, Adam Siepel, University of California */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <trees.h>
#include <getopt.h>

void print_usage() {
  fprintf(stderr, "\n\
PROGRAM:        draw_tree\n\
DESCRIPTION:    Produces a simple postscript rendering of a tree.\n\
USAGE:          draw_tree [-dbvs] <tree_fname>\n\
\n\
OPTIONS:\n\
    <tree_fname>    (Required) File name of tree (file must be in \n\
                    Newick format).\n\
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

  FILE *F;
  TreeNode *tree;
  char c;

  while ((c = getopt(argc, argv, "dbvs")) != -1) {
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
    case '?':
      print_usage();
      exit(1);
    }
  }

  if (optind >= argc) {
    print_usage();
    exit(1);
  }

  if ((F = fopen(argv[optind], "r")) == NULL) {
    fprintf(stderr, "ERROR: Cannot open %s.\n", argv[optind]);
    exit(1);
  }
  tree = tr_new_from_file(F);
  tr_print_ps(stdout, tree, show_branch_lens, square_branches, draw_to_scale, 
              horizontal);
  fclose(F);

  return 0;
}
