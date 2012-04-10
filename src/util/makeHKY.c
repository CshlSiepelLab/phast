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
#include <tree_model.h>
#include <prob_vector.h>
#include <subst_mods.h>
#include "makeHKY.help"

#define ALPHABET "ACGT"

int main(int argc, char *argv[]) {
  char c;
  int opt_idx;

  double gc = 0.4, t = 1, kappa;
  TreeNode *tree = NULL;
  Vector *pi = NULL;
  List *l;
  int i;
  TreeModel *mod;
  char tmpstr[50];

  struct option long_opts[] = {
    {"gc", 1, 0, 'g'},
    {"pi", 1, 0, 'p'},
    {"branch-length", 1, 0, 't'},
    {"tree", 1, 0, 'T'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "g:p:t:T:h", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'g':
      gc = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'p':
      l = get_arg_list_dbl(optarg);
      if (lst_size(l) != 4)
        die("ERROR: argument of --pi must be list of size 4.\n");
      pi = vec_new(4);
      for (i = 0; i < lst_size(l); i++)
        pi->data[i] = lst_get_dbl(l, i);
      pv_normalize(pi);
      break;
    case 't':
      t = get_arg_dbl_bounds(optarg, 0, INFTY);
      break;
    case 'T':
      tree = tr_new_from_file(phast_fopen(optarg, "r"));
      break;
    case 'h':
      printf("%s", HELP);
      exit(0);
    case '?':
      die("Bad argument.  Try 'program -h'.\n");
    }
  }

  if (optind != argc - 1) 
    die("Bad arguments.  Try 'makeHKY -h'.\n");

  set_seed(-1);
    
  kappa = get_arg_dbl_bounds(argv[optind], 0, INFTY);

  if (pi == NULL) {
    pi = vec_new(4);
    pi->data[0] = (1-gc)/2;
    pi->data[1] = gc/2;
    pi->data[2] = gc/2;
    pi->data[3] = (1-gc)/2;
  }

  if (tree == NULL) {
    sprintf(tmpstr, "(s1:%f,s2:%f)", t/2, t/2);
    tree = tr_new_from_string(tmpstr);
  }

  mod = tm_new(tree, NULL, pi, HKY85, ALPHABET, 1, 1, NULL, -1);
  tm_set_HKY_matrix(mod, kappa, -1);
  tm_scale_rate_matrix(mod);

  tm_print(stdout, mod);

  return 0;
}
