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
#include <trees.h>
#include <tree_model.h>

void usage(char *prog) {
  printf("\n\
PROGRAM:      %s\n\
\n\
DESCRIPTION:  Given a tree in Newick (*.nh) format, report distances\n\
              between all pairs of leaves.  If multiple files are given,\n\
              then distances are computed by averaging across models,\n\
              and statistics describing the errors in the estimates\n\
              are reported (can be useful for bootstrapping; see\n\
              'phyloBoot --dump-mods').\n\
\n\
USAGE:        %s <tree.nh> [<tree2.nh> <tree3.nh>...]\n\
\n\
OPTIONS:\n\
    --mod, -m\n\
        Read from tree model (*.mod) file(s) instead of Newick file.\n\
\n\
    --tree, -t <file>|<string>\n\
        Use leaf names from given tree.  Useful when primary files\n\
        use numbers rather than names.\n\
\n\
    --help, -h\n\
        Print this help message.\n\n", prog, prog);
  exit(0);
}

int main(int argc, char *argv[]) {
  char c;
  int i, j, t, opt_idx, ntrees, nleaves = -1;
  TreeNode *n, *node_i, *node_j, *lca, *nametree = NULL;
  TreeNode **tree;
  List *leaves, ***distance, *tree_fnames, *tot_dist;
  int mod = FALSE;
  char **leaf_name;
  String *trees_arg;
  FILE *F;

  struct option long_opts[] = {
    {"mod", 0, 0, 'm'},
    {"tree", 1, 0, 't'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "mt:h", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'm':
      mod = TRUE;
      break;
    case 't':
      if (optarg[0] == '(')
        nametree = tr_new_from_string(optarg);
      else 
        nametree = tr_new_from_file(phast_fopen(optarg, "r"));
      break;
    case 'h':
      usage(argv[0]);
    case '?':
      die("Bad argument.  Try '%s -h'.\n", argv[0]);
    }
  }

  if (optind > argc - 1) 
    die("Input filename required.  Try '%s -h'.\n", argv[0]);

  set_seed(-1);

  /* build a comma-delimited list and pass to get_arg_list; allows
     possibility of reading from file via '*' operator */
  trees_arg = str_new(1000);
  for (i = optind; i < argc; i++) {
    str_append_charstr(trees_arg, argv[i]);
    if (i < argc - 1) str_append_char(trees_arg, ',');
  }
  tree_fnames = get_arg_list(trees_arg->chars);

  ntrees = lst_size(tree_fnames);
  tree = smalloc(ntrees * sizeof(void*));

  /* read trees */
  for (t = 0; t < ntrees; t++) {
    String *fname = lst_get_ptr(tree_fnames, t);
    if (mod) {
      TreeModel *m = tm_new_from_file(F = phast_fopen(fname->chars, "r"), 1);
      tree[t] = tr_create_copy(m->tree);
      tm_free(m);
      phast_fclose(F);
    }
    else
      tree[t] = tr_new_from_file(phast_fopen(fname->chars, "r"));
  }

  /* initialization */
  nleaves = (tree[0]->nnodes + 1)/2;
  leaves = lst_new_ptr(nleaves);    
  distance = smalloc(nleaves * sizeof(void*));
  leaf_name = smalloc(nleaves * sizeof(void*));
  for (i = 0; i < nleaves; i++) {
    distance[i] = smalloc(nleaves * sizeof(void*));
    for (j = i+1; j < nleaves; j++) 
      distance[i][j] = lst_new_dbl(ntrees);
  }
  if (nametree == NULL) nametree = tree[0];
  for (i = 0, j = 0; i < lst_size(nametree->nodes); i++) {
    n = lst_get_ptr(nametree->nodes, i);
    if (n->lchild == NULL && n->rchild == NULL)
      leaf_name[j++] = n->name;
  }
  tot_dist = lst_new_dbl(ntrees);

  /* now compute distances */
  for (t = 0; t < ntrees; t++) {
    /* obtain list of leaves */
    lst_clear(leaves);
    for (i = 0; i < lst_size(tree[t]->nodes); i++) {
      n = lst_get_ptr(tree[t]->nodes, i);
      if (n->lchild == NULL && n->rchild == NULL)
        lst_push_ptr(leaves, n);
    }

    if (lst_size(leaves) != nleaves)
      die("ERROR: trees have different numbers of leaves.\n");

    /* look at all pairs */
    for (i = 0; i < nleaves; i++) {
      node_i = lst_get_ptr(leaves, i);
      for (j = i+1; j < nleaves; j++) {
        double dist = 0;
        node_j = lst_get_ptr(leaves, j);
        /* because ids are assigned in pre-order, the first ancestor of
           node j that has an id less than i is the LCA of i and j; we
           seek the sum of distances from both i and j to this node */
        for (n = node_j; n->id >= node_i->id; n = n->parent)
          dist += n->dparent;      
        lca = n;
        for (n = node_i; n != lca; n = n->parent)
          dist += n->dparent;            
        lst_push_dbl(distance[i][j], dist);
      }
    }
    lst_push_dbl(tot_dist, tr_total_len(tree[t]));
  }


  /* print distances and (optionally) stats */
  if (ntrees == 1) {
    for (i = 0; i < nleaves; i++) {
      for (j = i+1; j < nleaves; j++) {
        printf ("%s\t%s\t%f\n", leaf_name[i], leaf_name[j], 
                lst_get_dbl(distance[i][j], 0));
      }
    }
    printf ("%s\t%s\t%f\n", "(total)", "-", lst_get_dbl(tot_dist, 0));
  }
  else {
    double mean, stdev;
    double quantiles[] = {0, 0.025, 0.05, 0.5, 0.95, 0.975, 1};
    double quantile_vals[7]; 

    printf("%-15s %-15s %9s %9s %9s %9s %9s %9s %9s %9s %9s\n", "leaf1", 
           "leaf2", "mean", "stdev", "median", "min", "max", "95%_min", 
           "95%_max", "90%_min", "90%_max");

    for (i = 0; i < nleaves; i++) {
      for (j = i+1; j < nleaves; j++) {
        mean = lst_dbl_mean(distance[i][j]);
        stdev = lst_dbl_stdev(distance[i][j]);
        lst_qsort_dbl(distance[i][j], ASCENDING);
        lst_dbl_quantiles(distance[i][j], quantiles, 7, quantile_vals);

        printf("%-15s %-15s %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n", 
               leaf_name[i], leaf_name[j], mean, stdev, quantile_vals[3], quantile_vals[0], 
               quantile_vals[6], quantile_vals[1], quantile_vals[5], quantile_vals[2], 
               quantile_vals[4]);
      }
    }

    /* also do total branch len */
    mean = lst_dbl_mean(tot_dist);
    stdev = lst_dbl_stdev(tot_dist);
    lst_qsort_dbl(tot_dist, ASCENDING);
    lst_dbl_quantiles(tot_dist, quantiles, 7, quantile_vals);
    
    printf("%-15s %-15s %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n", 
	   "(total)", "-", mean, stdev, quantile_vals[3], quantile_vals[0], 
	   quantile_vals[6], quantile_vals[1], quantile_vals[5], quantile_vals[2], 
	   quantile_vals[4]);
  }

  return 0;
}
