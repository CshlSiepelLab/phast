/* prior distributions over trees and branch lengths, for use in variational inference */

#ifndef TREEPRIOR_H
#define TREEPRIOR_H

#include <stdio.h>
#include <phast/trees.h>
#include <phast/tree_model.h>
#include <phast/matrix.h>

struct cvdat;

/* prior under a simple Yule model with a relaxed local clock, defined
   by branchwise draws from a lognormal distribution with mean one */
typedef struct {
  double relclock_sig; /* stdev of log normal prior for relaxed clock */
  double relclock_sig_grad; /* gradient of same, for use in gradient ascent */
  Vector *nodetimes; /* absolute time of each node */
  Vector *nodetimes_grad; /* gradient of same */
} TreePrior;

/* package up a tree node and its distance from the root for sorting */
typedef struct {
  TreeNode *node;
  double dist_from_root;
} NodeDist;

TreePrior *tp_new();

double tp_compute_log_prior(TreeModel *mod, struct cvdat *data, Vector *branchgrad);

List *tp_sort_nodes(TreeNode *tree);

void tp_init_nodetimes(TreePrior *tp, TreeModel *mod);

#endif 
