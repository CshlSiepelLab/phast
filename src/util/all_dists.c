#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <misc.h>
#include <trees.h>

void usage(char *prog) {
  printf("\n\
PROGRAM:      %s\n\
DESCRIPTION:  Given a tree in Newick (*.nh) format, reports distances\n\
              between all pairs of leaves.\n\
\n\
USAGE:        %s <tree.nh>\n\
OPTIONS:\n\
    --help, -h\n\
        Print this help message.\n\n", prog, prog);
  exit(0);
}

int main(int argc, char *argv[]) {
  FILE *INF;
  char c;
  int i, j, opt_idx;
  TreeNode *tree, *n, *node_i, *node_j, *lca;
  List *leaves;


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
    
  tree = parse_nh_from_file(fopen_fname(argv[optind], "r"));

  /* obtain list of leaves */
  leaves = lst_new_ptr((tree->nnodes+1)/2);
  for (i = 0; i < lst_size(tree->nodes); i++) {
    n = lst_get_ptr(tree->nodes, i);
    if (n->lchild == NULL && n->rchild == NULL)
      lst_push_ptr(leaves, n);
  }

  /* look at all pairs */
  for (i = 0; i < lst_size(leaves); i++) {
    node_i = lst_get_ptr(leaves, i);
    for (j = i+1; j < lst_size(leaves); j++) {
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
      printf ("%-10s %-10s %f\n", node_i->name, node_j->name, dist);
    }
  }

  return 0;
}
