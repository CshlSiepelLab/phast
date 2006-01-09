/* $Id: trees.h,v 1.15 2006-01-09 21:53:58 acs Exp $
   Written by Adam Siepel, 2002
   Copyright 2002, Adam Siepel, University of California */

#ifndef TREES_H
#define TREES_H

#define MAX_TREESTR_LEN 10000
#define MAX_LINE_LEN 10000

#include <stdio.h>
#include <lists.h>
#include <stringsplus.h>

typedef struct tree_node TreeNode;

/* TODO: starting to accumulate a lot of data about the entire tree at
   each node; perhaps time to adopt a less purely recursive structure,
   with a higher level struct describing the tree as a whole. */
struct tree_node {
  TreeNode *parent;
  TreeNode *lchild, *rchild;
  double dparent;
  char name [STR_MED_LEN];
  void *data;                   /* allows generic data to be 
                                   associated with tree node */ 
  int id;
  int nnodes;                   /* number of nodes in subtree defined
                                   by this node */
  int height;                   /* height of this node in the tree
                                   (maximum distance to a leaf, in
                                   terms of number of edges) */

  /* the following attributes are only guaranteed to be defined for
     the TreeNode at the root of a tree; 'nodes' is defined upon
     initialization, but the others are constructed on demand; they
     should be accessed only via the functions tr_preorder,
     tr_inorder, tr_postorder.  Note that most of these auxiliary
     attributes (nnodes, height, nodes, preorder, inorder) are not
     guaranteed to remain correct if the structure of a tree is
     altered after initialization */

  List *nodes;                  /* List of nodes: ith element is a
                                   pointer to the node with id = i */
  List *preorder;               /* List of nodes in the order of a
                                   preorder traversal from the root */
  List *inorder;                /* List of nodes in the order of an
                                   inorder traversal from the root */
  List *postorder;              /* List of nodes in the order of a
                                   postorder traversal from the
                                   root */
};

TreeNode *tr_new_from_file(FILE *f);
TreeNode *tr_new_from_string(char *s);
TreeNode *tr_new_node();
void tr_add_child(TreeNode *parent, TreeNode *child);
void tr_print(FILE* f, TreeNode *root, int show_branch_lengths);
void tr_print_recur(FILE* f, TreeNode *n, int show_branch_lengths);
void tr_free(TreeNode *n);
void tr_set_nnodes(TreeNode *tree);
void tr_reset_id();
void tr_cpy(TreeNode *dest, TreeNode *src);
TreeNode *tr_create_copy(TreeNode *src);
void tr_node_cpy(TreeNode *dest, TreeNode *src);
void tr_print_ordered(FILE* f, TreeNode *root, int show_branch_lengths);
void tr_print_ordered_recur(FILE* f, TreeNode *n, int *left_right,
                            int show_branch_lengths);
List *tr_preorder(TreeNode *tr);
List *tr_inorder(TreeNode *tr);
List *tr_postorder(TreeNode *tr);
void tr_layout_xy(TreeNode *tree, int x0, int y0, int x1, int y1, 
                  int *x, int *y, int use_branch_lens, int horizontal);
void tr_print_ps(FILE *F, TreeNode *tree, int show_branch_lens,
                 int square_branches, int use_branch_lens, 
                 int horizontal_layout);
double tr_total_len(TreeNode *t);
double tr_total_len_subtree(TreeNode *sub_root);
TreeNode *tr_get_node(TreeNode *t, char *name);
void tr_scale(TreeNode *t, double scale_const);
void tr_scale_subtree(TreeNode *t, TreeNode *sub, double scale_const);
void tr_prune(TreeNode **t, List *names, int all_but);
TreeNode *tr_lca(TreeNode *tree, List *names);
TreeNode *tr_hybrid(TreeNode *sub, TreeNode *super);
double tr_scale_by_subtree(TreeNode *tree, TreeNode *sub);
void tr_partition_leaves(TreeNode *tree, TreeNode *sub, List *inside, 
                         List *outside);
void tr_partition_nodes(TreeNode *tree, TreeNode *sub, List *inside, 
			List *outside);
List *tr_leaf_names(TreeNode *tree);
void tr_name_ancestors(TreeNode *tree);
void tr_print_nodes(FILE *F, TreeNode *tree);
void tr_reroot(TreeNode *tree, TreeNode *newroot, int include_branch);
int* tr_in_subtree(TreeNode *t, TreeNode *sub);

#endif 
