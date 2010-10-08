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
#include "treeGen.help"

int num_rooted_topologies(int n);
TreeNode *tr_new_trivial(char *name1, char *name2);
void tr_add_leaf_internal(TreeNode *t, int branch, char *lname, int lgroup);
void tr_add_leaf_at_root(TreeNode *t, char *lname, int lgroup);


int main(int argc, char *argv[]) {
  char c;
  int opt_idx, i, j, k, N, nleaves;
  List *names, *treelist, *newlist, *tmpl, *groups = NULL;
  TreeNode *t, *tnew;
  int *used=NULL;

  struct option long_opts[] = {
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "h", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'h':
      printf("%s", HELP);
      exit(0);
    case '?':
      die("Bad argument.  Try 'treeGen -h'.\n");
    }
  }

  if (optind < argc - 2 || optind > argc - 1)
    die("ERROR: Wrong number of arguments.  Try 'treeGen -h'.\n");

  set_seed(-1);

  names = get_arg_list(argv[optind]);

  if (lst_size(names) <= 1)
    die("ERROR: must specify at least two species names.\n");

  if (optind == argc - 2) {
    groups = get_arg_list_int(argv[optind+1]);
    if (lst_size(names) != lst_size(groups))
      die("ERROR: name list and group list must be equal in length.\n");
  }

  nleaves = lst_size(names) - 1; /* excluding outgroup */

  N = num_rooted_topologies(nleaves); 

  if (groups != NULL) {
    int maxgroup = 0;
    for (i = 0; i < lst_size(groups); i++)
      if (lst_get_int(groups, i) > maxgroup)
          maxgroup = lst_get_int(groups, i);
    used = smalloc((maxgroup+1) * sizeof(int));
    for (i = 0; i <= maxgroup; i++)
      used[i] = FALSE;
  }

  /* FIXME: eventually need to consider constraints here */

  if (N > 1e9)
    fprintf(stderr, "WARNING: very large number of topologies expected (%d).  Program may not finish.\n", N);

  /* start with tree consisting of first two names */
  t = tr_new_trivial(((String*)lst_get_ptr(names, 0))->chars, 
                     ((String*)lst_get_ptr(names, 1))->chars);

  treelist = lst_new_ptr(1000);
  newlist = lst_new_ptr(1000);
  lst_push_ptr(treelist, t);

  if (groups != NULL) {         /* use branch lengths to encode group
                                   membership -- sort of an ugly hack
                                   but should be okay here */
    t->lchild->dparent = lst_get_int(groups, 0);
    t->rchild->dparent = lst_get_int(groups, 1);
    if (t->lchild->dparent == t->rchild->dparent)
      t->dparent = t->lchild->dparent;    
    used[lst_get_int(groups, 0)] = TRUE;
    used[lst_get_int(groups, 1)] = TRUE;
  }

  for (i = 2; i < nleaves; i++) {
    char *nextname = ((String*)lst_get_ptr(names, i))->chars;
    int nextgroup = groups != NULL ? lst_get_int(groups, i) : -1;
    lst_clear(newlist);

    for (j = 0; j < lst_size(treelist); j++) {
      t = lst_get_ptr(treelist, j);

      /* create copies and add leaf to each internal branch */
      for (k = 1; k < t->nnodes; k++) {
        TreeNode *n = lst_get_ptr(t->nodes, k);

        /* decide whether adding leaf to this branch is consistent
           with monophyletic groups */
        if (groups != NULL) {
          int branchgroup = n->dparent;
          int ancgroup = n->parent->dparent;
          if (nextgroup > 0 && used[nextgroup]) { 
                                /* group is represented in the tree */
            if (nextgroup != branchgroup) {
              continue;   /* can only add to the designated subtree */
            }
          }

          else {                  /* group is zero (background) or not
                                     yet represented in the tree */ 
            if (branchgroup != 0 && nextgroup != branchgroup && 
                branchgroup == ancgroup) {
              continue;             /* only prohibit adding inside
                                       another designated subtree
                                       (adding to leading branch is
                                       okay) */
            }
          }
        }

        tnew = tr_create_copy(t);
        tr_add_leaf_internal(tnew, k, nextname, nextgroup);
        lst_push_ptr(newlist, tnew);
    }

      /* now add leaf at root; this time reuse the original copy to
         avoid unnecessary memory reallocation */
      if (nextgroup <= 0 || !used[nextgroup] || t->dparent == nextgroup) {
        tr_add_leaf_at_root(t, nextname, nextgroup);
        lst_push_ptr(newlist, t);
      }
      else
        tr_free(t);
    }

    /* swap treelist and newlist */
    tmpl = treelist;
    treelist = newlist;
    newlist = tmpl;

    if (groups != NULL)
      used[nextgroup] = TRUE;
  }

  /* traverse list and add outgroup at root of each tree */
  if (nleaves > 1) {
    for (j = 0; j < lst_size(treelist); j++) {
      t = lst_get_ptr(treelist, j);
      tr_add_leaf_at_root(t, ((String*)lst_get_ptr(names, nleaves))->chars, 0);
    }
  }

  /* print trees */
  for (j = 0; j < lst_size(treelist); j++) {
    t = lst_get_ptr(treelist, j);
    tr_print(stdout, t, FALSE);
  }

  return 0;
}

/* compute number of rooted topologies for n species: (2n-3)!! */
int num_rooted_topologies(int n) {
  int multiplier = 2*n - 3, retval = 1;
  while (multiplier > 1) {
    retval *= multiplier;
    multiplier -= 2;
  }
  return retval;
}

/* create a trivial, two-leaf tree */
TreeNode *tr_new_trivial(char *name1, char *name2) {
  TreeNode *root;
  root = tr_new_node();
  root->lchild = tr_new_node();
  strcpy(root->lchild->name, name1);
  root->lchild->parent = root;
  root->rchild = tr_new_node();
  strcpy(root->rchild->name, name2);
  root->rchild->parent = root;

  /* bypass default handling of ids and nodes list */
  root->nnodes = 3;
  root->id = 0;
  root->lchild->id = 1;
  root->rchild->id = 2;  
  root->nodes = lst_new_ptr(root->nnodes);
  lst_push_ptr(root->nodes, root);
  lst_push_ptr(root->nodes, root->lchild);
  lst_push_ptr(root->nodes, root->rchild);

  return root;
}

/* add leaf with specified name to specified internal branch */
void tr_add_leaf_internal(TreeNode *t, int branch, char *lname, int lgroup) {
  TreeNode *oldnode, *newanc, *newleaf;

  oldnode = lst_get_ptr(t->nodes, branch); /* node beneath branch in question */
  if (oldnode == t)
    die("ERROR tr_add_leaf_internal: oldnode == t\n");

  newanc = tr_new_node();
  newleaf = tr_new_node();
  strcpy(newleaf->name, lname);
  newleaf->dparent = lgroup;

  newanc->rchild = newleaf;
  newleaf->parent = newanc;
  newanc->lchild = oldnode;
  newanc->parent = oldnode->parent; 

  if (oldnode->parent->lchild == oldnode)
    oldnode->parent->lchild = newanc;
  else 
    oldnode->parent->rchild = newanc;

  oldnode->parent = newanc;

  if (lgroup > 0 && lgroup == oldnode->dparent)
    newanc->dparent = lgroup;

  /* fix up ids and nodes list */
  lst_push_ptr(t->nodes, newanc);
  newanc->id = lst_size(t->nodes) - 1; /* circumvent normal id assignment */
  lst_push_ptr(t->nodes, newleaf);
  newleaf->id = lst_size(t->nodes) - 1;
  t->nnodes += 2;
}

/* add a leaf with specified name to root branch */
void tr_add_leaf_at_root(TreeNode *t, char *lname, int lgroup) {
  TreeNode *newanc, *newleaf;

  newanc = tr_new_node();
  newleaf = tr_new_node();
  strcpy(newleaf->name, lname);
  newleaf->dparent = lgroup;

  /* we don't want to change the identity of the root node, so will
     add the new node below it and rewire as necessary */
  newanc->lchild = t->lchild;
  newanc->rchild = t->rchild;
  t->lchild->parent = newanc;
  t->rchild->parent = newanc;
  t->lchild = newanc;
  t->rchild = newleaf;
  newanc->parent = t;
  newleaf->parent = t;

  newanc->dparent = t->dparent;

  if (lgroup == newanc->dparent) 
    t->dparent = lgroup;    
  else
    t->dparent = 0; 

  /* fix up ids and nodes list */
  lst_push_ptr(t->nodes, newanc);
  newanc->id = lst_size(t->nodes) - 1; /* circumvent normal id assignment */
  lst_push_ptr(t->nodes, newleaf);
  newleaf->id = lst_size(t->nodes) - 1;
  t->nnodes += 2;
}
