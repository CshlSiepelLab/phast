/* $Id: trees.c,v 1.7 2004-06-22 07:29:47 acs Exp $ 
   Written by Adam Siepel, 2002
   Copyright 2002, Adam Siepel, University of California */

/** \file trees.c 
  Functions for manipulating phylogenetic trees.  Includes functions
  for reading, writing, traversing, and printing.  Trees are
  represented as rooted binary trees with non-negative real branch
  lengths.  
  \ingroup base
*/

/* TODO: 
   - reasonable way to handle unrooted trees (e.g., have "virtual root", 
     understood to be ternary) 
   - better error checking in parsing function
   - make printing routines nonrecursive
   - function names: consistent prefix (e.g., tr_)
   - use tr_inorder, tr_preorder, tr_postorder where possible 
   - possibly new scheme for labels at leaves; need general way to go 
     between leaf labels and leaf numbers (perhaps wrt a list of names) 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "stacks.h"
#include "trees.h"
#include "misc.h"
#include "stringsplus.h"


/* NOTE: when parsed from NH, node ids correspond to a preorder
   traversal of the tree.  A number of useful properties result.  For
   example, if two nodes u and v are labeled such that v has the larger
   id, then the first ancestor of v that has an id smaller than that
   of u is the LCA of u and v */

static int idcounter = 0;

/* coords for postscript printing */
/* top-left x */
#define TL_X 15                 

/* top-left y */
#define TL_Y 10                 
                                
/* bottom-right x */
#define BR_X 500                
                                
/* bottom-right y */
#define BR_Y 700                
                                

/** Parse a single tree from a file in "New Hampshire" format */
TreeNode *parse_nh_from_file(FILE *f) { 
  String *s = str_new(STR_VERY_LONG_LEN);
  TreeNode *retval;

  str_slurp(s, f);

  str_double_trim(s);
  if (s->chars[0] != '(')
    die("ERROR: Can't parse tree file (Newick).\n");

  if (s->chars[s->length-1] == ';') 
    s->chars[--s->length] = '\0';

  retval = parse_nh_from_string(s->chars);
  str_free(s);
  return retval;
}

/** Parse a single New Hampshire-formatted tree from a string */
TreeNode *parse_nh_from_string(char *treestr) { 
  char diststr[STR_MED_LEN];
  TreeNode *root, *node, *newnode;
  int i, in_distance = 0, len = strlen(treestr);
  char c;
  char *currentname = NULL;

  idcounter = 0;                /* start at 0 for each tree */
  root = new_tree_node(); root->nnodes = 1;
  node = root;
  for (i = 0; i < len; i++) {
    c = treestr[i];

    if (in_distance) {
      if (isdigit(c) || c == '.' || c == '-') {
        strncat(diststr, &c, 1);
        continue;
      }
      else 
        node->dparent = atof(diststr);
      in_distance = 0;
    }

    if (c == '(') {
      addchild(node, newnode = new_tree_node());
      node = newnode;
      currentname = newnode->name;
      root->nnodes++;
    }
    else if (c == ',') {
      addchild(node->parent, newnode = new_tree_node());
      node = newnode;
      currentname = node->name;
      root->nnodes++;
    }
    else if (c == ')') {
      node = node->parent;
      currentname = NULL;
    }
    else if (c == ':') {
      diststr[0] = '\0';
      in_distance = 1;
    }
    else if (currentname != NULL) {
      if (!isspace(c) || strlen(currentname) > 0) /* avoid leading spaces */
        strncat(currentname, &c, 1);
    }
  }

  tr_set_nnodes(root);
  return root;
}

/* traverse tree to set nnodes at each node (postorder); also set
   height at each node, and create "nodes" list */
void tr_set_nnodes(TreeNode *tree) {
  Stack *stack;
  TreeNode *node;
  int i, j;

  tree->nodes = lst_new_ptr(tree->nnodes);
  for (i = 0; i < tree->nnodes; i++) lst_push_ptr(tree->nodes, NULL);
  stack = stk_new_ptr(tree->nnodes);
  stk_push_ptr(stack, tree);
  while ((node = stk_pop_ptr(stack)) != NULL) {
    assert((node->lchild == NULL && node->rchild == NULL) ||
           (node->lchild != NULL && node->rchild != NULL));

    if (node->lchild == NULL) {
      node->nnodes = 1;
      node->height = 0;
      if (node->id >= tree->nnodes) 
        for (j = tree->nnodes; j <= node->id; j++) 
          lst_push_ptr(tree->nodes, NULL); /* this hack necessary
                                              because original estimate
                                              of size of list may be
                                              wrong */
      lst_set_ptr(tree->nodes, node->id, node);
    }
    else if (node->lchild->nnodes != -1 && node->rchild->nnodes != -1) {
      node->nnodes = node->lchild->nnodes + node->rchild->nnodes + 1;
      node->height = max(node->lchild->height, node->rchild->height) + 1;

      for (j = tree->nnodes; j <= node->id; j++) 
        lst_push_ptr(tree->nodes, NULL); /* this hack necessary
                                            because original estimate
                                            of size of list may be
                                            wrong */
      lst_set_ptr(tree->nodes, node->id, node);
    }
    else {			/* internal node whose children have
                       not yet been visited */
      stk_push_ptr(stack, node);
      stk_push_ptr(stack, node->lchild);
      stk_push_ptr(stack, node->rchild);
    }
  }
  stk_free(stack);
}

/* Create and initialize a new tree node */
TreeNode *new_tree_node() {
  TreeNode *n = (TreeNode*)smalloc(sizeof(TreeNode));
  n->parent = n->lchild = n->rchild = NULL;
  n->data = NULL;
  n->id = idcounter++;
  n->name[0] = '\0';
  n->dparent = 0;
  n->nnodes = -1;
  n->height = 0;
  n->nodes = n->preorder = n->inorder = n->postorder = NULL;
  return(n);
}

/* Add specified child to specified parent, creating all requisite
   links.  If the parent already has two children, add a new node
   (to simulate an nary tree with a binary tree) */
void addchild(TreeNode *parent, TreeNode *child) {
  if (parent->lchild == NULL) {
    parent->lchild = child;
  }
  else if (parent->rchild == NULL) {
    parent->rchild = child;
  }
  else {
    /* add intermediate node to accommodate extra child */
    TreeNode *tmp = new_tree_node();
    tmp->lchild = parent->lchild;
    tmp->rchild = parent->rchild;
    tmp->parent = parent;
    tmp->lchild->parent = tmp;
    tmp->rchild->parent = tmp;
    parent->lchild = tmp;
    parent->rchild = child;
  }
  child->parent = parent;
}

/* Print tree in New Hampshire format */
void print_tree(FILE* f, TreeNode *root, int show_branch_lengths) {
  print_tree_recur(f, root, show_branch_lengths);
  fprintf(f, ";\n");
}

/* Recursive subroutine used by print_tree */
void print_tree_recur(FILE* f, TreeNode *n, int show_branch_lengths) {

  assert((n->lchild == NULL && n->rchild == NULL) || 
	 (n->lchild != NULL && n->rchild != NULL));

  if (n->lchild != NULL) {
    fprintf(f, "(");
    print_tree_recur(f, n->lchild, show_branch_lengths);
    fprintf(f, ",");
    print_tree_recur(f, n->rchild, show_branch_lengths);
    fprintf(f, ")");
  }
  else {
    fprintf(f, "%s", n->name);
  }

  if (show_branch_lengths && n->parent != NULL)
    fprintf(f, ":%f", n->dparent);
}

/* For debugging purposes */
void print_tree_debug(FILE *f, TreeNode *n) {
  assert((n->lchild == NULL && n->rchild == NULL) || 
	 (n->lchild != NULL && n->rchild != NULL));
  print_node_debug(f, n);
  if (n->lchild != NULL) {
    print_tree_debug(f, n->lchild);
    print_tree_debug(f, n->rchild);
  }  
}

void print_node_debug(FILE *f, TreeNode *n) {
  fprintf(f, "\nNODE %d\n", n->id);
  if (strlen(n->name) > 0)
    fprintf(f, "\tNAME = %s\n", n->name);
  fprintf(f, "\tNNODES = %d\n", n->nnodes);
  fprintf(f, "\tHEIGHT = %d\n", n->height);
  fprintf(f, "\tParent = %d\n", n->parent == NULL ? -1 : n->parent->id);
  if (n->lchild != NULL) 
    fprintf(f, "\tChildren = (%d, %d)\n", n->lchild->id, n->rchild->id);
}

void free_tree(TreeNode *tree) {
  Stack *stack;
  TreeNode *n;
  stack = stk_new_ptr(tree->nnodes);
  stk_push_ptr(stack, tree);
  while ((n = stk_pop_ptr(stack)) != NULL) {
    if (n->lchild != NULL) stk_push_ptr(stack, n->lchild);
    if (n->rchild != NULL) stk_push_ptr(stack, n->rchild);
    if (n->nodes != NULL) lst_free(n->nodes);
    if (n->preorder != NULL) lst_free(n->preorder);
    if (n->inorder != NULL) lst_free(n->inorder);
    if (n->postorder != NULL) lst_free(n->postorder);
    free(n);
  }
  stk_free(stack);
}

void tr_reset_id() {
  idcounter = 0;
}

void tr_cpy(TreeNode *dest, TreeNode *src) {
  Stack *stack, *nodes, *cpystack;
  TreeNode *n, *ncpy, *lcpy, *rcpy;

  /* first just flatten dest into a list; we won't try to match
   * corresponding nodes, point is just to reuse memory */
  stack = stk_new_ptr(src->nnodes);
  nodes = stk_new_ptr(src->nnodes);
  cpystack = stk_new_ptr(src->nnodes);
  stk_push_ptr(stack, dest); 
  while ((n = stk_pop_ptr(stack)) != NULL) {
    if (n != dest)              /* leave the root aside */
      stk_push_ptr(nodes, n);

    if (n->lchild != NULL)
      stk_push_ptr(stack, n->lchild);
    if (n->rchild != NULL)
      stk_push_ptr(stack, n->rchild);
  }

  /* now copy node by node */
  tr_node_cpy(dest, src);       /* copy root */
  stk_push_ptr(stack, src); 
  stk_push_ptr(cpystack, dest);
  while ((n = stk_pop_ptr(stack)) != NULL) {
    ncpy = stk_pop_ptr(cpystack); 

    if (n->lchild != NULL) {
      lcpy = stk_pop_ptr(nodes);
      tr_node_cpy(lcpy, n->lchild);
      addchild(ncpy, lcpy);
      stk_push_ptr(stack, n->lchild);
      stk_push_ptr(cpystack, lcpy);
    }
    if (n->rchild != NULL) {
      rcpy = stk_pop_ptr(nodes);
      tr_node_cpy(rcpy, n->rchild);
      addchild(ncpy, rcpy);
      stk_push_ptr(stack, n->rchild);
      stk_push_ptr(cpystack, rcpy);
    }
  }

  stk_free(stack);
  stk_free(cpystack);
  stk_free(nodes);
}

TreeNode *tr_create_copy(TreeNode *src) {
  Stack *stack, *cpystack;
  TreeNode *n, *ncpy, *lcpy, *rcpy, *dest;
  
  tr_reset_id();

  stack = stk_new_ptr(src->nnodes);
  cpystack = stk_new_ptr(src->nnodes);
  stk_push_ptr(stack, src); 
  dest = new_tree_node();
  tr_node_cpy(dest, src);
  stk_push_ptr(cpystack, dest);
  while ((n = stk_pop_ptr(stack)) != NULL) {
    ncpy = stk_pop_ptr(cpystack);
    if (n->lchild != NULL) {
      lcpy = new_tree_node();
      tr_node_cpy(lcpy, n->lchild);
      addchild(ncpy, lcpy);
      stk_push_ptr(stack, n->lchild);
      stk_push_ptr(cpystack, lcpy);
    }
    if (n->rchild != NULL) {
      rcpy = new_tree_node();
      tr_node_cpy(rcpy, n->rchild);
      addchild(ncpy, rcpy);
      stk_push_ptr(stack, n->rchild);
      stk_push_ptr(cpystack, rcpy);
    }
  }  
  stk_free(stack);
  stk_free(cpystack);

  /* leave inorder, preorder, postorder NULL at dest; will be created
     on demand (can't copy from src because pointers will be wrong) */

  dest->nnodes = src->nnodes;
  tr_set_nnodes(dest);
  return dest;
}


/* Copy contents of tree node (ignore pointers) */
void tr_node_cpy(TreeNode *dest, TreeNode *src) {
  dest->parent = dest->lchild = dest->rchild = NULL;
/*   dest->data = NULL; */ /* src->data; */       /* WARNING: only copying pointer! */
  dest->id = src->id;
  strcpy(dest->name, src->name); 
  dest->dparent = src->dparent;
/*   dest->nnodes = -1; */
/*   dest->height = src->height; */

  /* don't copy data, nnodes, height, preorder, inorder, postorder */
}


/* Print tree in New Hampshire format.  This version imposes an
 * ordering on the leaves (useful when comparing trees that have been
 * rearranged).  At every internal node, we store the name of the leaf
 * beneath it that comes first alphanumerically.  When recursively
 * printing the tree, at each internal node, we call its children in
 * the order of these names.  */
void print_tree_alph(FILE* f, TreeNode *root, int show_branch_lengths) {
  int *left_right, *mark;
  char **names;
  TreeNode *n;
  Stack *stack;
  int i;

  left_right = (int*)smalloc(root->nnodes * sizeof(int));
  mark = (int*)smalloc(root->nnodes * sizeof(int));
  names = (char**)smalloc(root->nnodes * sizeof(char*));
  for (i = 0; i < root->nnodes; i++) { 
    left_right[i] = 0; 
    mark[i] = 0; 
    names[i] = NULL;
  }

  stack = stk_new_ptr(root->nnodes);
  stk_push_ptr(stack, root);
  while ((n = stk_pop_ptr(stack)) != NULL) {
    assert((n->lchild == NULL && n->rchild == NULL) ||
           (n->lchild != NULL && n->rchild != NULL));
    if (n->lchild == NULL) {
      names[n->id] = n->name;
      mark[n->id] = 1;
    }
    else if (mark[n->lchild->id] == 1 && mark[n->rchild->id] == 1) {
      if (names[n->rchild->id] == NULL)
        left_right[n->id] = 1;
      else if (names[n->lchild->id] == NULL) 
        left_right[n->id] = 0;
      /* now we know neither is NULL */
      else if (strcmp(names[n->lchild->id], names[n->rchild->id]) <= 0)
        left_right[n->id] = 1;
      else
        left_right[n->id] = 0;

      names[n->id] = left_right[n->id] == 1 ? names[n->lchild->id] :
        names[n->rchild->id];
      
      mark[n->id] = 1;
    }
    else {                      /* push back on stack for post-order
                                 * traversal */
      stk_push_ptr(stack, n);
      stk_push_ptr(stack, n->rchild);
      stk_push_ptr(stack, n->lchild);
    }
  }

  print_tree_ordered_recur(f, root, left_right, show_branch_lengths);
  fprintf(f, ";\n");
  
  stk_free(stack);
  free(left_right);
  free(mark);
  free(names);
}

/* Recursive subroutine for printing trees that allows an arbitrary
 * left/right ordering of the children at every internal node. */
void print_tree_ordered_recur(FILE* f, TreeNode *n, int *left_right,
                              int show_branch_lengths) {

  assert((n->lchild == NULL && n->rchild == NULL) || 
	 (n->lchild != NULL && n->rchild != NULL));

  if (n->lchild != NULL) {
    fprintf(f, "(");
    if (left_right[n->id]) {
      print_tree_ordered_recur(f, n->lchild, left_right, show_branch_lengths);
      fprintf(f, ",");
      print_tree_ordered_recur(f, n->rchild, left_right, show_branch_lengths);
    }
    else {
      print_tree_ordered_recur(f, n->rchild, left_right, show_branch_lengths);
      fprintf(f, ",");
      print_tree_ordered_recur(f, n->lchild, left_right, show_branch_lengths);
    }
    fprintf(f, ")");
  }
  else {
    fprintf(f, "%s", n->name);
  }

  if (show_branch_lengths)
    fprintf(f, ":%f", n->dparent);
}

/* below routines provide lists that represent preorder, inorder, and
   postorder traversals of the specified tree.  Doing things this way
   simplifies calling code considerably, and also can substantially
   increase efficiency (avoid lots of alloc and srealloc of stacks;
   avoid two-pass problem with inorder and postorder) */
List *tr_preorder(TreeNode *tr) {
  if (tr->preorder == NULL) {   /* produce on demand */
    Stack *stack;
    TreeNode *n;

    tr->preorder = lst_new_ptr(tr->nnodes);
    stack = stk_new_ptr(tr->nnodes);
    stk_push_ptr(stack, tr);
    while ((n = stk_pop_ptr(stack)) != NULL) {
      assert((n->lchild == NULL && n->rchild == NULL) ||
             (n->lchild != NULL && n->rchild != NULL));
      if (n->lchild != NULL) {
        stk_push_ptr(stack, n->rchild);
        stk_push_ptr(stack, n->lchild);
      }
      lst_push_ptr(tr->preorder, n);
    }
    stk_free(stack);
  }

  return tr->preorder;
}

List *tr_inorder(TreeNode *tr) {
  if (tr->inorder == NULL) {    /* produce on demand */
    int i;
    int *mark;
    Stack *stack;
    TreeNode *n;

    tr->inorder = lst_new_ptr(tr->nnodes);
    mark = (int*)smalloc(tr->nnodes * sizeof(int));
    for (i = 0; i < tr->nnodes; i++) mark[i] = 0;
    stack = stk_new_ptr(tr->nnodes);
    stk_push_ptr(stack, tr);
    while ((n = stk_pop_ptr(stack)) != NULL) {
      assert((n->lchild == NULL && n->rchild == NULL) ||
             (n->lchild != NULL && n->rchild != NULL));
      if (n->lchild != NULL && mark[n->lchild->id] == 0) {
        stk_push_ptr(stack, n->rchild);
        stk_push_ptr(stack, n);
        stk_push_ptr(stack, n->lchild);
      }
      else {
        lst_push_ptr(tr->inorder, n);
        mark[n->id] = 1;
      }
    }
    stk_free(stack);
    free(mark);
  }

  return tr->inorder;
}

List *tr_postorder(TreeNode *tr) {
  if (tr->postorder == NULL) {  /* produce on demand */
    int i;
    int *mark;
    Stack *stack;
    TreeNode *n;

    tr->postorder = lst_new_ptr(tr->nnodes);
    mark = (int*)smalloc(tr->nnodes * sizeof(int));
    for (i = 0; i < tr->nnodes; i++) mark[i] = 0;
    stack = stk_new_ptr(tr->nnodes);
    stk_push_ptr(stack, tr);
    while ((n = stk_pop_ptr(stack)) != NULL) {
      assert((n->lchild == NULL && n->rchild == NULL) ||
             (n->lchild != NULL && n->rchild != NULL));
      if (n->lchild != NULL && mark[n->lchild->id] == 0) {
        stk_push_ptr(stack, n); /* order? */
        stk_push_ptr(stack, n->rchild);
        stk_push_ptr(stack, n->lchild);
      }
      else {
        lst_push_ptr(tr->postorder, n);
        mark[n->id] = 1;
      }
    }
    stk_free(stack);
    free(mark);
  }
  
  return tr->postorder;
}

/* Provide x-y coordinates for every node in a tree, given a bounding
   box defined by (x0,y0) [upper left] and (x1,y1) [lower right].  (It
   is assumed that x0 < x1 and y0 > y1).  The arrays x and y must be
   pre-allocated to dimension tree->nnodes.  After the routine is
   called, the coordinates for the node with id = i will be given by
   (x[i],y[i]).  If use_branch_lens == 1, the tree will be be laid out
   such that edges are proportional to the 'dparent' attributes of the
   TreeNodes.  If horizontal = 1, the tree will be laid out with the
   root on the left and the leaves on the right; otherwise, the root
   will be at top and the leaves at bottom.  */
void tr_layout_xy(TreeNode *tree, int x0, int y0, int x1, int y1, 
                  int *x, int *y, int use_branch_lens, int horizontal) {

  int delt_x, delt_y, i; 
  List *traversal; 
  TreeNode *n;
  double scale = 0;

  delt_x = x1 - x0; 
  delt_y = y0 - y1; 

  /* if scaling according to branch lens, we need the total height of
     the tree (in terms of branch lengths, not levels) */
  if (use_branch_lens) {
    double *total_height = (double*)smalloc(tree->nnodes * sizeof(double));
    traversal = tr_postorder(tree);
    for (i = 0; i < tree->nnodes; i++) {
      n = lst_get_ptr(traversal, i);
      if (n->lchild == NULL)
        total_height[n->id] = 0;
      else
        total_height[n->id] = 
          max(total_height[n->lchild->id] + n->lchild->dparent, 
              total_height[n->rchild->id] + n->rchild->dparent);
    }
    scale = (horizontal == 1 ? abs(delt_x)/total_height[tree->id] :
             abs(delt_y)/total_height[tree->id]);
    free(total_height);
  }

  /* set x coords (or y's if horizontal) by spacing evenly in an
     inorder traversal; if not scaling according to branch lens, can
     simultaneously set opposite coordinates, using the "height" of
     each node in the tree */
  traversal = tr_inorder(tree); 
  for (i = 0; i < tree->nnodes; i++) {  
    n = lst_get_ptr(traversal, i); 
    if (horizontal == 0) {
      x[n->id] = x0 + i*delt_x/tree->nnodes;     
      if (!use_branch_lens) y[n->id] = y1 + n->height * delt_y/tree->height;
    }
    else {
      y[n->id] = y1 + i*delt_y/tree->nnodes;     
      if (!use_branch_lens) x[n->id] = x1 - n->height * delt_x/tree->height;
    }
  } 

  /* if scaling according to branch lens, set y coords (or x's if
     horizontal) incrementally in a preorder traversal */
  if (use_branch_lens) {
    traversal = tr_preorder(tree);
    for (i = 0; i < tree->nnodes; i++) {  
      n = lst_get_ptr(traversal, i); 
      if (horizontal == 0) 
        y[n->id] = n->parent == NULL ? y0 :
          y[n->parent->id] - n->dparent*scale;
      else
        x[n->id] = n->parent == NULL ? x0 :
          x[n->parent->id] + n->dparent*scale;
    }
  }
} 

/* Print a simple postscript rendering of a tree.  Node coordinates
   are defined by tr_layout_xy (above).  Parameters specify whether to
   show branch lengths, whether to draw the tree horizontally or
   vertically, whether to scale proportionally to branch lengths, and
   whether to draw "square" (right-angled) branches or diagonal ones.
   To do: support scaling of entire tree in x or y dimension
   ... perhaps should do some of this automatically, to avoid really
   tall & skinny trees ... also may want to accept BoundingBox params
   as input */
void tr_print_ps(FILE *F, TreeNode *tree, int show_branch_lens,
                 int square_branches, int use_branch_lens, 
                 int horizontal_layout) {
  int i, xoffset, yoffset;
  int *x, *y;
  List *traversal;

  x = (int*)smalloc(tree->nnodes * sizeof(int));
  y = (int*)smalloc(tree->nnodes * sizeof(int));
  tr_layout_xy(tree, TL_X, TL_Y, BR_X, BR_Y, x, y, 
               use_branch_lens, horizontal_layout);

  /* print header */
  fprintf(F, "%%!\n\
1 setlinecap 1 setlinejoin 1 setlinewidth 0 setgray\n\
/basefont /Times-Roman findfont 12 scalefont def\n\
50 50 translate\n\
basefont setfont\n");

  traversal = tr_postorder(tree);
  for (i = 0; i < lst_size(traversal); i++) {
    TreeNode *n = lst_get_ptr(traversal, i);
    if (n->lchild == NULL) {    /* leaf */

      /* offsets for labels; should parameterize this a bit better,
         perhaps accepting font size as a function parameter */
      if (horizontal_layout) {
        xoffset = 10;
        yoffset = -6;
      }
      else {
        xoffset = -3 * (n->name != NULL ? strlen(n->name) : 0);
        yoffset = -18;
      }

      /* draw leaf label */
      fprintf(F, "%d %d moveto\n(%s) show\n", x[n->id]+xoffset, 
              y[n->id]+yoffset, 
              n->name != NULL ? n->name : "");
    }
    else {                      /* internal node */
      /* draw branches from node to each of its children; this is all
         a bit messy but it will have to do for now ... */
      if (square_branches) {
        if (horizontal_layout) {
          fprintf(F, "%d %d moveto\n%d %d lineto stroke\n", x[n->lchild->id], 
                  y[n->lchild->id], x[n->id], y[n->lchild->id]);
          fprintf(F, "%d %d moveto\n%d %d lineto stroke\n", x[n->id], 
                  y[n->lchild->id], x[n->id], y[n->id]);
          fprintf(F, "%d %d moveto\n%d %d lineto stroke\n", x[n->rchild->id], 
                  y[n->rchild->id], x[n->id], y[n->rchild->id]);
          fprintf(F, "%d %d moveto\n%d %d lineto stroke\n", x[n->id], 
                  y[n->rchild->id], x[n->id], y[n->id]);

          if (show_branch_lens) {
            fprintf(F, "%d %d moveto\n(%.4f) show\n", 
                    (x[n->lchild->id]-x[n->id])/2 + x[n->id], 
                    y[n->lchild->id] + 6, n->lchild->dparent);
            fprintf(F, "%d %d moveto\n(%.4f) show\n", 
                    (x[n->rchild->id]-x[n->id])/2 + x[n->id], 
                    y[n->rchild->id] + 6, n->rchild->dparent);
          }
        }
        else {
          fprintf(F, "%d %d moveto\n%d %d lineto stroke\n", x[n->lchild->id], 
                  y[n->lchild->id], x[n->lchild->id], y[n->id]);
          fprintf(F, "%d %d moveto\n%d %d lineto stroke\n", x[n->lchild->id], 
                  y[n->id], x[n->id], y[n->id]);
          fprintf(F, "%d %d moveto\n%d %d lineto stroke\n", x[n->rchild->id], 
                  y[n->rchild->id], x[n->rchild->id], y[n->id]);
          fprintf(F, "%d %d moveto\n%d %d lineto stroke\n", x[n->rchild->id], 
                  y[n->id], x[n->id], y[n->id]);

          if (show_branch_lens) {
            fprintf(F, "%d %d moveto\n(%.4f) show\n", 
                    x[n->lchild->id] + 6,
                    (y[n->id]-y[n->lchild->id])/2 + y[n->lchild->id],
                    n->lchild->dparent);
            fprintf(F, "%d %d moveto\n(%.4f) show\n", 
                    x[n->rchild->id] + 6,
                    (y[n->id]-y[n->rchild->id])/2 + y[n->rchild->id],
                    n->rchild->dparent);
          }
        }
      }
      else {                    /* diag branches */
        fprintf(F, "%d %d moveto\n%d %d lineto stroke\n", x[n->lchild->id], 
                y[n->lchild->id], x[n->id], y[n->id]);
        fprintf(F, "%d %d moveto\n%d %d lineto stroke\n", x[n->rchild->id], 
                y[n->rchild->id], x[n->id], y[n->id]);

        if (show_branch_lens) {
          if (horizontal_layout) {
            fprintf(F, "%d %d moveto\n(%.4f) show\n", 
                    (x[n->lchild->id]-x[n->id])/2 + x[n->id], 
                    (y[n->id]-y[n->lchild->id])/2 + y[n->lchild->id] + 3, 
                    n->lchild->dparent);
            fprintf(F, "%d %d moveto\n(%.4f) show\n", 
                    (x[n->rchild->id]-x[n->id])/2 + x[n->id], 
                    (y[n->rchild->id]-y[n->id])/2 + y[n->id] - 12, 
                    n->rchild->dparent);
          }
          else {
            fprintf(F, "%d %d moveto\n(%.4f) show\n", 
                    (x[n->lchild->id]-x[n->id])/2 + x[n->id] + 6,
                    (y[n->id]-y[n->lchild->id])/2 + y[n->lchild->id],
                    n->lchild->dparent);
            fprintf(F, "%d %d moveto\n(%.4f) show\n", 
                    (x[n->id]-x[n->rchild->id])/2 + x[n->rchild->id] + 6,
                    (y[n->id]-y[n->rchild->id])/2 + y[n->rchild->id],
                    n->rchild->dparent);
          }
        }
      }
    }
  }
  fprintf(F, "showpage\n");     /* complete PS file */

  free(x);
  free(y);
}

/* return sum of lengths at all edges */
double tr_total_len(TreeNode *t) {
  double retval = 0;
  int i;
  for (i = 0; i < t->nnodes; i++) {
    TreeNode *n = lst_get_ptr(t->nodes, i);
    if (n->parent != NULL) 
      retval += n->dparent;
  }
  return retval;
}

/* return node having specified name or NULL if none found.  */
TreeNode *tr_get_node(TreeNode *t, char *name) {
  int i;
  for (i = 0; i < t->nnodes; i++) {
    TreeNode *n = lst_get_ptr(t->nodes, i);
    if (n->name != NULL && !strcmp(n->name, name))
      return n;
  }
  return NULL;
}

void tr_number_leaves(TreeNode *t, char **names, int nnames) {
  int i, j;
  char *endptr;
  for (i = 0; i < t->nnodes; i++) {
    TreeNode *n = lst_get_ptr(t->nodes, i);
    if (strlen(n->name) > 0) {
      strtol(n->name, &endptr, 0);
      if (*endptr != '\0') {    /* doesn't parse as int -- assume
                                   name */
        for (j = 0; j < nnames; j++) {
          if (strcmp(names[j], n->name) == 0) {
            sprintf(n->name, "%d", j+1);
            break;
          }
        }
        if (j == nnames) die("ERROR: no match for name '%s' given in tree topology.\n", n->name);
      }
    }
  }
}
