/* $Id: trees.c,v 1.16 2004-08-11 20:47:05 acs Exp $ 
   Written by Adam Siepel, 2002
   Copyright 2002, Adam Siepel, University of California */

/** \file trees.c 
  Functions for manipulating phylogenetic trees.  Includes functions
  for reading, writing, traversing, and printing.  Trees are
  represented as rooted binary trees with non-negative real branch
  lengths.  
  \ingroup phylo
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


static int idcounter = 0;
/* NOTE: when tree is parsed from Newick file, node ids are assigned
   sequentially in a preorder traversal.  Some useful properties
   result.  For example, if two nodes u and v are such that v->id >
   u->id, then the first ancestor a of v s.t. a->id < u->id is the LCA
   of u and v */

/* coords for postscript printing */
/* top-left x */
#define TL_X 15                 

/* top-left y */
#define TL_Y 10                 
                                
/* bottom-right x */
#define BR_X 500                
                                
/* bottom-right y */
#define BR_Y 700                
                                

/** Parse a tree from a file in Newick (New Hampshire) format */
TreeNode *tr_new_from_file(FILE *f) { 
  String *s = str_new(STR_VERY_LONG_LEN);
  TreeNode *retval;

  str_slurp(s, f);

  str_double_trim(s);
  if (s->chars[0] != '(')
    die("ERROR: This doesn't look like a tree (Newick format): \"%s\".\n", 
        s->chars);

  if (s->chars[s->length-1] == ';') 
    s->chars[--s->length] = '\0';

  retval = tr_new_from_string(s->chars);
  str_free(s);
  return retval;
}

/** Parse a Newick-formatted tree from a character string */
TreeNode *tr_new_from_string(char *treestr) { 
  TreeNode *root, *node, *newnode;
  int i, in_distance = 0, len = strlen(treestr), nopen_parens = 0,
    nclose_parens = 0, already_allowed = FALSE;
  char c;
  char *currentname = NULL;
  String *diststr = str_new(STR_SHORT_LEN);

  tr_reset_id();
  root = tr_new_node(); root->nnodes = 1;
  node = root;
  for (i = 0; i < len; i++) {
    c = treestr[i];

    if (in_distance) {
      if (c != '(' && c != ',' && c != ')' && c != ':') {
        str_append_char(diststr, c);
        continue;
      }
      else {
        if (str_as_dbl(diststr, &node->dparent) != 0)
          die("ERROR: Can't parse distance in tree (\"%s\").\n", 
              diststr->chars);
      }
      in_distance = 0;
    }

    if (c == '(') {
      tr_add_child(node, newnode = tr_new_node());
      node = newnode;
      currentname = newnode->name;
      root->nnodes++;
      nopen_parens++;
    }
    else if (c == ',') {
      if (node->parent == NULL)
        die("ERROR: invalid binary tree (check parens).\n");
      if (node->parent->lchild != NULL && node->parent->rchild != NULL){
        if (node->parent == root && !already_allowed)
          already_allowed = TRUE;
        else
          die("ERROR (tree parser): invalid binary tree (too many children)\n");
                                /* we'll prohibit multinary
                                   branchings, except that we'll allow
                                   a single trinary branch immediately
                                   below the root (common with
                                   reversible models) */
      }

      tr_add_child(node->parent, newnode = tr_new_node());
      node = newnode;
      currentname = node->name;
      root->nnodes++;
    }
    else if (c == ')') {
      node = node->parent;
      currentname = NULL;
      nclose_parens++;
    }
    else if (c == ':') {
      str_clear(diststr);
      in_distance = 1;
    }
    else if (currentname != NULL) {
      if (!isspace(c) || currentname[0] != '\0') /* avoid leading spaces */
        strncat(currentname, &c, 1);
    }
  }

  if (nopen_parens != nclose_parens)
    die("ERROR: mismatching parens in tree.\n");

  tr_set_nnodes(root);
  str_free(diststr);
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
TreeNode *tr_new_node() {
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
void tr_add_child(TreeNode *parent, TreeNode *child) {
  if (parent->lchild == NULL) {
    parent->lchild = child;
  }
  else if (parent->rchild == NULL) {
    parent->rchild = child;
  }
  else {
    /* add intermediate node to accommodate extra child */
    TreeNode *tmp = tr_new_node();
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

/** Print tree in New Hampshire format. */
void tr_print(FILE* f, TreeNode *root, int show_branch_lengths) {
  /* It's simplest to do this recursively. */
  tr_print_recur(f, root, show_branch_lengths);
  fprintf(f, ";\n");
}

/* Recursive subroutine used by print_tree */
void tr_print_recur(FILE* f, TreeNode *n, int show_branch_lengths) {

  assert((n->lchild == NULL && n->rchild == NULL) || 
	 (n->lchild != NULL && n->rchild != NULL));

  if (n->lchild != NULL) {
    fprintf(f, "(");
    tr_print_recur(f, n->lchild, show_branch_lengths);
    fprintf(f, ",");
    tr_print_recur(f, n->rchild, show_branch_lengths);
    fprintf(f, ")");
  }
  else {
    fprintf(f, "%s", n->name);
  }

  if (show_branch_lengths && n->parent != NULL)
    fprintf(f, ":%f", n->dparent);
}

/** Free memory for tree */
void tr_free(TreeNode *tree) {
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

/** Copy tree */
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
      tr_add_child(ncpy, lcpy);
      stk_push_ptr(stack, n->lchild);
      stk_push_ptr(cpystack, lcpy);
    }
    if (n->rchild != NULL) {
      rcpy = stk_pop_ptr(nodes);
      tr_node_cpy(rcpy, n->rchild);
      tr_add_child(ncpy, rcpy);
      stk_push_ptr(stack, n->rchild);
      stk_push_ptr(cpystack, rcpy);
    }
  }

  stk_free(stack);
  stk_free(cpystack);
  stk_free(nodes);
}

/** Create a new tree that's a copy of another one */
TreeNode *tr_create_copy(TreeNode *src) {
  Stack *stack, *cpystack;
  TreeNode *n, *ncpy, *lcpy, *rcpy, *dest;
  
  tr_reset_id();

  stack = stk_new_ptr(src->nnodes);
  cpystack = stk_new_ptr(src->nnodes);
  stk_push_ptr(stack, src); 
  dest = tr_new_node();
  tr_node_cpy(dest, src);
  stk_push_ptr(cpystack, dest);
  while ((n = stk_pop_ptr(stack)) != NULL) {
    ncpy = stk_pop_ptr(cpystack);
    if (n->lchild != NULL) {
      lcpy = tr_new_node();
      tr_node_cpy(lcpy, n->lchild);
      tr_add_child(ncpy, lcpy);
      stk_push_ptr(stack, n->lchild);
      stk_push_ptr(cpystack, lcpy);
    }
    if (n->rchild != NULL) {
      rcpy = tr_new_node();
      tr_node_cpy(rcpy, n->rchild);
      tr_add_child(ncpy, rcpy);
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
  dest->id = src->id;
  strcpy(dest->name, src->name); 
  dest->dparent = src->dparent;
  /* don't copy data, nnodes, height, preorder, inorder, postorder */
}


/** Print tree in Newick format.  This version imposes an
   ordering on the leaves (useful when comparing trees that have been
   rearranged).  At every internal node, we store the name of the leaf
   beneath it that comes first alphanumerically.  When recursively
   printing the tree, at each internal node, we call its children in
   the order of these names.  */
void tr_print_ordered(FILE* f, TreeNode *root, int show_branch_lengths) {
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

  tr_print_ordered_recur(f, root, left_right, show_branch_lengths);
  fprintf(f, ";\n");
  
  stk_free(stack);
  free(left_right);
  free(mark);
  free(names);
}

/* Recursive subroutine for tr_print_ordered */
void tr_print_ordered_recur(FILE* f, TreeNode *n, int *left_right,
                            int show_branch_lengths) {

  assert((n->lchild == NULL && n->rchild == NULL) || 
	 (n->lchild != NULL && n->rchild != NULL));

  if (n->lchild != NULL) {
    fprintf(f, "(");
    if (left_right[n->id]) {
      tr_print_ordered_recur(f, n->lchild, left_right, show_branch_lengths);
      fprintf(f, ",");
      tr_print_ordered_recur(f, n->rchild, left_right, show_branch_lengths);
    }
    else {
      tr_print_ordered_recur(f, n->rchild, left_right, show_branch_lengths);
      fprintf(f, ",");
      tr_print_ordered_recur(f, n->lchild, left_right, show_branch_lengths);
    }
    fprintf(f, ")");
  }
  else {
    fprintf(f, "%s", n->name);
  }

  if (show_branch_lengths)
    fprintf(f, ":%f", n->dparent);
}

/** Obtain a list representing a preorder traversal of the tree.
    Allows for simplified loops (no stacks!) and improved
    efficiency */
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

/** Obtain a list representing an in-order traversal of the tree.
    Allows for simplified loops (no stacks!) and improved
    efficiency */
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

/** Obtain a list representing a postorder traversal of the tree.
    Allows for simplified loops (no stacks!) and improved
    efficiency */
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

/** Provide x-y coordinates for layout. */
void tr_layout_xy(TreeNode *tree, 
                  int x0,       /**< Upper left x bound */
                  int y0,       /**< Upper left y bound  */
                  int x1,       /**< Lower right x bound */
                  int y1,       /**< Lower right y bound  */
                  int *x,       /**< On return, will contain
                                   x-coordinates for nodes, in order
                                   of tree->nodes.  Must be
                                   preallocated. */
                  int *y,       /**< On return, will contain
                                   y-coordinates for nodes, in order
                                   of tree->nodes.  Must be
                                   preallocated. */
                  int use_branch_lens, 
                                /**< If TRUE, tree will be laid out
                                   such that edges are proportional to
                                   branch lengths (dparent
                                   attributes) */
                  int horizontal
                                /**< If TRUE, tree will be laid out
                                   with root on left and leaves on
                                   right; otherwise, root will be at
                                   top and leaves at bottom */
                  ) {

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

/** Print a (very basic!) postscript rendering of a tree. */
void tr_print_ps(FILE *F,       /**< Destination file */
                 TreeNode *tree, 
                                /**< Tree root */
                 int show_branch_lens,
                                /**< Whether to print branch lengths
                                   by edges */
                 int square_branches,
                                /**< If TRUE, branches will be
                                   right-angled, otherwise will be
                                   diagonal */
                 int use_branch_lens, 
                                /**< If TRUE, tree will be laid out
                                   such that edges are proportional to
                                   branch lengths (dparent
                                   attributes) */
                 int horizontal_layout
                                /**< If TRUE, tree will be laid out
                                   with root on left and leaves on
                                   right; otherwise, root will be at
                                   top and leaves at bottom */
                 ) {
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

/** Compute and return sum of lengths at all edges */
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

/** Compute and return sum of lengths of edges in subtree below
    given node */
double tr_total_len_subtree(TreeNode *sub_root) {
  TreeNode *n;
  Stack *stack = stk_new_ptr(sub_root->nnodes);
  double retval = 0;
  stk_push_ptr(stack, sub_root->lchild);
  stk_push_ptr(stack, sub_root->rchild);
  while ((n = stk_pop_ptr(stack)) != NULL) {
    if (n->lchild != NULL) {
      stk_push_ptr(stack, n->rchild);
      stk_push_ptr(stack, n->lchild);
    }
    retval += n->dparent;
  }
  stk_free(stack);
  return retval;
}

/** Return node having specified name or NULL if none found.  */
TreeNode *tr_get_node(TreeNode *t, char *name) {
  int i;
  for (i = 0; i < t->nnodes; i++) {
    TreeNode *n = lst_get_ptr(t->nodes, i);
    if (n->name[0] != '\0' && !strcmp(n->name, name))
      return n;
  }
  return NULL;
}

/** Scale all branch lengths by constant factor. */
void tr_scale(TreeNode *t, double scale_const) {
  int i;
  for (i = 0; i < t->nnodes; i++) {
    TreeNode *n = lst_get_ptr(t->nodes, i);
    if (n->parent != NULL) 
      n->dparent *= scale_const;
  }
}

/** Prune away all leaves whose names are in (or not in) the specified
    list.  Nodes will be removed and branches combined (branch lengths
    added) to restore as a proper binary tree.  */
void tr_prune(TreeNode **t,     /**< Tree to prune (may be altered
                                   because root can change) */
              List *names,      /**< List of names.  On return, will
                                   contain list of names of leaves
                                   that were pruned away.  */
              int all_but       /**< if FALSE, prune leaves *in*
                                   'names'; if TRUE, prune leaves *not
                                   in* 'names'  */
              ) {

  TreeNode *n;
  int i, new_nnodes = (*t)->nnodes;
  int *is_leaf;
  List *traversal, *pruned_leaves = lst_new_ptr((*t)->nnodes / 2);

  /* first identify original leaves; will need to distinguish them
     from leaves that are created by pruning */
  is_leaf = smalloc((*t)->nnodes * sizeof(int));
  for (i = 0; i < (*t)->nnodes; i++) {
    n = lst_get_ptr((*t)->nodes, i);
    is_leaf[i] = (n->lchild == NULL && n->rchild == NULL);
  }

  /* get rid of nodes, preorder, and inorder lists (do now because
     root may change) */
  if ((*t)->nodes != NULL) { lst_free((*t)->nodes); (*t)->nodes = NULL; }
  if ((*t)->preorder != NULL) { lst_free((*t)->preorder); (*t)->preorder = NULL; }
  if ((*t)->inorder != NULL) { lst_free((*t)->inorder); (*t)->inorder = NULL; }

  /* remove nodes and combine branches in postorder traversal */
  traversal = tr_postorder(*t);
  for (i = 0; i < lst_size(traversal); i++) {
    TreeNode *n = lst_get_ptr(traversal, i);
    if (n->lchild == NULL && n->rchild == NULL){ /* missing both children */
      String *s;
      int prune;
      
      if (!is_leaf[n->id])      /* if not originally a leaf, must be pruned */
        prune = TRUE;
      else {
        s = str_new_charstr(n->name);
        prune = str_in_list(s, names);
        if (all_but) prune = !prune;

        if (prune) lst_push_ptr(pruned_leaves, s);
        else str_free(s);
      }

      if (prune) {
        if (n->parent == NULL) 
          *t = NULL;            /* entire tree has been pruned away! */
        else {
          if (n == n->parent->lchild) n->parent->lchild = NULL;
          else n->parent->rchild = NULL;
        }
        free(n);
        new_nnodes--;
      }
    }
    else if (n->lchild == NULL) { /* missing left child only */
      if (n->parent == NULL) {
        assert(n == *t);        /* n must be root */
        n->rchild->parent = NULL; /* redefine root */
        *t = n->rchild;
        (*t)->dparent = 0;
      }
      else {                    /* mid-level node; remove and combine
                                   branch lengths */
        n->rchild->parent = n->parent;
        n->rchild->dparent += n->dparent;
        if (n == n->parent->lchild) n->parent->lchild = n->rchild;
        else n->parent->rchild = n->rchild;
      }
      free(n);
      new_nnodes--;
    }
    else if (n->rchild == NULL) { /* missing right child only */
      if (n->parent == NULL) {
        assert(n == *t);        /* n must be root */
        n->lchild->parent = NULL; /* redefine root */
        *t = n->lchild;         
        (*t)->dparent = 0;
      }
      else {                    /* mid-level node; remove and combine
                                   branch lengths */
        n->lchild->parent = n->parent;
        n->lchild->dparent += n->dparent;
        if (n == n->parent->lchild) n->parent->lchild = n->lchild;
        else n->parent->rchild = n->lchild;
      }
      free(n);
      new_nnodes--;
    }
  }

  /* finally, free postorder list */
  lst_free(traversal); 
  if (*t != NULL && (*t)->postorder != NULL) (*t)->postorder = NULL;

  /* reset ids, nodes, nnodes, heights */
  if (*t != NULL) {
    (*t)->nnodes = new_nnodes;
    traversal = tr_preorder(*t);
    for (i = 0; i < lst_size(traversal); i++) {
      n = lst_get_ptr(traversal, i);
      n->id = i;
      if (n != *t) n->nnodes = -1;
    }
    tr_set_nnodes(*t);
  }

  lst_free_strings(names);
  lst_clear(names);
  for (i = 0; i < lst_size(pruned_leaves); i++) 
    lst_push_ptr(names, lst_get_ptr(pruned_leaves, i));

  lst_free(pruned_leaves);
  free(is_leaf);
}

/** Return the LCA of the given species.  Assumes ids are numbered in
    preorder (a node's parent always has a smaller id than it does and
    left descendants have smaller ids than right descendants). */
TreeNode *tr_lca(TreeNode *tree, List *names) {
  int i, min = tree->nnodes, max = -1, idx;
  String *tmpstr = str_new(STR_MED_LEN);
  TreeNode *n;
  int *found = smalloc(lst_size(names) * sizeof(int));

  for (i = 0; i < lst_size(names); i++) found[i] = FALSE;

  for (i = 0; i < tree->nnodes; i++) {
    n = lst_get_ptr(tree->nodes, i);
    if (n->lchild == NULL && n->rchild == NULL && n->name[0] != '\0') {
      str_cpy_charstr(tmpstr, n->name);
      if (str_in_list_idx(tmpstr, names, &idx)) {
        found[idx] = TRUE;
        if (n->id < min) min = n->id;
        if (n->id > max) max = n->id;
      }
    }
  }

  for (i = 0; i < lst_size(names); i++)
    if (!found[i])
      die("ERROR: species name not found in tr_lca ('%s')\n",
          ((String*)lst_get_ptr(names, i))->chars);

  /* now the LCA must be the first ancestor of the node with the max
     id that has an id smaller than the min id */
  for (n = lst_get_ptr(tree->nodes, max); n->id > min; n = n->parent);

  str_free(tmpstr);
  free(found);
  return n;
}

/** Given two trees, one of which is a subtree of the other, create a
    hybrid tree composed of the smaller tree and a scaled version of
    the larger tree.  First, a copy of the larger tree will be created
    and scaled such that the total branch length in the subtree in
    question is equal to the total branch length of the smaller tree.
    Then, (a copy of) the smaller tree will be used in place of the
    overlapping subtree in the larger tree.  This function can be used
    to extrapolate from a small phylogeny for which accurate branch
    length estimation is possible (e.g., including eutherian mammals)
    to a larger phylogeny for which approximate branch length
    proportions are available, but absolute branch length estimates
    are not (e.g., including more distant vertebrates). */
TreeNode *tr_hybrid(TreeNode *sub, TreeNode *super) {
  TreeNode *retval, *n, *lca, *sub_copy;
  int i;
  double lfrac, sum;
  List *names = lst_new_ptr((sub->nnodes + 1) / 2);

  if (sub->nnodes < 3)
    die("ERROR: subtree must have at least two leaves in tr_hybrid.\n");

  /* copy supertree then find LCA corresponding to subtree */
  retval = tr_create_copy(super);
  for (i = 0; i < sub->nnodes; i++) {
    n = lst_get_ptr(sub->nodes, i);
    if (n->lchild == NULL && n->rchild == NULL)
      lst_push_ptr(names, str_new_charstr(n->name));
  }
  lca = tr_lca(retval, names);
  lst_free_strings(names);
  lst_free(names);

  sub_copy = tr_create_copy(sub);

  if (lca == super) {           /* rule out trivial case -- trees equal */
    tr_free(retval);
    return sub_copy;
  }

  /* scale supertree so that overlapping portions have equal total length */
  tr_scale(retval, tr_total_len(sub_copy) / tr_total_len_subtree(lca));

  /* now recombine */
  sub_copy->parent = lca->parent;
  sub_copy->dparent = lca->dparent;
  if (lca == lca->parent->lchild) lca->parent->lchild = sub_copy;
  else lca->parent->rchild = sub_copy;

  /* also ensure that the subtree is rooted with proportions equal to
     those of the original supertree (usually hard to get root right
     in subtree) */
  lfrac = lca->lchild->dparent / (lca->lchild->dparent + lca->rchild->dparent);
  sum = sub_copy->lchild->dparent + sub_copy->rchild->dparent;
  sub_copy->lchild->dparent = lfrac * sum;
  sub_copy->rchild->dparent = sum - sub_copy->lchild->dparent;

  tr_free(lca);                 /* this works recursively */

  return retval;
}

/** Partition leaves of tree at (branch above) given node.  All
    descendant leaves of 'sub' will be added to 'inside' list and all
    non-descendants will be added to 'outside' list.  Lists must be
    pre-allocated. */
void tr_partition_leaves(TreeNode *tree, TreeNode *sub, List *inside, 
                         List *outside) {
  int i;
  TreeNode *n;
  int *mark = smalloc(tree->nnodes * sizeof(int));
  Stack *stack = stk_new_ptr(sub->nnodes);

  for (i = 0; i < tree->nnodes; i++) mark[i] = FALSE;

  lst_clear(inside);
  lst_clear(outside);
  stk_push_ptr(stack, sub);
  while ((n = stk_pop_ptr(stack)) != NULL) {
    if (n->lchild == NULL) {
      lst_push_ptr(inside, n);
      mark[n->id] = TRUE;
    }
    else {
      stk_push_ptr(stack, n->rchild);
      stk_push_ptr(stack, n->lchild);
    }
  }
  for (i = 0; i < tree->nnodes; i++) {
    n = lst_get_ptr(tree->nodes, i);
    if (n->lchild == NULL && !mark[n->id])
      lst_push_ptr(outside, n);
  }
  stk_free(stack);
  free(mark);
}
