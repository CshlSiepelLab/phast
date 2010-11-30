/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/*****************************************************
rph_trees.c
The RPHAST handles to functions dealing with trees from
the phast package.

Melissa Hubisz
Last updated: 1/5/2010
*****************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <ctype.h>
#include <sufficient_stats.h>
#include <local_alignment.h>
#include <trees.h>
#include <misc.h>
#include <rph_util.h>

#include <Rdefines.h>


void rph_tree_protect(TreeNode *tr) {
  TreeNode *n;
  int i;
  if (tr == NULL) return;
  rph_mem_protect(tr);
  for (i=0; i < tr->nnodes; i++) {
    n = (TreeNode*)lst_get_ptr(tr->nodes, i);
    rph_mem_protect(n);
    if (n->nodes != NULL) 
      rph_lst_protect(n->nodes);
    if (n->preorder != NULL)
      rph_lst_protect(n->preorder);
    if (n->postorder != NULL)
      rph_lst_protect(n->postorder);
    if (n->inorder != NULL)
      rph_lst_protect(n->inorder);
    if (n->label != NULL)
      rph_mem_protect(n->label);
  }
}


TreeNode* rph_tree_new(SEXP treeStr) {
  TreeNode *tree = tr_new_from_string(CHARACTER_VALUE(treeStr));
  return tree;
}


/* read in a tree from a file.  Return character string 
   representing tree */
SEXP rph_tree_read(SEXP filename) {
  FILE *infile;
  char c;
  SEXP result;
  char *currStr, **strvec;
  int i, pos=0, currLen=10000, numparen=0, numtrees_alloc=1000, numtrees=0;
  TreeNode *tempTree;

  infile = fopen_fname(CHARACTER_VALUE(filename), "r");
  currStr = smalloc((currLen+2)*sizeof(char));
  strvec = smalloc(numtrees_alloc*sizeof(char*));
  while (1) {
    pos=0;
    numparen=0;
    while (';'!=(c=fgetc(infile))) {
      if (c==EOF) {
	if (pos==0) break;
	die("unexpected EOF in tree file.  Trees should be terminated by \";\"");
      }
      if (isspace(c)) continue;
      if (c=='(') numparen++;
      if (c==')') numparen--;
      if (numparen < 0) die("bad tree in tree file");
      if (pos==currLen) {
	currLen += 10000;
	currStr = srealloc(currStr, (currLen+2)*sizeof(char));
      }
      currStr[pos++] = c;
    }
    if (pos > 0) {
      if (numparen != 0) die("unbalanced parenthesis in tree file");
      currStr[pos++]=';';
      currStr[pos]='\0';
      if (numtrees == numtrees_alloc) {
	numtrees_alloc += 1000;
	strvec = srealloc(strvec, numtrees_alloc*sizeof(char*));
      }
      strvec[numtrees] = smalloc((strlen(currStr)+1)*sizeof(char));
      strcpy(strvec[numtrees], currStr);

      //check to make sure phast can read this tree
      tempTree = tr_new_from_string(strvec[numtrees]);
      numtrees++;
    }
    else break;
  }
  fclose(infile);
  PROTECT(result = NEW_CHARACTER(numtrees));
  for (i=0; i<numtrees; i++) 
    SET_STRING_ELT(result, i, mkChar(strvec[i]));
  UNPROTECT(1);
  return result;
}

SEXP rph_tree_numnodes(SEXP tree) {
  TreeNode *tr = rph_tree_new(tree);
  SEXP result;
  int *resultP;
  PROTECT(result = NEW_INTEGER(1));
  resultP = INTEGER_POINTER(result);
  resultP[0] = tr->nnodes;
  UNPROTECT(1);
  return result;
}


SEXP rph_tree_prune(SEXP treeStr, SEXP seqsP, SEXP allButP) {
  TreeNode *tr = rph_tree_new(treeStr);
  List *names = lst_new_ptr(LENGTH(seqsP));
  String *tempStr;
  char *temp;
  int i;
  SEXP result;
  for (i=0; i<LENGTH(seqsP); i++) {
    tempStr = str_new_charstr(CHAR(STRING_ELT(seqsP, i)));
    lst_push_ptr(names, tempStr);
  }
  tr_prune(&tr, names, INTEGER_VALUE(allButP));
  temp = tr_to_string(tr, 1);
  PROTECT(result = NEW_CHARACTER(1));
  SET_STRING_ELT(result, 0, mkChar(temp));
  UNPROTECT(1);
  return result;
}


SEXP rph_tree_name_ancestors(SEXP treeStr) {
  TreeNode *tr = rph_tree_new(treeStr);
  char *newTreeStr;
  SEXP result;
  tr_name_ancestors(tr);
  newTreeStr = tr_to_string(tr, 1);
  PROTECT(result = NEW_CHARACTER(1));
  SET_STRING_ELT(result, 0, mkChar(newTreeStr));
  UNPROTECT(1);
  return result;
}
  
  

SEXP rph_tree_subtree(SEXP treeStr, SEXP nodeStr) {
  TreeNode *tr = rph_tree_new(treeStr);
  TreeNode *n;
  char *newTreeStr;
  SEXP result;
  n = tr_get_node(tr, CHARACTER_VALUE(nodeStr));
  if (n == NULL) {
    tr_name_ancestors(tr);
    n = tr_get_node(tr, CHARACTER_VALUE(nodeStr));
    if (n == NULL)
      die("No node named %s", CHARACTER_VALUE(nodeStr));
  }
  tr_prune_supertree(&tr, n);
  newTreeStr = tr_to_string(tr, 1);
  PROTECT(result = NEW_CHARACTER(1));
  SET_STRING_ELT(result, 0, mkChar(newTreeStr));
  UNPROTECT(1);
  return result;
}


SEXP rph_tree_supertree(SEXP treeStr, SEXP nodeStr) {
  TreeNode *tr = rph_tree_new(treeStr);
  TreeNode *n;
  char *newTreeStr;
  SEXP result;

  n = tr_get_node(tr, CHARACTER_VALUE(nodeStr));
  if (n == NULL) {
    tr_name_ancestors(tr);
    n = tr_get_node(tr, CHARACTER_VALUE(nodeStr));
    if (n == NULL)
      die("No node named %s", CHARACTER_VALUE(nodeStr));
  }
  tr_prune_subtree(&tr, n);
  newTreeStr = tr_to_string(tr, 1);
  PROTECT(result = NEW_CHARACTER(1));
  SET_STRING_ELT(result, 0, mkChar(newTreeStr));
  UNPROTECT(1);
  return result;
}


SEXP rph_tree_scale(SEXP treeStr, SEXP scaleP, SEXP nodeStr,
		    SEXP includeLeadingP) {
  TreeNode *tr = rph_tree_new(treeStr);
  double scale = NUMERIC_VALUE(scaleP);
  char *newTreeStr;
  SEXP result;

  if (nodeStr != R_NilValue) {
    TreeNode *n;
    int includeLeading=LOGICAL_VALUE(includeLeadingP);
    n = tr_get_node(tr, CHARACTER_VALUE(nodeStr));
    if (n == NULL) {
      tr_name_ancestors(tr);
      n = tr_get_node(tr, CHARACTER_VALUE(nodeStr));
      if (n == NULL) 
	die("No node named %s in %s\n", CHARACTER_VALUE(nodeStr),
	    CHARACTER_VALUE(treeStr));
    }
    tr_scale_subtree(tr, n, scale, includeLeading);
  } 
  else tr_scale(tr, scale);
  newTreeStr = tr_to_string(tr, 1);
  PROTECT(result = NEW_CHARACTER(1));
  SET_STRING_ELT(result, 0, mkChar(newTreeStr));
  UNPROTECT(1);
  return result;    
}


SEXP rph_tree_rename(SEXP treeVec, SEXP oldNamesP, SEXP newNamesP) {
  int i, numtree = LENGTH(treeVec), treeIdx;
  TreeNode *tr, *n;
  SEXP result;
  Hashtable *hash = hsh_new(20);
  char *str;
  
  for (i=0; i<LENGTH(oldNamesP); i++) {
    str = smalloc((strlen(CHAR(STRING_ELT(newNamesP, i)))+1)*sizeof(char));
    strcpy(str, CHAR(STRING_ELT(newNamesP, i)));
    hsh_put(hash, CHAR(STRING_ELT(oldNamesP, i)), str);
  }

  PROTECT(result = NEW_CHARACTER(numtree));
  for (treeIdx=0; treeIdx < numtree; treeIdx++) {
    tr = rph_tree_new(STRING_ELT(treeVec, treeIdx));
    //    tr = tr_new_from_string(CHAR(STRING_ELT(treeVec, treeIdx)));
    for (i=0; i<tr->nnodes; i++) {
      n = lst_get_ptr(tr->nodes, i);
      if (n->name != NULL && n->name[0] != '\0' &&
	  (str = hsh_get(hash, n->name)) != (char*)-1)
	strcpy(n->name, str);
    }
    str = tr_to_string(tr, 1);
    SET_STRING_ELT(result, treeIdx, mkChar(str));
  }
  UNPROTECT(1);
  return result;
}


SEXP rph_tree_nodeName(SEXP treeP, SEXP idP) {
  TreeNode *tr, *n;
  int id;
  SEXP result;
  
  if (idP == R_NilValue || treeP == R_NilValue) return R_NilValue;
  id = INTEGER_VALUE(idP);
  tr = rph_tree_new(treeP);
  n = (TreeNode*)lst_get_ptr(tr->nodes, id);
  if (id != n->id) die("id-mixup in tree");
  PROTECT(result = NEW_CHARACTER(1));
  SET_STRING_ELT(result, 0, mkChar(n->name));
  UNPROTECT(1);
  return result;
}


SEXP rph_tree_isNode(SEXP treeP, SEXP nodeName) {
  TreeNode *tr, *n;
  SEXP result;
  int *resultP, i;
  tr = rph_tree_new(treeP);
  for (i=0; i<tr->nnodes; i++) {
    n = (TreeNode*)lst_get_ptr(tr->nodes, i);
    if (strcmp(n->name, CHARACTER_VALUE(nodeName))==0)
      break;
  }
  PROTECT(result = NEW_LOGICAL(1));
  resultP = LOGICAL_POINTER(result);
  resultP[0] = (i < tr->nnodes);
  UNPROTECT(1);
  return result;
}


SEXP rph_tree_branchlen(SEXP treeP) {
  TreeNode *tr = rph_tree_new(treeP);
  SEXP rv;

  PROTECT(rv = NEW_NUMERIC(1));
  REAL(rv)[0] = tr_total_len(tr);
  UNPROTECT(1);
  return rv;
}


SEXP rph_tree_depth(SEXP treeP, SEXP nodeP) {
  TreeNode *tr = rph_tree_new(treeP), *node;
  SEXP rv;
  
  node = tr_get_node(tr, CHARACTER_VALUE(nodeP));
  if (node == NULL) 
    die("no node named %s", CHARACTER_VALUE(nodeP));
  PROTECT(rv = NEW_NUMERIC(1));
  REAL(rv)[0] = tr_distance_to_root(node);
  UNPROTECT(1);
  return rv;
}

SEXP rph_tree_label_branches(SEXP treeP, SEXP nodep, SEXP labelP) {
  int i, j, numtree = LENGTH(treeP);
  TreeNode *tr;
  SEXP result;
  char *label;

  label = copy_charstr(CHARACTER_VALUE(labelP));
  PROTECT(result = NEW_CHARACTER(numtree));
  for (i=0; i < numtree; i++) {
    tr = rph_tree_new(STRING_ELT(treeP, i));
    for (j = 0; j < LENGTH(nodep); j++) 
      tr_label_node(tr, CHARACTER_VALUE(STRING_ELT(nodep, j)), label);
    SET_STRING_ELT(result, i, mkChar(tr_to_string(tr, 1)));
  }
  UNPROTECT(1);
  return result;
}
  
SEXP rph_tree_label_subtree(SEXP treeP, SEXP nodeP,
			    SEXP includeLeadingBranchP,
			    SEXP labelP) {
  int i, numtree = LENGTH(treeP), include_leading_branch=0;
  TreeNode *tr;
  SEXP result;
  char *node, *label;
  label = copy_charstr(CHARACTER_VALUE(labelP));
  node = copy_charstr(CHARACTER_VALUE(nodeP));
  include_leading_branch = LOGICAL_VALUE(includeLeadingBranchP);
  PROTECT(result = NEW_CHARACTER(numtree));
  for (i=0; i < numtree; i++) {
    tr = rph_tree_new(STRING_ELT(treeP, i));
    tr_label_subtree(tr, node, include_leading_branch, label);
    SET_STRING_ELT(result, i, mkChar(tr_to_string(tr, 1)));
  }
  UNPROTECT(1);
  return result;
}


SEXP rph_tree_summary_nodenames(SEXP treeP) {
  TreeNode *tr = rph_tree_new(treeP), *node;
  int i;
  List *nodes = tr_preorder(tr);
  SEXP result = PROTECT(NEW_CHARACTER(lst_size(nodes)));
  for (i=0; i < lst_size(nodes); i++) {
    node = (TreeNode*)lst_get_ptr(nodes, i);
    if (strlen(node->name)==0 || strcmp(node->name, ";")==0)
      SET_STRING_ELT(result, i, NA_STRING);
    else SET_STRING_ELT(result, i, mkChar(node->name));
  }
  UNPROTECT(1);
  return result;
}


SEXP rph_tree_summary_len(SEXP treeP) {
  TreeNode *tr = rph_tree_new(treeP), *node;
  int i;
  List *nodes = tr_preorder(tr);
  SEXP result = PROTECT(NEW_NUMERIC(lst_size(nodes)));
  double *d = NUMERIC_POINTER(result);
  for (i=0; i < lst_size(nodes); i++) {
    node = (TreeNode*)lst_get_ptr(nodes, i);
    if (node->parent == NULL) d[i] = -1;
    else d[i] = node->dparent;
  }
  UNPROTECT(1);
  return result;
}


SEXP rph_tree_summary_depth(SEXP treeP) {
  TreeNode *tr = rph_tree_new(treeP), *node;
  int i;
  List *nodes = tr_preorder(tr);
  SEXP result = PROTECT(NEW_NUMERIC(lst_size(nodes)));
  double *d = NUMERIC_POINTER(result);
  for (i=0; i < lst_size(nodes); i++) {
    node = (TreeNode*)lst_get_ptr(nodes, i);
    d[i] = tr_distance_to_root(node);
  }
  UNPROTECT(1);
  return result;
}


SEXP rph_tree_summary_parent(SEXP treeP) {
  TreeNode *tr = rph_tree_new(treeP), *node;
  int i, *parent, nnode, *idmap;
  List *nodes = tr_preorder(tr);
  SEXP result;

  nnode = lst_size(nodes);
  result = PROTECT(NEW_INTEGER(nnode));
  parent = INTEGER_POINTER(result);
  idmap = smalloc((nnode+1)*sizeof(int));
  for (i=0; i < nnode; i++) {
    node = (TreeNode*)lst_get_ptr(nodes, i);
    if (node->id > nnode || node->id < 0)
      die("invalid id (%i) in tree node\n", node->id);
    idmap[(int)node->id] = i;
  }
  for (i=0; i < nnode; i++) {
    node = (TreeNode*)lst_get_ptr(nodes, i);
    if (node->parent == NULL) 
      parent[idmap[node->id]] = -1;
    else parent[idmap[node->id]] = idmap[node->parent->id] + 1;
  }
  UNPROTECT(1);
  return result;
}


SEXP rph_tree_summary_lchild(SEXP treeP) {
  TreeNode *tr = rph_tree_new(treeP), *node;
  int i, *lchild, nnode, *idmap;
  List *nodes = tr_preorder(tr);
  SEXP result;

  nnode = lst_size(nodes);
  result = PROTECT(NEW_INTEGER(nnode));
  lchild = INTEGER_POINTER(result);
  idmap = smalloc((nnode+1)*sizeof(int));
  for (i=0; i < nnode; i++) {
    node = (TreeNode*)lst_get_ptr(nodes, i);
    if (node->id > nnode || node->id < 0)
      die("invalid id (%i) in tree node\n", node->id);
    idmap[(int)node->id] = i;
  }
  for (i=0; i < nnode; i++) {
    node = (TreeNode*)lst_get_ptr(nodes, i);
    if (node->lchild == NULL) 
      lchild[idmap[node->id]] = -1;
    else lchild[idmap[node->id]] = idmap[node->lchild->id] + 1;
  }
  UNPROTECT(1);
  return result;
}


SEXP rph_tree_summary_rchild(SEXP treeP) {
  TreeNode *tr = rph_tree_new(treeP), *node;
  int i, *rchild, nnode, *idmap;
  List *nodes = tr_preorder(tr);
  SEXP result;

  nnode = lst_size(nodes);
  result = PROTECT(NEW_INTEGER(nnode));
  rchild = INTEGER_POINTER(result);
  idmap = smalloc((nnode+1)*sizeof(int));
  for (i=0; i < nnode; i++) {
    node = (TreeNode*)lst_get_ptr(nodes, i);
    if (node->id > nnode || node->id < 0)
      die("invalid id (%i) in tree node\n", node->id);
    idmap[(int)node->id] = i;
  }
  for (i=0; i < nnode; i++) {
    node = (TreeNode*)lst_get_ptr(nodes, i);
    if (node->rchild == NULL) 
      rchild[idmap[node->id]] = -1;
    else rchild[idmap[node->id]] = idmap[node->rchild->id] + 1;
  }
  UNPROTECT(1);
  return result;
}


SEXP rph_tree_summary_label(SEXP treeP) {
  TreeNode *tr = rph_tree_new(treeP), *node;
  List *nodes = tr_preorder(tr);
  int i, nnodes = lst_size(nodes);
  SEXP result;
  for (i=0; i < nnodes; i++) {
    node = (TreeNode*)lst_get_ptr(nodes, i);
    if (node->label != NULL) break;
  }
  if (i == nnodes) return R_NilValue;
  
  PROTECT(result = NEW_CHARACTER(nnodes));
  for (i=0; i < nnodes; i++) {
    node = (TreeNode*)lst_get_ptr(nodes, i);
    if (node->label == NULL) SET_STRING_ELT(result, i, NA_STRING);
    else SET_STRING_ELT(result, i, mkChar(node->label));
  }
  UNPROTECT(1);
  return result;
}
