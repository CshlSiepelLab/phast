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
#include <hashtable.h>

void usage(char *prog) {
  printf("\n\
PROGRAM:      %s\n\
\n\
DESCRIPTION:  Scale, prune, merge, and otherwise tweak phylogenetic trees.\n\
              Expects input to be a tree model (.mod) file unless\n\
              filename ends with '.nh' or -n option is used, in which case \n\
              it will be expected to be a tree file in Newick format.\n\
\n\
USAGE:        %s [OPTIONS] <file.mod>|<file.nh>\n\
\n\
OPTIONS:\n\
    --prune, -p <list>\n\
        Remove all leaves whose names are included in the given list\n\
        (comma-separated), then remove nodes and combine branches\n\
        to restore as a complete binary tree (i.e., with each\n\
        node having zero children or two children).  This option is\n\
        applied *before* all other options.\n\
\n\
    --prune-all-but, -P <list>\n\
        Like --prune, but remove all leaves *except* the ones specified.\n\
\n\
    --get-subtree, -g <node_name>\n\
        Like --prune, but remove all leaves who are not descendants of \n\
        node.  (Note: implies --name-ancestors if given node not \n\
        explicitly named in input tree)\n\
\n\
    --rename, -r <mapping>\n\
        Rename leaves according to the given mapping.  The format of\n\
        <mapping> must be: \"oldname1 -> newname1 ; oldname2 ->\n\
        newname2 ; ...\".  This option is applied *after* all other\n\
        options (i.e., old names will be used for --prune, --merge,\n\
        etc.)\n\
\n\
    --scale, -s <factor>\n\
        Scale all branches by the specified factor.\n\
\n\
    --name-ancestors, -a\n\
        Ensure names are assigned to all ancestral nodes.  If a node\n\
        is unnamed, create a name by concatenating the names of a leaf\n\
        from its left subtree and a leaf from its right subtree.\n\
\n\
   --label-subtree, -L <node[+]:label>\n\
        Add a label to the subtree of the named node.  If the node name\n\
        is followed by a \"+\" sign, then the branch leading to that node\n\
        is included in the subtree.  This may be used multiple times to add\n\
        more than one label, though a single branch may have only one\n\
        label.  --label-subtree and --label-branches options are parsed in\n\
        the order given, so that later uses may override earlier ones.\n\
        Labels are applied *after* all pruning, re-rooting, and re-naming\n\
        options are applied.\n\
\n\
    --label-branches, -l <branch1,branch2,...:label>\n\
        Add a label to the branches listed.  Branches are named by the name\n\
        of the node which descends from that branch.  See --label-subtree\n\
        above for more information.\n\
\n\
    --tree-only, -t\n\
        Output tree only in Newick format rather than complete tree model.\n\
\n\
    --no-branchlen, -N\n\
        (Implies --tree-only).  Output only topology in Newick format.\n\
\n\
    --dissect, -d\n\
        In place of ordinary output, print a description of the id,\n\
        name, parent, children, and distance to parent for each node\n\
        of the tree.  Sometimes useful for debugging.  Can be used with\n\
        other options.\n\
\n\
    --branchlen, -b\n\
        In place of ordinary output, print the total branch length of\n\
        the tree that would have been printed.\n\
\n\
    --depth, -D <node_name>\n\
        In place of ordinary output, report distance from named node to\n\
        root\n\
\n\
    --reroot, -R <node_name>\n\
        Reroot tree at internal node with specified name.\n\
\n\
    --subtree, -S <node_name>\n\
        (for use with --scale)  Alter only the branches in the subtree \n\
        beneath the specified node.\n\
\n\
    --with-branch, -B <node_name>\n\
        (For use with --reroot or --subtree) include branch above specified\n\
        node with subtree beneath it.\n\
\n\
    --merge, -m <file2.mod> | <file2.nh>\n\
        Merge with another tree model or tree.  The primary model\n\
        (<file.mod>) must have a subset of the species (leaves) in the\n\
        secondary model (<file2.mod>), and the primary tree must be a\n\
        proper subtree of the secondary tree (i.e., the subtree of the\n\
        secondary tree beneath the LCA of the species in the primary\n\
        tree must equal the primary tree in terms of topology and\n\
        species names).  If a full tree model is given for the\n\
        secondary tree, only the tree will be considered.  The merged\n\
        tree model will have the rate matrix, equilibrium frequencies,\n\
        and rate distribution of the primary model, but a merged tree\n\
        that includes all species from both models.  The trees will be\n\
        merged by first scaling the secondary tree such that the\n\
        subtree corresponding to the primary tree has equal overall\n\
        length to the primary tree, then combining the primary tree\n\
        with the non-overlapping portion of the secondary tree.  The\n\
        names of matching species (leaves) must be exactly equal.\n\
\n\
    --extrapolate, -e <phylog.nh> | default\n\
        Extrapolate to a larger set of species based on the given\n\
        phylogeny (Newick-format).  The primary tree must be a subtree\n\
        of the phylogeny given in <phylog.nh>, but it need not be\n\
        a \"proper\" subtree (see --merge).  A copy will be created\n\
        of the larger phylogeny then scaled such that the total branch\n\
        length of the subtree corresponding to the primary tree equals\n\
        the total branch length of the primary tree; this new version\n\
        will then be used in place of the primary tree.  If the string\n\
        \"default\" is given instead of a filename, then a phylogeny\n\
        for 25 vertebrate species, estimated from sequence data for\n\
        Target 1 (CFTR) of the NISC Comparative Sequencing Program\n\
        (Thomas et al., Nature 424:788-793, 2003), will be assumed.\n\
        This option is similar to merge but differs in that the branch\n\
        length proportions of the output tree come completely from the\n\
        larger tree and the smaller tree doesn't have to be a proper\n\
        subset of the larger tree.\n\
\n\
    --newick,-n\n\
        The input file is in Newick format (necessary if file name does\n\
        not end in .nh)\n\
\n\
    --help, -h\n\
        Print this help message.\n\n", prog, prog);
  exit(0);
}

int main(int argc, char *argv[]) {
  /* variables for options, with defaults */
  TreeNode *tree = NULL, *merge_tree = NULL, *extrapolate_tree = NULL;
  Hashtable *rename_hash = NULL;
  double scale_factor = 1;
  List *prune_names = NULL, *label = NULL, *labelType = NULL;
  int prune_all_but = FALSE, tree_only = FALSE, dissect = FALSE,
    name_ancestors = FALSE, with_branch = FALSE, print_branchlen=FALSE,
    inNewick=FALSE, no_branchlen = FALSE, print_distance_to_root = FALSE;
  TreeModel *mod = NULL, *merge_mod = NULL;
  char *reroot_name = NULL, *subtree_name =NULL, *get_subtree_name = NULL,
    *node_distance_name = NULL;
  
  /* other variables */
  String *suffix,  *optstr;
  char c;
  int i, opt_idx;
  TreeNode *n;

  struct option long_opts[] = {
    {"scale", 1, 0, 's'},
    {"extrapolate", 1, 0, 'e'},
    {"prune", 1, 0, 'p'},
    {"prune-all-but", 1, 0, 'P'},
    {"get-subtree", 1, 0, 'g'},
    {"merge", 1, 0, 'm'},
    {"rename", 1, 0, 'r'},
    {"tree-only", 0, 0, 't'},
    {"no-branchlen", 0, 0, 'N'},
    {"dissect", 0, 0, 'd'},
    {"name-ancestors", 0, 0, 'a'},
    {"reroot", 1, 0, 'R'},
    {"with-branch", 1, 0, 'B'},
    {"subtree", 1, 0, 'S'},
    {"branchlen", 0, 0, 'b'},
    {"newick", 0, 0, 'n'},
    {"label-subtree", 1, 0, 'L'},
    {"label-branches", 1, 0, 'l'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "s:p:P:g:m:r:R:B:S:D:l:L:adtNbnh", 
                          long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 's':
      scale_factor = get_arg_dbl_bounds(optarg, 0, INFTY);
      break;
    case 'e':
      if (!strcmp(optarg, "default")) {
        optarg = smalloc(1000 * sizeof(char));
        #if defined(__MINGW32__)
          sprintf(optarg, "%s\\data\\exoniphy\\mammals\\cftr25_hybrid.nh",
		  PHAST_HOME);
        #else
          sprintf(optarg, "%s/data/exoniphy/mammals/cftr25_hybrid.nh", 
                  PHAST_HOME);
        #endif
      }
      extrapolate_tree = tr_new_from_file(phast_fopen(optarg, "r"));
      break;
    case 'p':
      prune_names = get_arg_list(optarg);
      break;
    case 'P':
      prune_names = get_arg_list(optarg);
      prune_all_but = TRUE;
      break;
    case 'g':
      get_subtree_name = optarg;
      break;
    case 'm':
      suffix = str_new_charstr(optarg);
      str_suffix(suffix, '.');
      if (str_equals_charstr(suffix, "nh"))
        merge_tree = tr_new_from_file(phast_fopen(optarg, "r"));
      else {
        merge_mod = tm_new_from_file(phast_fopen(optarg, "r"), 1);
        merge_tree = merge_mod->tree;
      }
      break;
    case 'r':
      rename_hash = make_name_hash(optarg);
      break;
    case 't':
      tree_only = TRUE;
      break;
    case 'N':
      no_branchlen = TRUE;
      tree_only = TRUE;
      break;
    case 'd':
      dissect = TRUE;
      break;
    case 'b':
      print_branchlen = TRUE;
      break;
    case 'D':
      print_distance_to_root = TRUE;
      node_distance_name = optarg;
      break;
    case 'R':
      reroot_name = optarg;
      break;
    case 'B':
      with_branch = TRUE;
      break;
    case 'a':
      name_ancestors = TRUE;
      break;
    case 'S':
      subtree_name = optarg;
      break;
    case 'n':
      inNewick=TRUE;
      break;
    case 'L':  //do the same for --label--subtree and --label-branches
    case 'l':
      if (label == NULL) {
	label = lst_new_ptr(1);
	labelType = lst_new_int(1);
      }
      optstr = str_new_charstr(optarg);
      lst_push_ptr(label, optstr);
      lst_push_int(labelType, (int)c);
      break;
    case 'h':
      usage(argv[0]);
    case '?':
      die("Bad argument.  Try '%s -h'.\n", argv[0]);
    }
  }

  if (optind != argc - 1) 
    die("Input filename required.  Try '%s -h'.\n", argv[0]);

  if (merge_tree != NULL && extrapolate_tree != NULL)
    die("ERROR: Can't use --merge and --extrapolate together");

  set_seed(-1);
    
  suffix = str_new_charstr(argv[optind]);
  str_suffix(suffix, '.');
  if (inNewick || str_equals_charstr(suffix, "nh")) {
    tree = tr_new_from_file(phast_fopen(argv[optind], "r"));
    tree_only = TRUE;           /* can't output tree model in this case */
  }
  else {
    mod = tm_new_from_file(phast_fopen(argv[optind], "r"), 1);
    tree = mod->tree;
  }

  if (prune_names != NULL) {
    tr_prune(&tree, prune_names, prune_all_but, NULL);
    if (mod != NULL) mod->tree = tree; /* root may have changed */
  }

  if (get_subtree_name != NULL) {
    n = tr_get_node(tree, get_subtree_name);
    if (n == NULL) {
      tr_name_ancestors(tree);
      n = tr_get_node(tree, get_subtree_name);
      if (n == NULL) {
	die("ERROR: no node named '%s'.\n", subtree_name);
      }
    }
    tr_prune_supertree(&tree, n);
    if (mod != NULL) mod->tree = tree;
  }

  if (merge_tree != NULL) {
    tree = tr_hybrid(tree, merge_tree);
    if (mod != NULL) mod->tree = tree;
  }

  else if (extrapolate_tree != NULL) {
    tr_scale_by_subtree(extrapolate_tree, tree);
    tree = extrapolate_tree;
    if (mod != NULL) mod->tree = tree;
  }

  if (scale_factor != 1) {
    if (subtree_name == NULL)
      tr_scale(tree, scale_factor);
    else {
      n = tr_get_node(tree, subtree_name);
      if (n == NULL) die("ERROR: no node named '%s'.\n", subtree_name);
      tr_scale_subtree(tree, n, scale_factor, with_branch);
    }
  }

  if (name_ancestors)
    tr_name_ancestors(tree);

  if (rename_hash != NULL) {
    char *newname;
    for (i = 0; i < tree->nnodes; i++) {
      n = lst_get_ptr(tree->nodes, i);
      if (n->name != NULL && n->name[0] != '\0' && 
          (newname = hsh_get(rename_hash, n->name)) != (char*)-1) {
        strcpy(n->name, newname);
      }
    }
  }

  if (reroot_name != NULL) {
    n = tr_get_node(tree, reroot_name);
    if (n == NULL) die("ERROR: no node named '%s'.\n", reroot_name);
    tr_reroot(tree, n, with_branch);
    if (mod != NULL) mod->tree = with_branch ? n->parent : n;
    tree = with_branch ? n->parent : n;
  }

  if (label != NULL) {
    for (i=0; i < lst_size(label); i++) {
      String *currstr = (String*)lst_get_ptr(label, i), *arg1, *labelVal;
      List *tmplst = lst_new_ptr(10);
      String *nodename;
      int j;
      str_split(currstr, ":", tmplst);
      if (lst_size(tmplst) != 2) 
	die("ERROR: bad argument to --label-branches or --label-subtree.\n");
      arg1 = lst_get_ptr(tmplst, 0);
      labelVal = lst_get_ptr(tmplst, 1);
      lst_clear(tmplst);
      if (lst_get_int(labelType, i) == (int)'l') {
	str_split(arg1, ",", tmplst);
	for (j=0; j < lst_size(tmplst); j++) {
	  nodename = (String*)lst_get_ptr(tmplst, j);
	  tr_label_node(tree, nodename->chars, labelVal->chars);
	}
	lst_free_strings(tmplst);
      } else if (lst_get_int(labelType, i) == (int)'L') {
	int include_leading_branch = FALSE;
	TreeNode *node;
	nodename = arg1;
	node = tr_get_node(tree, nodename->chars);
	if (node == NULL && nodename->chars[nodename->length-1] == '+') {
	  nodename->chars[--nodename->length] = '\0';
	  node = tr_get_node(tree, nodename->chars);
	  include_leading_branch = TRUE;
	}
	tr_label_subtree(tree, nodename->chars, include_leading_branch, 
			 labelVal->chars);
      } else die("ERROR got label_type %c\n", lst_get_int(labelType, (char)i));
      str_free(arg1);
      str_free(labelVal);
      lst_free(tmplst);
      str_free(currstr);
    }
    lst_free(label);
    lst_free(labelType);
  }

  if (dissect) 
    tr_print_nodes(stdout, tree);
  if (print_branchlen) 
    printf("TOTAL_TREE_LEN: %f\n", tr_total_len(tree));
  if (print_distance_to_root) {
    TreeNode *node = tr_get_node(tree, node_distance_name);
    if (node == NULL) 
      die("ERROR: no node named '%s'.\n", node_distance_name);
    printf("length(root-%s): %f\n", node_distance_name, 
	   tr_distance_to_root(node));
  }

  if (dissect==0 && print_branchlen==0 && print_distance_to_root==0) {
    if (tree_only)
      tr_print(stdout, tree, no_branchlen==FALSE);
    else
      tm_print(stdout, mod);
  }
  return 0;
}
