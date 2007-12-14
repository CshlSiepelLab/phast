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
              filename ends with '.nh', in which case it will be\n\
              expected to be a tree file in Newick format.\n\
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
    --tree-only, -t\n\
        Output tree only in Newick format rather than complete tree model.\n\
\n\
    --dissect, -d\n\
        In place of ordinary output, print a description of the id,\n\
        label (name), parent, children, and distance to parent for\n\
        each node of the tree.  Sometimes useful for debugging.  Can be\n\
        used with other options.\n\
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
    --help, -h\n\
        Print this help message.\n\n", prog, prog);
  exit(0);
}

int main(int argc, char *argv[]) {
  /* variables for options, with defaults */
  TreeNode *tree = NULL, *merge_tree = NULL, *extrapolate_tree = NULL;
  Hashtable *rename_hash = NULL;
  double scale_factor = 1;
  List *prune_names = NULL;
  int prune_all_but = FALSE, tree_only = FALSE, dissect = FALSE,
    name_ancestors = FALSE, with_branch = FALSE;
  TreeModel *mod = NULL, *merge_mod = NULL;
  char *reroot_name = NULL, *subtree_name = FALSE;
  
  /* other variables */
  String *suffix;
  char c;
  int i, opt_idx;
  TreeNode *n;

  struct option long_opts[] = {
    {"scale", 1, 0, 's'},
    {"extrapolate", 1, 0, 'e'},
    {"prune", 1, 0, 'p'},
    {"prune-all-but", 1, 0, 'P'},
    {"merge", 1, 0, 'm'},
    {"rename", 1, 0, 'r'},
    {"tree-only", 0, 0, 't'},
    {"dissect", 0, 0, 'd'},
    {"name-ancestors", 0, 0, 'a'},
    {"reroot", 1, 0, 'R'},
    {"with-branch", 1, 0, 'B'},
    {"subtree", 1, 0, 'S'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "s:p:P:m:r:R:B:S:adth", 
                          long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 's':
      scale_factor = get_arg_dbl_bounds(optarg, 0, INFTY);
      break;
    case 'e':
      if (!strcmp(optarg, "default")) {
        optarg = smalloc(1000 * sizeof(char));
        sprintf(optarg, "%s/data/exoniphy/mammals/cftr25_hybrid.nh", 
                PHAST_HOME);
      }
      extrapolate_tree = tr_new_from_file(fopen_fname(optarg, "r"));
      break;
    case 'p':
      prune_names = get_arg_list(optarg);
      break;
    case 'P':
      prune_names = get_arg_list(optarg);
      prune_all_but = TRUE;
      break;
    case 'm':
      suffix = str_new_charstr(optarg);
      str_suffix(suffix, '.');
      if (str_equals_charstr(suffix, "nh"))
        merge_tree = tr_new_from_file(fopen_fname(optarg, "r"));
      else {
        merge_mod = tm_new_from_file(fopen_fname(optarg, "r"));
        merge_tree = merge_mod->tree;
      }
      break;
    case 'r':
      rename_hash = make_name_hash(optarg);
      break;
    case 't':
      tree_only = TRUE;
      break;
    case 'd':
      dissect = TRUE;
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
    
  suffix = str_new_charstr(argv[optind]);
  str_suffix(suffix, '.');
  if (str_equals_charstr(suffix, "nh")) {
    tree = tr_new_from_file(fopen_fname(argv[optind], "r"));
    tree_only = TRUE;           /* can't output tree model in this case */
  }
  else {
    mod = tm_new_from_file(fopen_fname(argv[optind], "r"));
    tree = mod->tree;
  }

  if (prune_names != NULL) {
    tr_prune(&tree, prune_names, prune_all_but);
    if (mod != NULL) mod->tree = tree; /* root may have changed */
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
      tr_scale_subtree(tree, n, scale_factor);
      if (with_branch) n->dparent *= scale_factor;
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

  if (dissect) 
    tr_print_nodes(stdout, tree);
  else if (tree_only)
    tr_print(stdout, tree, TRUE);
  else
    tm_print(stdout, mod);
  
  return 0;
}
