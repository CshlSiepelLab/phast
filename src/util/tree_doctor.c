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
        node having zero children or two children).\n\
\n\
    --prune-all-but, -P <list>\n\
        Like --prune, but remove all leaves *except* the ones specified.\n\
\n\
    --merge, -m <file2.mod>|<file2.nh>\n\
        Merge with another tree model or tree.  The primary model\n\
        (<file.mod>) must have a subset of the species (leaves) in the\n\
        secondary model (<file2.mod>), and the primary tree must be a\n\
        subtree of the secondary tree (in terms of topology only).\n\
        If a full tree model is given, only the tree will be considered. \n\
        The merged tree model will have the rate matrix, equilibrium\n\
        frequencies, and rate distribution of the primary model, but a\n\
        merged tree that includes all species from both models.  The\n\
        trees will be merged by first scaling the secondary tree such\n\
        that the subtree corresponding to the primary tree has equal\n\
        overall length to the primary tree, then combining the primary\n\
        tree with the non-overlapping portion of the secondary tree.\n\
        The names of matching species (leaves) must be exactly equal.\n\
        The basic idea is to extrapolate from a small species set to a\n\
        larger one, using the branch-length proportions of a given\n\
        tree.\n\
\n\
    --rename, -r <mapping>\n\
        Rename leaves according to the given mapping.  The format of\n\
        <mapping> must be: \"oldname1 -> newname1 ; oldname2 ->\n\
        newname2 ; ...\".  This option is applied *after* all other\n\
        options (i.e., old names will be used for --prune, --merge,\n\
        etc.)\n\
\n\
    --scale, -s <factor>\n\
        Scale all branches by the specified factor.  This option is\n\
        applied *after* --merge.\n\
\n\
    --tree-only, -t\n\
        Output tree only in Newick format rather than complete tree model.\n\
\n\
    --help, -h\n\
        Print this help message.\n\n", prog, prog);
  exit(0);
}

/* parse string defining mapping from old names to new and store as
   hash */
Hashtable *make_name_hash(char *mapstr) {
  Hashtable *retval = hsh_new(20);
  Regex *map_re = str_re_new("^[[:space:]]*([A-Za-z0-9_]+)[[:space:]]*->[[:space:]]*([A-Za-z0-9_]+)[[:space:]]*");
  List *mappings = lst_new_ptr(20), *names = lst_new_ptr(3);
  String *s = str_new_charstr(mapstr);
  int i;

  str_split(s, ";", mappings);
  for (i = 0; i < lst_size(mappings); i++) {
    String *oldname, *newname;
    if (str_re_match(lst_get_ptr(mappings, i), map_re, names, 2) < 0)
      die("ERROR: cannot parse mapping ('%s')\n", lst_get_ptr(mappings, i));
    oldname = lst_get_ptr(names, 1);
    newname = lst_get_ptr(names, 2);
    hsh_put(retval, oldname->chars, strdup(newname->chars));
    lst_free_strings(names);
  }
  lst_free_strings(mappings);
  lst_free(mappings);
  lst_free(names);

  return retval;
}

int main(int argc, char *argv[]) {
  /* variables for options, with defaults */
  TreeNode *tree = NULL, *merge_tree = NULL;
  Hashtable *rename_hash = NULL;
  double scale_factor = 1;
  List *prune_names = NULL;
  int prune_all_but = FALSE, tree_only = FALSE;
  TreeModel *mod = NULL;

  /* other variables */
  String *suffix;
  char c;
  int i, opt_idx;

  struct option long_opts[] = {
    {"scale", 1, 0, 's'},
    {"prune", 1, 0, 'p'},
    {"prune-all-but", 1, 0, 'P'},
    {"merge", 1, 0, 'm'},
    {"rename", 1, 0, 'r'},
    {"tree-only", 0, 0, 't'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "s:p:P:m:r:th", 
                          long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 's':
      scale_factor = get_arg_dbl_bounds(optarg, 0, INFTY);
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
        TreeModel *tmp = tm_new_from_file(fopen_fname(optarg, "r"));
        merge_tree = tmp->tree;
        tm_free(tmp);
      }
      break;
    case 'r':
      rename_hash = make_name_hash(optarg);
      break;
    case 't':
      tree_only = TRUE;
      break;
    case 'h':
      usage(argv[0]);
    case '?':
      die("Bad argument.  Try '%s -h'.\n", argv[0]);
    }
  }

  if (optind != argc - 1) 
    die("Input filename required.  Try '%s -h'.\n", argv[0]);
    
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

  if (scale_factor != 1)
    tr_scale(tree, scale_factor);

  if (rename_hash != NULL) {
    char *newname;
    for (i = 0; i < tree->nnodes; i++) {
      TreeNode *n = lst_get_ptr(tree->nodes, i);
      if (n->name != NULL && n->name[0] != '\0' && 
          (newname = hsh_get(rename_hash, n->name)) != (char*)-1) {
        strcpy(n->name, newname);
      }
    }
  }

  if (tree_only)
    tr_print(stdout, tree, TRUE);
  else
    tm_print(stdout, mod);
  
  return 0;
}
