#include <tree_model.h>

void print_usage() {
  printf("\n\
PROGRAM: merge_tree_models\n\
\n\
USAGE: merge_tree_models <branchlens.mod> <subst.mod> [<norm.mod>]\n\
\n\
DESCRIPTION: Creates a new tree model with the branch lengths of one\n\
tree model and the substitution parameters (equilibrium frequencies,\n\
rate matrix, rate variation properties) of another.  Useful for\n\
approximating a richly parameterized tree model for a large set of\n\
species when there are insufficient sites to estimate all parameters.\n\
If a third tree model is specified (<norm.mod>), then the branches of\n\
<branchlens.mod> will be scaled by a factor equal to the total branch\n\
length of <subst.mod> divided by the total branch length of\n\
<norm.mod>.\n\n");

}

int main(int argc, char *argv[]) {
  TreeModel *bl, *subst, *norm;

  if (argc < 3 || argc > 4) { print_usage(); exit(1); }

  bl = tm_new_from_file(fopen_fname(argv[1], "r"));
  subst = tm_new_from_file(fopen_fname(argv[2], "r"));
  norm = argc > 3 ? 
    tm_new_from_file(fopen_fname(argv[3], "r")) : NULL;

  if (subst->tree == NULL) {    /* weight matrix */
    tm_print(stdout, subst);
    exit (0);
  }

  if (norm != NULL && norm->tree != NULL)
    tm_scale(bl, tr_total_len(subst->tree)/tr_total_len(norm->tree), 0);

  subst->tree = bl->tree;  
  tm_print(stdout, subst);

  return 0;
}
