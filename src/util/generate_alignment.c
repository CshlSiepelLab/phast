#include <stdlib.h>
#include <stdio.h>
#include "tree_model.h"
#include "markov_matrix.h"
#include "msa.h"

#define MSA_FNAME "alignment.ph"
#define GFF_FNAME "alignment.gff"

void print_usage() {
    fprintf(stderr, "USAGE: generate_alignment <ncolumns> <class-transition-matrix_fname>\n\t<class-1-model_fname> <class-2-model_fname> ...\n");
    fprintf(stderr, "\nNumber of class-specific models must equal number of states in transition\nmatrix.  Each model must have a tree with the same number N of taxa,\nlabeled 1...N.  Tree topologies need not be the same (but normally will be).\n");
}

int main(int argc, char* argv[]) {
  FILE* F;
  MarkovMatrix *classmat;
  TreeModel **classmods;
  MSA *msa;
  int ncols, i;
  int *labels;
  int *path_to_cat, *reverse_compl;
  CategoryMap *cm;
  GFF_Set *gff;

  if (argc < 4) {
    print_usage();
    exit(1);
  }

  ncols = atoi(argv[1]);

  /* read class transition matrix */
  F = fopen(argv[2], "r");
  classmat = mm_new_from_file(F, DISCRETE);
  fclose(F);

  if (classmat->size <= 0 || argc - 3 != classmat->size) {
    print_usage();
    exit(1);
  }

  /* read class-specific tree models */
  classmods = (TreeModel**)malloc(classmat->size * sizeof(TreeModel*));
  for (i = 0; i < classmat->size; i++) {
    F = fopen(argv[i+3], "r");
    classmods[i] = tm_new_from_file(F);
    fclose(F);
  }

  /* generate and print alignment */
  labels = (int*)malloc(ncols * sizeof(int));
  msa = tm_generate_msa(ncols, classmat, classmods, labels);
  F = fopen(MSA_FNAME, "w+");
  msa_print(F, msa, PHYLIP, 0);
  fclose(F);

  /* also print labels */
  cm = cm_read(fopen_fname("cats.cm", "r"));
  path_to_cat = smalloc(classmat->size * sizeof(int));
  reverse_compl = smalloc(classmat->size * sizeof(int));
  for (i = 0; i < classmat->size; i++) {
    path_to_cat[i] = i; 
    reverse_compl[i] = 0;
  }
  gff = cm_labeling_as_gff(cm, labels, path_to_cat, reverse_compl, msa->length, 
                           str_new_charstr("simulated"), 
                           str_new_charstr("generate_alignment"), 0, '+',
                           NULL, NULL);
  free(path_to_cat);
  free(reverse_compl);

  F = fopen(GFF_FNAME, "w+");
  gff_print_set(F, gff);
/*   for (i = 0; i < ncols; i++)  */
/*     fprintf(F, "%d\t%d\n", i+1, labels[i]); */
  fclose(F);
  gff_free_set(gff);
  
  fprintf(stderr, "Output written to %s and %s.\n", MSA_FNAME, GFF_FNAME);

  return 0;
}
