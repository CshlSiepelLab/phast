#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <phast/misc.h>
#include <phast/sparse_matrix.h>

int main(int argc, char *argv[]) {
  int i, j;
  double u;
  SparseMatrix *sm = spmat_new(1000,1000,10);
  for (i = 0; i < 1000; i++) {
    for (j = 0; j < 1000; j++) {
      u = unif_rand();
      if (u < 0.01) {
        spmat_set(sm, i, j, u);
        printf("setting %d %d to %f\n", i, j, u);
      }
    }
  }

  for (i = 0; i < 1000; i++) {
    for (j = 0; j < 1000; j++) {
      u = spmat_get(sm, i, j);
      if (u > 0)
        printf("value at %d %d is %f\n", i, j, u);
    }
  }
}
