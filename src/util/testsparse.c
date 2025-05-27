#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <phast/misc.h>
#include <phast/matrix.h>
#include <phast/sparse_matrix.h>

int main(int argc, char *argv[]) {
  int i, j;
  double u;
  SparseMatrix *sm = spmat_new(100000,100000,100);
  //  Matrix *sm = mat_new(10000,10000);
  for (i = 0; i < 100000; i++) {
    for (j = 0; j < 100000; j++) {
      u = unif_rand();
      if (u < 0.0001) {
        spmat_set_sorted(sm, i, j, u);
        //mat_set(sm, i, j, u);
        printf("setting %d %d to %f\n", i, j, u);
      }
    }
  }

  for (i = 0; i < 100000; i++) {
    for (j = 0; j < 100000; j++) {
      u = spmat_get(sm, i, j);
      //u = mat_get(sm, i, j);
      if (u > 0)
        printf("value at %d %d is %f\n", i, j, u);
    }
  }
}
