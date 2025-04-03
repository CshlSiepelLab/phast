#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <phast/misc.h>
#include <phast/nj.h>
#include <phast/mvn.h>

#define MVNDIM 2

int main(int argc, char *argv[]) {
  int i;
  Vector *x = vec_new(MVNDIM);
  MVN *mvn = mvn_new(MVNDIM, NULL, NULL, MVN_STD);
  /* check: what happens other types; is this sensible? */

  vec_set(mvn->mu, 0, 10);
  vec_set(mvn->mu, 1, -10);
  
  mat_set(mvn->sigma, 0, 0, 0.8);
  mat_set(mvn->sigma, 0, 1, -.2);
  mat_set(mvn->sigma, 1, 0, -.2);
  mat_set(mvn->sigma, 1, 1, 0.8);
  mvn->type = MVN_GEN;
  mvn_update_cholesky(mvn);
  
  for (i = 0; i < 1000; i++) {
    mvn_sample(mvn, x);

    vec_print(x, stdout);
  }
  return(0);
}
