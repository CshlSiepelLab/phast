#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <phast/misc.h>
#include <phast/nj.h>
#include <phast/mvn.h>

#define MVNDIM 3

int main(int argc, char *argv[]) {
  int i;
  double ldens, ldet;
  Vector *x = vec_new(MVNDIM);
  MVN *mvn = mvn_new(MVNDIM, NULL, NULL, MVN_STD);
  /* check: what happens other types; is this sensible? */

  for (i = 0; i < 1; i++) {
    mvn_sample_std(x);
    /* extend to other covariance matrices */
    vec_print(x, stdout);

    ldens = mvn_log_dens(mvn, x);
    ldet = mvn_log_det(mvn);

    printf("%f\t%f\n", ldens, ldet);
    vec_print(mvn->mu, stdout);
    mat_print(mvn->sigma, stdout);
  }
  return(0);
}
