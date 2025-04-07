#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <phast/misc.h>
#include <phast/nj.h>
#include <phast/mvn.h>

#define MVNDIM 5

int main(int argc, char *argv[]) {
  int i, j;
  Vector *x = vec_new(MVNDIM), *x_std = vec_new(MVNDIM);
  MVN *mvn = mvn_new(MVNDIM, NULL, NULL);
  /* check: what happens other types; is this sensible? */

  for (i = 0; i < MVNDIM; i++) {
    vec_set(mvn->mu, i, i);
    for (j = 0; j < MVNDIM; j++) {
      if (i == j)
        mat_set(mvn->sigma, i, j, 0.5);
      else
        mat_set(mvn->sigma, i, j, 0.1);          
    }
  }
  
  mvn_update_type(mvn);
  mvn_preprocess(mvn, TRUE);
  
  for (i = 0; i < 10; i++) {
    double ldens1, ldet1, ldens2, ldet2;
    mvn_sample(mvn, x);
    vec_print(x, stdout);

    mvn_rederive_std(mvn, x, x_std);
    vec_print(x_std, stdout);
    
    /* mvn_preprocess(mvn, TRUE); */
    ldens1 = mvn_log_dens(mvn, x);
    ldet1 = mvn_log_det(mvn);
    /*    mvn_print(mvn, stdout); */

    /* mvn_preprocess(mvn, TRUE); */
    ldens2 = mvn_log_dens(mvn, x);
    ldet2 = mvn_log_det(mvn);
    /*    mvn_print(mvn, stdout); */

    printf("%f\t%f\t%f\t%f\n", ldens1, ldens2, ldet1, ldet2);
  }
  return(0);
}
