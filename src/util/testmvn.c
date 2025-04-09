#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <phast/misc.h>
#include <phast/nj.h>
#include <phast/mvn.h>
#include <phast/multi_mvn.h>

#define MVN_N 3
#define MVN_D 3

int main(int argc, char *argv[]) {
  int i, j, d, d2, idx1, idx2;
  double val;
  Vector *x = vec_new(MVN_N * MVN_D), *x_std = vec_new(MVN_N * MVN_D);
  MVN *mvn;
  multi_MVN *mmvn;

  /* set up standard MVN */
  mvn = mvn_new(MVN_N * MVN_D, NULL, NULL);
  for (i = 0; i < MVN_N; i++) {
    for (d = 0; d < MVN_D; d++) {
      idx1 = i*MVN_D + d;
      vec_set(mvn->mu, idx1, idx1);
      for (j = 0; j < MVN_N; j++) {
        for (d2 = 0; d2 < MVN_D; d2++) {
          idx2 = j*MVN_D + d2;
          if (idx1 == idx2)
            val = 0.5;
          else if (d == d2)
            val = 0.2;
          else
            val = 0;
          mat_set(mvn->sigma, idx1, idx2, val);
        }
      }
    }
  }  
  mvn_update_type(mvn);
  mvn_preprocess(mvn, TRUE);  /* CHECK: when does this happen with a multi_mvn? */

  /* set up equivalent multi MVN */
  mmvn = mmvn_new(MVN_N, MVN_D, MVN_GEN);  
  mmvn_restore_mu(mmvn, mvn->mu);

  for (i = 0; i < MVN_N; i++) {
    for (j = 0; j < MVN_N; j++) {
      if (i == j)
        val = 0.5;
      else
        val = 0.2;
      mat_set(mmvn->mvn->sigma, i, j, val);
    }
  }  

  mvn_update_type(mmvn->mvn);  /* FIXME: is this happening in main prog? */
  mvn_preprocess(mmvn->mvn, TRUE); 
  
  /* try diagonal and other parameterization */

  for (i = 0; i < 10; i++) {
    double ldens1, ldet1, ldens2, ldet2, trace1, trace2, mu2_1, mu2_2;
    printf("Draw #%d:\n", i+1);
    mvn_sample(mvn, x);
    vec_print(x, stdout);

    mvn_rederive_std(mvn, x, x_std);
    vec_print(x_std, stdout);

    mmvn_rederive_std(mmvn, x, x_std);
    vec_print(x_std, stdout);
    /* note that, suprisingly, these two values of x_std are not
       guaranteed to be identical because of differences in sign in
       eigenvectors. However, they must correspond to the same
       density */
    
    ldens1 = mvn_log_dens(mvn, x);
    ldet1 = mvn_log_det(mvn);
    trace1 = mvn_trace(mvn);
    mu2_1 = mvn_mu2(mvn);
    
    ldens2 = mmvn_log_dens(mmvn, x);
    ldet2 = mmvn_log_det(mmvn);
    trace2 = mmvn_trace(mmvn);
    mu2_2 = mmvn_mu2(mmvn);
    
    printf("STATS: %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", ldens1, ldens2, ldet1, ldet2, trace1, trace2, mu2_1, mu2_2);
  }

  /* now do the same thing but sample from the mmvn */
  for (i = 0; i < 10; i++) {
    double ldens1, ldet1, ldens2, ldet2, trace1, trace2, mu2_1, mu2_2;
    printf("Draw #%d:\n", i+1);
    mmvn_sample(mmvn, x);
    vec_print(x, stdout);

    mvn_rederive_std(mvn, x, x_std);
    vec_print(x_std, stdout);

    mmvn_rederive_std(mmvn, x, x_std);
    vec_print(x_std, stdout);
    
    ldens1 = mvn_log_dens(mvn, x);
    ldet1 = mvn_log_det(mvn);
    trace1 = mvn_trace(mvn);
    mu2_1 = mvn_mu2(mvn);
    
    ldens2 = mmvn_log_dens(mmvn, x);
    ldet2 = mmvn_log_det(mmvn);
    trace2 = mmvn_trace(mmvn);
    mu2_2 = mmvn_mu2(mmvn);
    
    printf("STATS: %f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", ldens1, ldens2, ldet1, ldet2, trace1, trace2, mu2_1, mu2_2);
  }

  printf("MVN:\n");
  mvn_print(mvn, stdout);

  printf("MMVN:\n");
  mmvn_print(mmvn, stdout, FALSE, TRUE);

  return(0);
}
