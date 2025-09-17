#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <phast/misc.h>
#include <phast/nj.h>
#include <phast/radial_flow.h>
#include <phast/mvn.h>

#define MVN_N 10
#define MVN_D 3
#define EPS 1e-6

/* FIXME: implement scale invariance first */

/* dummy elbo: log p(y) + logdet */
double elbo(Vector *y, double logdet, MVN *mvn) {
  return mvn_log_dens(mvn, y) + logdet;
}

int main(int argc, char *argv[]) {
  int nd = MVN_N * MVN_D;
  Vector *x = vec_new(nd);
  Vector *y = vec_new(nd);
  Vector *gy = vec_new(nd);
  Vector *gx = vec_new(nd);
  Vector *gx_num = vec_new(nd);
  
  /* draw x from standard normal */
  mvn_sample_std(x);
  printf("x: ");
  vec_print(x, stdout);
  
  /* convert x to y using radial flow */
  RadialFlow *rf = rf_new(MVN_N, MVN_D);
  rf->a = 18;
  rf->b = 7;
  rf_update(rf);
  rf_rescale(rf, POINTSPAN_EUC/sqrt(2));
  rf_forward(rf, y, x);
  printf("y: ");
  vec_print(y, stdout);

  /* compute dL/dy = -y (simple in this case) */
  vec_copy(gy, y);
  vec_scale(gy, -1.0);

  /* compute dL/dx by backprop */
  rf_backprop(rf, x, gx, gy);
  printf("gx: ");
  vec_print(gx, stdout);
  
  /* now approximate dL/dx numerically by brute force */
  MVN *mvn = mvn_new(MVN_N * MVN_D, NULL, NULL);   /* standard MVN */
  Vector *x_tweak = vec_create_copy(x);
  Vector *y_tweak = vec_new(nd);
  for (int i = 0; i < nd; i++) {
    double xi = vec_get(x, i);

    vec_set(x_tweak, i, xi + EPS);
    double logdet_tweak = rf_forward(rf, y_tweak, x_tweak);
    double elbo_tweak = elbo(y_tweak, logdet_tweak, mvn);

    vec_set(x_tweak, i, xi - EPS);
    double logdet_tweak2 = rf_forward(rf, y_tweak, x_tweak);
    double elbo_tweak2 = elbo(y_tweak, logdet_tweak2, mvn);


    double gx_num_i = (elbo_tweak - elbo_tweak2) / (2*EPS);
    vec_set(gx_num, i, gx_num_i);
    
    vec_set(x_tweak, i, xi);
  }
  printf("gx_num: ");
  vec_print(gx_num, stdout);
  
  /* now check a and b grads */
  double a0 = rf->a;
  rf->a = a0 + EPS;
  rf_update(rf);
  double logdet_p = rf_forward(rf, y_tweak, x);
  double Lp = elbo(y_tweak, logdet_p, mvn);
  rf->a = a0 - EPS;
  rf_update(rf);
  double logdet_m = rf_forward(rf, y_tweak, x);
  double Lm = elbo(y_tweak, logdet_m, mvn);
  double a_grad_num = (Lp - Lm) / (2.0 * EPS);
  rf->a = a0; 

  double b0 = rf->b;
  rf->b = b0 + EPS;
  rf_update(rf);
  logdet_p = rf_forward(rf, y_tweak, x);
  Lp = elbo(y_tweak, logdet_p, mvn);
  rf->b = b0 - EPS;
  rf_update(rf);
  logdet_m = rf_forward(rf, y_tweak, x);
  Lm = elbo(y_tweak, logdet_m, mvn);
  double b_grad_num = (Lp - Lm) / (2.0 * EPS);

  printf("grad_a: %f\ngrad_a_num: %f\n", rf->a_grad, a_grad_num);
  printf("grad_b: %f\ngrad_b_num: %f\n", rf->b_grad, b_grad_num);
  
  /* FIXME: check center also? */
  /* FIXME: loop through various parameter values? */
}
