/* in SGA add a very small L2 prior to discourage extreme values, e.g., lambda * (||a||^2 + ||b||^2) where lambda is 1e-6 to 1e-4 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <float.h>
#include "phast/vector.h"
#include "phast/misc.h"
#include "phast/radial_flow.h"

RadialFlow *rf_new(int npoints, int ndim) {
  RadialFlow *rf = smalloc(sizeof(RadialFlow));
  rf->npoints = npoints;
  rf->ndim = ndim;
  rf->ctr = vec_new(ndim);
  vec_zero(rf->ctr);
  rf->ctr_grad = vec_new(ndim);
  vec_zero(rf->ctr_grad);
  rf->a = log(exp(1) - 1); /* initial values that cause near identity transformation */
  rf->b = -10;
  rf->a_grad = rf->b_grad = 0;
  rf_update(rf);
  return rf;
}

void rf_free(RadialFlow *rf) {
  vec_free(rf->ctr);
  sfree(rf);
}

/* update alpha and beta according to new value of raw parameters a and b */
void rf_update(RadialFlow *rf) {
  rf->alpha = softplus(rf->a) + RF_EPS;
  rf->beta = -rf->alpha + softplus(rf->b) + RF_EPS;
}

/* assumes points x, computes y = f(x) where x is an input set of
   points and y is transformed by the flow, assumes structured as n x
   d.  returns log determinant of Jacobian of transformation */
double rf_forward (RadialFlow *rf, Vector *y, Vector *x) {
  assert(y->size == x->size && y->size == rf->npoints * rf->ndim);
  double *diff = (double*)malloc(sizeof(double) * rf->ndim);
  double logdet = 0;
  int d;
  
  for (int i = 0; i < rf->npoints; i++) {
    int idx;
    double r2 = 0, r, h, h2, bh, A, B;
    for (d = 0; d < rf->ndim; d++) {
      idx = i*rf->ndim + d;
      diff[d] = vec_get(x, idx) - vec_get(rf->ctr, d);
      r2 += diff[d] * diff[d];
    }
    r = sqrt(r2);
    if (r < 1e-18) r = 1e-18;  /* avoid dividing by very small number */

    h = 1.0 / (rf->alpha + r);
    h2 = h * h;
    bh = rf->beta * h; 
    A  = 1.0 + bh;
    B  = A - rf->beta * h2 * r;

    for (d = 0; d < rf->ndim; d++) {
      idx = i*rf->ndim + d;
      vec_set(y, idx, vec_get(x, idx) + bh * diff[d]);
    }

    if (A <= 0.0) A = RF_EPS;
    if (B <= 0.0) B = RF_EPS;
    
    logdet += (rf->ndim - 1) * log(A) + log(B);
  }
  free(diff);
  return logdet;
}

/* param grad: first d components are derivs wrt ctr; then deriv wrt alpha and deriv wrt beta */
void rf_backprop(RadialFlow *rf, Vector *x, Vector *newgrad, Vector *origgrad,
                 Vector *paramgrad) {
  double *u = (double*)malloc(sizeof(double) * rf->ndim);

  assert(x->size == rf->npoints * rf->ndim && x->size == newgrad->size &&
         x->size == origgrad->size);
  assert(paramgrad->size == rf->ndim + 2);
  
  double alpha_grad = 0, beta_grad = 0;
  Vector *ctr_grad = vec_new(rf->ndim); vec_zero(ctr_grad);
  int d;
  
  for (int i = 0; i < rf->npoints; i++) {
    int idx;
    double r2 = 0, u_dot_g = 0, r, h, h2, h_prime, bh, A, B, dF_dr;
    
    for (d = 0; d < rf->ndim; d++) {
      idx = i*rf->ndim + d;
      u[d] = vec_get(x, idx) - vec_get(rf->ctr, idx);
      r2 += u[d] * u[d];
      u_dot_g += u[d] * vec_get(origgrad, idx);
    }

    r = sqrt(r2);  
    h = 1.0 / (rf->alpha + r);
    h2 = h*h;
    h_prime = -h2;
    bh = rf->beta * h;
    A = 1.0 + bh;
    B = A + rf->beta * h_prime * r;
    dF_dr = rf->beta * (-(rf->ndim - 1) * h2 / A + (-2.0 * h2 + 2.0 * r * h2 * h) / B);

    for (d = 0; d < rf->ndim; d++) {
      idx = i*rf->ndim + d;
      
      /* path contribution */
      double gx = A * vec_get(origgrad, idx) + rf->beta * h_prime / r * u[d] * u_dot_g;

      /* log determinant contribution */
      gx += dF_dr / r * u[d];
        
      vec_set(newgrad, idx, gx);
    }

    /* gradients with respect to parameters */
    alpha_grad += rf->beta * -h2 * u_dot_g +
      rf->beta * (-(rf->ndim - 1) * h2 / A + (-h2 + 2.0 * r * h2 * h) / B);

    beta_grad += h * u_dot_g + (rf->ndim - 1) * h / A + (h + h_prime * r) / B;

    /* coord-wise gradients for center */
    for (d = 0; d < rf->ndim; d++) {
      int idx = i*rf->ndim + d;
      double ctrg = -rf->beta * h * vec_get(origgrad, idx) + rf->beta * h2/r * u_dot_g * u[d]
        - dF_dr / r * u[d];
      vec_set(ctr_grad, d, vec_get(ctr_grad, d) + ctrg);
    }
  }
  free(u);

  /* finally populate paramgrad */
  int j = 0;
  for (d = 0; d < rf->ndim; d++)
    vec_set(paramgrad, j++, vec_get(ctr_grad, d));
    
  /* have to apply chain rule for a and b */
  double a_grad = alpha_grad * 1.0 / (1.0 + exp(-rf->a));
  double b_grad = beta_grad * 1.0 / (1.0 + exp(-rf->b));

  vec_set(paramgrad, j++, a_grad);
  vec_set(paramgrad, j++, b_grad);

  vec_free(ctr_grad);
}
