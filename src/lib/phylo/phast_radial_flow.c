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
  rf->a = log(exp(1) - 1); /* initial values that cause near identity transformation */
  rf->b = -10;
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
  double diff[rf->ndim];
  double logdet = 0;
  
  for (int i = 0; i < rf->npoints; i++) {
    int idx, d;
    double r2 = 0, r, h, bh;
    for (d = 0; d < rf->ndim; d++) {
      idx = i*rf->ndim + d;
      diff[d] = vec_get(x, idx) - vec_get(rf->ctr, idx);
      r2 += diff[d] * diff[d];
    }
    r = sqrt(r2);
    h = 1 / (rf->alpha + r);
    for (d = 0; d < rf->ndim; d++) {
      idx = i*rf->ndim + d;
      vec_set(y, idx, vec_get(x, idx) + h * diff[d]);
    }
    bh = rf->beta*h;
    logdet += (rf->ndim - 1) * log(1 + bh) + log(1 + bh - bh*h*r);
  }
  return logdet;
}

/* param grad: first d components are derivs wrt ctr; then deriv wrt alpha and deriv wrt beta */
void rf_backprop(RadialFlow *rf, Vector *x, Vector *newgrad, Vector *origgrad,
                 Vector *paramgrad) {
  double u[rf->ndim];
  assert(x->size == rf->npoints * rf->ndim && x->size == newgrad->size &&
         x->size == origgrad->size);

  for (int i = 0; i < rf->npoints; i++) {
    int idx, d;
    double r2 = 0, u_dot_g = 0, r, h, h2, h_prime, bh, A, B, dF_dr;
    for (d = 0; d < rf->ndim; d++) {
      idx = i*rf->ndim + d;
      u[d] = vec_get(x, idx) - vec_get(rf->ctr, idx);
      r2 += u[d] * u[d];
      u_dot_g += u[d] * vec_get(origgrad, idx);
    }

    r = sqrt(r2);  
    if (r < 1e-18) r = 1e-18;  /* avoid dividing by very small number */
      
    h = 1 / (rf->alpha + r);
    h2 = h*h;
    h_prime = -h2;

    bh = rf->beta * h;
    A = 1 + bh;
    B = 1 + bh + rf->beta * h_prime * r;
    dF_dr = rf->beta * (-(rf->ndim-1) * h2 / A + (-2*h2 + 2*r*h2*h) / B);
      
    for (d = 0; d < rf->ndim; d++) {
      idx = i*rf->ndim + d;
      
      /* path contribution */
      double newg = A * vec_get(origgrad, idx) + rf->beta * h_prime / r * u[d] * u_dot_g;

      /* log determinant contribution */
      newg += dF_dr / r * u[d];
        
      vec_set(newgrad, idx, newg);
    }

    /* populate paramgrad */
  }      
}
