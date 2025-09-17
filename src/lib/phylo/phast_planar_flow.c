/* this version written by ChatGPT with limited hand editing */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "phast/vector.h"
#include "phast/misc.h"
#include "phast/planar_flow.h"

PlanarFlow *pf_new(int npoints, int ndim) {
  PlanarFlow *pf = smalloc(sizeof(PlanarFlow));
  pf->npoints = npoints;
  pf->ndim = ndim;

  pf->u = vec_new(ndim);
  pf->w = vec_new(ndim);
  vec_zero(pf->u);
  vec_zero(pf->w);
  pf->b = 0.0;

  /* tiny random init so the layer is near-identity but not dead */
  for (int d = 0; d < ndim; d++) {
    vec_set(pf->u, d, norm_draw(0, 0.05));
    vec_set(pf->w, d, norm_draw(0, 0.05));
  }

  pf->u_grad = vec_new(ndim); vec_zero(pf->u_grad);
  pf->w_grad = vec_new(ndim); vec_zero(pf->w_grad);
  pf->b_grad = 0.0;

  return pf;
}

void pf_free(PlanarFlow *pf) {
  vec_free(pf->u);
  vec_free(pf->w);
  vec_free(pf->u_grad);
  vec_free(pf->w_grad);
  sfree(pf);
}

/* Forward pass:
   For each point i:
     s = w^T x_i + b
     t = tanh(s)
     y_i = x_i + u * t
     psi = (1 - t^2) * w
     logdet_i = log( 1 + u^T psi )
   Returns sum_i logdet_i
*/
double pf_forward(PlanarFlow *pf, Vector *y, Vector *x) {
  assert(y->size == x->size && x->size == pf->npoints * pf->ndim);
  const int D = pf->ndim;
  double logdet_sum = 0.0;

  for (int i = 0; i < pf->npoints; i++) {
    /* s = w^T x_i + b */
    double s = 0.0;
    for (int d = 0; d < D; d++) {
      int idx = i*D + d;
      s += vec_get(pf->w, d) * vec_get(x, idx);
    }
    s += pf->b;

    const double t  = tanh(s);
    const double dt = 1.0 - t*t;           /* derivative of tanh */

    /* y_i = x_i + u * t */
    for (int d = 0; d < D; d++) {
      int idx = i*D + d;
      double yi = vec_get(x, idx) + vec_get(pf->u, d) * t;
      vec_set(y, idx, yi);
    }

    /* log|det J| = log|1 + u^T psi|,  psi = dt * w */
    double u_dot_w = 0.0;
    for (int d = 0; d < D; d++) u_dot_w += vec_get(pf->u, d) * vec_get(pf->w, d);
    double det_term = 1.0 + dt * u_dot_w;

    /* For training, we expect det_term > 0. Add a tiny floor for stability. */
    if (det_term <= PF_EPS) det_term = PF_EPS;
    logdet_sum += log(det_term);
  }

  return logdet_sum;
}

/* Backprop:
   Given g_y (origgrad), return g_x (newgrad), and accumulate param grads.

   J = I + u psi^T,   psi = dt * w,  s = w^T x + b, t = tanh(s), dt = 1 - t^2
   J^T g = g + psi (u^T g)

   logdet = log(1 + u^T psi) = log(1 + dt * u^T w)
   d/dx logdet = [ (u^T w) * (-2 t dt) / (1 + dt * u^T w) ] * w
   d/du logdet = psi / (1 + dt * u^T w)
   d/dw logdet = [ dt * u + (-2 t dt) * (u^T w) * x ] / (1 + dt * u^T w)
   d/db logdet = [ (-2 t dt) * (u^T w) ] / (1 + dt * u^T w)

   Path (through y):
     dL/du += t * g
     dL/dw += (u^T g) * dt * x
     dL/db += (u^T g) * dt
     g_x   = g + psi * (u^T g)
*/
void pf_backprop(PlanarFlow *pf, Vector *x, Vector *newgrad, Vector *origgrad) {
  assert(x->size == newgrad->size && x->size == origgrad->size);
  assert(x->size == pf->npoints * pf->ndim);

  const int D = pf->ndim;
  vec_zero(pf->u_grad);
  vec_zero(pf->w_grad);
  pf->b_grad = 0.0;

  for (int i = 0; i < pf->npoints; i++) {
    /* s = w^T x_i + b */
    double s = 0.0;
    for (int d = 0; d < D; d++) {
      int idx = i*D + d;
      s += vec_get(pf->w, d) * vec_get(x, idx);
    }
    s += pf->b;

    const double t  = tanh(s);
    const double dt = 1.0 - t*t;

    /* u^T g_y and w^T g_y (the latter only for an alternative form; not needed) */
    double u_dot_gy = 0.0;
    /*    double w_dot_gy = 0.0; */
    for (int d = 0; d < D; d++) {
      int idx = i*D + d;
      double gy = vec_get(origgrad, idx);
      u_dot_gy += vec_get(pf->u, d) * gy;
      /*      w_dot_gy += vec_get(pf->w, d) * gy; */ /* not used below, but handy if you refactor */
    }

    /* psi = dt * w, u^T psi, denom for logdet grads */
    double u_dot_w = 0.0;
    for (int d = 0; d < D; d++) u_dot_w += vec_get(pf->u, d) * vec_get(pf->w, d);
    double denom = 1.0 + dt * u_dot_w;
    if (denom <= PF_EPS) denom = PF_EPS;

    /* g_x = J^T g_y + d/dx logdet */
    /* J^T g_y = g_y + psi * (u^T g_y) */
    for (int d = 0; d < D; d++) {
      int idx = i*D + d;
      double gy = vec_get(origgrad, idx);
      double psi_d = dt * vec_get(pf->w, d);
      double gx = gy + psi_d * u_dot_gy;

      /* + d/dx logdet */
      double coeff = (-2.0 * t * dt * u_dot_w) / denom; /* scalar */
      gx += coeff * vec_get(pf->w, d);

      vec_set(newgrad, idx, gx);
    }

    /* parameter grads: path + logdet contributions */
    for (int d = 0; d < D; d++) {
      int idx = i*D + d;

      /* dL/du_d: path t * g_y_d + logdet psi_d / denom */
      double gy = vec_get(origgrad, idx);
      double psi_d = dt * vec_get(pf->w, d);
      double du = t * gy + psi_d / denom;
      vec_set(pf->u_grad, d, vec_get(pf->u_grad, d) + du);

      /* dL/dw_d: path (u^T g_y) * dt * x_d + logdet term */
      double x_d = vec_get(x, idx);
      double dw = (u_dot_gy * dt) * x_d
                + ( (dt * vec_get(pf->u, d)) + (-2.0 * t * dt) * u_dot_w * x_d ) / denom;
      vec_set(pf->w_grad, d, vec_get(pf->w_grad, d) + dw);
    }

    /* dL/db: path (u^T g_y) * dt + logdet term */
    pf->b_grad += (u_dot_gy * dt) + ( (-2.0 * t * dt) * u_dot_w ) / denom;
  }
}
