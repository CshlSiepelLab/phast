#ifndef PLANAR_FLOW_H
#define PLANAR_FLOW_H

#include <stdio.h>
#include <math.h>
#include "phast/vector.h"

#define RF_EPS 1.0e-6

/* reuse the same numerically-stable helpers you had */
static inline double softplus(double x) { return log1p(exp(x)); }
static inline double inv_softplus(double sp) { return log(exp(sp) - 1.0); }

/* Planar flow:
   y = x + u * tanh(w^T x + b)
   J = I + u * psi^T,  psi = (1 - tanh(s)^2) * w,  s = w^T x + b
   log|det J| = log|1 + u^T psi|
*/
typedef struct {
  int npoints;      /* number of points (taxa) */
  int ndim;         /* embedding dimension per point */

  /* parameters (shared across points) */
  Vector *u;        /* length ndim */
  Vector *w;        /* length ndim */
  double b;         /* scalar bias */

  /* gradients */
  Vector *u_grad;   /* length ndim */
  Vector *w_grad;   /* length ndim */
  double b_grad;

} PlanarFlow;

PlanarFlow *pf_new(int npoints, int ndim);
void pf_free(PlanarFlow *pf);

/* optional “update” hook to keep parity with radial API (noop for planar) */
static inline void pf_update(PlanarFlow *pf) { (void)pf; }

/* forward: y = f(x); returns sum over points of log|det J| */
double pf_forward(PlanarFlow *pf, Vector *y, Vector *x);

/* backprop:
   - origgrad: dL/dy (length npoints*ndim)
   - newgrad:  populated with dL/dx
   - also fills pf->{u,w,b}_grad (accumulated across points)
*/
void pf_backprop(PlanarFlow *pf, Vector *x, Vector *newgrad, Vector *origgrad);

#endif /* PLANAR_FLOW_H */
