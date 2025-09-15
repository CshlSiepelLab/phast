/* FIXME: impose softmax on both alpha and beta; differentiate between
   raw and transformed */

/* function to create new version.   Initialize center */

/* fucntion to free memory including interior vector */

/* in SGA add a very small L2 prior to discourage extreme values, e.g., lambda * (||a||^2 + ||b||^2) where lambda is 1e-6 to 1e-4 */

/* function to update params based on a and b using softmax */

/* initialize with a about log (e - 1) and b negative */

/* set up an EPS local to this program */

/* need functions for gradients wrt a, b, and c including both main flow and log determinant */

/* then need a function to compute the Jacobian for the chain rule.
   Possibly don't want to express in full matrix form */

typedef struct {
  int npoints;
  int ndim;
  Vector *ctr;
  double a; /* raw parameter, unconstrained */
  double b; /* raw parameter, unconstrained */
  double alpha; /* alpha = softplus(a) + eps */
  double beta; /* beta = -alpha + softplus(b) + eps */
} RadialFlow;

/* assumes points x, computes y = f(x) where x is an input set of
   points and y is transformed by the flow, assumes structured as n x
   d.  returns log determinant of Jacobian of transformation */
double rf_forward (RadialFlow *rf, Vector *y, Vector *x) {
  assert(y->size == x->size && y->size == rf->npoints * rf->ndim);
  double diff[rf->ndim];
  
  for (int i = 0; i < rf->npoints; i++) {
    int idx, d;
    double r2 = 0, r, h, logdet = 0;
    for (d = 0; d < rf->ndim; d++) {
      idx = i*rf->ndim + d;
      diff[d] = vec_get(x, idx) - vec_get(rf->ctr, idx);
      r2 += (diff[d] * diff[d]);
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
}



