#ifndef RADIAL_FLOW_H
#define RADIAL_FLOW_H

#include <stdio.h>
#include <math.h>
#include "phast/vector.h"

#define RF_EPS 1.0e-6

static inline double softplus(double x) {
  return log(1 + exp(x));
}

typedef struct {
  int npoints;
  int ndim;
  Vector *ctr;
  double a; /* raw parameter, unconstrained */
  double b; /* raw parameter, unconstrained */
  double alpha; /* alpha = softplus(a) + eps */
  double beta; /* beta = -alpha + softplus(b) + eps */
} RadialFlow;

RadialFlow *rf_new(int npoints, int ndim);

void rf_free(RadialFlow *rf);

void rf_update(RadialFlow *rf);

double rf_forward (RadialFlow *rf, Vector *y, Vector *x);

void rf_backprop(RadialFlow *rf, Vector *x, Vector *newgrad,
                 Vector *origgrad, Vector *paramgrad);

#endif
