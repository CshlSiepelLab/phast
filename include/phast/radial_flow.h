#ifndef RADIAL_FLOW_H
#define RADIAL_FLOW_H

#include <stdio.h>
#include <math.h>
#include "phast/vector.h"

#define RF_EPS 1.0e-6

static inline double softplus(double x) {
  return log(1 + exp(x));
}

static inline double inv_softplus(double sp) {
  return log(exp(sp)-1);
}

typedef struct {
  int npoints; /* flow applies to npoints points, each in ndim dimensions */
  int ndim;
  Vector *ctr; /* center of flow, another point in ndim dimensions */
  double a; /* raw parameter, unconstrained */
  double b; /* raw parameter, unconstrained */
  double alpha; /* alpha = softplus(a) + eps */
  double beta; /* beta = softplus(b) + eps */
  Vector *ctr_grad; /* gradients */
  double a_grad;
  double b_grad;
  double r_med; /* used to keep scale-free */
  unsigned int center_update; /* can be used to turn off update of
                                 center; currently left on by
                                 default */
} RadialFlow;

RadialFlow *rf_new(int npoints, int ndim);

void rf_free(RadialFlow *rf);

void rf_update(RadialFlow *rf);

void rf_rescale(RadialFlow *rf, double scale);

double rf_forward (RadialFlow *rf, Vector *y, Vector *x);

void rf_backprop(RadialFlow *rf, Vector *x, Vector *newgrad,
                 Vector *origgrad);

#endif
