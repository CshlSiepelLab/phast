/* $Id: numerical_opt.h,v 1.1.1.1 2004-06-03 22:43:11 acs Exp $
   Written by Adam Siepel, 2002
   Copyright 2002, Adam Siepel, University of California */

#ifndef NEWRAPH_H
#define NEWRAPH_H

#include "gsl/gsl_matrix.h"
#include "math.h"
#include "assert.h"

typedef enum {
  OPT_DERIV_BACKWARD,
  OPT_DERIV_CENTRAL,
  OPT_DERIV_FORWARD
} opt_deriv_method;

typedef enum {
  OPT_LOWER_BOUND,
  OPT_UPPER_BOUND,
  OPT_NO_BOUND
} opt_at_bound_type;

/* level of precision to use in estimating parameters */
typedef enum {                  
  OPT_LOW_PREC, 
  OPT_MED_PREC, 
  OPT_HIGH_PREC,
  OPT_CRUDE_PREC
} opt_precision_type; 

void opt_gradient(gsl_vector *grad, double (*f)(gsl_vector*, void*), 
                  gsl_vector *params, void* data, opt_deriv_method method,
                  double reference_val, gsl_vector *lower_bounds, 
                  gsl_vector *upper_bounds);

int opt_bfgs(double (*f)(gsl_vector*, void*), gsl_vector *params, 
             void *data, double *retval, gsl_vector *lower_bounds, 
             gsl_vector *upper_bounds, FILE *logf,
             void (*compute_grad)(gsl_vector *grad, gsl_vector *params,
                                  void *data, gsl_vector *lb, gsl_vector *ub),
             opt_precision_type precision, gsl_matrix *inv_Hessian);

void opt_lnsrch(gsl_vector *xold, double fold, gsl_vector *g, gsl_vector *p, 
                gsl_vector *x, double *f, double stpmax, 
                int *check_convergence, double (*func)(gsl_vector*, void*), 
                void *data, int *nevals, double *lambda);

void opt_log(FILE *F, int header_only, double val, gsl_vector *params, 
             gsl_vector *derivs, int trunc, double lambda);


double opt_brent(double ax, double bx, double cx, double (*f)(double, void*), double tol,
                 double *xmin, void *data);

void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
            double (*func)(double, void*), void *data);

int opt_min_sigfig(gsl_vector *p1, gsl_vector *p2);

/***************************************************************************/
/* inline functions -- also defined in numerical_opt.c                     */
/***************************************************************************/

#define MAXSIGFIGS 6


/* given two vectors of consecutive parameter estimates, return the
   minimum number of shared significant figures */
extern inline
int opt_min_sigfig(gsl_vector *p1, gsl_vector *p2) {
  int i, sf, min = 99999;
  assert(p1->size == p2->size);
  for (i = 0; i < p1->size; i++) {
    double val1 = gsl_vector_get(p1, i);
    double val2 = gsl_vector_get(p2, i);
    double tmp = pow(10, floor(log10(val1)));
    val1 /= tmp; val2 /= tmp;
    for (sf = 0; sf < MAXSIGFIGS; sf++) {
      if (floor(val1) != floor(val2)) break;
      val1 *= 10;
      val2 *= 10;
    }    
    if (sf < min) min = sf;
  }
  if (min < 0) min = 0;
  return min;
}

#endif
