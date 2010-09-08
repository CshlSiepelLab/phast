/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: numerical_opt.h,v 1.6 2008-11-15 19:42:03 acs Exp $ */

#ifndef NEWRAPH_H
#define NEWRAPH_H

#include <matrix.h>
#include <vector.h>
#include <math.h>

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

void opt_gradient(Vector *grad, double (*f)(Vector*, void*), 
                  Vector *params, void* data, opt_deriv_method method,
                  double reference_val, Vector *lower_bounds, 
                  Vector *upper_bounds);

int opt_bfgs(double (*f)(Vector*, void*), Vector *params, 
             void *data, double *retval, Vector *lower_bounds, 
             Vector *upper_bounds, FILE *logf,
             void (*compute_grad)(Vector *grad, Vector *params,
                                  void *data, Vector *lb, Vector *ub),
             opt_precision_type precision, Matrix *inv_Hessian);

void opt_lnsrch(Vector *xold, double fold, Vector *g, Vector *p, 
                Vector *x, double *f, double stpmax, 
                int *check_convergence, double (*func)(Vector*, void*), 
                void *data, int *nevals, double *lambda, FILE *logf);

void opt_log(FILE *F, int header_only, double val, Vector *params, 
             Vector *derivs, int trunc, double lambda);


double opt_brent(double ax, double bx, double cx, 
                 double (*f)(double, void*), double tol,
                 double *xmin, void *data, FILE *logf);

void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
            double (*func)(double, void*), void *data, FILE *logf);

int opt_newton_1d(double (*f)(double, void*), double (*x), void *data, 
                  double *fx, int sigfigs, double lb, double ub, FILE *logf, 
                  double (*compute_deriv)(double x, void *data, double lb, 
                                          double ub),
                  double (*compute_deriv2)(double x, void *data, double lb, 
                                           double ub));

void opt_derivs_1d(double *deriv, double *deriv2, double x, double fx, 
                   double lb, double ub, double (*f)(double, void*), void *data,
                   double (*compute_deriv)(double x, void *data, double lb, 
                                           double ub),
                   double (*compute_deriv2)(double x, void *data, double lb, 
                                            double ub));

int opt_bfgs_1d(double (*f)(double, void*), double (*x), void *data, 
                  double *fx, int sigfigs, double lb, double ub, FILE *logf, 
                  double (*compute_deriv)(double x, void *data, double lb, 
                                          double ub));

void opt_lnsrch_1d(double direction, double xold, double fold, double *x, 
                   double *fx, double deriv, double (*func)(double, void*), 
                   void *data, int *nevals, double *final_lambda, FILE *logf);

int opt_min_sigfig(Vector *p1, Vector *p2);

int opt_sigfig(double val1, double val2);

#endif
