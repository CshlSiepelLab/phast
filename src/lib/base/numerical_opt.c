/* $Id: numerical_opt.c,v 1.2 2004-06-11 05:58:51 acs Exp $
   Written by Adam Siepel, 2002
   Copyright 2002, Adam Siepel, University of California */

#include <stdlib.h>
#include <numerical_opt.h>
#include <matrix.h>
#include <markov_matrix.h>
#include <math.h>
#include <misc.h>
#include <assert.h>
#include <sys/time.h>

/* Numerical optimization of multidimensional functions by the
   "variable metric" or "quasi-Newton" Broyden-Fletcher-Goldfarb-
   Shanno (BFGS) algorithm; see Numerical Recipes in C (Press et al),
   section 10.7 (Second Edition, 1992, as reprinted 2002). */

#define DERIV_EPSILON 1e-6      /* for numerical computation of
                                   derivatives */

#define EPS 1.11e-16            /* approx machine precision */

#define TOLX_HIGH (4*EPS)       /* convergence criterion for param */
#define TOLX_MED 1.0e-8         /* vals (high, medium, and low */
#define TOLX_LOW 1.0e-6         /* precision versions) */

#define TOLX(P) ( (P) == OPT_HIGH_PREC ? TOLX_HIGH : ( (P) == OPT_MED_PREC ? TOLX_MED : TOLX_LOW) )

/* below are definitions for some alternative convergence criteria.
   They are primarily useful for approximate (low-precision)
   estimation, for which the TOLX measure does not seem well suited */

/* convergence criteria in terms of minimum number of stable
   significant figures of parameter estimates */
#define SIGFIG_HIGH 4
#define SIGFIG_MED 3
#define SIGFIG_LOW 2

#define SIGFIG(P) ( (P) == OPT_HIGH_PREC ? SIGFIG_HIGH : ( (P) == OPT_MED_PREC ? SIGFIG_MED : SIGFIG_LOW) )

/* convergence criteria in terms of changes to function value.  These
   are tuned roughly for phylogenetic models and may need to be
   adjusted for other problems. */

#define DELTA_FUNC_HIGH 1e-9    /* (f(x2) - f(x1)) / f(x2) */
#define DELTA_FUNC_MED 1e-7
#define DELTA_FUNC_LOW 1e-6

#define DELTA_FUNC(P) ( (P) == OPT_HIGH_PREC ? DELTA_FUNC_HIGH : ( (P) == OPT_MED_PREC ? DELTA_FUNC_MED : DELTA_FUNC_LOW) )

#define LAMBDA_THRESHOLD 0.01   /* SIGFIG and DELTA_FUNC convergence
                                   criteria will not be applied when a
                                   line search in the Newton direction
                                   is radically truncated, as
                                   defined by this threshold (see
                                   opt_lnsrch, below)  */


#define ITMAX 200               /* maximum allowed number of
                                   iterations (opt_bfgs) */

#define STEP_SCALE 100          /* scale factor for maximum step size
                                   in line searches (opt_bfgs) */

#define ALPHA 1.0e-4            /* threshold for sufficient decrease
                                   in function value (lnsrch) */

#define GTOL 1.0e-5             /* convergence criterion for zeroing
                                   the gradient */

#define BOUNDARY_EPS 1.0e-3     /* FIXME: larger? */

#ifdef DEBUG
FILE *debugf = NULL;
#endif

/* Numerically compute the gradient for the specified function at the
   specified parameter values.  Vector "grad" must already be
   allocated.  Will pass on to the specified function the auxiliary
   parameter "data" as well as the adjusted parameter vector.  The
   parameter "method" may be "OPT_DERIV_CENTRAL", "OPT_DERIV_FORWARD",
   or "OPT_DERIV_BACKWARD"; in the latter two cases, the specified
   "reference_val" will be used, and one function evaluation per
   parameter will be saved.  If either of "lower_bounds" and
   "upper_bounds" is non-NULL, each parameter will be tested against
   the specified bounds, and the selected derivative method will be
   overridden as necessary, to avoid stepping "out of bounds". */
void opt_gradient(gsl_vector *grad, double (*f)(gsl_vector*, void*), 
                  gsl_vector *params, void* data, opt_deriv_method method,
                  double reference_val, gsl_vector *lower_bounds, 
                  gsl_vector *upper_bounds) {
  int i;
  double val1, val2;

  for (i = 0; i < params->size; i++) {
    double origparm = gsl_vector_get(params, i);
    double delta = 2 * DERIV_EPSILON;

    if (method == OPT_DERIV_FORWARD ||
        (lower_bounds != NULL && 
         origparm - gsl_vector_get(lower_bounds, i) < DERIV_EPSILON)) {
      delta = DERIV_EPSILON;
      val1 = reference_val;
    }
    else {
      gsl_vector_set(params, i, origparm - DERIV_EPSILON);
      val1 = f(params, data);
    }

    if (method == OPT_DERIV_BACKWARD || 
        (upper_bounds != NULL && 
         gsl_vector_get(upper_bounds, i) - origparm < DERIV_EPSILON)) {
      delta = DERIV_EPSILON;
      val2 = reference_val;
    }
    else {
      gsl_vector_set(params, i, origparm + DERIV_EPSILON);
      val2 = f(params, data);
    }

    gsl_vector_set(grad, i, (val2 - val1) / delta);
    gsl_vector_set(params, i, origparm);
  }
}

/* Test each parameter against specified bounds, and set "at_bounds"
   accordingly (every element will be given value "OPT_LOWER_BOUND",
   "OPT_UPPER_BOUND", or "OPT_NO_BOUND").  Either or both boundary
   vectors may be NULL, causing no parameter to be considered at the
   boundary.  If the gradient is non-NULL, then a parameter will be
   considered to be at a boundary *only if* the gradient is in the
   "wrong" direction (that is, a direction pushing the parameter out
   of bounds).  If the "only_add" parameter is set to 1, then
   parameters already at a boundary will be left alone.  The return
   value is equal to the number of params that have been determined to
   be at a bound (if only_add == 1, then only newly identified ones
   will be counted).  This is a support function for opt_bfgs.  */
inline int test_bounds(gsl_vector *params, gsl_vector *grad, 
                       gsl_vector *lower_bounds, gsl_vector *upper_bounds,
                       gsl_vector *at_bounds, int only_add) {
  int i;
  int retval = 0;
  assert((grad == NULL || params->size == grad->size) && 
         params->size == at_bounds->size && 
         (lower_bounds == NULL || lower_bounds->size == params->size) && 
         (upper_bounds == NULL || upper_bounds->size == params->size));

  /* short circuit if only_add == 1 and no bounds are defined */
  if (only_add && lower_bounds == NULL && upper_bounds == NULL)
    return 0;

  for (i = 0; i < params->size; i++) {
    double p;

    if (only_add == 1 && gsl_vector_get(at_bounds, i) != OPT_NO_BOUND)
      continue;

    p = gsl_vector_get(params, i);
    if (lower_bounds != NULL && 
        p - gsl_vector_get(lower_bounds, i) < BOUNDARY_EPS && 
        (grad == NULL || gsl_vector_get(grad, i) > 0)) {
/*         (grad == NULL || gsl_vector_get(grad, i) < 0)) { */
      gsl_vector_set(at_bounds, i, OPT_LOWER_BOUND);
      retval++;
#ifdef DEBUG
      fprintf(debugf, "\nParameter %d detected at boundary.\n", i);
#endif
    }
    else if (upper_bounds != NULL && 
             gsl_vector_get(upper_bounds, i) - p < BOUNDARY_EPS && 
             (grad == NULL || gsl_vector_get(grad, i) < 0)) {
/*              (grad == NULL || gsl_vector_get(grad, i) > 0)) { */
      gsl_vector_set(at_bounds, i, OPT_UPPER_BOUND);
      retval++;
#ifdef DEBUG
      fprintf(debugf, "Parameter %d at boundary.\n", i);
#endif
    }
    else
      gsl_vector_set(at_bounds, i, OPT_NO_BOUND);
  }
  return retval;
}

/* Simulate a projection of the specified matrix to a
   lower-dimensional space by zeroing all rows and columns
   corresponding to parameters at a boundary.  For use in opt_bfgs. */
inline void project_matrix(gsl_matrix *M, gsl_vector *at_bounds) {
  int i, j;
  assert(M->size1 == at_bounds->size && M->size2 == at_bounds->size);
  for (i = 0; i < at_bounds->size; i++) {
    if (gsl_vector_get(at_bounds, i) != OPT_NO_BOUND) {
      for (j = 0; j < M->size1; j++) {
        gsl_matrix_set(M, j, i, 0);
        gsl_matrix_set(M, i, j, 0);
      }
    }
  }
}

/* Simulate a projection of the specified vector to a
   lower-dimensional space by zeroing all elements corresponding to
   parameters at a boundary.  For use in opt_bfgs. */
inline void project_vector(gsl_vector *v, gsl_vector *at_bounds) {
  int i;
  assert(v->size == at_bounds->size);
  for (i = 0; i < at_bounds->size; i++) 
    if (gsl_vector_get(at_bounds, i) != OPT_NO_BOUND) 
      gsl_vector_set(v, i, 0);
}

/* Scale a line vector, if necessary, such that it will cause no
   parameter to be updated beyond a boundary; return the index of the
   limiting parameter, or -1. */
inline int scale_for_bounds(gsl_vector *linev, gsl_vector *params, 
                            gsl_vector *lower_bounds, 
                            gsl_vector *upper_bounds) {
  int i;
  double minscale = 1;
  int retval = -1;
  
  if (lower_bounds == NULL && upper_bounds == NULL) return -1;
  for (i = 0; i < params->size; i++) {
    double scale1 = 1, scale2 = 1;
    if (lower_bounds != NULL && 
        gsl_vector_get(params, i) + gsl_vector_get(linev, i) <
        gsl_vector_get(lower_bounds, i) && gsl_vector_get(linev, i) != 0)
      scale1 = (gsl_vector_get(params, i) - gsl_vector_get(lower_bounds, i) - EPS) /
        -gsl_vector_get(linev, i);
    if (upper_bounds != NULL && 
        gsl_vector_get(params, i) + gsl_vector_get(linev, i) >
        gsl_vector_get(upper_bounds, i) && gsl_vector_get(linev, i) != 0) 
      scale2 = (gsl_vector_get(upper_bounds, i) - gsl_vector_get(params, i) - EPS) /
        gsl_vector_get(linev, i);

    if (scale1 < minscale) { minscale = scale1; retval = i; }    
    if (scale2 < minscale) { minscale = scale2; retval = i; }
  }
  if (minscale < 1)
    gsl_vector_scale(linev, minscale);

  return retval;
}

#ifdef DEBUG
/* for debugging: check that H is positive definite and print to log file  */
inline void check_H(gsl_matrix *H, gsl_vector *at_bounds){
  int nfree = 0, pos_def, i, j, k, l;
  gsl_matrix *H_proj;
  gsl_vector *evals;

  /* first print H */
  fprintf(debugf, "\nUpdated inverse Hessian:\n");
  gsl_matrix_pretty_print(debugf, H);

  /* verify H has correct dimensionality */
  for (i = 0; i < H->size1; i++) {
    if (gsl_vector_get(at_bounds, i) == OPT_NO_BOUND) nfree++;
    for (j = i; j < H->size2; j++) 
      assert((gsl_vector_get(at_bounds, i) == OPT_NO_BOUND &&
              gsl_vector_get(at_bounds, j) == OPT_NO_BOUND) ||
             (gsl_matrix_get(H, i, j) == 0 && gsl_matrix_get(H, j, i) == 0));
  }
      
  /* create a "proper" projection of H */
  H_proj = gsl_matrix_alloc(nfree, nfree);
  k = 0;
  for (i = 0; i < H->size1; i++) {
    if (gsl_vector_get(at_bounds, i) != OPT_NO_BOUND) continue;
    l = 0;
    for (j = 0; j < H->size2; j++) {
      if (gsl_vector_get(at_bounds, j) != OPT_NO_BOUND) continue;
      gsl_matrix_set(H_proj, k, l, gsl_matrix_get(H, i, j));
      l++;
    }
    k++;
  }

  /* check to see if positive definite */
  evals = gsl_vector_alloc(H_proj->size1);
  mat_eigenvals(H_proj, evals); 
  pos_def = 1;
  fprintf(debugf, "\nEigenvalues: ");
  for (i = 0; i < evals->size; i++) {
    fprintf(debugf, "%f ", gsl_vector_get(evals, i));
    if (gsl_vector_get(evals, i) <= 0)
      pos_def = 0;
  }
  if (pos_def)
    fprintf(debugf, "\nH is positive definite\n\n");
  else
    fprintf(debugf, "\nH is *not* positive definite\n\n");

  gsl_matrix_free(H_proj);
  gsl_vector_free(evals);
}
#endif


/* Find a minimum of the specified function with the BFGS algorithm,
   starting at the specified parameter values.  The implementation
   here closely follows the one presented in Press et al.  The
   parameter "f" is the function to be minimized, "params" are the
   parameters to the function, "data" is auxiliary data to pass the
   function, "retval" will be assigned the function value at the
   minimum, upper_bounds and lower_bounds define boundaries for
   parameters (vector pointers of NULL indicate no boundary).  A log
   will be written to "logf" (set to NULL to disable logging).  On
   exit, "params" will hold the values that minimize the function.
   Returns a non-zero value on error.  

   NOTE: I have layered on top of Press et al.'s BFGS algorithm the
   capability to obey parameter bounds, piecing together a strategy
   from the PAML source code (which has a sparsely-documented
   implementation of BFGS with bounds) and from some reading on the
   web.  Instead of explicitly moving between the complete
   n-dimensional space (for n parameters) and a reduced space in which
   parameters at the boundary are fixed, I simulate the reduction in
   dimension by zeroing out certain rows and columns of matrices and
   elements of vectors.  This strategy avoids some complexity in
   coding.  */

/* NOTE: added optional gradient function, to be used instead of
   opt_gradient if non-NULL */
int opt_bfgs(double (*f)(gsl_vector*, void*), gsl_vector *params, 
             void *data, double *retval, gsl_vector *lower_bounds, 
             gsl_vector *upper_bounds, FILE *logf,
             void (*compute_grad)(gsl_vector *grad, gsl_vector *params,
                                  void *data, gsl_vector *lb, gsl_vector *ub),
             opt_precision_type precision, gsl_matrix *inv_Hessian) {
  
  int check, i, its, n = params->size, success = 0, nevals = 0, 
    params_at_bounds = 0, new_at_bounds, changed_dimension = 0,
    trunc, already_failed = 0, minsf;
  double den, fac, fae, fval, stpmax, temp, test, lambda, fval_old;
  gsl_vector *dg, *g, *hdg, *params_new, *xi, *at_bounds;
  gsl_matrix *H, *first_frac, *sec_frac, *bfgs_term;
  opt_deriv_method deriv_method = OPT_DERIV_FORWARD;
  struct timeval start_time, end_time;

  if (logf != NULL)
    gettimeofday(&start_time, NULL);

  g = gsl_vector_alloc(n);      /* gradient */
  xi = gsl_vector_alloc(n);     /* current direction along which to minimize */
  H = inv_Hessian != NULL ? inv_Hessian : gsl_matrix_alloc(n, n);   
                                /* inverse Hessian */
  dg = gsl_vector_alloc(n);     /* difference between old and new gradients */
  hdg = gsl_vector_alloc(n);    /* H * dg */
  params_new = gsl_vector_alloc(n); /* updated params vector */
  at_bounds = gsl_vector_alloc(n); /* record of whether each param is
                                      at a boundary */
  /* remainder are auxiliary matrices used in update of inverse
     Hessian */  
  first_frac = gsl_matrix_alloc(n, n);
  sec_frac = gsl_matrix_alloc(n, n);
  bfgs_term = gsl_matrix_alloc(n, n);

#ifdef DEBUG
  debugf = fopen("opt.debug", "w+");
#endif

  fval = f(params, data);       /* Calculate starting function value
                                   and gradient, */

  nevals++;

  if (compute_grad != NULL) {
    compute_grad(g, params, data, lower_bounds, upper_bounds);
    nevals++;                   /* here assume equiv of one function
                                   eval -- not necessarily accurate,
                                   but prob. okay approx. */
  }
  else {
    opt_gradient(g, f, params, data, deriv_method, fval, lower_bounds, 
                 upper_bounds);
    nevals += (deriv_method == OPT_DERIV_CENTRAL ? 2 : 1)*params->size;
  }

  /* test bounds of each parameter and set "at_bounds" appropriately */
  params_at_bounds = 
    test_bounds(params, NULL, lower_bounds, upper_bounds, at_bounds, 0);

  /* TODO: report an error if starting parameter actually *out* of bounds */

  /* if appropriate, print header and initial line to log file */
  if (logf != NULL) {
    opt_log(logf, 1, 0, params, g, -1, -1);
    opt_log(logf, 0, fval, params, g, -1, -1);
  }

  /* initialize inv Hessian and direction */
  if (inv_Hessian == NULL) gsl_matrix_set_identity(H);
  gsl_vector_memcpy(xi, g);
  gsl_vector_scale(xi, -1);

  /* if there are parameters at boundaries, reduce the dimensionality
     of xi accordingly */
  if (params_at_bounds) {
    project_vector(xi, at_bounds);
    changed_dimension = 1;
  }
/*   for (i = 0; i < at_bounds->size; i++)  */
/*     gsl_vector_set(at_bounds, i, OPT_NO_BOUND); */

  stpmax = STEP_SCALE * max(vect_norm(params), n);

  for (its = 0; its < ITMAX; its++) { /* main loop */

    /* see if any parameters are (newly) at a boundary, and update
       total number at boundary */
    /* FIXME: should this be here? */
    gsl_vector_scale(xi, -1);   /* temporary hack */
    if ((new_at_bounds = test_bounds(params, xi, lower_bounds, upper_bounds, 
                                     at_bounds, 1)) > 0) {
      params_at_bounds += new_at_bounds;
      project_vector(xi, at_bounds); 
      changed_dimension = 1;
    }
    gsl_vector_scale(xi, -1);

#ifdef DEBUG
    fprintf(debugf, "BFGS, iteration %d\nParameters at a boundary (prior to linesearch): ", its) ;
    if (params_at_bounds == 0) fprintf(debugf, "None\n\n");
    else {
      for (i = 0; i < at_bounds->size; i++)
        if (gsl_vector_get(at_bounds, i) != OPT_NO_BOUND)
          fprintf(debugf, "%d ", i);
      fprintf(debugf, "\n\n");
    }
#endif

    /* scale the direction vector such that the update will send no
       parameter out of bounds */
    trunc = scale_for_bounds(xi, params, lower_bounds, upper_bounds); 

    /* minimize along xi */
    opt_lnsrch(params, fval, g, xi, params_new, retval, stpmax, &check, 
               f, data, &nevals, &lambda); 
    /* function is evaluated in opt_lnsrch, value is returned in
       retval.  We'll ignore the value of "check" here (see Press, et
       al.) */
    fval_old = fval;
    fval = *retval;

    /* update line direction and current version of params */
    gsl_vector_memcpy(xi, params_new);
    gsl_vector_sub(xi, params);
    minsf = opt_min_sigfig(params, params_new);  /* first grab min
                                                    stable sig figs */
    gsl_vector_memcpy(params, params_new);

#ifdef DEBUG
    /* verify xi has correct dimensionality */
    for (i = 0; i < params->size; i++)
      assert(gsl_vector_get(at_bounds, i) == OPT_NO_BOUND ||
             gsl_vector_get(xi, i) == 0);
#endif

    /* test for convergence. "test" will be set to max_i
       abs(xi_i/max(abs(params_i, 1))) -- an adjusted version of the
       largest absolute value in xi (adjusted for large vals in
       params).  Convergence is taken to have occurred if "test" is
       smaller than the threshold TOLX */
    /* NOTE: here "xi" is the difference between the old and new
       parameter vectors */
    test = 0;                   
    for (i = 0; i < n; i++) {
      temp = fabs(gsl_vector_get(xi, i))/
        max(fabs(gsl_vector_get(params, i)), 1.0);
      if (temp > test) test = temp;
    }
    if (test < TOLX(precision)) {
      if (logf != NULL) fprintf(logf, "Convergence via TOLX\n");
      success = 1;
      break;
    }

    /* alternative tests for convergence */
    /* FIXME: also test for radical truncation due to bounds? */
    if (lambda > LAMBDA_THRESHOLD && 
        minsf >= SIGFIG(precision)) {
      if (logf != NULL) fprintf(logf, "Convergence via sigfigs (%d)\n", minsf);
      success = 1;
      break;
    }

    if (lambda > LAMBDA_THRESHOLD && 
        fabs((fval_old - fval) / fval) < DELTA_FUNC(precision)) {
      if (logf != NULL) fprintf(logf, "Convergence via delta func\n");
      success = 1;
      break;
    }

    /* save the old gradient and obtain the new gradient */
    /* first update the method for derivatives, if necessary */
    if (deriv_method == OPT_DERIV_FORWARD && trunc == -1 && lambda == 1)
      deriv_method = OPT_DERIV_CENTRAL; 
                                /* this heuristic seems pretty good:
                                   if we're taking complete Newton
                                   steps, then we're getting close to
                                   the minimum, so we'll start to
                                   obtain more precise (but expensive)
                                   gradient estimates */
    else if (deriv_method == OPT_DERIV_CENTRAL && (trunc != -1 || 
                                                   lambda < 1))
      deriv_method = OPT_DERIV_FORWARD;
                                /* we also need to be able to switch
                                   back, in case we have a lucky
                                   step early in the search */
    gsl_vector_memcpy(dg, g);
    if (compute_grad != NULL) {
      compute_grad(g, params, data, lower_bounds, upper_bounds);
      nevals++;
    }
    else {
      opt_gradient(g, f, params, data, deriv_method, fval, lower_bounds, 
                   upper_bounds);
      nevals += (deriv_method == OPT_DERIV_CENTRAL ? 2 : 1)*params->size;
    }

    if (logf != NULL) 
      opt_log(logf, 0, fval, params, g, trunc, lambda);

#ifdef DEBUG
    opt_log(debugf, 1, 0, params, g, -1, -1);
    opt_log(debugf, 0, fval, params, g, trunc, lambda);
#endif

    /* test for convergence via zero gradient.  here "test" is an
       adjusted version of the largest absolute value in g (the
       gradient) */
    test = 0;                 
    den = max(*retval, 1.0);
    for (i = 0; i < n; i++) {
      temp = fabs(gsl_vector_get(g, i)) * 
        max(fabs(gsl_vector_get(params, i)), 1.0) / den;
      if (temp > test) test = temp;
    }
    if (test < GTOL) {
      success = 1;
      break;
    }

    /* compute difference of gradients */
    gsl_vector_scale(dg, -1);
    gsl_vector_add(dg, g);    

    /* see if any parameters are (newly) at a boundary, and update
       total number at boundary */
    if ((new_at_bounds = test_bounds(params, g, lower_bounds, upper_bounds, 
                                     at_bounds, 1)) > 0) {
      params_at_bounds += new_at_bounds;
      project_vector(xi, at_bounds); 
      changed_dimension = 1;
    }

    /* see about relaxing some constraints (enlarging H) */
    if (params_at_bounds > 0) {
      double max_grad = 0;      /* maximum gradient component in the
                                   "right" direction, for parameters
                                   at the boundary (i.e., suggesting a
                                   move into the permitted region of
                                   the parameter space) */
      int max_idx;              /* corresponding index */
      for (i = 0; i < at_bounds->size; i++) {
        if (gsl_vector_get(at_bounds, i) == OPT_LOWER_BOUND && 
            -1 * gsl_vector_get(g, i) > max_grad) {
          /*             gsl_vector_get(xi, i) > max_grad) { */
          max_idx = i; 
          max_grad = -1 * gsl_vector_get(g, i);
          /*           max_grad = gsl_vector_get(xi, i); */
        }
        else if (gsl_vector_get(at_bounds, i) == OPT_UPPER_BOUND && 
                 gsl_vector_get(g, i) > max_grad) {
          /*                  -1 * gsl_vector_get(xi, i) > max_grad) { */
          max_idx = i; 
          max_grad = gsl_vector_get(g, i);
          /*           max_grad = -1 * gsl_vector_get(xi, i); */
        }            
      }
      /* enlarge in the dimension of max_grad, if its value is
         sufficiently large */
      if (max_grad > 0) {       /* FIXME: zero? */
#ifdef DEBUG
        fprintf(debugf, "\nParameter %d now free; will be included in update.\n", max_idx);
#endif
        gsl_matrix_set(H, max_idx, max_idx, 1);
        gsl_vector_set(at_bounds, max_idx, OPT_NO_BOUND); 
        params_at_bounds--;
        changed_dimension = 1;
      }
    }

    /* FIXME: make this more efficient, once working! */
    /*     if (changed_dimension) { */
    if (params_at_bounds > 0) {
      /* project H, dg, and xi, according to new representation of
         bounds */
      project_matrix(H, at_bounds);
      project_vector(dg, at_bounds);
      project_vector(xi, at_bounds); /* necessary? */
#ifdef DEBUG
      fprintf(debugf, "Dimensionality changed and params at bounds; re-projecting H, dg, and xi.\n");
#endif
    } 
    /* NOTE: if dimensionality has been *increased* such that there
       are no params currently at the bounds, then we need not do
       anything (H already taken care of, dg and xi can remain as
       they are).  */
    changed_dimension = 0;
    /*     } */

    /* compute product of current inv Hessian and difference in
       gradients.  NOTE: if one or more parameters are at a boundary,
       all of the following computations will take place in the
       lower-dimensional space, by virtue of H, dg, and xi having been
       "projected". */
    mat_vect_mult(hdg, H, dg); 

#ifdef DEBUG
    fprintf(debugf, "H:\n");
    gsl_matrix_pretty_print(debugf, H);
    fprintf(debugf, "g:\n");
    gsl_vector_fprintf(debugf, g, "%f");
    fprintf(debugf, "xi:\n");
    gsl_vector_fprintf(debugf, xi, "%f");
    fprintf(debugf, "dg:\n");
    gsl_vector_fprintf(debugf, dg, "%f");
    fprintf(debugf, "hdg:\n");
    gsl_vector_fprintf(debugf, hdg, "%f");
#endif 

    /* update the inv Hessian, using equation 10.7.8 (Numerical
       Recipes in C) */

    /* here xi holds the difference in the param vectors (the "x's"),
       dg holds the difference in the gradients, and hdg is the vector
       given by H*dg (where H is the inverse Hessian) */

    fac = vect_dot_prod(dg, xi); /* denom of first fraction */
    fae = vect_dot_prod(dg, hdg); /* denom of second fraction */

#ifdef DEBUG
    fprintf(debugf, "fac = %f, fae = %f\n", fac, fae);
#endif 
    
    if (fac > sqrt(EPS * vect_dot_prod(dg, dg) * vect_dot_prod(xi, xi))) { 
      /* skip update if fac not "sufficiently
         positive" */
      gsl_vector *u, *u2;

      /* first compute first and second fractional terms in 10.7.8 */
      vect_cross_prod(first_frac, xi, xi); 
      gsl_matrix_scale(first_frac, 1/fac);
      vect_cross_prod(sec_frac, hdg, hdg);
      gsl_matrix_scale(sec_frac, -1/fae);
      
      /* now compute the extra term for BFGS (eqs 10.7.9 and 10.7.10) */
      u = xi;                   /* rename for clarity; note that it is
                                   okay now to destroy xi (will soon
                                   recompute) */
      gsl_vector_scale(u, 1/fac); 
      u2 = hdg;                 /* similarly, use hdg for second term
                                   in u */
      gsl_vector_scale(u2, 1/fae);
      gsl_vector_sub(u, u2);
      vect_cross_prod(bfgs_term, u, u);
      gsl_matrix_scale(bfgs_term, fae);

      /* now compute the updated matrix */
      gsl_matrix_add(H, first_frac);
      gsl_matrix_add(H, sec_frac);
      gsl_matrix_add(H, bfgs_term);

#ifdef DEBUG
/*       check_H(H, at_bounds); */
#endif
    }
    else {
      if (logf != NULL) fprintf(logf, "WARNING: resetting H!\n");
      gsl_matrix_set_identity(H);
      project_matrix(H, at_bounds);
    }


    /* finally, update the direction vector */
    mat_vect_mult(xi, H, g);
    gsl_vector_scale(xi, -1);


    /* make sure direction of xi is okay; if not, reset H to identity;
       on second failure of this type, assume convergence */
    if (vect_dot_prod(g, xi) >= 0) {
      if (already_failed) {
        success = 1; 
        break;
      }
      else {
        if (logf != NULL) fprintf(logf, "WARNING: resetting H! (%s time)\n", 
                already_failed ? "2nd" : "1st");
        gsl_matrix_set_identity(H);
        project_matrix(H, at_bounds);
        mat_vect_mult(xi, H, g);
        gsl_vector_scale(xi, -1);
        already_failed = 1;
      }
    }
  } 

  if (logf != NULL) {
    opt_log(logf, 0, fval, params, g, trunc, lambda); /* final versions */
    gettimeofday(&end_time, NULL);
    fprintf(logf, "\nNumber of iterations: %d\nNumber of function evaluations: %d\nTotal time: %.4f sec.\n", 
            its, nevals, end_time.tv_sec - start_time.tv_sec + 
            (end_time.tv_usec - start_time.tv_usec)/1.0e6);
  }

  gsl_vector_free(dg);
  gsl_vector_free(g);
  gsl_vector_free(hdg);
  gsl_vector_free(params_new);
  gsl_vector_free(xi);
  gsl_vector_free(at_bounds);
  if (inv_Hessian == NULL) gsl_matrix_free(H);
  gsl_matrix_free(first_frac);
  gsl_matrix_free(sec_frac);
  gsl_matrix_free(bfgs_term);

  if (success == 0) {
    fprintf(stderr, 
            "ERROR: exceeded maximum number of iterations in opt_bfgs.\n");
    return 1;
  }
  
  return 0;    
}

/* Given a point "xold", the value of the function "f" and its gradient
   "g" at that point, and a direction "p", find a new point "x" along
   "p" from "xold" such that "f" is minimized or has "decreased
   sufficiently" (see Press et al., pp 383-386).  "p" is generally the
   Newton direction, or an approximation of it (as in the BFGS
   algorithm).  This routine ensures progress will occur toward
   convergence on every iteration, even far from a solution (unlike in
   a straightforward implementation of Newton's method).  The idea is
   to try the full Newton step (when close to the solution it will
   be optimal) but to require that any step taken must reduce
   the function value; if it does not, "backtrack" along "p" until an
   acceptable point is reached.  One key insight is that "f" need not
   be absolutely minimized along "p"; indeed, minimizing it can be
   extremely wasteful of function evaluations.  Instead, any new point
   along "p" is acceptable as long as two criteria are met: the
   average rate of decrease of the function must be at least a
   fraction alpha of the initial rate of decrease (given by the
   gradient), and the rate of decrease at the new point be greater
   than a fraction beta of the initial rate of decrease.  The first
   must be imposed explicitly but the second evidently falls out of the
   backtracking strategy.  See Press et al. for details.

   The parameter "stpmax" limits the length of the step (to help avoid
   sending the function into undefined regions).  The value
   "check_convergence" is true when the routine terminates due to a
   new point falling too close to the old one.  In a minimization
   algorithm this apparently tends to happen only at convergence, and
   thus the flag can be ignored; but in a zero-finding algorithm it
   should be examined.  The parameter "nevals" will incremented each
   time the function is evaluated.  Function returns the value of
   the final scaling factor for "p" (lambda). */
void opt_lnsrch(gsl_vector *xold, double fold, gsl_vector *g, gsl_vector *p, 
                gsl_vector *x, double *f, double stpmax, 
                int *check_convergence, double (*func)(gsl_vector*, void*), 
                void *data, int *nevals, double *final_lambda) {

  int i, n = xold->size;
  double a, lambda, lambda2, lamda_min, b, disc, f2, rhs1, rhs2, slope, sum,
    temp, test, tmplam;
  *check_convergence = 0;

  /* scale step if necessary */
  sum = vect_norm(p);
  if (sum > stpmax)  
    gsl_vector_scale(p, stpmax/sum);

  slope = vect_dot_prod(g, p);
  if (slope >= 0.0) {
    fprintf(stderr, "ERROR: positive slope in opt_lnsrch.  Roundoff error?\n");
    gsl_vector_memcpy(x, xold);
    *check_convergence = 1;
    *final_lambda = lambda;
    return;
  }

  test = 0.0;                   /* compute lambda_min */
  for (i = 0; i < n; i++) {
    temp = fabs(gsl_vector_get(p, i)) / 
      max(fabs(gsl_vector_get(xold, i)), 1.0);
    if (temp > test) test = temp;
  }
  lamda_min = TOLX(OPT_HIGH_PREC)/test;

  lambda = 1.0;                   /* try full Newton step first. */
  for (;;) {                    /* main loop */

    /* compute candidate "point" (parameter values), along direction p */
    gsl_vector_memcpy(x, p);
    gsl_vector_scale(x, lambda);
    gsl_vector_add(x, xold);

    /* function call */
    *f = func(x, data);
    (*nevals)++;

    if (lambda < lamda_min) {   /* if the loop is terminated this way,
                                   the calling program should verify
                                   convergence (mostly important for
                                   zero-finding; see Press et al.) */
      gsl_vector_memcpy(x, xold);
      *check_convergence = 1;
      *final_lambda = lambda;
      return;
    } 
    else if (*f <= fold + ALPHA*lambda*slope) {
      *final_lambda = lambda;
      return;                   /* sufficient decrease; see eq 9.7.7 */
    }

    else {                      /* backtrack */
      if (lambda == 1.0)        /* first time through loop */
        tmplam = -slope/(2.0*(*f-fold-slope)); 
      else {                    /* subsequent backtracks */
        rhs1 = *f - fold - lambda*slope;
        rhs2 = f2 - fold - lambda2*slope;
        a = (rhs1/(lambda*lambda) - rhs2/(lambda2*lambda2))/
          (lambda-lambda2);
        b = (-lambda2*rhs1/(lambda*lambda)+lambda*rhs2/(lambda2*lambda2))/
          (lambda-lambda2);
        if (a == 0.0) tmplam = -slope/(2.0*b);
        else {
          disc = b*b-3.0*a*slope;
          if (disc < 0.0) tmplam = 0.5*lambda;
          else if (b <= 0.0) tmplam = (-b+sqrt(disc))/(3.0*a);
          else tmplam = -slope/(b+sqrt(disc));
        }
        if (tmplam > 0.5*lambda)
          tmplam = 0.5*lambda;      
      }
    }
    lambda2 = lambda;
    f2 = *f;
    lambda = max(tmplam,0.1*lambda); 
  } 
}

/* Print a line to a log file that describes the state of the
   optimization procedure on a given iteration.  The value of the
   function is output, along with the values of all parameters, and
   the components of the gradient.  If "header_only == 1", an
   appropriate header is printed. */
void opt_log(FILE *logf, int header_only, double val, gsl_vector *params, 
             gsl_vector *grad, int trunc, double lambda) {
  int i;
  char tmp[30];
  if (header_only) {
    fprintf(logf, "%15s ", "f(x)");
    for (i = 0; i < params->size; i++) {
      sprintf(tmp, "x_%d", i);
      fprintf(logf, "%15s ", tmp);
    }
    for (i = 0; i < params->size; i++) {
      sprintf(tmp, "df/dx_%d", i);
      fprintf(logf, "%15s ", tmp);
    }
    fprintf(logf, "%15s %15s\n", "trunc", "lambda");
  }
  else {
    fprintf(logf, "%15.6f ", val);
    for (i = 0; i < params->size; i++) 
      fprintf(logf, "%15.6f ", gsl_vector_get(params, i));
    for (i = 0; i < params->size; i++) 
      fprintf(logf, "%15.6f ", gsl_vector_get(grad, i));
    fprintf(logf, "%15d %15f\n", trunc, lambda);
  }
  fflush(logf);
}

/***************************************************************************
 Implementation of Brent's method; very slightly adapted from
 Numerical Recipes in C.
****************************************************************************/

#define ITMAX_BRENT 100
#define CGOLD 0.3819660
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define ZEPS 1.0e-10
/* Here ITMAX_BRENT is the maximum allowed number of iterations; CGOLD is
   the golden ratio; ZEPS is a small number that protects against
   trying to achieve fractional accuracy for a minimum that happens to
   be exactly zero. */
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
double opt_brent(double ax, double bx, double cx, 
                 double (*f)(double, void*), double tol,
                 double *xmin, void *data, FILE *logf)
     /* Given a function f, and given a bracketing triplet of
        abscissas ax, bx, cx (such that bx is between ax and cx, and
        f(bx) is less than both f(ax) and f(cx)), this routine
        isolates the minimum to a fractional precision of about tol
        using Brent's method. The abscissa of the minimum is returned
        as xmin, and the minimum function value is returned as brent,
        the returned function value. */
{
  int iter;
  double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0;                  /*  This will be the distance moved on
                                    the step before last. */
  a=(ax < cx ? ax : cx);        /* a and b must be in ascending order, 
                                   but input abscissas need not be. */
  b=(ax > cx ? ax : cx);
  x=w=v=bx;                     /*  Initializations... */
  fw=fv=fx=(*f)(x, data);
  if (logf != NULL) 
    fprintf(logf, "opt_brent:\nStarting with x_a = %f, x_b = %f, x_c = %f, f(x_b) = %f\n", ax, bx, cx, fx);
  for (iter=1;iter<=ITMAX_BRENT;iter++) { /* Main program loop. */
    xm=0.5*(a+b);
    tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
    if (fabs(x-xm) <= (tol2-0.5*(b-a))) { /* Test for done here. */
      *xmin=x;
      if (logf != NULL) 
        fprintf(logf, "Returning x_min = %f, f(x_min) = %f\n", x, fx);
      return fx;
    }
    if (fabs(e) > tol1) {       /* Construct a trial parabolic fit. */
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q > 0.0) p = -p;
      q=fabs(q);
      etemp=e;
      e=d;
      if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
        d=CGOLD*(e=(x >= xm ? a-x : b-x));
      /* The above conditions determine the acceptability of the
         parabolic fit. Here we take the golden section step into the
         larger of the two segments. */
      else {
        d=p/q;                  /* Take the parabolic step. */
        u=x+d;
        if (u-a < tol2 || b-u < tol2)
          d=SIGN(tol1,xm-x);
      }
    } else {
      d=CGOLD*(e=(x >= xm ? a-x : b-x));
    }
    u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu=(*f)(u, data);
    if (logf != NULL) {
      fprintf(logf, "u = %f, f(u) = %f\n", u, fu);
      fflush(logf);
    }
    /* This is the one function evaluation per iteration. */
    if (fu <= fx) {             /* Now decide what to do with our
                                   function evaluation.  */
      if (u >= x) a=x; else b=x;
      SHFT(v,w,x,u)             /* Housekeeping follows: */
        SHFT(fv,fw,fx,fu)
        } else {
          if (u < x) a=u; else b=u;
          if (fu <= fw || w == x) {
            v=w;
            w=u;
            fv=fw;
            fw=fu;
          } else if (fu <= fv || v == x || v == w) {
            v=u;
            fv=fu;
          }
        } /* Done with housekeeping. Back for another iteration. */ 
  }
  die("ERROR: exceeded max iterations in brent.\n");
  return -1;                    /* never get here */
}

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
/* Here GOLD is the default ratio by which successive intervals are
   magnified; GLIMIT is the maximum magnification allowed for a
   parabolic-fit step. */
void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, 
            double *fc, double (*func)(double, void*), void *data,
            FILE *logf)
     /* Given a function func, and given distinct initial points ax
        and bx, this routine searches in the downhill direction
        (defined by the function as evaluated at the initial points)
        and returns new points ax, bx, cx that bracket a minimum of
        the function. Also returned are the function values at the
        three points, fa, fb, and fc. */
{
  double ulim,u,r,q,fu,dum;
  *fa=(*func)(*ax, data);
  *fb=(*func)(*bx, data);
  if (logf != NULL)
    fprintf(logf, "opt_mnbrak:\nx_a = %f, f(x_a) = %f\nx_b = %f, f(x_b) = %f\n", 
            *ax, *fa, *bx, *fb);
  if (*fb > *fa) {              /* Switch roles of a and b so that we
                                   can go downhill in the direction
                                   from a to b.*/
    SHFT(dum,*ax,*bx,dum) 
    SHFT(dum,*fb,*fa,dum)
      }
  *cx=(*bx)+GOLD*(*bx-*ax);     /* First guess for c. */
  *fc=(*func)(*cx, data);
  while (*fb > *fc) {           /* Keep returning here until we
                                   bracket. */
    if (logf != NULL)
      fprintf(logf, "x_a = %f, f(x_a) = %f\nx_b = %f, f(x_b) = %f\nx_c = %f, f(x_c) = %f\n", 
              *ax, *fa, *bx, *fb, *cx, *fc);
    r=(*bx-*ax)*(*fb-*fc);      /* Compute u by parabolic
                                   extrapolation from a, b, c. TINY is
                                   used to prevent any possible
                                   division by zero. */
    q=(*bx-*cx)*(*fb-*fa);
    u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
      (2.0*SIGN(max(fabs(q-r),TINY),q-r));
    ulim=(*bx)+GLIMIT*(*cx-*bx);
    /* We wont go farther than this. Test various possibilities: */
    if ((*bx-u)*(u-*cx) > 0.0) { /* Parabolic u is between b and c:
                                    try it. */
      fu=(*func)(u, data);
      if (fu < *fc) {           /* Got a minimum between b and c. */
        *ax=(*bx);
        *bx=u;
        *fa=(*fb);
        *fb=fu;
        return;
      } else if (fu > *fb) {    /* Got a minimum between between a and
                                   u. */
        *cx=u;
        *fc=fu;
        return;
      }
      u=(*cx)+GOLD*(*cx-*bx);   /* Parabolic fit was no use. Use
                                   default magnification. */
      fu=(*func)(u, data);
    } else if ((*cx-u)*(u-ulim) > 0.0) { /* Parabolic fit is between c
                                            and its allowed limit.  */
      fu=(*func)(u, data);
      if (fu < *fc) {
        SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
          SHFT(*fb,*fc,fu,(*func)(u, data))
          }
    } else if ((u-ulim)*(ulim-*cx) >= 0.0) { /* Limit parabolic u to
                                                maximum allowed
                                                value.  */
      u=ulim;
      fu=(*func)(u, data);
    } else {                    /* Reject parabolic u, use default
                                   magnification. */
      u=(*cx)+GOLD*(*cx-*bx);
      fu=(*func)(u, data);
    }
    SHFT(*ax,*bx,*cx,u)         /* Eliminate oldest point and
                                   continue. */
      SHFT(*fa,*fb,*fc,fu)
      }
  if (logf != NULL)
    fprintf(logf, "(final)\nx_a = %f, f(x_a) = %f\nx_b = %f, f(x_b) = %f\nx_c = %f, f(x_c) = %f\n", 
              *ax, *fa, *bx, *fb, *cx, *fc);
}

/***************************************************************************/
/* inline functions -- also defined in numerical_opt.h                     */
/***************************************************************************/

/* given two vectors of consecutive parameter estimates, return the
   minimum number of shared significant figures */
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
