/* $Id: fit_em.h,v 1.1.1.1 2004-06-03 22:43:11 acs Exp $
   Written by Adam Siepel, 2002-2004
   Copyright 2002-2004, Adam Siepel, University of California */

#ifndef FIT_EM_H
#define FIT_EM_H

#include <tree_model.h>
#include <msa.h>
#include <gsl/gsl_vector.h>
#include <numerical_opt.h>

int tm_fit_em(TreeModel *mod, MSA *msa, gsl_vector *params, int cat, 
              opt_precision_type precision, FILE *logf);

#endif
