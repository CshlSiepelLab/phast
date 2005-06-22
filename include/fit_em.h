/* $Id: fit_em.h,v 1.2 2005-06-22 07:11:20 acs Exp $
   Written by Adam Siepel, 2002-2004
   Copyright 2002-2004, Adam Siepel, University of California */

#ifndef FIT_EM_H
#define FIT_EM_H

#include <tree_model.h>
#include <msa.h>
#include <vector.h>
#include <numerical_opt.h>

int tm_fit_em(TreeModel *mod, MSA *msa, Vector *params, int cat, 
              opt_precision_type precision, FILE *logf);

#endif
