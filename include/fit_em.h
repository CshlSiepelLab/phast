/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: fit_em.h,v 1.3 2008-11-12 02:07:59 acs Exp $ */

#ifndef FIT_EM_H
#define FIT_EM_H

#include <tree_model.h>
#include <msa.h>
#include <vector.h>
#include <numerical_opt.h>

int tm_fit_em(TreeModel *mod, MSA *msa, Vector *params, int cat, 
              opt_precision_type precision, FILE *logf);

#endif
