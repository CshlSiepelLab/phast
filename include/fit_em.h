/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file fit_em.h
   Function for fitting tree models with EM
   @ingroup phylo
*/ 


#ifndef FIT_EM_H
#define FIT_EM_H

#include <tree_model.h>
#include <msa.h>
#include <vector.h>
#include <numerical_opt.h>

/** 
   Fit tree model using EM algorithm.
   @param mod Tree Model to fit to
   @param msa Multiple Sequence Alignment sequence data
   @param params Parameters to fit
   @param cat Site category in MSA
   @param precision Precision with which to fit to model
   @param logf output file to write to
   @param error_file If non-NULL, write estimate, variance, and 95% 
   confidence interval for each parameter to this file.
*/
int tm_fit_em(TreeModel *mod, MSA *msa, Vector *params, int cat, 
              opt_precision_type precision, int max_its, FILE *logf,
	      FILE *error_file);

#endif
