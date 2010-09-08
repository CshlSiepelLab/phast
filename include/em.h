/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: em.h,v 1.7 2008-11-12 02:07:59 acs Exp $ */

#ifndef EM_H
#define EM_H

#include <hmm.h>

#define EM_CONVERGENCE_THRESHOLD 0.1 /* TEMPORARY! */

double hmm_train_by_em(HMM *hmm, void *models, void *data, int nsamples, 
                       int *sample_lens, Matrix *pseudocounts, 
                       void (*compute_emissions)(double**, void**, int, void*, 
                                                 int, int), 
                       void (*estimate_state_models)(void**, int, void*, 
                                                     double**, int, FILE*),
                       void (*estimate_transitions)(HMM*, void*, double**),
                       int (*get_observation_index)(void*, int, int),
                       void (*log_function)(FILE*, double, HMM*, void*, int),
                       double **emissions_alloc, FILE *logf);

#endif
