/* $Id: em.h,v 1.5 2005-05-31 06:38:07 acs Exp $
   Written by Adam Siepel, 2002
   Copyright 2002, Adam Siepel, University of California */


#ifndef EM_H
#define EM_H

#include <hmm.h>

#define EM_CONVERGENCE_THRESHOLD 0.1 /* TEMPORARY! */

double hmm_train_by_em(HMM *hmm, void *models, void *data, int nsamples, 
                       int *sample_lens, gsl_matrix *pseudocounts, 
                       void (*compute_emissions)(double**, void**, int, void*, 
                                                 int, int), 
                       void (*estimate_state_models)(void**, int, void*, 
                                                     double**, int, FILE*),
                       void (*estimate_transitions)(HMM*, void*, double**),
                       int (*get_observation_index)(void*, int, int),
                       void (*log_function)(FILE*, double, HMM*, void*, int),
                       double **emissions_alloc, FILE *logf);

#endif
