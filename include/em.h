/* $Id: em.h,v 1.2 2004-06-16 06:21:17 acs Exp $
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
                                                     double**, int),
                       int (*get_observation_index)(void*, int, int),
                       FILE *logf);

void compute_emissions_phyhmm(double **emissions, void **models, int nmodels,
                              void *data, int sample, int length);
void estimate_state_models_phyhmm(void **models, int nmodels, void *data, 
                                  double **E, int nobs);
int get_observation_index_phyhmm(void *data, int sample, int position);


#endif
