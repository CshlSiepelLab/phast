/* $Id: subst_mods.h,v 1.2 2004-06-19 20:35:14 acs Exp $
   Written by Adam Siepel, 2002-2004
   Copyright 2002-2004, Adam Siepel, University of California */

#ifndef SUBST_MODS_H
#define SUBST_MODS_H

#include <markov_matrix.h>
#include <gsl/gsl_vector.h>
#include <stringsplus.h>

typedef enum {
  JC69,
  K80,
  F81,
  HKY85,
  REV,
  UNREST,
  HKY2,
  R2,
  U2,
  R2S,
  U2S,
  R3,
  R3S,
  U3,
  U3S,
  UNDEF_MOD
} subst_mod_type;

struct tm_struct;               /* use incomplete type because of
                                   reciprocal dependencies with tree_model.h */

void tm_set_probs_JC69(struct tm_struct *mod, MarkovMatrix *P, double t);
void tm_set_probs_F81(struct tm_struct *mod, MarkovMatrix *P, double scale, double t);
subst_mod_type tm_get_subst_mod_type(char *str);
char *tm_get_subst_mod_string(subst_mod_type type);
int tm_get_nratematparams(struct tm_struct *mod);
int tm_order(int subst_mod);
void tm_set_rate_matrix(struct tm_struct *mod, gsl_vector *params, int i);
void tm_rate_params_init(struct tm_struct *mod, gsl_vector *params, 
                         int params_idx, double kappa);
void tm_rate_params_init_from_model(struct tm_struct *mod, gsl_vector *params, 
                                    int params_idx);

#endif
