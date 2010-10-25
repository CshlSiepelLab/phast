/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: subst_mods.h,v 1.9.2.1 2009-03-18 19:35:57 mt269 Exp $ */

#ifndef SUBST_MODS_H
#define SUBST_MODS_H

#include <markov_matrix.h>
#include <vector.h>
#include <stringsplus.h>

typedef enum {
  JC69,
  K80,
  F81,
  HKY85,
  HKY85G,
  REV,
  SSREV,
  UNREST,
//  HKY2,
  R2,
  U2,
  R2S,
  U2S,
  R3,
  R3S,
  U3,
  U3S,
  GC,
  //  HB,
  HKY_CODON,
  REV_CODON,
  SSREV_CODON,
  UNDEF_MOD
} subst_mod_type;

struct tm_struct;               /* use incomplete type because of
                                   reciprocal dependencies with tree_model.h */

void tm_set_probs_JC69(struct tm_struct *mod, MarkovMatrix *P, double t);
void tm_set_probs_F81(Vector *backgd_freqs, MarkovMatrix *P, double scale, double t);
void tm_set_probs_independent(struct tm_struct *mod, MarkovMatrix *P);
subst_mod_type tm_get_subst_mod_type(const char *str);
char *tm_get_subst_mod_string(subst_mod_type type);
int tm_get_nratematparams(struct tm_struct *mod);
int tm_order(int subst_mod);
int subst_mod_is_codon_model(int subst_mod);
void tm_set_rate_matrix(struct tm_struct *mod, Vector *params, int i);
void tm_set_rate_matrix_sel_bgc(struct tm_struct *mod, Vector *params, int i,
				double selection, double bgc);
void tm_rate_params_init(struct tm_struct *mod, Vector *params, 
                         int params_idx, double kappa);
void tm_rate_params_init_from_model(struct tm_struct *mod, Vector *params, 
                                    int params_idx, 
				    double selection, double bgc);
void tm_set_HKY_matrix(struct tm_struct *mod, double kappa, int kappa_idx);
int tm_substmod_get_nratematparams(subst_mod_type submod, 
				   struct tm_struct *mod);
int tm_flag_subst_param_pos(struct tm_struct *mod, int *flag, 
			    String *param_name);
void tm_apply_selection_bgc(MarkovMatrix *mm, double sel, double bgc);
void tm_unapply_selection_bgc(MarkovMatrix *mm, double sel, double bgc);

#endif
