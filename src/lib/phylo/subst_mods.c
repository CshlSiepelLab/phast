/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: subst_mods.c,v 1.12.2.4 2009-03-20 19:42:26 mt269 Exp $ */

/* Handling of specific substitution models.  This needs reworking:
   was originally set up for small number of models but has become
   unwieldy... */

#include <subst_mods.h>
#include <tree_model.h>
#include <stringsplus.h>
#include <string.h>
#include <ctype.h>
#include <misc.h>

/* internal functions (model-specific) */
void tm_set_JC69_matrix(TreeModel *mod);
void tm_set_K80_matrix(TreeModel *mod, double kappa);
void tm_set_HKY_matrix(TreeModel *mod, double kappa, int kappa_idx);
void tm_set_HKYG_matrix(TreeModel *mod, Vector *params, int start_idx);
void tm_set_REV_matrix(TreeModel *mod, Vector *params, int start_idx);
void tm_set_SSREV_matrix(TreeModel *mod, Vector *params, int start_idx);
void tm_set_UNREST_matrix(TreeModel *mod, Vector *params, int start_idx);
void tm_set_R2_matrix(TreeModel *mod, Vector *params, int start_idx);
void tm_set_U2_matrix(TreeModel *mod, Vector *params, int start_idx);
void tm_set_R2S_matrix(TreeModel *mod, Vector *params, int start_idx);
void tm_set_U2S_matrix(TreeModel *mod, Vector *params, 
                       int start_idx);
void tm_set_R3_matrix(TreeModel *mod, Vector *params, int start_idx);
void tm_set_R3S_matrix(TreeModel *mod, Vector *params, int start_idx);
void tm_set_U3_matrix(TreeModel *mod, Vector *params, int start_idx);
void tm_set_U3S_matrix(TreeModel *mod, Vector *params, int start_idx);
void tm_set_GC_matrix(TreeModel *mod, double kappa, int kappa_idx, double alpha);
void tm_set_REV_GC_matrix(TreeModel *mod, Vector *params, int start_idx);
void tm_init_mat_REV(TreeModel *mod, Vector *params, int nbranches, 
                     double kappa);
void tm_init_mat_SSREV(TreeModel *mod, Vector *params, int nbranches,
		       double kappa);
void tm_init_mat_REV_GC(TreeModel *mod, Vector *params, int nbranches, double kappa);
void tm_init_mat_UNREST(TreeModel *mod, Vector *params, int nbranches, 
                        double kappa);
void tm_init_mat_R2(TreeModel *mod, Vector *params, 
                    int start_idx, double kappa);
void tm_init_mat_U2(TreeModel *mod, Vector *params, 
                    int start_idx, double kappa);
void tm_init_mat_R2S(TreeModel *mod, Vector *params, 
                     int start_idx, double kappa);
void tm_init_mat_U2S(TreeModel *mod, Vector *params, 
                     int start_idx, double kappa);
void tm_init_mat_R3(TreeModel *mod, Vector *params, 
                    int start_idx, double kappa);
void tm_init_mat_R3S(TreeModel *mod, Vector *params, 
                     int start_idx, double kappa);
void tm_init_mat_U3(TreeModel *mod, Vector *params, 
                    int start_idx, double kappa);
void tm_init_mat_U3S(TreeModel *mod, Vector *params, 
                     int start_idx, double kappa);
void tm_init_mat_from_model_REV(TreeModel *mod, Vector *params, 
                                int start_idx);
void tm_init_mat_from_model_SSREV(TreeModel *mod, Vector *params, 
                                int start_idx);
void tm_init_mat_from_model_REV_GC(TreeModel *mod, Vector *params, 
				    int start_idx);
void tm_init_mat_from_model_UNREST(TreeModel *mod, Vector *params, 
                                   int start_idx);
void tm_init_mat_from_model_R2(TreeModel *mod, Vector *params, 
                               int start_idx);
void tm_init_mat_from_model_U2(TreeModel *mod, Vector *params, 
                               int start_idx);
void tm_init_mat_from_model_R2S(TreeModel *mod, Vector *params, 
                                int start_idx);
void tm_init_mat_from_model_U2S(TreeModel *mod, Vector *params, 
                                int start_idx);
void tm_init_mat_from_model_R3(TreeModel *mod, Vector *params, 
                               int start_idx);
void tm_init_mat_from_model_R3S(TreeModel *mod, Vector *params, 
                                int start_idx);
void tm_init_mat_from_model_U3(TreeModel *mod, Vector *params, 
                               int start_idx);
void tm_init_mat_from_model_U3S(TreeModel *mod, Vector *params, 
                                int start_idx);
int tm_flag_subst_param_pos(TreeModel *mod, int *flag, 
			    String *param_name);


/* Return the substitution model (enum val) corresponding to the
   specified string */
subst_mod_type tm_get_subst_mod_type(const char *str) {
  String *subst_mod_str = str_new_charstr(str);
  subst_mod_type retval = UNDEF_MOD;
  if (str_equals_nocase_charstr(subst_mod_str, "JC69"))
    retval = JC69;
  else if (str_equals_nocase_charstr(subst_mod_str, "K80"))
    retval = K80;
  else if (str_equals_nocase_charstr(subst_mod_str, "F81"))
    retval = F81;
  else if (str_equals_nocase_charstr(subst_mod_str, "HKY85"))
    retval = HKY85;
  else if (str_equals_nocase_charstr(subst_mod_str, "HKY85+Gap"))
    retval = HKY85G;
  else if (str_equals_nocase_charstr(subst_mod_str, "REV"))
    retval = REV;
  else if (str_equals_nocase_charstr(subst_mod_str, "SSREV"))
    retval = SSREV;
  else if (str_equals_nocase_charstr(subst_mod_str, "REV+GC"))
    retval = REV_GC;
  else if (str_equals_nocase_charstr(subst_mod_str, "UNREST"))
    retval = UNREST;
  else if (str_equals_nocase_charstr(subst_mod_str, "R2"))
    retval = R2;
  else if (str_equals_nocase_charstr(subst_mod_str, "U2"))
    retval = U2;
  else if (str_equals_nocase_charstr(subst_mod_str, "R2S"))
    retval = R2S;
  else if (str_equals_nocase_charstr(subst_mod_str, "U2S"))
    retval = U2S;
  else if (str_equals_nocase_charstr(subst_mod_str, "R3"))
    retval = R3;
  else if (str_equals_nocase_charstr(subst_mod_str, "R3S"))
    retval = R3S;
  else if (str_equals_nocase_charstr(subst_mod_str, "U3"))
    retval = U3;
  else if (str_equals_nocase_charstr(subst_mod_str, "U3S"))
    retval = U3S;
  else if (str_equals_nocase_charstr(subst_mod_str, "GC"))
    retval = GC;
  else if (str_equals_nocase_charstr(subst_mod_str, "HB"))
    retval = HB;
  str_free(subst_mod_str);
  return retval;
}

/* Return a string description for the specified subst_mod_type  */
char *tm_get_subst_mod_string(subst_mod_type type) {
  switch(type) {
  case JC69:
    return "JC69";
  case K80:
    return "K80";
  case F81:
    return "F81";
  case HKY85:
    return "HKY85";
  case HKY85G:
    return "HKY85+Gap";
  case REV:
    return "REV";
  case SSREV:
    return "SSREV";
  case REV_GC:
    return "REV+GC";
  case UNREST:
    return "UNREST";
  case R2:
    return "R2";
  case U2:
    return "U2";
  case R2S:
    return "R2S";
  case U2S:
    return "U2S";
  case R3:
    return "R3";
  case R3S:
    return "R3S";
  case U3:
    return "U3";
  case U3S:
    return "U3S";
  case GC:
    return "GC";
  case HB:
    return "HB";
  default:
    return "(unknown model)";
  }
}


/* number of rate matrix params (not counting eq. freqs) */
int tm_get_nratematparams(TreeModel *mod) {

  //  if (mod->estimate_ratemat == FALSE) return 0;

  switch (mod->subst_mod) {
  case JC69:
  case F81:
    return 0;
  case K80:
  case HKY85:
    return 1;
  case HKY85G:
    return 2;
  case REV:
    return (mod->rate_matrix->size * mod->rate_matrix->size 
            - mod->rate_matrix->size) / 2; 
    /* allows use with alternative alphabets, e.g., {ACGT-} */
  case SSREV:
    return (mod->rate_matrix->size * mod->rate_matrix->size
	    - mod->rate_matrix->size)/2 - 
      mod->rate_matrix->size/2 -
      2*((mod->rate_matrix->size/2*2)!=mod->rate_matrix->size);
    /* subtract mod->rate_matrix->size/2 assuming that all but possibly one
       of the characters (ie, gap) has a complement.  The last line adds 2 
       if there is an odd number of bases*/
  case REV_GC:
    return (mod->rate_matrix->size * mod->rate_matrix->size
	    - mod->rate_matrix->size)/2 + 1;
  case UNREST:
    return (mod->rate_matrix->size * mod->rate_matrix->size 
            - mod->rate_matrix->size);     
    /* allows use with alternative alphabets, e.g., {ACGT-} */
  case R2:
    return 48; 
  case U2:
    return 96;
  case R2S:
    return 24; 
  case U2S:
    return 48;
  case R3:
    return 288;
  case R3S:
    return 148;
    /* every parameter covers 4 elements of the rate matrix, except
       for 8 parameters, which cover only two (the ones corresponding
       to subst. of the form X Y X_compl <-> X Y_compl X_compl); thus,
       the number of params required is (576 - 8*2)/4 + 8 = 148 */
  case U3:
    return 576;
  case U3S:
    return 288;
  case GC:
    return 2;
  case HB:
    return 0;
  default:
    die("ERROR: tm_get_nrateparams unknown substitution model\n");
  }
  return -1;
}

int tm_is_reversible(int subst_mod) {
  return !(subst_mod == UNREST || subst_mod == U2 || subst_mod == U2S || 
           subst_mod == U3 || subst_mod == U3S || subst_mod == GC);
}

int tm_order(int subst_mod) {
  if (subst_mod == R2 || subst_mod == U2 || 
      subst_mod == R2S || subst_mod == U2S)
    return 1; 
  else if (subst_mod == R3 || subst_mod == R3S || 
           subst_mod == U3 || subst_mod == U3S) 
    return 2;

  return 0;
}

/* Set rate matrix according to elements of parameter vector (i is
   starting index) */
void tm_set_rate_matrix(TreeModel *mod, Vector *params, int i) {
  switch(mod->subst_mod) {
  case JC69:
    tm_set_JC69_matrix(mod);
    break;
  case K80:
    tm_set_K80_matrix(mod, vec_get(params, i));
    break;
  case F81:
    tm_set_HKY_matrix(mod, 1, -1); 
    break;
  case HKY85:
    tm_set_HKY_matrix(mod, vec_get(params, i), i);
    break;
  case HKY85G:
    tm_set_HKYG_matrix(mod, params, i);
    break;
  case REV:
    tm_set_REV_matrix(mod, params, i);
    break;
  case SSREV:
    tm_set_SSREV_matrix(mod, params, i);
    break;
  case REV_GC:
    tm_set_REV_GC_matrix(mod, params, i);
    break;
  case UNREST:
    tm_set_UNREST_matrix(mod, params, i);
    break;
  case R2:
    tm_set_R2_matrix(mod, params, i);
    break;
  case U2:
    tm_set_U2_matrix(mod, params, i);
    break;
  case R2S:
    tm_set_R2S_matrix(mod, params, i);
    break;
  case U2S:
    tm_set_U2S_matrix(mod, params, i);
    break;
  case R3:
    tm_set_R3_matrix(mod, params, i);
    break;
  case R3S:
    tm_set_R3S_matrix(mod, params, i);
    break;
  case U3:
    tm_set_U3_matrix(mod, params, i);
    break;
  case U3S:
    tm_set_U3S_matrix(mod, params, i);
    break;
  case GC:
    tm_set_GC_matrix(mod, vec_get(params, i), i, vec_get(params, i+1));
    break;
  default:
    die("ERROR tm_set_rate_matrix: unknown substitution model\n");
  }

  if (mod->scale_during_opt && mod->subst_mod!=JC69 && mod->subst_mod != F81)
    tm_scale_rate_matrix(mod);

  /* NOTE: used to scale rate matrices here, but have removed scaling
     for more well-behaved derivatives.  */
}

/* initialize rate-matrix parameters in parameter vector, using an
   HKY-like strategy (with 'kappa' defining transition/transversion
   bias).  Starting index is 'params_idx'. */
void tm_rate_params_init(TreeModel *mod, Vector *params, 
                         int params_idx, double kappa) {
  switch(mod->subst_mod) {
  case JC69:
  case F81:
    break;                      /* do nothing */
  case K80:
  case HKY85:
    vec_set(params, params_idx, kappa);
    break;      
  case HKY85G:
    vec_set(params, params_idx, kappa);
    /* sigma often ends up near kappa, so this provides a reasonable
       initial guess. */
    vec_set(params, params_idx+1, kappa);
    break;    
  case GC:
    vec_set(params, params_idx, kappa);
    vec_set(params, params_idx+1, 1);
    break;    
  case REV:
    tm_init_mat_REV(mod, params, params_idx, kappa);
    break;      
  case SSREV:
    tm_init_mat_SSREV(mod, params, params_idx, kappa);
    break; 
  case REV_GC:
    tm_init_mat_REV_GC(mod, params, params_idx, kappa);
    break;
  case UNREST:
    tm_init_mat_UNREST(mod, params, params_idx, kappa);
    break;      
  case R2:
    tm_init_mat_R2(mod, params, params_idx, kappa);
    break;      
  case U2:
    tm_init_mat_U2(mod, params, params_idx, kappa);
    break;      
  case R2S:
    tm_init_mat_R2S(mod, params, params_idx, kappa);
    break;      
  case U2S:
    tm_init_mat_U2S(mod, params, params_idx, kappa);
    break;      
  case R3:
    tm_init_mat_R3(mod, params, params_idx, kappa);
    break;
  case R3S:
    tm_init_mat_R3S(mod, params, params_idx, kappa);
    break;
  case U3:
    tm_init_mat_U3(mod, params, params_idx, kappa);
    break;
  case U3S:
    tm_init_mat_U3S(mod, params, params_idx, kappa);
    break;
  default:
    die("ERROR tm_rate_params_init: unknown substitution model\n");
  }
}

/* initialize rate-matrix parameters in parameter vector, based on an
   existing model.  Starting index is 'params_idx'. */
void tm_rate_params_init_from_model(TreeModel *mod, Vector *params, 
                                    int params_idx) {
  double kappa, alpha;

  switch(mod->subst_mod) {
  case JC69:
  case F81:
    break;                      /* do nothing */
  case K80:
  case HKY85:
    /* infer kappa from rate matrix */
    kappa = 
      mm_get(mod->rate_matrix, 
             mod->rate_matrix->inv_states['G'], 
             mod->rate_matrix->inv_states['A']) /
      mm_get(mod->rate_matrix, 
             mod->rate_matrix->inv_states['C'], 
             mod->rate_matrix->inv_states['A']);

    vec_set(params, params_idx, kappa);
    break;      
  case HKY85G:
    /* infer kappa from rate matrix */
    kappa = 
      mm_get(mod->rate_matrix, 
             mod->rate_matrix->inv_states['G'], 
             mod->rate_matrix->inv_states['A']) /
      mm_get(mod->rate_matrix, 
             mod->rate_matrix->inv_states['C'], 
             mod->rate_matrix->inv_states['A']);
    vec_set(params, params_idx, kappa);
    vec_set(params, params_idx+1, kappa);
    break;    
  case GC:
    /* infer kappa from rate matrix */
    kappa = 
      mm_get(mod->rate_matrix, 
             mod->rate_matrix->inv_states['G'], 
             mod->rate_matrix->inv_states['A']) /
      mm_get(mod->rate_matrix, 
             mod->rate_matrix->inv_states['C'], 
             mod->rate_matrix->inv_states['A']);
    vec_set(params, params_idx, kappa);
    alpha = 
      mm_get(mod->rate_matrix,
	     mod->rate_matrix->inv_states['A'],
	     mod->rate_matrix->inv_states['C']) /
      mm_get(mod->rate_matrix,
	     mod->rate_matrix->inv_states['G'],
	     mod->rate_matrix->inv_states['C']);
    vec_set(params, params_idx+1, alpha); /* just use a 1 for gamma */
    break;
  case REV_GC:
    vec_set(params, params_idx++, 1.0);
  case REV:
    tm_init_mat_from_model_REV(mod, params, params_idx);
    break;      
  case SSREV:
    tm_init_mat_from_model_SSREV(mod, params, params_idx);
    break;  
  case UNREST:
    tm_init_mat_from_model_UNREST(mod, params, params_idx);
    break;      
  case R2:
    tm_init_mat_from_model_R2(mod, params, params_idx);
    break;      
  case U2:
    tm_init_mat_from_model_U2(mod, params, params_idx);
    break;      
  case R2S:
    tm_init_mat_from_model_R2S(mod, params, params_idx);
    break;      
  case U2S:
    tm_init_mat_from_model_U2S(mod, params, params_idx);
    break;      
  case R3:
    tm_init_mat_from_model_R3(mod, params, params_idx);
    break;      
  case U3:
    tm_init_mat_from_model_U3(mod, params, params_idx);
    break;      
  case U3S:
    tm_init_mat_from_model_U3S(mod, params, params_idx);
    break;      
  default:
    die("ERROR tm_rate_params_init_from_model: unknown substitution model\n");
  }
}

/* functions to compute subst. probability matrices directly (without
   diagonalization of a rate matrix), for some simple models */

void tm_set_probs_JC69(TreeModel *mod, MarkovMatrix *P, double t) {
  int i, j;
  double scale = mod->rate_matrix->size * 1.0/(mod->rate_matrix->size - 1);
  if (t < 0) die("ERROR tm_set_probs_JC69 t should be >=0 but is %f\n", t);
  for (i = 0; i < mod->rate_matrix->size; i++) {
    for (j = 0; j < mod->rate_matrix->size; j++) {
      if (i == j)
        mm_set(P, i, j, 1.0/mod->rate_matrix->size + 
               (1 - 1.0/mod->rate_matrix->size) * exp(-t * scale));
      else
        mm_set(P, i, j, 1.0/mod->rate_matrix->size - 
               1.0/mod->rate_matrix->size * exp(-t * scale));
    }
  }
}

void tm_set_probs_F81(Vector *backgd_freqs, MarkovMatrix *P, double scale, 
                      double t) {
  int i, j;
  for (i = 0; i < backgd_freqs->size; i++) {
    for (j = 0; j < backgd_freqs->size; j++) {
      if (i == j)
        mm_set(P, i, j, exp(-t * scale) + 
               vec_get(backgd_freqs, j) * (1 - exp(-t * scale)));
      else
        mm_set(P, i, j, vec_get(backgd_freqs, j) * 
               (1 - exp(-t * scale)));

    }
  }
}

/* set matrix such that element (i,j) has value pi_j, as for an
   infinitely long branch */
void tm_set_probs_independent(TreeModel *mod, MarkovMatrix *P) {
  int i, j;
  for (i = 0; i < mod->rate_matrix->size; i++) 
    for (j = 0; j < mod->rate_matrix->size; j++) 
      mm_set(P, i, j, vec_get(mod->backgd_freqs, j));
}


/***************************************************************************/
/* Functions to map from parameter vectors to rate matrices -- one         */
/* function per substitution model                                         */
/***************************************************************************/

void tm_set_JC69_matrix(TreeModel *mod) {
  int i, j;
  for (i = 0; i < mod->rate_matrix->size; i++) {
    for (j = 0; j < mod->rate_matrix->size; j++) {
      if (j == i) continue;
      mm_set(mod->rate_matrix, i, j, 1.0/mod->rate_matrix->size);
    }
    mm_set(mod->rate_matrix, i, i, -1);
  }
}

void tm_set_K80_matrix(TreeModel *mod, double kappa) {
  int i, j;
  for (i = 0; i < mod->rate_matrix->size; i++) {
    double rowsum = 0;
    for (j = 0; j < mod->rate_matrix->size; j++) {
      double val;
      if (j == i) continue;
      val = 1;
      if (is_transition(mod->rate_matrix->states[i], 
                        mod->rate_matrix->states[j]))
        val *= kappa;
      mm_set(mod->rate_matrix, i, j, val);
      rowsum += val;
    }
    mm_set(mod->rate_matrix, i, i, -1 * rowsum);
  }
}

void tm_set_HKY_matrix(TreeModel *mod, double kappa, int kappa_idx) {
  int i, j;
  int setup_mapping = 
    (kappa_idx >= 0 && lst_size(mod->rate_matrix_param_row[kappa_idx]) == 0);
  for (i = 0; i < mod->rate_matrix->size; i++) {
    double rowsum = 0;
    for (j = 0; j < mod->rate_matrix->size; j++) {
      double val;
      if (j == i) continue;
      val = vec_get(mod->backgd_freqs, j); 
      if (is_transition(mod->rate_matrix->states[i], 
                           mod->rate_matrix->states[j])) {
        val *= kappa;

        if (setup_mapping) {
          lst_push_int(mod->rate_matrix_param_row[kappa_idx], i);
          lst_push_int(mod->rate_matrix_param_col[kappa_idx], j);
        }
      }

      mm_set(mod->rate_matrix, i, j, val);
      rowsum += val;
    }
    mm_set(mod->rate_matrix, i, i, -1 * rowsum);
  }
}

void tm_set_HKYG_matrix(TreeModel *mod, Vector *params, int start_idx ) {
  int i, j;
  int setup_mapping = lst_size(mod->rate_matrix_param_row[start_idx]) == 0;
  int kappa_idx = start_idx;
  int sigma_idx = start_idx + 1;
  double kappa = vec_get(params, kappa_idx);
  double sigma = vec_get(params, sigma_idx);
  for (i = 0; i < mod->rate_matrix->size; i++) {
    double rowsum = 0;
    for (j = 0; j < mod->rate_matrix->size; j++) {
      double val;
      if (j == i) continue;
      val = vec_get(mod->backgd_freqs, j); 
      if (is_transition(mod->rate_matrix->states[i], 
                        mod->rate_matrix->states[j])) {
        val *= kappa;
        if (setup_mapping) {
          lst_push_int(mod->rate_matrix_param_row[kappa_idx], i);
          lst_push_int(mod->rate_matrix_param_col[kappa_idx], j);
        }
      }
      if (is_indel(mod->rate_matrix->states[i], 
                   mod->rate_matrix->states[j])) {
        val *= sigma;
        if (setup_mapping) {
          lst_push_int(mod->rate_matrix_param_row[sigma_idx], i);
          lst_push_int(mod->rate_matrix_param_col[sigma_idx], j);
        }
      }
      mm_set(mod->rate_matrix, i, j, val);
      rowsum += val;
    }
    mm_set(mod->rate_matrix, i, i, -1 * rowsum);
  }
}

void tm_set_GC_matrix(TreeModel *mod, double kappa, int kappa_idx, double alpha) {
  int i, j;
  char c;
  double sum = 0.0;
  int setup_mapping = 
    (kappa_idx >= 0 && lst_size(mod->rate_matrix_param_row[kappa_idx]) == 0);

  //first scale eq frequencies
  for (i=0; i<mod->rate_matrix->size; i++) {
    c = toupper(mod->rate_matrix->states[i]);
    if (c=='G' || c=='C')
      vec_set(mod->backgd_freqs, i, vec_get(mod->backgd_freqs, i)*alpha);
    sum += vec_get(mod->backgd_freqs, i);
  }
  for (i=0; i<mod->rate_matrix->size; i++)
    vec_set(mod->backgd_freqs, i, vec_get(mod->backgd_freqs, i)/sum);

  for (i = 0; i < mod->rate_matrix->size; i++) {
    double rowsum = 0;
    for (j = 0; j < mod->rate_matrix->size; j++) {
      double val;
      if (j == i) continue;
      val = vec_get(mod->backgd_freqs, j); 
      if (is_transition(mod->rate_matrix->states[i], 
                           mod->rate_matrix->states[j])) {
        val *= kappa;

        if (setup_mapping) {
          lst_push_int(mod->rate_matrix_param_row[kappa_idx], i);
          lst_push_int(mod->rate_matrix_param_col[kappa_idx], j);
        }
      }

      /*      if ((mod->rate_matrix->states[j] == 'G' || mod->rate_matrix->states[j] == 'C') &&
	      (mod->rate_matrix->states[i] != 'G' && mod->rate_matrix->states[i] != 'C')) {*/
      if ((mod->rate_matrix->states[j] == 'G' || mod->rate_matrix->states[j]=='C') &&
	  (mod->rate_matrix->states[i] == 'G' || mod->rate_matrix->states[i] == 'C')) {

	//        val *= alpha;
	val *= sum/alpha;

        if (setup_mapping) {
          lst_push_int(mod->rate_matrix_param_row[kappa_idx+1], i);
          lst_push_int(mod->rate_matrix_param_col[kappa_idx+1], j);
        } /* Check: is this right? */
      }

      mm_set(mod->rate_matrix, i, j, val);
      rowsum += val;
    }
    mm_set(mod->rate_matrix, i, i, -1 * rowsum);
  }
}

void tm_set_REV_matrix(TreeModel *mod, Vector *params, int start_idx) {
  int i, j;
  int setup_mapping = (lst_size(mod->rate_matrix_param_row[start_idx]) == 0);
  for (i = 0; i < mod->rate_matrix->size; i++) {
    double rowsum = 0;
    for (j = i+1; j < mod->rate_matrix->size; j++) {
      double val;
      val = vec_get(params, start_idx);
      mm_set(mod->rate_matrix, i, j, 
             val * vec_get(mod->backgd_freqs, j));
      mm_set(mod->rate_matrix, j, i, 
             val * vec_get(mod->backgd_freqs, i));
      rowsum += (val * vec_get(mod->backgd_freqs, j));

      /* experimental */
      if (setup_mapping) {
        lst_push_int(mod->rate_matrix_param_row[start_idx], i);
        lst_push_int(mod->rate_matrix_param_col[start_idx], j);
        lst_push_int(mod->rate_matrix_param_row[start_idx], j);
        lst_push_int(mod->rate_matrix_param_col[start_idx], i);
      }

      start_idx++;
    }
    for (j = 0; j < i; j++)
      rowsum += mm_get(mod->rate_matrix, i, j);
    mm_set(mod->rate_matrix, i, i, -1 * rowsum);
  }
}

void tm_set_SSREV_matrix(TreeModel *mod, Vector *params, int start_idx) {
  int i, j, compi, compj;
  int setup_mapping = (lst_size(mod->rate_matrix_param_row[start_idx]) == 0);
  int count=0;  //testing

  for (i = 0; i < mod->rate_matrix->size; i++) {
    double rowsum = 0;
    compi = mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->states[i])];
    for (j = i+1; j < mod->rate_matrix->size; j++) {
      double val;
      compj = mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->states[j])];

      //if this is true we have already computed this position
      if ((compi < compj && compi < i) ||
	  (compj < compi && compj < i)) continue;

      val = vec_get(params, start_idx);
      mm_set(mod->rate_matrix, i, j, 
             val * vec_get(mod->backgd_freqs, j));
      mm_set(mod->rate_matrix, j, i, 
             val * vec_get(mod->backgd_freqs, i));
      /* experimental */
      if (setup_mapping) {
        lst_push_int(mod->rate_matrix_param_row[start_idx], i);
        lst_push_int(mod->rate_matrix_param_col[start_idx], j);
        lst_push_int(mod->rate_matrix_param_row[start_idx], j);
        lst_push_int(mod->rate_matrix_param_col[start_idx], i);
      }
      if (compi!=j) {
	mm_set(mod->rate_matrix, compi, compj, 
	       val * vec_get(mod->backgd_freqs, compj));
	mm_set(mod->rate_matrix, compj, compi,
	       val * vec_get(mod->backgd_freqs, compi));
	/* experimental */
	if (setup_mapping) {
	  lst_push_int(mod->rate_matrix_param_row[start_idx], compi);
	  lst_push_int(mod->rate_matrix_param_col[start_idx], compj);
	  lst_push_int(mod->rate_matrix_param_row[start_idx], compj);
	  lst_push_int(mod->rate_matrix_param_col[start_idx], compi);
	}
      }
      
      start_idx++;
      count++;
    }
    rowsum = 0.0;
    for (j=0; j<mod->rate_matrix->size; j++)
      if (j!=i) rowsum += mm_get(mod->rate_matrix, i, j);
    mm_set(mod->rate_matrix, i, i, -1 * rowsum);
  }
  //testing
  if (count != tm_get_nratematparams(mod))
    die("ERROR tm_set_SSREV_matrix count (%i) != tm_get_nratematparams(mod) (%i)\n",
	count, tm_get_nratematparams(mod));
}


void tm_set_REV_GC_matrix(TreeModel *mod, Vector *params, int start_idx) {
  int i, j, gamma_idx = start_idx;
  int setup_mapping = (lst_size(mod->rate_matrix_param_row[start_idx])==0);
  double gamma, sum=0.0;
  char c1, c2;

  gamma = vec_get(params, start_idx++);

  //first scale the eq frequencies
  for (i=0; i<mod->rate_matrix->size; i++) {
    c1 = toupper(mod->rate_matrix->states[i]);
    if (c1 == 'C' || c1 == 'G')
      vec_set(mod->backgd_freqs, i, vec_get(mod->backgd_freqs, i)*gamma);
    sum += vec_get(mod->backgd_freqs, i);
  }
  for (i=0; i<mod->rate_matrix->size; i++)
    vec_set(mod->backgd_freqs, i, vec_get(mod->backgd_freqs, i)/sum);

  for (i = 0; i < mod->rate_matrix->size; i++) {
    double rowsum = 0;
    c1 = toupper(mod->rate_matrix->states[i]);
    for (j = i+1; j<mod->rate_matrix->size; j++) {
      double val;
      //set val to alpha{i,j}
      val = vec_get(params, start_idx);
      c2 = toupper(mod->rate_matrix->states[j]);
      if ((c1=='C' && c2=='G') ||
	  (c1=='G' && c2=='C')) {
	val *= sum/gamma;
	if (setup_mapping) {
	  lst_push_int(mod->rate_matrix_param_row[gamma_idx], i);
	  lst_push_int(mod->rate_matrix_param_col[gamma_idx], j);
	  lst_push_int(mod->rate_matrix_param_row[gamma_idx], j);
	  lst_push_int(mod->rate_matrix_param_col[gamma_idx], i);
	}
      }
      mm_set(mod->rate_matrix, i, j, val*vec_get(mod->backgd_freqs, j));
      mm_set(mod->rate_matrix, j, i, val*vec_get(mod->backgd_freqs, i));
      rowsum += (val * vec_get(mod->backgd_freqs, j));
      if (setup_mapping) {
        lst_push_int(mod->rate_matrix_param_row[start_idx], i);
        lst_push_int(mod->rate_matrix_param_col[start_idx], j);
        lst_push_int(mod->rate_matrix_param_row[start_idx], j);
        lst_push_int(mod->rate_matrix_param_col[start_idx], i);
      }
      start_idx++;
    }
    for (j=0; j<i; j++)
      rowsum += mm_get(mod->rate_matrix, i, j);
    mm_set(mod->rate_matrix, i, i, -1.0*rowsum);
  }
}


void tm_set_UNREST_matrix(TreeModel *mod, Vector *params, int start_idx) {
  int i, j;
  int setup_mapping = (lst_size(mod->rate_matrix_param_row[start_idx]) == 0);
  mat_zero(mod->rate_matrix->matrix);
  for (i = 0; i < mod->rate_matrix->size; i++) {
    double rowsum = 0;
    for (j = 0; j < mod->rate_matrix->size; j++) {
      double val;
      if (i == j) continue;
      mm_set(mod->rate_matrix, i, j, 
             val = vec_get(params, start_idx));
      rowsum += val;

      /* experimental */
      if (setup_mapping) {
        lst_push_int(mod->rate_matrix_param_row[start_idx], i);
        lst_push_int(mod->rate_matrix_param_col[start_idx], j);
      }

      start_idx++;
    }
    mm_set(mod->rate_matrix, i, i, -1 * rowsum);
  }
}

/* void tm_set_HKY2_matrix(TreeModel *mod, double kappa) { */
/*   int i, j; */
/*   int alph_size = strlen(mod->rate_matrix->states); */
/*   for (i = 0; i < mod->rate_matrix->size; i++) { */
/*     double rowsum = 0, val; */
/*     int b1_i, b2_i, b1_j, b2_j; */
/*     b1_i = i / alph_size; */
/*     b2_i = i % alph_size; */
/*     for (j = 0; j < mod->rate_matrix->size; j++) { */
/*       if (i == j) continue; */
/*       val = vec_get(mod->backgd_freqs, j); */
/*       b1_j = j / alph_size; */
/*       b2_j = j % alph_size; */
/*       if (b1_i == b1_j) { */
/*         if (is_transition(mod->rate_matrix->states[b2_i],  */
/*                              mod->rate_matrix->states[b2_j]))  */
/*           val *= kappa; */
/*       } */
/*       else if (b2_i == b2_j) { */
/*         if (is_transition(mod->rate_matrix->states[b1_i],  */
/*                              mod->rate_matrix->states[b1_j]))  */
/*           val *= kappa; */
/*       } */
/*       else  */
/*         val = 0; */

/*       mm_set(mod->rate_matrix, i, j, val); */
/*       rowsum += val; */
/*     } */
/*     mm_set(mod->rate_matrix, i, i, -1 * rowsum); */
/*   } */
/* } */

void tm_set_R2_matrix(TreeModel *mod, Vector *params, int start_idx) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  int setup_mapping = (lst_size(mod->rate_matrix_param_row[start_idx]) == 0);
  mat_zero(mod->rate_matrix->matrix);

  for (i = 0; i < mod->rate_matrix->size; i++) {
    double rowsum = 0, val;
    int b1_i, b2_i, b1_j, b2_j;
    b1_i = i / alph_size;
    b2_i = i % alph_size;
    for (j = i+1; j < mod->rate_matrix->size; j++) {
      b1_j = j / alph_size;
      b2_j = j % alph_size;
      if (b1_i != b1_j && b2_i != b2_j)
        continue;

      val = vec_get(params, start_idx);
      mm_set(mod->rate_matrix, i, j, 
             val * vec_get(mod->backgd_freqs, j));
      mm_set(mod->rate_matrix, j, i, 
             val * vec_get(mod->backgd_freqs, i));
      rowsum += (val * vec_get(mod->backgd_freqs, j));

      /* experimental */
      if (setup_mapping) {
        lst_push_int(mod->rate_matrix_param_row[start_idx], i);
        lst_push_int(mod->rate_matrix_param_col[start_idx], j);
        lst_push_int(mod->rate_matrix_param_row[start_idx], j);
        lst_push_int(mod->rate_matrix_param_col[start_idx], i);
      }

      start_idx++;
    }        
    for (j = 0; j < i; j++)
      rowsum += mm_get(mod->rate_matrix, i, j);
    mm_set(mod->rate_matrix, i, i, -1 * rowsum);
  }
  if (start_idx != params->size)
    die("ERROR tm_set_R2_matrix: start_idx (%i) != params->size (%i)\n",
	start_idx, params->size);
}

void tm_set_R2S_matrix(TreeModel *mod, Vector *params, 
                            int start_idx) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  int nstates = mod->rate_matrix->size;
  int done_row[nstates];
  int setup_mapping = (lst_size(mod->rate_matrix_param_row[start_idx]) == 0);
  double rowsum;
  mat_zero(mod->rate_matrix->matrix);
  
  for (i = 0; i < nstates; i++) done_row[i] = 0;

  for (i = 0; i < nstates; i++) {
    double val;
    int b1_i, b2_i, b1_j, b2_j, b1_comp_i, b2_comp_i, i_comp, b1_comp_j, 
      b2_comp_j, j_comp;

    if (done_row[i]) continue;

    b1_i = i / alph_size;       /* idx of first char of "from" dinuc */
    b2_i = i % alph_size;       /* idx of second char */
    b1_comp_i =                 /* idx of compl of first char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b1_i])];
    b2_comp_i =                 /* idx of compl of second char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b2_i])];
    i_comp = b2_comp_i * alph_size + b1_comp_i;
                                /* state idx of compl dinuc */

    for (j = i+1; j < nstates; j++) {

      b1_j = j / alph_size;     /* idx of first char of "to" dinuc */
      b2_j = j % alph_size;     /* idx of second char */
      b1_comp_j =               /* idx of compl of first char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b1_j])];
      b2_comp_j =               /* idx of compl of second char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b2_j])];
      j_comp = b2_comp_j * alph_size + b1_comp_j;
                                /* state idx of compl dinuc */

      if (done_row[j_comp] || (i == i_comp && j_comp < j)) continue;

      if (b1_i != b1_j && b2_i != b2_j)
        continue;

      val = vec_get(params, start_idx);

      /* set rates for i,j and j,i */
      mm_set(mod->rate_matrix, i, j, 
             val * vec_get(mod->backgd_freqs, j));
      mm_set(mod->rate_matrix, j, i, 
             val * vec_get(mod->backgd_freqs, i));

      /* set rates for their complements */
      mm_set(mod->rate_matrix, i_comp, j_comp, 
             val * vec_get(mod->backgd_freqs, j_comp));
      mm_set(mod->rate_matrix, j_comp, i_comp, 
             val * vec_get(mod->backgd_freqs, i_comp));

      /* experimental */
      if (setup_mapping) {
        lst_push_int(mod->rate_matrix_param_row[start_idx], i);
        lst_push_int(mod->rate_matrix_param_col[start_idx], j);
        lst_push_int(mod->rate_matrix_param_row[start_idx], j);
        lst_push_int(mod->rate_matrix_param_col[start_idx], i);
        lst_push_int(mod->rate_matrix_param_row[start_idx], i_comp);
        lst_push_int(mod->rate_matrix_param_col[start_idx], j_comp);
        lst_push_int(mod->rate_matrix_param_row[start_idx], j_comp);
        lst_push_int(mod->rate_matrix_param_col[start_idx], i_comp);
      }

      start_idx++;
    }        
    done_row[i] = done_row[i_comp] = 1;
  }

  for (i = 0; i < nstates; i++) {
    rowsum = 0;    
    for (j = 0; j < nstates; j++)
      if (j != i) rowsum += mm_get(mod->rate_matrix, i, j);
    mm_set(mod->rate_matrix, i, i, -1 * rowsum);
  }
  
  if (start_idx != params->size)
    die("ERROR: tm_set_R2S_matrix start_idx (%i) != params->size (%i)\n",
	start_idx, params->size);
}

void tm_set_U2_matrix(TreeModel *mod, Vector *params, int start_idx) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  int setup_mapping = (lst_size(mod->rate_matrix_param_row[start_idx]) == 0);
  mat_zero(mod->rate_matrix->matrix);
  for (i = 0; i < mod->rate_matrix->size; i++) {
    double rowsum = 0, val;
    int b1_i, b2_i, b1_j, b2_j;
    b1_i = i / alph_size;
    b2_i = i % alph_size;
    for (j = 0; j < mod->rate_matrix->size; j++) {
      if (i == j) continue;
      b1_j = j / alph_size;
      b2_j = j % alph_size;
      if (b1_i != b1_j && b2_i != b2_j)
        continue;

      val = vec_get(params, start_idx);
      mm_set(mod->rate_matrix, i, j, val);
      rowsum += val;

      /* experimental */
      if (setup_mapping) {
        lst_push_int(mod->rate_matrix_param_row[start_idx], i);
        lst_push_int(mod->rate_matrix_param_col[start_idx], j);
      }

      start_idx++;
    }
    mm_set(mod->rate_matrix, i, i, -1 * rowsum);
  }
  if (start_idx != params->size)
    die("ERROR tm_set_U2_matrix start_idx (%i) != params->size (%i)\n",
	start_idx, params->size);
}

void tm_set_U2S_matrix(TreeModel *mod, Vector *params, 
                               int start_idx) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  int nstates = mod->rate_matrix->size;
  int done_row[nstates];
  int setup_mapping = (lst_size(mod->rate_matrix_param_row[start_idx]) == 0);
  double rowsum;
  mat_zero(mod->rate_matrix->matrix);
  
  for (i = 0; i < nstates; i++) done_row[i] = 0;

  for (i = 0; i < nstates; i++) {
    double val;
    int b1_i, b2_i, b1_j, b2_j, b1_comp_i, b2_comp_i, i_comp, b1_comp_j, 
      b2_comp_j, j_comp;

    if (done_row[i]) continue;

    b1_i = i / alph_size;       /* idx of first char of "from" dinuc */
    b2_i = i % alph_size;       /* idx of second char */
    b1_comp_i =                 /* idx of compl of first char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b1_i])];
    b2_comp_i =                 /* idx of compl of second char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b2_i])];
    i_comp = b2_comp_i * alph_size + b1_comp_i;
                                /* state idx of compl dinuc */

    for (j = 0; j < nstates; j++) {

      if (j == i) continue;

      b1_j = j / alph_size;     /* idx of first char of "to" dinuc */
      b2_j = j % alph_size;     /* idx of second char */
      b1_comp_j =               /* idx of compl of first char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b1_j])];
      b2_comp_j =               /* idx of compl of second char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b2_j])];
      j_comp = b2_comp_j * alph_size + b1_comp_j;
                                /* state idx of compl dinuc */

      if (i == i_comp && j_comp < j) continue;
      if (b1_i != b1_j && b2_i != b2_j) continue;

      val = vec_get(params, start_idx);

      /* set rate for i,j */
      mm_set(mod->rate_matrix, i, j, val);
      
      /* set rate for its complement */
      mm_set(mod->rate_matrix, i_comp, j_comp, val);

      /* experimental */
      if (setup_mapping) {
        lst_push_int(mod->rate_matrix_param_row[start_idx], i);
        lst_push_int(mod->rate_matrix_param_col[start_idx], j);
        lst_push_int(mod->rate_matrix_param_row[start_idx], i_comp);
        lst_push_int(mod->rate_matrix_param_col[start_idx], j_comp);
      }

      start_idx++;
    }        
    done_row[i] = done_row[i_comp] = 1;
  }

  for (i = 0; i < nstates; i++) {
    rowsum = 0;    
    for (j = 0; j < nstates; j++)
      if (j != i) rowsum += mm_get(mod->rate_matrix, i, j);
    mm_set(mod->rate_matrix, i, i, -1 * rowsum);
  }
  if (start_idx != params->size)
    die("ERROR tm_set_U2S_matrix start_idx (%i) != params->size (%i)\n",
	start_idx, params->size);
}

/* FIXME: this should probably be generalized for a reversible matrix
   of arbitrary order */
void tm_set_R3_matrix(TreeModel *mod, Vector *params, int start_idx) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  int setup_mapping = (lst_size(mod->rate_matrix_param_row[start_idx]) == 0);
  mat_zero(mod->rate_matrix->matrix);

  for (i = 0; i < mod->rate_matrix->size; i++) {
    double rowsum = 0, val;
    int b1_i, b2_i, b3_i, b1_j, b2_j, b3_j;
    b1_i = i / (alph_size*alph_size);
    b2_i = (i % (alph_size*alph_size)) / alph_size;
    b3_i = i % alph_size;
    for (j = i+1; j < mod->rate_matrix->size; j++) {
      b1_j = j / (alph_size*alph_size);
      b2_j = (j % (alph_size*alph_size)) / alph_size;
      b3_j = j % alph_size;
      if ((b1_i != b1_j && b2_i != b2_j) || (b1_i != b1_j && b3_i != b3_j) || 
          (b2_i != b2_j && b3_i != b3_j))
        continue;

      val = vec_get(params, start_idx);
      mm_set(mod->rate_matrix, i, j, 
             val * vec_get(mod->backgd_freqs, j));
      mm_set(mod->rate_matrix, j, i, 
             val * vec_get(mod->backgd_freqs, i));
      rowsum += (val * vec_get(mod->backgd_freqs, j));

      /* experimental */
      if (setup_mapping) {
        lst_push_int(mod->rate_matrix_param_row[start_idx], i);
        lst_push_int(mod->rate_matrix_param_col[start_idx], j);
        lst_push_int(mod->rate_matrix_param_row[start_idx], j);
        lst_push_int(mod->rate_matrix_param_col[start_idx], i);
      }

      start_idx++;
    }        
    for (j = 0; j < i; j++)
      rowsum += mm_get(mod->rate_matrix, i, j);
    mm_set(mod->rate_matrix, i, i, -1 * rowsum);
  }
  if (start_idx != params->size)
    die("ERROR tm_set_RS_matrix start_idx (%i) != params->size (%i)\n",
	start_idx, params->size);
}

void tm_set_R3S_matrix(TreeModel *mod, Vector *params, int start_idx) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  int nstates = mod->rate_matrix->size;
  int setup_mapping = (lst_size(mod->rate_matrix_param_row[start_idx]) == 0);
  int done_row[nstates];
  double rowsum;
  mat_zero(mod->rate_matrix->matrix);

  for (i = 0; i < nstates; i++) done_row[i] = 0;

  for (i = 0; i < mod->rate_matrix->size; i++) {
    double val;
    int b1_i, b2_i, b3_i, b1_j, b2_j, b3_j, b1_comp_i, b2_comp_i, 
      b3_comp_i, i_comp, b1_comp_j, b2_comp_j, b3_comp_j, j_comp;

    if (done_row[i]) continue;

    b1_i = i / (alph_size*alph_size);
    b2_i = (i % (alph_size*alph_size)) / alph_size;
    b3_i = i % alph_size;

    b1_comp_i =                 /* idx of compl of first char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b1_i])];
    b2_comp_i =                 /* idx of compl of second char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b2_i])];
    b3_comp_i =                 /* idx of compl of second char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b3_i])];
    i_comp = b3_comp_i * alph_size * alph_size + 
      b2_comp_i * alph_size + b1_comp_i;
                                /* state idx of compl dinuc */

    for (j = i+1; j < mod->rate_matrix->size; j++) {
      b1_j = j / (alph_size*alph_size);
      b2_j = (j % (alph_size*alph_size)) / alph_size;
      b3_j = j % alph_size;

      b1_comp_j =               /* idx of compl of first char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b1_j])];
      b2_comp_j =               /* idx of compl of second char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b2_j])];
      b3_comp_j =               /* idx of compl of second char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b3_j])];
      j_comp = b3_comp_j * alph_size * alph_size +
        b2_comp_j * alph_size + b1_comp_j;
                                /* state idx of compl dinuc */

      if ((b1_i != b1_j && b2_i != b2_j) || (b1_i != b1_j && b3_i != b3_j) || 
          (b2_i != b2_j && b3_i != b3_j) || done_row[j_comp])
        continue;

      val = vec_get(params, start_idx);
      mm_set(mod->rate_matrix, i, j, 
             val * vec_get(mod->backgd_freqs, j));
      mm_set(mod->rate_matrix, j, i, 
             val * vec_get(mod->backgd_freqs, i));
      if (i != j_comp) {
        mm_set(mod->rate_matrix, i_comp, j_comp, 
               val * vec_get(mod->backgd_freqs, j_comp));
        mm_set(mod->rate_matrix, j_comp, i_comp, 
               val * vec_get(mod->backgd_freqs, i_comp));
      }

      /* experimental */
      if (setup_mapping) {
        lst_push_int(mod->rate_matrix_param_row[start_idx], i);
        lst_push_int(mod->rate_matrix_param_col[start_idx], j);
        lst_push_int(mod->rate_matrix_param_row[start_idx], j);
        lst_push_int(mod->rate_matrix_param_col[start_idx], i);
        if (i != j_comp) {
          lst_push_int(mod->rate_matrix_param_row[start_idx], i_comp);
          lst_push_int(mod->rate_matrix_param_col[start_idx], j_comp);
          lst_push_int(mod->rate_matrix_param_row[start_idx], j_comp);
          lst_push_int(mod->rate_matrix_param_col[start_idx], i_comp);
        }
      }

      start_idx++;
    }        
    done_row[i] = done_row[i_comp] = 1;
  }

  for (i = 0; i < nstates; i++) {
    rowsum = 0;    
    for (j = 0; j < nstates; j++)
      if (j != i) rowsum += mm_get(mod->rate_matrix, i, j);
    mm_set(mod->rate_matrix, i, i, -1 * rowsum);
  }
  if (start_idx != params->size)
    die("ERROR tm_set_R3S_matrix start_idx (%i) != params->size (%i)\n",
	start_idx, params->size);
}

void tm_set_U3_matrix(TreeModel *mod, Vector *params, int start_idx) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  int setup_mapping = (lst_size(mod->rate_matrix_param_row[start_idx]) == 0);
  mat_zero(mod->rate_matrix->matrix);

  for (i = 0; i < mod->rate_matrix->size; i++) {
    double rowsum = 0, val;
    int b1_i, b2_i, b3_i, b1_j, b2_j, b3_j;
    b1_i = i / (alph_size*alph_size);
    b2_i = (i % (alph_size*alph_size)) / alph_size;
    b3_i = i % alph_size;
    for (j = 0; j < mod->rate_matrix->size; j++) {
      if (i == j) continue;
      b1_j = j / (alph_size*alph_size);
      b2_j = (j % (alph_size*alph_size)) / alph_size;
      b3_j = j % alph_size;
      if ((b1_i != b1_j && b2_i != b2_j) || (b1_i != b1_j && b3_i != b3_j) || 
          (b2_i != b2_j && b3_i != b3_j))
        continue;

      val = vec_get(params, start_idx);
      mm_set(mod->rate_matrix, i, j, val);
      rowsum += val;

      /* experimental */
      if (setup_mapping) {
        lst_push_int(mod->rate_matrix_param_row[start_idx], i);
        lst_push_int(mod->rate_matrix_param_col[start_idx], j);
      }

      start_idx++;
    }        
    mm_set(mod->rate_matrix, i, i, -1 * rowsum);
  }
  if (start_idx != params->size)
    die("ERROR tm_set_US_matrix start_idx (%i) != params->size (%i)\n",
	start_idx, params->size);
}

void tm_set_U3S_matrix(TreeModel *mod, Vector *params, 
                               int start_idx) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  int nstates = mod->rate_matrix->size;
  int done_row[nstates];
  int setup_mapping = (lst_size(mod->rate_matrix_param_row[start_idx]) == 0);
  double rowsum;
  mat_zero(mod->rate_matrix->matrix);
  
  for (i = 0; i < nstates; i++) done_row[i] = 0;

  for (i = 0; i < nstates; i++) {
    double val;
    int b1_i, b2_i, b3_i, b1_j, b2_j, b3_j, b1_comp_i, b2_comp_i, 
      b3_comp_i, i_comp, b1_comp_j, b2_comp_j, b3_comp_j, j_comp;

    if (done_row[i]) continue;

    b1_i = i / (alph_size*alph_size);
    b2_i = (i % (alph_size*alph_size)) / alph_size;
    b3_i = i % alph_size;

    b1_comp_i =                 /* idx of compl of first char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b1_i])];
    b2_comp_i =                 /* idx of compl of second char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b2_i])];
    b3_comp_i =                 /* idx of compl of second char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b3_i])];
    i_comp = b3_comp_i * alph_size * alph_size + 
      b2_comp_i * alph_size + b1_comp_i;
                                /* state idx of compl dinuc */

    for (j = 0; j < nstates; j++) {

      if (j == i) continue;

      b1_j = j / (alph_size*alph_size);
      b2_j = (j % (alph_size*alph_size)) / alph_size;
      b3_j = j % alph_size;
      b1_comp_j =               /* idx of compl of first char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b1_j])];
      b2_comp_j =               /* idx of compl of second char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b2_j])];
      b3_comp_j =               /* idx of compl of second char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b3_j])];
      j_comp = b3_comp_j * alph_size * alph_size +
        b2_comp_j * alph_size + b1_comp_j;
                                /* state idx of compl dinuc */

      if ((b1_i != b1_j && b2_i != b2_j) || (b1_i != b1_j && b3_i != b3_j) || 
          (b2_i != b2_j && b3_i != b3_j))
        continue;

      val = vec_get(params, start_idx);

      /* set rate for i,j */
      mm_set(mod->rate_matrix, i, j, val);
      
      /* set rate for its complement */
      mm_set(mod->rate_matrix, i_comp, j_comp, val);

      /* experimental */
      if (setup_mapping) {
        lst_push_int(mod->rate_matrix_param_row[start_idx], i);
        lst_push_int(mod->rate_matrix_param_col[start_idx], j);
        lst_push_int(mod->rate_matrix_param_row[start_idx], i_comp);
        lst_push_int(mod->rate_matrix_param_col[start_idx], j_comp);
      }

      start_idx++;
    }        
    done_row[i] = done_row[i_comp] = 1;
  }

  for (i = 0; i < nstates; i++) {
    rowsum = 0;    
    for (j = 0; j < nstates; j++)
      if (j != i) rowsum += mm_get(mod->rate_matrix, i, j);
    mm_set(mod->rate_matrix, i, i, -1 * rowsum);
  }
  if (start_idx != params->size)
    die("ERROR tm_set_U3S_matrix start_idx (%i) != params->size (%i)\n",
	start_idx, params->size);
}

/* void tm_set_GY_matrix(TreeModel *mod, double kappa, double omega) { */
/*   int i, j; */
/*   int alph_size = strlen(mod->rate_matrix->states); */
/*   int setup_mapping = (lst_size(mod->rate_matrix_param_row[start_idx]) == 0); */
/*   char *aa = get_codon_mapping(mod->rate_matrix->states); */

/*   mat_zero(mod->rate_matrix->matrix); */

/*   for (i = 0; i < mod->rate_matrix->size; i++) { */
/*     double rowsum = 0, val; */
/*     int b1_i, b2_i, b3_i, b1_j, b2_j, b3_j; */
/*     b1_i = i / (alph_size*alph_size); */
/*     b2_i = (i % (alph_size*alph_size)) / alph_size; */
/*     b3_i = i % alph_size; */
/*     for (j = i+1; j < mod->rate_matrix->size; j++) { */
/*       b1_j = j / (alph_size*alph_size); */
/*       b2_j = (j % (alph_size*alph_size)) / alph_size; */
/*       b3_j = j % alph_size; */
/*       if ((b1_i != b1_j && b2_i != b2_j) || (b1_i != b1_j && b3_i != b3_j) ||  */
/*           (b2_i != b2_j && b3_i != b3_j)) */
/*         continue; */

/*       val = vec_get(params, start_idx); */
/*       mm_set(mod->rate_matrix, i, j,  */
/*              val * vec_get(mod->backgd_freqs, j)); */
/*       mm_set(mod->rate_matrix, j, i,  */
/*              val * vec_get(mod->backgd_freqs, i)); */
/*       rowsum += (val * vec_get(mod->backgd_freqs, j)); */

      /* experimental */
/*       if (setup_mapping) { */
/*         lst_push_int(mod->rate_matrix_param_row[start_idx], i); */
/*         lst_push_int(mod->rate_matrix_param_col[start_idx], j); */
/*         lst_push_int(mod->rate_matrix_param_row[start_idx], j); */
/*         lst_push_int(mod->rate_matrix_param_col[start_idx], i); */
/*       } */

/*       start_idx++; */
/*     }         */
/*     for (j = 0; j < i; j++) */
/*       rowsum += mm_get(mod->rate_matrix, i, j); */
/*     mm_set(mod->rate_matrix, i, i, -1 * rowsum); */
/*   } */
/*  if (start_idx != params_size) */
/*    die("ERROR tm_set_GY_matrix start_idx (%i) != params->size (%i)\n", */
/*	start_idx, params->size); */
/*   free(aa); */
/* } */

/***************************************************************************/
/* Functions to initialize rate-matrix parameters to reasonable
   values, one function per substitution model.  Matrices are either
   initialized using an HKY-like scheme, based on a
   transition/transversion bias parameter kappa, or are initialized
   from existing models. */
/***************************************************************************/

/* FIXME: the functions below should be able to have access to backgd
   freqs but currently they are uninitialized at the time of
   invocation */

/* initialize REV as if HKY */
void tm_init_mat_REV(TreeModel *mod, Vector *params, int parm_idx, 
                     double kappa) {
  int i, j;
  for (i = 0; i < mod->rate_matrix->size; i++) {
    for (j = i+1; j < mod->rate_matrix->size; j++) {
      double val = 1;
      if (is_transition(mod->rate_matrix->states[i], 
                        mod->rate_matrix->states[j])) 
        val *= kappa;
      vec_set(params, parm_idx++, val + .05 * unif_rand());
                                /* add a little noise to initial
                                   values to avoid making matrix
                                   undiagonalizable (becomes important
                                   with gap-as-fifth-char) */
    }
  }    
}

/* initialize SSREV as if HKY */
void tm_init_mat_SSREV(TreeModel *mod, Vector *params, int parm_idx, 
		       double kappa) {
  int i, j, compi, compj;
  int count=0;  //testing
  for (i = 0; i < mod->rate_matrix->size; i++) {
    compi=mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->states[i])];
    for (j = i+1; j < mod->rate_matrix->size; j++) {
      double val = 1;
      compj=mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->states[j])];

      if ((compi < compj && compi < i) ||
	  (compj < compi && compj < i)) continue;

      if (is_transition(mod->rate_matrix->states[i], 
                        mod->rate_matrix->states[j])) 
        val *= kappa;
      vec_set(params, parm_idx++, val + .05 * unif_rand());
      count++;
                                /* add a little noise to initial
                                   values to avoid making matrix
                                   undiagonalizable (becomes important
                                   with gap-as-fifth-char) */
    }
  }    
  if (count != tm_get_nratematparams(mod))
    die("ERROR tm_init_mat_SSERV count (%i) != tm_get_nratematparams(mod) (%i)\n",
	count, tm_get_nratematparams(mod));
}


void tm_init_mat_REV_GC(TreeModel *mod, Vector *params, int parm_idx,
			 double kappa) {
  vec_set(params, parm_idx, 2.0);  //init GC to 1?
  tm_init_mat_REV(mod, params, parm_idx+1, kappa);
}


/* initialize UNREST as if HKY */
void tm_init_mat_UNREST(TreeModel *mod, Vector *params, int parm_idx, 
                     double kappa) {
  int i, j;
  for (i = 0; i < mod->rate_matrix->size; i++) {
    for (j = 0; j < mod->rate_matrix->size; j++) {
      double val = 1;
      if (i == j) continue;
      if (is_transition(mod->rate_matrix->states[i], 
                           mod->rate_matrix->states[j])) 
        val *= kappa;
      vec_set(params, parm_idx++, val + .05 * unif_rand());
                                /* add a little noise to initial
                                   values to avoid making matrix
                                   undiagonalizable :) */
    }
  }    
}

/* initialize R2 as if HKY2 */
void tm_init_mat_R2(TreeModel *mod, Vector *params, int parm_idx, 
                      double kappa) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  for (i = 0; i < mod->rate_matrix->size; i++) {
    int b1_i, b2_i, b1_j, b2_j;
    b1_i = i / alph_size;
    b2_i = i % alph_size;
    for (j = i+1; j < mod->rate_matrix->size; j++) {
      double val = 1;
      b1_j = j / alph_size;
      b2_j = j % alph_size;
      if (b1_i != b1_j && b2_i != b2_j) continue;
      if (is_transition(mod->rate_matrix->states[b1_i], 
                           mod->rate_matrix->states[b1_j]) ||
          is_transition(mod->rate_matrix->states[b2_i], 
                           mod->rate_matrix->states[b2_j])) 
        val *= kappa;
      vec_set(params, parm_idx++, val);
    }
  }    
  if (parm_idx != params->size)
    die("ERROR tm_init_mat_R2: parm_idx (%i) != params->size (%i)\n",
	parm_idx, params->size);
}

/* initialize R2S as if HKY2 */
void tm_init_mat_R2S(TreeModel *mod, Vector *params, 
                          int parm_idx, double kappa) {
  int i;
  int alph_size = strlen(mod->rate_matrix->states);
  int nstates = mod->rate_matrix->size;
  int done_row[nstates];

  for (i = 0; i < nstates; i++) done_row[i] = 0;
  for (i = 0; i < nstates; i++) {
    double val = 1;
    int j, b1_i, b2_i, b1_j, b2_j, b1_comp_i, b2_comp_i, i_comp, b1_comp_j, 
      b2_comp_j, j_comp;

    if (done_row[i]) continue;

    b1_i = i / alph_size;       /* idx of first char of "from" dinuc */
    b2_i = i % alph_size;       /* idx of second char */
    b1_comp_i =                 /* idx of compl of first char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b1_i])];
    b2_comp_i =                 /* idx of compl of second char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b2_i])];
    i_comp = b2_comp_i * alph_size + b1_comp_i;
                                /* state idx of compl dinuc */

    for (j = i+1; j < nstates; j++) {

      b1_j = j / alph_size;     /* idx of first char of "to" dinuc */
      b2_j = j % alph_size;     /* idx of second char */
      b1_comp_j =               /* idx of compl of first char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b1_j])];
      b2_comp_j =               /* idx of compl of second char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b2_j])];
      j_comp = b2_comp_j * alph_size + b1_comp_j;
                                /* state idx of compl dinuc */

      if ((b1_i != b1_j && b2_i != b2_j) || done_row[j_comp] || 
          (i == i_comp && j_comp < j)) continue;

      val = 1;
      if (is_transition(mod->rate_matrix->states[b1_i], 
                           mod->rate_matrix->states[b1_j]) ||
          is_transition(mod->rate_matrix->states[b2_i], 
                           mod->rate_matrix->states[b2_j])) 
        val *= kappa;

      vec_set(params, parm_idx++, val);
    }        
    done_row[i] = done_row[i_comp] = 1;
  }
  if (parm_idx != params->size)
    die("ERROR tm_init_mat_R2S: parm_idx (%i) != params->size (%i)\n",
	parm_idx, params->size);
}

/* initialize U2 as if HKY2 */
void tm_init_mat_U2(TreeModel *mod, Vector *params, 
                         int parm_idx, double kappa) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  for (i = 0; i < mod->rate_matrix->size; i++) {
    int b1_i, b2_i, b1_j, b2_j;
    b1_i = i / alph_size;
    b2_i = i % alph_size;
    for (j = 0; j < mod->rate_matrix->size; j++) {
      double val = 1;
      b1_j = j / alph_size;
      b2_j = j % alph_size;
      if (i == j || (b1_i != b1_j && b2_i != b2_j)) continue;
      if (is_transition(mod->rate_matrix->states[b1_i], 
                           mod->rate_matrix->states[b1_j]) ||
          is_transition(mod->rate_matrix->states[b2_i], 
                           mod->rate_matrix->states[b2_j])) 
        val *= kappa;
      vec_set(params, parm_idx++, val + .05 * unif_rand());
                                /* add a little noise to initial
                                   values to avoid making matrix
                                   undiagonalizable :) */
    }
  }    
  if (parm_idx != params->size)
    die("ERROR tm_init_mat_U2: parm_idx (%i) != params->size (%i)\n",
	parm_idx, params->size);
}

/* initialize U2S as if HKY2 */
void tm_init_mat_U2S(TreeModel *mod, Vector *params, 
                             int parm_idx, double kappa) {
  int i;
  int alph_size = strlen(mod->rate_matrix->states);
  int nstates = mod->rate_matrix->size;
  int done_row[nstates];

  for (i = 0; i < nstates; i++) done_row[i] = 0;
  for (i = 0; i < nstates; i++) {
    double val = 1;
    int j, b1_i, b2_i, b1_j, b2_j, b1_comp_i, b2_comp_i, i_comp, b1_comp_j, 
      b2_comp_j, j_comp;

    if (done_row[i]) continue;

    b1_i = i / alph_size;       /* idx of first char of "from" dinuc */
    b2_i = i % alph_size;       /* idx of second char */
    b1_comp_i =                 /* idx of compl of first char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b1_i])];
    b2_comp_i =                 /* idx of compl of second char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b2_i])];
    i_comp = b2_comp_i * alph_size + b1_comp_i;
                                /* state idx of compl dinuc */

    for (j = 0; j < nstates; j++) {

      if (j == i) continue;

      b1_j = j / alph_size;     /* idx of first char of "to" dinuc */
      b2_j = j % alph_size;     /* idx of second char */
      b1_comp_j =               /* idx of compl of first char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b1_j])];
      b2_comp_j =               /* idx of compl of second char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b2_j])];
      j_comp = b2_comp_j * alph_size + b1_comp_j;
                                /* state idx of compl dinuc */

      if ((b1_i != b1_j && b2_i != b2_j) || 
          (i == i_comp && j_comp < j)) continue;

      val = 1;
      if (is_transition(mod->rate_matrix->states[b1_i], 
                           mod->rate_matrix->states[b1_j]) ||
          is_transition(mod->rate_matrix->states[b2_i], 
                           mod->rate_matrix->states[b2_j])) 
        val *= kappa;

      vec_set(params, parm_idx++, val + .05 * unif_rand());
                                /* add a little noise to initial
                                   values to avoid making matrix
                                   undiagonalizable :) */
    }        
    done_row[i] = done_row[i_comp] = 1;
  }
  if (parm_idx != params->size)
    die("ERROR tm_init_mat_U2S: parm_idx (%i) != params->size (%i)\n",
	parm_idx, params->size);
}

void tm_init_mat_R3(TreeModel *mod, Vector *params, int parm_idx, 
                    double kappa) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  for (i = 0; i < mod->rate_matrix->size; i++) {
    int b1_i, b2_i, b3_i, b1_j, b2_j, b3_j;
    b1_i = i / (alph_size*alph_size);
    b2_i = (i % (alph_size*alph_size)) / alph_size;
    b3_i = i % alph_size;
    for (j = i+1; j < mod->rate_matrix->size; j++) {
      double val = 1;
      b1_j = j / (alph_size*alph_size);
      b2_j = (j % (alph_size*alph_size)) / alph_size;
      b3_j = j % alph_size;
      if ((b1_i != b1_j && b2_i != b2_j) || (b1_i != b1_j && b3_i != b3_j) || 
          (b2_i != b2_j && b3_i != b3_j))
        continue;
      if (is_transition(mod->rate_matrix->states[b1_i], 
                           mod->rate_matrix->states[b1_j]) ||
          is_transition(mod->rate_matrix->states[b2_i], 
                           mod->rate_matrix->states[b2_j]) ||
          is_transition(mod->rate_matrix->states[b3_i], 
                           mod->rate_matrix->states[b3_j])) 
        val *= kappa;
      vec_set(params, parm_idx++, val);
    }
  }
  if (parm_idx != params->size)
    die("ERROR tm_init_mat_RS: parm_idx (%i) != params->size (%i)\n",
	parm_idx, params->size);
}

void tm_init_mat_R3S(TreeModel *mod, Vector *params, int parm_idx, 
                     double kappa) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  int nstates = mod->rate_matrix->size;
  int done_row[nstates];

  for (i = 0; i < nstates; i++) done_row[i] = 0;
  for (i = 0; i < mod->rate_matrix->size; i++) {
    int b1_i, b2_i, b3_i, b1_j, b2_j, b3_j;
    int b1_comp_i, b2_comp_i, b3_comp_i, b1_comp_j, b2_comp_j, b3_comp_j,
      i_comp, j_comp;

    if (done_row[i]) continue;

    b1_i = i / (alph_size*alph_size);
    b2_i = (i % (alph_size*alph_size)) / alph_size;
    b3_i = i % alph_size;

    b1_comp_i =                 /* idx of compl of first char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b1_i])];
    b2_comp_i =                 /* idx of compl of second char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b2_i])];
    b3_comp_i =                 /* idx of compl of third char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b3_i])];
    i_comp = b3_comp_i * alph_size * alph_size + 
      b2_comp_i * alph_size + b1_comp_i;
                                /* state idx of compl dinuc */

    for (j = i+1; j < mod->rate_matrix->size; j++) {
      double val = 1;
      b1_j = j / (alph_size*alph_size);
      b2_j = (j % (alph_size*alph_size)) / alph_size;
      b3_j = j % alph_size;

      b1_comp_j =               /* idx of compl of first char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b1_j])];
      b2_comp_j =               /* idx of compl of second char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b2_j])];
      b3_comp_j =               /* idx of compl of second char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b3_j])];
      j_comp = b3_comp_j * alph_size * alph_size + 
        b2_comp_j * alph_size + b1_comp_j;
                                /* state idx of compl dinuc */

      if ((b1_i != b1_j && b2_i != b2_j) || (b1_i != b1_j && b3_i != b3_j) || 
          (b2_i != b2_j && b3_i != b3_j) || done_row[j_comp])
        continue;

      if (is_transition(mod->rate_matrix->states[b1_i], 
                           mod->rate_matrix->states[b1_j]) ||
          is_transition(mod->rate_matrix->states[b2_i], 
                           mod->rate_matrix->states[b2_j]) ||
          is_transition(mod->rate_matrix->states[b3_i], 
                           mod->rate_matrix->states[b3_j])) 
        val *= kappa;
      vec_set(params, parm_idx++, val);
    }
    done_row[i] = done_row[i_comp] = 1;
  }
  if (parm_idx != params->size)
    die("ERROR tm_init_mat_R3S: parm_idx (%i) != params->size (%i)\n",
	parm_idx, params->size);
}

void tm_init_mat_U3(TreeModel *mod, Vector *params, 
                    int parm_idx, double kappa) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  for (i = 0; i < mod->rate_matrix->size; i++) {
    int b1_i, b2_i, b3_i, b1_j, b2_j, b3_j;
    b1_i = i / (alph_size*alph_size);
    b2_i = (i % (alph_size * alph_size)) / alph_size;
    b3_i = i % alph_size;
    for (j = 0; j < mod->rate_matrix->size; j++) {
      double val = 1.0 / mod->rate_matrix->size;
      b1_j = j / (alph_size*alph_size);
      b2_j = (j % (alph_size*alph_size)) / alph_size;
      b3_j = j % alph_size;
      if (i == j || (b1_i != b1_j && b2_i != b2_j) || 
          (b1_i != b1_j && b3_i != b3_j) || (b2_i != b2_j && b3_i != b3_j))
        continue;

      if (is_transition(mod->rate_matrix->states[b1_i], 
                           mod->rate_matrix->states[b1_j]) ||
          is_transition(mod->rate_matrix->states[b2_i], 
                           mod->rate_matrix->states[b2_j]) ||
          is_transition(mod->rate_matrix->states[b3_i], 
                           mod->rate_matrix->states[b3_j])) 
        val *= kappa;
      vec_set(params, parm_idx++, val + .001 * unif_rand());
                                /* add a little noise to initial
                                   values to avoid making matrix
                                   undiagonalizable :) */
    }
  }    
  if (parm_idx != params->size)
    die("ERROR tm_init_mat_US: parm_idx (%i) != params->size (%i)\n",
	parm_idx, params->size);
}

void tm_init_mat_U3S(TreeModel *mod, Vector *params, 
                     int parm_idx, double kappa) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  int nstates = mod->rate_matrix->size;
  int done_row[nstates];

  for (i = 0; i < nstates; i++) done_row[i] = 0;
  for (i = 0; i < nstates; i++) {
    int b1_i, b2_i, b3_i, b1_j, b2_j, b3_j;
    int b1_comp_i, b2_comp_i, b3_comp_i, b1_comp_j, b2_comp_j, b3_comp_j,
      i_comp, j_comp;

    if (done_row[i]) continue;

    b1_i = i / (alph_size*alph_size);
    b2_i = (i % (alph_size * alph_size)) / alph_size;
    b3_i = i % alph_size;

    b1_comp_i =                 /* idx of compl of first char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b1_i])];
    b2_comp_i =                 /* idx of compl of second char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b2_i])];
    b3_comp_i =                 /* idx of compl of third char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b3_i])];
    i_comp = b3_comp_i * alph_size * alph_size + 
      b2_comp_i * alph_size + b1_comp_i;
                                /* state idx of compl dinuc */

    for (j = 0; j < nstates; j++) {
      double val = 1.0 / nstates;
      if (j == i) continue;

      b1_j = j / (alph_size*alph_size);
      b2_j = (j % (alph_size*alph_size)) / alph_size;
      b3_j = j % alph_size;

      b1_comp_j =               /* idx of compl of first char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b1_j])];
      b2_comp_j =               /* idx of compl of second char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b2_j])];
      b3_comp_j =               /* idx of compl of second char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b3_j])];
      j_comp = b3_comp_j * alph_size * alph_size + 
        b2_comp_j * alph_size + b1_comp_j;
                                /* state idx of compl dinuc */

      if ((b1_i != b1_j && b2_i != b2_j) || (b1_i != b1_j && b3_i != b3_j) || 
          (b2_i != b2_j && b3_i != b3_j)) 
        continue;

      if (is_transition(mod->rate_matrix->states[b1_i], 
                           mod->rate_matrix->states[b1_j]) ||
          is_transition(mod->rate_matrix->states[b2_i], 
                           mod->rate_matrix->states[b2_j]) ||
          is_transition(mod->rate_matrix->states[b3_i], 
                           mod->rate_matrix->states[b3_j])) 
        val *= kappa;
      vec_set(params, parm_idx++, val + .001 * unif_rand());
                                /* add a little noise to initial
                                   values to avoid making matrix
                                   undiagonalizable :) */
    }
    done_row[i] = done_row[i_comp] = 1;
  }
  if (parm_idx != params->size)
    die("ERROR tm_init_mat_U3S: parm_idx (%i) != params->size (%i)\n",
	parm_idx, params->size);
}

void tm_init_mat_GY(TreeModel *mod, double kappa, double omega) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  char *aa = get_codon_mapping(mod->rate_matrix->states);
  for (i = 0; i < mod->rate_matrix->size; i++) {
    int b1_i, b2_i, b3_i, b1_j, b2_j, b3_j;
    b1_i = i / (alph_size*alph_size);
    b2_i = (i % (alph_size*alph_size)) / alph_size;
    b3_i = i % alph_size;
    for (j = i+1; j < mod->rate_matrix->size; j++) {
      double val = 1;
      b1_j = j / (alph_size*alph_size);
      b2_j = (j % (alph_size*alph_size)) / alph_size;
      b3_j = j % alph_size;
      if ((b1_i != b1_j && b2_i != b2_j) || (b1_i != b1_j && b3_i != b3_j) || 
          (b2_i != b2_j && b3_i != b3_j))
        continue;
      if (is_transition(mod->rate_matrix->states[b1_i], 
                           mod->rate_matrix->states[b1_j]) ||
          is_transition(mod->rate_matrix->states[b2_i], 
                           mod->rate_matrix->states[b2_j]) ||
          is_transition(mod->rate_matrix->states[b3_i], 
                           mod->rate_matrix->states[b3_j])) 
        val *= kappa;

      if (aa[i] != aa[j]) val *= omega;
    }
  }    
  free(aa);
}

void tm_init_mat_from_model_REV(TreeModel *mod, Vector *params, 
                                int start_idx) {
  int i, j;
  for (i = 0; i < mod->rate_matrix->size; i++) {
    for (j = i+1; j < mod->rate_matrix->size; j++) {
      double val = safediv(mm_get(mod->rate_matrix, i, j),
                           vec_get(mod->backgd_freqs, j));
      vec_set(params, start_idx++, val);
    }
  }
}


void tm_init_mat_from_model_SSREV(TreeModel *mod, Vector *params, 
				  int start_idx) {
  int i, j, compi, compj;
  double val;
  int count=0; //testing
  for (i = 0; i < mod->rate_matrix->size; i++) {
    compi=mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->states[i])];
    for (j = i+1; j < mod->rate_matrix->size; j++) {
      compj=mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->states[j])];
      if ((compi < compj && compi < i) ||
	  (compj < compi && compj < i)) continue;
      val = safediv(mm_get(mod->rate_matrix, i, j),
		    vec_get(mod->backgd_freqs, j));
      vec_set(params, start_idx++, val);
      count++;
    }
  }
  if (count != tm_get_nratematparams(mod))
    die("ERROR tm_init_mat_from_model_SSREV: count (%i) != tm_get_nratematparams(mod) (%i)\n",
	count, tm_get_nratematparams(mod));
}



void tm_init_mat_from_model_UNREST(TreeModel *mod, Vector *params, 
                                   int start_idx) {
  int i, j;
  for (i = 0; i < mod->rate_matrix->size; i++) {
    for (j = 0; j < mod->rate_matrix->size; j++) {
      double val;
      if (i == j) continue;
      val = mm_get(mod->rate_matrix, i, j);
      vec_set(params, start_idx++, val);
    }
  }
}

void tm_init_mat_from_model_R2(TreeModel *mod, Vector *params, 
                               int start_idx) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  for (i = 0; i < mod->rate_matrix->size; i++) {
    int b1_i, b2_i, b1_j, b2_j;
    b1_i = i / alph_size;
    b2_i = i % alph_size;
    for (j = i+1; j < mod->rate_matrix->size; j++) {
      double val;
      b1_j = j / alph_size;
      b2_j = j % alph_size;
      if (b1_i != b1_j && b2_i != b2_j) continue;
      val = safediv(mm_get(mod->rate_matrix, i, j),
                    vec_get(mod->backgd_freqs, j));
      vec_set(params, start_idx++, val);
    }
  }
  if (start_idx != params->size)
    die("ERROR tm_init_mat_from_model_R2: start_idx (%i) != params->size (%i)\n",
	start_idx, params->size);
}

void tm_init_mat_from_model_U2(TreeModel *mod, Vector *params, 
                               int start_idx) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  for (i = 0; i < mod->rate_matrix->size; i++) {
    int b1_i, b2_i, b1_j, b2_j;
    b1_i = i / alph_size;
    b2_i = i % alph_size;
    for (j = 0; j < mod->rate_matrix->size; j++) {
      double val;
      b1_j = j / alph_size;
      b2_j = j % alph_size;
      if (i == j || (b1_i != b1_j && b2_i != b2_j)) continue;
      val = mm_get(mod->rate_matrix, i, j);
      vec_set(params, start_idx++, val);
    }
  }
  if (start_idx != params->size)
    die("ERROR tm_init_mat_from_model_U2: start_idx (%i) != params->size (%i)\n",
	start_idx, params->size);
}

void tm_init_mat_from_model_R2S(TreeModel *mod, Vector *params, 
                                int start_idx) {
  int i;
  int alph_size = strlen(mod->rate_matrix->states);
  int nstates = mod->rate_matrix->size;
  int done_row[nstates];

  for (i = 0; i < nstates; i++) done_row[i] = 0;
  for (i = 0; i < nstates; i++) {
    double val;
    int j, b1_i, b2_i, b1_j, b2_j, b1_comp_i, b2_comp_i, i_comp, b1_comp_j, 
      b2_comp_j, j_comp;

    if (done_row[i]) continue;

    b1_i = i / alph_size;       /* idx of first char of "from" dinuc */
    b2_i = i % alph_size;       /* idx of second char */
    b1_comp_i =                 /* idx of compl of first char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b1_i])];
    b2_comp_i =                 /* idx of compl of second char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b2_i])];
    i_comp = b2_comp_i * alph_size + b1_comp_i;
                                /* state idx of compl dinuc */

    for (j = i+1; j < nstates; j++) {

      b1_j = j / alph_size;     /* idx of first char of "to" dinuc */
      b2_j = j % alph_size;     /* idx of second char */
      b1_comp_j =               /* idx of compl of first char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b1_j])];
      b2_comp_j =               /* idx of compl of second char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b2_j])];
      j_comp = b2_comp_j * alph_size + b1_comp_j;
                                /* state idx of compl dinuc */

      if ((b1_i != b1_j && b2_i != b2_j) || done_row[j_comp] || 
          (i == i_comp && j_comp < j)) continue;

      val = safediv(mm_get(mod->rate_matrix, i, j),
                    vec_get(mod->backgd_freqs, j));
      vec_set(params, start_idx++, val);
    }        
    done_row[i] = done_row[i_comp] = 1;
  }
  if (start_idx != params->size)
    die("ERROR tm_init_mat_from_model_R2S: start_idx (%i) != params->size (%i)\n",
	start_idx, params->size);
}

void tm_init_mat_from_model_U2S(TreeModel *mod, Vector *params, 
                                int start_idx) {
  int i;
  int alph_size = strlen(mod->rate_matrix->states);
  int nstates = mod->rate_matrix->size;
  int done_row[nstates];

  for (i = 0; i < nstates; i++) done_row[i] = 0;
  for (i = 0; i < nstates; i++) {
    double val;
    int j, b1_i, b2_i, b1_j, b2_j, b1_comp_i, b2_comp_i, i_comp, b1_comp_j, 
      b2_comp_j, j_comp;

    if (done_row[i]) continue;

    b1_i = i / alph_size;       /* idx of first char of "from" dinuc */
    b2_i = i % alph_size;       /* idx of second char */
    b1_comp_i =                 /* idx of compl of first char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b1_i])];
    b2_comp_i =                 /* idx of compl of second char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b2_i])];
    i_comp = b2_comp_i * alph_size + b1_comp_i;
                                /* state idx of compl dinuc */

    for (j = 0; j < nstates; j++) {
      if (j == i) continue;
      b1_j = j / alph_size;     /* idx of first char of "to" dinuc */
      b2_j = j % alph_size;     /* idx of second char */
      b1_comp_j =               /* idx of compl of first char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b1_j])];
      b2_comp_j =               /* idx of compl of second char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b2_j])];
      j_comp = b2_comp_j * alph_size + b1_comp_j;
                                /* state idx of compl dinuc */

      if ((b1_i != b1_j && b2_i != b2_j) || 
          (i == i_comp && j_comp < j)) continue;

      val = mm_get(mod->rate_matrix, i, j);
      vec_set(params, start_idx++, val);
    }        
    done_row[i] = done_row[i_comp] = 1;
  }
  if (start_idx != params->size)
    die("ERROR tm_init_mat_from_model_U2S: start_idx (%i) != params->size (%i)\n",
	start_idx, params->size);
}

void tm_init_mat_from_model_R3(TreeModel *mod, Vector *params, 
                               int start_idx) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  for (i = 0; i < mod->rate_matrix->size; i++) {
    int b1_i, b2_i, b3_i, b1_j, b2_j, b3_j;
    b1_i = i / (alph_size*alph_size);
    b2_i = (i % (alph_size*alph_size)) / alph_size;
    b3_i = i % alph_size;
    for (j = i+1; j < mod->rate_matrix->size; j++) {
      double val;
      b1_j = j / (alph_size*alph_size);
      b2_j = (j % (alph_size*alph_size)) / alph_size;
      b3_j = j % alph_size;
      if ((b1_i != b1_j && b2_i != b2_j) || (b1_i != b1_j && b3_i != b3_j) || 
          (b2_i != b2_j && b3_i != b3_j))
        continue;
      val = safediv(mm_get(mod->rate_matrix, i, j),
                    vec_get(mod->backgd_freqs, j));
      vec_set(params, start_idx++, val);
    }
  }
  if (start_idx != params->size)
    die("ERROR tm_init_mat_from_model_R3: start_idx (%i) != params->size (%i)\n",
	start_idx, params->size);
}

void tm_init_mat_from_model_R3S(TreeModel *mod, Vector *params, 
                               int start_idx) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  int nstates = mod->rate_matrix->size;
  int done_row[nstates];
  for (i = 0; i < nstates; i++) done_row[i] = 0;
  for (i = 0; i < mod->rate_matrix->size; i++) {
    int b1_i, b2_i, b3_i, b1_j, b2_j, b3_j;
    int b1_comp_i, b2_comp_i, b3_comp_i, b1_comp_j, b2_comp_j, b3_comp_j,
      i_comp, j_comp;

    if (done_row[i]) continue;

    b1_i = i / (alph_size*alph_size);
    b2_i = (i % (alph_size*alph_size)) / alph_size;
    b3_i = i % alph_size;

    b1_comp_i =                 /* idx of compl of first char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b1_i])];
    b2_comp_i =                 /* idx of compl of second char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b2_i])];
    b3_comp_i =                 /* idx of compl of third char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b3_i])];
    i_comp = b3_comp_i * alph_size * alph_size + 
      b2_comp_i * alph_size + b1_comp_i;

    for (j = i+1; j < mod->rate_matrix->size; j++) {
      double val;
      b1_j = j / (alph_size*alph_size);
      b2_j = (j % (alph_size*alph_size)) / alph_size;
      b3_j = j % alph_size;

      b1_comp_j =               /* idx of compl of first char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b1_j])];
      b2_comp_j =               /* idx of compl of second char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b2_j])];
      b3_comp_j =               /* idx of compl of second char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b3_j])];
      j_comp = b3_comp_j * alph_size * alph_size + 
        b2_comp_j * alph_size + b1_comp_j;
                                /* state idx of compl dinuc */

      if ((b1_i != b1_j && b2_i != b2_j) || (b1_i != b1_j && b3_i != b3_j) || 
          (b2_i != b2_j && b3_i != b3_j) || done_row[j_comp])
        continue;

      val = safediv(mm_get(mod->rate_matrix, i, j),
                    vec_get(mod->backgd_freqs, j));
      vec_set(params, start_idx++, val);
    }
    done_row[i] = done_row[i_comp] = 1;
  }
  if (start_idx != params->size)
    die("ERROR tm_init_mat_from_model_R3S: start_idx (%i) != params->size (%i)\n",
	start_idx, params->size);
}

void tm_init_mat_from_model_U3(TreeModel *mod, Vector *params, 
                               int start_idx) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  for (i = 0; i < mod->rate_matrix->size; i++) {
    int b1_i, b2_i, b3_i, b1_j, b2_j, b3_j;
    b1_i = i / (alph_size*alph_size);
    b2_i = (i % (alph_size*alph_size)) / alph_size;
    b3_i = i % alph_size;
    for (j = 0; j < mod->rate_matrix->size; j++) {
      double val;
      b1_j = j / (alph_size*alph_size);
      b2_j = (j % (alph_size*alph_size)) / alph_size;
      b3_j = j % alph_size;
      if (i == j || 
          (b1_i != b1_j && b2_i != b2_j) || (b1_i != b1_j && b3_i != b3_j) || 
          (b2_i != b2_j && b3_i != b3_j))
        continue;
      val = mm_get(mod->rate_matrix, i, j);
      vec_set(params, start_idx++, val);
    }
  }
  if (start_idx != params->size)
    die("ERROR tm_init_mat_from_model_U3: start_idx (%i) != params->size (%i)\n",
	start_idx, params->size);
}

void tm_init_mat_from_model_U3S(TreeModel *mod, Vector *params, 
                               int start_idx) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  int nstates = mod->rate_matrix->size;
  int done_row[nstates];

  for (i = 0; i < nstates; i++) done_row[i] = 0;
  for (i = 0; i < nstates; i++) {
    int b1_i, b2_i, b3_i, b1_j, b2_j, b3_j;
    int b1_comp_i, b2_comp_i, b3_comp_i, b1_comp_j, b2_comp_j, b3_comp_j,
      i_comp, j_comp;

    if (done_row[i]) continue;

    b1_i = i / (alph_size*alph_size);
    b2_i = (i % (alph_size * alph_size)) / alph_size;
    b3_i = i % alph_size;

    b1_comp_i =                 /* idx of compl of first char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b1_i])];
    b2_comp_i =                 /* idx of compl of second char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b2_i])];
    b3_comp_i =                 /* idx of compl of third char */
      mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                       states[b3_i])];
    i_comp = b3_comp_i * alph_size * alph_size + 
      b2_comp_i * alph_size + b1_comp_i;
                                /* state idx of compl dinuc */

    for (j = 0; j < nstates; j++) {
      double val;
      if (j == i) continue;

      b1_j = j / (alph_size*alph_size);
      b2_j = (j % (alph_size*alph_size)) / alph_size;
      b3_j = j % alph_size;

      b1_comp_j =               /* idx of compl of first char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b1_j])];
      b2_comp_j =               /* idx of compl of second char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b2_j])];
      b3_comp_j =               /* idx of compl of second char */
        mod->rate_matrix->inv_states[(int)msa_compl_char(mod->rate_matrix->
                                                         states[b3_j])];
      j_comp = b3_comp_j * alph_size * alph_size + 
        b2_comp_j * alph_size + b1_comp_j;
                                /* state idx of compl dinuc */

      if ((b1_i != b1_j && b2_i != b2_j) || (b1_i != b1_j && b3_i != b3_j) || 
          (b2_i != b2_j && b3_i != b3_j)) 
        continue;

      val = mm_get(mod->rate_matrix, i, j);
      vec_set(params, start_idx++, val);
   }
    done_row[i] = done_row[i_comp] = 1;
  }    
  if (start_idx != params->size)
    die("ERROR tm_init_mat_from_model_U3S: start_idx (%i) != params->size (%i)\n",
	start_idx, params->size);
}




int tm_flag_subst_param_pos(TreeModel *mod, int *flag, 
			    String *param_name) {
  int numpar, i;
  if (str_equals_nocase_charstr(param_name, RATEMAT_STR)) {
    numpar = tm_get_nratematparams(mod);
    if (numpar==0) return 0;

    for (i=0; i<numpar; i++)
      flag[i] = 1;
    return 1;
  }
  switch(mod->subst_mod) {
  case JC69:  //no named params: return error unless param_name is NULL
  case F81:
    return (param_name == NULL);
  case K80:
  case HKY85:
    if (str_equals_nocase_charstr(param_name, "kappa")) {
      flag[0] = 1;
      return 1;
    }
    return 0;
  case HKY85G:
    if (str_equals_nocase_charstr(param_name, "kappa")) {
      flag[0] = 1;
      return 1;
    }
    else if (str_equals_nocase_charstr(param_name, "gap_param")) {
      flag[1] = 1;
      return 1;
    }
    return 0;
  case REV:
  case  R2:
  case R2S:
  case R3:
  case R3S:
    return 0;
  case REV_GC:
    if (str_equals_nocase_charstr(param_name, "gc_param")) {
      flag[0] = 1;
      return 1;
    }
    else if (str_equals_nocase_charstr(param_name, "alphas")) {
      numpar = mod->rate_matrix->size*(mod->rate_matrix->size-1)/2;
      for (i=1; i<=numpar; i++)
	flag[i] = 1;
      return 1;
    }
    return 0;
  case UNREST:
  case U2:
  case U2S:
  case U3:
  case U3S:
    return 0;
  case GC:
    if (str_equals_nocase_charstr(param_name, "kappa")) {
      flag[0] = 1;
      return 1;
    }
    if (str_equals_nocase_charstr(param_name, "gc_param")) {
      flag[1] = 1;
      return 1;
    }
    return 0;
  default:
    return 0;
  }
}
