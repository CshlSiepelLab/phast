/* $Id: subst_mods.c,v 1.3 2004-07-26 15:02:42 acs Exp $
   Written by Adam Siepel, 2002-2004
   Copyright 2002-2004, Adam Siepel, University of California */

/* Handling of specific substitution models.  This needs reworking:
   was originally set up for small number of models but has become
   unwieldy... */

#include <subst_mods.h>
#include <tree_model.h>
#include <stringsplus.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>

/* internal functions (model-specific) */
void tm_set_JC69_matrix(TreeModel *mod);
void tm_set_K80_matrix(TreeModel *mod, double kappa);
void tm_set_HKY_matrix(TreeModel *mod, double kappa, int kappa_idx);
void tm_set_REV_matrix(TreeModel *mod, gsl_vector *params, int start_idx);
void tm_set_UNREST_matrix(TreeModel *mod, gsl_vector *params, int start_idx);
void tm_set_R2_matrix(TreeModel *mod, gsl_vector *params, int start_idx);
void tm_set_U2_matrix(TreeModel *mod, gsl_vector *params, int start_idx);
void tm_set_R2S_matrix(TreeModel *mod, gsl_vector *params, int start_idx);
void tm_set_U2S_matrix(TreeModel *mod, gsl_vector *params, 
                       int start_idx);
void tm_set_R3_matrix(TreeModel *mod, gsl_vector *params, int start_idx);
void tm_set_R3S_matrix(TreeModel *mod, gsl_vector *params, int start_idx);
void tm_set_U3_matrix(TreeModel *mod, gsl_vector *params, int start_idx);
void tm_set_U3S_matrix(TreeModel *mod, gsl_vector *params, int start_idx);
void tm_init_mat_REV(TreeModel *mod, gsl_vector *params, int nbranches, 
                     double kappa);
void tm_init_mat_UNREST(TreeModel *mod, gsl_vector *params, int nbranches, 
                        double kappa);
void tm_init_mat_R2(TreeModel *mod, gsl_vector *params, 
                    int start_idx, double kappa);
void tm_init_mat_U2(TreeModel *mod, gsl_vector *params, 
                    int start_idx, double kappa);
void tm_init_mat_R2S(TreeModel *mod, gsl_vector *params, 
                     int start_idx, double kappa);
void tm_init_mat_U2S(TreeModel *mod, gsl_vector *params, 
                     int start_idx, double kappa);
void tm_init_mat_R3(TreeModel *mod, gsl_vector *params, 
                    int start_idx, double kappa);
void tm_init_mat_R3S(TreeModel *mod, gsl_vector *params, 
                     int start_idx, double kappa);
void tm_init_mat_U3(TreeModel *mod, gsl_vector *params, 
                    int start_idx, double kappa);
void tm_init_mat_U3S(TreeModel *mod, gsl_vector *params, 
                     int start_idx, double kappa);
void tm_init_mat_from_model_REV(TreeModel *mod, gsl_vector *params, 
                                int start_idx);
void tm_init_mat_from_model_UNREST(TreeModel *mod, gsl_vector *params, 
                                   int start_idx);
void tm_init_mat_from_model_R2(TreeModel *mod, gsl_vector *params, 
                               int start_idx);
void tm_init_mat_from_model_U2(TreeModel *mod, gsl_vector *params, 
                               int start_idx);
void tm_init_mat_from_model_R2S(TreeModel *mod, gsl_vector *params, 
                                int start_idx);
void tm_init_mat_from_model_U2S(TreeModel *mod, gsl_vector *params, 
                                int start_idx);
void tm_init_mat_from_model_R3(TreeModel *mod, gsl_vector *params, 
                               int start_idx);
void tm_init_mat_from_model_R3S(TreeModel *mod, gsl_vector *params, 
                                int start_idx);
void tm_init_mat_from_model_U3(TreeModel *mod, gsl_vector *params, 
                               int start_idx);
void tm_init_mat_from_model_U3S(TreeModel *mod, gsl_vector *params, 
                                int start_idx);

/* Return the substitution model (enum val) corresponding to the
   specified string */
subst_mod_type tm_get_subst_mod_type(char *str) {
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
  else if (str_equals_nocase_charstr(subst_mod_str, "REV"))
    retval = REV;
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
  case REV:
    return "REV";
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
  default:
    return "(unknown model)";
  }
}

/* number of rate matrix params (not counting eq. freqs) */
int tm_get_nratematparams(TreeModel *mod) {
  switch (mod->subst_mod) {
  case JC69:
  case F81:
    return 0;
  case K80:
  case HKY85:
    return 1;
  case REV:
    return (mod->rate_matrix->size * mod->rate_matrix->size 
            - mod->rate_matrix->size) / 2; 
    /* allows use with alternative alphabets, e.g., {ACGT-} */
  case UNREST:
    return 12;
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
  default:
    assert(0);
  }
  return -1;
}

int tm_is_reversible(int subst_mod) {
  return !(subst_mod == UNREST || subst_mod == U2 || subst_mod == U2S || 
           subst_mod == U3 || subst_mod == U3S);
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
void tm_set_rate_matrix(TreeModel *mod, gsl_vector *params, int i) {
  switch(mod->subst_mod) {
  case JC69:
    tm_set_JC69_matrix(mod);
    break;
  case K80:
    tm_set_K80_matrix(mod, gsl_vector_get(params, i));
    break;
  case F81:
    tm_set_HKY_matrix(mod, 1, -1); 
    break;
  case HKY85:
    tm_set_HKY_matrix(mod, gsl_vector_get(params, i), i);
    break;
  case REV:
    tm_set_REV_matrix(mod, params, i);
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
  default:
    assert(0);
  }

  /* NOTE: used to scale rate matrices here, but have removed scaling
     for more well-behaved derivatives.  */
}

/* initialize rate-matrix parameters in parameter vector, using an
   HKY-like strategy (with 'kappa' defining transition/transversion
   bias).  Starting index is 'params_idx'. */
void tm_rate_params_init(TreeModel *mod, gsl_vector *params, 
                         int params_idx, double kappa) {
  switch(mod->subst_mod) {
  case JC69:
  case F81:
    break;                      /* do nothing */
  case K80:
  case HKY85:
    gsl_vector_set(params, params_idx, kappa);
    break;      
  case REV:
    tm_init_mat_REV(mod, params, params_idx, kappa);
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
    assert(0);
  }
}

/* initialize rate-matrix parameters in parameter vector, based on an
   existing model.  Starting index is 'params_idx'. */
void tm_rate_params_init_from_model(TreeModel *mod, gsl_vector *params, 
                                    int params_idx) {
  double kappa;

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

    gsl_vector_set(params, params_idx, kappa);
    break;      
  case REV:
    tm_init_mat_from_model_REV(mod, params, params_idx);
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
    assert(0);
  }
}

/* functions to compute subst. probability matrices directly (without
   diagonalization of a rate matrix), for some simple models */

void tm_set_probs_JC69(TreeModel *mod, MarkovMatrix *P, double t) {
  int i, j;
  double scale = 4.0/3;
  for (i = 0; i < mod->rate_matrix->size; i++) {
    for (j = 0; j < mod->rate_matrix->size; j++) {
      if (i == j)
        mm_set(P, i, j, 0.25 + 0.75 * exp(-t * scale));
      else
        mm_set(P, i, j, 0.25 - 0.25 * exp(-t * scale));
    }
  }
}

void tm_set_probs_F81(TreeModel *mod, MarkovMatrix *P, double scale, double t) {
  int i, j;
  for (i = 0; i < mod->rate_matrix->size; i++) {
    for (j = 0; j < mod->rate_matrix->size; j++) {
      if (i == j)
        mm_set(P, i, j, exp(-t * scale) + 
               gsl_vector_get(mod->backgd_freqs, j) * (1 - exp(-t * scale)));
      else
        mm_set(P, i, j, gsl_vector_get(mod->backgd_freqs, j) * 
               (1 - exp(-t * scale)));

    }
  }
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
      mm_set(mod->rate_matrix, i, j, (float)1/3);
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
      val = gsl_vector_get(mod->backgd_freqs, j); 
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

void tm_set_REV_matrix(TreeModel *mod, gsl_vector *params, int start_idx) {
  int i, j;
  int setup_mapping = (lst_size(mod->rate_matrix_param_row[start_idx]) == 0);
  for (i = 0; i < mod->rate_matrix->size; i++) {
    double rowsum = 0;
    for (j = i+1; j < mod->rate_matrix->size; j++) {
      double val;
      val = gsl_vector_get(params, start_idx);
      mm_set(mod->rate_matrix, i, j, 
             val * gsl_vector_get(mod->backgd_freqs, j));
      mm_set(mod->rate_matrix, j, i, 
             val * gsl_vector_get(mod->backgd_freqs, i));
      rowsum += (val * gsl_vector_get(mod->backgd_freqs, j));

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

void tm_set_UNREST_matrix(TreeModel *mod, gsl_vector *params, int start_idx) {
  int i, j;
  int setup_mapping = (lst_size(mod->rate_matrix_param_row[start_idx]) == 0);
  gsl_matrix_set_zero(mod->rate_matrix->matrix);
  for (i = 0; i < mod->rate_matrix->size; i++) {
    double rowsum = 0;
    for (j = 0; j < mod->rate_matrix->size; j++) {
      double val;
      if (i == j) continue;
      mm_set(mod->rate_matrix, i, j, 
             val = gsl_vector_get(params, start_idx));
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
/*       val = gsl_vector_get(mod->backgd_freqs, j); */
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

void tm_set_R2_matrix(TreeModel *mod, gsl_vector *params, int start_idx) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  int setup_mapping = (lst_size(mod->rate_matrix_param_row[start_idx]) == 0);
  gsl_matrix_set_zero(mod->rate_matrix->matrix);

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

      val = gsl_vector_get(params, start_idx);
      mm_set(mod->rate_matrix, i, j, 
             val * gsl_vector_get(mod->backgd_freqs, j));
      mm_set(mod->rate_matrix, j, i, 
             val * gsl_vector_get(mod->backgd_freqs, i));
      rowsum += (val * gsl_vector_get(mod->backgd_freqs, j));

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
  assert(start_idx == params->size);
}

void tm_set_R2S_matrix(TreeModel *mod, gsl_vector *params, 
                            int start_idx) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  int nstates = mod->rate_matrix->size;
  int done_row[nstates];
  int setup_mapping = (lst_size(mod->rate_matrix_param_row[start_idx]) == 0);
  double rowsum;
  gsl_matrix_set_zero(mod->rate_matrix->matrix);
  
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

      val = gsl_vector_get(params, start_idx);

      /* set rates for i,j and j,i */
      mm_set(mod->rate_matrix, i, j, 
             val * gsl_vector_get(mod->backgd_freqs, j));
      mm_set(mod->rate_matrix, j, i, 
             val * gsl_vector_get(mod->backgd_freqs, i));

      /* set rates for their complements */
      mm_set(mod->rate_matrix, i_comp, j_comp, 
             val * gsl_vector_get(mod->backgd_freqs, j_comp));
      mm_set(mod->rate_matrix, j_comp, i_comp, 
             val * gsl_vector_get(mod->backgd_freqs, i_comp));

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
  
  assert(start_idx == params->size);
}

void tm_set_U2_matrix(TreeModel *mod, gsl_vector *params, int start_idx) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  int setup_mapping = (lst_size(mod->rate_matrix_param_row[start_idx]) == 0);
  gsl_matrix_set_zero(mod->rate_matrix->matrix);
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

      val = gsl_vector_get(params, start_idx);
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
  assert(start_idx == params->size);
}

void tm_set_U2S_matrix(TreeModel *mod, gsl_vector *params, 
                               int start_idx) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  int nstates = mod->rate_matrix->size;
  int done_row[nstates];
  int setup_mapping = (lst_size(mod->rate_matrix_param_row[start_idx]) == 0);
  double rowsum;
  gsl_matrix_set_zero(mod->rate_matrix->matrix);
  
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

      val = gsl_vector_get(params, start_idx);

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
  
  assert(start_idx == params->size);
}

/* FIXME: this should probably be generalized for a reversible matrix
   of arbitrary order */
void tm_set_R3_matrix(TreeModel *mod, gsl_vector *params, int start_idx) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  int setup_mapping = (lst_size(mod->rate_matrix_param_row[start_idx]) == 0);
  gsl_matrix_set_zero(mod->rate_matrix->matrix);

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

      val = gsl_vector_get(params, start_idx);
      mm_set(mod->rate_matrix, i, j, 
             val * gsl_vector_get(mod->backgd_freqs, j));
      mm_set(mod->rate_matrix, j, i, 
             val * gsl_vector_get(mod->backgd_freqs, i));
      rowsum += (val * gsl_vector_get(mod->backgd_freqs, j));

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
  assert(start_idx == params->size);
}

void tm_set_R3S_matrix(TreeModel *mod, gsl_vector *params, int start_idx) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  int nstates = mod->rate_matrix->size;
  int setup_mapping = (lst_size(mod->rate_matrix_param_row[start_idx]) == 0);
  int done_row[nstates];
  double rowsum;
  gsl_matrix_set_zero(mod->rate_matrix->matrix);

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

      val = gsl_vector_get(params, start_idx);
      mm_set(mod->rate_matrix, i, j, 
             val * gsl_vector_get(mod->backgd_freqs, j));
      mm_set(mod->rate_matrix, j, i, 
             val * gsl_vector_get(mod->backgd_freqs, i));
      if (i != j_comp) {
        mm_set(mod->rate_matrix, i_comp, j_comp, 
               val * gsl_vector_get(mod->backgd_freqs, j_comp));
        mm_set(mod->rate_matrix, j_comp, i_comp, 
               val * gsl_vector_get(mod->backgd_freqs, i_comp));
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

  assert(start_idx == params->size);
}

void tm_set_U3_matrix(TreeModel *mod, gsl_vector *params, int start_idx) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  int setup_mapping = (lst_size(mod->rate_matrix_param_row[start_idx]) == 0);
  gsl_matrix_set_zero(mod->rate_matrix->matrix);

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

      val = gsl_vector_get(params, start_idx);
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
  assert(start_idx == params->size);
}

void tm_set_U3S_matrix(TreeModel *mod, gsl_vector *params, 
                               int start_idx) {
  int i, j;
  int alph_size = strlen(mod->rate_matrix->states);
  int nstates = mod->rate_matrix->size;
  int done_row[nstates];
  int setup_mapping = (lst_size(mod->rate_matrix_param_row[start_idx]) == 0);
  double rowsum;
  gsl_matrix_set_zero(mod->rate_matrix->matrix);
  
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

      val = gsl_vector_get(params, start_idx);

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
  
  assert(start_idx == params->size);
}

/* void tm_set_GY_matrix(TreeModel *mod, double kappa, double omega) { */
/*   int i, j; */
/*   int alph_size = strlen(mod->rate_matrix->states); */
/*   int setup_mapping = (lst_size(mod->rate_matrix_param_row[start_idx]) == 0); */
/*   char *aa = get_codon_mapping(mod->rate_matrix->states); */

/*   gsl_matrix_set_zero(mod->rate_matrix->matrix); */

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

/*       val = gsl_vector_get(params, start_idx); */
/*       mm_set(mod->rate_matrix, i, j,  */
/*              val * gsl_vector_get(mod->backgd_freqs, j)); */
/*       mm_set(mod->rate_matrix, j, i,  */
/*              val * gsl_vector_get(mod->backgd_freqs, i)); */
/*       rowsum += (val * gsl_vector_get(mod->backgd_freqs, j)); */

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
/*   assert(start_idx == params->size); */
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
void tm_init_mat_REV(TreeModel *mod, gsl_vector *params, int parm_idx, 
                     double kappa) {
  int i, j;
  for (i = 0; i < mod->rate_matrix->size; i++) {
    for (j = i+1; j < mod->rate_matrix->size; j++) {
      double val = 1;
      if (is_transition(mod->rate_matrix->states[i], 
                        mod->rate_matrix->states[j])) 
        val *= kappa;
      gsl_vector_set(params, parm_idx++, val);
    }
  }    
}

/* initialize UNREST as if HKY */
void tm_init_mat_UNREST(TreeModel *mod, gsl_vector *params, int parm_idx, 
                     double kappa) {
  int i, j;
  for (i = 0; i < mod->rate_matrix->size; i++) {
    for (j = 0; j < mod->rate_matrix->size; j++) {
      double val = 1;
      if (i == j) continue;
      if (is_transition(mod->rate_matrix->states[i], 
                           mod->rate_matrix->states[j])) 
        val *= kappa;
      gsl_vector_set(params, parm_idx++, val + .05 * rand()/(RAND_MAX + 1.0));
                                /* add a little noise to initial
                                   values to avoid making matrix
                                   undiagonalizable :) */
    }
  }    
}

/* initialize R2 as if HKY2 */
void tm_init_mat_R2(TreeModel *mod, gsl_vector *params, int parm_idx, 
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
      gsl_vector_set(params, parm_idx++, val);
    }
  }    
  assert(parm_idx == params->size);
}

/* initialize R2S as if HKY2 */
void tm_init_mat_R2S(TreeModel *mod, gsl_vector *params, 
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

      gsl_vector_set(params, parm_idx++, val);
    }        
    done_row[i] = done_row[i_comp] = 1;
  }
  assert(parm_idx == params->size);
}

/* initialize U2 as if HKY2 */
void tm_init_mat_U2(TreeModel *mod, gsl_vector *params, 
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
      gsl_vector_set(params, parm_idx++, val + .05 * rand()/(RAND_MAX + 1.0));
                                /* add a little noise to initial
                                   values to avoid making matrix
                                   undiagonalizable :) */
    }
  }    
  assert(parm_idx == params->size);
}

/* initialize U2S as if HKY2 */
void tm_init_mat_U2S(TreeModel *mod, gsl_vector *params, 
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

      gsl_vector_set(params, parm_idx++, val + .05 * rand()/(RAND_MAX + 1.0));
                                /* add a little noise to initial
                                   values to avoid making matrix
                                   undiagonalizable :) */
    }        
    done_row[i] = done_row[i_comp] = 1;
  }
  assert(parm_idx == params->size);
}

void tm_init_mat_R3(TreeModel *mod, gsl_vector *params, int parm_idx, 
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
      gsl_vector_set(params, parm_idx++, val);
    }
  }    
  assert(parm_idx == params->size);
}

void tm_init_mat_R3S(TreeModel *mod, gsl_vector *params, int parm_idx, 
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
      gsl_vector_set(params, parm_idx++, val);
    }
    done_row[i] = done_row[i_comp] = 1;
  }    
  assert(parm_idx == params->size);
}

void tm_init_mat_U3(TreeModel *mod, gsl_vector *params, 
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
      gsl_vector_set(params, parm_idx++, val + .001 * rand()/(RAND_MAX + 1.0));
                                /* add a little noise to initial
                                   values to avoid making matrix
                                   undiagonalizable :) */
    }
  }    
  assert(parm_idx == params->size);
}

void tm_init_mat_U3S(TreeModel *mod, gsl_vector *params, 
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
      gsl_vector_set(params, parm_idx++, val + .001 * rand()/(RAND_MAX + 1.0));
                                /* add a little noise to initial
                                   values to avoid making matrix
                                   undiagonalizable :) */
    }
    done_row[i] = done_row[i_comp] = 1;
  }    
  assert(parm_idx == params->size);
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

void tm_init_mat_from_model_REV(TreeModel *mod, gsl_vector *params, 
                                int start_idx) {
  int i, j;
  for (i = 0; i < mod->rate_matrix->size; i++) {
    for (j = i+1; j < mod->rate_matrix->size; j++) {
      double val = safediv(mm_get(mod->rate_matrix, i, j),
                           gsl_vector_get(mod->backgd_freqs, j));
      gsl_vector_set(params, start_idx++, val);
    }
  }
}

void tm_init_mat_from_model_UNREST(TreeModel *mod, gsl_vector *params, 
                                   int start_idx) {
  int i, j;
  for (i = 0; i < mod->rate_matrix->size; i++) {
    for (j = 0; j < mod->rate_matrix->size; j++) {
      double val;
      if (i == j) continue;
      val = mm_get(mod->rate_matrix, i, j);
      gsl_vector_set(params, start_idx++, val);
    }
  }
}

void tm_init_mat_from_model_R2(TreeModel *mod, gsl_vector *params, 
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
                    gsl_vector_get(mod->backgd_freqs, j));
      gsl_vector_set(params, start_idx++, val);
    }
  }
  assert(start_idx == params->size);
}

void tm_init_mat_from_model_U2(TreeModel *mod, gsl_vector *params, 
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
      gsl_vector_set(params, start_idx++, val);
    }
  }
  assert(start_idx == params->size);
}

void tm_init_mat_from_model_R2S(TreeModel *mod, gsl_vector *params, 
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
                    gsl_vector_get(mod->backgd_freqs, j));
      gsl_vector_set(params, start_idx++, val);
    }        
    done_row[i] = done_row[i_comp] = 1;
  }
  assert(start_idx == params->size);
}

void tm_init_mat_from_model_U2S(TreeModel *mod, gsl_vector *params, 
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
      gsl_vector_set(params, start_idx++, val);
    }        
    done_row[i] = done_row[i_comp] = 1;
  }
  assert(start_idx == params->size);
}

void tm_init_mat_from_model_R3(TreeModel *mod, gsl_vector *params, 
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
                    gsl_vector_get(mod->backgd_freqs, j));
      gsl_vector_set(params, start_idx++, val);
    }
  }
  assert(start_idx == params->size);
}

void tm_init_mat_from_model_R3S(TreeModel *mod, gsl_vector *params, 
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
                    gsl_vector_get(mod->backgd_freqs, j));
      gsl_vector_set(params, start_idx++, val);
    }
    done_row[i] = done_row[i_comp] = 1;
  }
  assert(start_idx == params->size);
}

void tm_init_mat_from_model_U3(TreeModel *mod, gsl_vector *params, 
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
      gsl_vector_set(params, start_idx++, val);
    }
  }
  assert(start_idx == params->size);
}

void tm_init_mat_from_model_U3S(TreeModel *mod, gsl_vector *params, 
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
      gsl_vector_set(params, start_idx++, val);
   }
    done_row[i] = done_row[i_comp] = 1;
  }    
  assert(start_idx == params->size);
}

