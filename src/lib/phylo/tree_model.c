/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2009 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: tree_model.c,v 1.41.2.4 2009-03-20 21:25:47 mt269 Exp $ */

#include <tree_model.h>
#include <subst_mods.h>
#include <stacks.h>
#include <stringsplus.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <numerical_opt.h>
#include <tree_likelihoods.h>
#include <time.h>
#include <sys/time.h>
#include <sufficient_stats.h>
#include <matrix.h>
#include <sys/types.h>
#include <unistd.h>
#include <dgamma.h>
#include <math.h>
#include <pthr.h>

#define ALPHABET_TAG "ALPHABET:"
#define BACKGROUND_TAG "BACKGROUND:"
#define ORDER_TAG "ORDER:"
#define SUBST_MOD_TAG "SUBST_MOD:"
#define RATE_MATRIX_TAG "RATE_MAT:"
#define TREE_TAG "TREE:"
#define EVEC_TAG "EVEC_MAT:"
#define EVAL_TAG "EVAL_MAT:"
#define INV_EVEC_TAG "INV_EVEC_MAT:"
#define LNL_TAG "TRAINING_LNL:"
#define ALPHA_TAG "ALPHA:"
#define NRATECATS_TAG "NRATECATS:"
#define RATE_CONSTS_TAG "RATE_CONSTS:"
#define RATE_WEIGHTS_TAG "RATE_WEIGHTS:"
#define ALT_MODEL_TAG "ALT_SUBST_MOD:"


/* internal functions */
double tm_likelihood_wrapper(Vector *params, void *data);

/* tree == NULL implies weight matrix (most other params ignored in
   this case) */
TreeModel *tm_new(TreeNode *tree, MarkovMatrix *rate_matrix, 
                  Vector *backgd_freqs, subst_mod_type subst_mod, 
                  char *alphabet, int nratecats, double alpha, 
                  List *rate_consts, int root_leaf_id) {
  TreeModel *tm = (TreeModel*)smalloc(sizeof(TreeModel));
  int i, j;

  tm->subst_mod = subst_mod;
  tm->order = tm_order(subst_mod);

  tm->backgd_freqs = backgd_freqs;
  tm->tree = tree;
  tm->alpha = alpha;
  tm->root_leaf_id = root_leaf_id;
  tm->empirical_rates = 0;      /* default (see below) */

  if (tree == NULL) {           /* weight matrix */
    tm->rate_matrix = rate_matrix != NULL ? rate_matrix :
      mm_new(int_pow(strlen(alphabet), tm->order+1), alphabet, CONTINUOUS);
                                /* dummy rate matrix; for now, this is
                                   the only way to store the
                                   alphabet */
    tm->P = NULL;
    tm->rK = tm->freqK = NULL;
    tm->nratecats = 1;
  }
  else {                        /* full tree model */
    tm->rate_matrix = rate_matrix != NULL ? rate_matrix :
      mm_new(int_pow(strlen(alphabet), tm->order+1), alphabet, CONTINUOUS);

    if (tm_is_reversible(subst_mod))
      mm_set_eigentype(tm->rate_matrix, REAL);
    /* can assume real eigenvalues and eigenvectors in this case */

    /* set up probability matrices and rate variation stuff */
    tm->nratecats = nratecats;
    tm->P = smalloc(tree->nnodes * sizeof(MarkovMatrix**));
    for (i = 0; i < tree->nnodes; i++) {
      tm->P[i] = smalloc(nratecats * sizeof(MarkovMatrix*));
      for (j = 0; j < nratecats; j++) tm->P[i][j] = NULL;
    }

    tm->rK = smalloc(nratecats * sizeof(double));
    tm->freqK = smalloc(nratecats * sizeof(double));

    if (rate_consts != NULL) {  /* empirical rate model */
      double interval_size, initalpha = (alpha > 0 ? alpha : 1);
      tm->empirical_rates = 1;
      if (nratecats != lst_size(rate_consts) )
        die("ERROR: number of explicitly defined rate constants must equal number of rate categories.\n");
      for (i = 0; i < nratecats; i++) {
        tm->rK[i] = lst_get_dbl(rate_consts, i);
        interval_size = tm->rK[i] - (i > 0 ? tm->rK[i-1] : 0);
        tm->freqK[i] = gamma_pdf(tm->rK[i], initalpha, 1/initalpha) * 
          interval_size; 
        /* init to approx gamma with shape param alpha. */
      }
      normalize_probs(tm->freqK, tm->nratecats);
    }
    else {                      /* dgamma or no rate variation; if
                                   dgamma, rates will be filled out
                                   as needed, according to alpha */
      tm->rK[0] = tm->freqK[0] = 1;
      for (i = 1; i < nratecats; i++) tm->rK[i] = tm->freqK[i] = 0;
    }
  }

  /* various attributes used when fitting a model to an alignment */
  tm->msa = NULL;
  tm->msa_seq_idx = NULL;
  tm->lnL = NULL_LOG_LIKELIHOOD;
  tm->tree_posteriors = NULL;
  tm->use_conditionals = 0;
  tm->category = -1;
  tm->allow_but_penalize_gaps = 0;
  tm->allow_gaps = 1;
  tm->inform_reqd = FALSE;
  tm->estimate_backgd = 0;
  tm->estimate_branchlens = TM_BRANCHLENS_ALL;
  tm->scale = 1;
  tm->subtree_root = NULL;
  tm->scale_sub = 1;
  tm->scale_sub_bound = NB;
  tm->in_subtree = NULL;
  tm->estimate_ratemat = TRUE;
  tm->ignore_branch = NULL;
  tm->alt_subst_mods = NULL;
  tm->alt_subst_mods_node = NULL;
  tm->all_params = NULL;
  tm->param_map = NULL;
  tm->scale_idx = tm->bl_idx = tm->backgd_idx = 
    tm->ratevar_idx = tm->ratematrix_idx = -1;
  tm->rate_matrix_param_row = tm->rate_matrix_param_col = NULL;
  tm->noopt_str = NULL;
  tm->eqfreq_sym = 0;
  return tm;
}

/** Re-initialize a tree model with an altered substitution model,
    number of rate categories, alpha, set of rate consts, or set of
    rate weights. */ 
void tm_reinit(TreeModel *tm,   /**< TreeModel object to reinitialize  */
               subst_mod_type new_subst_mod, 
                                /**< New substitution model */
               int new_nratecats, 
                                /**< New number of rate categories */
               double new_alpha, 
                                /**< New alpha parameter (ignored if
                                   new_nratecats == 1) */
               List *new_rate_consts, 
                                /**< New list of explicit rate
                                   constants.  If non-NULL, implies
                                   use of empirical rates. */
               List *new_rate_weights
                                /**< New list of explicit rate weights
                                   (mixing proportions).  Will be
                                   normalized automatically.  If
                                   new_rate_consts != NULL and
                                   new_rate_weights == NULL, weights
                                   will be initialized to approx gamma
                                   distrib. defined by alpha */
               ) {
  int i, j;
  int old_nratecats = tm->nratecats;
  assert(new_nratecats >= 1);
  tm_free_rmp(tm);
  tm->rate_matrix_param_row = tm->rate_matrix_param_col = NULL;
  if (tm_order(new_subst_mod) != tm->order)
    die("ERROR: Cannot initialize tree model from model of different order.\n");
  tm->subst_mod = new_subst_mod;
  tm->nratecats = new_nratecats;
  tm->alpha = (tm->nratecats > 1 ? new_alpha : 0);
  tm->empirical_rates = new_rate_consts != NULL ? 1 : 0;
  tm->rK = srealloc(tm->rK, new_nratecats * sizeof(double));
  tm->freqK = srealloc(tm->freqK, new_nratecats * sizeof(double));

  if (new_rate_consts != NULL) {  /* empirical rate model */
    double interval_size, initalpha = (new_alpha > 0 ? new_alpha : 1);
    if (new_nratecats != lst_size(new_rate_consts) ||
        (new_rate_weights != NULL && new_nratecats != lst_size(new_rate_weights)))
      die("ERROR: number of explicitly defined rate constants and/or rate weights must equal number of rate categories.\n");
    for (i = 0; i < new_nratecats; i++) {
      tm->rK[i] = lst_get_dbl(new_rate_consts, i);
      if (new_rate_weights != NULL)
        tm->freqK[i] = lst_get_dbl(new_rate_weights, i);
      else {
        interval_size = tm->rK[i] - (i > 0 ? tm->rK[i-1] : 0);
        tm->freqK[i] = gamma_pdf(tm->rK[i], initalpha, 1/initalpha) * 
          interval_size; 
        /* init to approx gamma with shape param alpha. */
      }
    }
    normalize_probs(tm->freqK, tm->nratecats);
  }
  else {                        /* dgamma or no rate variation */
    tm->rK[0] = tm->freqK[0] = 1;
    for (i = 1; i < tm->nratecats; i++) tm->rK[i] = tm->freqK[i] = 0;
  }

  for (i = 0; i < tm->tree->nnodes; i++) {
    tm->P[i] = srealloc(tm->P[i], new_nratecats * sizeof(MarkovMatrix*));
    for (j = old_nratecats; j < new_nratecats; j++) tm->P[i][j] = NULL;
  }

  if (tm_is_reversible(new_subst_mod) && tm->rate_matrix->eigentype == COMPLEX)
    mm_set_eigentype(tm->rate_matrix, REAL);
  else if (!tm_is_reversible(new_subst_mod) && tm->rate_matrix->eigentype == REAL)
    mm_set_eigentype(tm->rate_matrix, COMPLEX);    
}

void tm_free(TreeModel *tm) {
  int i, j;
  if (tm->tree != NULL) {
    if (tm->rate_matrix_param_row != NULL)
      tm_free_rmp(tm);
    for (i = 0; i < tm->tree->nnodes; i++) {
      for (j = 0; j < tm->nratecats; j++)
        if (tm->P[i][j] != NULL) mm_free(tm->P[i][j]);
      free(tm->P[i]);
    }
    if (tm->msa_seq_idx != NULL) free(tm->msa_seq_idx);
    free(tm->P);
    free(tm->rK);
    free(tm->freqK);
    tr_free(tm->tree);
  }
  if (tm->rate_matrix != NULL) mm_free(tm->rate_matrix);
  if (tm->backgd_freqs != NULL) vec_free(tm->backgd_freqs);
  if (tm->ignore_branch != NULL) free(tm->ignore_branch);
  if (tm->in_subtree != NULL) free(tm->in_subtree);
  if (tm->alt_subst_mods != NULL) {
    for (i=0; i<lst_size(tm->alt_subst_mods); i++) 
      tm_free_alt_subst_mod(lst_get_ptr(tm->alt_subst_mods, i));
    free(tm->alt_subst_mods);
  }
  if (tm->alt_subst_mods_node != NULL)
    free(tm->alt_subst_mods_node);
  if (tm->param_map != NULL) free(tm->param_map);
  if (tm->all_params != NULL) vec_free(tm->all_params);
  free(tm);
}

/* set up lists that allow each rate matrix parameter to be mapped to
   the rows and columns in which it appears */
void tm_init_rmp(TreeModel *tm) {
  int i, j;
  if (tm->rate_matrix_param_row != NULL) return;
  if (tm->rate_matrix == NULL || tm->subst_mod == UNDEF_MOD || 
      tm->tree == NULL) 
    tm->rate_matrix_param_row = tm->rate_matrix_param_col = NULL;
  else {
    int nparams = tm_get_nparams(tm);
    int nrmparams = tm_get_nratematparams(tm);
    tm->rate_matrix_param_row = (List**)smalloc(nparams * sizeof(List*));
    tm->rate_matrix_param_col = (List**)smalloc(nparams * sizeof(List*));
    for (i = 0; i < nparams; i++) 
      tm->rate_matrix_param_row[i] = tm->rate_matrix_param_col[i] = NULL;
                                /* NULL for non-rate-matrix params;
                                   keeps indexing simple */
    for (i=0; i<nrmparams; i++) {
      tm->rate_matrix_param_row[tm->ratematrix_idx + i] = lst_new_int(2);
      tm->rate_matrix_param_col[tm->ratematrix_idx + i] = lst_new_int(2);
    }

    if (tm->alt_subst_mods != NULL) {
      subst_mod_type tempmod = tm->subst_mod;
      for (j=0; j<lst_size(tm->alt_subst_mods); j++) {
	AltSubstMod *altmod = lst_get_ptr(tm->alt_subst_mods, j);
	tm->subst_mod = altmod->subst_mod;
	nrmparams = tm_get_nratematparams(tm);
	for (i=0; i<nrmparams; i++) {
	  tm->rate_matrix_param_row[altmod->ratematrix_idx+i] = lst_new_int(2);
	  tm->rate_matrix_param_col[altmod->ratematrix_idx+i] = lst_new_int(2);
	}
      }
      tm->subst_mod = tempmod;
    }
  }
}

void tm_free_rmp(TreeModel *tm) {
  int i;
  int nparams = tm_get_nparams(tm);
  if (tm->rate_matrix_param_row == NULL) return;
  for (i = 0; i < nparams; i++) {
    if (tm->rate_matrix_param_row[i] != NULL) {
      lst_free(tm->rate_matrix_param_row[i]);
      lst_free(tm->rate_matrix_param_col[i]);
    }
  }  
  free(tm->rate_matrix_param_row);
  free(tm->rate_matrix_param_col);
  tm->rate_matrix_param_row = NULL;
  tm->rate_matrix_param_col = NULL;
}

/* read tree model from file, with simple format specifying background
   rates, rate matrix, and tree (topology and branch lengths) */
TreeModel *tm_new_from_file(FILE *f) {
  char tag[STR_MED_LEN], alphabet[MAX_ALPH_SIZE]; 
  String *tmpstr = str_new(STR_LONG_LEN);
  Vector *backgd = NULL, *rate_weights = NULL;
  Matrix *rmat = NULL;
  MarkovMatrix *M = NULL;
  TreeNode *tree = NULL;
  int size = 0, order = 0, nratecats = -1, empty = TRUE;
  double alpha = 0;
  TreeModel *retval;
  subst_mod_type subst_mod = UNDEF_MOD;
  int i, j;
  List *rate_consts = NULL;

  while (fscanf(f, "%s", tag) != EOF) {
    empty = FALSE;
    if (!strcmp(tag, ALPHABET_TAG)) {
      str_readline(tmpstr, f);

      j = 0;
      for (i = 0; i < tmpstr->length; i++) 
        if (!isspace(tmpstr->chars[i])) alphabet[j++] = tmpstr->chars[i];
      alphabet[j] = '\0';
      size = j;
    }
    else if (!strcmp(tag, ORDER_TAG)) {
      str_readline(tmpstr, f);
      if (str_as_int(tmpstr, &order) == 1 || order < 0) {
        fprintf(stderr, "ERROR: bad ORDER line in tree model file.\n");
        exit(1);
      }
      size = int_pow(size, order+1);
    }
    else if (!strcmp(tag, SUBST_MOD_TAG)) {
      str_readline(tmpstr, f);
      str_double_trim(tmpstr);
      subst_mod = tm_get_subst_mod_type(tmpstr->chars); 
    }
    else if (!strcmp(tag, NRATECATS_TAG)) {
      str_readline(tmpstr, f);
      if (str_as_int(tmpstr, &nratecats) == 1 || nratecats < 0) {
        fprintf(stderr, "ERROR: bad NRATECATS line in tree model file.\n");
        exit(1);
      }
    }
    else if (!strcmp(tag, ALPHA_TAG)) {
      str_readline(tmpstr, f);
      if (str_as_dbl(tmpstr, &alpha) == 1 || alpha < 0) {
        fprintf(stderr, "ERROR: bad ALPHA line in tree model file.\n");
        exit(1);
      }
    }
    else if (!strcmp(tag, BACKGROUND_TAG)) {
      if (size == 0) {
        fprintf(stderr, "ERROR: ALPHABET line must precede BACKGROUND and RATE_MATRIX in tree model file.\n");
        exit(1);
      }
      backgd = vec_new_from_file(f, size);
    }
    else if (!strcmp(tag, RATE_MATRIX_TAG)) {
      if (size == 0) {
        fprintf(stderr, "ERROR: ALPHABET line must precede BACKGROUND and RATE_MAT in tree model file.\n");
        exit(1);
      }
      rmat = mat_new_from_file(f, size, size);
    }
    else if (!strcmp(tag, TREE_TAG)) {
      str_readline(tmpstr, f);
      str_double_trim(tmpstr);
      if (tmpstr->chars[tmpstr->length-1] == ';') 
	tmpstr->chars[--tmpstr->length] = '\0';
      tree = tr_new_from_string(tmpstr->chars);      
    }
    else if (strcmp(tag, LNL_TAG) == 0) 
      str_readline(tmpstr, f);  /* discard */
    else if (!strcmp(tag, RATE_CONSTS_TAG)) {
      Vector *tmpvect;
      if (nratecats < 0) 
        die("ERROR: NRATECATS must precede RATE_CONSTS in tree model file.\n");
      /* easiest to use vec_read and convert */
      tmpvect = vec_new_from_file(f, nratecats);
      rate_consts = lst_new_dbl(nratecats);
      for (i = 0; i < nratecats; i++) 
        lst_push_dbl(rate_consts, vec_get(tmpvect, i));
      vec_free(tmpvect);                     
    }
    else if (!strcmp(tag, RATE_WEIGHTS_TAG)) {
      Vector *rate_weights;
      if (nratecats < 0) 
        die("ERROR: NRATECATS must precede RATE_WEIGHTS in tree model file.\n");
      rate_weights = vec_new_from_file(f, nratecats);
    }
    else if (!strcmp(tag, ALT_MODEL_TAG)) {
      die("ERROR: reading ALT_MODELS from model file not yet supported.\n");
    }
    else {
      fprintf(stderr, "ERROR: unrecognized tag in model file (\"%s\").\n", 
	      tag);
      exit(1);
    }
  }

  if (empty) die("ERROR: empty tree model file.\n");

  if (rmat != NULL) 
    M = mm_new_from_matrix(rmat, alphabet, CONTINUOUS);
  else {
    if (size == 0) 
      die ("ERROR: ALPHABET line must precede RATE_MAT in tree model file.");
    M = mm_new(size, alphabet, CONTINUOUS);
  }

  if (nratecats < 0) nratecats = 1;

  retval = tm_new(tree, M, backgd, subst_mod, alphabet, nratecats, alpha, 
                  rate_consts, -1);	 

  if (retval->order != order) 
    die("ERROR: substitution model and order are inconsistent in tree model file.\n");

  /* if rate weights are explicitly defined, set freqK; otherwise
     default weights will be set up by tm_new */
  if (rate_weights != NULL) {
    if (!retval->empirical_rates) 
      die("ERROR: RATE_CONSTS required if RATE_WEIGHTS.\n");
    for (i = 0; i < nratecats; i++) 
      retval->freqK[i] = vec_get(rate_weights, i);
    normalize_probs(retval->freqK, retval->nratecats);
  }

  str_free(tmpstr);
  if (rate_consts != NULL) lst_free(rate_consts);
  if (rate_weights != NULL) vec_free(rate_weights);

  return retval;
}


void tm_print_altmodel(FILE *F, AltSubstMod *alt_model, TreeModel *tm) {
  int i;
  Vector *backgd_freqs;
  MarkovMatrix *rate_matrix;

  fprintf(F, "\n%s %s\n", ALT_MODEL_TAG, alt_model->nodename->chars);
  fprintf(F, "%s %s\n", SUBST_MOD_TAG, 
          tm_get_subst_mod_string(alt_model->subst_mod));
  if (alt_model->backgd_freqs == NULL)
    backgd_freqs = tm->backgd_freqs;
  else backgd_freqs = alt_model->backgd_freqs;
  if (alt_model->rate_matrix == NULL)
    rate_matrix = tm->rate_matrix;
  else rate_matrix = alt_model->rate_matrix;

  fprintf(F, "%s ", BACKGROUND_TAG);
  for (i = 0; i < backgd_freqs->size; i++) 
    fprintf(F, "%f ", vec_get(backgd_freqs, i));
  fprintf(F, "\n");

  fprintf(F, "%s\n", RATE_MATRIX_TAG);
  mat_print(rate_matrix->matrix, F);
}


void tm_print(FILE *F, TreeModel *tm) {
  int i;
  fprintf(F, "%s ", ALPHABET_TAG);
  for (i = 0; i < strlen(tm->rate_matrix->states); i++) 
    fprintf(F, "%c ", tm->rate_matrix->states[i]);
  fprintf(F, "\n");

  fprintf(F, "%s %d\n", ORDER_TAG, tm->order);

  fprintf(F, "%s %s\n", SUBST_MOD_TAG, 
          tm_get_subst_mod_string(tm->subst_mod));

  if (tm->nratecats > 1) {
    fprintf(F, "%s %d\n", NRATECATS_TAG, tm->nratecats);
    if (tm->empirical_rates) {
      fprintf(F, "%s ", RATE_CONSTS_TAG);
      for (i = 0; i < tm->nratecats; i++) fprintf(F, "%f ", tm->rK[i]);
      fprintf(F, "\n");
      fprintf(F, "%s ", RATE_WEIGHTS_TAG);
      for (i = 0; i < tm->nratecats; i++) fprintf(F, "%f ", tm->freqK[i]);
      fprintf(F, "\n");
    }
    else
      fprintf(F, "%s %f\n", ALPHA_TAG, tm->alpha);
  }

  if (tm->lnL != NULL_LOG_LIKELIHOOD)
    fprintf(F, "%s %f\n", LNL_TAG, tm->lnL);

  fprintf(F, "%s ", BACKGROUND_TAG);
  for (i = 0; i < tm->backgd_freqs->size; i++) 
    fprintf(F, "%f ", vec_get(tm->backgd_freqs, i));
  fprintf(F, "\n");

  if (tm->rate_matrix != NULL) {
    fprintf(F, "%s\n", RATE_MATRIX_TAG);
    mat_print(tm->rate_matrix->matrix, F);
  }

  if (tm->tree != NULL) {
    fprintf(F, "%s ", TREE_TAG);
    tr_print(F, tm->tree, 1);
  }

  if (tm->alt_subst_mods != NULL) {
    for (i=0; i<lst_size(tm->alt_subst_mods); i++)
      tm_print_altmodel(F, lst_get_ptr(tm->alt_subst_mods, i), tm);
  }
}

TreeModel *tm_create_copy(TreeModel *src) {
  TreeModel *retval;
  Vector *freqs = vec_new(src->backgd_freqs->size);
  int i, j;

  vec_copy(freqs, src->backgd_freqs);
  retval = tm_new(src->tree != NULL ? tr_create_copy(src->tree) : NULL, 
                  mm_create_copy(src->rate_matrix), freqs, 
                  src->subst_mod, NULL, 
                  src->nratecats, src->alpha, NULL, 
                  src->root_leaf_id);
  retval->msa = src->msa;
  retval->lnL = src->lnL;
  retval->use_conditionals = src->use_conditionals;
  retval->category = src->category;
  retval->allow_gaps = src->allow_gaps;
  retval->allow_but_penalize_gaps = src->allow_but_penalize_gaps;
  retval->inform_reqd = src->inform_reqd;
  retval->estimate_backgd = src->estimate_backgd;
  retval->estimate_branchlens = src->estimate_branchlens;
  retval->scale = src->scale;
  
  if (src->empirical_rates) {
    for (i = 0; i < src->nratecats; i++) {
      retval->rK[i] = src->rK[i];
      retval->freqK[i] = src->freqK[i];
    }
    retval->empirical_rates = 1;
  }

  retval->scale_idx = src->scale_idx;
  retval->bl_idx = src->bl_idx;
  retval->ratematrix_idx = src->ratematrix_idx;
  retval->backgd_idx = src->backgd_idx;
  retval->ratevar_idx = src->ratevar_idx;
  if (src->noopt_str != NULL)
    retval->noopt_str = str_new_charstr(src->noopt_str->chars);
  else retval->noopt_str = NULL;
  retval->eqfreq_sym = src->eqfreq_sym;
  
  if (src->all_params != NULL) {
    retval->all_params = vec_create_copy(src->all_params);
    if (src->param_map != NULL) {
      retval->param_map = malloc(src->all_params->size*sizeof(int));
      for (i=0; i<src->all_params->size; i++)
	retval->param_map[i] = src->param_map[i];
    }
    else retval->param_map = NULL;
  }

  /* Copy lineage-specific models */
  if (retval->alt_subst_mods != NULL) {
    AltSubstMod *currmod, *newmod;
    src->alt_subst_mods = lst_new(lst_size(retval->alt_subst_mods), sizeof(AltSubstMod*));
    for (i = 0; i < lst_size(retval->alt_subst_mods); i++) {
      currmod = (AltSubstMod*)lst_get(src->alt_subst_mods, i);
      newmod = malloc(sizeof(AltSubstMod));
      newmod->subst_mod = currmod->subst_mod;
      if (currmod->backgd_freqs != NULL)
	newmod->backgd_freqs = vec_create_copy(currmod->backgd_freqs);
      else newmod->backgd_freqs = NULL;
      if (currmod->rate_matrix != NULL)
	newmod->rate_matrix = mm_create_copy(currmod->rate_matrix);
      else newmod->rate_matrix = NULL;
      newmod->ratematrix_idx = currmod->ratematrix_idx;
      newmod->backgd_idx = currmod->backgd_idx;
      newmod->separate_model = currmod->separate_model;
      if (currmod->param_list != NULL) {
	String *currstr;
	newmod->param_list = lst_new_ptr(lst_size(currmod->param_list));
	for (j = 0; j < lst_size(currmod->param_list); j++) {
	  currstr = str_new_charstr(((String*)lst_get(currmod->param_list, j))->chars);
	  lst_push_ptr(newmod->param_list, currstr);
	}
      }
      else newmod->param_list = NULL;
      lst_push_ptr(retval->alt_subst_mods, (void*)newmod);
    }
  }
  else retval->alt_subst_mods = NULL;
  if (retval->alt_subst_mods_node != NULL) {
    List *traversal;
    TreeNode *n;
    assert(retval->tree != NULL);
    traversal = tr_preorder(retval->tree);
    for (i=0; i<lst_size(traversal); i++) {
      n = lst_get_ptr(traversal, i);
      if (retval->alt_subst_mods_node[n->id] != NULL) {
	/* Need to find the model for this lineage */
	for (j = 0; j<lst_size(src->alt_subst_mods); j++) {
	  if (lst_get_ptr(src->alt_subst_mods, j) == src->alt_subst_mods_node[n->id])
	    break;
	}
	assert(j < lst_size(src->alt_subst_mods));
	src->alt_subst_mods_node[n->id] = lst_get_ptr(retval->alt_subst_mods, j);
      }
      else src->alt_subst_mods_node[n->id] = NULL;
    }
  }
  else retval->alt_subst_mods_node = NULL;


  /* NOTE: ignoring params, tree_posteriors, etc. */

  return retval;
}


//add a lineage-specific model.   altmod_string should be in the
//form node[+]:MODNAME   or node[+]:param1,param2.
//then an alternative model will be set up for subtree of node, 
//including branch above if the + is present.
//string following the colon will be saved in altmod to be parsed later.
//if it is a model name, an entirely separate subst model will be used 
//for this lineage.  If it is a list of parameters, then same subst 
//model will be used, with these parameters estimated separately.
void tm_add_alt_mod(TreeModel *mod, String *altmod_str) {
  AltSubstMod *altmod = malloc(sizeof(AltSubstMod));
  String *nodestr, *modstr;
  List *templst, *traversal;
  TreeNode *n, *tempnode;
  int parentbranch=0, i;
  templst = lst_new_ptr(2);
  str_split(altmod_str, ":", templst);
  nodestr = (String*)lst_get_ptr(templst, 0);
  modstr = (String*) lst_get_ptr(templst, 1);
  if (lst_size(templst) > 2) 
    die("ERROR: can't parse alt-mod argument %s\n", altmod_str->chars);
  lst_free(templst);
  
  //parse nodes
 tm_add_alt_mod_getnode:
  traversal = tr_preorder(mod->tree);
  for (i=0; i<lst_size(traversal); i++) {
    n = lst_get_ptr(traversal, i);
    if (str_equals_nocase_charstr(nodestr, n->name))
      break;
  }
  if (i == lst_size(traversal)) {
    if (parentbranch==0 && nodestr->chars[nodestr->length-1]=='+') {
      parentbranch=1;
      nodestr->chars[--(nodestr->length)] = '\0';
      goto tm_add_alt_mod_getnode;
    }
    die("ERROR: can't find node %s from alt-mod argument\n", nodestr->chars);
  }
  if (parentbranch) 
    nodestr->chars[nodestr->length++] = '+';
  altmod->nodename = nodestr;

  //set up pointers to node
  if (mod->alt_subst_mods_node ==  NULL) {
    mod->alt_subst_mods_node = malloc(mod->tree->nnodes*sizeof(AltSubstMod*));
    for (i=0; i<mod->tree->nnodes; i++)
      mod->alt_subst_mods_node[i] = NULL;
  }
  traversal = tr_preorder(n);
  for (i=0; i<lst_size(traversal); i++) {
    tempnode = lst_get_ptr(traversal, i);
    if (parentbranch == 0 && tempnode==n) continue;
    assert(tempnode->id < mod->tree->nnodes);
    mod->alt_subst_mods_node[tempnode->id] = altmod;
  }

  //add altmod to mod
  if (mod->alt_subst_mods == NULL)
    mod->alt_subst_mods = lst_new_ptr(1);
  lst_push_ptr(mod->alt_subst_mods, altmod);



  //if modstr is empty, then this model is exactly the same as main model, and
  //no new parameters are estimated.  This may be useful to override a 
  //lineage-specific model in part of a subtree
  if (modstr ==  NULL) {
    altmod->subst_mod = mod->subst_mod;
    altmod->separate_model = 0;
    altmod->param_list = NULL;
  }
  else {
    //if modstr is a subst mod, then this is a separate model 
    //(optimize all params separately)
    altmod->subst_mod = tm_get_subst_mod_type(modstr->chars);
    if (altmod->subst_mod != UNDEF_MOD) {
      altmod->param_list = NULL;
      altmod->separate_model = 1;
    }
    else {
      int *boundpos, inparens=0, tempidx, count=0, maxidx, j;
      char tempchar='*';
      String *tempstr;
      //otherwise it is a list of parameters to be optimized.  Use the 
      //same subst model as main model, but optimize these parameters
      //separately.
      altmod->subst_mod = mod->subst_mod;
      altmod->param_list = lst_new_ptr(10);
      
      //since commas are used to separate boundaries and parameters, 
      //need to temporarily change boundary separators so that we
      //can split param_list.  Only want to keep commas that are not
      //within parenthesis.  This is a bit uglier than it probably should be.
      boundpos = malloc(modstr->length*sizeof(int));
      for (i=0; i<modstr->length; i++) {
	boundpos[i]=0;
	if (modstr->chars[i] == '(') {
	  inparens = 1;
	  count=0;
	}
	else if (modstr->chars[i] == ')')
	    inparens = 0;
	else if (modstr->chars[i] == ',' && inparens) {
	  count++;
	  if (count > 1)  {
	    for (i=0; i<modstr->length; i++) if (boundpos[i]) modstr->chars[i]=',';
	    die("ERROR: format error in --alt-mod argument %s\n", modstr);
	  }
	  boundpos[i]=1;
	  modstr->chars[i] = tempchar;
	}
      }
      if (inparens) die("ERROR: format error in --alt-mod argument (unbalanced parentheses)\n");
		 
      str_split(modstr, ",", altmod->param_list);
      altmod->separate_model = 0;

      tempidx=0;
      for (j=0; j<lst_size(altmod->param_list); j++) {
	tempstr = lst_get_ptr(altmod->param_list, j);
	maxidx = tempidx + tempstr->length;
	for (i=tempidx; i<maxidx; i++) {
	  if (boundpos[i]) {
	    assert(tempstr->chars[i-tempidx]==tempchar);
	    tempstr->chars[i-tempidx]=',';
	  }
	}
	tempidx += tempstr->length+1;
      }
      free(boundpos);
    }
    str_free(modstr);
  }
  altmod->rate_matrix = NULL;
  altmod->backgd_freqs = NULL;
  if (mod->rate_matrix_param_row != NULL) {
    tm_free_rmp(mod);
    tm_init_rmp(mod);
  }
}


void tm_set_subst_matrices(TreeModel *tm) {
  int i, j;
  double scaling_const, tmp, branch_scale;
  Vector *backgd_freqs = tm->backgd_freqs;
  subst_mod_type subst_mod = tm->subst_mod;
  MarkovMatrix *rate_matrix = tm->rate_matrix;
  TreeNode *n;

  scaling_const = -1;

  if (tm->estimate_branchlens != TM_SCALE_ONLY) 
    tm->scale = 1;
                                /* be sure scale factor has an effect
                                   only if estimating scale */

  /* if estimating scale factor for subtree, identify branches in subtree */
  if (tm->estimate_branchlens == TM_SCALE_ONLY && tm->subtree_root != NULL &&
      tm->in_subtree == NULL) 
    tm->in_subtree = tr_in_subtree(tm->tree, tm->subtree_root);

  /* need to compute a matrix scaling constant from the equilibrium
     freqs, in this case (see below) */
  if (subst_mod == F81) {   
    for (i = 0, tmp = 0; i < tm->rate_matrix->size; i++)
      tmp += vec_get(tm->backgd_freqs, i) *
        vec_get(tm->backgd_freqs, i);
    scaling_const = 1.0/(1 - tmp);
  }

  for (i = 0; i < tm->tree->nnodes; i++) {
    branch_scale = tm->scale;
    n = lst_get_ptr(tm->tree->nodes, i);

    if (n->parent == NULL) continue;

    /* special case where subst models differ by branch */
    if (tm->alt_subst_mods != NULL) {
      if (tm->alt_subst_mods_node[n->id] != NULL) {
	backgd_freqs = tm->alt_subst_mods_node[n->id]->backgd_freqs;
	if (backgd_freqs==NULL) backgd_freqs = tm->backgd_freqs;
        subst_mod = tm->alt_subst_mods_node[n->id]->subst_mod;
        rate_matrix = tm->alt_subst_mods_node[n->id]->rate_matrix;
	if (rate_matrix==NULL) rate_matrix = tm->rate_matrix;
      }
      else {
        backgd_freqs = tm->backgd_freqs;
        subst_mod = tm->subst_mod;
        rate_matrix = tm->rate_matrix;
      }
      if (subst_mod == F81) {   /* need branch-specific scale */
        for (j = 0, tmp = 0; j < rate_matrix->size; j++)
          tmp += vec_get(backgd_freqs, j) * vec_get(backgd_freqs, j);
        scaling_const = 1.0/(1 - tmp);
      }
    }

    if (tm->estimate_branchlens == TM_SCALE_ONLY && tm->subtree_root != NULL
        && tm->in_subtree[i]) 
      branch_scale *= tm->scale_sub;

    for (j = 0; j < tm->nratecats; j++) {
      if (tm->P[i][j] == NULL)
        tm->P[i][j] = mm_new(rate_matrix->size, rate_matrix->states, DISCRETE);

      if (tm->ignore_branch != NULL && tm->ignore_branch[i])  
                                /* treat as if infinitely long */
        tm_set_probs_independent(tm, tm->P[i][j]);

      /* for simple models, full matrix exponentiation is not necessary */
      else if (subst_mod == JC69)
        tm_set_probs_JC69(tm, tm->P[i][j], 
                          n->dparent * branch_scale * tm->rK[j]);
      else if (subst_mod == F81)
        tm_set_probs_F81(backgd_freqs, tm->P[i][j], scaling_const, 
                         n->dparent * branch_scale * tm->rK[j]);

      else {                     /* full matrix exponentiation */
        mm_exp(tm->P[i][j], rate_matrix, 
               n->dparent * branch_scale * tm->rK[j]);
      }
    }
  }
}

/* version of above that can be used with specified branch length and
   prob matrix */
void tm_set_subst_matrix(TreeModel *tm, MarkovMatrix *P, double t) {
  int i;
  double scaling_const = -1, tmp;

  assert(tm->alt_subst_mods == NULL);
  assert(tm->estimate_branchlens != TM_SCALE_ONLY);

  /* need to compute a matrix scaling constant from the equilibrium
     freqs, in this case (see below) */
  if (tm->subst_mod == F81) {   
    for (i = 0, tmp = 0; i < tm->rate_matrix->size; i++)
      tmp += vec_get(tm->backgd_freqs, i) * vec_get(tm->backgd_freqs, i);
    scaling_const = 1.0/(1 - tmp);
  }

  /* for simple models, full matrix exponentiation is not necessary */
  if (tm->subst_mod == JC69)
    tm_set_probs_JC69(tm, P, t);
  else if (tm->subst_mod == F81)
    tm_set_probs_F81(tm->backgd_freqs, P, scaling_const, t);
  else 
    mm_exp(P, tm->rate_matrix, t);
}

/* scale evolutionary rate by const factor (affects branch lengths
   only) */
void tm_scale(TreeModel *tm, double scale_const, int reset_subst_mats) {
  if (tm->tree == NULL) return;
  tr_scale(tm->tree, scale_const);
  if (reset_subst_mats) tm_set_subst_matrices(tm);
}

/* Generates an alignment according to set of Tree Models and a
   Markov matrix defing how to transition among them.  TreeModels must
   appear in same order as the states of the Markov matrix. 
   NOTE: call srandom externally. */
MSA *tm_generate_msa(int ncolumns, 
                     HMM *hmm,  /* if NULL, single tree model assumed */
                     TreeModel **classmods, 
                     int *labels /* if non-NULL, will be used to
                                    record state (model) responsible
                                    for generating each site; pass
                                    NULL if hmm is NULL */
                     ) {

  int i, class, nseqs, col, ntreenodes, idx;
  MSA *msa;
  Stack *stack;
  char *newchar;
  char **names, **seqs;

  int nclasses = hmm == NULL ? 1 : hmm->nstates;

  /* obtain number of sequences from tree models; ensure all have same
     number */
  ntreenodes = classmods[0]->tree->nnodes; 

  stack = stk_new_ptr(ntreenodes);
  nseqs = -1;
  for (i = 0; i < nclasses; i++) {
    /* count leaves in tree */
    int num = (classmods[i]->tree->nnodes + 1) / 2;

    assert(classmods[i]->nratecats == 1); /* assuming no rate variation */

    if (nseqs == -1) 
      nseqs = num;
    else if (nseqs != num) 
      die("ERROR in tm_generate_msa: model #%d has %d taxa, while a previous model had %d taxa.\n", i+1, num, nseqs);
  }

  /* create new MSA */
  names = (char**)smalloc(nseqs * sizeof(char*));
  seqs = (char**)smalloc(nseqs * sizeof(char*));
  for (i = 0; i < nseqs; i++) 
    seqs[i] = (char*)smalloc((ncolumns + 1) * sizeof(char));
  msa = msa_new(seqs, names, nseqs, ncolumns, 
                classmods[0]->rate_matrix->states);

  /* build sequence idx map; only need one for first model */
  /* FIXME: this assumes all tree models have the same topology; may
     want to relax... */
  classmods[0]->msa_seq_idx = smalloc(classmods[0]->tree->nnodes * 
                                      sizeof(int));
  for (i = 0, idx = 0; i < classmods[0]->tree->nnodes; i++) {
    TreeNode *n = lst_get_ptr(classmods[0]->tree->nodes, i);
    if (n->lchild == NULL && n->rchild == NULL) {
      classmods[0]->msa_seq_idx[i] = idx;
      names[idx] = strdup(n->name);
      idx++;
    }
    else classmods[0]->msa_seq_idx[i] = -1;
  }

  /* generate sequences, column by column */
  if (hmm != NULL && hmm->begin_transitions != NULL)
    class = draw_index(hmm->begin_transitions->data, hmm->nstates);
  else
    class = 0;
  newchar = (char*)smalloc(ntreenodes * sizeof(char));
  for (col = 0; col < ncolumns; col++) {
    List *traversal = tr_preorder(classmods[class]->tree);
    newchar[classmods[class]->tree->id] = 
      mm_sample_backgd(classmods[class]->rate_matrix->states, 
                       classmods[class]->backgd_freqs);
    for (i = 0; i < lst_size(traversal); i++) {
      TreeNode *n = lst_get_ptr(traversal, i);
      TreeNode *l = n->lchild;
      TreeNode *r = n->rchild;
      assert ((l == NULL && r == NULL) || (l != NULL && r != NULL));

      if (l == NULL) 
        msa->seqs[classmods[0]->msa_seq_idx[n->id]][col] = newchar[n->id];
      else {
        MarkovMatrix *lsubst_mat, *rsubst_mat;
        if (classmods[class]->P[l->id][0] == NULL)
          tm_set_subst_matrices(classmods[class]);
        lsubst_mat = classmods[class]->P[l->id][0];
        rsubst_mat = classmods[class]->P[r->id][0];
        newchar[l->id] = mm_sample_char(lsubst_mat, newchar[n->id]);
        newchar[r->id] = mm_sample_char(rsubst_mat, newchar[n->id]);
      }
    }
    if (labels != NULL) labels[col] = class;
    if (hmm != NULL)
      class = mm_sample_state(hmm->transition_matrix, class);
  }

  return msa;
}


/* Assumes *lower_bounds and *upper_bounds have not been allocated.
   Also assumes tm_setup_params has already been called.
   Sets each to NULL if no boundaries are needed */
void tm_set_boundaries(Vector **lower_bounds, Vector **upper_bounds,
		       int npar, TreeModel *mod) {
  int i, j;

  /* most params have lower bound of zero and no upper bound */
  *lower_bounds = vec_new(npar);
  vec_zero(*lower_bounds);
  *upper_bounds = NULL;

  /* however, in this case we don't want the eq freqs to go to zero */
  if (mod->estimate_backgd) {
    for (i = 0; i < mod->backgd_freqs->size; i++) {
      if (mod->param_map[mod->backgd_idx+i] >= 0)
	vec_set(*lower_bounds, mod->param_map[mod->backgd_idx+i], 0.001);
    }
  }
  /* Check eq freqs in lineage-specific models too */
  if (mod->alt_subst_mods != NULL) {
    for (j=0; j<lst_size(mod->alt_subst_mods); j++) {
      AltSubstMod *altmod = lst_get_ptr(mod->alt_subst_mods, j);
      if (altmod->backgd_freqs == NULL) continue;
      if (mod->param_map[altmod->backgd_idx] >= 0)
	for (i=0; i<altmod->backgd_freqs->size; i++)
	  vec_set(*lower_bounds, mod->param_map[altmod->backgd_idx+i], 
		  0.001);
    }
  }

  /* Also, in this case, we need to bound the scale of the subtree */
  if (mod->estimate_branchlens == TM_SCALE_ONLY && mod->subtree_root != NULL 
      && mod->scale_sub_bound != NB) {
    if (mod->scale_sub_bound == LB && mod->param_map[mod->scale_idx+1] >= 0) 
      vec_set(*lower_bounds, mod->param_map[mod->scale_idx+1], 1);
    if (mod->scale_sub_bound == UB && mod->param_map[mod->scale_idx+1] >= 0) {
      *upper_bounds = vec_new(npar);
      vec_set_all(*upper_bounds, INFTY);
      vec_set(*upper_bounds, mod->param_map[mod->scale_idx+1], 1);
    }
  }

}
		       


/* Given an MSA, a tree topology, and a substitution model, fit a tree
   model using a multidimensional optimization algorithm (BFGS).
   TreeModel 'mod' must already be allocated, and initialized with
   desired tree topology, substitution model, and (if appropriate)
   background frequencies.  The vector 'params' should define the
   initial values for the optimization procedure.  Fuction returns 0
   on success, 1 on failure.  */  
int tm_fit(TreeModel *mod, MSA *msa, Vector *params, int cat, 
           opt_precision_type precision, FILE *logf) {
  double ll;
  Vector *lower_bounds, *upper_bounds, *opt_params;
  int i, retval = 0, npar;

  if (msa->ss == NULL) {
    assert(msa->seqs != NULL);
    ss_from_msas(msa, mod->order+1, 0, NULL, NULL, NULL, -1);
  }

  if (mod->backgd_freqs == NULL) { 
    mod->backgd_freqs = vec_new(int_pow(strlen(msa->alphabet), 
					mod->order+1));
    if (mod->subst_mod == JC69 || mod->subst_mod == K80)
      vec_set_all(mod->backgd_freqs, 1.0/mod->backgd_freqs->size);
    else
      msa_get_base_freqs_tuples(msa, mod->backgd_freqs, mod->order + 1, cat);
    for (i=0; i<mod->backgd_freqs->size; i++)
      vec_set(params, mod->backgd_idx+i, vec_get(mod->backgd_freqs, i));
  }
  if (mod->alt_subst_mods != NULL) {
    AltSubstMod *altmod;
    for (i = 0; i < lst_size(mod->alt_subst_mods); i++) {
      altmod = lst_get_ptr(mod->alt_subst_mods, i);
      if (altmod->separate_model == 0 && 
	  altmod->param_list == NULL) continue;

      if (altmod->backgd_freqs == NULL) 
	altmod->backgd_freqs = vec_create_copy(mod->backgd_freqs);
      if (altmod->rate_matrix == NULL)
	altmod->rate_matrix = mm_create_copy(mod->rate_matrix);
    }
  }

  if (mod->tree == NULL) {      /* weight matrix */
    mod->lnL = tl_compute_log_likelihood(mod, msa, NULL, cat, NULL) * 
      log(2);
    return 0;
  }

  mod->msa = msa;               /* package with mod any data needed to
                                   compute likelihoods */
  mod->category = cat;

  npar=0;
  for (i=0; i<params->size; i++) {
    vec_set(mod->all_params, i, vec_get(params, i));
    if (mod->param_map[i] >= npar)
      npar = mod->param_map[i]+1;
  }
  assert(npar > 0);
  opt_params = vec_new(npar);
  for (i=0; i<params->size; i++)
    if (mod->param_map[i] >=0)
      vec_set(opt_params, mod->param_map[i], vec_get(params, i));
  
  tm_set_boundaries(&lower_bounds, &upper_bounds, npar, mod);

  retval = opt_bfgs(tm_likelihood_wrapper, opt_params, (void*)mod, &ll, 
                    lower_bounds, upper_bounds, logf, NULL, precision, NULL);

  mod->lnL = ll * -1 * log(2);  /* make negative again and convert to
                                   natural log scale */
  tm_unpack_params(mod, opt_params, -1);
  vec_copy(params, mod->all_params);
  vec_free(opt_params);

  if (mod->subst_mod != JC69 && mod->subst_mod != F81 && 
      mod->estimate_branchlens != TM_BRANCHLENS_NONE) {   
                                /* in this case, need to scale final
                                   model (JC69 and F81 are exceptions,
                                   because probability subst. matrices
                                   are not derived from the rate
                                   matrix) */
    double scale = tm_scale_rate_matrix(mod);
    tm_scale(mod, scale, 0);    /* branch lengths */
    tm_scale_params(mod, params, scale);
  }

  if (mod->estimate_branchlens == TM_SCALE_ONLY) { 
                                /* also, do final scaling of tree if
                                   estimating scale only */
    tm_scale(mod, mod->scale, 0);
    mod->scale = 1;
    if (mod->subtree_root != NULL) {  /* estimating subtree scale */
      tr_scale_subtree(mod->tree, mod->subtree_root, mod->scale_sub);
      mod->subtree_root->dparent *= mod->scale_sub; /* also need to scale 
                                                       leading branch */
      mod->scale_sub = 1;
    }
  }

  if (lower_bounds != NULL) vec_free(lower_bounds);
  if (upper_bounds != NULL) vec_free(upper_bounds);

  if (retval != 0) 
    fprintf(stderr, "WARNING: BFGS algorithm reached its maximum number of iterations.\n");

  return retval;
}

int void_str_equals_charstr(void *strptr, void *charptr) {
  return str_equals_charstr(*((String**)strptr),
			    (char*)charptr);
}


/* tm_setup_params: assigns parameter indices to mod->scale_idx,
   mod->bl_idx, mod->ratevar_idx, mod->ratematrix_idx which point
   to indices in mod->all_params where corresponding parameters are
   stored.  (as well as pointers in altmod structures).
   
   Then, it determines which parameters are to be optimized, and allocates
   and assigns values to mod->param_map, which maps indices in all_params
   to indices in the optimization vector, so that param_map[i]=j implies
   that all_params[i] is stored in opt_vec[j].  Only parameters which
   are to be optimized are stored in the optimization vector.  Parameters
   which are held constant have param_map[i]=-1.   
 */
void tm_setup_params(TreeModel *mod) {
  int i, opt_idx = 0, next_idx, alph_size, pos, numpar, *flag;
  List *noopt=NULL;

  //first assign indices in mod->all_params to give position of each
  //type of parameter
  mod->scale_idx = 0;
  mod->bl_idx = 2;
  assert(mod->tree != NULL);
  if (mod->backgd_freqs == NULL)
    alph_size = strlen(mod->rate_matrix->states);
  else alph_size = mod->backgd_freqs->size;

  //subtract 1 for root 
  mod->backgd_idx = mod->bl_idx + mod->tree->nnodes - 1;
  next_idx = mod->backgd_idx + alph_size;
  if (mod->nratecats > 1) {
    mod->ratevar_idx = next_idx;
    if (mod->empirical_rates)
      next_idx += mod->nratecats;
    else next_idx++;
  }
  else mod->ratevar_idx = -1;

  mod->ratematrix_idx = next_idx;
  next_idx += tm_get_nratematparams(mod);

  if (mod->alt_subst_mods != NULL) {
    AltSubstMod *altmod;
    subst_mod_type tempmod = mod->subst_mod;
    for (i = 0; i<lst_size(mod->alt_subst_mods); i++) {
      altmod = lst_get_ptr(mod->alt_subst_mods, i);
      mod->subst_mod = altmod->subst_mod;
      altmod->backgd_idx = next_idx;
      next_idx += alph_size;
      altmod->ratematrix_idx = next_idx;
      next_idx += tm_get_nratematparams(mod);
    }
    mod->subst_mod = tempmod;
  }

  //now, allocate all_params and param_map if necessary
  if (mod->all_params != NULL && mod->all_params->size != next_idx) {
    vec_free(mod->all_params);
    if (mod->param_map != NULL) {
      free(mod->param_map);
      mod->param_map = NULL;
    }
  }
  if (mod->all_params == NULL) 
    mod->all_params = vec_new(next_idx);
  if (mod->param_map == NULL)
    mod->param_map = malloc(next_idx*sizeof(int));

  //now assign positions in param_map so that
  //param_map[i] >= 0 if this parameter is optimized
  //param_map[i] = -1 if this parameter is held constant
  for (i=0; i<next_idx; i++)
    mod->param_map[i] = -1;

  if (mod->noopt_str != NULL) {
    noopt = lst_new_ptr(3);
    str_split(mod->noopt_str, ",", noopt);
    pos = lst_find_compare(noopt, BRANCHES_STR, void_str_equals_charstr);
    if (pos >= 0) {
      //      printf("holding branches const\n");
      mod->estimate_branchlens = TM_BRANCHLENS_NONE;
      str_free(lst_get_ptr(noopt, pos));
      lst_delete_idx(noopt, pos);
    }
    pos = lst_find_compare(noopt, BACKGD_STR, void_str_equals_charstr);
    if (pos >= 0) {
      //      printf("holding backgd const\n");
      mod->estimate_backgd = 0;
      str_free(lst_get_ptr(noopt, pos));
      lst_delete_idx(noopt, pos);
    }
    pos = lst_find_compare(noopt, RATEMAT_STR, void_str_equals_charstr);
    if (pos >= 0) {
      //      printf("holding ratematrix const\n");
      mod->estimate_ratemat = 0;
      str_free(lst_get_ptr(noopt, pos));
      lst_delete_idx(noopt, pos);
    }
  }

  //first do scale/branchlength params
  if (mod->estimate_branchlens == TM_SCALE_ONLY) {
    mod->param_map[mod->scale_idx] = opt_idx++;
    if (mod->subtree_root != NULL)
      mod->param_map[mod->scale_idx+1] = opt_idx++;
  }
  else if (mod->estimate_branchlens == TM_BRANCHLENS_ALL) {
    List *traversal = tr_preorder(mod->tree);
    int node_idx, root_idx = -1, temp_idx=0;
    int isreversible = tm_is_reversible(mod->subst_mod);
    TreeNode *n;
    for (node_idx=0; node_idx < lst_size(traversal); node_idx++) {
      n = lst_get_ptr(traversal, node_idx);
      if (n->parent == NULL) continue;  //skip root
      if (mod->root_leaf_id == node_idx) {
	temp_idx++;
	continue;  //don't optimize
      }
      if ((n == mod->tree->lchild || n == mod->tree->rchild) &&
	  isreversible) {
	if (root_idx == -1) root_idx = opt_idx++;
	mod->param_map[mod->bl_idx + temp_idx++] = root_idx;
      }
      else mod->param_map[mod->bl_idx + temp_idx++] = opt_idx++;
    }
  }
  //else TM_BRANCHLENS_NONE: keep everything at -1

  if (mod->estimate_backgd) {
    if (mod->eqfreq_sym==1) {
      int atpos=-1, gcpos=-1;
      char c;
      for (i=0; i<alph_size; i++) {
	c = tolower(mod->rate_matrix->states[i]);
	switch (c) {
	case 'g':
	case 'c':
	  if (gcpos == -1) gcpos = opt_idx++;
	  pos = gcpos;
	  break;
	case 'a':
	case 't':
	  if (atpos == -1) atpos = opt_idx++;
	  pos = atpos;
	  break;
	case '-':
	  pos = opt_idx++;
	  break;
	default:
	  die("ERROR: eqfreq_sym only defined for ACGT alphabet right now (got %c)\n", mod->rate_matrix->states[i]);
	}
	mod->param_map[mod->backgd_idx+i] = pos;
      }
    }
    else {
      for (i=0; i<alph_size; i++)
	mod->param_map[mod->backgd_idx+i] = opt_idx++;
    }
  }

  if (mod->nratecats > 1) {
    int est_rates = (noopt==NULL);
    if (noopt != NULL) {
      pos = lst_find_compare(noopt, RATEVAR_STR, void_str_equals_charstr);
      if (pos >= 0) {
	//	printf("holding ratevar const\n");
	est_rates = 0;
	str_free(lst_get_ptr(noopt, pos));
	lst_delete_idx(noopt, pos);
      }
    }
    if (est_rates) {
      if (mod->empirical_rates) {
	for (i=0; i<mod->nratecats; i++)
	  mod->param_map[mod->ratevar_idx+i] = opt_idx++;
      }
      else {
	mod->param_map[mod->ratevar_idx] = opt_idx++;
      }
    }
  }
  
  numpar = tm_get_nratematparams(mod);
  if (numpar > 0) 
    flag = malloc(numpar*sizeof(int));
  for (i=0; i<numpar; i++) flag[i]=0;
  if (noopt != NULL) {
    for (i=0; i<lst_size(noopt); i++) {
      if (0 == tm_flag_subst_param_pos(mod, flag, lst_get_ptr(noopt, i))) 
	die("ERROR: couldn't parse parameter name %s from no-opt argument for model %s\n", ((String*)lst_get_ptr(noopt, i))->chars, tm_get_subst_mod_string(mod->subst_mod));
      str_free(lst_get_ptr(noopt, i));
    }
    lst_free(noopt);
  }
  
  if (mod->estimate_ratemat) {
    for (i=0; i<numpar; i++)
      if (flag[i] == 0)
	mod->param_map[mod->ratematrix_idx+i] = opt_idx++;
  }
  if (numpar > 0) free(flag);

  if (mod->alt_subst_mods != NULL) {
    AltSubstMod* altmod;
    int j, *opt_par, opt_freq;
    subst_mod_type tempmod = mod->subst_mod;
    for (j=0; j<lst_size(mod->alt_subst_mods); j++) {
      altmod = lst_get_ptr(mod->alt_subst_mods, j);
      opt_freq = mod->estimate_backgd && altmod->separate_model;
      mod->subst_mod = altmod->subst_mod;
      altmod = lst_get_ptr(mod->alt_subst_mods, j);
      numpar = tm_get_nratematparams(mod);
      opt_par = malloc(numpar*sizeof(int));

      //if altmod->separate_model, optimize everything separately
      if (altmod->separate_model == 1) {
	for (i=0; i<numpar; i++)
	  opt_par[i] = 1;
      }
      //if sepatate_model==0 but param_list==NULL, this isn't really 
      //a separate model at all.  Optimize nothing.
      else if (altmod->param_list == NULL) {
	for (i=0; i<numpar; i++)
	  opt_par[i] = 0;
	if (altmod->backgd_freqs != NULL) {
	  vec_free(altmod->backgd_freqs);
	  altmod->backgd_freqs=NULL;
	}
	if (altmod->rate_matrix != NULL) {
	  mm_free(altmod->rate_matrix);
	  altmod->rate_matrix = NULL;
	}
      }
      else {
	//otherwise we only optimize pars in list
	for (i=0; i<numpar; i++)
	  opt_par[i] = 0;
	for (i=0; i<lst_size(altmod->param_list); i++) {
	  String *currparam;
	  int parenpos=-1, k;
	  currparam = lst_get_ptr(altmod->param_list, i);
	  for (k=0; k<currparam->length; k++)
	    if (currparam->chars[k]=='(') {
	      parenpos = k;
	      currparam->chars[parenpos]='\0';
	      currparam->length=parenpos;
	      break;
	    }
	  //	  printf("parsing string %s\n", currparam->chars);
	  if (str_equals_nocase_charstr(currparam, BACKGD_STR))
	    opt_freq = 1;
	  else  {
	    if (0==tm_flag_subst_param_pos(mod, opt_par, currparam)) {
	      die("altmod %s couldn't parse parameter name %s\n",
		  tm_get_subst_mod_string(altmod->subst_mod),
		  currparam->chars);
	    }
	  }
	  if (parenpos >= 0) {
	    currparam->chars[parenpos] = '(';
	    currparam->length = strlen(currparam->chars);
	  }
	}
      }
      
      if (opt_freq) {
	for (i=0; i<alph_size; i++)
	  mod->param_map[altmod->backgd_idx+i] = opt_idx++;
      }
      else {
	for (i=0; i<alph_size; i++)
	  mod->param_map[altmod->backgd_idx+i] = 
	    mod->param_map[mod->backgd_idx+i];
      }
      for (i=0; i<numpar; i++) {
	if (opt_par[i])
	  mod->param_map[altmod->ratematrix_idx+i] = opt_idx++;
	else mod->param_map[altmod->ratematrix_idx+i] = 
	       mod->param_map[mod->ratematrix_idx+i];
      }
      free(opt_par);
    }
    mod->subst_mod = tempmod;
  }
  tm_init_rmp(mod);
}


/* Wrapper for computation of likelihood */
double tm_likelihood_wrapper(Vector *params, void *data) {
  TreeModel *mod = (TreeModel*)data;
  tm_unpack_params(mod, params, -1);
  return -1 * tl_compute_log_likelihood(mod, mod->msa, NULL, mod->category,
                                        NULL);
}

/* Set specified TreeModel according to specified parameter vector
   (exact behavior depends on substitution model).  An index offset
   can be specified for cases in which vectors of parameters are
   nested within larger vectors of parameters (set to -1 if not
   needed) */
void tm_unpack_params(TreeModel *mod, Vector *params_in, int idx_offset) {
  TreeNode *n;
  int nodeidx, i, j;
  List *traversal;
  Vector *params = mod->all_params;

  if (idx_offset == -1) idx_offset = 0;

  for (j=0; j<params->size; j++) {
    if (mod->param_map[j] >= 0)
      vec_set(params, j, vec_get(params_in, mod->param_map[j] + idx_offset));
  }

  mod->scale = vec_get(params, mod->scale_idx);
  mod->scale_sub = vec_get(params, mod->scale_idx+1);
  traversal = tr_preorder(mod->tree);
  i=0;
  for (nodeidx=0; nodeidx < lst_size(traversal); nodeidx++) {
    n = lst_get_ptr(traversal, nodeidx);
    if (n->parent == NULL) continue;
    n->dparent = vec_get(params, mod->bl_idx + i);
    i++;
  }

  /* backgd freqs */
  double sum = 0;
  for (j = 0; j < mod->backgd_freqs->size; j++) 
    sum += vec_get(params, mod->backgd_idx+j);
  for (j = 0; j < mod->backgd_freqs->size; j++)
    vec_set(mod->backgd_freqs, j, vec_get(params, mod->backgd_idx+j)/sum);

  /* rate variation */
  if (mod->nratecats > 1) {
    if (mod->empirical_rates) {
      for (j = 0; j < mod->nratecats; j++)
        mod->freqK[j] = vec_get(params, mod->ratevar_idx+j);
                                /* assume normalized */
    }
    else {                      /* discrete gamma model */
      mod->alpha = vec_get(params, mod->ratevar_idx);
      DiscreteGamma(mod->freqK, mod->rK, mod->alpha, mod->alpha, 
                    mod->nratecats, 0); 
    }
  }

  /* lineage-specific models */ 
  if (mod->alt_subst_mods != NULL) {
    AltSubstMod *altmod;
    subst_mod_type temp_mod = mod->subst_mod;
    MarkovMatrix *temp_mm = mod->rate_matrix;
    Vector *temp_backgd = mod->backgd_freqs;

    for (j = 0; j < lst_size(mod->alt_subst_mods); j++) {
      altmod = lst_get_ptr(mod->alt_subst_mods, j);
      mod->subst_mod = altmod->subst_mod;

      if (altmod->backgd_freqs != NULL) {
	double sum = 0;
	for (i = 0; i < altmod->backgd_freqs->size; i++) 
	  sum += vec_get(params, altmod->backgd_idx+i);
	for (i = 0; i < altmod->backgd_freqs->size; i++)
	  vec_set(altmod->backgd_freqs, i, 
		  vec_get(params, altmod->backgd_idx+i)/sum);
	mod->backgd_freqs = altmod->backgd_freqs;
      }
      if (altmod->rate_matrix != NULL) {
	mod->rate_matrix = altmod->rate_matrix;
	tm_set_rate_matrix(mod, params, altmod->ratematrix_idx);
	if (altmod->subst_mod != JC69 && altmod->subst_mod != F81)
	  mm_diagonalize(altmod->rate_matrix);
      }
    }
    mod->subst_mod = temp_mod;
    mod->rate_matrix = temp_mm;
    mod->backgd_freqs = temp_backgd;
  }

  /* redefine substitution matrix */
  tm_set_rate_matrix(mod, params, mod->ratematrix_idx);
  
  /* diagonalize, if necessary */
  if (mod->subst_mod != JC69 && mod->subst_mod != F81)
    mm_diagonalize(mod->rate_matrix);

  /* set exponentiated version at each edge */
  tm_set_subst_matrices(mod); 
}


/* Scale the rate matrix of the specified tree model such that the
   expected rate of substitution at equilibrium is 1.  This is the
   standard convention; it allows the unit of measurement of branch
   lengths to be expected substitutions/site.  
   NOTE: now returns scaling factor, which can be used to adjust
   branch lengths if the rate matrix has not been scaled throughout the
   fitting procedure.

   NOTE NOTE!: how do we fix this to deal with alt models?  Just leaving
   it unchanged for now.
*/ 
double tm_scale_rate_matrix(TreeModel *mod) {
  double scale = 0;
  int i, j;
  for (i = 0; i < mod->rate_matrix->size; i++) {
    double rowsum = 0;
    for (j = 0; j < mod->rate_matrix->size; j++) {
      if (j != i) rowsum += mm_get(mod->rate_matrix, i, j);
    }
    scale += (vec_get(mod->backgd_freqs, i) * rowsum);
  }

  mm_scale(mod->rate_matrix, (mod->order + 1)/scale);
  return scale/(mod->order + 1);
}


/* scale a parameter vector according to a specified rate matrix scale
   factor.  Branch length params are multiplied by the specified factor
   and rate matrix params by its inverse */
/* Again, I (mt269) am not sure what to do with this in terms of 
   incorporating lineage-specific models.  Leave unchanged for now */
void tm_scale_params(TreeModel *mod, Vector *params, double scale_factor) {
  int i;
  int nbl = tm_get_nbranchlenparams(mod);
  int nrm = tm_get_nratematparams(mod);
  int np = tm_get_nparams(mod);     /* we won't assume anything about the
                                       number of params between the
                                       branchlen and ratemat params */
  for (i = 0; i < nbl; i++)
    vec_set(params, i, vec_get(params, i) * scale_factor);
  for (i = np - nrm; i < np; i++) 
    vec_set(params, i, vec_get(params, i) / scale_factor);
}


void tm_init_rootleaf(TreeModel *mod, Vector *params) {
  List *traversal = tr_preorder(mod->tree);
  int idx=0, i;
  TreeNode *n;
  if (mod->root_leaf_id < 0) return;
  for (i=0; i<lst_size(traversal); i++) {
    n = lst_get_ptr(traversal, i);
    if (n == mod->tree) continue;
    if (n->id == mod->root_leaf_id) {
      vec_set(params, mod->bl_idx + idx, 0.0);
      break;
    }
    idx++;
  }
  assert(i < lst_size(traversal));
  return;
}
  

/* initializes all branch lengths to designated constant; initializes
   rate-matrix params using specified kappa; kappa ignored for JC69;
   REV and UNREST initialized as if HKY.  Initializes alpha as
   specified, if dgamma.  In the case of empirical rates, uniform
   weights are used for initialization. */
Vector *tm_params_init(TreeModel *mod, double branchlen, double kappa,
                       double alpha) {
  int nparams = tm_get_nparams(mod);
  Vector *params = vec_new(nparams);
  int i, nbranches, alph_size;

  tm_setup_params(mod);

  /* initialize branch-length parameters */
  nbranches = tm_get_nbranchlenparams(mod);
  vec_set(params, mod->scale_idx, 1.0);
  vec_set(params, mod->scale_idx+1, 1.0);

  for (i = 0; i < nbranches; i++) 
    vec_set(params, mod->bl_idx+i, branchlen);
  tm_init_rootleaf(mod, params);

  if (mod->backgd_freqs == NULL) {
    alph_size = strlen(mod->rate_matrix->states);
    for (i = 0; i < alph_size; i++)
      vec_set(params, mod->backgd_idx+i, 1.0/(double)alph_size);
  }
  else {
    alph_size = mod->backgd_freqs->size;
    for (i = 0; i < alph_size; i++)
      vec_set(params, mod->backgd_idx+i, vec_get(mod->backgd_freqs, i));
  }

  if (mod->nratecats > 1) {
    if (mod->empirical_rates)
      for (i = 0; i < mod->nratecats; i++) 
        vec_set(params, mod->ratevar_idx+i, mod->freqK[i]);
    else
      vec_set(params, mod->ratevar_idx, alpha);
  }

  /* initialize rate-matrix parameters */
  tm_rate_params_init(mod, params, mod->ratematrix_idx, kappa);

  if (mod->alt_subst_mods != NULL) {
    AltSubstMod *altmod;
    subst_mod_type tempmod = mod->subst_mod;
    int j;
    for (i = 0; i<lst_size(mod->alt_subst_mods); i++) {
      altmod = lst_get_ptr(mod->alt_subst_mods, i);
      mod->subst_mod = altmod->subst_mod;
      tm_rate_params_init(mod, params, altmod->ratematrix_idx, kappa);
      
      if (mod->backgd_freqs == NULL) {
	for (j=0; j<alph_size; j++)
	  vec_set(params, altmod->backgd_idx+j, 1.0/(double)alph_size);
      }
      else {
	for (j=0; j<alph_size; j++)
	  vec_set(params, altmod->backgd_idx+j, vec_get(mod->backgd_freqs, j));
      }
    }
    mod->subst_mod = tempmod;
  }
  return params;
}


/* initializes branch lengths and rate matrix parameters
   randomly; can be used multiple times to ensure the MLE is real */
Vector *tm_params_init_random(TreeModel *mod) {
  int i;
  Vector *params = vec_new(tm_get_nparams(mod));
  int nbranches = tm_get_nbranchlenparams(mod);
  int nratematparams = tm_get_nratematparams(mod);

  assert(!mod->estimate_backgd && 
         mod->estimate_branchlens == TM_BRANCHLENS_ALL);

  srandom(time(NULL));

  tm_setup_params(mod);
  
  vec_set(params, mod->scale_idx, 1.0);
  vec_set(params, mod->scale_idx+1, 1.0);
  
  for (i = 0; i < nbranches; i++)
    vec_set(params, mod->bl_idx+i, 
                   0.01 + random() * (0.5 - 0.01) / RAND_MAX);
                                /* we'll use the interval from 0.01 to 0.5 */
  tm_init_rootleaf(mod, params);

  if (mod->nratecats > 1) {
    if (mod->empirical_rates) { /* empirical rate model (category weights) */
      double val, sum = 0;
      for (i = 0; i < mod->nratecats; i++) {
        vec_set(params, mod->ratevar_idx + i, 
                       val = 0.1 + random() * (1 - 0.5) / RAND_MAX);
                                /* we'll use the interval from 0.1 to 1 */
        sum += val;
      }
      for (i = 0; i < mod->nratecats; i++) /* normalize */
        vec_set(params, mod->ratevar_idx+i, 
                       vec_get(params, mod->ratevar_idx + i) / sum);
    }

    else                        /* discrete gamma (alpha) */                       
      vec_set(params, mod->ratevar_idx, 
                     0.5 + random() * (10 - 0.5) / RAND_MAX);
                                /* we'll use the interval from 0.5 to 10 */
  }

  for (i = 0; i < nratematparams; i++) 
    vec_set(params, mod->ratematrix_idx+i, 
                   0.1 + random() * (5 - 0.1) / RAND_MAX);
                                /* we'll use the interval from 0.1 to 5 */

  if (mod->alt_subst_mods != NULL) {
    AltSubstMod *altmod;
    subst_mod_type temp_mod = mod->subst_mod;
    int j;
    for (j=0; j<lst_size(mod->alt_subst_mods); j++) {
      altmod = lst_get_ptr(mod->alt_subst_mods, j);
      mod->subst_mod = altmod->subst_mod;
      nratematparams = tm_get_nratematparams(mod);
      for (i=0; i<nratematparams; i++)
	vec_set(params, altmod->ratematrix_idx+i, 0.1 + 
		random()*(5-0.1)/RAND_MAX);
      //don't worry about backgd_freq since it is not estimated
    }
    mod->subst_mod = temp_mod;
  }

  return params;
}


/* Functions to initialize a parameter vector from an existing tree model */
Vector *tm_params_new_init_from_model(TreeModel *mod) {
  Vector *params = vec_new(tm_get_nparams(mod));
  tm_params_init_from_model(mod, params, 0);
  return params;
}

void tm_params_init_from_model(TreeModel *mod, Vector *params, 
                               int start_idx) {
  int nodeidx, j, i;
  List *traversal;
  TreeNode *n;

  tm_setup_params(mod);

  vec_set(params, mod->scale_idx, mod->scale);
  vec_set(params, mod->scale_idx+1, mod->scale_sub);
  traversal = tr_preorder(mod->tree);
  i=0;
  for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
    n = lst_get_ptr(traversal, nodeidx);
    if (n->parent == NULL) continue;

    /* Note: if the model is reversible, then the distances from root
       to lchild equals distance from root to rchild.  Setting each to
       average distance as a convention.
    */
    if ((n == mod->tree->lchild || n == mod->tree->rchild) && 
	tm_is_reversible(mod->subst_mod)) 
      vec_set(params, mod->bl_idx+i, (mod->tree->lchild->dparent +
				      mod->tree->rchild->dparent)/2.0);
    else vec_set(params, mod->bl_idx+i, n->dparent);
    i++;
  }
  tm_init_rootleaf(mod, params);

  if (mod->backgd_freqs != NULL) {
    for (i=0; i<mod->backgd_freqs->size; i++)
      vec_set(params, mod->backgd_idx+i, vec_get(mod->backgd_freqs, i));
  }

  /* next parameters are for rate variation */
  if (mod->nratecats > 1 || mod->alpha < 0) {
    if (mod->empirical_rates) 
      for (j = 0; j < mod->nratecats; j++)
        vec_set(params, mod->ratevar_idx+j, mod->freqK[j]);
    else                        /* discrete gamma model */
      vec_set(params, mod->ratevar_idx, mod->alpha);
  }

  /* initialize rate-matrix parameters */
  //  if (mod->estimate_ratemat)
    tm_rate_params_init_from_model(mod, params, mod->ratematrix_idx);

  /* Now initialize alt_models */
  if (mod->alt_subst_mods != NULL) {
    subst_mod_type tempmod = mod->subst_mod;
    AltSubstMod *altmod;
    for (j=0; j<lst_size(mod->alt_subst_mods); j++) {
      altmod = lst_get_ptr(mod->alt_subst_mods, j);
      mod->subst_mod = altmod->subst_mod;
      tm_rate_params_init_from_model(mod, params, altmod->ratematrix_idx);
      if (mod->backgd_freqs != NULL)
	for (i=0; i<mod->backgd_freqs->size; i++)
	  vec_set(params, altmod->backgd_idx+i, vec_get(mod->backgd_freqs, i));
    }
    mod->subst_mod = tempmod;
  }

}

/* Given a codon model, create and return the induced amino acid model */
TreeModel *tm_induced_aa(TreeModel *codon_mod) {
  TreeModel *retval;
  char *codon_to_aa = get_codon_mapping(codon_mod->rate_matrix->states);
  Vector *aa_freqs = vec_new(strlen(AA_ALPHABET));
  MarkovMatrix *aa_mat = mm_new(strlen(AA_ALPHABET), AA_ALPHABET, CONTINUOUS);
  int i, j;
  int nstates = codon_mod->rate_matrix->size;

  assert(codon_mod->order == 2);
  vec_zero(aa_freqs);

  /* compute induced equilibrium freqs */
  for (i = 0; i < nstates; i++) {
    int aa_state_i = aa_mat->inv_states[(int)codon_to_aa[i]]; 
                                /* state number corresponding to amino
                                   acid which corresponds to codon i */
    vec_set(aa_freqs, aa_state_i, vec_get(aa_freqs, aa_state_i) + 
                   vec_get(codon_mod->backgd_freqs, i));
  }

  /* compute induced matrix */
  for (i = 0; i < nstates; i++) {
    int aa_state_i = aa_mat->inv_states[(int)codon_to_aa[i]]; 
    double cond_freq = safediv(vec_get(codon_mod->backgd_freqs, i),
                               vec_get(aa_freqs, aa_state_i)); 
                                /* freq of codon i given its aa */
    for (j = 0; j < nstates; j++) {
      int aa_state_j = aa_mat->inv_states[(int)codon_to_aa[j]]; 
      mm_set(aa_mat, aa_state_i, aa_state_j, 
             mm_get(aa_mat, aa_state_i, aa_state_j) + 
             cond_freq * mm_get(codon_mod->rate_matrix, i, j));
    }
  }

  retval = tm_new(tr_create_copy(codon_mod->tree), aa_mat, aa_freqs, 
                  UNDEF_MOD, NULL, codon_mod->nratecats, codon_mod->alpha, 
                  NULL, codon_mod->root_leaf_id);

  retval->msa = codon_mod->msa;
  for (i = 0; i < codon_mod->nratecats; i++) {
    retval->freqK[i] = codon_mod->freqK[i];
    retval->rK[i] = codon_mod->rK[i];
  }
  retval->empirical_rates = codon_mod->empirical_rates;

  /* NOTE: ignoring params, tree_posteriors, etc. */

  free(codon_to_aa);
  return retval;
}

/* return number of parameters for specified TreeModel (based on
   number of taxa and substitution model) */
int tm_get_nparams(TreeModel *mod) {
  int rv, neqfreq, nummod=1;
  
  /*old code:
  return (tm_get_nbranchlenparams(mod) + tm_get_neqfreqparams(mod) + 
          tm_get_nratevarparams(mod) + tm_get_nratematparams(mod));    
   */
  
  //2 for scale and subtree scale, subtract one for root
  rv = 2 + mod->tree->nnodes - 1 + tm_get_nratevarparams(mod) + 
    tm_get_nratematparams(mod);
  neqfreq = tm_get_neqfreqparams(mod);
  if (mod->alt_subst_mods != NULL) {
    subst_mod_type temp_mod = mod->subst_mod;
    AltSubstMod *altmod;
    int i;
    for (i=0; i<lst_size(mod->alt_subst_mods); i++) {
      altmod = lst_get_ptr(mod->alt_subst_mods, i);
      mod->subst_mod = altmod->subst_mod;
      rv += tm_get_nratematparams(mod);
      nummod++;
    }
    mod->subst_mod = temp_mod;
  }
  rv += nummod*neqfreq;
  return rv;
}

int tm_get_neqfreqparams(TreeModel *mod) {
  if (mod->backgd_freqs == NULL) return strlen(mod->rate_matrix->states);
  return mod->backgd_freqs->size;
}

int tm_get_nratevarparams(TreeModel *mod) {
  if (mod->nratecats > 1 || mod->alpha < 0) {
    if (!mod->empirical_rates)  /* discrete gamma (alpha parameter) */
      return 1;
    else if (mod->nratecats > 1) /* empirical rates */
      return mod->nratecats;
    else if (mod->alpha < 0)    /* empirical rates but temp. disabled */
      return -mod->alpha;       /* (here, nratecats is stored as -alpha) */
  }
  return 0;
}

/** Return number of branch length params */
int tm_get_nbranchlenparams(TreeModel *mod) {
  int retval;
  return mod->tree->nnodes - 1;


  if (mod->estimate_branchlens == TM_BRANCHLENS_NONE) return 0;
  else if (mod->estimate_branchlens == TM_SCALE_ONLY) {
    if (mod->subtree_root == NULL) return 1;
    else return 2;
  }
  retval = tm_is_reversible(mod->subst_mod) ? mod->tree->nnodes - 2 :
    mod->tree->nnodes - 1;
  return retval;
}

/** Build index of leaf ids to sequence indices in given alignment.
    Leaves not present in the alignment will be ignored.  Also, it's
    not required that there's a leaf for every sequence in the
    alignment. */
void tm_build_seq_idx(TreeModel *mod, MSA *msa) {
  int i, idx;
  mod->msa_seq_idx = smalloc(mod->tree->nnodes * sizeof(int));
  for (i = 0; i < mod->tree->nnodes; i++) {
    TreeNode *n = lst_get_ptr(mod->tree->nodes, i);
    mod->msa_seq_idx[i] = -1;
    if (n->lchild == NULL && n->rchild == NULL &&
        (idx = msa_get_seq_idx(msa, n->name)) >= 0)
      mod->msa_seq_idx[i] = idx;
  }
}

/** Prune away leaves in tree that don't correspond to sequences in a
    given alignment.  Warning: root of tree (value of mod->tree) may
    change. */
void tm_prune(TreeModel *mod,   /**< TreeModel whose tree is to be pruned  */
              MSA *msa,         /**< Alignment; all leaves whose names
                                    are not in msa->names will be
                                    pruned away */
              List *names       /**< Will contain names of deleted
                                   leaves on return.  Must be
                                   pre-allocated */
              ) {
  int i, j, old_nnodes = mod->tree->nnodes;

  assert(mod->tree->nnodes >= 3);

  lst_clear(names);
  for (i = 0; i < msa->nseqs; i++)
    lst_push_ptr(names, str_new_charstr(msa->names[i]));

  tr_prune(&mod->tree, names, TRUE);

  if (mod->tree == NULL)
    return;                     /* whole tree pruned away; special case */

  if (lst_size(names) > 0) {
    /* free memory for eliminated nodes */
    for (i = mod->tree->nnodes; i < old_nnodes; i++) {
      for (j = 0; j < mod->nratecats; j++)
        if (mod->P[i][j] != NULL) mm_free(mod->P[i][j]);
      free(mod->P[i]);
    }
  }

  /* necessary because parameter indices can change */
  if (mod->rate_matrix_param_row != NULL) {
    tm_free_rmp(mod);
    tm_init_rmp(mod);
  }
}

/** Extrapolate tree model and prune leaves not represented in
   alignment (see tr_scale_by_subtree).  Returns scale factor */
double tm_extrapolate_and_prune(TreeModel *mod, TreeNode *extrapolate_tree, 
                                MSA *msa, List *pruned_names) {
  int i;
  TreeNode *t = tr_create_copy(extrapolate_tree);
  double scale = tr_scale_by_subtree(t, mod->tree);
  for (i = 0; i < msa->nseqs; i++)
    lst_push_ptr(pruned_names, str_new_charstr(msa->names[i]));
  tr_prune(&t, pruned_names, TRUE);
  tm_reset_tree(mod, t);
  return scale;
}

/** Reset TreeModel with new or altered tree. */
void tm_reset_tree(TreeModel *mod,   /** TreeModel */
                   TreeNode *newtree /** New tree */
                   ) {
  /* merge this with tm_reinit? */
  int i, j;

  /* free P matrices */
  for (i = 0; i < mod->tree->nnodes; i++) {
    for (j = 0; j < mod->nratecats; j++)
      if (mod->P[i][j] != NULL) mm_free(mod->P[i][j]);
    free(mod->P[i]);
  }

  tm_free_rmp(mod);             /* necessary because parameter indices
                                   can change */
  tr_free(mod->tree);
  mod->tree = newtree;
  tm_init_rmp(mod);

  /* realloc P matrices */
  mod->P = srealloc(mod->P, mod->tree->nnodes * sizeof(void**));
  for (i = 0; i < mod->tree->nnodes; i++) {
    mod->P[i] = smalloc(mod->nratecats * sizeof(MarkovMatrix*));
    for (j = 0; j < mod->nratecats; j++) mod->P[i][j] = NULL;
  }
}

/* Set branches to be ignored in likelihood calculation and parameter
   estimation.  Argument 'ignore_branches' should be a list of Strings
   indicating nodes in the tree and the branches leading to those
   nodes. */
void tm_set_ignore_branches(TreeModel *mod, List *ignore_branches) {
  int j;
  assert(mod->ignore_branch == NULL);
  mod->ignore_branch = smalloc(mod->tree->nnodes * sizeof(int));
  for (j = 0; j < mod->tree->nnodes; j++) 
    mod->ignore_branch[j] = FALSE;
  for (j = 0; j < lst_size(ignore_branches); j++) {
    String *s = lst_get_ptr(ignore_branches, j);
    TreeNode *n  = tr_get_node(mod->tree, s->chars);
    if (n == NULL)
      die("ERROR: no node named '%s'.\n", s->chars);
    mod->ignore_branch[n->id] = TRUE;
  }
}

AltSubstMod *tm_new_alt_subst_mod(subst_mod_type subst_mod,
                                  Vector *backgd_freqs, 
                                  MarkovMatrix *rate_matrix) {
  AltSubstMod *am = smalloc(sizeof(AltSubstMod));
  am->subst_mod = subst_mod;
  am->backgd_freqs = backgd_freqs;
  am->rate_matrix = rate_matrix;
  return am;
}

void tm_free_alt_subst_mod(AltSubstMod *am) {
  int i;
  if (am->backgd_freqs != NULL)
    vec_free(am->backgd_freqs);
  if (am->rate_matrix != NULL)
    mm_free(am->rate_matrix);
  if (am->param_list != NULL) {
    for (i = 0; i < lst_size(am->param_list); i++) 
      str_free(lst_get_ptr(am->param_list, i));
    lst_free(am->param_list);
  }
  str_free(am->nodename);
  free(am);
}
