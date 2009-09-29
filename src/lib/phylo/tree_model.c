/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2009 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: tree_model.c,v 1.42 2009-03-20 18:16:48 agd27 Exp $ */

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
    tm->P = (MarkovMatrix***)smalloc(tree->nnodes * sizeof(MarkovMatrix**));
    for (i = 0; i < tree->nnodes; i++) {
      tm->P[i] = (MarkovMatrix**)smalloc(nratecats * sizeof(MarkovMatrix*));
      for (j = 0; j < nratecats; j++) tm->P[i][j] = NULL;
    }

    tm->rK = (double*)smalloc(nratecats * sizeof(double));
    tm->freqK = (double*)smalloc(nratecats * sizeof(double));

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
  tm_init_rmp(tm);

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
  if (tm_order(new_subst_mod) != tm->order)
    die("ERROR: Cannot initialize tree model from model of different order.\n");
  tm->subst_mod = new_subst_mod;
  tm->nratecats = new_nratecats;
  tm->alpha = (tm->nratecats > 1 ? new_alpha : 0);
  tm->empirical_rates = new_rate_consts != NULL ? 1 : 0;
  tm_init_rmp(tm);
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
      if (tm->alt_subst_mods != NULL && tm->alt_subst_mods[i] != NULL) 
        tm_free_alt_subst_mod(tm->alt_subst_mods[i]);
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
  if (tm->alt_subst_mods != NULL) free(tm->alt_subst_mods);
  free(tm);
}

/* set up lists that allow each rate matrix parameter to be mapped to
   the rows and columns in which it appears */
void tm_init_rmp(TreeModel *tm) {
  int i;
  if (tm->rate_matrix == NULL || tm->subst_mod == UNDEF_MOD || 
      tm->tree == NULL) 
    tm->rate_matrix_param_row = tm->rate_matrix_param_col = NULL;
  else {
    int nparams = tm_get_nparams(tm);
    int nrmparams = tm_get_nratematparams(tm);
    tm->rate_matrix_param_row = (List**)smalloc(nparams * sizeof(List*));
    tm->rate_matrix_param_col = (List**)smalloc(nparams * sizeof(List*));
    for (i = 0; i < nparams - nrmparams; i++) 
      tm->rate_matrix_param_row[i] = tm->rate_matrix_param_col[i] = NULL;
                                /* NULL for non-rate-matrix params;
                                   keeps indexing simple */
    for (; i < nparams; i++) {
      tm->rate_matrix_param_row[i] = lst_new_int(2); 
      tm->rate_matrix_param_col[i] = lst_new_int(2);
      /* usually 2 elements (general reversible models) */
    }
  }
}

void tm_free_rmp(TreeModel *tm) {
  int i;
  int nparams = tm_get_nparams(tm);
  int nrmparams = tm_get_nratematparams(tm);
  for (i = nparams - nrmparams; i < nparams; i++) {
    lst_free(tm->rate_matrix_param_row[i]);
    lst_free(tm->rate_matrix_param_col[i]);
  }  
  free(tm->rate_matrix_param_row);
  free(tm->rate_matrix_param_col);
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
    else {
      fprintf(stderr, "ERROR: unrecognized tag in model file (\"%s\").\n", 
	      tag);
      exit(1);
    }
  }

  if (empty) die("ERROR: empty tree model file.\n");

  if (rmat != NULL) {
    M = mm_new_from_matrix(rmat, alphabet, CONTINUOUS);
  }  else {
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
}

TreeModel *tm_create_copy(TreeModel *src) {
  TreeModel *retval;
  Vector *freqs = vec_new(src->backgd_freqs->size);
  int i;

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

  /* NOTE: ignoring params, tree_posteriors, etc. */

  return retval;
}

void tm_set_subst_matrices(TreeModel *tm) {
  int i, j;
  double scaling_const, tmp, branch_scale;
  Vector *backgd_freqs = tm->backgd_freqs;
  subst_mod_type subst_mod = tm->subst_mod;
  MarkovMatrix *rate_matrix = tm->rate_matrix;
  
  scaling_const = -1;

  if (tm->estimate_branchlens != TM_SCALE_ONLY && tm->scale != 1) 
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
    TreeNode *n = lst_get_ptr(tm->tree->nodes, i);

    if (n->parent == NULL) continue;

    /* special case where subst models differ by branch */
    if (tm->alt_subst_mods != NULL) {
      if (tm->alt_subst_mods[n->id] != NULL) {
        backgd_freqs = tm->alt_subst_mods[n->id]->backgd_freqs;
        subst_mod = tm->alt_subst_mods[n->id]->subst_mod;
        rate_matrix = tm->alt_subst_mods[n->id]->rate_matrix;
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

  int i, class, nseqs, col, ntreenodes, idx, ratecat;
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

    if (classmods[i]->nratecats > 1 && !classmods[i]->empirical_rates)
      DiscreteGamma(classmods[i]->freqK, classmods[i]->rK, 
		    classmods[i]->alpha, classmods[i]->alpha,
		    classmods[i]->nratecats, 0);
    //    assert(classmods[i]->nratecats == 1); /* assuming no rate variation */

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
    if (classmods[class]->nratecats > 1)
      ratecat = pv_draw_idx_arr(classmods[class]->freqK, classmods[class]->nratecats);
    else ratecat=0;

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
        if (classmods[class]->P[l->id][ratecat] == NULL)
          tm_set_subst_matrices(classmods[class]);
        lsubst_mat = classmods[class]->P[l->id][ratecat];
        rsubst_mat = classmods[class]->P[r->id][ratecat];
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

/* Generates a random alignment as above but uses a list of 
   scales and subtree scales to use.  (If subtreeName is NULL
   subtree scales must be 1).  */
MSA *tm_generate_msa_scaleLst(List *nsitesLst, List *scaleLst, 
			      List *subtreeScaleLst,
			      TreeModel *mod, char *subtreeName) {
  double scale, subtreeScale;
  int nsite, ncolumns=0, i, j, k, idx, col, nseqs, ratecat;
  char **seqs, **names, *newchar;
  MSA *msa;
  TreeNode *subtreeNode;
  List *traversal = tr_preorder(mod->tree);

  if (subtreeName == NULL) {
    for (i=0; i<lst_size(subtreeScaleLst); i++) {
      subtreeScale = lst_get_dbl(subtreeScaleLst, i);
      if (subtreeScale != 1.0) {
	fprintf(stderr, "Warning: ignoring subtreeScales because subtree not given\n");
	break;
      }
    }
  }
  else {
    subtreeNode = tr_get_node(mod->tree, subtreeName);
    if (subtreeNode == NULL) die("ERROR: no node with name %s\n", subtreeName);
  }
  for (i=0; i<lst_size(nsitesLst); i++) 
    ncolumns += lst_get_int(nsitesLst, i);
  nseqs = (mod->tree->nnodes+1)/2;
  
    /* create new MSA */
  names = (char**)smalloc(nseqs * sizeof(char*));
  seqs = (char**)smalloc(nseqs * sizeof(char*));
  for (i = 0; i < nseqs; i++) 
    seqs[i] = (char*)smalloc((ncolumns + 1) * sizeof(char));
  msa = msa_new(seqs, names, nseqs, ncolumns, 
                mod->rate_matrix->states);
  
  mod->msa_seq_idx = smalloc(mod->tree->nnodes * sizeof(int));
  for (i = 0, idx = 0; i < mod->tree->nnodes; i++) {
    TreeNode *n = (TreeNode*)lst_get_ptr(traversal, i);
    if (n->lchild == NULL && n->rchild == NULL) {
      mod->msa_seq_idx[i] = idx;
      names[idx] = strdup(n->name);
      idx++;
    }
    else mod->msa_seq_idx[i] = -1;
  }

  newchar = (char*)smalloc(mod->tree->nnodes * sizeof(char));
  col=0;
  for (i=0; i<lst_size(nsitesLst); i++) {
    nsite = lst_get_int(nsitesLst, i);
    scale = lst_get_dbl(scaleLst, i);
    subtreeScale = lst_get_dbl(subtreeScaleLst, i);
    tr_scale(mod->tree, scale);
    if (subtreeName != NULL) 
      tr_scale_subtree(mod->tree, subtreeNode, subtreeScale);
    tm_set_subst_matrices(mod);
    for (j=0; j<nsite; j++) {
      if (mod->nratecats > 1) 
	ratecat = pv_draw_idx_arr(mod->freqK, mod->nratecats);
      else ratecat=0;
      newchar[mod->tree->id] = 
	mm_sample_backgd(mod->rate_matrix->states, 
			 mod->backgd_freqs);
      
      for (k=0; k<lst_size(traversal); k++) {
	TreeNode *n = (TreeNode*)lst_get_ptr(traversal, k);
	TreeNode *l = n->lchild;
	TreeNode *r = n->rchild;
	assert( (l==NULL && r==NULL) || (l!=NULL && r!=NULL));
	if (l == NULL)
	  msa->seqs[mod->msa_seq_idx[n->id]][col] = newchar[n->id];
	else {
	  newchar[l->id] = mm_sample_char(mod->P[l->id][ratecat], newchar[n->id]);
	  newchar[r->id] = mm_sample_char(mod->P[r->id][ratecat], newchar[n->id]);
	}
      }
      col++;
    }
    tr_scale_subtree(mod->tree, subtreeNode, 1.0/subtreeScale);
    tr_scale(mod->tree, 1.0/scale);
  }
  tm_set_subst_matrices(mod);
  free(newchar);
  return msa;
}


/* Generates an alignment according to set of Tree Models and a
   Markov matrix defing how to transition among them.  TreeModels must
   appear in same order as the states of the Markov matrix. 
   NOTE: call srandom externally. */
MSA *tm_generate_msa_random_subtree(int ncolumns, TreeModel *mod,
				    TreeModel *subtreeMod, char *subtree,
				    double subtreeSwitchProb) {
  int i, nseqs, col, idx, ratecat;
  MSA *msa;
  char *newchar;
  char **names, **seqs;
  TreeNode *subtreeNode;
  List *traversal = tr_preorder(mod->tree);
  int *inSubtree;

  subtreeNode = tr_get_node(mod->tree, subtree);
  if (subtreeNode==NULL) {
    die("ERROR: no node with name %s\n", subtree);
  }
  nseqs = (mod->tree->nnodes+1)/2;
  inSubtree = tr_in_subtree(mod->tree, subtreeNode);
  inSubtree[subtreeNode->id] = TRUE;

  /* create new MSA */
  names = (char**)smalloc(nseqs * sizeof(char*));
  seqs = (char**)smalloc(nseqs * sizeof(char*));
  for (i = 0; i < nseqs; i++) 
    seqs[i] = (char*)smalloc((ncolumns + 1) * sizeof(char));
  msa = msa_new(seqs, names, nseqs, ncolumns, 
                mod->rate_matrix->states);

  /* build sequence idx map; only need one for first model */
  mod->msa_seq_idx = smalloc(mod->tree->nnodes * sizeof(int));
  for (i = 0, idx = 0; i < mod->tree->nnodes; i++) {
    TreeNode *n = (TreeNode*)lst_get_ptr(traversal, i);
    if (n->lchild == NULL && n->rchild == NULL) {
      mod->msa_seq_idx[i] = idx;
      names[idx] = strdup(n->name);
      idx++;
    }
    else mod->msa_seq_idx[i] = -1;
  }

  /* generate sequences, column by column */
  newchar = (char*)smalloc(mod->tree->nnodes * sizeof(char));
  for (col = 0; col < ncolumns; col++) {
    if (mod->nratecats > 1)
      ratecat = pv_draw_idx_arr(mod->freqK, mod->nratecats);
    else ratecat=0;

    newchar[mod->tree->id] = 
      mm_sample_backgd(mod->rate_matrix->states, 
                       mod->backgd_freqs);
    for (i = 0; i < lst_size(traversal); i++) {
      TreeNode *n = (TreeNode*)lst_get_ptr(traversal, i);
      TreeNode *l = n->lchild;
      TreeNode *r = n->rchild;
      int j, inSub[2];
      TreeModel *lmod, *rmod;
      assert ((l == NULL && r == NULL) || (l != NULL && r != NULL));

      if (l == NULL) 
        msa->seqs[mod->msa_seq_idx[n->id]][col] = newchar[n->id];
      else {
        MarkovMatrix *lsubst_mat, *rsubst_mat;
        inSub[0] = inSubtree[l->id];
        inSub[1] = inSubtree[r->id];

	for (j=0; j<2; j++) {
	  if (1.0*(double)random()/(double)RAND_MAX < subtreeSwitchProb) 
	    inSub[j]=!inSub[j];
	}
	if (inSub[0]) lmod = subtreeMod;
	else lmod = mod;
	if (inSub[1]) rmod = subtreeMod;
	else rmod=mod;
	
	if (lmod->P[l->id][ratecat]==NULL) {
	  //	  printf("setting subst matcies inSub=%i\n", inSub[0]);
	  tm_set_subst_matrices(lmod);
	}
	if (rmod->P[r->id][ratecat]==NULL) {
	  //	  printf("setting subst matrices inSub=%i\n", inSub[1]);
	  tm_set_subst_matrices(rmod);
	}
        lsubst_mat = lmod->P[l->id][ratecat];
        rsubst_mat = rmod->P[r->id][ratecat];
        newchar[l->id] = mm_sample_char(lsubst_mat, newchar[n->id]);
        newchar[r->id] = mm_sample_char(rsubst_mat, newchar[n->id]);
      }
    }
  }
  free(newchar);
  free(inSubtree);
  return msa;
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
  Vector *lower_bounds, *upper_bounds;
  int retval = 0;

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
  }

  if (mod->tree == NULL) {      /* weight matrix */
    mod->lnL = tl_compute_log_likelihood(mod, msa, NULL, cat, NULL) * 
      log(2);
    return 0;
  }

  mod->msa = msa;               /* package with mod any data needed to
                                   compute likelihoods */
  mod->category = cat;

  /* most params have lower bound of zero and no upper bound */
  lower_bounds = vec_new(params->size);
  vec_zero(lower_bounds);
  upper_bounds = NULL;

  /* however, in this case we don't want the eq freqs to go to zero */
  if (mod->estimate_backgd) {
    int i, offset = tm_get_nbranchlenparams(mod);
    for (i = 0; i < mod->backgd_freqs->size; i++)
      vec_set(lower_bounds, i + offset, 0.001);
  }

  /* Also, in this case, we need to bound the scale of the subtree */
  if (mod->estimate_branchlens == TM_SCALE_ONLY && mod->subtree_root != NULL 
      && mod->scale_sub_bound != NB) {
    if (mod->scale_sub_bound == LB) vec_set(lower_bounds, 1, 1);
    if (mod->scale_sub_bound == UB) {
      upper_bounds = vec_new(params->size);
      vec_set_all(upper_bounds, INFTY);
      vec_set(upper_bounds, 1, 1);
    }
  }

  retval = opt_bfgs(tm_likelihood_wrapper, params, (void*)mod, &ll, 
                    lower_bounds, upper_bounds, logf, NULL, precision, NULL);

  mod->lnL = ll * -1 * log(2);  /* make negative again and convert to
                                   natural log scale */
  tm_unpack_params(mod, params, -1);

  if (mod->subst_mod != JC69 && mod->subst_mod != F81) {   
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

  vec_free(lower_bounds);
  if (upper_bounds != NULL) vec_free(upper_bounds);

  if (retval != 0) 
    fprintf(stderr, "WARNING: BFGS algorithm reached its maximum number of iterations.\n");

  return retval;
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
void tm_unpack_params(TreeModel *mod, Vector *params, int idx_offset) {
  TreeNode *n;
  int nparams = tm_get_nparams(mod);
  int assigned = 0, nodeidx, i, j;
  List *traversal;

  if (idx_offset == -1) {       /* don't do if nested */
    assert(params->size == nparams && mod->tree->nnodes >= 3);
    /* this is a check on parameter values that helps catch errors
       before they've become untraceable.  Make sure all params are
       finite.  As of yet, all models use nonnegative parameters only,
       so also make sure all params are nonnegative. */
    for (i = 0; i < params->size; i++) {
      double p = vec_get(params, i);
      if (p < 0 && abs(p) < TM_IMAG_EPS) /* consider close enough to 0 */
        vec_set(params, i, p=0);
      if (p < 0) die("ERROR: parameter %d has become negative (%f).\n", i, p);
      if (!finite(p)) die("ERROR: parameter %d is no longer finite (%f).\n", 
                          i, p);
    }
    i = 0;
  }
  else 
    i = idx_offset;

  if (mod->estimate_branchlens == TM_SCALE_ONLY) {
    if (mod->empirical_rates && (mod->nratecats > 1 || mod->alpha < 0)) i++;
                                /* in this case, skip scale */
    else {
      mod->scale = vec_get(params, i++);
      if (mod->subtree_root) /* estimating subtree scale */
        mod->scale_sub = vec_get(params, i++);
    }
  }
  else if (mod->estimate_branchlens == TM_BRANCHLENS_CLOCK) { 
    /* molecular clock; set total height of each node from parameter
       vector; temporarily store in dparent */
    for (nodeidx = 0; nodeidx < mod->tree->nnodes; nodeidx++) {
      n = lst_get_ptr(mod->tree->nodes, nodeidx);
      if (n->lchild == NULL)  /* leaf */
        n->dparent = 0;
      else 
        n->dparent = vec_get(params, i++);
    }
    /* set branch lengths */
    traversal = tr_postorder(mod->tree);
    /* first ensure ordering constraints are obeyed; can be violated
       during parameter estimation */
    for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
      n = lst_get_ptr(traversal, nodeidx);
      if (n->parent != NULL && n->dparent > n->parent->dparent)
        n->parent->dparent = n->dparent;
    }
    /* now set branch lengths to differences in heights */
    for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
      n = lst_get_ptr(traversal, nodeidx);
      if (n->parent == NULL)
        n->dparent = 0;
      else
        n->dparent = n->parent->dparent - n->dparent;
    }      
  }
  else if (mod->estimate_branchlens == TM_BRANCHLENS_ALL) {
    /* no clock; first nnodes-2 elements define branch lengths */
    /* NOTE: order of traversal is not important as long as it is
       consistently obeyed */
    traversal = tr_preorder(mod->tree);
    for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
      n = lst_get_ptr(traversal, nodeidx);

      if (n->parent != NULL) {
        /* if the model is reversible, divide the first parameter evenly
           between the two branches from the root; this is a simple way
           to simulate an unrooted tree */
        if ((n == mod->tree->lchild || n == mod->tree->rchild) && 
            tm_is_reversible(mod->subst_mod)) {
          n->dparent = vec_get(params, 0)/2;
          if (!assigned) {
            i++;     /* only increment the first time */
            assigned = 1;
          }
        }
        else if (n->id == mod->root_leaf_id) 
          n->dparent = 0;
        else 
          n->dparent = vec_get(params, i++);
      }
    }
  }

  /* if estimating backgd, next params are backgd freqs */
  if (mod->estimate_backgd) {
    double sum = 0;
    for (j = 0; j < mod->backgd_freqs->size; j++) 
      sum += vec_get(params, i+j);
    for (j = 0; j < mod->backgd_freqs->size; j++)
      vec_set(mod->backgd_freqs, j, vec_get(params, i++)/sum);
  }

  /* next parameters are for rate variation */
  if (mod->nratecats > 1) {
    if (mod->empirical_rates) {
      for (j = 0; j < mod->nratecats; j++)
        mod->freqK[j] = vec_get(params, i++);
                                /* assume normalized */
    }
    else {                      /* discrete gamma model */
      mod->alpha = vec_get(params, i++);
      DiscreteGamma(mod->freqK, mod->rK, mod->alpha, mod->alpha, 
                    mod->nratecats, 0); 
    }
  }
  else if (mod->alpha < 0) {    /* code to ignore rate variation temporarily */
    if (mod->empirical_rates) i += -mod->alpha; 
                                /* original nratecats is stored as -alpha */
    else i++;
  }

  /* (remaining parameters define the subst matrix) */ 

  if (mod->estimate_ratemat) {
    /* redefine substitution matrix */
    tm_set_rate_matrix(mod, params, i);

    /* diagonalize, if necessary */
    if (mod->subst_mod != JC69 && mod->subst_mod != F81)
      mm_diagonalize(mod->rate_matrix);
  }

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

/* initializes all branch lengths to designated constant; initializes
   rate-matrix params using specified kappa; kappa ignored for JC69;
   REV and UNREST initialized as if HKY.  Initializes alpha as
   specified, if dgamma.  In the case of empirical rates, uniform
   weights are used for initialization. */
Vector *tm_params_init(TreeModel *mod, double branchlen, double kappa,
                       double alpha) {
  int nparams = tm_get_nparams(mod);
  Vector *params = vec_new(nparams);
  int i, nbranches, params_idx = 0;

  /* initialize branch-length parameters */
  nbranches = tm_get_nbranchlenparams(mod);

  if (mod->estimate_branchlens == TM_BRANCHLENS_CLOCK)  {     
                                /* set parameter to constant times
                                   height of corresponding node
                                   (parameters represent total branch
                                   length to ancestral nodes)  */
    for (i = 0; i < mod->tree->nnodes; i++) {
      TreeNode *n = lst_get_ptr(mod->tree->nodes, i);
      if (n->lchild != NULL)  /* non-leaf */
        vec_set(params, params_idx++, branchlen * n->height);
    }
  }
  else {                        /* simply set branch length to constant */
    for (params_idx = 0; params_idx < nbranches; params_idx++) 
      vec_set(params, params_idx, branchlen);
  }

  if (mod->estimate_backgd) {
    if (mod->backgd_freqs == NULL) {
      double alph_size = strlen(mod->rate_matrix->states);
      for (i = 0; i < alph_size; i++)
        vec_set(params, params_idx++, 1.0/alph_size);
    }
    else
      for (i = 0; i < mod->backgd_freqs->size; i++)
        vec_set(params, params_idx++, vec_get(mod->backgd_freqs, i));
  }

  if (mod->nratecats > 1) {
    if (mod->empirical_rates)
      for (i = 0; i < mod->nratecats; i++) 
        vec_set(params, params_idx++, mod->freqK[i]);
    else
      vec_set(params, params_idx++, alpha);
  }

  /* initialize rate-matrix parameters */
  tm_rate_params_init(mod, params, params_idx, kappa);

  return params;
}



/* initializes branch lengths and rate matrix parameters
   randomly; can be used multiple times to ensure the MLE is real */
Vector *tm_params_init_random(TreeModel *mod) {
  int i, params_idx = 0;
  Vector *params = vec_new(tm_get_nparams(mod));
  int nbranches = tm_get_nbranchlenparams(mod);
  int nratematparams = tm_get_nratematparams(mod);
  TreeNode *n;
  List *traversal;
  double *heights;

  assert(!mod->estimate_backgd && 
         (mod->estimate_branchlens == TM_BRANCHLENS_ALL ||
          mod->estimate_branchlens == TM_BRANCHLENS_CLOCK));

  srandom(time(NULL));

  if  (mod->estimate_branchlens == TM_BRANCHLENS_CLOCK) { 
    /* in this case, have to ensure ancestral nodes have larger values
       than their descendants */
    traversal = tr_postorder(mod->tree);
    heights = smalloc(mod->tree->nnodes * sizeof(double));
    for (i = 0; i < lst_size(traversal); i++) {
      n = lst_get_ptr(traversal, i);
      if (n->lchild == NULL)    /* leaf */
        heights[n->id] = 0;
      else
        heights[n->id] = 
          max(heights[n->lchild->id], heights[n->rchild->id]) +
          0.01 + random() * (0.5 - 0.01) / RAND_MAX;
    }
    for (i = 0; i < mod->tree->nnodes; i++) { /* has to follow index
                                                 order, excluding leaves */
      n = lst_get_ptr(mod->tree->nodes, i);
      if (n->lchild != NULL)
        vec_set(params, params_idx++, heights[n->id]);      
    }
    free(heights);
  }
  else {                        /* no clock */
    for (i = 0; i < nbranches; i++)
      vec_set(params, params_idx++, 
              0.01 + random() * (0.5 - 0.01) / RAND_MAX);
              /* we'll use the interval from 0.01 to 0.5 */
  }

  if (mod->nratecats > 1) {
    if (mod->empirical_rates) { /* empirical rate model (category weights) */
      double val, sum = 0;
      for (i = 0; i < mod->nratecats; i++) {
        vec_set(params, params_idx + i, 
                       val = 0.1 + random() * (1 - 0.5) / RAND_MAX);
                                /* we'll use the interval from 0.1 to 1 */
        sum += val;
      }
      for (i = 0; i < mod->nratecats; i++) /* normalize */
        vec_set(params, params_idx + i, 
                       vec_get(params, params_idx + i) / sum);
      params_idx += mod->nratecats;
    }

    else                        /* discrete gamma (alpha) */                       
      vec_set(params, params_idx++, 
                     0.5 + random() * (10 - 0.5) / RAND_MAX);
                                /* we'll use the interval from 0.5 to 10 */
  }

  for (i = 0; i < nratematparams; i++) 
    vec_set(params, params_idx++, 
                   0.1 + random() * (5 - 0.1) / RAND_MAX);
                                /* we'll use the interval from 0.1 to 5 */

  return params;
}


//assume that numMin and minState have already been set for leaf nodes
//recurse downward and set state to intersection if exists, 
//otherwise union
int tm_fitch_rec_down(TreeNode *tree, int *numMin, int **minState) {
  int sum;
  int i, j, child[2];

  if (tree->lchild==NULL) return 0;
  sum = tm_fitch_rec_down(tree->lchild, numMin, minState);
  sum += tm_fitch_rec_down(tree->rchild, numMin, minState);

  child[0] = tree->lchild->id;
  child[1] = tree->rchild->id;

  numMin[tree->id]=0;
  
  //first see if there's an intersection in minStates for children
  for (i=0; i<numMin[child[0]]; i++)
    for (j=0; j<numMin[child[1]]; j++)
      if (minState[child[0]][i] == minState[child[1]][j]) 
	minState[tree->id][numMin[tree->id]++] = minState[child[0]][i];
  if (numMin[tree->id] > 0) return sum;
  

  //otherwise take union and add 1 to cost
  for (i=0; i<2; i++) 
    for (j=0; j<numMin[child[i]]; j++)
      minState[tree->id][numMin[tree->id]++] = minState[child[i]][j];
  return sum+1;
}


void tm_fitch_rec_up(int *nodecost, TreeNode *tree, 
		     int *numMin, int **minState, int pState) {
  int node, state=-1, i;

  node = tree->id;
  nodecost[node]=0;
  
  //get state at node.  If only one minimum, we already have it.
  //Otherwise, if one of the minimums matches parent state, use that.
  //Otherwise pick one at random.
  if (tree->parent != NULL) {
    if (numMin[node]==1) state = minState[node][0];
    else {
      for (i=0; i<numMin[node]; i++)
	if (minState[node][i] == pState) {
	  state = pState;
	  break;
	}
      if (i == numMin[node]) 
	state = minState[node][(int)(random()/RAND_MAX*(double)numMin[node])];
    }      

    //now add cost to nodecost if state != parentState
    if (state != pState) nodecost[node]=1;
  } else state = pState;
  
  if (tree->lchild != NULL) {
    tm_fitch_rec_up(nodecost, tree->lchild, numMin, minState, state);
    tm_fitch_rec_up(nodecost, tree->rchild, numMin, minState, state);
  }
}



//initialize branchlengths to average number of mutations under parsimony.
//returns parsimony score
//If params is NULL, sets the branchlengths in mod->tree directly.
//Otherwise sets the parameter vector
double tm_params_init_branchlens_parsimony(Vector *params, TreeModel *mod, MSA *msa) {
  int i, numnode = mod->tree->nnodes;
  int numstate=strlen(msa->alphabet), **minState, *numMinState;
  int tupleIdx, spec, node, rootMinState, *nodecost, params_idx;
  double *brlen, weight, denom=0, totalCost=0;
  List *traversal;
  char currCh;
  TreeNode *currNode;

  if (msa->ss==NULL) 
    die("ERROR: tm_params_init_branchlens_parsimony needs msa->ss!=NULL\n");

  if (mod->order!=0) 
    die("ERROR: tm_params_init_branches currently only works for order 1 models\n");

  if  (mod->estimate_branchlens == TM_BRANCHLENS_CLOCK) { 
    die("Sorry, parsimony algorithm not implemented for molecular clock\n");
  }

  srandom(time(NULL));


  //this array keeps track of number of mutations on each branch
  brlen = malloc(numnode*sizeof(double));
  for (i=0; i<numnode; i++)
    brlen[i] = 0;

  //get sequence index if not already there
  if (mod->msa_seq_idx==NULL)
    tm_build_seq_idx(mod, msa);

  minState = malloc(numnode*sizeof(int*));
  for (i=0; i<numnode; i++)
    minState[i] = malloc(numstate*sizeof(int));
  numMinState = malloc(numnode*sizeof(int));
  nodecost = malloc(numnode*sizeof(int));
  traversal = tr_preorder(mod->tree);

  for (tupleIdx=0; tupleIdx<msa->ss->ntuples; tupleIdx++) {
    //get states at tips
    for (node=0; node<numnode; node++) {
      spec = mod->msa_seq_idx[node];
      if (spec==-1) continue;
      currCh = ss_get_char_tuple(msa, tupleIdx, spec, 0);
      if (msa->is_missing[(int)currCh] || currCh==GAP_CHAR) {
	numMinState[node] = numstate;
	for (i=0; i<numstate; i++) {
	  minState[node][i] = i;
	}
      } else {
	numMinState[node]=1;
	for (i=0; i<numstate; i++) 
	  if (msa->alphabet[i] == currCh) {
	    minState[node][0] = i;
	    break;
	  }
	assert(i < numstate);
      }
    }
    
    totalCost += (double)tm_fitch_rec_down(mod->tree, numMinState, minState);

    if (numMinState[mod->tree->id] > 0) 
      rootMinState = minState[mod->tree->id][(int)(random()/RAND_MAX*(double)numMinState[mod->tree->id])];
    else rootMinState = minState[mod->tree->id][0];

    tm_fitch_rec_up(nodecost, mod->tree, numMinState, minState, rootMinState);

    weight = msa->ss->counts[tupleIdx];
    denom += weight;
    for (node=0; node<numnode; node++)
      brlen[node] += nodecost[node]*weight;
  }

 // for (node=0; node<numnode; node++)
  //  brlen[node] /= denom;

 //for now try JC correction:
  for (node=0; node<numnode; node++) 
    if (brlen[node]!=0.0) {
      brlen[node] = -0.75*log(1.0-4.0*brlen[node]/(3.0*denom));
    }
  if (params == NULL) {
    for (node=0; node < numnode; node++)
      ((TreeNode*)(lst_get_ptr(mod->tree->nodes, node)))->dparent = brlen[node];
  }   
  else if (mod->estimate_branchlens == TM_BRANCHLENS_ALL) {
    params_idx = 0;
    for (node=0; node<lst_size(traversal); node++) {
      currNode = (TreeNode*)lst_get_ptr(traversal, node);
      if (currNode->parent == NULL) continue;
      if (currNode == mod->tree->lchild && tm_is_reversible(mod->subst_mod)) 
	vec_set(params, params_idx++, 
		brlen[currNode->id] + brlen[mod->tree->rchild->id]);
      else if (currNode != mod->tree->rchild || 
	       !tm_is_reversible(mod->subst_mod))
	vec_set(params, params_idx++, brlen[currNode->id]);
    }
  } else if (mod->estimate_branchlens==TM_SCALE_ONLY) {
    double sum=0.0;
    for (node=0; node<numnode; node++)
      sum += brlen[node];
    vec_set(params, 0, sum);
    if (mod->subtree_root != NULL) {
      sum=0.0;
      traversal = tr_preorder(mod->subtree_root);
      for (node=0; node<lst_size(traversal); node++) {
	currNode = (TreeNode*)lst_get_ptr(traversal,node);
	sum += brlen[currNode->id];
      }
      vec_set(params, 1, sum/vec_get(params, 0));
    }
  }
	
  free(brlen);
  for (i=0; i<numnode; i++)
    free(minState[i]);
  free(minState);
  free(numMinState);
  free(nodecost);
  //  printf("totalCost=%f\n", totalCost);
  return totalCost;
}


/* Functions to initialize a parameter vector from an existing tree model */
Vector *tm_params_new_init_from_model(TreeModel *mod) {
  Vector *params = vec_new(tm_get_nparams(mod));
  tm_params_init_from_model(mod, params, 0);
  return params;
}

void tm_params_init_from_model(TreeModel *mod, Vector *params, 
                               int start_idx) {
  int params_idx = start_idx, nodeidx, j;
  List *traversal;
  TreeNode *n;

  if (mod->estimate_branchlens == TM_SCALE_ONLY) {
    vec_set(params, params_idx++, mod->scale); 
    if (mod->subtree_root != NULL) /* estimating subtree scale */
      vec_set(params, params_idx++, mod->scale_sub);
  }
  else if (mod->estimate_branchlens == TM_BRANCHLENS_CLOCK) {
    /* compute height of each node */
    double *heights = smalloc(mod->tree->nnodes * sizeof(double));
    traversal = tr_postorder(mod->tree);
    for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
      n = lst_get_ptr(traversal, nodeidx);
      if (n->lchild == NULL)    /* leaf */
        heights[n->id] = 0;
      else 
        heights[n->id] = max(heights[n->lchild->id], heights[n->rchild->id]) + 
          n->dparent;
    }
    /* set params equal to heights, in order of ids (skipping leaves) */
    for (nodeidx = 0; nodeidx < mod->tree->nnodes; nodeidx++) {
      n = lst_get_ptr(mod->tree->nodes, nodeidx);
      if (n->lchild != NULL)
        vec_set(params, params_idx++, heights[n->id]);
    }
    free(heights);
  }
  else if (mod->estimate_branchlens == TM_BRANCHLENS_ALL) {
    traversal = tr_preorder(mod->tree);
    for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
      n = lst_get_ptr(traversal, nodeidx);
      if (n->parent == NULL) continue;

      /* Note: if the model is reversible, then the first parameter is
         the sum of the lengths of the two branches from the root */
      if (n == mod->tree->lchild && tm_is_reversible(mod->subst_mod))
        vec_set(params, params_idx++, n->dparent*2);
      else if (n != mod->tree->rchild || !tm_is_reversible(mod->subst_mod))
        vec_set(params, params_idx++, n->dparent);
    }
  }

  /* if estimating backgd, next params are backgd freqs */
  if (mod->estimate_backgd) 
    for (j = 0; j < mod->backgd_freqs->size; j++) 
      vec_set(params, params_idx++, 
                     vec_get(mod->backgd_freqs, j));

  /* next parameters are for rate variation */
  if (mod->nratecats > 1 || mod->alpha < 0) {
    if (mod->empirical_rates) 
      for (j = 0; j < mod->nratecats; j++)
        vec_set(params, params_idx++, mod->freqK[j]);
    else                        /* discrete gamma model */
      vec_set(params, params_idx++, mod->alpha);
  }

  /* initialize rate-matrix parameters */
  if (mod->estimate_ratemat)
    tm_rate_params_init_from_model(mod, params, params_idx);
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
  return (tm_get_nbranchlenparams(mod) + tm_get_neqfreqparams(mod) + 
          tm_get_nratevarparams(mod) + tm_get_nratematparams(mod));
}

int tm_get_neqfreqparams(TreeModel *mod) {
  if (mod->estimate_backgd) {
    if (mod->backgd_freqs == NULL) return strlen(mod->rate_matrix->states);
    else return mod->backgd_freqs->size;
  }
  return 0;
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


int tm_get_nleaf(TreeModel *mod) {
  int i,retval=0;
  for (i=0; i<mod->tree->nnodes; i++) {
    TreeNode *n = lst_get_ptr(mod->tree->nodes, i);
    if (n->lchild==NULL) retval++;
  }
  return retval;
}

/** Return number of branch length params */
int tm_get_nbranchlenparams(TreeModel *mod) {
  int retval;
  if (mod->estimate_branchlens == TM_BRANCHLENS_NONE) return 0;
  else if (mod->estimate_branchlens == TM_SCALE_ONLY) {
    if (mod->subtree_root == NULL) return 1;
    else return 2;
  }
  else if (mod->estimate_branchlens == TM_BRANCHLENS_CLOCK) {
    return ((mod->tree->nnodes + 1) / 2 - 1); /* number of ancestral nodes */
  }
  retval = tm_is_reversible(mod->subst_mod) ? mod->tree->nnodes - 2 :
    mod->tree->nnodes - 1;
  if (mod->root_leaf_id != -1) retval--;
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

  tm_free_rmp(mod);             /* necessary because parameter indices
                                   can change */
  tr_prune(&mod->tree, names, TRUE);
  tm_init_rmp(mod);

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
  AltSubstMod *am = (AltSubstMod*)smalloc(sizeof(AltSubstMod));
  am->subst_mod = subst_mod;
  am->backgd_freqs = backgd_freqs;
  am->rate_matrix = rate_matrix;
  return am;
}

void tm_free_alt_subst_mod(AltSubstMod *am) {
  vec_free(am->backgd_freqs);
  mm_free(am->rate_matrix);
  free(am);
}


void tm_variance(TreeModel *mod, MSA *msa, Vector *params, int cat, char *error_fname, int appendToFile) {
  FILE *outfile = fopen_fname(error_fname, appendToFile ? "a" : "w");
  double delta=1.0e-6, origParam, origLike, like1, like2, var, sd;
  int idx;
  
  tm_unpack_params(mod, params, -1);
  origLike = tl_compute_log_likelihood(mod, msa, NULL, cat, NULL);
  for (idx=0; idx < params->size; idx++) {
    origParam = vec_get(params, idx);
    vec_set(params, idx, origParam + 2*delta);
    tm_unpack_params(mod, params, -1);
    like1 = tl_compute_log_likelihood(mod, msa, NULL, cat, NULL);
    vec_set(params, idx, origParam + delta);
    tm_unpack_params(mod, params, -1);
    like2 = tl_compute_log_likelihood(mod, msa, NULL, cat, NULL);
    var = -(delta*delta)/(like1 - 2*like2 + origLike);
    sd = sqrt(var);
    fprintf(outfile, "%f\t%e\t%f\t%f\n", origParam, var, origParam - 1.96*sd, origParam + 1.96*sd);
    vec_set(params, idx, origParam);
  }
  tm_unpack_params(mod, params, -1);
  fclose(outfile);
}
