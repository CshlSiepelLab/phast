/* $Id: tree_model.c,v 1.2 2004-06-04 21:56:33 acs Exp $
   Written by Adam Siepel, 2002
   Copyright 2002, Adam Siepel, University of California */

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
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_randist.h>
#include <dgamma.h>
#include <math.h>

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

/* (acs, 3/03 [edited 5/04]) This package is getting so that it needs 
   rewriting.  Some things to do:
   - order, reversibility, and strand symmetry should be considered
   attributes of a substitution model.  The number of rate matrix
   parameters follows from these things.  Simple models such as JC,
   K2P, and HKY may still be treated as special cases.  When adding a
   new model, should be able simply to add this kind of info in
   tree_model.h (also, mapping between name and ordinal value of
   subst_mod).
   - matrix initialization should be rethought (e.g., handling of eq
   freqs).  The mapping from param indices to matrix coords should
   probably be set up at this stage.
   - need a sensible, consistent strategy for matrix scaling 
   - a single, simple (and much more efficient) function should be
   used to set up a matrix from a parameter vector, using the mapping
   from indices to matrix coords, and (if necessary), information on
   reversibility, strand symmetry, etc.
*/

/* internal functions */
double tm_likelihood_wrapper(gsl_vector *params, void *data);

/* tree == NULL implies weight matrix (most other params ignored in
   this case) */
TreeModel *tm_new(TreeNode *tree, MarkovMatrix *rate_matrix, 
                  gsl_vector *backgd_freqs, subst_mod_type subst_mod, 
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
      double interval_size, initalpha = alpha > 0 ? alpha : 1;
      tm->empirical_rates = 1;
      if (nratecats != lst_size(rate_consts) )
        die("ERROR: number of explicitly defined rate constants must equal number of rate categories.\n");
      for (i = 0; i < nratecats; i++) {
        tm->rK[i] = lst_get_dbl(rate_consts, i);
        interval_size = tm->rK[i] - (i > 0 ? tm->rK[i-1] : 0);
        tm->freqK[i] = gsl_ran_gamma_pdf(tm->rK[i], initalpha, 1/initalpha) *
          interval_size; 
        /* init to approx gamma with shape param alpha */
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
  tm->lnL = NULL_LOG_LIKELIHOOD;
  tm->max_samples = -1;
  tm->tree_posteriors = NULL;
  tm->tree_posteriors_marg = NULL;
  tm->use_conditionals = 0;
  tm->category = -1;
  tm->allow_but_penalize_gaps = 0;
  tm->allow_gaps = 1;
  tm->estimate_backgd = 0;
  tm->estimate_branchlens = TM_BRANCHLENS_ALL;
  tm->scale = 1;

  tm_init_rmp(tm);

  return tm;
}

/* re-initialize a tree model with an altered substitution model,
   number of rate categories, alpha, set of rate consts, etc.  */ 
void tm_reinit(TreeModel *tm, subst_mod_type new_subst_mod, 
               int new_nratecats, double new_alpha, 
               List *new_rate_consts) {
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
    if (new_nratecats != lst_size(new_rate_consts))
      die("ERROR: number of explicitly defined rate constants must equal number of rate categories.\n");
    for (i = 0; i < new_nratecats; i++) {
      tm->rK[i] = lst_get_dbl(new_rate_consts, i);
      tm->freqK[i] = exp(-tm->rK[i]); /* init to approx exponential
                                         distrib with mean 1 */
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
    free(tm->P);
    free(tm->rK);
    free(tm->freqK);
    free_tree(tm->tree);
  }
  if (tm->rate_matrix != NULL) mm_free(tm->rate_matrix);
  if (tm->backgd_freqs != NULL) gsl_vector_free(tm->backgd_freqs);
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
  gsl_vector *backgd = NULL, *rate_weights = NULL;
  gsl_matrix *rmat = NULL;
  MarkovMatrix *M = NULL;
  TreeNode *tree = NULL;
  int size = 0, order = 0, nratecats = -1;
  double alpha = 0;
  TreeModel *retval;
  subst_mod_type subst_mod = UNDEF_MOD;
  int i, j;
  List *rate_consts = NULL;

  while (fscanf(f, "%s", tag) != EOF) {
    if (!strcmp(tag, ALPHABET_TAG)) {
      str_readline(tmpstr, f);

      j = 0;
      for (i = 0; i < tmpstr->length; i++) 
        if (isalnum(tmpstr->chars[i])) alphabet[j++] = tmpstr->chars[i];
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
      subst_mod = tm_get_subst_mod_type(tmpstr); 
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
      backgd = gsl_vector_alloc(size);
      gsl_vector_fscanf(f, backgd);
    }
    else if (!strcmp(tag, RATE_MATRIX_TAG)) {
      if (size == 0) {
        fprintf(stderr, "ERROR: ALPHABET line must precede BACKGROUND and RATE_MAT in tree model file.\n");
        exit(1);
      }
      rmat = gsl_matrix_alloc(size, size);
      gsl_matrix_fscanf(f, rmat);
    }
    else if (!strcmp(tag, TREE_TAG)) {
      str_readline(tmpstr, f);
      tree = parse_nh_from_string(tmpstr);      
    }
    else if (strcmp(tag, LNL_TAG) == 0) 
      str_readline(tmpstr, f);  /* discard */
    else if (!strcmp(tag, RATE_CONSTS_TAG)) {
      gsl_vector *tmpvect;
      if (nratecats < 0) 
        die("ERROR: NRATECATS must precede RATE_CONSTS in tree model file.\n");
      /* easiest to use gsl_vector_fscanf and convert */
      tmpvect = gsl_vector_alloc(nratecats); 
      gsl_vector_fscanf(f, tmpvect);
      rate_consts = lst_new_dbl(nratecats);
      for (i = 0; i < nratecats; i++) 
        lst_push_dbl(rate_consts, gsl_vector_get(tmpvect, i));
      gsl_vector_free(tmpvect);                     
    }
    else if (!strcmp(tag, RATE_WEIGHTS_TAG)) {
      gsl_vector *rate_weights;
      if (nratecats < 0) 
        die("ERROR: NRATECATS must precede RATE_WEIGHTS in tree model file.\n");
      rate_weights = gsl_vector_alloc(nratecats);
      gsl_vector_fscanf(f, rate_weights);
    }
    else {
      fprintf(stderr, "Error: unrecognized tag in model file (\"%s\").\n", 
	      tag);
      exit(1);
    }
  }

  if (rmat != NULL) 
    M = mm_new_from_matrix(rmat, alphabet, CONTINUOUS);
  else 
    M = mm_new(size, alphabet, CONTINUOUS);

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
      retval->freqK[i] = gsl_vector_get(rate_weights, i);
    normalize_probs(retval->freqK, retval->nratecats);
  }

  str_free(tmpstr);
  if (rate_consts != NULL) lst_free(rate_consts);
  if (rate_weights != NULL) gsl_vector_free(rate_weights);

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
    fprintf(F, "%f ", gsl_vector_get(tm->backgd_freqs, i));
  fprintf(F, "\n");

  if (tm->rate_matrix != NULL) {
    fprintf(F, "%s\n", RATE_MATRIX_TAG);
    gsl_matrix_pretty_print(F, tm->rate_matrix->matrix);
  }

  if (tm->tree != NULL) {
    fprintf(F, "%s ", TREE_TAG);
    print_tree(F, tm->tree, 1);
  }
}

TreeModel *tm_create_copy(TreeModel *src) {
  TreeModel *retval;
  gsl_vector *freqs = gsl_vector_alloc(src->backgd_freqs->size);
  int i;

  gsl_vector_memcpy(freqs, src->backgd_freqs);
  retval = tm_new(src->tree != NULL ? tr_create_copy(src->tree) : NULL, 
                  mm_create_copy(src->rate_matrix), freqs, 
                  src->subst_mod, NULL, 
                  src->nratecats, src->alpha, NULL, 
                  src->root_leaf_id);
  retval->msa = src->msa;
  retval->lnL = src->lnL;
  retval->use_conditionals = src->use_conditionals;
  retval->max_samples = src->max_samples;
  retval->category = src->category;
  retval->allow_gaps = src->allow_gaps;
  retval->allow_but_penalize_gaps = src->allow_but_penalize_gaps;
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
  double scaling_const = -1, tmp;

  if (tm->estimate_branchlens != TM_SCALE_ONLY && tm->scale != 1) 
    tm->scale = 1;
                                /* be sure scale factor has an effect
                                   only if estimating scale */

  /* need to compute a matrix scaling constant from the equilibrium
     freqs, in this case (see below) */
  if (tm->subst_mod == F81) {   
    for (i = 0, tmp = 0; i < tm->rate_matrix->size; i++)
      tmp += gsl_vector_get(tm->backgd_freqs, i) *
        gsl_vector_get(tm->backgd_freqs, i);
    scaling_const = 1.0/(1 - tmp);
  }

  for (i = 0; i < tm->tree->nnodes; i++) {
    TreeNode *n = lst_get_ptr(tm->tree->nodes, i);
    if (n->parent == NULL) continue;
    for (j = 0; j < tm->nratecats; j++) {
      if (tm->P[i][j] == NULL)
        tm->P[i][j] = mm_new(tm->rate_matrix->size, tm->rate_matrix->states, 
                             DISCRETE);

      /* for simple models, full matrix exponentiation is not necessary */
      if (tm->subst_mod == JC69)
        tm_set_probs_JC69(tm, tm->P[i][j], 
                          n->dparent * tm->scale * tm->rK[j]);
      else if (tm->subst_mod == F81)
        tm_set_probs_F81(tm, tm->P[i][j], scaling_const, 
                         n->dparent * tm->scale * tm->rK[j]);
      else 
        mm_exp(tm->P[i][j], tm->rate_matrix, 
               n->dparent * tm->scale * tm->rK[j]);
    }
  }
}

/* scale evolutionary rate by const factor (affects branch lengths
   only) */
void tm_scale(TreeModel *tm, double scale_const, int reset_subst_mats) {
  TreeNode *n;
  Stack *s;

  if (tm->tree == NULL) return;

  s = stk_new_ptr(tm->tree->nnodes);
  stk_push_ptr(s, tm->tree);
  while ((n = stk_pop_ptr(s)) != NULL) {
    n->dparent *= scale_const;
    if (n->lchild != NULL) 
      stk_push_ptr(s, n->lchild);
    if (n->rchild != NULL) 
      stk_push_ptr(s, n->rchild);
  }
  stk_free(s);
  if (reset_subst_mats) tm_set_subst_matrices(tm);
}

/* Generates an alignment according to set of Tree Models and a
   Markov matrix defing how to transition among them.  TreeModels must
   appear in same order as the states of the Markov matrix . */
MSA *tm_generate_msa(int ncolumns, MarkovMatrix *classmat, 
                     TreeModel **classmods, int *labels) {

  int i, class, nseqs, col, ntreenodes;
  MSA *msa;
  Stack *stack;
  char *newchar;
  char **names, **seqs;

  int nclasses = classmat->size;
  assert(nclasses > 0);

  /* obtain number of sequences from tree models; ensure all have same
     number */
  ntreenodes = classmods[0]->tree->nnodes; 
                                /* all should be the same (binary trees) */
  stack = stk_new_ptr(ntreenodes);
  nseqs = -1;
  for (i = 0; i < classmat->size; i++) {
    /* count leaves in tree */
    int num = (classmods[i]->tree->nnodes + 1) / 2;

    assert(classmods[i]->nratecats == 1); /* assuming no rate variation */

    if (nseqs == -1) 
      nseqs = num;
    else if (nseqs != num) {
      fprintf(stderr, "ERROR in msa_generate: model #%d has %d taxa, while a previous model had %d taxa.\n", i+1, num, nseqs);
      assert(0);
    }
  }


  /* create new MSA */
  names = (char**)smalloc(nseqs * sizeof(char*));
  seqs = (char**)smalloc(nseqs * sizeof(char*));
  for (i = 0; i < nseqs; i++) {
    names[i] = (char*)smalloc(MAX_NAME_LEN * sizeof(char));
    seqs[i] = (char*)smalloc((ncolumns + 1) * sizeof(char));
    sprintf(names[i], "seq%d", i + 1);
  }
  msa = msa_new(seqs, names, nseqs, ncolumns, 
                classmods[0]->rate_matrix->states);

  /* generate sequences, column by column */
  class = 0;
  srand(time(NULL));
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

      if (l == NULL) {
        int seqno = atoi(n->name) - 1;
        msa->seqs[seqno][col] = newchar[n->id];
      }
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
    labels[col] = class;
    class = mm_sample_state(classmat, class);
  }

  return msa;
}

/* Given an MSA, a tree topology, and a substitution model, fit a tree
   model using a multidimensional optimization algorithm (BFGS).
   TreeModel 'mod' must already be allocated, and initialized with
   desired tree topology, substitution model, and (if appropriate)
   background frequencies.  The vector 'params' should define the
   initial values for the optimization procedure.  Fuction returns 0
   on success, 1 on failure.  */  
int tm_fit(TreeModel *mod, MSA *msa, gsl_vector *params, int cat, 
           int max_samples, opt_precision_type precision, FILE *logf) {
  double ll;
  gsl_vector *lower_bounds, *upper_bounds;
  int retval = 0;

  if (msa->ss == NULL) {
    assert(msa->seqs != NULL);
    ss_from_msas(msa, mod->order+1, 0, NULL, NULL, NULL, -1);
  }

  if (mod->backgd_freqs == NULL) { 
    mod->backgd_freqs = gsl_vector_calloc(int_pow(strlen(msa->alphabet), 
                                                  mod->order+1));
    if (mod->subst_mod == JC69 || mod->subst_mod == K80)
      gsl_vector_set_all(mod->backgd_freqs, 1.0/mod->backgd_freqs->size);
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
  mod->max_samples = max_samples;

  /* most params have lower bound of zero and no upper bound */
  lower_bounds = gsl_vector_calloc(params->size);
  upper_bounds = NULL;

  /* however, in this case we don't want the eq freqs to go to zero */
  if (mod->estimate_backgd) {
    int i, offset = tm_get_nbranchlenparams(mod);
    for (i = 0; i < mod->backgd_freqs->size; i++)
      gsl_vector_set(lower_bounds, i + offset, 0.001);
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
  }

  if (mod->estimate_branchlens == TM_SCALE_ONLY) { 
                                /* also, do final scaling of tree if
                                   estimating scale only */
    tm_scale(mod, mod->scale, 0);
    mod->scale = 1;
  }

  gsl_vector_free(lower_bounds);

  if (retval != 0) 
    fprintf(stderr, "WARNING: BFGS algorithm reached its maximum number of iterations.\n");

  return retval;
}


/* Wrapper for computation of likelihood, for use by nr_optimize */
double tm_likelihood_wrapper(gsl_vector *params, void *data) {
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
void tm_unpack_params(TreeModel *mod, gsl_vector *params, int idx_offset) {
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
      double p = gsl_vector_get(params, i);
      if (p < 0 && abs(p) < TM_IMAG_EPS) /* consider close enough to 0 */
        gsl_vector_set(params, i, p=0);
      if (p < 0) die("ERROR: parameter %d has become negative (%f).\n", i, p);
      if (!finite(p)) die("ERROR: parameter %d is no longer finite (%f).\n", 
                          i, p);
    }
    i = 0;
  }
  else 
    i = idx_offset;

  if (mod->estimate_branchlens == TM_SCALE_ONLY)
    mod->scale = gsl_vector_get(params, i++);
  else if (mod->estimate_branchlens == TM_BRANCHLENS_ALL) {
    /* first nnodes-2 elements define branch lengths */
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
          n->dparent = gsl_vector_get(params, 0)/2;
          if (!assigned) {
            i++;     /* only increment the first time */
            assigned = 1;
          }
        }
        else if (n->id == mod->root_leaf_id) 
          n->dparent = 0;
        else 
          n->dparent = gsl_vector_get(params, i++);
      }
    }
  }

  /* if estimating backgd, next params are backgd freqs */
  if (mod->estimate_backgd) {
    double sum = 0;
    for (j = 0; j < mod->backgd_freqs->size; j++) 
      sum += gsl_vector_get(params, i+j);
    for (j = 0; j < mod->backgd_freqs->size; j++)
      gsl_vector_set(mod->backgd_freqs, j, gsl_vector_get(params, i++)/sum);
  }

  /* next parameters are for rate variation */
  if (mod->nratecats > 1) {
    if (mod->empirical_rates) {
      for (j = 0; j < mod->nratecats; j++)
        mod->freqK[j] = gsl_vector_get(params, i++);
                                /* assume normalized */
    }
    else {                      /* discrete gamma model */
      mod->alpha = gsl_vector_get(params, i++);
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

  /* redefine substitution matrix */
  tm_set_rate_matrix(mod, params, i);

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
*/ 
double tm_scale_rate_matrix(TreeModel *mod) {
  double scale = 0;
  int i, j;
  for (i = 0; i < mod->rate_matrix->size; i++) {
    double rowsum = 0;
    for (j = 0; j < mod->rate_matrix->size; j++) {
      if (j != i) rowsum += mm_get(mod->rate_matrix, i, j);
    }
    scale += (gsl_vector_get(mod->backgd_freqs, i) * rowsum);
  }

  mm_scale(mod->rate_matrix, (mod->order + 1)/scale);
  return scale/(mod->order + 1);
}

/* scale a parameter vector according to a specified rate matrix scale
   factor.  Branch length params are multiplied by the specified factor
   and rate matrix params by its inverse */
void tm_scale_params(TreeModel *mod, gsl_vector *params, double scale_factor) {
  int i;
  int nbl = tm_get_nbranchlenparams(mod);
  int nrm = tm_get_nratematparams(mod);
  int np = tm_get_nparams(mod);     /* we won't assume anything about the
                                       number of params between the
                                       branchlen and ratemat params */
  for (i = 0; i < nbl; i++)
    gsl_vector_set(params, i, gsl_vector_get(params, i) * scale_factor);
  for (i = np - nrm; i < np; i++) 
    gsl_vector_set(params, i, gsl_vector_get(params, i) / scale_factor);
}

/* initializes all branch lengths to designated constant; initializes
   rate-matrix params using specified kappa; kappa ignored for JC69;
   REV and UNREST initialized as if HKY.  Initializes alpha as
   specified, if dgamma.  In the case of empirical rates, uniform
   weights are used for initialization. */
gsl_vector *tm_params_init(TreeModel *mod, double branchlen, double kappa,
                           double alpha) {
  int nparams = tm_get_nparams(mod);
  gsl_vector *params = gsl_vector_alloc(nparams);
  int i, nbranches, params_idx;

  /* initialize branch-length parameters */
  nbranches = tm_get_nbranchlenparams(mod);
  
  for (i = 0; i < nbranches; i++)
    gsl_vector_set(params, i, branchlen);

  params_idx = nbranches;

  if (mod->estimate_backgd) {
    if (mod->backgd_freqs == NULL) {
      double alph_size = strlen(mod->rate_matrix->states);
      for (i = 0; i < alph_size; i++)
        gsl_vector_set(params, params_idx++, 1.0/alph_size);
    }
    else
      for (i = 0; i < mod->backgd_freqs->size; i++)
        gsl_vector_set(params, params_idx++, gsl_vector_get(mod->backgd_freqs, i));
  }

  if (mod->nratecats > 1) {
    if (mod->empirical_rates)
      for (i = 0; i < mod->nratecats; i++) 
        gsl_vector_set(params, params_idx++, mod->freqK[i]);
    else
      gsl_vector_set(params, params_idx++, alpha);
  }

  /* initialize rate-matrix parameters */
  tm_rate_params_init(mod, params, params_idx, kappa);

  return params;
}

/* initializes branch lengths and rate matrix parameters
   randomly; can be used multiple times to ensure the MLE is real */
gsl_vector *tm_params_init_random(TreeModel *mod) {
  int i, params_idx = 0;
  gsl_vector *params = gsl_vector_alloc(tm_get_nparams(mod));
  int nbranches = tm_get_nbranchlenparams(mod);
  int nratematparams = tm_get_nratematparams(mod);

  assert(!mod->estimate_backgd && 
         mod->estimate_branchlens == TM_BRANCHLENS_ALL);

  srandom(time(NULL));
  
  for (i = 0; i < nbranches; i++)
    gsl_vector_set(params, params_idx++, 
                   0.01 + random() * (0.5 - 0.01) / RAND_MAX);
                                /* we'll use the interval from 0.01 to 0.5 */

  if (mod->nratecats > 1) {
    if (mod->empirical_rates) { /* empirical rate model (category weights) */
      double val, sum = 0;
      for (i = 0; i < mod->nratecats; i++) {
        gsl_vector_set(params, params_idx + i, 
                       val = 0.1 + random() * (1 - 0.5) / RAND_MAX);
                                /* we'll use the interval from 0.1 to 1 */
        sum += val;
      }
      for (i = 0; i < mod->nratecats; i++) /* normalize */
        gsl_vector_set(params, params_idx + i, 
                       gsl_vector_get(params, params_idx + i) / sum);
      params_idx += mod->nratecats;
    }

    else                        /* discrete gamma (alpha) */                       
      gsl_vector_set(params, params_idx++, 
                     0.5 + random() * (10 - 0.5) / RAND_MAX);
                                /* we'll use the interval from 0.5 to 10 */
  }

  for (i = 0; i < nratematparams; i++) 
    gsl_vector_set(params, params_idx++, 
                   0.1 + random() * (5 - 0.1) / RAND_MAX);
                                /* we'll use the interval from 0.1 to 5 */

  return params;
}


/* Functions to initialize a parameter vector from an existing tree model */
gsl_vector *tm_params_init_from_model(TreeModel *mod) {
  int nparams = tm_get_nparams(mod);
  gsl_vector *params = gsl_vector_alloc(nparams);
  int params_idx = 0, nodeidx, j;
  List *traversal;
  TreeNode *n;

  if (mod->estimate_branchlens == TM_SCALE_ONLY) 
    gsl_vector_set(params, params_idx++, mod->scale); 
  else if (mod->estimate_branchlens == TM_BRANCHLENS_ALL) {
    traversal = tr_preorder(mod->tree);
    for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
      n = lst_get_ptr(traversal, nodeidx);
      if (n->parent == NULL) continue;

      /* Note: if the model is reversible, then the first parameter is
         the sum of the lengths of the two branches from the root */
      if (n == mod->tree->lchild && tm_is_reversible(mod->subst_mod))
        gsl_vector_set(params, params_idx++, n->dparent*2);
      else if (n != mod->tree->rchild || !tm_is_reversible(mod->subst_mod))
        gsl_vector_set(params, params_idx++, n->dparent);
    }
  }

  /* if estimating backgd, next params are backgd freqs */
  if (mod->estimate_backgd) 
    for (j = 0; j < mod->backgd_freqs->size; j++) 
      gsl_vector_set(params, params_idx++, 
                     gsl_vector_get(mod->backgd_freqs, j));

  /* next parameters are for rate variation */
  if (mod->nratecats > 1 || mod->alpha < 0) {
    if (mod->empirical_rates) 
      for (j = 0; j < mod->nratecats; j++)
        gsl_vector_set(params, params_idx++, mod->freqK[j]);
    else                        /* discrete gamma model */
      gsl_vector_set(params, params_idx++, mod->alpha);
  }

  /* initialize rate-matrix parameters */
  tm_rate_params_init_from_model(mod, params, params_idx);

  return params;
}

/* Given a codon model, create and return the induced amino acid model */
TreeModel *tm_induced_aa(TreeModel *codon_mod) {
  TreeModel *retval;
  char *codon_to_aa = get_codon_mapping(codon_mod->rate_matrix->states);
  gsl_vector *aa_freqs = gsl_vector_calloc(strlen(AA_ALPHABET));
  MarkovMatrix *aa_mat = mm_new(strlen(AA_ALPHABET), AA_ALPHABET, CONTINUOUS);
  int i, j;
  int nstates = codon_mod->rate_matrix->size;

  assert(codon_mod->order == 2);

  /* compute induced equilibrium freqs */
  for (i = 0; i < nstates; i++) {
    int aa_state_i = aa_mat->inv_states[(int)codon_to_aa[i]]; 
                                /* state number corresponding to amino
                                   acid which corresponds to codon i */
    gsl_vector_set(aa_freqs, aa_state_i, gsl_vector_get(aa_freqs, aa_state_i) + 
                   gsl_vector_get(codon_mod->backgd_freqs, i));
  }

  /* compute induced matrix */
  for (i = 0; i < nstates; i++) {
    int aa_state_i = aa_mat->inv_states[(int)codon_to_aa[i]]; 
    double cond_freq = safediv(gsl_vector_get(codon_mod->backgd_freqs, i),
                               gsl_vector_get(aa_freqs, aa_state_i)); 
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

/* number of branch length params */
int tm_get_nbranchlenparams(TreeModel *mod) {
  int retval;
  if (mod->estimate_branchlens == TM_BRANCHLENS_NONE) return 0;
  else if (mod->estimate_branchlens == TM_SCALE_ONLY) return 1;
  retval = tm_is_reversible(mod->subst_mod) ? mod->tree->nnodes - 2 :
    mod->tree->nnodes - 1;
  if (mod->root_leaf_id != -1) retval--;
  return retval;
}

