/* $Id: tree_model.h,v 1.6 2004-08-05 07:15:04 acs Exp $
   Written by Adam Siepel, 2002
   Copyright 2002, Adam Siepel, University of California */

#ifndef TREE_MODEL_H
#define TREE_MODEL_H

#include <gsl/gsl_vector.h>
#include <markov_matrix.h> 
#include <trees.h> 
#include <stringsplus.h> 
#include <msa.h> 
#include <misc.h>
#include <numerical_opt.h>
#include <assert.h>
#include <subst_mods.h>

#define MAX_ALPH_SIZE 100

#define TM_IMAG_EPS 1e-6           
                                /* in cases where a function of
                                   complex numbers is supposed to be
                                   real, we'll consider it close
                                   enough if the imaginary component is
                                   smaller in magnitude than this
                                   threshold  */

/* convergence thresholds for EM, expressed as portions of function
   value: (f(xold) - f(xnew)) / f(xnew) */
#define TM_EM_CONV_HIGH 1e-8
#define TM_EM_CONV_MED 5e-7
#define TM_EM_CONV_LOW 1e-5
#define TM_EM_CONV_CRUDE 1e-4

#define TM_EM_CONV(P) ( (P) == OPT_HIGH_PREC ? TM_EM_CONV_HIGH : ( (P) == OPT_MED_PREC ? TM_EM_CONV_MED : ( (P) == OPT_CRUDE_PREC ? TM_EM_CONV_CRUDE : TM_EM_CONV_LOW) ))

typedef enum { 
  TM_BRANCHLENS_ALL, 
  TM_SCALE_ONLY, 
  TM_BRANCHLENS_NONE 
} blen_estim_type;

struct tp_struct;

struct tm_struct {
  TreeNode *tree;
  gsl_vector *backgd_freqs;
  MarkovMatrix *rate_matrix;
  subst_mod_type subst_mod;
  int *msa_seq_idx;             /**< (optional) Mapping from leaf
                                   indices to sequence indices in a
                                   given MSA; used in likelihood
                                   computations */
  MSA *msa;                     /* (optional) MSA from which TreeModel
                                   was estimated; for use in tm_fit */
  int category;                 /* (optional) site category in MSA or -1 */
  int order;
  double alpha;                 /* for gamma or discrete-gamma rate
                                   variation; set to 0 for constant rate */
  int nratecats;                /* number of rate categories, discrete gamma */
  double lnL;                   /* log likelihood from training */
  struct tp_struct *tree_posteriors; 
                                /* (optional) associated
                                   TreePosteriors object */
  int use_conditionals;         /* (optional) compute likelihood using
                                   conditional probabilities at each
                                   column; only relevant when order >
                                   0 */
  MarkovMatrix ***P;            /* probability matrices for edges,
                                   indexed by node id and rate category */
  double *rK, *freqK;           /* rate constants and freqs */
  List **rate_matrix_param_row, **rate_matrix_param_col;
  int root_leaf_id;             /* indicates id of leaf node to be
                                   interpretted as the root.  Must be
                                   a child of the actual root of the
                                   tree.  Ordinarily will have a null
                                   value (-1), but if non-negati, the
                                   distance to its parent will be
                                   constrained to be zero.  */
  int allow_gaps;
  int allow_but_penalize_gaps;
  int min_informative;          /* minimum number of informative bases
                                   for non-zero probability */
  int estimate_backgd;          /* estimate backgd freqs as free
                                   parameters in the optimization */
  blen_estim_type  estimate_branchlens; 
                                /* if value is TM_BRANCHLENS_ALL, then
                                   estimate all branch lengths; if
                                   value is TM_SCALE_ONLY, then
                                   estimate only a scale factor; and
                                   if value is TM_BRANCHLENS_NONE, do
                                   not estimate branch lengths */
  double scale;                 /* current scale factor */
  int empirical_rates;          /* indicates "empirical"
                                   (nonparameteric) model for
                                   rate-variation */
};

typedef struct tm_struct TreeModel;


TreeModel *tm_new(TreeNode *tree, MarkovMatrix *rate_matrix, 
                  gsl_vector *backgd_freqs, subst_mod_type subst_mod, 
                  char *alphabet, int nratecats, double alpha,
                  List *rate_consts, int root_leaf_id);

void tm_reinit(TreeModel *tm, subst_mod_type subst_mod,
               int new_nratecats, double new_alpha, 
               List *new_rate_consts, List *new_rate_weights);

TreeModel *tm_new_from_file(FILE *F);

void tm_init_rmp(TreeModel *tm);

void tm_free_rmp(TreeModel *tm);

void tm_free(TreeModel *tm);

void tm_print(FILE *F, TreeModel *tm);

void tm_cpy(TreeModel *dest, TreeModel *src);

TreeModel *tm_create_copy(TreeModel *src);

void tm_set_subst_matrices(TreeModel *tm);

void tm_scale(TreeModel *tm, double scale_const, int reset_subst_mats);

int tm_fit(TreeModel *mod, MSA *msa, gsl_vector *params, int cat, 
           opt_precision_type precision, FILE *logf);

void tm_unpack_params(TreeModel *mod, gsl_vector *params, int idx_offset);

double tm_scale_rate_matrix(TreeModel *mod);

void tm_scale_params(TreeModel *mod, gsl_vector *params, double scale_factor);

gsl_vector *tm_params_init(TreeModel *mod, double branchlen, double kappa,
                           double alpha);

gsl_vector *tm_params_init_random(TreeModel *mod);

gsl_vector *tm_params_init_from_model(TreeModel *mod);

int tm_get_nparams(TreeModel *mod);

int tm_get_neqfreqparams(TreeModel *mod);

int tm_get_nratevarparams(TreeModel *mod);

int tm_is_reversible(int subst_mod);

int tm_get_nbranchlenparams(TreeModel *mod);

MSA *tm_generate_msa(int ncolumns, MarkovMatrix *classmat, 
                     TreeModel **classmods, int *labels);

TreeModel *tm_induced_aa(TreeModel *codon_mod);

void tm_build_seq_idx(TreeModel *mod, MSA *msa);

void tm_prune(TreeModel *mod, MSA *msa, int warn);

#endif
