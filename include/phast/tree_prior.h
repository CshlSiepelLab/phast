/* prior distributions over trees and branch lengths, for use in
   variational inference */

#ifndef TREEPRIOR_H
#define TREEPRIOR_H

#include <stdio.h>
#include <phast/trees.h>
#include <phast/tree_model.h>
#include <phast/matrix.h>
#include <phast/bitset.h>

struct cvdat;

/* types of priors over trees and branch lengths */
enum tree_prior_type {YULE, GAMMA, NONE};

/* metadata for prior; allows for Yule or Gamma model for trees, and
   orthogonally, for a relaxed local clock model (or not) for absolute
   substitution rates along branches.  The clock model uses a
   lognormal prior (with mean 1) for relative substitution rates on
   each branch. The standard deviation (sigma) of this model is
   estimated from the data.  The prior for sigma may be defined as a
   normal prior on log sigma or an exponential prior on sigma (see
   below).  */
typedef struct {
  enum tree_prior_type type; /* YULE or GAMMA for trees and branch lengths,
                           or NONE to ignore tree prior */
  double gamma_shape;  /* hyperparameter: shape of gamma prior (if GAMMA) */
  double gamma_scale;  /* scale of gamma prior; estimated by empirical
                          Bayes from initial tree */
  unsigned int relclock; /* TRUE or FALSE; whether or not relaxed
                            local clock is active.  All relevant
                            parameter are ignored if FALSE */
  double relclock_sig; /* stdev of log normal prior for relaxed clock
                          (estimated).  Actual value of sigma is
                          obtained by applying softplus() to this raw
                          value */
  double relclock_sig_grad; /* gradient of raw sigma (before
                               softplus), for use in gradient
                               ascent */
  double relclock_lsig_mean; /* hyperparameter for mean of log sigma
                                or -1 to use exponential prior
                                instead */
  double relclock_lsig_sd; /* hyperparameter for sd of log sigma if
                              active */
  double relclock_sig_exp_mean; /* hyperparameter for mean of exp
                                   sigma if active (or -1) */
  Vector *nodetimes; /* absolute time of each node; estimated from
                        data. Stored as raw values to which softplus
                        is applied to obtain actual times */
  Vector *nodetimes_grad; /* gradient of raw nodetimes (before softplus) */
  double tau_beta;  /* scale parameter for scaled softplus
                       parameterization of nodetimes; set from branch
                       lengths on initialization */
  BSHash *bs2idx; /* bitset hash used to index internal nodes by sets of
                     descendant leaves; needed for persistance of
                     nodetimes as different trees are sampled */
} TreePrior;

/* shape parameter for Gamma prior; use moderately informative value */
#define GAMMA_SHAPE 2

/* constants for sigma, the standard deviation of the lognormal
   distribution for relative clock rates along branches */
/* initialization for sigma */
#define SIG_INIT 0.5
/* floor for sigma in estimation, to keep from collapsing */
#define SIG_FLOOR 0.4
/* hyperparameter for mean of log sigma */
#define LSIG_MEAN log(0.7)
/* hyperparameter for sd of log sigma */
#define LSIG_SD 0.1
/* hyperparameter for mean of sigma under exponential */
#define SIG_EXP_MEAN 2

TreePrior *tp_new();

double tp_compute_log_prior(TreeModel *mod, struct cvdat *data, Vector *branchgrad);

double tp_prior_noclock(TreeModel *mod, TreePrior *tp, Vector *branchgrad);

void tp_init_nodetimes(TreePrior *tp, TreeModel *mod, List *bs_by_id);

void tp_init_gamma_scale(TreePrior *tp, TreeModel *mod);

double tp_treelen(TreeModel *mod);

#endif 
