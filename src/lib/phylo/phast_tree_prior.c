/* prior distributions over trees and branch lengths, for use in variational inference */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include "phast/trees.h"
#include "phast/misc.h"
#include "phast/nj.h"
#include "phast/lists.h"
#include "phast/stacks.h"
#include "phast/tree_prior.h"

/* floor parameters to keep branch lengths (bl) and internode times (tau) from going to zero */
#define BL_EPS 1.0e-6
#define TAU_EPS 1.0e-6

/* helper functions for scaled softplus parameterization of nodetime deltas (tau) */

/* numerically safe clamp for exponent arguments */
static inline double clamp_exp_arg(double x) {
  if (x > 60.0) return 60.0;
  if (x < -60.0) return -60.0;
  return x;
}

/* scaled softplus:
   ssp(x;β) = (1/β) * log(1 + exp(βx))
   strictly positive, smooth, monotone increasing
*/
static inline double ssp(double x, double beta) {
  double z = clamp_exp_arg(beta * x);
  /* stable evaluation: for large positive z, log1p(exp(z)) ≈ z */
  if (z > 0.0)
    return (z + log1p(exp(-z))) / beta;
  else
    return log1p(exp(z)) / beta;
}

/* derivative wrt x: d/dx ssp(x;β) = sigmoid(βx)
*/
static inline double ssp_prime(double x, double beta) {
  return sigmoid(beta * x);
}

/* inverse of scaled softplus:
   given y ≥ 0, x = (1/β) * log(exp(βy) - 1)
   numerically stable for small y
*/
static inline double ssp_inv(double y, double beta) {
  double by = beta * (y > 0.0 ? y : 0.0);
  /* for very small by, expm1(by) ≈ by */
  if (by < 1e-8)
    return log(fmax(by, 1e-300)) / beta;
  double arg = expm1(by);
  if (arg < 1e-300) arg = 1e-300;
  return log(arg) / beta;
}

/* return a new TreePrior object with default settings; designed to be
   updated later with a particular TreeModel */
TreePrior *tp_new(enum tree_prior_type type, unsigned int relclock) {
  TreePrior *retval = smalloc(sizeof(TreePrior));

  retval->type = type;
  retval->relclock = relclock;

  retval->relclock_sig = inv_softplus(SIG_INIT);
  retval->relclock_sig_grad = 0;

  /* these will be initialized later, after an initial tree is available */
  retval->nodetimes = NULL; 
  retval->nodetimes_grad = NULL;
  retval->bs2idx = NULL;

  /* set these to defaults */
  retval->relclock_lsig_mean = LSIG_MEAN;
  retval->relclock_lsig_sd = LSIG_SD;
  retval->relclock_sig_exp_mean = -1; 
  /* set this instead to use exp prior:  retval->relclock_sig_exp_mean = SIG_EXP_MEAN; */

  retval->gamma_shape = GAMMA_SHAPE;
  retval ->gamma_scale = -1; /* will be set later if needed */
  
  return retval;
}

/* (helper for below) Get (or assign) a nodetime index for this bitset key.
   - used[] marks which indices in nodetimes are already claimed this pass.
   - Returns [0..K-1] where K = nodetimes->size. */
static int tp_get_or_assign_idx(TreePrior *tp, const BSet *key, int *used) {
  void *pv = bs_hash_get(tp->bs2idx, key);
  if (pv != NULL) {
    int idx = *(int*)pv;
    used[idx] = 1;
    return idx;
  }
  /* assign first unused slot */
  int K = tp->nodetimes->size;
  int pick = -1;
  for (int j = 0; j < K; ++j) if (!used[j]) { pick = j; break; }
  if (pick < 0) {
    /* Fallback: reuse 0 if all are marked (shouldn't happen if K==#internals) */
    pick = 0;
  }
  int *stored = (int*)smalloc(sizeof(int));
  *stored = pick;
  bs_hash_put(tp->bs2idx, key, stored);
  used[pick] = 1;
  return pick;
}

/* Huber functions with threshold kappa, for reduced sensitivity to
   outliers compared with standard quadratic */
static inline double huber_rho(double z, double kappa) {
  double a = fabs(z);
  if (a <= kappa) return 0.5 * z * z;
  return kappa * (a - 0.5 * kappa);
}

static inline double huber_psi(double z, double kappa) {  // derivative wrt z
  if (z >=  kappa) return  kappa;
  if (z <= -kappa) return -kappa;
  return z;
}

/* compute log prior for a tree and branch lengths under a simple Yule
   model with relaxed local clock */
double tp_compute_log_prior(TreeModel *mod, struct cvdat *data, Vector *branchgrad) {
  TreePrior *tp = data->treeprior;
  int nn = mod->tree->nnodes;
  int nbranches = nn - 1;

  /* reset grads */
  tp->relclock_sig_grad = 0.0;
  vec_zero(branchgrad);
  
  if (tp->type == GAMMA && tp->gamma_scale == -1)
    tp_init_gamma_scale(tp, mod);

  /* handle simplified case of no clock */
  if (tp->type == NONE && tp->relclock == FALSE)
    return 0;
  else if (tp->relclock == FALSE) /* much simpler case (see below) */
    return tp_prior_noclock(mod, tp, branchgrad);

  /* build per-node bitsets */
  List *bs_by_id = tr_set_leaf_bitsets(mod->tree);       /* size nn, BSet* indexed by node id */
  
  /* Ensure nodetimes allocated once */
  if (tp->nodetimes == NULL)
    tp_init_nodetimes(tp, mod, bs_by_id);
  double beta = tp->tau_beta;
  vec_zero(tp->nodetimes_grad);

  /* Track which nodetime slots are used by keys encountered this pass. */
  int K = tp->nodetimes->size;
  int *used = (int*)scalloc(K, sizeof(int));

  /* Working arrays per node id (for the two-pass tau propagation). */
  Vector *tau_base = vec_new(nn);   /* stores base term: diff/(sig2 * tau) (without softplus and Yule) */
  vec_zero(tau_base);

  /* sigma and constants */
  double raw_sigma = tp->relclock_sig;
  double sp        = softplus(raw_sigma);       // >= 0
  double sig       = SIG_FLOOR + sp;            // >= floor
 
  double sig2 = sig * sig, mu = -0.5 * sig2, lsig = log(sig);

  /* chain terms for backprop to the raw parameter */
  double dsig_draw = sigmoid(raw_sigma);  
  
  /* First pass: per-edge contributions and accumulate timesum */
  double retval = 0.0, timesum = 0.0;

  /* constant term per branch for UC lognormal */
  retval += nbranches * (-lsig - 0.5 * log(2 * M_PI));

  /* Iterate nodes */
  for (int it = 0; it < nn; ++it) {
    TreeNode *n = (TreeNode*)lst_get_ptr(mod->tree->nodes, it);
    if (n->parent == NULL) continue;  /* skip root */

    /* child time (leaf=0; internal=lookup by bitset) */
    double thistime = 0.0;
    int child_idx = -1;
    if (n->lchild != NULL) {  /* internal node */
      BSet *child_bs = (BSet*)lst_get_ptr(bs_by_id, n->id);
      child_idx = tp_get_or_assign_idx(tp, child_bs, used);
      thistime  = vec_get(tp->nodetimes, child_idx);
    }

    /* parent time */
    TreeNode *p = n->parent;
    BSet *parent_bs = (BSet*)lst_get_ptr(bs_by_id, p->id);
    int parent_idx = tp_get_or_assign_idx(tp, parent_bs, used);
    double partime = vec_get(tp->nodetimes, parent_idx);

    /* Edge measures */
    double bl = BL_EPS + n->dparent;

    /* tau via clamped zero-centered softplus on delta */
    double delta = partime - thistime;
    double tau = TAU_EPS + ssp(delta, beta);

    timesum += tau;

    /* UC lognormal */
    double zi    = log(bl / tau);
    double z     = zi - mu;

    /* Huber objective */
    const double KAPPA = 0.5;  
    double rho   = huber_rho(z, KAPPA);
    double psi   = huber_psi(z, KAPPA);

    retval += -log(bl) - (2.0 * rho) / (2.0 * sig2);  // i.e., -rho/sig2

    /* store base tau-term for later propagation */
    /* tau chain uses psi (not z) */
    double base  = psi / (sig2 * tau); 
    vec_set(tau_base, n->id, base);
    
    /* branch-length gradient */
    double dbl_raw = -1.0 / bl - psi / (sig2 * bl);
    vec_set(branchgrad, n->id, dbl_raw);
    
    /* sigma gradient per edge for -log σ - ρ(z)/σ^2 with z = log(bl/tau) - μ, μ = -0.5 σ^2
       d/dσ = -1/σ  - ψ(z)/σ + 2 ρ(z)/σ^3   (consistent with quadratic case when ψ=z, ρ=z^2/2)
    */
    tp->relclock_sig_grad += -1.0/sig - psi/sig + (2.0 * rho)/(sig*sig2);   // since sig*sig2 = σ^3
  }

  /* Yule (lambda integrated out) */
  if (tp->type == YULE)
    retval += -nbranches * log(timesum);
  else if (tp->type == GAMMA)
    retval += (tp->gamma_shape-1.0) * log(timesum) - timesum/tp->gamma_scale -
      tp->gamma_shape * log(tp->gamma_scale) - lgamma(tp->gamma_shape);
  else
    assert(tp->type == NONE); /* do nothing */

  /* sigma prior: depends on selected model */
  if (tp->relclock_sig_exp_mean == -1) { /* gaussian prior on log sigma */
    double lsigdif = lsig - tp->relclock_lsig_mean;
    double lsigvar =  tp->relclock_lsig_sd * tp->relclock_lsig_sd;
    retval += -0.5 * (lsigdif * lsigdif) / lsigvar;
    /* grad w.r.t. sigma: note this is for sigma, not raw param */
    tp->relclock_sig_grad += -lsigdif / lsigvar * (1.0 / sig);
  }
  else { /* exponential prior on sigma */
    retval += -log(tp->relclock_sig_exp_mean) - sig / tp->relclock_sig_exp_mean;
    tp->relclock_sig_grad += -1.0 / tp->relclock_sig_exp_mean; /* again, for sigma not raw */
  }

  assert(isfinite(retval));
  
  if (sig < SIG_FLOOR && tp->relclock_sig_grad < 0.0)
    tp->relclock_sig_grad = 0.0;   /* do not push below floor */

  /* final chain to raw param */
  tp->relclock_sig_grad *= dsig_draw;

  /* Iterate nodes again properly */
  for (int it = 0; it < nn; ++it) {
    TreeNode *n = (TreeNode*)lst_get_ptr(mod->tree->nodes, it);
    if (n->parent == NULL) continue;

    /* retrieve child/parent idx and thistime/partime again (cheap) */
    double thistime = 0.0;
    int child_idx = -1;
    if (n->lchild != NULL) {
      BSet *child_bs = (BSet*)lst_get_ptr(bs_by_id, n->id);
      child_idx = tp_get_or_assign_idx(tp, child_bs, used);   /* will hit existing mapping */
      thistime  = vec_get(tp->nodetimes, child_idx);
    }
    TreeNode *p = n->parent;
    BSet *parent_bs = (BSet*)lst_get_ptr(bs_by_id, p->id);
    int parent_idx = tp_get_or_assign_idx(tp, parent_bs, used);
    double partime = vec_get(tp->nodetimes, parent_idx);

    /* d softplus(delta) */
    double delta = partime - thistime;
    double dsoft = ssp_prime(delta, beta);

    /* tau derivative including Yule/Gamma piece */
    double base = vec_get(tau_base, n->id);     
    double gT = 0;    
    if (tp->type == YULE)
      gT = -(double)nbranches / timesum;
    else if (tp->type == GAMMA)
      gT = (tp->gamma_shape - 1.0) / timesum - 1.0 / tp->gamma_scale;
    else
      assert(tp->type == NONE); /* do nothing */
    
    double tau_d = (base + gT) * dsoft;

    /* accumulate to nodetimes_grad: + to parent, − to child (if internal) */
    vec_set(tp->nodetimes_grad, parent_idx,
            vec_get(tp->nodetimes_grad, parent_idx) + tau_d);
    if (child_idx >= 0)
      vec_set(tp->nodetimes_grad, child_idx,
              vec_get(tp->nodetimes_grad, child_idx) - tau_d);
  }

  /* gauge projection: center nodetimes_grad on its mean */
  double gmean = vec_sum(tp->nodetimes_grad) / tp->nodetimes_grad->size;
  for (int k = 0; k < tp->nodetimes_grad->size; k++) 
    vec_set(tp->nodetimes_grad, k,  vec_get(tp->nodetimes_grad, k) - gmean);
 
  for (int id = 0; id < nn; ++id) {
    BSet *bs = (BSet*)lst_get_ptr(bs_by_id, id);
    bs_free(bs);
  }
  lst_free(bs_by_id);
  vec_free(tau_base);
  sfree(used);

  return retval;
}

/* Special case when the relaxed clock is off:
   Use a tree-length prior on T = sum_b (BL_EPS + n->dparent).
   - YULE:    log p ∝ -B * log T
              d/d bl_i log p = -(B / T)   where B = nbranches
   - GAMMA(k, θ) on T: log p = (k-1)log T - T/θ - k log θ - lgamma(k)
                       d/d bl_i log p = (k-1)/T - 1/θ
   This version approximates Yule’s “sum of waiting times” using total
   tree length T.
*/
double tp_prior_noclock(TreeModel *mod, TreePrior *tp, Vector *branchgrad) {
  const int nn = mod->tree->nnodes;
  const int nbranches = nn - 1;

  double T = tp_treelen(mod);
  if (T < 1e-12) T = 1e-12;  /* avoid NaN */

  double lp = 0.0;
  double dlogp_dT = 0.0;

  if (tp->type == YULE) {
    lp = -(double)nbranches * log(T);
    dlogp_dT = -(double)nbranches / T;
  }
  else if (tp->type == GAMMA) {
    lp = (tp->gamma_shape - 1.0) * log(T) - T / tp->gamma_scale -
      tp->gamma_shape * log(tp->gamma_scale) - lgamma(tp->gamma_shape);
    dlogp_dT = (tp->gamma_shape - 1.0) / T - 1.0 / tp->gamma_scale;
  }

  /* gradient wrt each branch length: dT/d bl_i = 1 */
  vec_set_all(branchgrad, dlogp_dT);
  
  return lp;
}

/* Initialize nodetimes so that, approximately, bl/tau ≈ 1 on every edge.
   Uses bitsets + hash to map each internal node to a stable slot in tp->nodetimes.
   Sets scale parameter beta as a side-effect.
*/
void tp_init_nodetimes(TreePrior *tp, TreeModel *mod, List *bs_by_id) {
  TreeNode *root = mod->tree;
  int nn = root->nnodes;
  int nleaves = (nn + 1) / 2;
  int ninternal = nleaves - 1;

  /* allocate nodetimes vector (fixed size) */
  tp->nodetimes = vec_new(ninternal);
  tp->nodetimes_grad = vec_new(ninternal);
  vec_zero(tp->nodetimes);
  vec_zero(tp->nodetimes_grad);

  if (tp->bs2idx == NULL) tp->bs2idx = bs_hash_new(128);

 /* we need to set the scale appropriately for the softmax on nodetimes.  Choose
     scale-factor beta based on current tree: beta ≈ 3 / ave-branch-length */
  double avebl = tp_treelen(mod) / (nn-1.0);
  double alpha = 0.1;  /* controls offset fraction: can be roughly 0.05–0.2 */
  tp->tau_beta = 1.0 / (alpha * fmax(avebl, 1e-12));
  
  /* Track which nodetime slots we’ve assigned this init */
  int K = tp->nodetimes->size;               /* should equal ninternal */
  int *used = (int*)scalloc(K, sizeof(int));

  /* Raw times per node id (these are the *raw* parameters, not softplussed).
     Leaves start at 0. */
  double *t_raw = (double*)scalloc(nn, sizeof(double));

  /* Postorder traversal to compute t_raw (bottom-up) */
  List *post = tr_postorder(root);
  for (int i = 0; i < nn; ++i) {
    TreeNode *v = (TreeNode*)lst_get_ptr(post, i);

    if (v->lchild == NULL && v->rchild == NULL) 
      continue; /* already zero */

    /* Children contributions */
    TreeNode *L = v->lchild, *R = v->rchild;

    double blL = BL_EPS + L->dparent;  
    double blR = BL_EPS + R->dparent;

    /* target y so that TAU_EPS + ssp(Δ;beta) ≈ bl  =>  Δ ≈ ssp_inv(bl - TAU_EPS; beta) */
    double yL = fmax(blL - TAU_EPS, 1e-12);
    double yR = fmax(blR - TAU_EPS, 1e-12);

    double dL = ssp_inv(yL, tp->tau_beta);
    double dR = ssp_inv(yR, tp->tau_beta);

    double pL = t_raw[L->id] + dL;
    double pR = t_raw[R->id] + dR;

    /* choose parent raw time (avg or max; max is safest to keep Δ≥0 on both edges) */
    double tp_raw = (pL > pR ? pL : pR);

    /* tiny separation so parent > children */
    if (tp_raw < t_raw[L->id] + TAU_EPS) tp_raw = t_raw[L->id] + TAU_EPS;
    if (tp_raw < t_raw[R->id] + TAU_EPS) tp_raw = t_raw[R->id] + TAU_EPS;

    t_raw[v->id] = tp_raw;
  }

  /* Write raw times into the fixed nodetimes vector via bitset mapping */
  for (int i = 0; i < nn; ++i) {
    TreeNode *v = (TreeNode*)lst_get_ptr(post, i);
    if (v->lchild == NULL && v->rchild == NULL) continue;   /* skip leaves */

    BSet *key = (BSet*)lst_get_ptr(bs_by_id, v->id);
    int idx = tp_get_or_assign_idx(tp, key, used);
    vec_set(tp->nodetimes, idx, t_raw[v->id]);
  }

  /* clean up */
  sfree(used);
  sfree(t_raw);
}

/* initialize scale parameter of gamma prior based on starting tree,
   using empirical bayes estimator */
void tp_init_gamma_scale(TreePrior *tp, TreeModel *mod) {
  tp->gamma_scale = tp_treelen(mod) / tp->gamma_shape; /* ensure expected value equals totlen */
}

/* compute total tree length using same floor for branch lengths as
   used in prior computations */
double tp_treelen(TreeModel *mod) {
  double totlen = 0; /* need a version that uses same floor as tp_compute_log_prior */
  for (int i = 0; i < mod->tree->nnodes; i++) {
    TreeNode *n = lst_get_ptr(mod->tree->nodes, i);
    if (n->parent == NULL) continue;
    totlen += BL_EPS + n->dparent;
  }
  return totlen;
}
