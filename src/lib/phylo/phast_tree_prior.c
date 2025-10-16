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

#define RELCLOCK_SIG_INIT 1
#define BL_EPS 1.0e-6
/* #define TAU_EPS 1.0e-6 */ /* TEMPORARY */
#define TAU_EPS 5.0e-6
#define BETA_TAU 0.5
#define DELTA_MIN -8.0

/* return a new TreePrior object with default settings; designed to be
   updated later with a particular TreeModel */
TreePrior *tp_new() {
  TreePrior *retval = smalloc(sizeof(TreePrior));
  retval->relclock_sig = RELCLOCK_SIG_INIT;
  retval->relclock_sig_grad = 0;
  retval->nodetimes = NULL;
  retval->nodetimes_grad = NULL;
  return retval;
}

/* compute log prior for a tree and branch lengths under a simple Yule
   model with relaxed local clock */
double tp_compute_log_prior(TreeModel *mod, struct cvdat *data, Vector *branchgrad) {
  int i;  
  /* double sig = softplus(data->treeprior->relclock_sig); */
  double sig = 0.25 + softplus(data->treeprior->relclock_sig); /* TEMPORARY */
  double sig2 = sig * sig, mu = -0.5 * sig2, lsig = log(sig), twice_sig2 = 2 * sig2;
  double retval = 0, timesum = 0, thistime, partime, bl, tau, zi, diff, diff2, tau_d;
  NodeDist *nd;
  int nbranches = mod->tree->nnodes - 1, nleaves = (mod->tree->nnodes+1)/2;
  Vector *tau_deriv = vec_new(mod->tree->nnodes);

  if (data->treeprior->nodetimes == NULL)
    tp_init_nodetimes(data->treeprior, mod);
  
  data->treeprior->relclock_sig_grad = 0;
  vec_zero(data->treeprior->nodetimes_grad);
  
  List *sorted_nodelist = tp_sort_nodes(mod->tree);
  /* This ordering (from largest to smallest distance from the root)
     will implicitly define the meaning of data->treeprior->nodetimes.
     The first item in this list has the first nodetime, the second
     has the second, etc. */

  assert(data->treeprior->nodetimes->size == mod->tree->nnodes - nleaves);
  /* indexing: first nleaves items in the sorted nodelist (leaves of
     the tree) will not have counterparts in nodetimes because we only
     need to keep track of nodetimes of internal nodes */
  
  /* build mapping from ids of internal nodes to associated indices in
     nodetimes.  Has to accommodate an implicit offset of nleaves in
     sorted_nodelist */
  int *map = smalloc(lst_size(sorted_nodelist) * sizeof(int));
  for (i = 0; i < lst_size(sorted_nodelist); i++) {
    nd = lst_get_ptr(sorted_nodelist, i);
    if (nd->node->lchild == NULL) /* null mapping for leaves */
      map[nd->node->id] = -1;
    else
      map[nd->node->id] = i - nleaves;
  }
    
  /* lognormal clock rate: constant terms per branch */
  retval = nbranches * (-lsig - 0.5 * log (2 * M_PI));
  for (i = 0; i < lst_size(sorted_nodelist); i++) {
    nd = lst_get_ptr(sorted_nodelist, i);
    TreeNode *n = nd->node;

    if (n->parent == NULL)
      continue; /* ignore root here */

    /* node times for this node and its parent */
    if (n->lchild == NULL) /* force leaves to have nodetime of zero */
      thistime = 0;
    else
      thistime = vec_get(data->treeprior->nodetimes, i - nleaves); /* note idx offset */
    partime = vec_get(data->treeprior->nodetimes, map[n->parent->id]);
    
    /* branchlength and corresponding time difference */
    bl = BL_EPS + softplus(n->dparent);   /* avoid zeroes */
    /* tau = TAU_EPS + softplus(partime - thistime);   */ /* TEMPORARY */
    double delta = partime - thistime;
    double delta_c = (delta < DELTA_MIN ? DELTA_MIN : delta);
    tau = TAU_EPS + (1.0/BETA_TAU) * softplus(BETA_TAU * delta_c);    
    timesum += tau;

    /* precompute for below */
    zi = log(bl / tau); /* log ratio */
    diff = zi - mu; /* difference from mean */
    diff2 = diff*diff; /* difference squared */
    
    /* per-branch contribution from local clock */
    retval += -log(bl) - diff2 / twice_sig2 ;

    /* branch-length gradient */
    double bl_deriv = (-1.0/bl - diff / (sig2 * bl)) * sigmoid(n->dparent); 
    /* (sigmoid needed for chain rule of softplus) */
    vec_set(branchgrad, n->id, bl_deriv); 

    /* tau gradient, for use below in nodetimes */  
    vec_set(tau_deriv, i, diff / (sig2 * tau));

    /* contribution to sig gradient */ 
    data->treeprior->relclock_sig_grad += (-1.0/sig + diff2/(sig*sig2) - diff/sig);
  }

  /* log prior contribution from Yule process (rate lambda is
     implicitly integrated out) */
  retval -= nbranches * log(timesum);
  /* (ignores a constant that does not depend on params) */

  /* log prior contribution from sigma */
  retval += -log(2) - 0.5 * sig;

  /* contribution of exponential prior to gradient of sigma. */
  data->treeprior->relclock_sig_grad -= 0.5;
  data->treeprior->relclock_sig_grad *= sigmoid(data->treeprior->relclock_sig);

  /* finalize tau gradient and propagate to node times */
  for (i = 0; i < mod->tree->nnodes; i++) {
    nd = lst_get_ptr(sorted_nodelist, i);
    TreeNode *n = nd->node;
    if (n->parent == NULL) continue;
    thistime = (n->lchild == NULL ? 0 : vec_get(data->treeprior->nodetimes, i - nleaves));
    partime = vec_get(data->treeprior->nodetimes, map[n->parent->id]);
    
    /* contribution to node-time gradients from Yule process */
    /* tau_d = (vec_get(tau_deriv, i) - (double)nbranches/timesum) * sigmoid(partime-thistime); */ /* TEMPORARY */
    double delta = partime - thistime;
    double dsoft = (delta < DELTA_MIN ? 0.0 : sigmoid(BETA_TAU * delta));
    tau_d = (vec_get(tau_deriv, i) - (double)nbranches/timesum) * dsoft;
    /* (sigmoid for softplus) */
    
    /* now propagate to node times: + to parent, - to child */    
    vec_set(data->treeprior->nodetimes_grad, map[n->parent->id],
            vec_get(data->treeprior->nodetimes_grad, map[n->parent->id]) + tau_d);
    if (n->lchild != NULL) /* only if not a leaf */
      vec_set(data->treeprior->nodetimes_grad, i - nleaves,
              vec_get(data->treeprior->nodetimes_grad, i - nleaves) - tau_d);
  }

  /* free sorted nodelist */
  for (i = 0; i < lst_size(sorted_nodelist); i++)
    sfree(lst_get_ptr(sorted_nodelist, i));
  lst_free(sorted_nodelist);
  vec_free(tau_deriv);
  
  sfree(map);

  return retval;
}

/* comparator for sorting in descending order */
static inline
int node_dist_comp(const void* ptr1, const void* ptr2) {
  NodeDist *nd1 = *((NodeDist**)ptr1);
  NodeDist *nd2 = *((NodeDist**)ptr2);
  if (nd1->dist_from_root == nd2->dist_from_root)
    return 0;
  else if (nd2->dist_from_root > nd1->dist_from_root)
    return 1;
  return -1;
}

/* helper to make a new node */
static inline
NodeDist *new_node_dist(TreeNode *n, double dist) {
  NodeDist *nd = smalloc(sizeof(NodeDist));
  nd->node = n;
  nd->dist_from_root = dist;
  return nd;
}

/* sort nodes in tree by descending distance from root.  Leaves are
   forced to come first in the list (with dist INFTY).  Returns a list of
   NodeDist structs */
List *tp_sort_nodes(TreeNode *tree) {
  Stack *s = stk_new_ptr(tree->nnodes);
  List *retval = lst_new_ptr(tree->nnodes);
  NodeDist *nd;
  TreeNode *n;
  double d;
  
  nd = new_node_dist(tree, 0);  
  stk_push_ptr(s, nd);
  while ((nd = stk_pop_ptr(s)) != NULL) {
    n = nd->node;
    d = nd->dist_from_root;

    /* force leaves to front of sorted list */
    if (nd->node->lchild == NULL)
      nd->dist_from_root = INFTY;

    lst_push_ptr(retval, nd); /* put on final list */

    if (n->lchild != NULL) {
      assert(n->rchild != NULL);
      NodeDist *lnd = new_node_dist(n->lchild, n->lchild->dparent + d);
      NodeDist *rnd = new_node_dist(n->rchild, n->rchild->dparent + d);      
      stk_push_ptr(s, lnd);
      stk_push_ptr(s, rnd);      
    }
  }

  /* sort */
  lst_qsort(retval, node_dist_comp);
  stk_free(s);
 
  return retval;
}

/* set up nodetimes parameter vector and initialize based on starting
   tree.  Initialization is such that differences between parent/child
   nodetimes are about equal to starting branch lengths on average */
void tp_init_nodetimes(TreePrior *tp, TreeModel *mod) {
  int nleaves = (mod->tree->nnodes + 1) / 2;
  int ninternal = nleaves - 1;
  tp->nodetimes = vec_new(ninternal);
  tp->nodetimes_grad = vec_new(ninternal);
  vec_zero(tp->nodetimes); vec_zero(tp->nodetimes_grad);
  
  List *sorted = tp_sort_nodes(mod->tree); 
  int N = lst_size(sorted);

  /* temporary array of real heights (time units), aligned to sorted list indices */
  double *h = smalloc(N * sizeof(double));

  /* build mapping from node ids to associated indices in
     sorted list. Simultaneously set h=0 for leaves  */
  int *map = smalloc(lst_size(sorted) * sizeof(int));
  for (int i = 0; i < N; i++) {
    NodeDist *nd = lst_get_ptr(sorted, i);
    TreeNode *n  = nd->node;
    map[n->id] = i;
    if (n->lchild == NULL) 
      h[i] = 0.0;
  }

  /* set internal nodes in single upward pass */
  for (int i = nleaves; i < N; i++) {
    NodeDist *nd = lst_get_ptr(sorted, i);
    TreeNode *v  = nd->node;

    assert(v->lchild != NULL && v->rchild != NULL);
    
    TreeNode *cl = v->lchild, *cr = v->rchild;
    int il = map[cl->id], ir = map[cr->id]; /* indices of children */
    
    /* transformed branch lengths child->parent */ 
    double bl = BL_EPS + softplus(cl->dparent);
    double br = BL_EPS + softplus(cr->dparent);

    /* use ave to define height but maintain minimal separation */
    double tL = h[il] + bl;
    double tR = h[ir] + br;
    double avet = 0.5 * (tL + tR);
    double maxh = (h[il] > h[ir] ? h[il] : h[ir]);
    if (avet < maxh + TAU_EPS) avet = maxh + TAU_EPS;

    h[i] = avet;

    /* back-calculate raw value (before softplus) and store in param vector */
    double y = h[i] - TAU_EPS;
    if (y < 1e-12) y = 1e-12; 
    vec_set(tp->nodetimes, i-nleaves, inv_softplus(y));
  }

  sfree(h);
  sfree(map);
  
  /* free NodeDist* list from tp_sort_nodes */
  for (int i = 0; i < N; i++) sfree(lst_get_ptr(sorted, i));
  lst_free(sorted);
}
