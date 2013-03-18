/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/


/* "birth-death" phylo-HMM -- nonconserved state, fully conserved
   state, one state per birth event and one state per death event per
   branch of tree */

#include <bd_phylo_hmm.h>
#include <lists.h>
#include <sufficient_stats.h>
#include <numerical_opt.h>

/* create a new birth-death phylo-HMM based on parameter values */
BDPhyloHmm *bd_new(TreeModel *source_mod, double rho, double mu, 
                   double nu, double phi, double alpha_c, double beta_c, 
                   double tau_c, double alpha_n, double beta_n, 
                   double tau_n, int estim_gamma, int estim_omega, 
                   int estim_phi) {
  int i, state;
  int nleaves = (source_mod->tree->nnodes + 1) / 2; 
  int nstates = 2 + 2 * (2*nleaves - 3); /* number of states */
  HMM *hmm;
  TreeModel **models = smalloc(nstates * sizeof(void*));
  TreeNode *n;
  CategoryMap *cm;
  List *inside = lst_new_ptr(source_mod->tree->nnodes),
    *outside = lst_new_ptr(source_mod->tree->nnodes);
  BDPhyloHmm *bdphmm;
  double *alpha = smalloc(source_mod->tree->nnodes * sizeof(double)), 
    *beta = smalloc(source_mod->tree->nnodes * sizeof(double)), 
    *tau = smalloc(source_mod->tree->nnodes * sizeof(double));
  if (rho <= 0 || rho >= 1)
    die("ERROR bd_new: rho (%e) out of bound (0, 1)\n", rho);
  if (mu <= 0 || 2*mu >=1) 
    die("ERROR bd_new: mu (%e) out of bound (0, 0.5)\n", mu);
  if (nu <= 0 || nu >= 1)
    die("ERROR bd_new: nu (%e) out of bound (0, 1)\n", nu);

  bdphmm = smalloc(sizeof(BDPhyloHmm));
  bdphmm->state_to_branch = smalloc(nstates * sizeof(int));
  bdphmm->branch_to_state_birth = smalloc(source_mod->tree->nnodes * sizeof(int));
  bdphmm->branch_to_state_death = smalloc(source_mod->tree->nnodes * sizeof(int));
  bdphmm->rho = rho;
  bdphmm->mu = mu;
  bdphmm->nu = nu;
  bdphmm->phi = phi;
  bdphmm->estim_gamma = estim_gamma;
  bdphmm->estim_omega = estim_omega;
  bdphmm->estim_phi = estim_phi;

  /* set up mappings between states and branches.  Branches are
     indexed by ids of child nodes; "birth" means functional element
     arose on branch, "death" means element was lost on branch. */
  state = 0;
  for (i = 0; i < source_mod->tree->nnodes; i++) {
    n = lst_get_ptr(source_mod->tree->nodes, i);
    bdphmm->state_to_branch[state] = n->id;
    bdphmm->branch_to_state_death[n->id] = state;
    state++;
  }
  for (i = 0; i < source_mod->tree->nnodes; i++) {
    n = lst_get_ptr(source_mod->tree->nodes, i);

    if (n->parent == source_mod->tree) {
      bdphmm->branch_to_state_birth[n->id] = -1;
      continue; 
      /* branches beneath root ignored -- events can't be
         distinguished from death events on opposite branch */
    }

    bdphmm->state_to_branch[state] = n->id;
    bdphmm->branch_to_state_birth[n->id] = state;
    state++;
  }
  if  (state != nstates)
    die("ERROR bd_new state (%i) != nstates (%i)\n", state, nstates);
  /* note: this is structured so that state 0 is the nonconserved
     model (like element "died" prior to root) and state nnodes is the
     fully conserved model (like element was "born" prior to root) */

  /* now make copies of models and scale branches as needed */  
  models[0] = source_mod;
  for (state = 1; state < nstates; state++) { 
    /* skip state == 0 */
    models[state] = tm_create_copy(source_mod);
    n = lst_get_ptr(models[state]->tree->nodes, bdphmm->state_to_branch[state]);
    tr_partition_nodes(models[state]->tree, n, inside, outside);

    if (state < models[state]->tree->nnodes)  /* death */
      for (i = 0; i < lst_size(outside); i++) 
        ((TreeNode*)lst_get_ptr(outside, i))->dparent *= rho;
    else 			/* birth */
      for (i = 0; i < lst_size(inside); i++) 
        ((TreeNode*)lst_get_ptr(inside, i))->dparent *= rho;

    tm_reinit(models[state], models[state]->subst_mod, 
              models[state]->nratecats, models[state]->alpha, NULL, NULL);
  }	  

  /* Require informative sites for all but state 0 */
  for (state = 1; state < nstates; state++) 
    models[state]->inform_reqd = TRUE;

  /* set up HMM transitions */
  hmm = hmm_new_nstates(nstates, TRUE, FALSE);

  /* build category map */
  cm = cm_create_trivial(nstates - 1, NULL);
  for (state = 0; state < nstates; state++) {
    String *name = lst_get_ptr(cm->ranges[state]->feature_types, 0);

    if (state == 0)
      str_cpy_charstr(name, "nonconserved");
    else if (state == source_mod->tree->nnodes)
      str_cpy_charstr(name, "conserved");
    else {
      n = lst_get_ptr(source_mod->tree->nodes, bdphmm->state_to_branch[state]);
      str_cpy_charstr(name, n->name);
      if (state < source_mod->tree->nnodes)
        str_append_charstr(name, "-death");
      else 
        str_append_charstr(name, "-birth");      
    }
  }

  bdphmm->phmm = phmm_new(hmm, models, cm, NULL, MISSING_DATA);
  bd_set_transitions(bdphmm);

  /* set up indel model, if necessary */
  if (alpha_c > 0) {
    bdphmm->indel_mods = smalloc(nstates * sizeof(void*));
    for (state = 0; state < nstates; state++) {
      List *l;

      n = lst_get_ptr(models[state]->tree->nodes, 
                      bdphmm->state_to_branch[state]);
      tr_partition_nodes(models[state]->tree, n, inside, outside);

      /* initialize alpha, beta, tau with nonconserved params */
      for (i = 0; i < models[state]->tree->nnodes; i++) {
        alpha[i] = alpha_n; beta[i] = beta_n; tau[i] = tau_n;
      }

      /* now change to conserved params for conserved branches */
      if (state < models[state]->tree->nnodes)  
        l = outside;            /* death */
      else 
        l = inside;             /* birth; includes conserved */

      for (i = 0; i < lst_size(l); i++) {
        n = lst_get_ptr(l, i);
        alpha[n->id] = alpha_c;
        beta[n->id] = beta_c;
        tau[n->id] = tau_c;
      }

      bdphmm->indel_mods[state] = im_new(alpha, beta, tau, models[state]->tree);
    }
  }
  else 
    bdphmm->indel_mods = NULL;

  lst_free(inside);
  lst_free(outside);
  sfree(alpha);
  sfree(beta);
  sfree(tau);
  return bdphmm;
}

void bd_set_transitions(BDPhyloHmm *bdphmm) {
  int state;
  HMM *hmm = bdphmm->phmm->hmm;
  TreeNode *tree = bdphmm->phmm->mods[0]->tree;

  /* set transition probs based on mu, nu, and phi */
  for (state = 1; state < hmm->nstates; state++) {
    if (state == tree->nnodes)
      mm_set(hmm->transition_matrix, 0, state, bdphmm->nu * (1-bdphmm->phi));
    else 
      mm_set(hmm->transition_matrix, 0, state, bdphmm->nu * bdphmm->phi / 
             (hmm->nstates-2));
    mm_set(hmm->transition_matrix, state, 0, bdphmm->mu);
    mm_set(hmm->transition_matrix, state, state, 1-bdphmm->mu);
  }
  mm_set(hmm->transition_matrix, 0, 0, 1-bdphmm->nu);

  /* begin transitions */
  vec_set(hmm->begin_transitions, 0, 
          bdphmm->mu/(bdphmm->mu + bdphmm->nu));
  for (state = 1; state < hmm->nstates; state++) {
    if (state == tree->nnodes)
      vec_set(hmm->begin_transitions, state, bdphmm->nu * (1-bdphmm->phi)/
              (bdphmm->mu + bdphmm->nu));
    else 
      vec_set(hmm->begin_transitions, state, bdphmm->nu * bdphmm->phi /
              ((bdphmm->mu + bdphmm->nu) * (hmm->nstates-2)));
  }

  hmm_reset(hmm);
}

/* adjust emission probabilities to prevent birth/death events from
   spanning regions where they would be supported only by missing data */
void bd_handle_missing_data(BDPhyloHmm *bdphmm, MSA *msa) {
  List **uninform = smalloc(msa->ss->ntuples * sizeof(void*));
  TreeNode *n, *tree = bdphmm->phmm->mods[0]->tree;
  typedef enum {MISSING, DATA, DATA_LEFT, DATA_RIGHT} node_type;
  node_type *mark = smalloc(tree->nnodes * sizeof(node_type));
  List *traversal;
  int i, j;

  /* enumerate the states for which each tuple is "uninformative" */
  for (i = 0; i < msa->ss->ntuples; i++) {
    uninform[i] = NULL;
    for (j = 0; j < msa->nseqs; j++)
      if (msa->is_missing[(int)ss_get_char_tuple(msa, i, j, 0)])
        break;

    if (j == msa->nseqs) continue; /* no missing chars so no need to
                                      look further */

    /* identify branches separating subtrees of all missing data from
       the rest of the tree, using a simple recursive algorithm; the
       tuple is uninformative about birth/death events on these branches */
    uninform[i] = lst_new_int(5);
    traversal = tr_postorder(tree);
    for (j = 0; j < lst_size(traversal); j++) {
      n = lst_get_ptr(traversal, j);
      if (n->lchild == NULL) {	/* leaf */
        if (msa->is_missing[(int)ss_get_char_tuple(msa, i, bdphmm->phmm->mods[0]->msa_seq_idx[n->id], 0)])
          mark[n->id] = MISSING;
        else
          mark[n->id] = DATA;
      }
      else if (n == tree) {
        if (mark[n->lchild->id] == MISSING || mark[n->rchild->id] == MISSING)
          mark[n->id] = MISSING;
        else 
          mark[n->id] = DATA;
      }
      else {			/* internal node */
        if (mark[n->lchild->id] == MISSING && mark[n->rchild->id] == MISSING)
          mark[n->id] = MISSING;
        else if (mark[n->lchild->id] != MISSING && mark[n->rchild->id] != MISSING)
          mark[n->id] = DATA;
        else if (mark[n->lchild->id] != MISSING)
          mark[n->id] = DATA_LEFT;
        else
          mark[n->id] = DATA_RIGHT;
      }
    }
    traversal = tr_preorder(tree);
    for (j = 0; j < lst_size(traversal); j++) {
      n = lst_get_ptr(traversal, j);
      if (n == tree) continue;
      
      /* resolve ambiguity, if necessary */
      if (mark[n->id] == DATA_LEFT || mark[n->id] == DATA_RIGHT)
        mark[n->id] = mark[n->parent->id];

      /* now check if missing data in subtree; if so, the tuple is 
         uninformative wrt this branch */
      if (mark[n->id] == MISSING) {
        lst_push_int(uninform[i], bdphmm->branch_to_state_death[n->id]);
        if (bdphmm->branch_to_state_birth[n->id] != -1)
          lst_push_int(uninform[i], bdphmm->branch_to_state_birth[n->id]);
      }
    }
  }

  /* now zero out all associated emissions probs */
  for (i = 0; i < msa->length; i++) {
    if (uninform[msa->ss->tuple_idx[i]] == NULL) continue;
    for (j = 0; j < lst_size(uninform[msa->ss->tuple_idx[i]]); j++) {
      int mod = lst_get_int(uninform[msa->ss->tuple_idx[i]], j);
      bdphmm->phmm->emissions[mod][i] = NEGINFTY;
    }
  }

  for (i = 0; i < msa->ss->ntuples; i++) 
    if (uninform[i] != NULL) lst_free(uninform[i]);
  sfree(uninform);
  sfree(mark);  
}

/* compute log odds scores for predictions (wrt neutral), using
   previously computed emissions */
void bd_score_predictions(BDPhyloHmm *bdphmm, GFF_Set *predictions) {
  int i, j, state;
  GFF_Feature *f;

  for (i = 0; i < lst_size(predictions->features); i++) {
    f = lst_get_ptr(predictions->features, i);
    f->score = 0;
    f->score_is_null = FALSE;
    state = cm_get_category(bdphmm->phmm->cm, f->feature);
    for (j = f->start-1; j < f->end; j++) 
      f->score += (bdphmm->phmm->emissions[state][j] - 
                   bdphmm->phmm->emissions[0][j]);
  }
}

/* combine indel emissions with substitution-based emissions */
void bd_add_indel_emissions(BDPhyloHmm *bdphmm, IndelHistory *ih) {
  int state, i;
  double *col_logl = smalloc(ih->ncols * sizeof(double));
  for (state = 0; state < bdphmm->phmm->hmm->nstates; state++) {   
    im_column_logl(ih, bdphmm->indel_mods[state], col_logl);
    for (i = 0; i < ih->ncols; i++) 
      bdphmm->phmm->emissions[state][i] += col_logl[i];
  }
  sfree(col_logl);
}

/* these two functions for use by bd_estimate_transitions */
void unpack_params(Vector *params, BDPhyloHmm *bdphmm) {
  int params_idx = 0;
  if (bdphmm->estim_omega)
    bdphmm->mu = 1/vec_get(params, params_idx++);
  if (bdphmm->estim_gamma) {
    double gamma = vec_get(params, params_idx++);
    bdphmm->nu = gamma/(1-gamma) * bdphmm->mu;
  }
  if (bdphmm->estim_phi)
    bdphmm->phi = vec_get(params, params_idx++);
  bd_set_transitions(bdphmm);
}

double lnl_wrapper(Vector *params, void *data) {
  BDPhyloHmm *bdphmm = data;
  unpack_params(params, data);
  return -hmm_forward(bdphmm->phmm->hmm, bdphmm->phmm->emissions, 
                      bdphmm->phmm->alloc_len, bdphmm->phmm->forward);
}

/* estimate free parameters */
double bd_estimate_transitions(BDPhyloHmm *bdphmm, MSA *msa) {
  int i, nparams = 0;
  Vector *params, *lb, *ub;
  double retval;

  if (bdphmm->phmm->forward == NULL) {
    bdphmm->phmm->forward = smalloc(bdphmm->phmm->hmm->nstates * 
                                    sizeof(double*));
    for (i = 0; i < bdphmm->phmm->hmm->nstates; i++)
      bdphmm->phmm->forward[i] = smalloc(bdphmm->phmm->alloc_len * 
                                         sizeof(double));
  }

  if (bdphmm->estim_omega) nparams++;
  if (bdphmm->estim_gamma) nparams++;
  if (bdphmm->estim_phi) nparams++;

  params = vec_new(nparams);
  lb = vec_new(nparams);
  ub = vec_new(nparams);
  vec_set_all(lb, 1e-6);

  nparams = 0;
  if (bdphmm->estim_omega) {
    vec_set(params, nparams, 1/bdphmm->mu);
    vec_set(ub, nparams++, INFTY);
  }
  if (bdphmm->estim_gamma) {
    vec_set(params, nparams, bdphmm->nu / (bdphmm->mu + bdphmm->nu));
    vec_set(ub, nparams++, 0.5 /* 1-1e-6 */);
  }
  if (bdphmm->estim_phi) {
    vec_set(params, nparams, bdphmm->phi);
    vec_set(ub, nparams++, 0.7 /* 1-1e-6 */);
  }

  opt_bfgs(lnl_wrapper, params, bdphmm, &retval, lb, ub, stderr, NULL, 
           OPT_HIGH_PREC, NULL, NULL);

  unpack_params(params, bdphmm);

  vec_free(params);
  vec_free(lb);
  vec_free(ub);

  return retval;
}

