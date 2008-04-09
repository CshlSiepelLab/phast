
/* dmotif phylo-HMM */

#include <dmotif_phmm.h>
#include <lists.h>
#include <sufficient_stats.h>
#include <numerical_opt.h>
#include <tree_model.h>
#include <msa.h>

/* create a new birth-death phylo-HMM based on parameter values */
DMotifPhyloHmm *dm_new(TreeModel *source_mod, PSSM *m, double rho, double mu, 
                       double nu, double phi, double zeta, double alpha_c, 
                       double beta_c, double tau_c, double alpha_n, 
                       double beta_n, double tau_n, int estim_gamma, 
                       int estim_omega, int estim_phi, int estim_zeta) {
  int i, state, pos;
  int nleaves = (source_mod->tree->nnodes + 1) / 2; 
  int k = 2*nleaves - 3;                            /* no. branches */
  int nstates = (2 + 2*k) * (m->width+1); /* number of states */
  HMM *hmm;
  TreeModel **models = smalloc(nstates * sizeof(void*));
  TreeNode *n;
  CategoryMap *cm;
  List *inside = lst_new_ptr(source_mod->tree->nnodes),
    *outside = lst_new_ptr(source_mod->tree->nnodes);
  DMotifPhyloHmm *dm;
  double *alpha = smalloc(source_mod->tree->nnodes * sizeof(double)), 
    *beta = smalloc(source_mod->tree->nnodes * sizeof(double)), 
    *tau = smalloc(source_mod->tree->nnodes * sizeof(double));
  assert(rho > 0 && rho < 1);  
  assert(mu > 0 && 2*mu < 1);
  assert(nu > 0 && nu < 1);

  dm = smalloc(sizeof(DMotifPhyloHmm));
  dm->state_to_branch = smalloc(nstates * sizeof(int));
  dm->state_to_event = smalloc(nstates * sizeof(int));
  dm->state_to_motifpos = smalloc(nstates * sizeof(int));
  dm->branch_to_states = smalloc(source_mod->tree->nnodes * sizeof(void*));
  dm->rho = rho;
  dm->mu = mu;
  dm->nu = nu;
  dm->phi = phi;
  dm->zeta = zeta;
  dm->estim_gamma = estim_gamma;
  dm->estim_omega = estim_omega;
  dm->estim_phi = estim_phi;
  dm->estim_zeta = estim_zeta;
  dm->m = m;
  dm->k = k;

  /* set up mappings between states and branches.  Branches are
     indexed by ids of child nodes; "birth" means functional element
     arose on branch, "death" means element was lost on branch. */

  for (i = 0; i < source_mod->tree->nnodes; i++)
    dm->branch_to_states[i] = lst_new_int(m->width + 1);

  state = 0;
  for (pos = -1; pos < m->width; pos++) {
    for (i = 0; i < source_mod->tree->nnodes; i++) {
      /* note: i=0 implies root note and fully nonconserved element */
      n = lst_get_ptr(source_mod->tree->nodes, i);
      dm->state_to_branch[state] = n->id;
      lst_push_int(dm->branch_to_states[n->id], state);
      dm->state_to_motifpos[state] = pos;
      dm->state_to_event[state] = (i == 0 ? NEUT : DEATH);
      state++;
    }
    for (i = 0; i < source_mod->tree->nnodes; i++) {
      /* note: i=0 implies root note and fully conserved element */
      n = lst_get_ptr(source_mod->tree->nodes, i);

      if (n->parent == source_mod->tree) 
        continue; 
        /* branches beneath root ignored -- events can't be
           distinguished from death events on opposite branch */

      dm->state_to_branch[state] = n->id;
      lst_push_int(dm->branch_to_states[n->id], state);
      dm->state_to_motifpos[state] = pos;
      dm->state_to_event[state] = (i == 0 ? CONS : BIRTH);
      state++;
    }
  }
  assert(state == nstates);
  /* note: this is structured so that state 0 is the nonconserved
     model (like element "died" prior to root) and state nnodes is the
     fully conserved model (like element was "born" prior to root) */

  /* now create tree models and scale branches as needed */  
  for (state = 0; state < nstates; state++) { 
    pos = dm->state_to_motifpos[state];
    if (pos == -1 && dm->state_to_event[state] == NEUT)
      models[state] = source_mod;
    else if (pos == -1)
      models[state] = tm_create_copy(source_mod);
    else 
      models[state] = tm_new(tr_create_copy(source_mod->tree), NULL, 
                             m->probs[pos], F81, DEFAULT_ALPHABET, 
                             1, -1, NULL, -1);
    /* F81 model with motif freqs as equilibrium freqs */

    if (dm->state_to_event[state] == NEUT)
      continue;

    n = lst_get_ptr(models[state]->tree->nodes, dm->state_to_branch[state]);
    tr_partition_nodes(models[state]->tree, n, inside, outside);

    /* adjust branch lengths if BIRTH, DEATH, or CONS */
    if (dm->state_to_event[state] == DEATH) 
      for (i = 0; i < lst_size(outside); i++) 
        ((TreeNode*)lst_get_ptr(outside, i))->dparent *= rho;
    else if (dm->state_to_event[state] == BIRTH ||
             dm->state_to_event[state] == CONS)
      for (i = 0; i < lst_size(inside); i++) 
        ((TreeNode*)lst_get_ptr(inside, i))->dparent *= rho;

    /* if BIRTH or DEATH, use alternative subst model (the background
       model) for non-motif branches */
     if (dm->state_to_event[state] == DEATH)  
       dm_set_backgd_branches(models[state], source_mod, inside); 
     else if (dm->state_to_event[state] == BIRTH)  
       dm_set_backgd_branches(models[state], source_mod, outside); 
    
    /* if birth event, reroot tree.  This ensures that the correct
       stationary distribution is used at the root of the subtree in
       question */ 
     if (dm->state_to_event[state] == BIRTH) 
       tr_reroot(models[state]->tree, n, TRUE); 

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

    if (dm->state_to_event[state] == NEUT)
      str_cpy_charstr(name, "nonconserved");
    else if (dm->state_to_event[state] == CONS)
      str_cpy_charstr(name, "conserved");
    else {
      n = lst_get_ptr(source_mod->tree->nodes, dm->state_to_branch[state]);
      str_cpy_charstr(name, n->name);

      if (dm->state_to_event[state] == DEATH)
        str_append_charstr(name, "-death");
      else 
        str_append_charstr(name, "-birth");      
    }
    
    if (dm->state_to_motifpos[pos] == -1)
        str_append_charstr(name, "-background");      
    else {
        str_append_charstr(name, "-M");      
        str_append_int(name, dm->state_to_motifpos[pos]+1);      
    }
  }

  dm->phmm = phmm_new(hmm, models, cm, NULL, MISSING_DATA);
  dm_set_transitions(dm);

  /* set up indel model, if necessary */
  if (alpha_c > 0) {
    dm->indel_mods = smalloc(nstates * sizeof(void*));
    for (state = 0; state < nstates; state++) {
      List *l;

      n = lst_get_ptr(models[state]->tree->nodes, 
                      dm->state_to_branch[state]);
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

      dm->indel_mods[state] = im_new(alpha, beta, tau, models[state]->tree);
    }
  }
  else 
    dm->indel_mods = NULL;

  lst_free(inside);
  lst_free(outside);
  free(alpha);
  free(beta);
  free(tau);
  return dm;
}

void dm_set_transitions(DMotifPhyloHmm *dm) {
  int s1, s2;
  HMM *hmm = dm->phmm->hmm;

  /* set transition probs based on mu, nu, and phi */
  for (s1 = 0; s1 < hmm->nstates; s1++) {
    int p1 = dm->state_to_motifpos[s1],
      b1 = dm->state_to_branch[s1],
      e1 = dm->state_to_event[s1];
    double rowsum = 0;

    for (s2 = 0; s2 < hmm->nstates; s2++) {
      int p2 = dm->state_to_motifpos[s2],
        b2 = dm->state_to_branch[s2],
        e2 = dm->state_to_event[s2];
      double trans = 0;

      /* first set motif trans probs */
      if (p1 == -1 && p2 == 0)
        trans = dm->zeta;
      else if (p1 == -1 && p2 == -1)
        trans = (1 - dm->zeta);
      else if (p2 == p1 + 1 || (p2 == -1 && p1 == dm->m->width-1))
        trans = 1;
      else 
        continue;              

      /* now mult by DLESS trans probs */
      if (e1 == NEUT) {
        if (e2 == NEUT)
          trans *= (1 - dm->nu);
        else {                  /* e2 is CONS, BIRTH, or DEATH */
          if (e2 == CONS)
            trans *= dm->nu * (1 - dm->phi);
          else
            trans *= dm->nu * dm->phi / (2*dm->k);
        }
      }        
      else {                    /* e1 is CONS, BIRTH, or DEATH */
        if (e2 == NEUT)
          trans *= dm->mu;
        else if (e2 == e1 && b2 == b1)
          trans *= (1 - dm->mu);
        else continue;          /* has to be zero */
      }

      mm_set(hmm->transition_matrix, s1, s2, trans);      
    }

    /* check that sum of outgoing arcs is 1 */
    for (s2 = 0; s2 < hmm->nstates; s2++) 
      rowsum += mm_get(hmm->transition_matrix, s1, s2);      
    assert(fabs(rowsum-1) < 1e-6);
  }


  /* begin transitions */
  for (s1 = 0; s1 < hmm->nstates; s1++) {
    /* motif stationary prob */
    double prob;
    if (dm->state_to_motifpos[s1] == -1)
      prob = 1/dm->zeta;
    else
      prob = dm->zeta / dm->m->width;

    /* multiply by dless stationary prob */
    if (dm->state_to_event[s1] == NEUT)
      prob *= dm->mu/(dm->mu + dm->nu);
    else if (dm->state_to_event[s1] == CONS)
      prob *= dm->nu * (1-dm->phi) / (dm->mu + dm->nu);
    else 
      prob *= dm->nu * dm->phi / ((dm->mu + dm->nu) * dm->k);

    vec_set(hmm->begin_transitions, s1, prob);
  }

  hmm_reset(hmm);
}

/* adjust emission probabilities to prevent birth/death events from
   spanning regions where they would be supported only by missing data */
void dm_handle_missing_data(DMotifPhyloHmm *dm, MSA *msa) {
  List **uninform = smalloc(msa->ss->ntuples * sizeof(void*));
  TreeNode *n, *tree = dm->phmm->mods[0]->tree;
  typedef enum {MISSING, DATA, DATA_LEFT, DATA_RIGHT} node_type;
  node_type *mark = smalloc(tree->nnodes * sizeof(node_type));
  List *traversal;
  int i, j, l;

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
        if (msa->is_missing[(int)ss_get_char_tuple(msa, i, dm->phmm->mods[0]->msa_seq_idx[n->id], 0)])
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
        else if (mark[n->lchild->id] != MISSING && mark[n->rchild->id != MISSING])
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
        for (l = 0; l < lst_size(dm->branch_to_states[n->id]); l++)
          lst_push_int(uninform[i], lst_get_int(dm->branch_to_states[n->id], l));
      }
    }
  }

  /* now zero out all associated emissions probs */
  for (i = 0; i < msa->length; i++) {
    if (uninform[msa->ss->tuple_idx[i]] == NULL) continue;
    for (j = 0; j < lst_size(uninform[msa->ss->tuple_idx[i]]); j++) {
      int mod = lst_get_int(uninform[msa->ss->tuple_idx[i]], j);
      dm->phmm->emissions[mod][i] = NEGINFTY;
    }
  }

  for (i = 0; i < msa->ss->ntuples; i++) 
    if (uninform[i] != NULL) lst_free(uninform[i]);
  free(uninform);
  free(mark);  
}

/* compute log odds scores for predictions (wrt neutral), using
   previously computed emissions */
void dm_score_predictions(DMotifPhyloHmm *dm, GFF_Set *predictions) {
  int i, j, state, len;
  GFF_Feature *f;

  for (i = 0; i < lst_size(predictions->features); i++) {
    f = lst_get_ptr(predictions->features, i);
    f->score = 0;
    f->score_is_null = FALSE;
    len = f->end - f->start + 1;
    state = cm_get_category(dm->phmm->cm, f->feature);
    for (j = f->start-1; j < f->end; j++) 
      f->score += (dm->phmm->emissions[state][j] - 
                   dm->phmm->emissions[0][j]);
  }
}

/* combine indel emissions with substitution-based emissions */
void dm_add_indel_emissions(DMotifPhyloHmm *dm, IndelHistory *ih) {
  int state, i;
  double *col_logl = smalloc(ih->ncols * sizeof(double));
  for (state = 0; state < dm->phmm->hmm->nstates; state++) {   
    im_column_logl(ih, dm->indel_mods[state], col_logl);
    for (i = 0; i < ih->ncols; i++) 
      dm->phmm->emissions[state][i] += col_logl[i];
  }
  free(col_logl);
}

/* these two functions for use by dm_estimate_transitions */
void unpack_params(Vector *params, DMotifPhyloHmm *dm) {
  int params_idx = 0;
  if (dm->estim_omega)
    dm->mu = 1/vec_get(params, params_idx++);
  if (dm->estim_gamma) {
    double gamma = vec_get(params, params_idx++);
    dm->nu = gamma/(1-gamma) * dm->mu;
  }
  if (dm->estim_phi)
    dm->phi = vec_get(params, params_idx++);
  if (dm->estim_zeta)
    dm->zeta = vec_get(params, params_idx++);
  dm_set_transitions(dm);
}

double lnl_wrapper(Vector *params, void *data) {
  DMotifPhyloHmm *dm = data;
  unpack_params(params, data);
  return -hmm_forward(dm->phmm->hmm, dm->phmm->emissions, 
                      dm->phmm->alloc_len, dm->phmm->forward);
}

/* estimate free parameters */
double dm_estimate_transitions(DMotifPhyloHmm *dm, MSA *msa) {
  int i, nparams = 0;
  Vector *params, *lb, *ub;
  double retval;

  if (dm->phmm->forward == NULL) {
    dm->phmm->forward = smalloc(dm->phmm->hmm->nstates * 
                                    sizeof(double*));
    for (i = 0; i < dm->phmm->hmm->nstates; i++)
      dm->phmm->forward[i] = smalloc(dm->phmm->alloc_len * 
                                         sizeof(double));
  }

  if (dm->estim_omega) nparams++;
  if (dm->estim_gamma) nparams++;
  if (dm->estim_phi) nparams++;
  if (dm->estim_zeta) nparams++;

  params = vec_new(nparams);
  lb = vec_new(nparams);
  ub = vec_new(nparams);
  vec_set_all(lb, 1e-6);

  nparams = 0;
  if (dm->estim_omega) {
    vec_set(params, nparams, 1/dm->mu);
    vec_set(ub, nparams++, INFTY);
  }
  if (dm->estim_gamma) {
    vec_set(params, nparams, dm->nu / (dm->mu + dm->nu));
    vec_set(ub, nparams++, 0.5 /* 1-1e-6 */);
  }
  if (dm->estim_phi) {
    vec_set(params, nparams, dm->phi);
    vec_set(ub, nparams++, 0.7 /* 1-1e-6 */);
  }
  if (dm->estim_zeta) {
    vec_set(params, nparams, dm->zeta);
    vec_set(ub, nparams++, 0.1);
  }

  opt_bfgs(lnl_wrapper, params, dm, &retval, lb, ub, stderr, NULL, 
           OPT_HIGH_PREC, NULL);

  unpack_params(params, dm);

  vec_free(params);
  vec_free(lb);
  vec_free(ub);

  return retval;
}

void dm_set_backgd_branches(TreeModel *tm, TreeModel *backgd_mod, 
                            List *nodelist) {
  int i;
  tm->alt_subst_mods = smalloc(tm->tree->nnodes * sizeof(void*));
  for (i = 0; i < tm->tree->nnodes; i++) tm->alt_subst_mods[i] = NULL;

  for (i = 0; i < lst_size(nodelist); i++) {
    TreeNode *n = lst_get_ptr(nodelist, i);
    tm->alt_subst_mods[n->id] = 
      tm_new_alt_subst_mod(backgd_mod->subst_mod,
                           vec_create_copy(backgd_mod->backgd_freqs),
                           mm_create_copy(backgd_mod->rate_matrix));
  }
}
