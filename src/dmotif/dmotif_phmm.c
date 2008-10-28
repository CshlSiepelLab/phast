 
/* dmotif phylo-HMM */
#include <dmotif_phmm.h>
#include <lists.h>
#include <sufficient_stats.h>
#include <numerical_opt.h>
#include <tree_model.h>
#include <msa.h>
#include <multi_msa.h>
#include "category_map.h"
#include <dmotif_indel_mod.h>

/* create a new birth-death phylo-HMM based on parameter values */
DMotifPhyloHmm *dm_new(TreeModel *source_mod, PSSM *m, double rho, double mu, 
                       double nu, double phi, double zeta, double alpha_c, 
                       double beta_c, double tau_c, double epsilon_c, 
		       double alpha_n, double beta_n, double tau_n,
		       double epsilon_n, int estim_gamma, int estim_omega, 
		       int estim_phi, int estim_zeta) {
  int i, state, pos, e, b;
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
    *tau = smalloc(source_mod->tree->nnodes * sizeof(double)),
    *epsilon = smalloc(source_mod->tree->nnodes * sizeof(double));
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
  for (i = 0; i < source_mod->tree->nnodes; i++){
    for (pos = -1; pos < m->width; pos++) {
      /* note: i=0 implies root note and fully nonconserved element */
      n = lst_get_ptr(source_mod->tree->nodes, i);
      dm->state_to_branch[state] = n->id;
      lst_push_int(dm->branch_to_states[n->id], state);
      dm->state_to_motifpos[state] = pos;
      dm->state_to_event[state] = (i == 0 ? NEUT : DEATH);
      state++;
    }
    for (pos = -1; pos < m->width; pos++) {
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
    else {
      models[state] = tm_new(tr_create_copy(source_mod->tree), NULL, 
                             m->probs[pos], F81, DEFAULT_ALPHABET, 
                             1, -1, NULL, -1);
      //      fprintf(stderr, "state %d pos %d ", state, pos);
      //      vec_print(m->probs[pos], stderr);
    }
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
  cm = dm_create_catmap(dm, source_mod, NULL);

  dm->phmm = phmm_new(hmm, models, cm, NULL, MISSING_DATA);

  dm_set_transitions(dm);

  /* set up indel model, if necessary */
  if (alpha_c > 0) {
    dm->indel_mods = smalloc(nstates * sizeof(void*));
    for (state = 0; state < nstates; state++) {
      List *l;
      /* Need event, motif position and branch to apply appropriate model 
	 params for motif and non-motif states */
      e = dm->state_to_event[state];
      pos = dm->state_to_motifpos[state];
      b = dm->state_to_branch[state];

      n = lst_get_ptr(models[state]->tree->nodes, b);
      tr_partition_nodes(models[state]->tree, n, inside, outside);

      /* initialize alpha, beta, tau with nonconserved params */
      for (i = 0; i < models[state]->tree->nnodes; i++) {
        alpha[i] = alpha_n; beta[i] = beta_n; tau[i] = tau_n;
	epsilon[i] = epsilon_n;
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
	epsilon[n->id] = epsilon_c;
      }

      dm->indel_mods[state] = dmih_new(alpha, beta, tau, epsilon,
				       models[state]->tree, e, l,
				       (pos == -1 ? 0 : 1));
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
  //  TreeNode *nn;

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

      /* Explicitly set params for each type of allowed transition */
      if ( e1 == NEUT ) { /* Transitions from NEUT states */
	if ( e2 == NEUT ) {
	  if (p1 == -1 && p2 == -1) /* Nb self-transition */
	    trans = ((1 - dm->zeta) * (1 - dm->nu));
	  else if (p1 == -1 && p2 == 0) /* first motif position */
	    trans = (dm->zeta * (1 - dm->nu));
	  else if (p2 == (p1 + 1) || (p1 == dm->m->width-1 && p2 == -1))
	    /* transitions between motif states and from NMw to Nb*/
	    trans = 1;
	  else continue;
	} else if ( e2 == CONS ) {
	  if (p1 == -1 && p2 == -1) /* Nb to Cb */
	    trans = (dm->nu * (1 - dm->zeta) * (1 - dm->phi));
	  else if (p1 == -1 && p2 == 0) /* first motif position */
	    trans = (dm->nu * dm->zeta * (1 - dm->phi));
	  else continue;
	} else { /* e2 == BIRTH or DEATH */
	  if (p1 == -1 && p2 == -1 ) /* lineage-specific background */
	    trans = (dm->nu * dm->phi * (1 - dm->zeta)) / (2 * dm->k);
	  else if (p1 == -1 && p2 == 0 ) /* lineage-specific M1 */
	    trans = (dm->nu * dm->phi * dm->zeta) / (2 * dm->k);
	  else continue;
	}
      } else { /* Transitions from CONS, BIRTH and DEATH states */
        if (e2 == e1 && b2 == b1) {
	  if (p1 == -1 && p2 == -1) /* background self-transitions */
	    trans = ((1 - dm->zeta) * (1 - dm->mu));
	  else if (p1 == -1 && p2 == 0) /* background to M1 */
	    trans = (dm->zeta * (1 - dm->mu));
	  else if (p1 == dm->m->width-1 && p2 == -1) /* Mw to B */
	    trans = (1 - dm->mu);
	  else if (p2 == (p1 + 1) && p2 < dm->m->width) /* transitions between 
							motif positions */
	    trans = 1;
	  else continue;
	} else if (e2 == NEUT && (p1 == -1 || p1 == dm->m->width-1) 
		   && p2 == -1) /* Cb or CMw to Nb */
	  trans = dm->mu;
	else continue;
      }
      mm_set(hmm->transition_matrix, s1, s2, trans);
      //      printf("s1 %d, s2 %d, p1 %d, p2 %d, trans %f\n", s1, s2, p1, p2, trans);
    }

    /* check that sum of outgoing arcs is 1 */
    for (s2 = 0; s2 < hmm->nstates; s2++)
      rowsum += mm_get(hmm->transition_matrix, s1, s2);
/*      printf("%d %f (%d, ", s1, rowsum, p1); */
/*      if (e1 == NEUT)  */
/*        printf("NEUT, ");  */
/*      else if (e1 == CONS)  */
/*        printf("CONS, ");  */
/*      else if (e1 == BIRTH)  */
/*        printf("BIRTH, ");  */
/*      else  */
/*        printf("DEATH, ");  */
  
/*      nn = lst_get_ptr(dm->phmm->mods[s1]->tree->nodes, dm->state_to_branch[s1]); */
/*      printf("%s)\n", nn->name); */

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
  int i, cbstate;
  String *cbname;
  GFF_Feature *f;

  cbname = str_new(STR_SHORT_LEN);
  str_append_charstr(cbname, "conserved-background");
  cbstate = cm_get_category(dm->phmm->cm, cbname);
  str_free(cbname);
  for (i = 0; i < lst_size(predictions->features); i++) {
    f = lst_get_ptr(predictions->features, i);
    str_append_charstr(f->attribute, "; ");
    dm_score_feat(dm, f, cbstate);
  }
}

/* Compute a log-odds score for a GFF feature */
void dm_score_feat(DMotifPhyloHmm *dm, GFF_Feature *f, int cbstate) {
  int j, bbstate, fstate;
  double s0, s1, s2;  

  s0 = s1 = s2 = 0; /* reset all scores, to be safe */  
  f->score = 0;
  f->score_is_null = FALSE;
  fstate = cm_get_category(dm->phmm->cm, f->feature);
  bbstate = (fstate - (dm->state_to_motifpos[fstate] + 1));
  for (j = f->start-1; j < f->end; j++) {
    s0 += (dm->phmm->emissions[fstate][j] -
	   dm->phmm->emissions[0][j]);
    s1 += (dm->phmm->emissions[fstate][j] -
	   dm->phmm->emissions[cbstate][j]);
    s2 += (dm->phmm->emissions[fstate][j] -
	   dm->phmm->emissions[bbstate][j]);
  }
  /* Set f->score to the most conservative score */
  f->score = (s0 < s1) ? ((s0 < s2) ? s0 : s2) : ((s1 < s2) ? s1 : s2);
  str_append_charstr(f->attribute, (s0 < s1) ?
		     ((s0 < s2) ? "BGDIST 0" : "BGDIST 1") :
		     ((s1 < s2) ? "BGDIST 1" : "BGDIST 2"));
}

/* Print a simple list of window scores for all motif windows in the msa. This
   can be used to produce a null score distribution for p-value computation */
void dm_print_motif_scores(DMotifPhyloHmm *dm) {
  int i, j, k, pstate, cbstate, bbstate;
  double s0, s1, s2;
  String *cbname;
  cbname = str_new(STR_SHORT_LEN);
  str_append_charstr(cbname, "conserved-background");
  cbstate = cm_get_category(dm->phmm->cm, cbname);


  for (i = 0; i <= ((dm->phmm->alloc_len - dm->m->width) + 1) ; i++) {
    for (j = 0; j < dm->phmm->hmm->nstates; j++) {
      if (dm->state_to_motifpos[j] != 0) { /* Only compute scores for motifs */
	if (dm->state_to_motifpos[j] == -1) {
	  bbstate = j;
	  s0 = s1 = s2 = 0;
	}
	continue;
      }
      pstate = j;
      for (k = i; k < i+dm->m->width; k++) {
	s0 += (dm->phmm->emissions[pstate][k] -
	       dm->phmm->emissions[0][k]);
	s1 += (dm->phmm->emissions[pstate][k] -
	       dm->phmm->emissions[cbstate][k]);
	s2 += (dm->phmm->emissions[pstate][k] -
	       dm->phmm->emissions[bbstate][k]);
	pstate++;
      }
      printf("%.3f\n", 
	     (s0 < s1) ? ((s0 < s2) ? s0 : s2) : ((s1 < s2) ? s1 : s2) );
    }
  }

}

/* combine indel emissions with substitution-based emissions */
void dm_add_indel_emissions(DMotifPhyloHmm *dm, IndelHistory *ih) {
  int state, i;
  double *col_logl = smalloc(ih->ncols * sizeof(double));
  for (state = 0; state < dm->phmm->hmm->nstates; state++) {   
    dmih_column_logl(ih, dm->indel_mods[state], col_logl);
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

/** Create a CategoryMap with features for each type of motif, with ranges
    defined by the width of the motif, and for each non-motif (background)
    state, with ranges size 1. */
CategoryMap* dm_create_catmap(DMotifPhyloHmm *dm, TreeModel *source_mod,
			      char *feature_prefix) {
  int i, j, k, nstates;
  TreeNode *node;
  nstates = (2 + 2*dm->k) * (dm->m->width+1); /* number of states */
/*   ncats = (4 * dm->k + 4)-1; */
  CategoryMap *retval = cm_new(nstates-1);
  CategoryRange *cr;
  for (i = 0; i < nstates; i++) {
    int p = dm->state_to_motifpos[i],
      b = dm->state_to_branch[i],
      e = dm->state_to_event[i];
    String *type;
    node = lst_get_ptr(source_mod->tree->nodes, b);
    if ( p == -1 ) /* background state */
      j = i;
    else if ( p == 0 ) /* First motif position */
      j = i + dm->m->width - 1;
    else  /* positions within motifs -- lumped into category with first pos. */
      continue;
    assert(j <= nstates);
    type = str_new(STR_SHORT_LEN);
    if (feature_prefix != NULL) str_append_charstr(type, feature_prefix);
    if ( e == NEUT ) { 
      str_append_charstr(type, "nonconserved");
      if ( p == -1 ) { /* 0th category is background */
      	retval->ranges[0] = 
	  cm_new_category_range(str_new_charstr(BACKGD_CAT_NAME), 0, 0);
 	continue;
      }
      else 
      	str_append_charstr(type, "-motif"); 
    } 
    else if ( e == CONS ) {
      str_append_charstr(type, "conserved");
      if ( p == -1 )
       	str_append_charstr(type, "-background");
      else 
    	str_append_charstr(type, "-motif");
    } 
    else if ( e == DEATH ) {
      str_append_charstr(type, "death-");
      str_append_charstr(type, node->name); 
      if ( p == -1 )
    	str_append_charstr(type, "-background");
      else
    	str_append_charstr(type, "-motif");
    }
    else { /* e == BIRTH */
      str_append_charstr(type, "birth-");
      str_append_charstr(type, node->name);  
      if ( p == -1 )
    	str_append_charstr(type, "-background");
      else  
    	str_append_charstr(type, "-motif");
    } 
/*      printf("i %d j %d nstates %d e %d b %d p %d (", i, j, nstates, e, b, p); */
/*      printf("%s)\n", type->chars); */
/*      printf("%d\n", ((int)((i-1) / dm->m->width))); */

    cr = cm_new_category_range(type, i, j);
    for (k = i ; k <= j; k++) 
      retval->ranges[k] = cr;
  }
  return retval;
}

/** Run the Viterbi algorithm and return a set of predictions.
    Emissions must have already been computed (see
    phmm_compute_emissions) */
GFF_Set* dm_phmm_predict_viterbi(DMotifPhyloHmm *dm, 
			      /**< DmotifPhyloHmm object */
                              char *seqname,
			      /**< seqname for feature set (e.g.,
				 "chr1") */
                              char *grouptag,
			      /**< tag to use for groups (e.g.,
                                   "exon_id", "transcript_id"); if
                                   NULL, a default will be used */
                              char *idpref,
			      /**< prefix for assigned ids (e.g.,
				 "chr1.15") (may be NULL) */
                              List *frame
			      /**< names of features for which to obtain
				 frame (NULL to ignore) */
			      ) {
  PhyloHmm *phmm = dm->phmm;
  int *path = (int*)smalloc(phmm->alloc_len * sizeof(int));
  GFF_Set *retval;

  if (phmm->emissions == NULL)
    die("ERROR: emissions required for phmm_viterbi_features.\n");
          
  hmm_viterbi(phmm->hmm, phmm->emissions, phmm->alloc_len, path);

  retval = dm_labeling_as_gff(phmm->cm, path, phmm->alloc_len, dm->m->width,
                              phmm->state_to_cat,
			      dm->state_to_motifpos,
                              phmm->reverse_compl, seqname, "DMOTIF",  
                              frame, grouptag, idpref);
  free(path);

  return retval;
}

/** Create a GFF_Set from a sequence of category/state numbers, using
   a specified category map and mapping from raw state numbers to
   category numbers.  */
GFF_Set *dm_labeling_as_gff(CategoryMap *cm, 
			    /**< CategoryMap to use in mapping */
                            int *path, 
			    /**< Raw sequence of state/category numbers  */
                            int length, 
			    /**< Length of sequence */
			    int w,
			    /**< Length of motif */
                            int *path_to_cat, 
			    /**< Mapping from raw numbers to
			       category numbers */
			    int *state_to_motifpos,
			    /**< Mapping of states to motif positions*/
                            int *reverse_compl, 
			    /**< Array of boolean values indicating
                                   whether each raw-sequence value
                                   corresponds to the reverse strand  */
                            char *seqname, 
			    /**< char string to use as 'seqname' in
			       generated GFF_Set  */
                            char *source, 
			    /**< char string to use as 'source' in
			       generated GFF_Set  */
                            List *frame_cats, 
			    /**< Categories for which to obtain frame
			       information (by name) */
                            char *grouptag,
                            /**< Tag to use to define groups in
                               GFF_Set (e.g., "transcript_id") */
                            char *idpref
			    /**< Prefix for ids of predicted
                                   elements (may be NULL).  Can be
                                   used to ensure ids are unique. */
                            ) {
  int beg, end, i, cat, lastcat, frame, groupno, mpos;
  GFF_Set *gff = gff_new_set_init("PHAST", PHAST_VERSION);
  int do_frame[cm->ncats+1];
  char strand;
  char groupstr[STR_SHORT_LEN];
  int ignore_0 = str_equals_charstr(cm_get_feature(cm, 0), BACKGD_CAT_NAME);
  /* ignore category 0 if background  */

  if (length <= 0) return gff;

  for (i = 0; i <= cm->ncats; i++) do_frame[i] = 0;
  if (frame_cats != NULL)
    for (i = 0; i < lst_size(frame_cats); i++) {
      cat = cm_get_category(cm, lst_get_ptr(frame_cats, i));
      if (cat != 0)             /* ignore background or unrecognized name */
        do_frame[cat] = 1;
    }

  groupno = 1;
  if (idpref != NULL)
    sprintf(groupstr, "%s \"%s.%d\"", grouptag != NULL ? grouptag : "id", 
	    idpref, groupno);
  else
    sprintf(groupstr, "%s \"%d\"", grouptag != NULL ? grouptag : "id", groupno);
  
  i = 0;
  cat = -1;
  while (i < length) {
    lastcat = cat;
    cat = cm->ranges[path_to_cat[path[i]]]->start_cat_no;
/*     printf("%d %d %d\n", path[i], path_to_cat[33], cat); */
    mpos = state_to_motifpos[path[i]];
    strand = reverse_compl[path[i]] ? '-' : '+';
    frame = do_frame[cat] ? path_to_cat[path[i]] - cat : GFF_NULL_FRAME;
    
    /* scan ahead until enter new category range (or reach end of seq) */
    beg = i + 1;                /* begin of feature (GFF coords) */
    for (i++; i < length && 
	   cm->ranges[path_to_cat[path[i]]]->start_cat_no == cat; i++);
    end = i;                    /* end of feature (GFF coords) */
    //    printf("mpos %d, start %d, end %d\n", mpos, beg, end);


    /* This is the wrong way of handling this -- ignores any motifs that have
     deletions in the reference species, even if they are strongly supported
     in other species. The indel model will do a better job of deciding whether
     these represent real motifs. */
/*     /\* If this is a motif feature and the start and end coordinates add up to */
/*        less than the motif length (i.e., for a partial motif lying at the end */
/*        of an alignment block or bordering a gap), skip this feature -- this  */
/*        will probably have to be modified once the indel model is in place. *\/ */
/*     if ( mpos >= 0 && (end-beg) < (w-1) ) */
/*       continue; */
    
    /* if minus strand, adjust frame to reflect end */
    if (strand == '-' && do_frame[cat]) 
      frame = path_to_cat[path[i-1]] - cat;
    
    /* if legitimate feature (non-background), then incorp into GFF_Set */
    if ((cat != 0 || !ignore_0) && mpos != -1)  /* create new feature and add */
      lst_push_ptr(gff->features, 
		   gff_new_feature(str_new_charstr(seqname), 
				   str_new_charstr(source), 
				   str_dup(cm_get_feature(cm, cat)), 
				   beg, end, 0, strand, frame, 
				   str_new_charstr(groupstr), TRUE));
    
/*     fprintf(stderr, "cat %d, beg %d, lastcat %d\n", cat, beg, lastcat); */
    if ((cat == 0 && beg > 1) || cat != lastcat) {
      groupno++;                /* increment group number each time a
				   new category or a 0 is encountered  */
      if (idpref != NULL)
	sprintf(groupstr, "%s \"%s.%d\"", grouptag != NULL ? grouptag : "id", 
		idpref, groupno);
      else
	sprintf(groupstr, "%s \"%d\"", grouptag != NULL ? grouptag : "id", 
		groupno);
    }
  }
  
  return gff;
}

void dms_sample_paths(DMotifPhyloHmm *dm, Multi_MSA *blocks,
		      double ***emissions, int bsamples, int nsamples,
		      int sample_interval, Hashtable *path_counts, 
		      int **priors, FILE *log, GFF_Set *reference, 
		      int ref_as_prior, int force_priors) {
  int i, j, k, l, **path, **trans, s1, s2, p1, p2, e1, e2, *counts,
    **ref_paths, *tmp_path, *stats, nwins;
  double ***forward_scores, llh, fpr, spc, ppv;
  char hsh_key[STR_SHORT_LEN], param[STR_SHORT_LEN];
  MSA *msa;
  GFF_Set *query_gff, *tmp_gff;
  GFF_Feature *f;

  /* Allocate space for the forward scores and paths */
  forward_scores = smalloc(blocks->nblocks * sizeof(double*));
  path = smalloc(blocks->nblocks * sizeof(int*));
  for (i = 0; i < blocks->nblocks; i++) {
    forward_scores[i] = smalloc(dm->phmm->hmm->nstates * sizeof(double*));
    msa = blocks->blocks[i];
    path[i] = smalloc(msa->length * sizeof(int));
    for (j = 0; j < dm->phmm->hmm->nstates; j++)
      forward_scores[i][j] = smalloc(msa->length * sizeof(double));
  }

  /* Allocate space for the transition counts and set to prior values */
  trans = smalloc(4 * sizeof(int*));
  for (i = 0; i < 4; i++) {
    trans[i] = smalloc(2 * sizeof(int));
    trans[i][0] = priors[i][0];
    trans[i][1] = priors[i][1];
  }

  /* If comparing to a reference GFF for logging purposes, allocate the
     stats array and compute the number of motif windows in the dataset. */
  if (reference != NULL && log != NULL) {
    stats = smalloc(4 * sizeof(int*));
    for (i = 0; i < blocks->nblocks; i++)
      nwins += (blocks->blocks[i]->length - (dm->m->width -1));      
  }
  /* If using the reference set as prior information, reconstruct a set of
     state paths from the GFF_Set. */
  if (reference != NULL && ref_as_prior == TRUE) {
    ref_paths = dm_gff_to_paths(dm, reference, blocks, -1);
  }

  for (i = 0; i < (bsamples + nsamples); i++) {
    fprintf(stderr, "*");

    /* Draw a set of params */
    dm->mu = beta_draw(trans[0][0], trans[0][1]);
    dm->nu = beta_draw(trans[1][0], trans[1][1]);
    dm->phi = beta_draw(trans[2][0], trans[2][1]);
    dm->zeta = beta_draw(trans[3][0], trans[3][1]);
    dm_set_transitions(dm);

    llh = 0;
    for (j = 0; j < blocks->nblocks; j++) {
      dm->phmm->emissions = emissions[j];
      msa = blocks->blocks[j]; 
      llh += hmm_forward(dm->phmm->hmm, emissions[j], msa->length,
			 forward_scores[j]);
      hmm_stochastic_traceback(dm->phmm->hmm, forward_scores[j], msa->length,
			       path[j]);
    }
 
    if (i % sample_interval) {
      continue;
    } else {
      /* If using reference set as priors, build a composite path containing
	 all known and predicted features and use this as the path */
      if (ref_as_prior == TRUE) {
	for (j = 0; j < blocks->nblocks; j++) {
	  tmp_path = dms_composite_path(dm, ref_paths[j], path[j], 
					msa->length, 0, force_priors);
	  for (k = 0; k < msa->length; k++)
	    path[j][k] = tmp_path[k];
	  free(tmp_path);
	}
      }
      
      /* Reset transition counts */
      for (j = 0; j < 4; j++) {
	trans[j][0] = priors[j][0];
	trans[j][1] = priors[j][1];
      }
      
      for (j = 0; j < blocks->nblocks; j++) {
	msa = blocks->blocks[j];
	s1 = -1;
	for (k = 0; k < msa->length; k++) {
	  s2 = path[j][k];
	  
/* 	  fprintf(stderr, "i %d, j %d, len %d, k %d, s2 %d\n", i, j,  */
/* 		  msa->length, k, s2); */
	  
	  p1 = (s1 == -1 ? -1 : dm->state_to_motifpos[s1]);
	  e1 = (s1 == -1 ? NEUT : dm->state_to_event[s1]);
	  p2 = dm->state_to_motifpos[s2];
	  e2 = dm->state_to_event[s2];
	  
	  /* Track different types of transitions for param estimation */
	  if (p1 == -1) { /* Transitions from bg states */
	    if (e1 == NEUT) { 
	      if (e2 != NEUT) { /* Affects nu */
		trans[1][0]++;
		if (e2 == BIRTH || e2 == DEATH) { /* Affects phi */
		  trans[2][0]++;
		} else {
		  trans[2][1]++;
		}
	      } else {
		trans[1][1]++;
	      }
	    } else { /* Source state is under selection */
	      if (e2 == NEUT) { /* Affects mu */
		trans[0][0]++;
	      } else {
		trans[0][1]++;
	      }
	    }
	    
	    if (p2 == 0) { /* Into a motif -- affects zeta */
	      trans[3][0]++;
	    } else { /* into another background category */
	      trans[3][1]++;
	    }
	    
	  } else if (p1 == (dm->m->width - 1)
		     && e1 != NEUT) { /* Transitions from motif end */
	    if (e2 == NEUT) { /* Affects mu */
	      trans[0][0]++;
	    } else {
	      trans[0][1]++;
	    }
	  }

	  /* Track entries into motifs for path sampling (after burn-in)*/
	  if (i >= bsamples) { 
	    if (p2 == 0) {
	      sprintf(hsh_key, "%d_%d", j, k);
	      counts = hsh_get(path_counts, hsh_key);
	      if (counts == -1) { /* First motif at this position */
		counts = smalloc(dm->phmm->hmm->nstates * sizeof(int));
		for (l = 0; l < dm->phmm->hmm->nstates; l++)
		  counts[l] = 0;
		counts[s2] = 1;
		hsh_put(path_counts, hsh_key, counts);
	      } else {
		counts[s2]++;
	      }
	    }
	    /* 	    fprintf(stderr, "key %s state %d count %i\n", hsh_key, s2, */
	    /* 		    ((int*)hsh_get(path_counts, hsh_key))[s2]); */
	  }
	  s1 = s2;
	}	
      }
    }
    
    /* Write stats to log file, if specified */
    if (log != NULL) {
      fprintf(log, "Sample = %d, total log likelihood = %f\n", i, llh);
      
      /* If using a reference set, compare known and predicted features and 
	 keep track of the matches, mismatches, etc. in the stats structure */
      if (reference != NULL) {
	fpr = spc = ppv = 0;
	query_gff = gff_new_set();	
	for (j = 0; j < 4; j++) /* Reset counters in stats array */
	  stats[j] = 0;
	for (j = 0; j < blocks->nblocks; j++) {
	  tmp_gff = gff_new_set();	
	  tmp_gff = dm_labeling_as_gff(dm->phmm->cm, path[j], msa->length,
				       dm->m->width,
				       dm->phmm->state_to_cat,
				       dm->state_to_motifpos,
				       dm->phmm->reverse_compl,
				       ((String*)lst_get(blocks->seqnames, j))->chars,
				       "DMSAMPLE", NULL, NULL, NULL);
	  for (k = 0; k < lst_size(tmp_gff->features); k++) {
	    f = lst_get_ptr(tmp_gff->features, k);
	    lst_push_ptr(query_gff->features, gff_new_feature_copy(f));
	  }
	  gff_free_set(tmp_gff);
	}
	dms_compare_gffs(reference, query_gff, stats, 0, NULL, NULL, NULL, 
			 NULL);
	fpr = (double)(stats[1] + stats[3])
	  / (double)(nwins - lst_size(reference->features));
	spc = (double)(nwins - lst_size(query_gff->features))
	  / (double)((nwins - lst_size(query_gff->features)) + stats[3] + stats[1]);
	if (stats[0])
	  ppv = (double)(stats[0]) / (double)(lst_size(query_gff->features));
	else
	  ppv = 0;

	fprintf(log, "True Positive Features: %d of %d (Sensitivity = %f)\n", 
		stats[0], lst_size(reference->features),
		((double)stats[0] / (double)lst_size(reference->features)));
	fprintf(log, "False Positive Features: \n\t%d total (FPR = %f)\n\t%d misidentified\n\t%d in locations with no annotated feature\n",
		(stats[3] + stats[1]), fpr, stats[1], stats[3]);
	fprintf(log, "False Negative Features: %d\n", stats[2]);
	fprintf(log, "Specificity: %f\n", spc);
	fprintf(log, "PPV: %f\n", ppv);
      }
      
      for (j = 0; j < 4; j++) {
	if (j == 0)
	  sprintf(param, "mu");
	else if (j == 1)
	  sprintf(param, "nu");
	else if (j == 2)
	  sprintf(param, "phi");
	else /* (j == 3) */
	  sprintf(param, "zeta");
	
	fprintf(log, "%s transitions: alpha = %d, beta = %d\n", param, 
		trans[j][0], trans[j][1]);
      }

      fprintf(log, "mu = %f, nu = %f, phi = %f, zeta = %f\n\n", dm->mu, 
	      dm->nu, dm->phi, dm->zeta);
      gff_free_set(query_gff);
    }
  }
  
  fprintf(stderr, "\nDone.\n");
  
  for (i = 0; i < blocks->nblocks; i++) {
    for (j = 0; j < dm->phmm->hmm->nstates; j++)
      free(forward_scores[i][j]);
    free(forward_scores[i]);
    free(path[i]);
  }
  free(forward_scores);
  free(path); 
}

void dms_read_priors(int **priors, FILE *prior_f) {
  Regex *mu_re = str_re_new("#[[:space:]]*MU[[:space:]]*([0-9]+[[:space:]]*[0-9]+)");
  Regex *nu_re = str_re_new("#[[:space:]]*NU[[:space:]]*([0-9]+[[:space:]]*[0-9]+)");
  Regex *phi_re = str_re_new("#[[:space:]]*PHI[[:space:]]*([0-9]+[[:space:]]*[0-9]+)");
  Regex *zeta_re = str_re_new("#[[:space:]]*ZETA[[:space:]]*([0-9]+[[:space:]]*[0-9]+)");

  String *line = str_new(STR_MED_LEN);
  List *matches = lst_new_ptr(2);
  List *counts = lst_new_ptr(2);
  int count;

  while (str_readline(line, prior_f) != EOF) {
    lst_clear(matches);
    lst_clear(counts);
    str_trim(line);
    if (line->length == 0) continue;
   
    if (str_re_match(line, mu_re, matches, 1) >= 0) {
      str_split(lst_get_ptr(matches, 1), NULL, counts);
      str_as_int(lst_get_ptr(counts, 0), &count);
      priors[0][0] = count;
      str_as_int(lst_get_ptr(counts, 1), &count);
      priors[0][1] = count;
    } else if (str_re_match(line, nu_re, matches, 1) >= 0) {
      str_split(lst_get_ptr(matches, 1), NULL, counts);
      str_as_int(lst_get_ptr(counts, 0), &count);
      priors[1][0] = count;
      str_as_int(lst_get_ptr(counts, 1), &count);
      priors[1][1] = count;
    } else if (str_re_match(line, phi_re, matches, 1) >= 0) {
      str_split(lst_get_ptr(matches, 1), NULL, counts);
      str_as_int(lst_get_ptr(counts, 0), &count);
      priors[2][0] = count;
      str_as_int(lst_get_ptr(counts, 1), &count);
      priors[2][1] = count;
    } else  if (str_re_match(line, zeta_re, matches, 1) >= 0) {
      str_split(lst_get_ptr(matches, 1), NULL, counts);
      str_as_int(lst_get_ptr(counts, 0), &count);
      priors[3][0] = count;
      str_as_int(lst_get_ptr(counts, 1), &count);
      priors[3][1] = count;
    } else { /* something unexpected */
      die("ERROR: Too many lines in priors file.\n");
    }
  }
  str_re_free(mu_re);
  str_re_free(nu_re);
  str_re_free(phi_re);
  str_re_free(zeta_re);
  str_free(line);
  lst_free(matches);
  lst_free(counts);
}

GFF_Feature* dms_motif_as_gff_feat(DMotifPhyloHmm *dm, double ***emissions,
				   Multi_MSA *blocks, char *key, int *counts, 
				   int nsamples, int sample_interval,
				   int cbstate) {
  int i, cat, seqnum, start, end, best,
    mtf_total, b_total, d_total, c_total, n_total;
  char *strand, delim[2], post[STR_LONG_LEN], *seqname;
  CategoryMap *cm = dm->phmm->cm;
  String *attributes, *key_string;
  List *split = lst_new(2, sizeof(String));
  GFF_Feature *f;

  key_string = str_new_charstr(key);
  attributes = str_new(STR_MED_LEN);

  sprintf(delim, "_");
  str_split(key_string, delim, split);
  str_free(key_string);
  str_as_int(lst_get_ptr(split, 0), &seqnum);
  str_as_int(lst_get_ptr(split, 1), &start);
  start++; /* Correct 0-based coordinate to 1-based */
  end = ((start + dm->m->width - 1) < blocks->blocks[seqnum]->length) ?
    (start + dm->m->width - 1) : (blocks->blocks[seqnum]->length - 1);
  seqname = ((char*)((String*)lst_get(blocks->seqnames, seqnum))->chars);

/*   fprintf(stderr, "key %s, seqnum %d, start %d, end %d, seqname %s\n", key, */
/* 	  seqnum, start, end, seqname); */
  
  sprintf(post, "ID \"%s.%d\"; ", seqname, start);
  str_append_charstr(attributes, post);

  /* First pass counts the total number of motifs predicted in this position
     and totals for conserved, gain, loss and neutral categories of motifs,
     along with finding the "best" category -- that is, the one with the
     highest posterior probability -- that will be used as the GFF feature
     type. */
  best = mtf_total = c_total = b_total = d_total = n_total = 0;
  for (i = 0; i < dm->phmm->hmm->nstates; i++) {
    if (counts[i] > 0) {
      mtf_total += counts[i];
      if (dm->state_to_event[i] == CONS)
	c_total += counts[i];
      else if (dm->state_to_event[i] == BIRTH)
	b_total += counts[i];
      else if (dm->state_to_event[i] == DEATH)
	d_total += counts[i];
      else /* dm->state_to_event[i] == NEUT */
	n_total += counts[i];
	
      if (counts[i] > counts[best]) {
	best = i;
      }
    }
  }

/*   fprintf(stderr, "nsamples %d, mtf_total %d, c_total %d, b_total %d, d_total %d, n_total %d\n", nsamples, mtf_total, c_total, b_total, d_total, n_total); */

  /* Compute marginal posteriors for total motifs and motifs of each flavor */
  sprintf(post, "MP %s %f; MP %s %f; MP %s %f; MP %s %f; MP %s %f; ",
	  "MOTIF", ((double)mtf_total / (double)(nsamples / sample_interval)),
	  "CONS", ((double)c_total / (double)mtf_total),
	  "BIRTH", ((double)b_total / (double)mtf_total),
	  "DEATH", ((double)d_total / (double)mtf_total),
	  "NEUT", ((double)n_total / (double)mtf_total));
  str_append_charstr(attributes, post);
  
  /* Second pass computes posterior probabilities of being in each flavor of 
     motif given that there is a motif in this position. */
  for (i = 0; i < dm->phmm->hmm->nstates; i++) {
    if (counts[i] > 0) {
      cat = cm->ranges[dm->phmm->state_to_cat[i]]->start_cat_no;
      strand = dm->phmm->reverse_compl[i] ? "-" : "+";
      sprintf(post, "PP %s %f; ", ((String*)cm_get_feature(cm, cat))->chars,
	      ((double)counts[i] / (double)mtf_total));
      str_append_charstr(attributes, post);
/*       printf("%s %d %s %s %d %d %s %s\n", key, counts[i], seqname,  */
/* 	     ((String*)cm_get_feature(cm, cat))->chars, start, end, strand,  */
/* 	     attributes->chars); */
    }
  }
  
/*   fprintf(stderr, "seqname: %s, attributes: %s\n", seqname, attributes->chars); */

  cat = cm->ranges[dm->phmm->state_to_cat[best]]->start_cat_no;
  strand = dm->phmm->reverse_compl[best] ? "-" : "+";  

  f = gff_new_feature(str_new_charstr(seqname),
		      str_new_charstr("DMSAMPLE"),
		      str_dup(cm_get_feature(cm, cat)),
		      start, end, 0, *strand, GFF_NULL_FRAME,
		      attributes, TRUE);
/*   fprintf(stderr, "%d %p %p \n", seqnum, emissions, */
/* 	  dm->phmm->emissions); */
  dm->phmm->emissions = emissions[seqnum];
  dm_score_feat(dm, f, cbstate);

/*   gff_print_feat(stderr, f); */

  lst_free(split);  

  return f;
}

/* Compare motif features in a query set to a target set and return stats on
   numbers of matches, mismatches, features unique to target list and features
   unique to query list. */
void dms_compare_gffs(GFF_Set *target_gff, GFF_Set *query_gff, int *stats,
		      int offset, GFF_Set *matches, GFF_Set *mismatches,
		      GFF_Set *unique_to_query, GFF_Set *unique_to_target) {
  /* offset can be used to adjust coordinates of the query to match the 
     coordinate base of the target set */
  int i, j, k, nmatches, *tmatch_idx, *qmatch_idx;
  GFF_Feature *ft, *fq;
  tmatch_idx = smalloc(lst_size(target_gff->features) * sizeof(int));
  qmatch_idx = smalloc(lst_size(query_gff->features) * sizeof(int));

  if (matches != NULL && (mismatches == NULL || unique_to_query == NULL ||
			  unique_to_target == NULL))
    die("ERROR in dms_compare_gffs: Not enough GFF_Sets for format!\n");

  nmatches = 0;
  for (i = 0; i < lst_size(query_gff->features); i++) {
    fq = lst_get_ptr(query_gff->features, i);
    for (j = 0; j < lst_size(target_gff->features); j++) {
      ft = lst_get_ptr(target_gff->features, j);
      if (str_compare(fq->seqname, ft->seqname) != 0) /* Different sequences */
	continue;
      if ((fq->start + offset) == ft->start) { 
	/* Motif found in same position */
	tmatch_idx[nmatches] = j;
	qmatch_idx[nmatches] = i;
	nmatches++;
	if (str_compare(fq->feature, ft->feature) == 0) { /* feature types
							     match */
	  stats[0]++;
	  if (matches != NULL)
	    lst_push_ptr(matches->features, gff_new_feature_copy(fq));
	} else { /* Correct position, wrong flavor of motif */
	  stats[1]++;
	  if (matches != NULL) {
	    lst_push_ptr(mismatches->features, gff_new_feature_copy(fq));
	    lst_push_ptr(mismatches->features, gff_new_feature_copy(ft));
	  }
	}
      }
    }
  }
  /* Number features unique to target */
  stats[2] = (lst_size(target_gff->features) - nmatches);
  /* Number features unique to query */
  stats[3] = (lst_size(query_gff->features) - nmatches);

  if (matches != NULL) { /* Build gff's for feats unique to target and query */
    /* unique_to_query set */
    j = lst_size(query_gff->features)-1;
    for (i = (nmatches-1); i >= -1; i--) {
      if (i == -1) {
	for (k = j; k >= 0; k--) {
	  fq = lst_get_ptr(query_gff->features, k);
	  lst_push_ptr(unique_to_query->features, gff_new_feature_copy(fq));
	}
      } else {
	for (k = j; k > qmatch_idx[i]; k--) {
	  fq = lst_get_ptr(query_gff->features, k);
	  lst_push_ptr(unique_to_query->features, gff_new_feature_copy(fq)); 
	}
      }
      j = (qmatch_idx[i] - 1);
    }
    /* unique_to_target set */
    j = lst_size(target_gff->features)-1;
    for (i = (nmatches-1); i >= -1; i--) {
      if (i == -1) {
	for (k = j; k >= 0; k--) {
	  fq = lst_get_ptr(target_gff->features, k);
	  lst_push_ptr(unique_to_target->features, gff_new_feature_copy(fq));
	}
      } else {
	for (k = j; k > tmatch_idx[i]; k--) {
	  fq = lst_get_ptr(target_gff->features, k);
	  lst_push_ptr(unique_to_target->features, gff_new_feature_copy(fq));
	}
      }
      j = (tmatch_idx[i] - 1);
    }
  }
}


/* Merge two GFF_Sets for motif features into a GFF containing all features
   from both sets. Sums posterior probabilities for positionally-matched
   features, but does not do any averaging or readjustment of scores. I.e.,
   the posteriors in the output may not be valid posteriors! Feature type is
   ignored when combining features. */
void dms_combine_gffs(GFF_Set *target_gff, GFF_Set *query_gff) {
    int i, j, k, nmatches, *tmatch_idx, *qmatch_idx;
  GFF_Feature *ft, *fq;
  tmatch_idx = smalloc(lst_size(target_gff->features) * sizeof(int));
  qmatch_idx = smalloc(lst_size(query_gff->features) * sizeof(int));

  nmatches = 0;
  for (i = 0; i < lst_size(query_gff->features); i++) {
    fq = lst_get_ptr(query_gff->features, i);
    for (j = 0; j < lst_size(target_gff->features); j++) {
      ft = lst_get_ptr(target_gff->features, j);
      if (str_compare(fq->seqname, ft->seqname) != 0) /* Different sequences */
	continue;
      if (fq->start == ft->start) { 
	/* Motif found in same position */
	tmatch_idx[nmatches] = j;
	qmatch_idx[nmatches] = i;
	nmatches++;

	/* Go through attributes fields, summing the posteriors for matching
	   features and simply appending unique types */

      }
    }
  }

  /* Now fold in all the features unique to the query set. */
  j = lst_size(query_gff->features)-1;
  for (i = (nmatches-1); i >= -1; i--) {
    if (i == -1) {
      for (k = j; k >= 0; k--) {
	fq = lst_get_ptr(query_gff->features, k);
	lst_push_ptr(target_gff->features, gff_new_feature_copy(fq));
      }
    } else {
      for (k = j; k > qmatch_idx[i]; k--) {
	fq = lst_get_ptr(query_gff->features, k);
	lst_push_ptr(target_gff->features, gff_new_feature_copy(fq)); 
      }
    }
    j = (qmatch_idx[i] - 1);
  }
  free(tmatch_idx);
  free(qmatch_idx);
}

/* Combine two paths into a single composite path containing features from
   both paths. E.g., for combining a path with predicted features with a
   reference path containing known features. The first path (path1) is
   considered the reference path -- motif features only will be
   used, while background features will be ignored. Background features in the
   composite path will be assigned based on the query path (path2). In cases
   where motif features are in identical or overlapping positions, the query
   path will be given priority (this ensures complete sampling of the posterior
   distribution, even when motifs are of known location and mode of selection).
*/
int* dms_composite_path(DMotifPhyloHmm *dm, int *path1, int *path2, 
			int seqlen, int offset, int force_priors) {
  int i, j, s1, s2, m1, m2, mc, *composite_path;

  composite_path = smalloc(seqlen * sizeof(int));
  for (i = 0; i < seqlen; i++) {
    s1 = path1[i];
    s2 = path2[1];
    m1 = dm->state_to_motifpos[s1];
    m2 = dm->state_to_motifpos[s2];
    
    if (s1 == 0 && s2 == 0) { /* Neutral BG in both paths */
      composite_path[i] = 0;
    } else if (s1 == 0 &&  s2 > 0) {
      /* Motif or selected BG in path2 only */
      composite_path[i] = path2[i];
    } else if (s1 > 0 && s2 == 0) {
      /* Motif or selected BG in path1 only */
      composite_path[i] = path1[i];
    } else { /* (s1 > 0 && s2 > 0) */
      /* Motif or selected BG in both paths -- priority based on invocation */
      composite_path[i] = force_priors ? path1[i] : path2[i];
    }
    /* Check for illegal transitions induced by path-flattening. This can be
       a result of either overlapping motif predictions or overlapping 
       predictions for background under different forms of selection. */
    if (mm_get(dm->phmm->hmm->transition_matrix, composite_path[i-1], 
	       composite_path[i]) == 0) {

      /* This is broken -- setting it up this way effectively forces the
	 prior path to be used if an overlapping feature begins earlier
	 in the sequence and negates prior_priority if the reference feature
	 begins later in the sequence */


      /* transition in composite_path is illegal -- use the state for [i]
	 from the path used to set [i-1], as this will always be legal */
      composite_path[i] = (composite_path[i-1] == path1[i-1]) ?
	path1[i] : path2[i];
    }
  }
  
  return(composite_path);
}

/* Reconstruct a state path from a set of GFF features. If the GFF contains only
   motifs, neutral background will be used in all spacer regions. */
int* dm_gff_to_path(DMotifPhyloHmm *dm, GFF_Set *gff, int seqlen, 
		     int offset) {
  int i, j, *path;
  GFF_Feature *f;
  
  /* Cycle through the indeces from 0 to seqlen and reconstruct the state path
     based on the presence/absence of a GFF feature. Motif features require
     special handling because of the constraint that all motif states must be
     cycled through sequentially. Where there is no motif or background
     feature specified, fill the path with 0's. */

  return(path);
}

/* Generate a set of state paths based on a GFF_Set and a Multi_MSA. */
int** dm_gff_to_paths(DMotifPhyloHmm *dm, GFF_Set *gff, Multi_MSA *blocks,
		      int offset) {
  int i, j, **paths;
  GFF_Set *tmp_gff;
  GFF_Feature *f;
  String *seqname;

/*   paths = smalloc(blocks->nblocks * sizeof(int*)); */
/*   for (i = 0; i < blocks->nblocks; i++) */
/*     paths[i] = smalloc(blocks->blocks[i]->length * sizeof(int)); */
  
  /* Cycle through (Multi_MSA)blocks->seqnames and stuff tmp_gff with all
     features from the gff that reside in a given sequence. Next call
     dm_gff_to_path to reconstruct the path and stuff the path into the
     set of paths for output. */

  return(paths);
}

/* Create a hashtable from a stored representation on disk */
Hashtable* dms_read_hash(FILE *hash_f, int nstates, int* nsamples) {
  int i, nelements, *counts;
  Hashtable *retval;
  Regex *ne_re = str_re_new("#[[:space:]]*NELEMENTS[[:space:]]*([0-9]+)");
  Regex *ns_re = str_re_new("#[[:space:]]*NSAMPLES[[:space:]]*([0-9]+)");
  String *line = str_new(STR_LONG_LEN);
  List *matches = lst_new_ptr(nstates+1);
  
  i = 0;
  while (str_readline(line, hash_f) != EOF) {
    lst_clear(matches);
    str_trim(line);
    if (line->length == 0) continue;
    
    if (str_re_match(line, ne_re, matches, 1) >= 0) {
      str_as_int(lst_get_int(matches, 1), &nelements);
      retval = hsh_new(nelements);
    } else if (str_re_match(line, ns_re, matches, 1) >= 0) {
      str_as_int(lst_get_ptr(matches, 1), nsamples);
    } else {
      str_split(line, NULL, matches);
      counts = smalloc(nstates * sizeof(int));
      for (i = 0; i < nstates; i++)
	str_as_int((String*)lst_get_ptr(matches, i), &counts[i]);
      hsh_put(retval, ((String*)lst_get_ptr(matches, 0))->chars, counts);
    }
  }
  str_re_free(ne_re);
  lst_free(matches);
  str_free(line);
  return(retval);
}

/* Write the contents of a hashtable to a file */
void dms_write_hash(Hashtable *path_counts, FILE *hash_f, int nstates, 
		    int nsamples) {
  int i, j, *counts;
  char *key;
  List *keys;
  
  keys = hsh_keys(path_counts);
  fprintf(hash_f, "# NELEMENTS %d\n", lst_size(keys));
  fprintf(hash_f, "# NSAMPLES %d\n", nsamples);
  for (i = 0; i < lst_size(keys); i++) {
    key = lst_get_ptr(keys, i);
    counts = hsh_get(path_counts, key);
    fprintf(hash_f, "%s\t", key);
    for (j = 0; j < nstates-1; j++)
      fprintf(hash_f, "%d\t", counts[j]);
    fprintf(hash_f, "%d\n", counts[nstates-1]);
  }
  lst_free(keys);
}

/* Fold features of two hashtables into a single unified hashtable */
void dms_combine_hashes(Hashtable *target, Hashtable *query, int nstates) {
  int i, j, *counts_t, *counts_q;
  List *keys;
  char *key;
  
  keys = hsh_keys(query);
  for (i = 0; i < lst_size(keys); i++) {
    key = lst_get_ptr(keys, i);
    if (hsh_get(target, key) != -1) {
      counts_t = hsh_get(target, key);
      counts_q = hsh_get(query, key);
      for (j = 0; j < nstates; j++) {
/* 	fprintf(stderr, "i %d, key %s, j %d, tc %d, qc %d, ", i, key, j,  */
/* 		counts_t[j], counts_q[j]); */
	counts_t[j] += counts_q[j];
/* 	fprintf(stderr, "tc_new %d\n", counts_t[j]); */
      }
    } else {
/*       fprintf(stderr, "i %d, key %s\n", i, key); */
      hsh_put(target, key, hsh_get(query, key));
    }
  }
}
