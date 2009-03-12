 
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
  int k = 2*nleaves - 3;                  /* no. branches (unrooted tree) */
  int nstates = (2 + 2*k) * (m->width+1); /* number of states */
  HMM *hmm;
  TreeModel **models = smalloc(nstates * sizeof(void*));
  TreeNode *n;
  CategoryMap *cm;
  List *inside = lst_new_ptr(source_mod->tree->nnodes),
    *outside = lst_new_ptr(source_mod->tree->nnodes);
  DMotifPhyloHmm *dm;
  double *alpha = (double*)smalloc(source_mod->tree->nnodes 
				   * sizeof(double)),
    *beta = (double*)smalloc(source_mod->tree->nnodes 
			     * sizeof(double)),
    *tau = (double*)smalloc(source_mod->tree->nnodes 
			    * sizeof(double)),
    *epsilon = (double*)smalloc(source_mod->tree->nnodes 
				* sizeof(double));
  assert(rho > 0 && rho < 1);
  assert(mu > 0 && 2*mu < 1);
  assert(nu > 0 && nu < 1);

  dm = smalloc(sizeof(DMotifPhyloHmm));
  dm->state_to_branch = (int*)smalloc(nstates 
				      * sizeof(int));
  dm->state_to_event = (int*)smalloc(nstates 
				     * sizeof(int));
  dm->state_to_motifpos = (int*)smalloc(nstates 
					* sizeof(int));
  dm->state_to_motif = (int*)smalloc(nstates 
				     * sizeof(int));
  dm->branch_to_states = (List**)smalloc(source_mod->tree->nnodes * sizeof(void*));
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
  dm->indel_mods = NULL;

  /* set up mappings between states and branches.  Branches are
     indexed by ids of child nodes; "birth" means functional element
     arose on branch, "death" means element was lost on branch. */
  
  for (i = 0; i < source_mod->tree->nnodes; i++)
    dm->branch_to_states[i] = lst_new_int(2 * m->width + 2);
  
  state = 0;
  for (i = 0; i < source_mod->tree->nnodes; i++){
    n = lst_get_ptr(source_mod->tree->nodes, i);

    for (pos = -1; pos < m->width; pos++) {
      /* note: i=0 implies root note and fully nonconserved element */
      dm->state_to_branch[state] = n->id;
      lst_push_int(dm->branch_to_states[n->id], state);
      dm->state_to_motifpos[state] = pos;
      dm->state_to_event[state] = (i == 0 ? NEUT : DEATH);
      state++;
    }

/*     fprintf(stderr, "mapping second half of model\n"); */

    for (pos = -1; pos < m->width; pos++) {
      /* note: i=0 implies root note and fully conserved element */
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
    
/*     fprintf(stderr, "done with second half of model\n"); */
  }


  assert(state == nstates);
			     
  /* Set up state to motif mapping. Maps states to first motif positions.
     States corresponding to internal motif or background positions are set
     to -1. */
  i = 0;
  for (state = 0; state < nstates; state++ ) {
    if (dm->state_to_motifpos[state] == 0)
      dm->state_to_motif[state] = i++;
    else
      dm->state_to_motif[state] = -1;
  }

 
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
    else { /* F81 model with motif freqs as equilibrium freqs */
      models[state] = tm_new(tr_create_copy(source_mod->tree), NULL,
                             vec_create_copy(m->probs[pos]), F81, 
			     DEFAULT_ALPHABET, 1, -1, NULL, -1);
      //      fprintf(stderr, "state %d pos %d ", state, pos);
      //      vec_print(m->probs[pos], stderr);
    }
    
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
  cm_free(cm); /* A copy of this is made within phmm_new, so it's safe to
		  free this one. */

  dm_set_transitions(dm);

  /* set up indel model, if necessary */
 if (alpha_c > 0) {
   dm->indel_mods = (DMotifIndelModel**)smalloc(nstates * sizeof(void*));
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
      if (e == DEATH)
        l = outside;            /* death */
      else
        l = inside;             /* birth, conserved or neutral */
      
      /* Set the params to conserved values unless we're at a completely
	 neutral state */
      if (e != NEUT) {
	for (i = 0; i < lst_size(l); i++) {
	  n = lst_get_ptr(l, i);
	  alpha[n->id] = alpha_c;
	  beta[n->id] = beta_c;
	  tau[n->id] = tau_c;
	  epsilon[n->id] = epsilon_c;
	}
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
  free(epsilon);
  return dm;
}

void dm_set_transitions(DMotifPhyloHmm *dm) {
  int s1, s2, p1, p2, b1, b2, e1, e2;
  double trans, rowsum;
  HMM *hmm = dm->phmm->hmm;
  //  TreeNode *nn;

  /* set transition probs based on mu, nu, and phi */
  for (s1 = 0; s1 < hmm->nstates; s1++) {
    p1 = dm->state_to_motifpos[s1];
    b1 = dm->state_to_branch[s1];
    e1 = dm->state_to_event[s1];
    rowsum = 0;

    for (s2 = 0; s2 < hmm->nstates; s2++) {
      p2 = dm->state_to_motifpos[s2];
      b2 = dm->state_to_branch[s2];
      e2 = dm->state_to_event[s2];
      trans = 0;

      /* Explicitly set params for each type of allowed transition */
      if ( e1 == NEUT ) { /* Transitions from NEUT states */
	if ( e2 == NEUT ) {
	  if (p1 == -1 && p2 == -1) /* Nb self-transition */
	    trans = ((1 - dm->zeta) * (1 - dm->nu));
	  else if (p1 == -1 && p2 == 0) /* Nb to Neutral-motif */
	    trans = (dm->zeta * (1 - dm->nu));
	  else if (p2 == (p1 + 1) || (p1 == dm->m->width-1 && p2 == -1))
	    /* transitions between motif states and from NMw to Nb*/
	    trans = 1;
	  else continue;
	} else if ( e2 == CONS ) {
	  if (p1 == -1 && p2 == -1) /* Nb to Cb */
	    trans = (dm->nu * (1 - dm->zeta) * (1 - dm->phi));
	  else if (p1 == -1 && p2 == 0) /* first cons-mtf position, from Nb */
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
        if (e2 == e1 && b2 == b1) { /* same selection and branch */
	  if (p1 == -1 && p2 == -1) /* background self-transitions */
	    trans = ((1 - dm->zeta) * (1 - dm->mu));
	  else if (p1 == -1 && p2 == 0) /* background to M1 */
	    trans = (dm->zeta * (1 - dm->mu));
	  else if (p1 == dm->m->width-1 && p2 == -1) /* Mw to B, same branch */
	    trans = (1 - dm->mu);
	  else if (p2 == (p1 + 1) && p2 < dm->m->width) /* transitions between
							motif positions */
	    trans = 1;
	  else continue;
	} else if (e2 == NEUT && (p1 == -1 || p1 == dm->m->width-1)
		   && p2 == -1) /* Last motif position or BG under selection
				   to Nb */
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

  /* Begin transitions */
  rowsum = 0;
  for (s1 = 0; s1 < hmm->nstates; s1++) {
    /* motif stationary probs */
    if (dm->state_to_motifpos[s1] == -1) /* background state */
      trans = 1 / (1 + dm->zeta);
    /*       trans = 1; */
    else if (dm->state_to_motifpos[s1] == 0) /* first motif position */
      trans = dm->zeta / (1 + dm->zeta);
    else /* internal motif position -- don't allow starts here! */
      trans = 0;
    
    /* multiply by dless stationary probs*/
    if (dm->state_to_event[s1] == NEUT)
      trans *= dm->mu/(dm->mu + dm->nu);
    
    else if (dm->state_to_event[s1] == CONS)
      trans *= dm->nu * (1-dm->phi) / (dm->mu + dm->nu);

    else
      trans *= dm->nu * dm->phi / ((dm->mu + dm->nu) * 2 * dm->k);

    vec_set(hmm->begin_transitions, s1, trans);
    rowsum += trans;
/*     fprintf(stderr, "s1 %d, trans %f, rowsum %f\n", s1, trans, rowsum); */
  }
  assert(fabs(rowsum-1) < 1e-6);  
  
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
  int i, j, l, mod;

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
/* 	fprintf(stderr, "i %d, n->id %d, char tuple %d, msa->is_missing %d\n", */
/* 		i, n->id, (int)ss_get_char_tuple(msa, i, dm->phmm->mods[0]->msa_seq_idx[n->id], 0), msa->is_missing[(int)ss_get_char_tuple(msa, i, dm->phmm->mods[0]->msa_seq_idx[n->id], 0)]); */
	
        if (msa->is_missing[(int)ss_get_char_tuple(msa, i, dm->phmm->mods[0]->msa_seq_idx[n->id], 0)])

          mark[n->id] = MISSING;
        else
          mark[n->id] = DATA;
      }

      else if (n == tree) { /* root */
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
    } /* End postorder traversal */

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
    } /* End preorder traversal */
  } /* End loop over tuples */

  /* now zero out all associated emissions probs */
  for (i = 0; i < msa->length; i++) {
    if (uninform[msa->ss->tuple_idx[i]] == NULL) continue;
    for (j = 0; j < lst_size(uninform[msa->ss->tuple_idx[i]]); j++) {
      mod = lst_get_int(uninform[msa->ss->tuple_idx[i]], j);
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
void dm_add_indel_emissions(DMotifPhyloHmm *dm, double **emissions, 
			    IndelHistory *ih) {
  int state, i;
  double *col_logl = (double*)smalloc(ih->ncols 
				      * sizeof(double));
  for (state = 0; state < dm->phmm->hmm->nstates; state++) {
    dmih_column_logl(ih, dm->indel_mods[state], col_logl);
    for (i = 0; i < ih->ncols; i++) {
/*       fprintf(stderr, "state %d, col %d, emissions[state][col] %f, col_logl[col] %f, ", */
/* 	      state, i, emissions[state][i], col_logl[i]); */
      emissions[state][i] += col_logl[i];
/*       fprintf(stderr, "emissions[state][col] %f\n", emissions[state][i]); */
    }
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
    dm->phmm->forward = (double**)smalloc(dm->phmm->hmm->nstates *
                                    sizeof(double*));
    for (i = 0; i < dm->phmm->hmm->nstates; i++)
      dm->phmm->forward[i] = (double*)smalloc(dm->phmm->alloc_len 
					      * sizeof(double));
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
  TreeNode *n;
  tm->alt_subst_mods = (AltSubstMod**)smalloc(tm->tree->nnodes
					      * sizeof(void*));
  for (i = 0; i < tm->tree->nnodes; i++) tm->alt_subst_mods[i] = NULL;

  for (i = 0; i < lst_size(nodelist); i++) {
    n = lst_get_ptr(nodelist, i);
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
  int i, j, k, nstates, p, e, b;
  TreeNode *node;
  String *type;
  CategoryMap *retval;
  CategoryRange *cr;

  nstates = (2 + 2*dm->k) * (dm->m->width+1); /* number of states */
  /*int ncats = (4 * dm->k + 4)-1; */
  retval = cm_new(nstates-1);

  for (i = 0; i < nstates; i++) {
    p = dm->state_to_motifpos[i];
    b = dm->state_to_branch[i];
    e = dm->state_to_event[i];
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
/*   fprintf(stderr, "i %d, j %d, nstates %d\n", i, j, nstates); */
  assert(j == nstates-1);
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
    snprintf(groupstr, STR_SHORT_LEN, "%s \"%s.%d\"", grouptag != NULL ? grouptag : "id",
	    idpref, groupno);
  else
    snprintf(groupstr, STR_SHORT_LEN, "%s \"%d\"", grouptag != NULL ? grouptag : "id", groupno);
  
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
	snprintf(groupstr, STR_SHORT_LEN, "%s \"%s.%d\"", grouptag != NULL ? grouptag : "id",
		idpref, groupno);
      else
	snprintf(groupstr, STR_SHORT_LEN, "%s \"%d\"", grouptag != NULL ? grouptag : "id",
		groupno);
    }
  }
  
  return gff;
}

List* dms_sample_paths(DMotifPhyloHmm *dm, PooledMSA *blocks,
		       double **tuple_scores, IndelHistory **ih,
		       List *seqnames, int max_seqlen, int bsamples,
		       int nsamples, int sample_interval,
		       int **priors, FILE *log,
		       GFF_Set *reference, int ref_as_prior,
		       int force_priors,
		       int quiet, char *cache_fname, int cache_int) {

  int i, j, k, l, *path, **trans, **ref_paths = NULL, nwins, reinit;
  unsigned long int seed;
  FILE *devrandom;
  double **forward_scores, llh;
  char *cache_out;
  GFF_Set *query_gff = NULL;
  Hashtable *path_counts;
  List *cache_files;

  path_counts = hsh_new(max((10 * lst_size(blocks->source_msas)), 10000));
  cache_files = lst_new_ptr(nsamples/cache_int);

  /* Allocate space for the forward scores and path. Set to the
     length of the longest sequence -- will reuse for each sequnce. */
  forward_scores = (double**)smalloc(dm->phmm->hmm->nstates * sizeof(double*));
  for (i = 0; i < dm->phmm->hmm->nstates; i++) {
    forward_scores[i] = (double*)smalloc(max_seqlen * sizeof(double));
  }
  path = (int*)smalloc(max_seqlen * sizeof(int));

  /* Allocate space for the transition counts and set to prior values */
  trans = (int**)smalloc(4 * sizeof(int*));
  for (i = 0; i < 4; i++) {
    trans[i] = (int*)smalloc(2 * sizeof(int));
    trans[i][0] = priors[i][0];
    trans[i][1] = priors[i][1];
  }

  /* If comparing to a reference GFF for logging purposes, allocate the
     stats array and compute the number of motif windows in the dataset. */
  nwins = 0;
  if (reference != NULL && log != NULL) {
    for (i = 0; i < lst_size(blocks->source_msas); i++)
      nwins += (blocks->lens[i] - (dm->m->width -1));
  }
  /* If using the reference set as prior information, reconstruct a set of
     state paths from the GFF_Set. NOT FULLY IMPLEMENTED!! */
  /*   if (reference != NULL && ref_as_prior == TRUE) { */
  /*     ref_paths = dm_gff_to_paths(dm, reference, blocks, -1); */
  /*   } */

  k = 0;
  l = 1;
  cache_out = (char*)smalloc(STR_MED_LEN * sizeof(char));

  /* Seed the random number generator based on the value of /dev/random */
  devrandom = fopen("/dev/random", "r");
  fread(&seed, sizeof(seed), 1, devrandom);
  srandom(abs(seed));
  fclose(devrandom);

  for (i = 0; i < (bsamples + nsamples); i++) {
    if (!quiet) fprintf(stderr, "*");
    /* Clean things up for logging motif stats for individual samples, if
       needed */
    if (log != NULL && reference != NULL) {
      if (query_gff != NULL)
	gff_free_set(query_gff);
      query_gff = gff_new_set();
    }

    /* Draw a set of params */
    dm->mu = beta_draw((double)trans[0][0], (double)trans[0][1]);
    dm->nu = beta_draw((double)trans[1][0], (double)trans[1][1]);
    dm->phi = beta_draw((double)trans[2][0], (double)trans[2][1]);
    dm->zeta = beta_draw((double)trans[3][0], (double)trans[3][1]);
    dm_set_transitions(dm);

    /* Reset transition counts to prior values */
    for (j = 0; j < 4; j++) {
      trans[j][0] = priors[j][0];
      trans[j][1] = priors[j][1];
    }

    llh = 0;
    for (j = 0; j < lst_size(blocks->source_msas); j++) {
      /* Zero out the emissions, forward and path structs to be safe */

      /* Fill in emissions table for this sequence */
      /*       fprintf(stderr, "j %d, blocks->lens[j] %d, ih[j] %p, forward_scores %p\
	       \n", j, */
      /*            blocks->lens[j], ih[j], forward_scores); */
      /*       fprintf(stderr, "seqname %s\n", ((String*)lst_get_ptr(seqnames, j))->c\
	       hars); */
      dms_lookup_emissions(dm, tuple_scores, dm->phmm->emissions,
			   blocks, j, blocks->lens[j],
			   (ih == NULL ? NULL : ih[j]));
      
      /* Run the forward algorithm */
      llh += hmm_forward(dm->phmm->hmm, dm->phmm->emissions, blocks->lens[j],
			 forward_scores);

      /* Run stochastic traceback to get the path */
      hmm_stochastic_traceback(dm->phmm->hmm, forward_scores, blocks->lens[j],
			       path);

      /* Determine if need to take a sample this iteration, else continue */
      if (i % sample_interval) {
	continue;
      } else {
	dms_count_transitions(dm, path, trans, blocks->lens[j],
			      ref_paths == NULL ? NULL : ref_paths[j],
			      force_priors);
	if (i >= bsamples) /* Past burn-in -- sample motifs */
	  dms_count_motifs(dm, path, blocks->lens[j], path_counts, j);
	if (log != NULL && reference != NULL) /* Prepare GFF for this sequence
						 and fold into query_gff for 
						 comparison to reference at 
						 end of sample */
	  dms_path_log(dm, path, blocks->lens[j], lst_get_ptr(seqnames, j),
		       query_gff);
      }
    }

    if (log != NULL) /* Write stats for this sample to the log file */
      dms_write_log(log, dm, trans, i, llh, query_gff, reference, nwins);

    /* Cache and reinitialize the hash every 100 samples after burn-in */
    if (i >= bsamples) {
      if (l == cache_int || i == nsamples+bsamples-1) {
	snprintf(cache_out, STR_MED_LEN, "%s.%d.tmp", cache_fname, k++);
	lst_push_ptr(cache_files, (void*)str_new_charstr(cache_out));
	/*      fprintf(stderr, "i %d, k %d, l %d, cache_file %s\n", i, k-1, l, cache_out); */
	if (i == nsamples+bsamples-1)
	  reinit = -1;
	else
	  reinit = max(10*lst_size(blocks->source_msas), 10000);
	path_counts = dms_cache_hash(path_counts, cache_out, (2*dm->k)+2, l,
				     reinit);
	l = 1;
      } else {
	l++;
      }
    }
  }

  fprintf(stderr, "\nDone.\n");

  /* Clean up our area */
  for (i = 0; i < dm->phmm->hmm->nstates; i++) {
    free(forward_scores[i]);
  }
  free(forward_scores);
  free(path);
  if (query_gff != NULL)
    gff_free_set(query_gff);
  free(cache_out);
  
  return cache_files;
}

List* dms_sample_paths_pthr(DMotifPhyloHmm *dm, PooledMSA *blocks,
			    double **tuple_scores, IndelHistory **ih,
			    List *seqnames, int bsamples, int nsamples, 
			    int sample_interval, int **priors, FILE *log,
			    GFF_Set *reference, int ref_as_prior, 
			    int force_priors, int quiet, char *cache_fname,
			    int cache_int, ThreadPool *pool, int nthreads) {
  
  int i, j, k, l, t, **trans, ***thread_trans, nwins, reinit, threads;
  unsigned long int seed;
  FILE *devrandom;
/*   int **ref_paths = NULL; */
  double llh, *thread_llh;
  char *cache_out;
  GFF_Set *query_gff = NULL;
  Hashtable *path_counts, **thread_counts;  /**< path_counts is the global hash
					       for all sequences. thread_counts
					       contains a separate hash for
					       each thread to ensure thread-
					       safe operations -- each of these
					       are folded into path_counts for
					       caching and final reporting. */
  List *cache_files, *work;
  DMsamplingThreadData *data;

  
  /* Allocate lists for the cache file names and thread work lists */
  cache_files = lst_new_ptr(max(nsamples/cache_int, 1));
  work = lst_new_ptr(lst_size(blocks->source_msas));
  lst_clear(work); /* to be safe */
  
  /* Simplifies looping for case where thread pool size == 0 */
  if (nthreads == 0)
    threads = 1;
  else
    threads = nthreads;
  
  /* Global hash initial allocation size */
  reinit = max(100*lst_size(blocks->source_msas), 1000);
  
  /* Allocate the global hashtable */
  path_counts = hsh_new(reinit);
  
  /* Allocate the thread_specific hashtables */
  thread_counts = (Hashtable**)smalloc(threads * sizeof(Hashtable*));
  for (i = 0; i < threads; i++) {
    thread_counts[i] = hsh_new(reinit/threads);
  }
  
  /* Allocate space for the global transition counts and set to prior values */
  trans = (int**)smalloc(4 * sizeof(int*));
  for (i = 0; i < 4; i++) {
    trans[i] = (int*)smalloc(2 * sizeof(int));
    trans[i][0] = priors[i][0];
    trans[i][1] = priors[i][1];
  }

  /* Allocate space for the single-thread transition counts. Initialize these
     to 0 as pseudocounts are added into the global struct! */
  thread_trans = (int***)smalloc(threads * sizeof(int**));
  for (i = 0; i < threads; i++) {
    thread_trans[i] = (int**)smalloc(4 * sizeof(int*));
    for (j = 0; j < 4; j++) {
      thread_trans[i][j] = (int*)smalloc(2 * sizeof(int));
      thread_trans[i][j][0] = thread_trans[i][j][1] = 0;
    }
  }
  
  /* Allocate space for the thread-specific log likelihoods and initialize
     to 0. */
  thread_llh = (double*)smalloc(threads * sizeof(double));
  for (i = 0; i < threads; i++) {
    thread_llh[i] = 0;
  }
    
  /* If comparing to a reference GFF for logging purposes, allocate the
     stats array and compute the number of motif windows in the dataset. */
  nwins = 0;
  if (reference != NULL && log != NULL) {
    for (i = 0; i < lst_size(blocks->source_msas); i++)
      nwins += (blocks->lens[i] - (dm->m->width -1));
  }
  /* If using the reference set as prior information, reconstruct a set of
     state paths from the GFF_Set. NOT FULLY IMPLEMENTED!! */
/*   if (reference != NULL && ref_as_prior == TRUE) { */
/*     ref_paths = dm_gff_to_paths(dm, reference, blocks, -1); */
/*   } */
  
  k = 0;
  l = 1;
  cache_out = (char*)smalloc(STR_MED_LEN * sizeof(char));

  /* Seed the random number generator based on the value of /dev/random */
  devrandom = fopen("/dev/random", "r");
  fread(&seed, sizeof(seed), 1, devrandom);
/*   fprintf(stderr, "seed %d\n", abs(seed)); */
  srandom(abs(seed));
  fclose(devrandom);
  
  /* Initialize the static "*l" Lists in hmm_max_or_sum_pthread */
  hmm_max_or_sum_pthread(dm->phmm->hmm, NULL, NULL, NULL, 0, 0, FORWARD,
			 nthreads, -1);
  
  for (i = 0; i < (bsamples + nsamples); i++) {
    if (!quiet) fprintf(stderr, "*");
    /* Clean things up for logging motif stats for individual samples, if
       needed */
    if (log != NULL && reference != NULL) {
      if (query_gff != NULL)
	gff_free_set(query_gff);
      query_gff = gff_new_set();
    }
    
    /* Draw a set of params */
    dm->mu = beta_draw((double)trans[0][0], (double)trans[0][1]);
    dm->nu = beta_draw((double)trans[1][0], (double)trans[1][1]);
    dm->phi = beta_draw((double)trans[2][0], (double)trans[2][1]);
    dm->zeta = beta_draw((double)trans[3][0], (double)trans[3][1]);
    dm_set_transitions(dm);

    /* Set transition_score matrix in the hmm */
    hmm_set_transition_score_matrix(dm->phmm->hmm);

    /* Reset transition counts to prior values */
    for (j = 0; j < 4; j++) {
      trans[j][0] = priors[j][0];
      trans[j][1] = priors[j][1];
    }

    /* Set up the work list for multithreading over sequences. First make sure
       the existing list is empty and all old work items have been freed. */
    if (!lst_empty(work)) {
      for (j = 0; j < lst_size(blocks->source_msas); j++) {
	data = lst_get_ptr(work, j);
	free(data);
      }
      lst_clear(work);
    }
    for (j = 0; j < lst_size(blocks->source_msas); j++) {
      data = (DMsamplingThreadData*)smalloc(sizeof(DMsamplingThreadData));
      
      data->dm = dm;
      data->blocks = blocks;
      data->ih = (ih == NULL ? NULL : ih[j]);
      data->tuple_scores = tuple_scores;

      data->thread_llh = thread_llh;
      data->thread_trans = thread_trans;
      data->thread_counts = thread_counts;
      
      data->do_sample = (i < bsamples ? 0 : (i % sample_interval ? 0 : 1));
      data->seqnum = j;
      data->seqname = ((char*)((String*)lst_get_ptr(seqnames, j))->chars);
      data->log = log;
      data->query_gff = query_gff;
      data->do_reference = (reference == NULL ? 0 : (log == NULL ? 0 : 1));
      data->nthreads = nthreads;
      data->p = pool;
      data->sample = i;
      
/*       fprintf(stderr, "i %d, j %d, t %d, data->llh %p, data->trans %p, data->path_counts %p, data->do_sample %d\n", */
/* 	      i, j, t, data->llh, data->trans, data->path_counts, */
/* 	      data->do_sample); */

      lst_push(work, (void*)(&data));
    }
    
    thr_foreach(pool, work, dms_launch_sample_thread);

    /* Combine transition counts from all threads into global counts and
       reset thread-specific counters. Also total up global llh and reset
       thread-specific vals to 0. */
    llh = 0;
    for (t = 0; t < threads; t++) {
      for (j = 0; j < 4; j++) {
	trans[j][0] += thread_trans[t][j][0];
	trans[j][1] += thread_trans[t][j][1];
	thread_trans[t][j][0] = thread_trans[t][j][1] = 0;
      }
      llh += thread_llh[t];
      thread_llh[t] = 0;
    }
    
    if (log != NULL) /* Write stats for this sample to the log file */
      dms_write_log(log, dm, trans, i, llh, query_gff, reference, nwins);
    
    /* Cache and reinitialize the hash periodically after burn-in */
    if (i >= bsamples) {
      if (l == cache_int || i == nsamples+bsamples-1) {
	/* Fold together hashes from individual threads */
	for (t = 0; t < threads; t++) {
	  dms_combine_hashes(path_counts, thread_counts[t], (2*dm->k)+2);
	  hsh_clear(thread_counts[t]);
	}
	/* build and store a temp file name to hold the cached data */
	snprintf(cache_out, STR_MED_LEN, "%s.%d.tmp", cache_fname, k++);
	lst_push_ptr(cache_files, (void*)str_new_charstr(cache_out));
/* 	fprintf(stderr, "i %d, k %d, l %d, cache_file %s\n", i, k-1, l, cache_out); */
	if (i == nsamples+bsamples-1) /* Don't reinit the hash if we're 
					 done with sampling. */
	  reinit = -1;
	path_counts = dms_cache_hash(path_counts, cache_out, (2*dm->k)+2, l, 
				     reinit);
	l = 1;
      } else {
	l++;
      }
    }
  }

  fprintf(stderr, "\n\tSampling complete.\n");

  /* Clean up our area */
  if (query_gff != NULL)
    gff_free_set(query_gff);

  for (i = 0; i < lst_size(blocks->source_msas); i++) {
    data = lst_get_ptr(work, i);
    free(data);
  }
  lst_free(work);

  for (i = 0; i < threads; i++) {
    for (j = 0; j < 4; j++) {
      free(thread_trans[i][j]);
    }
    free(thread_trans[i]);
    hsh_free_with_vals(thread_counts[i]);
  }
  free(thread_trans);
  free(thread_counts);
  free(thread_llh);
  free(cache_out);
  
  return cache_files;
}

/* This will be used to run individual sampling threads, called by
   dms_launch_sample_thread. */
void dms_sample_path(DMotifPhyloHmm *dm, PooledMSA *blocks, IndelHistory *ih, 
		     double **tuple_scores, double *thread_llh, 
		     int ***thread_trans, Hashtable **thread_counts, 
		     int do_sample, int seqnum, char *seqname, FILE *log, 
		     GFF_Set *query_gff, int do_reference, int nthreads,
		     ThreadPool *pool, int sample) {

  int i,*path, thread_idx, **trans;
  double **emissions, **forward_scores, *llh;
  MSA *msa = lst_get_ptr(blocks->source_msas, seqnum);
  Hashtable *path_counts;

  /* Get the index of the current thread and appropriate vectors in thread-
     specific structures. */
  thread_idx = thr_index(pool);
  
  if (thread_idx == -1)
    die("ERROR: Attempt to access out-of-bounds thread in dms_sample_path!\n");

  trans = thread_trans[thread_idx];
  llh = &(thread_llh[thread_idx]);
  path_counts = thread_counts[thread_idx];
  
  /* Allocate space for the forward scores, emissions and path. Set to the
     length of the longest sequence -- will reuse for each sequnce. */
  forward_scores = (double**)smalloc(dm->phmm->hmm->nstates * sizeof(double*));
  emissions = (double**)smalloc(dm->phmm->hmm->nstates * sizeof(double*));
  for (i = 0; i < dm->phmm->hmm->nstates; i++) {
    forward_scores[i] = (double*)smalloc(msa->length * sizeof(double));
    emissions[i] = (double*)smalloc(msa->length * sizeof(double));
  }
  path = (int*)smalloc(msa->length * sizeof(int));
  
  dms_lookup_emissions(dm, tuple_scores, emissions, blocks, seqnum,
		       msa->length, ih);
  
  /* Run the forward algorithm */
  *llh += hmm_forward_pthread(dm->phmm->hmm, emissions, msa->length,
			      forward_scores, nthreads, thread_idx);
  /* Run stochastic traceback to get the path */
/*   fprintf(stderr, "Thread: %d, Sequence: %s, *l %p\n",  */
/* 	  thread_id, seqname, l); */
  hmm_stochastic_traceback(dm->phmm->hmm, forward_scores, msa->length, path);
/*   fprintf(stderr, "Thread: %d done with sequence %s\n", */
/* 	  thread_id, seqname); */
  
  dms_count_transitions(dm, path, trans, msa->length,
			NULL, FALSE);
  
  if (do_sample == 1) { /* Past burn-in -- sample motifs */
    dms_count_motifs(dm, path, msa->length, path_counts, seqnum);
    if (log != NULL && do_reference != 0) /* Prepare GFF for this sequence
					     and fold into query_gff for
					     comparison to reference at
					     end of sample */
      dms_path_log(dm, path, msa->length, seqname, query_gff);
/*     dms_dump_sample_data(sample, thread_idx, seqname, msa->length, path, trans, */
/* 			 path_counts, stderr, 2*dm->k+2); */
  }
  
  /* Clean up our area */
  for (i = 0; i < dm->phmm->hmm->nstates; i++) {
    free(forward_scores[i]);
    free(emissions[i]);
  }
  free(forward_scores);
  free(emissions);
  free(path);  
}

/* Worker function to launch sampling threads */
void dms_launch_sample_thread(void *data) {
  DMsamplingThreadData *thread = *((void**)data);

  dms_sample_path(thread->dm, thread->blocks, thread->ih, 
		  thread->tuple_scores, thread->thread_llh, 
		  thread->thread_trans, thread->thread_counts,
		  thread->do_sample, thread->seqnum, thread->seqname, 
		  thread->log, thread->query_gff, thread->do_reference, 
		  thread->nthreads, thread->p, thread->sample);  
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

GFF_Feature* dms_motif_as_gff_feat(DMotifPhyloHmm *dm, PooledMSA *blocks,
				   List *seqnames, char *key, int *counts,
				   int nsamples, int sample_interval,
				   int refidx) {
  int i, cat, seqnum, start, end, best, mtf_total, b_total, d_total,
    c_total, n_total, m, width;
/*   int j, all_missing; */
  char strand[2], delim[2], post[STR_LONG_LEN], *seqname;
/*   char *candidate; */
  CategoryMap *cm = dm->phmm->cm;
  String *attributes, *key_string, *tmp_str;
  List *split = lst_new_ptr(3);
  GFF_Feature *f;
  MSA *msa; 
  /* MSA *sub_msa; */
  Regex *rc_re = str_re_new("[A-Za-z0-9.]+_rc");

  key_string = str_new_charstr(key);
  attributes = str_new(STR_LONG_LEN);

  snprintf(delim, 2, "_");
  str_split(key_string, delim, split);
  str_free(key_string);
  str_as_int(lst_get_ptr(split, 0), &seqnum);
  str_as_int(lst_get_ptr(split, 1), &start);
  str_as_int(lst_get_ptr(split, 2), &end);
  if (end >= blocks->lens[seqnum])
    end = (blocks->lens[seqnum] - 1);
  seqname = ((char*)((String*)lst_get_ptr(seqnames, seqnum))->chars);

/*   fprintf(stderr, "key %s, seqnum %d, start %d, end %d, &seqname %p, seqname %s\n",  */
/* 	  key, seqnum, start, end, lst_get_ptr(seqnames, seqnum), seqname); */
  
  snprintf(post, STR_LONG_LEN, "ID \"%s.%d\"; ", seqname, start);
  str_append_charstr(attributes, post);

  /* First pass counts the total number of motifs predicted in this position
     and totals for conserved, gain, loss and neutral categories of motifs,
     along with finding the "best" category -- that is, the one with the
     highest posterior probability -- that will be used as the GFF feature
     type. */
  mtf_total = c_total = b_total = d_total = n_total = 0;
  best = 1;
  for (i = 0; i < dm->phmm->hmm->nstates; i++) {
    m = dm->state_to_motif[i];

    if (m == -1 || counts[m] == 0)
      continue;
    
    mtf_total += counts[m];
    if (dm->state_to_event[i] == CONS)
      c_total += counts[m];
    else if (dm->state_to_event[i] == BIRTH)
      b_total += counts[m];
    else if (dm->state_to_event[i] == DEATH)
      d_total += counts[m];
    else /* dm->state_to_event[i] == NEUT */
      n_total += counts[m];
    
    if (counts[m] > counts[ (dm->state_to_motif[best] == -1 ? 
			     0 : dm->state_to_motif[best]) ]) {
      best = i;      
    }
  }
  
/*   fprintf(stderr, "nsamples %d, mtf_total %d, c_total %d, b_total %d, d_total %d, n_total %d\n", nsamples, mtf_total, c_total, b_total, d_total, n_total); */

  /* Compute marginal posteriors for total motifs and motifs of each flavor */
  snprintf(post, STR_LONG_LEN, "MP %s %f; MP %s %f; MP %s %f; MP %s %f; MP %s %f; ",
	  "MOTIF", ((double)mtf_total / (double)(nsamples / sample_interval)),
	  "CONS", ((double)c_total / (double)mtf_total),
	  "BIRTH", ((double)b_total / (double)mtf_total),
	  "DEATH", ((double)d_total / (double)mtf_total),
	  "NEUT", ((double)n_total / (double)mtf_total));
  str_append_charstr(attributes, post);
  
  /* Second pass computes posterior probabilities of being in each flavor of
     motif given that there is a motif in this position and appends all to the
     GFF "attributes" field. */
  for (i = 0; i < dm->phmm->hmm->nstates; i++) {
    m = dm->state_to_motif[i];

    if (m == -1 || counts[m] == 0)
      continue;

    cat = cm->ranges[dm->phmm->state_to_cat[i]]->start_cat_no;
    snprintf(post, STR_LONG_LEN, "PP %s %f; ",
	     ((String*)cm_get_feature(cm, cat))->chars,
	     ((double)counts[m] / (double)mtf_total));
    str_append_charstr(attributes, post);
    /*       printf("%s %d %s %s %d %d %s %s\n", key, counts[i], seqname,  */
    /* 	     ((String*)cm_get_feature(cm, cat))->chars, start, end, strand,  */
    /* 	     attributes->chars); */
  }
  
  /*   fprintf(stderr, "seqname: %s, attributes: %s\n", seqname, attributes->chars); */
  
  /* Get the category for the feature with the highest PP */
  cat = cm->ranges[dm->phmm->state_to_cat[best]]->start_cat_no;

  /* Get the strand for the feature -- this is always + if we're not dealing
     with reverese complemented sequences. Also adjust coordinates if we're
     dealing with a feature on the - strand */
  tmp_str = str_new_charstr(seqname);
  if (str_re_match(tmp_str, rc_re, NULL, 1) >= 0) {
    sprintf(strand, "%s", "-");
    msa = lst_get_ptr(blocks->source_msas, seqnum);
    width = end - start;
    start = (msa->length - end) + 1;
    end = start + width;
    seqname = ((char*)((String*)lst_get_ptr(seqnames, seqnum-1))->chars);
  } else {
    sprintf(strand, "%s", "+");
  }
  str_free(tmp_str);

  /* Get the subalignment representing the current motif. This will need seq
     strings rebuilt if using SS as the input format. Padd with 5 bases up and
     downstream. -- NOTE: this is now handled by the browser display code! */
/*   msa = lst_get_ptr(blocks->source_msas, seqnum); */
/*   sub_msa = msa_sub_alignment(msa, NULL, 0, (start-5 < 0 ? 0 : start-5), */
/* 			      (end+5 > msa->length ? msa->length : end+5)); */
/*   if (sub_msa->seqs == NULL) */
/*     ss_to_msa(sub_msa); */

/*   /\* Stringify the sequences. Use comma delimiting. *\/ */
/*   tmp_str = str_new(STR_LONG_LEN); */
/*   str_append_charstr(tmp_str, "SEQS \""); */
/*   for (i = 0; i < sub_msa->nseqs; i++) { */
/*     /\* Do not include strings of all missing data characters *\/ */
/*     all_missing = 1; */
/*     candidate = sub_msa->seqs[i]; */
/*     for (j = 0; j < strlen(candidate); j++) { */
/*       if (!msa->is_missing[(int)candidate[j]]) { */
/* 	all_missing = 0; */
/* 	break; */
/*       } */
/*     } */
    
/*     if (all_missing) */
/*       continue; */
    
/*     str_append_charstr(tmp_str, sub_msa->names[i]); */
/*     str_append_char(tmp_str, ':'); */
/*     str_append_charstr(tmp_str, sub_msa->seqs[i]); */
/*     str_append_char(tmp_str, ','); */
/*   } */
/*   str_append_charstr(tmp_str, "\"; "); */
/*   str_append(attributes, tmp_str); */
/*   str_free(tmp_str); */
/*   msa_free(sub_msa); */

  /* Construct the GFF feature */
  f = gff_new_feature(str_new_charstr(seqname),
		      str_new_charstr("DMSAMPLE"),
		      str_dup(cm_get_feature(cm, cat)), start, end,
		      ((double)mtf_total /
		       (double)(nsamples / sample_interval)),
		      *strand, GFF_NULL_FRAME,
		      attributes, TRUE);
  
/*   gff_print_feat(stderr, f); */
  /* Adjust coordinates to the reference sequence, if needed */
  if (refidx != 0)
    dms_map_gff_coords(blocks, seqnum, f, 0, refidx);

  lst_free(split);
  str_re_free(rc_re);
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
  tmatch_idx = (int*)smalloc(lst_size(target_gff->features) 
			     * sizeof(int));
  qmatch_idx = (int*)smalloc(lst_size(query_gff->features) 
			     * sizeof(int));

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
  free(tmatch_idx);
  free(qmatch_idx);
}

/* Create a hashtable from a stored representation on disk */
Hashtable* dms_read_hash(FILE *hash_f, /* 2*dm->k+2 */ int nstates, 
			 int* nsamples) {
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
      str_as_int(lst_get_ptr(matches, 1), &nelements);
      retval = hsh_new(nelements);
    } else if (str_re_match(line, ns_re, matches, 1) >= 0) {
      str_as_int(lst_get_ptr(matches, 1), nsamples);
    } else {
      str_split(line, NULL, matches);
      counts = (int*)smalloc(nstates * sizeof(int));
      /* Count entries in matches are offset by 1 due to the key occupying
	 index 0 */
      for (i = 1; i <= nstates; i++)
	str_as_int((String*)lst_get_ptr(matches, i), &counts[i-1]);
      hsh_put(retval, ((String*)lst_get_ptr(matches, 0))->chars, counts);
    }
  }
  str_re_free(ne_re);
  lst_free(matches);
  str_free(line);
  return(retval);
}

/* Write the contents of a hashtable to a file */
void dms_write_hash(Hashtable *path_counts, FILE *hash_f, 
		    /* 2*dm->k+2 */ int nstates,
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

/* Fold features of two hashtables into a single unified hashtable, using
   the target hash as the destination for the combined hash. Values in the
   query hash will be unchanged. */
void dms_combine_hashes(Hashtable *target, Hashtable *query,
			/* 2*dm->k+2 */ int nstates) {
  int i, j, *counts_t, *counts_q;
  List *keys;
  char *key;
  
  keys = hsh_keys(query);
  for (i = 0; i < lst_size(keys); i++) {
    key = lst_get_ptr(keys, i);
    counts_q = hsh_get(query, key);
    counts_t = hsh_get(target, key);

    if (counts_t == (void*)-1) { /* Query key not found in destination hash. 
				    Allocate and zero storage for the counts
				    and insert in the destination hash. */
      counts_t = (int*)smalloc(nstates * sizeof(int));
      for (j = 0; j < nstates; j++)
	counts_t[j] = 0;
      hsh_put(target, key, counts_t);
    }

    /* Loop over states and add counts from the query values for this key to
       the destination values for this key. */
    for (j = 0; j < nstates; j++) {
/* 	fprintf(stderr, "i %d, key %s, j %d, tc %d, qc %d, ", i, key, j, */
/* 		counts_t[j], counts_q[j]); */
      counts_t[j] += counts_q[j];
/* 	fprintf(stderr, "tc_new %d\n", counts_t[j]); */
    }
  }
}

/** Compute emissions for given PhyloHmm and MSA using posix threads. Based on
    phmm_compute_emissions, but does not allocate any memory, so all must be 
    allocated externally. */
void dms_compute_emissions(PhyloHmm *phmm,
                                /**< Initialized PhyloHmm */
			   MSA *msa,
                                /**< Source alignment */
			   int quiet,
			   /**< Determins whether progress is
                                   reported to stderr */
			   ThreadPool *pool,
			   int nthreads
                            ) {
  int i, mod;
  DMemissionsThreadData *data;
  List *work = lst_new_ptr(phmm->hmm->nstates);
  /* build the list of calls to dm_compute_log_likelihood to thread out */
  for (i = phmm->hmm->nstates-1; i >= 0; i--) {
    data = (DMemissionsThreadData*)smalloc(sizeof(DMemissionsThreadData));

    mod = phmm->state_to_mod[i];
    data->mod = phmm->mods[mod];
    data->msa = msa;
    data->emissions_vec = phmm->emissions[i];
    data->cat = -1;
    data->state = i+1;
    data->nstates = phmm->hmm->nstates;
    data->quiet = quiet;
    lst_push(work, (void*)(&data));

/*     fprintf(stderr, "data %p %p  %p  %p  %d\n", data, data->mod, data->msa, data->emissions_vec, data->cat); */
    
/*     data = *((void**)lst_get(work, i)); */
/*     fprintf(stderr, "data %p %p  %p  %p  %d\n", data, data->mod, data->msa, data->emissions_vec, data->cat); */
  }
  
  if (nthreads == 0)
    lst_reverse(work);
  
  thr_foreach(pool, work, dms_do_emissions_row);
  for (i = 0; i < phmm->hmm->nstates; i++) {
    data = lst_get_ptr(work, i);
    free(data);
  }
  lst_free(work);
}

/* Launch a thread to compute a single emissions row */
void dms_do_emissions_row(void *data) {
  DMemissionsThreadData *thread = *((void**)data);
  
  if (!thread->quiet) {
    fprintf(stderr, "\tState %d  (%d remaining)\n", 
	    thread->state, (thread->nstates - thread->state));
  }
/*   fprintf(stderr, "thread id %d, data %p, mod %p, msa %p, emissions row %p, cat %d\n", */
/* 	  (int)pthread_self(), thread, thread->mod, thread->msa,  */
/* 	  thread->emissions_vec, thread->cat); */
  dm_compute_log_likelihood(thread->mod, thread->msa, thread->emissions_vec,
			  thread->cat);
}

/* Fill in an emissions matrix for a single sequence based on a precomputed
   table of values for each tuple. */
void dms_lookup_emissions(DMotifPhyloHmm *dm, double **tuple_scores,
			  double **emissions, PooledMSA *blocks, 
			  int seqnum, int seqlen, IndelHistory *ih) {
  int i, j, tuple_idx, tuple;
  MSA *smsa = lst_get_ptr(blocks->source_msas, seqnum);
  
  for (i = 0; i < dm->phmm->hmm->nstates; i++) {
    for (j = 0; j < seqlen; j++) {
      tuple = smsa->ss->tuple_idx[j];
      tuple_idx = blocks->tuple_idx_map[seqnum][tuple];
/*       fprintf(stderr, "seqnum %d, seqlen %d, i %d, j %d, tuple %d, tuple_idx %d\n", */
/* 	      seqnum, seqlen, i, j, tuple, tuple_idx); */
      emissions[i][j] = tuple_scores[i][tuple_idx];
/*       fprintf(stderr, "i %d, j %d, emissions[i][j] %f\n", */
/* 	      i, j, emissions[i][j]); */
    }
  }
  /*   Adjust for indels, if needed */
  if (ih != NULL)
    dm_add_indel_emissions(dm, emissions, ih);
}

/* Count transitions between states for parameter draws */
void dms_count_transitions(DMotifPhyloHmm *dm, int *path, int **trans,
			   int seqlen, int *ref_path, int force_priors) {
  int i, s1, s2, p1, p2, e1, e2/*, *tmp_path */;

  /* If using reference set as priors, build a composite path containing
     all known and predicted features and use this as the path */
  /* NOT YET IMPLEMENTED */
/*   if (ref_path != NULL) { */
/*     tmp_path = dms_composite_path(dm, ref_path, path, */
/* 				  seqlen, 0, force_priors); */
/*     for (i = 0; i < seqlen; i++) */
/*       path[i] = tmp_path[i]; */
/*     free(tmp_path); */
/*   } */
  
  s1 = BEGIN_STATE;
  for (i = 0; i < seqlen; i++) {
    s2 = path[i];
    
    /* 	  fprintf(stderr, "i %d, len %d, s2 %d\n", i, seqlen, s2); */

    /* 2-11-2009 -- After discussions with ACS, seems the best thing to do is
       ignore begin transitions -- entries into these states are poorly modeled
       in this context and don't give much information as to the actual param
       values. */
    if (s1 == BEGIN_STATE) {
      s1 = s2;
      continue;
    }

    /* Begin transitions are counted as if they are from the neutral
       background state */
    p1 = (s1 == BEGIN_STATE ? -1 : dm->state_to_motifpos[s1]);
    e1 = (s1 == BEGIN_STATE ? NEUT : dm->state_to_event[s1]);
    p2 = dm->state_to_motifpos[s2];
    e2 = dm->state_to_event[s2];
    
    /* Track different types of transitions for param estimation */
    if (p1 == -1) { /* Transitions from bg states */
      if (e1 == NEUT) { /* Neutral BG to...  */
	if (e2 != NEUT) { /* selected state -- affects alpha count for nu */
	  trans[1][0]++;
	  if (e2 == BIRTH || e2 == DEATH) { /* lineage-specific state --
					       also affects alpha count 
					       for phi */
	    trans[2][0]++;
	  } else { /* conserved state -- also affects beta count for phi */
	    trans[2][1]++;
	  }
	} else { /* Transition from neutral bg to neutral state (bg or motif)
		    -- affects beta count for nu */
	  trans[1][1]++;
	}
      } else { /* Source state is under selection (can be cons, birth 
		  or death) -- count transitions into... */
	if (e2 == NEUT) { /* neutral states -- Affects alpha count for mu */
	  trans[0][0]++;
	} else { /* selected states -- Affects beta count for mu */
	  trans[0][1]++;
	}
      }
      /* Motif entries can be counted without regard to selection -- they only
	 affect zeta.... */
      if (p2 == 0) { /* BG into a motif -- affects alpha count for zeta */
	trans[3][0]++;
      } else { /* into another background category -- affects beta count for
		  zeta. */
	trans[3][1]++;
      }
      
    } else if (p1 == (dm->m->width - 1)
	       && e1 != NEUT) { /* Transitions from motif end */
      if (e2 == NEUT) { /* into neutral background -- affects alpha count for 
			   mu */
	trans[0][0]++;
      } else { /* into the background state associated with the current branch
		  -- affects beta count for mu */
	trans[0][1]++;
      }
    }
    s1 = s2;
  }
}


void dms_count_motifs(DMotifPhyloHmm *dm, int *path, int seqlen, 
		      Hashtable *path_counts, int seqnum) {
  int i, j, k, p, s, *counts, m, w, end;
  char hsh_key[STR_SHORT_LEN];
  m = (2*dm->k)+2;
  w = dm->m->width;
  
  for (i = 0; i < seqlen; i++) {
    end = -1;
    p = dm->state_to_motif[path[i]];
    
    /* Position 0 in the path needs special handling to cover the case where
       a path starts with an internal motif position. */
    if (i == 0 && p == -1) {
      s = dm->state_to_motifpos[path[i]];
      if (s == -1) /* background state -- no need to proceed */
	continue;
      else {
	for (k = i+1; k < (i + w) && s >= 0; k++) {
	  s = dm->state_to_motifpos[path[k]];
	}
	end = k;
      }
    }

    assert((end - i) <= w);
    
    if (p != -1 || end != -1) {
      if (end == -1)
	end = i + w;
      snprintf(hsh_key, STR_SHORT_LEN, "%d_%d_%d", seqnum, i, end);
      counts = hsh_get(path_counts, hsh_key);
      if (counts == (void*)-1) { /* First motif at this position */
	counts = (int*)smalloc(m * sizeof(int));
	for (j = 0; j < m; j++)
	  counts[j] = 0;
	counts[p] = 1;
	hsh_put(path_counts, hsh_key, counts);
      } else {
	counts[p]++;
      }
    }
    /* 	    fprintf(stderr, "key %s state %d count %i\n", hsh_key, path[i], */
    /* 		    ((int*)hsh_get(path_counts, hsh_key))[s]); */
  }
}

void dms_path_log(DMotifPhyloHmm *dm, int *path, int seqlen, char *seqname,
		  GFF_Set *motifs) {
  int i;
  GFF_Set *tmp_gff;
  GFF_Feature *f;

  tmp_gff = gff_new_set();	
  tmp_gff = dm_labeling_as_gff(dm->phmm->cm, path, seqlen,
			       dm->m->width,
			       dm->phmm->state_to_cat,
			       dm->state_to_motifpos,
			       dm->phmm->reverse_compl, seqname, "DMSAMPLE", 
			       NULL, NULL, NULL);
  for (i = 0; i < lst_size(tmp_gff->features); i++) {
    f = lst_get_ptr(tmp_gff->features, i);
    lst_push_ptr(motifs->features, gff_new_feature_copy(f));
  }
  gff_free_set(tmp_gff);  
}

void dms_write_log(FILE *log, DMotifPhyloHmm *dm, int **trans, int sample, 
		   double llh, GFF_Set *query_gff, GFF_Set *reference, 
		   int nwins) { 
  int i, *stats;
  double fpr, spc, ppv;
  char param[STR_SHORT_LEN];
 
  fprintf(log, "Sample = %d, total log likelihood = %f\n", sample, llh);
  
  /* If using a reference set, compare known and predicted features and 
     keep track of the matches, mismatches, etc. in the stats structure */
  if (reference != NULL) {
    fpr = spc = ppv = 0;
    stats = (int*)smalloc(4 * sizeof(int));
    for (i = 0; i < 4; i++)
      stats[i] = 0;
    
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
  
  for (i = 0; i < 4; i++) {
    if (i == 0)
      snprintf(param, STR_SHORT_LEN, "mu");
    else if (i == 1)
      snprintf(param, STR_SHORT_LEN, "nu");
    else if (i == 2)
      snprintf(param, STR_SHORT_LEN, "phi");
    else /* (i == 3) */
      snprintf(param, STR_SHORT_LEN, "zeta");
    
    fprintf(log, "%s transitions: alpha = %d, beta = %d\n", param, 
	    trans[i][0], trans[i][1]);
  }
  
  fprintf(log, "mu = %f, nu = %f, phi = %f, zeta = %f\n\n", dm->mu, 
	  dm->nu, dm->phi, dm->zeta);
}

DMotifPmsaStruct *dms_read_alignments(FILE *F, int do_ih, int quiet, 
				      int revcomp) {
  
  Regex *blocks_re = str_re_new("#[[:space:]]*BLOCKS[[:space:]]*=[[:space:]]*([0-9]+)");
  Regex *alph_re = str_re_new("#[[:space:]]*ALPHABET[[:space:]]*=[[:space:]]*([A-Z]+)");
  Regex *format_re = str_re_new("#[[:space:]]*FORMAT[[:space:]]*=[[:space:]]*([A-Z]+)");
  
  int i, j, nblocks, ncats, tuple_size;
  char msa_fname[STR_MED_LEN], ih_fname[STR_MED_LEN];
  FILE *msa_f, *ih_f;
  msa_format_type format;
  List *msas = NULL, *matches = lst_new_ptr(4);
  MSA *msa, *msa_rc;
  String *line = str_new(STR_MED_LEN), *alphabet = NULL, 
    *fname = str_new(STR_MED_LEN);
  DMotifPmsaStruct *dmpmsa;

  dmpmsa = (DMotifPmsaStruct*)smalloc(sizeof(DMotifPmsaStruct));
  
  nblocks = i = dmpmsa->max_seqlen = 0;

  while (str_readline(line, F) != EOF) {

    /* Free all substrings in the matches list or there will be a memory 
       leak */
    for (j = 0; j < lst_size(matches); j++)
      str_free((String*)lst_get_ptr(matches, j));

    str_trim(line);
    if (line->length == 0) continue; /* ignore blank lines */
    
    if (msas == NULL) {
      if (str_re_match(line, blocks_re, matches, 1) >= 0) { 
	str_as_int(lst_get_ptr(matches, 1), &nblocks);
	if (revcomp == TRUE)
	  nblocks *= 2;
	i++;
      } else if (str_re_match(line, alph_re, matches, 1) >= 0) {
	alphabet = str_dup(lst_get_ptr(matches, 1));
	i++;
      } else if (str_re_match(line, format_re, matches, 1) >= 0) {
	format = msa_str_to_format(((String*)lst_get_ptr(matches, 1))->chars);
	i++;
      } else {
	die("Bad header in alignment list file\n");
      }
      
      if (nblocks != 0 && alphabet != NULL && i == 3) {
	msas = lst_new_ptr(nblocks);
	dmpmsa->seqnames = lst_new_ptr(nblocks);
	if (do_ih) { /* Allocate space for indel histories, if
			called for -- this will need to be filled
			in explicitly later */
	  dmpmsa->ih = smalloc(nblocks * sizeof(void*));
	} else {
	  dmpmsa->ih = NULL;
	}
	i = 0;
      }
    } else {
      if ((i >= nblocks && !revcomp) || (i >= (nblocks / 2) && revcomp))
	  die("Too many alignment files for format\n");
      
      if (str_split(line, NULL, matches) > 1) {
	strncpy(msa_fname, ((String*)lst_get_ptr(matches, 0))->chars,
		STR_MED_LEN);
	strncpy(ih_fname, ((String*)lst_get_ptr(matches, 1))->chars,
		STR_MED_LEN);
      } else {
	strncpy(msa_fname, line->chars, STR_MED_LEN);
      }
      
      if (!quiet)
	fprintf(stderr, "\t%s (%d of %d)\n", msa_fname, (i+1), nblocks);
      msa_f = fopen_fname(msa_fname, "r");
      msa = msa_new_from_file(msa_f, format, (char*)(alphabet->chars));
      lst_push_ptr(msas, msa);
      fclose(msa_f);

      str_cpy_charstr(fname, msa_fname);
      str_remove_path(fname);
      str_root(fname, '.');
      lst_push_ptr(dmpmsa->seqnames, str_dup(fname));

      if (msa->length > dmpmsa->max_seqlen)
	dmpmsa->max_seqlen = msa->length;

      /* Need to process the input alignments individually to be sure all have
	 SS structures in place */
      if (msa_alph_has_lowercase(msa))
	msa_toupper(msa);
      msa_remove_N_from_alph(msa);
      
      if (msa->ss == NULL) {
	if (!quiet)
	  fprintf(stderr, "\t\tExtracting sufficient statistics...\n");
	ss_from_msas(msa, 1, TRUE, NULL, NULL, NULL, -1);
      }
      else if (msa->ss->tuple_idx == NULL)
	die("ERROR: ordered representation of alignment required.\n");

      if (i == 0) {
	ncats = msa->ncats;
	tuple_size = msa->ss->tuple_size;
      } else if (ncats != msa->ncats || tuple_size != msa->ss->tuple_size) {   
       	die("ERROR: All MSA's must have same ncats and tuple_size.\n");
      }
      
      if (do_ih) {
	if (strlen(ih_fname) == 0)
	  die("ERROR: --indel-model requires an indel history for all alignment files! Try 'dmsample -h'\n");
	if (!quiet)
	  fprintf(stderr, "\tReading indel history for %s from %s\n", 
		  msa_fname, ih_fname);
	ih_f = fopen_fname(ih_fname, "r");
	dmpmsa->ih[i] = ih_new_from_file(ih_f);
	fclose(ih_f);
      }
      
      /* Create reverse complement of sequence i, if called for */
      if (revcomp == TRUE) {
	if (!quiet)
	  fprintf(stderr, "\tProducing reverse complement of %s\n",
		  msa_fname);
	msa_rc = msa_create_copy(msa, (format == SS ? 1 : 0));
	msa_reverse_compl(msa_rc);
	lst_push_ptr(msas, msa_rc);
        str_append_charstr(fname, "_rc");
	lst_push_ptr(dmpmsa->seqnames, str_dup(fname));
	/* TO DO: There is a way to create an indel history for the reverse 
	   seq from the indel history for the forward seq -- this needs to be
           implemented before using indel mod! */
      }

      i++;
    }
  }
  if ((i < (nblocks - 1) && !revcomp) || (i < (nblocks / 2)-1 && revcomp))
    die("ERROR: Not enough files in alignments list!\n");

  /* Create the PooledMSA structure from the list of MSA's */
  if (!quiet)
    fprintf(stderr, "\tBuilding pooled SS...\n");
/*   fprintf(stderr, "blocks %p, ", blocks); */
  dmpmsa->pmsa = ss_pooled_from_msas(msas, tuple_size, ncats, NULL);
/*   fprintf(stderr, "blocks_alloc %p, blocks->source_msas %p, n %d\n", blocks, */
/* 	  blocks->source_msas, lst_size(blocks->source_msas)); */
/*   fprintf(stderr, "msas[i] %p, blocks->source_msas[0] %p\n", */
/* 	  lst_get_ptr(msas, 0), lst_get_ptr(blocks->source_msas, 0)); */

  str_re_free(blocks_re);
  str_re_free(alph_re);
  str_re_free(format_re);
  lst_free(matches);
  str_free(alphabet);
  str_free(line);
  str_free(fname);
  return dmpmsa;
}

void dms_free_dmpmsa_struct(DMotifPmsaStruct *dmpmsa) {
  int i;
  List *msas;
  msas = dmpmsa->pmsa->source_msas;

  for (i = 0; i < lst_size(msas); i++) {
    msa_free(lst_get_ptr(msas, i));
    str_free(lst_get_ptr(dmpmsa->seqnames, i));
  }

  if(dmpmsa->ih != NULL) {
    for (i = 0; i < lst_size(msas); i++)
      ih_free(dmpmsa->ih[i]);
    free(dmpmsa->ih);
  }

  lst_free(dmpmsa->seqnames);
  ss_free_pooled_msa(dmpmsa->pmsa);
  free(dmpmsa);
  lst_free(msas);
}

/* Pared down version of tL_compute_log_likelihood with improved memory
   performance */
double dm_compute_log_likelihood(TreeModel *mod, MSA *msa, 
                                 double *col_scores, int cat) {

  int i, j, k, nodeidx, tupleidx, defined, nstates, alph_size,
    *partial_match, skip_fels, thisseq, ninform, observed_state;
  double retval, total_prob, **inside_joint, *tuple_scores, totl, totr;
  TreeNode *n;
  List *traversal;
  MarkovMatrix *lsubst_mat, *rsubst_mat;

  inside_joint = NULL;
  retval = 0;
  nstates = mod->rate_matrix->size;
  alph_size = strlen(mod->rate_matrix->states);

  /* allocate memory */
  partial_match = (int*)smalloc(alph_size * sizeof(int));
  inside_joint = (double**)smalloc(nstates * sizeof(double*));
  for (j = 0; j < nstates; j++) 
    inside_joint[j] = (double*)smalloc((mod->tree->nnodes+1)
				       * sizeof(double));
  
  if (col_scores != NULL) {
    tuple_scores = (double*)smalloc(msa->ss->ntuples 
				    * sizeof(double));
    for (tupleidx = 0; tupleidx < msa->ss->ntuples; tupleidx++)
      tuple_scores[tupleidx] = 0;
  }

  /* set up leaf to sequence mapping, if necessary */
  if (mod->msa_seq_idx == NULL)
    tm_build_seq_idx(mod, msa);

  /* set up prob matrices, if any are undefined */
  for (i = 0, defined = TRUE; defined && i < mod->tree->nnodes; i++) {
    if (((TreeNode*)lst_get_ptr(mod->tree->nodes, i))->parent == NULL)
      continue;  		/* skip root */
    for (j = 0; j < mod->nratecats; j++)
      if (mod->P[i][j] == NULL) defined = FALSE;
  }
  if (!defined) {
    tm_set_subst_matrices(mod);
  }

  /* Compute log probability for each tuple in the dataset */
  for (tupleidx = 0; tupleidx < msa->ss->ntuples; tupleidx++) {
    skip_fels = FALSE;
    total_prob = 0;

    /* check for gaps and whether column is informative, if necessary */
    if (!mod->allow_gaps)
      for (j = 0; !skip_fels && j < msa->nseqs; j++) 
        if (ss_get_char_tuple(msa, tupleidx, j, 0) == GAP_CHAR) 
          skip_fels = TRUE;
    if (!skip_fels && mod->inform_reqd) {
      ninform = 0;
      for (j = 0; j < msa->nseqs; j++) {
        if (msa->is_informative != NULL && !msa->is_informative[j])
          continue;
        else if (!msa->is_missing[(int)ss_get_char_tuple(msa, tupleidx, j, 0)])
          ninform++;
      }
      if (ninform < 2) skip_fels = TRUE;
    }
          
    if (!skip_fels) {
      traversal = tr_postorder(mod->tree);
      for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
	n = lst_get_ptr(traversal, nodeidx);      
	if (n->lchild == NULL) { 
	  /* leaf: base case of recursion */	  
	  assert(n->name != NULL);
	  thisseq = mod->msa_seq_idx[n->id];
	  observed_state = mod->rate_matrix->
	    inv_states[(int)ss_get_char_tuple(msa, tupleidx, 
					      thisseq, 0)];
	  
	  for (i = 0; i < alph_size; i++) {
	    if (observed_state < 0 || i == observed_state) 
	      partial_match[i] = 1;
	    else
	      partial_match[i] = 0; 
	  }
	  
	  /* now find the intersection of the partial matches */
	  for (i = 0; i < nstates; i++)
	    inside_joint[i][n->id] = (double)partial_match[i]; 
	}
	else {                    
	  /* general recursive case */
	  lsubst_mat = mod->P[n->lchild->id][0];
	  rsubst_mat = mod->P[n->rchild->id][0];
	  for (i = 0; i < nstates; i++) {
	    totl = totr = 0;
	    for (j = 0; j < nstates; j++) 
	      totl += inside_joint[j][n->lchild->id] *
		mm_get(lsubst_mat, i, j);
	    
	    for (k = 0; k < nstates; k++) 
	      totr += inside_joint[k][n->rchild->id] *
		mm_get(rsubst_mat, i, k);
	    
	    inside_joint[i][n->id] = totl * totr;
	  }
	}
      }

      for (i = 0; i < nstates; i++) {
	total_prob += vec_get(mod->backgd_freqs, i) * 
	  inside_joint[i][mod->tree->id] * mod->freqK[0];
      }
    } /* if skip_fels */
    
    total_prob = log2(total_prob);
    tuple_scores[tupleidx] = total_prob; /* Unweighted log prob! */
    total_prob *= msa->ss->counts[tupleidx]; /* log space */
    retval += total_prob;     /* log space */        
  } /* for tupleidx */
  
  for (j = 0; j < nstates; j++)
    free(inside_joint[j]);
  free(inside_joint);
  if (col_scores != NULL) {
    for (i = 0; i < msa->length; i++)
      col_scores[i] = tuple_scores[msa->ss->tuple_idx[i]];
    free(tuple_scores);
  }
  free(partial_match);
/*   free(mod->msa_seq_idx); */
/*   mod->msa_seq_idx = NULL; */
/*   dm_free_subst_matrices(mod); */
/*   free(mod->in_subtree); */
/*   mod->in_subtree = NULL; */
/*   lst_free(traversal); */
/*   mod->tree->postorder = NULL; */
  return(retval);
}

void dm_set_subst_mods(DMotifPhyloHmm *dm) {
  int i;
  TreeModel *mod;

  for (i = 0; i < dm->phmm->hmm->nstates; i++) {
    mod = dm->phmm->mods[ (dm->phmm->state_to_mod[i]) ];
    tm_set_subst_matrices(mod);
  }
}

void dm_free_subst_matrices(TreeModel *tm) {
  int i, j;
  TreeNode *n;
  for (i = 0; i < tm->tree->nnodes; i++) {
    n = lst_get_ptr(tm->tree->nodes, i);
    if (n->parent == NULL) continue;
    for (j = 0; j < tm->nratecats; j++) {
      mm_free(tm->P[i][j]);
      tm->P[i][j] = NULL;
    }
    free(tm->P[i]);
    tm->P[i] = NULL;
  }
  free(tm->P);
  tm->P = NULL;
}

/* Given an MSA without indels and a dmotif indel model, produce an MSA with
   indels and a corresponding indel history according to the dmotif indel
   model. */
MSA *dm_indel_mask(DMotifPhyloHmm *dm, MSA *msa, IndelHistory *ih,
		   int *path) {
  
  int i, j, k, state, c;
  unsigned long int seed;
  FILE *devrandom;
  double trans;
  TreeNode *n, *ns;
  List *outside;
  DMotifIndelModel *im;
  MSA *indel_msa;

  /* Seed the random number generator based on the value of /dev/random */
  devrandom = fopen("/dev/random", "r");
  fread(&seed, sizeof(seed), 1, devrandom);
  srandom(abs(seed));
  fclose(devrandom);

  outside = lst_new_ptr(dm->phmm->mods[0]->tree->nnodes);
  
  /* Do a preorder traversal to overlay insertions first */
  for (i = 0; i < ih->tree->nnodes; i++) {
    n = lst_get_ptr(ih->tree->nodes, i);
    if (n == ih->tree)
      continue;
    
    for (j = 0; j < msa->length; j++) {
      state = path[j];
      im = dm->indel_mods[state];

/*       fprintf(stderr, "i %d, j %d, state %d, n->id %d, name %s, ", */
/* 	      i, j, state, n->id, n->name); */

      if (ih->indel_strings[n->id][j] == INS)
	continue;

      if (j == 0) { /* Use begin probs if first seq. column */
	trans = vec_get(im->branch_mods[n->id]->beg_probs, CHILDINS);
      } else {
	c = ih->indel_strings[n->id][j-1];
	if (c == BASE) {
	  trans = mm_get(im->branch_mods[n->id]->probs, MATCH, CHILDINS);
	} else if (c == DEL) {
	  trans = mm_get(im->branch_mods[n->id]->probs, CHILDDEL, CHILDINS);
	} else /* c == INS */ {
	  trans = mm_get(im->branch_mods[n->id]->probs, CHILDINS, CHILDINS);
	}
      }
      
      if (bn_draw(1, trans)) {
	ih->indel_strings[n->id][j] = INS;
	/* Insertions always propagate through the supertree above a node */
	tr_partition_nodes(ih->tree, n, NULL, outside);
	for (k = 0; k < lst_size(outside); k++) {
	  ns = lst_get_ptr(outside, k);
	  ih->indel_strings[ns->id][j] = INS;
	}
      }
    }
  }
      
  /* Do a preorder traversal to overlay deletions. */
  for (i = 0; i < ih->tree->nnodes; i++) {
    n = lst_get_ptr(ih->tree->nodes, i);
      
    if (n == ih->tree) /* skip root node */
      continue;

    for (j = 0; j < msa->length; j++) {
      state = path[j];
      im = dm->indel_mods[state];
      
      if (ih->indel_strings[n->id][j] == INS)
	continue;
      else if (ih->indel_strings[n->parent->id][j] == DEL)
	ih->indel_strings[n->id][j] = DEL; /* Deletion characters can only give
					      rise to deletion characters */
      else {
	if (j == 0) { /* Use begin probs if first seq. column */
	  trans = vec_get(im->branch_mods[n->id]->beg_probs, CHILDDEL);
	} else {
	  c = ih->indel_strings[n->id][j-1];
	  
	  if (c == BASE)
	    trans = mm_get(im->branch_mods[n->id]->probs, MATCH, CHILDDEL);
	  else if (c == DEL)
	    trans = mm_get(im->branch_mods[n->id]->probs, CHILDDEL, CHILDDEL);
	  else /* c == INS */
	    trans = mm_get(im->branch_mods[n->id]->probs, CHILDINS, CHILDDEL);
	}
	
	if (bn_draw(1, trans))
	  ih->indel_strings[n->id][j] = DEL;
      }
    }
  }
    
/*   for (i = 0; i < ih->tree->nnodes; i++) { */
    
/*     n = lst_get_ptr(ih->tree->nodes, i); */
/*     if (n == ih->tree) */
/*       continue; */
    
/*     fprintf(stderr, "i %d, n->name %s\n", i, n->name); */
/*     for (j = 0; j < msa->length; j++) { */
/*       fprintf(stderr, "%d", ih->indel_strings[n->id][j]); */
/*     } */
/*     fprintf(stderr, "\n"); */
/*   } */
  
  indel_msa = dmih_as_alignment(ih, msa);
  return indel_msa;
}

/* Incrementally cache the motifs hash to disk and reset the working hash */
Hashtable* dms_cache_hash(Hashtable *path_counts, char *cache_fname, 
			  int nstates, int nsamples, int reinit) {
  Hashtable *tmp_hash;
  FILE *hash_f;

  hash_f = fopen_fname(cache_fname, "w");
  dms_write_hash(path_counts, hash_f, nstates, nsamples);
  fclose(hash_f);
  hsh_free_with_vals(path_counts);
  if (reinit == -1) /* Do not reinitialize -- return null pointer */
    return NULL;
  tmp_hash = hsh_new(reinit);
  return tmp_hash;
}

Hashtable* dms_uncache(List *cache_files, int init_size, int nstates, 
		       int *nsamples, int quiet) {
  int i, csamples;
  char fname[STR_MED_LEN];
  Hashtable *path_counts, *tmp_hash;
  FILE *hash_f;
  path_counts = hsh_new(init_size);
  
  *nsamples = csamples = 0;

  for (i = 0; i < lst_size(cache_files); i++) {
    snprintf(fname, STR_MED_LEN, "%s", ((String*)lst_get_ptr(cache_files, i))->chars);
    
    if (!quiet)
      fprintf(stderr, "Reading sampling data from file %s...\n", fname);
    
    hash_f = fopen_fname(fname, "r");
    tmp_hash = dms_read_hash(hash_f, nstates, &csamples);
    *nsamples += csamples;
    fclose(hash_f);
    dms_combine_hashes(path_counts, tmp_hash, nstates);

    /* For debugging */
/*     sprintf(fname, "test_%i.1.hash", i); */
/*     hash_f = fopen_fname(fname, "w"); */
/*     dms_write_hash(tmp_hash, hash_f, nstates, csamples); */
/*     fclose(hash_f); */
/*     sprintf(fname, "test_%i.2.hash", i); */
/*     hash_f = fopen_fname(fname, "w"); */
/*     dms_write_hash(path_counts, hash_f, nstates, *nsamples); */
/*     fclose(hash_f); */
    /* End debug */

    hsh_free_with_vals(tmp_hash);
  }
  return path_counts;
}

/* Free memory associated with a DMotifPhyloHmm object */
void dm_free(DMotifPhyloHmm *dm) {
  int i;
  int nnodes = dm->phmm->mods[0]->tree->nnodes;
  int nstates = dm->phmm->hmm->nstates;     /* number of states */

  free(dm->state_to_branch);
  free(dm->state_to_event);
  free(dm->state_to_motifpos);
  free(dm->state_to_motif);
  /* free dm->branch_to_states */
  for (i = 0; i < nnodes; i++)
    lst_free(dm->branch_to_states[i]);
  free(dm->branch_to_states);
  /* free dm->indel_mods, if used */
  if (dm->indel_mods != NULL) {
    for (i = 0; i < nstates; i++)
      dmih_free(dm->indel_mods[i]);
    free(dm->indel_mods);
  }
  phmm_free(dm->phmm);
  mot_free(dm->m);
  free(dm);
}

List* dms_read_tmp_from_file(FILE *tmp_lst_f) {
  Regex *nfiles_re = str_re_new("#[[:space:]]*NFILES[[:space:]]*=[[:space:]]*([0-9]+)");
  int i, nfiles;
  String *line = str_new(STR_MED_LEN);
  List *matches = lst_new_ptr(2), *ret;

  i = 0;
  while (str_readline(line, tmp_lst_f) != EOF) {
    lst_clear(matches);
    str_trim(line);
    if (line->length == 0) continue;
    
    if (str_re_match(line, nfiles_re, matches, 1) >= 0) {
      str_as_int(lst_get_ptr(matches, 1), &nfiles);
      ret = lst_new_ptr(nfiles);
    } else {
      lst_push_ptr(ret, (void*)str_dup(line));
      i++;
    }
  }
  if (i != nfiles) die("ERROR: Wrong number of lines in temp file list!\n");
  return ret;
}

/* Dump the sampling data to an output file for debugging purposes */
void dms_dump_sample_data(int sample, int thread_id, char *seqname, 
			  int seqlen, int *path, int **trans, 
			  Hashtable *path_counts, FILE *out, int dim) {
  int i;

  /* Print the sample and thread numbers */
  fprintf(out, "Sample %d, thread %d\n", sample, thread_id);

  /* Print the sequence name */
  fprintf(out, "%s\n\n", seqname);
  
  /* Print out the path */
  fprintf(out, "Path from hmm_stochastic_traceback:\n");
  for (i = 0; i < seqlen; i++)
    fprintf(out, "%d ", path[i]);
  fprintf(out, "\n\n");

  /* Print the transition counts */
  fprintf(out, "Cumulative transition counts for thread %d:\n", thread_id);
  for (i = 0; i < 4; i++) {
    if (i == 0)
      fprintf(out, "Mu trans: ");
    else if (i == 1)
      fprintf(out, "Nu trans: ");
    else if (i == 2)
      fprintf(out, "Phi trans: ");
    else
      fprintf(out, "Zeta trans: ");

    fprintf(out, "alpha: %d, beta: %d\n", trans[i][0], trans[i][1]);
  }

  /* Print out the hashtable */
  fprintf(out, "\nMotifs Hash:\n");
  dms_write_hash(path_counts, out, dim, 1);
}

/* Adjust coordinates of a motif GFF feature to the coordinate frame of the
   reference sequence. Assumes all features will be found within the coordinate
   frame of the given MSA -- behavior is undefined if this is not true! */
void dms_map_gff_coords(PooledMSA *blocks, int seqidx, GFF_Feature *f,
			int from_seq, int to_seq) {

  int i, j, s, e, offset, fseq, tseq;
  static msa_coord_map ***all_maps = NULL;
  msa_coord_map **maps = NULL, *from_map = NULL, *to_map = NULL;
  MSA *msa;
  String *prev_name = NULL;
  
  /* Check for static maps and allocate if needed */
  if (all_maps == NULL) {
    all_maps = (msa_coord_map***)smalloc(lst_size(blocks->source_msas) *
				     sizeof(msa_coord_map**));
    for (i = 0; i < lst_size(blocks->source_msas); i++) {
      msa = lst_get_ptr(blocks->source_msas, i);
      all_maps[i] = (msa_coord_map**)smalloc((msa->nseqs + 1) *
					     sizeof(msa_coord_map*));
      for (j = 0; j <= msa->nseqs; j++)
	all_maps[i][j] = NULL;
    }      
  }
  
  /* Get appropriate pointers for the sequence containing this feature */
  msa = lst_get_ptr(blocks->source_msas, seqidx);
  offset = msa->idx_offset;
  maps = all_maps[seqidx];
  fseq = from_seq;
  tseq = to_seq;

  if (from_seq == to_seq) {
    f->start += offset;
    f->end += offset;
  }
  else if (from_seq == -1) {
    if (str_equals_nocase_charstr(f->seqname, "MSA"))
      fseq = 0;
    else if (prev_name == NULL || !str_equals(prev_name, f->seqname)) {
      /* generally all seqs will have the same name; take advantage of
	 this property */
      if ((fseq = msa_get_seq_idx(msa, f->seqname->chars)) == -1)
	die("ERROR: name %s not present in MSA.\n", f->seqname->chars);
      fseq++;                 /* need 1-based index */
      prev_name = f->seqname;
    }
  }
  else if (to_seq == -1) {
    if (str_equals_nocase_charstr(f->seqname, "MSA"))
      tseq = 0;
    else if (prev_name == NULL || !str_equals(prev_name, f->seqname)) {
      if ((tseq = msa_get_seq_idx(msa, f->seqname->chars)) == -1)
	die("ERROR: name %s not present in MSA.\n", f->seqname->chars);
      prev_name = f->seqname;
      tseq++;                 /* need 1-based index */
    }
  }
  
  if ((from_map = maps[fseq]) == NULL && fseq > 0) 
    from_map = maps[fseq] = msa_build_coord_map(msa, fseq);
  
  if ((to_map = maps[tseq]) == NULL && tseq > 0) 
    to_map = maps[tseq] = msa_build_coord_map(msa, tseq);
  
/*   fprintf(stderr, "f->start %d, f->end %d, to_map->msa_len %d, to_map->seq_len %d, offset %d, ", */
/* 	  f->start, f->end, to_map->msa_len, to_map->seq_len, offset); */
  
  /* from_map, to_map will be NULL iff fseq, to_seq are 0 */
  s = msa_map_seq_to_seq(from_map, to_map, f->start);
  e = msa_map_seq_to_seq(from_map, to_map, f->end);
    
  f->start = (s < 0 ? 1 : s) + offset;

  if (e < 0)
    f->end = (to_map != NULL ? to_map->seq_len : msa->length) + offset;
  else
    f->end = e + offset;
  
/*   fprintf(stderr, "f->start new %d, f->end new %d\n", f->start, f->end); */
}
