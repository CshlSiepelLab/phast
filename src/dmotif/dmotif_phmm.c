 
/* dmotif phylo-HMM */
#include <dmotif_phmm.h>
#include <stacks.h>
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
                       double nu, double phi, double zeta, double xi,
		       int xi_mode, double alpha_c, double beta_c, 
		       double tau_c, double epsilon_c, double alpha_n, 
		       double beta_n, double tau_n, double epsilon_n,
		       int estim_gamma, int estim_omega, int estim_phi, 
		       int estim_zeta, int estim_xi, subst_mod_type mmod_type,
		       int scale_by_branch) {

  int i, state, pos, e, b, nleaves, k, nstates;
  HMM *hmm;
  TreeModel **models;
  TreeNode *n;
  CategoryMap *cm;
/*   CategoryMap *cm_copy; */
  List *inside, *outside;
  double *alpha, *beta, *tau, *epsilon;
  DMotifPhyloHmm *dm;

  /* Some sanity checks */
  assert(rho > 0 && rho < 1);
  assert(mu > 0 && 2*mu < 1);
  assert(nu > 0 && nu < 1);
  assert(zeta > 0 && zeta < 1);
  assert(xi > 0 && xi < 1);

  /* Compute the number of leaves, branches and states */
  nleaves = (source_mod->tree->nnodes + 1) / 2;
  k = 2*nleaves - 3;                  /* no. branches (unrooted tree) */
  nstates = (2 + 2*k) * (m->width+1); /* number of states */

  /* Allocate memory for temporary items */
  inside = lst_new_ptr(source_mod->tree->nnodes);
  outside = lst_new_ptr(source_mod->tree->nnodes);
  alpha = (double*)smalloc(source_mod->tree->nnodes * sizeof(double));
  beta = (double*)smalloc(source_mod->tree->nnodes * sizeof(double));
  tau = (double*)smalloc(source_mod->tree->nnodes * sizeof(double));
  epsilon = (double*)smalloc(source_mod->tree->nnodes * sizeof(double));

  /* Allocate memory that will be retained in the DMotifPhyloHmm object */
  models = (TreeModel**)smalloc(nstates * sizeof(TreeModel*));

  /* Set up the DMotifPhyloHmm object */
  dm = (DMotifPhyloHmm*)smalloc(sizeof(DMotifPhyloHmm));
  dm->state_to_branch = (int*)smalloc(nstates 
				      * sizeof(int));
  dm->state_to_event = (int*)smalloc(nstates 
				     * sizeof(int));
  dm->state_to_motifpos = (int*)smalloc(nstates 
					* sizeof(int));
  dm->state_to_motif = (int*)smalloc(nstates 
				     * sizeof(int));
  dm->branch_to_states = (List**)smalloc(source_mod->tree->nnodes *
					 sizeof(List*));
  dm->rho = rho;
  dm->mu = mu;
  dm->nu = nu;
  dm->phi = phi;
  dm->zeta = zeta;
  dm->xi = xi;
  dm->estim_gamma = estim_gamma;
  dm->estim_omega = estim_omega;
  dm->estim_phi = estim_phi;
  dm->estim_zeta = estim_zeta;
  dm->estim_xi = estim_xi;
  dm->m = m;
  dm->k = k;
  dm->indel_mods = NULL;

  /* set up mappings between states and branches.  Branches are
     indexed by ids of child nodes; "birth" means functional element
     arose on branch, "death" means element was lost on branch. */
  for (i = 0; i < source_mod->tree->nnodes; i++)
    dm->branch_to_states[i] = lst_new_int(2 * m->width + 2);
  
  state = 0;
  for (i = 0; i < source_mod->tree->nnodes; i++) {
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
      /* Neutral background state */
      models[state] = source_mod;
    else if (pos == -1) /* Background state */
      models[state] = tm_create_copy(source_mod);
    else { /* Motif state */

      if (mmod_type == F81 || dm->state_to_event[state] == NEUT) {
	/* F81 model with motif freqs as equilibrium freqs */
	models[state] = tm_new(tr_create_copy(source_mod->tree), NULL,
			       vec_create_copy(m->probs[pos]), F81, 
			       DEFAULT_ALPHABET, 1, -1, NULL, -1);
	//      fprintf(stderr, "state %d pos %d ", state, pos);
	//      vec_print(m->probs[pos], stderr);
	
      } else {
	/* Halpern-Bruno model -- use REV bg model as the source */
	models[state] = tm_create_copy(source_mod);
      }
    }
    
    if (dm->state_to_event[state] == NEUT)
      continue;
    
    n = lst_get_ptr(models[state]->tree->nodes, dm->state_to_branch[state]);
    tr_partition_nodes(models[state]->tree, n, inside, outside);

    /* adjust branch lengths if BIRTH, DEATH, or CONS and either background
       state or using F81 as the motif model type */
    if (pos == -1 || mmod_type == F81) {
      if (dm->state_to_event[state] == DEATH)      
	for (i = 0; i < lst_size(outside); i++)
	  ((TreeNode*)lst_get_ptr(outside, i))->dparent *= rho;
      else if (dm->state_to_event[state] == BIRTH ||
	       dm->state_to_event[state] == CONS)
	for (i = 0; i < lst_size(inside); i++)
	  ((TreeNode*)lst_get_ptr(inside, i))->dparent *= rho;
    }

    /* if BIRTH or DEATH, and mmod_type is F81, use alternative subst model 
       (the background model) for non-motif branches */
    if (pos != -1) {
      if (mmod_type == F81) {
	if (dm->state_to_event[state] == DEATH)
	  dm_set_backgd_branches(models[state], source_mod, inside);
	else if (dm->state_to_event[state] == BIRTH)
	  dm_set_backgd_branches(models[state], source_mod, outside);
	
      } else { /* Using Halpern-Bruno -- use alternative subst model for motif
		  branches. Conserved model is handled specially within
		  dm_set_motif_branches_hb and does not use the alt_subst_mods
		  attribute. */
/*       if (pos == -1) /\* BG states handled above *\/ */
/* 	continue; */
/* 	fprintf(stderr, "state %d, pos %d\n", state, pos); */
	if (dm->state_to_event[state] == BIRTH ||
	    dm->state_to_event[state] == CONS)
	  dm_set_motif_branches_hb(models[state], m->probs[pos], inside);
	else if (dm->state_to_event[state] == DEATH)
	  dm_set_motif_branches_hb(models[state], m->probs[pos], outside);  
      }
    }

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
  hmm = hmm_new_nstates(nstates, TRUE, TRUE);

  /* build category map */
  cm = dm_create_catmap(dm, source_mod, NULL);
  dm->phmm = phmm_new(hmm, models, cm, NULL, MISSING_DATA);
/*   cm_copy = dm->phmm->cm; */
/*   dm->phmm->cm = cm; */
/*   cm_free(cm_copy); */
  cm_free(cm); /* This is copied within phmm_new, so it's safe to free. */
  dm_set_transitions(dm, xi_mode, scale_by_branch);

  /* set up indel model, if necessary */
  if (alpha_c > 0) {
    dm->indel_mods = (DMotifIndelModel**)smalloc(nstates * 
						 sizeof(DMotifIndelModel*));
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
  
  /* Free all temporary memory items */
  lst_free(inside);
  lst_free(outside);
  free(alpha);
  free(beta);
  free(tau);
  free(epsilon);
  
  return dm;
}

void dm_set_transitions(DMotifPhyloHmm *dm, int xi_mode, int scale_by_branch) {
  int i, s1, s2, p1, p2, b1, b2, e1, e2;
  double trans, rowsum, treelen, *scale;
/*   double bsum; */
  HMM *hmm = dm->phmm->hmm;
  TreeNode *nn;

  if (scale_by_branch == TRUE) {
    /* Adjusting lineage-specific transistion probs by branch length instead 
       of uniformly distributing across all branches -- need total tree length
       to do this. */
    scale = smalloc(dm->phmm->mods[0]->tree->nnodes * sizeof(double));
    treelen = 0;
    for (i = 0; i < dm->phmm->mods[0]->tree->nnodes; i++) {
      nn = lst_get_ptr(dm->phmm->mods[0]->tree->nodes, i);
      scale[nn->id] = 0;
/*       if (nn->parent == dm->phmm->mods[0]->tree) */
/* 	continue; */
      treelen += nn->dparent;
    }
    
/*     bsum = 0; */
    for (i = 0; i < dm->phmm->mods[0]->tree->nnodes; i++) {
      nn = lst_get_ptr(dm->phmm->mods[0]->tree->nodes, i);
      /* Nodes immediately below the root need special handling, since we
	 cannot distinguish gain on one branch from loss on the other. We have
	 already accounted for this in setting up the HMM, by assigning gains
	 to one branch and losses to the opposite branch. Since the underlying
	 tree structure is unrooted, we need to count both branches under the
	 root toward the branch length for these gain/loss events. */
      if (nn->parent == dm->phmm->mods[0]->tree)
	scale[nn->id] = ((nn->parent->lchild->dparent +
			  nn->parent->rchild->dparent) / treelen) / 2;
      else
	scale[nn->id] = (nn->dparent / treelen) / 2;
/*       bsum += scale[nn->id]; */
    }
/*     fprintf(stderr, "treelen %f, bsum %f\n", treelen, bsum); */
  }
  
  /* set transition probs for each state */
  for (s1 = 0; s1 < hmm->nstates; s1++) {
    p1 = dm->state_to_motifpos[s1];
    b1 = dm->state_to_branch[s1];
    e1 = dm->state_to_event[s1];
    rowsum = 0;
/*     bsum = 0; */

    for (s2 = 0; s2 < hmm->nstates; s2++) {
      p2 = dm->state_to_motifpos[s2];
      b2 = dm->state_to_branch[s2];
      e2 = dm->state_to_event[s2];
      trans = 0;
      
      nn = lst_get_ptr(dm->phmm->mods[s2]->tree->nodes,
		       dm->state_to_branch[s2]);
      
      /* Set params for each type of legal transition */
      
      /* Transitions from neutral background state */
      if (e1 == NEUT && p1 == -1) {
	/* Into another neutral state (bg or motif) */
	if (e2 == NEUT) {
	  if (p2 == -1) { /* Self-transition */
	    if (xi_mode == FALSE)
	      trans = ((1 - dm->nu) * (1 - dm->zeta));
	    else
	      trans = ((1 - dm->nu) * (1 - dm->xi));
	  } else if (p2 == 0) { /* Motif position 1 */
	    if (xi_mode == FALSE)
	      trans = ((1 - dm->nu) * dm->zeta);
	    else
	      trans = ((1 - dm->nu) * dm->xi);
	  } else { /* Internal motif position */
	    continue;
	  }
	}
	
	/* Into a conserved state (bg or motif) */
	else if (e2 == CONS) {
	  if (p2 == -1) /* bg state */
	    trans = (dm->nu * (1 - dm->zeta) * (1 - dm->phi));
	  else if (p2 == 0) /* Motif position 1 */
	    trans = (dm->nu * dm->zeta * (1 - dm->phi));
	  else { /* Internal motif position */
	    continue;
	  }
	}

	/* Into a lineage-specific state (bg or motif) */
	else if (e2 == BIRTH || e2 == DEATH) {
	  if (p2 == -1) { /* bg state */
	    if (scale_by_branch == TRUE) {
	      trans = (dm->nu *  (1 - dm->zeta) * dm->phi) * scale[nn->id];
	    } else {
	      trans = (dm->nu *  (1 - dm->zeta) * dm->phi) / (2 * dm->k);
	    }
	  } else if (p2 == 0) { /* Motif position 1 */
            if (scale_by_branch == TRUE) {
	      trans = (dm->nu * dm->phi * dm->zeta) * scale[nn->id];
/* 	      bsum += scale[nn->id]; */
	    } else {
	      trans = (dm->nu * dm->phi * dm->zeta) / (2 * dm->k);
	    }
	  } else { /* Internal motif position */
	    continue;
	  }
	}
/* 	fprintf(stderr, "s1 %d, s2, %d, p1 %d, p2 %d, node %s, blen %f, treelen %f, scale %f, trans %.9f, rowsum %.9f, bsum %.9f\n", */
/* 		s1, s2, p1, p2, nn->name, n0->dparent, treelen, scale[nn->id], */
/* 		trans, (rowsum + trans), bsum); */	
      }

      /* Transitions between neutral motif states and other neutral motif 
	 states or neutral background state */
      else if (e1 == NEUT && p1 >= 0) {
	if (e2 != NEUT)
	  continue;
	if ((p1 < dm->m->width-1 && p2 == p1+1) ||
	    (p1 == dm->m->width-1 && p2 == -1))
	  trans = 1;
	else
	  continue;
      }
      
      /* Transitions from CONS, BIRTH and DEATH states */
      else if (e1 == CONS || e1 == BIRTH || e1 == DEATH) {
	
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
		   && p2 == -1) /* Last motif position or non-neutral
				   background to Nb */
	  trans = dm->mu;
	else continue;
      } else continue;
      
      rowsum += trans;
/*       fprintf(stderr, "s1 %d, s2, %d, p1 %d, p2 %d, node %s, blen %f, treelen %f, scale %f, trans %.9f, rowsum %.9f, bsum %.9f\n", */
/* 	      s1, s2, p1, p2, nn->name, n0->dparent, treelen, scale[nn->id], */
/* 	      trans, rowsum, bsum); */
      mm_set(hmm->transition_matrix, s1, s2, trans);
      //      printf("s1 %d, s2 %d, p1 %d, p2 %d, trans %f\n", s1, s2, p1, p2, trans);
    }
    
    /* check that sum of outgoing arcs is 1 */
    rowsum = 0;
    for (s2 = 0; s2 < hmm->nstates; s2++)
      rowsum += mm_get(hmm->transition_matrix, s1, s2);
/*     fprintf(stderr, "%d %f (%d, ", s1, rowsum, p1); */
/*     if (e1 == NEUT) */
/*       fprintf(stderr, "NEUT, "); */
/*     else if (e1 == CONS) */
/*       fprintf(stderr, "CONS, "); */
/*     else if (e1 == BIRTH) */
/*       fprintf(stderr, "BIRTH, "); */
/*     else */
/*       fprintf(stderr, "DEATH, "); */
    
/*     nn = lst_get_ptr(dm->phmm->mods[s1]->tree->nodes, dm->state_to_branch[s1]); */
/*     printf("%s)\n", nn->name); */
    
    assert(fabs(rowsum-1) < 1e-6);
  }

  /* Begin transitions */
  rowsum = 0;
  for (s1 = 0; s1 < hmm->nstates; s1++) {
    /* motif stationary probs */

    if (xi_mode == FALSE) {

      if (dm->state_to_motifpos[s1] == -1) /* background state */
	trans = 1 / (1 + dm->zeta);
      else if (dm->state_to_motifpos[s1] == 0) /* first motif position */
	trans = dm->zeta / (1 + dm->zeta);      
      else /* internal motif position -- don't allow starts here! */
	trans = 0;

    } else { /* using xi for transitions into neutral motifs */

      if (dm->state_to_motifpos[s1] == -1) {/* background state */
	if (dm->state_to_event == NEUT)
	  trans = 1 / (1 + dm->xi);
	else
	  trans = 1 / (1 + dm->zeta);
      } else if (dm->state_to_motifpos[s1] == 0) { /* first motif position */
	if (dm->state_to_event == NEUT)
	  trans = dm->xi / (1 + dm->xi);
	else
	  trans = dm->zeta / (1 + dm->zeta);	
      } else { /* internal motif position -- don't allow starts here! */
	trans = 0;
      }
    }

    /* multiply by dless stationary probs*/
    nn = lst_get_ptr(dm->phmm->mods[s1]->tree->nodes,
		     dm->state_to_branch[s1]);
    
    if (dm->state_to_event[s1] == NEUT) {
      trans *= dm->mu/(dm->mu + dm->nu);
      
    } else if (dm->state_to_event[s1] == CONS) {
      trans *= dm->nu * (1-dm->phi) / (dm->mu + dm->nu);
      
    } else {
      if (scale_by_branch == FALSE) {
	trans *= dm->nu * dm->phi / ((dm->mu + dm->nu) * 2 * dm->k);
      } else {
	trans *= dm->nu * dm->phi / ((dm->mu + dm->nu) / scale[nn->id]);
      }
    }

    vec_set(hmm->begin_transitions, s1, trans);
    rowsum += trans;
/*     fprintf(stderr, "s1 %d, trans %f, rowsum %f\n", s1, trans, rowsum); */
  }
  assert(fabs(rowsum-1) < 1e-6);

  /* End transitions -- we need these to prevent predictions of partial motifs
     at the end of a ChIP fragment */
  rowsum = 0;
  for (s1 = 0; s1 < hmm->nstates; s1++) {

    if (xi_mode == FALSE) {

      if (dm->state_to_motifpos[s1] == -1) /* background state */
	trans = 1 / (1 + dm->zeta);
      else if (dm->state_to_motifpos[s1] == dm->m->width-1) /* motif end pos */
	trans = dm->zeta / (1 + dm->zeta);
      else /* first or internal motif position -- don't allow ends from here */
	trans = 0;

    } else { /* using xi for transitions into neutral motifs */

      if (dm->state_to_motifpos[s1] == -1) {/* background state */
	if (dm->state_to_event == NEUT)
	  trans = 1 / (1 + dm->xi);
	else
	  trans = 1 / (1 + dm->zeta);
      
      } else if (dm->state_to_motifpos[s1] == dm->m->width-1) { /* motif end */
	if (dm->state_to_event == NEUT)
	  trans = dm->xi / (1 + dm->xi);
	else
	  trans = dm->zeta / (1 + dm->zeta);
	
      } else { /* internal motif position -- don't allow starts here! */
	trans = 0;
      }
    }
    
    /* multiply by dless stationary probs*/
    nn = lst_get_ptr(dm->phmm->mods[s1]->tree->nodes,
		     dm->state_to_branch[s1]);

    if (dm->state_to_event[s1] == NEUT) {
      trans *= dm->mu/(dm->mu + dm->nu);
    } else if (dm->state_to_event[s1] == CONS) {
      trans *= dm->nu * (1-dm->phi) / (dm->mu + dm->nu);
    } else {
      if (scale_by_branch == FALSE)
	trans *= dm->nu * dm->phi / ((dm->mu + dm->nu) * 2 * dm->k);
      else
	trans *= dm->nu * dm->phi / ((dm->mu + dm->nu) / scale[nn->id]);
    }

    vec_set(hmm->end_transitions, s1, trans);
    rowsum += trans;
/*     fprintf(stderr, "s1 %d, trans %f, rowsum %f\n", s1, trans, rowsum); */
  }
  assert(fabs(rowsum-1) < 1e-6);

  if (scale_by_branch == TRUE)
    free(scale);
  
  hmm_reset(hmm);

/*   mm_pretty_print(stderr, hmm->transition_matrix); */
}

/* adjust emission probabilities to prevent birth/death events from
   spanning regions where they would be supported only by missing data */
void dm_handle_missing_data(DMotifPhyloHmm *dm, MSA *msa) {
  List **uninform = (List**)smalloc(msa->ss->ntuples * sizeof(List*));
  TreeNode *n, *tree = dm->phmm->mods[0]->tree, *ns;
  typedef enum {MISSING, DATA, DATA_LEFT, DATA_RIGHT} node_type;
  node_type *mark = (node_type*)smalloc(tree->nnodes * sizeof(node_type));
  List *traversal;
  int i, j, l, state, mod;

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
      
      /* resolve ambiguity, if necessary. */
      if (mark[n->id] == DATA_LEFT || mark[n->id] == DATA_RIGHT) {
        mark[n->id] = mark[n->parent->id];
	/* If the parent ID is missing, then there is no information regarding
	   losses on the branch leading to the subtree with data. Zero loss
	   states on this branch. */
	if (mark[n->parent->id] == MISSING) {
	  if (mark[n->lchild->id] != MISSING)
	    ns = n->lchild;
	  else
	    ns = n->rchild;
	  for (l = 0; l < lst_size(dm->branch_to_states[ns->id]); l++) {
	    state = lst_get_int(dm->branch_to_states[ns->id], l);
/* 	    if (dm->state_to_event[state] == DEATH) */
	    lst_push_int(uninform[i], state);
	  }
	}
      }
      /* now check if missing data in subtree; if so, the tuple is
         uninformative wrt this branch */
      if (mark[n->id] == MISSING) {
        for (l = 0; l < lst_size(dm->branch_to_states[n->id]); l++)
          lst_push_int(uninform[i],
		       lst_get_int(dm->branch_to_states[n->id], l));
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
void unpack_params(Vector *params, DMotifPhyloHmm *dm, int xi_mode,
		   int scale_by_branch) {
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
  if (dm->estim_xi)
    dm->xi = vec_get(params, params_idx++);
  dm_set_transitions(dm, xi_mode, scale_by_branch);
}

double lnl_wrapper(Vector *params, void *data, int xi_mode, 
		   int scale_by_branch) {
  DMotifPhyloHmm *dm = data;
  unpack_params(params, data, xi_mode, scale_by_branch);
  return -hmm_forward(dm->phmm->hmm, dm->phmm->emissions,
                      dm->phmm->alloc_len, dm->phmm->forward);
}

/* estimate free parameters */
double dm_estimate_transitions(DMotifPhyloHmm *dm, MSA *msa, int xi_mode,
			       int scale_by_branch) {
  int i, nparams = 0;
  Vector *params, *lb, *ub;
  double retval;

  /* This method for parameter estimation does not work well and will not be
     supported. */
  die("ERROR: Parameter estimation is no longer supported in dmotif. Please use dmsample. (See dmsample -h)\n");

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
  if (dm->estim_xi) nparams++;

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
  if (dm->estim_xi) {
    vec_set(params, nparams, dm->xi);
    vec_set(ub, nparams++, 0.1);
  }

/*   opt_bfgs(lnl_wrapper, params, dm, &retval, lb, ub, stderr, NULL, */
/*            OPT_HIGH_PREC, NULL); */

  unpack_params(params, dm, xi_mode, scale_by_branch);

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
					      * sizeof(AltSubstMod*));
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
	str_free(type);
	type = str_new_charstr(BACKGD_CAT_NAME);
      	retval->ranges[0] =
	  cm_new_category_range(type, 0, 0);
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
    beg = i;                /* begin of feature */
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
		       int quiet, char *cache_fname, int cache_int,
		       int xi_mode, int scale_by_branch) {

  int i, j, k, l, *path, **trans, **ref_paths = NULL, nwins;
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
  trans = (int**)smalloc((xi_mode == FALSE ? 4 : 5) * sizeof(int*));
  for (i = 0; i < 5; i++) {
    if (xi_mode == FALSE && i == 4)
      continue;
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
    if (xi_mode == TRUE)
      dm->xi = beta_draw((double)trans[4][0], (double)trans[4][1]);

    dm_set_transitions(dm, xi_mode, scale_by_branch);

    /* Reset transition counts to prior values */
    for (j = 0; j < 5; j++) {
      if (xi_mode == FALSE && j == 4)
	continue;
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
			      force_priors, xi_mode);
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
      dms_write_log(log, dm, trans, i, llh, query_gff, reference, nwins, 
		    xi_mode);

    /* Cache and reinitialize the hash every 100 samples after burn-in */
    if (i >= bsamples) {
      if (l == cache_int || i == nsamples+bsamples-1) {
	snprintf(cache_out, STR_MED_LEN, "%s.%d.tmp", cache_fname, k++);
	lst_push_ptr(cache_files, (void*)str_new_charstr(cache_out));
	/*      fprintf(stderr, "i %d, k %d, l %d, cache_file %s\n", i, k-1, l, cache_out); */
	dms_cache_hash(path_counts, cache_out, (2*dm->k)+2, l);
	hsh_clear_with_vals(path_counts);
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
			    List *seqnames, int max_seqlen, int bsamples,
			    int nsamples, int sample_interval, int **priors,
			    FILE *log, GFF_Set *reference, int ref_as_prior, 
			    int force_priors, int quiet, char *cache_fname,
			    int cache_int, ThreadPool *pool, int nthreads,
			    List **zeroed_states, int xi_mode, 
			    int scale_by_branch) {
  
  int i, j, k, l, t, **trans, ***thread_trans, nwins, threads, **thread_path;
  unsigned long int seed;
  FILE *devrandom;
/*   int **ref_paths = NULL; */
  double llh, *thread_llh, ***thread_emissions, ***thread_forward;
  char *cache_out;
  GFF_Set *query_gff = NULL, **thread_gff = NULL;
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
  lst_clear(cache_files);
  
  /* Simplifies looping for case where thread pool size == 0 */
  if (nthreads == 0)
    threads = 1;
  else
    threads = nthreads;

  /* Allocate the global hashtable */
  path_counts = hsh_new(max(100*lst_size(blocks->source_msas), 1000));
  
  /* Allocate the thread_specific hashtables */
  thread_counts = (Hashtable**)smalloc(threads * sizeof(Hashtable*));
  for (i = 0; i < threads; i++) {
    thread_counts[i] = hsh_new(max(100*lst_size(blocks->source_msas), 1000) /
			       threads);
  }
  
  /* Allocate space for the global transition counts and set to prior values */
  trans = (int**)smalloc((xi_mode == FALSE ? 4 : 5) * sizeof(int*));
  for (i = 0; i < 5; i++) {
    if (xi_mode == FALSE && i == 4)
      continue;
    trans[i] = (int*)smalloc(2 * sizeof(int));
    trans[i][0] = priors[i][0];
    trans[i][1] = priors[i][1];
  }

  /* Allocate space for the single-thread transition counts. Initialize these
     to 0 as pseudocounts are added into the global struct! */
  thread_trans = (int***)smalloc(threads * sizeof(int**));
  for (i = 0; i < threads; i++) {
    thread_trans[i] = (int**)smalloc((xi_mode == FALSE ? 4 : 5) * 
				     sizeof(int*));
    for (j = 0; j < 5; j++) {
      if (xi_mode == FALSE && j == 4)
	continue;
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

  /* Allocate space for thread-specific emissions and forward matrices, and
     paths, sized for the longest sequence in the dataset. */
  thread_emissions = (double***)smalloc(threads * sizeof(double**));
  thread_forward = (double***)smalloc(threads * sizeof(double**));
  thread_path = (int**)smalloc(threads * sizeof(int*));
  for (i = 0; i < threads; i++) {
    thread_emissions[i] = (double**)smalloc(dm->phmm->hmm->nstates *
					    sizeof(double*));
    thread_forward[i] = (double**)smalloc(dm->phmm->hmm->nstates *
					  sizeof(double*));
    thread_path[i] = (int*)smalloc(max_seqlen * sizeof(int));
    for (j = 0; j < dm->phmm->hmm->nstates; j++) {
      thread_emissions[i][j] = (double*)smalloc(max_seqlen *
					      sizeof(double));
      thread_forward[i][j] = (double*)smalloc(max_seqlen *
					    sizeof(double));
    }
  }
    
  /* If comparing to a reference GFF for logging purposes, allocate the
     stats array and compute the number of motif windows in the dataset. */
  nwins = 0;
  if (reference != NULL && log != NULL) {
    for (i = 0; i < lst_size(blocks->source_msas); i++)
      nwins += (blocks->lens[i] - (dm->m->width -1));
  
    /* Also set up GFF Sets for tracking of true and false positives */
    query_gff = gff_new_set();
    thread_gff = (GFF_Set**)smalloc(threads * sizeof(GFF_Set*));
    for (i = 0; i < threads; i++)
      thread_gff[i] = gff_new_set();
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
    if (query_gff != NULL) {
	gff_clear_set(query_gff);
	for (j = 0; j <  threads; j++)
	  gff_clear_set(thread_gff[j]);
    }
    
    /* Draw a set of params */
    dm->mu = beta_draw((double)trans[0][0], (double)trans[0][1]);
    dm->nu = beta_draw((double)trans[1][0], (double)trans[1][1]);
    dm->phi = beta_draw((double)trans[2][0], (double)trans[2][1]);
    dm->zeta = beta_draw((double)trans[3][0], (double)trans[3][1]);
    if (xi_mode == TRUE)
      dm->xi = beta_draw((double)trans[4][0], (double)trans[4][1]);
    
    dm_set_transitions(dm, xi_mode, scale_by_branch);

    /* Set transition_score matrix in the hmm */
    hmm_set_transition_score_matrix(dm->phmm->hmm);

    /* Reset transition counts to prior values */
    for (j = 0; j < 5; j++) {
      if (xi_mode == FALSE && j == 4)
	continue;
      trans[j][0] = priors[j][0];
      trans[j][1] = priors[j][1];
    }

    /* Set up the work list for multithreading over sequences. First make sure
       the existing list is empty. Once this is set up, we can reuse the
       work items in subsequent samples by just changing two values. */
    if (i == 0) {
      if (!lst_empty(work)) {
	for (j = 0; j < lst_size(blocks->source_msas); j++) {
	  data = lst_get_ptr(work, j);
	  free(data);
	}
	lst_clear(work);
      }
      /* Set up the individual work items */
      for (j = 0; j < lst_size(blocks->source_msas); j++) {
	data = (DMsamplingThreadData*)smalloc(sizeof(DMsamplingThreadData));
	
	data->dm = dm;
	data->blocks = blocks;
	data->thread_path = thread_path;
	data->thread_emissions = thread_emissions;
	data->thread_forward = thread_forward;
	data->ih = (ih == NULL ? NULL : ih[j]);
	data->zeroed_states = (zeroed_states == NULL ? NULL : zeroed_states[j]);
	data->tuple_scores = tuple_scores;
	
	data->thread_llh = thread_llh;
	data->thread_trans = thread_trans;
	data->thread_counts = thread_counts;
	
	data->do_sample = (i < bsamples ? 0 : (i % sample_interval ? 0 : 1));
	data->seqnum = j;
	data->seqname = ((char*)((String*)lst_get_ptr(seqnames, j))->chars);
	data->log = log;
	data->thread_gff = thread_gff;
	data->do_reference = (reference == NULL ? 0 : (log == NULL ? 0 : 1));
	data->nthreads = nthreads;
	data->p = pool;
	data->sample = i;
	data->xi_mode = xi_mode;
	
	/*       fprintf(stderr, "i %d, j %d, t %d, data->llh %p, data->trans %p, data->path_counts %p, data->do_sample %d\n", */
	/* 	      i, j, t, data->llh, data->trans, data->path_counts, */
	/* 	      data->do_sample); */
	
	lst_push(work, (void*)(&data));
      }
    } else {
      /* Only a couple things need to be changed in the work list from sample
	 to sample -- only changing things that need to be changed instead of
	 freeing and reallocating/rebuilding the work list on each sample
	 should save a lot of memory overhead. -- 3 June 2009 */
      for (j = 0; j < lst_size(blocks->source_msas); j++) {
	data = (DMsamplingThreadData*)lst_get_ptr(work, j);
	data->do_sample = (i < bsamples ? 0 : (i % sample_interval ? 0 : 1));
	data->sample = i;
/* 	fprintf(stderr, "i %d, j %d, seqnum %d, seqname %s, do_sample %d, sample %d\n", */
/* 		i, j, data->seqnum, data->seqname, data->do_sample, */
/* 		data->sample); */
      }
    }

    thr_foreach(pool, work, dms_launch_sample_thread);

    /* Combine transition counts from all threads into global counts and
       reset thread-specific counters. Also total up global llh and reset
       thread-specific vals to 0. */
    llh = 0;
    for (t = 0; t < threads; t++) {
      for (j = 0; j < 5; j++) {
	if (xi_mode == FALSE && j == 4)
	  continue;
	trans[j][0] += thread_trans[t][j][0];
	trans[j][1] += thread_trans[t][j][1];
	thread_trans[t][j][0] = thread_trans[t][j][1] = 0;
      }
      llh += thread_llh[t];
      thread_llh[t] = 0;
    }
    
    if (log != NULL) { /* Write stats for this sample to the log file */
      if (query_gff != NULL) {
	for (t = 0; t < threads; t++)
	  dms_merge_thread_gffs(query_gff, thread_gff[t]);
      }
      dms_write_log(log, dm, trans, i, llh, query_gff, reference, nwins, 
		    xi_mode);    
    }

    /* Cache and reinitialize the hash periodically after burn-in */
    if (i >= bsamples) {

      if (l == cache_int || i == nsamples+bsamples-1) {
	/* Fold together hashes from individual threads */
	for (t = 0; t < threads; t++) {
	  dms_combine_hashes(path_counts, thread_counts[t], (2*dm->k)+2);
	  hsh_clear_with_vals(thread_counts[t]);
	}
	/* build and store a temp file name to hold the cached data */
	snprintf(cache_out, STR_MED_LEN, "%s.%d.tmp", cache_fname, k++);
	lst_push_ptr(cache_files, (void*)str_new_charstr(cache_out));
/* 	fprintf(stderr, "i %d, k %d, l %d, cache_file %s\n", i, k-1, l, cache_out); */
	dms_cache_hash(path_counts, cache_out, (2*dm->k)+2, l);
	hsh_clear_with_vals(path_counts);
	l = 1;
      } else {
	l++;
      }
    }
  }

  fprintf(stderr, "\n\tSampling complete.\n");

  /* Clean up our area */
  if (query_gff != NULL) {
    gff_free_set(query_gff);
    for (i = 0; i < threads; i++)
      gff_free_set(thread_gff[i]);
  }
  free(thread_gff);

  for (i = 0; i < lst_size(blocks->source_msas); i++) {
    data = lst_get_ptr(work, i);
    free(data);
  }
  lst_free(work);

  for (i = 0; i < 5; i++) {
    if (xi_mode == FALSE && i == 4)
      continue;
    free(trans[i]);
  }
  free(trans);

  for (i = 0; i < threads; i++) {
    for (j = 0; j < 5; j++) {
      if (xi_mode == FALSE && j == 4)
	continue;
      free(thread_trans[i][j]);
    }
    free(thread_trans[i]);
    
    for (j = 0; j < dm->phmm->hmm->nstates; j++) {
/*       fprintf(stderr, "i %d, j %d\n", i, j); */
      free(thread_emissions[i][j]);
      free(thread_forward[i][j]);
    }
    free(thread_emissions[i]);
    free(thread_forward[i]);
    hsh_free_with_vals(thread_counts[i]);
    free(thread_path[i]);
  }
  hsh_free_with_vals(path_counts);
  free(thread_trans);
  free(thread_counts);
  free(thread_llh);
  free(thread_path);
  free(thread_emissions);
  free(thread_forward);
  free(cache_out);
  
  return cache_files;
}

/* This will be used to run individual sampling threads, called by
   dms_launch_sample_thread. */
void dms_sample_path(DMotifPhyloHmm *dm, PooledMSA *blocks, IndelHistory *ih, 
		     double **tuple_scores, double *thread_llh,
		     int **thread_path,
		     double ***thread_emissions, double ***thread_forward,
		     int ***thread_trans, Hashtable **thread_counts, 
		     int do_sample, int seqnum, char *seqname, FILE *log, 
		     GFF_Set **thread_gff, int do_reference, int nthreads,
		     ThreadPool *pool, int sample, List *zeroed_states,
		     int xi_mode) {

  int *path, thread_idx, **trans;
  double **emissions, **forward_scores, *llh;
  MSA *msa;
  Hashtable *path_counts;
  GFF_Set *query_gff;

  /* Get a pointer to the msa we're looking at */
  msa = lst_get_ptr(blocks->source_msas, seqnum);

  /* Get the index of the current thread and appropriate vectors in thread-
     specific structures. */
  thread_idx = thr_index(pool);
  
  if (thread_idx == -1)
    die("ERROR: Attempt to access out-of-bounds thread in dms_sample_path!\n");

  trans = thread_trans[thread_idx];
  llh = &(thread_llh[thread_idx]);
  path_counts = thread_counts[thread_idx];
  emissions = thread_emissions[thread_idx];
  forward_scores = thread_forward[thread_idx];
  path = thread_path[thread_idx];
  if (thread_gff != NULL)
    query_gff = thread_gff[thread_idx];
  else
    query_gff = NULL;
 
  dms_lookup_emissions(dm, tuple_scores, emissions, blocks, seqnum,
		       msa->length, ih);
  if (zeroed_states != NULL)
    dms_zero_states(msa, emissions, zeroed_states);
  
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
			NULL, FALSE, xi_mode);
  
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
}

/* Worker function to launch sampling threads */
void dms_launch_sample_thread(void *data) {
  DMsamplingThreadData *thread = *((void**)data);

  dms_sample_path(thread->dm, thread->blocks, thread->ih, 
		  thread->tuple_scores, thread->thread_llh, 
		  thread->thread_path, thread->thread_emissions,
		  thread->thread_forward,
		  thread->thread_trans, thread->thread_counts,
		  thread->do_sample, thread->seqnum, thread->seqname, 
		  thread->log, thread->thread_gff, thread->do_reference, 
		  thread->nthreads, thread->p, thread->sample,
		  thread->zeroed_states, thread->xi_mode);  
}

void dms_read_priors(int **priors, FILE *prior_f, int xi_mode) {
  Regex *mu_re = str_re_new("#[[:space:]]*MU[[:space:]]*([0-9]+[[:space:]]*[0-9]+)");
  Regex *nu_re = str_re_new("#[[:space:]]*NU[[:space:]]*([0-9]+[[:space:]]*[0-9]+)");
  Regex *phi_re = str_re_new("#[[:space:]]*PHI[[:space:]]*([0-9]+[[:space:]]*[0-9]+)");
  Regex *zeta_re = str_re_new("#[[:space:]]*ZETA[[:space:]]*([0-9]+[[:space:]]*[0-9]+)");
  Regex *xi_re = str_re_new("#[[:space:]]*XI[[:space:]]*([0-9]+[[:space:]]*[0-9]+)");

  String *line = str_new(STR_MED_LEN);
  List *matches = lst_new_ptr(2);
  List *counts = lst_new_ptr(2);
  int count, xi_found = FALSE;

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
    } else  if (str_re_match(line, xi_re, matches, 1) >= 0) {
      if (xi_mode == FALSE) {
	lst_free_strings(matches);
	continue;
      }
      str_split(lst_get_ptr(matches, 1), NULL, counts);
      str_as_int(lst_get_ptr(counts, 0), &count);
      priors[4][0] = count;
      str_as_int(lst_get_ptr(counts, 1), &count);
      priors[4][1] = count;
      xi_found = TRUE;
    } else { /* something unexpected */
      lst_free_strings(matches);
      lst_free_strings(counts);
      die("ERROR: Too many lines in priors file.\n");
    }
    lst_free_strings(matches);
    lst_free_strings(counts);
  }

  if (xi_mode == TRUE && xi_found == FALSE)
    die("ERROR: Priors file must contain a line for xi when xi parameter is used.\n");

  str_re_free(mu_re);
  str_re_free(nu_re);
  str_re_free(phi_re);
  str_re_free(zeta_re);
  str_re_free(xi_re);
  str_free(line);
  lst_free_strings(matches);
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
  Regex *rc_re = str_re_new("_rc");

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
  if (str_re_search(tmp_str, rc_re, tmp_str->length-3, NULL, 1) >= 0) {
    sprintf(strand, "%s", "-");
    msa = lst_get_ptr(blocks->source_msas, seqnum);
    width = (end - start) + 1;
    start = msa->length - (end + 1);
    end = start + width - 1;
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

  if (matches != NULL) {
    tmatch_idx = (int*)smalloc(lst_size(target_gff->features) == 0 ? 1 :
			       lst_size(target_gff->features)
			       * sizeof(int));
    qmatch_idx = (int*)smalloc(lst_size(query_gff->features) == 0 ? 1 :
			       lst_size(query_gff->features)
			       * sizeof(int));
  
    if (mismatches == NULL || unique_to_query == NULL ||
	unique_to_target == NULL)
      die("ERROR in dms_compare_gffs: Not enough GFF_Sets for format!\n");
  }

  nmatches = 0;
  for (i = 0; i < lst_size(query_gff->features); i++) {
    fq = lst_get_ptr(query_gff->features, i);
    for (j = 0; j < lst_size(target_gff->features); j++) {
      ft = lst_get_ptr(target_gff->features, j);
      if (str_compare(fq->seqname, ft->seqname) != 0) /* Different sequences */
	continue;
      if ((fq->start + offset) == ft->start) {
	/* Motif found in same position */
	nmatches++;
	if (matches != NULL) {
	  tmatch_idx[nmatches] = j;
	  qmatch_idx[nmatches] = i;
	}
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
  if (matches != NULL) {
    free(tmatch_idx);
    free(qmatch_idx);
  }
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
    lst_free_strings(matches);
  }
  str_re_free(ne_re);
  str_re_free(ns_re);
  lst_free_strings(matches);
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
  lst_free(keys);
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
  dm_compute_log_likelihood(thread->mod, thread->msa, thread->emissions_vec);
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
			   int seqlen, int *ref_path, int force_priors,
			   int xi_mode) {
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
       and don't give much information as to the param values. Using the
       stationary frequencies of each state as the begin transition probs,
       though probably reasonable, is a bit arbitrary. */
    if (s1 == BEGIN_STATE) {
      s1 = s2;
      continue;
    }

    p1 = dm->state_to_motifpos[s1] ;
    e1 = dm->state_to_event[s1];
    p2 = dm->state_to_motifpos[s2];
    e2 = dm->state_to_event[s2];

    /* Track different types of transitions for param estimation */

    /* Transitions that affect mu -- i.e., from non-neutral bg states and 
       motif end positions */
    if (e1 != NEUT && (p1 == -1 || p1 == dm->m->width - 1)) {
      if (e2 == NEUT)
	trans[0][0]++;
      else
	trans[0][1]++;
    }

    /* Transitions that affect nu -- all of these are from NB state */
    if (e1 == NEUT && p1 == -1) {
      if (e2 != NEUT)
	trans[1][0]++;
      else
	trans[1][1]++;
    }
      
    /* Transitions that affect phi -- all are from NB. Must pay attention to
       whether destination state is conserved or not. */
    if (e1 == NEUT && e2 != NEUT) {
      if (e2 != CONS)
	trans[2][0]++;
      else
	trans[2][1]++;
    }

    /* Transitions that affect zeta or xi. All these are from bg states, which
       may be either neutral, conserved or lineage-specific. If we're using
       xi, also pay attention to whether we're starting a neutral motif or a
       conserved/gain/loss motif. */
    if (p1 == -1) {
      if (p2 == 0) { /* Started a motif -- add an alpha count */

	if (xi_mode == TRUE) {
	  if (e2 == NEUT)
	    trans[4][0]++;
	  else
	    trans[3][0]++;
	} else {
	   trans[3][0]++;
	}

      } else { /* Stayed in background -- add a beta count */

	if (xi_mode == TRUE) {
	  if (e2 == NEUT)
	     trans[4][1]++;
	  else
	     trans[3][1]++;
	} else {
	   trans[3][1]++;
	}
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
	end = i + w - 1;
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
		   int nwins, int xi_mode) { 
  int i;
  static int *stats = NULL;
  double fpr, spc, ppv;
  char param[STR_SHORT_LEN];
 
  fprintf(log, "Sample = %d, total log likelihood = %f\n", sample, llh);
  
  /* If using a reference set, compare known and predicted features and 
     keep track of the matches, mismatches, etc. in the stats structure */
  if (reference != NULL) {
    fpr = spc = ppv = 0;

    if (stats == NULL) {
      stats = (int*)smalloc(4 * sizeof(int));
    }
    for (i = 0; i < 4; i++)
      stats[i] = 0;
    
    dms_compare_gffs(reference, query_gff, stats, 0, NULL, NULL, NULL, 
		     NULL);
    fpr = (double)(stats[1] + stats[3])
       / (double)(nwins - lst_size(reference->features));
    spc = (double)(nwins - lst_size(query_gff->features))
      / (double)((nwins - lst_size(query_gff->features)) + stats[3] + stats[1]);
    if (stats[0] != 0)
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
  
  for (i = 0; i < 5; i++) {
    if (xi_mode == FALSE && i == 4)
      continue;
    if (i == 0)
      snprintf(param, STR_SHORT_LEN, "mu");
    else if (i == 1)
      snprintf(param, STR_SHORT_LEN, "nu");
    else if (i == 2)
      snprintf(param, STR_SHORT_LEN, "phi");
    else if (i == 3)
      snprintf(param, STR_SHORT_LEN, "zeta");
    else /* (i == 4) */
      snprintf(param, STR_SHORT_LEN, "xi");
    
    fprintf(log, "%s transitions: alpha = %d, beta = %d\n", param, 
	    trans[i][0], trans[i][1]);
  }
  
  if (xi_mode == FALSE) {
    fprintf(log, "mu = %f, nu = %f, phi = %f, zeta = %f\n\n", dm->mu, 
	    dm->nu, dm->phi, dm->zeta);
  } else {
    fprintf(log, "mu = %f, nu = %f, phi = %f, zeta = %f, xi = %f\n\n", dm->mu, 
	    dm->nu, dm->phi, dm->zeta, dm->xi);
  }
}

DMotifPmsaStruct *dms_read_alignments(FILE *F, int do_ih, int quiet, 
				      int revcomp, int do_zeroed,
				      FILE *cond_spec_f) {
  
  Regex *blocks_re = str_re_new("#[[:space:]]*BLOCKS[[:space:]]*=[[:space:]]*([0-9]+)");
  Regex *alph_re = str_re_new("#[[:space:]]*ALPHABET[[:space:]]*=[[:space:]]*([A-Z]+)");
  Regex *format_re = str_re_new("#[[:space:]]*FORMAT[[:space:]]*=[[:space:]]*([A-Z]+)");
  Regex *comment_re = str_re_new("#");
  
  int i, nblocks, ncats, tuple_size;
  char msa_fname[STR_MED_LEN], ih_fname[STR_MED_LEN],
    zeroed_fname[STR_MED_LEN];
  FILE *msa_f, *ih_f, *zeroed_f;
  msa_format_type format;
  List *msas = NULL, *matches = lst_new_ptr(5);
  MSA *msa, *msa_rc;
  String *line = str_new(STR_MED_LEN), *alphabet = NULL, 
    *fname = str_new(STR_MED_LEN);
  DMotifPmsaStruct *dmpmsa;

  dmpmsa = (DMotifPmsaStruct*)smalloc(sizeof(DMotifPmsaStruct));
  
  nblocks = i = dmpmsa->max_seqlen = 0;

  while (str_readline(line, F) != EOF) {
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
      } else if (str_re_match(line, comment_re, matches, 1) >= 0) {
	continue;
      } else {
	die("Bad header in alignment list file\n");
      }
      
      if (nblocks != 0 && alphabet != NULL && i == 3) {
	msas = lst_new_ptr(nblocks);
	dmpmsa->seqnames = lst_new_ptr(nblocks);
	if (do_ih) { /* Allocate space for indel histories, if
			called for -- this will need to be filled
			in explicitly later */
	  dmpmsa->ih = smalloc(nblocks * sizeof(IndelHistory*));
	} else {
	  dmpmsa->ih = NULL;
	}
	if (do_zeroed) {
	  dmpmsa->zeroed_states = smalloc((nblocks+1) * sizeof(List*));
	  if (cond_spec_f != NULL) {
	    /* zeroed_states[nblocks] is reserved for zeroed rows that affect
 	       all seqs */
	    dmpmsa->zeroed_states[nblocks] = dms_read_zeroed_states(cond_spec_f);
	  } else {
	    dmpmsa->zeroed_states[nblocks] = NULL;
	  }
	} else {
	  dmpmsa->zeroed_states = NULL;
	}
	i = 0;
      }
    } else {
      if (i >= nblocks)
	  die("Too many alignment files for format\n");
      
      if (str_split(line, NULL, matches) > 1) {
	strncpy(msa_fname, ((String*)lst_get_ptr(matches, 0))->chars,
		STR_MED_LEN);
	strncpy(ih_fname, ((String*)lst_get_ptr(matches, 1))->chars,
		STR_MED_LEN);
	if (lst_size(matches) > 2)
	  strncpy(zeroed_fname, ((String*)lst_get_ptr(matches, 2))->chars,
		  STR_MED_LEN);
      } else {
	strncpy(msa_fname, line->chars, STR_MED_LEN);
      }
      
      if (!quiet)
	fprintf(stderr, "\tReading alignment from %s  (%d of %d)\n", 
		msa_fname, (i+1), nblocks);
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

      if (do_zeroed && lst_size(matches) > 2) {
	if (!quiet)
	  fprintf(stderr, "\tReading conditioning data from %s...\n",
		  zeroed_fname);
	if (lst_size(matches) < 3)
	  die("ERROR: --cond-on-subs and --cond-on-species both require a zeroed-states file created with dmcondition. See 'dmsample -h' and 'dmcondition -h' for details.\n");
	zeroed_f = fopen_fname(zeroed_fname, "r");
	dmpmsa->zeroed_states[i] = dms_read_zeroed_states(zeroed_f);
	fclose(zeroed_f);
      }
      
      /* Create reverse complement of sequence i, if called for */
      if (revcomp == TRUE) {
	if (!quiet)
	  fprintf(stderr, "\tReverse complementing  %s  (%d of %d)\n",
		  msa_fname, (i+2), nblocks);
	msa_rc = msa_create_copy(msa, (format == SS ? 1 : 0));
	msa_reverse_compl(msa_rc);
	lst_push_ptr(msas, msa_rc);
        str_append_charstr(fname, "_rc");
	lst_push_ptr(dmpmsa->seqnames, str_dup(fname));
	if (do_zeroed) {
	  dmpmsa->zeroed_states[i+1] =
	    dms_reverse_zeroed_states(dmpmsa->zeroed_states[i], msa);
	}
	if (do_ih)
	  dmpmsa->ih[i+1] = dms_reverse_ih(dmpmsa->ih[i]);
	i++;
      }

      i++;
    }
    lst_free_strings(matches);
  }
  if (i < (nblocks - 1))
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

  /* Compact the zeroed_states lists, if needed */
  if (cond_spec_f != NULL)
    dms_compact_zeroed_states(dmpmsa->zeroed_states, nblocks);
  
  str_re_free(blocks_re);
  str_re_free(alph_re);
  str_re_free(format_re);
  str_re_free(comment_re);
  lst_free_strings(matches);
  lst_free(matches);
  str_free(alphabet);
  str_free(line);
  str_free(fname);
  return dmpmsa;
}

void dms_free_dmpmsa_struct(DMotifPmsaStruct *dmpmsa) {
  int i, j;
  List *msas;
  DMzeroedState *z;
  msas = dmpmsa->pmsa->source_msas;

  ss_free_pooled_msa(dmpmsa->pmsa);

  if(dmpmsa->ih != NULL) {
    for (i = 0; i < lst_size(msas); i++) {
      ih_free(dmpmsa->ih[i]);
    }
    free(dmpmsa->ih);
  }

  if (dmpmsa->zeroed_states != NULL) {
    for (i = 0; i <= lst_size(msas); i++) {
      if (dmpmsa->zeroed_states[i] == NULL)
	continue;
      for (j = 0; j < lst_size(dmpmsa->zeroed_states[i]); j++) {
	z = lst_get_ptr(dmpmsa->zeroed_states[i], j);
	if (z->do_row == 1 && i < lst_size(msas))
	  continue;
	dms_free_zeroed_state(z);
      }
      lst_free(dmpmsa->zeroed_states[i]);
    }
  }
  
  for (i = 0; i < lst_size(msas); i++) {
    msa_free(lst_get_ptr(msas, i));
    str_free(lst_get_ptr(dmpmsa->seqnames, i));
  }
  
  lst_free(dmpmsa->seqnames);
  free(dmpmsa);
  lst_free(msas);
}

/* Pared down version of tL_compute_log_likelihood with improved memory
   performance */
double dm_compute_log_likelihood(TreeModel *mod, MSA *msa, 
				 double *col_scores) {
  int i, j, k, nodeidx, tupleidx, defined, nstates, alph_size,
    *partial_match, skip_fels, thisseq, ninform, observed_state;
  double retval, total_prob, **inside_joint, *tuple_scores, totl, totr;
  TreeNode *n;
  List *traversal;
  MarkovMatrix *lsubst_mat, *rsubst_mat;

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

/*   fprintf(stderr, "mod->allow_gaps %d, mod->inform_reqd %d\n",  */
/* 	  mod->allow_gaps, mod->inform_reqd); */

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
	
	} else {                    
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
    } /* if !skip_fels */
    
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
/*     fprintf(stderr, "state %d, pos %d\n", i, dm->state_to_motifpos[i]); */
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
void dms_cache_hash(Hashtable *path_counts, char *cache_fname, 
			  int nstates, int nsamples) {
  FILE *hash_f;
  
  hash_f = fopen_fname(cache_fname, "w");
  dms_write_hash(path_counts, hash_f, nstates, nsamples);
  fclose(hash_f);
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
      fprintf(stderr, "\tReading sampling data from file %s...\n", fname);
    
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
  Regex *comment_re = str_re_new("#");
  int i, nfiles = -1;
  String *line = str_new(STR_MED_LEN);
  List *matches = lst_new_ptr(2), *ret;

  i = 0;
  while (str_readline(line, tmp_lst_f) != EOF) {
    lst_clear(matches);
    str_trim(line);
    if (line->length == 0) continue;
    
    if (str_re_match(line, nfiles_re, matches, 1) >= 0) {
      str_as_int(lst_get_ptr(matches, 1), &nfiles);
      if (nfiles > 0)
	ret = lst_new_ptr(nfiles);
      else
	die("ERROR: Bad header in temp file list!\n");
    } else if (nfiles == -1) {
      die("ERROR: Bad or mising header in temp file list!\n");
    } else if (str_re_match(line, comment_re, matches, 1) >= 0) {
      continue;
    } else {
      if (i == nfiles)
	die("ERROR: Too many lines in temp file list!\n");
      lst_push_ptr(ret, (void*)str_dup(line));
      i++;
    }
    lst_free_strings(matches);
  }
  if (i != nfiles) die("ERROR: Not enough lines in temp file list!\n");
  str_re_free(nfiles_re);
  str_re_free(comment_re);
  str_free(line);
  lst_free_strings(matches);
  lst_free(matches);
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

/* Zero out emissions for combiations of states and sequence positions where
   no valid motifs are possible */
void dms_zero_states(MSA *msa, double **emissions, List *zeroed_states) {
  int i, j, k, pos, len;
  DMzeroedState *z;

  for (i = 0; i < lst_size(zeroed_states); i++) {
    z = lst_get_ptr(zeroed_states, i);
    if (z->do_row == 1) {
/*       fprintf(stderr, "Zeroing emissions row for state %d\n", z->state); */
      dms_zero_emissions_row(emissions[z->state], msa->length);
    } else {
/*       fprintf(stderr, "state %d, do_row %d\n", z->state, z->do_row); */
      for (j = 0; j < lst_size(z->starts); j++) {
	pos = lst_get_int(z->starts, j);
	len = lst_get_int(z->lengths, j);
/* 	fprintf(stderr, "z->state %d, pos %d, len %d, pos + len %d, msa->length %d\n", */
/* 		z->state, pos, len, (pos + len), msa->length); */
	for (k = pos; k < (pos + len); k++) {
/* 	  if (k >= msa->length) */
/* 	    fprintf(stderr, "\tk %d, pos %d, len %d, pos + len %d\n", */
/* 		    k, pos, len, (pos + len)); */
/* 	  fprintf(stderr, "%.3f  ", emissions[z->state][k]); */
	  emissions[z->state][k] = -INFTY;
/* 	  fprintf(stderr, "%.3f\n", emissions[z->state][k]); */
	}
      }	
    } 
  }  
}

/* Zero out an entire row in the emissions matrix */
void dms_zero_emissions_row(double *row, int seqlen) {
  int i;
  for (i = 0; i < seqlen; i++) {
/*     fprintf(stderr, "%.3f  ", row[i]); */
    row[i] = -INFTY;
/*     fprintf(stderr, "%.3f\n", row[i]); */
  }
}

/* Create an array of DMzeroedState objects to emulate conditioning predictions
   on presence of a site in a given species. Right now this will only work for
   a single species, but should be possible to modify for conditioning on a
   site's presence in multiple species. */ 
DMzeroedState **dms_condition_on_species(DMotifPhyloHmm *dm, TreeNode *tree,
					 int cond_on /* Node ID of (leaf) 
							sequence believed to 
							contain a site */
					 ) {
  int i, j, state, nodes_visited, states_visited, *mark;
  TreeNode *n;
  DMzeroedState *z, **retval;

  retval = (DMzeroedState**)smalloc(dm->phmm->hmm->nstates *
				    sizeof(DMzeroedState*));
  for (i = 0; i < dm->phmm->hmm->nstates; i++) {
    retval[i] = NULL;
  }

  /* Create an array of ints to indicate whether to zero births or deaths for
     a given node id. 0 indicates births will be zeroed, 1 indicates deaths
     will be zeroed. */
  mark = smalloc(tree->nnodes * sizeof(int));
  for (i = 0; i < tree->nnodes; i++)
    mark[i] = 0;


  /* Starting at the leaf sequence we're conditioning on, climb the tree and
     mark nodes that are predecessors of the node we're conditioning on */
  n = lst_get_ptr(tree->nodes, cond_on);
  assert(n->lchild == NULL && n->rchild == NULL);
  while (n != tree) {
/*     fprintf(stderr, "n->id %d, n->name %s, parent->id %d\n",  */
/* 	    n->id, n->name, n->parent->id); */
    mark[n->id] = 1;
    n = n->parent;
  }
/*   fprintf(stderr, "\n"); */

  /* Traverse the nodes in the tree in arbitrary order and zero birth or death
     states according to the value or mark[n->id] */
  nodes_visited = states_visited = 0;
  for (i = 0; i < tree->nnodes; i++) {
    n = lst_get_ptr(tree->nodes, i);
    nodes_visited++;
    if (n == tree) {
      continue;
    } else if (mark[n->id] == 1) { /* Zero death states */
/*       fprintf(stderr, "n->id %d, n->name %s, zeroing deaths\n", */
/* 	      n->id, n->name); */
      for (j = 0; j < lst_size(dm->branch_to_states[n->id]); j++) {
	states_visited++;
	state = lst_get_int(dm->branch_to_states[n->id], j);
	if (dm->state_to_event[state] == DEATH &&
	    dm->state_to_motifpos[state] != -1) {
/* 	  fprintf(stderr, "\tstate %d (death, zeroed)\n", state); */
	  z = dms_new_zeroed_state(state, 1);
	  retval[state] = z;
	} else {
/* 	  fprintf(stderr, "\tstate %d (%s)\n",  */
/* 		  state, dm->state_to_motifpos[state] == -1 ? "bg" : "birth"); */
	}
      }
    } else { /* mark[n->id] == 0 -- zero birth states */
/*       fprintf(stderr, "n->id %d, n->name %s, zeroing births\n", */
/* 	      n->id, n->name); */
      for (j = 0; j < lst_size(dm->branch_to_states[n->id]); j++) {
	states_visited++;
	state = lst_get_int(dm->branch_to_states[n->id], j);
	if (dm->state_to_event[state] == BIRTH &&
	    dm->state_to_motifpos[state] != -1) {
/* 	  fprintf(stderr, "\tstate %d (birth, zeroed)\n", state); */
	  z = dms_new_zeroed_state(state, 1);
	  retval[state] = z;
	} else {
/* 	  fprintf(stderr, "\tstate %d (%s)\n", */
/* 		  state, dm->state_to_motifpos[state] == -1 ? "bg" : "death"); */
	}
      }
    }
  }

/*   fprintf(stderr, "nodes visited: %d of %d\n", nodes_visited, tree->nnodes); */
/*   fprintf(stderr, "states visited: %d of %d\n", states_visited, dm->phmm->hmm->nstates); */
  assert(states_visited == (dm->phmm->hmm->nstates - (2 * dm->m->width + 2)));
  assert(nodes_visited == tree->nnodes);
  free(mark);
  return retval;
}

/* Build a list of DMzeroedState objects to condition on presence
   of substitutions by zeroing chains of motif states for which there are no
   substitutions to support a gain/loss prediction. */
void dms_condition_on_subs(DMotifPhyloHmm *dm, TreeModel *mod, MSA *msa,
			   DMzeroedState **zeroed_states, int nosubs) {
  int i, j, k, n_id, state, start, lastpos, len, offset, **starts_ar;
  List **zeroed_nodes;
  DMzeroedState *z;

  /* Used to keep the starts lists as compact as possible */
  starts_ar = smalloc(dm->phmm->hmm->nstates * sizeof(int*));
  for (i = 0; i < dm->phmm->hmm->nstates; i++) {
    starts_ar[i] = smalloc(msa->length * sizeof(int));
    for (j = 0; j < msa->length; j++)
      starts_ar[i][j] = 0;
  }
  
  zeroed_nodes = dms_label_subst_nodes(msa, mod, dm->m, nosubs);

  for (i = 0; i < (msa->length - dm->m->width) + 1; i++) {
    if (zeroed_nodes[i] != NULL) {
      for (j = 0; j < lst_size(zeroed_nodes[i]); j++) {
	n_id = lst_get_int(zeroed_nodes[i], j);

	offset = 0;
	for (k = 0; k < lst_size(dm->branch_to_states[n_id]); k++) {
	  state = lst_get_int(dm->branch_to_states[n_id], k);
	  /* Only zero birth and death states. Here we DO NOT pay attention
	     to whether the state is a motif state; we don't want to count
	     any substitutionless gain/loss events within background states
	     because these may skew phi. Background states will need special
	     handling. */
	  if (dm->state_to_event[state] != BIRTH &&
	      dm->state_to_event[state] != DEATH)
	    continue;
	  
	  /* Check for an existing DMzeroedState object with a starts list and
	     create one if not found or continue if the whole row is already
	     being zeroed. */
	  z = zeroed_states[state];
	  if (z != NULL && z->do_row == 1) {
	    continue;
	  } else if (z == NULL) {
	    z = zeroed_states[state] = dms_new_zeroed_state(state, 0);
	    z->starts = lst_new_int(100); /* Wild guess at size needed */
	    z->lengths = lst_new_int(100);
	  } else {
	    if (z->starts == NULL)
	      die("ERROR: Expected starts list for state %d but found NULL!\n",
		  state);
	  }
	  
	  if (dm->state_to_motifpos[state] == -1) { /* bg state */
	    for (offset = 0; offset < dm->m->width; offset++)
	      starts_ar[state][i + offset] = 1;
	    offset = 0;
	  } else {
	    starts_ar[state][i + offset] = 1;
	    if (offset < dm->m->width - 1)
	      offset++;
	    else
	      offset = 0;
	  }
	}
      }
    }
  }
  
  /* Compact the states_ar arrays into the z->starts lists */
  for (i = 0; i < dm->phmm->hmm->nstates; i++) {
    z = zeroed_states[i];
    if (z == NULL || z->do_row == 1)
      continue;
    
    lastpos = 0;
    for (j = 0; j <= msa->length; j++) {

      if (j < msa->length) {
	if (starts_ar[i][j] == 1) {
	  if (lastpos == 0) { /* Starting a new block of zeroed pos's */
	    start = j;
	    len = 1;
	    lastpos = 1;
	  } else {
	    len++;
	  }
	} else { /* (starts_ar[i][j] == 0)  */
	  if (lastpos == 1) { /* Ending a block of zeroed pos's */
	    lst_push_int(z->starts, start);
	    lst_push_int(z->lengths, len);
	    lastpos = 0;
	  } else {
	    continue;
	  }
	}
      } else { /* (j == msa->length) */
	if (lastpos == 1) { /* Ending a block of zeroed pos's */
	  lst_push_int(z->starts, start);
	  lst_push_int(z->lengths, len);
	}
      }
    }
  }
  
  for (i = 0; i < dm->phmm->hmm->nstates; i++) {
    if (starts_ar[i] != NULL)
      free(starts_ar[i]);
  }
  free(starts_ar);
}

/* Create a new DMzeroedState object */
DMzeroedState *dms_new_zeroed_state(int state, int do_row) {
  DMzeroedState *z = smalloc(sizeof(DMzeroedState));
  z->state = state;
  z->do_row = do_row;
  z->starts = NULL;
  z->lengths = NULL;
  return z;
}

/* Free a DMzeroedState object */
void dms_free_zeroed_state(DMzeroedState *z) {
  if (z->starts != NULL) {
    lst_free(z->starts);
    lst_free(z->lengths);
  }
  free(z);
}

void dms_write_zeroed_states(FILE *f, List *zeroed_states) {
  int i, j, start, len;
  DMzeroedState *z;

  fprintf(f, "# NZEROED_STATES = %d\n", lst_size(zeroed_states));
  for (i = 0; i < lst_size(zeroed_states); i++) {
    z = lst_get_ptr(zeroed_states, i);
    fprintf(f, "%d\t%d\t", z->state, z->do_row);
    if (z->starts != NULL) {
      for (j = 0; j < lst_size(z->starts); j++) {
	start = lst_get_int(z->starts, j);
	fprintf(f, "%d,", start);
      }
      fprintf(f, "\t");
      for (j = 0; j < lst_size(z->lengths); j++) {
	len = lst_get_int(z->lengths, j);
	fprintf(f, "%d,", len);
      }
    }
    fprintf(f, "\n");
  }
}

List *dms_read_zeroed_states(FILE *f) {
  Regex *nzeroed_states_re = str_re_new("#[[:space:]]*NZEROED_STATES[[:space:]]*=[[:space:]]*([0-9]+)");
  int i, j, state, do_row, start, len, nzeroed = -1;
  char delim[2];
  List *retval, *starts, *lengths, *matches, *tmpl;
  String *line = str_new(STR_LONG_LEN), *starts_str, *lens_str;
  DMzeroedState *z;

  matches = lst_new_ptr(4);
  tmpl = lst_new_ptr(5000); /* Connservative guess at the number of starts,
				to avoid having to realloc this list, which
				seems to be a major contributor to bus error
				crashes! */
  snprintf(delim, 2, ",");

  i = 0;
  while (str_readline(line, f) != EOF) {
    lst_clear(matches);
    str_trim(line);
    if (line->length == 0) continue;

    if (str_re_match(line, nzeroed_states_re, matches, 1) >= 0) {
      str_as_int(lst_get_ptr(matches, 1), &nzeroed);
      if (nzeroed > 0)
	retval = lst_new_ptr(nzeroed);
      else
	die("ERROR: Bad header in zeroed_states file!\n");
    } else if (nzeroed == -1) {
      die("ERROR: Bad or mising header in zeroed_states file!\n");
    } else {
      if (i == nzeroed)
	die("ERROR: Too many lines in zeroed_states file!\n");
      str_split(line, NULL, matches);

      str_as_int((String*)lst_get_ptr(matches, 0), &state);
      str_as_int((String*)lst_get_ptr(matches, 1), &do_row);

      if (do_row == 0) {
	starts_str = lst_get_ptr(matches, 2);
	str_split(starts_str, delim, tmpl);
	starts = lst_new_int(lst_size(tmpl));
	for (j = 0; j < lst_size(tmpl); j++) {
	  str_as_int(lst_get_ptr(tmpl, j), &start);
	  lst_push_int(starts, start);
	}
	lens_str = lst_get_ptr(matches, 3);
	str_split(lens_str, delim, tmpl);
	lengths = lst_new_int(lst_size(tmpl));
	for (j = 0; j < lst_size(tmpl); j++) {
	  str_as_int(lst_get_ptr(tmpl, j), &len);
	  lst_push_int(lengths, len);
	}
      } else {
	starts = NULL;
	lengths = NULL;
      }

      z = dms_new_zeroed_state(state, do_row);
      z->starts = starts;
      z->lengths = lengths;
      lst_push_ptr(retval, z);

      i++;
    }
    lst_free_strings(matches);
  }

  if (i != nzeroed) die("ERROR: Not enough lines in zeroed_states file!\n");
  lst_free_strings(matches);
  lst_free(matches);
  lst_free(tmpl);
  str_free(line);
  str_re_free(nzeroed_states_re);
  return retval;
}

/* Convert an array of DMzeroedState objects into a List. Copies pointers to
   DMzeroedState objects directly from the array -- does not make copies of the
   objects themselves, so NOT safe to free the objects after calling this
   function, but can free zeroed_states array safely. */
List *dms_zeroed_states_array_to_list(DMzeroedState **zeroed_states,
				      int nstates, int init_size) {
  int i;
  List *retval = lst_new_ptr(init_size);

  for (i = 0; i < nstates; i++) {
    if (zeroed_states[i] != NULL)
      lst_push_ptr(retval, zeroed_states[i]);
  }

  return retval;
}

/* Given a dm object, a tree and a zereoed_states list, print the states,
   branches and events for zeroedState objects in the list to the given
   stream */
void dms_print_zeroed_states(FILE *f, DMotifPhyloHmm *dm, TreeNode *tree,
			     List *zeroed_states) {
  int i, j;
  TreeNode *n;
  DMzeroedState *z;

  for (i = 0; i < lst_size(zeroed_states); i++) {
    z = lst_get_ptr(zeroed_states, i);
    n = lst_get_ptr(tree->nodes, dm->state_to_branch[z->state]);
    fprintf(f, "i = %d, state = %d, motifpos %d, n->id = %d, n->name = %s, event = %s, do_row = %d\n",
	    i, z->state, dm->state_to_motifpos[z->state], n->id, n->name,
	    dm->state_to_event[z->state] == BIRTH ? "birth" :
	    dm->state_to_event[z->state] == DEATH ? "death" : "other",
	    z->do_row);
    if (z->starts != NULL) {
      fprintf(f, "\tstarts: ");
      for (j = 0; j < lst_size(z->starts); j++) {
	fprintf(f, "%d,", lst_get_int(z->starts, j));
      }
      fprintf(f, "\n");
      fprintf(f, "\tlengths: ");
      for (j = 0; j < lst_size(z->lengths); j++) {
	fprintf(f, "%d,", lst_get_int(z->lengths, j));
      }
      fprintf(f, "\n");
    }
  }
}

/* Reverse POSITIONS ONLY in a List of DMzeroedState objects */
List *dms_reverse_zeroed_states(List *zeroed_states, MSA *msa) {
  int i, j, start, len;
  DMzeroedState *z, *z_copy;
  List *retval, *starts_copy = NULL, *lens_copy = NULL;

  retval = lst_new_ptr(lst_size(zeroed_states));

  for (i = 0; i < lst_size(zeroed_states); i++) {
    z = lst_get_ptr(zeroed_states, i);
    z_copy = dms_new_zeroed_state(z->state, z->do_row);
    if (z->starts != NULL) {
      starts_copy = lst_new_int(lst_size(z->starts));
      lens_copy = lst_new_int(lst_size(z->lengths));
      for (j = lst_size(z->starts) - 1; j >= 0; j--) {
	start = lst_get_int(z->starts, j);
	len = lst_get_int(z->lengths, j);
	lst_push_int(starts_copy, (msa->length - start - len));
	lst_push_int(lens_copy, len);
      }
      z_copy->starts = starts_copy;
      z_copy->lengths = lens_copy;
    }
    lst_push_ptr(retval, z_copy);
  }

  return retval;
}

/* Given an MSA, Tree and Motif models, label the nodes of the tree that lack
   substitutions for zeroing of emissions as a way of conditioning gain and
   loss events on presence of substitutions to support the predictions. Uses
   a simple Fitch parsimony algorithm. */
List **dms_label_subst_nodes(MSA *msa, TreeModel *mod, PSSM *m, int nosubs) {
  int i, j, k, tuple_idx, seq_idx, overlap, p, *mark, all_missing, *zeroed;
  char *key;
  MSA *sub_msa;
  Hashtable *tuple_hash;
  TreeNode *n;
  List *traversal, *lchild_seqs, *rchild_seqs, *seqs, **pre, **post, **retval;
  
  /* The return value is an array of pointers to List objects describing nodes
     that need to be zeroed for each overlapping motif window in the dataset
     or null pointers for windows that need no nodes zeroed. */
  retval = smalloc(((msa->length - m->width) + 1) * sizeof(List*));
  for (i = 0; i < (msa->length - m->width) + 1; i++)
    retval[i] = NULL;
  zeroed = (int*)smalloc(mod->tree->nnodes * sizeof(int));
  mark = smalloc(msa->nseqs * sizeof(int));
  tuple_hash = hsh_new(mod->tree->nnodes);
  pre = smalloc(mod->tree->nnodes * sizeof(List*));
  post = smalloc(mod->tree->nnodes * sizeof(List*));
  for (i = 0; i < mod->tree->nnodes; i++) {
    pre[i] = lst_new_int(mod->tree->nnodes);
    post[i] =lst_new_int(mod->tree->nnodes);
  }
  
  if (mod->msa_seq_idx == NULL)
    tm_build_seq_idx(mod, msa);

  /* Find nodes for all overlapping motif windows that do not have subsitutions
     to support gain/loss predictions. */
  for (i = 0; i < (msa->length - m->width) + 1; i++) {
    
    /* First encode each leaf string as a single character */
    hsh_clear(tuple_hash);
    tuple_idx = 0;
    /* Initialize mark to all missing data */
    for (j = 0; j < msa->nseqs; j++)
      mark[j] = -1;
    /* Get the alignment for the current window */
    sub_msa = msa_sub_alignment(msa, NULL, 0, i, (i + m->width));
    if (sub_msa->seqs == NULL)
      ss_to_msa(sub_msa);
    
    for (j = 0; j < mod->tree->nnodes; j++) {
      n = lst_get_ptr(mod->tree->nodes, j);
      if (n->lchild != NULL) /* internal node */
	continue;
      seq_idx = mod->msa_seq_idx[n->id];

      key = sub_msa->seqs[seq_idx];

      /* Exclude species with all missing data in this window */
      all_missing = 1;
      for (k = 0; k < m->width; k++) {
	if (!msa->is_missing[(int)key[k]] /* && key[k] != GAP_CHAR */) {
	  /* FIXME: It is not clear if windows of all gap characters should
	     be considered missing and zeroed -- there are (at least) two
	     things these could be: deletions of entire motifs, which could be
	     considered legitimate losses, or regions of poor alignment
	     quality. Quality-masked MAF's may resolve the latter issue. For
	     the moment, do not zero these -- instead use quality-masked MAF
	     files. */
	  all_missing = 0;
	  break;
	}
      }
      if (all_missing == 1) {
	mark[seq_idx] = -1;
/* 	fprintf(stderr, "n->name %s, key %s\n", n->name, key); */
	continue;
      }
/*       fprintf(stderr, "pos %d, spec %d, key %s\n", i, j, key); */
      if (hsh_get(tuple_hash, key) == (void*)-1) {
	mark[seq_idx] = tuple_idx;
	hsh_put(tuple_hash, key, (void*)&mark[seq_idx]);
	tuple_idx++;
/*         fprintf(stderr, "mark[%d] %d\n", j, mark[j]); */
      } else {
	mark[seq_idx] = *(int*)hsh_get(tuple_hash, key);
/* 	fprintf(stderr, "mark[%d] %d\n", j, mark[j]); */
      }
    }
    msa_free(sub_msa);

/*     for (j = 0; j < msa->nseqs; j++) */
/*       fprintf(stderr, "mark[%d] %d, ", j, mark[j]); */
/*     fprintf(stderr, "\n"); */
    
    /* Do the postorder traversal step of the Fitch parsimony algorithm */
    traversal = tr_postorder(mod->tree);
    for (j = 0; j < lst_size(traversal); j++) {
      n = lst_get_ptr(traversal, j);
      seqs = pre[n->id];
      lst_clear(seqs);
      if (n->lchild == NULL) { /* base case: leaf node */
	seq_idx = mod->msa_seq_idx[n->id];
	if (mark[seq_idx] == -1) 
	  continue;
	lst_push_int(seqs, mark[seq_idx]);
      } else { /* internal node or root */
	lchild_seqs = pre[n->lchild->id];
	rchild_seqs = pre[n->rchild->id];
	dms_intersect_or_union(seqs, lchild_seqs, rchild_seqs);

/* 	fprintf(stderr, "i %d, j %d\n", i, j); */
/* 	fprintf(stderr, "\tlchild_seqs: "); */
/* 	for (k = 0; k < lst_size(lchild_seqs); k++) */
/* 	  fprintf(stderr, "%d,", lst_get_int(lchild_seqs, k)); */
/* 	fprintf(stderr, "\n"); */
/* 	fprintf(stderr, "\trchild_seqs: "); */
/* 	for (k = 0; k < lst_size(rchild_seqs); k++) */
/* 	  fprintf(stderr, "%d,", lst_get_int(rchild_seqs, k)); */
/* 	fprintf(stderr, "\n"); */
/* 	fprintf(stderr, "\tseqs: "); */
/* 	for (k = 0; k < lst_size(seqs); k++) */
/* 	  fprintf(stderr, "%d,", lst_get_int(seqs, k)); */
/* 	fprintf(stderr, "\n\n"); */
      }
    }
    
    /* As we go down the tree, we want to zero out scenarios where there
       are NO parsimonious scenarios that REQUIRE a substitution. This is much
       simpler than it sounds: because in the downward pass, lists at all child
       nodes will intersect. */
    traversal = tr_preorder(mod->tree);
    for (j = 0; j < lst_size(traversal); j++) {
      n = lst_get_ptr(traversal, j);
      lst_clear(post[n->id]);

      if (n == mod->tree) { /* Root node -- base case */
	for (k = 0; k < lst_size(pre[n->id]); k++) {
	  p = lst_get_int(pre[n->id], k);
	  lst_push_int(post[n->id], p);
	}
	continue;
      } else {
	lst_clear(post[n->id]);
	overlap = dms_intersect_or_union(post[n->id], post[n->parent->id],
					 pre[n->id]);
	/* Since dms_intersect_or_union will take the union if there is no
	   intersection between parent and child, which may cause incorrect
	   behavior at internal nodes. dms_intersect_or_union returns 0 in
	   these cases, so we can take corrective measures to avoid any
	   problems. */
	if (overlap == 0) { /* union was taken */
	  lst_clear(post[n->id]);
	  for (k = 0; k < lst_size(pre[n->id]); k++) {
	    p = lst_get_int(pre[n->id], k);
	    lst_push_int(post[n->id], p);
	  }
	}
	
	if (nosubs == 1) {
	  /* Use comparison method that zeroes branches only when no
	     substitution is possible along a branch */
	  zeroed[n->id] = dms_compare_lists_nosubs(post[n->parent->id],
						   post[n->id]);
	} else {
	  /* Use comparison method that zeroes branches when parent and
	     child nodes contain sets of characters that are identical in
	     size and content. */
	  zeroed[n->id] = dms_compare_lists_identity(post[n->parent->id],
						     post[n->id]);
	}
	
	/* Check for leaf nodes with all missing data and zero branches leading
	   to these nodes also. */
	if (n->lchild == NULL) {
	  seq_idx =  mod->msa_seq_idx[n->id];
	  if (mark[seq_idx] == -1)
	    zeroed[n->id] = -1;
	}

/* 	fprintf(stderr, "i %d, j %d\n", i, j); */
/* 	fprintf(stderr, "\tpost[%d] (%s): ", n->id, n->name); */
/* 	for (k = 0; k < lst_size(post[n->id]); k++) */
/* 	  fprintf(stderr, "%d,", lst_get_int(post[n->id], k)); */
/* 	fprintf(stderr, "\n"); */
/* 	fprintf(stderr, "\tparent (%d) (%s): ", n->parent->id, n->parent->name); */
/* 	for (k = 0; k < lst_size(post[n->parent->id]); k++) */
/* 	  fprintf(stderr, "%d,", lst_get_int(post[n->parent->id], k)); */
/* 	fprintf(stderr, "\n"); */
/* 	fprintf(stderr, "\tpre[%d] (%s): ", n->id, n->name); */
/* 	for (k = 0; k < lst_size(pre[n->id]); k++) */
/* 	  fprintf(stderr, "%d,", lst_get_int(pre[n->id], k)); */
/* 	fprintf(stderr, "\n"); */
/* 	fprintf(stderr, "zero_branch %d\n", zero_branch); */
/* 	if (zero_branch == 1) { */
/* 	  fprintf(stderr, " Zeroed: "); */
/* 	  for (k = 0; k < lst_size(retval[i]); k++) */
/* 	    fprintf(stderr, "%d,", lst_get_int(retval[i], k)); */
/* 	} */
/* 	fprintf(stderr, "\n"); */

      }
    }
    
    /* Go through the zeroed array and check for subtrees of all missing data
       -- zero branches leading to these subtrees. */
    traversal = tr_postorder(mod->tree);
    for (j = 0; j < lst_size(traversal); j++) {
      n = lst_get_ptr(traversal, j);
      if (n->lchild == NULL)
	continue;
      if (zeroed[n->lchild->id] == -1 && zeroed[n->rchild->id] == -1)
	zeroed[n->id] = -1;
    }
    
    /* Push zeroed node id's onto the return list */
    for (j = 0; j < lst_size(traversal); j++) {
      n = lst_get_ptr(traversal, j);
      if (retval[i] == NULL)
	retval[i] = lst_new_int(mod->tree->nnodes);
      if (zeroed[n->id] != 0)
	lst_push_int(retval[i], n->id);
    }
    
  }
  
  for (i = 0; i < mod->tree->nnodes; i++) {
    if (pre[i] != NULL)
      lst_free(pre[i]);
    if (post[i] != NULL)
      lst_free(post[i]);
  }
  free(pre);
  free(post);
  free(mark);
  free(zeroed);
  hsh_free(tuple_hash);
  return retval;
}

/* Create a list that represents the intersection of two int lists or, if there
   is no overlap between the two lists, the union of both lists. *dest may be
   NULL, in which case the function only returns the size of the overlap
   or 0 if there is none. */
int dms_intersect_or_union(List *dest, List *left, List *right) {
  int i, j, l, r, found;
  List *l1, *l2;

  assert((dest == NULL || lst_empty(dest)));
  
  if (lst_size(left) > lst_size(right)) {
    l1 = left;
    l2 = right;
  } else {
    l1 = right;
    l2 = left;
  }

  /* Handle simple cases in which no comparison is needed */
  if (lst_size(l1) == lst_size(l2) && lst_size(l2) == 0) {
    return 0;
  } else if (lst_size(l2) == 0) {
    if (dest == NULL)
      return 0;

    for (i = 0; i < lst_size(l1); i++) {
      l = lst_get_int(l1, i);
      lst_push_int(dest, l);
    }
    return 0;
  }
  
  /* Search for and keep track of the number of elements in the intersection */
  found = 0;
  for (i = 0; i < lst_size(l1); i++) {
    l = lst_get_int(l1, i);

    for (j = 0; j < lst_size(l2); j++) {
      r = lst_get_int(l2, j);
      if (l == r) {
	if (dest != NULL)
	  lst_push_int(dest, l);
	found++;
	break;
      }
    }
  }

  /* No intersection -- take the union of the two lists */
  if (found == 0 && dest != NULL) {
    for (i = 0; i < lst_size(l1); i++) {
      l = lst_get_int(l1, i);
      lst_push_int(dest, l);
    }
    for (i = 0; i < lst_size(l2); i++) {
      r = lst_get_int(l2, i);
      lst_push_int(dest, r);
    }
  }
  return found;
}

/* Given an array of lists of DMzeroedState objects, make all DMzeroedState
   pointers to states with entire rows zeroed out (i.e., do_row == 1) within
   List** indeces above 0 pointers to the corresponding object in List**
   index 0.

   WARNING: This assumes that zeroed_states List items for each emissions
   row needing zeroing are in the same order in all Lists!!! */
void dms_compact_zeroed_states(List **zeroed_states, int nblocks) {
  int i, j, k;
  DMzeroedState *z, *z_row;

  for (i = 0; i < nblocks; i++) {
    k = 0;
    for (j = 0; j < lst_size(zeroed_states[i]); j++) {
      z = lst_get_ptr(zeroed_states[i], j);
      if (z->do_row == 1) {
	z_row = lst_get_ptr(zeroed_states[nblocks], k);
/* 	fprintf(stderr, "k %d, z->state %d, z_row->state %d\n", */
/* 		k, z->state, z_row->state); */
	dms_free_zeroed_state(z);
	lst_set_ptr(zeroed_states[i], j, z_row);
	k++;
      }
    }
    assert(k == lst_size(zeroed_states[nblocks]));
  }
}

/* If assume that the only branches we want to zero are those where there
   is no possibility of substitution along that branch,
   the only cases where we need to zero out a branch will be
   indicated by a single character in the parent that is identical
   to the character at its child after the intersection at the child
   node is taken -- all other cases allow/require substitutions on
   the connecting branch. */
int dms_compare_lists_nosubs(List *parent, List *child) {
  int p, c;

  /* We only need to pay attention if the parent list is size 1, otherwise
     there exist scenarios where substitutions may occur along the branch */
  if (lst_size(parent) == 1) {
    if (lst_size(child) > 1) { /* There can be no intersection with parent or
				  the child list would be size 1. */
      return 0;
    } else { /* Need to see if parent char matches child char */
      p = lst_get_int(parent, 0);
      c = lst_get_int(child, 0);
      if (p == c)
	return 1;
    }
  }
  return 0;
}

/* A less restrictive case is to zero branches where the parent and child
   have identical sets of characters. */
int dms_compare_lists_identity(List *parent, List *child) {
  int overlap;

  if (lst_size(parent) != lst_size(child)) {
    /* lists cannot be identical */
    return 0;
  } else if (lst_size(parent) == 1) {
    /* use simple comparison */
    return dms_compare_lists_nosubs(parent, child);
  } else {
    overlap = dms_intersect_or_union(NULL, parent, child);
    if (overlap == lst_size(parent)) {
      /* the two lists are identical -- they have 100% overlap */
      return 1;
    }
  }
  return 0;
}

/* Print out zeroed regions in bed format. Skips zeroed rows. */
void dms_zeroed_as_bed(FILE *f, DMotifPhyloHmm *dm, List *zeroed_states,
		       char *fname) {
  int i, j, start, end, len, cat, genstart;
  char delim[2], chrom[10];
  List *split;
  String *fname_str, *name, *chrom_str = str_new(STR_MED_LEN);
  CategoryMap *cm = dm->phmm->cm;
  DMzeroedState *z;

  snprintf(delim, 2, ".");
  split = lst_new_ptr(4);
  fname_str = str_new_charstr(fname);
  str_split(fname_str, delim, split);
  str_as_int((String*)lst_get_ptr(split, 1), &genstart);
  str_append(chrom_str, (String*)lst_get_ptr(split, 0));
  str_remove_path(chrom_str);
  snprintf(chrom, 10, chrom_str->chars);

  for (i = 0; i < lst_size(zeroed_states); i++) {
    z = lst_get_ptr(zeroed_states, i);
    if (z->do_row == 1 || dm->state_to_motifpos[z->state] > 0)
      continue;

    name = str_new_charstr(fname);
    str_remove_path(name);
    str_append_charstr(name, "_");
    cat = cm->ranges[dm->phmm->state_to_cat[z->state]]->start_cat_no;
     str_append(name, (String*)cm_get_feature(cm, cat));
    
    for (j = 0; j < lst_size(z->starts); j++) {
      start = lst_get_int(z->starts, j);
      len = lst_get_int(z->lengths, j);
      if (dm->state_to_motifpos[z->state] == 0)
	end = (start + len);
      else
	end = (start + len);
      
      fprintf(f, "%s\t%d\t%d\t%s\n", chrom, (start + genstart - 1),
	      (end + genstart - 1), name->chars);
    }
    str_free(name);
  }
  str_free(fname_str);
  str_free(chrom_str);
  lst_free(split);
}

/* Condition transition probabilities in the hmm to require site presence in
   a given species, based on a List of DMzeroedState objects. This needs to
   be done only when simulating data conditional on site presence in a species
   -- this is automatically accounted for in the stochastic traceback
   recurrence under normal sampling. */
void dms_condition_transitions(DMotifPhyloHmm *dm, List *zeroed_states) {
  int i, s1, s2;
  double norm, rowsum, trans;
  DMzeroedState *z;
  HMM *hmm = dm->phmm->hmm;
  
  for (i = 0; i < lst_size(zeroed_states); i++) {
    z = lst_get_ptr(zeroed_states, i);
    if (z->do_row == 0) {
      continue;
    } else {
      s2 = z->state;
      for (s1 = 0; s1 < hmm->nstates; s1++) {
	mm_set(hmm->transition_matrix, s1, s2, 0);
	mm_set(hmm->transition_matrix, s2, s1, 0);
      }
      vec_set(hmm->begin_transitions, s2, 0);
    }
  }
  /* Renormalize the transitions matrix */
  mm_renormalize(hmm->transition_matrix);

  /* Renormalize the begin transitions vector */
  rowsum = norm = 0;
  for (i = 0; i < hmm->nstates; i++)
    rowsum += vec_get(hmm->begin_transitions, i);
  norm = 1 / rowsum;
  for (i = 0; i < hmm->nstates; i++) {
    trans = vec_get(hmm->begin_transitions, i);
    trans /= norm;
    vec_set(hmm->begin_transitions, i, trans);
  }

  hmm_reset(hmm);
}

/* Generate a MSA using a DMotifPhyloHmm and a set of tree models,
   requiring gain/loss elements to contain at least one substitution on 
   the branch where the gain/loss has occured. Based on tm_generate_msa. */
MSA *dm_generate_msa(int ncolumns, 
                     DMotifPhyloHmm *dm,
                     TreeModel **classmods, 
                     int *labels, /* if non-NULL, will be used to
				     record state (model) responsible
				     for generating each site; pass
				     NULL if hmm is NULL */
		     int keep_ancestral,
		     int fixed_path
                     ) {

  int i, j, class, nextclass, nseqs, col, ntreenodes, idx, nclasses, has_subs;
  char *newchar, **names, **seqs, **motif, p, c;
  MSA *msa;
  Stack *stack;
  HMM *hmm = dm->phmm->hmm;
  TreeModel *mod, *nextmod;
  TreeNode *n;
/*   TreeNode *n2; */
  
  nclasses = hmm->nstates;
  
  /* obtain number of sequences from tree models; ensure all have same
     number */
  ntreenodes = classmods[0]->tree->nnodes; 
  
  stack = stk_new_ptr(ntreenodes);
  nseqs = -1;
  for (i = 0; i < nclasses; i++) {
    /* count leaves in tree */
    int num = (classmods[i]->tree->nnodes + 1) / 2;

    assert(classmods[i]->nratecats == 1); /* assuming no rate variation */

    if (nseqs == -1) 
      nseqs = num;
    else if (nseqs != num) 
      die("ERROR in tm_generate_msa: model #%d has %d taxa, while a previous model had %d taxa.\n", i+1, num, nseqs);
  }

  /* create new MSA */
  names = (char**)smalloc(nseqs * sizeof(char*));
  seqs = (char**)smalloc(nseqs * sizeof(char*));
  for (i = 0; i < nseqs; i++) {
    seqs[i] = (char*)smalloc((ncolumns + 1) * sizeof(char));
    seqs[i][ncolumns] = '\0';
    names[i] = (char*)smalloc(STR_MED_LEN * sizeof(char));
  }
  msa = msa_new(seqs, names, nseqs, ncolumns, 
                classmods[0]->rate_matrix->states);

  /* build sequence idx maps for each model */
  for (i = 0; i < hmm->nstates; i++) {
    classmods[i]->msa_seq_idx = smalloc(classmods[0]->tree->nnodes * 
					sizeof(int));
    
    for (j = 0, idx = 0; j < classmods[i]->tree->nnodes; j++) {
      n = lst_get_ptr(classmods[i]->tree->nodes, j);
      if (n->lchild == NULL && n->rchild == NULL) {
	classmods[i]->msa_seq_idx[j] = idx;
	names[idx] = strdup(n->name);
	idx++;
      }
      else classmods[i]->msa_seq_idx[j] = -1;
    }
  }
  
  motif = (char**)smalloc(classmods[0]->tree->nnodes * sizeof(char*));
  for (i = 0; i < classmods[0]->tree->nnodes; i++)
    motif[i] = (char*)smalloc(dm->m->width * sizeof(char));
  
  /* generate sequences, column by column */
  if (fixed_path == TRUE) {
    class = labels[0];				\
  } else {
    if (hmm != NULL && hmm->begin_transitions != NULL)
      class = draw_index(hmm->begin_transitions->data, hmm->nstates);
    else
      class = 0;
  }
  
  newchar = (char*)smalloc(ntreenodes * sizeof(char));
  
  for (col = 0; col < ncolumns; /* incrementing done inside loop */) {
    mod = classmods[class];
/*     fprintf(stderr, "col %d, class %d\n", col, class); */
    
    if ( (dm->state_to_event[class] == BIRTH || 
	  dm->state_to_event[class] == DEATH) &&
	 dm->state_to_motifpos[class] != -1) { 
      /* Motif gain or loss -- must make sure motif window has
	 at least one substitution on the gain/loss branch. */
      n = lst_get_ptr(classmods[class]->tree->nodes, 
		      dm->state_to_branch[class]);
      has_subs = FALSE;

      j = 0;
      while (has_subs == FALSE) {
	/* Generate a candidate motif alignment */
	nextclass = class;
	for (i = 0; i < dm->m->width; i++) {
/* 	  fprintf(stderr, "%d ", nextclass); */

	  nextmod = classmods[nextclass];
	  dm_sample_char_col(motif, nextmod, newchar, nextclass, i, TRUE);
	  if (fixed_path == TRUE && (col + i + 1) < ncolumns) {
	    nextclass = labels[(col + i + 1)];
	  } else {
	    nextclass = mm_sample_state(hmm->transition_matrix, nextclass);
	  }
	}
	/* Check that the gain/loss branch has substitution(s) */
	for (i = 0; i < dm->m->width; i++) {
	  p = motif[n->parent->id][i];
	  c = motif[n->id][i];
	  if (p != c) {
	    has_subs = TRUE;
	    break;
	  }
	}
/*         fprintf(stderr, "Motif attempt %d (branch %d, parent %d):\n", j++, n->id, n->parent->id); */
/* 	for (i = 0; i < classmods[class]->tree->nnodes; i++) { */
/* 	  n2 = lst_get_ptr(classmods[class]->tree->nodes, i); */
/* 	  fprintf(stderr, "\tn->id %d: %s\n", n2->id, motif[n2->id]); */
/* 	} */
      }

      /* Once we have a motif alignment with appropriately places subs, push
	 the leaf sequences onto the MSA */
      for (i = 0; i < dm->m->width && col < ncolumns; i++) {
	for (j = 0; j < classmods[class]->tree->nnodes; j++) {
	  n = lst_get_ptr(classmods[class]->tree->nodes, j);
	  if (n->lchild == NULL) {
/* 	    fprintf(stderr, "i %d, j %d, col %d, ncols %d, n->id %d\n", i, j, col, ncolumns, n->id); */
	    msa->seqs[classmods[class]->msa_seq_idx[n->id]][col] = 
	      motif[n->id][i];
	  }
	}
	if (labels != NULL && fixed_path == FALSE)
	  labels[col] = class;
	if (fixed_path == TRUE && (col + 1) < ncolumns) {
	  class = labels[(col + 1)];
	} else {
	  class = mm_sample_state(hmm->transition_matrix, class);
	}
	col++;
      }
      
    } else {
      dm_sample_char_col(msa->seqs, mod, newchar, class, col,
			 keep_ancestral);

      if (labels != NULL && fixed_path == FALSE)
	labels[col] = class;
      if (fixed_path == TRUE && (col + 1) < ncolumns) {
	class = labels[(col + 1)];
      } else {
	class = mm_sample_state(hmm->transition_matrix, class);
      }
      col++;
    }
/*     fprintf(stderr, "%d ", labels[col-1]); */
  }
/*   fprintf(stderr, "\n"); */

  for (i = 0; i < classmods[0]->tree->nnodes; i++)
    free(motif[i]);
  free(motif);

  return msa;
}

/* Sample a single column of characters when simulating an MSA */
void dm_sample_char_col(char **seqs, TreeModel *mod, char *newchar, 
			int class, int col, int keep_ancestral) {
  int i;
  List *traversal;
  TreeNode *n, *l, *r;
  MarkovMatrix *lsubst_mat, *rsubst_mat;

  /* Check for substitution matrices */
  if (mod->P[0][0] == NULL)
    tm_set_subst_matrices(mod);
    
  /* Sample a character at the root */
  if (mod->alt_subst_mods != NULL &&
      mod->alt_subst_mods[mod->tree->id] != NULL) {
    newchar[mod->tree->id] = 
      mm_sample_backgd(mod->alt_subst_mods[mod->tree->id]->rate_matrix->states, 
		       mod->alt_subst_mods[mod->tree->id]->backgd_freqs);
  } else {
    newchar[mod->tree->id] = 
      mm_sample_backgd(mod->rate_matrix->states, 
		       mod->backgd_freqs);
  }
/*   fprintf(stderr, "n->id %d, newchar[root] %c\n", mod->tree->id, newchar[mod->tree->id]); */

  traversal = tr_preorder(mod->tree);  
  for (i = 0; i < lst_size(traversal); i++) {
    n = lst_get_ptr(traversal, i);
    l = n->lchild;
    r = n->rchild;
    assert ((l == NULL && r == NULL) || (l != NULL && r != NULL));
    
    if (l != NULL) { /* internal or root node -- draw chars for children */
      lsubst_mat = mod->P[l->id][0];
      rsubst_mat = mod->P[r->id][0];
      newchar[l->id] = mm_sample_char(lsubst_mat, newchar[n->id]);
      newchar[r->id] = mm_sample_char(rsubst_mat, newchar[n->id]);
/*       fprintf(stderr, "i %d, n->id %d, newchar[n] %c, newchar[lchild] %c, newchar[rchild] %c\n", i, n->id, newchar[n->id], newchar[l->id], newchar[r->id]); */
    } else { /* Leaf node */
      if (!keep_ancestral) /* not storing ancestral sequences */
	seqs[mod->msa_seq_idx[n->id]][col] = newchar[n->id];
    }
  }
  if (keep_ancestral) {
    for (i = 0; i < lst_size(traversal); i++) {
      n = lst_get_ptr(traversal, i);
      seqs[n->id][col] = newchar[n->id];
    }
  }
}
  
/* Set Halpern-Bruno model alt_subst_mod objects for motif branches within a
   TreeModel object. If nodelist size is equal to the number of nodes in the
   source tree, conserved motif is assumed and alt_subst_mods will not be
   used -- instead, rate_matrix, background_freqs and subst_mod in the
   original model will be reset and originals freed. */
void dm_set_motif_branches_hb(TreeModel *tm, Vector *motif_freqs, 
			      List *nodelist) {
  int i, j;
  double rab, pia, pib, pab, pba, rowsum;
  MarkovMatrix *hb_matrix;
  TreeNode *n;

  if (tm->subst_mod != REV)
    die("ERROR: Halpern-Bruno motif model requires use of an REV background model!\n");

  /* First compute the Halpern-Bruno substitution matrix based on the equations
     in Halpern and Bruno, 1998. MBE. 15(7):910-917. Note that the scaling
     factor, k, used in Halpern-Bruno is not needed when using REV as the 
     background model. */
  hb_matrix = mm_new(motif_freqs->size, tm->rate_matrix->states, CONTINUOUS);
  for (i = 0; i < motif_freqs->size; i++) {
    rowsum = 0;
    for (j = 0; j < motif_freqs->size; j++) {
      if (i == j)
	continue;
      else { /* off-diagonal element -- compute rab and add to the row sum */
	pia = vec_get(motif_freqs, i); /* Starting nucleotide frequency */ 
	pib = vec_get(motif_freqs, j); /* Ending nucleotide frequency */
	pab = mm_get(tm->rate_matrix, i, j); /* REV rate param for a->b mut. */
	pba = mm_get(tm->rate_matrix, j, i); /* REV rate param for b->a mut. */

	if ((pib * pba) / (pia * pab) != 1) {
	  rab = pab * ( (log((pib * pba) / (pia * pab))) / 
			(1 - ((pia * pab) / (pib * pba)))
			);
	} else {
	  rab = pab;
	}

	mm_set(hb_matrix, i, j, rab);
	rowsum += rab;
/* 	fprintf(stderr, "i %d, j %d, pia %f, pib %f, pab %f, pba %f, rab %f\n", */
/* 		i, j, pia, pib, pab, pba, rab); */
      }
    }
    mm_set(hb_matrix, i, i, -rowsum);
/*     for (j = 0; j < motif_freqs->size; j++) */
/*       fprintf(stderr, "%f\t", mm_get(hb_matrix, i, j)); */
/*     fprintf(stderr, "\n"); */
  }
/*   fprintf(stderr, "\n"); */
/*   mm_pretty_print(stderr, hb_matrix); */
  
  /* Now set up the alternate subst. mods for the motif branches in the source
     model, according to the nodelist. Check to see if the entire tree is
     contained within the nodes list -- this indicates the entire tree is under
     selection as a motif (i.e., conserved motif) and we can simply substitute
     the freqs and rate_matrix within the entire mod instead of using the
     alt_subst_mods. */
/*   if (lst_size(nodelist) == tm->tree->nnodes) { */
/*     mm_free(tm->rate_matrix); */
/*     vec_free(tm->backgd_freqs); */
/*     tm->rate_matrix = mm_create_copy(hb_matrix); */
/*     tm->backgd_freqs = vec_create_copy(motif_freqs); */
/*     tm->subst_mod = HB; */
/*   } else { */
    tm->alt_subst_mods = (AltSubstMod**)smalloc(tm->tree->nnodes
						* sizeof(AltSubstMod*));

    for (i = 0; i < tm->tree->nnodes; i++) tm->alt_subst_mods[i] = NULL;
    
    for (i = 0; i < lst_size(nodelist); i++) {
      n = lst_get_ptr(nodelist, i);
      tm->alt_subst_mods[n->id] =
	tm_new_alt_subst_mod(HB, vec_create_copy(motif_freqs),
			     mm_create_copy(hb_matrix));
/*       fprintf(stderr, "mod[%d]:\n", n->id); */
/*       fprintf(stderr, "\tbg: "); */
/*       vec_print(tm->alt_subst_mods[n->id]->backgd_freqs, stderr); */
/*       fprintf(stderr, "rate_matrix:\n"); */
/*       mm_pretty_print(stderr, tm->alt_subst_mods[n->id]->rate_matrix); */
/*       fprintf(stderr, "\n"); */
    }
/*   } */
  /* Free up the source hb_matrix and we're done */
  mm_free(hb_matrix);
}

dmevent_t dm_get_event_type(char *ename) {
  if (strcmp(ename, "NEUT") == 0)
    return NEUT;
  else if (strcmp(ename, "CONS") == 0)
    return CONS;
  else if (strcmp(ename, "BIRTH") == 0)
    return BIRTH;
  else if (strcmp(ename, "DEATH") == 0)
    return DEATH;
  else
    return -1;
}

/* Produce an indel history for the reverse strand */
IndelHistory *dms_reverse_ih(IndelHistory *ih) {
  int i, j;
  CompactIndelHistory *cih_f, *cih_r;
  IndelHistory *retval;
  Indel *indel_f, *indel_r;

  cih_f = ih_compact(ih);

  cih_r = ih_new_compact(tr_create_copy(cih_f->tree), cih_f->ncols);
  for (i = 0; i < cih_f->tree->nnodes; i++) {
    for (j = 0; j < lst_size(cih_f->indels[i]); j++) {
      indel_f = lst_get_ptr(cih_f->indels[i], j);
      indel_r = (Indel*)smalloc(sizeof(Indel));
      indel_r->type = indel_f->type;
      indel_r->len = indel_f->len;
      indel_r->start = (cih_f->ncols -
			((indel_f->start + indel_f->len - 1) + 1));
      lst_push_ptr(cih_r->indels[i], indel_r);
    }
  }
  retval = ih_expand(cih_r);
  ih_free_compact(cih_f);
  ih_free_compact(cih_r);
  return retval;
}

/* Fold non-identical motif features from the second gff set into the first
   gff set, which will then contain the union of both sets */
void dms_combine_gffs(GFF_Set *target_gff, GFF_Set *query_gff) {
  int i, j, found;
  GFF_Feature *tf, *qf;
  GFF_Set *unique = gff_new_set();

  for (i = 0; i < lst_size(query_gff->features); i++) {
    qf = lst_get_ptr(query_gff->features, i);
    found = 0;
    for (j = 0; j < lst_size(target_gff->features); j++) {
      tf = lst_get_ptr(target_gff->features, j);
      if (qf->seqname == tf->seqname &&
	  qf->feature == tf->feature &&
	  qf->start == tf->start &&
	  qf->end == tf-> end &&
	  qf->strand == tf->strand) { /* Features are identical */
	found = 1;
	break;
      }
    }
    if (found == 0) /* query feature not found in target set */
      lst_push_ptr(unique->features, gff_new_feature_copy(qf));
  }
  for (i = 0; i < lst_size(unique->features); i++) {
    qf = lst_get_ptr(unique->features, i);
    lst_push_ptr(target_gff->features, gff_new_feature_copy(qf));
  }
  gff_free_set(unique);
}

/* Merge two GFF sets without comparing the contents -- assumes no overlap
   between target and query! This is valid for single thread gff_sets from
   dms_sample_path because each thread samples paths in a unique set of
   sequences, thus gff_sets from different threads should never overlap. */
void dms_merge_thread_gffs(GFF_Set *target_gff, GFF_Set *query_gff) {
  int i;
  GFF_Feature *f;
  for (i = 0; i < lst_size(query_gff->features); i++) {
    f = lst_get_ptr(query_gff->features, i);
    lst_push_ptr(target_gff->features, gff_new_feature_copy(f));
  }
}
