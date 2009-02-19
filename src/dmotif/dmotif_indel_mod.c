/* $Id: dmotif_indel_mod.c,v 1.5 2009-02-19 19:38:48 agd27 Exp $
   Written by Adam Siepel, 2005
   Copyright 2005, Adam Siepel, University of California */

/* simple model of insertions and deletions, assumes given indel history */

#include <misc.h>
#include <assert.h>
#include <dmotif_indel_mod.h>
#include <numerical_opt.h>

DMotifBranchIndelModel *dmih_new_branch(double alpha, double beta, double tau,
					double epsilon, double t, int motif) {
  DMotifBranchIndelModel *bim =
    (DMotifBranchIndelModel*)smalloc(sizeof(DMotifBranchIndelModel));
  bim->probs = mm_new(NINDEL_STATES, NULL, DISCRETE);
  bim->log_probs = mat_new(NINDEL_STATES, NINDEL_STATES);
  bim->beg_probs = vec_new(NINDEL_STATES);
  bim->beg_log_probs = vec_new(NINDEL_STATES);
  dmih_set_branch(bim, alpha, beta, tau, epsilon, t, motif);
  return bim;
}


/* To Do: This needs to be retooled to pay attention to motif vs. bg states */
DMotifIndelModel *dmih_new_all(double alpha, double beta, double tau, 
			       double epsilon, TreeNode *tree) {
  int i, mnode;
  /* Figure out how mnode should be used!! */
  DMotifIndelModel *im = (DMotifIndelModel*)smalloc(sizeof(DMotifIndelModel));
  mnode = 0;
  im->training_lnl = 0;
  im->branch_mods = smalloc(tree->nnodes * sizeof(void*));
  for (i = 0; i < tree->nnodes; i++) {
    TreeNode *n = lst_get_ptr(tree->nodes, i);
    if (n == tree) { im->branch_mods[i] = NULL; continue; }
    im->branch_mods[i] = dmih_new_branch(alpha, beta, tau, epsilon, 
					 n->dparent, mnode);
  }
  dmih_set_all(im, alpha, beta, tau, epsilon, tree);
  return im;
}

DMotifIndelModel *dmih_new(double *alpha, double *beta, double *tau, 
			   double *epsilon, TreeNode *tree, int e, List *l,
			   int motif) {
  int i, j, mnode;
  DMotifIndelModel *im = (DMotifIndelModel*)smalloc(sizeof(DMotifIndelModel));
  im->training_lnl = 0;
  im->branch_mods = (DMotifBranchIndelModel**)smalloc(tree->nnodes 
						      * sizeof(void*));
  
  for (i = 0; i < tree->nnodes; i++) {
    TreeNode *n = lst_get_ptr(tree->nodes, i);
    mnode = 0;
    if (n == tree) {
      im->branch_mods[i] = NULL; 
      continue;
    }
    if (motif) { /* See if we're at a motif node */
      if (e == CONS || e == NEUT) { /* All nodes are motif nodes */
	mnode = 1;
      } else { /* Gain or Loss -- search the motif nodes list (l) for the
		  current node id */
	for (j = 0; j < lst_size(l); j++) {
	  if (i == ((TreeNode*)(lst_get_ptr(l, j)))->id) {
	    mnode = 1;
	    break;
	  }
	}
      }
    }    
    im->branch_mods[i] = dmih_new_branch(alpha[i], beta[i], tau[i], epsilon[i],
					 n->dparent, mnode);
  }
  im->tree = tree;
  im->alpha = im->beta = im->tau = im->epsilon = -1;
  /* This call seems to be redundant when the above two lines are set here! */
/*   dmih_set(im, alpha, beta, tau, epsilon, tree, e, l, motif); */
  return im;
}

void dmih_free_branch(DMotifBranchIndelModel *bim) {
  mm_free(bim->probs);
  mat_free(bim->log_probs);
  vec_free(bim->beg_probs);
  vec_free(bim->beg_log_probs);
  free(bim);
}

 void dmih_free(DMotifIndelModel *im) {
  int i;
  for (i = 0; i < im->tree->nnodes; i++)
    if (im->branch_mods[i] != NULL)
      dmih_free_branch(im->branch_mods[i]);
  free(im->branch_mods);
  free(im);
}

void dmih_set_branch(DMotifBranchIndelModel *bim, double alpha, double beta, 
		     double tau, double epsilon, double t, int motif) {
  int i, j;

  bim->alpha = alpha;
  bim->beta = beta;
  bim->tau = tau;
  bim->epsilon = epsilon;
  bim->t = t;

  if (motif) { /* Motif state -- use epsilon as the indel param wherever we
		  have a motif (supertree for losses, subtree for gains,
		  everywhere for cons or neutral) */
    /* transition probabilities */
    mm_set(bim->probs, CHILDINS, CHILDINS, 1 - tau - epsilon * t);
    mm_set(bim->probs, CHILDINS, CHILDDEL, epsilon * t);
    mm_set(bim->probs, CHILDINS, MATCH, tau);
    
    mm_set(bim->probs, CHILDDEL, CHILDINS, epsilon * t);
    mm_set(bim->probs, CHILDDEL, CHILDDEL, 1 - tau - epsilon * t);
    mm_set(bim->probs, CHILDDEL, MATCH, tau);
    
    mm_set(bim->probs, MATCH, CHILDINS, epsilon * t);
    mm_set(bim->probs, MATCH, CHILDDEL, epsilon * t);
    mm_set(bim->probs, MATCH, MATCH, 1 - 2 * epsilon * t);
    
    /* set begin probabilities to stationary distribution */
    vec_set(bim->beg_probs, CHILDINS, 
	    epsilon * t / (epsilon * t + tau));
    vec_set(bim->beg_probs, CHILDDEL, 
	    epsilon * t / (epsilon * t + tau));
    vec_set(bim->beg_probs, MATCH, 
	    tau / (epsilon * t + tau));
  } else { /* Background state -- use default dless model */
    /* transition probabilities */
    mm_set(bim->probs, CHILDINS, CHILDINS, 1 - tau - beta * t);
    mm_set(bim->probs, CHILDINS, CHILDDEL, beta * t);
    mm_set(bim->probs, CHILDINS, MATCH, tau);
    
    mm_set(bim->probs, CHILDDEL, CHILDINS, alpha * t);
    mm_set(bim->probs, CHILDDEL, CHILDDEL, 1 - tau - alpha * t);
    mm_set(bim->probs, CHILDDEL, MATCH, tau);
    
    mm_set(bim->probs, MATCH, CHILDINS, alpha * t);
    mm_set(bim->probs, MATCH, CHILDDEL, beta * t);
    mm_set(bim->probs, MATCH, MATCH, 1 - alpha * t - beta * t);
    
    /* set begin probabilities to stationary distribution */
    vec_set(bim->beg_probs, CHILDINS, 
	    alpha * t / (alpha * t + beta * t + tau));
    vec_set(bim->beg_probs, CHILDDEL, 
	    beta * t / (alpha * t + beta * t + tau));
    vec_set(bim->beg_probs, MATCH, 
	    tau / (alpha * t + beta * t + tau));
  }

  /* set log probs; also check that everything is a probability */
  for (i = 0; i < NINDEL_STATES; i++) {
    double rowsum = 0;
    for (j = 0; j < NINDEL_STATES; j++) {
      double prob = mm_get(bim->probs, i, j);
      rowsum += prob;
      if (prob < 0 || prob > 1) 
        die("ERROR: invalid indel probability.  Alpha, beta, tau and/or epsilon\nare probably too large given branch lengths.\n");
      mat_set(bim->log_probs, i, j, log2(prob));
    }
    assert(fabs(1-rowsum) < 1e-6);
    vec_set(bim->beg_log_probs, i, log2(vec_get(bim->beg_probs, i)));
  }
}

/* To Do: This needs to be retooled to pay attention to motif vs. bg states */
void dmih_set_all(DMotifIndelModel *im, double alpha, double beta, double tau, 
		  double epsilon, TreeNode *tree) {
  int i, mnode = 0; /* To do check value */
  im->alpha = alpha;
  im->beta = beta;
  im->tau = tau;
  im->tree = tree;

  for (i = 0; i < tree->nnodes; i++) {
    TreeNode *n = lst_get_ptr(tree->nodes, i);
    if (n != tree)
      dmih_set_branch(im->branch_mods[i], alpha, beta, tau, epsilon,
		      n->dparent, mnode);
  }
}

void dmih_set(DMotifIndelModel *im, double *alpha, double *beta, double *tau, 
	      double *epsilon, TreeNode *tree, int e, List *l, int motif) {
  int i, j, mnode;
  im->alpha = im->beta = im->tau = im->epsilon = -1;
  im->tree = tree;

  for (i = 0; i < tree->nnodes; i++) {
    TreeNode *n = lst_get_ptr(tree->nodes, i);
    mnode = 0;
    
    if (motif) { /* See if we're at a motif node */
      if (e == CONS || e == NEUT) { /* All nodes are motif nodes */
	mnode = 1;
      } else { /* Gain or Loss -- search the motif nodes list (l) for the
		  current node id */
	for (j = 0; j < lst_size(l); j++) {
	  if (i == ((TreeNode*)(lst_get_ptr(l, j)))->id) {
	    mnode = 1;
	    break;
	  }
	}
      }
    }

    if (n != tree)
      dmih_set_branch(im->branch_mods[i], alpha[i], beta[i], tau[i],
		      epsilon[i], n->dparent, mnode);
  }
}

static inline
col_type get_col_type(IndelHistory *ih, int child_id, int parent_id, int col) {
  if (ih->indel_strings[parent_id][col] == BASE && 
      ih->indel_strings[child_id][col] == BASE)
    return MATCH;
  else if (ih->indel_strings[parent_id][col] == INS && 
           ih->indel_strings[child_id][col] == BASE)
    return CHILDINS;
  else if (ih->indel_strings[parent_id][col] == BASE && 
           ih->indel_strings[child_id][col] == DEL)
    return CHILDDEL;
  else if (ih->indel_strings[parent_id][col] == INS && 
           ih->indel_strings[child_id][col] == INS)
    return SKIP;
  else if (ih->indel_strings[parent_id][col] == DEL && 
           ih->indel_strings[child_id][col] == DEL)
    return SKIP;
  else 
    return ERROR;
}

double dmih_branch_column_logl(IndelHistory *ih, DMotifBranchIndelModel *bim, 
			       int child, double *col_logl) {
  int i;
  col_type this_type, last_type = SKIP;
  double logl = 0;
  int parent_id = ((TreeNode*)lst_get_ptr(ih->tree->nodes, child))->parent->id;

/*   for (i = 0; i < ih->ncols; i++) { */
/*     last_type = get_col_type(ih, child, parent_id, i); */
/*     if (last_type == SKIP)  */
/*       col_logl[i] = 0; */
/*     else { */
/*       col_logl[i] = vec_get(bim->beg_log_probs, last_type); */
/*       break; */
/*     } */
/*   } */
/*   logl = col_logl[i++]; */
/*   /\* Note that in DLESS, i is NOT advanced before restarting the loop -- this */
/*      seems to result in the standard transition prob. replacing the begin prob. */
/*      in the col_logl vector and a double-counted transition. *\/ */
/*   for (; i < ih->ncols; i++) { */
/*     this_type = get_col_type(ih, child, parent_id, i); */

/*     assert(this_type != ERROR); */

/*     if (this_type == SKIP) { */
/*       col_logl[i] = 0; */
/*       continue; */
/*     } */

/*     col_logl[i] = bim->log_probs->data[last_type][this_type]; */
/*     logl += col_logl[i]; */
/*     last_type = this_type; */
/*   } */

  for (i = 0; i < ih->ncols; i++) {
    this_type = get_col_type(ih, child, parent_id, i);
    
    assert(this_type != ERROR);
    
    if (this_type == SKIP) {
      col_logl[i] = 0;
      continue;
    }
    if (last_type == SKIP)
      col_logl[i] = vec_get(bim->beg_log_probs, this_type);
    else
      col_logl[i] = bim->log_probs->data[last_type][this_type];
    
    logl += col_logl[i];
    last_type = this_type;
  }
  
  return logl;
}

/* returns total log likelihood */
double dmih_column_logl(IndelHistory *ih, DMotifIndelModel *im, 
			double *col_logl) {
  int i, j;
  double logl = 0;
  double *branch_col_logl = smalloc(ih->ncols * sizeof(double));
  TreeNode *n;

  for (i = 0; i < ih->ncols; i++) col_logl[i] = 0;
  for (i = 0; i < im->tree->nnodes; i++) {
    n = lst_get_ptr(im->tree->nodes, i);
    if (n == im->tree) continue;
    logl += dmih_branch_column_logl(ih, im->branch_mods[n->id], i, 
				    branch_col_logl);
    for (j = 0; j < ih->ncols; j++) col_logl[j] += branch_col_logl[j];
  }
  free(branch_col_logl);

  /* NOTE: it looks like we're combining joint probabilities over
     indel strings here, when we should be using conditional
     probabilities, but what we're really combining are the
     probabilities of indel *events* along branches, which are
     independent given the indel history */

  return logl;
}

double dmih_likelihood(DMotifIndelModel *im, DMotifIndelSuffStats *iss) {
  double logl = 0;
  int i, j, k;
  for (i = 0; i < im->tree->nnodes; i++) {
    TreeNode *n = lst_get_ptr(im->tree->nodes, i);
    if (n == im->tree) continue;
    for (j = 0; j < NINDEL_STATES; j++) {
      logl += iss->branch_counts[i]->beg_counts->data[j] * 
        im->branch_mods[i]->beg_log_probs->data[j];
      for (k = 0; k < NINDEL_STATES; k++)
        logl += iss->branch_counts[i]->trans_counts->data[j][k] * 
          im->branch_mods[i]->log_probs->data[j][k];
    }
  }
  
  return logl;
}

DMotifBranchIndelSuffStats *dmih_suff_stats_branch(IndelHistory *ih, int child_id) {
  int i, j;
  char c;
  col_type this_type, last_type = SKIP;
  int parent_id = ((TreeNode*)lst_get_ptr(ih->tree->nodes, child_id))->parent->id;
  DMotifBranchIndelSuffStats *ss = smalloc(sizeof(DMotifBranchIndelSuffStats));
  ss->trans_counts = mat_new(NINDEL_STATES, NINDEL_STATES);
  ss->beg_counts = vec_new(NINDEL_STATES);
  mat_zero(ss->trans_counts);
  vec_zero(ss->beg_counts);

  for (i = 0; last_type == SKIP && i < ih->ncols; i++)
    last_type = get_col_type(ih, child_id, parent_id, i);
  ss->beg_counts->data[last_type]++;
  for (; i < ih->ncols; i++) {
    this_type = get_col_type(ih, child_id, parent_id, i);

    if (this_type == ERROR) {
      fprintf(stderr, "ERROR at column %d of indel history:\n", i);
      for (j = 0; j < ih->tree->nnodes; j++) {
        if (ih->indel_strings[j][i] == BASE)
          c = 'b';
        else if (ih->indel_strings[j][i] == INS)
          c = '^';
        else
          c = '.';
        fprintf(stderr, "%25s %c\n", 
                ((TreeNode*)lst_get_ptr(ih->tree->nodes, j))->name, c);
      }
      assert(0);
    }

    else if (this_type == SKIP) continue;

    ss->trans_counts->data[last_type][this_type]++;
    last_type = this_type;
  }

  return ss;  
}

DMotifIndelSuffStats *dmih_suff_stats(IndelHistory *ih) {
  int i;
  DMotifIndelSuffStats *iss = smalloc(sizeof(DMotifIndelSuffStats));
  iss->tree = ih->tree;
  iss->branch_counts = smalloc(ih->tree->nnodes * sizeof(void*));
  for (i = 0; i < ih->tree->nnodes; i++) {
    TreeNode *n = lst_get_ptr(ih->tree->nodes, i);
    if (n == ih->tree) {
      iss->branch_counts[i] = NULL;
      continue;
    }
    iss->branch_counts[i] = dmih_suff_stats_branch(ih, n->id);
  }
  return iss;
}

void dmih_free_suff_stats(DMotifIndelSuffStats *iss) {
  int i;
  for (i = 0; i < iss->tree->nnodes; i++) {
    if (iss->branch_counts[i] != NULL) {
      mat_free(iss->branch_counts[i]->trans_counts);
      vec_free(iss->branch_counts[i]->beg_counts);
      free(iss->branch_counts[i]);
    }
  }
  free(iss->branch_counts);
  free(iss);
}

/* double dmih_simulate_history(DMotifIndelModel *tim, int ncols) { */
/* } */

/* package of data used below */
struct likelihood_data {
  DMotifIndelModel *im;
  DMotifIndelSuffStats *ss;
};

/* wrapper for likelihood function, for use in numerical maximization */
double dmih_likelihood_wrapper(Vector *params, void *data) {
  struct likelihood_data *d = data;
  dmih_set_all(d->im, vec_get(params, 0), vec_get(params, 1), 
	       vec_get(params, 2), vec_get(params,3), d->im->tree);
  return -dmih_likelihood(d->im, d->ss);
}

/* computes gradient of (negative) likelihood, for use in numerical
   maximization */
void dmih_likelihood_gradient(Vector *grad, Vector *params, void *data,
                            Vector *lb, Vector *ub) {
  struct likelihood_data *d = data;
  double alpha = vec_get(params, 0), beta = vec_get(params, 1),
    tau = vec_get(params, 2);
  double alpha_deriv = 0, beta_deriv = 0, tau_deriv = 0;
  int i;

  for (i = 0; i < d->im->tree->nnodes; i++) {
    Matrix *c;
    Vector *begc;
    double sumc, denom;
    TreeNode *n = lst_get_ptr(d->im->tree->nodes, i);

    if (n == d->im->tree) continue;

    c = d->ss->branch_counts[i]->trans_counts;

    alpha_deriv +=
      (c->data[MATCH][CHILDINS] + c->data[CHILDDEL][CHILDINS]) / alpha -
      (c->data[MATCH][MATCH] * n->dparent) /
      (1 - alpha * n->dparent - beta * n->dparent) -
      (c->data[CHILDDEL][CHILDDEL] * n->dparent) /
      (1 - tau - alpha * n->dparent);

    beta_deriv +=
      (c->data[MATCH][CHILDDEL] + c->data[CHILDINS][CHILDDEL]) / beta -
      (c->data[MATCH][MATCH] * n->dparent) /
      (1 - alpha * n->dparent - beta * n->dparent) -
      (c->data[CHILDINS][CHILDINS] * n->dparent) /
      (1 - tau - beta * n->dparent);

    tau_deriv +=
      (c->data[CHILDINS][MATCH] + c->data[CHILDDEL][MATCH]) / tau -
      c->data[CHILDINS][CHILDINS] / (1 - tau - beta * n->dparent) -
      c->data[CHILDDEL][CHILDDEL] / (1 - tau - alpha * n->dparent);

    /* also consider begins */
    begc = d->ss->branch_counts[i]->beg_counts;
    sumc = begc->data[MATCH] + begc->data[CHILDINS] + begc->data[CHILDDEL];
    denom = tau + alpha * n->dparent + beta * n->dparent;
    alpha_deriv += begc->data[CHILDINS] / alpha - sumc * n->dparent / denom;
    beta_deriv += begc->data[CHILDDEL] / beta - sumc * n->dparent / denom;
    tau_deriv += begc->data[MATCH] / tau - sumc / denom;
  }

  vec_set(grad, 0, alpha_deriv);
  vec_set(grad, 1, beta_deriv);
  vec_set(grad, 2, tau_deriv);
  vec_scale(grad, -1/log(2));    /* because min rather than max and
                                    because im_likelihood uses log
                                    base 2 */
}

/* estimate alpha, beta, and tau from an indel history by
   maximum likelihood */
void dmih_estimate(DMotifIndelModel *im, IndelHistory *ih, 
		   DMotifIndelSuffStats *ss, FILE *logf) {
  Vector *params = vec_new(3), *lb = vec_new(3), *ub = vec_new(3);
  struct likelihood_data *d = smalloc(sizeof(struct likelihood_data));
  double neglogl;

  d->im = im;
  d->ss = ss;
  vec_set(params, 0, im->alpha);
  vec_set(params, 1, im->beta);
  vec_set(params, 2, im->tau);
  vec_set_all(lb, 1e-6);
  vec_set_all(ub, 0.5);

  opt_bfgs(dmih_likelihood_wrapper, params, d, &neglogl, lb, ub, logf,
           dmih_likelihood_gradient, OPT_HIGH_PREC, NULL);

  dmih_set_all(im, vec_get(params, 0), vec_get(params, 1),
             vec_get(params, 2), vec_get(params, 3), im->tree);
  im->training_lnl = -neglogl * log(2);
  
  vec_free(params);
  vec_free(lb);
  dmih_free_suff_stats(d->ss);
  free(d);
}

/* collect sufficient stats for a branch, considering only sites in
   the specified category */
DMotifBranchIndelSuffStats *dmih_suff_stats_branch_cat(IndelHistory *ih, 
						 int child_id, int *categories,
						 int do_cat) {
  int i, j;
  char c;
  col_type this_type, last_type;
  int parent_id = ((TreeNode*)lst_get_ptr(ih->tree->nodes, child_id))->parent->id;
  DMotifBranchIndelSuffStats *ss = smalloc(sizeof(DMotifBranchIndelSuffStats));
  ss->trans_counts = mat_new(NINDEL_STATES, NINDEL_STATES);
  ss->beg_counts = vec_new(NINDEL_STATES);
  mat_zero(ss->trans_counts);
  vec_zero(ss->beg_counts);

  /* scan to first non-SKIP in category of interest */
  for (i = 0; i < ih->ncols; i++) {
    if (categories[i] != do_cat) continue;
    if ((this_type = get_col_type(ih, child_id, parent_id, i)) != SKIP)
      break;
  }
  if (i == ih->ncols) return ss;

  ss->beg_counts->data[this_type]++;
  last_type = this_type;
  for (; i < ih->ncols; i++) {
    this_type = get_col_type(ih, child_id, parent_id, i);

    if (this_type == ERROR) {
      fprintf(stderr, "ERROR at column %d of indel history:\n", i);
      for (j = 0; j < ih->tree->nnodes; j++) {
        if (ih->indel_strings[j][i] == BASE)
          c = 'b';
        else if (ih->indel_strings[j][i] == INS)
          c = '^';
        else
          c = '.';
        fprintf(stderr, "%25s %c\n", 
                ((TreeNode*)lst_get_ptr(ih->tree->nodes, j))->name, c);
      }
      assert(0);
    }
    else if (this_type == SKIP) continue;

    if (categories[i] == do_cat) 
      ss->trans_counts->data[last_type][this_type]++;

    last_type = this_type;      /* need to set last_type even if not
                                   in category; will use if next site
                                   is in category  */
  }

  return ss;  
}

/* collect sufficient stats, considering only sites in a particular
   category */
DMotifIndelSuffStats *dmih_suff_stats_cat(IndelHistory *ih, int *categories, 
                                  int do_cat) {
  int i;
  DMotifIndelSuffStats *iss = smalloc(sizeof(DMotifIndelSuffStats));
  iss->tree = ih->tree;
  iss->branch_counts = smalloc(ih->tree->nnodes * sizeof(void*));
  for (i = 0; i < ih->tree->nnodes; i++) {
    TreeNode *n = lst_get_ptr(ih->tree->nodes, i);
    if (n == ih->tree) {
      iss->branch_counts[i] = NULL;
      continue;
    }
    iss->branch_counts[i] = dmih_suff_stats_branch_cat(ih, n->id, 
                                                     categories, do_cat);
  }
  return iss;
}

/* convert to an alignment, including sequences for ancestral nodes as
   well as leaf nodes, and with '^' characters in place of '-' for
   insertions and '.' characters in place of '-' for deletions.
   Useful for debugging */
MSA *dmih_as_alignment(IndelHistory *ih, MSA *msa) {
  int i, j, k, s, ins;
  char **seqs = smalloc(ih->tree->nnodes * sizeof(char*));
  char **names = smalloc(ih->tree->nnodes * sizeof(char*));
  List *inside, *outside;
  TreeNode *n, *n2;

  inside = lst_new_ptr(10);
  outside = lst_new_ptr(10);

  for (i = 0; i < ih->tree->nnodes; i++) {
    n = lst_get_ptr(ih->tree->nodes, i);
    names[i] = strdup(n->name);
    seqs[i] = smalloc((ih->ncols+1) * sizeof(char));
  }

  for (i = 0; i < ih->tree->nnodes; i++) {
    n = lst_get_ptr(ih->tree->nodes, i);
    /* initialize with actual bases if available or 'N's otherwise */
    if (n->lchild == NULL) {    /* leaf */
      if (msa != NULL) {
        if ((s = msa_get_seq_idx(msa, n->name)) < 0)
          die("ERROR: no match for leaf \"%s\" in alignment.\n", n->name);
      }
      for (j = 0; j < ih->ncols; j++) {
        if (ih->indel_strings[i][j] == BASE) 
          seqs[i][j] = msa == NULL ? 'N' : msa_get_char(msa, s, j);
        else {
	  if (ih->indel_strings[i][j] == INS) { /* Insertion */
	    /* Find the node below the branch where the insertion happened */
	    for (k = 0; 
		 k < ih->tree->nnodes && ih->indel_strings[k][j] != BASE;
		 k++){
	      if (k == 0 || i == ih->tree->nnodes)
		ins = -1;
	      else 
		ins = k;
	    }
	    /* Nodes in the subtree under the node receiving the insertion
	       will all have a non-insertion char while those in the supertree
	       will all have an insertion character */
	    tr_partition_nodes(ih->tree, lst_get_ptr(ih->tree->nodes, ins), 
			       inside, outside);
	    for (k = 0; k < lst_size(inside); k++) {
	      n2 = lst_get_ptr(inside, k);
	      s = msa_get_seq_idx(msa, n2->name);
	      seqs[n2->id][j] = (n2->lchild == NULL) ?
		msa_get_char(msa, s, j) : 'N';
	    }
	    for (k = 0; k < lst_size(outside); k++) {
	      n2 = lst_get_ptr(outside, k);
	      seqs[n2->id][j] = '^';
	    }
	  } else { /* Deletion */
	    seqs[i][j] = '.';
	  }
	}
      }
    }


    else {                      /* ancestor */
      for (j = 0; j < ih->ncols; j++) {
        if (ih->indel_strings[i][j] == BASE) 
          seqs[i][j] = 'N';
        else
          seqs[i][j] = ih->indel_strings[i][j] == INS ? '^' : '.';
      }
    }

    seqs[i][ih->ncols] = '\0';
  }
  lst_free(inside);
  lst_free(outside);
  return msa_new(seqs, names, ih->tree->nnodes, ih->ncols, "ACGTN-^.");
}
