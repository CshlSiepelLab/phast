/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: indel_mod.c,v 1.9 2009-02-19 22:16:16 acs Exp $ */

/* simple model of insertions and deletions, assumes given indel history */

#include <misc.h>
#include <indel_mod.h>
#include <numerical_opt.h>
#include <external_libs.h>

BranchIndelModel *im_new_branch(double alpha, double beta, double tau,
                                double t) {
  BranchIndelModel *bim = smalloc(sizeof(BranchIndelModel));
  bim->probs = mm_new(NINDEL_STATES, NULL, DISCRETE);
  bim->log_probs = mat_new(NINDEL_STATES, NINDEL_STATES);
  bim->beg_probs = vec_new(NINDEL_STATES);
  bim->beg_log_probs = vec_new(NINDEL_STATES);
  im_set_branch(bim, alpha, beta, tau, t);
  return bim;
}

IndelModel *im_new_all(double alpha, double beta, double tau, 
                       TreeNode *tree) {
  int i;
  IndelModel *im = smalloc(sizeof(IndelModel));
  im->training_lnl = 0;
  im->branch_mods = smalloc(tree->nnodes * sizeof(void*));
  for (i = 0; i < tree->nnodes; i++) {
    TreeNode *n = lst_get_ptr(tree->nodes, i);
    if (n == tree) { im->branch_mods[i] = NULL; continue; }
    im->branch_mods[i] = im_new_branch(alpha, beta, tau, n->dparent);
  }
  im_set_all(im, alpha, beta, tau, tree);
  return im;
}

IndelModel *im_new(double *alpha, double *beta, double *tau, 
                   TreeNode *tree) {
  int i;
  IndelModel *im = smalloc(sizeof(IndelModel));
  im->training_lnl = 0;
  im->branch_mods = smalloc(tree->nnodes * sizeof(void*));
  for (i = 0; i < tree->nnodes; i++) {
    TreeNode *n = lst_get_ptr(tree->nodes, i);
    if (n == tree) { im->branch_mods[i] = NULL; continue; }
    im->branch_mods[i] = im_new_branch(alpha[i], beta[i], tau[i], n->dparent);
  }
  im_set(im, alpha, beta, tau, tree);
  return im;
}

void im_free_branch(BranchIndelModel *bim) {
  mm_free(bim->probs);
  mat_free(bim->log_probs);
  vec_free(bim->beg_probs);
  vec_free(bim->beg_log_probs);
  sfree(bim);
}

void im_free(IndelModel *im) {
  int i;
  for (i = 0; i < im->tree->nnodes; i++)
    if (im->branch_mods[i] != NULL)
      im_free_branch(im->branch_mods[i]);
  sfree(im->branch_mods);
  sfree(im);
}

void im_set_branch(BranchIndelModel *bim, double alpha, 
                   double beta, double tau, double t) {
  int i, j;

  bim->alpha = alpha;
  bim->beta = beta;
  bim->tau = tau;
  bim->t = t;

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

  /* set log probs; also check that everything is a probability */
  for (i = 0; i < NINDEL_STATES; i++) {
    for (j = 0; j < NINDEL_STATES; j++) {
      double prob = mm_get(bim->probs, i, j);
      if (prob < 0 || prob > 1) 
        die("ERROR: invalid indel probability.  Alpha, beta, and/or tau are probably too\nlarge given branch lengths.\n");
      mat_set(bim->log_probs, i, j, log2(prob));
    }
    vec_set(bim->beg_log_probs, i, log2(vec_get(bim->beg_probs, i)));
  }
}

void im_set_all(IndelModel *im, double alpha, double beta, double tau, 
            TreeNode *tree) {
  int i;
  im->alpha = alpha;
  im->beta = beta;
  im->tau = tau;
  im->tree = tree;

  for (i = 0; i < tree->nnodes; i++) {
    TreeNode *n = lst_get_ptr(tree->nodes, i);
    if (n != tree)
      im_set_branch(im->branch_mods[i], alpha, beta, tau, n->dparent);
  }
}

void im_set(IndelModel *im, double *alpha, double *beta, double *tau,
            TreeNode *tree) {
  int i;
  im->alpha = im->beta = im->tau = -1;
  im->tree = tree;

  for (i = 0; i < tree->nnodes; i++) {
    TreeNode *n = lst_get_ptr(tree->nodes, i);
    if (n != tree)
      im_set_branch(im->branch_mods[i], alpha[i], beta[i], tau[i], n->dparent);
  }
}

static PHAST_INLINE
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

double im_branch_column_logl(IndelHistory *ih, BranchIndelModel *bim, 
                             int child, double *col_logl) {
  int i;
  col_type this_type, last_type = SKIP;
  double logl = 0;
  int parent_id = ((TreeNode*)lst_get_ptr(ih->tree->nodes, child))->parent->id;

  for (i = 0; i < ih->ncols; i++) {
    last_type = get_col_type(ih, child, parent_id, i);
    if (last_type == SKIP) 
      col_logl[i] = 0;
    else {
      col_logl[i] = vec_get(bim->beg_log_probs, last_type);
      break;
    }
  }
  logl = col_logl[i];
  i++;

  for (; i < ih->ncols; i++) {
    this_type = get_col_type(ih, child, parent_id, i);
    
    if (this_type == ERROR)
      die("ERROR im_branch_column_logl\n");

    if (this_type == SKIP) {
      col_logl[i] = 0;
      continue;
    }

    col_logl[i] = bim->log_probs->data[last_type][this_type];
    logl += col_logl[i];
    last_type = this_type;
  }

  return logl;
}

/* returns total log likelihood */
double im_column_logl(IndelHistory *ih, IndelModel *im, double *col_logl) {
  int i, j;
  double logl = 0;
  double *branch_col_logl = smalloc(ih->ncols * sizeof(double));
  for (i = 0; i < ih->ncols; i++) col_logl[i] = 0;
  for (i = 0; i < im->tree->nnodes; i++) {
    TreeNode *n = lst_get_ptr(im->tree->nodes, i);
    if (n == im->tree) continue;
    logl += im_branch_column_logl(ih, im->branch_mods[n->id], i, 
                                  branch_col_logl);
    for (j = 0; j < ih->ncols; j++) col_logl[j] += branch_col_logl[j];
  }
  sfree(branch_col_logl);

  /* NOTE: it looks like we're combining joint probabilities over
     indel strings here, when we should be using conditional
     probabilities, but what we're really combining are the
     probabilities of indel *events* along branches, which are
     independent given the indel history */

  return logl;
}

double im_likelihood(IndelModel *im, IndelSuffStats *iss) {
  double logl = 0;
  int i, j, k;
  checkInterrupt();
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

BranchIndelSuffStats *im_suff_stats_branch(IndelHistory *ih, int child_id) {
  int i, j;
  char c;
  col_type this_type, last_type = SKIP;
  int parent_id = ((TreeNode*)lst_get_ptr(ih->tree->nodes, child_id))->parent->id;
  BranchIndelSuffStats *ss = smalloc(sizeof(BranchIndelSuffStats));
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
      die("ERROR im_suff_stats_branch\n");
    }

    else if (this_type == SKIP) continue;

    ss->trans_counts->data[last_type][this_type]++;
    last_type = this_type;
  }

  return ss;  
}

IndelSuffStats *im_suff_stats(IndelHistory *ih) {
  int i;
  IndelSuffStats *iss = smalloc(sizeof(IndelSuffStats));
  iss->tree = ih->tree;
  iss->branch_counts = smalloc(ih->tree->nnodes * sizeof(void*));
  for (i = 0; i < ih->tree->nnodes; i++) {
    TreeNode *n = lst_get_ptr(ih->tree->nodes, i);
    if (n == ih->tree) {
      iss->branch_counts[i] = NULL;
      continue;
    }
    iss->branch_counts[i] = im_suff_stats_branch(ih, n->id);
  }
  return iss;
}

void im_free_suff_stats(IndelSuffStats *iss) {
  int i;
  for (i = 0; i < iss->tree->nnodes; i++) {
    if (iss->branch_counts[i] != NULL) {
      mat_free(iss->branch_counts[i]->trans_counts);
      vec_free(iss->branch_counts[i]->beg_counts);
      sfree(iss->branch_counts[i]);
    }
  }
  sfree(iss->branch_counts);
  sfree(iss);
}

/* double im_simulate_history(IndelModel *tim, int ncols) { */
/* } */

/* package of data used below */
struct likelihood_data {
  IndelModel *im;
  IndelSuffStats *ss;
};

/* wrapper for likelihood function, for use in numerical maximization */
double im_likelihood_wrapper(Vector *params, void *data) {
  struct likelihood_data *d = data;
  im_set_all(d->im, vec_get(params, 0), vec_get(params, 1), vec_get(params, 2), 
             d->im->tree);
  return -im_likelihood(d->im, d->ss);
}

/* computes gradient of (negative) likelihood, for use in numerical
   maximization */
void im_likelihood_gradient(Vector *grad, Vector *params, void *data, 
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
void im_estimate(IndelModel *im, IndelHistory *ih, IndelSuffStats *ss, 
                 FILE *logf) {
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

  opt_bfgs(im_likelihood_wrapper, params, d, &neglogl, lb, ub, logf,  
           im_likelihood_gradient, OPT_HIGH_PREC, NULL);  

  im_set_all(im, vec_get(params, 0), vec_get(params, 1), 
             vec_get(params, 2), im->tree);
  im->training_lnl = -neglogl * log(2);
  
  vec_free(params);
  vec_free(lb);
  im_free_suff_stats(d->ss);
  sfree(d);
}

/* collect sufficient stats for a branch, considering only sites in
   the specified category */
BranchIndelSuffStats *im_suff_stats_branch_cat(IndelHistory *ih, int child_id,
                                               int *categories, int do_cat) {
  int i, j;
  char c;
  col_type this_type=SKIP, last_type;
  int parent_id = ((TreeNode*)lst_get_ptr(ih->tree->nodes, child_id))->parent->id;
  BranchIndelSuffStats *ss = smalloc(sizeof(BranchIndelSuffStats));
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
    checkInterruptN(i, 1000);
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
      die("ERROR im_suff_stats_branch_cat\n");
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
IndelSuffStats *im_suff_stats_cat(IndelHistory *ih, int *categories, 
                                  int do_cat) {
  int i;
  IndelSuffStats *iss = smalloc(sizeof(IndelSuffStats));
  iss->tree = ih->tree;
  iss->branch_counts = smalloc(ih->tree->nnodes * sizeof(void*));
  for (i = 0; i < ih->tree->nnodes; i++) {
    TreeNode *n = lst_get_ptr(ih->tree->nodes, i);
    if (n == ih->tree) {
      iss->branch_counts[i] = NULL;
      continue;
    }
    iss->branch_counts[i] = im_suff_stats_branch_cat(ih, n->id, 
                                                     categories, do_cat);
  }
  return iss;
}

