/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: fit_column.c,v 1.24 2009-02-19 17:23:17 acs Exp $ */

/* Functions to compute likelihoods for individual alignment columns,
   estimate column-by-column scale factors by maximum likelihood,
   perform single-base LRTs, score tests, phyloP, etc. */

#include <stdlib.h>
#include <fit_column.h>
#include <sufficient_stats.h>
#include <tree_likelihoods.h>
#include <time.h>

#define DERIV_EPSILON 1e-6
/* for numerical computation of derivatives */

#define NSAMPLES_FIM 50
/* number of samples to use in estimating FIM */

#define SIGFIGS 4
/* number of significant figures to which to estimate column scale
   parameters (currently affects 1d parameter estimation only) */

/* Compute and return the log likelihood of a tree model with respect
   to a single column tuple in an alignment.  This is a pared-down
   version of tl_compute_log_likelihood for use in estimation of
   base-by-base scale factors.  It assumes a 0th order model,
   leaf-to-sequence mapping already available, prob matrices computed,
   sufficient stats already available.  Note that this function uses
   natural log rather than log2.  This function does allow for rate
   variation. */
double col_compute_likelihood(TreeModel *mod, MSA *msa, int tupleidx,
                                  double **scratch) {

  int i, j, k, nodeidx, rcat;
  int nstates = mod->rate_matrix->size;
  TreeNode *n;
  double total_prob = 0;
  List *traversal = tr_postorder(mod->tree);
  double **pL = NULL;

  if (msa->ss->tuple_size != 1)
    die("ERROR col_compute_likelihood: need tuple size 1, got %i\n",
	msa->ss->tuple_size);
  if (mod->order != 0)
    die("ERROR col_compute_likelihood: got mod->order of %i, expected 0\n",
	mod->order);
  if (!mod->allow_gaps)
    die("ERROR col_compute_likelihood: need mod->allow_gaps to be TRUE\n");

  /* allocate memory or use scratch if avail */
  if (scratch != NULL)
    pL = scratch;
  else {
    pL = smalloc(nstates * sizeof(double*));
    for (j = 0; j < nstates; j++)
      pL[j] = smalloc((mod->tree->nnodes+1) * sizeof(double));
  }

  for (rcat = 0; rcat < mod->nratecats; rcat++) {
    for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
      n = lst_get_ptr(traversal, nodeidx);
      if (n->lchild == NULL) {
        /* leaf: base case of recursion */
        int state = mod->rate_matrix->
          inv_states[(int)ss_get_char_tuple(msa, tupleidx,
                                            mod->msa_seq_idx[n->id], 0)];
        for (i = 0; i < nstates; i++) {
          if (state < 0 || i == state)
            pL[i][n->id] = 1;
          else
            pL[i][n->id] = 0;
        }
      }
      else {
        /* general recursive case */
        MarkovMatrix *lsubst_mat = mod->P[n->lchild->id][rcat];
        MarkovMatrix *rsubst_mat = mod->P[n->rchild->id][rcat];
        for (i = 0; i < nstates; i++) {
          double totl = 0, totr = 0;
          for (j = 0; j < nstates; j++)
            totl += pL[j][n->lchild->id] *
              mm_get(lsubst_mat, i, j);

          for (k = 0; k < nstates; k++)
            totr += pL[k][n->rchild->id] *
              mm_get(rsubst_mat, i, k);

          pL[i][n->id] = totl * totr;
        }
      }
    }

    /* termination (for each rate cat) */
    for (i = 0; i < nstates; i++)
      total_prob += vec_get(mod->backgd_freqs, i) *
        pL[i][mod->tree->id] * mod->freqK[rcat];
  }

  if (scratch == NULL) {
    for (j = 0; j < nstates; j++) sfree(pL[j]);
    sfree(pL);
  }

  return(total_prob);
}



/* See col_compute_likelihood above for notes.
   Note that this function uses natural log rather than log2
 */
double col_compute_log_likelihood(TreeModel *mod, MSA *msa, int tupleidx,
                                  double **scratch) {
  return log(col_compute_likelihood(mod, msa, tupleidx, scratch));
}


/* version of col_scale_derivs_subst that allows for the general case
   of complex eigenvalues and eigenvectors */
void col_scale_derivs_subst_complex(ColFitData *d) {
  Zmatrix *S = d->mod->rate_matrix->evec_matrix_z,
    *Sinv = d->mod->rate_matrix->evec_matrix_inv_z;
  MarkovMatrix *Q = d->mod->rate_matrix;
  int size = Q->size;
  int rcat, nid, i;

  if (S==NULL)
    die("ERROR col_scale_derivs_subst_complex: got S==NULL\n");
  if (Sinv==NULL)
    die("ERROR col_scale_derivs_subst_complex: got Sinv==NULL\n");
  if (d->mod->alt_subst_mods != NULL)
    die("ERROR col_scale_derivs_subst_complex cannot handle lineage-specific models");

  for (rcat = 0; rcat < d->mod->nratecats; rcat++) {
    for (nid = 1; nid < d->mod->tree->nnodes; nid++) { /* skip root */

      double t = ((TreeNode*)lst_get_ptr(d->mod->tree->nodes, nid))->dparent;
      double l1 = d->mod->scale;
      double l2 = (d->stype == SUBTREE && d->mod->in_subtree[nid] ?
                   d->mod->scale_sub : 1);

      /* set up exponentiated diagonal matrix */
      for (i = 0; i < size; i++) {
        double r = t * l1 * l2 * d->mod->rK[rcat];
        zvec_set(d->expdiag_z, i, z_exp(z_mul_real(zvec_get(Q->evals_z, i), r)));
      }

      /* PP */
      zvec_copy(d->vec_scratch1_z, Q->evals_z);
      zvec_scale(d->vec_scratch1_z, t * l2);
      zvec_had_prod(d->vec_scratch2_z, d->vec_scratch1_z, d->expdiag_z);
      zmat_mult_real_diag(d->PP[nid][rcat], S, d->vec_scratch2_z, Sinv,
                          d->mat_scratch_z);

      if (d->second_derivs) {
        /* PPP */
        zvec_had_prod(d->vec_scratch2_z, d->vec_scratch1_z, d->vec_scratch1_z);
        zvec_had_prod(d->vec_scratch1_z, d->vec_scratch2_z, d->expdiag_z);
        zmat_mult_real_diag(d->PPP[nid][rcat], S, d->vec_scratch1_z, Sinv,
                            d->mat_scratch_z);
      }

      if (d->stype == SUBTREE && d->mod->in_subtree[nid]) {
        /* if not in subtree, leave all of these equal to 0 (as
           initialized) */

        /* QQ */
        zvec_copy(d->vec_scratch1_z, Q->evals_z);
        zvec_scale(d->vec_scratch1_z, t * l1);
        zvec_had_prod(d->vec_scratch2_z, d->vec_scratch1_z, d->expdiag_z);
        zmat_mult_real_diag(d->QQ[nid][rcat], S, d->vec_scratch2_z, Sinv,
                            d->mat_scratch_z);

        if (d->second_derivs) {
          /* QQQ */
          zvec_had_prod(d->vec_scratch2_z, d->vec_scratch1_z, d->vec_scratch1_z);
          zvec_had_prod(d->vec_scratch1_z, d->vec_scratch2_z, d->expdiag_z);
          zmat_mult_real_diag(d->QQQ[nid][rcat], S, d->vec_scratch1_z, Sinv,
                              d->mat_scratch_z);

          /* RRR */
          zvec_copy(d->vec_scratch1_z, Q->evals_z);
          zvec_scale(d->vec_scratch1_z, t);
          zvec_had_prod(d->vec_scratch2_z, Q->evals_z, Q->evals_z);
          zvec_scale(d->vec_scratch2_z, t * t * l1 * l2);
          zvec_plus_eq(d->vec_scratch2_z, d->vec_scratch1_z);
          zvec_had_prod(d->vec_scratch1_z, d->vec_scratch2_z, d->expdiag_z);
          zmat_mult_real_diag(d->RRR[nid][rcat], S, d->vec_scratch1_z, Sinv,
                              d->mat_scratch_z);
        }
      }
    }
  }
}

/* version of col_scale_derivs_subst that is optimized for the case in
   which eigenvalues and eigenvectors can be assumed to be real */
void col_scale_derivs_subst_real(ColFitData *d) {
  Matrix *S = d->mod->rate_matrix->evec_matrix_r,
    *Sinv = d->mod->rate_matrix->evec_matrix_inv_r;
  MarkovMatrix *Q = d->mod->rate_matrix;
  int size = Q->size;
  int rcat, nid, i;

  if (S==NULL)
    die("ERROR col_scale_derivs_subst_real: got S==NULL\n");
  if (Sinv==NULL)
    die("ERROR col_scale_derivs_subst_real: got Sinv==NULL\n");
  if (d->mod->alt_subst_mods != NULL)
    die("ERROR col_scale_derivs_subst_real: cannot handle lineage-specific models");

  for (rcat = 0; rcat < d->mod->nratecats; rcat++) {
    for (nid = 1; nid < d->mod->tree->nnodes; nid++) { /* skip root */

      double t = ((TreeNode*)lst_get_ptr(d->mod->tree->nodes, nid))->dparent;
      double l1 = d->mod->scale;
      double l2 = (d->stype == SUBTREE && d->mod->in_subtree[nid] ?
                   d->mod->scale_sub : 1);

      /* set up exponentiated diagonal matrix */
      for (i = 0; i < size; i++)
        d->expdiag_r->data[i] = exp(Q->evals_r->data[i] * t *
                                    l1 * l2 * d->mod->rK[rcat]);

      /* PP */
      for (i = 0; i < size; i++) {
        d->vec_scratch1_r->data[i] = Q->evals_r->data[i] * t * l2;
        d->vec_scratch2_r->data[i] = d->vec_scratch1_r->data[i] *
          d->expdiag_r->data[i];
      }
      mat_mult_diag(d->PP[nid][rcat], S, d->vec_scratch2_r, Sinv);

      if (d->second_derivs) {
        /* PPP */
        for (i = 0; i < size; i++) {
          d->vec_scratch2_r->data[i] = d->vec_scratch1_r->data[i] *
            d->vec_scratch1_r->data[i];
          d->vec_scratch1_r->data[i] = d->vec_scratch2_r->data[i] *
            d->expdiag_r->data[i];
        }
        mat_mult_diag(d->PPP[nid][rcat], S, d->vec_scratch1_r, Sinv);
      }

      if (d->stype == SUBTREE && d->mod->in_subtree[nid]) {
        /* if not in subtree, leave all of these equal to 0 (as
           initialized) */

        /* QQ */
        for (i = 0; i < size; i++) {
          d->vec_scratch1_r->data[i] = Q->evals_r->data[i] * t * l1;
          d->vec_scratch2_r->data[i] = d->vec_scratch1_r->data[i] *
            d->expdiag_r->data[i];
        }
        mat_mult_diag(d->QQ[nid][rcat], S, d->vec_scratch2_r, Sinv);

        if (d->second_derivs) {
          /* QQQ */
          for (i = 0; i < size; i++) {
            d->vec_scratch2_r->data[i] = d->vec_scratch1_r->data[i] *
              d->vec_scratch1_r->data[i];
            d->vec_scratch1_r->data[i] = d->vec_scratch2_r->data[i] *
              d->expdiag_r->data[i];
          }
          mat_mult_diag(d->QQQ[nid][rcat], S, d->vec_scratch1_r, Sinv);

          /* RRR */
          for (i = 0; i < size; i++)
            d->vec_scratch1_r->data[i] = d->expdiag_r->data[i] *
              (Q->evals_r->data[i] * Q->evals_r->data[i] * t * t * l1 * l2 +
               Q->evals_r->data[i] * t);
          mat_mult_diag(d->RRR[nid][rcat], S, d->vec_scratch1_r, Sinv);
        }
      }
    }
  }
}

/* Compute 1st and 2nd derivs wrt scale params of substitution
   matrices for each branch of the tree (and each rate category).
   These are used in the recursive computation of derivatives of the
   column likelihoods */
void col_scale_derivs_subst(ColFitData *d) {
  /* now broken into real and complex cases for efficiency */
  if (d->mod->rate_matrix->eigentype == REAL_NUM)
    col_scale_derivs_subst_real(d);
  else
    col_scale_derivs_subst_complex(d);
}

/* Compute the first and (optionally) second derivatives with respect
   to the scale parameter for the single-column log likelihood
   function (col_compute_log_likelihood).  This version assumes a
   single scale parameter; see below for the subtree case.  Return
   value is log likelihood, which is computed as a by-product.  Derivs
   will be stored in *first_deriv and *second_deriv.  If second_deriv
   == NULL, it will not be computed (saves some time).  */
double col_scale_derivs(ColFitData *d, double *first_deriv,
                        double *second_deriv, double ***scratch) {

  int i, j, k, nodeidx, rcat;
  int nstates = d->mod->rate_matrix->size;
  TreeNode *n;
  double total_prob = 0;
  List *traversal = tr_postorder(d->mod->tree);
  double **L=NULL;                   /* partial likelihoods */
  double **LL=NULL;                  /* 1st deriv of partial likelihoods wrt
                                   scale param */
  double **LLL=NULL;                 /* 2nd deriv of partial likelihoods
                                   wrt scale param */
  if (d->msa->ss->tuple_size != 1)
    die("ERROR col_scale_derivs: need tuple size 1, got %i\n",
	d->msa->ss->tuple_size);
  if (d->mod->order != 0)
    die("ERROR col_scale_derivs: got mod->order of %i, expected 0\n",
	d->mod->order);
  if (!d->mod->allow_gaps)
    die("ERROR col_scale_derivs: need mod->allow_gaps to be TRUE\n");

  *first_deriv = 0;
  if (second_deriv != NULL) *second_deriv = 0;

  /* allocate memory or use scratch if available */
  if (scratch == NULL) {
    L = smalloc(nstates * sizeof(double*));
    LL = smalloc(nstates * sizeof(double*));
    if (second_deriv != NULL)
      LLL = smalloc(nstates * sizeof(double*));
    for (j = 0; j < nstates; j++) {
      L[j] = smalloc((d->mod->tree->nnodes+1) * sizeof(double));
      LL[j] = smalloc((d->mod->tree->nnodes+1) * sizeof(double));
      if (second_deriv != NULL)
        LLL[j] = smalloc((d->mod->tree->nnodes+1) * sizeof(double));
    }
  }
  else {
    L = scratch[0];
    LL = scratch[1];
    LLL = scratch[2];
  }

  col_scale_derivs_subst(d);

  for (rcat = 0; rcat < d->mod->nratecats; rcat++) {
    for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
      n = lst_get_ptr(traversal, nodeidx);
      if (n->lchild == NULL) {
        /* leaf: base case of recursion */
        int state = d->mod->rate_matrix->
          inv_states[(int)ss_get_char_tuple(d->msa, d->tupleidx,
                                            d->mod->msa_seq_idx[n->id], 0)];
        for (i = 0; i < nstates; i++) {
          if (state < 0 || i == state)
            L[i][n->id] = 1;
          else
            L[i][n->id] = 0;

          LL[i][n->id] = 0;
          if (second_deriv != NULL) LLL[i][n->id] = 0;
        }
      }
      else {
        /* general recursive case */
        MarkovMatrix *lsubst_mat = d->mod->P[n->lchild->id][rcat];
        MarkovMatrix *rsubst_mat = d->mod->P[n->rchild->id][rcat];
        for (i = 0; i < nstates; i++) {
          double totl = 0, totr = 0, A = 0, B = 0, E = 0, F = 0;
          for (j = 0; j < nstates; j++) {
            totl += L[j][n->lchild->id] * mm_get(lsubst_mat, i, j);

            A += (L[j][n->lchild->id] * d->PP[n->lchild->id][rcat]->data[i][j]) +
              (LL[j][n->lchild->id] * mm_get(lsubst_mat, i, j));
          }

          for (k = 0; k < nstates; k++) {
            totr += L[k][n->rchild->id] * mm_get(rsubst_mat, i, k);

            B += (L[k][n->rchild->id] * d->PP[n->rchild->id][rcat]->data[i][k]) +
              (LL[k][n->rchild->id] * mm_get(rsubst_mat, i, k));

          }

          L[i][n->id] = totl * totr;
          LL[i][n->id] = totr*A + totl*B;

          if (second_deriv != NULL) {
            for (j = 0; j < nstates; j++)
              E += L[j][n->lchild->id] * d->PPP[n->lchild->id][rcat]->data[i][j] +
                2 * LL[j][n->lchild->id] * d->PP[n->lchild->id][rcat]->data[i][j] +
                LLL[j][n->lchild->id] * mm_get(lsubst_mat, i, j);

            for (k = 0; k < nstates; k++)
              F += L[k][n->rchild->id] * d->PPP[n->rchild->id][rcat]->data[i][k] +
                2 * LL[k][n->rchild->id] * d->PP[n->rchild->id][rcat]->data[i][k] +
                LLL[k][n->rchild->id] * mm_get(rsubst_mat, i, k);

            LLL[i][n->id] = totr*E + 2*A*B + totl*F;
          }
        }
      }
    }

    /* termination (for each rate cat) */
    for (i = 0; i < nstates; i++) {
      total_prob += L[i][d->mod->tree->id] * vec_get(d->mod->backgd_freqs, i) *
        d->mod->freqK[rcat];

      *first_deriv += LL[i][d->mod->tree->id] * vec_get(d->mod->backgd_freqs, i) *
        d->mod->freqK[rcat];

      if (second_deriv != NULL)
        *second_deriv += LLL[i][d->mod->tree->id] * vec_get(d->mod->backgd_freqs, i) *
          d->mod->freqK[rcat];
    }
  }

  /* convert to log space */
  if (second_deriv != NULL)
    *second_deriv = (*second_deriv)/total_prob -
      ((*first_deriv)/total_prob) * ((*first_deriv)/total_prob);
                                /* deriv of log followed by quotient
                                   rule; rearrange terms to avoid
                                   underflow */

  *first_deriv = *first_deriv / total_prob; /* deriv of log */
  total_prob = log(total_prob);

  if (scratch == NULL) {
    for (j = 0; j < nstates; j++) {
      sfree(L[j]);
      sfree(LL[j]);
      if (second_deriv != NULL) sfree(LLL[j]);
    }
    sfree(L);
    sfree(LL);
    if (second_deriv != NULL) sfree(LLL);
  }

  return(total_prob);
}

/* Compute the first and (optionally) second derivatives with respect
   to the scale parameters for the single-column log likelihood
   function (col_compute_log_likelihood).  This version assumes scale
   parameters for the whole tree and for the subtree.  Return value is
   log likelihood, which is computed as a by-product.  Derivs will be
   stored in *gradient and *hessian.  If hessian == NULL,
   it will not be computed (saves some time).  */
double col_scale_derivs_subtree(ColFitData *d, Vector *gradient,
                                Matrix *hessian, double ***scratch) {
  int i, j, k, nodeidx, rcat;
  int nstates = d->mod->rate_matrix->size;
  TreeNode *n;
  double total_prob = 0;
  List *traversal = tr_postorder(d->mod->tree);
  double **L=NULL;                   /* partial likelihoods */
  double **LL=NULL;                  /* 1st deriv of partial likelihoods wrt
                                   1st scale param */
  double **LLL=NULL;                 /* 2nd deriv of partial likelihoods
                                   wrt 1st scale param */
  double **MM=NULL;                  /* 1st deriv of partial likelihoods
                                   wrt 2nd scale param */
  double **MMM=NULL;                 /* 2nd deriv of partial likelihoods
                                   wrt 2nd scale param */
  double **NNN=NULL;                 /* 2nd cross deriv (off diagonal in
                                   Hessian) of partial likelihoods */

  double *pd = gradient->data;  /* 1st partial derivatives */
  double **pd2 = (hessian == NULL ? NULL : hessian->data);
                                /* 2nd partial derivatives; because of
                                   symmetry, only pd2[0][0],
                                   pd2[1][1], and pd2[1][0] need to be
                                   considered during computation */
  if (d->msa->ss->tuple_size != 1)
    die("ERROR col_scale_derivs_subtree: need tuple size 1, got %i\n",
	d->msa->ss->tuple_size);
  if (d->mod->order != 0)
    die("ERROR col_scale_derivs_subtree: got mod->order of %i, expected 0\n",
	d->mod->order);
  if (!d->mod->allow_gaps)
    die("ERROR col_scale_derivs_subtree: need mod->allow_gaps to be TRUE\n");

  pd[0] = pd[1] = 0;
  if (pd2 != NULL)
    pd2[0][0] = pd2[1][1] = pd2[0][1] = pd2[1][0] = 0;

  /* allocate memory or use scratch if available */
  if (scratch == NULL) {
    L = smalloc(nstates * sizeof(double*));
    LL = smalloc(nstates * sizeof(double*));
    MM = smalloc(nstates * sizeof(double*));
    if (pd2 != NULL) {
      LLL = smalloc(nstates * sizeof(double*));
      MMM = smalloc(nstates * sizeof(double*));
      NNN = smalloc(nstates * sizeof(double*));
    }
    for (j = 0; j < nstates; j++) {
      L[j] = smalloc((d->mod->tree->nnodes+1) * sizeof(double));
      LL[j] = smalloc((d->mod->tree->nnodes+1) * sizeof(double));
      MM[j] = smalloc((d->mod->tree->nnodes+1) * sizeof(double));
      if (pd2 != NULL) {
        LLL[j] = smalloc((d->mod->tree->nnodes+1) * sizeof(double));
        MMM[j] = smalloc((d->mod->tree->nnodes+1) * sizeof(double));
        NNN[j] = smalloc((d->mod->tree->nnodes+1) * sizeof(double));
      }
    }
  }
  else {
    L = scratch[0];
    LL = scratch[1];
    LLL = scratch[2];
    MM = scratch[3];
    MMM = scratch[4];
    NNN = scratch[5];
  }

  col_scale_derivs_subst(d);

  for (rcat = 0; rcat < d->mod->nratecats; rcat++) {
    for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
      n = lst_get_ptr(traversal, nodeidx);
      if (n->lchild == NULL) {
        /* leaf: base case of recursion */
        int state = d->mod->rate_matrix->
          inv_states[(int)ss_get_char_tuple(d->msa, d->tupleidx,
                                            d->mod->msa_seq_idx[n->id], 0)];
        for (i = 0; i < nstates; i++) {
          if (state < 0 || i == state)
            L[i][n->id] = 1;
          else
            L[i][n->id] = 0;

          LL[i][n->id] = MM[i][n->id] = 0;
          if (pd2 != NULL)
            LLL[i][n->id] = MMM[i][n->id] = NNN[i][n->id] = 0;
        }
      }
      else {
        /* general recursive case */
        MarkovMatrix *lsubst_mat = d->mod->P[n->lchild->id][rcat];
        MarkovMatrix *rsubst_mat = d->mod->P[n->rchild->id][rcat];
        for (i = 0; i < nstates; i++) {
          double totl = 0, totr = 0, A = 0, B = 0, C = 0, D = 0, E = 0,
            F = 0, G = 0, H = 0, I = 0, J = 0;
          for (j = 0; j < nstates; j++) {
            totl += L[j][n->lchild->id] * mm_get(lsubst_mat, i, j);

            A += (L[j][n->lchild->id] * d->PP[n->lchild->id][rcat]->data[i][j]) +
              (LL[j][n->lchild->id] * mm_get(lsubst_mat, i, j));

            C += (L[j][n->lchild->id] * d->QQ[n->lchild->id][rcat]->data[i][j]) +
              (MM[j][n->lchild->id] * mm_get(lsubst_mat, i, j));
          }

          for (k = 0; k < nstates; k++) {
            totr += L[k][n->rchild->id] * mm_get(rsubst_mat, i, k);

            B += (L[k][n->rchild->id] * d->PP[n->rchild->id][rcat]->data[i][k]) +
              (LL[k][n->rchild->id] * mm_get(rsubst_mat, i, k));

            D += (L[k][n->rchild->id] * d->QQ[n->rchild->id][rcat]->data[i][k]) +
              (MM[k][n->rchild->id] * mm_get(rsubst_mat, i, k));

          }

          L[i][n->id] = totl * totr;
          LL[i][n->id] = totr*A + totl*B;
          MM[i][n->id] = totr*C + totl*D;

          if (pd2 != NULL) {
            for (j = 0; j < nstates; j++) {
              E += L[j][n->lchild->id] * d->PPP[n->lchild->id][rcat]->data[i][j] +
                2 * LL[j][n->lchild->id] * d->PP[n->lchild->id][rcat]->data[i][j] +
                LLL[j][n->lchild->id] * mm_get(lsubst_mat, i, j);
              G += L[j][n->lchild->id] * d->QQQ[n->lchild->id][rcat]->data[i][j] +
                2 * MM[j][n->lchild->id] * d->QQ[n->lchild->id][rcat]->data[i][j] +
                MMM[j][n->lchild->id] * mm_get(lsubst_mat, i, j);
              I += L[j][n->lchild->id] * d->RRR[n->lchild->id][rcat]->data[i][j] +
                MM[j][n->lchild->id] * d->PP[n->lchild->id][rcat]->data[i][j] +
                LL[j][n->lchild->id] * d->QQ[n->lchild->id][rcat]->data[i][j] +
                NNN[j][n->lchild->id] * mm_get(lsubst_mat, i, j);
            }

            for (k = 0; k < nstates; k++) {
              F += L[k][n->rchild->id] * d->PPP[n->rchild->id][rcat]->data[i][k] +
                2 * LL[k][n->rchild->id] * d->PP[n->rchild->id][rcat]->data[i][k] +
                LLL[k][n->rchild->id] * mm_get(rsubst_mat, i, k);
              H += L[k][n->rchild->id] * d->QQQ[n->rchild->id][rcat]->data[i][k] +
                2 * MM[k][n->rchild->id] * d->QQ[n->rchild->id][rcat]->data[i][k] +
                MMM[k][n->rchild->id] * mm_get(rsubst_mat, i, k);
              J += L[k][n->rchild->id] * d->RRR[n->rchild->id][rcat]->data[i][k] +
                MM[k][n->rchild->id] * d->PP[n->rchild->id][rcat]->data[i][k] +
                LL[k][n->rchild->id] * d->QQ[n->rchild->id][rcat]->data[i][k] +
                NNN[k][n->rchild->id] * mm_get(rsubst_mat, i, k);
            }

            LLL[i][n->id] = totr*E + 2*A*B + totl*F;
            MMM[i][n->id] = totr*G + 2*C*D + totl*H;
            NNN[i][n->id] = totr*I + A*D + B*C + totl*J;
          }
        }
      }
    }

    /* termination (for each rate cat) */
    for (i = 0; i < nstates; i++) {
      total_prob += L[i][d->mod->tree->id] * vec_get(d->mod->backgd_freqs, i) *
        d->mod->freqK[rcat];

      pd[0] += LL[i][d->mod->tree->id] * vec_get(d->mod->backgd_freqs, i) *
        d->mod->freqK[rcat];
      pd[1] += MM[i][d->mod->tree->id] * vec_get(d->mod->backgd_freqs, i) *
        d->mod->freqK[rcat];

      if (pd2 != NULL) {
        pd2[0][0] += LLL[i][d->mod->tree->id] * vec_get(d->mod->backgd_freqs, i) *
          d->mod->freqK[rcat];
        pd2[1][1] += MMM[i][d->mod->tree->id] * vec_get(d->mod->backgd_freqs, i) *
          d->mod->freqK[rcat];
        pd2[1][0] += NNN[i][d->mod->tree->id] * vec_get(d->mod->backgd_freqs, i) *
          d->mod->freqK[rcat];
      }
    }
 }

  /* convert to log space */
  if (pd2 != NULL) {
    /* deriv of log and quotient rule */
    pd2[0][0] =  pd2[0][0]/total_prob - (pd[0]/total_prob)*(pd[0]/total_prob);
    pd2[1][1] =  pd2[1][1]/total_prob - (pd[1]/total_prob)*(pd[1]/total_prob);
    pd2[1][0] =  pd2[1][0]/total_prob - (pd[1]/total_prob)*(pd[0]/total_prob);
    pd2[0][1] =  pd2[1][0];
  }
  pd[0] = pd[0] / total_prob; /* deriv of log */
  pd[1] = pd[1] / total_prob;
  total_prob = log(total_prob);

  if (scratch == NULL) {
    for (j = 0; j < nstates; j++) {
      sfree(L[j]);
      sfree(LL[j]);
      sfree(MM[j]);
      if (pd2 != NULL) {
        sfree(LLL[j]);
        sfree(MMM[j]);
        sfree(NNN[j]);
      }
    }
    sfree(L);
    sfree(LL);
    sfree(MM);
    if (pd2 != NULL) {
      sfree(LLL);
      sfree(MMM);
      sfree(NNN);
    }
  }

  return(total_prob);
}


/* Wrapper for likelihood function for use in parameter estimation */
double col_likelihood_wrapper(Vector *params, void *data) {
  ColFitData *d = (ColFitData*)data;

  d->mod->scale = vec_get(params, 0);
  if (d->stype == SUBTREE)
    d->mod->scale_sub = vec_get(params, 1);

  /* reestimate subst models on edges */
  tm_set_subst_matrices(d->mod);

  return -1 * col_compute_log_likelihood(d->mod, d->msa, d->tupleidx,
                                         d->fels_scratch[0]);
}

/* Wrapper for likelihood function for use in parameter estimation;
   version for use with opt_newton_1d */
double col_likelihood_wrapper_1d(double x, void *data) {
  ColFitData *d = (ColFitData*)data;
  if (d->stype == SUBTREE)
    die("ERROR col_likelihood_wrapper_1d: d->stype cannot be SUBTREE\n");

  d->mod->scale = x;

  /* reestimate subst models on edges */
  tm_set_subst_matrices(d->mod);

  return -1 * col_compute_log_likelihood(d->mod, d->msa, d->tupleidx,
                                         d->fels_scratch[0]);
}

/* Wrapper for gradient function for use in parameter estimation */
void col_grad_wrapper(Vector *grad, Vector *params, void *data,
                      Vector *lb, Vector *ub) {
  ColFitData *d = (ColFitData*)data;
  double deriv;

  if (d->stype == ALL) {
    col_scale_derivs(d, &deriv, NULL, d->fels_scratch);
    vec_set(grad, 0, -deriv);   /* because working with neg lnl */
  }
  else {
    col_scale_derivs_subtree(d, grad, NULL, d->fels_scratch);
    vec_scale(grad, -1);
  }
}

/* Wrapper for gradient function for use in parameter estimation;
   version for use with opt_newton_1d */
double col_grad_wrapper_1d(double x, void *data, double lb, double ub) {
  double deriv, deriv2;
  ColFitData *d = (ColFitData*)data;
  if (d->stype != ALL)
    die("ERROR col_grad_wrapper_1d: d->stype must be ALL\n");
  col_scale_derivs(d, &deriv, &deriv2, d->fels_scratch);
  d->deriv2 = -deriv2;           /* store for use by wrapper below */
  return -deriv; /* because working with neg lnl */
}

/* Wrapper for second derivative function for use in parameter
   estimation; version for use with opt_newton_1d.  Simply returns
   value computed in col_grad_wrapper_1d (1st and 2nd derivs are
   computed simultaneously) */
double col_deriv2_wrapper_1d(double x, void *data, double lb, double ub) {
  ColFitData *d = (ColFitData*)data;
  return d->deriv2;
}

/* Perform a likelihood ratio test for each column tuple in an
   alignment, comparing the given null model with an alternative model
   that has a free scaling parameter for all branches.  Assumes a 0th
   order model, leaf-to-sequence mapping already available, prob
   matrices computed, sufficient stats available.  Computes p-values
   based using the chi-sq distribution and stores them in tuple_pvals.
   Will optionally store the individual scale factors in tuple_scales
   and raw log likelihood ratios in tuple_llrs if these variables are
   non-NULL.  Must define mode as CON (for 0 <= scale <= 1), ACC
   (for 1 <= scale), NNEUT (0 <= scale), or CONACC (0 <= scale) */
void col_lrts(TreeModel *mod, MSA *msa, mode_type mode, double *tuple_pvals,
              double *tuple_scales, double *tuple_llrs, FILE *logf) {
  int i;
  ColFitData *d;
  double null_lnl, alt_lnl, delta_lnl, this_scale = 1;

  /* init ColFitData */
  d = col_init_fit_data(mod, msa, ALL, mode, FALSE);

  /* iterate through column tuples */
  for (i = 0; i < msa->ss->ntuples; i++) {
    checkInterruptN(i, 100);

    /* first check for actual substitution data in column; if none,
       don't waste time computing likelihoods */
    if (!col_has_data(mod, msa, i)) {
      delta_lnl = 0;
      this_scale = 1;
    }

    else {                      /* compute null and alt lnl */
      mod->scale = 1;
      tm_set_subst_matrices(mod);

      /* compute log likelihoods under null and alt hypotheses */
      null_lnl = col_compute_log_likelihood(mod, msa, i, d->fels_scratch[0]);

      vec_set(d->params, 0, d->init_scale);
      d->tupleidx = i;

      opt_newton_1d(col_likelihood_wrapper_1d, &d->params->data[0], d,
                    &alt_lnl, SIGFIGS, d->lb->data[0], d->ub->data[0],
                    logf, NULL, NULL);
      /* turns out to be faster (roughly 15% in limited experiments)
         to use numerical rather than exact derivatives */

      alt_lnl *= -1;
      this_scale = d->params->data[0];

      delta_lnl = alt_lnl - null_lnl;
      if (delta_lnl <= -0.01)
	die("ERROR col_lrts: delta_lnl = %e < -0.01\n", delta_lnl);
      if (delta_lnl < 0) delta_lnl = 0;
    } /* end estimation of delta_lnl */

    /* compute p-vals via chi-sq */
    if (tuple_pvals != NULL) {
      if (mode == NNEUT || mode == CONACC)
        tuple_pvals[i] = chisq_cdf(2*delta_lnl, 1, FALSE);
      else
        tuple_pvals[i] = half_chisq_cdf(2*delta_lnl, 1, FALSE);
        /* assumes 50:50 mix of chisq and point mass at zero, due to
           bounding of param */

      if (tuple_pvals[i] < 1e-20)
        tuple_pvals[i] = 1e-20;
      /* approx limit of eval of tail prob; pvals of 0 cause problems */

      if (mode == CONACC && this_scale > 1)
          tuple_pvals[i] *= -1; /* mark as acceleration */
    }

    /* store scales and log likelihood ratios if necessary */
    if (tuple_scales != NULL) tuple_scales[i] = this_scale;
    if (tuple_llrs != NULL) tuple_llrs[i] = delta_lnl;
  }

  col_free_fit_data(d);
}

/* Subtree version of LRT */
void col_lrts_sub(TreeModel *mod, MSA *msa, mode_type mode,
                  double *tuple_pvals, double *tuple_null_scales,
                  double *tuple_scales, double *tuple_sub_scales,
                  double *tuple_llrs, FILE *logf) {
  int i;
  ColFitData *d, *d2;
  double null_lnl, alt_lnl, delta_lnl;
  TreeModel *modcpy;
  List *inside=NULL, *outside=NULL;

  modcpy = tm_create_copy(mod);   /* need separate copy of tree model
                                     with different internal scaling
                                     data for supertree/subtree case */
  modcpy->subtree_root = NULL;

  /* init ColFitData -- one for null model, one for alt */
  d = col_init_fit_data(modcpy, msa, ALL, NNEUT, FALSE);
  d2 = col_init_fit_data(mod, msa, SUBTREE, mode, FALSE);
                                /* mod has the subtree info, modcpy
                                   does not */

  /* prepare lists of leaves inside and outside root, for use in
     checking for informative substitutions */
  if (mod->subtree_root != NULL) {
    inside = lst_new_ptr(mod->tree->nnodes);
    outside = lst_new_ptr(mod->tree->nnodes);
    tr_partition_leaves(mod->tree, mod->subtree_root, inside, outside);
  }

  /* iterate through column tuples */
  for (i = 0; i < msa->ss->ntuples; i++) {
    checkInterruptN(i, 100);

    /* first check for informative substitution data in column; if none,
       don't waste time computing likeihoods */
    if (!col_has_data_sub(mod, msa, i, inside, outside)) {
      delta_lnl = 0;
      d->params->data[0] = d2->params->data[0] = d2->params->data[1] = 1;
    }

    else {
      /* compute log likelihoods under null and alt hypotheses */
      d->tupleidx = i;
      vec_set(d->params, 0, d->init_scale);
      opt_newton_1d(col_likelihood_wrapper_1d, &d->params->data[0], d,
                    &null_lnl, SIGFIGS, d->lb->data[0], d->ub->data[0],
                    logf, NULL, NULL);

      //      opt_bfgs(col_likelihood_wrapper, d->params, d, &null_lnl, d->lb,
      //	       d->ub, logf, NULL, OPT_HIGH_PREC, NULL, NULL);

      /* turns out to be faster (roughly 15% in limited experiments)
         to use numerical rather than exact derivatives */
      null_lnl *= -1;

      d2->tupleidx = i;
      vec_set(d2->params, 0, max(0.05, d->params->data[0]));
      /* init to previous estimate to save time, but don't init to
         value at boundary */
      vec_set(d2->params, 1, d2->init_scale_sub);

      if (opt_bfgs(col_likelihood_wrapper, d2->params, d2, &alt_lnl, d2->lb,
                   d2->ub, logf, NULL, OPT_HIGH_PREC, NULL, NULL) != 0)
        ;                         /* do nothing; nonzero exit typically
                                     occurs when max iterations is
                                     reached; a warning is printed to
                                     the log */
      alt_lnl *= -1;

      delta_lnl = alt_lnl - null_lnl;
      if (delta_lnl <= -0.1)
	die("ERROR col_lrts_sub: delta_lnl = %e <= -0.1\n", delta_lnl);
      if (delta_lnl < 0) delta_lnl = 0;
    }

    /* compute p-vals via chi-sq */
    if (tuple_pvals != NULL) {
      if (mode == NNEUT || mode == CONACC)
        tuple_pvals[i] = chisq_cdf(2*delta_lnl, 1, FALSE);
      else
        tuple_pvals[i] = half_chisq_cdf(2*delta_lnl, 1, FALSE);
        /* assumes 50:50 mix of chisq and point mass at zero, due to
           bounding of param */

      if (tuple_pvals[i] < 1e-20)
        tuple_pvals[i] = 1e-20;
      /* approx limit of eval of tail prob; pvals of 0 cause problems */

      if (mode == CONACC && d2->params->data[1] > 1)
        tuple_pvals[i] *= -1;    /* mark as acceleration */
    }

    /* store scales and log likelihood ratios if necessary */
    if (tuple_null_scales != NULL)
      tuple_null_scales[i] = d->params->data[0];
    if (tuple_scales != NULL)
      tuple_scales[i] = d2->params->data[0];
    if (tuple_sub_scales != NULL)
      tuple_sub_scales[i] = d2->params->data[1];
    if (tuple_llrs != NULL)
      tuple_llrs[i] = delta_lnl;
  }

  col_free_fit_data(d);
  col_free_fit_data(d2);
  modcpy->estimate_branchlens = TM_BRANCHLENS_ALL;
                                /* have to revert for tm_free to work
                                   correctly */
  tm_free(modcpy);
  if (inside != NULL) lst_free(inside);
  if (outside != NULL) lst_free(outside);
}

/* Score test */
void col_score_tests(TreeModel *mod, MSA *msa, mode_type mode,
                     double *tuple_pvals, double *tuple_derivs,
                     double *tuple_teststats) {
  int i;
  ColFitData *d;
  double first_deriv, teststat, fim;

  /* init ColFitData */
  d = col_init_fit_data(mod, msa, ALL, NNEUT, FALSE);

  /* precompute FIM */
  fim = col_estimate_fim(mod);

  if (fim < 0)
    die("ERROR: negative fisher information in col_score_tests\n");

  /* iterate through column tuples */
  for (i = 0; i < msa->ss->ntuples; i++) {
    checkInterruptN(i, 1000);

    /* first check for actual substitution data in column; if none,
       don't waste time computing score */
    if (!col_has_data(mod, msa, i)) {
      first_deriv = 0;
      teststat = 0;
    }

    else {
      d->tupleidx = i;

      col_scale_derivs(d, &first_deriv, NULL, d->fels_scratch);

      teststat = first_deriv*first_deriv / fim;

      if ((mode == ACC && first_deriv < 0) ||
          (mode == CON && first_deriv > 0))
        teststat = 0;             /* derivative points toward boundary;
                                     truncate at 0 */
    }

    if (tuple_pvals != NULL) {
      if (mode == NNEUT || mode == CONACC)
        tuple_pvals[i] = chisq_cdf(teststat, 1, FALSE);
      else
        tuple_pvals[i] = half_chisq_cdf(teststat, 1, FALSE);
        /* assumes 50:50 mix of chisq and point mass at zero */

      if (tuple_pvals[i] < 1e-20)
        tuple_pvals[i] = 1e-20;
      /* approx limit of eval of tail prob; pvals of 0 cause problems */

      if (mode == CONACC && first_deriv > 0)
        tuple_pvals[i] *= -1; /* mark as acceleration */
    }

    /* store scales and log likelihood ratios if necessary */
    if (tuple_derivs != NULL) tuple_derivs[i] = first_deriv;
    if (tuple_teststats != NULL) tuple_teststats[i] = teststat;
  }

  col_free_fit_data(d);
}

/* Subtree version of score test */
void col_score_tests_sub(TreeModel *mod, MSA *msa, mode_type mode,
                         double *tuple_pvals, double *tuple_null_scales,
                         double *tuple_derivs, double *tuple_sub_derivs,
                         double *tuple_teststats, FILE *logf) {
  int i;
  ColFitData *d, *d2;
  Vector *grad = vec_new(2);
  Matrix *fim;
  double lnl, teststat;
  FimGrid *grid;
  List *inside=NULL, *outside=NULL;
  TreeModel *modcpy = tm_create_copy(mod); /* need separate copy of tree model
                                              with different internal scaling
                                              data for supertree/subtree case */
  modcpy->subtree_root = NULL;

  /* init ColFitData -- one for null model, one for alt */
  d = col_init_fit_data(modcpy, msa, ALL, NNEUT, FALSE);
  d2 = col_init_fit_data(mod, msa, SUBTREE, NNEUT, FALSE);
                                /* mod has the subtree info, modcpy
                                   does not */

  /* precompute Fisher information matrices for a grid of scale values */
  grid = col_fim_grid_sub(mod);

  /* prepare lists of leaves inside and outside root, for use in
     checking for informative substitutions */
  if (mod->subtree_root != NULL) {
    inside = lst_new_ptr(mod->tree->nnodes);
    outside = lst_new_ptr(mod->tree->nnodes);
    tr_partition_leaves(mod->tree, mod->subtree_root, inside, outside);
  }

  /* iterate through column tuples */
  for (i = 0; i < msa->ss->ntuples; i++) {
    checkInterruptN(i, 100);

    /* first check for informative substitution data in column; if none,
       don't waste time computing score */
    if (!col_has_data_sub(mod, msa, i, inside, outside)) {
      teststat = 0;
      vec_zero(grad);
      d->params->data[0] = 1.0;
    }

    else {
      d->tupleidx = i;
      vec_set(d->params, 0, d->init_scale);

      opt_newton_1d(col_likelihood_wrapper_1d, &d->params->data[0], d,
                    &lnl, SIGFIGS, d->lb->data[0], d->ub->data[0],
                    logf, NULL, NULL);
      /* turns out to be faster (roughly 15% in limited experiments)
         to use numerical rather than exact derivatives */

      d2->tupleidx = i;
      d2->mod->scale = d->params->data[0];
      d2->mod->scale_sub = 1;
      tm_set_subst_matrices(d2->mod);
      col_scale_derivs_subtree(d2, grad, NULL, d2->fels_scratch);

      fim = col_get_fim_sub(grid, d2->mod->scale);

      teststat = grad->data[1]*grad->data[1] /
        (fim->data[1][1] - fim->data[0][1]*fim->data[1][0]/fim->data[0][0]);

      if (teststat < 0) {
        fprintf(stderr, "WARNING: teststat < 0 (%f\t%f\t%f\t%f\t%f\t%f)\n",
                teststat, fim->data[0][0], fim->data[0][1],
                fim->data[1][0], fim->data[1][1],
                fim->data[0][1]*fim->data[1][0]/fim->data[0][0]);
        teststat = 0;
      }
      mat_free(fim);

      if ((mode == ACC && grad->data[1] < 0) ||
          (mode == CON && grad->data[1] > 0))
        teststat = 0;             /* derivative points toward boundary;
                                     truncate at 0 */
    }

    if (tuple_pvals != NULL) {
      if (mode == NNEUT || mode == CONACC)
        tuple_pvals[i] = chisq_cdf(teststat, 1, FALSE);
      else
        tuple_pvals[i] = half_chisq_cdf(teststat, 1, FALSE);
      /* assumes 50:50 mix of chisq and point mass at zero */

      if (tuple_pvals[i] < 1e-20)
        tuple_pvals[i] = 1e-20;
      /* approx limit of eval of tail prob; pvals of 0 cause problems */

      if (mode == CONACC && grad->data[1] > 0)
        tuple_pvals[i] *= -1; /* mark as acceleration */
    }

    /* store scales and log likelihood ratios if necessary */
    if (tuple_null_scales != NULL) tuple_null_scales[i] = d->params->data[0];
    if (tuple_derivs != NULL) tuple_derivs[i] = grad->data[0];
    if (tuple_sub_derivs != NULL) tuple_sub_derivs[i] = grad->data[1];
    if (tuple_teststats != NULL) tuple_teststats[i] = teststat;
  }

  col_free_fit_data(d);
  col_free_fit_data(d2);
  vec_free(grad);
  modcpy->estimate_branchlens = TM_BRANCHLENS_ALL;
                                /* have to revert for tm_free to work
                                   correctly */
  tm_free(modcpy);
  if (inside != NULL) lst_free(inside);
  if (outside != NULL) lst_free(outside);
  col_free_fim_grid(grid);
}

/* Create object with metadata and scratch memory for fitting scale
   factors */
ColFitData *col_init_fit_data(TreeModel *mod, MSA *msa, scale_type stype,
                              mode_type mode, int second_derivs) {
  ColFitData *d = smalloc(sizeof(ColFitData));
  int size = mod->rate_matrix->size, nrcats = mod->nratecats,
    nnodes = mod->tree->nnodes;
  int nid, rcat, i, j, dim;

  d->mod = mod;
  d->msa = msa;
  d->stype = stype;
  d->mode = mode;
  d->second_derivs = second_derivs;
  d->tupleidx = -1;         /* will be set as needed */

  d->mod->estimate_branchlens = TM_SCALE_ONLY;
  if (stype == SUBTREE) {
    if (!(mod->subtree_root != NULL || mod->in_subtree!=NULL))
      die("ERROR col_init_fit_data: mod->subtree_root or mod->in_subtree must not be NULL in SUBTREE mode\n");
  }
  else {
    if (mod->subtree_root != NULL)
      die("ERROR col_init_fit_data: mod->subtree_root must be NULL if not in subtree mode\n");
    mod->scale_sub = 1;
  }
  if (mod->msa_seq_idx == NULL)
    tm_build_seq_idx(mod, msa);
  tm_set_subst_matrices(mod);

  dim = (stype == ALL ? 1 : 2);
  d->params = vec_new(dim);
  d->lb = vec_new(dim);
  d->ub = vec_new(dim);
  vec_set(d->lb, 0, 0);
  vec_set(d->ub, 0, INFTY);
  d->init_scale = d->init_scale_sub = 1;

  if (stype == ALL) {
    if (mode == CON) {
      vec_set(d->ub, 0, 1);
      d->init_scale = 0.9;      /* don't start on boundary but avoid
                                   strong initialization bias */
    }
    else if (mode == ACC) {
      vec_set(d->lb, 0, 1);
      d->init_scale = 1.1;      /* don't start on boundary but avoid
                                   strong initialization bias */
    }
  }
  else {                        /* stype == SUBTREE */
    vec_set(d->lb, 1, 0);
    vec_set(d->ub, 1, INFTY);
    vec_set(d->lb, 0, 1e-6);    /* can't let this param go quite to
                                   zero, because other param becomes
                                   undefined */

    if (mode == CON) {
      vec_set(d->ub, 1, 1);
      d->init_scale_sub = 0.9;  /* don't start on boundary but avoid
                                   strong initialization bias */
    }
    else if (mode == ACC) {
      vec_set(d->lb, 1, 1);
      d->init_scale_sub = 1.1;  /* don't start on boundary but avoid
                                   strong initialization bias */
    }
  }


  d->PP = smalloc(nnodes * sizeof(void*));
  d->PPP = smalloc(nnodes * sizeof(void*));
  d->QQ = smalloc(nnodes * sizeof(void*));
  d->QQQ = smalloc(nnodes * sizeof(void*));
  d->RRR = smalloc(nnodes * sizeof(void*));

  for (nid = 0; nid < nnodes; nid++) {
    d->PP[nid] = smalloc(nrcats * sizeof(void*));
    d->PPP[nid] = smalloc(nrcats * sizeof(void*));
    d->QQ[nid] = smalloc(nrcats * sizeof(void*));
    d->QQQ[nid] = smalloc(nrcats * sizeof(void*));
    d->RRR[nid] = smalloc(nrcats * sizeof(void*));

    for (rcat = 0; rcat < nrcats; rcat++) {
      d->PP[nid][rcat] = mat_new(size, size);
      d->PPP[nid][rcat] = mat_new(size, size);
      d->QQ[nid][rcat] = mat_new(size, size);
      d->QQQ[nid][rcat] = mat_new(size, size);
      d->RRR[nid][rcat] = mat_new(size, size);

      mat_zero(d->PP[nid][rcat]);
      mat_zero(d->PPP[nid][rcat]);
      mat_zero(d->QQ[nid][rcat]);
      mat_zero(d->QQQ[nid][rcat]);
      mat_zero(d->RRR[nid][rcat]);
    }
  }

  d->expdiag_z = zvec_new(size);
  d->expdiag_r = vec_new(size);
  d->nfels_scratch = stype == SUBTREE ? 6 : 3;
  d->fels_scratch = smalloc(d->nfels_scratch * sizeof(void*));
  for (i = 0; i < d->nfels_scratch; i++) {
    d->fels_scratch[i] = smalloc(size * sizeof(void*));
    for (j = 0; j < size; j++)
      d->fels_scratch[i][j] = smalloc((nnodes+1) * sizeof(double));
  }
  d->mat_scratch_z = zmat_new(size, size);
  d->vec_scratch1_z = zvec_new(size);
  d->vec_scratch2_z = zvec_new(size);
  d->vec_scratch1_r = vec_new(size);
  d->vec_scratch2_r = vec_new(size);
  return d;
}

/* Free metadata and memory for fitting scale factors */
void col_free_fit_data(ColFitData *d) {
  int nid, rcat, i, j;

  vec_free(d->params);
  vec_free(d->lb);
  vec_free(d->ub);

  for (nid = 0; nid < d->mod->tree->nnodes; nid++) {
    for (rcat = 0; rcat < d->mod->nratecats; rcat++) {
      mat_free(d->PP[nid][rcat]);
      mat_free(d->PPP[nid][rcat]);
      mat_free(d->QQ[nid][rcat]);
      mat_free(d->QQQ[nid][rcat]);
      mat_free(d->RRR[nid][rcat]);
    }
    sfree(d->PP[nid]);
    sfree(d->PPP[nid]);
    sfree(d->QQ[nid]);
    sfree(d->QQQ[nid]);
    sfree(d->RRR[nid]);
  }
  sfree(d->PP);
  sfree(d->PPP);
  sfree(d->QQ);
  sfree(d->QQQ);
  sfree(d->RRR);

  zvec_free(d->expdiag_z);
  vec_free(d->expdiag_r);
  for (i = 0; i < d->nfels_scratch; i++) {
    for (j = 0; j < d->mod->rate_matrix->size; j++)
      sfree(d->fels_scratch[i][j]);
    sfree(d->fels_scratch[i]);
  }
  sfree(d->fels_scratch);
  zmat_free(d->mat_scratch_z);
  zvec_free(d->vec_scratch1_z);
  zvec_free(d->vec_scratch2_z);
  vec_free(d->vec_scratch1_r);
  vec_free(d->vec_scratch2_r);

  sfree(d);
}

/* Perform a GERP-like computation for each tuple.  Computes expected
   number of subst. under neutrality (tuple_nneut), expected number
   after rescaling by ML (tuple_nobs), expected number of rejected
   substitutions (tuple_nrejected), and number of species with data
   (tuple_nspecies).  If any arrays are NULL, values will not be
   retained.  Gaps and missing data are handled by working with
   induced subtree.  */
void col_gerp(TreeModel *mod, MSA *msa, mode_type mode, double *tuple_nneut,
              double *tuple_nobs, double *tuple_nrejected,
              double *tuple_nspec, FILE *logf) {
  int i, j, nspec = 0;
  double nneut, scale, lnl;
  int *has_data = smalloc(mod->tree->nnodes * sizeof(int));
  ColFitData *d;

  /* init ColFitData */
  d = col_init_fit_data(mod, msa, ALL, NNEUT, FALSE);

  /* iterate through column tuples */
  for (i = 0; i < msa->ss->ntuples;i++) {
    checkInterruptN(i, 1000);
    col_find_missing_branches(mod, msa, i, has_data, &nspec);

    if (nspec < 3)
      nneut = scale = 0;
    else {
      vec_set(d->params, 0, d->init_scale);
      d->tupleidx = i;

      opt_newton_1d(col_likelihood_wrapper_1d, &d->params->data[0], d,
                    &lnl, SIGFIGS, d->lb->data[0], d->ub->data[0],
                    logf, NULL, NULL);
      /* turns out to be faster (roughly 15% in limited experiments)
         to use numerical rather than exact derivatives */

      scale = d->params->data[0];
      for (j = 1, nneut = 0; j < mod->tree->nnodes; j++)  /* node 0 is root */
        if (has_data[j])
          nneut += ((TreeNode*)lst_get_ptr(mod->tree->nodes, j))->dparent;
    }

    if (tuple_nspec != NULL) tuple_nspec[i] = (double)nspec;
    if (tuple_nneut != NULL) tuple_nneut[i] = nneut;
    if (tuple_nobs != NULL) tuple_nobs[i] = scale * nneut;
    if (tuple_nrejected != NULL) {
      tuple_nrejected[i] = nneut * (1 - scale);
      if (mode == ACC) tuple_nrejected[i] *= -1;
      else if (mode == NNEUT) tuple_nrejected[i] = fabs(tuple_nrejected[i]);
    }
  }
  col_free_fit_data(d);
  sfree(has_data);
}

/* Identify branches wrt which a given column tuple is uninformative,
   in the sense that all leaves beneath these branches having missing
   data.  Will set (preallocated) array has_data[i] = I(branch above
   node i is informative).  Will also set *nspec equal to number of
   leaves that have data. */
void col_find_missing_branches(TreeModel *mod, MSA *msa, int tupleidx,
                               int *has_data, int *nspec) {
  int i;
  List *traversal = tr_postorder(mod->tree);
  *nspec = 0;
  for (i = 0; i < lst_size(traversal); i++) {
    TreeNode *n = lst_get_ptr(traversal, i);
    if (!((n->lchild == NULL && n->rchild == NULL) ||
	  (n->lchild != NULL && n->rchild != NULL)))
      die("ERROR: col_find_missing_branches: lchild and rchild must be either both NULL or both non-NULL\n");
    if (n->parent == NULL)      /* root */
      has_data[n->id] = FALSE;
    else if (n->lchild == NULL) {    /* leaf */
      if (mod->rate_matrix->
          inv_states[(int)ss_get_char_tuple(msa, tupleidx,
                                            mod->msa_seq_idx[n->id], 0)] >= 0) {
        has_data[n->id] = TRUE;
        (*nspec)++;
      }
      else
        has_data[n->id] = FALSE;
    }
    else {                      /* non-root ancestral node */
      if (has_data[n->lchild->id] || has_data[n->rchild->id])
        has_data[n->id] = TRUE;
      else
        has_data[n->id] = FALSE;
    }
  }
}

/* Numerically compute first and second derivatives of single-column
   log likelihood function.  For debugging */
void col_scale_derivs_num(ColFitData *d, double *first_deriv,
                          double *second_deriv) {
  double lnl1, lnl2, lnl3;
  double orig_scale = d->mod->scale;
  d->mod->scale += (2*DERIV_EPSILON);
  tm_set_subst_matrices(d->mod);
  lnl1 = col_compute_log_likelihood(d->mod, d->msa, d->tupleidx,
                                    d->fels_scratch[0]);
  d->mod->scale = orig_scale + DERIV_EPSILON;
  tm_set_subst_matrices(d->mod);
  lnl2 = col_compute_log_likelihood(d->mod, d->msa, d->tupleidx,
                                    d->fels_scratch[0]);
  d->mod->scale = orig_scale;
  tm_set_subst_matrices(d->mod);
  lnl3 = col_compute_log_likelihood(d->mod, d->msa, d->tupleidx,
                                    d->fels_scratch[0]);
  *first_deriv = (lnl2 - lnl3) / DERIV_EPSILON;
  *second_deriv = (lnl1 - 2*lnl2 + lnl3) / (DERIV_EPSILON * DERIV_EPSILON);
}

/* Numerically compute first and second derivatives of single-column
   log likelihood function.  For debugging */
void col_scale_derivs_subtree_num(ColFitData *d, Vector *gradient,
                                  Matrix *hessian) {
  double lnl00, lnl11, lnl0, lnl1, lnl01, lnl;
  double orig_scale = d->mod->scale, orig_scale_sub = d->mod->scale_sub;

  d->mod->scale += (2*DERIV_EPSILON);
  tm_set_subst_matrices(d->mod);
  lnl00 = col_compute_log_likelihood(d->mod, d->msa, d->tupleidx,
                                    d->fels_scratch[0]);

  d->mod->scale = orig_scale + DERIV_EPSILON;
  tm_set_subst_matrices(d->mod);
  lnl0 = col_compute_log_likelihood(d->mod, d->msa, d->tupleidx,
                                    d->fels_scratch[0]);

  d->mod->scale = orig_scale + DERIV_EPSILON;
  d->mod->scale_sub = orig_scale_sub + DERIV_EPSILON;
  tm_set_subst_matrices(d->mod);
  lnl01 = col_compute_log_likelihood(d->mod, d->msa, d->tupleidx,
                                    d->fels_scratch[0]);

  d->mod->scale = orig_scale;
  d->mod->scale_sub = orig_scale_sub + DERIV_EPSILON;
  tm_set_subst_matrices(d->mod);
  lnl1 = col_compute_log_likelihood(d->mod, d->msa, d->tupleidx,
                                    d->fels_scratch[0]);

  d->mod->scale = orig_scale;
  d->mod->scale_sub = orig_scale_sub + (2*DERIV_EPSILON);
  tm_set_subst_matrices(d->mod);
  lnl11 = col_compute_log_likelihood(d->mod, d->msa, d->tupleidx,
                                    d->fels_scratch[0]);

  d->mod->scale = orig_scale;
  d->mod->scale_sub = orig_scale_sub;
  tm_set_subst_matrices(d->mod);
  lnl = col_compute_log_likelihood(d->mod, d->msa, d->tupleidx,
                                    d->fels_scratch[0]);

  gradient->data[0] = (lnl0 - lnl) / DERIV_EPSILON;
  gradient->data[1] = (lnl1 - lnl) / DERIV_EPSILON;

  hessian->data[0][0] = (lnl00 - 2*lnl0 + lnl) /
    (DERIV_EPSILON * DERIV_EPSILON);
  hessian->data[1][1] = (lnl11 - 2*lnl1 + lnl) /
    (DERIV_EPSILON * DERIV_EPSILON);
  hessian->data[0][1] = hessian->data[1][0] =
    (lnl01 - lnl0 - lnl1 + lnl) /
    (DERIV_EPSILON * DERIV_EPSILON);
}

/* Estimate 2x2 Fisher Information Matrix (expected value of the
   negative Hessian) for the subtree case, based on a particular value
   of the scale parameter (set in calling code).  Estimation is done
   by sampling. */
Matrix *col_estimate_fim_sub(TreeModel *mod) {
  Vector *grad = vec_new(2);
  Matrix *hessian = mat_new(2, 2), *fim = mat_new(2, 2);
  MSA *msa = tm_generate_msa(NSAMPLES_FIM, NULL, &mod, NULL);
  ColFitData *d = col_init_fit_data(mod, msa, SUBTREE, NNEUT, TRUE);
  int i;

  ss_from_msas(msa, 1, TRUE, NULL, NULL, NULL, -1, 0);
  mat_zero(fim);
  for (i = 0; i < msa->ss->ntuples; i++) {
    d->tupleidx = i;
    col_scale_derivs_subtree(d, grad, hessian, d->fels_scratch);
    mat_scale(hessian, -1 * msa->ss->counts[i]);
                                /* now (observed) Fisher matrix
                                   weighted by count */
    mat_plus_eq(fim, hessian);   /* add to running total */
  }
  mat_scale(fim, 1.0/NSAMPLES_FIM);   /* convert total to sample mean */

  msa_free(msa);
  col_free_fit_data(d);
  vec_free(grad);
  mat_free(hessian);
  return (fim);
}

/* Precompute estimates of FIM for a grid of possible scale params
   (subtree case) */
FimGrid *col_fim_grid_sub(TreeModel *mod) {
  int i;
  FimGrid *g = smalloc(sizeof(FimGrid));

  g->ngrid1 = (int)(1.0/GRIDSIZE1);
  g->ngrid2 = (int)((1.0 * GRIDMAXLOG / GRIDSIZE2) + 1);
  g->ngrid = g->ngrid1 + g->ngrid2;
  g->scales = smalloc(g->ngrid * sizeof(double));

  mod->scale_sub = 1;

  for (i = 0; i < g->ngrid1; i++)
    g->scales[i] = i * GRIDSIZE1;

  for (i = 0; i < g->ngrid2; i++)
    g->scales[g->ngrid1 + i] = exp(i * GRIDSIZE2);

  g->fim = smalloc(g->ngrid * sizeof(void*));
  for (i = 0; i < g->ngrid; i++) {
    mod->scale = g->scales[i];
    tm_set_subst_matrices(mod);
    g->fim[i] = col_estimate_fim_sub(mod);
  }

  return g;
}

/* free FimGrid object */
void col_free_fim_grid(FimGrid *g) {
  int i;
  for (i = 0; i < g->ngrid; i++)
    mat_free(g->fim[i]);
  sfree(g->fim);
  sfree(g->scales);
}

/* Estimate scale Fisher Information Matrix for the non-subtree case.
   This version does not depend on any free parameters, so no grid is
   required.  Estimation is done by sampling, as above */
double col_estimate_fim(TreeModel *mod) {
  double deriv1, deriv2, retval = 0;
  MSA *msa = tm_generate_msa(NSAMPLES_FIM, NULL, &mod, NULL);
  ColFitData *d = col_init_fit_data(mod, msa, ALL, NNEUT, FALSE);
  int i;

  ss_from_msas(msa, 1, TRUE, NULL, NULL, NULL, -1, 0);

  for (i = 0; i < msa->ss->ntuples; i++) {
    d->tupleidx = i;
    col_scale_derivs(d, &deriv1, &deriv2, d->fels_scratch);
    retval += (-deriv2 * msa->ss->counts[i]); /* add (observed) Fisher matrix
                                                 weighted by count */
  }
  retval /= NSAMPLES_FIM;   /* convert total to sample mean */

  msa_free(msa);
  col_free_fit_data(d);
  return (retval);
}

/* Retrieve estimated FIM for given scale function; uses linear
   interpolation from precomputed grid */
Matrix *col_get_fim_sub(FimGrid *g, double scale) {
  int idx;
  double frac;
  Matrix *retval;

  if (scale < 0)
    die("ERROR col_get_fix_sub: scale should be >= 0 but is %e\n", scale);

  if (scale < 1)
    idx = (int)floor(scale / GRIDSIZE1);
  else
    idx = g->ngrid1 + (int)floor(log(scale) / GRIDSIZE2);

  if (idx >= g->ngrid - 1)
    retval = mat_create_copy(g->fim[g->ngrid - 1]);
                                /* just use last one in this case */

  else {
    if (!(g->scales[idx] <= scale && g->scales[idx+1] > scale))
      die("ERROR col_get_fim_sub: g->scales[%i]=%e should be <= %e and g->scales[%i+1]=%e should be > %e\n", idx, g->scales[idx], scale, idx+1, g->scales[idx+1], scale);

    frac = (scale - g->scales[idx]) / (g->scales[idx+1] - g->scales[idx]);

    if (frac < 0.1)
      retval = mat_create_copy(g->fim[idx]);
    else if (frac > 0.9)
      retval = mat_create_copy(g->fim[idx+1]);
    else {  /* interpolate */
      retval = mat_new(2,2);
      mat_linear_comb(retval, g->fim[idx], frac, g->fim[idx+1], 1-frac);
    }
  }
  return retval;
}

/* returns TRUE if column has two or more actual bases (not gaps or
   missing data), otherwise returns FALSE */
int col_has_data(TreeModel *mod, MSA *msa, int tupleidx) {
  int i, nbases = 0;
  for (i = 0; i < msa->nseqs && nbases < 2; i++) {
    int state = mod->rate_matrix->
      inv_states[(int)ss_get_char_tuple(msa, tupleidx, i, 0)];
    if (state >= 0)
      nbases++;
  }
  return(nbases >= 2);
}

/* returns TRUE if column has at least one base in the subtree of
   interest, at least one in the supertree of interest, and at least
   three bases total (the minimum required for a meaningful subtree
   test), otherwise returns FALSE.
   If inside and outside are both NULL, then returns TRUE if there
   are at least three bases total.
 */
int col_has_data_sub(TreeModel *mod, MSA *msa, int tupleidx, List *inside,
                     List *outside) {
  int i, nbases = 0, state;
  TreeNode *n;

  if (inside == NULL && outside == NULL) {
    for (i=0; i<mod->tree->nnodes; i++) {
      n = lst_get_ptr(mod->tree->nodes, i);
      if (n->lchild == NULL) {
	state = mod->rate_matrix->
	  inv_states[(int)ss_get_char_tuple(msa, tupleidx,
					    mod->msa_seq_idx[n->id], 0)];
	if (state >=0) nbases++;
	if (nbases==2) return TRUE;
      }
    }
    return FALSE;
  }
  for (i = 0; i < lst_size(inside) && nbases < 2; i++) {
    n = lst_get_ptr(inside, i);
    state = mod->rate_matrix->
      inv_states[(int)ss_get_char_tuple(msa, tupleidx,
                                        mod->msa_seq_idx[n->id], 0)];
    if (state >= 0) nbases++;
  }

  if (nbases == 0)              /* no leaves in subtree */
    return FALSE;
  /* NOTE: nbases at most two here */

  for (i = 0; i < lst_size(outside) && nbases < 3; i++) {
    n = lst_get_ptr(outside, i);
    state = mod->rate_matrix->
      inv_states[(int)ss_get_char_tuple(msa, tupleidx,
                                        mod->msa_seq_idx[n->id], 0)];
    if (state >= 0) nbases++;
  }

  if (nbases == 3)              /* has to be at least one in each
                                   partition and a third from one of
                                   the two partitions */
    return TRUE;

  return FALSE;
}
