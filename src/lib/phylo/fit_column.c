/* $Id: fit_column.c,v 1.8 2008-03-18 14:30:57 acs Exp $
   Written by Adam Siepel, 2008
*/

/* Functions to compute likelihoods for individual alignment columns
   and estimate column-by-column scale factors by maximum likelhood.
   For use with single-base LRTs, score tests, phyloP, etc. */

#include <fit_column.h>
#include <sufficient_stats.h>

/* Compute and return the likelihood of a tree model with respect to a
   single column tuple in an alignment.  This is a pared-down version
   of tl_compute_log_likelihood for use in estimation of base-by-base
   scale factors.  It assumes a 0th order model, leaf-to-sequence
   mapping already available, prob matrices computed, sufficient stats
   already available.  Note that this function uses natural log rather
   than log2.  This function does allow for rate variation. */
double col_compute_log_likelihood(TreeModel *mod, MSA *msa, int tupleidx,
                                  double **scratch) {

  int i, j, k, nodeidx, rcat;
  int nstates = mod->rate_matrix->size;
  TreeNode *n;
  double total_prob = 0;
  List *traversal = tr_postorder(mod->tree);  
  double **pL = NULL;

  assert(msa->ss->tuple_size == 1);
  assert(mod->order == 0);
  assert(mod->allow_gaps == TRUE);

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
    for (j = 0; j < nstates; j++) free(pL[j]);
    free(pL);
  }

  return(log(total_prob));
}

/* Compute 1st and 2nd derivs wrt scale params of substitution
   matrices for each branch of the tree (and each rate category).
   These are used in the recursive computation of derivatives of the
   column likelihoods */
void col_scale_derivs_subst(ColFitData *d) {
  Zmatrix *S = d->mod->rate_matrix->evec_matrix, 
    *Sinv = d->mod->rate_matrix->evec_matrix_inv;
  MarkovMatrix *Q = d->mod->rate_matrix;
  int size = Q->size;
  int rcat, nid, i;

  assert(S != NULL && Sinv != NULL);

  for (rcat = 0; rcat < d->mod->nratecats; rcat++) {
    for (nid = 1; nid < d->mod->tree->nnodes; nid++) { /* skip root */

      double t = ((TreeNode*)lst_get_ptr(d->mod->tree->nodes, nid))->dparent;
      double l1 = d->mod->scale;
      double l2 = (d->stype == SUBTREE && d->mod->in_subtree[nid] ? 
                   d->mod->scale_sub : 1);

      /* set up exponentiated diagonal matrix */
      for (i = 0; i < size; i++) {
        double r = t * l1 * l2 * d->mod->rK[rcat];
        zvec_set(d->expdiag, i, z_exp(z_mul_real(zvec_get(Q->evals, i), r)));
      }

      /* PP */
      zvec_copy(d->vec_scratch1, Q->evals);
      zvec_scale(d->vec_scratch1, t * l2);
      zvec_had_prod(d->vec_scratch2, d->vec_scratch1, d->expdiag);
      zmat_mult_real_diag(d->PP[nid][rcat], S, d->vec_scratch2, Sinv, 
                          d->mat_scratch);

      if (d->second_derivs) {
        /* PPP */
        zvec_had_prod(d->vec_scratch2, d->vec_scratch1, d->vec_scratch1);
        zvec_had_prod(d->vec_scratch1, d->vec_scratch2, d->expdiag);
        zmat_mult_real_diag(d->PPP[nid][rcat], S, d->vec_scratch1, Sinv, 
                            d->mat_scratch);
      }

      if (d->stype == SUBTREE && d->mod->in_subtree[nid]) {
        /* if not in subtree, leave all of these equal to 0 (as
           initialized) */

        /* QQ */
        zvec_copy(d->vec_scratch1, Q->evals);
        zvec_scale(d->vec_scratch1, t * l1);
        zvec_had_prod(d->vec_scratch2, d->vec_scratch1, d->expdiag);
        zmat_mult_real_diag(d->QQ[nid][rcat], S, d->vec_scratch2, Sinv, 
                            d->mat_scratch);

        if (d->second_derivs) {
          /* QQQ */
          zvec_had_prod(d->vec_scratch2, d->vec_scratch1, d->vec_scratch1);
          zvec_had_prod(d->vec_scratch1, d->vec_scratch2, d->expdiag);
          zmat_mult_real_diag(d->QQQ[nid][rcat], S, d->vec_scratch1, Sinv, 
                              d->mat_scratch);

          /* RRR */        
          zvec_copy(d->vec_scratch1, Q->evals);
          zvec_scale(d->vec_scratch1, t);
          zvec_had_prod(d->vec_scratch2, Q->evals, Q->evals);
          zvec_scale(d->vec_scratch2, t * l1 * l2);
          zvec_plus_eq(d->vec_scratch2, d->vec_scratch1);
          zvec_had_prod(d->vec_scratch1, d->vec_scratch2, d->expdiag);      
          zmat_mult_real_diag(d->RRR[nid][rcat], S, d->vec_scratch1, Sinv, 
                              d->mat_scratch);
        }
      }
    }
  }
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
  double **L;                   /* partial likelihoods */
  double **LL;                  /* 1st deriv of partial likelihoods wrt
                                   scale param */ 
  double **LLL;                 /* 2nd deriv of partial likelihoods
                                   wrt scale param */

  assert(d->msa->ss->tuple_size == 1);
  assert(d->mod->order == 0);
  assert(d->mod->allow_gaps == TRUE);

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
    *second_deriv = (total_prob * (*second_deriv) - ((*first_deriv) * (*first_deriv)))
      / (total_prob * total_prob); /* deriv of log followed by quotient rule */
  *first_deriv = *first_deriv / total_prob; /* deriv of log */
  total_prob = log(total_prob);

  if (scratch == NULL) {
    for (j = 0; j < nstates; j++) {
      free(L[j]);
      free(LL[j]);
      if (second_deriv != NULL) free(LLL[j]);
    }
    free(L);
    free(LL);
    if (second_deriv != NULL) free(LLL);
  }

  return(total_prob);
}

/* Compute the first and (optionally) second derivatives with respect
   to the scale parameter for the single-column log likelihood
   function (col_compute_log_likelihood).  This version assumes scale
   parameters for the whole tree and for the subtree.  Return value is
   log likelihood, which is computed as a by-product.  Derivs will be
   stored in *first_deriv and *second_deriv.  If second_deriv == NULL,
   it will not be computed (saves some time).  */
double col_scale_derivs_subtree(ColFitData *d, Vector *gradient, 
                                Matrix *hessian, double ***scratch) {
  int i, j, k, nodeidx, rcat;
  int nstates = d->mod->rate_matrix->size;
  TreeNode *n;
  double total_prob = 0;
  List *traversal = tr_postorder(d->mod->tree);  
  double **L;                   /* partial likelihoods */
  double **LL;                  /* 1st deriv of partial likelihoods wrt
                                   1st scale param */ 
  double **LLL;                 /* 2nd deriv of partial likelihoods
                                   wrt 1st scale param */
  double **MM;                  /* 1st deriv of partial likelihoods
                                   wrt 2nd scale param */
  double **MMM;                 /* 2nd deriv of partial likelihoods
                                   wrt 2nd scale param */
  double **NNN;                 /* 2nd cross deriv (off diagonal in
                                   Hessian) of partial likelihoods */

  double *pd = gradient->data;  /* 1st partial derivatives */
  double **pd2 = (hessian == NULL ? NULL : hessian->data);
                                /* 2nd partial derivatives; because of
                                   symmetry, only pd2[0][0],
                                   pd2[1][1], and pd2[1][0] need to be
                                   considered during computation */

  assert(d->msa->ss->tuple_size == 1);
  assert(d->mod->order == 0);
  assert(d->mod->allow_gaps == TRUE);

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
    pd2[0][0] =  (total_prob*pd2[0][0] - pd[0]*pd[0]) / (total_prob*total_prob);
    pd2[1][1] =  (total_prob*pd2[1][1] - pd[1]*pd[1]) / (total_prob*total_prob);
    pd2[1][0] =  (total_prob*pd2[1][0] - pd[1]*pd[0]) / (total_prob*total_prob);
    pd2[0][1] =  pd2[1][0];
  }
  pd[0] = pd[0] / total_prob; /* deriv of log */
  pd[1] = pd[1] / total_prob; 
  total_prob = log(total_prob);

  if (scratch == NULL) {
    for (j = 0; j < nstates; j++) {
      free(L[j]);
      free(LL[j]);
      free(MM[j]);
      if (pd2 != NULL) {
        free(LLL[j]);
        free(MMM[j]);
        free(NNN[j]);
      }
    }
    free(L);
    free(LL);
    free(MM);
    if (pd2 != NULL) {
      free(LLL);
      free(MMM);
      free(NNN);
    }
  }

  return(total_prob);
}

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

/* Perform a likelihood ratio test for each column tuple in an
   alignment, comparing the given null model with an alternative model
   that has a free scaling parameter for all branches.  Assumes a 0th
   order model, leaf-to-sequence mapping already available, prob
   matrices computed, sufficient stats available.  Computes p-values
   based using the chi-sq distribution and stores them in tuple_pvals.
   Will optionally store the individual scale factors in tuple_scales
   and raw log likelihood ratios in tuple_llrs if these variables are
   non-NULL.  Must define mode as CON (for 0 <= scale <= 1), ACC
   (for 1 <= scale), or NNEUT (0 <= scale) */ 
void col_lrts(TreeModel *mod, MSA *msa, mode_type mode, double *tuple_pvals, 
              double *tuple_scales, double *tuple_llrs, FILE *logf) {
  int i;
  ColFitData *d;
  double null_lnl, alt_lnl, delta_lnl;

  /* init ColFitData */
  d = col_init_fit_data(mod, msa, ALL, mode, FALSE);

  /* iterate through column tuples */
  for (i = 0; i < msa->ss->ntuples; i++) {
    mod->scale = 1;
    tm_set_subst_matrices(mod);

    /* compute log likelihoods under null and alt hypotheses */
    null_lnl = col_compute_log_likelihood(mod, msa, i, d->fels_scratch[0]);

    vec_set(d->params, 0, d->init_scale);
    d->tupleidx = i;
    if (opt_bfgs(col_likelihood_wrapper, d->params, d, &alt_lnl, d->lb, 
                 d->ub, logf, col_grad_wrapper, OPT_HIGH_PREC, NULL) != 0) 
      ;                         /* do nothing; nonzero exit typically
                                   occurs when max iterations is
                                   reached; a warning is printed to
                                   the log */

    alt_lnl *= -1;

    delta_lnl = alt_lnl - null_lnl;
    assert(delta_lnl > -0.01);
    if (delta_lnl < 0) delta_lnl = 0;

    /* compute p-vals via chi-sq */
    if (tuple_pvals != NULL) {
      if (mode == NNEUT) 
        tuple_pvals[i] = chisq_cdf(2*delta_lnl, 1, FALSE);
      else
        tuple_pvals[i] = half_chisq_cdf(2*delta_lnl, 1, FALSE);
    /* assumes 50:50 mix of chisq and point mass at zero, due to
       bounding of param */
    }

    /* store scales and log likelihood ratios if necessary */
    if (tuple_scales != NULL) tuple_scales[i] = vec_get(d->params, 0);
    if (tuple_llrs != NULL) tuple_llrs[i] = delta_lnl;
  }
  
  col_free_fit_data(d);
}

void col_lrts_sub(TreeModel *mod, MSA *msa, mode_type mode, 
                  double *tuple_pvals, double *tuple_null_scales, 
                  double *tuple_scales, double *tuple_sub_scales, 
                  double *tuple_llrs, FILE *logf) {
  int i;
  ColFitData *d, *d2;
  double null_lnl, alt_lnl, delta_lnl;
  TreeModel *modcpy;

  modcpy = tm_create_copy(mod);   /* need separate copy of tree model
                                     with different internal scaling
                                     data for supertree/subtree case */

  /* init ColFitData -- one for null model, one for alt */
  d = col_init_fit_data(modcpy, msa, ALL, NNEUT, FALSE);
  d2 = col_init_fit_data(mod, msa, SUBTREE, mode, FALSE); 
                                /* mod has the subtree info, modcpy
                                   does not */

  /* iterate through column tuples */
  for (i = 0; i < msa->ss->ntuples; i++) {
    /* compute log likelihoods under null and alt hypotheses */
    d->tupleidx = i;
    vec_set(d->params, 0, d->init_scale);
    if (opt_bfgs(col_likelihood_wrapper, d->params, d, &null_lnl, d->lb, 
                 d->ub, logf, col_grad_wrapper, OPT_HIGH_PREC, NULL) != 0)
      die("ERROR in estimation of scale for tuple %d.\n", i);
    null_lnl *= -1;

    d2->tupleidx = i;
    vec_set(d2->params, 0, d2->init_scale);
    vec_set(d2->params, 1, d2->init_scale_sub);
    if (opt_bfgs(col_likelihood_wrapper, d2->params, d2, &alt_lnl, d2->lb, 
                 d2->ub, logf, col_grad_wrapper, OPT_HIGH_PREC, NULL) != 0)
      die("ERROR in estimation of supertree/subtree scale for tuple %d.\n", i);
    alt_lnl *= -1;

    delta_lnl = alt_lnl - null_lnl;
    assert(delta_lnl > -0.1);
    if (delta_lnl < 0 || delta_lnl < 0.001) delta_lnl = 0;
    /* within tolerance of optimization */

    /* compute p-vals via chi-sq */
    if (tuple_pvals != NULL) {
      if (mode == NNEUT) 
        tuple_pvals[i] = chisq_cdf(2*delta_lnl, 1, FALSE);
      else
        tuple_pvals[i] = half_chisq_cdf(2*delta_lnl, 1, FALSE);
      /* assumes 50:50 mix of chisq and point mass at zero, due to
         bounding of param */
    }

    /* store scales and log likelihood ratios if necessary */
    if (tuple_null_scales != NULL) 
      tuple_null_scales[i] = vec_get(d->params, 0);
    if (tuple_scales != NULL) 
      tuple_scales[i] = vec_get(d2->params, 0);
    if (tuple_sub_scales != NULL) 
      tuple_sub_scales[i] = vec_get(d2->params, 1);
    if (tuple_llrs != NULL) 
      tuple_llrs[i] = delta_lnl;
  }
  
  col_free_fit_data(d);
  col_free_fit_data(d2);
  modcpy->estimate_branchlens = TM_BRANCHLENS_ALL; 
                                /* have to revert for tm_free to work
                                   correctly */
  tm_free(modcpy);
}

void col_score_tests(TreeModel *mod, MSA *msa, double *tuple_pvals, 
                     double *tuple_derivs, double *tuple_teststats) {
  int i;
  ColFitData *d;
  double first_deriv, second_deriv, teststat;

  /* init ColFitData */
  d = col_init_fit_data(mod, msa, ALL, NNEUT, FALSE);
  /* FIXME: is a one-sided test possible? */

  /* iterate through column tuples */
  for (i = 0; i < msa->ss->ntuples; i++) {
    d->tupleidx = i;
    col_scale_derivs(d, &first_deriv, &second_deriv, d->fels_scratch);

    teststat = -first_deriv*first_deriv / second_deriv;
    if (teststat < 0) teststat = 0;

    if (tuple_pvals != NULL)
      tuple_pvals[i] = chisq_cdf(teststat, 1, FALSE);

    /* store scales and log likelihood ratios if necessary */
    if (tuple_derivs != NULL) tuple_derivs[i] = first_deriv;
    if (tuple_teststats != NULL) tuple_teststats[i] = teststat;
  }

  col_free_fit_data(d);
}

void col_score_tests_sub(TreeModel *mod, MSA *msa, double *tuple_pvals, 
                         double *tuple_null_scales, double *tuple_derivs,
                         double *tuple_sub_derivs, double *tuple_teststats,
                         FILE *logf) {
  int i;
  ColFitData *d, *d2;
  Vector *grad = vec_new(2);
  Matrix *hessian = mat_new(2, 2);
  double det, A, B, C, D, E, F, lnl, teststat;
  TreeModel *modcpy = tm_create_copy(mod); /* need separate copy of tree model
                                              with different internal scaling
                                              data for supertree/subtree case */

  /* init ColFitData -- one for null model, one for alt */
  d = col_init_fit_data(modcpy, msa, ALL, NNEUT, FALSE);
  d2 = col_init_fit_data(mod, msa, SUBTREE, NNEUT, TRUE); 
                                /* mod has the subtree info, modcpy
                                   does not */

  /* FIXME: is a one-sided test possible? */

  /* iterate through column tuples */
  for (i = 0; i < msa->ss->ntuples; i++) {
    d->tupleidx = i;
    vec_set(d->params, 0, d->init_scale);
    if (opt_bfgs(col_likelihood_wrapper, d->params, d, &lnl, d->lb, 
                 d->ub, logf, col_grad_wrapper, OPT_HIGH_PREC, NULL) != 0)
      die("ERROR in estimation of scale for tuple %d.\n", i);

    d2->tupleidx = i;
    d2->mod->scale = d->params->data[0];
    tm_set_subst_matrices(d2->mod);
    col_scale_derivs_subtree(d2, grad, hessian, d2->fels_scratch);

    det = hessian->data[0][0] * hessian->data[1][1]
      - hessian->data[0][1] * hessian->data[1][0];
    assert(det != 0);
    A = hessian->data[1][1] / det; /* cell 0,0 of inverse Fisher matrix */
    B = -hessian->data[0][1] / det; /* cell 0,1 */
    C = -hessian->data[1][0] / det; /* cell 1,0 */
    D = hessian->data[0][0] / det;  /* cell 1,1 */
    E = grad->data[0];
    F = grad->data[1];

    teststat = E*E*A + E*F*C + E*B*F + F*F*D;
    /* grad' * inv_fish * grad */

    if (teststat < 0) teststat = 0;

    if (tuple_pvals != NULL)
      tuple_pvals[i] = chisq_cdf(teststat, 1, FALSE);

    /* store scales and log likelihood ratios if necessary */
    if (tuple_null_scales != NULL) tuple_null_scales[i] = d->params->data[0];
    if (tuple_derivs != NULL) tuple_derivs[i] = E;
    if (tuple_sub_derivs != NULL) tuple_sub_derivs[i] = F;
    if (tuple_teststats != NULL) tuple_teststats[i] = teststat;
  }

  col_free_fit_data(d);
  col_free_fit_data(d2);
  vec_free(grad);
  mat_free(hessian);
  modcpy->estimate_branchlens = TM_BRANCHLENS_ALL; 
                                /* have to revert for tm_free to work
                                   correctly */
  tm_free(modcpy);
}

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
  if (stype == SUBTREE) 
    assert(mod->subtree_root != NULL);
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

  d->expdiag = zvec_new(size); 
  d->nfels_scratch = stype == SUBTREE ? 6 : 3;
  d->fels_scratch = smalloc(d->nfels_scratch * sizeof(void*));
  for (i = 0; i < d->nfels_scratch; i++) {
    d->fels_scratch[i] = smalloc(size * sizeof(void*));
    for (j = 0; j < size; j++) 
      d->fels_scratch[i][j] = smalloc((nnodes+1) * sizeof(double)); 
  }
  d->mat_scratch = zmat_new(size, size);   
  d->vec_scratch1 = zvec_new(size);
  d->vec_scratch2 = zvec_new(size);
  return d;
}

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
    free(d->PP[nid]);
    free(d->PPP[nid]);
    free(d->QQ[nid]);
    free(d->QQQ[nid]);
    free(d->RRR[nid]);
  }
  free(d->PP);
  free(d->PPP);
  free(d->QQ);
  free(d->QQQ);
  free(d->RRR);

  zvec_free(d->expdiag); 
  for (i = 0; i < d->nfels_scratch; i++) {
    for (j = 0; j < d->mod->rate_matrix->size; j++) 
      free(d->fels_scratch[i][j]);
    free(d->fels_scratch[i]);
  }
  free(d->fels_scratch);
  zmat_free(d->mat_scratch); 
  zvec_free(d->vec_scratch1);
  zvec_free(d->vec_scratch2);

  free(d);
}
