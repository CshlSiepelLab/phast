/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: tree_likelihoods.c,v 1.14 2008-11-12 02:07:59 acs Exp $ */

#include <tree_likelihoods.h>
#include <subst_mods.h>
#include <markov_matrix.h>
#include <subst_mods.h>
#include <dgamma.h>
#include <sufficient_stats.h>

/* Computation of likelihoods for columns of a given multiple
   alignment, according to a given tree model.  */

/* NOTE: routines currently assume that the names at the leaves of the
   tree (as defined in the TreeModel object) are numbers from 1 to
   nseqs.  */

/* FIXME: inside and outside computations should be in log space --
   probably only an issue when the number of leaves is large */

int tuple_index_missing_data(char *tuple, int *inv_alph, int *is_missing,
                             int alph_size);



/* Compute the likelihood of a tree model with respect to an
   alignment.  Optionally retain column-by-column likelihoods,
   optionally compute posterior probabilities.  If 'post' is NULL, no
   posterior probabilities (or related quantities) will be computed.
   If 'post' is non-NULL each of its attributes must either be NULL or
   previously allocated to the required size. */
double tl_compute_log_likelihood(TreeModel *mod, MSA *msa, 
                                 double *col_scores, int cat,
                                 TreePosteriors *post) {

  int i, j;
  double retval = 0;
  int nstates = mod->rate_matrix->size;
  int alph_size = strlen(mod->rate_matrix->states); 
  int npasses = (mod->order > 0 && mod->use_conditionals == 1 ? 2 : 1); 
  int pass, col_offset, k, nodeidx, rcat, /* colidx, */ tupleidx, defined;
  TreeNode *n;
  double total_prob, marg_tot;
  List *traversal;
  double **inside_joint = NULL, **inside_marginal = NULL, 
    **outside_joint = NULL, **outside_marginal = NULL, 
    ****subst_probs = NULL;
  double *tuple_scores=NULL;
  double rcat_prob[mod->nratecats];
  double tmp[nstates];

  checkInterrupt();

  /* allocate memory */
  inside_joint = (double**)smalloc(nstates * sizeof(double*));
  for (j = 0; j < nstates; j++) 
    inside_joint[j] = (double*)smalloc((mod->tree->nnodes+1) * 
                                       sizeof(double)); 
  outside_joint = (double**)smalloc(nstates * sizeof(double*));
  for (j = 0; j < nstates; j++) 
    outside_joint[j] = (double*)smalloc((mod->tree->nnodes+1) * 
                                        sizeof(double)); 
  /* only needed if post != NULL? */
  if (mod->order > 0) {
    inside_marginal = (double**)smalloc(nstates * sizeof(double*));
    for (j = 0; j < nstates; j++) 
      inside_marginal[j] = (double*)smalloc((mod->tree->nnodes+1) * 
                                            sizeof(double));
  }
  if (mod->order > 0 && post != NULL) {
    outside_marginal = (double**)smalloc(nstates * sizeof(double*));
    for (j = 0; j < nstates; j++) 
      outside_marginal[j] = (double*)smalloc((mod->tree->nnodes+1) * 
                                             sizeof(double));
  }
  if (post != NULL) {
    subst_probs = (double****)smalloc(mod->nratecats * sizeof(double***));
    for (rcat = 0; rcat < mod->nratecats; rcat++) {
      subst_probs[rcat] = (double***)smalloc(nstates * sizeof(double**));
      for (j = 0; j < nstates; j++) {
        subst_probs[rcat][j] = (double**)smalloc(nstates * sizeof(double*));
        for (k = 0; k < nstates; k++)
          subst_probs[rcat][j][k] = (double*)smalloc(mod->tree->nnodes * sizeof(double));          
      }
    }
  }

  /* create IUPAC mapping if needed */
  if (mod->iupac_inv_map == NULL)
    mod->iupac_inv_map = build_iupac_inv_map(mod->rate_matrix->inv_states, 
                                             alph_size);

  
  if (cat > msa->ncats)
    die("ERROR tl_compute_log_likelihood: cat (%i) > msa->ncats (%i)\n", cat, msa->ncats);

  if (!(cat < 0 || col_scores == NULL || msa->categories != NULL))
    die("ERROR tl_compute_log_likelihood: cat=%i, col_scores==NULL=%i, msa->categories==NULL=%i\n", cat, col_scores==NULL, msa->categories==NULL);
  /* if using categories and col-by-col
     scoring, then must have col-by-col
     categories */

  /* obtain sufficient statistics, if necessary */
  if (msa->ss != NULL){ 
    if (msa->ss->tuple_size <= mod->order)
      die("ERROR tl_compute_log_likelihood: tuple_size (%i) must be greater than mod->order (%i)\n",
	  msa->ss->tuple_size, mod->order);
  }
  else 
    ss_from_msas(msa, mod->order+1, col_scores == NULL ? 0 : 1, 
                 NULL, NULL, NULL, -1, subst_mod_is_codon_model(mod->subst_mod));

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
/*   printf("mod %d\n", mod); */
  if (!defined) {
    tm_set_subst_matrices(mod);
  }
  if (col_scores != NULL) {
    tuple_scores = (double*)smalloc(msa->ss->ntuples * sizeof(double));
    for (tupleidx = 0; tupleidx < msa->ss->ntuples; tupleidx++)
      tuple_scores[tupleidx] = 0;
  }

  if (post != NULL && post->expected_nsubst_tot != NULL) {
    for (rcat = 0; rcat < mod->nratecats; rcat++)
      for (i = 0; i < nstates; i++) 
        for (j = 0; j < nstates; j++) 
          for (k = 0; k < mod->tree->nnodes; k++)
            post->expected_nsubst_tot[rcat][i][j][k] = 0;
  }
  if (post != NULL && post->rcat_expected_nsites != NULL)
    for (rcat = 0; rcat < mod->nratecats; rcat++)
      post->rcat_expected_nsites[rcat] = 0;

  for (tupleidx = 0; tupleidx < msa->ss->ntuples; tupleidx++) {
    int skip_fels = FALSE;

    if ((cat >= 0 && msa->ss->cat_counts[cat][tupleidx] == 0) || 
        (cat < 0 && msa->ss->counts[tupleidx] == 0))
      continue;
    checkInterruptN(tupleidx, 1000);

    total_prob = 0;
    marg_tot = NULL_LOG_LIKELIHOOD;

    /* check for gaps and whether column is informative, if necessary */
    if (!mod->allow_gaps)
      for (j = 0; !skip_fels && j < msa->nseqs; j++) 
        if (ss_get_char_tuple(msa, tupleidx, j, 0) == GAP_CHAR) 
          skip_fels = TRUE;
    if (!skip_fels && mod->inform_reqd) {
      int ninform = 0;
      for (j = 0; j < msa->nseqs; j++) {
        if (msa->is_informative != NULL && !msa->is_informative[j])
          continue;
        else if (!msa->is_missing[(int)ss_get_char_tuple(msa, tupleidx, j, 0)])
          ninform++;
      }
      if (ninform < 2) skip_fels = TRUE;
    }
          
    if (!skip_fels) {
      for (pass = 0; pass < npasses; pass++) {
        double **pL = (pass == 0 ? inside_joint : inside_marginal);
        double **pLbar = (pass == 0 ? outside_joint : outside_marginal);
        /*         TreePosteriors *postpass = (pass == 0 ? post : postmarg); */

        if (pass > 0)
          marg_tot = 0;         /* will need to compute */

        for (rcat = 0; rcat < mod->nratecats; rcat++) {
          traversal = tr_postorder(mod->tree);      
          for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
            int partial_match[mod->order+1][alph_size];
            n = lst_get_ptr(traversal, nodeidx);      
            if (n->lchild == NULL) { 
              /* leaf: base case of recursion */
              int thisseq;

	      if (n->name == NULL)
		die("ERROR tl_compute_log_likelihood: n->name is NULL\n");
              thisseq = mod->msa_seq_idx[n->id];

              /* first figure out whether there is a match for each
                 character in each position; we'll call this the record of
                 "partial_matches". */
              for (col_offset = -1*mod->order; col_offset <= 0; col_offset++) {
                int observed_state = -1;
                int *iupac_prob = NULL;

                if (pass == 0 || col_offset < 0) {
                  char thischar = (int)ss_get_char_tuple(msa, tupleidx, 
                                                         thisseq, col_offset);
                  observed_state = mod->rate_matrix->inv_states[(int)thischar];
                  if (observed_state < 0)
                    iupac_prob = mod->iupac_inv_map[(int)thischar];
                }
 
                /* otherwise, we're on a second pass and looking the
                   current base, so we want to use the "missing
                   information" principle */

                if (iupac_prob != NULL) {
                  for (i = 0; i < alph_size; i++) 
                    partial_match[mod->order+col_offset][i] = iupac_prob[i];
                }
                else {
                  for (i = 0; i < alph_size; i++) {
                    if (observed_state < 0 || i == observed_state) 
                      partial_match[mod->order+col_offset][i] = 1;
                    else
                      partial_match[mod->order+col_offset][i] = 0; 
                  }
                }
              }

              /* now find the intersection of the partial matches */
              for (i = 0; i < nstates; i++) {
                if (mod->order == 0)  /* handle 0th order model as special
                                         case, for efficiency.  In this case
                                         the partial match *is* the total
                                         match */
                  pL[i][n->id] = partial_match[0][i];
                else {
                  int total_match = 1;
                  /* figure out the "projection" of state i in the dimension
                     of each position, and see whether there is a
                     corresponding partial match. */
                  /* NOTE: mod->order is approx equal to log nstates
                     (prob no more than 2) */
                  for (col_offset = -1*mod->order; col_offset <= 0 && total_match; 
                       col_offset++) {
                    int projection = (i / int_pow(alph_size, -1 * col_offset)) % 
                      alph_size;

                    if (!partial_match[mod->order+col_offset][projection]) 
                      total_match = 0; /* must have partial matches in all
                                          dimensions for a total match */
                  }
                  pL[i][n->id] = total_match;
                }
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

          if (post != NULL && pass == 0) {
            MarkovMatrix *subst_mat;
            double this_total, denom;

            /* do outside calculation */
            traversal = tr_preorder(mod->tree);
            for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
              n = lst_get_ptr(traversal, nodeidx);
              if (n->parent == NULL) { /* base case */ 
                for (i = 0; i < nstates; i++) 
                  pLbar[i][n->id] = vec_get(mod->backgd_freqs, i);
              }
              else {            /* recursive case */
                TreeNode *sibling = (n == n->parent->lchild ? 
                                     n->parent->rchild : n->parent->lchild);
                MarkovMatrix *par_subst_mat = mod->P[n->id][rcat];
                MarkovMatrix *sib_subst_mat = mod->P[sibling->id][rcat];

                /* breaking this computation into two parts as follows
                   reduces its complexity by a factor of nstates */

                for (j = 0; j < nstates; j++) { /* parent state */
                  tmp[j] = 0;
                  for (k = 0; k < nstates; k++) { /* sibling state */
                    tmp[j] += pLbar[j][n->parent->id] *
                      pL[k][sibling->id] * mm_get(sib_subst_mat, j, k);
                  }
                }

                for (i = 0; i < nstates; i++) { /* child state */
                  pLbar[i][n->id] = 0;
                  for (j = 0; j < nstates; j++) { /* parent state */
                    pLbar[i][n->id] += 
                      tmp[j] * mm_get(par_subst_mat, j, i);
                  }
                }
              }
              

              /* compute total probability based on current node, to
                 avoid numerical errors */
              this_total = 0;
              for (i = 0; i < nstates; i++)
                this_total += pL[i][n->id] * pLbar[i][n->id];

              if (post->expected_nsubst != NULL && n->parent != NULL)
                post->expected_nsubst[rcat][n->id][tupleidx] = 1;

              subst_mat = mod->P[n->id][rcat];
              for (i = 0; i < nstates; i++) {
                /* compute posterior prob of base (tuple) i at node n */
                if (post->base_probs != NULL) {
                  post->base_probs[rcat][i][n->id][tupleidx] = 
                    safediv(pL[i][n->id] * pLbar[i][n->id], this_total);
                }

                if (n->parent == NULL) continue;

                /* (intermediate computation used for subst probs) */
                denom = 0;
                for (k = 0; k < nstates; k++)
                  denom += pL[k][n->id] * mm_get(subst_mat, i, k);

                for (j = 0; j < nstates; j++) { 
                  /* compute posterior prob of a subst of base j at
                     node n for base i at node n->parent */
                  subst_probs[rcat][i][j][n->id] = 
                    safediv(pL[i][n->parent->id] * pLbar[i][n->parent->id], 
                            this_total) *
                    pL[j][n->id] * mm_get(subst_mat, i, j);
                  subst_probs[rcat][i][j][n->id] = 
                    safediv(subst_probs[rcat][i][j][n->id], denom);
                    
                  if (post->subst_probs != NULL)
                    post->subst_probs[rcat][i][j][n->id][tupleidx] = 
                      subst_probs[rcat][i][j][n->id];
                  
                  if (post->expected_nsubst != NULL && j == i)
                    post->expected_nsubst[rcat][n->id][tupleidx] -= 
                      subst_probs[rcat][i][j][n->id];

                }
              }
            }
          }

          if (pass == 0) {
            rcat_prob[rcat] = 0;
            for (i = 0; i < nstates; i++) {
              rcat_prob[rcat] += vec_get(mod->backgd_freqs, i) * 
                inside_joint[i][mod->tree->id] * mod->freqK[rcat];
            }
            total_prob += rcat_prob[rcat];
          }
          else { 
            for (i = 0; i < nstates; i++) 
              marg_tot += vec_get(mod->backgd_freqs, i) * 
                inside_marginal[i][mod->tree->id] * mod->freqK[rcat]; 
          }
        } /* for rcat */
      } /* for pass */
    } /* if skip_fels */

      /* compute posterior prob of each rate cat and related quantities */
    if (post != NULL) {
      if (skip_fels) die("ERROR: tl_compute_log_likelihood: skip_fels should be 0 but is %i\n", skip_fels);
      for (rcat = 0; rcat < mod->nratecats; rcat++) {
        double rcat_post_prob = safediv(rcat_prob[rcat], total_prob);
        if (post->rcat_probs != NULL) 
          post->rcat_probs[rcat][tupleidx] = rcat_post_prob;
        if (post->rcat_expected_nsites != NULL)
          post->rcat_expected_nsites[rcat] += rcat_post_prob * 
            (cat >= 0 ? msa->ss->cat_counts[cat][tupleidx] : 
             msa->ss->counts[tupleidx]);
        if (post->expected_nsubst_tot != NULL) {
          for (nodeidx = 0; nodeidx < mod->tree->nnodes; nodeidx++) {
            n = lst_get_ptr(mod->tree->nodes, nodeidx);
            if (n->parent == NULL) continue;
            for (i = 0; i < nstates; i++) 
              for (j = 0; j < nstates; j++) 
                post->expected_nsubst_tot[rcat][i][j][n->id] += 
                  subst_probs[rcat][i][j][n->id] * 
                  (cat >= 0 ? msa->ss->cat_counts[cat][tupleidx] : 
                   msa->ss->counts[tupleidx]) *
                  rcat_post_prob;
          }
        }
	if (post->expected_nsubst_col != NULL) {
	  for (nodeidx = 0; nodeidx < mod->tree->nnodes; nodeidx++) {
            n = lst_get_ptr(mod->tree->nodes, nodeidx);
            if (n->parent == NULL) continue;
            for (i = 0; i < nstates; i++) 
              for (j = 0; j < nstates; j++) 
                post->expected_nsubst_col[rcat][n->id][tupleidx][i][j] = 
                  subst_probs[rcat][i][j][n->id] * rcat_post_prob;
          }
        }
      }
    }

    if (mod->order > 0 && mod->use_conditionals == 1 && !skip_fels) 
      total_prob /= marg_tot; 

    total_prob = log2(total_prob);

    if (col_scores != NULL && 
        (cat < 0 || msa->ss->cat_counts[cat][tupleidx] > 0))
      tuple_scores[tupleidx] = total_prob;
    /* NOTE: tuple_scores contains the
     *unweighted* (log) probabilities. */

    total_prob *= (cat >= 0 ? msa->ss->cat_counts[cat][tupleidx] : 
                   msa->ss->counts[tupleidx]); /* log space */

    retval += total_prob;     /* log space */
        
  } /* for tupleidx */

  for (j = 0; j < nstates; j++) {
    sfree(inside_joint[j]);
    sfree(outside_joint[j]);
    if (mod->order > 0) sfree(inside_marginal[j]);
    if (mod->order > 0 && post != NULL) sfree(outside_marginal[j]);
  }
  sfree(inside_joint);
  sfree(outside_joint);
  if (mod->order > 0) sfree(inside_marginal);
  if (mod->order > 0 && post != NULL) sfree(outside_marginal);
  if (col_scores != NULL) {
    if (cat >= 0) 
      for (i = 0; i < msa->length; i++)
        col_scores[i] = msa->categories[i] == cat ?
          tuple_scores[msa->ss->tuple_idx[i]] :
          NEGINFTY;
    else
      for (i = 0; i < msa->length; i++)
        col_scores[i] = tuple_scores[msa->ss->tuple_idx[i]];
    sfree(tuple_scores);
  }
  if (post != NULL) {
    for (rcat = 0; rcat < mod->nratecats; rcat++) {
      for (j = 0; j < nstates; j++) {
        for (k = 0; k < nstates; k++)
          sfree(subst_probs[rcat][j][k]);
        sfree(subst_probs[rcat][j]);
      }
      sfree (subst_probs[rcat]);
    }
    sfree(subst_probs);
  }
  return(retval);
}

/* this is retained for possible use in the future; not using weight
   matrices for much anymore */
void tl_compute_log_likelihood_weight_matrix(TreeModel *mod, MSA *msa, 
                                             double *col_scores, int cat) {
  int i, seq, idx, alph_size = strlen(msa->alphabet);
  double retval = 0;
  char tuple[mod->order + 2];
  Vector *margfreqs = 
    get_marginal_eq_freqs(mod->rate_matrix->states, mod->order+1,
                          mod->backgd_freqs);
  int col_by_col = (col_scores != NULL || msa->ss == NULL);
                                /* evaluate the alignment
                                   column-by-column if either
                                   col-by-col scores are required or
                                   the sufficient stats are not
                                   available */
                                /* NOTE: !col_by_col -> msa->ss != NULL */

  checkInterrupt();
  tuple[mod->order+1] = '\0';

  if (mod->tree != NULL)
    die("ERROR tl_compute_log_likelihood_weight_matrix: mod->tree should be NULL\n");
  if (msa->ss == NULL && msa->seqs == NULL)
    die("ERROR tl_compute_log_likelihood_weight_matrix: mod->ss and mod->seqs are both NULL\n");

  if (cat >= 0) {
    if (col_by_col) {
      if (msa->categories == NULL)
	die("ERROR tl_compute_log_likelihood_weight_matrix: msa->categories is NULL\n");
    }
    /* if using categories and col-by-col
       scoring, then must have col-by-col
       categories */
    else if (msa->ss->cat_counts == NULL)
      die("ERROR tl_compute_log_likelihood_weight_matrix: msa->ss->cat_counts is NULL\n");
    /* if using categories and unordered
       sufficient statistics, must have
       category-by-category counts */
  }

  if (col_by_col)
    if (msa->seqs == NULL && msa->ss->tuple_idx == NULL)
      die("ERROR tl_compute_log_likelihood requires ordered alignment\n");
  /* if using col-by-col scoring, must
     have ordered representation */
    
  retval = 0;
  for (idx = 0; 
       idx < (col_by_col ? msa->length : msa->ss->ntuples);
       idx++) {
    int thisstate, col, tupleidx;
    double col_val = 0, prob = 0;

    /* NOTE: when evaluating col-by-col, idx is a column, but otherwise
       idx is a tuple index.  Let's be clear about this ... */
    col = (col_by_col ? idx : -1); 
    if (msa->ss == NULL) tupleidx = -1; 
    else if (col_by_col) tupleidx = msa->ss->tuple_idx[col];
    else tupleidx = idx;      /* NOTE: tupleidx will be defined
                                 whenever msa->ss != NULL */

    if (cat < 0 || (col_by_col && msa->categories[col] == cat) ||
        (!col_by_col && msa->ss->cat_counts[cat][tupleidx] > 0)) {

      for (seq = 0; seq < msa->nseqs; seq++) {

        for (i = -mod->order; i <= 0; i++) {
          if (msa->ss != NULL)
            tuple[mod->order+i] = ss_get_char_tuple(msa, tupleidx, seq, i); 
          else if (col + i >= 0)
            tuple[mod->order+i] = msa->seqs[seq][col+i];
          else 
            tuple[mod->order+i] = msa->missing[0];
        }

        if (!mod->allow_gaps && 
            msa->inv_alphabet[(int)tuple[mod->order]] < 0 &&
            !msa->is_missing[(int)tuple[mod->order]]) {
            
          col_val = NEGINFTY; /* we want to apply the strict penalty
                                 iff there is an unrecognized
                                 character in *this* (the rightmost)
                                 col; missing data is a special case -- don't
                                 penalize even if not in alphabet */

                                /* FIXME: seems too complicated --
                                   just check for gap? */
          break; 
        }
        else if (mod->allow_but_penalize_gaps &&
                 msa->inv_alphabet[(int)tuple[mod->order]] < 0) { 
          /* temporary */
          double tmp_prob;
          prob = 1;
          for (i = 0; i < alph_size; i++) {
            tuple[mod->order] = msa->alphabet[i];
            tmp_prob = vec_get(margfreqs, tuple_index_missing_data(tuple, msa->inv_alphabet, msa->is_missing, alph_size));
            if (tmp_prob < prob) prob = tmp_prob;
          }
          if (prob == 0) prob = 0.01;
        }
        else {
          thisstate = tuple_index_missing_data(tuple, msa->inv_alphabet, 
                                               msa->is_missing, alph_size);
          prob = vec_get(margfreqs, thisstate);
        }

        if (prob == 0) { col_val = NEGINFTY; break; }

        col_val += log2(prob);

        if (mod->use_conditionals && mod->order > 0) {
          tuple[mod->order] = msa->missing[0];
          col_val -= log2(vec_get(margfreqs, tuple_index_missing_data(tuple, msa->inv_alphabet, msa->is_missing, alph_size)));
        }
      }
    }
    if (!col_by_col)   /* tuple-by-tuple scoring */
      col_val *= (cat >= 0 ? msa->ss->cat_counts[cat][tupleidx] : 
                  msa->ss->counts[tupleidx]);
    retval += col_val;
    if (col_scores != NULL) col_scores[col] = col_val;
  }
  if (retval < NEGINFTY) retval = NEGINFTY; 
  /* must be true if any of the columns
     considered had prob NEGINFTY */
  vec_free(margfreqs);
}


TreePosteriors *tl_new_tree_posteriors(TreeModel *mod, MSA *msa, int do_bases, 
                                       int do_substs, int do_expected_nsubst, 
                                       int do_expected_nsubst_tot,
				       int do_expected_nsubst_col,
                                       int do_rate_cats, int do_rate_cats_exp) {
  int i, j, k, r, ntuples, nnodes, nstates;
  TreePosteriors *tp = (TreePosteriors*)smalloc(sizeof(TreePosteriors));

  if (mod->tree ==  NULL)
    die("ERROR tl_new_tree_posteriors: mod->tree is NULL\n");
  if (msa->ss == NULL)
    die("ERROR tl_new_tree_posteriors: msa->ss is NULL\n");

  ntuples = msa->ss->ntuples;
  nnodes = mod->tree->nnodes;
  nstates = mod->rate_matrix->size;

  if (do_bases) {
    tp->base_probs = (double****)smalloc(mod->nratecats * sizeof(double***));
    for (r = 0; r < mod->nratecats; r++) {
      tp->base_probs[r] = (double***)smalloc(nstates * sizeof(double**));
      for (i = 0; i < nstates; i++) {
        tp->base_probs[r][i] = (double**)smalloc(nnodes * sizeof(double*));
        for (j = 0; j < nnodes; j++) {
          tp->base_probs[r][i][j] = (double*)smalloc(ntuples * sizeof(double));
        }
      }
    }
  }
  else tp->base_probs = NULL;

  if (do_substs) {
    tp->subst_probs = (double*****)smalloc(mod->nratecats * sizeof(double****));
    for (r = 0; r < mod->nratecats; r++) {
      checkInterrupt();
      tp->subst_probs[r] = (double****)smalloc(nstates * sizeof(double***));
      for (i = 0; i < nstates; i++) {
        tp->subst_probs[r][i] = (double***)smalloc(nstates * sizeof(double**));
        for (j = 0; j < nstates; j++) {
          tp->subst_probs[r][i][j] = (double**)smalloc(nnodes * sizeof(double*));
          for (k = 0; k < nnodes; k++) 
            if (k != mod->tree->id) /* don't need one for the root */
              tp->subst_probs[r][i][j][k] = (double*)smalloc(ntuples * 
                                                            sizeof(double));
            else 
              tp->subst_probs[r][i][j][k] = NULL;
        }
      }
    }
  }
  else tp->subst_probs = NULL;

  if (do_expected_nsubst) {
    tp->expected_nsubst = (double***)smalloc(mod->nratecats * sizeof(double**));
    for (r = 0; r < mod->nratecats; r++) {
      tp->expected_nsubst[r] = (double**)smalloc(nnodes * sizeof(double*));
      for (i = 0; i < nnodes; i++) {
        if (i != mod->tree->id)
          tp->expected_nsubst[r][i] = (double*)smalloc(ntuples * sizeof(double));
        else
          tp->expected_nsubst[r][i] = NULL;
      }
    }
  }
  else tp->expected_nsubst = NULL;

  if (do_expected_nsubst_tot) {
    tp->expected_nsubst_tot = (double****)smalloc(mod->nratecats * sizeof(double***));
    for (r = 0; r < mod->nratecats; r++) {
      tp->expected_nsubst_tot[r] = (double***)smalloc(nstates * sizeof(double**));
      for (i = 0; i < nstates; i++) {
        tp->expected_nsubst_tot[r][i] = (double**)smalloc(nstates * 
                                                         sizeof(double*));
        for (j = 0; j < nstates; j++) 
          tp->expected_nsubst_tot[r][i][j] = (double*)smalloc(nnodes * 
                                                             sizeof(double));
      }
    }
  }
  else tp->expected_nsubst_tot = NULL;

  if (do_expected_nsubst_col) {
    tp->expected_nsubst_col = (double*****)smalloc(mod->nratecats * sizeof(double****));
    for (r=0; r < mod->nratecats; r++) {
      tp->expected_nsubst_col[r] = (double****)smalloc(nnodes * sizeof(double***));
      for (i=0; i < nnodes; i++) {
	tp->expected_nsubst_col[r][i] = (double***)smalloc(ntuples * sizeof(double**));
	for (j=0; j < ntuples; j++) {
	  tp->expected_nsubst_col[r][i][j] = (double**)smalloc(nstates * sizeof(double*));
	  for (k=0; k < nstates; k++) 
	    tp->expected_nsubst_col[r][i][j][k] = (double*)smalloc(nstates * sizeof(double));
	}
      }
    }
  }
  else tp->expected_nsubst_col = NULL;

  if (do_rate_cats) {
    tp->rcat_probs = (double**)smalloc(mod->nratecats * sizeof(double*));
    for (i = 0; i < mod->nratecats; i++)
      tp->rcat_probs[i] = (double*)smalloc(ntuples * sizeof(double));
  }
  else tp->rcat_probs = NULL;

  if (do_rate_cats_exp) 
    tp->rcat_expected_nsites = (double*)smalloc(mod->nratecats * sizeof(double));
  else tp->rcat_expected_nsites = NULL;

  return tp;
}

void tl_free_tree_posteriors(TreeModel *mod, MSA *msa, TreePosteriors *tp) {
  int i, j, k, r, ntuples, nnodes, nstates;

  if (mod->tree == NULL)
    die("ERROR tl_free_tree_posteriors: mod->tree is NULL\n");
  if (msa->ss == NULL)
    die("ERROR tl_free_tree_posteriors: msa->ss is NULL\n");
  ntuples = msa->ss->ntuples;
  nnodes = mod->tree->nnodes;
  nstates = mod->rate_matrix->size;

  if (tp->base_probs != NULL) {
    for (r = 0; r < mod->nratecats; r++) {
      for (i = 0; i < nstates; i++) {
        for (j = 0; j < nnodes; j++) 
          if (tp->base_probs[r][i][j] != NULL)
            sfree(tp->base_probs[r][i][j]);
        sfree(tp->base_probs[r][i]);
      }
      sfree(tp->base_probs[r]);
    }
    sfree(tp->base_probs);
  }
  if (tp->subst_probs != NULL) {
    for (r = 0; r < mod->nratecats; r++) {
      for (i = 0; i < nstates; i++) {
        for (j = 0; j < nstates; j++) {
          for (k = 0; k < nnodes; k++) 
            if (k != mod->tree->id) 
              sfree(tp->subst_probs[r][i][j][k]);
          sfree(tp->subst_probs[r][i][j]);
        }
        sfree(tp->subst_probs[r][i]);
      }
      sfree(tp->subst_probs[r]);
    }
    sfree(tp->subst_probs);
  }
  if (tp->expected_nsubst != NULL) {
    for (r = 0; r < mod->nratecats; r++) {
      for (i = 0; i < nnodes; i++) 
        if (i != mod->tree->id)
          sfree(tp->expected_nsubst[r][i]);
      sfree(tp->expected_nsubst[r]);
    }
    sfree(tp->expected_nsubst);
  }
  if (tp->expected_nsubst_tot != NULL) {
    for (r = 0; r < mod->nratecats; r++) {
      for (i = 0; i < nstates; i++) {
        for (j = 0; j < nstates; j++) 
          sfree(tp->expected_nsubst_tot[r][i][j]);
        sfree(tp->expected_nsubst_tot[r][i]);
      }
      sfree(tp->expected_nsubst_tot[r]);
    }
    sfree(tp->expected_nsubst_tot);
  }
  if (tp->expected_nsubst_col != NULL) {
    for (r = 0; r < mod->nratecats; r++) {
      for (i=0; i < nnodes; i++) {
	for (j=0; j < ntuples; j++) {
	  for (k=0; k < nstates; k++) 
	    sfree(tp->expected_nsubst_col[r][i][j][k]);
	  sfree(tp->expected_nsubst_col[r][i][j]);
	}
	sfree(tp->expected_nsubst_col[r][i]);
      }
      sfree(tp->expected_nsubst_col[r]);
    }
    sfree(tp->expected_nsubst_col);
  }
  if (tp->rcat_probs != NULL) {
    for (i = 0; i < mod->nratecats; i++)
      sfree(tp->rcat_probs[i]);
    sfree(tp->rcat_probs);           
  }
  if (tp->rcat_expected_nsites != NULL) {
    sfree(tp->rcat_expected_nsites);           
  }

  sfree(tp);
}

/* compute the expected (posterior) complete log likelihood of a tree
   model based on a TreePosteriors object.  Equilibrium frequencies
   are not considered. */
double tl_compute_partial_ll_suff_stats(TreeModel *mod, TreePosteriors *post) {
  double retval = 0;
  int i, j, k, rcat;
  TreeNode *n;
  int nstates = mod->rate_matrix->size;

  for (rcat = 0; rcat < mod->nratecats; rcat++) {
    for (i = 0; i < mod->tree->nnodes; i++) {
      MarkovMatrix *subst_mat;
      if (i == mod->tree->id) continue; /* skip root */
      n = lst_get_ptr(mod->tree->nodes, i);
      subst_mat = mod->P[n->id][rcat];
      for (j = 0; j < nstates; j++) { /* from tuple */
        for (k = 0; k < nstates; k++) { /* to tuple */
          retval += (post->expected_nsubst_tot[rcat][j][k][i] *
                     log2(mm_get(subst_mat, j, k)));
        }
      }
    }
  }
  return retval;
}

/* The functions below are used for computing likelihoods with
   weight-matrix models.  They should possibly be moved to misc.c.  */

/* given an alphabet, a tuple size, and a vector of equilibrium
   frequences, create a new vector of marginal equilibrium
   frequencies describing the space of "meta-tuples", which contain
   actual characters *or* missing data characters.  Each meta-tuple is
   given an equilibrium frequency equal to the sum of the frequencies
   of all "matching" ordinary tuples.  Missing data characters are
   assumed to be gap characters or Ns. */
Vector *get_marginal_eq_freqs (char *alphabet, int tuple_size, 
			       Vector *eq_freqs) {
  int alph_size = strlen(alphabet);
  int ntuples = int_pow(alph_size, tuple_size);
  int i;
  Vector *retval = vec_new(int_pow(alph_size+1, tuple_size));
  vec_zero(retval);

  /* loop through the ordinary (non-meta) tuples */
  for (i = 0; i < ntuples; i++) {
    int digits[tuple_size];
    int j, k, remainder;
    
    /* first decompose the tuple into its "digits" */
    remainder = i;
    for (j = 0; j < tuple_size; j++) { /* from least sig. to most */
      digits[j] = remainder % alph_size;
      remainder /= alph_size;
    }

    /* now consider every pattern of missing-data characters that can
       be overlaid on it.  The equilibrium frequency of the tuple
       contributes to the marginal frequency corresponding to every
       such pattern.  There are 2^tuple_size of them to consider. */
    for (k = 0; k < (1 << tuple_size); k++) {
      int newtuple = 0, base = 1;
      for (j = 0; j < tuple_size; j++) {
        if (k & (1 << j)) 
          newtuple += alph_size * base;
        else 
          newtuple += digits[j] * base;
        base *= (alph_size + 1);
      }
      vec_set(retval, newtuple, 
                     vec_get(retval, newtuple) + 
                     vec_get(eq_freqs, i));
    }
  }
  return retval;
}

/* given a tuple consisting of actual characters and/or missing data,
   return the corresponding state number in the "meta-tuple" space.
   Returns -1 for unallowed tuples */
int tuple_index_missing_data(char *tuple, int *inv_alph, int *is_missing,
                             int alph_size) {
  int retval = 0, i;
  int tuple_size = strlen(tuple);
  for (i = 0; i < tuple_size; i++) {
    int charidx = inv_alph[(int)tuple[tuple_size-i-1]];
    if (charidx < 0) {
      if (tuple[tuple_size-i-1] == GAP_CHAR || 
          is_missing[(int)tuple[tuple_size-i-1]])
        charidx = alph_size;
      else return -1;
    }
    retval += charidx * int_pow(alph_size+1, i);
                                /* i == 0 => least sig. dig; 
                                   i == tuple_size-1 => most sig. dig */
  }
  return retval;
}

