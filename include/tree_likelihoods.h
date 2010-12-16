/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: tree_likelihoods.h,v 1.5 2008-11-12 02:07:59 acs Exp $ */

#ifndef TREE_LIK_H
#define TREE_LIK_H

#include <tree_model.h>
#include <msa.h>
#include <math.h>
#include <misc.h>

/* Structure for information related to posterior probability of tree
   model wrt an alignment.  Each array is indexed as appropriate for
   rate categories, bases in the alignment (indexed according to a
   model's inv_states; with higher order models, actually tuples of
   bases), nodes or edges in the tree (indexed by node ids; a node is
   associated with the edge that connects it to its parent), and
   column tuples in a "sufficient statistics" representation of the
   alignment (all quantities will be the same for all instances of a
   column tuple). */
struct tp_struct {
  double ****base_probs;        /* posterior probability of each base
                                   given a node, a column
                                   tuple, and a rate category.  First
                                   index is rate category, second is
                                   base, third is node, fourth is
                                   column tuple */
  double *****subst_probs;      /* posterior probability of a
                                   substitution of each base for each
                                   other, given a branch, column
                                   tuple, and rate category.  First
                                   index is rate category, second is
                                   original base, third is replacement
                                   base, fourth is branch, fifth is
                                   column tuple */
  double ***expected_nsubst;    /* expected number of subst for each
                                   branch x column tuple, given a rate
                                   category (conditioned on rate
                                   category in case posterior
                                   probabilities of rate categories
                                   depend on an HMM or similar).
                                   First index is rate category,
                                   second is branch, third is column
                                   tuple */ 
  double ****expected_nsubst_tot; 
                                /* total expected number of
                                   substitutions of each type along
                                   each branch for each rate category,
                                   summed over all column tuples
                                   (considering the number of
                                   instances of each tuple).  These
                                   are the sufficient statistics for
                                   computing the likelihood of a tree
                                   model.  Note that they are based on
                                   *joint* probabilities with rate
                                   categories, rather than being
                                   conditioned on rate categories (the
                                   posterior probability of each rate
                                   at each site in incorporated).
                                   First index is rate category,
                                   second is original base, third is
                                   replacement base, fourth is
                                   branch  */
  double *****expected_nsubst_col;
                                /* expected number of substitutions of each
                                   type along each branch for each rate 
				   category, for each tuple column.
				   First index is rate category, second
				   is branch, third is tuple, fourth is
				   original base, fourth is replacement base 
				*/
  double **rcat_probs;          /* posterior probability of each rate
                                   category for each column tuple.
                                   First index is rate category,
                                   second is column tuple */
  double *rcat_expected_nsites; /* expected number of sites in each
                                   rate category */
};

typedef struct tp_struct TreePosteriors;
                                /* see incomplete type in tree_model.h */

#define NULL_LOG_LIKELIHOOD 1   /* safe val for null when dealing with
                                   log likelihoods (should always be <= 0) */

void tl_dump_matrices(TreeModel *mod, double **inside_vals, 
                      double **outside_vals, double **posterior_probs);
double tl_compute_log_likelihood(TreeModel *mod, MSA *msa, 
                                 double *col_scores, int cat,
                                 TreePosteriors *post);
TreePosteriors *tl_new_tree_posteriors(TreeModel *mod, MSA *msa, int do_bases, 
                                       int do_substs, int do_expected_nsubst, 
                                       int do_expected_nsubst_tot,
				       int do_expected_nsubst_col,
                                       int do_rate_cats, int do_rate_cats_exp);
void tl_free_tree_posteriors(TreeModel *mod, MSA *msa, TreePosteriors *tp);
double tl_compute_partial_ll_suff_stats(TreeModel *mod, TreePosteriors *post);
double tl_compute_ll_suff_stats(TreeModel *mod, MSA *msa, TreePosteriors *post);
Vector *get_marginal_eq_freqs (char *alphabet, int tuple_size, 
                                   Vector *eq_freqs);


#endif
