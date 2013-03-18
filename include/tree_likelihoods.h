/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file tree_likelihoods.h
    Computation of likelihoods for columns of a given multiple
    alignment, according to a given tree model.
    @ingroup phylo
 */

#ifndef TREE_LIK_H
#define TREE_LIK_H

#include <tree_model.h>
#include <msa.h>
#include <math.h>
#include <misc.h>

/** Structure for information related to posterior probability of tree
   model wrt an alignment.  
     Each array is indexed as appropriate for
   rate categories, bases in the alignment (indexed according to a
   model's inv_states; with higher order models, actually tuples of
   bases), nodes or edges in the tree (indexed by node ids; a node is
   associated with the edge that connects it to its parent), and
   column tuples in a "sufficient statistics" representation of the
   alignment (all quantities will be the same for all instances of a
   column tuple). */
struct tp_struct {
  double ****base_probs;        /**< Posterior probability of each base
                                   given a node, a column
                                   tuple, and a rate category.  
					- First index is rate category
					- Second is base
					- Third is node
					- Fourth is column tuple
				 */
  double *****subst_probs;      /**< Posterior probability of a
                                   substitution of each base for each
                                   other, given a branch, column
                                   tuple, and rate category.  
					- First index is rate category
					- Second is original base
					- Third is replacement base
					- Fourth is branch
					- Fifth is column tuple 
				*/
  double ***expected_nsubst;    /**< Expected number of substitutions for each
                                   branch x column tuple, given a rate
                                   category (conditioned on rate
                                   category in case posterior
                                   probabilities of rate categories
                                   depend on an HMM or similar).
                                   	- First index is rate category
					- Second is branch 
					- Third is column tuple 
				*/ 
  double ****expected_nsubst_tot; 
                                /**< Total expected number of
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
                                   	- First index is rate category
					- Second is original base
					- Third is replacement base
					- Fourth is branch 
 */
  double *****expected_nsubst_col;
                                /**< Expected number of substitutions of each
                                   type along each branch for each rate 
				   category, for each tuple column.
				   	- First index is rate category
					- Second is branch 
					- Third is tuple
					- Fourth is original base
					- Fifth is replacement base 
				*/
  double **rcat_probs;          /**< Posterior probability of each rate
                                   category for each column tuple.
                                    	- First index is rate category
					- Second is column tuple 
				*/
  double *rcat_expected_nsites; /**< Expected number of sites in each
                                   rate category */
};

typedef struct tp_struct TreePosteriors;
                                /* see incomplete type in tree_model.h */

#define NULL_LOG_LIKELIHOOD 1   /** Safe value for null when dealing with
                                   log likelihoods (should always be <= 0) FIXME? */

/* does not appear to be implemented */
void tl_dump_matrices(TreeModel *mod, double **inside_vals, 
                      double **outside_vals, double **posterior_probs);

/** Compute the likelihood of a tree model with respect to an
   alignment; Optionally retain column-by-column likelihoods and/or posterior probabilities.  
   @param[in] mod Tree Model to compute likelihood for
   @param[in] msa Multiple Alignment containing data related to tree model
   @param[out] col_scores (Optional) Log likelihood score per column
   @param[out] tuple_scores (Optional) Log likelihood score per tuple
   @param[in] cat Whether to use categories
   @param[out] post (Optional) Computed posterior probabilities; If NULL, no
   posterior probabilities (or related quantities) will be computed.
   If non-NULL each of its attributes must either be NULL or
   previously allocated to the required size. 
   @result Log likelihood of entire tree model specified
*/
double tl_compute_log_likelihood(TreeModel *mod, MSA *msa, 
                                 double *col_scores, 
				 double *tuple_scores, 
				 int cat,
                                 TreePosteriors *post);

/** Create a new TreePosteriors object.
    @param mod Tree Model of which the posterior probabilities are calculated
    @param msa Multiple Alignment
    @param do_bases Whether to allocate space for base probabilities
    @param do_subst Whether to allocate space for substitution probabilities
    @param do_expected_nsubst Whether to allocate space for expected number of substitutions matrix
    @param do_expected_nsubst_tot Whether to allocate space for total expected number of substitutions
    @param do_expected_nsubst_col Whether to allocate space for expected number of substitutions per column
    @param do_rate_cats Whether to allocate space for rate categories
    @param do_rate_cats_exp Whether to allocate space for expected rate categories
    @result Newly allocated TreePosteriors object
*/
TreePosteriors *tl_new_tree_posteriors(TreeModel *mod, MSA *msa, int do_bases, 
                                       int do_substs, int do_expected_nsubst, 
                                       int do_expected_nsubst_tot,
				       int do_expected_nsubst_col,
                                       int do_rate_cats, int do_rate_cats_exp);

/** Free TreePosteriors object
   @param mod Tree model of which posterior are calculated
   @param msa Multiple Alignment
   @param tp TreePosteriors object to free
 */
void tl_free_tree_posteriors(TreeModel *mod, MSA *msa, TreePosteriors *tp);

/** Compute the expected (posterior) complete log likelihood of a tree
   model based on a TreePosteriors object.  
   @param[in] mod Tree Model
   @param[in] post Pre-calculated posterior probabilities
   @note Equilibrium frequencies are not considered
   @result Log Likelihood of tree
*/
double tl_compute_partial_ll_suff_stats(TreeModel *mod, TreePosteriors *post);

/* Could not find implementation */
double tl_compute_ll_suff_stats(TreeModel *mod, MSA *msa, TreePosteriors *post);

/** Given an alphabet, a tuple size, and a vector of equilibrium
   frequencies, create a new vector of marginal equilibrium
   frequencies describing the space of "meta-tuples", which contain
   actual characters *or* missing data characters.  
   Each meta-tuple is
   given an equilibrium frequency equal to the sum of the frequencies
   of all "matching" ordinary tuples.  
    Missing data characters are
   assumed to be gap characters or Ns. 
   @param alphabet List of possible characters
   @param tuple_size Size of tuples
   @param eq_freqs Equilibrium frequencies
   @param New vector of marginal equilibrium frequencies
*/
Vector *get_marginal_eq_freqs (char *alphabet, int tuple_size, 
                                   Vector *eq_freqs);


#endif
