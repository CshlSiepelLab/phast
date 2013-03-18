/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file subst_mods.h
    Functions to setup substitution models
    @ingroup phylo
*/

#ifndef SUBST_MODS_H
#define SUBST_MODS_H

#include <markov_matrix.h>
#include <vector.h>
#include <stringsplus.h>

/** Substitution models available */
typedef enum {
  JC69,	/**< Jukes Cantor 1969 */
  K80,  /**< Kimura 1980 */
  F81,  /**< Felsenstein 1981 */
  HKY85, /**< Hasegawa, Kishino and Yano 1985 */
  HKY85G, /**< Hasegawa, Kishino and Yano (Gamma) 1985 */
  REV,   /**< Reversible model*/
  SSREV,  /**< Strand-symmetric Reversible  */
  UNREST, /**< Unrest */
//  HKY2,
  R2,	  /**< Reversible di-nucleotide  */
  U2,	  /**< Unrestricted di-nucleotide */
  R2S,    /**< Reversible di-nucleotide strand-symmetric */
  U2S,	  /**< Unrestricted di-nucleotide strand-symmetric */
  R3,	  /**< Reversible tri-nucleotide */
  R3S,    /**< Reversible tri-nucleotide strand-symmetric */
  U3,     /**< Unrestricted tri-nucleotide */
  U3S,    /**< Unrestricted tri-nucleotide strand-symmetric */
  GC,	  /**< GC */
  //  HB,
  HKY_CODON, /**< Hasegawa, Kishino and Yano (Codon)*/
  REV_CODON, /**< Reversible (Codon) */
  SSREV_CODON,  /**< Strand-symmetric reversible (Codon) */
  UNDEF_MOD   /**< No Model */
} subst_mod_type;

struct tm_struct;               /* use incomplete type because of
                                   reciprocal dependencies with tree_model.h */
/** \name Setup probability matrix
\{ */
/** Setup probability matrix for JC69 
    @param mod Tree Model
    @param p Markov Matrix to set the probabilities for
    @param t T parameter (dparent * branch_scale * Tree Models rK)
    @note For an example of how to calculate 't', see tm_set_subst_matrices
    @see tm_set_subst_matrices
*/
void tm_set_probs_JC69(struct tm_struct *mod, MarkovMatrix *P, double t);

/** Setup probability matrix for F81 
    @param backgd_freqs Background frequencies (Usually from a Tree Model)
    @param P Markov Matrix to set probabilities for
    @param scale Scaling constant (1/(1-sum(backgd_freqs^2))) 
    @param t T parameter  (dparent * branch_scale * Tree Models rK)
    @note For an example of how to calculate 'scale' and 't', see tm_set_subst_matrices
    @see tm_set_subst_matrices
*/
void tm_set_probs_F81(Vector *backgd_freqs, MarkovMatrix *P, double scale, double t);

/** Setup probability matrix by copying an existing probability matrix.
    Set matrix such that element (i,j) has value pi_j, as for an
     infinitely long branch
    @param mod Tree Model to setup probability matrix for
    @param P Existing markov matrix to copy probabilities from
 */
void tm_set_probs_independent(struct tm_struct *mod, MarkovMatrix *P);

/** Return the substitution model (enum val) corresponding to the
   specified string.
   @param str Substitution Model as string
   @result Substitution Model as enumerated value of type subst_mod_type
 */
subst_mod_type tm_get_subst_mod_type(const char *str);

/** \} */

/**  Return a string description for the specified subst_mod_type.
   @param type Substitution Model type
   @result Name of substitution model as a string
 */
char *tm_get_subst_mod_string(subst_mod_type type);

/** Return number of rate matrix parameters (not counting eq. freqs) 
   @param mod Tree Model with a defined substitution model and rate_matrix
   @result Number of rate matrix parameters 
   @note Some substitution models do not need to have an allocated rate matrix for this function to work
*/
int tm_get_nratematparams(struct tm_struct *mod);

/** Get the order of a substitution model
    @param subst_mod subst_mod_type Substitution Model i.e. R2 or U2S
    @result Order of the substitution model specified
    @note Although codon models (such as HKY_CODON, REV_CODON, and SSREV_CODON) are technically 0th order models representing 3 bases, tm_order returns 2 for these models.
 */
int tm_order(int subst_mod);

/** Test whether the substitution model specified handles codons
    @param subst_mod subst_mod_type Substitution Model i.e. R2 or U2S
    @result 1 if subst_mod supports codons, otherwise 0
 */
int subst_mod_is_codon_model(int subst_mod);

/** \name Initialize rate parameters
 \{ */
/** Initialize rate-matrix parameters in parameter vector, using an
   HKY-like strategy 
   @param mod Tree Model containing rate matrix to initialize
   @param params Parameter vector
   @param kappa Defines transition/transversion bias
   @param params_idx Starting index of vector params
   @see tm_set_rate_matrix
*/
void tm_rate_params_init(struct tm_struct *mod, Vector *params, 
                         int params_idx, double kappa);
/** Initialize rate-matrix parameters in parameter vector, based on an
   existing model. 
   @param mod Tree Model containing rate matrix to initialize
   @param params Parameter vector
   @param params_idx Starting index of vector params
   @param selection Selection factor 
   @param bgc Bias gene conversion factor
   @see tm_set_rate_matrix_sel_bgc
 */
void tm_rate_params_init_from_model(struct tm_struct *mod, Vector *params, 
                                    int params_idx, 
				    double selection, double bgc);
 /** \} \name Set rate matrix 
\{ */

/** Initialize rate-matrix parameters in Models rate-matrix
    @param mod Tree Model containing rate-matrix to initialize
    @param kappa Interacts with background frequency in case of transition to set rate-matrix values
    @param kappa_idx Index of rate_matrix_param_row to setup mapping at
 */
void tm_set_HKY_matrix(struct tm_struct *mod, double kappa, int kappa_idx);


/** Set rate matrix according to elements of parameter vector
    @param mod Tree Model containing rate matrix to set
    @param params Parameter vector containing elements used to set rate matrix
    @param i Starting index
    @note Neither JC69 nor F81 use the parameters 'params' and 'i'
   starting index)
*/
void tm_set_rate_matrix(struct tm_struct *mod, Vector *params, int i);

/** Set rate matrix according to elements of parameter vector; Then set bias gene conversion and selection factors
    @param mod Tree Model containing rate matrix to set
    @param params Parameter vector containing elements used to set rate matrix
    @param i Starting index
    @param selection Selection factor
    @param bgc Bias gene conversion 
    @note Neither JC69 nor F81 use the parameters 'params' and 'i'
   starting index)
*/
void tm_set_rate_matrix_sel_bgc(struct tm_struct *mod, Vector *params, int i,
				double selection, double bgc);

 /** \} */

/* Couldn't find implementation */
int tm_substmod_get_nratematparams(subst_mod_type submod, 
				   struct tm_struct *mod);
/** Find the position(s) of substitution parameter used by model.
    @param mod Tree Model with substitution specified
    @param flag Array of int size of maximum number of matrix parameters
    @param param_name Parameter name to find position for.
    @result 1 on success, 0 on error
*/
int tm_flag_subst_param_pos(struct tm_struct *mod, int *flag, 
			    String *param_name);

/** \name Apply / Remove selection & bias gene conversion (bgc) factors
\{ */

/** Apply selection factor and bias gene conversion factor to either 4state or codon Markov Matrix 
    @param mm Markov Matrix to apply selection and bgc factors to
    @param sel Selection Factor
    @param bgc Bias gene conversion factor
*/
void tm_apply_selection_bgc(MarkovMatrix *mm, double sel, double bgc);

/** Remove selection factor and bias gene conversion factor on either 4state or codon Markov MAatrix
    @param mm Markov Matrix to remove selection and bgc factors from
    @param sel Selection Factor to remove
    @param bgc Bias gene conversion factor to remove
*/
void tm_unapply_selection_bgc(MarkovMatrix *mm, double sel, double bgc);
 /** \} */
/** \name Misc
\{ */

/** Get the substitution model with same meaning/parameterization but which corresponds to codons
    @param subst_mod A nucleotide substitution model
    @return The "codon version" of subst_mod, which has the same parameterization but has been expanded from 4x4 to 64x64, or UNDEF_MOD if none exists
    @note if subst_mod is a codon model, returns subst_mod
    @note Most substitution models do not have codon version.  Currently only HKY85, REV, and SSREV do.  Prints a warning if result is UNDEF_MOD.
 */
subst_mod_type tm_codon_version(subst_mod_type subst_mod);


/** Get the substitution model with the same meaning/parameterization but which corresponds to nucleotides
   @param subst_mod A codon substitution model
   @return The "nucleotide version" of subst_mod, which has the same parameterization but uses a 4x4 matrix rather than 64x64
   @note if subst_mod is a nucleotide model, returns itself
   @note Prints a warning if result is UNDEF_MOD.
 */
subst_mod_type tm_nucleotide_version(subst_mod_type subst_mod);
/** \} */
#endif
