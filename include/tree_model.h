/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file tree_model.h
    @brief A tree model represents a phylogenetic tree, substitution rate matrix,
    and background frequencies.  
    The model allows for rate variation and
    also for varying substitution models on different branches.  If the
    tree model is optimized by maximum likelihood, the tree model object 
    contains data which indicate which parameters to hold constant, and 
    which to optimize, as well as boundary conditions.
    @ingroup phylo
*/

#ifndef TREE_MODEL_H
#define TREE_MODEL_H

#include <vector.h>
#include <markov_matrix.h> 
#include <trees.h> 
#include <stringsplus.h> 
#include <msa.h> 
#include <misc.h>
#include <numerical_opt.h>
#include <subst_mods.h>
#include <hmm.h>

/** Maximum alphabet size */
#define MAX_ALPH_SIZE 100

/** In cases where a function of complex numbers is supposed to be
    real, we'll consider it close enough if the imaginary component is
    smaller in magnitude than this threshold 
*/
#define TM_IMAG_EPS 1e-6           
                           

/* convergence thresholds for EM, expressed as portions of function
   value: (f(xold) - f(xnew)) / f(xnew) */
/** High convergence threshold for EM  (as portion of function value (f(xold) - f(xnew)) / f(xnew))*/
#define TM_EM_CONV_HIGH 1e-8
/** Medium convergence threshold for EM  (as portion of function value (f(xold) - f(xnew)) / f(xnew))*/
#define TM_EM_CONV_MED 5e-7
/** Low convergence threshold for EM  (as portion of function value (f(xold) - f(xnew)) / f(xnew))*/
#define TM_EM_CONV_LOW 1e-5
/** Crude convergence threshold for EM  (as portion of function value (f(xold) - f(xnew)) / f(xnew))*/
#define TM_EM_CONV_CRUDE 1e-4

/** Test of conversion */
#define TM_EM_CONV(P) ( (P) == OPT_HIGH_PREC ? TM_EM_CONV_HIGH : ( (P) == OPT_MED_PREC ? TM_EM_CONV_MED : ( (P) == OPT_CRUDE_PREC ? TM_EM_CONV_CRUDE : TM_EM_CONV_LOW) ))

#define BACKGD_STR "backgd"
#define RATEMAT_STR "ratematrix"
#define RATEVAR_STR "ratevar"
#define BRANCHES_STR "branches"
#define SCALE_STR "scale"
#define SCALE_SUB_STR "scale_sub"
#define SELECTION_STR "sel"
#define BGC_STR "bgc"

/** Type of branch length estimation */
typedef enum { 
  TM_BRANCHLENS_ALL, 		/**< Estimate branch lengths */
  TM_BRANCHLENS_CLOCK,          /**< Molecular clock */
  TM_SCALE_ONLY, 		/**< Scale factor only */
  TM_BRANCHLENS_NONE 		/**< Do not estimate branch lengths */
} blen_estim_type;

/** Type of parameter bound; used in bounded estimation of subtree scale */
typedef enum {
  NB,                           /**< No bound */
  LB,                           /**< Lower bound -- subtree scale at
                                   least as large as scale of whole
                                   tree */
  UB                            /**< Upper bound -- subtree scale at
                                   most as large as scale of whole
                                   tree */
} scale_bound_type; 

struct tp_struct;


/** Defines alternative substitution model for a particular branch */
typedef struct alt_subst_mod {
  subst_mod_type subst_mod;    /**< Substitution model */
  Vector *backgd_freqs;        /**< Eq freqs (set to NULL if separate_backgd=0) */
  MarkovMatrix *rate_matrix;   /**< Rate matrix (set to NULL if separate_model=0) */
/* Indices in main model parameters where lineage-specific paramers start */  
  int ratematrix_idx,		/**< Indices of rate matrix model parameters */ 
  backgd_idx, 			/**< Index of first background frequency parameter */
  selection_idx, 		/**< Index of selection parameter */
  bgc_idx;			/**< Index of bias gene conversion parameter */
  double selection,		/**< Optional selection parameters for this model (only used if selection_idx >=0 ) .  Total selection on the model is the sum of this parameter and the selection value in the main model.*/
  bgc;       			/**< Optional bgc parameters for this model (only used if bgc_idx >=0 ) */
  int separate_model;   /**< ==1 if no parameters shared with main model */
  int separate_backgd;  /**< ==1 if no background frequencies shared with non-alternate substitution model */
  List *param_list;     /**< List of string arguments giving which parameters to
			   estimate separately (only if separate_model=0) */
  String *defString;    /**< Name of alternative Model. This is the argument given to phyloFit
			   to define the alternative model */
  String *noopt_arg;    /**< Which parameters to hold constant (represented as a command line argument) */
  struct alt_subst_mod *share_sel, /**< If not NULL, share selection parameter with another altSubstMod */
    *share_bgc;  /**< If not NULL, share bgc parameter with another altSubstMod */
  /** Note: could implement sharing with other parameters, hasn't been needed yet */
} AltSubstMod;



/** Tree model object */
struct tm_struct {
  TreeNode *tree;		/**< Root node of tree (used to traverse tree node by node) */
  Vector *backgd_freqs;		/**< Equilibrium frequencies */
  MarkovMatrix *rate_matrix;	/**< Rate matrix */
  subst_mod_type subst_mod;	/**< Substitution model being used for this tree model */
  int *msa_seq_idx;             /**< (Optional) Mapping from leaf
                                   indices to sequence indices in a
                                   given MSA; used in likelihood
                                   computations */
  MSA *msa;                     /**< (Optional) MSA from which TreeModel
                                   was estimated; for use in tm_fit */
  int category;                 /**< (Optional) site category in MSA or -1 */
  int order;                    /**< Order of model */
  double alpha;                 /**< For gamma or discrete-gamma rate
                                   variation; set to 0 for constant rate */
  int nratecats;                /**< Number of rate categories, discrete gamma */
  double lnL;                   /**< Log likelihood from training */
  struct tp_struct *tree_posteriors; 
                                /**< (Optional) associated
                                   TreePosteriors object */
  int use_conditionals;         /**< (Optional) compute likelihood using
                                   conditional probabilities at each
                                   column; only relevant when order >
                                   0 */
  MarkovMatrix ***P;            /**< Probability matrices for edges,
                                   indexed by node id and rate category */
  double *rK, *freqK;           /**< Rate constants and frequencies */
  List **rate_matrix_param_row, /**< Rate matrix parameter row */
  **rate_matrix_param_col;	/**< Rate matrix parameter column */
  int root_leaf_id;             /**< Indicates id of leaf node to be
                                   Interpreted as the root.  Must be
                                   a child of the actual root of the
                                   tree.  Ordinarily will have a null
                                   value (-1), but if non-negative, the
                                   distance to its parent will be
                                   constrained to be zero.  */
  int allow_gaps;		/**< If TRUE, gaps are not allowed to be
				   taken into account for a model */
  int allow_but_penalize_gaps;  /**< If TRUE, gaps are allowed but are
				   penalized NOTE: not used */
  int inform_reqd;              /**< If TRUE, only "informative" sites
                                   will be given non-zero probability */
  int estimate_backgd;          /**< Estimate background frequencies as free
                                   parameters in the optimization */
  blen_estim_type estimate_branchlens; 
                                /**< If value is TM_BRANCHLENS_ALL, then
                                   estimate all branch lengths; if
                                   value is TM_SCALE_ONLY, then
                                   estimate only a scale factor; and
                                   if value is TM_BRANCHLENS_NONE, do
                                   not estimate branch lengths */
  double scale;                 /**< Current scale factor */
  double scale_sub;             /**< Scale factor for subtree, when
                                   estimating separate scale factors
                                   for subtree and supertree */
  double selection;
  int *in_subtree;              /**< Array indicating whether each
                                   branch is in designated subtree */
  scale_bound_type scale_sub_bound;
                                /**< Bound on scale of subtree */
  TreeNode *subtree_root;       /**< Node defining subtree */
  int *ignore_branch;           /**< If ignore_branch is non-NULL and
                                   ignore_branch[i] == TRUE, then
                                   branch i will be ignored (treated
                                   as infinitely long) in likelihood
                                   calculations */
  int empirical_rates;          /**< Indicates "empirical"
                                   (non-parametric) model for
                                   rate-variation */
  int site_model;               /* TRUE if this is a Nielsen-Yang site model; indicates parameterization for site categories */
  int estimate_ratemat;         /**< Indicates whether rate-matrix
                                   parameters should be estimated */
  AltSubstMod ***alt_subst_mods_ptr; /**< Pointer to alt_subst_mod for 
				       each branch and each rate category */

  List *alt_subst_mods;         /**< List of relevant AltSubstMods for 
                                   this tree.  alt_subst_mods_ptr above are 
				   usually pointers to elements in this 
				   list */
  Vector *all_params;           /**< All parameters relevant to this model */
  //  Vector *lowbound, *upbound;
  int *param_map;               /**< Map of parameter indices used by 
				   phyloFit.  If param_map[i]==-1, then 
				   parameter i is held constant.  Otherwise,
				   if param_map[i]==j, then it is the j'th
				   element in the vector of parameters which
				   are optimized */
 /* These are the indices in all_params that show where scale, branchlens, rateMatrix,
    backgd_freqs, and rate variation parameters are stored*/
  int scale_idx,		/**< Starting index of scale parameters in parameter vector*/ 
  bl_idx, 			/**< Starting index of branch length parameters in parameter vector */
  ratematrix_idx, 		/**< Starting index of rate matrix info in parameter vector */
  backgd_idx, 			/**< Starting index of background frequency info in parameter vector */
  ratevar_idx,			/**< Starting index of rate variation parameters in parameter vector */
  selection_idx;		/**< Starting index of selection parameters in parameter vector */
  List *bound_arg;              /**< Used by phyloFit, this is a copy of the
				   command-line argument(s) which specify
				   boundaries on parameters */
  String *noopt_arg;            /**< Used by phyloFit, this is a copy of the
				   command-line argument which specifies
				   which parameters to hold constant */
  int eqfreq_sym;               /**< Use symmetrical equilibrium frequencies? */
  int scale_during_opt;       /**< Whether to scale rate matrix during optimization.
				 Normally 0, but 1 if TM_BRANCHLENS_NONE, or
				 if TM_SCALE and alt_subst_mods!=NULL */
  int **iupac_inv_map;          /**< Inverse map for IUPAC ambiguity characters */
};

typedef struct tm_struct TreeModel;

/** \name Tree Model allocation/cleanup functions 
\{ */
/** Create new tree model object.*/
TreeModel *tm_new(TreeNode *tree,  /**< Phylogenetic tree */
		  MarkovMatrix *rate_matrix, /**< Matrix of transition/transversion rates */
                  Vector *backgd_freqs, /**< Background frequencies */
		  subst_mod_type subst_mod,  /**< Substitution Model i.e. REV (reversable)*/
                  char *alphabet,  /**< List of possible characters in the sequence data i.e. 'ATGC' */
		  int nratecats, 		/**< Number of rate categories */
		  double alpha,			/**< Alpha parameter */
                  List *rate_consts, /**< List defining rate constants using
                  	  	  	  	  	  non-parametric mixture model for rates instead of gamma dist. */
		  int root_leaf_id     /**< ID number of the root node */
		  );

/** Re-initialize a tree model with an altered substitution model,
    number of rate categories, alpha, set of rate consts, or set of
    rate weights.
    @param tm Tree Model to re-initialize
    @param subst_mod New substitution model
    @param new_nratecats New number of rate categories
    @param new_alpha  New alpha parameter (ignored if new_nratecats == 1)
    @param new_rate_consts New list of explicit rate
       constants.  If non-NULL, implies use of empirical rates.
    @param new_rate_weights  New list of explicit rate weights
       (mixing proportions).  Will be normalized automatically.  If
       new_rate_consts != NULL and new_rate_weights == NULL, weights
       will be initialized to approx gamma distrib. defined by alpha
 */
void tm_reinit(TreeModel *tm, subst_mod_type subst_mod,
               int new_nratecats, double new_alpha, 
               List *new_rate_consts, List *new_rate_weights);


/**  Set up lists that allow each rate matrix parameter to be mapped to
   the rows and columns in which it appears.
   @param tm Tree Model that contains rate matrix
    */
void tm_init_rmp(TreeModel *tm);


/** Free rate matrix parameter
   @param tm Tree Model that contains rate matrix parameter */
void tm_free_rmp(TreeModel *tm);

/** Free tree model
   @param tm Tree Model to free */
void tm_free(TreeModel *tm);



/** \} */
/** \name Tree Model file access functions 
\{ */

/** Read tree model from file.
   @param F File descriptor of file containing tree model
   @param discard_likelihood Whether to discard the log likelihood saved in the file
   @result Newly allocated TreeModel from file
   @note Simple format specifies: Background
   rates, rate matrix, and tree (topology and branch lengths) */
TreeModel *tm_new_from_file(FILE *F, int discard_likelihood);

/** Tree Model to save to file
   @param F File descriptor to save Tree Model to
   @param tm Tree model to save to file
*/
void tm_print(FILE *F, TreeModel *tm);

/** \} */

/** Create a copy of TreeModel with an alternative substitution model.
   @param tm Tree Model to copy from
   @param altmod_str Alternate substitution model */
AltSubstMod* tm_add_alt_mod(TreeModel *tm, String *altmod_str);

/** Tree Model to save to file
   @param F File descriptor to save Tree Model to
   @param tm Tree model to save to file
*/
void tm_print(FILE *F, TreeModel *tm);

/** Copy a Tree Model to an existing tree model
    @param dest Destination to copy tree model to
    @param src Source to copy tree model from
*/
void tm_cpy(TreeModel *dest, TreeModel *src);

/** Copy a tree model into a new TreeModel object.
   @param src Tree Model source to copy contents from
   @result Newly allocated TreeModel containing data from src */
TreeModel *tm_create_copy(TreeModel *src);

/** \name Tree Model substitution matrix functions 
\{ */

/** Setup the substitution matrices on a Tree Model.
   @param tm Tree Model to setup substitution matrix for
*/
void tm_set_subst_matrices(TreeModel *tm);

/** Setup the substitution matrices on a Tree Model with custom probability matrix and branch length.
   @param P Probability matrix to use
   @param t Branch length to use
 */
void tm_set_subst_matrix(TreeModel *tm, MarkovMatrix *P, double t);


/** Create an alternate substitution model with new background frequencies and probabilities.
   @param subst_mod Existing Substitution Model to base new model on
   @param backgd_freqs Background frequencies for new model
   @param rate_matrix Markov Matrix defining  substitution probabilities
   @result Newly created, alternate substitution method
*/
AltSubstMod* tm_new_alt_subst_mod(subst_mod_type subst_mod,
                                  Vector *backgd_freqs, 
                                  MarkovMatrix *rate_matrix);

/** Free an alternate substitution model.
   @param am Alternate substitution model to free
*/
void tm_free_alt_subst_mod(AltSubstMod *am); 


/** \} \name Tree Model scale functions 
\{ */

/** Scales rate matrix and Optionally parameters, branch lengths, substitution matrices.
   @param tm Tree Model containing rate matrix to scale
   @param params (Optional) Parameters to scale
   @param scale_blens Whether to scale branch lengths as well
   @param reset_subst_matrices Whether to reset substitution matrices as well
   @warning Works with lineage-specific models (although this is not necessarily the proper way to calculate scale factor!)
   */
void tm_scale_model(TreeModel *tm, Vector *params, int scale_blens,
		    int reset_subst_matrices);


/** Allocate (if necessary) and initialize backgd frequencies based on observed frequencies in MSA.
    If model is symmetric (SSREV) then makes backgd frequencies symmetric.
    If mod->backgd_freqs is already allocated, it is not reallocated and assumed to be the correct size.
    @param mod A Tree model object.  mod->backgd_freqs will be set by this function.
    @param msa An MSA object
    @param cat The category of the msa to use to compute the backgd
 */
void tm_init_backgd(TreeModel *mod, MSA *msa, int cat);


/** Modifies equilibrium frequency of a model in such a way that
   reversibility is maintained.
   @param tm Tree Model containing equilibrium frequency to change
   @param newfreqs New equilibrium frequencies for tree model
*/
void tm_mod_freqs(TreeModel *tm, Vector *newfreqs);

/** Scale evolutionary rate (only effects branch lengths)
  @param tm Tree Model containing branch lengths to scale
  @param scale_const Constant to scale by
  @param reset_subst_mats Whether to reset substitution matrices
  @note to scale rate matrix and more see tm_scale_model
  @see tm_scale_model
*/
void tm_scale_branchlens(TreeModel *tm, double scale_const, int reset_subst_mats);

/** Scale the rate matrix of the specified tree model such that the
   expected rate of substitution at equilibrium is 1.  This is the
   standard convention; it allows the unit of measurement of branch
   lengths to be expected substitutions/site.
   @param mod Tree Model to scale
   @result Scaling factor
   @note Scaling factor can be used to adjust
   branch lengths if the rate matrix has not been scaled throughout the
   fitting procedure. */
double tm_scale_rate_matrix(TreeModel *mod);

/** Scale a parameter vector according to a specified rate matrix scale
   factor.  Branch length params are multiplied by the specified factor
   and rate matrix params by its inverse.
   @param mod Tree Model
   @param params Parameters for scaling
   @param scale_factor Factor to scale parameters by
*/
void tm_scale_params(TreeModel *mod, Vector *params, double scale_factor);


/** \} */

/** Fit a tree model to data using BFGS (multidimensional optimization algorithm)
   @param mod Tree Model containing desired tree topology to fit, substitution model, and (if appropriate) background frequencies
   @param params Initial values for optimization procedure
   @param cat MSA category
   @param precision Precision describing BFGS convergence criteria
   @param logf File descriptor of log file
   @param quiet Whether to report progress to stderr
   @param error_file If non-NULL, write estimate, variance, and 95% 
   confidence interval for each parameter to this file.
   @returns 0 on success, 1 on failure
 */
int tm_fit(TreeModel *mod, MSA *msa, Vector *params, int cat, 
           opt_precision_type precision, FILE *logf, int quiet,
	   FILE *error_file);


/** Fit several tree models (which share parameters) to data using BFGS
    @param mod Array of tree models
    @param nmod Length of mod array
    @param msa Array of msas.  There should be one for each mod, or a single msa with nmod categories.  In that case mod[i] will apply to category i+1 (i in 0,...,nmod-1).
    @param nmsa Length of msa array.
    @param precision Precision describing BFGS convergence criteria
    @param logf log file
    @param quiet Whether to report progress to stderr
    @returns 0 on success, 1 on failure
 */
int tm_fit_multi(TreeModel **mod, int nmod, MSA **msa, int nmsa,
		 opt_precision_type precision,
		 FILE *logf, int quiet);

/** Set specified TreeModel according to specified parameter vector.
   Exact behavior depends on substitution model.
   @param mod Tree Model to adjust parameter vector for
   @param params Parameter vector values
   @param idx_offset (Optional) Index offset for cases in which vectors of parameters are
   nested within larger vectors of parameters (set to -1 if not
   needed) */
void tm_unpack_params(TreeModel *mod, Vector *params, int idx_offset);

/** \name Tree Model parameter initialization 
\{ */

/** Initializes all branch lengths and rate matrix parameters.
   Initializes branch lengths to designated constant; initializes
   rate-matrix params using specified kappa; kappa ignored for JC69;
   REV and UNREST initialized as if HKY.  Initializes alpha as
   specified, if dgamma.  In the case of empirical rates, uniform
   weights are used for initialization.
   @param mod Tree Model to initialize parameters for
   @param branchlen Constant to use as initial branch length
   @param kappa Rate matrix initialization parameter
   @param alpha Initial value for alpha (if dgamma)
   @result Parameter vector
   */
Vector *tm_params_init(TreeModel *mod, double branchlen, double kappa,
                           double alpha);
/** Initializes all branch lengths and rate matrix parameters using random values.
   @param mod Tree Model to initialize parameters for
   @result Parameter vector
*/
Vector *tm_params_init_random(TreeModel *mod);

/** Functions to initialize a parameter vector from an existing tree model
   @param mod Tree Model used to generate parameter vector
   @result Parameter vector
 */
Vector *tm_params_new_init_from_model(TreeModel *mod);

/** Functions to initialize a parameter vector from an existing tree model
   @param[in] mod Tree Model used to generate parameter vector
   @param[out] params Parameter vector
   @param[in] setup_params If TRUE, calls tm_setup_params to set scale_idx, bl_idx, backgd_idx, etc to determine which parameters are stored where in the parameter vector.  Otherwise assume this has already been done.
 */
void tm_params_init_from_model(TreeModel *mod, Vector *params, int setup_params);

/** Initialize branch lengths to average number of mutations under parsimony, and apply Jukes Cantor correction.
  @param params (Optional) If params is NULL, sets the branch lengths in mod->tree directly. Otherwise sets the parameter vector
  @param mod Tree Model to init branch lengths for
  @param cat Category associated with tree model
  @result Parsimony Score
 */
double tm_params_init_branchlens_parsimony(Vector *params, 
					   TreeModel *mod, MSA *msa, int cat);

/** Setup parameter indexes (mod->scale_idx, mod->ratevar_idx, etc.) and initialize optimization vector. 
   Assigns parameter indices to mod->scale_idx,
   mod->bl_idx, mod->ratevar_idx, mod->ratematrix_idx which point
   to indices in mod->all_params where corresponding parameters are
   stored.  (as well as pointers in altmod structures).
   
   Then, it determines which parameters are to be optimized, and allocates
   and assigns values to mod->param_map, which maps indices in all_params
   to indices in the optimization vector, so that param_map[i]=j implies
   that all_params[i] is stored in opt_vec[j].  Only parameters which
   are to be optimized are stored in the optimization vector.  Parameters
   which are held constant have param_map[i]=-1.  
   @param mod Tree Model
   @param offset The first position to use for parameter indices (usually 0).
   @result The number of parameters in this model
 */
int tm_setup_params(TreeModel *mod, int offset);

/** \} \name Tree Model get counts
\{  */

/** Return number of parameters for specified TreeModel
   @param mod Tree Model
   @note Based on number of taxa and substitution model
   @result Nubmer of parameters for specified tree model
*/
int tm_get_nparams(TreeModel *mod);

/** Get the number of equilibrium frequency parameters.
   @param mod Tree Model
   @result Number of equilibrium frequencies for specified model
*/
int tm_get_neqfreqparams(TreeModel *mod);

/** Get the number of rate variable parameters.
   @param mod Tree Model
   @result Number of rate variable parameters
*/
int tm_get_nratevarparams(TreeModel *mod);

/** Return the number of leaves in a Tree
   @param mod Tree Model to count leafs for
   @result Number of leaves in tree model
 */
int tm_get_nleaf(TreeModel *mod);

/** Return the number of branch length parameters
    @param mod Tree Model
    @result Number of branch length parameters
*/
int tm_get_nbranchlenparams(TreeModel *mod);


/** \} \name Tree Model check reversibility 
\{ */

/** Get whether the model being used is reversible
  @param tm Tree model to test
  @result 1 if reversible, 0 otherwise
*/
int tm_is_reversible(TreeModel *tm);

/** Get whether a node within a tree model is reversible
    @param tm Tree model to test
    @param node Node with alternate substitution model to test
    @result 1 if reversible, 0 otherwise
*/
int tm_node_is_reversible(TreeModel *tm, TreeNode *node);

/** Get whether a substitution model is reversible
    @param subst_mod Substitution Model of type subst_mod_type
    @result 1 if reversible, 0 otherwise
 */
int subst_mod_is_reversible(int subst_mod);

/** \} \name Tree Model generate MSA
\{ */


/** Generates an alignment according to set of Tree Models and a
   Markov matrix defining how to translate among them.
   @pre Call srandom externally
   @pre TreeModels must appear in same order as the states of the Markov matrix.
   @param ncolumns Number of columns in MSA
   @param hmm (Optional) Markov matrix of probabilities if NULL, single tree model assumed
   @param classmods Array of tree models with different classes, one per HMM state (or 1 total if HMM is NULL)
   @param labels (Optional) Used to record state (model) responsible for generating each site; pass NULL if hmm is NULL
   @result Multiple Alignment generated from HMM and Tree Model
   @note Only appropriate for order 0 models
*/
MSA *tm_generate_msa(int ncolumns, HMM *hmm, 
                     TreeModel **classmods, int *labels);

/** Generates a random alignment using a list of scales and subtree scales.
   @pre Call srandom externally
   @param nsitesList List of int with each entry indicating number of sites per scale
   @param scaleLst List of double with each entry indicating a scale
   @param subtreeScaleLst List of double with each entry indicating scale of a subtree
   @param model Tree Model
   @param subtreeName (Optional) Name of node identifying subtree root
   @result Multiple Alignment generated using a list of scales and subtree scales
   @note If subtreeName is NULL subtree scales must be 1
*/
MSA *tm_generate_msa_scaleLst(List *nsitesList, List *scaleLst,
			      List *subtreeScaleLst, TreeModel *model,
			      char *subtreeName);

/** Generates an alignment according to a Tree Model, subtree, and and probability of subtree switching.
   @pre Call srandom externally.
   @param ncolumns Number of columns in MSA
   @param mod Tree Model
   @param subtreeMod Tree Model of the subtree
   @param subtree Name of node identifying subtree root
   @param subtreeSwitchProb Probability that a node switches its state of being in or out of the subtree
   @result Multiple Alignment generated using a Tree Model, subtree, and probability of subtree switching
*/
MSA *tm_generate_msa_random_subtree(int ncolumns, TreeModel *mod,
				    TreeModel *subtreeMod, 
				    char *subtree, double subtreeSwitchProb);

/** \} */

/** Given a codon model, create and return the induced amino acid model
  @param codon_mod Codon model
  @param Amino Acid model
*/
TreeModel *tm_induced_aa(TreeModel *codon_mod);

/** Build index of leaf ids to sequence indices in given alignment.
    @param mod Tree Model containing leaves
    @param msa Multiple Alignment containing sequence data
    @note Leaves not present in the alignment will be ignored.
    @note its not required that there's a leaf for every sequence in the
    alignment.
*/
void tm_build_seq_idx(TreeModel *mod, MSA *msa);

/**  Prune away leaves in tree that don't correspond to sequences in a
    given alignment.
    @param[in,out] mod Tree Model to prune
    @param[in] msa Multiple Alignment; all leaves whose names are not in msa->names will be pruned away
    @param[out] names Will contain names of deleted leaves on return.  Must be pre-allocated
    @warning Root of tree (value of mod->tree) may change.
*/
void tm_prune(TreeModel *mod, MSA *msa, List *names);

/** Extrapolate tree model and prune leaves not represented in
   alignment.
   @param[in] mod Tree Model to prune
   @param[in] extrapolate_tree Super tree of mod->tree to scale
   @param[in] msa Multiple Alignment; all leaves whose names are not in msa->names will be pruned away
   @param[out] pruned_names Will contain names of deleted leaves on return.  Must be pre-allocated
   @result scale factor
   @see tr_scale_by_subtree   */
double tm_extrapolate_and_prune(TreeModel *mod, TreeNode *extrapolate_tree, 
                                MSA *msa, List *pruned_names);

/** Reset TreeModel with new or altered tree.
  @param[in,out] mod Tree model to modify
  @param[in] newtree New tree to model
*/
void tm_reset_tree(TreeModel *mod, TreeNode *newtree);

/** Set branches to be ignored in likelihood calculation and parameter
   estimation.
   @param mod Tree Model containing branches to be ignored
   @param ignore_branches List of Strings indicating nodes in the tree, and the branches leading to those nodes.
*/
void tm_set_ignore_branches(TreeModel *mod, List *ignore_branches);

/** For each parameter, report estimate, variance, and 95% confidence
        interval, printed to given filename, one parameter per line.
    @param mod Tree Model to determine error for
    @param msa Multiple Alignment
    @param params Vector of parameters to test on model
    @param cat Whether using categories
    @param error_fname Filename to save to
    @param appendToFile Whether to append to file or overwrite
*/
void tm_free_alt_subst_mods(TreeModel *tm);

void tm_variance(FILE *outfile, TreeModel *mod, MSA *msa, Vector *params, int cat);


/** Calculate Lower/Upper bounds for each parameter or if not empty, get them from parsing mod->bound_arg
   @pre *lower_bounds and *upper_bounds have not been allocated
   @pre tm_setup_params has already been called
   @param lower_bounds Unallocated vector pointer
   @param upper_bounds Unallocated vector pointer
   @param npar Number of parameters
   @param mod Tree Model to estimate bounds for
   @param allocate_default If TRUE, allocate space for lower_bounds and upper_bounds even if the boundaries are negative and positive infinity, respectively.  If FALSE, set lower_bounds and upper_bounds to NULL to indicate these "default" boundaries.
 */
void tm_new_boundaries(Vector **lower_bounds, 
		       Vector **upper_bounds, 
		       int npar, TreeModel *mod,
		       int allocate_default);


/** Set lower/upper boundaries for each parameter describing a tree model
    @param lower_bounds A vector describing the lower bound for each parameter
    @param upper_bounds A vector describing the upper bound for each parameter
    @param mod A tree model object
 */
void tm_set_boundaries(Vector *lower_bounds, Vector *upper_bounds,
		       TreeModel *mod);

struct likelihood_wrapper_struct {
  MSA *align;
  TreeModel *mod;
  void (*unpackFunction)(Vector *params, TreeModel *mod);
  double (*likelihoodFunction)(Vector *params, void *data);
};


void tm_site_model_set_ml_weights(TreeModel *mod, Vector *params, 
				  double *counts);
void tm_setup_site_model(TreeModel *mod, const char *foreground, int bgc, 
			 int alt_hypothesis, double selNeg, double selPlus,
			 double initBgc, double *initWeights);

List *tm_setup_bgc_model_hmm(TreeModel *mod, const char *foreground, double selNeg,
			     double selPos, double initBgc, double *initWeights);
#endif
