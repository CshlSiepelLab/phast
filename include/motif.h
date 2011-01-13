/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** 
   @file motif.h
   Finding motifs in sequences
   @ingroup motif

*/

#ifndef MOTIF_H
#define MOTIF_H

#include "msa.h"
#include "vector.h"
#include "tree_model.h"

/** Epsilon value for motifs */
#define MTF_EPSILON 0.001
/** Threshold met for convergence when finding a motif */
#define MTF_EM_CONVERGENCE_THRESHOLD 0.1
/** URL to use when creating links on output HTML webpages */
#define HGTRACKS_URL "http://hgwdev-acs.cse.ucsc.edu/cgi-bin/hgTracks?db=hg16"

/** Structure for a motif (single- or multi-sequence) */
typedef struct {
  int motif_size;               /**< Width of motif */
  int multiseq;                 /**< Whether or not multi-sequence */
  char *alphabet;               /**< Alphabet for the motif */
  int alph_size;                /**< Size of alphabet */
  Vector **freqs;           /**< array of position-specific base
                                   frequencies */
  TreeModel **ph_mods;          /**< Array of position-specific
                                   phylogenetic models (NULL if
                                   multiseq == 0). Element 0 represents
                                   background */
  void *training_data;          /**< A pointer to a PooledMSA (multiseq == 1) 
                                   or SeqSet (multiseq == 0) */
  int refseq;                   /**< Reference sequence (-1 if multiseq == 0) */
  double *postprob;             /**< Posterior probability that each
                                   sample in the training data has the
                                   motif */
  int *bestposition;            /**< Most likely starting position in
                                   each sample */
  double *samplescore;          /**< Score of best motif in each sample  */
  double score;                 /**< Maximized score of whole motif wrt
                                   the training data */
  double *has_motif;            /**< Used in discriminative training */
  msa_coord_map **coord_maps;   /**< Used for translation to coord
                                   system of reference sequence */
} Motif;

/** Structure for set of individual sequences */
typedef struct {
  MSA *set;                     /**< Under the hood, just use an MSA
                                   object */
  int *lens;                    /**< Lengths of sequences (array of size
                                   set->nseqs) */
} SeqSet;

/** \name Motif Allocation/Cleanup functions 
\{ */

/** Create a new Motif object from parameters.
    @param motif_size Size of the new motif
    @param multiseq If multiple sequences then = 1
    @param freqs Array of possible position-specific base frequencies
    @param training_data A pointer to a PooledMSA (multiseq == 1) or SeqSet (multiseq == 0)
    @param backgd_phmod (Multi-Seq Only) Array of position-specific phylogenetic models (NULL if multiseq == 0). Element 0 represents background.
    @param scale_factor (Multi-Seq Only) Scale branch lengths by this amount
    @result Newly created Motif object populated with data
*/
Motif* mtf_new(int motif_size, int multiseq, Vector **freqs, 
               void *training_data, TreeModel *backgd_phmod, 
               double scale_factor);

/** Free a Motif object.
   @param m Motif object to free
*/
void mtf_free(Motif *m);

/** \} */

/** Find motifs in one or more sequences.
   Uses either using EM or discriminative training.  
   @param data Sequence data to look for motif in.  If 'multiseq' == 1 then 'data' must be a PooledMSA object; otherwise, it must be a SeqSet object.   
   @param nmotif Upper limit on number of motifs to find.
   @param motif_size Size of Motif trying to find
   @param backgd Background model; a TreeModel (if multiseq == 1), otherwise a Vector (only used for initialization)
   @param has_modif If non-NULL, then motifs are learned by
   discriminative training, with 'has_motif' interpreted as an array
   of values indicating whether each training example should be
   considered a positive (1) or negative (0) example, or somewhere in
   between.  
   @param prior Initial value for the prior probability that a motif instance appears in each sequence (used
   with EM only).  
   @param init_list List of initial tuples, mtf_sample_ntuples or mtf_get_common_ntuples works well to generate this
   @param sample_parms If == 1 Sample
                     parameters from a Dirichlet distribution defined by the
                     pseudocounts. If == 0 Do a deterministic
                     initialization based on a consensus sequence
   @param npseudocounts Number of Pseudo counts for consensus bases
   @result List of Motif objects. 
*/
List* mtf_find(void *data, int multiseq, int motif_size, int nmotifs, 
               TreeNode *tree, void *backgd, double *has_motif, double prior, 
               int nrestarts, List *init_list, int sample_parms, 
               int npseudocounts);

/** This is the function that is optimized in discriminative training;
   see Segal et al., RECOMB '02 */
double mtf_compute_conditional(Vector *params, void *data);

/** Find a single motif by EM, given a pre-initialized set of models. 
   Functions must be provided for computing "emission" probabilities
   under all models, for updating model parameters given posterior
   expected counts, and for indexing observations.  It is assumed that
   the motif appears at most once in each sequence of observations.
   @param models  Array of position-specific phylogenetic models i.e. motif->ph_mods
   @param data Sequence data as either PooledMSA (multiple MSAs) or SeqSet (single MSA)
   @param nsamples Number of data samples
   @param sample_lens Array of size nsamples, each element indicates the length of a sample in data
   @param width Width of the motif, models are also assumed to be width+1
   @param motif_prior Prior probability that each sequence has a motif.  
   @param compute_emissions Function to compute emissions
   @param estimate_state_models Function to estimate state models
   @param get_observation_index Function to get observation index
   @param postprob (Optional) Array of size nsamples to be populated with Posterior probabilities that a motif appears
   @param bestposition (Optional) Array of size nsamples to be populated with starting position of the best instance of the motif
   @result Maximized log likelihood.  
   @note This function can be used with phylogenetic models or ordinary multinomial models.  
   @note The first model is assumed to represent the background distribution and its parameter
   @warning The models that are passed in are updated, and at convergence, represent the (apparent) m.l.e. rs will not be updated.
*/
double mtf_em(void *models, void *data, int nsamples, 
              int *sample_lens, int width, double motif_prior,
              void (*compute_emissions)(double**, void**, int, void*, 
                                        int, int), 
              void (*estimate_state_models)(void**, int, void*, 
                                            double**, int),
              int (*get_observation_index)(void*, int, int),
              double *postprob, int *bestposition);

/** Estimate a (multinomial) background model from a set of sequences.
   @param[in] s Set of sequences to estimate background model from
   @param[out] model Estimated derive a consensus sequence from a motif (majority each position)
   Assumes DNA rather than amino-acid alphabet. background model
*/
void mtf_estim_backgd_mn(SeqSet *s, Vector *model);

/** Draw the parameters of a multinomial distribution (from a Dirichlet distrib).
   @param[out] v Scale of distribution 
   @param[in] alpha Array of distribution parameters of. Array alpha is same size as v
*/
void mtf_draw_multinomial(Vector *v, double *alpha);

/** \name Create list of tuples 
  \{  */

/** Scan a set of DNA sequences for the most prevalent n-tuples of bases.
    @param[in] s Set of DNA sequences to scan through
    @param[out] tuples List of most common tuples
    @param[in] tuple_size Size of tuple(s) to find
    @param[in] number Number of tuples to put into parameter 'tuples'
*/
void mtf_get_common_ntuples(SeqSet *s, List *tuples, int tuple_size, 
                            int number);

/** Randomly sample n-tuples from Sequence Set.
   @param[in] s Set of sequences to scan through
   @param[out] tuples List of the randomly sampled tuples
   @param[in] tuple_size Size of tuple(s) to find
   @param[in] number Number of tuples to put into parameter 'tuples'
 */
void mtf_sample_ntuples(SeqSet *s, List *tuples, int tuple_size, int number);
/** \} */

/** Perform a "soft" initialization of a series of multinomial models
   from a consensus sequence. 
   @param consensus Consensus sequence to perform initialization from
   @param mods Models
   @param inv_alpha Inverted Alphabet
   @param npseudocounts Number of pseudo counts
   @param probalistic Whether treated as probabilistic (draw parameters from Dirichlet) or deterministic (treat pseudocounts as counts)
   @param target_size Size sequences generated for initialization
   @note If target size is larger than consensus size, flanking models will be added with uniform
   distributions (or draws from uniform Dirichlet) */
void mtf_init_from_consensus(String *consensus, Vector **mods, 
                             int *inv_alph, int npseudocounts, 
                             int probabilistic, int target_size);

/**  Subset of consensus sequences representing starting models.
   @param data Models PooledMSA (multiple MSAs) or SeqSet (single MSA)
   @param[in] origseqs List of strings, containing the original sequences
   @param[in] ntochoose Number of sequences to choose
   @param[out] bestseqs List of original sequences that survive score thresholding
   @param[in] multiseq Whether there are multiple MSAs (multiseq==1) or a single sequence set (multiseq==0)
   @param[in] motif_size Size of motif to look for
   @param tree NOT USED
   @param[in] backgd Vector of background frequencies (If multiseq==1 is instead a pointer to a TreeModel with background frequencies residing in TreeModel->backgd_freqs)
   @param has_motif Array indicating whether each sequences has a motif or not
   @result Subset that looks promising based on unmaximized scores
   @note This is similar to the heuristic used by MEME
*/
void mtf_winnow_starts(void *data, List *origseqs, int ntochoose, 
                       List *bestseqs, int multiseq, int motif_size, 
                       TreeNode *tree, void *backgd, double *has_motif);

/** Derive a consensus sequence from a motif (majority each position).
   @pre Sequences must be DNA alphabet
   @param[in] m Motif to derive consensus sequence from
   @param[out] consensus Consensus Sequence
*/
void mtf_get_consensus(Motif *m, char *consensus);

/** \name Motif Write to file functions 
\{ */

/** Save Motif to file. 
    @param File descriptor for file to save to
    @param m Motif to save to file
*/
void mtf_print(FILE *F, Motif *m);

/** Save Motif to file as an HTML document.
    @param File descriptor for file to save to
    @param m Motif to save to file
*/
void mtf_print_html(FILE *F, Motif *m);
/** Save a summary of motifs to a file.
   @param F File descriptor for file to save to
   @param motifs List of motifs to print summary about
   @param prefix First part of filename for motif specific documents i.e. 'phastm' makes links to phastm.1.html, phastm.2.html, etc.
 */
void mtf_print_summary_html(FILE *F, List *motifs, String *prefix);

/** \} */

/** Create a set of individual (gapless) sequences from a set of
   multiple alignments. 
   @param msas List of Alignments
   @param refseq The new set can either include only the
   reference sequence from each alignment (set 'refseq' to desired
   one-based index) or all sequences (set 'refseq' to -1).  
   @param min_allowable_size Sequences smaller than this length are not included in result
   @result New SeqSet object.
*/
SeqSet *mtf_get_seqset(List *msas, int refseq, int min_allowable_size);

/** Predict the best motif in each sample of a data set.  
   @param m If m->multiseq, data should be a PooledMSA, otherwise it should be a
   SeqSet.  
   @param bestposition The starting position
   of the best motifs for each sample i is stored in bestposition[i]
   @param score The score of the best motifs for each sample i is stored in score[i].
   @param has_motif If has_motif != NULL, then predictions will be made only
   for samples i such that has_motif[i] >= 0.5.  
*/
void mtf_predict(Motif *m, void *data, int *bestposition, double *score, 
                 double *has_motif);

/** Add a feature to a gff for each predicted instance of a motif in
   the training set.
   @param m Motif containing predicted motifs
   @param gff Feature Set to add new feature(s) to
 */
void mtf_add_features(Motif *m, GFF_Set *gff);

#endif
