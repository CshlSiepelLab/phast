/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#ifndef MOTIF_H
#define MOTIF_H

#include "msa.h"
#include "vector.h"
#include "tree_model.h"

#define MTF_EPSILON 0.001
#define MTF_EM_CONVERGENCE_THRESHOLD 0.1
#define HGTRACKS_URL "http://hgwdev-acs.cse.ucsc.edu/cgi-bin/hgTracks?db=hg16"

/* structure for a motif (single- or multi-sequence) */
typedef struct {
  int motif_size;               /* width of motif */
  int multiseq;                 /* whether or not multi-sequence */
  char *alphabet;               /* alphabet for the motif */
  int alph_size;                /* size of alphabet */
  Vector **freqs;           /* array of position-specific base
                                   frequencies */
  TreeModel **ph_mods;          /* array of position-specific
                                   phylogenetic models (NULL if
                                   multiseq == 0). Element 0 represents
                                   background */
  void *training_data;          /* a pointer to a PooledMSA (multiseq == 1) 
                                   or SeqSet (multiseq == 0) */
  int refseq;                   /* reference sequence (-1 if multiseq == 0) */
  double *postprob;             /* posterior probability that each
                                   sample in the training data has the
                                   motif */
  int *bestposition;            /* most likely starting position in
                                   each sample */
  double *samplescore;          /* score of best motif in each sample  */
  double score;                 /* maximized score of whole motif wrt
                                   the training data */
  double *has_motif;            /* used in discriminative training */
  msa_coord_map **coord_maps;   /* used for translation to coord
                                   system of reference sequence */
} Motif;

/* structure for set of individual sequences */
typedef struct {
  MSA *set;                     /* under the hood, just use an MSA
                                   object */
  int *lens;                    /* lengths of sequences (array of size
                                   set->nseqs) */
} SeqSet;

Motif* mtf_new(int motif_size, int multiseq, Vector **freqs, 
               void *training_data, TreeModel *backgd_phmod, 
               double scale_factor);
void mtf_free(Motif *m);
List* mtf_find(void *data, int multiseq, int motif_size, int nmotifs, 
               TreeNode *tree, void *backgd, double *has_motif, double prior, 
               int nrestarts, List *init_list, int sample_parms, 
               int npseudocounts);
double mtf_compute_conditional(Vector *params, void *data);
double mtf_em(void *models, void *data, int nsamples, 
              int *sample_lens, int width, double motif_prior,
              void (*compute_emissions)(double**, void**, int, void*, 
                                        int, int), 
              void (*estimate_state_models)(void**, int, void*, 
                                            double**, int),
              int (*get_observation_index)(void*, int, int),
              double *postprob, int *bestposition);
void mtf_estim_backgd_mn(SeqSet *s, Vector *model);
void mtf_draw_multinomial(Vector *v, double *alpha);
void mtf_get_common_ntuples(SeqSet *s, List *tuples, int tuple_size, 
                            int number);
void mtf_sample_ntuples(SeqSet *s, List *tuples, int tuple_size, int number);
void mtf_init_from_consensus(String *consensus, Vector **mods, 
                             int *inv_alph, int npseudocounts, 
                             int probabilistic, int target_size);
void mtf_winnow_starts(void *data, List *origseqs, int ntochoose, 
                       List *bestseqs, int multiseq, int motif_size, 
                       TreeNode *tree, void *backgd, double *has_motif);
void mtf_get_consensus(Motif *m, char *consensus);
void mtf_print(FILE *F, Motif *m);
void mtf_print_html(FILE *F, Motif *m);
void mtf_print_summary_html(FILE *F, List *motifs, String *prefix);
SeqSet *mtf_get_seqset(List *msas, int refseq, int min_allowable_size);
void mtf_predict(Motif *m, void *data, int *bestposition, double *score, 
                 double *has_motif);
void mtf_add_features(Motif *m, GFF_Set *gff);

#endif
