#ifndef PHMM_H
#define PHMM_H

#include <category_map.h>
#include <hmm.h>
#include <tree_model.h>
#include <lists.h>

typedef struct {
  CategoryMap *cm;              /**< category map */
  HMM *hmm;                     /**< final HMM, after reflection and
                                   allowance for rate categories  */
  HMM *functional_hmm, *autocorr_hmm;
                                /**< original HMMs used to create cross
                                   product (NULL if no rate
                                   categories) */
  TreeModel **mods;             /**< array of tree models, after
                                   allowance for rate categories  */
  int nmods;                    /**< number of tree models (length of
                                   array mods) */
  int *state_to_mod;            /**< mapping of HMM state number to tree
                                   model number */
  int *state_to_cat;            /**< mapping of HMM state number to (spooled)
                                   category number */
  int *reverse_compl;           /**< array of length hmm->nstates
                                   with value 1 for each state
                                   that corresponds to the reverse
                                   strand and value 0 otherwise */
  List **cat_to_states;         /**< one to many mapping */
  int *state_to_pattern;        /**< gap pattern associated with each
                                   state, when modeling indels (-1 for
                                   no gap pattern) */
  int nratecats;                /**< number of rate categories (for
                                   cross-product constructions) */
  int reflected;                /**< whether "reflected" for reverse compl */
  double **emissions;           /**< values computed by
                                   phmm_compute_emissions */
  double **forward;             /**< forward scores */
  int alloc_len;                /**< length for which emissions and/or
                                   forward are (or are to be) allocated */
  int *state_pos, *state_neg;   /**< contain tracking data for emissions */
} PhyloHmm;

PhyloHmm *phmm_new(HMM *hmm, TreeModel **tree_models, CategoryMap *cm, 
                   List *pivot_cats, List *indel_cats, int nseqs);
void phmm_reflect_hmm(PhyloHmm *phmm, List *pivot_cats);
void phmm_create_autocorr_hmm(HMM *hmm, double lambda);
void phmm_rates_cross(PhyloHmm *phmm, int nratecats, double lambda, 
                      int expand_cats);
void phmm_update_cross_prod(PhyloHmm *phmm, double lambda);
void phmm_free(PhyloHmm *phmm);
void phmm_compute_emissions(PhyloHmm *phmm, MSA *msa, int quiet);
double phmm_fit_lambda(PhyloHmm *phmm, double *lambda, FILE *logf);
void phmm_update_cross_prod(PhyloHmm *phmm, double lambda);
GFF_Set* phmm_predict_viterbi(PhyloHmm *phmm, char *seqname, List *frame);
GFF_Set* phmm_predict_viterbi_cats(PhyloHmm *phmm, List *cats, char *seqname,
                                     List *frame, char *new_type );
double phmm_lnl(PhyloHmm *phmm);
double phmm_postprobs(PhyloHmm *phmm, double **post_probs);
double* phmm_postprobs_cats(PhyloHmm *phmm, List *cats, double *lnl);
void phmm_score_predictions(PhyloHmm *phmm, GFF_Set *preds, List *score_cats, 
                            List *helper_cats, List *null_cats,
                            int score_only_score_cats);
void phmm_rates_cut(PhyloHmm *phmm, int nrates, int cut_idx, double p, double q);
double phmm_fit_rates_cut(PhyloHmm *phmm, double *p, double *q, FILE *logf);
double phmm_fit_rates_cut_em(PhyloHmm *phmm, double *p, double *q, FILE *logf);

#endif
