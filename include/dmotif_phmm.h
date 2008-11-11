/* "birth-death" phylo-HMM -- nonconserved state, fully conserved
   state, one state per birth event and one state per death event per
   branch of tree */

#ifndef DMPHMM
#define DMPHMM

#include <tree_model.h>
#include <phylo_hmm.h>
#include <dmotif_indel_mod.h>
#include <indel_history.h>
#include <gff.h>
#include <pssm.h>
#include <tree_model.h>
#include <category_map.h>
#include <hashtable.h>
#include <multi_msa.h>
#include <sufficient_stats.h>

/* typedef enum dmevent_t; /\* Defined in dmotif_indel_mod.h -- declare as an */
/* 			   incomplete type here to avoid errors with  */
/* 			   reciprocal dependencies. *\/ */

typedef struct {
  PhyloHmm *phmm;               /* phylo-HMM */
  PSSM *m;                      /* underlying motif model */
  int *state_to_branch;         /* mapping from states to branches on
                                   which birth/death events occur */
  int *state_to_event;
  int *state_to_motifpos;
  List **branch_to_states;       /* mapping from branches to states
                                    indicating birth or death events */
  DMotifIndelModel **indel_mods; /* indel models, one per state */
  double rho, mu, nu, phi, zeta;
  int estim_gamma, estim_omega, estim_phi, estim_zeta;
  int k;
} DMotifPhyloHmm;

/* Convenience structure for reading in PooledMSA, indel histories and seqnames
   for a list of msa's. Used by dmsample. */
typedef struct {
  PooledMSA *pmsa;
  IndelHistory **ih;
  List *seqnames;
  int max_seqlen;
} DMotifPmsaStruct;


DMotifPhyloHmm *dm_new(TreeModel *source_mod, PSSM *m, double rho, double mu, 
                       double nu, double phi, double zeta, double alpha_c, 
                       double beta_c, double tau_c, double epsilon_c, 
		       double alpha_n, double beta_n, double tau_n,
		       double epsilon_n, int estim_gamma, int estim_omega, 
		       int estim_phi, int estim_zeta);
void dm_set_transitions(DMotifPhyloHmm *dm);
void dm_handle_missing_data(DMotifPhyloHmm *dm, MSA *msa);
void dm_score_predictions(DMotifPhyloHmm *dm, GFF_Set *predictions);
void dm_score_feat(DMotifPhyloHmm *dm, GFF_Feature *f, int cbstate);
void dm_add_indel_emissions(DMotifPhyloHmm *dm, IndelHistory *ih);
double dm_estimate_transitions(DMotifPhyloHmm *dm, MSA *msa);
void dm_set_backgd_branches(TreeModel *tm, TreeModel *backgd_mod, 
                            List *nodelist);
CategoryMap* dm_create_catmap(DMotifPhyloHmm *dm, TreeModel *source_mod,
			      char *feature_prefix);
GFF_Set* dm_phmm_predict_viterbi(DMotifPhyloHmm *dm, char *seqname, 
				 char *grouptag, char *idpref, List *frame);
GFF_Set *dm_labeling_as_gff(CategoryMap *cm, int *path, int length, int w, 
			    int *path_to_cat, int *state_to_motifpos, 
			    int *reverse_compl, char *seqname, char *source, 
			    List *frame_cats, char *grouproot, char *idpref);
void dm_print_motif_scores(DMotifPhyloHmm *dm);
void dms_sample_paths(DMotifPhyloHmm *dm, PooledMSA *blocks, 
		      double **tuple_scores, IndelHistory **ih,
		      List *seqnames, int max_seqlen, int bsamples, 
		      int nsamples, int sample_interval, 
		      Hashtable *path_counts, int **priors, FILE *log, 
		      GFF_Set *reference, int ref_as_prior, int force_priors);
void dms_read_priors(int **priors, FILE *prior_f);
GFF_Feature* dms_motif_as_gff_feat(DMotifPhyloHmm *dm, PooledMSA *blocks, 
				   List *seqnames, char *key, int *counts, 
				   int nsamples, int sample_interval);
void dms_compare_gffs(GFF_Set *reference, GFF_Set *query, int *stats, 
		      int offset, /*< Can be used to adjust the coordinate base
				    of the query set to match the target set
				    (e.g., from 0-based to 1-based) */
		      GFF_Set *matches, GFF_Set *mismatches,
		      GFF_Set *unique_to_query, GFF_Set *unique_to_target
		      );
int* dms_composite_path(DMotifPhyloHmm *dm, int *path1, int *path2, 
			int seqlen, int offset, int force_priors);
int* dm_gff_to_path(DMotifPhyloHmm *dm, GFF_Set *gff, int seqlen,
		    int offset);
int** dm_gff_to_paths(DMotifPhyloHmm *dm, GFF_Set *gff, Multi_MSA *blocks,
		      int offset);
void dms_combine_gffs(GFF_Set *target_gff, GFF_Set *query_gff);
Hashtable* dms_read_hash(FILE *hash_f, int nstates, int* nsamples);
void dms_write_hash(Hashtable *path_counts, FILE *hash_f, int nstates, 
		    int nsamples);
void dms_combine_hashes(Hashtable *target, Hashtable *query, int nstates);
void dms_compute_emissions(PhyloHmm *phmm, MSA *pmsa, int quiet);
void dms_lookup_emissions(DMotifPhyloHmm *dm, double **tuple_scores, 
			  PooledMSA *blocks, int seqnum, int seqlen, 
			  IndelHistory *ih);
void dms_count_transitions(DMotifPhyloHmm *dm, int *path, int **trans, 
			   int **priors, int seqlen, int *ref_path,
			   int force_priors);
void dms_count_motifs(DMotifPhyloHmm *dm, int *path, int seqlen,
		      Hashtable *path_counts, int seqnum);
void dms_path_log(DMotifPhyloHmm *dm, int *path, int seqlen, char *seqname,
		  GFF_Set *motifs);
void dms_write_log(FILE *log, DMotifPhyloHmm *dm, int **trans, int sample, 
		   double llh, GFF_Set *query_gff, GFF_Set *reference, 
		   int nwins);
DMotifPmsaStruct *dms_read_alignments(FILE *F, int do_ih);
double dm_compute_log_likelihood(TreeModel *mod, MSA *msa, double *col_scores,
				 int cat);
void dm_free_subst_matrices(TreeModel *tm);

#endif
