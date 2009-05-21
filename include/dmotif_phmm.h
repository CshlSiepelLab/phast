/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2009 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

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
#include <pthr.h>

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
  int *state_to_motif;           /* Mapping of state to motif/branch comb. */
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
  List **zeroed_states;
  int max_seqlen;
} DMotifPmsaStruct;

/* Convenience structure for emissions multithreading */
typedef struct {
  TreeModel *mod;
  MSA *msa;
  double *emissions_vec;
  int cat;
  int state;
  int nstates;
  int quiet;
} DMemissionsThreadData;

/* Convenience structure for sampling multithreading */
typedef struct {
  DMotifPhyloHmm *dm;
  PooledMSA *blocks;
  IndelHistory *ih;
  List *zeroed_states;
  double **tuple_scores;
  int **thread_path;
  double ***thread_emissions;
  double ***thread_forward;
  double *thread_llh;
  int ***thread_trans;
  Hashtable **thread_counts;
  int do_sample;
  int seqnum;
  char *seqname;
  FILE *log;
  GFF_Set *query_gff;
  int do_reference;
  int nthreads;
  ThreadPool *p;
  int sample;
} DMsamplingThreadData;

/* Structure to specify states and positions to zero out when conditioning on
   presence of substitutions for gain/loss states and presence of a site in a
   specific specie(s) */
typedef struct {
  int state;
  int do_row;
  List *starts;
  List *lengths;
} DMzeroedState;

DMotifPhyloHmm *dm_new(TreeModel *source_mod, PSSM *m, double rho, double mu, 
                       double nu, double phi, double zeta, double alpha_c, 
                       double beta_c, double tau_c, double epsilon_c, 
		       double alpha_n, double beta_n, double tau_n,
		       double epsilon_n, int estim_gamma, int estim_omega, 
		       int estim_phi, int estim_zeta);
void dm_free(DMotifPhyloHmm *dm);
void dm_set_transitions(DMotifPhyloHmm *dm);
void dm_handle_missing_data(DMotifPhyloHmm *dm, MSA *msa);
void dm_score_predictions(DMotifPhyloHmm *dm, GFF_Set *predictions);
void dm_score_feat(DMotifPhyloHmm *dm, GFF_Feature *f, int cbstate);
void dm_add_indel_emissions(DMotifPhyloHmm *dm, double **emissions,
			    IndelHistory *ih);
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
List* dms_sample_paths(DMotifPhyloHmm *dm, PooledMSA *blocks, 
		       double **tuple_scores, IndelHistory **ih,
		       List *seqnames, int max_seqlen, int bsamples, 
		       int nsamples, int sample_interval, 
		       int **priors, FILE *log, 
		       GFF_Set *reference, int ref_as_prior, 
		       int force_priors,
		       int quiet, char *cache_fname, int cache_int);
/* Multithreaded version of dms_sample_paths */
List* dms_sample_paths_pthr(DMotifPhyloHmm *dm, PooledMSA *blocks,
			    double **tuple_scores, IndelHistory **ih,
			    List *seqnames, int max_seqlen, int bsamples,
			    int nsamples, int sample_interval, int **priors,
			    FILE *log, GFF_Set *reference, int ref_as_prior,
			    int force_priors, int quiet, char *cache_fname,
			    int cache_int, ThreadPool *pool, int nthreads,
			    List **zeroed_states);
void dms_sample_path(DMotifPhyloHmm *dm, PooledMSA *blocks, IndelHistory *ih,
		     double **tuple_scores, double *thread_llh,
		     int **thread_path,
		     double ***thread_emissions, double ***thread_forward,
		     int ***thread_trans, Hashtable **thread_counts, 
		     int do_sample, int seqnum, char *seqname, FILE *log, 
		     GFF_Set *query_gff, int do_reference, int nthreads,
		     ThreadPool *p, int sample, List *zeroed_states);
void dms_launch_sample_thread(void *data);
void dms_read_priors(int **priors, FILE *prior_f);
GFF_Feature* dms_motif_as_gff_feat(DMotifPhyloHmm *dm, PooledMSA *blocks, 
				   List *seqnames, char *key, int *counts, 
				   int nsamples, int sample_interval,
				   int refidx);
void dms_compare_gffs(GFF_Set *reference, GFF_Set *query, int *stats, 
		      int offset, /*< Can be used to adjust the coordinate base
				    of the query set to match the target set
				    (e.g., from 0-based to 1-based) */
		      GFF_Set *matches, GFF_Set *mismatches,
		      GFF_Set *unique_to_query, GFF_Set *unique_to_target);
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
Hashtable* dms_cache_hash(Hashtable *path_counts, char *hash_fname,
			  int nstates, int nsamples, int reinit_val);
void dms_combine_hashes(Hashtable *target, Hashtable *query, int nstates);
void dms_compute_emissions(PhyloHmm *phmm, MSA *pmsa, int quiet, 
			   ThreadPool *pool, int nthreads);
void dms_lookup_emissions(DMotifPhyloHmm *dm, double **tuple_scores,
			  double **emissions, PooledMSA *blocks, int seqnum,
			  int seqlen, IndelHistory *ih);
void dms_count_transitions(DMotifPhyloHmm *dm, int *path, int **trans, 
			   int seqlen, int *ref_path, int force_priors);
void dms_count_motifs(DMotifPhyloHmm *dm, int *path, int seqlen,
		      Hashtable *path_counts, int seqnum);
void dms_path_log(DMotifPhyloHmm *dm, int *path, int seqlen, char *seqname,
		  GFF_Set *motifs);
void dms_write_log(FILE *log, DMotifPhyloHmm *dm, int **trans, int sample, 
		   double llh, GFF_Set *query_gff, GFF_Set *reference, 
		   int nwins);
DMotifPmsaStruct *dms_read_alignments(FILE *F, int do_ih, int quiet,
				      int revcomp, int do_zeroed,
				      FILE *cond_spec_f);
double dm_compute_log_likelihood(TreeModel *mod, MSA *msa, double *col_scores,
				 int cat);
void dm_free_subst_matrices(TreeModel *tm);
MSA *dm_indel_mask(DMotifPhyloHmm *dm, MSA *msa, IndelHistory *ih,
		   int *path);
Hashtable* dms_uncache(List *cache_files, int init_size, int nstates,
		       int* nsamples, int quiet);
void dms_do_emissions_row(void *data);
void dm_set_subst_mods(DMotifPhyloHmm *dm);
List* dms_read_tmp_from_file(FILE *tmp_lst_f);
void dms_dump_sample_data(int sample, int thread_id, char *seqname, 
			  int seqlen, int *path, int **trans, 
			  Hashtable *path_counts, FILE *out, int dim);
void dms_free_dmpmsa_struct(DMotifPmsaStruct *dmpmsa);
void dms_map_gff_coords(PooledMSA *blocks, int seqidx, GFF_Feature *f,
			int from_seq, int to_seq);

/* Functions related to conditioning on presence of sites and substitutions */
void dms_zero_states(MSA *msa, double **emissions, List *zeroed_states);
void dms_zero_emissions_row(double *row, int seqlen);
DMzeroedState **dms_condition_on_species(DMotifPhyloHmm *dm, TreeNode *tree,
					 int cond_on /* Node ID of (leaf) 
							sequence believed to 
							contain a site */
			       );
void dms_condition_on_subs(DMotifPhyloHmm *dm, TreeModel *mod, MSA *msa, 
			   DMzeroedState **zeroed_states, int nosubs);
DMzeroedState *dms_new_zeroed_state(int state, int do_row);
void dms_free_zeroed_state(DMzeroedState *z);
void dms_write_zeroed_states(FILE *f, List *zeroed_states);
List *dms_read_zeroed_states(FILE *F);
List *dms_zeroed_states_array_to_list(DMzeroedState **zeroed_states,
				      int nstates, int init_size);
void dms_print_zeroed_states(FILE *f, DMotifPhyloHmm *dm, TreeNode *tree,
			     List *zeroed_states);
List *dms_reverse_zeroed_states(List *zeroed_states, MSA *msa);

/* Create an array of Lists that describe nodes to zero emissions for at
   each motif window. Nosubs toggles between choosing branches to zero
   based on whether parent and child character sets are identical (nosubs ==
   0) or whether scenarios exist that allow substitutions on a branch (nosubs
   == 1) */
List **dms_label_subst_nodes(MSA *msa, TreeModel *mod, PSSM *m, int nosubs);

/* Create a list that represents the intersection of two int lists or, if there
   is no overlap between the two lists, the union of both lists. Returns the
   size of the intersection or 0 if dest represents the union. */
int dms_intersect_or_union(List *dest, List *left, List *right);
void dms_compact_zeroed_states(List **zeroed_states, int nblocks);
int dms_compare_lists_identity(List *parent, List *child);
int dms_compare_lists_nosubs(List *parent, List *child);
void dms_zeroed_as_bed(FILE *f, DMotifPhyloHmm *dm, List *zeroed_states,
		       char *fname);
void dms_condition_transitions(DMotifPhyloHmm *dm, List *zeroed_states);
MSA *dm_generate_msa(int ncolumns, 
                     DMotifPhyloHmm *dm,
                     TreeModel **classmods, 
                     int *labels /* if non-NULL, will be used to
                                    record state (model) responsible
                                    for generating each site; pass
                                    NULL if hmm is NULL */
                     );
void dm_sample_char_col(char **seqs, TreeModel *mod, char *newchar, 
			int class, int col, int keep_ancestral);

#endif
