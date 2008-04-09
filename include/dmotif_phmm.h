/* "birth-death" phylo-HMM -- nonconserved state, fully conserved
   state, one state per birth event and one state per death event per
   branch of tree */

#ifndef DMPHMM
#define DMPHMM

#include <tree_model.h>
#include <phylo_hmm.h>
#include <indel_mod.h>
#include <indel_history.h>
#include <gff.h>
#include <pssm.h>
#include <tree_model.h>

typedef enum {NEUT, CONS, DEATH, BIRTH} dmevent_t;

typedef struct{
  PhyloHmm *phmm;               /* phylo-HMM */
  PSSM *m;                      /* underlying motif model */
  int *state_to_branch;         /* mapping from states to branches on
                                   which birth/death events occur */
  int *state_to_event;
  int *state_to_motifpos;
  List **branch_to_states;       /* mapping from branches to states
                                    indicating birth or death events */
  IndelModel **indel_mods;      /* indel models, one per state */
  double rho, mu, nu, phi, zeta;
  int estim_gamma, estim_omega, estim_phi, estim_zeta;
  int k;
} DMotifPhyloHmm;


DMotifPhyloHmm *dm_new(TreeModel *source_mod, PSSM *m, double rho, double mu, 
                       double nu, double phi, double zeta, double alpha_c, 
                       double beta_c, double tau_c, double alpha_n, 
                       double beta_n, double tau_n, int estim_gamma, 
                       int estim_omega, int estim_phi, int estim_zeta);
void dm_set_transitions(DMotifPhyloHmm *dm);
void dm_handle_missing_data(DMotifPhyloHmm *dm, MSA *msa);
void dm_score_predictions(DMotifPhyloHmm *dm, GFF_Set *predictions);
void dm_add_indel_emissions(DMotifPhyloHmm *dm, IndelHistory *ih);
double dm_estimate_transitions(DMotifPhyloHmm *dm, MSA *msa);
void dm_set_backgd_branches(TreeModel *tm, TreeModel *backgd_mod, 
                            List *nodelist);

#endif
