/* 
 * migration.h
 *
 * Header file for migration module.  Generalizes tissue migration of
 * BEAM model to any kind of change in cell state.
 */


#ifndef MIGRATION_H
#define MIGRATION_H

#include <stdio.h>
#include <phast/lists.h>
#include <phast/crispr.h>
#include <phast/hashtable.h>

typedef struct {
  int ncells;
  int nstates;
  int nparams; /* number of free parameters in migration model */
  List *cellnames;
  List *states; /* parallel list of corresponding state indices */
  Hashtable *statehash; /* maps state label to state index */
  List *statenames;     /* maps state index to state label */
  Vector *gtr_params; /* GTR parameters */
  Vector *deriv_gtr;  /* gradient w.r.t. GTR parameters */
  List *Pt;       /* list of migration probability matrices for each branch */
  Vector *backgd_freqs; /* equilibrium frequencies of states */
  MarkovMatrix *rate_matrix; /* rate matrix for migration model */
  List **rate_matrix_param_row; /* for each free param, list of row
                                   indices it affects */
  List **rate_matrix_param_col; /* for each free param, list of col
                                   indices it affects */
} MigTable;

MigTable *mig_new();

void mig_free(MigTable *M);

MigTable *mig_read_table(FILE *F);

void mig_update_states(MigTable *M);

void mig_check_table(MigTable *mg, CrisprMutTable *mm);

double mig_compute_log_likelihood(TreeModel *mod, MigTable *mg,
                                  CrisprMutModel *cprmod, Vector *branchgrad);

void mig_set_REV_matrix(MigTable *mg, Vector *params, int start_idx);

void mig_grad_REV_dr(MigTable *mg, List *dP_dr_lst, double t);

void mig_grad_REV_dt(MigTable *mg, Matrix *grad, double t);

void mig_update_subst_matrices(TreeModel *mod, MigTable *mg);

#endif /* MIGRATION_H */
