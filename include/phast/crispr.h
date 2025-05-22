/* handling of crispr mutation models for phylogeny reconstruction */
/* reads tab-delimited mutation data, calculates substitution probabilities, custom pruning algorithm for better efficiency */

#ifndef CPR_H
#define CPR_H

#include <stdio.h>
#include <phast/tree_model.h>
#include <phast/vector.h>

typedef struct {
  int nsites;
  int ncells;
  int nstates;
  List *sitenames;
  List *cellnames;
  List *cellmuts;
  Vector *eqfreqs;
} CrisprMutTable;

CrisprMutTable *cpr_new_table();

CrisprMutTable *cpr_read_table(FILE *F);

void cpr_free_table(CrisprMutTable *M);

void cpr_print_table(CrisprMutTable *M, FILE *F);

int cpr_get_mut(CrisprMutTable *M, int cell, int site);

void cpr_set_mut(CrisprMutTable *M, int cell, int site, int val);

void cpr_renumber_states(CrisprMutTable *M);

double cpr_compute_log_likelihood(TreeModel *mod, CrisprMutTable *M,
                                  Vector *branchgrad);

Matrix *cpr_compute_dist(CrisprMutTable *M);

double cpr_compute_pw_dist(CrisprMutTable *M, int i, int j);

void cpr_set_subst_matrices(TreeModel *mod, List *Pt, Vector *eqfreqs);

void cpr_set_branch_matrix(MarkovMatrix *P, double t, Vector *eqfreqs);

void cpr_branch_grad(Matrix *grad, double t, Vector *eqfreqs);

Vector *cpr_estim_mut_rates(CrisprMutTable *M, unsigned int ignore_silent);

void cpr_build_seq_idx(TreeModel *mod, CrisprMutTable *M);

#endif
