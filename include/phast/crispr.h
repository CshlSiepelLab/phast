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

#endif
