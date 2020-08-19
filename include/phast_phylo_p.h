#ifndef PHYLO_P_H
#define PHYLO_P_H

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <phast_misc.h>
#include <phast_msa.h>
#include <phast_maf.h>
#include <phast_tree_model.h>
#include <phast_sufficient_stats.h>
#include <phast_subst_distrib.h>
#include <phast_prob_vector.h>
#include <phast_prob_matrix.h>
#include "phast_list_of_lists.h"
#include "phast_phylo_p_print.h"
#include "phast_fit_column.h"
#include "phast_fit_feature.h"

/* default values for epsilon; can tolerate larger value with --wig-scores or
   --base-by-base */
#define DEFAULT_EPSILON 1e-10
#define DEFAULT_EPSILON_BASE_BY_BASE 1e-6

typedef enum{SPH, LRT, SCORE, GERP} method_type;

/* maximum size of matrix for which to do explicit convolution of
   joint prior; beyond this size an approximation is used.
   Computational complexity is proportional to square of this number.
   This only comes into play when --features and --subtree are used
   together */
#define MAX_CONVOLVE_SIZE 22500


struct phyloP_struct {
  MSA *msa;
  int prior_only, post_only, quantiles_only;
  int output_wig, output_gff;
  int nsites, fit_model, base_by_base, default_epsilon, refidx, refidx_feat;
  double ci, epsilon;
  char *subtree_name, *chrom;
  List *branch_name;
  GFF_Set *feats;
  method_type method;
  mode_type mode;
  FILE *outfile, *logf;
  TreeModel *mod;
  List *cats_to_do;
  CategoryMap *cm;
  char *help, *mod_fname, *msa_fname;
  ListOfLists *results;
  int no_prune;
};

struct phyloP_struct *phyloP_struct_new(int rphast);
void phyloP(struct phyloP_struct *p);

#endif
