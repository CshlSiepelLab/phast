/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: phylo_fit.h,v 1.17 2008-11-12 02:07:59 acs Exp $ */

#ifndef PHYLO_FIT_H
#define PHYLO_FIT_H


/* default minimum number of informative sites (see -I) */
#define DEFAULT_NSITES_THRESHOLD 50

/* default starting alpha for dgamma */
#define DEFAULT_ALPHA 1

#include <tree_model.h>
#include <msa.h>
#include <category_map.h>
#include <stringsplus.h>
#include <lists.h>
#include <gff.h>


struct phyloFit_struct {
  MSA *msa;
  char *output_fname_root, 
    *log_fname, *reverse_group_tag,
    *root_seqname, *subtree_name, *error_fname,
    *see_for_help, *parsimony_cost_fname,
    *msa_fname;  //note: msa_fname will not be read by run_phyloFit.  
                 //The msa should already be loaded.  msa_fname can be
                 //optionally set for producing informative messages
  int subst_mod, quiet, nratecats, use_em, 
    window_size, window_shift, use_conditionals,
    precision, likelihood_only, do_bases, 
    do_expected_nsubst, do_expected_nsubst_tot, 
    random_init, estimate_backgd, estimate_scale_only,
    do_column_probs, nonoverlapping, gaps_as_bases,
    no_freqs, no_rates, assume_clock, 
    init_parsimony, parsimony_only, no_branchlens,
    label_categories, symfreq;
  unsigned int nsites_threshold;
  TreeNode *tree;
  CategoryMap *cm;
  String *nooptstr;
  List *cats_to_do_str,  *window_coords, 
    *ignore_branches, *alt_mod_str,
    *bound_arg, *rate_consts;
  double alpha;
  GFF_Set *gff;
  TreeModel *input_mod;

  //results go in these if not-null
  List *estimated_models;
  List *model_labels;
};


struct phyloFit_struct* phyloFit_struct_new();
int run_phyloFit(struct phyloFit_struct *pf);

#endif
