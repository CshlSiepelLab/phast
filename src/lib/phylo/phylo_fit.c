/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: trees.c,v 1.25 2008-11-12 02:07:59 acs Exp $ */

/*
  Contains phyloFit function, the main engine behind the phyloFit program.
*/

#include <stdlib.h>
#include <stdio.h>
#include <lists.h>
#include <stringsplus.h>
#include <msa.h>
#include <gff.h>
#include <category_map.h>
#include <getopt.h>
#include <tree_model.h>
#include <fit_em.h>
#include <subst_mods.h>
#include <local_alignment.h>
#include <ctype.h>
#include <tree_likelihoods.h>
#include <numerical_opt.h>
#include <sufficient_stats.h>
#include <maf.h>
#include <phylo_fit.h>
#include <stacks.h>
#include <trees.h>
#include <misc.h>

/* initialize phyloFit options to defaults (slightly different
   for rphast).
 */
struct phyloFit_struct* phyloFit_struct_new(int rphast) {
  struct phyloFit_struct *pf = smalloc(sizeof(struct phyloFit_struct));
  pf->msa = NULL;
  pf->output_fname_root = rphast ? NULL : "phyloFit";
  pf->logf = NULL;
  pf->reverse_group_tag = NULL;
  pf->root_seqname = NULL;
  pf->subtree_name = NULL;
  pf->error_fname = NULL;
  pf->see_for_help = rphast ? "phyloFit" : "'phyloFit -h'";
  pf->parsimony_cost_fname = NULL;
  pf->msa_fname = NULL;
  pf->subst_mod = UNDEF_MOD;
  pf->quiet = FALSE; //probably want to switch to TRUE for rphast after debugging
  pf->nratecats = -1;
  pf->use_em = FALSE;
  pf->window_size = -1;
  pf->window_shift = -1;
  pf->use_conditionals = FALSE;
  pf->precision = OPT_HIGH_PREC;
  pf->likelihood_only = FALSE;
  pf->do_bases = FALSE;
  pf->do_expected_nsubst = FALSE;
  pf->do_expected_nsubst_tot = FALSE;
  pf->do_expected_nsubst_col = FALSE;
  pf->random_init = FALSE;
  pf->estimate_backgd = rphast ? TRUE : FALSE;  //in rphast mode we use no.opt to specify no estimate backgd
  pf->estimate_scale_only = FALSE;
  pf->do_column_probs = FALSE;
  pf->nonoverlapping = FALSE;
  pf->gaps_as_bases = FALSE;
  pf->no_freqs = FALSE;
  pf->init_backgd_from_data = TRUE;
  pf->no_rates = FALSE;
  pf->assume_clock = FALSE;
  pf->init_parsimony = FALSE;
  pf->parsimony_only = FALSE;
  pf->no_branchlens = FALSE;
  pf->label_categories = TRUE;  //if false, assume MSA already has
                                //categories labelled with correct gff
  pf->nsites_threshold = DEFAULT_NSITES_THRESHOLD;
  pf->tree = NULL;
  pf->cm = NULL;
  pf->symfreq = FALSE;
  pf->nooptstr = NULL;
  pf->cats_to_do_str = NULL;
  pf->window_coords=NULL;
  pf->ignore_branches = NULL;
  pf->alt_mod_str = NULL;
  pf->label_str = NULL;
  pf->label_type = NULL;
  pf->bound_arg = NULL;
  pf->rate_consts = NULL;
  pf->alpha = DEFAULT_ALPHA;
  pf->gff = NULL;
  pf->input_mod = NULL;
  pf->use_selection = 0;
  pf->selection = 0.0;
  pf->max_em_its = -1;
  
  pf->results = rphast ? lol_new(2) : NULL;
  return pf;
}


void set_output_fname(String *fname, char *root, int cat, char *suffix) {
  str_cpy_charstr(fname, root);
  if (cat != -1) {
    str_append_charstr(fname, ".");
    str_append_int(fname, cat);
  }
  str_append_charstr(fname, suffix);
}



char **get_ratecat_names(TreeModel *mod, int *len) {
  int i;
  char tempstr2[100];
  char **rv = smalloc(mod->nratecats * sizeof(char*));
  for (i=0; i < mod->nratecats; i++) {
    sprintf(tempstr2, "rate.cat.%i", i);
    rv[i] = smalloc((strlen(tempstr2)+1)*sizeof(char));
    strcpy(rv[i], tempstr2);
  }
  *len = mod->nratecats;
  return rv;
}


char **get_node_names(TreeModel *mod, int include_leaf_branches,
		      int include_root, int add_count_column,
		      int *len) {
  int node, idx=0;
  TreeNode *n;
  char **rv;

  *len = mod->tree->nnodes;
  if (add_count_column) (*len)++;
  if (!include_root) (*len)--;
  if (!include_leaf_branches) (*len) -= (mod->tree->nnodes+1)/2;
  rv = smalloc((*len)*sizeof(char*));
  if (add_count_column) {
    rv[idx] = smalloc(6*sizeof(char));
    strcpy(rv[idx++], "nsite");
  }
  tr_name_ancestors(mod->tree);
  for (node = 0; node < mod->tree->nnodes; node++) {
    n = lst_get_ptr(mod->tree->nodes, node);
    if (n == mod->tree && !include_root) continue;
    if ((n->lchild == NULL || n->rchild == NULL) && !include_leaf_branches)
      continue;
    rv[idx] = smalloc((strlen(n->name)+1)*sizeof(char));
    strcpy(rv[idx++], n->name);
  }
  return rv;
}

char **get_state_names(TreeModel *mod, const char *prefix, int *len) {
  int state, prefixlen;
  char **rv = smalloc(mod->rate_matrix->size*sizeof(char*));
  *len = mod->rate_matrix->size;
  if (prefix == NULL) prefixlen = 0;
  else prefixlen = (int)strlen(prefix);
  for (state=0; state < mod->rate_matrix->size; state++) {
    rv[state] = smalloc((mod->order + 2 + prefixlen)*sizeof(char));
    rv[state][mod->order + 1 + prefixlen] = '\0';
    if (prefix != NULL) strcpy(rv[state], prefix);
    get_tuple_str(rv[state]+prefixlen, state, mod->order + 1, 
		  mod->rate_matrix->states);
  }
  return rv;
}


char **get_site_names(MSA *msa, int *len) {
  char **rv;
  char tempch[100];
  int i;
  *len = (int)msa->length;
  rv = smalloc((*len) * sizeof(char*));
  for (i=0; i < *len; i++) {
    sprintf(tempch, "%i", i+1);
    rv[i] = smalloc((strlen(tempch)+1)*sizeof(char));
    strcpy(rv[i], tempch);
  }
  return rv;
}


char **get_tuple_names(TreeModel *mod, MSA *msa, int cat, int *len) {
  int tup, idx=0;
  char **rv;
  *len = 0;
  for (tup = 0; tup < msa->ss->ntuples; tup++) {
    if ((cat >= 0 && msa->ss->cat_counts[cat][tup] == 0) ||
	msa->ss->counts[tup] == 0) continue;
    (*len)++;
  }
  rv = smalloc((*len)*sizeof(char*));
  for (tup = 0; tup < msa->ss->ntuples; tup++) {
    checkInterruptN(tup, 1000);
    if ((cat >= 0 && msa->ss->cat_counts[cat][tup] == 0) ||
	msa->ss->counts[tup] == 0) continue;
    rv[idx] = smalloc((msa->nseqs+1)*sizeof(char));
    tuple_to_string_pretty(rv[idx++], msa, tup);
  }
  return rv;
}


/* Compute and output statistics based on posterior probabilities,
   including (optionally) the post prob of each tuple of bases at each
   ancestral node at each site (do_bases), the expected total number
   of substs per site (do_expected_nsubst), and the expected number of
   substitutions of each type on each edge across all sites
   (do_expected_nsubst_tot).  A separate file is output for each
   selected option, with an appropriate filename suffix (".postprob",
   ".expsub", and ".exptotsub", respectively).    */
void print_post_prob_stats(TreeModel *mod, MSA *msa, char *output_fname_root, 
                           int do_bases, int do_expected_nsubst, 
                           int do_expected_nsubst_tot, 
			   int do_expected_nsubst_col,
			   int do_every_site,
			   int cat, int quiet,
			   ListOfLists *results) {
  String *fname = str_new(STR_MED_LEN);
  FILE *POSTPROBF, *EXPSUBF, *EXPTOTSUBF;
  int i, j, tup, node, state, state2;
  TreeNode *n;
  char tuplestr[mod->order+2];
  char coltupstr[msa->nseqs+1];
  int ratecat;
  char ***dimnames;

  if (msa->ss == NULL) 
    die("Error: print_post_prob_stats needs sufficient statistics");
  if (do_every_site && msa->ss->tuple_idx == NULL)
    die("Error in print_post_prob_stats: do_every_site option requires ordered sufficient statistics");

  tuplestr[mod->order+1] = '\0';
  coltupstr[msa->nseqs] = '\0';

  /* FIXME: rate variation!     need rate post probs! */
  if (mod->nratecats != 1)
    die("ERROR print_post_prob_stats nratecats should be 1 but is %i\n",
	mod->nratecats);

  /* compute desired stats */
  if (mod->tree_posteriors != NULL)
    die("ERROR: mod->tree_posteriors should be NULL\n");
  if (!quiet) 
    fprintf(stderr, "Computing posterior probabilities and/or related stats ...\n");
  mod->tree_posteriors = tl_new_tree_posteriors(mod, msa, do_bases, 0, 
                                                do_expected_nsubst, 
                                                do_expected_nsubst_tot, 
						do_expected_nsubst_col,
						0, 0);
  tl_compute_log_likelihood(mod, msa, NULL, NULL, cat, mod->tree_posteriors);
  tr_name_ancestors(mod->tree);

  if (do_bases) {
    if (results != NULL) {
      double ****arr;
      int dimsize[4];
      dimnames = smalloc(4*sizeof(char**));
      dimnames[0] = get_ratecat_names(mod, &dimsize[0]);
      dimnames[1] = get_node_names(mod, 0, 1, 0, &dimsize[1]);
      dimnames[2] = get_state_names(mod, NULL, &dimsize[2]);
      if (do_every_site)
	dimnames[3] = get_site_names(msa, &dimsize[3]);
      else dimnames[3] = get_tuple_names(mod, msa, cat, &dimsize[3]);

      arr = (double****)alloc_n_dimensional_array(4, dimsize, sizeof(double));
      for (ratecat=0; ratecat < mod->nratecats; ratecat++) {
	int node_idx = 0;
	for (node = 0; node < mod->tree->nnodes; node++) {
	  n = (TreeNode*)lst_get_ptr(mod->tree->nodes, node);
	  if (n->lchild == NULL || n->rchild == NULL) continue;
	  if (do_every_site) {
	    for (i=0; i < msa->length; i++) {
	      for (state=0; state < mod->rate_matrix->size; state++) {
		arr[ratecat][node_idx][state][i] = mod->tree_posteriors->base_probs[ratecat][state][n->id][msa->ss->tuple_idx[i]];
	      }
	    }
	  } else {
	    for (state = 0 ; state < mod->rate_matrix->size; state++) {
	      int tup_idx=0;
	      for (tup=0; tup < msa->ss->ntuples; tup++) {
		if ((cat >=0 && msa->ss->cat_counts[cat][tup] == 0) ||
		    msa->ss->counts[tup] == 0) continue;
		arr[ratecat][node_idx][state][tup_idx++] = mod->tree_posteriors->base_probs[ratecat][state][n->id][tup];
	      }
	    }
	  }
	  node_idx++;
	}
      }
      lol_push_dbl_array(results, arr, "bybase", 4, dimsize, dimnames);
      free_n_dimensional_array(arr, 4, dimsize);
      for (i=0; i < 4; i++) {
	for (j=0; j < dimsize[i]; j++)
	  sfree(dimnames[i][j]);
	sfree(dimnames[i]);
      }
      sfree(dimnames);
    }
    if (output_fname_root != NULL) {
      set_output_fname(fname, output_fname_root, cat, ".postprob");
      if (!quiet) 
	fprintf(stderr, "Writing posterior probabilities to %s ...\n", 
		fname->chars);

      POSTPROBF = phast_fopen(fname->chars, "w+");
      
      /* print header */
      fprintf(POSTPROBF, "%-6s ", "#");
      for (i = 0; i < msa->nseqs; i++) fprintf(POSTPROBF, " ");
      fprintf(POSTPROBF, "    ");
      for (node = 0; node < mod->tree->nnodes; node++) {
	n = lst_get_ptr(mod->tree->nodes, node);
	if (n->lchild == NULL || n->rchild == NULL) continue;
	for (state = 0; state < mod->rate_matrix->size; state++) {
	  if (state == mod->rate_matrix->size/2)
	    fprintf(POSTPROBF, "node %-2d", n->id);
	  else
	    fprintf(POSTPROBF, "%6s ", "");
	}
      }
      fprintf(POSTPROBF, "\n%-6s ", "#");
      for (state = 0; state < msa->nseqs-5; state++) fprintf(POSTPROBF, " ");
      fprintf(POSTPROBF, "tuple     ");
      for (node = 0; node < mod->tree->nnodes; node++) {
	n = lst_get_ptr(mod->tree->nodes, node);
	if (n->lchild == NULL || n->rchild == NULL) continue;
	for (state = 0; state < mod->rate_matrix->size; state++) {
	  get_tuple_str(tuplestr, state, mod->order + 1, 
			mod->rate_matrix->states);
	  fprintf(POSTPROBF, "%6s ", tuplestr);
	}
      }
      fprintf(POSTPROBF, "\n");
      
      /* print post probs */
      for (tup = 0; tup < msa->ss->ntuples; tup++) {
	checkInterruptN(tup, 1000);
	
	if ((cat >= 0 && msa->ss->cat_counts[cat][tup] == 0) ||
	    msa->ss->counts[tup] == 0) continue;
	
	tuple_to_string_pretty(coltupstr, msa, tup);
	fprintf(POSTPROBF, "%-6d %5s      ", tup, coltupstr);
	for (node = 0; node < mod->tree->nnodes; node++) {
	  n = lst_get_ptr(mod->tree->nodes, node);
	  if (n->lchild == NULL || n->rchild == NULL) continue;
	  for (state = 0; state < mod->rate_matrix->size; state++) 
	    fprintf(POSTPROBF, "%6.4f ", 
		    mod->tree_posteriors->base_probs[0][state][n->id][tup]);
	}                 
	fprintf(POSTPROBF, "\n");
      }
      phast_fclose(POSTPROBF);
    }
  }

  if (do_expected_nsubst) {
    if (results != NULL) {
      int dimsize[3], node_idx, tup_idx;
      double ***arr;

      dimnames = smalloc(3*sizeof(char**));
      dimnames[0] = get_ratecat_names(mod, &dimsize[0]);
      dimnames[1] = get_node_names(mod, 1, 0, 1, &dimsize[1]);
      dimnames[2] = get_tuple_names(mod, msa, cat, &dimsize[2]);

      arr = alloc_n_dimensional_array(3, dimsize, sizeof(double));
      for (ratecat=0; ratecat < mod->nratecats; ratecat++) {
	tup_idx = 0;
	for (tup=0; tup < msa->ss->ntuples; tup++) {
	  if ((cat >= 0 && msa->ss->cat_counts[cat][tup] == 0) ||
	      msa->ss->counts[tup] == 0) continue;
	  arr[ratecat][0][tup_idx++] = cat >= 0 ? msa->ss->cat_counts[cat][tup] : msa->ss->counts[tup];
	}
	node_idx = 1;
	for (node=0; node < mod->tree->nnodes; node++) {
	  tup_idx = 0;
	  n = lst_get_ptr(mod->tree->nodes, node);
	  if (n == mod->tree) continue;
	  for (tup=0; tup < msa->ss->ntuples; tup++) {
	    checkInterruptN(tup, 1000);
	    if ((cat >= 0 && msa->ss->cat_counts[cat][tup] == 0) ||
		msa->ss->counts[tup] == 0) continue;
	    arr[ratecat][node_idx][tup_idx++] = mod->tree_posteriors->expected_nsubst[ratecat][n->id][tup];
	  }
	  node_idx++;
	}
      }
      lol_push_dbl_array(results, arr, "exp.nsub", 3, dimsize, dimnames);
      free_n_dimensional_array(arr, 3, dimsize);
      for (i=0; i < 3; i++) {
	for (j=0; j < dimsize[i]; j++)
	  sfree(dimnames[i][j]);
	sfree(dimnames[i]);
      }
      sfree(dimnames);
    }
    if (output_fname_root != NULL) {
      set_output_fname(fname, output_fname_root, cat, ".expsub");
      if (!quiet) 
	fprintf(stderr, "Writing expected numbers of substitutions to %s ...\n", 
		fname->chars);
      EXPSUBF = phast_fopen(fname->chars, "w+");
      
      fprintf(EXPSUBF, "%-3s %10s %7s ", "#", "tuple", "count");
      for (node = 0; node < mod->tree->nnodes; node++) {
	n = lst_get_ptr(tr_postorder(mod->tree), node);
	if (n == mod->tree) continue;
	fprintf(EXPSUBF, " node_%-2d", n->id);
      }
      fprintf(EXPSUBF, "    total\n");
      for (tup = 0; tup < msa->ss->ntuples; tup++) {
	double total = 0;
	checkInterruptN(tup, 1000);
	
	if ((cat >= 0 && msa->ss->cat_counts[cat][tup] == 0) ||
	    msa->ss->counts[tup] == 0) continue;
	
	tuple_to_string_pretty(coltupstr, msa, tup);
	fprintf(EXPSUBF, "%-3d %10s %.0f ", tup, coltupstr, msa->ss->counts[tup]);
	for (node = 0; node < mod->tree->nnodes; node++) {
	  n = lst_get_ptr(tr_postorder(mod->tree), node);
	  if (n == mod->tree) continue;
	  fprintf(EXPSUBF, "%7.4f ", 
		  mod->tree_posteriors->expected_nsubst[0][n->id][tup]);
	  total += mod->tree_posteriors->expected_nsubst[0][n->id][tup];
	}                 
	fprintf(EXPSUBF, "%7.4f\n", total);
      }
      phast_fclose(EXPSUBF);
    }
  }

  if (do_expected_nsubst_col) {
    if (results != NULL) {
      int dimsize[5];
      double *****arr;
      dimnames = smalloc(5*sizeof(char**));
      dimnames[0] = get_ratecat_names(mod, &dimsize[0]);
      dimnames[1] = get_node_names(mod, 1, 0, 1, &dimsize[1]);
      dimnames[2] = get_tuple_names(mod, msa, cat, &dimsize[2]);
      dimnames[3] = get_state_names(mod, "from.", &dimsize[3]);
      dimnames[4] = get_state_names(mod, "to.", &dimsize[4]);
      
      arr = alloc_n_dimensional_array(5, dimsize, sizeof(double));
      for (ratecat=0; ratecat < mod->nratecats; ratecat++) {
	int node_idx = 1;
	for (node = 0; node < mod->tree->nnodes; node++) {
	  int tuple_idx=0;
	  n = lst_get_ptr(mod->tree->nodes, node);
	  if (n == mod->tree) continue;
	  for (tup=0; tup < msa->ss->ntuples; tup++) {
	    if ((cat >= 0 && msa->ss->cat_counts[cat][tup] == 0) ||
		msa->ss->counts[tup] == 0) continue;
	    for (state=0; state < mod->rate_matrix->size; state++) {
	      for (state2=0; state2 < mod->rate_matrix->size; state2++) {
		arr[ratecat][0][tuple_idx][state][state2] = cat >= 0 ? msa->ss->cat_counts[cat][tup] : msa->ss->counts[tup];
		arr[ratecat][node_idx][tuple_idx][state][state2] = mod->tree_posteriors->expected_nsubst_col[ratecat][n->id][tup][state][state2];
	      }
	    }
	    tuple_idx++;
	  }
	  node_idx++;
	}
      }
      lol_push_dbl_array(results, arr, "col.exp.nsub", 5, dimsize, dimnames);
      free_n_dimensional_array(arr, 5, dimsize);
      for (i=0; i < 5; i++) {
	for (j=0; j<dimsize[i]; j++)
          sfree(dimnames[i][j]);
	sfree(dimnames[i]);
      }
      sfree(dimnames);
    }
    if (output_fname_root != NULL) {
      set_output_fname(fname, output_fname_root, cat, ".expcolsub");
      if (!quiet) 
	fprintf(stderr, "Writing expected numbers of substitutions per site to %s ...\n",
		fname->chars);
      EXPSUBF = phast_fopen(fname->chars, "w+");

      /* print header */
      fprintf(EXPSUBF, "#tuple\tcount\tbranch");
      for (state=0; state < mod->rate_matrix->size; state++) {
	for (state2=0; state2 < mod->rate_matrix->size; state2++) {
	  get_tuple_str(tuplestr, state, mod->order+1, mod->rate_matrix->states);
	  fprintf(EXPSUBF, "\t%s", tuplestr);
	  get_tuple_str(tuplestr, state2, mod->order+1, mod->rate_matrix->states);
	  fprintf(EXPSUBF, "->%s", tuplestr);
	}
      }
      fprintf(EXPSUBF, "\n");

      for (tup=0; tup < msa->ss->ntuples; tup++) {
	if ((cat >= 0 && msa->ss->cat_counts[cat][tup] == 0) ||
	    msa->ss->counts[tup] == 0) continue;
	tuple_to_string_pretty(coltupstr, msa, tup);
	for (node = 0; node < mod->tree->nnodes; node++) {
	  n = lst_get_ptr(mod->tree->nodes, node);
	  if (n->parent == NULL) continue;
	  fprintf(EXPSUBF, "%s\t%.0f\t%s", coltupstr, 
		  cat >=0 ? msa->ss->cat_counts[cat][tup] : 
		  msa->ss->counts[tup], 
		  n->name);
	  for (state=0; state < mod->rate_matrix->size; state++) {
	    for (state2=0; state2 < mod->rate_matrix->size; state2++) {
	      fprintf(EXPSUBF, "\t%g", mod->tree_posteriors->expected_nsubst_col[0][n->id][tup][state][state2]);
	    }
	  }
	  fprintf(EXPSUBF, "\n");
	}
      }
      phast_fclose(EXPSUBF);
    }
  }

  if (do_expected_nsubst_tot) {
    if (results != NULL) {
      double ****arr;
      int dimsize[4];

      dimnames = smalloc(4*sizeof(char**));
      dimnames[0] = get_ratecat_names(mod, &dimsize[0]);
      dimnames[1] = get_node_names(mod, 1, 0, 0, &dimsize[1]);
      dimnames[2] = get_state_names(mod, "from.", &dimsize[2]);
      dimnames[3] = get_state_names(mod, "to.", &dimsize[3]);
      
      arr = alloc_n_dimensional_array(4, dimsize, sizeof(double));
      for (ratecat=0; ratecat < mod->nratecats; ratecat++) {
	int node_idx=0;
	for (node = 0 ; node < mod->tree->nnodes; node++) {
	  n = lst_get_ptr(mod->tree->nodes, node);
	  if (n == mod->tree) continue;
	  for (state=0; state < mod->rate_matrix->size; state++) {
	    for (state2=0; state2 < mod->rate_matrix->size; state2++) 
	      arr[ratecat][node_idx][state][state2] = mod->tree_posteriors->expected_nsubst_tot[ratecat][state][state2][n->id];
	  }
	  node_idx++;
	}
      }
      lol_push_dbl_array(results, arr, "tot.exp.nsub", 4, dimsize, dimnames);
      free_n_dimensional_array(arr, 4, dimsize);
      for (i=0; i < 4; i++) {
	for (j=0; j < dimsize[i]; j++)
	  sfree(dimnames[i][j]);
	sfree(dimnames[i]);
      }
      sfree(dimnames);
    }
    if (output_fname_root != NULL) {
      set_output_fname(fname, output_fname_root, cat, ".exptotsub");
      if (!quiet) 
	fprintf(stderr, "Writing total expected numbers of substitutions to %s ...\n", 
		fname->chars);
      EXPTOTSUBF = phast_fopen(fname->chars, "w+");
      
      fprintf(EXPTOTSUBF, "\n\
A separate matrix of expected numbers of substitutions is shown for each\n\
branch of the tree.     Nodes of the tree are visited in a postorder traversal,\n\
and each node is taken to be representative of the branch between itself and\n\
its parent.     Starting bases or tuples of bases appear on the vertical axis\n\
of each matrix, and destination bases or tuples of bases appear on the\n\
horizontal axis.\n\n");
      
      for (node = 0; node < mod->tree->nnodes; node++) {
	n = lst_get_ptr(tr_postorder(mod->tree), node);
	if (n == mod->tree) continue;
	
	fprintf(EXPTOTSUBF, "Branch above node %d", n->id);
	if (n->name != NULL && strlen(n->name) > 0) 
	  fprintf(EXPTOTSUBF, " (leaf labeled '%s')", n->name);
	fprintf(EXPTOTSUBF, ":\n\n");
	
	/* print header */
	fprintf(EXPTOTSUBF, "%-4s ", "");
	for (state2 = 0; state2 < mod->rate_matrix->size; state2++) {
	  get_tuple_str(tuplestr, state2, mod->order + 1, 
			mod->rate_matrix->states);
	  fprintf(EXPTOTSUBF, "%12s ", tuplestr);
	}
	fprintf(EXPTOTSUBF, "\n");
	for (state = 0; state < mod->rate_matrix->size; state++) {
	  get_tuple_str(tuplestr, state, mod->order + 1, 
			mod->rate_matrix->states);
	  fprintf(EXPTOTSUBF, "%-4s ", tuplestr);
	  for (state2 = 0; state2 < mod->rate_matrix->size; state2++) 
	    fprintf(EXPTOTSUBF, "%12.2f ", 
		    mod->tree_posteriors->expected_nsubst_tot[0][state][state2][n->id]);
	  fprintf(EXPTOTSUBF, "\n");
	}
	fprintf(EXPTOTSUBF, "\n\n");
      }
      phast_fclose(EXPTOTSUBF);
    }
  }

  tl_free_tree_posteriors(mod, msa, mod->tree_posteriors);
  mod->tree_posteriors = NULL;
  str_free(fname);
}


void print_window_summary(FILE* WINDOWF, List *window_coords, int win, 
                          int cat, TreeModel *mod, double *gc,
                          int ninf_sites, int nseqs, int header_only) {
  int j, i;
  if (header_only) {
    if (WINDOWF != NULL) {
      fprintf(WINDOWF, "%5s %8s %8s %4s", "win", "beg", "end", "cat");
      for (i=0; i < nseqs; i++) {
	char temp[100];
	sprintf(temp, "GC%i", i);
	fprintf(WINDOWF, " %6s", temp);
      }
      fprintf(WINDOWF, " %7s", "ninf");
      fprintf(WINDOWF, " %7s\n", "t");
    }
  }
  else {
    if (WINDOWF != NULL) {
      fprintf(WINDOWF, "%5d %8d %8d %4d", win/2+1, 
	      lst_get_int(window_coords, win), 
	      lst_get_int(window_coords, win+1), cat);
      fprintf(WINDOWF, " %6.4f", 
	      vec_get(mod->backgd_freqs, 
		      mod->rate_matrix->inv_states[(int)'G']) + 
	      vec_get(mod->backgd_freqs, 
		      mod->rate_matrix->inv_states[(int)'C']));
      for (j = 0; j < nseqs; j++) 
	fprintf(WINDOWF, " %6.4f", gc==NULL ? -1.0 : gc[j]);
      fprintf(WINDOWF, " %7d", ninf_sites);
      fprintf(WINDOWF, " %7.4f\n", tr_total_len(mod->tree));
    }
  }
}



int run_phyloFit(struct phyloFit_struct *pf) {
  FILE *F, *WINDOWF=NULL;
  int i, j, win, root_leaf_id = -1;
  String *mod_fname;
  MSA *source_msa;
  String *tmpstr = str_new(STR_SHORT_LEN);
  List *cats_to_do=NULL;
  double *gc=NULL;
  char tmpchstr[STR_MED_LEN];
  FILE *parsimony_cost_file = NULL;
  int free_cm = FALSE, free_cats_to_do_str=FALSE, free_tree=FALSE,
    free_window_coords = FALSE;

  //copy some heavily used variables directly from pf
  MSA *msa = pf->msa;
  int subst_mod = pf->subst_mod;
  TreeNode *tree = pf->tree;
  GFF_Set *gff = pf->gff;
  int quiet = pf->quiet;
  TreeModel *input_mod = pf->input_mod;
  FILE *error_file=NULL;

  if (pf->no_freqs)
    pf->init_backgd_from_data = FALSE;

  if (pf->parsimony_cost_fname != NULL)
    parsimony_cost_file = phast_fopen(pf->parsimony_cost_fname, "w");

  if (pf->use_conditionals && pf->use_em) 
    die("ERROR: Cannot use --markov with --EM.    Type %s for usage.\n",
	pf->see_for_help);
  
  if (pf->likelihood_only && input_mod == NULL)  
    die("ERROR: --lnl requires --init-model.  Type '%s' for usage.\n",
	pf->see_for_help);

  if (input_mod != NULL && tree != NULL)
    die("ERROR: --tree is not allowed with --init-model.\n");

  if (subst_mod == UNDEF_MOD) {
    if (pf->input_mod != NULL)
      subst_mod = pf->input_mod->subst_mod;
    else subst_mod = REV;
  }

  if (pf->gaps_as_bases && subst_mod != JC69 && subst_mod != F81 && 
      subst_mod != HKY85G && subst_mod != REV && 
      subst_mod != UNREST && subst_mod != SSREV)
    die("ERROR: --gaps-as-bases currently only supported with JC69, F81, HKY85+Gap, REV, SSREV, and UNREST.\n");
                                /* with HKY, yields undiagonalizable matrix */
  if ((pf->no_freqs || pf->no_rates) && input_mod == NULL)
    die("ERROR: --init-model required with --no-freqs and/or --no-rates.\n");

  if (pf->no_freqs && pf->estimate_backgd)
    die("ERROR: can't use both --no-freqs and --estimate-freqs.\n");

  if (gff != NULL && pf->cm == NULL) {
    pf->cm = cm_new_from_features(gff);
    free_cm = TRUE;
  }

  if (pf->subtree_name != NULL && pf->estimate_scale_only == FALSE) {
    if (!quiet) 
      fprintf(stderr, "warning: specifying subtree implies scale_only\n");
    pf->estimate_scale_only = TRUE;
  }

  if (pf->rate_consts != NULL) {
    lst_qsort_dbl(pf->rate_consts, ASCENDING);
    if (lst_size(pf->rate_consts) < 2 || lst_get_dbl(pf->rate_consts, 0) <= 0) 
      die("ERROR: must be >= 2 rate constants and all must be positive.\n");
    if (pf->nratecats != lst_size(pf->rate_consts))
      die("ERROR: rate_consts must have length nratecats");
  }

  /* internally, --non-overlapping is accomplished via --do-cats */
  if (pf->nonoverlapping) {
    if (pf->cats_to_do_str != NULL)
      die("ERROR: cannot use --do-cats with nonoverlapping");
    pf->cats_to_do_str = lst_new_ptr(1);
    lst_push_ptr(pf->cats_to_do_str, str_new_charstr("1"));
    free_cats_to_do_str = TRUE;
  }

  if (tree == NULL) {
    if (input_mod != NULL) tree = input_mod->tree;
    else if (msa->nseqs == 2) {
      sprintf(tmpchstr, "(%s,%s)", msa->names[0], msa->names[1]);
      tree = tr_new_from_string(tmpchstr);
      free_tree = TRUE;
    }
    else if (msa->nseqs == 3 && subst_mod_is_reversible(subst_mod) && pf->alt_mod_str == NULL) {
      sprintf(tmpchstr, "(%s,(%s,%s))", msa->names[0], msa->names[1], 
              msa->names[2]);
      tree = tr_new_from_string(tmpchstr);
      free_tree = TRUE;
    }
    else die("ERROR: --tree required.\n");
  }


  /* allow for specified ancestor */
  if (pf->root_seqname != NULL) {
    TreeNode *rl;
    if (tree == NULL || subst_mod_is_reversible(subst_mod)) 
      die("ERROR: --ancestor requires --tree and a non-reversible model.\n");
    rl = tr_get_node(tree, pf->root_seqname);     
    if (rl == NULL || rl->parent != tree) 
      die("ERROR: Sequence specified by --ancestor must be a child of the root.\n");
    root_leaf_id = rl->id;
  }

  if (pf->label_str != NULL || pf->label_type != NULL) {
    if (pf->label_str == NULL || pf->label_type == NULL ||
	lst_size(pf->label_str) != lst_size(pf->label_type))
      //shouldn't happen unless bug
      die("label_str and label_type should both be lists of same length\n");
    for (i=0; i < lst_size(pf->label_str); i++) {
      String *currstr = (String*)lst_get_ptr(pf->label_str, i), *arg1, *label;
      List *tmplst = lst_new_ptr(10);
      String *nodename;
      str_split(currstr, ":", tmplst);
      if (lst_size(tmplst) != 2) 
	die("ERROR: bad argument to --label-branches or --label-subtree.\n");
      arg1 = lst_get_ptr(tmplst, 0);
      label = lst_get_ptr(tmplst, 1);
      lst_clear(tmplst);
      if (lst_get_int(pf->label_type, i) == BRANCH_TYPE) {
	str_split(arg1, ",", tmplst);
	for (j=0; j < lst_size(tmplst); j++) {
	  nodename = (String*)lst_get_ptr(tmplst, j);
	  tr_label_node(tree, nodename->chars, label->chars);
	}
	lst_free_strings(tmplst);
      } else if (lst_get_int(pf->label_type, i) == SUBTREE_TYPE) {
	int include_leading_branch = FALSE;
	TreeNode *node;
	nodename = arg1;
	node = tr_get_node(tree, nodename->chars);
	if (node == NULL && nodename->chars[nodename->length-1] == '+') {
	  nodename->chars[--nodename->length] = '\0';
	  node = tr_get_node(tree, nodename->chars);
	  include_leading_branch = TRUE;
	}
	tr_label_subtree(tree, nodename->chars, include_leading_branch, 
			 label->chars);
      } else die("ERROR got label_type %i\n", lst_get_int(pf->label_type, i));
      str_free(arg1);
      str_free(label);
      lst_free(tmplst);
    }
  }

  
  if (msa_alph_has_lowercase(msa)) msa_toupper(msa); 
  msa_remove_N_from_alph(msa);    /* for backward compatibility */

  /* set up for categories */
  /* first label sites, if necessary */
  if (gff != NULL && pf->label_categories) {
    if (msa->idx_offset > 0) {
      /* if these's an offset, we'll just subtract it from all features */
      for (i = 0; i < lst_size(gff->features); i++) {
        GFF_Feature *f = lst_get_ptr(gff->features, i);
        f->start -= msa->idx_offset;
        f->end -= msa->idx_offset;
      }
      msa->idx_offset = 0;
    }

    /* convert GFF to coordinate frame of alignment */
    msa_map_gff_coords(msa, gff, 1, 0, 0);

    /* reverse complement segments of MSA corresponding to features on
       reverse strand (if necessary) */
    if (pf->reverse_group_tag != NULL) {
      gff_group(gff, pf->reverse_group_tag);
      msa_reverse_compl_feats(msa, gff, NULL);
    }
    
    /* label categories */
    if (!quiet) fprintf(stderr, "Labeling alignment sites by category ...\n");
    msa_label_categories(msa, gff, pf->cm);
  }
  else if (pf->nonoverlapping && pf->label_categories) {
                                /* (already taken care of if MAF) */
    int cycle_size = tm_order(subst_mod) + 1;
    if (!(msa->seqs != NULL && msa->ss == NULL))
      die("ERROR run_phyloFit: need explicit sequences, not sufficient statistics\n");
    msa->categories = smalloc(msa->length * sizeof(int));
    for (i = 0; i < msa->length; i++) 
      msa->categories[i] = (i % cycle_size) + 1;
    msa->ncats = cycle_size;
  }
  /* at this point, we have msa->ncats > 0 iff we intend to do
     category-by-category estimation */

  /* now set up list of categories to process.    There are several
     cases to consider */
  if (msa->ncats < 0) {            
    if (pf->cats_to_do_str != NULL)
      fprintf(stderr, "WARNING: ignoring --do-cats; no category information.\n");
    cats_to_do = lst_new_int(1);
    lst_push_int(cats_to_do, -1);
    /* no categories -- pool all sites */
  }
  else if (pf->cats_to_do_str == NULL) {
    cats_to_do = lst_new_int(msa->ncats + 1);
    for (i = 0; i <= msa->ncats; i++) lst_push_int(cats_to_do, i);
                                /* have categories but no --do-cats --
                                   process all categories */
  }
  else if (pf->cm != NULL) 
    cats_to_do = cm_get_category_list(pf->cm, pf->cats_to_do_str, 0);
                                /* have --do-cats and category map;
                                   use cm_get_category_list (allows
                                   use of names as well as numbers) */
  else if (pf->cats_to_do_str != NULL)
    cats_to_do = str_list_as_int(pf->cats_to_do_str);
                                /* have --do-cats but no category map;
                                   use literal numbers */
  /* set up windows, if necessary */
  if (pf->window_size != -1) {
    if (pf->window_coords != NULL) 
      die("ERROR: cannot use both --windows and --windows-explicit.\n");
    pf->window_coords = lst_new_int(msa->length/pf->window_shift + 1);
    for (i = 1; i < msa->length; i += pf->window_shift) {
      lst_push_int(pf->window_coords, i);
      lst_push_int(pf->window_coords, 
                   min(i + pf->window_size - 1, msa->length));
    }
    free_window_coords = TRUE;
  }
  if (pf->window_coords != NULL) {
    /* set up summary file */
    String *sumfname=NULL;
    msa_coord_map *map;
    if (pf->output_fname_root != NULL) {
      sumfname = str_new_charstr(pf->output_fname_root);
      str_append_charstr(sumfname, ".win-sum");
      WINDOWF = phast_fopen(sumfname->chars, "w+");
      str_free(sumfname);
    } 
    print_window_summary(WINDOWF, NULL, 0, 0, NULL, NULL, 0, 0, TRUE);
    
    /* map to coord frame of alignment */
    map = msa_build_coord_map(msa, 1);
    for (i = 0; i < lst_size(pf->window_coords); i += 2) {
      lst_set_int(pf->window_coords, i, 
                  msa_map_seq_to_msa(map, lst_get_int(pf->window_coords, i)));
      lst_set_int(pf->window_coords, i+1, 
                  msa_map_seq_to_msa(map, lst_get_int(pf->window_coords, i+1)));
    }
    msa_map_free(map);
  }

  if (pf->error_fname != NULL)
    error_file = phast_fopen(pf->error_fname, "w");
  
  /* now estimate models (window by window, if necessary) */
  mod_fname = str_new(STR_MED_LEN);
  source_msa = msa;
  for (win = 0; 
       win < (pf->window_coords == NULL ? 1 : lst_size(pf->window_coords)); 
       win += 2) {
    int win_beg, win_end;

    if (pf->window_coords != NULL) {
      win_beg = lst_get_int(pf->window_coords, win);
      win_end = lst_get_int(pf->window_coords, win+1);
      if (win_beg < 0 || win_end < 0) continue;
      
      /* note: msa_sub_alignment uses a funny indexing system (see docs) */
      msa = msa_sub_alignment(source_msa, NULL, 0, win_beg-1, win_end);
    }

    /* process each category */
    for (i = 0; i < lst_size(cats_to_do); i++) {
      TreeModel *mod;
      Vector *params = NULL;
      List *pruned_names;
      int old_nnodes, cat = lst_get_int(cats_to_do, i);
      unsigned int ninf_sites;

      if (input_mod == NULL) 
        mod = tm_new(tr_create_copy(tree), NULL, NULL, subst_mod, 
                     msa->alphabet, pf->nratecats == -1 ? 1 : pf->nratecats, 
		     pf->alpha, pf->rate_consts, root_leaf_id);
      else if (pf->likelihood_only)
        mod = input_mod;
      else {
	List *rate_consts, *freq;
	double alpha;
	int nratecats;

	if (pf->nratecats != -1) {
	  nratecats = pf->nratecats;
	  alpha = pf->alpha;
	  rate_consts = pf->rate_consts;
	  freq = NULL;
	} else {
	  nratecats = input_mod->nratecats;
	  alpha = input_mod->alpha;
	  if (input_mod->rK != NULL) {
	    rate_consts = lst_new_dbl(input_mod->nratecats);
	    for (j=0; j < input_mod->nratecats; j++)
	      lst_push_dbl(rate_consts, input_mod->rK[j]);
	  } else rate_consts = NULL;
	  if (input_mod->freqK != NULL) {
	    freq = lst_new_dbl(input_mod->nratecats);
	    for (j=0; j < input_mod->nratecats; j++)
	      lst_push_dbl(freq, input_mod->freqK[j]);
	  } else freq = NULL;
	}
        mod = input_mod;
        tm_reinit(mod, subst_mod, nratecats, alpha, 
		  rate_consts, freq);
	if (rate_consts != pf->rate_consts)
	  lst_free(rate_consts);
	if (freq != NULL) 
	  lst_free(freq);
      }

      if (pf->use_selection) {
	mod->selection_idx = 0;
	mod->selection = pf->selection;
      }

      mod->noopt_arg = pf->nooptstr == NULL ? NULL : str_new_charstr(pf->nooptstr->chars);
      mod->eqfreq_sym = pf->symfreq || subst_mod == SSREV;
      if (pf->bound_arg != NULL) {
	mod->bound_arg = lst_new_ptr(lst_size(pf->bound_arg));
	for (j=0; j < lst_size(pf->bound_arg); j++) {
	  String *tmp = lst_get_ptr(pf->bound_arg, j);
	  lst_push_ptr(mod->bound_arg, str_new_charstr(tmp->chars));
	}
      } else mod->bound_arg = NULL;

      mod->use_conditionals = pf->use_conditionals;

      if (pf->estimate_scale_only || 
	  pf->estimate_backgd || 
	  pf->no_rates || 
	  pf->assume_clock) {
        if (pf->estimate_scale_only) {
          mod->estimate_branchlens = TM_SCALE_ONLY;

          if (pf->subtree_name != NULL) { /* estimation of subtree scale */
            String *s1 = str_new_charstr(pf->subtree_name), 
              *s2 = str_new_charstr(pf->subtree_name);
            str_root(s1, ':'); str_suffix(s2, ':'); /* parse string */
            mod->subtree_root = tr_get_node(mod->tree, s1->chars);
            if (mod->subtree_root == NULL) {
	      tr_name_ancestors(mod->tree);
	      mod->subtree_root = tr_get_node(mod->tree, s1->chars);
	      if (mod->subtree_root == NULL)
		die("ERROR: no node named '%s'.\n", s1->chars);
	    }
            if (s2->length > 0) {
              if (str_equals_charstr(s2, "loss")) 
		mod->scale_sub_bound = LB;
              else if (str_equals_charstr(s2, "gain")) 
		mod->scale_sub_bound = UB;
              else die("ERROR: unrecognized suffix '%s'\n", s2->chars);
            }
            str_free(s1); str_free(s2);
          }
        }
	
        else if (pf->assume_clock)
          mod->estimate_branchlens = TM_BRANCHLENS_CLOCK;
        
        if (pf->no_rates)
          mod->estimate_ratemat = FALSE;

        mod->estimate_backgd = pf->estimate_backgd;
      }
      
      if (pf->no_branchlens)
	mod->estimate_branchlens = TM_BRANCHLENS_NONE;

      if (pf->ignore_branches != NULL) 
        tm_set_ignore_branches(mod, pf->ignore_branches);

      old_nnodes = mod->tree->nnodes;
      pruned_names = lst_new_ptr(msa->nseqs);
      tm_prune(mod, msa, pruned_names);
      if (lst_size(pruned_names) == (old_nnodes + 1) / 2)
        die("ERROR: no match for leaves of tree in alignment (leaf names must match alignment names).\n");
      if (!quiet && lst_size(pruned_names) > 0) {
        fprintf(stderr, "WARNING: pruned away leaves of tree with no match in alignment (");
        for (j = 0; j < lst_size(pruned_names); j++)
          fprintf(stderr, "%s%s", ((String*)lst_get_ptr(pruned_names, j))->chars, 
                  j < lst_size(pruned_names) - 1 ? ", " : ").\n");
      }
      lst_free_strings(pruned_names);
      lst_free(pruned_names);

      if (pf->alt_mod_str != NULL) {
	for (j = 0 ; j < lst_size(pf->alt_mod_str); j++) 
	  tm_add_alt_mod(mod, (String*)lst_get_ptr(pf->alt_mod_str, j));
      }

      str_clear(tmpstr);
      
      if  (pf->msa_fname != NULL)
	str_append_charstr(tmpstr, pf->msa_fname);
      else str_append_charstr(tmpstr, "alignment");
      
      if (cat != -1 || pf->window_coords != NULL) {
	str_append_charstr(tmpstr, " (");
	if (cat != -1) {
	  str_append_charstr(tmpstr, "category ");
	  str_append_int(tmpstr, cat);
	}
	
	if (pf->window_coords != NULL) {
	  if (cat != -1) str_append_charstr(tmpstr, ", ");
	  str_append_charstr(tmpstr, "window ");
	  str_append_int(tmpstr, win/2 + 1);
	}
	
	str_append_char(tmpstr, ')');
      }

      ninf_sites = msa_ninformative_sites(msa, cat);
      if (ninf_sites < pf->nsites_threshold) {
        if (input_mod == NULL) tm_free(mod);
        fprintf(stderr, "Skipping %s; insufficient informative sites ...\n", 
                tmpstr->chars);
        continue;
      }

      if (pf->init_parsimony) {
	double parsimony_cost = tm_params_init_branchlens_parsimony(NULL, mod, msa, cat);
        if (parsimony_cost_file != NULL) 
           fprintf(parsimony_cost_file, "%f\n", parsimony_cost);
        if (pf->parsimony_only) continue;
      }

      if (pf->likelihood_only) {
        double *col_log_probs = pf->do_column_probs ? 
          smalloc(msa->length * sizeof(double)) : NULL;
        String *colprob_fname;
        if (!quiet) 
          fprintf(stderr, "Computing likelihood of %s ...\n", tmpstr->chars);
        tm_set_subst_matrices(mod);
        if (pf->do_column_probs && msa->ss != NULL && msa->ss->tuple_idx == NULL) {
          msa->ss->tuple_idx = smalloc(msa->length * sizeof(int));
          for (j = 0; j < msa->length; j++)
            msa->ss->tuple_idx[j] = j;
        }
        mod->lnL = tl_compute_log_likelihood(mod, msa, col_log_probs, NULL, cat, NULL) * 
          log(2);
        if (pf->do_column_probs) {
	  //we don't need to implement this in RPHAST because there is
	  //already a msa.likelihood function
	  if (pf->output_fname_root == NULL)
	    die("ERROR: currently do_column_probs requires output file");
          colprob_fname = str_new_charstr(pf->output_fname_root);
          str_append_charstr(colprob_fname, ".colprobs");
          if (!quiet) 
            fprintf(stderr, "Writing column probabilities to %s ...\n", 
                    colprob_fname->chars);
          F = phast_fopen(colprob_fname->chars, "w+");
          for (j = 0; j < msa->length; j++)
            fprintf(F, "%d\t%.6f\n", j, col_log_probs[j]);
          phast_fclose(F);
          str_free(colprob_fname);
          sfree(col_log_probs);
        }
      }
      else {                    /* fit model */

        if (msa->ss == NULL) {    /* get sufficient stats if necessary */
          if (!quiet)
            fprintf(stderr, "Extracting sufficient statistics ...\n");
          ss_from_msas(msa, mod->order+1, 0, 
                       pf->cats_to_do_str != NULL ? cats_to_do : NULL, 
                       NULL, NULL, -1, subst_mod_is_codon_model(mod->subst_mod));
          /* (sufficient stats obtained only for categories of interest) */
      
          if (msa->length > 1000000) { /* throw out original data if
                                          very large */
            for (j = 0; j < msa->nseqs; j++) sfree(msa->seqs[j]);
            sfree(msa->seqs);
            msa->seqs = NULL;
          }
        }
        if (pf->random_init) 
          params = tm_params_init_random(mod);
        else if (input_mod != NULL) 
          params = tm_params_new_init_from_model(input_mod);
	else 
          params = tm_params_init(mod, .1, 5, pf->alpha);     

	if (pf->init_parsimony)
	  tm_params_init_branchlens_parsimony(params, mod, msa, cat);

        if (input_mod != NULL && mod->backgd_freqs != NULL && !pf->no_freqs && pf->init_backgd_from_data) {
          /* in some cases, the eq freqs are needed for
             initialization, but now they should be re-estimated --
             UNLESS user specifies --no-freqs */
	  vec_free(mod->backgd_freqs);
	  mod->backgd_freqs = NULL;
        }


        if (i == 0) {
          if (!quiet) fprintf(stderr, "Compacting sufficient statistics ...\n");
          ss_collapse_missing(msa, !pf->gaps_as_bases);
                                /* reduce number of tuples as much as
                                   possible */
        }

        if (!quiet) {
          fprintf(stderr, "Fitting tree model to %s using %s%s ...\n",
                  tmpstr->chars, tm_get_subst_mod_string(subst_mod),
                  mod->nratecats > 1 ? " (with rate variation)" : "");
        }

        if (pf->use_em)
          tm_fit_em(mod, msa, params, cat, pf->precision, pf->max_em_its, pf->logf, error_file);
        else
          tm_fit(mod, msa, params, cat, pf->precision, pf->logf, pf->quiet, error_file);
      }

      if (pf->output_fname_root != NULL) 
	str_cpy_charstr(mod_fname, pf->output_fname_root);
      else str_clear(mod_fname);
      if (pf->window_coords != NULL) {
	if (mod_fname->length != 0) 
	  str_append_char(mod_fname, '.');
        str_append_charstr(mod_fname, "win-");
        str_append_int(mod_fname, win/2 + 1);
      }
      if (cat != -1 && pf->nonoverlapping == FALSE) {
	if (mod_fname->length != 0) 
	  str_append_char(mod_fname, '.');
        if (pf->cm != NULL)  
          str_append(mod_fname, cm_get_feature_unique(pf->cm, cat));
        else 
          str_append_int(mod_fname, cat);
      }
      if (pf->output_fname_root != NULL)
	str_append_charstr(mod_fname, ".mod");

      if (pf->output_fname_root != NULL) {
	if (!quiet) fprintf(stderr, "Writing model to %s ...\n", 
			    mod_fname->chars);
	F = phast_fopen(mod_fname->chars, "w+");
	tm_print(F, mod);
	phast_fclose(F);
      }
      if (pf->results != NULL)
	lol_push_treeModel(pf->results, mod, mod_fname->chars);

      /* output posterior probabilities, if necessary */
      if (pf->do_bases || pf->do_expected_nsubst || 
	  pf->do_expected_nsubst_tot || pf->do_expected_nsubst_col) {
	print_post_prob_stats(mod, msa, pf->output_fname_root, 
			      pf->do_bases, pf->do_expected_nsubst, 
			      pf->do_expected_nsubst_tot, 
			      pf->do_expected_nsubst_col, 0,
			      cat, quiet, NULL);
      }

      /* print window summary, if window mode */
      if (pf->window_coords != NULL) {
	int i, j, total=0;
	char c;
	if (gc == NULL)
	  gc = smalloc(msa->nseqs*sizeof(double));
	for (i=0; i < msa->nseqs; i++) {
	  total=0;
	  gc[i]=0;
	  for (j=0; j<msa->length; j++) {
	    c = msa_get_char(msa, i, j);
	    if ((!msa->is_missing[(int)c]) && c != GAP_CHAR) {
	      total++;
	      if (c=='C' || c=='G') gc[i]++;
	    }
	  }
	  gc[i] /= (double)total;
	}
        print_window_summary(WINDOWF, pf->window_coords, win, cat, mod, gc, 
                             ninf_sites, msa->nseqs, FALSE);
      }
      
      if (input_mod == NULL) tm_free(mod);
      if (params != NULL) vec_free(params);
    }
    if (pf->window_coords != NULL) 
      msa_free(msa);
  }
  if (error_file != NULL) phast_fclose(error_file);
  if (parsimony_cost_file != NULL) phast_fclose(parsimony_cost_file); 
  str_free(mod_fname);
  str_free(tmpstr);
  if (free_cm) {
    cm_free(pf->cm);
    pf->cm = NULL;
  }
  if (free_cats_to_do_str) {
    lst_free_strings(pf->cats_to_do_str);
    lst_free(pf->cats_to_do_str);
    pf->cats_to_do_str = NULL;
  }
  if (free_tree) 
    tr_free(tree);
  if (cats_to_do != NULL) lst_free(cats_to_do);
  if (free_window_coords) {
    lst_free(pf->window_coords);
    pf->window_coords = NULL;
  }
  if (gc != NULL)
    sfree(gc);
  return 0;
}
