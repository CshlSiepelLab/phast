/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <misc.h>
#include "pbsScoreMatrix.help"
#include <pbs_code.h>
#include <tree_model.h>

int main(int argc, char *argv[]) {
  char c;
  int opt_idx, bl, i, j, k, l, alph_size;
  PbsCode *code = NULL;
  TreeModel *mod;
  List *branchlens = lst_new_dbl(50);
  MarkovMatrix *P;
  Matrix *S = NULL, *S2 = NULL, *expS = NULL;

  struct option long_opts[] = {
    {"branch-length", 1, 0, 't'},    
    {"half-pbs", 0, 0, 'H'},
    {"no-pbs", 0, 0, 'N'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  /* argument variables and defaults */
  enum {FULL, HALF, NONE} pbs_mode = FULL;

  while ((c = getopt_long(argc, argv, "a:b:h", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 't':
      lst_push_dbl(branchlens, get_arg_dbl_bounds(optarg, 0, INFTY));
      break;
    case 'H':
      pbs_mode = HALF;
      break;
    case 'N':
      pbs_mode = NONE;
      break;
    case 'h':
      printf(HELP);
      exit(0);
    case '?':
      die("Bad argument.  Try 'pbsScoreMatrix -h'.\n");
    }
  }

  if (pbs_mode != NONE && optind != argc - 2) 
    die("Two input filenames required.  Try 'pbsScoreMatrix -h'.\n");
  else if (pbs_mode == NONE && optind != argc - 1) 
    die("One input filenames required.  Try 'pbsScoreMatrix -h'.\n");
  
  set_seed(-1);
    
  mod = tm_new_from_file(fopen_fname(argv[optind], "r"));
  
  if (pbs_mode != NONE)
    code = pbs_new_from_file(fopen_fname(argv[optind+1], "r"));

  alph_size = mod->rate_matrix->size;
  P = mm_new(alph_size, NULL, DISCRETE);

  if (pbs_mode == FULL) {
    expS = mat_new(alph_size, alph_size);
    S2 = mat_new(code->code_size, code->code_size);
    mat_zero(S2);
  }
  else if (pbs_mode == HALF) {
    expS = mat_new(alph_size, alph_size);
    S2 = mat_new(code->code_size, alph_size);
    mat_zero(S2);
  }
  else 
    S = mat_new(alph_size, alph_size);

  if (lst_size(branchlens) == 0) { /* no branch length specified; use
				      all in tree */
    for (i = 0; i < mod->tree->nnodes; i++) {
      TreeNode *n = lst_get_ptr(mod->tree->nodes, i);
      if (n->parent != NULL)
	lst_push_dbl(branchlens, n->dparent);
    }
  }

  /* print header */
  printf("# Score matrix/matrices generated by pbsScoreMatrix based on model \"%s\"", argv[optind]);
  if (pbs_mode != NONE)
    printf("\n# and code file \"%s\"", argv[optind+1]);
  printf("\n# Element (i, j) of matrix is log-odds score for aligning ith %s with jth %s\n", pbs_mode == NONE ? "base" : "code index", pbs_mode == FULL ? "code index" : "base");
  if (pbs_mode == NONE) {
    printf("# Bases are ordered as in model: ");
    for (i = 0; i < alph_size; i++) printf("%c ", mod->rate_matrix->states[i]);
    printf("\n");
  }

  for (bl = 0; bl < lst_size(branchlens); bl++) {
    tm_set_subst_matrix(mod, P, lst_get_dbl(branchlens, bl));

    /* set up S or expS, whichever is needed */
    for (i = 0; i < alph_size; i++) {
      for (j = 0; j < alph_size; j++) {
	double odds = mm_get(P, i, j) / vec_get(mod->backgd_freqs, j);
	if (pbs_mode == NONE)
	  mat_set(S, i, j, log2(odds));
	else 
	  mat_set(expS, i, j, odds); 
      }
    }

    /* now derive output matrix S2 from S */
    if (pbs_mode == FULL) {
      for (k = 0; k < code->code_size; k++) {
	for (l = 0; l < code->code_size; l++) {
	  for (i = 0; i < alph_size; i++) 
	    for (j = 0; j < alph_size; j++) 
	      S2->data[k][l] += vec_get(code->rp[k], i) * 
		vec_get(code->rp[l], j) * mat_get(expS, i, j);
	  S2->data[k][l] = log2(S2->data[k][l]);
	}
      }
    }
    else if (pbs_mode == HALF) {
      for (k = 0; k < code->code_size; k++) {
	for (j = 0; j < alph_size; j++) {
	  for (i = 0; i < alph_size; i++) 
	    S2->data[k][j] += vec_get(code->rp[k], i) * mat_get(expS, i, j);
	  S2->data[k][j] = log2(S2->data[k][j]);	  
	}
      }
    }
    else 
      S2 = S;
    
    /* now print */
    printf("\n# t=%f\n", lst_get_dbl(branchlens, bl));	     
    mat_print(S2, stdout);
  }
  
  return 0;
}
