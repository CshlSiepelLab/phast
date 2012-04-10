/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* phyloFit - fit phylogenetic model(s) to a multiple alignment */
 
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
#include "phyloFit.help"


int main(int argc, char *argv[]) {
  char *msa_fname = NULL, *alph = "ACGT";
  msa_format_type input_format = UNKNOWN_FORMAT;
  char c;
  int opt_idx, seed=-1;
  String *optstr;
  List *tmplist = NULL; 
  struct phyloFit_struct *pf;
  FILE *infile;
  
  struct option long_opts[] = {
    {"msa", 1, 0, 'm'},
    {"tree", 1, 0, 't'},
    {"subst-mod", 1, 0, 's'},
    {"msa-format", 1, 0, 'i'},
    {"nrates", 1, 0, 'k'},
    {"alpha", 1, 0, 'a'},
    {"features", 1, 0, 'g'},
    {"catmap", 1, 0, 'c'},
    {"log", 1, 0, 'l'},
    {"out-root", 1, 0, 'o'},
    {"EM", 0, 0, 'E'},
    {"error", 1, 0, 'e'},
    {"precision", 1, 0, 'p'},
    {"do-cats", 1, 0, 'C'},
    {"non-overlapping", 0, 0, 'V'},
    {"markov", 0, 0, 'N'},
    {"reverse-groups", 1, 0, 'R'},
    {"init-model", 1, 0, 'M'},
    {"init-random", 0, 0, 'r'},
    {"init-parsimony", 0, 0, 'y'},
    {"print-parsimony", 1, 0, 'Y'},
    {"lnl", 0, 0, 'L'},
    {"scale-only", 0, 0, 'B'},
    {"scale-subtree", 1, 0, 'S'},
    {"estimate-freqs", 0, 0, 'F'},
    {"sym-freqs", 0, 0, 'W'},
    {"no-freqs", 0, 0, 'f'},
    {"no-rates", 0, 0, 'n'},
    {"no-opt", 1, 0, 'O'},
    {"min-informative", 1, 0, 'I'},
    {"gaps-as-bases", 0, 0, 'G'},     
    {"quiet", 0, 0, 'q'},
    {"help", 0, 0, 'h'},
    {"windows", 1, 0, 'w'},
    {"windows-explicit", 1, 0, 'v'},
    {"ancestor", 1, 0, 'A'},
    {"post-probs", 0, 0, 'P'},
    {"expected-subs", 0, 0, 'X'},
    {"expected-total-subs", 0, 0, 'Z'},
    {"expected-subs-col", 0, 0, 'J'},
    {"column-probs", 0, 0, 'U'},
    {"rate-constants", 1, 0, 'K'},
    {"ignore-branches", 1, 0, 'b'},
    {"clock", 0, 0, 'z'},
    {"alt-model", 1, 0, 'd'},
    {"label-branches", 1, 0, 0},
    {"label-subtree", 1, 0, 0},
    {"selection", 1, 0, 0},
    {"bound", 1, 0, 'u'},
    {"seed", 1, 0, 'D'},
    {0, 0, 0, 0}
  };

  // NOTE: remaining shortcuts left: HjQx

  pf = phyloFit_struct_new(0);

  while ((c = getopt_long(argc, argv, "m:t:s:g:c:C:i:o:k:a:l:w:v:M:p:A:I:K:S:b:d:O:u:Y:e:D:GVENRqLPXZUBFfnrzhWyJ", long_opts, &opt_idx)) != -1) {
    switch(c) {
    case 'm':
      msa_fname = optarg;
      break;
    case 't':
      if (optarg[0] == '(')        /* in this case, assume topology given
                                   at command line */
        pf->tree = tr_new_from_string(optarg);
      else 
        pf->tree = tr_new_from_file(phast_fopen(optarg, "r"));
      break;
    case 's':
      pf->subst_mod = tm_get_subst_mod_type(optarg);
      if (pf->subst_mod == UNDEF_MOD) 
        die("ERROR: illegal substitution model.     Type \"phyloFit -h\" for usage.\n");
      break;
    case 'g':
      pf->gff = gff_read_set(phast_fopen(optarg, "r"));
      break;
    case 'c':
      pf->cm = cm_new_string_or_file(optarg);
      break;
    case 'C':
      pf->cats_to_do_str = get_arg_list(optarg);
      break;
    case 'V':
      pf->nonoverlapping = TRUE;
      break;
    case 'o':
      pf->output_fname_root = optarg;
      break;
    case 'k':
      pf->nratecats = get_arg_int_bounds(optarg, 0, INFTY);
      break;
    case 'a':
      pf->alpha = get_arg_dbl(optarg);
      break;
    case 'R':
      pf->reverse_group_tag = optarg;
      break;
    case 'i':
      input_format = msa_str_to_format(optarg);
      if (input_format == UNKNOWN_FORMAT)
        die("ERROR: unrecognized alignment format.    Type 'phyloFit -h' for usage.\n");
      break;
    case 'l':
      if (!strcmp(optarg, "-"))
	pf->logf = stderr;
      else pf->logf = phast_fopen(optarg, "w+");
      break;
    case 'N':
      pf->use_conditionals = 1;
      break;
    case 'w':
      tmplist = get_arg_list(optarg);
      if (lst_size(tmplist) != 2 ||
          str_as_int(lst_get_ptr(tmplist, 0), &(pf->window_size)) != 0 ||
          str_as_int(lst_get_ptr(tmplist, 1), &(pf->window_shift)) != 0) 
        die("ERROR: illegal arguments to --windows.\n");
      lst_free_strings(tmplist);
      lst_free(tmplist);
      break;
    case 'v':
      tmplist = get_arg_list(optarg);
      if (lst_size(tmplist) % 2 != 0) 
        die("ERROR: argument to --windows-explicit must be a list of even length.\n");
      pf->window_coords = str_list_as_int(tmplist);
      lst_free(tmplist);
      break;
    case 'E':
      pf->use_em = TRUE;
      break;
    case 'e':
      pf->error_fname=optarg;
      break;
    case 'p':
      if (!strcmp(optarg, "LOW")) pf->precision = OPT_LOW_PREC;
      else if (!strcmp(optarg, "MED")) pf->precision = OPT_MED_PREC;
      else if (!strcmp(optarg, "HIGH")) pf->precision = OPT_HIGH_PREC;
      else if (!strcmp(optarg, "VERY_HIGH")) pf->precision = OPT_VERY_HIGH_PREC;
      else die("ERROR: --precision must be LOW, MED, or HIGH.\n\n");
      break;
    case 'M':
      pf->input_mod = tm_new_from_file(phast_fopen(optarg, "r"), 1);
      break;
    case 'r':
      pf->random_init = TRUE;
      break;
    case 'y':
      pf->init_parsimony = TRUE;
      break;
    case 'Y':
      pf->init_parsimony = TRUE;
      pf->parsimony_cost_fname = optarg;
      pf->parsimony_only=TRUE;
      break; 
    case 'L':
      pf->likelihood_only = TRUE;
      break;
    case 'q':
      pf->quiet = TRUE;
      break;
    case 'P':
      pf->do_bases = TRUE;
      break;
    case 'X':
      pf->do_expected_nsubst = TRUE;
      break;
    case 'Z':
      pf->do_expected_nsubst_tot = TRUE;
      break;
    case 'J':
      pf->do_expected_nsubst_col = TRUE;
      break;
    case 'U':
      pf->likelihood_only = TRUE;        /* force -L */
      pf->nsites_threshold = 0;        /* also force this; typical use is
                                   with small number of tuples, no
                                   tuple_idx */
      pf->do_column_probs = TRUE;
      break;
    case 'A':
      pf->root_seqname = optarg;
      break;
    case 'I':
      pf->nsites_threshold = get_arg_int(optarg);
      break;
    case 'G':
      pf->gaps_as_bases = TRUE;
      alph = "ACGT-";
      break;
    case 'B':
      pf->estimate_scale_only = TRUE;
      break;
    case 'S':
      pf->subtree_name = optarg;
      break;       
    case 'F':
      pf->estimate_backgd = TRUE;
      break;
    case 'W':
      pf->estimate_backgd = TRUE;
      pf->symfreq = TRUE;
      break;
    case 'f':
      pf->no_freqs = TRUE;
      break;
    case 'n':
      pf->no_rates = TRUE;
      break;
    case 'K':
      tmplist = get_arg_list(optarg);
      pf->rate_consts = str_list_as_dbl(tmplist);
      pf->nratecats = lst_size(pf->rate_consts);
      pf->use_em = 1;
      lst_free_strings(tmplist); lst_free(tmplist);
      break;
    case 'b':
      pf->ignore_branches = get_arg_list(optarg);
      break;
    case 'z':
      pf->assume_clock = TRUE;
      break;
    case 'O':
      if (pf->nooptstr == NULL) 
	pf->nooptstr = str_new_charstr(optarg);
      else die("ERROR: no-opt argument can only be used once!  parameters can be comma-separated list.");
      break;
    case 'd':
      if (pf->alt_mod_str == NULL) {
	pf->alt_mod_str = lst_new_ptr(1);
      }
      optstr = str_new_charstr(optarg);
      lst_push_ptr(pf->alt_mod_str, optstr);
      break;
    case 0:
      if (strcmp(long_opts[opt_idx].name, "label-branches") == 0 ||
	  strcmp(long_opts[opt_idx].name, "label-subtree") == 0) {
	optstr = str_new_charstr(optarg);
	if (pf->label_str == NULL) {
	  pf->label_str = lst_new_ptr(3);
	  pf->label_type = lst_new_int(3);
	}
	lst_push_ptr(pf->label_str, optstr);
	lst_push_int(pf->label_type, 
		     strcmp(long_opts[opt_idx].name, "label-branches") == 0 ? 
		     BRANCH_TYPE : SUBTREE_TYPE);
      }
      else if (strcmp(long_opts[opt_idx].name, "selection") == 0) {
	pf->selection = get_arg_dbl(optarg);
	pf->use_selection = TRUE;
      }
      else {
	die("ERROR: unknown option.  Type 'phyloFit -h' for usage.\n");
      }
      break;
    case 'u':
      if (pf->bound_arg == NULL) 
	pf->bound_arg = lst_new_ptr(1);
      optstr = str_new_charstr(optarg);
      lst_push_ptr(pf->bound_arg, optstr);
      break;
    case 'D':
      seed = get_arg_int_bounds(optarg, 1, INFTY);
      break;
    case 'h':
      printf("%s", HELP);
      exit(0);
    case '?':
      die("ERROR: illegal argument.     Type 'phyloFit -h' for usage.\n");
    }
  }

  set_seed(seed);

  if (msa_fname == NULL) {
    if (optind >= argc) 
      die("ERROR: missing alignment filename.  Type 'phyloFit -h' for usage.\n");
    msa_fname = argv[optind];
    pf->msa_fname = msa_fname;
  }

  infile = phast_fopen(msa_fname, "r");

  if (input_format == UNKNOWN_FORMAT)
    input_format = msa_format_for_content(infile, 1);

  if (pf->nonoverlapping && (pf->use_conditionals || pf->gff != NULL || 
			     pf->cats_to_do_str || input_format == SS))
    die("ERROR: cannot use --non-overlapping with --markov, --features,\n--msa-format SS, or --do-cats.\n");


  /* read alignment */
  if (!pf->quiet) fprintf(stderr, "Reading alignment from %s ...\n", msa_fname);
  if (input_format == MAF) {
    pf->msa = maf_read(infile, NULL, 
		       tm_order(pf->subst_mod) + 1, 
		       NULL, pf->gff, pf->cm, 
		       pf->nonoverlapping ? tm_order(pf->subst_mod) + 1 : -1, 
		       FALSE, pf->reverse_group_tag, NO_STRIP, FALSE);
    if (pf->gaps_as_bases) 
      msa_reset_alphabet(pf->msa, alph);
  }
  else 
    pf->msa = msa_new_from_file_define_format(infile, 
				input_format, alph);

  /* set up for categories */
  /* first label sites, if necessary */
  pf->label_categories = (input_format != MAF);

  run_phyloFit(pf);

  if (pf->logf != NULL && pf->logf != stderr && pf->logf != stdout)
    phast_fclose(pf->logf);
  if (!pf->quiet) fprintf(stderr, "Done.\n");
  sfree(pf);
  
  return 0;
}
