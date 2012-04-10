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
#include <phylo_hmm.h>
#include <em.h>
#include <sufficient_stats.h>
#include <bed.h>
#include <dgamma.h>
#include <tree_likelihoods.h>
#include <maf.h>
#include "phast_cons.h"
#include "phastCons.help"


int main(int argc, char *argv[]) {
  struct phastCons_struct *p = phastCons_struct_new(0);
  struct option long_opts[] = {
    {"states", 1, 0, 'S'},
    {"hmm", 1, 0, 'H'},
    {"viterbi", 1, 0, 'V'},
    {"most-conserved", 1, 0, 'V'}, /* same as --viterbi */
    {"no-post-probs", 0, 0, 'n'},
    {"msa-format", 1, 0, 'i'},
    {"FC", 0, 0, 'X'},
    {"lambda", 1, 0, 'l'},
    {"target-coverage", 1, 0, 'C'},
    {"transitions", 1, 0, 't'},
    {"expected-length", 1, 0, 'E'},
    {"expected-lengths", 1, 0, 'E'}, /* for backward compatibility */
    {"estimate-trees", 1, 0, 'T'},
    {"estimate-rho", 1, 0, 'O'},
    {"rho", 1, 0, 'R'},
    {"gc", 1, 0, 'G'},
    {"ignore-missing", 0, 0, 'z'},
    {"nrates", 1, 0, 'k'},
    {"log", 1, 0, 'g'},
    {"refidx", 1, 0, 'r'},
    {"suppress-missing", 0, 0, 'x'}, /* for backward compatibility */
    {"reflect-strand", 1, 0, 'U'},
    {"catmap", 1, 0, 'c'},
    {"extrapolate", 1, 0, 'e'},
    {"indels", 0, 0, 'I'},
    {"max-micro-indel", 1, 0, 'Y'},
    {"indel-params", 1, 0, 'D'},
    {"min-informative-types", 1, 0, 'M'}, /* for backward compatibility */
    {"require-informative", 1, 0, 'M'},
    {"not-informative", 1, 0, 'F'},
    {"lnl", 1, 0, 'L'},
    {"seqname", 1, 0, 'N'},
    {"idpref", 1, 0, 'P'},
    {"score", 0, 0, 's'},
    {"coding-potential", 0, 0, 'p'},
    {"indels-only", 0, 0, 'J'},
    {"alias", 1, 0, 'A'},
    {"quiet", 0, 0, 'q'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  /* other vars */
  FILE *infile;
  char *msa_fname;
  char c;
  int opt_idx, i, coding_potential=FALSE;
  List *tmpl = NULL;
  String *tmpstr;
  char *mods_fname = NULL;
  List *mod_fname_list;
  msa_format_type msa_format = UNKNOWN_FORMAT;

  while ((c = getopt_long(argc, argv, 
			  "S:H:V:ni:k:l:C:G:zt:E:R:T:O:r:xL:sN:P:g:U:c:e:IY:D:JM:F:pA:Xqh", 
                          long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'S':
      p->states = get_arg_list(optarg);
      break;
    case 'H':
      p->hmm = hmm_new_from_file(phast_fopen(optarg, "r"));
      p->two_state = FALSE;
      break;
    case 'V':
      p->viterbi_f = phast_fopen(optarg, "w+");
      tmpstr = str_new_charstr(optarg);
      if (str_ends_with_charstr(tmpstr, ".gff")) 
	p->gff = TRUE;
      str_free(tmpstr);
      break;
    case 'n':
      p->post_probs = FALSE;
      break;
    case 'i':
      msa_format = msa_str_to_format(optarg);
      if (msa_format == UNKNOWN_FORMAT) 
	die("ERROR: bad argument to --msa-format\n");
      break;
    case 'X':
      p->FC = TRUE;
      p->two_state = FALSE;
      break;
    case 'l':
      if (optarg[0] != '~') 
	p->estim_lambda = FALSE;
      else optarg = &optarg[1];
      p->lambda = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'C':
      p->gamma = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'G':
      p->gc = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 't':
      p->set_transitions = TRUE;
      if (optarg[0] != '~') 
	p->estim_transitions = FALSE;
      else optarg = &optarg[1];
      tmpl = get_arg_list_dbl(optarg);
      if (lst_size(tmpl) != 2) 
        die("ERROR: bad argument to --transitions.\n");
      p->mu = lst_get_dbl(tmpl, 0);
      p->nu = lst_get_dbl(tmpl, 1);
      if (p->mu <= 0 || p->mu >= 1 || p->nu <= 0 || p->nu >= 1)
        die("ERROR: bad argument to --transitions.\n");
      lst_free(tmpl);
      break;
    case 'E':
      if (optarg[0] != '~') 
	p->estim_transitions = FALSE;
      else optarg = &optarg[1];
      p->omega = get_arg_dbl_bounds(optarg, 1, INFTY);
      p->mu = 1/p->omega;
      break;
    case 'T':
      p->estim_trees = TRUE;
      p->estim_trees_fname_root = optarg;
      break;
    case 'O':
      p->estim_rho = TRUE;
      p->estim_trees_fname_root = optarg;
      break;
    case 'z':
      p->ignore_missing = TRUE;
      break;
    case 'k':
      tmpl = get_arg_list_int(optarg);
      if (lst_size(tmpl) > 2) 
        die("ERROR: too many arguments with --nrates.\n");
      p->nrates = lst_get_int(tmpl, 0);
      if (p->nrates <= 0) 
        die("ERROR: bad argument to --nrates (%d).\n", p->nrates);
      if (lst_size(tmpl) == 2) {
        p->nrates2 = lst_get_int(tmpl, 1);
        if (p->nrates2 <= 0) 
          die("ERROR: bad argument to --nrates (%d).\n", p->nrates2);
      }
      lst_free(tmpl);
      break;
    case 'R':
      p->rho = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'g':
      if (!strcmp(optarg, "-")) 
	p->log_f = stderr;
      else p->log_f = phast_fopen(optarg, "w+");
      break;
    case 'r':
      p->refidx = get_arg_int_bounds(optarg, 0, INFTY);
      break;
    case 'x':
      /* do nothing; left in for backward compatibility */
      break;
    case 'U':
      p->pivot_states = get_arg_list(optarg); /* we want strings not ints
						 for phmm_new */
      break;
    case 'e':
      p->extrapolate_tree_fname = optarg;
      break;
    case 'I':
      p->indels = TRUE;
      break;
    case 'Y':
      p->max_micro_indel = get_arg_int_bounds(optarg, 1, INFTY);
      break;
    case 'D':
      if (optarg[0] != '~')
	p->estim_indels = FALSE;
      else optarg = &optarg[1];
      tmpl = get_arg_list_dbl(optarg);
      if (lst_size(tmpl) != 6) die("ERROR: bad argument to --indel-params.\n");
      p->alpha_0 = lst_get_dbl(tmpl, 0);
      p->beta_0 = lst_get_dbl(tmpl, 1);
      p->tau_0 = lst_get_dbl(tmpl, 2);
      p->alpha_1 = lst_get_dbl(tmpl, 3);
      p->beta_1 = lst_get_dbl(tmpl, 4);
      p->tau_1 = lst_get_dbl(tmpl, 5);
      if (p->alpha_0 < 0 || p->beta_0 < 0 || p->tau_0 < 0 || 
          p->alpha_1 < 0 || p->beta_1 < 0 || p->tau_1 < 0)
        die("ERROR: bad argument to --indel-params.\n");
      lst_free(tmpl);
      break;
    case 'J':
      p->indels_only = TRUE;
      p->two_state = FALSE;
      p->indels = TRUE;
      p->post_probs = FALSE;
      break;
    case 'M':
      p->inform_reqd = get_arg_list(optarg);
      break;
    case 'F':
      p->not_informative = get_arg_list(optarg);
      break;
    case 'c':
      p->cm = cm_new_string_or_file(optarg);
      break;
    case 'L':
      p->lnl_f = phast_fopen(optarg, "w+");
      break;
    case 'N':
      p->seqname = optarg;
      break;
    case 'P':
      p->idpref = optarg;
      break;
    case 's':
      p->score = TRUE;
      break;
    case 'p':
      coding_potential = TRUE;
      break;
    case 'A':
      p->alias_hash = make_name_hash(optarg);
      break;
    case 'q':
      p->results_f = NULL;
      break;
    case 'h':
      printf("%s", HELP);
      exit(0);
    case '?':
      die("Bad argument.  Try '%s -h'.\n", argv[0]);
    }
  }

  if ((!coding_potential && optind != argc - 2) ||
      (coding_potential && optind != argc - 2 && optind != argc - 1))
    die("ERROR: extra or missing arguments.  Try '%s -h'.\n", argv[0]);

  set_seed(-1);

  if (p->extrapolate_tree_fname != NULL &&
      !strcmp(p->extrapolate_tree_fname, "default")) {
    p->extrapolate_tree_fname = smalloc((strlen(PHAST_HOME)+100)*sizeof(char));
    #if defined(__MINGW32__)
      sprintf(p->extrapolate_tree_fname,
	      "%s\\data\\exoniphy\\mammals\\cftr25_hybrid.nh", PHAST_HOME);
    #else
      sprintf(p->extrapolate_tree_fname, 
              "%s/data/exoniphy/mammals/cftr25_hybrid.nh", PHAST_HOME);
    #endif
  }
  if (p->extrapolate_tree_fname != NULL)
    p->extrapolate_tree = tr_new_from_file(phast_fopen(p->extrapolate_tree_fname, "r"));

  mods_fname = (optind == argc - 2 ? argv[argc - 1] : NULL);
  /* if there are two args, mods are the second one; otherwise will
     use default mods for coding potential (see below) */
  
  /* set defaults for coding-potential mode */
  if (coding_potential) {
    char tmp[5000];
    p->two_state = FALSE;
    if (p->cm == NULL) 
      p->cm = cm_new_string_or_file("NCATS=4; CNS 1; CDS 2-4");
    if (p->hmm == NULL) {
      #if defined(__MINGW32__)
        sprintf(tmp, "%s\\data\\phastCons\\%s", PHAST_HOME,
                p->indels ? "simple-coding-indels.hmm" : "simple-coding.hmm");
      #else
        sprintf(tmp, "%s/data/phastCons/%s", PHAST_HOME,
                p->indels ? "simple-coding-indels.hmm" : "simple-coding.hmm");
      #endif
      if (p->results_f!=NULL) 
	fprintf(p->results_f, "Reading HMM from %s...\n", tmp);
      p->hmm = hmm_new_from_file(phast_fopen(tmp, "r"));
    }
    if (mods_fname == NULL) {
      #if defined(__MINGW32__)
        sprintf(tmp, "%s\\data\\exoniphy\\mammals\\r3.ncns.mod, %s\\data\\exoniphy\\mammals\\r3.cns.mod, %s\\data\\exoniphy\\mammals\\r3.cds-1.mod, %s\\data\\exoniphy\\mammals\\r3.cds-2.mod, %s\\data\\exoniphy\\mammals\\r3.cds-3.mod",  PHAST_HOME, PHAST_HOME, PHAST_HOME, PHAST_HOME, PHAST_HOME);
      #else
      sprintf(tmp, "\
%s/data/exoniphy/mammals/r3.ncns.mod,\
%s/data/exoniphy/mammals/r3.cns.mod,\
%s/data/exoniphy/mammals/r3.cds-1.mod,\
%s/data/exoniphy/mammals/r3.cds-2.mod,\
%s/data/exoniphy/mammals/r3.cds-3.mod", 
              PHAST_HOME, PHAST_HOME, PHAST_HOME, PHAST_HOME, PHAST_HOME);
      #endif
      mods_fname = tmp;
    }
    if (p->states == NULL) 
      p->states = get_arg_list("CDS");
    if (p->pivot_states == NULL) 
      p->pivot_states = get_arg_list("background,CNS");
  }

   /* read tree models */
  mod_fname_list = get_arg_list(mods_fname);
  p->nummod = lst_size(mod_fname_list);
  p->mod = (TreeModel**)smalloc(sizeof(TreeModel*) * p->nummod);
  for (i = 0; i < p->nummod; i++) {
    String *fname = lst_get_ptr(mod_fname_list, i);

    if (p->results_f != NULL)
      fprintf(p->results_f, "Reading tree model from %s...\n", fname->chars);
    p->mod[i] = tm_new_from_file(phast_fopen(fname->chars, "r"), 1);
    p->mod[i]->use_conditionals = 1;     
  }

  /* read alignment */
  msa_fname = argv[optind];
  infile = phast_fopen(msa_fname, "r");

  if (msa_format == UNKNOWN_FORMAT)
    msa_format = msa_format_for_content(infile, 1);
  if (p->results_f != NULL)
    fprintf(p->results_f, "Reading alignment from %s...\n", msa_fname);
  if (msa_format == MAF) {
    List *keepSeqs = tr_leaf_names(p->mod[0]->tree);
    p->msa = maf_read_cats_subset(infile, NULL, 1, NULL, NULL, 
				  NULL, -1, TRUE, NULL, NO_STRIP, FALSE, NULL, keepSeqs, 1);
    lst_free_strings(keepSeqs);
    lst_free(keepSeqs);
  }
  else
    p->msa = msa_new_from_file_define_format(infile, msa_format, NULL);

  /* use file name root for default seqname */
  if (p->viterbi_f != NULL && (p->seqname == NULL || p->idpref == NULL)) {
    String *tmp = str_new_charstr(msa_fname);
    if (!str_equals_charstr(tmp, "-")) {
      str_remove_path(tmp);
      str_root(tmp, '.');
      if (p->idpref == NULL) p->idpref = copy_charstr(tmp->chars);
      str_root(tmp, '.');         /* apply one more time for double suffix */
      if (p->seqname == NULL) p->seqname = tmp->chars;    
    }
  }

  phastCons(p);

  return 0;
}

