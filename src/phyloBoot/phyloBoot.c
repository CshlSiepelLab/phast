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
#include <msa.h>
#include <maf.h>
#include <sufficient_stats.h>
#include <numerical_opt.h>
#include <tree_model.h>
#include <fit_em.h>
#include <time.h>
#include "phyloBoot.help"

/* attempt to provide a brief description of each estimated parameter,
   based on a given TreeModel definition */
void set_param_descriptions(char **descriptions, TreeModel *mod) {
  List *traversal;
  String *str = str_new(STR_MED_LEN);
  int nrv_params = tm_get_nratevarparams(mod);
  int nrm_params = tm_get_nratematparams(mod);
  int i, j, idx;
  char *tempstr = smalloc((mod->order+2)*sizeof(char));

  if (mod->estimate_branchlens != TM_BRANCHLENS_ALL)
    die("ERROR set_param_descriptions: mod->estimate_branchlens != TM_BRANCHLENS_ALL\n");
  if (mod->estimate_backgd)
    die("ERROR set_param_descriptions: mod->estimate_backgd is TRUE\n");
  if (mod->alt_subst_mods != NULL)
    die("ERROR set_param_descriptions: mod->alt_subst_mods is not NULL\n");

  /* Tree scale descriptions */
  if (mod->estimate_branchlens != TM_BRANCHLENS_ALL) {
    sprintf(descriptions[mod->scale_idx], "tree scale");
    if (mod->subtree_root != NULL)
      sprintf(descriptions[mod->scale_idx+1], "subtree scale");
  }

  /* Branch length descriptions */
  traversal = tr_preorder(mod->tree);
  idx=0;
  for (i = 0; i < lst_size(traversal); i++) {
    TreeNode *n = lst_get_ptr(traversal, i);
    if (n->parent == NULL) continue;
    /* if the model is reversible, then the first parameter is
       the sum of the lengths of the two branches from the root */
    if ((n == mod->tree->lchild ||
	 n == mod->tree->rchild) && tm_is_reversible(mod))
      sprintf(descriptions[mod->bl_idx+idx], "branch (spans root)"); 
    else {
      if (strlen(n->name) > 0)
        sprintf(descriptions[mod->bl_idx+idx], 
		"branch (lf_%s->anc_%d)", n->name, n->parent->id);
      else 
        sprintf(descriptions[mod->bl_idx+idx], 
		"branch (anc_%d->anc_%d)", n->id, n->parent->id);
    }
    idx++;
  }

  /* Background frequency descriptions */
  if (mod->backgd_idx >= 0) {
   for (i=0; i < mod->backgd_freqs->size; i++) {
      get_tuple_str(tempstr, i, mod->order+1, mod->rate_matrix->states);
      tempstr[mod->order+1] = '\0';
      sprintf(descriptions[mod->backgd_idx + i], "background freq(%s)", tempstr);
   } 
  }

  /* Rate variation descriptions */
  for (i = 0; i < nrv_params; i++) {
    if (nrv_params == 1) 
      strcpy(descriptions[mod->ratevar_idx+i], "alpha");
    else
      sprintf(descriptions[mod->ratevar_idx+i], "rate var #%d", i+1);
  }

  /* Rate matrix value descriptions */
  if (nrm_params == 1)
  {
    strcpy(descriptions[mod->ratematrix_idx], "kappa");
  }
  else {
    for (i = 0; i < nrm_params; i++) {
      List *rows = mod->rate_matrix_param_row[mod->ratematrix_idx+i];
      List *cols = mod->rate_matrix_param_col[mod->ratematrix_idx+i];
      str_cpy_charstr(str, "rmatrix");
      for (j = 0; j < lst_size(rows); j++) {
        char tmp[STR_SHORT_LEN];
        sprintf(tmp, " (%d,%d)", lst_get_int(rows, j) + 1, 
                lst_get_int(cols, j) + 1);
        str_append_charstr(str, tmp);
      }
      strcpy(descriptions[mod->ratematrix_idx+i], str->chars);
    }
  }
  str_free(str);
  free(tempstr);
}

int main(int argc, char *argv[]) {
  
  /* variables for args with default values */
  int default_nsites = 1000, nsites = -1, nreps = 100, input_format = FASTA, 
    subst_mod = REV, nrates = 1, precision = OPT_HIGH_PREC, dump_format = FASTA;
  int quiet = FALSE, use_em = FALSE, random_init = FALSE, parametric = FALSE, 
    do_estimates = TRUE;
  TreeNode *tree = NULL;
  TreeModel *model = NULL, *init_mod = NULL;
  MSA *msa=NULL;
  char *dump_mods_root = NULL, *dump_msas_root = NULL, *ave_model = NULL;
  TreeModel **input_mods = NULL;

  /* other variables */
  FILE *INF, *F;
  char c;
  int i, j, opt_idx, nparams = -1, seed = -1;
  String *tmpstr;
  List **estimates=NULL;
  double *p = NULL;
  int *tmpcounts=NULL;
  char **descriptions = NULL;
  List *tmpl;
  char fname[STR_MED_LEN];
  char tmpchstr[STR_MED_LEN];
  TreeModel *repmod = NULL;
  double subtreeScale=1.0, subtreeSwitchProb=0.0, scale=1.0;
  char *subtreeName=NULL, *scaleFileName=NULL;
  TreeModel *subtreeModel=NULL;
  List *scaleLst=NULL, *subtreeScaleLst=NULL, *nsitesLst=NULL;
  FILE *scaleFile;

  struct option long_opts[] = {
    {"nsites", 1, 0, 'L'},
    {"nreps", 1, 0, 'n'},
    {"msa-format", 1, 0, 'i'},
    {"dump-mods", 1, 0, 'd'},
    {"alignments-only", 1, 0, 'a'},
    {"dump-samples", 1, 0, 'm'},
    {"dump-format", 1, 0, 'o'},
    {"no-estimates", 0, 0, 'x'}, /* for backward compatibility */
    {"read-mods", 1, 0, 'R'},
    {"output-average", 1, 0, 'A'},
    {"quiet", 0, 0, 'q'},
    {"help", 0, 0, 'h'},
    {"tree", 1, 0, 't'},
    {"subst-mod", 1, 0, 's'},
    {"nrates", 1, 0, 'k'},
    {"EM", 0, 0, 'E'},
    {"precision", 1, 0, 'p'},
    {"init-model", 1, 0, 'M'},
    {"init-random", 0, 0, 'r'},
    {"subtree", 1, 0, 'S'},
    {"subtree-switch", 1, 0, 'w'},
    {"subtree-scale", 1, 0, 'l'},
    {"scale", 1, 0, 'P'},
    {"scale-file", 1, 0, 'F'},
    {"seed", 1, 0, 'D'},
    {0, 0, 0, 0}
  };
  
  while ((c = getopt_long(argc, argv, "L:n:i:d:a:m:o:xR:qht:s:k:Ep:M:S:w:l:P:F:D:r", 
                          long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'L':
      nsites = get_arg_int_bounds(optarg, 10, INFTY);
      break;
    case 'n':
      if (input_mods != NULL) die("ERROR: Can't use --nreps with --read-mods.\n");
      nreps = get_arg_int_bounds(optarg, 1, INFTY);
      break;
    case 'i':
      input_format = msa_str_to_format(optarg);
      if (input_format == -1)
        die("ERROR: unrecognized alignment format.  Type 'phyloBoot -h' for usage.\n");
      break;
    case 'd':
      dump_mods_root = optarg;
      break;
    case 'a':
      dump_msas_root = optarg;
      do_estimates = FALSE;
      break;
    case 'm':
      dump_msas_root = optarg;
      break;
    case 'o':
      dump_format = msa_str_to_format(optarg);
      if (dump_format == -1)
        die("ERROR: unrecognized dump format.  Type 'phyloBoot -h' for usage.\n");      
      break;
    case 'x':                   /* for backward compatibility */
      do_estimates = FALSE;
      break;
    case 'R':
      tmpl = get_arg_list(optarg);
      nreps = lst_size(tmpl);
      input_mods = smalloc(nreps * sizeof(void*));
      for (i = 0; i < lst_size(tmpl); i++) {
        FILE *F = phast_fopen(((String*)lst_get_ptr(tmpl, i))->chars, "r");
        input_mods[i] = tm_new_from_file(F, 1);
        phast_fclose(F);
      }
      lst_free_strings(tmpl); lst_free(tmpl);
      break;
    case 'A':
      ave_model = optarg;
      break;
    case 'q':
      quiet = 1;
      break;
    case 'P':
      scale=atof(optarg);
      break;
    case 'S':
      subtreeName=optarg;
      break;
    case 'F':
      scaleFileName = optarg;
      break;
    case 'w':
      subtreeSwitchProb=atof(optarg);
      if (subtreeSwitchProb > 1.0 || subtreeSwitchProb < 0.0) 
	die("ERROR: --subtree-switch argument should be between 0 and 1\n");
      break;
    case 'l':
      subtreeScale = atof(optarg);
      if (subtreeScale < 0.0)
	die("ERROR: --subtree-scale argument should be >=0.0\n");
      break;
    case 'h':
      printf("%s", HELP);
      exit(0);
    case 't':
      if (optarg[0] == '(')     /* in this case, assume topology given
                                   at command line */
        tree = tr_new_from_string(optarg);
      else 
        tree = tr_new_from_file(phast_fopen(optarg, "r"));
      break;
    case 's':
      subst_mod = tm_get_subst_mod_type(optarg);
      if (subst_mod == UNDEF_MOD) 
        die("ERROR: illegal substitution model.  Type \"phyloBoot -h\" for usage.\n");
      break;
    case 'k':
      nrates = get_arg_int_bounds(optarg, 1, INFTY);
      break;
    case 'E':
      use_em = 1;
      break;
    case 'p':
      if (!strcmp(optarg, "LOW")) precision = OPT_LOW_PREC;
      else if (!strcmp(optarg, "MED")) precision = OPT_MED_PREC;
      else if (!strcmp(optarg, "HIGH")) precision = OPT_HIGH_PREC;
      else die("ERROR: --precision must be LOW, MED, or HIGH.\n\n");
      break;
    case 'M':
      init_mod = tm_new_from_file(phast_fopen(optarg, "r"), 1);
      break;
    case 'r':
      random_init = 1;
      break;
    case 'D':
      seed = get_arg_int_bounds(optarg, 1, INFTY);
      break;
    case '?':
      die("Bad argument.  Try '%s -h'.\n", argv[0]);
    }
  }

  set_seed(seed);

  if ((subtreeScale!=1.0 || subtreeSwitchProb!=0.0) &&
      subtreeName==NULL)
    die("ERROR: need to use --subtree with --subtree-scale or --subtree-switch\n");

  if (input_mods == NULL) {     /* only do if models aren't given */
    if (optind != argc - 1) 
      die("Input filename required.  Try '%s -h'.\n", argv[0]);

    INF = phast_fopen(argv[optind], "r");

    tmpstr = str_new_charstr(argv[optind]);
    if (str_ends_with_charstr(tmpstr, ".mod")) {
      parametric = TRUE;
      model = tm_new_from_file(INF, 1);
      tree = model->tree;
      if (scale != 1.0)
	tm_scale_branchlens(model, scale, 0);
      if (subtreeName != NULL) {
	subtreeModel = tm_create_copy(model);
	tm_scale_branchlens(subtreeModel, subtreeScale, 0);
      }
    }
    else {
      if (subtreeName != NULL) 
	die("subtree option only works if .mod file given (parametric mode)\n");
      if (input_format == -1)
	input_format = msa_format_for_content(INF, 1);
      if (input_format == MAF)
        msa = maf_read(INF, NULL, 1, NULL, NULL, NULL, -1, FALSE, NULL, NO_STRIP, FALSE);
      else
        msa = msa_new_from_file_define_format(INF, input_format, NULL);

      /* represent as SS and get rid of seqs (important below) */
      if (msa->ss == NULL) {
        ss_from_msas(msa, tm_order(subst_mod) + 1, FALSE, NULL, NULL, NULL, -1,
		     subst_mod_is_codon_model(subst_mod));
        for (i = 0; i < msa->nseqs; i++)
          sfree(msa->seqs[i]);
        sfree(msa->seqs);
        msa->seqs = NULL;
      }
    }

    /* general set up -- different for parametric and non-parametric cases */
    if (!parametric) {
      if  (scaleFileName != NULL) die("ERROR: --scale-file only works in parametric mode\n");
      if (tree == NULL) {
        if (msa->nseqs == 2) {
          sprintf(tmpchstr, "(%s,%s)", msa->names[0], msa->names[1]);
          tree = tr_new_from_string(tmpchstr);
        }
        else if (msa->nseqs == 3 && subst_mod_is_reversible(subst_mod)) {
          sprintf(tmpchstr, "(%s,(%s,%s))", msa->names[0], msa->names[1], 
                  msa->names[2]);
          tree = tr_new_from_string(tmpchstr);
        }
        else if (do_estimates)
          die("ERROR: must specify tree topology.\n");
      }
      else if (msa->nseqs * 2 - 1 != tree->nnodes)
        die("ERROR: Tree must have 2n-1 nodes, where n is the number of sequences in the\nalignment.  Even with a reversible model, specify a rooted tree; the root\nwill be ignored in the optimization procedure.\n");

      if (msa_alph_has_lowercase(msa)) msa_toupper(msa); 
      msa_remove_N_from_alph(msa); /* for backward compatibility */

      if (nsites == -1) nsites = msa->length;

      /* define probability vector from tuple counts */
      if (msa->ss == NULL)
        ss_from_msas(msa, tm_order(subst_mod) + 1, FALSE, NULL, NULL, NULL, -1,
		     subst_mod_is_codon_model(subst_mod));
      p = smalloc(msa->ss->ntuples * sizeof(double));
      for (i = 0; i < msa->ss->ntuples; i++) p[i] = msa->ss->counts[i];
      normalize_probs(p, msa->ss->ntuples);
      tmpcounts = smalloc(msa->ss->ntuples * sizeof(int));
    }
    else {                        /* parametric */
      if (scaleFileName != NULL) {
	double tempScale, tempSubtreeScale;
	int tempNumSite;
	if (subtreeSwitchProb!=0.0) die("ERROR: Cannot use --subtree-switch with --scale-file (not implemented)\n");
	scaleFile = phast_fopen(scaleFileName, "r");
	scaleLst = lst_new_dbl(100);
	subtreeScaleLst = lst_new_dbl(100);
	nsitesLst = lst_new_int(100);
	if (nsites != -1) fprintf(stderr, "Warning: number of simulated sites will be determined by --scale-file.  Ignoring --nsites argument\n");
	if (subtreeScale!=1.0) fprintf(stderr, "Warning: --subtree-file overrides --subtree-scale.  Ignoring --subtree-scale argument\n");
	subtreeScale = 1.0;
	nsites=0;
	while (EOF != fscanf(scaleFile, "%i %lf %lf", 
			     &tempNumSite, &tempScale, &tempSubtreeScale)) {
	  lst_push_int(nsitesLst, tempNumSite);
	  lst_push_dbl(scaleLst, tempScale);
	  lst_push_dbl(subtreeScaleLst, tempSubtreeScale);
	  nsites += tempNumSite;
	}
	phast_fclose(scaleFile);
      }
      if (nsites == -1) nsites = default_nsites;
    }
  } /* if input_mods == NULL */

  for (i = 0; i < nreps; i++) {
    Vector *params=NULL;
    TreeModel *thismod=NULL;

    /* generate alignment */
    if (input_mods == NULL) {   /* skip if models given */
      if (parametric) {
	if (scaleLst != NULL)
	  msa = tm_generate_msa_scaleLst(nsitesLst, scaleLst, subtreeScaleLst,
					 model, subtreeName);
	else if (subtreeName!=NULL && (subtreeScale!=1.0 || subtreeSwitchProb!=0.0)) 
	  msa = tm_generate_msa_random_subtree(nsites, model, subtreeModel, 
					       subtreeName, subtreeSwitchProb);
	else msa = tm_generate_msa(nsites, NULL, &model, NULL);
      }
      else {

        double sum=0;
        for (j = 0; j < msa->ss->ntuples; j++) sum+=p[j];

        mn_draw(nsites, p, msa->ss->ntuples, tmpcounts);
                                /* here we simply redraw numbers of
                                   tuples from multinomial distribution
                                   defined by orig alignment */
        for (j = 0; j < msa->ss->ntuples; j++) msa->ss->counts[j] = tmpcounts[j];
                                /* (have to convert from int to double) */
        msa->length = nsites;
      }

      if (dump_msas_root != NULL) {
        sprintf(fname, "%s.%d.%s", dump_msas_root, i+1, 
                msa_suffix_for_format(dump_format));
        if (!quiet) fprintf(stderr, "Dumping alignment to %s...\n", fname);
        F = phast_fopen(fname, "w+");

        if (dump_format == SS) { /* output ss */
          if (msa->ss == NULL)   /* (only happens in parametric case) */
            ss_from_msas(msa, tm_order(subst_mod) + 1, FALSE, NULL, NULL, NULL, -1, subst_mod_is_codon_model(subst_mod));
          ss_write(msa, F, FALSE);
        }
        else {                  /* output actual seqs */
          if (!parametric) {   /* only have SS; need to create seqs */
            ss_to_msa(msa);            
            msa_permute(msa);
          }
          msa_print(F, msa, dump_format, FALSE);
          if (!parametric) {   /* need to get rid of seqs because msa
                                   object reused */
            for (j = 0; j < msa->nseqs; j++) sfree(msa->seqs[j]);
            sfree(msa->seqs);
            msa->seqs = NULL;
          }
        }
        phast_fclose(F);
      }
    }

    /* now estimate model parameters */
    if (input_mods == NULL && do_estimates) {
      if (init_mod == NULL) 
        thismod = tm_new(tr_create_copy(tree), NULL, NULL, subst_mod, 
                         msa->alphabet, nrates, 1, NULL, -1);
      else {
        thismod = tm_create_copy(init_mod);  
        tm_reinit(thismod, subst_mod, nrates, thismod->alpha, NULL, NULL);
      }

      if (random_init) 
        params = tm_params_init_random(thismod);
      else if (init_mod != NULL)
        params = tm_params_new_init_from_model(init_mod);
      else
        params = tm_params_init(thismod, .1, 5, 1);    

      if (init_mod != NULL && thismod->backgd_freqs != NULL) {
        vec_free(thismod->backgd_freqs);
        thismod->backgd_freqs = NULL; /* force re-estimation */
      }

      if (!quiet) 
        fprintf(stderr, "Estimating model for replicate %d of %d...\n", i+1, nreps);

      if (use_em)
        tm_fit_em(thismod, msa, params, -1, precision, -1, NULL, NULL);
      else
        tm_fit(thismod, msa, params, -1, precision, NULL, quiet, NULL);

      if (dump_mods_root != NULL) {
        sprintf(fname, "%s.%d.mod", dump_mods_root, i+1);
        if (!quiet) fprintf(stderr, "Dumping model to %s...\n", fname);
        F = phast_fopen(fname, "w+");
        tm_print(F, thismod);
        phast_fclose(F);
      }
    } 

    else if (input_mods != NULL) { 
      /* in this case, we need to set up a parameter vector from
         the input model */
      thismod = input_mods[i];
      params = tm_params_new_init_from_model(thismod);
      if (nparams > 0 && params->size != nparams)
        die("ERROR: input models have different numbers of parameters.\n");
      if (repmod == NULL) repmod = thismod; /* keep around one representative model */
    }

    /* collect parameter estimates */
    if (do_estimates) {
      /* set up record of estimates; easiest to init here because number
         of parameters not always known above */
      if (nparams <= 0) {
        nparams = params->size;
        estimates = smalloc(nparams * sizeof(void*));
        descriptions = smalloc(nparams * sizeof(char*));
        for (j = 0; j < nparams; j++) {
          estimates[j] = lst_new_dbl(nreps);
          descriptions[j] = smalloc(STR_MED_LEN * sizeof(char));
	  descriptions[j][0] = '\0';
        }
        set_param_descriptions(descriptions, thismod);
      }

      /* record estimates for this replicate */
      for (j = 0; j < nparams; j++)
        lst_push_dbl(estimates[j], vec_get(params, j));
    }

    if (input_mods == NULL && do_estimates) {
      if (repmod == NULL) repmod = thismod; /* keep around one representative model */
      else tm_free(thismod);
    }
    if (do_estimates) vec_free(params);
    if (parametric) msa_free(msa);
  }

  /* finally, compute and print stats */
  if (do_estimates) {
    Vector *ave_params = vec_new(nparams); 
    printf("%-7s %-25s %9s %9s %9s %9s %9s %9s %9s %9s %9s\n", "param", 
           "description", "mean", "stdev", "median", "min", "max", "95%_min", 
           "95%_max", "90%_min", "90%_max");
    int param_num = 0;
    for (j = 0; j < nparams; j++) {
      if (strlen(descriptions[j]) < 1)
        continue;
      double mean = lst_dbl_mean(estimates[j]);
      double stdev = lst_dbl_stdev(estimates[j]);
      double quantiles[] = {0, 0.025, 0.05, 0.5, 0.95, 0.975, 1};
      double quantile_vals[7]; 
      lst_qsort_dbl(estimates[j], ASCENDING);
      lst_dbl_quantiles(estimates[j], quantiles, 7, quantile_vals);

      printf("%-7d %-25s %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n", 
             param_num, descriptions[j], mean, stdev, quantile_vals[3], quantile_vals[0], 
             quantile_vals[6], quantile_vals[1], quantile_vals[5], quantile_vals[2], 
             quantile_vals[4]);
      vec_set(ave_params, j, mean);
      param_num++;
    }

    if (ave_model != NULL) {
      for (i=0; i < repmod->all_params->size; i++) repmod->param_map[i]=i;
      tm_unpack_params(repmod, ave_params, -1);
      if (!quiet) fprintf(stderr, "Writing average model to %s...\n", ave_model);
      tm_print(phast_fopen(ave_model, "w+"), repmod);
    }
    vec_free(ave_params);
  }

  if (!quiet) fprintf(stderr, "Done.\n");
  return 0;
}
