/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: exoniphy.c,v 1.44 2008-11-12 02:07:59 acs Exp $ */

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <phylo_hmm.h>
#include <gff.h>
#include <category_map.h>
#include <sufficient_stats.h>
#include <stringsplus.h>
#include <maf.h>
#include "exoniphy.help"

/* default background feature types; used when scoring predictions and
   reflecting HMM */
#define DEFAULT_BACKGD_TYPES "background,CNS"

/* default "cds" and "signal" feature tupes */
#define DEFAULT_CDS_TYPES "CDS"
#define DEFAULT_SIGNAL_TYPES "start_codon,stop_codon,5'splice,3'splice,prestart,cds5'ss,cds3'ss"

/* categories to be "absorbed" into CDS (want coords to be included in
   CDS) and categories to be "invisible" in output. For now, these are
   fixed; all of this should become simpler with generalized HMM
   architecture */
#define CDS_ABSORB_TYPES "start_codon,cds5'ss,cds3'ss"
#define INVISIBLE_TYPES "prestart,5'splice,3'splice,cds5'ss,cds3'ss"

/* parameters controlling evaluation of Sn/Sp tradeoff (see -Y option) */
#define SCALE_RANGE_MIN -20
#define SCALE_RANGE_MAX 10
#define NSENS_SPEC_TRIES 10

/* thresholds defining G+C ranges (-1 indicates end) */
double GC_THRESHOLDS[] = {0.40, 0.45, 0.50, 0.55, -1};

int main(int argc, char* argv[]) {

  /* variables for options, with defaults */
  int msa_format = UNKNOWN_FORMAT;
  int quiet = FALSE, reflect_hmm = FALSE, score = FALSE, indels = FALSE, 
    no_cns = FALSE;
  double bias = NEGINFTY;
  char *seqname = NULL, *grouptag = "transcript_id", *sens_spec_fname_root = NULL,
    *idpref = NULL, *extrapolate_tree_fname = NULL, *newname;
  List *model_fname_list = NULL, *no_gaps_str = NULL, *inform_reqd = NULL,
    *not_informative = NULL,
    *backgd_types = get_arg_list(DEFAULT_BACKGD_TYPES), 
    *cds_types = get_arg_list(DEFAULT_CDS_TYPES), 
    *signal_types = get_arg_list(DEFAULT_SIGNAL_TYPES),
    *cds_absorb_types = get_arg_list(CDS_ABSORB_TYPES),
    *invisible_types = get_arg_list(INVISIBLE_TYPES);
  TreeNode *extrapolate_tree = NULL;
  Hashtable *alias_hash = NULL;

  struct option long_opts[] = {
    {"hmm", 1, 0, 'H'},
    {"tree-models", 1, 0, 'm'},
    {"catmap", 1, 0, 'c'},
    {"msa-format", 1, 0, 'i'},
    {"data-path", 1, 0, 'D'}, 
    {"score", 0, 0, 'S'},
    {"seqname", 1, 0, 's'},
    {"idpref", 1, 0, 'p'},
    {"grouptag", 1, 0, 'g'},
    {"no-cns", 0, 0, 'x'},
    {"reflect-strand", 0, 0, 'U'},
    {"bias", 1, 0, 'b'},
    {"sens-spec", 1, 0, 'Y'},
    {"cds-types", 1, 0, 'C'},
    {"backgd-types", 1, 0, 'B'},
    {"signal-types", 1, 0, 'L'},
    {"indels", 0, 0, 'I'},
    {"no-gaps", 1, 0, 'W'},
    {"require-informative", 1, 0, 'N'},
    {"not-informative", 1, 0, 'n'},
    {"extrapolate", 1, 0, 'e'},
    {"alias", 1, 0, 'A'},
    {"quiet", 0, 0, 'q'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  /* other variables */
  FILE* F;
  FILE *infile;
  PhyloHmm *phmm;
  MSA *msa;
  TreeModel **mod;
  HMM *hmm = NULL;
  CategoryMap *cm = NULL;
  GFF_Set *predictions;
  String *data_path=NULL;
  char c;
  int i, j, ncats, trial, ntrials, opt_idx, gc_cat;
  double gc;
  char tmpstr[STR_LONG_LEN];
  char *msa_fname = NULL;
  String *fname_str = str_new(STR_LONG_LEN), *str;

  while ((c = getopt_long(argc, argv, "i:D:c:H:m:s:p:g:B:T:L:F:IW:N:n:b:e:A:xSYUhq", 
                          long_opts, &opt_idx)) != -1) {
    switch(c) {
    case 'i':
      msa_format = msa_str_to_format(optarg);
      if (msa_format == UNKNOWN_FORMAT) die("ERROR: bad alignment format.\n");
      break;
    case 'D':
      data_path = str_new_charstr(optarg);
      break;
    case 'c':
      cm = cm_new_string_or_file(optarg);
      break;
    case 'm':
      model_fname_list = get_arg_list(optarg);
      break;
    case 'H':
      hmm = hmm_new_from_file(phast_fopen(optarg, "r"));
      break;
    case 'x':
      no_cns = TRUE;
      break;
    case 'U':
      reflect_hmm = TRUE;
      break;
    case 'B':
      lst_free_strings(backgd_types); lst_free(backgd_types);
      backgd_types = get_arg_list(optarg);
      break;
    case 'T':
      lst_free_strings(cds_types); lst_free(cds_types); 
      cds_types = get_arg_list(optarg);
      break;
    case 'L':
      lst_free_strings(signal_types); lst_free(signal_types); 
      signal_types = get_arg_list(optarg);
      break;
    case 'S':
      score = TRUE;
      break;
    case 'b':
      bias = get_arg_dbl(optarg);
      break;
    case 'Y':
      sens_spec_fname_root = optarg;
      break;
    case 'W':
      no_gaps_str = get_arg_list(optarg);
      break;
    case 'N':
      inform_reqd = get_arg_list(optarg);
      break;
    case 'n':
      not_informative = get_arg_list(optarg);
      break;
    case 'I':
      indels = TRUE;
      break;
    case 's':
      seqname = optarg;
      break;
    case 'p':
      idpref = optarg;
      break;
    case 'g':
      grouptag = optarg;
      break;
    case 'e':
      extrapolate_tree_fname = optarg;
      break;
    case 'A':
      alias_hash = make_name_hash(optarg);
      break;
    case 'q':
      quiet = TRUE;
      break;
    case 'h':
      printf("%s", HELP);
      exit(0);
    case '?':
      die("ERROR: unrecognized option.  Try 'exoniphy -h' for help.\n");
    }
  }

  if (optind != argc - 1) {
      die("ERROR: alignment filename is required argument.  Try 'exoniphy -h' for help.\n");
  }

  set_seed(-1);

  if(data_path == NULL) {
    data_path = str_new(strlen(PHAST_HOME)+strlen("/data")+1);
    str_append_charstr(data_path, PHAST_HOME);
    #if defined(__MINGW32__)
      str_append_charstr(data_path, "\\data"); 
    #else
      str_append_charstr(data_path, "/data");
    #endif
  }

  if (sens_spec_fname_root != NULL && bias != NEGINFTY)
    die("ERROR: can't use --bias and --sens-spec together.\n");

  if (extrapolate_tree_fname != NULL &&
      !strcmp(extrapolate_tree_fname, "default")) {
    extrapolate_tree_fname = smalloc(1000 * sizeof(char));
    #if defined(__MINGW32__)
      sprintf(extrapolate_tree_fname,
	      "%s\\exoniphy\\mammals\\cftr25_hybrid.nh", data_path->chars);
    #else
      sprintf(extrapolate_tree_fname, 
              "%s/exoniphy/mammals/cftr25_hybrid.nh", data_path->chars);
    #endif
  }
  if (extrapolate_tree_fname != NULL)
    extrapolate_tree = tr_new_from_file(phast_fopen(extrapolate_tree_fname, "r"));

  /* read alignment */
  msa_fname = argv[optind];
  infile = phast_fopen(msa_fname, "r");

  if (msa_format == UNKNOWN_FORMAT)
    msa_format = msa_format_for_content(infile, 1);
  if (!quiet)
    fprintf(stderr, "Reading alignment from %s...\n",  msa_fname);

  if (msa_format == MAF)
    msa = maf_read(infile, NULL, 1, NULL, NULL, 
                   NULL, -1, TRUE, NULL, NO_STRIP, FALSE);
  else
    msa = msa_new_from_file_define_format(infile, msa_format, NULL);

  phast_fclose(infile);

  if (msa_alph_has_lowercase(msa)) msa_toupper(msa); 
  msa_remove_N_from_alph(msa);

  if (msa_format == SS && msa->ss->tuple_idx == NULL) 
    die("ERROR: Ordered representation of alignment required.\n");
  if (not_informative != NULL)
    msa_set_informative(msa, not_informative);

  /* rename if aliases are defined */
  if (alias_hash != NULL) {
    for (i = 0; i < msa->nseqs; i++) {
      if ((newname = hsh_get(alias_hash, msa->names[i])) != (char*)-1) {
        sfree(msa->names[i]);
        msa->names[i] = copy_charstr(newname);
      }
    }
  }

  /* use filename root for default seqname and/or idpref */
  if (seqname == NULL || idpref == NULL) {
    String *tmp = str_new_charstr(msa_fname);
    if (!str_equals_charstr(tmp, "-")) {
      str_remove_path(tmp);
      str_root(tmp, '.');
      if (idpref == NULL) idpref = copy_charstr(tmp->chars);
      str_root(tmp, '.');         /* apply one more time for double suffix */
      if (seqname == NULL) seqname = tmp->chars;    
    }
    else if (seqname == NULL) seqname = "refseq";
  }

  /* set hmm, tree models, and category map to defaults, if not
     already specified */
  if (hmm == NULL) {
    char *default_hmm = "default.hmm";
    if (indels) {
      if (no_cns) default_hmm = "default-indels-no-cns.hmm";
      else default_hmm = "default-indels.hmm";
    }
    else if (no_cns) default_hmm = "default-no-cns.hmm";
    #if defined(__MINGW32__)
      sprintf(tmpstr, "%s\\exoniphy\\mammals\\%s", data_path->chars, default_hmm);
    #else
      sprintf(tmpstr, "%s/exoniphy/mammals/%s", data_path->chars, default_hmm);
    #endif
    if (!quiet) fprintf(stderr, "Reading default HMM from %s...\n", tmpstr);
    hmm = hmm_new_from_file(phast_fopen(tmpstr, "r"));
    reflect_hmm = TRUE;
    if (bias == NEGINFTY) bias = -3.33;
  }

  if (model_fname_list == NULL) {
    /* Figure out which set of models to use based on G+C content */
    Vector *f = msa_get_base_freqs(msa, -1, -1); 
    gc = vec_get(f, msa->inv_alphabet[(int)'G']) +
      vec_get(f, msa->inv_alphabet[(int)'C']);    
    for (gc_cat = 0; 
         GC_THRESHOLDS[gc_cat] != -1 && gc >= GC_THRESHOLDS[gc_cat]; 
         gc_cat++);
    gc_cat++;                   /* use 1-based index */
    if (!quiet) 
      fprintf(stderr, "(G+C content is %.1f%%; using tree models for G+C category %d)\n",
              gc*100, gc_cat);
    model_fname_list = lst_new_ptr(30);
    #if defined(__MINGW32__)
      sprintf(tmpstr, "%s\\exoniphy\\%s-gc%d", data_path->chars,
	      no_cns ? "models-no-cns" : "models", gc_cat);
    #else
      sprintf(tmpstr, "%s/exoniphy/%s-gc%d", data_path->chars, 
              no_cns ? "models-no-cns" : "models", gc_cat);
    #endif
    str_slurp(fname_str, phast_fopen(tmpstr, "r"));
    str_split(fname_str, NULL, model_fname_list);
    if (!quiet) 
      #if defined(__MINGW32__)
        fprintf(stderr, "Reading default tree models from %s\\exoniphy\\mammals\\*.mod...\n", data_path->chars);
      #else
        fprintf(stderr, "Reading default tree models from %s/exoniphy/mammals/*.mod...\n", data_path->chars);
      #endif
    for (i = 0; i < lst_size(model_fname_list); i++) {
      str = lst_get_ptr(model_fname_list, i);
      #if defined(__MINGW32__)
        sprintf(tmpstr, "%s\\exoniphy\\mammals\\%s", data_path->chars, str->chars);
      #else
        sprintf(tmpstr, "%s/exoniphy/mammals/%s", data_path->chars, str->chars);
      #endif

      str_cpy_charstr(str, tmpstr);
    }
    vec_free(f);
  }

  if (cm == NULL) {
    #if defined(__MINGW32__)
      sprintf(tmpstr, "%s\\exoniphy\\%s", data_path->chars,
	      no_cns ? "default-no-cns.cm" : "default.cm");
    #else
      sprintf(tmpstr, "%s/exoniphy/%s", data_path->chars, 
              no_cns ? "default-no-cns.cm" : "default.cm");
    #endif
    if (!quiet) fprintf(stderr, "Reading default category map from %s...\n", tmpstr);
    if ((cm = cm_read(phast_fopen(tmpstr, "r"))) == NULL)
      die("Unable to open default category map %s...\n", tmpstr);

    if (no_gaps_str == NULL) 
      no_gaps_str = get_arg_list("10,11,20,21,cds5\'ss,cds3\'ss,start_codon,stop_codon");
    if (inform_reqd == NULL) 
      inform_reqd = get_arg_list("10,11,20,21,cds5\'ss,cds3\'ss,start_codon,stop_codon,CDS");
  }

  ncats = cm->ncats + 1;

  /* read tree models */
  if (lst_size(model_fname_list) != ncats) 
    die("ERROR: number of tree models must equal number of site categories.\n");
    
  mod = (TreeModel**)smalloc(sizeof(TreeModel*) * ncats);
  for (i = 0; i < ncats; i++) {
    String *fname = (String*)lst_get_ptr(model_fname_list, i);
    List *pruned_names = lst_new_ptr(msa->nseqs);
    int old_nnodes;

    F = phast_fopen(fname->chars, "r");
    mod[i] = tm_new_from_file(F, 1);
    mod[i]->use_conditionals = 1;
    phast_fclose(F);

    old_nnodes = mod[i]->tree->nnodes;
    /* extrapolate tree and/or prune away extra species */
    if (extrapolate_tree != NULL) {
      double scale = tm_extrapolate_and_prune(mod[i], extrapolate_tree, 
                                              msa, pruned_names);
      if (!quiet) 
        fprintf(stderr, "Extrapolating based on %s (scale=%f)...\n", 
                extrapolate_tree_fname, scale);
    }
    else
      tm_prune(mod[i], msa, pruned_names);

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
  }

  /* disallow gaps, if necessary */
  if (no_gaps_str != NULL) {
    List *l = cm_get_category_list(cm, no_gaps_str, 0);
    for (i = 0; i < lst_size(l); i++) 
      mod[lst_get_int(l, i)]->allow_gaps = FALSE;
    lst_free(l);
  }      

  /* set min informative bases, if necessary */
  if (inform_reqd != NULL) {
    List *l = cm_get_category_list(cm, inform_reqd, 0);
    for (i = 0; i < lst_size(l); i++) 
      mod[lst_get_int(l, i)]->inform_reqd = TRUE;
    lst_free(l);
  }      

  phmm = phmm_new(hmm, mod, cm, reflect_hmm ? backgd_types : NULL, 
                  indels ? NONPARAMETERIC : MISSING_DATA);
                                /* FIXME: allow nonparameteric also */

  /* add bias, if necessary */
  if (bias != NEGINFTY) {
    if (!quiet) fprintf(stderr, "Applying coding bias of %f...\n", bias);
    phmm_add_bias(phmm, backgd_types, bias);
  }

  /* compute emissions */
  phmm_compute_emissions(phmm, msa, quiet);

  /* now produce predictions.  Need to do this in a loop because
     of sens-spec mode */
  if (sens_spec_fname_root != NULL) {    
    phmm_add_bias(phmm, backgd_types, SCALE_RANGE_MIN);
    ntrials = NSENS_SPEC_TRIES;
  }
  else ntrials = 1;

  for (trial = 0; trial < ntrials; trial++) {
        
    if (ntrials > 1 && !quiet)
      fprintf(stderr, "(Sensitivity/specificity trial #%d)\n", trial+1);

    /* run Viterbi */
    if (!quiet)
      fprintf(stderr, "Executing Viterbi algorithm...\n");
    predictions = phmm_predict_viterbi(phmm, seqname, grouptag, idpref,
                                       cds_types);

    /* score predictions */
    if (score) {
      if (!quiet) fprintf(stderr, "Scoring predictions...\n");            
      phmm_score_predictions(phmm, predictions, cds_types, 
                             signal_types, backgd_types, TRUE);
    }

    /* adjust GFF -- absorb helper features, filter out unwanted
       types, add group_id tag */
    gff_group(predictions, grouptag);
    gff_absorb_helpers(predictions, cds_types, cds_absorb_types);
    gff_filter_by_type(predictions, invisible_types, TRUE, NULL);
    gff_group(predictions, grouptag); /* will be ungrouped by gff_filter_by_type */
    gff_add_gene_id(predictions);

    /* convert to coord frame of reference sequence and adjust for
       idx_offset.  FIXME: make clear in help page assuming refidx 1 */
    msa_map_gff_coords(msa, predictions, 0, 1, msa->idx_offset);

    /* output predictions */
    if (sens_spec_fname_root != NULL) { 
      char tmpstr[STR_MED_LEN];
      sprintf(tmpstr, "%s.v%d.gff", sens_spec_fname_root, trial+1);
      gff_print_set(phast_fopen(tmpstr, "w+"), predictions);      
      if (trial < ntrials - 1)     /* also set up for next iteration */
        phmm_add_bias(phmm, backgd_types, (SCALE_RANGE_MAX - SCALE_RANGE_MIN)/
                      (NSENS_SPEC_TRIES-1));
    }
    else                        /* just output to stdout */
      gff_print_set(stdout, predictions);

  } /* ntrials loop */

  if (!quiet)
    fprintf(stderr, "Done.\n");

  return 0;
}
