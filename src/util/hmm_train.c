/* hmm_train - estimation of HMM transition probabilities from labeled
   training data */

/* $Id: hmm_train.c,v 1.3 2004-07-27 17:59:06 acs Exp $
   Written by Adam Siepel, 2002-2004
   Copyright 2002-2004, Adam Siepel, University of California 
*/


#include <stdlib.h>
#include <stdio.h>
#include <msa.h>
#include <hmm.h>
#include <assert.h>
#include <getopt.h>
#include <string.h>
#include <gff.h>
#include <category_map.h>
#include <sufficient_stats.h>
#include <stringsplus.h>
#include <gap_patterns.h>

/* categories for which complex gap patterns are prohibited;
   temporarily hardwired */
#define NO_COMPLEX "cds,cds5'ss,cds3'ss"

void print_usage() {
  printf("\n\
PROGRAM: hmm_train\n\
\n\
USAGE: hmm_train -m <msa_fname_list> -c <category_map_fname> \\\n\
              -g <gff_fname_list> [OPTIONS] > out.hmm \n\
\n\
DESCRIPTION: \n\
    Estimate the transition probabilities of an HMM, based on multiple\n\
    alignments, sequence annotations, and a category map.\n\
\n\
OPTIONS:\n\
\n\
 (required options)\n\
    -m <msa_fname_list>\n\
        List of multiple sequence alignment files.\n\
        Currently, in testing mode, the list must be of length one.\n\
\n\
    -c <category_map_fname>\n\
        File defining mapping of feature types to category\n\
        numbers.\n\
\n\
    -g <gff_fname_list>\n\
        Files in GFF defining sequence\n\
        features to be used in labeling sites.   Frame of reference of \n\
        feature indices is determined feature-by-feature according to \n\
         'seqname' attribute.  Filenames must correspond in number and order\n\
        to the elements of <msa_fname_list>. \n\
\n\
 (alignment options)\n\
    -M <msa_length_list>\n\
        (Mutually exclusive with -m) Assume alignments\n\
        of the specified lengths (comma-separated list) and do not not\n\
        attempt to map the coordinates in the specified GFFs (assume\n\
        they are in the desired coordinate frame).  This option allows\n\
        an HMM to be trained directly from GFFs, without alignments.\n\
        Not permitted with -I.\n\
\n\
    -i PHYLIP|FASTA|MPM|SS \n\
        (default SS) Alignment format.\n\
\n\
 (indel options)
    -I <indel_cat_list>\n\
        Model indels for specified categories.  To have\n\
        nonzero probability for the states corresponding to a\n\
        specified category range, indels must be \"clean\"\n\
        (nonoverlapping), must be assignable by parsimony to a single\n\
        branch in the phylogenetic tree, and must have lengths that\n\
        are exact multiples of the category range size.  Avoid -G with\n\
        this option.  If used in training mode, requires -T.\n\
\n\
    -t <tree_fname>\n\
        Use the specified tree topology when training\n\
        for indels. \n\
\n\
    -n <nseqs> \n\
        Train an indel model for <nseqs>\n\
        sequences, despite that the training alignment has a different\n\
        number.  All (non-trivial) gap patterns are assumed to be\n\
        equally frequent.\n\
\n\
 (other options)\n\
    -q \n\
        Proceed quietly (without updates to stderr).\n\
\n\
    -h\n\
        Print this help message and exit.\n\n");
}

int main(int argc, char* argv[]) {
  FILE* F;
  MSA *msa;
  int *msa_gap_patterns = NULL;
  HMM *hmm = NULL;
  TreeNode *tree = NULL;
  int i, input_format = SS, msa_idx, quiet_mode = FALSE, 
    ncats, nmsas, ncats_unspooled, indel_nseqs = -1;
  String *msa_fname, *gff_fname;
  List *gff_fname_list = NULL, *msa_fname_list = NULL, 
    *msa_length_list = NULL, *model_indels_str = NULL;
  gsl_matrix *traincounts = NULL;
  gsl_vector *begcounts = NULL, *statecounts = NULL;
  CategoryMap *cm = NULL;
  char c;
  GapPatternMap *gpm = NULL;
  GFF_Set *gff;

  while ((c = getopt(argc, argv, "i:g:c:m:M:I:n:t:P:G:qh")) != -1) {
    switch(c) {
    case 'i':
      input_format = msa_str_to_format(optarg);
      if (input_format == -1)
        die("ERROR: bad alignment format.\n");
      break;
    case 'g':
      gff_fname_list = get_arg_list(optarg);
      break;
    case 'c':
      cm = cm_new_string_or_file(optarg);
      break;
    case 'm':
      msa_fname_list = get_arg_list(optarg);
      break;
    case 'M':
      msa_length_list = str_list_as_int(get_arg_list(optarg));
      break;
    case 'I':
      model_indels_str = get_arg_list(optarg);
      break;
    case 'n':
      indel_nseqs = get_arg_int(optarg);
      break;
    case 't':
      if (optarg[0] == '(')     /* in this case, assume topology given
                                   at command line */
        tree = tr_new_from_string(optarg);
      else 
        tree = tr_new_from_file(fopen_fname(optarg, "r"));
      break;
    case 'q':
      quiet_mode = TRUE;
      break;
    case 'h':
      print_usage();
      exit(0);
    case '?':
      die("ERROR: unrecognized option.\n\nType 'hmm_train -h' for usage.\n");
    }
  }

  if (gff_fname_list == NULL) 
    die("ERROR: -g required in training mode.  Type 'hmm_train -h' for usage.\n");
  if (msa_length_list != NULL && msa_fname_list != NULL) 
    die("ERROR: -m and -M are mutually exclusive.  Type 'hmm_train -h' for usage.\n");
  if (model_indels_str != NULL && tree == NULL)
    die("ERROR: -I requires -t.  Type 'hmm_train -h' for usage.\n");
  if (cm == NULL) 
    die("ERROR: category map required.\n");
  
  ncats = cm->ncats + 1;
  ncats_unspooled = cm->unspooler != NULL ? cm->unspooler->nstates_unspooled : 
    ncats;
  nmsas = (msa_length_list != NULL ? lst_size(msa_length_list) : 
           lst_size(msa_fname_list));

  if (model_indels_str != NULL) {
    gpm = gp_create_gapcats(cm, model_indels_str, 
                            indel_nseqs > 0 ? indel_nseqs :
                            (tree->nnodes+1)/2); /* note: no alignment yet */
    ncats = cm->ncats + 1;    /* numbers will change */ 
    ncats_unspooled = cm->unspooler == NULL ? ncats : 
      cm->unspooler->nstates_unspooled;
  }

  /* allocate memory for storage of "training paths" */
  traincounts = gsl_matrix_calloc(ncats_unspooled, ncats_unspooled);
  statecounts = gsl_vector_calloc(ncats_unspooled);
  begcounts = gsl_vector_calloc(ncats_unspooled);
    
  /* create skeleton of new HMM. */
  hmm = hmm_new_nstates(ncats_unspooled, 0, 0);

  /* Main loop: consider each MSA in turn */
  for (msa_idx = 0; msa_idx < nmsas; msa_idx++) {
    if (msa_fname_list != NULL) {
      msa_fname = (String*)lst_get_ptr(msa_fname_list, msa_idx);
      F = fopen_fname(msa_fname->chars, "r");
      if (!quiet_mode)
        fprintf(stderr, "Reading alignment from %s ...\n", 
                F == stdin ? "stdin" : msa_fname->chars);
      msa = msa_new_from_file(F, input_format, NULL);
      fclose(F);

    }
    else {                      /* only lengths of alignments specified */
      msa = msa_new(NULL, NULL, 0, lst_get_int(msa_length_list, msa_idx), NULL);
      /* just a shell in this case */
    }

    gff_fname = (String*)lst_get_ptr(gff_fname_list, msa_idx);
    if (!quiet_mode)
      fprintf(stderr, "Reading annotations from %s ...\n", gff_fname->chars);
    gff_read_set(fopen_fname(gff_fname->chars, "r"));

    /* convert GFF to coordinate frame of alignment */
    if (msa_length_list == NULL) {
      if (!quiet_mode)
        fprintf(stderr, "Mapping annotations to alignment ...\n");
      msa_map_gff_coords(msa, gff, -1, 0, 0, NULL);
    }

    if (model_indels_str != NULL) {
      if (!quiet_mode)
        fprintf(stderr, "Obtaining gap patterns ...\n");
      msa_gap_patterns = smalloc(msa->length * sizeof(int));
      gp_set_phylo_patterns(msa_gap_patterns, msa, tree);
    }

    /* at this point, we don't actually need the alignment anymore;
       if using ordered suff stats (likely with large data sets),
       can free them now, to avoid running out of memory */
    if (msa->ss != NULL) { ss_free(msa->ss); msa->ss = NULL; }
      
    if (!quiet_mode)
      fprintf(stderr, "Labeling sites by category ...\n");       
    msa_label_categories(msa, gff, cm);

    gff_free_set(gff);

    if (model_indels_str != NULL) {
      if (!quiet_mode)
        fprintf(stderr, "Remapping categories according to gap patterns ...\n");

      if (indel_nseqs > 0 && indel_nseqs != msa->nseqs) {
        /* in this case, we'll simply reassign non-trivial gap
           patterns randomly.  This will achieve the desired
           effect with minimal coding, as long as the number of
           sites is not too small (the indel model is probably
           useless anyway if the number is small) */
        int pat, newpat;
        int npatterns = 4 * indel_nseqs - 5;
        int complex_allowed[cm->ncats+1];
        List *no_complex_names, *no_complex_nums;

        if (!quiet_mode)
          fprintf(stderr, "(target number of sequences: %d)\n", indel_nseqs);

        /* set up index indicating by cat no. whether complex gaps
           are allowed */
        for (i = 0; i < ncats; i++) complex_allowed[i] = 1;
        no_complex_names = lst_new_ptr(10);
        str_split(str_new_charstr(NO_COMPLEX), ",", no_complex_names);
        no_complex_nums = cm_get_category_list(cm, no_complex_names, 1);
        for (i = 0; i < lst_size(no_complex_nums); i++)
          complex_allowed[lst_get_int(no_complex_nums, i)] = 0;
        lst_free(no_complex_nums);
        lst_free_strings(no_complex_names); lst_free(no_complex_names);

        /* now reassign all non-null numbers */
        for (i = 0; i < msa->length; ) {
          if ((pat = msa_gap_patterns[i]) != 0) {
            if (complex_allowed[msa->categories[i]])
              newpat = 1 + (random() * (double)npatterns / (RAND_MAX + 1.0));
            /* random number in interval [1, npatterns] */
            else 
              newpat = 1 + (random() * (double)(npatterns-1) / (RAND_MAX + 1.0));
            /* random number in interval [1,npatterns-1] 
               (excludes complex gap pattern) */
            for (; i < msa->length && msa_gap_patterns[i] == pat; i++)
              msa_gap_patterns[i] = newpat; /* change for whole sequence */
          }
          else i++;
        }
      }

      /* obtain gapped category number for each site */
      for (i = 0; i < msa->length; i++) 
        if (gpm->cat_x_pattern_to_gapcat[msa->categories[i]] != NULL)
          msa->categories[i] = gpm->cat_x_pattern_to_gapcat[msa->categories[i]][msa_gap_patterns[i]];
    }

    if (!quiet_mode)
      fprintf(stderr, "Unspooling categories ...\n");
    cm_spooled_to_unspooled(cm, msa->categories, msa->length);
      
    if (!quiet_mode)
      fprintf(stderr, "Collecting training data ...\n");
    hmm_train_update_counts(traincounts, statecounts, begcounts, 
                            msa->categories, msa->length, 
                            ncats_unspooled);

    if (msa_gap_patterns != NULL) free(msa_gap_patterns);
    msa_free(msa);
  }

  /* now train HMM, using cumulative data */
  hmm_train_from_counts(hmm, traincounts, NULL, statecounts, NULL, 
                        begcounts, NULL);

  /* if modeling indels, adjust begin transitions so probability is
     distributed among different "gap pattern" states that all
     correspond to the same ungapped state (category); this helps
     avoid problems that occur when training on a few large sequences
     (e.g., whole chromosomes) and then testing on many shorter ones */
  if (model_indels_str != NULL) {
    double tprob[gpm->ncats]; 
    int nst[gpm->ncats];  /* total prob and number of states per
                             spooled, ungapped category */ 
    for (i = 0; i < gpm->ncats; i++) tprob[i] = nst[i] = 0;
    for (i = 0; i < hmm->nstates; i++) {
      if (gsl_vector_get(hmm->begin_transitions, i) > 0) 
        /* have to go from unspooled space to spooled space, then to
           ungapped space (HMM states correspond to unspooled,
           gapped categories).  Note that states with nonzero begin
           probs shouldn't be conditioned on other states. */
        tprob[gpm->gapcat_to_cat[cm_unspooled_to_spooled_cat(cm, i)]] += 
          gsl_vector_get(hmm->begin_transitions, i);
      nst[gpm->gapcat_to_cat[cm_unspooled_to_spooled_cat(cm, i)]]++;
    }
    for (i = 0; i < hmm->nstates; i++) 
      if (tprob[gpm->gapcat_to_cat[cm_unspooled_to_spooled_cat(cm, i)]] > 0)
        gsl_vector_set(hmm->begin_transitions, i, 
                       tprob[gpm->gapcat_to_cat[cm_unspooled_to_spooled_cat(cm, i)]] / 
                       nst[gpm->gapcat_to_cat[cm_unspooled_to_spooled_cat(cm, i)]]);
    /* (uniform prior) */
  }

  /* write trained HMM */
  hmm_print(stdout, hmm);

  if (!quiet_mode) fprintf(stderr, "Done.\n");

  return 0;
}
