/* label - label the columns of alignment(s) by category */

/* $Id: exoniphy.c,v 1.10 2004-06-29 03:18:03 acs Exp $
   Written by Adam Siepel, 2002-2004
   Copyright 2002-2004, Adam Siepel, University of California 
*/


#include <stdlib.h>
#include <stdio.h>
#include <phylo_hmm.h>
#include <assert.h>
#include <getopt.h>
#include <string.h>
#include <tree_likelihoods.h>
#include <gff.h>
#include <bed.h>
#include <category_map.h>
#include <dgamma.h>
#include <numerical_opt.h>
#include <sufficient_stats.h>
#include <stringsplus.h>
#include <gap_patterns.h>

/* default background feature types; used when scoring predictions and
   reflecting HMM */
#define DEFAULT_BACKGD_CATS "background,CNS"

/* default "cds" and "signal" feature tupes */
#define DEFAULT_CDS_CATS "CDS,start_codon"
#define DEFAULT_SIGNAL_CATS "stop_codon,5'splice,3'splice"
                                /* cat names that aren't present will
                                   be ignored */

#define DEFAULT_CATMAP "PHASTHOME/data/exoniphy/default.cm"

/* parameters controlling evaluation of Sn/Sp tradeoff (see -Y option) */
#define SCALE_RANGE_MIN -20
#define SCALE_RANGE_MAX 10
#define NSENS_SPEC_TRIES 10

void print_usage() {
    printf("
PROGRAM: exoniphy

USAGE: exoniphy --hmm <file> --tree-models <list> <msa_fname> > predictions.gff

FIXME: need default HMM and tree models options (mammals), --tree option
what about --tree for indel model?  


DESCRIPTION: 
    

OPTIONS:
    --msa-format, -i PHYLIP|FASTA|MPM|SS 
        (default FASTA) Alignment format.
 
    --catmap, -c <category_map_fname>
        (required) File defining mapping of feature types to category
        numbers.

    --tree-models, -m <model_fname_list>
        List of files defining a tree model
        for each functional category.  Order of models must correspond
        to order of states in HMM.  Tree model
        files may be produced with phyloFit.  Note: this option
        is given a special interpretation with -D (see below).

    --hmm, -H <hmm_fname>
        (for use with -d; indicates testing mode)  Name of HMM file,
        defining the probability of transition from each functional
        category to each other.  Generally a file is used that was
        produced by this same program, running in training mode.

    --seqname, -s <name>
        (For use with -f, -b, or -Q)  Use specified string as the \"seqname\"
        field in GFF output (e.g., chr22).  By default, the root of the input file name is used.

    --quiet, -q 
        Proceed quietly (without updates to stderr).

    --help -h
        Print this help message.


Auxiliary options

    --reflect-strand, -U 
        Given an HMM describing the forward
        strand, create a larger HMM that allows for features on both
        strands by \"reflecting\" the HMM about all states associated with background categories
        (see -V).  The new HMM will be used for predictions on both strands.

    --score, -S
        Assign log-odds scores to predictions, equal to their log total
        probability under an exon model minus their log total probability
        under a background model.  The exon model can be altered using the --cds-types and --signal-types options and the background model can be altered using the --backgd-types option.

    --bias, -b <val>         /* FIXME: use log */
        Predict with a \"coding bias\" equal to the specified value.  The specified value is added to the log transition probabilities from background states to 
        non-background states (see -V), then all transition probabilities are renormalized.  If the value is positive, more predictions will tend to be made and sensitivity will tend to increase, at some cost to specificity; if the value is negative, fewer predictions will tend to be made, and specificity will tend to improve, at some cost to sensitivity.

    --sens-spec, -Y <fname-root
        Make predictions for a range of different values of coding biase (see --bias), and write results to files with given filename root.  The range is fixed at %d to %d, and %d different sets of predictions are produced.  Allows analysis of sensitivity/specificity tradeoff.

    --cds-types, -C <list>
        Feature types that represent protein-coding regions (default value: \"%s\").  Used when scoring predictions and filling out 'frame' field in GFF output.

    --backgd-types, -B <list>
        Feature types to be considered \"background\" (default value: \"%s\").  Associated
        states are considered \"background states\".  Affects --reflect-strand, --score, and --bias.  

    --signal-types, -L <list>
        (for use with --score)  Types of features to be considered \"signals\" during scoring (default value: \"%s\").  One score is produced for each CDS feature (as defined by --cds-types) and adjacent signal features; the score is then assigned to the CDS feature.  

    --indel-types, -I <list>
        Model indels for features of the specified types.  To have
        nonzero probability for the states corresponding to a
        specified category range, indels must be \"clean\"
        (nonoverlapping), must be assignable by parsimony to a single
        branch in the phylogenetic tree, and must have lengths that
        are exact multiples of the category range size.

    --gc-ranges, -D <range-cutoffs>
        (Changes interpretation of -d) Use different sets of tree
        models, depending on the G+C content of the input alignment.
        The list <gc-thresholds> must consist of x ordered values in
        (0,1), defining x+1 G+C classes.  The argument to -d must then
        consist of the names of x+1 files, each of which contains a list
        of tree-model filenames.

    --no-gaps, -W <cat_list>
        Prohibit gaps in the specified categories (gaps result in
        emission probabilities of zero).  By default, gaps are treated
        as missing data.\n\n", SCALE_RANGE_MIN, 
           SCALE_RANGE_MAX, NSENS_SPEC_TRIES);
}

/* FIXME: put this in phylo_hmm.c; convert to use phmm */
/* multiply all transition probabilities from designated background
   categories to non-background categories by the specified factor,
   then renormalize.  This provides a simple "knob" for the
   sensitivity/specificity tradeoff.  Mathematically, it's a way of
   changing the expected number of predicted features.  Category 0 is
   assumed as a background category */
void scale_trans_from_backgd(PhyloHmm *phmm, List *backgd_cat_names, 
                             double scale_factor) {
  /* FIXME: check cat nos */
  double is_backgd_cat[phmm->cm->ncats+1];
  int i, j;

  is_backgd_cat[0] = 1;
  for (i = 1; i <= phmm->cm->ncats; i++) is_backgd_cat[i] = 0;
  if (backgd_cat_names != NULL) {
    List *backgd_cat_nos = cm_get_category_list(phmm->cm, backgd_cat_names, 0); 
    for (i = 0; i < lst_size(backgd_cat_nos); i++) 
      is_backgd_cat[lst_get_int(backgd_cat_nos, i)] = 1;
    lst_free(backgd_cat_nos);
  }

  for (i = 0; i < phmm->hmm->nstates; i++)
    if (is_backgd_cat[phmm->state_to_cat[i]])
      for (j = 0; j < phmm->hmm->nstates; j++)
        if (!is_backgd_cat[phmm->state_to_cat[j]] && 
            mm_get(phmm->hmm->transition_matrix, i, j) != 0) 
          mm_set(phmm->hmm->transition_matrix, i, j, 
                 scale_factor * 
                 mm_get(phmm->hmm->transition_matrix, i, j));

  hmm_renormalize(phmm->hmm);
}

int main(int argc, char* argv[]) {

  /* variables for options, with defaults */
  int msa_format = FASTA;
  int quiet = FALSE, reflect_hmm = FALSE, score = FALSE;
  double mult_trans_probs = -1;
  char *seqname = NULL, *sens_spec_fname_root = NULL;
  List *model_fname_list = NULL, *no_gaps_str = NULL, *model_indels_str = NULL,
    *penalize_gaps_str = NULL, *gc_thresholds = NULL;
  List *backgd_cats = get_arg_list(DEFAULT_BACKGD_CATS), 
    *cds_cats = get_arg_list(DEFAULT_CDS_CATS), 
    *signal_cats = get_arg_list(DEFAULT_SIGNAL_CATS);

  /* other variables */
  FILE* F;
  PhyloHmm *phmm;
  MSA *msa;
  TreeModel **mod;
  HMM *hmm = NULL;
  CategoryMap *cm = NULL;
  TreeNode *tree = NULL;
  GapPatternMap *gpm = NULL;
  GFF_Set *predictions;
  int *msa_gap_patterns = NULL;
  char c;
  int i, j, ncats, ncats_unspooled,  trial, ntrials;

  while ((c = getopt(argc, argv, "i:c:H:m:s:B:T:L:I:W:b:D:SYUhq")) != -1) {
    switch(c) {
    case 'i':
      msa_format = msa_str_to_format(optarg);
      if (msa_format == -1) die("ERROR: bad alignment format.\n");
      break;
    case 'c':
      cm = cm_new_string_or_file(optarg);
      break;
    case 'm':
      model_fname_list = get_arg_list(optarg);
      break;
    case 'H':
      if (!quiet) fprintf(stderr, "Reading HMM from %s...", optarg);
      hmm = hmm_new_from_file(fopen_fname(optarg, "r"));
      break;
    case 'U':
      reflect_hmm = TRUE;
      break;
    case 'B':
      lst_free_strings(backgd_cats); lst_free(backgd_cats); /* free defaults */
      backgd_cats = get_arg_list(optarg);
      break;
    case 'T':
      lst_free_strings(cds_cats); lst_free(cds_cats); /* free defaults */
      cds_cats = get_arg_list(optarg);
      break;
    case 'L':
      lst_free_strings(signal_cats); lst_free(signal_cats); /* free defaults */
      signal_cats = get_arg_list(optarg);
      break;
    case 'S':
      score = TRUE;
      break;
    case 'b':
      mult_trans_probs = get_arg_dbl(optarg);
      break;
    case 'Y':
      sens_spec_fname_root = optarg;
      break;
    case 'W':
      no_gaps_str = get_arg_list(optarg);
      break;
    case 'I':
      model_indels_str = get_arg_list(optarg);
      break;
    case 'D':
      gc_thresholds = str_list_as_dbl(get_arg_list(optarg));
      for (i = 0; i < lst_size(gc_thresholds); i++)
        if (lst_get_dbl(gc_thresholds, i) <= 0 || 
            lst_get_dbl(gc_thresholds, i) >= 1) 
          die("ERROR: Bad G+C threshold.\n");
      lst_qsort_dbl(gc_thresholds, ASCENDING);
      break; 
    case 's':
      seqname = optarg;
      break;
    case 'q':
      quiet = TRUE;
      break;
    case 'h':
      print_usage();
      exit(0);
    case '?':
      fprintf(stderr, "ERROR: unrecognized option.  Try 'exoniphy -h' for help.\n");
      exit(1);
    }
  }

  /* Check validity of arguments */
  if (optind != argc - 1) 
    die("ERROR: alignment filename is required argument.  Try 'exoniphy -h' for help.\n");

  if (model_fname_list == NULL) 
    die("ERROR: --tree-models is a required argument.  Try 'exoniphy -h' for help.\n");

  /* read alignment */
  if (!quiet)
    fprintf(stderr, "Reading alignment from %s ...\n", 
            !strcmp(argv[optind], "-") ? "stdin" : argv[optind]);
  
  msa = msa_new_from_file(fopen_fname(argv[optind], "r"), msa_format, NULL);
  msa_remove_N_from_alph(msa);
  if (msa_format == SS && msa->ss->tuple_idx == NULL) 
    die("ERROR: Ordered representation of alignment required.\n");

  /* use filename root for default seqname */
  if (seqname == NULL) {
    String *tmp = str_new_charstr(argv[optind]);
    str_remove_path(tmp);
    str_root(tmp, '.');
    seqname = tmp->chars;
  }

  /* get default cat map if not specified */
  if (cm == NULL)
    cm = cm_read(fopen_fname(DEFAULT_CATMAP, "r"));

  ncats = cm->ncats + 1;
  ncats_unspooled = cm->unspooler != NULL ? cm->unspooler->nstates_unspooled : 
    ncats;

  if (gc_thresholds != NULL) {
    double gc;
    String *models_fname;
    char tmpstr[STR_SHORT_LEN];
    gsl_vector *f = msa_get_base_freqs(msa, -1, -1); 
    /* note: does not consider categories */
    assert(msa->inv_alphabet[(int)'G'] >= 0 && 
           msa->inv_alphabet[(int)'C'] >= 0);
    gc = gsl_vector_get(f, msa->inv_alphabet[(int)'G']) +
      gsl_vector_get(f, msa->inv_alphabet[(int)'C']);
    assert(lst_size(model_fname_list) == lst_size(gc_thresholds) + 1);
    for (i = 0; models_fname == NULL && i < lst_size(gc_thresholds); i++)
      if (gc < lst_get_dbl(gc_thresholds, i))
        models_fname = lst_get_ptr(model_fname_list, i);
    if (models_fname == NULL) 
      models_fname = lst_get_ptr(model_fname_list, i); /* last element */
    
    if (!quiet) 
      fprintf(stderr, "G+C content is %.1f%%; using models for partition %d (%s) ...\n", gc*100, i, models_fname->chars);
    
    /* this trick makes it as if the correct set of models had been
       specified directly */
    sprintf(tmpstr, "*%s", models_fname->chars);
    lst_free_strings(model_fname_list); lst_clear(model_fname_list);
    model_fname_list = get_arg_list(tmpstr);
    
    gsl_vector_free(f);
  }

  /* read tree models */
  if (lst_size(model_fname_list) != ncats) 
    die("ERROR: number of tree models must equal number of site categories.\n");
    
  mod = (TreeModel**)smalloc(sizeof(TreeModel*) * ncats);
  for (i = 0; i < ncats; i++) {
    String *fname = (String*)lst_get_ptr(model_fname_list, i);
    if (!quiet)
      fprintf(stderr, "Reading tree model from %s ...\n", fname->chars);
    F = fopen_fname(fname->chars, "r");
    mod[i] = tm_new_from_file(F); 
    mod[i]->use_conditionals = 1;
    if (mod[i]->tree != NULL) { /* FIXME: put this in phylo_hmm.c */
      if (tree == NULL)
        tree = mod[i]->tree;
      else if (mod[i]->tree->nnodes != tree->nnodes) 
        die("ERROR: trees of models have different numbers of nodes.\n");
    }
    fclose(F);
  }

  if (tree != NULL && msa->nseqs * 2 - 1 != tree->nnodes) 
    die("ERROR: number of leaves in trees must equal number of sequences in alignment.\n");

  /* set up gap cats */
  /* FIXME: put in phylo-hmm.c? */
  if (no_gaps_str != NULL) {
    List *l = cm_get_category_list(cm, no_gaps_str, 0);
    for (i = 0; i < lst_size(l); i++) 
      mod[lst_get_int(l, i)]->allow_gaps = 0;
    lst_free(l);
  }      
  else if (penalize_gaps_str != NULL) {
    List *l = cm_get_category_list(cm, penalize_gaps_str, 0);
    for (i = 0; i < lst_size(l); i++) 
      mod[lst_get_int(l, i)]->allow_but_penalize_gaps = 1;
    lst_free(l);
  }      

  phmm = phmm_new(hmm, mod, cm, reflect_hmm ? backgd_cats : NULL, 
                  model_indels_str, msa->nseqs);

  /* if modeling indels, obtain gap pattern for each site in the alignment */
  /* FIXME: put some of this in phylo_hmm.c? */
  if (model_indels_str != NULL) {
    if (!quiet)
      fprintf(stderr, "Obtaining gap patterns ...\n");
    msa_gap_patterns = smalloc(msa->length * sizeof(int));
    gp_set_phylo_patterns(msa_gap_patterns, msa, tree);
  }

  /* multiply transition probs, if necessary */
  if (mult_trans_probs > 0) 
    scale_trans_from_backgd(phmm, backgd_cats, mult_trans_probs);

  phmm_compute_emissions(phmm, msa, quiet);

  /* redefine emissions for indel states, if necessary */
  /* FIXME: put in phylo_hmm.c ! */
  if (model_indels_str != NULL) {
    if (!quiet)
      fprintf(stderr, "Adjusting emission probs according to gap patterns ...\n");
    for (i = phmm->hmm->nstates - 1; i >= 0; i--) {
                                /* by going backwards, we ensure that
                                   the "base" state (with gap pattern
                                   == 0) is visited last (see
                                   puzzler.c and gap_patterns.c) */
      if (phmm->state_to_pattern[i] >= 0) {
        double *orig_emissions = phmm->emissions[i];
        if (phmm->state_to_pattern[i] > 0)
          phmm->emissions[i] = smalloc(msa->length * sizeof(double));
                                /* otherwise, use the array already
                                   allocated */
        for (j = 0; j < msa->length; j++) 
          phmm->emissions[i][j] = 
            (msa_gap_patterns[j] == phmm->state_to_pattern[i] ? 
             orig_emissions[j] : NEGINFTY);                                 
      }
    }
  }

  /* need to do the part below in a loop in the case of sens-spec
     mode */
  if (sens_spec_fname_root != NULL) {    
    scale_trans_from_backgd(phmm, backgd_cats, exp(SCALE_RANGE_MIN));
    ntrials = NSENS_SPEC_TRIES;
  }
  else ntrials = 1;

  for (trial = 0; trial < ntrials; trial++) {
        
    if (ntrials > 1 && !quiet)
      fprintf(stderr, "(Sensitivity/specificity trial #%d)\n", trial+1);

    /* run Viterbi */
    if (!quiet)
      fprintf(stderr, "Executing Viterbi algorithm ...\n");
    predictions = phmm_predict_viterbi(phmm, seqname, cds_cats);

    /* score predictions */
    if (score) {
      if (!quiet) fprintf(stderr, "Scoring predictions...\n");            
      phmm_score_predictions(phmm, predictions, cds_cats, 
                             signal_cats, backgd_cats, TRUE);
    }

    /* convert to coord frame of reference sequence and adjust for
       idx_offset.  FIXME: make clear in help page assuming refidx 1 */
    msa_map_gff_coords(msa, predictions, 0, 1, msa->idx_offset, NULL);

    /* output predictions */
    if (sens_spec_fname_root != NULL) { 
      char tmpstr[STR_MED_LEN];
      sprintf(tmpstr, "%s.v%d.gff", sens_spec_fname_root, trial+1);
      gff_print_set(fopen_fname(tmpstr, "w+"), predictions);      
      if (trial < ntrials - 1)     /* also set up for next iteration */
        scale_trans_from_backgd(phmm, backgd_cats, 
                                exp((SCALE_RANGE_MAX - SCALE_RANGE_MIN)/
                                    (NSENS_SPEC_TRIES-1)));
    }
    else                        /* just output to stdout */
      gff_print_set(stdout, predictions);

  } /* ntrials loop */

  if (!quiet)
    fprintf(stderr, "Done.\n");

  return 0;
}
