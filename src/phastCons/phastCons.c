#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <misc.h>
#include <phylo_hmm.h>
#include <sufficient_stats.h>
#include <bed.h>

#define MIN_BLOCK_SIZE 30
/* defines block size for -x option */

void usage(char *prog) {
  printf("\n\
PROGRAM:     %s

DESCRIPTION: Identify conserved elements or produce conservation
             scores (posterior probabilities of conserved states),
             given a multiple alignment and a phylo-HMM.  The
             state-transition structure of the phylo-HMM may be
             specified explicitly (see --hmm) or implicitly (see
             --auto-dgamma and --cut-rates).  Phylogenetic models
             (*.mod files) must be provided in either case (they can
             be produced with 'phyloFit').  By default, the posterior
             probability of the selected states in the HMM (see
             --states) is written to stdout in a simple two-column
             format (position and probability, tab separated).
             Prediction of discrete elements using the Viterbi
             algorithm is available via the --viterbi option.
             Predicted elements may optionally be assigned log-odds
             scores based on the model (see --score).

             This program is written in a general enough way that it
             can be useful for other things besides analyzing rates of
             substitution and evolutionary conservation.  E.g., by
             giving it a simple gene-finding phylo-HMM (e.g., with a
             state for non-coding regions and three states for the
             three codon positions), and specifying the coding states
             with --states, you can obtain a \"coding potential\"
             score, or even a set of crude exon predictions (using
             --viterbi).  (The program 'exoniphy' is a more
             full-featured exon predictor.)

USAGE:       %s [OPTIONS] <msa_fname> <mod_fname_list>

             <mod_fname_list> must be a comma- or whitespace-
             delimited list of *.mod files, as produced by 'phyloFit';
             <msa_fname> should be the name of the multiple alignment
             file, which can use any of several possible file formats
             (see --msa-format).

OPTIONS:

    --states, -S <state_list>
        States of interest in the phylo-HMM, specified by number
        (indexing starts with 1).  Default value is 1.
        This option defines the meaning of the program's output.  For
        example, --states \"1,2,3\" causes output of the sum of the
        posterior probabilities for states 1, 2, and 3, and/or of
        regions in which the Viterbi path coincides with (any of)
        states 1, 2, or 3 (see --viterbi).

    --hmm, -H <hmm_fname>
        Name of HMM file, explicitly defining the probabilities of all
        state transitions.  States in the file must correspond in
        number and order to phylogenetic models in <mod_fname_list>.
        Expected file format is as produced by 'hmm_train.'

    --viterbi, -V <fname>
        Compute the Viterbi (maximum likelihood) path and write
        start and end coordinates of predicted elements to specified
        file.  Output is in BED format, unless <fname> has suffix
        \".gff\", in which case output is in GFF.  The features
        output will represent maximal segments in which the Viterbi
        path remains in the selected set of states (see --states).

    --lnl, -L <fname>
        Compute total log likelihood using the forward algorithm and
        write it to the specified file.

    --no-post-probs, -n
        Suppress output of posterior probabilities.  Useful if only
        Viterbi path or likelihood is of interest (saves computation
        time and disk space).

    --msa-format, -i PHYLIP|FASTA|PSU|SS|LAV|MAF
        (Default SS) Alignment file format.

    --auto-dgamma, -k <nratecats>
        (Alternative to --hmm; specify only one *.mod file with this
        option) Use an HMM with <nratecats> states, corresponding to
        \"scaled\" versions of the given phylogenetic model.  The
        scaling (rate) constants will be computed using Yang's (1994)
        discrete gamma method, and all transition probabilities will
        be defined by a single autocorrelation parameter, lambda, as
        in Felsenstein and Churchill (1996).  (See Siepel & Haussler,
        2003 for details on how these methods are combined.)  The
        parameter lambda will be estimated by maximum likelihood.  The
        specified *.mod file must have been produced by phyloFit using
        the -k option (an estimate of the shape parameter 'alpha' is
        required).

FIXME: actually, can do this with empirical rates as well

    --lambda, -l <lambda>
        (Optionally use with --auto-dgamma) Fix lambda at the
        specified value rather than estimating it by maximum likelihood.
        Allowable range is 0-1.  With k rate categories, the
        transition probability between state i and state j will be
        lambda * I(i == j) + (1-lambda)/k, where I is the indicator
        function.  Thus, lambda = 0 implies no autocorrelation and
        lambda = 1 implies perfect autocorrelation.

    --cut-rates, -c <cut_idx>
        (Alternative to --hmm and --auto-dgamma; specify only one
        phylogenetic model) Define a simple HMM with a conserved state
        (state 1) and a non-conserved state (state 2), such that rate
        categories 1-<cut_idx> from the specified *.mod file are
        associated with the conserved state, and the remaining rate
        categories are associated with the non-conserved state.  The
        *.mod file must allow for rate variation, via either the
        discrete gamma model (-k option to phyloFit) or a
        non-parametric alternative (-K option to phyloFit).  The
        new phylogenetic model associated with each state will be a
        mixture model with one or more mixture components, whose
        (unnormalized) mixing proportions are given by the *.mod file,
        either implicitly (discrete gamma case) or explicitly
        (non-parameteric case).  The transition probabilities of the
        HMM will be estimated by maximum likelihood.

    --cut-params, -p <alpha, beta>
        (Optionally use with --cut-rates) Fix the transition
        probabilities of the HMM at the specified values, rather than
        estimating them by maximum likelihood.  Here, alpha is the
        probability of transitioning from the conserved to the
        non-conserved state, and beta is the probability of the reverse
        transition (probabilities of self transitions are thus 1-alpha
        and 1-beta).

    --refidx, -r <refseq_idx> 
        Use coordinate frame of specified sequence (the value 0
        indicates the frame of the entire multiple alignment).
        Default value is 1 (first sequence assumed reference).

    --suppress-missing, -x 
        Suppress posterior probabilities where cross-species alignment
        data does not appear to be available.  The heuristic used is
        to consider blocks of %d or more sites in which only the
        reference sequence is present (as determined by --refidx) to
        be regions of \"missing data\".

    --reflect-strand, -U <pivot_states>
        (Optionally use with --hmm) Given an HMM describing the
        forward strand, create a larger HMM that allows for features
        on both strands by \"reflecting\" the original HMM about the
        specified \"pivot\" states.  The new HMM will be used for
        prediction on both strands.  For example, suppose you use -hmm
        to specify a four-state gene-finding HMM (state 1 noncoding,
        states 2, 3, and 4 coding).  You can use --reflect-strand 1 to
        dynamically create a new HMM with two versions of states 2, 3,
        and 4, one for each strand, and with appropriately defined
        transition probabilities.  Then you can use --states 2,3,4 to
        obtain posterior probabilities and/or a Viterbi path
        describing coding regions on either strand.

    --score, -s
        (Optionally use with --viterbi) Assign a log-odds score to
        each predictions, equal to the log total probability of the
        region in question under the portion of the model defined by
        --states minus its log total probability under the remainder
        of the model.  To compute these scores, the states of the
        model are partitioned and transition probabilities are
        re-normalized.

    --quiet, -q 
        Proceed quietly (without updates to stderr).

    --help, -h
        Print this help message.


REFERENCES:

    J. Felsenstein and G. Churchill.  1996. A hidden Markov model
      approach to variation among sites in rate of evolution. 
      Mol. Biol. Evol., 13:93-104. 

    A. Siepel and D. Haussler.  2003.  Combining phylogenetic and
      hidden Markov models in biosequence analysis.  RECOMB '03.

    Z. Yang. 1994. Maximum likelihood phylogenetic estimation from
      DNA sequences with variable rates over sites: approximate
      methods. J. Mol. Evol., 39:306-314.\n\n", prog, prog, MIN_BLOCK_SIZE);

  exit(0);
}

/* attempts to merge cats by incorporating together into ranges */
void collapse_cats(CategoryMap *cm, List *cats_to_merge) {
  /* assumes all ranges are of size one initially */
  int i, beg, end = -INFTY;
  lst_qsort_int(cats_to_merge, ASCENDING);
  for (i = 0; i < lst_size(cats_to_merge); i++) {
    int cat = lst_get_int(cats_to_merge, i);
    if (cat == end + 1) {
      cm_free_category_range(cm->ranges[cat]);
      cm->ranges[cat] = cm->ranges[beg];
      end = cat;
      cm->ranges[beg]->end_cat_no = end;
    }
    else beg = end = cat;
  }
}

int main(int argc, char *argv[]) {

  /* arguments and defaults */
  int post_probs = TRUE, no_missing = FALSE, score = FALSE, quiet = FALSE, 
    gff = FALSE;
  int dgamma_nrates = -1, cut_rate_idx = -1, refidx = 1;
  double lambda = -1, alpha = -1, beta = -1;
  msa_format_type msa_format = SS;
  FILE *viterbi_f = NULL, *lnl_f = NULL;
  List *states = NULL, *pivot_states = NULL;
  HMM *hmm = NULL;

  struct option long_opts[] = {
    {"states", 1, 0, 'S'},
    {"hmm", 1, 0, 'H'},
    {"viterbi", 1, 0, 'V'},
    {"no-post-probs", 0, 0, 'n'},
    {"msa-format", 1, 0, 'i'},
    {"auto-dgamma", 1, 0, 'k'},
    {"lambda", 1, 0, 'l'},
    {"cut-rates", 1, 0, 'c'},
    {"cut-params", 1, 0, 'p'},
    {"refidx", 1, 0, 'r'},
    {"suppress-missing", 0, 0, 'x'},
    {"reflect-strand", 0, 0, 'U'},
    {"lnl", 1, 0, 'L'},
    {"score", 1, 0, 's'},
    {"quiet", 1, 0, 'q'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  /* other vars */
  char c;
  int opt_idx, i;
  List *tmpl = NULL;
  MSA *msa = NULL;
  double lnl = INFTY;
  String *tmpstr;
  TreeModel **mod;
  PhyloHmm *phmm;

  while ((c = getopt_long(argc, argv, "S:H:V:ni:k:l:c:p:r:xL:s:qh", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'S':
      states = get_arg_list_int(optarg);
      break;
    case 'H':
      if (!quiet) fprintf(stderr, "Reading HMM from %s...", optarg);
      hmm = hmm_new_from_file(fopen_fname(optarg, "r"));
      break;
    case 'V':
      viterbi_f = fopen_fname(optarg, "w+");
      tmpstr = str_new_charstr(optarg);
      if (str_ends_with_charstr(tmpstr, ".gff")) gff = TRUE;
      str_free(tmpstr);
      break;
    case 'n':
      post_probs = FALSE;
      break;
    case 'i':
      msa_format = msa_str_to_format(optarg);
      if (msa_format == -1) die("ERROR: bad argument to --msa-format\n");
      break;
    case 'k':
      dgamma_nrates = get_arg_int_bounds(optarg, 2, 100);
      break;
    case 'l':
      lambda = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'c':
      cut_rate_idx = get_arg_int_bounds(optarg, 1, 100);
      break;
    case 'p':
      tmpl = get_arg_list_dbl(optarg);
      if (lst_size(tmpl) != 2) die("ERROR: bad argument to --cut-params.\n");
      alpha = lst_get_dbl(tmpl, 0);
      beta = lst_get_dbl(tmpl, 1);
      if (alpha <= 0 || alpha >= 1 || beta <= 0 || beta >= 1)
        die("ERROR: bad argument to --cut-params.\n");
      lst_free(tmpl);
      break;
    case 'r':
      refidx = get_arg_int_bounds(optarg, 0, INFTY);
      break;
    case 'x':
      no_missing = TRUE;
      break;
    case 'U':
      pivot_states = get_arg_list(optarg); /* we want strings not ints
                                             for phmm_new */
      break;
    case 'L':
      lnl_f = fopen_fname(optarg, "w+");
      break;
    case 's':
      score = TRUE;
      break;
    case 'q':
      quiet = TRUE;
      break;
    case 'h':
      usage(argv[0]);
    case '?':
      die("Bad argument.  Try '%s -h'.\n", argv[0]);
    }
  }

  if ((hmm != NULL && dgamma_nrates > 1) ||
      (hmm != NULL && cut_rate_idx != -1) ||
      (dgamma_nrates > 1 && cut_rate_idx != -1))
    die("ERROR: --hmm, --auto-dgamma, and --cut-rates are mutually exclusive.\n");

  if (optind != argc - 2) 
    die("ERROR: missing required arguments.  Try '%s -h'.\n", argv[0]);

  /* read these first (msa may take a while) */
  tmpl = get_arg_list(argv[optind+1]);
  mod = (TreeModel**)smalloc(sizeof(TreeModel*) * lst_size(tmpl));
  for (i = 0; i < lst_size(tmpl); i++) {
    String *fname = lst_get_ptr(tmpl, i);
    if (!quiet)
      fprintf(stderr, "Reading tree model from %s ...\n", fname->chars);
    mod[i] = tm_new_from_file(fopen_fname(fname->chars, "r"));
    mod[i]->use_conditionals = 1; /* FIXME: necessary? */
  }
  
  if (!quiet)
    fprintf(stderr, "Reading alignment from %s ...\n", argv[optind]);
  msa = msa_new_from_file(fopen_fname(argv[optind], "r"), msa_format, NULL);

  if (msa_format == SS) {
    if (msa->ss->tuple_idx == NULL) 
      die("ERROR: Ordered representation of alignment required.\n");
    msa_remove_N_from_alph(msa);
  }

  if (states == NULL) {
    states = lst_new_int(1);
    lst_push_int(states, 0);
  }
  else {
    for (i = 0; i < lst_size(states); i++) {
      int state = lst_get_int(states, i);
      if (state < 1 || state > hmm->nstates) 
        die("ERROR: selected state out of range.\n");
      lst_set_int(states, i, state - 1); /* internally use 0-based
                                            indexing */
    }
  }

  /* set up PhyloHmm */
  phmm = phmm_new(hmm, mod, NULL, pivot_states, NULL, -1);

  if (dgamma_nrates != -1) {
    if (!quiet) 
      fprintf(stderr, "Creating %d scaled versions of model...\n", 
              dgamma_nrates);
    phmm_rates_cross(phmm, dgamma_nrates, lambda, TRUE);
  }

  /* FIXME: something similar for cut-rates version */

  /* compute emissions */
  phmm_compute_emissions(phmm, msa, quiet);

  /* fit lambda, if necessary */
  if (dgamma_nrates > 1 && lambda == -1) {
    if (!quiet) fprintf(stderr, "Finding MLE for lambda ...");
    lambda = phmm_fit_lambda(phmm, msa->length);
    if (!quiet) fprintf(stderr, "lambda = %f\n", lambda);
    phmm_update_cross_prod(phmm, lambda);
  }
    
  /* Viterbi */
  if (viterbi_f != NULL) {
    GFF_Set *predictions;

    if (lst_size(states) > 1) 
      /* if possible, we want to merge categories of specified states,
         so that the "union" of states is automatically considered by
         phmm_viterbi_features.  Reduces potential for generation of
         huge numbers of features (could be as many as one per site) */
      collapse_cats(phmm->cm, states);

    if (!quiet) fprintf(stderr, "Running Viterbi algorithm...\n");
    predictions = phmm_predict_viterbi_cats(phmm, states, "myseq", NULL, 
                                            "phastCons_predicted");
                                /* FIXME: need to propagate chromosome
                                   name -- add option */
   
    /* score predictions, if necessary */
    if (score) { 
      if (!quiet) fprintf(stderr, "Scoring predictions...\n");            
      phmm_score_predictions(phmm, predictions, states, NULL, NULL);
    }

    /* convert GFF to coord frame of reference sequence and adjust
       coords by idx_offset, if necessary  */
    if (refidx != 0 || msa->idx_offset != 0)
      msa_map_gff_coords(msa, predictions, 0, refidx, msa->idx_offset, NULL);

    /* now output predictions */
    if (gff)
      gff_print_set(viterbi_f, predictions);
    else                        /* BED format */
      gff_print_bed(viterbi_f, predictions, NULL, NULL); 
                                /* FIXME: need to group by "exon_id"
                                   or something similar */
  }

  /* posterior probs */
  if (post_probs) {
    int j, k;
    double *postprobs;
    int *missing;

    if (!quiet) fprintf(stderr, "Computing posterior probabilities...\n");

    postprobs = phmm_postprobs_states(phmm, states, &lnl);

    /* check for missing data, if necessary */
    if (no_missing) {
      missing = smalloc(msa->length * sizeof(int));
      msa_find_noaln(msa, refidx, MIN_BLOCK_SIZE, missing);
    }

    /* print to stdout */
    for (j = 0, k = 0; j < msa->length; j++) {
      if (refidx == 0 || msa_get_char(msa, refidx-1, j) != GAP_CHAR) {
        if (!no_missing || !missing[j]) 
          printf("%d\t%.4f\n", k + msa->idx_offset + 1, postprobs[j]);
        k++;
      }
    }
  }

  /* likelihood */
  if (lnl_f != NULL) {
    if (lnl > 0) {              /* may have already been computed */
      if (!quiet) fprintf(stderr, "Computing total log likelihood ...\n");
      lnl = phmm_lnl(phmm); 
    }
    fprintf(lnl_f, "lnL = %.4f\n", lnl); 
    if (dgamma_nrates > 1) fprintf(lnl_f, "(lambda = %f)\n", lambda);
  }

  if (!quiet)
    fprintf(stderr, "Done.\n");

  return 0;
}
