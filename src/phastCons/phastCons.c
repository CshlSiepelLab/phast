#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <misc.h>
#include <phylo_hmm.h>
#include <sufficient_stats.h>
#include <bed.h>
#include <dgamma.h>

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
             --rates-cross and --rates-cut).  Phylogenetic models
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

    --rates-cross, -X 
        (Alternative to --hmm; specify only one *.mod file with this
        option) Create and use an HMM with a state for every rate
        category in the given phylogenetic model, and transition
        probabilities defined by an autocorrelation parameter lambda
        (as described by Felsenstein and Churchill, 1996).  A rate
        constant for each state (rate category) will be multiplied by
        the branch lengths of the phylogenetic model, to create a
        \"scaled\" version of the model for that state.  If the
        phylogenetic model was estimated using Yang's discrete gamma
        method (-k option to phyloFit), then the rate constants will
        be defined according to the estimated shape parameter 'alpha',
        as described by Yang (1994).  Otherwise, a nonparameteric
        model of rate variation must have been used (-K option to
        phyloFit), and the rate constants will be as defined
        (explicitly) in the *.mod file.  By default, the parameter
        lambda will be estimated by maximum likelihood (see --lambda).

    --lambda, -l <lambda>
        (Optionally use with --rates-cross) Fix lambda at the
        specified value rather than estimating it by maximum likelihood.
        Allowable range is 0-1.  With k rate categories, the
        transition probability between state i and state j will be
        lambda * I(i == j) + (1-lambda)/k, where I is the indicator
        function.  Thus, lambda = 0 implies no autocorrelation and
        lambda = 1 implies perfect autocorrelation.

    --rates-cut, -c <cut_idx>
        (Alternative to --hmm and --rates-cross; specify only one
        phylogenetic model) Define a simple HMM with a conserved state
        (state 1) and a non-conserved state (state 2), such that the
        conserved state is defined by rate categories 1-<cut_idx> from
        the specified *.mod file, and the non-conserved state is
        defined by the remaining rate categories.  The *.mod file must
        allow for rate variation, via either the discrete gamma 
        (-k option to phyloFit) or non-parametric (-K option) method.
        The new phylogenetic model associated with each state will be
        a mixture model of rates, whose (unnormalized) mixing
        proportions are given by the *.mod file, either implicitly
        (discrete gamma case) or explicitly (non-parameteric case).
        The transition probabilities of the HMM will be estimated by
        maximum likelihood.

    --cut-params, -p <p>,<q>
        (Optionally use with --rates-cut) Fix the transition
        probabilities of the HMM at the specified values, rather than
        estimating them by maximum likelihood.  Here, <p> is the
        probability of transitioning from the conserved to the
        non-conserved state, and <q> is the probability of the reverse
        transition (probabilities of self transitions are thus 1-<p>
        and 1-<q>).

    --nrates, -k <nrates>
        (Optionally use with --rates-cross or --rates-cut and a
        discrete-gamma model) Assume the specified number of rate
        categories, instead of the number given in the *.mod file.
        The shape parameter 'alpha' will be as given in the *.mod file.

    --log, -g <log_fname>
        (Optionally use with --rates-cross or --rates-cut) Write log
        of optimization procedure to specified file.

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

    --seqname, -N <name>
        (Optionally use with --viterbi) Use specified string for
        'seqname' (GFF) or 'chrom' field in output file.  By default,
        the filename root of <msa_fname> will be used.

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
    gff = FALSE, rates_cross = FALSE;
  int nrates = -1, rates_cut_idx = -1, refidx = 1;
  double lambda = -1, p = -1, q = -1;
  msa_format_type msa_format = SS;
  FILE *viterbi_f = NULL, *lnl_f = NULL, *log_f = NULL;
  List *states = NULL, *pivot_states = NULL;
  char *seqname = NULL;
  HMM *hmm = NULL;

  struct option long_opts[] = {
    {"states", 1, 0, 'S'},
    {"hmm", 1, 0, 'H'},
    {"viterbi", 1, 0, 'V'},
    {"no-post-probs", 0, 0, 'n'},
    {"msa-format", 1, 0, 'i'},
    {"rates-cross", 0, 0, 'X'},
    {"lambda", 1, 0, 'l'},
    {"rates-cut", 1, 0, 'c'},
    {"cut-params", 1, 0, 'p'},
    {"nrates", 1, 0, 'k'},
    {"log", 1, 0, 'g'},
    {"refidx", 1, 0, 'r'},
    {"suppress-missing", 0, 0, 'x'},
    {"reflect-strand", 0, 0, 'U'},
    {"lnl", 1, 0, 'L'},
    {"seqname", 1, 0, 'N'},
    {"score", 0, 0, 's'},
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

  while ((c = getopt_long(argc, argv, "S:H:V:ni:k:l:c:p:r:xL:s:N:g:Xqh", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'S':
      states = get_arg_list(optarg);
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
    case 'X':
      rates_cross = TRUE;
      break;
    case 'l':
      lambda = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'c':
      rates_cut_idx = get_arg_int_bounds(optarg, 1, 100);
      break;
    case 'p':
      tmpl = get_arg_list_dbl(optarg);
      if (lst_size(tmpl) != 2) die("ERROR: bad argument to --cut-params.\n");
      p = lst_get_dbl(tmpl, 0);
      q = lst_get_dbl(tmpl, 1);
      if (p <= 0 || p >= 1 || q <= 0 || q >= 1)
        die("ERROR: bad argument to --cut-params.\n");
      lst_free(tmpl);
      break;
    case 'k':
      nrates = get_arg_int_bounds(optarg, 2, 100);
      break;
    case 'g':
      log_f = fopen_fname(optarg, "w+");
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
    case 'N':
      seqname = optarg;
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

  if ((hmm != NULL && rates_cross) ||
      (hmm != NULL && rates_cut_idx != -1) ||
      (rates_cross && rates_cut_idx != -1))
    die("ERROR: --hmm, --rates-cross, and --rates-cut are mutually exclusive.\n");
  
  if (hmm == NULL && !rates_cross && rates_cut_idx == -1)
    die("ERROR: must specify one of --hmm, --rates-cross, and --rates-cut.\n");

  if (optind != argc - 2) 
    die("ERROR: missing required arguments.  Try '%s -h'.\n", argv[0]);

  /* read tree models first (alignment may take a while) */
  tmpl = get_arg_list(argv[optind+1]);

  if ((rates_cross || rates_cut_idx != -1) && lst_size(tmpl) != 1)
    die("ERROR: only one tree model allowed with --rates-cross or --rates-cut.\n");
  else if (hmm != NULL && hmm->nstates != lst_size(tmpl)) 
    die("ERROR: number of states in HMM must equal number of tree models.\n");
    
  mod = (TreeModel**)smalloc(sizeof(TreeModel*) * lst_size(tmpl));
  for (i = 0; i < lst_size(tmpl); i++) {
    String *fname = lst_get_ptr(tmpl, i);
    if (!quiet)
      fprintf(stderr, "Reading tree model from %s...\n", fname->chars);
    mod[i] = tm_new_from_file(fopen_fname(fname->chars, "r"));
    mod[i]->use_conditionals = 1; /* FIXME: necessary? */
  }

  /* check rates-cross and rates-cut options vis-a-vis the tree model */
  if (rates_cross || rates_cut_idx) {
    if (mod[0]->nratecats <= 1)
      die("ERROR: --rates-cross and --rates-cut require tree model allowing for rate variation.\n");
    if (nrates != -1 && mod[0]->empirical_rates)
      die("ERROR: can't use --nrates with nonparameteric rate model.\n");
    if (rates_cut_idx != -1 && rates_cut_idx > mod[0]->nratecats)
      die("ERROR: --rates-cut arg must be <= NRATECATS from *.mod file.\n");
    if (nrates == -1) nrates = mod[0]->nratecats;
  }

  /* read alignment */
  if (!quiet)
    fprintf(stderr, "Reading alignment from %s...\n", argv[optind]);
  msa = msa_new_from_file(fopen_fname(argv[optind], "r"), msa_format, NULL);
  /* use file name root for default seqname */

  if (seqname == NULL) {
    String *tmp = str_new_charstr(argv[optind]);
    str_remove_path(tmp);
    str_root(tmp, '.');
    seqname = tmp->chars;
  }

  /* tweak alphabet, if necessary */
  if (msa_format == SS) {
    if (msa->ss->tuple_idx == NULL) 
      die("ERROR: Ordered representation of alignment required.\n");
    msa_remove_N_from_alph(msa);
  }

  /* set up states */
  if (states == NULL) {
    states = lst_new_ptr(1);
    lst_push_ptr(states, str_new_charstr("0"));
  }
  else {
    for (i = 0; i < lst_size(states); i++) {
      int state_no;
      String *state = lst_get_ptr(states, i);
      if (str_as_int(state, &state_no) != 0 || state_no < 1) 
        die("ERROR: illegal state '%s'.\n", state->chars);
      str_clear(state);
      str_append_int(state, state_no-1); /* internally use 0-based
                                            indexing */
    }
  }

  /* set up PhyloHmm */
  phmm = phmm_new(hmm, mod, NULL, pivot_states, NULL, -1);

  if (rates_cross) {
    if (!quiet) 
      fprintf(stderr, "Creating %d scaled versions of tree model...\n", nrates);
    phmm_rates_cross(phmm, nrates, lambda, TRUE);
  }

  else if (rates_cut_idx != -1) {
    if (!quiet) 
      fprintf(stderr, "Partitioning at rate category %d to create 'conserved' and 'nonconserved' states...\n", rates_cut_idx);

    phmm_rates_cut(phmm, nrates, rates_cut_idx, p, q);
  }

  /* compute emissions */
  phmm_compute_emissions(phmm, msa, quiet);

  /* fit lambda, if necessary */
  if (rates_cross && lambda == -1) {
    if (!quiet) fprintf(stderr, "Finding MLE for lambda ...");
    lnl = phmm_fit_lambda(phmm, &lambda, log_f);
    if (!quiet) fprintf(stderr, " (lambda = %f)\n", lambda);
    phmm_update_cross_prod(phmm, lambda);
  }

  /* fit p and q, if necessary */
  else if (rates_cut_idx != -1 && (p == -1 || q == -1)) {
    if (!quiet) fprintf(stderr, "Finding MLE for 'p' and 'q'...");
    p = 0.01; q = 0.001;        /* initial values */
    lnl = phmm_fit_rates_cut(phmm, &p, &q, log_f);
    if (!quiet) fprintf(stderr, " (p = %f. q = %f)\n", p, q);
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
    predictions = phmm_predict_viterbi_cats(phmm, states, seqname, NULL,
                                            "phastCons_predicted");
    /* note that selected state numbers are also cat numbers  */
   
    /* score predictions, if necessary */
    if (score) { 
      if (!quiet) fprintf(stderr, "Scoring predictions...\n");            
      phmm_score_predictions(phmm, predictions, states, NULL, NULL, FALSE);
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
  }

  /* posterior probs */
  if (post_probs) {
    int j, k;
    double *postprobs;
    int *missing;

    if (!quiet) fprintf(stderr, "Computing posterior probabilities...\n");

    postprobs = phmm_postprobs_cats(phmm, states, &lnl);
    /* note that selected state numbers are also cat numbers  */

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
    if (rates_cross) fprintf(lnl_f, "(lambda = %f)\n", lambda);
    else if (rates_cut_idx != -1) fprintf(lnl_f, "(p = %f, q = %f)\n", p, q);
  }

  if (!quiet)
    fprintf(stderr, "Done.\n");

  return 0;
}
