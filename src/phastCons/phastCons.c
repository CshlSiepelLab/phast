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

/* default value of lambda, used with --rates-cross */
#define DEFAULT_LAMBDA 0.9
/* default values of p and q, used with --rates-cut */
#define DEFAULT_P 0.01
#define DEFAULT_Q 0.01

void usage(char *prog) {
  printf("\n\
PROGRAM: %s\n\
\n\
DESCRIPTION: \n\
\n\
    Identify conserved elements or produce conservation scores\n\
    (posterior probabilities of conserved states), given a multiple\n\
    alignment and a phylo-HMM.  The state-transition structure of the\n\
    phylo-HMM may be specified explicitly (see --hmm) or implicitly\n\
    (see --rates-cross and --rates-cut).  Phylogenetic models (*.mod\n\
    files) must be provided in either case (they can be produced with\n\
    'phyloFit').  By default, the posterior probability of the\n\
    selected states in the HMM (see --states) is written to stdout in\n\
    a simple two-column format (position and probability, tab\n\
    separated).  Prediction of discrete elements using the Viterbi\n\
    algorithm is available via the --viterbi option.  Predicted\n\
    elements may optionally be assigned log-odds scores based on the\n\
    model (see --score).\n\
\n\
    This program is written in a general enough way that it can be\n\
    useful for other things besides analyzing rates of substitution\n\
    and evolutionary conservation.  E.g., by giving it a simple\n\
    gene-finding phylo-HMM (e.g., with a state for non-coding regions\n\
    and three states for the three codon positions), and specifying\n\
    the coding states with --states, you can obtain a \"coding\n\
    potential\" score, or even a set of crude exon predictions (using\n\
    --viterbi).  (The program 'exoniphy' is a more full-featured exon\n\
    predictor.)\n\
\n\
USAGE: %s [OPTIONS] <msa_fname> <mod_fname_list>\n\
\n\
    <mod_fname_list> must be a comma- or whitespace- delimited list of\n\
    *.mod files, as produced by 'phyloFit'; <msa_fname> should be the\n\
    name of the multiple alignment file, which can use any of several\n\
    possible file formats (see --msa-format).\n\
\n\
EXAMPLES:\n\
\n\
    (coming soon)\n\
\n\
OPTIONS:\n\
\n\
    --states, -S <state_list>\n\
        States of interest in the phylo-HMM, specified by number\n\
        (indexing starts with 1).  Default value is 1.\n\
        This option defines the meaning of the program's output.  For\n\
        example, --states \"1,2,3\" causes output of the sum of the\n\
        posterior probabilities for states 1, 2, and 3, and/or of\n\
        regions in which the Viterbi path coincides with (any of)\n\
        states 1, 2, or 3 (see --viterbi).\n\
\n\
    --hmm, -H <hmm_fname>\n\
        Name of HMM file, explicitly defining the probabilities of all\n\
        state transitions.  States in the file must correspond in\n\
        number and order to phylogenetic models in <mod_fname_list>.\n\
        Expected file format is as produced by 'hmm_train.'\n\
\n\
    --viterbi, -V <fname>\n\
        Compute the Viterbi (maximum likelihood) path and write\n\
        start and end coordinates of predicted elements to specified\n\
        file.  Output is in BED format, unless <fname> has suffix\n\
        \".gff\", in which case output is in GFF.  The features\n\
        output will represent maximal segments in which the Viterbi\n\
        path remains in the selected set of states (see --states).\n\
\n\
    --lnl, -L <fname>\n\
        Compute total log likelihood using the forward algorithm and\n\
        write it to the specified file.\n\
\n\
    --no-post-probs, -n\n\
        Suppress output of posterior probabilities.  Useful if only\n\
        Viterbi path or likelihood is of interest (saves computation\n\
        time and disk space).\n\
\n\
    --msa-format, -i PHYLIP|FASTA|MPM|SS|MAF\n\
        (Default SS) Alignment file format.\n\
\n\
    --rates-cross, -X \n\
        (Alternative to --hmm; specify only one *.mod file with this\n\
        option) Create and use an HMM with a state for every rate\n\
        category in the given phylogenetic model, and transition\n\
        probabilities defined by an autocorrelation parameter lambda\n\
        (as described by Felsenstein and Churchill, 1996).  A rate\n\
        constant for each state (rate category) will be multiplied by\n\
        the branch lengths of the phylogenetic model, to create a\n\
        \"scaled\" version of the model for that state.  If the\n\
        phylogenetic model was estimated using Yang's discrete gamma\n\
        method (-k option to phyloFit), then the rate constants will\n\
        be defined according to the estimated shape parameter 'alpha',\n\
        as described by Yang (1994).  Otherwise, a nonparameteric\n\
        model of rate variation must have been used (-K option to\n\
        phyloFit), and the rate constants will be as defined\n\
        (explicitly) in the *.mod file.  By default, the parameter\n\
        lambda will be estimated by maximum likelihood (see --lambda).\n\
\n\
    --lambda, -l [~]<lambda>\n\
        (Optionally use with --rates-cross) Fix lambda at the\n\
        specified value rather than estimating it by maximum\n\
        likelihood.  Alternatively, if first character is '~',\n\
        estimate but initialize at specified value.  Allowable range\n\
        is 0-1.  With k rate categories, the transition probability\n\
        between state i and state j will be lambda * I(i == j) +\n\
        (1-lambda)/k, where I is the indicator function.  Thus, lambda\n\
        = 0 implies no autocorrelation and lambda = 1 implies perfect\n\
        autocorrelation.\n\
\n\
    --rates-cut, -c <cut_idx>\n\
        (Alternative to --hmm and --rates-cross; specify only one\n\
        phylogenetic model) Define a simple HMM with a conserved state\n\
        (state 1) and a non-conserved state (state 2), such that the\n\
        conserved state is defined by rate categories 1-<cut_idx> from\n\
        the specified *.mod file, and the non-conserved state is\n\
        defined by the remaining rate categories.  The *.mod file must\n\
        allow for rate variation, via either the discrete gamma \n\
        (-k option to phyloFit) or non-parametric (-K option) method.\n\
        The new phylogenetic model associated with each state will be\n\
        a mixture model of rates, whose (unnormalized) mixing\n\
        proportions are given by the *.mod file, either implicitly\n\
        (discrete gamma case -- mixing proportions uniform) or explicitly\n\
         (non-parameteric case).  The transition probabilities of the HMM\n\
        will be estimated by maximum likelihood using an EM algorithm.\n\
\n\
    --cut-params, -p [~]<p>,<q>\n\
        (Optionally use with --rates-cut) Fix the transition\n\
        probabilities of the HMM at the specified values, rather than\n\
        estimating them by maximum likelihood.  Alternatively, if\n\
        first character of argument is '~', estimate parameters, but\n\
        initialize with specified values.  Here, <p> is the\n\
        probability of transitioning from the conserved to the\n\
        non-conserved state, and <q> is the probability of the reverse\n\
        transition (probabilities of self transitions are thus 1-<p>\n\
        and 1-<q>).\n\
\n\
    --nrates, -k <nrates>\n\
        (Optionally use with --rates-cross or --rates-cut and a\n\
        discrete-gamma model) Assume the specified number of rate\n\
        categories, instead of the number given in the *.mod file.\n\
        The shape parameter 'alpha' will be as given in the *.mod file.\n\
\n\
    --log, -g <log_fname>\n\
        (Optionally use with --rates-cross or --rates-cut) Write log\n\
        of optimization procedure to specified file.\n\
\n\
    --refidx, -r <refseq_idx> \n\
        Use coordinate frame of specified sequence (the value 0\n\
        indicates the frame of the entire multiple alignment).\n\
        Default value is 1 (first sequence assumed reference).\n\
\n\
    --suppress-missing, -x \n\
        Suppress posterior probabilities where cross-species alignment\n\
        data does not appear to be available.  The heuristic used is\n\
        to consider blocks of %d or more sites in which only the\n\
        reference sequence is present (as determined by --refidx) to\n\
        be regions of \"missing data\".\n\
\n\
    --reflect-strand, -U <pivot_states>\n\
        (Optionally use with --hmm) Given an HMM describing the\n\
        forward strand, create a larger HMM that allows for features\n\
        on both strands by \"reflecting\" the original HMM about the\n\
        specified \"pivot\" states.  The new HMM will be used for\n\
        prediction on both strands.  For example, suppose you use -hmm\n\
        to specify a four-state gene-finding HMM (state 1 noncoding,\n\
        states 2, 3, and 4 coding).  You can use --reflect-strand 1 to\n\
        dynamically create a new HMM with two versions of states 2, 3,\n\
        and 4, one for each strand, and with appropriately defined\n\
        transition probabilities.  Then you can use --states 2,3,4 to\n\
        obtain posterior probabilities and/or a Viterbi path\n\
        describing coding regions on either strand.\n\
\n\
    --seqname, -N <name>\n\
        (Optionally use with --viterbi) Use specified string for\n\
        'seqname' (GFF) or 'chrom' field in output file.  By default,\n\
        the filename root of <msa_fname> will be used.\n\
\n\
    --score, -s\n\
        (Optionally use with --viterbi) Assign a log-odds score to\n\
        each predictions, equal to the log total probability of the\n\
        region in question under the portion of the model defined by\n\
        --states minus its log total probability under the remainder\n\
        of the model.  To compute these scores, the states of the\n\
        model are partitioned and transition probabilities are\n\
        re-normalized.\n\
\n\
    --quiet, -q \n\
        Proceed quietly (without updates to stderr).\n\
\n\
    --help, -h\n\
        Print this help message.\n\
\n\
\n\
REFERENCES:\n\
\n\
    J. Felsenstein and G. Churchill.  1996. A hidden Markov model\n\
      approach to variation among sites in rate of evolution. \n\
      Mol. Biol. Evol., 13:93-104. \n\
\n\
    A. Siepel and D. Haussler.  2003.  Combining phylogenetic and\n\
      hidden Markov models in biosequence analysis.  RECOMB '03.\n\
\n\
    Z. Yang. 1994. Maximum likelihood phylogenetic estimation from\n\
      DNA sequences with variable rates over sites: approximate\n\
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
    gff = FALSE, rates_cross = FALSE, estim_lambda = TRUE, 
    estim_cut_params = TRUE;
  int nrates = -1, rates_cut_idx = -1, refidx = 1;
  double lambda = DEFAULT_LAMBDA, p = DEFAULT_P, q = DEFAULT_Q;
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
    {"quiet", 0, 0, 'q'},
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
      if (optarg[0] != '~') estim_lambda = FALSE;
      else optarg = &optarg[1];
      lambda = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'c':
      rates_cut_idx = get_arg_int_bounds(optarg, 1, 100);
      break;
    case 'p':
      if (optarg[0] != '~') estim_cut_params = FALSE;
      else optarg = &optarg[1];
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

  /* estimate lambda, if necessary */
  if (rates_cross && estim_lambda) {
    if (!quiet) fprintf(stderr, "Finding MLE for lambda ...");
    lnl = phmm_fit_lambda(phmm, &lambda, log_f);
    if (!quiet) fprintf(stderr, " (lambda = %f)\n", lambda);
    phmm_update_cross_prod(phmm, lambda);
  }

  /* estimate p and q, if necessary */
  else if (rates_cut_idx != -1 && estim_cut_params) {
    if (!quiet) fprintf(stderr, "Finding MLE for 'p' and 'q'...");
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
