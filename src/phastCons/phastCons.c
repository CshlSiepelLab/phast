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
DESCRIPTION:\n\
\n\
    Identify conserved elements or produce conservation scores, given\n\
    a multiple alignment and a phylo-HMM.  By default, a phylo-HMM\n\
    consisting of a single conserved state and a single non-conserved\n\
    state is assumed, and its transition probabilities are estimated\n\
    from the data -- i.e., no HMM structure needs to be defined by the\n\
    user.  (The transition probabilities can be set a priori if\n\
    desired; see --transitions.)  It is also possible to use a k-state\n\
    phylo-HMM of the kind described by Felsenstein and Churchill\n\
    (1996) (see --rates-cross), or to define the state-transition\n\
    structure of the HMM explicitly (see --hmm).\n\
\n\
    In the (default) two-state HMM and --rates-cross cases, a single\n\
    phylogenetic model (*.mod file) must be given, and this model must\n\
    allow for variation in evolutionary rate (it can be produced using\n\
    the -k/-K options to phyloFit).  In the --hmm case, a model file\n\
    must be given for each state in the HMM.  These files must be\n\
    given in the same order as the states in the HMM file.\n\
\n\
    By default, the program computes the posterior probability at each\n\
    site of the *first* state in the HMM.  In the cases of the default\n\
    two-state HMM and the Felsenstein-Churchill HMM, this is the most\n\
    conserved state, and these probabilities can be interpreted as\n\
    conservation scores.  They are written to stdout in a simple\n\
    tab-separated two-column format (position and probability).  The\n\
    set of states whose total (marginal) posterior probability is\n\
    reported can be changed using the --states option.  In addition,\n\
    discrete elements can be predicted using the --viterbi option, and\n\
    they can be assigned log-odds scores using the --score option.\n\
    The set of states considered when predicting discrete elements is\n\
    also controlled by --states.\n\
\n\
    This program is written in a general enough way that it can be\n\
    useful for other things besides analyzing rates of substitution\n\
    and evolutionary conservation.  E.g., by giving it a simple\n\
    gene-finding phylo-HMM (e.g., with a state for non-coding regions\n\
    and three states for the three codon positions), and specifying\n\
    the coding states via --states, you can obtain posterior\n\
    probabilities that can be interpreted as a measure of \"coding\n\
    potential.\"\n\
\n\
USAGE: %s [OPTIONS] <msa_fname> <mod_fname_list>\n\
\n\
    <mod_fname_list> must be a comma-delimited list of *.mod files, as\n\
    produced by 'phyloFit'; <msa_fname> should be the name of the\n\
    multiple alignment file, which can use any of several possible\n\
    file formats (see --msa-format).\n\
\n\
EXAMPLES:\n\
\n\
    1. Fit a phylogenetic model to a data set, using the discrete\n\
    gamma model for rate variation, then produce conservation scores\n\
    via a two-state phylo-HMM.  Use an input file in SS format (see\n\
    msa_view).\n\
\n\
        phyloFit --tree mytree.nh --subst-mod REV --nrates 5 \\\n\
            mydata.ss --msa-format SS --out-root rev-dg\n\
        phastCons mydata.ss rev-dg.mod --nrates 20 > cons.dat\n\
\n\
    (The --nrates 20 option to phastCons causes the estimated\n\
    distribution of evolutionary rates to be partitioned into 20 parts\n\
    of equal probability mass; as a result, the partition having\n\
    the smallest rate [by default, the one used for the conserved\n\
    state of the HMM] is representative of the most conserved 1 / 20 =\n\
    5%% of sites in the data set.  Thus, the reported scores represent\n\
    the posterior probability that each site is among the 5%% most\n\
    conserved sites in the data set, modulo the \"smoothing\" imposed\n\
    by the HMM.)\n\
\n\
    2. As in (1), but also predict discrete conserved elements and\n\
    score them using log-odds scoring.  Write predictions in BED format\n\
    to elements.bed.\n\
\n\
        phastCons mydata.ss rev-dg.mod --nrates 20 --viterbi elements.bed \\\n\
            --score > cons.dat\n\
\n\
    (if output file were \"elements.gff,\" then output would be in GFF\n\
    instead)\n\
\n\
    3. As in (1), but bypass finding the maximum likelihood estimate of\n\
    the state transition probabilities, and instead just use the values\n\
    given (much faster).\n\
\n\
        phastCons mydata.ss rev-dg.mod --transitions 0.005,0.001 \\\n\
            --nrates 20 > cons.dat\n\
\n\
    4. As in (1), but this time use a k-state HMM with transition\n\
    probabilities defined by an autocorrelation parameter lambda, as\n\
    described by Felsenstein and Churchill (1996).  Use k = 10 states\n\
    and estimate the parameter lambda from the data.  Report posterior\n\
    probabilities not just of the single most conserved state, but of\n\
    the two most conserved states.\n\
\n\
        phastCons mydata.ss rev-dg.mod --rates-cross --nrates 10 \\\n\
            --states 1,2 > cons.dat\n\
\n\
    5. As in (4), but fix lambda at 0.9 rather than estimating it from\n\
    the data.\n\
\n\
        phastCons mydata.ss rev-dg.mod --rates-cross --nrates 10 \\\n\
            --states 1,2 --lambda 0.9 > cons.dat\n\
\n\
    6. Compute a coding potential score, using a simple gene-finding\n\
    HMM and phylogenetic models for the three codon positions and\n\
    non-coding sites.  Allow for genes on either strand.\n\
\n\
        phastCons mydata.ss noncoding.mod,codon1.mod,codon2.mod,codon3.mod \\\n\
            --hmm simple-4state.hmm --reflect-strand 1 --states 2,3,4 \\\n\
            > coding-potential.dat\n\
\n\
\n\
OPTIONS:\n\
\n\
 (HMM structure and transition probabilities)\n\
    --hmm, -H <hmm_fname>\n\
        Name of HMM file, explicitly defining the probabilities of all\n\
        state transitions.  States in the file must correspond in\n\
        number and order to phylogenetic models in <mod_fname_list>.\n\
        Expected file format is as produced by 'hmm_train.'\n\
\n\
    --reflect-strand, -U <pivot_states>\n\
        (Optionally use with --hmm) Given an HMM describing the\n\
        forward strand, create a larger HMM that allows for features\n\
        on both strands by \"reflecting\" the original HMM about the\n\
        specified \"pivot\" states.  The new HMM will be used for\n\
        prediction on both strands.\n\
\n\
    --rates-cross, -X\n\
        (Alternative to --hmm; specify only one *.mod file with this\n\
        option) Use an HMM with a state for every rate\n\
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
    --cut-at, -c <cut_idx>\n\
        (For use with default two-state HMM) Use rate categories\n\
        1-<cut_idx> for the conserved state (state 1) and the\n\
        remaining rate categories for the non-conserved state.  The\n\
        given phylogenetic model must allow for rate variation, via\n\
        either the discrete gamma (-k option to phyloFit) or\n\
        non-parametric (-K option) method.  The new phylogenetic model\n\
        associated with each state will be a mixture model of rates,\n\
        whose (unnormalized) mixing proportions are given by the\n\
        original *.mod file, either implicitly (discrete gamma case;\n\
        mixing proportions uniform) or explicitly (non-parameteric\n\
        case).  By default, the transition probabilities of the HMM\n\
        will be estimated by maximum likelihood using an EM algorithm\n\
        (see --transitions).\n\
\n\
    --transitions, -p [~]<p>,<q>\n\
        (Optionally use with default two-state HMM) Fix the transition\n\
        probabilities of the two-state HMM as specified, rather than\n\
        estimating them by maximum likelihood.  Alternatively, if\n\
        first character of argument is '~', estimate parameters, but\n\
        initialize with specified values.  The argument <p> is the\n\
        probability of transitioning from the conserved to the\n\
        non-conserved state, and <q> is the probability of the reverse\n\
        transition (probabilities of self transitions are thus 1-<p>\n\
        and 1-<q>).\n\
\n\
    --nrates, -k <nrates>\n\
        (Optionally use with a discrete-gamma model) Assume the\n\
        specified number of rate categories, instead of the number\n\
        given in the *.mod file.  The shape parameter 'alpha' will be\n\
        as given in the *.mod file.\n\
\n\
 (Output)\n\
    --states, -S <state_list>\n\
        States of interest in the phylo-HMM, specified by number\n\
        (indexing starts with 1).  Default value is 1.\n\
        Choosing --states \"1,2,3\" will cause output of the sum of the\n\
        posterior probabilities for states 1, 2, and 3, and/or of\n\
        regions in which the Viterbi path coincides with (any of)\n\
        states 1, 2, or 3 (see --viterbi).\n\
\n\
    --viterbi, -V <fname>\n\
        Compute the Viterbi (maximum likelihood) path and write\n\
        start and end coordinates of predicted elements to specified\n\
        file.  Output is in BED format, unless <fname> has suffix\n\
        \".gff\", in which case output is in GFF.  The features\n\
        output will represent maximal segments in which the Viterbi\n\
        path remains in the selected set of states (see --states).\n\
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
    --lnl, -L <fname>\n\
        Compute total log likelihood using the forward algorithm and\n\
        write it to the specified file.\n\
\n\
    --no-post-probs, -n\n\
        Suppress output of posterior probabilities.  Useful if only\n\
        Viterbi path or likelihood is of interest (saves computation\n\
        time and disk space).\n\
\n\
    --suppress-missing, -x\n\
        Suppress posterior probabilities where cross-species alignment\n\
        data does not appear to be available.  The heuristic used is\n\
        to consider blocks of %d or more sites in which only the\n\
        reference sequence is present (as determined by --refidx) to\n\
        be regions of \"missing data\".\n\
\n\
    --log, -g <log_fname>\n\
        (Optionally use when estimating transition probabilities)\n\
        Write log of optimization procedure to specified file.\n\
\n\
    --refidx, -r <refseq_idx>\n\
        Use coordinate frame of specified sequence (the value 0\n\
        indicates the frame of the entire multiple alignment) in\n\
        output.  Default value is 1 (first sequence assumed\n\
        reference).\n\
\n\
    --seqname, -N <name>\n\
        (Optionally use with --viterbi) Use specified string for\n\
        'seqname' (GFF) or 'chrom' field in output file.  By default,\n\
        the filename root of <msa_fname> will be used.\n\
\n\
 (Other)\n\
    --msa-format, -i PHYLIP|FASTA|MPM|SS|MAF\n\
        (Default SS) Alignment file format.\n\
\n\
    --quiet, -q\n\
        Proceed quietly (without updates to stderr).\n\
\n\
    --help, -h\n\
        Print this help message.\n\
\n\
\n\
REFERENCES:\n\
\n\
    J. Felsenstein and G. Churchill.  1996. A hidden Markov model\n\
      approach to variation among sites in rate of evolution.\n\
      Mol. Biol. Evol., 13:93-104.\n\
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
    estim_transitions = TRUE, two_state = TRUE;
  int nrates = -1, rates_cut_idx = 1, refidx = 1;
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
    {"cut-at", 1, 0, 'c'},
    {"transitions", 1, 0, 't'},
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

  while ((c = getopt_long(argc, argv, "S:H:V:ni:k:l:c:t:r:xL:s:N:g:Xqh", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'S':
      states = get_arg_list(optarg);
      break;
    case 'H':
      if (!quiet) fprintf(stderr, "Reading HMM from %s...", optarg);
      hmm = hmm_new_from_file(fopen_fname(optarg, "r"));
      two_state = FALSE;
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
      two_state = FALSE;
      break;
    case 'l':
      if (optarg[0] != '~') estim_lambda = FALSE;
      else optarg = &optarg[1];
      lambda = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'c':
      rates_cut_idx = get_arg_int_bounds(optarg, 1, 100);
      break;
    case 't':
      if (optarg[0] != '~') estim_transitions = FALSE;
      else optarg = &optarg[1];
      tmpl = get_arg_list_dbl(optarg);
      if (lst_size(tmpl) != 2) die("ERROR: bad argument to --transitions.\n");
      p = lst_get_dbl(tmpl, 0);
      q = lst_get_dbl(tmpl, 1);
      if (p <= 0 || p >= 1 || q <= 0 || q >= 1)
        die("ERROR: bad argument to --transitions.\n");
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

  if ((hmm != NULL && rates_cross))
    die("ERROR: --hmm and --rates-cross are mutually exclusive.\n");
  
  if (optind != argc - 2) 
    die("ERROR: missing required arguments.  Try '%s -h'.\n", argv[0]);

  /* read tree models first (alignment may take a while) */
  tmpl = get_arg_list(argv[optind+1]);

  if ((rates_cross || two_state) && lst_size(tmpl) != 1)
    die("ERROR: only one tree model allowed with --rates-cross or --cut-at.\n");
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

  /* check rates-cross and cut-at options vis-a-vis the tree model */
  if (rates_cross || two_state) {
    if (mod[0]->nratecats <= 1)
      die("ERROR: a tree model allowing for rate variation is required.\n");
    if (nrates != -1 && mod[0]->empirical_rates)
      die("ERROR: can't use --nrates with nonparameteric rate model.\n");
    if (nrates == -1) nrates = mod[0]->nratecats;
    if (rates_cut_idx > nrates)
      die("ERROR: --cut-at arg must be <= nrates.\n");
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

  else if (two_state) {
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
  else if (two_state && estim_transitions) {
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
      gff_print_bed(viterbi_f, predictions, FALSE); 
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
    else if (two_state) fprintf(lnl_f, "(p = %f, q = %f)\n", p, q);
  }

  if (!quiet)
    fprintf(stderr, "Done.\n");

  return 0;
}
