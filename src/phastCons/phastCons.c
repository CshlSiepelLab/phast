#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <misc.h>
#include <phylo_hmm.h>
#include <sufficient_stats.h>
#include <bed.h>
#include <dgamma.h>

/* default value of lambda, used with --rates-cross */
#define DEFAULT_LAMBDA 0.9
/* default values of p and q, used with two-state model */
#define DEFAULT_P 0.01
#define DEFAULT_Q 0.01

void usage(char *prog) {
  printf("\n\
PROGRAM: %s\n\
\n\
USAGE: %s [OPTIONS] <msa_fname> <mod_fname_list>\n\
\n\
    <msa_fname> must be the name of a multiple alignment file,\n\
    which can use any of several possible file formats (see\n\
    --msa-format); <mod_fname_list> must be a comma-delimited list of\n\
    *.mod files, as produced by 'phyloFit.'\n\
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
    site of the *first* state (index 0) in the HMM.  In the cases of the\n\
    default two-state HMM and the Felsenstein-Churchill HMM, this is the\n\
    most conserved state, and these probabilities can be interpreted as\n\
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
    potential (see --coding-potential).\"\n\
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
            --states 0,1 > cons.dat\n\
\n\
    5. As in (4), but fix lambda at 0.9 rather than estimating it from\n\
    the data.\n\
\n\
        phastCons mydata.ss rev-dg.mod --rates-cross --nrates 10 \\\n\
            --states 0,1 --lambda 0.9 > cons.dat\n\
\n\
    6. Compute a coding potential score, using a simple gene-finding\n\
    HMM and phylogenetic models for the three codon positions and\n\
    non-coding sites.  Allow for genes on either strand.\n\
\n\
        phastCons mydata.ss noncoding.mod,codon1.mod,codon2.mod,codon3.mod \\\n\
            --hmm simple-4state.hmm --reflect-strand 0 --states 1,2,3 \\\n\
            > coding-potential.dat\n\
\n\
    7. Compute a coding potential score using a default phylo-HMM,\n\
    a simplified version of the one used in exoniphy.  Allows for\n\
    conserved non-coding sequences and makes use of the different\n\
    patterns of indels seen in coding and non-coding regions.\n\
    Currently assumes a human/mouse/rat alignment.\n\
\n\
        phastCons --coding-potential human-mouse-rat.ss > cp.dat\n\
\n\
\n\
OPTIONS:\n\
\n\
 (HMM structure and transition probabilities)\n\
    --hmm, -H <hmm_fname>\n\
        Name of HMM file explicitly defining the probabilities of all\n\
        state transitions.  States in the file must correspond in\n\
        number and order to phylogenetic models in <mod_fname_list>.\n\
        Expected file format is as produced by 'hmm_train.'\n\
\n\
    --catmap, -c <fname>|<string>\n\
        (Optionally use with --hmm)  Mapping of feature types to category\n\
        numbers.  Can give either a filename or an \"inline\" description\n\
        of a simple category map, e.g., --catmap \"NCATS = 3 ; CDS 1-3\".\n\
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
        1-<cut_idx> for the conserved state (state 0) and the\n\
        remaining rate categories for the non-conserved state (default\n\
        value is 1).  The given phylogenetic model must allow for rate\n\
        variation, via either the discrete gamma (-k option to\n\
        phyloFit) or non-parametric (-K option) method.  The new\n\
        phylogenetic model associated with each state will be a\n\
        mixture model of rates, whose (unnormalized) mixing\n\
        proportions are given by the original *.mod file, either\n\
        implicitly (discrete gamma case; mixing proportions uniform)\n\
        or explicitly (non-parameteric case).  By default, the\n\
        transition probabilities of the HMM will be estimated by\n\
        maximum likelihood using an EM algorithm (see --transitions).\n\
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
 (Indels, forward/reverse strands, missing data, and coding potential)\n\
    --indels, -I\n\
        (Optionally use with --hmm) Expand HMM state space to model\n\
        indels as described in Siepel & Haussler (2004).\n\
\n\
    --reflect-strand, -U <pivot_states>\n\
        (Optionally use with --hmm) Given an HMM describing the\n\
        forward strand, create a larger HMM that allows for features\n\
        on both strands by \"reflecting\" the original HMM about the\n\
        specified \"pivot\" states.  The new HMM will be used for\n\
        prediction on both strands.  States can be specified by number\n\
        (indexing starts with 0), or if --catmap, by category name.\n\
\n\
    --min-informative-types, -M <list>\n\
        Require a minimum number of \"informative\" bases (i.e.,\n\
        non-missing-data characters) for the specified HMM states.\n\
        (Specify by state number [indexing starts with 0] or, if\n\
        --catmap, by category name.)  Columns not meeting the minimum\n\
        number (see --min-informative-bases) will be given emission\n\
        probabilities of zero.\n\
\n\
    --min-informative-bases, -m <number>\n\
        Minimum number of informative bases used with --min-informative-types\n\
        (default is 2).\n\
\n\
    --coding-potential, -p\n\
        Use parameter settings that cause output to be interpretable\n\
        as a coding potential score.  By default, a simplified version\n\
        of exoniphy's phylo-HMM is used, with a noncoding\n\
        (background) state, a conserved non-coding (CNS) state, and\n\
        states for the three codon positions.  This option implies\n\
        --catmap \"NCATS=4; CNS 1; CDS 2-4\" --hmm <default-HMM-file>\n\
        --states CDS --indels --reflect-strand background,CNS\n\
        --min-informative-types CDS, plus a set of default *.mod files.\n\
        (All of these options can be overridden.)\n\
\n\
    --indels-only, -J\n\
        Like --indels but force the use of a single-state HMM.  This\n\
        option allows the effect of the indel model in isolation to be\n\
        observed.  Implies --no-post-probs.  Use with --lnl.\n\
\n\
 (Output)\n\
    --states, -S <state_list>\n\
        States of interest in the phylo-HMM, specified by number\n\
        (indexing starts with 0), or if --catmap, by category name.\n\
        Default value is 1.  Choosing --states \"0,1,2\" will cause\n\
        output of the sum of the posterior probabilities for states 0,\n\
        1, and 2, and/or of regions in which the Viterbi path\n\
        coincides with (any of) states 0, 1, or 2 (see --viterbi).\n\
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
        Viterbi path or likelihood is of interest.\n\
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
        'seqname' (GFF) or 'chrom' field in output file.  Default\n\
        is obtained from input file name (double filename\n\
        root, e.g., \"chr22\" if input file is \"chr22.35.ss\").\n\
\n\
    --idpref, -P <name>\n\
        (Optionally use with --viterbi) Use specified string as\n\
        prefix of generated ids in output file.  Can be used to\n\
        ensure ids are unique.  Default is obtained from input\n\
        file name (single filename root, e.g., \"chr22.35\" if\n\
        input file is \"chr22.35.ss\").\n\
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
    A. Siepel and D. Haussler.  2004.  Computational identification of\n\
      evolutionarily conserved exons.  Proc. 8th Annual Int'l Conf.\n\
      on Research in Computational Biology (RECOMB '04), pp. 177-186.\n\
\n\
    Z. Yang. 1994. Maximum likelihood phylogenetic estimation from\n\
      DNA sequences with variable rates over sites: approximate\n\
      methods. J. Mol. Evol., 39:306-314.\n\n", prog, prog);

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

void setup_rates_cut(HMM **hmm, TreeModel ***mods, CategoryMap **cm, int nrates, 
                     int cut_idx, double p, double q) {

  double freq1 = 0;
  int i;
  int dgamma = !(*mods)[0]->empirical_rates; /* whether using discrete
                                                gamma model */

  List *rconsts, *rweights;
  double tmp_rates[nrates], tmp_weights[nrates];

  if (nrates <= 1) die("ERROR: must have nrates > 1.\n");
  if (cut_idx < 1 || cut_idx > nrates) 
    die("ERROR: must have 1 <= cut_idx <= nrates.\n");

  *hmm = hmm_new_nstates(2, TRUE, FALSE);

  if (dgamma) 
    /* if using dgamma, need to compute rate consts and weights -- may
       not have been computed yet */
    DiscreteGamma(tmp_weights, tmp_rates, (*mods)[0]->alpha, 
                  (*mods)[0]->alpha, nrates, 0); 
  else {
    for (i = 0; i < nrates; i++) {
      tmp_weights[i] = (*mods)[0]->freqK[i];
      tmp_rates[i] = (*mods)[0]->rK[i];
    }
  }

  /* set HMM transitions according to p and q */
  mm_set((*hmm)->transition_matrix, 0, 0, 1-p);
  mm_set((*hmm)->transition_matrix, 0, 1, p);
  mm_set((*hmm)->transition_matrix, 1, 0, q);
  mm_set((*hmm)->transition_matrix, 1, 1, 1-q);

  /* set HMM begin transitions according to weights */
  for (i = 0; i < cut_idx; i++) freq1 += tmp_weights[i];
  gsl_vector_set((*hmm)->begin_transitions, 0, freq1);
  gsl_vector_set((*hmm)->begin_transitions, 1, 1 - freq1);

  hmm_reset(*hmm);

  /* create 2nd tree model, then update rate categories in both */
  (*mods) = srealloc(*mods, 2 * sizeof(void*));
  (*mods)[1] = tm_create_copy((*mods)[0]);

  rconsts = lst_new_dbl(nrates);
  rweights = lst_new_dbl(nrates);
  for (i = 0; i < cut_idx; i++) {
    lst_push_dbl(rweights, tmp_weights[i]);
    lst_push_dbl(rconsts, tmp_rates[i]);
  }

  if (cut_idx == 1) {
    tm_reinit((*mods)[0], (*mods)[0]->subst_mod, 1, 0, NULL, NULL);
    tm_scale((*mods)[0], lst_get_dbl(rconsts, 0), TRUE);
                                /* in this case, have to by-pass rate
                                   variation machinery and just scale
                                   tree; tree model code won't do
                                   rate variation with single rate
                                   category */
  }
  else
    tm_reinit((*mods)[0], (*mods)[0]->subst_mod, cut_idx, 
              (*mods)[0]->alpha, rconsts, rweights);
                                /* note that dgamma model will be
                                   redefined as empirical rates mod */

  lst_clear(rconsts); lst_clear(rweights);
  for (i = cut_idx; i < nrates; i++) {
    lst_push_dbl(rweights, tmp_weights[i]);
    lst_push_dbl(rconsts, tmp_rates[i]);
  }

  if (cut_idx == nrates-1) {    /* unlikely but possible */
    tm_reinit((*mods)[1], (*mods)[1]->subst_mod, 1, 0, NULL, NULL);
    tm_scale((*mods)[1], lst_get_dbl(rconsts, nrates-1), TRUE);
  }
  else
    tm_reinit((*mods)[1], (*mods)[1]->subst_mod, nrates - cut_idx, 
              (*mods)[1]->subst_mod, rconsts, rweights);
  lst_free(rconsts); lst_free(rweights);

  /* define two-category category map */
  *cm = cm_create_trivial(1, "cons_");
}

/** Estimate the parameters 'p' and 'q' that define the two-state
    "rates-cut" model using an EM algorithm.  Returns ln likelihood. */
double fit_rates_cut(PhyloHmm *phmm, 
                          double *p, 
                          double *q, 
                          FILE *logf
                          ) {
  double retval;

  mm_set(phmm->functional_hmm->transition_matrix, 0, 0, 1-*p);
  mm_set(phmm->functional_hmm->transition_matrix, 0, 1, *p);
  mm_set(phmm->functional_hmm->transition_matrix, 1, 0, *q);
  mm_set(phmm->functional_hmm->transition_matrix, 1, 1, 1-*q);
                                /* note that phmm->functional_hmm ==
                                   phmm->hmm if no indel model */
  phmm_reset(phmm); 

  retval = phmm_fit_em(phmm, NULL, logf);

  *p = mm_get(phmm->functional_hmm->transition_matrix, 0, 1);
  *q = mm_get(phmm->functional_hmm->transition_matrix, 1, 0);

  return retval;
}
 
int main(int argc, char *argv[]) {

  /* arguments and defaults */
  int post_probs = TRUE, score = FALSE, quiet = FALSE, 
    gff = FALSE, rates_cross = FALSE, estim_lambda = TRUE, 
    estim_transitions = TRUE, two_state = TRUE, indels = FALSE,
    coding_potential = FALSE, indels_only = FALSE;
  int nrates = -1, rates_cut_idx = 1, refidx = 1, min_inform_bases = 2;
  double lambda = DEFAULT_LAMBDA, p = DEFAULT_P, q = DEFAULT_Q;
  msa_format_type msa_format = SS;
  FILE *viterbi_f = NULL, *lnl_f = NULL, *log_f = NULL;
  List *states = NULL, *pivot_states = NULL, *min_inform_str = NULL, 
    *mod_fname_list;
  char *seqname = NULL, *idpref = NULL;
  HMM *hmm = NULL;
  
  struct option long_opts[] = {
    {"states", 1, 0, 'S'},
    {"hmm", 1, 0, 'H'},
    {"viterbi", 1, 0, 'V'},
    {"no-post-probs", 0, 0, 'n'},
    {"msa-format", 1, 0, 'i'},
    {"rates-cross", 0, 0, 'X'},
    {"lambda", 1, 0, 'l'},
    {"cut-at", 1, 0, 'C'},
    {"transitions", 1, 0, 't'},
    {"nrates", 1, 0, 'k'},
    {"log", 1, 0, 'g'},
    {"refidx", 1, 0, 'r'},
    {"suppress-missing", 0, 0, 'x'}, /* for backward compatibility */
    {"reflect-strand", 1, 0, 'U'},
    {"catmap", 1, 0, 'c'},
    {"indels", 0, 0, 'I'},
    {"min-informative-types", 1, 0, 'M'},
    {"min-informative-bases", 1, 0, 'm'},
    {"lnl", 1, 0, 'L'},
    {"seqname", 1, 0, 'N'},
    {"idpref", 1, 0, 'P'},
    {"score", 0, 0, 's'},
    {"coding-potential", 0, 0, 'p'},
    {"indels-only", 0, 0, 'J'},
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
  CategoryMap *cm = NULL;
  char *mods_fname = NULL;
  indel_mode_type indel_mode;

  while ((c = getopt_long(argc, argv, "S:H:V:ni:k:l:C:t:r:xL:s:N:P:g:U:c:IJM:m:pXqh", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'S':
      states = get_arg_list(optarg);
      break;
    case 'H':
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
    case 'C':
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
      if (!strcmp(optarg, "-")) log_f = stderr;
      else log_f = fopen_fname(optarg, "w+");
      break;
    case 'r':
      refidx = get_arg_int_bounds(optarg, 0, INFTY);
      break;
    case 'x':
      /* do nothing; left in for backward compatibility */
      break;
    case 'U':
      pivot_states = get_arg_list(optarg); /* we want strings not ints
                                             for phmm_new */
      break;
    case 'I':
      indels = TRUE;
      break;
    case 'J':
      indels_only = TRUE;
      two_state = FALSE;
      indels = TRUE;
      post_probs = FALSE;
      break;
    case 'M':
      min_inform_str = get_arg_list(optarg);
      break;
    case 'm':
      min_inform_bases = get_arg_int_bounds(optarg, 0, INFTY);
      break;
    case 'c':
      cm = cm_new_string_or_file(optarg);
      break;
    case 'L':
      lnl_f = fopen_fname(optarg, "w+");
      break;
    case 'N':
      seqname = optarg;
      break;
    case 'P':
      idpref = optarg;
      break;
    case 's':
      score = TRUE;
      break;
    case 'p':
      coding_potential = TRUE;
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

  if (indels_only && (hmm != NULL || rates_cross))
    die("ERROR: --indels-only cannot be used with --hmm or --rates-cross.\n");

  if (cm != NULL && hmm == NULL) 
    die("ERROR: --catmap can only be used with --hmm.\n");
  
  if (indels == TRUE && rates_cross)
    die("ERROR: --indels cannot be used with --rates-cross.\n");

  if ((!coding_potential && optind != argc - 2) ||
      (coding_potential && optind != argc - 2 && optind != argc - 1))
    die("ERROR: extra or missing arguments.  Try '%s -h'.\n", argv[0]);

  mods_fname = (optind == argc - 2 ? argv[argc - 1] : NULL);
  /* if there are two args, mods are the second one; otherwise will
     use default mods for coding potential (see below) */
  
  /* set defaults for coding-potential mode */
  if (coding_potential) {
    char tmp[5000];
    two_state = FALSE;
    if (cm == NULL) cm = cm_new_string_or_file("NCATS=4; CNS 1; CDS 2-4");
    if (hmm == NULL) {
      sprintf(tmp, "%s/data/phastCons/simple-coding.hmm", PHAST_HOME);
      hmm = hmm_new_from_file(fopen_fname(tmp, "r"));
    }
    if (mods_fname == NULL) {
      sprintf(tmp, "\
%s/data/exoniphy/mammals/r3.ncns.mod,\
%s/data/exoniphy/mammals/r3.cns.mod,\
%s/data/exoniphy/mammals/r3.cds-1.mod,\
%s/data/exoniphy/mammals/r3.cds-2.mod,\
%s/data/exoniphy/mammals/r3.cds-3.mod", 
              PHAST_HOME, PHAST_HOME, PHAST_HOME, PHAST_HOME, PHAST_HOME);
      mods_fname = tmp;
    }
    if (states == NULL) states = get_arg_list("CDS");
    if (pivot_states == NULL) pivot_states = get_arg_list("background,CNS");
    if (min_inform_str == NULL) min_inform_str = get_arg_list("CDS");
    indels = TRUE;
  }
  
  /* read tree models first (alignment may take a while) */
  mod_fname_list = get_arg_list(mods_fname);

  if ((rates_cross || two_state || indels_only) && lst_size(mod_fname_list) != 1)
    die("ERROR: only one tree model allowed unless --hmm.\n");
    
  mod = (TreeModel**)smalloc(sizeof(TreeModel*) * lst_size(mod_fname_list));
  for (i = 0; i < lst_size(mod_fname_list); i++) {
    String *fname = lst_get_ptr(mod_fname_list, i);
    if (!quiet)
      fprintf(stderr, "Reading tree model from %s...\n", fname->chars);
    mod[i] = tm_new_from_file(fopen_fname(fname->chars, "r"));
    mod[i]->use_conditionals = 1; 
  }

  /* set min informative bases, if necessary */
  if (min_inform_str != NULL) {
    List *l = cm_get_category_list(cm, min_inform_str, 0);
    for (i = 0; i < lst_size(l); i++) 
      mod[lst_get_int(l, i)]->min_informative = min_inform_bases;
    lst_free(l);
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
  msa_remove_N_from_alph(msa);  /* for backward compatibility */
  if (msa_format == SS && msa->ss->tuple_idx == NULL) 
    die("ERROR: Ordered representation of alignment required.\n");

  for (i = 0; i < lst_size(mod_fname_list); i++)
    tm_prune(mod[i], msa, !quiet);

  /* use file name root for default seqname */
  if (viterbi_f != NULL && (seqname == NULL || idpref == NULL)) {
    String *tmp = str_new_charstr(argv[optind]);
    if (!str_equals_charstr(tmp, "-")) {
      str_remove_path(tmp);
      str_root(tmp, '.');
      if (idpref == NULL) idpref = strdup(tmp->chars);
      str_root(tmp, '.');         /* apply one more time for double suffix */
      if (seqname == NULL) seqname = tmp->chars;    
    }
  }

  /* set up states */
  if (states == NULL) {
    states = lst_new_ptr(1);
    lst_push_ptr(states, str_new_charstr("0"));
  }

  if (two_state) {
    if (!quiet) 
      fprintf(stderr, "Partitioning at rate category %d to create 'conserved' and 'nonconserved' states...\n", rates_cut_idx);
    setup_rates_cut(&hmm, &mod, &cm, nrates, rates_cut_idx, p, q);
  }
  else if (cm == NULL)
    cm = cm_create_trivial(lst_size(mod_fname_list)-1, NULL);

  /* set up PhyloHmm */
  if (!indels) indel_mode = MISSING_DATA;
  else if (hmm == NULL || hmm->nstates == cm->ncats + 1)
    indel_mode = PARAMETERIC;
  else indel_mode = NONPARAMETERIC;

  phmm = phmm_new(hmm, mod, cm, pivot_states, indel_mode);

  if (rates_cross) {
    if (!quiet) 
      fprintf(stderr, "Creating %d scaled versions of tree model...\n", nrates);
    phmm_rates_cross(phmm, nrates, lambda, TRUE);
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
    lnl = fit_rates_cut(phmm, &p, &q, log_f);
    if (!quiet) fprintf(stderr, " (p = %f. q = %f)\n", p, q);
  }

  /* estimate indel parameters only, if necessary */
  else if (indels_only) {
    if (!quiet) fprintf(stderr, "Estimating parameters for indel model...");
    lnl = phmm_fit_em(phmm, msa, log_f);
    if (!quiet) fprintf(stderr, "...\n");
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
                                            idpref, NULL, "phastCons_predicted");
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

    if (!quiet) fprintf(stderr, "Computing posterior probabilities...\n");

    postprobs = phmm_postprobs_cats(phmm, states, &lnl);

    /* print to stdout */
    for (j = 0, k = 0; j < msa->length; j++) {
      if (refidx == 0 || msa_get_char(msa, refidx-1, j) != GAP_CHAR) {
        if (!msa_missing_col(msa, refidx, j))
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
