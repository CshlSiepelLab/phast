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

#define DEFAULT_CONSERVED_SCALE 0.3

/* functions implemented below and used internally */
void setup_rates_cut(HMM **hmm, CategoryMap **cm, double p, double q);

double fit_rates_cut(PhyloHmm *phmm, MSA *msa, int estim_func, int estim_indels,
                     int estim_trees, double *p, double *q, 
                     double *alpha_0, double *beta_0, double *omega_0, 
                     double *alpha_1, double *beta_1, double *omega_1, 
                     double conserved_scale, double target_coverage, FILE *logf);

void reestimate_trees(void **models, int nmodels, void *data, 
                      double **E, int nobs, FILE *logf);

void phmm_estim_trans_em_coverage(HMM *hmm, void *data, double **A);

void usage(char *prog) {
  printf("\n\
PROGRAM: %s\n\
\n\
USAGE: %s [OPTIONS] <msa_fname> <mod_fname_list>\n\
    The alignment file, <msa_fname>, can be in any of several file\n\
    formats (see --msa-format); SS format is assumed by default (see\n\
    msa_view).  The list of *.mod files, <mod_fname_list>, should be\n\
    comma delimited.  (Files should be as output by phyloFit.)\n\
\n\
DESCRIPTION:\n\
    Identify conserved elements or produce conservation scores, given\n\
    a multiple alignment and a phylo-HMM.  By default, a phylo-HMM\n\
    consisting of two states is assumed: a \"conserved\" state and a\n\
    \"non-conserved\" state.  Separate phylogenetic models can be\n\
    specified for these two states, e.g.,\n\
\n\
        phastCons myfile.ss conserved.mod,non-conserved.mod\n\
\n\
    or a single model can be given for the non-conserved state, e.g.,\n\
\n\
        phastCons myfile.ss --conserved-scale 0.5 non-conserved.mod\n\
\n\
    in which case the model for the conserved state will be obtained\n\
    by multiplying all branch lengths by the specified factor (must be\n\
    between 0 and 1; if not specified, the default value of 0.3 will\n\
    be used).\n\
\n\
    By default, the phylogenetic models will be left unaltered, but if\n\
    the --estimate-trees option is used, e.g.,\n\
\n\
        phastCons myfile.ss init.mod --estimate-trees newtree\n\
\n\
    then the phylogenetic models for the two states will be estimated\n\
    from the data and the given model (there must be only one in this\n\
    case) will be used for initialization only.  (The estimated models\n\
    for the two states will differ only by a scale factor, which will\n\
    be estimated from the data.)\n\
\n\
    The transition probabilities for the HMM will be estimated from\n\
    the data using an EM algorithm.  This behavior can be overridden\n\
    using the --transitions option (see below).  Also, the transition\n\
    probabilities can be constrained such that the portion of sites\n\
    predicted to be conserved will be approximately equal to some\n\
    target value.  For example, to estimate parameters consistent with\n\
    identifying the 5%% most conserved sites in a genome, use\n\
    --target-coverage 0.05.  (You must also consider the coverage of\n\
    your alignments, however.  For example, if only 40%% of the genome\n\
    of interest is aligned to other genomes and you want to identify\n\
    the 5%% most conserved sites, use a target coverage of 0.05/0.4 =\n\
    0.125.)\n\
\n\
    It's also possible to use a a k-state phylo-HMM of the kind\n\
    described by Felsenstein and Churchill (1996) (see --rates-cross),\n\
    or to define the state-transition structure of the HMM explicitly\n\
    (see --hmm).  In the Felsenstein-Churchill case, a single\n\
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
    tab-delimited two-column format (position and probability).  The\n\
    set of states whose total (marginal) posterior probability is\n\
    reported can be changed using the --states option.  In addition,\n\
    discrete elements can be predicted using the --viterbi option, and\n\
    they can be assigned log-odds scores using the --score option.\n\
    The set of states considered when predicting discrete elements is\n\
    also defined by --states.\n\
\n\
    This program is written in a general enough way that it can be\n\
    useful for other things besides analyzing rates of substitution\n\
    and evolutionary conservation.  E.g., by giving it a simple\n\
    gene-finding phylo-HMM (e.g., with a state for non-coding regions\n\
    and three states for the three codon positions), and specifying\n\
    the coding states via --states, you can obtain posterior\n\
    probabilities that can be interpreted as a measure of \"coding\n\
    potential.\" (The --coding-potential option sets appropriate\n\
    defaults for such an application).\n\
\n\
EXAMPLES:\n\
    1. Fit a phylogenetic model to a data set, then produce\n\
    conservation scores via a two-state phylo-HMM.  Use an input file\n\
    in SS format.\n\
\n\
        phyloFit --tree mytree.nh mydata.ss --msa-format SS \\\n\
            --out-root init\n\
        phastCons mydata.ss init.mod --target-coverage 0.10 \\\n\
            --estimate-trees newtree > cons.dat\n\
\n\
    2. As in (1), but also predict discrete conserved elements and\n\
    score them using log-odds scoring.  This time, use the tree models\n\
    that were estimated above.  Write predictions in BED format to\n\
    elements.bed.\n\
\n\
        phastCons mydata.ss newtree.cons.mod,newtree.noncons.mod \\\n\
             --target-coverage 0.10 --viterbi elements.bed \\\n\
             --score > cons.dat\n\
\n\
    (if output file were \"elements.gff,\" then output would be in GFF\n\
    instead)\n\
\n\
    3. As above, but bypass finding the maximum likelihood estimate of\n\
    the state transition probabilities, and instead just use the value\n\
    given (much faster).\n\
\n\
        phastCons mydata.ss newtree.cons.mod,newtree.noncons.mod \\\n\
            --transitions 0.01 --target-coverage 0.10 \\\n\
            --viterbi elements.bed --score > cons.dat\n\
\n\
    (Notice that two transition parameters are required without\n\
    --target-coverage and only one with --target-coverage; see details\n\
    below)\n\
\n\
    4. Similar to above, but this time use a k-state HMM with transition\n\
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
    --transitions, -t [~]<p>,<q> | [~]<p>\n\
        (Optionally use with default two-state HMM) Fix the transition\n\
        probabilities of the two-state HMM as specified, rather than\n\
        estimating them by maximum likelihood.  Alternatively, if\n\
        first character of argument is '~', estimate parameters, but\n\
        initialize with specified values.  The argument <p> is the\n\
        probability of transitioning from the conserved to the\n\
        non-conserved state, and <q> is the probability of the reverse\n\
        transition (probabilities of self transitions are thus 1-<p>\n\
        and 1-<q> and expected lengths of conserved and nonconserved\n\
        elements are (1-p)/p and (1-q)/q, respectively).  If --target-\n\
        coverage is used (see below), then <p> may be specified alone.\n\
\n\
    --target-coverage, -C <coverage>\n\
        (Optionally use with default two-state HMM) Constrain\n\
        transition parameters such that the expected portion of\n\
        conserved sites is <coverage> (a number between 0 and 1).\n\
        This is a *prior* rather than *posterior* expectation and\n\
        assumes stationarity of the state-transition process.  Adding\n\
        this constraint causes the ratio p/q to be fixed at\n\
        (1-coverage)/coverage.  Therefore, any value of <q> specified\n\
        via --transitions will be ignored; only <p> will be used.\n\
\n\
    --expected-lengths, -E [~]<len1>,<len2> | [~]<len1>\n\
        (Alternative to --transitions)  Set transition probabilities\n\
        such that the expected lengths of conserved and non-conserved\n\
        regions are <len1> and <len2>, respectively.  Using this\n\
        option is equivalent to using --transitions with arguments of\n\
        1/(1+<len1>) and 1/(1+<len2>).  As with --transitions, the\n\
        second argument is optional (and will be ignored) if\n\
        --target-coverage is used.\n\
\n\
    --ignore-missing, -z\n\
        (For use when estimating transition probabilities) Ignore\n\
        regions of missing data (i.e., in all sequences but the\n\
        reference sequence, excluding those specified by\n\
        --not-informative) when estimating transition probabilities.\n\
        Can help avoid too-low estimates of <p> and <q> or too-high\n\
        estimates of <lambda>.  Warning: this option should be used\n\
        with --no-post-probs and without --viterbi because coordinates\n\
        in output will be unrecognizable.\n\
\n\
    --nrates, -k <nrates> | <nrates_conserved,nrates_nonconserved>\n\
        (Optionally use with a discrete-gamma model) Assume the\n\
        specified number of rate categories, instead of the number\n\
        given in the *.mod file.  The shape parameter 'alpha' will be\n\
        as given in the *.mod file.  In the case of the default\n\
        two-state HMM, two values can be specified, for the numbers of\n\
        rates for the conserved and the nonconserved states, resp.\n\
\n\
 (Tree models)\n\
    --estimate-trees, -T <fname_root>\n\
        (Optionally use with default two-state HMM) Re-estimate tree\n\
        model parameters, in the context of the two-state HMM, and\n\
        write new models to files with given file-name root (filenames\n\
        will be given suffixes \".cons.mod\" and \".noncons.mod\").  \n\
\n\
    --conserved-scale, -R <scale_factor>\n\
        (Optionally use with default two-state HMM) Set the *scale*\n\
        (overall evolutionary rate) of the model for the conserved\n\
        state to be <scale_factor> times that of the model for the\n\
        non-conserved state.  The parameter <scale_factor> must be\n\
        between 0 and 1; default is 0.3.  If used with\n\
        --estimate-trees, the specified value will be used for\n\
        initialization only (the scale factor will be estimated).\n\
        This option is ignored if two tree models are given.\n\
\n\
    --gc, -G <val>\n\
        (Optionally use with --estimate-trees) Assume equilibrium base\n\
        frequencies consistent with the given average G+C content when\n\
        estimating tree models.  (The frequencies of G and C will be\n\
        set to <val>/2 and the frequencies of A and T will be set to\n\
        (1-<val>)/2.)  This option overrides the default behavior of\n\
        estimating the equilibrium frequencies to be equal to the\n\
        relative frequencies observed in the data.  It can be useful\n\
        when model parameters are to be estimated separately, e.g.,\n\
        for fragments of genome-wide alignments, then combined into a\n\
        single set of estimates which are to be applied globally.\n\
        The argument <val> should be between 0 and 1.\n\
\n\
    --extrapolate, -e <phylog.nh> | default\n\
        Extrapolate to a larger set of species based on the given\n\
        phylogeny (Newick-format).  The trees in the given tree models\n\
        (*.mod files) must be subtrees of the larger phylogeny.  For\n\
        each tree model M, a copy will be created of the larger\n\
        phylogeny, then scaled such that the total branch length of\n\
        the subtree corresponding to M's tree equals the total branch\n\
        length of M's tree; this new version will then be used in\n\
        place of M's tree.  (Any species name present in this tree but\n\
        not in the data will be ignored.)  If the string \"default\"\n\
        is given instead of a filename, then a phylogeny for 25\n\
        vertebrate species, estimated from sequence data for Target 1\n\
        (CFTR) of the NISC Comparative Sequencing Program (Thomas et\n\
        al., 2003), will be assumed.\n\
\n\
 (Indels, forward/reverse strands, missing data, and coding potential)\n\
    --indels, -I\n\
        (Optionally use with --hmm) Expand HMM state space to model\n\
        indels as described in Siepel & Haussler (2004).\n\
\n\
    --max-micro-indel, -Y <length> \n\
        (Optionally use with --indels) Maximum length of an alignment\n\
        gap to be considered a \"micro-indel\" and therefore\n\
        addressed by the indel model.  Gaps longer than this threshold\n\
        will be treated as missing data.  Default value is 20.\n\
\n\
    --indel-params, -D [~]<alpha_0,beta_0,omega_0,alpha_1,beta_1,omega_1>\n\
        (Optionally use with --indels and default two-state HMM) Fix\n\
        the indel parameters at (alpha_0, beta_0, omega_0) for the\n\
        conserved state and at (alpha_1, beta_1, omega_1) for the\n\
        non-conserved state, rather than estimating them by maximum\n\
        likelihood.  Alternatively, if first character of argument is\n\
        '~', estimate parameters, but initialize with specified\n\
        values.  Alpha_j is the rate of insertion events per\n\
        substitution per site in state j (typically ~0.05), beta_j is\n\
        the rate of deletion events per substitution per site in state\n\
        j (typically ~0.05), and omega_j is approximately the inverse\n\
        of the expected indel length in state j (typically 0.2-0.5).\n\
\n\
    --reflect-strand, -U <pivot_states>\n\
        (Optionally use with --hmm) Given an HMM describing the\n\
        forward strand, create a larger HMM that allows for features\n\
        on both strands by \"reflecting\" the original HMM about the\n\
        specified \"pivot\" states.  The new HMM will be used for\n\
        prediction on both strands.  States can be specified by number\n\
        (indexing starts with 0), or if --catmap, by category name.\n\
\n\
    --require-informative, -M <list>\n\
        Require \"informative\" columns (i.e., columns with\n\
        more than two non-missing-data characters, excluding sequences\n\
        specified by --not-informative) in the specified HMM states.\n\
        States can be specified by number (indexing starts with 0) or,\n\
        if --catmap is used, by category name.  Non-informative\n\
        columns will be given emission probabilities of zero.  This\n\
        option can be used to eliminate false positive predictions in\n\
        regions of missing data.\n\
\n\
    --not-informative, -F <list>\n\
        Do not consider the specified sequences (listed by name) when\n\
        deciding whether a column is informative.  This option can be\n\
        useful when sequences are present that are very close to the\n\
        reference sequence and thus do not contribute much in the way\n\
        of phylogenetic information.  E.g., one might use\n\
        \"--not-informative chimp\" with a human-referenced multiple\n\
        alignment including chimp sequence.\n\
\n\
    --coding-potential, -p\n\
        Use parameter settings that cause output to be interpretable\n\
        as a coding potential score.  By default, a simplified version\n\
        of exoniphy's phylo-HMM is used, with a noncoding\n\
        (background) state, a conserved non-coding (CNS) state, and\n\
        states for the three codon positions.  This option implies\n\
        --catmap \"NCATS=4; CNS 1; CDS 2-4\" --hmm <default-HMM-file>\n\
        --states CDS --reflect-strand background,CNS\n\
        --require-informative CDS, plus a set of default *.mod files.\n\
        (All of these options can be overridden; --coding potential can\n\
        be used with or without --indels.)\n\
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
    --alias, -A <alias_def>\n\
        Alias names in input alignment according to given definition,\n\
        e.g., \"hg17=human; mm5=mouse; rn3=rat\".  Useful with default\n\
        *.mod files, e.g., with --coding-potential.  (Default models\n\
        use generic common names such as \"human\", \"mouse\", and\n\
        \"rat\".  This option allows a mapping to be established\n\
        between the leaves of trees in these files and the sequences\n\
        of an alignment that uses an alternative naming convention.)\n\
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
    J. Thomas et al.  2003.  Comparative analyses of multi-species\n\
      sequences from targeted genomic regions.  Nature 424:788-793.\n\
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

/* initialize equilibrium freqs for tree model; either make consistent
   with given G+C content or estimate from alignment */
void init_eqfreqs(TreeModel *mod, MSA *msa, double gc) {
  if (gc != -1) {               /* gc specified */
    if (strlen(mod->rate_matrix->states) != 4 || 
        mod->rate_matrix->inv_states[(int)'A'] < 0 ||
        mod->rate_matrix->inv_states[(int)'C'] < 0 ||
        mod->rate_matrix->inv_states[(int)'G'] < 0 ||
        mod->rate_matrix->inv_states[(int)'T'] < 0)
      die("ERROR: Four-character DNA alphabet required with --gc.\n");
    assert(gc > 0 && gc < 1);
    gsl_vector_set(mod->backgd_freqs, 
                   mod->rate_matrix->inv_states[(int)'G'], gc/2);
    gsl_vector_set(mod->backgd_freqs, 
                   mod->rate_matrix->inv_states[(int)'C'], gc/2);
    gsl_vector_set(mod->backgd_freqs, 
                   mod->rate_matrix->inv_states[(int)'A'], (1-gc)/2);
    gsl_vector_set(mod->backgd_freqs, 
                   mod->rate_matrix->inv_states[(int)'T'], (1-gc)/2);
  }
  else {                        /* estimate from alignment */
    if (mod->subst_mod == JC69 || mod->subst_mod == K80)
      gsl_vector_set_all(mod->backgd_freqs, 1.0/mod->backgd_freqs->size);
    else
      msa_get_base_freqs_tuples(msa, mod->backgd_freqs, mod->order+1, -1);
  }
}

int main(int argc, char *argv[]) {

  /* arguments and defaults */
  int post_probs = TRUE, score = FALSE, quiet = FALSE, 
    gff = FALSE, rates_cross = FALSE, estim_lambda = TRUE, 
    estim_transitions = TRUE, two_state = TRUE, indels = FALSE,
    coding_potential = FALSE, indels_only = FALSE, estim_indels = TRUE,
    estim_trees = FALSE, ignore_missing = FALSE;
  int nrates = -1, nrates2 = -1, refidx = 1, max_micro_indel = 20;
  double lambda = 0.9, p = 0.01, q = 0.01, alpha_0 = 0.05, beta_0 = 0.05, 
    omega_0 = 0.45, alpha_1 = 0.05, beta_1 = 0.05, omega_1 = 0.2, gc = -1,
    target_coverage = -1, conserved_scale = DEFAULT_CONSERVED_SCALE;
  msa_format_type msa_format = SS;
  FILE *viterbi_f = NULL, *lnl_f = NULL, *log_f = NULL;
  List *states = NULL, *pivot_states = NULL, *inform_reqd = NULL, 
    *mod_fname_list = NULL, *not_informative = NULL;
  char *seqname = NULL, *idpref = NULL, *estim_trees_fname_root = NULL,
    *extrapolate_tree_fname;
  HMM *hmm = NULL;
  Hashtable *alias_hash = NULL;
  TreeNode *extrapolate_tree = NULL;

  struct option long_opts[] = {
    {"states", 1, 0, 'S'},
    {"hmm", 1, 0, 'H'},
    {"viterbi", 1, 0, 'V'},
    {"no-post-probs", 0, 0, 'n'},
    {"msa-format", 1, 0, 'i'},
    {"rates-cross", 0, 0, 'X'},
    {"lambda", 1, 0, 'l'},
    {"target-coverage", 1, 0, 'C'},
    {"transitions", 1, 0, 't'},
    {"expected-lengths", 1, 0, 'E'},
    {"estimate-trees", 1, 0, 'T'},
    {"conserved-scale", 1, 0, 'R'},
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
    {"min-informative-types", 1, 0, 'M'}, /* for backward compatibilitye */
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
  char c;
  int opt_idx, i;
  List *tmpl = NULL;
  MSA *msa = NULL;
  double lnl = INFTY;
  String *tmpstr;
  TreeModel **mod;
  PhyloHmm *phmm;
  CategoryMap *cm = NULL;
  char *mods_fname = NULL, *newname;
  indel_mode_type indel_mode;

  while ((c = getopt_long(argc, argv, "S:H:V:ni:k:l:C:G:zt:E:R:T:r:xL:s:N:P:g:U:c:IY:D:JM:F:pA:Xqh", 
                          long_opts, &opt_idx)) != -1) {
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
      target_coverage = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'G':
      gc = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 't':
      if (optarg[0] != '~') estim_transitions = FALSE;
      else optarg = &optarg[1];
      tmpl = get_arg_list_dbl(optarg);
      if (lst_size(tmpl) > 2) 
        die("ERROR: bad argument to --transitions.\n");
      p = lst_get_dbl(tmpl, 0);
      if (lst_size(tmpl) == 2) q = lst_get_dbl(tmpl, 1);
      if (p <= 0 || p >= 1 || q <= 0 || q >= 1)
        die("ERROR: bad argument to --transitions.\n");
      lst_free(tmpl);
      break;
    case 'E':
      if (optarg[0] != '~') estim_transitions = FALSE;
      else optarg = &optarg[1];
      tmpl = get_arg_list_dbl(optarg);
      if (lst_size(tmpl) > 2) 
        die("ERROR: bad argument to --expected-lengths.\n");
      p = 1.0/(1.0+lst_get_dbl(tmpl, 0));
      if (lst_size(tmpl) == 2) q = 1.0/(1.0+lst_get_dbl(tmpl, 1));
      if (p <= 0 || p >= 1 || q <= 0 || q >= 1)
        die("ERROR: bad argument to --expected-lengths.\n");
      lst_free(tmpl);
      break;
    case 'T':
      estim_trees = TRUE;
      estim_trees_fname_root = optarg;
      break;
    case 'z':
      ignore_missing = TRUE;
      break;
    case 'k':
      tmpl = get_arg_list_int(optarg);
      if (lst_size(tmpl) > 2) 
        die("ERROR: too many arguments with --nrates.\n");
      nrates = lst_get_int(tmpl, 0);
      if (nrates <= 0) 
        die("ERROR: bad argument to --nrates (%d).\n", nrates);
      if (lst_size(tmpl) == 2) {
        nrates2 = lst_get_int(tmpl, 1);
        if (nrates2 <= 0) 
          die("ERROR: bad argument to --nrates (%d).\n", nrates2);
      }
      lst_free(tmpl);
      break;
    case 'R':
      conserved_scale = get_arg_dbl_bounds(optarg, 0, 1);
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
    case 'e':
      extrapolate_tree_fname = optarg;
      break;
    case 'I':
      indels = TRUE;
      break;
    case 'Y':
      max_micro_indel = get_arg_int_bounds(optarg, 1, INFTY);
      break;
    case 'D':
      if (optarg[0] != '~') estim_indels = FALSE;
      else optarg = &optarg[1];
      tmpl = get_arg_list_dbl(optarg);
      if (lst_size(tmpl) != 6) die("ERROR: bad argument to --indel-params.\n");
      alpha_0 = lst_get_dbl(tmpl, 0);
      beta_0 = lst_get_dbl(tmpl, 1);
      omega_0 = lst_get_dbl(tmpl, 2);
      alpha_1 = lst_get_dbl(tmpl, 3);
      beta_1 = lst_get_dbl(tmpl, 4);
      omega_1 = lst_get_dbl(tmpl, 5);
      if (alpha_0 < 0 || beta_0 < 0 || omega_0 < 0 || 
          alpha_1 < 0 || beta_1 < 0 || omega_1 < 0)
        die("ERROR: bad argument to --indel-params.\n");
      lst_free(tmpl);
      break;
    case 'J':
      indels_only = TRUE;
      two_state = FALSE;
      indels = TRUE;
      post_probs = FALSE;
      break;
    case 'M':
      inform_reqd = get_arg_list(optarg);
      break;
    case 'F':
      not_informative = get_arg_list(optarg);
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
    case 'A':
      alias_hash = make_name_hash(optarg);
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

  if ((estim_trees || target_coverage != -1) && !two_state)
    die("ERROR: --estimate-trees, and --target-coverage can only be used with default two-state HMM.\n");

  if (cm != NULL && hmm == NULL) 
    die("ERROR: --catmap can only be used with --hmm.\n");
  
  if (indels == TRUE && rates_cross)
    die("ERROR: --indels cannot be used with --rates-cross.\n");

  if (nrates != -1 && hmm != NULL)
    die("ERROR: --nrates currently can't be used with --hmm.\n");
  
  if ((!coding_potential && optind != argc - 2) ||
      (coding_potential && optind != argc - 2 && optind != argc - 1))
    die("ERROR: extra or missing arguments.  Try '%s -h'.\n", argv[0]);

  if (!indels) estim_indels = FALSE;

  if (!strcmp(extrapolate_tree_fname, "default")) {
    extrapolate_tree_fname = smalloc(1000 * sizeof(char));
    sprintf(extrapolate_tree_fname, 
            "%s/data/exoniphy/mammals/cftr25_hybrid.nh", PHAST_HOME);
  }
  if (extrapolate_tree_fname != NULL)
    extrapolate_tree = tr_new_from_file(fopen_fname(extrapolate_tree_fname, "r"));

  mods_fname = (optind == argc - 2 ? argv[argc - 1] : NULL);
  /* if there are two args, mods are the second one; otherwise will
     use default mods for coding potential (see below) */
  
  /* set defaults for coding-potential mode */
  if (coding_potential) {
    char tmp[5000];
    two_state = FALSE;
    if (cm == NULL) cm = cm_new_string_or_file("NCATS=4; CNS 1; CDS 2-4");
    if (hmm == NULL) {
      sprintf(tmp, "%s/data/phastCons/%s", PHAST_HOME,
              indels ? "simple-coding-indels.hmm" : "simple-coding.hmm");
      if (!quiet) fprintf(stderr, "Reading HMM from %s...\n", tmp);
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
    if (inform_reqd == NULL) inform_reqd = get_arg_list("CDS");
  }
  
  /* read alignment */
  if (!quiet)
    fprintf(stderr, "Reading alignment from %s...\n", argv[optind]);
  msa = msa_new_from_file(fopen_fname(argv[optind], "r"), msa_format, NULL);
  msa_remove_N_from_alph(msa);  /* for backward compatibility */
  if (msa_format == SS && msa->ss->tuple_idx == NULL) 
    die("ERROR: Ordered representation of alignment required.\n");
  if (msa->ss == NULL) ss_from_msas(msa, 1, TRUE, NULL, NULL, NULL, -1);
                                /* SS assumed below */

  /* rename if aliases are defined */
  if (alias_hash != NULL) {
    for (i = 0; i < msa->nseqs; i++) {
      if ((newname = hsh_get(alias_hash, msa->names[i])) != (char*)-1) {
        free(msa->names[i]);
        msa->names[i] = strdup(newname);
      }
    }
  }

  /* mask out macro-indels, if necessary */
  if (indels) {
    /* this little hack allows gaps in refseq to be restored before
       output (needed for proper coord conversion) */
    if (msa->seqs == NULL) { ss_to_msa(msa); ss_free(msa->ss); msa->ss = NULL; }
    assert(strlen(msa->missing) >= 2);
    for (i = 0; i < msa->length; i++) 
      if (msa->is_missing[(int)msa->seqs[0][i]]) msa->seqs[0][i] = msa->missing[1];
                                /* msa->missing[0] is used in msa_mask_macro_indels */

    msa_mask_macro_indels(msa, max_micro_indel);
  }

  /* Set up array indicating which seqs are informative, if necessary */
  if (not_informative != NULL)
    msa_set_informative(msa, not_informative);

  /* strip missing columns, if necessary */
  if (ignore_missing)
    ss_strip_missing(msa, refidx);

  /* read tree models */
  mod_fname_list = get_arg_list(mods_fname);

  if ((rates_cross || indels_only) && lst_size(mod_fname_list) != 1)
    die("ERROR: only one tree model allowed with --rates-cross and --indels-only.\n");

  if (two_state && lst_size(mod_fname_list) > 2)
    die("ERROR: must specify either one or two tree models with default two-state model.\n");
    
  mod = (TreeModel**)smalloc(sizeof(TreeModel*) * lst_size(mod_fname_list));
  for (i = 0; i < lst_size(mod_fname_list); i++) {
    String *fname = lst_get_ptr(mod_fname_list, i);
    if (!quiet)
      fprintf(stderr, "Reading tree model from %s...\n", fname->chars);
    mod[i] = tm_new_from_file(fopen_fname(fname->chars, "r"));
    mod[i]->use_conditionals = 1;     

    if (extrapolate_tree != NULL) {
      TreeNode *t = tr_create_copy(extrapolate_tree);
      double scale = tr_scale_by_subtree(t, mod[i]->tree);
      if (!quiet) 
        fprintf(stderr, "Extrapolating based on %s (scale=%f)...\n", 
                extrapolate_tree_fname, scale);
      tr_free(mod[i]->tree);
      mod[i]->tree = t;
    }

    tm_prune(mod[i], msa, !quiet); /* prune away extra species, if possible */
  }

  /* initial checks and setup of tree models with two-state and
     rates-cross options */
  if (two_state) {
    if (lst_size(mod_fname_list) == 2 && estim_trees)
      die("ERROR: If re-estimating tree models, pass in only one model for initialization.\n");
    if (mod[0]->empirical_rates || 
        (lst_size(mod_fname_list) == 2 && mod[1]->empirical_rates))
      die("ERROR: nonparameteric rate variation not allowed with default two-state HMM.\n");

    /* set equilibrium frequencies if estimating tree models */
    if (estim_trees) init_eqfreqs(mod[0], msa, gc);

    if (lst_size(mod_fname_list) == 1) { /* create 2nd tree model &
                                            rescale first */
      mod = srealloc(mod, 2 * sizeof(void*));
      mod[1] = tm_create_copy(mod[0]);
      tm_scale(mod[0], conserved_scale, TRUE);
    }
    if (nrates != -1 && nrates != mod[0]->nratecats) 
      tm_reinit(mod[0], mod[0]->subst_mod, nrates, mod[0]->alpha, NULL, NULL);
    if (nrates2 != -1 && nrates2 != mod[1]->nratecats) 
      tm_reinit(mod[1], mod[1]->subst_mod, nrates2, mod[1]->alpha, NULL, NULL);
  }
  else if (rates_cross) {
    if (mod[0]->nratecats <= 1)
      die("ERROR: a tree model allowing for rate variation is required.\n");
    if (nrates != -1 && mod[0]->empirical_rates)
      die("ERROR: can't use --nrates with nonparameteric rate model.\n");
    if (nrates == -1) nrates = mod[0]->nratecats;
  }

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
    else if (seqname == NULL) seqname = "refseq";
  }

  /* set up states */
  if (states == NULL) {
    states = lst_new_ptr(1);
    lst_push_ptr(states, str_new_charstr("0"));
  }

  if (two_state) {
    if (!quiet) 
      fprintf(stderr, "Creating 'conserved' and 'nonconserved' states in HMM...\n");
    if (target_coverage != -1) {
      q = target_coverage/(1-target_coverage) * p;
      if (q >= 1) 
        die("ERROR: p=%f and coverage=%f imply q >= 1.\n", p, target_coverage);
    }
    setup_rates_cut(&hmm, &cm, p, q);
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

  /* set min informative bases, if necessary.  This has to be done
     *after* the set of models is expanded (two-state or
     rates-cross) */
  if (inform_reqd != NULL) {
    List *l = cm_get_category_list(cm, inform_reqd, 0);
    for (i = 0; i < lst_size(l); i++) {
      int modno = lst_get_int(l, i);
      if (modno < 0 || modno >= phmm->nmods) 
        die("ERROR: illegal argument to --require-informative.\n");
      phmm->mods[modno]->inform_reqd = TRUE;
    }
    lst_free(l);
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

  /* estimate p and q and indel params, if necessary */
  else if (two_state && (estim_transitions || estim_indels || estim_trees)) {
    if (!quiet) {
      fprintf(stderr, "Finding MLE for (");
      if (estim_transitions) 
        fprintf(stderr, "p, q%s", estim_indels || estim_trees ? ", " : "");
      if (estim_indels) 
        fprintf(stderr, "alpha_0, beta_0, omega_0, alpha_1, beta_1, omega_1%s",
                estim_trees ? ", " : "");
      if (estim_trees)
        fprintf(stderr, "[tree models]");
      fprintf(stderr, ")...\n");
    }
    lnl = fit_rates_cut(phmm, msa, estim_transitions, estim_indels, estim_trees,
                        &p, &q, &alpha_0, &beta_0, &omega_0, 
                        &alpha_1, &beta_1, &omega_1, conserved_scale,
                        target_coverage, log_f);
    if (!quiet && (estim_transitions || estim_indels)) {      
      fprintf(stderr, "(");
      if (estim_transitions)
        fprintf(stderr, "p = %f. q = %f%s", p, q, estim_indels ? ", " : "");
      if (estim_indels)
        fprintf(stderr, "alpha_0 = %f, beta_0 = %f, omega_0 = %f, alpha_1 = %f, beta_1 = %f, omega_1 = %f", alpha_0, beta_0, omega_0, alpha_1, beta_1, omega_1);
      fprintf(stderr, ")\n");
    }

    if (estim_trees) {
      char cons_fname[STR_MED_LEN], noncons_fname[STR_MED_LEN];
      sprintf(cons_fname, "%s.cons.mod", estim_trees_fname_root);
      sprintf(noncons_fname, "%s.noncons.mod", estim_trees_fname_root);
      if (!quiet)
        fprintf(stderr, "Writing re-estimated tree models to %s and %s...\n", 
                cons_fname, noncons_fname);
      tm_print(fopen_fname(cons_fname, "w+"), phmm->mods[0]);
      tm_print(fopen_fname(noncons_fname, "w+"), phmm->mods[1]);
    }
  }

  /* estimate indel parameters only, if necessary */
  else if (indels_only) {
    if (!quiet) fprintf(stderr, "Estimating parameters for indel model...");
    lnl = phmm_fit_em(phmm, msa, TRUE, FALSE, log_f);
    if (!quiet) fprintf(stderr, "...\n");
  }

  /* still have to set indel params if not estimating */
  else if (indel_mode == PARAMETERIC) {
    phmm->alpha[0] = alpha_0; phmm->beta[0] = beta_0; phmm->omega[0] = omega_0;
    phmm->alpha[1] = alpha_1; phmm->beta[1] = beta_1; phmm->omega[1] = omega_1;
    phmm_reset(phmm);
  }

  /* before output, have to restore gaps in reference sequence, for
     proper coord conversion */
  if (indels && (post_probs || viterbi_f != NULL)) {
    ss_free(msa->ss); msa->ss = NULL; /* msa->seqs must already exist */
    for (i = 0; i < msa->length; i++) 
      if (msa->seqs[0][i] == msa->missing[0]) msa->seqs[0][i] = GAP_CHAR;
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
    else if (two_state && (estim_transitions || estim_indels)) {
      fprintf(lnl_f, "(");
      if (estim_transitions)
        fprintf(lnl_f, "p = %f, q = %f%s", p, q, estim_indels ? ", " : "");
      if (estim_indels)
        fprintf(lnl_f, "alpha_0 = %f, beta_0 = %f, omega_0 = %f, alpha_1 = %f, beta_1 = %f, omega_1 = %f", alpha_0, beta_0, omega_0, alpha_1, beta_1, omega_1);
      fprintf(lnl_f, ")\n");
    }
  }

  if (!quiet)
    fprintf(stderr, "Done.\n");

  return 0;
}

/* Set up HMM and category map for two-state case */
void setup_rates_cut(HMM **hmm, CategoryMap **cm, double p, double q) {

  *hmm = hmm_new_nstates(2, TRUE, FALSE);

  /* set HMM transitions according to p and q */
  mm_set((*hmm)->transition_matrix, 0, 0, 1-p);
  mm_set((*hmm)->transition_matrix, 0, 1, p);
  mm_set((*hmm)->transition_matrix, 1, 0, q);
  mm_set((*hmm)->transition_matrix, 1, 1, 1-q);

  /* just use stationary distribution for begin transitions */
  gsl_vector_set((*hmm)->begin_transitions, 0, q/(p+q));
  gsl_vector_set((*hmm)->begin_transitions, 1, p/(p+q));

  hmm_reset(*hmm);

  /* define two-category category map */
  *cm = cm_create_trivial(1, "cons_");
}

/* Estimate parameters for the two-state model using an EM algorithm.
   Any or all of the parameters 'p' and 'q', the indel parameters, and
   the tree models themselves may be estimated.  Returns ln
   likelihood. */
double fit_rates_cut(PhyloHmm *phmm, MSA *msa, int estim_func, int estim_indels,
                     int estim_trees, double *p, double *q, 
                     double *alpha_0, double *beta_0, double *omega_0, 
                     double *alpha_1, double *beta_1, double *omega_1, 
                     double conserved_scale, double target_coverage, FILE *logf) {
  double retval;

  mm_set(phmm->functional_hmm->transition_matrix, 0, 0, 1-*p);
  mm_set(phmm->functional_hmm->transition_matrix, 0, 1, *p);
  mm_set(phmm->functional_hmm->transition_matrix, 1, 0, *q);
  mm_set(phmm->functional_hmm->transition_matrix, 1, 1, 1-*q);
                                /* note that phmm->functional_hmm ==
                                   phmm->hmm if no indel model */

  phmm->em_data = smalloc(sizeof(EmData));
  phmm->em_data->msa = msa;
  phmm->em_data->fix_functional = !estim_func;
  phmm->em_data->fix_indel = !estim_indels;
  phmm->em_data->conserved_scale = conserved_scale;
  phmm->em_data->target_coverage = target_coverage;
  phmm->em_data->H = NULL;      /* will be defined as needed */

  if (phmm->indel_mode == PARAMETERIC) {
    phmm->alpha[0] = *alpha_0;
    phmm->beta[0] = *beta_0;
    phmm->omega[0] = *omega_0;
    phmm->alpha[1] = *alpha_1;
    phmm->beta[1] = *beta_1;
    phmm->omega[1] = *omega_1;
  }

  phmm_reset(phmm); 

  if (estim_trees) {
    msa->ncats = phmm->nmods - 1;   /* ?? */
    if (msa->ss == NULL) 
      ss_from_msas(msa, phmm->mods[0]->order+1, TRUE, NULL, NULL, NULL, -1);
    else if (msa->ss->cat_counts == NULL)
      ss_realloc(msa, msa->ss->tuple_size, msa->ss->ntuples, TRUE, TRUE);

    /* force re-initialization of tree models with rate variation;
       this is a hack that helps keep the parameterization simple */
    if (phmm->mods[0]->nratecats == 1) {
      tm_reinit(phmm->mods[0], phmm->mods[0]->subst_mod, 2, 
                phmm->mods[0]->alpha, NULL, NULL);
      phmm->mods[0]->nratecats = 1; phmm->mods[0]->alpha = -1; /* ignore rate variation */
      phmm->mods[0]->freqK[0] = phmm->mods[0]->rK[0] = 1;
    }
    if (phmm->mods[1]->nratecats == 1) {
      tm_reinit(phmm->mods[1], phmm->mods[1]->subst_mod, 2, 
                phmm->mods[1]->alpha, NULL, NULL);
      phmm->mods[1]->nratecats = 1; phmm->mods[1]->alpha = -1;
      phmm->mods[1]->freqK[0] = phmm->mods[1]->rK[0] = 1;
    }

    retval = hmm_train_by_em(phmm->hmm, phmm->mods, phmm, 1, &phmm->alloc_len, NULL, 
                             phmm_compute_emissions_em, reestimate_trees,
                             target_coverage > 0 ? 
                             phmm_estim_trans_em_coverage : phmm_estim_trans_em, 
                             phmm_get_obs_idx_em, 
                             phmm_log_em, logf);
    retval *= log(2);

    /* have to do final rescaling of tree models to get units of subst/site */
    if (phmm->mods[0]->subst_mod != JC69 && phmm->mods[0]->subst_mod != F81) {   
                                /* JC69 and F81 are exceptions */
      double scale = tm_scale_rate_matrix(phmm->mods[0]); 
      tm_scale_rate_matrix(phmm->mods[1]); /* will be the same */
      tm_scale(phmm->mods[0], scale, 0); 
      tm_scale(phmm->mods[1], scale, 0);
    }

    phmm->mods[0]->lnL = phmm->mods[1]->lnL = retval;
  }

  else {                        /* not estimating tree models */
    retval = hmm_train_by_em(phmm->hmm, phmm->mods, phmm, 1, 
                             &phmm->alloc_len, NULL, 
                             phmm_compute_emissions_copy_em, NULL,
                             target_coverage > 0 ? 
                             phmm_estim_trans_em_coverage : phmm_estim_trans_em, 
                             NULL, phmm_log_em, logf);
    retval *= log(2);
  }

  *p = mm_get(phmm->functional_hmm->transition_matrix, 0, 1);
  *q = mm_get(phmm->functional_hmm->transition_matrix, 1, 0);

  if (phmm->indel_mode == PARAMETERIC) {
    *alpha_0 = phmm->alpha[0];
    *beta_0 = phmm->beta[0];
    *omega_0 = phmm->omega[0];
    *alpha_1 = phmm->alpha[1];
    *beta_1 = phmm->beta[1];
    *omega_1 = phmm->omega[1];
  }

  return retval;
}

/* Special-purpose unpack function, adapted from tm_unpack_params */
void unpack_params_mod(TreeModel *mod, gsl_vector *params) {
  TreeNode *n;
  int assigned = 0, nodeidx, i;
  List *traversal;

  assert(!mod->estimate_backgd && !mod->empirical_rates);

  /* check parameter values */
  for (i = 0; i < params->size; i++) {
    double p = gsl_vector_get(params, i);
    if (p < 0 && abs(p) < TM_IMAG_EPS) /* consider close enough to 0 */
      gsl_vector_set(params, i, p=0);
    if (p < 0) die("ERROR: parameter %d has become negative (%f).\n", i, p);
    if (!finite(p)) die("ERROR: parameter %d is no longer finite (%f).\n", i, p);
  }
  i = 0;

  if (mod->estimate_branchlens == TM_SCALE_ONLY) 
    mod->scale = gsl_vector_get(params, i++);
  else if (mod->estimate_branchlens == TM_BRANCHLENS_ALL) {
    traversal = tr_preorder(mod->tree);
    for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
      n = lst_get_ptr(traversal, nodeidx);

      if (n->parent != NULL) {
        if ((n == mod->tree->lchild || n == mod->tree->rchild) && 
            tm_is_reversible(mod->subst_mod)) {
          n->dparent = gsl_vector_get(params, 0)/2;
          if (!assigned) {
            i++;     /* only increment the first time */
            assigned = 1;
          }
        }
        else if (n->id == mod->root_leaf_id) 
          n->dparent = 0;
        else 
          n->dparent = gsl_vector_get(params, i++);
      }
    }
  }

  /* next parameters are for rate variation */
  if (mod->nratecats > 1) 
    mod->alpha = gsl_vector_get(params, i++);
  else i++;                     /* there's always a placeholder
                                   here */

  tm_set_rate_matrix(mod, params, i);

  /* diagonalize, if necessary */
  if (mod->subst_mod != JC69 && mod->subst_mod != F81)
    mm_diagonalize(mod->rate_matrix);
}

/* Unpack all params for two-state HMM.  Used during M step of EM */
void unpack_params_phmm(PhyloHmm *phmm, gsl_vector *params) {
  unpack_params_mod(phmm->mods[0], params);
  unpack_params_mod(phmm->mods[1], params);
  phmm->em_data->conserved_scale = gsl_vector_get(params, params->size - 1);
  tm_scale(phmm->mods[0], gsl_vector_get(params, params->size-1), FALSE);
  
  if (phmm->mods[0]->nratecats > 1) 
    DiscreteGamma(phmm->mods[0]->freqK, phmm->mods[0]->rK, phmm->mods[0]->alpha, 
                  phmm->mods[0]->alpha, phmm->mods[0]->nratecats, 0); 
                                /* mods[0]->alpha will already be set */
  if (phmm->mods[1]->nratecats > 1) {
    phmm->mods[1]->alpha = gsl_vector_get(params, params->size - 2);
    DiscreteGamma(phmm->mods[1]->freqK, phmm->mods[1]->rK, phmm->mods[1]->alpha, 
                  phmm->mods[1]->alpha, phmm->mods[1]->nratecats, 0); 
  }
  tm_set_subst_matrices(phmm->mods[0]);
  tm_set_subst_matrices(phmm->mods[1]);
}
 
/* Wrapper for computation of likelihood, for use by reestimate_trees (below) */
double likelihood_wrapper(gsl_vector *params, void *data) {
  PhyloHmm *phmm = (PhyloHmm*)data;
  double retval0, retval1;

  unpack_params_phmm(phmm, params);

  retval0 = -tl_compute_log_likelihood(phmm->mods[0], phmm->em_data->msa, NULL, 0, NULL);
  retval1 = -tl_compute_log_likelihood(phmm->mods[1], phmm->em_data->msa, NULL, 1, NULL);
  
  return retval0 + retval1;
                                /* FIXME: what happens when not one to
                                   one cats and mods? */
}

/* Re-estimate phylogenetic model based on expected counts (M step of EM) */
void reestimate_trees(void **models, int nmodels, void *data, 
                      double **E, int nobs, FILE *logf) {

  PhyloHmm *phmm = (PhyloHmm*)data;
  int k, obsidx;
  gsl_vector *params, *lower_bounds, *upper_bounds;
  double ll;

  /* FIXME: what about when multiple states per model?  Need to
     collapse sufficient stats.  Could probably be done generally...
     need to use state_to_cat, etc. in deciding which categories to
     use */

  for (k = 0; k < phmm->nmods; k++) 
    for (obsidx = 0; obsidx < nobs; obsidx++) 
      phmm->em_data->msa->ss->cat_counts[k][obsidx] = E[k][obsidx];

  params = gsl_vector_alloc(tm_get_nparams(phmm->mods[1]) + 2);
  tm_params_init_from_model(phmm->mods[1], params, 0); /* unscaled branch lens */

  /* special initialization of rate-variation parameters and conserved_scale */
  gsl_vector_set(params, tm_get_nbranchlenparams(phmm->mods[0]), 
                 phmm->mods[0]->nratecats > 1 ? phmm->mods[0]->alpha : 0);
  gsl_vector_set(params, params->size - 2,  
                 phmm->mods[1]->nratecats > 1 ? phmm->mods[1]->alpha : 0);
  gsl_vector_set(params, params->size - 1, phmm->em_data->conserved_scale);

  lower_bounds = gsl_vector_calloc(params->size);
  upper_bounds = gsl_vector_alloc(params->size);
  gsl_vector_set_all(upper_bounds, INFTY);
  gsl_vector_set(upper_bounds, params->size - 1, 1); /* 0 < conserved_scale < 1 */

  if (logf != NULL)
    fprintf(logf, "\nRE-ESTIMATION OF TREE MODEL:\n");

  /* keep Hessian arround so it can be used from one iteration to the
     next */
  if (phmm->em_data->H == NULL) {
    phmm->em_data->H = gsl_matrix_alloc(params->size, params->size);
    gsl_matrix_set_identity(phmm->em_data->H);
  }

  if (opt_bfgs(likelihood_wrapper, params, phmm, &ll, lower_bounds, 
               NULL, logf, NULL, OPT_MED_PREC, phmm->em_data->H) != 0)
    die("ERROR returned by opt_bfgs.\n");

  unpack_params_phmm(phmm, params);

  if (phmm->indel_mode == PARAMETERIC)
    phmm_set_branch_len_factors(phmm);

  gsl_vector_free(params); 
  gsl_vector_free(lower_bounds);
  gsl_vector_free(upper_bounds);
}

/* Maximize HMM transition parameters subject to constrain implied by
   target coverage (M step of EM).  For use with two-state HMM.  This
   function is passed to hmm_train_by_em in phmm_fit_em */
void phmm_estim_trans_em_coverage(HMM *hmm, void *data, double **A) {

  PhyloHmm *phmm = data;
  IndelEstimData *ied = NULL;
  double **C;

  if (phmm->em_data->fix_functional && phmm->em_data->fix_indel) return;

  if (phmm->indel_mode == PARAMETERIC) {
    ied = phmm_new_ied(phmm, A);
    C = ied->fcounts;
  }
  else C = A;

  /* estimate transition probs for functional cats subject to
     constraint on coverage */
  if (!phmm->em_data->fix_functional) {
    double a, b, c, p, q, q1, q2, z, tmp;
    /* if you take the first derivative wrt q of the expression inside
       the argmax and set it to zero, you get a quadratic eqn which
       can be solved using the quadratic formula */
    z = (1-phmm->em_data->target_coverage)/phmm->em_data->target_coverage;
    a = z * (C[0][0] + C[0][1] + C[1][0] + C[1][1]);
    b = -C[0][1] - C[1][0] - C[1][1] - z * (C[0][0] + C[0][1] + C[1][0]);
    c = C[0][1] + C[1][0];

    tmp = b*b - 4*a*c;
    assert (tmp >= 0);
    tmp = sqrt(tmp);
    q1 = (-b + tmp) / (2*a);
    q2 = (-b - tmp) / (2*a);
    /* only one root can be valid */
    if (q1 < 1e-10 || z * q1 > 1 - 1e-10)                                 
      q = q2;                   /* (allow for rounding errors) */
    else q = q1;

    /* double check that derivative is really zero */
    assert(fabs(-z*C[0][0]/(1-z*q) + (C[0][1] + C[1][0])/q - C[1][1]/(1-q)) < 1e-6);

    p = z * q;
    assert(q >= 0 && q <= 1 && p >= 0 && p <= 1);

    mm_set(phmm->functional_hmm->transition_matrix, 0, 0, 1-p);
    mm_set(phmm->functional_hmm->transition_matrix, 0, 1, p);
    mm_set(phmm->functional_hmm->transition_matrix, 1, 0, q);
    mm_set(phmm->functional_hmm->transition_matrix, 1, 1, 1-q);

    /* use stationary distribution for begin transitions */
    gsl_vector_set(phmm->functional_hmm->begin_transitions, 0, q/(p+q));
    gsl_vector_set(phmm->functional_hmm->begin_transitions, 1, p/(p+q));
  }

  if (phmm->indel_mode == PARAMETERIC) {
    if (!phmm->em_data->fix_indel) phmm_em_estim_indels(phmm, ied);
    phmm_free_ied(ied);
  }

  phmm_reset(phmm);
}
