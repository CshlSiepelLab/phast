/* phyloFit - fit phylogenetic model(s) to a multiple alignment
   
   $Id: phyloFit.c,v 1.2 2004-06-04 21:56:33 acs Exp $
   Written by Adam Siepel, 2002-2004
   Copyright 2002-2004, Adam Siepel, University of California 

   Things to add:
        - if -g but not -c, then use default category map of one category for 
        each type of feature 
*/

#include <stdlib.h>
#include <stdio.h>
#include "lists.h"
#include "stringsplus.h"
#include "msa.h"
#include "gff.h"
#include "category_map.h"
#include <getopt.h>
#include "tree_model.h"
#include "fit_em.h"
#include "subst_mods.h"
#include <string.h>
#include "local_alignment.h"
#include <assert.h>
#include <ctype.h>
#include "tree_likelihoods.h"
#include "numerical_opt.h"
#include "sufficient_stats.h"

#define DEFAULT_NSITES_THRESHOLD 50

/* default starting alpha for dgamma */
#define DEFAULT_ALPHA 1

void print_usage() {
  printf("\n\
PROGRAM: phyloFit\n\
\n\
DESCRIPTION: \n\
    Fits one or more tree models to a multiple alignment by maximum \n\
    likelihood, using the specified tree topology and substitution\n\
    model.  If functional categories are defined via --features and\n\
    --catmap (see below), then a separate model will be fitted for each\n\
    category.  A description of each fitted model will be written to a\n\
    separate file, with the suffix \".mod\".  Each such file will\n\
    include a substitution rate matrix, a tree with branch lengths,\n\
    and estimates of nucleotide equilibrium frequencies.\n\
\n\
USAGE: phyloFit -m <msa_fname> -t <tree_fname> [OPTIONS]\n\
\n\
NOTE:\n\
    Use --EM and --precision MED (at most) for context-dependent models\n\
    (see below).\n\
\n\
OPTIONS:\n\
    --msa, -m <msa_fname>\n\
        (required) Name of file containing multiple sequence alignment,\n\
        in either PHYLIP format or a specified alternative (see\n\
        --msa-format).\n\
\n\
    --tree, -t <tree_fname>|<tree_string>\n\
        (required) Name of file *or* literal string defining tree\n\
        topology.  In either case, tree must be in New Hampshire (NH)\n\
        format, with the label at each leaf equal to the index of the\n\
        corresponding sequence in the alignment (indexing begins with\n\
        1).  Example: --tree \"(1,(2,3))\".  Currently, the topology\n\
        must be rooted.  When a reversible substitution model is used,\n\
        the root is ignored during the optimization procedure.\n\
\n\
    --subst-mod, -s JC69|F81|HKY85|REV|UNREST|R2|R2S|U2|U2S|R3|R3S|U3|U3S\n\
        (default REV).  Nucleotide substitution model.  JC69, F81, HKY85\n\
        REV, and UNREST have the usual meanings (see, e.g., Yang, \n\
        Goldman, and Friday, 1994).  The others (all considered \"context-\n\
        dependent\") are as defined in Siepel and Haussler, 2004.\n\
\n\
    --msa-format, -i PHYLIP|FASTA|PSU|SS|LAV\n\
        (default PHYLIP) Alignment format.  PSU is the \"raw\" text\n\
        format used by several tools developed by Webb Miller and\n\
        colleagues at Penn. State University.  SS is a\n\
        self-explanatory representation of a multiple alignment in\n\
        terms of its sufficient statistics.  LAV is the format used by\n\
        BLASTZ to represent local pairwise alignments.  If it is\n\
        selected, the alignment will be treated like a global\n\
        alignment, with unaligned portions of the target sequence\n\
        replaced by gaps.  The original sequences must be accessible\n\
        via the filenames specified in the LAV file.\n\
\n\
    --nrates, -k <nratecats>\n\
        (default 1).  Number of rate categories to use.  Specifying a\n\
        value of greater than one causes the discrete gamma model for\n\
        rate variation to be used (Yang, 1994).\n\
\n\
    --alpha, -a <alpha>\n\
        (for use with --nrates).  Initial value for alpha, the shape\n\
        parameter of the gamma distribution.  Default is 1.\n\
\n\
    --features, -g <gff_fname>\n\
        (must use with --catmap) File in GFF describing features on one\n\
        or more of the sequences in the alignment.\n\
\n\
    --catmap, -c <cat_map_fname>\n\
        (must use with --features) File defining mapping of sequence\n\
        features to category numbers.  If this option and the -g option\n\
        are selected, sites of the alignment will be labeled as specified,\n\
        and a separate model will be fitted to the sites of each\n\
        category.\n\
\n\
    --log, -l <log_fname>\n\
        Write log file to <log_fname>, describing details of the\n\
        optimization procedure.\n\
\n\
    --out-root, -o <output_fname_root>\n\
        (default is \"ftm\").  Use specified string as root filename\n\
        for all files created.\n\
\n\
    --output-tree, -T\n\
        Output a tree in NH format for each model, in addition to the\n\
        one that appears in the model file (each tree will appear in a\n\
        separate file).\n\
\n\
    --EM, -E \n\
        Fit model(s) using EM rather than with the BFGS quasi-Newton\n\
        algorithm (the default).\n\
\n\
    --precision, -p HIGH|MED|LOW\n\
        (default HIGH) Level of precision to use in estimating model\n\
        parameters.  Affects convergence criteria for iterative\n\
        algorithms: higher precision means more iterations and longer\n\
        execution time.\n\
\n\
    --do-cats, -C <cat_list>\n\
        (optionally use with --features and --catmap) Fit models for\n\
        only the specified categories (comma-delimited list).  Default\n\
        is to fit a model for every category.\n\
\n\
    --cats-cycle, -Y <cycle_size>\n\
        (alternative to --features and --catmap) Assign site categories in\n\
        cycles of the specified size, e.g., as 1,2,3,...,1,2,3 (for\n\
        cycle_size == 3).  Useful for coding sequence or, with higher\n\
        order models, to partition a data set into nonoverlapping\n\
        tuples of columns (can be used with --do-cats).\n\
\n\
    --markov-dependence, -N\n\
        (for use with context-dependent substitutions models and not\n\
        available with --EM.)  Assume Markov dependence of alignment\n\
        columns, and compute the conditional probability of each\n\
        column given its N-1 predecessors using the two-pass algorithm\n\
        described by Siepel and Haussler (2004).  (Here, N is the\n\
        \"order\" of the model, as defined by --subst-mod; e.g., N=1\n\
        for REV, N=2 for U2S, N=3 for U3S.) The alternative (the\n\
        default) is simply to work with joint probabilities of tuples\n\
        of columns.  (You can ensure that these tuples are\n\
        nonoverlapping by using --cats-cycle and --do-cats.)  The use\n\
        of joint probabilities during parameter estimation allows the\n\
        use of the --EM option and can be much faster; in addition, it\n\
        appears to produce nearly equivalent estimates.  If desired,\n\
        parameters can be estimated without --markov-dependence, and\n\
        then the likelihood can be evaluated using --lnl and\n\
        --markov-dependence together.  This gives a lower bound on the\n\
        likelihood of the Markov-dependent model.\n\
\n\
    --gap-strip, -G ALL|ANY|<seqno>\n\
        Strip columns in alignment containing all gaps, any gaps, or \n\
        a gap in the specified sequence (<seqno>; indexing starts with\n\
        one).  Default is not to strip any columns.\n\
\n\
    --weight-matrix-cats, -W <wmat_cats>\n\
        (optionally use with --features and --catmap) Argument should be a\n\
        comma-delimited list of strings, or the single string \"all\".\n\
        Fit sites corresponding to designated category names (as they\n\
        appear in the category map and GFF) as \"weight matrices\" rather\n\
        than as full tree models (that is, consider equilibrium\n\
        frequencies only).  Useful for sites corresponding to\n\
        \"signals\", such as start codons or splice sites, when training\n\
        a phylo-HMM for gene prediction.  If strings are numbers, they\n\
        will be assumed to represent category numbers rather than names.\n\
\n\
    --no-reverse, -R\n\
        (for use with --features and --catmap) Allow segments of alignment\n\
        corresponding to features on reverse strand to remain as they\n\
        are (by default, they are reverse complemented).  Strandedness\n\
        is obtained from the GFF file.\n\
\n\
    --init-model, -M <mod_fname>\n\
        Initialize with the specified tree model.  By choosing good\n\
        starting values for parameters, it is possible to reduce\n\
        execution time dramatically.  If this option is chosen, -t is\n\
        not required.  Note: currently only one mod_fname may be\n\
        specified; it will be used for all categories.\n\
\n\
    --init-random, -r\n\
        Initialize parameters randomly.  Can be used multiple times to test\n\
        whether the m.l.e. is real.\n\
\n\
    --lnl, -L\n\
        (for use with --init-model) Simply evaluate the log likelihood of\n\
        the specified tree model, without performing any further\n\
        optimization.  Can be used with --post-probs, --expected-subs, and\n\
        --expected-total-subs.\n\
\n\
    --scale-only, -B\n\
        (for use with --init-mod) Estimate only the scale of the tree,\n\
        rather than individual branch lengths (branch proportions fixed).\n\
\n\
    --estimate-freqs, -F\n\
        Estimate equilibrium frequencies by maximum likelihood, rather\n\
        than approximating them by the relative frequencies in the data.\n\
\n\
    --min-informative, -I <ninf_sites>\n\
        Require at least <ninf_sites> \"informative\" sites -- i.e., \n\
        sites at which at least two non-gap characters are present.  Default\n\
        is %d.\n\
\n\
    --quiet, -q\n\
        Proceed quietly.\n\
\n\
    --help, -h\n\
        Print this help message.\n\
\n\
\n\
 (Experimental options)\n\
    --windows, -w <size,shift>\n\
        Apply a sliding window to the alignment, and fit a separate\n\
        tree to each window.  Arguments specify size of window and\n\
        amount by which to shift it on each iteration, both in number\n\
        of columns.  Separate versions of all output files will be\n\
        created for each window (file name will include window number).\n\
\n\
    --windows-explicit, -v <window_coord_list>\n\
        Like --windows, except that all start and end coordinates must\n\
        be explicitly specified.  Each successive pair of numbers is\n\
        interpreted as defining the start and end of a window.\n\
\n\
    --coord-frame, -d <coord_frame>\n\
        Sequence defining frame of reference for all coordinates.  Default \n\
        is 1.  Use 0 for the entire multiple alignment.\n\
\n\
    --feature-counts, -n <nclasses>\n\
        (for use with --windows or --windows-explicit and --features.)\n\
        Count the number of sites in each window corresponding to each\n\
        class of feature, as designated in the GFF.  The total number\n\
        of classes must be specified as an argument.  Counts will be\n\
        included in summary file.  Only informative sites will be\n\
        counted (sites with at least two non-gap characters).\n\
\n\
    --window-map, -y <lav_fname>\n\
        (for use with --windows or --windows-explicit).  Transform\n\
        window coordinates using the specified local alignment (LAV\n\
        format).  Raw coordinates are assumed to refer to the query\n\
        sequence, and coordinates will be transformed to refer to the\n\
        target sequence, which should be present in the multiple\n\
        alignment being analyzed (see --coord-frame).  With --windows,\n\
        raw coordinates will be generated to span the entire query\n\
        sequence.\n\
\n\
    --max-samples, -S <max_samples>\n\
        Maximum number of sites to consider when fitting any model.\n\
        If the number of eligible columns (considering category\n\
        designations) is greater than this number, max_samples of them\n\
        will be selected randomly.\n\
\n\
    --ancestor, -A <seqname>\n\
        Treat specified sequence as the root of the tree.  The tree\n\
        topology must define this sequence to be a child of the root\n\
        (in practice, the branch from the root to the specified\n\
        sequence will be retained, but will be constrained to have\n\
        length zero).\n\
\n\
    --post-probs, -P\n\
        Output posterior probabilities of all bases at all ancestral \n\
        nodes.  Output will be to auxiliary file(s) with suffix \n\
        \".postprobs\".\n\
\n\
    --expected-subs, -X\n\
        Output posterior expected number of substitutions on each branch at\n\
        each site, summed across all types of substitutions. \n\
        Output will be to auxiliary file(s) with suffix \".expsub\".\n\
\n\
    --expected-total-subs, -Z\n\
        Output posterior expected number of substitutions of each type \n\
        on each branch, summed across all sites.  Output will be to \n\
        auxiliary file(s) with suffix \".exptotsub\".\n\
\n\
    --column-probs, -U\n\
        (for use with -init-model; implies --lnl)  Output a separate log\n\
        probability for each type of column in the input.  Output will\n\
        be to a file with suffix \".colprobs\".  Values are log base 2.\n\
\n\
    --rate-constants, -K <rate_consts>\n\
        Use a non-parameteric mixture model for rates, instead of\n\
        assuming a gamma distribution.  The argument <rate_consts>\n\
        must be a comma-delimited list explicitly defining the rate\n\
        constants to be used.  The \"weight\" (mixing proportion)\n\
        associated with each rate constant will be estimated by EM\n\
        (this option implies --EM).  If --alpha is used with\n\
        this option, then the mixing proportions will be initialized\n\
        to reflect a gamma distribution with the specified shape\n\
        parameter.\n\
\n\
\n\
REFERENCES:\n\
\n\
    A. Siepel and D. Haussler.  2004.  Phylogenetic estimation of\n\
      context-dependent substitution rates by maximum likelihood.\n\
      Mol. Biol. Evol., 21:468-488.\n\
\n\
    Z. Yang, N. Goldman, and A. Friday.  1994. Comparison of models for\n\
      nucleotide substution used in maximum likelihood phylogenetic\n\
      estimation. Mol. Biol. Evol., 11:316-324.\n\
\n\
    Z. Yang. 1994. Maximum likelihood phylogenetic estimation from\n\
      DNA sequences with variable rates over sites: approximate\n\
      methods. J. Mol. Evol., 39:306-314.\n\n", DEFAULT_NSITES_THRESHOLD);
}

void set_output_fname(String *fname, char *root, int cat, char *suffix) {
  str_cpy_charstr(fname, root);
  if (cat != -1) {
    str_append_charstr(fname, ".");
    str_append_int(fname, cat);
  }
  str_append_charstr(fname, suffix);
}

/* Compute and output statistics based on posterior probabilities,
   including (optionally) the post prob of each tuple of bases at each
   ancestral node at each site (do_bases), the expected total number
   of substs per site (do_expected_nsubst), and the expected number of
   substitutions of each type on each edge across all sites
   (do_expected_nsubst_tot).  A separate file is output for each
   selected option, with an appropriate filename suffix (".postprob",
   ".expsub", and ".exptotsub", respectively).  */
void print_post_prob_stats(TreeModel *mod, MSA *msa, char *output_fname_root, 
                           int do_bases, int do_expected_nsubst, 
                           int do_expected_nsubst_tot, int cat, int quiet) {
  String *fname = str_new(STR_MED_LEN);
  FILE *POSTPROBF, *EXPSUBF, *EXPTOTSUBF;
  int i, tup, node, state, state2;
  TreeNode *n;
  char tuplestr[mod->order+2];
  char coltupstr[msa->nseqs+1];
  tuplestr[mod->order+1] = '\0';
  coltupstr[msa->nseqs] = '\0';

  /* FIXME: rate variation!  need rate post probs! */
  assert(mod->nratecats == 1);

  /* compute desired stats */
  assert(mod->tree_posteriors == NULL); 
  mod->tree_posteriors = tl_new_tree_posteriors(mod, msa, do_bases, 0, 
                                                do_expected_nsubst, 
                                                do_expected_nsubst_tot, 0, 0);
  tl_compute_log_likelihood(mod, msa, NULL, cat, mod->tree_posteriors);

  if (do_bases) {
    set_output_fname(fname, output_fname_root, cat, ".postprob");
    if (!quiet) 
      fprintf(stderr, "Writing posterior probabilities to %s ...\n", 
              fname->chars);
    POSTPROBF = fopen_fname(fname->chars, "w+");

    /* print header */
    fprintf(POSTPROBF, "%-6s ", "#");
    for (i = 0; i < msa->nseqs; i++) fprintf(POSTPROBF, " ");
    fprintf(POSTPROBF, "    ");
    for (node = 0; node < mod->tree->nnodes; node++) {
      n = lst_get_ptr(mod->tree->nodes, node);
      if (n->lchild == NULL || n->rchild == NULL) continue;
      for (state = 0; state < mod->rate_matrix->size; state++) {
        if (state == mod->rate_matrix->size/2)
          fprintf(POSTPROBF, "node %-2d", n->id);
        else
          fprintf(POSTPROBF, "%6s ", "");
      }
    }
    fprintf(POSTPROBF, "\n%-6s ", "#");
    for (state = 0; state < msa->nseqs-5; state++) fprintf(POSTPROBF, " ");
    fprintf(POSTPROBF, "tuple    ");
    for (node = 0; node < mod->tree->nnodes; node++) {
      n = lst_get_ptr(mod->tree->nodes, node);
      if (n->lchild == NULL || n->rchild == NULL) continue;
      for (state = 0; state < mod->rate_matrix->size; state++) {
        get_tuple_str(tuplestr, state, mod->order + 1, 
                      mod->rate_matrix->states);
        fprintf(POSTPROBF, "%6s ", tuplestr);
      }
    }
    fprintf(POSTPROBF, "\n");

    /* print post probs */
    for (tup = 0; tup < msa->ss->ntuples; tup++) {

      if ((cat >= 0 && msa->ss->cat_counts[cat][tup] == 0) ||
          msa->ss->counts[tup] == 0) continue;

      tuple_to_string_pretty(coltupstr, msa, tup);
      fprintf(POSTPROBF, "%-6d %5s    ", tup, coltupstr);
      for (node = 0; node < mod->tree->nnodes; node++) {
        n = lst_get_ptr(mod->tree->nodes, node);
        if (n->lchild == NULL || n->rchild == NULL) continue;
        for (state = 0; state < mod->rate_matrix->size; state++) 
          fprintf(POSTPROBF, "%6.4f ", 
                  mod->tree_posteriors->base_probs[0][state][node][tup]);
      }              
      fprintf(POSTPROBF, "\n");
    }
    fclose(POSTPROBF);
  }

  if (do_expected_nsubst) {
    set_output_fname(fname, output_fname_root, cat, ".expsub");
    if (!quiet) 
      fprintf(stderr, "Writing expected numbers of substitutions to %s ...\n", 
              fname->chars);
    EXPSUBF = fopen_fname(fname->chars, "w+");

    fprintf(EXPSUBF, "%-3s %10s ", "#", "tuple");
    for (node = 0; node < mod->tree->nnodes; node++) {
      n = lst_get_ptr(tr_postorder(mod->tree), node);
      if (n == mod->tree) continue;
      fprintf(EXPSUBF, " node_%-2d", n->id);
    }
    fprintf(EXPSUBF, "  total\n");
    for (tup = 0; tup < msa->ss->ntuples; tup++) {
      double total = 0;

      if ((cat >= 0 && msa->ss->cat_counts[cat][tup] == 0) ||
          msa->ss->counts[tup] == 0) continue;

      tuple_to_string_pretty(coltupstr, msa, tup);
      fprintf(EXPSUBF, "%-3d %10s ", tup, coltupstr);
      for (node = 0; node < mod->tree->nnodes; node++) {
        n = lst_get_ptr(tr_postorder(mod->tree), node);
        if (n == mod->tree) continue;
        fprintf(EXPSUBF, "%7.4f ", 
                mod->tree_posteriors->expected_nsubst[0][n->id][tup]);
        total += mod->tree_posteriors->expected_nsubst[0][n->id][tup];
      }              
      fprintf(EXPSUBF, "%7.4f\n", total);
    }
    fclose(EXPSUBF);
  }

  if (do_expected_nsubst_tot) {
    set_output_fname(fname, output_fname_root, cat, ".exptotsub");
    if (!quiet) 
      fprintf(stderr, "Writing total expected numbers of substitutions to %s ...\n", 
              fname->chars);
    EXPTOTSUBF = fopen_fname(fname->chars, "w+");

    fprintf(EXPTOTSUBF, "\n\
A separate matrix of expected numbers of substitutions is shown for each\n\
branch of the tree.  Nodes of the tree are visited in a postorder traversal,\n\
and each node is taken to be representative of the branch between itself and\n\
its parent.  Starting bases or tuples of bases appear on the vertical axis\n\
of each matrix, and destination bases or tuples of bases appear on the\n\
horizontal axis.\n\n");

    for (node = 0; node < mod->tree->nnodes; node++) {
      n = lst_get_ptr(tr_postorder(mod->tree), node);
      if (n == mod->tree) continue;

      fprintf(EXPTOTSUBF, "Branch above node %d", n->id);
      if (n->name != NULL && strlen(n->name) > 0) 
        fprintf(EXPTOTSUBF, " (leaf labeled '%s')", n->name);
      fprintf(EXPTOTSUBF, ":\n\n");

      /* print header */
      fprintf(EXPTOTSUBF, "%-4s ", "");
      for (state2 = 0; state2 < mod->rate_matrix->size; state2++) {
        get_tuple_str(tuplestr, state2, mod->order + 1, 
                      mod->rate_matrix->states);
        fprintf(EXPTOTSUBF, "%12s ", tuplestr);
      }
      fprintf(EXPTOTSUBF, "\n");
      for (state = 0; state < mod->rate_matrix->size; state++) {
        get_tuple_str(tuplestr, state, mod->order + 1, 
                      mod->rate_matrix->states);
        fprintf(EXPTOTSUBF, "%-4s ", tuplestr);
        for (state2 = 0; state2 < mod->rate_matrix->size; state2++) 
          fprintf(EXPTOTSUBF, "%12.2f ", 
                  mod->tree_posteriors->expected_nsubst_tot[0][state][state2][n->id]);
        fprintf(EXPTOTSUBF, "\n");
      }
      fprintf(EXPTOTSUBF, "\n\n");
    }
    fclose(EXPTOTSUBF);
  }

  tl_free_tree_posteriors(mod, msa, mod->tree_posteriors);
  mod->tree_posteriors = NULL;
  str_free(fname);
}

int main(int argc, char *argv[]) {
  char *msa_fname = NULL, *gff_fname = NULL,
    *cat_map_fname = NULL, *output_fname_root = "ftm", 
    *log_fname = NULL, *coord_aln_fname = NULL, *input_mod_fname = NULL;
  int output_trees = 0, input_format = PHYLIP, subst_mod = REV, quiet = 0,
    gap_strip_mode = NO_STRIP, 
    nratecats = 1, max_samples = -1, reverse_complement = 1, ncats = 0, 
    use_em = 0, window_size = -1, window_shift = -1, coord_frame = 1,
    nclasses = -1, use_conditionals = 0, 
    precision = OPT_HIGH_PREC, cycle_size = -1,
    likelihood_only = 0, do_bases = 0, do_expected_nsubst = 0, 
    do_expected_nsubst_tot = 0, nsites_threshold = DEFAULT_NSITES_THRESHOLD,
    random_init = 0, estimate_backgd = 0, estimate_scale_only = 0,
    do_column_probs;
  char c;
  FILE *F, *WINDOWF;
  TreeNode *tree = NULL;
  CategoryMap *cm = NULL;
  int i, j, win, tbases, opt_idx;
  String *s, *mod_fname, *out_tree_fname, *root_seqname = NULL;
  MSA *msa, *source_msa;
  FILE *logf = NULL;
  String *tmpstr = str_new(STR_SHORT_LEN);
  List *cats_to_do = NULL, *weight_matrix_list = NULL,
    *tmplist = NULL, *window_coords = NULL, *match_list = NULL,
    *cats_to_do_str = NULL;
  int *weight_matrix, *nsites, *totbases;
  double *gc, cpg;
  double alpha = DEFAULT_ALPHA;
  LocalPwAlignment *lpwa = NULL;
  GFF_Set *gff = NULL;
  TreeModel *input_mod = NULL;
  int root_leaf_id = -1;
  Regex *class_re = str_re_new(".*category[[:space:]]+[A-Za-z0-9_]+[[:space:]]+\\[([[:digit:]]+)\\]");
  List *rate_consts = NULL;

  struct option long_opts[] = {
    {"msa", 1, 0, 'm'},
    {"tree", 1, 0, 't'},
    {"subst-mod", 1, 0, 's'},
    {"msa-format", 1, 0, 'i'},
    {"nrates", 1, 0, 'k'},
    {"alpha", 1, 0, 'a'},
    {"features", 1, 0, 'g'},
    {"catmap", 1, 0, 'c'},
    {"log", 1, 0, 'l'},
    {"out-root", 1, 0, 'o'},
    {"output-tree", 0, 0, 'T'},
    {"EM", 0, 0, 'E'},
    {"precision", 1, 0, 'p'},
    {"do-cats", 1, 0, 'C'},
    {"cats-cycle", 1, 0, 'Y'},
    {"markov-dependence", 0, 0, 'N'},
    {"gap-strip", 1, 0, 'G'},
    {"weight-matrix-cats", 1, 0, 'W'},
    {"no-reverse", 0, 0, 'R'},
    {"init-model", 1, 0, 'M'},
    {"init-random", 0, 0, 'r'},
    {"lnl", 0, 0, 'L'},
    {"scale-only", 0, 0, 'B'},
    {"estimate-freqs", 0, 0, 'F'},
    {"min-informative", 1, 0, 'I'},
    {"quiet", 0, 0, 'q'},
    {"help", 0, 0, 'h'},
    {"windows", 1, 0, 'w'},
    {"windows-explicit", 1, 0, 'v'},
    {"coord-frame", 1, 0, 'd'},
    {"feature-counts", 1, 0, 'n'},
    {"window-map", 1, 0, 'y'},
    {"max-samples", 1, 0, 'S'},
    {"ancestor", 1, 0, 'A'},
    {"post-probs", 1, 0, 'P'},
    {"expected-subs", 1, 0, 'X'},
    {"expected-total-subs", 1, 0, 'Z'},
    {"column-probs", 1, 0, 'U'},
    {"rate-constants", 1, 0, 'K'},
  };

  while ((c = getopt_long(argc, argv, "m:t:s:g:c:C:Y:G:i:o:k:a:l:S:W:w:v:y:d:n:M:p:A:I:K:EeNDRTqLPXZUBFrh", long_opts, &opt_idx)) != -1) {
    switch(c) {
    case 'm':
      msa_fname = optarg;
      break;
    case 't':
      if (optarg[0] == '(')     /* in this case, assume topology given
                                   at command line */
        tree = parse_nh_from_string(str_new_charstr(optarg));
      else {
        if (!quiet) fprintf(stderr, "Reading tree from %s ...\n", optarg);
        tree = parse_nh_from_file(fopen_fname(optarg, "r"));
      }
      break;
    case 's':
      subst_mod = tm_get_subst_mod_type(s = str_new_charstr(optarg));
      if (subst_mod == UNDEF_MOD) 
        die("ERROR: illegal substitution model.  Type \"phyloFit -h\" for usage.\n");
      break;
    case 'g':
      gff_fname = optarg;
      break;
    case 'c':
      cat_map_fname = optarg;
      break;
    case 'C':
      cats_to_do_str = get_arg_list(optarg);
      break;
    case 'Y':
      cycle_size = atoi(optarg);
      break;
    case 'o':
      output_fname_root = optarg;
      break;
    case 'T':
      output_trees = 1;
      break;
    case 'k':
      nratecats = atoi(optarg);
      if (nratecats <= 0) 
        die("ERROR: number of rate categories must be >= 1.\n");
      break;
    case 'a':
      alpha = atof(optarg);
      break;
    case 'R':
      reverse_complement = 0;
      break;
    case 'G':
      if (!strcmp(optarg, "ALL")) gap_strip_mode = STRIP_ALL_GAPS;
      else if (!strcmp(optarg, "ANY")) gap_strip_mode = STRIP_ANY_GAPS;
      else gap_strip_mode = atoi(optarg);
      break;        
    case 'i':
      if (!strcmp(optarg, "PSU")) input_format = PSU;
      else if (!strcmp(optarg, "FASTA")) input_format = FASTA;
      else if (!strcmp(optarg, "SS")) input_format = SS;
      else if (!strcmp(optarg, "LAV")) input_format = LAV;
      else if (strcmp(optarg, "PHYLIP") != 0) 
        die("ERROR: unrecognized alignment format.  Type 'phyloFit -h' for usage.\n");
      break;
    case 'l':
      log_fname = optarg;
      if (!strcmp(log_fname, "-"))
        logf = stderr;
      logf = fopen_fname(log_fname, "w+");
      break;
    case 'N':
      use_conditionals = 1;
      break;
    case 'S':
      max_samples = atoi(optarg);
      if (max_samples <= 0) 
        die("ERROR: maximum number of samples (--max-samples) must be greater than zero.\n");
      break;
    case 'W':
      weight_matrix_list = get_arg_list(optarg);
      break;
    case 'w':
      tmplist = get_arg_list(optarg);
      if (lst_size(tmplist) != 2 ||
          str_as_int(lst_get_ptr(tmplist, 0), &window_size) != 0 ||
          str_as_int(lst_get_ptr(tmplist, 1), &window_shift) != 0) 
        die("ERROR: illegal arguments to --windows.\n");
      lst_free_strings(tmplist);
      lst_free(tmplist);
      break;
    case 'v':
      tmplist = get_arg_list(optarg);
      if (lst_size(window_coords) % 2 != 0) 
        die("ERROR: argument to -v must be a list of even length.\n");
      window_coords = lst_new_int(lst_size(tmplist));
      for (i = 0; i < lst_size(tmplist); i++) {
        int tmp;
        if (str_as_int(lst_get_ptr(tmplist, i), &tmp) != 0) 
          die("ERROR: illegal argument to --windows-explicit.\n");
        lst_push_int(window_coords, tmp);
        str_free(lst_get_ptr(tmplist, i));
      }
      lst_free(tmplist);
      break;
    case 'n':
      nclasses = atoi(optarg);
      break;
    case 'd':
      coord_frame = atoi(optarg);
      assert(0);                /* not yet implemented! */
      break;
    case 'y':
      coord_aln_fname = optarg;
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
      input_mod_fname = optarg;
      break;
    case 'r':
      random_init = 1;
      break;
    case 'L':
      likelihood_only = 1;
      break;
    case 'q':
      quiet = 1;
      break;
    case 'P':
      do_bases = 1;
      break;
    case 'X':
      do_expected_nsubst = 1;
      break;
    case 'Z':
      do_expected_nsubst_tot = 1;
      break;
    case 'U':
      likelihood_only = 1;      /* force -L */
      nsites_threshold = 0;     /* also force this; typical use is
                                   with small number of tuples, no
                                   tuple_idx */
      do_column_probs = 1;
      break;
    case 'A':
      root_seqname = str_new_charstr(optarg);
      break;
    case 'I':
      nsites_threshold = atoi(optarg);
      break;
    case 'B':
      estimate_scale_only = 1;
      break;
    case 'F':
      estimate_backgd = 1;
      break;
    case 'K':
      tmplist = get_arg_list(optarg);
      rate_consts = str_list_as_dbl(tmplist);
      lst_qsort_dbl(rate_consts, ASCENDING);
      if (lst_size(rate_consts) < 2 || lst_get_dbl(rate_consts, 0) <= 0) 
        die("ERROR: must be >= 2 rate constants and all must be positive.\n");
      nratecats = lst_size(rate_consts);
      use_em = 1;
      lst_free_strings(tmplist); lst_free(tmplist);
      break;
    case 'h':
      print_usage();
      exit(0);
    case '?':
      die("ERROR: illegal argument.  Type 'phyloFit -h' for usage.\n");
    }
  }

  if (msa_fname == NULL || (tree == NULL && input_mod_fname == NULL)) 
    die("ERROR: must specify --msa and either --tree or --init-model.  Type 'phyloFit -h' for usage.\n");

  if (gff_fname != NULL && cat_map_fname == NULL) 
    die("ERROR: --features requires --catmap.  Type 'phyloFit -h' for usage.\n");

  if (use_conditionals && use_em) 
    die("ERROR: Cannot use --markov-dependence with --EM.  Type 'phyloFit -h' for usage.\n");

  if (cycle_size != -1 && gff_fname != NULL) 
    die("ERROR: Cannot use both --cats-cycle and --features.  Type 'phyloFit -h' for usage.\n");

  if (likelihood_only && input_mod_fname == NULL) 
    die("ERROR: --lnl requires --init-model.  Type 'phyloFit -h' for usage.\n");

  /* read input model */
  if (input_mod_fname != NULL) {
    if (!quiet) fprintf(stderr, "Reading tree model from %s ...\n", 
                        input_mod_fname);
    input_mod = tm_new_from_file(fopen_fname(input_mod_fname, "r"));
  }
  
  /* read alignment */
  if (!quiet) fprintf(stderr, "Reading alignment from %s ...\n", msa_fname);
  msa = msa_new_from_file(fopen_fname(msa_fname, "r"), input_format, "ACGT");
                                /* FIXME: need better strategy for Ns */

  /* make sure alignment and tree topology consistent */
  if (msa->nseqs * 2 - 1 != 
      (input_mod == NULL ? tree->nnodes : input_mod->tree->nnodes)) 
    die("ERROR: Tree must have 2n-1 nodes, where n is the number of sequences in the\nalignment.  Even with a reversible model, specify a rooted tree; the two\nbranches adjoining the root will simply be assigned a single parameter.\n");

  /* allow for specified ancestor */
  if (root_seqname != NULL) {
    char tmpstr[5];
    int idx;
    TreeNode *rl = NULL;
    if (tree == NULL || tm_is_reversible(subst_mod)) 
      die("ERROR: --ancestor requires --tree and a non-reversible model.\n");
    idx = msa_get_seq_idx(msa, root_seqname);
    sprintf(tmpstr, "%d", idx);
    rl = tr_get_node(tree, tmpstr);
    
    if (rl == NULL || rl->parent != tree) 
      die("ERROR: Sequence specified by --ancestor must be a child of the root.\n");
    
    root_leaf_id = rl->id;
  }

  /* temporary: in case of SS or LAV format, remove 'N' from alphabet
     if necessary (model fitting code will treat Ns as missing data;
     we don't want to set up with 'N' as a legitimate character in the
     alphabet) */
  if (input_format == SS || input_format == LAV) 
    msa_remove_N_from_alph(msa);

  /* read coordinate-transforming alignment, if necessary */
  if (coord_aln_fname != NULL) 
    lpwa = la_read_lav(fopen_fname(coord_aln_fname, "r"), 0);
  
  if (cat_map_fname != NULL) {
    /* read category map */
    if (!quiet) fprintf(stderr, "Reading category map from %s ...\n", 
                        cat_map_fname);
    cm = cm_read(fopen_fname(cat_map_fname, "r"));
    ncats = cm->ncats;
  }

  if (gff_fname != NULL) {      /* GFF specified: label categories */
    /* read gff */
    if (!quiet) fprintf(stderr, "Reading annotations from %s ...\n", 
                        gff_fname);
    gff = gff_read_set(fopen_fname(gff_fname, "r"));

    /* convert GFF to coordinate frame of alignment */
    msa_map_gff_coords(msa, gff, -1, 0, 0, NULL);

    /* reverse complement segments of MSA corresponding to features on
       reverse strand (if necessary) */
    if (reverse_complement) 
      msa_reverse_compl_gff(msa, gff, NULL);

    /* label categories */
    if (!quiet) fprintf(stderr, "Labeling alignment sites by category ...\n");
    msa_label_categories(msa, gff, cm);
  }
  else if (msa->ss != NULL && msa->ncats > 0)
    ncats = msa->ncats;
                                /* in this case, categories are
                                   specified indirectly, via
                                   category-specific counts in an SS
                                   file */
  else if (cycle_size != -1) {
    assert(cycle_size > 0);
    msa->categories = (int*)smalloc(msa->length * sizeof(int));
    for (i = 0; i < msa->length; i++) 
      msa->categories[i] = (i % cycle_size) + 1;
    ncats = msa->ncats = cycle_size;
  }

  /* strip gaps, if necessary.  Category labels will be adjusted
     appropriately */
  if (gap_strip_mode != NO_STRIP)
    strip_gaps(msa, gap_strip_mode);

  mod_fname = str_new(STR_MED_LEN);
  if (output_trees) out_tree_fname = str_new(STR_MED_LEN);

  /* set up list of categories to process */
  if (cm == NULL && cycle_size == -1 && ncats == 0) {
    cats_to_do = lst_new_int(1);
    lst_push_int(cats_to_do, -1);
  }
  else if (cats_to_do_str == NULL) {
    cats_to_do = lst_new_int(ncats+1);
    for (i = 0; i <= ncats; i++)
      lst_push_int(cats_to_do, i);
  }
  else {
    if (cm == NULL) {           /* e.g., with SS file and cat-specific
                                   counts */
      cats_to_do = lst_new_int(lst_size(cats_to_do_str));
      for (i = 0; i < lst_size(cats_to_do_str); i++) {
        int tmpint;
        if (str_as_int(lst_get_ptr(cats_to_do_str, i), &tmpint) != 0) 
          die("ERROR: categories list must consist of integers.\n");
        lst_push_int(cats_to_do, tmpint);
      }
    }
    else 
      cats_to_do = cm_get_category_list(cm, cats_to_do_str, 0);
  }

  /* identify weight-matrix categories */
  weight_matrix = (int*)smalloc((ncats+1) * sizeof(int));
  for (i = 0; i <= ncats; i++) weight_matrix[i] = 0;
  if (weight_matrix_list != NULL) {
    if (lst_size(weight_matrix_list) == 1 && 
        str_equals_nocase_charstr(lst_get_ptr(weight_matrix_list, 0), "all")) 
      for (i = 0; i <= ncats; i++) weight_matrix[i] = 1;
    else {
      List *weight_matrix_cats = 
        cm_get_category_list(cm, weight_matrix_list, 0);
      for (i = 0; i < lst_size(weight_matrix_cats); i++)
        weight_matrix[lst_get_int(weight_matrix_cats, i)] = 1;
      lst_free(weight_matrix_cats);
    }
    lst_free_strings(weight_matrix_list);
    lst_free(weight_matrix_list);
  }

  /* set up windows, if necessary */
  if (window_size != -1) {
    int len;
    if (window_coords != NULL) 
      die("ERROR: cannot use both --windows and --windows-explicit.\n");
    len = (lpwa != NULL ? lpwa->query_len : msa->length);
    window_coords = lst_new_int(len/window_shift + 1);
    for (i = 1; i < len; i += window_shift) {
      lst_push_int(window_coords, i);
      lst_push_int(window_coords, 
                   min(i + window_size - 1, len));
    }
  }

  if (window_coords != NULL) {
    /* set up summary file */
    String *sumfname = str_new_charstr(output_fname_root);
    char tmpstr[10];

    str_append_charstr(sumfname, ".win-sum");
    WINDOWF = fopen_fname(sumfname->chars, "w+");
    str_free(sumfname);

    fprintf(WINDOWF, "%5s %8s %8s %4s", "win", "beg", "end", "cat");
/*     for (j = 0; j < strlen(msa->alphabet); j++) */
/*       fprintf(WINDOWF, " %6c", msa->alphabet[j]); */

    fprintf(WINDOWF, " %6s", "GC");
    for (j = 0; nclasses > 0 && j < msa->nseqs; j++) {
      sprintf(tmpstr, "GC%d", j+1);
      fprintf(WINDOWF, " %6s", tmpstr);
    }
    fprintf(WINDOWF, " %8s", "CpG");
    fprintf(WINDOWF, " %7s", "ninf");
    for (j = 0; j < nclasses; j++) { /* note: will only be executed with -n */
      sprintf(tmpstr, "ninf%d", j+1);
      fprintf(WINDOWF, " %7s", tmpstr);
    }
    fprintf(WINDOWF, " %7s\n", "t");

    /* transform coordinates, if necessary */
    if (lpwa != NULL) {
      for (i = 0; i < lst_size(window_coords); i += 2) {
        int win_beg = la_get_target_coord(lpwa, lst_get_int(window_coords, i), 
                                          ADJUSTRIGHT);
        int win_end = la_get_target_coord(lpwa, lst_get_int(window_coords, 
                                                            i+1), ADJUSTLEFT);
        if (win_beg == -1 || win_end == -1 || win_beg >= win_end)
          win_beg = win_end = -1;
        lst_set_int(window_coords, i, win_beg);
        lst_set_int(window_coords, i+1, win_end);
      }
      /* FIXME: need a way to free lpwa! */
    }

    if (coord_frame != 0) { 
      msa_coord_map *map = msa_build_coord_map(msa, coord_frame);
      for (i = 0; i < lst_size(window_coords); i += 2) {
        lst_set_int(window_coords, i, 
                    msa_map_seq_to_msa(map, lst_get_int(window_coords, i)));
        lst_set_int(window_coords, i+1, 
                    msa_map_seq_to_msa(map, lst_get_int(window_coords, i+1)));
      }
      msa_map_free(map);
    }
  }

  source_msa = msa;

  if (nclasses > 0) {           /* inits for -n option */
    match_list = lst_new_ptr(2);
    nsites = (int*)smalloc(nclasses * sizeof(int));
    gc = (double*)smalloc(msa->nseqs * sizeof(double));
    totbases = (int*)smalloc(msa->nseqs * sizeof(int));
    gff_sort(gff);
  }

  /* process each window, if necessary */
  for (win = 0; 
       win < (window_coords == NULL ? 1 : lst_size(window_coords)); 
       win += 2) {
    int win_beg, win_end;

    if (window_coords != NULL) {
      win_beg = lst_get_int(window_coords, win);
      win_end = lst_get_int(window_coords, win+1);
      if (win_beg < 0 || win_end < 0) continue;

      /* note: msa_sub_alignment uses a funny indexing system (see docs) */
      msa = msa_sub_alignment(source_msa, NULL, 0, win_beg-1, win_end);
    }

    /* process each category */
    for (i = 0; i < lst_size(cats_to_do); i++) {
      TreeModel *mod;
      gsl_vector *params = NULL;
      int cat = lst_get_int(cats_to_do, i);
      int wm_cat = cat >= 0 ? cat : 0;
      int ninf_sites;

      /* count number of aligned sites in window belonging to each
         class of feature (if necessary) */
      if (nclasses > 0) {
        assert(gff != NULL && cm != NULL && window_coords != NULL);
        for (j = 0; j < nclasses; j++) nsites[j] = 0;
        for (j = 0; j < source_msa->nseqs; j++) gc[j] = totbases[j] = 0;
        cpg = tbases = 0; 
        for (j = 0; j < lst_size(gff->features); j++) {
          int class_num, beg, end, k;
          GFF_Feature *feat = lst_get_ptr(gff->features, j);

          if (cm_get_category(cm, feat->feature) != cat) continue;

          if (feat->end < win_beg) continue;
          else if (feat->start > win_end) 
            break;

          /* obtain class number from feature */
          if (str_re_match(feat->attribute, class_re, match_list, 1) <= 0)
            continue;
          str_as_int(lst_get_ptr(match_list, 1), &class_num);
          str_free(lst_get_ptr(match_list, 0));
          str_free(lst_get_ptr(match_list, 1));

          /* now count the number of aligned bases */
          /* at the same time, obtain the GC content of each seq, and 
             the number of CpGs */
          beg = max(feat->start, win_beg);
          end = min(feat->end, win_end);
          for (k = beg-1; k < end; k++) { /* check index */
            int ninf = 0, seq;
            for (seq = 0; seq < source_msa->nseqs; seq++) {
              if (source_msa->seqs[seq][k] != GAP_CHAR) {
                ninf++;
                totbases[seq]++;
                tbases++;
              }
              if (toupper(source_msa->seqs[seq][k]) == 'G' || 
                  toupper(source_msa->seqs[seq][k]) == 'C')
                gc[seq]++;              
              if (toupper(source_msa->seqs[seq][k]) == 'G' && k > 0 && 
                  toupper(source_msa->seqs[seq][k-1] == 'C'))
                cpg++;
            }
            assert(class_num >= 1 && class_num <= nclasses);
            if (ninf >= 2)
              nsites[class_num-1]++;
          }
        }
        for (j = 0; j < source_msa->nseqs; j++) gc[j] /= totbases[j];
        cpg /= tbases;
      }

      if (input_mod == NULL) {
        mod = tm_new(weight_matrix[wm_cat] == 0 ? tr_create_copy(tree) : NULL, 
                     NULL, NULL, subst_mod, msa->alphabet, nratecats, alpha, 
                     rate_consts, root_leaf_id);
      }
      else if (likelihood_only)
        mod = input_mod;
      else {
        double newalpha = 
          (input_mod->nratecats > 1 && alpha == DEFAULT_ALPHA ? 
           input_mod->alpha : alpha);
                                /* if the input_mod has a meaningful
                                   alpha and a non-default alpha has
                                   not been specified, then use
                                   input_mod's alpha  */
        mod = input_mod;
        tm_reinit(mod, subst_mod, nratecats, newalpha, rate_consts);
      }

      mod->use_conditionals = use_conditionals;

      if (estimate_scale_only || estimate_backgd) {
        tm_free_rmp(mod);
        if (estimate_scale_only) mod->estimate_branchlens = TM_SCALE_ONLY;
        mod->estimate_backgd = estimate_backgd;
        tm_init_rmp(mod);       /* necessary because number of
                                   parameters changes */
      }

      if (!quiet) {
        str_clear(tmpstr);

        str_append_charstr(tmpstr, msa_fname);
        if (cat != -1 || window_coords != NULL) {
          str_append_charstr(tmpstr, " (");
          if (cat != -1) {
            str_append_charstr(tmpstr, "category ");
            str_append_int(tmpstr, cat);
          }
          if (window_coords != NULL) {
            if (cat != -1) str_append_charstr(tmpstr, ", ");
            str_append_charstr(tmpstr, "window ");
            str_append_int(tmpstr, win/2 + 1);
          }
          str_append_char(tmpstr, ')');
        }
      }

      ninf_sites = msa_ninformative_sites(msa, cat);
      if (!weight_matrix[wm_cat] && ninf_sites < nsites_threshold) {
        tm_free(mod);
        fprintf(stderr, "Skipping %s; insufficient informative sites ...\n", 
                tmpstr->chars);
        continue;
      }

      if (likelihood_only) {
        double *col_log_probs = do_column_probs ? 
          smalloc(msa->length * sizeof(double)) : NULL;
        String *colprob_fname;
        if (!quiet) 
          fprintf(stderr, "Computing likelihood of %s ...\n", tmpstr->chars);
        tm_set_subst_matrices(mod);
        if (do_column_probs && msa->ss != NULL && msa->ss->tuple_idx == NULL) {
          msa->ss->tuple_idx = smalloc(msa->length * sizeof(int));
          for (j = 0; j < msa->length; j++)
            msa->ss->tuple_idx[j] = j;
        }
        mod->lnL = tl_compute_log_likelihood(mod, msa, col_log_probs, cat, NULL) * 
          log(2);
        if (do_column_probs) {
          colprob_fname = str_new_charstr(output_fname_root);
          str_append_charstr(colprob_fname, ".colprobs");
          if (!quiet) 
            fprintf(stderr, "Writing column probabilities to %s ...\n", 
                    colprob_fname->chars);
          F = fopen_fname(colprob_fname->chars, "w+");
          for (j = 0; j < msa->length; j++)
            fprintf(F, "%d\t%.6f\n", j, col_log_probs[j]);
          fclose(F);
          str_free(colprob_fname);
          free(col_log_probs);
        }
      }
      else {                    /* fit model */
        if (weight_matrix[wm_cat] == 0) {
          if (random_init) 
            params = tm_params_init_random(mod);
          else if (input_mod != NULL)
            params = tm_params_init_from_model(input_mod);
          else
            params = tm_params_init(mod, .1, 5, alpha);    
          if (input_mod != NULL && mod->backgd_freqs != NULL) {
            /* in some cases, these are needed for initialization, but
               now they should be re-estimated */
            gsl_vector_free(mod->backgd_freqs);
            mod->backgd_freqs = NULL;
          }
        }

        if (!quiet) {
          fprintf(stderr, "Fitting tree model to %s using %s%s ...\n",
                  tmpstr->chars, weight_matrix[wm_cat] ? "weight matrix" : 
                  tm_get_subst_mod_string(subst_mod),
                  !weight_matrix[wm_cat] && mod->nratecats > 1 ? " (with rate variation)" : "");
          if (log_fname != NULL)
            fprintf(stderr, "(writing log to %s)\n", log_fname);
        }

        if (msa->ss == NULL) 
          /* ensures that sufficient stats are obtained only for
             categories of interest */
          ss_from_msas(msa, mod->order+1, 0, 
                       cats_to_do_str != NULL ? cats_to_do : NULL, 
                       NULL, NULL, -1);

        if (use_em && !weight_matrix[wm_cat])
          tm_fit_em(mod, msa, params, cat, precision, logf);
        else
          tm_fit(mod, msa, params, cat, max_samples, precision, logf);
      }


      str_cpy_charstr(mod_fname, output_fname_root);
      if (window_coords != NULL) {
        str_append_charstr(mod_fname, ".win-");
        str_append_int(mod_fname, win/2 + 1);
      }
      if (cat != -1) {
        str_append_char(mod_fname, '.');
        if (cm != NULL) 
          str_append(mod_fname, cm_get_feature_unique(cm, cat));
        else 
          str_append_int(mod_fname, cat);
      }
      str_append_charstr(mod_fname, ".mod");

      if (!quiet) fprintf(stderr, "Writing model to %s ...\n", 
                          mod_fname->chars);
      F = fopen_fname(mod_fname->chars, "w+");
      tm_print(F, mod);
      fclose(F);

      /* output posterior probabilities, if necessary */
      if (do_bases || do_expected_nsubst || do_expected_nsubst_tot) 
        print_post_prob_stats(mod, msa, output_fname_root, do_bases, 
                              do_expected_nsubst, do_expected_nsubst_tot, 
                              cat, quiet);

      /* also print tree, if requested */
      if (!weight_matrix[wm_cat] && output_trees) {
        /* first create a copy, with leaves labeled with sequence names */
        TreeNode *trcpy = tr_create_copy(mod->tree);
        for (j = 0; j < lst_size(trcpy->nodes); j++) {
          int id;
          TreeNode *n = lst_get_ptr(trcpy->nodes, j);
          if (n->lchild != NULL) continue; /* internal node */
          id = atoi(n->name);
          strcpy(n->name, msa->names[id-1]);
        }

        /* now print altered tree */
        str_cpy_charstr(out_tree_fname, output_fname_root);
        if (cat != -1) {
          str_append_charstr(out_tree_fname, ".");
          str_append_int(out_tree_fname, cat);
        }
        str_append_charstr(out_tree_fname, ".nh");
        if (!quiet) fprintf(stderr, "Writing tree to %s ...\n", 
                            out_tree_fname->chars);
        F = fopen_fname(out_tree_fname->chars, "w+");
        print_tree(F, trcpy, 1);
        fclose(F);
        free_tree(trcpy);
      }

      /* print window summary, if window mode */
      if (window_coords != NULL) {
        fprintf(WINDOWF, "%5d %8d %8d %4d", win/2+1, lst_get_int(window_coords, win), lst_get_int(window_coords, win+1), cat);
/*         for (j = 0; j < mod->backgd_freqs->size; j++) */
/*           fprintf(WINDOWF, " %6.4f", gsl_vector_get(mod->backgd_freqs, j)); */
        fprintf(WINDOWF, " %6.4f", 
                gsl_vector_get(mod->backgd_freqs, 
                               mod->rate_matrix->inv_states[(int)'G']) + 
                gsl_vector_get(mod->backgd_freqs, 
                               mod->rate_matrix->inv_states[(int)'C']));
        for (j = 0; j < msa->nseqs; j++)
          fprintf(WINDOWF, " %6.4f", gc[j]);
        fprintf(WINDOWF, " %8.6f", cpg);
        fprintf(WINDOWF, " %7d", ninf_sites);
        for (j = 0; j < nclasses; j++) /* note: will only be executed with -n */
          fprintf(WINDOWF, " %7d", nsites[j]);
        fprintf(WINDOWF, " %7.4f\n", tr_total_len(mod->tree));
      }

      if (input_mod == NULL) tm_free(mod);
      if (params != NULL) gsl_vector_free(params);
    }
    if (window_coords != NULL) 
      msa_free(msa);
  }

  if (gff != NULL) {
    gff_free_set(gff);
    cm_free(cm);
  }
  if (!quiet) fprintf(stderr, "Done.\n");

  str_free(mod_fname);
  if (output_trees) str_free(out_tree_fname);
  msa_free(source_msa);
  if (input_mod != NULL)
    tm_free(input_mod);
  else
    free_tree(tree);
  if (logf != NULL) fclose(logf);
  free(weight_matrix);
  if (window_coords != NULL) {
    lst_free(window_coords);
    fclose(WINDOWF);
  }
  str_free(tmpstr);
  if (nclasses > 0) {
    lst_free(match_list);
    free(nsites);
    free(gc);
    free(totbases);
  }
  lst_free(cats_to_do);
  if (cats_to_do_str != NULL) 
    lst_free_strings(cats_to_do_str);

  return 0;
}
