/* phyloFit - fit phylogenetic model(s) to a multiple alignment
   
   $Id: phyloFit.c,v 1.17 2004-07-26 18:07:16 acs Exp $
   Written by Adam Siepel, 2002-2004
   Copyright 2002-2004, Adam Siepel, University of California 
*/

#include <stdlib.h>
#include <stdio.h>
#include <lists.h>
#include <stringsplus.h>
#include <msa.h>
#include <gff.h>
#include <category_map.h>
#include <getopt.h>
#include <tree_model.h>
#include <fit_em.h>
#include <subst_mods.h>
#include <string.h>
#include <local_alignment.h>
#include <assert.h>
#include <ctype.h>
#include <tree_likelihoods.h>
#include <numerical_opt.h>
#include <sufficient_stats.h>
#include <maf.h>

/* default minimum number of informative sites (see -I) */
#define DEFAULT_NSITES_THRESHOLD 50

/* default starting alpha for dgamma */
#define DEFAULT_ALPHA 1

void print_usage() {
  printf("\n\
PROGRAM: phyloFit\n\
\n\
DESCRIPTION: \n\
\n\
    Fits one or more tree models to a multiple alignment of DNA\n\
    sequences by maximum likelihood, using the specified tree topology\n\
    and substitution model.  If categories of sites are defined via\n\
    --features and --catmap (see below), then a separate model will be\n\
    estimated for each category.  A description of each model will\n\
    be written to a separate file, with the suffix \".mod\".  These\n\
    .mod files minimally include a substitution rate matrix, a tree with\n\
    branch lengths, and estimates of nucleotide equilibrium\n\
    frequencies.  They may also include information about parameters\n\
    for modeling rate variation.  In addition to the .mod file, a\n\
    description of each each estimated tree will be written to a file with\n\
    the suffix \".nh\" (Newick, a.k.a. New Hampshire, format).\n\
\n\
USAGE: phyloFit [OPTIONS] <msa_fname>\n\
\n\
    <msa_fname> should be a multiple alignment in FASTA format or\n\
    one of several alternative formats (see --msa-format).  For\n\
    backward compatibility, this argument may be preceded by '-m' or\n\
    '--msa'.  Note that --tree is required in most cases.  By default,\n\
    all output files will have the prefix \"phyloFit\" (see\n\
    --out-root).\n\
\n\
EXAMPLES:\n\
\n\
    (If you're like me, you want some basic examples first, and a list\n\
    of all options later.)\n\
\n\
    1. Compute the distance between two aligned sequences (in FASTA file\n\
    pair.fa) under the REV model.\n\
\n\
        phyloFit pair.fa\n\
\n\
    (output is to phyloFit.mod; distance in substitutions per site\n\
    appears in the TREE line in the output file)\n\
\n\
    2. Fit a phylogenetic model to an alignment of human, chimp, mouse,\n\
    and rat sequences.  Use the HKY85 substitution model.  Write output\n\
    to files with prefix \"myfile\".\n\
\n\
        phyloFit --tree \"((human,chimp),(mouse,rat))\" --subst-mod HKY85\n\
            --out-root myfile primate-rodent.fa\n\
\n\
    3. As above, but use the discrete-gamma model for rate variation,\n\
    with 4 rate categories.\n\
\n\
        phyloFit --tree \"((human,chimp),(mouse,rat))\" --subst-mod HKY85\n\
            --out-root myfile --nrates 4 primate-rodent.fa\n\
\n\
    4. As above, but use genome-wide data, stored in the compact\n\
    \"sufficient-statistics\" format (can be produced with \"msa_view\n\
    -o SS\").\n\
\n\
        phyloFit --tree \"((human,chimp),(mouse,rat))\" --subst-mod HKY85\n\
            --out-root myfile --nrates 4 --msa-format SS \n\
            primate-rodent.ss\n\
\n\
    5. Fit a context-dependent phylogenetic model (U2S) to an\n\
    alignment of human, mouse, and rat sequences.  Use\n\
    an EM algorithm for parameter optimization and relax the\n\
    convergence criteria a bit (recommended with context-dependent\n\
    models).  Write a log file for the optimization procedure.\n\
    Consider only non-overlapping pairs of sites.\n\
\n\
        phyloFit --tree \"(human,(mouse,rat))\" --subst-mod U2S --EM\n\
            --precision MED --non-overlapping --log u2s.log --out-root\n\
            hmr-u2s hmr.fa\n\
\n\
    6. As above, but allow overlapping pairs of sites, and compute\n\
    likelihoods by assuming Markov-dependence of columns (see Siepel &\n\
    Haussler, 2004).  The EM algorithm can no longer be used\n\
    (optimization will be much slower).\n\
\n\
        phyloFit --tree \"(human,(mouse,rat))\" --subst-mod U2S\n\
            --precision MED --log u2s-markov.log --markov hmr.fa\n\
\n\
    7. Compute a likelihood using parameter estimates obtained in (5)\n\
    and an assumption of Markov dependence.  This provides a lower\n\
    bound on the likelihood of the Markov-dependent model.\n\
\n\
        phyloFit --init-model hmr-u2s.mod --lnl --markov hmr.fa\n\
\n\
    8. Given an alignment of several mammalian sequences (mammals.fa), a\n\
    tree topology (tree.nh), and a set of gene annotations in GFF\n\
    (genes.gff), fit a separate model to sites in 1st, 2nd, and 3rd\n\
    codon positions.  Use the REV substitution model.  Assume coding\n\
    regions have feature type 'CDS'.\n\
\n\
        phyloFit --tree tree.nh --features genes.gff --out-root mammals-rev\n\
            --catmap \"NCATS = 3; CDS 1-3\" --do-cats 1,2,3 mammals.fa\n\
\n\
    (output will be to mammals-rev.cds-1.mod, mammals-rev.cds-2.mod, and \n\
    mammals-rev.cds-3.mod; first one is for non-coding sites, others are for\n\
    three codon positions)\n\
\n\
\n\
OPTIONS:\n\
\n\
    --tree, -t <tree_fname>|<tree_string>\n\
        (Required if more than three species, or more than two species\n\
        and a non-reversible substitution model, e.g., UNREST, U2, U3)\n\
        Name of file or literal string defining tree topology.  Tree\n\
        must be in Newick format, with the label at each leaf equal to\n\
        the index or name of the corresponding sequence in the alignment\n\
        (indexing begins with 1).  Examples: --tree \"(1,(2,3))\", \n\
        --tree \"(human,(mouse,rat))\".  Currently, the topology must be\n\
        rooted.  When a reversible substitution model is used, the root\n\
        is ignored during the optimization procedure.\n\
\n\
    --subst-mod, -s JC69|F81|HKY85|REV|UNREST|R2|R2S|U2|U2S|R3|R3S|U3|U3S\n\
        (default REV).  Nucleotide substitution model.  JC69, F81, HKY85\n\
        REV, and UNREST have the usual meanings (see, e.g., Yang, \n\
        Goldman, and Friday, 1994).  The others (all considered \"context-\n\
        dependent\") are as defined in Siepel and Haussler, 2004.  The\n\
        options --EM and --precision MED are recommended with context-\n\
        dependent models (see below).\n\
\n\
    --msa-format, -i FASTA|PHYLIP|MPM|MAF|SS\n\
        (default FASTA) Alignment format.  FASTA is as usual.  PHYLIP\n\
        is compatible with the formats used in the PHYLIP and PAML\n\
        packages.  MPM is the format used by the MultiPipMaker aligner\n\
        and some other of Webb Miller's older tools.  MAF (\"Multiple\n\
        Alignment Format\") is used by MULTIZ/TBA and the UCSC Genome\n\
        Browser.  SS is a simple format describing the sufficient\n\
        statistics for phylogenetic inference (distinct columns or\n\
        tuple of columns and their counts).  Note that the program\n\
        \"msa_view\" can be used for file conversion.\n\
\n\
    --out-root, -o <output_fname_root>\n\
        (default \"phyloFit\").  Use specified string as root filename\n\
        for all files created.\n\
\n\
    --min-informative, -I <ninf_sites>\n\
        Require at least <ninf_sites> \"informative\" sites -- i.e., \n\
        sites at which at least two non-gap and non-missing-data ('N'\n\
        or '*') characters are present.  Default is %d.\n\
\n\
    --gaps-as-bases, -G\n\
        Treat alignment gap characters ('-') like ordinary bases.  By\n\
        default, they are treated as missing data.\n\
\n\
    --quiet, -q\n\
        Proceed quietly.\n\
\n\
    --help, -h\n\
        Print this help message.\n\
\n\
\n\
 (Options for controlling and monitoring the optimization procedure)\n\
\n\
    --lnl, -L\n\
        (for use with --init-model) Simply evaluate the log likelihood of\n\
        the specified tree model, without performing any further\n\
        optimization.  Can be used with --post-probs, --expected-subs, and\n\
        --expected-total-subs.\n\
\n\
    --EM, -E \n\
        Fit model(s) using EM rather than the BFGS quasi-Newton\n\
        algorithm (the default).\n\
\n\
    --precision, -p HIGH|MED|LOW\n\
        (default HIGH) Level of precision to use in estimating model\n\
        parameters.  Affects convergence criteria for iterative\n\
        algorithms: higher precision means more iterations and longer\n\
        execution time.\n\
\n\
    --log, -l <log_fname>\n\
        Write log to <log_fname> describing details of the optimization\n\
        procedure.\n\
\n\
    --init-model, -M <mod_fname>\n\
        Initialize with specified tree model.  By choosing good\n\
        starting values for parameters, it is possible to reduce\n\
        execution time dramatically.  If this option is chosen, --tree\n\
        is not required.  Note: currently only one mod_fname may be\n\
        specified; it will be used for all categories.\n\
\n\
    --init-random, -r\n\
        Initialize parameters randomly.  Can be used multiple times to test\n\
        whether the m.l.e. is real.\n\
\n\
    --scale-only, -B\n\
        (for use with --init-model) Estimate only the scale of the tree,\n\
        rather than individual branch lengths (branch proportions fixed).\n\
\n\
    --estimate-freqs, -F\n\
        Estimate equilibrium frequencies by maximum likelihood, rather\n\
        than approximating them by the relative frequencies in the data.\n\
\n\
    --ancestor, -A <seqname>\n\
        Treat specified sequence as the root of the tree.  The tree\n\
        topology must define this sequence to be a child of the root\n\
        (in practice, the branch from the root to the specified\n\
        sequence will be retained, but will be constrained to have\n\
        length zero).\n\
\n\
\n\
 (Options for modeling rate variation)\n\
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
 (Options for separate handling of sites in different annotation categories)\n\
\n\
    --features, -g <fname>\n\
        Annotations file (GFF or BED format) describing features on\n\
        one or more sequences in the alignment.  Together with a\n\
        category map (see --catmap), will be taken to define site\n\
        categories, and a separate model will be estimated for each\n\
        category.  If no category map is specified, a category will be\n\
        assumed for each type of feature, and they will be numbered in\n\
        the order of appearance of the features.  Features are assumed\n\
        to use the coordinate frame of the first sequence in the\n\
        alignment and should be non-overlapping (see 'refeature\n\
        --unique').\n\
\n\
    --catmap, -c <fname>|<string>\n\
        (optionally use with --features) Mapping of feature types to\n\
        category numbers.  Can either give a filename or an \"inline\"\n\
        description of a simple category map, e.g., --catmap \"NCATS =\n\
        3 ; CDS 1-3\" or --catmap \"NCATS = 1 ; UTR 1\".  Note that\n\
        category 0 is reserved for \"background\" (everything that is\n\
        not described by a defined feature type).\n\
\n\
    --do-cats, -C <cat_list>\n\
        (optionally use with --features) Estimate models for only the\n\
        specified categories (comma-delimited list categories, by name\n\
        or numbera).  Default is to fit a model for every category.\n\
\n\
    --reverse-groups, -R <tag>\n\
        (optionally use with --features) Group features by <tag> (e.g.,\n\
        \"transcript_id\" or \"exon_id\") and reverse complement\n\
        segments of the alignment corresponding to groups on the\n\
        reverse strand.  Groups must be non-overlapping (see refeature\n\
        --unique).  Useful with categories corresponding to\n\
        strand-specific phenomena (e.g., codon positions).\n\
\n\
\n\
 (Options for context-dependent substitution models)\n\
\n\
    --markov, -N\n\
        (for use with context-dependent substitutions models and not\n\
        available with --EM.)  Assume Markov dependence of alignment\n\
        columns, and compute the conditional probability of each\n\
        column given its N-1 predecessors using the two-pass algorithm\n\
        described by Siepel and Haussler (2004).  (Here, N is the\n\
        \"order\" of the model, as defined by --subst-mod; e.g., N=1\n\
        for REV, N=2 for U2S, N=3 for U3S.) The alternative (the\n\
        default) is simply to work with joint probabilities of tuples\n\
        of columns.  (You can ensure that these tuples are\n\
        non-overlapping with the --non-overlapping option.)  The use\n\
        of joint probabilities during parameter estimation allows the\n\
        use of the --EM option and can be much faster; in addition, it\n\
        appears to produce nearly equivalent estimates.  If desired,\n\
        parameters can be estimated without --markov, and\n\
        then the likelihood can be evaluated using --lnl and\n\
        --markov together.  This gives a lower bound on the\n\
        likelihood of the Markov-dependent model.\n\
\n\
    --non-overlapping, -V\n\
        (for use with context-dependent substitution models; not\n\
        compatible with --markov, --features, or\n\
        --msa-format SS) Avoid using overlapping tuples of sites\n\
        in parameter estimation.  If a dinucleotide model is selected,\n\
        every other tuple will be considered, and if a nucleotide\n\
        triplet model is selected, every third tuple will be\n\
        considered.  This option cannot be used with an alignment\n\
        represented only by unordered sufficient statistics.\n\
\n\
\n\
 (Options for posterior probabilities)\n\
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
\n\
 (Options for estimation in sliding window)\n\
\n\
    --windows, -w <size,shift>\n\
        Apply a sliding window to the alignment, and fit a separate\n\
        tree to each window.  Arguments specify size of window and\n\
        amount by which to shift it on each iteration, both in bases\n\
        of the first sequence in the alignment (assumed to be the\n\
        reference sequence).  Separate versions of all output files\n\
        will be created for each window.\n\
\n\
    --windows-explicit, -v <window_coord_list>\n\
        Like --windows, except that all start and end coordinates must\n\
        be explicitly specified.  Each successive pair of numbers is\n\
        interpreted as defining the start and end of a window.  Can be\n\
        used with a two-column file and the '*' operator, e.g.,\n\
        --windows-explicit '*mycoords'.\n\
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

void print_window_summary(FILE* WINDOWF, List *window_coords, int win, 
                          int cat, TreeModel *mod, double *gc, double cpg, 
                          int ninf_sites, int nseqs, int header_only) {
  int j;
  if (header_only) {
    fprintf(WINDOWF, "%5s %8s %8s %4s", "win", "beg", "end", "cat");
    fprintf(WINDOWF, " %6s", "GC");
    fprintf(WINDOWF, " %8s", "CpG");
    fprintf(WINDOWF, " %7s", "ninf");
    fprintf(WINDOWF, " %7s\n", "t");
  }
  else {
    fprintf(WINDOWF, "%5d %8d %8d %4d", win/2+1, 
            lst_get_int(window_coords, win), 
            lst_get_int(window_coords, win+1), cat);
    fprintf(WINDOWF, " %6.4f", 
            gsl_vector_get(mod->backgd_freqs, 
                           mod->rate_matrix->inv_states[(int)'G']) + 
            gsl_vector_get(mod->backgd_freqs, 
                           mod->rate_matrix->inv_states[(int)'C']));
    for (j = 0; j < nseqs; j++) fprintf(WINDOWF, " %6.4f", gc[j]);
    fprintf(WINDOWF, " %8.6f", cpg);
    fprintf(WINDOWF, " %7d", ninf_sites);
    fprintf(WINDOWF, " %7.4f\n", tr_total_len(mod->tree));
  }
}

int main(int argc, char *argv[]) {
  char *msa_fname = NULL, *output_fname_root = "phyloFit", 
    *log_fname = NULL, *reverse_group_tag = NULL, *alph = "ACGT";
  int output_trees = TRUE, subst_mod = REV, quiet = FALSE,
    nratecats = 1, use_em = FALSE, window_size = -1, 
    window_shift = -1, use_conditionals = FALSE, 
    precision = OPT_HIGH_PREC, 
    likelihood_only = FALSE, do_bases = FALSE, do_expected_nsubst = FALSE, 
    do_expected_nsubst_tot = FALSE, 
    random_init = FALSE, estimate_backgd = FALSE, estimate_scale_only = FALSE,
    do_column_probs = FALSE, nonoverlapping = FALSE, gaps_as_bases = FALSE;
  unsigned int nsites_threshold = DEFAULT_NSITES_THRESHOLD;
  msa_format_type input_format = FASTA;
  char c;
  FILE *F, *WINDOWF;
  TreeNode *tree = NULL;
  CategoryMap *cm = NULL;
  int i, j, win, opt_idx;
  String *mod_fname, *out_tree_fname, *root_seqname = NULL;
  MSA *msa, *source_msa;
  FILE *logf = NULL;
  String *tmpstr = str_new(STR_SHORT_LEN);
  List *cats_to_do = NULL, *tmplist = NULL, *window_coords = NULL, 
    *cats_to_do_str = NULL;
  double *gc;
  double cpg, alpha = DEFAULT_ALPHA;
  GFF_Set *gff = NULL;
  TreeModel *input_mod = NULL;
  int root_leaf_id = -1;
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
    {"non-overlapping", 0, 0, 'V'},
    {"markov", 0, 0, 'N'},
    {"reverse-groups", 1, 0, 'R'},
    {"init-model", 1, 0, 'M'},
    {"init-random", 0, 0, 'r'},
    {"lnl", 0, 0, 'L'},
    {"scale-only", 0, 0, 'B'},
    {"estimate-freqs", 0, 0, 'F'},
    {"min-informative", 1, 0, 'I'},
    {"gaps-as-bases", 0, 0, 'G'},    
    {"quiet", 0, 0, 'q'},
    {"help", 0, 0, 'h'},
    {"windows", 1, 0, 'w'},
    {"windows-explicit", 1, 0, 'v'},
    {"ancestor", 1, 0, 'A'},
    {"post-probs", 1, 0, 'P'},
    {"expected-subs", 1, 0, 'X'},
    {"expected-total-subs", 1, 0, 'Z'},
    {"column-probs", 1, 0, 'U'},
    {"rate-constants", 1, 0, 'K'},
  };

  while ((c = getopt_long(argc, argv, "m:t:s:g:c:C:i:o:k:a:l:w:v:M:p:A:I:K:GVEeNDRTqLPXZUBFrh", long_opts, &opt_idx)) != -1) {
    switch(c) {
    case 'm':
      msa_fname = optarg;
      break;
    case 't':
      if (optarg[0] == '(')     /* in this case, assume topology given
                                   at command line */
        tree = tr_new_from_string(optarg);
      else 
        tree = tr_new_from_file(fopen_fname(optarg, "r"));
      break;
    case 's':
      subst_mod = tm_get_subst_mod_type(optarg);
      if (subst_mod == UNDEF_MOD) 
        die("ERROR: illegal substitution model.  Type \"phyloFit -h\" for usage.\n");
      break;
    case 'g':
      gff = gff_read_set(fopen_fname(optarg, "r"));
      break;
    case 'c':
      cm = cm_new_string_or_file(optarg);
      break;
    case 'C':
      cats_to_do_str = get_arg_list(optarg);
      break;
    case 'V':
      nonoverlapping = TRUE;
      break;
    case 'o':
      output_fname_root = optarg;
      break;
    case 'T':                   /* this is now done by default, but
                                   we'll leave the option in for
                                   backward compatibility */
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
      reverse_group_tag = optarg;
      break;
    case 'i':
      input_format = msa_str_to_format(optarg);
      if (input_format == -1)
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
      if (lst_size(tmplist) % 2 != 0) 
        die("ERROR: argument to --windows-explicit must be a list of even length.\n");
      window_coords = str_list_as_int(tmplist);
      lst_free(tmplist);
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
      input_mod = tm_new_from_file(fopen_fname(optarg, "r"));
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
    case 'G':
      gaps_as_bases = TRUE;
      alph = "ACGT-";
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

  if (msa_fname == NULL) {
    if (optind >= argc) 
      die("ERROR: missing alignment filename.  Type 'phyloFit -h' for usage.\n");
    msa_fname = argv[optind];
  }

  if (use_conditionals && use_em) 
    die("ERROR: Cannot use --markov with --EM.  Type 'phyloFit -h' for usage.\n");

  if (likelihood_only && input_mod == NULL) 
    die("ERROR: --lnl requires --init-model.  Type 'phyloFit -h' for usage.\n");

  if (nonoverlapping && (use_conditionals || gff != NULL || 
                         cats_to_do_str || input_format == SS))
    die("ERROR: cannot use --non-overlapping with --markov, --features,\n--msa-format SS, or --do-cats.\n");

  if (gaps_as_bases && subst_mod != JC69 && subst_mod != F81 && 
      subst_mod != REV && subst_mod != UNREST)
    die("ERROR: --gaps-as-bases currently only supported with JC69, F81, REV, and UNREST.\n");
                                /* with HKY, yields undiagonalizable matrix */
  
  if (gff != NULL && cm == NULL) cm = cm_new_from_features(gff);

  /* internally, --non-overlapping is accomplished via --do-cats */
  if (nonoverlapping) {
    cats_to_do_str = lst_new_ptr(1);
    lst_push_ptr(cats_to_do_str, str_new_charstr("1"));
  }

  /* read alignment */
  if (!quiet) fprintf(stderr, "Reading alignment from %s ...\n", msa_fname);
  if (input_format == MAF)
    msa = maf_read(fopen_fname(msa_fname, "r"), NULL, tm_order(subst_mod) + 1, 
                   gff, cm, nonoverlapping ? tm_order(subst_mod) + 1 : -1, 
                   FALSE, reverse_group_tag, NO_STRIP, FALSE);
  else 
    msa = msa_new_from_file(fopen_fname(msa_fname, "r"), input_format, alph);

  if (tree == NULL) {
    if (msa->nseqs == 2)
      tree = tr_new_from_string("(1,2)");
    else if (msa->nseqs == 3 && tm_is_reversible(subst_mod))
      tree = tr_new_from_string("(1,(2,3))");
    else die("ERROR: --tree required.\n");
  }
  else 
    tr_number_leaves(tree, msa->names, msa->nseqs);
                         /* convert from names to numbers, if
                            necessary */

  /* make sure alignment and tree topology are consistent */
  if (msa->nseqs * 2 - 1 != 
      (input_mod == NULL ? tree->nnodes : input_mod->tree->nnodes)) 
    die("ERROR: Tree must have 2n-1 nodes, where n is the number of sequences in the\nalignment.  Even with a reversible model, specify a rooted tree; the root\nwill be ignored in the optimization procedure.\n");

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
  
  msa_remove_N_from_alph(msa);  /* for backward compatibility */

  /* set up for categories */
  /* first label sites, if necessary */
  if (gff != NULL && input_format != MAF) {
    /* convert GFF to coordinate frame of alignment */
    msa_map_gff_coords(msa, gff, 1, 0, 0, NULL);

    /* reverse complement segments of MSA corresponding to features on
       reverse strand (if necessary) */
    if (reverse_group_tag != NULL) {
      gff_group(gff, reverse_group_tag);
      msa_reverse_compl_feats(msa, gff, NULL);
    }

    /* label categories */
    if (!quiet) fprintf(stderr, "Labeling alignment sites by category ...\n");
    msa_label_categories(msa, gff, cm);
  }
  else if (nonoverlapping && input_format != MAF) {
                                /* (already taken care of if MAF) */
    int cycle_size = tm_order(subst_mod) + 1;
    assert(msa->seqs != NULL && msa->ss == NULL);  /* need explicit seqs */
    msa->categories = smalloc(msa->length * sizeof(int));
    for (i = 0; i < msa->length; i++) 
      msa->categories[i] = (i % cycle_size) + 1;
    msa->ncats = cycle_size;
  }
  /* at this point, we have msa->ncats > 0 iff we intend to do
     category-by-category estimation */

  /* now set up list of categories to process.  There are several
     cases to consider */
  if (msa->ncats < 0) {         
    if (cats_to_do_str != NULL)
      fprintf(stderr, "WARNING: ignoring --do-cats; no category information.\n");
    cats_to_do = lst_new_int(1);
    lst_push_int(cats_to_do, -1);
                                /* no categories -- pool all sites */
  }
  else if (cats_to_do_str == NULL) {
    cats_to_do = lst_new_int(msa->ncats + 1);
    for (i = 0; i <= msa->ncats; i++) lst_push_int(cats_to_do, i);
                                /* have categories but no --do-cats --
                                   process all categories */
  }
  else if (cm != NULL) 
    cats_to_do = cm_get_category_list(cm, cats_to_do_str, 0);
                                /* have --do-cats and category map;
                                   use cm_get_category_list (allows
                                   use of names as well as numbers) */
  else if (cats_to_do_str != NULL)
    cats_to_do = str_list_as_int(cats_to_do_str);
                                /* have --do-cats but no category map;
                                   use literal numbers */

  /* set up windows, if necessary */
  if (window_size != -1) {
    if (window_coords != NULL) 
      die("ERROR: cannot use both --windows and --windows-explicit.\n");
    window_coords = lst_new_int(msa->length/window_shift + 1);
    for (i = 1; i < msa->length; i += window_shift) {
      lst_push_int(window_coords, i);
      lst_push_int(window_coords, 
                   min(i + window_size - 1, msa->length));
    }
  }
  if (window_coords != NULL) {
    /* set up summary file */
    String *sumfname = str_new_charstr(output_fname_root);
    msa_coord_map *map;

    str_append_charstr(sumfname, ".win-sum");
    WINDOWF = fopen_fname(sumfname->chars, "w+");
    print_window_summary(WINDOWF, NULL, 0, 0, NULL, NULL, 0, 0, 0, TRUE);
    
    /* map to coord frame of alignment */
    map = msa_build_coord_map(msa, 1);
    for (i = 0; i < lst_size(window_coords); i += 2) {
      lst_set_int(window_coords, i, 
                  msa_map_seq_to_msa(map, lst_get_int(window_coords, i)));
      lst_set_int(window_coords, i+1, 
                  msa_map_seq_to_msa(map, lst_get_int(window_coords, i+1)));
    }
    msa_map_free(map);
  }

  /* now estimate models (window by window, if necessary) */
  mod_fname = str_new(STR_MED_LEN);
  if (output_trees) out_tree_fname = str_new(STR_MED_LEN);
  source_msa = msa;
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
      unsigned int ninf_sites;

      if (input_mod == NULL) 
        mod = tm_new(tr_create_copy(tree), NULL, NULL, subst_mod, 
                     msa->alphabet, nratecats, alpha, rate_consts, 
                     root_leaf_id);
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
        tm_reinit(mod, subst_mod, nratecats, newalpha, rate_consts, NULL);
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
      if (ninf_sites < nsites_threshold) {
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

        if (!quiet) {
          fprintf(stderr, "Fitting tree model to %s using %s%s ...\n",
                  tmpstr->chars, tm_get_subst_mod_string(subst_mod),
                  mod->nratecats > 1 ? " (with rate variation)" : "");
          if (log_fname != NULL)
            fprintf(stderr, "(writing log to %s)\n", log_fname);
        }

        if (msa->ss == NULL) 
          /* ensures that sufficient stats are obtained only for
             categories of interest */
          ss_from_msas(msa, mod->order+1, 0, 
                       cats_to_do_str != NULL ? cats_to_do : NULL, 
                       NULL, NULL, -1);

        if (i == 0) ss_collapse_missing(msa, !gaps_as_bases);
                                /* reduce number of tuples as much as
                                   possible */

        if (use_em)
          tm_fit_em(mod, msa, params, cat, precision, logf);
        else
          tm_fit(mod, msa, params, cat, precision, logf);
      }

      str_cpy_charstr(mod_fname, output_fname_root);
      if (window_coords != NULL) {
        str_append_charstr(mod_fname, ".win-");
        str_append_int(mod_fname, win/2 + 1);
      }
      if (cat != -1 && nonoverlapping == FALSE) {
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
      if (output_trees) {
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
        tr_print(F, trcpy, 1);
        fclose(F);
        tr_free(trcpy);
      }

      /* print window summary, if window mode */
      if (window_coords != NULL) 
        print_window_summary(WINDOWF, window_coords, win, cat, mod, gc, 
                             cpg, ninf_sites, msa->nseqs, FALSE);

      if (input_mod == NULL) tm_free(mod);
      if (params != NULL) gsl_vector_free(params);
    }
    if (window_coords != NULL) 
      msa_free(msa);
  }

  if (!quiet) fprintf(stderr, "Done.\n");

  return 0;
}
