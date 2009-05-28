/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2009 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* phyloFit - fit phylogenetic model(s) to a multiple alignment
   
   $Id: phyloFit.c,v 1.38 2008-11-12 02:07:58 acs Exp $ */

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
#include "phyloFit.help"

/* default minimum number of informative sites (see -I) */
#define DEFAULT_NSITES_THRESHOLD 50

/* default starting alpha for dgamma */
#define DEFAULT_ALPHA 1

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
   ".expsub", and ".exptotsub", respectively).    */
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

  /* FIXME: rate variation!     need rate post probs! */
  assert(mod->nratecats == 1);

  /* compute desired stats */
  assert(mod->tree_posteriors == NULL); 
  if (!quiet) 
    fprintf(stderr, "Computing posterior probabilities and/or related stats ...\n");
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
    fprintf(POSTPROBF, "tuple     ");
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
      fprintf(POSTPROBF, "%-6d %5s      ", tup, coltupstr);
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

    fprintf(EXPSUBF, "%-3s %10s %7s ", "#", "tuple", "count");
    for (node = 0; node < mod->tree->nnodes; node++) {
      n = lst_get_ptr(tr_postorder(mod->tree), node);
      if (n == mod->tree) continue;
      fprintf(EXPSUBF, " node_%-2d", n->id);
    }
    fprintf(EXPSUBF, "    total\n");
    for (tup = 0; tup < msa->ss->ntuples; tup++) {
      double total = 0;

      if ((cat >= 0 && msa->ss->cat_counts[cat][tup] == 0) ||
          msa->ss->counts[tup] == 0) continue;

      tuple_to_string_pretty(coltupstr, msa, tup);
      fprintf(EXPSUBF, "%-3d %10s %.0f ", tup, coltupstr, msa->ss->counts[tup]);
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
branch of the tree.     Nodes of the tree are visited in a postorder traversal,\n\
and each node is taken to be representative of the branch between itself and\n\
its parent.     Starting bases or tuples of bases appear on the vertical axis\n\
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
            vec_get(mod->backgd_freqs, 
                           mod->rate_matrix->inv_states[(int)'G']) + 
            vec_get(mod->backgd_freqs, 
                           mod->rate_matrix->inv_states[(int)'C']));
    for (j = 0; j < nseqs; j++) fprintf(WINDOWF, " %6.4f", gc[j]);
    fprintf(WINDOWF, " %8.6f", cpg);
    fprintf(WINDOWF, " %7d", ninf_sites);
    fprintf(WINDOWF, " %7.4f\n", tr_total_len(mod->tree));
  }
}

int main(int argc, char *argv[]) {
  char *msa_fname = NULL, *output_fname_root = "phyloFit", 
    *log_fname = NULL, *reverse_group_tag = NULL, *alph = "ACGT", 
    *root_seqname = NULL, *subtree_name = NULL;
  int subst_mod = REV, quiet = FALSE, nratecats = 1, use_em = FALSE, 
    window_size = -1, window_shift = -1, use_conditionals = FALSE, 
    precision = OPT_HIGH_PREC, likelihood_only = FALSE, do_bases = FALSE, 
    do_expected_nsubst = FALSE, do_expected_nsubst_tot = FALSE, 
    random_init = FALSE, estimate_backgd = FALSE, estimate_scale_only = FALSE,
    do_column_probs = FALSE, nonoverlapping = FALSE, gaps_as_bases = FALSE,
    no_freqs = FALSE, no_rates = FALSE;
  unsigned int nsites_threshold = DEFAULT_NSITES_THRESHOLD;
  msa_format_type input_format = FASTA;
  char c;
  FILE *F, *WINDOWF;
  TreeNode *tree = NULL;
  CategoryMap *cm = NULL;
  int i, j, win, opt_idx;
  String *mod_fname;
  MSA *msa, *source_msa;
  FILE *logf = NULL;
  String *tmpstr = str_new(STR_SHORT_LEN);
  List *cats_to_do = NULL, *tmplist = NULL, *window_coords = NULL, 
    *cats_to_do_str = NULL, *ignore_branches = NULL;
  double *gc;
  double cpg, alpha = DEFAULT_ALPHA;
  GFF_Set *gff = NULL;
  TreeModel *input_mod = NULL;
  int root_leaf_id = -1;
  List *rate_consts = NULL;
  char tmpchstr[STR_MED_LEN];

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
    {"scale-subtree", 1, 0, 'S'},
    {"estimate-freqs", 0, 0, 'F'},
    {"no-freqs", 0, 0, 'f'},
    {"no-rates", 0, 0, 'n'},
    {"min-informative", 1, 0, 'I'},
    {"gaps-as-bases", 0, 0, 'G'},     
    {"quiet", 0, 0, 'q'},
    {"help", 0, 0, 'h'},
    {"windows", 1, 0, 'w'},
    {"windows-explicit", 1, 0, 'v'},
    {"ancestor", 1, 0, 'A'},
    {"post-probs", 0, 0, 'P'},
    {"expected-subs", 0, 0, 'X'},
    {"expected-total-subs", 0, 0, 'Z'},
    {"column-probs", 0, 0, 'U'},
    {"rate-constants", 1, 0, 'K'},
    {"ignore-branches", 1, 0, 'b'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "m:t:s:g:c:C:i:o:k:a:l:w:v:M:p:A:I:K:S:b:GVEeNDRTqLPXZUBFfnrh", long_opts, &opt_idx)) != -1) {
    switch(c) {
    case 'm':
      msa_fname = optarg;
      break;
    case 't':
      if (optarg[0] == '(')        /* in this case, assume topology given
                                   at command line */
        tree = tr_new_from_string(optarg);
      else 
        tree = tr_new_from_file(fopen_fname(optarg, "r"));
      break;
    case 's':
      subst_mod = tm_get_subst_mod_type(optarg);
      if (subst_mod == UNDEF_MOD) 
        die("ERROR: illegal substitution model.     Type \"phyloFit -h\" for usage.\n");
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
    case 'T': 
      fprintf(stderr, "WARNING: --output-tree (-T) is deprecated; leaf names now appear in .mod\nfile.    (If necessary, use 'tree_doctor --tree-only' to extract tree.)\n");
      break;
    case 'k':
      nratecats = get_arg_int_bounds(optarg, 0, INFTY);
      break;
    case 'a':
      alpha = get_arg_dbl(optarg);
      break;
    case 'R':
      reverse_group_tag = optarg;
      break;
    case 'i':
      input_format = msa_str_to_format(optarg);
      if (input_format == -1)
        die("ERROR: unrecognized alignment format.    Type 'phyloFit -h' for usage.\n");
      break;
    case 'l':
      log_fname = optarg;
      if (!strcmp(log_fname, "-"))
        logf = stderr;
      else
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
      use_em = TRUE;
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
      random_init = TRUE;
      break;
    case 'L':
      likelihood_only = TRUE;
      break;
    case 'q':
      quiet = TRUE;
      break;
    case 'P':
      do_bases = TRUE;
      break;
    case 'X':
      do_expected_nsubst = TRUE;
      break;
    case 'Z':
      do_expected_nsubst_tot = TRUE;
      break;
    case 'U':
      likelihood_only = TRUE;        /* force -L */
      nsites_threshold = 0;        /* also force this; typical use is
                                   with small number of tuples, no
                                   tuple_idx */
      do_column_probs = TRUE;
      break;
    case 'A':
      root_seqname = optarg;
      break;
    case 'I':
      nsites_threshold = get_arg_int(optarg);
      break;
    case 'G':
      gaps_as_bases = TRUE;
      alph = "ACGT-";
      break;
    case 'B':
      estimate_scale_only = TRUE;
      break;
    case 'S':
      subtree_name = optarg;
      break;       
    case 'F':
      estimate_backgd = TRUE;
      break;
    case 'f':
      no_freqs = TRUE;
      break;
    case 'n':
      no_rates = TRUE;
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
    case 'b':
      ignore_branches = get_arg_list(optarg);
      break;
    case 'h':
      printf("%s", HELP);
      exit(0);
    case '?':
      die("ERROR: illegal argument.     Type 'phyloFit -h' for usage.\n");
    }
  }

  if (msa_fname == NULL) {
    if (optind >= argc) 
      die("ERROR: missing alignment filename.  Type 'phyloFit -h' for usage.\n");
    msa_fname = argv[optind];
  }

  if (use_conditionals && use_em) 
    die("ERROR: Cannot use --markov with --EM.    Type 'phyloFit -h' for usage.\n");

  if (likelihood_only && input_mod == NULL) 
    die("ERROR: --lnl requires --init-model.  Type 'phyloFit -h' for usage.\n");

  if (input_mod != NULL && tree != NULL)
    die("ERROR: --tree is not allowed with --init-model.\n");

  if (nonoverlapping && (use_conditionals || gff != NULL || 
                         cats_to_do_str || input_format == SS))
    die("ERROR: cannot use --non-overlapping with --markov, --features,\n--msa-format SS, or --do-cats.\n");

  if (gaps_as_bases && subst_mod != JC69 && subst_mod != F81 && 
      subst_mod != HKY85G && subst_mod != REV && subst_mod != UNREST &&
      subst_mod != SSREV)
    die("ERROR: --gaps-as-bases currently only supported with JC69, F81, HKY85+Gap, REV, SSREV, and UNREST.\n");
                                /* with HKY, yields undiagonalizable matrix */

  if ((no_freqs || no_rates) && input_mod == NULL)
    die("ERROR: --init-model required with --no-freqs and/or --no-rates.\n");

  if (no_freqs && estimate_backgd)
    die("ERROR: can't use both --no-freqs and --estimate-freqs.\n");
  
  if (gff != NULL && cm == NULL) cm = cm_new_from_features(gff);

  /* internally, --non-overlapping is accomplished via --do-cats */
  if (nonoverlapping) {
    cats_to_do_str = lst_new_ptr(1);
    lst_push_ptr(cats_to_do_str, str_new_charstr("1"));
  }

  /* read alignment */
  if (!quiet) fprintf(stderr, "Reading alignment from %s ...\n", msa_fname);
  if (input_format == MAF) {
    msa = maf_read(fopen_fname(msa_fname, "r"), NULL, tm_order(subst_mod) + 1, 
                   NULL, gff, cm, 
                   nonoverlapping ? tm_order(subst_mod) + 1 : -1, 
                   FALSE, reverse_group_tag, NO_STRIP, FALSE);
    if (gaps_as_bases) msa_reset_alphabet(msa, alph);
  }
  else 
    msa = msa_new_from_file(fopen_fname(msa_fname, "r"), input_format, alph);

  if (tree == NULL) {
    if (input_mod != NULL) tree = input_mod->tree;
    else if (msa->nseqs == 2) {
      sprintf(tmpchstr, "(%s,%s)", msa->names[0], msa->names[1]);
      tree = tr_new_from_string(tmpchstr);
    }
    else if (msa->nseqs == 3 && tm_is_reversible(subst_mod)) {
      sprintf(tmpchstr, "(%s,(%s,%s))", msa->names[0], msa->names[1], 
              msa->names[2]);
      tree = tr_new_from_string(tmpchstr);
    }
    else die("ERROR: --tree required.\n");
  }

  /* allow for specified ancestor */
  if (root_seqname != NULL) {
    TreeNode *rl;
    if (tree == NULL || tm_is_reversible(subst_mod)) 
      die("ERROR: --ancestor requires --tree and a non-reversible model.\n");
    rl = tr_get_node(tree, root_seqname);     
    if (rl == NULL || rl->parent != tree) 
      die("ERROR: Sequence specified by --ancestor must be a child of the root.\n");
    root_leaf_id = rl->id;
  }
  
  if (msa_alph_has_lowercase(msa)) msa_toupper(msa); 
  msa_remove_N_from_alph(msa);    /* for backward compatibility */

  /* set up for categories */
  /* first label sites, if necessary */
  if (gff != NULL && input_format != MAF) {
    if (msa->idx_offset > 0) {
      /* if these's an offset, we'll just subtract it from all features */
      for (i = 0; i < lst_size(gff->features); i++) {
        GFF_Feature *f = lst_get_ptr(gff->features, i);
        f->start -= msa->idx_offset;
        f->end -= msa->idx_offset;
      }
      msa->idx_offset = 0;
    }

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

  /* now set up list of categories to process.    There are several
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
      Vector *params = NULL;
      List *pruned_names;
      int old_nnodes, cat = lst_get_int(cats_to_do, i);
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

      if (estimate_scale_only || estimate_backgd || no_rates) {
        tm_free_rmp(mod);

        if (estimate_scale_only) {
          mod->estimate_branchlens = TM_SCALE_ONLY;

          if (subtree_name != NULL) { /* estimation of subtree scale */
            String *s1 = str_new_charstr(subtree_name), 
              *s2 = str_new_charstr(subtree_name);
            str_root(s1, ':'); str_suffix(s2, ':'); /* parse string */
            mod->subtree_root = tr_get_node(mod->tree, s1->chars);
            if (mod->subtree_root == NULL)
              die("ERROR: no node named '%s'.\n", s1->chars);
            if (s2->length > 0) {
              if (str_equals_charstr(s2, "loss")) mod->scale_sub_bound = LB;
              else if (str_equals_charstr(s2, "gain")) mod->scale_sub_bound = UB;
              else die("ERROR: unrecognized suffix '%s'\n", s2->chars);
            }
            str_free(s1); str_free(s2);
          }
        }
        
        if (no_rates)
          mod->estimate_ratemat = FALSE;

        mod->estimate_backgd = estimate_backgd;
        tm_init_rmp(mod);        /* necessary because number of
                                    parameters changes */
      }

      if (ignore_branches != NULL) 
        tm_set_ignore_branches(mod, ignore_branches);

      old_nnodes = mod->tree->nnodes;
      pruned_names = lst_new_ptr(msa->nseqs);
      tm_prune(mod, msa, pruned_names);
      if (lst_size(pruned_names) == (old_nnodes + 1) / 2)
        die("ERROR: no match for leaves of tree in alignment (leaf names must match alignment names).\n");
      if (!quiet && lst_size(pruned_names) > 0) {
        fprintf(stderr, "WARNING: pruned away leaves of tree with no match in alignment (");
        for (j = 0; j < lst_size(pruned_names); j++)
          fprintf(stderr, "%s%s", ((String*)lst_get_ptr(pruned_names, j))->chars, 
                  j < lst_size(pruned_names) - 1 ? ", " : ").\n");
      }
      lst_free_strings(pruned_names);
      lst_free(pruned_names);

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
          params = tm_params_new_init_from_model(input_mod);
        else
          params = tm_params_init(mod, .1, 5, alpha);     
        if (input_mod != NULL && mod->backgd_freqs != NULL && !no_freqs) {
          /* in some cases, the eq freqs are needed for
             initialization, but now they should be re-estimated --
             UNLESS user specifies --no-freqs */
          vec_free(mod->backgd_freqs);
          mod->backgd_freqs = NULL;
        }

        if (msa->ss == NULL) {    /* get sufficient stats if necessary */
          if (!quiet)
            fprintf(stderr, "Extracting sufficient statistics ...\n");
          ss_from_msas(msa, mod->order+1, 0, 
                       cats_to_do_str != NULL ? cats_to_do : NULL, 
                       NULL, NULL, -1);
          /* (sufficient stats obtained only for categories of interest) */
      
          if (msa->length > 1000000) { /* throw out original data if
                                          very large */
            for (j = 0; j < msa->nseqs; j++) free(msa->seqs[j]);
            free(msa->seqs);
            msa->seqs = NULL;
          }
        }

        if (i == 0) {
          if (!quiet) fprintf(stderr, "Compacting sufficient statistics ...\n");
          ss_collapse_missing(msa, !gaps_as_bases);
                                /* reduce number of tuples as much as
                                   possible */
        }

        if (!quiet) {
          fprintf(stderr, "Fitting tree model to %s using %s%s ...\n",
                  tmpstr->chars, tm_get_subst_mod_string(subst_mod),
                  mod->nratecats > 1 ? " (with rate variation)" : "");
          if (log_fname != NULL)
            fprintf(stderr, "(writing log to %s)\n", log_fname);
        }

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

      /* print window summary, if window mode */
      if (window_coords != NULL) 
        print_window_summary(WINDOWF, window_coords, win, cat, mod, gc, 
                             cpg, ninf_sites, msa->nseqs, FALSE);

      if (input_mod == NULL) tm_free(mod);
      if (params != NULL) vec_free(params);
    }
    if (window_coords != NULL) 
      msa_free(msa);
  }
  
  if (!quiet) fprintf(stderr, "Done.\n");
  
  return 0;
}
