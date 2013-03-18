/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <misc.h>
#include <tree_likelihoods.h>
#include <sufficient_stats.h>
#include <maf.h>
#include <pbs_code.h>
#include "prequel.help"

void do_indels(MSA *msa, TreeModel *mod);

int main(int argc, char *argv[]) {
  char c;
  int opt_idx, node;
  FILE *out_f = NULL, *msa_f, *mod_f;
  char *out_root;
  TreeModel *mod;
  MSA *msa;
  char out_fname[STR_MED_LEN];

  struct option long_opts[] = {
    {"refseq", 1, 0, 'r'},
    {"msa-format", 1, 0, 'i'},
    {"seqs", 1, 0, 's'},
    {"exclude", 0, 0, 'x'},
    {"no-probs", 0, 0, 'n'},
    {"suff-stats", 0, 0, 'S'},
    {"encode", 1, 0, 'e'},
    {"keep-gaps", 0, 0, 'k'},
    {"gibbs", 1, 0, 'G'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  /* arguments and defaults for options */
  FILE *refseq_f = NULL;
  msa_format_type msa_format = UNKNOWN_FORMAT;
  int suff_stats = FALSE, exclude = FALSE, keep_gaps = FALSE, do_probs = TRUE;
  List *seqlist = NULL;
  PbsCode *code = NULL;
  int gibbs_nsamples = -1;

  while ((c = getopt_long(argc, argv, "r:i:s:e:knxSh", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'r':
      refseq_f = phast_fopen(optarg, "r");
      break;
    case 'i':
      msa_format = msa_str_to_format(optarg);
      if (msa_format == UNKNOWN_FORMAT)
	die("ERROR: unrecognized alignment format.\n");
      break;
    case 'S':
      suff_stats = TRUE;
      break;
    case 'e':
      code = pbs_new_from_file(phast_fopen(optarg, "r"));
      break;
    case 's':
      seqlist = get_arg_list(optarg);
      break;
    case 'x':
      exclude = TRUE;
      break;
    case 'n':
      do_probs = FALSE;
      break;
    case 'k':
      keep_gaps = TRUE;
      break;
    case 'G':
      gibbs_nsamples = get_arg_int_bounds(optarg, 1, INFTY);
      break;
    case 'h':
      printf("%s", HELP);
      exit(0);
    case '?':
      die("Bad argument.  Try 'prequel -h'.\n");
    }
  }

  if (optind != argc - 3)
    die("Three arguments required.  Try 'prequel -h'.\n");

  set_seed(-1);

  if (!do_probs && (suff_stats || code != NULL))
    die("ERROR: --no-probs can't be used with --suff-stats or --encode.\n");

  msa_f = phast_fopen(argv[optind], "r");
  if (msa_format == UNKNOWN_FORMAT)
    msa_format = msa_format_for_content(msa_f, 1);
  fprintf(stderr, "Reading alignment from %s...\n", argv[optind]);
  if (msa_format == MAF) {
    msa = maf_read(msa_f, refseq_f, 1, NULL, NULL, NULL, -1, !suff_stats, NULL,
                   NO_STRIP, FALSE);
    /* (no need to store order if suff_stats mode) */
  }
  else 
    msa = msa_new_from_file_define_format(msa_f, msa_format, NULL);

  if (msa->ss == NULL) {
    fprintf(stderr, "Extracting sufficient statistics...\n");
    ss_from_msas(msa, 1, TRUE, NULL, NULL, NULL, -1, 0);
  }
  else if (msa->ss->tuple_idx == NULL && !suff_stats)
    die("ERROR: ordered representation of alignment required unless --suff-stats.\n");

  mod_f = phast_fopen(argv[optind+1], "r");
  out_root = argv[optind+2];

  mod = tm_new_from_file(mod_f, 1);

  /* MH prune just like in phastcons */
  int old_nnodes = mod->tree->nnodes;
  List *pruned_names = lst_new_ptr(msa->nseqs);
  tm_prune(mod, msa, pruned_names);  
  if (lst_size(pruned_names) == (old_nnodes + 1) / 2)
    die("ERROR: no match for leaves of tree in alignment (leaf names must match alignment names).\n");
  if (lst_size(pruned_names) > 0) {
    fprintf(stderr, "WARNING: pruned away leaves of tree with no match in alignment (");
    int j;
    for (j = 0; j < lst_size(pruned_names); j++)
      fprintf(stderr, "%s%s", ((String*)lst_get_ptr(pruned_names, j))->chars, j < lst_size(pruned_names) - 1 ? ", " : ").\n");
  }
  lst_free_strings(pruned_names);  


  tr_name_ancestors(mod->tree);

  if (mod->order != 0)
    die("ERROR: Only single nucleotide models are supported.\n");

  if (mod->nratecats > 1)
    die("ERROR: Rate variation not supported.\n");


  mod->tree_posteriors = tl_new_tree_posteriors(mod, msa, TRUE, FALSE, 
                                                FALSE, FALSE, FALSE, FALSE,
						FALSE);

  fprintf(stderr, "Computing posterior probabilities...\n");

  if (gibbs_nsamples > 0)
    die("ERROR: --gibbs not implemented yet.");
  /*     gb_sample_ancestral_seqs(mod, msa, mod->tree_posteriors, gibbs_nsamples); */
  else
    tl_compute_log_likelihood(mod, msa, NULL, NULL, -1, mod->tree_posteriors);

  fprintf(stderr, "Reconstructing indels by parsimony...\n");
  do_indels(msa, mod);

  for (node = 0; node < mod->tree->nnodes; node++) {
    int i, j;
    TreeNode *n = lst_get_ptr(mod->tree->nodes, node);

    if (n->lchild == NULL || n->rchild == NULL) continue;

    if (seqlist != NULL) {
      int in_list = str_in_list_charstr(n->name, seqlist);
      if ((in_list && exclude) || (!in_list && !exclude))
        continue;
    }

    fprintf(stderr, "Writing output for ancestral node '%s'...\n", 
            n->name);

    if (suff_stats) {
      if (out_f == NULL) {
        sprintf(out_fname, "%s.stats", out_root);
        out_f = phast_fopen(out_fname, "w+");

        fprintf(out_f, "#count\t");
        for (j = 0; j < mod->rate_matrix->size; j++) 
          fprintf(out_f, "p(%c)%c", mod->rate_matrix->states[j], 
                  j == mod->rate_matrix->size - 1 ? '\n' : '\t');
      }

      for (i = 0; i < msa->ss->ntuples; i++) {
        if (mod->tree_posteriors->base_probs[0][0][node][i] == -1)
          continue;		/* no base this node */
        fprintf(out_f, "%.0f\t", msa->ss->counts[i]);
        for (j = 0; j < mod->rate_matrix->size; j++) {
          fprintf(out_f, "%f%c", 
                  mod->tree_posteriors->base_probs[0][j][node][i], 
                  j == mod->rate_matrix->size - 1 ? '\n' : '\t');
        }
      }
    }

    else if (code == NULL && do_probs) {	/* ordinary sequence-by-sequence 
						   output */
      sprintf(out_fname, "%s.%s.probs", out_root, n->name);
      out_f = phast_fopen(out_fname, "w+");

      fprintf(out_f, "#");
      for (j = 0; j < mod->rate_matrix->size; j++) 
        fprintf(out_f, "p(%c)%c", mod->rate_matrix->states[j], 
                j == mod->rate_matrix->size - 1 ? '\n' : '\t');

      for (i = 0; i < msa->length; i++) {
        if (mod->tree_posteriors->base_probs[0][0][node][msa->ss->tuple_idx[i]] == -1) {
          /* no base */
          if (keep_gaps) fprintf(out_f, "-\n"); 
          /* otherwise do nothing */
        }
        else 
          for (j = 0; j < mod->rate_matrix->size; j++) 
            fprintf(out_f, "%f%c", 
                    mod->tree_posteriors->base_probs[0][j][node][msa->ss->tuple_idx[i]], 
                    j == mod->rate_matrix->size - 1 ? '\n' : '\t');
      }

      phast_fclose(out_f);
    }

    else if (code == NULL && !do_probs) {	/* write point estimates
						   to FASTA file */
      char *outseq = smalloc((msa->length + 1) * sizeof(char));
      int len = 0;

      for (i = 0; i < msa->length; i++) {
        if (mod->tree_posteriors->base_probs[0][0][node][msa->ss->tuple_idx[i]] == -1) {
          /* no base */
          if (keep_gaps) outseq[len++] = GAP_CHAR;
          /* otherwise do nothing */
        }
        else {
          double maxprob = 0;
          int maxidx = -1;
          for (j = 0; j < mod->rate_matrix->size; j++) {
            if (mod->tree_posteriors->base_probs[0][j][node][msa->ss->tuple_idx[i]] > maxprob) {
              maxprob = mod->tree_posteriors->base_probs[0][j][node][msa->ss->tuple_idx[i]];
              maxidx = j;
            }
          }
          outseq[len++] = mod->rate_matrix->states[maxidx];
        }
      }
      outseq[len] = '\0';

      /* print in FASTA format */
      sprintf(out_fname, "%s.%s.fa", out_root, n->name);
      out_f = phast_fopen(out_fname, "w+");
      print_seq_fasta(out_f, outseq, n->name, len);
      phast_fclose(out_f);
      sfree(outseq); 
    }

    else {			/* encoded sequence-by-sequence
				   output */
      double error, tot_error = 0;
      int ngaps = 0;
      Vector *v;
      unsigned *encoded;

      /* first encode tuple by tuple */
      v = vec_new(mod->rate_matrix->size);
      encoded = smalloc(msa->ss->ntuples * sizeof(unsigned));
      for (i = 0; i < msa->ss->ntuples; i++) {
        if (mod->tree_posteriors->base_probs[0][0][node][i] == -1) {
          encoded[i] = code->gap_code;
          ngaps += msa->ss->counts[i];
        }
        else {
          for (j = 0; j < mod->rate_matrix->size; j++) 
            vec_set(v, j, mod->tree_posteriors->base_probs[0][j][node][i]);
          encoded[i] = pbs_get_index(code, v, &error); 
          tot_error += error * msa->ss->counts[i];	 
        }
      }
      vec_free(v);

      /* now write site by site */
      sprintf(out_fname, "%s.%s.bin", out_root, n->name);
      out_f = phast_fopen(out_fname, "w+");
      for (i = 0; i < msa->length; i++) {
        if (keep_gaps || encoded[msa->ss->tuple_idx[i]] != code->gap_code)
          pbs_write_binary(code, encoded[msa->ss->tuple_idx[i]], out_f);
      }

      fprintf(stderr, "Average approximation error ('%s'): %f bits\n",
              n->name, tot_error/(msa->length - ngaps));

      sfree(encoded);
    }
  }

  fprintf(stderr, "Done.\n");
  return 0;
}

/* reconstruct indels by parsimony and assign all base probs to -1
   where ancestral bases are inferred not to have been present */
void do_indels(MSA *msa, TreeModel *mod) {
  int s, tup, i, j;
  TreeNode *n, *lca;
  char c;
  typedef enum {IGNORE, GAP, BASE, MISSING, AMBIG} label_type;
  List *postorder;

  label_type *label = smalloc(mod->tree->nnodes * sizeof(label_type));
  List *inside = lst_new_ptr(mod->tree->nnodes), 
    *outside = lst_new_ptr(mod->tree->nnodes),
    *ambig_cases = lst_new_ptr(mod->tree->nnodes);
  int *seq_to_leaf = smalloc(msa->nseqs * sizeof(int));

  /* build mapping from seqs to leaf indices in tree */
  for (s = 0; s < msa->nseqs; s++) {
    TreeNode *n = tr_get_node(mod->tree, msa->names[s]);
    if (n == NULL)
      die("ERROR: no match for sequence \"%s\" in tree.\n", msa->names[s]);
    seq_to_leaf[s] = n->id;
  }    

  if (mod->msa_seq_idx == NULL)
    tm_build_seq_idx(mod, msa);

  postorder = tr_postorder(mod->tree);

  for (tup = 0; tup < msa->ss->ntuples; tup++) {
    int min = mod->tree->nnodes, max = -1, ngaps = 0, skip_root = FALSE;

    /* find min and max ids of seqs that actually have bases (non-gaps) */
    for (s = 0; s < msa->nseqs; s++) {
      if (ss_get_char_tuple(msa, tup, s, 0) == GAP_CHAR) {
        ngaps++;
        continue;
      }
      if (seq_to_leaf[s] < min) min = seq_to_leaf[s];
      if (seq_to_leaf[s] > max) max = seq_to_leaf[s];

      /* NOTE: missing data being handled like bases here; in some
         cases, a base may be inferred at an ancestral node, when the
         only evidence for it is missing data in the leaves.  There
         are ambiguous cases; we'll err on the side of predicting
         bases rather than indels */
    }

    if (ngaps <= 1) continue;	/* short cut -- impossible to infer
                                   gaps in ancestors */

    else if (ngaps >= msa->nseqs - 1) {
      /* in this case, all ancestors must be gaps */
      for (i = 0; i < mod->tree->nnodes; i++) {
        n = lst_get_ptr(mod->tree->nodes, i);
        if (n->lchild == NULL || n->rchild == NULL) 
          continue;               /* ignore leaves */
        for (j = 0; j < mod->rate_matrix->size; j++)
          mod->tree_posteriors->base_probs[0][j][n->id][tup] = -1;
	/* mark as gap */
      }
      continue;
    }

    if (min < 0) die("prequel.c: min = %e < 0\n", min);
    if (max < min) die("prequel.c: max (%e) < min (%e)", max, min);

    /* the LCA of all leaves with non-gaps must be the first ancestor of
       the node with the max id that has an id smaller than the min
       id.  This is based on the assumption that node ids are assigned
       sequentially in a preorder traversal of the tree, which will be
       true as long as the tree is read from a Newick file by the code
       in trees.c */
    for (lca = lst_get_ptr(mod->tree->nodes, max); lca->id > min; 
         lca = lca->parent);

    /* by parsimony, the base was inserted on the branch to the LCA,
       and all ancestral nodes outside the subtree rooted at the LCA
       did not have bases */

    if (lca == mod->tree->lchild || lca == mod->tree->rchild)
      skip_root = TRUE;        /* don't mark root as gap in this case:
                                  can't distinguish insertion from
                                  deletion so assume deletion */

    /* mark ancestral bases outside subtree beneath LCA as gaps */
    tr_partition_nodes(mod->tree, lca, inside, outside);
    for (i = 0; i < mod->tree->nnodes; i++) label[i] = BASE;
    for (i = 0; i < lst_size(outside); i++) {
      n = lst_get_ptr(outside, i);
      label[n->id] = IGNORE;
      if (n->lchild == NULL || n->rchild == NULL) 
        continue;               /* skip leaves */
      if (n == mod->tree && skip_root) 
        continue;               /* skip root if condition above */
      for (j = 0; j < mod->rate_matrix->size; j++)
        mod->tree_posteriors->base_probs[0][j][n->id][tup] = -1;
      /* mark as gap */
    }

    /* check for gaps in subtree; if there's at most one, we can go
       on; otherwise have to use parsimony to infer history in subtree */
    ngaps = 0;
    for (i = 0; i < lst_size(inside); i++) {
      n = lst_get_ptr(inside, i);
      if (n->lchild == NULL &&
          ss_get_char_tuple(msa, tup, mod->msa_seq_idx[n->id], 0) == GAP_CHAR)
        ngaps++;
    }
    if (ngaps <= 1) continue;

    /* use Dollo parsimony to infer the indel history of the subtree
       beneath the LCA.  Use the fact that every base must have a
       chain of bases to the LCA, because, assuming the alignment is
       correct, no insertions are possible beneath the LCA */
    lst_clear(ambig_cases);
    for (i = 0; i < lst_size(postorder); i++) {
      n = lst_get_ptr(postorder, i);
      if (label[n->id] == IGNORE) continue; /* outside subtree */

      /* MISSING means all leaves beneath node have missing data */
      /* AMBIG means combination of gaps and missing data beneath node */

      else if (n->lchild == NULL) {  /* leaf in subtree */
        c = ss_get_char_tuple(msa, tup, mod->msa_seq_idx[n->id], 0);
        if (c == GAP_CHAR)
          label[n->id] = GAP;
        else if (msa->is_missing[(int)c]) 
          label[n->id] = MISSING;
        else
          label[n->id] = BASE;
      }
      else {                    /* internal node in subtree */
        if (label[n->lchild->id] == BASE || label[n->rchild->id] == BASE)
          label[n->id] = BASE;  /* by Dollo parsimony */
        else if ((label[n->lchild->id] == GAP || label[n->lchild->id] == AMBIG) &&
                 (label[n->rchild->id] == GAP || label[n->rchild->id] == AMBIG))
          label[n->id] = GAP;   /* gaps from both sides and no bases -- must be gap */
        else if (label[n->lchild->id] == MISSING && label[n->rchild->id] == MISSING)
          label[n->id] = MISSING;
        else {              /* must be GAP/MISSING or AMBIG/MISSING */
          label[n->id] = AMBIG;
          lst_push_ptr(ambig_cases, n);
        }
      }
    }

    /* now resolve any ambiguities, by giving each ambiguous node the same
       label as its parent; traversing ambig_cases in reverse order
       ensures that parents are visited before children  */

    /* first make sure root of subtree has a base */
    if (label[lca->id] == MISSING || label[lca->id] == AMBIG)
      label[lca->id] = BASE;
    /* in this case there is all missing data and gaps beneath the LCA;
       hard to know what is right, but let's force a base and err on
       the side of bases rather than gaps */

    for (i = lst_size(ambig_cases) - 1; i >= 0; i--) {
      n = lst_get_ptr(ambig_cases, i);
      if (n == lca) continue;
      else label[n->id] = label[n->parent->id];
    }

    /* now mark gaps inside subtree, as needed */
    for (i = 0; i < lst_size(inside); i++) {
      n = lst_get_ptr(inside, i);
      if (n->lchild == NULL || n->rchild == NULL) continue;
      if (label[n->id] == GAP) 
        for (j = 0; j < mod->rate_matrix->size; j++)
          mod->tree_posteriors->base_probs[0][j][n->id][tup] = -1;
    }
  }

  lst_free(inside);
  lst_free(outside);
  lst_free(ambig_cases);
  sfree(seq_to_leaf);
  sfree(label);
}

