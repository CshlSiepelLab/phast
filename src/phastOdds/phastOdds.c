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
#include <tree_model.h>
#include <hmm.h>
#include <msa.h>
#include <maf.h>
#include <sufficient_stats.h>
#include <gff.h>
#include <bed.h>
#include <tree_likelihoods.h>
#include "phastOdds.help"

#define MIN_BLOCK_SIZE 30
/* used in identifying regions of missing data in a reference-sequence
   alignment */

int main(int argc, char *argv[]) {
  char c;
  List *l;
  int i, j, strand, bed_output = 0, backgd_nmods = -1, feat_nmods = -1, 
    winsize = -1, verbose = 0, max_nmods, memblocksize, old_nleaves,
    refidx = 1, base_by_base = FALSE, windowWig = FALSE;
  TreeModel **backgd_mods = NULL, **feat_mods = NULL;
  HMM *backgd_hmm = NULL, *feat_hmm = NULL;
  msa_format_type inform = UNKNOWN_FORMAT;
  GFF_Set *features = NULL;
  MSA *msa, *msa_compl=NULL;
  double **backgd_emissions, **feat_emissions, **mem, **dummy_emissions,
    *winscore_pos=NULL, *winscore_neg=NULL;
  int *no_alignment=NULL;
  List *pruned_names;
  char *msa_fname;
  FILE *infile;

  int opt_idx;
  struct option long_opts[] = {
    {"background-mods", 1, 0, 'b'},
    {"background-hmm", 1, 0, 'B'},
    {"feature-mods", 1, 0, 'f'},
    {"feature-hmm", 1, 0, 'F'},
    {"features", 1, 0, 'g'},
    {"window", 1, 0, 'w'},
    {"window-wig", 1, 0, 'W'},
    {"base-by-base", 0, 0, 'y'},
    {"msa-format", 1, 0, 'i'},
    {"refidx", 1, 0, 'r'},
    {"output-bed", 0, 0, 'd'},
    {"verbose", 0, 0, 'v'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "B:b:F:f:r:g:w:W:i:ydvh", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'B':
      backgd_hmm = hmm_new_from_file(phast_fopen(optarg, "r"));
      break;
    case 'b':
      l = get_arg_list(optarg);
      backgd_nmods = lst_size(l);
      backgd_mods = smalloc(backgd_nmods * sizeof(void*));
      for (i = 0; i < backgd_nmods; i++) 
        backgd_mods[i] = tm_new_from_file(phast_fopen(((String*)lst_get_ptr(l, i))->chars, "r"), 1);
      lst_free_strings(l); lst_free(l);
      break;
    case 'F':
      feat_hmm = hmm_new_from_file(phast_fopen(optarg, "r"));
      break;
    case 'f':
      l = get_arg_list(optarg);
      feat_nmods = lst_size(l);
      feat_mods = smalloc(feat_nmods * sizeof(void*));
      for (i = 0; i < feat_nmods; i++) 
        feat_mods[i] = tm_new_from_file(phast_fopen(((String*)lst_get_ptr(l, i))->chars, "r"), 1);
      lst_free_strings(l); lst_free(l);
      break;
    case 'g':
      features = gff_read_set(phast_fopen(optarg, "r"));
      break;
    case 'w':
      winsize = get_arg_int(optarg);
      if (winsize <= 0) die("ERROR: window size must be positive.\n");
      break;
    case 'W':
      winsize = get_arg_int(optarg);
      if (winsize <= 0) die("ERROR: window size must be positive.\n");
      windowWig = TRUE;
      break;      
    case 'y':
      base_by_base = TRUE;
      break;
    case 'i':
      inform = msa_str_to_format(optarg);
      if (inform == UNKNOWN_FORMAT) die("Bad argument to -i.\n");
      break;
    case 'r':
      refidx = get_arg_int_bounds(optarg, 0, INFTY);
      break;
    case 'd':
      bed_output = 1;
      break;
    case 'h':
      printf("%s", HELP);
      exit(0);
    case 'v':
      verbose = 1;
      break;
    case '?':
      die("Bad argument.  Try '%s -h'.\n", argv[0]);
    }
  }

  set_seed(-1);

  if (backgd_mods == NULL || feat_mods == NULL) 
    die("ERROR: -b and -f required.  Try '%s -h'.\n", argv[0]);

  if (backgd_nmods == 1 && backgd_hmm == NULL) 
    backgd_hmm = hmm_create_trivial();
  else if (backgd_hmm == NULL)
    die("ERROR: -B required.  Try '%s -h'.\n", argv[0]);

  if (feat_nmods == 1 && feat_hmm == NULL) 
    feat_hmm = hmm_create_trivial();
  else if (feat_hmm == NULL)
    die("ERROR: -F required.  Try '%s -h'.\n", argv[0]);

  if ((winsize == -1 && features == NULL && !base_by_base) || 
      (winsize != -1 && features != NULL) || 
      (winsize != -1 && base_by_base) || 
      (features != NULL && base_by_base))
    die("ERROR: must specify exactly one of -g, -w, and -y.  Try '%s -h'.\n", argv[0]);

  if (backgd_hmm->nstates != backgd_nmods) 
    die("ERROR: number of states must equal number of tree models for background.\n");

  if (feat_hmm->nstates != feat_nmods) 
    die("ERROR: number of states must equal number of tree models for features.\n");

  if (features != NULL && lst_size(features->features) == 0)
    die("ERROR: empty features file.\n");

  if (base_by_base && (backgd_nmods > 1 || feat_nmods > 1))
    die("ERROR: only single phylogenetic models (not HMMs) are supported with --base-by-base.\n");

  if (optind != argc - 1) 
    die("ERROR: too few arguments.  Try '%s -h'.\n", argv[0]);

  if (verbose) fprintf(stderr, "Reading alignment ...\n");
  msa_fname = argv[optind];
  infile = phast_fopen(msa_fname, "r");
  if (inform == UNKNOWN_FORMAT)
    inform = msa_format_for_content(infile, 1);
  if (inform == MAF)
    msa = maf_read(infile, NULL, 1, NULL, NULL, 
                   NULL, -1, TRUE, NULL, NO_STRIP, FALSE);
  else
    msa = msa_new_from_file_define_format(infile, inform, NULL);
  if (msa_alph_has_lowercase(msa)) msa_toupper(msa); 
  msa_remove_N_from_alph(msa);

  /* need ordered representation of alignment */
  if (msa->seqs == NULL && (msa->ss == NULL || msa->ss->tuple_idx == NULL) )
    die("ERROR: ordered sufficient statistics are required.\n");

  pruned_names = lst_new_ptr(msa->nseqs);
  for (i = 0; i < backgd_nmods; i++) {
    old_nleaves = (backgd_mods[i]->tree->nnodes + 1) / 2;
    tm_prune(backgd_mods[i], msa, pruned_names);
    if (lst_size(pruned_names) >= old_nleaves)
      die("ERROR: no match for leaves of tree in alignment (background model #%d)\n", i+1);
    else if (lst_size(pruned_names) > 0) {
      fprintf(stderr, "WARNING: pruned away leaves in background model (#%d) with no match in alignment (", i+1);
      for (j = 0; j < lst_size(pruned_names); j++)
        fprintf(stderr, "%s%s", ((String*)lst_get_ptr(pruned_names, j))->chars, 
                j < lst_size(pruned_names) - 1 ? ", " : ").\n");
    }
    lst_free_strings(pruned_names);
  }
  for (i = 0; i < feat_nmods; i++) {
    old_nleaves = (feat_mods[i]->tree->nnodes + 1) / 2;
    tm_prune(feat_mods[i], msa, pruned_names);
    if (lst_size(pruned_names) >= old_nleaves)
      die("ERROR: no match for leaves of tree in alignment (features model #%d)\n", i+1);
    else if (lst_size(pruned_names) > 0) {
      fprintf(stderr, "WARNING: pruned away leaves in features model (#%d) with no match in alignment (", i+1);
      for (j = 0; j < lst_size(pruned_names); j++)
        fprintf(stderr, "%s%s", ((String*)lst_get_ptr(pruned_names, j))->chars, 
                j < lst_size(pruned_names) - 1 ? ", " : ").\n");
    }
    lst_free_strings(pruned_names);
  }
  lst_free(pruned_names);

  /* first have to subtract offset from features, if necessary */
  if (msa->idx_offset != 0 && features != NULL) {
    for (i = 0; i < lst_size(features->features); i++) {
      GFF_Feature *f = lst_get_ptr(features->features, i);
      f->start -= msa->idx_offset;
      f->end -= msa->idx_offset;
    }
  }

  /* convert to coord frame of alignment */
  if (features != NULL && refidx != 0) {
    if (verbose) fprintf(stderr, "Mapping coordinates ...\n");
    msa_map_gff_coords(msa, features, refidx, 0, 0); 
    if (lst_size(features->features) == 0)
      die("ERROR: no features within coordinate range of alignment.\n");
  }

  /* Make a reverse complemented copy of the alignment.  The two
     strands will be processed separately, to avoid problems with
     overlapping features, etc. */
  if (!base_by_base) {          /* skip in base by base case */
    if (verbose) fprintf(stderr, "Creating reverse complemented alignment ...\n");
    msa_compl = msa_create_copy(msa, 0);
    /* temporary workaround: make sure reverse complement not based on
       sufficient stats */
    if (msa_compl->seqs == NULL) ss_to_msa(msa_compl);
    if (msa_compl->ss != NULL) {
      ss_free(msa_compl->ss);
      msa_compl->ss = NULL;
    }
    msa_reverse_compl(msa_compl);
  }

  /* allocate memory for computing scores */
  backgd_emissions = smalloc(backgd_nmods * sizeof(void*));
  for (i = 0; i < backgd_nmods; i++) 
    backgd_emissions[i] = smalloc(msa->length * sizeof(double));
  feat_emissions = smalloc(feat_nmods * sizeof(void*));
  for (i = 0; i < feat_nmods; i++) 
    feat_emissions[i] = smalloc(msa->length * sizeof(double));
  max_nmods = max(backgd_nmods, feat_nmods);
  dummy_emissions = smalloc(max_nmods * sizeof(void*));
  mem = smalloc(max_nmods * sizeof(void*));
  /* memory for forward algorithm -- each block must be as large as
     the largest feature */
  if (features != NULL) {
    for (i = 0, memblocksize = -1; i < lst_size(features->features); i++) {
      GFF_Feature *f = lst_get_ptr(features->features, i);
      if (f->end - f->start + 1 > memblocksize) 
        memblocksize = f->end - f->start + 1;
    }
  }
  else memblocksize = winsize;  /* -1 if base-by-base mode */

  if (memblocksize > 0)
    for (i = 0; i < max_nmods; i++)
      mem[i] = smalloc(memblocksize * sizeof(double));

  if (winsize != -1) {
    winscore_pos = smalloc(msa->length * sizeof(double));
    winscore_neg = smalloc(msa->length * sizeof(double));
    no_alignment = smalloc(msa->length * sizeof(int));

    for (i = 0; i < msa->length; i++) {
      winscore_pos[i] = winscore_neg[i] = NEGINFTY; 
      if (refidx == 0)
        no_alignment[i] = FALSE;
      else
        no_alignment[i] = msa_missing_col(msa, refidx, i);
    }
  }

  /* the rest will be repeated for each strand */
  for (strand = 1; strand <= 2; strand++) {
    MSA *thismsa = strand == 1 ? msa : msa_compl;
    double *winscore = strand == 1 ? winscore_pos : winscore_neg;

    if (base_by_base && strand == 2) break; /* don't do second pass in
                                               base_by_base case */

    if (verbose) fprintf(stderr, "Processing %c strand ...\n",
                         strand == 1 ? '+' : '-');

    /* set up dummy categories array, so that emissions are only
       computed where needed */
    thismsa->categories = smalloc(thismsa->length * sizeof(int));
    thismsa->ncats = 1;
    if (winsize != -1) {
      if (strand == 1)
        for (i = 0; i < thismsa->length; i++) 
          thismsa->categories[i] = no_alignment[i] ? 0 : 1;
      else
        for (i = 0; i < thismsa->length; i++) 
          thismsa->categories[i] = no_alignment[thismsa->length - i - 1] ? 0 : 1;
    }
    else if (features != NULL) {
      for (i = 0; i < thismsa->length; i++) thismsa->categories[i] = 0;
      for (i = 0; i < lst_size(features->features); i++) {
        GFF_Feature *f = lst_get_ptr(features->features, i);
        if (f->start <= 0 || f->end <= 0) {
          fprintf(stderr, "WARNING: feature out of range ('");
          gff_print_feat(stderr, f);
          fprintf(stderr, "')\n");
          continue;
        }

        if (strand == 1 && f->strand != '-') 
          for (j = f->start - 1; j < f->end; j++)
            thismsa->categories[j] = 1;
        else if (strand == 2 && f->strand == '-')
          for (j = thismsa->length - f->end; 
               j < thismsa->length - f->start + 1; j++)
            thismsa->categories[j] = 1;
      }
    }
    else {                      /* base-by-base scores */
      for (i = 0; i < thismsa->length; i++) thismsa->categories[i] = 1;
    }
    if (thismsa->ss != NULL) ss_update_categories(thismsa);

    /* compute emissions */
    for (i = 0; i < backgd_nmods; i++) {
      if (verbose) 
        fprintf(stderr, "Computing emissions for background model #%d ...\n", i+1);
      tl_compute_log_likelihood(backgd_mods[i], thismsa, 
                                backgd_emissions[i], NULL, 1, NULL);
    }
    for (i = 0; i < feat_nmods; i++) {
      if (verbose) 
        fprintf(stderr, "Computing emissions for features model #%d ...\n", i+1);
      tl_compute_log_likelihood(feat_mods[i], thismsa, 
                                feat_emissions[i], NULL, 1, NULL);
    }

    /* now compute scores */
    if (winsize != -1) {        /* windows case */
      int winstart;
      if (verbose) fprintf(stderr, "Computing scores ...\n");

      for (winstart = 0; winstart <= thismsa->length - winsize; winstart++) {
        int centeridx = winstart + winsize/2;

        if (strand == 2) centeridx = thismsa->length - centeridx - 1;

        if (no_alignment[centeridx]) continue;

        for (j = 0; j < feat_nmods; j++)
          dummy_emissions[j] = &(feat_emissions[j][winstart]);
        winscore[centeridx] = hmm_forward(feat_hmm, dummy_emissions, 
                                          winsize, mem);

        if (winscore[centeridx] <= NEGINFTY) {
          winscore[centeridx] = NEGINFTY;
          continue;
        }

        for (j = 0; j < backgd_nmods; j++)
          dummy_emissions[j] = &(backgd_emissions[j][winstart]);
        winscore[centeridx] -= hmm_forward(backgd_hmm, dummy_emissions, 
                                           winsize, mem);

        if (winscore[centeridx] < NEGINFTY) winscore[centeridx] = NEGINFTY;
      }
    }
    else if (features != NULL) { /* features case */
      if (verbose) fprintf(stderr, "Computing scores ...\n");
      for (i = 0; i < lst_size(features->features); i++) {
        GFF_Feature *f = lst_get_ptr(features->features, i);
        int s, e;

        if ((strand == 1 && f->strand == '-') || 
            (strand == 2 && f->strand != '-') ||
            f->start <= 0 || f->end <= 0 || f->end - f->start < 0)
          continue;
        
        /* effective coords */
        if (f->strand == '-') {   
          s = thismsa->length - f->end + 1;
          e = thismsa->length - f->start + 1;
        }
        else { s = f->start; e = f->end; }
        
        f->score_is_null = 0;
        
        for (j = 0; j < feat_nmods; j++)
          dummy_emissions[j] = &(feat_emissions[j][s-1]);
        f->score = hmm_forward(feat_hmm, dummy_emissions, e - s + 1, mem);
        
        if (f->score <= NEGINFTY) {
          f->score = NEGINFTY;
          continue;
        }
        
        for (j = 0; j < backgd_nmods; j++)
          dummy_emissions[j] = &(backgd_emissions[j][s-1]);
        f->score -= hmm_forward(backgd_hmm, dummy_emissions, e - s + 1, mem);

        if (f->score < NEGINFTY) f->score = NEGINFTY;
      }
    }
  }

  if (verbose) fprintf(stderr, "Generating output ...\n");
  
  if (winsize != -1 && windowWig == FALSE) { /* standard windows output */
    for (i = 0, j = 0; i < msa->length; i++) {
      if (no_alignment[i] == FALSE)
        printf("%d\t%.3f\t%.3f\n", j + msa->idx_offset + 1, winscore_pos[i], 
               winscore_neg[i]);
      if (ss_get_char_pos(msa, i, 0, 0) != GAP_CHAR) j++;
    }
  }
  else if (windowWig == TRUE) { /* windows with wig output */
    int last = NEGINFTY;
    for (i = 0, j = 0; i < msa->length; i++) {
      if (refidx == 0 || msa_get_char(msa, refidx-1, i) != GAP_CHAR) {
        if (no_alignment[i] == FALSE && winscore_pos[i] > NEGINFTY) {
          if (j > last + 1) 
            printf("fixedStep chrom=%s start=%d step=1\n", 
                   refidx > 0 ? msa->names[refidx-1] : "alignment",
                   j + msa->idx_offset + 1);
          printf("%.3f\n", winscore_pos[i]);
          last = j;
        }
        j++;
      }
    }
  }
  else if (features != NULL) {  /* features output */
    /* return to coord frame of reference seq (also, replace offset) */
    if (refidx != 0)
      msa_map_gff_coords(msa, features, 0, refidx, msa->idx_offset); 
    else if (msa->idx_offset != 0) {
      for (i = 0; i < lst_size(features->features); i++) {
        GFF_Feature *f = lst_get_ptr(features->features, i);
        f->start += msa->idx_offset;
        f->end += msa->idx_offset;
      }
    }

    if (bed_output) 
      gff_print_bed(stdout, features, FALSE);
    else
      gff_print_set(stdout, features);
  }
  else {           /* base-by-base scores */
    /* in this case, we can just output the difference between the emissions */
    printf("fixedStep chrom=%s start=%d step=1\n", 
           refidx > 0 ? msa->names[refidx-1] : "alignment",
           msa->idx_offset + 1);
    for (i = 0, j = 0; i < msa->length; i++) {
      if (refidx == 0 || msa_get_char(msa, refidx-1, i) != GAP_CHAR) {
        printf("%.3f\n", feat_emissions[0][i] - backgd_emissions[0][i]);
        j++;
      }
    }
  }

  if (verbose) fprintf(stderr, "\nDone.\n");

  return 0;
}
