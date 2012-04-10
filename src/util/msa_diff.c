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
#include <msa.h>
#include <maf.h>
#include "msa_diff.help" 

int main(int argc, char *argv[]) {
  char c;
  int opt_idx;
  MSA *msa1, *msa2;
  int i, j, ncommon, len, ndiffs, lastsame, same, diffstart, col;
  char **common_names;
  int *common_to_msa1, *common_to_msa2, *mark;

  int ignore_base_id = FALSE, ignore_gap_type = FALSE;
  char *alphabet = "ACGT";
  msa_format_type format1 = UNKNOWN_FORMAT, format2 = UNKNOWN_FORMAT;
  FILE *infile1, *infile2;

  struct option long_opts[] = {
    {"ignore-base-id", 0, 0, 'b'},
    {"ignore-gap-type", 0, 0, 'g'},
    {"alphabet", 1, 0, 'a'},
    {"format1", 1, 0, 'i'},
    {"format2", 1, 0, 'j'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "bga:i:j:h", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'b':
      ignore_base_id = TRUE;
      break;
    case 'g':
      ignore_gap_type = TRUE;
      break;
    case 'a':
      alphabet = optarg;
      break;
    case 'i':
      format1 = msa_str_to_format(optarg);
      if (format1 == UNKNOWN_FORMAT) 
        die("ERROR: bad input format.\n");
      break;      
    case 'j':
      format2 = msa_str_to_format(optarg);
      if (format2 == UNKNOWN_FORMAT) 
        die("ERROR: bad input format.\n");
      break;      
    case 'h':
      printf("%s", HELP); 
      exit(0);
    case '?':
      die("Bad argument.  Try 'program -h'.\n");
    }
  }

  if (optind != argc - 2) 
    die("Two filenames required.  Try 'msa_diff -h'.\n");

  set_seed(-1);
  infile1 = phast_fopen(argv[optind], "r");
  if (format1 == UNKNOWN_FORMAT)
    format1 = msa_format_for_content(infile1, 1);
  if (format1 == MAF) 
    msa1 = maf_read(phast_fopen(argv[optind], "r"), NULL, 1, alphabet,
                    NULL, NULL, -1, TRUE, NULL, NO_STRIP, FALSE);
  else
    msa1 = msa_new_from_file_define_format(phast_fopen(argv[optind], "r"), 
                             format1, alphabet);

  infile2 = phast_fopen(argv[optind+1], "r");
  if (format2 == UNKNOWN_FORMAT)
    format2 = msa_format_for_content(infile2, 1);
  if (format2 == MAF)
    msa2 = maf_read(phast_fopen(argv[optind+1], "r"), NULL, 1, alphabet,
                    NULL, NULL, -1, TRUE, NULL, NO_STRIP, FALSE);
  else 
    msa2 = msa_new_from_file_define_format(phast_fopen(argv[optind+1], "r"), 
                             format2, alphabet);

  common_names = smalloc(msa1->nseqs * sizeof(void*));
  common_to_msa1 = smalloc(msa1->nseqs * sizeof(int));
  common_to_msa2 = smalloc(msa1->nseqs * sizeof(int));
  mark = smalloc(msa2->nseqs * sizeof(int));
  for (j = 0; j < msa2->nseqs; j++) mark[j] = FALSE;

  ncommon = 0;
  for (i = 0; i < msa1->nseqs; i++) {
    j = msa_get_seq_idx(msa2, msa1->names[i]);
    if (j < 0)
      printf("Sequence '%s' from msa1 not found in msa2\n", msa1->names[i]);
    else {
      common_names[ncommon] = msa1->names[i];
      common_to_msa1[ncommon] = i;
      common_to_msa2[ncommon] = j;
      mark[j] = TRUE;
      ncommon++;
    }
  }
  for (j = 0; j < msa2->nseqs; j++)
    if (!mark[j])
      printf("Sequence '%s' from msa2 not found in msa1\n", msa2->names[j]);

  len = min(msa1->length, msa2->length);
  if (msa1->length != msa2->length) 
    printf("Lengths differ (msa1: %d, msa2: %d); comparing common part only\n",
           msa1->length, msa2->length);

  ndiffs = 0;
  printf("Columns with differences:\n");
  lastsame = TRUE;
  diffstart = 0;
  for (col = 0; col < len; col++) {
    same = TRUE;
    for (i = 0; same && i < ncommon; i++) {
      char c1 = msa_get_char(msa1, common_to_msa1[i], col);
      char c2 = msa_get_char(msa2, common_to_msa2[i], col);
      if (ignore_base_id) {
        if (isalpha(c1)) c1 = 'b';
        if (isalpha(c2)) c2 = 'b';
      }
      if (ignore_gap_type) {
        if (c1 == '-' || c1 == '.' || c1 == '^') c1 = '-';
        if (c2 == '-' || c2 == '.' || c2 == '^') c2 = '-';
      }
      if (c1 != c2) same = FALSE;
    } 

    if (same == FALSE && lastsame == TRUE) {
      ndiffs++;
      diffstart = col;
    }
    else if (same == TRUE && lastsame == FALSE)
      printf("%d - %d\n", diffstart, col-1);

    lastsame = same;
  }
  if (lastsame == FALSE)
    printf("%d - %d\n", diffstart, col-1);  

  if (ndiffs == 0) printf("none\n");

  return 0;
}
