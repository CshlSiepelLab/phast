/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/


#include "bgc_hmm.h"
#include "phastBias.help"

/* Basic idea: 
   Read in an MSA and a features file describing which regions contain
   codons.  For now assume codons are clean, and always length three.
   Everything else is flanking.  Set up 6 state HMM.  2 states for coding,
   4 states for non-coding.
   The coding states are bgc and non-bgc, each use the site model, with 4
   categories of sites as in Nielsen Yang model.
   The noncoding states are bgc+sel, bgc, sel, neutral.  The rates to/from
   bgc states should be the same everywhere.

   Parameters that describe this HMM:
   tree scale
   noncoding conserved tree scale (hold at 0.3?)
   noncoding transition rate cons->neutral (use phastCons estimates?)
   noncoding transition rate neutral->cons (use phastCons estimates?)
   rate bgc->no_bgc
   rate no_bgc->bgc
   2 weight parameters from nielsen/yang model
   negative selection
   positive selection
   bgc 

   We can use EM to estimate rates between bgc and no bgc states, and
   estimate the rest of the parameters with tm_fit_multi (which can
   optimize parameters shared across multiple models).

   Need to be careful about codon vs noncoding regions, make sure rates
   are on the same scale even though codons are in units of 3 bases.
 */


int main(int argc, char *argv[]) {
  FILE *infile;
  char c;
  int opt_idx;
  msa_format_type msa_format;
  struct bgchmm_struct *b = bgchmm_struct_new(0);

  struct option long_opts[] = {
    {"bgc", 1, 0, 'B'},
    {"estimate-bgc", 1, 0, 'b'},
    {"bgc-exp-length", 1, 0, 'L'},
    {"estimate-bgc-exp-length", 1, 0, 'l'},
    {"bgc-target-coverage", 1, 0, 'C'},
    {"estimate-bgc-target-coverage", 1, 0, 'c'},
    {"rho", 1, 0, 'R'},
    {"cons-exp-length", 1, 0, 'E'},
    {"cons-target-coverage", 1, 0, 'T'},
    {"scale", 1, 0, 'S'},
    {"estimate-scale", 1, 0, 's'},
    {"eqfreqs-from-msa", 1, 0, 'f'},
    {"output-tracts", 1, 0, 'g'},
    {"posteriors", 1, 0, 'p'},
    {"output-mods", 1, 0, 'm'},
    {"help", 0, 0, 'h'},
    {0,0,0,0}};

  while ((c = getopt_long(argc, argv, "B:b:L:l:C:c:R:E:T:S:s:f:g:p:m:Wh", long_opts, &opt_idx))
	 != -1) {
    switch (c) {
    case 'B':
      b->bgc = get_arg_dbl_bounds(optarg, 0, INFTY);
      break;
    case 'b':
      b->estimate_bgc = get_arg_int_bounds(optarg, 0, 1);
      break;
    case 'L':
      b->bgc_expected_length = get_arg_dbl_bounds(optarg, 0, INFTY);
      break;
    case 'l':
      b->estimate_bgc_expected_length = get_arg_int_bounds(optarg, 0, 1);
      break;
    case 'C':
      b->bgc_target_coverage = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'c':
      b->estimate_bgc_target_coverage = get_arg_int_bounds(optarg, 0, 1);
      break;
    case 'R':
      b->rho = get_arg_dbl_bounds(optarg, 0, INFTY);
      if (b->rho > 1.0) phast_warning("Warning: rho is a scale for conserved states and is usually less than 1, got %f", b->rho);
      break;
    case 'E':
      b->cons_expected_length = get_arg_dbl_bounds(optarg, 0, INFTY);
      break;
    case 'T':
      b->cons_target_coverage = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'S':
      b->scale = get_arg_dbl_bounds(optarg, 0, INFTY);
      break;
    case 's':
      b->estimate_scale = get_arg_int_bounds(optarg, 0, 1);
      break;
    case 'f':
      b->eqfreqs_from_msa = get_arg_int_bounds(optarg, 0, 1);
      break;
    case 'g':
      b->tract_fn = optarg;
      break;
    case 'p':
      if (strcasecmp(optarg, "none")==0)
	b->post_probs = NONE;
      else if (strcasecmp(optarg, "wig")==0)
	b->post_probs = WIG;
      else if (strcasecmp(optarg, "full")==0)
	b->post_probs = FULL;
      else die("--posteriors option expects either none, wig, or full, got %s", optarg);
      break;
    case 'm':
      b->mods_fn = optarg;
      break;
    case 'W':
      b->post_probs = NONE;
      break;
    case 'h':
      printf("%s", HELP);
      exit(0);
    case '?':
    default:
      die("Bad argument.  Try '%s -h'.\n", argv[0]);
    }
  }
   
  if (optind != argc - 3)
    die("ERROR: extra or missing arguments.  Try '%s -h'.\n", argv[0]);


  /* read the MSA and make sure there are no sufficient stats (we want to get these
     later depending on whether coding/noncoding, tuple size will be 3 or 1)
   */
  infile = phast_fopen(argv[optind], "r");
  msa_format = msa_format_for_content(infile, 1);
  if (msa_format == MAF) 
    b->msa = maf_read(infile, NULL, 1, NULL, NULL, NULL, 0, 1, NULL, NO_STRIP, 0);
  else b->msa = msa_new_from_file_define_format(infile, msa_format, NULL);
  fclose(infile);

  infile = phast_fopen(argv[optind+1], "r");
  b->mod = tm_new_from_file(infile, 1);
  fclose(infile);

  b->foregd_branch = argv[optind+2];
  
  bgcHmm(b);
  
  return 0;
}

