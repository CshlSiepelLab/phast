/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#include "phylo_p.h"
#include "phyloP.help"
#include <misc.h>


int main(int argc, char *argv[]) {
  struct phyloP_struct *p = phyloP_struct_new(0);
  char c;
  FILE *msa_f = NULL;
  msa_format_type msa_format = UNKNOWN_FORMAT;

  /* other variables */
  int opt_idx, seed = -1;
  List *cats_to_do_str=NULL;
  struct timeval now;

  struct option long_opts[] = {
    {"method", 1, 0, 'm'},
    {"mode", 1, 0, 'o'},
    {"msa-format", 1, 0, 'i'},
    {"null", 1, 0, 'n'},
    {"posterior", 0, 0, 'p'},
    {"confidence-interval", 1, 0, 'c'},
    {"subtree", 1, 0, 's'},
    {"branch", 1, 0, 'B'},
    {"features", 1, 0, 'f'},
    {"fit-model", 0, 0, 'F'},
    {"epsilon", 1, 0, 'e'},
    {"quantiles", 0, 0, 'q'},
    {"wig-scores", 0, 0, 'w'},
    {"base-by-base", 0, 0, 'b'},
    {"refidx", 1, 0, 'r'},
    {"chrom", 1, 0, 'N'},
    {"log", 1, 0, 'l'},
    {"gff-scores", 0, 0, 'g'},
    {"do-cats", 1, 0, 'C'},
    {"catmap", 1, 0, 'M'},
    {"no-prune", 0, 0, 'P'},
    {"seed", 1, 0, 'd'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

#ifndef RPHAST
  /* set seed for sampling */
  gettimeofday(&now, NULL);
  srandom(now.tv_usec);
#endif

  while ((c = getopt_long(argc, argv, "m:o:i:n:pc:s:f:Fe:l:r:B:d:qwgbPN:h", 
                          long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'm':
      if (!strcmp(optarg, "SPH"))
        p->method = SPH;
      else if (!strcmp(optarg, "LRT"))
        p->method = LRT;
      else if (!strcmp(optarg, "SCORE"))
        p->method = SCORE;
      else if (!strcmp(optarg, "GERP"))
        p->method = GERP;
      else die("ERROR: bad argument to --method (-m).\n");
      break;
    case 'o':
      if (!strcmp(optarg, "CON"))
        p->mode = CON;
      else if (!strcmp(optarg, "ACC"))
        p->mode = ACC;
      else if (!strcmp(optarg, "NNEUT"))
        p->mode = NNEUT;
      else if (!strcmp(optarg, "CONACC"))
        p->mode = CONACC;
      else die("ERROR: bad argument to --mode (-o).\n");
      break;
    case 'i':
      msa_format = msa_str_to_format(optarg);
      if (msa_format == UNKNOWN_FORMAT)
        die("ERROR: unrecognized alignment format.\n");
      break;
    case 'g':
      p->output_gff = TRUE;
      break;
    case 'n':
      p->nsites = get_arg_int_bounds(optarg, 1, INFTY);
      p->prior_only = TRUE;
      break;
    case 'p':
      p->post_only = TRUE;
      break;
    case 'c':
      p->ci = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 's':
      p->subtree_name = optarg;
      break;
    case 'B':
      p->branch_name = get_arg_list(optarg);
      break;
    case 'f':
      p->feats = gff_read_set(phast_fopen(optarg, "r"));
      break;
    case 'F':
      p->fit_model = TRUE;
      break;
    case 'e':
      p->epsilon = get_arg_dbl_bounds(optarg, 0, 0.1);
      break;
    case 'q':
      p->quantiles_only = TRUE;
      break;
    case 'w':
      p->base_by_base = TRUE;
      p->output_wig = TRUE;
      break;
    case 'b':
      p->base_by_base = TRUE;
      break;
    case 'N':
      p->chrom = optarg;
      break;
    case 'r':
      p->refidx = get_arg_int_bounds(optarg, 0, INFTY);
      break;
    case 'l':
      if (!strcmp(optarg, "-"))
        p->logf = stderr;
      else
        p->logf = phast_fopen(optarg, "w+");
      break;
    case 'C':
      cats_to_do_str = get_arg_list(optarg);
      break;
    case 'M':
      p->cm = cm_new_string_or_file(optarg);
      break;
    case 'd':
      seed = get_arg_int_bounds(optarg, 1, INFTY);
      break;
    case 'P':
      p->no_prune = TRUE;
      break;
    case 'h':
      printf("%s", HELP);
      exit(0);
    case '?':
      die("Bad argument.  Try 'phyloP -h'.\n");
    }
  }

  set_seed(seed);

  if ((p->prior_only && optind > argc - 1) || 
      (!p->prior_only && optind != argc - 2))
    die("ERROR: bad arguments.  Try 'phyloP -h'.\n");
  p->mod_fname = argv[optind];

  p->mod = tm_new_from_file(phast_fopen(p->mod_fname, "r"), 1);

  if (cats_to_do_str != NULL) {
    if (p->cm == NULL) die("ERROR: --cats-to-do requires --catmap option\n");
    p->cats_to_do = cm_get_category_list(p->cm, cats_to_do_str, FALSE);
  }

  if (!p->prior_only) {
    p->msa_fname = argv[optind+1];
    msa_f = phast_fopen(p->msa_fname, "r");
    if (msa_format == UNKNOWN_FORMAT)
      msa_format = msa_format_for_content(msa_f, 1);
    if (msa_format == MAF) 
      p->msa = maf_read_cats(msa_f, NULL, 1, NULL, 
			     p->cats_to_do==NULL ? NULL : p->feats, p->cm, -1, 
			     (p->feats == NULL && p->base_by_base==0) ? FALSE : TRUE, /* --features requires order */
			     NULL, NO_STRIP, FALSE, p->cats_to_do); 
    else 
      p->msa = msa_new_from_file_define_format(msa_f, msa_format, NULL);
    phast_fclose(msa_f);

    /* if base_by_base and undefined chrom, use filename root as chrom */
    if (p->base_by_base && p->chrom == NULL) {
      String *tmpstr = str_new_charstr(p->msa_fname);
      if (str_equals_charstr(tmpstr, "-")) p->chrom = "NA";
      else {
        str_remove_path(tmpstr);
        str_shortest_root(tmpstr, '.');
        p->chrom = tmpstr->chars;    
      }
    }
  }
  
  phyloP(p);    
  return 0;
}

