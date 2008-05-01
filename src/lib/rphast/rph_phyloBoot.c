/*****************************************************
rph_phyloBoot.c
The RPHAST handles to phyloBoot, a part
of the phast package.

Alexandra Denby
Last updated: 4/22/08
*****************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <misc.h>
#include <msa.h>
#include <maf.h>
#include <sufficient_stats.h>
#include <numerical_opt.h>
#include <tree_model.h>
#include <fit_em.h>
#include <time.h>
#include <rph_util.h>

/******************functions defined herein******************/


void rph_tm_generate_msa(double* modAddress, int* numSites, double* msaAddress, int* numberSpecies, int* length, char** alphabet, int* error, char** errstr){
  MSA* msa;
  int nsites=*numSites;
  TreeModel* model;
  int i;

  model=(TreeModel*)ad2ptr(*modAddress);

  msa = tm_generate_msa(nsites, NULL, &model, NULL);

  *msaAddress=ptr2ad(msa);
  *numberSpecies=msa->nseqs;
  *length=msa->length;
  *alphabet=msa->alphabet;

  *error=rphast_errno;
  *errstr=rphast_errmsg;

}

