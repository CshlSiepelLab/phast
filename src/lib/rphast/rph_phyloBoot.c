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

/*********DELETE ME LATER*********/
int rphast_errno;
char* rphast_errmsg;

/******************functions defined herein******************/


void rph_tm_generate_msa(int* modAddress, int* numSites, int* msaAddress, int* numberSpecies, int* length, char** alphabet, char** species, int* error, char** errstr){
  MSA* msa;
  int nsites=numSites[0];
  TreeModel* model;
  int i;

  model=(TreeModel*)modAddress[0];

  msa = tm_generate_msa(nsites, NULL, &model, NULL);

  *msaAddress=(unsigned int)msa;
  *numberSpecies=msa->nseqs;
  *length=msa->length;
  *alphabet=msa->alphabet;

  *error=rphast_errno;
  *errstr=rphast_errmsg;

}
