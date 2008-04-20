/*****************************************************
rph_phyloBoot.c
The RPHAST handles to function from phyloBoot, a part
of the phast package.

Alexandra Denby
4/12/08
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



void rph_tm_generate_msa(int* modAddress, int* numSites, int* msaAddress, int* numberSpecies, int* length, char** alphabet, char** species, int* error){
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
  for(i=0; i<*numberSpecies; i++){
    species[i]=msa->names[i];
  }


}
