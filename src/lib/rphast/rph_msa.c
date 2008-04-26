/*****************************************************
rph_msa.c
The RPHAST handles to functions dealing with multiple
sequence alignment functions from the phast package.

Alexandra Denby
Last updated: 4/22/08
*****************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <msa.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include <ctype.h>
#include <sufficient_stats.h>
#include <local_alignment.h>
#include <maf.h>
#include <rph_util.h>

/******************functions defined herein******************/
void rph_msa_read(char** fname, char** format, double* address, int* numberSpecies, int* length, char** alphabet, int* error, char** errstr);
void readSpecies(double* address, int* numberSpecies, char** species, int* error, char** errstr);
void rph_msa_print(char** fname, char** format, double* address, int* error, char** errstr);
void rph_msa_free(double* address, int* error, char** errstr);


/*******************************************************
rph_msa_read
reads an alignment from a file.
*******************************************************/
void rph_msa_read(char** fname, char** format, double* address, int* numberSpecies, int* length, char** alphabet, int* error, char** errstr){
  MSA* msa;
  FILE* file;
  msa_format_type input_format;

  input_format = msa_str_to_format(format[0]);
  file=fopen_fname(*fname,"r");
  if (input_format == MAF) {
    msa = maf_read(file, NULL, 1, 
                   NULL, NULL, NULL, -1, 1, NULL, NO_STRIP, FALSE);
  }
  else{
    msa = msa_new_from_file(fopen_fname(fname[0], "r"), input_format, NULL);
  }

  *address=ptr2ad(msa);
  *numberSpecies=msa->nseqs;
  *length=msa->length;
  *alphabet=msa->alphabet;
  fclose(file);

  *error=rphast_errno;
  *errstr=rphast_errmsg;
  
}

void readSpecies(double* address, int* numberSpecies, char** species, int* error, char** errstr){
  MSA* msa=(MSA*)ad2ptr(*address);
  int i;

  for(i=0; i<*numberSpecies; i++){
    species[i]=msa->names[i];
  }

  *error=rphast_errno;
  *errstr=rphast_errmsg;


}

/*******************************************************
rph_msa_print
writes an alignment to a file.
*******************************************************/
void rph_msa_print(char** fname, char** format, double* address, int* error, char** errstr){
  FILE* outfile=fopen_fname(fname[0],"w");
  msa_format_type output_format;

  output_format=msa_str_to_format(format[0]);
  msa_print(outfile, (MSA*)ad2ptr(*address), output_format, FALSE);
  fclose(outfile);

  *error=rphast_errno;
  *errstr=rphast_errmsg;

}

/*******************************************************
rph_msa_free
frees a multiple sequence alignment.
*******************************************************/
void rph_msa_free(double* address, int* error, char** errstr){ 
  msa_free((MSA*)ad2ptr(*address));

  *error=rphast_errno;
  *errstr=rphast_errmsg;

}
