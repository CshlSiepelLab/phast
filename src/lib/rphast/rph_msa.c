/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2009 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

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
void rph_msa_read(char** fname, char** format, double* address, double* gffAddress, int* fourD, int* SS, int* numberSpecies, int* length, char** alphabet, int* error, char** errstr);
void readSpecies(double* address, int* numberSpecies, char** species, int* error, char** errstr);
void rph_msa_print(char** fname, char** format, double* address, int* error, char** errstr);
void rph_msa_free(double* address, int* error, char** errstr);
void reduce_to_4d(MSA *msa, CategoryMap *cm);
void rph_msa_concatenate(double* baseAddress, double* newAddress, int* error, char** errstr);
void rph_msa_new(int* nseqs, char** spec, char** seqs, double* address, int* numberSpecies, int* length, char** alphabet, int* error, char** errstr);

/*******************************************************
rph_msa_read
reads an alignment from a file.
*******************************************************/
void rph_msa_read(char** fname, char** format, double* address, double* gffAddress, int* fourD, int* SS, int* numberSpecies, int* length, char** alphabet, int* error, char** errstr){
  MSA* msa;
  CategoryMap *cm;
  GFF_Set *gff;
  FILE* file=fopen_fname(*fname, "r");
  msa_format_type input_format;
  int tuple_size=1;
  char* reverse_groups_tag=NULL;

  input_format = msa_str_to_format(*format);
  file=fopen_fname(*fname,"r");
  
  if (*gffAddress==0){
    gff=NULL;
    cm=NULL;
  }
  else{
    gff=(GFF_Set*)ad2ptr(*gffAddress);
    cm=cm_new_from_features(gff);
  }

  if (*fourD){
    if (gff==NULL){
      *error=1;
      strcpy(*errstr, "Error: 4D requires gff\n");
      return;
    }
    tuple_size=3;
    reverse_groups_tag="transcript_id";
    cm = cm_new_string_or_file("NCATS=3; CDS 1-3");
  }

  if(input_format==MAF){
    msa=maf_read(file, NULL,tuple_size,NULL, gff, cm, -1, FALSE, reverse_groups_tag, NO_STRIP, FALSE);
  }
  else{
    msa = msa_new_from_file(file, input_format, NULL);
  }

  if (*SS){
    if (msa->ss==NULL){
      ss_from_msas(msa, tuple_size, FALSE, NULL, NULL, NULL, -1);
    }
  }

  if (*fourD){
    reduce_to_4d(msa, cm);
  }

  *address=ptr2ad(msa);
  *numberSpecies=msa->nseqs;
  *length=msa->length;
  *alphabet=msa->alphabet;
  fclose(file);

  *error=rphast_errno;
  *errstr=rphast_errmsg;
  
}


void rph_msa_new(int* nseqs, char** spec, char** seqs, double* address, int* numberSpecies, int* length, char** alphabet, int* error, char** errstr){

  int len=0;
  char* al=NULL; //alphabet

  MSA* msa=msa_new(seqs, spec, *nseqs, len, al);

  *address=ptr2ad(msa);
  *numberSpecies=msa->nseqs;
  *length=msa->length;
  
  *alphabet=msa->alphabet;

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




void rph_msa_concatenate(double* baseAddress, double* newAddress, int* error, char** errstr){
  MSA* base;
  MSA* new;

  base=(MSA*)ad2ptr(*baseAddress);
  new=(MSA*)ad2ptr(*newAddress);

  msa_concatenate(base, new);

  *error=rphast_errno;
  *errstr=rphast_errmsg;

}






/* reduce SS representation of alignment with nucleotide triples to 4d
   sites only (supports --4d option) */
void reduce_to_4d(MSA *msa, CategoryMap *cm) {
  String *tmpstr = str_new_charstr("CDS");
  int cat_pos3 = cm->ranges[cm_get_category(cm, tmpstr)]->end_cat_no;
  int i, j, len = 0;
  if (cat_pos3 == 0)
    die("ERROR: no match for 'CDS' feature type (required with --4d).\n");
  assert(msa->ss->cat_counts != NULL && msa->ncats >= cat_pos3);
  for (i = 0; i < msa->ss->ntuples; i++) {
    if (ss_is_4d(msa, i)) {
      msa->ss->counts[i] = msa->ss->cat_counts[cat_pos3][i];
      len += msa->ss->counts[i];
    }
    else 
      msa->ss->counts[i] = 0;
  }
  for (j = 0; j <= msa->ncats; j++) free(msa->ss->cat_counts[j]);
  free(msa->ss->cat_counts);
  msa->ss->cat_counts = NULL;
  msa->ncats = -1;
  msa->length = len;
  ss_remove_zero_counts(msa);
  str_free(tmpstr);
}
