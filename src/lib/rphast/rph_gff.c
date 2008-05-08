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
#include <gff.h>
#include <rph_util.h>

/******************functions defined herein******************/
void rph_gff_read(char** fname, double* gffAddress, int* error, char** errstr);
void rph_gff_print_set(double* address, char** fname, int* error, char** errstr);
void rph_gff_print_features(double* address, char** fname, int* error, char** errstr);
void rph_gff_free(double* address, int* error, char** errstr);

/***********************************************************
rph_gff_read
Reads a gff file and returns the address of the GFF_Set
**********************************************************/
void rph_gff_read(char** fname, double* gffAddress, int* error, char** errstr){
  GFF_Set* gff;
  CategoryMap *cm;
  FILE* file=fopen_fname(*fname, "r");

  gff=gff_read_set(file);

  *gffAddress=ptr2ad(gff);

  fclose(file);

  *error=rphast_errno;
  *errstr=rphast_errmsg;

}

void rph_gff_print_set(double* address, char** fname, int* error, char** errstr){
  GFF_Set* gff;
  FILE* file=fopen_fname(*fname, "w");

  gff=(GFF_Set*)ad2ptr(*address);
  gff_print_set(file, gff);

  fclose(file);

  *error=rphast_errno;
  *errstr=rphast_errmsg;

}

void rph_gff_print_features(double* address, char** fname, int* error, char** errstr){
  GFF_Set* gff;
  FILE* file;
  GFF_Feature* feat;
  int i;
  char filename[1000];

  gff=(GFF_Set*)ad2ptr(*address);
  for(i=0; i<lst_size(gff->features); i++){
    feat=lst_get_ptr(gff->features, i);
    sprintf(filename, "%s_%d", *fname, i);
    file=fopen_fname(filename,"w");
    gff_print_feat(file, feat);
    fclose(file);
  }

  *error=rphast_errno;
  *errstr=rphast_errmsg;

}

void rph_gff_free(double* address, int* error, char** errstr){
  GFF_Set* gff;
  gff=(GFF_Set*)ad2ptr(*address);
  gff_free_set(gff);

  *error=rphast_errno;
  *errstr=rphast_errmsg;

}
