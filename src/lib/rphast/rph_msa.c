/* $Id: rph_msa.c,v 1.1 2008-04-01 15:19:32 acs Exp $
   Written by Adam Siepel, 2002
   Copyright 2002, Adam Siepel, University of California */

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

/* minimum number of codons required for -L */
#define MIN_NCODONS 10

/* number of bases considered to border each indel and minimum number
   of gapless bases required for -I */
#define INDEL_BORDER 3
#define MIN_NBASES 15

/* new functions */
void init();
void readCommandLineArgs(int* argc, char** argv);
void readAlignment();
void print_usage();
void freeAlignment();
void writeAlignment();
void printSequence();

/* global variables, used for keeping around during access from R */
  MSA *msa, *sub_msa;
  msa_format_type input_format, output_format;
  List *seqlist_str, *l;
  char *infname, *clean_seqname, *rseq_fname,
    *reverse_groups_tag, *alphabet;
  int opt_idx, startcol, endcol, include, gap_strip_mode,
    pretty_print, refseq, tuple_size, ordered_stats, 
    indel_clean_nseqs, cats_done, rand_perm, reverse_compl, 
    stats_only, win_size, cycle_size, maf_keep_overlapping, 
    collapse_missing, fourD, mark_missing_maxsize, 
    missing_as_indels, unmask;
  List *cats_to_do, *aggregate_list, *msa_fname_list, 
    *order_list, *fill_N_list;
  msa_coord_map *map;
  GFF_Set *gff;
  CategoryMap *cm;

/* initialize our globals */
void init(){
  msa = NULL;
  sub_msa = NULL;
  input_format = FASTA;
  output_format = FASTA;
  seqlist_str = NULL;
  l = NULL;
  infname = NULL;
  clean_seqname = NULL; 
  rseq_fname = NULL;
  reverse_groups_tag = NULL; 
  alphabet = NULL;
  startcol = 1; 
  endcol = -1; 
  include = 1; 
  gap_strip_mode = NO_STRIP;
  pretty_print = FALSE;  
  refseq = 0; 
  tuple_size = 1; 
  ordered_stats = TRUE; 
  indel_clean_nseqs = -1; 
  cats_done = FALSE;
  rand_perm = FALSE; 
  reverse_compl = FALSE; 
  stats_only = FALSE; 
  win_size = -1; 
  cycle_size = -1; 
  maf_keep_overlapping = FALSE; 
  collapse_missing = FALSE;
  fourD = FALSE; 
  mark_missing_maxsize = -1; 
  missing_as_indels = FALSE;
  unmask = FALSE;
  cats_to_do = NULL; 
  aggregate_list = NULL;
  msa_fname_list = NULL; 
  order_list = NULL; 
  fill_N_list = NULL;
  map = NULL;
  gff = NULL;
  cm = NULL;

}

void readAlignment(char** filename, char** format, int* error, int* address, int* numberSpecies, int* length, char** alphabet, char** species){
  MSA* align;
  int i;

  input_format = msa_str_to_format(format[0]);
  align = msa_new_from_file(fopen_fname(filename[0], "r"), input_format, NULL);
  if (align == NULL) {
    *error=-1;
    return;
  }
  *address=(unsigned int)align;
  *numberSpecies=align->nseqs;
  *length=align->length;
  *alphabet=align->alphabet;
  for(i=0; i<*numberSpecies; i++){
    species[i]=align->names[i];
  }
}

void extractSS(int *address, int *storeOrder) {
  MSA *msa = (MSA*)(*address);
  ss_from_msas(msa, 1, *storeOrder, NULL, NULL, NULL, -1);
}

void writeAlignment(char** filename, char** format, int* error, int* address){
  FILE* outfile=fopen_fname(filename[0],"w");
  if (outfile==NULL){
    *error=-1;
    return;
  }
  output_format=msa_str_to_format(format[0]);
  msa_print(outfile, (MSA*)address[0], output_format, FALSE);
  fclose(outfile);
}

void freeAlignment(int* address){ 
  msa_free((MSA*)address[0]);
}


