/* $Id: pasteCds.c,v 1.1 2006-05-28 15:54:13 acs Exp $:
   Amit Indap indapa@gmail.com
   Bustamante Lab C Source File
   Copyright 2006, Cornell University


   Description: This program pastes together CDS portions of clean_genes
   output. It takes 3 arguments: a ss file containing the MSA, a gff file
   which contains the output of a clean_genes process and the basename for
   the files which output is written to.

   Output is written in phylip format. Each "transcript_id" i.e. NM_ Refseq
   accession is written to its own *.phy format. 

   This program was written with help of  Adam Siepel and makes extensive use
   of his PHAST library. 

   Usage: pasteCds [OPTIONS] <ss> <gff> <outroot>

   History:
   $Log: not supported by cvs2svn $


*/


#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <misc.h>
#include <msa.h>
#include <gff.h>
#include <sufficient_stats.h>

void usage(char *prog) {
  printf("\n\
PROGRAM:      pasteCds\n\
DESCRIPTION:  \n\
USAGE:        pasteCds [OPTIONS] <ss> <gff> <outroot>\n\
OPTIONS:\n\
    --help, -h\n\
        Print this help message.\n\n");
  exit(0);
}

int main(int argc, char *argv[]) {
  FILE *SSF, *GFFF, *OUTF;
  MSA *msa;
  GFF_Set *gff;
  char c;
  int opt_idx;
  char *outroot = NULL;
  int i, j;
  char outname[50];
   List *l = lst_new_ptr(1);
  /*  List *l = lst_new_ptr(2);*/

  struct option long_opts[] = {
    {"a_long", 1, 0, 'a'},
    {"b_long", 1, 0, 'b'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "a:b:h", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'a':
      break;
    case 'b':
      break;
    case 'h':
      usage(argv[0]);
    case '?':
      die("Bad argument.  Try '%s -h'.\n", argv[0]);
    }
  }

  if (optind != argc - 3) 
    die("Three arguments required.  Try '%s -h'.\n", argv[0]);
    
  SSF = fopen_fname(argv[optind], "r");
  GFFF = fopen_fname(argv[optind+1], "r");
  outroot = argv[optind+2];

  fprintf(stderr, "Reading alignment...\n");
  msa = msa_new_from_file(SSF, SS, NULL);

  fprintf(stderr, "Reading GFF...\n");
  gff = gff_read_set(GFFF);

  fprintf(stderr, "Mapping coordinates...\n");
  msa_map_gff_coords(msa, gff, 1, 0, 0, NULL);

  lst_push_ptr(l, str_new_charstr("CDS"));
  /*lst_push_ptr(l, str_new_charstr("stop_codon"));*/
  gff_filter_by_type(gff, l, FALSE, NULL);
  gff_group(gff, "transcript_id");

  /* iterate throught the gff groups */ 
  for (i = 0; i < lst_size(gff->groups); i++) {
    MSA *gene_msa = NULL;
    GFF_FeatureGroup *group = lst_get_ptr(gff->groups, i);
    char strand = '.';
    fprintf(stderr, "Processing group '%s'...\n", group->name->chars);
    for (j = 0; j < lst_size(group->features); j++) {
      GFF_Feature *f = lst_get_ptr(group->features, j); 
      MSA *subaln = msa_sub_alignment(msa, NULL, FALSE, f->start-1, f->end);
      ss_to_msa(subaln);
      ss_free(subaln->ss);
      subaln->ss = NULL;
     
      if (j == 0)               /* initialize gene_msa */
        gene_msa = subaln;
      else {
        msa_concatenate(gene_msa, subaln);
        msa_free(subaln); 
      }
      if (strand == '.') strand = f->strand;
      else if (strand != f->strand) 
        die("ERROR: features in group don't have same strand.\n");
    }

    if (strand == '-')
      msa_reverse_compl(gene_msa);

    sprintf(outname, "%s.%s.ph", outroot, group->name->chars);
    OUTF = fopen_fname(outname, "w+");
    msa_print(OUTF, gene_msa, PHYLIP, FALSE);
    msa_free(gene_msa);
  }

  return 0;
}
