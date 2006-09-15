/* $Id: findMafDiffs.c,v 1.1 2006-09-15 20:16:33 amit Exp $ */
/* $Log: not supported by cvs2svn $ */
/* this program compares the human, chimp, rhesus alleles at a given chromosome position in a maf file
   and prints out a message stating where bp differences occur. It takes 3 arguments:
   argc[1]: a maf file
   argc[2]: a reference sequence
   argc[3]: a space delimited  coordinate file of the form <chr> <start> <end> giving the ranges of the 
            of the reference sequence to get bp information for. 
   Output is written to STDOUT.
   Sample usage: findMafDiffs  chr22.maf /usr/data/hg18/chromFa/chr22.fa  test.txt
   Sample out: 
*/
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <getopt.h>
#include <misc.h>
#include <msa.h>
#include <maf.h>
#include <gff.h>
#include <sufficient_stats.h>

#define HUMAN 0
#define CHIMP 1
#define RHESUS 2

void getAlleles ( MSA *msa, int position, char *alleles);


int main (int argc, char *argv[]) {


  FILE *MSAF, *REFSEQF, *COORDF;
  MSA *msa;

  int i;
  int j;
  char chrom[16]; 
  char inputline[3000];
  int start, end, range;

  /* open the msa, refseq, coordinate  file handles */
  MSAF = fopen_fname(argv[1], "r"); 
  REFSEQF = fopen_fname(argv[2], "r");
  COORDF = fopen_fname(argv[3], "r");

  /* read the maf into memory */
    printf("reading %s file ....\n", argv[1]);
    msa = maf_read(MSAF, REFSEQF, 1, NULL, NULL, NULL, -1, 1, NULL, 1, FALSE);

    printf("alphabet %s\n", msa->alphabet);
    printf("src size %d\n", msa->length);

  /* array of char to hold alleles at a particular bp 
     this is the array passed to getAlleles function */
    char alleles[msa->nseqs];
  
    while (fgets(inputline, 3000, COORDF) ) { /* read the whole line */
      sscanf(inputline,"%s %d %d\n", chrom, &start, &end); /* now parse the line into its seperate fields */
      range = (end - start) + 1;
      printf("%s:%d..%d range:%d\n", chrom, start, end, range);
      
      /* get the alleles in the chromosome range between start .. end */

      for ( j= start; j < end; j++ ) { 
	getAlleles(msa, j, alleles);

	/* compare the human and chimp alleles */
	if ( (alleles[HUMAN] != alleles[CHIMP]) && alleles[CHIMP] != '*' && alleles[CHIMP] != 'N' && alleles[CHIMP] != '-' ) 
	  printf("%s %c differs from %s %c at %d\n", msa->names[HUMAN], alleles[HUMAN], msa->names[CHIMP], alleles[CHIMP], j);

	/*if ( alleles[HUMAN] != alleles[RHESUS] && alleles[RHESUS] != '*' ) 
	  printf("%s %c differs from %s %c at %d\n", msa->names[HUMAN], alleles[HUMAN], msa->names[RHESUS], alleles[RHESUS], j);
	*/
      }
      printf("\n");
    }
    
}


/* return a pointer to *char containing the chracters at the specified position in the maf */
void  getAlleles ( MSA *msa, int position, char *alleles) {

  int i;

  for (i=0; i < msa->nseqs; i++) 
    *(alleles+i) = ss_get_char_pos(msa, position,i,0);

}
