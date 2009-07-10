/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2009 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: maf_parse.c,v 1.42 2009-01-09 22:01:00 mt269 Exp $ */

#include <stdlib.h>
#include <stdio.h>
#include <msa.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include <string.h>
#include <ctype.h>
#include <sufficient_stats.h>
#include <local_alignment.h>
#include <maf.h>
#include <maf_block.h>

void print_usage() {
    printf("\n\
USAGE: maf_parse [OPTIONS] <infile>\n\
\n\
DESCRIPTION:\n\
\n\
    Reads a MAF file and perform various operations on it.\n\
    Performs parsing operations block-by-block whenever possible,\n\
    rather than storing entire alignment in memory.  \n\
    Can extract a sub-alignment from an alignment (by row\n \
    or by column).  Can extract features given GFF, BED, or \n\
    genepred file.  Can also extract sub-features such as CDS1,2,3\n\
    or 4d sites.  Can perform various functions such as gap\n\
    stripping or re-ordering of sequences.  Capable of reading and\n\
    writing in a few common formats, but will not load input or output\n\
    alignments into memory if output format is MAF.\n\
\n\
OPTIONS:\n\
\n\
 (Obtaining sub-alignments and re-ordering rows)\n\
    --start, -s <start_col>\n\
        Starting column of sub-alignment (indexing starts with 1).\n\
        Default is 1.\n\
\n\
    --end, -e <end_col>\n\
        Ending column of sub-alignment.  Default is length of\n\
        alignment.\n\
\n\
    --seqs, -l <seq_list>\n\
        Comma-separated list of sequences to include (default)\n\
        exclude (if --exclude).  Indicate by sequence number or name\n\
        (numbering starts with 1 and is evaluated *after* --order is\n\
        applied).\n\
\n\
    --exclude, -x\n\
        Exclude rather than include specified sequences.\n\
\n\
    --order, -O <name_list>\n\
        Change order of rows in alignment to match sequence names\n\
        specified in name_list.  The first name in the alignment becomes\n\
        the reference sequence.\n\
\n\
    --no-refseq, -n\n\
        Do not assume first sequence in MAF is refseq.  Instead, use coordinates\n\
        given by absolute position in alignment (starting from 1).\n\
\n\
(Splitting into multiple MAFs by length)\n\
    --split, -S length \n\
        Split MAF into pieces by length, and puts output in outRootX.maf, where\n\
        X=1,2,...,numPieces.  outRoot can be modified with --out-root, and the\n\
        minimum number of digits in X can be modified with --out-root-digits.\n\
        Splits between blocks, so that each output file does not exceed specified\n\
        length.  By default, length is counted by distance spanned in alignment\n\
        by refseq, unless --no-refseq is specified.\n\
\n\
   --out-root, -r <name>\n\
        Filename root for output files produced by --split (default \"maf_parse\").\n\
\n\
   --out-root-digits, -d <numdigits>\n\
        (for use with --split).  The minimum number of digits used to index each\n\
        output file produced by split.\n\
\n\
 (Other)\n\
    --help, -h\n\
        Print this help message.\n\n");
}

int main(int argc, char* argv[]) {
  char *maf_fname = NULL, *out_root_fname = "maf_parse";
  String *refseq = NULL, *currRefseq;
  int opt_idx, startcol = 1, endcol = -1, include = 1, splitInterval = -1;
  char c, outfilename[1000], splitFormat[100]="%s%.1i.maf";
  List *order_list = NULL, *seqlist_str = NULL;
  MafBlock *block;
  FILE *mfile, *outfile=NULL;
  int useRefseq=TRUE, currLen=-1, blockIdx=0, currSize, sortWarned=0;
  int lastIdx = 0, currStart=0;

  struct option long_opts[] = {
    {"start", 1, 0, 's'},
    {"end", 1, 0, 'e'},
    {"seqs", 1, 0, 'l'},
    {"exclude", 0, 0, 'x'},
    {"order", 1, 0, 'O'},
    {"split", 1, 0, 'S'},
    {"out-root", 1, 0, 'r'},
    {"out-root-digits", 1, 0, 'd'},
    {"no-refseq", 0, 0, 'n'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "s:e:l:O:r:S:d:nxh", long_opts, &opt_idx)) != -1) {
    switch(c) {
    case 's':
      startcol = get_arg_int(optarg);
      break;
    case 'e':
      endcol = get_arg_int(optarg);
      break;
    case 'l':
      seqlist_str = get_arg_list(optarg);
      break;
    case 'O':
      order_list = get_arg_list(optarg);
      break;
    case 'x':
      include = FALSE;
      break;
    case 'S':
      splitInterval = atoi(optarg);
      break;
    case 'r':
      out_root_fname = optarg;
      break;
    case 'd':
      sprintf(splitFormat, "%%s%%.%si.maf", optarg);
      break;
    case 'n':
      useRefseq = FALSE;
      break;
    case 'h':
      print_usage();
      exit(0);
    case '?':
      fprintf(stderr, "Bad argument.  Try 'msa_view -h' for help.\n");
      exit(1); 
    }
  }

  if (optind >= argc) 
    die("Missing alignment filename.  Try 'msa_view -h' for help.\n");
  else if (optind == argc - 1) 
    maf_fname = argv[optind];
  else 
    die("ERROR: Too many arguments.  Try 'msa_view -h' for help.\n");

  if (startcol < 1 || (endcol != -1 && endcol < startcol))
    die("ERROR: must have 1 <= start <= end <= [msa_length]\n");

  /* Check to see if --do-cats names a feature which is length 1. 
     If so, set output_format to SS ? or FASTA ? */
  
  mfile = fopen_fname(maf_fname, "r");
  block = mafBlock_read_next(mfile, NULL, NULL);

  if (splitInterval == -1) {
    //do we want to copy header from original MAF in this case?
    mafBlock_open_file(NULL);
  }

  while (block != NULL) {
    if (order_list != NULL)
      mafBlock_reorder(block, order_list);
    if (seqlist_str != NULL)
      mafBlock_subSpec(block, seqlist_str, include);
    if (mafBlock_numSpec(block)==0 || mafBlock_all_gaps(block)) 
      goto get_next_block;

    if (useRefseq) {  //get refseq and check that it is consistent in MAF file
      currRefseq = mafBlock_get_refSpec(block);
      if (refseq == NULL) 
	refseq = str_new_charstr(currRefseq->chars);
      else if (str_compare(refseq, currRefseq)!=0)
	die("Error: refseq not consistent in MAF (got %s, %s)\n",
	    refseq->chars, currRefseq->chars);
    }
    
    if (startcol != 1 || endcol != -1) 
      if (0 == mafBlock_trim(block, startcol, endcol, refseq, useRefseq ? 0 : lastIdx))
	goto get_next_block;

    currSize = mafBlock_get_size(block, refseq);
    if (useRefseq) {
      currStart = mafBlock_get_start(block, refseq);
      if (currStart < lastIdx && sortWarned == 0) {
	fprintf(stderr, "Warning: input MAF not sorted with respect to refseq.  Output files may not represent contiguous alignments. (%i, %i)\n", lastIdx, currStart);
	sortWarned = 1;
      }
    }
    else currStart = lastIdx;
    
    lastIdx = currStart + currSize;

    //split by length
    if (splitInterval != -1) {
      if (currLen == -1 || currLen+currSize > splitInterval) {
	sprintf(outfilename, splitFormat, out_root_fname, ++blockIdx);
	if (outfile != NULL) {
	  mafBlock_close_file(outfile);
	}
	outfile = mafBlock_open_file(outfilename);
	currLen = 0;
      }
      currLen += currSize;
    }
    else outfile = stdout;

    mafBlock_print(outfile, block);

  get_next_block:
    mafBlock_free(block);
    block = mafBlock_read_next(mfile, NULL, NULL);
  }
  mafBlock_close_file(outfile);

  fclose(mfile);


  return 0;
}
