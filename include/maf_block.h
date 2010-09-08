/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: maf_block.h,v 1.7 2009-01-09 22:01:00 mt269 Exp $ */

/** \file maf_block.h 
    Reading of alignments from MAF ("Multiple Alignment Format")
    files, as produced by MULTIZ and TBA.
    (See http://www.bx.psu.edu/miller_lab.)  These functions read
    a MAF file block by block, and perform operations on the block
    without the need to convert to sufficient statistic format.
    The MAF_BLOCK format is not currently supported by most phast
    programs, except maf_read.c.
    \ingroup msa
*/

#ifndef MAF_BLOCK_H
#define MAF_BLOCK_H

#include "stdio.h"
#include "msa.h"
#include "hashtable.h"

typedef struct {
  String *seq;
  String *src, *specName;  //specName is part of src before first '.'
  int start;
  int size;
  char strand;
  int srcSize;
  int numLine;     //number of lines corresponding to this 
                   //species in this block.
  char lineType[4];  //type of line i, either 's', 'q', 'i', 'e'
  char iStatus[2];  //leftStatus and rightStatus; only defined if
                    //one of lineType[0..(numLine-1)] is 'i'
  int iCount[2];    //leftCount and rightCount for i-line
  char eStatus;
  String *quality;  //array of quality scores of length seqlen.  Only
                  //defined if one of lineType is 'q'
} MafSubBlock;

typedef struct MAFBLOCK {
  String *aLine;
  Hashtable *specMap;  //hash that shows which element of data refers to
                       //a given species
  int seqlen;
  List *data;  //List of pointers to type MafSubBlock.  One for each species
  struct MAFBLOCK *prev, *next;
} MafBlock;
  
             
void mafBlock_print(FILE *outfile, MafBlock *block, int pretty_print);
MafBlock *mafBlock_read_next(FILE *mfile, Hashtable *specHash, int *numspec);

MafSubBlock *mafSubBlock_copy(MafSubBlock *src);
MafBlock *mafBlock_copy(MafBlock *src);

//reorder rows of maf block given list with names of species in desired order
void mafBlock_reorder(MafBlock *block, List *specNameOrder);

//removes species not in specNameList if include=0,
//otherwise removes species in specNameList
void mafBlock_subSpec(MafBlock *block, List *specNameList, int include);

//opens new MAF file for writing and prints minimal header for MAF.  
//If fn==NULL uses stdout.  If argc > 0 then it uses this to 
//print second comment line of MAF giving parameters
//this file
FILE *mafBlock_open_outfile(char *fn, int argc, char *argv[]);

//prints #eof and closes file if outfile != stdout
void mafBlock_close_outfile(FILE *outfile);

//returns species name of first species in block
String *mafBlock_get_refSpec(MafBlock *block);

//returns the start idx for species given
//returns -1 if specName not in block
int mafBlock_get_start(MafBlock *block, String *specName);

//returns the number of bases (non-gaps) in alignment for species
//given.  Returns total alignment length if specName == NULL.  Returns
//-1 if specName not in block
int mafBlock_get_size(MafBlock *block, String *specName);

void mafBlock_free_data(MafBlock *block);

void mafBlock_free(MafBlock *block);

int mafBlock_numSpec(MafBlock *block);

//returns 1 if block is entirely gaps, 0 otherwise
int mafBlock_all_gaps(MafBlock *block);

//returns # of species with data for block (including e-lines)
int mafBlock_numSpec(MafBlock *block);

//trims the block (if necessary) so that it only contains columns
//from startcol to endcol (1-based, endcol inclusive).  If refseq
//is NULL use alignment column coordinates (1st column is 1).  
//Otherwise use refseq coordinates.  Add offset to all coordinates in
//block.  Returns 0 if resulting block is empty, otherwise returns 1.
//If endcol is -1, it is ignored (all columns >= startcol are included).
int mafBlock_trim(MafBlock *block, int startcol, int endcol,
		  String *refseq, int offset);

//sets block to the sub-alignment from alignment column start to end 
//(1-based, in the reference frame of entire alignment).  start should be
//<= end and both should be in the range [1,block->seqlen]. 
void mafBlock_subAlign(MafBlock *block, int start, int end);

//strip iLines
void mafBlock_strip_iLines(MafBlock *block);

//strip eLines
void mafBlock_strip_eLines(MafBlock *block);

//strip i- and e-Lines
void mafBlock_strip_ieLines(MafBlock *block);

void mafBlock_mask_bases(MafBlock *block, int cutoff);

//void mafBlock_mask_indels(MafBlock *block, int cutoff, FILE *mfile);
#endif
