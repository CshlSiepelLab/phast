/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file maf_block.h
    Reading of blocks from MAF ("Multiple Alignment Format")
    files, without generating Sufficient Statistics.
    (See http://www.bx.psu.edu/miller_lab.)  These functions read
    a MAF file block by block, and perform operations on the block
    without the need to convert to sufficient statistic format.
    The MAF_BLOCK format is not currently supported by most phast
    programs, except maf_read.c.
    @ingroup msa
*/

#ifndef MAF_BLOCK_H
#define MAF_BLOCK_H

#include "stdio.h"
#include "msa.h"
#include "hashtable.h"

/** Holds per-species data for a Maf Block */
typedef struct {
  String *seq;  /**< Sequence data of sub block */
  String *src,   /**< Source of sequence data */
    *specName;  /**< Part of src before the first '.' */
  long start;   /**< Starting column */
  int size;    /**< Length of the block */
  char strand; /**< Type of strand of the sequence*/
  long srcSize; /**< Size of the source */
  int numLine;  /**< Number of lines corresponding to this 
			species in this block. */
  char lineType[4];  /**< type of line i, either 's', 'q', 'i', 'e' */
  char iStatus[2];  /**< leftStatus and rightStatus; only defined if
                       one of lineType[0..(numLine-1)] is 'i' */
  int iCount[2];    /**< leftCount and rightCount for i-line */
  char eStatus;     /**< ??? */
  String *quality;  /**< Array of quality scores of length seqlen.  Only
                      defined if one of lineType is 'q' */
} MafSubBlock;

/** Holds data about a Maf Block */
typedef struct MAFBLOCK {
  String *aLine;
  Hashtable *specMap;  /**< Hash that shows which element of data refers to
                         a given species */
  int seqlen;   /**< Length of aligned sequences */
  List *data;  /**< List of pointers to type MafSubBlock.  One for each species */
  struct MAFBLOCK *prev, *next; /**< Pointers to other Maf blocks in a MAF file */
} MafBlock;

/** \name MAF block read/write file functions 
 \{ */

/** Save a Maf Block to a file.
    @param outfile File to save Maf Block to
    @param block Maf Block to save
    @param pretty_print Use identical sign '.' when matches refseq
*/
void mafBlock_print(FILE *outfile, MafBlock *block, int pretty_print);

/**  Read next block from Maf File.
     @pre If optional parameter specHash is not NULL make sure it is initialized
     @param mfile MAF file to read next MAF block from
     @param specHash  (Optional) Any new species encountered added to this hash
     @param numspec   (Optional) Number of species added to specHash
     @result MafBlock read from MAF file, OR NULL if EOF
     @warning If you use specHash, you also must use numspec.
*/
MafBlock *mafBlock_read_next(FILE *mfile, Hashtable *specHash, int *numspec);

/** Opens new MAF file for writing and prints minimal header for MAF.  
    @param fn File to write MAF to, NULL == stdout
    @param argc Number of comment lines (in argv below) to write to MAF 
    @param argv Comment lines of MAF giving parameters this file
 */
FILE *mafBlock_open_outfile(char *fn, int argc, char *argv[]);

/** Finish writing MAF file; Prints #eof and closes file if outfile != stdout
    @param outfile File to close writing for
*/
void mafBlock_close_outfile(FILE *outfile);

/** \name MAF block cleanup functions 
 \{ */

/** Free data within a Maf Block but not the object itself 
    @param block Maf Block containing data to free
*/
void mafBlock_free_data(MafBlock *block);

/** Free entire Maf Block and its data.
    @param block Maf Block to free
*/
void mafBlock_free(MafBlock *block);

/** \} \name MAF block copy functions 
   \{ */

/** Create a copy of a Maf Sub Block.
    @param src Maf Sub Block to be copied
    @return New copy of Maf Sub Block
*/
MafSubBlock *mafSubBlock_copy(MafSubBlock *src);

/** Create a copy of a Maf Block.
    @param src Maf Block to be copied
    @result New copy of Maf Block
*/
MafBlock *mafBlock_copy(MafBlock *src);


/** \name MAF block modify functions 
 \{ */

/** Reorder rows of maf block given list with names of species in desired order.
    @param block Maf Block to re-order
    @param specNameOrder List of species names in the order desired
 */
void mafBlock_reorder(MafBlock *block, List *specNameOrder);

/** Remove or Keep only species specified in list of names from MafBlock.
    @param block Maf Block to remove species from
    @param specNameList List of names to either be removed or kept
    @param if include=0 remove species in specNameList, else keep only species in specNameList
 */
void mafBlock_subSpec(MafBlock *block, List *specNameList, int include);


/** Trims a Maf Block based on start and end columns.
    @param block Maf Block to trim
    @param statcol Starting column used to crop sequences
    @param endcol Ending column used to crop sequences, ignored if == -1
    @param refseq (Optional) If NULL use alignment column coords (1st column is 1), otherwise use refseq coords.
    @param offset Add an offset to all coords in block
    @result 0 if block is empty, otherwise 1
    @note (1-based, endcol inclusive).
*/
int mafBlock_trim(MafBlock *block, int startcol, int endcol,
		  String *refseq, int offset);

/** Sets block to the sub-alignment from alignment column start to end.
    @param block Maf Block to perform sub alignment on
    @param start Column in alignment defining the beginning of the sub block to align
    @param end Column in alignment defining the ending of the sub block to align
    @note (1-based, in the reference frame of entire alignment)
    @note Both start and end should be in the range [1,block->seqlen]. 
*/
void mafBlock_subAlign(MafBlock *block, int start, int end);

/** Remove any line in the Maf Block that is of type 'i'
    @param block Maf Block to remove 'i' type lines from
*/
void mafBlock_strip_iLines(MafBlock *block);

/** Remove any line in the Maf Block that is of type 'e'
    @param block Maf Block to remove 'e' type lines from
*/
void mafBlock_strip_eLines(MafBlock *block);

/** Remove any line in the Maf Block that is of type 'i' or 'e'
    @param block Maf Block to remove 'i' or 'e' type lines from
*/
void mafBlock_strip_ieLines(MafBlock *block);

/** Mask a region of the alignment.
    @param block MafBlock to mask
    @param mask_feats features defining regions of alignment to mask
    @param speclist A list of character vectors identifying species to mask
 */
void mafBlock_mask_region(MafBlock *block, GFF_Set *mask_feats, List *speclist);

/**  Threshold all bases based on quality score.
     Change all bases with quality score <= cutoff to N
     @param block Maf Block to threshold
     @param cutoff Threshold value; each base's quality score must be greater than cutoff, otherwise it is changed to 'N'
     @param outfile File to output coordinates of masked bases to, if not NULL.  Coordinates will be relative to refseq, and the final column gives the name of the species which was masked.
     @note If outfile is not null, note that masked bases which align to a gap in refseq may not be indicated in the outfile.
*/
void mafBlock_mask_bases(MafBlock *block, int cutoff, FILE *outfile);

/** \} \name MAF block get info functions 
   \{ */

/** Get species name of first species in block 
    @param block Maf Block to get first species from
    @result Species name of first species in block
*/
String *mafBlock_get_refSpec(MafBlock *block);

/** Get the starting index for a given species 
    @param block Maf Block containing the species
    @param specName Name of the species to get starting index of
    @result Starting index of species, OR -1 if species does not exist in given$
*/
long mafBlock_get_start(MafBlock *block, String *specName);

/** Get number of bases (non-gaps) in alignment for a given species. 
    @param block Maf Block containing species specName sequence data
    @param specName (Optional) Species Name to count bases for
    @result Number of bases (non-gaps) for a given species, OR total alignment $
*/
int mafBlock_get_size(MafBlock *block, String *specName);

/** Count number of species in a Maf Block 
    @param block Maf Block containing species to count
*/
int mafBlock_numSpec(MafBlock *block);

/** Test if block is entirely gaps.
    @param block Maf Block to test
    @result 1 if block is entirely gaps, 0 otherwise
*/
int mafBlock_all_gaps(MafBlock *block);


/** \} */

#endif
