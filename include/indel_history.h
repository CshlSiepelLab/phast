/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file indel_history.h
    Functions and structs to hold and interrogate the history of insertions / deletions
    @ingroup phylo
*/

#ifndef IND_HIST
#define IND_HIST

#include <stdio.h>
#include <trees.h>
#include <msa.h>

/** Number of Indel types i.e. insert, delete, none */
#define NINDEL_CHARS 3
/** What type of change, insertion, deletion, none */
typedef enum {INS, /**< Insertion */
	      DEL, /**< Deletion */
	      BASE /**< No insertion or deletion */
		} indel_char; /* note: order is used in places */

/** Information about an indel */
typedef struct {
  indel_char type;  /**< What type of change, insertion, deletion, none */
  int start;  	    /**< Site the indel starts at */
  int len;	    /**< Length of indel */
} Indel;

/** Compact Indel history of a Multiple Sequence Alignment */
typedef struct {
  TreeNode *tree; /**< Tree describing structure of data */
  int ncols;   /**< Number of sites from data */
  List **indels; /**< List of each indel */
} CompactIndelHistory;

/** Indel history of a Multiple Sequence Alignment */
typedef struct {
  TreeNode *tree; /**< Tree describing structure of data */
  int ncols;      /**< Number of sites from data */
  char **indel_strings;    /**< list of strings describing each indel */        /* make bin vector later */
} IndelHistory;

/** \name Indel History allocation functions 
 \{ */

/** Create an Indel History object for a dataset 
  @param tree Representing dataset structure
  @param ncols Number of columns in dataset
  @result New indel history object
*/
IndelHistory *ih_new(TreeNode *tree, int ncols);

/** Create a compact Indel history object for a dataset 
  @param tree Representing dataset structure
  @param ncols Number of columns in dataset
  @result New compact indel history object
*/
CompactIndelHistory *ih_new_compact(TreeNode *tree, int ncols);

/** Given an alignment and tree structure, extract indels
   @param msa Multiple Sequence Alignment sequence data
   @param tree Tree structure
   @result New Indel history object from sequence and tree data
   @note Includes sequences for ancestral nodes and leaves
*/
IndelHistory *ih_extract_from_alignment(MSA *msa, TreeNode *tree);

/**  Reconstruct an indel history by parsimony from an alignment, given a tree 
  @param msa Multiple Sequence Alignment sequence data
  @param tree Tree structure
  @result Indel history object from sequence and tree data
*/
IndelHistory *ih_reconstruct(MSA *msa, TreeNode *tree);

/** \} \name Indel History cleanup functions 
 \{ */

/** Free an Indel History object
 @param ih Indel History object to free */
void ih_free(IndelHistory *ih);

/** Free a Compact Indel History object
 @param cih Compact Indel History object to free
 */
void ih_free_compact(CompactIndelHistory *cih);

/** \} \name Indel History convert between compact and normal functions 
 \{ */

/** Expand a compact indel history.
  @param cih Compact indel history
  @result Indel history
*/
IndelHistory *ih_expand(CompactIndelHistory *cih);

/** Compact an indel history
  @param ih Indel history
  @result Compact indel history
*/
CompactIndelHistory *ih_compact(IndelHistory *ih);

/** \} \name Indel History read/write file access functions 
 \{ */

/** Save an indel history to a file
    @param ih Indel History to write to file
    @param outf File descriptor to write Indel History into
    @param msa_name Name of the Multiple Sequence Alignment
    @param prog_name Name of program that generated indel history
    @warning Make sure ancestral nodes have been labeled before calling
  */
void ih_print(IndelHistory *ih, FILE *outf, char *msa_name, char *prog_name);

/** Save a compact indel history to a file
    @param cih Compact Indel History to write to file
    @param outf File descriptor to write compact Indel History into
    @param msa_name Name of the Multiple Sequence Alignment
    @param prog_name Name of program that generated indel history
    @warning Make sure ancestral nodes have been labeled before calling
*/
void ih_print_compact(CompactIndelHistory *cih, FILE *outf, char *msa_name, 
                      char *prog_name);


/** Read a compact indel history object from a file.
   @param inf Input file to read compact indel history from
   @result Compact indel history data read from file
*/
CompactIndelHistory *ih_read_compact(FILE *inf);

/** Read a indel history from a file.
    @param inf Input file to read indel history from
    @result Indel history data read from file
*/
IndelHistory *ih_new_from_file(FILE* inf);

/** \} \name Indel History Misc. functions 
 \{ */

/**  Convert an indel history into an alignment.
   Alignment includes sequences for ancestral nodes as
   well as leaf nodes, and with '^' characters in place of '-' for
   insertions and '.' characters in place of '-' for deletions.
   @param ih Indel History to use as template to create modified MSA
   @param msa Multiple Sequence Alignment 
   @result MSA modified at sites where indels occurred
*/
MSA *ih_as_alignment(IndelHistory *ih, MSA *msa);


/** Convert names in an alignment from the convention used by
   inferAncestors to the convention used in PHAST, based on a given
   tree.
   @param[in,out] msa Multiple Sequence Alignment sequence data
   @param[in] tree Tree structure
 */
void ih_convert_ia_names(MSA *msa, TreeNode *tree);

#endif
