/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/


/** @file tfbs.h
   Transcription Factor binding sites functions.

   Identitifes Transcription Factor Binding Sites (TFBS) given a set of sequences
   that might contain TFBS and a set of Transcription Factors to search for.
   Transcription Factors are specified as MEME files containing Position Weight Matrices (PWMs)

   Steps to identify TFBS
    1. Read in PWM
          pwm = read_pwm(
    2. Read in multiple sequences (ms) file
	ms = ms_read("input.fas", "ACGT")
    3. (Optional) Split sequences based on windows
	ms =
    4. Group sequences based on % GC content
	groups = ms_group_by_gc(
    5. Build MarkovModels (one for each group)
	 models = mm_build(
    6. Calculate scores for sequences
	rScores = mm_comute_scores(
    7. Simulate sequences (one per group) based on Markov Models
	simGroups = mm_simulate_seq(
    8. Calculate Scores for simulated sequences
	sScores = mm_compute_scores(
    9. Identify threshold
    10. Output binding sites

   @ingroup motif
*/ 

#ifndef TFBS_H
#define TFBS_H

/** Read a MEME text format file and extract Position Weight Matrix (PWM)
    @param filename Full Path to MEME formated file containing at least one PWM
    @result List of PWMs found in file
*/
List *read_pwm(const char *filename);

/** Read in non-aligned Multiple Sequences.
    @param filename Full path to file of type FASTA
    @param alphabet Characters allowed in sequences (normally 'ACGT')
    @result List of MSA objects with one sequence per MSA.
*/
List *ms_read(const char *filename, const char *alphabet);


/** Helper function for building Markov Matrix of order N 
    @param group GC content group of sequences
    @param norder Specifies order of the markov matrix to construct
    @param lowerOrderMM Markov Matrix of order (norder-1), or NULL if order=0
    @result Markov Matrix of order 'norder'
 */
Matrix *mm_build_helper(ListOfLists *group, int norder, Matrix *lowerOrderMM);


/** Calculate Markov Models for all groups 
    @param msGroupListP List of GC groups containing sequences
    @param norderP Order of the Markov Model to construct
	@result List of Markov Models (1 per group)
	@note Each Markov Model contains norder Markov Matrices i.e. norder=3 matrices of norder 1, 2, 3
*/
List *mm_build(ListOfLists *msGroupListP, int norderP);


/** Identifies column of a PWM or Markov Matrix corrisponding to a single base.
    @param base Single character i.e. 'A' or 'T'
    @result Integer 0 to 3 specifying column that contains data for that base
*/
int basetocol(char base);

/** Lookup row of Markov Matrix corrisponding to bases i.e 'AT' == 3
    @param norderP Order of Markov Matrix represented
    @param alph_size Alphabet size of data represented in Markov Matrix
    @param Integer specifying row of Markov Matrix mapped to given bases
*/
int basesToRow(int *previousBases, int norderP, int alph_size);

/** Scores sequences for matches to motif represented as PWM.
    @param ms Single sequence to scan for matches to PWM
    @param mm Markov Model of GC content group sequence belongs to
    @param pwm Position Weight Matrix representing a Motif to scan for
    @param mmOrder Order of the Markov Model
    @result List of scores, one per base
*/
List *mm_compute_scores(MSA *ms, List *MarkovMatrices, Matrix *pwm, int mmOrder);
#endif
