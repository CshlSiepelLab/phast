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

   @ingroup motif
*/ 

#ifndef TFBS_H
#define TFBS_H
#include <ms.h>

/** Convert a double to a charstring
    @param num Double to convert to charstring
    @result Charstring containing double
*/
char *dtoa(double num);

/** Create a new MS object to hold sequences that are not nececssairly aligned
    @param seqs 2D array, one entry per sequence containing the bases.
    @param names 2D array, one entry per sequence containing name
    @param nseqs number of sequences in this MS
    @param lengths 2D array, one entry per sequence containing length of sequence
    @param alphabet characters thare are valid bases
    @param rangeLow Lower bound of %GC content this bin will hold (between 0 and 1)
    @param rangeHigh Upper bound of %GC content this bin will hold (between 0 and 1) 
    @result New GC bin to hold sequences
*/
MS *ms_new(char **seqs, char **names, int nseqs, const char *alphabet, double rangeLow, double rangeHigh);

/** Free a MS object 
    @param toFree the MS object to free
*/
void ms_free(MS *toFree);

/** Read a MEME text format file and extract Position Weight Matrix (PWM)
    @param filename Full Path to MEME formated file containing at least one PWM
    @result List of PWMs found in file
*/
List *pwm_read(const char *filename);

/** Read in non-aligned Multiple Sequences.
    @param filename Full path to file of type FASTA
    @param alphabet Characters allowed in sequences (normally 'ACGT')
    @result List of MSA objects with one sequence per MSA.
*/
MS *ms_read(const char *filename, const char *alphabet);

/** Helper function for building Markov Matrix of order N 
    @param inputMS MS object containing sequences used to make Markov Model
    @param norder Specifies order of the markov MATRIX to construct
    @param pseudoCount Number of counts to add for each possible entry (smoothing parameter)
    @param considerReverse If TRUE also take into account reverse complement frequencies
    @result Markov Matrix of order 'norder'
 */
Matrix *mm_build_helper(MS *inputMS, int norder, int pseudoCount, int considerReverse);


/** Reverse complement a DNA PWM
    @param m PWM to reverse complement
    @result Reverse complement of m
*/
Matrix *mat_reverse_complement(Matrix *m);


/** Calculate Markov Models for all groups 
    @param inputMS Multiple Sequences object containing sequences used to make Markov Model
    @param norderP Order of the Markov MODEL to construct
    @param pseudoCount Number of counts to add for each possible entry (smoothing parameter)
    @param considerReverse If TRUE also take into account reverse complement frequencies
	@result Markov Model (List of Markov Matrices)
	@note Each Markov Model contains norder Markov Matrices i.e. norder=3 matrices of norder 1, 2, 3
*/
List *mm_build(MS *inputMS, int norder, int pseudoCount, int considerReverse);


/** Identifies column of a PWM or Markov Matrix corrisponding to a single base.
    @param base Single character i.e. 'A' or 'T'
    @result Integer 0 to 3 specifying column that contains data for that base
*/
int basetocol(char base);

/** Lookup row of Markov Matrix corrisponding to bases i.e 'AT' == 3
    @param norder Order of Markov Matrix represented
    @param alph_size Alphabet size of data represented in Markov Matrix
    @param Integer specifying row of Markov Matrix mapped to given bases
*/
int basesToRow(int *previousBases, int norder, int alph_size);

/** Prints MS to file in FASTA format.  
   @param F File pointer where to output FASTA file
   @param ms MS object containing sequences & names to write out
*/
void ms_print(FILE *F, MS *ms);

void ms_print_fasta(FILE *F, MS *ms);

/** Prints MS to file in FASTA format.  (wrapper for ms_print)
   @param filename Full Path to where file should be written
   @param ms MS object containing sequences & names to write out
*/
void ms_print_to_file(const char *filename, MS *ms);

/** Compute GC content of each sequence in an ms
    @param ms MS object
    @return A vector of the same length as the number of sequences in ms, 
    containing the GC fraction of each sequence.
 */
Vector *ms_gc_content(MS *ms);


/** Scores sequences for matches to motif represented as PWM.
    @param seqName Name of the sequence being scored
    @param seqData Sequence data (bases) of the sequence being scored
    @param seqLen Length of the sequence being scored
    @param seqIdxOff Index offset of the sequence being scored
    @param seqAlphLen Length of the alphabet for the sequence being scored
    @param MarkovMatrices Markov Model (list of Markov Matrices) of GC content group sequence belongs to
    @param pwm Position Weight Matrix representing a Motif to scan for
    @param reverseCmpPwm Reverse complemented version of pwm parameter
    @param conservative If == 1 and encounters an 'N' base, the site gets a -Inf score
    @param threshold Score threshold that any score must be above to be part of the returned feature set
    @param strand Which strands to score and which results to return ("best", "both", "+", "-") more info in R docs
    @result List of scores as a Feature Set
*/
GFF_Set *ms_score(char *seqName, char *seqData, int seqLen, int seqIdxOff, int seqAlphLen, List *MarkovMatrices, Matrix *pwm, Matrix *reverseCmpPWM, int conservative, double threshold, char *strand); 
/** Simulate a sequence given a Markov Model
    @param mm Markov Model containing probabilities used to generate sequence
    @param norder Order of Markov Model mm
    @param alph_size Alphabet size of data used to generate mm
    @param length Length of sequence to generate
    @result Char string containing generated sequence i.e. "ATAGCGC" for length 7
 */
char *ms_simulate(List *mm, int norder, int alph_size, int length);
 
#endif
