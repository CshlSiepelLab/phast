/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/*****************************************************
   rph_tfbs.c
   The RPHAST handles to functions that identify
   Transcription factor binding sites

   Nick Peterson
*****************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <ctype.h>
#include <misc.h>
#include <sufficient_stats.h>
#include <local_alignment.h>
#include <trees.h>
#include <rph_util.h>
#include <Rdefines.h>
#include <R_ext/Random.h>
#include <tfbs.h>
#include <R.h>
#include <Rinternals.h>


//ncols and nrows have a name collision with Matrix->ncols, Matrix->nrows
#undef ncols
#undef nrows

/** Free an MS object 
    @param msP SEXP pointer to MS object to free
*/
void rph_ms_free(SEXP msP) {
  MS *ms;
  ms = (MS*)EXTPTR_PTR(msP);
  phast_unregister_protected(ms);
  ms_free(ms);
}

/** Create a new MS object
    @param seqsP SEXP pointer to list of strings, each string contains a sequence
    @param namesP SEXP pointer to list of strings, each string contains a name
    @param nseqsP SEXP pointer to Integer specifying how many sequences we have (probably unnecessary since Length(seqsP) would work)
    @param alphabetP SEXP pointer to string containing the valid non-missing characters that make up the sequences
    @param idxOffsets SEXP pointer to list of integer, each integer specifying the index offset of a sequence
    @return External pointer to newly created MS object in C.
*/
SEXP rph_ms_new(SEXP seqsP, SEXP namesP, SEXP nseqsP, 
                SEXP alphabetP, SEXP idxOffsetsP) {
  char **seqs=NULL, **names=NULL, *alphabet=NULL; 
  int nseqs=0, i, numProtect=0;
  MS *ms=NULL;
  SEXP result;

  nseqs = INTEGER_VALUE(nseqsP);


  if (namesP != R_NilValue) {
    PROTECT(namesP = AS_CHARACTER(namesP));
    numProtect++;
    names = smalloc(nseqs*sizeof(char*));
    for (i=0; i<nseqs; i++) {
      names[i] = smalloc((strlen(CHAR(STRING_ELT(namesP, i)))+1)*sizeof(char));
      strcpy(names[i], CHAR(STRING_ELT(namesP, i)));
    }
  }
  if (seqsP != R_NilValue) {
    PROTECT(seqsP = AS_CHARACTER(seqsP));
    numProtect++;
    seqs = smalloc(nseqs*sizeof(char*));
    for (i=0; i<nseqs; i++) {
      seqs[i] = smalloc((strlen(CHAR(STRING_ELT(seqsP, i)))+1)*sizeof(char));
      strcpy(seqs[i], CHAR(STRING_ELT(seqsP, i)));
    }
  }
  if (alphabetP != R_NilValue) {
    PROTECT(alphabetP = AS_CHARACTER(alphabetP));
    numProtect++;
    alphabet = smalloc((strlen(CHARACTER_VALUE(alphabetP))+1)*sizeof(char));
    strcpy(alphabet, CHARACTER_VALUE(alphabetP));
  }

  ms = ms_new(seqs, names, nseqs, alphabet, 0, 1);
  if (idxOffsetsP != R_NilValue) {
    PROTECT(idxOffsetsP = AS_INTEGER(idxOffsetsP));
    numProtect++;
    ms->idx_offsets = smalloc(nseqs*sizeof(int));
    for (i=0; i<nseqs; i++) {
      ms->idx_offsets[i] = INTEGER_POINTER(idxOffsetsP)[i];
    }
  }

  //don't free seqs, names, or alphabet because they are used by MS object
    
  PROTECT(result = rph_ms_new_extptr(ms));
  numProtect++;

  UNPROTECT(numProtect);
  return result;
}

/** Make an external pointer (to return to R) from an MS object stored in C.
    @param ms MS object stored in C to return by reference to R
    @result External pointer to MS object
*/
SEXP rph_ms_new_extptr(MS *ms) {
  SEXP result;
  ms_register_protect(ms);
  PROTECT(result=R_MakeExternalPtr((void*)ms, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(result, rph_ms_free, 1);
  UNPROTECT(1);
  return result;
}

/** Makes an external pointer (to return to R) from an MS object stored in C and optinally the GC range it represents.
    @param ms MS object in C
    @param includeGcRange Whether to include the GC range that the MS represents
    @return External Pointer to MS object in C
*/
SEXP group_to_SEXP(MS *ms, int includeGcRange)
{
  ListOfLists  *lolResult, *lolSequences;
  double *range = smalloc(2*sizeof(double));

  lolResult = lol_new(2);
  lolSequences = lol_new(1);
  //As long as the MS is not empty
  if (ms->nseqs != 0)
	{
	  lol_push_ms_ptr(lolSequences, ms, "");

      if(includeGcRange) {
        lol_push_lol(lolResult, lolSequences, "seqs");
        range[0] = ms->rangeLow;
        range[1] = ms->rangeHigh;
        lol_push_dbl(lolResult, range, 2, "GcRange");
        return rph_listOfLists_to_SEXP(lolResult);
      }  
	}
  return rph_listOfLists_to_SEXP(lolSequences);
}	

/** Transforms a matrix in C into a Matrix in R
    @param mat Matrix in C
    @result Matrix in R (SEXP)
*/
SEXP Matrix_to_SEXP(Matrix *mat)
{
  SEXP result;
  int col, row;

  PROTECT(result = allocMatrix(REALSXP, mat->nrows, mat->ncols));

  for (row = 0; row < mat->nrows; row++) {
    for (col = 0; col < mat->ncols; col++) {
      REAL(result)[col * mat->nrows + row] = mat_get(mat, row, col);
    }
  }
	
  UNPROTECT(1);
  return result;
}

/** Read multiple sequences from file.
    @param[in] filenameP File path to Multiple Sequence file
    @param[in] alphabetP Alphabet of characters used in sequences
    @result External pointer to MS object 
*/
SEXP rph_ms_read(SEXP filenameP, SEXP alphabetP)
{
  MS *sequences;

  if (filenameP == R_NilValue)
    die("ERROR: No sequences filename was provided");

  if (alphabetP == R_NilValue)
    sequences = ms_read(CHARACTER_VALUE(filenameP), NULL);
  else
    sequences = ms_read(CHARACTER_VALUE(filenameP), CHARACTER_VALUE(alphabetP));

  return group_to_SEXP(sequences, FALSE); 
}


/** Given an externalPointer (SEXP) to MS object in C, return a normal C pointer to the MS object.
    @param msP SEXP pointer to MS object in C
    @param Pointer to MS object in C
*/
MS *SEXP_to_group(SEXP msP) //Returns a single group of sequences only NOT ENTIRE LIST OF GROUPS
{
  MS *ms = (MS*)EXTPTR_PTR(msP);
  return ms;
}

/** Sum lengths of all sequences in MS.
    @param SEXP pointer to MS object
    @return Integer specifying total # of bases in MS object
*/
SEXP rph_ms_totalSeqLengths(SEXP msP)
{
  SEXP result;
  MS *inputMS = SEXP_to_group(msP);
  
  int i, *resultP;
  long totalSeqLengths = 0;
  for(i=0; i < inputMS->nseqs; i++)
    {
      totalSeqLengths += strlen(inputMS->seqs[i]);
    }
  
  PROTECT(result = NEW_INTEGER(1));
  resultP = INTEGER_POINTER(result);
  resultP[0] = totalSeqLengths;
  UNPROTECT(1);
  return result;
}

/** Build a Markov Model based on sequences for a given MS
    @param ms SEXP pointer to MS object containing at least one sequence
    @param norderP SEXP pointer to an integer containing the order of the markov model to create
    @param pseudoCountP SEXP pointer to an integer containing the pseudoCount to add 
    @param considerReverseP SEXP pointer to logical whether to take into account reverse complement freqs
    @result R list of Matrices that make up the Markov Model
*/
SEXP rph_mm_build(SEXP ms, SEXP norderP, SEXP pseudoCountP, SEXP considerReverseP)
{
  int nOrder, j, numProtect = 0, pseudoCount, considerReverse;
  nOrder = INTEGER_VALUE(norderP);
  pseudoCount = INTEGER_VALUE(pseudoCountP);
  considerReverse = INTEGER_VALUE(considerReverseP);
  if ((considerReverse != 0) && (considerReverse != 1)) 
    die("considerReverse must be a logical value either 0 or 1");
  List *markovModel;  //Matrices from order 0 to norder that make up the Markov Model
  SEXP SMarkovModel;
  MS *group = SEXP_to_group(ms);

  //Build the Markov Model
  markovModel = mm_build(group, nOrder, pseudoCount, considerReverse);

  //Return a list of Matrices to R
  PROTECT(SMarkovModel = allocVector(VECSXP, lst_size(markovModel)));
  numProtect++;
  for (j = 0; j < lst_size(markovModel); j++)
    SET_VECTOR_ELT(SMarkovModel, j, Matrix_to_SEXP(lst_get_ptr(markovModel, j)));

  UNPROTECT(numProtect);
  return SMarkovModel;
}


/** Return a subset of the MS sequences.
    @param msP SEXP pointer to MS object
    @param rowsP List of integers indicating the sequences to return
    @result subset of input MS sequences
*/
SEXP rph_ms_square_brackets(SEXP msP, SEXP rowsP) {
  MS *ms, *newMs;
  char **names, **seqs;
  int *rows=NULL, i, spec, nrow, *offsets, numprotect=0;

  ms = (MS*)EXTPTR_PTR(msP);
  ms_register_protect(ms);
  
  if (rowsP != R_NilValue) {
    nrow = LENGTH(rowsP);
    PROTECT(rowsP = AS_INTEGER(rowsP));
    rows = INTEGER_POINTER(rowsP);
    numprotect++;
  } else nrow = ms->nseqs;


  names = smalloc(nrow*sizeof(char*));
  seqs = smalloc(nrow*sizeof(char*));
  offsets = smalloc(nrow*sizeof(int));
  for (i=0; i < nrow; i++) {
    checkInterrupt();
    if (rows == NULL) spec = i;
    else spec = rows[i]-1; //convert to 0-based numbers from R indices
    if (spec < 0 || spec >= (ms->nseqs)) { 
      names[i] = copy_charstr("");
      seqs[i] = copy_charstr("");
      offsets[i] = 0;
    } else {
      names[i] = copy_charstr(ms->names[spec]);
      seqs[i] = copy_charstr(ms->seqs[spec]);
      offsets[i] = ms->idx_offsets[spec];
    }
}

  newMs = ms_new(seqs, names, nrow, ms->alphabet, 0,1);
  newMs->idx_offsets = offsets;
  if (numprotect > 0) UNPROTECT(numprotect);
  return rph_ms_new_extptr(newMs);
}


SEXP rph_ms_gc_content(SEXP sequencesP) {
  Vector *gc_content;
  ListOfLists *result;
  gc_content = ms_gc_content(SEXP_to_group(sequencesP));
  result = lol_new(1);
  lol_push_dbl(result, gc_content->data, gc_content->size, "gc.content");
  return rph_listOfLists_to_SEXP(result);
}


/** RPHAST split sequences given a fixed window length
    @param sequencesP SEXP pointer to MS object containing sequences
    @param windowSize Size of the window we will use to chunk each sequence (i.e. split every 700 bases)
    @result MS object containing sequences split up.
*/
SEXP rph_ms_split_size(SEXP sequencesP, SEXP windowSizeP)
{
  int  inputSeqNum, outputSeqNum, nextCutAt, windowSize, outputNseqs, inputSeqLen;
  MS *outputMS;
  MS *inputMS;
  
  windowSize = INTEGER_VALUE(windowSizeP);
  
  inputMS = SEXP_to_group(sequencesP);
  
  //Calculate number of sequences we will have in the split up MS
  for(inputSeqNum=0, outputNseqs=0; inputSeqNum < inputMS->nseqs; inputSeqNum++)
    {
      outputNseqs += ceil(strlen(inputMS->seqs[inputSeqNum]) / (double)windowSize);
    }
  
  //Setup output MS specifying number of sequences and copying alphabet from inputMS
  outputMS = ms_new(NULL, NULL, outputNseqs, inputMS->alphabet, 0, 1);
  //Allocate memory for pointers to each sequence
  outputMS->seqs = (char**)smalloc(outputNseqs * sizeof(char*));
  outputMS->names = (char**)smalloc(outputNseqs * sizeof(char*));
  outputMS->idx_offsets = (int*)smalloc(outputNseqs * sizeof(int));

  //For each input sequence (unsplit)
  for(inputSeqNum=0, outputSeqNum=0; inputSeqNum < inputMS->nseqs; inputSeqNum++)
    {
      nextCutAt = 0;
      //Break the sequence into sub sequences until it is smaller than the window size and add to default group
      while ((strlen(inputMS->seqs[inputSeqNum])-nextCutAt) > windowSize)  
        {
          //Copy name of sequence
          outputMS->names[outputSeqNum] = (char*)smalloc((strlen(inputMS->names[inputSeqNum])+1) * sizeof(char));
          strncpy(outputMS->names[outputSeqNum], inputMS->names[inputSeqNum], strlen(inputMS->names[inputSeqNum]));
          outputMS->names[outputSeqNum][strlen(inputMS->names[inputSeqNum])] = '\0';
	   
          //Copy subset of sequence bases
          outputMS->seqs[outputSeqNum] = (char*)smalloc((windowSize + 1) * sizeof(char));
          strncpy(outputMS->seqs[outputSeqNum], inputMS->seqs[inputSeqNum] + nextCutAt, windowSize);
          outputMS->seqs[outputSeqNum][windowSize] = '\0';
	   
          //Set index offset of new sequence
          outputMS->idx_offsets[outputSeqNum] = inputMS->idx_offsets[inputSeqNum] + nextCutAt;
	   
          nextCutAt += windowSize;
          outputSeqNum++;
        } 
    
      //Make a copy of the sequence that is smaller than the window size and add to outputMS
      if ((strlen(inputMS->seqs[inputSeqNum])-nextCutAt) > 0)
        {
          //Copy name of sequence
          outputMS->names[outputSeqNum] = (char*)smalloc((strlen(inputMS->names[inputSeqNum])+1) * sizeof(char*));
          strncpy(outputMS->names[outputSeqNum], inputMS->names[inputSeqNum], strlen(inputMS->names[inputSeqNum]));
          outputMS->names[outputSeqNum][strlen(inputMS->names[inputSeqNum])] = '\0';
	   
          //Copy subset of sequence bases 
          inputSeqLen = strlen(inputMS->seqs[inputSeqNum]);
          outputMS->seqs[outputSeqNum] = (char*)smalloc((inputSeqLen-nextCutAt + 1) * sizeof(char*));
          strncpy(outputMS->seqs[outputSeqNum], inputMS->seqs[inputSeqNum] + nextCutAt, inputSeqLen-nextCutAt);
          outputMS->seqs[outputSeqNum][inputSeqLen-nextCutAt] = '\0';
	   
          //Set index offset of the new sequence
          outputMS->idx_offsets[outputSeqNum] = inputMS->idx_offsets[inputSeqNum] + nextCutAt;
	   
          nextCutAt += windowSize;
          outputSeqNum++;
        }
    }

  return group_to_SEXP(outputMS, FALSE);
}


/** Print sequences in MS object to FASTA file.
    @param msP SEXP pointer to MS object
    @param fileP SEXP pointer to string specifying location of where to save object (stdout if NULL)
*/
SEXP rph_ms_printSeq_fasta(SEXP msP, SEXP fileP) {
  MS *ms;
  ms = (MS*)EXTPTR_PTR(msP);
  ms_register_protect(ms);

  if (fileP != R_NilValue)
    ms_print_to_file(CHARACTER_VALUE(fileP), 
                     ms);
  else
    ms_print_fasta(stdout, ms);
  return R_NilValue;
}


/** Print sequences in MS object to screen
    @param msP SEXP pointer to MS object
*/
SEXP rph_ms_printSeq(SEXP msP) {
  MS *ms;
  ms = (MS*)EXTPTR_PTR(msP);
  ms_register_protect(ms);

  ms_print(stdout, ms);

  return R_NilValue;
}


/** RPHAST split sequences given a list of locations
    @param sequencesP SEXP pointer to MS object
    @param featuresP SEXP Pointer to GFF_Set specifying window ranges
    @result List of pointers to MSAs
*/
SEXP rph_ms_split_gff(SEXP sequencesP, SEXP featuresP)
{
  int outputSeqNum, currentFeature, currentSequence, lengthOfSubSequence, outputNseqs;
  GFF_Feature *feature;
  GFF_Set *gff;
  MS *outputMS;
  MS *inputMS;

  gff = (GFF_Set*)EXTPTR_PTR(featuresP);
  inputMS = SEXP_to_group(sequencesP);
  outputNseqs = lst_size(gff->features);
	
  //Setup output MS specifying number of output sequences and the alphabet to use
  outputMS = ms_new(NULL, NULL, outputNseqs,  inputMS->alphabet, 0, 1);
  //Allocate memory for pointers to each sequence
  outputMS->seqs = (char**)smalloc(outputNseqs * sizeof(char*));
  outputMS->names = (char**)smalloc(outputNseqs * sizeof(char*));
  outputMS->idx_offsets = (int*)smalloc(outputNseqs * sizeof(int));
	
  //For each window
  for (currentFeature = 0, outputSeqNum=0; currentFeature < lst_size(gff->features); currentFeature++) {
    feature = (GFF_Feature*)lst_get_ptr(gff->features, currentFeature);		
    
    //For each sequence
    for (currentSequence = 0; currentSequence < inputMS->nseqs; currentSequence++) {
			
      //fprintf(stderr, "seqname %s, feature %s\n", ms->names[0], feature->seqname->chars); //DEBUG
      //Find the sequence with the same name as the feature
      if (strcmp(inputMS->names[currentSequence], feature->seqname->chars) == 0) {
				
        //Copy the name of the sequence for inclusion in outputMS
        //        if(feature->start <= inputMS->idx_offsets[currentSequence])
        feature->start = feature->start - inputMS->idx_offsets[currentSequence];
        feature->end = feature->end - inputMS->idx_offsets[currentSequence];
  
        printf("Start %d, end %d.\n", feature->start, feature->end);
        outputMS->names[outputSeqNum] = (char*)smalloc(feature->seqname->length + 1 * sizeof(char*));
        outputMS->names[outputSeqNum] = copy_charstr(feature->seqname->chars);

        if (feature->start < 1)
          feature->start = 1;

        lengthOfSubSequence = (feature->end - feature->start)+1;

        //Check if window specified is within bounds of sequence
       
        if ((feature->start > 0) && ((feature->start-1+ lengthOfSubSequence) <= strlen(inputMS->seqs[currentSequence]))) {
					
          //Copy substring from sequence and create new MSA to hold subsequence specified by window
          outputMS->seqs[outputSeqNum] = (char*)smalloc((lengthOfSubSequence + 1) * sizeof(char));
          strncpy(outputMS->seqs[outputSeqNum], inputMS->seqs[currentSequence] + (feature->start-1), lengthOfSubSequence);
          outputMS->seqs[outputSeqNum][lengthOfSubSequence] = '\0';
					
          //Set index offset
          outputMS->idx_offsets[outputSeqNum] = feature->start-1 + inputMS->idx_offsets[currentSequence];
          outputSeqNum++;
        }else
          Rf_warning("Feature %s from %d to %d is outside sequence boundaries\n", feature->seqname->chars, feature->start, feature->end);
        break;
      }
    }
  }
  return group_to_SEXP(outputMS, FALSE);
}


/** RPHAST Read Position Weight Matrix (PWM) from file
    @param filenameP Full path to MEME file containing at least one PWM
    @result List of Matrices (PWMs)
*/
SEXP rph_pwm_read(SEXP filenameP)
{
  int i;
  SEXP result;
  List *pwms;

  if (filenameP == R_NilValue)
    die("ERROR: No positition weight matrix filename was provided");

  pwms = pwm_read(CHARACTER_VALUE(filenameP));

  PROTECT(result = allocVector(VECSXP, lst_size(pwms)));
  for (i = 0; i < lst_size(pwms); i++)
    SET_VECTOR_ELT(result, i, Matrix_to_SEXP(lst_get_ptr(pwms, i)));
  UNPROTECT(1);
  return result;
}

/** RPHAST return the number of sequences in an MS object
    @param msP SEXP pointer to MS object
    @result Integer containing number of sequences in MS object
*/
SEXP rph_ms_nseq(SEXP msP) {
  MS *ms = (MS*)EXTPTR_PTR(msP);
  SEXP result;
  int *resultP;
  PROTECT(result = NEW_INTEGER(1));
  resultP = INTEGER_POINTER(result);
  resultP[0] = ms->nseqs;
  UNPROTECT(1);
  return result;
}

/** Translate a R SEXP object to a C Matrix object
    @param matrixP SEXP object representing a matrix
    @result C matrix containing data represented in SEXP object
*/
Matrix *SEXP_to_Matrix(SEXP matrixP)
{
  int row, col, nrows, ncols;
  Matrix *result;
  double MatrixElement;

  nrows = INTEGER(GET_DIM(matrixP))[0];
  ncols = INTEGER(GET_DIM(matrixP))[1];
  //printf("SEXP_to_matrix dimensions nrows=%d, ncols=%d\n", nrows, ncols);
  result = mat_new(nrows, ncols);
  for (row = 0; row < nrows; row++) {
    for (col = 0; col < ncols; col++) {
      MatrixElement = REAL(matrixP)[col * nrows + row];
      // printf("Matrix element [row=%d, col=%d ] = %f\n", row, col, MatrixElement);
      mat_set(result, row, col, MatrixElement);
    }
  }
  return result;
}


// CGD: ADDED 17-October-2013
// Function modified from ms_score to return a SEXP which can be passed into R.
// List structure represesnts: scores on the plus strand and scores on the minus strand, for both the motif and background models. 
/** RPHAST Compute scores for a given sequence.
    @param seqName Name of the sequence to scan
    @param seqData Char representation of the DNA sequence
    @param seqLen Length of the DNA sequence
    @param seqIdxOff ?? Offsets for the sequence 
    @param seqAlphLen ?? Length of the seqeunce alphabet
    @param MarkovMatrices Pointer to Markov Model used for the background
    @param pwm Pointer to matrix to be used as Position Weight Matrix (PWM) to score for matches against
    @param reverseCmpPWM Pointer to the reverse complement matrix
    @param conservative Pointer to logical whether to treat regions containing 'N' as possible binding sites
    @result An R list structure representing the posterior probability of every base in every strand in both backgorund and motif states.
*/
SEXP ms_posterior_list(char *seqName, char *seqData, int seqLen, int seqIdxOff, int seqAlphLen, List *MarkovMatrices, Matrix *pwm, Matrix *reverseCmpPWM, int conservative) {
  int i, k,j,l,col;
  double MMprob, PWMprob=0, ReversePWMprob=0;
//  GFF_Set *scores = gff_new_set(); // CGD: Changed return type.
  double *MMprobs = (double*)smalloc((pwm->nrows+1) * sizeof(double));    //Sliding window of mmOrder previous MM probabilities

  // CGD: Create R list type.
  SEXP scores, motif, motif_plus, motif_minus, background;
  PROTECT(scores = allocVector(VECSXP, 2));
  PROTECT(motif = allocVector(VECSXP, 2));
  PROTECT(background = allocVector(REALSXP, seqLen-(pwm->nrows)+1));
  PROTECT(motif_plus = allocVector(REALSXP, seqLen-(pwm->nrows)+1));
  PROTECT(motif_minus = allocVector(REALSXP, seqLen-(pwm->nrows)+1));
  SET_VECTOR_ELT(motif, 0, motif_plus);
  SET_VECTOR_ELT(motif, 1, motif_minus);
  SET_VECTOR_ELT(scores, 0, motif);
  SET_VECTOR_ELT(scores, 1, background);

  // Set list names.
  SEXP model_COL_Names, score_COL_Names;
  PROTECT(model_COL_Names = NEW_CHARACTER(2));
  SET_STRING_ELT(model_COL_Names, 0, mkChar("Forward"));
  SET_STRING_ELT(model_COL_Names, 1, mkChar("Reverse"));
  setAttrib(motif, R_NamesSymbol, model_COL_Names);

  PROTECT(score_COL_Names = NEW_CHARACTER(2));
  SET_STRING_ELT(score_COL_Names, 0, mkChar("MotifModel"));
  SET_STRING_ELT(score_COL_Names, 1, mkChar("Background"));
  setAttrib(scores, R_NamesSymbol, score_COL_Names);


  double *mmPlusScore= REAL(motif_plus);
  double *mmMinusScore= REAL(motif_minus);
  double *backgroundScore= REAL(background);

  if ((conservative != 0) && (conservative != 1))
    die("ERROR: Conserverative (boolean) value must be 0 or 1");

  if (seqLen < pwm->nrows)  //Check to see if the sequence is shorter than the pwm
    return(scores);

  for (i = 0; i <= pwm->nrows; i++)                                                     //Calculate MM scores from sites 0 to pwm->nrows
    if (i < seqLen)
      MMprobs[i] = calcMMscore(seqData, i, MarkovMatrices, conservative);

  for (i = 0; i <= seqLen-(pwm->nrows); i++) {                          //For each base in the sequence
    PWMprob = 0; MMprob = 0; ReversePWMprob = 0;

    for (k = 0, j = i; k < pwm->nrows; k++, j++) {              //Sum PWM, ReversePWM, MM probabilities for score calculation
      col = basetocol(seqData[j]);
      if (col >= 0)
        {
          PWMprob += mat_get(pwm, k, col);
          ReversePWMprob += mat_get(reverseCmpPWM, k, col);
          MMprob += MMprobs[k];
        }
      else {
        if (conservative)
          {
            PWMprob = log(0);                   //If we get something other than the expected language (i.e. A,C,T,G) i.e. N, then our probability is -Inf
            ReversePWMprob = log(0);
            break;
          }
        else
          {
            PWMprob = 0;
            ReversePWMprob = 0;
          }
      }
    }

    if (i < (seqLen - pwm->nrows)) { //Only if there are more bases in this sequence to test
      for (l = 0; l < pwm->nrows; l++)          //Shift probs left to make room for next
        MMprobs[l] = MMprobs[l + 1];

      MMprobs[pwm->nrows-1] = calcMMscore(seqData, i+pwm->nrows,  //Calculate MM probability for site at (i+pwm->nrows)
                                          MarkovMatrices, conservative);
    }

    // CGD: Record values here.
    backgroundScore[i] = MMprob;
    mmPlusScore[i] = PWMprob;
    mmMinusScore[i] = ReversePWMprob;

  }
  sfree(MMprobs);
  unprotect(7); // CGD: Must unprotect.
  return(scores);
}



// CGD (dankoc@gmail.com): ADDED 17-October-2013
// Intended to return vectors of the posterior probabilies under the motif and background models into R
// Formats as a list() structure.
/** RPHAST Compute scores for each base in every sequence in every group
    @param inputMSP SEXP pointer to MS object containing sequences
    @param pwmP SEXP pointer to R matrix to be used as Position Weight Matrix (PWM) to score for matches against
    @param markovModelP SEXP pointer to Markov Model (List of Matrices in R) generated by build.mm for provided MS object
    @param nOrderP  SEXP pointer to Integer specifying order of Markov Model used
    @param conservativeP SEXP pointer to logical whether to treat regions containing 'N' as possible binding sites
    @result An R list structure representing the posterior probability of every base in every sequence in every group
*/
SEXP rph_ms_posterior(SEXP inputMSP, SEXP pwmP, SEXP markovModelP, SEXP nOrderP, SEXP conservativeP)
{
  int i, currentSequence, conservative;
  
  Matrix *mm, *pwm, *reverseCompPWM;
  List *MarkovMatrices;

  MS *inputMS;
  inputMS = SEXP_to_group(inputMSP);

  // CGD: Return types.
  SEXP scores, locus_score;
  PROTECT(scores = allocVector(VECSXP, inputMS->nseqs));
  conservative = asLogical(conservativeP);
  pwm = SEXP_to_Matrix(pwmP);
  reverseCompPWM = mat_reverse_complement(pwm);

  MarkovMatrices = lst_new_ptr(length(markovModelP));
  for (i = 0; i < length(markovModelP); i++) {
    mm = SEXP_to_Matrix(VECTOR_ELT(markovModelP, i));
    lst_push_ptr(MarkovMatrices, mm);
  }

  for (currentSequence = 0; currentSequence < inputMS->nseqs; currentSequence++) { //For each sequence in the inputMS
    locus_score = ms_posterior_list(inputMS->names[currentSequence], inputMS->seqs[currentSequence],
                      strlen(inputMS->seqs[currentSequence]), inputMS->idx_offsets[currentSequence],
                      strlen(inputMS->alphabet), MarkovMatrices, pwm, reverseCompPWM,
                      conservative);
    SET_VECTOR_ELT(scores, currentSequence, locus_score);
  }

  unprotect(1);
  return(scores);
}

/** RPHAST Compute scores for each base in every sequence in every group
    @param inputMSP SEXP pointer to MS object containing sequences
    @param pwmP SEXP pointer to R matrix to be used as Position Weight Matrix (PWM) to score for matches against
    @param markovModelP SEXP pointer to Markov Model (List of Matrices in R) generated by build.mm for provided MS object
    @param nOrderP  SEXP pointer to Integer specifying order of Markov Model used
    @param conservativeP SEXP pointer to logical whether to treat regions containing 'N' as possible binding sites
    @param thresholdP Minimum score a binding site must receive to be returned
    @param strandP Which strands to look through and which results to return
    @result Pointer to GFF_Set object in C, a Features object in R containing the locations & scores of binding sites
*/
SEXP rph_ms_score(SEXP inputMSP, SEXP pwmP, SEXP markovModelP, SEXP nOrderP, SEXP conservativeP, SEXP thresholdP, SEXP strandP)
{
  int site, i, currentSequence, conservative;
  double threshold;
  char *strand;
  GFF_Feature *score;
  Matrix *mm, *pwm, *reverseCompPWM;
  List *MarkovMatrices;
  GFF_Set *groupScores, *scores;
  MS *inputMS;
  ListOfLists *result;

  threshold = NUMERIC_VALUE(thresholdP);
  conservative = asLogical(conservativeP);
  strand = (char*)translateChar(STRING_ELT(strandP, 0));

  pwm = SEXP_to_Matrix(pwmP);
  reverseCompPWM = mat_reverse_complement(pwm);

  inputMS = SEXP_to_group(inputMSP);
	
  result = lol_new(1);

  MarkovMatrices = lst_new_ptr(length(markovModelP));
  for (i = 0; i < length(markovModelP); i++) {
    mm = SEXP_to_Matrix(VECTOR_ELT(markovModelP, i));
    lst_push_ptr(MarkovMatrices, mm);
  }

  groupScores = gff_new_set();
	
  //For each sequence calculate score each site
  for (currentSequence = 0; currentSequence < inputMS->nseqs; currentSequence++) { //For each sequence in the inputMS
    scores = ms_score(inputMS->names[currentSequence], inputMS->seqs[currentSequence], 
                      strlen(inputMS->seqs[currentSequence]), inputMS->idx_offsets[currentSequence],
                      strlen(inputMS->alphabet), MarkovMatrices, pwm, reverseCompPWM, 
                      conservative, threshold, strand);

    //Add scores to list of scores for current group
    for (site = 0; site < lst_size(scores->features); site++) {
      score = (GFF_Feature*)lst_get_ptr(scores->features, site);
      lst_push_ptr(groupScores->features, score);
    }
    scores->features = NULL;
    gff_free_set(scores);
  }
  lol_push_gff(result, groupScores, "scores");

  //printf("Finished with compute Scores\n");
  return rph_listOfLists_to_SEXP(result);
}


/** RPHAST Simulate a sequence given a Markov Model
    @param mmP Markov Model created with build.mm()
    @param norderP Order of Markov Matrices specified in mmP
    @param alph_sizeP Alphabet size of data used to generate mmP
    @param lengthP Length of sequence to generate
    @result SEXP pointer to MS object containing a single simulated sequence
*/
SEXP rph_ms_simulate(SEXP mmP, SEXP norderP, SEXP alph_sizeP, SEXP lengthP) //don't need to pass the order (size of mmP-1) or the alphabet size (num cols in mmP of order 0)
{
  MS *outputMS;
  char *seq;
  char *name = (char*)smalloc(100*sizeof(char));
  List *MarkovMatrices;
  int norder, alph_size, mmNum;
  unsigned int *length;
  GetRNGstate();	//Get RNG state from R
    
  norder = INTEGER_VALUE(norderP);
  alph_size = INTEGER_VALUE(alph_sizeP);
  length = (unsigned int*)INTEGER(lengthP);//(unsigned int)INTEGER_VALUE(lengthP);
  int nSeqs = Rf_nrows(lengthP);

  MarkovMatrices = lst_new_ptr(length(mmP));
  for (mmNum = 0; mmNum < length(mmP); mmNum++)
    lst_push_ptr(MarkovMatrices, SEXP_to_Matrix(VECTOR_ELT(mmP, mmNum)));

  outputMS = ms_new(NULL, NULL, nSeqs,  NULL, 0, 1);

  //Setup output MS specifying number of sequences and copying alphabet from inputMS
  //Allocate memory for pointers to each sequence
  outputMS->seqs = (char**)smalloc(sizeof(char*)*nSeqs);
  outputMS->names = (char**)smalloc(sizeof(char*)*nSeqs);
  outputMS->idx_offsets = (int*)smalloc(sizeof(int)*nSeqs);

  for(int i = 0; i < nSeqs; i++) {
    seq = ms_simulate(MarkovMatrices, norder, alph_size, length[i]);
    //printf("simulated seq %s\n", seq);

    snprintf(name, 10, "S%d", i);
    outputMS->names[i] = (char*)smalloc((strlen(name)+1) * sizeof(char));
    strncpy(outputMS->names[i], name, strlen(name));
    outputMS->names[i][strlen(name)] = '\0';

    //Copy subset of sequence bases
    outputMS->seqs[i] = (char*)smalloc((length[i] + 1) * sizeof(char));
    strncpy(outputMS->seqs[i], seq, length[i]);
    outputMS->seqs[i][length[i]] = '\0';
    sfree(seq);
	   
    //Set length & index offset of new sequence
    outputMS->idx_offsets[i] = 0;
  }
  PutRNGstate();
	   
  return group_to_SEXP(outputMS, FALSE);
}

/****** Accessors functions**********/

/** Return list of sequences contained in the MS object
    @param SEXP pointer to MS object
    @return List of strings, each string containing the bases for a sequence
*/
SEXP rph_ms_seqs(SEXP msP) {
  SEXP result;
  MS *ms = (MS*)EXTPTR_PTR(msP);
  int seq;
  
  PROTECT(result = NEW_CHARACTER(ms->nseqs));
  ms_register_protect(ms);

  for (seq = 0; seq < ms->nseqs; seq++) { 
    SET_STRING_ELT(result, seq, mkChar(ms->seqs[seq]));
  }
  
  UNPROTECT(1);
  return result;
}

/** Return list of names for each sequence contained in the MS object.
    @param SEXP pointer to MS object
    @return LIst of string, each string containing the name of a sequence
*/
SEXP rph_ms_seqNames(SEXP msP) {
  MS *ms = (MS*)EXTPTR_PTR(msP);
  SEXP result;
  int i;
  
  if (ms->names==NULL) return R_NilValue;
  PROTECT(result = NEW_CHARACTER(ms->nseqs));
  for (i=0; i<ms->nseqs; i++) 
    SET_STRING_ELT(result, i, mkChar(ms->names[i]));
  UNPROTECT(1);
  return result;
}

/** Return the alphabet (valid non-missing characters) used to make the sequences in an MS object.
    @param SEXP pointer to MS object
    @return String containing characters that make up the alphabet (valid non-missing characters)
*/
SEXP rph_ms_alphabet(SEXP msP) {
  MS *ms = (MS*)EXTPTR_PTR(msP);
  SEXP result;
  
  if (ms->alphabet==NULL) return R_NilValue;
  PROTECT(result = NEW_CHARACTER(1));
  SET_STRING_ELT(result, 0, mkChar(ms->alphabet));
  UNPROTECT(1);
  return result;
}

/** Return the Index offsets for each sequence in an MS object.
    @param SEXP pointer to MS object
    @return List of integers, one for each sequence specifying its index offset.
*/
SEXP rph_ms_idxOffsets(SEXP msP) {
  MS *ms = (MS*)EXTPTR_PTR(msP);
  SEXP result;
  int *resultP, i;
  
  PROTECT(result = NEW_INTEGER(ms->nseqs));
  resultP = INTEGER_POINTER(result);
  for (i=0; i<ms->nseqs; i++)  //TODO: can probably skip this loop and return ms->idx_offsets
    resultP[i] = ms->idx_offsets[i];
  UNPROTECT(1);
  return result;
}

/** Return the lengths for each sequence in an MS object.
    @param SEXP pointer to MS object.
    @return List of integers, one for each sequence specifying its length
*/
SEXP rph_ms_lengths(SEXP msP) {
  MS *ms = (MS*)EXTPTR_PTR(msP);
  SEXP result;
  int *resultP, i;
  
  PROTECT(result = NEW_INTEGER(ms->nseqs));
  resultP = INTEGER_POINTER(result);
  for (i=0; i<ms->nseqs; i++)
    resultP[i] = strlen(ms->seqs[i]);
  UNPROTECT(1);
  return result;
}
