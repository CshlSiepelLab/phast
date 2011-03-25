
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

Melissa Hubisz
Nick Peterson
Last updated: 2/15/2011
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

//ncols and nrows have a name collision with Matrix->ncols, Matrix->nrows
#undef ncols
#undef nrows

SEXP *ms_group_by_gc(List *MSs, int ngroups);
char *mm_simulate_seq(List *mmP, int norderP, int alph_size, int lengthP);


SEXP Matrix_to_SEXP(Matrix *mat)
{
  SEXP result;
  int col, row;
  PROTECT(result = allocMatrix(REALSXP, mat->nrows, mat->ncols));

  for(row=0; row< mat->nrows; row++)
  {
    for(col=0; col< mat->ncols; col++)
    {
      REAL(result)[col*mat->nrows+row] = mat_get(mat, row, col);
     // printf("Setting [row=%d, col=%d] = %f\n", row, col, mat_get(mat, row, col));
      
    }
  }
  UNPROTECT(1);

  return result;
}

/** Read multiple sequences from file.
    @param[in] filenameP File path to Multiple Sequence file
    @param[in] alphabetP Alphabet of characters used in sequences
    @result Multiple sequence object 
    @note Each sequence gets its own MS
*/
SEXP rph_ms_read(SEXP filenameP, SEXP alphabetP) {
  int i, numProtect=0;
  MSA *msa;
  List *mss;
  SEXP result_list, Rmsa;

  if (filenameP == R_NilValue)
    die("ERROR: No sequences filename was provided");

  if (alphabetP == R_NilValue)
    mss = ms_read(CHARACTER_VALUE(filenameP), NULL);
  else
    mss = ms_read(CHARACTER_VALUE(filenameP), CHARACTER_VALUE(alphabetP));

  PROTECT(result_list = allocVector(VECSXP, lst_size(mss)));
  numProtect++;
  for(i=0; i< lst_size(mss); i++)
  {
    msa = (MSA*)lst_get_ptr(mss, i);
    PROTECT(Rmsa = rph_msa_new_extptr(msa));
    numProtect++;
    SET_VECTOR_ELT(result_list, i, Rmsa);
  }
  UNPROTECT(numProtect);
  return result_list; 
}



SEXP rph_ms_seqNames(SEXP msListP) {
  int numMS, i; 
  List *msList;
  MSA *ms;
  SEXP result;

  numMS = length(msListP);
  //msList = lst_new_ptr(numMS);
  PROTECT(result = NEW_CHARACTER(numMS));
  for (i = 0; i < numMS; i++) {
    ms = (MSA*)EXTPTR_PTR(VECTOR_ELT(msListP, i));
    //lst_push_ptr(msList, ms);
    //msa_registenew_characprotect(ms);
    //printf("length of ms %d\n", ms->length);
    SET_STRING_ELT(result, i, mkChar(ms->names[0]));
  }
  UNPROTECT(1);
  
  return result;
}


ListOfLists *SEXP_to_groups(SEXP Slol)
{
  int numGroups = length(Slol);
  int numMS;
  int i,j;
  ListOfLists *result = lol_new(numGroups);
  ListOfLists *group;
  SEXP groupNames = getAttrib(Slol, R_NamesSymbol);
  SEXP Sgroup;
  SEXP sequenceNames;
  MSA *ms;
  char *sequenceName, *quantileName;
  for(i=0; i< numGroups; i++)
  {
    //quantileName = CHAR(STRING_ELT(groupNames, i));
    quantileName = "";
    Sgroup = VECTOR_ELT(Slol, i);
    numMS = length(Sgroup);
    group = lol_new(numMS);
    lol_push_lol(result, group, quantileName);
    sequenceNames = getAttrib(Sgroup, R_NamesSymbol);
    for(j=0; j<numMS; j++)
    {
      //sequenceName = CHAR(STRING_ELT(sequenceNames, j));
      sequenceName = "";
      ms = ((MSA*)EXTPTR_PTR(VECTOR_ELT(Sgroup, j)));
      lol_push_msa_ptr(group, ms, sequenceName);
    }
  }
  return result;
}

SEXP groups_to_SEXP(ListOfLists *groups)
{
  int i, j, k, numGroups, numMS, numProtected, temp;
  SEXP result;
  List *group;
  MSA *ms;
  
  numProtected =0;
  numGroups = lst_size(groups->lst);
  
  for(i=0, temp=numGroups; i< numGroups; i++)
  {
    if(lst_size((List*)(((ListOfLists*)lst_get_ptr(groups->lst, i))->lst)) == 0)
      temp--;
  }
  PROTECT(result = allocVector(VECSXP, temp));
  numProtected++;
  for(i=0, k=0; i< numGroups; i++)
  {
    SEXP Sgroup;
    group = (List*)(((ListOfLists*)lst_get_ptr(groups->lst, i))->lst);
    numMS = lst_size(group);
    if (numMS > 0)
    {
      PROTECT(Sgroup = allocVector(VECSXP, numMS));
	  numProtected++;
	  for(j=0; j<numMS; j++)
	  {
	    ms = (MSA*)lst_get_ptr(((ListOfLists*)lst_get_ptr(group, j))->lst, 0);
	    SEXP Rmsa;
	    PROTECT(Rmsa = rph_msa_new_extptr(ms));
	    numProtected++;
	    SET_VECTOR_ELT(Sgroup, j, Rmsa);
	  }
      SET_VECTOR_ELT(result, k, Sgroup);
      k++;
    }
  }
  UNPROTECT(numProtected);
  return result;
}

SEXP rph_build_mm(SEXP msListP, SEXP norderP)
{
  int mmOrder, i, j, numProtect =0;
  List *mms, *MarkovModel;
  ListOfLists *msList = SEXP_to_groups(msListP);
  mmOrder = INTEGER_VALUE(norderP);
  mms = mm_build(msList, mmOrder);
  SEXP result;
  
  PROTECT(result = allocVector(VECSXP, lst_size(mms)));
  numProtect++;
  for(i=0; i< lst_size(mms); i++)
  {
    MarkovModel = (List*)lst_get_ptr(mms, i);
    SEXP SMarkovModel;
    PROTECT(SMarkovModel = allocVector(VECSXP, lst_size(MarkovModel)));
    numProtect++;
    for(j=0; j< lst_size(MarkovModel); j++)
    {
      SET_VECTOR_ELT(SMarkovModel, j, Matrix_to_SEXP(lst_get_ptr(MarkovModel, j)));
    }
    SET_VECTOR_ELT(result, i, SMarkovModel); 
  }
  UNPROTECT(numProtect);
  return result;
}

/** RPHAST group sequences by their GC content
    @param msListP List of pointers to MSAs
    @param ngroupsP Pointer to number of groups to group into
    @result List of groups which each contain a list of pointers to MSAs
*/
SEXP rph_ms_group_by_gc_content(SEXP msListP, SEXP ngroupsP)
{
  int i;
  int ngroups = INTEGER_VALUE(ngroupsP);
  int numMS = length(msListP);
  List *mss = lst_new_ptr(numMS);
  MSA *ms;
  
  if (ngroups < 1) 
   die("Number of gc content groups must be at least 1");  
  
  
  for(i=0; i < numMS; i++)
  {
    ms = (MSA*)EXTPTR_PTR(VECTOR_ELT(msListP, i));
    lst_push_ptr(mss, ms);
  }
  return ms_group_by_gc(mss, ngroups);

}

/** RPHAST split sequences given a list of windows
    @param msListP List of pointers to MSAs
    @param windowsP Pointer to GFF_Set specifying window ranges
    @result List of pointers to MSAs 
*/
SEXP rph_ms_clip_to_supplied_windows(SEXP msListP, SEXP windowsP)
{
  int i, j, currentFeature, currentSequence, numMS, lengthOfSubSequence, numProtect=0;
  GFF_Feature *feature;
  int gcCount;
  GFF_Set *gff;
  MSA *ms, *subSequence;
  char **name, **sequence;
  SEXP result_list, Rmsa;
  List *clipedSubSequences, *mss;


  gff = (GFF_Set*)EXTPTR_PTR(windowsP);  
  
  numMS = length(msListP);
  
  clipedSubSequences = lst_new_ptr(lst_size(gff->features));

  //For each window
  for(currentFeature=0; currentFeature < lst_size(gff->features); currentFeature++)
  {
    feature = (GFF_Feature*)lst_get_ptr(gff->features, currentFeature);
    fprintf(stderr, "feature name %s, start %d, end %d\n", feature->seqname->chars, feature->start, feature->end); //DEBUG
    //For each sequence
    for(currentSequence=0; currentSequence < numMS; currentSequence++)
    {   
      ms = (MSA*)EXTPTR_PTR(VECTOR_ELT(msListP, currentSequence));
      fprintf(stderr, "seqname %s, feature %s\n", ms->names[0], feature->seqname->chars); //DEBUG
      //Find the sequence with the same name as the feature
      if(strcmp(ms->names[0], feature->seqname->chars) == 0)
      {
	//Copy the name of the sequence for making a new MSA
        name = smalloc(feature->seqname->length+1*sizeof(char*));
        name[0] = copy_charstr(feature->seqname->chars);
	
	lengthOfSubSequence = (feature->end - feature->start);
	sequence = smalloc((2) * sizeof(char*));
	//Check if window specified is within bounds of sequence
	if ((feature->start > 0) && (feature->start+lengthOfSubSequence <= feature->end)) 
	{  
	  //Copy substring from sequence and create new MSA to hold subsequence specified by window
	  sequence[0] = (char*)smalloc((lengthOfSubSequence +1) * sizeof(char));
	  strncpy(sequence[0], ms->seqs[0]+feature->start, lengthOfSubSequence);
	  sequence[0][lengthOfSubSequence] = '\0';
	  fprintf(stderr, "Clipped sequence %s from %d to %d, seq %s\n", feature->seqname->chars, feature->start, feature->end, sequence[0]);
	  
	  subSequence = msa_new(sequence, name, 1, lengthOfSubSequence, ms->alphabet);
	  lst_push_ptr(clipedSubSequences, subSequence); 

	}else
          Rf_warning("Feature %s from %d to %d is outside sequence boundaries\n", feature->seqname->chars, feature->start, feature->end);
	break;
      }
    }
  }

  PROTECT(result_list = allocVector(VECSXP, lst_size(clipedSubSequences)));
  numProtect++;
  for(i=0; i< lst_size(clipedSubSequences); i++)
  {
    ms = (MSA*)lst_get_ptr(clipedSubSequences, i);
    PROTECT(Rmsa = rph_msa_new_extptr(ms));
    numProtect++;
    SET_VECTOR_ELT(result_list, i, Rmsa);
  }
  UNPROTECT(numProtect);
  return result_list; 

}


/** RPHAST Lookup row of Markov Matrix corrisponding to bases i.e 'AT' == 3
    @param norderP Order of Markov Matrix represented
    @param alph_size Alphabet size of data represented in Markov Matrix
    @param Integer specifying row of Markov Matrix mapped to given bases
*/
int basesToRow(int *previousBases, int norderP, int alph_size)
{
  int col = 0;
  int i;
  for (i=norderP-1; i > 0; i--)
  {
    col += int_pow(alph_size, i) * previousBases[i]; 
  }

  return col;
}

/** RPHAST Read Position Weight Matrix (PWM) from file
    @param filenameP Full path to MEME file containing at least one PWM
    @result List of Matrices (PWMs) 
*/
SEXP rph_pwm_read(SEXP filenameP)
{
  int i, numProtect=0;
  SEXP result;
  Matrix *pwm;
  ListOfLists *pwms;

  if (filenameP == R_NilValue)
    die("ERROR: No positition weight matrix filename was provided");

  pwms = pwm_read(CHARACTER_VALUE(filenameP));
  
  PROTECT(result = allocVector(VECSXP, lst_size(pwms)));
  for(i=0; i< lst_size(pwms); i++)
  {
    SET_VECTOR_ELT(result, i, Matrix_to_SEXP(lst_get_ptr(pwms, i)));
  }
  UNPROTECT(1);
  return result;
}

/** Translate a R SEXP object to a C Matrix object
    @param Smatrix SEXP object representing a matrix
    @result C matrix containing data represented in SEXP object
*/
Matrix *SEXP_to_Matrix(SEXP Smatrix) 
{
  int row, col, nrows, ncols;	
  Matrix *result;
  double MatrixElement;
  nrows = INTEGER(GET_DIM(Smatrix))[0];
  ncols = INTEGER(GET_DIM(Smatrix))[1];
  printf("SEXP_to_matrix dimensions nrows=%d, ncols=%d\n", nrows, ncols);
  result = mat_new(nrows, ncols);
  for(row=0; row < nrows; row++)
  {
    for(col=0; col < ncols; col++)
    {
      MatrixElement = REAL(Smatrix)[col*nrows + row];
     // printf("Matrix element [row=%d, col=%d ] = %f\n", row, col, MatrixElement);
      mat_set(result, row, col, MatrixElement);
    }
  }
  return result;
}

/** RPHAST Compute scores for each base in every sequence in every group
    @param msGroupListP List of groups containing pointers to MSAs
    @param pwmP Position Weight Matrix (PWM) to score for matches against
    @param mmListP List of Markov Matrices (one per group)
    @param mmOrder Order of Markov Matrix specified in mmListP
*/
SEXP rph_compute_scores(SEXP msGroupListP, SEXP pwmP, SEXP mmListP, SEXP mmOrderP)
{
  int numGroups, numPositiveScores=0,  site, currentGroupNum, i, numMS, currentSequence, numProtected =0;
  char currentStrand = "-";
  MSA *ms;
  double score;
  Matrix *mm;
  Matrix *pwm;
  List *scores, *MarkovMatrices;
  ListOfLists *groups;
  ListOfLists *currentGroup;
  SEXP result;
  SEXP groupScores;
  SEXP scoreEntry;
  int mmOrder = INTEGER_VALUE(mmOrderP);
  pwm = SEXP_to_Matrix(pwmP); 
  
  //Groups of sequences
  groups = SEXP_to_groups(msGroupListP);
  
  numGroups = lst_size(groups->lst);
  if (numGroups != length(mmListP))
    die("ERROR: mmListP has a different number of groups than mmListP, you must have one markov model per group");

  //PROTECT(result = allocVector(VECSXP, numGroups));
  PROTECT(result = allocVector(VECSXP, numGroups));
  numProtected++;
  for(currentGroupNum=0; currentGroupNum < numGroups; currentGroupNum++)
  {
    SEXP MarkovModel = VECTOR_ELT(mmListP, currentGroupNum);
    MarkovMatrices = lst_new_ptr(length(MarkovModel));
    for(i=0; i<length(MarkovModel); i++)
    {
       mm = SEXP_to_Matrix(VECTOR_ELT(MarkovModel, i));
       lst_push_ptr(MarkovMatrices, mm); 
    }
    numPositiveScores =0;
    currentGroup = lst_get_ptr(groups->lst, currentGroupNum);
    numMS = lst_size(currentGroup->lst);
    PROTECT(groupScores = allocVector(VECSXP, numMS));
    numProtected++;
    for(currentSequence=0; currentSequence < numMS; currentSequence++) //For each sequence in the group
    {
      ms = (MSA*)lst_get_ptr(((ListOfLists*)(lst_get_ptr(currentGroup->lst, currentSequence)))->lst, 0); //Sequence #currentSequence in group #currentGroup
      scores = mm_compute_scores(ms, MarkovMatrices, pwm, mmOrder);
      for(site=0; site<lst_size(scores); site++)
      {
        score = lst_get_dbl(scores, site);
        if (score > 0) //Note, we can throw away negative scores???
        {
          PROTECT(scoreEntry = allocVector(VECSXP, 5));
      	  numProtected++;
          
          //Object structure
          //      Score | Strand | Base | Seq Name | Seq Length 
          // i.e.  70.4     +       554    Chimp.5     48
          SEXP Sscore;
          PROTECT(Sscore = allocVector(REALSXP, 1));
          numProtected++;
          REAL(Sscore)[0] = score;
          SET_VECTOR_ELT(scoreEntry, 0, Sscore);
          
          SEXP Sstrand;
          PROTECT(Sstrand = allocVector(STRSXP, 1));
          numProtected++;
          if (currentStrand == "-")
            SET_STRING_ELT(Sstrand, 0, mkChar("-"));
          else
            SET_STRING_ELT(Sstrand, 0, mkChar("+"));
          SET_VECTOR_ELT(scoreEntry, 1, Sstrand);
          
          SEXP Sbase;
          PROTECT(Sbase = allocVector(INTSXP, 1));
          numProtected++;
          INTEGER(Sbase)[0] = site;
          SET_VECTOR_ELT(scoreEntry, 2, Sbase);
          
          SEXP SseqName;
          PROTECT(SseqName = allocVector(STRSXP, 1));
          numProtected++;
          SET_STRING_ELT(SseqName, 0, mkChar(ms->names[0]));
          SET_VECTOR_ELT(scoreEntry, 3, SseqName);
          
          SEXP SseqLen;
          PROTECT(SseqLen = allocVector(INTSXP, 1));
          numProtected++;
          INTEGER(SseqLen)[0] = ms->length;
          SET_VECTOR_ELT(scoreEntry, 4, SseqLen);
          
          SET_VECTOR_ELT(groupScores, numPositiveScores, scoreEntry);
          numPositiveScores++;
        }
      }
     
    }
    SET_VECTOR_ELT(result, currentGroupNum, groupScores);
    
  }
  UNPROTECT(numProtected);
  printf("Finished with compute Scores\n");
  return result;
}

SEXP rph_alph_size(SEXP msListP)
{
  SEXP ans;
  PROTECT(ans = allocVector(REALSXP, 1));
  REAL(ans)[0] = strlen(((MSA*)EXTPTR_PTR(VECTOR_ELT(VECTOR_ELT(msListP, 0), 0)))->alphabet);
  UNPROTECT(1);
  return ans;
}

/** RPHAST Generate a sequence given a Markov Model
    @param mmP List groups with each element being a list of Markov Models from (0...nth order)
    @param norderP Order of Markov Matrices specified in mmP
    @param alph_sizeP Alphabet size of data used to generate mmP
    @param lengthP Length of sequence to generate
    @result List of pointers to MSAs (one per model) containing generated sequences    
*/
SEXP rph_mm_simulate_seq(SEXP mmP, SEXP norderP, SEXP alph_sizeP, SEXP lengthP)
{
  MSA *ms;
  Matrix *mm;
  char *seq;
  char **sequence;
  char **name;
  List *MarkovMatrices;
  int norder, alph_size, length, numMarkovMatrices, groupNum, numGroups, mmNum, numProtected;
  SEXP result, SMarkovMatrices, Sms, Sseqs;
  numProtected=0;
  
  norder = INTEGER_VALUE(norderP);
  alph_size = INTEGER_VALUE(alph_sizeP);
  length = INTEGER_VALUE(lengthP);
  numGroups = length(mmP);
  PROTECT(result = allocVector(VECSXP, numGroups));
  numProtected++;
  
  for(groupNum=0; groupNum < numGroups; groupNum++)
  {
    SMarkovMatrices = VECTOR_ELT(mmP, groupNum);
    MarkovMatrices = lst_new_ptr(length(SMarkovMatrices));
    for(mmNum=0; mmNum < length(SMarkovMatrices); mmNum++)
    {
      lst_push_ptr(MarkovMatrices, SEXP_to_Matrix(VECTOR_ELT(SMarkovMatrices, mmNum)));
    }
    
    seq = mm_simulate_seq(MarkovMatrices, norder, alph_size, length);
    printf("simulated seq %s\n", seq);
    sequence = smalloc((length+1)*sizeof(char*));
    sequence[0] = copy_charstr(seq);
    name = smalloc(100*sizeof(char*));
    name[0] = copy_charstr("GC group");
  
  
    ms = msa_new(sequence, name, 1, length, NULL);
    PROTECT(Sseqs = allocVector(VECSXP, 1));
    numProtected++;
    PROTECT(Sms = rph_msa_new_extptr(ms));
    numProtected++;
    SET_VECTOR_ELT(Sseqs, 0, Sms);
    SET_VECTOR_ELT(result, groupNum, Sseqs);
  }
  UNPROTECT(numProtected);
  return result;
}

/** Generate a sequence given a Markov Model 
    @param mm Markov Model containing probabilities used to generate sequence
    @param norderP Order of Markov Model mm
    @param alph_size Alphabet size of data used to generate mm
    @param lenthP Length of sequence to generate
    @result Char string containing generated sequence i.e. "ATAGCGC" for length 7
*/
char *mm_simulate_seq(List *mmP, int norderP, int alph_size, int lengthP)
{
  printf("Called mm_simulate_seq\n");
   int currentMMorder, i, j, l, k, base, a=0, t=0, g=0, c=0;
   char *result = (char*)smalloc((lengthP+1) * sizeof(char));
   int *previousBases = smalloc(norderP-1* sizeof(int));
   Matrix *mm;
   srand ( time(NULL) );

   for(base=0; base<lengthP; base++)
   {
     //j = column corrisponding to last two bases
     
     if(base+1 >= norderP)
       currentMMorder = norderP;
     else
       currentMMorder = base+1;
     
     mm = (Matrix*)lst_get_ptr(mmP, currentMMorder-1);
     
     j = basesToRow(previousBases, currentMMorder-1, alph_size);
	
     double r = (random() % 100)/(double)100;
   
     for (i =0; i< mm->ncols; i++)
     {
       r = r - mat_get(mm, j, i);
       if (r <= 0)
         break; 	
     }
     if (i >= mm->ncols)
       i = mm->ncols-1;

     //Shift characters in string left 1
     for(k=0; k<norderP-1; k++)
     {
        previousBases[k] = previousBases[k+1];
     }
     previousBases[norderP-2] = i;
     switch(i)
     {
	  case 0: result[base] = 'A'; break;
	  case 1: result[base] = 'C'; break;
	  case 2: result[base] = 'G'; break;
	  case 3: result[base] = 'T'; break;	
	  //case 0: a++; break;
	  //case 1: c++; break;
	  //case 2: g++; break;
	  //case 3: t++; break;
     }
   }
   //printf("a=%d, c=%d, g=%d, t=%d\n", a, c, g, t);
   //printf("\n");
   result[base]= '\0';
   return result;
}

/** Convert floating point number to Char string
    @param num Floating point number
    @result String containing floating point number
*/
char *dtoa(double num)
{
  char *str = smalloc(STR_MED_LEN);
  sprintf(str,"%f", num);
  return str;
}

/** RPHAST group sequences by their GC content level
    @param MSs List of pointers to MSAs
    @param ngroups Number of GC content level groups
    @result List of groups each containing list of pointers to MSAs each of which has one sequence
*/
SEXP *ms_group_by_gc(List *MSs, int ngroups)
{
  int numMS = lst_size(MSs);
  int i, j, gcCount;
  List *gcContentScores = lst_new_dbl(numMS);
  List *gcContentScoresUnsorted = lst_new_dbl(numMS);
  ListOfLists *groups = lol_new(ngroups);
  MSA *ms;
  SEXP result;
  
  
  for(i=0; i<numMS; i++)
  {
    ms = (MSA*)lst_get_ptr(MSs, i);
    //GC calc
    gcCount = 0;
    //Caclulate GC content
    for(j=0; j< ms->length; j++)
    {
      if ((ms->seqs[0][j] == 'G') || (ms->seqs[0][j] == 'C'))
        gcCount = gcCount +1;
    }
    fprintf(stderr, " gc found %d, gc score %f\n", gcCount, ((double)gcCount/(double)ms->length));

    lst_push_dbl(gcContentScores, ((double)gcCount/(double)ms->length));
    lst_push_dbl(gcContentScoresUnsorted, ((double)gcCount/(double)ms->length));
  }

 //Sort list of scores via quick sort
  lst_qsort_dbl(gcContentScores, ASCENDING);

  //Calculate parameter for quantile calculation
  double q[ngroups];
  for(i=1; i<= ngroups; i++)
  {
    q[i-1] = i*((double)1/(double)ngroups);
    printf("q[%d] = %f\n", i, q[i]);
  }
  
  //Create array for quantile boundaries
  double quantile_vals[ngroups+1]; 

  //Calculate quantiles
  lst_dbl_quantiles(gcContentScores, q, ngroups+1, quantile_vals);

  //Print quantiles  SET_VECTOR_ELT(result, i group);
  for(i=0; i < ngroups; i++)
  {
    printf("%f ", quantile_vals[i]);
  }

  //Create groups per quantile
  for(i=0; i< ngroups; i++)
  {
    lol_push_lol(groups, lol_new(1), dtoa(quantile_vals[i])); 
  }

  //Group sequences by quantiles
  for(i=0; i < numMS; i++)
  {
    ms = (MSA*)lst_get_ptr(MSs, i);
    for(j=ngroups-1; j>=0; j--)
    {
	  if (lst_get_dbl(gcContentScoresUnsorted, i) >= quantile_vals[j])
	  {
	    ListOfLists *lol = lol_find_lol(groups, dtoa(quantile_vals[j]));	
	    lol_push_msa_ptr(lol, ms, ms->names[0]);
        
	    printf("quantile %f, val %f , %d, %d\n", quantile_vals[j],lst_get_dbl(gcContentScoresUnsorted, i), i, j);
	    break;
	  }
    }
     if (j == -1) //If it wasn't larger than any of the quantile values it belongs in the lowest one
	  {
	    ListOfLists *lol = lol_find_lol(groups, dtoa(quantile_vals[0]));
	    lol_push_msa_ptr(lol, ms, ms->names[0]);
	     printf("quantile %f, val %f , %d, %d\n", quantile_vals[j],lst_get_dbl(gcContentScoresUnsorted, i), i, 0);
	  }
  }
  result = groups_to_SEXP(groups);
  return result;
}



