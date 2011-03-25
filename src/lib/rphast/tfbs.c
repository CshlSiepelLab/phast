/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/*****************************************************
tfbs.c
Transcription Factor binding sites functions

Nick Peterson
Last updated: 3/18/2011
*****************************************************/
#include <msa.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <stdio.h>

#include <list_of_lists.h>
#include <lists.h>
#include <stacks.h>
#include <misc.h>
#include <gff.h>
#include <category_map.h>
#include <hashtable.h>
#include <sufficient_stats.h>
#include <local_alignment.h>
#include <indel_history.h>
#include <tfbs.h>


//////////////////////////////////////////
List *read_pwm(const char *filename)
{ 
  List *result;
  Matrix *pwm = NULL;
  int i, currBase, nBases = 0;
  FILE * F;
  char *motifName;
  String *line = str_new(STR_MED_LEN);
  List *l = lst_new_ptr(3);
  List *probabilitiesStr = lst_new_ptr(4);
  List *probabilitiesDbl;
  Regex *pssm_re = NULL;
  Regex *motif_name_re = NULL;
  int alphabetLength;

  result = lst_new_ptr(1);
  //letter-probability matrix: alength= 4 w= 8 nsites= 2 E= 1.5e+004 

  pssm_re = str_re_new("^letter-probability matrix: alength= ([0-9]+) w= ([0-9]+)");
  motif_name_re = str_re_new("^MOTIF[[:space:]]+(.+?)[[:space:]].*");
  //open file
  if( (F = fopen(filename, "r")) == NULL)
   die("ERROR: Unable to open MEME file containing PWM \n");
  currBase = 0;
  nBases = -1;
  //For each line in the MEME file
  while ((str_readline(line, F)) != EOF) {

    //If line matches Motif name
    if(str_re_match(line, motif_name_re, l, 1) > 0) {
      motifName = copy_charstr(((String*)lst_get_ptr(l, 1))->chars);
      printf("motifName=%s\n", motifName);
    }
    //If line matches beginning of a probability matrix
    else if (str_re_match(line, pssm_re, l, 2) > 0) {
      //Extract the alphabet size & number of bases in matrix
      
      if(str_as_int((String*)lst_get_ptr(l, 1), &alphabetLength) != 0)
      {
        die("ERROR: Unable to parse 'alength=' from MEME file, expected integer, read %s", ((String*)lst_get_ptr(l, 1))->chars);
      }
      if (str_as_int((String*)lst_get_ptr(l, 2), &nBases) != 0)
      {
        die("ERROR: Unable to parse 'w=' from MEME file, expected integer, read %s ", ((String*)lst_get_ptr(l, 2))->chars);
      }
      currBase = 0;
      pwm = mat_new(nBases, alphabetLength);
      continue;
    //If this row contains matrix data
    } else if(currBase < nBases)
    {
      //Parse row of probabilities
      str_double_trim(line);
      str_split(line, NULL, probabilitiesStr);
      probabilitiesDbl = str_list_as_dbl(probabilitiesStr);
      for(i=0; i< lst_size(probabilitiesDbl); i++)
      {
        mat_set(pwm, currBase, i, lst_get_dbl(probabilitiesDbl, i));
      }
      currBase++;
    } else if((currBase == nBases) && (pwm != NULL))
    {
      //Push full matrix
      lst_push_ptr(result, pwm);
      pwm = NULL;
    }
  }
  return result;
}


//////////////////////////////////////////////////
List *ms_read(const char *filename, const char *alphabet)
{
  
  FILE * F;
  if( (F = fopen(filename, "r")) == NULL)
   die("ERROR: Unable to open FASTA file \n");

  List *mss = lst_new_ptr(10);
  static Regex *descrip_re = NULL;
  int maxlen, i, nseqs, j, do_toupper, line_no;
  String *line = str_new(STR_MED_LEN);
  List *l = lst_new_ptr(2);
  String *new_str = str_new(STR_MED_LEN);
  char **sequence, **name;
  MSA *msa;

  if (descrip_re == NULL) 
    descrip_re = str_re_new("[[:space:]]*>[[:space:]]*([^[:space:]]+)");

  line_no=1;
  while ((str_readline(line, F)) != EOF) {
    if (str_re_match(line, descrip_re, l, 1) > 0) {
      if(new_str->length != 0)
      {

        sequence = smalloc((new_str->length+1)*sizeof(char*));
	sequence[0] = copy_charstr(new_str->chars);
        msa = msa_new(sequence, name, 1, new_str->length, copy_charstr(alphabet));
        lst_push_ptr(mss, msa);
        
      }
      name = smalloc(((String*)lst_get_ptr(l, 1))->length*sizeof(char*));
      name[0] = copy_charstr(((String*)lst_get_ptr(l, 1))->chars);
      str_free((String*)lst_get_ptr(l, 0));

      new_str = str_new(STR_MED_LEN);
      continue;
    }

    str_double_trim(line);
    if (line->length == 0) continue;

    if (new_str == NULL) 
      die("ERROR in FASTA file: non-blank line preceding first description ('>') line.\n");

    str_append(new_str, line);
    checkInterruptN(line_no++, 1000);
  }
 
  if (new_str->length > 0)
  {
    sequence = smalloc((new_str->length+1)*sizeof(char*));
    sequence[0] = copy_charstr(new_str->chars);
    msa = msa_new(sequence, name, 1, new_str->length, copy_charstr(alphabet));
    lst_push_ptr(mss, msa);
  }
  

  if (lst_size(mss) == 0)
    die("ERROR: empty FASTA file.\n");

  /* upcase chars unless there are lowercase characters in the alphabet */
  do_toupper = 1;    
 
  //lst_free(l);
  str_free(line);

  for (i = 0; i < lst_size(mss); i++) {
    String *n, *s;

    msa = (MSA*)lst_get_ptr(mss, i);
   // scan chars and adjust if necessary 
    for (j = 0; j < strlen(msa->seqs[0]); j++) {
      msa->seqs[0][j] = do_toupper ? toupper(msa->seqs[0][j]) : msa->seqs[0][j];
      if (isalpha(msa->seqs[0][j]) && msa->inv_alphabet[(int)msa->seqs[0][j]] == -1 && get_iupac_map()[(int)msa->seqs[0][j]] == NULL) 
        msa->seqs[0][j] = 'N';   // assume 'N' if unrecognized letter 
    }

  }  

  return mss;
}

///////////////////////////////////////////////////////////////////////////////
Matrix *mm_build_helper(ListOfLists *group, int norder, Matrix *lowerOrderMM)
{
	int numMS, alph_size, i, endOfSeqs, j, ignore, tup_idx, l, alph_idx, s;
	MSA *ms;
	double sum = 0, val;
	double lowerOrderMMval = 1;   
	Vector *freqs = vec_new(int_pow(4, norder)); 
          
    vec_zero(freqs);  
    numMS = lst_size(group->lst);
    
       alph_size = strlen(((MSA*)lst_get_ptr(((ListOfLists*)lst_get_ptr(group->lst, 0))->lst, 0))->alphabet);

      for (i = 0, endOfSeqs = FALSE; endOfSeqs == FALSE; i++) {
        checkInterruptN(i, 10000);
        endOfSeqs = TRUE;
      
        for (j = 0; j < numMS; j++) {
          ignore = 0;
          tup_idx = 0;
        
          for (l = 0; !ignore && l < norder; l++) {
            ms = (MSA*)lst_get_ptr(((ListOfLists*)lst_get_ptr(group->lst, j))->lst, 0);
            if (ms->length <= (i+l))
               continue;
            endOfSeqs = FALSE;
            char c = ms->seqs[0][i+l];
            if ((alph_idx = ms->inv_alphabet[(int)c]) == -1)
              ignore = 1;
            else 
              tup_idx += alph_idx * int_pow(alph_size, (norder-l)-1); 
          }
          
          if (!ignore) {
            vec_set(freqs, tup_idx, 
                           vec_get(freqs, tup_idx) + 1); 
          }
        }
    }  
       
       //Normalize
       sum = 0;
       for (i = 0; i < freqs->size; i++) sum += vec_get(freqs, i);
       vec_scale(freqs, 1.0/sum);

       for(s=0; s<freqs->size; s++)
       {
	     printf("vec %f ", vec_get(freqs, s));
       }
       printf("\n");

       Matrix *mm = mat_new(int_pow(alph_size, (norder-1)), alph_size);

       for (i=0; i< freqs->size; i=i+alph_size)
	   { 
	     sum = 0;
 	     for(j=0; j<alph_size; j++)
	     {
	       sum += vec_get(freqs, i+j);
	     }
	     if (sum == 0) //To avoid a divide by 0
	       sum = 1; //If the sum is 0, then in the next loop, the sum can be anything but zero.
	     for(j=0; j<alph_size; j++)
	     {
	      // if (lowerOrderMM != NULL) //If lowerOrderMM is NULL it defaults to 1, not modifying freqs just counted
	      //   lowerOrderMMval = mat_get(lowerOrderMM, floor(i / (alph_size*alph_size)), (i+j % alph_size));
	       val = vec_get(freqs, i+j)/sum;// * lowerOrderMMval;
	       
	       mat_set(mm, i/alph_size, j, val);
	     }
	   }
	   
	   return mm;
}

////////////////////////////////////////////////////////////////////
List *mm_build(ListOfLists *msGroupListP, int norderP)
{
   int i, j, k, s, currentGroup, numMS, endOfSeqs;
   int numGroups = lst_size(msGroupListP->lst);
   MSA *ms;
   ListOfLists *group;
   List *MarkovModelList;
   Matrix *mm = NULL;
   
   k=norderP; //tuple size

   printf("called mm_build\n");
   MarkovModelList = lst_new_ptr(numGroups);
   for(currentGroup=0; currentGroup< numGroups; currentGroup++)
   {
     List *MatrixList = lst_new_ptr(norderP);
       /*double sum = 0;    
       int  ignore, tup_idx, l, alph_idx;
       vec_zero(freqs);  
       int alph_size;*/

     group = (ListOfLists*)lst_get_ptr(msGroupListP->lst, currentGroup);
     
	for(i=1; i <= norderP; i++)
	{
     mm = mm_build_helper(group, i, mm);
	 lst_push_ptr(MatrixList, mm);
	}
	lst_push_ptr(MarkovModelList, MatrixList);
/*
     if (numMS > 0)
     {
       alph_size = strlen(((MSA*)lst_get_ptr(((ListOfLists*)lst_get_ptr(group->lst, 0))->lst, 0))->alphabet);

      for (i = 0, endOfSeqs = FALSE; endOfSeqs == FALSE; i++) {
        checkInterruptN(i, 10000);
        endOfSeqs = TRUE;
      
        for (j = 0; j < numMS; j++) {
          ignore = 0;
          tup_idx = 0;
        
          for (l = 0; !ignore && l < k; l++) {
            ms = (MSA*)lst_get_ptr(((ListOfLists*)lst_get_ptr(group->lst, j))->lst, 0);
            if (ms->length <= (i+l))
               continue;
            endOfSeqs = FALSE;
            char c = ms->seqs[0][i+l];
            if ((alph_idx = ms->inv_alphabet[(int)c]) == -1)
              ignore = 1;
            else 
              tup_idx += alph_idx * int_pow(alph_size, (k-l)-1); 
          }
          
          if (!ignore) {
            vec_set(freqs, tup_idx, 
                           vec_get(freqs, tup_idx) + 1); 
          }
        }
    }  
       
       //Normalize
       sum = 0;
       for (i = 0; i < freqs->size; i++) sum += vec_get(freqs, i);
       vec_scale(freqs, 1.0/sum);

       for(s=0; s<freqs->size; s++)
       {
	     printf("vec %f ", vec_get(freqs, s));
       }
       printf("\n");

       Matrix *mm = mat_new(int_pow(alph_size, (k-1)), alph_size);

       for (i=0; i< freqs->size; i=i+alph_size)
	   { 
	     sum = 0;
 	     for(j=0; j<alph_size; j++)
	     {
	       sum += vec_get(freqs, i+j);
	     }
	     if (sum == 0) //To avoid a divide by 0
	       sum = 1; //If the sum is 0, then in the next loop, the sum can be anything but zero.
	     for(j=0; j<alph_size; j++)
	     {
	       mat_set(mm, i/alph_size, j, vec_get(freqs, i+j)/sum);
	     }
	   }
	   lst_push_ptr(mmList, mm);
     }
     */

   /** @todo calculate probabilities for sites 0 to norderP -1 and save as list.. or something */  
   }
   printf("done build mm\n");
   return MarkovModelList;
}

//////////////////////////////////////////////////////
int basetocol(char base)
{
 //Put next base from sequence in opening at end of array (opened by shift above)
    switch(base)
    {
      case 'A': return 0;
      case 'C': return 1;
      case 'G': return 2;
      case 'T': return 3;
      default: fprintf(stderr, "Encountered unknown base %c\n", base); return -1;
    }
}

////////////////////////////////////////////////////////////////////
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

//////////////////////////////////////////////////////////////////////////////////
List *mm_compute_scores(MSA *ms, List *MarkovMatrices, Matrix *pwm, int mmOrder)
{
  //ms = Current sequence 
  //mm = List of Markov Matrices (order 1..n) built for current sequence's group
  //pwm = Position Weight matrix from MEME file
  int i, k, j, alph_size, col, currentMMorder;
  double val;
  Matrix *mm;
  double MMprob, PWMprob;
  List *scores = lst_new_dbl(ms->length);
  int *previousMMbases = (int*)smalloc(mmOrder * sizeof(int));  //Sliding window of mmOrder previous bases
  double *previousMMprob = (double*)smalloc(mmOrder * sizeof(double)); //Sliding window of mmOrder previous MM probabilities
  alph_size = strlen(ms->alphabet);
  MMprob =0;
  PWMprob =0;

  //Zero out previous bases
  for(i=0; i< mmOrder; i++)
  {
    previousMMbases[i] = 0;
  }
  
  for(i=0; i<mmOrder; i++)
  {
    previousMMprob[i] = 0;
  }
    
  //For each base in the sequence
  for(i=0; i< ms->length; i++)
  {
     PWMprob =0;
     MMprob = 0;
     
     //Choose appropreate Markov Matrix for this site
     if(i+1 >= mmOrder)
       currentMMorder = mmOrder;
     else
       currentMMorder = i;
       
     //Shift characters in string left 1
     for(k=0; k<mmOrder-1; k++)
     {
        previousMMbases[k] = previousMMbases[k+1];
     }
     
     //Shift probs left to make room for next
     for(k=0; k<mmOrder-1; k++)
     {
        previousMMprob[k] = previousMMprob[k+1];
     }
     
     //Put next base from sequence in opening at end of array (opened by shift above)
    previousMMbases[mmOrder-1] = basetocol(ms->seqs[0][i]);
    
    //Calculate PWM match probability
    if ((i+pwm->nrows) > ms->length) //If PWM is larger than window, mark as not possible match and move on.
      PWMprob = log(0); // -inf probability
    else
    {
      for(k=0, j=i; k < pwm->nrows; k++, j++)
      {
           col = basetocol(ms->seqs[0][j]);
          if (col >= 0)
            val = mat_get(pwm, k, col);
          else
            val = 0; //If we get something other than the expected language (i.e. A,C,T,G) 
            			//such as N, then the probability of TFBS @ this base is -inf, log(0) = -inf
            			//TODO: Optimize this so we don't keep checking other bases
          PWMprob += log(val); //log(p1*p2*...pn) 
      }
    
      //Determine row in MM that corrisponds to previous mmOrder windows
	  mm =  lst_get_ptr(MarkovMatrices, currentMMorder-1);
      j = basesToRow(previousMMbases, currentMMorder-1, alph_size);
      if (j >= 0) 
        val = mat_get(mm, j, previousMMbases[mmOrder-1]);
      else
        val = log(0);
      
      //Append new MM probability to previousMMprobs
      previousMMprob[mmOrder-1] = log(val);
    
      //Add site MM log probabilities 
      for(k=0; k < mmOrder; k++)
      {
         MMprob += previousMMprob[k];
      }
    }
    
    if (PWMprob == log(0))
      lst_push_dbl(scores, log(0));
     else
     {
      lst_push_dbl(scores, (PWMprob - MMprob));
      //printf("PWMprob=%f, MMprob=%f\n", PWMprob, MMprob);
      }
    

  }
  
  return scores; // log(pwm/mm)
}

