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
MS *ms_new(char **seqs, char **names, int nseqs, const char *alphabet, double rangeLow, double rangeHigh) {
  int i;
  MS *ms;

  //Ranges are from 0 to 1 representing GC content for this group
  if ((rangeLow < 0) || (rangeLow > 1) || (rangeHigh < 0) || (rangeHigh > 1))
    die("ERROR: Creating a new GC group for sequences: range values must be between 0 and 1");
  
  ms = (MS*)smalloc(sizeof(MS));
  ms->rangeLow = rangeLow;
  ms->rangeHigh = rangeHigh;
  ms->seqs = seqs;
  ms->names = names;
  ms->nseqs = nseqs;
  if (alphabet != NULL) {
    ms->alphabet = (char*)smalloc((strlen(alphabet) + 1) * sizeof(char));
    strcpy(ms->alphabet, alphabet);
  }
  else {
    ms->alphabet = (char*)smalloc((strlen(DEFAULT_ALPHABET)+1) * sizeof(char));
    strcpy(ms->alphabet, DEFAULT_ALPHABET);
  }
  ms->missing = DEFAULT_MDATA_CHARS;
  
  for (i = 0; i < NCHARS; i++) { 
    ms->inv_alphabet[i] = -1;
    ms->is_missing[i] = 0;
  }
  for (i = 0; ms->alphabet[i] != '\0'; i++)
    ms->inv_alphabet[(int)ms->alphabet[i]] = i;
  for (i = 0; ms->missing[i] != '\0'; i++)
    ms->is_missing[(int)ms->missing[i]] = 1;

  return ms;
}

////////////////////////////////////////
void ms_free(MS *ms) {
  int i;
  
  for (i = 0; i < ms->nseqs; i++) {
    if (ms->names != NULL && ms->names[i] != NULL) 
      sfree(ms->names[i]);
    if (ms->seqs != NULL && ms->seqs[i] != NULL) 
      sfree(ms->seqs[i]);
  }
  if (ms->names != NULL) sfree(ms->names);
  if (ms->seqs != NULL) sfree(ms->seqs);
  if (ms->alphabet != NULL) sfree(ms->alphabet);
  sfree(ms);

}

///////////////////////////////////////////////////////////////
void ms_print(FILE *F, MS *ms) {
  int i;
  
  for (i = 0; i < ms->nseqs; i++) {
    checkInterrupt();
    fprintf(F, "  Name    %s\n", ms->names[i]);
    fprintf(F, "  Offset  %d\n", ms->idx_offsets[i]);
    fprintf(F, "  Seq     %s\n", ms->seqs[i]);
    if(i != (ms->nseqs-1))
      fprintf(F, "\n");
  }
}

//////////////////////////////////////////////////////////
void ms_print_fasta(FILE *F, MS *ms) {
  int i, j, k, seqLen;
  
  for (i = 0; i < ms->nseqs; i++) {
    checkInterrupt();
    fprintf(F, ">%s\n", ms->names[i]);
    seqLen = strlen(ms->seqs[i]);
    for (j = 0; j < seqLen; j += OUTPUT_LINE_LEN) {
      checkInterruptN(j, 100);
      for (k = 0; k < OUTPUT_LINE_LEN && j + k < seqLen; k++) 
        fprintf(F, "%c", ms->seqs[i][j+k]);
      fprintf(F, "\n");
    }
  }
}

///////////////////////////////////////
void ms_print_to_file(const char *filename, MS *ms) {
  FILE *outfile = phast_fopen(filename, "w");
  ms_print_fasta(outfile, ms);
  phast_fclose(outfile);
}



//////////////////////////////////////////
List *pwm_read(const char *filename) {
  List *result;
  Matrix *pwm = NULL;
  int i, currBase, nBases = 0;
  FILE * F;
  //  char *motifName;
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
  //open PWM file
  F = phast_fopen(filename, "r");
  currBase = 0;
  nBases = -1;
  //For each line in the MEME file
  while ((str_readline(line, F)) != EOF) {
    //If line matches Motif name
    if (str_re_match(line, motif_name_re, l, 1) > 0) {
      //      motifName = copy_charstr(((String*)lst_get_ptr(l, 1))->chars);
      //printf("motifName=%s\n", motifName);
    }
    //If line matches beginning of a probability matrix
    else if (str_re_match(line, pssm_re, l, 2) > 0) {
      //Extract the alphabet size & number of bases in matrix

      if (str_as_int((String*)lst_get_ptr(l, 1), &alphabetLength) != 0)
        die("ERROR: Unable to parse 'alength=' from MEME file, expected integer, read %s", ((String*)lst_get_ptr(l, 1))->chars);
      if (str_as_int((String*)lst_get_ptr(l, 2), &nBases) != 0)
        die("ERROR: Unable to parse 'w=' from MEME file, expected integer, read %s ", ((String*)lst_get_ptr(l, 2))->chars);
      currBase = 0;
      if (nBases <= 0) //We must have at least one base in the PWM
        die("ERROR: No Position Weight Matrices were detected in the provided PWM file");
      if (alphabetLength <= 0) //We must have a positive alphabet length
        die("ERROR: Alphabet lengh specified in PWM file must be greater than zero");
      pwm = mat_new(nBases, alphabetLength);
      mat_set_all(pwm, -1);
      continue;
      //If this row contains matrix data
    } else if (currBase < nBases) {
      //Parse row of probabilities
      str_double_trim(line);
      str_split(line, NULL, probabilitiesStr);
      probabilitiesDbl = str_list_as_dbl(probabilitiesStr);
      for (i = 0; i < lst_size(probabilitiesDbl); i++)
        mat_set(pwm, currBase, i, log(lst_get_dbl(probabilitiesDbl, i)));
      currBase++;
    } else if ((currBase == nBases) && (pwm != NULL)) {
      //Push full matrix
      lst_push_ptr(result, pwm);
      pwm = NULL;
    }
  }
  if (currBase == nBases && pwm != NULL) 
    lst_push_ptr(result, pwm);
  else if (pwm != NULL) 
    die("Premature end of PWM file\n");
  str_re_free(motif_name_re);
  str_re_free(pssm_re);
  phast_fclose(F);
  return result;
}


//////////////////////////////////////////////////


/* return TRUE if alphabet has lowercase letters, FALSE otherwise */
int ms_alph_has_lowercase(MS *ms) {
  int i;
  
  for (i = 0; ms->alphabet[i] != '\0'; i++)
    if (ms->alphabet[i] >= 'a' && ms->alphabet[i] <= 'z')
      return TRUE;
  return FALSE;
}

MS *ms_read(const char *filename, const char *alphabet) {
  List *names = lst_new_ptr(10);
  List *seqs = lst_new_ptr(10);
  static Regex *descrip_re = NULL;
  int i, nseqs, j, do_toupper, line_no;
  String *line = str_new(STR_MED_LEN);
  List *l = lst_new_ptr(2);
  String *n, *s, *new_str = NULL;
  MS *ms;
  FILE * F;

  F = phast_fopen(filename, "r");

  if (descrip_re == NULL) 
    descrip_re = str_re_new("[[:space:]]*>[[:space:]]*(.+)");

  line_no=1;
  while ((str_readline(line, F)) != EOF) {
    if (str_re_match(line, descrip_re, l, 1) > 0) {
      lst_push_ptr(names, lst_get_ptr(l, 1));
      str_free((String*)lst_get_ptr(l, 0));

      new_str = str_new(STR_MED_LEN);
      lst_push_ptr(seqs, new_str);
      continue;
    }

    str_double_trim(line);
    if (line->length == 0) continue;

    if (new_str == NULL) 
      die("ERROR in FASTA file: non-blank line preceding first description ('>') line.\n");

    str_append(new_str, line);
    checkInterruptN(line_no++, 1000);
  }

  if (lst_size(seqs) == 0)
    die("ERROR: empty FASTA file.\n");


  // now create MS 
  nseqs = lst_size(names);
  if (nseqs != lst_size(seqs))
    die("ERROR ms_read: nseqs (%i) != lst_size(seqs) (%i)\n",
        nseqs, lst_size(seqs));

  ms = ms_new(NULL, NULL, nseqs, NULL,  0, 1); 
  ms->names = (char**)smalloc(nseqs * sizeof(char*));
  ms->seqs = (char**)smalloc(nseqs * sizeof(char*));
  ms->idx_offsets = (int*)smalloc(nseqs * sizeof(int));

  // upcase chars unless there are lowercase characters in the alphabet 
  do_toupper = !ms_alph_has_lowercase(ms);    

  for (i = 0; i < nseqs; i++) {
    n = (String*)lst_get_ptr(names, i);
    ms->names[i] = (char*)smalloc((n->length + 1) * sizeof(char));
    strcpy(ms->names[i], n->chars);
    str_free(n);

    s = (String*)lst_get_ptr(seqs, i);
    ms->seqs[i] = (char*)smalloc((s->length + 1) * sizeof(char));

	ms->idx_offsets[i] = 0;
	
    // scan chars and adjust if necessary 
    for (j = 0; j < s->length; j++) {
      ms->seqs[i][j] = do_toupper ? toupper(s->chars[j]) : s->chars[j];
      if (ms->seqs[i][j] == '.' && ms->inv_alphabet[(int)'.'] == -1) 
        ms->seqs[i][j] = ms->missing[0]; // interpret '.' as missing
      //   data; maybe no longer
      //   necessary 
      if (isalpha(ms->seqs[i][j]) && ms->inv_alphabet[(int)ms->seqs[i][j]] == -1 && get_iupac_map()[(int)ms->seqs[i][j]] == NULL) 
        ms->seqs[i][j] = 'N';   // assume 'N' if unrecognized letter
    }
    ms->seqs[i][s->length] = '\0';

    str_free(s);
    //    str_free(n);
  }  

  lst_free(names);
  lst_free(seqs);
  lst_free(l);
  str_free(line);
  phast_fclose(F);
  return ms;
}

///////////////////////////////////////////////////////////////////////////////
Matrix *mm_build_helper(MS *inputMS, int norder, int pseudoCount, int considerReverse) {
  int alph_size, i, j, ignore, tup_idx, l, alph_idx, seqLen;
  double sum = 0, val;
  Vector *freqs = vec_new(int_pow(4, norder+1));
  Matrix *mm; 
  char c;
	
  if(inputMS == NULL)
    die("ERROR: GC% group passed to mm_build_helper was null");
  if (inputMS->nseqs <= 0) //Must have at least one sequence
    die("ERROR: At least one sequence must be present to build a markov model");
  if (norder < 0)//Order of Markov matrix must be positive
    die("Order of markov model to create must be zero or greater");
  
  vec_zero(freqs);
  
  alph_size = strlen(inputMS->alphabet);
  
  //Apply Pseudo-Counts
  vec_set_all(freqs, pseudoCount);
  
  //For each sequence
  for (j = 0; j < inputMS->nseqs; j++) { 
    seqLen = strlen(inputMS->seqs[j]);
    //For each site
    for (i = 0; i < seqLen; i++) { 
      checkInterruptN(i, 10000);
      ignore = 0;
      tup_idx = 0;
      //For each base in the tuple
      for (l = 0; !ignore && l <= norder; l++) {
        
        c = inputMS->seqs[j][i + l];
        if ((alph_idx = inputMS->inv_alphabet[(int)c]) == -1) //If we get an unknown base
          ignore = 1;
        else
          tup_idx += alph_idx * int_pow(alph_size, (norder - l));
      }
      
      if (!ignore)
        vec_set(freqs, tup_idx,
                vec_get(freqs, tup_idx) + 1);
    }
  }
  
  //Take into account reverse complement frequencies
  if(considerReverse == 1)
	{
      //For each sequence
      for (j = 0; j < inputMS->nseqs; j++) { 
        
        //For each site
        for (i = strlen(inputMS->seqs[j]); i >= 0 ; i--) { 
          checkInterruptN(i, 10000);
          ignore = 0;
          tup_idx = 0;
          //For each base in the tuple
          for (l = 0; !ignore && l <= norder; l++) {
            c = inputMS->seqs[j][i - l];
            switch(c)
              {
              case 'A': 
                c = 'T';
                break;
              case 'C':
                c = 'G';
                break;
              case 'G':
                c = 'C';
                break;
              case 'T':
                c = 'A';
                break;
              }
            if ((alph_idx = inputMS->inv_alphabet[(int)c]) == -1) //If we get an unknown base
              ignore = 1;
            else
              tup_idx += alph_idx * int_pow(alph_size, (norder - l));
          }
          
          if (!ignore)
            vec_set(freqs, tup_idx,
                    vec_get(freqs, tup_idx) + 1);
        }
      }
    }
  mm = mat_new(int_pow(alph_size, norder), alph_size);
  //Transform count vector into Markov Matrix of order norder
  for (i = 0; i < freqs->size; i = i + alph_size) {
    sum = 0;
    for (j = 0; j < alph_size; j++) //Calculate sum i.e. for AA = count(AAA) + count(AAC) + count(AAG) + count(AAT)
      sum += vec_get(freqs, i + j);
    for (j = 0; j < alph_size; j++) { //For each base in alphabet
      if (sum == 0) //Handle unknown <prefix, base> tuples
        val = 1.0/alph_size;
      else
        val = vec_get(freqs, i + j) / sum; //i.e. for AAT, AA identified by i, T defined by j; #AAT/#AA
      if((val < 0) || (val > 1)) //Value should be a probability between 0 and 1
        die("ERROR: Generating Markov Models, generated probability must be between 0 and 1");
      mat_set(mm, i / alph_size, j, val);
    }
  }
  
  return mm;
}

////////////////////////////////////////////////////////////////////
List *mm_build(MS *inputMS, int norder, int pseudoCount, int considerReverse) {
  int i;
  Matrix *mm = NULL;
  List *MatrixList; 
	
  //testBaseToRow();

  if (norder < 0) //Must have a positive order to build markov Model
    die("ERROR: Order of Markov Models must be zero or greater");
	
  MatrixList = lst_new_ptr(norder+1);

  //Build a Markov Model (list of Matrix order 0 -> norder)
  for (i = 0; i <= norder; i++) {
    mm = mm_build_helper(inputMS, i, pseudoCount, considerReverse); //Build MarkovMatrix of order i
    lst_push_ptr(MatrixList, mm);
  }

  return MatrixList;
}

//////////////////////////////////////////////////////
int basetocol(char base) {
  //Put next base from sequence in opening at end of array (opened by shift above)
  switch (base) {
  case 'A': return 0;
  case 'C': return 1;
  case 'G': return 2;
  case 'T': return 3;
  case 'N': return -1;
  default: fprintf(stderr, "Encountered unknown base %c\n", base); return -1;
  }
}

//////////////////////////////////////////////////
void testBaseToRow() {
  int c,d,e;
  char bases[4] = {'A','C','G','T'};
  int previousBases[3];
  for(c=0; c <= 3; c++)
    {
      previousBases[0] = basetocol(bases[c]);
      for(d=0; d <= 3; d++)
        {
          previousBases[1] = basetocol(bases[d]);
          for(e=0; e <= 3; e++)
            {
              previousBases[2] = basetocol(bases[e]);
              printf("%c%c%c = ", bases[c],bases[d],bases[e]);
              printf("%d\n", basesToRow(previousBases, 3, 4));
            }
        }
    }

}

////////////////////////////////////////////
char *dtoa(double num)
{
  char *str = smalloc(STR_MED_LEN);
  sprintf(str, "%f", num);
  return str;
}

///////////////////////////////////////////
Matrix *mat_reverse_complement(Matrix *m) {
  int i;
  Matrix *result;
  if (m == NULL) //Matrix must be initialized
    die("ERROR: PWM to reverse complement is NULL");
  result = mat_new(m->nrows, m->ncols);
  //Reverse complement
  for(i = 0; i < m->nrows; i++) 
    {
      result->data[(m->nrows-1)-i][0] = (double)(m->data[i][3]);
      result->data[(m->nrows-1)-i][2] = (double)(m->data[i][1]);
      result->data[(m->nrows-1)-i][1] = (double)(m->data[i][2]);
      result->data[(m->nrows-1)-i][3] = (double)(m->data[i][0]);
    }
  return result;
} 

////////////////////////////////////////////////////////////////////
int basesToRow(int *previousBases, int norder, int alph_size) {
  int col = 0;
  int i;
	
  if(norder < 0) //Order must be positive
    die("Order of Markov Matrix must be zero or greater");
  if (alph_size <= 0) //Alphabet size must be at least 1
    die("Alphabet size must be at least 1");
	
  for (i = norder-1; i >= 0; i--)
    col += int_pow(alph_size, i) * (previousBases[norder-(i+1)]);

  return col;
}

/////////////////////////////////////////////////////////////////////////
char *ms_simulate(List *mModel, int norder, int alph_size, int length) {
  if (norder < 0) //Must have positive order of Markov Matrices
    die("Order of Markov Model must be zero or greater");
  if(alph_size < 1) //Must have at least one character in alphabet
    die("Alphabet size must be at least 1");
  if (length <= 0) //The length of the sequence to generate must be positive
    die("Length of sequence to generate must be at least 1");
  //printf("Called mm_simulate_seq\n");
  int currentMMorder, i, j, k, base, a = 0, t = 0, g = 0, c = 0;
  char *result = (char*)smalloc((length + 1) * sizeof(char));
  int *previousBases = smalloc((norder +1) * sizeof(int));
  double probability, r;
  Matrix *mm;


  //For length of the simulated sequence
  for (base = 0; base < length; base++) {
    checkInterruptN(base, 1000);
    if (base >= norder) 	//If at a site with less than norder previous bases, use a lower order Markov Matrix
      currentMMorder = norder;
    else
      currentMMorder = base;

    mm = (Matrix*)lst_get_ptr(mModel, currentMMorder);

    j = basesToRow(previousBases, currentMMorder, alph_size); //Map previous bases to row i.e. "AAA"->0, "AAG"->2, "TCA"->53

    r = unif_rand();  //Get random double 0..1

    //Using frequencies of Markov Matrix as weights, determine which base is next in the sequence
    for (i = 0; i < mm->ncols; i++) {
      probability = mat_get(mm, j, i);
      if ((probability < 0) || (probability > 1)) //Probabilities should be between 0 and 1
        die("ERROR: Simulating sequence, probability must be between 0 and 1");
      r = r - probability;
      if (r <= 0)
        break;
    }
		
    if (i >= mm->ncols)  //Counters the final increment of i in previous for loop if it chooses 'T'
      i = mm->ncols - 1;

    //Shift characters in sliding window of previous bases left by 1 and append newly picked base to list
    if( currentMMorder == norder) {
      for (k = 0; k < norder; k++)
        previousBases[k] = previousBases[k + 1];
      previousBases[norder-1] = i;
    }
    else{
      previousBases[currentMMorder] = i;
    }
    switch (i) {  //Convert the number we picked to its alpha equivilent
    case 0: result[base] = 'A'; break;
    case 1: result[base] = 'C'; break;
    case 2: result[base] = 'G'; break;
    case 3: result[base] = 'T'; break;
		
    }
    switch (i) {
    case 0: a++; break;
    case 1: c++; break;
    case 2: g++; break;
    case 3: t++; break;
    }
  }
  //printf("a=%d, c=%d, g=%d, t=%d\n", a, c, g, t);
  //printf("\n");
  result[base] = '\0'; //Add an end to the newly generated sequence string
  return result;
}


double calcMMscore(char *seqData, int base, List *MarkovMatrices, int conservative) {
  int i, baseAsNum, j;
  double val;
  int mmOrder = lst_size(MarkovMatrices)-1;
  Matrix *mm;
  int previousMMbases[mmOrder];
    
  //If there aren't mmOrder previous bases @ base, then adjust mmOrder to take advantage of however many we have
  if (base < mmOrder)
    mmOrder = base;
      
  //If we run into any unknown "N" characters, adjust the mmOrder accordingly
  for(i=mmOrder; i>0; i--)
    {
      baseAsNum = basetocol(seqData[base-i]);
      if (baseAsNum < 0)
        mmOrder = i-1;
      else
        previousMMbases[mmOrder-i] = baseAsNum;
    }
   	
  //Get score from Markov Matrix
  mm =  lst_get_ptr(MarkovMatrices, mmOrder);
  j = basesToRow(previousMMbases, mmOrder, mm->ncols);
  if (j >= 0)
    val = log(mat_get(mm, j, basetocol(seqData[base])));
  else
	{
      if (conservative == 1)
        val = log(0);	//If it is an unknown base, probability is 0, in log space =inf
      else
        val = 0; //If it is an unknown base probability is 1, in log space log(1)=0
	}
  return val;
}

//////////////////////////////////////////////////////////////////////////////////
GFF_Set *ms_score(char *seqName, char *seqData, int seqLen, int seqIdxOff, int seqAlphLen, List *MarkovMatrices, Matrix *pwm, Matrix *reverseCmpPWM, int conservative, double threshold, char *strand) { 
  int i, k,j,l,col;
  double MMprob, PWMprob=0, ReversePWMprob=0;
  GFF_Set *scores = gff_new_set();
  double *MMprobs = (double*)smalloc((pwm->nrows+1) * sizeof(double));    //Sliding window of mmOrder previous MM probabilities
		
  if ((conservative != 0) && (conservative != 1))
    die("ERROR: Conserverative (boolean) value must be 0 or 1");
	
  if (seqLen < pwm->nrows)  //Check to see if the sequence is shorter than the pwm
    return scores;

  for (i = 0; i <= pwm->nrows; i++)							//Calculate MM scores from sites 0 to pwm->nrows
    if (i < seqLen)
      MMprobs[i] = calcMMscore(seqData, i, MarkovMatrices, conservative);
		
  for (i = 0; i <= seqLen-(pwm->nrows); i++) {				//For each base in the sequence
    PWMprob = 0; MMprob = 0; ReversePWMprob = 0;
		
    for (k = 0, j = i; k < pwm->nrows; k++, j++) {		//Sum PWM, ReversePWM, MM probabilities for score calculation
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
            PWMprob = log(0);			//If we get something other than the expected language (i.e. A,C,T,G) i.e. N, then our probability is -Inf
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
      for (l = 0; l < pwm->nrows; l++)		//Shift probs left to make room for next
	MMprobs[l] = MMprobs[l + 1];

      MMprobs[pwm->nrows-1] = calcMMscore(seqData, i+pwm->nrows,  //Calculate MM probability for site at (i+pwm->nrows)
                                          MarkovMatrices, conservative);
    }

    if (((PWMprob - MMprob) > threshold) && ((strcmp(strand, "+") == 0) || (strcmp(strand, "both") == 0) || ((strcmp(strand, "best") == 0) && ((PWMprob - MMprob) >= (ReversePWMprob - MMprob))))) {			//If we have a positive score add it to the list of scores
      GFF_Feature *feat = gff_new_feature(str_new_charstr(seqName), str_new_charstr(""), 
                                          str_new_charstr(""), seqIdxOff+i+1, 
                                          seqIdxOff+i+pwm->nrows, (PWMprob - MMprob), '+', 
                                          0, str_new_charstr(""), 0);
      lst_push_ptr(scores->features, feat);
    }

    if (((ReversePWMprob - MMprob) > threshold) && ((strcmp(strand, "-") == 0) || (strcmp(strand, "both") == 0) || ((strcmp(strand, "best") == 0) && ((ReversePWMprob - MMprob) > (PWMprob - MMprob))))) {
      GFF_Feature *feat = gff_new_feature(str_new_charstr(seqName), str_new_charstr(""), 
                                          str_new_charstr(""), seqIdxOff+i+1, 
                                          seqIdxOff+i+pwm->nrows, (ReversePWMprob - MMprob), '-', 
                                          0, str_new_charstr(""), 0);
      lst_push_ptr(scores->features, feat);
    }
  }
  sfree(MMprobs);
  return scores; 
}


Vector *ms_gc_content(MS *ms) {
  Vector *rv = vec_new(ms->nseqs);
  int gcCount, atCount, acgtCount, i, j;
  for (i = 0; i < ms->nseqs; i++) {
    //GC calc
    gcCount = 0;
    atCount = 0;
    //Caclulate GC content
    for (j = 0; j < strlen(ms->seqs[i]); j++)
      {
        if ((ms->seqs[i][j] == 'G') || (ms->seqs[i][j] == 'C'))
          gcCount = gcCount + 1;
        else if ((ms->seqs[i][j] == 'A') || (ms->seqs[i][j] == 'T'))
          atCount = atCount + 1;
      }
    //fprintf(stderr, " gc found %d, gc score %f\n", gcCount, ((double)gcCount / (double)ms->length));

    acgtCount = atCount + gcCount;
    vec_set(rv, i, ((double)gcCount/(double)acgtCount));
  }
  return rv;
}


