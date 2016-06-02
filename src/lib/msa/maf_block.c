 /***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#include <misc.h>
#include <sufficient_stats.h>
#include <msa.h>
#include <maf_block.h>
#include <hashtable.h>
#include <ctype.h>
#include <assert.h>

MafBlock *mafBlock_new() {
  MafBlock *block = smalloc(sizeof(MafBlock));
  block->aLine = NULL;
  block->specMap = hsh_new(100);
  block->seqlen = -1;
  block->data = lst_new_ptr(20);
  block->prev = block->next = NULL;
  return block;
}

//allocates new subBlock and initializes
MafSubBlock *mafBlock_new_subBlock()  {
  MafSubBlock *sub = smalloc(sizeof(MafSubBlock));
  sub->seq = NULL;
  sub->quality = NULL;
  sub->src = NULL;
  sub->specName = NULL;
  sub->numLine = 0;
  return sub;
}


MafSubBlock *mafSubBlock_copy(MafSubBlock *src) {
  MafSubBlock *sub = smalloc(sizeof(MafSubBlock));
  int i, j;
  if (src->seq == NULL) sub->seq = NULL;
  else sub->seq = str_new_charstr(src->seq->chars);
  if (src->src == NULL) sub->src = NULL;
  else sub->src = str_new_charstr(src->src->chars);
  if (src->specName == NULL) sub->specName = NULL;
  else sub->specName = str_new_charstr(src->specName->chars);
  sub->start = src->start;
  sub->size = src->size;
  sub->strand = src->strand;
  sub->srcSize = src->srcSize;
  sub->numLine = src->numLine;
  for (i=0; i<src->numLine; i++) {
    sub->lineType[i] = src->lineType[i];
    if (src->lineType[i]=='i') {
      for (j=0; j<2; j++) {
	sub->iStatus[j] = src->iStatus[j];
	sub->iCount[j]  = src->iCount[j];
      }
    }
    if (src->lineType[i]=='e') 
      sub->eStatus = src->eStatus;
  }
  if (src->quality == NULL) sub->quality = NULL;
  else sub->quality = str_new_charstr(src->quality->chars);
  return sub;
}

MafBlock* mafBlock_copy(MafBlock *src) {
  MafBlock *block = smalloc(sizeof(MafBlock));
  MafSubBlock *sub;
  int i;
  if (src->aLine == NULL) block->aLine = NULL;
  else block->aLine = str_new_charstr(src->aLine->chars);
  if (src->specMap == NULL) block->specMap = NULL;
  else block->specMap = hsh_copy(src->specMap);
  block->seqlen = src->seqlen;
  if (src->data==NULL) block->data = NULL;
  else {
    block->data = lst_new_ptr(lst_size(src->data));
    for (i=0; i<lst_size(src->data); i++) {
      sub = mafSubBlock_copy((MafSubBlock*)lst_get_ptr(src->data, i));
      lst_push_ptr(block->data, (void*)sub);
    }
  }
  return block;
}

//parses a line from maf block starting with 'e' or 's' and returns a new MafSubBlock 
//object. 
MafSubBlock *mafBlock_get_subBlock(String *line) {
  int i;
  List *l = lst_new_ptr(7);
  String *str;
  MafSubBlock *sub;

  if (7 != str_split(line, NULL, l)) 
    die("Error: mafBlock_get_subBlock expected seven fields in MAF line starting "
	"with %s\n",
	((String*)lst_get_ptr(l, 0))->chars);
  
  sub = mafBlock_new_subBlock();
  
  //field 0: should be 's' or 'e'
  str = (String*)lst_get_ptr(l, 0);
  if (str_compare_charstr(str, "s")==0)
    sub->lineType[0]='s';
  else if (str_compare_charstr(str, "e")==0)
    sub->lineType[0]='e';
  else die("ERROR: mafBlock_get_subBlock expected first field 's' or 'e' (got %s)\n",
	   str->chars);

  //field 1: should be src.  Also set specName
  sub->src = (String*)lst_get_ptr(l, 1);
  sub->specName = str_new_charstr(sub->src->chars);
  str_shortest_root(sub->specName, '.');

  //field 2: should be start
  sub->start = atol(((String*)lst_get_ptr(l, 2))->chars);
  
  //field 3: should be length
  sub->size = atoi(((String*)lst_get_ptr(l, 3))->chars);

  //field 4: should be strand
  str = (String*)lst_get_ptr(l, 4);
  if (str_compare_charstr(str, "+")==0)
    sub->strand = '+';
  else if (str_compare_charstr(str, "-")==0)
    sub->strand = '-';
  else die("ERROR: got strand %s\n", str->chars);
  
  //field 5: should be srcSize
  sub->srcSize = atol(((String*)lst_get_ptr(l, 5))->chars);

  //field 6: sequence if sLine, eStatus if eLine.
  str = (String*)lst_get_ptr(l, 6);
  if (sub->lineType[0]=='s')
    sub->seq = str;
  else {
    if (sub->lineType[0] != 'e')
      die("ERROR mafBlock_get_subBlock: bad lineType (expected 'e', got %c)\n",
	  sub->lineType[0]);
    if (str->length != 1)
      die("ERROR: e-Line with status %s in MAF block\n", str->chars);
    sub->eStatus = str->chars[0];
    //note: don't know what status 'T' means (it's not in MAF documentation), but
    //it is in the 44-way MAFs
    if (sub->eStatus != 'C' && sub->eStatus != 'I' && sub->eStatus != 'M' &&
	sub->eStatus != 'n' && sub->eStatus != 'T')
      die("ERROR: e-Line has illegal status %c\n", sub->eStatus);
  }
  sub->numLine = 1;
  //free all strings except field 1 and field 6 when lineType=='s'
  for (i=0; i<7; i++)
    if (i!=1 && (i!=6 || sub->lineType[0]!='s')) 
      str_free((String*)lst_get_ptr(l, i));
  lst_free(l);
  return sub;
}

void mafBlock_add_iLine(String *line, MafSubBlock *sub) {
  List *l = lst_new_ptr(6);
  String *str;
  int i;

  if (sub->numLine<1 || sub->lineType[0]!='s') 
    die("ERROR: got i-Line without preceding s-Line in MAF block\n");
  
  if (6 != str_split(line, NULL, l))
    die("ERROR: expected six fields in MAF line starting with 'i' (got %i)\n",
	lst_size(l));

  //field[0] should be 'i'
  if (!(str_compare_charstr((String*)lst_get_ptr(l, 0), "i")==0))
    die("ERROR: mafBlock_add_iLine: field[0] should be 'i', got %s\n",
	((String*)lst_get_ptr(l, 0))->chars);

  //field[1] should be src, and should match src already set in sub
  if (str_compare((String*)lst_get_ptr(l, 1), sub->src) != 0)
    die("iLine sourceName does not match preceding s-Line (%s, %s)\n", 
	((String*)lst_get_ptr(l, 1))->chars, sub->src->chars);

  for (i=0; i<2; i++) {

    //field[2,4] should be leftStatus, rightStauts
    str = (String*)lst_get_ptr(l, i*2+2);
    if (str->length != 1) die("ERROR: i-Line got illegal %sStatus = %s\n",
			      i==0 ? "left": "right", str->chars);
    sub->iStatus[i] = str->chars[0];
    if (sub->iStatus[i] != 'C' && sub->iStatus[i] != 'I' &&
	sub->iStatus[i] != 'N' && sub->iStatus[i] != 'n' &&
	sub->iStatus[i] != 'M' && sub->iStatus[i] != 'T')
      die("ERROR: i-Line got illegal %sStatus = '%c'\n",
	  i==0 ? "left" : "right", sub->iStatus[i]);

    //field 3,5 should be leftCount, rightCount
    str = (String*)lst_get_ptr(l, i*2+3);
    sub->iCount[i] = atoi(str->chars);
  }
  
  for (i=0; i<6; i++) str_free((String*)lst_get_ptr(l, i));
  lst_free(l);
  if (sub->numLine >= 4) die("Error: bad MAF file");
  sub->lineType[sub->numLine++] = 'i';


}


void mafBlock_add_qLine(String *line, MafSubBlock *sub) {
  List *l = lst_new_ptr(3);
  String *str;
  int i;

  if (sub->numLine<1 || sub->lineType[0]!='s') 
    die("ERROR: got q-Line without preceding s-Line in MAF block\n");

  if (3 != str_split(line, NULL, l))
    die("ERROR: expected three fields in q-Line of maf file, got %i\n", lst_size(l));
  
  //field[0] should be 'q'
  if (!(str_compare_charstr((String*)lst_get_ptr(l, 0), "q")==0))
    die("ERROR mafBlock_add_qLine expected 'q' got %s\n",
	((String*)lst_get_ptr(l, 0))->chars);
  
  //field[1] should be src, and should match src already set in sub
  if (str_compare((String*)lst_get_ptr(l, 1), sub->src) != 0)
    die("iLine sourceName does not match preceding s-Line (%s, %s)\n", 
	((String*)lst_get_ptr(l, 1))->chars, sub->src->chars);

  //field[2] should be quality
  if (sub->seq == NULL)
    die("ERROR mafBlock_add_qLine: sub->seq is NULL\n");
  str = (String*)lst_get_ptr(l, 2);
  if (sub->seq->length != str->length) 
    die("ERROR: length of q-line does not match sequence length\n");
  sub->quality = str;
  for (i=0; i<sub->quality->length; i++) {
    if (sub->seq->chars[i] == '-') {
      if (sub->quality->chars[i] != '-') 
	die("ERROR: got quality score where alignment char is gap\n");
    } else {
      if (sub->quality->chars[i] != 'F' && sub->quality->chars[i] < '0' &&
	  sub->quality->chars[i] > '9')
	die("ERROR: Illegal quality score '%c' in MAF block\n", 
	    sub->quality->chars[i]);
    }
  }
   
  for (i=0; i<2; i++) str_free((String*)lst_get_ptr(l, i));
  lst_free(l);
  if (sub->numLine >= 4) die("Error: bad MAF file");
  sub->lineType[sub->numLine++] = 'q';

}


//read next block in mfile and return MafBlock object or NULL if EOF.
//specHash and numSpec are not used, but if specHash is not NULL,
//it should be initialized, and any new species encountered will be added
//to the hash, with numSpec increased accordingly.  If specHash is NULL,
//numSpec will not be used or modified.
MafBlock *mafBlock_read_next(FILE *mfile, Hashtable *specHash, int *numSpec) {
  int i;
  char firstchar;
  String *currLine = str_new(1000);
  MafBlock *block=NULL;
  MafSubBlock *sub=NULL;

  if (specHash != NULL && numSpec==NULL) 
    die("ERROR: mafBlock_read_next: numSpec cannot be NULL "
	"if specHash is not NULL\n");

  while (EOF != str_readline(currLine, mfile)) {
    str_trim(currLine);
    if (currLine->length==0) {  //if blank line, it is either first or last line
      if (block == NULL) continue;
      else break;
    }
    firstchar = currLine->chars[0];
    if (firstchar == '#') continue;  //ignore comments
    if (block == NULL) {
      if (firstchar != 'a') 
	die("ERROR: first line of MAF block should start with 'a'\n");
      block = mafBlock_new();
      block->aLine = str_new_charstr(currLine->chars);
    }
    //if 's' or 'e', then this is first line of data for this species
    else if (firstchar == 's' || firstchar == 'e') {
      sub = mafBlock_get_subBlock(currLine);
      if (hsh_get_int(block->specMap, sub->src->chars) != -1) 
	die("ERROR: mafBlock has two alignments with same srcName (%s)\n", 
	    sub->src->chars);
      hsh_put_int(block->specMap, sub->src->chars, lst_size(block->data));
      hsh_put_int(block->specMap, sub->specName->chars, lst_size(block->data));
      lst_push_ptr(block->data, (void*)sub);
      if (specHash != NULL) {
	if (-1 == hsh_get_int(specHash, sub->specName->chars)) {
	  hsh_put_int(specHash, sub->specName->chars, *numSpec);
	  (*numSpec)++;
	}
      }
    }
    else {
      if (firstchar == 'i')
	mafBlock_add_iLine(currLine, sub);
      else if (firstchar == 'q')
	mafBlock_add_qLine(currLine, sub);
      else die("ERROR: found line in MAF block starting with '%c'\n", firstchar);
    }
  }
  str_free(currLine);
  if (block == NULL) return NULL;

  //set seqlen and make sure all seq arrays agree
  for (i=0; i<lst_size(block->data); i++) {
    sub = (MafSubBlock*)lst_get_ptr(block->data, i);
    if (sub->lineType[0]=='e') continue;
    if (block->seqlen == -1) block->seqlen = sub->seq->length;
    else if (sub->seq->length != block->seqlen) {
      die("ERROR: lengths of sequences in MAF block do not agree (%i, %i)\n",
	  block->seqlen, sub->seq->length);
    }
  }
  return block;
}

//returns 1 if the block entirely consists of gaps.
int mafBlock_all_gaps(MafBlock *block) {
  MafSubBlock *sub;
  int i, j;
  for (i=0; i<lst_size(block->data); i++) {
    sub = (MafSubBlock*)lst_get_ptr(block->data, i);
    if (sub->lineType[0]=='e') continue;
    for (j=0; j<block->seqlen; j++)
      if (sub->seq->chars[j] != '-') return 0;
  }
  return 1;
}


/* Checks for columns which are only gaps (may happen after processing 
   such as removal of certain species or indel masking).  Deletes these
   columns from the block.  If all bases in a species are gaps, then
   convert to e-line.
 */
void mafBlock_remove_gap_cols(MafBlock *block) {
  int *isGap = NULL, nonGapSpecies=0, i, j, k, gapCount, start, pos=-1;
  MafSubBlock *sub;

  for (i=0; i<lst_size(block->data); i++) {
    sub = (MafSubBlock*)lst_get_ptr(block->data, i);
    if (sub->lineType[0]=='e') continue;
    if (isGap==NULL) {
      isGap = smalloc(block->seqlen*sizeof(int));
      for (j=0; j<block->seqlen; j++) isGap[j]=1;
    }
    gapCount = 0;
    for (j=0; j<block->seqlen; j++) {
      if (sub->seq->chars[j] != '-') 
	isGap[j]=0;
      else gapCount++;
    }
    if (gapCount == block->seqlen) {
      //in this case, convert to e-line
      for (k=1; k<sub->numLine; k++) {
	if (sub->lineType[k] == 'q') {
	  str_free(sub->quality);
	  sub->quality = NULL;
	}
      }
      sub->numLine = 1;
      sub->lineType[0] = 'e';
      str_free(sub->seq);
      sub->seq = NULL;
      sub->eStatus = 'C';  //FIXME? this seems like the right status in the case
                           //that bases were deleted and turned to gaps,
                           //but not completely sure.
    }
    else nonGapSpecies++;
  }
  if (isGap==NULL) return;

  for (start=0; start<block->seqlen; start++)
    if (isGap[start]) break;

  if (nonGapSpecies > 0 && start != block->seqlen) {
    for (i=0; i<lst_size(block->data); i++) {
      sub = (MafSubBlock*)lst_get_ptr(block->data, i);
      if (sub->lineType[0] == 'e') continue;
      pos = start;
      for (j=start; j<block->seqlen; j++) {
	if (isGap[j]==0) {
	  if (pos != j)
	    sub->seq->chars[pos] = sub->seq->chars[j];
	  pos++;
	}
      }
      sub->seq->chars[pos]='\0';
      if (sub->quality != NULL) {
	pos = start;
	for (j=start; j<block->seqlen; j++) {
	  if (isGap[j] == 0) {
	    if (pos != j)
	      sub->quality->chars[pos] = sub->quality->chars[j];
	    pos++;
	  }
	}
	sub->quality->chars[pos]='\0';
      }
    }
    if (pos <= 0)
      die("ERROR mafBlock_remove_gap_cols: pos=%i, should be >0\n", pos);
    block->seqlen = pos;
  } 
  else if (nonGapSpecies==0) 
    block->seqlen = 0;
  sfree(isGap);
  return;
}


//sets fieldSize[i] to maximum length of field i in MAF, so that block
//can be printed with nice formatting
void mafBlock_get_fieldSizes(MafBlock *block, int fieldSize[6]) {
  int i;
  MafSubBlock *sub;
  char tempstr[1000];
  for (i=0; i<6; i++) fieldSize[i] = 0;

  fieldSize[0] = 1;  //this is always one character
  fieldSize[4] = 1;  //this is always one character (strand)
  for (i=0; i<lst_size(block->data); i++) {
    sub = (MafSubBlock*)lst_get_ptr(block->data, i);

    //field[1] is src
    if (sub->src->length > fieldSize[1])
      fieldSize[1] = sub->src->length;

    //field[2] is start
    sprintf(tempstr, "%li", sub->start);
    if (strlen(tempstr) > fieldSize[2])
      fieldSize[2] = (int)strlen(tempstr);

    //field[3] is size
    sprintf(tempstr, "%i", sub->size);
    if (strlen(tempstr) > fieldSize[3])
      fieldSize[3] = (int)strlen(tempstr);
    
    //field[4] is strand... skip
    
    //field[5] is srcSize
    sprintf(tempstr, "%li", sub->srcSize);
    if (strlen(tempstr) > fieldSize[5])
      fieldSize[5] = (int)strlen(tempstr);

    //don't worry about size of lastField since it just goes to end-of-line
  }
}


void mafBlock_print(FILE *outfile, MafBlock *block, int pretty_print) {
  int i, j, k, numSpace;
  int fieldSize[6];  //maximum # of characters in the first 6 fields of block
  MafSubBlock *sub;
  char firstChar, formatstr[1000];
  char *firstseq=NULL;

  //if processing has reduced the number of species with data to zero, or has
  //reduced the block to all gaps, don't print
  if (lst_size(block->data) == 0 ||
      mafBlock_all_gaps(block)) return;
  mafBlock_remove_gap_cols(block);
  mafBlock_get_fieldSizes(block, fieldSize);

  fprintf(outfile, "%s\n", block->aLine->chars);
  for (i=0; i<lst_size(block->data); i++) {
    sub = (MafSubBlock*)lst_get_ptr(block->data, i);
    for (j=0; j<sub->numLine; j++) {
      firstChar = sub->lineType[j];
      if (firstChar == 's' || firstChar == 'e') {
	sprintf(formatstr, "%%c %%-%is %%%ii %%%ii %%c %%%ii ",
		fieldSize[1], fieldSize[2], fieldSize[3], fieldSize[5]);
	fprintf(outfile, formatstr, firstChar, sub->src->chars,
		sub->start, sub->size, sub->strand, sub->srcSize);
	if (firstChar == 's') {
	  if (firstseq == NULL) {
	    fprintf(outfile, "%s\n", sub->seq->chars);
	    if (pretty_print) firstseq = sub->seq->chars;
	  }
	  else {
	    for (k=0; k<block->seqlen; k++)
	      fputc(tolower(sub->seq->chars[k])==tolower(firstseq[k]) ? 
		    '.' : sub->seq->chars[k], 
		    outfile);
	  }
	}
	else fprintf(outfile, "%c\n", sub->eStatus);
      } else if (firstChar=='i') {
	sprintf(formatstr, "i %%-%is %%c %%i %%c %%i",
		fieldSize[1]);
	fprintf(outfile, formatstr, sub->src->chars,
		sub->iStatus[0], sub->iCount[0],
		sub->iStatus[1], sub->iCount[1]);
	fputc('\n', outfile);
      } else {
	if (firstChar != 'q')
	  die("ERROR mafBlock_print: firstChar should be q, got %c\n", firstChar);
	sprintf(formatstr, "q %%-%is", fieldSize[1]);
	fprintf(outfile, formatstr, sub->src->chars);
	numSpace = 6 + fieldSize[2] + fieldSize[3] + fieldSize[5];
	for (k=0; k<numSpace; k++) fputc(' ', outfile);
	fprintf(outfile, "%s\n", sub->quality->chars);
      }
    }
  }
  fputc('\n', outfile);  //blank line to mark end of block
  //  fflush(outfile);
}

void mafSubBlock_free(MafSubBlock *sub) {
  if (sub->seq != NULL) {
    str_free(sub->seq);
    sub->seq = NULL;
  }
  if (sub->src != NULL) {
    str_free(sub->src);
    sub->src = NULL;
  }
  if (sub->specName != NULL) {
    str_free(sub->specName);
    sub->specName = NULL;
  }
  if (sub->quality != NULL) {
    str_free(sub->quality);
    sub->quality = NULL;
  }
  sfree(sub);
}


void mafBlock_remove_lines(MafBlock *block, int *keep) {
  int i, oldSize = lst_size(block->data), newSize=0;
  MafSubBlock *sub, *testSub;
  for (i=0; i<oldSize; i++) {
    if (keep[i]) {
      if (i != newSize) {
	sub = (MafSubBlock*)lst_get_ptr(block->data, i);
	hsh_reset_int(block->specMap, sub->src->chars, newSize);
	hsh_reset_int(block->specMap, sub->specName->chars, newSize);
	testSub = (MafSubBlock*)lst_get_ptr(block->data, newSize);
	if (testSub != NULL)
	  die("ERROR: mafBlock_remove_lines: testSub should be NULL\n");
	lst_set_ptr(block->data, newSize, (void*)sub);
	lst_set_ptr(block->data, i, NULL);
      }
      newSize++;
    } else {
      sub = (MafSubBlock*)lst_get_ptr(block->data, i);
      hsh_reset_int(block->specMap, sub->src->chars, -1);
      hsh_reset_int(block->specMap, sub->specName->chars, -1);
      mafSubBlock_free(sub);
      lst_set_ptr(block->data, i, NULL);
    }
  }
  for (i=oldSize-1; i>=newSize; i--)
    lst_delete_idx(block->data, i);
}

//if exclude==0, removes all species not in list.
//if exclude==1, removes all species in list
void mafBlock_subSpec(MafBlock *block, List *specNameList, int include) {
  String *str;
  int i, idx, *keep, oldSize = lst_size(block->data);

  keep = smalloc(oldSize*sizeof(int));
  for (i=0; i<oldSize; i++) keep[i]=(include==0);

  for (i=0; i<lst_size(specNameList); i++) {
    str = (String*)lst_get_ptr(specNameList, i);
    idx = hsh_get_int(block->specMap, str->chars);
    if (idx != -1) keep[idx] = !(include==0);
  }
  mafBlock_remove_lines(block, keep);
  sfree(keep);
  return;
}


void mafBlock_reorder(MafBlock *block, List *specNameOrder) {
  String *str;
  MafSubBlock *sub;
  List *newData;
  Hashtable *newSpecMap;
  int i, idx, *found, oldSize = lst_size(block->data), newSize = lst_size(specNameOrder);

  found = smalloc(oldSize*sizeof(int));
  for (i=0; i<oldSize; i++) found[i]=0;

  newData = lst_new_ptr(oldSize);
  newSpecMap = hsh_new(100);

  for (i=0; i<newSize; i++) {
    str = (String*)lst_get_ptr(specNameOrder, i);
    idx = hsh_get_int(block->specMap, str->chars);
    if (idx != -1) {
      if (found[idx]==1) die("ERROR: species %s appears twice in reorder list\n", 
			     str->chars);
      sub = (MafSubBlock*)lst_get_ptr(block->data, idx);
      hsh_put_int(newSpecMap, sub->src->chars, lst_size(newData));
      hsh_put_int(newSpecMap, sub->specName->chars, lst_size(newData));
      lst_push_ptr(newData, (void*)sub);
      found[idx] = 1;
    }
  }
  for (i=0; i<oldSize; i++) {
    if (found[i]==0) {
      sub = (MafSubBlock*)lst_get_ptr(block->data, i);
      mafSubBlock_free(sub);
    }
  }
  hsh_free(block->specMap);
  lst_free(block->data);
  block->specMap = newSpecMap;
  block->data = newData;
  sfree(found);
}


void mafBlock_free_data(MafBlock *block) {
  MafSubBlock *sub;
  int i;
  
  if (block->data != NULL) {
    for (i=0; i<lst_size(block->data); i++) {
      sub = (MafSubBlock*)lst_get_ptr(block->data, i);
      mafSubBlock_free(sub);
    }
    lst_free(block->data);
    block->data = NULL;
  }
  block->seqlen = 0;
}


//frees all elements of block except prev and next pointers
void mafBlock_free(MafBlock *block) {
  if (block->aLine != NULL) {
    str_free(block->aLine);
    block->aLine = NULL;
  }
  if (block->specMap != NULL) {
    hsh_free(block->specMap);
    block->specMap = NULL;
  }
  mafBlock_free_data(block);
  sfree(block);
}


FILE *mafBlock_open_outfile(char *fn, int argc, char *argv[]) {
  FILE *outfile;
  int i;
  if (fn != NULL) 
    outfile = phast_fopen_no_exit(fn, "w");
  else outfile = stdout;
  if (outfile == NULL) return NULL;
  fprintf(outfile, "##maf version=1\n#");
  for (i=0; i<argc; i++) 
    fprintf(outfile, " %s", argv[i]);
  fputc('\n', outfile);
  return outfile;
}

void mafBlock_close_outfile(FILE *outfile) {
  fprintf(outfile, "#eof\n");
  if (outfile != stdout) phast_fclose(outfile);
}

String *mafBlock_get_refSpec(MafBlock *block) {
  MafSubBlock *sub = (MafSubBlock*)lst_get_ptr(block->data, 0);
  if (sub==NULL) return NULL;
  return sub->specName;
}

long mafBlock_get_start(MafBlock *block, String *specName) {
  int idx=0;
  if (specName != NULL) 
    idx = hsh_get_int(block->specMap, specName->chars);
  if (idx == -1 || idx >= lst_size(block->data)) return -1;
  return ((MafSubBlock*)lst_get_ptr(block->data, idx))->start;
}

int mafBlock_get_size(MafBlock *block, String *specName) {
  int idx=0;
  MafSubBlock *sub;
  if (specName == NULL) return block->seqlen;
    idx = hsh_get_int(block->specMap, specName->chars);
  if (idx == -1 || idx >= lst_size(block->data)) return -1;
  sub = (MafSubBlock*)lst_get_ptr(block->data, idx);
  if (sub->lineType[0]=='s') return sub->size;
  if (sub->lineType[0] != 'e')
    die("ERROR mafBlock_get_size, expected line type 'e', got %c\n",
	sub->lineType[0]);
  return 0;
}

int mafBlock_numSpec(MafBlock *block) {
  return lst_size(block->data);
}

void mafSubBlock_strip_iLine(MafSubBlock *sub) {
  int i, j;
  for (i=0; i<sub->numLine; i++) 
    if (sub->lineType[i]=='i') break;
  if (i < sub->numLine) {
    assert(i < 4);
    for (j=i+1; j<sub->numLine && j<4; j++) {
      sub->lineType[j-1] = sub->lineType[j];
      if (sub->lineType[j] == 'i')
	die("ERROR mafSubBlock_strip_iLine: sub->lineType[%i]=%c\n",
	    j, sub->lineType[j]);
    }
    sub->numLine--;
  }
}


void mafBlock_subAlign(MafBlock *block, int start, int end) {
  int i, j, oldSeqlen = block->seqlen;
  MafSubBlock *sub;
  String *str;
  if (start > end || start <= 0 || start > oldSeqlen ||
      end <= 0 || end > oldSeqlen) 
    die("ERROR: mafBlock_subAlign got start=%i, end=%i, seqlen=%i\n", 
	start, end, oldSeqlen);
  if (end==oldSeqlen && start==1) return;  //nothing to do

  start--; //convert to zero-based coords
  block->seqlen = end-start;

  for (i=0; i<lst_size(block->data); i++) {
    sub = (MafSubBlock*)lst_get_ptr(block->data, i);
    if (sub->lineType[0]=='e') continue;  //e-lines remain unchanged

    if (sub->lineType[0] != 's')
      die("ERROR mafBlock_sub_align: expected lineType 's', got %c\n",
	  sub->lineType[0]);
    for (j=0; j<start; j++)
      if (sub->seq->chars[j]!='-') sub->start++;
    sub->size = 0;
    for (j=start; j<end; j++)
      if (sub->seq->chars[j] != '-') sub->size++;
    
    mafSubBlock_strip_iLine(sub);

    //get rid of i-line if exists
    mafSubBlock_strip_iLine(sub);

    //trim seq and quality scores
    str = str_new(end-start);
    str_substring(str, sub->seq, start, end-start);
    str_free(sub->seq);
    sub->seq = str;

    if (sub->quality != NULL) {
      str = str_new(end-start);
      str_substring(str, sub->quality, start, end-start);
      str_free(sub->quality);
      sub->quality = str;
    }
  }
}


//trim mafblock to only keep columns with indcies[startcol..endcol] wrt
//refseq.  If refseq is null use frame of entire alignment.  If endcol is -1
//then keep everything with index >= startcol.
int mafBlock_trim(MafBlock *block, int startcol, int endcol, String *refseq,
		  int offset) {
  MafSubBlock *sub=NULL;
  int i, specIdx, first=-1, last=-1, keep, length;
  long startIdx, lastIdx, idx;
  if (block->seqlen == 0) return 0;
  if (refseq == NULL) {
    startIdx = 1;
    length = block->seqlen;
  }
  else {
    specIdx = hsh_get_int(block->specMap, refseq->chars);
    if (specIdx == -1)
      die("Error: mafBlock_trim got specIdx -1\n");
    sub = (MafSubBlock*)lst_get_ptr(block->data, specIdx);
    startIdx = sub->start + 1;
    length = sub->size;
  }
  startIdx += offset;
  lastIdx = startIdx + length - 1;
  idx = startIdx;
  if (refseq != NULL && sub->seq->chars[0]=='-') startIdx--;

  if (startcol != 1 && endcol != -1 && startcol > endcol) 
    die("ERROR: startcol > endcol\n");
  if (startcol > lastIdx || (endcol != -1 && endcol < startIdx)) {
    mafBlock_free_data(block);
    return 0;
  }
  if (startcol <= startIdx && 
      (endcol   == -1 || endcol  >= lastIdx)) return 1;
  
  //now we know we have to do some trimming
  for (i=0; i<block->seqlen; i++) {
    if (refseq != NULL && sub->seq->chars[i]=='-') idx--;

    keep = (idx >= startcol &&
	    (idx <= endcol   || endcol == -1));
    if (first == -1 && keep) first = i+1;
    if (keep) last=i+1;
    idx++;
  }
  mafBlock_subAlign(block, first, last);
  return 1;
}


void mafBlock_strip_iLines(MafBlock *block) {
  MafSubBlock *sub;
  int i;
  for (i=0; i<lst_size(block->data); i++) {
    sub = (MafSubBlock*)lst_get_ptr(block->data, i);
    mafSubBlock_strip_iLine(sub);
  }
}


void mafBlock_strip_eLines(MafBlock *block) {
  int i, *keep = smalloc(lst_size(block->data)*sizeof(int));
  for (i=0; i<lst_size(block->data); i++) 
    keep[i] = (((MafSubBlock*)lst_get_ptr(block->data, i))->lineType[0] != 'e');
  mafBlock_remove_lines(block, keep);
  sfree(keep);
}

//strip both i- and e-lines
void mafBlock_strip_ieLines(MafBlock *block) {
  mafBlock_strip_iLines(block);
  mafBlock_strip_eLines(block);
}


void mafBlock_mask_region(MafBlock *block, GFF_Set *mask_feats, List *speclist) {
  MafSubBlock *refblock, *maskblock;
  int i, j, spec_idx;
  GFF_Set *feat;
  GFF_Feature *f, *prevf=NULL;
  int next_feat_idx = 1;
  char **maskseq;
  int num_mask_seq=0;
  long coord;
  if (mask_feats == NULL || lst_size(mask_feats->features) == 0L) return;
  maskseq = smalloc(lst_size(speclist)*sizeof(char*));
  for (i=0; i < lst_size(speclist); i++) {
    spec_idx = hsh_get_int(block->specMap, ((String*)lst_get_ptr(speclist, i))->chars);
    if (spec_idx == -1) continue;
    maskblock = lst_get_ptr(block->data, spec_idx);
    if (maskblock->seq == NULL) continue;
    maskseq[num_mask_seq++] = maskblock->seq->chars;
  }
  if (num_mask_seq == 0) {
    sfree(maskseq);
    return;
  }
  feat = gff_copy_set_no_groups(mask_feats);
  gff_flatten_mergeAll(feat);
  f = lst_get_ptr(feat->features, 0);

  refblock = lst_get_ptr(block->data, 0);
  coord = refblock->start;
  for (i=0; i < block->seqlen; i++) {
    if (refblock->seq->chars[i] != '-') coord++;  //this is 1-based coordinate
    if (coord > f->end) {
      if (next_feat_idx == lst_size(feat->features))
	break;
      prevf = f;
      f = lst_get_ptr(feat->features, next_feat_idx++);
      if (f->start <= prevf->end) {
	die("Error: feats not sorted in mafBlock_mask_region");  //shouldn't happen
      }
    }
    if (coord >= f->start && coord <= f->end) {
      for (j=0; j < num_mask_seq; j++)
	if (maskseq[j][i] != '-') maskseq[j][i] = 'N';
    }
  }
  gff_free_set(feat);
  sfree(maskseq);
}

//change all bases with quality score <= cutoff to N
void mafBlock_mask_bases(MafBlock *block, int cutoff, FILE *outfile) {
  MafSubBlock *sub;
  int i, j, firstMasked;
  char *refseq=NULL, *refseqName;
  long firstCoord, lastCoord=-1, *coord;

  sub = (MafSubBlock*)lst_get_ptr(block->data, 0);
  refseq = sub->seq->chars;
  coord = smalloc(block->seqlen*sizeof(long));
  firstCoord = sub->start;
  refseqName = sub->src->chars;
  for (i=0; i < block->seqlen; i++) {
    if (refseq[i] != '-')
      coord[i] = firstCoord++;
    else coord[i] = -1;
  }
  lastCoord = firstCoord;

  for (i=0; i<lst_size(block->data); i++) {
    sub = (MafSubBlock*)lst_get_ptr(block->data, i);
    if (sub->quality==NULL) continue;
    firstMasked=-1;
    for (j=0; j<block->seqlen; j++) {
      if (sub->quality->chars[j] == '-' && refseq[j]=='-') continue;
      if ((sub->quality->chars[j] != '-' && sub->quality->chars[j] != 'F')
	  && sub->quality->chars[j] - '0' <= cutoff) {
	sub->seq->chars[j] = 'N';
	if (firstMasked == -1 && refseq[j]!='-')
	  firstMasked = j;
      }
      else if (firstMasked != -1) {
	if (outfile != NULL) {
	  fprintf(outfile, "%s\t%li\t%li\t%s\n", 
		  refseqName, coord[firstMasked], coord[j], sub->src->chars);
	}
	firstMasked = -1;
      }
    }
    if (outfile != NULL && firstMasked != -1) 
      fprintf(outfile, "%s\t%li\t%li\t%s\n", refseqName, coord[firstMasked], lastCoord, sub->src->chars);
  }
  sfree(coord);
}

/* mask any indels that start in this block.  An indel will be "masked" if
   it is lineage-specific and has a quality score <= cutoff.  An insertion
   is masked by deleting the inserted bases (changing them to gaps).  A 
   deletion is masked by changing the gaps to N's.
   Treat indels of length > max_indel_length as missing data (don't mask).  
   Reads additional blocks from mfile as necessary to check contiguity and
   deal with indels which span multiple blocks.
   flank_size determines the number of neighbors used to determine quality
   score of indel.  If deletion, then it is the minimum score across this
   many neighbors.  If insertion, it is minimum score across insertion plus
   this many neighbors.
   WARNING:  masking indels changes the coordinate system.  The coordinates 
   will no longer be consistent with each other in the resulting block!  
   High-quality species should still be OK, since we will not mask there. */
/*void mafBlock_mask_indels(MafBlock *block, int cutoff, int flank_size, 
			  int max_indel_length, TreeNode *tree, 
			  FILE *mfile) {
  
			  }*/
