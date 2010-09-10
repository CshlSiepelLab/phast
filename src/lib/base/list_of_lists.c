#include <stdlib.h>
#include <misc.h>
#include <tree_likelihoods.h>
#include <list_of_lists.h>

//make a new list of lists (LOL)
ListOfLists *lol_new(int approx_size) {
  ListOfLists *lol = smalloc(sizeof(ListOfLists));
  lol->lst = lst_new_ptr(approx_size);
  lol->lstName = lst_new_ptr(approx_size);
  lol->lstType = lst_new_int(approx_size);
  lol->class = NULL;
  return lol;
}


void lol_set_class(ListOfLists *lol, char *class) {
  if (lol->class != NULL) {
    if (strcmp(lol->class, class)==0) return;
    phast_warning("warning: changing class of list from %s to %s\n", 
		  lol->class, class);
    free(lol->class);
  }
  lol->class = copy_charstr(class);
}


//adds a pointer to list to end of LOL object.
void lol_push(ListOfLists *lol, void *data, const char *name, 
	      list_element_type listType) {
  char *tempname;
  lst_push_ptr(lol->lst, (void*)data);
  if (name == NULL)
    tempname = NULL;
  else {
    tempname = smalloc((strlen(name)+1)*sizeof(char));
    strcpy(tempname, name);
  }
  lst_push_ptr(lol->lstName, (void*)tempname);
  lst_push_int(lol->lstType, listType);
}


void lol_push_list(ListOfLists *lol, List *data, const char *name,
		   list_element_type listType) {
  lol_push(lol, (void*)data, name, listType);
}

void lol_push_lol(ListOfLists *lol, ListOfLists *data,
				  const char *name) {
  lol_push(lol, (void*)data, name, LIST_LIST);
}


//add a list of doubles to end of LOL object.  Copies all values.
void lol_push_dbl(ListOfLists *lol, double *vals, int len,
		  const char *name) {
  List *lst = lst_new_dbl(len);
  int i;
  for (i=0; i<len;  i++) 
    lst_push_dbl(lst, vals[i]);
  lol_push(lol, lst, name, DBL_LIST);
}



//add a list of ints to end of LOL object.  Copies all values.
void lol_push_int(ListOfLists *lol, int *vals, int len,
		  const char *name) {
  List *lst = lst_new_int(len);
  int i;
  for (i=0; i<len;  i++) 
    lst_push_int(lst, vals[i]);
  lol_push(lol, lst, name, INT_LIST);
}


//add a list of char* to end of LOL object.  Copies all values.
void lol_push_charvec(ListOfLists *lol, char **vals, int len,
		      const char *name) {
  List *lst = lst_new_ptr(len);
  char *tempstr;
  int i;
  for (i=0; i<len; i++) {
    tempstr = smalloc((strlen(vals[i])+1)*sizeof(char));
    strcpy(tempstr, vals[i]);
    lst_push_ptr(lst, tempstr);
  }
  lol_push(lol, lst, name, CHAR_LIST);
}


void lol_push_matrix(ListOfLists *lol, Matrix *mat,
		     const char *name) {
  ListOfLists *matList =lol_new(1 + mat->nrows);
  List *rowNames = lst_new_ptr(mat->ncols);
  Vector *tmpVec;
  char *tmpstr = smalloc(10*sizeof(char));
  int i;

  for (i=0; i < mat->nrows; i++) {
    sprintf(tmpstr, "%i", i);
    tmpVec = mat_get_row(mat, i);
    lol_push_dbl(matList, tmpVec->data, tmpVec->size, tmpstr);
    vec_free(tmpVec);
  }
  free(tmpstr);
  for (i=0; i<mat->ncols; i++) {
    tmpstr = smalloc(10*sizeof(char));
    sprintf(tmpstr, "%i", i);
    lst_push_ptr(rowNames, (void*)tmpstr);
  }
  lol_push_list(matList, rowNames, "row.names", CHAR_LIST);
  lol_set_class(matList, "matrix");
  lol_push_lol(lol, matList, name);
}


void lol_push_treeModel(ListOfLists *lol, TreeModel *tm,
			const char *name) {
  char *str;
  ListOfLists *tmList = lol_new(11);
  if (tm->rate_matrix->states != NULL)
    lol_push_charvec(tmList, &tm->rate_matrix->states, 1, "alphabet");
  if (tm->backgd_freqs != NULL)
    lol_push_dbl(tmList, tm->backgd_freqs->data, tm->backgd_freqs->size,
			 "backgd");
  if (tm->rate_matrix != NULL && tm->rate_matrix->matrix != NULL) 
    lol_push_matrix(tmList, tm->rate_matrix->matrix, "rate.matrix");
  str = copy_charstr(tm_get_subst_mod_string(tm->subst_mod));
  lol_push_charvec(tmList, &str, 1, "subst.mod");
  free(str);
  if (tm->lnL != NULL_LOG_LIKELIHOOD)
    lol_push_dbl(tmList, &(tm->lnL), 1, "likelihood");
  if (tm->alpha != 0.0)
    lol_push_dbl(tmList, &(tm->alpha), 1, "alpha");
  lol_push_int(tmList, &(tm->nratecats), 1, "nratecats");
  if (tm->rK != NULL)
    lol_push_dbl(tmList, tm->rK, tm->nratecats, "rate.consts");
  if (tm->freqK != NULL)
    lol_push_dbl(tmList, tm->freqK, tm->nratecats, "rate.weights");
  if (tm->tree != NULL) {
    str = tr_to_string(tm->tree, 1);
    lol_push_charvec(tmList, &str, 1, "tree");
    free(str);
  }
  if (tm->root_leaf_id != -1)
    lol_push_int(tmList, &(tm->root_leaf_id), 1, "root.leaf");
  lol_set_class(tmList, "tm");
  lol_push_lol(lol, tmList, name);
}


void lol_push_gff(ListOfLists *lol, GFF_Set *gff, const char *name) {
  ListOfLists *gffList = lol_new(9);
  int gffLen = lst_size(gff->features);
  char **names, **src, **feature, **strand, **attribute, tempStrand[2];
  int *start, *end, *frame;
  double *score;
  int haveScore=0, haveStrand=0, haveFrame=0, haveAttribute=0, i;
  GFF_Feature *feat;

  names = smalloc(gffLen*sizeof(char*));
  src = smalloc(gffLen*sizeof(char*));
  feature = smalloc(gffLen*sizeof(char*));
  start = smalloc(gffLen*sizeof(int));
  end=smalloc(gffLen*sizeof(int));
  score = smalloc(gffLen*sizeof(double));
  strand = smalloc(gffLen*sizeof(char*));
  frame = smalloc(gffLen*sizeof(int));
  attribute = smalloc(gffLen*sizeof(char*));
  tempStrand[1]='\0';

  for (i=0; i<gffLen; i++) {
    feat = (GFF_Feature*)lst_get_ptr(gff->features, i);
    names[i] = copy_charstr(feat->seqname->chars);
    src[i] = copy_charstr(feat->source->chars);
    feature[i] = copy_charstr(feat->feature->chars);
    start[i] = feat->start;
    end[i] = feat->end;
    if (feat->score_is_null)  {
      if (haveScore != 0)
	die("ERROR lol_push_gff haveScore is not zero but score is NULL\n");
    }
    else {
      if (!(i==0 || haveScore))
	die("ERROR haveScore should be same for all features\n");
      haveScore=1;
      score[i] = feat->score;
    }
    if (feat->strand != '.') haveStrand=1;
    tempStrand[0] = feat->strand;
    strand[i] = copy_charstr(tempStrand);
    if (feat->frame == GFF_NULL_FRAME)
      frame[i] = -1;
    else {
      haveFrame = 1;
      if (feat->frame == 0)
	frame[i] = 0;
      else if (feat->frame == 1)
	frame[i] = 2;
      else if (feat->frame == 2)
	frame[i] = 1;
      else die("invalid frame %i in feature\n", feat->frame);
    }
    if (feat->attribute == NULL || feat->attribute->length==0)
      attribute[i] = copy_charstr(".");
    else {
      haveAttribute=1;
      attribute[i] = copy_charstr(feat->attribute->chars);
    }
  }

  lol_push_charvec(gffList, names, gffLen, "seqname");
  lol_push_charvec(gffList, src, gffLen, "src");
  lol_push_charvec(gffList, feature, gffLen, "feature");
  lol_push_int(gffList, start, gffLen, "start");
  lol_push_int(gffList, end, gffLen, "end");
  if (haveScore) lol_push_dbl(gffList, score, gffLen, "score");
  if (haveStrand) lol_push_charvec(gffList, strand, gffLen, "strand");
  if (haveFrame) lol_push_int(gffList, frame, gffLen, "frame");
  if (haveAttribute) lol_push_charvec(gffList, attribute, gffLen, "attribute");
  lol_set_class(gffList, "feat");
  lol_push_lol(lol, gffList, name);

  for (i=0; i<gffLen; i++) {
    free(names[i]);
    free(src[i]);
    free(feature[i]);
    free(strand[i]);
    free(attribute[i]);
  }
  free(names);
  free(src);
  free(feature);
  free(start);
  free(end);
  free(score);
  free(strand);
  free(frame);
  free(attribute);
}


//free a lol and all associated memory (if type is char free the strings)
void lol_free(ListOfLists *lol) {
  List *currlst;
  int i, j;
  list_element_type currtype;
  char *currname=NULL, *currstr;

  for (i=0; i<lst_size(lol->lst); i++) {
    currlst = (List*)lst_get_ptr(lol->lst, i);
    currtype = lst_get_int(lol->lstType, i);
    currname = (char*)lst_get_ptr(lol->lstName, i);
    if (currname != NULL)
      free(currname);
    if (currtype == LIST_LIST) 
      lol_free((ListOfLists*)currlst);
    else {
      if (currtype == CHAR_LIST) {
	for (j=0; j<lst_size(currlst); j++) {
	  currstr = (char*)lst_get_ptr(currlst, j);
	  if ((void*)currstr != NULL)
	    free(currstr);
	}
      }
      lst_free(currlst);
    }
  }
  lst_free(lol->lstType);
  lst_free(lol->lstName);
  lst_free(lol->lst);
  if (lol->class != NULL) free(lol->class);
  free(lol);
}
