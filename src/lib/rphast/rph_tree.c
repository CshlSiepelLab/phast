/*****************************************************
rph_tree.c
The RPHAST handles to functions dealing with tree
model functions from the phast package.

Alexandra Denby
Last updated: 4/22/08
*****************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <misc.h>
#include <tree_model.h>
#include <hashtable.h>

/*********DELETE ME LATER*********/
int rphast_errno;
char* rphast_errmsg;


/******************functions defined herein******************/
void rph_tm_read(char** fname, int* treeAddress, int* modAddress, int* error, char** errstr);
void rph_tr_scale(int* modAddress, int* treeAddress, double* scale, int* error, char** errstr);
void rph_tm_print(char** fname, int* address, int* error, char** errstr);
void rph_tm_free(int* modAddress, int* error, char** errstr);

/*******************************************************
rph_tm_read
reads a tree model from a file.
*******************************************************/
void rph_tm_read(char** fname, int* treeAddress, int* modAddress, int* error, char** errstr){
  TreeNode* tree;
  TreeModel* mod;
  FILE* file;

  mod=tm_new_from_file(file=fopen_fname(fname[0], "r"));
  tree=mod->tree;

  treeAddress[0]=(unsigned int)tree;
  modAddress[0]=(unsigned int)mod;
  fclose(file);

  *error=rphast_errno;
  *errstr=rphast_errmsg;

}

/*******************************************************
rph_tr_scale
Takes the address of a tree model, and returns a new
address of a new tree model with tree scaled by "scale"
*******************************************************/
void rph_tr_scale(int* modAddress, int* treeAddress, double* scale, int* error, char** errstr){
  TreeModel *origModel, *scaledModel;
  TreeNode *tree;
  double scale_factor=scale[0];
  
  origModel=(TreeModel*)modAddress[0];
  scaledModel=tm_create_copy(origModel);
  tr_scale(scaledModel->tree, scale_factor);
  tree=scaledModel->tree;

  modAddress[0]=(unsigned int)scaledModel;
  treeAddress[0]=(unsigned int)tree;

  *error=rphast_errno;
  *errstr=rphast_errmsg;

}

/*******************************************************
rph_tm_print
Prints a tree model to a file
*******************************************************/
void rph_tm_print(char** fname, int* address, int* error, char** errstr){
  FILE* outfile=fopen_fname(fname[0],"w");
  TreeModel* mod;

  if (outfile==NULL){
    error[0]=-1;
  }

  mod=(TreeModel*)address[0];

  tm_print(outfile, mod);
  fclose(outfile);

  *error=rphast_errno;
  *errstr=rphast_errmsg;

}

/*******************************************************
rph_tm_free
Frees a tree model
*******************************************************/
void rph_tm_free(int* modAddress, int* error, char** errstr){
  TreeModel *mod;

  mod=(TreeModel*)modAddress[0];
  tm_free(mod);

  *error=rphast_errno;
  *errstr=rphast_errmsg;

}

