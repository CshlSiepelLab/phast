/*****************************************************
rph_tree_doctor.c
The RPHAST handles to function from tree_doctor, a part
of the phast package.

Alexandra Denby
4/12/08
*****************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <misc.h>
#include <tree_model.h>
#include <hashtable.h>


void rph_tm_read(char** fname, int* treeAddress, int* modAddress, int* error){
  String* suffix;
  TreeNode* tree;
  TreeModel* mod;

  mod=tm_new_from_file(fopen_fname(fname[0], "r"));
  tree=mod->tree;

  treeAddress[0]=(unsigned int)tree;
  modAddress[0]=(unsigned int)mod;

}

void rph_tr_scale(int* modAddress, int* treeAddress, double* scale, int* error){
  TreeModel *origModel, *scaledModel;
  TreeNode *tree;
  double scale_factor=scale[0];
  
  origModel=(TreeModel*)modAddress[0];
  scaledModel=tm_create_copy(origModel);
  tr_scale(scaledModel->tree, scale_factor);
  tree=scaledModel->tree;

  modAddress[0]=(unsigned int)scaledModel;
  treeAddress[0]=(unsigned int)tree;
}

void rph_tm_print(char** fname, int* address, int* error){
  FILE* outfile=fopen_fname(fname[0],"w");
  TreeModel* mod;

  if (outfile==NULL){
    error[0]=-1;
  }

  mod=(TreeModel*)address[0];

  tm_print(outfile, mod);
}

void rph_tm_free(int* modAddress){
  TreeModel *mod;

  mod=(TreeModel*)modAddress[0];
  tm_free(mod);
}

