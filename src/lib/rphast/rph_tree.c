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
#include <rph_util.h>


/******************functions defined herein******************/
void rph_tm_read(char** fname, double* treeAddress, double* modAddress, int* error, char** errstr);
void rph_tr_scale(double* modAddress, double* treeAddress, double* scale, int* error, char** errstr);
void rph_tm_print(char** fname, double* address, int* error, char** errstr);
void rph_tm_free(double* modAddress, int* error, char** errstr);


/*******************************************************
rph_tm_read
reads a tree model from a file.
*******************************************************/
void rph_tm_read(char** fname, double* treeAddress, double* modAddress, int* error, char** errstr){
  TreeNode* tree;
  TreeModel* mod;
  FILE* file;

  file=fopen_fname(*fname, "r");

  mod=tm_new_from_file(file);
  tree=mod->tree;

  *treeAddress=ptr2ad(tree);
  *modAddress=ptr2ad(mod);
  fclose(file);

  *error=rphast_errno;
  *errstr=rphast_errmsg;

}

/*******************************************************
rph_tr_scale
Takes the address of a tree model, and returns a new
address of a new tree model with tree scaled by "scale"
*******************************************************/
void rph_tr_scale(double* modAddress, double* treeAddress, double* scale, int* error, char** errstr){
  TreeModel *origModel, *scaledModel;
  TreeNode *tree;
  double scale_factor=*scale;
  
  origModel=(TreeModel*)ad2ptr(*modAddress);
  scaledModel=tm_create_copy(origModel);
  tr_scale(scaledModel->tree, scale_factor);
  tree=scaledModel->tree;

  *modAddress=ptr2ad(scaledModel);
  *treeAddress=ptr2ad(tree);

  *error=rphast_errno;
  *errstr=rphast_errmsg;

}

/*******************************************************
rph_tm_print
Prints a tree model to a file
*******************************************************/
void rph_tm_print(char** fname, double* address, int* error, char** errstr){
  FILE* outfile=fopen_fname(fname[0],"w");
  TreeModel* mod;

  mod=(TreeModel*)ad2ptr(*address);

  tm_print(outfile, mod);
  fclose(outfile);

  *error=rphast_errno;
  *errstr=rphast_errmsg;

}

/*******************************************************
rph_tm_free
Frees a tree model
*******************************************************/
void rph_tm_free(double* modAddress, int* error, char** errstr){
  TreeModel *mod;

  mod=(TreeModel*)ad2ptr(*modAddress);
  tm_free(mod);

  *error=rphast_errno;
  *errstr=rphast_errmsg;

}

void rph_tm_fit(double* modAddress, double* msaAddress, char** substModel, char**treestr, int* prec, int* em, double* fittedAddress, int* error, char** errstr){

  TreeModel *init_mod, *mod;
  TreeNode* tree;
  MSA* msa;
  Vector* params;
  int precision;
  int use_em=*em;
  FILE* logf=NULL;
  subst_mod_type subst_mod;

  if (*prec==1){
    precision=OPT_LOW_PREC;
  }
  else if (*prec==2){
    precision=OPT_MED_PREC;
  }
  else{
    precision=OPT_HIGH_PREC;
  }

  msa=(MSA*)ad2ptr(*msaAddress);

  //from R, we either pass an initial model, or a string to create an 
  //initial model from
  if(*modAddress==0){
    tree = tr_new_from_string(*treestr);
  }
  else{
    init_mod=(TreeModel*)ad2ptr(*modAddress);
  }

  subst_mod = tm_get_subst_mod_type(*substModel);
  if (subst_mod==UNDEF_MOD){
    *error=1;
    strcpy(*errstr, "ERROR: Undefined Model Type\n");
    return;
  }

  if (init_mod == NULL){
    mod = tm_new(tr_create_copy(tree), NULL, NULL, subst_mod, 
		     msa->alphabet, 1, 1, NULL, -1);
    params = tm_params_init(mod, .1, 5, 1);  
  }  
  else {
    mod = tm_create_copy(init_mod);  
    tm_reinit(mod, subst_mod, 1, mod->alpha, NULL, NULL);
    params = tm_params_new_init_from_model(init_mod);
  }

  if (init_mod != NULL && mod->backgd_freqs != NULL) {
    vec_free(mod->backgd_freqs);
    mod->backgd_freqs = NULL; /* force re-estimation */
  }

  if (use_em)
    tm_fit_em(mod, msa, params, -1, precision, logf);
  else
    tm_fit(mod, msa, params, -1, precision, logf);

  *fittedAddress=ptr2ad(mod);

  *error=rphast_errno;
  *errstr=rphast_errmsg;

}
