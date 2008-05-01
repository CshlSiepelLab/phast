/*************************************************
rph_phyloP.c
Alexandra Denby
4/29/08

The RPHAST handles to phyloP, a part of the
phast package

NEEDS TO BE FIXED: pruneTree currently does something
funny when you try to call it twice.  Need to find
a workaround.

*************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <misc.h>
#include <msa.h>
#include <maf.h>
#include <tree_model.h>
#include <sufficient_stats.h>
#include <subst_distrib.h>
#include <prob_vector.h>
#include <prob_matrix.h>
#include "phyloP.h"
#include "fit_column.h"
#include <rph_util.h>


/***********************functions*********************************/

void preprocessMSA(double* msaAddress, int* error, char** errstr);
void pruneTree(MSA* msa, TreeModel* mod);
void toAlignmentSpace(double* align_array, double* tuple_array, MSA* msa);
mode_type getMode(int modeNum);
void phyloP_SPH(double* msaAddress, double* tmAddress, int* fm, double* eps, int* stats, int* modeNum, double* pVals, double* postMeans, double* postVars, int* error, char** errstr);
void phyloP_LRT(double* msaAddress, double* tmAddress, int* stats, int* modeNum, double* pvals, double* llrs, double* scales, int* error, char** errstr);
void phyloP_SCORE(double* msaAddress, double* tmAddress, int* stats, int* modeNum, double* pvals, double* teststats, double* derivs, int* error, char** errstr);

void preprocessMSA(double* msaAddress, int* error, char** errstr){

  MSA* msa=(MSA*)ad2ptr(*msaAddress);

  if (msa->ss == NULL)
    ss_from_msas(msa, 1, TRUE, NULL, NULL, NULL, -1);
  
  if (msa_alph_has_lowercase(msa)) msa_toupper(msa);     
  
  msa_remove_N_from_alph(msa);
  
  *error=rphast_errno;
  *errstr=rphast_errmsg;

}

void pruneTree(MSA* msa, TreeModel* mod){
  List* pruned_names;
  int old_nleaves;
  int j;

  pruned_names = lst_new_ptr(msa->nseqs);
  old_nleaves = (mod->tree->nnodes + 1) / 2;
  tm_prune(mod, msa, pruned_names);

  if (lst_size(pruned_names) >= old_nleaves){
    rphast_errno=1;
    strcpy(rphast_errmsg,"ERROR: no match for leaves of tree in alignment.\n");
  }
  else if (lst_size(pruned_names) > 0) {
    printf("WARNING: pruned away leaves with no match in alignment (");
    for (j = 0; j < lst_size(pruned_names); j++)
      printf("%s%s", ((String*)lst_get_ptr(pruned_names, j))->chars, 
	      j < lst_size(pruned_names) - 1 ? ", " : ").\n");
  }
  
}

/* Assigns values in tuple space from tuple array into preallocated
   align_array in alignment space, using msa as it's basis.
Sufficient stats should already have been calculated*/
void toAlignmentSpace(double* align_array, double* tuple_array, MSA* msa){
  int i;
  for(i=0; i<msa->length; i++){
    align_array[i]=tuple_array[msa->ss->tuple_idx[i]];
  }
}

mode_type getMode(int modeNum){
  if (modeNum==1){
    return(CON);
  }
  else if (modeNum==2){
    return(NNEUT);
  }
  else{
    return(ACC);
  }

}

void phyloP_SPH(double* msaAddress, double* tmAddress, int* fm, double* eps, int* stats, int* modeNum, double* pVals, double* postMeans, double* postVars, int* error, char** errstr){

  JumpProcess *jp,*jp_post;
  int nsites;
  double epsilon;
  MSA* msa;
  TreeModel* mod;
  int fit_model;
  double *tuple_pvals = NULL;
  double *tuple_post_means = NULL, *tuple_post_vars = NULL;
  double prior_mean, prior_var;
  FILE* logf = NULL;
  mode_type mode;

  msa=(MSA*)ad2ptr(*msaAddress);
  mod=(TreeModel*)ad2ptr(*tmAddress);
  nsites = msa->length;
  epsilon=*eps;
  fit_model=*fm;
  mode=getMode(*modeNum);

  //pruneTree(msa,mod);

  jp = sub_define_jump_process(mod, epsilon, tr_total_len(mod->tree));
  if (fit_model){
    jp_post = sub_define_jump_process(mod, epsilon, 
                                  10 * tr_max_branchlen(mod->tree));
  }
  else{
    jp_post = jp;
  }

  tuple_pvals = smalloc(msa->ss->ntuples * sizeof(double));

  if (*stats){
    tuple_post_means = smalloc(msa->ss->ntuples * sizeof(double));
    tuple_post_vars = smalloc(msa->ss->ntuples * sizeof(double));
  }

  sub_pval_per_site(jp, msa, mode, fit_model, &prior_mean, &prior_var, 
		    tuple_pvals, tuple_post_means, tuple_post_vars,
		    logf);

  toAlignmentSpace(pVals, tuple_pvals, msa);

  if(*stats){
    toAlignmentSpace(postMeans, tuple_post_means, msa);
    toAlignmentSpace(postVars, tuple_post_vars, msa);
  }

  *error=rphast_errno;
  *errstr=rphast_errmsg;
  

}

void phyloP_LRT(double* msaAddress, double* tmAddress, int* stats, int* modeNum, double* pvals, double* llrs, double* scales, int* error, char** errstr){
  MSA* msa=(MSA*)ad2ptr(*msaAddress);
  TreeModel* mod=(TreeModel*)ad2ptr(*tmAddress);
  double *tuple_pvals = NULL;
  double *tuple_llrs = NULL, *tuple_scales = NULL;
  mode_type mode;
  FILE* logf = NULL;

  mode=getMode(*modeNum);
  
  //pruneTree(msa,mod);
  
  tuple_pvals = smalloc(msa->ss->ntuples * sizeof(double));
  if (*stats){
        tuple_llrs = smalloc(msa->ss->ntuples * sizeof(double));
        tuple_scales = smalloc(msa->ss->ntuples * sizeof(double));
  }


  col_lrts(mod, msa, mode, tuple_pvals, tuple_scales, tuple_llrs, logf);

  toAlignmentSpace(pvals, tuple_pvals, msa);

  if(*stats){
    toAlignmentSpace(llrs, tuple_llrs, msa);
    toAlignmentSpace(scales, tuple_scales, msa);
  }

  *error=rphast_errno;
  *errstr=rphast_errmsg;


}

void phyloP_SCORE(double* msaAddress, double* tmAddress, int* stats, int* modeNum, double* pvals, double* teststats, double* derivs, int* error, char** errstr){

  MSA* msa=(MSA*)ad2ptr(*msaAddress);
  TreeModel* mod=(TreeModel*)ad2ptr(*tmAddress);
  double *tuple_pvals = NULL;
  double *tuple_teststats = NULL, *tuple_derivs = NULL;
  mode_type mode;

  mode=getMode(*modeNum);
  
  //pruneTree(msa,mod);

  tuple_pvals = smalloc(msa->ss->ntuples * sizeof(double));
  if (*stats){
    tuple_teststats = smalloc(msa->ss->ntuples * sizeof(double));
    tuple_derivs = smalloc(msa->ss->ntuples * sizeof(double));
  }
  
  col_score_tests(mod, msa, mode, tuple_pvals, tuple_derivs,                         tuple_teststats);

 toAlignmentSpace(pvals, tuple_pvals, msa);

  if(*stats){
    toAlignmentSpace(teststats, tuple_teststats, msa);
    toAlignmentSpace(derivs, tuple_derivs, msa);
  }

  *error=rphast_errno;
  *errstr=rphast_errmsg;
  
}

void phyloP_GERP(double* msaAddress, double* tmAddress, int* stats, int* modeNum, double* nrejected, double* nneut, double* nobs, double* nspec, int* error, char** errstr){

  MSA* msa=(MSA*)ad2ptr(*msaAddress);
  TreeModel* mod=(TreeModel*)ad2ptr(*tmAddress);
  double *tuple_nrejected = NULL;
  double *tuple_nneut = NULL, *tuple_nobs = NULL, *tuple_nspec=NULL;
  mode_type mode;
  FILE* logf = NULL;

  mode=getMode(*modeNum);
  
  //pruneTree(msa,mod);

  tuple_nrejected = smalloc(msa->ss->ntuples * sizeof(double));

  if (*stats){
    tuple_nneut = smalloc(msa->ss->ntuples * sizeof(double));
    tuple_nobs = smalloc(msa->ss->ntuples * sizeof(double));
    tuple_nspec = smalloc(msa->ss->ntuples * sizeof(double));
  }

  col_gerp(mod, msa, tuple_nneut, tuple_nobs, tuple_nrejected,              
	   tuple_nspec, logf);


  toAlignmentSpace(nrejected, tuple_nrejected, msa);

  if(*stats){
    toAlignmentSpace(nneut, tuple_nneut, msa);
    toAlignmentSpace(nobs, tuple_nobs, msa);
    toAlignmentSpace(nspec, tuple_nspec, msa);
  }

  *error=rphast_errno;
  *errstr=rphast_errmsg;

}
