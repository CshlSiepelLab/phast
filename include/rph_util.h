//NOTE: this file is only included by files in src/lib/rphast
#include <Rdefines.h>
#undef Matrix
#undef nrows
#include <misc.h>
#include <matrix.h>
#include <vector.h>
#include <list_of_lists.h>
#include <stringsplus.h>
#include <complex_matrix.h>
#include <gap_patterns.h>

SEXP rph_free_all();
Vector *rph_get_vector(SEXP doubleP);
Matrix *rph_get_matrix(SEXP matP);
SEXP rph_listOfLists_to_SEXP(ListOfLists *lol);

void rph_mem_protect(void *ptr);
void rph_lst_protect(List *l);
void rph_str_protect(String *s);
void rph_vec_protect(Vector *v);
void rph_zvec_protect(Zvector *v);
void rph_mat_protect(Matrix *m);
void rph_zmat_protect(Zmatrix *m);
void rph_msa_protect(MSA *msa);
void rph_gff_protect(GFF_Set *gff);
void rph_tree_protect(TreeNode *tr);
void rph_mm_protect(MarkovMatrix *mm);
void rph_tm_protect(TreeModel *tm);
void rph_hmm_protect(HMM *hmm);
void rph_cm_protect(CategoryMap *cm);
void rph_gp_protect(GapPatternMap *gpm);

void rph_register_protected_object(void *ptr, void (*protect_function)(void*));
void rph_unregister_protected(void *ptr);
void rph_gff_register_protect(GFF_Set *gff);
void rph_msa_register_protect(MSA *msa);
void rph_tm_register_protect(TreeModel *tm);
void rph_hmm_register_protect(HMM *hmm);

SEXP rph_gff_new_extptr(GFF_Set *gff);
SEXP rph_msa_new_extptr(MSA *msa);
SEXP rph_tm_new_extptr(TreeModel *tm);
SEXP rph_hmm_new_extptr(HMM *hmm);
