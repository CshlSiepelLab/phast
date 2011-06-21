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
#include <memory_handler.h>

Vector *rph_get_vector(SEXP doubleP);
Matrix *rph_get_matrix(SEXP matP);
SEXP rph_listOfLists_to_SEXP(ListOfLists *lol);

SEXP rph_gff_new_extptr(GFF_Set *gff);
SEXP rph_msa_new_extptr(MSA *msa);
SEXP rph_ms_new_extptr(MS *ms);
SEXP rph_tm_new_extptr(TreeModel *tm);
SEXP rph_hmm_new_extptr(HMM *hmm);
