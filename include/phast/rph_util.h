//NOTE: this file is only included by files in src/lib/rphast
#include <Rdefines.h>
#undef Matrix
#undef nrows
#include <phast/misc.h>
#include <phast/matrix.h>
#include <phast/vector.h>
#include <phast/list_of_lists.h>
#include <phast/stringsplus.h>
#include <phast/complex_matrix.h>
#include <phast/gap_patterns.h>
#include <phast/memory_handler.h>

Vector *rph_get_vector(SEXP doubleP);
Matrix *rph_get_matrix(SEXP matP);
SEXP rph_listOfLists_to_SEXP(ListOfLists *lol);

SEXP rph_gff_new_extptr(GFF_Set *gff);
SEXP rph_msa_new_extptr(MSA *msa);
SEXP rph_ms_new_extptr(MS *ms);
SEXP rph_tm_new_extptr(TreeModel *tm);
SEXP rph_hmm_new_extptr(HMM *hmm);
