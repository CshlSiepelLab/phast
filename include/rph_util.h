//#ifdef RPHAST

#include <Rdefines.h>
#undef Matrix
#include <matrix.h>
#include <vector.h>
#include <list_of_lists.h>
#include <stringsplus.h>

void rph_protect_mem(void *ptr);
SEXP rph_free_all();
void rph_lst_protect(List *l);
void rph_string_protect(String *s);
Vector *rph_get_vector(SEXP doubleP);
Matrix *rph_get_matrix(SEXP matP);
SEXP rph_listOfLists_to_SEXP(ListOfLists *lol);
void rph_msa_protect(MSA *msa);
void rph_gff_protect(GFF_Set *gff);
//#endif
