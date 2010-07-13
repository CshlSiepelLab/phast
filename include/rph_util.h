//#ifdef RPHAST

#include <Rdefines.h>
#undef Matrix
#include <matrix.h>
#include <vector.h>
#include <list_of_lists.h>

Vector *rph_get_vector(SEXP doubleP);
Matrix *rph_get_matrix(SEXP matP);
SEXP rph_listOfLists_to_SEXP(ListOfLists *lol);
//#endif
