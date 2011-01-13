/** discrete encodings of probabilistic biological sequences */

#ifndef PBS_CODE
#define PBS_CODE

#include <simplex_grid.h>
#include <lists.h>
#include <stringsplus.h>

#define MAX_NBYTES 4
/* needs to be <= sizeof(int) */

typedef struct {
  SimplexGrid *sg;		/** Simplex grid */
  Vector **rp;		        /** array of representative points, by code
				    index (array has code_size elements) */
  List **codes_by_region;	/** code indices by simplex region;
				    codes_by_region[i] is a list of
				    the code indices associated with
				    simplex region i (array has
				    sg->nregs elements) */
  unsigned nbytes;		/** bytes of storage per encoded
				    "letter"; between 1 and MAX_NBYTES */
  unsigned max_size;		/** maximum code size: 256^nbytes - 1  (because of
				    reserved code for gaps) */
  unsigned code_size;	        /** actual code size; at most
				    max_size */
  unsigned gap_code;		/** reserved code for gap characters */
} PbsCode;

/* package of data used in code estimation */
typedef struct {
  PbsCode *code;		/** code being estimated */
  List *prob_vectors;		/** raw training data: probability vectors */
  List *counts;			/** raw training data: counts  */
  List **vectors_by_region;	/** element i is list of vectors for
				    simplex region i */
  List **counts_by_region;	/** element i is list of counts for
				    simplex region i */
  List **vectors_by_code;	/** element i is list of vectors for
				    code index i */
  List **counts_by_code;  	/** element i is list of counts for
				    code index i */
  double *error_by_region;	/** element i is contribution of
				    region i to total approximation
				    error */
  double *error_by_code;	/** element i is contribution of code
				    index i to total approximation
				    error */
} PbsCodeTrainingData;

typedef enum {FULL, NO_GREEDY, NO_TRAIN} training_mode;

PbsCode *pbs_new(int dim, int nrows, int nbytes);
PbsCode *pbs_new_shell(int dim, int nrows, int nbytes, int code_size);
void pbs_free(PbsCode *code);
PbsCode *pbs_new_from_file(FILE *F);
void pbs_write(PbsCode *c, FILE *F, char *comment);
void pbs_assign_points(PbsCode *c);
unsigned pbs_get_index(PbsCode *code, Vector *p, double *error);
double pbs_estimate_from_data(PbsCode *code, List *prob_vectors, List *counts,
			      FILE *logf, training_mode mode);
void pbs_write_binary(PbsCode *code, unsigned code_idx, FILE *F);
int pbs_read_binary(PbsCode *code, unsigned *code_idx, FILE *F);

#endif
