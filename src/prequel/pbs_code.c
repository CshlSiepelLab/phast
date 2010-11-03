/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** discrete encodings of probabilistic biological sequences */


#include <misc.h>
#include <pbs_code.h>
#include <time.h>

PbsCode *pbs_new(int dim, int nrows, int nbytes) {
  int i;
  PbsCode *retval = smalloc(sizeof(PbsCode));
  retval->max_size = ~(~0 << (8*nbytes)); /* e.g., 255 for nbytes = 1, 65535
					    for nbytes = 2; note that
					    max_size itself will be a
					    reserved code, for gaps */
  retval->sg = sxg_build_grid(dim, nrows);
  retval->rp = smalloc(retval->max_size * sizeof(void*));
  retval->nbytes = nbytes;
  retval->code_size = retval->sg->nregs;
  retval->gap_code = retval->max_size;

  if (retval->code_size >= retval->max_size)
    die("pbs_new: retval->code_size %i >= retval->max_size %i", 
	retval->code_size, retval->max_size);
  if (nbytes > MAX_NBYTES)
    die("pbs_new: nbytes (%i) <= %i", nbytes, MAX_NBYTES);

  /* initialize representative points to centroids of simplex regions */
  retval->codes_by_region = smalloc(retval->sg->nregs * sizeof(void*));
  for (i = 0; i < retval->sg->nregs; i++) {
    retval->rp[i] = vec_create_copy(retval->sg->sr[i]->centroid);
    retval->codes_by_region[i] = lst_new_int(1);
    lst_push_int(retval->codes_by_region[i], i);
  }

  return retval;
}

/* like pbs_new, but take given code_size and don't initialize based
   on simplex; for use when code is known */
PbsCode *pbs_new_shell(int dim, int nrows, int nbytes, int code_size) {
  PbsCode *retval = smalloc(sizeof(PbsCode));
  int i;
  retval->max_size = ~(~0 << (8*nbytes)); /* e.g., 255 for nbytes = 1, 65535
					    for nbytes = 2; note that
					    max_size itself will be a
					    reserved code, for gaps */
  if (code_size > retval->max_size)
    die("pbs_new_shell: code_size (%i) > retval->max_size (%i)\n",
	code_size, retval->max_size);
  if (nbytes > MAX_NBYTES)
    die("pbs_new_shell: nbytes (%i) > MAX_NBYTES\n", nbytes, MAX_NBYTES);

  retval->sg = sxg_build_grid(dim, nrows);
  retval->rp = smalloc(code_size * sizeof(void*));
  for (i = 0; i < code_size; i++) retval->rp[i] = NULL;
  retval->nbytes = nbytes;
  retval->code_size = code_size;
  retval->gap_code = retval->max_size;
  retval->codes_by_region = smalloc(retval->sg->nregs * sizeof(void*));
  for (i = 0; i < retval->sg->nregs; i++) 
    retval->codes_by_region[i] = lst_new_int(1);

  return retval;
}

void pbs_free(PbsCode *code) {
  int i;
  for (i = 0; i < code->sg->nregs; i++)
    lst_free(code->codes_by_region[i]);
  sfree(code->codes_by_region);
  sxg_free_grid(code->sg);
  for (i = 0; i < code->code_size; i++)
    vec_free(code->rp[i]);
  sfree(code->rp);
  sfree(code);
}

PbsCode *pbs_new_from_file(FILE *F) {
  Regex *nrows_re = str_re_new("##NROWS[[:space:]]*=[[:space:]]*([0-9]+)"),
    *dimension_re = str_re_new("##DIMENSION[[:space:]]*=[[:space:]]*([0-9]+)"),
    *nbytes_re = str_re_new("##NBYTES[[:space:]]*=[[:space:]]*([0-9]+)"),
    *codesize_re = str_re_new("##CODESIZE[[:space:]]*=[[:space:]]*([0-9]+)");
  String *line = str_new(STR_MED_LEN);
  List *fields = lst_new_ptr(50);
  int nrows = -1, dimension = -1, nbytes = -1, codesize = -1;
  PbsCode *code = NULL;

  while (str_readline(line, F) != EOF) {
    str_trim(line);
    if (line->length == 0) continue;

    if (code == NULL) {		/* still reading header */
      if (str_re_match(line, nrows_re, fields, 1) >= 0) 
	str_as_int(lst_get_ptr(fields, 1), &nrows);
      else if (str_re_match(line, dimension_re, fields, 1) >= 0) 
	str_as_int(lst_get_ptr(fields, 1), &dimension);
      else if (str_re_match(line, nbytes_re, fields, 1) >= 0) 
	str_as_int(lst_get_ptr(fields, 1), &nbytes);
      else if (str_re_match(line, codesize_re, fields, 1) >= 0) 
	str_as_int(lst_get_ptr(fields, 1), &codesize);
      else if (line->chars[0] == '#') continue; /* comment between
						   header lines */
      else die("ERROR: malformed header in code file.\n");

      if (nrows >= 0 && dimension >= 0 && nbytes >= 0 && codesize >= 0) 
	code = pbs_new_shell(dimension, nrows, nbytes, codesize);
    }

    else if (line->chars[0] == '#') continue;

    else {			/* code lines */
      Vector *p;
      int i, codeidx;
      double tmpprob;

      str_split(line, NULL, fields);
      if (lst_size(fields) != dimension + 1)
	die("ERROR: code lines in code file must have dimension + 1 columns.\n");

      if (str_as_int(lst_get_ptr(fields, 0), &codeidx) != 0 ||
	  codeidx < 0 || codeidx >= codesize)
	die ("ERROR: bad index in code file ('%s')\n", lst_get_ptr(fields, 0));

      p = vec_new(dimension);
      for (i = 0; i < dimension; i++) {
	if (str_as_dbl(lst_get_ptr(fields, i+1), &tmpprob) != 0 ||
	    tmpprob < 0 || tmpprob > 1)
	  die("ERROR: bad probability in code file ('%s')\n", lst_get_ptr(fields, i+1));
	vec_set(p, i, tmpprob);
      }

      if (code->rp[codeidx] != NULL) 
	die("ERROR: nonunique code index in code file (%d)\n", codeidx);

      code->rp[codeidx] = p;
    }

    lst_free_strings(fields);
  }

  pbs_assign_points(code);

  str_free(line);
  lst_free(fields);
  str_re_free(nrows_re);
  str_re_free(dimension_re);
  str_re_free(nbytes_re);
  str_re_free(codesize_re);

  return code;
}

void pbs_write(PbsCode *c, FILE *F, char *comment) {
  int i, j;

  fprintf(F, "##NROWS = %d\n", c->sg->nrows);
  fprintf(F, "##DIMENSION = %d\n", c->sg->d);
  fprintf(F, "##NBYTES = %d\n", c->nbytes);
  fprintf(F, "##CODESIZE = %d\n\n", c->code_size);

  if (comment != NULL) 		/* line formatting and '#' prefix must be
				   taken care of in calling code  */
    fprintf(F, "%s\n", comment);

  fprintf(F, "# Each index of the code is shown below with its representative probability\n\
# vector (p1, p2, ..., pd).\n\n");

  fprintf(F, "#code_index p1 p2 ...\n");

  for (i = 0; i < c->code_size; i++) {
    fprintf(F, "%d\t", i);
    for (j = 0; j < c->sg->d; j++) 
      fprintf(F, "%f%s", c->rp[i]->data[j], j == c->sg->d - 1? "\n" : "\t");
  }
}

/* assign representative points to simplex regions */
void pbs_assign_points(PbsCode *c) {
  int i;
  for (i = 0; i < c->sg->nregs; i++) lst_clear(c->codes_by_region[i]); 
  for (i = 0; i < c->code_size; i++) {
    SimplexRegion *r = sxg_get_region(c->sg, c->rp[i]);
    lst_push_int(c->codes_by_region[r->idx], i);
  }    
}

/* get code index for probability vector; if 'errorVal' is non-null, it
   will be set equal to the symmetric KL divergence between the vector
   and the representative point */
unsigned pbs_get_index(PbsCode *code, Vector *p, double *errorVal) {
  unsigned retval=-1;
  double min_d = INFTY + 1;	/* because min distance could be INFTY */
  int i;
  SimplexRegion *r = sxg_get_region(code->sg, p);
  int ncodes = lst_size(code->codes_by_region[r->idx]);

  if (ncodes == 0)
    die("ERROR: no representative points for simplex region.\n");
  else if (ncodes == 1 && errorVal == NULL)
    return lst_get_int(code->codes_by_region[r->idx], 0);

  for (i = 0; i < ncodes; i++) {
    int idx = lst_get_int(code->codes_by_region[r->idx], i);
    double d = sym_rel_entropy(p->data, code->rp[idx]->data, p->size);
    if (d < min_d) {
      retval = idx;
      min_d = d;
    }
  }

  if (errorVal != NULL) *errorVal = min_d;

  return retval;
}

/* save a copy of the representative points for a given simplex
   region (used below) */
void save_points(PbsCode *code, int region_idx, Vector **copy) {
  int i;
  for (i = 0; i < lst_size(code->codes_by_region[region_idx]); i++) 
    vec_copy(copy[i], code->rp[lst_get_int(code->codes_by_region[region_idx], i)]);
}

/* restore saved representative points (used below) */
void restore_points(PbsCode *code, int region_idx, Vector **copy) {
  int i;
  for (i = 0; i < lst_size(code->codes_by_region[region_idx]); i++) 
    vec_copy(code->rp[lst_get_int(code->codes_by_region[region_idx], i)], copy[i]);
}


/* for a given region, assign vectors to representative points and
   return total error. */
double assign_vectors(PbsCodeTrainingData *td, int region_idx) {
  int i, code;
  double error, tot_error = 0;  

  /* clear previous assignment */
  for (i = 0; i < lst_size(td->code->codes_by_region[region_idx]); i++) {
    code = lst_get_int(td->code->codes_by_region[region_idx], i);
    lst_clear(td->vectors_by_code[code]);
    lst_clear(td->counts_by_code[code]);
    td->error_by_code[code] = 0;
  }

  for (i = 0; i < lst_size(td->vectors_by_region[region_idx]); i++) {
    code = pbs_get_index(td->code, lst_get_ptr(td->vectors_by_region[region_idx], i), 
			 &error);
    error *= lst_get_int(td->counts_by_region[region_idx], i);
    tot_error += error;
    td->error_by_code[code] += error;
    lst_push_ptr(td->vectors_by_code[code], 
		 lst_get_ptr(td->vectors_by_region[region_idx], i));
    lst_push_int(td->counts_by_code[code], 
		 lst_get_int(td->counts_by_region[region_idx], i));
  }

  td->error_by_region[region_idx] = tot_error;
  return tot_error;
}

/* refine the representative points within a given simplex region
   using a variant of the k-means algorithm.  Initialization must be
   external */
double k_means(PbsCodeTrainingData *td, int region_idx) {
  int i, code, niterations = 0;
  double err, olderr = INFTY;

  while (TRUE) {
    err = assign_vectors(td, region_idx);

    if (err >= olderr) {
      err = olderr;		/* rp will still reflect last
				   iteration; vectors_by_code
				   etc. will be out of sync but that's
				   okay */
      break;
    }

    for (i = 0; i < lst_size(td->code->codes_by_region[region_idx]); i++) {
      code = lst_get_int(td->code->codes_by_region[region_idx], i);
      if (lst_size(td->vectors_by_code[code]) > 0)
	vec_ave(td->code->rp[code], td->vectors_by_code[code], 
		td->counts_by_code[code]);
    }
    
    niterations++;

    if (niterations == 50) break;

    olderr = err;
  }

  return err;
}

/* attempt to optimize codes for given region; uses k-means algorithm
   with multiple restarts */
double pbs_optimize_region(PbsCodeTrainingData *td, int region_idx,
			   FILE *log_f) {
  int i, trial, idx, 
    ncodes = lst_size(td->code->codes_by_region[region_idx]);
  unsigned long nchoices;	/* can get large */
  Vector *freqs = vec_new(lst_size(td->vectors_by_region[region_idx]));
  double error, best_error;
  Vector **best_rp = smalloc(ncodes * sizeof(void*));

  /* find initial assignment and error, save initial representative
     points */
  best_error = assign_vectors(td, region_idx);
  if (log_f != NULL)
    fprintf(log_f, "Before k-means: %f\n", best_error);
  for (i = 0; i < ncodes; i++) best_rp[i] = vec_new(td->code->sg->d);
  save_points(td->code, region_idx, best_rp);

  /* now try several initializations of k-means and take the best one.
     Use elements in region_vectors as starting points.  If the number
     of vectors is small enough, try all combinations; otherwise, make
     several random draws, weighting the vectors by their counts */

  nchoices = combinations(lst_size(td->vectors_by_region[region_idx]), ncodes);

  if (nchoices < 50) {
    int *index = smalloc(ncodes * sizeof(int));
    index[0] = -1;		/* used by next_comb */
    while (next_comb(lst_size(td->vectors_by_region[region_idx]), ncodes, index)) {
      for (i = 0; i < ncodes; i++) 
	vec_copy(td->code->rp[lst_get_int(td->code->codes_by_region[region_idx], i)], 
		 lst_get_ptr(td->vectors_by_region[region_idx], index[i]));
      error = k_means(td, region_idx);
      if (error < best_error) {
	save_points(td->code, region_idx, best_rp);
	best_error = error;
      }
    }
    sfree(index);
  }
  else {
    for (trial = 0; trial < 10; trial++) {
      /* randomly draw ncodes starting points from region_vectors,
	 without replacement */
      for (i = 0; i < lst_size(td->vectors_by_region[region_idx]); i++) 
	vec_set(freqs, i, lst_get_int(td->counts_by_region[region_idx], i));
      for (i = 0; i < ncodes; i++) {
	normalize_probs(freqs->data, freqs->size);
	idx = draw_index(freqs->data, freqs->size);
	vec_copy(td->code->rp[lst_get_int(td->code->codes_by_region[region_idx], i)], 
		 lst_get_ptr(td->vectors_by_region[region_idx], idx));
	freqs->data[idx] = 0;	/* ensures won't be drawn again */
      }
      error = k_means(td, region_idx);
      if (error < best_error) {
	save_points(td->code, region_idx, best_rp);
	best_error = error;
      }
    }
  }

  restore_points(td->code, region_idx, best_rp);
  best_error = assign_vectors(td, region_idx);

  if (log_f != NULL)
    fprintf(log_f, "After k-means: %f\n", best_error);

  for (i = 0; i < ncodes; i++) vec_free(best_rp[i]);
  sfree(best_rp);
  vec_free(freqs);

  return best_error;
}

PbsCodeTrainingData *pbs_new_training_data(PbsCode *code, List *prob_vectors, 
					   List *counts) {
  int i, init_size;
  PbsCodeTrainingData *td;

  if (lst_size(prob_vectors) != lst_size(counts))
    die("ERROR: pbs_new_training_data: prob_vectors of different size (%i) than counts (%i)\n", lst_size(prob_vectors), lst_size(counts));
  if (lst_size(prob_vectors) <= 0)
    die("ERROR: pbs_new_training_data: prob_vectors must have size > 0 (has size %i)\n", lst_size(prob_vectors));

  td = smalloc(sizeof(PbsCodeTrainingData));
  td->code = code;
  td->prob_vectors = prob_vectors;
  td->counts = counts;
  td->vectors_by_region = smalloc(code->sg->nregs * sizeof(void*));
  td->counts_by_region = smalloc(code->sg->nregs * sizeof(void*));
  td->vectors_by_code = smalloc(code->max_size * sizeof(void*));
  td->counts_by_code = smalloc(code->max_size * sizeof(void*));
  td->error_by_region = smalloc(code->sg->nregs * sizeof(double));
  td->error_by_code = smalloc(code->max_size * sizeof(double));
  init_size = max(5, lst_size(prob_vectors) / code->sg->nregs);

  for (i = 0; i < code->sg->nregs; i++) {
    td->vectors_by_region[i] = lst_new_ptr(init_size);
    td->counts_by_region[i] = lst_new_int(init_size);
    td->error_by_region[i] = 0;
  }
  for (i = 0; i < code->max_size; i++) {
    td->vectors_by_code[i] = lst_new_ptr(init_size);
    td->counts_by_code[i] = lst_new_int(init_size);
    td->error_by_code[i] = 0;
  }

  /* vectors/counts by region can be defined up front and won't
     change */
  for (i = 0; i < lst_size(prob_vectors); i++) {
    Vector *v = lst_get_ptr(prob_vectors, i);
    SimplexRegion *reg = sxg_get_region(code->sg, v);
    lst_push_ptr(td->vectors_by_region[reg->idx], v);
    lst_push_int(td->counts_by_region[reg->idx], lst_get_int(counts, i));
  }

  return td;
}

void pbs_free_training_data(PbsCodeTrainingData *td) {
  int i;
  for (i = 0; i < td->code->sg->nregs; i++) {
    lst_free(td->vectors_by_region[i]);
    lst_free(td->counts_by_region[i]);
  }
  for (i = 0; i < td->code->max_size; i++) {
    lst_free(td->vectors_by_code[i]);
    lst_free(td->counts_by_code[i]);
  }
  sfree(td->vectors_by_region);
  sfree(td->counts_by_region);
  sfree(td->vectors_by_code);
  sfree(td->counts_by_code);
  sfree(td->error_by_region);
  sfree(td->error_by_code);
  sfree(td);
}

/* returns average training error */
/* works with any initial set of representative points */
double pbs_estimate_from_data(PbsCode *code, List *prob_vectors, 
			      List *counts, FILE *logf, 
			      training_mode mode) {
  int i, j, tot_count = 0;
  unsigned idx;
  double tot_error = 0;
  PbsCodeTrainingData *td = pbs_new_training_data(code, prob_vectors, counts);

  if (lst_size(prob_vectors) != lst_size(counts))
    die("ERROR: pbs_estimate_from_data: prob_vectors of different size (%i) than counts (%i)\n", lst_size(prob_vectors), lst_size(counts));
  if (lst_size(prob_vectors) <= 0)
    die("ERROR: pbs_estimate_from_data: prob_vectors must have size > 0 (has size %i)\n", lst_size(prob_vectors));

  for (i = 0; i < lst_size(counts); i++) 
    tot_count += lst_get_int(counts, i);

  /* initialize by setting representative point for each code index to
     pointwise average of assigned vectors; this is a first order
     optimization */
  for (i = 0; i < code->sg->nregs; i++) { /* process by region */
    assign_vectors(td, i);
    for (j = 0; j < lst_size(code->codes_by_region[i]); j++) {
      idx = lst_get_int(code->codes_by_region[i], j);
      if (lst_size(td->vectors_by_code[idx]) > 0)
	vec_ave(td->code->rp[idx], td->vectors_by_code[idx], 
		td->counts_by_code[idx]);
    }
    assign_vectors(td, i);	/* needed to update error_by_code */
  }

  /* output to log */
  if (logf != NULL) {
    for (i = 0; i < code->code_size; i++) {
      fprintf(logf, "%5d ", i);
      for (j = 0; j < code->rp[i]->size; j++)
	fprintf(logf, "%7.3f ", code->rp[i]->data[j]);
      fprintf(logf, "%6d %9.1f\n", lst_size(td->vectors_by_code[i]), 
	      fabs(td->error_by_code[i]));
    }
  }      

  /* now add code indices greedily until all have been used */
  if (mode == FULL) {
    for (i = code->code_size; i < code->max_size; i++) {
      double max_error = -1;
      int worst = -1;

      /* identify region with worst error */
      for (j = 0; j < code->sg->nregs; j++) {
	if (lst_size(code->codes_by_region[j]) > 8) 
	  continue;		/* don't add more than 8 codes to a region */ 
	if (td->error_by_region[j] > max_error) {
	  max_error = td->error_by_region[j];
	  worst = j;
	}
      }

      if (max_error == 0) break; /* possible if code size is as large
				    as number of input vectors */

      if (logf != NULL) {
	fprintf(logf, "\nWorst region: %d\n", worst);
	if (lst_size(code->codes_by_region[worst]) == 8)
	  fprintf(logf, "(Now max no. representative points)\n");
      }

      code->rp[i] = vec_create_copy(code->sg->sr[worst]->centroid); 
      /* (arbitrary initialization) */
      lst_push_int(code->codes_by_region[worst], i);
      code->code_size++;

      pbs_optimize_region(td, worst, logf);
    }
  }

   if (logf != NULL) { 
     int k; 
     fprintf(logf, "\n\n"); 
     for (i = 0; i < code->code_size; i++) { 
       fprintf(logf, "%3d ", i); 
       for (j = 0; j < code->rp[i]->size; j++) 
	 fprintf(logf, "%7.3f ", code->rp[i]->data[j]); 
       fprintf(logf, "\n"); 
       for (j = 0; j < lst_size(td->vectors_by_code[i]); j++) { 
	 Vector *v = lst_get_ptr(td->vectors_by_code[i], j); 
	 fprintf(logf, "    "); 
	 for (k = 0; k < v->size; k++) 
	   fprintf(logf, "%7.3f ", v->data[k]); 
	 fprintf(logf, "%8d\n", lst_get_int(td->counts_by_code[i], j)); 
       } 
     } 
   } 

  tot_error = 0; 
  for (i = 0; i < code->code_size; i++)
    tot_error += td->error_by_code[i];

  pbs_free_training_data(td);
  
  return tot_error / tot_count;
}

/* write code index in binary form, allowing for variable nbytes */
void pbs_write_binary(PbsCode *code, unsigned code_idx, FILE *F) {
  unsigned char bytes[MAX_NBYTES];
  int i;

  /* as discussed in Kernighan & Pike (The Practice of Programming;
     pp. 206-207), write in canonical format, to avoid
     incompatibilities between little endian and big endian
     architectures.  We'll go from high- to low-order */

  for (i = code->nbytes - 1; i >= 0; i--) {
    bytes[i] = code_idx & ~(~0 << 8); /* see Kernighan and Ritchie, p. 49 */
    if (i > 0) code_idx >>= 8;
  }

  fwrite(bytes, code->nbytes, 1, F);
}

/* read code index in binary form, allowing for variable nbytes;
   returns EOF when end of file is reached */
int pbs_read_binary(PbsCode *code, unsigned *code_idx, FILE *F) {
  unsigned char bytes[MAX_NBYTES];
  int i;
  *code_idx = 0;
  if (fread(bytes, code->nbytes, 1, F) == 0) 
    return EOF;
  for (i = 0; i < code->nbytes; i++) {
    if (i > 0) *code_idx <<= 8;
    *code_idx |= bytes[i];
  }
  return 0;
}

