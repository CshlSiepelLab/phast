/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#include "misc.h"
#include "sufficient_stats.h"
#include "maf.h"
#include "queues.h"

#define MAX_NTUPLE_ALLOC 100000
                                /* maximum number of tuples to
                                   accommodate initially */

/* Given a multiple alignment object, create a representation based on
   its sufficient statistics -- i.e., the distinct columns that it
   includes, the number of times each one appears, and (if store_order
   == 1), the order in which they appear.  An 'MSA_SS' object will be
   created and linked to the provided alignment object.  If the
   argument 'source_msa' is non-NULL, then the stats for the specified
   alignment will be *added* to the MSA_SS object of 'msa'.  This
   option can be used in repeated calls to create aggregate
   representations of sets of alignments (see ss_pooled_from_msas and
   ss_aggregate_from_files below).  In this case, a running hash table
   should also be passed in (existing_hash).  The argument
   'tuple_size' determines what size of column tuple to consider
   (e.g., if tuple_size==1, then each column is considered
   individually, and if tuple_size==2, then each is considered wrt its
   predecessor).  If store_order==1, then the order of the columns
   will be stored.

   If source_msa != NULL, then any sequences in the 'msa' object will
   be ignored, but 'msa' must be initialized with the appropriate
   sequence names, alphabet, and number of sequences (length will be
   adjusted as needed).  A specified source alignment must have the
   same number of sequences as 'msa' and the sequences must appear in
   corresponding order.  

   If msa->ncats > 0 and msa->categories != NULL then
   category-by-category counts will be maintained.  If in addition
   source_msas != NULL, then each source alignment will be expected to
   have non-NULL 'categories' vectors containing values of at most
   msa->ncats.  The optional list 'cats_to_do' (expected to be a list
   of integers if non-NULL) allows consideration to be limited to
   certain categories.  All of this does not apply when source
   alignments are represented by (unordered) sufficient statistics
   only, in which case category-by-category counts are maintained iff
   they exist for source alignments.  

   If non_overlapping == TRUE and tuple_size > 1, then will only
   convert non_overlapping tuples into the SS structure.  The first
   tuple will start with position 1 and be of length tuple_size, the
   second will start at tuple_size+1 and go to tuple_size*2, etc.
   (This is useful for converting coding sequence into codons).

*/

/* ADDENDUM: added an 'idx_offset' parameter for use when storing
   order and source alignments refer to local segments of a long
   alignment.  Will be added to coordinates in source_msa when setting
   tuple_idx for msa.  Set to -1 when not in use.  If idx_offset is
   non-negative, then it must be true that store_order == 1 and
   source_msa != NULL.  If idx_offset is non-negative, then msa->ss is
   assumed to be pre-allocated (offsets complicate reallocs). */

void ss_from_msas(MSA *msa, int tuple_size, int store_order, 
                  List *cats_to_do, MSA *source_msa, 
                  Hashtable *existing_hash,
                  int idx_offset, int non_overlapping) {
  int i, j, do_cats, idx, upper_bound;
  int max_tuples;
  MSA_SS *main_ss, *source_ss = NULL;
  Hashtable *tuple_hash = NULL;
  int *do_cat_number = NULL;
  char key[msa->nseqs * tuple_size + 1];
  MSA *smsa;
  int effective_offset = (idx_offset < 0 ? 0 : idx_offset); 

  if (source_msa == NULL && 
      (msa->seqs == NULL || msa->length <= 0 || msa->ss != NULL))
    die("ERROR: with no separate source alignment, ss_from_msas expects sequences of positive length and no SS object.\n");
  if (idx_offset >= 0) 
    if (!(store_order && source_msa != NULL))
      die("ERROR ss_from_msas: idx_offset=%i store_order=%i, source_msa=NULL=%i\n",
	  idx_offset, store_order, source_msa==NULL);
                                /* this is a little clumsy but it
                                   allows idx_offset both to signal
                                   the mode of usage and to specify
                                   the amount of the offset */
  if (store_order && source_msa != NULL && source_msa->seqs == NULL &&
      source_msa->ss->tuple_idx == NULL)
    die("ERROR: (ss_from_msas) must have ordered representation in source alignment for ordered destination alignment.\n");

  if (source_msa != NULL) {
    if (msa->nseqs != source_msa->nseqs)
      die("ERROR: (ss_from_msas) numbers of sequences must be equal in source and destination alignments.\n");
    if (msa->ncats >= 0 && source_msa->ncats >= 0 && 
        msa->ncats != source_msa->ncats)
      die("ERROR: (ss_from_msas) numbers of categories must be equal in source and destination alignments.\n");
  }

  do_cats = (msa->ncats >= 0);
  key[msa->nseqs * tuple_size] = '\0';
  if (do_cats && cats_to_do != NULL) {
    do_cat_number = (int*)smalloc((msa->ncats + 1) * sizeof(int));
    for (i = 0; i <= msa->ncats; i++) do_cat_number[i] = 0;
    for (i = 0; i < lst_size(cats_to_do); i++)
      do_cat_number[lst_get_int(cats_to_do, i)] = 1;
  }

  if (msa->ss == NULL) {
    if (source_msa != NULL && source_msa->ss != NULL)
      upper_bound = source_msa->ss->ntuples;
    else if (source_msa != NULL) upper_bound = source_msa->length;
    else upper_bound = min(msa->length, MAX_NTUPLE_ALLOC);
    max_tuples = pow_bounded(strlen(msa->alphabet)+ (int)strlen(msa->missing) + 1,
                         msa->nseqs * tuple_size,   
                     upper_bound);
    if (max_tuples < 0) max_tuples = upper_bound;  //protect against overflow
    ss_new(msa, tuple_size, max_tuples, do_cats, store_order);
    if (source_msa != NULL) msa->length = source_msa->length;
  }
  else {  //if (idx_offset < 0) {    
    // msa->ss != NULL so  we know source_msa != NULL
    /* if storing order based on source
                                   alignments and offset, then assume
                                   proper preallocation; otherwise
                                   (this case), realloc to accommodate
                                   new source msa */
    //    int newlen = effective_offset + msa->length + source_msa->length;
    int newlen = effective_offset + source_msa->length;
    msa_realloc(msa, newlen, newlen+100000, 
		do_cats, store_order);
    if (source_msa->ss != NULL) 
      upper_bound = msa->ss->ntuples + source_msa->ss->ntuples;
    else
      upper_bound = msa->ss->ntuples + source_msa->length;
      max_tuples = pow_bounded(strlen(msa->alphabet) + (int)strlen(msa->missing) + 1,
                         msa->nseqs * tuple_size,
                     upper_bound);
    if (max_tuples < 0) max_tuples = upper_bound;
    ss_realloc(msa, tuple_size, max_tuples, do_cats, store_order);
  }


  main_ss = msa->ss;
  tuple_hash = existing_hash != NULL ? existing_hash : hsh_new((int)((double)max_tuples/3));

  if (source_msa != NULL && source_msa->ss != NULL)
    source_ss = source_msa->ss;
  /* otherwise, source_ss will remain NULL */

  if (source_ss != NULL && !store_order && !non_overlapping) { /* in this case, just use
                                              existing suff stats */
    for (i = 0; i < source_ss->ntuples; i++) {
      checkInterruptN(i, 1000);
/*       fprintf(stderr, "col_tuple %d: %s\n", i, source_ss->col_tuples[i]); */

      if ((idx = ss_lookup_coltuple(source_ss->col_tuples[i], tuple_hash, msa)) == -1) {
 	idx = main_ss->ntuples++;
	ss_add_coltuple(source_ss->col_tuples[i], int_to_ptr(idx), tuple_hash, msa);
        main_ss->col_tuples[idx] = (char*)smalloc((tuple_size * msa->nseqs +1) 
						  * sizeof(char));
	main_ss->col_tuples[idx][msa->nseqs * tuple_size] = '\0';
        strncpy(main_ss->col_tuples[idx], source_ss->col_tuples[i],
		(msa->nseqs * tuple_size + 1));
                                /* NOTE: here main_ss must have
                                   sufficient size; we don't need to
                                   worry about reallocating */
      }

      main_ss->counts[idx] += source_ss->counts[i];
      if (do_cats && source_ss->cat_counts != NULL) {
        for (j = 0; j <= source_msa->ncats; j++) 
          main_ss->cat_counts[j][idx] += source_ss->cat_counts[j][i];
      }
    }
  }
  else {                        /* suff stats not available or storing
                                   order; in either case, go col-by-col */
    smsa = source_msa != NULL ? source_msa : msa; 
                                /* smsa is the source of the sequence
                                   data, whether msa itself or a
                                   separate source msa */

    for (i = 0; i < smsa->length; i++) { 
      checkInterruptN(i, 1000);
      if (non_overlapping &&  ((i+1) % tuple_size != 0)) continue;
      if (do_cats && cats_to_do != NULL && 
          do_cat_number[smsa->categories[i]] == 0) {
        if (store_order) main_ss->tuple_idx[i + effective_offset] = -1;
        continue;
      }

      if (smsa->seqs != NULL)
        col_to_string(key, smsa, i, tuple_size);
      else                      /* NOTE: must have ordered suff stats */
        strncpy(key, smsa->ss->col_tuples[smsa->ss->tuple_idx[i]], 
		(msa->nseqs * tuple_size + 1));

      if ((idx = ss_lookup_coltuple(key, tuple_hash, msa)) == -1) {
                                /* column tuple has not been seen
                                   before */
        idx = main_ss->ntuples++;
	ss_add_coltuple(key, int_to_ptr(idx), tuple_hash, msa);

        if (main_ss->ntuples > main_ss->alloc_ntuples) 
                                /* possible if allocated only for
                                   MAX_NTUPLE_ALLOC */
          ss_realloc(msa, tuple_size, main_ss->ntuples, do_cats, store_order);

        main_ss->col_tuples[idx] = (char*)smalloc((tuple_size * msa->nseqs + 1)
						  * sizeof(char));
	main_ss->col_tuples[idx][tuple_size * msa->nseqs] = '\0';
        strncpy(main_ss->col_tuples[idx], key, msa->nseqs * tuple_size);
      }

      main_ss->counts[idx]++;
      if (do_cats && smsa->categories != NULL) {
        if (!(smsa->categories[i] >= 0 && smsa->categories[i] <= msa->ncats))
	  die("ERROR ss_from_msas: smsa->categories[i]=%i should be in [0,%i]\n",
	      smsa->categories[i], msa->ncats);
        main_ss->cat_counts[smsa->categories[i]][idx]++;
      }
      if (store_order)
        main_ss->tuple_idx[i + effective_offset] = idx;
    }
  }

  if (existing_hash == NULL) {
    ss_compact(main_ss);        /* only compact if it looks like this
                                   function is not being called
                                   repeatedly */
    hsh_free(tuple_hash);
  }

  if (do_cats) sfree(do_cat_number);
}

/* creates a new sufficient statistics object and links it to the
   specified alignment.  Allocates sufficient space for max_ntuples
   distinct column tuples (see ss_compact, below).  The optional
   argument 'col_tuples' specifies an already existing array, so that
   a representation of the individual tuples can be shared (use NULL
   to create a new one; note that individual char* entries are
   initialized to NULL and must be allocated as needed).  */
void ss_new(MSA *msa, int tuple_size, int max_ntuples, int do_cats, 
            int store_order/* , char **col_tuples */) {
  int i, j;
  MSA_SS *ss;

  ss = (MSA_SS*)smalloc(sizeof(MSA_SS));
  msa->ss = ss;
  ss->msa = msa;
  ss->tuple_size = tuple_size;
  ss->ntuples = 0;
  ss->tuple_idx = NULL;
  ss->cat_counts = NULL;
  ss->alloc_len = max(1000, msa->length);
  if (store_order) {
    ss->tuple_idx = (int*)smalloc(ss->alloc_len * sizeof(int));
    for (i=0; i < ss->alloc_len; i++)
      ss->tuple_idx[i] = -1;
  }
  ss->col_tuples = (char**)smalloc(max_ntuples * sizeof(char*));
  for (i = 0; i < max_ntuples; i++) ss->col_tuples[i] = NULL;
  ss->counts = (double*)smalloc(max_ntuples * sizeof(double));
  for (i = 0; i < max_ntuples; i++) ss->counts[i] = 0; 
  if (do_cats) {
    if (msa->ncats < 0)
      die("ERROR ss_new: msa->ncats=%i (should be >=0)\n", msa->ncats);
    ss->cat_counts = (double**)smalloc((msa->ncats+1) * sizeof(double*));
    for (j = 0; j <= msa->ncats; j++) {
      ss->cat_counts[j] = (double*)smalloc(max_ntuples * sizeof(double));
      for (i = 0; i < max_ntuples; i++) 
        ss->cat_counts[j][i] = 0;
    }
  }
  ss->alloc_ntuples = max_ntuples;
}

/* ensures a suff stats object has enough memory allocated to
   accommodate new msa->length and max_ntuples */
void ss_realloc(MSA *msa, int tuple_size, int max_ntuples, int do_cats, 
                int store_order) {

  int i, j, cat_counts_done = FALSE, old_alloc_len;
  MSA_SS *ss = msa->ss;
  if (store_order && msa->length > ss->alloc_len) {
    old_alloc_len = ss->alloc_len;
    ss->alloc_len = max(ss->alloc_len * 2, msa->length);
    ss->tuple_idx = (int*)srealloc(ss->tuple_idx, ss->alloc_len * sizeof(int));
    for (i = old_alloc_len; i < ss->alloc_len; i++)
      ss->tuple_idx[i] = -1;
    
  }
  if (do_cats && ss->cat_counts == NULL) {
    if (msa->ncats < 0)
      die("ERROR ss_realloc: msa->ncats=%i (should be >=0)\n", msa->ncats);
    ss->cat_counts = smalloc((msa->ncats+1) * sizeof(double*));
    for (j = 0; j <= msa->ncats; j++) {
      ss->cat_counts[j] = smalloc(max_ntuples * sizeof(double));
      for (i = 0; i < max_ntuples; i++) ss->cat_counts[j][i] = 0;
    }
    cat_counts_done = TRUE;
  }
  if (max_ntuples > ss->alloc_ntuples) {
    int new_alloc_ntuples = max(max_ntuples, ss->alloc_ntuples*2);
/*     fprintf(stderr, "Realloc: max_ntuples = %d, alloc_ntuples = %d; realloc to %d\n", max_ntuples, ss->alloc_ntuples, new_alloc_ntuples); */
    ss->col_tuples = (char**)srealloc(ss->col_tuples, new_alloc_ntuples * 
                                     sizeof(char*));
    for (i = ss->alloc_ntuples; i < new_alloc_ntuples; i++) 
      ss->col_tuples[i] = NULL;

    ss->counts = (double*)srealloc(ss->counts, new_alloc_ntuples * 
                                  sizeof(double));
    for (i = ss->alloc_ntuples; i < new_alloc_ntuples; i++) 
      ss->counts[i] = 0; 
    if (do_cats && !cat_counts_done) {
      for (j = 0; j <= msa->ncats; j++) {
        ss->cat_counts[j] = 
          (double*)srealloc(ss->cat_counts[j], new_alloc_ntuples * 
                           sizeof(double));
        for (i = ss->alloc_ntuples; i < new_alloc_ntuples; i++) 
          ss->cat_counts[j][i] = 0;
      }
    }
    ss->alloc_ntuples = new_alloc_ntuples;
  }
}

/* Create a PooledMSA from a list of source msas.  Note that the same
   source_msas list will be used in the new object as was passed in.
   Also, this function assumes all msas have same names, nseqs, and
   alphabet (it uses those of the first) */
PooledMSA *ss_pooled_from_msas(List *source_msas, int tuple_size, int ncats, 
                               List *cats_to_do, int non_overlapping) {
  int i, j;
  MSA *rep_msa;
  PooledMSA *pmsa = (PooledMSA*)smalloc(sizeof(PooledMSA));
  Hashtable *tuple_hash = hsh_new(100000); /* wild guess on size.  Big
                                              enough?  Preview files? */
  char *key;

  if (lst_size(source_msas) <= 0)
    die("ERROR ss_pooled_from_msas: lst_size(source_msas)=%i, should be >0\n",
	lst_size(source_msas));
  rep_msa = (MSA*)lst_get_ptr(source_msas, 0);

  pmsa->pooled_msa = msa_new(NULL, NULL, rep_msa->nseqs, 0, rep_msa->alphabet);
  pmsa->source_msas = source_msas;
  pmsa->pooled_msa->names = (char**)smalloc(rep_msa->nseqs * sizeof(char*));
  for (i = 0; i < rep_msa->nseqs; i++) 
    pmsa->pooled_msa->names[i] = copy_charstr(rep_msa->names[i]);
  if (ncats >= 0) pmsa->pooled_msa->ncats = ncats;
  pmsa->lens = smalloc(lst_size(source_msas) * sizeof(int));

  pmsa->tuple_idx_map = smalloc(lst_size(source_msas) * sizeof(int*));
  key = smalloc((rep_msa->nseqs * tuple_size + 1) * sizeof(char));
  key[rep_msa->nseqs * tuple_size] = '\0';
  for (i = 0; i < lst_size(source_msas); i++) {
    MSA *smsa = (MSA*)lst_get_ptr(source_msas, i);
    checkInterrupt();
    if (smsa->nseqs != rep_msa->nseqs)
      die("ERROR: All MSA's must contain the same number of species! The offending sequence is: %s\n", pmsa->pooled_msa->names[i]);
    if (smsa->ss == NULL)
      ss_from_msas(smsa, tuple_size, 1, cats_to_do, NULL, NULL, -1, non_overlapping);
                                /* assume we want ordered suff stats
                                   for source alignments */
    ss_from_msas(pmsa->pooled_msa, tuple_size, 0, cats_to_do, 
                 smsa, tuple_hash, -1, non_overlapping);
    pmsa->lens[i] = smsa->length;

    /* keep around a mapping from the tuple indices of each source
       alignment to those of the pooled alignment, to allow a unified
       indexing scheme to be used */
    pmsa->tuple_idx_map[i] = smalloc(smsa->ss->ntuples * sizeof(int));
    for (j = 0; j < smsa->ss->ntuples; j++) {
      strncpy(key, smsa->ss->col_tuples[j], (smsa->nseqs * tuple_size + 1));
/*       fprintf(stderr, "i %d, j %d, key %s, tuple_idx[i][j] %d\n", */
/* 	      i, j, key, (int)hsh_get(tuple_hash, key)); */
      pmsa->tuple_idx_map[i][j] = ss_lookup_coltuple(key, tuple_hash, smsa);
      if (pmsa->tuple_idx_map[i][j] < 0)
	die("ERROR ss_pooled_from_msas: pmsa->tuple_idx_map[%i][%i]=%i, should be >= 0\n",
	    i, j, pmsa->tuple_idx_map[i][j]);
    }
  }
  hsh_free(tuple_hash);
  sfree(key);
  return pmsa;
}

/* Free a PooledMSA.  Does not free source alignments. */
void ss_free_pooled_msa(PooledMSA *pmsa) {
  int i;
  msa_free(pmsa->pooled_msa);
  for (i = 0; i < lst_size(pmsa->source_msas); i++)
    sfree(pmsa->tuple_idx_map[i]);
  sfree(pmsa->tuple_idx_map);
  sfree(pmsa->lens);
  sfree(pmsa);
}

/* Create an aggregate MSA from a list of MSA filenames, a list of
   sequence names, and an alphabet.  The list 'seqnames' will be used
   to define the order and contents of the sequences in the aggregate
   MSA (missing sequences will be replaced with missing data).  All source
   MSAs must share the same alphabet, and each must contain a subset
   of the names in 'seqnames'.  This function differs from the one
   above in that no direct representation is retained of the source
   MSAs.  Also, no information about tuple order is retained.  If
   cycle_size > 0, site categories will be labeled
   1,2,...,<cycle_size>,...,1,2,...,<cycle_size>. */
/* TODO: support collection of ordered sufficient stats -- possible
   now using idx_offset arg of ss_from_msas.  See maf_read in maf.c
   and warning message in msa_view.c. */
MSA *ss_aggregate_from_files(List *fnames,
                             List *seqnames, int tuple_size, List *cats_to_do, 
                             int cycle_size) {

  MSA *retval;
  Hashtable *tuple_hash = hsh_new(100000);
  int nseqs = lst_size(seqnames);
  MSA *source_msa = NULL;
  int i, j;
  FILE *F;
  msa_format_type format;
  char **names = (char**)smalloc(nseqs * sizeof(char*));

  /* set up names for new MSA object */
  for (i = 0; i < nseqs; i++) {
    String *s = lst_get_ptr(seqnames, i);
    names[i] = copy_charstr(s->chars);
  }

  retval = msa_new(NULL, names, nseqs, 0, NULL);
  retval->ncats = cycle_size > 0 ? cycle_size : -1;

  for (i = 0; i < lst_size(fnames); i++) {
    String *fname = lst_get_ptr(fnames, i);
    checkInterrupt();
    fprintf(stderr, "Reading alignment from %s ...\n", fname->chars);

    F = phast_fopen(fname->chars, "r");
    format = msa_format_for_content(F, 1);
    if (format == MAF)
      source_msa = maf_read_cats(F, NULL, tuple_size, NULL, NULL, NULL, cycle_size, 
                            FALSE, NULL, NO_STRIP, FALSE, cats_to_do); 
                                /* note: assuming unordered, not allowing
                                   overlapping blocks */
    else 
      source_msa = msa_new_from_file_define_format(F, format, NULL);
    phast_fclose(F);

    if (source_msa->seqs == NULL && 
        source_msa->ss->tuple_size != tuple_size) 
      die("ERROR: tuple size of input file '%s' (%d) does not match desired tuple size (%d).\n", fname->chars, source_msa->ss->tuple_size, tuple_size);

    if (cycle_size > 0 && format != MAF) { /* already done if MAF */
      source_msa->categories = smalloc(source_msa->length * sizeof(int));
      source_msa->ncats = cycle_size;
      for (j = 0; j < source_msa->length; j++)
        source_msa->categories[j] = (j % cycle_size) + 1;
    }

    if (source_msa->ncats != retval->ncats) { /* only an issue with SS */
      if (i == 0) retval->ncats = source_msa->ncats;
      else die("ERROR: input alignments have different numbers of categories.\n");
    }

    /* reorder the seqs and names; add seqs of missing data as necessary */
    msa_reorder_rows(source_msa, seqnames);

    /* now add the source MSA to the aggregate */
    ss_from_msas(retval, tuple_size, 0, cats_to_do, source_msa, 
                 tuple_hash, -1, 0);

    msa_free(source_msa);
  }

  hsh_free(tuple_hash);
  return retval;
}


char *ss_get_one_seq(MSA *msa, int spec) {
  char *seq, c;
  int i, j, col;
  
  seq = (char*)smalloc((msa->length+1)*sizeof(char));
  seq[msa->length] = '\0';
  if (msa->ss->tuple_idx == NULL) { /*unordered sufficient stats */
    col = 0;
    for (i=0; i<msa->ss->ntuples; i++) {
      checkInterruptN(i, 1000);
      c = col_string_to_char(msa, msa->ss->col_tuples[i], spec,
			     msa->ss->tuple_size, 0);
      while (col + msa->ss->counts[i] > msa->length) {  
	//this shouldn't happen, but the length isn't necessarily initialized
	//when SS is created
	msa->length *= 2;
	seq = srealloc(seq, (msa->length*sizeof(char)));
      }
      for (j=0; j<msa->ss->counts[i]; j++)
	seq[col++] = c;
    }
    seq[col] = '\0';
    if (col < msa->length)
      seq = srealloc(seq, (col+1)*sizeof(char));
  } else { /* ordered sufficient stats */
    for (col = 0; col < msa->length; col++) {
      seq[col] = col_string_to_char(msa, 
				    msa->ss->col_tuples[msa->ss->tuple_idx[col]], 
				    spec, msa->ss->tuple_size, 0);
    }
  }
  return seq;
}
				   
    
  

/* Reconstruct sequences from suff stats.  Only uses right-most col in
   tuple; requires length and nseqs to be correct, seqs to be NULL,
   categories to be NULL; does not set category labels (impossible
   from sufficient stats as they are currently defined) */
/* FIXME: this will be wrong if suff stats only describe certain
   categories and a higher-order model is used */
void ss_to_msa(MSA *msa) {
  int i;

  if (!(msa->seqs == NULL && msa->categories == NULL))
    die("ERROR ss_to_msa: msa->seqs and msa->categories should be NULL\n");
  msa->seqs = (char**)smalloc(msa->nseqs*sizeof(char*));
  for (i = 0; i < msa->nseqs; i++) 
    msa->seqs[i] = ss_get_one_seq(msa, i);
  msa->alloc_len = msa->length;
}


/* NOTE: assumes one file per species ... is this reasonable? */
/* this function should eventually be moved to msa.c */
void msa_read_AXT(MSA *msa, List *axt_fnames) {
  FILE *F;
  String *line, *ref, *targ;
  int i, j, k, start, line_no;
  List *fields;

  msa->nseqs = lst_size(axt_fnames)+1;
  msa->names = (char**)srealloc(msa->names, msa->nseqs * sizeof(char*));
  msa->seqs = (char**)srealloc(msa->seqs, msa->nseqs * sizeof(char*));

  line = str_new(STR_MED_LEN);
  ref = str_new(STR_VERY_LONG_LEN);
  targ = str_new(STR_VERY_LONG_LEN);
  fields = lst_new_ptr(9);
  for (i = 0; i < lst_size(axt_fnames); i++) {
    String *axtfname = (String*)lst_get_ptr(axt_fnames, i);

    msa->names[i+1] = (char*)smalloc(STR_MED_LEN * sizeof(char));
    msa->seqs[i+1] = (char*)smalloc((msa->length+1) * sizeof(char));
    for (j = 0; j < msa->length; j++) msa->seqs[i+1][j] = GAP_CHAR;
    msa->seqs[i+1][msa->length] = '\0';
    strcpy(msa->names[i+1], axtfname->chars); /* ?? */

    F = phast_fopen(axtfname->chars, "r");

    line_no=0;
    /* FIXME: need to deal with strand!  Also, soft masking ... */
    while (str_readline(line, F) != EOF) {
      checkInterruptN(line_no, 1000);
      line_no++;
      str_trim(line);

      if (line->length == 0) continue;

      str_split(line, NULL, fields);
      str_as_int(lst_get_ptr(fields, 2), &start); /* FIXME check error */
      for (j = 0; j < lst_size(fields); j++) 
        str_free((String*)lst_get_ptr(fields, j));

      if (str_readline(ref, F) == EOF || str_readline(targ, F) == EOF)
        break;

      str_trim(ref); str_trim(targ);

      k = start;
      for (j = 0; j < ref->length; j++) {
        if (ref->chars[j] != GAP_CHAR) {
          msa->seqs[i+1][k] = targ->chars[j];
          k++;
        }
      }
    }
  }
  str_free(line);
  str_free(ref);
  str_free(targ);
  lst_free(fields);
}

/* write a representation of an alignment in terms of its sufficient
   statistics */
void ss_write(MSA *msa, FILE *F, int show_order) {
  MSA_SS *ss = msa->ss;
  char tmp[ss->tuple_size * msa->nseqs + (ss->tuple_size-1) + 1];
  int i, j;
  String *namestr;
  tmp[ss->tuple_size * msa->nseqs + (ss->tuple_size-1) + 1] = '\0';

  namestr = str_new(STR_MED_LEN);
  for (i = 0; i < msa->nseqs; i++) {
    str_append_charstr(namestr, msa->names[i]);
    if (i < msa->nseqs-1) str_append_charstr(namestr, ",");
  }
                    
  fprintf(F, "NSEQS = %d\nLENGTH = %u\nTUPLE_SIZE = %d\nNTUPLES = %d\nNAMES = %s\nALPHABET = %s\n", 
          msa->nseqs, msa->length, ss->tuple_size, ss->ntuples, 
          namestr->chars, msa->alphabet);
                                /* NOTE: it's important to print
                                   LENGTH as unsigned; can overflow
                                   with whole mammalian genomes */
  if (msa->idx_offset != 0) fprintf(F, "IDX_OFFSET = %d\n", msa->idx_offset);
  fprintf(F, "NCATS = %d\n\n", msa->ncats);
  str_free(namestr);  

  for (i = 0; i < ss->ntuples; i++) {
    checkInterruptN(i, 100);
    tuple_to_string_pretty(tmp, msa, i);
    fprintf(F, "%d\t%s\t%.0f", i, tmp, ss->counts[i]);
    if (msa->ncats > 0 && ss->cat_counts != NULL) 
      for (j = 0; j <= msa->ncats; j++) 
        fprintf(F, "\t%.0f", ss->cat_counts[j][i]);
    fprintf(F, "\n");
  }
  if (show_order && ss->tuple_idx != NULL) {
    fprintf(F, "\nTUPLE_IDX_ORDER:\n");
    for (i = 0; i < msa->length; i++) {
      checkInterruptN(i, 100);
      fprintf(F, "%d\n", ss->tuple_idx[i]);
    }
  }
}

/* make reading order optional?  alphabet argument overrides alphabet
   in file (use NULL to use version in file) */
MSA* ss_read(FILE *F, char *alphabet) {
  Regex *nseqs_re, *length_re, *tuple_size_re, *ntuples_re, *tuple_re, 
    *names_re, *alph_re, *ncats_re, *order_re, *offset_re;
  String *line, *alph = NULL;
  int nseqs, length, tuple_size, ntuples, i, ncats = -99, header_done = 0, 
    idx_offset = 0, idx, offset, line_no=0;
  MSA *msa = NULL;
  List *matches;
  char **names = NULL;

  nseqs_re = str_re_new("NSEQS[[:space:]]*=[[:space:]]*([0-9]+)");
  length_re = str_re_new("LENGTH[[:space:]]*=[[:space:]]*([0-9]+)");
  tuple_size_re = str_re_new("TUPLE_SIZE[[:space:]]*=[[:space:]]*([0-9]+)");
  ntuples_re = str_re_new("NTUPLES[[:space:]]*=[[:space:]]*([0-9]+)");
  names_re = str_re_new("NAMES[[:space:]]*=[[:space:]]*(.*)");
  alph_re = str_re_new("ALPHABET[[:space:]]*=[[:space:]]*([-.^A-Za-z ]+)");
  ncats_re = str_re_new("NCATS[[:space:]]*=[[:space:]]*([-0-9]+)");
  offset_re = str_re_new("IDX_OFFSET[[:space:]]*=[[:space:]]*([-0-9]+)");
  tuple_re = str_re_new("^([0-9]+)[[:space:]]+([-.^A-Za-z ]+)[[:space:]]+([0-9.[:space:]]+)");
  order_re = str_re_new("TUPLE_IDX_ORDER:");

  line = str_new(STR_MED_LEN);
  matches = lst_new_ptr(3);
  nseqs = length = tuple_size = ntuples = -1;

  while (str_readline(line, F) != EOF) {
    checkInterruptN(line_no, 1000);
    line_no++;
    str_trim(line);
    if (line->length == 0) continue;
    if (line->chars[0]=='#') continue;

    if (!header_done) {
      if (str_re_match(line, nseqs_re, matches, 1) >= 0) {
        str_as_int(lst_get_ptr(matches, 1), &nseqs);
      }
      else if (str_re_match(line, length_re, matches, 1) >= 0) {
        str_as_int(lst_get_ptr(matches, 1), &length);
      }
      else if (str_re_match(line, tuple_size_re, matches, 1) >= 0) {
        str_as_int(lst_get_ptr(matches, 1), &tuple_size);
      }
      else if (str_re_match(line, ntuples_re, matches, 1) >= 0) {
        str_as_int(lst_get_ptr(matches, 1), &ntuples);
      }
      else if (str_re_match(line, alph_re, matches, 1) >= 0) {
        alph = str_dup(lst_get_ptr(matches, 1));
        str_remove_all_whitespace(alph);
      }
      else if (str_re_match(line, ncats_re, matches, 1) >= 0) {
        str_as_int(lst_get_ptr(matches, 1), &ncats);
        if (ncats < -1) ncats = -1;
      }
      else if (str_re_match(line, offset_re, matches, 1) >= 0) {
        str_as_int(lst_get_ptr(matches, 1), &idx_offset);
        if (idx_offset < -1) idx_offset = -1;
      }
      else if (str_re_match(line, names_re, matches, 1) >= 0) {
        List *names_list = lst_new_ptr(nseqs > 0 ? nseqs : 20);
        str_split(lst_get_ptr(matches, 1), ", ", names_list);
        names = (char**)smalloc(lst_size(names_list) * sizeof(char*));
        for (i = 0; i < lst_size(names_list); i++) {
          String *s = (String*)lst_get_ptr(names_list, i);
          names[i] = copy_charstr(s->chars);
          str_free(s);
        }
        lst_free(names_list);
      }      
      else 
        die("ERROR: unrecognized line in sufficient statistics file.  Is your header information complete?\nOffending line: '%s'\n", line->chars);

      if (nseqs > 0 && length >= 0 && tuple_size > 0 && ntuples > 0 && 
          alph != NULL && names != NULL && ncats != -99) {
        msa = msa_new(NULL, names, nseqs, length, alphabet != NULL ? alphabet : alph->chars);
                                /* allow alphabet from file to be overridden */
        if (ncats > 0) msa->ncats = ncats;
        msa->idx_offset = idx_offset;
        ss_new(msa, tuple_size, ntuples, (ncats > 0 ? 1 : 0), 0);
        msa->ss->ntuples = ntuples;
        /* in this case, we can preallocate the col_tuples, because we
           know exactly how many there will be */
/* 	fprintf(stderr, "ntuples %d, nseqs %d, tuple_size %d, alloc_size %d\n", */
/* 		ntuples, nseqs, tuple_size, (nseqs * tuple_size + 1)); */
        for (i = 0; i < ntuples; i++) {
          msa->ss->col_tuples[i] = (char*)smalloc((nseqs * tuple_size + 1) *
						   sizeof(char));
	  msa->ss->col_tuples[i][nseqs * tuple_size] = '\0';
	}

        str_free(alph);
        header_done = 1;
      }
    }
    
    else if (str_split(line, NULL, matches) >= 3) {
      String *tmpstr;

      tmpstr = lst_get_ptr(matches, 0);
      str_as_int(tmpstr, &idx);
      if (idx < 0 || idx >= ntuples) 
        die("ERROR: tuple line has index out of bounds.\nOffending line is, \"%s\"\n", line->chars);

      for (offset = -1 * (msa->ss->tuple_size-1); offset <= 0; offset++) {
        tmpstr = lst_get_ptr(matches, offset + msa->ss->tuple_size);
        if (tmpstr->length != msa->nseqs) 
          die("ERROR: length of column tuple does not match NSEQS.\nOffending line is, \"%s\"\n", line->chars);

/* 	fprintf(stderr, "tuple %d from input msa: %s\n", idx, tmpstr->chars); */
        for (i = 0; i < msa->nseqs; i++) {
          set_col_char_in_string(msa, msa->ss->col_tuples[idx], i,
                                 msa->ss->tuple_size, offset,
                                 tmpstr->chars[i]);
        }
/* 	fprintf(stderr, "tuple from msa in memory: %s\n", msa->ss->col_tuples[idx]); */
      }

      for (i = -1; i <= ncats; i++) { 
                                /* i == -1 -> global count; other
                                   values correspond to particular
                                   categories */
        tmpstr = lst_get_ptr(matches, i+msa->ss->tuple_size+2);
        str_as_dbl(tmpstr, i == -1 ? &msa->ss->counts[idx] : 
                   &msa->ss->cat_counts[i][idx]);
      }
    }
    else if (str_re_match(line, order_re, NULL, 0) >= 0) {
      msa->ss->tuple_idx = smalloc(msa->length * sizeof(int));
      for (i = 0; str_readline(line, F) != EOF && i < msa->length; i++) {
        str_trim(line);
        if (line->length == 0) continue;
        if (str_as_int(line, &msa->ss->tuple_idx[i]) != 0) 
          die("ERROR: bad integer in TUPLE_IDX_ORDER list.\n");
      }
      if (i < msa->length) 
        die("ERROR: too few numbers in TUPLE_IDX_ORDER list.\n");
    }

    for (i = 0; i < lst_size(matches); i++)
      str_free(lst_get_ptr(matches, i));
    lst_clear(matches);
  }

  if (!header_done || msa == NULL)
    die("ERROR: Missing or incomplete header in SS file.\n");

  lst_free(matches);
  str_re_free(nseqs_re);
  str_re_free(length_re);
  str_re_free(tuple_size_re);
  str_re_free(ntuples_re);
  str_re_free(names_re);
  str_re_free(alph_re);
  str_re_free(ncats_re);
  str_re_free(offset_re);
  str_re_free(tuple_re);
  str_re_free(order_re);
  str_free(line);
  
/*   for (idx = 0; idx < ntuples; idx++) */
/*     fprintf(stderr, "Tuple %d in msa: %s\n", idx, msa->ss->col_tuples[idx]); */
  
  return msa;
}

void ss_free_categories(MSA_SS *ss) {
  int j;
  if (ss->cat_counts != NULL) {
    for (j = 0; j < ss->msa->ncats; j++)
      sfree(ss->cat_counts[j]);
    sfree(ss->cat_counts);
    ss->cat_counts = NULL;
  }
}


/* free all memory associated with a sufficient stats object */
void ss_free(MSA_SS *ss) {
  int j;
  for (j = 0; j < ss->alloc_ntuples; j++)
    sfree(ss->col_tuples[j]);
  sfree(ss->col_tuples);
  ss_free_categories(ss);
  if (ss->counts != NULL) sfree(ss->counts);
  if (ss->tuple_idx != NULL) sfree(ss->tuple_idx);
  sfree(ss);
}

/* update category counts, according to 'categories' attribute of MSA
   object.  Requires ordered sufficient stats.  Will allocate and
   initialize cat_counts if necessary. */
void ss_update_categories(MSA *msa) {
  MSA_SS *ss = msa->ss;
  int i, j;
  if (!(msa->ncats >= 0 && msa->categories != NULL && 
	ss->tuple_idx != NULL))
    die("ERROR ss_update_categories: msa->ncats=%i msa->categories==NULL=%i, ss->tuple_idx==NULL=%i\n", msa->ncats, msa->categories==NULL, ss->tuple_idx==NULL);
  if (ss->cat_counts == NULL) {
    ss->cat_counts = (double**)smalloc((msa->ncats+1) * sizeof(double*));
    for (i = 0; i <= msa->ncats; i++) 
      ss->cat_counts[i] = (double*)smalloc(ss->ntuples * sizeof(double));
  }
  /* otherwise assume allocated correctly (reasonable?) */
  for (i = 0; i <= msa->ncats; i++) 
    for (j = 0; j < ss->ntuples; j++) 
      ss->cat_counts[i][j] = 0;
  for (i = 0; i < msa->length; i++) {
    checkInterruptN(i, 10000);
    if (msa->categories[i] > msa->ncats)
      die("ERROR ss_update_categories: msa->categories[%i]=%i, should be <= msa->ncats (%i)\n", i, msa->categories[i], msa->ncats);
    ss->cat_counts[msa->categories[i]][ss->tuple_idx[i]]++;
  }
}

/* Shrinks arrays to size ss->ntuples. */
void ss_compact(MSA_SS *ss) {
  int j;
  ss->col_tuples = (char**)srealloc(ss->col_tuples, 
                                    ss->ntuples*sizeof(char*));
  ss->counts = (double*)srealloc(ss->counts, 
                                ss->ntuples*sizeof(double));
  for (j = 0; ss->cat_counts != NULL && j <= ss->msa->ncats; j++)
    ss->cat_counts[j] = (double*)srealloc(ss->cat_counts[j], 
                                         ss->ntuples*sizeof(double));
  ss->alloc_ntuples = ss->ntuples;
}

/* given an MSA (with or without suff stats), create an alternative
   representation with sufficient statistics of a different tuple
   size.  The new alignment will share the seqs and names of the old
   alignment, but will have a new sufficient statistics object.  The
   'col_offset' argument allows category labels to be shifted *left*
   by the specified amount, as is convenient when estimating
   conditional probabilities with a two-pass approach.  */
MSA* ss_alt_msa(MSA *orig_msa, int new_tuple_size, int store_order, 
                int col_offset) {
  int i;
  MSA *new_msa = NULL;

  if (col_offset < 0)
    die("ERROR ss_alt_msa: col_offset=%i, should be >=0\n", col_offset);
  if (new_tuple_size <= 0)
    die("ERROR ss_alt_msa: new_tuple_size=%i, should be >=0\n", new_tuple_size);

  /* for now, must have sequences, so reconstruct them if all we have
     are suff stats.  This could be done more efficiently, but it's
     currently not worth the trouble */
  if (orig_msa->seqs == NULL) ss_to_msa(orig_msa);

  new_msa = msa_new(orig_msa->seqs, orig_msa->names, orig_msa->nseqs,
                    orig_msa->length, orig_msa->alphabet);
  new_msa->ncats = orig_msa->ncats;

  if (orig_msa->categories != NULL) {
    new_msa->categories = (int*)smalloc(orig_msa->length * sizeof(int));
    for (i = 0; i < orig_msa->length - col_offset; i++) 
      new_msa->categories[i] = orig_msa->categories[i+col_offset];
  }

  ss_from_msas(new_msa, new_tuple_size, store_order, NULL, NULL, NULL, -1, 0);
  return new_msa;
}

/* Extract a sub-alignment from an alignment, in terms of ordered
   sufficient statistics (refer to msa_sub_alignment in msa.c).  The
   new alignment will represent the interval [start_col, end_col). */
MSA *ss_sub_alignment(MSA *msa, char **new_names, List *include_list, 
                      int start_col, int end_col) {
  MSA *retval;
  int do_cats = (msa->ncats >= 0 && msa->categories != NULL);
  int i, offset, seqidx, tupidx, sub_ntuples, sub_tupidx, cat;
  int *full_to_sub;
  int unordered_seqs = (msa->ss->tuple_idx == NULL && 
                        start_col == 0 && end_col == msa->length);
                                /* special case: taking subset of seqs only
                                   with unordered sufficient stats */
  MSA_SS *ss;

  if (msa->ss == NULL)
    die("ERROR: sufficient stats required in ss_sub_alignment.\n");

  if (!unordered_seqs && msa->ss->tuple_idx == NULL) 
    die("ERROR: ordered sufficient statistics required in ss_sub_alignment.\n");

  if (msa->ncats >= 0 && !do_cats)
    fprintf(stderr, "WARNING: ss_sub_alignment can't handle site categories with categories vector.  Ignoring category-specific counts.\n");

  retval = msa_new(NULL, new_names, lst_size(include_list), 
                   end_col - start_col, msa->alphabet);
  if (do_cats) {
    retval->ncats = msa->ncats;
    retval->categories = smalloc(retval->length * sizeof(int));
  }

  /* mapping of original tuple numbers to tuple numbers in the
     sub-alignment */
  full_to_sub = smalloc(msa->ss->ntuples * sizeof(int));
  if (unordered_seqs) {         /* in this case, we know we still have
                                   all tuples.  (Some may be
                                   redundant, but that will be handled
                                   below.) */
    for (tupidx = 0; tupidx < msa->ss->ntuples; tupidx++) 
      full_to_sub[tupidx] = 0; /* placeholder */
    sub_ntuples = msa->ss->ntuples;
  }
  else {
    for (tupidx = 0; tupidx < msa->ss->ntuples; tupidx++) 
      full_to_sub[tupidx] = -1; /* indicates absent from subalignment */
    sub_ntuples = 0;
    for (i = 0; i < retval->length; i++) {
      checkInterruptN(i, 1000);
      if (!(msa->ss->tuple_idx[i+start_col] >= 0 && 
	    msa->ss->tuple_idx[i+start_col] < msa->ss->ntuples))
	die("ERROR: ss_sub_alignment: msa->ss->tuple_idx[%i]=%i, should be in [0, %i)\n", i+start_col, msa->ss->tuple_idx[i+start_col], msa->ss->ntuples);
      if (full_to_sub[msa->ss->tuple_idx[i+start_col]] == -1) {
        full_to_sub[msa->ss->tuple_idx[i+start_col]] = 0; /* placeholder */
        sub_ntuples++;
      }
    }
  }

  ss_new(retval, msa->ss->tuple_size, sub_ntuples, do_cats, 
         msa->ss->tuple_idx != NULL);
  ss = retval->ss;
  ss->ntuples = sub_ntuples;

  /* copy column tuples for specified seqs */
  for (tupidx = sub_tupidx = 0; tupidx < msa->ss->ntuples; tupidx++) {
    checkInterruptN(tupidx, 1000);
    if (full_to_sub[tupidx] == -1) continue;

    ss->col_tuples[sub_tupidx] = smalloc((retval->nseqs * ss->tuple_size + 1) * 
                                         sizeof(char));
    for (offset = -(ss->tuple_size-1); offset <= 0; offset++) {
      for (i = 0; i < lst_size(include_list); i++) {
        seqidx = lst_get_int(include_list, i);
	ss->col_tuples[sub_tupidx][ss->tuple_size*i + ss->tuple_size-1 + offset] =
	  msa->ss->col_tuples[tupidx][msa->ss->tuple_size*seqidx + msa->ss->tuple_size-1 + offset];
      }
    }
    full_to_sub[tupidx] = sub_tupidx++;
  }
  
  if (sub_tupidx != sub_ntuples)
    die("ERROR ss_sub_alignment: sub_tupidx (%i) != sub_ntuples (%i)\n",
	sub_tupidx, sub_ntuples);

  /* copy ordering info for specified columns and recompute counts.
     Also take care of category labels and category-specific counts */

  if (unordered_seqs) {         /* in this case, just copy counts
                                   directly tuple by tuple */
    for (i = 0; i < msa->ss->ntuples; i++) {
      checkInterruptN(i, 1000);
      ss->counts[full_to_sub[i]] = msa->ss->counts[i];
      if (do_cats) {
        for (cat = 0; cat <= msa->ncats; cat++)
          ss->cat_counts[cat][full_to_sub[i]] = msa->ss->cat_counts[cat][i];
      }
    }
  }
  else {                        /* go site by site */
    for (i = 0; i < retval->length; i++) {
      checkInterruptN(i, 1000);
      ss->tuple_idx[i] = full_to_sub[msa->ss->tuple_idx[i+start_col]];
      if (ss->tuple_idx[i] < 0) 
	die("ERROR ss_sub_alignment: ss->tuple_idx[%i]=%i, should be >=0\n", 
	    i, ss->tuple_idx);
      ss->counts[ss->tuple_idx[i]]++;
      if (do_cats) {
        retval->categories[i] = msa->categories[i+start_col];
        ss->cat_counts[retval->categories[i]][ss->tuple_idx[i]]++;
      }
    }  
  }

  if (lst_size(include_list) != msa->nseqs) 
    ss_unique(retval);

  sfree(full_to_sub);
  return retval;
}


/* adjust sufficient statistics to reflect the reverse complement of
   an alignment.  Refer to msa_reverse_compl */
void ss_reverse_compl(MSA *msa) {
  int i, j, k, offset1, offset2, midpt, do_cats, idx1, idx2;
  char c1, c2;
  MSA_SS *ss;
  char *first_tuple;
  char new_tuple[msa->ss->tuple_size * msa->nseqs + 1];
  Queue *overwrites = que_new_int(msa->ss->tuple_size);

  if (msa->ss == NULL || msa->ss->tuple_idx == NULL)
    die("ERROR ss_reverse_compl: Need ordered sufficient statistics\n");
  ss = msa->ss;

  if (msa->categories == NULL && ss->cat_counts != NULL)
    fprintf(stderr, "WARNING: ss_reverse_compl cannot address category-specific counts without a\ncategories vector.  Ignoring category counts.  They will be wrong!\n");

  do_cats = (msa->categories != NULL && ss->cat_counts != NULL);

  /* adjust counts for first few columns; these can't be reverse
     complemented */
  for (i = 0; i < ss->tuple_size - 1; i++) {
    checkInterruptN(i, 1000);
    ss->counts[ss->tuple_idx[i]]--;
    if (do_cats) ss->cat_counts[msa->categories[i]][ss->tuple_idx[i]]--;
    if (ss->counts[ss->tuple_idx[i]] == 0) 
      que_push_int(overwrites, ss->tuple_idx[i]); 
    /* it's likely that these leading tuples will not be present
       elsewhere in the alignment, so we'll make them available for
       overwriting (see below) */
  }

  /* reverse complement column tuples; counts remain unaltered */
  midpt = (int)ceil(ss->tuple_size/2.0);
  for (i = 0; i < ss->ntuples; i++) {
    checkInterruptN(i, 1000);
    for (j = 0; j < msa->nseqs; j++) {
      for (k = 0; k < midpt; k++) {
        offset1 = -(ss->tuple_size-1) + k;
        offset2 = -k;

        c1 = msa_compl_char(col_string_to_char(msa, ss->col_tuples[i],
                                               j, ss->tuple_size, offset1));

        c2 = msa_compl_char(col_string_to_char(msa, ss->col_tuples[i],
                                               j, ss->tuple_size, offset2));

        set_col_char_in_string(msa, ss->col_tuples[i], j, ss->tuple_size, 
                               offset1, c2);

        if (offset1 != offset2) /* true iff at midpoint */
          set_col_char_in_string(msa, ss->col_tuples[i], j, ss->tuple_size, 
                                 offset2, c1);          
      }
    }
  }

  /* reverse order */
  midpt = (int)ceil(msa->length/2.0);
  idx1 = ss->tuple_size - 1; /* note offset due to edge effect */
  idx2 = msa->length - 1;
  for (; idx1 < idx2; idx1++, idx2--) {
    int tmp = ss->tuple_idx[idx2];
    ss->tuple_idx[idx2] = ss->tuple_idx[idx1];
    ss->tuple_idx[idx1] = tmp;
  }    

  /* also reverse categories, if necessary */
  if (do_cats) {
    idx1 = 0; idx2 = msa->length - 1;
    for (; idx1 < idx2; idx1++, idx2--) {
      int tmp = msa->categories[idx2];
      msa->categories[idx2] = msa->categories[idx1];
      msa->categories[idx1] = tmp;
    }
  }

  /* now add representation of initial columns of reverse compl */
  first_tuple = ss->col_tuples[ss->tuple_idx[ss->tuple_size-1]];
  new_tuple[ss->tuple_size * msa->nseqs] = '\0';
  for (i = 0; i < ss->tuple_size-1; i++) {
    int new_tuple_idx;
    checkInterruptN(i, 10000);
    for (j = 0; j < ss->tuple_size * msa->nseqs; j++) new_tuple[j] = GAP_CHAR;
    for (offset2 = -i; offset2 <= 0; offset2++) { /* offset in new_tuple */
      offset1 = offset2 + i - (ss->tuple_size - 1); /* offset in first_tuple */
      for (j = 0; j < msa->nseqs; j++)
        set_col_char_in_string(msa, new_tuple, j, ss->tuple_size, offset2, 
                               col_string_to_char(msa, first_tuple, j,
                                                  ss->tuple_size, offset1));
    }

    /* if there's a zero-count tuple left over from above, overwrite
       it.  Otherwise, let's assume new_tuple is not already present.
       In the worst case, there will be redundant tuples in the suff
       stats, which is probably not too serious. */
    if (!que_empty(overwrites)) new_tuple_idx = que_pop_int(overwrites);
    else {
      ss_realloc(msa, ss->tuple_size, ss->ntuples + 1, do_cats, 1);
      ss->col_tuples[ss->ntuples] = smalloc((ss->tuple_size * msa->nseqs+1)*
                                           sizeof(char));
      new_tuple_idx = ss->ntuples++;
    }

    strncpy(ss->col_tuples[new_tuple_idx], new_tuple, 
            (msa->nseqs * ss->tuple_size + 1));
    ss->counts[new_tuple_idx]++;
    if (do_cats) ss->cat_counts[msa->categories[i]][new_tuple_idx]++;

  }

  que_free(overwrites);
}


/* change sufficient stats to reflect reordered rows of an alignment --
   see msa_reorder_rows.  */
void ss_reorder_rows(MSA *msa, int *new_to_old, int new_nseqs) {
  int ts = msa->ss->tuple_size;
  char tmp[msa->nseqs * ts];
  int col_offset, j, tup;
  for (tup = 0; tup < msa->ss->ntuples; tup++) {
    checkInterruptN(tup, 10000);
    strncpy(tmp, msa->ss->col_tuples[tup], msa->nseqs * ts);
    if (new_nseqs > msa->nseqs) 
      msa->ss->col_tuples[tup] = srealloc(msa->ss->col_tuples[tup], 
                                          new_nseqs * ts);
    for (col_offset = -ts+1; col_offset <= 0; col_offset++) {
      for (j = 0; j < new_nseqs; j++) {
        if (new_to_old[j] >= 0)
	  msa->ss->col_tuples[tup][ts*j + ts - 1 + col_offset] = 
	    tmp[ts*new_to_old[j] + ts - 1 + col_offset];
        else
	  msa->ss->col_tuples[tup][ts*j + ts-1 + col_offset] = msa->missing[0];
      }
    }
  }
}

/** Eliminate all tuples with counts of zero.  Main counts must
    reflect cat counts (cat counts won't be checked).  Tuples with
    counts of zero are assumed not to appear in tuple_idx.  */
void ss_remove_zero_counts(MSA *msa) {
  int i, cat, new_ntuples = 0;
  int *old_to_new = smalloc(msa->ss->ntuples * sizeof(int));

  for (i = 0; i < msa->ss->ntuples; i++) {
    checkInterruptN(i, 10000);
    if (msa->ss->counts[i] > 0) {
      if (new_ntuples != i) {
        msa->ss->col_tuples[new_ntuples] = msa->ss->col_tuples[i];
                                /* note new_ntuples <= i, so
                                   tuple with idx new_tuples has
                                   already been visited (data was
                                   transferred or memory was freed) */
        msa->ss->counts[new_ntuples] = msa->ss->counts[i];
        if (msa->ss->cat_counts != NULL) 
          for (cat = 0; cat <= msa->ncats; cat++)
            msa->ss->cat_counts[cat][new_ntuples] = msa->ss->cat_counts[cat][i];
      }
      old_to_new[i] = new_ntuples++;
    }
    else { sfree(msa->ss->col_tuples[i]); msa->ss->col_tuples[i] = NULL; }
  }

  if (msa->ss->tuple_idx != NULL)
    for (i = 0; i < msa->length; i++)
      msa->ss->tuple_idx[i] = old_to_new[msa->ss->tuple_idx[i]];

  msa->ss->ntuples = new_ntuples;
  ss_compact(msa->ss);  
  sfree(old_to_new);
}

/** Ensure all tuples are unique.  Combine counts and remap indices as
    necessary. */
void ss_unique(MSA *msa) {
  char key[msa->nseqs * msa->ss->tuple_size + 1];
  Hashtable *hash = hsh_new(msa->ss->ntuples);
  int i, idx, cat;
  int *old_to_new = smalloc(msa->ss->ntuples * sizeof(int));
  key[msa->nseqs * msa->ss->tuple_size] = '\0';

  for (i = 0; i < msa->ss->ntuples; i++) {
    checkInterruptN(i, 10000);
    strncpy(key, msa->ss->col_tuples[i], (msa->nseqs * msa->ss->tuple_size + 1));
    if ((idx = hsh_get_int(hash, key)) == -1) { /* tuple not seen before */
      hsh_put_int(hash, key, i);
      old_to_new[i] = i;
    }
    else { /* seen before: combine counts with prev, zero out this version */
      msa->ss->counts[idx] += msa->ss->counts[i];
      if (msa->ss->cat_counts != NULL) 
        for (cat = 0; cat <= msa->ncats; cat++)
          msa->ss->cat_counts[cat][idx] += msa->ss->cat_counts[cat][i];
      msa->ss->counts[i] = 0;   /* no need to zero cat_counts */
      old_to_new[i] = idx;
    }
  }

  if (msa->ss->tuple_idx != NULL)
    for (i = 0; i < msa->length; i++)
      msa->ss->tuple_idx[i] = old_to_new[msa->ss->tuple_idx[i]];

  ss_remove_zero_counts(msa);
  sfree(old_to_new);
  hsh_free(hash);
}

/** Convert all missing data characters to the default missing data
    character (msa->missing[0]).  Optionally also convert gap
    characters.  Can be useful in reducing number of tuples */
void ss_collapse_missing(MSA *msa, int do_gaps) {
  int i, j, len = msa->nseqs * msa->ss->tuple_size;
  int changed_missing = FALSE, changed_gaps = FALSE, exists_missing = FALSE;
  for (i = 0; i < msa->ss->ntuples; i++) {
    checkInterruptN(i, 10000);
    for (j = 0; j < len; j++) {
      char c = msa->ss->col_tuples[i][j];
      if (!exists_missing && c == msa->missing[0]) exists_missing = TRUE;
      else if (c != msa->missing[0] && msa->is_missing[(int)c]) {
        msa->ss->col_tuples[i][j] = msa->missing[0];
        if (!changed_missing) changed_missing = TRUE;
      }
      else if (do_gaps && c == GAP_CHAR) {
        msa->ss->col_tuples[i][j] = msa->missing[0];
	if (!changed_gaps) changed_gaps = TRUE;
      }
    }
  }

  if (changed_missing || (exists_missing && changed_gaps))
      ss_unique(msa);
      /* avoid expensive call to ss_unique if there were no missing
	 data characters and all we did was relabel gaps as missing
	 data */
}


/* return the number of non-gap characters for a particular sequnce */
int ss_seqlen(MSA *msa, int seqidx) {
  int i, count=0;
  for (i=0; i<msa->ss->ntuples; i++)
    if (ss_get_char_tuple(msa, i, seqidx, 0) != GAP_CHAR)
      count += (int)msa->ss->counts[i];
  return count;
}


/** Like msa_strip_gaps, but in terms of suffient statistics.  If
   gap_strip_mode is STRIP_ALL_GAPS or STRIP_ANY_GAPS, removes all
   columns with ALL or ANY gaps, respectively.  Otherwise, assumes a
   *projection* is desired onto the sequence whose index is
   gap_strip_mode (indexing starts with 1).  Changes are made to
   original alignment.  Gaps are expected to be represented by
   GAP_CHAR.  If msa->categories is non-NULL, it will be adjusted
   accordingly. */
void ss_strip_gaps(MSA *msa, int gap_strip_mode) {
  int i, j;
  int newlen = msa->length;
  for (i = 0; i < msa->ss->ntuples; i++) {
    int strip;
    checkInterruptN(i, 10000);
    if (gap_strip_mode > 0)    /* project on refseq */
      strip = (ss_get_char_tuple(msa, i, gap_strip_mode-1, 0) == GAP_CHAR);
    else {                      /* stip columns with all or any gaps */
      strip = (gap_strip_mode == STRIP_ALL_GAPS);
      for (j = 0; j < msa->nseqs; j++) {
        char c = ss_get_char_tuple(msa, i, j, 0);
        if (gap_strip_mode == STRIP_ANY_GAPS && c == GAP_CHAR) {
          strip = TRUE; break;
        }
        else if (gap_strip_mode == STRIP_ALL_GAPS && c != GAP_CHAR) {
          strip = FALSE; break;
        }        
      }
    }
    if (strip) {
      newlen -= msa->ss->counts[i];
      msa->ss->counts[i] = 0;
    }
  }
  
  if (msa->ss->tuple_idx != NULL) {
    for (i = 0, j = 0; i < msa->length; i++) {
      checkInterruptN(j, 10000);
      if (msa->ss->counts[msa->ss->tuple_idx[i]] > 0) {
        msa->ss->tuple_idx[j] = msa->ss->tuple_idx[i];
        if (msa->categories != NULL)
          msa->categories[j] = msa->categories[i];
        j++;
      }
    }
    /* above assumes any tuple is to be removed that is present in
       tuple_idx but has a count of zero */

    if (newlen != j)
      die("ERROR ss_strip_gaps: newlen (%i) != j (%i)\n", newlen, j);
    msa->length = newlen;
  } else 
    msa_update_length(msa);

  ss_remove_zero_counts(msa);
}

/** Strip columns consisting of missing data in all
   sequences but the reference sequence.  If msa->categories is
   non-NULL, it will be adjusted accordingly */
/* FIXME: what about category counts?  what happens in ss_strip_gaps? */
void ss_strip_missing(MSA *msa, /* Input alignment; will be altered */
                      int refseq /* Index of reference sequence
                                    (indexing starts with 1) */
                      ) { 
  int i, j;
  int newlen = msa->length;

  for (i = 0; i < msa->ss->ntuples; i++) {
    int strip = TRUE;
    checkInterruptN(i, 1000);
    for (j = 0; j < msa->nseqs && strip; j++) {
      if (j == refseq-1 || 
          (msa->is_informative != NULL && !msa->is_informative[j]))
        continue;
      else if (!msa->is_missing[(int)ss_get_char_tuple(msa, i, j, 0)])
        strip = FALSE; 
    }
    if (strip) {
      newlen -= msa->ss->counts[i];
      msa->ss->counts[i] = 0;
    }
  }
  
  if (msa->ss->tuple_idx != NULL) {
    for (i = 0, j = 0; i < msa->length; i++) {
      checkInterruptN(i, 10000);
      if (msa->ss->counts[msa->ss->tuple_idx[i]] > 0) {
        msa->ss->tuple_idx[j] = msa->ss->tuple_idx[i];
        if (msa->categories != NULL)
          msa->categories[j] = msa->categories[i];
        j++;
      }
    }
    /* above assumes any tuple is to be removed that is present in
       tuple_idx but has a count of zero */

    if (newlen != j)
      die("ERROR ss_strip_missing: newlen (%i) != j (%i)\n", newlen, j);
    msa->length = newlen;
  }

  ss_remove_zero_counts(msa);
}

/* return 1 if a column triple represents a 4d site and 0 otherwise */
int ss_is_4d(MSA *msa, int tuple) {
  char base[2];
  int offset, seq;
  if (msa->ss->tuple_size != 3)
    die("ERROR ss_is_4d need tuple_size 3, got %i\n", msa->ss->tuple_size);
  for (offset = -2; offset < 0; offset++) {
    base[offset+2] = '\0';
    for (seq = 0; seq < msa->nseqs; seq++) {
      char c = ss_get_char_tuple(msa, tuple, seq, offset);
      if (msa->is_missing[(int)c]) continue;
      else if (c == GAP_CHAR || msa->inv_alphabet[(int)c] < 0) return FALSE;
      else if (base[offset+2] == '\0') base[offset+2] = c;
      else if (base[offset+2] != c) return FALSE;
    }
    if (base[offset+2] == '\0') return FALSE; /* in case column has
                                                 only missing data */
  }
  for (seq = 0; seq < msa->nseqs; seq++)
    if (ss_get_char_tuple(msa, tuple, seq, 0) == GAP_CHAR)
      return FALSE;             /* also prohibit gap in 3rd pos. */

  if (((base[0] == 'A' || base[0] == 'T') && base[1] == 'C') ||
      ((base[0] == 'C' || base[0] == 'G') && base[1] != 'A'))
    return TRUE;                /* first two bases defining 4d codons in
                                   standard genetic code */
  return FALSE;
}

/* reduce to smaller tuple size representation */
/* NOTE: this only works because we know that ss are represented in
   strings in such a way that the old tuples won't be over-written by 
   the new ones. */
void ss_reduce_tuple_size(MSA *msa, int new_tuple_size) {
  int i, j, k, newlen;
  if (new_tuple_size >= msa->ss->tuple_size)
    die("ERROR: new tuple size must be smaller than old in ss_reduce_tuple_size.\n");
  newlen = msa->nseqs * new_tuple_size;
  for (i = 0; i < msa->ss->ntuples; i++)  {
    checkInterruptN(i, 10000);
    for (j = 0; j < msa->nseqs; j++) 
      for (k= -new_tuple_size+1; k<=0; k++) {
	set_col_char_in_string(msa, msa->ss->col_tuples[i], j, new_tuple_size, k,
			       col_string_to_char(msa, msa->ss->col_tuples[i], j, msa->ss->tuple_size, k));
      }
    msa->ss->col_tuples[i][newlen] = '\0';
    msa->ss->col_tuples[i] = srealloc(msa->ss->col_tuples[i], (newlen+1)*sizeof(char));
  }
  msa->ss->tuple_size = new_tuple_size;
  ss_unique(msa);
}

int ss_lookup_coltuple(char *coltuple_str, Hashtable *tuple_hash, MSA *msa) {
  int allgap[msa->ss->tuple_size], i, j, rv, tuple_size, len;
  char tempchar;
  tuple_size = msa->ss->tuple_size;
  len = tuple_size * msa->nseqs;
  for (i=0; i < tuple_size; i++) {
    for (j=0; j<msa->nseqs; j++)
      if (coltuple_str[j*tuple_size + i]!=GAP_CHAR &&
	  coltuple_str[j*tuple_size + i]!=msa->missing[0]) break;
    allgap[i] = (j == msa->nseqs);
  }
  for (i=len-1; i>=0; i--) {
/*     fprintf(stderr, "i %d, coltuple_str %s, coltuple_str[i] %c, msa->missing[0] %c, allgap[i tuple_size] %d\n", */
/* 	    i, coltuple_str, coltuple_str[i], msa->missing[0], */
/* 	    allgap[i%tuple_size]); */
    if (coltuple_str[i] != msa->missing[0] && allgap[i%tuple_size]==0)
      break;
  }
  i++;
  while (i%msa->ss->tuple_size != 0) i++;
  tempchar = coltuple_str[i];
  coltuple_str[i] = '\0';
  rv = hsh_get_int(tuple_hash, coltuple_str);
  coltuple_str[i] = tempchar;
  return rv;
}

void ss_add_coltuple(char *coltuple_str, void *val, Hashtable *tuple_hash, 
		     MSA *msa) {
  int allgap[msa->ss->tuple_size], i, j, tuple_size, len;
  char tempchar;
  tuple_size = msa->ss->tuple_size;
  len = tuple_size * msa->nseqs;

  for (i=0; i < tuple_size; i++) {
    for (j=0; j<msa->nseqs; j++)
      if (coltuple_str[j*tuple_size + i]!=GAP_CHAR &&
	  coltuple_str[j*tuple_size + i]!=msa->missing[0]) break;
    allgap[i] = (j==msa->nseqs);
  }

  for (i=len-1; i>=0; i--)
    if (coltuple_str[i] != msa->missing[0] && allgap[i%tuple_size]==0)
      break;
  i++;

  while (i%msa->ss->tuple_size != 0) i++;
  /* Until 1-27-2009, tempchar was not being stored! This bug took three months
     to find and wreaked havoc with dmsample! --agd27 */
  tempchar = coltuple_str[i];
  coltuple_str[i] = '\0';
  hsh_put(tuple_hash, coltuple_str, val);
  coltuple_str[i] = tempchar;
}


/* impose an artificial ordering on tuples if they aren't already ordered */
void ss_make_ordered(MSA *msa) {
  int i, j, idx=0, count;
  if (msa->ss == NULL || msa->ss->tuple_idx != NULL) return;
  msa_update_length(msa);
  msa->ss->tuple_idx = smalloc(msa->length*sizeof(int));
  for (i=0; i < msa->ss->ntuples; i++) {
    checkInterruptN(i, 10000);
    count = (int)msa->ss->counts[i];
    if (fabs(msa->ss->counts[i] - count) > 0.00001)
      die("can't impose order on alignment with non-integral counts");
    for (j=0; j < count; j++)
      msa->ss->tuple_idx[idx++] = i;
  }
  if (idx != msa->length)
    die("ERROR ss_make_ordered: idx (%i) != msa->length (%i)\n", 
	idx, msa->length);
}
