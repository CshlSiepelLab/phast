/* $Id: sufficient_stats.c,v 1.21 2005-05-28 21:31:06 acs Exp $
   Written by Adam Siepel, 2002 and 2003
   Copyright 2002, 2003, Adam Siepel, University of California */

#include "sufficient_stats.h"
#include "maf.h"
#include "assert.h"
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
   they exist for source alignments.  */

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
                  int idx_offset) {
  int i, j, do_cats, idx, upper_bound;
  long int max_tuples;
  MSA_SS *main_ss, *source_ss = NULL;
  Hashtable *tuple_hash = NULL;
  int *do_cat_number = NULL;
  char key[msa->nseqs * tuple_size + 1];
  MSA *smsa;
  int effective_offset = (idx_offset < 0 ? 0 : idx_offset); 
  if (source_msa == NULL && 
      (msa->seqs == NULL || msa->length <= 0 || msa->ss != NULL))
    die("ERROR: with no separate source alignment, ss_from_msas expects sequences of positive length and no SS object.\n");
  if (idx_offset >= 0) assert(store_order && source_msa != NULL);
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
    max_tuples = min(pow(strlen(msa->alphabet)+ strlen(msa->missing) + 1, 
                         msa->nseqs * tuple_size),
                     upper_bound);
    ss_new(msa, tuple_size, max_tuples, do_cats, store_order);
    if (source_msa != NULL) msa->length = source_msa->length;
  }
  else if (idx_offset < 0) {    /* if storing order based on source
                                   alignments and offset, then assume
                                   proper preallocation; otherwise
                                   (this case), realloc to accommodate
                                   new source msa */
    msa->length += source_msa->length;
    if (source_msa->ss != NULL) 
      upper_bound = msa->ss->ntuples + source_msa->ss->ntuples;
    else
      upper_bound = msa->ss->ntuples + source_msa->length;
    max_tuples = min(pow(strlen(msa->alphabet) + strlen(msa->missing) + 1, 
                         msa->nseqs * tuple_size),
                     upper_bound);
    ss_realloc(msa, tuple_size, max_tuples, do_cats, store_order);
  }

  main_ss = msa->ss;
  tuple_hash = existing_hash != NULL ? existing_hash : hsh_new(max_tuples/3);

  if (source_msa != NULL && source_msa->ss != NULL)
    source_ss = source_msa->ss;
  /* otherwise, source_ss will remain NULL */

  if (source_ss != NULL && !store_order) { /* in this case, just use
                                              existing suff stats */
    for (i = 0; i < source_ss->ntuples; i++) {
      strncpy(key, source_ss->col_tuples[i], msa->nseqs * tuple_size);
      if ((idx = (int)hsh_get(tuple_hash, key)) == -1) {
                                /* column tuple has not been seen
                                   before */
        idx = main_ss->ntuples++;
        hsh_put(tuple_hash, key, (void*)idx);
        main_ss->col_tuples[idx] = (char*)smalloc(tuple_size * msa->nseqs * 
                                                  sizeof(char));
        strncpy(main_ss->col_tuples[idx], key, msa->nseqs * tuple_size);
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
      if (do_cats && cats_to_do != NULL && 
          do_cat_number[smsa->categories[i]] == 0) {
        if (store_order) main_ss->tuple_idx[i + effective_offset] = -1;
        continue;
      }

      if (smsa->seqs != NULL)
        col_to_string(key, smsa, i, tuple_size);
      else                      /* NOTE: must have ordered suff stats */
        strncpy(key, smsa->ss->col_tuples[smsa->ss->tuple_idx[i]], 
                msa->nseqs * tuple_size);

      if ((idx = (int)hsh_get(tuple_hash, key)) == -1) {
                                /* column tuple has not been seen
                                   before */
        idx = main_ss->ntuples++;
        hsh_put(tuple_hash, key, (void*)idx);

        if (main_ss->ntuples > main_ss->alloc_ntuples) 
                                /* possible if allocated only for
                                   MAX_NTUPLE_ALLOC */
          ss_realloc(msa, tuple_size, main_ss->ntuples, do_cats, store_order);

        main_ss->col_tuples[idx] = (char*)smalloc(tuple_size * msa->nseqs * 
                                                 sizeof(char));
        strncpy(main_ss->col_tuples[idx], key, msa->nseqs * tuple_size);
      }

      main_ss->counts[idx]++;
      if (do_cats && smsa->categories != NULL) {
        assert(smsa->categories[i] >= 0 && smsa->categories[i] <= msa->ncats);
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

  if (do_cats) free(do_cat_number);
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
  if (store_order)
    ss->tuple_idx = (int*)smalloc(msa->length * sizeof(int));
  ss->col_tuples = (char**)smalloc(max_ntuples * sizeof(char*));
  for (i = 0; i < max_ntuples; i++) ss->col_tuples[i] = NULL;
  ss->counts = (double*)smalloc(max_ntuples * sizeof(double));
  for (i = 0; i < max_ntuples; i++) ss->counts[i] = 0; 
  if (do_cats) {
    assert(msa->ncats >= 0);
    ss->cat_counts = (double**)smalloc((msa->ncats+1) * sizeof(double*));
    for (j = 0; j <= msa->ncats; j++) {
      ss->cat_counts[j] = (double*)smalloc(max_ntuples * sizeof(double));
      for (i = 0; i < max_ntuples; i++) 
        ss->cat_counts[j][i] = 0;
    }
  }
  ss->alloc_len = msa->length;  /* see ss_realloc */
  ss->alloc_ntuples = max_ntuples;
}

/* ensures a suff stats object has enough memory allocated to
   accommodate new msa->length and max_ntuples */
void ss_realloc(MSA *msa, int tuple_size, int max_ntuples, int do_cats, 
                int store_order) {

  int i, j, cat_counts_done = FALSE;
  MSA_SS *ss = msa->ss;
  if (store_order && msa->length > ss->alloc_len) {
    ss->alloc_len = max(ss->alloc_len * 2, msa->length);
    ss->tuple_idx = (int*)srealloc(ss->tuple_idx, ss->alloc_len * sizeof(int));
  }
  if (do_cats && ss->cat_counts == NULL) {
    assert(msa->ncats >= 0);
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
                               List *cats_to_do) {
  int i, j;
  MSA *rep_msa;
  PooledMSA *pmsa = (PooledMSA*)smalloc(sizeof(PooledMSA));
  Hashtable *tuple_hash = hsh_new(100000); /* wild guess on size.  Big
                                              enough?  Preview files? */
  char *key;

  assert(lst_size(source_msas) > 0);
  rep_msa = (MSA*)lst_get_ptr(source_msas, 0);

  pmsa->pooled_msa = msa_new(NULL, NULL, rep_msa->nseqs, 0, rep_msa->alphabet);
  pmsa->source_msas = source_msas;
  pmsa->pooled_msa->names = (char**)smalloc(rep_msa->nseqs * sizeof(char*));
  for (i = 0; i < rep_msa->nseqs; i++) 
    pmsa->pooled_msa->names[i] = strdup(rep_msa->names[i]);
  if (ncats >= 0) pmsa->pooled_msa->ncats = ncats;
  pmsa->lens = smalloc(lst_size(source_msas) * sizeof(int));

  pmsa->tuple_idx_map = smalloc(lst_size(source_msas) * sizeof(int));
  key = smalloc((rep_msa->nseqs * tuple_size + 1) * sizeof(char));
  key[rep_msa->nseqs * tuple_size] = '\0';
  for (i = 0; i < lst_size(source_msas); i++) {
    MSA *smsa = (MSA*)lst_get_ptr(source_msas, i);
    assert(smsa->nseqs == rep_msa->nseqs);
    if (smsa->ss == NULL)
      ss_from_msas(smsa, tuple_size, 1, cats_to_do, NULL, NULL, -1);
                                /* assume we want ordered suff stats
                                   for source alignments */
    ss_from_msas(pmsa->pooled_msa, tuple_size, 0, cats_to_do, 
                 smsa, tuple_hash, -1);
    pmsa->lens[i] = smsa->length;

    /* keep around a mapping from the tuple indices of each source
       alignment to those of the pooled alignment, to allow a unified
       indexing scheme to be used */
    pmsa->tuple_idx_map[i] = smalloc(smsa->ss->ntuples * sizeof(int));
    for (j = 0; j < smsa->ss->ntuples; j++) {
      strncpy(key, smsa->ss->col_tuples[j], smsa->nseqs * tuple_size);      
      pmsa->tuple_idx_map[i][j] = (int)hsh_get(tuple_hash, key);
      assert(pmsa->tuple_idx_map[i][j] >= 0);
    }
  }
  hsh_free(tuple_hash);
  free(key);
  return pmsa;
}

/* Free a PooledMSA.  Does not free source alignments. */
void ss_free_pooled_msa(PooledMSA *pmsa) {
  int i;
  msa_free(pmsa->pooled_msa);
  for (i = 0; i < lst_size(pmsa->source_msas); i++)
    free(pmsa->tuple_idx_map[i]);
  free(pmsa->tuple_idx_map);
  free(pmsa->lens);
  free(pmsa);
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
MSA *ss_aggregate_from_files(List *fnames, msa_format_type format,
                             List *seqnames, int tuple_size, List *cats_to_do, 
                             int cycle_size) {

  MSA *retval;
  Hashtable *tuple_hash = hsh_new(100000);
  int nseqs = lst_size(seqnames);
  MSA *source_msa = NULL;
  int i, j;
  FILE *F;
  char **names = (char**)smalloc(nseqs * sizeof(char*));

  /* set up names for new MSA object */
  for (i = 0; i < nseqs; i++) {
    String *s = lst_get_ptr(seqnames, i);
    names[i] = strdup(s->chars);
  }

  retval = msa_new(NULL, names, nseqs, 0, NULL);
  retval->ncats = cycle_size > 0 ? cycle_size : -1;

  for (i = 0; i < lst_size(fnames); i++) {
    String *fname = lst_get_ptr(fnames, i);
    fprintf(stderr, "Reading alignment from %s ...\n", fname->chars);

    F = fopen_fname(fname->chars, "r");
    if (format == MAF)
      source_msa = maf_read(F, NULL, tuple_size, NULL, NULL, cycle_size, 
                            FALSE, NULL, NO_STRIP, FALSE); 
                                /* note: assuming unordered, not allowing
                                   overlapping blocks */
    else 
      source_msa = msa_new_from_file(F, format, NULL);
    fclose(F);

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
                 tuple_hash, -1);

    msa_free(source_msa);
  }

  hsh_free(tuple_hash);
  return retval;
}

/* Reconstruct sequences from suff stats.  Only uses right-most col in
   tuple; requires length and nseqs to be correct, seqs to be NULL,
   categories to be NULL; does not set category labels (impossible
   from sufficient stats as they are currently defined) */
/* FIXME: this will be wrong if suff stats only describe certain
   categories and a higher-order model is used */
void ss_to_msa(MSA *msa) {
  int i, j, k, col;
  char *colstr;

  assert(msa->seqs == NULL && msa->categories == NULL);
  msa->seqs = (char**)smalloc(msa->nseqs*sizeof(char*));
  for (i = 0; i < msa->nseqs; i++) {
    msa->seqs[i] = (char*)smalloc((msa->length+1) * sizeof(char));
    msa->seqs[i][msa->length] = '\0';
  }  

  if (msa->ss->tuple_idx == NULL) { /* unordered sufficient stats */
    col = 0;
    for (i = 0; i < msa->ss->ntuples; i++) {
      colstr = msa->ss->col_tuples[i];
      for (j = 0; j < msa->ss->counts[i]; j++) {
        for (k = 0; k < msa->nseqs; k++) {
          msa->seqs[k][col] = col_string_to_char(msa, colstr, k, 
                                                 msa->ss->tuple_size, 0);
        }
        col++;
      }
    }
  }
  else {                        /* ordered sufficient stats */
    for (col = 0; col < msa->length; col++) {
      colstr = msa->ss->col_tuples[msa->ss->tuple_idx[col]];
      for (k = 0; k < msa->nseqs; k++) {
        msa->seqs[k][col] = col_string_to_char(msa, colstr, k, 
                                               msa->ss->tuple_size, 0);
      }
    }
  }
}

/* NOTE: assumes one file per species ... is this reasonable? */
/* this function should eventually be moved to msa.c */
void msa_read_AXT(MSA *msa, List *axt_fnames) {
  FILE *F;
  String *line, *ref, *targ;
  int i, j, k, start;
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

    if ((F = fopen(axtfname->chars, "r")) == NULL) {
      fprintf(stderr, "ERROR: unable to open %s\n", axtfname->chars);
      exit(1);
    }

    /* FIXME: need to deal with strand!  Also, soft masking ... */
    while (str_readline(line, F) != EOF) {
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
/*           assert(msa->seqs[0][k] == ref->chars[j]); */
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

/* Produce a string representation of an alignment column tuple, given
   the model order; str must be allocated externally to size
   msa->nseqs * (tuple_size) + 1 */
void col_to_string(char *str, MSA *msa, int col, int tuple_size) {
  int col_offset, j;
  for (col_offset = -1 * (tuple_size-1); col_offset <= 0; col_offset++) 
  for (j = 0; j < msa->nseqs; j++) 
    str[msa->nseqs*(tuple_size-1 + col_offset) + j] = 
      (col+col_offset >= 0 ? msa->seqs[j][col+col_offset] : GAP_CHAR);
  str[msa->nseqs*tuple_size] = '\0';
}

/* Produce a printable representation of the specified tuple.  Strings
   representing each column will be separated by spaces */
void tuple_to_string_pretty(char *str, MSA *msa, int tupleidx) {
  int stridx = 0, offset, j;
  for (offset = -1 * (msa->ss->tuple_size-1); offset <= 0; offset++) {
    for (j = 0; j < msa->nseqs; j++) {
      str[stridx++] = col_string_to_char(msa, msa->ss->col_tuples[tupleidx], 
                                         j, msa->ss->tuple_size, offset);
    }
    if (offset < 0) str[stridx++] = ' ';
  }
  str[stridx] = '\0';
}

/* Given a string representation of a column tuple, return the
   character corresponding to the specified sequence and column.  Here
   sequence indexing begins with 0.  The index 'col_offset' is defined
   relative to the last column in the tuple.  If col_offset == 0 then
   the last column is considered, if col_offset == -1 then the
   preceding one is considered, and so on. */
char col_string_to_char(MSA *msa, char *str, int seqidx, int tuple_size, int col_offset) {
  return str[msa->nseqs*(tuple_size-1+col_offset) + seqidx];
}

void set_col_char_in_string(MSA *msa, char *str, int seqidx, 
                            int tuple_size, int col_offset, char c) {
  str[msa->nseqs*(tuple_size-1+col_offset) + seqidx] = c;
}

/* return character for specified sequence given index of tuple */
char ss_get_char_tuple(MSA *msa, int tupleidx, int seqidx, 
                       int col_offset) {
  return col_string_to_char(msa, msa->ss->col_tuples[tupleidx], seqidx, 
                            msa->ss->tuple_size, col_offset);
}

/* return character for specified sequence at specified alignment
   position; requires representation of column order */
char ss_get_char_pos(MSA *msa, int position, int seqidx,
                     int col_offset) {
  assert(msa->ss->tuple_idx != NULL);
  return col_string_to_char(msa, 
                            msa->ss->col_tuples[msa->ss->tuple_idx[position]], 
                            seqidx, msa->ss->tuple_size, col_offset);
}

/* fill out 'tuplestr' with tuple of characters present in specified
   sequence and column tuple; tuplestr must be allocated to at least
   tuple_size (null terminator will not be added, but if present will
   be left unchanged). */
void ss_get_tuple_of_chars(MSA *msa, int tupleidx, int seqidx, 
                           char *tuplestr) {
  int offset;
  for (offset = -1 * (msa->ss->tuple_size-1); offset <= 0; offset++) {
    tuplestr[msa->ss->tuple_size + offset - 1] = 
      col_string_to_char(msa, msa->ss->col_tuples[tupleidx], 
                         seqidx, msa->ss->tuple_size, offset);
  }
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
    tuple_to_string_pretty(tmp, msa, i);
    fprintf(F, "%d\t%s\t%.0f", i, tmp, ss->counts[i]);
    if (msa->ncats > 0 && ss->cat_counts != NULL) 
      for (j = 0; j <= msa->ncats; j++) 
        fprintf(F, "\t%.0f", ss->cat_counts[j][i]);
    fprintf(F, "\n");
  }
  if (show_order && ss->tuple_idx != NULL) {
    fprintf(F, "\nTUPLE_IDX_ORDER:\n");
    for (i = 0; i < msa->length; i++)
      fprintf(F, "%d\n", ss->tuple_idx[i]);
  }
}

/* make reading order optional?  alphabet argument overrides alphabet
   in file (use NULL to use version in file) */
MSA* ss_read(FILE *F, char *alphabet) {
  Regex *nseqs_re, *length_re, *tuple_size_re, *ntuples_re, *tuple_re, 
    *names_re, *alph_re, *ncats_re, *order_re, *offset_re;
  String *line, *alph = NULL;
  int nseqs, length, tuple_size, ntuples, i, ncats = -99, header_done = 0, 
    idx_offset = 0;
  MSA *msa = NULL;
  List *matches;
  char **names = NULL;

  nseqs_re = str_re_new("NSEQS[[:space:]]*=[[:space:]]*([0-9]+)");
  length_re = str_re_new("LENGTH[[:space:]]*=[[:space:]]*([0-9]+)");
  tuple_size_re = str_re_new("TUPLE_SIZE[[:space:]]*=[[:space:]]*([0-9]+)");
  ntuples_re = str_re_new("NTUPLES[[:space:]]*=[[:space:]]*([0-9]+)");
  names_re = str_re_new("NAMES[[:space:]]*=[[:space:]]*(.*)");
  alph_re = str_re_new("ALPHABET[[:space:]]*=[[:space:]]*([A-Z ]+)");
  ncats_re = str_re_new("NCATS[[:space:]]*=[[:space:]]*([-0-9]+)");
  offset_re = str_re_new("IDX_OFFSET[[:space:]]*=[[:space:]]*([-0-9]+)");
  tuple_re = str_re_new("^([0-9]+)[[:space:]]+([-A-Z ]+)[[:space:]]+([0-9.[:space:]]+)");
  order_re = str_re_new("TUPLE_IDX_ORDER:");
  line = str_new(STR_MED_LEN);
  matches = lst_new_ptr(3);
  nseqs = length = tuple_size = ntuples = -1;

  while (str_readline(line, F) != EOF) {
    str_trim(line);
    if (line->length == 0) continue;

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
          names[i] = strdup(s->chars);
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
        for (i = 0; i < ntuples; i++)
          msa->ss->col_tuples[i] = (char*)smalloc(nseqs * tuple_size *
                                                 sizeof(char));

        str_free(alph);
        header_done = 1;
      }
    }
    
    else if (str_split(line, NULL, matches) >= 3) {
      int idx, offset;
      String *tmpstr;

      tmpstr = lst_get_ptr(matches, 0);
      str_as_int(tmpstr, &idx);
      if (idx < 0 || idx >= ntuples) {
        fprintf(stderr, "ERROR: tuple line has index out of bounds.\nOffending line is, \"%s\"\n", line->chars);
        exit(1);
      }

      for (offset = -1 * (msa->ss->tuple_size-1); offset <= 0; offset++) {
        tmpstr = lst_get_ptr(matches, offset + msa->ss->tuple_size);
        if (tmpstr->length != msa->nseqs) {
          fprintf(stderr, "ERROR: length of column tuple does not match NSEQS.\nOffending line is, \"%s\"\n", line->chars);
          exit(1);
        }          
        for (i = 0; i < msa->nseqs; i++) {
          set_col_char_in_string(msa, msa->ss->col_tuples[idx], i, 
                                 msa->ss->tuple_size, offset, 
                                 tmpstr->chars[i]);
        }
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
        if (str_as_int(line, &msa->ss->tuple_idx[i]) != 0) {
          fprintf(stderr, "ERROR: bad integer in TUPLE_IDX_ORDER list.\n");
          exit(1);
        }        
      }
      if (i < msa->length) {
        fprintf(stderr, "ERROR: too few numbers in TUPLE_IDX_ORDER list.\n");
        exit(1);
      }
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
  str_re_free(tuple_re);
  return msa;
}

/* free all memory associated with a sufficient stats object */
void ss_free(MSA_SS *ss) {
  int j;
    for (j = 0; j < ss->alloc_ntuples; j++)
      free(ss->col_tuples[j]);
    free(ss->col_tuples);
  if (ss->cat_counts != NULL) {
    for (j = 0; j <= ss->msa->ncats; j++) 
      free(ss->cat_counts[j]); 
    free(ss->cat_counts); 
  } 
  if (ss->counts != NULL) free(ss->counts);
  if (ss->tuple_idx != NULL) free(ss->tuple_idx);
  free(ss);
}

/* update category counts, according to 'categories' attribute of MSA
   object.  Requires ordered sufficient stats.  Will allocate and
   initialize cat_counts if necessary. */
void ss_update_categories(MSA *msa) {
  MSA_SS *ss = msa->ss;
  int i, j;
  assert(msa->ncats >= 0 && msa->categories != NULL && 
         ss->tuple_idx != NULL);
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
    assert(msa->categories[i] <= msa->ncats);
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

  assert(col_offset >= 0 && new_tuple_size > 0);

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

  ss_from_msas(new_msa, new_tuple_size, store_order, NULL, NULL, NULL, -1);
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
      assert(msa->ss->tuple_idx[i+start_col] >= 0 && 
             msa->ss->tuple_idx[i+start_col] < msa->ss->ntuples);
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
    if (full_to_sub[tupidx] == -1) continue;

    ss->col_tuples[sub_tupidx] = smalloc(retval->nseqs * ss->tuple_size * 
                                         sizeof(char));
    for (offset = -(ss->tuple_size-1); offset <= 0; offset++) {
      for (i = 0; i < lst_size(include_list); i++) {
        seqidx = lst_get_int(include_list, i);
        ss->col_tuples[sub_tupidx][retval->nseqs*(ss->tuple_size-1 + offset) + i] = 
          msa->ss->col_tuples[tupidx][msa->nseqs*(msa->ss->tuple_size-1 + offset) + seqidx];
      }
    }
    full_to_sub[tupidx] = sub_tupidx++;
  }

  assert(sub_tupidx == sub_ntuples);

  /* copy ordering info for specified columns and recompute counts.
     Also take care of category labels and category-specific counts */

  if (unordered_seqs) {         /* in this case, just copy counts
                                   directly tuple by tuple */
    for (i = 0; i < msa->ss->ntuples; i++) {
      ss->counts[full_to_sub[i]] = msa->ss->counts[i];
      if (do_cats) {
        for (cat = 0; cat <= msa->ncats; cat++)
          ss->cat_counts[cat][full_to_sub[i]] = msa->ss->cat_counts[cat][i];
      }
    }
  }
  else {                        /* go site by site */
    for (i = 0; i < retval->length; i++) {
      ss->tuple_idx[i] = full_to_sub[msa->ss->tuple_idx[i+start_col]];
      assert(ss->tuple_idx[i] >= 0);
      ss->counts[ss->tuple_idx[i]]++;
      if (do_cats) {
        retval->categories[i] = msa->categories[i+start_col];
        ss->cat_counts[retval->categories[i]][ss->tuple_idx[i]]++;
      }
    }  
  }

  if (lst_size(include_list) != msa->nseqs) 
    ss_unique(retval);

  free(full_to_sub);
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

  assert(msa->ss != NULL && msa->ss->tuple_idx != NULL);
  ss = msa->ss;

  if (msa->categories == NULL && ss->cat_counts != NULL)
    fprintf(stderr, "WARNING: ss_reverse_compl cannot address category-specific counts without a\ncategories vector.  Ignoring category counts.  They will be wrong!\n");

  do_cats = (msa->categories != NULL && ss->cat_counts != NULL);

  /* adjust counts for first few columns; these can't be reverse
     complemented */
  for (i = 0; i < ss->tuple_size - 1; i++) {
    ss->counts[ss->tuple_idx[i]]--;
    if (do_cats) ss->cat_counts[msa->categories[i]][ss->tuple_idx[i]]--;
    if (ss->counts[ss->tuple_idx[i]] == 0) 
      que_push_int(overwrites, ss->tuple_idx[i]); 
    /* it's likely that these leading tuples will not be present
       elsewhere in the alignment, so we'll make them available for
       overwriting (see below) */
  }

  /* reverse complement column tuples; counts remain unaltered */
  midpt = ceil(ss->tuple_size/2.0);
  for (i = 0; i < ss->ntuples; i++) {
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
  midpt = ceil(msa->length/2.0);
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
      ss->col_tuples[ss->ntuples] = smalloc(ss->tuple_size * msa->nseqs * 
                                           sizeof(char));
      new_tuple_idx = ss->ntuples++;
    }

    strncpy(ss->col_tuples[new_tuple_idx], new_tuple, 
            msa->nseqs * ss->tuple_size);
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
    strncpy(tmp, msa->ss->col_tuples[tup], msa->nseqs * ts);
    if (new_nseqs > msa->nseqs) 
      msa->ss->col_tuples[tup] = srealloc(msa->ss->col_tuples[tup], 
                                          new_nseqs * ts);
    for (col_offset = -ts+1; col_offset <= 0; col_offset++) {
      for (j = 0; j < new_nseqs; j++) {
        if (new_to_old[j] >= 0)
          msa->ss->col_tuples[tup][new_nseqs * (ts - 1 + col_offset) + j] = 
            tmp[msa->nseqs * (ts - 1 + col_offset) + new_to_old[j]];
        else
          msa->ss->col_tuples[tup][new_nseqs * (ts-1 + col_offset) + j] = 
           msa->missing[0];
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
    else { free(msa->ss->col_tuples[i]); msa->ss->col_tuples[i] = NULL; }
  }

  if (msa->ss->tuple_idx != NULL)
    for (i = 0; i < msa->length; i++)
      msa->ss->tuple_idx[i] = old_to_new[msa->ss->tuple_idx[i]];

  msa->ss->ntuples = new_ntuples;
  ss_compact(msa->ss);  
  free(old_to_new);
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
    strncpy(key, msa->ss->col_tuples[i], msa->nseqs * msa->ss->tuple_size);
    if ((idx = (int)hsh_get(hash, key)) == -1) { /* tuple not seen before */
      hsh_put(hash, key, (void*)i);
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
  free(old_to_new);
  hsh_free(hash);
}

/** Convert all missing data characters to the default missing data
    character (msa->missing[0]).  Optionally also convert gap
    characters.  Can be useful in reducing number of tuples */
void ss_collapse_missing(MSA *msa, int do_gaps) {
  int i, j, len = msa->nseqs * msa->ss->tuple_size;
  int changed_missing = FALSE, changed_gaps = FALSE, exists_missing = FALSE;
  for (i = 0; i < msa->ss->ntuples; i++) {
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
    for (i = 0, j = 0; i < msa->length; i++)
      if (msa->ss->counts[msa->ss->tuple_idx[i]] > 0) {
        msa->ss->tuple_idx[j] = msa->ss->tuple_idx[i];
        if (msa->categories != NULL)
          msa->categories[j] = msa->categories[i];
        j++;
      }
    /* above assumes any tuple is to be removed that is present in
       tuple_idx but has a count of zero */

    assert(newlen == j);
    msa->length = newlen;
  }

  ss_remove_zero_counts(msa);
}

/** Strip columns consisting of missing data in all
   sequences but the reference sequence.  If msa->categories is
   non-NULL, it will be adjusted accordingly */
/* FIXME: what about category counts?  what happens in ss_strip_gaps? */
void ss_strip_missing(MSA *msa, /**< Input alignment; will be altered */
                      int refseq /**< Index of reference sequence
                                    (indexing starts with 1) */
                      ) { 
  int i, j;
  int newlen = msa->length;

  for (i = 0; i < msa->ss->ntuples; i++) {
    int strip = TRUE;
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
    for (i = 0, j = 0; i < msa->length; i++)
      if (msa->ss->counts[msa->ss->tuple_idx[i]] > 0) {
        msa->ss->tuple_idx[j] = msa->ss->tuple_idx[i];
        if (msa->categories != NULL)
          msa->categories[j] = msa->categories[i];
        j++;
      }
    /* above assumes any tuple is to be removed that is present in
       tuple_idx but has a count of zero */

    assert(newlen == j);
    msa->length = newlen;
  }

  ss_remove_zero_counts(msa);
}

/* return 1 if a column triple represents a 4d site and 0 otherwise */
int ss_is_4d(MSA *msa, int tuple) {
  char base[2];
  int offset, seq;
  assert(msa->ss->tuple_size == 3);
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
void ss_reduce_tuple_size(MSA *msa, int new_tuple_size) {
  int i, j, len, offset;
  if (new_tuple_size >= msa->ss->tuple_size)
    die("ERROR: new tuple size must be smaller than old in ss_reduce_tuple_size.\n");
  len = new_tuple_size * msa->nseqs;
  offset = (msa->ss->tuple_size - new_tuple_size) * msa->nseqs;
  for (i = 0; i < msa->ss->ntuples; i++) 
    for (j = 0; j < len; j++) 
      msa->ss->col_tuples[i][j] = msa->ss->col_tuples[i][j+offset];
                                /* because of the way the tuples are
                                   stored, we can just shift everthing
                                   back offset bytes; could realloc
                                   each col_tuple but probably not
                                   worth the cost */
  msa->ss->tuple_size = new_tuple_size;
  ss_unique(msa);
}
