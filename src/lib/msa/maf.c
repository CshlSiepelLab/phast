/* $Id: maf.c,v 1.18 2005-08-30 16:32:26 acs Exp $
   Written by Adam Siepel, 2003
   Copyright 2003, Adam Siepel, University of California */

/** \file maf.c
    Reading of alignments from MAF ("Multiple Alignment Format")
    files, as produced by MULTIZ and TBA.  See maf.h for details.
    \ingroup msa
*/

#include <sufficient_stats.h>
#include <msa.h>
#include <maf.h>
#include <ctype.h>

/** Read an alignment from a MAF file.  The alignment won't be
   constructed explicitly; instead, a sufficient-statistics
   representation will be extracted directly from the MAF.  Blocks
   corresponding to overlapping segments of the reference sequence are
   permitted, but all except the first one will be discarded.  */
MSA *maf_read(FILE *F,          /**< MAF file */
              FILE *REFSEQF,    /**< optional reference sequence.  If
                                   non-NULL, the indicated file will
                                   be used to define bases in regions
                                   of no alignment (not represented in
                                   the MAF).  File format is expected
                                   to be FASTA.  Ignored if
                                   store_order == FALSE.  If NULL and
                                   store_order == TRUE, then bases in
                                   reference seq not present in MAF
                                   are represented as Ns  */
              int tuple_size,   /**< tuple size for sufficient
                                   statistics */
              GFF_Set *gff,     /**< optional GFF_Set.  If non-NULL,
                                   category-specific counts will be
                                   collected (cm must be non-NULL
                                   also).  The gff is assumed to use
                                   the indexing system of the
                                   reference sequence (sequence 1).
                                   Currently, a non-NULL gff implies
                                   gap_strip_mode == 1 (projection
                                   onto reference sequence).  */
              CategoryMap *cm,  /**< Used for category-specific
                                   counts, ignored otherwise */
              int cycle_size,   /**< Label site categories
                                   12...<cycle_size>...12...<cycle_size>
                                   instead of using gff and cm.
                                   Useful when stats are to be
                                   collected for non-overlapping
                                   tuples.  Use -1 to ignore. */
              int store_order,  /**< Whether to store order in which
                                   tuples occur.  Storing order
                                   requires more memory and larger
                                   files, and is not necessary in many
                                   cases. */
              char *reverse_groups, 
                                /**< Tag defining groups in gff;
                                   indicates groups on negative strand
                                   should be reverse complemented.
                                   Ignored if NULL.  Useful when
                                   collecting counts for
                                   strand-specific categories.  Can't
                                   be used if store_order == TRUE */
              int gap_strip_mode,
                                /**< Gap stripping mode.  Currently,
                                   if store_order == 1, may only have
                                   value NO_STRIP or 1 (indicating
                                   projection onto sequence 1, the
                                   reference sequence).  This is
                                   simply to avoid some complexity in
                                   coordinate mapping. */
              int keep_overlapping
                                /**< If TRUE, keep overlapping blocks,
                                   otherwise keep only first instance.
                                   Must be FALSE if store_order ==
                                   TRUE or gff != NULL or 
                                   cycle_size != -1 */
              ) {

/* NOTE: for now, if a GFF is defined, then all blocks are projected
   onto the reference seq (gff != NULL -> gap_strip_mode == 1).  This
   simplifies things somewhat, and it's rare that you want
   category-specific counts without projecting (because gaps in the
   reference sequence make it difficult to assign sites to categories
   rationally).  */

  int i, start_idx, length, max_tuples, block_no, rbl_idx, refseqlen = -1;
  Hashtable *tuple_hash;
  Hashtable *name_hash = hsh_new(25);
  MSA *msa, *mini_msa;
  GFF_Set *mini_gff = NULL;
  int gff_idx = 0;
  List *redundant_blocks = lst_new_int(100);
  msa_coord_map *map = NULL;

  if (gff != NULL) gap_strip_mode = 1; /* for now, automatically
                                          project if GFF (see comment
                                          above) */

  if ((gff == NULL && cm != NULL) || (gff != NULL && cm == NULL)) 
    die("ERROR: maf_read should be passed either both a set of features and a category map, or neither one.\n");

  if (gff != NULL && cycle_size > 0)
    die("ERROR: gff and cycle_size mutually exclusive in maf_read.\n");

  if (store_order) {
    if (reverse_groups != NULL) 
      die("ERROR: Can't reverse complement if storing order in maf_read.\n");
    if (gap_strip_mode != NO_STRIP && gap_strip_mode != 1)
      die("ERROR: Gap strip mode must be either NO_STRIP or 1 if storing order in maf_read.\n");
  }

  if (keep_overlapping && (gff != NULL || store_order || cycle_size > 0))
    die("ERROR: Can't keep overlapping blocks if storing order or collecting stats for site categories (maf_read).\n");

  /* a coordinate map is necessary only if storing order AND not
     projecting on the reference sequence */
  if (store_order && gap_strip_mode == NO_STRIP)
    map = smalloc(sizeof(msa_coord_map));
                                /* inner lists will be allocated by maf_peek */

  /* scan MAF file for total number of sequences and their names, and
     initialize msa accordingly.  Simultaneously build coordinate map,
     if necessary */
  msa = msa_new(NULL, NULL, -1, 0, NULL);
  maf_peek(F, &msa->names, name_hash, &msa->nseqs, map, redundant_blocks, 
           keep_overlapping, &refseqlen);
  /* NOTE: it seems as if this could be avoided when store_order == 0,
     but things would become quite a lot more complicated; e.g., if a
     new seq was encountered midway in the file, all previously
     encountered tuples would have to be redefined  */

  /* init MSA object to be used for individual blocks */
  mini_msa = msa_new(NULL, msa->names, msa->nseqs, -1, NULL);
                                /* note that names are shared */

  if (cm != NULL) msa->ncats = mini_msa->ncats = cm->ncats;
  else if (cycle_size > 0) msa->ncats = mini_msa->ncats = cycle_size;
  else msa->ncats = mini_msa->ncats = -1;

  mini_msa->seqs = smalloc(mini_msa->nseqs * sizeof(char*));
  for (i = 0; i < mini_msa->nseqs; i++) mini_msa->seqs[i] = NULL;
                                /* set up array of seqs; actual alloc
                                   of char strs will occur as needed in
                                   maf_read_block */

  if (gff != NULL) {            /* set up for category labeling */
    gff_sort(gff);
    mini_gff = gff_new_set();
  }

  if (store_order) {
    msa->length = map != NULL ? map->msa_len : refseqlen;
    max_tuples = min(msa->length,
                     pow(strlen(msa->alphabet)+strlen(msa->missing)+1, msa->nseqs * tuple_size));
    if (max_tuples > 1000000) max_tuples = 1000000; 
  }
  else {
    msa->length = 0;
    max_tuples = min(50000,
                     pow(strlen(msa->alphabet)+strlen(msa->missing)+1, msa->nseqs * tuple_size));
  }

  tuple_hash = hsh_new(max_tuples); 
  ss_new(msa, tuple_size, max_tuples, gff != NULL || cycle_size > 0 ? 1 : 0, 
         store_order); 

  if (store_order)
    for (i = 0; i < msa->length; i++) msa->ss->tuple_idx[i] = -1;

  /* process MAF one block at a time */
  block_no = 0;
  rbl_idx = 0;
  while (maf_read_block(F, mini_msa, name_hash, &start_idx, 
                        &length) != EOF) {
    int idx_offset;

    /* ignore if block is marked as redundant */
    if (lst_size(redundant_blocks) > rbl_idx &&
        lst_get_int(redundant_blocks, rbl_idx) == block_no++) {
      rbl_idx++;
      continue;
    }

    if (gap_strip_mode != NO_STRIP) 
      msa_strip_gaps(mini_msa, gap_strip_mode);

    if (gff != NULL) {
      /* extract subset of features in GFF corresponding to block */
      lst_clear(mini_gff->features);
      maf_block_sub_gff(mini_gff, gff, start_idx + 1, start_idx + length - 1, 
                        &gff_idx, cm, reverse_groups != NULL, tuple_size); 
                                /* coords in GFF are 1-based */

      /* if we're not using a global coordinate map, we need to map the
         mini_gff to the coords of the mini_msa */
      /* NOTE: not necessary because automatically projecting */
/*       if (map == NULL && lst_size(mini_gff->features) > 0)  */
/*         msa_map_gff_coords(mini_msa, mini_gff, 1, 0, 0, cm); */

      if (reverse_groups != NULL && lst_size(mini_gff->features) > 0) {
          gff_group(mini_gff, reverse_groups);
          msa_reverse_compl_feats(mini_msa, mini_gff, NULL);
      }

      /* now label categories of mini_msa accordingly */
      msa_label_categories(mini_msa, mini_gff, cm);   
    }
    else if (cycle_size > 0)
      for (i = 0; i < mini_msa->length; i++)
        mini_msa->categories[i] = (i % cycle_size) + 1;

    /* fold new block into aggregate representation */
    /* first map starting coordinate */
    if (map != NULL) {
      idx_offset = msa_map_seq_to_msa(map, start_idx + 1) - 1;

      /* when the reference sequence begins with gaps, 
         start_idx will actually map to the first *non-gap*
         character; we have to adjust accordingly */
      for (i = 0; mini_msa->seqs[0][i] == GAP_CHAR; i++) idx_offset--;

      assert(idx_offset >= 0);
    }

    else if (store_order) idx_offset = start_idx; 

    else idx_offset = -1;           /* no offset */

    /* extract the suff stats from the mini alignment and fold them
       into the new msa */
    ss_from_msas(msa, tuple_size, store_order, NULL, mini_msa, 
                 tuple_hash, idx_offset);

    if (gff != NULL) {          /* free features and clear list */
      for (i = 0; i < lst_size(mini_gff->features); i++)
        gff_free_feature(lst_get_ptr(mini_gff->features, i));
      lst_clear(mini_gff->features);
    }
  }

  /* if necessary, read reference sequence, make sure consistent with
     alignments, fill in remaining tuples */
  if (store_order) {
    char tuple_str[msa->nseqs * tuple_size + 1];
    int offset, tuple_idx, msa_idx, map_idx, fasthash_idx;
    String *refseq;
    int alph_size = strlen(msa->alphabet), nreftuples = int_pow(alph_size, tuple_size);
    int *fasthash = smalloc(nreftuples * sizeof(int));
    char reftuple[tuple_size + 1];

    for (i = 0; i < nreftuples; i++) fasthash[i] = -1;

    if (REFSEQF != NULL)
      refseq = msa_read_seq_fasta(REFSEQF);
    else {
      /* in this case, create a dummy sequence (see below) */
      refseq = str_new(map == NULL ? refseqlen : map->seq_len);
      refseq->length = map == NULL ? refseqlen : map->seq_len;
    }

    if ((map == NULL && refseq->length != refseqlen) ||
        (map != NULL && refseq->length != map->seq_len)) 
      die("ERROR: reference sequence length (%d) does not match description in MAF file (%d).\n", 
          refseq->length, map == NULL ? refseqlen : map->seq_len);

    for (offset = -1 * (tuple_size-1); offset <= 0; offset++) 
      for (i = 1; i < msa->nseqs; i++)
        tuple_str[msa->nseqs*(tuple_size-1 + offset) + i] = msa->missing[0];
    tuple_str[msa->nseqs * tuple_size] = '\0';

    /* look at each site in the reference sequence */
    map_idx = 0;
    reftuple[tuple_size] = '\0';
    for (i = 0, msa_idx = 0; i < refseq->length; i++, msa_idx++) {

      /* use the coord map but avoid a separate lookup at each position */
      if (map != NULL) {
        if (lst_get_int(map->seq_list, map_idx) - 1 == i) 
          msa_idx = lst_get_int(map->msa_list, map_idx++) - 1;
      }
      else msa_idx = i;

      assert(msa_idx >= 0);

      /* simple hack to handle the case where order is stored but 
         refseq is not available: use the char from the alignment if
         available or a missing-data char otherwise */
      if (REFSEQF == NULL) 
        refseq->chars[i] = msa->ss->tuple_idx[msa_idx] == -1 ? 
          msa->missing[0] :
          ss_get_char_pos(msa, msa_idx, 0, 0);

      refseq->chars[i] = toupper(refseq->chars[i]);

      if (msa->inv_alphabet[(int)refseq->chars[i]] < 0 &&
          refseq->chars[i] != GAP_CHAR &&
          !msa->is_missing[(int)refseq->chars[i]] &&
          isalpha(refseq->chars[i]))
        refseq->chars[i] = msa->missing[0];
      /* (for now, assume ambiguity character and treat as missing data) */

      if (msa->ss->tuple_idx[msa_idx] == -1) { /* nothing known about this
                                                  position; set tuple_idx from refseq */
        tuple_idx = fasthash_idx = -1;

        /* try shortcut based on fact that most of the time we'll have
           tuple_size characters from the alphabet */
        if (i >= tuple_size - 1) {
          strncpy(reftuple, &refseq->chars[i-tuple_size+1], tuple_size);
          fasthash_idx = tuple_index(reftuple, msa->inv_alphabet, alph_size);
          if (fasthash_idx != -1 && fasthash[fasthash_idx] != -1)
            tuple_idx = fasthash[fasthash_idx];
        }

        if (tuple_idx == -1) {  /* need full hash lookup */
          for (offset = -1 * (tuple_size-1); offset <= 0; offset++) 
            tuple_str[msa->nseqs*(tuple_size-1 + offset)] = 
              i+offset >= 0 ? refseq->chars[i+offset] : msa->missing[0];

          if ((tuple_idx = (int)hsh_get(tuple_hash, tuple_str)) == -1) {
                                /* tuple isn't in hash yet; have to add */
            tuple_idx = msa->ss->ntuples++;
            hsh_put(tuple_hash, tuple_str, (void*)tuple_idx);
            msa->ss->col_tuples[tuple_idx] = smalloc(tuple_size * msa->nseqs * sizeof(char));
            strncpy(msa->ss->col_tuples[tuple_idx], tuple_str, msa->nseqs * tuple_size);
            if (fasthash_idx != -1) fasthash[fasthash_idx] = tuple_idx;
          }
        }

        msa->ss->counts[tuple_idx]++;
        msa->ss->tuple_idx[msa_idx] = tuple_idx;

        /* NOTE: some context will be lost for sites in the reference
           sequence immediately *following* an alignment block
           (non-reference sequence context will be lost, as will gaps
           in the reference sequence).  However, these sites are no
           important in most analyses (they contain no information
           about substitutions), and the error will have no effect on
           coordinate mapping, format conversion, or anything else I
           can think of; it seems best just to leave things as they
           are.  */
      }

      else if (refseq->chars[i] != ss_get_char_pos(msa, msa_idx, 0, 0) &&
               ss_get_char_pos(msa, msa_idx, 0, 0) != GAP_CHAR) {
        /* here a tuple is available for this position but it does not
           match the sequence */

        /* (make an exception if both chars are not in the alphabet,
           to allow for differences in ambiguity characters) */
        if (msa->inv_alphabet[(int)ss_get_char_pos(msa, msa_idx, 0, 0)] == -1 &&
            msa->inv_alphabet[(int)refseq->chars[i]] == -1)
          ;                     /* okay */
        else 
          die("ERROR: character '%c' at position %d of reference sequence does not match character '%c' given in MAF file.\n", refseq->chars[i], i, ss_get_char_pos(msa, msa_idx, 0, 0));
      }
    }
    
    str_free(refseq);
    free(fasthash);
  }

  mini_msa->names = NULL;       /* will prohibit names from being
                                   freed (they are shared) */
  msa_free(mini_msa);
  if (mini_gff != NULL) gff_free_set(mini_gff);

  hsh_free(tuple_hash);
  hsh_free(name_hash);
  lst_free(redundant_blocks);
  if (map != NULL) msa_map_free(map);

  return msa;
}

/** Read a block from an MAF file and store it as a "mini-msa" using
   the provided object.  Allocates memory for sequences if they are
   NULL (as with first block).  Reads to next "a" line or EOF.
   Returns EOF when no more alignments are available.  Sets start
   coord and length of reference sequence if non-NULL
   pointers are provided.  Uses provided hash to map sequence names to
   sequence indices (prefix of name wrt '.' character); sequences not
   present in a block will be represented by missing-data characters. */
int maf_read_block(FILE *F, MSA *mini_msa, Hashtable *name_hash,
                   int *start_idx, int *length) {

  int seqidx, more_blocks = 0, i, j;
  String *this_seq, *linebuffer = str_new(STR_VERY_LONG_LEN);
  List *l = lst_new_ptr(7);
  String *this_name = str_new(STR_SHORT_LEN);
  int mark[mini_msa->nseqs];

  mini_msa->length = -1;
  for (i = 0; i < mini_msa->nseqs; i++) mark[i] = 0;
  while (str_readline(linebuffer, F) != EOF) {
    if (str_starts_with_charstr(linebuffer, "#")) continue;
    else if (str_starts_with_charstr(linebuffer, "a")) {
      if (mini_msa->length == -1) continue;   /* assume first block (?) */
      more_blocks = 1;          /* want to distinguish a new block
                                   from an EOF */
      break;
    }
    str_trim(linebuffer);
    if (linebuffer->length == 0) continue;

    /* if we get here, linebuffer should contain a sequence line */
    str_split(linebuffer, NULL, l);    
    if (lst_size(l) != 7 || !str_equals_charstr(lst_get_ptr(l, 0), "s")) {
      fprintf(stderr, "ERROR: bad sequence line in MAF file --\n\t\"%s\"\n", linebuffer->chars);
      exit(1);
    }
    str_cpy(this_name, lst_get_ptr(l, 1));
    str_root(this_name, '.');
    this_seq = lst_get_ptr(l, 6);

    /* if this is the reference sequence, also grab start_idx and
       length and check strand */
    if (mini_msa->length == -1 && 
        ((start_idx != NULL && str_as_int(lst_get_ptr(l, 2), start_idx) != 0) ||
        (length != NULL && str_as_int(lst_get_ptr(l, 3), length) != 0) ||
        ((String*)lst_get_ptr(l, 4))->chars[0] != '+')) {
      fprintf(stderr, "ERROR: bad integers or strand in MAF (strand must be + for reference sequence) --\n\t\"%s\"\n", linebuffer->chars);
      exit(1);
    }

    /* ensure lengths of all seqs are consistent */
    if (mini_msa->length == -1) mini_msa->length = this_seq->length;
    else if (this_seq->length != mini_msa->length) {
      fprintf(stderr, "ERROR: sequence lengths do not match in MAF block -- \n\tsee line \"%s\"\n", linebuffer->chars);
      exit(1);
    }

    /* enlarge allocated sequence lengths as necessary */
    if (this_seq->length > mini_msa->alloc_len) {
      mini_msa->alloc_len = this_seq->length;
      for (i = 0; i < mini_msa->nseqs; i++)
        mini_msa->seqs[i] = 
          srealloc(mini_msa->seqs[i], (mini_msa->alloc_len+1) * sizeof(char));
      if (mini_msa->ncats >= 0) 
        mini_msa->categories = 
          srealloc(mini_msa->categories, mini_msa->alloc_len * sizeof(int)); 
    }

    /* obtain index of seq */
    seqidx = (int)hsh_get(name_hash, this_name->chars);
    if (seqidx == -1) {
      fprintf(stderr, "ERROR: unexpected sequence name '%s' --\n\tsee line \"%s\"\n", this_name->chars, linebuffer->chars);
      exit(1);
    }
    assert(str_equals_charstr(this_name, mini_msa->names[seqidx]));

    for (i = 0; i < this_seq->length; i++) {
      mini_msa->seqs[seqidx][i] = toupper(this_seq->chars[i]);
      if (mini_msa->seqs[seqidx][i] == '.' && mini_msa->inv_alphabet[(int)'.'] == -1) 
        mini_msa->seqs[seqidx][i] = mini_msa->missing[0];
      if (mini_msa->seqs[seqidx][i] != GAP_CHAR && 
          !mini_msa->is_missing[(int)mini_msa->seqs[seqidx][i]] &&
          mini_msa->inv_alphabet[(int)mini_msa->seqs[seqidx][i]] == -1) {
        if (isalpha(mini_msa->seqs[seqidx][i]))
          mini_msa->seqs[seqidx][i] = 'N';
        else 
          die("ERROR: unrecognized character in sequence in MAF block ('%c')\n",
              mini_msa->seqs[seqidx][i]);
      }
    }
    mini_msa->seqs[seqidx][this_seq->length] = '\0';
    mark[seqidx] = 1;

    for (i = 0; i < lst_size(l); i++) str_free(lst_get_ptr(l, i));
  }

  lst_free(l);
  str_free(linebuffer);

  if (mini_msa->length == -1 && !more_blocks) 
    return EOF;                 /* in this case, an EOF must have been
                                   encountered before any alignment
                                   blocks were found */


  /* pad unmarked seqs with missing-data characters */
  for (i = 0; i < mini_msa->nseqs; i++) {
    if (!mark[i]) {
      for (j = 0; j < mini_msa->length; j++) 
        mini_msa->seqs[i][j] = mini_msa->missing[0];
      mini_msa->seqs[i][mini_msa->length] = '\0';
    }
  }

  return 0;
}

/** Reads a block from a MAF file and returns it as a new MSA object.
    Reads to next "a" line or EOF.  Returns NULL when no more
    alignments are available.  This is a slightly less efficient but
    simpler and more flexible version of maf_read_block.  */
MSA *maf_read_next_msa(FILE *F) {

  int seqidx = 0, i, firstcall = 0;
  String *this_seq, *this_name, *linebuffer = str_new(STR_VERY_LONG_LEN);
  List *l = lst_new_ptr(7);

  MSA *msa = msa_new(NULL, NULL, 0, -1, NULL);
  while (str_readline(linebuffer, F) != EOF) {
    if (str_starts_with_charstr(linebuffer, "#")) {
      firstcall = 1;            /* probably but not guaranteed */
      continue;
    }
    else if (str_starts_with_charstr(linebuffer, "a")) {
      if (firstcall) { firstcall = 0; continue; }
      break;
    }
    str_trim(linebuffer);
    if (linebuffer->length == 0) continue;

    /* if we get here, linebuffer should contain a sequence line */
    str_split(linebuffer, NULL, l);    
    if (lst_size(l) != 7 || !str_equals_charstr(lst_get_ptr(l, 0), "s")) {
      fprintf(stderr, "ERROR: bad sequence line in MAF file --\n\t\"%s\"\n", linebuffer->chars);
      exit(1);
    }
    this_name = lst_get_ptr(l, 1);
    this_seq = lst_get_ptr(l, 6);

    /* ensure lengths of all seqs are consistent */
    if (msa->length == -1) msa->length = this_seq->length;
    else if (this_seq->length != msa->length) {
      fprintf(stderr, "ERROR: sequence lengths do not match in MAF block -- \n\tsee line \"%s\"\n", linebuffer->chars);
      exit(1);
    }

    /* add new sequence and name */
    msa->nseqs++;
    msa->names = srealloc(msa->names, msa->nseqs * sizeof(char*));
    msa->seqs = srealloc(msa->seqs, msa->nseqs * sizeof(char*));
    msa->names[seqidx] = strdup(this_name->chars);
    msa->seqs[seqidx] = smalloc((msa->length+1) * sizeof(char));
    for (i = 0; i < this_seq->length; i++) {
      msa->seqs[seqidx][i] = toupper(this_seq->chars[i]);
      if (msa->seqs[seqidx][i] == '.')
        msa->seqs[seqidx][i] = msa->missing[0];
      if (msa->seqs[seqidx][i] != GAP_CHAR && 
          !msa->is_missing[(int)msa->seqs[seqidx][i]] &&
          msa->inv_alphabet[(int)msa->seqs[seqidx][i]] == -1) {
        fprintf(stderr, "ERROR: unrecognized character in sequence in MAF block ('%c')\n",
                msa->seqs[seqidx][i]);
        exit(1);
      }
    }
    msa->seqs[seqidx][this_seq->length] = '\0';
    seqidx++;

    for (i = 0; i < lst_size(l); i++) str_free(lst_get_ptr(l, i));
  }

  lst_free(l);
  str_free(linebuffer);

  if (msa->length == -1) {
    msa_free(msa);
    return NULL; 
  }

  return msa;
}


/* these are used in the function below */
struct gap_pair {
  int idx;
  int len;
};

int gap_pair_compare(const void* ptr1, const void* ptr2) {
  struct gap_pair *gp1 = *(struct gap_pair**)ptr1;
  struct gap_pair *gp2 = *(struct gap_pair**)ptr2;
  return gp1->idx - gp2->idx;
}

/* Scan an MAF file for the complete list of sequence names that
   appear (only roots of names are considered).  In addition, fill a
   hashtable that maps each name to the corresponding sequence index
   (must be preallocated), and add the indices of blocks that overlap
   previous blocks to the 'redundant_blocks' list.  If map is
   non-NULL, then construct a coordinate map for the reference
   sequence (map object assumed to be preallocated).  */
void maf_peek(FILE *F, char ***names, Hashtable *name_hash, 
              int *nseqs, msa_coord_map *map, List *redundant_blocks,
              int keep_overlapping, int *refseqlen) {
  String *line = str_new(STR_VERY_LONG_LEN);
  int count = 0, seqidx = 0, gaplen, tmp, startidx, i, 
    length, last_endidx = -1, block_no = 0, skip = 0, endidx;
  String *fullname = str_new(STR_SHORT_LEN), *name = str_new(STR_SHORT_LEN);
  String *s;
  List *gp_list = (map != NULL ? lst_new_ptr(10000) : NULL);
  List *l = lst_new_ptr(7), *block_starts = lst_new_int(1000), 
    *block_ends = lst_new_int(1000);
  fpos_t pos;

  *refseqlen = -1;

  if (fgetpos(F, &pos) != 0)
    die("ERROR: Currently, MAF input stream must be seekable (can't be stdin).\n");

  while (str_readline(line, F) != EOF) {
    if (line->chars[0] == 'a') seqidx = 0;
    else if (line->chars[0] == 's') {
      /* avoid calling str_split for efficiency */
      for (i = 1; i < line->length && isspace(line->chars[i]); i++);
      str_clear(fullname);
      for (; i < line->length && !isspace(line->chars[i]); i++)
        str_append_char(fullname, line->chars[i]);
      str_cpy(name, fullname);
      str_root(name, '.');
      assert(name->length > 0); /* must be a non-empty name */
      if ((int)hsh_get(name_hash, name->chars) == -1) {
        hsh_put(name_hash, name->chars, (void*)count);
        *names = srealloc(*names, (count+1) * sizeof(char*));
        (*names)[count] = strdup(name->chars);
        count++;
      }

      if (seqidx == 0 && !keep_overlapping) { /* reference sequence */
        str_split(line, NULL, l);
        if (lst_size(l) != 7 || 
            str_as_int(lst_get_ptr(l, 2), &startidx) != 0 ||
            str_as_int(lst_get_ptr(l, 3), &length) != 0 ||
            str_as_int(lst_get_ptr(l, 5), &tmp) != 0) 
          die("ERROR: bad line in MAF file --\n\t\"%s\"\n", line->chars);

        if (*refseqlen == -1) *refseqlen = tmp;

        skip = 0;
        if (tmp != *refseqlen || 
            ((String*)lst_get_ptr(l, 4))->chars[0] != '+') {
          fprintf(stderr, "WARNING: Reference sequence length or strand does not match previous blocks, ignoring block starting with:\n\t%s", line->chars);
          lst_push_int(redundant_blocks, block_no);
          skip = 1;
        }

        /* Check whether new block overlaps a previously seen block.
           Most of the time, it will be downstream of the previous
           block; check for this first */
        endidx = startidx + length - 1;
        if (!skip && startidx > last_endidx) {
          lst_push_int(block_starts, startidx);
          lst_push_int(block_ends, endidx);
          last_endidx = endidx;
        }
        else if (!skip) {       /* have to search list */
          int block_list_idx = lst_bsearch_int(block_starts, startidx);
          int prev_end = block_list_idx >= 0 ? 
            lst_get_int(block_ends, block_list_idx) : -1;
          int next_start = block_list_idx+1 < lst_size(block_starts) ?
            lst_get_int(block_starts, block_list_idx+1) : endidx+1;         
          if (prev_end >= startidx || next_start <= endidx) {
/*             fprintf(stderr, "WARNING: MAF block (%d-%d in ref. seq.) overlaps a previous block -- ignoring.\n", startidx, endidx); */
            lst_push_int(redundant_blocks, block_no);
            skip = 1;
          }
          else {
            lst_insert_idx_int(block_starts, block_list_idx, startidx);
            lst_insert_idx_int(block_ends, block_list_idx, endidx);
            if (endidx > last_endidx) last_endidx = endidx;
          }
        }

        /* collect info on gaps for coordinate map */
        if (map != NULL && !skip) {
          s = lst_get_ptr(l, 6);
          gaplen = 0;
          for (i = 0; i < s->length; i++) {
            if (s->chars[i] == GAP_CHAR) gaplen++;
            else {
              if (gaplen > 0) {
                struct gap_pair *gp = smalloc(sizeof(struct gap_pair));
                gp->idx = startidx;
                gp->len = gaplen;
                lst_push_ptr(gp_list, gp);

                gaplen = 0;
              }
              startidx++;
            }
          }
          /* there may be a gap at the end of the block */
          if (gaplen > 0) {
            struct gap_pair *gp = smalloc(sizeof(struct gap_pair));
            gp->idx = startidx;
            gp->len = gaplen;
            lst_push_ptr(gp_list, gp);
          }
        }
        for (i = 0; i < lst_size(l); i++) str_free(lst_get_ptr(l, i));
        block_no++;
      }
      seqidx++;
    }
  }
  fsetpos(F, &pos);
  str_free(line); str_free(fullname); str_free(name);
  lst_free(l);
  *nseqs = count;

  /* now build coordinate map, if necessary */
  if (map != NULL) {
    int partial_gap_sum = 0;
    struct gap_pair *gp, *nextgp;
    
    lst_qsort(gp_list, gap_pair_compare);    

    map->seq_list = lst_new_int(lst_size(gp_list) + 1);
    map->msa_list = lst_new_int(lst_size(gp_list) + 1);

    /* "prime" coord map */
    lst_push_int(map->seq_list, 1); 
    lst_push_int(map->msa_list, 1);

    /* build coord map from gap list */
    for (i = 0; i < lst_size(gp_list); i++) {
      gp = lst_get_ptr(gp_list, i);
      nextgp = (i != lst_size(gp_list) - 1 ? lst_get_ptr(gp_list, i+1) : NULL);

      partial_gap_sum += gp->len;

      /* if there is a gap prior to the beginning of the reference seq,
         then the first element of map->msa_list has to be reset */
      if (i == 0 && gp->idx == 0) {
        lst_set_int(map->msa_list, 0, partial_gap_sum + 1);
        continue;
      }

      /* if gaps occur at the end of one block and at the beginning of an
         immediate successor, then they have to be merged */
      if (nextgp != NULL && nextgp->idx == gp->idx) continue;

      lst_push_int(map->seq_list, gp->idx + 1);
      lst_push_int(map->msa_list, gp->idx + partial_gap_sum + 1);
                                /* note: coord map uses 1-based indexing */
      free(gp);
    }
    map->seq_len = *refseqlen;
    map->msa_len = *refseqlen + partial_gap_sum;
    lst_free(gp_list); lst_free(block_starts); lst_free(block_ends);
  }
}

/* Extracts features from gff relevant to the interval [start_idx,
   end_idx] and stores them in sub_gff (assumed to be allocated but to
   have an empty feature list).  Truncates overlapping features if
   possible (see details below).  Shifts all coords such that
   start_idx is position 1.  Assumes main gff is sorted.  Designed for
   repeated calls. */
void maf_block_sub_gff(GFF_Set *sub_gff, GFF_Set *gff, int start_idx, 
                       int end_idx, int *gff_idx, CategoryMap *cm,
                       int reverse_compl, int tuple_size) {
  static int last_nonoverlap = -1;
  GFF_Feature *feat;
  for (; *gff_idx < lst_size(gff->features) && 
         ((GFF_Feature *)lst_get_ptr(gff->features, *gff_idx))->end < 
         start_idx;
       (*gff_idx)++);            /* ignore upstream features */
  for (; *gff_idx < lst_size(gff->features) && 
         (feat = lst_get_ptr(gff->features, *gff_idx))->start <= end_idx;
       (*gff_idx)++) {           /* look at all that overlap */
    GFF_Feature *featcpy;
    if (feat->end < start_idx) continue; /* possible because sorted by
                                            start index */

    /* address overlapping features */

    /* truncate the feature and keep it if it is does not correspond
       to a "range category" or if it has frame information, but
       otherwise throw it out, because truncating it will likely
       result in an incorrect category labeling.  This strategy is
       possible because frame is sufficient for correct labeling of
       range cats */
    if (feat->start < start_idx || feat->end > end_idx) {
      int cat = cm_get_category(cm, feat->feature);
      if (cm->ranges[cat]->start_cat_no != cm->ranges[cat]->end_cat_no &&
          feat->frame == GFF_NULL_FRAME)
        continue;
    }

    featcpy = gff_new_feature_copy(feat);
    if (featcpy->start < start_idx) {
      if (featcpy->strand == '+' && featcpy->frame != GFF_NULL_FRAME) 
        featcpy->frame = (featcpy->frame + start_idx - featcpy->start) % 3;
      featcpy->start = start_idx;
    }
    if (featcpy->end > end_idx) {
      int effective_end = end_idx;
      if (featcpy->strand == '-' && reverse_compl) 
        effective_end -= (tuple_size - 1);
                                /* if we truncate a feature that is to
                                   be reverse complemented, we have to
                                   be careful not to introduce
                                   artificial context */
      if (featcpy->strand == '-' && featcpy->frame != GFF_NULL_FRAME) 
        featcpy->frame = (featcpy->frame + featcpy->end - effective_end) % 3;
      featcpy->end = effective_end;
    }

    if (feat->end < end_idx) last_nonoverlap = *gff_idx;

    /* shift coords and add feature */
    featcpy->start -= (start_idx - 1);
    featcpy->end -= (start_idx - 1);

    lst_push_ptr(sub_gff->features, featcpy);
  }

  if (last_nonoverlap != -1) *gff_idx = last_nonoverlap;
                                /* this allows features that overlap a
                                   block to be considered for the next
                                   block */
}
