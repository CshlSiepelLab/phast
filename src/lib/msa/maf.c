 /***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#include <sufficient_stats.h>
#include <msa.h>
#include <maf.h>
#include <ctype.h>
#include <maf_block.h>
#include <misc.h>


/** Read An Alignment from a MAF file.  The alignment won't be
   constructed explicitly; instead, a sufficient-statistics
   representation will be extracted directly from the MAF.  The MAF file must
   be sorted with respect to the reference sequence.  Any blocks falling out
   of order, or which are redundant with previous blocks, will be discarded.
   This allows the MAF to be read in one pass.*/
MSA *maf_read_cats_subset(FILE *F, FILE *REFSEQF, int tuple_size, 
   char *alphabet, GFF_Set *gff, CategoryMap *cm, int cycle_size, 
   int store_order, char *reverse_groups, int gap_strip_mode, 
   int keep_overlapping, List *cats_to_do, List *seqnames, int seq_keep ) {

  int i, start_idx, length, max_tuples, block_no,  
    refseqlen = -1, do_toupper, last_refseqpos = -1;
  Hashtable *tuple_hash;
  Hashtable *name_hash = hsh_new(25);
  MSA *msa, *mini_msa;
  GFF_Set *mini_gff = NULL;
  int gff_idx = 0, refseq_sorted = 1;
  msa_coord_map *map = NULL;
  List *block_starts = lst_new_int(1000), *block_ends = lst_new_int(1000);
  int last_gap_start = -1;
  int idx_offset, end_idx, gap_sum=0;
  int block_list_idx, prev_end, next_start;
  int first_idx=-1, last_idx=-1, free_cm=0;

  if (gff != NULL) gap_strip_mode = 1; /* for now, automatically
                                          project if GFF (see comment
                                          above) */

  if (tuple_size != 1) 
    phast_warning("Warning: reading in MAF with tuple_size > 1 loses information that crosses blocks");

  if (gff ==  NULL && cm != NULL) 
    die("ERROR: maf_read got non-null category map without a set of features");
  if (gff != NULL && cm == NULL) {
    cm = cm_new_from_features(gff);
    free_cm = 1;
  }
  if (gff != NULL && cycle_size > 0)
    die("ERROR: gff and cycle_size mutually exclusive in maf_read.\n");
  
  if (cats_to_do != NULL && (gff == NULL || cm == NULL)) 
    die("ERROR: cats_to_do requires gff and cm != NULL in maf_read.\n");

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
  if (store_order && gap_strip_mode == NO_STRIP) {
    map = smalloc(sizeof(msa_coord_map));
    map->seq_list = lst_new_int(1);
    map->msa_list = lst_new_int(1);
    /* "prime" coord map */
    /* Note: re-prime map->seq_list later if msa->idx_offset > 0 */
    lst_push_int(map->seq_list, 1);
    lst_push_int(map->msa_list, 1);
  }
                                /* inner lists will be allocated by maf_peek */

  msa = msa_new(NULL, NULL, -1, 0, alphabet);

  if (seqnames != NULL) {
    String *currname;
    if (seq_keep) {
      msa->names = smalloc(lst_size(seqnames)*sizeof(char*));
      for (i=0; i<lst_size(seqnames); i++) {
	currname = (String*)lst_get_ptr(seqnames, i);
	hsh_put_int(name_hash, currname->chars, i);
	msa->names[i] = copy_charstr(currname->chars);
      }
      msa->nseqs = lst_size(seqnames);
      maf_quick_peek(F, &msa->names, name_hash, NULL, &refseqlen, 0);
    } else {
      for (i=0; i<lst_size(seqnames); i++) {
	currname = (String*)lst_get_ptr(seqnames, i);
	hsh_put_int(name_hash, currname->chars, -2);
      }
      maf_quick_peek(F, &msa->names, name_hash, &msa->nseqs, &refseqlen, 1);
    }
  } else {
    //look at first block to get initial sequence names and refseqlen
    maf_quick_peek(F, &msa->names, name_hash, &msa->nseqs, &refseqlen, 1);
  }
  if (msa->nseqs == 0 || refseqlen==-1) 
    die("ERROR: got invalid maf file\n");
  if (map != NULL)
    map->seq_len = refseqlen;

  /* upcase chars unless there are lowercase characters in the alphabet */
  do_toupper = !msa_alph_has_lowercase(msa);    

  /* init MSA object to be used for individual blocks */
  mini_msa = msa_new(NULL, msa->names, msa->nseqs, -1, alphabet);
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

  msa->length = 0;
  if (store_order) {
    if (REFSEQF != NULL) 
      msa->alloc_len = msa->length = refseqlen;  //this may still not be big enough because of gaps in refseq
    else 
      msa->alloc_len = 50000;
    max_tuples =  max(1000000, (int)pow(strlen(msa->alphabet)+strlen(msa->missing)+1, 
				   2 * msa->nseqs * tuple_size));
    if (max_tuples > 10000000 || max_tuples < 0) max_tuples = 10000000;
    if (max_tuples < 1000000) max_tuples = 1000000;
  }
  else 
    max_tuples = min(1000000,
		     (int)pow(strlen(msa->alphabet)+strlen(msa->missing)+1, 2 * msa->nseqs * tuple_size));

  tuple_hash = hsh_new(max_tuples); 
  ss_new(msa, tuple_size, max_tuples, gff != NULL || cycle_size > 0 ? 1 : 0, 
         store_order); 

  if (store_order) {
    for (i = 0; i < msa->length; i++) 
      msa->ss->tuple_idx[i] = -1;
  }

  /* process MAF one block at a time */
  block_no = 0;
  while (maf_read_block_addseq(F, mini_msa, name_hash, &start_idx,
			       &length, do_toupper, seqnames != NULL && seq_keep) != EOF) {
    checkInterruptN(block_no++, 1000);

    //sequence may have been added in maf_read_block so reset numseqs
    //do not have to reset msa->names since they are shared with mini_msa
    if (msa->nseqs < mini_msa->nseqs) {
      msa_add_seq_ss(msa, mini_msa->nseqs);
      msa->nseqs = mini_msa->nseqs;
    }
    end_idx = start_idx + length - 1;

    /* if creating a map, require MAF to be sorted wrt reference sequence, otherwise skip block */
    if (map != NULL && start_idx <= last_refseqpos) {
      if (refseq_sorted) {
	phast_warning("warning: maf_read: MAF file must be sorted with respect to reference" \
		      " sequence if store_order=TRUE.  Ignoring out-of-order blocks\n");
	/* only print the warning once */
	refseq_sorted=0;
      }
      continue;
    }
    /* ignore if redundant block: if start_idx < last_refseqpos, need to check list to
       see if region is redundant */
    if (!keep_overlapping && start_idx <= last_refseqpos) {
       block_list_idx = lst_bsearch_int(block_starts, start_idx);
       prev_end = block_list_idx >=0 ? lst_get_int(block_ends, block_list_idx) : -1;
       next_start = block_list_idx + 1 < lst_size(block_starts) ?
	lst_get_int(block_starts, block_list_idx+1) : end_idx + 1;
       if (prev_end >= start_idx || next_start <= end_idx)  { //redundant
	 continue;
       }
    }
    /* also ignore if block size is less than tuple size */
    if (length < tuple_size)  {
      continue; 
    }

    if (first_idx == -1) {
      first_idx = start_idx;
      if (store_order && REFSEQF == NULL) {
        msa->idx_offset = first_idx < 0 ? 0 : first_idx;
        /* reprime map->seq_list if necessary */
        if (map != NULL && first_idx != 0)
          lst_set_int(map->seq_list, 0, msa->idx_offset + 1);
      }
    }
    if (start_idx + length > last_idx)
      last_idx = start_idx + length;

    /* add block to list to check for redundant blocks later */
    lst_push_int(block_starts, start_idx);
    lst_push_int(block_ends, end_idx);

    last_refseqpos = end_idx;

    /* collect info on gaps for coordinate map */
    if (map != NULL) {
      int idx = start_idx, gaplen = 0, gapsum_block=0;
      for (i = 0; i < mini_msa->length; i++) {
	if (mini_msa->seqs[0][i] == GAP_CHAR) gaplen++;
	else {
	  if (gaplen > 0) {
	    gap_sum += gaplen;
	    if (idx == msa->idx_offset) 
	      lst_set_int(map->msa_list, 0, gap_sum + 1);
	    else if (idx == last_gap_start) 
	      lst_set_int(map->msa_list, lst_size(map->msa_list)-1, 
			  idx + gap_sum + 1 - msa->idx_offset);
	    else {
	      lst_push_int(map->msa_list, idx + gap_sum + 1 - msa->idx_offset);
	      lst_push_int(map->seq_list, idx + 1);
	    }
	    last_gap_start = idx;
	  }
	  gapsum_block += gaplen;
	  gaplen = 0;
	  idx++;
	}
      }
      /* there may be a gap at the end of the block */
      if (gaplen > 0) {
	gapsum_block += gaplen;
	gap_sum += gaplen;
	if (idx == last_gap_start) 
	  lst_set_int(map->msa_list, lst_size(map->msa_list)-1, 
		      idx + gap_sum + 1 - msa->idx_offset);
	else {
	  lst_push_int(map->msa_list, idx + gap_sum + 1 - msa->idx_offset);
	  lst_push_int(map->seq_list, idx + 1);
	}
	last_gap_start = idx;
      }
      /*      msa->length += gapsum_block;
      if (msa->length > msa->alloc_len) {
      	msa_realloc(msa, msa->length, max(2 * msa->length, idx + gap_sum), 0,
		    store_order);
	
		    }*/
    } /* end coordinate map section */
    
    if (gap_strip_mode != NO_STRIP) 
      msa_strip_gaps(mini_msa, gap_strip_mode);

    if (gff != NULL) {
      /* extract subset of features in GFF corresponding to block */
      lst_clear(mini_gff->features);
      maf_block_sub_gff(mini_gff, gff, start_idx + 1, start_idx + length, 
                        &gff_idx, cm, reverse_groups != NULL, tuple_size); 
                                /* coords in GFF are 1-based */

      /* if we're not using a global coordinate map, we need to map the
         mini_gff to the coords of the mini_msa */
      /* NOTE: not necessary because automatically projecting */
/*       if (map == NULL && lst_size(mini_gff->features) > 0)  */
/*         msa_map_gff_coords(mini_msa, mini_gff, 1, 0, 0); */

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
      if (idx_offset < 0)
	die("ERROR maf_read_subset: invalid idx_offset %i\n", idx_offset);

      /* when the reference sequence begins with gaps, 
         start_idx will actually map to the first *non-gap*
         character; we have to adjust accordingly */
      for (i = 0; mini_msa->seqs[0][i] == GAP_CHAR; i++) idx_offset--;

      if (idx_offset < 0)
	die("ERROR maf_read_subset: invalid idx_offset2 %i\n", idx_offset);
    }

    else if (store_order) idx_offset = start_idx - msa->idx_offset; 

    else idx_offset = -1;           /* no offset */

    /* extract the suff stats from the mini alignment and fold them
       into the new msa */
    ss_from_msas(msa, tuple_size, store_order, cats_to_do, mini_msa, 
                 tuple_hash, idx_offset, 0);

    if (gff != NULL) {          /* free features and clear list */
      for (i = 0; i < lst_size(mini_gff->features); i++)
        gff_free_feature(lst_get_ptr(mini_gff->features, i));
      lst_clear(mini_gff->features);
    }
  }
  if (map != NULL)
    map->msa_len = map->seq_len + gap_sum;

  /* if necessary, read reference sequence, make sure consistent with
     alignments, fill in remaining tuples */
  if (store_order) {
    char tuple_str[msa->nseqs * tuple_size + 1];
    int offset, tuple_idx, msa_idx, map_idx, fasthash_idx;
    String *refseq;
    int alph_size = (int)strlen(msa->alphabet), nreftuples = int_pow(alph_size, tuple_size);
    int *fasthash = smalloc(nreftuples * sizeof(int));
    char reftuple[tuple_size + 1];

    for (i = 0; i < nreftuples; i++) fasthash[i] = -1;

    if (REFSEQF != NULL) {
      refseq = msa_read_seq_fasta(REFSEQF);
      if (refseq->length != refseqlen) 
	die("ERROR: reference sequence length (%d) does not match description in MAF file (%d).\n", 
	    refseq->length, refseqlen);
    }
    else {
      /* in this case, create a dummy sequence (see below) */
      /* update: if store_order is true but there is no REFSEQF, only store parts of
       the alignment that fall inside the range of coordinates found in MAF file.  So
       we still create a dummy sequence (in order to fill parts of alignment between blocks),
       but it is length last_idx-first_idx.*/
      msa->length = last_idx - first_idx + gap_sum;
      msa_realloc(msa, msa->length, msa->length, 0, store_order);
      refseq = str_new(last_idx - first_idx);
      refseq->length = last_idx - first_idx;
    }

    for (offset = -1 * (tuple_size-1); offset <= 0; offset++) 
      for (i = 1; i < msa->nseqs; i++)
	tuple_str[tuple_size*i + tuple_size -1 + offset] = msa->missing[0];
    tuple_str[msa->nseqs * tuple_size] = '\0';

    /* look at each site in the reference sequence */
    map_idx = 0;
    reftuple[tuple_size] = '\0';


    for (i = 0, msa_idx = 0; i < refseq->length; i++, msa_idx++) {

      /* use the coord map but avoid a separate lookup at each position */
      if (map != NULL) {
        if (lst_get_int(map->seq_list, map_idx) - 1 == i + msa->idx_offset) 
          msa_idx = lst_get_int(map->msa_list, map_idx++) - 1;
      }
      else msa_idx = i;

      if (msa_idx >= msa->length)
	msa_realloc(msa, msa_idx+1, msa_idx + 10000, 0, store_order);
      
      if (msa_idx < 0)
	die("ERROR maf_read_subset: msa_idx=%i, should be >=0\n",
	    msa_idx);

      /* simple hack to handle the case where order is stored but 
         refseq is not available: use the char from the alignment if
         available or a missing-data char otherwise */
      if (REFSEQF == NULL) 
        refseq->chars[i] = msa->ss->tuple_idx[msa_idx] == -1 ? 
          msa->missing[1] :
          ss_get_char_pos(msa, msa_idx, 0, 0);

      if (do_toupper)
        refseq->chars[i] = (char)toupper(refseq->chars[i]);

      if (msa->inv_alphabet[(int)refseq->chars[i]] < 0 &&
          refseq->chars[i] != GAP_CHAR &&
          !msa->is_missing[(int)refseq->chars[i]] &&
	  get_iupac_map()[(int)refseq->chars[i]] == NULL &&
          isalpha(refseq->chars[i]))
        refseq->chars[i] = msa->missing[1];
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
	    tuple_str[tuple_size-1 + offset] =
              i+offset >= 0 ? refseq->chars[i+offset] : msa->missing[1];

	  if ((tuple_idx = ss_lookup_coltuple(tuple_str, tuple_hash, msa)) == -1) {
                                /* tuple isn't in hash yet; have to add */
            tuple_idx = msa->ss->ntuples++;
	    ss_add_coltuple(tuple_str, int_to_ptr(tuple_idx), tuple_hash, msa);
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

      else {
	if (refseq->chars[i] != ss_get_char_pos(msa, msa_idx, 0, 0) &&
	    ss_get_char_pos(msa, msa_idx, 0, 0) != GAP_CHAR) {
	  /* here a tuple is available for this position but it does not
	     match the sequence */
	  
	  /* (make an exception if both chars are not in the alphabet,
           to allow for differences in ambiguity characters) */
        if (msa->inv_alphabet[(int)ss_get_char_pos(msa, msa_idx, 0, 0)] == -1 &&
            msa->inv_alphabet[(int)refseq->chars[i]] == -1)
          ;                     /* okay */
        /* (also make an exception if in the alphabet but missing) */
        else if (msa->is_missing[(int)ss_get_char_pos(msa, msa_idx, 0, 0)] &&
            msa->is_missing[(int)refseq->chars[i]])
          ;                     /* okay */
        /* (also make an exception if one is softmasked and the other isn't) */
        else if (!do_toupper && toupper(refseq->chars[i]) == 
            toupper(ss_get_char_pos(msa, msa_idx, 0, 0)))
          ;
        else 
          die("ERROR: character '%c' at position %d of reference sequence does not match character '%c' given in MAF file.\n", refseq->chars[i], i, ss_get_char_pos(msa, msa_idx, 0, 0));
      }
      }
    }
    
    str_free(refseq);
    sfree(fasthash);
  }

  msa->names = mini_msa->names;
  mini_msa->names = NULL;       /* will prohibit names from being
                                   freed (they are shared) */
  msa_free(mini_msa);
  if (mini_gff != NULL) gff_free_set(mini_gff);

  hsh_free(tuple_hash);
  hsh_free(name_hash);
  lst_free(block_starts);
  lst_free(block_ends);
  if (map != NULL) msa_map_free(map);
  if (free_cm) cm_free(cm);
  return msa;
}


MSA *maf_read_cats(FILE *F, FILE *REFSEQF, int tuple_size,  
	char *alphabet,  GFF_Set *gff,   CategoryMap *cm, int cycle_size,
	int store_order,  char *reverse_groups, int gap_strip_mode,
	int keep_overlapping,List *cats_to_do) {
  return maf_read_cats_subset(F, REFSEQF, tuple_size, alphabet, gff, cm, cycle_size,
			      store_order, reverse_groups, gap_strip_mode, 
			      keep_overlapping, cats_to_do, NULL, 0);
}

MSA *maf_read(FILE *F, FILE *REFSEQF, int tuple_size, char *alphabet,
	      GFF_Set *gff, CategoryMap *cm, int cycle_size, int store_order, 
	      char *reverse_groups, int gap_strip_mode, int keep_overlapping) {
  return maf_read_cats(F, REFSEQF, tuple_size, alphabet, gff, cm, cycle_size, store_order,
		reverse_groups, gap_strip_mode, keep_overlapping, NULL);
}

/* Read An Alignment from a MAF file which is not necessarily sorted wrt the
    reference sequence.  The alignment won't be
   constructed explicitly; instead, a sufficient-statistics
   representation will be extracted directly from the MAF.  Blocks
   corresponding to overlapping segments of the reference sequence are
   permitted, but all except the first one will be discarded.  */  
MSA *maf_read_unsorted(FILE *F, FILE *REFSEQF, int tuple_size, char *alphabet,
	 GFF_Set *gff, CategoryMap *cm, int cycle_size, int store_order,
	 char *reverse_groups, int gap_strip_mode, int keep_overlapping, 
	 List *cats_to_do ) {

/* NOTE: for now, if a GFF is defined, then all blocks are projected
   onto the reference seq (gff != NULL -> gap_strip_mode == 1).  This
   simplifies things somewhat, and it's rare that you want
   category-specific counts without projecting (because gaps in the
   reference sequence make it difficult to assign sites to categories
   rationally).  */

/* Note: this is the original version of maf_read which does NOT require the
   MAF file to be sorted with respect to the reference sequence.  As a result,
   it has to read through the MAF file twice to build a map, which can be 
   quite slow with large MAF files */

  int i, start_idx, length, max_tuples, block_no, rbl_idx, 
    refseqlen = -1, do_toupper;
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
  msa = msa_new(NULL, NULL, -1, 0, alphabet);
  maf_peek(F, &msa->names, name_hash, &msa->nseqs, map, redundant_blocks, 
           keep_overlapping, &refseqlen);
  /* NOTE: it seems as if this could be avoided when store_order == 0,
     but things would become quite a lot more complicated; e.g., if a
     new seq was encountered midway in the file, all previously
     encountered tuples would have to be redefined  */

  /* upcase chars unless there are lowercase characters in the alphabet */
  do_toupper = !msa_alph_has_lowercase(msa);    

  /* init MSA object to be used for individual blocks */
  mini_msa = msa_new(NULL, msa->names, msa->nseqs, -1, alphabet);
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
                     (int)pow(strlen(msa->alphabet)+strlen(msa->missing)+1, msa->nseqs * tuple_size));
    if (max_tuples > 1000000) max_tuples = 1000000; 
  }
  else {
    msa->length = 0;
    max_tuples = min(50000,
                     (int)pow(strlen(msa->alphabet)+strlen(msa->missing)+1, msa->nseqs * tuple_size));
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
                        &length, do_toupper) != EOF) {
    int idx_offset;
    checkInterruptN(block_no, 1000);

    /* ignore if block is marked as redundant */
    if (lst_size(redundant_blocks) > rbl_idx &&
        lst_get_int(redundant_blocks, rbl_idx) == block_no++) {
      rbl_idx++;
      continue;
    }

    /* also ignore if block size is less than tuple size */
    if (length < tuple_size) 
      continue; 

    if (gap_strip_mode != NO_STRIP) 
      msa_strip_gaps(mini_msa, gap_strip_mode);

    if (gff != NULL) {
      /* extract subset of features in GFF corresponding to block */
      lst_clear(mini_gff->features);
      maf_block_sub_gff(mini_gff, gff, start_idx + 1, start_idx + length, 
                        &gff_idx, cm, reverse_groups != NULL, tuple_size); 
                                /* coords in GFF are 1-based */

      /* if we're not using a global coordinate map, we need to map the
         mini_gff to the coords of the mini_msa */
      /* NOTE: not necessary because automatically projecting */
/*       if (map == NULL && lst_size(mini_gff->features) > 0)  */
/*         msa_map_gff_coords(mini_msa, mini_gff, 1, 0, 0); */

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

      if (idx_offset < 0)
	die("ERROR maf_read_unsorted: idx_offset=%i\n", idx_offset);
    }

    else if (store_order) idx_offset = start_idx; 

    else idx_offset = -1;           /* no offset */

    /* extract the suff stats from the mini alignment and fold them
       into the new msa */
    ss_from_msas(msa, tuple_size, store_order, cats_to_do, mini_msa, 
                 tuple_hash, idx_offset, 0);

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
    int alph_size = (int)strlen(msa->alphabet), nreftuples = int_pow(alph_size, tuple_size);
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
	tuple_str[i*tuple_size + tuple_size - 1 + offset] = msa->missing[0];
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

      if (msa_idx < 0)
	die("ERROR maf_read_unsorted: msa_idx=%i\n", msa_idx);

      /* simple hack to handle the case where order is stored but 
         refseq is not available: use the char from the alignment if
         available or a missing-data char otherwise */
      if (REFSEQF == NULL) 
        refseq->chars[i] = msa->ss->tuple_idx[msa_idx] == -1 ? 
          msa->missing[1] :
          ss_get_char_pos(msa, msa_idx, 0, 0);

      if (do_toupper)
        refseq->chars[i] = (char)toupper(refseq->chars[i]);

      if (msa->inv_alphabet[(int)refseq->chars[i]] < 0 &&
          refseq->chars[i] != GAP_CHAR &&
          !msa->is_missing[(int)refseq->chars[i]] &&
	  get_iupac_map()[(int)refseq->chars[i]] == NULL &&
          isalpha(refseq->chars[i]))
        refseq->chars[i] = msa->missing[1];
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
	    tuple_str[tuple_size - 1 + offset] = 
              i+offset >= 0 ? refseq->chars[i+offset] : msa->missing[1];
	  if ((tuple_idx = ss_lookup_coltuple(tuple_str, tuple_hash, msa)) == -1) {
                                /* tuple isn't in hash yet; have to add */
            tuple_idx = msa->ss->ntuples++;
	    ss_add_coltuple(tuple_str, int_to_ptr(tuple_idx), tuple_hash, msa);
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
        /* (also make an exception if in the alphabet but missing) */
        else if (msa->is_missing[(int)ss_get_char_pos(msa, msa_idx, 0, 0)] &&
            msa->is_missing[(int)refseq->chars[i]])
          ;                     /* okay */
        /* (also make an exception if one is softmasked and the other isn't) */
        else if (!do_toupper && toupper(refseq->chars[i]) == 
            toupper(ss_get_char_pos(msa, msa_idx, 0, 0)))
          ;
        else 
          die("ERROR: character '%c' at position %d of reference sequence does not match character '%c' given in MAF file.\n", refseq->chars[i], i, ss_get_char_pos(msa, msa_idx, 0, 0));
      }
    }
    
    str_free(refseq);
    sfree(fasthash);
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


MSA *maf_read_old(FILE *f, FILE *REFSEQF, int tuple_size, char *alphabet,
		  GFF_Set *gff, CategoryMap *cm, int cycle_size, int store_order,
		  char *reverse_groups, int gap_strip_mode, int keep_overlapping) {
  return maf_read_unsorted(f, REFSEQF, tuple_size, alphabet, gff, cm, cycle_size,
		    store_order, reverse_groups, gap_strip_mode, keep_overlapping, NULL);
}


/* Read a block from an MAF file and store it as a "mini-msa" using
   the provided object.  Allocates memory for sequences if they are
   NULL (as with first block).  Reads to next "a" line or EOF.
   Returns EOF when no more alignments are available.  Sets start
   coord and length of reference sequence if non-NULL
   pointers are provided.  Uses provided hash to map sequence names to
   sequence indices (prefix of name wrt '.' character); sequences not
   present in a block will be represented by missing-data characters. */
int maf_read_block_addseq(FILE *F, MSA *mini_msa, Hashtable *name_hash, 
			  int *start_idx, int *length, int do_toupper,
			  int skip_new_species) {

  int seqidx, more_blocks = 0, i, j;
  String *this_seq, *linebuffer = str_new(STR_VERY_LONG_LEN);
  List *l = lst_new_ptr(7);
  String *this_name = str_new(STR_SHORT_LEN);
  int *mark;

  mini_msa->length = -1;
  mark = smalloc(mini_msa->nseqs*sizeof(int));
  for (i = 0; i < mini_msa->nseqs; i++) mark[i] = 0;
  while (str_readline(linebuffer, F) != EOF) {
    if (str_starts_with_charstr(linebuffer, "#") ||
        str_starts_with_charstr(linebuffer, "i ") ||
        str_starts_with_charstr(linebuffer, "e ") ||
        str_starts_with_charstr(linebuffer, "q ")) 
      continue;                 /* ignore i, e, and q lines for now */
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
    if (lst_size(l) != 7 || !str_equals_charstr(lst_get_ptr(l, 0), "s")) 
      die("ERROR: bad sequence line in MAF file --\n\t\"%s\"\n", linebuffer->chars);
    str_cpy(this_name, lst_get_ptr(l, 1));
    str_shortest_root(this_name, '.');
    this_seq = lst_get_ptr(l, 6);

    /* if this is the reference sequence, also grab start_idx and
       length and check strand */
    if (mini_msa->length == -1 && 
        ((start_idx != NULL && str_as_int(lst_get_ptr(l, 2), start_idx) != 0) ||
        (length != NULL && str_as_int(lst_get_ptr(l, 3), length) != 0) ||
        ((String*)lst_get_ptr(l, 4))->chars[0] != '+')) {
      die("ERROR: bad integers or strand in MAF (strand must be + for reference sequence) --\n\t\"%s\"\n", linebuffer->chars);
    }

    /* ensure lengths of all seqs are consistent */
    if (mini_msa->length == -1) 
      mini_msa->length = this_seq->length;
    else if (this_seq->length != mini_msa->length) {
      die("ERROR: sequence lengths do not match in MAF block -- \n\tsee line \"%s\"\n", linebuffer->chars);
    }

    /* obtain index of seq */
    seqidx = hsh_get_int(name_hash, this_name->chars);
    if (seqidx == -2 || (seqidx == -1 && !skip_new_species)) {
      seqidx = msa_add_seq(mini_msa, this_name->chars);
      hsh_put_int(name_hash, this_name->chars, seqidx);
      mark = srealloc(mark, mini_msa->nseqs*sizeof(int));
    } else if (seqidx == -1) {
      goto msa_read_block_addseq_free_loop;
    }
    if (!(str_equals_charstr(this_name, mini_msa->names[seqidx])))
      die("ERROR: maf_read_block_addseq: %s != %s\n",
	  this_name->chars, mini_msa->names[seqidx]);


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

    for (i = 0; i < this_seq->length; i++) {
      mini_msa->seqs[seqidx][i] = do_toupper ? (char)toupper(this_seq->chars[i]) : 
        this_seq->chars[i];
      if (mini_msa->seqs[seqidx][i] == '.' && mini_msa->inv_alphabet[(int)'.'] == -1) 
        mini_msa->seqs[seqidx][i] = mini_msa->missing[0];
      if (mini_msa->seqs[seqidx][i] != GAP_CHAR && 
          !mini_msa->is_missing[(int)mini_msa->seqs[seqidx][i]] &&
          mini_msa->inv_alphabet[(int)mini_msa->seqs[seqidx][i]] == -1 &&
	  get_iupac_map()[(int)mini_msa->seqs[seqidx][i]] == NULL) {
        if (isalpha(mini_msa->seqs[seqidx][i]))
          mini_msa->seqs[seqidx][i] = 'N';
        else 
          die("ERROR: unrecognized character in sequence in MAF block ('%c')\n",
              mini_msa->seqs[seqidx][i]);
      }
    }
    fflush(stdout);
    mini_msa->seqs[seqidx][this_seq->length] = '\0';
    mark[seqidx] = 1;
  msa_read_block_addseq_free_loop:
    for (i = 0; i < lst_size(l); i++) str_free(lst_get_ptr(l, i));
  }

  lst_free(l);
  str_free(linebuffer);
  str_free(this_name);

  if (mini_msa->length == -1 && !more_blocks) {
    sfree(mark);
    return EOF;}                 /* in this case, an EOF must have been
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
  sfree(mark);
  return 0;
}




/* Read a block from an MAF file and store it as a "mini-msa" using
   the provided object.  Allocates memory for sequences if they are
   NULL (as with first block).  Reads to next "a" line or EOF.
   Returns EOF when no more alignments are available.  Sets start
   coord and length of reference sequence if non-NULL
   pointers are provided.  Uses provided hash to map sequence names to
   sequence indices (prefix of name wrt '.' character); sequences not
   present in a block will be represented by missing-data characters. */
int maf_read_block(FILE *F, MSA *mini_msa, Hashtable *name_hash,
                   int *start_idx, int *length, int do_toupper) {

  int seqidx, more_blocks = 0, i, j;
  String *this_seq, *linebuffer = str_new(STR_VERY_LONG_LEN);
  List *l = lst_new_ptr(7);
  String *this_name = str_new(STR_SHORT_LEN);
  int mark[mini_msa->nseqs];

  mini_msa->length = -1;
  for (i = 0; i < mini_msa->nseqs; i++) mark[i] = 0;
  while (str_readline(linebuffer, F) != EOF) {
    if (str_starts_with_charstr(linebuffer, "#") ||
        str_starts_with_charstr(linebuffer, "i ") ||
        str_starts_with_charstr(linebuffer, "e ") ||
        str_starts_with_charstr(linebuffer, "q ")) 
      continue;                 /* ignore i, e, and q lines for now */
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
      die("ERROR: bad sequence line in MAF file --\n\t\"%s\"\n", linebuffer->chars);
    }
    str_cpy(this_name, lst_get_ptr(l, 1));
    str_shortest_root(this_name, '.');
    this_seq = lst_get_ptr(l, 6);

    /* if this is the reference sequence, also grab start_idx and
       length and check strand */
    if (mini_msa->length == -1 && 
        ((start_idx != NULL && str_as_int(lst_get_ptr(l, 2), start_idx) != 0) ||
        (length != NULL && str_as_int(lst_get_ptr(l, 3), length) != 0) ||
        ((String*)lst_get_ptr(l, 4))->chars[0] != '+'))
      die("ERROR: bad integers or strand in MAF (strand must be + for reference sequence) --\n\t\"%s\"\n", linebuffer->chars);

    /* ensure lengths of all seqs are consistent */
    if (mini_msa->length == -1) mini_msa->length = this_seq->length;
    else if (this_seq->length != mini_msa->length) 
      die("ERROR: sequence lengths do not match in MAF block -- \n\tsee line \"%s\"\n", linebuffer->chars);

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
    seqidx = hsh_get_int(name_hash, this_name->chars);
    if (seqidx == -1) 
      die("ERROR: unexpected sequence name '%s' --\n\tsee line \"%s\"\n", this_name->chars, linebuffer->chars);
    if (!(str_equals_charstr(this_name, mini_msa->names[seqidx])))
      die("ERROR: maf_read_block: %s != %s\n", this_name->chars, mini_msa->names[seqidx]);

    for (i = 0; i < this_seq->length; i++) {
      mini_msa->seqs[seqidx][i] = do_toupper ? (char)toupper(this_seq->chars[i]) : 
        this_seq->chars[i];
      if (mini_msa->seqs[seqidx][i] == '.' && mini_msa->inv_alphabet[(int)'.'] == -1) 
        mini_msa->seqs[seqidx][i] = mini_msa->missing[0];
      if (mini_msa->seqs[seqidx][i] != GAP_CHAR && 
          !mini_msa->is_missing[(int)mini_msa->seqs[seqidx][i]] &&
          mini_msa->inv_alphabet[(int)mini_msa->seqs[seqidx][i]] == -1 &&
	  get_iupac_map()[(int)mini_msa->seqs[seqidx][i]] == NULL) {
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


/* Scan the first block of a MAF file to get a partial list of 
   sequence names (only roots of names are considered), and fill
   the hashtable with these names and a corresponding sequence
   index.  Also get the length of refseq.  If add_seqs==0 will
   not add any new seqs to names or name_hash (but still may
   re-order names to put refseq first) */
void maf_quick_peek(FILE *F, char ***names, Hashtable *name_hash, int *nseqs, int *refseqlen, int add_seqs) {
  String *line = str_new(STR_VERY_LONG_LEN);
  int count = 0, seqidx = 0, tmp, startidx, i, j, length, linenum=0;
  String *fullname = str_new(STR_SHORT_LEN), *name = str_new(STR_SHORT_LEN);
  List *l = lst_new_ptr(7);
  fpos_t pos;
  int hash_val, inParens;

  *refseqlen = -1;

  if (fgetpos(F, &pos) != 0)
    die("ERROR: Currently, MAF input stream must be seekable (can't be stdin).\n");

  while (str_readline(line, F) != EOF) {
    linenum++;
    //if second line is comment, assume it contains parameters as 
    //described by UCSC specification.  If a tree is found, use it
    if (linenum==2 && line->chars[0]=='#') {
      inParens=0;
      for (i=1; i<line->length; i++) {
	if (line->chars[i]=='(') inParens++;
	else if (line->chars[i]==')') inParens--;
	else if (inParens && !isspace(line->chars[i])) {
	  str_clear(fullname);
	  for (j=i; j <line->length; j++) {
	    if (line->chars[j]==')' || isspace(line->chars[j])) break;
	    str_append_char(fullname, line->chars[j]);
	  }
	  str_cpy(name, fullname);
	  str_shortest_root(name, '.');
	  if (name->length <= 0)
	    die("ERROR maf_quick_peek name->length=%i\n", name->length);
	  if (hsh_get_int(name_hash, name->chars) == -1) {
	    if (add_seqs) {
	      hsh_put_int(name_hash, name->chars, count);
	      *names = srealloc(*names, (count+1) * sizeof(char*));
	      (*names)[count] = copy_charstr(name->chars);
	    }
	    count++;
	  }
	  i=j-1;
	}
      }
    }
    if (line->chars[0] == 'a' && seqidx > 0) break;  //break after first block is read
    if (line->chars[0] == 's') {
      /* avoid calling str_split for efficiency */
      for (i = 1; i < line->length && isspace(line->chars[i]); i++);
      str_clear(fullname);
      for (; i < line->length && !isspace(line->chars[i]); i++)
        str_append_char(fullname, line->chars[i]);
      str_cpy(name, fullname);
      str_shortest_root(name, '.');
      if (name->length <= 0)
	die("ERROR maf_quick_peek2: name->length=%i\n", name->length);

      if (hsh_get_int(name_hash, name->chars) == -1 && add_seqs) {
	hsh_put_int(name_hash, name->chars, count);
	*names = srealloc(*names, (count+1) * sizeof(char*));
	(*names)[count] = copy_charstr(name->chars);
	count++;
      }

      if (seqidx == 0) { /* reference sequence */
	if (hsh_get_int(name_hash, name->chars) == -1) {
	  die("cannot disregard reference species (%s) when reading MAF file", name->chars);
	}
        str_split(line, NULL, l);
        if (lst_size(l) != 7 || 
            str_as_int(lst_get_ptr(l, 2), &startidx) != 0 ||
            str_as_int(lst_get_ptr(l, 3), &length) != 0 ||
            str_as_int(lst_get_ptr(l, 5), &tmp) != 0) 
          die("ERROR: bad line in MAF file --\n\t\"%s\"\n", line->chars);

        if (*refseqlen == -1) *refseqlen = tmp;

	//make sure reference sequence is first sequence.  Swap if necessary.
	hash_val = hsh_get_int(name_hash, name->chars);
	if (hash_val != 0) {
	  char *tempCharPtr = (*names)[hash_val];
	  (*names)[hash_val] = (*names)[0];
	  (*names)[0] = tempCharPtr;
	  hsh_reset_int(name_hash, (*names)[hash_val], hash_val);
	  hsh_reset_int(name_hash, name->chars, 0);
	}
      }
      seqidx++;
    }
  }
  fsetpos(F, &pos);
  str_free(line); str_free(fullname); str_free(name);
  lst_free_strings(l); lst_free(l);
  if (nseqs != NULL) *nseqs = count;
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
      str_shortest_root(name, '.');
      if (name->length <= 0)  /* must be a non-empty name */
	die("ERROR: maf_peek: name->length=%i\n", name->length);
      if (hsh_get_int(name_hash, name->chars) == -1) {
        hsh_put_int(name_hash, name->chars, count);
        *names = srealloc(*names, (count+1) * sizeof(char*));
        (*names)[count] = copy_charstr(name->chars);
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
      sfree(gp);
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
  int first_extend = -1;
  GFF_Feature *feat;
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
      if (featcpy->strand != '-' && featcpy->frame != GFF_NULL_FRAME) 
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

    if (first_extend == -1 && feat->end > end_idx) first_extend = *gff_idx;

    /* shift coords and add feature */
    featcpy->start -= (start_idx - 1);
    featcpy->end -= (start_idx - 1);

    lst_push_ptr(sub_gff->features, featcpy);
  }

  if (first_extend != -1) *gff_idx = first_extend;
                                /* this allows features that extend beyond a
                                   block to be considered for the next
                                   block */
}


