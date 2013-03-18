/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/*
   To do:
	msas->ih = NULL;
      - more general coordinate mapping?  Most functions use MSA
      coords, but some can convert using coordinate maps.

      - clean up indexing.  Usually 0-based indexing is used, but
      not always.

      - some obscure aspects of PHYLIP format not supported -- e.g.,
      problem when sequence name is not separated from sequence by
      whitespace.

      - should automatically recognize alignment format when reading.

      - unified handling of suff stats; they've become a special case
      that nearly every function has to consider.  Perhaps should just
      always do *everything* in terms of sufficient statistics, and
      eliminate explicit representation of alignments.
*/

#include <stdlib.h>
#include <ctype.h>
#include <time.h>

#include <lists.h>
#include <stacks.h>
#include <msa.h>
#include <multi_msa.h>
#include <misc.h>
#include <gff.h>
#include <category_map.h>
#include <hashtable.h>
#include <sufficient_stats.h>
#include <local_alignment.h>
#include <indel_history.h>

/* Creates a new multi-fasta object containing alignments read from files
    specified in a config file. Config file contains the number of alignments,
    format, alphabet and list of files to be read. If "alphabet" is
    NULL, default alphabet for DNA will be used.  This routine will
    abort if the sequence contains a character not in the alphabet. */
Multi_MSA *multimsa_new_from_files(FILE *F) {

  Regex *blocks_re = str_re_new("#[[:space:]]*BLOCKS[[:space:]]*=[[:space:]]*([0-9]+)");
  Regex *alph_re = str_re_new("#[[:space:]]*ALPHABET[[:space:]]*=[[:space:]]*([A-Z]+)");
  Regex *format_re = str_re_new("#[[:space:]]*FORMAT[[:space:]]*=[[:space:]]*([A-Z]+)");
  
  int i, num_msa, line_no=0;
  char *msa_fname;
  FILE *msa_f;
  msa_format_type format=0;
  Multi_MSA *msas = NULL;
  MSA *msa;
  List *matches = lst_new_ptr(4);
  String *line = str_new(STR_MED_LEN), *alphabet = NULL, 
    *fname = str_new(STR_MED_LEN);
  num_msa = i = 0;

  while (str_readline(line, F) != EOF) {
    str_trim(line);
    checkInterruptN(line_no, 1000);
    line_no++;
    if (line->length == 0) continue; /* ignore blank lines */
    
    if (msas == NULL) {
      if (str_re_match(line, blocks_re, matches, 1) >= 0) { 
	str_as_int(lst_get_ptr(matches, 1), &num_msa);
	i++;
      } else if (str_re_match(line, alph_re, matches, 1) >= 0) {
	alphabet = str_dup(lst_get_ptr(matches, 1));
	i++;
      } else if (str_re_match(line, format_re, matches, 1) >= 0) {
	format = msa_str_to_format(((String*)lst_get_ptr(matches, 1))->chars);
	i++;
      } else {
	die("Bad header in alignment list file\n");
      }
      
      if (num_msa != 0 && alphabet != NULL && i == 3) {
	msas = smalloc(sizeof(Multi_MSA));
	msas->num_msa = num_msa;
	msas->msa = smalloc(num_msa * sizeof(MSA));
	msas->seqnames = lst_new(num_msa, sizeof(String));
	msas->concat_msa = NULL;
	msas->start_coord = msas->end_coord = NULL;
	msas->type = NULL;
	i = 0;
      }
    } else {
      if (i == msas->num_msa)
	die("Too many alignment files for format\n");

      msa_fname = line->chars;

      fprintf(stderr, "\t%s (%d of %d)\n", msa_fname, (i+1), msas->num_msa);
      msa_f = phast_fopen(msa_fname, "r");
      msa = msa_new_from_file_define_format(msa_f, format, (char*)alphabet->chars);
      str_cpy_charstr(fname, msa_fname);
      str_remove_path(fname);
      str_root(fname, '.');
      lst_push(msas->seqnames, str_dup(fname));
      msas->msa[i] = msa;
      phast_fclose(msa_f);

      i++;
    }
  }
  if (i < (num_msa - 1))
    die("Not enough files in alignments list\n");
  lst_free(matches);
  str_free(alphabet);
  str_free(line);
  str_free(fname);
  return msas;
}


/* Create a multi MSA from an array of MSAs.   Does NOT copy MSAs.
 */
Multi_MSA *multimsa_new_from_msas(MSA **msa, int nmsa) {
  Multi_MSA *rv = smalloc(sizeof(Multi_MSA));
  int i;
  rv->num_msa = nmsa;
  rv->msa = smalloc(nmsa * sizeof(MSA*));
  for (i=0; i < rv->num_msa; i++)
    rv->msa[i] = msa[i];
  rv->seqnames = NULL;
  rv->concat_msa = NULL;
  rv->start_coord = rv->end_coord = NULL;
  rv->type = NULL;
  return rv;
}


void multimsa_make_concat(Multi_MSA *mmsa) {
  int i;
  if (mmsa->concat_msa != NULL) return;  /* assume this one is correct? */
  if (mmsa->num_msa == 0) {
    mmsa->concat_msa = NULL;
    return;
  }
  mmsa->concat_msa = msa_create_copy(mmsa->msa[0], 0);
  for (i=1; i < mmsa->num_msa; i++) 
    msa_concatenate(mmsa->concat_msa, mmsa->msa[i]);
}


void multimsa_add_msa(Multi_MSA *mmsa, MSA *msa) {
  if (mmsa->num_msa == 0) 
    mmsa->msa = smalloc(sizeof(MSA*));
  else mmsa->msa = srealloc(mmsa->msa, (mmsa->num_msa+1)*sizeof(MSA*));
  mmsa->msa[mmsa->num_msa] = msa;
  if (mmsa->concat_msa != NULL) 
    msa_concatenate(mmsa->concat_msa, mmsa->msa[mmsa->num_msa]);
  mmsa->num_msa++;
}


