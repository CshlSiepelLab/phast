/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: multi_msa.c,v 1.2 2008-11-12 02:07:59 acs Exp $ */

/** \file msa.c
   Multiple sequence alignments.
   \ingroup msa
*/

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
#include <string.h>
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

/** Creates a new multi-fasta object containing alignments read from files
    specified in a config file. Config file contains the number of alignments,
    format, alphabet and list of files to be read. If "alphabet" is
    NULL, default alphabet for DNA will be used.  This routine will
    abort if the sequence contains a character not in the alphabet. */
Multi_MSA *msa_multimsa_new(FILE *F, int do_ih) {

  Regex *blocks_re = str_re_new("#[[:space:]]*BLOCKS[[:space:]]*=[[:space:]]*([0-9]+)");
  Regex *alph_re = str_re_new("#[[:space:]]*ALPHABET[[:space:]]*=[[:space:]]*([A-Z]+)");
  Regex *format_re = str_re_new("#[[:space:]]*FORMAT[[:space:]]*=[[:space:]]*([A-Z]+)");
  
  int i, nblocks, line_no=0;
  char *msa_fname, *ih_fname = NULL;
  FILE *msa_f, *ih_f;
  msa_format_type format=0;
  Multi_MSA *msas = NULL;
  MSA *msa;
  List *matches = lst_new_ptr(4);
  String *line = str_new(STR_MED_LEN), *alphabet = NULL, 
    *fname = str_new(STR_MED_LEN);
  nblocks = i = 0;

  while (str_readline(line, F) != EOF) {
    str_trim(line);
    checkInterruptN(line_no, 1000);
    line_no++;
    if (line->length == 0) continue; /* ignore blank lines */
    
    if (msas == NULL) {
      if (str_re_match(line, blocks_re, matches, 1) >= 0) { 
	str_as_int(lst_get_ptr(matches, 1), &nblocks);
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
      
      if (nblocks != 0 && alphabet != NULL && i == 3) {
	msas = smalloc(sizeof(Multi_MSA));
	msas->nblocks = nblocks;
	msas->blocks = smalloc(nblocks * sizeof(MSA));
	msas->seqnames = lst_new(nblocks, sizeof(String));
	if (do_ih) { /* Allocate space for indel histories, if
			called for -- this will need to be filled
			in explicitly later */
	  msas->ih = smalloc(nblocks * sizeof(IndelHistory));
	} else {
	  msas->ih = NULL;
	}
	i = 0;
      }
    } else {
      if (i == msas->nblocks)
	die("Too many alignment files for format\n");

      if (str_split(line, NULL, matches) > 1) {
	msa_fname = ((String*)lst_get_ptr(matches, 0))->chars;
	ih_fname = ((String*)lst_get_ptr(matches, 1))->chars;
	
      } else {
	msa_fname = line->chars;
      }

      fprintf(stderr, "\t%s (%d of %d)\n", msa_fname, (i+1), msas->nblocks);
      msa_f = fopen_fname(msa_fname, "r");
      msa = msa_new_from_file(msa_f, format, (char*)alphabet->chars);
      str_cpy_charstr(fname, msa_fname);
      str_remove_path(fname);
      str_root(fname, '.');
      lst_push(msas->seqnames, str_dup(fname));
      msas->blocks[i] = msa;
      fclose(msa_f);

      if (do_ih) {
	if (ih_fname == NULL)
	  die("ERROR: --indel-model requires an indel history for all alignment files! Try 'dmsample -h'\n");
	fprintf(stderr, "\tReading indel history for %s from %s\n", msa_fname, 
		ih_fname);
	ih_f = fopen_fname(ih_fname, "r");
	msas->ih[i] = ih_new_from_file(ih_f);
	fclose(ih_f);
      }
      i++;
    }
  }
  if (i < (nblocks - 1))
    die("Not enough files in alignments list\n");
  lst_free(matches);
  str_free(alphabet);
  str_free(line);
  str_free(fname);
  return msas;
}
