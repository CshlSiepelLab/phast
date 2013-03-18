/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: clean_genes.c,v 1.42 2008-11-12 02:07:58 acs Exp $ */

#include <stdlib.h>
#include <stdio.h>
#include <msa.h>
#include <gff.h>
#include <sufficient_stats.h>
#include <getopt.h>
#include <maf.h>
#include <external_libs.h>
#include "clean_genes.help"

/* types of features examined */
#define SPLICE_5 "5'splice"
#define SPLICE_3 "3'splice"
#define SPLICE_5_UTR "5'splice_utr"
#define SPLICE_3_UTR "3'splice_utr"
/* NOTE: SPLICE_5_UTR and SPLICE_3_UTR must have SPLICE_5 and SPLICE_3
   as prefixes (respectively).  This property is used to simplify some
   of the case handling below */ 

/* possible values for status of each group of features */
/* NOTE: these are now more important as individual problem types than
   as a global status for a whole gene */
typedef enum {OKAY, BAD_REF, BAD_REF_START, BAD_REF_STOP, BAD_REF_5_SPLICE,
              BAD_REF_3_SPLICE, BAD_REF_ORF, BAD_REF_INDEL_STRICT_FAIL,
              NO_ALN, BAD_START, BAD_STOP, BAD_5_SPLICE, BAD_3_SPLICE, 
              BAD_5_SPLICE_UTR, BAD_3_SPLICE_UTR, NONSENSE, FSHIFT, 
              BAD_INTRON, TOO_MANY_Ns, WARN_FSHIFT, WARN_Ns, NTYPES} 
status_type;                    /* NTYPES marks the total number */

/* possible gap types for cds exon -- frame-shift gaps (FSHIFT_BAD),
   no apparent frame shift (FSHIFT_OK), all gaps multiple of 3 in len
   (CLN_GAPS), gaps mult. of 3 and nonoverlapping (NOVRLP_CLN_GAPS),
   and no gaps present (NGAPS).  Defined in order of least to most
   stringent criterion (ordering gets used) */
typedef enum {FSHIFT_BAD, FSHIFT_OK, CLN_GAPS, 
              NOVRLP_CLN_GAPS, NGAPS, NGAP_TYPES} 
cds_gap_type;                   /* NGAP_TYPES marks the total number */

/* for use in adjusting cds's to exclude stop codons (see below) */
typedef enum {NO_ADJUST, ADJUST_END, ADJUST_START} adjust_type;

/* minimum number of sites required for a gapless alignment block,
   when checking for frame-shift indels -- approximately the same as
   the maximum number of bases that can separate compensatory
   indels (see function 'is_fshift_okay') */
#define MIN_GAPLESS_BLOCK_SIZE 15

/* in is_fshift_okay, maximum size of a "gappy region" of an alignment
   that is tolerated, even if there's no net frame shift.  Used to
   catch cases where frame has been conserved by random chance */
#define MAX_GAPPY_BLOCK_SIZE 45

/* "signal" feature types, used for 'indel-strict' mode (special-case
   stuff related to exoniphy training) */
#define SIGNALS "start_codon,stop_codon,5'splice,3'splice,cds5'ss,cds3'ss,prestart"

/* description of stats kept, written to bottom of stats file */
char *STATS_DESCRIPTION = "#\n\
# *** PASS/FAIL STATS ***\n\
# total         total number of entries in input file\n\
# nbad_ref      number failing test against reference sequence\n\
# nconsid       number considered in other tests (total - nbad_ref)\n\
# nkept         number passing all tests\n\
# nno_aln       number failing alignment test (completely missing in at least one species)\n\
# nbad_starts   number failing start codon test\n\
# (out of)      (number tested)\n\
# nbad_stops    number failing stop codon test\n\
# (out of)      (number tested)\n\
# nbad_5spl     number failing 5' splice site test (not designated UTR)\n\
# (out of)      (number tested)\n\
# nbad_3spl     number failing 3' splice site test (not designated UTR)\n\
# (out of)      (number tested)\n\
# nbad_5utr     number failing 5' splice site test (designated UTR)\n\
# (out of)      (number tested)\n\
# nbad_3utr     number failing 3' splice site test (designated UTR)\n\
# (out of)      (number tested)\n\
# nbad_intron   number of pairs of splice sites failing consistency test (not GT-AG, GC-AG, AT-AC)\n\
# nnons         number of cds's failing nonsense mutation test (not counting ones that fail above tests)\n\
# nfshifts      number of cds's failing frameshift test (not counting ones that fail above tests)\n\
# nNs           number of cds's exceeding limit on number of Ns\n\
# *** ALIGNMENT GAP STATS (for use at individual exon level) ***\n\
# ncons_exons   number of conserved exons (same as nkept if --conserved)\n\
# nce_ngaps     number with no alignment gaps\n\
# nce_nov_cln   number with only nonoverlapping \"clean\" gaps (multiple of three lengths)\n\
# nce_clean     number with clean gaps that overlap\n\
# nce_fshftok   number with compensatory frame-shifting gaps as allowed by --fshift\n";

/* description of a problem */
typedef struct Problem {
  GFF_Feature *feat;
  status_type status;
  int start;
  int end;
  cds_gap_type cds_gap;  /* if status if FSHIFT */
} Problem;

/* create a new problem.  feat can be null for whole gene.
 * start and end are < 0, then they are filled in from feat   */
Problem *problem_new(GFF_Feature *feat, status_type status,
                     int start, int end) {
  Problem *p = smalloc(sizeof(Problem));
  p->feat = feat;
  p->status = status;
  if ((feat != NULL) && (start < 0) && (end < 0)) {
    p->start = feat->start;
    p->end = feat->end;
  } else {
    p->start = start;
    p->end = end;
  }
  return p;
}

/* create a new problem, and add to the list */
Problem *problem_add(List *problems, GFF_Feature *feat, status_type status,
                     int start, int end) {
  Problem *p = problem_new(feat, status, start, end);
  lst_push_ptr(problems, p);
  return p;
}

/* free a problem */
void problem_free(Problem *p) {
  if (p != NULL) {
    sfree(p);
  }
}

/* Reset a problem list to the empty state */
void problems_clear(List *problems) {
  int i;
  for (i = 0; i < lst_size(problems); i++) {
    problem_free(lst_get_ptr(problems, i));
  }
  lst_clear(problems);
}

/* free list of problem objects */
void problems_free(List *problems) {
  problems_clear(problems);
  lst_free(problems);
}

static PHAST_INLINE 
int is_signal(char *str, int n_sig, char **signals, int *is_missing) {
  /* determine if str equals one of the strings in signals,
   * possibly with replacing some chars with missing values*/
  int i, j;
  int ok;

  for (i = 0; i < n_sig; i++) {
    ok = 1;
    for (j = 0; j < strlen(signals[i]); j++) {
      if (! is_missing[(int)str[j]] && str[j] != signals[i][j]) { 
	ok = 0;
	break;
      }
    }
    if (ok) return 1;
  }
  return 0;
}


static PHAST_INLINE 
int is_conserved_start(GFF_Feature *feat, MSA *msa) {
  char tuplestr[4];
  int j;
  int start = feat->start - 1;  /* base 0 indexing */
  tuplestr[3] = '\0';
  char * start_signal[1] = { "ATG" };

  for (j = 0; j < msa->nseqs; j++) {
    tuplestr[0] = ss_get_char_tuple(msa, msa->ss->tuple_idx[start], j, 0);
    tuplestr[1] = ss_get_char_tuple(msa, msa->ss->tuple_idx[start+1], j, 0);
    tuplestr[2] = ss_get_char_tuple(msa, msa->ss->tuple_idx[start+2], j, 0);
    if (feat->strand == '-') msa_reverse_compl_seq(tuplestr, 3);    
    if(! is_signal(tuplestr, 1, start_signal, msa->is_missing)) return 0;
  }
  return 1;
}

static PHAST_INLINE 
int is_stop_codon(char *str) {
 return (strncmp(str, "TAA", 3) == 0 || strncmp(str, "TAG", 3) == 0 ||
         strncmp(str, "TGA", 3) == 0);
 /* strncmp allows testing of codons in the middle of larger strings */
}


static PHAST_INLINE 
int is_conserved_stop(GFF_Feature *feat, MSA *msa) {
  char tuplestr[4];
  int j;
  int start = feat->start - 1;
  tuplestr[3] = '\0';
  char * stop_signals[3] = { "TAA", "TAG", "TGA" };
  for (j = 0; j < msa->nseqs; j++) {
    tuplestr[0] = ss_get_char_tuple(msa, msa->ss->tuple_idx[start], j, 0);
    tuplestr[1] = ss_get_char_tuple(msa, msa->ss->tuple_idx[start+1], j, 0);
    tuplestr[2] = ss_get_char_tuple(msa, msa->ss->tuple_idx[start+2], j, 0);
    if (feat->strand == '-') msa_reverse_compl_seq(tuplestr, 3);
    if (!is_signal(tuplestr, 3, stop_signals, msa->is_missing)) return 0;
  }
  return 1;
}

/* returns 1 if GT, GC or AT, and 0 otherwise.  If splice_strict,
   returns 0 unless GT. */
static PHAST_INLINE 
int is_valid_5splice(char *str, int splice_strict) {
  if (strncmp(str, "GT", 2) == 0) return 1;
  else if (splice_strict) return 0;
  else if (strncmp(str, "GC", 2) == 0 || strncmp(str, "AT", 2) == 0)
    return 1;
  return 0;
}

static PHAST_INLINE 
int is_conserved_5splice(GFF_Feature *feat, MSA *msa, int offset5,
			 int splice_strict) {
  char tuplestr[3];
  int j, start = feat->start - 1;  /* base 0 indexing */
  if (feat->strand == '-') start += offset5;
  tuplestr[2] = '\0';
  char * splice_signals[3] = {"GT", "GC", "AT"};

  for (j = 0; j < msa->nseqs; j++) {
    tuplestr[0] = ss_get_char_tuple(msa, msa->ss->tuple_idx[start], j, 0);
    tuplestr[1] = ss_get_char_tuple(msa, msa->ss->tuple_idx[start+1], j, 0);
    if (feat->strand == '-') msa_reverse_compl_seq(tuplestr, 2);
    if (!is_signal(tuplestr, (splice_strict?1:3), splice_signals, msa->is_missing)) 
      return 0;
  }
  return 1;
}

/* returns 1 for AG or AC, 0 otherwise.  If splice_strict, returns 0
   unless AG. */
static PHAST_INLINE 
int is_valid_3splice(char *str, int splice_strict) {
  if (strncmp(str, "AG", 2) == 0) return 1;
  else if (splice_strict) return 0;
  else if (strncmp(str, "AC", 2) == 0) return 1;
  return 0;
}

static PHAST_INLINE 
int is_conserved_3splice(GFF_Feature *feat, MSA *msa, int offset3,
			 int splice_strict) {
  char tuplestr[3];
  int j, start = feat->start - 1;  /* base 0 indexing */
  if (feat->strand == '+') start += offset3;
  tuplestr[2] = '\0';
  char * splice_signals[2] = {"AG", "AC"};
  for (j = 0; j < msa->nseqs; j++) {
    tuplestr[0] = ss_get_char_tuple(msa, msa->ss->tuple_idx[start], j, 0);
    tuplestr[1] = ss_get_char_tuple(msa, msa->ss->tuple_idx[start+1], j, 0);
    if (feat->strand == '-') msa_reverse_compl_seq(tuplestr, 2);
    if (!is_signal(tuplestr, (splice_strict?1:2), splice_signals, msa->is_missing)) 
      return 0;
  }
  return 1;
}

int is_nonsense_clean(GFF_Feature *feat, MSA *msa, List *problems) {
  int i, j, len;
  char seq[feat->end - feat->start + 2];
  for (j = 1; j < msa->nseqs; j++) { /* no need to check reference seq */
    /* first copy entire sequence without gaps */
    for (i = feat->start - 1, len = 0; i < feat->end; i++) 
      if (ss_get_char_pos(msa, i, j, 0) != GAP_CHAR)
        seq[len++] = ss_get_char_pos(msa, i, j, 0);
    seq[len] = '\0';

    if (feat->strand == '-') msa_reverse_compl_seq(seq, len);

    /* now scan for stop codons */
    for (i = (3 - feat->frame) % 3; i <= len - 3; i += 3) 
      if (is_stop_codon(&seq[i])) {
	int problem_start;
	if(feat->strand == '+') problem_start = feat->start+i;
	else problem_start = feat->end-i-2;

        problem_add(problems, feat, NONSENSE, problem_start,
                    problem_start+2);
        return 0;
      }
  }

  return 1;
}


/* for use in scan_for_gaps */
struct gap { int start; int end; };
int gap_compare(const void *ptr1, const void* ptr2) {
  struct gap *g1 = *(struct gap**)ptr1;
  struct gap *g2 = *(struct gap**)ptr2;
  return g1->start - g2->start;
}
/* scans a cds for gaps.  Returns CLN_GAPS, NOVRLP_CLN_GAPS, NO_GAPS,
   or FSHIFT_BAD; doesn't try to check for compensatory indels, which
   is more complicated (this is left for the special-purpose function
   below) */
int scan_for_gaps(GFF_Feature *feat, MSA *msa, Problem **problem) {
  int msa_start = feat->start - 1;
  int msa_end = feat->end - 1;
  int i, j;
  int near_boundary = 0;
  cds_gap_type retval = NGAPS;
  List *gaps = lst_new_ptr(10);

  for (j = 0; retval != FSHIFT_BAD && j < msa->nseqs; j++) {
    for (i = msa_start; i <= msa_end; i++) {
      if (ss_get_char_pos(msa, i, j, 0) == GAP_CHAR) {
        int gap_start, gap_end;
        struct gap *g;

        for (gap_start = i-1; gap_start >= msa_start && 
               ss_get_char_pos(msa, gap_start, j, 0) == GAP_CHAR; gap_start--);
        gap_start++;            /* inclusive */
        for (gap_end = i+1; gap_end <= msa_end && 
               ss_get_char_pos(msa, gap_end, j, 0) == GAP_CHAR; gap_end++);
        gap_end--;              /* inclusive */

        if ((gap_end - gap_start + 1) % 3 != 0) {
          retval = FSHIFT_BAD;
          *problem = problem_new(feat, FSHIFT, gap_start, gap_end);
          (*problem)->cds_gap = FSHIFT_BAD;
          break;
        }

        /* note whether gaps occur near a cds boundary (within 3 sites) */
        if (gap_start <= msa_start + 3 || gap_end >= msa_end - 3)
          near_boundary = 1;
        
        if (retval == NGAPS) retval = CLN_GAPS;
        g = smalloc(sizeof(struct gap));
        g->start = gap_start;
        g->end = gap_end;
        lst_push_ptr(gaps, g);

        i = gap_end;
      }
    }
  }

  if (retval == CLN_GAPS) {     /* now check for overlaps */
    lst_qsort(gaps, gap_compare);
    retval = NOVRLP_CLN_GAPS;
    for (i = 1; i < lst_size(gaps); i++) {
      struct gap *g1 = lst_get_ptr(gaps, i-1);
      struct gap *g2 = lst_get_ptr(gaps, i);
      if (g2->start <= g1->end && 
          (g2->start != g1->start || g2->end != g1->end)) {
        retval = CLN_GAPS;
        break;
      }
    }
    if (retval == NOVRLP_CLN_GAPS && near_boundary)
      retval = CLN_GAPS;        /* note that the boundary criterion is
                                   being confounded with the overlap
                                   criterion.  Doesn't seem worth
                                   fixing at the moment ...  */
  }

  for (i = 0; i < lst_size(gaps); i++) sfree(lst_get_ptr(gaps, i));
  lst_free(gaps);
  return retval;
}

/* look for frame-shift gaps using a slightly more sophisticated
   algorithm, which allows for compensatory indels.  The strategy here
   is to identify maximal gapless blocks of greater than
   MIN_GAPLESS_BLOCK_SIZE sites, then to make sure that in the gappy portions
   between them, each sequence has a total number of gaps that equals
   the total number for the reference sequence, modulo 3.  Returns 1 if
   all gaps look okay (no net frame shift) and 0 otherwise. */
int is_fshift_okay(GFF_Feature *feat, MSA *msa) {
  int *ngaps = smalloc(msa->nseqs * sizeof(int));
  int i, j, blk_beg, blk_end, start_gappy_reg;

  for (j = 0; j < msa->nseqs; j++) ngaps[j] = 0;
  start_gappy_reg = 0;

  for (i = feat->start - 1; i < feat->end; ) {
    /* find next gapless column, simultaneously keeping track of the
       number of gaps encountered in each sequence */
    for (; i < feat->end; i++) {
      int gapless_col = 1;  
      for (j = 0; j < msa->nseqs; j++) {
        if (ss_get_char_pos(msa, i, j, 0) == GAP_CHAR) {
          ngaps[j]++;
          gapless_col = 0;
        }
      }
      if (gapless_col == 1)
        break;     
    }

    blk_beg = i;                /* inclusive */

    /* find next col with gap */
    for (i++; i < feat->end; i++) {
      for (j = 0; 
           j < msa->nseqs && ss_get_char_pos(msa, i, j, 0) != GAP_CHAR; 
           j++);
      if (j != msa->nseqs) break;
    }
    blk_end = i;                /* exclusive */

    if (blk_end - blk_beg >= MIN_GAPLESS_BLOCK_SIZE ||
        blk_beg == feat->end || /* gaps at end of aln */
        blk_end == feat->end) { /* short block at end of aln */
      /* check total number of gaps since last retained block or
         beginning of alignment; must be same as reference sequence,
         mod 3 */
      for (j = 0; j < msa->nseqs && ngaps[j] % 3 == ngaps[0] % 3; j++);
      /* reject alignment if mod 3 test fails OR if the total length
         of the gappy region exceeds MAX_GAPPY_BLOCK_SIZE */
      if (j != msa->nseqs || blk_beg - start_gappy_reg > MAX_GAPPY_BLOCK_SIZE) {
        sfree(ngaps);
        return 0;
      }

      /* reset ngaps (note: done only if block exceeds size
         threshold) */
      for (j = 0; j < msa->nseqs; j++) ngaps[j] = 0;
      start_gappy_reg = blk_end;
    }
  }
  sfree(ngaps);
  return 1;
}

/* Given a cds feature, determine whether it has no gaps (NGAPS),
   "clean" gaps (all multiples of 3 in length; CLEAN_GAPS)
   non-overlapping clean gaps (NOVRLP_CLN_GAPS), "okay" gaps (only
   temporary frame shifts, corrected by compensatory indels;
   FSHIFT_OK), or real frame-shift gaps (FSHIFT_BAD) */
cds_gap_type get_cds_gap_type(GFF_Feature *feat, MSA *msa, List *problems) {
  Problem *problem = NULL;
  cds_gap_type retval = scan_for_gaps(feat, msa, &problem);

  if (retval == FSHIFT_BAD && is_fshift_okay(feat, msa)) {
    retval = FSHIFT_OK;
                                /* most of the time the call to
                                   is_fshift_okay won't be
                                   necessary */
    problem->status = WARN_FSHIFT;
    problem->cds_gap = FSHIFT_OK;
  }
  if (problem != NULL) {
    lst_push_ptr(problems, problem);
    /* FIXME: It's possible that the single problem identified in
       scan_for_gaps is actually okay, but there's a frameshift
       without compensation downstream.  In this case, the status will
       be correct but the problem will point to the wrong place */
  }
  return retval;
}

/* returns TRUE if the alignment is incomplete (i.e., contains '*'
   characters) or if any row consists completely of gaps in the region
   of the specified feature, and returns FALSE otherwise. */
int is_incomplete_alignment(GFF_Feature *feat, MSA *msa) {
  int i, j;
  for (i = 1; i < msa->nseqs; i++) { /* don't check reference seq */
    int row_all_gaps = TRUE;
    for (j = feat->start - 1; row_all_gaps && j < feat->end; j++) {
      char c = ss_get_char_tuple(msa, msa->ss->tuple_idx[j], i, 0); 
      if (msa->is_missing[(int)c] && c != 'N') /* N is special case; ignore */
        return TRUE;
      if (c != GAP_CHAR)
        row_all_gaps = FALSE;
    }
    if (row_all_gaps) return TRUE;
    /* NOTE: the test for a row of all gaps is not really necessary
       any longer (generally '*'s will appear rather than gaps in such
       a case) but it doesn't really hurt to leave it in and it may be
       useful in cases where information about missing alignments has
       been lost */
  }
  return FALSE;
}

/* comparators used below in is_intron_okay */
int feature_comparator_ascending(const void* ptr1, const void* ptr2) {
  GFF_Feature *feat1 = *((GFF_Feature**)ptr1);
  GFF_Feature *feat2 = *((GFF_Feature**)ptr2);
  if (feat1->start != feat2->start) 
    return (feat1->start - feat2->start);
  return (feat1->end - feat2->end);
}
int feature_comparator_descending(const void* ptr1, const void* ptr2) {
  GFF_Feature *feat1 = *((GFF_Feature**)ptr1);
  GFF_Feature *feat2 = *((GFF_Feature**)ptr2);
  if (feat1->start != feat2->start) 
    return (feat2->start - feat1->start);
  return (feat2->end - feat1->end);
}


/* given a list of 5' and 3' splice sites extracted from a group,
   check whether they form valid pairs in all species */
int are_introns_okay(List *intron_splice,  MSA *msa, List *problems,
                     int offset5, int offset3) {
  int i, j, start1, start2;
  char str1[3], str2[3], str12[5];
  char strand;
  int retval = 1;
  char * splice_pairs[3] = {"GTAG", "GCAG", "ATAC"};

  str1[2] = '\0'; str2[2] = '\0';

  if (lst_size(intron_splice) < 2) return 1;

  strand = ((GFF_Feature*)lst_get_ptr(intron_splice, 0))->strand;
                                /* assume all same strand */

  if (strand == '+')
    lst_qsort(intron_splice, feature_comparator_ascending); 
  else
    lst_qsort(intron_splice, feature_comparator_descending); 

  for (i = 0; i < lst_size(intron_splice) - 1; i++) {
    /* assume every 5' splice and immediately following 3' splice
       form a pair */
    GFF_Feature *f1 = lst_get_ptr(intron_splice, i);
    GFF_Feature *f2 = lst_get_ptr(intron_splice, i+1);
    if (str_starts_with_charstr(f1->feature, SPLICE_5) &&
        str_starts_with_charstr(f2->feature, SPLICE_3)) {
      start1 = f1->start - 1 + (strand == '-' ? offset5 : 0);
      start2 = f2->start - 1 + (strand == '+' ? offset3 : 0);
      for (j = 0; j < msa->nseqs; j++) {
        str1[0] = ss_get_char_tuple(msa, msa->ss->tuple_idx[start1], j, 0);
        str1[1] = ss_get_char_tuple(msa, msa->ss->tuple_idx[start1+1], j, 0);
        str2[0] = ss_get_char_tuple(msa, msa->ss->tuple_idx[start2], j, 0);
        str2[1] = ss_get_char_tuple(msa, msa->ss->tuple_idx[start2+1], j, 0);
        if (strand == '-') {
          msa_reverse_compl_seq(str1, 2);
          msa_reverse_compl_seq(str2, 2);
        }
	strcpy(str12, str1); strcat(str12, str2);
        if (!is_signal(str12, 3, splice_pairs, msa->is_missing)) {
          problem_add(problems, f1, BAD_INTRON, -1, -1);
          problem_add(problems, f2, BAD_INTRON, -1, -1);
          retval = 0;
          break;
        }
      }
      i++;                      /* no need to look at next feature */
    }
  }
  return retval;
}

/* add counts of all (non-gap) CDS bases and all Ns to running total
   for a CDS feature */
void get_N_counts(int *countNs, int *countCDSs, GFF_Feature *feat, MSA *msa) {
  int i, j;
  for (i = 0; i < msa->nseqs; i++) {
    for (j = feat->start-1; j < feat->end; j++) {
      char c = ss_get_char_tuple(msa, msa->ss->tuple_idx[j], i, 0);
      if (c != GAP_CHAR)
        countCDSs[i]++;
      if (c == 'N')
        countNs[i]++;
    }
  }
}

/* dump sub-alignment corresponding to feature */
void dump_aln(FILE *F, GFF_Feature *feat, MSA *msa, int show_frame) {
  int i, j;
  MSA *sub = msa_sub_alignment(msa, NULL, 0, feat->start - 1, feat->end);
  if (sub->seqs == NULL && sub->ss != NULL) ss_to_msa(sub);
  if (feat->strand == '-') msa_reverse_compl(sub);
  for (i = 0; i < sub->nseqs; i++) {
    fprintf(F, "%-12s ", sub->names[i]);
    for (j = 0; j < sub->length; j++) {
      if (show_frame && (j + feat->frame) % 3 == 0) {
        if (is_stop_codon(&(sub->seqs[i][j]))) fprintf(F, "*");
        /* this may not be right when there are frameshift gaps (e.g.,
           temporary ones with compensatory indels) or gaps
           interrupting stop codons, but it will be good enough to see
           what's going on most of the time  */
        else fprintf(F, ".");
      }
      fprintf(F, "%c", sub->seqs[i][j]);
    }
    fprintf(F, "\n");
  }
  msa_free(sub);
}

/* write log entry for discarded feature */
void write_log(FILE *logf, GFF_FeatureGroup *group, status_type status, 
               List *problems, MSA *msa, msa_coord_map *map) {

  int i;

  /* special case */
  if (status == NO_ALN) {
    fprintf(logf, "****\nSkipped '%s' -- lack of alignment across species.\n", 
            group->name->chars);
    return;
  }

  if (status != OKAY) 
    fprintf(logf, "****\nDiscarded '%s'\n", group->name->chars);
  else
    fprintf(logf, "****\nWarnings for '%s'\n", group->name->chars);

  for (i = 0; i < lst_size(problems); i++) {
    struct Problem *problem = lst_get_ptr(problems, i);
    char *reason=NULL;
    int print_alignment = TRUE;

    switch (problem->status) {
    case BAD_REF_START:
      reason = "Bad start codon in reference sequence";
      break;
    case BAD_REF_STOP:
      reason = "Bad stop codon in reference sequence";
      break;
    case BAD_REF_5_SPLICE:
      reason = "Bad 5' splice in reference sequence";
      break;
    case BAD_REF_3_SPLICE:
      reason = "Bad 3' splice in reference sequence";
      break;
    case BAD_REF_ORF:
      reason = "Premature stop codon in reference sequence";
      break;
    case BAD_REF_INDEL_STRICT_FAIL:
      reason = "Indel-strict failure in reference sequence";
      break;
    case BAD_START:
      reason = "Start not conserved";
      break;
    case BAD_STOP:
      reason = "Stop not conserved";
      break;
    case BAD_5_SPLICE:
      reason = "5' splice not canonical or not conserved";
      break;
    case BAD_3_SPLICE:
      reason = "3' splice not canonical or not conserved";
      break;
    case BAD_5_SPLICE_UTR:
      reason = "5' splice not canonical or not conserved (UTR)";
      break;
    case BAD_3_SPLICE_UTR:
      reason = "3' splice not canonical or not conserved (UTR)";
      break;
    case NONSENSE:
      reason = "Nonsense mutation";
      break;
    case FSHIFT:
      {
        if (problem->cds_gap ==  FSHIFT_OK) 
          reason = "Frame-shift gap [gaps not clean]";
        else if (problem->cds_gap == CLN_GAPS) 
          reason = "Frame-shift gap [gaps clean but overlapping/near boundary]";
        else
          reason = "Frame-shift gap";
      }
      break;
    case BAD_INTRON:
      reason = "Adjacent 5' and 3' splice sites don't match (not GT-AG, GC-AG, or AT-AC)";
      break;
    case NO_ALN:
      reason = "Incomplete alignment for cds";
      break;
    case TOO_MANY_Ns:
      reason = "Too many Ns";
      print_alignment = FALSE;
      break;
    case WARN_FSHIFT:
      reason = "Frame shift gap (with compensatory gap)";
      break;
    case WARN_Ns:
      reason = "CDS contains Ns";
      break;
    default:
      die("ERROR in write_log: unknown problem->status %i\n", problem->status);
    }

    fprintf(logf, "%s (%d-%d):\n", reason,
            msa_map_msa_to_seq(map, problem->feat->start), 
            msa_map_msa_to_seq(map, problem->feat->end));

    if (print_alignment)
      dump_aln(logf, problem->feat, msa, problem->status == NONSENSE);
  }
}

/* convert a status_type to a string */
char *status_type_str(status_type status) {
  switch(status) {
  case BAD_REF_START:
    return "bad_ref_start";
  case BAD_REF_STOP:
    return "bad_ref_stop";
  case BAD_REF_5_SPLICE:
    return "bad_ref_5_splice";
  case BAD_REF_3_SPLICE:
    return "bad_ref_3_splice";
  case BAD_REF_ORF:
    return "bad_ref_orf";
  case BAD_REF_INDEL_STRICT_FAIL:
    return "bad_ref_indel_strict";
  case NO_ALN:
    return "no_alignment";
  case OKAY:
    return "okay";
  case BAD_START:
    return "bad_start";
  case BAD_STOP:
    return "bad_stop";
  case BAD_5_SPLICE:
    return "bad_5_splice";
  case BAD_3_SPLICE:
    return "bad_3_splice";
  case BAD_5_SPLICE_UTR:
    return "bad_5_splice_utr";
  case BAD_3_SPLICE_UTR:
    return "bad_3_splice_utr";
  case NONSENSE:
    return "nonsense";
  case FSHIFT:
    return "frameshift";
  case BAD_INTRON:
    return "bad_intron";
  case TOO_MANY_Ns:
    return "too_many_Ns";
  case WARN_FSHIFT:
    return "warning_fshift";
  case WARN_Ns:
    return "warning_Ns";
  default:
    die("ERROR: status_type_str unknown status %i\n", status);
    return "unknown status";
  }
} 

/* write one problem to machine log */
void write_machine_problem(FILE *mlogf, GFF_FeatureGroup *group, 
                           Problem *problem, msa_coord_map *map) {
  char *featName;
  int start, end;
  if (problem->feat == NULL) {
    featName = ((GFF_Feature*)lst_get_ptr(group->features, 0))->seqname->chars;
  } else {
    featName = problem->feat->seqname->chars;
  }
  if (problem->start >= 0) {
    start = problem->start;
    end = problem->end;
  } else {
    start = group->start;
    end = group->end;
  }
  fprintf(mlogf, "%s\t%s\t%d\t%d\t%s\t%d\t%d\n",
          group->name->chars, featName,
          msa_map_msa_to_seq(map, start)-1, 
          msa_map_msa_to_seq(map, end),
          status_type_str(problem->status),
          msa_map_msa_to_seq(map, group->start)-1, 
          msa_map_msa_to_seq(map, group->end));
}

/* write machine-readable log entry for discarded feature */
void write_machine_log(FILE *mlogf, GFF_FeatureGroup *group, List *problems,
                       msa_coord_map *map) {
  int i;
  for (i = 0; i < lst_size(problems); i++) {
    write_machine_problem(mlogf, group, lst_get_ptr(problems, i), map);
  }
}


/* checks to see if reference sequence looks okay wrt a given
   list of features */
int ref_seq_okay(List *features, MSA *msa, int offset3, 
                 int indel_strict, int splice_strict, List *problems) {
  List *signals = NULL;
  char *seq = NULL;
  int seqalloc = 0;
  int idx, retval = TRUE;
  GFF_Feature *feat, *lastfeat_helper = NULL;

  if (indel_strict) {
    signals = lst_new_ptr(10);
    str_split(str_new_charstr(SIGNALS), ",", signals);
  }

  for (idx = 0; idx < lst_size(features); idx++) {
    int i, j, len, has_gaps = 0; 

    feat = lst_get_ptr(features, idx);

    if (seqalloc <= feat->end - feat->start + 2) {
      seqalloc = (feat->end - feat->start) * 2; 
      seq = srealloc(seq, seqalloc * sizeof(char));
    }

    for (i = feat->start - 1, len = 0; i < feat->end; i++) {
      if (ss_get_char_pos(msa, i, 0, 0) != GAP_CHAR)
        seq[len++] = ss_get_char_pos(msa, i, 0, 0);
      else if (!has_gaps) has_gaps = 1;
    }
    seq[len] = '\0';
    if (feat->strand == '-') msa_reverse_compl_seq(seq, len);

    if (str_equals_charstr(feat->feature, GFF_START_TYPE) && strcmp(seq, "ATG") != 0) {
      problem_add(problems, feat, BAD_REF_START, -1, -1);
      retval = FALSE;
    }
    else if (str_equals_charstr(feat->feature, GFF_STOP_TYPE) && 
             (feat->frame != 0 || !is_stop_codon(seq))) {
      problem_add(problems, feat, BAD_REF_STOP, -1, -1);
      retval = FALSE;
    }
    else if (str_starts_with_charstr(feat->feature, SPLICE_5) && 
             !is_valid_5splice(seq, splice_strict)) {
      problem_add(problems, feat, BAD_REF_5_SPLICE, -1, -1);
      retval = FALSE;
    }
    else if (str_starts_with_charstr(feat->feature, SPLICE_3) &&
             !is_valid_3splice(&seq[offset3], splice_strict)) {
      problem_add(problems, feat, BAD_REF_3_SPLICE, -1, -1);
      retval = FALSE;
    }
    else if (str_equals_charstr(feat->feature, GFF_CDS_TYPE)) {
      for (i = (3 - feat->frame) % 3; i <= len - 3; i += 3) {
        if (is_stop_codon(&seq[i])) {
          problem_add(problems, feat, BAD_REF_ORF, -1, -1);
          retval = FALSE;
          break;
        }
      }
    }

    if (indel_strict) {
      int strict_okay = TRUE;
      List *signals = lst_new_ptr(10);
      str_split(str_new_charstr(SIGNALS), ",", signals);

      if (str_in_list(feat->feature, signals)) {
        /* reject any signal feature with gaps in the ref seq, unless they
           appear in a non-critical part of a splice site or in a
           "prestart" feature  */
        if (has_gaps) {          
          if (str_starts_with_charstr(feat->feature, SPLICE_5)) {
            if (ss_get_char_pos(msa, feat->start-1, 0, 0) == GAP_CHAR ||
                ss_get_char_pos(msa, feat->start, 0, 0) == GAP_CHAR)
              strict_okay = FALSE;
          }
          else if (str_starts_with_charstr(feat->feature, SPLICE_3)) {
            if (ss_get_char_pos(msa, feat->end-1, 0, 0) == GAP_CHAR ||
                ss_get_char_pos(msa, feat->end-2, 0, 0) == GAP_CHAR)
              strict_okay = FALSE;
          }
          else if (!str_equals_charstr(feat->feature, "prestart"))
            strict_okay = FALSE;
        }
        /* in addition, if two signals occur consec. with gaps and
           only gaps between them, assume a violation of
           --indel-strict */
        if (lastfeat_helper != NULL && lastfeat_helper->end < feat->start-1) {
          int allgaps = 1;
          for (j = lastfeat_helper->end; allgaps && j < feat->start-1; j++) 
                                /* note indexing: -1+1 for end and -1
                                   for start  */
            if (ss_get_char_pos(msa, j, 0, 0) != GAP_CHAR) allgaps = 0;
          if (allgaps) 
            strict_okay = FALSE;
        }
        lastfeat_helper = feat;
      }
      else lastfeat_helper = NULL;
    
      /* also exclude CDS exons of length less than 6 in indel_strict
         case -- these cause problems in exoniphy training because
         start_codon is adjacent to cds5ss */
      if (str_equals_charstr(feat->feature, GFF_CDS_TYPE) && len <= 6)
        strict_okay = FALSE;

      if (!strict_okay) {
        problem_add(problems, feat, BAD_REF_INDEL_STRICT_FAIL, -1, -1);
        retval = FALSE;
      }
      lst_free_strings(signals);
      lst_free(signals);
    }
  }
  if (seq != NULL) sfree(seq);
  return retval;
}

/* Exclude stop codons from all CDS in a group, as necessary.  Record
   any features that are changed, so they can be changed back before
   data is output */
void exclude_stops(GFF_FeatureGroup *group, List *starts_adjusted, 
                   List *ends_adjusted) {
  int j, k;
  List *stops = lst_new_ptr(1), *gfeatures = group->features;
  GFF_Feature *feat;
  lst_clear(stops); lst_clear(ends_adjusted); lst_clear(starts_adjusted);
  for (j = 0; j < lst_size(gfeatures); j++) { /* first grab all stops.  We 
                                                 expect at most one, but more 
                                                 are possible */
    feat = lst_get_ptr(gfeatures, j);
    if (str_equals_charstr(feat->feature, GFF_STOP_TYPE)) lst_push_ptr(stops, feat);
  }
  for (j = 0; j < lst_size(gfeatures); j++) { /* now look at CDSs */
    feat = lst_get_ptr(gfeatures, j);
    if (str_equals_charstr(feat->feature, GFF_CDS_TYPE)) {
      for (k = 0; k < lst_size(stops); k++) { /* check stops */
        GFF_Feature *stop = lst_get_ptr(stops, k);
        if (feat->strand == '+' && stop->strand == '+' && 
            feat->end == stop->end) {
          feat->end -= 3; 
          lst_push_ptr(ends_adjusted, feat);
        }
        else if (feat->strand == '-' && stop->strand == '-' && 
                 feat->start == stop->start) {
          feat->start += 3; 
          lst_push_ptr(starts_adjusted, feat);
        }
      }
    }
  }
  lst_free(stops);
}

/* Restore cds coords to include stop codons, as necessary */
void restore_stops(GFF_FeatureGroup *group, List *starts_adjusted,
                   List *ends_adjusted) {
  int j;
  if (lst_size(ends_adjusted) == 0 && lst_size(starts_adjusted) == 0)
    return;
  for (j = 0; j < lst_size(group->features); j++) {
    GFF_Feature *feat = lst_get_ptr(group->features, j);
    if (str_equals_charstr(feat->feature, GFF_CDS_TYPE)) {
      if (lst_find_ptr(ends_adjusted, feat) != -1) feat->end += 3;
      else if (lst_find_ptr(starts_adjusted, feat) != -1) feat->start -= 3;
    }
  }
}

int main(int argc, char *argv[]) {

  int check_start = 0, check_stop = 0, check_splice = 0, check_nonsense = 0,
    offset5 = 0, offset3 = 0, opt_idx, i, j, indel_strict = 0, no_output = 0,
    check_alignment = 0, splice_strict = 0;
  int ncons_tested, nkept, nconserved_exons;
  int nce_gap_type[NGAP_TYPES], nconsid[NTYPES], nfail[NTYPES];
  double Nfrac = 0.05;
  char c;
  MSA *msa;
  GFF_Set *gff;
  msa_format_type msa_format = UNKNOWN_FORMAT;
  List *keepers, *problems = lst_new_ptr(10), 
    *ends_adjusted = lst_new_ptr(1), *starts_adjusted = lst_new_ptr(1), 
    *discards=NULL, *intron_splice = lst_new_ptr(10);
  char *rseq_fname = NULL;
  FILE *logf = NULL, *mlogf = NULL, *statsf = NULL, *discardf = NULL;
  cds_gap_type fshift_mode = FSHIFT_BAD;
  char *groupby = "transcript_id";
  msa_coord_map *map;
  int *countNs, *countCDSs;
  FILE *infile;
  char *msa_fname;

  struct option long_opts[] = {
    {"start", 0, 0, 's'},
    {"stop", 0, 0, 't'},
    {"splice", 0, 0, 'l'},
    {"nonsense", 0, 0, 'n'},
    {"fshift", 0, 0, 'f'},
    {"conserved", 0, 0, 'c'},
    {"N-limit", 1, 0, 'N'},
    {"clean-gaps", 0, 0, 'e'},
    {"indel-strict", 0, 0, 'I'},
    {"splice-strict", 0, 0, 'C'},
    {"groupby", 1, 0, 'g'},
    {"msa-format", 1, 0, 'i'},
    {"refseq", 1, 0, 'r'},
    {"offset5", 1, 0, 'o'},
    {"offset3", 1, 0, 'p'},
    {"no-output", 0, 0, 'x'},
    {"discards", 1, 0, 'd'},
    {"log", 1, 0, 'L'},
    {"machine-log", 1, 0, 'M'},
    {"stats", 1, 0, 'S'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "N:i:r:L:M:S:g:d:stlnfceICxh", 
                          long_opts, &opt_idx)) != -1) {
    switch(c) {
    case 's':
      check_alignment = check_start = 1;
      break;
    case 't':
      check_alignment = check_stop = 1;
      break;
    case 'l':
      check_alignment = check_splice = 1;
      break;
    case 'n':
      check_alignment = check_nonsense = 1;
      break;
    case 'f':
      check_alignment = 1;
      fshift_mode = FSHIFT_OK;
      break;
    case 'c':
      check_alignment = check_start = check_stop = check_splice = check_nonsense = 1;
      if (fshift_mode < FSHIFT_OK) fshift_mode = FSHIFT_OK;
      break;
    case 'N':
      Nfrac = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'e':
      check_alignment = 1;
      if (fshift_mode < CLN_GAPS) fshift_mode = CLN_GAPS;
      break;
    case 'I':
      check_alignment = 1;
      fshift_mode = NOVRLP_CLN_GAPS;
      indel_strict = 1;
      break;
    case 'C':
      check_alignment = check_splice = splice_strict = 1;
      break;
    case 'g':
      groupby = optarg;
      break;
    case 'i':
      msa_format = msa_str_to_format(optarg);
      if (msa_format == UNKNOWN_FORMAT) die("Bad alignment format.\n");
      break;
    case 'r':
      rseq_fname = optarg;
      break;
    case 'o':
      offset5 = get_arg_int(optarg);
      break;
    case 'p':
      offset3 = get_arg_int(optarg);
      break;
    case 'L':
      logf = phast_fopen(optarg, "w+");
      break;
    case 'M':
      mlogf = phast_fopen(optarg, "w+");
      break;
    case 'S':
      statsf = phast_fopen(optarg, "w+");
      break;
    case 'd':
      discardf = phast_fopen(optarg, "w+");
      break;
    case 'x':
      no_output = 1;
      break;
    case 'h':
      printf("%s", HELP);
      exit(0);
    case '?':
      die("ERROR: Bad argument.  Try the --help option.\n");
    }
  }

  if (optind + 1 >= argc ) {
    die("ERROR:  Missing required arguments.  Try the --help option.\n");
  }
  
  set_seed(-1);

  gff = gff_read_set(phast_fopen(argv[optind], "r"));
  msa_fname = argv[optind+1];
  infile = phast_fopen(msa_fname, "r");
  if (msa_format == UNKNOWN_FORMAT)
    msa_format = msa_format_for_content(infile, 1);
  if (msa_format == MAF) {
    msa = maf_read(infile, 
                   rseq_fname == NULL ? NULL : phast_fopen(rseq_fname, "r"), 
                   1, NULL, NULL, NULL, -1, TRUE, NULL, NO_STRIP, FALSE); 
  }
  else {
    msa = msa_new_from_file_define_format(infile,
                            msa_format, NULL); 
    if (msa->ss == NULL) 
      ss_from_msas(msa, 1, 1, NULL, NULL, NULL, -1, 0);
  }
  if (!msa->ss->tuple_idx)
    die("ERROR: need ordered tuples\n");
  msa_remove_N_from_alph(msa);  /* for backward compatibility (old SS files) */

  if (msa->idx_offset != 0) {   /* avoids offset problem */
    for (i = 0; i < lst_size(gff->features); i++) {
      GFF_Feature *f = lst_get_ptr(gff->features, i);
      f->start -= msa->idx_offset;
      f->end -= msa->idx_offset;
    }
  }

  /* set up coordinate map; assume GFF is for sequence 1 */
  map = msa_build_coord_map(msa, 1);

  /* convert all features */
  for (i = 0; i < lst_size(gff->features); i++) {
    GFF_Feature *f = lst_get_ptr(gff->features, i);
    int newstart, newend;
 
    if (f->start < 0 || f->end < f->start)
      die("ERROR: bad feature in GFF (start=%d, end=%d).\n",
          f->start, f->end);

    newstart = msa_map_seq_to_msa(map, f->start);
    newend = msa_map_seq_to_msa(map, f->end);

    if (newstart < 0 || newend < newstart)
      die("ERROR: unable to map coordinates for feature (start=%d, end=%d).\n",
          f->start, f->end);

    f->start = newstart;
    f->end = newend;
  }

  gff_group(gff, groupby);	/* do this after coord conversion, or
                               group coords and feature coords
                               will be out of sync */

  keepers = lst_new_ptr(lst_size(gff->features));
  if (discardf != NULL) discards = lst_new_ptr(lst_size(gff->features));

  ncons_tested = nkept = nconserved_exons = 0;
  for (i = 0; i < NTYPES; i++) nconsid[i] = 0;
  for (i = 0; i < NTYPES; i++) nfail[i] = 0;
  for (i = 0; i < NGAP_TYPES; i++) nce_gap_type[i] = 0;  

  countNs = smalloc(msa->nseqs * sizeof(int));
  countCDSs = smalloc(msa->nseqs * sizeof(int));

  for (i = 0; i < lst_size(gff->groups); i++) {
    GFF_FeatureGroup *group = lst_get_ptr(gff->groups, i);
    List *gfeatures = group->features;
    GFF_Feature *feat;
    status_type status = OKAY;
    cds_gap_type gt = FSHIFT_BAD;
    problems_clear(problems);

    /* make sure have frame info for CDSs */
    for (j = 0; j < lst_size(gfeatures); j++) {
      feat = lst_get_ptr(gfeatures, j);
      if (str_equals_charstr(feat->feature, GFF_CDS_TYPE) && 
          feat->frame == GFF_NULL_FRAME)
        die("ERROR: Missing frame info for CDS.\n");
    }

    /* First, exclude stop codons from cds's, if necessary (simplifies
       the detection of nonsense mutations). */
    exclude_stops(group, starts_adjusted, ends_adjusted);

    /* In all cases, discard any group for which the reference sequence
       doesn't have valid splice sites or start/stop codons, or has a
       premature stop codon */
    if (!ref_seq_okay(gfeatures, msa, offset3, indel_strict, splice_strict,
                      problems)) {
      status = BAD_REF;
      nfail[BAD_REF]++;
    }
    else
      /* Everything else counts as a potentially valid group */
      ncons_tested++;

    if (status == OKAY && check_alignment) {      
                                /* only bother with below if
                                   interested in cross-species
                                   conservation */

      /* Check first to make sure there's alignment across species in
         the cds; if not, there's no need to look at individual
         features. */
      for (j = 0; j < lst_size(gfeatures); j++) { 
        feat = lst_get_ptr(gfeatures, j);
        if (str_equals_charstr(feat->feature, GFF_CDS_TYPE) &&
            is_incomplete_alignment(feat, msa)) {
          status = NO_ALN;
          nfail[NO_ALN]++;
          problem_add(problems, feat, NO_ALN, -1, -1);
          break;
        }
      }

      if (status == OKAY) {     /* we have alignment and agreement
                                   with the ref seq; now check feature
                                   by feature  */

        lst_clear(intron_splice);
        for (j = 0; j < msa->nseqs; j++) countNs[j] = countCDSs[j] = 0;

        for (j = 0; j < lst_size(gfeatures); j++) {
          feat = lst_get_ptr(gfeatures, j);

          if (feat->end - 1 >= msa->length) 
            die("ERROR: feature extends beyond alignment (%d >= %d).\n",
                feat->end - 1, msa->length);
        
          if (check_start && str_equals_charstr(feat->feature, GFF_START_TYPE)) {

            nconsid[BAD_START]++;

            if (!is_conserved_start(feat, msa)) {
              status = BAD_START;
              problem_add(problems, feat, BAD_START, -1, -1);
            }
          }

          else if (check_stop && str_equals_charstr(feat->feature, GFF_STOP_TYPE)) {

            nconsid[BAD_STOP]++;

            if (!is_conserved_stop(feat, msa)) {
              status = BAD_STOP;
              problem_add(problems, feat, BAD_STOP, -1, -1);
            }
          }

          else if (check_splice && 
                   str_equals_charstr(feat->feature, SPLICE_5)) {

            nconsid[BAD_5_SPLICE]++;

            if (!is_conserved_5splice(feat, msa, offset5, splice_strict)) {
              status = BAD_5_SPLICE;
              problem_add(problems, feat, BAD_5_SPLICE, -1, -1);
            }
            else lst_push_ptr(intron_splice, feat);
          }

          else if (check_splice && 
                   str_equals_charstr(feat->feature, SPLICE_5_UTR)) {

            nconsid[BAD_5_SPLICE_UTR]++;

            if (!is_conserved_5splice(feat, msa, offset5, splice_strict)) {
              status = BAD_5_SPLICE_UTR;
              problem_add(problems, feat, BAD_5_SPLICE_UTR, -1, -1);
            }
            else lst_push_ptr(intron_splice, feat);
          }

          else if (check_splice && str_equals_charstr(feat->feature, SPLICE_3)) {


            nconsid[BAD_3_SPLICE]++;

            if (!is_conserved_3splice(feat, msa, offset3, splice_strict)) {
              status = BAD_3_SPLICE;
              problem_add(problems, feat, BAD_3_SPLICE, -1, -1);
            }
            else lst_push_ptr(intron_splice, feat);
          }

          else if (check_splice && str_equals_charstr(feat->feature, SPLICE_3)) {

            nconsid[BAD_3_SPLICE_UTR]++;

            if (!is_conserved_3splice(feat, msa, offset3, splice_strict)) {
              status = BAD_3_SPLICE_UTR;
              problem_add(problems, feat, BAD_3_SPLICE_UTR, -1, -1);
            }
            else lst_push_ptr(intron_splice, feat);
          }

          else if (str_equals_charstr(feat->feature, GFF_CDS_TYPE)) {
 
            if (fshift_mode > FSHIFT_BAD 
		&& (gt = get_cds_gap_type(feat, msa, problems)) < fshift_mode) {
              if (status == OKAY || status == NONSENSE) status = FSHIFT;
            }

            if (check_nonsense && !is_nonsense_clean(feat, msa, problems)) {
              if (status == OKAY) status = NONSENSE;
            }

            if (Nfrac < 1) 
              get_N_counts(countNs, countCDSs, feat, msa);
          }
        } /* end loop through features in group */

        /* still have to make sure splice sites are paired correctly
           (GT-AG, GC-AG, AT-AC) */
        if (status == OKAY && !splice_strict && lst_size(intron_splice) >= 2 &&
            !are_introns_okay(intron_splice, msa, problems, offset5, offset3)) 
          status = BAD_INTRON;

        /* also check fraction of Ns */
        if (Nfrac < 1) {
          enum {MY_OKAY, MY_FAIL, MY_WARN} Nstatus = MY_OKAY;
          for (j = 0; j < msa->nseqs; j++) {
            if ((double)countNs[j] / countCDSs[j] > Nfrac) Nstatus = MY_FAIL;
            if (Nstatus == MY_OKAY && countNs[j] > 0) Nstatus = MY_WARN;
          }
          if (Nstatus == MY_FAIL) {
            problem_add(problems, NULL, TOO_MANY_Ns, -1, -1);
            if (status == OKAY) status = TOO_MANY_Ns;
          }
          else if (Nstatus == MY_WARN) 
            problem_add(problems, NULL, WARN_Ns, -1, -1);
        }

        /* if collecting stats, record counts for failures */
        if (statsf != NULL) {
          if (status != OKAY) {
            for (j = 0; j < lst_size(problems); j++) {
              struct Problem *problem = lst_get_ptr(problems, j);
              status_type ftype = problem->status;
              if ((ftype == FSHIFT || ftype == NONSENSE) && 
                  status != FSHIFT && status != NONSENSE)
                continue;       /* don't count secondary frame shifts
                                   and nonsense mutations */ 

              if (ftype == BAD_INTRON && j % 2 == 0)
                continue;       /* only count one of every pair of these */

              nfail[ftype]++;
            }
          }

          /* also keep track of the total number of "conserved exons", and
             the number having each kind of gap */
          if ((status == OKAY || (status == FSHIFT && gt >= FSHIFT_OK))) {
            nconserved_exons++;
            nce_gap_type[gt]++;     /* number of conserved exons having
                                       given type of gaps */
          }
        }
      } /* end if (status == OKAY) [checks for conserved features] */
    } /* end if (status == OKAY && check_alignment) [all cross-species
         checks] */

    /* now we have looked at the whole group; we just need to do some
       final accounting and logging */

    if (status == OKAY) {
      nkept++;
      if (!no_output) {
        restore_stops(group, starts_adjusted, ends_adjusted);
        for (j = 0; j < lst_size(gfeatures); j++)
          lst_push_ptr(keepers, lst_get_ptr(gfeatures, j));
      }
      if (logf != NULL && lst_size(problems) > 0) /* warnings only */
        write_log(logf, group, status, problems, msa, map);
      if (mlogf != NULL) {
        /* no problem, need to add an okay status to log */
        problem_add(problems, NULL, OKAY, -1, -1);
        write_machine_log(mlogf, group, problems, map); /* may include
                                                           warnings */
      }
    }
    else {
      if (discardf != NULL) {
        restore_stops(group, starts_adjusted, ends_adjusted);
        for (j = 0; j < lst_size(gfeatures); j++) 
          lst_push_ptr(discards, lst_get_ptr(gfeatures, j));
      }
      if (logf != NULL) 
        write_log(logf, group, status, problems, msa, map);
      if (mlogf != NULL)
        write_machine_log(mlogf, group, problems, map);
    }
  } /* end loop over groups */

  /* write main output and discards */
  if (!no_output || discardf != NULL) {
    /* first map features back to coord frame of reference seq. */
    for (i = 0; i < lst_size(gff->features); i++) {
      GFF_Feature *f = lst_get_ptr(gff->features, i);
      f->start = msa_map_msa_to_seq(map, f->start) + msa->idx_offset;
      f->end = msa_map_msa_to_seq(map, f->end) + msa->idx_offset;
    }

    if (!no_output) {
      gff->features = keepers;
      gff_print_set(stdout, gff);
    }

    if (discardf != NULL) {
      gff->features = discards;
      gff_print_set(discardf, gff);
    }
  }


  /* dump counts to stats file */
  if (statsf != NULL) {
    fprintf(statsf, "#%11s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n", 
            "total", "nbad_ref", "nconsid", "nkept", "nno_aln", 
            "nbad_starts", "(out of)", "nbad_stops", "(out of)", 
            "nbad_5spl", "(out of)", "nbad_3spl", "(out of)", 
            "nbad_5utr", "(out of)", "nbad_3utr", "(out of)", 
            "nbad_intron", "nnons", "nfshifts", "nNs", "ncons_exons", 
            "nce_ngaps", "nce_nov_cln", "nce_clean", "nce_fshftok");
    fprintf(statsf, "%12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d\n", 
            nfail[BAD_REF]+ncons_tested, nfail[BAD_REF], ncons_tested, nkept, 
            nfail[NO_ALN], nfail[BAD_START], nconsid[BAD_START], 
            nfail[BAD_STOP], nconsid[BAD_STOP], nfail[BAD_5_SPLICE], 
            nconsid[BAD_5_SPLICE], nfail[BAD_3_SPLICE], nconsid[BAD_3_SPLICE],
            nfail[BAD_5_SPLICE_UTR], nconsid[BAD_5_SPLICE_UTR],
            nfail[BAD_3_SPLICE_UTR], nconsid[BAD_3_SPLICE_UTR], 
            nfail[BAD_INTRON], nfail[NONSENSE], nfail[FSHIFT], 
            nfail[TOO_MANY_Ns], nconserved_exons, nce_gap_type[NGAPS], 
            nce_gap_type[NOVRLP_CLN_GAPS], nce_gap_type[CLN_GAPS], 
            nce_gap_type[FSHIFT_OK]);
    fprintf(statsf, "%s", STATS_DESCRIPTION);
  }

  if (logf != NULL) phast_fclose(logf);
  if (mlogf != NULL) phast_fclose(mlogf);
  if (statsf != NULL) phast_fclose(statsf);
  if (discardf != NULL) phast_fclose(discardf);

  return 0;
}


