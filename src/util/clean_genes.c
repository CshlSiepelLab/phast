/* $Id: clean_genes.c,v 1.4 2004-06-05 06:28:55 acs Exp $
   Written by Adam Siepel, 2003-2004
   Copyright 2003-2004, Adam Siepel, University of California */

#include <stdlib.h>
#include <stdio.h>
#include <msa.h>
#include <gff.h>
#include <sufficient_stats.h>
#include <getopt.h>
#include <maf.h>

/* types of features examined */
#define START "start"
#define STOP "stop"
#define SPLICE_5 "5'splice"
#define SPLICE_3 "3'splice"
#define SPLICE_5_UTR "5'splice_utr"
#define SPLICE_3_UTR "3'splice_utr"
#define CDS "cds"
/* NOTE: SPLICE_5_UTR and SPLICE_3_UTR must have SPLICE_5 and SPLICE_3
   as prefixes (respectively).  This property is used to simplify some
   of the case handling below */ 

/* possible values for status of each group of features */
typedef enum {OKAY, BAD_REF, NO_ALN, BAD_START, BAD_STOP, BAD_5_SPLICE, 
              BAD_3_SPLICE, BAD_5_SPLICE_UTR, BAD_3_SPLICE_UTR, 
              NONSENSE, FSHIFT, BAD_INTRON, NTYPES} 
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
#define MIN_BLOCK_SIZE 30

/* "signal" feature types, used for 'indel-strict' mode (special-case
   stuff related to exoniphy training) */
#define SIGNALS "start,stop,5'splice,3'splice,cds5'ss,cds3'ss,prestart"

void print_usage() {
    printf("\n\
PROGRAM:        clean_genes\n\
\n\
DESCRIPTION:    Given a GFF describing a set of genes and a corresponding \n\
                multiple alignment, output a new GFF with only those\n\
                genes that meet certain \"cleanliness\" criteria. The\n\
                coordinates in the GFF are assumed to correspond to\n\
                the reference sequence in the alignment, which is\n\
                assumed to be the first one listed.  Default behavior\n\
                is simply to require that all annotated start/stop codons and\n\
                splice sites are valid in the reference sequence (GT-AG, \n\
                GC-AG, and AT-AC splice sites are allowed).  This can\n\
                be used with an \"alignment\" consisting of a single\n\
                sequence to filter out incorrect annotations.  Options\n\
                are available to impose additional criteria as well,\n\
                mostly having to do with conservation across species\n\
                (see the '--conserved' option in particular).\n\
\n\
USAGE:          clean_genes [options] <gff_fname> <msa_fname>\n\
\n\
OPTIONS:        \n\
\n\
    --start, -s\n\
        Require conserved start codons (all species)\n\
\n\
    --stop, -t\n\
        Require conserved stop codons (all species)\n\
\n\
    --splice, -l    \n\
        Require conserved splice sites (all species).  By default,\n\
        only GT-AG, GC-AG, and AT-AC splice sites are allowed (see also\n\
        --splice-strict)\n\
\n\
    --fshift, -h     \n\
        Require that no frame-shift gap is present in any species.  Frame \n\
        shifts are evaluated with respect to the reference sequence.  Gaps \n\
        that have non-multiple-of-three lengths are allowed if  \n\
        compensatory gaps occur nearby (see source code for details).\n\
\n\
    --nonsense, -n  \n\
        Require that no premature stop codon is present in any species.\n\
\n\
    --conserved, -c\n\
        Implies --start, --stop, --splice, --fshift, and --nonsense.  
        Recommended option for cross-species analysis.\n\
\n\
    --clean-gaps, -e\n\
        Require all cds gaps to be multiples of three in length.  Can be \n\
        used with --conserved.\n\
\n\
    --indel-strict, -I\n\
        Implies --clean_gaps, usually used with --conserved.  Prohibits\n\
        overlapping cds gaps in different sequences, gaps near cds \n\
        boundaries, and gaps in the reference sequence within and between\n\
        flanking features (splice sites, etc.; see code for details).\n\
        Designed for use in training a phylo-HMM with an indel model.\n\
\n\
    --splice-strict, -C\n\
        Implies --splice.  Allow only GT-AG canonical splice sites.  Useful\n\
        when training a gene finder with a simple model for splice sites.\n\
\n\
    --groupby, -g <tag>\n\
        Group features according to specified tag (default\n\
        \"transcript_id\").  If any feature within a group fails, the\n\
        entire group will be discarded.  By choosing to group features\n\
        according to different criteria, you can make the program\n\
        \"clean\" the data set at different levels.  For example, to\n\
        clean at the level of individual exons, add a tag like\n\
        \"exon_id\" to indicate exons (see the program \"refeature\"),\n\
        and then invoke clean_genes with \"--groupby exon_id\".\n\
\n\
    --msa-format, -i PHYLIP|FASTA|PSU|SS|LAV|MAF\n\
        (Default PHYLIP) Alignment file format. \n\
\n\
    --refseq, -r <seqfile.fa>\n\
        (Required with --msa-format MAF)  Complete reference \n\
        sequence for alignment (FASTA format).\n\
\n\
    --offset5, -o <n>\n\
        (Default 0)  Offset of canonical \"GT\" with respect to boundary \n\
        on *intron side* of annotated 5' splice sites.  Useful with\n\
        annotations that describe a window around the canonical splice site.\n\
\n\
    --offset3, -p <n>\n\
        (Default 0)  Offset of canonical \"AG\" with respect to boundary \n\
        on intron side of annotated 3' splice sites.\n\
\n\
    --log, -L <fname>\n\
        Write human-readable log to specified file.\n\
\n\
    --machine-log, -M <fname>\n\
        Like --log, but produces more concise, machine-readable log.\n\
\n\
    --stats, -S <fname>\n\
        Write statistics on retained and discarded features to specified file.\n\
\n\
    --discards, -d <fname>\n\
        Write discarded features to specified file.\n\
\n\
    --no-output, -x\n\
        Suppress output of \"cleaned\" features to stdout.  Useful if only\n\
        log file and/or stats are of interest.\n\
\n\
    --help, -h\n\
        Print this help message.\n\
\n\
NOTES:  Currently, feature types are defined at compile time as follows.\n\
\n\
               coding exon    <-> \"%s\"\n\
               start codon    <-> \"%s\"\n\
               stop codon     <-> \"%s\"\n\
               5' splice site <-> \"%s\"\n\
               3' splice site <-> \"%s\"\n\
\n\
        In addition, splice sites in UTR can be separately designated as\n\
        \"%s\" and \"%s\".  Errors in these sites will be given a different\n\
        code in the log files, which can be useful for tracking purposes.\n\
\n\
        If evaluation is done at the level of individual exons (see\n\
        --groupby), then splice sites are considered independently\n\
        rather than in the context of introns.  Thus, the exons flanking\n\
        a GT-AC or AT-AG intron might be considered \"clean\".\n\
\n\
        With --fshift and --nonsense, it is possible for entries\n\
        to pass through that have stop codons in the frame of the\n\
        *reference* sequence, although they do not have any in their\n\
        own frame.  Use --clean-gaps instead to guarantee that no stop\n\
        codons occur in any sequence in the frame of the reference\n\
        sequence.\n\n", 
           CDS, START, STOP, SPLICE_5, SPLICE_3, SPLICE_5_UTR, SPLICE_3_UTR);
}

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
# nbad_intron   number splice sites incorrectly paired at intron\n\
# nnons         number failing nonsense mutation test (not counting ones that fail above tests)\n\
# nfshifts      number failing frameshift test (not counting ones that fail above tests)\n\
# *** ALIGNMENT GAP STATS ***\n\
# ncons_exons   number of conserved exons (same as nkept if --conserved)\n\
# nce_ngaps     number with no alignment gaps\n\
# nce_nov_cln   number with only nonoverlapping \"clean\" gaps (multiple of three lengths)\n\
# nce_clean     number with clean gaps that overlap\n\
# nce_fshftok   number with compensatory frame-shifting gaps as allowed by --fshift\n";

inline int is_conserved_start(GFF_Feature *feat, MSA *msa) {
  char tuplestr[4];
  int j;
  int start = feat->start - 1;  /* base 0 indexing */
  tuplestr[3] = '\0';
  for (j = 0; j < msa->nseqs; j++) {
    tuplestr[0] = ss_get_char_tuple(msa, msa->ss->tuple_idx[start], j, 0);
    tuplestr[1] = ss_get_char_tuple(msa, msa->ss->tuple_idx[start+1], j, 0);
    tuplestr[2] = ss_get_char_tuple(msa, msa->ss->tuple_idx[start+2], j, 0);
    if (strcmp(tuplestr, "NNN") == 0) continue;
    if (feat->strand == '-') msa_reverse_compl_seq(tuplestr, 3);
    if (strcmp(tuplestr, "ATG") != 0) return 0;
  }
  return 1;
}

inline int is_stop_codon(char *str) {
 return (strncmp(str, "TAA", 3) == 0 || strncmp(str, "TAG", 3) == 0 ||
         strncmp(str, "TGA", 3) == 0);
 /* strncmp allows testing of codons in the middle of larger strings */
}

inline int is_conserved_stop(GFF_Feature *feat, MSA *msa) {
  char tuplestr[4];
  int j;
  int start = feat->start - 1;
  tuplestr[3] = '\0';
  for (j = 0; j < msa->nseqs; j++) {
    tuplestr[0] = ss_get_char_tuple(msa, msa->ss->tuple_idx[start], j, 0);
    tuplestr[1] = ss_get_char_tuple(msa, msa->ss->tuple_idx[start+1], j, 0);
    tuplestr[2] = ss_get_char_tuple(msa, msa->ss->tuple_idx[start+2], j, 0);
    if (strcmp(tuplestr, "NNN") == 0) continue;
    if (feat->strand == '-') msa_reverse_compl_seq(tuplestr, 3);
    if (!is_stop_codon(tuplestr)) return 0;
  }
  return 1;
}

/* returns 1 if GT, GC or AT, and 0 otherwise.  If splice_strict,
   returns 0 unless GT. */
inline int is_valid_5splice(char *str, int splice_strict) {
  if (strncmp(str, "GT", 2) == 0) return 1;
  else if (splice_strict) return 0;
  else if (strncmp(str, "GC", 2) == 0 || strncmp(str, "AT", 2) == 0)
    return 1;
  return 0;
}

inline int is_conserved_5splice(GFF_Feature *feat, MSA *msa, int offset5,
                                int splice_strict) {
  char tuplestr[3];
  int j, start = feat->start - 1;  /* base 0 indexing */
  if (feat->strand == '-') start += offset5;
  tuplestr[2] = '\0';
  for (j = 0; j < msa->nseqs; j++) {
    tuplestr[0] = ss_get_char_tuple(msa, msa->ss->tuple_idx[start], j, 0);
    tuplestr[1] = ss_get_char_tuple(msa, msa->ss->tuple_idx[start+1], j, 0);
    if (strcmp(tuplestr, "NN") == 0) continue;
    if (feat->strand == '-') msa_reverse_compl_seq(tuplestr, 2);
    if (!is_valid_5splice(tuplestr, splice_strict)) return 0;
  }
  return 1;
}

/* returns 1 for AG or AC, 0 otherwise.  If splice_strict, returns 0
   unless AG. */
inline int is_valid_3splice(char *str, int splice_strict) {
  if (strncmp(str, "AG", 2) == 0) return 1;
  else if (splice_strict) return 0;
  else if (strncmp(str, "AC", 2) == 0) return 1;
  return 0;
}

inline int is_conserved_3splice(GFF_Feature *feat, MSA *msa, int offset3,
                                int splice_strict) {
  char tuplestr[3];
  int j, start = feat->start - 1;  /* base 0 indexing */
  if (feat->strand == '+') start += offset3;
  tuplestr[2] = '\0';
  for (j = 0; j < msa->nseqs; j++) {
    tuplestr[0] = ss_get_char_tuple(msa, msa->ss->tuple_idx[start], j, 0);
    tuplestr[1] = ss_get_char_tuple(msa, msa->ss->tuple_idx[start+1], j, 0);
    if (strcmp(tuplestr, "NN") == 0) continue;
    if (feat->strand == '-') msa_reverse_compl_seq(tuplestr, 2);
    if (!is_valid_3splice(tuplestr, splice_strict)) return 0;
  }
  return 1;
}

int is_nonsense_clean(GFF_Feature *feat, MSA *msa) {
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
      if (is_stop_codon(&seq[i])) return 0;
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
int scan_for_gaps(GFF_Feature *feat, MSA *msa) {
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

  for (i = 0; i < lst_size(gaps); i++) free(lst_get_ptr(gaps, i));
  lst_free(gaps);
  return retval;
}

/* look for frame-shift gaps using a slightly more sophisticated
   algorithm, which allows for compensatory indels.  The strategy here
   is to identify maximal gapless blocks of greater than
   MIN_BLOCK_SIZE sites, then to make sure that in the gappy portions
   between them, each sequence has a total number of gaps that equals
   the total number for the reference sequence, modulo 3.  Returns 1 if
   all gaps look okay (no net frame shift) and 0 otherwise. */
int is_fshift_okay(GFF_Feature *feat, MSA *msa) {
  int *ngaps = smalloc(msa->nseqs * sizeof(int));
  int i, j, blk_beg, blk_end;

  for (j = 0; j < msa->nseqs; j++) ngaps[j] = 0;
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

    if (blk_end - blk_beg >= MIN_BLOCK_SIZE ||
        blk_beg == feat->end || /* gaps at end of aln */
        blk_end == feat->end) { /* short block at end of aln */
      /* check total number of gaps since last retained block or
         beginning of alignment; must be same as reference sequence,
         mod 3 */
      for (j = 0; j < msa->nseqs && ngaps[j] % 3 == ngaps[0] % 3; j++);
      if (j != msa->nseqs) {   /* reject alignment */
        free(ngaps);
        return 0;
      }

      /* reset ngaps (note: done only if block exceeds size
         threshold) */
      for (j = 0; j < msa->nseqs; j++) ngaps[j] = 0;
    }
  }
  free(ngaps);
  return 1;
}

/* Given a cds feature, determine whether it has no gaps (NGAPS),
   "clean" gaps (all multiples of 3 in length; CLEAN_GAPS)
   non-overlapping clean gaps (NOVRLP_CLN_GAPS), "okay" gaps (only
   temporary frame shifts, corrected by compensatory indels;
   FSHIFT_OK), or real frame-shift gaps (FSHIFT_BAD) */
cds_gap_type get_cds_gap_type(GFF_Feature *feat, MSA *msa) {
  cds_gap_type retval = scan_for_gaps(feat, msa);

  if (retval == FSHIFT_BAD && is_fshift_okay(feat, msa))
    retval = FSHIFT_OK;
                                /* most of the time the call to
                                   is_fshift_okay won't be
                                   necessary */
  return retval;
}

/* returns 1 if any row of the alignment consists completely of gaps in
   the region of the specified feature, and returns 0 otherwise.
   For use in testing whether a cds is present in the alignment for all
   species */
inline int is_incomplete_alignment(GFF_Feature *feat, MSA *msa) {
  int i, j;
  for (i = 1; i < msa->nseqs; i++) { /* no need to check reference seq */
    int row_all_gaps = 1;
    for (j = feat->start - 1; row_all_gaps && j < feat->end; j++)
      if (ss_get_char_tuple(msa, msa->ss->tuple_idx[j], i, 0) != GAP_CHAR)
        row_all_gaps = 0;
    if (row_all_gaps) return 1;
  }
  return 0;
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

/* returns 1 if a 5' and a 3' splice site form a valid pair */
inline int is_valid_splice_pair(char *ss5, char *ss3) {
  if ((strcmp(ss5, "GT") == 0 && strcmp(ss3, "AG") == 0) ||
      (strcmp(ss5, "GC") == 0 && strcmp(ss3, "AG") == 0) ||
      (strcmp(ss5, "AT") == 0 && strcmp(ss3, "AC") == 0)) 
    return 1;
  return 0;
}

/* given a list of 5' and 3' splice sites extracted from a group,
   check whether they form valid pairs in all species */
int are_introns_okay(List *intron_splice,  MSA *msa, GFF_Feature **bad_feat,
                     int offset5, int offset3) {
  int i, j, start1, start2;
  char str1[3], str2[3];
  char strand;
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
        if (strcmp(str1, "NN") == 0 || strcmp(str2, "NN") == 0) continue;
        if (strand == '-') {
          msa_reverse_compl_seq(str1, 2);
          msa_reverse_compl_seq(str2, 2);
        }
        if (!is_valid_splice_pair(str1, str2)) {
          *bad_feat = f2;
          return 0;
        }
      }
      i++;                      /* no need to look at next one */
    }
  }
  return 0;
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
               cds_gap_type gt, List *badfeats, List *failure_types, MSA *msa,
               msa_coord_map *map) {

  int i;

  /* special cases: we don't have info on individual features */
  if (status == BAD_REF) {
    fprintf(logf, "****\nSkipped '%s' -- features not consistent with reference sequence.\n", group->name->chars);
    return;
  }
  else if (status == NO_ALN) {
    fprintf(logf, "****\nSkipped '%s' -- lack of alignment across species.\n", 
            group->name->chars);
    return;
  }

  fprintf(logf, "****\nDiscarded '%s'\n", group->name->chars);

  for (i = 0; i < lst_size(failure_types); i++) {
    status_type ftype = lst_get_int(failure_types, i);
    GFF_Feature *badfeat = lst_get_ptr(badfeats, i);
    char *reason;

    switch (ftype) {
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
        if (gt ==  FSHIFT_OK) 
          reason = "Frame-shift gap [gaps not clean]";
        else if (gt == CLN_GAPS) 
          reason = "Frame-shift gap [gaps clean but overlapping]";
        else
          reason = "Frame-shift gap";
      }
      break;
    case BAD_INTRON:
      reason = "3' splice site doesn't pair correctly with 5' splice site (not GT-AG, GC-AG, or AT-AC)";
      break;
    case NO_ALN:
      reason = "Incomplete alignment for cds";
      break;
    default: 
      assert(0);
    }

    fprintf(logf, "%s (%d-%d):\n", reason, msa_map_msa_to_seq(map, badfeat->start), 
            msa_map_msa_to_seq(map, badfeat->end));

    dump_aln(logf, badfeat, msa, ftype == NONSENSE);
  }
}

/* get range of group; used for machine log */
void get_range(GFF_FeatureGroup *group, int *start, int *end,
               String **seqname) {
  int i;
  *start = INFTY; *end = -1; *seqname = NULL;
  for (i = 0; i < lst_size(group->features); i++) {
    GFF_Feature *f = lst_get_ptr(group->features, i);
    if (f->start < *start) *start = f->start;
    if (f->end > *end) *end = f->end;
    if (*seqname == NULL) *seqname = f->seqname;
    else if (!str_equals(*seqname, f->seqname)) 
      die("ERROR: Members of group '%s' have inconsistent sequence names (chromosomes).\n");
  }
}

/* write machine-readable log entry for discarded feature */
void write_machine_log(FILE *mlogf, GFF_FeatureGroup *group, status_type status, 
                       cds_gap_type gt, List *badfeats, List *failure_types, 
                       msa_coord_map *map) {

  int i;
  char *reason;

  /* special cases: no info on individual features */
  if (status == BAD_REF || status == NO_ALN || status == OKAY) {
    int s, e;
    String *seqname;
    get_range(group, &s, &e, &seqname);

    switch(status) {
    case BAD_REF:
      reason = "bad_ref";
      break;
    case NO_ALN:
      reason = "no_alignment";
      break;
    case OKAY:
      reason = "okay";
      break;
    default:
      assert(0);
    }

    fprintf(mlogf, "%s\t%s\t%d\t%d\t%s\n", group->name->chars, 
            seqname->chars, msa_map_msa_to_seq(map, s), 
            msa_map_msa_to_seq(map, e), reason);

    return;
  }

  for (i = 0; i < lst_size(failure_types); i++) {
    status_type ftype = lst_get_int(failure_types, i);
    GFF_Feature *badfeat = lst_get_ptr(badfeats, i);

    switch (ftype) {
    case BAD_START:
      reason = "bad_start";
      break;
    case BAD_STOP:
      reason = "bad_stop";
      break;
    case BAD_5_SPLICE:
      reason = "bad_5_splice";
      break;
    case BAD_3_SPLICE:
      reason = "bad_3_splice";
      break;
    case BAD_5_SPLICE_UTR:
      reason = "bad_5_splice_utr";
      break;
    case BAD_3_SPLICE_UTR:
      reason = "bad_3_splice_utr";
      break;
    case NONSENSE:
      reason = "nonsense";
      break;
    case FSHIFT:
      reason = "frameshift";
      break;
    case BAD_INTRON:
      reason = "bad_intron";
      break;
    case NO_ALN:
      reason = "no_alignment";
      break;
    default: 
      assert(0);
    }

    fprintf(mlogf, "%s\t%s\t%d\t%d\t%s\n", group->name->chars, 
            badfeat->seqname->chars,
            msa_map_msa_to_seq(map, badfeat->start), 
            msa_map_msa_to_seq(map, badfeat->end),
            reason);
  }
}

/* checks to see if reference sequence looks okay wrt a given
   list of features */
int ref_seq_okay(List *features, MSA *msa, int offset3, int indel_strict,
                 int splice_strict) {
  static List *signals = NULL;
  static char *seq = NULL;
  static int seqalloc = 0;
  int idx;
  GFF_Feature *feat, *lastfeat_helper = NULL;

  if (indel_strict && signals == NULL) {
    signals = lst_new_ptr(10);
    str_split(str_new_charstr(SIGNALS), ",", signals);
  }

  for (idx = 0; idx < lst_size(features); idx++) {
    int i, j, len, has_gaps = 0; 

    feat = lst_get_ptr(features, idx);

    if (seqalloc <= feat->end - feat->start + 2) {
      seqalloc = (feat->end - feat->start) * 2; 
      seq = realloc(seq, seqalloc * sizeof(char));
    }

    for (i = feat->start - 1, len = 0; i < feat->end; i++) {
      if (ss_get_char_pos(msa, i, 0, 0) != GAP_CHAR)
        seq[len++] = ss_get_char_pos(msa, i, 0, 0);
      else if (!has_gaps) has_gaps = 1;
    }
    seq[len] = '\0';
    if (feat->strand == '-') msa_reverse_compl_seq(seq, len);

    if (str_equals_charstr(feat->feature, START) && strcmp(seq, "ATG") != 0) 
      return 0;
    else if (str_equals_charstr(feat->feature, STOP) && !is_stop_codon(seq)) 
      return 0;
    else if (str_starts_with_charstr(feat->feature, SPLICE_5) && 
             !is_valid_5splice(seq, splice_strict))
      return 0;
    else if (str_starts_with_charstr(feat->feature, SPLICE_3) &&
             !is_valid_3splice(seq, splice_strict))
      return 0;
    else if (str_equals_charstr(feat->feature, CDS)) {
      for (i = (3 - feat->frame) % 3; i <= len - 3; i += 3) 
        if (is_stop_codon(&seq[i])) return 0;
    }

    if (indel_strict) {
      if (str_in_list(feat->feature, signals)) {
        /* reject any signal feature with gaps in the ref seq, unless they
           appear in a non-critical part of a splice site or in a
           "prestart" feature  */
        if (has_gaps) {          
          if (str_starts_with_charstr(feat->feature, SPLICE_5)) {
            if (ss_get_char_pos(msa, feat->start-1, 0, 0) == GAP_CHAR ||
                ss_get_char_pos(msa, feat->start, 0, 0) == GAP_CHAR)
              return 0;
          }
          else if (str_starts_with_charstr(feat->feature, SPLICE_3)) {
            if (ss_get_char_pos(msa, feat->end-1, 0, 0) == GAP_CHAR ||
                ss_get_char_pos(msa, feat->end-2, 0, 0) == GAP_CHAR)
              return 0;
          }
          else if (!str_equals_charstr(feat->feature, "prestart"))
            return 0;
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
          if (allgaps) return 0;
        }
        lastfeat_helper = feat;
      }
      else lastfeat_helper = NULL;
    }
  }
  return 1;
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
    if (str_equals_charstr(feat->feature, STOP)) lst_push_ptr(stops, feat);
  }
  for (j = 0; j < lst_size(gfeatures); j++) { /* now look at CDSs */
    feat = lst_get_ptr(gfeatures, j);
    if (str_equals_charstr(feat->feature, CDS)) {
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
    if (str_equals_charstr(feat->feature, CDS)) {
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
  char c;
  MSA *msa;
  GFF_Set *gff;
  msa_format_type msa_format = PHYLIP;
  List *keepers, *badfeats = lst_new_ptr(10), *failure_types = lst_new_int(10),
    *ends_adjusted = lst_new_ptr(1), *starts_adjusted = lst_new_ptr(1), 
    *discards, *intron_splice = lst_new_ptr(10);
  char *rseq_fname = NULL;
  FILE *logf = NULL, *mlogf = NULL, *statsf = NULL, *discardf = NULL;
  cds_gap_type fshift_mode = FSHIFT_BAD;
  String *groupby = str_new_charstr("transcript_id");
  msa_coord_map *map;

  struct option long_opts[] = {
    {"start", 0, 0, 's'},
    {"stop", 0, 0, 't'},
    {"splice", 0, 0, 'l'},
    {"nonsense", 0, 0, 'n'},
    {"fshift", 0, 0, 'f'},
    {"conserved", 0, 0, 'c'},
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

  while ((c = getopt_long(argc, argv, "i:r:L:M:S:g:d:stlnfceICxh", 
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
      groupby = str_new_charstr(optarg);
      break;
    case 'i':
      if (!strcmp(optarg, "PSU")) msa_format = PSU;
      else if (!strcmp(optarg, "FASTA")) msa_format = FASTA;
      else if (!strcmp(optarg, "SS")) msa_format = SS;
      else if (!strcmp(optarg, "LAV")) msa_format = LAV;
      else if (!strcmp(optarg, "MAF")) msa_format = MAF;
      else if (strcmp(optarg, "PHYLIP") != 0) { 
        fprintf(stderr, "Bad alignment format.\n");
        exit(1); 
      }
      break;
    case 'r':
      rseq_fname = optarg;
      break;
    case 'o':
      offset5 = atoi(optarg);
      break;
    case 'p':
      offset3 = atoi(optarg);
      break;
    case 'L':
      logf = fopen_fname(optarg, "w+");
      break;
    case 'M':
      mlogf = fopen_fname(optarg, "w+");
      break;
    case 'S':
      statsf = fopen_fname(optarg, "w+");
      break;
    case 'd':
      discardf = fopen_fname(optarg, "w+");
      break;
    case 'x':
      no_output = 1;
      break;
    case 'h':
      print_usage();
      exit(0);
    case '?':
      fprintf(stderr, "ERROR: Bad argument.  Try the --help option.\n");
      print_usage();
      exit(1);
    }
  }

  if (optind + 1 >= argc ) {
    fprintf(stderr, "ERROR:  Missing required arguments.  Try the --help option.\n");
    exit(1);
  }

  gff = gff_read_from_fname(argv[optind]);
  if (msa_format == MAF) {
    msa = maf_read(fopen_fname(argv[optind+1], "r"), 
                   rseq_fname == NULL ? NULL : fopen_fname(rseq_fname, "r"), 
                   NULL, 1, NULL, NULL, -1, 1, 0, NO_STRIP); 
  }
  else {
    msa = msa_new_from_file(fopen_fname(argv[optind+1], "r"),
                            msa_format, NULL); 
    if (msa->ss == NULL) 
      ss_from_msas(msa, 1, 1, NULL, NULL, NULL, -1);
  }
  assert(msa->ss->tuple_idx);

  gff_group(gff, groupby);

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
    f->start = msa_map_seq_to_msa(map, f->start);
    f->end = msa_map_seq_to_msa(map, f->end);
  }

  keepers = lst_new_ptr(lst_size(gff->features));
  if (discardf != NULL) discards = lst_new_ptr(lst_size(gff->features));

  ncons_tested = nkept = nconserved_exons = 0;
  for (i = 0; i < NTYPES; i++) nconsid[i] = 0;
  for (i = 0; i < NTYPES; i++) nfail[i] = 0;
  for (i = 0; i < NGAP_TYPES; i++) nce_gap_type[i] = 0;  

  for (i = 0; i < lst_size(gff->groups); i++) {
    GFF_FeatureGroup *group = lst_get_ptr(gff->groups, i);
    List *gfeatures = group->features;
    GFF_Feature *feat;
    status_type status = OKAY;
    cds_gap_type gt = FSHIFT_BAD;
    int no_alignment;
    lst_clear(failure_types);
    lst_clear(badfeats);

    /* First, exclude stop codons from cds's, if necessary (simplifies
       the detection of nonsense mutations). */
    exclude_stops(group, starts_adjusted, ends_adjusted);

    /* In all cases, discard any group for which the reference sequence
       doesn't have valid splice sites or start/stop codons, or has a
       premature stop codon */
    if (!ref_seq_okay(gfeatures, msa, offset3, indel_strict, splice_strict)) {
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
      for (j = 0, no_alignment = 0; j < lst_size(gfeatures); j++) { 
        feat = lst_get_ptr(gfeatures, j);
        if (str_equals_charstr(feat->feature, CDS) &&
            is_incomplete_alignment(feat, msa)) {
          status = NO_ALN;
          nfail[NO_ALN]++;
          break;
        }
      }

      if (status == OKAY) {     /* we have alignment and agreement
                                   with the ref seq; now check feature
                                   by feature  */

        lst_clear(intron_splice);
        for (j = 0; j < lst_size(gfeatures); j++) {
          feat = lst_get_ptr(gfeatures, j);

          if (feat->end - 1 >= msa->length) 
            die("ERROR: feature extends beyond alignment (%d >= %d).\n",
                feat->end - 1, msa->length);
        
          if (check_start && str_equals_charstr(feat->feature, START)) {

            nconsid[BAD_START]++;

            if (!is_conserved_start(feat, msa)) {
              status = BAD_START;
              lst_push_int(failure_types, BAD_START);
              lst_push_ptr(badfeats, feat);
            }
          }

          else if (check_stop && str_equals_charstr(feat->feature, STOP)) {

            nconsid[BAD_STOP]++;

            if (!is_conserved_stop(feat, msa)) {
              status = BAD_STOP;
              lst_push_int(failure_types, BAD_STOP);
              lst_push_ptr(badfeats, feat);
            }
          }

          else if (check_splice && 
                   str_equals_charstr(feat->feature, SPLICE_5)) {

            nconsid[BAD_5_SPLICE]++;

            if (!is_conserved_5splice(feat, msa, offset5, splice_strict)) {
              status = BAD_5_SPLICE;
              lst_push_int(failure_types, BAD_5_SPLICE);
              lst_push_ptr(badfeats, feat);
            }
            else lst_push_ptr(intron_splice, feat);
          }

          else if (check_splice && 
                   str_equals_charstr(feat->feature, SPLICE_5_UTR)) {

            nconsid[BAD_5_SPLICE_UTR]++;

            if (!is_conserved_5splice(feat, msa, offset5, splice_strict)) {
              status = BAD_5_SPLICE_UTR;
              lst_push_int(failure_types, BAD_5_SPLICE_UTR);
              lst_push_ptr(badfeats, feat);
            }
            else lst_push_ptr(intron_splice, feat);
          }

          else if (check_splice && str_equals_charstr(feat->feature, SPLICE_3)) {


            nconsid[BAD_3_SPLICE]++;

            if (!is_conserved_3splice(feat, msa, offset3, splice_strict)) {
              status = BAD_3_SPLICE;
              lst_push_int(failure_types, BAD_3_SPLICE);
              lst_push_ptr(badfeats, feat);
            }
            else lst_push_ptr(intron_splice, feat);
          }

          else if (check_splice && str_equals_charstr(feat->feature, SPLICE_3)) {

            nconsid[BAD_3_SPLICE_UTR]++;

            if (!is_conserved_3splice(feat, msa, offset3, splice_strict)) {
              status = BAD_3_SPLICE_UTR;
              lst_push_int(failure_types, BAD_3_SPLICE_UTR);
              lst_push_ptr(badfeats, feat);
            }
            else lst_push_ptr(intron_splice, feat);
          }

          else if (str_equals_charstr(feat->feature, CDS)) {
 
            if ((gt = get_cds_gap_type(feat, msa)) < fshift_mode) {
              if (status == OKAY || status == NONSENSE) status = FSHIFT;
              /* status records most basic type of failure; frame shifts
                 and nonsense mutations are often secondary */
              lst_push_int(failure_types, FSHIFT);
              lst_push_ptr(badfeats, feat);
            }

            if (check_nonsense && !is_nonsense_clean(feat, msa)) {
              if (status == OKAY) status = NONSENSE;
              lst_push_int(failure_types, NONSENSE);
              lst_push_ptr(badfeats, feat);
            }
          }
        } /* end loop through features in group */

        /* still have to make sure splice sites are paired correctly
           (GT-AG, GC-AG, AT-AC) */
        if (status == OKAY && !splice_strict && lst_size(intron_splice) >= 2 &&
            !are_introns_okay(intron_splice, msa, &feat, offset5, offset3)) {
          status = BAD_INTRON;
          lst_push_int(failure_types, BAD_INTRON);
          lst_push_ptr(badfeats, feat);
        }

        /* if collecting stats, record counts for failures */
        if (statsf != NULL) {
          if (status != OKAY) {
            for (j = 0; j < lst_size(failure_types); j++) {
              status_type ftype = lst_get_int(failure_types, j);

              if ((ftype == FSHIFT || ftype == NONSENSE) && 
                  status != FSHIFT && status != NONSENSE)
                continue;           /* don't count secondary frame shifts
                                       and nonsense mutations */ 

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
      if (mlogf != NULL) 
        write_machine_log(mlogf, group, status, gt, badfeats, 
                          failure_types, map);
    }
    else {
      if (discardf != NULL) {
        restore_stops(group, starts_adjusted, ends_adjusted);
        for (j = 0; j < lst_size(gfeatures); j++) 
          lst_push_ptr(discards, lst_get_ptr(gfeatures, j));
      }
      if (logf != NULL) 
        write_log(logf, group, status, gt, badfeats, failure_types, msa, map);
      if (mlogf != NULL)
        write_machine_log(mlogf, group, status, gt, badfeats, failure_types, map);
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
    fprintf(statsf, "#%11s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n", 
            "total", "nbad_ref", "nconsid", "nkept", "nno_aln", 
            "nbad_starts", "(out of)", "nbad_stops", "(out of)", 
            "nbad_5spl", "(out of)", "nbad_3spl", "(out of)", 
            "nbad_5utr", "(out of)", "nbad_3utr", "(out of)", 
            "nbad_intron", "nnons", "nfshifts", "ncons_exons", 
            "nce_ngaps", "nce_nov_cln", "nce_clean", "nce_fshftok");
    fprintf(statsf, "%12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d\n", 
            nfail[BAD_REF]+ncons_tested, nfail[BAD_REF], ncons_tested, nkept, 
            nfail[NO_ALN], nfail[BAD_START], nconsid[BAD_START], 
            nfail[BAD_STOP], nconsid[BAD_STOP], nfail[BAD_5_SPLICE], 
            nconsid[BAD_5_SPLICE], nfail[BAD_3_SPLICE], nconsid[BAD_3_SPLICE],
            nfail[BAD_5_SPLICE_UTR], nconsid[BAD_5_SPLICE_UTR],
            nfail[BAD_3_SPLICE_UTR], nconsid[BAD_3_SPLICE_UTR], 
            nfail[BAD_INTRON], nfail[NONSENSE], nfail[FSHIFT], 
            nconserved_exons, nce_gap_type[NGAPS], 
            nce_gap_type[NOVRLP_CLN_GAPS], nce_gap_type[CLN_GAPS], 
            nce_gap_type[FSHIFT_OK]);
    fprintf(statsf, STATS_DESCRIPTION);
  }

  if (logf != NULL) fclose(logf);
  if (mlogf != NULL) fclose(mlogf);
  if (statsf != NULL) fclose(statsf);
  if (discardf != NULL) fclose(discardf);

  return 0;
}


