/* $Id: clean_genes.c,v 1.1.1.1 2004-06-03 22:43:12 acs Exp $
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
#define CDS "cds"

/* "helper" feature types, used for 'indel-strict' mode */
#define HELPERS "start,stop,5'splice,3'splice,cds5'ss,cds3'ss,prestart"

/* possible values for status of each group of features */
typedef enum {OKAY, BAD_START, BAD_STOP, BAD_5_SPLICE, BAD_3_SPLICE, 
              NONSENSE_MUTATION, FSHIFT, NO_ALN} status_type;

/* possible gap types for cds exon -- frame-shift gaps (FSHIFT_BAD),
   no apparent frame shift (FSHIFT_OK), all gaps multiple of 3 in len
   (CLN_GAPS), gaps mult. of 3 and nonoverlapping (NOVRLP_CLN_GAPS),
   and no gaps present (NGAPS).  Defined in order of least to most
   stringent criterion */
typedef enum {FSHIFT_BAD, FSHIFT_OK, CLN_GAPS, 
              NOVRLP_CLN_GAPS, NGAPS} cds_gap_type;

/* for use in adjusting cds's to exclude stop codons (see below) */
typedef enum {NO_ADJUST, ADJUST_END, ADJUST_START} adjust_type;

/* minimum number of sites required for a gapless alignment block,
   when checking for frame-shift indels -- approximately the same as
   the maximum number of bases that can separate compensatory
   indels (see function 'is_fshift_okay') */
#define MIN_BLOCK_SIZE 30

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
                splice sites are valid in the reference sequence.\n\
                This can be used with an \"alignment\" consisting of a\n\
                single sequence to filter out incorrect annotations.\n\
                Options are available to impose additional criteria\n\
                as well, mostly having to do with conservation across species\n\
                (see the '--conserved' option in particular).\n\
\n\
USAGE:          clean_genes [options] <gff_fname> <msa_fname>\n\
\n\
OPTIONS:        \n\
\n\
    --start, -s\n\
        Require conserved start codons \n\
\n\
    --stop, -t\n\
        Require conserved stop codons \n\
\n\
    --splice, -l    \n\
        Require conserved canonical splice sites \n\
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
        Implies --start, --stop, --splice, --fshift, and --nonsense.\n\
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
NOTE:   Currently, feature types are defined at compile time as follows.\n\
\n\
               coding exon    <-> \"%s\"\n\
               start codon    <-> \"%s\"\n\
               stop codon     <-> \"%s\"\n\
               5' splice site <-> \"%s\"\n\
               3' splice site <-> \"%s\"\n\
\n\
        NOTE: with --fshift and --nonsense, it is possible for entries\n\
        to pass through that have stop codons in the frame of the\n\
        *reference* sequence, although they do not have any in their\n\
        own frame.  Use --clean-gaps instead to guarantee that no stop\n\
        codons occur in any sequence in the frame of the reference\n\
        sequence.\n\n", 
           CDS, START, STOP, SPLICE_5, SPLICE_3);
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
# (out of)      (number tested for start codons)\n\
# nbad_stops    number failing stop codon test\n\
# (out of)      (number tested for stop codons)\n\
# nbad_5spl     number failing 5' splice site test\n\
# (out of)      (number tested for 5' splice sites)\n\
# nbad_3spl     number failing 3' splice site test\n\
# (out of)      (number tested for 3' splice sites)\n\
# nnons         number failing nonsense mutation test (not counting ones that fail above tests)\n\
# nfshifts      number failing frameshift test (not counting ones that fail above tests)\n\
# *** ALIGNMENT GAP STATS ***\n\
# ncons_exons   number of conserved exons (same as nkept if --conserved)\n\
# nce_ngaps     number with no alignment gaps\n\
# nce_nov_cln   number with only nonoverlapping \"clean\" gaps (multiple of three lengths)\n\
# nce_clean     number with clean gaps that overlap\n\
# nce_fshftok   number with compensatory frame-shifting gaps as allowed by --fshift\n";

/* return 1 if all seqs for the specified column tuple and offset
   have the specified char; otherwise returns 0 */
inline int all_chars(MSA *msa, int tupleidx, char c, int offset) {
  int i;
  for (i = 0; i < msa->nseqs; i++) 
    if (ss_get_char_tuple(msa, tupleidx, i, offset) != c &&
        ss_get_char_tuple(msa, tupleidx, i, offset) != 'N') 
      return 0;
  return 1;
}

inline int is_conserved_start(GFF_Feature *feat, MSA *msa) {
  int start = feat->start - 1;  /* base 0 indexing */
  return ((feat->strand == '+' &&
           all_chars(msa, msa->ss->tuple_idx[start], 'A', 0) &&
           all_chars(msa, msa->ss->tuple_idx[start+1], 'T', 0) &&
           all_chars(msa, msa->ss->tuple_idx[start+2], 'G', 0)) ||
          (feat->strand == '-' &&
           all_chars(msa, msa->ss->tuple_idx[start+2], 'T', 0) &&
           all_chars(msa, msa->ss->tuple_idx[start+1], 'A', 0) &&
           all_chars(msa, msa->ss->tuple_idx[start], 'C', 0)));
}

inline int is_stop_codon(char *str) {
 return (strncmp(str, "TAA", 3) == 0 || strncmp(str, "TAG", 3) == 0 ||
         strncmp(str, "TGA", 3) == 0);
                                /* use of strncmp allows codons in the
                                   middle of larger strings to be
                                   tested */
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
    if (feat->strand == '-') msa_reverse_compl_seq(tuplestr, 3);
    if (!is_stop_codon(tuplestr) && strcmp(tuplestr, "NNN") != 0) return 0;
  }
  return 1;
}

inline int is_conserved_5splice(GFF_Feature *feat, MSA *msa, int offset5) {
  int start = feat->start - 1;  /* base 0 indexing */
  if (feat->strand == '-') start += offset5;

  return ((feat->strand == '+' && 
           all_chars(msa, msa->ss->tuple_idx[start], 'G', 0) &&
           all_chars(msa, msa->ss->tuple_idx[start+1], 'T', 0)) ||
          (feat->strand == '-' && 
           all_chars(msa, msa->ss->tuple_idx[start+1], 'C', 0) &&
           all_chars(msa, msa->ss->tuple_idx[start], 'A', 0)));
}

inline int is_conserved_3splice(GFF_Feature *feat, MSA *msa, int offset3) {
  int start = feat->start - 1;  /* base 0 indexing */
  if (feat->strand == '+') start += offset3;
  return ((feat->strand == '+' && 
           all_chars(msa, msa->ss->tuple_idx[start], 'A', 0) &&
           all_chars(msa, msa->ss->tuple_idx[start+1], 'G', 0)) ||
          (feat->strand == '-' && 
           all_chars(msa, msa->ss->tuple_idx[start+1], 'T', 0) &&
           all_chars(msa, msa->ss->tuple_idx[start], 'C', 0)));
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
    case NONSENSE_MUTATION:
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
    case NO_ALN:
      reason = "Incomplete alignment for cds";
      break;
    default: 
      assert(0);
    }

    fprintf(logf, "%s (%d-%d):\n", reason, msa_map_msa_to_seq(map, badfeat->start), 
            msa_map_msa_to_seq(map, badfeat->end));

    dump_aln(logf, badfeat, msa, ftype == NONSENSE_MUTATION);
  }
}

/* write machine-readable log entry for discarded feature */
void write_machine_log(FILE *mlogf, GFF_FeatureGroup *group, status_type status, 
                       cds_gap_type gt, List *badfeats, List *failure_types, 
                       msa_coord_map *map) {

  int i;

  for (i = 0; i < lst_size(failure_types); i++) {
    status_type ftype = lst_get_int(failure_types, i);
    GFF_Feature *badfeat = lst_get_ptr(badfeats, i);
    char *reason;

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
    case NONSENSE_MUTATION:
      reason = "nonsense";
      break;
    case FSHIFT:
      reason = "frameshift";
      break;
    case NO_ALN:
      reason = "no_alignment";
      break;
    default: 
      assert(0);
    }

    fprintf(mlogf, "%s\t%d\t%d\t%s\n", group->name->chars, 
            msa_map_msa_to_seq(map, badfeat->start), 
            msa_map_msa_to_seq(map, badfeat->end),
            reason);
  }
}

/* checks to see if reference sequence looks okay wrt a given
   list of features */
int ref_seq_okay(List *features, MSA *msa, int offset3, int indel_strict) {
  static List *helpers = NULL;
  static char *seq = NULL;
  static int seqalloc = 0;
  int idx;
  GFF_Feature *feat, *lastfeat_helper = NULL;

  if (indel_strict && helpers == NULL) {
    helpers = lst_new_ptr(10);
    str_split(str_new_charstr(HELPERS), ",", helpers);
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
    else if (str_equals_charstr(feat->feature, SPLICE_5) && 
             strncmp(seq, "GT", 2) != 0) 
      return 0;
    else if (str_equals_charstr(feat->feature, SPLICE_3) &&
             strncmp(&seq[offset3], "AG", 2) != 0) 
      return 0;
    else if (str_equals_charstr(feat->feature, CDS)) {
      for (i = (3 - feat->frame) % 3; i <= len - 3; i += 3) 
        if (is_stop_codon(&seq[i])) return 0;
    }

    if (indel_strict) {
      if (str_in_list(feat->feature, helpers)) {
        /* reject any helper feature with gaps in the ref seq, unless they
           appear in a non-critical part of a splice site or in a
           "prestart" feature  */
        if (has_gaps) {          
          if (str_equals_charstr(feat->feature, SPLICE_5)) {
            if (ss_get_char_pos(msa, feat->start-1, 0, 0) == GAP_CHAR ||
                ss_get_char_pos(msa, feat->start, 0, 0) == GAP_CHAR)
              return 0;
          }
          else if (str_equals_charstr(feat->feature, SPLICE_3)) {
            if (ss_get_char_pos(msa, feat->end-1, 0, 0) == GAP_CHAR ||
                ss_get_char_pos(msa, feat->end-2, 0, 0) == GAP_CHAR)
              return 0;
          }
          else if (!str_equals_charstr(feat->feature, "prestart"))
            return 0;
        }
        /* in addition, if two helpers occur consec. with gaps and
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

/* Exclude stop codons from all CDS in a group, as necessary.  Record any features
   that are changed, so they can be changed back before data is
   output */
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

/* get range of group */
inline void get_range(GFF_FeatureGroup *group, int *start, int *end) {
  int i;
  *start = INFTY; *end = -1;
  for (i = 0; i < lst_size(group->features); i++) {
    GFF_Feature *f = lst_get_ptr(group->features, i);
    if (f->start < *start) *start = f->start;
    if (f->end > *end) *end = f->end;
  }
}

int main(int argc, char *argv[]) {

  int check_start = 0, check_stop = 0, check_splice = 0, check_nonsense = 0,
    offset5 = 0, offset3 = 0, opt_idx, i, j, indel_strict = 0, no_output = 0,
    check_alignment = 0;
  int nconsidered, nkept, nconserved_exons,
    nstarts_considered, nstops_considered, nsplice_5_considered, 
    nsplice_3_considered, nbad_ref, s, e;
  int nce_gap_type[NGAPS+1], nfail[NO_ALN+1];
  char c;
  MSA *msa;
  GFF_Set *gff;
  msa_format_type msa_format = PHYLIP;
  List *keepers, *badfeats = lst_new_ptr(10), *failure_types = lst_new_int(10),
    *ends_adjusted = lst_new_ptr(1), *starts_adjusted = lst_new_ptr(1), *discards;
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

  while ((c = getopt_long(argc, argv, "i:r:L:M:S:g:d:stlnfceIxh", 
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

  nconsidered = nkept = nfail[BAD_START] = nfail[BAD_STOP] = nfail[BAD_5_SPLICE] =
    nfail[BAD_3_SPLICE] = nfail[NONSENSE_MUTATION] = nfail[FSHIFT] = 
    nfail[NO_ALN] = nconserved_exons = 
    nce_gap_type[FSHIFT_BAD] = nce_gap_type[FSHIFT_OK] = 
    nce_gap_type[CLN_GAPS] = nce_gap_type[NOVRLP_CLN_GAPS] = 
    nce_gap_type[NGAPS] = nstarts_considered = nstops_considered =
    nsplice_5_considered = nsplice_3_considered = nbad_ref = 0;

  for (i = 0; i < lst_size(gff->groups); i++) {
    GFF_FeatureGroup *group = lst_get_ptr(gff->groups, i);
    List *gfeatures = group->features;
    GFF_Feature *feat;
    status_type status = OKAY;
    cds_gap_type gt = FSHIFT_BAD;
    int no_alignment;

    /* First, exclude stop codons from cds's, if necessary (simplifies
       the detection of nonsense mutations). */
    exclude_stops(group, starts_adjusted, ends_adjusted);

    /* In all cases, discard any group for which the reference sequence
       doesn't have valid splice sites or start/stop codons, or has a
       premature stop codon */
    if (!ref_seq_okay(gfeatures, msa, offset3, indel_strict)) {
      if (logf != NULL)
        fprintf(logf, "****\nSkipped '%s' -- features not consistent with reference sequence.\n", group->name->chars);
      if (mlogf != NULL) {
        get_range(group, &s, &e);
        fprintf(mlogf, "%s\t%d\t%d\tbad_ref\n", group->name->chars, 
                msa_map_msa_to_seq(map, s), msa_map_msa_to_seq(map, e));
      }
      nbad_ref++;
      if (discardf != NULL) {
        restore_stops(group, starts_adjusted, ends_adjusted);
        for (j = 0; j < lst_size(gfeatures); j++) 
          lst_push_ptr(discards, lst_get_ptr(gfeatures, j));
      }
      continue;
    }

    /* Everything else counts as a potentially valid group */
    nconsidered++;

    /* If not interested in cross-species conservation at all, can
       simply continue */
    if (!check_alignment) {
      nkept++;
      for (j = 0; j < lst_size(gfeatures); j++)
        lst_push_ptr(keepers, lst_get_ptr(gfeatures, j));
      continue;
    }

    /* Otherwise (we are interested in cross-species conservation),
       check first to make sure there's alignment across species in
       the cds; if not, there's no need to look at individual
       features. */
    for (j = 0, no_alignment = 0; j < lst_size(gfeatures) && !no_alignment; j++) { 
      feat = lst_get_ptr(gfeatures, j);
      if (str_equals_charstr(feat->feature, CDS) &&
          is_incomplete_alignment(feat, msa)) 
        no_alignment = 1;
    }
    if (no_alignment) {
      if (logf != NULL) 
        fprintf(logf, "****\nSkipped '%s' -- lack of alignment across species.\n", 
                group->name->chars);
      if (mlogf != NULL) {
        get_range(group, &s, &e);
        fprintf(mlogf, "%s\t%d\t%d\tno_alignment\n", group->name->chars, 
                msa_map_msa_to_seq(map, s), msa_map_msa_to_seq(map, e));
      }
      nfail[NO_ALN]++;
      if (discardf != NULL) {
        restore_stops(group, starts_adjusted, ends_adjusted);
        for (j = 0; j < lst_size(gfeatures); j++) 
          lst_push_ptr(discards, lst_get_ptr(gfeatures, j));
      }
      continue;
    }

    /* Now check feature by feature */
    lst_clear(failure_types);
    lst_clear(badfeats);
    for (j = 0; j < lst_size(gfeatures); j++) {
      feat = lst_get_ptr(gfeatures, j);

      if (feat->end - 1 >= msa->length) 
        die("ERROR: feature extends beyond alignment (%d >= %d).\n",
            feat->end - 1, msa->length);
        
      if (check_start && status != NO_ALN &&
          str_equals_charstr(feat->feature, START)) {

        nstarts_considered++;

        if (!is_conserved_start(feat, msa)) {
          status = BAD_START;
          lst_push_int(failure_types, BAD_START);
          lst_push_ptr(badfeats, feat);
        }
      }

      else if (check_stop && status != NO_ALN &&
               str_equals_charstr(feat->feature, STOP)) {

        nstops_considered++;

        if (!is_conserved_stop(feat, msa)) {
          status = BAD_STOP;
          lst_push_int(failure_types, BAD_STOP);
          lst_push_ptr(badfeats, feat);
        }
      }

      else if (check_splice && status != NO_ALN &&
               str_equals_charstr(feat->feature, SPLICE_5)) {

        nsplice_5_considered++;

        if (!is_conserved_5splice(feat, msa, offset5)) {
          status = BAD_5_SPLICE;
          lst_push_int(failure_types, BAD_5_SPLICE);
          lst_push_ptr(badfeats, feat);
        }
      }

      else if (check_splice && status != NO_ALN &&
               str_equals_charstr(feat->feature, SPLICE_3)) {

        nsplice_3_considered++;

        if (!is_conserved_3splice(feat, msa, offset3)) {
          status = BAD_3_SPLICE;
          lst_push_int(failure_types, BAD_3_SPLICE);
          lst_push_ptr(badfeats, feat);
        }
      }

      else if (str_equals_charstr(feat->feature, CDS)) {
 
        if ((gt = get_cds_gap_type(feat, msa)) < fshift_mode) {
          if (status == OKAY || status == NONSENSE_MUTATION) status = FSHIFT;
                                /* status records most basic type
                                   of failure; frame shifts and
                                   nonsense mutations are often
                                   secondary */
          lst_push_int(failure_types, FSHIFT);
          lst_push_ptr(badfeats, feat);
        }

        if (check_nonsense && !is_nonsense_clean(feat, msa)) {
          if (status == OKAY) status = NONSENSE_MUTATION;
          lst_push_int(failure_types, NONSENSE_MUTATION);
          lst_push_ptr(badfeats, feat);
        }
      }
    }

    /* if collecting stats, record counts for failures */
    if (statsf != NULL) {
      if (status != OKAY) {
        for (j = 0; j < lst_size(failure_types); j++) {
          status_type ftype = lst_get_int(failure_types, j);

          if ((ftype == FSHIFT || ftype == NONSENSE_MUTATION) && 
              status != FSHIFT && status != NONSENSE_MUTATION)
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

    if (status == OKAY) {
      nkept++;
      if (!no_output) {
        restore_stops(group, starts_adjusted, ends_adjusted);
        for (j = 0; j < lst_size(gfeatures); j++)
          lst_push_ptr(keepers, lst_get_ptr(gfeatures, j));
      }
      if (mlogf != NULL) {
        get_range(group, &s, &e);
        fprintf(mlogf, "%s\t%d\t%d\tokay\n", group->name->chars, 
                msa_map_msa_to_seq(map, s), msa_map_msa_to_seq(map, e));
      }

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


  if (statsf != NULL) {
    fprintf(statsf, "#%11s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\n", "total", "nbad_ref", "nconsid", "nkept", "nno_aln", "nbad_starts", "(out of)", "nbad_stops", "(out of)", "nbad_5spl", "(out of)", "nbad_3spl", "(out of)", "nnons", "nfshifts", "ncons_exons", "nce_ngaps", "nce_nov_cln", "nce_clean", "nce_fshftok");
    fprintf(statsf, "%12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d %12d\n", nbad_ref+nconsidered, nbad_ref, nconsidered, nkept, nfail[NO_ALN], nfail[BAD_START], nstarts_considered, nfail[BAD_STOP], nstops_considered, nfail[BAD_5_SPLICE], nsplice_5_considered, nfail[BAD_3_SPLICE], nsplice_3_considered, nfail[NONSENSE_MUTATION], nfail[FSHIFT], nconserved_exons, nce_gap_type[NGAPS], nce_gap_type[NOVRLP_CLN_GAPS], nce_gap_type[CLN_GAPS], nce_gap_type[FSHIFT_OK]);
    fprintf(statsf, STATS_DESCRIPTION);
  }

  if (logf != NULL) fclose(logf);
  if (mlogf != NULL) fclose(mlogf);
  if (statsf != NULL) fclose(statsf);
  if (discardf != NULL) fclose(discardf);

  return 0;
}


