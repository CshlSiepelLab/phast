/* label - label the columns of alignment(s) by category */

/* $Id: exoniphy.c,v 1.6 2004-06-15 22:33:57 acs Exp $
   Written by Adam Siepel, 2002 and 2003
   Copyright 2002, Adam Siepel, University of California 

   NOTE: This program needs a rewrite!  Start by separating training
   and testing modes into two programs ...
*/


#include <stdlib.h>
#include <stdio.h>
#include <msa.h>
#include <hmm.h>
#include <assert.h>
#include <getopt.h>
#include <string.h>
#include <tree_likelihoods.h>
#include <gff.h>
#include <bed.h>
#include <category_map.h>
#include <dgamma.h>
#include <numerical_opt.h>
#include <sufficient_stats.h>
#include <stringsplus.h>
#include <puzzler.h>
#include <gap_patterns.h>
/* #include <em.h> */

#define STARTING_LAMBDA 0.9     /* used for initialization when
                                   autofitting lambda */

/* feature tupes used for scoring predictions; hardwired for now */
#define SCORE_CATS "cds,cds1,cds2,cds5'ss,cds3'ss,start"
#define HELPER_CATS "stop,5'splice,3'splice"
#define NULL_CATS "background,cnc,icnc"
                                /* cat names that aren't present will
                                   be ignored */

/* categories for which complex gap patterns are prohibited;
   temporarily hardwired */
#define NO_COMPLEX "cds,cds5'ss,cds3'ss"

/* features for which frame will be inferred; also temporarily
   hardwired */
#define FRAME_CATS "cds,cds1,cds2"

/* parameters controlling evaluation of Sn/Sp tradeoff (see -Y option) */
#define SCALE_RANGE_MIN -20
#define SCALE_RANGE_MAX 10
#define NSENS_SPEC_VERSIONS 10

#define MIN_BLOCK_SIZE 30
/* defines block size for -x option */

typedef enum {TRAINING, TESTING  /* , EM */} mode_type;

void print_usage() {
    printf("\n\
PROGRAM: label\n\
\n\
USAGE: \n\
    (training) \n\
        label -m <msa_fname_list> -c <category_map_fname> \n\
              -g <gff_fname_list> -t <out_hmm_fname> [OPTIONS] \n\
\n\
    (testing)  \n\
        label -m <msa_fname> -c <category_map_fname> \n\
              -d <model_fname_list> [OPTIONS] \n\
\n\
DESCRIPTION: \n\
    Labels the columns of alignment(s) by category.  Runs in two\n\
    modes: a TRAINING mode, in which column labels are obtained from\n\
    explicit sequence annotations, and are used to train an\n\
    HMM; and a TESTING mode, in which labels are predicted according to\n\
    a combined phylogenetic and hidden Markov model (specified as\n\
    input).  Note that training here refers only to the transition\n\
    probabilities of the HMM; the emission probabilities are defined\n\
    by phylogenetic models, which must be fitted separately (see the\n\
    program fit_tree_model).  More than one alignment may be specified\n\
    in training mode, but currently only one may be specified in testing\n\
    mode.\n\
\n\
OPTIONS:\n\
    -m <msa_fname_list>\n\
        (required) List of multiple sequence alignment files.\n\
        Currently, in testing mode, the list must be of length one.\n\
\n\
    -i PHYLIP|FASTA|MPM|SS \n\
        (default FASTA) Alignment format.\n\
 \n\
    -c <category_map_fname>\n\
        (required) File defining mapping of feature types to category\n\
        numbers.\n\
\n\
    -g <gff_fname_list>\n\
        (indicates training mode)  Files in GFF defining sequence\n\
        features to be used in labeling sites.   Frame of reference of \n\
        feature indices is determined feature-by-feature according to \n\
         'seqname' attribute.  Filenames must correspond in number and order\n\
        to the elements of <msa_fname_list>. \n\
\n\
    -d <model_fname_list>\n\
        (indicates testing mode)  List of files defining a tree model\n\
        for each functional category.  Order of models must correspond\n\
        to order of states in HMM.  If -d is specified and -H is not,\n\
        then a trivial HMM is assumed (a single state, with\n\
        probability 1 of transitioning to itself).  In such a case,\n\
        <model_fname_list> must have only one element.  Tree model\n\
        files may be produced with fit_tree_model.  Note: this option\n\
        is given a special interpretation with -D (see below).\n\
\n\
    -H <hmm_fname>\n\
        (for use with -d; indicates testing mode)  Name of HMM file,\n\
        defining the probability of transition from each functional\n\
        category to each other.  Generally a file is used that was\n\
        produced by this same program, running in training mode.\n\
\n\
    -R\n\
        Reverse complement.  In training mode, causes reverse\n\
        complementation of segments of alignment corresponding to\n\
        groups of related features on the reverse strand, before\n\
        transition probabilities are estimated.  In testing mode,\n\
        causes predictions to be made for the reverse\n\
        complement of the entire alignment.\n\
\n\
    -G  ALL|ANY|<seqno>\n\
        Strip columns in alignment(s) containing all gaps, any gaps, or \n\
        a gap in the specified sequence (<seqno>; indexing starts with\n\
        one).  Default is not to strip any columns.\n\
\n\
    -o <output_fname_root>\n\
        Use specified filename root for output files, rather than \n\
        \"label-out\" (the default).\n\
\n\
    -f\n\
        Write output (training labels or Viterbi path) in GFF format.\n\
        This option will be ignored if -l or -F is selected.  \n\
\n\
    -b <bed_feature_list>\n\
        Write output also in BED format, for specified feature types.\n\
        This option will be ignored if -l or -F is selected.  It may\n\
        be used with -f.\n\
\n\
    -s <source_name>\n\
        (For use with -f, -b, or -Q)  Use specified string as the \"source\"\n\
        field in generated GFF, BED, or samples file (e.g., chr22).\n\
\n\
    -q \n\
        Proceed quietly (without updates to stderr).\n\
\n\
    -h\n\
        Print this help message.\n\
\n\
\n\
Auxiliary options, training mode:\n\
    -t <out_hmm_fname>\n\
        Train an HMM based on observed labels.  Training will reflect\n\
        the aggregate of all specified alignments and annotation\n\
        files.  New model will be written to <out_hmm_fname>.\n\
\n\
    -M <msa_length_list>\n\
        (For use with -t; mutually exclusive with -m) Assume alignments\n\
        of the specified lengths (comma-separated list) and do not not\n\
        attempt to map the coordinates in the specified GFFs (assume\n\
        they are in the desired coordinate frame).  This option allows\n\
        an HMM to be trained directly from GFFs, without alignments.\n\
\n\
\n\
Auxiliary options, testing mode:\n\
    -p <state_number_list>\n\
        Compute and output posterior probabilities for the specified states.\n\
        A row for each position and a column for each state will be written\n\
        to a file with suffix \"postprob\".  Argument can be the string \"all\".\n\
\n\
    -j <refseq_idx>\n\
        (For use with -p and -Q)  Index posterior probabilities according \n\
        to coordinate frame of specified sequence (1-based indexing; \n\
        the value 0 indicates the frame of the entire multiple alignment, \n\
        which is the default).\n\
\n\
    -x \n\
        (For use with -p and -Q)  Don't report posterior probabilities for\n\
        blocks of sites at which only the reference sequence is present \n\
        (as determined by -j).  Threshold for block size is %d.\n\
\n\
    -l\n\
        (Cannot be used with -F) Output total log likelihood only,\n\
        rather than site-by-site labels.  Likelihoods are computed\n\
        using the forward algorithm.\n\
\n\
    -r <scaling_const_list>\n\
        List of rate constants by which to scale each tree model.  If\n\
        m rates are specified, each state in the original HMM\n\
        (specified by -h, or consisting of a single state, if -h is\n\
        not specified) will be represented by m states in a new HMM,\n\
        each of which represents a functional-category/rate-category\n\
        pair.  Transition probabilities between any two new states\n\
        will be defined as the product of the transition probability\n\
        between corresponding functional categories (defined by the\n\
        original HMM) and rate categories (defined by the\n\
        autocorrelation parameter; see below).  See Siepel and\n\
        Haussler, 2003 for a more detailed discussion.\n\
\n\
    -L <autocorr_param>\n\
        (For use with -r)  Autocorrelation parameter (\"lambda\" in the\n\
        literature), defining the tendency of each column to obey the\n\
        same evolutionary rate as the previous column.  Allowable\n\
        range is 0-1.  Default is 0 (no autocorrelation).  With m rate\n\
        categories, the probability of changing from any one rate\n\
        category to any other (when moving from a column its\n\
        successor) is (1-lambda)/m, and the probability of remaining\n\
        in a rate category is lambda + (1-lambda)/m.  See Siepel and\n\
        Haussler, 2003.\n\
\n\
    -A \n\
        Automatically determine rate constants and autocorrelation\n\
        parameter, to maximize likelihood.  Model must have been fitted\n\
        using the discrete gamma approximation and the parameter\n\
        \"alpha\" must have been estimated.  This approach does not\n\
        precisely find the MLE, but approximates it, by choosing the\n\
        rate constants according to alpha, then choosing lambda to\n\
        maximize the likelihood.  Generally, use this option instead\n\
        of -r/-L.  If it is used *with* L, then the rate constants\n\
        will be chosen according to alpha, but lambda will be fixed as\n\
        defined.  See Siepel and Haussler, 2003.\n\
\n\
    -k <nratecats>\n\
        (For use with -A).  Use <nratecats> rate categories instead of the\n\
        number given in the specified tree model files.\n\
\n\
    -U \n\
        (Alternative to -R).  Given an HMM describing the forward\n\
        strand, create a larger HMM that allows for features on both\n\
        strands by \"reflecting\" the HMM about all background states\n\
        (see -V).  The new HMM will be used for predictions on both strands.\n\
\n\
    -S\n\
        Assign log-odds scores to predictions, equal to their log total\n\
        probability under the given model minus their log total probability\n\
        under a null model (which consist of only background states -- see -V).\n\
\n\
    -X <pred_factor>\n\
        Multiply transition probabilities from background states into \n\
        non-background states (see -V) by the specified prediction factor,\n\
        then renormalize.  Useful when an HMM has been trained on a sparse \n\
        set of features.\n\
\n\
    -V <backgd_cat_names>\n\
        Categories to be considered \"background\", in addition to the\n\
        standard background category (category number 0).  Associated\n\
        states are considered \"background states\".  Affects -U, -S, \n\
        and -X.\n\
\n\
\n\
 (Experimental options)\n\
    -I <indel_cat_list>\n\
        Model indels for specified categories.  To have\n\
        nonzero probability for the states corresponding to a\n\
        specified category range, indels must be \"clean\"\n\
        (nonoverlapping), must be assignable by parsimony to a single\n\
        branch in the phylogenetic tree, and must have lengths that\n\
        are exact multiples of the category range size.  Avoid -G with\n\
        this option.  If used in training mode, requires -T.\n\
\n\
  (training)\n\
    -T <tree_fname>\n\
        (For use with -I).  Use the specified tree topology when training\n\
        for indels. \n\
\n\
    -n <nseqs> \n\
        (For use with -I).  Train an indel model for <nseqs>\n\
        sequences, despite that the training alignment has a different\n\
        number.  All (non-trivial) gap patterns are assumed to be\n\
        equally frequent.\n\
\n\
    -P <pseudocount_matrix_fname>\n\
        (For use with -t).  When training, use pseudocounts for\n\
        category-to-category transitions as defined in\n\
        <pseudocount_matrix_fname> (one row per line, columns\n\
        separated by whitespace; rows and columns must correspond in\n\
        order to category indices).  FIXME: not yet implemented!\n\
\n\
  (testing)\n\
    -D <gc-thresholds>\n\
        (Changes interpretation of -d) Use different sets of tree\n\
        models, depending on the G+C content of the input alignment.\n\
        The list <gc-thresholds> must consist of x ordered values in\n\
        (0,1), defining x+1 G+C classes.  The argument to -d must then\n\
        consist of the names of x+1 files, each of which contains a list\n\
        of tree-model filenames.\n\
\n\
    -F <label_fname_list>\n\
        (Cannot be used with -l) Output total log likelihood, for the\n\
        path defined by the labels in the specified files (in other\n\
        words, the joint probability of the path and the alignment\n\
        given the tree models).  \n\
\n\
    -C <gff_fname>\n\
        Use a priori labels, as specified in <gff_fname>.  Compute\n\
        likelihoods simply by computing the likelihood at each site\n\
        according to the designated categories, then sum\n\
        across all sites (the HMM will be ignored).\n\
\n\
    -Y\n\
        Make predictions for a range of different prediction factors, on\n\
        an exponential scale.  Allows analysis of the sensitivity/specificity\n\
        tradeoff.  Currently the range is hardcoded to be from exp(%d) to \n\
        exp(%d), with %d different sets predicted.\n\
\n\
    -W <cat_list>\n\
        Prohibit gaps in the specified categories (gaps result in\n\
        emission probabilities of zero).  By default, gaps are treated\n\
        as missing data.\n\
\n\
    -Z <cat_list>\n\
        (Similar to -W).  Allow gaps, but penalize them by assigning\n\
        them a probability equal to that of the least likely character\n\
        in the alphabet.  Currently only works for weight-matrix\n\
        categories.\n\
\n\
    -Q <state_number_list>\n\
        (Experimental, for use with -p)  Output a \"samples\" file for the \n\
        posterior probabilities of the specified states.  Can be used to \n\
        create a \"wiggle\" track in the UCSC browser.\n\
\n\
\n\
NOTE ON OPTIONS: \n\
    An arguments denoted as a list (e.g., <scaling_const_list>) may be\n\
    specified either directly, as a list of comma-separated strings, or\n\
    indirectly, by a filename preceded by an asterisk.  If\n\
    an asterisk is used, the elements of a list are read from the\n\
    specified filename.  For example, -r\"1,2,3\" and -r\"*numbers.txt\" \n\
    are equivalent, if the file \"numbers.txt\" contains 1, 2, and 3 \n\
    (separated by whitespace).  (The \"*\" can be thought of as like the \n\
    \"dereferencing\" operator in C.)\n\n", MIN_BLOCK_SIZE, SCALE_RANGE_MIN, 
           SCALE_RANGE_MAX, NSENS_SPEC_VERSIONS);
}

/* Defunct options: */
/*     -E  \n\ */
/*         Train by EM, rather than from labeled data.  This option\n\ */
/*         produces a sort of hybrid of training and testing modes.  As\n\ */
/*         in testing mode, option -g is not required, and -H and -d are\n\ */
/*         required (for initialization).  As in training mode, more than\n\ */
/*         one multiple alignment is allowed, as are options -t and -P.\n\ */
/*         Currently none of the auxiliary options of testing mode are\n\ */
/*         allowed.\n\ */

/*     -e  \n\ */
/*         Compute and output expected number of substitutions per site.\n\ */
/*         These are posterior estimates, which consider the forward and\n\ */
/*         backward values for the HMM.  They will appear as an\n\ */
/*         additional column in the output file.\n\ */
/* \n\ */


/* FIXME: if -r and -L are specified with -F, should sum over all
   possible assignments of rate categories, using the forward alg. */


typedef struct {                /* "package" for data needed by
                                   log_likelihood_wrapper (see below) */
  double **forward_scores;
  double **emission_scores;
  int msa_len;
  PhyloHMM_Puzzler *puzzler;
} AutoratesData;

double log_likelihood_wrapper(double lambda, void *data) {
  AutoratesData *ad = (AutoratesData*)data;
  if (lambda < 0 || lambda > 1) return INFTY;
  puz_update_cross_prod(ad->puzzler, lambda);
  return -hmm_forward(ad->puzzler->hmm, ad->emission_scores, ad->msa_len, 
                      ad->forward_scores);
}

double fit_lambda(PhyloHMM_Puzzler *puz, double **emission_scores, 
                  int msa_len) {
  AutoratesData ad;
  double lambda, final_score, ax, bx, cx, fa, fb, fc;
  int i;

  /* allocate memory for forward alg */
  ad.forward_scores = smalloc(puz->hmm->nstates * sizeof(double*));
  for (i = 0; i < puz->hmm->nstates; i++)
    ad.forward_scores[i] = smalloc(msa_len * sizeof(double));

  ad.emission_scores = emission_scores;
  ad.puzzler = puz;
  ad.msa_len = msa_len;

  ax = .80; bx = .97;           /* FIXME -- parameterize */
  mnbrak(&ax, &bx, &cx, &fa, &fb, &fc, log_likelihood_wrapper, &ad, NULL);
  final_score = opt_brent(ax, bx, cx, log_likelihood_wrapper, 5e-3, &lambda, &ad, NULL);
/*   fprintf(stderr, "Returned from opt_brent; freeing forward_scores ...\n"); */

  for (i = 0; i < puz->hmm->nstates; i++) free(ad.forward_scores[i]);
  free(ad.forward_scores);

  return lambda;
}

/* score a set of predicted features using log odds scoring.  gff
   contains features to score, score_cats are the names of the
   categories to score, helper cats are secondary categories to be
   included in scoring if adjacent to a primary category (e.g., splice
   sites), HMM and emissions should be as used for prediction, cm and
   puz are as usual, offset is from msa */
void score_predictions(GFF_Set *gff, List *score_cats, List *helper_cats,  
                       List *null_cats, double **emissions, 
                       CategoryMap *cm, PhyloHMM_Puzzler *puz, int offset)  {
  int i, j, cat, ncats;
  List *cats, *score_states, *null_states;
  List **cat_to_states;

  /* convert lists of cat names to lists of states */

  /* first create inverse mapping from category numbers to lists of
     states.  FIXME: should have this available in puzzler (partially
     implemented, I believe)  */
  ncats = cm->ncats + 1;
  cat_to_states = smalloc(ncats * sizeof(List*));
  for (i = 0; i < ncats; i++) 
    cat_to_states[i] = lst_new_int(2*puz->hmm->nstates/ncats);
  for (i = 0; i < puz->hmm->nstates; i++) 
    lst_push_int(cat_to_states[puz->state_to_cat[i]], i);

  /* now fill in state lists */
  score_states = lst_new_int((lst_size(score_cats) + 
                              lst_size(helper_cats)) * 10);
  null_states = lst_new_int(20);

  cats = cm_get_category_list(cm, score_cats, 1);
  for (i = 0; i < lst_size(cats); i++) {
    cat = lst_get_int(cats, i);
    for (j = 0; j < lst_size(cat_to_states[cat]); j++)
      lst_push_int(score_states, lst_get_int(cat_to_states[cat], j));
  }
  lst_free(cats);

  cats = cm_get_category_list(cm, helper_cats, 1);
  for (i = 0; i < lst_size(cats); i++) {
    cat = lst_get_int(cats, i);
    for (j = 0; j < lst_size(cat_to_states[cat]); j++)
      lst_push_int(score_states, lst_get_int(cat_to_states[cat], j));
  }
  lst_free(cats);

  /* null states; first assume background states */
  for (j = 0; j < lst_size(cat_to_states[0]); j++)
    lst_push_int(null_states, lst_get_int(cat_to_states[0], j));
  /* now check for others */
  if (null_cats != NULL) {
    cats = cm_get_category_list(cm, null_cats, 1);
    for (i = 0; i < lst_size(cats); i++) {
      cat = lst_get_int(cats, i);
      for (j = 0; j < lst_size(cat_to_states[cat]); j++)
        lst_push_int(null_states, lst_get_int(cat_to_states[cat], j));
    }
    lst_free(cats);
  }

  /* now score each feature */
  for (i = 0; i < lst_size(gff->features); i++) {
    GFF_Feature *feat = lst_get_ptr(gff->features, i);
    if (str_in_list(feat->feature, score_cats)) {
      int start = feat->start;
      int end = feat->end;

      /* extend range as far as possible using helper cats */
      for (j = i-1; j >= 0; j--) {
        GFF_Feature *prev_feat = lst_get_ptr(gff->features, j);
        if (str_in_list(prev_feat->feature, helper_cats) && 
            prev_feat->end == start - 1) 
          start = prev_feat->start;
        else break;
      }
      for (j = i+1; j < lst_size(gff->features); j++) {
        GFF_Feature *next_feat = lst_get_ptr(gff->features, j);
        if (str_in_list(next_feat->feature, helper_cats) && 
            next_feat->start == end + 1) 
          end = next_feat->end;
        else break;
      }

      /* score from start to end */
      feat->score = 
        hmm_log_odds_subset(puz->hmm, emissions, score_states, null_states, 
                            start - 1 - offset, end - start + 1);

      feat->score_is_null = 0;
    }
  }

  lst_free(score_states);
  lst_free(null_states);
  for (i = 0; i < ncats; i++) lst_free(cat_to_states[i]);
  free(cat_to_states);

}

/* post-processing test for a predicted feature.  This version is
   designed for cds predictions and is primarily concerned with
   alignment gaps (necessary until better gap handling can be
   integrated into the model).  */
int valid_feat(GFF_Feature *pred_feat, MSA *msa) {
  int i, j, msa_start, msa_end;

  /* first, special case for exon pairs: throw out any exon < 30 bp in
     length (compensates for "filler exons") */
  if ((str_equals_charstr(pred_feat->feature, "cds1") ||
       str_equals_charstr(pred_feat->feature, "cds2")) &&
      pred_feat->end - pred_feat->start + 1 < 30)
    return 0;

  /* now, if introns are present, discard any that are most likely too short */
  if (str_equals_charstr(pred_feat->feature, "intron") &&
      pred_feat->end - pred_feat->start + 1 < 50)
    return 0;  

  if (!str_equals_charstr(pred_feat->feature, "cds") &&
      !str_equals_charstr(pred_feat->feature, "cds1") &&
      !str_equals_charstr(pred_feat->feature, "cds2"))
    return 1;

  /* finally, we'll look at gaps in cds features: discard any such
     feature with *any* gaps in *any* sequence, unless only perfectly
     in-frame, internal gaps of relatively short length are present */

  /* assume ordered sufficient statistics */
  if (msa->ss == NULL) ss_from_msas(msa, 1, 1, NULL, NULL, NULL, -1);
  assert(msa->ss->tuple_idx != NULL);

  msa_start = pred_feat->start - 1; /*  - msa->idx_offset; (now handled after) */
  msa_end = pred_feat->end - 1; /*  - msa->idx_offset; */

  for (i = msa_start; i <= msa_end ; i++) {
    for (j = 0; j < msa->nseqs; j++) {
      if (ss_get_char_pos(msa, i, j, 0) == GAP_CHAR) {
        int gap_start, gap_end;
        for (gap_start = i-1; gap_start >= msa_start && 
               ss_get_char_pos(msa, gap_start, j, 0) == GAP_CHAR; gap_start--);
        gap_start++;
        for (gap_end = i+1; gap_end <= msa_end && 
               ss_get_char_pos(msa, gap_end, j, 0) == GAP_CHAR; gap_end++);
        gap_end--;

        if (gap_start == msa_start || gap_end == msa_end ||
            (gap_end - gap_start + 1) > 0.2 * (msa_end - msa_start + 1) ||
            (gap_end - gap_start + 1) % 3 != 0)
                                               
          return 0;
      }
    }
  }
  return 1;
}

/* filter at group level.  This version simply calls a group invalid
   if any of its feature are invalid */
int valid_group(GFF_Set *group, MSA *msa) {
  int i;
  for (i = 0; i < lst_size(group->features); i++)
    if (!valid_feat(lst_get_ptr(group->features, i), msa))
      return 0;
  return 1;
}

/* filter a set of predictions according to 'valid_group' */
/* void filter_predictions(GFF_Set *gff, MSA *msa) { */
/*   int i, j; */
/*   List *groups = lst_new_ptr(max(lst_size(gff->features) / 10, 5)); */
/*   List *newfeats = lst_new_ptr(lst_size(gff->features)); */
/*   gff_partition_by_group(gff, groups); */
/*   for (i = 0; i < lst_size(groups); i++) { */
/*     GFF_Set *g = lst_get_ptr(groups, i); */
/*     int valid = valid_group(g, msa); */
/*     for (j = 0; j < lst_size(g->features); j++) { */
/*       GFF_Feature *f = lst_get_ptr(g->features, j); */
/*       if (!valid) */
/*         str_append_char(f->feature, 'X'); */
/*       lst_push_ptr(newfeats, f);       */
/*     } */
/*     lst_clear(g->features); */
/*     gff_free_set(g); */
/*   } */
/*   for (i = 0; i < lst_size(gff->features); i++) */
/*     gff_free_feature(lst_get_ptr(gff->features, i)); */
/*   lst_free(gff->features); */
/*   gff->features = newfeats; */
/*   lst_free(groups); */
/* } */

/* multiply all transition probabilities from designated background
   categories to non-background categories by the specified factor,
   then renormalize.  This provides a simple "knob" for the
   sensitivity/specificity tradeoff.  Mathematically, it's a way of
   changing the expected number of predicted features.  Category 0 is
   assumed as a background category */
void scale_trans_from_backgd(PhyloHMM_Puzzler *puz, List *backgd_cat_names, 
                             double scale_factor) {
  /* FIXME: check cat nos */
  double is_backgd_cat[puz->cm->ncats+1];
  int i, j;

  is_backgd_cat[0] = 1;
  for (i = 1; i <= puz->cm->ncats; i++) is_backgd_cat[i] = 0;
  if (backgd_cat_names != NULL) {
    List *backgd_cat_nos = cm_get_category_list(puz->cm, backgd_cat_names, 0); 
    for (i = 0; i < lst_size(backgd_cat_nos); i++) 
      is_backgd_cat[lst_get_int(backgd_cat_nos, i)] = 1;
    lst_free(backgd_cat_nos);
  }

  for (i = 0; i < puz->hmm->nstates; i++)
    if (is_backgd_cat[puz->state_to_cat[i]])
      for (j = 0; j < puz->hmm->nstates; j++)
        if (!is_backgd_cat[puz->state_to_cat[j]] && 
            mm_get(puz->hmm->transition_matrix, i, j) != 0) 
          mm_set(puz->hmm->transition_matrix, i, j, 
                 scale_factor * 
                 mm_get(puz->hmm->transition_matrix, i, j));

  hmm_renormalize(puz->hmm);
}

int main(int argc, char* argv[]) {
  FILE* F;
  MSA *msa, *msa_compl = NULL;
  TreeModel **mod;
  double **posterior_probs, **forward_scores, **emissions;
/* *exp_vals_mod, *exp_vals_tot,     *smoothed_exp_vals, ; */
  int *path, *state_pos, *state_neg, *msa_gap_patterns = NULL,
    *missing = NULL;
  HMM *hmm = NULL;
  TreeNode *tree = NULL;
  int i, j, input_format = FASTA,
    likelihood_only = 0, msa_idx, quiet_mode = 0, gff_mode = 0, 
    reverse_complement = 0, ncats, nmsas,
    ncats_unspooled, nratecats = 1, autorates = 0, 
    gap_strip_mode = NO_STRIP, explicit_nratecats = -1,
    reflect_hmm = 0, sens_spec_mode = 0, score_mode = 0, indel_nseqs = -1,
    postprob_ref = 0, hide_missing = 0;
  double lambda = -1, loglikelihood = 0, mult_trans_probs = -1;
  char *catmap_fname = NULL, *hmm_fname = NULL, *hmm_train_fname = NULL,
    *pseudocount_matrix_fname = NULL, *output_fname_root = "label-out",
    *apriori_label_fname = NULL, *source_name = "chr?", *tree_fname = NULL;
  String *msa_fname, *gff_fname;
  List *gff_fname_list = NULL, *msa_fname_list = NULL, 
    *model_fname_list = NULL, *scaling_const_list = NULL, 
    *fixed_labels_list = NULL, /* *em_training_msas = NULL,  */
    *msa_length_list = NULL, *bed_feature_list = NULL, 
    *backgd_cat_names = NULL, *no_gaps_str = NULL, *model_indels_str = NULL,
    *penalize_gaps_str = NULL, *samples_states_str = NULL, 
    *postprob_states_str = NULL, *gc_thresholds = NULL;
  gsl_matrix *traincounts = NULL;
  gsl_vector *begcounts = NULL, *statecounts = NULL;
  CategoryMap *cm = NULL;
  char c;
  mode_type mode = TRAINING;
  String *outfname = str_new(STR_SHORT_LEN);
  PhyloHMM_Puzzler *puz;
  GapPatternMap *gpm = NULL;

  while ((c = getopt(argc, argv, "p:i:g:c:H:d:m:t:lr:L:C:P:F:G:o:M:b:s:k:V:I:W:Z:T:w:X:n:Q:j:D:xSYUEARfhq")) != -1) {
    switch(c) {
    case 'p':
      postprob_states_str = get_arg_list(optarg);
      break;
    case 'Q':
      samples_states_str = get_arg_list(optarg);
      break;
    case 'j':
      postprob_ref = atoi(optarg);
      break;
    case 'x':
      hide_missing = 1;
      break;
    case 'i':
      input_format = msa_str_to_format(optarg);
      if (input_format == -1) die("ERROR: bad alignment format.\n");
      break;
    case 'g':
      gff_fname_list = get_arg_list(optarg);
      break;
    case 'c':
      catmap_fname = optarg;
      break;
    case 'm':
      msa_fname_list = get_arg_list(optarg);
      break;
    case 'd':
      model_fname_list = get_arg_list(optarg);
      break;
    case 'H':
      hmm_fname = optarg;
      break;
    case 't':
      hmm_train_fname = optarg;
      break;
    case 'l':
      likelihood_only = 1;
      break;
    case 'r':
      scaling_const_list = get_arg_list(optarg);
      break;
    case 'L':
      lambda = atof(optarg);
      break;
    case 'A':
      autorates = 1;
      break;
    case 'k':
      explicit_nratecats = atoi(optarg);
      break;
    case 'C':
      apriori_label_fname = optarg;
      break;
    case 'M':
      msa_length_list = str_list_as_int(get_arg_list(optarg));
      break;
    case 'R':
      reverse_complement = 1;
      break;
    case 'U':
      reflect_hmm = 1;
      break;
    case 'V':
      backgd_cat_names = get_arg_list(optarg);
      break;
    case 'S':
      score_mode = 1;
      break;
    case 'X':
      mult_trans_probs = atof(optarg);
      break;
    case 'Y':
      sens_spec_mode = 1;
      break;
    case 'W':
      no_gaps_str = get_arg_list(optarg);
      break;
    case 'Z':
      penalize_gaps_str = get_arg_list(optarg);
      break;
    case 'I':
      model_indels_str = get_arg_list(optarg);
      break;
    case 'n':
      indel_nseqs = atoi(optarg);
      break;
    case 'T':
      tree_fname = optarg;
      break;
    case 'P':
      pseudocount_matrix_fname = optarg;
      break;
    case 'F':
      fixed_labels_list = get_arg_list(optarg);
      break;
    case 'G':
      if (!strcmp(optarg, "ALL")) gap_strip_mode = STRIP_ALL_GAPS;
      else if (!strcmp(optarg, "ANY")) gap_strip_mode = STRIP_ANY_GAPS;
      else gap_strip_mode = atoi(optarg);
      break;
    case 'o':
      output_fname_root = optarg;
      break;
    case 'f':
      gff_mode = 1;
      break;
    case 'b':
      bed_feature_list = get_arg_list(optarg);
      break;
    case 'D':
      gc_thresholds = str_list_as_dbl(get_arg_list(optarg));
      for (i = 0; i < lst_size(gc_thresholds); i++)
        if (lst_get_dbl(gc_thresholds, i) <= 0 || 
            lst_get_dbl(gc_thresholds, i) >= 1) 
          die("ERROR: Bad G+C threshold.\n");
      lst_qsort_dbl(gc_thresholds, ASCENDING);
      break; 
    case 's':
      source_name = optarg;
      break;
    case 'q':
      quiet_mode = 1;
      break;
    case 'h':
      print_usage();
      exit(0);
    case '?':
      fprintf(stderr, "ERROR: unrecognized option.\n\nType label -h for usage.\n");
      exit(1);
    }
  }

  /* Check validity of arguments */
  if (gff_fname_list == NULL && hmm_train_fname == NULL)
    mode = TESTING;             /* otherwise remain TRAINING */

  if (msa_fname_list == NULL && (mode != TRAINING || msa_length_list == NULL)) {
    fprintf(stderr, "ERROR: -m is required.\n\nType label -h for usage.\n");
    exit(1);
  }

  if (mode == TRAINING) {
    if (gff_fname_list == NULL || hmm_train_fname == NULL) {
        fprintf(stderr, "ERROR: -g and -t required in training mode.  Type label -h for usage.\n");
        exit(1);
    }
    if (hmm_fname != NULL || model_fname_list != NULL || 
        postprob_states_str != NULL || samples_states_str ||
        likelihood_only == 1 || fixed_labels_list != NULL || 
        scaling_const_list != NULL || lambda != -1 || reflect_hmm ||
        no_gaps_str != NULL || penalize_gaps_str != NULL) {
      fprintf(stderr, "ERROR: -H, -d, -p, -l, -F, -r, -L, -Q, -U, -W, and -Z not valid in training mode.\n\nType label -h for usage.\n");
      exit(1);
    }
    if (msa_length_list != NULL && msa_fname_list != NULL) {
      fprintf(stderr, "ERROR: -m and -M are mutually exclusive.  Type label -h for usage.\n");
      exit(1);
    }
    if (model_indels_str != NULL && tree_fname == NULL) {
      fprintf(stderr, "ERROR: -I requires -T in training mode.  Type label -h for usage.\n");
      exit(1);
    }
  }
  else if (mode == TESTING) {
    if (lst_size(msa_fname_list) != 1) {
      fprintf(stderr, "ERROR: only one MSA allowed in testing mode.\n");
      exit(1);
    }
    if (model_fname_list == NULL) {
      fprintf(stderr, "ERROR: -d required in testing mode.\n\nType label -h for usage.\n");
      exit(1);
    }
    if (hmm_train_fname != NULL || pseudocount_matrix_fname != NULL || 
        msa_length_list != NULL) {
      fprintf(stderr, "ERROR: -t, -P, and -M not valid in testing mode.\n\nType label -h for usage.\n");
      exit(1);
    }
    if (lambda != -1 && scaling_const_list == NULL && autorates == 0) {
      fprintf(stderr, "ERROR: -L requires -r or -A.\n\nType label -h for usage.\n");
      exit(1);
    }
    if (autorates == 1 && scaling_const_list != NULL) {
      fprintf(stderr, "ERROR: -r and -A are mutually exclusive.\n\nType label -h for usage.\n");
      exit(1);
    }
    if (reverse_complement && reflect_hmm) {
      fprintf(stderr, "ERROR: -U and -R are mutually exclusive.\n\nType label -h for usage.\n");
      exit(1);
    }
    if (lambda != -1 && (lambda < 0 || lambda > 1)) {
      fprintf(stderr, "ERROR: argument of -L must be between 0 and 1 (inclusive.\n\nType label -h for usage.\n");
      exit(1);
    }
  }

  /* General set-up */
  
  /* initialize category map */
  if (catmap_fname == NULL) {
    if (!quiet_mode)
      fprintf(stderr, "Creating trivial (single-category) category map ...\n");
    cm = cm_new(0);
    cm->ranges[0] = 
      cm_new_category_range(str_new_charstr(BACKGD_CAT_NAME), 0, 0);
  }
  else {
    if (!quiet_mode)
      fprintf(stderr, "Reading category map from %s ...\n", catmap_fname);

    if ((F = fopen(catmap_fname, "r")) == NULL ||
        (cm = cm_read(F)) == NULL) {
      fprintf(stderr, "ERROR reading from %s.\n", catmap_fname);
      exit(1);
    }      
    fclose(F);
  }
  
  ncats = cm->ncats + 1;
  ncats_unspooled = cm->unspooler != NULL ? cm->unspooler->nstates_unspooled : 
    ncats;
  nmsas = (msa_length_list != NULL ? lst_size(msa_length_list) : 
           lst_size(msa_fname_list));

  if (mode == TRAINING) {
    if (model_indels_str != NULL) {
      tree = parse_nh_from_file(fopen_fname(tree_fname, "r"));
      if (tree == NULL) { 
        fprintf(stderr, "ERROR: invalid tree in %s.\n", tree_fname); 
        exit (1); 
      }
      gpm = gp_create_gapcats(cm, model_indels_str, 
                              indel_nseqs > 0 ? indel_nseqs :
                              (tree->nnodes+1)/2); /* note: no alignment yet */
      ncats = cm->ncats + 1;    /* numbers will change */ 
      ncats_unspooled = cm->unspooler == NULL ? ncats : 
        cm->unspooler->nstates_unspooled;
    }

    if (hmm_train_fname != NULL) {
      /* allocate memory for storage of "training paths" */
      traincounts = gsl_matrix_calloc(ncats_unspooled, ncats_unspooled);
      statecounts = gsl_vector_calloc(ncats_unspooled);
      begcounts = gsl_vector_calloc(ncats_unspooled);
    
      /* create skeleton of new HMM. */
      hmm = hmm_new_nstates(ncats_unspooled, 0, 0);
    }
  }
  else {                        /* mode == TESTING */

    double *scaling_consts = NULL;
    int nrc;

    if (likelihood_only || fixed_labels_list != NULL || 
        apriori_label_fname != NULL) {
      if (likelihood_only && fixed_labels_list != NULL) {
        fprintf(stderr, "ERROR: cannot specify both -l and -F.\n\nType label -h for usage.\n");
        exit(1);
      }
      if (fixed_labels_list != NULL && 
          lst_size(fixed_labels_list) != lst_size(msa_fname_list)) {
        fprintf(stderr, "ERROR: list of fixed-labels files (-F) must be equal in length to list of MSA filenames.\n");
        exit(1);
      }
    }

    /* read HMM */
    if (hmm_fname == NULL) {
      if (!quiet_mode)
        fprintf(stderr, "Creating trivial (single-state) HMM ...\n");
      hmm = hmm_create_trivial();
    }
    else {
      if ((F = fopen(hmm_fname, "r")) == NULL) {        fprintf(stderr, "ERROR: cannot open %s.\n", hmm_fname);
        exit(1);
      }      

      if (!quiet_mode)
        fprintf(stderr, "Reading HMM from %s ...\n", hmm_fname);

      hmm = hmm_new_from_file(F); 
    }

    /* read alignment */
    msa_fname = (String*)lst_get_ptr(msa_fname_list, 0);
    if (!quiet_mode)
      fprintf(stderr, "Reading alignment from %s ...\n", 
              str_equals_charstr(msa_fname, "-") ? "stdin" :
              msa_fname->chars);
    msa = msa_new_from_file(fopen_fname(msa_fname->chars, "r"), 
                            input_format, "ACGT");

    if (input_format == SS) {
      if (msa->ss->tuple_idx == NULL) {
        fprintf(stderr, "ERROR: Ordered representation of alignment required.\n");
        exit(1);
      }

      /* temporary: remove 'N' from alphabet, in this case */
      for (i = 0, j = 0; i < strlen(msa->alphabet); i++) {
        if (msa->alphabet[i] != 'N') {
          msa->alphabet[j++] = msa->alphabet[i];
          msa->inv_alphabet[(int)'N'] = -1;
        }
      }
      msa->alphabet[j] = '\0';
    }

    if (gc_thresholds != NULL) {
      double gc;
      String *models_fname;
      char tmpstr[STR_SHORT_LEN];
      gsl_vector *f = msa_get_base_freqs(msa, -1, -1); 
                                /* note: does not consider categories */
      assert(msa->inv_alphabet[(int)'G'] >= 0 && 
             msa->inv_alphabet[(int)'C'] >= 0);
      gc = gsl_vector_get(f, msa->inv_alphabet[(int)'G']) +
        gsl_vector_get(f, msa->inv_alphabet[(int)'C']);
      assert(lst_size(model_fname_list) == lst_size(gc_thresholds) + 1);
      for (i = 0; models_fname == NULL && i < lst_size(gc_thresholds); i++)
        if (gc < lst_get_dbl(gc_thresholds, i))
          models_fname = lst_get_ptr(model_fname_list, i);
      if (models_fname == NULL) 
        models_fname = lst_get_ptr(model_fname_list, i); /* last element */

      if (!quiet_mode) 
        fprintf(stderr, "G+C content is %.1f%%; using models for partition %d (%s) ...\n", gc*100, i, models_fname->chars);

      /* this trick makes it as if the correct set of models had been
         specified directly */
      sprintf(tmpstr, "*%s", models_fname->chars);
      lst_free_strings(model_fname_list); lst_clear(model_fname_list);
      model_fname_list = get_arg_list(tmpstr);

      gsl_vector_free(f);
    }

    /* read tree models */
    if (lst_size(model_fname_list) != ncats) {
      fprintf(stderr, "ERROR: number of tree models must equal number of site categories.\n");
      exit(1);
    }
    
    mod = (TreeModel**)smalloc(sizeof(TreeModel*) * ncats);
    for (i = 0, j = 0; i < ncats; i++) {
      String *fname = (String*)lst_get_ptr(model_fname_list, i);
      if (!quiet_mode)
        fprintf(stderr, "Reading tree model from %s ...\n", fname->chars);
      if ((F = fopen(fname->chars, "r")) == NULL) {
        fprintf(stderr, "ERROR: cannot open %s.\n", fname->chars);
        exit(1);
      }
      mod[i] = tm_new_from_file(F); 
      mod[i]->use_conditionals = 1;
      if (mod[i]->tree != NULL) {
        if (tree == NULL)
          tree = mod[i]->tree;
        else if (mod[i]->tree->nnodes != tree->nnodes) {
          fprintf(stderr, "ERROR: trees of models have different numbers of nodes.\n");
          exit(1);
        }
      }
      fclose(F);
      j++;
    }

    if (model_indels_str != NULL && tree == NULL) { 
      /* possible if all models are weight matrices */
      if (tree_fname == NULL || 
          (tree = parse_nh_from_file(fopen_fname(tree_fname, "r"))) == NULL) { 
        fprintf(stderr, "ERROR: indel mode requires a tree; you have not specified one, or you have specified one that is invalid.\n"); 
        exit (1); 
      }
    }

    /* set up gap cats */
    if (no_gaps_str != NULL) {
      List *l = cm_get_category_list(cm, no_gaps_str, 0);
      for (i = 0; i < lst_size(l); i++) 
        mod[lst_get_int(l, i)]->allow_gaps = 0;
      lst_free(l); lst_free_strings(no_gaps_str); lst_free(no_gaps_str);
    }      
    else if (penalize_gaps_str != NULL) {
      List *l = cm_get_category_list(cm, penalize_gaps_str, 0);
      for (i = 0; i < lst_size(l); i++) 
        mod[lst_get_int(l, i)]->allow_but_penalize_gaps = 1;
      lst_free(l); lst_free_strings(penalize_gaps_str); 
      lst_free(penalize_gaps_str);
    }      

    if (tree != NULL && msa->nseqs * 2 - 1 != tree->nnodes) {
      fprintf(stderr, "ERROR: number of leaves in trees must equal number of sequences in alignment.\n");
      exit(1);
    }

    if (gap_strip_mode != NO_STRIP)
      strip_gaps(msa, gap_strip_mode);

    if (apriori_label_fname != NULL) {
      GFF_Set *gff;
      if (!quiet_mode)
        fprintf(stderr, "Reading annotations from %s ...\n", 
                apriori_label_fname);

      if ((F = fopen(apriori_label_fname, "r")) == NULL || 
          (gff = gff_read_set(F)) == NULL) {
        fprintf(stderr, "ERROR reading from %s.\n", apriori_label_fname);
        exit(1);
      }
      fclose(F);

      if (reverse_complement) {
        gff_group(gff, "exon_id");
        msa_reverse_compl_feats(msa, gff, NULL);
      }

      msa_label_categories(msa, gff, cm);

      gff_free_set(gff);
    }    
    else if (reverse_complement) {
      if (!quiet_mode)
        fprintf(stderr, "Reverse complementing alignment ...\n");
      msa_reverse_compl(msa);        
    }
    else if (reflect_hmm) {
      int idx1, idx2, tupsize = -1;

      msa_compl = msa_create_copy(msa, 0);

      /* temporary -- work-around for problem with context being wrong
         in suff stats at boundaries of MAF blocks; reverse complement
         using the complete alignment, not just the suff stats */
      if (msa_compl->ss != NULL && msa_compl->ss->tuple_size > 1) {
        tupsize = msa_compl->ss->tuple_size;
        if (msa_compl->seqs == NULL) ss_to_msa(msa_compl);
        ss_free(msa_compl->ss);
        msa_compl->ss = NULL;
      }
      /* end temporary */

      msa_reverse_compl(msa_compl);

      /* temporary */
      if (tupsize != -1) {
        ss_from_msas(msa_compl, tupsize, 1, NULL, NULL, NULL, -1);
        /* get rid of the sequences! they'll be wrong! */
        for (i = 0; i < msa_compl->nseqs; i++) free(msa_compl->seqs[i]);
        free(msa_compl->seqs);
        msa_compl->seqs = NULL;
      }
      /* end temporary */
      
      /* we want the reverse complemented alignment, but the indexing
         of the forward strand; to save code, we'll just *reverse* the
         reverse complement */
      for (idx1 = 0, idx2 = msa_compl->length-1; 
           idx1 < idx2; idx1++, idx2--) {
        int tmp = msa_compl->ss->tuple_idx[idx2];
        msa_compl->ss->tuple_idx[idx2] = msa_compl->ss->tuple_idx[idx1];
        msa_compl->ss->tuple_idx[idx1] = tmp;
      }
    }

    /* if modeling indels, obtain gap pattern for each site in the alignment */
    if (model_indels_str != NULL) {
      if (!quiet_mode)
        fprintf(stderr, "Obtaining gap patterns ...\n");
      msa_gap_patterns = smalloc(msa->length * sizeof(int));
      gp_set_phylo_patterns(msa_gap_patterns, msa, tree);
    }

    /* set up puzzler */
    nrc = 1;
    if (autorates) nrc = explicit_nratecats;
    else if (scaling_const_list != NULL) { 
      nrc = lst_size(scaling_const_list);
      scaling_consts = smalloc(nrc * sizeof(double));
      for (i = 0; i < nratecats; i++) {
        if (str_as_dbl((String*)lst_get_ptr(scaling_const_list, i), 
                       &scaling_consts[i]) != 0) {
          fprintf(stderr, "ERROR: cannot parse number in rate-constant list.\n"); 
          exit(1);
        }
      }
    }
    puz = puz_new(cm, mod, hmm, reflect_hmm, backgd_cat_names, nrc, 
                  scaling_consts, lambda, model_indels_str, msa->nseqs);
    if (scaling_consts != NULL) free(scaling_consts);
    hmm = NULL;                 /* must use puzzler's hmm now */

    /* multiply transition probs, if necessary */
    if (mult_trans_probs > 0) 
      scale_trans_from_backgd(puz, backgd_cat_names, mult_trans_probs);

    /* now compute emissions */
    emissions = smalloc(puz->hmm->nstates * sizeof(double*));

    /* set up mapping from model/strand to first associated state
       (allows emissions to be computed only once for each
       model/strand pair) */
    state_pos = smalloc(puz->nmods * sizeof(int));
    state_neg = smalloc(puz->nmods * sizeof(int));
    for (i = 0; i < puz->nmods; i++) state_pos[i] = state_neg[i] = -1;

    for (i = 0; i < puz->hmm->nstates; i++) {
      if (!quiet_mode)
        fprintf(stderr, "Computing emission probs (state %d,  cat %d, mod %d, pattern %d, strand %c) ...\n", i, puz->state_to_cat[i], puz->state_to_mod[i], puz->state_to_pattern[i], puz->reverse_compl[i] ? '-' : '+');

      if (apriori_label_fname != NULL) {
        emissions[i] = smalloc(msa->length * sizeof(double));
        loglikelihood += 
          tl_compute_log_likelihood(puz->mods[puz->state_to_mod[i]], 
                                    puz->reverse_compl[i] ? msa : msa_compl,
                                    emissions[i], puz->state_to_cat[i], NULL);
      }
      else {
        /* reuse already computed values if possible */
        int mod = puz->state_to_mod[i];
        if (!puz->reverse_compl[i] && state_pos[mod] != -1)
          emissions[i] = emissions[state_pos[mod]]; /* saves memory */
        else if (puz->reverse_compl[i] && state_neg[mod] != -1)
          emissions[i] = emissions[state_neg[mod]];
        else {
          emissions[i] = smalloc(msa->length * sizeof(double));
          tl_compute_log_likelihood(puz->mods[mod], 
                                    puz->reverse_compl[i] ? msa_compl : msa,
                                    emissions[i], -1, NULL);
          if (!puz->reverse_compl[i]) state_pos[mod] = i;
          else state_neg[mod] = i;            
        }
      }
    }

    /* redefine emissions for indel states, if necessary */
    if (model_indels_str != NULL) {
      if (!quiet_mode)
        fprintf(stderr, "Adjusting emission probs according to gap patterns ...\n");
      for (i = puz->hmm->nstates - 1; i >= 0; i--) {
                                /* by going backwards, we ensure that
                                   the "base" state (with gap pattern
                                   == 0) is visited last (see
                                   puzzler.c and gap_patterns.c) */
        if (puz->state_to_pattern[i] >= 0) {
          double *orig_emissions = emissions[i];
          if (puz->state_to_pattern[i] > 0)
            emissions[i] = smalloc(msa->length * sizeof(double));
                                /* otherwise, use the array already
                                   allocated */
          for (j = 0; j < msa->length; j++) 
            emissions[i][j] = 
              (msa_gap_patterns[j] == puz->state_to_pattern[i] ? 
              orig_emissions[j] : NEGINFTY);
                                 
        }
      }

      /* in addition, states corresponding to complex gap patterns must have a fixed, small probability of emission for all columns matching their pattern */
      /* WAIT: is this necessary? */
/*       for (j = 0; j < msa->length; j++) { */
/*         if (msa_gap_patterns[j] == puz->complex_gap_pattern)  */
/*           for (i = 0; i < puz->hmm->nstates */
/*       } */
    }

    /* fit lambda, if necessary */
    if (autorates == 1 && lambda == -1) {
      if (!quiet_mode) fprintf(stderr, "Finding MLE for lambda ...");
      lambda = fit_lambda(puz, emissions, msa->length);
      if (!quiet_mode) fprintf(stderr, "lambda = %f\n", lambda);
      puz_update_cross_prod(puz, lambda);
    }
    
  } /* end testing mode stuff */

  /* Main loop: consider each MSA in turn */
  for (msa_idx = 0; msa_idx < nmsas; msa_idx++) {
    if (msa_fname_list != NULL) {
      msa_fname = (String*)lst_get_ptr(msa_fname_list, msa_idx);

      if (mode == TRAINING) {
        if (str_equals_charstr(msa_fname, "-")) F = stdin;
        else if ((F = fopen(msa_fname->chars, "r")) == NULL) {
          fprintf(stderr, "ERROR: cannot open %s.\n", msa_fname->chars);
          exit(1);
        }
        if (!quiet_mode)
          fprintf(stderr, "Reading alignment from %s ...\n", 
                  F == stdin ? "stdin" : msa_fname->chars);
        msa = msa_new_from_file(F, input_format, NULL);
        fclose(F);

        if (gap_strip_mode != NO_STRIP)
          strip_gaps(msa, gap_strip_mode);
      }
    }
    else {                      /* only lengths of alignments specified */
      msa = msa_new(NULL, NULL, 0, lst_get_int(msa_length_list, msa_idx), NULL);
                                /* just a shell in this case */
    }

    if (mode == TRAINING) {
      GFF_Set *gff;

      gff_fname = (String*)lst_get_ptr(gff_fname_list, msa_idx);

      if (!quiet_mode)
        fprintf(stderr, "Reading annotations from %s ...\n", 
                gff_fname->chars);

      if ((F = fopen(gff_fname->chars, "r")) == NULL || 
          (gff = gff_read_set(F)) == NULL) {
        fprintf(stderr, "ERROR reading from %s.\n", gff_fname->chars);
        exit(1);
      }
      fclose(F);

      /* convert GFF to coordinate frame of alignment */
      if (msa_length_list == NULL) {
        if (!quiet_mode)
          fprintf(stderr, "Mapping annotations to alignment ...\n");
        msa_map_gff_coords(msa, gff, -1, 0, 0, NULL);
      }

/*       if (!quiet_mode)  */
/*         fprintf(stderr, "Removing overlapping annotations ...\n"); */
/*       gff_remove_overlaps_by_group(gff); */
      
      /* NOTE: be sure nonoverlapping input (put in USAGE) */

      if (model_indels_str != NULL) {
        if (!quiet_mode)
          fprintf(stderr, "Obtaining gap patterns ...\n");
        msa_gap_patterns = smalloc(msa->length * sizeof(int));
        gp_set_phylo_patterns(msa_gap_patterns, msa, tree);
      }

      /* at this point, we don't actually need the alignment anymore;
         if using ordered suff stats (likely with large data sets),
         can free them now, to avoid running out of memory */
      if (msa->ss != NULL) { ss_free(msa->ss); msa->ss = NULL; }

      if (reverse_complement) {
        if (!quiet_mode)
          fprintf(stderr, "Reverse complementing features on negative strand ...\n");
        /* we don't need to reverse complement the whole alignment --
           just the gff and possibly the gap pattern array (pass a
           NULL msa) */        
        gff_group(gff, "exon_id"); /* FIXME: tag? */
        msa_reverse_compl_feats(NULL, gff, msa_gap_patterns);
      }

      if (!quiet_mode)
        fprintf(stderr, "Labeling sites by category ...\n");       
      msa_label_categories(msa, gff, cm);

      gff_free_set(gff);

      if (hmm_train_fname != NULL) { 
        if (model_indels_str != NULL) {
          if (!quiet_mode)
            fprintf(stderr, "Remapping categories according to gap patterns ...\n");

          if (indel_nseqs > 0 && indel_nseqs != msa->nseqs) {
            /* in this case, we'll simply reassign non-trivial gap
               patterns randomly.  This will achieve the desired
               effect with minimal coding, as long as the number of
               sites is not too small (the indel model is probably
               useless anyway if the number is small) */
            int pat, newpat;
            int npatterns = 4 * indel_nseqs - 5;
            int complex_allowed[cm->ncats+1];
            List *no_complex_names, *no_complex_nums;

            if (!quiet_mode)
              fprintf(stderr, "(target number of sequences: %d)\n", indel_nseqs);

            /* set up index indicating by cat no. whether complex gaps
               are allowed */
            for (i = 0; i < ncats; i++) complex_allowed[i] = 1;
            no_complex_names = lst_new_ptr(10);
            str_split(str_new_charstr(NO_COMPLEX), ",", no_complex_names);
            no_complex_nums = cm_get_category_list(cm, no_complex_names, 1);
            for (i = 0; i < lst_size(no_complex_nums); i++)
              complex_allowed[lst_get_int(no_complex_nums, i)] = 0;
            lst_free(no_complex_nums);
            lst_free_strings(no_complex_names); lst_free(no_complex_names);

            /* now reassign all non-null numbers */
            for (i = 0; i < msa->length; ) {
              if ((pat = msa_gap_patterns[i]) != 0) {
                if (complex_allowed[msa->categories[i]])
                  newpat = 1 + (random() * (double)npatterns / (RAND_MAX + 1.0));
                                /* random number in interval [1, npatterns] */
                else 
                  newpat = 1 + (random() * (double)(npatterns-1) / (RAND_MAX + 1.0));
                                /* random number in interval [1,npatterns-1] 
                                   (excludes complex gap pattern) */
                for (; i < msa->length && msa_gap_patterns[i] == pat; i++)
                  msa_gap_patterns[i] = newpat; /* change for whole sequence */
              }
              else i++;
            }
          }

          /* obtain gapped category number for each site */
          for (i = 0; i < msa->length; i++) 
            if (gpm->cat_x_pattern_to_gapcat[msa->categories[i]] != NULL)
              msa->categories[i] = gpm->cat_x_pattern_to_gapcat[msa->categories[i]][msa_gap_patterns[i]];
        }

        if (!quiet_mode)
          fprintf(stderr, "Unspooling categories ...\n");
        cm_spooled_to_unspooled(cm, msa->categories, msa->length);

        if (!quiet_mode)
          fprintf(stderr, "Collecting training data ...\n");
        hmm_train_update_counts(traincounts, statecounts, begcounts, 
                                msa->categories, msa->length, 
                                ncats_unspooled);
      }
    }

    else {                      /* mode == TESTING */
      int try, ntries;
      char *orig_output_fname_root;
      char tmpstr[STR_MED_LEN];
      

      /* NOTE: matrix of emission scores will have already been
         created */

      /* TEMPORARY!!! */
      ntries = 1;
      if (sens_spec_mode) {    /* FIXME: add param  */
        scale_trans_from_backgd(puz, backgd_cat_names, exp(SCALE_RANGE_MIN));
        ntries = NSENS_SPEC_VERSIONS;
        orig_output_fname_root = output_fname_root;
        output_fname_root = tmpstr;
      }

      for (try = 0; try < ntries; try++) {
        
        if (ntries > 1) {
          if (!quiet_mode)
            fprintf(stderr, "(Sens-spec tradeoff #%d)\n", try+1);
          sprintf(output_fname_root, "%s.v%d", orig_output_fname_root, try+1);
        }

        /* run Viterbi */
        if (gff_mode || bed_feature_list != NULL) {
          if (!quiet_mode)
            fprintf(stderr, "Executing Viterbi algorithm ...\n");
          path = (int*)smalloc(msa->length * sizeof(int));
          hmm_viterbi(puz->hmm, emissions, msa->length, path);
        }

        /* obtain posterior probs */
        if (postprob_states_str != NULL || samples_states_str) {

          posterior_probs = (double**)smalloc(puz->hmm->nstates * sizeof(double*));

          /* only allocate memory for states of interest; NULLs for
             the others will cause hmm_posterior_probs to ignore
             them */
          if (postprob_states_str != NULL &&
              lst_size(postprob_states_str) == 1 && 
              str_equals_charstr(lst_get_ptr(postprob_states_str, 0), "all")) 
            for (i = 0; i < puz->hmm->nstates; i++) 
              posterior_probs[i] = (double*)smalloc(msa->length * sizeof(double));
          else {
            int state;
            for (i = 0; i < puz->hmm->nstates; i++) posterior_probs[i] = NULL;
            if (postprob_states_str != NULL) {
              for (i = 0; i < lst_size(postprob_states_str); i++) {
                if (str_as_int(lst_get_ptr(postprob_states_str, i), &state) != 0 || 
                    state < 0 || state >= puz->hmm->nstates) {
                  fprintf(stderr, "ERROR: illegal state number (\"%s\")\n", 
                          ((String*)lst_get_ptr(postprob_states_str, i))->chars);
                  exit(1);
                }
                if (posterior_probs[i] == NULL)
                  posterior_probs[i] = (double*)smalloc(msa->length * sizeof(double));
              }
            }
            if (samples_states_str != NULL) {
              for (i = 0; i < lst_size(samples_states_str); i++) {
                if (str_as_int(lst_get_ptr(samples_states_str, i), &state) != 0 || 
                    state < 0 || state >= puz->hmm->nstates) {
                  fprintf(stderr, "ERROR: illegal state number (\"%s\")\n", 
                          ((String*)lst_get_ptr(samples_states_str, i))->chars);
                  exit(1);
                }
                if (posterior_probs[i] == NULL)
                  posterior_probs[i] = (double*)smalloc(msa->length * sizeof(double));
              }
            }
          }

          if (!quiet_mode)
            fprintf(stderr, "Computing posterior probabilities ...\n");
          hmm_posterior_probs(puz->hmm, emissions, msa->length, posterior_probs);
        }
    
        if (likelihood_only && apriori_label_fname == NULL) {
          forward_scores = (double**)smalloc(puz->hmm->nstates * sizeof(double*));
          for (i = 0; i < puz->hmm->nstates; i++)
            forward_scores[i] = (double*)smalloc(msa->length * sizeof(double));
          if (!quiet_mode)
            fprintf(stderr, "Computing total log likelihood ...\n");
          loglikelihood = hmm_forward(puz->hmm, emissions, msa->length, forward_scores);
        }
        if (fixed_labels_list != NULL) {
          String *labels_fname = (String*)lst_get_ptr(fixed_labels_list, msa_idx);
          if ((F = fopen(labels_fname->chars, "r")) == NULL || 
              msa_read_category_labels(msa, F) != 0) {
            fprintf(stderr, "ERROR: Unable to read category labels from %s.\n", 
                    labels_fname->chars);
            exit(1);
          }
          if (!quiet_mode)
            fprintf(stderr, "Computing likelihood of path defined in %s ...\n", 
                    labels_fname->chars);
          loglikelihood = hmm_path_likelihood(puz->hmm, emissions, msa->length, 
                                              msa->categories);
        }

        /* print output */
        if (likelihood_only || fixed_labels_list != NULL || 
            apriori_label_fname != NULL) {
          if (scaling_const_list != NULL || autorates == 1)
            printf("lambda = %f\n", lambda);
          printf("lnL = %.4f\n", loglikelihood * log(2)); /* output natural log */
        }

        /* check for missing data, if -x option */
        if (hide_missing && 
            (postprob_states_str != NULL || samples_states_str != NULL)) {
          missing = smalloc(msa->length * sizeof(int));
          msa_find_noaln(msa, postprob_ref, MIN_BLOCK_SIZE, missing);
        }

        if (postprob_states_str != NULL) {
          int k;

          str_cpy_charstr(outfname, output_fname_root);
          str_append_charstr(outfname, ".postprob");
          if ((F = fopen(outfname->chars, "w+")) == NULL) {
            fprintf(stderr, "ERROR: cannot write to %s.\n", outfname->chars);
            exit(1);
          }

          if (!quiet_mode)
            fprintf(stderr, "Writing output to %s ...\n", outfname->chars);

          for (j = 0, k = 0; j < msa->length; j++) {
            if (postprob_ref == 0 || 
                (msa->seqs != NULL && msa->seqs[postprob_ref-1][j] != GAP_CHAR) ||
                (ss_get_char_pos(msa, j, postprob_ref-1, 0)) != GAP_CHAR) {

              if (!hide_missing || !missing[j]) {
                fprintf(F, "%d\t", k + msa->idx_offset + 1);
                for (i = 0; i < puz->hmm->nstates; i++) 
                  if (posterior_probs[i] != NULL)
                    fprintf(F, "%.4f\t", posterior_probs[i][j]);
                fprintf(F, "\n");
              }
              k++;
            }
          }
          fclose(F);
        }
      
        if (samples_states_str != NULL) {
          for (i = 0; i < lst_size(samples_states_str); i++) {
            int state;
            str_as_int(lst_get_ptr(samples_states_str, i), &state);
                                /* at this point we know the state is valid */
            str_cpy_charstr(outfname, output_fname_root);
            str_append_charstr(outfname, ".");
            str_append_int(outfname, state);
            str_append_charstr(outfname, ".samples");
            
            if (!quiet_mode) 
              fprintf(stderr, "Writing output to %s ...\n", outfname->chars);

            F = fopen_fname(outfname->chars, "w+");
            msa_scores_as_samples(msa, F, posterior_probs[state], 
                                  source_name, "postprob", 1000, 0, 
                                  postprob_ref, 0);
                                /* for now, assume first seq is reference */
            fclose(F);
          }
        }

        if (missing != NULL) free(missing);

        if (gff_mode || bed_feature_list != NULL) {
          GFF_Set *gff;
          List *framelist = lst_new_ptr(3);
          
          String *group_root = str_new(STR_SHORT_LEN);

          if (!quiet_mode)
            fprintf(stderr, "Reconstructing GFF from category labels ...\n");

          if (reverse_complement) {
            int idx1 = 0, idx2 = msa->length - 1;
            for (; idx1 < idx2; idx1++, idx2--) {
              int tmp = path[idx2];
              path[idx2] = path[idx1];
              path[idx1] = tmp;
            }
          }

          str_split(str_new_charstr(FRAME_CATS), ",", framelist);

          str_cpy_charstr(group_root, output_fname_root);
          str_remove_path(group_root);
          gff = cm_labeling_as_gff(cm, mode == TRAINING ? msa->categories : path, 
                                   msa->length, puz->state_to_cat, puz->reverse_compl,
                                   source_name, 
                                   "PHAST", 
                                   reverse_complement && mode == TESTING ? '-' : '+', 
                                   framelist, group_root->chars);

          /* score predictions */
          if (score_mode) { 
            List *score_cats = lst_new_ptr(10);  
            List *helper_cats = lst_new_ptr(10); 
            str_split(str_new_charstr(SCORE_CATS), ",", score_cats);
            str_split(str_new_charstr(HELPER_CATS), ",", helper_cats);

            if (!quiet_mode)
              fprintf(stderr, "Scoring predictions ...\n");
            
            score_predictions(gff, score_cats, helper_cats, backgd_cat_names, 
                              emissions, cm, puz, 0 /* msa->idx_offset */);

        /* TEMPORARY -- validate predictions by postprocessing */
/*         if (!quiet_mode) */
/*           fprintf(stderr, "Filtering predictions ...\n"); */
/*         filter_predictions (gff, msa); */

            lst_free_strings(score_cats); lst_free_strings(helper_cats); 
            lst_free(score_cats); lst_free(helper_cats); 
          }

          /* we generally want the GFF to be in the coord frame of the
             reference sequence.  Note that we handle the idx_offset
             here, after scoring.  FIXME: parameterize this? */
          msa_map_gff_coords(msa, gff, 0, 1, msa->idx_offset, NULL);

          /* now output predictions */
          if (gff_mode) {
            str_cpy_charstr(outfname, output_fname_root);
            if (nmsas > 1) {
              str_append_charstr(outfname, ".");
              str_append_int(outfname, msa_idx + 1);
            }
            str_append_charstr(outfname, ".gff");
            if ((F = fopen(outfname->chars, "w+")) == NULL) {
              fprintf(stderr, "ERROR: cannot write to %s.\n", outfname->chars);
              exit(1);
            }

            if (!quiet_mode)
              fprintf(stderr, "Writing output to %s ...\n", outfname->chars);

            gff_print_set(F, gff);
            fclose(F);
          }
          if (bed_feature_list != NULL) {
            str_cpy_charstr(outfname, output_fname_root);
            if (nmsas > 1) {
              str_append_charstr(outfname, ".");
              str_append_int(outfname, msa_idx + 1);
            }
            str_append_charstr(outfname, ".bed");
            if ((F = fopen(outfname->chars, "w+")) == NULL) {
              fprintf(stderr, "ERROR: cannot write to %s.\n", outfname->chars);
              exit(1);
            }

            if (!quiet_mode)
              fprintf(stderr, "Writing output to %s ...\n", outfname->chars);

            gff_print_bed(F, gff, NULL, bed_feature_list); /* FIXME: need to group by "exon_id" or something similar */
            fclose(F);
          }
          gff_free_set(gff);
          lst_free_strings(framelist);
          lst_free(framelist);
          str_free(group_root);
        }

        if (ntries > 1)
          scale_trans_from_backgd(puz, backgd_cat_names, 
                                  exp((SCALE_RANGE_MAX - SCALE_RANGE_MIN)/
                                      (NSENS_SPEC_VERSIONS-1)));
      } /* ntries loop (sens_spec_mode) */

      /* free stuff specific to testing mode (note that there can only
         be one MSA in testing mode */
      if (gff_mode || bed_feature_list != NULL)
        free(path);
      if (postprob_states_str != NULL || samples_states_str != NULL) {
        for (i = 0; i < puz->hmm->nstates; i++) 
          if (posterior_probs[i] != NULL) free(posterior_probs[i]);
        free(posterior_probs);
      }
      if (likelihood_only && apriori_label_fname == NULL) {
        for (i = 0; i < puz->hmm->nstates; i++) free(forward_scores[i]);
        free(forward_scores);
      }
      for (i = 0; i < puz->hmm->nstates; i++) 
        if (state_pos[puz->state_to_mod[i]] == i ||
            state_neg[puz->state_to_mod[i]] == i || 
            puz->state_to_pattern[i] >= 0)
          free(emissions[i]);
      free(emissions); free(state_pos); free(state_neg);
      puz_free(puz); 
    }

    /* free memory allocated per msa */
    if (msa_gap_patterns != NULL) free(msa_gap_patterns);
    msa_free(msa);
  }

  /* if training an HMM, do so now, using cumulative data */
  if (hmm_train_fname != NULL) {

    hmm_train_from_counts(hmm, traincounts, NULL, statecounts, NULL, 
                          begcounts, NULL);


    /* if modeling indels, adjust begin transitions so probability is
       distributed among different "gap pattern" states that all
       correspond to the same ungapped state (category); this helps
       avoid problems that occur when training on a few large sequences
       (e.g., whole chromosomes) and then testing on many shorter ones */
    if (model_indels_str != NULL) {
      double tprob[gpm->ncats]; 
      int nst[gpm->ncats];  /* total prob and number of states per
                               spooled, ungapped category */ 
      for (i = 0; i < gpm->ncats; i++) tprob[i] = nst[i] = 0;
      for (i = 0; i < hmm->nstates; i++) {
        if (gsl_vector_get(hmm->begin_transitions, i) > 0) 
          /* have to go from unspooled space to spooled space, then to
             ungapped space (HMM states correspond to unspooled,
             gapped categories).  Note that states with nonzero begin
             probs shouldn't be conditioned on other states. */
          tprob[gpm->gapcat_to_cat[cm_unspooled_to_spooled_cat(cm, i)]] += 
            gsl_vector_get(hmm->begin_transitions, i);
        nst[gpm->gapcat_to_cat[cm_unspooled_to_spooled_cat(cm, i)]]++;
      }
      for (i = 0; i < hmm->nstates; i++) 
        if (tprob[gpm->gapcat_to_cat[cm_unspooled_to_spooled_cat(cm, i)]] > 0)
          gsl_vector_set(hmm->begin_transitions, i, 
                         tprob[gpm->gapcat_to_cat[cm_unspooled_to_spooled_cat(cm, i)]] / 
                         nst[gpm->gapcat_to_cat[cm_unspooled_to_spooled_cat(cm, i)]]);
                                /* (uniform prior) */
    }

    /* write trained HMM */
    if ((F = fopen(hmm_train_fname, "w+")) == NULL) {
      fprintf(stderr, "ERROR writing to %s.\n", hmm_train_fname);
      exit(1);
    }        
    if (!quiet_mode)
      fprintf(stderr, "Saving trained HMM to %s ...\n", hmm_train_fname);
    hmm_print(F, hmm);
    fclose(F);
  }

  /* free remaining memory */
  cm_free(cm);
  if (outfname != NULL) str_free(outfname);
  if (gff_fname_list != NULL) 
    { lst_free_strings(gff_fname_list); lst_free(gff_fname_list); }
  if (msa_fname_list != NULL) 
    { lst_free_strings(msa_fname_list); lst_free(msa_fname_list); }
  if (model_fname_list != NULL) 
    { lst_free_strings(model_fname_list); lst_free(model_fname_list); }
  if (scaling_const_list != NULL) 
    { lst_free_strings(scaling_const_list); lst_free(scaling_const_list); }
  if (fixed_labels_list != NULL)
    { lst_free_strings(fixed_labels_list); lst_free(fixed_labels_list); }
  if (msa_length_list != NULL) lst_free(msa_length_list);
  if (gpm != NULL) gp_free_map(gpm);

  if (!quiet_mode)
    fprintf(stderr, "Done.\n");

  return 0;
}
