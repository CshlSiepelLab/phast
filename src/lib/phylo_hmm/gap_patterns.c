/* $Id: gap_patterns.c,v 1.1.1.1 2004-06-03 22:43:12 acs Exp $
   Written by Adam Siepel, 2003
   Copyright 2003, Adam Siepel, University of California */

/* Code relating to the patterns of gaps and nongaps that occur at
   each site in an alignment, which are the basis of the phylo-hmm's
   indel model */

#include <gap_patterns.h>
#include <category_map.h>
#include <sufficient_stats.h>

/* expands the set of categories in the specified cat map such that
   each designated "indel category" (to be given by feature name) is
   replaced by ngap_patterns categories, where ngap_patterns is a
   function of the number of sequences to be considered (nseqs).  In
   addition, for category ranges of size greater than one, categories
   corresponding to non-empty gap patterns are "conditioned on"
   categories corresponding to gapless sites, and the cat maps's
   unspooler is redefined accordingly (allows for "cyclic" gaps, as in
   coding regions).  The end result will be that each category in the
   original *unspooled* space will be replaced by ngap_patterns *
   range_size categories in the new unspooled space.  Return value is
   an object that defines the key mappings between the original
   (spooled) categories and the new "gapped categories", defined as
   category x gap pattern pairs.  */
GapPatternMap *gp_create_gapcats(CategoryMap *cm, List *indel_cats, int nseqs) {
  int i, j, k, cat, gapcat, range_size;
  List *indel_cat_nos;
  int new_dependencies = 0;

  GapPatternMap *gpm = smalloc(sizeof(GapPatternMap));
  gpm->ncats = cm->ncats + 1;
  gpm->ngap_patterns = 2 * (2 * nseqs - 3) + 2;
                                /* phylogenetic gap patterns: two per
                                   branch in the unrooted tree plus
                                   one corresponding to no gaps plus
                                   one for "complex" gap patterns; see
                                   function gp_set_phylo_patterns
                                   below */

  /* this is a simple way to scan for the total number of categories
     implied by the given list of names, which we need to determine
     the total number of gapcats */
  {
    List *catnos = cm_get_category_list(cm, indel_cats, 0);
    gpm->ngap_cats = gpm->ncats + (gpm->ngap_patterns-1) * lst_size(catnos);
    lst_free(catnos);
  }

  /* init mapping from gap cats to (spooled) cats */
  gpm->gapcat_to_cat = smalloc(gpm->ngap_cats * sizeof(int)); 
  for (i = 0; i < gpm->ngap_cats; i++) gpm->gapcat_to_cat[i] = i;

  /* init mapping from gap cats to gap patterns */
  gpm->gapcat_to_pattern = smalloc(gpm->ngap_cats * sizeof(int));
  for (i = 0; i < gpm->ngap_cats; i++) gpm->gapcat_to_pattern[i] = -1;

  /* init mapping from categories x gap patterns to gap cats */
  gpm->cat_x_pattern_to_gapcat = smalloc(gpm->ncats * sizeof(int*));
  for (i = 0; i < gpm->ncats; i++) gpm->cat_x_pattern_to_gapcat[i] = NULL;

  /* the mapping to gap cats has to be independent of the order
     in which the indel cats have been specified; we'll sort by base
     category number */
  indel_cat_nos = lst_new_int(lst_size(indel_cats));
  for (i = 0; i < lst_size(indel_cats); i++)
    lst_push_int(indel_cat_nos, cm_get_category(cm, lst_get_ptr(indel_cats, i)));
  lst_qsort_int(indel_cat_nos, ASCENDING);

  /* now enlarge category map and fill in mappings */
  cm_realloc(cm, gpm->ngap_cats-1);
  gapcat = gpm->ncats;
  for (i = 0; i < lst_size(indel_cat_nos); i++) {
    String *catname;
    int cyclic = 0;

    cat = lst_get_int(indel_cat_nos, i);
    catname = cm_get_feature(cm, cat);
    range_size = cm->ranges[cat]->end_cat_no - cm->ranges[cat]->start_cat_no + 1;

    if (range_size > 1) {
      cyclic = 1;
      if (str_equals_charstr(catname, "5'splice") /* TEMPORARY!  FOR TESTING */
          || str_equals_charstr(catname, "3'splice") 
          || str_equals_charstr(catname, "prestart")) 
        cyclic = 0;
    }

    /* initialize per-gap-pattern mapping (non-NULL only for indel
       categories).  Also init mappings for pattern 0 (no gaps) */
    assert(gpm->cat_x_pattern_to_gapcat[cat] == NULL); /* check redund. */
    for (k = 0; k < range_size; k++) {
      gpm->cat_x_pattern_to_gapcat[cat+k] = 
        smalloc(gpm->ngap_patterns*sizeof(int));
      gpm->cat_x_pattern_to_gapcat[cat+k][0] = cat+k; 
      gpm->gapcat_to_pattern[cat+k] = 0;
    }

    /* now handle the non-trivial gap patterns, which correspond to new cats */
    for (j = 1; j < gpm->ngap_patterns; j++) {
      String *newname;

      /* set up new category range */
      newname = str_dup(catname);
      str_append_charstr(newname, "_gp");
      str_append_int(newname, j);
      cm->ranges[gapcat] = cm_new_category_range(newname, gapcat, 
                                                 gapcat + range_size - 1);

      if (cyclic) {

        cm->conditioned_on[gapcat] = lst_new_int(range_size);
        for (k = 0; k < range_size; k++) {
          if (k > 0) {
            cm->ranges[gapcat+k] = cm->ranges[gapcat];
            cm->conditioned_on[gapcat+k] = cm->conditioned_on[gapcat];
          }
          lst_push_int(cm->conditioned_on[gapcat], cat+k);
        }
        /* NOTE: it's not possible to inherit the "conditioned on"
           categories of the base cat, because it will lead to cycles in
           dependencies.  Problems could result when indel categories
           are conditioned on other categories. */
        if (cm->conditioned_on[cat] != NULL)
          fprintf(stderr, "WARNING: modeling indels for a category that is \"conditioned on\" other categories; check the transition probs of your HMM carefully.\n");
        new_dependencies = 1;
      }
      else 
        for (k = 1; k < range_size; k++) 
          cm->ranges[gapcat+k] = cm->ranges[gapcat];

      /* fill in mappings */
      for (k = 0; k < range_size; k++) {
        gpm->gapcat_to_cat[gapcat] = cat + k;
        gpm->gapcat_to_pattern[gapcat] = j;
        gpm->cat_x_pattern_to_gapcat[cat+k][j] = gapcat;
        gapcat++;
      }
    }
  }

  /* rebuild unspooler, if necessary */
  if (cm->unspooler != NULL || new_dependencies) {
    if (cm->unspooler != NULL) cm_free_unspooler(cm->unspooler);
    cm->unspooler = cm_create_unspooler(cm->ncats + 1, cm->conditioned_on);
  }

  lst_free(indel_cat_nos);
  return gpm;
}

void gp_free_map(GapPatternMap *gpm) {
  int i;
  free(gpm->gapcat_to_cat);
  free(gpm->gapcat_to_pattern);
  for (i = 0; i < gpm->ncats; i++) 
    if (gpm->cat_x_pattern_to_gapcat[i] != NULL)
      free(gpm->cat_x_pattern_to_gapcat[i]);
  free(gpm->cat_x_pattern_to_gapcat);
  free(gpm);
}

/* given a multiple alignment and a tree, fill an array with indices
   describing a "phylogenetic gap pattern" at each position.  The
   value for each column will be an integer i such that i == 0 if
   there are no gaps, 1 <= i <= 2N-3 (where N is the number of leaves
   in the tree) if gaps occur that can be explained by a deletion on
   branch i, 2N-2 <= i <= 4N-6 if gaps occur that can be explained by
   an insertion on branch i, and i == 4N-5 if gaps occur that require
   more than one indel event to explain (so-called "complex" gap
   patterns).  The branches of the tree are numbered from 1 to 2N-3 in
   an in-order traversal.  The branch between the root and its right
   child is ignored, because insertions on one branch below the root
   cannot be distinguished from deletions on the other.  The array
   'patterns' is expected to be allocated to at least size
   msa->length.  It is assumed that no column of the alignment
   consists only of gap characters */
void gp_set_phylo_patterns(int *patterns, MSA *msa, TreeNode *tree) {

  int i, idx, tup;
  List *traversal;
  int *gap_code, *node_to_branch, *leaf_to_seq, *tuple_patterns;
  TreeNode *n;
  String *namestr = str_new(STR_SHORT_LEN);
  int nbranches = tree->nnodes - 2; /* unrooted */
  int complex = nbranches*2 + 1;

  /* require ordered sufficient statistics representation */
  if (msa->ss == NULL)
    ss_from_msas(msa, 1, 1, NULL, NULL, NULL, -1);
  assert(msa->ss->tuple_idx != NULL);

  /* set up mappings of node ids to branch numbers (nodes correspond
     to branches above them), and node ids to sequence indices */
  node_to_branch = smalloc(tree->nnodes * sizeof(int));
  leaf_to_seq = smalloc(tree->nnodes * sizeof(int));
  for (i = 0; i < tree->nnodes; i++) node_to_branch[i] = -1;
  traversal = tr_inorder(tree);
  for (i = 0, idx = 1; i < lst_size(traversal); i++) {
    n = lst_get_ptr(traversal, i);
    if (n != tree && n != tree->rchild) 
      node_to_branch[n->id] = idx++;
    if (n->lchild == NULL) {
      str_cpy_charstr(namestr, n->name);
      if (str_as_int(namestr, &leaf_to_seq[n->id]) != 0 || 
          leaf_to_seq[n->id] < 1 || leaf_to_seq[n->id] > (tree->nnodes+1)/2) {
        fprintf(stderr, "ERROR: cannot map tree leaf '%s' to sequence index.\n", n->name);
        exit(1);
      }
      leaf_to_seq[n->id]--;     /* 0-based indexing */
    }
  }

  gap_code = smalloc(tree->nnodes * sizeof(int));
  /* code is as follows: 
        gap_code[n->id] == 0 -> no gap at node n
        gap_code[n->id] == 1 -> gap at node n 
        gap_code[n->id] == 2 -> ambiguous (implies at least
                                one indel event beneath n) */

  traversal = tr_postorder(tree);
  tuple_patterns = smalloc(msa->ss->ntuples * sizeof(int));
  for (tup = 0; tup < msa->ss->ntuples; tup++) {

    TreeNode *indel_node = NULL; /* node beneath which indel event is
                                    postulated to occur */
    int nchanges = 0;
    tuple_patterns[tup] = 0;
    for (i = 0; i < lst_size(traversal); i++) {
      n = lst_get_ptr(traversal, i);
      if (n->lchild == NULL)
        gap_code[n->id] = 
          (ss_get_char_tuple(msa, tup, leaf_to_seq[n->id], 0) == GAP_CHAR);
      else {                    /* ancestral node */
        if (gap_code[n->lchild->id] == gap_code[n->rchild->id]) {
          if (gap_code[n->lchild->id] == 2) {
                                /* both children have ambiguous codes;
                                   there is no scenario involving only
                                   one indel event */
            tuple_patterns[tup] = complex;
            break;        
          }

          else gap_code[n->id] = gap_code[n->lchild->id];
                                /* children have the same code and it
                                   is unambiguous -- this will be the
                                   most common case */
        }
        else {                  /* children have different codes, one
                                   of which may be ambiguous */

          if (gap_code[n->lchild->id] == 2) {
            gap_code[n->id] = gap_code[n->lchild->id] = gap_code[n->rchild->id];
                                /* lchild ambiguous; the only
                                   allowable scenario has the lchild
                                   and parent equal to the rchild,
                                   with an indel event below the
                                   lchild */

            /* if lchild and one of *its* children now have different
               codes, then an indel event must have occurred between
               them.  Note that lchild must have children; otherwise
               it could not have an indel code of 2 */
            if (gap_code[n->lchild->lchild->id] != 
                gap_code[n->lchild->id] ||
                gap_code[n->lchild->rchild->id] != 
                gap_code[n->lchild->id])
              indel_node = n->lchild;
          }
          else if (gap_code[n->rchild->id] == 2) {
                                /* (symmetric case with rchild) */

            gap_code[n->id] = gap_code[n->rchild->id] = gap_code[n->lchild->id];

            if (gap_code[n->rchild->lchild->id] != 
                gap_code[n->rchild->id] ||
                gap_code[n->rchild->rchild->id] != 
                gap_code[n->rchild->id])
              indel_node = n->rchild;
          }
          else {                /* different and unambiguous */
            if (++nchanges > 1) {
              tuple_patterns[tup] = complex;
              break;
            } 
            if (n == tree) {    /* special case: indel on branch
                                   between subtrees of root; assume
                                   the event occured on the branch to
                                   the left subtree */
              indel_node = tree;
              gap_code[n->id] = gap_code[n->rchild->id];
            }
            else
              gap_code[n->id] = 2;
         }
        }
      }        
    }

    if (tuple_patterns[tup] == complex || nchanges == 0) 
                                /* either zero or more than one indel
                                   events; do nothing (tuple_patterns[tup]
                                   already correct) */
      continue;

    else {                      /* single indel event immediately beneath
                                   indel_node; identify the event and
                                   assign the pattern appropriately */

      if (gap_code[indel_node->id] == 0) {
        /* deletion between indel_node and its left child */
        if (gap_code[indel_node->lchild->id] == 1)
          tuple_patterns[tup] = node_to_branch[indel_node->lchild->id];
        /* deletion between indel_node and its right child */
        else {
          assert(gap_code[indel_node->rchild->id] == 1);
          tuple_patterns[tup] = node_to_branch[indel_node->rchild->id];
        }
      }
      else {                    /* gap_code[indel_node->id] == 1 */
        /* insertion between indel_node and its left child */
        if (gap_code[indel_node->lchild->id] == 0)
          tuple_patterns[tup] = node_to_branch[indel_node->lchild->id] +
            nbranches;
        /* insertion between indel_node and its right child */
        else {
          assert(gap_code[indel_node->rchild->id] == 0);
          tuple_patterns[tup] = node_to_branch[indel_node->rchild->id] +
            nbranches;
        }
      }
    }
  }

  /* finally, assign column patterns from tuple patterns */
  for (i = 0; i < msa->length; i++)
    patterns[i] = tuple_patterns[msa->ss->tuple_idx[i]];

  free(node_to_branch);
  free(leaf_to_seq);
  free(gap_code);
  free(tuple_patterns);
  str_free(namestr);
}
