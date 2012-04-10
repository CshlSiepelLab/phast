/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <misc.h>
#include <hmm.h>
#include <category_map.h>
#include <gap_patterns.h>

void usage(char *prog) {
  printf("\n\
PROGRAM:       %s\n\
\n\
DESCRIPTION:   Alter transition probabilities in an HMM definition file.\n\
               After specified operations are performed, transition\n\
               probabilities are renormalized and the adjusted file is\n\
               written to standard out.  This program may be used\n\
               multiple times in a pipe.\n\
\n\
USAGE:         %s [OPTIONS] <file.hmm> <cmap.cm>\n\
\n\
OPTIONS:\n\
    -f <cats>  Operate on transitions *from* states corresponding to the \n\
               specified category names (default all)\n\
    -t <cats>  Operate on transitions *to* states corresponding to the \n\
               specified category names (default all)\n\
    -m <fact>  Multiply transition probabilities by the specified factor.\n\
    -a <const> Add the specified constant to transition probabilities.\n\
    -e <val>   Set transition probabilities equal to the specified value.\n\
    -i <icats> Assume a phylo-HMM indel model for states corresponding to \n\
               the specified category names.\n\
    -u <tree>  (Required with -i) Assume given tree topology (.nh file).\n\
    -F <gps>   (For use with -i) Operate on transitions from states corresp.\n\
               to the specified gap-pattern numbers (ANDed with -f).\n\
    -T <gps>   (For use with -i) Operate on transitions to states corresp.\n\
               to the specified gap-pattern numbers (ANDed with -t).\n\
    -z         Equalize transition probabilities.  Set all transition\n\
               probabilities indicated by -f/-t/-F/-T to their overall\n\
               average value.  Options -m and/or -a can be used to adjust\n\
               this average value.\n\
    -R         Restrict to successive transitions within a category range.\n\
    -y         Like -z, except compute separate averages for five classes\n\
               of transitions, based on the gap patterns of the states\n\
               involved: between null gap patterns, between equal\n\
               non-null gap patterns, from null to non-null gap\n\
               patterns, from non-null to null gap patterns, and all\n\
               others.  Useful with the indel model when training data\n\
               is sparse (e.g., for splice-site states).  Options -m and -a\n\
               will be applied to transitions of the 3rd and 5th classes\n\
               described.\n\
    -h         Print this help message.\n\n", prog, prog);
  exit(0);
}

typedef enum {BOTH_NULL, FROM_NULL, TO_NULL, EQUAL, OTHER} trans_type;

trans_type gp_transition(int fgap_pat, int tgap_pat) {
  if (fgap_pat == 0) {
    if (fgap_pat == tgap_pat) return BOTH_NULL;
    else return FROM_NULL;
  }
  else if (tgap_pat == 0) return TO_NULL;
  else {
    if (fgap_pat == tgap_pat) return EQUAL;
    else return OTHER;
  }
}

/* NOTE: might be useful to have options to operate on begin or end
   transitions or to select transitions by gap pattern */

int main(int argc, char *argv[]) {
  double mult_fact = 1, add_const = 0, new_val = -1;
  int i, j, within_range = 0, equalize = 0, equalize_gp = 0;
  List *from_cats = NULL, *to_cats = NULL, *indel_cats = NULL, 
    *from_patterns = NULL, *to_patterns = NULL;
  TreeNode *tree = NULL;
  HMM *hmm;
  CategoryMap *cm;
  int *tcat, *fcat, *tpat=NULL, *fpat=NULL, *isolated;
  GapPatternMap *gpm = NULL;
  char c;
  double sum = 0;
  int count = 0, npasses = 1, pass;
  double gp_sum[5] = {0, 0, 0, 0, 0};
  int gp_count[5] = {0, 0, 0, 0, 0};

  while ((c = getopt(argc, argv, "m:a:e:f:t:i:u:F:T:zyRh")) != -1) {
    switch (c) {
    case 'm':
      mult_fact = get_arg_dbl(optarg);
      break;
    case 'a':
      add_const = get_arg_dbl(optarg);
      break;
    case 'e':
      new_val = get_arg_dbl_bounds(optarg, 0, 1);
      break;
    case 'f':
      from_cats = get_arg_list(optarg);
      break;
    case 't':
      to_cats = get_arg_list(optarg);
      break;
    case 'i':
      indel_cats = get_arg_list(optarg);
      break;
    case 'u':
      tree = tr_new_from_file(phast_fopen(optarg, "r"));
      break;
    case 'F':
      from_patterns = str_list_as_int(get_arg_list(optarg));
      break;
    case 'T':
      to_patterns = str_list_as_int(get_arg_list(optarg));
      break;
    case 'R':
      within_range = 1;
      break;
    case 'z':
      equalize = 1;
      break;
    case 'y':
      equalize_gp = 1;
      break;
    case 'h':
      usage(argv[0]);
    case '?':
      die("Bad argument.  Try '%s -h'.\n", argv[0]);
    }
  }

  if (equalize && equalize_gp) die("Only one of -y and -z may be specified.\n");

  if (optind != argc - 2) 
    die("Input filename required.  Try '%s -h'.\n", argv[0]);

  set_seed(-1);
    
  hmm = hmm_new_from_file(phast_fopen(argv[optind], "r"));
  cm = cm_read(phast_fopen(argv[optind+1], "r"));

  fcat = smalloc((cm->ncats+1) * sizeof(int));
  if (from_cats != NULL) {
    List *l = cm_get_category_list(cm, from_cats, 0);
    for (i = 0; i <= cm->ncats; i++) fcat[i] = 0;
    for (i = 0; i < lst_size(l); i++)
      fcat[lst_get_int(l, i)] = 1;
    lst_free(l);
  }
  else 
    for (i = 0; i <= cm->ncats; i++) fcat[i] = 1;

  tcat = smalloc((cm->ncats+1) * sizeof(int));
  if (to_cats != NULL) {
    List *l = cm_get_category_list(cm, to_cats, 0);
    for (i = 0; i <= cm->ncats; i++) tcat[i] = 0;
    for (i = 0; i < lst_size(l); i++)
      tcat[lst_get_int(l, i)] = 1;
    lst_free(l);
  }
  else 
    for (i = 0; i <= cm->ncats; i++) tcat[i] = 1;

  if (indel_cats != NULL) {
    if (tree == NULL) die("ERROR: must use -u with -i.\n");
    gpm = gp_create_gapcats(cm, indel_cats, tree, FALSE);  

    fpat = smalloc(gpm->ngap_patterns * sizeof(int));
    if (from_patterns != NULL) {
      for (i = 0; i < gpm->ngap_patterns; i++) fpat[i] = 0;
      for (i = 0; i < lst_size(from_patterns); i++) {
        int pat = lst_get_int(from_patterns, i);
        if (pat < 0 || pat > gpm->ngap_patterns) 
          die("ERROR: gap pattern must be in range [0, %d].\n", 
              gpm->ngap_patterns);
        fpat[pat] = 1;
      }
    }
    else
      for (i = 0; i < gpm->ngap_patterns; i++) fpat[i] = 1;

    tpat = smalloc(gpm->ngap_patterns * sizeof(int));
    if (to_patterns != NULL) {
      for (i = 0; i < gpm->ngap_patterns; i++) tpat[i] = 0;
      for (i = 0; i < lst_size(to_patterns); i++) {
        int pat = lst_get_int(to_patterns, i);
        if (pat < 0 || pat > gpm->ngap_patterns) 
          die("ERROR: gap pattern must be in range [0, %d].\n", 
              gpm->ngap_patterns);
        tpat[pat] = 1;
      }
    }
    else
      for (i = 0; i < gpm->ngap_patterns; i++) tpat[i] = 1;
    
  }
  else if (from_patterns != NULL || to_patterns != NULL)
    die("ERROR: -T and -F not valid without -i.\n");

  if (hmm->nstates != (cm->unspooler == NULL ? cm->ncats + 1 : 
                       cm->unspooler->nstates_unspooled)) {
    die("ERROR: number of states in HMM must equal number of site categories (unspooled).\n");
  }

  /* identify isolated states; need to eliminate their default
     self-transitions if they become linked to other states */
  isolated = smalloc(hmm->nstates * sizeof(int));
  for (i = 0; i < hmm->nstates; i++) isolated[i] = 1;
  for (j = 0; j < hmm->nstates; j++)
    for (i = 0; isolated[j] && i < hmm->nstates; i++)
      if (i != j && mm_get(hmm->transition_matrix, i, j) > 0)
        isolated[j] = 0;

  /* need to make two passes if equalizing */
  if (equalize || equalize_gp) npasses = 2;

  for (pass = 1; pass <= npasses; pass++) {
    
    if (pass == 2 && equalize) 
      new_val = sum/count * mult_fact + add_const;

    for (i = 0; i < hmm->nstates; i++) {
      int fbasecat, fgap_pat;
      int fgapped_cat = cm_unspooled_to_spooled_cat(cm, i);
      if (gpm != NULL) {
        fbasecat = gpm->gapcat_to_cat[fgapped_cat];
        fgap_pat = gpm->gapcat_to_pattern[fgapped_cat];
      }
      else {
        fbasecat = fgapped_cat;
        fgap_pat = 0;
      }
    
      if (!fcat[fbasecat] || !fpat[fgap_pat]) continue;

      for (j = 0; j < hmm->nstates; j++) {
        int tbasecat, tgap_pat;
        int tgapped_cat = cm_unspooled_to_spooled_cat(cm, j);
        if (gpm != NULL) {
          tbasecat = gpm->gapcat_to_cat[tgapped_cat];
          tgap_pat = gpm->gapcat_to_pattern[tgapped_cat];
        }
        else {
          tbasecat = tgapped_cat;
          tgap_pat = 0;
        }
        if (!tcat[tbasecat] || !tpat[tgap_pat]) continue;

        if (within_range && (tbasecat != fbasecat + 1 || 
                             cm->ranges[fbasecat]->start_cat_no != 
                             cm->ranges[tbasecat]->start_cat_no)) 
          continue;               /* FIXME: for cyclic ranges, should
                                     allow wraparound (need notion of
                                     cyclic vs. noncyclic) */

        if (equalize && pass == 1)
          sum += mm_get(hmm->transition_matrix, i, j);
        else if (equalize_gp && pass == 1) {
          trans_type type = gp_transition(fgap_pat, tgap_pat);
          gp_sum[type] += mm_get(hmm->transition_matrix, i, j);
          gp_count[type]++;
        }
        else {                  /* either pass == 2 and (equalize or
                                   equalize_gp), OR pass == 1 and
                                   !(equalize or equalize_gp) */

          if (equalize_gp) {
            trans_type type = gp_transition(fgap_pat, tgap_pat);
            new_val = gp_sum[type]/gp_count[type];
            if (type == FROM_NULL || type == OTHER)
              new_val = new_val * mult_fact + add_const;
          }

          if (new_val != -1) 
            mm_set(hmm->transition_matrix, i, j, new_val);
          else
            mm_set(hmm->transition_matrix, i, j, 
                   mm_get(hmm->transition_matrix, i, j) * mult_fact + add_const);

          if (i != j && isolated[i] && 
              mm_get(hmm->transition_matrix, i, j) > 0) {
            mm_set(hmm->transition_matrix, i, i, 0);
            isolated[i] = 0;
          }
                                /* eliminate self-transition for
                                   isolated state if it is no longer
                                   isolated (creates trouble in
                                   normalization) */
        }
      }
    }
  }

  hmm_renormalize(hmm);
  hmm_print(stdout, hmm);

  return 0;
}
