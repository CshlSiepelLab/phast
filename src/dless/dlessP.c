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
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <getopt.h>
#include <misc.h>
#include <bd_phylo_hmm.h>
#include <msa.h>
#include <maf.h>
#include <sufficient_stats.h>
#include <tree_model.h>
#include <subst_distrib.h>
#include "dlessP.help"

/* maximum size of matrix for which to do explicit convolution of
   joint prior; beyond this size an approximation is used.
   Computational complexity is proportional to square of this
   number. */
#define MAX_CONVOLVE_SIZE 22500

void do_p_values(BDPhyloHmm *bdphmm, GFF_Set *predictions,
                 MSA *msa, msa_coord_map *map, char *htmldir, 
                 FILE *timing_f);
void write_html(char *htmldir, GFF_Feature *feat, void *stats,
                event_t type, char *subtree_root_name, String *id);
void write_stats(FILE *F, GFF_Feature *feat, void *stats,
                 event_t type, char *subtree_root_name, String *id);

int main(int argc, char *argv[]) {
  char c;
  int i, opt_idx, old_nnodes;
  FILE *msa_f;
  TreeModel *mod;
  GFF_Set *predictions;
  MSA *msa;
  List *pruned_names = lst_new_ptr(5);
  BDPhyloHmm *bdphmm;
  struct stat st;
  msa_coord_map *map = NULL;

  /* arguments and defaults for options */
  FILE *refseq_f = NULL, *timing_f = NULL;
  msa_format_type msa_format = UNKNOWN_FORMAT;
  int refidx = 1;
  char *htmldir = NULL;

  struct option long_opts[] = {
    {"refseq", 1, 0, 'M'},
    {"msa-format", 1, 0, 'i'},
    {"refidx", 1, 0, 'r'},
    {"timing", 1, 0, 't'},
    {"html", 1, 0, 'H'},
    {"help", 0, 0, 'h'},
    {0, 0, 0, 0}
  };

  while ((c = getopt_long(argc, argv, "r:M:i:t:h", long_opts, &opt_idx)) != -1) {
    switch (c) {
    case 'r':
      refidx = get_arg_int_bounds(optarg, 0, INFTY);
      break;
    case 'M':
      refseq_f = phast_fopen(optarg, "r");
      break;
    case 'i':
      msa_format = msa_str_to_format(optarg);
      if (msa_format == UNKNOWN_FORMAT)
        die("ERROR: unrecognized alignment format.\n");
      break;
    case 't':
      timing_f = phast_fopen(optarg, "w+");;
      break;
    case 'H':
      htmldir = optarg;
      break;
    case 'h':
      printf("%s", HELP);
      exit(0);
    case '?':
      die("Bad argument.  Try 'dlessP -h'.\n");
    }
  }

  if (optind != argc - 3) 
    die("Three arguments required.  Try 'dlessP -h'.\n");

  set_seed(-1);

  /* read input files; do alignment last, because it may be slow */

  /* read tree model */
  mod = tm_new_from_file(phast_fopen(argv[optind+1], "r"), 1);

  /* read predictions */
  predictions = gff_read_set(phast_fopen(argv[optind+2], "r"));

  /* read alignment */
  msa_f = phast_fopen(argv[optind], "r");
  if (msa_format == UNKNOWN_FORMAT)
    msa_format = msa_format_for_content(msa_f, 1);
  fprintf(stderr, "Reading alignment from %s...\n", argv[optind]);
  if (msa_format == MAF) {
    msa = maf_read(msa_f, refseq_f, 1, NULL, NULL, NULL, -1, TRUE, NULL, 
                   NO_STRIP, FALSE); 
  }
  else 
    msa = msa_new_from_file_define_format(msa_f, msa_format, NULL);

  if (msa_alph_has_lowercase(msa)) msa_toupper(msa); 
  msa_remove_N_from_alph(msa);
      
  if (msa->ss == NULL) {
    fprintf(stderr, "Extracting sufficient statistics...\n");
    ss_from_msas(msa, 1, TRUE, NULL, NULL, NULL, -1, 0);
  }
  else if (msa->ss->tuple_idx == NULL)
    die("ERROR: ordered representation of alignment required unless --suff-stats.\n");

  /* prune tree, if necessary */
  old_nnodes = mod->tree->nnodes;
  tm_prune(mod, msa, pruned_names);

  if (lst_size(pruned_names) == (old_nnodes + 1) / 2)
    die("ERROR: no match for leaves of tree in alignment (leaf names must match alignment names).\n");
  if (lst_size(pruned_names) > 0) {
    fprintf(stderr, "WARNING: pruned away leaves of tree with no match in alignment (");
    for (i = 0; i < lst_size(pruned_names); i++)
      fprintf(stderr, "%s%s", ((String*)lst_get_ptr(pruned_names, i))->chars, 
              i < lst_size(pruned_names) - 1 ? ", " : ").\n");
  }

  /* also name ancestors */
  tr_name_ancestors(mod->tree);

  /* set up html dir, if necessary */
  if (htmldir != NULL) {
    fprintf(stderr, "Setting up directory '%s'...\n", htmldir);
    if (stat(htmldir, &st) != 0) {
      if (errno == ENOENT) {	/* missing; create dir */
	#if defined(__MINGW32__)
	    if (mkdir(htmldir) != 0)
    #else
        if (mkdir(htmldir, S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH) != 0)
	#endif
          die("ERROR: cannot create directory '%s'\n", htmldir);
      }
      else 			/* stat returned some other error */
        die("ERROR: bad directory ('%s')\n", htmldir);
    }
    else if (!S_ISDIR(st.st_mode))
      die("ERROR: '%s' exists and is not a directory\n", htmldir);
  }

  bdphmm = bd_new(mod, 0.3, 0.01, 0.01, 0.5, 0.01, 0.01, 0.01,
                  0.01, 0.01, 0.01, FALSE, FALSE, FALSE);
  /* (only the state_to_branch mapping and the category map will
     actually be used) */

  /* need to remap to coord frame of alignment and add offset */
  if (msa->idx_offset != 0) {
    for (i = 0; i < lst_size(predictions->features); i++) {
      GFF_Feature *f = lst_get_ptr(predictions->features, i);
      f->start -= msa->idx_offset;
      f->end -= msa->idx_offset;
    }
  }
  if (refidx != 0) {
    msa_map_gff_coords(msa, predictions, refidx, 0, 0);
    map = msa_build_coord_map(msa, refidx);
  }

  fprintf(stderr, "Computing p-values and generating HTML...\n");
  do_p_values(bdphmm, predictions, msa, map, htmldir, timing_f);

  fprintf(stderr, "Done.\n");

  return 0;
}

void do_p_values(BDPhyloHmm *bdphmm, GFF_Set *predictions, 
                 MSA *msa, msa_coord_map *map, char *htmldir, 
                 FILE *timing_f) {
  int i, j, state;
  int nnodes = bdphmm->phmm->mods[0]->tree->nnodes;
  JumpProcess *jp;
  List *types = lst_new_ptr(nnodes * 2), *type_lists = lst_new_ptr(nnodes * 2);
  TreeModel *mod = bdphmm->phmm->mods[0]; /* nonconserved */
  Regex *id_re = str_re_new(".*id \"([^\"]*)\"");
  String *id = str_new(STR_SHORT_LEN);
  List *l = lst_new_ptr(1);

  /* partition predictions by type */
  gff_partition_by_type(predictions, types, type_lists);

  /* write header of output file */
  write_stats(stdout, NULL, NULL, -1, NULL, NULL);

  for (i = 0; i < lst_size(types); i++) {
    List *feats_this_type = lst_get_ptr(type_lists, i);
    TreeNode *orig_tree, *subtree_root=NULL, *tmp;
    event_t event_type;
    p_value_stats *stats_cons = NULL;
    p_value_joint_stats *stats_bd = NULL;

    fprintf(stderr, "  (category '%s' [%d features])...\n",
            ((String*)lst_get_ptr(types, i))->chars, lst_size(feats_this_type));

    state = cm_get_category(bdphmm->phmm->cm, lst_get_ptr(types, i));

    if (state == nnodes) {    /* fully conserved state */
      jp = sub_define_jump_process(mod, 1e-10, tr_total_len(mod->tree));
      stats_cons = sub_p_value_many(jp, msa, feats_this_type, -1);
      event_type = CONS;
      sub_free_jump_process(jp);
    }
    else {                    /* birth or death state */
      event_type = state < nnodes ? DEATH : BIRTH;

      /* first reroot tree, using copy */
      orig_tree = mod->tree;
      mod->tree = tr_create_copy(orig_tree);

      subtree_root = lst_get_ptr(mod->tree->nodes, 
                                 bdphmm->state_to_branch[state]);

      /* check for special case of a subtree just below the root and a
         supertree consisting of a single leaf node.  This case is
         meaningless if the tree is unrooted and the branch in
         question is considered part of the subtree, as here -- the
         supertree has zero total branch length and a trivial
         distribution of numbers of substitutions.  */ 
      if ((mod->tree->lchild == subtree_root && 
           mod->tree->rchild->lchild == NULL) ||
          (mod->tree->rchild == subtree_root && 
           mod->tree->lchild->lchild == NULL)) {
        fprintf(stderr, "WARNING: ignoring type '%s'; supertree consists of a single leaf node.\n", ((String*)lst_get_ptr(types, i))->chars);
        lst_free(feats_this_type);
        continue;
      }

      tr_reroot(mod->tree, subtree_root, TRUE); /* include branch above node */
      mod->tree = subtree_root->parent;

      /* swap left and right children.  This is necessary because
         routines for computing joint distrib assume branch to right
         has length zero, but because branch is included, tr_reroot
         will put zero length branch on left */
      tmp = mod->tree->lchild;
      mod->tree->lchild = mod->tree->rchild;
      mod->tree->rchild = tmp;

      /* now get stats */
      jp = sub_define_jump_process(mod, 1e-10, tr_total_len(mod->tree));
      stats_bd = sub_p_value_joint_many(jp, msa, feats_this_type, -1, 
                                        MAX_CONVOLVE_SIZE, timing_f);

      tr_free(mod->tree);
      sub_free_jump_process(jp);
      mod->tree = orig_tree;
    }

    for (j = 0; j < lst_size(feats_this_type); j++) {
      GFF_Feature *f = lst_get_ptr(feats_this_type, j);
        
      /* grab id from attribute */
      if (! (str_re_match(f->attribute, id_re, l, 1) >= 0)) 
	die("ERROR in do_pvalues: could not parse %s\n", f->attribute);
      str_cpy(id, lst_get_ptr(l, 1));
      lst_free_strings(l);

      /* map coords back to orig, if necessary */
      if (map != NULL) {
        f->start = msa_map_msa_to_seq(map, f->start);
        f->end = msa_map_msa_to_seq(map, f->end);
      }
      f->start += msa->idx_offset;
      f->end += msa->idx_offset;

      write_stats(stdout, f, 
                  event_type == CONS ? (void*)(&stats_cons[j]) : (void*)(&stats_bd[j]), 
                  event_type, subtree_root->name, id);

      if (htmldir != NULL)
        write_html(htmldir, f, 
                   event_type == CONS ? (void*)(&stats_cons[j]) : (void*)(&stats_bd[j]), 
                   event_type, subtree_root->name, id);
    }

    if (event_type == CONS)
      sfree(stats_cons);
    else
      sfree(stats_bd);

    lst_free(feats_this_type);
  }

  lst_free(types);
  lst_free(type_lists);
  str_re_free(id_re);
  str_free(id);
  lst_free(l);
}

void write_html(char *htmldir, GFF_Feature *feat, void *stats,
                event_t type, char *subtree_root_name, String *id) {

  FILE *F;
  char tmpstr[STR_MED_LEN];
  p_value_stats *s_cons = NULL;
  p_value_joint_stats *s_bd = NULL;

  if (type == CONS)
    s_cons = stats;
  else 
    s_bd = stats;

  sprintf(tmpstr, "%s/%s.html", htmldir, id->chars); 
  F = phast_fopen(tmpstr, "w+");

  fprintf(F, "\
 <!DOCTYPE HTML PUBLIC \"-//IETF//DTD HTML//EN\">\n\
 <html> <head>\n\
 <title>DLESS output for %s</title>\n\
 </head>\n\
 \n\
 <body>\n\
 <h3>DLESS output for %s</h3>\n", id->chars, id->chars);

  if (type == CONS) 
    fprintf(F, "\
 Prediction: conserved in all species<br>\n\
 Log-odds score wrt neutral: %.1f bits<br>\n\
 P-value of conservation: %e<br>\n\
<br>\n\
 Numbers of substitutions:\n\
 <ul>\n\
 <li>Null distribution: mean = %f, var = %f, 95%% c.i. = [%d, %d]\n\
 <li>Posterior distribution: mean = %f, var = %f\n\
 </ul>\n", 
            feat->score, s_cons->p_cons, s_cons->prior_mean, 
            s_cons->prior_var, s_cons->prior_min, s_cons->prior_max,
            s_cons->post_mean, s_cons->post_var);
  else {
    fprintf(F, "\
 Prediction: %s of element on branch above node labeled \"%s\"<br>\n\
 Log-odds score wrt neutral: %.1f bits<br>\n\
 P-value of conservation in subtree: %e<br>\n\
 P-value of conservation in rest of tree: %e<br>\n\
 P-value of conservation in subtree given total: %e%s<br>\n\
 P-value of conservation in rest of tree given total: %e%s<br>\n\
 <br>\n\
 Numbers of substitutions in subtree beneath event:\n\
 <ul>\n\
 <li>Null distribution: mean = %f, var = %f, 95%% c.i. = [%d, %d]\n\
 <li>Posterior distribution: mean = %f, var = %f\n\
 </ul>\n\
 Numbers of substitutions in rest of tree:\n\
 <ul>\n\
 <li>Null distribution: mean = %f, var = %f, 95%% c.i. = [%d, %d]\n\
 <li>Posterior distribution: mean = %f, var = %f\n\
 </ul>\n", type == BIRTH ? "gain" : "loss", 
            subtree_root_name, feat->score, s_bd->p_cons_left,
            s_bd->p_cons_right, 
            s_bd->cond_p_cons_left, s_bd->cond_p_approx ? "*" : "",
            s_bd->cond_p_cons_right, s_bd->cond_p_approx ? "*" : "",
            s_bd->prior_mean_left, s_bd->prior_var_left, s_bd->prior_min_left, 
            s_bd->prior_max_left, s_bd->post_mean_left, s_bd->post_var_left, 
            s_bd->prior_mean_right, s_bd->prior_var_right, s_bd->prior_min_right, 
            s_bd->prior_max_right, s_bd->post_mean_right, s_bd->post_var_right);
    if (s_bd->cond_p_approx)
      fprintf(F, "* = Approximate p-value (usually conservative)<br>\n");
  }

  fprintf(F, "</body> </html>\n");
    
  phast_fclose(F);
}

void write_stats(FILE *F, GFF_Feature *feat, void *stats,
                 event_t type, char *subtree_root_name, String *id) {

  p_value_stats *s_cons = NULL;
  p_value_joint_stats *s_bd = NULL;

  if (type == -1) {             /* code to print header */
    fprintf(F, "#chr\tstart\tend\tid\tscore\ttype\tbranch\tpConsSub\tpConsSup\tpConsSubCond\tpConsSupCond\tcondApprox\tpriorMeanSub\tpriorVarSub\tpriorMinSub\tpriorMaxSub\tpostMeanSub\tpostVarSub\tpriorMeanSup\tpriorVarSup\tpriorMinSup\tpriorMaxSup\tpostMeanSup\tpostVarSup\n");
    return;
  }

  if (type == CONS)
    s_cons = stats;
  else 
    s_bd = stats;

  if (type == CONS) 
    fprintf(F, "%s\t%d\t%d\t%s\t%.1f\t%s\t%s\t%.3e\t%.3e\t%.3e\t%.3e\t%s\t%.2f\t%.2f\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%d\t%.2f\t%.2f\n",
            feat->seqname->chars, feat->start-1, feat->end, id->chars, 
            feat->score, "conserved", "NA", 
            s_cons->p_cons, (double)1, (double)1, (double)1, "exact", 
            s_cons->prior_mean, s_cons->prior_var, 
            s_cons->prior_min, s_cons->prior_max, 
            s_cons->post_mean, s_cons->post_var, 
            (double)0, (double)0, 0, 0, (double)0, (double)0);
  else 
    fprintf(F, "%s\t%d\t%d\t%s\t%.1f\t%s\t%s\t%.3e\t%.3e\t%.3e\t%.3e\t%s\t%.2f\t%.2f\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%d\t%.2f\t%.2f\n", 
            feat->seqname->chars, feat->start-1, feat->end, id->chars, 
            feat->score, type == BIRTH ? "gain" : "loss", subtree_root_name, 
            s_bd->p_cons_left, s_bd->p_cons_right, 
            s_bd->cond_p_cons_left, s_bd->cond_p_cons_right, 
            s_bd->cond_p_approx ? "approx" : "exact",
            s_bd->prior_mean_left, s_bd->prior_var_left, 
            s_bd->prior_min_left, s_bd->prior_max_left, 
            s_bd->post_mean_left, s_bd->post_var_left, 
            s_bd->prior_mean_right, s_bd->prior_var_right, 
            s_bd->prior_min_right, s_bd->prior_max_right, 
            s_bd->post_mean_right, s_bd->post_var_right);
}

