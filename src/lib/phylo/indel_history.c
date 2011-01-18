/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: indel_history.c,v 1.10 2008-12-10 18:09:17 agd27 Exp $ */

#include <trees.h>
#include <msa.h>
#include <sufficient_stats.h>
#include <indel_history.h>
#include <misc.h>

/* note: when there are nested indels a (compact) indel history is
   *not* a most parsimonious description of indel events; however it
   *is* an complete and unambiguous description of the insertion and
   deletion "characters" in all sequences (i.e., whether each position
   in each sequence is a base, has been deleted, or is "padding"
   required for an insertion) */

/* create new indel history */
IndelHistory *ih_new(TreeNode *tree, int ncols) {
  int i, j;
  IndelHistory *ih = smalloc(sizeof(IndelHistory));
  ih->tree = tree;
  ih->ncols = ncols;
  ih->indel_strings = smalloc(tree->nnodes * sizeof(void*));

  for (i = 0; i < tree->nnodes; i++) {
    ih->indel_strings[i] = smalloc(ncols * sizeof(char));
    /* initialize with base characters */
    for (j = 0; j < ncols; j++)
      ih->indel_strings[i][j] = BASE;
  }
  return ih;
}

/* free indel history */
void ih_free(IndelHistory *ih) {
  int i;
  for (i = 0; i < ih->tree->nnodes; i++)
    sfree(ih->indel_strings[i]);
  sfree(ih->indel_strings);
  sfree(ih);
}

/* create new compact indel history based on alignment and tree */
CompactIndelHistory *ih_new_compact(TreeNode *tree, int ncols) {
  int i;
  CompactIndelHistory *cih = smalloc(sizeof(CompactIndelHistory));
  cih->tree = tree;
  cih->ncols = ncols;
  cih->indels = smalloc(tree->nnodes * sizeof(void*));
  for (i = 0; i < tree->nnodes; i++) 
    cih->indels[i] = lst_new_ptr(max(100, ncols / 10));
  return cih;
}

/* free compact indel history */
void ih_free_compact(CompactIndelHistory *cih) {
  int i, j;
  for (i = 0; i < cih->tree->nnodes; i++) {
    for (j = 0; j < lst_size(cih->indels[i]); j++)
      sfree(lst_get_ptr(cih->indels[i], j));
    lst_free(cih->indels[i]);
  }
  sfree(cih);
}

/* create indel history from compact indel history  */
IndelHistory *ih_expand(CompactIndelHistory *cih) {
  int i, j, k, node_idx;
  List *inside = lst_new_ptr(cih->tree->nnodes), 
    *outside = lst_new_ptr(cih->tree->nnodes);
  TreeNode *n;
  IndelHistory *ih = ih_new(cih->tree, cih->ncols);
  for (i = 0; i < cih->tree->nnodes; i++) {
    tr_partition_nodes(cih->tree, lst_get_ptr(cih->tree->nodes, i), 
                       inside, outside);

    for (j = 0; j < lst_size(cih->indels[i]); j++) {
      Indel *indel = lst_get_ptr(cih->indels[i], j);
      if (indel->type == DEL) {
        for (node_idx = 0; node_idx < lst_size(inside); node_idx++) {
          n = lst_get_ptr(inside, node_idx);
          for (k = 0; k < indel->len; k++)
            ih->indel_strings[n->id][indel->start + k] = DEL;
        }
      }
      else {                    /* indel->type == INS */
        for (node_idx = 0; node_idx < lst_size(outside); node_idx++) {
          n = lst_get_ptr(outside, node_idx);
          for (k = 0; k < indel->len; k++)
            ih->indel_strings[n->id][indel->start + k] = INS;
        }
      }
    }
  }

  lst_free(inside);
  lst_free(outside);

  return ih;
}

/* create compact indel history from indel history */
CompactIndelHistory *ih_compact(IndelHistory *ih) {
  int i, j, k;
  CompactIndelHistory *cih = ih_new_compact(ih->tree, ih->ncols);
  int *ins = smalloc(ih->ncols * sizeof(int));
  List *postorder;
  Indel *indel;
  TreeNode *n;

  /* eliminate all deletions whose parents are deletions -- these are
     implicit */
  postorder = tr_postorder(ih->tree);
  for (i = 0; i < lst_size(postorder); i++) {
    n = lst_get_ptr(postorder, i);
    if (n == ih->tree) continue;
    for (j = 0; j < ih->ncols; j++) 
      if (ih->indel_strings[n->id][j] == DEL && 
          ih->indel_strings[n->parent->id][j] == DEL)
        ih->indel_strings[n->id][j] = BASE + 1;
  }

  /* find the branch of the single insertion event corresponding
     to all insertion gaps.  This will be the branch above the node of
     smallest id that has a base, because ids are assigned in
     preorder */
  for (j = 0; j < ih->ncols; j++) {
    for (i = 0; i < ih->tree->nnodes && ih->indel_strings[i][j] != BASE; )
      i++;
    
    if (i == 0 || i == ih->tree->nnodes)
      ins[j] = -1;
    else 
      ins[j] = i; 
  }

  /* summarize remaining deletions with Indel objects */
  for (i = 0; i < ih->tree->nnodes; i++) {
    for (j = 0; j < ih->ncols; ) {
      if (ih->indel_strings[i][j] == DEL) {
        for (k = 0; 
             j + k < ih->ncols && 
               ih->indel_strings[i][j+k] == DEL;
             k++);
        indel = smalloc(sizeof(Indel));
        indel->type = DEL;
        indel->start = j;
        indel->len = k;
        lst_push_ptr(cih->indels[i], indel);
        j += k;
      }
      else j++;
    }
  }

  /* summarize insertions with Indel objects */
  for (j = 0; j < ih->ncols; ) {
    if (ins[j] > 0) {
      for (k = 0; 
           j + k < ih->ncols && ins[j+k] == ins[j];
           k++);
      indel = smalloc(sizeof(Indel));
      indel->type = INS;
      indel->start = j;
      indel->len = k;
      lst_push_ptr(cih->indels[ins[j]], indel);
      j += k;
    }
    else j++;
  }

  /* restore deletions */
  for (i = 1; i < ih->tree->nnodes; i++) {
    for (j = 0; j < ih->ncols; j++) 
      if (ih->indel_strings[i][j] == BASE + 1)
        ih->indel_strings[i][j] = DEL;
  }

  sfree(ins);
  return cih;
}

void ih_print(IndelHistory *ih, FILE *outf, char *msa_name, char *prog_name) {
  CompactIndelHistory *cih = ih_compact(ih);
  ih_print_compact(cih, outf, msa_name, prog_name);
  ih_free_compact(cih);
}

/* print indel history to file; make sure ancestral nodes have been
   labeled before calling */
void ih_print_compact(CompactIndelHistory *cih, FILE *outf, char *msa_name, char *prog_name) {
  int i, j;
  TreeNode *n;
  fprintf(outf, "## indel history for %s generated by %s\n", 
          msa_name, prog_name);
  fprintf(outf, "## tree: ");
  tr_print(outf, cih->tree, FALSE);
  fprintf(outf, "## ncols: %d\n\n", cih->ncols);
  for (i = 0; i < cih->tree->nnodes; i++) {
    n = lst_get_ptr(cih->tree->nodes, i);

    if (lst_size(cih->indels[i]) == 0) continue;

    fprintf(outf, "s %s\n", n->name);
    for (j = 0; j < lst_size(cih->indels[i]); j++) {
      Indel *indel = lst_get_ptr(cih->indels[i], j);
      fprintf(outf, "%c %d %d\n", indel->type == INS ? 'i' : 'd',
              indel->start, indel->len);
    }
  }    
}

/* convert to an alignment, including sequences for ancestral nodes as
   well as leaf nodes, and with '^' characters in place of '-' for
   insertions and '.' characters in place of '-' for deletions.
   Useful for debugging */
MSA *ih_as_alignment(IndelHistory *ih, MSA *msa) {
  int i, j, k, s=-1, ins=-1;
  char **seqs = smalloc(ih->tree->nnodes * sizeof(char*));
  char **names = smalloc(ih->tree->nnodes * sizeof(char*));
  List *inside, *outside;
  TreeNode *n, *n2;

  inside = lst_new_ptr(10);
  outside = lst_new_ptr(10);

  for (i = 0; i < ih->tree->nnodes; i++) {
    n = lst_get_ptr(ih->tree->nodes, i);
    names[i] = copy_charstr(n->name);
    seqs[i] = smalloc((ih->ncols+1) * sizeof(char));
  }

  for (i = 0; i < ih->tree->nnodes; i++) {
    n = lst_get_ptr(ih->tree->nodes, i);
    /* initialize with actual bases if available or 'N's otherwise */
    if (n->lchild == NULL) {    /* leaf */
      if (msa != NULL) {
        if ((s = msa_get_seq_idx(msa, n->name)) < 0)
          die("ERROR: no match for leaf \"%s\" in alignment.\n", n->name);
      }
      for (j = 0; j < ih->ncols; j++) {
        if (ih->indel_strings[i][j] == BASE) 
          seqs[i][j] = msa == NULL ? 'N' : msa_get_char(msa, s, j);
        else {
	  if (ih->indel_strings[i][j] == INS) { /* Insertion */
	    /* Find the node below the branch where the insertion happened */
	    for (k = 0; 
		 k < ih->tree->nnodes && ih->indel_strings[k][j] != BASE;
		 k++){
	      if (k == 0 || i == ih->tree->nnodes)
		ins = -1;
	      else 
		ins = k;
	    }
	    /* Nodes in the subtree under the node receiving the insertion
	       will all have a non-insertion char while those in the supertree
	       will all have an insertion character */
	    tr_partition_nodes(ih->tree, lst_get_ptr(ih->tree->nodes, ins), 
			       inside, outside);
	    for (k = 0; k < lst_size(inside); k++) {
	      n2 = lst_get_ptr(inside, k);
	      s = msa_get_seq_idx(msa, n2->name);
	      seqs[n2->id][j] = (n2->lchild == NULL) ?
		msa_get_char(msa, s, j) : 'N';
	    }
	    for (k = 0; k < lst_size(outside); k++) {
	      n2 = lst_get_ptr(outside, k);
	      seqs[n2->id][j] = '^';
	    }
	  } else { /* Deletion */
	    seqs[i][j] = '.';
	  }
	}
      }
    }


    else {                      /* ancestor */
      for (j = 0; j < ih->ncols; j++) {
        if (ih->indel_strings[i][j] == BASE) 
          seqs[i][j] = 'N';
        else
          seqs[i][j] = ih->indel_strings[i][j] == INS ? '^' : '.';
      }
    }

    seqs[i][ih->ncols] = '\0';
  }
  lst_free(inside);
  lst_free(outside);
  return msa_new(seqs, names, ih->tree->nnodes, ih->ncols, "ACGTN-^.");
}

CompactIndelHistory *ih_read_compact(FILE *inf) {
  TreeNode *node=NULL, *tree = NULL;
  String *line = str_new(STR_MED_LEN);
  CompactIndelHistory *cih = NULL;
  List *l = lst_new_ptr(3);
  int ncols = -1;

  while (str_readline(line, inf) != EOF) {
    str_trim(line);
    if (line->length == 0) continue;

    if (str_starts_with_charstr(line, "## tree:")) {
      if (line->chars[line->length-1] == ';') 
        line->chars[--line->length] = '\0';
     
      tree = tr_new_from_string(&line->chars[8]);
    }

    else if (str_starts_with_charstr(line, "## ncols:")) {
      str_split(line, NULL, l);
      if (str_as_int(lst_get_ptr(l, 2), &ncols) != 0 || ncols <= 0)
        die("ERROR: bad 'ncols' in indel history file.\n");
      cih = ih_new_compact(tree, ncols);
      lst_free_strings(l);
    }

    else if (line->chars[0] == '#') 
      continue;

    else {
      if (cih == NULL || ncols < 0)
        die("ERROR: malformed header in indel history file.\n");

      str_split(line, NULL, l);

      if (line->chars[0] == 's') {
        node = tr_get_node(tree, ((String*)lst_get_ptr(l, 1))->chars);
        if (node == NULL)
          die("ERROR: no match for node \"%s\" in tree.\n", 
              ((String*)lst_get_ptr(l, 1))->chars);
      }

      else {
        Indel *indel = smalloc(sizeof(Indel));
        if ((line->chars[0] != 'i' && line->chars[0] != 'd') ||
            (str_as_int(lst_get_ptr(l, 1), &indel->start) != 0) ||
            (str_as_int(lst_get_ptr(l, 2), &indel->len) != 0) || 
            node == NULL)
          die("ERROR: bad indel line in history file ('%s')\n",
              line->chars);
        indel->type = line->chars[0] == 'i' ? INS : DEL;

        lst_push_ptr(cih->indels[node->id], indel);
      }

      lst_free_strings(l);
    }
  }

  str_free(line);
  lst_free(l);

  return cih;
}

IndelHistory *ih_new_from_file(FILE* inf) { 
  CompactIndelHistory *cih = ih_read_compact(inf);
  IndelHistory *ih = ih_expand(cih);
  ih_free_compact(cih);
  return ih;
} 

/* extract an indel history from an augmented alignment, including
   sequences for ancestral nodes as well as leaves */
IndelHistory *ih_extract_from_alignment(MSA *msa, TreeNode *tree) {
  int i, j;
  TreeNode *n;
  List *preorder;
  IndelHistory *ih = ih_new(tree, msa->length);
  int *done = smalloc(tree->nnodes * sizeof(int));

  /* first record all gaps as insertions */
  for (i = 0; i < tree->nnodes; i++) done[i] = FALSE;
  for (i = 0; i < msa->nseqs; i++) {
    n = tr_get_node(tree, msa->names[i]);

    if (n == NULL)
      die("ERROR: no match for sequence \"%s\" in tree.\n", msa->names[i]);    

    for (j = 0; j < msa->length; j++) {
      char c = msa_get_char(msa, i, j);
      if (c == GAP_CHAR || c == '^' || c == '.') 
        ih->indel_strings[n->id][j] = INS;
    }

    done[n->id] = TRUE;
  }

  /* make sure all nodes were covered */
  for (i = 0; i < tree->nnodes; i++) 
    if (!done[i]) 
      die("ERROR: no match for node \"%s\" in alignment.\n", 
          ((TreeNode*)lst_get_ptr(tree->nodes, i))->name);

  /* now change gaps that derive from bases to deletions */
  preorder = tr_preorder(tree);
  for (i = 0; i < lst_size(preorder); i++) {
    n = lst_get_ptr(preorder, i);
    if (n == tree) continue;
    for (j = 0; j < msa->length; j++) {
      if (ih->indel_strings[n->id][j] == INS && 
          ih->indel_strings[n->parent->id][j] != INS)
        ih->indel_strings[n->id][j] = DEL;

      /* also check for violation of rule that bases cannot derive
         from deletions */
      else if (ih->indel_strings[n->id][j] == BASE && 
               ih->indel_strings[n->parent->id][j] == DEL)
        die("ERROR: illegal history in column %d; deletions cannot re-emerge as aligned bases.\n", j);
    }
  }

  /* special case: columns of all indels are handled as deletions */
  for (j = 0; j < msa->length; j++) {
    int has_bases = FALSE;
    for (i = 0; !has_bases && i < tree->nnodes; i++) 
      if (ih->indel_strings[i][j] == BASE) has_bases = TRUE;
    if (!has_bases)
      for (i = 0; !has_bases && i < tree->nnodes; i++) 
        ih->indel_strings[i][j] = DEL;
  }

  sfree(done);
  return ih;
}

/* reconstruct an indel history by parsimony from an alignment, given
   a tree */
IndelHistory *ih_reconstruct(MSA *msa, TreeNode *tree) {
  int s, tup, i, j;
  TreeNode *n, *lca;
  char c;
  typedef enum {IGNORE, GAP, OBS_BASE, MISSING, AMBIG} label_type;
  List *postorder;
  IndelHistory *ih = ih_new(tree, msa->length);

  label_type *label = smalloc(tree->nnodes * sizeof(label_type));
  List *inside = lst_new_ptr(tree->nnodes), 
    *outside = lst_new_ptr(tree->nnodes),
    *ambig_cases = lst_new_ptr(tree->nnodes);
  int *seq_to_leaf = smalloc(msa->nseqs * sizeof(int)),
    *leaf_to_seq = smalloc(tree->nnodes * sizeof(int));
  char **tup_hist = smalloc(msa->ss->ntuples * sizeof(void*));

  if (!(msa->ss != NULL && msa->ss->tuple_idx != NULL))
    die("ERROR ih_reconstruct: Need ordered sufficient statistics\n");

  /* build mappings between seqs and leaf indices in tree */
  for (s = 0; s < msa->nseqs; s++) {
    n = tr_get_node(tree, msa->names[s]);
    if (n == NULL)
      die("ERROR: no match for sequence \"%s\" in tree.\n", msa->names[s]);
    seq_to_leaf[s] = n->id;
  }    
  for (i = 0; i < tree->nnodes; i++) {
    n = lst_get_ptr(tree->nodes, i);
    leaf_to_seq[i] = -1;
    if (n->lchild == NULL && n->rchild == NULL) {
      if ((s = msa_get_seq_idx(msa, n->name)) < 0)
        die("ERROR: no match for leaf \"%s\" in alignment.\n", n->name);
      leaf_to_seq[i] = s;
    }
  }

  /* obtain an indel history for each column tuple */
  postorder = tr_postorder(tree);
  for (tup = 0; tup < msa->ss->ntuples; tup++) {
    int min = tree->nnodes, max = -1, ngaps = 0, nmissing = 0,
      skip_root = FALSE;
    checkInterruptN(tup, 1000);

    /* initialize tuple history to all bases */
    tup_hist[tup] = smalloc(tree->nnodes * sizeof(char));
    for (i = 0; i < tree->nnodes; i++) tup_hist[tup][i] = BASE;

    /* find min and max ids of seqs that actually have bases (non-gaps
       and non-missing-data) */
    for (s = 0; s < msa->nseqs; s++) {
      c = ss_get_char_tuple(msa, tup, s, 0);
      if (c == GAP_CHAR) {
        ngaps++;
        continue;
      }
      if (msa->is_missing[(int)c]) {
        nmissing++;
        continue;
      }
      if (seq_to_leaf[s] < min) min = seq_to_leaf[s];
      if (seq_to_leaf[s] > max) max = seq_to_leaf[s];
    }

    /* several special cases allow short cutting */

    if (ngaps == 0) 
      /* impossible to infer gaps in ancestors */
      continue;

    else if (ngaps == 1) {
      /* single base must be deletion, leave others as bases */
      for (i = 0; i < tree->nnodes; i++) 
        if (leaf_to_seq[i] >= 0 &&
            ss_get_char_tuple(msa, tup, leaf_to_seq[i], 0) == GAP_CHAR)
          tup_hist[tup][i] = DEL;
      continue;
    }

    else if (ngaps == msa->nseqs - 1) {
      /* single base must be insertion, so make all others insertion chars */
      for (i = 0; i < tree->nnodes; i++) 
        if (leaf_to_seq[i] == -1 ||
            ss_get_char_tuple(msa, tup, leaf_to_seq[i], 0) == GAP_CHAR)
          tup_hist[tup][i] = INS;
      continue;
    }

    else if (nmissing + ngaps == msa->nseqs) {
      /* all must be deletions */
      for (i = 0; i < tree->nnodes; i++) tup_hist[tup][i] = DEL;
      continue;
    }

    if (!(min >= 0 && max >= min))
      die("ERROR ih_reconstruct min=%e should be >=0 and <= max=%e\n", 
	  min, max);

    /* the LCA of all leaves with bases must be the first ancestor of
       the node with the max id that has an id smaller than the min
       id.  This is based on the assumption that node ids are assigned
       sequentially in a preorder traversal of the tree, which will be
       true as long as the tree is read from a Newick file by the code
       in trees.c */
    for (lca = lst_get_ptr(tree->nodes, max); lca->id > min; 
         lca = lca->parent);

    /* by parsimony, the base was inserted on the branch to the LCA,
       and all ancestral nodes outside the subtree rooted at the LCA
       did not have bases */

    if (lca == tree->lchild || lca == tree->rchild)
      skip_root = TRUE;        /* don't mark root as indel in this case:
                                  can't distinguish insertion from
                                  deletion so assume deletion */

    /* mark ancestral bases outside subtree beneath LCA as insertions
       (or as deletions if skip_root) */
    tr_partition_nodes(tree, lca, inside, outside);
    for (i = 0; i < tree->nnodes; i++) label[i] = OBS_BASE;
    for (i = 0; i < lst_size(outside); i++) {
      n = lst_get_ptr(outside, i);
      label[n->id] = IGNORE;
      if (n == tree && skip_root) 
        continue;               /* skip root if condition above */
      tup_hist[tup][n->id] = skip_root ? DEL : INS;
    }

    /* check for gaps in subtree; if there's at most one, we can take
       a shortcut; otherwise have to use parsimony to infer history in
       subtree */
    ngaps = 0;
    for (i = 0; i < lst_size(inside); i++) {
      n = lst_get_ptr(inside, i);
      if (n->lchild == NULL &&
          ss_get_char_tuple(msa, tup, leaf_to_seq[n->id], 0) == GAP_CHAR)
        ngaps++;
    }
    if (ngaps == 0) 
      continue;
    else if (ngaps == 1) {
      for (i = 0; i < lst_size(inside); i++) {
        n = lst_get_ptr(inside, i);
        if (leaf_to_seq[n->id] >= 0 &&
            ss_get_char_tuple(msa, tup, leaf_to_seq[n->id], 0) == GAP_CHAR)
          tup_hist[tup][n->id] = DEL;
      }
      continue;
    }

    /* use Dollo parsimony to infer the indel history of the subtree
       beneath the LCA.  Use the fact that every base must have a
       chain of bases to the LCA, because, assuming the alignment is
       correct, no insertions are possible beneath the LCA */
    lst_clear(ambig_cases);
    for (i = 0; i < lst_size(postorder); i++) {
      n = lst_get_ptr(postorder, i);
      if (label[n->id] == IGNORE) continue; /* outside subtree */

      /* MISSING means all leaves beneath node have missing data */
      /* AMBIG means combination of gaps and missing data beneath node */

      else if (n->lchild == NULL) {  /* leaf in subtree */
        c = ss_get_char_tuple(msa, tup, leaf_to_seq[n->id], 0);
        if (c == GAP_CHAR)
          label[n->id] = GAP;
        else if (msa->is_missing[(int)c]) 
          label[n->id] = MISSING;
        else
          label[n->id] = OBS_BASE;
      }
      else {                    /* internal node in subtree */
        if (label[n->lchild->id] == OBS_BASE || label[n->rchild->id] == OBS_BASE)
          label[n->id] = OBS_BASE;  /* by Dollo parsimony */
        else if ((label[n->lchild->id] == GAP || label[n->lchild->id] == AMBIG) &&
                 (label[n->rchild->id] == GAP || label[n->rchild->id] == AMBIG))
          label[n->id] = GAP;   /* gaps from both sides and no bases -- must be gap */
        else if (label[n->lchild->id] == MISSING && label[n->rchild->id] == MISSING)
          label[n->id] = MISSING;
        else {              /* must be GAP/MISSING or AMBIG/MISSING */
          label[n->id] = AMBIG;
          lst_push_ptr(ambig_cases, n);
        }
      }
    }

    /* now resolve any ambiguities, by giving each ambiguous node the same
       label as its parent; traversing ambig_cases in reverse order
       ensures that parents are visited before children  */
    if (label[lca->id] != OBS_BASE)
      die("ERROR ih_reconstruct label[%i] (%i) != OBS_BASE (%i)\n",
	  lca->id, label[lca->id], OBS_BASE);
    for (i = lst_size(ambig_cases) - 1; i >= 0; i--) {
      n = lst_get_ptr(ambig_cases, i);
      if (n == lca) continue;
      else label[n->id] = label[n->parent->id];
    }

    /* now mark all gaps inside subtree as deletions */
    for (i = 0; i < lst_size(inside); i++) {
      n = lst_get_ptr(inside, i);
      if (label[n->id] == GAP) 
        tup_hist[tup][n->id] = DEL;
    }
  }

  /* finally, fill out indel history using tuple histories */
  for (i = 0; i < tree->nnodes; i++) {
    for (j = 0; j < msa->length; j++) {
      if (tup_hist[msa->ss->tuple_idx[j]][i] != BASE) {
	if (leaf_to_seq[i] >= 0) 
        {
	  c = ss_get_char_tuple(msa, msa->ss->tuple_idx[j], leaf_to_seq[i], 0);
	  if (!(c==GAP_CHAR || msa->is_missing[(int)c])) die("ERROR reconstructing history in indel_history.c \n");
        }
        ih->indel_strings[i][j] = tup_hist[msa->ss->tuple_idx[j]][i];
      }
    }
  }

  for (tup = 0; tup < msa->ss->ntuples; tup++)
    sfree(tup_hist[tup]);
  sfree(tup_hist);
  lst_free(inside);
  lst_free(outside);
  lst_free(ambig_cases);
  sfree(seq_to_leaf);
  sfree(label);

  return ih;
}

/* convert names in an alignment from the convention used by
   inferAncestors to the convention used in PHAST, based on a given
   tree */
void ih_convert_ia_names(MSA *msa, TreeNode *tree) {
  int i;
  char *newname;
  String **ia_names = smalloc(tree->nnodes * sizeof(void*));
  Hashtable *name_map = hsh_new(tree->nnodes);
  List *postorder = tr_postorder(tree);

  /* create a mapping from Mathieu's names to ours */
  for (i = 0; i < lst_size(postorder); i++) {
    TreeNode *n = lst_get_ptr(postorder, i);
    if (n->lchild == NULL) {
      ia_names[n->id] = str_new_charstr(n->name);
      str_toupper(ia_names[n->id]);
      str_append_char(ia_names[n->id], '+');
    }
    else {
      ia_names[n->id] = str_dup(ia_names[n->lchild->id]);
      str_append(ia_names[n->id], ia_names[n->rchild->id]);
    }
    hsh_put(name_map, ia_names[n->id]->chars, n->name);
  }

  /* now rename */
  for (i = 0; i < msa->nseqs; i++) {
    if ((newname = hsh_get(name_map, msa->names[i])) != (char*)-1) {
      sfree(msa->names[i]);
      msa->names[i] = copy_charstr(newname);
    }
    else 
      die("ERROR: can't convert name '%s'\n", msa->names[i]);
  }

  for (i = 0; i < tree->nnodes; i++)
    str_free(ia_names[i]);
  sfree(ia_names);
  hsh_free(name_map);
}
