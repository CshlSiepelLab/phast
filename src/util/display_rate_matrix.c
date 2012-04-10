/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: display_rate_matrix.c,v 1.6 2008-11-12 02:07:58 acs Exp $ */

/* simple program to display the contents of a rate matrix in various
   ways (see USAGE) */

#include <stdlib.h>
#include <stdio.h>
#include "tree_model.h"
#include <getopt.h>
#include <stringsplus.h>
#include <ctype.h>

void print_usage() {
  fprintf(stdout, "PROGRAM: display_rate_matrix\n\
USAGE: display_rate_matrix [OPTIONS] <model_fname>\n\
OPTIONS:\n\
    -t <t>: Output P(t) = exp(Qt) instead of Q.  Requires t >= 0.  \n\
            Use \"-t A\" to output a matrix for each branch of the tree.\n\
    -f:     Show equilibrium frequencies as an additional table row.  \n\
            In list node they are shown with first tuple being -.\n\
    -e:     Show \"exchangeabilities\" instead of raw matrix elements \n\
            (that is, divide each element by the equilibrium frequency \n\
            of its column).  Not available with -t.\n\
    -d:     Suppress printing of elements on main diagonal.\n\
    -L:     Format table for typesetting with LATEX.  Incompatible with -l.\n\
    -l:     Show matrix elements as a list rather than as a table.  \n\
            When -t is not specified (rate matrix case), only off-diagonal\n\
            elements will be printed.\n\
    -i:     (For use with -l only) Report whether each substitution is\n\
            a transition or a transversion.\n\
    -z:     (For use with -l) Report elements equal to zero (omitted by \n\
            default, except with -t).  Implied by -a.\n\
    -S:     (For use with -l)  Assume a symmetric matrix and report half \n\
            as many lines.  Useful with -e.\n\
    -E:     (for use with -l) Print rates and probabilities \n\
            in scientific notation (format %%e instead of %%f).\n\
    -a:     (Requires a model of order 3).  Replace a matrix of codon\n\
            substitution rates with the induced matrix of amino acid\n\
            substitution rates, according to the universal genetic\n\
            code.  See Yang, Nielsen, and Hasegawa, 1998.\n\
    -s:     (For use with -a)  Include stop codons (by default suppressed).\n\
    -M <f>: (For use with -l only; implies -a) Read an amino-acid\n\
            substitution matrix from file <f> and report values from\n\
            this matrix with the induced amino acid substitution rates.\n\
            Matrix should be in the format used by BLAST (as\n\
            produced by the NCBI \"pam\" program) \n\
    -N <f>  Like -M but for matrices in the format used by the PAML \n\
            package for amino acid substitution and rate matrices.\n\
    -A <f>: (For use with -l only and not with -M/-N)  Read alternative\n\
            substitution scores from file <f> and report values in\n\
            output.  File <f> should have three columns: a \"from\"\n\
            tuple, a \"to\" tuple, and a real-valued score.\n\
            Substitutions not listed will be given null scores and\n\
            reported as \"NA\".\n\
    -B <f>  Like -A but compares to rates of a single-nucleotide model \n\
            (order 1).  File <f> should be a standard tree model (.mod) file.\n\
    -C      Report context-dependent transition/transversion rates, as \n\
            shown in Tables 2 and 3 of Morton et al., JME 45:227-231, 1997. \n\
            Requires a model of order 3 with a DNA alphabet.\n\
    -h      Print this help message.\n\n");
}

/* given a model, return the tuple of characters corresponding to the
   specified state */
void get_state_tuple(TreeModel *mod, char *tuple, int state) {
  int alph_size = strlen(mod->rate_matrix->states);
  int tuple_idx;
  for (tuple_idx = -1*mod->order; tuple_idx <= 0; tuple_idx++) {
    int projection = (state / int_pow(alph_size, -1 * tuple_idx)) % 
      alph_size;
    tuple[tuple_idx + mod->order] = mod->rate_matrix->states[projection];
  }
  tuple[mod->order+1] = '\0';
}

/* Read substitution scores from specified file and return as a kind
   of pseudo substitution matrix.  All nonspecified elements in matrix
   will be equal to NEGINFTY, which is to be interpretted as "NA" */
Matrix* read_subst_scores(TreeModel *mod, FILE *F) {
  Matrix *retval = mat_new(mod->rate_matrix->size,
                                        mod->rate_matrix->size);
  String *line = str_new(STR_MED_LEN), *tuple1, *tuple2;
  List *l = lst_new_ptr(3);
  int alph_size = strlen(mod->rate_matrix->states);
  int *inv_alph = mod->rate_matrix->inv_states;
  double val;
  mat_set_all(retval, NEGINFTY);
  while (str_readline(line, F) != EOF) {
    str_double_trim(line);
    if (str_starts_with_charstr(line, "#") || line->length == 0) 
      continue;
    str_split(line, NULL, l);
    if (lst_size(l) < 3) {
      die("ERROR: wrong number of columns in subst. score file.\n");
    }
    tuple1 = lst_get_ptr(l, 0);
    tuple2 = lst_get_ptr(l, 1);
    if (str_as_dbl(lst_get_ptr(l, 2), &val) != 0) {
      die("ERROR: bad value in subst. score file.\n");
    }
    mat_set(retval, tuple_index(tuple1->chars, inv_alph, alph_size),
                   tuple_index(tuple2->chars, inv_alph, alph_size), val);
    str_free(tuple1); str_free(tuple2); str_free(lst_get_ptr(l, 2));
  }
  lst_free(l);
  str_free(line);
  return retval;
}

/* Read an amino acid rate matrix in the format used by PAML.  Reorder
   the rows and columns to match 'alph'.  Warning: the ordering in the
   file is assumed to match that used in the files in the PAML
   distribution (alphabetical order of 3-letter codes), which is also
   the order of AA_ALPHABET (therefore AA_ALPHABET may not be
   changed!).  Equilibrium frequencies are ignored.  */ 
Matrix *read_paml_matrix(FILE *F, char *alph) {
  char *paml_alph = "ARNDCQEGHILKMFPSTWYV$";
  int size = strlen(paml_alph);
  Matrix *retval = mat_new(size, size);
  List *fields = lst_new_ptr(100);
  String *line = str_new(STR_MED_LEN);
  int i, j;
  if (strcmp(alph, paml_alph) != 0)
    die("ERROR read_paml_matrix (alph (%s) != paml_alph (%s))\n",
	alph, paml_alph);
  mat_zero(retval);

  for (i = 1; i < size-1 && str_readline(line, F) != EOF; ) {
    /* NOTE: size of matrix allows for stop, but stop not included in
       file; therefore, only read size-1 lines */
    str_double_trim(line);
    if (line->length == 0) continue;
    str_split(line, NULL, fields);
    if (lst_size(fields) != i) {
      die("ERROR: row %d of matrix must have %d columns.\n",
	  i+1, i);
    }
    for (j = 0; j < lst_size(fields); j++) {
      double val;

      if (str_as_dbl(lst_get_ptr(fields, j), &val) != 0) {
        die("ERROR: non-numeric matrix element in subst. matrix ('%s')\n", 
	    ((String*)lst_get_ptr(fields, j+1))->chars);
      }
      str_free(lst_get_ptr(fields, j));

      if (j >= size)
	die("ERROR read_paml_matrix j (%i) should be < size (%i)\n", j, size);
      mat_set(retval, i, j, val);
      mat_set(retval, j, i, val);
    }
    i++;
  }

  if (i != size - 1) {
    die("ERROR: too few rows in subst. matrix.\n");
  }
  
  lst_free(fields);
  str_free(line);
  return retval;
}

Matrix* unproject_rates(TreeModel *mod_tuples, TreeModel *mod_single) {
  int dim = mod_tuples->rate_matrix->size;
  int alph_size = strlen(mod_tuples->rate_matrix->states);
  char tuple_i[mod_tuples->order+1], tuple_j[mod_tuples->order+1];
  int position, i, j;
  Matrix *retval = mat_new(dim, dim);
  mat_zero(retval);
  for (i = 0; i < dim; i++) {
    get_tuple_str(tuple_i, i, mod_tuples->order+1, 
                  mod_tuples->rate_matrix->states);
    for (j = 0; j < dim; j++) {
      if (i == j || mm_get(mod_tuples->rate_matrix, i, j) == 0) continue;
      /* WARNING: we'll ignore any rate matrix elements that have
         *actually been estimated* to be zero */

      get_tuple_str(tuple_j, j, mod_tuples->order+1, 
                    mod_tuples->rate_matrix->states);
      position = mod_tuples->order - floor(log(abs(i - j))/log(alph_size));
      mat_set(retval, i, j, 
                     mm_get(mod_single->rate_matrix, 
                            mod_single->rate_matrix->inv_states[(int)tuple_i[position]],
                            mod_single->rate_matrix->inv_states[(int)tuple_j[position]]));
    }
  }
  return retval;
}

/* this function implements the -C option */
void do_context_dependent_ti_tv(TreeModel *mod) {
  char *alph = mod->rate_matrix->states;
  int alph_size = strlen(alph);
  char tuple_i[mod->order+2], tuple_j[mod->order+2];
  double context_ti[alph_size][alph_size], context_tv[alph_size][alph_size],
    ti_AT_5_pyrim[3][3], tv_AT_5_pyrim[3][3],
    all_ti, all_tv, expected_rate;
  int mid_src, mid_targ, first, last, at_states[2], gc_states[2], i, j, 
    num_5_pyrim;
  if (mod->order != 2)
    die("ERROR do_context_dependent_ti_tv: mod->order (%i) should be 2\n",
	mod->order);
  if (alph_size != 4)
    die("ERROR do_contect_dependent_ti_tv: alph_size (%i) should be 4\n", 
	alph_size);

  tuple_i[mod->order+1] = tuple_j[mod->order+1] = '\0';

  /* We only care about substitutions at the middle position */
  for (first = 0; first < alph_size; first++) {
    for (last = 0; last < alph_size; last++) {
      context_ti[first][last] = context_tv[first][last] = 0;
      for (mid_src = 0; mid_src < alph_size; mid_src++) {
        for (mid_targ = 0; mid_targ < alph_size; mid_targ++) {
          if (mid_src == mid_targ) continue;

          i = last + mid_src*alph_size + first*alph_size*alph_size;
          j = last + mid_targ*alph_size + first*alph_size*alph_size;
          expected_rate = mm_get(mod->rate_matrix, i, j) * 
            vec_get(mod->backgd_freqs, i);

/*           get_tuple_str(tuple_i, i, mod->order+1,  */
/*                         mod->rate_matrix->states); */
/*           get_tuple_str(tuple_j, j, mod->order+1,  */
/*                         mod->rate_matrix->states); */

/*           if ((tuple_i[1] == 'C' && tuple_i[2] == 'G') || tuple_i[0] == 'C' && tuple_i[1] == 'G') continue; */

          if (is_transition(alph[mid_src], alph[mid_targ]))
            context_ti[first][last] += expected_rate;
          else if (i != j)
            context_tv[first][last] += expected_rate;
        }
      }
    }
  }

  /* first print equivalent of table 2 */

  /* get "states" associated with A+T bases and G+C bases */
  at_states[0] = mod->rate_matrix->inv_states['A'];
  at_states[1] = mod->rate_matrix->inv_states['T'];
  gc_states[0] = mod->rate_matrix->inv_states['C'];
  gc_states[1] = mod->rate_matrix->inv_states['G'];

  /* header */
  printf("%5s %5s %5s %8s %8s %8s\n", "A+T", "5'", "3'", "Ts", "Tv", "Tv/(Tv+Ts)");

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      ti_AT_5_pyrim[i][j] = tv_AT_5_pyrim[i][j] = 0;

  /* AT content 0 (GC either side) */
  all_ti = all_tv = 0;
  for (first = 0; first < 2; first++) {
    for (last = 0; last < 2; last++) {
      printf("%5d %5c %5c %8.4f %8.4f %8.4f\n", 0, 
             alph[(int)gc_states[first]], 
             alph[(int)gc_states[last]], 
             context_ti[gc_states[first]][gc_states[last]],
             context_tv[gc_states[first]][gc_states[last]],
             context_tv[gc_states[first]][gc_states[last]] /
             (context_ti[gc_states[first]][gc_states[last]] +
              context_tv[gc_states[first]][gc_states[last]]));
      all_ti += context_ti[gc_states[first]][gc_states[last]];
      all_tv += context_tv[gc_states[first]][gc_states[last]];
      
      num_5_pyrim = IS_PYRIMIDINE(alph[(int)gc_states[first]]);
      num_5_pyrim += IS_PYRIMIDINE(msa_compl_char(alph[(int)gc_states[last]]));
      ti_AT_5_pyrim[0][num_5_pyrim] += 
        context_ti[gc_states[first]][gc_states[last]];
      tv_AT_5_pyrim[0][num_5_pyrim] += 
        context_tv[gc_states[first]][gc_states[last]];
    }
  }
  printf("%5d %5s %5s %8.4f %8.4f %8.4f\n", 0, "all", "all", all_ti, 
         all_tv, all_tv / (all_ti + all_tv));
      
  /* AT content 1 (AT one side, GC other side) */
  all_ti = all_tv = 0;
  for (first = 0; first < 2; first++) {
    for (last = 0; last < 2; last++) {
      printf("%5d %5c %5c %8.4f %8.4f %8.4f\n", 1, 
             alph[(int)at_states[first]], 
             alph[(int)gc_states[last]], 
             context_ti[at_states[first]][gc_states[last]],
             context_tv[at_states[first]][gc_states[last]],
             context_tv[at_states[first]][gc_states[last]] /
             (context_ti[at_states[first]][gc_states[last]] +
              context_tv[at_states[first]][gc_states[last]]));
      printf("%5d %5c %5c %8.4f %8.4f %8.4f\n", 1, 
             alph[(int)gc_states[first]], 
             alph[(int)at_states[last]], 
             context_ti[gc_states[first]][at_states[last]],
             context_tv[gc_states[first]][at_states[last]],
             context_tv[gc_states[first]][at_states[last]] /
             (context_ti[gc_states[first]][at_states[last]] +
              context_tv[gc_states[first]][at_states[last]]));
      all_ti += context_ti[at_states[first]][gc_states[last]] +
        context_ti[gc_states[first]][at_states[last]];
      all_tv += context_tv[at_states[first]][gc_states[last]] +
        context_tv[gc_states[first]][at_states[last]];

      num_5_pyrim = IS_PYRIMIDINE(alph[(int)at_states[first]]);
      num_5_pyrim += IS_PYRIMIDINE(msa_compl_char(alph[(int)gc_states[last]]));
      ti_AT_5_pyrim[1][num_5_pyrim] += 
        context_ti[at_states[first]][gc_states[last]];
      tv_AT_5_pyrim[1][num_5_pyrim] += 
        context_tv[at_states[first]][gc_states[last]];

      num_5_pyrim = IS_PYRIMIDINE(alph[(int)gc_states[first]]);
      num_5_pyrim += IS_PYRIMIDINE(msa_compl_char(alph[(int)at_states[last]]));
      ti_AT_5_pyrim[1][num_5_pyrim] += 
        context_ti[gc_states[first]][at_states[last]];
      tv_AT_5_pyrim[1][num_5_pyrim] += 
        context_tv[gc_states[first]][at_states[last]];
    }
  }
  printf("%5d %5s %5s %8.4f %8.4f %8.4f\n", 1, "all", "all", all_ti, 
         all_tv, all_tv / (all_ti + all_tv));

  /* AT content 2 (AT both sides) */
  all_ti = all_tv = 0;
  for (first = 0; first < 2; first++) {
    for (last = 0; last < 2; last++) {
      printf("%5d %5c %5c %8.4f %8.4f %8.4f\n", 2, 
             alph[(int)at_states[first]], 
             alph[(int)at_states[last]], 
             context_ti[at_states[first]][at_states[last]],
             context_tv[at_states[first]][at_states[last]],
             context_tv[at_states[first]][at_states[last]] /
             (context_ti[at_states[first]][at_states[last]] +
              context_tv[at_states[first]][at_states[last]]));
      all_ti += context_ti[at_states[first]][at_states[last]];
      all_tv += context_tv[at_states[first]][at_states[last]];

      num_5_pyrim = IS_PYRIMIDINE(alph[(int)at_states[first]]);
      num_5_pyrim += IS_PYRIMIDINE(msa_compl_char(alph[(int)at_states[last]]));
      ti_AT_5_pyrim[2][num_5_pyrim] += 
        context_ti[at_states[first]][at_states[last]];
      tv_AT_5_pyrim[2][num_5_pyrim] += 
        context_tv[at_states[first]][at_states[last]];
    }
  }
  printf("%5d %5s %5s %8.4f %8.4f %8.4f\n", 2, "all", "all", all_ti, 
         all_tv, all_tv / (all_ti + all_tv));

  /* now print equivalent of table 3 */
  printf("\n\n%5s %5s %8s %8s %8s\n", "A+T", "5'Y", "Ts", "Tv", "Tv/(Tv+Ts)");
  for (i = 0; i < 3; i++) 
    for (j = 0; j < 3; j++) 
      printf("%5d %5d %8.4f %8.4f %8.4f\n", i, j, ti_AT_5_pyrim[i][j], tv_AT_5_pyrim[i][j], tv_AT_5_pyrim[i][j]/(tv_AT_5_pyrim[i][j] + ti_AT_5_pyrim[i][j]));
  for (j = 0; j < 3; j++) {
    double all_ti = ti_AT_5_pyrim[0][j] + ti_AT_5_pyrim[1][j] + 
      ti_AT_5_pyrim[2][j];
    double all_tv = tv_AT_5_pyrim[0][j] + tv_AT_5_pyrim[1][j] + 
      tv_AT_5_pyrim[2][j];
    printf("%5s %5d %8.4f %8.4f %8.4f\n", "all", j, all_ti, all_tv, all_tv/(all_tv + all_ti));
  }
}

int main(int argc, char* argv[]) {
  FILE* F;
  TreeModel *model;
  int i, j, k, alph_size, nstates, do_eqfreqs = 0, exch_mode = 0, 
    list_mode = 0, latex_mode = 0, suppress_diag = 0, ti_tv = 0, 
    scientific_mode = 0,
    induced_aa = 0, do_stop_codons = 0, do_zeroes = 0, symmetric = 0, 
    context_ti_tv = 0, all_branches = 0;
  int startcol, endcol, ncols, branch_no = 0, matrix_idx = 0;
/*   int aa_inv[256]; */
  double t = -1, total_ti = 0, total_tv = 0, rho_s = 0, cpg_ti = 0, 
    cpg_tv = 0, non_cpg_ti = 0, non_cpg_tv = 0, cpg_eqfreq = 0;
  char *rate_format_string = "%8.6f";
  MarkovMatrix *M;
  char c;
  char tuple[5], tuple2[5]; /* , aa_alph[50]; */
  char *subst_mat_fname = NULL, *subst_score_fname = NULL, 
    *subst_mat_fname_paml = NULL, *order1_mod_fname = NULL;
  Matrix *subst_mat = NULL;
  List *matrix_list = lst_new_ptr(20), *traversal = NULL;

  while ((c = getopt(argc, argv, "t:fedlLiM:N:A:B:aszSECh")) != -1) {
   switch(c) {
    case 't':
      if (optarg[0] == 'A') all_branches = 1;
      else t = get_arg_dbl_bounds(optarg, 0, INFTY);
      break;
    case 'f':
      do_eqfreqs = 1;
      break;
    case 'e':
      exch_mode = 1;
      break;
    case 'd':
      suppress_diag = 1;
      break;
    case 'l':
      list_mode = 1;
      break;
    case 'L':
      latex_mode = 1;
      break;
    case 'i':
      ti_tv = 1;
      break;
    case 'M':
      subst_mat_fname = optarg;
      induced_aa = 1;
      break;
    case 'N':
      subst_mat_fname_paml = optarg;
      induced_aa = 1;
      break;
    case 'A':
      subst_score_fname = optarg;
      break;
    case 'B':
      order1_mod_fname = optarg;
      break;
    case 'a':
      induced_aa = 1;
      do_zeroes = 1;
      break;
    case 's':
      do_stop_codons = 1;
      break;
    case 'z':
      do_zeroes = 1;
      break;
    case 'S':
      symmetric = 1;
      break;
    case 'E':
      scientific_mode = 1;
      rate_format_string = "%13.6e";
      break;
    case 'C':
      context_ti_tv = 1;
      break;
    case 'h':
      print_usage();
      exit(0);
    case '?':
      die("Unrecognized option.  Try \"display_rate_matrix -h\" for help.\n");
    }
  }

  set_seed(-1);

  if ((t >= 0 && exch_mode) || (latex_mode && list_mode) || 
      ((ti_tv || subst_mat_fname != NULL || subst_score_fname != NULL || 
        subst_mat_fname_paml != NULL || scientific_mode) && !list_mode) || 
      (subst_mat_fname != NULL && subst_score_fname != NULL) || 
      (subst_score_fname != NULL && subst_mat_fname_paml != NULL) || 
      (subst_mat_fname != NULL && subst_mat_fname_paml != NULL) || 
      optind != argc - 1) {
    die("ERROR: missing required arguments or illegal combination of arguments.\nTry \"display_rate_matrix -h\" for help.\n");
  }

  F = phast_fopen(argv[optind], "r");
  model = tm_new_from_file(F, 1);

  if (context_ti_tv) {
    /* this option requires completely different handling from the others */
    if (model->order != 2) { 
      die("ERROR: -C requires a model of order 3.\n");
    }
    do_context_dependent_ti_tv(model);
    exit(0);
  }

  if (induced_aa) {
    TreeModel *aa_model = tm_induced_aa(model);
    char *codon_to_aa = get_codon_mapping(model->rate_matrix->states);

    /* before freeing model, grab the expected rate of synonymous
       subst, rho_s */
    for (i = 0; i < model->rate_matrix->size; i++)
      for (j = 0; j < model->rate_matrix->size; j++)
        if (i != j && codon_to_aa[i] == codon_to_aa[j])
          rho_s += mm_get(model->rate_matrix, i, j) * 
            vec_get(model->backgd_freqs, i);

    sfree(codon_to_aa);

    tm_free(model);
    model = aa_model;
  }

  if (all_branches) {
    traversal = tr_inorder(model->tree);
    for (matrix_idx = 0; matrix_idx < lst_size(traversal); matrix_idx++) {
      TreeNode *n = lst_get_ptr(traversal, matrix_idx);
      if (n->parent == NULL) { lst_push_ptr(matrix_list, NULL); continue; }
      M = mm_new(model->rate_matrix->size, model->rate_matrix->states, DISCRETE);
      mm_exp(M, model->rate_matrix, n->dparent);
      lst_push_ptr(matrix_list, M);      
    }
  }
  else if (t >= 0) {
    M = mm_new(model->rate_matrix->size, model->rate_matrix->states, DISCRETE);
    mm_exp(M, model->rate_matrix, t);
    lst_push_ptr(matrix_list, M);
  }
  else 
    lst_push_ptr(matrix_list, model->rate_matrix);

  alph_size = strlen(model->rate_matrix->states);
  nstates = model->rate_matrix->size;

  if (subst_mat_fname != NULL) {
    if ((F = fopen(subst_mat_fname, "r")) == NULL) {
      die("ERROR: Can't open %s.\n", subst_mat_fname);
    }    
    subst_mat = read_subst_mat(F, AA_ALPHABET); 
  }
  else if (subst_mat_fname_paml != NULL) {
    if ((F = fopen(subst_mat_fname_paml, "r")) == NULL) {
      die("ERROR: Can't open %s.\n", subst_mat_fname_paml);
    }    
    subst_mat = read_paml_matrix(F, AA_ALPHABET); 
  }
  else if (subst_score_fname != NULL) {
    if ((F = fopen(subst_score_fname, "r")) == NULL) {
      die("ERROR: Can't open %s.\n", subst_score_fname);
    }    
    subst_mat = read_subst_scores(model, F);
  }
  else if (order1_mod_fname != NULL) {
    if ((F = fopen(order1_mod_fname, "r")) == NULL) {
      die("ERROR: Can't open %s.\n", order1_mod_fname);
    }    
    subst_mat = unproject_rates(model, tm_new_from_file(F, 1));
  }

  /* loop through matrices to print */
  for (matrix_idx = 0; matrix_idx < lst_size(matrix_list); matrix_idx++) {
    M = lst_get_ptr(matrix_list, matrix_idx);

    if (all_branches) {
      if (M == NULL) continue;  /* root */
      printf("BRANCH %d (t = %.6f)\n", ++branch_no,
             ((TreeNode*)lst_get_ptr(traversal, matrix_idx))->dparent);
    }

  /* print no more than 16 columns at a time (except with -a) */
  ncols = (induced_aa ? nstates : 16);
  for (startcol = 0; startcol < nstates; startcol += ncols) {
    endcol = min(nstates, startcol+ncols);

    /* table header */
    if (! list_mode) {
      if (latex_mode) {
        printf("\\begin{tabular}{|c|");
        for (i = startcol; i < endcol; i++) printf("r");
        printf("|}\n\\hline\n");
      }
      printf("%-5s ", "");
      if (latex_mode) printf("& ");
      for (i = startcol; i < endcol; i++) {
        get_state_tuple(model, tuple, i);
        if (latex_mode) {
          printf("{\\bf %s}", tuple);
          if (i < endcol-1) printf("& ");
        }
        else printf("%8s ", tuple);
    }
      if (latex_mode) printf("\\\\\n\\hline\n");
      else printf("\n");
    }

    /* table or list contents */
    for (i = 0; i < nstates; i++) {
      if (induced_aa && AA_ALPHABET[i] == '$' && !do_stop_codons) continue;
      get_state_tuple(model, tuple, i);

      /* get total eq freq of tuples containing CpG dinucs */
      for (k = 0; k < model->order; k++) {
        if (tuple[k] == 'C' && tuple[k+1] == 'G') {
          cpg_eqfreq += vec_get(model->backgd_freqs, i);
/*           printf("***CPG***"); */
          break;
        }
      }

      if (latex_mode) printf("{\\bf %s}& ", tuple);
      else if (!list_mode) printf("%-5s ", tuple);
      for (j = startcol; j < endcol; j++) {
        if (induced_aa && AA_ALPHABET[j] == '$' && !do_stop_codons) continue;
        if (latex_mode) printf("$");
        if (list_mode) {
          if (symmetric && j <= i) continue;
          else if ((t < 0 && ! all_branches) 
		   && (i == j || (!do_zeroes && mm_get(M, i, j) == 0))) 
            continue;
          get_state_tuple(model, tuple2, j);
          printf("%-5s %-5s ", tuple, tuple2);
        }
        if (i == j && suppress_diag && !list_mode) printf("%-7s", "-");
        else { 
	  /* get rate or probability */
	  double val = exch_mode == 0 ? mm_get(M, i, j) : 
	    safediv(mm_get(M, i, j), vec_get(model->backgd_freqs,j));
	  /* print value in format %8.6f or %13.6e */
	  printf(rate_format_string, val); 
	  printf(" ");
	}
        if (latex_mode) {
          printf("$");
          if (j < endcol-1) printf("& ");
        }
        else if (list_mode) {
          int ti, is_cpg;
          if (ti_tv) {
            ti = -1;
            is_cpg = 0;
            for (k = 0; k <= model->order; k++) {
              int dig_i = (i % int_pow(alph_size, k+1)) / int_pow(alph_size, k);
              int dig_j = (j % int_pow(alph_size, k+1)) / int_pow(alph_size, k);
              char next_char = '\0', prev_char = '\0';
              if (dig_i != dig_j) {
                ti = is_transition(M->states[dig_i], M->states[dig_j]);
                if (k != model->order)
                  prev_char = M->states[(i % int_pow(alph_size, k+2)) / 
                                        int_pow(alph_size, k+1)];
                if (k != 0)
                  next_char = M->states[(i % int_pow(alph_size, k)) / 
                                        int_pow(alph_size, k-1)];
                if ((M->states[dig_i] == 'C' && next_char == 'G') || 
                    (M->states[dig_i] == 'G' && prev_char == 'C')) 
                  is_cpg = 1;
              }
            }
	    if (ti == -1)
	      die("ERROR ti=-1\n");
            printf("%5s ", ti ? "ti" : "tv");
/*             printf("%5s ", is_cpg ? "CPG" : "-"); */
            if (ti) {
              total_ti += mm_get(M, i, j) * 
                vec_get(model->backgd_freqs, i);
              if (is_cpg) 
                cpg_ti += mm_get(M, i, j) * 
                  vec_get(model->backgd_freqs, i);
              else non_cpg_ti += mm_get(M, i, j) * 
                     vec_get(model->backgd_freqs, i);
            }
            else {
              total_tv += mm_get(M, i, j) * 
                vec_get(model->backgd_freqs, i);
              if (is_cpg)
                cpg_tv += mm_get(M, i, j) * 
                  vec_get(model->backgd_freqs, i);
              else non_cpg_tv += mm_get(M, i, j) * 
                     vec_get(model->backgd_freqs, i);
            }
          }
          if (subst_mat != NULL) {
            if (mat_get(subst_mat, i, j) == NEGINFTY) 
              printf("%8s", "-"); 
            else printf("%8.4f", mat_get(subst_mat, i, j)); 
          }
          printf("\n");
        }
      }
      if (latex_mode) printf("\\\\\n");
      else if (!list_mode) printf("\n");
    }
    
    /* equilibrium freqs (table case only) */
    if (do_eqfreqs && ! list_mode) {
      if (latex_mode) 
        printf("\\hline\n$\\boldsymbol{\\mathbf{\\pi}}$&");
      else 
        printf("%-5s ", "pi");
      for (i = startcol; i < endcol; i++) {
        if (latex_mode) 
          printf("$%8.4f$ ", vec_get(model->backgd_freqs, i));      
        else 
          printf("%8.4f ", vec_get(model->backgd_freqs, i));      
        if (latex_mode && i < endcol-1) printf("& ");
      }
      if (latex_mode) printf("\\\\\n");
      else printf("\n");
    }

    if (latex_mode) printf("\\hline\n\\end{tabular}\n\n");
  }

  /* equilibrium freqs (list case only) */
  if (do_eqfreqs &&  list_mode) {
    for (i = 0; i < nstates; i++) {
      get_state_tuple(model, tuple, i);
      printf("%-5s %-5s ", "-", tuple); //!!
      printf(rate_format_string, vec_get(model->backgd_freqs, i)); 
      printf("\n");
    }
  }
  
  if (ti_tv && list_mode) {
    printf("\n#Total ti/tv = %.4f\n", total_ti/total_tv);
    printf("#CpG ti ratio = %.4f, CpG tv ratio = %.4f\n", 
           cpg_ti/non_cpg_ti /* * (1 - cpg_eqfreq) */ / cpg_eqfreq, 
           cpg_tv/non_cpg_tv /* * (1 - cpg_eqfreq) */ / cpg_eqfreq);
  }
  else if (induced_aa) 
    printf("\n#Total rho_s/rho_v = %.4f\n", rho_s/(3-rho_s));

  if (all_branches == 1) printf("\n\n");
  }

  tm_free(model);
  lst_free(matrix_list);

  return 0;
}
