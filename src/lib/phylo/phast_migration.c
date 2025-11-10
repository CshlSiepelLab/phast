/* handling of migration models for phylogeny reconstruction */

#include <stdio.h>
#include <stdlib.h>
#include <phast/misc.h>
#include <phast/stringsplus.h>
#include <phast/hashtable.h>
#include <phast/crispr.h>
#include <phast/migration.h>
#include <phast/multiDAG.h>

MigTable *mig_new() {
  MigTable *M = smalloc(sizeof(MigTable));
  M->cellnames = lst_new_ptr(200);
  M->states = lst_new_int(50);
  M->statehash = hsh_new(50);
  M->statenames = lst_new_ptr(50);
  M->ncells = 0;
  M->nstates = 0;
  M->nparams = 0;
  M->gtr_params = NULL;
  M->deriv_gtr = NULL;
  M->Pt = NULL;
  M->backgd_freqs = NULL;
  M->rate_matrix = NULL;
  M->rate_matrix_param_row = NULL;
  M->rate_matrix_param_col = NULL;
  return M;
}

void mig_free(MigTable *M) {
  lst_free_strings(M->cellnames);
  lst_free(M->cellnames);
  lst_free(M->states);
  hsh_free(M->statehash);
  lst_free_strings(M->statenames);
  lst_free(M->statenames);
  sfree(M);
}

/* read migration table from file; expects two comma-delimited
   columns: cell name, state name */
MigTable *mig_read_table(FILE *F) {
  MigTable *M = mig_new();
  String *line = str_new(STR_MED_LEN);
  List *cols = lst_new_ptr(2);
  int lineno = 0;
  while (str_readline(line, F) != EOF) {
    lineno++;

    if (str_starts_with_charstr(line, "#")) /* comment line */
      continue;
    
    str_split(line, ",", cols);
    if (lst_size(cols) != 2)
      die("ERROR in line %d of input file: each line must have two columns (comma-delimited)\n", lineno);
    lst_push_ptr(M->cellnames, lst_get_ptr(cols, 0));

    String *statename = lst_get_ptr(cols, 1);
    int stateno = hsh_get_int(M->statehash, statename->chars);
    if (stateno == -1) { /* new state */
      stateno = M->nstates++;
      hsh_put_int(M->statehash, statename->chars, stateno);
      lst_push_ptr(M->statenames, str_dup(statename));      
    }
    
    lst_push_int(M->states, stateno);
  }

  lst_free(cols);
  str_free(line);

  M->ncells = lst_size(M->cellnames);
  mig_update_states(M);
  return M;
}

void mig_update_states(MigTable *M) {
  M->nparams = (M->nstates * (M->nstates - 1)) / 2;
  M->gtr_params = vec_new(M->nparams);
  vec_set_random(M->gtr_params, 1.0, 0.1);
  M->deriv_gtr = vec_new(M->nparams);
  vec_zero(M->deriv_gtr);
  M->backgd_freqs = vec_new(M->nstates); /* has to be done before below */
  vec_set_all(M->backgd_freqs, 1.0 / M->nstates);
  M->rate_matrix_param_row = (List**)smalloc(M->nparams * sizeof(List*));
  M->rate_matrix_param_col = (List**)smalloc(M->nparams * sizeof(List*));
  for (int i = 0; i < M->nparams; i++) {
    M->rate_matrix_param_row[i] = lst_new_int(2);
    M->rate_matrix_param_col[i] = lst_new_int(2);
  }
  M->rate_matrix = mm_new(M->nstates, NULL, DISCRETE);
  mm_set_eigentype(M->rate_matrix, REAL_NUM);
  /* FIXME: send in state labels */
  mig_set_REV_matrix(M, M->gtr_params);
}

/* check that migration table contains the same list of cellnames as a
   mutation table. Also sort migration table cells to match mutation
   table */
void mig_check_table(MigTable *mg, CrisprMutTable *mm) {
  int i;
  if (mg->ncells != mm->ncells)
    die("ERROR: migration table and mutation table have different numbers of "
        "cells.\n");

  /* build hashtable for mutation table cellnames */
  Hashtable *mutnamehash = hsh_new(mg->ncells * 2);
  for (i = 0; i < mm->ncells; i++) {
    String *s = lst_get_ptr(mm->cellnames, i);
    hsh_put_int(mutnamehash, s->chars, i);    
  }

  /* now check that all migration table cellnames are in mutation table */
  List *new_cellnames = lst_new_ptr(mg->ncells);
  List *new_states = lst_new_int(mg->ncells);
  for (i = 0; i < mg->ncells; i++) {
    String *s = lst_get_ptr(mg->cellnames, i);
    int idx = hsh_get_int(mutnamehash, s->chars);
    if (idx == -1)
      die("ERROR: cell '%s' in migration table not found in mutation table.\n",
          s->chars);
    lst_push_ptr(new_cellnames, str_dup(s));
    lst_push_int(new_states, lst_get_int(mg->states, i));
  }
  
  /* replace cellnames and states with sorted versions */
  lst_free_strings(mg->cellnames);
  lst_free(mg->cellnames);
  lst_free(mg->states);
  mg->cellnames = new_cellnames;
  mg->states = new_states;
}

/* compute log likelihood of migration model based on a given tree
   model and migration table for tips.  If branchgrad is non-null, it
   will be populated with the gradient with respect to the individual
   branches of the tree, in post-order */
double mig_compute_log_likelihood(TreeModel *mod, MigTable *mg, 
                                  CrisprMutModel *cprmod, Vector *branchgrad) {

  int i, j, k, nodeidx, cell, state;
  int nstates = mg->nstates; 
  TreeNode *n, *sibling;
  double total_prob = 0;
  List *traversal, *pre_trav;
  double **pL = NULL, **pLbar = NULL;
  double scaling_threshold = DBL_MIN * 1.0e10;  /* need some padding */
  double lscaling_threshold = log(scaling_threshold), ll = 0;
  double tmp[nstates], root_eqfreqs[nstates];
  Matrix **grad_mat = NULL;
  List **grad_mat_P;
  MarkovMatrix *par_subst_mat, *sib_subst_mat, *leading_Pt, *lsubst_mat, *rsubst_mat; ;
  Vector *lscale, *lscale_o; /* inside and outside versions */
  Vector *this_deriv_gtr = NULL;
  unsigned int rescale;

  /* set up "inside" probability matrices for pruning algorithm */
  pL = smalloc(nstates * sizeof(double*));
  for (j = 0; j < nstates; j++)
    pL[j] = smalloc((cprmod->mod->tree->nnodes+1) * sizeof(double));

  /* we also need to keep track of the log scale of every node for
     underflow purposes */
  lscale = vec_new(mod->tree->nnodes+1); 
  lscale_o = vec_new(mod->tree->nnodes+1); 
  vec_zero(lscale); vec_zero(lscale_o);
  
  if (branchgrad != NULL) {
    if (branchgrad->size != mod->tree->nnodes-1) /* rooted */
      die("ERROR in mig_compute_log_likelihood: size of branchgrad must be 2n-2\n");
    vec_zero(branchgrad);
    /* set up complementary "outside" probability matrices */
    pLbar = smalloc(nstates * sizeof(double*));
    for (j = 0; j < nstates; j++)
      pLbar[j] = smalloc((mod->tree->nnodes+1) * sizeof(double));

    vec_zero(mg->deriv_gtr);
    grad_mat = malloc(mod->tree->nnodes * sizeof(Matrix*));
    for (j = 0; j < mod->tree->nnodes; j++)
      grad_mat[j] = mat_new(nstates, nstates);
    /* each node of the tree needs a list of gradient
       matrices, one for each free GTR parameter */
    grad_mat_P = malloc(mod->tree->nnodes * sizeof(void*));
    for (j = 0; j < mod->tree->nnodes; j++) {
      grad_mat_P[j] = lst_new_ptr(mg->gtr_params->size);
      for (int jj = 0; jj < mg->gtr_params->size; jj++)
        lst_push_ptr(grad_mat_P[j],
                     mat_new(nstates, nstates));
    }
    this_deriv_gtr = vec_new(mg->gtr_params->size);
  }

  mig_update_subst_matrices(mod, mg); /* compute all necessary migration probability
                                       matrices */
  if (cprmod->mod->msa_seq_idx == NULL)
    cpr_build_seq_idx(cprmod->mod, cprmod->mut);

  traversal = tr_postorder(cprmod->mod->tree);
    
  /* this model allows a leading branch to the root of the tree but
     forces the unedited state at the start of that branch.  We can
     simulate this behavior by setting the root eq freqs equal to
     the conditional distribution at the end of the branch given the
     unedited state at the start */
  leading_Pt = lst_get_ptr(mg->Pt, mod->tree->id); 
  for (i = 0; i < nstates; i++)
    root_eqfreqs[i] = mm_get(leading_Pt, 0, i);
    
  for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
    n = lst_get_ptr(traversal, nodeidx);

    if (n->lchild == NULL) {
      /* leaf: base case of recursion */
      cell = cprmod->mod->msa_seq_idx[n->id]; /* CHECK.  okay to use this version? */
      
      if (cell == -1)
        die("ERROR in mig_compute_log_likelihood: leaf '%s' not found in "
            "migration table.\n", n->name);

      state = lst_get_int(mg->states, cell);
      assert(state >= 0 && state < nstates);
      
      for (i = 0; i < nstates; i++) {
        if (i == state)
          pL[i][n->id] = 1;
        else
          pL[i][n->id] = 0;
      }
    }
    else {
      /* general recursive case */
      lsubst_mat = lst_get_ptr(mg->Pt, n->lchild->id);
      rsubst_mat = lst_get_ptr(mg->Pt, n->rchild->id);

      rescale = FALSE;
      for (i = 0; i < nstates; i++) {
        double totl = 0, totr = 0;
        for (j = 0; j < nstates; j++)
          totl += pL[j][n->lchild->id] *
            mm_get(lsubst_mat, i, j);
          
        for (k = 0; k < nstates; k++)
          totr += pL[k][n->rchild->id] *
            mm_get(rsubst_mat, i, k);

        pL[i][n->id] = totl * totr;
        if (totl > 0 && totr > 0 && pL[i][n->id] < scaling_threshold)
          rescale = TRUE;
      }

      /* deal with nodewise scaling */
      vec_set(lscale, n->id, vec_get(lscale, n->lchild->id) +
              vec_get(lscale, n->rchild->id));
      if (rescale == TRUE) { /* have to rescale for all states */
        vec_set(lscale, n->id, vec_get(lscale, n->id) + lscaling_threshold);
        for (i = 0; i < nstates; i++) 
          pL[i][n->id] /= scaling_threshold;
      }
    }
  }
  /* termination */
  total_prob = 0;
  for (i = 0; i < nstates; i++)
    total_prob += root_eqfreqs[i] *
      pL[i][mod->tree->id] * root_eqfreqs[i];
    
  ll += (log(total_prob) + vec_get(lscale, mod->tree->id));

  /* to compute gradients efficiently, need to make a second pass
     across the tree to compute "outside" probabilities */
  if (branchgrad != NULL) {
    double expon;
    pre_trav = tr_preorder(mod->tree);

    for (nodeidx = 0; nodeidx < lst_size(pre_trav); nodeidx++) {
      n = lst_get_ptr(pre_trav, nodeidx);

      if (n->parent == NULL) { /* base case */
        for (i = 0; i < nstates; i++)
          pLbar[i][n->id] = root_eqfreqs[i];
      }
      else {            /* recursive case */
        sibling = (n == n->parent->lchild ? n->parent->rchild : n->parent->lchild);
        par_subst_mat = lst_get_ptr(mg->Pt, n->id);
        sib_subst_mat = lst_get_ptr(mg->Pt, sibling->id);
        rescale = FALSE;
          
        for (j = 0; j < nstates; j++) { /* parent state */
          tmp[j] = 0;
          for (k = 0; k < nstates; k++)  /* sibling state */
            tmp[j] += pLbar[j][n->parent->id] *
              pL[k][sibling->id] * mm_get(sib_subst_mat, j, k);
        }
          
        for (i = 0; i < nstates; i++) { /* child state */
          pLbar[i][n->id] = 0;
          for (j = 0; j < nstates; j++) { /* parent state */
            pLbar[i][n->id] +=
              tmp[j] * mm_get(par_subst_mat, j, i);
            if (tmp[j] > 0 && mm_get(par_subst_mat, j, i) > 0 &&
                pLbar[i][n->id] < scaling_threshold)
              rescale = TRUE;
          }
        }
        vec_set(lscale_o, n->id, vec_get(lscale_o, n->parent->id) +
                vec_get(lscale, sibling->id));
        if (rescale == TRUE) { /* rescale for all states */
          vec_set(lscale_o, n->id, vec_get(lscale_o, n->id) + lscaling_threshold);
          for (i = 0; i < nstates; i++)
            pLbar[i][n->id] /= scaling_threshold;     
        }
      }
    }

    /* TEMPORARY: check inside/outside */
    /* for (nodeidx = 0; nodeidx < lst_size(mod->tree->nodes); nodeidx++) { */
    /*   double pr = 0; */
    /*   n = lst_get_ptr(mod->tree->nodes, nodeidx); */
    /*   assert(vec_get(lscale, n->id) == 0 && vec_get(lscale_o, n->id) == 0); */
    /*   for (j = 0; j < nstates; j++) */
    /*     pr += pL[j][n->id] * pLbar[j][n->id]; */
    /*   printf("Tuple %d, node %d: %f (%f)\n", tupleidx, nodeidx, log(pr), log(total_prob)); */
    /* } */

    /* now compute branchwise derivatives in a final pass */
    for (nodeidx = 0; nodeidx < lst_size(mod->tree->nodes); nodeidx++) {
      TreeNode *par, *sibling;
      double base_prob = total_prob, deriv;
        
      n = lst_get_ptr(mod->tree->nodes, nodeidx);
      par = n->parent;
	
      if (par == NULL) 
        continue;
       
      sibling = (n == n->parent->lchild ?
                 n->parent->rchild : n->parent->lchild);

      sib_subst_mat = lst_get_ptr(mg->Pt, sibling->id);

      /* this part is just a constant to propagate through to the
         derivative */
      for (i = 0; i < nstates; i++) {  /* parent */
        tmp[i] = 0;
        for (k = 0; k < nstates; k++)  /* sibling */
          tmp[i] += pL[k][sibling->id] * mm_get(sib_subst_mat, i, k);
      }

      if (n != mod->tree->rchild) { /* skip branch to right of root because unrooted */
        /* calculate derivative analytically */
        deriv = 0;
        mig_grad_REV_dt(mg, grad_mat[n->id], n->dparent); /* FIXME: customize */
          
        for (i = 0; i < nstates; i++)   
          for (j = 0; j < nstates; j++)    
            deriv +=  tmp[i] * pLbar[i][par->id] * pL[j][n->id] * mat_get(grad_mat[n->id], i, j);

        /* adjust for all relevant scale terms; do everything in log space */
        expon = -vec_get(lscale, mod->tree->id)
          + vec_get(lscale, sibling->id) + vec_get(lscale_o, par->id)
          + vec_get(lscale, n->id) - log(base_prob);
        /* note division by base_prob because we need deriv of log P */

        /* avoid overflow */
        if (expon > 700.0) expon = 700.0;
        if (expon < -745.0) expon = -745.0;
          
        deriv *= exp(expon);
        assert(isfinite(deriv));
                  
        vec_set(branchgrad, n->id, vec_get(branchgrad, n->id) + deriv);
      }

      /* we need partial derivatives for migration rates also;
         they have to be aggregated across all branches */
      mig_grad_REV_dr(mg, grad_mat_P[n->id], n->dparent);  /* FIXME: customize */
      /* loop over rate parameters */
      for (int pidx = 0; pidx < mg->gtr_params->size; pidx++) {
        double pderiv = 0; /* partial deriv wrt this param */
        Matrix *dP_dr = lst_get_ptr(grad_mat_P[n->id], pidx);
        for (int i = 0; i < nstates; i++) {
          for (int j = 0; j < nstates; j++) {
            pderiv += tmp[i] * pLbar[i][par->id] * pL[j][n->id] *
              mat_get(dP_dr, i, j);
          }
        }
        vec_set(this_deriv_gtr, pidx, pderiv);
      }
      /* adjust for all relevant scale terms */
      vec_scale(this_deriv_gtr, exp(expon));
      vec_plus_eq(mg->deriv_gtr, this_deriv_gtr);
    }
  }
  
  for (j = 0; j < nstates; j++)
    sfree(pL[j]);
  sfree(pL);
  
  if (branchgrad != NULL) {
    for (j = 0; j < nstates; j++)
      sfree(pLbar[j]);
    sfree(pLbar);
    for (j = 0; j < mod->tree->nnodes; j++)      
      mat_free(grad_mat[j]);
    free(grad_mat);
    for (j = 0; j < mod->tree->nnodes; j++) {
      List *gmats = grad_mat_P[j];
      for (int jj = 0; jj < lst_size(gmats); jj++)
        mat_free(lst_get_ptr(gmats, jj));
      lst_free(gmats);
    }
    free(grad_mat_P);
    vec_free(this_deriv_gtr);
  }

  vec_free(lscale);
  vec_free(lscale_o);
  
  return ll;
}

/* set P = exp(Qt) for each branch in the tree model */
void mig_update_subst_matrices(TreeModel *mod, MigTable *mg) {
  if (mg->Pt == NULL) { /* first time through */
    mg->Pt = lst_new_ptr(mod->tree->nnodes);
    for (int nodeidx = 0; nodeidx < mod->tree->nnodes; nodeidx++) {
      MarkovMatrix *P = mm_new(mg->nstates, NULL, DISCRETE);
      mm_set_eigentype(P, REAL_NUM);
      lst_push_ptr(mg->Pt, P);
    }
  }
  for (int nodeidx = 0; nodeidx < mod->tree->nnodes; nodeidx++) {
    TreeNode *n = lst_get_ptr(mod->tree->nodes, nodeidx);
    MarkovMatrix *P = lst_get_ptr(mg->Pt, nodeidx);
    mm_exp(P, mg->rate_matrix, n->dparent);  /* CHECK: diagonalization, rate_matrix updated, etc */
  }
}

/* helper function to sample states at internal nodes.  Given an array of
   doubles of dimension nstates, sample in proportion to those values and return
   an index between 0 and nstates-1.  Original probabilities do not need to be
   normalized */
static inline
int mig_sample_state(double *probs, int nstates) {
  double total = 0, cumprob = 0;
  for (int i = 0; i < nstates; i++)
    total += probs[i];
  assert(total > 0.0);
  double r = unif_rand() * total;
  for (int i = 0; i < nstates; i++) {
    cumprob += probs[i];
    if (r <= cumprob)
      return i;
  }
  return nstates - 1; /* should not get here */
}

/* given a tree model and migration table, sample states at internal
   nodes of the tree.  Populates a list parallel to mod->tree->nodes
   whose elements are state indices */
void mig_sample_states(TreeModel *mod, MigTable *mg, 
                       CrisprMutModel *cprmod, List *state_samples) {

  int i, j, k, nodeidx, cell, state;
  int nstates = mg->nstates; 
  TreeNode *n;
  List *traversal, *pre_trav;
  double **pL = NULL;
  double scaling_threshold = DBL_MIN * 1.0e10;  /* need some padding */
  double lscaling_threshold = log(scaling_threshold);
  double root_eqfreqs[nstates];
  MarkovMatrix *par_subst_mat, *leading_Pt, *lsubst_mat, *rsubst_mat; ;
  Vector *lscale; /* inside and outside versions */
  unsigned int rescale;

  /* set up "inside" probability matrices */
  pL = smalloc(nstates * sizeof(double*));
  for (j = 0; j < nstates; j++)
    pL[j] = smalloc((cprmod->mod->tree->nnodes+1) * sizeof(double));

  /* we also need to keep track of the log scale of every node for
     underflow purposes */
  lscale = vec_new(mod->tree->nnodes+1); 
  vec_zero(lscale);
  
  mig_update_subst_matrices(mod, mg); /* compute all necessary migration probability
                                         matrices */
  if (cprmod->mod->msa_seq_idx == NULL)
    cpr_build_seq_idx(cprmod->mod, cprmod->mut);

  traversal = tr_postorder(cprmod->mod->tree);
    
  /* this model allows a leading branch to the root of the tree but
     forces the unedited state at the start of that branch.  We can
     simulate this behavior by setting the root eq freqs equal to
     the conditional distribution at the end of the branch given the
     unedited state at the start */
  leading_Pt = lst_get_ptr(mg->Pt, mod->tree->id); 
  for (i = 0; i < nstates; i++)
    root_eqfreqs[i] = mm_get(leading_Pt, 0, i);
    
  for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
    n = lst_get_ptr(traversal, nodeidx);

    if (n->lchild == NULL) {
      /* leaf: base case of recursion */
      cell = cprmod->mod->msa_seq_idx[n->id];
      
      if (cell == -1)
        die("ERROR in mig_compute_log_likelihood: leaf '%s' not found in "
            "migration table.\n", n->name);

      state = lst_get_int(mg->states, cell);
      assert(state >= 0 && state < nstates);
      
      for (i = 0; i < nstates; i++) {
        if (i == state)
          pL[i][n->id] = 1;
        else
          pL[i][n->id] = 0;
      }

      /* we can set the state samples for leaves now */
      lst_set_int(state_samples, n->id, state);
    }
    else {
      /* general recursive case */
      lsubst_mat = lst_get_ptr(mg->Pt, n->lchild->id);
      rsubst_mat = lst_get_ptr(mg->Pt, n->rchild->id);

      double maxP = 0;
      for (i = 0; i < nstates; i++) {
        double totl = 0, totr = 0;
        for (j = 0; j < nstates; j++)
          totl += pL[j][n->lchild->id] *
            mm_get(lsubst_mat, i, j);
          
        for (k = 0; k < nstates; k++)
          totr += pL[k][n->rchild->id] *
            mm_get(rsubst_mat, i, k);

        pL[i][n->id] = totl * totr;
        if (maxP < pL[i][n->id])
          maxP = pL[i][n->id];
      }
      rescale = FALSE;
      if (maxP > 0 && maxP < scaling_threshold)
        rescale = TRUE;

      /* deal with nodewise scaling */
      vec_set(lscale, n->id, vec_get(lscale, n->lchild->id) +
              vec_get(lscale, n->rchild->id));
      if (rescale == TRUE) { /* have to rescale for all states */
        vec_set(lscale, n->id, vec_get(lscale, n->id) + lscaling_threshold);
        for (i = 0; i < nstates; i++) 
          pL[i][n->id] /= scaling_threshold;
      }
    }
  }

  /* Now pass from root to leaves and sample notes based on smpled parent and
     inside probabilities */
  pre_trav = tr_preorder(mod->tree);
  double *sampdens = smalloc(nstates * sizeof(double));
  for (nodeidx = 0; nodeidx < lst_size(pre_trav); nodeidx++) {
    n = lst_get_ptr(pre_trav, nodeidx);

    if (n->parent == NULL) { /* base case */
      for (i = 0; i < nstates; i++)
        sampdens[i] = root_eqfreqs[i] * pL[i][n->id];
      state = mig_sample_state(sampdens, nstates);
      lst_set_int(state_samples, n->id, state);
    }
    else if (n->lchild == NULL)  /* leaf: already handled */
      continue;
    else { /* recursive case */
      int parstate = lst_get_int(state_samples, n->parent->id);
      par_subst_mat = lst_get_ptr(mg->Pt, n->id);

      for (i = 0; i < nstates; i++)
        sampdens[i] = pL[i][n->id] * mm_get(par_subst_mat, parstate, i);;
      state = mig_sample_state(sampdens, nstates);

      lst_set_int(state_samples, n->id, state);
    }
  }
  
  for (j = 0; j < nstates; j++)
    sfree(pL[j]);
  sfree(pL);
  
  vec_free(lscale);
  sfree(sampdens);
}

/* based on a tree model and list of states at all nodes, obtain a
   multigraph representing migration events and their times */
struct mdag *mig_get_graph(TreeModel *mod, MigTable *mg, List *state_samples) {
  MultiDAG *g = mdag_new(mg);

  /* first compute height of each node */
  List *traversal = tr_postorder(mod->tree);
  Vector *heights = vec_new(mod->tree->nnodes);
  vec_set_all(heights, 0.0);
  for (int nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
    TreeNode *n = lst_get_ptr(traversal, nodeidx);
    if (n->lchild == NULL) /* leaf */
      vec_set(heights, n->id, 0.0);
    else { /* mostly concerned with ultrametric trees here but we'll
              take the max of both children to be safe */
      double lht = vec_get(heights, n->lchild->id) + n->lchild->dparent;
      double rht = vec_get(heights, n->rchild->id) + n->rchild->dparent;
      vec_set(heights, n->id, (lht > rht ? lht : rht));
    }
  }  
 
  for (int nodeidx = 0; nodeidx < mod->tree->nnodes; nodeidx++) {
    TreeNode *n = lst_get_ptr(mod->tree->nodes, nodeidx);
    if (n->parent == NULL)
      continue;
    int childstate = lst_get_int(state_samples, n->id);
    int parstate = lst_get_int(state_samples, n->parent->id);
    if (childstate != parstate) {
      double start_time = vec_get(heights, n->parent->id);
      double end_time = vec_get(heights, n->id);
      mdag_add_edge(g, parstate, childstate, start_time, end_time); 
    }
  }
  return g;
}
  
/***************************************************************************
 Functions adapted from phast_subst_mods.c to handle migration models 
 ***************************************************************************/
void mig_set_REV_matrix(MigTable *mg, Vector *params) {
  int i, j, start_idx = 0;
  if (mg->backgd_freqs == NULL)
    die("mig_set_REV_matrix: mg->backgd_freqs is NULL\n");
  for (i = 0; i < mg->rate_matrix->size; i++) {
    double rowsum = 0;
    for (j = i+1; j < mg->rate_matrix->size; j++) {
      double val;
      val = vec_get(params, start_idx);
      mm_set(mg->rate_matrix, i, j,
             val * vec_get(mg->backgd_freqs, j));
      mm_set(mg->rate_matrix, j, i,
             val * vec_get(mg->backgd_freqs, i));
      rowsum += (val * vec_get(mg->backgd_freqs, j));
      lst_clear(mg->rate_matrix_param_row[start_idx]);
      lst_clear(mg->rate_matrix_param_col[start_idx]);
      lst_push_int(mg->rate_matrix_param_row[start_idx], i);
      lst_push_int(mg->rate_matrix_param_col[start_idx], j);
      lst_push_int(mg->rate_matrix_param_row[start_idx], j);
      lst_push_int(mg->rate_matrix_param_col[start_idx], i);

      start_idx++;
    }
    for (j = 0; j < i; j++)
      rowsum += mm_get(mg->rate_matrix, i, j);
    mm_set(mg->rate_matrix, i, i, -1 * rowsum);
  }
  mig_scale_rate_matrix(mg);
  mm_diagonalize(mg->rate_matrix);
}

void mig_scale_rate_matrix(MigTable *mg) {
  double scale = 0;
  for (int i = 0; i < mg->rate_matrix->size; i++) {
    double rowsum = 0;
    for (int j = 0; j < mg->rate_matrix->size; j++) 
      if (j != i) rowsum += mm_get(mg->rate_matrix, i, j);
    scale += (vec_get(mg->backgd_freqs, i) * rowsum);
  }
  mm_scale(mg->rate_matrix, 1.0/scale);
}

/* dP/dt for REV */
void mig_grad_REV_dt(MigTable *mg, Matrix *grad, double t) {
  double g_ij, s_ik, sprime_kj, lambda_k;
  MarkovMatrix *Q = mg->rate_matrix;

  if (t == 0) {
    mat_copy(grad, Q->matrix); /* at t=0, dP/dt = Q */    
    return;
  }
  
  if (Q->evec_matrix_r == NULL || Q->evals_r == NULL ||
      Q->evec_matrix_inv_r == NULL) {
    mm_diagonalize(Q);
    if (Q->diagonalize_error) {
      mat_print(Q->matrix, stderr);
      die("ERROR in mig_grad_REV_dt: rate matrix could not be diagonalized.\n");
    }
  }
  for (int i = 0; i < Q->size; i++) {
    for (int j = 0; j < Q->size; j++) {
      g_ij = 0;
      for (int k = 0; k < Q->size; k++) {
        s_ik = mat_get(Q->evec_matrix_r, i, k);
        sprime_kj = mat_get(Q->evec_matrix_inv_r, k, j);
        lambda_k = vec_get(Q->evals_r, k);
        g_ij += s_ik * lambda_k * exp(lambda_k * t) * sprime_kj;
        /* see Siepel & Hausser 2004 Eq C.10 but in this case we omit
           the denominator and avoid summing over elements (handled by
           calling code) */
      }
      mat_set(grad, i, j, g_ij);
    }
  }  
}

/* dP/dr for REV, where r is a free parameter corresponding in the
   rate matrix.  Fills dQ_dr, a list of matrices, the ith element of
   which corresponds to the ith free parameter */
void mig_grad_REV_dr(MigTable *mg, List *dP_dr_lst, double t) {

  /* this code is adapted from compute_grad_em_exact in
     phast_fit_em.c */
  int i, j, k, l, m, idx, lidx, orig_size;
  int nstates = mg->nstates;
  List *erows = lst_new_int(4), *ecols = lst_new_int(4), 
    *distinct_rows = lst_new_int(2);

  double **tmpmat = smalloc(nstates * sizeof(double *));
  double **sinv_dq_s = smalloc(nstates * sizeof(double *));
  double **dq = smalloc(nstates * sizeof(double *));
  double **f = smalloc(nstates * sizeof(double *));
  for (i = 0; i < nstates; i++) {
    dq[i] = smalloc(nstates * sizeof(double));
    tmpmat[i] = smalloc(nstates * sizeof(double));
    sinv_dq_s[i] = smalloc(nstates * sizeof(double));
    f[i] = smalloc(nstates * sizeof(double));
  }

  MarkovMatrix *Q = mg->rate_matrix;

  if (Q->evec_matrix_r == NULL || Q->evals_r == NULL ||
      Q->evec_matrix_inv_r == NULL) {
    mm_diagonalize(Q);
    if (Q->diagonalize_error) {
      mat_print(Q->matrix, stderr);
      die("ERROR in mig_grad_REV_dr: rate matrix could not be diagonalized.\n");
    }
  }
  
  for (idx = 0; idx < lst_size(dP_dr_lst); idx++) {

    Matrix *dP_dr = (Matrix *)lst_get_ptr(dP_dr_lst, idx);
    mat_zero(dP_dr);
    
    for (i = 0; i < nstates; i++) 
      for (j = 0; j < nstates; j++) 
        dq[i][j] = tmpmat[i][j] = sinv_dq_s[i][j] = 0;

    /* element coords (rows/col pairs) at which current param appears in Q */
    lst_cpy(erows, mg->rate_matrix_param_row[idx]);
    lst_cpy(ecols, mg->rate_matrix_param_col[idx]);
    if (lst_size(erows) != lst_size(ecols))
      die("ERROR mig_grad_REV_dr: size of erows (%i) does not match size of "
          "ecols (%i)\n",
          lst_size(erows), lst_size(ecols));

    /* set up dQ, the partial deriv of Q wrt the current param */
    lst_clear(distinct_rows);
    for (i = 0, orig_size = lst_size(erows); i < orig_size; i++) {
      l = lst_get_int(erows, i); 
      m = lst_get_int(ecols, i);

      if (dq[l][m] != 0)    /* row/col pairs should be unique */
        die("ERROR tm_grad_REV_dr: dq[%i][%i] should be zero but is %f\n",
	    l, m, dq[l][m]);

      dq[l][m] = vec_get(mg->backgd_freqs, m);      
      if (dq[l][m] == 0) continue; 

      /* keep track of distinct rows and cols with non-zero entries */
      /* also add diagonal elements to 'rows' and 'cols' lists, as
         necessary */
      if (dq[l][l] == 0) {      /* new row */
        lst_push_int(distinct_rows, l);
        lst_push_int(erows, l);
        lst_push_int(ecols, l);
      }

      dq[l][l] -= dq[l][m]; /* row sums to zero */
    }

    /* compute S^-1 dQ S */
    for (lidx = 0; lidx < lst_size(erows); lidx++) {
      i = lst_get_int(erows, lidx);
      k = lst_get_int(ecols, lidx);
      for (j = 0; j < nstates; j++)
        tmpmat[i][j] += mat_get(Q->evec_matrix_r, k, j) * dq[i][k];
    }

    for (lidx = 0; lidx < lst_size(distinct_rows); lidx++) {
      k = lst_get_int(distinct_rows, lidx);
      for (i = 0; i < nstates; i++) {
        for (j = 0; j < nstates; j++) {
          sinv_dq_s[i][j] += mat_get(Q->evec_matrix_inv_r, i, k) * tmpmat[k][j];
        }
      }
    }

    /* set up Schadt and Lange's F matrix */
    for (i = 0; i < nstates; i++) {
      for (j = 0; j < nstates; j++) {
        if (fabs(vec_get(Q->evals_r, i) - vec_get(Q->evals_r, j)) < 1e-8)
          f[i][j] = exp((vec_get(Q->evals_r, i)) * t) * t;
        else
          f[i][j] = (exp((vec_get(Q->evals_r, i)) * t) 
                     - exp((vec_get(Q->evals_r, j)) * t)) /
            ((vec_get(Q->evals_r, i)) - (vec_get(Q->evals_r, j)));
      }
    }

    /* compute (F o S^-1 dQ S) S^-1 */
    for (i = 0; i < nstates; i++) {
      for (j = 0; j < nstates; j++) {
        tmpmat[i][j] = 0;
        for (k = 0; k < nstates; k++) 
          tmpmat[i][j] += f[i][k] * sinv_dq_s[i][k] *
            mat_get(Q->evec_matrix_inv_r, k, j);
      }
    }

    /* compute S (F o S^-1 dQ S) S^-1 */
    for (i = 0; i < nstates; i++) {
      for (j = 0; j < nstates; j++) {
        double partial_p = 0;        
        for (k = 0; k < nstates; k++) 
          partial_p += mat_get(Q->evec_matrix_r, i, k) * tmpmat[k][j];
        mat_set(dP_dr, i, j, partial_p);
      }
    }
  }

  /* free allocated memory */
  for (i = 0; i < nstates; i++) {
    sfree(dq[i]); sfree(tmpmat[i]); sfree(sinv_dq_s[i]); sfree(f[i]);
  }
  sfree(dq); sfree(tmpmat); sfree(sinv_dq_s); sfree(f);
  lst_free(erows); lst_free(ecols); lst_free(distinct_rows);
}
