/* handling of crispr mutation models for phylogeny reconstruction */
/* reads tab-delimited mutation data, calculates substitution probabilities, custom pruning algorithm for better efficiency */

/* still to do
   x make header file
   x get it to compile
   x test reading and writing, renumbering
   - implement mapping between cell names and leaf ids
   x need to compute distance matrix for init
   - figure out background freqs
   - estimate or pre-estimate relative rates?
   - handle silent states correctly
   - more efficient state traversal
   - exponentiation of rate matrix; make sure freeing any old subst matrices and replacing them.  Use the markov-matrix machinery but replace the exponentiation with something here.  maybe better to keep separate from tree model
   - implement derivative
   - don't forget to come back and clean up scale factor
   - need to handle new subst model and call appropriate special case functions through varPHAST (also update help msg)
   - clean up comments
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <float.h>
#include <phast/stringsplus.h>
#include <phast/crispr.h>

CrisprMutTable *cpr_new_table() {
  CrisprMutTable *retval = malloc(sizeof(CrisprMutTable));
  retval->sitenames = lst_new_ptr(30); /* will resize as needed */
  retval->cellnames = lst_new_ptr(100);
  retval->cellmuts = lst_new_ptr(100);
  retval->nsites = 0;
  retval->ncells = 0;
  retval->nstates = 0;
  retval->eqfreqs = NULL;
  return retval;
}

CrisprMutTable *cpr_read_table(FILE *F) {
  int i, state, lineno = 0;
  List *cols, *muts;
  String *line = str_new(STR_MED_LEN);
  CrisprMutTable *M = cpr_new_table();
  
  if (str_readline(line, F) == EOF)
    die("ERROR: cannot read from file.\n");

  cols = lst_new_ptr(M->nsites + 1);
  str_split(line, NULL, cols);
  str_free(lst_get_ptr(cols, 0)); /* blank */
  for (i = 1; i < lst_size(cols); i++) 
    lst_push_ptr(M->sitenames, lst_get_ptr(cols, i));
  M->nsites = lst_size(M->sitenames); 

  while (str_readline(line, F) != EOF) {
    lineno++;
    str_split(line, NULL, cols);
    if (lst_size(cols) != M->nsites + 1)
      die("ERROR in line %d of input file: number of mutations does not match header\n", lineno);
    lst_push_ptr(M->cellnames, lst_get_ptr(cols, 0));
    muts = lst_new_int(M->nsites);
    lst_push_ptr(M->cellmuts, muts);
    for (i = 1; i< M->nsites + 1; i++) {
      str_as_int(lst_get_ptr(cols, i), &state);
      lst_push_int(muts, state);
      if (state >= M->nstates)
        M->nstates = state + 1;
      str_free(lst_get_ptr(cols, i));
    }
  }

  M->ncells = lst_size(M->cellmuts);
  lst_free(cols);
  str_free(line);

  return M;
}

void cpr_free_table(CrisprMutTable *M) {
  int i;
  for (i = 0; i < lst_size(M->sitenames); i++)
    str_free(lst_get_ptr(M->sitenames, i));
  for (i = 0; i < lst_size(M->cellnames); i++)
    str_free(lst_get_ptr(M->cellnames, i));
  for (i = 0; i < lst_size(M->cellmuts); i++) 
    lst_free(lst_get_ptr(M->cellmuts, i));
  lst_free(M->sitenames);
  lst_free(M->cellnames);
  lst_free(M->cellmuts);
}

void cpr_print_table(CrisprMutTable *M, FILE *F) {
  int i, j;
  for (j = 0; j < M->nsites; j++)
    fprintf(F, "\t%s", ((String*)lst_get_ptr(M->sitenames, j))->chars);
  fprintf(F, "\n");
  for (i = 0; i < M->ncells; i++) {
    fprintf(F, "%s\t", ((String*)lst_get_ptr(M->cellnames, i))->chars);
    for (j = 0; j < M->nsites; j++) 
      fprintf(F, "%d%c", cpr_get_mut(M, i, j),
              j == M->nsites - 1 ? '\n' : '\t');      
  }
}

int cpr_get_mut(CrisprMutTable *M, int cell, int site) {
  List *muts;
  
  if (cell >= M->ncells || site >= M->nsites)
    die("ERROR in cpr_get_mut: index out of range.\n");

  muts = lst_get_ptr(M->cellmuts, cell);
  return lst_get_int(muts, site);
}

void cpr_set_mut(CrisprMutTable *M, int cell, int site, int val) {
  List *muts;
  
  if (cell >= M->ncells || site >= M->nsites)
    die("ERROR in cpr_get_mut: index out of range.\n");
  if (val >= M->nstates)
    die("ERROR in cpr_get_mut: value out of range.\n");
  
  muts = lst_get_ptr(M->cellmuts, cell);
  lst_set_int(muts, site, val);
}

/* renumber states so they fall densely between 0 and nstates - 1,
   with -1 set aside to represent silent states [CHECK] */
void cpr_renumber_states(CrisprMutTable *M) {
  int i, j, newnstates;
  int *statemap = malloc(M->nstates * sizeof(int));
  for (i = 0; i < M->nstates; i++) statemap[i] = -1;
  newnstates = 1;  /* assume 0; ignore and don't count -1 */
  for (i = 0; i < M->ncells; i++) {
    for (j = 0; j < M->nsites; j++) {
      int stateij = cpr_get_mut(M, i, j);
      if (stateij == -1 || stateij == 0) continue;
      else if (stateij >= M->nstates)
        die("ERROR: state number out of range.\n");
      if (statemap[stateij] == -1) 
        statemap[stateij] = newnstates++;
      cpr_set_mut(M, i, j, statemap[stateij]);
    }
  }
  M->nstates = newnstates;
  free(statemap);
}

/* Compute and return the log likelihood of a tree model with respect
   to a CRISPR mutation table.  This is derived from
   nj_compute_log_likelihood but customized for the irreversible
   CRISPR mutation model of Seidel and Stadler (Proc. R. Soc. B
   289:20221844, 2022). If branchgrad is non-null, it will be
   populated with the gradient of the log likelihood with respect to
   the individual branches of the tree, in post-order.  */
double cpr_compute_log_likelihood(TreeModel *mod, CrisprMutTable *M, Vector *branchgrad) {

  int i, j, k, nodeidx, site, cell, state;
  int nstates = M->nstates;
  TreeNode *n, *sibling;
  double total_prob = 0;
  List *traversal, *Pt;
  double **pL = NULL, **pLbar = NULL;
  double log_scale;
  double scaling_threshold = DBL_MIN;
  double ll = 0;
  double tmp[nstates];
  Matrix *grad_mat = NULL;
  MarkovMatrix *par_subst_mat, *sib_subst_mat;
    
  /* set up "inside" probability matrices for pruning algorithm */
  pL = smalloc(nstates * sizeof(double*));
  for (j = 0; j < nstates; j++)
    pL[j] = smalloc((mod->tree->nnodes+1) * sizeof(double));

  /* set up branchwise substitution probability matrices */
  Pt = lst_new_ptr(mod->tree->nnodes);
  for (nodeidx = 0; nodeidx < mod->tree->nnodes; nodeidx++)
    lst_push_ptr(Pt, mm_new(nstates, NULL, DISCRETE));

  /* compute equilibrium freqs if not yet available */
  if (M->eqfreqs == NULL)
    M->eqfreqs = cpr_estim_mut_rates(M, FALSE);
  
  if (branchgrad != NULL) {
    if (branchgrad->size != mod->tree->nnodes-1) /* rooted */
      die("ERROR in nj_compute_log_likelihood: size of branchgrad must be 2n-2\n");
    vec_zero(branchgrad);
    /* set up complementary "outside" probability matrices */
    pLbar = smalloc(nstates * sizeof(double*));
    for (j = 0; j < nstates; j++)
      pLbar[j] = smalloc((mod->tree->nnodes+1) * sizeof(double));
    grad_mat = mat_new(nstates, nstates);
  }

  /* FIXME: set up a mapping from cell numbers to leaf numbers */
  
  cpr_set_subst_matrices(mod, Pt, M->eqfreqs); /* set up P matrices for each edge;
                                               we'll keep this separate from the
                                               TreeModel in this case */

  if (mod->msa_seq_idx == NULL)
    cpr_build_seq_idx(mod, M);
  
  traversal = tr_postorder(mod->tree);
  for (site = 0; site < M->nsites; site++) {
    log_scale = 0;
    for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
      n = lst_get_ptr(traversal, nodeidx);
      if (n->lchild == NULL) {
        /* leaf: base case of recursion */
        cell = mod->msa_seq_idx[n->id]; /* -1 if not found */
        state = cell >= 0 ? cpr_get_mut(M, site, cell) : -1;  /* FIXME: silent states */
        for (i = 0; i < nstates; i++) {
          if (state < 0 || i == state) 
            pL[i][n->id] = 1;
          else
            pL[i][n->id] = 0;
        }
      }
      else {
        /* general recursive case */
        MarkovMatrix *lsubst_mat = lst_get_ptr(Pt, n->lchild->id);
        MarkovMatrix *rsubst_mat = lst_get_ptr(Pt, n->rchild->id);
        for (i = 0; i < nstates; i++) {
          double totl = 0, totr = 0;
          for (j = 0; j < nstates; j++)
            totl += pL[j][n->lchild->id] *
              mm_get(lsubst_mat, i, j);
          
          for (k = 0; k < nstates; k++)
            totr += pL[k][n->rchild->id] *
              mm_get(rsubst_mat, i, k);
          
          if (totl * totr < scaling_threshold) {
            pL[i][n->id] = (totl / scaling_threshold) * totr;
            log_scale -= log(scaling_threshold);
          }
          else {
            pL[i][n->id] = totl * totr;
          }
        }
      }
    }
  
    /* termination */
    total_prob = 0;
    for (i = 0; i < nstates; i++) 
      total_prob += vec_get(M->eqfreqs, i) * pL[i][mod->tree->id];
    
    ll += (log(total_prob) - log_scale);

    assert(isfinite(ll));

    /* to compute gradients efficiently, need to make a second pass
       across the tree to compute "outside" probabilities */
    if (branchgrad != NULL) {
      traversal = tr_preorder(mod->tree);
      assert(log_scale == 0); /* not yet generalized to handle scale factor */

      for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
        n = lst_get_ptr(traversal, nodeidx);

        if (n->parent == NULL) { /* base case */
          for (i = 0; i < nstates; i++)
            pLbar[i][n->id] = vec_get(M->eqfreqs, i);
        }
        else {            /* recursive case */
          sibling = (n == n->parent->lchild ?
                               n->parent->rchild : n->parent->lchild);
          par_subst_mat = lst_get_ptr(Pt, n->id);
          sib_subst_mat = lst_get_ptr(Pt, sibling->id);
          	  
          for (j = 0; j < nstates; j++) { /* parent state */
            tmp[j] = 0;
            for (k = 0; k < nstates; k++)  /* sibling state */
              tmp[j] += pLbar[j][n->parent->id] *
                pL[k][sibling->id] * mm_get(sib_subst_mat, j, k);
          }
          
          for (i = 0; i < nstates; i++) { /* child state */
            pLbar[i][n->id] = 0;
            for (j = 0; j < nstates; j++)  /* parent state */
              pLbar[i][n->id] +=
                tmp[j] * mm_get(par_subst_mat, j, i);
          }
        }
      }

      /* now compute branchwise derivatives in a final pass */
      traversal = mod->tree->nodes;
      for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
        TreeNode *par;
        double base_prob = total_prob, deriv = 0;
        
        n = lst_get_ptr(traversal, nodeidx);
        par = n->parent;
	
        if (par == NULL || n == mod->tree->rchild) 
          continue;
       
        sibling = (n == n->parent->lchild ?
                   n->parent->rchild : n->parent->lchild);

        sib_subst_mat = lst_get_ptr(Pt, sibling->id);
        
        /* this part is just a constant to propagate through to the
           derivative */
        for (i = 0; i < nstates; i++) {  /* parent */
          tmp[i] = 0;
          for (k = 0; k < nstates; k++)  /* sibling */
            tmp[i] += pL[k][sibling->id] * mm_get(sib_subst_mat, i, k);
        }
        
        /* calculate derivative analytically */
        cpr_branch_grad(grad_mat, n->dparent, M->eqfreqs); 
        for (i = 0; i < nstates; i++)   
          for (j = 0; j < nstates; j++)    
            deriv +=  tmp[i] * pLbar[i][par->id] * pL[j][n->id] * mat_get(grad_mat, i, j);  

        deriv *= 1.0 / base_prob; /* because need deriv of log P */

        vec_set(branchgrad, nodeidx, vec_get(branchgrad, nodeidx) + deriv );
      }
    }
  }
  
  for (j = 0; j < nstates; j++)
    sfree(pL[j]);
  sfree(pL);

  for (nodeidx = 0; nodeidx < mod->tree->nnodes; nodeidx++)
    mm_free(lst_get_ptr(Pt, nodeidx));
  lst_free(Pt);
                 
  if (branchgrad != NULL) {
    for (j = 0; j < nstates; j++)
      sfree(pLbar[j]);
    sfree(pLbar);
    mat_free(grad_mat);
  }

  return ll;
}

/* build and return an upper triangular distance matrix for the cells
   in a CrisprMutTable using a simple Poisson-type distance measure */
Matrix *cpr_compute_dist(CrisprMutTable *M) {
  int i, j;
  Matrix *retval = mat_new(M->ncells, M->ncells);

  mat_zero(retval);
  
  for (i = 0; i < M->ncells; i++) 
    for (j = i+1; j < M->ncells; j++) 
      mat_set(retval, i, j,
              cpr_compute_pw_dist(M, i, j));

  return retval;  
}

/* compute pairwise distance between two cells using Poisson-type
   measure */
double cpr_compute_pw_dist(CrisprMutTable *M, int i, int j) {
  int k, diff = 0, n = 0;
  double d;
  for (k = 0; k < M->nsites; k++) {
    int typei = cpr_get_mut(M, i, k),
      typej = cpr_get_mut(M, j, k);

    if (typei == -1 || typej == -1)
      continue;
    n++;
    if (typei != typej)
      diff++;
  }
  d = -log(1.0 - diff*1.0/n);
  /* assumes mutations arise by a Poisson process with rate one; this
     is the mle for the time elapsed in units of expected mutations
     per site */

  if (n == 0 || d > 3)
    d = 3; /* set a max to keep the initialization reasonable */

  return d;
}

/* set a substitution matrix for each edge based on current branch
   lengths and sitewise rate parameters */
void cpr_set_subst_matrices(TreeModel *mod, List *Pt, Vector *eqfreqs) {
  for (int nodeidx = 0; nodeidx < mod->tree->nnodes; nodeidx++) {
    TreeNode *n = lst_get_ptr(mod->tree->nodes, nodeidx);
    MarkovMatrix *P = lst_get_ptr(Pt, nodeidx);
    if (n->parent == NULL) continue;
    cpr_set_branch_matrix(P, n->dparent, eqfreqs); 
  }    
}

void cpr_set_branch_matrix(MarkovMatrix *P, double t, Vector *eqfreqs) {
  /* FIXME: make sure other code expects silent state last */
  /* should already have zeroes for all but the defined rates */

  /* what does the scale factor mean? how to compute? */

  /* silent state is the last one */
  int j, silst = P->size - 1;
  double scale = 0.0; /* FIXME: pass in */
  
  /* substitution probabilities from 0 (unedited) state to all edited
     (and not silent) states */
  for (j = 1; j < silst; j++)
    mm_set(P, 0, j, vec_get(eqfreqs, j) * exp(-t * scale) * (1 - exp(-t)));
  mm_set(P, 0, 0, exp(-t*(1+scale)));
  
  /* substitution probabilities from edited states to themselves */
  for (j = 1; j < silst; j++)
    mm_set(P, j, j, exp(-t * scale));

  /* substitution probabilities to silent state */
  for (j = 0; j < silst; j++)
    mm_set(P, j, silst, 1 - exp(-t * scale));
  mm_set(P, silst, silst, 1); /* absorbing state */
}

/* compute gradients of elements of substitution matrix with respect
   to branch length */
void cpr_branch_grad(Matrix *grad, double t, Vector *eqfreqs) {
  int j, silst = grad->ncols - 1;
  double scale = 0.0; /* FIXME: pass in */
  mat_zero(grad);
  
  /* derivatives of substitution probabilities from 0 (unedited) state
     to all edited (and not silent) states */
  for (j = 1; j < silst; j++)
    mat_set(grad, 0, j, vec_get(eqfreqs, j) *
            (-scale * exp(-t*scale) * (1 - exp(-t)) +
             exp(-t * (1+scale))));
  mat_set(grad, 0, 0, -(1+scale) * exp(-t*(1+scale)));
  
  /* derivatives of substitution probabilities from edited states to
     themselves */
  for (j = 1; j < silst; j++)
    mat_set(grad, j, j, -scale * exp(-t * scale));

  /* derivatives of substitution probabilities to silent state */
  for (j = 0; j < silst; j++)
    mat_set(grad, j, silst, scale * exp(-t * scale));
}

/* estimate relative mutation rates based on relative frequencies in
   data set.  Returns a vector of size M->nstates + 1, with the last
   element representing the relative frequency of the silent state.
   If ignore_silent == TRUE, that frequency will be forced to zero.
   To ensure all other states have nonzero values, it may be helpful to
   preprocess with cpr_renumber_states */
Vector *cpr_estim_mut_rates(CrisprMutTable *M, unsigned int ignore_silent) {
  Vector *retval = vec_new(M->nstates + 1);
  int silentstate = M->nstates;
  vec_zero(retval);
  for (int i = 0; i < M->ncells; i++) {
    for (int j = 0; j < M->ncells; j++) {
      int mut = cpr_get_mut(M, i, j);
      if (mut == -1 && ignore_silent == FALSE)
        vec_set(retval, silentstate, vec_get(retval, silentstate) + 1.0);
      else
        vec_set(retval, mut, vec_get(retval, mut) + 1.0);
    }
  }
  pv_normalize(retval);
  return retval;
}

/*  Build index of leaf ids to cell indices based on matching names.
    Leaves not present in the alignment will be ignored.  Also, it's
    not required that there's a leaf for every sequence in the
    alignment. This is a version of tm_build_seq_idx adapted to work
    with a CrisprMutTable */
void cpr_build_seq_idx(TreeModel *mod, CrisprMutTable *M) {
  int i, idx;  
  mod->msa_seq_idx = smalloc(mod->tree->nnodes * sizeof(int));
  /* let's just reuse this even though it's misnamed for the purpose */
  
  for (i = 0; i < mod->tree->nnodes; i++) {
    TreeNode *n = lst_get_ptr(mod->tree->nodes, i);
    mod->msa_seq_idx[i] = -1;
    if (n->lchild == NULL && n->rchild == NULL) {
      String *namestr = str_new_charstr(n->name);
      if (str_in_list_idx(namestr, M->cellnames, &idx) == 1)
        mod->msa_seq_idx[i] = idx;
      str_free(namestr);
    }
  }
}
