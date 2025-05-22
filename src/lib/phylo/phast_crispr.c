/* handling of crispr mutation models for phylogeny reconstruction */
/* reads tab-delimited mutation data, calculates substitution probabilities, custom pruning algorithm for better efficiency */

/* still to do
   x make header file
   x get it to compile
   x test reading and writing, renumbering
   - implement mapping between cell names and leaf ids
   - need to compute distance matrix for init
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
      die("ERROR in line %d of input file: number of mutations does not match header\n");
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
  TreeNode *n;
  double total_prob = 0;
  List *traversal;
  double **pL = NULL, **pLbar = NULL;
  double log_scale;
  double scaling_threshold = DBL_MIN;
  double ll = 0;
  double tmp[nstates];
  Matrix *grad_mat = NULL;

  /* set up "inside" probability matrices for pruning algorithm */
  pL = smalloc(nstates * sizeof(double*));
  for (j = 0; j < nstates; j++)
    pL[j] = smalloc((mod->tree->nnodes+1) * sizeof(double));

  /* set up branchwise substitution probability matrices */
  Pt = lst_new_ptr(mod->tree->nnodes);
  for (nodeidx = 0; nodeidx < mod->tree->nnodes; nodeidx++)
    lst_push_ptr(Pt, mm_new(nstates, NULL, DISCRETE));
  
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
  
  cpr_set_subst_matrices(mod, Pt); /* set up P matrices for each edge;
                                      we'll keep this separate from the
                                      TreeModel in this case */

  traversal = tr_postorder(mod->tree);
  for (site = 0; site < M->nsites; site++) {
    log_scale = 0;
    for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
      n = lst_get_ptr(traversal, nodeidx);
      if (n->lchild == NULL) {
        /* leaf: base case of recursion */
        cell = 5; /*XXXXXX();*/ /* FIXME: map from leaf index to cell number */
        state = cpr_get_mut(M, site, cell);
        for (i = 0; i < nstates; i++) {
          if (i == state)  /* FIXME: what to do with silent states? -1 */
            pL[i][n->id] = 1;
          else
            pL[i][n->id] = 0;
        }
      }
      else {
        /* general recursive case */
        MarkovMatrix *lsubst_mat = mod->P[n->lchild->id];
        MarkovMatrix *rsubst_mat = mod->P[n->rchild->id];
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
    for (i = 0; i < nstates; i++)  /* FIXME: how to handle background freqs? */
      total_prob += vec_get(mod->backgd_freqs, i) *
        pL[i][mod->tree->id] * mod->freqK[0];
    
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
            pLbar[i][n->id] = vec_get(mod->backgd_freqs, i);   /* FIXME */
        }
        else {            /* recursive case */
          TreeNode *sibling = (n == n->parent->lchild ?
                               n->parent->rchild : n->parent->lchild);
          MarkovMatrix *par_subst_mat = mod->P[n->id];
          MarkovMatrix *sib_subst_mat = mod->P[sibling->id];
          	  
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
        TreeNode *par, *sib;
        double base_prob = total_prob, deriv = 0;
        
        n = lst_get_ptr(traversal, nodeidx);
        par = n->parent;
	
        if (par == NULL || n == mod->tree->rchild) 
          continue;
       
        sib = (n == n->parent->lchild ?
               n->parent->rchild : n->parent->lchild);

        /* this part is just a constant to propagate through to the
           derivative */
        for (i = 0; i < nstates; i++) {  /* parent */
          tmp[i] = 0;
          for (k = 0; k < nstates; k++)  /* sibling */
            tmp[i] += pL[k][sib->id] * mm_get(mod->P[sib->id], i, k);
        }
        
        /* calculate derivative analytically */
        cpr_branch_grad(grad_mat, n->dparent, eqfreqs); 
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

  if (d > 3)
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
  
  /* substitution probabilities from 0 (unedited) state to all edited
     (and not silent) states */
  for (j = 1; j < silst; j++)
    mm_set(0, j, vec_get(eqfreq, j) * exp(-t * scale) * (1 - exp(-t)));
  mm_set(0, 0, exp(-t*(1+scale)));
  
  /* substitution probabilities from edited states to themselves */
  for (j = 1; j < silst; j++)
    mm_set(j, j, exp(-t * scale));

  /* substitution probabilities to silent state */
  for (j = 0, j < silst, j++)
    mm_set(j, silst, 1 - exp(-t * scale));
  mm_set(silst, silst, 1); /* absorbing state */
}

/* compute gradients of elements of substitution matrix with respect
   to branch length */
cpr_branch_grad(Matrix *grad, double t, Vector *eqfreqs) {
  int j, silst = P->size - 1;
  mat_zero(grad);
  
  /* derivatives of substitution probabilities from 0 (unedited) state
     to all edited (and not silent) states */
  for (j = 1; j < silst; j++)
    mat_set(grad, 0, j, vec_get(eqfreq, j) *
            (-scale * exp(-t*scale) * (1 - exp(-t)) +
             exp(-t * (1+scale))));
  mm_set(0, 0, -(1+scale) * exp(-t*(1+scale)));
  
  /* derivatives of substitution probabilities from edited states to
     themselves */
  for (j = 1; j < silst; j++)
    mm_set(j, j, -scale * exp(-t * scale));

  /* derivatives of substitution probabilities to silent state */
  for (j = 0, j < silst, j++)
    mm_set(j, silst, scale * exp(-t * scale));
}
