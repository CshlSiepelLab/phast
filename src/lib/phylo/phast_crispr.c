/* handling of crispr mutation models for phylogeny reconstruction */
/* reads tab-delimited mutation data, calculates substitution probabilities, custom pruning algorithm for better efficiency */


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
  retval->sitewise_nstates = NULL;
  return retval;
}

CrisprMutTable *cpr_copy_table(CrisprMutTable *orig) {
  CrisprMutTable *copy = cpr_new_table();
  int i, j;

  copy->nsites = orig->nsites;
  copy->ncells = orig->ncells;
  copy->nstates = orig->nstates;
  
  for (i = 0; i < lst_size(orig->sitenames); i++) 
    lst_push_ptr(copy->sitenames, str_dup(lst_get_ptr(orig->sitenames, i)));
  for (i = 0; i < lst_size(orig->cellnames); i++) 
    lst_push_ptr(copy->cellnames, str_dup(lst_get_ptr(orig->cellnames, i)));
  for (i = 0; i < orig->ncells; i++) {
    List *muts = lst_new_int(orig->nsites);
    lst_push_ptr(copy->cellmuts, muts);
    for (j = 0; j < orig->nsites; j++) 
      lst_push_int(muts, cpr_get_mut(orig, i, j));
  }  

  return copy;
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
  fprintf(F, "cell");
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
   with -1 set aside to represent silent states */
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
   to a CRISPR mutation table.  This function is derived from
   nj_compute_log_likelihood but is customized for the irreversible
   CRISPR mutation model of Seidel and Stadler (Proc. R. Soc. B
   289:20221844, 2022). If branchgrad is non-null, it will be
   populated with the gradient of the log likelihood with respect to
   the individual branches of the tree, in post-order.  */
double cpr_compute_log_likelihood(CrisprMutModel *cprmod, Vector *branchgrad) {

  int i, j, k, nodeidx, site, cell, state;
  int nstates;
  TreeNode *n, *sibling;
  double total_prob = 0;
  List *traversal, *pre_trav;
  double **pL = NULL, **pLbar = NULL;
  double scaling_threshold = DBL_MIN * 1.0e10;  /* need some padding */
  double lscaling_threshold = log(scaling_threshold), ll = 0;
  double tmp[cprmod->nstates+1], root_eqfreqs[cprmod->nstates+1];
  Matrix *grad_mat = NULL;
  MarkovMatrix *par_subst_mat, *sib_subst_mat;
  static CrisprAncestralStateSets *ancsets = NULL;
  Vector *lscale, *lscale_o; /* inside and outside versions */
  unsigned int rescale;
  List *par_states, *lchild_states, *rchild_states, *child_states, *sib_states;
  int pstate, lcstate, rcstate, cstate, sstate;
      
  /* set up "inside" probability matrices for pruning algorithm */
  pL = smalloc((cprmod->nstates+1) * sizeof(double*));
  for (j = 0; j < (cprmod->nstates+1); j++)
    pL[j] = smalloc((cprmod->mod->tree->nnodes+1) * sizeof(double));

  /* we also need to keep track of the log scale of every node for
     underflow purposes */
  lscale = vec_new(cprmod->mod->tree->nnodes+1); 
  lscale_o = vec_new(cprmod->mod->tree->nnodes+1); 
  
  if (branchgrad != NULL) {
    if (branchgrad->size != cprmod->mod->tree->nnodes-1) /* rooted */
      die("ERROR in cpr_compute_log_likelihood: size of branchgrad must be 2n-2\n");
    vec_zero(branchgrad);
    /* set up complementary "outside" probability matrices */
    pLbar = smalloc((cprmod->nstates+1) * sizeof(double*));
    for (j = 0; j < (cprmod->nstates+1); j++)
      pLbar[j] = smalloc((cprmod->mod->tree->nnodes+1) * sizeof(double));
    cprmod->deriv_sil = 0.0;
  }

  cprmod->mod->tree->dparent = cprmod->leading_t; /* update length of leading branch */
  
  cpr_update_model(cprmod); /* compute all necessary mutation probability matrices */

  if (cprmod->mod->msa_seq_idx == NULL)
    cpr_build_seq_idx(cprmod->mod, cprmod->mut);

  /* also set up ancestral state sets if not already available */
  if (ancsets == NULL)
    ancsets = cpr_new_state_sets(cprmod->mod->tree->nnodes);
  if (cprmod->mod->tree->nnodes > ancsets->nnodes) /* in case number of nodes grows */
    cpr_state_sets_resize(ancsets, cprmod->mod->tree->nnodes);
  
  traversal = tr_postorder(cprmod->mod->tree);
  for (site = 0; site < cprmod->nsites; site++) {
    int silst;
    List *Pt = lst_get_ptr(cprmod->Pt, site);
    MarkovMatrix *leading_Pt;
    double this_deriv_sil;
    
    nstates = cprmod->mut->sitewise_nstates[site] + 1; /* have to allow for silent state */
    silst = nstates - 1; /* silent state will always be last */

    /* first zero out all pL values because with the smart
       algorithm, we won't visit most elements in the matrix */
    for (nodeidx = 0; nodeidx < cprmod->mod->tree->nnodes; nodeidx++) {
      ancsets->nodetypes[nodeidx] = -99; /* also initialize these */
      for (i = 0; i < nstates; i++) 
        pL[i][nodeidx] = 0;
    }

    /* also reset scale */
    vec_zero(lscale); vec_zero(lscale_o);

    /* this model allows a leading branch to the root of the tree but
       forces the unedited state at the start of that branch.  We can
       simulate this behavior by setting the root eq freqs equal to
       the conditional distribution at the end of the branch given the
       unedited state at the start */
    leading_Pt = lst_get_ptr(Pt, cprmod->mod->tree->id); 
    for (i = 0; i < nstates; i++)
      root_eqfreqs[i] = mm_get(leading_Pt, 0, i);
    
    for (nodeidx = 0; nodeidx < lst_size(traversal); nodeidx++) {
      int mut;

      n = lst_get_ptr(traversal, nodeidx);

      if (n->lchild == NULL) {
        /* leaf: base case of recursion */
        cell = cprmod->mod->msa_seq_idx[n->id];
        if (cell == -1)
          die("ERROR in cpr_compute_log_likelihood: leaf '%s' not found in mutation table.\n",
              n->name);

        mut = cpr_get_mut(cprmod->mut, cell, site);
        state = (mut == -1 ? silst : mut);
        assert(state >= 0 && state <= silst);
        pL[state][n->id] = 1;

        /* also update nodetype */
        if (state < silst) /* ancestral (0) or derived edit */
          ancsets->nodetypes[n->id] = state;
        else /* unrestricted if silent */
          ancsets->nodetypes[n->id] = ancsets->NORESTRICT;
      }
      else {
        /* general recursive case */
        MarkovMatrix *lsubst_mat = lst_get_ptr(Pt, n->lchild->id);
        MarkovMatrix *rsubst_mat = lst_get_ptr(Pt, n->rchild->id);
        int lchildtype, rchildtype, thistype;
        
        /* first set nodetype based on nodetypes of children */
        lchildtype = ancsets->nodetypes[n->lchild->id];
        rchildtype = ancsets->nodetypes[n->rchild->id];
        
        if (lchildtype == 0 || rchildtype == 0) /* if either child is
                                                   unedited, parent
                                                   must be unedited */
          thistype = 0;
        else if (lchildtype != ancsets->NORESTRICT &&
                 rchildtype != ancsets->NORESTRICT &&
                 lchildtype != rchildtype) /* if children have
                                              different edits, parent
                                              must be unedited */
          thistype = 0;
        else if (lchildtype == rchildtype && 
                 lchildtype != ancsets->NORESTRICT) /* if children have
                                                       same edits,
                                                       parent must
                                                       share it or be
                                                       unedited */
          
          thistype = lchildtype;
        else if (lchildtype != ancsets->NORESTRICT &&
                 rchildtype == ancsets->NORESTRICT) /* if one edited
                                                       child, parent
                                                       must have same
                                                       edit or be
                                                       unedited */
          thistype = lchildtype;
        
        else if (rchildtype != ancsets->NORESTRICT &&
                 lchildtype == ancsets->NORESTRICT) /* converse case */
          thistype = rchildtype;

        else /* otherwise we have to consider all possible states */
          thistype = ancsets->NORESTRICT;

        ancsets->nodetypes[n->id] = thistype;

        /* now get corresponding sets of eligible states */
        par_states = cpr_get_state_set(ancsets, n, nstates);
        lchild_states = cpr_get_state_set(ancsets, n->lchild, nstates);
        rchild_states = cpr_get_state_set(ancsets, n->rchild, nstates);
        
        rescale = FALSE;
        for (i = 0; i < lst_size(par_states); i++) {
          double totl = 0, totr = 0;
          pstate = lst_get_int(par_states, i);
          for (j = 0; j < lst_size(lchild_states); j++) {
            lcstate = lst_get_int(lchild_states, j);
            totl += pL[lcstate][n->lchild->id] *
              mm_get(lsubst_mat, pstate, lcstate);
          }
          for (k = 0; k < lst_size(rchild_states); k++) {
            rcstate = lst_get_int(rchild_states, k);
            totr += pL[rcstate][n->rchild->id] *
              mm_get(rsubst_mat, pstate, rcstate);
          }
        
          pL[pstate][n->id] = totl * totr;
          if (totl > 0 && totr > 0 && pL[pstate][n->id] < scaling_threshold) 
            rescale = TRUE;
        }

        /* deal with nodewise scaling */
        vec_set(lscale, n->id, vec_get(lscale, n->lchild->id) +
                vec_get(lscale, n->rchild->id));
        if (rescale == TRUE) { /* have to rescale for all states */
          vec_set(lscale, n->id, vec_get(lscale, n->id) - lscaling_threshold);
          for (i = 0; i < lst_size(par_states); i++) {
            pstate = lst_get_int(par_states, i);
            pL[pstate][n->id] /= scaling_threshold;
          }
        }
      }
    }
  
    /* termination */
    par_states = cpr_get_state_set(ancsets, cprmod->mod->tree, nstates);
    total_prob = 0;
    for (i = 0; i < lst_size(par_states); i++) {
      int rstate = lst_get_int(par_states, i);
      total_prob += root_eqfreqs[rstate] * pL[rstate][cprmod->mod->tree->id];
    }    
    ll += (log(total_prob) - vec_get(lscale, cprmod->mod->tree->id));

    if (!isfinite(ll)) break;  /* this is possible with the crispr
                                  model (total_prob == 0); in this
                                  case we'll bail out and return
                                  -inf */
  
    /* to compute gradients efficiently, need to make a second pass
       across the tree to compute "outside" probabilities */
    if (branchgrad != NULL) {  
      pre_trav = tr_preorder(cprmod->mod->tree);

      for (nodeidx = 0; nodeidx < lst_size(pre_trav); nodeidx++) {
        n = lst_get_ptr(pre_trav, nodeidx);

        if (n->parent == NULL) { /* base case */
          par_states = cpr_get_state_set(ancsets, n, nstates);
          for (i = 0; i < lst_size(par_states); i++) {
            pstate = lst_get_int(par_states, i);
            pLbar[pstate][n->id] = root_eqfreqs[pstate];
          }
        }
        else {            /* recursive case */
          sibling = (n == n->parent->lchild ?
                     n->parent->rchild : n->parent->lchild);
          par_subst_mat = lst_get_ptr(Pt, n->id);
          sib_subst_mat = lst_get_ptr(Pt, sibling->id);

          rescale = FALSE;
          par_states = cpr_get_state_set(ancsets, n->parent, nstates);
          child_states = cpr_get_state_set(ancsets, n, nstates);
          sib_states = cpr_get_state_set(ancsets, sibling, nstates);
          for (j = 0; j < lst_size(par_states); j++) { /* parent state */
            pstate = lst_get_int(par_states, j);
            tmp[pstate] = 0;
            for (k = 0; k < lst_size(sib_states); k++) { /* sibling state */
              sstate = lst_get_int(sib_states, k);
              tmp[pstate] += pLbar[pstate][n->parent->id] *
                pL[sstate][sibling->id] * mm_get(sib_subst_mat, pstate, sstate);
            }
          }
          
          for (i = 0; i < lst_size(child_states); i++) { /* child state */
            cstate = lst_get_int(child_states, i);
            pLbar[cstate][n->id] = 0;
            for (j = 0; j < lst_size(par_states); j++) { /* parent state */
              pstate = lst_get_int(par_states, j);
              pLbar[cstate][n->id] +=
                tmp[pstate] * mm_get(par_subst_mat, pstate, cstate);
              if (tmp[pstate] > 0 && mm_get(par_subst_mat, pstate, cstate) > 0 &&
                  pLbar[cstate][n->id] < scaling_threshold)
                rescale = TRUE;
            }
          }

          vec_set(lscale_o, n->id, vec_get(lscale_o, n->parent->id) +
                  vec_get(lscale, sibling->id));
          if (rescale == TRUE) { /* rescale for all states */
            vec_set(lscale_o, n->id, vec_get(lscale_o, n->id) - lscaling_threshold);
            for (i = 0; i < lst_size(child_states); i++)
              pLbar[lst_get_int(child_states, i)][n->id] /= scaling_threshold;     
          }
        }
      }

      /* TEMPORARY: check inside/outside */
      /* for (nodeidx = 0; nodeidx < lst_size(cprmod->mod->tree->nodes); nodeidx++) { */
      /*   double pr = 0; */
      /*   n = lst_get_ptr(cprmod->mod->tree->nodes, nodeidx); */
      /*   assert(vec_get(lscale, n->id) == 0 && vec_get(lscale_o, n->id) == 0); */
      /*   for (j = 0; j < nstates; j++) */
      /*     pr += pL[j][n->id] * pLbar[j][n->id]; */
      /*   printf("Site %d, node %d: %f (%f)\n", site, nodeidx, log(pr), log(total_prob)); */
      /* } */
        
      
      /* now compute branchwise derivatives in a final pass */
      grad_mat = mat_new(nstates, nstates);
      for (nodeidx = 0; nodeidx < lst_size(cprmod->mod->tree->nodes); nodeidx++) {
        TreeNode *par;
        double base_prob = total_prob, deriv;
        
        n = lst_get_ptr(cprmod->mod->tree->nodes, nodeidx);
        par = n->parent;
	
        if (par == NULL) 
          continue;
       
        sibling = (n == n->parent->lchild ?
                   n->parent->rchild : n->parent->lchild);

        sib_subst_mat = lst_get_ptr(Pt, sibling->id);

        /* get corresponding sets of eligible states */
        par_states = cpr_get_state_set(ancsets, par, nstates);
        child_states = cpr_get_state_set(ancsets, n, nstates);
        sib_states = cpr_get_state_set(ancsets, sibling, nstates);
        
        /* this part is just a constant to propagate through to the
           derivative */
        for (i = 0; i < lst_size(par_states); i++) {  /* parent */
          pstate = lst_get_int(par_states, i);
          tmp[pstate] = 0;
          for (k = 0; k < lst_size(sib_states); k++) { /* sibling */
            sstate = lst_get_int(sib_states, k);
            tmp[pstate] += pL[sstate][sibling->id] * mm_get(sib_subst_mat, pstate, sstate);
          }
        }

        if (n != cprmod->mod->tree->rchild) { /* skip this for right branch from root because unrooted */
          /* calculate derivative analytically */
          deriv = 0;
          cpr_branch_grad(grad_mat, n->dparent, cprmod->sil_rate,
                          lst_get_ptr(cprmod->sitewise_mutrates, site));
          for (i = 0; i < lst_size(par_states); i++) {
            pstate = lst_get_int(par_states, i);
            for (j = 0; j < lst_size(child_states); j++) {
              cstate = lst_get_int(child_states, j);
              deriv +=  tmp[pstate] * pLbar[pstate][par->id] * pL[cstate][n->id] *
                mat_get(grad_mat, pstate, cstate);
            }
          }

          /* adjust for all relevant scale terms; do everything in log space */
          deriv *= exp(vec_get(lscale, cprmod->mod->tree->id)
                       - vec_get(lscale, sibling->id) - vec_get(lscale_o, par->id)
                       - vec_get(lscale, n->id) - log(base_prob));
          /* note division by base_prob because we need deriv of log P */
        
          vec_set(branchgrad, n->id, vec_get(branchgrad, n->id) + deriv);
        }
        
        /* now do the same for the derivative wrt the silent
           rate; has to be aggregated across branches */
        this_deriv_sil = 0;
        cpr_silent_rate_grad(grad_mat, n->dparent, cprmod->sil_rate,
                             lst_get_ptr(cprmod->sitewise_mutrates, site));
        for (i = 0; i < lst_size(par_states); i++) {
          pstate = lst_get_int(par_states, i);
          for (j = 0; j < lst_size(child_states); j++) {
            cstate = lst_get_int(child_states, j);
            this_deriv_sil +=  tmp[pstate] * pLbar[pstate][par->id] * pL[cstate][n->id] *
              mat_get(grad_mat, pstate, cstate);
          }
        }

        /* adjust for all relevant scale terms; do everything in log space */
        this_deriv_sil *= exp(vec_get(lscale, cprmod->mod->tree->id)
                              - vec_get(lscale, sibling->id) - vec_get(lscale_o, par->id)
                              - vec_get(lscale, n->id) - log(base_prob));
        cprmod->deriv_sil += this_deriv_sil;        
      }


      /* also compute gradient for leading branch */
      child_states = cpr_get_state_set(ancsets, cprmod->mod->tree, nstates);
      cpr_branch_grad(grad_mat, cprmod->mod->tree->dparent, cprmod->sil_rate,
                      lst_get_ptr(cprmod->sitewise_mutrates, site));
      cprmod->deriv_leading_t = 0;
      for (j = 0; j < lst_size(child_states); j++) {
        cstate = lst_get_int(child_states, j);
        cprmod->deriv_leading_t += pL[cstate][cprmod->mod->tree->id]
          * mat_get(grad_mat, 0, cstate);
      }
      /* rescale */
      cprmod->deriv_leading_t *= exp(vec_get(lscale, cprmod->mod->tree->id) - log(total_prob)); 

      /* leading branch also contributes to derivative of silent rate */
      this_deriv_sil = 0;
      cpr_silent_rate_grad(grad_mat, cprmod->mod->tree->dparent, cprmod->sil_rate,
                           lst_get_ptr(cprmod->sitewise_mutrates, site));
      for (j = 0; j < lst_size(child_states); j++) {
        cstate = lst_get_int(child_states, j);
        this_deriv_sil +=  pL[cstate][cprmod->mod->tree->id]
          * mat_get(grad_mat, 0, cstate);
      }
      this_deriv_sil *= exp(vec_get(lscale, cprmod->mod->tree->id) - log(total_prob));
      cprmod->deriv_sil += this_deriv_sil;        
      
      mat_free(grad_mat);
    }
  }
  
  for (j = 0; j < cprmod->nstates+1; j++)
    sfree(pL[j]);
  sfree(pL);

  if (branchgrad != NULL) {
    for (j = 0; j < cprmod->nstates+1; j++)
      sfree(pLbar[j]);
    sfree(pLbar);
  }

  vec_free(lscale);
  vec_free(lscale_o);
  
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
void cpr_set_subst_matrices(TreeModel *mod, double silent_rate,
                            List *Pt, Vector *mutrates) {
  for (int nodeidx = 0; nodeidx < mod->tree->nnodes; nodeidx++) {
    TreeNode *n = lst_get_ptr(mod->tree->nodes, nodeidx);
    MarkovMatrix *P = lst_get_ptr(Pt, nodeidx);
    cpr_set_branch_matrix(P, n->dparent, silent_rate, mutrates); 
  }    
}

/* set P = exp(Qt) matrix for branch length t, using parameterization
   of Mai, Chu, and Raphael, doi:10.1101/2024.03.05.583638 */  
void cpr_set_branch_matrix(MarkovMatrix *P, double t, double silent_rate, Vector *mutrates) {  
  int j, silst = P->size - 1; /* silent state is the last one */
  mat_zero(P->matrix);
  
  /* substitution probabilities from 0 (unedited) state to all edited
     (and not silent) states */
  for (j = 1; j < silst; j++)
    mm_set(P, 0, j, vec_get(mutrates, j) * exp(-t * silent_rate) * (1 - exp(-t)));
  mm_set(P, 0, 0, exp(-t*(1+silent_rate)));
  
  /* substitution probabilities from edited states to themselves */
  for (j = 1; j < silst; j++)
    mm_set(P, j, j, exp(-t * silent_rate));

  /* substitution probabilities to silent state */
  for (j = 0; j < silst; j++)
    mm_set(P, j, silst, 1 - exp(-t * silent_rate));
  mm_set(P, silst, silst, 1); /* absorbing state */
}

/* compute gradients of elements of substitution matrix with respect
   to branch length */
void cpr_branch_grad(Matrix *grad, double t, double silent_rate, Vector *mutrates) {
  int j, silst = grad->nrows - 1; 
  mat_zero(grad);
         
  /* derivatives of substitution probabilities from 0 (unedited) state
     to all edited (and not silent) states */
  for (j = 1; j < silst; j++)
    mat_set(grad, 0, j, vec_get(mutrates, j) *
            (-silent_rate * exp(-t*silent_rate) * (1 - exp(-t)) +
             exp(-t * (1+silent_rate))));
  mat_set(grad, 0, 0, -(1+silent_rate) * exp(-t*(1+silent_rate)));
  
  /* derivatives of substitution probabilities from edited states to
     themselves */
  for (j = 1; j < silst; j++)
    mat_set(grad, j, j, -silent_rate * exp(-t * silent_rate));

  /* derivatives of substitution probabilities to silent state */
  for (j = 0; j < silst; j++)
    mat_set(grad, j, silst, silent_rate * exp(-t * silent_rate));
}

/* compute gradients of elements of substitution matrix with respect
   to the silencing rate */
void cpr_silent_rate_grad(Matrix *grad, double t, double silent_rate, Vector *mutrates) {
  int j, silst = grad->nrows - 1;
  double Es = exp(-silent_rate * t);
  double E1 = exp(-t);    
  mat_zero(grad);
         
  /* derivatives of substitution probabilities from 0 (unedited) state
     to all edited (and not silent) states */
  for (j = 1; j < silst; j++)
    mat_set(grad, 0, j, vec_get(mutrates, j) * -t * Es * (1.0-E1));
  mat_set(grad, 0, 0, -t * exp(-t*(1.0+silent_rate)));
  
  /* derivatives of substitution probabilities from edited states to
     themselves */
  for (j = 1; j < silst; j++)
    mat_set(grad, j, j, -t * Es);

  /* derivatives of substitution probabilities to silent state */
  for (j = 0; j < silst; j++)
    mat_set(grad, j, silst, t * Es);
}

/* estimate relative mutation rates based on relative frequencies in
   data set.  Returns a vector of size M->nstates.  To ensure all
   other states have nonzero values, it may be helpful to preprocess
   with cpr_renumber_states */
Vector *cpr_estim_mutrates(CrisprMutTable *M,
                            enum crispr_mutrates_type type) {
  int i, j;
  Vector *retval = vec_new(M->nstates);
  vec_zero(retval);

  if (type == UNIF) { /* uniform distrib over non-zero states */
    for (i = 1; i < M->nstates; i++)
      vec_set(retval, i, 1.0/(M->nstates-1.0));
    return retval;
  }

  /* otherwise compute empirically */
  for (i = 0; i < M->ncells; i++) {
    for (j = 0; j < M->nsites; j++) {
      int mut = cpr_get_mut(M, i, j);
      if (mut > 0)
        vec_set(retval, mut, vec_get(retval, mut) + 1.0);
    }
  }
  pv_normalize(retval);
  return retval;
}

/* estimate relative mutation rates based on relative frequencies in
   data set, but separately for each site.  Returns a list of
   vectors. Assumes mutation matrix has been re-indexed using
   cpr_new_sitewise_table */
List *cpr_estim_sitewise_mutrates(CrisprMutTable *M,
                                  enum crispr_mutrates_type type) {
  int i, j, mut;
  List *sitewise_mutrates = lst_new_ptr(M->nsites);
  for (j = 0; j < M->nsites; j++) {
    Vector *f = vec_new(M->sitewise_nstates[j]);
    vec_zero(f);

    if (type == UNIF) {
      for (mut = 1; mut < M->sitewise_nstates[j]; mut++)
      vec_set(f, mut, 1.0/(M->sitewise_nstates[j]-1.0));
    }

    else {
      for (i = 0; i < M->ncells; i++) {
        int mut = cpr_get_mut(M, i, j);
        if (mut > 0)
          vec_set(f, mut, vec_get(f, mut) + 1.0);
      }
      pv_normalize(f);
    }
    
    lst_push_ptr(sitewise_mutrates, f);
  }
  
  return sitewise_mutrates;
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

CrisprAncestralStateSets *cpr_new_state_sets(int nnodes) {
  CrisprAncestralStateSets *retval = malloc(sizeof(CrisprAncestralStateSets));
  retval->nnodes = nnodes;
  retval->NORESTRICT = -1; 
  retval->nodetypes = smalloc(nnodes * sizeof(int));
  retval->unr_lists = lst_new_ptr(100); /* starting size; will realloc */
  retval->restr_lists = lst_new_ptr(100); 
  retval->sil_lists = lst_new_ptr(100); 
  return retval;
}

void cpr_state_sets_resize(CrisprAncestralStateSets *sets, int newnnodes) {
  sets->nodetypes = srealloc(sets->nodetypes, newnnodes * sizeof(int));
  sets->nnodes = newnnodes;
}

/* efficiently retrieve ancestral state set for a given node type and
   total number of states, using caching.  Specifically, if nodetype
   is NORESTRICT, returns a list of integers from 0 to nstates - 1
   inclusive.  If nodetype has another (restricted) value, returns a
   list of two integers consisting of 0 and that value.  If nodetype
   is 0, returns a list containing zero only. */
List *cpr_get_state_set(CrisprAncestralStateSets *set, TreeNode *n, int nstates) {
  int i, j;

  int nodetype = set->nodetypes[n->id];
  
  if (nodetype == set->NORESTRICT) {

    /* special case: node is a leaf but has NORESTRICT type; can
       only arise when node is observed to have silent state, in
       which case that is the only state we need to consider */
    if (n->lchild == NULL) {
        if (nstates >= lst_size(set->sil_lists)) {
          /* fill in all of the ones up to the requested size and cache them */
          for (i = lst_size(set->sil_lists); i <= nstates; i++) {
            List *l = NULL;
            if (i > 1) {
              l = lst_new_int(1);
              lst_push_int(l, i-1);
            }
            lst_push_ptr(set->sil_lists, l); /* the ith element of
                                                sil_lists will be a list
                                                consisting of just i-1 (or
                                                NULL if i <= 1) */
          }
        }
        return lst_get_ptr(set->sil_lists, nstates);
    }

    /* otherwise return the appropriate unrestricted list */
    if (nstates >= lst_size(set->unr_lists)) {
      /* fill in all of the ones up to the requested size and cache them */
      for (i = lst_size(set->unr_lists); i <= nstates; i++) {
        List *l = lst_new_int(i);
        for (j = 0; j < i; j++)
          lst_push_int(l, j);
        lst_push_ptr(set->unr_lists, l); /* the ith element of
                                            unr_lists will be a list
                                            of integers from 0 to i,
                                            inclusive */
      }
    }
    return lst_get_ptr(set->unr_lists, nstates);
  }
  else { /* restricted list; includes the 0 case */
    if (nodetype >= lst_size(set->restr_lists)) {
      /* fill in all of the ones up to the requested size and cache them */
      for (i = lst_size(set->restr_lists); i <= nodetype; i++) {
        List *l = lst_new_int(2);
        lst_push_int(l, 0);
        if (i > 0) 
          lst_push_int(l, i);
        lst_push_ptr(set->restr_lists, l); /* the ith element of
                                              restr_lists will be a list
                                              consisting of 0 and i (or
                                              just 0 if i == 0) */
      }
    }
    return lst_get_ptr(set->restr_lists, nodetype);
  }
}

void cpr_free_state_sets(CrisprAncestralStateSets *sets) {
  int i;
  List *l;
  free(sets->nodetypes);
  for (i = 0; i < lst_size(sets->unr_lists); i++) {
    l = lst_get_ptr(sets->unr_lists, i);
    if (l == NULL) break; /* non-NULL must be contiguous */
    lst_free(l);
  }
  for (i = 0; i < lst_size(sets->restr_lists); i++) {
    l = lst_get_ptr(sets->restr_lists, i);
    if (l == NULL) break;
    lst_free(l);
  }
  for (i = 0; i < lst_size(sets->sil_lists); i++) {
    l = lst_get_ptr(sets->restr_lists, i);
    if (l != NULL) 
      lst_free(l);
  }
  free(sets);
}

/* renumber mutation states so they are dense for each site; needed
   for sitewise mutation matrices */
CrisprMutTable *cpr_new_sitewise_table(CrisprMutTable *origM) {
  int i, j, k, state, newstate;
  CrisprMutTable *M = cpr_copy_table(origM);
  int map[origM->nstates]; 

  M->sitewise_nstates = smalloc(origM->nsites * sizeof(int));
  
  for (j = 0; j < origM->nsites; j++) {
    /* create mapping for col j */
    for (k = 0; k < origM->nstates; k++) map[k] = -1;
    M->sitewise_nstates[j] = 1;
    
    for (i = 0; i < origM->ncells; i++) {
      state = cpr_get_mut(origM, i, j);
      assert(state < origM->nstates);  
      
      if (state == -1 || state == 0)
        newstate = state;
      else {
        if (map[state] == -1)
          map[state] = M->sitewise_nstates[j]++;
        
        newstate = map[state];
      }
      cpr_set_mut(M, i, j, newstate);
    }
  }
  return M;
}

/* create and return a new model based on a given mutation matrix and tree model */
CrisprMutModel *cpr_new_model(CrisprMutTable *M, TreeModel *mod,
                              enum crispr_model_type modtype,
                              enum crispr_mutrates_type mrtype) {
  CrisprMutModel *retval = smalloc(sizeof(CrisprMutModel));
  retval->model_type = modtype;
  retval->mod = mod;
  retval->mut = M;
  retval->nsites = M->nsites;
  retval->ncells = M->ncells;
  retval->nstates = M->nstates;
  retval->sil_rate = CPR_SIL_RATE_INIT;
  retval->leading_t = 0.05;
  retval->mutrates = NULL;
  retval->sitewise_mutrates = NULL;
  retval->Pt = NULL;
  retval->mutrates_type = mrtype;
  return retval;
}
  
/* preprocessing steps for likelihood calculations -- compute
   equilibrium frequencies, initialize substitution models, etc. This
   function allocates new memory and should be called only once prior
   to repeated likelihood calculations */
void cpr_prep_model(CrisprMutModel *cprmod) {
  int j, nodeidx;
  List *thisPt;
  
  if (cprmod->model_type == SITEWISE) {
    cprmod->mut = cpr_new_sitewise_table(cprmod->mut); /* replace with pointer to new table */
    
    /* compute equilibrium frequencies */
    cprmod->sitewise_mutrates = cpr_estim_sitewise_mutrates(cprmod->mut,
                                                            cprmod->mutrates_type);  
 
    /* allocate memory for sitewise, branchwise substitution matrices */
    cprmod->Pt = lst_new_ptr(cprmod->nsites);
    for (j = 0; j < cprmod->nsites; j++) {
      thisPt = lst_new_ptr(cprmod->mod->tree->nnodes);
      for (nodeidx = 0; nodeidx < cprmod->mod->tree->nnodes; nodeidx++)
        lst_push_ptr(thisPt, mm_new(cprmod->mut->sitewise_nstates[j]+1, NULL, DISCRETE));
      lst_push_ptr(cprmod->Pt, thisPt);
    }
  }
  else {
    /* no need to alter the mutation table in this case; just build
       global eq freqs */
    cprmod->mutrates = cpr_estim_mutrates(cprmod->mut,
                                          cprmod->mutrates_type);

    /* make all sitewise eqfreqs point to the global eqfreqs */
    cprmod->sitewise_mutrates = lst_new_ptr(cprmod->nsites);
    for (j = 0; j < cprmod->nsites; j++) 
      lst_push_ptr(cprmod->sitewise_mutrates, cprmod->mutrates);
    
    /* build one list of substitution matrices and make
       all sitewise models point to it */
    cprmod->Pt = lst_new_ptr(cprmod->nsites);
    thisPt = lst_new_ptr(cprmod->mod->tree->nnodes);
    for (nodeidx = 0; nodeidx < cprmod->mod->tree->nnodes; nodeidx++)
      lst_push_ptr(thisPt, mm_new(cprmod->nstates+1, NULL, DISCRETE));
    for (j = 0; j < cprmod->nsites; j++) 
      lst_push_ptr(cprmod->Pt, thisPt);
  }
}

/* dump a CrisprMutModel object to a file.  For debugging */
void cpr_print_model(CrisprMutModel *cprmod, FILE *F) {
  int j, nodeidx;
  fprintf(F, "CrisprMutModel:\nmodel_type = %s\nnsites = %d\nncells = %d\nnstates = %d\nsilencing_rate = %f\n",
          (cprmod->model_type == SITEWISE ? "SITEWISE" : "GLOBAL"), cprmod->nsites,
          cprmod->ncells, cprmod->nstates, cprmod->sil_rate);
  
  for (j = 0; j < cprmod->nsites; j++) {
    List *thisPt = lst_get_ptr(cprmod->Pt, j);
    Vector *mutrates = lst_get_ptr(cprmod->sitewise_mutrates, j);
    fprintf(F, "Model for site %d:\n", j);
    fprintf(F, "Mutation rates:\n");
    vec_print(mutrates, F);
    for (nodeidx = 0; nodeidx < cprmod->mod->tree->nnodes; nodeidx++) {
      TreeNode *n = lst_get_ptr(cprmod->mod->tree->nodes, nodeidx);
      MarkovMatrix *mm = lst_get_ptr(thisPt, nodeidx);
      fprintf(F, "Node %d (dparent %f):\nPt:\n", n->id, n->dparent);
      mm_pretty_print(F, mm);
    }
  }
}

/* free memory allocated by cpr_prep_model */
void cpr_free_model(CrisprMutModel *cprmod) {
  int j, nodeidx;
  List *l;
  if (cprmod->model_type == SITEWISE) {
    for (j = 0; j < cprmod->nsites; j++) {
      vec_free(lst_get_ptr(cprmod->sitewise_mutrates, j));
      l = lst_get_ptr(cprmod->Pt, j);
      for (nodeidx = 0; nodeidx < cprmod->mod->tree->nnodes; nodeidx++)
        mm_free(lst_get_ptr(l, nodeidx));
      lst_free(l);
    }
  }
  else {
    l = lst_get_ptr(cprmod->Pt, 0);
    for (nodeidx = 0; nodeidx < cprmod->mod->tree->nnodes; nodeidx++)
      mm_free(lst_get_ptr(l, nodeidx));
    lst_free(l);
  }
  if (cprmod->mutrates != NULL)
    vec_free(cprmod->mutrates);
  if (cprmod->sitewise_mutrates != NULL)
    lst_free(cprmod->sitewise_mutrates);
  if (cprmod->Pt != NULL)
    lst_free(cprmod->Pt);
}

/* update substitution matrices for a new set of branch lengths */
void cpr_update_model(CrisprMutModel *cprmod) {
  int j;
  if (cprmod->model_type == SITEWISE) {
    for (j = 0; j < cprmod->nsites; j++) 
      cpr_set_subst_matrices(cprmod->mod, cprmod->sil_rate,
                             lst_get_ptr(cprmod->Pt, j),
                             lst_get_ptr(cprmod->sitewise_mutrates, j));                                              
  }
  else 
    cpr_set_subst_matrices(cprmod->mod, cprmod->sil_rate,
                           lst_get_ptr(cprmod->Pt, 0), cprmod->mutrates);
}
