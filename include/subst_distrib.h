#ifndef SUBST_DISTRIB
#define SUBST_DISTRIB

#include <vector.h>
#include <msa.h>
#include <tree_model.h>

typedef struct {
  int njumps_max;
  double lambda;
  Matrix *R;
  TreeModel *mod;
  Matrix **A;
  Matrix ***B;
} JumpProcess;


JumpProcess *sub_define_jump_process(TreeModel *mod, int njumps_max);
void sub_free_jump_process(JumpProcess *jp);
Vector *sub_distrib_branch(JumpProcess *jp, double t);
Matrix **sub_distrib_branch_conditional(JumpProcess *jp, double t);
Vector *sub_prior_distrib_site(JumpProcess *jp);
Vector *sub_posterior_distrib_site(JumpProcess *jp, MSA *msa, int tuple_idx);
Vector *sub_prior_distrib_alignment(JumpProcess *jp, int nsites);
Vector *sub_posterior_distrib_alignment(JumpProcess *jp, MSA *msa);
void sub_posterior_stats_alignment(JumpProcess *jp, MSA *msa, 
                                   double *mean, double *variance);
Matrix *sub_joint_distrib_site(JumpProcess *jp, MSA *msa, int tuple_idx);
Matrix *sub_prior_joint_distrib_alignment(JumpProcess *jp, int nsites);
Matrix *sub_posterior_joint_distrib_alignment(JumpProcess *jp, MSA *msa);
void sub_posterior_joint_stats_alignment(JumpProcess *jp, MSA *msa, 
                                         double *mean_tot, double *var_tot,
                                         double *mean_left, double *var_left,
                                         double *mean_right, double *var_right);
#endif
