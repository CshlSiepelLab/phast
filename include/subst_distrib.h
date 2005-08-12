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

typedef struct {
  double prior_mean, prior_var, post_mean, post_var, p_cons, p_anti_cons;
  int prior_min, prior_max, post_min, post_max;
} p_value_stats;

typedef struct {
  double prior_mean_sup, prior_mean_sub, prior_var_sup, prior_var_sub, 
    post_mean_sup, post_mean_sub, post_mean_tot, post_var_sup, 
    post_var_sub, post_var_tot, p_cons_sub, p_cons_sup, p_anti_cons_sup, 
    p_anti_cons_sub, cond_p_cons_sub, cond_p_anti_cons_sub;
  int prior_min_sup, prior_min_sub, prior_max_sup, prior_max_sub, 
    post_min_sup, post_min_sub, post_min_tot, post_max_sup, post_max_sub,
    post_max_tot;
} p_value_joint_stats;

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
p_value_stats *sub_p_value_many(JumpProcess *jp, MSA *msa, GFF_Set *feats, 
                                double ci);
p_value_joint_stats* sub_p_value_joint_many(JumpProcess *jp, MSA *msa, 
                                            GFF_Set *feats, double ci);
#endif
