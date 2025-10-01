#ifndef VARRES_H
#define VARRES_H

#include <stdio.h>
#include <limits.h>
#include <phast/tree_model.h>
#include <phast/mvn.h>
#include <phast/multi_mvn.h>
#include <phast/nj.h>

List *nj_importance_sample(int nsamples, List *trees, Vector *logdens,
                           TreeModel *mod, CovarData *data, FILE *logf);

List *nj_var_sample_rejection(int nsamples, multi_MVN *mmvn,
                              CovarData *data, TreeModel *mod,
                              FILE *logf);

List *nj_var_sample_importance(int nsamples, multi_MVN *mmvn,
                               CovarData *data, TreeModel *mod,
                               FILE *logf);

List *nj_var_sample_importance_psis(int nsamples, multi_MVN *mmvn,
                                    CovarData *data, TreeModel *mod,
                                    FILE *logf);

List *nj_var_sample_importance_moc(int nsamples, multi_MVN *mmvn,
                                   CovarData *data, TreeModel *mod,
                                   FILE *logf);

List *nj_var_sample_smc(int nsamples, multi_MVN *mmvn,
                        CovarData *data, TreeModel *mod,
                        FILE *logf);

List *nj_var_sample_pcn(int nsamples, multi_MVN *mmvn,
                        CovarData *data, TreeModel *mod,
                        FILE *logf /* unused; we print to stderr */);

List *nj_var_sample_pcn_blocked(int nsamples, multi_MVN *mmvn,
                               CovarData *data, TreeModel *mod,
                               FILE *logf /* unused; diagnostics -> stderr */);

List *nj_var_sample_pcn_pt(int nsamples, multi_MVN *mmvn,
                           CovarData *data, TreeModel *mod,
                           FILE *logf /* unused; diagnostics -> stderr */);
#endif
