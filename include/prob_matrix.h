#ifndef PROB_MATRIX
#define PROB_MATRIX

#include <vector.h>
#include <matrix.h>

#define PM_EPS 1e-10

void pm_mean(Matrix *p, double *mean_x, double *mean_y);
Vector *pm_marg_x(Matrix *p);
Vector *pm_marg_y(Matrix *p);
Vector *pm_marg_tot(Matrix *p);
Vector *pm_x_given_tot(Matrix *p, int tot);
Vector *pm_y_given_tot(Matrix *p, int tot);
void pm_normalize(Matrix *p);
Matrix *pm_convolve(Matrix *p, int n);
Matrix *pm_convolve_many(Matrix **p, int *counts, int n);

#endif
