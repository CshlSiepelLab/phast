/* $Id: eigen.h,v 1.1 2005-06-22 07:11:20 acs Exp $
   Written by Adam Siepel, 2002-2005
   Copyright 2002-2005, Adam Siepel, University of California 
*/

/* \file eigen.h
   Eigenvalue computation and matrix diagonalization.  CLAPACK is used.
   \ingroup base
*/

#ifndef EIG_H
#define EIG_H

#include <matrix.h>
#include <complex_vector.h>
#include <complex_matrix.h>

int mat_diagonalize(Matrix *M, Zvector *eval, Zmatrix *revect, Zmatrix *levect);
int mat_eigenvals(Matrix *M, Zvector *evals);

#endif
