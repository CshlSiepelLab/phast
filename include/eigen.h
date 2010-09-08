/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/* $Id: eigen.h,v 1.2 2008-11-12 02:07:59 acs Exp $ */

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
