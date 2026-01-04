/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/*
  
  This file contains some definitions which are useful for linking to
  external libraries (R, clapack, veclib)

 */

#ifndef _PHAST_EXTERNAL_LIBS_
#define _PHAST_EXTERNAL_LIBS_

#include <stdlib.h>

#ifdef RPHAST
#include <Rconfig.h>
#define PHAST_INLINE R_INLINE
#else
#define PHAST_INLINE inline
#endif

/* ---------------- LAPACK backend selection ----------------
 *
 * Exactly one of the following should be defined:
 *   - R_LAPACK
 *   - VECLIB
 *   - PHAST_USE_SYSTEM_LAPACK
 *   - PHAST_USE_CLAPACK   (legacy)
 */

#if defined(R_LAPACK)

  #include <R_ext/Lapack.h>
  typedef int    LAPACK_INT;
  typedef double LAPACK_DOUBLE;

#elif defined(VECLIB)

  #include <Accelerate/Accelerate.h>

  #if defined(ACCELERATE_NEW_LAPACK)
    /* New Accelerate LAPACK interface: use plain C types */
    typedef int    LAPACK_INT;
    typedef double LAPACK_DOUBLE;
  #else
    /* Old Accelerate CLAPACK-style API */
    typedef __CLPK_integer    LAPACK_INT;
    typedef __CLPK_doublereal LAPACK_DOUBLE;
  #endif

#elif defined(PHAST_USE_SYSTEM_LAPACK)

  /* System LAPACK via Fortran ABI (OpenBLAS/netlib/MKL).
     We call dgesv_(), dgeev_(), etc. directly; no clapack.h needed. */
  typedef int    LAPACK_INT;
  typedef double LAPACK_DOUBLE;

#elif defined(PHAST_USE_CLAPACK)

  /* Legacy CLAPACK + f2c (only if explicitly enabled) */
  #include <f2c.h>
  #include <clapack.h>
  typedef integer    LAPACK_INT;
  typedef doublereal LAPACK_DOUBLE;

#else

  #error "No LAPACK backend selected. Define one of: R_LAPACK, VECLIB, PHAST_USE_SYSTEM_LAPACK (or PHAST_USE_CLAPACK)."

#endif

/* ------------------------------------------------------------------ */
/* Fortran LAPACK entry points (LP64).                                  */
/* We only declare these when using system LAPACK via Fortran ABI.      */
/* R_LAPACK provides its own declarations via R_ext/Lapack.h, and       */
/* Accelerate provides declarations via Accelerate headers.             */
/* ------------------------------------------------------------------ */

#if defined(PHAST_USE_SYSTEM_LAPACK) && !defined(R_LAPACK) && !defined(VECLIB)

void dgebal_(const char *job, const LAPACK_INT *n,
             LAPACK_DOUBLE *a, const LAPACK_INT *lda,
             LAPACK_INT *ilo, LAPACK_INT *ihi,
             LAPACK_DOUBLE *scale, LAPACK_INT *info);

void dgebak_(const char *job, const char *side, const LAPACK_INT *n,
             const LAPACK_INT *ilo, const LAPACK_INT *ihi,
             const LAPACK_DOUBLE *scale, const LAPACK_INT *m,
             LAPACK_DOUBLE *v, const LAPACK_INT *ldv, LAPACK_INT *info);

void dgesv_(const LAPACK_INT *n, const LAPACK_INT *nrhs,
            LAPACK_DOUBLE *a, const LAPACK_INT *lda,
            LAPACK_INT *ipiv,
            LAPACK_DOUBLE *b, const LAPACK_INT *ldb,
            LAPACK_INT *info);

void dgetrf_(const LAPACK_INT *m, const LAPACK_INT *n,
             LAPACK_DOUBLE *a, const LAPACK_INT *lda,
             LAPACK_INT *ipiv, LAPACK_INT *info);

void dgetri_(const LAPACK_INT *n,
             LAPACK_DOUBLE *a, const LAPACK_INT *lda,
             const LAPACK_INT *ipiv,
             LAPACK_DOUBLE *work, const LAPACK_INT *lwork,
             LAPACK_INT *info);

void dpotrf_(const char *uplo, const LAPACK_INT *n,
             LAPACK_DOUBLE *a, const LAPACK_INT *lda,
             LAPACK_INT *info);

#endif /* PHAST_USE_SYSTEM_LAPACK */

#endif /* _PHAST_EXTERNAL_LIBS_ */

