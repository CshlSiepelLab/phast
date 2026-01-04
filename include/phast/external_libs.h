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

#endif /* _PHAST_EXTERNAL_LIBS_ */

