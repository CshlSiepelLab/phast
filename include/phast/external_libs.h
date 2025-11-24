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

#ifdef R_LAPACK
   #include <R_ext/Lapack.h>
   #define LAPACK_INT int
   #define LAPACK_DOUBLE double

#elif defined(VECLIB)
   #include <Accelerate/Accelerate.h>
   #ifdef ACCELERATE_NEW_LAPACK
       /* New Accelerate LAPACK interface: use plain C types */
       typedef int    LAPACK_INT;
       typedef double LAPACK_DOUBLE;
   #else
     /* Old CLAPACK-style API (kept for older macOS) */
     #define LAPACK_INT    __CLPK_integer
     #define LAPACK_DOUBLE __CLPK_doublereal

   #endif  /* ACCELERATE_NEW_LAPACK */

#elif !defined(SKIP_LAPACK) && defined(PHAST_LAPACK_GENERIC)
   /* Generic LAPACK: use C types and Fortran symbols (dgesv_, dsyev_, etc.) */
   #define LAPACK_INT    int
   #define LAPACK_DOUBLE double

#elif !defined(SKIP_LAPACK)
    /* Legacy CLAPACK + f2c path */
    #include <f2c.h>
    #include <clapack.h>
    #define LAPACK_INT    integer
    #define LAPACK_DOUBLE doublereal

#else
   /* SKIP_LAPACK defined: LAPACK is disabled */

#endif /* R_LAPACK, VECLIB, SKIP_LAPACK, etc. */

#endif  /* ifndef PHAST_EXTERNAL_LIBS */
