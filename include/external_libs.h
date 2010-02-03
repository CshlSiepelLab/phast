/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2009 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/*
  
  This file contains some definitions which are useful for linking to
  external libraries (R, clapack, veclib)

 */

#ifndef PHAST_MISC_DEFS
#define PHAST_MISC_DEFS

#include <stdlib.h>

#ifdef RPHAST
//#include <R_ext/Error.h>
#include <Rconfig.h>
#define PHAST_INLINE R_INLINE
#else
#define PHAST_INLINE inline
#endif

#ifdef R_LAPACK
#include <R_ext/Lapack.h>
#define LAPACK_INT int
#define LAPACK_DOUBLE double
#else

#ifdef VECLIB
#include <vecLib/clapack.h>
#define LAPACK_INT __CLPK_integer
#define LAPACK_DOUBLE __CLPK_doublereal
#else

#ifndef SKIP_LAPACK 
#include <f2c.h>  
#include <clapack.h> 
#define LAPACK_INT integer
#define LAPACK_DOUBLE doublereal
#endif  /*ifndef SKIP_LAPACK */

#endif  /*ifdef VECLIB */

#endif  /* ifdef R_LAPACK */

#endif  /* ifndef PHAST_MISC_DEFS */
