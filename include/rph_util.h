/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2009 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#ifndef RPHASTH
#define RPHASTH

extern int rphast_errno;
extern char rphast_errmsg[1000];

void* ad2ptr(double address);
double ptr2ad(void* ptr);

#endif
