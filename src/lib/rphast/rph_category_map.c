/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2009 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/*****************************************************
rph_category_map.c
The RPHAST handles to functions dealing with
category maps in the PHAST package.

Melissa Hubisz
Last updated: 1/5/2010
*****************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <msa.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>
#include <ctype.h>
#include <sufficient_stats.h>
#include <local_alignment.h>
#include <category_map.h>

#include <Rdefines.h>


SEXP rph_cm_new_from_gff(SEXP gff) {
  return R_MakeExternalPtr(cm_new_from_features((GFF_Set*)EXTPTR_PTR(gff)),
			   R_NilValue, R_NilValue);
}
