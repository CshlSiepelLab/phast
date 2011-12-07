/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "phast.help"

int main(int argc, char *argv[]) {
  printf("PHAST %s\n", PHAST_VERSION);
  printf("%s", HELP);
  return 0;
}
