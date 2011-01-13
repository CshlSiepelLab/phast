/***************************************************************************
 * PHAST: PHylogenetic Analysis with Space/Time models
 * Copyright (c) 2002-2005 University of California, 2006-2010 Cornell 
 * University.  All rights reserved.
 *
 * This source code is distributed under a BSD-style license.  See the
 * file LICENSE.txt for details.
 ***************************************************************************/

/** @file simplex_grid.h
   Grid for d-dimensional probability simplex: partitions into regions
   that intersect hypercubes in d-dimensional space 
   @ingroup prequel
*/

#ifndef SIMPLEX_GRID
#define SIMPLEX_GRID

#include <vector.h>

/** Region within simplex -- intersection of simplex with hypercube in
    d-dimensional space. */
typedef struct {
  unsigned d;			/**< Dimension */
  unsigned idx;			/**< Index of region, in [0, nregs) */
  unsigned idx_d;		/**< D-dimensional index of region, in
				   [0, nrows^d) */
  double *lb;			/**< Lower bounds (d dimensions) */
  double *ub;			/**< Upper bounds (d dimensions) */
  Vector *centroid;		/**< Center of the region in d dimensions*/
} SimplexRegion;

/** Simplex Grid made up of multiple Simplex Regions*/
typedef struct {
  unsigned d;			/**< Dimension */
  unsigned nrows;		/**< Number of rows per dimension */
  unsigned nregs;		/**< Number of regions in grid */
  SimplexRegion **sr;		/**< Pointers to regions, by index */
  SimplexRegion **sr_d;		/**< Pointers to regions, by
				   d-dimensional index */
} SimplexGrid;

/** Create a new Simplex region.
  @param dim Dimension the simplex region should be created in
  @param nrows Number of rows within the simplex region
  @param coord Coordinates for new region in dim dimensions represented as Array of size dim
*/
SimplexRegion *sxg_new_region(int dim, int nrows, int *coord);

/** Free a simplex region object
  @param sr Simplex region to free
*/
void sxg_free_region(SimplexRegion *sr);

/** Create a new Simplex grid.
  @param dim Dimension of simplex grid
  @param nrows Number of rows in simplex grid  
 */
SimplexGrid *sxg_build_grid(int dim, int nrows);

/** Free Simplex Grid object.
   @param g Simplex Grid object to free
*/
void sxg_free_grid(SimplexGrid *g);

/** Return a specified Simplex region containing given point from Simplex Grid.
  @param g Simplex Grid containing region to retrieve
  @param p Coordinates of simplex region to retrieve (may be multi dimensional)
  @result Simplex Region 
*/
SimplexRegion *sxg_get_region(SimplexGrid *g, Vector *p);

/** Number of regions in a simplex grid having dimension dim and number
   of rows equal to nrows
   @param dim Only count regions with dimension == dim
   @param nrows Only count regions with number of rows == nrows 
   @result Number of regions that matched query
*/
unsigned sxg_nregions(int dim, int nrows);

/** Return maximum number of rows a simplex grid of dimension d can have
   if the number of regions can be at most nregs 
   @param dim Dimension of simplex grid
   @param nregs Maximum amount of regions for which to get max number of rows
   @result Max number of rows in grid that matched query
*/
unsigned sxg_max_nrows(int dim, int nregs);

#endif






