/* grid for d-dimensional probability simplex: partitions into regions
   that intersect hypercubes in d-dimensional space */

#ifndef SIMPLEX_GRID
#define SIMPLEX_GRID

#include <vector.h>

/** Region within simplex -- intersection of simplex with hypercube in
    d-dimensional space */
typedef struct {
  unsigned d;			/* dimension */
  unsigned idx;			/* index of region, in [0, nregs) */
  unsigned idx_d;		/* d-dimensional index of region, in
				   [0, nrows^d) */
  double *lb;			/* lower bounds (d dimensions) */
  double *ub;			/* upper bounds (d dimensions) */
  Vector *centroid;
} SimplexRegion;

typedef struct {
  unsigned d;			/* dimension */
  unsigned nrows;		/* number of rows per dimension */
  unsigned nregs;		/* number of regions in grid */
  SimplexRegion **sr;		/* pointers to regions, by index */
  SimplexRegion **sr_d;		/* pointers to regions, by
				   d-dimensional index */
} SimplexGrid;

SimplexRegion *sxg_new_region(int dim, int nrows, int *coord);
void sxg_free_region(SimplexRegion *sr);
SimplexGrid *sxg_build_grid(int dim, int nrows);
void sxg_free_grid(SimplexGrid *g);
SimplexRegion *sxg_get_region(SimplexGrid *g, Vector *p);
unsigned sxg_nregions(int dim, int nrows);
unsigned sxg_max_nrows(int dim, int nregs);

#endif
