/* grid for d-dimensional probability simplex: partitions into regions
   that intersect hypercubes in d-dimensional space */

#include <simplex_grid.h>
#include <misc.h>
#include <external_libs.h>

/* return single index from coords in d dimensions */
static PHAST_INLINE
int int_idx(int *coord, int dim, int nrows) {
  int i, retval = 0;
  for (i = 0; i < dim; i++) {
    retval += coord[i];
    if (i < dim-1) retval *= nrows;
  }
  return retval;
}

/* inverse of above: obtain d-dimensional coords from single index */
static PHAST_INLINE
void get_coords(int *coord, int dim, int nrows, int idx_d) {
  int i;
  for (i = dim-1; i >= 0; i--) {
    coord[i] = idx_d % nrows;
    if (i > 0) idx_d /= nrows;    
  }
}

SimplexRegion *sxg_new_region(int dim, int nrows, int *coord) {
  SimplexRegion *sr = smalloc(sizeof(SimplexRegion));
  int i, alpha = nrows;

  sr->d = dim;
  sr->lb = smalloc(dim * sizeof(double));
  sr->ub = smalloc(dim * sizeof(double));
  sr->centroid = vec_new(dim);

  for (i = 0; i < dim; i++) alpha -= coord[i];

  if (alpha < 1 || alpha > dim - 1)
    die("ERROR: region does not intersect with simplex.\n");

  for (i = 0; i < dim; i++) {
    sr->lb[i] = (double) coord[i] / nrows;
    sr->ub[i] = ((double) coord[i] + 1) / nrows;
    vec_set(sr->centroid, i, sr->lb[i] + ((double) alpha) / (dim * nrows));
  }

  return sr;
}

void sxg_free_region(SimplexRegion *sr) {
  sfree(sr->lb);
  sfree(sr->ub);
  vec_free(sr->centroid);
  sfree(sr);
}

/* create and return a SimplexGrid for a simplex of the given
   dimension; the grid will be defined by nrows rows per
   dimension */
SimplexGrid *sxg_build_grid(int dim, int nrows) {
  int i, j, idx, alpha;
  int coord[256];
  int maxsize = int_pow(nrows, dim);
  SimplexGrid *g = smalloc(sizeof(SimplexGrid));

  g->nregs = sxg_nregions(dim, nrows);
  
  if (dim >= 256) die("ERROR sxg_build_grid: dim must be < 256, but is %i\n",
		      dim);
  g->d = dim;
  g->nrows = nrows;
  g->sr_d = smalloc(maxsize * sizeof(void*));
  g->sr = smalloc(g->nregs * sizeof(void*));

  idx = 0;			
  for (i = 0; i < maxsize; i++) {
    get_coords(coord, dim, nrows, i);
    alpha = nrows;
    for (j = 0; j < dim; j++) alpha -= coord[j];
    if (alpha >= 1 && alpha <= dim - 1) {
      g->sr_d[i] = sxg_new_region(dim, nrows, coord);
      g->sr_d[i]->idx_d = i;
      g->sr[idx] = g->sr_d[i];
      g->sr_d[i]->idx = idx;
      idx++;
    }
    else 
      g->sr_d[i] = NULL;
  }
  if (idx != g->nregs)
    die("ERROR sxg_build_grid: idx (%i) should equal g->nregs (%i)\n",
	idx, g->nregs);

  return g;
}

void sxg_free_grid(SimplexGrid *g) {
  int i;
  for (i = 0; i < g->nregs; i++)
    sxg_free_region(g->sr[i]);
  sfree(g->sr_d);
  sfree(g->sr);
  sfree(g);
}

/* return SimplexRegion containing given point */
SimplexRegion *sxg_get_region(SimplexGrid *g, Vector *p) {
  int i, alpha = g->nrows;
  int coord[256];
  if (g->d >= 256)
    die("ERROR sxg_get_region: g->d should be < 256 but is %i\n", g->d);
  for (i = 0; i < g->d; i++) {
    if (vec_get(p, i) == 1) coord[i] = g->nrows - 1;
    else coord[i] = vec_get(p, i) * g->nrows;
    alpha -= coord[i];
  }
  if (alpha < 1 || alpha > g->d-1)
    die("ERROR sxg_get_region: alpha should be between 1 and %i, but is %i\n",
	g->d-1, alpha);
  return g->sr_d[int_idx(coord, g->d, g->nrows)];
}

/* number of regions in a simplex grid having dimension dim and number
   of rows equal to nrows */
unsigned sxg_nregions(int dim, int nrows) {
  return combinations(nrows + dim - 1, dim) - combinations(nrows, dim);
}

/* return maximum number of rows a simplex grid of dimension d can have
   if the number of regions can be at most nregs */
unsigned sxg_max_nrows(int dim, int nregs) {
  /* The number of regions is upper bounded by (nrows + dim - 1)^dim /
     dim!.  Therefore, we know nrows can be at least floor((nregs *
     dim!)^(1/dim) - dim + 1)  */
  int nrows = floor(pow(nregs * permutations(dim), 1.0/dim) - dim + 1);

  while (sxg_nregions(dim, nrows+1) < nregs)
    nrows++;

  return nrows;
}

