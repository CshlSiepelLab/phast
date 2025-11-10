/* simple representations of multiDAGs for cell state migration histories.
 * Allows timing information to be associated with edges.  Designed to support
 * many variations on the same graph */

#ifndef MULTIDAG_H
#define MULTIDAG_H

#include <phast/stringsplus.h>
#include <phast/hashtable.h>
#include <phast/crispr.h>
#include <phast/migration.h>
#include <phast/lists.h>

typedef struct {
  int from_state;
  int to_state;
  double start_time; 
  double end_time; 
} Edge;

typedef struct mdag {
  int id; /* unique identifier for graph */
  List *edges; /* list of Edge* */
  int nedges;
  MigTable *migtable; /* associated migration table */
  int nstates; /* for convenience; matches MigTable */   
} MultiDAG;

/* set of related multiDAGs, typically sampled from posterior */
typedef struct {
  List *dags; /* list of MultiDAG* */
  int ndags;
} MultiDAGSet;

MultiDAG *mdag_new(MigTable *mg);
void mdag_free(MultiDAG *G);
Edge *mdag_add_edge(MultiDAG *G, int from_state, int to_state,
                    double start_time, double end_time);
void mdag_add_to_set(MultiDAGSet *S, MultiDAG *G);
void mdag_set_free(MultiDAGSet *S);
void mdag_print_dot(MultiDAG *G, FILE *F);
void mdg_set_print_dot(MultiDAGSet *S, FILE *F);


#endif /* MULTIDAG_H */
