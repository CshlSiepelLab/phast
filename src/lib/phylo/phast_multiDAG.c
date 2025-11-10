/* simple representations of multiDAGs for cell state migration histories.
 * Allows timing information to be associated with edges.  Designed to support
 * many variations on the same graph */

#include <stdio.h>
#include <stdlib.h>
#include <phast/stringsplus.h>
#include <phast/hashtable.h>
#include <phast/crispr.h>
#include <phast/migration.h>
#include <phast/multiDAG.h>

static int idcounter = 0;

MultiDAG *mdag_new(MigTable *mg) {
  MultiDAG *G = smalloc(sizeof(MultiDAG));
  G->edges = lst_new_ptr(100);
  G->nedges = 0;
  G->migtable = mg;
  G->nstates = 0;
  G->id = idcounter++;
  return G;
}

void mdag_free(MultiDAG *G) {
  for (int i = 0; i < lst_size(G->edges); i++) {
    Edge *e = lst_get_ptr(G->edges, i);
    sfree(e);
  }
  lst_free(G->edges);
  sfree(G);
}

Edge *mdag_add_edge(MultiDAG *G, int from_state, int to_state,
                    double start_time, double end_time) {
  Edge *e = smalloc(sizeof(Edge));
  e->from_state = from_state;
  e->to_state = to_state;
  e->start_time = start_time;
  e->end_time = end_time;
  lst_push_ptr(G->edges, e);
  G->nedges++;
  return e;
}

void mdag_add_to_set(MultiDAGSet *S, MultiDAG *G) {
  lst_push_ptr(S->dags, G);
  S->ndags++;
}

void mdag_set_free(MultiDAGSet *S) {
  for (int i = 0; i < lst_size(S->dags); i++) {
    MultiDAG *G = lst_get_ptr(S->dags, i);
    mdag_free(G);
  }
  lst_free(S->dags);
  sfree(S);
}

/* output a DAG as a single line in dot format. Floating point label
   is midpoint of start and end times.  The labels for the nodes (such
   as A, B, C) are obtained from the associated MigTable */
void mdag_print_dot(MultiDAG *S, FILE *F) {
  fprintf(F, "digraph G%d { rankdir=LR;\n", S->id);
  for (int i = 0; i < lst_size(S->edges); i++) {
    Edge *e = lst_get_ptr(S->edges, i);
    String *from_name = lst_get_ptr(S->migtable->statenames, e->from_state);
    String *to_name = lst_get_ptr(S->migtable->statenames, e->to_state);
    double mid_time = (e->start_time + e->end_time) / 2.0;
    fprintf(F, "  %s -> %s [label=\"%.4f\"];\n",
            from_name->chars, to_name->chars, mid_time);
  }
  fprintf(F, "}\n");
}

void mdag_print_set_dot(MultiDAGSet *S, FILE *F) {
  for (int i = 0; i < lst_size(S->dags); i++) {
    MultiDAG *G = lst_get_ptr(S->dags, i);
    mdag_print_dot(G, F);
  }
}
