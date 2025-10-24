/* simple scheduler for Adam algorithm used in VINE */

#ifndef ASCHED
#define ASCHED

#include <stdio.h>

/* config & state */
typedef struct {
  int    N_sites;           /* total alignment length */
  int    m0;                /* initial sites per minibatch (e.g., 256) */
  int    inc_every;         /* grow minibatch every 'inc_every' steps */
  double lr0;               /* base learning rate */
  double noise_const;       /* keep lr * m ~= noise_const (set to lr0*m0) */
  int    keep_noise_const;  /* 1: keep lr*m constant; 0: independent lr schedule */
  double clip_max_norm;     /* 0 disables clipping; else L2 cap (e.g., 5.0) */
  int    adaptive_clip;     /* 1: clip threshold tracks EMA of grad norm */
  double clip_beta;         /* EMA factor (e.g., 0.95) */
  int    persist_k;         /* reuse same subset of sites for k updates (e.g., 5) */
  int    fullgrad_every;    /* 0 disables; else compute full grad every K steps */
  int    T_total;           /* total intended steps (for optional cosine ramps) */
} Scheduler;

typedef struct {
  int    t;                 /* step counter */
  int    m_t;               /* current minibatch size */
  int    persist_left;      /* steps left before resampling sites */
  double ema_gnorm;         /* EMA of gradient norm */
  double lr_t;              /* current learning rate */
  double lambda_full;       /* blend weight for full grad [0..1] */
} SchedState;

/* outputs per step */
typedef struct {
  double lr;            /* learning rate to use *this* step */
  int    m;             /* how many sites to sample *this* step */
  double clip_norm;     /* clip threshold to apply (0 => no clip) */
  int    resample_sites;/* 1 => pick a fresh subset; 0 => reuse previous */
  double lambda_full;   /* optional: blend weight if anchor to full grad */
} SchedDirectives;

/* simple telemetry fed back each step */
typedef struct {
  double grad_norm;     /* L2 norm of the *raw* minibatch gradient (pre-clip) */
} SchedMetrics;

/* API */
Scheduler* sched_new(int N_sites, int init_subsample, int inc_every,
                     double init_lr, int persist_k, int fullgrad_every);
SchedState* sched_new_state(const Scheduler *cfg);
void sched_next(const Scheduler *cfg, SchedState *st,
                const SchedMetrics *metrics, SchedDirectives *out);
double sched_clip_scale(double *g, int n, double clip_norm); 

#endif
