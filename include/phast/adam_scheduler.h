/* simple scheduler for Adam algorithm used in VINE */

#ifndef ASCHED
#define ASCHED

#include <stdio.h>

/* config & state */
typedef struct {
  int    N_sites;           /* total alignment length */
  int    m0;                /* initial sites per minibatch */
  int    inc_every;         /* grow minibatch every 'inc_every' steps */
  double lr_full;           /* target learning rate at end (full gradients) */
  double lr_alpha;          /* exponent for lr ramp-down (e.g., 1.0
                               linear, 0.5 sqrt) */
  double tau_full;          /* fraction of data for "full" gradients
                               (e.g., 0.95) */
  double clip_max_norm;     /* 0 disables clipping; else L2 cap (e.g., 5.0) */
  int    adaptive_clip;     /* 1: clip threshold tracks EMA of grad norm */
  int    clip_warmup;       /* steps before adaptive clip kicks in */
  double clip_factor;       /* factor for adaptive clipping (e.g., 2.0) */
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
} SchedState;

/* outputs per step */
typedef struct {
  double lr;            /* learning rate to use *this* step */
  int    m;             /* how many sites to sample *this* step */
  double clip_norm;     /* clip threshold to apply (0 => no clip) */
  int resample_sites;   /* 1 => pick a fresh subset; 0 => reuse previous */
  int full_grad_now;    /* 1 => compute full gradient this step */
} SchedDirectives;

/* simple telemetry fed back each step */
typedef struct {
  double grad_norm;     /* L2 norm of the *raw* minibatch gradient (pre-clip) */
} SchedMetrics;

/* API */
Scheduler *sched_new(int N_sites, int init_subsample, int inc_every,
                     double target_lr, int persist_k, int fullgrad_every,
                     int clip_warmup);
SchedState* sched_new_state(const Scheduler *cfg);
void sched_next(const Scheduler *cfg, SchedState *st,
                const SchedMetrics *metrics, SchedDirectives *out);

#endif
