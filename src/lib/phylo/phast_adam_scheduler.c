#include <math.h>
#include <string.h>
#include <phast/misc.h>
#include <phast/adam_scheduler.h>

/* set up a scheduler with appropriate defaults */
Scheduler* sched_new(int N_sites, int init_subsample, int inc_every,
                     double target_lr, int persist_k, int fullgrad_every) {
  Scheduler *s = smalloc(sizeof(Scheduler));
  s->N_sites = N_sites;
  s->m0 = init_subsample;
  s->inc_every = inc_every;
  s->lr_full = target_lr;
  s->lr_alpha = 0.5; /* sqrt decay */
  s->tau_full = 0.95; 
  s->clip_max_norm = 0; /* CHECK.  Should we use this? */
  s->adaptive_clip = FALSE; /* CHECK */
  s->persist_k = persist_k;
  s->fullgrad_every = fullgrad_every;
  s->T_total = 1000; /* CHECK */
  return s;
}

/* init state with defaults derived from cfg */
SchedState* sched_new_state(const Scheduler *cfg) {
  SchedState *st = smalloc(sizeof(SchedState));
  st->t = 0;
  st->m_t = (cfg->m0 > 0 ? cfg->m0 : (cfg->N_sites < 256 ? cfg->N_sites : 256));
  st->persist_left = cfg->persist_k > 0 ? cfg->persist_k : 1;
  st->ema_gnorm = 0.0;
  st->lr_t = cfg->lr_full; /* will be overwritten on first step */
  return st;
}

/* One scheduler “tick”. Call at the *start* of step t (before
 * computing next gradient), passing the metrics from step t-1 (or
 * NULL on t=0). */
void sched_next(const Scheduler *cfg, SchedState *st,
                const SchedMetrics *metrics, SchedDirectives *out) {
  /* 1) Update EMAs / adaptive stuff based on *previous* step’s gradient */
  if (metrics && metrics->grad_norm > 0) {
    if (st->ema_gnorm == 0.0) st->ema_gnorm = metrics->grad_norm;
    else st->ema_gnorm = cfg->clip_beta * st->ema_gnorm +
                         (1.0 - cfg->clip_beta) * metrics->grad_norm;
  }

  /* 2) Progressive batch growth (geometric-ish, smooth in practice) */
  int grow = (cfg->inc_every > 0 && (st->t > 0) && (st->t % cfg->inc_every == 0));
  if (grow) {
    int new_m = st->m_t * 2;                          /* double */
    if (new_m > cfg->N_sites) new_m = st->m_t + (cfg->N_sites - st->m_t)/2; /* taper */
    if (new_m < st->m_t + 1) new_m = st->m_t + 1;
    if (new_m > cfg->N_sites) new_m = cfg->N_sites;
    st->m_t = new_m;
  }

  /* 2.5) Check full-batch mode */
  double tau_full = (cfg->tau_full > 0.0 ? cfg->tau_full : 0.95);
  double frac = (double)st->m_t / (double)cfg->N_sites;
  if (frac < 1e-12) frac = 1e-12;
  int in_full_mode = (frac >= tau_full) || (st->m_t >= cfg->N_sites);
  if (in_full_mode) {
    st->m_t = cfg->N_sites;    /* pin to full */
    frac    = 1.0;
  }

  /* 3) Learning-rate policy  */
  double alpha = (cfg->lr_alpha > 0 ? cfg->lr_alpha : 0.5); // 0.5 good for Adam
  if (in_full_mode) {
    st->lr_t = cfg->lr_full;   /* exact tuned LR in full mode */
  } else {
    st->lr_t = cfg->lr_full * pow(frac, alpha);
    /* gentle global decay (keeps end behavior intact) */
    if (cfg->T_total > 0) {
      double mult = 0.9 + 0.1 * (0.5 * (1.0 + cos(M_PI * (double)st->t / (double)cfg->T_total)));
      st->lr_t *= mult;
    }
  }

  /* 4) Site subset persistency: reuse same subset for k updates.
   * In full mode, do NOT resample and do NOT decrement persist_left. */
  int resample = 0;
  if (!in_full_mode) {
    if (--st->persist_left <= 0) {
      resample = 1;
      st->persist_left = (cfg->persist_k > 0 ? cfg->persist_k : 1);
    }
  } else {
    resample = 0;
  }

  /* 5) Clip threshold */
  double clip = 0.0;
  if (cfg->clip_max_norm > 0.0) {
    clip = cfg->clip_max_norm;
    if (cfg->adaptive_clip && st->ema_gnorm > 0.0) {
      /* clip tracks EMA, e.g., 2x the typical norm */
      clip = fmax(clip, 2.0 * st->ema_gnorm);
    }
  }

  /* 6) Emit directives */
  out->lr             = st->lr_t;
  out->m              = st->m_t;
  out->clip_norm      = clip;                 /* 0 => no clip */
  out->resample_sites = in_full_mode ? 0 : resample;

  /* full-gradient scheduling:
   * - in full mode: always compute full gradient
   * - otherwise: periodic anchor if configured
   */
  if (in_full_mode) {
    out->full_grad_now = 1;
  } else if (cfg->fullgrad_every > 0 && (st->t % cfg->fullgrad_every) == 0) {
    out->full_grad_now = 1;
  } else {
    out->full_grad_now = 0;
  }

  st->t += 1;
}
