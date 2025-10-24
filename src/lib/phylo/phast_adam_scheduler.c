#include <math.h>
#include <string.h>
#include <phast/misc.h>
#include <phast/adam_scheduler.h>

/* set up a scheduler with appropriate defaults */
Scheduler* sched_new(int N_sites, int init_subsample, int inc_every,
                     double init_lr, int persist_k, int fullgrad_every) {
  Scheduler *s = smalloc(sizeof(Scheduler));
  s->N_sites = N_sites;
  s->m0 = init_subsample;
  s->inc_every = inc_every;
  s->lr0 = init_lr;
  s->noise_const = s->lr0 * s->m0;
  s->keep_noise_const = TRUE;
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
  st->lr_t = cfg->lr0;
  st->lambda_full = 0.0;
  return st;
}

/* cosine helper for smooth ramps [0,1] over total steps */
static inline double cosine01(int t, int T) {
  if (T <= 0) return 0.0;
  double x = (double)t / (double)T;
  if (x < 0) x = 0; if (x > 1) x = 1;
  return 0.5 * (1.0 - cos(M_PI * x));
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

  /* 3) Learning-rate policy tied to batch to control noise scale */
  if (cfg->keep_noise_const) {
    /* keep lr * m ≈ constant => lr_t = noise_const / m_t */
    st->lr_t = cfg->noise_const / (double)st->m_t;
  } else {
    /* independent; example: mild cosine decay over run */
    double decay = 0.1 + 0.9 * (1.0 - cosine01(st->t, cfg->T_total));
    st->lr_t = cfg->lr0 * decay;
  }

  /* 4) Blending weight for occasional full-gradient anchoring */
  if (cfg->fullgrad_every > 0) {
    /* Option A: ramp lambda_full over the latter half of training */
    double ramp = cosine01(st->t - cfg->T_total/2, cfg->T_total/2);
    if (ramp < 0) ramp = 0; if (ramp > 1) ramp = 1;
    st->lambda_full = ramp;
  } else {
    st->lambda_full = 0.0;
  }

  /* 5) Site subset persistency: reuse same subset for k steps */
  int resample = 0;
  if (--st->persist_left <= 0) {
    resample = 1;
    st->persist_left = (cfg->persist_k > 0 ? cfg->persist_k : 1);
  }

  /* 6) Clip threshold */
  double clip = 0.0;
  if (cfg->clip_max_norm > 0.0) {
    clip = cfg->clip_max_norm;
    if (cfg->adaptive_clip && st->ema_gnorm > 0.0) {
      /* clip tracks EMA, e.g., 2x the typical norm */
      clip = fmax(clip, 2.0 * st->ema_gnorm);
    }
  }

  /* 7) Emit directives */
  out->lr            = st->lr_t;
  out->m             = st->m_t;
  out->clip_norm     = clip;          // 0 => do nothing
  out->resample_sites= resample;
  out->lambda_full   = st->lambda_full;

  st->t += 1;
}

/* Global-norm clipping (L2). Returns scale factor applied (<=1). */
double sched_clip_scale(double *g, int n, double clip_norm) {
  if (clip_norm <= 0.0) return 1.0;
  double s = 0.0;
  for (int i = 0; i < n; ++i) s += g[i]*g[i];
  double gn = sqrt(s);
  if (gn <= clip_norm || gn == 0.0) return 1.0;
  double scale = clip_norm / gn;
  for (int i = 0; i < n; ++i) g[i] *= scale;
  return scale;
}
