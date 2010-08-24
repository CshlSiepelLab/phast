#ifdef PHAST_PTHREAD

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <misc.h>
#include <pthr.h>

#define WRK_FOREACH 0
#define WRK_FORK 1

struct thrpool_ {
  int work_type;

  /* foreach work */
  wrk_func work_func;
  List * work_list;
  int next_unit;
  int unfinished;
  
  /* fork work */
  void * work;
  wrk_func * next_func;
  
  pthread_t * threads;
  int n_threads;
  
  /* locking */
  pthread_cond_t cond_wait_work;
  pthread_mutex_t mtx_wait_work;
  
  pthread_cond_t cond_wait_finish;
  pthread_mutex_t mtx_wait_finish;

  pthread_attr_t attr;
  int pool_quit;
};

static int has_registered = 0;

static void thread_final_exit(void) {
  pthread_exit(NULL);
}

static void * thread_worker(void * ptr) {
  ThreadPool * pool = (ThreadPool*) ptr;
  void * current_work_data = NULL;
  
  /* get work & check for quit */
  pthread_mutex_lock(&(pool->mtx_wait_work));

  /* inform that init has reached wait condition */
  pthread_mutex_lock(&(pool->mtx_wait_finish));
  --(pool->unfinished);
  if (pool->unfinished == 0) {
    pthread_cond_signal(&(pool->cond_wait_finish));
  }    
  pthread_mutex_unlock(&(pool->mtx_wait_finish));

  /* wait for work */
wait_state:
  pthread_cond_wait(&(pool->cond_wait_work),&(pool->mtx_wait_work));
  
  if (pool->pool_quit) {
    pthread_mutex_unlock(&(pool->mtx_wait_work));
    pthread_exit(NULL);
    return NULL;
  }
  
  /* check for work */
check_work:
  if (pool->work_type == WRK_FOREACH) {
    if (pool->next_unit < 0)
      goto wait_state;

    current_work_data = lst_get(pool->work_list, pool->next_unit);
    --(pool->next_unit);

    pthread_mutex_unlock(&(pool->mtx_wait_work));
    
    /* do some work */
    (*(pool->work_func))(current_work_data);
  } else if (pool->work_type == WRK_FORK) {
    wrk_func current_func;
    
    if (*(pool->next_func) == NULL)
      goto wait_state;

    current_func = *(pool->next_func);
    ++(pool->next_func);

    pthread_mutex_unlock(&(pool->mtx_wait_work));
    
    /* do some work */
    (*current_func)(pool->work);
  } else
    fprintf(stderr, "error: invalid thread pool work type.\n");
  
  /* decrement work counter */
  pthread_mutex_lock(&(pool->mtx_wait_finish));
  --(pool->unfinished);
  if (pool->unfinished == 0) {
    pthread_mutex_lock(&(pool->mtx_wait_work)); /* Must acquire this lock before
                                                   sending the signal, or else
                                                   there is a race condition
                                                   between getting to the wait_state
                                                   and the main thread issuing the
                                                   signal that the thread should
                                                   terminate.
                                                 */
    pthread_cond_signal(&(pool->cond_wait_finish));
    pthread_mutex_unlock(&(pool->mtx_wait_finish));

    goto wait_state;
  }    
  pthread_mutex_unlock(&(pool->mtx_wait_finish));
      
  /* go check for next work chunk */
  pthread_mutex_lock(&(pool->mtx_wait_work));
  goto check_work;
}

static void serial_worker(List * work, wrk_func func) {
  int idx;
  int size = lst_size(work);
  
  for (idx = 0; idx < size; ++idx)
    (*func)(lst_get(work, idx));
}

ThreadPool * thr_pool_init(int n_threads) {
  int i;
  ThreadPool * pool = (ThreadPool*) smalloc(sizeof(ThreadPool));
  
  /* clear memory (phast does not have a calloc implementation) */
  memset(pool, 0, sizeof(ThreadPool));
  
  if (n_threads <= 0) {
    pool->n_threads = 0;
    return pool;
  }
  pool->n_threads = n_threads;

  /* clean-up context */
  if (!has_registered)
    atexit(thread_final_exit);

  /* setup finish count */
  pool->unfinished = n_threads;
  
  /* set up locks */
  pthread_mutex_init(&(pool->mtx_wait_work), NULL);
  pthread_cond_init (&(pool->cond_wait_work), NULL);

  pthread_mutex_init(&(pool->mtx_wait_finish), NULL);
  pthread_cond_init (&(pool->cond_wait_finish), NULL);
  
  /* create threads */
  pool->threads = (pthread_t*) smalloc(n_threads*sizeof(pthread_t));
  pthread_attr_init(&(pool->attr));
  pthread_attr_setdetachstate(&(pool->attr), PTHREAD_CREATE_JOINABLE);
  for (i = 0; i < n_threads; ++i) {
    int rc = pthread_create(&(pool->threads[i]), &(pool->attr), thread_worker, pool);
    if (rc) {
      printf("ERROR; return code from pthread_create() is %d\n", rc);
      exit(-1);
    }
  }

  /* wait for threads to finish init */
  pthread_mutex_lock(&(pool->mtx_wait_finish));
  while (pool->unfinished > 0)
    pthread_cond_wait(&(pool->cond_wait_finish),&(pool->mtx_wait_finish));
  pthread_mutex_unlock(&(pool->mtx_wait_finish));
  
  return pool;
}

void thr_pool_free(ThreadPool * pool) {
  int i;
  
  if (pool->n_threads > 0) {
    /* terminate threads */
    /* send termination signal */
    pthread_mutex_lock(&(pool->mtx_wait_work));
    pool->pool_quit = 1;
    pthread_cond_broadcast(&(pool->cond_wait_work));
    pthread_mutex_unlock(&(pool->mtx_wait_work));
    
    /* wait for all threads to finish */
    pthread_attr_destroy(&(pool->attr));
    for (i = 0; i < pool->n_threads; ++i) {
      void *status;
      int rc = pthread_join(pool->threads[i], &status);
      if (rc) {
        printf("ERROR; return code from pthread_join() is %d\n", rc);
        exit(-1);
      }
    }
    
    /* free resources */
    pthread_cond_destroy(&(pool->cond_wait_work));
    pthread_mutex_destroy(&(pool->mtx_wait_work));
    pthread_cond_destroy(&(pool->cond_wait_finish));
    pthread_mutex_destroy(&(pool->mtx_wait_finish));
    
    free(pool->threads);
  }
  free(pool);
}

/* Parallel foreach
  
  Given a list of work units, execute them in parallel using one thread
  for each with a maximum of as many simultaneous threads as those available
  in the thread pool.
  
  Synchronized wait for the work to finish. Function will only return when
  all the work has finished.
*/
void thr_foreach(ThreadPool * pool, List * work, wrk_func func) {
  if (pool->n_threads == 0) {
    serial_worker(work, func);
    return;
  }
  
  /* broadcast work */
  pthread_mutex_lock(&(pool->mtx_wait_work));
  {
    int n_items = lst_size(work);
    
    pool->work_type = WRK_FOREACH;

    /* fill pool data */
    pool->work_list = work;
    pool->next_unit = n_items - 1;
    pool->unfinished = n_items;
    pool->work_func = func;
  }
  pthread_cond_broadcast(&(pool->cond_wait_work));
  pthread_mutex_unlock(&(pool->mtx_wait_work));
  
  /* wait for work to finish */
  pthread_mutex_lock(&(pool->mtx_wait_finish));
  while (pool->next_unit >= 0 || pool->unfinished > 0)
    pthread_cond_wait(&(pool->cond_wait_finish),&(pool->mtx_wait_finish));
  pthread_mutex_unlock(&(pool->mtx_wait_finish));
}

/* Index of current thread within thread pool.

  Given a thread pool, compare the current thread with each of the ThreadPool
  threads and return an index (from 0 to n_threads - 1).
  
  Returns -1 if the current thread does not belong to the thread pool and
  the number of threads is > 0. If the number of threads is zero, it will always
  return 0.
*/
int thr_index(ThreadPool * pool) {
  int i;
  pthread_t current;
  
  if (pool->n_threads == 0)
    return 0;
  
  current = pthread_self();
  
  for (i = 0; i < pool->n_threads; ++i)
    if (pthread_equal(current, pool->threads[i]))
      return i;

  return -1;
}

static void serial_fork(void * work, wrk_func * funcs) {
  while (*funcs != NULL) {
    (*(*funcs))(work);
    ++funcs;
  }
}

/* Thread fork using thread pool.

  Given a null-terminated array of functions and a data pointer, apply each 
  function to the data. Each function runs in its separate thread, with the
  maximum number of simultaneous threads determined by the size of the thread
  pool.
  
  \warning Functions cannot call thr_foreach with the same pool.
  
  Synchronized wait for all the work to finish.
*/
void thr_fork(ThreadPool * pool, void * work, wrk_func * funcs) {
  if (pool->n_threads == 0) {
    serial_fork(work, funcs);
    return;
  }

  /* broadcast work */
  pthread_mutex_lock(&(pool->mtx_wait_work));
  {
    wrk_func * ptr;

    pool->work_type = WRK_FORK;
    
    /* fill pool data */
    pool->next_func = funcs;
    pool->work = work;
    
    /* count units */
    pool->unfinished = 0;
    for (ptr = funcs; *ptr != NULL; ++ptr)
      ++(pool->unfinished);
  }
  pthread_cond_broadcast(&(pool->cond_wait_work));
  pthread_mutex_unlock(&(pool->mtx_wait_work));
  
  /* wait for work to finish */
  pthread_mutex_lock(&(pool->mtx_wait_finish));
  while (*(pool->next_func) != NULL || pool->unfinished > 0)
    pthread_cond_wait(&(pool->cond_wait_finish),&(pool->mtx_wait_finish));
  pthread_mutex_unlock(&(pool->mtx_wait_finish));
}

#endif
