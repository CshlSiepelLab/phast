#ifndef PTHR_H
#define PTHR_H

/** \file pthr.h
    Utilities for parallelizing code based on POSIX Threads.

    A pool of worker threads are used in the functions to carry out
    the requests. The threads remain asleep between function calls
    and can be reused many times in the program.

    \note Thread pools are not reentrant and so cannot be reused from
          within worker functions.
    \ingroup base
*/

#include <pthread.h>
#include <lists.h>

/** Worker function prototype. The argument is the value returned by
    lst_get() function and will require conversion to the appropriate
    type. See lists.h for how to do this.
*/  
typedef void (*wrk_func) (void * work_unit);

/** Object representing the thread pool. */
typedef struct thrpool_ ThreadPool;

/** Create a new thread pool. This function creates a new thread
    pool and initializes the specified number of threads.

    If the specified number of threads is zero, then no threads are
    created and running code is executed in a serial fashion. May be
    useful for debugging.
*/
ThreadPool * thr_pool_init(int n_threads);

/** Terminate thread pool and release resources. */
void thr_pool_free(ThreadPool * pool);

/** Parallel foreach using thread pool.
  
  Given a list of work units, execute them in parallel using one thread
  for each with a maximum of as many simultaneous threads as those available
  in the thread pool.
  
  Synchronized wait for the work to finish. Function will only return when
  all the work has finished.

  When the thread pool has zero threads, this is equivalent to a for loop
  on the list.
*/
void thr_foreach(ThreadPool * pool, List * work, wrk_func func);

/** Index of current thread within thread pool.

  Given a thread pool, compare the current thread with each of the ThreadPool
  threads and return an index (from 0 to n_threads - 1).
  
  Returns -1 if the current thread does not belong to the thread pool and
  the number of threads is > 0. If the number of threads is zero, it will always
  return 0.
*/
int thr_index(ThreadPool * pool);

/** Thread fork using thread pool.

  Given a null-terminated array of functions and a data pointer, apply each 
  function to the data. Each function runs in its separate thread, with the
  maximum number of simultaneous threads determined by the size of the thread
  pool.
  
  \warning Worker functions cannot call thr_foreach with the same pool.
  
  Synchronized wait for all the work to finish.
*/
void thr_fork(ThreadPool * pool, void * work, wrk_func * funcs);

#endif
