#ifndef PTHR_H
#define PTHR_H

/** \file pthr.h
    Utilities for parallelizing code based on POSIX Threads.

    A pool of worker threads are used in the functions to carry out
    the requests. The threads remain asleep between function calls
    and can be reused many times in the program.

    Note: Thread pools are not reentrant and so cannot be called from
          different threads.
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

#endif
