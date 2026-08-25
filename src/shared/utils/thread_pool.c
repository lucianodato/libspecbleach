/*
libspecbleach - A spectral processing library

Copyright 2022 Luciano Dato <lucianodato@gmail.com>

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include "thread_pool.h"

#include <stdatomic.h>
#include <stdlib.h>

#include "simd_utils.h"

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <process.h>
#include <windows.h>
#elif defined(__APPLE__)
// macOS does not implement unnamed POSIX semaphores (sem_init fails with
// ENOSYS); Grand Central Dispatch semaphores are the native equivalent.
#include <dispatch/dispatch.h>
#include <pthread.h>
#include <sched.h>
#else
#include <pthread.h>
#include <sched.h>
#include <semaphore.h>
#endif

// Each worker owns a private semaphore: every dispatch posts exactly one
// token per participating worker, so wakeups are exactly-once and a worker
// can never steal another worker's dispatch.
typedef struct {
  SbThreadPool* pool;
  uint32_t participant; // Range index: caller is 0, workers are 1..N
#ifdef _WIN32
  HANDLE sem;
#elif defined(__APPLE__)
  dispatch_semaphore_t sem;
#else
  sem_t sem;
#endif
} SbWorker;

typedef struct {
  void (*task)(void* arg, uint32_t start, uint32_t count);
  void* arg;
  uint32_t num_items;
  uint32_t base_count;   // Minimum range size per participant
  uint32_t remainder;    // First `remainder` participants get one extra item
  uint32_t participants; // Ranges executed per dispatch (caller included)
} SbWorkItem;

struct SbThreadPool {
  atomic_uint remaining; // Workers that have not finished the current dispatch
  atomic_bool shutdown;  // Set once on destruction
  SbWorkItem work;       // Published to workers via semaphore synchronization

  uint32_t worker_count;
  uint32_t threads_created;
  uint32_t semaphores_initialized; // Workers with a live semaphore
  SbWorker* workers;
#ifdef _WIN32
  HANDLE* threads;
#else
  pthread_t* threads;
#endif
};

static void sb_yield_cpu(void) {
#ifdef _WIN32
  SwitchToThread();
#else
  sched_yield();
#endif
}

static inline void sb_wake_worker(SbWorker* worker) {
#ifdef _WIN32
  ReleaseSemaphore(worker->sem, 1, NULL);
#elif defined(__APPLE__)
  dispatch_semaphore_signal(worker->sem);
#else
  sem_post(&worker->sem);
#endif
}

static inline void sb_worker_wait(SbWorker* worker) {
#ifdef _WIN32
  WaitForSingleObject(worker->sem, INFINITE);
#elif defined(__APPLE__)
  dispatch_semaphore_wait(worker->sem, DISPATCH_TIME_FOREVER);
#else
  while (sem_wait(&worker->sem) != 0) {
    // Retry: sem_wait only fails here with EINTR
  }
#endif
}

static inline uint32_t sb_range_start(uint32_t p, uint32_t base_count,
                                      uint32_t remainder) {
  return p * base_count + (p < remainder ? p : remainder);
}

static inline uint32_t sb_range_count(uint32_t base_count, uint32_t has_extra) {
  return base_count + has_extra;
}

static void sb_worker_execute(SbWorker* worker) {
  SbThreadPool* pool = worker->pool;
  const uint32_t count =
      sb_range_count(pool->work.base_count,
                     (worker->participant < pool->work.remainder) ? 1U : 0U);
  const uint32_t start = sb_range_start(
      worker->participant, pool->work.base_count, pool->work.remainder);
  if (count > 0U && start + count <= pool->work.num_items) {
    pool->work.task(pool->work.arg, start, count);
  }
  atomic_fetch_sub_explicit(&pool->remaining, 1U, memory_order_acq_rel);
}

#ifdef _WIN32
static unsigned __stdcall sb_worker_main(void* raw) {
#else
static void* sb_worker_main(void* raw) {
#endif
  SbWorker* worker = (SbWorker*)raw;

  // Dedicated library threads run with denormals flushed for their lifetime.
  sb_simd_enable_ftz_daz();

  for (;;) {
    sb_worker_wait(worker);
    if (atomic_load_explicit(&worker->pool->shutdown, memory_order_acquire)) {
      break;
    }
    sb_worker_execute(worker);
  }

#ifdef _WIN32
  return 0U;
#else
  return NULL;
#endif
}

static inline int sb_worker_spawn(SbThreadPool* pool, uint32_t index) {
#ifdef _WIN32
  const uintptr_t handle =
      _beginthreadex(NULL, 0U, sb_worker_main, &pool->workers[index], 0U, NULL);
  pool->threads[index] = (HANDLE)handle;
  return pool->threads[index] ? 0 : -1;
#else
  return pthread_create(&pool->threads[index], NULL, sb_worker_main,
                        &pool->workers[index]);
#endif
}

SbThreadPool* sb_thread_pool_create(uint32_t num_workers) {
  if (num_workers == 0U) {
    return NULL;
  }

  SbThreadPool* pool = (SbThreadPool*)calloc(1U, sizeof(SbThreadPool));
  if (!pool) {
    return NULL;
  }

  atomic_init(&pool->remaining, 0U);
  atomic_init(&pool->shutdown, false);
  pool->worker_count = num_workers;

  pool->workers = (SbWorker*)calloc(num_workers, sizeof(SbWorker));
  pool->threads = (void*)calloc(num_workers, sizeof(*pool->threads));
  if (!pool->workers || !pool->threads) {
    sb_thread_pool_free(pool);
    return NULL;
  }

  for (uint32_t i = 0U; i < num_workers; i++) {
    SbWorker* worker = &pool->workers[i];
    worker->pool = pool;
    worker->participant = i + 1U;
#ifdef _WIN32
    worker->sem = CreateSemaphore(NULL, 0, 1, NULL);
    if (!worker->sem) {
      sb_thread_pool_free(pool);
      return NULL;
    }
#elif defined(__APPLE__)
    worker->sem = dispatch_semaphore_create(0);
    if (!worker->sem) {
      sb_thread_pool_free(pool);
      return NULL;
    }
#else
    if (sem_init(&worker->sem, 0, 0U) != 0) {
      sb_thread_pool_free(pool);
      return NULL;
    }
#endif
    pool->semaphores_initialized++;
  }

  for (uint32_t i = 0U; i < num_workers; i++) {
    if (sb_worker_spawn(pool, i) != 0) {
      sb_thread_pool_free(pool);
      return NULL;
    }
    pool->threads_created++;
  }

  return pool;
}

uint32_t sb_thread_pool_num_workers(const SbThreadPool* pool) {
  return pool ? pool->worker_count : 0U;
}

void sb_thread_pool_parallel_for(SbThreadPool* pool, uint32_t num_items,
                                 void (*task)(void* arg, uint32_t start,
                                              uint32_t count),
                                 void* arg) {
  if (!pool || num_items == 0U || !task) {
    if (task && num_items > 0U) {
      task(arg, 0U, num_items);
    }
    return;
  }

  // Caller plus workers; never more participants than items.
  uint32_t participants = pool->worker_count + 1U;
  if (participants > num_items) {
    participants = num_items;
  }
  const uint32_t workers = participants - 1U;
  const uint32_t base_count = num_items / participants;
  const uint32_t remainder = num_items - base_count * participants;

  pool->work.task = task;
  pool->work.arg = arg;
  pool->work.num_items = num_items;
  pool->work.base_count = base_count;
  pool->work.remainder = remainder;
  pool->work.participants = participants;

  atomic_store_explicit(&pool->remaining, workers, memory_order_relaxed);

  // Semaphore synchronization publishes the work item above to every woken
  // worker. Each worker owns a private token, so exactly `workers` distinct
  // workers execute exactly once.
  for (uint32_t i = 0U; i < workers; i++) {
    sb_wake_worker(&pool->workers[i]);
  }

  // The caller executes the first range, then joins the workers.
  const uint32_t caller_count =
      sb_range_count(base_count, (remainder > 0U) ? 1U : 0U);
  if (caller_count > 0U) {
    task(arg, 0U, caller_count);
  }

  while (atomic_load_explicit(&pool->remaining, memory_order_acquire) != 0U) {
    sb_yield_cpu();
  }
}

void sb_thread_pool_free(SbThreadPool* pool) {
  if (!pool) {
    return;
  }

  if (pool->threads && pool->workers) {
    atomic_store_explicit(&pool->shutdown, true, memory_order_release);
    for (uint32_t i = 0U; i < pool->threads_created; i++) {
      sb_wake_worker(&pool->workers[i]);
    }
    for (uint32_t i = 0U; i < pool->threads_created; i++) {
#ifdef _WIN32
      WaitForSingleObject(pool->threads[i], INFINITE);
      CloseHandle(pool->threads[i]);
#else
      pthread_join(pool->threads[i], NULL);
#endif
    }
  }

  if (pool->workers) {
    for (uint32_t i = 0U; i < pool->semaphores_initialized; i++) {
#ifdef _WIN32
      if (pool->workers[i].sem) {
        CloseHandle(pool->workers[i].sem);
      }
#elif defined(__APPLE__)
      if (pool->workers[i].sem) {
        dispatch_release(pool->workers[i].sem);
      }
#else
      sem_destroy(&pool->workers[i].sem);
#endif
    }
  }

  free(pool->workers);
  free(pool->threads);
  free(pool);
}
