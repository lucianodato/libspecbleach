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

#ifndef SHARED_UTILS_THREAD_POOL_H
#define SHARED_UTILS_THREAD_POOL_H

#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Fixed-size worker pool for deterministic fork/join parallel loops.
 *
 * Workers are created once at pool creation and block between dispatches, so
 * dispatching from a real-time audio thread performs no allocation, no lock
 * acquisition, and no thread management. Work is split into contiguous static
 * ranges, which keeps the per-element accumulation order identical to
 * single-threaded execution.
 */
typedef struct SbThreadPool SbThreadPool;

/**
 * @brief Creates a pool with the specified number of worker threads.
 *
 * The calling thread always participates in work execution, so the effective
 * parallelism of a dispatch is `num_workers + 1`. Returns NULL on failure or
 * if `num_workers` is zero; callers should fall back to sequential execution.
 */
SbThreadPool* sb_thread_pool_create(uint32_t num_workers);

/**
 * @brief Returns the number of worker threads in the pool (excludes the
 * calling thread).
 */
uint32_t sb_thread_pool_num_workers(const SbThreadPool* pool);

/**
 * @brief Runs `task(arg, start, count)` over `[0, num_items)` split into
 * contiguous ranges.
 *
 * The calling thread executes the first range synchronously and joins the
 * workers before returning, so on return all items are fully processed. When
 * the pool is NULL, `num_items` is zero, or there is a single range, the task
 * runs sequentially on the calling thread.
 */
void sb_thread_pool_parallel_for(SbThreadPool* pool, uint32_t num_items,
                                 void (*task)(void* arg, uint32_t start,
                                              uint32_t count),
                                 void* arg);

/**
 * @brief Stops all workers and releases the pool. Passing NULL is a no-op.
 */
void sb_thread_pool_free(SbThreadPool* pool);

#ifdef __cplusplus
}
#endif

#endif /* SHARED_UTILS_THREAD_POOL_H */
