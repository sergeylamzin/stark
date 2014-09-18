/*
    Copyright (C) 2014  Sergey Lamzin, https://github.com/sergeylamzin/stark

    This file is part of the StarK genome assembler.

    StarK is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    StarK is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

*/


// #include "threadpool.h"
#include <pthread.h>
// #include <semaphore.h>
#include "main.h"
#include "threadpool.h"
#include <errno.h>
#include <string.h>
#include "stark_status.h"
#include "stark_assemble_hierarchical.h"

// #ifdef UV_SYSTEM
// #include "uv_system.h"
// #endif

// sem_t *sem_open(const char *name, int oflag);


size_t thread_pool_join(struct thread_pool_s* thread_pool) {
	size_t jobs_done = 0;
	pthread_mutex_lock(&thread_pool->mutex);
	if (thread_pool->threads_length <= thread_pool->numthreads) {
		void* tmp = realloc(thread_pool->threads, sizeof(pthread_t) * thread_pool->threads_length * 2);
		if (!tmp) 
			return 0;
			
		thread_pool->threads = tmp;
		thread_pool->threads_length *= 2;
	}
	thread_pool->threads[thread_pool->numthreads++] = pthread_self();
	// thread_pool->numthreads++;
	
	// stark_status_register_thread();
	
	// #ifdef UV_SYSTEM
	// if (uv_init_atomic()) {
	// 	CRITICAL_ERROR("Cannot initialize UV GRUs on thread 0x%lX\n", pthread_self());
	// }
	// pthread_cleanup_push(uv_free_atomic, NULL);
	// #endif
	
	job_t job;
	for (;;) {
		// int value; 
		// sem_getvalue(&thread_pool->semaphore, &value);
		// DEBUG_MSG("Waiting for Job. size = %lu, sem = %d", thread_pool->size, value);
		// sem_wait(&thread_pool->semaphore);
		// DEBUG_MSG("thread_pool->size = %zd, thread_pool->running_threads = %zd, thread_pool->numthreads = %zd", thread_pool->size, thread_pool->running_threads, thread_pool->numthreads);
		if (!(thread_pool->size || thread_pool->running_threads))
			pthread_cond_broadcast(&thread_pool->pool_idle);
		
		// DEBUG_MSG("thread_pool->size = %d, thread_pool->running_threads = %d, thread_pool->numthreads = %d, pthread_self() = 0x%lx", thread_pool->size, thread_pool->running_threads, thread_pool->numthreads, pthread_self());
		while (!(thread_pool->size || thread_pool->shutdown)) {
			#if defined(STARK_STATUS_H)
			stark_set_status_action("Waiting for jobs");
			#endif
			pthread_cond_wait(&thread_pool->jobs_available,&thread_pool->mutex);
		}
			
		// DEBUG_MSG("Executing for Job. size = %lu, next = %d", thread_pool->size, thread_pool->next);
		if (thread_pool->shutdown) {
			#if defined(STARK_STATUS_H)
			stark_set_status_action("Shutting down");
			#endif
			if (thread_pool->numthreads <= 1) {
				// im the last thread to shut down
				// destroy the pool;
				pthread_cond_broadcast(&thread_pool->pool_idle);
				
				// DEBUG_MSG("size = %lu, im the last thread to shut down, destroying the pool.", thread_pool->size);
				pthread_mutex_unlock(&thread_pool->mutex);
				
				while (pthread_cond_destroy(&thread_pool->jobs_available) == EBUSY) {
				        pthread_cond_broadcast(&thread_pool->jobs_available);
				        sched_yield();
				}
				
				while (pthread_cond_destroy(&thread_pool->pool_idle) == EBUSY) {
				        pthread_cond_broadcast(&thread_pool->pool_idle);
				        sched_yield();
				}
				
				while (pthread_mutex_destroy(&thread_pool->mutex) == EBUSY) {
				        sched_yield();
				}
				
				free((void*)thread_pool->jobs);
				
				// pthread_mutex_destroy(&thread_pool->mutex);
				
				// pthread_cond_destroy(&thread_pool->jobs_available);
				// pthread_cond_destroy(&thread_pool->pool_idle);
				// sem_destroy(&thread_pool->semaphore);
				free(thread_pool);
			} else {
				thread_pool->numthreads--;
				pthread_mutex_unlock(&thread_pool->mutex);
			}
			
			break;
		}
		job = thread_pool->jobs[thread_pool->next];
		
		if (!(job.flags & THREAD_POOL_JOB_BARRIER || job.flags & THREAD_POOL_JOB_BROADCAST)) {
			thread_pool->next = (thread_pool->next + 1) % thread_pool->length;
			thread_pool->size--;
		}	
		thread_pool->running_threads++;
		// DEBUG_MSG("Thread %lu picked up job", pthread_self());
		pthread_mutex_unlock(&thread_pool->mutex);
		#if defined(STARK_STATUS_H)
		stark_set_status_action("Executing job");
		// stark_status_update_last_cpu();
		#endif
		if (job.function)
			job.function(job.arg);
		// if (job.destroy)
		// 	job.destroy(job.arg);
		jobs_done++;
		
		
			
		
		if (job.flags & THREAD_POOL_JOB_BARRIER || job.flags & THREAD_POOL_JOB_BROADCAST) {
			pthread_mutex_lock(&thread_pool->mutex);
			thread_pool->jobs[thread_pool->next].numthreads_waiting++;
			if (thread_pool->jobs[thread_pool->next].numthreads_waiting < thread_pool->numthreads) {
				// wait for the barrier
				// DEBUG_MSG("Thread %lu waiting for barrier", pthread_self());
				// 
				#if defined(STARK_STATUS_H)
				stark_set_status_action("Waiting for barrier");
				#endif
				
				pthread_cond_wait(job.barrier_cond, &thread_pool->mutex);
				// skip destroying the job
				goto end_job;
			} else {
				// I am the last one to complete the broadcast barrier job
				// Remove the boradcast job
				thread_pool->next = (thread_pool->next + 1) % thread_pool->length;
				thread_pool->size--;
				// Free everyone and destroy the condition
				// DEBUG_MSG("Thread 0x%lX releasing barrier 0x%lX", pthread_self(), job.barrier_cond);
				pthread_mutex_unlock(&thread_pool->mutex);
				do {
					pthread_cond_broadcast(job.barrier_cond);
					sched_yield();
				} while (pthread_cond_destroy(job.barrier_cond) == EBUSY);
				free(job.barrier_cond);
				// pthread_mutex_lock(&thread_pool->mutex);
			}
		}
		// destroy argument
		if (job.destroy && (job.flags & THREAD_POOL_JOB_PROTECT_DESTRUCTOR)) {
			pthread_mutex_lock(&thread_pool->mutex);
			job.destroy(job.arg);
		} else {
			if (job.destroy) {
				job.destroy(job.arg);
			}
			pthread_mutex_lock(&thread_pool->mutex);
		}
		
		end_job:
			
		// DEBUG_MSG("Thread %lu finished job", pthread_self());
		thread_pool->running_threads--;
	}
	
	// #ifdef UV_SYSTEM
	// pthread_cleanup_pop(1);
	// #endif
	#if defined(STARK_STATUS_H)
	stark_set_status_action("dead");
	#endif
	
	return jobs_done;
}

struct thread_pool_s* thread_pool_create(int thread_count) {
	struct thread_pool_s * thread_pool = calloc(1, sizeof(struct thread_pool_s));
	if (!thread_pool)
		CRITICAL_ERROR("calloc(1, sizeof(struct thread_pool_s)) failed");
	thread_pool->jobs = calloc(8, sizeof(job_t));
	thread_pool->length = 8;
	thread_pool->running_threads = 0;
	thread_pool->numthreads = 0;
	// int error;
	// if(error = sem_init(&thread_pool->semaphore, 0, 0)) {
	// 	CRITICAL_ERROR("error initializing semaphore. errno = %d\n%s", errno, strerror(errno));
	// }
	pthread_mutex_init(&thread_pool->mutex, NULL);
	pthread_cond_init(&thread_pool->jobs_available, NULL);
	pthread_cond_init(&thread_pool->pool_idle, NULL);
	// pthread_cond_init(&thread_pool->barrier, NULL);
	thread_pool->threads = calloc(thread_count+1, sizeof(pthread_t));
	thread_pool->threads_length = thread_count+1;
	pthread_mutex_lock(&thread_pool->mutex);
	int i;
	pthread_t thread;
	for (i = 0; i < thread_count; i++) {
		pthread_create(&thread, 0, (void *(*)(void*))thread_pool_join, thread_pool);
	}
	pthread_mutex_unlock(&thread_pool->mutex);
	DEBUG_MSG("Thread Pool Initialized with %d threads.",thread_count);
	return thread_pool;
}

static inline int thread_pool_dispatch_nolock(struct thread_pool_s* thread_pool, void* (*function)(void*), void* arg, void (*destroy)(void*), int flags) {
	if (thread_pool->size >= thread_pool->length) {
		void* tmp = realloc((void*)thread_pool->jobs, sizeof(job_t) * thread_pool->length * 2);
		if (tmp) {
			thread_pool->jobs = tmp;
			memmove ((void *)(thread_pool->jobs + thread_pool->length), (void*)thread_pool->jobs, thread_pool->next * sizeof(job_t) );
			thread_pool->length *= 2;
		}
		else
			return -1;
	}
	size_t newpos = (thread_pool->next + thread_pool->size) % thread_pool->length;
	thread_pool->jobs[newpos].function = function;
	thread_pool->jobs[newpos].arg = arg;
	thread_pool->jobs[newpos].destroy = destroy;
	thread_pool->jobs[newpos].flags = flags;
	
	thread_pool->size++;
	
	if (flags & THREAD_POOL_JOB_BARRIER || flags & THREAD_POOL_JOB_BROADCAST) {
		thread_pool->jobs[newpos].flags |= THREAD_POOL_JOB_PROTECT_DESTRUCTOR;
		pthread_cond_init((thread_pool->jobs[newpos].barrier_cond = calloc(1, sizeof(pthread_cond_t))), NULL);
		// pthread_mutex_init(&thread_pool->jobs[newpos].barrier_mutex, NULL);
		thread_pool->jobs[newpos].numthreads_waiting = 0;
		pthread_cond_broadcast(&thread_pool->jobs_available);
	}
	else
		pthread_cond_signal(&thread_pool->jobs_available);
	
	return 0;
}

int thread_pool_dispatch(struct thread_pool_s* thread_pool, void* (*function)(void*), void* arg, void (*destroy)(void*), int flags) {
	
	pthread_mutex_lock(&thread_pool->mutex);
	
	thread_pool_dispatch_nolock(thread_pool, function, arg, destroy, flags);
	
	pthread_mutex_unlock(&thread_pool->mutex);
	
	return 0;
}

// static void do_nothing() {
// 	__asm__ __volatile__ ("");
// }

int thread_pool_dispatch_barrier(struct thread_pool_s* thread_pool, void* (*function)(void*), void* arg, void (*destroy)(void*)) {
	pthread_mutex_lock(&thread_pool->mutex);
	
	thread_pool_dispatch_nolock(thread_pool, function, arg, destroy, THREAD_POOL_JOB_BARRIER | THREAD_POOL_JOB_BROADCAST);
	
	pthread_mutex_unlock(&thread_pool->mutex);
	
	return 0;
}

void* thread_pool_shutdown_real(struct thread_pool_s* thread_pool) {
	pthread_mutex_lock(&thread_pool->mutex);
	thread_pool->shutdown = 1;
	pthread_cond_broadcast(&thread_pool->jobs_available);
	pthread_mutex_unlock(&thread_pool->mutex);
	
	return 0;
}

// void * thread_pool_shutdown_thread(struct thread_pool_s* thread_pool) {
// 	pthread_mutex_lock(&thread_pool->mutex);
// 	if (thread_pool->numthreads > 1) {
// 		thread_pool_dispatch_nolock(thread_pool, (void* (*)(void*))thread_pool_shutdown_thread, thread_pool, NULL);
// 		thread_pool->numthreads--;
// 		pthread_mutex_unlock(&thread_pool->mutex);
// 
// 		pthread_exit(NULL);
// 	} else { // I am the last thread to shut down
// 		pthread_mutex_unlock(&thread_pool->mutex);
// 		
// 		while (pthread_cond_destroy(&thread_pool->jobs_available) == EBUSY) {
// 		        pthread_cond_broadcast(&thread_pool->jobs_available);
// 		        sched_yield();
// 		}
// 		
// 		while (pthread_cond_destroy(&thread_pool->pool_idle) == EBUSY) {
// 		        pthread_cond_broadcast(&thread_pool->pool_idle);
// 		        sched_yield();
// 		}
// 		
// 		while (pthread_mutex_destroy(&thread_pool->mutex) == EBUSY) {
// 		        sched_yield();
// 		}
// 		
// 		free((void*)thread_pool->jobs);
// 		
// 		// pthread_mutex_destroy(&thread_pool->mutex);
// 		
// 		// pthread_cond_destroy(&thread_pool->jobs_available);
// 		// pthread_cond_destroy(&thread_pool->pool_idle);
// 		// sem_destroy(&thread_pool->semaphore);
// 		free(thread_pool);
// 	}
// 	
// 	
// }

int thread_pool_shutdown(struct thread_pool_s* thread_pool) {
	thread_pool_dispatch(thread_pool, (void* (*)(void*))thread_pool_shutdown_real, thread_pool, NULL, 0);
	// pthread_mutex_lock(&thread_pool->mutex);
	// thread_pool->shutdown = 1;
	// pthread_cond_broadcast(&thread_pool->jobs_available);
	// pthread_mutex_unlock(&thread_pool->mutex);
	
	return 0;
}

int thread_pool_purge(struct thread_pool_s* thread_pool) {
	pthread_mutex_lock(&thread_pool->mutex);
	job_t job;
	while(thread_pool->size) {
		job = thread_pool->jobs[thread_pool->next];
		thread_pool->next = (thread_pool->next + 1) % thread_pool->length;
		thread_pool->size--;
		if (job.destroy)
			job.destroy(job.arg);
	}
	// thread_pool->size = 0;
	// thread_pool->next = 0;
	pthread_mutex_unlock(&thread_pool->mutex);
	// pthread_mutex_destroy(&thread_pool->mutex);
	
	return 0;
}

int thread_pool_destroy(struct thread_pool_s* thread_pool);

void thread_pool_wait(struct thread_pool_s* thread_pool) {
	pthread_mutex_lock(&thread_pool->mutex);
	// DEBUG_MSG("thread_pool->size = %zd, thread_pool->running_threads = %zd, thread_pool->numthreads = %zd", thread_pool->size, thread_pool->running_threads, thread_pool->numthreads);
	while (thread_pool->size || thread_pool->running_threads)
		pthread_cond_wait(&thread_pool->pool_idle, &thread_pool->mutex);
	
	pthread_mutex_unlock(&thread_pool->mutex);
}


