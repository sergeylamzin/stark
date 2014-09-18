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


#include "pthread.h"
#include "main.h"

#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#define THREAD_POOL_JOB_BROADCAST 3
#define THREAD_POOL_JOB_BARRIER 3
#define THREAD_POOL_JOB_PROTECT_DESTRUCTOR 4

typedef struct {
	size_t job_id;
	void* (*function)(void*);
	void* arg;
	void (*destroy)(void*);
	int flags;
	pthread_cond_t * barrier_cond;
	// pthread_mutex_t barrier_mutex;
	volatile size_t numthreads_waiting;
	volatile size_t numjobs_completed_after_broadcast;
} job_t;

struct thread_pool_s {
	// sem_t semaphore;
	pthread_cond_t jobs_available;
	pthread_cond_t pool_idle;
	pthread_mutex_t mutex;
	volatile job_t* jobs;
	size_t length;
	volatile size_t next;
	volatile size_t size;
	size_t numthreads;
	volatile size_t running_threads;
	// volatile size_t barrier_waiting_threads;
	pthread_t* threads;
	size_t threads_length;
	volatile bool shutdown;
};

typedef struct thread_pool_s* thread_pool_t;

thread_pool_t thread_pool_create(int thread_count);
int thread_pool_dispatch(struct thread_pool_s* thread_pool, void* (*function)(void*), void* arg, void (*destroy)(void*), int flags); 
// int thread_pool_dispatch_barrier(thread_pool_t thread_pool, void* (*function)(void*), void* arg, void (*destroy)(void*));

void thread_pool_wait(thread_pool_t thread_pool);
int thread_pool_shutdown(thread_pool_t thread_pool); 
int thread_pool_purge(thread_pool_t thread_pool); 
// size_t thread_pool_join(thread_pool_t thread_pool);

#endif
