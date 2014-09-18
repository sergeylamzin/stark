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


#ifndef STARK_STATUS_H
#define STARK_STATUS_H
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <sys/syscall.h>

struct stark_status_s {
	pthread_mutex_t mutex;
	// thread_pool_t threadpool;
	struct distributor_s * distributor;
	volatile struct stark_thread_status_s * stark_thread_status;
	struct stark_hierarchical_assembler_s * hierarchical_assembler;
	// stark_t * stark;
	size_t reads_total;
};

struct stark_thread_status_s {
	volatile struct stark_thread_status_s * next;
	pthread_t thread_id;
	size_t reads_inserted;
	// const char* current_action;
	
	struct {
		char * seq4;
		size_t length;
		size_t maxk;
	} current_insert;
	
	size_t total_insert_mis;
	size_t total_insertions;
	int last_cpu;
	pid_t linux_tid;
	struct stark_coverage_mytoken_s * mytoken;
	struct stark_hierarchical_assembler_test_and_merge_groups_dispatch_s * volatile dispatch;
	char current_action[128];
};

int stark_status_init(uintptr_t port);
void stark_status_kill();
void stark_status_register_thread();

// const char stark_status_action_waiting_for_jobs[] = "Waiting for Jobs.";
// const char stark_status_action_waiting_for_barrier[] = "Waiting on Barrier Job.";
extern const char stark_status_action_inserting_read[];

extern volatile struct stark_status_s stark_status;




// #ifndef __MACH__



#define get_stark_thread_status() stark_thread_status

// #else
// static inline volatile struct stark_thread_status_s * get_stark_thread_status() {
// 	volatile struct stark_thread_status_s * next;
// 	for (next = stark_status.stark_thread_status; next; next = next->next) {
// 		if (next->thread_id == pthread_self())
// 			return next;
// 	}
// 	return NULL;
// }
// #endif

// #ifdef SYS_getcpu
// static inline int getcpu() {
// 	int cpu, status;
// 	status = syscall(SYS_getcpu, &cpu, NULL, NULL);
// 	return (status == -1) ? status : cpu;
// }
// 
// static inline void stark_status_update_last_cpu() {
// 	if (get_stark_thread_status())
// 		get_stark_thread_status()->last_cpu = getcpu();
// }
// #else
// #define stark_status_update_last_cpu()
// #endif

// #define stark_set_status_action(...) if (get_stark_thread_status()) 
	// snprintf((char *)(get_stark_thread_status()->current_action), sizeof(get_stark_thread_status()->current_action), __VA_ARGS__)
	
#ifndef STARK_MAIN_STATUS_C
#define STARK_MAIN_STATUS_C extern
#endif

STARK_MAIN_STATUS_C struct stark_thread_status_s * stark_thread_status;
#pragma omp threadprivate(stark_thread_status)

#include <stdarg.h>

static inline void stark_set_status_action(const char *fmt, ...) {
	if (get_stark_thread_status())
	{
		va_list args;
		va_start(args, fmt);
		vsnprintf((char *)(get_stark_thread_status()->current_action), sizeof(get_stark_thread_status()->current_action), fmt, args);
		va_end(args);
	}
}
// #endif

// static inline void stark_set_status_action(const char * status) {
// 	if (get_stark_thread_status()) {
// 		
// 	}
// 		get_stark_thread_status()->current_action = status;
// }

#endif
