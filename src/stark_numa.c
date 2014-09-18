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


// #define _GNU_SOURCE
#include <pthread.h>
// #include <sched.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <sched.h>
#include <unistd.h>
#include <sys/syscall.h>

#ifdef USE_LIBNUMA
#include <numa.h>
#include <numaif.h>
#endif

// CPU_SETSIZE

static int initialized = 0;
cpu_set_t allcpus = {0};
int master_node = -1;
// nodemask_t master_nodemask;


#ifndef CPU_COUNT

int CPU_COUNT(cpu_set_t * set) {
	int count = 0;
	size_t i;
	for (i = 0; i < (sizeof(set->__bits) / sizeof(set->__bits[0])); i++) {
		__cpu_mask l;
		// #if defined(__x86_64__)
		// __asm__ ("popcnt %1, %0"
		// 		: "=r" (l)
		// 		: "g" (set->__bits[i]));
		// #else
			l = set->__bits[i];
			if (!l)
				continue;
			
			# if LONG_BIT > 32
				l = (l & 0x5555555555555555ul) + ((l >> 1) & 0x5555555555555555ul);
				l = (l & 0x3333333333333333ul) + ((l >> 2) & 0x3333333333333333ul);
				l = (l & 0x0f0f0f0f0f0f0f0ful) + ((l >> 4) & 0x0f0f0f0f0f0f0f0ful);
				l = (l & 0x00ff00ff00ff00fful) + ((l >> 8) & 0x00ff00ff00ff00fful);
				l = (l & 0x0000ffff0000fffful) + ((l >> 16) & 0x0000ffff0000fffful);
				l = (l & 0x00000000fffffffful) + ((l >> 32) & 0x00000000fffffffful);
			# else
				l = (l & 0x55555555ul) + ((l >> 1) & 0x55555555ul);
				l = (l & 0x33333333ul) + ((l >> 2) & 0x33333333ul);
				l = (l & 0x0f0f0f0ful) + ((l >> 4) & 0x0f0f0f0ful);
				l = (l & 0x00ff00fful) + ((l >> 8) & 0x00ff00fful);
				l = (l & 0x0000fffful) + ((l >> 16) & 0x0000fffful);
			# endif
		// #endif
		count += l;
	}
	return count;
}

#endif 

#include <omp.h>

int stark_shed_init() {
	CPU_ZERO(&allcpus);
	sched_getaffinity(0, sizeof(cpu_set_t), &allcpus);
	initialized = 1;
	
	int num_cpus = CPU_COUNT(&allcpus);
	
	fprintf(stderr, "stark_shed_init initialized\n");
	fprintf(stderr, "%d cpus available for scheduling\nAvailable CPUs: ", num_cpus);
	
	size_t i;
	for (i = 0; i < CPU_SETSIZE; i++) {
		if (CPU_ISSET(i, &allcpus))
			fprintf(stderr, "%zu, ", i);
	}
	fprintf(stderr, "\n");
	
	return num_cpus;
}

size_t get_num_available_cpus() {
	if (!initialized)
		stark_shed_init();
	return CPU_COUNT(&allcpus);
}

// numa_run_on_node

int stark_autopin_to_cpu(int threadnum) {
	size_t i;
	for (i = 0; i < CPU_SETSIZE; i++) {
		if (CPU_ISSET(i, &allcpus)) {
			if (threadnum)
				threadnum--;
			else {
				cpu_set_t cpu_set;
				CPU_ZERO(&cpu_set);
				CPU_SET(i, &cpu_set);
				return pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpu_set);
			}
		}
	}
	
	return -1;
}

int stark_run_on_cpu(int cpu) {
	if (!initialized)
		stark_shed_init();
	
	if (CPU_ISSET(cpu, &allcpus)) {
		cpu_set_t cpu_set;
		CPU_ZERO(&cpu_set);
		CPU_SET(cpu, &cpu_set);
		return pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpu_set);
	}
	
	return -1;
}

// int numa_node_of_cpu(int cpu);

#ifdef USE_LIBNUMA

#if !defined(SYS_getcpu) && defined(__NR_get_cpu)
#define SYS_getcpu __NR_get_cpu
#endif

// has to be called from master/main thread
void stark_numa_lock_master_thread() {
	int cpu, node, status;
	#ifdef SYS_getcpu
	
	status = syscall(SYS_getcpu, &cpu, &node, NULL);
	
	if (status == -1) 
	#endif
		return; // error
	
		// master_node = numa_node_of_cpu(curcpu);
	numa_run_on_node(master_node = node);
	
	// struct bitmask * bindmask = numa_bitmask_alloc(numa_num_possible_nodes());
	// numa_bitmask_clearall(bindmask);
	// numa_bitmask_setbit(bindmask, master_node);
	// copy_bitmask_to_nodemask(bindmask, &master_nodemask);
	// numa_bitmask_free(bindmask);
}

#endif

