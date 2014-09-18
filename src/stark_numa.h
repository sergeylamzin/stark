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



#ifndef STARK_NUMA_H
#define STARK_NUMA_H

// #define _GNU_SOURCE
#include <sched.h>

#ifdef __CPU_SETSIZE
int stark_shed_init();
size_t get_num_available_cpus();
int stark_autopin_to_cpu(int threadnum);
#else
#define stark_shed_init(a) 0
#define get_num_available_cpus() 1
#define stark_autopin_to_cpu(a) 0
#endif

#ifdef USE_LIBNUMA
#include <numa.h>
#include <numaif.h>

extern int master_node;
extern nodemask_t master_nodemask;

// has to be called from master/main thread
void stark_numa_lock_master_thread();

static inline void stark_numa_migrate_pages_to_master_node(void * addr, unsigned long len) {
	if (master_node < 0)
		return;
	
	if ( mbind(	  addr
				, len
				, MPOL_BIND
	            , master_nodemask.n
				, numa_max_node()
				, MPOL_MF_MOVE)) {
		perror("mbind");
		fprintf(stderr, "stark_numa_migrate_pages_to_master_node(%p, %lu)\n", addr, len);
	}
	
}
#else

// make them noops
#define stark_numa_lock_master_thread()
#define stark_numa_migrate_pages_to_master_node(addr, len)

#endif

#endif
