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


// #include "stark_alloc.h"
#include <unistd.h>
#include <sys/types.h>
#include <sys/mman.h>
#include "stark.h"
#include "stark_numa.h"
#include "main.h"

// sysconf(_SC_PAGE_SIZE)

#ifndef MAP_ANONYMOUS
#define MAP_ANONYMOUS MAP_ANON
#endif

#ifndef MAP_NORESERVE
#define MAP_NORESERVE 0
#endif

#define stark_alloc_pages(num_pages, flags) ({\
	void *tmp = mmap(NULL, (num_pages) * sysconf(_SC_PAGE_SIZE), PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE | (flags), -1, 0);\
	if (tmp == MAP_FAILED) { perror(__FILE__ ":" TOSTRING(__LINE__) ": " "mmap"); fflush(stderr); exit(-1); }\
	tmp; })

#define stark_free_pages(addr, num_pages) munmap((addr), (num_pages) * sysconf(_SC_PAGE_SIZE))
;

static inline void stark_malloc_starknode_depth_array(stark_t *stark, depth_t depth) {
	size_t allocated_kmers = 0;
	if (sizeof(offset_t) < 6) {
		
		#ifdef MAXMAP
		allocated_kmers = MAXMAP;
		#else
		allocated_kmers = OFFSET_MAX; // maximum adressable
		#endif
		
		if (depth < 30) {
			size_t max_kmers_possible = (UINT64_C(1) << (depth * 2)) + 1;
			if (max_kmers_possible < allocated_kmers)
				allocated_kmers = max_kmers_possible;
		}
		
		size_t memory_allocating  = allocated_kmers * sizeof(starknode_t);
		
		#ifdef DEBUG_MEMORY_CAP
		if (memory_allocating > DEBUG_MEMORY_CAP)
			memory_allocating = DEBUG_MEMORY_CAP;
		#endif
		
		// increase to page size;
		size_t page_mod = memory_allocating % sysconf(_SC_PAGE_SIZE);
		if (page_mod)
			memory_allocating += sysconf(_SC_PAGE_SIZE) - page_mod;
		
		stark->maxsize[depth] = (size_t)(memory_allocating / sizeof(starknode_t));
		void *tmp = mmap(NULL, memory_allocating, PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE | MAP_NORESERVE, -1, 0);
		if (tmp == MAP_FAILED) {
			perror(__FILE__ ":" TOSTRING(__LINE__) ": " "mmap");
			DEBUG_MSG("stark_malloc_starknode_depth_array: depth = %d, memory_allocating = %zu", depth, memory_allocating);
			fflush(stderr);
			exit(-1);
		}
		
		// fprintf(stderr, "Allocated %zu bytes for %zu %d-mers (%zu bytes each) at address %p\n", memory_allocating, allocated_kmers, depth, sizeof(starknode_t), tmp);
		
		stark->level[depth] = tmp;
		stark->memory_allocated[depth] = memory_allocating;
		
	} else {
		CRITICAL_ERROR("Huston, we have a memory allocation Problem!\nLarge Offsets are not supported in MMAP mode.\nNeigther are extravagant depths.");
	}
}

static inline void stark_malloc_starknode_depths(stark_t *stark) {
	depth_t depth;
	for (depth = 1; depth < HARD_MAXDEPTH; depth++) {
		stark_malloc_starknode_depth_array(stark, depth);
	}
	
	// exit(0);
}

static inline void stark_free_starknode_depth_array(stark_t *stark, depth_t depth) {
	if(munmap(stark->level[depth], stark->memory_allocated[depth])) {
		perror(__FILE__ ":" TOSTRING(__LINE__) ": " "munmap");
		fflush(stderr);
	}
	else {
		stark->level[depth] = NULL;
		stark->memory_allocated[depth] = 0;
		stark->maxsize[depth] = 0;
	}
}

static inline void stark_shrink_starknode_depth_array(stark_t *stark, depth_t depth) {
	size_t used_memory = (stark->size[depth] + 1) * sizeof(starknode_t);
	// increase to page size;
	size_t page_mod = used_memory % sysconf(_SC_PAGE_SIZE);
	if (page_mod)
		used_memory += sysconf(_SC_PAGE_SIZE) - page_mod;
	
	;
	
	if(stark->memory_allocated[depth] - used_memory && mmap(((char*)(stark->level[depth])) + used_memory, stark->memory_allocated[depth] - used_memory, PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE | MAP_NORESERVE | MAP_FIXED, -1, 0) == MAP_FAILED) {
	// if(madvise(((char*)(stark->level[depth])) + used_memory, stark->memory_allocated[depth] - used_memory, MADV_DONTNEED)) {
		perror(__FILE__ ":" TOSTRING(__LINE__) ": " "madvise");
		fflush(stderr);
	}
}

/*
#ifdef USE_LIBNUMA
// doesn't work at the moment
static inline void stark_move_starknode_depth_array_to_numa_master(stark_t *stark, depth_t depth) {
	size_t used_memory = (stark->size[depth] + 1) * sizeof(starknode_t);
	// increase to page size;
	size_t page_mod = used_memory % sysconf(_SC_PAGE_SIZE);
	if (page_mod)
		used_memory += sysconf(_SC_PAGE_SIZE) - page_mod;
	
	stark_numa_migrate_pages_to_master_node(stark->level[depth], used_memory);
}


static inline void stark_move_all_to_master_memory(stark_t* const stark) {
	depth_t depth;
	
	for (depth = 1; stark->level[depth] && depth < HARD_MAXDEPTH; depth++) {
		stark_move_starknode_depth_array_to_numa_master(stark, depth);
	}
	
}

#else */
// nop
#define stark_move_starknode_depth_array_to_numa_master(stark, depth)
#define stark_move_all_to_master_memory(stark)
// #endif

static inline offset_t stark_create_node(stark_t* stark, depth_t depth) {
	offset_t offset = __sync_add_and_fetch(stark->size + depth, 1);
	
	if (offset <= 0 || offset == OFFSET_MAX) { // offset overflow, memory error
		CRITICAL_ERROR("offset_t overflow, you have too many k-mers in your job. recompile with a larger offset_t.");
	}
	
	// memset(stark->level[depth] + offset, 0, sizeof(stark->level[depth][offset]));
	
	// #define MADVISE_WILLNEED_CHUNKSIZE 0x400000
	// 	starknode_t * node = stark->level[depth] + offset;
	// 	if (
	// 		(
	// 		(uintptr_t)(node)
	// 		^
	// 		(uintptr_t)(node - 1) 
	// 		)
	// 		& MADVISE_WILLNEED_CHUNKSIZE
	// 	) {
	// 		madvise(
	// 		(void *)(((uintptr_t)(node) & ~(MADVISE_WILLNEED_CHUNKSIZE - 1)) + MADVISE_WILLNEED_CHUNKSIZE)
	// 		, MADVISE_WILLNEED_CHUNKSIZE
	// 		, MADV_WILLNEED
	// 		);
	// 	}
	
	// #define PAGE_PREFETCH_CHUNKSIZE 0x1000
	// starknode_t * node = stark->level[depth] + offset;
	// char * node_nextpage = (int *)node + 0x1000;
	// if (
	// 	(
	// 	(uintptr_t)(node)
	// 	^
	// 	(uintptr_t)(node - 1) 
	// 	)
	// 	& PAGE_PREFETCH_CHUNKSIZE
	// ) {
	// 	*node_nextpage = 0;
	// }
	
	// __builtin_prefetch((void *)((uintptr_t)((char*)(stark->level[depth] + offset) + 0x1000) & ~0xFFF), 1, 0);
	
	return offset;
}


