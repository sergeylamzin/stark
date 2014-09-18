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


// this header should never be loaded directly, but only through stark.h
// this whole thing is not supported on MACH (OS X) targets
#ifndef STARK_INCREMENT_CACHE_H
#define STARK_INCREMENT_CACHE_H
#if defined(AMO_CACHE)

// #define _GNU_SOURCE
#include <unistd.h>
#include <stdint.h>
#include "stark_alloc.h"
#include <sys/types.h>
#include <sys/mman.h>
#include "stark_status.h"
#include "list.h"

#ifndef MAP_ANONYMOUS
#define MAP_ANONYMOUS MAP_ANON
#endif

#ifdef UV_SYSTEM
#include "uv_system.h"
#endif


union stark_coverage_counter_cache_u {
	struct {
		char boffset[5];
		uint8_t counter;
		uint16_t bdepth;
	} b;
	uint64_t q;
};

#ifndef STARK_COVERAGE_CACHE_L1_MEMORY_PAGES
#define STARK_COVERAGE_CACHE_L1_MEMORY_PAGES 1024
#endif

#define STARK_COVERAGE_CACHE_L2_MEMORY 0x400000
#define STARK_COVERAGE_CACHE_L2_MODULUS(value) (((STARK_COVERAGE_CACHE_L2_MEMORY >> 3)-1) & value)


struct stark_coverage_mytoken_s {
	uintmax_t mytoken;
	// coverage_t ** cache;
	// 	size_t index;
	// 	size_t maxsize;
	list_t(coverage_t *) cache;
};


#ifdef AMO_CACHE_MAIN_C

coverage_t * stark_coverage_counter_cache_full[32];
depth_t stark_coverage_counter_cache_full_maxdepth = 0;
void * stark_coverage_counter_page;

#ifdef STARK_COVERAGE_CACHE_TOKEN
uintmax_t stark_coverage_token_mask = 0;
volatile uintmax_t stark_coverage_token[0x10] = {[0 ... 0xF] = 1};
volatile uintmax_t stark_coverage_token_round[0x10] = {0};
struct stark_coverage_mytoken_s stark_coverage_mytoken[0x10] = {0};
#endif

#else

extern coverage_t * stark_coverage_counter_cache_full[32];
extern depth_t stark_coverage_counter_cache_full_maxdepth;
extern void * stark_coverage_counter_page;

#ifdef STARK_COVERAGE_CACHE_TOKEN
extern uintmax_t stark_coverage_token_mask;
extern volatile uintmax_t stark_coverage_token[0x10];
extern volatile uintmax_t stark_coverage_token_round[0x10];
extern struct stark_coverage_mytoken_s stark_coverage_mytoken[0x10];
#endif

#endif

#pragma omp threadprivate(stark_coverage_counter_cache_full, stark_coverage_counter_page)

// AMO_CACHE_MAIN_C coverage_t * stark_coverage_counter_cache_full[32];
// AMO_CACHE_MAIN_C depth_t stark_coverage_counter_cache_full_maxdepth;
// AMO_CACHE_MAIN_C void * stark_coverage_counter_page;


#ifdef STARK_COVERAGE_CACHE_TOKEN


#pragma omp threadprivate(stark_coverage_counter_cache_full_maxdepth, stark_coverage_mytoken)

// AMO_CACHE_MAIN_C uintmax_t stark_coverage_token_mask;
// AMO_CACHE_MAIN_C volatile uintmax_t stark_coverage_token[0x10];
// AMO_CACHE_MAIN_C volatile uintmax_t stark_coverage_token_round[0x10];

// AMO_CACHE_MAIN_C struct stark_coverage_mytoken_s stark_coverage_mytoken[0x10];


static inline void stark_coverage_token_shift(int_fast32_t token_index) {
	if (!stark_coverage_token_mask)
		return;
	
	volatile uintmax_t * token_addr = stark_coverage_token + token_index;
	
	register uintmax_t token = *token_addr;
	do {
		#if defined(__x86_64__)
		asm (
			  "rol %1, %0"
			: "+r" (token)
			: "I"  (1)
		);
		#else
		if (!(token <<= 1))
			token = 1;
		#endif
	} while (!(token & stark_coverage_token_mask));
	
	stark_coverage_token_round[token_index]++;
	__sync_synchronize();
	*token_addr = token;
}

static inline int stark_coverage_token_init_tryset(const uintmax_t mytoken) {
	const uintmax_t lastmask = __sync_fetch_and_or(&stark_coverage_token_mask, mytoken);
	if (!(~lastmask))
		return -1; // mask is full
	if (lastmask & mytoken)
		return 0;
	
	return 1;
}

static inline void stark_coverage_token_init() {
	if (!~stark_coverage_token_mask) {
		fprintf(stderr, "too many threads, currently maxmimum supported thread with STARK_COVERAGE_CACHE_TOKEN is %zu\n", sizeof(stark_coverage_token_mask) * 8);
		exit(-1);
		return;
	}
	
	stark_coverage_mytoken->mytoken = 1;
	
	int status;
	size_t i;
	for (i = 0; i < (sizeof(stark_coverage_token_mask) * 8); i++) {
		if (!(status = stark_coverage_token_init_tryset(stark_coverage_mytoken->mytoken))) {
			stark_coverage_mytoken->mytoken <<= 1;
		} else if (status < 0) {
			stark_coverage_mytoken->mytoken = 0;
			break;
		} else
			break;
	}
	
	for (i = 0; i < (sizeof(stark_coverage_mytoken) / sizeof(*stark_coverage_mytoken)); i++) {
		list_init_size(&stark_coverage_mytoken[i].cache, 0x1000);
	}
	
	#if defined(STARK_STATUS_H)
	if (get_stark_thread_status())
		get_stark_thread_status()->mytoken = stark_coverage_mytoken;
	#endif
}


static inline void stark_increment_coverage_cached_tokenized_free(size_t i) {
	list_free(stark_coverage_mytoken[i].cache);
}

static inline void stark_coverage_token_release() {
	__sync_fetch_and_and(&stark_coverage_token_mask, ~stark_coverage_mytoken->mytoken);
	size_t i;
	for (i = 0; i < (sizeof(stark_coverage_token) / sizeof(*stark_coverage_token)); i++) {
		if (stark_coverage_token[i] == stark_coverage_mytoken->mytoken) {
			// if (stark_coverage_mytoken[i].cache.list)
				// stark_increment_coverage_cached_tokenized_free(i);
			stark_coverage_token_shift(i);
		}
	}
	
	for (i = 0; i < (sizeof(stark_coverage_mytoken) / sizeof(*stark_coverage_mytoken)); i++) {
		list_free(stark_coverage_mytoken[i].cache);
	}
	
	#if defined(STARK_STATUS_H)
	get_stark_thread_status()->mytoken = NULL;
	#endif
}

static inline void stark_increment_coverage_cached_token_wait(volatile uintmax_t * token_addr ) {
	while (*token_addr != stark_coverage_mytoken->mytoken);
}

static inline void stark_increment_coverage_cached_tokenized_tryflush() {
	size_t i;
	for (i = 0; i < (sizeof(stark_coverage_token) / sizeof(*stark_coverage_token)); i++) {
		if (stark_coverage_token[i] == stark_coverage_mytoken->mytoken) {
			stark_set_status_action("Got token, flushing cache.");
			size_t cache_size = stark_coverage_mytoken[i].cache.size;
			while (cache_size) {
				coverage_t * address = stark_coverage_mytoken[i].cache.list[--cache_size];
				*address += 1;
			}
			stark_coverage_mytoken[i].cache.size = cache_size;
			stark_coverage_token_shift(i);
		}
	}
	
}

static inline void stark_increment_coverage_cached_tokenized_waitflush() {
	size_t i, total = 0;
	for (i = 0; i < (sizeof(stark_coverage_token) / sizeof(*stark_coverage_token)); i++) {
		total += stark_coverage_mytoken[i].cache.size;
	}
	
	while (total) {
		for (i = 0; i < (sizeof(stark_coverage_token) / sizeof(*stark_coverage_token)); i++) {
			stark_set_status_action("Waiting for token.");
			if (stark_coverage_token[i] == stark_coverage_mytoken->mytoken) {
				stark_set_status_action("Got token, flushing cache.");
				size_t cache_size = stark_coverage_mytoken[i].cache.size;
				while (cache_size) {
					coverage_t * address = stark_coverage_mytoken[i].cache.list[--cache_size];
					*address += 1;
					total--;
				}
				stark_coverage_mytoken[i].cache.size = cache_size;
				stark_coverage_token_shift(i);
				stark_set_status_action("Waiting for token.");
			}
		}
	}
	
	
}

static inline void stark_increment_coverage_cached_tokenized(coverage_t * address) {
	size_t bucket_num = ((uintptr_t)address >> 12) & 0xF;
	
	list_insert(&(stark_coverage_mytoken[bucket_num].cache), address);
	
	// if ((stark_coverage_mytoken[bucket_num].index * sizeof(*stark_coverage_mytoken->cache)) >= stark_coverage_mytoken[bucket_num].maxsize) {
	// 	if (stark_coverage_mytoken[bucket_num].cache) {
	// 		#ifdef MREMAP_MAYMOVE
	// 		// #warning "Notice: Using Linux specific mremap call."
	// 		void * tmp = mremap(stark_coverage_mytoken[bucket_num].cache
	// 							, stark_coverage_mytoken[bucket_num].maxsize
	// 							, stark_coverage_mytoken[bucket_num].maxsize << 1
	// 		                    , MREMAP_MAYMOVE);
	// 		if (tmp == MAP_FAILED) {
	// 			perror("mremap");
	// 			exit(0);
	// 		}
	// 		
	// 		#else
	// 		void * tmp = mmap(NULL
	// 						, stark_coverage_mytoken[bucket_num].maxsize << 1
	// 						, PROT_READ | PROT_WRITE
	// 						, MAP_PRIVATE | MAP_ANONYMOUS
	// 						, -1
	// 						, 0);
	// 		if (tmp == MAP_FAILED) {
	// 			perror("mmap");
	// 			exit(0);
	// 		}
	// 		
	// 		memcpy(tmp, stark_coverage_mytoken[bucket_num].cache, stark_coverage_mytoken[bucket_num].maxsize);
	// 		stark_increment_coverage_cached_tokenized_free(bucket_num);
	// 		
	// 		#endif
	// 		
	// 		stark_coverage_mytoken[bucket_num].cache = tmp;
	// 		stark_coverage_mytoken[bucket_num].maxsize <<= 1;
	// 		
	// 	} else {
	// 		stark_coverage_mytoken[bucket_num].cache = mmap(NULL
	// 										, sysconf(_SC_PAGE_SIZE)
	// 										, PROT_READ | PROT_WRITE
	// 										, MAP_PRIVATE | MAP_ANONYMOUS
	// 										, -1
	// 										, 0);
	// 		if (stark_coverage_mytoken[bucket_num].cache == MAP_FAILED) {
	// 			perror("mmap");
	// 			exit(0);
	// 		}
	// 		
	// 		stark_coverage_mytoken[bucket_num].maxsize = sysconf(_SC_PAGE_SIZE);
	// 	}
	// }
	// 
	// stark_coverage_mytoken[bucket_num].cache[stark_coverage_mytoken[bucket_num].index++] = address;
	
}

#endif

// #ifdef AMO_CACHE_L2
// extern __thread union stark_coverage_counter_cache_u * stark_coverage_counter_cache_L2;
// #endif
// 
static inline void stark_coverage_counter_cache_init() {
	if (stark_coverage_counter_page)
		return;
	
	coverage_t * page = stark_alloc_pages(STARK_COVERAGE_CACHE_L1_MEMORY_PAGES, 0);
	stark_coverage_counter_page = page;
	
	// #ifdef AMO_CACHE_L2
	// stark_coverage_counter_cache_L2 = mmap(NULL, STARK_COVERAGE_CACHE_L2_MEMORY, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
	// if (stark_coverage_counter_cache_L2 == MAP_FAILED) {
	// 	perror("mmap");
	// 	exit(-1);
	// }
	// #endif
	
	ssize_t num_counters = (STARK_COVERAGE_CACHE_L1_MEMORY_PAGES * sysconf(_SC_PAGE_SIZE)) / sizeof(coverage_t);
	
	size_t i;
	for (i = 1; i < 32; i++) {
		num_counters -= k_mer_universe_size(i) +1;
		if (num_counters >= 0) {
			stark_coverage_counter_cache_full[i] = page;
			page += k_mer_universe_size(i) +1;
			// DEBUG_MSG("Allocated counter cache for level %zu, %zu+1 counters for thread 0x%lu", i, k_mer_universe_size(i), pthread_self());
		}
		else
			break;
	}
	
	stark_coverage_counter_cache_full_maxdepth = i-1;
	
	// for (; i < 16; i++)
	// 	stark_coverage_counter_cache[i] = NULL;
	
}

#ifdef UV_SYSTEM
	#define stark_increment_coverage_real(addr, increment, depth) (depth < 5 ? uv_gru_atomic_async(addr, increment, GAMER_AMO_ADD_32_OPWORD) : __sync_fetch_and_add(addr, increment))
#else
	#define stark_increment_coverage_real(addr, increment, depth) __sync_fetch_and_add(addr, increment)
#endif

static inline void stark_increment_coverage_cached_L2(stark_t * stark, coverage_t * addr, depth_t depth, offset_t offset) {
	
}

static inline void stark_increment_coverage_flush_L2_entry(stark_t * stark, union stark_coverage_counter_cache_u * cachefield_addr) {
	size_t old_depth = cachefield_addr->b.bdepth;
	size_t old_offset = 0xFFFFFFFFFFUL & cachefield_addr->q;
	__sync_fetch_and_add(&(stark->level[old_depth][old_offset].coverage), cachefield_addr->b.counter);
	// cachefield_addr->b.counter = 0;
}



static inline void stark_increment_coverage_cached(stark_t * stark, coverage_t * addr, depth_t depth, offset_t offset) {
	if (depth > stark_coverage_counter_cache_full_maxdepth) {
		// #ifdef AMO_CACHE_L2
		// // go to L2
		// union stark_coverage_counter_cache_u * cachefield_addr = stark_coverage_counter_cache_L2 + (STARK_COVERAGE_CACHE_L2_MODULUS(depth * offset));
		// uint64_t cachefield_id = ((uint64_t)offset & 0xFFFFFFFFFFUL) | (((uint64_t)depth & 0xFFFF) << 48);
		// if ((cachefield_addr->q & 0xFFFF00FFFFFFFFFFUL) != cachefield_id) {
		// 	// cachefield in use
		// 	// flush out the old value
		// 	stark_increment_coverage_flush_L2_entry(stark, cachefield_addr);
		// 	cachefield_addr->q = cachefield_id | 0x10000000000;
		// 	// cachefield_addr->b.counter = 1;
		// } else if (!(cachefield_addr->b.counter += 1)) {
		// 	__sync_fetch_and_add(&(stark->level[depth][offset].coverage), (coverage_t)0x100);
		// }
		// #el
		#ifdef STARK_COVERAGE_CACHE_TOKEN
		stark_increment_coverage_cached_tokenized(addr);
		#elif UV_SYSTEM
		stark_uv_gru_increment_coverage_w_fallback(addr, 1);
		#else
		__sync_fetch_and_add(addr, (coverage_t)1);
		#endif
	}
	else if (!(++stark_coverage_counter_cache_full[depth][offset])) {
		fprintf(stderr, "Coverage counter overflow on depth %d, offset %lu", depth, (uint64_t)offset);
	}
}

void stark_coverage_counter_cache_flush_free(stark_t * stark);

#else

void nop();
#define stark_coverage_counter_cache_init nop
#define stark_coverage_counter_cache_flush_free nop
#define stark_increment_coverage_cached(stark, addr, depth, offset) __sync_fetch_and_add(addr, (coverage_t)1)

#endif
#endif
