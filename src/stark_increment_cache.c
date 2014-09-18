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


#if defined(AMO_CACHE)

#define AMO_CACHE_MAIN_C

#include "main.h"
#include "stark.h"

#include "stark_alloc.h"
#include "stark_increment_cache.h"

/*
void stark_coverage_counter_cache_free() {
	if (stark_coverage_counter_page)
		stark_free_pages(stark_coverage_counter_page, 1);
}

void stark_coverage_counter_cache_flush(stark_t * stark) {
	if (!stark_coverage_counter_page)
		return;
	
	depth_t depth;
	offset_t offset;
	
	for(depth = 1; depth < 16; depth++) {
		if (stark_coverage_counter_cache[depth]) {
			for(offset = 0; offset < k_mer_universe_size(i) +1; offset++) {
				coverage_t c;
				if ((c = stark_coverage_counter_cache[depth][offset])) {
					stark->level[depth][offset].coverage += c;
					stark_coverage_counter_cache[depth][offset] = 0;
				}
			}
		} else
			break;
	}
}
*/

// #ifdef STARK_COVERAGE_CACHE_TOKEN
// 
// uintmax_t stark_coverage_token_mask = 0;
// volatile uintmax_t stark_coverage_token[0x10] = {[0 ... 0xF] = 1};
// volatile uintmax_t stark_coverage_token_round[0x10] = {[0 ... 0xF] = 0};
// struct stark_coverage_mytoken_s stark_coverage_mytoken[0x10] = {0};
// 
// #endif
// 
// coverage_t * stark_coverage_counter_cache_full[32];
// depth_t stark_coverage_counter_cache_full_maxdepth = 0;
// void * stark_coverage_counter_page;
// __thread union stark_coverage_counter_cache_u * stark_coverage_counter_cache_L2;

// #pragma omp threadprivate(stark_coverage_counter_cache_full, stark_coverage_counter_cache_full_maxdepth, stark_coverage_counter_page, stark_coverage_mytoken)

void stark_coverage_counter_cache_flush_free(stark_t * stark) {
	if (!stark_coverage_counter_page)
		return;
	
	depth_t depth;
	offset_t offset;
	
	for(depth = 1; depth <= stark_coverage_counter_cache_full_maxdepth; depth++) {
		// if (stark_coverage_counter_cache[depth]) {
		for(offset = 0; offset < k_mer_universe_size(depth) +1; offset++) {
			coverage_t c;
			if ((c = stark_coverage_counter_cache_full[depth][offset])) {
				__sync_fetch_and_add(&(stark->level[depth][offset].coverage), c);
				stark_coverage_counter_cache_full[depth][offset] = 0;
			}
		}
		// } else
		// 	break;
	}
	
	stark_free_pages(stark_coverage_counter_page, STARK_COVERAGE_CACHE_L1_MEMORY_PAGES);
	
	// #ifdef AMO_CACHE_L2
	// // DEBUG_MSG("Clearing thread 0x%lX L2 AMO cache.", pthread_self());
	// 
	// size_t i = STARK_COVERAGE_CACHE_L2_MEMORY >> 3;
	// union stark_coverage_counter_cache_u * cachefield_addr = stark_coverage_counter_cache_L2;
	// do {
	// 	// union stark_coverage_counter_cache_u * cachefield_addr = stark_coverage_counter_cache_L2 + (--i);
	// 	if (cachefield_addr->b.counter) {
	// 		stark_increment_coverage_flush_L2_entry(stark, cachefield_addr);
	// 		cachefield_addr->q = 0;
	// 	}
	// 	cachefield_addr++;
	// } while (--i);
	// munmap(stark_coverage_counter_cache_L2, STARK_COVERAGE_CACHE_L2_MEMORY);
	// #endif
	
}

#else

void nop() {}

#endif
