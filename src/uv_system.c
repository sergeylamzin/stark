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


#ifdef UV_SYSTEM
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include "uv_system.h"
#include <uv/gru/gru_alloc.h>
#include "main.h"

// __thread gru_cookie_t gru_cookie[UV_GRU_CONTEXTS_PER_THREAD];
__thread gru_control_block_t *gru_control_block[UV_GRU_MAX_CBS_PER_THREAD];
// static __thread gru_segment_t* gru_gruseg[UV_GRU_CONTEXTS_PER_THREAD];


static __thread int uv_gru_initialized = 0;


size_t gru_cbs_per_cpu = 0;
__thread size_t uv_gru_num_context = 0;

__thread size_t times_waited = 0;
__thread size_t total_cycles_waited = 0;
__thread size_t times_called;
__thread size_t total_cycles_between;
__thread size_t last_call;

int uv_init_global() {
	gru_alloc_thdata_t thdata;
	if (gru_temp_reserve_try(&thdata)) {
		perror("gru_temp_reserve_try");
		return -1;
	}
	
	gru_temp_release();
	
	return (gru_cbs_per_cpu = thdata.cb_cnt);
}

int uv_init_atomic() {
	if (uv_gru_initialized)
		return 0;
	
	size_t gru_cbs_per_thread = gru_cbs_per_cpu >> 3;
	if (!gru_cbs_per_thread)
		return -1;
	
	if (gru_cbs_per_thread > UV_GRU_MAX_CBS_PER_THREAD)
		gru_cbs_per_thread = UV_GRU_MAX_CBS_PER_THREAD;
	
	
	
	// for (num_context = 0; num_context < UV_GRU_CONTEXTS_PER_THREAD; num_context++) {
		void* dseg;
		gru_segment_t* gseg;
		int cbnum;
		if(gru_pallocate(gru_cbs_per_thread, 0, &gseg, &cbnum, &dseg)) {
			perror("gru_pallocate");
			return -1;
		}
		
		for (uv_gru_num_context = 0; uv_gru_num_context < gru_cbs_per_thread; uv_gru_num_context++) {
			gru_control_block[uv_gru_num_context] = gru_get_cb_pointer(gseg, cbnum + uv_gru_num_context);
		}
		
		/*
		if (gru_create_context(gru_cookie + num_context, NULL,
		             1, 0,
		             1, GRU_OPT_MISS_DEFAULT)) {
			perror("gru_create_context");
			return -1;
		}
		
		gru_segment_t *gseg;
		gseg = gru_get_thread_gru_segment(gru_cookie[num_context], 0);
		if (!gseg) {
			perror("gru_get_thread_gru_segment");
			return -1;
		}
		
		gru_control_block[num_context] = gru_get_cb_pointer(gseg, 0);
		*/
	// }
	
	last_call = rdtsc();
	
	uv_gru_initialized = 1;
	
	return 0;
}

void uv_free_atomic() {
	DEBUG_MSG("Freeing GRU resources.");
	
	if (!uv_gru_initialized)
		return;
	
	gru_alloc_reset_state();
	
	if (times_called) {
		DEBUG_MSG("GRU AMO (add dword) on thread 0x%lX was called %zu times with %zu cycles on average between calls.", pthread_self(), times_called, total_cycles_between / times_called);
	}
	
	if (times_waited) {
		DEBUG_MSG("Thread 0x%lX called gru_wait_abort() a total of %zu times taking %zu cycles on average.", pthread_self(), times_waited, total_cycles_waited / times_waited);
	}
	else {
		DEBUG_MSG("Thread 0x%lX never waited for the GRU.", pthread_self());
	}
		
}




#endif // UV_SYSTEM
