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


#define stark_alloc_pages(num_pages, flags) calloc(num_pages, sysconf(_SC_PAGE_SIZE))
#define stark_free_pages(addr, num_pages) free(addr)

static inline void stark_malloc_starknode_depths(stark_t *stark) {
	if (!(stark->level[1] = calloc(3,sizeof(starknode_t)))) {
		CRITICAL_ERROR("Out of Memory!");
	}
	stark->memory_allocated[1] = 3 * sizeof(starknode_t);
}

static inline void stark_free_starknode_depth_array(stark_t *stark, depth_t depth) {
	if(stark->level[depth]) {
		free(stark->level[depth]);
		stark->level[depth] = NULL;
		stark->memory_allocated[depth] = 0;
		stark->maxsize[depth] = 0;
		stark->size[depth] = 0;
	}
}

static inline void stark_shrink_starknode_depth_array(stark_t *stark, depth_t depth) {

	if (stark->size[depth]) {
		void* tmp = realloc(stark->level[depth], (stark->size[depth] + 1) * sizeof(starknode_t));
		if (tmp) {
			stark->level[depth] = tmp;
			stark->maxsize[depth] = stark->size[depth] + 1;
			stark->memory_allocated[depth] = (stark->size[depth] + 1) * sizeof(starknode_t);
		}
	}
	else
		stark_free_starknode_depth_array(stark, depth);
}

static inline offset_t stark_create_node(stark_t* stark, depth_t depth) {
	if (stark->maxsize[depth] <= stark->size[depth] +1 ) {
		starknode_t* newstark;
		intptr_t sizeincrease;
		//long long maxneededsize = (1 << ((2 * level) - 1)) + 1;
		
		if (stark->maxsize[depth] * sizeof(starknode_t) < MAXALLOC)
			sizeincrease = stark->maxsize[depth] / sizeof(starknode_t);
		else
			sizeincrease = MAXALLOC / sizeof(starknode_t);
			
		if (!sizeincrease)
			sizeincrease = 2;
		
		while (sizeincrease >= (MINALLOC / sizeof(starknode_t)) && !(newstark = realloc(stark->level[depth],sizeof(starknode_t) * (stark->maxsize[depth] + sizeincrease))))
			sizeincrease >>= 1;
			
		if (!newstark) {
			//fprintf(stderr,"Allocation of %ld MB failed. Out of Memroy!\n",sizeincrease >> 20);
			fprintf(stderr,"Allocation of %ld bytes failed. Out of Memroy!\n",sizeincrease);
			exit(-1);
		} else {
			memset(newstark + stark->maxsize[depth], 0, sizeincrease * sizeof(starknode_t));
			stark->level[depth] = newstark;
			stark->maxsize[depth] += sizeincrease;
			stark->memory_allocated[depth] = sizeof(starknode_t) * (stark->maxsize[depth] + sizeincrease);
		}
		
	}
	//stark->level[stark->size].order = stark->size;
	return ++(stark->size[depth]);
}

#define stark_move_starknode_depth_array_to_numa_master(stark, depth)

