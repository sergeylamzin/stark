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
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifndef HASHMAP_H
#define HASHMAP_H

#include <sys/types.h>
#include <sys/mman.h>
#include <unistd.h>
#include <stdint.h>


#define hashmap_t(type) struct { \
	type * map; \
	ssize_t size; \
	size_t maxlength; \
	size_t modulus; \
	unsigned char * incore; \
	ssize_t incore_length; \
} 


#ifndef map_hashfunction
#define map_hashfunction(_object) _object
#endif


#define map_free(_map) (free((_map)->incore), munmap((_map)->map, (_map)->incore_length * sysconf(_SC_PAGESIZE)))

typedef hashmap_t(uint64_t) hashmap_uint64_t;

#define map_init(_map) map_init_size(_map, 0x200)
static inline void map_init_size( hashmap_uint64_t * const map, size_t initsize) {
	// if (objectsize & 0x3) {
	// 	fprintf(stderr, "map_init_size_alloc: Object size has to be a multiple of 4\n");
	// 	abort();
	// }
	
	map->size = 0;
	
	if (initsize & (initsize - 1)) {
		int bits = 1;
		while (initsize) {
			bits ++;
			initsize >>=1;
		}

		initsize = ((size_t)1) << bits;
	}
	
	map->modulus = initsize - 1;
	map->maxlength = initsize;
	
	size_t memsize = initsize * sizeof(map->map[0]);
	if (memsize & (sysconf(_SC_PAGESIZE) - 1)) {
		memsize += sysconf(_SC_PAGESIZE);
		memsize &= ~(sysconf(_SC_PAGESIZE) - 1);
	}
	
	void * tmp = mmap(
			  NULL
			, memsize
			, PROT_READ | PROT_WRITE
			, MAP_PRIVATE | MAP_ANON
			, -1
			, 0
		);
	if (tmp == MAP_FAILED) {
		perror("mmap");
		abort();
	}
	
	map->map = tmp;
	
	map->incore_length = memsize / sysconf(_SC_PAGESIZE);
	if (!(map->incore = calloc(map->incore_length, 1))) {
		perror("calloc");
		abort();
	}
}

static inline void map_shrink(hashmap_uint64_t * const map, size_t initsize) {
	if (initsize & (initsize - 1)) {
		int bits = 1;
		while (initsize) {
			bits ++;
			initsize >>=1;
		}

		initsize = ((size_t)1) << bits;
	}
	
	size_t memsize = initsize * sizeof(map->map[0]);
	if (memsize & (sysconf(_SC_PAGESIZE) - 1)) {
		memsize += sysconf(_SC_PAGESIZE);
		memsize &= ~(sysconf(_SC_PAGESIZE) - 1);
	}
	
	if (memsize < map->maxlength * sizeof(map->map[0])) {
		if (munmap(((char*)map->map) + memsize, (map->maxlength * sizeof(map->map[0])) - memsize)) {
			perror("munmap");
			abort();
		}
		
		map->modulus = initsize - 1;
		map->maxlength = initsize;
		
		map->incore_length = memsize / sysconf(_SC_PAGESIZE);
		if (!(map->incore = realloc(map->incore, map->incore_length))) {
			perror("realloc");
			abort();
		}
	}
	
}

static inline size_t map_hash_real(hashmap_uint64_t * const map, uint64_t object, size_t hash) {
		
	hash &= map->modulus;
	
	for (;;) {
		// if (!hash)
		// 	hash = 2305843009213693951 & map->modulus;
		
		if (!map->map[hash]
			|| (map->map[hash] == object)
		)
			return hash;
		
		hash = (hash + 1) & map->modulus;
	}
}

static inline int map_insert(hashmap_uint64_t * const map, uint64_t object) {
	
	if (__builtin_expect((map->size << 1) >= map->maxlength, 0)) {
		uint64_t * old = map->map;
		free(map->incore);
		
		size_t memsize_old = map->maxlength * sizeof(map->map[0]);
		if (memsize_old & (sysconf(_SC_PAGESIZE) - 1)) {
			memsize_old += sysconf(_SC_PAGESIZE);
			memsize_old &= ~(sysconf(_SC_PAGESIZE) - 1);
		}
		
		size_t memsize = 2 * map->maxlength * sizeof(map->map[0]);
		if (memsize & (sysconf(_SC_PAGESIZE) - 1)) {
			memsize += sysconf(_SC_PAGESIZE);
			memsize &= ~(sysconf(_SC_PAGESIZE) - 1);
		}
		void * tmp = mmap(
				  NULL
				, memsize
				, PROT_READ | PROT_WRITE
				, MAP_PRIVATE | MAP_ANON
				, -1
				, 0
			);
		if (tmp == MAP_FAILED) {
			perror("mmap");
			abort();
		}
		
		map->maxlength <<= 1;
		map->modulus = (map->modulus << 1) |  1;
		map->map = tmp;
		map->size = 0;
		
		map->incore_length = memsize / sysconf(_SC_PAGESIZE);
		if (!(map->incore = calloc(map->incore_length, 1))) {
			perror("calloc");
			abort();
		}
		
		size_t i;
		
		// rehash
		for (i = 0 ; i < (map->maxlength >> 1); i++) {
			if (old[i])
				map_insert(map, old[i]);
		}
		
		munmap(old, memsize_old);		
	}
	
	size_t realhash = map_hash_real(map, object, map_hashfunction(object));
	uint64_t * bucket = map->map + realhash;
	if (__builtin_expect(!*bucket, 1)) {
		*bucket = object;
		map->size++;
		map->incore[sizeof(object) * realhash / sysconf(_SC_PAGESIZE)] = 1;
		return 0;
	} else
		return 1;
}

static inline void map_zero(hashmap_uint64_t * const map) {
	char *array = (void *)map->map;
	size_t i;
	
	// #pragma omp for schedule(static, 0x100)
	for (i = 0; i < map->incore_length; i++) {
		if (map->incore[i])
			memset(array + (i * sysconf(_SC_PAGESIZE)), 0, sysconf(_SC_PAGESIZE));
	}
	
	// #pragma omp master
	// {
	memset(map->incore, 0, map->incore_length);
	map->size = 0;
	// }
}

// static inline size_t map_hashfunction_default(void * object, size_t objectsize) {
// 	size_t hash = 2305843009213693951;
// 	
// 	int i;
// 	for (i = 0; i < objectsize; i++) {
// 		hash ^= ((char *)object)[i];
// 		hash *= 2305843009213693951;
// 	}
// }

#define map_get(_map, object) ((_map)->map[map_hash_real(_map, object, (object))])

#ifdef HASHMAP_DEBUG_TEST

int main(int argc, char ** argv) {
	
	
	hashmap_uint64_t map;
	map_init_size(&map, 0x200);
	
	srand(0);
	
	size_t i;
	for (i = 0; i < atoi(argv[1]); i++) {
		int a = rand() % 0x10000;
		a++;
		if (map_get(&map, a)) {
			printf("Element %d already in map.\n", a);
		} else {
			map_insert(&map, a);
		}
	}
	
	
	printf("sysconf(_SC_PAGESIZE) = %zu\n", sysconf(_SC_PAGESIZE));
	
	map_zero(&map);
	
	for (i = 0 ; i < map.maxlength; i++) {
		if (map.map[i])
			printf("Zero didn't work : map.map[%zu] = %zu \n", i , map.map[i]);
	}
	
	map_free(&map);
	
	return 0;
}

#endif

#endif
