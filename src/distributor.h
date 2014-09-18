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


// #include <pthread.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <unistd.h>
#include <sys/mman.h>

#ifndef DISTRIBUTOR_H
#define DISTRIBUTOR_H

#ifndef MAP_ANONYMOUS
#define MAP_ANONYMOUS MAP_ANON
#endif

struct distributor_s {
	void* pages;
	size_t page_count;
	size_t chunk_size;
	size_t modulus;
	// #ifdef _POSIX_THREADS
	// pthread_mutex_t mutex;
	// #endif
	volatile size_t read_offset;
	volatile size_t write_offset;
	void * volatile free_pages[0];
};

typedef struct distributor_s page_distributor_t;

// #ifdef DISTRIBUTOR_TRUST
static inline void distributor_put(struct distributor_s * const __restrict distributor, void* const page ) {
	if (page) {
		distributor->free_pages[distributor->write_offset & distributor->modulus] = page;
		distributor->write_offset++;
	}
}

// #define STRINGIFY(x) #x
// #define TOSTRING(x) STRINGIFY(x)

/*
#if defined(__x86_64__) && defined(ASM_DOES_NOT_WORK_YET)
static inline void distributor_put_sync(struct distributor_s * const __restrict distributor, void* const page ) {
	if (page) {
		size_t write_offset;
		void * a;
		__asm__ (
			// "xor %%rdx, %%rdx"
			"1:\n\t"
			"mov %[dw], %[write_offset]\n\t"
			"mov %[offset](%[distributor], %[write_offset], 8), %[a]\n\t"
			"test %[a], %[a]\n\t"
			"jz 1b\n\t"
			"lock cmpxchg %[page], %[offset](%[distributor], %[write_offset], 8)\n\t"
			"jnz 1b\n\t"
			// "inc %[write_offset]\n\t"
			"inc %[write_offset]\n\t"
			"cmp %[page_count], %[write_offset]\n\t"
			"cmovg $0, %[write_offset]\n\t"
			"mov %[write_offset], %[dw]"
			: [write_offset] "+d" (write_offset)
			, [dw] "+m" (distributor->write_offset)
			, [addr] "=m" (distributor->free_pages[write_offset])
			, [a] "=a" (a)
			: [distributor] "r" (distributor)
			, [page] "r" (page)
			, [page_count] "rm" (distributor->page_count)
			, [offset] "N" (sizeof(struct distributor_s))
			: "memory", "cc"
		);
	}
}
#else
*/
/*
static inline void distributor_put_sync(struct distributor_s * const __restrict distributor, void* const page ) {
	if (page) {
		size_t write_offset;
		for (;;) {
			write_offset = distributor->write_offset & distributor->modulus;
			if (!distributor->free_pages[write_offset])
				if (__sync_bool_compare_and_swap((uintptr_t*)&(distributor->free_pages[write_offset]), (uintptr_t)0, (uintptr_t)page))
					break;
		}
		
		distributor->write_offset++;
	}
}
*/
// #endif

static inline void* distributor_get(struct distributor_s * const __restrict distributor) {
	size_t read_offset = distributor->read_offset & distributor->modulus;
	if (distributor->write_offset != distributor->read_offset) {
		void* const page = distributor->free_pages[read_offset];
		distributor->free_pages[read_offset] = NULL;
		// read_offset++;
		// if (read_offset > distributor->page_count)
		// 	read_offset = 0;
		distributor->read_offset++;
		return page;
	} else
		return NULL;
}

/*
static inline void* distributor_get_sync(struct distributor_s * const __restrict distributor) {
	size_t read_offset;
	void * volatile page;
	for (;;) {
		
		read_offset = distributor->read_offset;
		// current_size = ;
		if (distributor->write_offset == read_offset)
			return NULL;
		read_offset &= distributor->modulus;
		page = distributor->free_pages[read_offset];
		if (page)
			if (__sync_bool_compare_and_swap((uintptr_t*)&(distributor->free_pages[read_offset]), (uintptr_t)page, (uintptr_t)NULL))
				break;
	}
	// read_offset++;
	// if (read_offset > distributor->page_count)
	// 	read_offset = 0;
	distributor->read_offset++;
	return page;
}
*/

static inline struct distributor_s * distributor_create(size_t num_pages, size_t chunk_size) {
	
	struct distributor_s * distributor;
	size_t modulus;
	size_t temp = num_pages + 1;
	do {
		modulus = temp;
		temp--;
		temp &= modulus;
	} while (temp);
	modulus <<= 1;
	
	size_t distributor_mem = sizeof(*distributor) + ( sizeof(distributor->free_pages[0]) * modulus);
	// modulus--;
	
	if (distributor_mem % sysconf(_SC_PAGE_SIZE)) distributor_mem += sysconf(_SC_PAGE_SIZE) - (distributor_mem % sysconf(_SC_PAGE_SIZE));
	
	distributor = mmap(NULL, distributor_mem, PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE , -1, 0);
	if (distributor == MAP_FAILED)
		return NULL;
	
	distributor->pages = mmap(NULL, sysconf(_SC_PAGE_SIZE) * num_pages * chunk_size, PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE , -1, 0);
	if (distributor->pages == MAP_FAILED)
		return NULL;
	
	distributor->chunk_size = chunk_size;
	distributor->modulus = modulus-1;
	distributor->read_offset = 0;
	distributor->write_offset = 0;
	distributor->page_count = num_pages;
	// #ifdef _POSIX_THREADS
	// pthread_mutex_init(&distributor->mutex, NULL);
	// #endif
	memset((void*)distributor->free_pages, 0, sizeof(distributor->free_pages[0]) * modulus);
	
	size_t i;
	for (i = 0; i < num_pages; i++) {
		distributor_put(distributor, (char*)distributor->pages + (i * sysconf(_SC_PAGE_SIZE) * chunk_size) );
	}
	
	return distributor;
}

static inline int distributor_destroy(struct distributor_s * const __restrict distributor) {
	size_t distributor_mem = sizeof(*distributor) + ( sizeof((*distributor).free_pages[0]) * (distributor->page_count + 1));
	if (distributor_mem % sysconf(_SC_PAGE_SIZE)) distributor_mem += sysconf(_SC_PAGE_SIZE) - (distributor_mem % sysconf(_SC_PAGE_SIZE));
	
	if (munmap(distributor->pages, sysconf(_SC_PAGE_SIZE) * distributor->page_count * distributor->chunk_size))
		return -1;
	
	return munmap(distributor, distributor_mem);
}

#endif
