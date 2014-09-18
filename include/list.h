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
#include <stdio.h>
#include <string.h>

#ifndef LIST_H
#define LIST_H

#include <sys/types.h>
#include <sys/mman.h>
#include <unistd.h>

#ifdef __MACH__
#include <mach/mach.h>
#else
#include <pthread.h>
#include <fcntl.h>
#endif

// #define list_access(listarg, index) ((((listarg).list))[index])
// #define list_pointer(listarg, index) ((((listarg).list)) + index)

#define list_qsort(listarg, comparator) qsort((listarg).list, (listarg).size, sizeof((listarg).list[0]), comparator)

// typedef struct list_s {
// 	void* list;
// 	size_t size;
// 	size_t maxlength;
// 	size_t objectsize;
// 	void* (*allocator) ( size_t ) ;
// } list_t;

// type##_list_t 
#define list_t(type) struct { \
	type *list; \
	ssize_t size; \
	size_t maxlength; \
	int shmfd; \
	void * oldmem; \
} 

// typedef list_t(void) void_list_t;

// __typeof__
#define list_init(listarg) list_init_size( (listarg), 1024)
#define list_init_flags(listarg, flags) list_init_size_flags( (listarg), 1024, (flags))
#define list_init_size( _list, initsize) list_init_size_flags((_list), (initsize), 0)
#define list_init_size_flags( _list, initsize, flags) list_init_size_alloc((_list), (initsize), sizeof((_list)->list[0]), (flags))

// inline void list_free( list_t(void) list);
// #define list_free(listarg) list_free_real(&(listarg), sizeof((listarg).list[0]))

#define list_free(listarg) ((listarg).maxlength * sizeof((listarg).list[0]) >= sysconf(_SC_PAGESIZE)) ? munmap((listarg).list, (listarg).maxlength * sizeof((listarg).list[0])) : free((listarg).list), (listarg).oldmem ? munmap((listarg).oldmem, ((listarg).maxlength * sizeof((listarg).list[0])) / 2) : 0



// static inline void list_free_real(void * const _list, size_t objectsize) {
// 	if (list->maxlength * objectsize >= sysconf(_SC_PAGESIZE)) {
// 		// use mmap
// 		munmap(list->list, list->maxlength * objectsize);
// 	} else {
// 		// use malloc
// 		free(list->list);
// 	}
// }
// 
// #ifdef MREMAP_MAYMOVE
// static void * list_allocator_mmap(size_t bytes) {
// 	void * tmp = mmap(NULL, bytes, PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE, -1, 0);
// 	if (tmp == MAP_FAILED)
// 		return NULL;
// 	return tmp;
// }
// 
// static void * list_allocator_mremap(void * old_address, size_t old_size, size_t new_size) {
// 	void * tmp = mremap(old_address, old_size, new_size, MREMAP_MAYMOVE);
// 	if (tmp == MAP_FAILED)
// 		return NULL;
// 	return tmp;
// }
// #endif

static inline void list_init_size_alloc( void * const _list, size_t initsize, size_t objectsize, int flags) {
	list_t(void) * list = _list;
	list->size = 0;
	list->shmfd = -1;
	list->oldmem = NULL;
	// list->objectsize = objectsize;
	
	if (flags & 1) { // use thread safe shared memory expansion
		size_t memsize = initsize * objectsize;
		
		if (memsize & (sysconf(_SC_PAGESIZE) - 1)) {
			fprintf(stderr, "list_init_size_alloc: initsize * objectsize = %zu * %zu must be page aligned to 0x%lx\n", initsize, objectsize, sysconf(_SC_PAGESIZE));
			abort();
		}
			
		
		#ifdef __MACH__
		{
			kern_return_t result = vm_allocate(mach_task_self(), (vm_address_t *)&list->list, memsize, VM_FLAGS_ANYWHERE);
			if ( result != ERR_SUCCESS ) {
				fprintf(stderr, "Error allocating VM space\n");
				abort();
			}
		}
		list->shmfd = 0;
		#else
		char shmname[256];
		sprintf(shmname, "list_t_%d_%zx", getpid(), pthread_self());

		list->shmfd = shm_open(shmname, O_RDWR | O_CREAT | O_EXCL, 0);
		if (list->shmfd < 0) {
			perror("shm_open");
			abort();
		}
		if (list->shmfd == 0) {
			list->shmfd = dup(list->shmfd);
			if (list->shmfd < 0) {
				perror("dup");
				abort();
			}
		}

		shm_unlink(shmname);

		if (ftruncate(list->shmfd, memsize)) {
			perror("ftruncate");
			abort();
		}

		list->list = mmap(NULL
			, memsize
			, PROT_READ | PROT_WRITE
			, MAP_SHARED | MAP_FILE
			, list->shmfd
			, 0
		);

		if (list->list == MAP_FAILED) {
			perror("mmap");
			abort();
		}
		#endif
	}
	else if (initsize * objectsize >= sysconf(_SC_PAGESIZE)) {
		// round up init size
		// size_t memsize = initsize * objectsize;
		// if (memsize & (sysconf(_SC_PAGESIZE) - 1))
		
		// use mmap
		void * tmp = mmap(
				  NULL
				, initsize * objectsize
				, PROT_READ | PROT_WRITE
				, MAP_PRIVATE | MAP_ANON
				, -1
				, 0
			);
		if (tmp == MAP_FAILED) {
			perror("mmap");
			abort();
		}
		list->list = tmp;
	} else {
		// use malloc
		if(!(list->list = calloc(initsize, objectsize))) {
			perror("malloc");
			abort();
		}
		
	}
	
	list->maxlength = initsize;
}

static inline void * remap_shm_shared_copy(void * old_address, size_t old_size, size_t new_size, int fd) {
	#ifdef __MACH__
	// OS X gets special treatment
	
	mach_port_t mytask = mach_task_self();
	kern_return_t result;
	
	vm_address_t new_address = (uintptr_t)old_address + old_size;
	
	// first try to just extend the VM mapping
	result = vm_allocate(mytask, &new_address, new_size - old_size, VM_FLAGS_FIXED);
	
	if (result == ERR_SUCCESS )
		return old_address;
	
	int attempts;
	result = ERR_SUCCESS;
	
	// ;
	for (attempts = 0; attempts < 3; attempts++, result = vm_deallocate(mytask, new_address + old_size, new_size - old_size)) {
		result = vm_allocate(mytask, &new_address, new_size, VM_FLAGS_ANYWHERE);

		if ( result != ERR_SUCCESS ) {
			fprintf(stderr, "Error allocating new VM space\n");
			abort();
		}

		result = vm_deallocate(mytask, new_address, old_size);

		if ( result != ERR_SUCCESS ) {
			fprintf(stderr, "%d: Error deallocating part of new VM space\n", __LINE__);
			abort();
		}
		
		// there is a race condition here, so we need to have more then one attempt

		vm_prot_t cur_prot, max_prot;

		result = vm_remap(mytask,
			&new_address,   // mirror target
			old_size,    // size of mirror
			0,                 // auto alignment
			0,                 // force remapping to mirrorAddress
			mytask,  // same task
			(uintptr_t)old_address,     // mirror source
			0,                 // MAP READ-WRITE, NOT COPY
			&cur_prot,         // unused protection struct
			&max_prot,         // unused protection struct
			VM_INHERIT_DEFAULT);

		if ( result == ERR_SUCCESS )
			return (void*) new_address;
	}
	
	fprintf(stderr, "attempts = %d\n", attempts);
	
	
	fprintf(stderr, "Error remapping old vm space\n");
	abort();
	return NULL;
	
	#else
	if (ftruncate(fd, new_size)) {
		perror("ftruncate");
		abort();
	}
	
	#ifdef MREMAP_MAYMOVE
	// try to remap first
	if(mremap(old_address, old_size, new_size, 0) != MAP_FAILED) {
		return old_address;
	}
	#endif
	
	void * tmp = mmap(NULL
		, new_size
		, PROT_READ | PROT_WRITE
		, MAP_SHARED
		, fd
		, 0
	);
	
	if (tmp == MAP_FAILED) {
		perror("mmap");
		abort();
	}
	
	return tmp;
	#endif
}


static inline int list_check_length_real(void * const _list, size_t objectsize) {
	list_t(void) * list = _list;
	if (list->size >= list->maxlength) { // enlarge the array
		if (list->shmfd > 0) {
			// we assume that all old pointers have expired by now
			if (list->oldmem) {
				munmap(list->oldmem, (list->maxlength * objectsize) / 2);
				list->oldmem = NULL;
			}
			
			void * oldmem = list->list;
			list->list = remap_shm_shared_copy(list->list, list->maxlength * objectsize, list->maxlength * objectsize * 2, list->shmfd);
			list->maxlength *=2;
			
			if (oldmem != list->list) {
				list->oldmem = oldmem;
			}
			
		} else if (list->maxlength * 2 * objectsize >= sysconf(_SC_PAGESIZE)) {
			if (list->maxlength * objectsize < sysconf(_SC_PAGESIZE)) {
				// switch from malloc to mmap
				void * tmp = mmap(
						  NULL
						, list->maxlength * objectsize * 2
						, PROT_READ | PROT_WRITE
						, MAP_PRIVATE | MAP_ANON
						, -1
						, 0
					);
				if (tmp == MAP_FAILED) {
					perror("mmap");
					return -1;
				}
				memcpy(tmp, list->list, list->maxlength * objectsize);
				free(list->list);
				list->list = tmp;
				list->maxlength *=2;
			} else {
				// stay with mmap
				void * tmp;
				#ifdef MREMAP_MAYMOVE
				tmp = mremap(list->list, list->maxlength * objectsize, list->maxlength * objectsize * 2, MREMAP_MAYMOVE);
				if (tmp == MAP_FAILED)
					return -1;
				#else
				tmp = mmap(
						  NULL
						, list->maxlength * objectsize * 2
						, PROT_READ | PROT_WRITE
						, MAP_PRIVATE | MAP_ANON
						, -1
						, 0
					);
				if (tmp == MAP_FAILED) {
					perror("mmap");
					return -1;
				}
				memcpy(tmp, list->list, list->maxlength * objectsize);
				munmap(list->list, list->maxlength * objectsize);
				#endif
				list->list = tmp;
				list->maxlength *=2;
			}
		} else {
			// stay with malloc
			void *tmp = realloc(list->list, list->maxlength * 2 * objectsize);
			if (tmp) {
			   	list->list = tmp; 
				list->maxlength *=2;
			} else {
				return -1;
			}
		}
	}
	return 0;
}

#define list_check_length(_list) list_check_length_real(_list, sizeof((_list)->list[0]))

#define list_new_empty( _list) (list_check_length(_list) ? -1 : (memset((_list)->list + (_list)->size, 0, sizeof((_list)->list[0])) , (_list)->size++))

#define list_insert_fast( _list, object) (list_check_length(_list) ? -1 : (((_list)->list[(_list)->size] = (object)), (_list)->size++))
#define list_insert list_insert_fast
#define list_push list_insert_fast

#define list_pop( _list) (_list)->list[--((_list)->size)]

static inline void* list_compact_real(void * const _list, size_t objectsize) {
	list_t(void) * list = _list;
	if (list->maxlength * objectsize >= sysconf(_SC_PAGESIZE))
		return list->list;
	
	if (list->size) {
		void *tmp = realloc(list->list, list->size * objectsize);
		if (tmp) {
		   	list->list = tmp; 
			list->maxlength = list->size;
		}
	} else if (list->list) {
		free(list->list);
		list->list = NULL;
		list->maxlength = 0;
	}
	
	return list->list;
}

#define list_compact(_list) list_compact_real(_list, sizeof((_list)->list[0]))

// heap functions

static inline void swap_memory(void * data1, void * data2, size_t size) {
	char temp[size];
	memcpy(temp, data1, size);
	memcpy(data1, data2, size);
	memcpy(data2, temp, size);
}

static inline void list_shiftUp(size_t size, char (*data)[size], size_t start, size_t end, int (*compare)(const void *, const void *)) {
	size_t parent;
	
	while (start && compare(data[parent = (start-1)/2], data[start]) > 0) {
		swap_memory(data[parent], data[start], size);
		start = parent;
	}
}

static inline void list_shiftDown(size_t size, char (*data)[size], size_t start, size_t end, int (*compare)(const void *, const void *)) {
	size_t root = start;
	size_t child;
	
	while ((child = root * 2 +1) <= end) {
		size_t swap = root;
		
		if (compare(data[swap], data[child]) < 0)
			swap = child;
		
		if (child < end && compare(data[swap], data[child+1]) < 0)
			swap = child + 1;
		
		if (swap != root) {
			swap_memory(data[root], data[swap], size);
			root = swap;
		} else
			return;
	}
}

static inline void list_heapify_real(void * const _list, size_t objectsize, int (*compare)(const void *, const void *)) {
	list_t(void) * list = _list;
	ssize_t start;
	for (start = (list->size - 2) / 2; start >= 0; start--)
		list_shiftDown(objectsize, list->list, start, list->size - 1, compare);
}

#endif
