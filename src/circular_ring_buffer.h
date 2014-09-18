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


#ifndef CIRCULAR_RING_BUFFER_H
#define CIRCULAR_RING_BUFFER_H
#ifdef __MACH__
#include <mach/mach.h>
#else
#include <pthread.h>
#include <fcntl.h>
#endif
#include <sys/mman.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <stdio.h>


#define CIRCULAR_RING_BUFFER_WRITE_WHOLE 0x01
#define CIRCULAR_RING_BUFFER_READ_WHOLE 0x01

struct circular_ring_buffer_s {
	void* address;
	
	size_t length;
	// volatile size_t current_size;
	volatile size_t read_offset;
	volatile size_t write_offset;
	volatile int eof;
};

typedef struct circular_ring_buffer_s circular_ring_buffer;

// int circular_ring_buffer_init(struct circular_ring_buffer_s * circular_ring_buffer, size_t length);
// void circular_ring_buffer_free(struct circular_ring_buffer_s * circular_ring_buffer);

static inline void circular_ring_buffer_free(struct circular_ring_buffer_s * const circular_ring_buffer) {
	#ifdef __MACH__
	vm_deallocate(mach_task_self(), (vm_address_t)circular_ring_buffer->address, circular_ring_buffer->length << 1);
	#else
	munmap(circular_ring_buffer->address, circular_ring_buffer->length << 1);
	#endif
	
	circular_ring_buffer->address = NULL;
}

static inline int circular_ring_buffer_active(const struct circular_ring_buffer_s * const circular_ring_buffer) {
	return (circular_ring_buffer->write_offset != circular_ring_buffer->read_offset) || !(circular_ring_buffer->eof);
}

#define circular_ring_buffer_write_eof(circular_ring_buffer) ((circular_ring_buffer)->eof = 1)

static inline char* circular_ring_buffer_get_write_ptr(const struct circular_ring_buffer_s * const circular_ring_buffer) {
	return (char *)circular_ring_buffer->address + (circular_ring_buffer->write_offset & (circular_ring_buffer->length -1));
}

static inline size_t circular_ring_buffer_get_current_size(const struct circular_ring_buffer_s * const circular_ring_buffer) {
	size_t current_size = circular_ring_buffer->write_offset - circular_ring_buffer->read_offset;
	current_size &= circular_ring_buffer->length -1;
	return current_size;
}

static inline size_t circular_ring_buffer_get_free_space(const struct circular_ring_buffer_s * const circular_ring_buffer) {
	size_t free_space = circular_ring_buffer->read_offset - circular_ring_buffer->write_offset - 1;
	free_space &= circular_ring_buffer->length -1;
	
	return free_space;
}

// static inline void circular_ring_buffer_flush_write(struct circular_ring_buffer_s * const circular_ring_buffer, size_t bytes_written) {
// 	circular_ring_buffer->write_offset += bytes_written;
// }

#define circular_ring_buffer_flush_write(crb, bytes_written) ((crb)->write_offset += (bytes_written))

static inline int circular_ring_buffer_write(struct circular_ring_buffer_s * const circular_ring_buffer, const void* const source, const size_t length, const int flags) {
	size_t free_space = circular_ring_buffer_get_free_space(circular_ring_buffer);
	
	if (free_space >= (flags & CIRCULAR_RING_BUFFER_WRITE_WHOLE ? length : 1)) {
		size_t bytes_written = length;
		if (free_space < bytes_written)
			bytes_written = free_space;
		memcpy(circular_ring_buffer_get_write_ptr(circular_ring_buffer), source, bytes_written);
		circular_ring_buffer_flush_write(circular_ring_buffer, bytes_written);
		return bytes_written;
	} else
		return 0;
}

static inline char* circular_ring_buffer_get_read_ptr(const struct circular_ring_buffer_s * const circular_ring_buffer) {
	return (char *)circular_ring_buffer->address + (circular_ring_buffer->read_offset & (circular_ring_buffer->length -1));
}

// static inline void circular_ring_buffer_flush_read(struct circular_ring_buffer_s * const circular_ring_buffer, size_t bytes_read) {
// 	circular_ring_buffer->read_offset += bytes_read;
// }

#define circular_ring_buffer_flush_read(crb, bytes_read) ((crb)->read_offset += (bytes_read))

static inline int circular_ring_buffer_read(struct circular_ring_buffer_s * const circular_ring_buffer, void* const target, const size_t length, const int flags) {
	size_t current_size = circular_ring_buffer_get_current_size(circular_ring_buffer);
	
	if (current_size >= (flags & CIRCULAR_RING_BUFFER_READ_WHOLE ? length : 1)) {
		size_t bytes_read = length;
		if (current_size < bytes_read)
			bytes_read = current_size;
		memcpy(target, circular_ring_buffer_get_read_ptr(circular_ring_buffer), bytes_read);
		circular_ring_buffer_flush_read(circular_ring_buffer, bytes_read);
		return bytes_read;
	} else
		return 0;
}

static inline int circular_ring_buffer_init(struct circular_ring_buffer_s * circular_ring_buffer, size_t length) {
	// round to page size
	if (length < sysconf(_SC_PAGE_SIZE))
		length = sysconf(_SC_PAGE_SIZE);
	else if (length & (length -1)) {
		while (length & (length -1))
			length = length & (length -1);
		
		length <<= 1;
	}
	
	if (!length || !(length << 1)) {
		fprintf(stderr, "Requested crb length too long\n");
		abort();
	}
	
	#ifdef __MACH__
	// is MACH (OS X)
	
	// fprintf(stderr, "allocating 0x%X memory.\n", length);
	
	int attempts = 3;
	do {
		vm_address_t tmp_address;
		// allocate twice the buffer size
		kern_return_t result = vm_allocate(mach_task_self(), &tmp_address, length << 1, VM_FLAGS_ANYWHERE);

		if ( result != ERR_SUCCESS ) {
	        if ( attempts == 1 ) {
				fprintf(stderr, "Error allocating initial 0x%zu memory.\n", length << 1);
	            return -1;
	        }
	        // Try again
	        continue;
	    }
	
		// deallocate second half of address apce
		result = vm_deallocate(mach_task_self(), tmp_address + length, length);
        if ( result != ERR_SUCCESS ) {
            if ( attempts == 1 ) {
				fprintf(stderr, "Error deallocating second half of address space.\n");
	            return -1;
	        }
			
            vm_deallocate(mach_task_self(), tmp_address, length);
            continue;
        }
		
		vm_address_t mirrorAddress = tmp_address + length;
        vm_prot_t cur_prot, max_prot;
        result = vm_remap(mach_task_self(),
                          &mirrorAddress,   // mirror target
                          length,    // size of mirror
                          0,                 // auto alignment
                          0,                 // force remapping to mirrorAddress
                          mach_task_self(),  // same task
                          tmp_address,     // mirror source
                          0,                 // MAP READ-WRITE, NOT COPY
                          &cur_prot,         // unused protection struct
                          &max_prot,         // unused protection struct
                          VM_INHERIT_DEFAULT);

        if ( result != ERR_SUCCESS ) {
            if ( attempts == 1 ) {
				perror("vm_remap");
				fprintf(stderr, "Error remapping pages.\n");
	            return -1;
	        }
	
            vm_deallocate(mach_task_self(), tmp_address, length);
            continue;
        }
		
		if ( mirrorAddress != tmp_address + length ) {
            if ( attempts == 1 ) {
				fprintf(stderr, "Error remapped pages are at wrong virtual address.\n");
	            return -1;
	        }

            vm_deallocate(mach_task_self(), tmp_address, length);
            vm_deallocate(mach_task_self(), mirrorAddress, length);
            continue;
        }
		
		circular_ring_buffer->address = (void*)tmp_address;
		break;
	} while (attempts--);
	
	
	#else
	// assume Linux
	// Attempt to allocate Anonymous shared non-linear mapping
	{
		
		void *main_address, *tmp_address, *mirrorAddress;
		
		tmp_address = mmap (NULL, length << 1, PROT_READ | PROT_WRITE,
		                        MAP_ANONYMOUS | MAP_SHARED, -1, 0);

		if (tmp_address == MAP_FAILED) {
			return -1;
		}
		
		mirrorAddress = (char *)tmp_address + length;
		
		
		int status = -1;
		
		// depending on the Linux kernel this may fail or succeed. Different versions of Linux support this differently
		#ifdef _GNU_SOURCE
		status = remap_file_pages(mirrorAddress, length, 0/*prot*/, 0/*offset*/, 0/*flags*/);
		#endif
		
		if (status) {
			// remap failed, unmap the region and try shared file mapping instead
			munmap(tmp_address, length << 1);
			
			int fd = -1;
			
			char shmname[256];
			sprintf(shmname, "crb_t_%d_%zx", getpid(), pthread_self());

			fd = shm_open(shmname, O_RDWR | O_CREAT | O_EXCL, 0);
			if (fd < 0) {
				perror("shm_open");
				abort();
			}

			shm_unlink(shmname);

			status = ftruncate (fd, length);
			if (status) {
				return -1;
			}

			tmp_address = mmap (NULL, length << 1, PROT_NONE,
			                        MAP_ANONYMOUS | MAP_PRIVATE, -1, 0);

			if (tmp_address == MAP_FAILED) {
				return -1;
			}

			main_address = mmap (tmp_address, length, PROT_READ | PROT_WRITE,
			        MAP_FIXED | MAP_SHARED, fd, 0);

			if (main_address != tmp_address)	{
				munmap(tmp_address, length << 1);
				return -1;
			}

			mirrorAddress = mmap ((char *)tmp_address + length,
			                length, PROT_READ | PROT_WRITE,
			                MAP_FIXED | MAP_SHARED, fd, 0);

			if (mirrorAddress != (char *)tmp_address + length) {
				munmap(tmp_address, length << 1);
				return -1;
			}
						
			status = close (fd);
			if (status) {
				return -1;
			}
		}

		circular_ring_buffer->address = tmp_address;
	}

	#endif
	
	circular_ring_buffer->read_offset = 0;
	circular_ring_buffer->write_offset = 0;
	circular_ring_buffer->eof = 0;
	circular_ring_buffer->length = length;
	
	return 0;
}


#endif
