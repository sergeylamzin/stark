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

#include <pthread.h>
#include <stdlib.h>

#define STARK_MAIN_STATUS_C

// struct stark_thread_status_s * stark_thread_status = NULL;
// #pragma omp threadprivate(stark_thread_status)

#include "main.h"
#include "stark.h"
#include "stark_status.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <sys/socket.h>
#include <sys/un.h>
#include <netdb.h>
#include <sched.h>
#include "threadpool.h"
#include "distributor.h"
#include <sys/syscall.h>
#include "stark_status_xsl.h"
#include <sys/un.h>
#include "stark_assemble_hierarchical.h"
#include <dlfcn.h>
#include <signal.h>

#ifdef PBS_STATUS
#include <pbs_error.h>
#include <pbs_ifl.h>

static int pbs_server_sockfd = -1;
#endif

const char stark_status_action_inserting_read[] = "Inserting read.";
volatile struct stark_status_s stark_status;
// #ifndef __MACH__

// #endif

#ifndef MSG_NOSIGNAL
#define MSG_NOSIGNAL 0
#endif

static pthread_t status_thread;
#define STATUS_BUFFER_MAXSIZE 0x1000000
static char* status_buffer = NULL;
static size_t status_buffer_size = 0;

#define stark_status_append(...) (status_buffer_size += snprintf(status_buffer + status_buffer_size,  STATUS_BUFFER_MAXSIZE - status_buffer_size, __VA_ARGS__))

struct socklist_s {
	int sockfd;
	struct socklist_s * next;
};
static struct socklist_s * listeners = NULL;

void stark_status_socket_thread(uintptr_t port);

int stark_status_init(uintptr_t port) {
	memset((void*)&stark_status, 0, sizeof(stark_status));
	if (!status_buffer) {
		void * tmp = mmap(NULL, STATUS_BUFFER_MAXSIZE, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS | MAP_NORESERVE, -1, 0);
		if (tmp == MAP_FAILED) {
			perror("mmap");
			return -1;
		}
		status_buffer = tmp;
	}
	
	pthread_mutex_init((pthread_mutex_t*)&stark_status.mutex, NULL);
	if (port)
		pthread_create(&status_thread, NULL, (void * (*)(void *))stark_status_socket_thread, (void*)port);
	else
		status_thread = (pthread_t)0;
		
	return 0;
}

void stark_status_kill() {
	if (status_thread)
		pthread_cancel(status_thread);
	struct socklist_s * listeners_i;
	for (listeners_i = listeners; listeners_i; listeners_i = listeners_i->next) {
		// if (FD_ISSET(listeners_i->sockfd, &exceptfds)) {
			close(listeners_i->sockfd);
			listeners_i->sockfd = -1;
		// }
	}
	
	#ifdef PBS_STATUS
	if (pbs_server_sockfd >= 0)
		pbs_disconnect(pbs_server_sockfd);
	#endif
	
	munmap(status_buffer, STATUS_BUFFER_MAXSIZE);
}

void stark_status_register_thread() {
	// #ifdef __MACH__
	// cannot use TLS
	// struct stark_thread_status_s * stark_thread_status;
	// #else
	if (stark_thread_status)
		return;
	// #endif
	stark_thread_status = calloc(sizeof(struct stark_thread_status_s), 1);
	stark_thread_status->thread_id = pthread_self();
	
	// DEBUG_MSG("Registering Thread 0x%lX for stark_status", (uint64_t)pthread_self());
	
	pthread_mutex_lock((pthread_mutex_t*)&(stark_status.mutex));
	
	if (stark_status.stark_thread_status) {
		volatile struct stark_thread_status_s * next = stark_status.stark_thread_status;
		while (next->next)
			next = next->next;

		// stark_thread_status->next = stark_status.stark_thread_status;
		next->next = stark_thread_status;
	} else
		stark_status.stark_thread_status = stark_thread_status;
	
	// stark_status_update_last_cpu();
	stark_thread_status->last_cpu = -1;
	
	#ifdef SYS_gettid
	stark_thread_status->linux_tid = syscall(SYS_gettid);
	#else
	stark_thread_status->linux_tid = -1;
	#endif
	
	#if defined(STARK_COVERAGE_CACHE_TOKEN)
	stark_thread_status->mytoken = stark_coverage_mytoken;
	#endif
	
	// volatile struct stark_thread_status_s * first;
	// do {
	// 	first = stark_status.stark_thread_status;
	// 	stark_thread_status->next = first;
	// } while (!__sync_bool_compare_and_swap(&stark_status.stark_thread_status, first, stark_thread_status));
	
	pthread_mutex_unlock((pthread_mutex_t*)&(stark_status.mutex));
}

/*
void stark_status_unregister_thread() {
	
	#ifdef __MACH__
	struct stark_thread_status_s * stark_thread_status = stark_status.stark_thread_status;
	if (stark_thread_status->thread_id == pthread_self()) {
	#else
	if (stark_status.stark_thread_status == stark_thread_status) {
	#endif
		if (__sync_bool_compare_and_swap(&stark_status.stark_thread_status, stark_thread_status, stark_thread_status->next))
			return;
	}
	
	struct stark_thread_status_s * prev;
	for (prev = stark_status.stark_thread_status; prev , prev = prev->next) {
		#ifdef __MACH__
		struct stark_thread_status_s * stark_thread_status = prev->next;
		if (stark_thread_status->thread_id == pthread_self())
		#else
		if (prev->next == stark_thread_status)
		#endif
			if (__sync_bool_compare_and_swap(&prev->next, stark_thread_status, stark_thread_status->next))
				return;
	}
}
*/

/*
struct stark_thread_status_s {
	struct stark_thread_status_s * next;
	pthread_t thread_id;
	
	struct {
		char * seq4;
		size_t length;
		size_t maxk;
	} current_insert;
	
	size_t total_insert_mis;
	size_t total_insertions;
	
};

*/

size_t stark_status_proc_thread(pid_t tid);

size_t stark_write_thread_status(FILE* fp, volatile struct stark_thread_status_s * stark_thread_status) {
	// size_t status_buffer_size = 0;
	stark_status_append("<stark_thread_status pthread_id=\"0x%lX\">\n", (uintptr_t)stark_thread_status->thread_id);
	
	stark_status_append("<current_action>%s</current_action>\n", stark_thread_status->current_action);
	stark_status_append("<reads_inserted>%zu</reads_inserted>\n", stark_thread_status->reads_inserted);
	if (stark_thread_status->total_insertions)
		stark_status_append("<avg_insert_time><ns>%zu</ns></avg_insert_time>\n", (stark_thread_status->total_insert_mis * 1000) / stark_thread_status->total_insertions);
	
	#ifdef __CPU_SETSIZE
	// sched_getcpu()
	stark_status_append("<sched_getcpu>%d</sched_getcpu>\n", stark_thread_status->last_cpu);
	#endif
	
	if (stark_thread_status->dispatch) {
		stark_status_append("<ha_dispatch><size>%ld</size></ha_dispatch>\n", circular_ring_buffer_get_current_size(&stark_thread_status->dispatch->buffer) / sizeof(struct stark_hierarchical_assembler_neighbours_global_s));
	}
	
	#if defined(STARK_COVERAGE_CACHE_TOKEN)
	if (stark_thread_status->mytoken) {
		stark_status_append("<stark_coverage_cache_mytoken>\n");
		stark_status_append("<Token>0x%llx</Token>\n", stark_thread_status->mytoken->mytoken);
		stark_status_append("<Index>%zu</Index>\n", stark_thread_status->mytoken->cache.size);
		stark_status_append("</stark_coverage_cache_mytoken>\n");
	}
	
	#endif
	
	/*
	if (stark_thread_status->current_action == stark_status_action_inserting_read) {
		stark_status_append("<current_insert>\n");
		// char buffer[2048];
		// int i;
		// for (i = 0; i < stark_thread_status->current_insert.length; i++) {
		// 	buffer[buflen++] = nucleotides[stark_thread_status->current_insert.seq4[i]];
		// }
		// buffer[i] = '\0';
		// stark_status_append("<sequence>Currently unsafe operation</sequence>\n", buffer);
		stark_status_append("<length>%zu</length>\n", stark_thread_status->current_insert.length);
		stark_status_append("<maxk>%zu</maxk>\n", stark_thread_status->current_insert.maxk);
		stark_status_append("</current_insert>\n");
	}
	*/
	
	if (stark_thread_status->linux_tid  > 0)
		stark_status_proc_thread(stark_thread_status->linux_tid);
		
	stark_status_append("</stark_thread_status>\n");
	
	return status_buffer_size;
}

// size_t stark_status_threadpool(FILE* fp) {
// 	// size_t status_buffer_size = 0;
// 	if (stark_status.threadpool) {
// 		stark_status_append("<threadpool>\n");
// 		
// 		stark_status_append("<num_threads>%zu</num_threads>\n", stark_status.threadpool->numthreads);
// 		stark_status_append("<numjobs>%zu</numjobs>\n", stark_status.threadpool->size);
// 		stark_status_append("<running_threads>%zu</running_threads>\n", stark_status.threadpool->running_threads);
// 		int i;
// 		for (i = 0; i < stark_status.threadpool->numthreads; i++)
// 			stark_status_append("<thread>0x%lX</thread>\n", (uintptr_t)stark_status.threadpool->threads[i]);
// 		
// 		
// 		stark_status_append("</threadpool>\n");
// 	}
// 	
// 	return status_buffer_size;
// }

/*

struct distributor_s {
	void* pages;
	size_t page_count;
	size_t chunk_size;
	// #ifdef _POSIX_THREADS
	// pthread_mutex_t mutex;
	// #endif
	volatile size_t read_offset;
	volatile size_t write_offset;
	void * volatile free_pages[0];
};

*/

size_t stark_status_distributor(FILE* fp) {
	// size_t status_buffer_size = 0;
	if (stark_status.distributor) {
		stark_status_append("<distributor>\n");
		
		stark_status_append("<page_count>%zu</page_count>\n", stark_status.distributor->page_count);
		stark_status_append("<chunk_size>%zu</chunk_size>\n", stark_status.distributor->chunk_size);
		stark_status_append("<read_offset>%zu</read_offset>\n", stark_status.distributor->read_offset);
		stark_status_append("<write_offset>%zu</write_offset>\n", stark_status.distributor->write_offset);
		stark_status_append("<modulus>%zu</modulus>\n", stark_status.distributor->modulus);
		stark_status_append("<free_pages>\n");
		int i;
		for (i = 0; i <= stark_status.distributor->modulus; i++)
			stark_status_append("<page index=\"%d\" %s %s>0x%lX</page>\n", i
				, (i == (stark_status.distributor->read_offset & stark_status.distributor->modulus) ? "read=\"here\"" : "")
				, (i == (stark_status.distributor->write_offset & stark_status.distributor->modulus) ? "write=\"here\"" : "")
				, (uintptr_t)stark_status.distributor->free_pages[i]);
		
		stark_status_append("</free_pages>\n");
		stark_status_append("</distributor>\n");
	}
	
	return status_buffer_size;
}

size_t stark_status_proc() {	
	FILE * procfp = fopen("/proc/self/status", "r");
	if (!procfp)
		return status_buffer_size;
	
	stark_status_append("<proc_status>\n");
	
	// char key[256];
	// key[255] = '\0';
	char line[4096];
	// value[4095] = '\0';
	// int result;
	
	
	while (fgets(line, 4096, procfp)) {
		char * key, *value;
		
		for (key = line; isspace(*key); key++);
		
		if ((value = strchr( key, ':'))) {
			*value = '\0';
			
			for (value++; isspace(*value); value++);
			
			char * end = value + strlen(value);
			for (end--; isspace(*end); end--) *end = '\0';
			
			stark_status_append("<%s>%s</%s>\n", key, value, key);
		}
	}
	
	fclose(procfp);
	
	stark_status_append("</proc_status>\n");
	
	return status_buffer_size;
}

size_t stark_status_proc_thread(pid_t tid) {
	char filename[256];
	snprintf(filename, sizeof(filename), "/proc/%d/task/%d/status", getpid(), tid);
	FILE * procfp = fopen(filename, "r");
	if (!procfp)
		return status_buffer_size;
	
	stark_status_append("<proc_status thread_id=\"%d\">\n", tid);
	
	snprintf(filename, sizeof(filename), "/proc/%d/task/%d/stat", getpid(), tid);
	FILE* statfp = fopen(filename, "r");
	if (statfp) {
		int cpuid = -1;
		fscanf(statfp, 
			"%*s " "%*s " "%*s " "%*s " "%*s " "%*s " "%*s " "%*s " // 8
			"%*s " "%*s " "%*s " "%*s " "%*s " "%*s " "%*s " "%*s " // 16
			"%*s " "%*s " "%*s " "%*s " "%*s " "%*s " "%*s " "%*s " // 24
			"%*s " "%*s " "%*s " "%*s " "%*s " "%*s " "%*s " "%*s " // 32
			"%*s " "%*s " "%*s " "%*s " "%*s " "%*s " "%d " "%*s " // 40
			, &cpuid
		);
		stark_status_append("<sched_getcpu>%d</sched_getcpu>\n", cpuid);
		fclose(statfp);
	}
	
	// char key[256];
	// key[255] = '\0';
	char line[4096];
	// value[4095] = '\0';
	// int result;
	
	
	while (fgets(line, 4096, procfp)) {
		char * key, *value;
		
		for (key = line; isspace(*key); key++);
		
		if ((value = strchr( key, ':'))) {
			*value = '\0';
			
			for (value++; isspace(*value); value++);
			
			char * end = value + strlen(value) -1;
			for (; isspace(*end); end--) *end = '\0';
			
			stark_status_append("<%s>%s</%s>\n", key, value, key);
		}
	}
	
	fclose(procfp);
	
	stark_status_append("</proc_status>\n");
	
	return status_buffer_size;
}

size_t stark_status_PBS() {
	// PBS_JOBID
	char * PBS_JOBID = getenv("PBS_JOBID");
	if (!PBS_JOBID)
		return status_buffer_size;
	
	
	
	#ifdef PBS_STATUS
	struct batch_status * pbs_status = pbs_statjob(pbs_server_sockfd, PBS_JOBID, NULL, NULL);
	
	if (pbs_status) {
		
		stark_status_append("<PBS_job id=\"%s\" xmlns:jsdl-hpcpa=\"http://schemas.ggf.org/jsdl/2006/07/jsdl-hpcpa\">\n", PBS_JOBID);
		
		struct attrl * job_attribute;
		
		char resource[256];
		resource[0] = '\0';
		for (job_attribute = pbs_status->attribs; job_attribute; job_attribute = job_attribute->next) {
			if (job_attribute->resource && *(job_attribute->resource)) {
				if (strcmp(job_attribute->name, resource)) {
					if (resource[0])
						stark_status_append("</%s>\n", resource);
					memcpy(resource, job_attribute->name, strlen(job_attribute->name) +1);
					stark_status_append("<%s>\n", job_attribute->name);
				}
				stark_status_append("<%s>%s</%s>\n", job_attribute->resource, job_attribute->value, job_attribute->resource);
			} else {
				if (resource[0]) {
					stark_status_append("</%s>\n", resource);
					resource[0] = '\0';
				}
				stark_status_append("<%s>%s</%s>\n", job_attribute->name, job_attribute->value, job_attribute->name);
			}
		}
		
		if (resource[0]) {
			stark_status_append("</%s>\n", resource);
			resource[0] = '\0';
		}
		
		pbs_statfree(pbs_status);
	} else
	#endif
	{
		stark_status_append("<PBS_job id=\"%s\">\n", PBS_JOBID);
		char * PBS_JOBNAME = getenv("PBS_JOBNAME");
		if (PBS_JOBNAME)
			stark_status_append("<name>%s</name>\n", PBS_JOBNAME);
	}
	
	
	
	stark_status_append("</PBS_job>\n");
	
	return status_buffer_size;
}

size_t stark_status_hierarchical_assembler(struct stark_hierarchical_assembler_s * ha) {
	stark_status_append("<hierarchical_assembler mink=\"%d\">\n", ha->minK);
	
	stark_status_append("<groups>\n");
	stark_status_append("<size>%zu</size>\n", ha->groups.size);
	stark_status_append("</groups>\n");
		
	stark_status_append("</hierarchical_assembler>\n", ha->minK);
}

size_t stark_write_status(FILE* fp) {
	// size_t status_buffer_size = 0;
	volatile struct stark_thread_status_s * stark_thread_status = stark_status.stark_thread_status;
	
	stark_status_append("<?xml version=\"1.0\"?>\n");
	stark_status_append("<?xml-stylesheet type=\"text/xsl\" href=\"/stark_status.xsl\"?>\n");
	stark_status_append("<stark_status>\n");
	stark_status_append("<version>" STARK_VERSION_STRING "</version>\n");
	
	size_t reads_inserted = 0;
	
	#if defined(STARK_COVERAGE_CACHE_TOKEN)
	if (stark_coverage_token_mask) {
		stark_status_append("<stark_coverage_cache_token mask=\"0x%llx\">\n", stark_coverage_token_mask);
		size_t i;
		for (i = 0; i < (sizeof(stark_coverage_token) / sizeof(*stark_coverage_token)); i++) {
			stark_status_append("<token id=\"%zu\">0x%llx</token>\n", i, stark_coverage_token[i]);
		}
		stark_status_append("</stark_coverage_cache_token>\n");
	}
	#endif
	// stark_threads
	stark_status_append("<stark_threads>\n");
	for (stark_thread_status = stark_status.stark_thread_status; stark_thread_status; stark_thread_status = stark_thread_status->next) {
		stark_write_thread_status(fp, stark_thread_status);
		reads_inserted += stark_thread_status->reads_inserted;
	}
	if (stark_status.reads_total) {
		stark_status_append("<reads total=\"%zu\" done=\"%zu\" />\n", stark_status.reads_total, reads_inserted);
	}
	stark_status_append("</stark_threads>\n");
	// stark_status_threadpool(fp);
	stark_status_proc();
	stark_status_PBS();
	stark_status_distributor(fp);
	stark_status_append("</stark_status>\n");
	
	return status_buffer_size;
}

struct sockaddr_un sockaddr_un;

static void unlink_unix_socket(void) {
	unlink(sockaddr_un.sun_path);
}

void stark_status_socket_thread(uintptr_t port) {
	{
		struct sigaction sa;
		memset(&sa, 0, sizeof(sa));
		sa.sa_handler = SIG_IGN;
		sigaction(SIGPIPE, &sa, NULL);	
	}
	
	#ifdef PBS_STATUS
	if (getenv("PBS_JOBID"))
		pbs_server_sockfd = pbs_connect(NULL);
	#endif
	
	struct addrinfo hints, *res, *p;
	
	memset(&hints, 0, sizeof hints);
	hints.ai_family = AF_UNSPEC;  // use IPv4 or IPv6, whichever
	hints.ai_socktype = SOCK_STREAM;
	hints.ai_flags = AI_PASSIVE;     // fill in my IP for me
	
	char portstr[32];
	sprintf(portstr, "%d", (uint16_t)port);
	getaddrinfo(NULL, portstr, &hints, &res);

	// make a socket:
	for(p = res;p != NULL; p = p->ai_next) {
        int sockfd = socket(res->ai_family, res->ai_socktype, res->ai_protocol);
		
		if (sockfd < 0)
			continue;
		
		int yes=1;

		// lose the pesky "Address already in use" error message
		if (setsockopt(sockfd,SOL_SOCKET,SO_REUSEADDR,&yes,sizeof(int)) == -1) {
		    perror("setsockopt");
			continue;
		}
		
		if(bind(sockfd, res->ai_addr, res->ai_addrlen))
			continue;
		
		if (listen(sockfd, 10))
			continue;
		
		struct socklist_s * socklist = calloc(sizeof(struct socklist_s), 1);
		socklist->sockfd = sockfd;
		socklist->next = listeners;
		listeners = socklist;
    }

	freeaddrinfo(res); // free the linked list
	
	#ifdef AF_UNIX
	char * tempdir = getenv("TMPDIR");
	if (!tempdir)
		tempdir = getenv("TMP");
	if (!tempdir)
		tempdir = getenv("TEMP");
	if (!tempdir)
		tempdir = getenv("TEMPDIR");
	
	if (tempdir) {
		sockaddr_un.sun_family = AF_UNIX;
		
		int sockfd = socket(AF_UNIX, SOCK_STREAM, 0);
		
		if (sockfd >= 0) {
			socklen_t len = sizeof(sockaddr_un.sun_family) +
			snprintf(sockaddr_un.sun_path, sizeof(sockaddr_un.sun_path), "%s/stark_status_%d.sock", tempdir, getpid());
			unlink(sockaddr_un.sun_path);
			
			if (!bind(sockfd, (struct sockaddr *)&sockaddr_un, len) && !listen(sockfd, 10)) {
				struct socklist_s * socklist = calloc(sizeof(struct socklist_s), 1);
				socklist->sockfd = sockfd;
				socklist->next = listeners;
				listeners = socklist;

				atexit(unlink_unix_socket);
				
				DEBUG_MSG("Created UNIX status socket on %s\n", sockaddr_un.sun_path);
			} else {
				perror("bind/listen");
			}
		} else {
			perror("socket");
		}
	}
	
	#endif
	

	if (!listeners) {
		fprintf(stderr, "Could not establish listener sockets on port %s\n", portstr);
		return;
	}

	fd_set readfds, exceptfds;
	
	struct timeval timeout;
	
	for (;;) {
		
		FD_ZERO(&readfds);
		FD_ZERO(&exceptfds);
		int nfds = 0;
		
		struct socklist_s * listeners_i;
		for (listeners_i = listeners; listeners_i; listeners_i = listeners_i->next) {
			if (listeners_i->sockfd < 0)
				continue;
			FD_SET(listeners_i->sockfd, &readfds);
			FD_SET(listeners_i->sockfd, &exceptfds);
			if (nfds < listeners_i->sockfd)
				nfds = listeners_i->sockfd + 1;
		}
		
		// DEBUG_MSG("selecting on listener sockets");
		if (status_buffer_size) {
			// if the buffer is non-empty select with a timeout
			timeout.tv_sec = 15;
			timeout.tv_usec = 0;
			nfds = select(nfds, &readfds, NULL, &exceptfds, &timeout);
		} else
			nfds = select(nfds, &readfds, NULL, &exceptfds, NULL);
		
		
		if(nfds < 0) {
			perror("select");
			return;
		} else if (nfds == 0) {
			// return due to timeout
			// free memory associated with the buffer
			madvise(status_buffer, STATUS_BUFFER_MAXSIZE, MADV_DONTNEED);
			status_buffer_size = 0;
			continue;
		}
		
		// DEBUG_MSG("select returned");
		
		for (listeners_i = listeners; listeners_i; listeners_i = listeners_i->next) {
			if (FD_ISSET(listeners_i->sockfd, &exceptfds)) {
				close(listeners_i->sockfd);
				listeners_i->sockfd = -1;
			}
		}
		
		for (listeners_i = listeners; listeners_i; listeners_i = listeners_i->next) {
			if (FD_ISSET(listeners_i->sockfd, &readfds)) {
				struct sockaddr addr;
				socklen_t addrlen = sizeof(struct sockaddr);
				int connectsockfd = accept(listeners_i->sockfd, &addr, &addrlen);
				if (connectsockfd < 0) {
					perror("accept");
					continue;
				}
				// DEBUG_MSG("connection accepted sockfd = %d", connectsockfd);
				// FILE* fp = fdopen(connectsockfd, "w");
				status_buffer_size = snprintf(status_buffer,  STATUS_BUFFER_MAXSIZE, "HTTP/1.1 200 OK\r\n"
							"Connection: close\r\n"
							"Content-Type: text/xml\r\n"
							"Pragma: no-cache\r\n"
							"\r\n");
				// fflush(fp);
				status_buffer_size = 0;
				stark_write_status(NULL);
	
				// truncate tailing '\0'
				status_buffer_size--;
				// puts(status_buffer);
				
				char headers[4096];
				size_t headers_size;
				size_t bytes_sent;
				// DEBUG_MSG("got a new connection");
				
				struct timeval timeout_read;
				timeout_read.tv_sec = 2;
				timeout_read.tv_usec = 0;
				fd_set thisfds;
				FD_ZERO(&thisfds);
				FD_SET(connectsockfd, &thisfds);
				if (select(connectsockfd + 1, &thisfds, NULL, NULL, &timeout_read)) {
					ssize_t bytes_got = recv(connectsockfd, headers, sizeof(headers) - 1, MSG_DONTWAIT);
					// DEBUG_MSG("Client sent %zu bytes", bytes_got);
					if (bytes_got >= 0) {

						headers[bytes_got] = '\0';

						// DEBUG_MSG("Client sent:\n%s", headers);

						if(bytes_got && strstr(headers, "stark_status.xsl")) {
							// send stylesheet
							while (recv(connectsockfd, headers, sizeof(headers), MSG_DONTWAIT) > 0) ;

							headers_size = snprintf(headers,  sizeof(headers), "HTTP/1.1 200 OK\r\n"
										"Connection: close\r\n"
										"Content-Type: text/xsl\r\n"
										"Content-Length: %u\r\n"
										"Server: StarK/" STARK_VERSION_STRING "\r\n"
										"\r\n", stark_status_xsl_len);

							bytes_sent = 0;
							while (bytes_sent < headers_size) {
								ssize_t sent;
								if ((sent = send(connectsockfd, headers + bytes_sent, headers_size - bytes_sent, MSG_NOSIGNAL)) < 0) {
									perror("send");
									break;
								} else
									bytes_sent += sent;
							}

							while (recv(connectsockfd, headers, sizeof(headers), MSG_DONTWAIT) > 0) ;

							bytes_sent = 0;
							while (bytes_sent < stark_status_xsl_len) {
								ssize_t sent;
								if ((sent = send(connectsockfd, stark_status_xsl + bytes_sent, stark_status_xsl_len - bytes_sent, MSG_NOSIGNAL)) < 0) {
									perror("send");
									break;
								} else
									bytes_sent += sent;
							}

						} else {
							// send xml status
							
							if (bytes_got) {
								while (recv(connectsockfd, headers, sizeof(headers), MSG_DONTWAIT) > 0) ;
								headers_size = snprintf(headers,  sizeof(headers), "HTTP/1.1 200 OK\r\n"
											"Connection: close\r\n"
											"Content-Type: text/xml\r\n"
											"Content-Length: %zu\r\n"
											"Pragma: no-cache\r\n"
											"Server: StarK/" STARK_VERSION_STRING "\r\n"
											"\r\n", status_buffer_size);

								bytes_sent = 0;
								while (bytes_sent < headers_size) {
									ssize_t sent;
									if ((sent = send(connectsockfd, headers + bytes_sent, headers_size - bytes_sent, MSG_NOSIGNAL)) < 0) {
										perror("send");
										break;
									} else
										bytes_sent += sent;
								}
							}

							while (recv(connectsockfd, headers, sizeof(headers), MSG_DONTWAIT) > 0) ;

							bytes_sent = 0;
							while (bytes_sent < status_buffer_size) {
								ssize_t sent;
								if ((sent = send(connectsockfd, status_buffer + bytes_sent, status_buffer_size - bytes_sent, MSG_NOSIGNAL)) < 0) {
									perror("send");
									break;
								} else
									bytes_sent += sent;
							}
						}
					}
				}
				
				// DEBUG_MSG("status written, flushing");
				// fflush(fp);
				shutdown(connectsockfd, SHUT_WR);
				// char buffer[4096];
				while (recv(connectsockfd, headers, sizeof(headers), MSG_DONTWAIT) > 0) ;
				// DEBUG_MSG("status flushed");
				close(connectsockfd);
			}
		}
		
	}
	
}




