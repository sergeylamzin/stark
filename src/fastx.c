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


#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <ctype.h>
#include <string.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <unistd.h>
#include <sys/stat.h>
#include "fastx.h"

#ifndef MAP_ANONYMOUS
#define MAP_ANONYMOUS MAP_ANON
#endif

#ifndef MAP_NORESERVE
#define MAP_NORESERVE 0
#endif

#ifdef __MACH__
char fastx_read_error_message_str[512];
#else
__thread char fastx_read_error_message_str[512];
#endif


#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

static inline size_t roundup(size_t num, size_t modulus) {
	size_t ret = num % modulus;
	if (ret)
		num += modulus - ret;
	return num;
}

ssize_t fasta_read_scan(FILE* fp, struct fastx_read_s * dest, int flags) {
	// read multiline FASTA
	size_t seqlen = 0;
	for (;;) {
		
		int next;
		do {
			next = fgetc(fp);
		} while (isspace(next));
		
		if ((next < 'A' || next > 'z')) {
			if (next != EOF && ungetc(next, fp) == EOF)
				{
					perror("ungetc");
					return -1;
				}
			
			break;
		}
		
		char * line = fgets(dest->data + dest->offset.seq + seqlen, 0x1000, fp);
		if (!line)
			break;
				
		size_t linelen = strlen(line);
		for (; linelen && isspace(line[linelen-1]); linelen--)
			line[linelen-1] = '\0';

		seqlen += linelen;
	}
	
	dest->length.seq = seqlen;
	
	return seqlen;
}

// read 4-line FASTQ, no multiline FASTQ allowed
ssize_t fastq_read_scan(FILE* fp, struct fastx_read_s * dest, int flags) {
	size_t seqlen = 0;
	char * seq = dest->data + dest->offset.seq;

	for (;;) {
		char * line = fgets(seq + seqlen, 0x1000, fp);
		if (!line) {
			perror("fgets");
			return -1;
		}
		
		size_t linelen = strlen(line);
		
		seqlen += linelen;
		
		if (line[linelen-1] == '\n')
			break;
		
		
	}
	
	if (flags & FASTQ_NOSEQ) {
		seqlen = 0;
	} else {
		for (; seqlen && isspace(seq[seqlen-1]); seqlen--)
			seq[seqlen-1] = '\0';
	}
	
	
	dest->length.seq = seqlen;
			
	dest->offset.qual = dest->offset.seq + seqlen;
	char * qual = dest->data + dest->offset.qual;
	
	for (;;) {
		char * line = fgets(qual, 0x1000, fp);
		if (!line) {
			perror("fgets");
			return -1;
		}
		
		size_t linelen = strlen(line);
		
		if (line[linelen-1] == '\n')
			break;
	}
	
	size_t quallen = 0;
	
	for (;;) {
		char * line = fgets(qual + quallen, 0x1000, fp);
		if (!line) {
			perror("fgets");
			return -1;
		}
		
		size_t linelen = strlen(line);
		
		quallen += linelen;
		
		if (line[linelen-1] == '\n')
			break;
	}
	
	if (flags & FASTQ_NOQUAL){
		quallen = 0;
	}
	else
		for (; quallen && isspace(qual[quallen-1]); quallen--)
			qual[quallen-1] = '\0';
	
	dest->length.qual = quallen;
	
	return seqlen + quallen;
}

ssize_t fastx_read_scan(FILE* fp, struct fastx_read_s * dest, int flags) {
	size_t bufsize = sizeof(struct fastx_read_s);
	// int status;
		
	int format;
	
	do {
		format = fgetc(fp);
	} while (isspace(format));
		
	switch (format) {
		case EOF:
		return 0;
		case '>': 
		case '@': break;
		default: sprintf(fastx_read_error_message_str, "scan_fastx_read: unrecognized format\n"); return -2;
	}
	
	char * name;
	
	if (!(name = fgets(dest->data, 0x1000, fp))) {
		if (!feof(fp)) {
			// perror(__FILE__ ":" TOSTRING(__LINE__) ": " "fgets");
			return -1;
		}
		return 0;
	}
	
	// chomp(name);
	
	size_t name_len = 0;
	ssize_t seqlen = 0;
	
	if (!(flags & FASTX_NONAME)) {
		name_len = strlen(name);
		for (; name_len && isspace(name[name_len-1]); name_len--)
			name[name_len-1] = '\0';

		bufsize += name_len + 1;

		dest->offset.seq = name_len + 1;
	} else {
		dest->offset.seq = 0;
	}
	
	dest->offset.name = 0;
	dest->length.name = name_len;
	
	
	switch (format) {
		case '>': 
		seqlen = fasta_read_scan(fp, dest, flags);
		if (seqlen < 0) {
			fprintf(stderr, "error reading FASTA\n");
			return -1;
		}
		return bufsize + seqlen;
		case '@': 
		seqlen = fastq_read_scan(fp, dest, flags);
		if (seqlen < 0) {
			fprintf(stderr, "error reading FASTQ\n");
			return -1;
		}
		return bufsize + seqlen;
		default: sprintf(fastx_read_error_message_str, "scan_fastx_read: unrecognized format\n"); return -2;
	}
}


/*
ssize_t fastx_read_scan(FILE* fp, struct fastx_read_s * dest) {
	
	// char debug_message[4096] = {'\0'};
	fastx_read_error_message_str[0] = '\0';
	
	size_t bufsize = sizeof(struct fastx_read_s);
	int status;
		
	char format;
	
	status = fscanf(fp, "%c", &format);
	if (status == EOF)
		return 0;
	else if (!status) {
		// perror(__FILE__ ":" TOSTRING(__LINE__) ": " "fscanf");
		return -1;
	}
	
	switch (format) {
		case '>': 
		case '@': break;
		default: sprintf(fastx_read_error_message_str, "scan_fastx_read: unrecognized format\n"); return -2;
	}
	
	char * name;
	
	if (!(name = fgets(dest->data, 0x1000, fp))) {
		if (!feof(fp)) {
			// perror(__FILE__ ":" TOSTRING(__LINE__) ": " "fgets");
			return -1;
		}
		return 0;
	}
	
	// chomp(name);
	
	size_t name_len = strlen(name);
	for (; name_len && isspace(name[name_len-1]); name_len--)
		name[name_len-1] = '\0';
	
	dest->offset.name = 0;
	dest->length.name = name_len;
	
	bufsize += name_len + 1;
	// buflen-= name_len + 1;
	
	// fill in name
	
	size_t seqlen = 0;
	dest->offset.seq = name_len + 1;
	char * seq = dest->data + name_len + 1;
	
	int before, after;
	status = fscanf(fp, " %n%s%n ", &before, seq, &after);
	if (status == EOF) {
		sprintf(fastx_read_error_message_str, "unrecognized format, expected sequence, got end of file\n");
		return 0;
	} if (!status) {
		// perror(__FILE__ ":" TOSTRING(__LINE__) ": " "fscanf");
		return -1;
	}
	
	
		
	for (;;) {
		seqlen += after - before;

		bufsize += seqlen;
		// buflen-= seqlen;
		
		status = fscanf(fp, "%c", seq + seqlen);
		if (status == EOF)
			break;
		else if (!status) {
			// perror(__FILE__ ":" TOSTRING(__LINE__) ": " "fscanf");
			break;
		}
		
		ungetc(seq[seqlen], fp);
		// fprintf(stderr, __FILE__ ":" TOSTRING(__LINE__) ": " "next line begins with %c\n", seq[seqlen]);
		
		if (seq[seqlen] == '+') {
			if (format == '>') {
				sprintf(fastx_read_error_message_str, "Encountered weird FASTA read \"%s\" with quality scores, assuming FASTQ\n", dest->data);
			}
			
			dest->length.seq = seqlen;
			char * qual = seq + seqlen;
			dest->offset.qual = dest->offset.seq + seqlen;
			// get the rest of the quality head line
			if (!fgets(qual, 0x1000, fp)) {
				if (!feof(fp)) {
					// perror(__FILE__ ":" TOSTRING(__LINE__) ": " "fgets");
				}
				else
					sprintf(fastx_read_error_message_str, "Encountered end of file where quality scores expected in read \"%s\"\n", dest->data);
				break;
			}
			
			// read uality scores
			status = fscanf(fp, " %n%s%n ", &before, qual, &after);
			if (status == EOF)
				break;
			else if (!status) {
				sprintf(fastx_read_error_message_str, "scan_fastx_read: unrecognized format, expected quality scores, got %s\n", feof(fp) ? "end of file" : "read error");
				
				break;
			}
			
			seqlen = 0;
			
			for (;;) {
				seqlen += after - before;

				bufsize += seqlen;
				// buflen-= seqlen;
				
				status = fscanf(fp, "%c", qual + seqlen);
				if (status == EOF)
					break;
				if (!status) {
					// perror(__FILE__ ":" TOSTRING(__LINE__) ": " "fscanf");
					break;
				}
				
				ungetc(qual[seqlen], fp);
				// sprintf(debug_message, __FILE__ ":" TOSTRING(__LINE__) ": " "next line begins with %c\n", qual[seqlen]);
				
				if (qual[seqlen] == '@' || qual[seqlen] == '>') {
					break;
				}
				
				status = fscanf(fp, " %n%s%n ", &before, qual + seqlen, &after);
				if (status == EOF)
					break;
				if (!status) {
					// perror(__FILE__ ":" TOSTRING(__LINE__) ": " "fscanf");
					break;
				}
			}
			
			dest->length.qual = seqlen;
			if (dest->length.seq != seqlen) {
				sprintf(fastx_read_error_message_str, "Unexpected inequality between sequence and quality string lengths: %zu, %zu\n", dest->length.seq, seqlen);
			}
			return bufsize;
		} else if (seq[seqlen] == '>' || seq[seqlen] == '@') {
			if (format == '@') {
				sprintf(fastx_read_error_message_str, "Encountered next read where quality scores expected in FASTQ read \"%s\"\n", dest->data);
			}
			break;
		} 
		
		status = fscanf(fp, " %n%s%n ", &before, seq + seqlen, &after);
		if (status == EOF) {
			if (format != '>') {
				sprintf(fastx_read_error_message_str, "Encountered end of file where quality scores expected in read \"%s\"\n", dest->data);
			}
		}
		if (!status) {
			// perror(__FILE__ ":" TOSTRING(__LINE__) ": " "fscanf");
			break;
		}
	}
	
	
	dest->length.seq = seqlen;
	dest->length.qual = 0;
	return bufsize;
		
}
*/

void fastx_read_trim(struct fastx_read_s * dest, struct fastx_read_s * fastx_read, int trimN, int trimQ) {
	*dest = *fastx_read;
	
	if (!fastx_read->length.seq || (fastx_read->length.qual && fastx_read->length.seq != fastx_read->length.qual))
		return;
	
	const char * seq = fastx_read_get_seq(fastx_read);
	size_t read_len = fastx_read->length.seq;
	
	
	if (fastx_read->length.qual && trimQ >= 0) {
		char * qual = fastx_read_get_qual(fastx_read);
		for (; read_len && ((trimN & (seq[read_len-1] == 'N')) || (qual[read_len-1] <= trimQ)); read_len--) ;

		for (; read_len && ((trimN & (*seq == 'N')) || (*qual <= trimQ)); read_len--) {
			seq++;
			qual++;
			dest->offset.seq++;
			dest->offset.qual++;
		}
		
		dest->length.seq = read_len;
		dest->length.qual = read_len;
		
	} else if (trimN) {
		for (; read_len && (seq[read_len-1] == 'N'); read_len--) ;

		for (; read_len && (*seq == 'N'); read_len--) {
			seq++;
			dest->offset.seq++;
		}
		dest->length.seq = read_len;
	}
}

struct fastx_read_list_s * fastx_read_list_create(FILE * fp1, FILE * fp2) {
	size_t memory_struct = roundup(sizeof(struct fastx_read_list_s) + (sizeof(void*) * (1 << 28)), sysconf(_SC_PAGE_SIZE));
	struct fastx_read_list_s * fastx_read_list = mmap(NULL
			, memory_struct
			, PROT_NONE
			, MAP_ANONYMOUS | MAP_PRIVATE | MAP_NORESERVE
			, -1
			, 0
		);
	
	if (fastx_read_list == MAP_FAILED) {
		perror(__FILE__ ":" TOSTRING(__LINE__) ": " "mmap");
		return NULL;
	}
	
	if (mprotect(fastx_read_list, MEMOPRY_PROTECT_CHUNK_SIZE, PROT_READ | PROT_WRITE)) {
		perror(__FILE__ ":" TOSTRING(__LINE__) ": " "mprotect");
		munmap(fastx_read_list, memory_struct);
		return NULL;
	}
	
	fastx_read_list->memory_this = memory_struct;
	fastx_read_list->numreads = 0;
	fastx_read_list->memory_open = 0;
	fastx_read_list->memory_used = 0;
	fastx_read_list->memory_total = 1UL << 38;

	struct stat fpstat;
	int fd;
	size_t filesize = 0;
	if (fp1 && (fd = fileno(fp1)) >= 0 && !fstat(fd, &fpstat)) {
		filesize += fpstat.st_size;
	}
	if (fp2 && (fd = fileno(fp2)) >= 0 && !fstat(fd, &fpstat)) {
		filesize += fpstat.st_size;
	}
	
	void * structs;
	
	structs =  mmap(NULL
	#ifdef DEBUG_MEMORY_CAP
			, DEBUG_MEMORY_CAP
	#else
			, 1UL << 38 // 256 GB
	#endif
			, PROT_NONE
			, MAP_ANONYMOUS | MAP_PRIVATE | MAP_NORESERVE
			, -1
			, 0
		);
		
	if (structs == MAP_FAILED) {
		perror(__FILE__ ":" TOSTRING(__LINE__) ": " "mmap");
		return NULL;
	}
	
	if (filesize) {
		filesize = roundup(filesize, MEMOPRY_PROTECT_CHUNK_SIZE);
		if (mprotect( structs
					, filesize
					, PROT_READ | PROT_WRITE)) {
			perror(__FILE__ ":" TOSTRING(__LINE__) ": " "mprotect");
			return NULL;
		}
		fastx_read_list->memory_open = filesize;
	}
	
	fastx_read_list->structs = structs;
	
	
	return fastx_read_list;
}

int fastx_read_list_scan_read(struct fastx_read_list_s * fastx_read_list, FILE * fp, int flags) {
	struct fastx_read_s * const read = (struct fastx_read_s *)(fastx_read_list->structs + fastx_read_list->memory_used);
	
	// size_t new_memory = roundup(sizeof(*read) + len0 + len1, sizeof(void*));
	if (fastx_read_list->memory_open - fastx_read_list->memory_used < MEMOPRY_PROTECT_CHUNK_SIZE) {
		if (mprotect( fastx_read_list->structs + fastx_read_list->memory_open
					, MEMOPRY_PROTECT_CHUNK_SIZE
					, PROT_READ | PROT_WRITE)) {
			perror(__FILE__ ":" TOSTRING(__LINE__) ": " "mprotect");
			return -1;
		}
		fastx_read_list->memory_open += MEMOPRY_PROTECT_CHUNK_SIZE;
	}
	
	ssize_t read_size = fastx_read_scan(fp, read, flags);
	if (read_size <= 0) {
		return read_size;
	}
	
	if (!((sizeof(*fastx_read_list) + (fastx_read_list->numreads * sizeof(void *))) & (MEMOPRY_PROTECT_CHUNK_SIZE - 1)) || 
		((sizeof(*fastx_read_list) + (fastx_read_list->numreads * sizeof(void *))) & (MEMOPRY_PROTECT_CHUNK_SIZE - 1)) + sizeof(void *) > MEMOPRY_PROTECT_CHUNK_SIZE ) {
		if (mprotect( ((char*)fastx_read_list) + roundup(sizeof(*fastx_read_list) + (fastx_read_list->numreads * sizeof(void *)) , MEMOPRY_PROTECT_CHUNK_SIZE)
					, MEMOPRY_PROTECT_CHUNK_SIZE
					, PROT_READ | PROT_WRITE)) {
			perror(__FILE__ ":" TOSTRING(__LINE__) ": " "mprotect");
			return -1;
		}
	}
	fastx_read_list->memory_used += roundup(read_size, 8);
	fastx_read_list->pointers[fastx_read_list->numreads] = read;
	fastx_read_list->numreads++;
	
	return 1;
}

int pair_comparator_qsort(struct fastx_read_s * const pair1[2], struct fastx_read_s * const pair2[2]) {
	int result;
	if (!(result = fastx_read_compare(pair1[0], pair2[0]))) {
		return fastx_read_compare(pair1[1], pair2[1]);
	}
	return result;
}

static void faxtq_read_join_quality(struct fastx_read_s * read1, const struct fastx_read_s * read2, int_fast8_t qual_adjust) {
	size_t length = read1->length.qual;
	if (!length || (length != read2->length.qual))
		return;
	
	char * read1_qual = fastx_read_get_qual(read1);
	const char * read2_qual = fastx_read_get_qual(read2);
	do {
		uint_fast8_t newqual = read1_qual[length - 1] + read2_qual[length - 1] - qual_adjust;
		if (newqual > '~')
			newqual = '~';
		
		read1_qual[length - 1] = newqual;
	} while (--length);
}


int fastx_pair_unique(struct fastx_read_list_s * fastx_read_list, int_fast8_t qual_adjust) {
	if (!fastx_read_list->numreads || (fastx_read_list->numreads & 0x1))
		return -2;
	// sort first
	// size_t numpairs = fastx_read_list->numreads >> 1;
	qsort(fastx_read_list->pointers, fastx_read_list->numreads >> 1, 2 * sizeof(struct fastx_read_s*), (int (*)(const void *, const void *))pair_comparator_qsort);
	
	size_t i, j, memory_for_reads = roundup(fastx_read_sizeof(fastx_read_list->pointers[0]), sizeof(void*)) + roundup(fastx_read_sizeof(fastx_read_list->pointers[1]), sizeof(void*));
	
	size_t backtrace = 2;
	j = 2;
	for (i = 2; i < fastx_read_list->numreads; i+=2) {
		if (!pair_comparator_qsort(fastx_read_list->pointers + i - backtrace, fastx_read_list->pointers + i )) {
			if (qual_adjust > 0)
				faxtq_read_join_quality(fastx_read_list->pointers[i - backtrace], fastx_read_list->pointers[i], qual_adjust);
			
			fastx_read_list->pointers[i] = NULL;
			fastx_read_list->pointers[i+1] = NULL;
			backtrace+=2;
		} else {
			if (backtrace > 2) {
				// fprintf(stderr, "%zu instances of a duplicate read deleted\n", backtrace - 1);
			}
			fastx_read_list->pointers[i - backtrace]->dupecount = backtrace >> 1;
			fastx_read_list->pointers[i - backtrace + 1]->dupecount = backtrace >> 1;
			memory_for_reads += roundup(fastx_read_sizeof(fastx_read_list->pointers[i]), sizeof(void*)) + roundup(fastx_read_sizeof(fastx_read_list->pointers[i+1]), sizeof(void*));
			backtrace = 2;
			fastx_read_list->pointers[j] = fastx_read_list->pointers[i];
			fastx_read_list->pointers[j+1] = fastx_read_list->pointers[i+1];
			j+=2;
		}
	}
	fastx_read_list->pointers[i - backtrace]->dupecount = backtrace >> 1;
	fastx_read_list->pointers[i - backtrace + 1]->dupecount = backtrace >> 1;
	
	// for (j = 0; j < fastx_read_list->numpairs && fastx_read_list->pointers[j]; j++);
	// 
	// for (i = j+1; i < fastx_read_list->numpairs; i++) {
	// 	if (fastx_read_list->pointers[i])
	// 		fastx_read_list->pointers[j++] = fastx_read_list->pointers[i];
	// }
	fastx_read_list->numreads = j;
	
	size_t useful_memory = roundup(sizeof(*fastx_read_list) + (sizeof(void *) * fastx_read_list->numreads), sysconf(_SC_PAGE_SIZE));
	munmap(((char *)fastx_read_list) + useful_memory, fastx_read_list->memory_this - useful_memory);
	fastx_read_list->memory_this = useful_memory;
	
	size_t memory_total = roundup(memory_for_reads, sysconf(_SC_PAGE_SIZE));
	char * structs =  mmap(NULL
			, memory_total
			, PROT_READ | PROT_WRITE
			, MAP_ANONYMOUS | MAP_PRIVATE
			, -1
			, 0
		);
	
	if (structs == MAP_FAILED) {
		perror(__FILE__ ":" TOSTRING(__LINE__) ": " "mmap");
		return -1;
	}
	
	size_t memory_used = 0;
	for (i = 0; i < fastx_read_list->numreads; i++) {
		size_t read_memory = roundup(fastx_read_sizeof(fastx_read_list->pointers[i]), sizeof(void*));
		void * fastx_read = (struct fastx_read_s *)(structs + memory_used);
		memcpy(fastx_read, fastx_read_list->pointers[i], read_memory);
		fastx_read_list->pointers[i] = fastx_read;
		memory_used += read_memory;
	}
	
	munmap(fastx_read_list->structs, fastx_read_list->memory_total);
	
	fastx_read_list->structs = structs;
	fastx_read_list->memory_used = memory_used;
	fastx_read_list->memory_open = memory_total;
	fastx_read_list->memory_total = memory_total;
	
	return 0;
}

/*


#include <printf.h>

int fastx_format (FILE *stream, const struct printf_info *info, const void *const *args) {
	struct fastx_read_s * const read = (struct fastx_read_s *)args[0];
	
	if (info->is_long_double || info->is_long) {
		return fastq_read_fprint(stream, read);
	} else {
		return fasta_read_fprint(stream, read);
	}
}

int fastx_arginfo (const struct printf_info *info, size_t n, int *argtypes) {
	if (n > 0)
	         argtypes[0] = PA_POINTER;
	return 1;
}

*/



