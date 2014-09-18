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


#ifndef FASTX_H_GHLHCTHM
#define FASTX_H_GHLHCTHM

#ifndef MEMOPRY_PROTECT_CHUNK_SIZE
#define MEMOPRY_PROTECT_CHUNK_SIZE 0x100000
#endif

#define FASTX_NONAME 1
#define FASTQ_NOQUAL 2
#define FASTQ_NOSEQ  4

struct fastx_read_s {
	struct {
		size_t name, seq, qual;
	} offset;
	struct {
		size_t name, seq, qual;
	} length;
	
	unsigned short dupecount;
	
	char data[0];
};

struct fastx_read_list_s {
	char * structs;
	size_t memory_used;
	size_t memory_open;
	size_t memory_total;
	size_t numreads;
	size_t memory_this;
	struct fastx_read_s* pointers[0];
};

#define fastx_read_sizeof(read) ({ \
	typeof(read) _read = (read); \
	sizeof(struct fastx_read_s) + _read->offset.qual + _read->length.qual; \
})
// returns size in bytes of the estruct + it's data

#define fastx_read_get_name(read) ({ \
	typeof(read) _read = (read); \
	_read->data + _read->offset.name; \
})
// null-terminated string

#define fastx_read_get_seq(read) ({ \
	typeof(read) _read = (read); \
	_read->data + _read->offset.seq; \
})
// non-null-terminated string
// print with printf("%.*s", read->length.seq, fastx_read_get_seq(read))

#define fastx_read_get_qual(read) ({ \
	typeof(read) _read = (read); \
	_read->data + _read->offset.qual; \
})
// non-null-terminated string
// print with printf("%.*s", read->length.qual, fastx_read_get_qual(read))

ssize_t fastx_read_scan(FILE* fp, struct fastx_read_s * dest, int flags);
/*
	Parses one FASTA or FASTQ read from fp into dest. The user should ensure that *dest is a large enough buffer to hold sizeof(struct fastx_read_s), the read name, sequence and quality strings (if appropriate).
	On success the amount of bytes written to *dest including the header is returned. On failure -1 is returned and *dest is undefined. On EOF 0 is returned and *dest is undefined.
	In certain cases a positive value may be returned even if an error was encountered indicating that up until the error a valid read was read in.
	Consult fastx_read_error_message() for errors.
*/

/*
#define fastx_read_fprint_format(fp, format, read) \
	fprintf(fp	, format \
						, fastx_read_get_name(read) \
						, (int)read->length.seq \
						, fastx_read_get_seq(read) \
						, fastx_read_get_name(read) \
						, (int)read->length.qual \
						, fastx_read_get_qual(read) \
					)

;
*/



static inline size_t fasta_read_fprint(FILE * fp, const struct fastx_read_s * read) { 
	// return fastx_read_fprint_format(fp, ">%sn%.*sn", read); 
	return fprintf(fp	, ">%s\n%.*s\n"
						, fastx_read_get_name(read)
						, (int)read->length.seq
						, fastx_read_get_seq(read)
					);
}
static inline size_t fastq_read_fprint(FILE * fp, const struct fastx_read_s * read) { 
	return fprintf(fp	, "@%s\n%.*s\n+%s\n%.*s\n"
						, fastx_read_get_name(read)
						, (int)read->length.seq
						, fastx_read_get_seq(read)
						, fastx_read_get_name(read)
						, (int)read->length.qual
						, fastx_read_get_qual(read)
					);
}

static inline size_t fastx_read_fprint(FILE * fp, const struct fastx_read_s * read) {
	return read->length.qual ? fastq_read_fprint(fp, read) : fasta_read_fprint(fp, read); 
}
/*
	Prints the read to fp.
	FASTQ is printed if a quality string is available.
	FASTA otherwise.
*/

/*
#define fastx_read_fprint(fp, read) ({ \
	const struct fastx_read_s *_read = (read); \
	fastx_read_fprint_format(fp, _read->length.qual ? "@%s\n%.*s\n+%s\n%.*s\n" : ">%s\n%.*s\n", _read); \
})

#define fasta_read_fprint(fp, read) ({ \
	const struct fastx_read_s *_read = (read); \
	fastx_read_fprint_format(fp, ">%s\n%.*s\n", _read); \
})
#define fastq_read_fprint(fp, read) ({ \
	const struct fastx_read_s *_read = (read); \
	fastx_read_fprint_format(fp, "@%s\n%.*s\n+%s\n%.*s\n", _read); \
})
*/

struct fastx_read_list_s * fastx_read_list_create(FILE * fp1, FILE * fp2);
int fastx_read_list_scan_read(struct fastx_read_list_s * fastx_read_list, FILE * fp, int flags);


static inline int fastx_read_compare(const struct fastx_read_s * read1, const struct fastx_read_s * read2) {
	size_t len = read1->length.seq;
	if (len > read2->length.seq)
		len = read2->length.seq;
	
	int result;
	if ((result = strncmp(
		  fastx_read_get_seq(read1)
		, fastx_read_get_seq(read2)
		, len
		)) == 0)
	{
		return read1->length.seq - read2->length.seq;
	} else
		return result;
}

static inline size_t fastx_read_list_scan_se(struct fastx_read_list_s * fastx_read_list, FILE * fp, int flags) {
	size_t numreads = 0;
	while (fastx_read_list_scan_read(fastx_read_list, fp, flags) > 0) numreads++;
	
	return numreads;
}

static inline ssize_t fastx_read_list_scan_pe(struct fastx_read_list_s * fastx_read_list, FILE * fp[2], int flags) {
	size_t numpairs = 0;
	while (fastx_read_list_scan_read(fastx_read_list, fp[0], flags) > 0) {
		if (fastx_read_list_scan_read(fastx_read_list, fp[1], flags) <= 0)
			return -1;
		
		struct fastx_read_s * read1 = fastx_read_list->pointers[fastx_read_list->numreads -2];
		struct fastx_read_s * read2 = fastx_read_list->pointers[fastx_read_list->numreads -1];
		
		if (fastx_read_compare(read1, read2) > 0)
		{
			fastx_read_list->pointers[fastx_read_list->numreads -2] = read2;
			fastx_read_list->pointers[fastx_read_list->numreads -1] = read1;
		}
		
		numpairs++;
	}
	
	return numpairs;
}

int fastx_pair_unique(struct fastx_read_list_s * fastx_read_list, int_fast8_t qual_adjust);

static inline size_t fastx_read_list_print_se(struct fastx_read_list_s * fastx_read_list, FILE * fp) {
	size_t numreads;
	size_t bytes_printed = 0;
	for (numreads = 0; numreads < fastx_read_list->numreads; numreads++) {
		bytes_printed += fastx_read_fprint(fp, fastx_read_list->pointers[numreads]);
	}
	
	return bytes_printed;
}

static inline size_t fastx_read_list_print_pe(struct fastx_read_list_s * fastx_read_list, FILE * fp[2]) {
	if (fastx_read_list->numreads & 0x1)
		return -1;
	size_t numreads;
	size_t bytes_printed = 0;
	for (numreads = 0; numreads < fastx_read_list->numreads; numreads+=2) {
		bytes_printed += fastx_read_fprint(fp[0], fastx_read_list->pointers[numreads]);
		bytes_printed += fastx_read_fprint(fp[1], fastx_read_list->pointers[numreads+1]);
	}
	
	return bytes_printed;
}

void fastx_read_trim(struct fastx_read_s * dest, struct fastx_read_s * read, int trimN, int trimQ);

/*
static inline size_t fastx_read_list_trim_and_print_se(struct fastx_read_list_s * fastx_read_list, FILE * fp, int_fast8_t trimN, int_fast8_t trimQ, ssize_t minlength, FILE* discardedfp) {
	size_t numreads;
	size_t bytes_printed = 0;
	for (numreads = 0; numreads < fastx_read_list->numreads; numreads++) {
		struct fastx_read_s trimmed_meta_read;
		fastx_read_trim(&trimmed_meta_read, fastx_read_list->pointers[numreads], trimN, trimQ);
		size_t minthislen = minlength >= 0 ? minlength : (fastx_read_list->pointers[numreads]->length.seq >> 1);
		
		if (trimmed_meta_read.length.seq < minthislen) {
			if (discardedfp)
				bytes_printed += fastx_read_fprint(discardedfp, fastx_read_list->pointers[numreads]);
		} else {
			*fastx_read_list->pointers[numreads] = trimmed_meta_read;
			bytes_printed += fastx_read_fprint(fp, fastx_read_list->pointers[numreads]);
		}
	}

	return bytes_printed;
}

*/

static inline size_t fastx_read_list_trim_and_print_pe(struct fastx_read_list_s * fastx_read_list, FILE * fp[2], int trimN, int trimQ, ssize_t minlength, FILE* discardedfp) {
	if (fastx_read_list->numreads & 0x1)
		return -1;
	size_t numreads;
	size_t bytes_printed = 0;
	for (numreads = 0; numreads < fastx_read_list->numreads; numreads+=2) {
		struct fastx_read_s trimmed_meta_read[2];
		
		fastx_read_trim(trimmed_meta_read, fastx_read_list->pointers[numreads], trimN, trimQ);
		fastx_read_trim(trimmed_meta_read + 1, fastx_read_list->pointers[numreads+1], trimN, trimQ);
		
		size_t minthislen = minlength >= 0 ? minlength : (fastx_read_list->pointers[numreads]->length.seq >> 1);
		
		if (trimmed_meta_read[0].length.seq < minthislen) {
			if (discardedfp) {
				bytes_printed += fastx_read_fprint(discardedfp, fastx_read_list->pointers[numreads]);
			}
		} else {
			*fastx_read_list->pointers[numreads] = trimmed_meta_read[0];
			bytes_printed += fastx_read_fprint(fp[0], fastx_read_list->pointers[numreads]);
		}
		
		minthislen = minlength >= 0 ? minlength : (fastx_read_list->pointers[numreads+1]->length.seq >> 1);
		
		if (trimmed_meta_read[1].length.seq < minthislen) {
			if (discardedfp) {
				bytes_printed += fastx_read_fprint(discardedfp, fastx_read_list->pointers[numreads+1]);
			}
		} else {
			*fastx_read_list->pointers[numreads+1] = trimmed_meta_read[1];
			bytes_printed += fastx_read_fprint(fp[1], fastx_read_list->pointers[numreads+1]);
		}
	}
	
	return bytes_printed;
}

static inline void fastx_read_list_free(struct fastx_read_list_s * fastx_read_list) {
	munmap(fastx_read_list->structs, fastx_read_list->memory_total);
	munmap(fastx_read_list, fastx_read_list->memory_this);
}

// #ifdef __MACH__
// extern char fastx_read_error_message_str[512];
// #else
// extern __thread char fastx_read_error_message_str[512];
// #endif
// 
// #define fastx_read_error_message() fastx_read_error_message_str

#endif /* end of include guard: FASTX_H_GHLHCTHM */
