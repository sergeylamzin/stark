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


#ifdef THREADED_IMPORT
#define THREADED
#endif

#if defined(THREADED_IMPORT) && defined(STARK_ALLOC_MALLOC)
#error THREADED_IMPORT and STARK_ALLOC_MALLOC are mutually exclusive. Sorry.
#endif

#ifdef THREADED
#include <pthread.h>
#include "threadpool.h"
#endif

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"
#include "list.h"
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <stdint.h>

#define ___RESTRICT __restrict

#ifdef UBIGRAPH
#include "UbigraphThreaded.h"
#endif

#ifndef STARK_H
#define STARK_H

// #define TRIE_INDEX(i) (foursymbol_index[i])
#define MINALLOC 2//(1L << 28)
#define MAXALLOC (1L << 34)

extern int STARK_MINOVERLAP;
extern const char CONTAINER_DEF[];

// #define STARK_NODE_IS_PALINDROME(node) ((node).flags & STARK_FLAG_PALINDROME)

#define STARK_NODE_UPLINK(node, upnum) ((node).uplink[upnum])
#define STARK_NODE_CHILD(node, direction, childnum) ((node).child[direction][upnum])


#define STARK_VERSION 0
#define STARK_SUBVERSION 0
#define STARK_SUBSUBVERSION 1

#define STARK_VERSION_STRING TOSTRING(STARK_VERSION) "." TOSTRING(STARK_SUBVERSION) "." TOSTRING(STARK_SUBSUBVERSION)

// list_t(char);

#include "stark_t.h"

// struct advancement {
// 	offset_t offset;
// 	depth_t depth;
// 	char direction;
// };

void stark_init(stark_t* const stark);
void stark_free(stark_t* const stark);

void insert_sequence(stark_t* const ___RESTRICT stark, const char* const sequence, const size_t length, depth_t maxdepth);

// __inline__ offset_t create_node(stark_t* const stark, unsigned char level);

// starknode_t* find (starknode_t* node, char* seq);

void stark_print_depth(stark_t* const stark, depth_t depth);
void stark_print(stark_t* const stark, depth_t maxdepth);
// size_t stark_print_node(char* buffer, starknode_t* node);

size_t stark_unambiguos_clean(stark_t* const stark, depth_t minK);

size_t stark_extract_sequence(char* seq, stark_t* stark, depth_t startdepth, offset_t startoffset);

// void reverse_complement_onspot(char* seq, size_t length);

void stark_compact(stark_t* const stark);

// void stark_save(stark_t* const stark, int fd);
// 
// int stark_load(stark_t* const stark, int fd);

int stark_serialize_fp(FILE* fp, stark_t* stark);
int stark_serialize_file(char* filename, stark_t* stark);

int stark_unserialize_file(char* filename, stark_t* stark);
int stark_unserialize_fp(FILE* fp, stark_t* stark);

int stark_find(depth_t* depth_r, offset_t* offset_r, stark_t* stark, char* sequence);

void stark_extract_all(stark_t* const stark);

#ifdef ENCHAIN_SEQUENCES
size_t stark_start_chain(char* seq, stark_t* stark, depth_t startdepth, offset_t startoffset);
#endif

size_t stark_memory_usage(stark_t* const stark);

void stark_assemble_maxcov(struct simple_assembly_s *assembly, stark_t* stark, depth_t startdepth, offset_t startoffset, int32_t minK);

void stark_assemble_single(stark_t* const stark, depth_t minK, depth_t depth, offset_t offset);

void stark_assemble_all(stark_t* const stark, depth_t minK);

size_t stark_print_node_info(char* buffer, starknode_t * node);

struct coverage_histogram_s* stark_statistics(stark_t* const stark, depth_t maxdepth);

// __inline__ void stark_statistics_log_depth(struct coverage_log_histogram_s* target, stark_t* stark, depth_t depth);
struct coverage_log_statistics_s* stark_statistics_log(stark_t* const stark);
void stark_print_log_statistics(FILE* stats_fp, struct coverage_log_statistics_s* log_statistics);

#ifdef THREADED_IMPORT
void stark_insert_sequence_job(struct stark_insert_threaded_s* ___RESTRICT arg);
#endif

#ifdef UBIGRAPHAPI_H
void stark_ubigraph_rebuild(stark_t* const stark);

#endif

static inline int64_t tuple(const int32_t x, const int32_t y) {
	int_fast32_t min, max;
	if (x < y) {
		min = x;
		max = y;
	} else {
		min = y;
		max = x;
	}
	return (max * (max + 1))/2 + min;
	// return x * y;
	//return ((x + y)*(x + y + 1))/2 + y;
}

// #universe for k odd = 4^k / 2 (every k-mer uniquely maps to exactly one reverse complement)
// #universe for k even = 4^k / 2 + 4^(k/2) / 2 (palindromes are unique)
static inline size_t k_mer_universe_size(int_fast32_t k) {
	if (k >= 32)
		return SIZE_MAX;
	size_t universe_size = UINT64_C(1) << ((k << 1) - 1);
	if (!(k & 1))
		universe_size += UINT64_C(1) << (k - 1);
	
	return universe_size;
}


struct uplevelcache {
	int8_t symbol_index;
	int8_t reverse;
	int8_t palindrome;
	// int8_t type;
	// depth_t depth;
	//starknode_t* node_p;
	int32_t compartment;
	// void * node;
	int64_t offset;
};

struct stark_read_coverage_statistics_s {
	size_t count;
	double sum, square_sum;
	size_t hist[1000];
};

struct delcache {
	//starknode_t* node_p;
	offset_t offset;
	intmax_t amount;
	intmax_t shared_coverage;
	depth_t depth;
	int8_t symbol_index;
	int8_t reverse;
	int8_t palindrome;
	int8_t forwardchar;
	int8_t backchar;
};

// #ifdef DEBUG
// extern size_t total_insert_mis;
// extern size_t total_insertions;
// #endif

static inline offset_t stark_get_neighbour(stark_t* const stark, offset_t offset, depth_t depth, int_fast8_t direction, int_fast8_t ch) {
	starknode_t * node = stark->level[depth] + offset;
	
	// char buffer[1024];
	// stark_print_node_info(buffer, node);
	// DEBUG_MSG("starting at node %s", buffer);
	
	starknode_t * parent = stark->level[depth-1] + node->uplink[1 ^ direction];
	// stark_print_node_info(buffer, parent);
	
	// DEBUG_MSG("going to parent %s", buffer);
	
	direction ^= node->flags & (direction ? STARK_FLAG_PARENT1_REVERSE : STARK_FLAG_PARENT2_REVERSE) ? 1 : 0;
	offset_t result = parent->child[direction][ch];
	// DEBUG_MSG("returning parent->child[%d][%d] = %d", direction, ch, result);
	
	
	return result;
}


static inline void get_sequence(char* ___RESTRICT sequence, stark_t* const ___RESTRICT stark, depth_t depth, offset_t offset) {
	#ifndef DEBUG_SEQUENCES
	sequence[depth] = 0;
	
	unsigned char flags, direction = 0;
	
	const starknode_t * ___RESTRICT node; //, *parent;
	
	while (depth) {
		node = stark->level[depth] + offset;
		// offset_t parent_offset 
		offset = node->uplink[direction];
		flags = node->flags;
		// parent = stark->level[--depth] + parent_offset;
	
	
		sequence[--depth] = parent_Link_to_char[flags & (direction ? STARK_PARENT2_LINK_MASK : STARK_PARENT1_LINK_MASK)];
	
		if (flags & (direction ? STARK_FLAG_PARENT2_REVERSE : STARK_FLAG_PARENT1_REVERSE))
			direction ^= 1;

	}

	#else
	memcpy(sequence, stark->level[depth][offset].sequence, depth +1);
	#endif
}


static inline void reverse_complement(char* const ___RESTRICT rev, const char* const ___RESTRICT seq, const size_t length) {
	size_t i,j = length-1;
	// char s,d;
	// int ret = 0;
	for (i = 0; i < length; i++) {
		// s = seq[i];
		// rev[j] = complement_index[s];
		rev[j] = complement_index[seq[i]];
		// d = seq[i] - complement_index[seq[j]];
		// if (d < 0)
		// 	ret = -1; // continue forward
		// else if (d > 0) {
		// 	ret = 1; // reverse
		// }
		j--;
	}
	// return ret;
}

static inline void reverse_complement_onspot(char* const ___RESTRICT seq, const size_t length) {
	size_t i,j = length-1;
	char temp;
	for (i = 0; i < j; i++) {
		temp = seq[j];
		seq[j] = complement_index[seq[i]];
		seq[i] = complement_index[temp];
		j--;
	}
	if (i == j)
		seq[i] = complement_index[seq[i]];
}


static inline void stark_node_init(starknode_t* node) {
	memset(node, 0, sizeof(starknode_t));
}


static inline double stark_link_strength (coverage_t coverage[2]) {
	intmax_t diff = coverage[0];
	diff -= coverage[1];
	intmax_t max = coverage[0];
	if (max < coverage[1])
		max = coverage[1];
	
	if (!max)
		max = 1;
	
	double link_strength = diff;
	// if (diff >= 0) {
	// 	link_strength /= max;
	// } else
	// 	link_strength /= -max;
	
	return link_strength / max;
}

size_t stark_print_node(char* buffer, stark_t* const stark, depth_t depth, offset_t offset);


#include "stark_increment_cache.h"

void stark_substract_coverage_single(stark_t* const stark, depth_t depth, offset_t offset, coverage_t amount);

typedef int (*sequence_process_f)(void *, stark_t * const, struct uplevelcache * const ___RESTRICT, const depth_t);
int stark_read_link_statistics_subf(void *, stark_t* const stark, struct uplevelcache* const ___RESTRICT uplevel, const depth_t depth);
void process_sequence(stark_t* const ___RESTRICT stark, const char* const sequence, const size_t length, depth_t maxdepth, sequence_process_f process_func, void * aux);
void process_sequence_incomplete(stark_t* const ___RESTRICT stark, const char* const sequence, const size_t length, depth_t maxdepth, sequence_process_f process_func, void * aux);

#include "stark_phase1.2.h"

#endif
