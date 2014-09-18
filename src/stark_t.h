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


// #include <arpa/inet.h>
#ifdef THREADED_IMPORT
#include <pthread.h>
#endif
#include "main.h"
#include "distributor.h"

#ifndef STARK_T_H
#define STARK_T_H


// #define STARK_FLAG_PALINDROME 0x20
#define STARK_NODE_IS_PALINDROME(node) ( (((node).flags & (STARK_FLAG_PARENT1_REVERSE | STARK_PARENT1_LINK_MASK)) >> 4) == (STARK_FLAG_PARENT2_REVERSE ^ ((node).flags & (STARK_FLAG_PARENT2_REVERSE | STARK_PARENT2_LINK_MASK))) && (node).uplink[0] == (node).uplink[1])

#define STARK_NODE_GET_CHILD_OFFSET(node, child_direction, child_index) ((node).child[child_direction][child_index])
#define STARK_NODE_VALID(node) (node).depth

/*
	starknode_t.flags

	mask 0x80: STARK_FLAG_NODELETE
	
	mask 0x40: STARK_FLAG_PARENT1_REVERSE
	mask 0x30: STARK_PARENT1_LINK 
	
	mask 0x08: STARK_FLAG_EXTRACTED
	
	mask 0x04: STARK_FLAG_PARENT2_REVERSE
	mask 0x03: STARK_PARENT2_LINK

*/

#define STARK_FLAG_NODELETE  0x80
#define STARK_FLAG_EXTRACTED 0x08

// #define STARK_FLAG_PARENT2_REVERSE 0x10

#define FOURSYMBOL_COMPLEMENT(idx) ((idx) ^ 0x3)

#define STARK_FLAG_PARENT1_REVERSE 0x40

#define STARK_PARENT1_LINK_MASK 0x30
#define STARK_PARENT1_LINK_A 0x00
#define STARK_PARENT1_LINK_C 0x10
#define STARK_PARENT1_LINK_G 0x20
#define STARK_PARENT1_LINK_T 0x30


#define STARK_FLAG_PARENT2_REVERSE 0x04

#define STARK_PARENT2_LINK_MASK 0x03
#define STARK_PARENT2_LINK_A 0x00
#define STARK_PARENT2_LINK_C 0x01
#define STARK_PARENT2_LINK_G 0x02
#define STARK_PARENT2_LINK_T 0x03

#define STARK_EDGES_FORWARD_MASK 0xF0
#define STARK_EDGES_REVERSE_MASK 0x0F

#define STARK_CHAIN_FLAG_REVERSE 0x80

extern const char starkinit[];

extern const char foursymbol_index[];

extern const char parent_Link_index[];

extern const char complement_index[];

extern const unsigned char edge_index[2]['U'+1];

extern const unsigned char edgeid_single_index[];

extern const unsigned char edgeid_single_char[];

extern const char parent_Link_to_char[];

extern const unsigned char nucleotides[];

#ifndef STARK_MD5
#define STARK_MD5 0
#endif

#ifdef STARK_VERYLONG
	#define HARD_MAXDEPTH UINT32_MAX
	typedef uint32_t depth_t;
	#define htondepth_t(value) htonl(value)
	#define ntohdepth_t(value) ntohl(value)
#elif STARK_LONG
	#define HARD_MAXDEPTH UINT16_MAX
	typedef uint16_t depth_t;
	#define htondepth_t(value) htons(value)
	#define ntohdepth_t(value) ntohs(value)
#else
	#define HARD_MAXDEPTH UINT8_MAX
	typedef uint8_t depth_t;
	#define htondepth_t(value) value
	#define ntohdepth_t(value) value
#endif

#if OFFSET_L < 64
	#if OFFSET_L < 32
		#if OFFSET_L < 16
			#define OFFSET_T int8_t
			#define htonoffset_t(value) value
			#define ntohoffset_t(value) value
			#define PRI_OFFSET PRId8
			#define OFFSET_MAX INT8_MAX
			#define locked_offset INT8_MIN
		#else
			#define OFFSET_T int16_t
			#define htonoffset_t(value) htons(value)
			#define ntohoffset_t(value) ntohs(value)
			#define PRI_OFFSET PRId16
			#define OFFSET_MAX INT16_MAX
			#define locked_offset INT16_MIN
		#endif
	#else
		#define OFFSET_T int32_t
		#define htonoffset_t(value) htonl(value)
		#define ntohoffset_t(value) ntohl(value)
		#define PRI_OFFSET PRId32
		#define OFFSET_MAX INT32_MAX
		#define locked_offset INT32_MIN
	#endif
#else
	#define OFFSET_T int64_t
	#define PRI_OFFSET PRId64
	#define htonoffset_t(value) htonll(value)
	#define ntohoffset_t(value) ntohll(value)
	#define OFFSET_MAX INT64_MAX
	#define locked_offset INT64_MIN
#endif

#ifdef DEBUG_MAXDEPTH
#define HARD_MAXDEPTH DEBUG_MAXDEPTH
#endif

typedef OFFSET_T offset_t;

#define OFFSET_VALID(offset) ((offset) > 0)
#define OFFSET_CONTAINS_SHARED_COVBERAGE(offset) (offset < 0)

typedef uint32_t coverage_t;

// list_t(coverage_t);

typedef uint8_t flags_t;

#define htoncoverage_t(value) htons(value)
#define ntohcoverage_t(value) ntohs(value)

#ifdef ENCHAIN_SEQUENCES
typedef uint32_t chainnode_offset_t;

#define STARK_CHAINNODE_FLAG_MULTIPREV 0x80
#define STARK_CHAINNODE_FLAG_MULTINEXT 0x40
#define STARK_CHAINNODE_FLAG_MULTIMASK 0xC0
#define STARK_CHAINNODE_FLAG_MULTIFIRST STARK_CHAINNODE_FLAG_MULTINEXT

#define STARK_CHAINNODE_MULTINODE 0x80
#define STARK_CHAINNODE_FLAG_HAS_MASTER 0x40

// #define STARK_CHAINNODE_MULTINODE

typedef struct {
	struct chainmultinode_size_t {
		unsigned char next;
		unsigned char previous;
	} size;
	chainnode_offset_t offset[2];
} chainmultinode_t;

typedef struct chainnode_t {
	unsigned char flags;
	char character;
	coverage_t coverage;
	union {
		struct {
			chainnode_offset_t next;
			chainnode_offset_t previous;
		};
		chainnode_offset_t link[2];
		chainmultinode_t* multinode;
		chainnode_offset_t master;
	};
} chainnode_t;

#define STARK_CHAIN_FLAG_REROUTE 0x80

typedef struct {
	// size_t reroute;
	unsigned char flags;
	chainnode_offset_t start;
	// size_t chain_id;
} chain_t;
#endif

typedef int16_t stark_lock_t;

typedef struct {
	
	// union {
		offset_t child[2][4];
	// };
	
	#ifdef ENCHAIN_SEQUENCES
	struct sequence_s {
		unsigned char flags;
		size_t chain_id;
		union {
			struct {
				chainnode_offset_t start, end;
			};
			chainnode_offset_t endpoint[2];
		};
		// chainnode_offset_t start, end;
		} sequence;
	#endif
	
	// offset_t child[2][4];
	
	offset_t uplink[2];
	
	coverage_t coverage;
	
	uint8_t edges;
	flags_t flags;
	
	depth_t depth;


	#ifdef DEBUG_SEQUENCES
	char sequence[DEBUG_SEQUENCES];
	#endif
	#ifdef UBIGRAPHAPI_H
	struct  {
		vertex_id_t id;
		edge_id_t edge[2];
		char visible;
	} ubigraph;
	#endif
	
	#if HIERARCHIAL_ASSEMBLY
	int64_t hgroup;
	// int32_t contig_position;
	// stark_lock_t lock[2];
	#endif
} starknode_t;

struct coverage_histogram_s {
	intmax_t nodes_counted;
	intmax_t sum_coverage;
	coverage_t max_coverage;
	coverage_t avg_coverage;
	offset_t hist[0];
};

struct coverage_log_histogram_s {
	uintmax_t nodes_counted;
	uintmax_t sum_coverage;
	coverage_t max_coverage;
	coverage_t avg_coverage;
	coverage_t estimated_sample_mean;
	double estimated_standard_deviation;
	int log_mass_hook;
	offset_t loghist[8 * sizeof(uintptr_t)];
};

struct coverage_log_statistics_s {
	uintmax_t nodes_counted;
	uintmax_t sum_coverage;
	coverage_t max_coverage;
	coverage_t avg_coverage;
	
	depth_t maxdepth;
	struct coverage_log_histogram_s hist[0];
};

struct simple_assembly_s {
	size_t length;
	char* sequence;
	coverage_t* coverage;
};


typedef struct {
	starknode_t root;
	starknode_t* level[HARD_MAXDEPTH];
	
	size_t size[HARD_MAXDEPTH];
	size_t deleted[HARD_MAXDEPTH];
	size_t maxsize[HARD_MAXDEPTH];
	size_t memory_allocated[HARD_MAXDEPTH];
	
	#ifdef THREADED
	thread_pool_t thread_pool;
	#endif
	#if defined(THREADED_IMPORT) && defined(STARK_ALLOC_MALLOC)
	pthread_mutex_t mutex[HARD_MAXDEPTH];
	#endif
	
	depth_t maxdepth;
	depth_t longest_palindrome;
	// unsigned short junctions;
	
	// list_t sequences;
	
	#ifdef ENCHAIN_SEQUENCES
	list_t(chain_t) chains;
	list_t(chainnode_t) chainnodes;
	#endif
	
	struct coverage_histogram_s* coverage_histogram;
	
	#ifdef UBIGRAPHAPI_H
	unsigned char level_shown[HARD_MAXDEPTH];
	vertex_id_t ubicount;
	#endif
	
	struct stark_phase1_s * phase1;
	
} stark_t;

// #ifdef THREADED_IMPORT
struct stark_insert_threaded_s {
	stark_t* stark;
	page_distributor_t* distributor;
	size_t length;
	depth_t maxdepth;
	char sequence[0];
};
// #endif


struct starknode_id_s {
	stark_t* stark;
	depth_t depth;
	offset_t offset;
	unsigned char direction;
	char transition_char;
	starknode_t *node_p;
	starknode_t node;
	double score;
};


#endif
