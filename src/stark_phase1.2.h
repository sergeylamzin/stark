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
#include <stdint.h>
#include "list.h"
#include "stark.h"

#ifndef STARK_PHASE1_2_H_TV4QJMGR
#define STARK_PHASE1_2_H_TV4QJMGR

struct stark_node_phase1_small_s {
	union {
		uint8_t offset_u[2];
		uint16_t relocation_offset_h;
		volatile uint16_t lock;
	} offsets_u;
	// uint8_t flags;
	uint8_t flags;
	uint8_t edges;
	union {
		uint16_t offset_l[2];
		uint32_t relocation_offset_l;
	} offsets_l;
	
	#ifdef DEBUG_P1_SEQ
	char seq[256];
	#endif
};

union stark_node_phase1_small_u {
	struct stark_node_phase1_small_s node;
	volatile uint64_t int64;
};

struct stark_node_phase1_extension_s // __attribute__((packed, aligned(4)))
{
	uint8_t offset_u[4];
	uint16_t offset_l[2][4];
};


#define STARK_PHASE1_MAX_SMALL_NODES_IN_COMPARTMENT 0x100000
#define STARK_PHASE1_MAX_COMPARTMENTS 0x100000
#define STARK_PHASE1_RELOCATION_OFFSET_MAX 0x800000000000

#define STARK_PHASE1_PARENTS_REVERSED 0x80
#define STARK_PHASE1_NODELETE 0x40

// #ifndef HARD_MAXDEPTH
// #define HARD_MAXDEPTH 0xFE
// typedef uint8_t depth_t;
// typedef void * stark_t;
// #endif

struct stark_phase1_s {
	int coverage_t_bytes;
	list_t(struct stark_node_phase1_extension_s) extensions;
	struct {
		union stark_node_phase1_small_u (*all_compartments)[STARK_PHASE1_MAX_SMALL_NODES_IN_COMPARTMENT];
		size_t compartments_modulus;
		// void * coverage_array;
		list_t(char) coverage;
		// list_t(list_t(union stark_node_phase1_small_u)) compartments;
		list_t(struct {uint32_t size;}) compartments;
	} level[HARD_MAXDEPTH];
};


// int stark_node_bucket_mask = 0xFFF;

#define TOKENPASTE2(x, y) x ## y
#define TOKENPASTE(x, y) TOKENPASTE2(x, y)

#define stark_node_phase1_small_get_relocation_offset(_node) (((int64_t)(_node).offsets_u.relocation_offset_h << 32) | (_node).offsets_l.relocation_offset_l)
#define stark_node_phase1_small_set_relocation_offset(_node, _offset) ( \
	  ((_node).offsets_u.relocation_offset_h = (_offset) >> 32)\
	, (_node).offsets_l.relocation_offset_l = (_offset))

#define stark_node_phase1_small_is_locked(_node_u) ((_node_u).node.offsets_u.lock == 0xFFFF)
#define edges_count(v) (((uint64_t)v) * 0x200040008001ULL & 0x111111111111111ULL) % 0xf

#define stark_node_phase1_small_child_exists(_node, _child_direction, _child_transition_char) ((_node).edges & (1 << (((_child_direction) << 2) | (_child_transition_char))))

#define stark_node_phase1_small_write_child(_node, _child_num, _child_direction, _child_transition_char, _child_offset) (\
	(_node).offsets_u.offset_u[_child_num] = (((_child_direction) & 0x1) << 6) \
											| (((_child_transition_char) & 0x3) << 4) \
											| (((_child_offset) >> 16) & 0xF) \
	, \
	((_node).offsets_l.offset_l[_child_num] = _child_offset & 0xFFFF))

#define stark_node_phase1_small_iterate_children(_node, _child_direction, _child_transition_char, _child_offset) \
	int_fast8_t _child_direction = ((_node).offsets_u.offset_u[0] >> 6) & 0x1, \
				_child_transition_char = ((_node).offsets_u.offset_u[0] >> 4) & 0x3 ; \
	int64_t _child_offset = ((int64_t)((_node).offsets_u.offset_u[0] & 0xF) << 16) | (int64_t)(_node).offsets_l.offset_l[0] ; \
	int_fast8_t TOKENPASTE(_stark_node_phase1_small_iterator_num_children_, __LINE__) = edges_count((_node).edges); \
	for (; TOKENPASTE(_stark_node_phase1_small_iterator_num_children_, __LINE__) -- ; \
		  _child_direction = ((_node).offsets_u.offset_u[1] >> 6) & 0x1 \
		, _child_transition_char = ((_node).offsets_u.offset_u[1] >> 4) & 0x3 \
		, _child_offset = ((int64_t)((_node).offsets_u.offset_u[1] & 0xF) << 16) | (int64_t)(_node).offsets_l.offset_l[1] \
	)

#define stark_node_phase1_extension_get_offset(_extension, _child_direction, _child_transition_char) (\
	  (((int64_t)(_extension).offset_u[_child_transition_char] << (16 - ((_child_direction) << 2))) & 0xF0000) \
	| (int64_t)(_extension).offset_l[_child_direction][_child_transition_char])

#define stark_node_phase1_extension_set_offset(_extension, _child_direction, _child_transition_char, _offset) (\
	  ((_extension).offset_u[_child_transition_char] |= ((_offset) & 0xF0000) >> (16 - ((_child_direction) << 2))) \
	, ((_extension).offset_l[_child_direction][_child_transition_char] = (_offset)))


int stark_phase1_insert_one_sequence(void *_target, stark_t* const stark, struct uplevelcache* const uplevel, const depth_t depth);
void stark_phase1_init(struct stark_phase1_s * phase1, int maxK, size_t _max_memory, int coverage_t_bytes);


#endif /* end of include guard: STARK_PHASE1_2_H_TV4QJMGR */

