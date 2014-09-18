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


#include "stark_phase1.2.h"

static int64_t stark_phase1_get_or_create_offset(struct stark_phase1_s * phase1, depth_t depth, union stark_node_phase1_small_u * backparent_p, int_fast8_t uplevel0reverse, int_fast8_t uplevel0palindrome, int_fast8_t forwardchar_index, int64_t compartment, int64_t *_offset, int flags) {
	union stark_node_phase1_small_u backparent;
	int64_t offset = -1;
	
	for (;;) {
		// this is hopefully an atomic read
		backparent.int64 = backparent_p->int64;
		if (stark_node_phase1_small_is_locked(backparent))
			continue;
		
		if (stark_node_phase1_small_child_exists(backparent.node, uplevel0reverse, forwardchar_index)) {
			if (edges_count(backparent.node.edges) > 2) {
				// this node has extension
				int64_t relocation_offset = stark_node_phase1_small_get_relocation_offset(backparent.node);
				struct stark_node_phase1_extension_s * backparent_extension = phase1->extensions.list + relocation_offset;
				offset = stark_node_phase1_extension_get_offset(*backparent_extension, uplevel0reverse, forwardchar_index);
				#ifndef NOCHECKS
				break;
				#else
				*_offset = offset; 
				return 0;
				#endif
			} else {
				stark_node_phase1_small_iterate_children(backparent.node, child_direction, child_transition_char, child_offset) {
					if (child_direction == uplevel0reverse && child_transition_char == forwardchar_index) {
						offset = child_offset;
						#ifndef NOCHECKS
						break;
						#else
						*_offset = offset; 
						return 0;
						#endif
					}
				}
			}
			#ifndef NOCHECKS
			break;
			#else
			*_offset = offset; 
			return 0;
			#endif
		} // else create the node
		
		// lock
		
		if (!__sync_bool_compare_and_swap(&(backparent_p->node.offsets_u.lock), backparent.node.offsets_u.lock, 0xFFFF))
			continue; // repeat the whole thing
				
		// I obtained the lock
		// Obtain a new lower offset for a new node
		if (*_offset >= 0)
			offset = *_offset;
		else {
			offset = __sync_fetch_and_add(&(phase1->level[depth].compartments.list[compartment & phase1->level[depth].compartments_modulus].size), 1);
			memset(phase1->level[depth].all_compartments[compartment & phase1->level[depth].compartments_modulus] + offset, 0, sizeof(backparent));
		}
		
		if (offset >= STARK_PHASE1_MAX_SMALL_NODES_IN_COMPARTMENT) {
			CRITICAL_ERROR("Offset out of bounds. Too few comparetments.\nRestart with a lower maxK");
		}
		
		phase1->level[depth].all_compartments[compartment & phase1->level[depth].compartments_modulus][offset].node.flags |= flags;
		
		// is the node extended?
		int num_edges;
		if ((num_edges = edges_count(backparent.node.edges)) < (2 - uplevel0palindrome)) {
			// no
			stark_node_phase1_small_write_child(backparent.node, num_edges, uplevel0reverse, forwardchar_index, offset);
			num_edges++;
			if (uplevel0palindrome) {
				stark_node_phase1_small_write_child(backparent.node, num_edges, uplevel0reverse ^ 1, forwardchar_index, offset);
				num_edges++;
			}
		} else {
			// yes
			struct stark_node_phase1_extension_s * backparent_extension;
			
			if (num_edges == 2 || (uplevel0palindrome && num_edges == 1)) { // i need to extend it
				int64_t relocation_offset = __sync_fetch_and_add(&(phase1->extensions.size), 1);
				if (relocation_offset >= phase1->extensions.maxlength || relocation_offset >= STARK_PHASE1_RELOCATION_OFFSET_MAX) {
					CRITICAL_ERROR("Too many extensions");
				}
				backparent_extension = phase1->extensions.list + relocation_offset;
			
				// transfer nodes
				stark_node_phase1_small_iterate_children(backparent.node, child_direction, child_transition_char, child_offset) {
					stark_node_phase1_extension_set_offset(*backparent_extension, child_direction, child_transition_char, child_offset);
				}
				// set relocation
				stark_node_phase1_small_set_relocation_offset(backparent.node, relocation_offset);
			} else { // it is already extended
				backparent_extension = phase1->extensions.list + stark_node_phase1_small_get_relocation_offset(backparent.node);
			}
			
			// set new offset
			stark_node_phase1_extension_set_offset(*backparent_extension, uplevel0reverse, forwardchar_index, offset);
			// num_edges++;
			if (uplevel0palindrome) {
				stark_node_phase1_extension_set_offset(*backparent_extension, uplevel0reverse ^ 1, forwardchar_index, offset);
			}
		}
		
		// enter the new edge
		backparent.node.edges |= 1 << ((uplevel0reverse << 2) | forwardchar_index);
		if (uplevel0palindrome) {
			backparent.node.edges |= 1 << (((uplevel0reverse ^ 1) << 2) | forwardchar_index);
		}
		
		// unlock and write back backparent
		backparent_p->int64 = backparent.int64;
		
		#ifndef NOCHECKS
		break;
		#else
		*_offset = offset; 
		return 1;
		#endif
	} 
	
	if (*_offset >= 0 && *_offset != offset) {
		CRITICAL_ERROR("Data inconsitency");
	}
	
	*_offset = offset; 
	
	return 1;
	
}


int stark_phase1_insert_one_sequence(void *_target, stark_t* const stark, struct uplevelcache* const uplevel, const depth_t depth)
{
	
	struct stark_phase1_s * phase1 = _target;
	
	// {
	// 	char buffer[256];
	// 	int i;
	// 	for (i = 0; i < depth; i++) {
	// 		buffer[i] = nucleotides[uplevel[i].symbol_index];
	// 	}
	// 	DEBUG_MSG("inserting %.*s", i, buffer);
	// }
	
	int_fast8_t palindrome = 1;
	int_fast8_t rindex = 0;

	{
		int_fast8_t d;
		int_fast32_t length = depth & 1;
		length += depth >> 1;
		const struct uplevelcache * seqr = uplevel;
		const struct uplevelcache * sh = uplevel + depth - 1;

		for (;;) {
			if (__builtin_expect(!length, 0)) {
				// palindrome = 1;
				break;
			}
			length--;
			if (__builtin_expect(d = ((seqr++)->symbol_index ^ 0x3) - (sh--)->symbol_index, 1)) {
				rindex = (d < 0);
				palindrome = 0;
				break;
			}
		}
	}
	
	int64_t offset = -1;
	int64_t compartment = uplevel[0].compartment;
	
	if (uplevel[0].offset - uplevel[1].offset <= 0) {
		if (stark_phase1_get_or_create_offset(phase1
										, depth
										, phase1->level[depth-1].all_compartments[0] + uplevel[0].offset
										, uplevel[0].reverse
										, uplevel[0].palindrome
										, uplevel[depth -1].symbol_index
										, compartment
										, &offset
										, 0))
		{
			int flags = 0;
			if (uplevel[0].palindrome | uplevel[1].palindrome) {
				flags = STARK_PHASE1_NODELETE;
			}
			if (uplevel[rindex].offset > uplevel[rindex ^ 1].offset) {
				flags |= 0x40;
			}
			stark_phase1_get_or_create_offset(phase1
											, depth
											, phase1->level[depth-1].all_compartments[0] + uplevel[1].offset
											, uplevel[1].reverse ^ 1
											, uplevel[1].palindrome
											, FOURSYMBOL_COMPLEMENT(uplevel[0].symbol_index)
											, compartment
											, &offset
											, flags);
		}
	} else {
		if (stark_phase1_get_or_create_offset(phase1
										, depth
										, phase1->level[depth-1].all_compartments[0] + uplevel[1].offset
										, uplevel[1].reverse ^ 1
										, uplevel[1].palindrome
										, FOURSYMBOL_COMPLEMENT(uplevel[0].symbol_index)
										, compartment
										, &offset
										, 0))
		{
			int flags = 0;
			if (uplevel[0].palindrome | uplevel[1].palindrome) {
				flags = STARK_PHASE1_NODELETE;
			}
			if (uplevel[rindex].offset > uplevel[rindex ^ 1].offset) {
				flags |= 0x40;
			}
			stark_phase1_get_or_create_offset(phase1
											, depth
											, phase1->level[depth-1].all_compartments[0] + uplevel[0].offset
											, uplevel[0].reverse
											, uplevel[0].palindrome
											, uplevel[depth -1].symbol_index
											, compartment
											, &offset
											, flags);
			}
	}
	
	offset += ((compartment & phase1->level[depth].compartments_modulus) * STARK_PHASE1_MAX_SMALL_NODES_IN_COMPARTMENT);
	
	uplevel->compartment = uplevel[1].offset;
	uplevel->offset =  + offset;
	
	if (phase1->coverage_t_bytes){
		void * coverage_p = phase1->level[depth].coverage.list + (offset * phase1->coverage_t_bytes);
		switch (phase1->coverage_t_bytes) {
			case 1:
			__sync_fetch_and_add((uint8_t *)coverage_p, 1);
			break;
			case 2:
			__sync_fetch_and_add((uint16_t *)coverage_p, 1);
			break;
			case 4:
			__sync_fetch_and_add((uint32_t *)coverage_p, 1);
			break;
			case 8:
			__sync_fetch_and_add((uint64_t *)coverage_p, 1);
			break;
			default:
			break;
		}
	}
	
	#ifdef DEBUG_P1_SEQ
	{
		
		
		int i;
		if (rindex) {
			for (i = 0; i < depth; i++) {
				node->seq[i] = nucleotides[FOURSYMBOL_COMPLEMENT(uplevel[depth - 1 - i].symbol_index)];
			}
		} else {
			for (i = 0; i < depth; i++) {
				node->seq[i] = nucleotides[uplevel[i].symbol_index];
			}
		}
		// node->seq[i] = '\0';
	}
	#endif
	
	;
	uplevel->reverse = rindex;
	uplevel->palindrome = palindrome;
	// uplevel->node = phase1->level[depth].all_compartments[0] + offset;
	
	return 0;
	
}

int stark_phase1_print_node(struct stark_phase1_s * phase1, depth_t depth, int64_t compartment, int64_t offset, void * aux) {
	FILE * fp = stdout;
	fprintf(fp, "%d:[%llu][%llu] -", depth, compartment, offset);
	
	
	struct stark_node_phase1_small_s * node = &phase1->level[depth].all_compartments[compartment][offset].node;
	
	#ifdef DEBUG_P1_SEQ
	char complement[256];
	reverse_complement(complement, node->seq, depth);
	fprintf(fp, " %.*s/%.*s -", (int)depth, node->seq, (int)depth, complement);
	#endif
	
	if (edges_count(node->edges) > 2) {
		struct stark_node_phase1_extension_s * extension = phase1->extensions.list + stark_node_phase1_small_get_relocation_offset(*node);
		
		int i;
		for (i = 0; i < 8; i++) {
			if (node->edges & (1 << i)) {
				int_fast8_t child_direction = i >> 2;
				int_fast8_t transition_char = i & 0x3;
				int64_t child_offset = stark_node_phase1_extension_get_offset(*extension, child_direction, transition_char);
				
				fprintf(fp, " %c:%c:[%llu][%llu]"
					, nucleotides[transition_char]
					, child_direction ? 'r' : 'f'
					, offset & phase1->level[depth].compartments_modulus
					, child_offset);
				
			}
		}
		
	} else {
		stark_node_phase1_small_iterate_children(*node, child_direction, transition_char, child_offset) {
			fprintf(fp, " %c:%c:[%llu][%llu]"
				, nucleotides[transition_char]
				, child_direction ? 'r' : 'f'
				, offset & phase1->level[depth].compartments_modulus
				, child_offset);
		}
	}
	fprintf(fp, "\n");
	
	return 0;
}

void stark_phase1_forall(struct stark_phase1_s * phase1
	, int (*callback)(struct stark_phase1_s *, depth_t, int64_t, int64_t, void *)
	, void * aux)
{
	depth_t depth;
	for (depth = 1; depth < HARD_MAXDEPTH; depth++) {
		int64_t compartment;
		for (compartment = 0; compartment < phase1->level[depth].compartments.size; compartment++) {
			int64_t offset;
			for (offset = 0; offset < phase1->level[depth].compartments.list[compartment].size; offset++) {
				if(callback(phase1, depth, compartment, offset, aux))
					return;
			}
		}
	}
}

void stark_phase1_print_all(struct stark_phase1_s * phase1) {
	stark_phase1_forall(phase1, stark_phase1_print_node, NULL);
}

void stark_phase1_compartment_fill_status(struct stark_phase1_s * phase1)
{
	
	int pages_used[STARK_PHASE1_MAX_SMALL_NODES_IN_COMPARTMENT / sysconf(_SC_PAGESIZE)];
	memset(pages_used, 0, sizeof pages_used);
	
	depth_t depth;
	for (depth = 1; depth < HARD_MAXDEPTH; depth++) {
		int64_t compartment;
		for (compartment = 0; compartment < phase1->level[depth].compartments.size; compartment++) {
			if (phase1->level[depth].compartments.list[compartment].size)
				pages_used[phase1->level[depth].compartments.list[compartment].size * sizeof(phase1->level[depth].all_compartments[0][0]) / sysconf(_SC_PAGESIZE)]++;
		}
	}
	
	int i;
	
	for (i = 0; i < STARK_PHASE1_MAX_SMALL_NODES_IN_COMPARTMENT / sysconf(_SC_PAGESIZE); i++) {
		printf("%d\t%d\n", i, pages_used[i]);
	}
}


#include <sys/mman.h>

static size_t mapmax(int flags, void ** target, size_t size) {
	do {
		void * map;
		map = mmap(
			  NULL
			, size
			, PROT_WRITE | PROT_READ
			, flags
			, -1
			, 0
			);
	
		if (map != MAP_FAILED) {
			*target = map;
			return size;
		}
	} while (size >>= 1);
	
	return 0;
}

void stark_phase1_init(struct stark_phase1_s * phase1, int maxK, size_t _max_memory, int coverage_t_bytes) {
	
	DEBUG_MSG("sizeof(union stark_node_phase1_small_u) = %d", sizeof(union stark_node_phase1_small_u));
	DEBUG_MSG("sizeof(struct stark_node_phase1_extension_s) = %d", sizeof(struct stark_node_phase1_extension_s));
	
	size_t max_memory = 1;
	
	while (max_memory < _max_memory)
		max_memory <<=1;
		
	max_memory <<= 2;
	
	memset(phase1, 0, sizeof(*phase1));
	
	
	
	int depth = 1;
	
	phase1->level[depth].compartments_modulus = 0;
	list_init_size(&phase1->level[depth].compartments, 1);
	list_new_empty(&phase1->level[depth].compartments);
	// list_init_size(&phase1->level[depth].compartments.list[0], 3);
	phase1->level[depth].compartments.list[0].size = 3;
	phase1->level[depth].all_compartments = calloc(sizeof phase1->level[depth].all_compartments[0][0], 3);  // phase1->level[depth].compartments.list[0].list;
	if (coverage_t_bytes) {
		list_init_size(&phase1->level[depth].coverage, coverage_t_bytes * 3);
	}
	
	for (depth = 2; depth <= maxK; depth++) {
		size_t max_kmers_possible = (UINT64_C(1) << (depth * 2));
		
		size_t memory_for_kmers = max_kmers_possible * sizeof(union stark_node_phase1_small_u);
		if (memory_for_kmers < sysconf(_SC_PAGESIZE) && max_kmers_possible <= STARK_PHASE1_MAX_SMALL_NODES_IN_COMPARTMENT) {
			phase1->level[depth].compartments_modulus = 0;
			list_init_size(&phase1->level[depth].compartments, 1);
			list_new_empty(&phase1->level[depth].compartments);
			// list_init_size(&phase1->level[depth].compartments.list[0], max_kmers_possible);
			phase1->level[depth].all_compartments = calloc(sizeof phase1->level[depth].all_compartments[0][0], max_kmers_possible); // phase1->level[depth].compartments.list[0].list;
			if (coverage_t_bytes) {
				list_init_size(&phase1->level[depth].coverage, coverage_t_bytes * max_kmers_possible);
			}
		} else {
			break;
		}
	}
	
	union stark_node_phase1_small_u (*all_compartments)[STARK_PHASE1_MAX_SMALL_NODES_IN_COMPARTMENT];
	size_t memory_got = mapmax(MAP_PRIVATE | MAP_ANON | MAP_NORESERVE, (void **)&all_compartments, max_memory);
	// size_t memory_used = 0;
	
	{
		size_t excess_memory = memory_got % sizeof(*all_compartments);
		excess_memory &= ~(sysconf(_SC_PAGESIZE) -1);
		if (excess_memory) {
			munmap((char*)all_compartments + memory_got - excess_memory, excess_memory);
			memory_got -= excess_memory;
		}
	}
	
	size_t num_compartments = memory_got / sizeof(*all_compartments);
	size_t compartments_used = 0;
	size_t max_num_compartments_per_level = num_compartments / (maxK - depth + 1);
	
	for (; depth <= maxK; depth++) {
		size_t max_kmers_possible = (UINT64_C(1) << (depth * 2));
		
		size_t this_depth_num_compartments = 2;
		while (this_depth_num_compartments <= max_num_compartments_per_level && ((depth > 30) || this_depth_num_compartments * STARK_PHASE1_MAX_SMALL_NODES_IN_COMPARTMENT < max_kmers_possible))
			this_depth_num_compartments <<= 1;
		
		this_depth_num_compartments >>=1;
		
		// DEBUG_MSG("depth = %3d, this_depth_num_compartments = 0x%zx", depth, this_depth_num_compartments);
		
		size_t i;
		phase1->level[depth].all_compartments = all_compartments + compartments_used;
		phase1->level[depth].compartments_modulus = this_depth_num_compartments - 1;
		list_init_size(&phase1->level[depth].compartments, this_depth_num_compartments);
		phase1->level[depth].compartments.size = this_depth_num_compartments;
		memset(phase1->level[depth].compartments.list, 0, sizeof(phase1->level[depth].compartments.list[i]) * this_depth_num_compartments);
		compartments_used += this_depth_num_compartments;
		
		if (coverage_t_bytes) {
			size_t memsize = coverage_t_bytes * this_depth_num_compartments * STARK_PHASE1_MAX_SMALL_NODES_IN_COMPARTMENT;
			list_init_size(&phase1->level[depth].coverage, memsize);
		}
		
		// for (i = 0; i < this_depth_num_compartments; i++) {
		// 	phase1->level[depth].compartments.list[i].maxlength = STARK_PHASE1_MAX_SMALL_NODES_IN_COMPARTMENT;
		// 	phase1->level[depth].compartments.list[i].list = all_compartments[compartments_used++];
		// }
		
		
	}
	
	if (compartments_used > num_compartments) {
		CRITICAL_ERROR("Used more compartments then were allocated");
	}
	
	munmap(all_compartments + compartments_used, memory_got - (compartments_used * sizeof(*all_compartments)));
	
	// extensions
	
	struct stark_node_phase1_extension_s * extensions;
	memory_got = mapmax(MAP_PRIVATE | MAP_ANON | MAP_NORESERVE, (void **)&extensions, STARK_PHASE1_RELOCATION_OFFSET_MAX * sizeof(*extensions));
	{
		size_t excess_memory = memory_got % sizeof(*extensions);
		excess_memory &= ~(sysconf(_SC_PAGESIZE) -1);
		if (excess_memory) {
			munmap((char*)extensions + memory_got - excess_memory, excess_memory);
			memory_got -= excess_memory;
		}
	}
	
	if (!memory_got) {
		CRITICAL_ERROR("Cannot allocate memory for extensions\n");
	}
	
	phase1->extensions.maxlength = memory_got / sizeof(*extensions);
	phase1->extensions.list = extensions;
	
}



















