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


#include "main.h"
#include "stark.h"
// #include "stark_alloc_mmap.h"
#include "stark_navigate.h"
#include "stark_assemble.h"
#include <math.h>

// #ifndef __INTEL_COMPILER
// #if GCC_VERSION < 430 && defined(__x86_64__)
// #define __builtin_bswap32(value) ({__asm__ ("bswapl %0" : "+r" ((uint32_t)value)); value;})
// #define __builtin_bswap64(value) ({__asm__ ("bswapq %0" : "+r" ((uint64_t)value)); value;})
// #endif
// #endif

static inline int stark_node_is_on_surface (starknode_t * node) {
	size_t i;
	for (i = 0; i < 8; i++) {
		if (OFFSET_VALID(node->child[0][i]))
			return 0;
	}
	
	return 1;
}

static inline int stark_redundant_node(stark_t * stark, depth_t depth, offset_t offset) {
	starknode_t * node = stark->level[depth] + offset;
	
	flags_t flags = node->flags;
	if (flags & STARK_FLAG_NODELETE || !stark_node_is_on_surface(node))
		return 0;
	
	size_t i;
	// first check that it has no children
	// for (i = 0; i < 8; i++) {
	// 	if (OFFSET_VALID(node->child[0][i]))
	// 		return 0;
	// }
	
	// if it has no children check if it's parents don't have any other children in the same direction
	starknode_t* forward_parent = stark->level[depth-1] + node->uplink[0];
	starknode_t* reverse_parent = stark->level[depth-1] + node->uplink[1];
	
	
	offset_t * forward_parent_offsets = forward_parent->child[flags & STARK_FLAG_PARENT1_REVERSE ? 1 : 0];
	for (i = 0; i < 4; i++) {
		if (i == ((flags & STARK_PARENT1_LINK_MASK) >> 4))
			continue;
		
		if (OFFSET_VALID(forward_parent_offsets[i]))
			return 0;
	}
	
	offset_t * reverse_parent_offsets = reverse_parent->child[flags & STARK_FLAG_PARENT2_REVERSE ? 0 : 1];
	for (i = 0; i < 4; i++) {
		if (i == (flags & STARK_PARENT2_LINK_MASK))
			continue;
		
		if (OFFSET_VALID(reverse_parent_offsets[i]))
			return 0;
	}
	
	return 1;
	
}

static inline void get_sequence4(char* ___RESTRICT sequence, stark_t* const ___RESTRICT stark, depth_t depth, offset_t offset) {
	sequence[depth] = 0;
	
	unsigned char flags, direction = 0;
	
	const starknode_t * ___RESTRICT node; //, *parent;
	
	while (depth) {
		node = stark->level[depth] + offset;
		// offset_t parent_offset 
		offset = node->uplink[direction];
		flags = node->flags;
		// parent = stark->level[--depth] + parent_offset;
	
		char ch = flags & (direction ? STARK_PARENT2_LINK_MASK : STARK_PARENT1_LINK_MASK);
		if (!direction)
			ch >>= 4;
		
		sequence[--depth] = ch;
	
		if (flags & (direction ? STARK_FLAG_PARENT2_REVERSE : STARK_FLAG_PARENT1_REVERSE))
			direction ^= 1;

	}
}

int stark_assembler_hook (struct starknode_navigator_s * navigator, struct stark_assembler_s * const assembler) {
	// size_t i;
	
	stark_t * stark = assembler->stark;
	
	depth_t depth;
	
	struct starknode_navigator_s * candidates = assembler->candidates;
	
	
	for (depth = HARD_MAXDEPTH - 1; depth > 1 && depth >= assembler->minK; depth--) {
		if (stark->size[depth] > stark->deleted[depth]) {
			offset_t offset;
			flags_t flags;
			int log_mass_hook = assembler->log_statistics->hist[depth].log_mass_hook;
			for (offset = 1; offset <= stark->size[depth]; offset++) {
				starknode_t * node = stark->level[depth] + offset;
				if (STARK_NODE_VALID(*node)) {
					if (stark_redundant_node(stark, depth, offset) || !node->coverage) {
						// node has depleted it's coverage or has become redundant due to deletion of another
						starknode_t* forward_parent = assembler->stark->level[depth-1] + node->uplink[0];
						starknode_t* reverse_parent = assembler->stark->level[depth-1] + node->uplink[1];
						forward_parent->child[flags & STARK_FLAG_PARENT1_REVERSE ? 1 : 0][(flags & STARK_PARENT1_LINK_MASK) >> 4] = 0;
						reverse_parent->child[flags & STARK_FLAG_PARENT2_REVERSE ? 0 : 1][flags & STARK_PARENT2_LINK_MASK] = 0;
						stark_node_init(node);
						__sync_fetch_and_add(stark->deleted + depth, 1);
					}
					else {
						if (!stark_node_is_on_surface(node))
							continue;
						
						int logdiff = fast_log2(node->coverage) - log_mass_hook;
						if (logdiff < 0)
							logdiff = 1 - logdiff;
						if (logdiff == 0) {
							navigator->stark = assembler->stark;
							navigator->depth = depth;
							navigator->offset = offset;
							navigator->node = node;
							return logdiff;
						}

						if (!candidates[logdiff].depth)	{
							candidates[logdiff].stark = assembler->stark;
							candidates[logdiff].depth = depth;
							candidates[logdiff].offset = offset;
							candidates[logdiff].node = node;
						}
					} 
				}
				
				dont_take:
				;
			}
		}
	}
	
	int logdiff = 1;
	for (logdiff = 0; logdiff < 70; logdiff++) {
		if (candidates[logdiff].depth) {
			*navigator = candidates[logdiff];
			return logdiff;
		}
	}
	
	return -1;
}

int stark_assemble_contig(struct starknode_navigator_s * hook, struct stark_assembler_s * const assembler) {
	
	hook->direction = 0;
	
	char * sequence = stark_alloc_pages(0x1000000, MAP_NORESERVE);
	
	get_sequence4(sequence, hook->stark, hook->depth, hook->offset);
	
	int turn = 1;
	
	size_t contig_length = hook->depth;
	
	struct {
		size_t num_kmers_used;
		size_t sum_coverage;
		size_t sum_sq_coverage;
	} statistic = {1, hook->node->coverage, hook->node->coverage * hook->node->coverage};
	
	struct starknode_navigator_s candidate[4], navigator = *hook;
	
	// stark_substract_coverage_single(navigator.stark, navigator.depth, navigator.offset, navigator.node->coverage);
	
	coverage_t hook_initial_coverage = hook->node->coverage;
	coverage_t compare_coverage[2] = {navigator.node->coverage};
	
	stark_substract_coverage_single(navigator.stark, navigator.depth, navigator.offset, navigator.node->coverage);
	
	for (;;) {
		int wentup = -1;

		starknode_t * compare_nodes[2] = {navigator.node};

		// int num_acceptable_extensions = 0;
		int num_possible_extensions = 0;


		int edgemask = 0xf;
		while (!num_possible_extensions && edgemask && navigator.depth > assembler->minK) {
			// num_possible_extensions = 0;
			size_t i;
			for (i = 0; i < 4; i++) {
				offset_t child_offset;
				if (((1 << i) & edgemask) && OFFSET_VALID(child_offset = navigator.node->child[navigator.direction][i])) {
					edgemask &= ~(1 << i);
					candidate[num_possible_extensions] = navigator;
					stark_update_starknode_navigator_down(candidate + num_possible_extensions, child_offset);
					candidate[0].transition_char = i;

					compare_coverage[1] = candidate[num_possible_extensions].node->coverage;

					double link_strength = stark_link_strength(compare_coverage);

					if (link_strength * link_strength < assembler->link_strength)
						num_possible_extensions++;
				}
			}

			wentup++;
			stark_navigate_node_go_up(&navigator);
		}

		// stark_substract_coverage_single(navigator.stark, navigator.depth, navigator.offset, navigator.node->coverage);

		if (num_possible_extensions == 1) {
			sequence[contig_length++] = candidate[0].transition_char;
			navigator = candidate[0];
			statistic.num_kmers_used++;
			statistic.sum_coverage += navigator.node->coverage;
			statistic.sum_sq_coverage += navigator.node->coverage * navigator.node->coverage;
			compare_coverage[0] = navigator.node->coverage;
			stark_substract_coverage_single(navigator.stark, navigator.depth, navigator.offset, navigator.node->coverage);
		} else {
			if (num_possible_extensions > 1) {
				DEBUG_MSG("Encountered ambiguous contig continuation with %d possible extensions.", num_possible_extensions);
				DEBUG_MSG("Feature not yet implemented.");
			} else if (num_possible_extensions == 0) {
				// DEBUG_MSG("Encountered no strong link at kmer, assuming end of contig.");
			}
			
			if (turn) {
				turn = 0;
				#if defined(__x86_64__)
				if (contig_length >= 8) {
					uint64_t * seq64i = (uint64_t *)sequence;
					uint64_t * seq64j = (uint64_t *)(sequence + contig_length);
					while (seq64i < seq64j) {
						seq64j--;
						register uint64_t temp = *seq64i & 0x0303030303030303;
						__asm__ ("bswapq %0" : "+r" (temp));
						register uint64_t temp2 = *seq64j & 0x0303030303030303;
						__asm__ ("bswapq %0" : "+r" (temp2));
						*seq64i = temp2;
						*seq64j = temp;
						seq64i++;
					}
				} else
				#endif
				{
					size_t i,j = contig_length-1;
					for (i = 0; i <= j; i++) {
						char temp = FOURSYMBOL_COMPLEMENT(sequence[i]);
						sequence[i] = FOURSYMBOL_COMPLEMENT(sequence[j]);
						sequence[j] = temp;
						j--;
					}
				}
				
				
				navigator = *hook;
				navigator.direction ^= 1;
				compare_coverage[0] = hook_initial_coverage;
			} else {
				break;
			}
		} 
	}
	
	size_t i;
	for (i = 0; i < contig_length; i++) {
		sequence[i] = parent_Link_to_char[sequence[i]];
	}
	
	sequence[i] = '\0';
	
	if (contig_length > assembler->garbage_length)
	fprintf(contig_length >= assembler->long_length ? assembler->longout : assembler->shortout, 
		">Contig:hook[%d,%d]:length[%zu]:coverage[%zu,%.2f]\n%*s\n"
		, hook->depth
		, hook->offset
		, contig_length
		, statistic.sum_coverage / statistic.num_kmers_used
		, sqrt(((double)(statistic.sum_sq_coverage) - ((double)statistic.sum_coverage / statistic.num_kmers_used)) / statistic.num_kmers_used)
		, (int)contig_length
		, sequence
	);
	else {
		DEBUG_MSG("Trashed garbage contig of length %zu", contig_length);
	}
	
	stark_free_pages(sequence, 0x1000000);
	
	return 0;
	
}

static inline void stark_deep_mark_nodes_extracted_link_strength(stark_t * stark, depth_t depth, starknode_t * node, struct stark_assembler_s * const assembler) {
	// starknode_t * node = stark->level[depth] + offset;
		
	if (depth == 2 || (node->flags & STARK_FLAG_EXTRACTED))
		return;
	
	node->flags |= STARK_FLAG_EXTRACTED;
	
	coverage_t coverage[2] = {node->coverage};
	
	starknode_t * parent_node;
	
	depth--;
	parent_node = stark->level[depth] + node->uplink[0];
	
	coverage[1] = parent_node->coverage;
	double link_strength = stark_link_strength(coverage);
	
	if (link_strength * link_strength < assembler->link_strength)
		stark_deep_mark_nodes_extracted_link_strength(stark, depth, parent_node, assembler);
	
	parent_node = stark->level[depth] + node->uplink[1];
	
	coverage[1] = parent_node->coverage;
	link_strength = stark_link_strength(coverage);
	
	if (link_strength * link_strength < assembler->link_strength)
		stark_deep_mark_nodes_extracted_link_strength(stark, depth, parent_node, assembler);
	
	
}

static inline int stark_node_is_coverage_unambigous(struct starknode_navigator_s * navigator, struct stark_assembler_s * const assembler) {
	starknode_t * node = navigator->node;
	coverage_t coverage[2] = {node->coverage};
	
	int num_forward_matching_children = 0;
	int num_back_matching_children = 0;
	
	size_t i;
	depth_t depth = navigator->depth +1;
	
	for (i = 0; i < 4; i++) {
		offset_t child_offset = node->child[0][i];
		if (OFFSET_VALID(child_offset)) {
			coverage[1] = navigator->stark->level[depth][child_offset].coverage;
			
			double link_strength = stark_link_strength(coverage);
			
			if (link_strength * link_strength < assembler->link_strength) {
				if (num_forward_matching_children++)
					return 0;
			}
		}
		
		child_offset = node->child[1][i];
		if (OFFSET_VALID(child_offset)) {
			coverage[1] = navigator->stark->level[depth][child_offset].coverage;
			
			double link_strength = stark_link_strength(coverage);
			
			if (link_strength * link_strength < assembler->link_strength) {
				if (num_back_matching_children++)
					return 0;
			}
		}
	}
	
	return 1;
}

int stark_assembler_hook2(struct starknode_navigator_s * navigator, struct stark_assembler_s * const assembler) {
	// size_t i;
	
	stark_t * stark = assembler->stark;
	
	depth_t depth;
	
	struct starknode_navigator_s * candidates = assembler->candidates;
	
	
	for (depth = assembler->minK; depth < HARD_MAXDEPTH && stark->level[depth]; depth++) {
		if (stark->size[depth] > stark->deleted[depth]) {
			offset_t offset;
			flags_t flags;
			int log_mass_hook = assembler->log_statistics->hist[depth].log_mass_hook;
			for (offset = 1; offset <= stark->size[depth]; offset++) {
				starknode_t * node = stark->level[depth] + offset;
				if (STARK_NODE_VALID(*node) && !((flags = node->flags) & STARK_FLAG_EXTRACTED)) {
					// if (stark_redundant_node(stark, depth, offset) || !node->coverage) {
					// 	// node has depleted it's coverage or has become redundant due to deletion of another
					// 	starknode_t* forward_parent = assembler->stark->level[depth-1] + node->uplink[0];
					// 	starknode_t* reverse_parent = assembler->stark->level[depth-1] + node->uplink[1];
					// 	forward_parent->child[flags & STARK_FLAG_PARENT1_REVERSE ? 1 : 0][(flags & STARK_PARENT1_LINK_MASK) >> 4] = 0;
					// 	reverse_parent->child[flags & STARK_FLAG_PARENT2_REVERSE ? 0 : 1][flags & STARK_PARENT2_LINK_MASK] = 0;
					// 	stark_node_init(node);
					// 	__sync_fetch_and_add(stark->deleted + depth, 1);
					// }
					// else {
						int logdiff = fast_log2(node->coverage) - log_mass_hook;
						if (logdiff)
							continue;
						
						navigator->stark = assembler->stark;
						navigator->depth = depth;
						navigator->offset = offset;
						navigator->node = node;
						navigator->direction = 0;
						
						if (!stark_node_is_coverage_unambigous(navigator, assembler))
							continue;
						
						// struct starknode_navigator_s downnav = *navigator;
						
						size_t i;
						do {
							coverage_t coverage[2] = {navigator->node->coverage};
							for (i = 0; i < 4; i++) {
								offset_t child_offset = node->child[navigator->direction][i];
								if (OFFSET_VALID(child_offset)) {
									coverage[1] = navigator->stark->level[navigator->depth +1][child_offset].coverage;

									double link_strength = stark_link_strength(coverage);

									if (link_strength * link_strength < assembler->link_strength) {
										stark_update_starknode_navigator_down(navigator, child_offset);
										break;
									}
								}
							}
						} while (i < 4);
						
						return 0;
						// if (!candidates[logdiff].depth)	{
						// 	candidates[logdiff].stark = assembler->stark;
						// 	candidates[logdiff].depth = depth;
						// 	candidates[logdiff].offset = offset;
						// 	candidates[logdiff].node = node;
						// }
					// } 
				}
			}
		}
	}
	
	// int logdiff = 1;
	// for (logdiff = 0; logdiff < 70; logdiff++) {
	// 	if (candidates[logdiff].depth) {
	// 		*navigator = candidates[logdiff];
	// 		return logdiff;
	// 	}
	// }
	
	return -1;
}

int stark_assembler_hook3(struct starknode_navigator_s * navigator, struct stark_assembler_s * const assembler) {
	// size_t i;
	
	stark_t * stark = assembler->stark;
	
	depth_t depth;
	
	struct starknode_navigator_s * candidates = assembler->candidates;
	
	
	for (depth = assembler->maxdepth; depth > 1 && depth >= assembler->minK; depth--) {
		if (stark->size[depth] > stark->deleted[depth]) {
			offset_t offset;
			if (!OFFSET_VALID(offset = assembler->nextoffset))
				offset = 1;
			flags_t flags;
			int log_mass_hook = assembler->log_statistics->hist[depth].log_mass_hook;
			for (; offset <= stark->size[depth]; offset++) {
				starknode_t * node = stark->level[depth] + offset;
				if (STARK_NODE_VALID(*node) && !((flags = node->flags) & STARK_FLAG_EXTRACTED)) {
					// if (stark_redundant_node(stark, depth, offset) || !node->coverage) {
					// 	// node has depleted it's coverage or has become redundant due to deletion of another
					// 	starknode_t* forward_parent = assembler->stark->level[depth-1] + node->uplink[0];
					// 	starknode_t* reverse_parent = assembler->stark->level[depth-1] + node->uplink[1];
					// 	forward_parent->child[flags & STARK_FLAG_PARENT1_REVERSE ? 1 : 0][(flags & STARK_PARENT1_LINK_MASK) >> 4] = 0;
					// 	reverse_parent->child[flags & STARK_FLAG_PARENT2_REVERSE ? 0 : 1][flags & STARK_PARENT2_LINK_MASK] = 0;
					// 	stark_node_init(node);
					// 	__sync_fetch_and_add(stark->deleted + depth, 1);
					// }
					// else {
						// int logdiff = fast_log2(node->coverage) - log_mass_hook;
						// 						if (logdiff)
						// 							continue;
						// 						
						navigator->stark = assembler->stark;
						navigator->depth = depth;
						navigator->offset = offset;
						navigator->node = node;
						navigator->direction = 0;
						
						// if (!stark_node_is_coverage_unambigous(navigator))
						// 	continue;
						
						// struct starknode_navigator_s downnav = *navigator;
						
						// size_t i;
						// do {
						// 	coverage_t coverage[2] = {navigator->node->coverage};
						// 	for (i = 0; i < 4; i++) {
						// 		offset_t child_offset = node->child[navigator->direction][i];
						// 		if (OFFSET_VALID(child_offset)) {
						// 			coverage[1] = navigator->stark->level[navigator->depth +1][child_offset].coverage;
						// 
						// 			double link_strength = stark_link_strength(coverage);
						// 
						// 			if (link_strength * link_strength < 0.1) {
						// 				stark_update_starknode_navigator_down(navigator, child_offset);
						// 				break;
						// 			}
						// 		}
						// 	}
						// } while (i < 4);
						
						return 0;
						// if (!candidates[logdiff].depth)	{
						// 	candidates[logdiff].stark = assembler->stark;
						// 	candidates[logdiff].depth = depth;
						// 	candidates[logdiff].offset = offset;
						// 	candidates[logdiff].node = node;
						// }
					// } 
				}
				assembler->nextoffset++;
			}
		}
		assembler->maxdepth--;
		assembler->nextoffset = 1;
	}
	
	// int logdiff = 1;
	// for (logdiff = 0; logdiff < 70; logdiff++) {
	// 	if (candidates[logdiff].depth) {
	// 		*navigator = candidates[logdiff];
	// 		return logdiff;
	// 	}
	// }
	
	return -1;
}

int stark_assemble_contig2(struct starknode_navigator_s * hook, struct stark_assembler_s * const assembler) {
	
	hook->direction = 0;
	
	char * sequence = stark_alloc_pages(0x1000000, MAP_NORESERVE);
	
	get_sequence4(sequence, hook->stark, hook->depth, hook->offset);
	
	int turn = 1;
	
	size_t contig_length = hook->depth;
	
	struct {
		size_t num_kmers_used;
		size_t sum_coverage;
		size_t sum_sq_coverage;
	} statistic = {1, hook->node->coverage, hook->node->coverage * hook->node->coverage};
	
	struct starknode_navigator_s candidate[4], navigator = *hook;
	
	// stark_substract_coverage_single(navigator.stark, navigator.depth, navigator.offset, navigator.node->coverage);
	
	coverage_t hook_initial_coverage = hook->node->coverage;
	coverage_t compare_coverage[2] = {navigator.node->coverage};
	
	stark_deep_mark_nodes_extracted_link_strength(navigator.stark, navigator.depth, navigator.node, assembler);
	
	int ambiguous_continuation = 0;
	
	for (;;) {
		int wentup = -1;

		starknode_t * compare_nodes[2] = {navigator.node};

		// int num_acceptable_extensions = 0;
		int num_possible_extensions = 0;


		int edgemask = 0xf;
		while (!num_possible_extensions && edgemask && navigator.depth > assembler->minK) {
			// num_possible_extensions = 0;
			size_t i;
			for (i = 0; i < 4; i++) {
				offset_t child_offset;
				if (((1 << i) & edgemask) && OFFSET_VALID(child_offset = navigator.node->child[navigator.direction][i])) {
					candidate[num_possible_extensions] = navigator;
					stark_update_starknode_navigator_down(candidate + num_possible_extensions, child_offset);
					
					if (candidate[num_possible_extensions].node->flags & STARK_FLAG_EXTRACTED) {
						continue;
					}
					
					edgemask &= ~(1 << i);
					
					candidate[0].transition_char = i;

					compare_coverage[1] = candidate[num_possible_extensions].node->coverage;

					double link_strength = stark_link_strength(compare_coverage);

					if (link_strength * link_strength < assembler->link_strength)
						num_possible_extensions++;
				}
			}

			wentup++;
			stark_navigate_node_go_up(&navigator);
		}

		// stark_substract_coverage_single(navigator.stark, navigator.depth, navigator.offset, navigator.node->coverage);

		if (num_possible_extensions == 1) {
			sequence[contig_length++] = candidate[0].transition_char;
			navigator = candidate[0];
			statistic.num_kmers_used++;
			statistic.sum_coverage += navigator.node->coverage;
			statistic.sum_sq_coverage += navigator.node->coverage * navigator.node->coverage;
			compare_coverage[0] = navigator.node->coverage;
			stark_deep_mark_nodes_extracted_link_strength(navigator.stark, navigator.depth, navigator.node, assembler);
			// stark_substract_coverage_single(navigator.stark, navigator.depth, navigator.offset, navigator.node->coverage);
		} else {
			if (num_possible_extensions > 1) {
				ambiguous_continuation = 1;
				DEBUG_MSG("Encountered ambiguous contig continuation with %d possible extensions.", num_possible_extensions);
				DEBUG_MSG("Feature not yet implemented.");
				if (compare_coverage[0] > 1000) {
					fprintf(stderr, "Big contig with ambiguous continuation\n");
				}
			} else if (num_possible_extensions == 0) {
				// DEBUG_MSG("Encountered no strong link at kmer, assuming end of contig.");
			}
			
			if (turn) {
				turn = 0;
				#if defined(__x86_64__)
				if (contig_length >= 8) {
					uint64_t * seq64i = (uint64_t *)sequence;
					uint64_t * seq64j = (uint64_t *)(sequence + contig_length);
					while (seq64i < seq64j) {
						seq64j--;
						register uint64_t temp = *seq64i & 0x0303030303030303;
						__asm__ ("bswapq %0" : "+r" (temp));
						register uint64_t temp2 = *seq64j & 0x0303030303030303;
						__asm__ ("bswapq %0" : "+r" (temp2));
						*seq64i = temp2;
						*seq64j = temp;
						seq64i++;
					}
				} else
				#endif
				{
					size_t i,j = contig_length-1;
					for (i = 0; i <= j; i++) {
						char temp = FOURSYMBOL_COMPLEMENT(sequence[i]);
						sequence[i] = FOURSYMBOL_COMPLEMENT(sequence[j]);
						sequence[j] = temp;
						j--;
					}
				}
				
				
				navigator = *hook;
				navigator.direction ^= 1;
				compare_coverage[0] = hook_initial_coverage;
			} else {
				break;
			}
		} 
	}
	
	size_t i;
	for (i = 0; i < contig_length; i++) {
		sequence[i] = parent_Link_to_char[sequence[i]];
	}
	
	sequence[i] = '\0';
	
	if (contig_length > assembler->garbage_length)
	fprintf(contig_length >= assembler->long_length ? assembler->longout : assembler->shortout, 
		">Contig:hook[%d,%d]:length[%zu]:coverage[%zu,%.2f]\n%*s\n"
		, hook->depth
		, hook->offset
		, contig_length
		, statistic.sum_coverage / statistic.num_kmers_used
		, sqrt(((double)(statistic.sum_sq_coverage) - ((double)statistic.sum_coverage / statistic.num_kmers_used)) / statistic.num_kmers_used)
		, (int)contig_length
		, sequence
	);
	else {
		DEBUG_MSG("Trashed garbage contig of length %zu", contig_length);
	}
	
	stark_free_pages(sequence, 0x1000000);
	
	return 0;
	
}

struct stark_fragment_s {
	struct stark_assembler_s * assembler;
	
	struct {
		struct starknode_navigator_s nav;
		struct {
			int_fast8_t ch;
			int_fast8_t level;
			double link_strength;
			struct stark_fragment_s * fragment;
		} continuation[4];
	} endpoint[2];
	
	struct {
		size_t surface;
		size_t subsize;
		size_t sequence;
	} size;
};

static inline signed char link_strength_level(double link_strength, double reference_strength) {
	signed char level = 1;
	while (link_strength >= reference_strength && level < INT8_MAX ) {
		link_strength *= link_strength;
		level++;
	}
	
	return level;
}

int stark_assemble_fragment1(struct stark_fragment_s * fragment, struct starknode_navigator_s * hook, struct stark_assembler_s * const assembler) {
	
	// hook->direction = 0;
	
	// char * sequence = stark_alloc_pages(0x1000000, MAP_NORESERVE);
	
	// get_sequence4(sequence, hook->stark, hook->depth, hook->offset);
	
	stark_t * stark = assembler->stark;
	
	int turn = 1;
	
	struct stark_fragment_s fragment = {.assembler = assembler};
	
	// size_t contig_length = hook->depth;
	
	struct {
		size_t num_kmers_used;
		size_t sum_coverage;
		size_t sum_sq_coverage;
	} statistic = {1, hook->node->coverage, hook->node->coverage * hook->node->coverage};
	
	struct starknode_navigator_s candidate[4], navigator, current = *hook;
	
	signed char candidate_strength_level[4] = {INT8_MAX, INT8_MAX, INT8_MAX, INT8_MAX};
	
	// stark_substract_coverage_single(navigator.stark, navigator.depth, navigator.offset, navigator.node->coverage);
	
	coverage_t hook_initial_coverage = hook->node->coverage;
	coverage_t compare_coverage[3] = {navigator.node->coverage, 0, 0};
	
	stark_deep_mark_nodes_extracted_link_strength(navigator.stark, navigator.depth, navigator.node, assembler);
	
	int ambiguous_continuation = 0;
	
	for (;;) {
		int wentup = -1;

		// starknode_t * compare_nodes[2] = {navigator.node};

		// int num_acceptable_extensions = 0;
		int num_possible_extensions = 0;


		int edgemask = 0xf;
		
		navigator =  current;
		compare_coverage[0] = navigator.node->coverage;
		
		// depth_t depth = navigator.depth;
		
		while (navigator.depth > assembler->minK) {
			size_t i;
			for (i = 0; i < 4; i++) {
				offset_t child_offset;
				if (OFFSET_VALID(child_offset = navigator.node->child[navigator.direction][i])) {
					
					struct starknode_navigator_s child_nav = navigator;
					stark_update_starknode_navigator_down(&child_nav, child_offset);
					
					compare_coverage[1] = child_nav.node->coverage;
					
					double link_strength = stark_link_strength(compare_coverage);
					signed char strength_level;
					
					if ((strength_level = link_strength_level(link_strength * link_strength, assembler->link_strength)) < candidate_strength_level[i]
					 	&& strength_level < 4 ) {
						candidate[i] = child_nav;
						candidate_strength_level[i] = strength_level;
					}
					
				}
			}
			
			stark_navigate_node_go_up(&current);
			compare_coverage[1] = current.node->coverage;
			if (link_strength * link_strength < assembler->link_strength) {
				navigator =  current;
				compare_coverage[0] = compare_coverage[1];
			}
		}
		
		int num_direct_extensions = 0;
		int_fast8_t direct_extension_ch = -1;
		size_t i;
		for (i = 0; i < 4; i++) {
			if (candidate_strength_level[i] == 1){
				num_direct_extensions++;
				direct_extension_ch = i;
			}
		}
		
		if (num_direct_extensions == 1) {
			current = candidate[direct_extension_ch];
			stark_deep_mark_nodes_extracted_link_strength(current.stark, current.depth, current.node, assembler);
			continue;
		}
		
		
		if (num_direct_extensions > 1) {
			
		}

		// stark_substract_coverage_single(navigator.stark, navigator.depth, navigator.offset, navigator.node->coverage);

		if (num_possible_extensions == 1) {
			sequence[contig_length++] = candidate[0].transition_char;
			navigator = candidate[0];
			statistic.num_kmers_used++;
			statistic.sum_coverage += navigator.node->coverage;
			statistic.sum_sq_coverage += navigator.node->coverage * navigator.node->coverage;
			compare_coverage[0] = navigator.node->coverage;
			stark_deep_mark_nodes_extracted_link_strength(navigator.stark, navigator.depth, navigator.node, assembler);
			// stark_substract_coverage_single(navigator.stark, navigator.depth, navigator.offset, navigator.node->coverage);
		} else {
			if (num_possible_extensions > 1) {
				ambiguous_continuation = 1;
				DEBUG_MSG("Encountered ambiguous contig continuation with %d possible extensions.", num_possible_extensions);
				DEBUG_MSG("Feature not yet implemented.");
				if (compare_coverage[0] > 1000) {
					fprintf(stderr, "Big contig with ambiguous continuation\n");
				}
			} else if (num_possible_extensions == 0) {
				// DEBUG_MSG("Encountered no strong link at kmer, assuming end of contig.");
			}
			
			if (turn) {
				turn = 0;
				#if defined(__x86_64__)
				if (contig_length >= 8) {
					uint64_t * seq64i = (uint64_t *)sequence;
					uint64_t * seq64j = (uint64_t *)(sequence + contig_length);
					while (seq64i < seq64j) {
						seq64j--;
						register uint64_t temp = *seq64i & 0x0303030303030303;
						__asm__ ("bswapq %0" : "+r" (temp));
						register uint64_t temp2 = *seq64j & 0x0303030303030303;
						__asm__ ("bswapq %0" : "+r" (temp2));
						*seq64i = temp2;
						*seq64j = temp;
						seq64i++;
					}
				} else
				#endif
				{
					size_t i,j = contig_length-1;
					for (i = 0; i <= j; i++) {
						char temp = FOURSYMBOL_COMPLEMENT(sequence[i]);
						sequence[i] = FOURSYMBOL_COMPLEMENT(sequence[j]);
						sequence[j] = temp;
						j--;
					}
				}
				
				
				navigator = *hook;
				navigator.direction ^= 1;
				compare_coverage[0] = hook_initial_coverage;
			} else {
				break;
			}
		} 
	}
	
	size_t i;
	for (i = 0; i < contig_length; i++) {
		sequence[i] = parent_Link_to_char[sequence[i]];
	}
	
	sequence[i] = '\0';
	
	if (contig_length > assembler->garbage_length)
	fprintf(contig_length >= assembler->long_length ? assembler->longout : assembler->shortout, 
		">Contig:hook[%d,%d]:length[%zu]:coverage[%zu,%.2f]\n%*s\n"
		, hook->depth
		, hook->offset
		, contig_length
		, statistic.sum_coverage / statistic.num_kmers_used
		, sqrt(((double)(statistic.sum_sq_coverage) - ((double)statistic.sum_coverage / statistic.num_kmers_used)) / statistic.num_kmers_used)
		, (int)contig_length
		, sequence
	);
	else {
		DEBUG_MSG("Trashed garbage contig of length %zu", contig_length);
	}
	
	stark_free_pages(sequence, 0x1000000);
	
	return 0;
	
}

struct stark_component_s {
	int id;
	
};

struct stark_assemble_component1_builder {
	struct starknode_navigator_s candidate[4], current;
	int numcandidates;
	double candidate_link_strength[4];
};

static inline void stark_assemble_component1_get_candidates(struct stark_assemble_component1_builder * builder, struct stark_assembler_s * const assembler) {
	struct starknode_navigator_s navigator = builder->current;
	int edgemask = 0xf;
	builder->numcandidates = 0;
	coverage_t compare_coverage[2] = {navigator.node->coverage, 0};
	while (edgemask && navigator.depth > assembler->minK) {
		// num_possible_extensions = 0;
		size_t i;
		for (i = 0; i < 4; i++) {
			offset_t child_offset;
			if (((1 << i) & edgemask) && OFFSET_VALID(child_offset = navigator.node->child[navigator.direction][i])) {
				builder->candidate[builder->numcandidates] = navigator;
				stark_update_starknode_navigator_down(builder->candidate + builder->numcandidates, child_offset);
				
				if (builder->candidate[builder->numcandidates].node->flags & STARK_FLAG_EXTRACTED) {
					continue;
				}
				
				edgemask &= ~(1 << i);
				
				builder->candidate[builder->numcandidates].transition_char = i;

				compare_coverage[1] = builder->candidate[builder->numcandidates].node->coverage;

				builder->candidate_link_strength[builder->numcandidates] = stark_link_strength(compare_coverage);

				builder->numcandidates++;
			}
		}

		// wentup++;
		stark_navigate_node_go_up(&navigator);
	}
}

void stark_assemble_component1(struct stark_component_s * component, struct starknode_navigator_s * hook, struct stark_assembler_s * const assembler) {
	stark_t * stark = assembler->stark;
	
	double link_strength_cap = 0;
	
	
	struct starknode_navigator_s navigator;
	
	struct stark_assemble_component1_builder builder[2];
	
	memset(builder, 0, sizeof(builder));
	
	builder[0].current = *hook;
	builder[1].current = *hook;
	
	stark_assemble_component1_get_candidates(builder);
	stark_assemble_component1_get_candidates(builder +1);
	
}


void stark_assembler_test(stark_t * stark, struct coverage_log_statistics_s* log_statistics, double link_strength) {
	struct stark_assembler_s assembler;
	memset(&assembler, 0, sizeof(assembler));
	assembler.stark = stark;
	assembler.minK = 30;
	assembler.log_statistics = log_statistics;
	assembler.garbage_length = 150;
	assembler.long_length = 500;
	assembler.longout = stdout;
	assembler.shortout = fopen("/dev/shm/shortcontigs.fa", "w");
	assembler.maxdepth = HARD_MAXDEPTH -1;
	assembler.nextoffset = 0;
	assembler.link_strength = link_strength * link_strength;
	
	struct starknode_navigator_s hook;
	
	memset(&hook, 0 , sizeof(hook));
	
	int hook_status;
	while ((hook_status = stark_assembler_hook3(&hook, &assembler)) >= 0) {
		
		#ifdef DEBUG
		// char buffer[0x1000];
		// get_sequence(buffer, hook.stark, hook.depth, hook.offset);
		// stark_print_node(buffer, hook.stark, hook.depth, hook.offset);
		
		DEBUG_MSG("Obtained contig hook at [%d, %d] with status %d"
			, hook.depth, hook.offset
			, hook_status
			// , hook.depth
			// , buffer
		);
		#endif
		
		stark_assemble_contig2(&hook, &assembler);
	}
	
}










