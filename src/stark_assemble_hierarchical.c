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


#ifdef HIERARCHIAL_ASSEMBLY
#include "stark_assemble_hierarchical.h"
#include <stdint.h>
#include <math.h>

#ifdef __clang__
// clang does not have support for openmp at the moment
#define omp_get_thread_num() 0
#define omp_get_team_size(x) 1
#define omp_get_level() 0
#define omp_get_num_threads() 1
#else
#include <omp.h>
#endif


// #if !defined(__MACH__) 
static hashmap_uint64_t thread_seen_map = {0};
static stark_navigate_dfs_stack_t thread_dfs_stack = {0};

#pragma omp threadprivate(thread_seen_map, thread_dfs_stack)

// #else
// #include <omp.h>
// static hashmap_uint64_t * mach_seen_map = NULL;
// #define thread_seen_map (mach_seen_map || mach_seen_map = calloc(omp_get_num_threads(), sizeof(*mach_seen_map)), mach_seen_map[omp_thread_num()])
// static hashmap_uint64_t thread_seen_map = {0};
// static stark_navigate_dfs_stack_t thread_dfs_stack = {0};
// #endif

#define navigator_mark_seen(_map, _nav) map_insert(_map, (uint64_t)(_nav).node ^ ((uint64_t)(_nav).direction))


struct stark_hierarchical_assembler_grow_group_dfs_cb_s {
	struct stark_hierarchical_assembler_s * hierarchical_assembler;
	struct stark_hierarchical_assembler_group_s * group;
	double current_link_stregth_cutoff;
	hirarchial_assembler_group_t group_orig;
	// struct starknode_navigator_s hook_orig;
	list_t(struct starknode_navigator_s) todo_nodes;
};



static inline void stark_hierarchical_assembler_init_add_node_to_group(struct starknode_navigator_s * nav, struct stark_hierarchical_assembler_group_s * group, double link_strength) {
	nav->node->hgroup = group->group_id;
	
	group->num_nodes++;
	if (group->internal_link_strength < link_strength) {
		group->internal_link_strength = link_strength;
	}

	group->mass += nav->node->coverage;
	
	if (!group->endpoint->depth || group->endpoint->depth < nav->depth) {
		*group->endpoint = *nav;
	}
	
}

static inline void stark_hierarchical_assembler_recycle_group(struct stark_hierarchical_assembler_s * hierarchical_assembler, struct stark_hierarchical_assembler_group_s * group) {
	hirarchial_assembler_group_t gid = group->group_id;
	
	#pragma omp critical (hierarchical_assembler)
	memset(group, 0, sizeof(struct stark_hierarchical_assembler_group_s));
	{
		list_insert(&hierarchical_assembler->recycle_groups, gid);
	}
}

static inline hirarchial_assembler_group_t stark_hierarchical_assembler_new_group(struct stark_hierarchical_assembler_s * hierarchical_assembler) {
	hirarchial_assembler_group_t gid;
	#pragma omp critical (hierarchical_assembler)
	{
		if (hierarchical_assembler->recycle_groups.size) {
			gid = hierarchical_assembler->recycle_groups.list[--(hierarchical_assembler->recycle_groups.size)];

		} else {
			gid = list_new_empty(&(hierarchical_assembler->groups));
		}
	}
	return gid;
}


size_t stark_hierarchical_assembler_print_group(FILE * fp, struct stark_hierarchical_assembler_group_s * group) {
	size_t bytes_printed = 0;
	
	bytes_printed += fprintf(fp, "Group %llu %s - num_nodes = %zu, mass = %.0lf, internal_link_strength = %0.4lf, contig length = %zu\n"
		, group->group_id
		, group->contains_palindrome ? "(palindromic)" : ""
		, group->num_nodes
		, group->mass
		, group->internal_link_strength
		, group->contig_length
		);
	
	char buffer[0x400];
	if (group->endpoint[0].depth) {
		stark_print_node(buffer, group->endpoint[0].stark, group->endpoint[0].depth, group->endpoint[0].offset);
		bytes_printed += fprintf(fp, "\t Endpoint 0: %s\n", buffer);
	}
	if (group->endpoint[1].depth) {
		stark_print_node(buffer, group->endpoint[1].stark, group->endpoint[1].depth, group->endpoint[1].offset);
		bytes_printed += fprintf(fp, "\t Endpoint 1: %s\n", buffer);
	}
	
	// int j;
	// for (j = 0; j < group->neighbours.size; j++) {
	// 	bytes_printed += fprintf(fp, "\tNeighbour %llu, link_strength = %0.4lf\n"
	// 		, group->neighbours.list[j].id
	// 		, group->neighbours.list[j].link_strength
	// 		);
	// }
	
	return bytes_printed;
}

size_t stark_hierarchical_assembler_print_neighbours(FILE * fp, struct stark_hierarchical_assembler_s * hierarchical_assembler) {
	size_t bytes_printed = 0;
	
	int j;
	for (j = 0; j < hierarchical_assembler->neighbours.size; j++) {
		bytes_printed += fprintf(fp, "\tNeighbours %llu and %llu, link_strength = %0.4lf\n"
			, hierarchical_assembler->neighbours.list[j].ids[0]
			, hierarchical_assembler->neighbours.list[j].ids[0]
			, hierarchical_assembler->neighbours.list[j].link_strength
			);
	}
	
	return bytes_printed;
}

ssize_t stark_print_all_groups(FILE * fp, struct stark_hierarchical_assembler_s * hierarchical_assembler) {
	size_t i;
	for (i = 1; i < hierarchical_assembler->groups.size; i++) {
		struct stark_hierarchical_assembler_group_s *group = hierarchical_assembler->groups.list + i;
		
		if (group->group_id)		
			stark_hierarchical_assembler_print_group(fp, group);
		else {
			DEBUG_MSG("Group %zu was merged with group %zu", i, group->num_nodes);
		}
	}
	
	return 0;
}

typedef list_t(char) char_list_t;


ssize_t stark_hierarchical_assembler_run_contig_choice(
	struct starknode_navigator_s * _navigator, depth_t mink, hirarchial_assembler_group_t hgroup, size_t hgroup_alt, // input
	struct starknode_navigator_s * endpoint, // output
	int  (*choice_callback) (struct starknode_navigator_s * from, struct starknode_navigator_s * candidates, int numcandidates, void * aux),
	void (*callback_taken) (struct starknode_navigator_s *, int, void *), void * cb_aux // callback
	) {
	
	struct starknode_navigator_s navigator = *_navigator;
	// if (endpoint)
	// 	*endpoint = *_navigator;
	
	// {
	// 	char buffer[0x400];
	// 	stark_print_node(buffer, _navigator->stark, _navigator->depth, _navigator->offset);
	// 	DEBUG_MSG("Hooked node \n%c - %s",_navigator->direction ? 'r' : 'f' , buffer);
	// }
	
	// find endpoint 1
	
	ssize_t num_transitions = navigator.depth;
	// map_zero(&thread_seen_map);
	struct starknode_navigator_s current, previous = navigator, candidate = {0}; // = *_navigator;
	offset_t child_offset;
	// stark_lock_t lock = navigator.node->lock[1] ^ 1;
	
	int breaknext = 0;
	int edges = 0xF;
	if (navigator_mark_seen(&thread_seen_map, navigator)) {
		return -num_transitions;
		// goto stark_hierarchical_assembler_run_contig_loop_error;
	}

	// navigator = current;
	for (;;) {
		
		
		for (;;) {
			candidate.stark = NULL;
			int ch;
			
			struct starknode_navigator_s candidates[4];
			int numcandidates = 0;
			
			for (ch = 0; ch < 4; ch++, edges <<= 1) {
				if ((edges & 0x8) && OFFSET_VALID(child_offset = navigator.node->child[navigator.direction][ch])) {
					// extensions[num_extensions] = current;
					candidates[numcandidates] = navigator;
					stark_update_starknode_navigator_down(candidates + numcandidates, child_offset);
					candidates[numcandidates].transition_char = ch;
					numcandidates++;
					// choice_callback(  &previous  // from
					// 					, &current   // to
					// 					, &candidate // previous_candidate
					// 					, cb_aux);
				}
			}
			
			int approved_candidate;
			
			if (numcandidates && (approved_candidate = choice_callback(&previous, candidates, numcandidates, cb_aux)) >= 0) {
				
			}
			
			if (candidate.stark) {
				if (navigator_mark_seen(&thread_seen_map, candidate)) {
					return -num_transitions;
					// goto stark_hierarchical_assembler_run_contig_loop_error;
				}
			
				num_transitions++;
				if (callback_taken)
					callback_taken(&candidate, ch, cb_aux);
				
				previous = candidate;
				navigator = candidate;
				// debug_print_node(navigator.stark, navigator.depth, navigator.offset);
				edges = navigator.node->edges;
				if (!navigator.direction)
					edges >>= 4;
				
				breaknext = 0;
				if (endpoint)
					*endpoint = navigator;
				
				continue;
			}
			break;
		} 
		
		if (__builtin_expect(breaknext // || STARK_NODE_IS_PALINDROME(navigator.node[0])
			, 0)
			) {
			// TODO palindrome handling
			// for now just abort
			// if (endpoint)
			// 	*endpoint = navigator;

			return num_transitions;
		}
		
		edges = navigator.node->edges;
		if (!navigator.direction)
			edges >>= 4;
		
		// current = navigator;
		stark_navigate_node_go_up(&navigator);
		// {
		// 	char buffer[0x400];
		// 	stark_print_node(buffer, navigator.stark, navigator.depth, navigator.offset);
		// 	DEBUG_MSG("Went up to node \n%c - %s",navigator.direction ? 'r' : 'f' , buffer);
		// }
		
		if (choice_callback(&previous, &navigator, 0, cb_aux) < 0)
			breaknext = 1;
		
		if (navigator_mark_seen(&thread_seen_map, navigator)) {
			return -num_transitions;
			// goto stark_hierarchical_assembler_run_contig_loop_error;
		}
		
		
	}
	
	// stark_hierarchical_assembler_run_contig_loop_error:
	// {
	// 	char buffer[0x400];
	// 	// stark_print_node(buffer, navigator.stark, navigator.depth, navigator.offset);
	// 	// DEBUG_MSG("Ran into a loop in hgroup %llu at node\n%s", hgroup, buffer);
	// 	// list_free(*seq);
	// 	return -num_transitions;
	// }
}



ssize_t stark_hierarchical_assembler_run_contig(
	struct starknode_navigator_s * _navigator, depth_t mink, hirarchial_assembler_group_t hgroup, size_t hgroup_alt, // inpit
	struct starknode_navigator_s * endpoint, // output
	void (*callback) (struct starknode_navigator_s *, int, void *), void * cb_aux // callback
	) {
	
	struct starknode_navigator_s navigator = *_navigator;
	// if (endpoint)
	// 	*endpoint = *_navigator;
	
	// {
	// 	char buffer[0x400];
	// 	stark_print_node(buffer, _navigator->stark, _navigator->depth, _navigator->offset);
	// 	DEBUG_MSG("Hooked node \n%c - %s",_navigator->direction ? 'r' : 'f' , buffer);
	// }
	
	// find endpoint 1
	
	ssize_t num_transitions = navigator.depth;
	// map_zero(&thread_seen_map);
	struct starknode_navigator_s current; // = *_navigator;
	offset_t child_offset;
	// stark_lock_t lock = navigator.node->lock[1] ^ 1;
	
	int breaknext = 0;
	int edges = 0xF;
	if (navigator_mark_seen(&thread_seen_map, navigator)) {
		goto stark_hierarchical_assembler_run_contig_loop_error;
	}

	// navigator = current;
	for (;;) {
		
		
		{
			int ch;
			for (ch = 0; ch < 4; ) {
				if ((edges & 0x8) && OFFSET_VALID(child_offset = navigator.node->child[navigator.direction][ch])) {
					// extensions[num_extensions] = current;
					current = navigator;
					stark_update_starknode_navigator_down(&current, child_offset);
			
					if (current.node->hgroup == hgroup || current.node->hgroup == hgroup_alt) {
					
						// {
						// 	char buffer[0x400];
						// 	stark_print_node(buffer, navigator.stark, current.depth, current.offset);
						// 	DEBUG_MSG("Went down to node \n%c - %s",current.direction ? 'r' : 'f' , buffer);
						// }
					
						if (navigator_mark_seen(&thread_seen_map, current)) {
							goto stark_hierarchical_assembler_run_contig_loop_error;
						}
					
						num_transitions++;
						if (callback)
							callback(&current, ch, cb_aux);
					
						ch = 0;
						navigator = current;
						// debug_print_node(navigator.stark, navigator.depth, navigator.offset);
						edges = navigator.node->edges;
						if (!navigator.direction)
							edges >>= 4;
						
						breaknext = 0;
						if (endpoint)
							*endpoint = navigator;
							
						continue;
					}
			
				}
				ch++;
				edges <<= 1;
			}
		}
		
		// if (breaknext) {
		// 	// if (endpoint)
		// 	// 	*endpoint = navigator;
		// 
		// 	return num_transitions;
		// }
		
		if (__builtin_expect(breaknext || STARK_NODE_IS_PALINDROME(navigator.node[0]), 0)) {
			// TODO palindrome handling
			// for now just abort
			// if (endpoint)
			// 	*endpoint = navigator;

			return num_transitions;
		}
		
		edges = navigator.node->edges;
		if (!navigator.direction)
			edges >>= 4;
		
		// current = navigator;
		stark_navigate_node_go_up(&navigator);
		// {
		// 	char buffer[0x400];
		// 	stark_print_node(buffer, navigator.stark, navigator.depth, navigator.offset);
		// 	DEBUG_MSG("Went up to node \n%c - %s",navigator.direction ? 'r' : 'f' , buffer);
		// }
		// if we walked into a palindrome
		
		
		if (!((navigator.node->hgroup == hgroup || navigator.node->hgroup == hgroup_alt) && navigator.depth >= mink)) {
			// try neighbours first
			breaknext = 1;
			// goto stark_hierarchical_assembler_run_contig_children;
			
		} else {
			if (navigator_mark_seen(&thread_seen_map, navigator)) {
				goto stark_hierarchical_assembler_run_contig_loop_error;
			}

			// navigator = current;
		}
		
		
	}
	
	stark_hierarchical_assembler_run_contig_loop_error:
	{
		char buffer[0x400];
		// stark_print_node(buffer, navigator.stark, navigator.depth, navigator.offset);
		// DEBUG_MSG("Ran into a loop in hgroup %llu at node\n%s", hgroup, buffer);
		// list_free(*seq);
		return -num_transitions;
	}
}


struct stark_hierarchical_assembler_group_get_contig_seq_and_stats_s {
	list_t(char) seq;
	// list_t(coverage_t) coverage;
	double sum_coverage;
};

void stark_hierarchical_assembler_group_get_contig_cb(struct starknode_navigator_s * navigator, int ch, void * aux) {
	struct stark_hierarchical_assembler_group_get_contig_seq_and_stats_s * seq_and_stats = aux;
	// char_list_t * seq = aux;
	list_insert(&seq_and_stats->seq, nucleotides[ch]);
	// list_insert(&seq_and_stats->coverage, navigator->node->coverage);
	seq_and_stats->sum_coverage += navigator->node->coverage;
}


void stark_hierarchical_assembler_group_get_contig(
	  struct stark_hierarchical_assembler_s * hierarchical_assembler
	, struct stark_hierarchical_assembler_group_s * group
	, struct stark_hierarchical_assembler_group_get_contig_seq_and_stats_s * seq_and_stats) {
		
	// char_list_t seq;
	list_init_size(&seq_and_stats->seq, group->contig_length + 1);
	// list_init_size(&seq_and_stats->coverage, group->contig_length + 1);
	seq_and_stats->sum_coverage = group->endpoint->node->coverage * group->endpoint->depth;
	
	get_sequence(seq_and_stats->seq.list, group->endpoint[0].stark, group->endpoint[0].depth, group->endpoint[0].offset);
	if (group->endpoint[0].direction)
		reverse_complement_onspot(seq_and_stats->seq.list, group->endpoint[0].depth);
	
	seq_and_stats->seq.size = group->endpoint[0].depth;
	
	map_zero(&thread_seen_map);
	ssize_t length = stark_hierarchical_assembler_run_contig(
		  group->endpoint						// start point
		, hierarchical_assembler->minK			// minK
		, group->group_id						// hgroup
		, hierarchical_assembler->groups.size	// invalid hgroup_alt
		, NULL									// output endpoint
		, stark_hierarchical_assembler_group_get_contig_cb
		, seq_and_stats							// callback aux
		);
	
	if (length < 0) {
		CRITICAL_ERROR("Error reading contig from group %llu", group->group_id);
	}
}

// int stark_hierarchical_assembler_init_split_group_dfs_test_cb(
// 		  struct starknode_navigator_s navigator[2]
// 		, void * _stark_hierarchical_assembler_control) {
// 	
// }

void stark_hierarchical_assembler_init_split_group (
	  struct stark_hierarchical_assembler_s * hierarchical_assembler
	, struct stark_hierarchical_assembler_group_s * group
	, ssize_t contig_length);

static inline void stark_hierarchical_assembler_init_group_sanity_check(
	  struct stark_hierarchical_assembler_s * hierarchical_assembler
	, struct stark_hierarchical_assembler_group_s * group) {
	
	group->endpoint[1] = group->endpoint[0];
	
	ssize_t length1, length0;
	
	map_zero(&thread_seen_map);
	// this one saves directional nodes in the map
	if ((length0 = stark_hierarchical_assembler_run_contig(
		  group->endpoint						// start point
		, hierarchical_assembler->minK			// minK
		, group->group_id						// hgroup
		, (hirarchial_assembler_group_t)-1		// invalid hgroup_alt
		, group->endpoint + 1					// output endpoint
		, NULL									// no callback
		, NULL									// callback aux
		)) < 0) {
		// for now we just attempt to split the group
		stark_hierarchical_assembler_init_split_group(hierarchical_assembler, group, length0);
		return;
		CRITICAL_ERROR("Cannot determine endpoint for group %llu.", group->group_id);
	}
	
	group->endpoint[1].direction ^= 1;
	map_zero(&thread_seen_map);
	if ((length1 = stark_hierarchical_assembler_run_contig(
		  group->endpoint + 1					// start point
		, hierarchical_assembler->minK			// minK
		, group->group_id						// hgroup
		, (hirarchial_assembler_group_t)-1		// invalid hgroup_alt
		, group->endpoint					// output endpoint
		, NULL									// no callback
		, NULL									// callback aux
		)) < 0) {
		// for now we just attempt to split the group
		stark_hierarchical_assembler_init_split_group(hierarchical_assembler, group, length1);
		return;
		CRITICAL_ERROR("Cannot determine endpoint for group %llu.", group->group_id);
	}
	// length1 += group->endpoint[1].depth;
	
	group->endpoint[0].direction ^= 1;
	map_zero(&thread_seen_map);
	if((length0 = stark_hierarchical_assembler_run_contig(
		  group->endpoint						// start point
		, hierarchical_assembler->minK			// minK
		, group->group_id						// hgroup
		, (hirarchial_assembler_group_t)-1		// invalid hgroup_alt
		, NULL									// output endpoint
		, NULL									// no callback
		, NULL									// callback aux
		)) < 0) {
		// for now we just attempt to split the group
		stark_hierarchical_assembler_init_split_group(hierarchical_assembler, group, length0);
		return;
		CRITICAL_ERROR("Cannot determine endpoint for group %llu.", group->group_id);
	}
	
	// length0 += group->endpoint[0].depth;
	
	group->contig_length = length0;
	
	// fprintf(stdout, "Groupo %zu contig length is %zd, %zd\n",group->group_id , length0, length1);
	
	if (length0 != length1) {
		// for now we just attempt to split the group
		stark_hierarchical_assembler_init_split_group(hierarchical_assembler, group, 1);
		return;
		DEBUG_MSG("Contig lengths mismatch in group %llu. Lengths: %zd != %zd",group->group_id , length0, length1);
	}
	
	
	#ifdef DEBUG
	stark_hierarchical_assembler_print_group(stdout, group);
	#endif
	
}


int stark_hierarchical_assembler_init_dfs_test_callback(
		  struct starknode_navigator_s navigator[2]
		, void * _stark_hierarchical_assembler_control) {
	struct stark_hierarchical_assembler_grow_group_dfs_cb_s * stark_hierarchical_assembler_control = _stark_hierarchical_assembler_control;
	
	if (!navigator[0].offset  || !navigator[1].offset || !navigator[0].node->depth || !navigator[1].node->depth) {
		CRITICAL_ERROR("Zero navigator offset, something went wrong");
	}
	
	// {
	// 	if (navigator[0].depth == 86 && navigator[0].offset == 60194) {
	// 		debug_print_node(navigator[0].stark, navigator[0].depth, navigator[0].offset);
	// 	}
	// }
	
	#if defined(UBIGRAPHAPI_H) || defined(UBIGRAPITHREADED_H)
	char color[10];
	sprintf(color, "%lld", navigator[0].node->hgroup);
	stark_ubigraph_set_vertex_attribute(navigator[0].node->ubigraph.id, "color", color);
	#endif
	
	struct stark_hierarchical_assembler_group_s * group = stark_hierarchical_assembler_control->group;
	coverage_t coverage[2] = {navigator[0].node->coverage, navigator[1].node->coverage};

	double link_strength = stark_link_strength(coverage);
	link_strength = fabs(link_strength);
	
	if (   navigator[1].node->hgroup != navigator[0].node->hgroup
	    && !(stark_hierarchical_assembler_find_group_by_id(stark_hierarchical_assembler_control->hierarchical_assembler, navigator[1].node->hgroup)->contains_palindrome)
	 	&& link_strength <= stark_hierarchical_assembler_control->current_link_stregth_cutoff) {
		if (navigator[1].node->hgroup != stark_hierarchical_assembler_control->group_orig) {
			stark_hierarchical_assembler_print_group(stderr, group);
			stark_hierarchical_assembler_print_group(stderr, stark_hierarchical_assembler_find_group_by_id(stark_hierarchical_assembler_control->hierarchical_assembler, navigator[1].node->hgroup));
			CRITICAL_ERROR("Thread Collision groups %lld and %lld, this code if not yet written", navigator[0].node->hgroup, navigator[1].node->hgroup);
		}
		
		// if (navigator[1].node == stark_hierarchical_assembler_control->hook_orig.node) {
		// 	stark_hierarchical_assembler_control->hook_orig.node = NULL;
		// }
		stark_hierarchical_assembler_init_add_node_to_group(navigator + 1, group, link_strength);
		return 1;
	} else if (stark_hierarchical_assembler_control->group_orig && navigator[1].node->hgroup == stark_hierarchical_assembler_control->group_orig && !map_insert(&thread_seen_map, (hirarchial_assembler_group_t)navigator[1].node) ) {
		list_insert(&stark_hierarchical_assembler_control->todo_nodes, navigator[1]);
	}
	
	return 0;
	
}

void stark_hierarchical_assembler_init_split_group (
	  struct stark_hierarchical_assembler_s * hierarchical_assembler
	, struct stark_hierarchical_assembler_group_s * group
	, ssize_t contig_length)
{
	
	if (group->internal_link_strength == 0.0) {
		if (contig_length < 0) {
			stark_hierarchical_assembler_group_mark_cyclic(*group);
			return;
		}
		// #ifdef DEBUG
		// stark_hierarchical_assembler_print_group(stderr, group);
		// fprintf(stderr, "For now will continue with uneven endpoint lengths\n");
		// #endif
		return;
		// CRITICAL_ERROR("Splitting group %zu is no longer possible", group->group_id);
	}
	if (group->internal_link_strength < (double)0.000000000232831) {
		group->internal_link_strength = 0.0;
	}
	
	struct starknode_navigator_s hook = *group->endpoint;
	
	// DEBUG_MSG("Splitting group %llu", group->group_id);
	
	int new_groups = 0;
	
	struct stark_hierarchical_assembler_grow_group_dfs_cb_s control_struct = {
		.current_link_stregth_cutoff = (group->internal_link_strength / 2),
		.hierarchical_assembler = hierarchical_assembler,
		// .group = group_alt,
		.group_orig = group->group_id,
	};
	list_init(&(control_struct.todo_nodes));
	
	while (hook.node) {
		new_groups++;
		
		hirarchial_assembler_group_t hgroup_alt = stark_hierarchical_assembler_new_group(hierarchical_assembler);
		struct stark_hierarchical_assembler_group_s * group_alt = hierarchical_assembler->groups.list + hgroup_alt;

		group_alt->group_id = hgroup_alt;

		control_struct.group = group_alt;

		stark_hierarchical_assembler_init_add_node_to_group(&hook, group_alt, 0.0);
		
		thread_dfs_stack.size = 0;
		list_insert(&thread_dfs_stack, hook);
		
		map_zero(&thread_seen_map);
		stark_navigate_dfs(
			  &thread_dfs_stack
			, hierarchical_assembler->minK
			// , &hook
			, stark_hierarchical_assembler_init_dfs_test_callback
			, &control_struct
		);
		
		group->num_nodes -= group_alt->num_nodes;
		group->mass -= group_alt->mass;
		
		hook.node = NULL;
		
		while (control_struct.todo_nodes.size) {
			struct starknode_navigator_s * potential_new_hook = control_struct.todo_nodes.list + --control_struct.todo_nodes.size;
			if (potential_new_hook->node->hgroup == group->group_id) {
				hook = *potential_new_hook;
				break;
			}
		}
		
		stark_hierarchical_assembler_init_group_sanity_check(
			  hierarchical_assembler
			, group_alt
		);
	}
	
	list_free(control_struct.todo_nodes);
	
	// DEBUG_MSG("Split group %llu into %d new groups", group->group_id, new_groups);
	
	stark_hierarchical_assembler_recycle_group(hierarchical_assembler, group);
	
}


int stark_hierarchical_assembler_init_forall_cb(struct starknode_navigator_s * current, void * _hierarchical_assembler) {
	struct stark_hierarchical_assembler_s * hierarchical_assembler = _hierarchical_assembler;
	// stark_lock_t lock;
	if (!current->node->hgroup
	 	// && !STARK_NODE_IS_PALINDROME(current->node[0]) // not a palindrome
		// && !(current->node->flags & STARK_FLAG_NODELETE) // not part of a longer palindrome
														 // TODO: remove edges leading to palindromes
	) {
		hirarchial_assembler_group_t hgroup = stark_hierarchical_assembler_new_group(hierarchical_assembler);
		struct stark_hierarchical_assembler_group_s * group = hierarchical_assembler->groups.list + hgroup;
		
		group->group_id = hgroup;
		// group->endpoint->depth = HARD_MAXDEPTH;
		// group->internal_link_strength = 0.0;
		
		// list_init(&(group->neighbours));
		
		struct stark_hierarchical_assembler_grow_group_dfs_cb_s control_struct = {
			.current_link_stregth_cutoff = 0.05,
			.hierarchical_assembler = hierarchical_assembler,
			.group = group,
			.group_orig = (hirarchial_assembler_group_t)0
		};
		
		map_zero(&thread_seen_map);
		
		map_insert(&thread_seen_map, (size_t)current->node);
		stark_hierarchical_assembler_init_add_node_to_group(current, group, 0.0);
				
		thread_dfs_stack.size = 0;
		list_insert(&thread_dfs_stack, *current);
		
		stark_navigate_dfs( &thread_dfs_stack,
			  hierarchical_assembler->minK
			// , current
			, stark_hierarchical_assembler_init_dfs_test_callback
			, &control_struct
		);
		
	}
	
	return 0;
}

int stark_hierarchical_assembler_init_palindromes_dfs_test_callback(
		  struct starknode_navigator_s navigator[2]
		, void * _stark_hierarchical_assembler_control) {
			
	if (navigator[0].depth > navigator[1].depth) { // new one is a parent
		return stark_hierarchical_assembler_init_dfs_test_callback(navigator, _stark_hierarchical_assembler_control);
	} else
		return 0;
}

int stark_hierarchical_assembler_init_palindromes_forall_cb(struct starknode_navigator_s * current, void * _hierarchical_assembler) {
	struct stark_hierarchical_assembler_s * hierarchical_assembler = _hierarchical_assembler;
	if (!current->node->hgroup
	 	&& STARK_NODE_IS_PALINDROME(current->node[0]) // not a palindrome
		&& current->depth >= hierarchical_assembler->minK
	) {
		hirarchial_assembler_group_t hgroup = stark_hierarchical_assembler_new_group(hierarchical_assembler);
		struct stark_hierarchical_assembler_group_s * group = hierarchical_assembler->groups.list + hgroup;
		
		group->group_id = hgroup;
		group->contains_palindrome = 1;
		// group->endpoint->depth = HARD_MAXDEPTH;
		// group->internal_link_strength = 0.0;
		
		// list_init(&(group->neighbours));
		
		struct stark_hierarchical_assembler_grow_group_dfs_cb_s control_struct = {
			.current_link_stregth_cutoff = 1.0,
			.hierarchical_assembler = hierarchical_assembler,
			.group = group,
			.group_orig = (hirarchial_assembler_group_t)0
		};
		
		map_zero(&thread_seen_map);
		
		map_insert(&thread_seen_map, (size_t)current->node);
		stark_hierarchical_assembler_init_add_node_to_group(current, group, 0.0);
		
		thread_dfs_stack.size = 0;
		list_insert(&thread_dfs_stack, *current);
		
		stark_navigate_dfs( &thread_dfs_stack,
			  hierarchical_assembler->minK
			// , current
			, stark_hierarchical_assembler_init_palindromes_dfs_test_callback
			, &control_struct
		);
		
	}
	
	return 0;
}


int stark_hierarchical_assembler_neighbour_comparator(const void * a, const void * b) {
	double l1 = ((struct stark_hierarchical_assembler_neighbours_global_s *)a)->link_strength;
	double l2 = ((struct stark_hierarchical_assembler_neighbours_global_s *)b)->link_strength;
	
	int ret = 0;
	
	if (l1 < l2)
		ret = -1;
	if (l1 > l2)
		ret = 1;
	
	return ret;
}

int stark_hierarchical_assembler_merge_groups_neighbour_id_comparator(const void * a, const void * b) {
	size_t ida = ((struct stark_hierarchical_assembler_neighbour_s *)a)->id;
	size_t idb = ((struct stark_hierarchical_assembler_neighbour_s *)b)->id;
	
	if (!ida && !idb)
		return 0;
	if (!ida)
		return 1;
	if (!idb)
		return -1;
	
	ssize_t diff = ida - idb;
	
	if (diff > 0)
		return 1;
	if (diff < 0)
		return -1;
	
	return stark_hierarchical_assembler_neighbour_comparator(a, b);
}

typedef list_t(struct stark_hierarchical_assembler_neighbour_s) stark_hierarchical_assembler_neighbours_list_t;

int stark_hierarchical_assembler_init_group_neighbours_dfs_test_cb_thread_safe(
		  struct starknode_navigator_s navigator[2]
		, void * arg)
{
	// struct stark_hierarchical_assembler_s * hierarchical_assembler = _stark_hierarchical_assembler;
	
	
	hirarchial_assembler_group_t group1 = navigator[0].node->hgroup;
	hirarchial_assembler_group_t group2 = navigator[1].node->hgroup;
			
	if (group1 && group2) {
		
		if (group1 != group2) {
			if (group1 < group2) {
				coverage_t coverage[2] = {navigator[0].node->coverage, navigator[1].node->coverage};

				double link_strength = stark_link_strength(coverage);
				link_strength = fabs(link_strength);
				
				struct stark_hierarchical_assembler_neighbour_s new_neighbour = {
					  .link_strength = link_strength
					, .id = group2
				};
				
				stark_hierarchical_assembler_neighbours_list_t * neighbours_list = arg;
				
				list_insert(neighbours_list, new_neighbour);
			}

		} else if (!map_insert(&thread_seen_map, (uint64_t)navigator[0].node)) {
			return 1;
		}
		
		
	}
		
	
	return 0;
}

void stark_hierarchical_assembler_init_group_neighbours_threadsafe(struct stark_hierarchical_assembler_s * hierarchical_assembler) {
	
	#pragma omp single
	{
		hierarchical_assembler->neighbours.size = 0;
	}
	
	#pragma omp barrier
	
	#if defined(STARK_STATUS_H)
	stark_set_status_action("Initialising neighbours");
	#endif
	
	stark_hierarchical_assembler_neighbours_list_t neighbours_list;
	list_init(&neighbours_list);
	
	size_t i;
	const size_t numgroups = hierarchical_assembler->groups.size;
	#pragma omp for schedule(dynamic,0x100)
	for (i = 1; i < numgroups; i++) {
		struct stark_hierarchical_assembler_group_s *group = hierarchical_assembler->groups.list + i;
		
		// it could be a recycled or merged group
		if (group->group_id) {
			
			neighbours_list.size = 0;
			
			
			thread_dfs_stack.size = 0;
			list_insert(&thread_dfs_stack, group->endpoint[0]);
			
			map_zero(&thread_seen_map);
			stark_navigate_dfs( &thread_dfs_stack,
				  hierarchical_assembler->minK
				// , group->endpoint
				, stark_hierarchical_assembler_init_group_neighbours_dfs_test_cb_thread_safe
				, &neighbours_list
			);
			
			list_qsort(neighbours_list, stark_hierarchical_assembler_merge_groups_neighbour_id_comparator);
			
			hirarchial_assembler_group_t last_group_id = 0;
			size_t j;
			#pragma omp critical (merge_neighbours)
			for (j = 0; j < neighbours_list.size; j++) {
				if (last_group_id != neighbours_list.list[j].id) {
					last_group_id = neighbours_list.list[j].id;
					struct stark_hierarchical_assembler_neighbours_global_s neighbours_entry = {
						  .link_strength = neighbours_list.list[j].link_strength
						, .ids = {
							i, neighbours_list.list[j].id
						}
					};

					list_insert(&hierarchical_assembler->neighbours, neighbours_entry);
					
				}
			}

			
		}
	}
	
	list_free(neighbours_list);
	
	#pragma omp single
	{
		list_qsort(hierarchical_assembler->neighbours, stark_hierarchical_assembler_neighbour_comparator);
	}
	
	#pragma omp barrier
	
}

void stark_hierarchical_assembler_init (stark_t * stark, struct stark_hierarchical_assembler_s * hierarchical_assembler, depth_t minK) {
	
	#pragma omp master
	{
		memset(hierarchical_assembler, 0, sizeof(*hierarchical_assembler));
		hierarchical_assembler->stark = stark;
		hierarchical_assembler->minK = minK;
		list_init_flags(&(hierarchical_assembler->groups), 1);
		list_init(&(hierarchical_assembler->recycle_groups));
		list_init(&(hierarchical_assembler->neighbours));

		list_new_empty(&(hierarchical_assembler->groups));
	}
	
	if (!thread_seen_map.map)
		map_init(&thread_seen_map);
	if (!thread_dfs_stack.list)
		list_init(&thread_dfs_stack);
	
	#pragma omp master
	{
		struct starknode_navigator_s start = {
			  .stark = stark
			, .depth = HARD_MAXDEPTH - 1
			, .offset = 1
		};
		
		#if defined(STARK_STATUS_H)
		stark_set_status_action("Initialising HA groups");
		#endif
	
		// palindromes get their own groups
		stark_navigate_loop_forall(&start, -1, stark_hierarchical_assembler_init_palindromes_forall_cb, hierarchical_assembler);
	
		start.depth = minK;
		start.offset = 1;
		stark_navigate_loop_forall(&start, 1, stark_hierarchical_assembler_init_forall_cb, hierarchical_assembler);
	
		// stark_print(stark, HARD_MAXDEPTH);
	}
	
	#pragma omp barrier
	
	#if defined(STARK_STATUS_H)
	stark_set_status_action("Running sanity check on HA groups");
	#endif
	
	size_t i;
	size_t numgroups = *(volatile size_t *)&hierarchical_assembler->groups.size;
	
	#pragma omp barrier
	
	#pragma omp for schedule(dynamic,0x100) nowait
	for (i = 1; i < numgroups; i++) {
		// #if defined(STARK_STATUS_H)
		// stark_set_status_action("Sanity check on HA group %zu out of %zu", i, numgroups);
		// #endif
		
		struct stark_hierarchical_assembler_group_s *group = hierarchical_assembler->groups.list + i;
	
		if (group->group_id) { // is not recycled
			stark_hierarchical_assembler_init_group_sanity_check(hierarchical_assembler, group);
		}
	
		// list_qsort(group->neighbours, stark_hierarchical_assembler_neighbour_comparator);
	}
	
	
	// #pragma omp single
	// {
	// 	#ifdef DEBUG
	// 	stark_print_all_groups(stdout, hierarchical_assembler);
	// 	#endif
	// }
	
	#if defined(STARK_STATUS_H)
	stark_set_status_action("Sleeping");
	#endif
	
	#pragma omp barrier
}

int stark_hierarchical_assembler_merge_groups_dfs_test_callback(
		  struct starknode_navigator_s navigator[2]
		, void * _hgroups) {
	hirarchial_assembler_group_t hgroup_primary = ((hirarchial_assembler_group_t *)_hgroups)[0];
	hirarchial_assembler_group_t hgroup_secondary = ((hirarchial_assembler_group_t *)_hgroups)[1];
	
	// if (!navigator[0].offset  || !navigator[1].offset || !navigator[0].node->depth || !navigator[1].node->depth) {
	// 	CRITICAL_ERROR("Zero navigator offset, something went wrong");
	// }
	
	if (!map_insert(&thread_seen_map, (size_t)navigator[1].node)) {
		
		if (navigator[1].node->hgroup == hgroup_secondary) { // node has no group yet
			
			navigator[1].node->hgroup = hgroup_primary;
			
			#if defined(UBIGRAPHAPI_H) || defined(UBIGRAPITHREADED_H)
			char color[10];
			sprintf(color, "%llu", hgroup_primary);
			stark_ubigraph_set_vertex_attribute(navigator[0].node->ubigraph.id, "color", color);
			#endif
			
			return 1;
			
		}
				
		// if (navigator[1].node->hgroup == hgroup_primary) {
		// 			return 0;
		// 		}
		
		return 0;
	}
	
	return 0;
}


struct stark_hierarchical_assembler_contig_candidate_s {
	struct starknode_navigator_s endpoint[2];
	ssize_t length;
};


void stark_hierarchical_assembler_merge_groups (
	  struct stark_hierarchical_assembler_s * hierarchical_assembler
	, struct stark_hierarchical_assembler_group_s * group[2]
	, struct stark_hierarchical_assembler_contig_candidate_s * contig_candidate
	, double neighbours_link_strength
) {
	int primary = 0;
	if (group[1]->group_id < group[0]->group_id)
		primary = 1;
	
	// set new contig length and endpoints
	memcpy(group[primary]->endpoint, contig_candidate->endpoint, sizeof(contig_candidate->endpoint));
	group[primary]->contig_length = contig_candidate->length;
	
	// transfer all nodes to primary group
	
	hirarchial_assembler_group_t hgroups[2] = {
		  group[primary]->group_id
		, group[primary ^ 1]->group_id
	};
	
	map_zero(&thread_seen_map);
	map_insert(&thread_seen_map, (hirarchial_assembler_group_t)group[primary ^ 1]->endpoint->node);
	group[primary ^ 1]->endpoint->node->hgroup = group[primary]->group_id;
	
	#if defined(UBIGRAPHAPI_H) || defined(UBIGRAPITHREADED_H)
	char color[10];
	sprintf(color, "%llu", group[primary]->group_id);
	stark_ubigraph_set_vertex_attribute(group[primary ^ 1]->endpoint->node->ubigraph.id, "color", color);
	#endif
	
	
	thread_dfs_stack.size = 0;
	list_insert(&thread_dfs_stack, group[primary ^ 1]->endpoint[0]);
	
	stark_navigate_dfs( &thread_dfs_stack,
		  hierarchical_assembler->minK
		// , group[primary ^ 1]->endpoint
		, stark_hierarchical_assembler_merge_groups_dfs_test_callback
		, &hgroups
	);
	
	// add nodes
	group[primary]->num_nodes += group[primary ^ 1]->num_nodes;
	
	// merge mass
	group[primary]->mass += group[primary ^ 1]->mass;
	
	// adjust internal_link_strength
	if (group[primary]->internal_link_strength < group[primary ^ 1]->internal_link_strength) {
		group[primary]->internal_link_strength = group[primary ^ 1]->internal_link_strength;
	}
	if (group[primary]->internal_link_strength < neighbours_link_strength) {
		group[primary]->internal_link_strength = neighbours_link_strength;
	}
	
	// reset secondary group
	// list_free(group[primary ^ 1]->neighbours);
	memset(group[primary ^ 1], 0, sizeof(*(group[primary ^ 1])));
	
	// indicate group redirection
	group[primary ^ 1]->num_nodes = group[primary]->group_id; 
	
}

static inline void stark_hierarchical_assembler_test_groups_for_merging_loop (struct stark_hierarchical_assembler_contig_candidate_s * possible_contig, struct stark_hierarchical_assembler_neighbours_global_s * pairing, int minK){
	map_zero(&thread_seen_map);
	possible_contig->length = stark_hierarchical_assembler_run_contig(
		  possible_contig->endpoint						// start point
		, minK									// minK
		, pairing->ids[0]						// hgroup
		, pairing->ids[1]						// hgroup_alt
		, possible_contig->endpoint + 1		// output endpoint
		, NULL									// callback
		, NULL									// callback aux
		);

	possible_contig->endpoint[1].direction ^= 1;

	
}

int stark_hierarchical_assembler_test_groups_for_merging(
	  struct stark_hierarchical_assembler_s * hierarchical_assembler
	// , hirarchial_assembler_group_t group_ids[2]
	, struct stark_hierarchical_assembler_neighbours_global_s * pairing
) {
	
	struct stark_hierarchical_assembler_group_s * group[2] = {
		  stark_hierarchical_assembler_find_group_by_id(hierarchical_assembler, pairing->ids[0])
		, stark_hierarchical_assembler_find_group_by_id(hierarchical_assembler, pairing->ids[1])
	};
	struct stark_hierarchical_assembler_contig_candidate_s possible_contigs[4];
	
	possible_contigs[0].endpoint[0] = group[0]->endpoint[0];
	possible_contigs[1].endpoint[0] = group[0]->endpoint[1];
	possible_contigs[2].endpoint[0] = group[1]->endpoint[0];
	possible_contigs[3].endpoint[0] = group[1]->endpoint[1];
	
	int i;
	int better_contig = -1;
	
	int largegroup = 0;
	if (group[0]->num_nodes + group[1]->num_nodes > 0x1000)
		largegroup = 1;
	
	for (i = 0; i < 4; i++) {
		// int taskit = largegroup && (hierarchical_assembler->threads_idle > 0);
		// #pragma omp task if(taskit != 0) default(none) firstprivate(i) shared(possible_contigs, pairing, hierarchical_assembler)
		{
			stark_hierarchical_assembler_test_groups_for_merging_loop(possible_contigs + i, pairing, hierarchical_assembler->minK);
		}
	}
	
	#pragma omp taskwait
	
	for (i = 0; i < 4; i++) {
		// stark_hierarchical_assembler_test_groups_for_merging_loop(possible_contigs + i, pairing);
		if (
			   possible_contigs[i].length > group[0]->contig_length
			&& possible_contigs[i].length > group[1]->contig_length
			&& (better_contig < 0 || possible_contigs[i].length > possible_contigs[better_contig].length)
		) {
			better_contig = i;
		}
	}
	
	
	if (better_contig >= 0) {
		stark_hierarchical_assembler_merge_groups (
			  hierarchical_assembler
			, group
			, possible_contigs + better_contig
			, pairing->link_strength
		);
		
		return 1;
	} else
		return 0;
	
}

static inline int stark_hierarchical_assembler_resolve_pairing_used_up_groups(struct stark_hierarchical_assembler_s * hierarchical_assembler, struct stark_hierarchical_assembler_neighbours_global_s * pairing)
{
	
	while (pairing->ids[0] && !hierarchical_assembler->groups.list[pairing->ids[0]].group_id) { // this group no longer exists
		pairing->ids[0] = hierarchical_assembler->groups.list[pairing->ids[0]].num_nodes;
	}
	while (pairing->ids[1] && !hierarchical_assembler->groups.list[pairing->ids[1]].group_id) { // this group no longer exists
		pairing->ids[1] = hierarchical_assembler->groups.list[pairing->ids[1]].num_nodes;
	}

	if (!pairing->ids[0] || !pairing->ids[1]) {
		fprintf(stderr, "Got an empty group id again...\n");
	}

	if (pairing->ids[0] == pairing->ids[1]) {
		return 1;
	}
	
	if (!pairing->ids[0] || !pairing->ids[1]) {
		// FILE * missing_pairings_fp = fopen("missing_pairings.txt", "a");
		// uint64_t id = hierarchical_assembler->neighbours.list[i].ids[0];
		// fprintf(missing_pairings_fp, "pairing: %.4lf\n%llu", pairing->link_strength, id);
		// while (id && !hierarchical_assembler->groups.list[id].group_id) {
		// 	id = hierarchical_assembler->groups.list[id].num_nodes;
		// 	fprintf(missing_pairings_fp, " -> %llu", id);
		// }
		// 
		// id = hierarchical_assembler->neighbours.list[i].ids[1];
		// fprintf(missing_pairings_fp, "\n%llu", id);
		// while (id && !hierarchical_assembler->groups.list[id].group_id) {
		// 	id = hierarchical_assembler->groups.list[id].num_nodes;
		// 	fprintf(missing_pairings_fp, " -> %llu", id);
		// }
		// fprintf(missing_pairings_fp, "\n");
		// 
		// fclose(missing_pairings_fp);

		// CRITICAL_ERROR("Missing groups, something went wrong");
		fprintf(stderr, "Missing groups, something went wrong\n");
		// hierarchical_assembler->neighbours.list[i].ids[0] = ((hirarchial_assembler_group_t)-1);
		return 1;
	}

	if (hierarchical_assembler->groups.list[pairing->ids[1]].contains_palindrome || hierarchical_assembler->groups.list[pairing->ids[0]].contains_palindrome) {
		// hierarchical_assembler->neighbours.list[i].ids[0] = ((hirarchial_assembler_group_t)-1);
		return 1; // skip palindromic groups
	}

	
	
	
	return 0;

}

static void map_zero_parallel(hashmap_uint64_t * const map) {
	char *array = (void *)map->map;
	size_t i, size = map->incore_length;
	
	#pragma omp for schedule(dynamic, 0x10) nowait
	for (i = 0; i < size; i++) {
		if (map->incore[i]) {
			memset(array + (i * sysconf(_SC_PAGESIZE)), 0, sysconf(_SC_PAGESIZE));
			map->incore[i] = 0;
		}
	}
	
	
	#pragma omp master
	{
		// memset(map->incore, 0, map->incore_length);
		map->size = 0;
	}
}


size_t stark_hierarchical_assembler_test_and_merge_groups_openmp(struct stark_hierarchical_assembler_s * hierarchical_assembler, size_t max_rounds) {
	int master = -1;
	int team_size;
	
	struct stark_hierarchical_assembler_test_and_merge_groups_dispatch_s * dispatch;
	
	#if defined(STARK_STATUS_H)
	stark_set_status_action("HA Slave");
	#endif
	
	#pragma omp master
	{
		
		#if defined(STARK_STATUS_H)
		stark_set_status_action("HA master");
		#endif
		
		master = omp_get_thread_num();
		team_size = omp_get_num_threads();
		dispatch = calloc(sizeof(*dispatch), team_size);
		hierarchical_assembler->dispatch = dispatch;
		map_init_size(&hierarchical_assembler->locked_groups_map, hierarchical_assembler->groups.size);
		
	}
	
	#pragma omp barrier
	
	struct stark_hierarchical_assembler_test_and_merge_groups_dispatch_s * my_dispatch = hierarchical_assembler->dispatch + omp_get_thread_num();

	circular_ring_buffer_init(&my_dispatch->buffer, sizeof(*hierarchical_assembler->neighbours.list) * 0x100);
	map_init(&my_dispatch->dispatch_map);
	
	#pragma omp barrier
	
	size_t rounds = 0;
	// size_t merged;
	
	do {
		stark_hierarchical_assembler_init_group_neighbours_threadsafe(hierarchical_assembler);
		
		#if defined(STARK_STATUS_H)
		stark_set_status_action("HA Slave");
		stark_thread_status->dispatch = my_dispatch;
		#endif
		
		rounds++;
		// merged = 0;
		
		struct rusage start, end;
		
		#pragma omp master
		{
			hierarchical_assembler->merged = 0;
			getrusage(RUSAGE_SELF, &start);
			#if defined(STARK_STATUS_H)
			stark_set_status_action("HA dispatch master round %zu", rounds);
			#endif
		}
		
		// size_t unhandled_pairings;
		// size_t nextstart = 0;
		// size_t pairings_in_round = hierarchical_assembler->neighbours.size;
		
		
		
		size_t pairings_done = 0;
		size_t skipfrom = 0;
		
		do {
			
			if (my_dispatch->dispatch_map.size)
				map_zero(&my_dispatch->dispatch_map);
			
			my_dispatch->buffer.eof = 0;
			
			// if (hierarchical_assembler->locked_groups_map.size)
			map_zero_parallel(&hierarchical_assembler->locked_groups_map);
			
			#pragma omp atomic
			hierarchical_assembler->threads_idle++;
			
			#pragma omp barrier
			
			#pragma omp atomic
			hierarchical_assembler->threads_idle--;
			
			
			#pragma omp master
			{			
				size_t i;
				size_t write_next_pairing = 0;
				
				for (i = 0;; i++) {
					if (i == skipfrom)
						i += pairings_done;
					
					if (i >= hierarchical_assembler->neighbours.size)
						break;
				
					struct stark_hierarchical_assembler_neighbours_global_s pairing = hierarchical_assembler->neighbours.list[i];
				
					if (stark_hierarchical_assembler_resolve_pairing_used_up_groups(hierarchical_assembler, &pairing)) {
						pairings_done++;
						continue;
					}
					
					uint64_t in_map[2] = {
						  map_get(&hierarchical_assembler->locked_groups_map, pairing.ids[0])
						, map_get(&hierarchical_assembler->locked_groups_map, pairing.ids[1])
					};
				
					int used_id = -1;
					int best_dispatch = -1;
					
					int idle_threads = 0;
					
					int j;
					for (j = 1; j < team_size; j++) {
						struct stark_hierarchical_assembler_test_and_merge_groups_dispatch_s * his_dispatch = dispatch + j;
						if (
							   (!in_map[0] || map_get(&his_dispatch->dispatch_map, pairing.ids[0]))
							&& (!in_map[1] || map_get(&his_dispatch->dispatch_map, pairing.ids[1]))
							&& (
								best_dispatch < 0
							 	|| circular_ring_buffer_get_free_space(&dispatch[best_dispatch].buffer) < circular_ring_buffer_get_free_space(&his_dispatch->buffer)
								)
						) { // it fits here and is a better fit then the one before
							best_dispatch = j;
						}
						
						if (write_next_pairing > 0x1000 && his_dispatch->dispatch_map.size > 0x100 && !circular_ring_buffer_get_current_size(&his_dispatch->buffer)) {
							idle_threads++;
						}
					}
					
					map_insert(&hierarchical_assembler->locked_groups_map, pairing.ids[0]);
					map_insert(&hierarchical_assembler->locked_groups_map, pairing.ids[1]);
									
					if (best_dispatch >= 0) {
						if (!circular_ring_buffer_write(&dispatch[best_dispatch].buffer, &pairing, sizeof(pairing), CIRCULAR_RING_BUFFER_WRITE_WHOLE)) {
							// TODO
							// buffer is full
							// rewind
							break;
						}
						pairings_done++;
						// *map_location[0] = pairing.ids[0];
						// *map_location[1] = pairing.ids[1];
						
						map_insert(&dispatch[best_dispatch].dispatch_map, pairing.ids[0]);
						map_insert(&dispatch[best_dispatch].dispatch_map, pairing.ids[1]);
					} else {
						// didnt find a free dispatch
						hierarchical_assembler->neighbours.list[write_next_pairing++] = pairing;
					}
				
					
					if ((idle_threads << 1) > team_size)
						break;
					
				}
				
				skipfrom = write_next_pairing;
				
				#if defined(STARK_STATUS_H)
				stark_set_status_action("HA master round %zu, %.0f%%", rounds, ((100.0 * pairings_done) / hierarchical_assembler->neighbours.size));
				#endif
				
				if (pairings_done == hierarchical_assembler->neighbours.size)
					hierarchical_assembler->neighbours.size = 0;
				
				map_shrink(&hierarchical_assembler->locked_groups_map, hierarchical_assembler->locked_groups_map.size << 1);
				
				int j;
				for (j = 1; j < team_size; j++) {
					circular_ring_buffer_write_eof(&dispatch[j].buffer);
				}
				
			}
			
			if (master != omp_get_thread_num()) {
				size_t thread_merged = 0;
				while (circular_ring_buffer_active(&my_dispatch->buffer)) {
					struct stark_hierarchical_assembler_neighbours_global_s * pairing;
					if (circular_ring_buffer_get_current_size(&my_dispatch->buffer) >= sizeof(*pairing)) {
						pairing = (void *)circular_ring_buffer_get_read_ptr(&my_dispatch->buffer);
						if (!stark_hierarchical_assembler_resolve_pairing_used_up_groups(hierarchical_assembler, pairing)
							&& stark_hierarchical_assembler_test_groups_for_merging(hierarchical_assembler, pairing)
						)
							thread_merged++;
						
						circular_ring_buffer_flush_read(&my_dispatch->buffer, sizeof(*pairing));
					}
					// else {
					// 	#pragma omp taskwait
					// }
				}
				
				#pragma omp atomic
				hierarchical_assembler->merged += thread_merged;
				
			}
			
			
		} while (hierarchical_assembler->neighbours.size);
		
		#pragma omp barrier
		
		#pragma omp master
		{
			fprintf(stderr, "Merged %zu groups in Round %zu\n", hierarchical_assembler->merged, rounds);
			getrusage(RUSAGE_SELF, &end);
			char timebuffer[128];
			print_timediff(timebuffer, &end, &start);
			fprintf(stderr, "Assembly round %zu took %s\n", rounds, timebuffer);
		}
		
	} while (hierarchical_assembler->merged && rounds < max_rounds);
	
	#if defined(STARK_STATUS_H)
	stark_thread_status->dispatch = NULL;
	#endif

	circular_ring_buffer_free(&my_dispatch->buffer);
	
	#pragma omp barrier
	
	map_free(&my_dispatch->dispatch_map);
	
	
	#pragma omp master
	{
		map_free(&hierarchical_assembler->locked_groups_map);
		hierarchical_assembler->dispatch = NULL;
		free(dispatch);
	}
	
	return rounds;
	
}



size_t stark_hierarchical_assembler_test_and_merge_groups_openmp_v2(struct stark_hierarchical_assembler_s * hierarchical_assembler, size_t max_rounds) {
	int master = -1;
	int team_size;
	
	struct stark_hierarchical_assembler_test_and_merge_groups_dispatch_s * dispatch;
	
	#if defined(STARK_STATUS_H)
	stark_set_status_action("HA Slave");
	#endif
	
	#pragma omp master
	{
		
		#if defined(STARK_STATUS_H)
		stark_set_status_action("HA master");
		#endif
		
		master = omp_get_thread_num();
		team_size = omp_get_num_threads();
		dispatch = calloc(sizeof(*dispatch), team_size);
		hierarchical_assembler->dispatch = dispatch;
		map_init_size(&hierarchical_assembler->locked_groups_map, hierarchical_assembler->groups.size);
		
	}
	
	#pragma omp barrier
	
	struct stark_hierarchical_assembler_test_and_merge_groups_dispatch_s * my_dispatch = hierarchical_assembler->dispatch + omp_get_thread_num();

	circular_ring_buffer_init(&my_dispatch->buffer, sizeof(*hierarchical_assembler->neighbours.list) * 0x40);
	map_init(&my_dispatch->dispatch_map);
	
	#pragma omp barrier
	
	size_t rounds = 0;
	// size_t merged;
	
	do {
		stark_hierarchical_assembler_init_group_neighbours_threadsafe(hierarchical_assembler);
		
		#if defined(STARK_STATUS_H)
		stark_set_status_action("HA Slave");
		stark_thread_status->dispatch = my_dispatch;
		#endif
		
		rounds++;
		// merged = 0;
		
		struct rusage start, end;
		
		#pragma omp master
		{
			hierarchical_assembler->merged = 0;
			getrusage(RUSAGE_SELF, &start);
			#if defined(STARK_STATUS_H)
			stark_set_status_action("HA dispatch master round %zu", rounds);
			#endif
		}
		
		// size_t unhandled_pairings;
		// size_t nextstart = 0;
		// size_t pairings_in_round = hierarchical_assembler->neighbours.size;
		
		
		
		size_t pairings_done = 0;
		size_t skipfrom = 0;
		
		#pragma omp master
		do {
			
			// if (my_dispatch->dispatch_map.size)
				
			
			// my_dispatch->buffer.eof = 0;
			
			// if (hierarchical_assembler->locked_groups_map.size)
			
			
			// #pragma omp atomic
			// hierarchical_assembler->threads_idle++;
			
			// #pragma omp barrier
			
			// #pragma omp atomic
			// hierarchical_assembler->threads_idle--;
			
			
			// #pragma omp master
			{			
				size_t i;
				size_t write_next_pairing = 0;
				
				for (i = 0;; i++) {
					if (i == skipfrom)
						i += pairings_done;
					
					if (i >= hierarchical_assembler->neighbours.size)
						break;
				
					struct stark_hierarchical_assembler_neighbours_global_s pairing = hierarchical_assembler->neighbours.list[i];
				
					if (stark_hierarchical_assembler_resolve_pairing_used_up_groups(hierarchical_assembler, &pairing)) {
						pairings_done++;
						continue;
					}
					
					uint64_t in_map[2] = {
						  map_get(&hierarchical_assembler->locked_groups_map, pairing.ids[0])
						, map_get(&hierarchical_assembler->locked_groups_map, pairing.ids[1])
					};
				
					int used_id = -1;
					int best_dispatch = -1;
					
					int idle_threads = 0;
					
					int j;
					for (j = 1; j < team_size; j++) {
						struct stark_hierarchical_assembler_test_and_merge_groups_dispatch_s * his_dispatch = dispatch + j;
						if (
							   (!in_map[0] || map_get(&his_dispatch->dispatch_map, pairing.ids[0]))
							&& (!in_map[1] || map_get(&his_dispatch->dispatch_map, pairing.ids[1]))
							&& (
								best_dispatch < 0
							 	|| circular_ring_buffer_get_free_space(&dispatch[best_dispatch].buffer) < circular_ring_buffer_get_free_space(&his_dispatch->buffer)
								)
						) { // it fits here and is a better fit then the one before
							best_dispatch = j;
						}
						
						if (write_next_pairing > 0x1000 && his_dispatch->dispatch_map.size > 0x100 && !circular_ring_buffer_get_current_size(&his_dispatch->buffer)) {
							idle_threads++;
						}
					}
					
					map_insert(&hierarchical_assembler->locked_groups_map, pairing.ids[0]);
					map_insert(&hierarchical_assembler->locked_groups_map, pairing.ids[1]);
									
					if (best_dispatch >= 0) {
						if (!circular_ring_buffer_write(&dispatch[best_dispatch].buffer, &pairing, sizeof(pairing), CIRCULAR_RING_BUFFER_WRITE_WHOLE)) {
							// TODO
							// buffer is full
							// rewind
							break;
						}
						pairings_done++;
						// *map_location[0] = pairing.ids[0];
						// *map_location[1] = pairing.ids[1];
						
						map_insert(&dispatch[best_dispatch].dispatch_map, pairing.ids[0]);
						map_insert(&dispatch[best_dispatch].dispatch_map, pairing.ids[1]);
					} else {
						// didnt find a free dispatch
						hierarchical_assembler->neighbours.list[write_next_pairing++] = pairing;
					}
				
					
					if ((idle_threads << 1) > team_size)
						break;
					
				}
				
				skipfrom = write_next_pairing;
				
				#if defined(STARK_STATUS_H)
				stark_set_status_action("HA master round %zu, %.0f%%", rounds, ((100.0 * pairings_done) / hierarchical_assembler->neighbours.size));
				#endif
				
				if (pairings_done == hierarchical_assembler->neighbours.size)
					hierarchical_assembler->neighbours.size = 0;
				else {
					map_shrink(&hierarchical_assembler->locked_groups_map, hierarchical_assembler->locked_groups_map.size << 1);
					// map_zero(&my_dispatch->dispatch_map);
					map_zero(&hierarchical_assembler->locked_groups_map);
					
					int j;
					for (j = 1; j < team_size; j++) {
						map_zero(&dispatch[j].dispatch_map);
						
						size_t read_offset = dispatch[j].buffer.read_offset;
						struct stark_hierarchical_assembler_neighbours_global_s * pairing = (char *)dispatch[j].buffer.address + (read_offset & (dispatch[j].buffer.length - 1));
						size_t num_pairings = dispatch[j].buffer.write_offset - read_offset;
						num_pairings &= dispatch[j].buffer.length - 1;
						num_pairings /= sizeof(*pairing);
						
						int k;
						for (k = 0; k < num_pairings; k++) {
							map_insert(&hierarchical_assembler->locked_groups_map, pairing[k].ids[0]);
							map_insert(&hierarchical_assembler->locked_groups_map, pairing[k].ids[1]);
							map_insert(&dispatch[j].dispatch_map, pairing[k].ids[0]);
							map_insert(&dispatch[j].dispatch_map, pairing[k].ids[1]);
						}
					}
					
				}
				
			}
			
			
			
			
		} while (hierarchical_assembler->neighbours.size);
		
		#pragma omp master
		{
			int j;
			for (j = 1; j < team_size; j++) {
				circular_ring_buffer_write_eof(&dispatch[j].buffer);
			}
		}
		
		
		if (master != omp_get_thread_num()) {
			size_t thread_merged = 0;
			while (circular_ring_buffer_active(&my_dispatch->buffer)) {
				struct stark_hierarchical_assembler_neighbours_global_s * pairing;
				if (circular_ring_buffer_get_current_size(&my_dispatch->buffer) >= sizeof(*pairing)) {
					pairing = (void *)circular_ring_buffer_get_read_ptr(&my_dispatch->buffer);
					if (!stark_hierarchical_assembler_resolve_pairing_used_up_groups(hierarchical_assembler, pairing)
						&& stark_hierarchical_assembler_test_groups_for_merging(hierarchical_assembler, pairing)
					)
						thread_merged++;
					
					circular_ring_buffer_flush_read(&my_dispatch->buffer, sizeof(*pairing));
				}
				// else {
				// 	#pragma omp taskwait
				// }
			}
			
			#pragma omp atomic
			hierarchical_assembler->merged += thread_merged;
			
		}
		
		my_dispatch->buffer.eof = 0;
		
		#pragma omp barrier
		
		#pragma omp master
		{
			fprintf(stderr, "Merged %zu groups in Round %zu\n", hierarchical_assembler->merged, rounds);
			getrusage(RUSAGE_SELF, &end);
			char timebuffer[128];
			print_timediff(timebuffer, &end, &start);
			fprintf(stderr, "Assembly round %zu took %s\n", rounds, timebuffer);
		}
		
	} while (hierarchical_assembler->merged && rounds < max_rounds);
	
	#if defined(STARK_STATUS_H)
	stark_thread_status->dispatch = NULL;
	#endif

	circular_ring_buffer_free(&my_dispatch->buffer);
	
	#pragma omp barrier
	
	map_free(&my_dispatch->dispatch_map);
	
	
	#pragma omp master
	{
		map_free(&hierarchical_assembler->locked_groups_map);
		hierarchical_assembler->dispatch = NULL;
		free(dispatch);
	}
	
	return rounds;
	
}



size_t stark_hierarchical_assembler_test_and_merge_groups_openmp_slave (struct stark_hierarchical_assembler_s * hierarchical_assembler, struct stark_hierarchical_assembler_neighbours_global_s * pairings, size_t num_pairings) {
	size_t merged = 0;
	size_t i;
	for (i = 0; i < num_pairings; i++) {
		
		if (stark_hierarchical_assembler_resolve_pairing_used_up_groups(hierarchical_assembler, pairings + i))
			continue;
		
		if (stark_hierarchical_assembler_test_groups_for_merging(hierarchical_assembler, pairings + i)) {
			merged++;
		}
	}
	
	return merged;
}

size_t stark_hierarchical_assembler_test_and_merge_groups_openmp_master(struct stark_hierarchical_assembler_s * hierarchical_assembler, size_t max_rounds) {
	
	hashmap_uint64_t locked_groups_map = {0};
	hashmap_uint64_t this_dispatch_map = {0};
	
	#pragma omp master
	{
		
		#if defined(STARK_STATUS_H)
		stark_set_status_action("HA master");
		#endif
	}
	
	size_t max_pairings = 0x1000;
	size_t max_nodes_dispatch = 0x100000;
	size_t rounds = 0;
	size_t numdispatches = SIZE_MAX;
	
	do {
		stark_hierarchical_assembler_init_group_neighbours_threadsafe(hierarchical_assembler);
		rounds++;
		
		struct rusage start, end;
		
		#pragma omp master
		{
			if (!locked_groups_map.map && hierarchical_assembler->neighbours.size)
			{
				map_init_size(&locked_groups_map, hierarchical_assembler->neighbours.size);
			}
			
			hierarchical_assembler->merged = 0;
			getrusage(RUSAGE_SELF, &start);
			#if defined(STARK_STATUS_H)
			stark_set_status_action("HA disparch master round %zu", rounds);
			#endif
		}
		
		// size_t unhandled_pairings;
		// size_t nextstart = 0;
		size_t pairings_in_round = hierarchical_assembler->neighbours.size;
		
		
		size_t write_next_pairing;
		#pragma omp master
		do {
			
			struct {
				struct stark_hierarchical_assembler_neighbours_global_s * pairings;
				size_t size;
				size_t num_nodes;
			} dispatch_pairings;
			
			// max_pairings = hierarchical_assembler->neighbours.size / 0x1000;
			// 
			// if (max_pairings > 0x1000)
			// 	max_pairings = 0x1000;
			// if (!max_pairings)
			// 	max_pairings = 1;
			
			if (numdispatches < 0x400)
				max_pairings >>= 1;
			
			if (max_pairings < 0x10)
				max_pairings = 0x10;
				
			numdispatches = 0;
			
			
			map_init_size(&this_dispatch_map, max_pairings << 4);
			
			// unhandled_pairings = 0;
			map_zero(&locked_groups_map);
		
			memset(&dispatch_pairings, 0, sizeof(dispatch_pairings));
			dispatch_pairings.pairings = calloc(sizeof(*dispatch_pairings.pairings), max_pairings);
		
			size_t i;
			write_next_pairing = 0;
			for (i = 0; i < hierarchical_assembler->neighbours.size; i++) {
				struct stark_hierarchical_assembler_neighbours_global_s pairing = hierarchical_assembler->neighbours.list[i];
				
				if (stark_hierarchical_assembler_resolve_pairing_used_up_groups(hierarchical_assembler, &pairing))
					continue;

				if (
					   (!map_insert(&locked_groups_map, pairing.ids[0]) || map_get(&this_dispatch_map, pairing.ids[0]))
				 	&& (!map_insert(&locked_groups_map, pairing.ids[1]) || map_get(&this_dispatch_map, pairing.ids[1]))
				) {
					// disable this pairing
					// hierarchical_assembler->neighbours.list[i].ids[0] = ((hirarchial_assembler_group_t)-1);
					
					map_insert(&this_dispatch_map, pairing.ids[0]);
					map_insert(&this_dispatch_map, pairing.ids[1]);
				
				
					dispatch_pairings.pairings[dispatch_pairings.size++] = pairing;
					dispatch_pairings.num_nodes += hierarchical_assembler->groups.list[pairing.ids[0]].num_nodes;
					dispatch_pairings.num_nodes += hierarchical_assembler->groups.list[pairing.ids[1]].num_nodes;
				
					if (dispatch_pairings.num_nodes >= max_nodes_dispatch || dispatch_pairings.size == max_pairings) {
						#pragma omp task firstprivate(dispatch_pairings)
						{
							#if defined(STARK_STATUS_H)
							stark_set_status_action("HA slave");
							#endif
							size_t task_merged = stark_hierarchical_assembler_test_and_merge_groups_openmp_slave (hierarchical_assembler, dispatch_pairings.pairings, dispatch_pairings.size);
							free(dispatch_pairings.pairings);
							
							#pragma omp atomic
							hierarchical_assembler->merged += task_merged;
							
							#if defined(STARK_STATUS_H)
							stark_set_status_action("HA Sleeping");
							#endif
						}
						#if defined(STARK_STATUS_H)
						stark_set_status_action("HA master round %zu, %.0f%%", rounds, 100.0 - (100.0 * (hierarchical_assembler->neighbours.size) / pairings_in_round));
						#endif
						
						numdispatches++;
						memset(&dispatch_pairings, 0, sizeof(dispatch_pairings));
						map_zero(&this_dispatch_map);
						dispatch_pairings.pairings = calloc(sizeof(*dispatch_pairings.pairings), max_pairings);
					}
					
				
				} else {
					hierarchical_assembler->neighbours.list[write_next_pairing++] = pairing;
				}
			
			}
			
			// truncate pairing list
			hierarchical_assembler->neighbours.size = write_next_pairing;
			
			
			
			if (dispatch_pairings.num_nodes) {
				size_t task_merged = stark_hierarchical_assembler_test_and_merge_groups_openmp_slave (hierarchical_assembler, dispatch_pairings.pairings, dispatch_pairings.size);
				free(dispatch_pairings.pairings);
				
				#pragma omp atomic
				hierarchical_assembler->merged += task_merged;
			}
			
			map_free(&this_dispatch_map);
			
			#if defined(STARK_STATUS_H)
			stark_set_status_action("HA master round %zu, %.0f%%", rounds, 100.0 - (100.0 * (hierarchical_assembler->neighbours.size) / pairings_in_round));
			#endif
			
			#pragma omp taskwait
			
		} while (write_next_pairing);
		
		#pragma omp taskwait
		#pragma omp barrier
		
		#pragma omp master
		{
			fprintf(stderr, "Merged %zu groups in Round %zu\n", hierarchical_assembler->merged, rounds);
			getrusage(RUSAGE_SELF, &end);
			char timebuffer[128];
			print_timediff(timebuffer, &end, &start);
			fprintf(stderr, "Assembly round %zu took %s\n", rounds, timebuffer);
		}
		
	} while (hierarchical_assembler->merged && rounds < max_rounds);
	
	
	
	// #ifdef DEBUG
	// #pragma omp master
	// {
	// 	stark_print_all_groups(stderr, hierarchical_assembler);
	// }
	// #endif
	
	return rounds;
}


size_t stark_hierarchical_assembler_test_and_merge_groups(struct stark_hierarchical_assembler_s * hierarchical_assembler, size_t max_rounds) {
	
	
	size_t i;
	size_t merged = 0;
	size_t rounds = 0;
	
	do {
		
		stark_hierarchical_assembler_init_group_neighbours_threadsafe(hierarchical_assembler);
		rounds++;
		merged = 0;
		for (i = 0; i < hierarchical_assembler->neighbours.size; i++) {
			struct stark_hierarchical_assembler_neighbours_global_s pairing = hierarchical_assembler->neighbours.list[i];
		
			// DEBUG_MSG("Found worthy linkkage groups %llu, %llu link %.4lf", pairing->ids[0], pairing->ids[1], pairing->link_strength);
		
			while (pairing.ids[0] && !hierarchical_assembler->groups.list[pairing.ids[0]].group_id) { // this group no longer exists
				pairing.ids[0] = hierarchical_assembler->groups.list[pairing.ids[0]].num_nodes;
			}
			while (pairing.ids[1] && !hierarchical_assembler->groups.list[pairing.ids[1]].group_id) { // this group no longer exists
				pairing.ids[1] = hierarchical_assembler->groups.list[pairing.ids[1]].num_nodes;
			}
			
			if (!pairing.ids[0] || !pairing.ids[1]) {
				fprintf(stderr, "Got an empty group id again...\n");
			}
			
			if (pairing.ids[0] == pairing.ids[1])
				continue;
		
			if (!pairing.ids[0] || !pairing.ids[1]) {
				FILE * missing_pairings_fp = fopen("missing_pairings.txt", "a");
				uint64_t id = hierarchical_assembler->neighbours.list[i].ids[0];
				fprintf(missing_pairings_fp, "pairing: %.4lf\n%llu", pairing.link_strength, id);
				while (id && !hierarchical_assembler->groups.list[id].group_id) {
					id = hierarchical_assembler->groups.list[id].num_nodes;
					fprintf(missing_pairings_fp, " -> %llu", id);
				}
				
				id = hierarchical_assembler->neighbours.list[i].ids[1];
				fprintf(missing_pairings_fp, "\n%llu", id);
				while (id && !hierarchical_assembler->groups.list[id].group_id) {
					id = hierarchical_assembler->groups.list[id].num_nodes;
					fprintf(missing_pairings_fp, " -> %llu", id);
				}
				fprintf(missing_pairings_fp, "\n");
				
				fclose(missing_pairings_fp);
				
				// CRITICAL_ERROR("Missing groups, something went wrong");
				fprintf(stderr, "Missing groups, something went wrong\n");
				continue;
			}
		
			if (hierarchical_assembler->groups.list[pairing.ids[1]].contains_palindrome || hierarchical_assembler->groups.list[pairing.ids[0]].contains_palindrome)
				continue; // skip palindromic groups
		
			// test and merge
			if (!stark_hierarchical_assembler_test_groups_for_merging(hierarchical_assembler, &pairing)) {
				// DEBUG_MSG("Did not merge groups %llu and %llu at link strength %.4lf", pairing.ids[0], pairing.ids[1], pairing.link_strength);
			} else {
				// DEBUG_MSG("Merged groups %llu and %llu at link strength %.4lf", pairing.ids[0], pairing.ids[1], pairing.link_strength);
				merged++;
			}
		
		}
	} while (merged && rounds < max_rounds);
	
	
	
	// #ifdef DEBUG
	// stark_print_all_groups(stderr, hierarchical_assembler);
	// #endif
	
	return rounds;
}


size_t stark_hierarchical_assembler_print_all_group_contigs(FILE * fp, struct stark_hierarchical_assembler_s * hierarchical_assembler, size_t bounds[2]) {
	// char_list_t seq;
	struct stark_hierarchical_assembler_group_get_contig_seq_and_stats_s seq_and_stats;
	
	size_t i, j, bytes_printed = 0;
	for (i = 1; i < hierarchical_assembler->groups.size; i++) {
		struct stark_hierarchical_assembler_group_s *group = hierarchical_assembler->groups.list + i;
		
		if (group->group_id && group->contig_length > bounds[0] && group->contig_length <= bounds[1]) {
			stark_hierarchical_assembler_group_get_contig(hierarchical_assembler, group, &seq_and_stats);

			// list_insert(&seq, '\0');
			
			bytes_printed += fprintf(fp, ">CONTIG_%zu%s%s length:%zu average coverage:%.2lf\n"
				, i
				, group->contains_palindrome ? " palindromic" : ""
				, group->flags & STARK_HIERARCHIAL_ASSEMBLER_GROUP_FLAG_CYCLIC ? " cyclic" : ""
				, seq_and_stats.seq.size
				, seq_and_stats.sum_coverage / seq_and_stats.seq.size);
			
			size_t bytes_to_print = seq_and_stats.seq.size;
			for (j = 0; j < seq_and_stats.seq.size && bytes_to_print > 80; j+= 80, bytes_to_print -= 80) 
				bytes_printed += fprintf(fp, "%.*s\n", 80, seq_and_stats.seq.list + j);
			if (bytes_to_print)
				bytes_printed += fprintf(fp, "%.*s\n", (int)bytes_to_print, seq_and_stats.seq.list + j);
			// puts(seq.list);

			list_free(seq_and_stats.seq);
			// list_free(seq_and_stats.coverage);
		}
	}
	
	return bytes_printed;	
}



#endif
