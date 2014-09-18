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

#ifndef STARK_NAVIGATE_H_TNO9P00H
#define STARK_NAVIGATE_H_TNO9P00H


struct starknode_navigator_s {
	stark_t* stark;
	starknode_t *node;
	offset_t offset;
	depth_t depth;
	int_fast8_t direction;
	int_fast8_t transition_char;
	// starknode_t node;
	// size_t score;
};

static inline void stark_update_starknode_navigator_down(struct starknode_navigator_s* const __restrict current, offset_t child_offset) {
	current->depth++;
		
	current->node = current->stark->level[current->depth] + child_offset;

	if (current->node->flags & ((current->node->uplink[0] == current->offset) ? STARK_FLAG_PARENT1_REVERSE : STARK_FLAG_PARENT2_REVERSE))
		current->direction ^= 1;

	current->offset = child_offset;
	// current->transition_char = child_offset;
}
// 
// static inline void stark_update_starknode_navigator_down_char(struct starknode_navigator_s* const __restrict current, int_fast8_t ch) {
// 	offset_t offset = current->node->child[current->direction][ch];
// 	
// 	current->depth++;
// 		
// 	current->node = current->stark->level[current->depth] + child_offset;
// 
// 	if (current->node->flags & ((current->node->uplink[0] == current->offset) ? STARK_FLAG_PARENT1_REVERSE : STARK_FLAG_PARENT2_REVERSE))
// 		current->direction ^= 1;
// 
// 	current->offset = child_offset;
// 	current->transition_char = child_offset;
// }

static inline void stark_navigate_node_go_up(struct starknode_navigator_s* const __restrict current) {
	current->offset = current->node->uplink[current->direction ^ 1];
	
	if (current->node->flags & (current->direction ? STARK_FLAG_PARENT1_REVERSE : STARK_FLAG_PARENT2_REVERSE)) {
		current->direction ^= 1;
	}
	
	current->node = current->stark->level[--current->depth] + current->offset;
	
	// if (STARK_NODE_IS_PALINDROME(*(current->node))) {
	// 	// There is only one way out of a oalindrome
	// 	current->direction = 0;
	// }
}

/*

i tcga
i AGCA
i tgac

print

test AGC f A

test CGA r A

test GCA r T

*/

static inline int stark_navigate_neighbour(struct starknode_navigator_s* current, int_fast8_t ch) {
	stark_navigate_node_go_up(current);
		
	offset_t child_offset;
	
	if (!OFFSET_VALID(child_offset = current->node->child[current->direction][ch])) {
		return -1;
	}
	
	stark_update_starknode_navigator_down(current, child_offset);
	
	return 0;
}

static inline int stark_navigate_shared_link_pair(offset_t * result[2], struct starknode_navigator_s* const __restrict current, int_fast8_t ch) {
	result[0] = current->node->child[current->direction] + ch;
	
	int_fast8_t char_to_me = current->node->flags & (current->direction ? STARK_PARENT1_LINK_MASK : STARK_PARENT2_LINK_MASK);
	if (current->direction) {
		char_to_me &= STARK_PARENT1_LINK_MASK;
		char_to_me >>= 4;
	} else
		char_to_me &= STARK_PARENT2_LINK_MASK;
	
	if (stark_navigate_neighbour(current, ch))
		return -1;
	
	result[1] = current->node->child[current->direction ^ 1] + char_to_me;
	
	return 0;
}

static inline void stark_navigate_advance(struct starknode_navigator_s* const __restrict current, int_fast8_t ch) {
	//char ch = edgeid_single_index[edge_mask];
	offset_t child_offset;

	while (!OFFSET_VALID(child_offset = current->node->child[current->direction][ch])) {
		
		stark_navigate_node_go_up(current);
		
	}

	stark_update_starknode_navigator_down(current, child_offset);
	
}

static inline offset_t * stark_navigate_parent_offset_to_self(stark_t * const stark, depth_t depth, offset_t offset, int_fast8_t direction) {
	starknode_t * node = stark->level[depth] + offset;
	// offset_t parent_offset 
	offset = node->uplink[direction];
	flags_t flags = node->flags;
	// parent = stark->level[--depth] + parent_offset;


	int_fast8_t ch = flags;
	if (!direction) {
		ch &= STARK_PARENT1_LINK_MASK;
		ch >>= 4;
	} else
		ch &= STARK_PARENT2_LINK_MASK;

	if ((flags & (direction ? STARK_FLAG_PARENT2_REVERSE : STARK_FLAG_PARENT1_REVERSE)))
		direction ^= 1;
	
	return stark->level[depth-1][offset].child[direction] + ch;
}

static inline void stark_navigate_loop_forall(struct starknode_navigator_s* const __restrict start, int direction, int (*callback) (struct starknode_navigator_s *, void *), void * aux) {
	start->direction = 0;
	for (; start->depth && start->depth < HARD_MAXDEPTH; start->depth+=direction, start->offset = 1) {
		for (; start->offset <= start->stark->size[start->depth]; start->offset++) {
			if (OFFSET_VALID(start->offset) && ((start->node = start->stark->level[start->depth] + start->offset)->coverage || start->depth == 1)) {
				// start->node = start->stark->level[start->depth] + start->offset;
				if (callback(start, aux))
					return;
			}
		}
	}
}

/*
static inline void stark_navigate_loop_forall_parallel_depth(struct starknode_navigator_s* const __restrict start, int (*callback) (struct starknode_navigator_s *, void *), void * aux) {
	
	#pragma omp for schedule(static,1)
	for (; start->depth < HARD_MAXDEPTH; start->depth+=start->direction, start->offset = 1) {
		for (; start->offset <= start->stark->size[start->depth]; start->offset+=start->direction) {
			if (OFFSET_VALID(start->offset)) {
				start->node = start->stark->level[start->depth] + start->offset;
				if (callback(start, aux))
					return;
			}
		}
	}
}
*/

// 
// static inline void stark_navigate_dfs(struct starknode_navigator_s* const result
// 	, struct starknode_navigator_s* const start
// 	, int advance_selector(starknode_t *, starknode_t *, void *)
// 	, void * advance_selector_arg
// 	, int return) {
// 	
// }


typedef int (*stark_navigate_dfs_test_callback_t) (struct starknode_navigator_s navigator[2], void * aux);

// static inline void stark_navigate_dfs_forward(
// 		  depth_t mink
// 		, struct starknode_navigator_s * navigator
// 		, stark_navigate_dfs_test_callback_t test_callback
// 		// , stark_navigate_dfs_use_callback_t use_callback
// 		, void * callback_aux) {
// 	struct starknode_navigator_s nav[2];
// 	struct starknode_navigator_s parent;
// 	
// 	nav[0] = *navigator;
// 	
// 	if (nav[0].depth >= mink) {
// 		nav[1] = nav[0];
// 		stark_navigate_node_go_up(nav +1);
// 
// 		if (nav[1].depth >= mink && test_callback(nav, callback_aux)) {
// 			list_insert_fast(&stack, nav[1]);
// 		}
// 		
// 		// neighbours
// 		parent = nav[1];
// 		int i;
// 		for (i = 0; i < 4; i++) {
// 			offset_t child_offset;
// 			if (OFFSET_VALID(child_offset = parent.node->child[parent.direction][i])) {
// 				nav[1] = parent;
// 				stark_update_starknode_navigator_down(nav + 1, child_offset);
// 
// 				if (test_callback(nav, callback_aux)) {
// 					list_insert_fast(&stack, nav[1]);
// 				}
// 			}
// 		}
// 	}
// 
// 	// children
// 	int i;
// 	for (i = 0; i < 4; i++) {
// 		offset_t child_offset;
// 		if (OFFSET_VALID(child_offset = nav->node->child[nav[0].direction][i])) {
// 			nav[1] = nav[0];
// 
// 			stark_update_starknode_navigator_down(nav + 1, child_offset);
// 
// 			if (test_callback(nav, callback_aux)) {
// 				list_insert_fast(&stack, nav[1]);
// 			}
// 		}
// 
// 	}
// }

typedef list_t(struct starknode_navigator_s) stark_navigate_dfs_stack_t;

static inline size_t stark_navigate_dfs(
		  stark_navigate_dfs_stack_t * stack
		, depth_t mink
		// , struct starknode_navigator_s * navigator
		, stark_navigate_dfs_test_callback_t test_callback
		// , stark_navigate_dfs_use_callback_t use_callback
		, void * callback_aux) {
	
	struct starknode_navigator_s nav[2]; // = {*navigator, *navigator};
	
	// if (!test_callback(nav, callback_aux))
	// 	return 0;
	
	struct starknode_navigator_s parent;
	size_t nodes_used = 1;
	
	// struct starknode_navigator_s * stack = stark_alloc_pages(8, 0);
	
	// list_t(struct starknode_navigator_s) stack;
	// 
	// list_init(&stack);
	
	// stack->size = 0;
	
	// if (navigator) {
	// 	list_insert_fast(stack, *navigator);
	// }
	
	while (stack->size) {
		nav[0] = stack->list[--(stack->size)];
		// if (use_callback)
		// 	use_callback(&nav, callback_aux);
		nodes_used++;
		
		
		// children
		int i;
		for (i = 0; i < 8; i++) {
			offset_t child_offset;
			if (OFFSET_VALID(child_offset = nav->node->child[0][i])) {
				nav[1] = nav[0];
				nav[1].direction = i >> 2;

				stark_update_starknode_navigator_down(nav + 1, child_offset);

				if (test_callback(nav, callback_aux)) {
					list_insert_fast(stack, nav[1]);
				}
			}

		}
		
		// parents
		if (nav[0].depth >= mink) {
			parent = nav[0];
			stark_navigate_node_go_up(&parent);
			// neighbours
			// int i;
			for (i = 0; i < 4; i++) {
				offset_t child_offset;
				if (OFFSET_VALID(child_offset = parent.node->child[parent.direction][i])) {
					nav[1] = parent;
					stark_update_starknode_navigator_down(nav + 1, child_offset);

					if (test_callback(nav, callback_aux)) {
						list_insert_fast(stack, nav[1]);
					}
				}
			}
			
			parent = nav[0];
			parent.direction ^= 1;
			stark_navigate_node_go_up(&parent);
			// neighbours
			// int i;
			for (i = 0; i < 4; i++) {
				offset_t child_offset;
				if (OFFSET_VALID(child_offset = parent.node->child[parent.direction][i])) {
					nav[1] = parent;
					stark_update_starknode_navigator_down(nav + 1, child_offset);

					if (test_callback(nav, callback_aux)) {
						list_insert_fast(stack, nav[1]);
					}
				}
			}
			
			if (parent.depth >= mink) {
				nav[1] = parent;

				if (test_callback(nav, callback_aux)) {
					list_insert_fast(stack, nav[1]);
				}

				nav[1] = nav[0];
				stark_navigate_node_go_up(nav +1);

				if (test_callback(nav, callback_aux)) {
					list_insert_fast(stack, nav[1]);
				}
			}
			
			
		}

	}
	
	// list_free(stack);
	
	return nodes_used;
	
}



#endif /* end of include guard: STARK_NAVIGATE_H_TNO9P00H */
