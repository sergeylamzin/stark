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
// #error "This file has to be compiled with -DHIERARCHIAL_ASSEMBLY"

#ifndef STARK_ASSEMBLE_HIERARCHIAL_H_SYNHZO4C
#define STARK_ASSEMBLE_HIERARCHIAL_H_SYNHZO4C

#include <string.h>
#include "stark.h"
#include "stark_alloc.h"
#include "stark_navigate.h"
#include "list.h"
#include "hashmap.h"
#include <sys/time.h>
#include <sys/resource.h>

#include "circular_ring_buffer.h"

// 
// struct stark_hierarchical_assembler_group_s {
// 	size_t id;
// 	double internal_link_strength;
// 	
// 	struct {
// 		size_t id;
// 		double link_strength;
// 	} * neighbours;
// 	
// };
// 

// 
// struct pairing_s {
// 	struct pairing_s * next;
// 	double link_strength;
// 	size_t group[2];
// 	struct starknode_navigator_s hook;
// };
// 
// static inline void add_pairing(struct pairing_s * start, struct pairing_s * new_pairing) {
// 	
// 	while (start->next && start->next->link_strength < new_pairing->link_strength) 
// 		start = start->next;
// 	
// 	new_pairing->next = start->next;
// 	start->next = new_pairing;
// }

// #define stark_hierarchical_assembler_broken_group

typedef uint64_t hirarchial_assembler_group_t;

struct stark_hierarchical_assembler_neighbour_s {
	double link_strength;
	hirarchial_assembler_group_t id;
};

struct stark_hierarchical_assembler_neighbours_global_s {
	double link_strength;
	hirarchial_assembler_group_t ids[2];
};

#define STARK_HIERARCHIAL_ASSEMBLER_GROUP_FLAG_CYCLIC 0x2
#define stark_hierarchical_assembler_group_mark_cyclic(_group) ((_group).flags |= STARK_HIERARCHIAL_ASSEMBLER_GROUP_FLAG_CYCLIC)

struct stark_hierarchical_assembler_group_s {
	hirarchial_assembler_group_t group_id;
	size_t num_nodes;
	double internal_link_strength;
	double mass;
	ssize_t contig_length;
	struct starknode_navigator_s endpoint[2];
	// list_t(struct stark_hierarchical_assembler_neighbour_s) neighbours;
	signed char contains_palindrome;
	signed char flags;
};


struct stark_hierarchical_assembler_test_and_merge_groups_dispatch_s {
	struct circular_ring_buffer_s buffer;
	hashmap_uint64_t dispatch_map;
};

struct stark_hierarchical_assembler_s {
	stark_t * stark;
	depth_t minK;
	volatile list_t(struct stark_hierarchical_assembler_group_s) groups;
	volatile list_t(hirarchial_assembler_group_t) recycle_groups;
	list_t(struct stark_hierarchical_assembler_neighbours_global_s) neighbours;
	volatile size_t merged;
	struct stark_hierarchical_assembler_test_and_merge_groups_dispatch_s * volatile dispatch;
	hashmap_uint64_t locked_groups_map;
	int threads_idle;
	// struct pairing_s pairings;
};



// int stark_hierarchical_assembler_init_group_neighbours_forall_cb(struct starknode_navigator_s * current, void * _hierarchical_assembler) ;
// 
// static inline void stark_hierarchical_assembler_init_group_neighbours(struct stark_hierarchical_assembler_s * hierarchical_assembler) {
// 	struct starknode_navigator_s start = {
// 		  .stark = hierarchical_assembler->stark
// 		, .depth = hierarchical_assembler->minK
// 		, .offset = 1
// 	};
// 	
// 	stark_navigate_loop_forall(&start, 1, stark_hierarchical_assembler_init_group_neighbours_forall_cb, hierarchical_assembler);
// }


void stark_hierarchical_assembler_init_group_neighbours_threadsafe(struct stark_hierarchical_assembler_s * hierarchical_assembler);


// void stark_hierarchical_assembler_get_sequences(stark_t * stark, depth_t mink);
void stark_hierarchical_assembler_init (stark_t * stark, struct stark_hierarchical_assembler_s * hierarchical_assembler, depth_t minK);
size_t stark_hierarchical_assembler_test_and_merge_groups(struct stark_hierarchical_assembler_s * hierarchical_assembler, size_t max_rounds);
size_t stark_hierarchical_assembler_print_all_group_contigs(FILE * fp, struct stark_hierarchical_assembler_s * hierarchical_assembler, size_t bounds[2]);
size_t stark_hierarchical_assembler_test_and_merge_groups_openmp_master(struct stark_hierarchical_assembler_s * hierarchical_assembler, size_t max_rounds);
size_t stark_hierarchical_assembler_test_and_merge_groups_openmp(struct stark_hierarchical_assembler_s * hierarchical_assembler, size_t max_rounds);
size_t stark_hierarchical_assembler_test_and_merge_groups_openmp_v2(struct stark_hierarchical_assembler_s * hierarchical_assembler, size_t max_rounds);

#define stark_hierarchical_assembler_find_group_by_id(_ha, _id) ((_ha)->groups.list + (_id))

// static inline void stark_hierarchical_assembler_group_add_neighbour(struct stark_hierarchical_assembler_group_s * group, size_t neighbour, double link_strength) {
// 	size_t i;
// 	for (i = 0; i < group->neighbours.size; i++) {
// 		if (group->neighbours.list[i].id == neighbour) {
// 			if (group->neighbours.list[i].link_strength > link_strength)
// 				group->neighbours.list[i].link_strength = link_strength;
// 			return;
// 		}
// 	}
// 	struct stark_hierarchical_assembler_neighbour_s * neighbour_s = group->neighbours.list + list_new_empty(&(group->neighbours));
// 	neighbour_s->id = neighbour;
// 	neighbour_s->link_strength = link_strength;
// }


static inline size_t print_timediff(char * buffer, struct rusage *a, struct rusage *b) {
	struct timeval t;
	timeradd(&a->ru_utime, &a->ru_stime, &t);
	timersub(&t, &b->ru_utime, &t);
	timersub(&t, &b->ru_stime, &t);
	
	size_t sec = t.tv_sec;
	int h, min;
	if (sec) {
		if ((h = (sec / (60 * 60)))) {
			min = (sec % (60 * 60)) / 60;
			sec %= 60;
			return sprintf(buffer, "%d:%d:%zu", h, min, sec);
		} else if ((min = sec / 60)) {
			sec %= 60;
			return sprintf(buffer, "%d:%zu", min, sec);
		} else 
			return sprintf(buffer, "%zusec", sec);
	} else {
		return sprintf(buffer, "%lums", (unsigned long)t.tv_usec / 1000);
	}
}


#define stark_hierarchical_assembler_group_add_neighbours(_ha, _n) list_insert(&((_ha)->neighbours), *(_n))


#endif /* end of include guard: STARK_ASSEMBLE_HIERARCHIAL_H_SYNHZO4C */

#endif
