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


/*
static void insert_one_sequence_phase1(void * aux, stark_t* const stark, struct uplevelcache* const ___RESTRICT uplevel, const depth_t depth) {
	
	const int_fast8_t forwardchar_index = uplevel[depth-1].symbol_index;
	const int_fast8_t backchar_complement_index = FOURSYMBOL_COMPLEMENT(uplevel->symbol_index);
	const int_fast8_t uplevel0reverse = uplevel->reverse;
	const int_fast8_t uplevel1reverse = uplevel[1].reverse;
		
	unsigned char rindex, palindrome;
	
	rindex = 0;
	palindrome = 0;
	{
		size_t i,j = depth-1;
		char d;
		for (i = 0; i < depth-1; i++) {
			d = uplevel[i].symbol_index - FOURSYMBOL_COMPLEMENT(uplevel[j].symbol_index);
			if (d < 0) {
				goto reversing_end;
			}
			else if (d > 0) {
				rindex = 1; // reverse
				goto reversing_end;
			}
			j--;
		}

		palindrome = 1; // sequence is a palindrome
		
		
	}
	reversing_end:
	;
	
	starknode_t * backparent, * forwardparent;
	{
		starknode_t* const tmp = stark->level[depth-1];
		backparent = tmp + uplevel->offset;
		forwardparent = tmp + uplevel[1].offset;
	}
	
	#if defined(UBIGRAPHAPI_H) || defined(UBIGRAPITHREADED_H)
	{
		char forwardchar = seq4[depth-1];
		char backchar = uplevel->symbol_index;
		
		if (backparent->edges & edge_index[uplevel0reverse][forwardchar] || 
		forwardparent->edges & edge_index[uplevel1reverse ^ 1][complement_index[backchar]])
			;
		else {
			//DEBUG_MSG("Inserting edge between %s and %s", backparent->sequence, forwardparent->sequence);
			stark_ubigraph_new_edge_stark(backparent->ubigraph.id, forwardparent->ubigraph.id, 2 * (depth - 1));
			// if (edge > 0)
			// 	stark_ubigraph_change_edge_style(edge, (2 * depth) - 1);
		}
	}
	
	#endif
	
	
	
	
	offset_t offset;
	
	// int offset_aquisition_method = 0; // direct
	
	starknode_t* ___RESTRICT newnode;
	// synched node creation
	
	{
		offset_t * offsets[4] = {
			&(backparent->child[uplevel0reverse][forwardchar_index]),
			&(forwardparent->child[uplevel1reverse ^ 1][backchar_complement_index]),
			uplevel->palindrome ? &(backparent->child[uplevel0reverse ^ 1][forwardchar_index]) : NULL,
			uplevel[1].palindrome ? &(forwardparent->child[uplevel1reverse][backchar_complement_index]) : NULL
			};
		
		volatile offset_t* lowptr = offsets[0];
		
		if ((uintptr_t)offsets[1] < (uintptr_t)lowptr) 
			lowptr = offsets[1];
		
		if (offsets[2] && (uintptr_t)offsets[2] < (uintptr_t)lowptr) 
			lowptr = offsets[2];
		
		if (offsets[3] && (uintptr_t)offsets[3] < (uintptr_t)lowptr) 
			lowptr = offsets[3];
		
		offset = *lowptr;
	
		if (!offset || offset == locked_offset) {
		
			
			// {
			// 	fprintf(stderr,"creating new node, its offset targets are: ");
			// 	int i;
			// 	for (i = 0; i < 4; i++)
			// 		fprintf(stderr, "%p, ", offsets[i]);
			// 	fprintf(stderr,"\n lowptr = %p\n", lowptr);
			// }
		
			// offset_t locked_offset = 0;
			// 		locked_offset--;
			if (__sync_bool_compare_and_swap(lowptr, (offset_t)0, locked_offset)) {
				offset = stark_create_node(stark, depth);
				
				// char buffer[1024];
				// size_t buffer_length = 0;
				// 
				// buffer_length += sprintf(buffer, "Created new node [%d,%" PRI_OFFSET "]", depth, offset);
			
				if ((uintptr_t)offsets[0] != (uintptr_t)lowptr) {
					// buffer_length += sprintf(buffer + buffer_length, ", %p = %" PRI_OFFSET, offsets[0], offset);
					*(offsets[0]) = offset;
				}
					
			
				if ((uintptr_t)offsets[1] != (uintptr_t)lowptr) {
					// buffer_length += sprintf(buffer + buffer_length, ", %p = %" PRI_OFFSET, offsets[1], offset);
					*(offsets[1]) = offset;
				}
				
				if (offsets[2] && (uintptr_t)offsets[2] != (uintptr_t)lowptr) {
					// buffer_length += sprintf(buffer + buffer_length, ", %p = %" PRI_OFFSET, offsets[2], offset);
					*(offsets[2]) = offset;
				}
				
				if (offsets[3] && (uintptr_t)offsets[3] != (uintptr_t)lowptr) {
					// buffer_length += sprintf(buffer + buffer_length, ", %p = %" PRI_OFFSET, offsets[3], offset);
					*(offsets[3]) = offset;
				}
							
				// printf("%s, %p = %" PRI_OFFSET "\n", buffer, lowptr, offset);
				*lowptr = offset;
				
				// init new node
				
				// offset_aquisition_method = 1; // I created this Node
			
				newnode = stark->level[depth] + offset;
				newnode->depth = depth;
				newnode->uplink[rindex] = uplevel->offset;
				newnode->uplink[rindex ^ 1] = uplevel[1].offset;

				{
					flags_t flags = (forwardchar_index << (rindex ? 0 : 4)) | (backchar_complement_index << (rindex ? 4 : 0));

					if (uplevel0reverse != rindex)
						flags |= rindex ? STARK_FLAG_PARENT2_REVERSE : STARK_FLAG_PARENT1_REVERSE;
					if (uplevel1reverse != rindex)
						flags |= rindex ? STARK_FLAG_PARENT1_REVERSE : STARK_FLAG_PARENT2_REVERSE;

					newnode->flags = flags;
				}


				#if defined(UBIGRAPHAPI_H) || defined(UBIGRAPITHREADED_H)
				{
					vertex_id_t vertex_id;
					stark_ubigraph_new_vertex_w_id(vertex_id = (++stark->ubicount));
					//DEBUG_MSG("Created vertex with ID %d",vertex_id);
					newnode->ubigraph.id = vertex_id;
					stark_ubigraph_change_vertex_style(vertex_id, depth);
					char graphsequence[256];
					// memcpy(graphsequence, sequence, depth);
					int i;
					for (i = 0; i < depth; i++)
						graphsequence[i] = nucleotides[uplevel[i].symbol_index];
					graphsequence[depth] = 0;
					stark_ubigraph_set_vertex_attribute(vertex_id, "label", graphsequence);

					if (palindrome)
						stark_ubigraph_set_vertex_attribute(vertex_id, "shape", "torus");

					newnode->ubigraph.edge[0] = stark_ubigraph_new_edge_stark(backparent->ubigraph.id, vertex_id, (2 * depth) - 1);
					newnode->ubigraph.edge[1] = stark_ubigraph_new_edge_stark(forwardparent->ubigraph.id, vertex_id, (2 * depth) - 1);
				}
				#endif
			
			
			} else {
				// we have a synchronisation problem, someone else was faster
			
				// current solution: spinlock
				// DEBUG_MSG("Thread %lu Going Spinlock on adress %p", pthread_self(), lowptr);
				while ((offset = *lowptr) == locked_offset);
				// DEBUG_MSG("Thread %lu Exiting Spinlock", pthread_self());
				
				newnode = stark->level[depth] + offset;
			}
		
		} else
			newnode = stark->level[depth] + offset;
	}
	
		
	#ifndef NOCHECKS
	if (offset != *(volatile offset_t *)(forwardparent->child[uplevel1reverse ^ 1] + backchar_complement_index) || offset != *(volatile offset_t *)(backparent->child[uplevel0reverse] + forwardchar_index) ) { // data inconsistency
		char forwardchar = nucleotides[uplevel[depth-1].symbol_index];
		char backchar = nucleotides[uplevel->symbol_index];
		
		char buffer[4048 * 2] = "was inserting k-mer ";
		int buflen = sizeof("was inserting k-mer ") - 1;
		// memcpy(buffer+buflen, sequence, depth);
		int i;
		for (i = 0; i < depth; i++) {
			buffer[buflen++] = nucleotides[uplevel[i].symbol_index];
		}
		// buflen += depth;
		// buffer[buflen++] = '\n';
		buflen += sprintf(buffer + buflen, "\n[%d,%d] ", depth-1, uplevel[0].offset);
		buflen += stark_print_node(buffer + buflen, stark, depth-1, uplevel[0].offset);
		buflen += sprintf(buffer + buflen, "\n[%d,%d] ", depth-1, uplevel[1].offset);
		buflen += stark_print_node(buffer + buflen, stark, depth-1, uplevel[1].offset);
		// fprintf(stderr, "%s", buffer);
		
		stark_print(stark, 128);
		
		// buflen += sprintf(buffer + buflen, "\noffset_aquisition_method = %d\nlowptr = %p\n", offset_aquisition_method, lowptr);
		
		
		// int i;
		// for (i = 0; i < 4; i++)
		// 	buflen += sprintf(buffer + buflen, "offsets[%d] = %p\n", i, offsets[i]);

		CRITICAL_ERROR("%sData Inconsistency.\noffset = %" PRI_OFFSET "\nforwardparent->child[uplevel1reverse ^ 1][backchar_complement_index] = %" PRI_OFFSET "\nbackparent->child[uplevel0reverse][forwardchar_index] = %" PRI_OFFSET "\n%s child on node [%d, %" PRI_OFFSET "] for %c does not match %s child on node [%d, %" PRI_OFFSET "] for %c.", buffer, offset, forwardparent->child[uplevel1reverse ^ 1][backchar_complement_index], backparent->child[uplevel0reverse][forwardchar_index], uplevel0reverse ? "Reverse" : "Forward", depth-1, uplevel->offset, forwardchar,(uplevel1reverse ^ 1) ? "reverse" : "forward", depth -1,uplevel[1].offset, complement_index[backchar]);
	}
	#endif
	
	{
		// for backparent
		uint8_t new_edges;
		if (uplevel->palindrome) {
			new_edges = 1 << (3 - forwardchar_index);
			new_edges |= new_edges << 4;
		}
		else
			new_edges = 1 << ((uplevel0reverse ? 3 : 7) - forwardchar_index); 
		
		#ifdef THREADED_IMPORT
		if ((backparent->edges & new_edges) != new_edges) 
			__sync_fetch_and_or(&(backparent->edges), new_edges);
		#else
		backparent->edges |= new_edges;
		#endif
	}
	
	{
		// for forwardparent
		uint8_t new_edges;
		if (uplevel[1].palindrome) {
			new_edges = 1 << (3 - backchar_complement_index);
			new_edges |= new_edges << 4;
		}
		else
			new_edges = 1 << ((uplevel1reverse ? 7 : 3) - backchar_complement_index); 
		
		#ifdef THREADED_IMPORT
		if ((forwardparent->edges & new_edges) != new_edges)
			__sync_fetch_and_or(&(forwardparent->edges), new_edges);
		#else
		forwardparent->edges |= new_edges;
		#endif
	}

	
	#ifdef THREADED_IMPORT
	stark_increment_coverage_cached(stark, &(newnode->coverage), depth, offset);
	#else
	newnode->coverage++;
	#endif
	
	uplevel->offset = offset;
	uplevel->reverse = rindex;
	uplevel->palindrome = palindrome;
	
}


static void insert_one_sequence_phase2(void * aux, stark_t* const stark, struct uplevelcache* const ___RESTRICT uplevel, const depth_t depth) {
	
	// const int_fast8_t forwardchar_index = uplevel[depth-1].symbol_index;
	// const int_fast8_t backchar_complement_index = FOURSYMBOL_COMPLEMENT(uplevel->symbol_index);
	// const int_fast8_t uplevel0reverse = uplevel->reverse;
	// const int_fast8_t uplevel1reverse = uplevel[1].reverse;
		
	unsigned char rindex = 0; //, palindrome;
	
	// rindex = 0;
	// palindrome = 0;
	{
		size_t i,j = depth-1;
		char d;
		for (i = 0; i < depth-1; i++) {
			d = uplevel[i].symbol_index - FOURSYMBOL_COMPLEMENT(uplevel[j].symbol_index);
			if (d < 0) {
				break;
			}
			else if (d > 0) {
				rindex = 1; // reverse
				break;
			}
			j--;
		}

		// palindrome = 1; // sequence is a palindrome
		
		
	}
	
	// const starknode_t * const backparent = stark->level[depth-1] + uplevel->offset;
	
	offset_t offset = stark->level[depth-1][uplevel->offset].child[uplevel->reverse][uplevel[depth-1].symbol_index];
	
	// int offset_aquisition_method = 0; // direct
	
	// synched node creation
	
	if (offset) {
		starknode_t* ___RESTRICT newnode = newnode = stark->level[depth] + offset;
		#ifdef THREADED_IMPORT
		stark_increment_coverage_cached(stark, &(newnode->coverage), depth, offset);
		#else
		newnode->coverage++;
		#endif
	}
	

	
	uplevel->offset = offset;
	uplevel->reverse = rindex;
	// uplevel->palindrome = palindrome;
	
}
*/

int stark_read_link_statistics_subf(void *_target, stark_t* const stark, struct uplevelcache* const ___RESTRICT uplevel, const depth_t depth) {
	struct stark_read_coverage_statistics_s *target = _target;
	starknode_t * backparent, * forwardparent;
	// if (!OFFSET_VALID(uplevel[1].offset))
	// 	return 1;
	// else {
	// 	if (!OFFSET_VALID(uplevel->offset))
	// 		return 0;
		
		starknode_t* const tmp = stark->level[depth-1];
		backparent = tmp + uplevel->offset;
		forwardparent = tmp + uplevel[1].offset;
		
		if (!backparent->coverage || !forwardparent->coverage)
			return 0;
	// } 
	
	const int_fast8_t forwardchar_index = uplevel[depth-1].symbol_index;
	// const int_fast8_t backchar_complement_index = FOURSYMBOL_COMPLEMENT(uplevel->symbol_index);
	const int_fast8_t uplevel0reverse = uplevel->reverse;
	// const int_fast8_t uplevel1reverse = uplevel[1].reverse;
		
	unsigned char rindex; //, palindrome;
	
	rindex = 0;
	// palindrome = 0;
	{
		size_t i,j = depth-1;
		char d;
		for (i = 0; i < depth-1; i++) {
			d = uplevel[i].symbol_index - FOURSYMBOL_COMPLEMENT(uplevel[j].symbol_index);
			if (d < 0) {
				break;
			}
			else if (d > 0) {
				rindex = 1; // reverse
				break;
			}
			j--;
		}

		// palindrome = 1;
		// sequence is a palindrome
		
		
	}
		
	offset_t offset = backparent->child[uplevel0reverse][forwardchar_index];
	
	
	
	// uplevel->palindrome = palindrome;
	
	intmax_t cov_diff = backparent->coverage;
	cov_diff -= forwardparent->coverage;
	double normalized_diff = cov_diff;
	if (cov_diff >= 0) {
		normalized_diff /= backparent->coverage;
	} else {
		normalized_diff /= -(forwardparent->coverage);
	}
	
	if (normalized_diff >= 1.0) {
		// static int amount = 0;
		DEBUG_MSG("backparent->coverage: %d, forwardparent->coverage %d, cov_diff: %zd, normalized_diff: %f\n"
			, backparent->coverage, forwardparent->coverage
			, cov_diff, normalized_diff);
			
		// if(++amount > 20) {
			exit(-1);
		// }
	}
	
	if (normalized_diff < 0.9) {
		target[depth-1].count++;
		target[depth-1].sum += normalized_diff;
		target[depth-1].square_sum += normalized_diff * normalized_diff;
	}
	
	target[depth-1].hist[(size_t)(normalized_diff * ( sizeof(target[depth-1].hist) / sizeof(target[depth-1].hist[0])))]++;
	
	if (!OFFSET_VALID(offset))
		return 1;
	
	uplevel->offset = offset;
	uplevel->reverse = rindex;
		
	return 0;
}


// struct uplevelcache {
// 	int8_t symbol_index;
// 	int8_t reverse;
// 	int8_t palindrome;
// 	// int8_t type;
// 	// depth_t depth;
// 	//starknode_t* node_p;
// 	int32_t compartment;
// 	void * node;
// 	int64_t offset;
// };


#define UPLEVELCACHE_DEFAULT {.symbol_index = 0, .reverse = 0, .palindrome = 0, .offset = 1, .compartment = 0}

static struct uplevelcache uplevelcache_presets[] = {
	[0 ... 0xFF] = UPLEVELCACHE_DEFAULT,
	[0]   = UPLEVELCACHE_DEFAULT,
	[1]   = {.symbol_index = 1, .reverse = 0, .palindrome = 0, .offset = 2, .compartment = 0},
	[2]   = {.symbol_index = 2, .reverse = 1, .palindrome = 0, .offset = 2, .compartment = 0},
	[3]   = {.symbol_index = 3, .reverse = 1, .palindrome = 0, .offset = 1, .compartment = 0},
	['C'] = {.symbol_index = 1, .reverse = 0, .palindrome = 0, .offset = 2, .compartment = 0},
	['G'] = {.symbol_index = 2, .reverse = 1, .palindrome = 0, .offset = 2, .compartment = 0},
	['T'] = {.symbol_index = 3, .reverse = 1, .palindrome = 0, .offset = 1, .compartment = 0},
	['U'] = {.symbol_index = 3, .reverse = 1, .palindrome = 0, .offset = 1, .compartment = 0},
	['c'] = {.symbol_index = 1, .reverse = 0, .palindrome = 0, .offset = 2, .compartment = 0},
	['g'] = {.symbol_index = 2, .reverse = 1, .palindrome = 0, .offset = 2, .compartment = 0},
	['t'] = {.symbol_index = 3, .reverse = 1, .palindrome = 0, .offset = 1, .compartment = 0},
	['u'] = {.symbol_index = 3, .reverse = 1, .palindrome = 0, .offset = 1, .compartment = 0},
};


void process_sequence_rec(stark_t* const ___RESTRICT stark, struct uplevelcache* uplevel, depth_t depth, ssize_t length, depth_t maxdepth, sequence_process_f process_func, void * aux) {
	for (; depth <= maxdepth && depth <= length; depth++) {
		
		ssize_t offset;
		
		for (offset = 0; offset <= length - depth; offset++ ) {
			if(process_func(aux, stark, uplevel + offset, depth)) {
				if (offset - 1 > length - offset) {
					// struct uplevelcache savedcache = uplevel[offset + 1];
					process_sequence_rec(stark, uplevel + offset + 1, depth, length - offset - 1, maxdepth, process_func, aux);
					length = offset + depth;
					// uplevel[offset + 1] = savedcache;
					break;
				} else {
					process_sequence_rec(stark, uplevel, depth + 1, offset + depth, maxdepth, process_func, aux);
					uplevel += offset + 1;
					length -= offset + 1;
					offset = -1;
				}
			}
		}
		

	}
}



void process_sequence_incomplete(stark_t* const ___RESTRICT stark, const char* const sequence, const size_t length, depth_t maxdepth, sequence_process_f process_func, void * aux) {
	size_t i;
	
	struct uplevelcache* uplevel;
	
	const size_t stackcachelength = (length <= 1024 ? length : 0);
	// stackcachelength = 0;
	struct uplevelcache stack_uplevel[stackcachelength];
	
	if (stackcachelength) {
		uplevel = stack_uplevel;
		memset(stack_uplevel, 0, stackcachelength * sizeof(struct uplevelcache));
	}
	else
		uplevel = calloc(length,sizeof(struct uplevelcache));

	for (i = 0; i < length; i++) {
		uplevel[i] = uplevelcache_presets[sequence[i]];
	}
	
	process_sequence_rec(stark, uplevel, 2, length, maxdepth, process_func, aux);
	
	if(!stackcachelength)
		free(uplevel);
	
	//return 0;
}

void process_sequence(stark_t* const ___RESTRICT stark, const char* const sequence, const size_t length, depth_t maxdepth, sequence_process_f process_func, void * aux) {
	depth_t depth;
	size_t i;
	

	
	struct uplevelcache* uplevel;
	
	const size_t stackcachelength = (length <= 1024 ? length : 0);
	// stackcachelength = 0;
	struct uplevelcache stack_uplevel[stackcachelength];
	
	if (stackcachelength) {
		uplevel = stack_uplevel;
		memset(stack_uplevel, 0, stackcachelength * sizeof(struct uplevelcache));
	}
	else
		uplevel = calloc(length,sizeof(struct uplevelcache));

	for (i = 0; i < length; i++) {
		uplevel[i] = uplevelcache_presets[sequence[i]];
	}
	
	
	if (maxdepth > length)
		maxdepth = length;
	

	
	for (depth = 2; depth <= maxdepth; depth++) {
		
		ssize_t offset;
		//DEBUG_MSG("Uplevel: %s , %s , %s , %s , %s", stark->root[uplevel[0].offset].sequence, stark->root[uplevel[1].offset].sequence, stark->root[uplevel[2].offset].sequence, stark->root[uplevel[3].offset].sequence, stark->root[uplevel[4].offset].sequence);
		
		for (offset = 0; offset <= length - depth; offset++ ) {
			offset += process_func(aux, stark, uplevel + offset, depth);
		}
		

	}
	
	if(!stackcachelength)
		free(uplevel);
	
	//return 0;
}


