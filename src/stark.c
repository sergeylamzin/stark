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

#include "stark_alloc.h"
#include "stark_status.h"
#include "stark_navigate.h"

#ifdef UV_SYSTEM
#include "uv_system.h"
#endif


// void* inline check_memory(void* p, size_t base_size, size_t array_size, size_t required_size) {
// 	if ()
// }


// typedef starknode node;


/*

	A: 65: bx01000001
	C: 67: bx01000011
	G: 71: bx01000111
	T: 84: bx01010100
	U: 85: bx01010101

	A^T :  bx00010101
	A^U :  bx00010100
	C^G :  bx00000100
	
	x & bx00000010
	
	x ^ (bx00000100 | (x & bx00000010 : ))

*/

int STARK_MINOVERLAP = 7;

void get_sequence(char* sequence, stark_t* stark, depth_t depth, offset_t offset);
size_t stark_print_node(char* buffer, stark_t* stark, depth_t depth, offset_t offset);

#define swap_int(a, b) {a ^= b; b ^= a; a ^= b;}


#ifndef __INTEL_COMPILER
#if GCC_VERSION < 430 && defined(__x86_64__)
#define __builtin_bswap32(value) ({__asm__ ("bswapl %0" : "+r" ((uint32_t)value)); value;})
#define __builtin_bswap64(value) ({__asm__ ("bswapq %0" : "+r" ((uint64_t)value)); value;})
#endif
#endif


inline char floor_log2(uint64_t n) {	
	uint_fast16_t l = 0;	
	while (n >>= 1)
		l++;
	
	return l;
}

inline uint64_t reverse_complement_uint64(uint64_t n) {  // 5*5+1 ops
	//n = ((n >> 1) & 0x55555555) | ((n << 1) & 0xaaaaaaaa);
	n = ((n >> 2) & 0x3333333333333333) | ((n << 2) & 0xcccccccccccccccc);
	n = ((n >> 4) & 0x0f0f0f0f0f0f0f0f) | ((n << 4) & 0xf0f0f0f0f0f0f0f0);
	n = ((n >> 8) & 0x00ff00ff00ff00ff) | ((n << 8) & 0xff00ff00ff00ff00);
	n = ((n >> 16) & 0x0000ffff0000ffff) | ((n << 16) & 0xffff0000ffff0000);
	n = ((n >> 32) & 0x00000000ffffffff) | ((n << 32) & 0xffffffff00000000);
	
	return ~n;
}

static inline int_fast8_t  reversing(const char* const ___RESTRICT seq, const size_t length) {
	size_t i,j = length-1;
	char d;
	for (i = 0; i < length; i++) {
		d = seq[i] - complement_index[seq[j]];
		if (d < 0)
			return 1; // continue forward
		else if (d > 0) {
			return -1; // reverse
		}
		j--;
	}
	
	return 0; // sequence is a palindrome
}


// starknode_t* down(stark_t* const stark, starknode_t* node, char direction, char ci) {
// 	long offeset;
// 	if (!(offeset = node->child[direction][ci])) { // node does not exist
// 		offeset = stark_create_node(stark);
// 	}
// 	
// }

static inline void kpm_table(const char* const pattern, size_t length, size_t * const table) { // Knuth-Morris-Patt
	size_t pos = 1;
	size_t cnd = 0;
	
	table[0] = -1;
	table[1] = 0;
	length--;
	
	while (pos < length) {
		if (pattern[pos] == pattern[cnd])
			table[++pos] = ++cnd;
		else if (cnd > 0)
			cnd = table[cnd];
		else 
			table[++pos] = 0;
	}
}

size_t check_cyclic_read(const char* const sequence, const size_t length) { // Adopted KMP
	size_t table[length];
	kpm_table(sequence, length, table);
	
	size_t t, i = 0;
	size_t m = 0;
	
	while (i + m < length) {
		if (sequence[i] == sequence[m + i]) {
			i++;
		} else {
			t = table[i];
			m += i - t;
			i = t > 0 ? t : 0;
		}
	}
	
	return i;
}


#if defined(UBIGRAPHAPI_H) && defined(UBIGRAPITHREADED_H)
static inline edge_id_t stark_ubigraph_new_edge_stark(const vertex_id_t x, const vertex_id_t y, const style_id_t style) {
	if (x == y)
		return -1;
	edge_id_t edge = tuple(x,y);
	// if (
	stark_ubigraph_new_edge_w_id(edge, x, y);
	// )
	// 	return -1;
	// else {
		stark_ubigraph_change_edge_style(edge, style);
		return edge;
	// }
}
#endif


static inline int reversing2(const char* const ___RESTRICT seq4, const size_t length) {
	
	// unsigned char rindex, palindrome;
	
	int_fast32_t i,j = length-1;
	int_fast8_t d;
	for (i = 0; i < length; i++) {
		d = seq4[i] - FOURSYMBOL_COMPLEMENT(seq4[j]);
		if (d < 0)
			return 0; // continue forward
		else if (d > 0) {
			return 1; // reverse
		}
		j--;
	}
	
	return 2; // sequence is a palindrome
	
}

static void insert_one_sequence(stark_t* const stark, const char* const ___RESTRICT _deprecated_sequence, struct uplevelcache* const ___RESTRICT uplevel, const depth_t depth, const char* const ___RESTRICT seq4) {
	
	__builtin_prefetch(&(stark->level[depth-1]), 0, 1);
	
	starknode_t * backparent, * forwardparent;
	{
		starknode_t* const tmp = stark->level[depth-1];
		backparent = tmp + uplevel->offset;
		__builtin_prefetch(backparent, 0, 1);
		forwardparent = tmp + uplevel[1].offset;
		__builtin_prefetch(forwardparent, 0, 1);
	}
	
	// const int_fast8_t forwardchar_index = foursymbol_index[sequence[depth-1]];
	// const int_fast8_t backchar_complement_index = FOURSYMBOL_COMPLEMENT(foursymbol_index[*sequence]);
	const int_fast8_t forwardchar_index = seq4[depth-1];
	const int_fast8_t backchar_complement_index = FOURSYMBOL_COMPLEMENT(*seq4);
	const int_fast8_t uplevel0reverse = uplevel->reverse;
	const int_fast8_t uplevel1reverse = uplevel[1].reverse;
		
	int_fast8_t palindrome = 1;
	int_fast8_t rindex = 0;
	
	{
		int_fast8_t d;
		int_fast32_t length = depth & 1;
		length += depth >> 1;
		const char * seqr = seq4;
		const char * sh = seq4 + depth - 1;

		for (;;) {
			if (__builtin_expect(!length, 0)) {
				// palindrome = 1;
				break;
			}
			length--;
			if (__builtin_expect(d = (*(seqr++) ^ 0x3) - *(sh--), 1)) {
				rindex = (d < 0);
				palindrome = 0;
				break;
			}
		}
	}
	
	// starknode_t * const backparent = stark->level[depth-1] + uplevel->offset;
	// starknode_t * const forwardparent = stark->level[depth-1] + uplevel[1].offset;
	

	
	/*{
		char backseq[256], forwardseq[256];
		get_sequence(backseq, stark, depth-1, uplevel->offset); 
		get_sequence(forwardseq, stark, depth-1, uplevel[1].offset);
		DEBUG_MSG("Local Uplevel: %c%s , %c%s, rindex = %d", uplevel0reverse ? 'r' : 'f', backseq, uplevel1reverse ? 'r' : 'f', forwardseq, rindex);
	}*/
	
	#if defined(UBIGRAPHAPI_H) || defined(UBIGRAPITHREADED_H)
	{
		char forwardchar = seq4[depth-1];
		char backchar = *seq4;
		
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
	
	starknode_t* ___RESTRICT newnode; // = stark->level[depth] + offset;
	// synched node creation
	
	{
		
		volatile offset_t* lowptr;
		offset_t * offset_1 = &(backparent->child[uplevel0reverse][forwardchar_index]);
		__builtin_prefetch(offset_1);
		
		lowptr = offset_1;
		
		offset_t * offset_2 = &(forwardparent->child[uplevel1reverse ^ 1][backchar_complement_index]);
		__builtin_prefetch(offset_2);
		
		if ((uintptr_t)offset_2 > (uintptr_t)lowptr)
			lowptr = offset_2;
		
		offset_t * offset_3 = NULL;
		offset_t * offset_4 = NULL;
		
		if (uplevel->palindrome) {
			offset_3 = &(backparent->child[uplevel0reverse ^ 1][forwardchar_index]);
		}
		if ((uintptr_t)offset_3 > (uintptr_t)lowptr)
			lowptr = offset_3;
		
		if (uplevel[1].palindrome) {
			offset_4 = &(forwardparent->child[uplevel1reverse][backchar_complement_index]);
		}
		if ((uintptr_t)offset_4 > (uintptr_t)lowptr)
			lowptr = offset_4;
		
		
		
		/*
		offset_t * offsets[4] = {
			&(backparent->child[uplevel0reverse][forwardchar_index]),
			&(forwardparent->child[uplevel1reverse ^ 1][backchar_complement_index]),
			uplevel->palindrome ? 
				&(backparent->child[uplevel0reverse ^ 1][forwardchar_index])
				: NULL,
			uplevel[1].palindrome ? 
				&(forwardparent->child[uplevel1reverse][backchar_complement_index])
				: NULL
			};
		
		__builtin_prefetch((void *)(lowptr = offsets[0]));
		
		if ((uintptr_t)offsets[1] < (uintptr_t)lowptr) 
			__builtin_prefetch((void *)(lowptr = offsets[1]));
		
		if (__builtin_expect(offsets[2] && (uintptr_t)offsets[2] < (uintptr_t)lowptr, 0)) 
			lowptr = offsets[2];
		
		if (__builtin_expect(offsets[3] && (uintptr_t)offsets[3] < (uintptr_t)lowptr, 0)) 
			lowptr = offsets[3];
		*/
		
		offset = *lowptr;
	
		if (__builtin_expect(!OFFSET_VALID(offset), 0)) {
		
			
			// {
			// 	fprintf(stderr,"creating new node, its offset targets are: ");
			// 	int i;
			// 	for (i = 0; i < 4; i++)
			// 		fprintf(stderr, "%p, ", offsets[i]);
			// 	fprintf(stderr,"\n lowptr = %p\n", lowptr);
			// }
		
			// offset_t locked_offset = 0;
			// 		locked_offset--;
			if (__builtin_expect(__sync_bool_compare_and_swap(lowptr, (offset_t)0, locked_offset), 1)) {
				offset = stark_create_node(stark, depth);
				
				// char buffer[1024];
				// size_t buffer_length = 0;
				// 
				// buffer_length += sprintf(buffer, "Created new node [%d,%" PRI_OFFSET "]", depth, offset);
				
				/*
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
				*/
				
				
				if ((uintptr_t)lowptr == (uintptr_t)offset_1)
					offset_1 = NULL;
				if ((uintptr_t)lowptr == (uintptr_t)offset_2)
					offset_2 = NULL;
				if ((uintptr_t)lowptr == (uintptr_t)offset_3)
					offset_3 = NULL;
				if ((uintptr_t)lowptr == (uintptr_t)offset_4)
					offset_4 = NULL;

				
				if (offset_1)
					*offset_1 = offset;
				if (offset_2)
					*offset_2 = offset;
				if (offset_3)
					*offset_3 = offset;
				if (offset_4)
					*offset_4 = offset;
				
											
				// printf("%s, %p = %" PRI_OFFSET "\n", buffer, lowptr, offset);
				*lowptr = offset;
				// __sync_synchronize();
				// init new node
				
				// offset_aquisition_method = 1; // I created this Node
			
				newnode = stark->level[depth] + offset;

				{
					flags_t flags = (forwardchar_index << (rindex ? 0 : 4)) | (backchar_complement_index << (rindex ? 4 : 0));

					if (uplevel0reverse != rindex)
						flags |= rindex ? STARK_FLAG_PARENT2_REVERSE : STARK_FLAG_PARENT1_REVERSE;
					if (uplevel1reverse != rindex)
						flags |= rindex ? STARK_FLAG_PARENT1_REVERSE : STARK_FLAG_PARENT2_REVERSE;

					newnode->flags = flags;
				}
				
				newnode->depth = depth;
				newnode->uplink[rindex] = uplevel->offset;
				newnode->uplink[rindex ^ 1] = uplevel[1].offset;
			

				#ifdef DEBUG_SEQUENCES
				#error DEBUG_SEQUENCES no longer supported
				/*
				if (depth < DEBUG_SEQUENCES-1) {
					if (rindex) {
						reverse_complement(newnode->sequence, sequence, depth);
					} else {
						memcpy(newnode->sequence, sequence, depth);
					}
					newnode->sequence[depth] = 0;
				}
				*/
				#endif

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
						graphsequence[i] = nucleotides[seq4[i]];
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
				while ((offset = *lowptr) == locked_offset) {
					#ifdef STARK_COVERAGE_CACHE_TOKEN
					stark_increment_coverage_cached_tokenized_tryflush();
					#endif
				}
				// DEBUG_MSG("Thread %lu Exiting Spinlock", pthread_self());
				// offset = *lowptr;
				// offset_aquisition_method = 2; // i wated for someone else to create it
				
				newnode = stark->level[depth] + offset;
			}
		
		} else
			newnode = stark->level[depth] + offset;
	}
	
		
	#ifndef NOCHECKS
	if (__builtin_expect(offset != *(volatile offset_t *)(forwardparent->child[uplevel1reverse ^ 1] + backchar_complement_index) || offset != *(volatile offset_t *)(backparent->child[uplevel0reverse] + forwardchar_index) , 0)) { // data inconsistency
		char forwardchar = nucleotides[seq4[depth-1]];
		char backchar = nucleotides[*seq4];
		
		char buffer[4048 * 2] = "was inserting k-mer ";
		int buflen = sizeof("was inserting k-mer ") - 1;
		// memcpy(buffer+buflen, sequence, depth);
		int i;
		for (i = 0; i < depth; i++) {
			buffer[buflen++] = nucleotides[seq4[i]];
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
		if (__builtin_expect(uplevel->palindrome, 0)) {
			new_edges = 1 << (3 - forwardchar_index);
			new_edges |= new_edges << 4;
		}
		else
			new_edges = 1 << ((uplevel0reverse ? 3 : 7) - forwardchar_index); 
		
		#ifdef THREADED_IMPORT
		if (__builtin_expect((backparent->edges & new_edges) != new_edges, 0)) 
			__sync_fetch_and_or(&(backparent->edges), new_edges);
		#else
		backparent->edges |= new_edges;
		#endif
	}
	
	{
		// for forwardparent
		uint8_t new_edges;
		if (__builtin_expect(uplevel[1].palindrome, 0)) {
			new_edges = 1 << (3 - backchar_complement_index);
			new_edges |= new_edges << 4;
		}
		else
			new_edges = 1 << ((uplevel1reverse ? 7 : 3) - backchar_complement_index); 
		
		#ifdef THREADED_IMPORT
		if (__builtin_expect((forwardparent->edges & new_edges) != new_edges, 0))
			__sync_fetch_and_or(&(forwardparent->edges), new_edges);
		#else
		forwardparent->edges |= new_edges;
		#endif
	}
		

	
	#ifdef THREADED_IMPORT
	/*
	#ifdef UV_SYSTEM
		uv_gru_atomic_async(&(newnode->coverage), 1, GAMER_AMO_ADD_32_OPWORD);
	#else
		__sync_fetch_and_add((coverage_t *)&(newnode->coverage), (coverage_t)1);
	#endif
	*/
	stark_increment_coverage_cached(stark, &(newnode->coverage), depth, offset);
	#else
	newnode->coverage++;
	#endif
	
	uplevel->offset = offset;
	uplevel->reverse = rindex;
	uplevel->palindrome = palindrome;
	//uplevel->node_p = newnode;
	
	
}

/*
struct uplevelcache uplevel_bases[256] = {
	[0] = {.offset = 1, .reverse = 0, .palindrome = 0},
	[1] = {.offset = 2, .reverse = 1, .palindrome = 0},
	[2] = {.offset = 2, .reverse = 0, .palindrome = 0},
	[3] = {.offset = 1, .reverse = 1, .palindrome = 0},
	
	['A'] = {.offset = 1, .reverse = 0, .palindrome = 0},
	['C'] = {.offset = 2, .reverse = 1, .palindrome = 0},
	['G'] = {.offset = 2, .reverse = 0, .palindrome = 0},
	['T'] = {.offset = 1, .reverse = 1, .palindrome = 0},
	['U'] = {.offset = 1, .reverse = 1, .palindrome = 0},
	
	['a'] = {.offset = 1, .reverse = 0, .palindrome = 0},
	['c'] = {.offset = 2, .reverse = 1, .palindrome = 0},
	['g'] = {.offset = 2, .reverse = 0, .palindrome = 0},
	['t'] = {.offset = 1, .reverse = 1, .palindrome = 0},
	['u'] = {.offset = 1, .reverse = 1, .palindrome = 0},
}
*/


#ifdef DEBUG
#include <sys/time.h>
#include <sys/resource.h>

#ifndef RUSAGE_THREAD
#define RUSAGE_THREAD   1
#endif

size_t total_insert_mis = 0;
size_t total_insertions = 0;
#endif

void insert_sequence(stark_t* const ___RESTRICT stark, const char* const sequence, const size_t length, depth_t maxdepth) {
	
	#if defined(STARK_STATUS_H)
	stark_set_status_action(stark_status_action_inserting_read);
	#endif
	
	depth_t depth;
	size_t i;
	
	// // trim sequence N's
	// 
	// while (*sequence == 'N') {
	// 	sequence++;
	// 	length--;
	// }
	// 
	// while (sequence[length - 1] == 'N') 
	// 	length--;
	
	// proofread the sequence
	
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
	
	char seq4[length+1];
	seq4[length] = '\0';
	// size_t stackcachelength = 1;
	// struct uplevelcache uplevel[length];
	// bzero(uplevel, length * sizeof(struct uplevelcache));
	// 
	// char buffer[length+1];
	// memcpy(buffer, sequence, length);
	// buffer[length] = '\0';
	// DEBUG_MSG("Sequence: %s", sequence);

	
	for (i = 0; i < length; i++) {
		switch(sequence[i]) {
			/*
			case 'a':
			case 'A':
			sequence[i] = 'A';
			uplevel[i].offset = 1;
			break;
			*/
			case 'c':
			case 'C':
			// sequence[i] = 'C';
			uplevel[i].offset = 2;
			seq4[i] = 0x1;
			break;
			case 'g':
			case 'G':
			// sequence[i] = 'G';
			uplevel[i].offset = 2;
			uplevel[i].reverse = 1;
			seq4[i] = 0x2;
			break;
			case 't':
			case 'T':
			case 'u':
			case 'U':
			// sequence[i] = 'T';
			uplevel[i].offset = 1;
			uplevel[i].reverse = 1;
			seq4[i] = 0x3;
			break;
			default: // unknown symbols are treated as Adenine
			// sequence[i] = 'A';
			uplevel[i].offset = 1;
			seq4[i] = 0x0;
			break;
		}
		
		// DEBUG_MSG("sequence[i] = %c, Uplevel[%d] = {.offset = %d, .reverse = %d};", sequence[i], i, uplevel[i].offset, uplevel[i].reverse );
	}
	
	
	if (maxdepth > length)
		maxdepth = length;
	
	// for (i = 0; i <= length; i++ )
	// 	uplevel[i].node_p = &(stark->root);
	
	#if defined(THREADED_IMPORT) && defined(STARK_ALLOC_MALLOC)
	// DEBUG_MSG("Picking up Mutex %d", 1);
	pthread_mutex_lock(stark->mutex + 1);
	#endif
	
	#if defined(STARK_STATUS_H) && defined(DEBUG)
	struct rusage usage_start, usage_end;
	size_t insertions = 0;
	getrusage(RUSAGE_THREAD, &usage_start);
	#endif
	
	for (depth = 2; depth <= maxdepth; depth++) {
		
		offset_t offset;
		//DEBUG_MSG("Uplevel: %s , %s , %s , %s , %s", stark->root[uplevel[0].offset].sequence, stark->root[uplevel[1].offset].sequence, stark->root[uplevel[2].offset].sequence, stark->root[uplevel[3].offset].sequence, stark->root[uplevel[4].offset].sequence);
		
		#if defined(THREADED_IMPORT) && defined(STARK_ALLOC_MALLOC)
		// DEBUG_MSG("Picking up Mutex %d", depth);
		pthread_mutex_lock(stark->mutex + depth);
		#endif
		for (offset = 0; offset <= length - depth; offset++ ) {
			insert_one_sequence(stark, NULL, uplevel + offset, depth, seq4 + offset);
		}
		#if defined(THREADED_IMPORT) && defined(STARK_ALLOC_MALLOC)
		// DEBUG_MSG("Releasing up Mutex %d", depth - 1);
		pthread_mutex_unlock(stark->mutex + depth - 1);
		#endif
		#if defined(STARK_STATUS_H) && defined(DEBUG)
		insertions += length - depth;
		#endif
	}
	
	#if defined(STARK_STATUS_H)
	#if defined(DEBUG)
	getrusage(RUSAGE_THREAD, &usage_end);
	get_stark_thread_status()->total_insert_mis += ((usage_end.ru_utime.tv_sec - usage_start.ru_utime.tv_sec)* 1000000) + usage_end.ru_utime.tv_usec - usage_start.ru_utime.tv_usec;
	get_stark_thread_status()->total_insertions += insertions;
	#endif
	get_stark_thread_status()->reads_inserted++;
	#endif
	
	#if defined(THREADED_IMPORT) && defined(STARK_ALLOC_MALLOC)
	// DEBUG_MSG("Releasing up Mutex %d", depth - 1);
	pthread_mutex_unlock(stark->mutex + depth - 1);
	#endif
	
	if(!stackcachelength)
		free(uplevel);
	
	//return 0;
}

#ifdef THREADED_IMPORT
void stark_insert_sequence_job(struct stark_insert_threaded_s* ___RESTRICT arg) {
	stark_coverage_counter_cache_init();
	while (arg->stark) {
		size_t size = sizeof(struct stark_insert_threaded_s) + arg->length;
		insert_sequence(arg->stark, arg->sequence, arg->length, arg->maxdepth);
		if (size & 0x7) {
			size &= ~0x7;
			size += 0x8;
		}
		arg = (struct stark_insert_threaded_s*)(((char *)arg) + size);
	}
	
}
#endif

/*
#define advance_func void advance(char ch) { \
	node *child_node, *node = stark->level[depth] + offset; \
	offset_t child_offset; \
	\
	while (!(child_offset = node->child[direction][ch])) { \
		offset = node->uplink[direction ^ 1]; \
		if (node->flags & (direction ? STARK_FLAG_PARENT1_REVERSE : STARK_FLAG_PARENT2_REVERSE)) { \
			direction ^= 1; \
		} \
		depth--; \
	} \
	 \
	child_node = stark->level[a->depth = (depth +1)] + child_offset; \
	 \
	if (child_node->flags & ((child_node->uplink[0] == offset) ? STARK_FLAG_PARENT1_REVERSE : STARK_FLAG_PARENT2_REVERSE)) \
		direction ^= 1; \
	 \
	offset = child_offset; \
	depth++; \
}



void advance(struct advancement* a, stark_t* stark, char ch) {
	
	depth_t depth = a->depth;
	offset_t offset = a->offset;
	char direction = a->direction;
	
	node *child_node, *node = stark->level[depth] + offset;
	offset_t child_offset;
	
	while (!(child_offset = node->child[direction][ch])) {
		offset = node->uplink[direction ^ 1];
		if (node->flags & (direction ? STARK_FLAG_PARENT1_REVERSE : STARK_FLAG_PARENT2_REVERSE)) {
			direction ^= 1;
		}
		depth--;
	}
	
	a->offset = child_offset;
	a->depth = depth+1;
	// starknode_t* child_node;
	child_node = stark->level[a->depth = (depth +1)] + child_offset;
	
	
	if (child_node->flags & ((child_node->uplink[0] == offset) ? STARK_FLAG_PARENT1_REVERSE : STARK_FLAG_PARENT2_REVERSE)) {
		a->direction = direction ^ 1;
	} else
		a->direction = direction;
	
	// return;
	
}
*/


#define PRINT_TRAVERSE(format, ...) // fprintf(stderr, format, __VA_ARGS__)

void stark_deep_flag_extracted(stark_t* const stark, depth_t depth, offset_t offset
#if defined(UBIGRAPHAPI_H) || defined(UBIGRAPITHREADED_H)
, const char* color
#endif
) {
	starknode_t* node_p = stark->level[depth] + offset;
	if (depth && !(node_p->flags & STARK_FLAG_EXTRACTED)) {
		node_p->flags |= STARK_FLAG_EXTRACTED;
		
		// struct coverage_histogram_s* coverage_histogram = stark->coverage_histogram;
		// coverage_t coverage = node_p->coverage;
		// 
		// coverage_histogram->nodes_counted--;
		// coverage_histogram->sum_coverage -= coverage;
		// coverage_histogram->hist[coverage]--;
		
		#if defined(UBIGRAPHAPI_H) || defined(UBIGRAPITHREADED_H)
		stark_ubigraph_set_vertex_attribute(node_p->ubigraph.id, "color", color);
		#endif
		stark_deep_flag_extracted(stark, depth-1,  node_p->uplink[0]
		#if defined(UBIGRAPHAPI_H) || defined(UBIGRAPITHREADED_H)
		, "#ffffaa"
		#endif
		);
		stark_deep_flag_extracted(stark, depth-1,  node_p->uplink[1]
		#if defined(UBIGRAPHAPI_H) || defined(UBIGRAPITHREADED_H)
		, "#ffffaa"
		#endif
		);
	}
}

size_t build_seq(unsigned short* junctions, char* seq, stark_t* stark, offset_t current_offset, depth_t current_depth, char run) {

	starknode_t currentnode;
	char current_direction = run;
	size_t seqlen = 0;
	
	unsigned char edges;
	
	currentnode = stark->level[current_depth][current_offset];
	
	while (1) {

		stark_deep_flag_extracted(stark, current_depth, current_offset
		#if defined(UBIGRAPHAPI_H) || defined(UBIGRAPITHREADED_H)
		, "#ffffff"
		#endif
		);

		
		edges = currentnode.edges & (current_direction ? STARK_EDGES_REVERSE_MASK : STARK_EDGES_FORWARD_MASK);
		
		// DEBUG_MSG("0x%02X", edges);
		
		
		char buffer[1024];
		get_sequence(buffer, stark, current_depth, current_offset);
		if (current_direction ^ run ^ 1) {
			reverse_complement_onspot(buffer, current_depth);
		}
		char format[256];
		if (!run)
			snprintf(format, sizeof(format),".%c%zds%cc",'%',50 - seqlen + current_depth,'%');
		else
			snprintf(format, sizeof(format),".%c%zds%cc",'%',50 + seqlen,'%');
		printf(format,buffer,'\n');
		
		
		while (!edges && current_depth >= STARK_MINOVERLAP) {
			current_offset = currentnode.uplink[current_direction ^ 1];
			if (currentnode.flags & (current_direction ? STARK_FLAG_PARENT1_REVERSE : STARK_FLAG_PARENT2_REVERSE)) {
				current_direction ^= 1;
			}
			current_depth--;
			
			currentnode = stark->level[current_depth][current_offset];
			
			if (STARK_NODE_IS_PALINDROME(currentnode))
				break;
			
			edges = currentnode.edges & (current_direction ? STARK_EDGES_REVERSE_MASK : STARK_EDGES_FORWARD_MASK);
			
		}

		if (!edges || (edges & (edges - 1)))
			break;


		char ch = edgeid_single_index[edges];
		// DEBUG_MSG("ch = %d", ch);
		seq[seqlen++] = edgeid_single_char[edges];
		
		
		
		// stark->level[current_depth][current_offset].flags |= STARK_FLAG_EXTRACTED;
		// stark_deep_flag_extracted(stark, current_depth, current_offset
		// #if defined(UBIGRAPHAPI_H) || defined(UBIGRAPITHREADED_H)
		// , "#ffffff"
		// #endif
		// );
		
		// Advance
		{
			offset_t child_offset;

			while (!OFFSET_VALID(child_offset = currentnode.child[current_direction][ch])) {
				current_offset = currentnode.uplink[current_direction ^ 1];
				if (currentnode.flags & (current_direction ? STARK_FLAG_PARENT1_REVERSE : STARK_FLAG_PARENT2_REVERSE)) {
					current_direction ^= 1;
				}
				current_depth--;
				
				currentnode = stark->level[current_depth][current_offset];
				
				
				
				PRINT_TRAVERSE("(%d,%" PRI_OFFSET "), ", current_depth, current_offset);
				
				// stark->level[current_depth][current_offset].flags |= STARK_FLAG_EXTRACTED;
				// #if defined(UBIGRAPHAPI_H) || defined(UBIGRAPITHREADED_H)
				// stark_ubigraph_set_vertex_attribute(stark->level[current_depth][current_offset].ubigraph.id, "color", "#aaffff");
				// #endif
				stark_deep_flag_extracted(stark, current_depth, current_offset
				#if defined(UBIGRAPHAPI_H) || defined(UBIGRAPITHREADED_H)
				, "#aaffff"
				#endif
				);
			}

			current_depth++;
			
			PRINT_TRAVERSE("[%d,%" PRI_OFFSET "], ", current_depth, child_offset);
			
			currentnode = stark->level[current_depth][child_offset];


			if (currentnode.flags & ((currentnode.uplink[0] == current_offset) ? STARK_FLAG_PARENT1_REVERSE : STARK_FLAG_PARENT2_REVERSE))
				current_direction = current_direction ^ 1;

			current_offset = child_offset;
			
			if (currentnode.flags & STARK_FLAG_EXTRACTED) {
				DEBUG_MSG("Encountered an already extracted node at %d %" PRI_OFFSET, current_depth, current_offset);
				return seqlen;
			}
			
		}
		
	}
	
	// if (currentnode.flags & STARK_FLAG_EXTRACTED) {
	// 	DEBUG_MSG("Encountered an already extracted node at %d %" PRI_OFFSET, current_depth, current_offset);
	// }
	
	if (edges) {
		// stark->level[current_depth][current_offset].flags |= STARK_FLAG_JUNCTION;
		// stark->level[current_depth][current_offset].junction_id = (junctions[run] = ++(stark->junctions));
	}
	
	return seqlen;
}

size_t stark_extract_sequence(char* seq, stark_t* stark, depth_t startdepth, offset_t startoffset) {
	size_t seqlen = startdepth;
	
	PRINT_TRAVERSE("[%d,%" PRI_OFFSET "], ", startdepth, startoffset);
	get_sequence(seq, stark, startdepth, startoffset);
	
	unsigned short junctions[2];
	
	seqlen += build_seq(junctions, seq + seqlen, stark, startoffset, startdepth, 0);
	
	seq[seqlen] = 0;
	
	// DEBUG_MSG("%s", seq);
	
	reverse_complement_onspot(seq, seqlen);
	
	// DEBUG_MSG("%s", seq);
	
	seqlen += build_seq(junctions, seq + seqlen, stark, startoffset, startdepth, 1);
	
	seq[seqlen] = 0;
	return seqlen;
	
}

/*
struct assembly_support_s {
	char* sequence;
	coverage_t* coverage;
	size_t memory;
	size_t length;
	stark_t* stark;
	
	struct starknode_id_s current;
	
	struct {
		depth_t depth;
		offset_t offset;
	} start;
	
	int32_t minK;
};
*/

static inline void free_simple_assembly_s(struct simple_assembly_s *s) {
	
	if (s->sequence)
		free(s->sequence);
	if (s->coverage)
		free(s->coverage);
}

static inline void stark_node_go_up(struct starknode_id_s* current) {
	current->offset = current->node.uplink[current->direction ^ 1];
	
	if (current->node.flags & (current->direction ? STARK_FLAG_PARENT1_REVERSE : STARK_FLAG_PARENT2_REVERSE)) {
		current->direction ^= 1;
	}
	current->depth--;
	
	current->node = *(current->node_p = current->stark->level[current->depth] + current->offset);
}

static inline void stark_advance(struct starknode_id_s* current, char ch) {
	//char ch = edgeid_single_index[edge_mask];
	offset_t child_offset;

	while (!OFFSET_VALID(child_offset = current->node.child[current->direction][ch])) {
		// 
		// if (current->node.depth < minK) {
		// 	memset(current, 0, sizeof(*current));
		// 	return;
		// }
		
		stark_node_go_up(current);
		
	}

	current->depth++;
		
	current->node = *(current->node_p = current->stark->level[current->depth] + child_offset);

	if (current->node.flags & ((current->node.uplink[0] == current->offset) ? STARK_FLAG_PARENT1_REVERSE : STARK_FLAG_PARENT2_REVERSE))
		current->direction = current->direction ^ 1;

	current->offset = child_offset;
	
}


// This choice method attempts to choose each nucleotide once (longest k-mer) and evalueattes based on coverage, parent depth, target depth
char stark_choose_check_all_chars_once(struct starknode_id_s *current, int32_t minK) {
	char next_char = 0;
	
	unsigned char edges, edge_mask, edges_tmp, edges_done = 0;
	depth_t current_candidate_parent_depth = 0;
	
	struct starknode_id_s current_candidate;
	current_candidate.node.coverage = 0;
	current_candidate.depth = 0;
	
	// printf("edges = 0x%x\n", edges);
	
	for (;;) {
		
		edges = (current->node.edges & (current->direction ? STARK_EDGES_REVERSE_MASK : STARK_EDGES_FORWARD_MASK));

		if (!current->direction)
			edges >>= 4; 

		edges = edges & (~edges_done);
		
		if (edges) {
			
			// extract edge masks one by one
			for (; edges ; edges = edges_tmp) {
				edges_tmp = edges & (edges - 1);
				edge_mask = edges ^ edges_tmp;
				
				edges_done |= edge_mask;

				// printf("edge_mask = 0x%x\n", edge_mask);

				struct starknode_id_s candidate = *current;

				// printf("Trying Edge %c\n", edgeid_single_char[edge_mask]);

				stark_advance(&candidate, edgeid_single_index[edge_mask]);
				// stark_advance(&candidate, ch, edge_index[current->direction][nucleotides[ch]] & edges ? 1 : minK);

				// char buffer[256];
				// get_sequence(buffer, candidate.stark, candidate.depth , candidate.offset);
				// if (candidate.direction)
				// 	reverse_complement_onspot(buffer, candidate.depth);
				// printf("Trying through [%d,%" PRI_OFFSET "]:%s\n", candidate.depth , candidate.offset, buffer);

				if ( !(candidate.node.flags & STARK_FLAG_EXTRACTED) && (candidate.node.coverage * candidate.depth * current->depth) > (current_candidate.node.coverage * current_candidate.depth * current_candidate_parent_depth)) {
					current_candidate = candidate;
					current_candidate_parent_depth = current->depth;
					next_char = edgeid_single_char[edge_mask];
					// next_char = nucleotides[ch];
				}

			}
		}
		
		if (current->depth >= minK && edges_done != STARK_EDGES_REVERSE_MASK)
			stark_node_go_up(current);
		else
			break;
						
	}
	
	*current = current_candidate;
	
	return next_char;
}



char* stark_assemble_maxcov_direction(struct simple_assembly_s* simple_assembly, stark_t* const stark, depth_t startdepth, offset_t startoffset, int32_t minK, char direction) {
	list_t(char) sequence_l;
	list_init(&sequence_l);
	
	list_t(coverage_t) coverage_l;
	list_init(&coverage_l);
	
	struct starknode_id_s current;
	
	current.stark = stark;
	current.depth = startdepth;
	current.offset = startoffset;
	current.direction = direction;
	current.node = *(current.node_p = stark->level[current.depth] + current.offset);
				
	// starknode_t currentnode;
	
	char next_char;
	
	for (;;) {
				
		// for (;;) {
		// 	edges = (current.node.edges & (current.direction ? STARK_EDGES_REVERSE_MASK : STARK_EDGES_FORWARD_MASK));
		// 	if (edges || current.depth < minK)
		// 		break;
		// 	stark_node_go_up(&current);
		// }
		// 
		// if (!edges) {
		// 	break;
		// }
		
		
		// DEBUG
		// char buffer[256];
		// get_sequence(buffer, current.stark, current.depth , current.offset);
		// if (current.direction)
		// 	reverse_complement_onspot(buffer, current.depth);
		// printf("Going through [%d,%" PRI_OFFSET "]:%s\n", current.depth , current.offset, buffer);
		
		
		stark_deep_flag_extracted(stark, current.depth , current.offset
		#if defined(UBIGRAPHAPI_H) || defined(UBIGRAPITHREADED_H)
		, "#ffffff"
		#endif
		);
		
		
		// offset_t flag_offset = current.offset;
		// depth_t flag_depth = current.depth;
		
		//struct starknode_id_s current_candidate;
		next_char = stark_choose_check_all_chars_once(&current, minK);
		/*
		{
			unsigned char edge_mask, edges_tmp;
			depth_t current_candidate_parent_depth = 0;
			
			unsigned char edges_done = 0;
			
			// printf("edges = 0x%x\n", edges);
			
			for (;;) {
				
				edges = (current.node.edges & (current.direction ? STARK_EDGES_REVERSE_MASK : STARK_EDGES_FORWARD_MASK));

				if (!current.direction)
					edges >>= 4; 

				edges = edges & (~edges_done);
				
				if (edges) {
					
					// extract edge masks one by one
					for (; edges ; edges = edges_tmp) {
						edges_tmp = edges & (edges - 1);
						edge_mask = edges ^ edges_tmp;
						
						edges_done |= edge_mask;

						// printf("edge_mask = 0x%x\n", edge_mask);

						struct starknode_id_s candidate = current;

						// printf("Trying Edge %c\n", edgeid_single_char[edge_mask]);

						stark_advance(&candidate, edgeid_single_index[edge_mask]);
						// stark_advance(&candidate, ch, edge_index[current.direction][nucleotides[ch]] & edges ? 1 : minK);

						// char buffer[256];
						// get_sequence(buffer, candidate.stark, candidate.depth , candidate.offset);
						// if (candidate.direction)
						// 	reverse_complement_onspot(buffer, candidate.depth);
						// printf("Trying through [%d,%" PRI_OFFSET "]:%s\n", candidate.depth , candidate.offset, buffer);

						if ( !(candidate.node.flags & STARK_FLAG_EXTRACTED) && (candidate.node.coverage * candidate.depth * current.depth) > (current_candidate.node.coverage * current_candidate.depth * current_candidate_parent_depth)) {
							current_candidate = candidate;
							current_candidate_parent_depth = current.depth;
							next_char = edgeid_single_char[edge_mask];
							// next_char = nucleotides[ch];
						}

					}
				}
				
				if (current.depth >= minK && edges_done != STARK_EDGES_REVERSE_MASK)
					stark_node_go_up(&current);
				else
					break;
								
			}
				
		}
		*/
		
		// char buffer[256];
		// get_sequence(buffer, current_candidate.stark, current_candidate.depth , current_candidate.offset);
		// if (current_candidate.direction)
		// 	reverse_complement_onspot(buffer, current_candidate.depth);
		// printf("Going through [%d,%" PRI_OFFSET "]:%s\n", current_candidate.depth , current_candidate.offset, buffer);
		
		if (next_char) {
			
			// current = current_candidate;
			
			if (sequence_l.size >= sequence_l.maxlength)
				list_insert_fast(&sequence_l, next_char);
			else
				sequence_l.list[sequence_l.size++] = next_char;
			
			//coverage_t coverage = current.coverage;
			if (coverage_l.size >= coverage_l.maxlength)
				list_insert_fast(&coverage_l, current.node.coverage);
			else
				coverage_l.list[coverage_l.size++] = current.node.coverage;
			
			int64_t i, backtrack = coverage_l.size - current.depth - 1;
			if (backtrack < 0)
				backtrack = 0;
			for (i = coverage_l.size - 2; i >= backtrack; i--) {
				if (coverage_l.list[i] < current.node.coverage)
					coverage_l.list[i] = current.node.coverage;
			}
			
		
		} else {
			break;
		}
	}
	
	
	next_char = 0;	
	list_insert_fast(&sequence_l, next_char);
	
	char* sequence = list_compact(&sequence_l);
	if (simple_assembly) {
		simple_assembly->length = sequence_l.size - 1;
		simple_assembly->sequence = sequence;
		simple_assembly->coverage = list_compact(&coverage_l);
	}
	
	return sequence;
}

void stark_assemble_maxcov(struct simple_assembly_s *assembly, stark_t* const stark, depth_t startdepth, offset_t startoffset, int32_t minK) {
	
	// char[startdepth] node_seq;
	// get_sequence(&node_seq, stark, startdepth, startoffset);
	
	struct simple_assembly_s forward_assembly;
	stark_assemble_maxcov_direction(&forward_assembly, stark, startdepth, startoffset, minK, 0);
	// size_t forward_length = forward_assembly.length;
	
	struct simple_assembly_s reverse_assembly;
	stark_assemble_maxcov_direction(&reverse_assembly, stark, startdepth, startoffset, minK, 1);
	// size_t reverse_length = strlen(reverse_assembly);
	
	// reverse_complement_onspot(reverse_assembly.sequence, reverse_assembly.length);
	
	struct simple_assembly_s joined_assembly;
	joined_assembly.length = startdepth + forward_assembly.length + reverse_assembly.length;
	joined_assembly.sequence = malloc(joined_assembly.length + 1);
	joined_assembly.coverage = malloc(joined_assembly.length * sizeof(coverage_t));
	
	// char* assembly = 
	
	reverse_complement(joined_assembly.sequence, reverse_assembly.sequence, reverse_assembly.length);
	// memcpy(assembly, reverse_assembly, reverse_length);
	get_sequence(joined_assembly.sequence + reverse_assembly.length, stark, startdepth, startoffset);
	memcpy(joined_assembly.sequence + reverse_assembly.length + startdepth, forward_assembly.sequence, forward_assembly.length);
	joined_assembly.sequence[reverse_assembly.length + startdepth + forward_assembly.length] = 0;
	
	int64_t i,j = 0;
	for (i = reverse_assembly.length - 1; i >= 0; i--) {
		joined_assembly.coverage[j++] = reverse_assembly.coverage[i];
	}
	coverage_t startcov = stark->level[startdepth][startoffset].coverage;
	for (i = 0; i < startdepth; i++)
		joined_assembly.coverage[j++] = startcov;
	memcpy(joined_assembly.coverage + j, forward_assembly.coverage, forward_assembly.length * sizeof(coverage_t));
	
	free_simple_assembly_s(&forward_assembly);
	free_simple_assembly_s(&reverse_assembly);
	
	*assembly = joined_assembly;
	
	// return assembly;
}

inline void stark_assemble_single(stark_t* const stark, depth_t minK, depth_t depth, offset_t offset) {
	// DEBUG_MSG("Attempting to assemble starting from %d %" PRI_OFFSET ".", depth, offset);
	struct simple_assembly_s assembly;
	stark_assemble_maxcov(&assembly, stark, depth, offset, minK);
	size_t i;
	coverage_t maxcov = 1;
	for (i = 0; i < assembly.length; i++)
		if (assembly.coverage[i] > maxcov)
			maxcov = assembly.coverage[i];
	
	char* output_quality = malloc(assembly.length + 1);
	for (i = 0; i < assembly.length; i++)
		output_quality[i] = 33 + (int8_t)(93 * assembly.coverage[i] / maxcov);
	output_quality[i] = 0;			
	
	printf("@CONTIG length = %zd\n%s\n+CONTIG max coverage = %d\n%s\n", assembly.length, assembly.sequence, maxcov, output_quality);
	free(output_quality);
	
	// printf("@CONTIG length = %d\n%s\n+CONTIG max coverage = %d\n", assembly.length, assembly.sequence, maxcov);
	// for (i = 0; i < assembly.length; i++)
	// 	printf("%d,", assembly.coverage[i]);
	// printf("\n");
	
	free_simple_assembly_s(&assembly);
}

void stark_assemble_all(stark_t* const stark, depth_t minK) {
	DEBUG_MSG("Assembling all.");
	
	struct coverage_histogram_s* coverage_histogram = stark_statistics(stark, HARD_MAXDEPTH);
	
	depth_t depth;
	
	offset_t offset;
	
	bool done = 0;
	
	while (!done) {
		done = 1;
		coverage_t avg_coverage = coverage_histogram->sum_coverage / (coverage_histogram->nodes_counted ? coverage_histogram->nodes_counted : 1);
		for (depth = minK; (stark->size[depth] - stark->deleted[depth]) > 0; depth++) {
			size_t levelsize = stark->size[depth];
			for (offset = 1; offset <= levelsize ; offset++) {
				if (stark->level[depth][offset].depth && !(stark->level[depth][offset].flags & STARK_FLAG_EXTRACTED) && stark->level[depth][offset].coverage >= avg_coverage) {
					stark_assemble_single(stark, minK, depth, offset);
					done = 0;
				}
			}
		}
	}

}

void stark_extract_all(stark_t* const stark) {
	DEBUG_MSG("Extracting all.");
	
	depth_t depth;
	for (depth = 1; (stark->size[depth] - stark->deleted[depth]) > 0; depth++);
	
	offset_t offset;
	char* buffer = malloc(1048576);
	
	for (depth--; depth; depth--) {
		size_t levelsize = stark->size[depth];
		for (offset = 1; offset <= levelsize ; offset++) {
			if (stark->level[depth][offset].depth && !(stark->level[depth][offset].flags & STARK_FLAG_EXTRACTED)) {
				DEBUG_MSG("Attempting to extract %d %" PRI_OFFSET ".", depth, offset);
				stark_extract_sequence(buffer, stark, depth, offset);
				printf("%s\n", buffer);
			}
		}
	}
	free(buffer);
}

size_t stark_print_node_info(char* buffer, starknode_t * node) {
	//char* buffer = malloc(1024);
	
	size_t string_offset = 0;
	
	string_offset += sprintf(buffer,"cov %d :\t up %c %" PRI_OFFSET " %c %c %" PRI_OFFSET " %c :\tedges :\t",
		node->coverage, node->flags & STARK_FLAG_PARENT1_REVERSE ? '^' : '.',
		node->uplink[0],
		parent_Link_to_char[node->flags & STARK_PARENT1_LINK_MASK],
		node->flags & STARK_FLAG_PARENT2_REVERSE ? '^' : '.',
		node->uplink[1],
		parent_Link_to_char[node->flags & STARK_PARENT2_LINK_MASK]);
	
	#ifdef ENCHAIN_SEQUENCES
	if (node->flags & STARK_FLAG_EXTRACTED)
		string_offset += sprintf(buffer + string_offset,"chain %c%3d -> %3d :\t", (node->sequence.flags & STARK_CHAIN_FLAG_REVERSE) ? '^' : ' ', node->sequence.start, node->sequence.end);
	#endif
	
	int_fast8_t i;
	unsigned char mask = 1 << 7;
	for (i = 0; i < 4; i++) {
		if (node->edges & mask) {
			char c = 'N';
			switch(i) {
				case 0: c = 'A'; break;
				case 1: c = 'C'; break;
				case 2: c = 'G'; break;
				case 3: c = 'T'; break;
			}
			string_offset += sprintf(buffer + string_offset,"f%c.%" PRI_OFFSET " ", c, node->child[0][i]);
		}
		mask >>= 1;
	}
	for (i = 0; i < 4; i++) {
		if (node->edges & mask) {
			char c = 'N';
			switch(i) {
				case 0: c = 'A'; break;
				case 1: c = 'C'; break;
				case 2: c = 'G'; break;
				case 3: c = 'T'; break;
			}
			string_offset += sprintf(buffer + string_offset,"r%c.%" PRI_OFFSET " ", c, node->child[1][i]);
		}
		mask >>= 1;
	}
	
	#ifdef HIERARCHIAL_ASSEMBLY
	string_offset += sprintf(buffer + string_offset,"\thgroup: %lld", node->hgroup);
	
	// if (node->contig_position)
	// 	string_offset += sprintf(buffer + string_offset,"\tpos: %d", node->contig_position);
	#endif
		
	return string_offset;
}

size_t stark_print_node(char* buffer, stark_t* const stark, depth_t depth, offset_t offset) {
	starknode_t* node = stark->level[depth] + offset;
	size_t string_offset = 0;
	
	char sequence[256];
	char reverse[256];
	if (!depth) {
		memcpy(sequence, "root", 5);
		memcpy(reverse, "root", 5);
	} else {
		get_sequence(sequence, stark, depth, offset);
		reverse_complement(reverse, sequence, depth);
		reverse[depth] = 0;
	}
	
	string_offset += sprintf(buffer, "%s %c %s :\t", sequence, node->flags & STARK_FLAG_EXTRACTED ? '*' : '/', reverse);
	string_offset += stark_print_node_info(buffer + string_offset, node);
		
	return string_offset;
	
	/*
	string_offset += sprintf(buffer,"%s %c %s :\tcov %d :\t up %s%" PRI_OFFSET "%c %s%" PRI_OFFSET "%c :\tedges :\t", sequence, node->flags & STARK_FLAG_EXTRACTED ? '*' : '/', reverse, node->coverage, node->flags & STARK_FLAG_PARENT1_REVERSE ? "^" : "", node->uplink[0], parent_Link_to_char[node->flags & STARK_PARENT1_LINK_MASK], node->flags & STARK_FLAG_PARENT2_REVERSE ? "^" : "", node->uplink[1], parent_Link_to_char[node->flags & STARK_PARENT2_LINK_MASK]);
	
	#ifdef ENCHAIN_SEQUENCES
	if (node->flags & STARK_FLAG_EXTRACTED)
		string_offset += sprintf(buffer + string_offset,"chain %c%3d -> %3d :\t", (node->sequence.flags & STARK_CHAIN_FLAG_REVERSE) ? '^' : ' ', node->sequence.start, node->sequence.end);
	#endif
	
	int_fast8_t i;
	unsigned char mask = 1 << 7;
	for (i = 0; i < 4; i++) {
		if (node->edges & mask) {
			char c = 'N';
			switch(i) {
				case 0: c = 'A'; break;
				case 1: c = 'C'; break;
				case 2: c = 'G'; break;
				case 3: c = 'T'; break;
			}
			string_offset += sprintf(buffer + string_offset,"f%c.%" PRI_OFFSET " ", c, node->child[0][i]);
		}
		mask >>= 1;
	}
	for (i = 0; i < 4; i++) {
		if (node->edges & mask) {
			char c = 'N';
			switch(i) {
				case 0: c = 'A'; break;
				case 1: c = 'C'; break;
				case 2: c = 'G'; break;
				case 3: c = 'T'; break;
			}
			string_offset += sprintf(buffer + string_offset,"r%c.%" PRI_OFFSET " ", c, node->child[1][i]);
		}
		mask >>= 1;
	}
	return offset;
	*/
}

// int stark_print_debug(stark_t* const stark, int maxdepth, depth_t depth) {
// 	printf("depth = %d, maxdepth = %d, stark->size[depth] = %d, stark->deleted[depth] = %d\n",depth, maxdepth, stark->size[depth], stark->deleted[depth]);
// 	return 1;
// }

void debug_print_node(stark_t* const stark, depth_t depth, offset_t offset) {
	char buffer[1024];
	stark_print_node(buffer, stark, depth, offset);
	puts(buffer);
}

void stark_print_depth(stark_t* const stark, depth_t depth) {
	printf("Level %d (size = %ld, deleted = %ld):\n",depth, stark->size[depth], stark->deleted[depth]);
	char buffer[1024];
	offset_t offset;
	for (offset = 1; offset <= stark->size[depth]; offset++) {
		if (stark->level[depth][offset].depth) {
			stark_print_node(buffer, stark, depth, offset);
			printf("N%d : %s \n",offset,buffer);
		}
		else
			printf("N%d : marked deleted\n", offset);
		
	}
}

void stark_print(stark_t* const stark, depth_t maxdepth) {
	printf("Root:\n");
	char buffer[1024];
	stark_print_node(buffer, stark, 0, 0);
	printf("%s\n",buffer);
	
	depth_t depth;
	// for (depth = 1; stark->size[depth] && depth <= maxdepth; depth++) {
	for (depth = 1; stark->size[depth] && (stark->size[depth] - stark->deleted[depth] > 0) && depth <= maxdepth; depth++) {
		stark_print_depth(stark, depth);
	}
	
}


struct coverage_histogram_s* stark_statistics(stark_t* const stark, depth_t maxdepth) {
	
	depth_t depth;
	coverage_t coverage;
	
	struct coverage_histogram_s* coverage_histogram = calloc(1, sizeof(struct coverage_histogram_s) + (1024 * sizeof(int32_t)));
	uint32_t coverage_histogram_memory = 1024;
	coverage_t max_coverage = 0;
	
	// for (depth = 1; stark->size[depth] && depth <= maxdepth; depth++) {
		
	//int i;
	offset_t offset;
	for (depth = 1; stark->size[depth] && (stark->size[depth] - stark->deleted[depth] > 0) && depth <= maxdepth; depth++) {
		for (offset = 1; offset <= stark->size[depth]; offset++) {
			if (stark->level[depth][offset].depth) {
				coverage = stark->level[depth][offset].coverage;
				
				if (coverage >= coverage_histogram_memory) {
					int new_memory_size = coverage_histogram_memory;
					while (coverage >= new_memory_size)
						new_memory_size *= 2;
					
					void* tmp = realloc(coverage_histogram, sizeof(struct coverage_histogram_s) + (new_memory_size * sizeof(int32_t)));
					if (tmp) {
						memset(((char*)tmp) + sizeof(struct coverage_histogram_s) + (coverage_histogram_memory * sizeof(int32_t)), 0, (new_memory_size - coverage_histogram_memory) * sizeof(int32_t));
						coverage_histogram_memory = new_memory_size;
						coverage_histogram = tmp;
					}
					else
					{
						DEBUG_MSG("CRITICAL ERROR");
						perror("realloc");
						exit(-1);
					}
				}
				
				coverage_histogram->hist[coverage]++;
				
				if (max_coverage < coverage)
					max_coverage = coverage;
				
			}
			
		}
	}
	
	// DEBUG_MSG("stats table populated.");
	
	coverage_histogram->max_coverage = max_coverage++;
	
	void* tmp = realloc(coverage_histogram, sizeof(struct coverage_histogram_s) + (max_coverage * sizeof(int32_t)));
	if (tmp)
		coverage_histogram = tmp;
		
	// DEBUG_MSG("realloc completed");
	
	uintmax_t sum_coverage = 0;
	OFFSET_T nodes_counted = 0;
	offset_t count;
	for (coverage = 0; coverage < max_coverage; coverage++) {
		count = coverage_histogram->hist[coverage];
		sum_coverage += count * coverage;
		nodes_counted += count;
	}
	
	
	coverage_histogram->sum_coverage = sum_coverage;
	coverage_histogram->nodes_counted = nodes_counted;
	coverage_histogram->avg_coverage = sum_coverage / (nodes_counted ? nodes_counted : 1);
	
	if (stark->coverage_histogram)
		free(stark->coverage_histogram);
		
	stark->coverage_histogram = coverage_histogram;
	
	return coverage_histogram;
	
}

static inline void stark_compute_coverage_hooks(struct coverage_log_histogram_s * statistics) {
	int log_mass_hook = 0;
	intmax_t current_value = 0;
	
	size_t i;
	for (i = 1; i < (sizeof(statistics->loghist) / sizeof(statistics->loghist[0])); i++) {
		intmax_t new_val = fast_log2(statistics->loghist[i]) * i;
		if (new_val > current_value) {
			log_mass_hook = i;
			current_value = new_val;
		}
	}
	
	statistics->log_mass_hook = log_mass_hook;
}

/*
static inline void stark_statistics_log_depth(struct coverage_log_histogram_s* target, stark_t* const stark, depth_t depth) {
	
	coverage_t coverage, max_coverage = 0;
	offset_t offset;
	uintmax_t sum_coverage = 0;
	
	for (offset = 1; offset <= stark->size[depth]; offset++) {
		if (stark->level[depth][offset].depth) {
			coverage = stark->level[depth][offset].coverage;
			sum_coverage += coverage;
			
			target->loghist[fast_log2(coverage)]++;
			
			if (max_coverage < coverage)
				max_coverage = coverage;
			
		}
		
	}
	
	target->max_coverage = max_coverage;
	
	OFFSET_T nodes_counted = 0;
	int_fast16_t maxlog = fast_log2(max_coverage) + 1;
	for (coverage = 0; coverage < maxlog; coverage++) {
		nodes_counted += target->loghist[coverage];
	}
	
	
	target->sum_coverage = sum_coverage;
	target->nodes_counted = nodes_counted;
	target->avg_coverage = sum_coverage / (nodes_counted ? nodes_counted : 1);
	
	coverage_t sweep, estimated_mean = 0;
	uintmax_t cov_top = 0;
	for (sweep = maxlog - 1; sweep > 1; sweep--) {
		if (target->loghist[sweep] > cov_top) {
			estimated_mean = sweep;
		} else if (target->loghist[sweep] < cov_top && target->loghist[sweep-1] <= cov_top) {
			break;
		}
	}
	
	target->estimated_sample_mean = estimated_mean;
	
	stark_compute_coverage_hooks(target);
}

struct coverage_log_statistics_s* stark_statistics_log(stark_t* const stark) {
	
	depth_t depth, maxdepth;
	
	for (maxdepth = 1; stark->size[maxdepth] && stark->size[maxdepth] - stark->deleted[maxdepth]; maxdepth++);
	
	struct coverage_log_statistics_s* statistics = malloc(sizeof(statistics) + maxdepth * sizeof(statistics->hist[0]));
	memset(statistics, 0, sizeof(statistics) + maxdepth * sizeof(statistics->hist[0]));
	statistics->maxdepth = maxdepth - 1;
	//struct coverage_log_histogram_s* statistics_a = calloc(maxdepth, sizeof(struct coverage_log_histogram_s));
	
	for (depth = 1; depth < maxdepth; depth++) {
		stark_statistics_log_depth(statistics->hist + depth, stark, depth);
	}
	
	return statistics;
}
*/

int stark_statistics_log_loopcb(struct starknode_navigator_s * nav, struct coverage_log_statistics_s * statistics) {
	struct coverage_log_histogram_s* target = statistics->hist + nav->depth;
	
	coverage_t cov = nav->node->coverage;
	target->sum_coverage += cov;
	target->nodes_counted++;
	
	if (target->max_coverage < cov) {
		target->max_coverage = cov;
	}
	
	int covlog = fast_log2(cov);
	
	target->loghist[covlog]++;
	
	return 0;
}

struct coverage_log_statistics_s* stark_statistics_log(stark_t* const stark) {
	
	depth_t depth, maxdepth;
	
	for (maxdepth = 1; stark->size[maxdepth] && stark->size[maxdepth] - stark->deleted[maxdepth] && maxdepth < HARD_MAXDEPTH; maxdepth++);
	
	struct coverage_log_statistics_s* statistics = calloc(1, sizeof(statistics) + maxdepth * sizeof(statistics->hist[0]));
	// memset(statistics, 0, sizeof(statistics) + maxdepth * sizeof(statistics->hist[0]));
	statistics->maxdepth = maxdepth - 1;
	//struct coverage_log_histogram_s* statistics_a = calloc(maxdepth, sizeof(struct coverage_log_histogram_s));
	// 
	// for (depth = 1; depth < maxdepth; depth++) {
	// 	stark_statistics_log_depth(statistics->hist + depth, stark, depth);
	// }
	// 
	struct starknode_navigator_s start = {
		  .stark = stark
		, .depth = 2
		, .offset = 1
		, .direction = 1
	};
	
	stark_navigate_loop_forall(&start, 1, (int (*) (struct starknode_navigator_s *, void *))stark_statistics_log_loopcb, statistics);
	
	for (depth = 2; depth < maxdepth; depth++) {
		statistics->nodes_counted += statistics->hist[depth].nodes_counted;
		statistics->sum_coverage += statistics->hist[depth].sum_coverage;
		if (statistics->max_coverage < statistics->hist[depth].max_coverage)
			statistics->max_coverage += statistics->hist[depth].max_coverage;
	}
	
	statistics->avg_coverage = statistics->sum_coverage / statistics->nodes_counted;
	
	return statistics;
}

/*

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

*/



void stark_print_log_statistics(FILE* stats_fp, struct coverage_log_statistics_s* log_statistics) {
	int_fast16_t i;
	fprintf(stats_fp,"k\tnodes\tsum\tavg");
	for (i = 0; i < 8 * sizeof(coverage_t); i++) {
		fprintf(stats_fp,"\t%" PRIdFAST16, i);
	}
	fprintf(stats_fp,"\n");
			
	depth_t depth;
	struct coverage_log_histogram_s* coverage_log_histogram;
	for (depth = 1; depth <= log_statistics->maxdepth; depth++) {
		coverage_log_histogram = log_statistics->hist + depth;
		fprintf(stats_fp,"%d\t%ld\t%ld\t%d", depth, coverage_log_histogram->nodes_counted, coverage_log_histogram->sum_coverage, coverage_log_histogram->avg_coverage);
		for (i = 0; i < 8 * sizeof(coverage_t); i++) {
			fprintf(stats_fp,"\t%d", coverage_log_histogram->loghist[i]);
		}
		fprintf(stats_fp,"\n");
	}
}


#if defined(UBIGRAPHAPI_H) || defined(UBIGRAPITHREADED_H)

void stark_ubigraph_init_real(void* arg) {
	
	//int i;
	depth_t depth;
	char buffer[256];
	for (depth = 1; depth < HARD_MAXDEPTH; depth++) {
		ubigraph_new_vertex_style_w_id(depth,0);
		//DEBUG_MSG("style_id = %d, i = %d",style_id,i);
		snprintf(buffer, sizeof(buffer),"%d",depth+1);
		ubigraph_set_vertex_style_attribute(depth, "color", buffer);
		if (depth % 2)
			ubigraph_set_vertex_style_attribute(depth, "shape", "sphere");
	}

	
	ubigraph_set_edge_style_attribute(0, "strength", "0.5");
	
	ubigraph_new_edge_style_w_id(-1, 0);
	ubigraph_set_edge_style_attribute(-1, "strength", "0.8");
	ubigraph_set_edge_style_attribute(-1, "width", "1.5");
	
	for (depth = 1; depth < HARD_MAXDEPTH; depth++) {
		style_id_t style = depth * 2;
		ubigraph_new_edge_style_w_id(style, -1);
		// sprintf(buffer,"%d",i +1);
		// ubigraph_set_edge_style_attribute(style, "color", buffer);
		
		style--;
		ubigraph_new_edge_style_w_id(style, -1);
		ubigraph_set_edge_style_attribute(style, "oriented", "true");
	}
	
	ubigraph_new_edge_style_w_id(-2, 0);
	ubigraph_set_edge_style_attribute(-2, "width", "2");
	ubigraph_set_edge_style_attribute(-2, "arrow", "true");
	ubigraph_set_edge_style_attribute(-2, "arrow_radius", "2.5");
	ubigraph_set_edge_style_attribute(-2, "strength", "0");
}

void stark_ubigraph_init() {
	stark_ubigraph_clear();
	thread_pool_dispatch(thread_pool, (void* (*)(void*))stark_ubigraph_init_real, NULL, NULL, 0);
}

void stark_ubigraph_rebuild(stark_t* const stark) {
	stark_ubigraph_init();
	
	DEBUG_MSG("Ubigraph initialized");
	
	vertex_id_t idnum;
	// edge_id_t edgeid;
	stark->ubicount = 0;
	// DEBUG_MSG("(idnum = ++(stark->ubicount)) = %d", (idnum = ++(stark->ubicount)));
	// DEBUG_MSG("stark->ubicount = %d", stark->ubicount);
	
	stark_ubigraph_new_vertex_w_id((idnum = ++(stark->ubicount)));
	// idnum = stark_ubigraph_new_vertex();
	// DEBUG_MSG("idnum = %d", idnum);
	stark_ubigraph_set_vertex_attribute(idnum, "label", "root");
	stark_ubigraph_set_vertex_attribute(idnum, "shape", "dodecahedron");
	stark_ubigraph_set_vertex_attribute(idnum, "color", "1");
	
	// DEBUG_MSG("Root node created", "");
	
	stark->root.ubigraph.id = idnum;
	
	
	
	depth_t depth;
	offset_t offset;
	for (depth = 1; stark->size[depth]; depth++) {
		for (offset = 1; offset <= stark->size[depth]; offset++) {
			
			// DEBUG_MSG("Creating node for depth = %d, offset = %d", depth, i);
			
			if (stark->level[depth][offset].depth) {
				stark_ubigraph_new_vertex_w_id((idnum = ++stark->ubicount));
				stark->level[depth][offset].ubigraph.id = idnum;
				stark_ubigraph_change_vertex_style(idnum, depth);
				char sequence[265];
				get_sequence(sequence, stark, depth, offset);
				stark_ubigraph_set_vertex_attribute(idnum, "label", sequence);
				
				if (STARK_NODE_IS_PALINDROME(stark->level[depth][offset]))
					stark_ubigraph_set_vertex_attribute(idnum, "shape", "torus");
				
				stark_ubigraph_new_edge_stark(stark->level[depth - 1][stark->level[depth][offset].uplink[0]].ubigraph.id, idnum, (depth * 2) - 1);
				stark_ubigraph_new_edge_stark(stark->level[depth - 1][stark->level[depth][offset].uplink[1]].ubigraph.id, idnum, (depth * 2) - 1);
				
				stark_ubigraph_new_edge_stark(stark->level[depth - 1][stark->level[depth][offset].uplink[0]].ubigraph.id, stark->level[depth - 1][stark->level[depth][offset].uplink[1]].ubigraph.id, 2 * (depth - 1));
				
				// stark_ubigraph_new_edge_w_id(edgeid = tuple(stark->level[depth - 1][stark->level[depth][i].uplink[0]].ubigraph.id, idnum),stark->level[depth - 1][stark->level[depth][i].uplink[0]].ubigraph.id, idnum);
				// 
				// stark_ubigraph_change_edge_style(edgeid, (depth * 2) - 1);
				// 
				// stark_ubigraph_new_edge_w_id(edgeid = tuple(stark->level[depth - 1][stark->level[depth][i].uplink[1]].ubigraph.id, idnum), stark->level[depth - 1][stark->level[depth][i].uplink[1]].ubigraph.id, idnum);
				// 
				// stark_ubigraph_change_edge_style(edgeid, (depth * 2) - 1);
				
				// stark->level[depth][i].ubigraph.edge[0] = stark_ubigraph_new_edge(stark->level[depth - 1][stark->level[depth][i].uplink[0]].ubigraph.id, idnum);
				// 				stark->level[depth][i].ubigraph.edge[1] = stark_ubigraph_new_edge(stark->level[depth - 1][stark->level[depth][i].uplink[1]].ubigraph.id, idnum);
			}
			
		}
	}
	
	for (depth = 0; depth < HARD_MAXDEPTH; depth++) {
		stark->level_shown[depth] = 1;
	}
	
}

void stark_ubigraph_hide_vertex_style(style_id_t style) {
	stark_ubigraph_set_vertex_style_attribute(style, "visible", "false");
}
#endif

#ifdef ENCHAIN_SEQUENCES
inline size_t stark_new_chain(stark_t* const stark) {
	return list_new_empty(&stark->chains);
}
#endif

void stark_free(stark_t* const stark) {
	depth_t depth;
	for (depth = 1; depth < HARD_MAXDEPTH; depth++) {
		if (stark->level[depth])
			stark_free_starknode_depth_array(stark, depth);
		
		#if defined(THREADED_IMPORT) && defined(STARK_ALLOC_MALLOC)
		pthread_mutex_destroy(stark->mutex + depth);
		#endif
	}
		
			
	// for (i = 0; i < stark->chains.size; i++)
	// 	list_free(stark->chains.list[i].nodes);

	#ifdef ENCHAIN_SEQUENCES
	list_free(stark->chains);
	list_free(stark->chainnodes);
	#endif
}

void stark_init(stark_t* const stark) {
	
	#ifdef THREADED
	{
	thread_pool_t thread_pool = stark->thread_pool;
	#endif
	
	memset(stark, 0, sizeof(stark_t));
	
	#ifdef THREADED
	thread_pool = stark->thread_pool = thread_pool;
	}
	#endif
	
	stark->maxdepth = HARD_MAXDEPTH;
	
	#ifdef ENCHAIN_SEQUENCES
	list_init( &stark->chains );
	list_init( &stark->chainnodes );
	#endif
	
	stark->root.child[0][0] = 1;
	stark->root.child[0][1] = 2;
	stark->root.child[0][2] = 2;
	stark->root.child[0][3] = 1;
	stark->root.child[1][0] = 1;
	stark->root.child[1][1] = 2;
	stark->root.child[1][2] = 2;
	stark->root.child[1][3] = 1;
	
	stark->root.edges = 0xFF;
	
	// memset(stark->level, 0, HARD_MAXDEPTH * sizeof(starknode_t*));
	// memset(stark->size, 0, HARD_MAXDEPTH * sizeof(intptr_t));
	// memset(stark->maxsize, 0, HARD_MAXDEPTH * sizeof(intptr_t));
	// memset(stark->deleted, 0, HARD_MAXDEPTH * sizeof(intptr_t));
	
	stark->level[0] = &(stark->root);
	
	stark_malloc_starknode_depths(stark);
		
	stark->size[1] = 2;
	stark->level[1][1].depth = 1;
	stark->level[1][1].flags |= STARK_FLAG_PARENT2_REVERSE | STARK_PARENT1_LINK_A | STARK_PARENT2_LINK_T;
	stark->level[1][2].depth = 1;
	stark->level[1][2].flags |= STARK_FLAG_PARENT2_REVERSE | STARK_PARENT1_LINK_C | STARK_PARENT2_LINK_G;
	#ifdef DEBUG_SEQUENCES
	*(stark->level[1][1].sequence) = 'A';
	*(stark->level[1][2].sequence) = 'C';
	#endif
	
	#if defined(THREADED_IMPORT) && defined(STARK_ALLOC_MALLOC)
	depth_t depth;
	for (depth = 1; depth < HARD_MAXDEPTH; depth++)
		pthread_mutex_init(stark->mutex + depth, NULL);
	#endif
	
	
	DEBUG_MSG("Stark Initialized.");
	
	#if defined(UBIGRAPHAPI_H) || defined(UBIGRAPITHREADED_H)
	
	stark_ubigraph_rebuild(stark);
	
	#endif
}


/*
#ifdef __INTEL_COMPILER
// icc intrinsic: int _mm_popcnt_u32(unsigned int a);
static inline bool stark_unambiguos_edges(uint8_t edges) {
	return ((_mm_popcnt_u32(edges & 0xF0) <= 1) && (_mm_popcnt_u32(edges & 0x0F) <= 1));
}
// #define stark_unambiguos_edges(edges) ((_mm_popcnt_u32(edges & 0xF0) > 1) || (_mm_popcnt_u32(edges & 0x0F) > 1))
#else
*/

static inline int_fast8_t stark_unambiguos_edges(int_fast8_t edges) {
	uint_fast8_t edgesforward = edges & 0xF0;
	uint_fast8_t edgesreserse = edges & 0x0F;
	
	return !(edgesforward & (edgesforward - 1)) && !(edgesreserse & (edgesreserse - 1));
	
	// if (!(edgesforward & (edgesforward - 1)) && !(edgesreserse & (edgesreserse - 1)))
	// 	return 1;
	// else
	// 	return 0;
	
}
// #endif

// static const char revser_edge_table[16] = {
// 	[0] = -1, [1] = 3, [1 << 1] = 2, [1 << 2] = 1, [1 << 3] = 0
// };

/*
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
*/

static inline coverage_t stark_substract_shared_coverage(stark_t* const stark, offset_t offset, depth_t depth, coverage_t amount, int_fast8_t direction, int_fast8_t link_char) {
	starknode_t * node = stark->level[depth] + offset;
	
	if (!OFFSET_CONTAINS_SHARED_COVBERAGE(node->child[direction][link_char]))
		return amount;
	
	struct starknode_navigator_s nav = {.stark = stark, .depth = depth, .offset = offset, .direction = direction, .node = node};
		
	offset_t * offsets[2];
	
	if(stark_navigate_shared_link_pair(offsets, &nav, link_char)) {
		CRITICAL_ERROR("Found no share link neighbour where one expected");
	}
	
	if (offsets[0][0] != offsets[1][0]) {
		CRITICAL_ERROR("Shared coverage has to be equal");
	}
	
	long int rest = amount + offsets[0][0];
	if (rest < 0) {
		offsets[0][0] += amount;
		offsets[1][0] += amount;
		
		return 0;
	} else {
		// coverage_t rest = amount + offsets[0][0];
		
		offsets[0][0] = 0;
		offsets[1][0] = 0;
		
		return rest;
	}
	
}


void stark_substract_coverage_seq(stark_t* const stark, struct delcache * uplevel, size_t length) {
	ssize_t i;
	
	int more;
	
	do {
		more = 0;
		starknode_t * node;
		i = length-1;
		// struct {
		// 	intmax_t own;
		// 	intmax_t shared;
		// } left, right;
		 /*
		if (uplevel[i].depth > 1) {
			node = stark->level[uplevel[i].depth] + uplevel[i].offset;
			// attempt to resolve shares
			
			right.own = uplevel[i].amount;
			
			int_fast8_t i;
			for (i = 0; right.own && i < 4; i++) {
				right.own = stark_substract_shared_coverage(stark, uplevel[i].offset, uplevel[i].depth, right.own, uplevel[i].reverse, i);
			}
			right.shared = uplevel[i].amount - right.own;
		}
		*/
		for (; i >= 0; i--) {
			depth_t depth;
			if ((depth = uplevel[i].depth) > 1) {
				more = 1;
				node = stark->level[uplevel[i].depth] + uplevel[i].offset;
				
				flags_t flags = node->flags;
				
				// if (i) {
				// 	
				// }
				
				
				// char buffer[1024];
				// stark_print_node(buffer
				// 	+
				// 	sprintf(buffer, "%c [%c%d, %c%d] "
				// 		, uplevel[i].reverse ? '^' : ' '
				// 		, parent_Link_to_char[uplevel[i].backchar]
				// 		, node->child[uplevel[i].reverse ^ 1][uplevel[i].backchar]
				// 		, parent_Link_to_char[uplevel[i].forwardchar]
				// 		, node->child[uplevel[i].reverse][uplevel[i].forwardchar]
				// 		)
				// 	, stark, uplevel[i].depth, uplevel[i].offset);
				// puts(buffer);
				
				coverage_t old_cov = node->coverage;
				node->coverage -= uplevel[i].amount;
				
				//normalize excessive coverage removal
				// this solution is not ideal, but temporary
				{
					size_t j;
					offset_t child_offset;
					coverage_t min_forward = 0;
					for (j = 0; j < 4; j++)
						if (OFFSET_VALID(child_offset = node->child[0][j]))
							min_forward += stark->level[depth+1][child_offset].coverage;
				
					coverage_t min_back = 0;
					for (j = 0; j < 4; j++)
						if (OFFSET_VALID(child_offset = node->child[1][j]))
							min_back += stark->level[depth+1][child_offset].coverage;
				
					if (node->coverage < min_forward)
						node->coverage = min_forward;
					if (node->coverage < min_back)
						node->coverage = min_back;
					
					if (node->coverage == old_cov) {
						memset(uplevel + i, 0, sizeof(uplevel[i]));
						continue;
					}
				}
				
				

				uplevel[i+1].depth = --uplevel[i].depth;
				if (!uplevel[i].reverse) {
					uplevel[i+1].offset = node->uplink[1];
					uplevel[i+1].reverse = flags & STARK_FLAG_PARENT2_REVERSE ? 1 : 0;
					uplevel[i+1].backchar = flags & STARK_PARENT2_LINK_MASK;
					uplevel[i].offset = node->uplink[0];
					uplevel[i].reverse = flags & STARK_FLAG_PARENT1_REVERSE ? 1 : 0;
					uplevel[i].forwardchar = (flags & STARK_PARENT1_LINK_MASK) >> 4;
				} else {
					uplevel[i].offset = node->uplink[1];
					uplevel[i].reverse = flags & STARK_FLAG_PARENT2_REVERSE ? 0 : 1;
					uplevel[i].forwardchar = flags & STARK_PARENT2_LINK_MASK;
					uplevel[i+1].offset = node->uplink[0];
					uplevel[i+1].reverse = flags & STARK_FLAG_PARENT1_REVERSE ? 0 : 1;
					uplevel[i+1].backchar = (flags & STARK_PARENT1_LINK_MASK) >> 4;
				}
				
				if (uplevel[i+1].amount < uplevel[i].amount)
					uplevel[i+1].amount = uplevel[i].amount;
				
				// delete the node if it's empty
				// if (!node->coverage && !(flags & STARK_FLAG_NODELETE)) {
				// 	starknode_t* forward_parent = stark->level[depth-1] + node->uplink[0];
				// 	starknode_t* reverse_parent = stark->level[depth-1] + node->uplink[1];
				// 	forward_parent->child[flags & STARK_FLAG_PARENT1_REVERSE ? 1 : 0][(flags & STARK_PARENT1_LINK_MASK) >> 4] = - (uplevel[i].amount);
				// 	reverse_parent->child[flags & STARK_FLAG_PARENT2_REVERSE ? 0 : 1][flags & STARK_PARENT2_LINK_MASK] = - (uplevel[i].amount);
				// 	stark_node_init(node);
				// 	__sync_fetch_and_add(stark->deleted + depth, 1);
				// }
				
				/*
				if (uplevel[i].reverse)
					buffer[0] = '^';
				else
					buffer[0] = ' ';
				stark_print_node(buffer+1, stark, uplevel[i].depth, uplevel[i].offset);
				puts(buffer);
				if (uplevel[i+1].reverse)
					buffer[0] = '^';
				else
					buffer[0] = ' ';
				stark_print_node(buffer+1, stark, uplevel[i+1].depth, uplevel[i+1].offset);
				puts(buffer);
				*/

				if (i == length-1)
					length++;
			}

		}
		// puts("\n");
	} while (more);
	
	
}

void stark_substract_coverage_single(stark_t* const stark, depth_t depth, offset_t offset, coverage_t amount) {
	struct delcache uplevel[depth + 1];
	memset(uplevel, 0 , sizeof(uplevel));
	
	uplevel[0].depth = depth;
	uplevel[0].offset = offset;
	uplevel[0].amount = amount;
	
	stark_substract_coverage_seq(stark, uplevel, 1);
}


void stark_substract_coverage(stark_t* const stark, offset_t offset, depth_t depth, coverage_t amount, int_fast8_t share_forward, int_fast8_t share_back) {
	
	stark_print(stark, HARD_MAXDEPTH);
	fflush(stdout);
	DEBUG_MSG("subbing %d coverage from [%d,%d]", amount, depth, offset);
	
	starknode_t * node = stark->level[depth] + offset;
	
	if (STARK_NODE_IS_PALINDROME(*node)) {
		DEBUG_MSG("Palindromes canot have coverage substracted. Sorry, feature unsupported.");
		return;
	}
	
	if (node->coverage < amount) {
		CRITICAL_ERROR("You cannot substract more coverage then available");
	}
	
	coverage_t sub_forward = amount;
	coverage_t sub_back = amount;
	int_fast8_t i;
	for (i = 0; i < 4; i++) {
		sub_forward = stark_substract_shared_coverage(stark, offset, depth, sub_forward, 0, i);
		sub_back = stark_substract_shared_coverage(stark, offset, depth, sub_back, 1, i);
	}
		
	node->coverage -= amount;
	
	if (!node->coverage) {
		// delete the mode and add shared coverage below
		*(stark_navigate_parent_offset_to_self(stark, depth, offset, 0)) = -amount;
		*(stark_navigate_parent_offset_to_self(stark, depth, offset, 1)) = -amount;
		// stark_node_init(node);
		__sync_fetch_and_add(stark->deleted + depth, 1);
	}
	
	flags_t flags = node->flags;
	if (sub_back > 0) {
		

		int_fast8_t ch = (flags & STARK_PARENT1_LINK_MASK) >> 4;

		if (depth > 2)
		
		stark_substract_coverage(stark, node->uplink[0], depth - 1, sub_back, 
			(flags & STARK_FLAG_PARENT2_REVERSE) ? ch : share_forward, 
			(flags & STARK_FLAG_PARENT2_REVERSE) ? share_forward : ch);
	}
	
	if (sub_forward > 0) {

		int_fast8_t ch = (flags & STARK_PARENT2_LINK_MASK);

		if (depth > 2)
		
		stark_substract_coverage(stark, node->uplink[1], depth - 1, sub_forward, 
			(flags & STARK_FLAG_PARENT1_REVERSE) ? share_back : ch, 
			(flags & STARK_FLAG_PARENT1_REVERSE) ? ch : share_back );
	}
	
}



size_t stark_unambiguos_clean(stark_t* const stark, depth_t minK) {
	depth_t depth;
	offset_t offset;
	size_t erased = 0;
	for (depth = HARD_MAXDEPTH -1; depth > minK; depth--) {
		
		size_t deleted = 0;
		#pragma omp for schedule(guided,0x1000)
		for (offset = 1; offset <= stark->size[depth]; offset++) {
			starknode_t* __restrict current = stark->level[depth] + offset;
			starknode_t* forward_parent = stark->level[depth-1] + current->uplink[0];
			starknode_t* reverse_parent = stark->level[depth-1] + current->uplink[1];
			flags_t flags = current->flags;
			if ((flags & STARK_FLAG_NODELETE) || STARK_NODE_IS_PALINDROME(*forward_parent) || STARK_NODE_IS_PALINDROME(*reverse_parent)) {
				forward_parent->flags |= STARK_FLAG_NODELETE;
				reverse_parent->flags |= STARK_FLAG_NODELETE;
				if (stark->longest_palindrome < depth-1) {
					stark->longest_palindrome = depth-1; // should be synchronous
				}
					
			} else {
				// int_fast8_t forward_parent_edges = forward_parent->edges;
				// int_fast8_t reverse_parent_edges = reverse_parent->edges;
				if (stark_unambiguos_edges(forward_parent->edges) && stark_unambiguos_edges(reverse_parent->edges)) { // this is a redundant child
					// DEBUG_MSG("Removing Node %d.%" PRI_OFFSET,depth,offset);
					
					
					// NEW IN 0.4.1, put in shared coverage
					forward_parent->child[flags & STARK_FLAG_PARENT1_REVERSE ? 1 : 0][(flags & STARK_PARENT1_LINK_MASK) >> 4] = - (current->coverage);
					reverse_parent->child[flags & STARK_FLAG_PARENT2_REVERSE ? 0 : 1][flags & STARK_PARENT2_LINK_MASK] = - (current->coverage);
					
					// old variant
					// forward_parent->child[flags & STARK_FLAG_PARENT1_REVERSE ? 1 : 0][(flags & STARK_PARENT1_LINK_MASK) >> 4] = 0;
					// reverse_parent->child[flags & STARK_FLAG_PARENT2_REVERSE ? 0 : 1][flags & STARK_PARENT2_LINK_MASK] = 0;
					
					
					#if defined(UBIGRAPHAPI_H) || defined(UBIGRAPITHREADED_H)
					stark_ubigraph_set_vertex_attribute(stark->level[depth][offset].ubigraph.id,"visible","false");
					// stark_ubigraph_set_edge_attribute(stark->level[depth][i].ubigraph.edge[0],"visible","false");
					// stark_ubigraph_set_edge_attribute(stark->level[depth][i].ubigraph.edge[1],"visible","false");
					
					depth_t depth2;
					offset_t j;
					for (depth2 = depth-1; depth2 <= depth +1; depth2++) {
						// printf("Level %d:\n",depth);
						for (j = 1; j <= stark->size[depth2]; j++) {
							starknode_t* node2 = stark->level[depth2] + j;
							
							if (node2->depth){
								stark_ubigraph_set_edge_attribute(tuple(node2->ubigraph.id, stark->level[depth][offset].ubigraph.id), "visible", "false");
							}
							
						}
					}
					#endif
					
					
					// erase node
					stark_node_init(current);
					
					// stark->deleted[depth]++;
					deleted++;
					
					erased++;
				} else {
					forward_parent->flags |= STARK_FLAG_NODELETE;
					reverse_parent->flags |= STARK_FLAG_NODELETE;
				}
			}
		}
		
		__sync_fetch_and_add(stark->deleted + depth, deleted);
	}
	
	// fix broken shared coverage info
	for (depth = HARD_MAXDEPTH -1; depth; depth--) {
		
		#pragma omp for schedule(guided,0x10000)
		for (offset = 1; offset <= stark->size[depth]; offset++) {
			starknode_t* __restrict current = stark->level[depth] + offset;
			
			if (!STARK_NODE_VALID(*current))
				continue;
			
			int_fast8_t i;
			for (i = 0; i < 4; i++) {
				if (OFFSET_CONTAINS_SHARED_COVBERAGE(current->child[0][i])) {
					struct starknode_navigator_s nav = {.stark = stark, .depth = depth, .offset = offset, .direction = 0, .node = current};
					offset_t * offsets[2];
					
					if(stark_navigate_shared_link_pair(offsets, &nav, i)) {
						current->child[0][i] = 0;
					}
				}
				if (OFFSET_CONTAINS_SHARED_COVBERAGE(current->child[1][i])) {
					struct starknode_navigator_s nav = {.stark = stark, .depth = depth, .offset = offset, .direction = 1, .node = current};
					offset_t * offsets[2];
					
					if(stark_navigate_shared_link_pair(offsets, &nav, i)) {
						current->child[1][i] = 0;
					}
				}
			}
		}
	}
	
	
	return erased;
}

void stark_compact_task(stark_t * stark, size_t resize_current, offset_t * repointers_upper, offset_t * revrepointers_current, offset_t * repointers_lower, depth_t depth, int * releasecount) {
	#if defined(STARK_STATUS_H)
	stark_set_status_action("Compacting Task");
	#endif
	
	starknode_t node;
	offset_t newoffset;
	for (newoffset = 1; newoffset <= resize_current; newoffset++) {
		// offset = revrepointers_current[newoffset];
		// repointers_current[offset] = newoffset;
		// revrepointers[depth][resizes[depth]] = offset;
		// printf("%d.%" PRI_OFFSET " = %" PRI_OFFSET "\n",depth, newoffset, revrepointers_current[newoffset]);

		node = stark->level[depth][revrepointers_current[newoffset]];

		node.uplink[0] = repointers_lower[node.uplink[0]];
		node.uplink[1] = repointers_lower[node.uplink[1]];

		size_t i;
		for (i = 0; i < 8; i++) {
			offset_t offset = ((offset_t*)(node.child))[i];
			if (OFFSET_VALID(offset))
				((offset_t*)(node.child))[i] = repointers_upper[offset];
		}

		stark->level[depth][newoffset] = node;

	}
	
	

	stark->size[depth] = resize_current;
	stark->deleted[depth] = 0;

	if (resize_current) {
		stark_shrink_starknode_depth_array(stark, depth);
		stark_move_starknode_depth_array_to_numa_master(stark, depth);
	}
	else 
		stark_free_starknode_depth_array(stark, depth);
	
	// while (depth < HARD_MAXDEPTH -2 && *((((volatile size_t *)(stark->deleted))) + depth + 2));
	
	if (__sync_fetch_and_add(releasecount +depth-2, -1) == 0) {
		free(repointers_lower);
	}

	if (__sync_fetch_and_add(releasecount +depth, -1) == 0) {
		free(repointers_upper);
	}
	
	free(revrepointers_current);
	
	#if defined(STARK_STATUS_H)
	stark_set_status_action("Sleeping");
	#endif
}

void stark_compact(stark_t* stark) {
	
	
	depth_t depth, maxdepth = 1;
	offset_t offset;
	// starknode_t node;
	offset_t* repointers_upper, *repointers_current, *repointers_lower;
	offset_t* revrepointers_current, *revrepointers_lower;
	size_t resize_current, resize_lower;
	
	for (depth = 1; stark->level[depth] && depth < HARD_MAXDEPTH; depth++) {
		if ((stark->size[depth] - stark->deleted[depth]))
			maxdepth = depth;
		else 
			stark_shrink_starknode_depth_array(stark, depth);
	}
	depth = maxdepth;
	

	repointers_upper = calloc(1, sizeof(offset_t));
	
	
	resize_current = 0;
	repointers_current = calloc(stark->size[depth] + 1, sizeof(offset_t));
	revrepointers_current = calloc(stark->size[depth] - stark->deleted[depth] + 1, sizeof(offset_t));
	for (offset = 1; offset <= stark->size[depth]; offset++) {
		if (!stark->level[depth][offset].depth) 
			continue;
		
		repointers_current[offset] = ++resize_current;
		revrepointers_current[resize_current] = offset;
		
	}
	
	int releasecount[depth];
		
	for (; depth > 1; depth--) {
		
		releasecount[depth-2] = 1;
		
		{
			depth--;
			repointers_lower = calloc(stark->size[depth] + 1, sizeof(offset_t));
			revrepointers_lower = calloc(stark->size[depth] - stark->deleted[depth] + 1, sizeof(offset_t));
		
			resize_lower = 0;
			for (offset = 1; offset <= stark->size[depth]; offset++) {
				if (stark->level[depth][offset].depth) {
					repointers_lower[offset] = ++resize_lower;
					revrepointers_lower[resize_lower] = offset;
				}
			
			}
			depth++;
		}
		
		#pragma omp task untied
		stark_compact_task(stark, resize_current, repointers_upper, revrepointers_current, repointers_lower, depth, releasecount);
		
		repointers_upper = repointers_current;
		repointers_current = repointers_lower;
		revrepointers_current = revrepointers_lower;
		resize_current = resize_lower;
	}
	
	#pragma omp taskwait
	
	free(repointers_upper);
	free(repointers_current);
	free(revrepointers_current);
	

}


/*
#include <md5.h>

inline void md5(const unsigned char* data, size_t length, unsigned char* digest) {
	md5_state_t pms;
	
	md5_init(&pms);
	md5_append(&pms, data, length);
	
	md5_finish(&pms, digest);
}

int md5_check(const char* data, size_t length, char* hash) {
	char digest[16];
	
	md5(data, length, digest);
	
	int i;
	for (i = 0; i < (16 / sizeof(intptr_t)); i++) {
		if(((intptr_t*)hash)[i] != ((intptr_t*)digest)[i])
			return 1;
	}
	return 0;
}
*/


int stark_find(depth_t* depth_r, offset_t* offset_r, stark_t* const stark, char* sequence) {
	
	depth_t depth = 0;
	offset_t offset = 0;
	char direction = 0;
	
	starknode_t* node = stark->level[depth] + offset;
	
	char c;
	
	while ((c = sequence[depth])) {
		if (OFFSET_VALID(offset = node->child[direction][foursymbol_index[c]])) {
			depth++;
			node = stark->level[depth] + offset;
			if (reversing(sequence, depth) == -1)
				direction = 1;
			else
				direction = 0;
		} else
			return 0;
	}
	
	*depth_r = depth;
	*offset_r = offset;
	
	return 1;
}

// void stark_swap_nodes(stark_t * stark, depth_t depth, offset_t offsets[2]) {
// 	starknode_t temp = 
// }

#include "getusage.h"

size_t stark_memory_usage(stark_t* const stark) {
	
	
	struct pstat result;
	if(!get_usage(&result)) {
		return result.rss;
	}
	
	
	size_t memory = sizeof(stark_t);
	
	depth_t depth;
	
	for (depth = 1; stark->level[depth]; depth++) {
		memory += sizeof(starknode_t) * stark->maxsize[depth];
	}
	
	fprintf(stderr, "depth = %d\n", depth);
	
	return memory;
}




