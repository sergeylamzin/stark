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
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>

#ifndef MAIN_H
#define MAIN_H

#if defined(__INTEL_COMPILER)
//
// Disable ICC's remark #869: "Parameter" was never referenced warning.
// This is legal ANSI C code so we disable the remark that is turned on with -Wall
//
// #pragma warning ( disable : 869 )

//
// Disable ICC's remark #1418: external function definition with no prior declaration.
// This is legal ANSI C code so we disable the remark that is turned on with /W4
//
#pragma warning ( disable : 1418 )

//
// Disable ICC's remark #1419: external declaration in primary source file
// This is legal ANSI C code so we disable the remark that is turned on with /W4
//
#pragma warning ( disable : 1419 )

//
// Disable ICC's remark #2259: non-pointer conversion from "int" to "... may lose significant bits
// This is legal ANSI C code so we disable the remark that is turned on with /W4
//
#pragma warning ( disable : 2259 )

//
// Disable ICC's remark #593: "Variable" was set but never used.
// This is legal ANSI C code so we disable the remark that is turned on with /W4
//
//#pragma warning ( disable : 593 )

#endif

#define CONCATENATE_2(Older, Small) Older##Small
#define CONCATENATE(Older, Small) CONCATENATE_2(Older, Small)

typedef int bool;
#define true 1
#define false 0

#define rdtsc64() ({uint64_t eax,edx; __asm__ ("rdtsc": "=a" (eax), "=d" (edx)); (eax | (edx << 32));})
#define rdtsc32() ({uint32_t eax; __asm__ ("rdtsc": "=a" (eax) : : "edx"); eax;})

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#define CRITICAL_ERROR(...) (\
	fprintf(stderr, __FILE__ ": Critical Error!\n" \
	 				__FILE__ ":" TOSTRING(__LINE__) ": " __VA_ARGS__), \
	fprintf(stderr, "\n"), \
	DEBUG_PRINT_STACK_TRACE(), \
	exit(-1) \
) 

#ifdef DEBUG
#define DEBUG_MSG(...) (\
	fprintf(stderr, __FILE__ ":" TOSTRING(__LINE__) ": " __VA_ARGS__), \
	fprintf(stderr, "\n"), \
	fflush(stderr) \
)

#include <execinfo.h>

static inline void DEBUG_PRINT_STACK_TRACE(void) {
	void *array[10];
    size_t size;

    size = backtrace (array, 10);
    backtrace_symbols_fd (array, size, 2);
}

#else
#define DEBUG_MSG(...) 0
#define DEBUG_PRINT_STACK_TRACE() 0
#endif

// #define DEBUG_MSG(format, ...) fprintf(stderr, __FILE__ ":%d: " format "\n", __LINE__, __VA_ARGS__)


// #define min(a, b) ((a < b) ? a : b)

#if defined(__INTEL_COMPILER)
#define fast_log2(ul) _bit_scan_reverse(ul)
#elif defined(__x86_64__)

static inline int fast_log2(uint64_t const value) {
	uint64_t dst;
	__asm__ (
		"bsrq %1, %0\n\
		cmove %2, %0"
		: "=r" (dst)
		: "mr" (value) , "mr" (-1L)
		);
	return dst;
}

#elif defined(__GNUC__)
#define fast_log2(ull) (63 - __builtin_clzll(ull))
#else
static inline char fast_log2(uint64_t n) {
	
	#define INT1_MASK  0x1
	#define INT2_MASK  0x3
	#define INT4_MASK  0xF
	#define INT8_MASK  0xFF
	#define INT16_MASK 0xFFFF
	#define INT32_MASK 0xFFFFFFFF
	#define INT64_MASK 0xFFFFFFFFFFFFFFFF
	
	uint_fast16_t l = 0;
	
	if (n & (INT32_MASK << 32)) {
		l += 32;
		n >>= 32;
	}
	
	if (n & (INT16_MASK << 16)) {
		l += 16;
		n >>= 16;
	}
	
	if (n & (INT16_MASK << 16)) {
		l += 16;
		n >>= 16;
	}
	
	if (n & (INT8_MASK << 8)) {
		l += 8;
		n >>= 8;
	}
	
	if (n & (INT8_MASK << 8)) {
		l += 8;
		n >>= 8;
	}
	
	if (n & (INT4_MASK << 4)) {
		l += 4;
		n >>= 4;
	}
	
	if (n & (INT2_MASK << 2)) {
		l += 2;
		n >>= 2;
	}
	
	if (n & (INT1_MASK << 1)) {
		l += 1;
		// n >>= 1;
	}
	
	// while (n >>= 1)
	// 	l++;
	
	return l;
}
#endif


#endif
