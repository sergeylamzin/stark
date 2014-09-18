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


#include "stark.h"

#ifdef UBIGRAPHAPI_H
#define UBISTRING ":UBI"
#else
#define UBISTRING
#endif

const char starkinit[] = "STARK:O" TOSTRING(OFFSET_L) ":HARD_MAXDEPTH:" TOSTRING(HARD_MAXDEPTH) UBISTRING ":MD5:" TOSTRING(STARK_MD5) "\n";

const char foursymbol_index[] = {['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3, ['U'] = 3};
//char foursymbol_index[2]['U'+1] = {{['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3, ['U'] = 3},{['A'] = 3, ['C'] = 2, ['G'] = 1, ['T'] = 0, ['U'] = 0}};

const char parent_Link_index[] = {['A'] = STARK_PARENT2_LINK_A, ['C'] = STARK_PARENT2_LINK_C, ['G'] = STARK_PARENT2_LINK_G, ['T'] = STARK_PARENT2_LINK_T, ['U'] = STARK_PARENT2_LINK_T};

const char parent_Link_to_char[] = {
	[STARK_PARENT1_LINK_A] = 'A',
	[STARK_PARENT1_LINK_C] = 'C',
	[STARK_PARENT1_LINK_G] = 'G',
	[STARK_PARENT1_LINK_T] = 'T',
	[STARK_PARENT2_LINK_A] = 'A',
	[STARK_PARENT2_LINK_C] = 'C',
	[STARK_PARENT2_LINK_G] = 'G',
	[STARK_PARENT2_LINK_T] = 'T'};

const char complement_index[] = {['A'] = 'T', ['C'] = 'G', ['G'] = 'C', ['T'] = 'A', ['U'] = 'A'};

const unsigned char edge_index[2]['U'+1] = {
	{
		[0] = (1<<(3+4))
		, [1] = (1<<(2+4))
		, [2] = (1<<(1+4))
		, [3] = (1<<4)
		,['A'] = (1<<(3+4))
		, ['C'] = (1<<(2+4))
		, ['G'] = (1<<(1+4))
		, ['T'] = (1<<4)
		, ['U'] = (1<<4)
	}
	, {
		[0] = (1<<3)
		, [1] = (1<<2)
		, [2] = (1<<1)
		, [3] = (1)
		,['A'] = (1<<3)
		, ['C'] = (1<<2)
		, ['G'] = (1<<1)
		, ['T'] = (1)
		, ['U'] = (1)
	}
};


const unsigned char edgeid_single_index[] = {
	[1 << 7] = 0,
	[1 << 6] = 1,
	[1 << 5] = 2,
	[1 << 4] = 3,
	[1 << 3] = 0,
	[1 << 2] = 1,
	[1 << 1] = 2,
	[1] = 3
};

const unsigned char edgeid_single_char[] = {
	[1 << 7] = 'A',
	[1 << 6] = 'C',
	[1 << 5] = 'G',
	[1 << 4] = 'T',
	[1 << 3] = 'A',
	[1 << 2] = 'C',
	[1 << 1] = 'G',
	[1] = 'T'
};

const unsigned char nucleotides[] = {'A', 'C', 'G', 'T'};

