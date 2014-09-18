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
#include "stark_navigate.h"

struct stark_assembler_segment_s {
	ssize_t id;
	
	size_t length;
	
	struct {
		offset_t offset;
		depth_t depth;
	} start, end;
	
	struct {
		coverage_t max;
		coverage_t min;
		size_t sum;
		size_t sum_square;
		double factor;
	} coverage;
	
};

struct stark_assembler_s {
	stark_t * stark;
	depth_t minK;
	struct coverage_log_statistics_s* log_statistics;
	
	size_t garbage_length;
	size_t long_length;
	
	depth_t maxdepth;
	offset_t nextoffset;
	
	FILE * longout, *shortout;
	
	double link_strength;
	
	struct {
		size_t size;
	} contigs;
		
	struct {
		struct starknode_navigator_s * list;
		size_t size;
		size_t num_pages;
	} surface[64];
	
	struct starknode_navigator_s candidates[65];
};


void stark_assembler_test(stark_t * stark, struct coverage_log_statistics_s* log_statistics, double link_strength);

