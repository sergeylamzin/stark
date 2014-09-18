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

#ifndef STARK_R_H_IKKODAS
#define STARK_R_H_IKKODAS

void init_R();
void stark_generate_coverage_contour_plot_R(struct coverage_log_statistics_s* log_statistics, char * Rbasefilename);

#endif /* end of include guard: STARK_R_H_IKKODAS */
