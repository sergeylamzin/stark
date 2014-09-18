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

#include <unistd.h>
#include <fcntl.h>
#include <termios.h>
#include "main.h"
#include "stark.h"
#include "stark_assemble_hierarchical.h"

#ifndef TERMINAL_H
#define TERMINAL_H


struct terminal_session {
	stark_t* stark;
	FILE* example_fp;
	char* example_next_line;
	struct stark_hierarchical_assembler_s *hierarchical_assembler;
};

// int terminal_setup(int fd);
int terminal_get_command(char* buffer, int maxlength, FILE* fp);
int process_command_line(struct terminal_session* session, char* input);


#endif 

