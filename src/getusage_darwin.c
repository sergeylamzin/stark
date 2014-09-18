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


#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/sysctl.h>
#include <mach/task.h>
#include <mach/mach_init.h>

// void getres(task_t task, unsigned int *rss, unsigned int *vs)
// {
//     struct task_basic_info t_info;
//     mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
// 
//     task_info(task, TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count);
//     *rss = t_info.resident_size;
//     *vs = t_info.virtual_size;
// }
// 
// void memcheck()
// {
//     unsigned int rss, vs, psize;
//     task_t task = MACH_PORT_NULL;
// 
//     if (task_for_pid(current_task(), getpid(), &task) != KERN_SUCCESS)
//         abort();
//     getres(task, &rss, &vs);
//     psize = getpagesize();
//     fprintf(stderr, "RSS: %u KiB, VS: %u KiB.\n", rss, vs);
// }

#include "getusage.h"


int get_usage(struct pstat* result) {
	struct task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

    task_info(mach_task_self(), TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count);
    result->rss = t_info.resident_size;
    result->vsize = t_info.virtual_size;

	// result->utime_ticks = t_info.user_time;
	
	// result->stime_ticks = t_info.system_time;
	
	return 0;
}