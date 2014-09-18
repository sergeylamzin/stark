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


//Last modified: 03/03/12 00:44:18(CET) by Fabian Holler
// Source: https://raw.github.com/fho/code_snippets/master/c/getusage.c
#include <stdlib.h> 
#include <sys/types.h>
#include <stdio.h>
#include <strings.h>
#include <string.h>
#include <unistd.h>
#include "getusage.h"

/*
 * read /proc data into the passed struct pstat
 * returns 0 on success, -1 on error
*/
int get_usage(struct pstat* result) {
	char stat_filepath[30];

	snprintf(stat_filepath, sizeof(stat_filepath), "/proc/%d/stat", getpid());

    //Open /proc/stat and /proc/$pid/stat fds successive(dont want that cpu
    //ticks increases too much during measurements)
    //TODO: open /proc dir, to lock all files and read the results from the
    //same timefragem
    FILE *fpstat = fopen(stat_filepath, "r");
    if(fpstat == NULL){
        // perror("FOPEN ERROR ");
        return -1;
    }

    FILE *fstat = fopen("/proc/stat", "r");
    if(fstat == NULL){
        // perror("FOPEN ERROR ");
        fclose(fstat);
        return -1;
    }

    //read values from /proc/pid/stat
    memset(result, 0, sizeof(struct pstat));
    long int rss;
    if(fscanf(fpstat, "%*d %*s %*c %*d %*d %*d %*d %*d %*u %*u %*u %*u %*u %lu"
                "%lu %ld %ld %*d %*d %*d %*d %*u %lu %ld",
                &result->utime_ticks, &result->stime_ticks,
                &result->cutime_ticks, &result->cstime_ticks, &result->vsize,
                &rss) == EOF) {
        fclose(fpstat);
        return -1;
    }
    fclose(fpstat);
    result->rss = rss * sysconf(_SC_PAGE_SIZE);

    //read+calc cpu total time from /proc/stat
    long unsigned int cpu_time[10];
    bzero(cpu_time, sizeof(cpu_time));
    if(fscanf(fstat, "%*s %lu %lu %lu %lu %lu %lu %lu %lu %lu %lu",
                &cpu_time[0], &cpu_time[1], &cpu_time[2], &cpu_time[3],
                &cpu_time[4], &cpu_time[5], &cpu_time[6], &cpu_time[7],
                &cpu_time[8], &cpu_time[9]) == EOF) {
        fclose(fstat);
        return -1;
    }

    fclose(fstat);

	int i;
    for(i=0; i < 10;i++){
        result->cpu_total_time += cpu_time[i];
    }

    return 0;
}

/*
* calculates the actual CPU usage(cur_usage - last_usage) in percent
* cur_usage, last_usage: both last measured get_usage() results
* ucpu_usage, scpu_usage: result parameters: user and sys cpu usage in %
*/
void calc_cpu_usage(const struct pstat* cur_usage, const struct pstat*
                    last_usage, double* ucpu_usage, double* scpu_usage){
    const long unsigned int total_time_diff = cur_usage->cpu_total_time -
                                              last_usage->cpu_total_time;

    *ucpu_usage = 100 * (((cur_usage->utime_ticks + cur_usage->cutime_ticks)
                    - (last_usage->utime_ticks + last_usage->cutime_ticks))
                    / (double) total_time_diff);

    *scpu_usage = 100 * ((((cur_usage->stime_ticks + cur_usage->cstime_ticks)
                    - (last_usage->stime_ticks + last_usage->cstime_ticks))) /
                    (double) total_time_diff);
}
