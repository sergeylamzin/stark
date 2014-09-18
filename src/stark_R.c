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
#include "stark_contour_log_stats_plot.h"
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

static inline int which_query(const char * query, const char * sub_path, int len, char * outpath, size_t outpath_len) {
	struct stat s;
	char filename[0x1000];
	size_t length = snprintf(filename, sizeof(filename), "%.*s/%s", len, sub_path, query);
	if (!stat(filename, &s)) {
		if (outpath)
			snprintf(outpath, outpath_len, "%s", filename);
		return 0;
	}
	return -1;
}

int which(const char * query, const char * path, const char * const * path_a, char * outpath, size_t outpath_len) {
	if (path) {
		const char * sub_path = path;
		const char * end;
		while ((end = strchr(sub_path, ':'))) {
			if (!which_query(query, sub_path, end - sub_path, outpath, outpath_len))
				return 0;
			sub_path = end + 1;
		}
		if (!which_query(query, sub_path, strlen(sub_path), outpath, outpath_len))
			return 0;
	}
	while (path_a) {
		if (!which_query(query, *path_a, strlen(*path_a), outpath, outpath_len))
			return 0;
		path_a++;
	}
	return -1;
}

static int useR = 1;

void init_R() {

	useR = which("Rscript", getenv("PATH"), NULL, NULL, 0) != 0;

	// pid_t pid = fork();
	// if (pid) {
	// 	int status = 1;
		
	// 	if (pid > 0)
	// 		waitpid(pid, &status, 0);
		
	// 	if (status) {
	// 		fprintf(stderr,"Rscript not found, will generate an .r file instead.\n");
	// 		useR = 1;
	// 	} else {
	// 		fprintf(stderr,"Rscript executable found.\n");
	// 		useR = 0;
	// 	}
	// } else {
	// 	dup2(2,1);
	// 	char* argv1[] = {"which", "Rscript", 0};
	// 	execvp(argv1[0],argv1);
	// 	perror("execvp");
	// 	exit(-1);
	// }
}


void stark_generate_coverage_contour_plot_R(struct coverage_log_statistics_s* log_statistics, char * Rbasefilename) {
	FILE* rfp;
	if (useR) {
		char fname[strlen(Rbasefilename)+3];
		sprintf(fname,"%s.r",Rbasefilename);
		rfp = fopen(fname,"w");
	} else {
		int pd[2];
		if(pipe(pd)){
			perror("pipe");
			return;
		}
		
		rfp = fdopen(pd[1],"w");
		
		pid_t child_pid = fork();
		
		if (child_pid < 0) {
			perror("fork");
			return;
		}
		
		if (!child_pid) {
			close(pd[1]);
			dup2(pd[0],0);
			dup2(2,1);
			char* argv1[] = {"Rscript", "--silent", "-", 0};
			execvp(argv1[0],argv1);
			perror("execvp");
			exit(-1);
		}
		close(pd[0]);
	}
	
	fprintf(rfp, "data <- scan()\n");
	
	size_t maxcov = 0;
	depth_t depth;
	size_t i;
	struct coverage_log_histogram_s* coverage_log_histogram;
	for (depth = 2; depth <= log_statistics->maxdepth; depth++) {
		coverage_log_histogram = log_statistics->hist + depth;
		for (i = 0; i < 8 * sizeof(coverage_t); i++) {
			if (maxcov < i && coverage_log_histogram->loghist[i])
				maxcov = i;
		}
	}
	for (depth = 2; depth <= log_statistics->maxdepth; depth++) {
		coverage_log_histogram = log_statistics->hist + depth;
		for (i = 0; i < maxcov; i++) {
			fprintf(rfp,"%d ", coverage_log_histogram->loghist[i] ? coverage_log_histogram->loghist[i] : 1);
		}
		fprintf(rfp,"\n");
	}
	fprintf(rfp,"\n\nx <- scan()\n");
	for (depth = 2; depth <= log_statistics->maxdepth; depth++) {
		fprintf(rfp,"%d ", depth);
	}
	
	fprintf(rfp,"\n\ny <- scan()\n");
	for (i = 0; i < maxcov; i++) {
		fprintf(rfp,"%lu ", 1UL << i);
	}
	

	fprintf(rfp,"\n\n"
	"z <- t(matrix(log(data, 10),ncol=%d))\n"
	"pdf(file=\"%s.pdf\",height=8.3, width=11.7,onefile = TRUE)\n"
	, log_statistics->maxdepth -1
	, Rbasefilename
	);
	
	fwrite(stark_contour_log_stats_plot_R, stark_contour_log_stats_plot_R_len, 1, rfp);
	
	fclose(rfp);
	
}




