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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "main.h"
#include "stark.h"
#include "stark_status.h"
#include "terminal.h"


#include <readline/readline.h>
#include <readline/history.h>


int test(int* i) {
	return (*i)++;
}

// typedef struct { 
// 	int A : 1;
// 	int C : 1;
// 	int G : 1;
// 	int T : 1;
// 	 } bfield;
// 
// typedef union {
// struct { 
// 	int A : 1;
// 	int C : 1;
// 	int G : 1;
// 	int T : 1;
// 	 };
// 
// 	int mask : 4;
// } edgemask;

typedef struct  {
	char direction;
} testu;

#if defined(__INTEL_COMPILER)
//
// Disable ICC's remark #111: statement is unreachable
// This is legal ANSI C code so we disable the remark that is turned on with /W4
//
#pragma warning ( disable : 111 )

#pragma warning ( disable : 181 )

#endif


extern const char* ubigraph_url;

// void test_f() {
// 	typedef int test_t;
// }

struct twoints {
	int a[2];
};

void test_f(struct twoints a) {
	a.a[0] = 1;
	a.a[1] = 2;
}

#ifdef UBIGRAPHAPI_H
int ubigraph_start_attempt() {
	
	// attempt to find out if ubigraph_server is already running
	int pd[2];
	if(pipe(pd))
		return -1;
	
	int nullfd = open("/dev/null", O_RDWR);
	if (nullfd < 0)
		return -1;
	
	pid_t chldpid;
	if ((chldpid = fork())) {
		close(pd[1]);
		FILE *fp = fdopen(pd[0], "r");
		char buffer[1024];
		int ubigraph_is_running = 0;
		DEBUG_MSG("reading from ps (pid = %d)", chldpid);
		while (fgets(buffer, 1024, fp)) {
			if (strstr(buffer, "ubigraph_server"))
				return 0;
		}
		
		close(pd[0]);
		DEBUG_MSG("waiting for ps (pid = %d)", chldpid);
		
		int status;
		waitpid(chldpid, &status, 0);
		if (status) {
			return -1;
		}
	} else {
		dup2(nullfd,0);
		dup2(pd[1],1);
		// dup2(pd[1],2);
		char* argv1[] = {"ps", "-A", 0};
		execvp(argv1[0],argv1);
		perror("execvp");
		exit(-1);
	}
	
	close(pd[1]);
	close(pd[0]);
	
	// attempt to start ubigraph server
	if ((chldpid = fork())) {
		usleep(1000000);
		return 1;
	} else {
		dup2(nullfd,0);
		dup2(nullfd,1);
		dup2(nullfd,2);
		char* argv1[] = {"ubigraph_server", "-quiet", 0};
		execvp(argv1[0],argv1);
		perror("execvp");
		exit(-1);
	}
	
}
#endif

int main(int argc, char** argv) {
	
	// int i = 0;
	// printf("test() = %d\n",test(&i));
	// printf("i = %d\n",i);
	// 
	// 
	
	char* filename = NULL;
	
	for (++argv;*argv;argv++) {
		if (**argv == '-') {
			// switch((*argv)[1]) {
			// 	case 'k': minK = atoi(*(++argv)) -2; break;
			// 	case 'd': makedump = 1; break;
			// 	case 'h': help("stark"); exit(0); break;
			// 	case '-': if (!memcmp((*argv) + 2,"help",4)) {
			// 		help("stark");
			// 		exit(0);
			// 		} break;
			// 	case 0: {if (!filename)
			// 		filename = *argv;
			// 	} break;
			// }
		} else if (!filename)
			filename = *argv;
		else {
			fprintf(stderr,"Unknown option %s\n",*argv);
			// help("stark");
			exit(EXIT_FAILURE);
		}
	}
	
	fprintf(stderr,"sizeof(pthread_mutex_t) = %zx\n",sizeof(pthread_mutex_t));
	
	starknode_t node;
	
	fprintf(stderr,"StarK version %s\n", STARK_VERSION_STRING);
	
	// printf("MAX_CHAININS = %d\n", MAX_CHAININS);
	
	fprintf(stderr,"sizeof(starknode_t) = %zx\n",sizeof(starknode_t));
	#ifdef UBIGRAPHAPI_H
	fprintf(stderr,"sizeof(starknode_t.ubigraph) = %zx\n",sizeof(node.ubigraph));
	#endif
	//printf("sizeof(starknode_t.sequence) = %p\n",sizeof(node.sequence));
	
	// chainmultinode_t *testmultinode = calloc(1,sizeof(testmultinode->size) + (3* sizeof(chainnode_offset_t)));
	// testmultinode->size.next = 2;
	// testmultinode->size.previous = 5;
	// testmultinode->offset[0] = 1;
	// testmultinode->offset[1] = 2;
	// testmultinode->offset[2] = 3;
	// 
	// printf("testmultinode.size.next = %d\n",testmultinode->size.next);
	// printf("testmultinode.size.previous = %d\n",testmultinode->size.previous);
	// printf("testmultinode.offset = %p\n",testmultinode->offset - (unsigned long long)testmultinode);
	// printf("testmultinode.offset[0] = %d\n",testmultinode->offset[0]);
	// printf("testmultinode.offset[1] = %d\n",testmultinode->offset[1]);
	// printf("testmultinode.offset[2] = %d\n",testmultinode->offset[2]);
	
	#ifdef UBIGRAPHAPI_H
	printf("%s\n",ubigraph_url);
	// *((char**)(&ubigraph_url)) = "http://127.0.0.1:20738/RPC2";
	ubigraph_url = "http://127.0.0.1:20738/RPC2";
	printf("%s\n",ubigraph_url);
	ubigraph_start_attempt();
	#endif
	
	stark_t t;
	
	stark_init(&t);
	//getc(stdin);
	
	stark_status_init(0);
	stark_status_register_thread();
	
	if (filename) {
		stark_unserialize_file(filename, &t);
	}
	
	#ifdef INTERACTIVE
	
	// struct terminal_session session = {.stark = &t};
	// char buffer[1024];
	// while(1) {
	// 	terminal_get_command(buffer, 1024, stdin);
	// 	process_command_line(&session, buffer);
	// }
	
	char* input, *shell_prompt = "stark> ";
	struct terminal_session session = {.stark = &t};
	
	for (;;) {
		input = readline(shell_prompt);
		
		if (!input)
		   break;
		
		add_history(input);
		
		process_command_line(&session, input);
		
		// free(input);
	}
	
	#endif
	
	
	return 0;
}




