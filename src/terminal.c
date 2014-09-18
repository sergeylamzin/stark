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


#include "terminal.h"
#include "stark_assemble_hierarchical.h"
#include <stdint.h>

int parse_args(char** argv, char* data) {
	char c;
	int argc = 0;
	// int length = strlen(data);
	// char* args = malloc(length +1);
	// memcpy(args, data, length +1);
	// char** argv = calloc(length +1,1);
	
	if (argv)
		*argv = data;
	
	int quoteopen = 0;
	int waitingnextarg = 1;
	for (;(c = *data); data++) {
		switch(c) {
			case '\n':
			case '\t':
			case ' ':
				if (!quoteopen) {
					if (argv)
						*data = 0;
					waitingnextarg = 1;
				}
			break;
			case '"':
				if (argv)
					*data = 0;
				if (quoteopen) {
					waitingnextarg = 1;
				}
				quoteopen ^= 1;
			break;
			default:
				if (waitingnextarg) {
					waitingnextarg = 0;
					if (argv)
						argv[argc] = data;
					argc++;
				}
			break;
		}
	}
	
	return argc;
}


typedef int (*terminal_func_t)(struct terminal_session*, int, char**);

struct command {
	const char* name;
	const char* alias;
	terminal_func_t function;
	const char* description;
};

const struct command commands[];

terminal_func_t select_str(char*);
int terminal_example_next_line(struct terminal_session* , int , char** );



int terminal_get_command(char* buffer, int maxlength, FILE* fp) {
	if (isatty(fileno(fp)))
		printf("> ");
	fgets(buffer, maxlength, fp);
	return EXIT_SUCCESS;
}

int process_command_line(struct terminal_session* session, char* input) {	
	// char** argv;
	int argc = parse_args(0, input);
	char* argv[argc+1];
	parse_args(argv, input);
	
	// printf("argc = %d\n",argc);
	// int i;
	// for (i = 0; i < argc; i++)	{
	// 	printf("'%s'\n",argv[i]);
	// }
	
	if (argc) {
		terminal_func_t func = select_str(argv[0]);
		if (func)
			return func(session,--argc, argv+1);
	
		printf("Unknown command %s. Type ? for help.\n",input);
	} else if (session->example_fp) {
		terminal_example_next_line(session, 0, NULL);
	}
	return EXIT_SUCCESS;
}

#ifdef UBIGRAPHAPI_H
int stark_showhide_level (struct terminal_session* session, char show, int level) {
	const char* visibility = show ? "true" : "false";
	
	session->stark->level_shown[level] = show;
	
	if (level) {
		ubigraph_set_vertex_style_attribute(level, "visible", visibility);
		ubigraph_set_edge_style_attribute(level * 2, "visible", visibility);
	}
	else
		ubigraph_set_vertex_attribute(session->stark->root.ubigraph.id, "visible", visibility);
		
	if (level && (!show || session->stark->level_shown[level - 1]))
		ubigraph_set_edge_style_attribute(level * 2 - 1, "visible", visibility);
	if (!show || session->stark->level_shown[level + 1])
		ubigraph_set_edge_style_attribute(level * 2 + 1, "visible", visibility);
	
	return EXIT_SUCCESS;
}

// inline int64_t tuple(int32_t x, int32_t y);

int terminal_stark_showhide_wrap(struct terminal_session* session, char* data, char show) {
	// char c,*pos;
	// int i, end;
	// for (pos = data; (c = *pos); pos++) {
	// 	if (c == ' ' || (end = (c == '\0'))) {
	// 		*pos = '\0';
	// 		if (!strcmp(data, "all")) {
	// 			for (i = 0; i < HARD_MAXDEPTH; i++)
	// 				stark_showhide_level(session, show, i);
	// 		} else {
	// 			int level = atoi(data);
	// 			stark_showhide_level(session, show, level);
	// 		}
	// 		if (!end) {
	// 			data = pos+1;
	// 		}
	// 	}
	// }
	
	int i;
	
	if (!strcmp(data, "all")) {
		// int i;
		for (i = 0; i < HARD_MAXDEPTH; i++)
			stark_showhide_level(session, show, i);
	} else if (!strcmp(data, "chain")) {
		ubigraph_set_edge_style_attribute(-2, "visible", show ? "true" : "false");
	} else if (!strcmp(data, "cleaned")) {
		ubigraph_set_vertex_attribute(session->stark->root.ubigraph.id, "visible", "false");
		ubigraph_set_edge_style_attribute(1, "visible", "false");
		int depth, depth2, j;
		
		for (i = 0; i < HARD_MAXDEPTH; i++)
			stark_showhide_level(session, 0, i);
		
		for (depth = 1; session->stark->size[depth]; depth++) {
			// printf("Level %d:\n",depth);
			for (i = 1; i <= session->stark->size[depth]; i++) {
				starknode_t* node = session->stark->level[depth] + i;
				
				node->ubigraph.visible = 0;
				
				if (node->depth) {
					if (!(node->child[0][0] || node->child[0][1] || node->child[0][2] || node->child[0][3]) || !(node->child[1][0] || node->child[1][1] || node->child[1][2] || node->child[1][3])) {
						node->ubigraph.visible = 1;
						ubigraph_set_vertex_attribute(node->ubigraph.id, "visible", "true");
					} // else {
					// 						ubigraph_set_vertex_attribute(node->ubigraph.id, "visible", "false");
					// 						for (depth2 = depth-1; depth2 <= depth +1; depth2++) {
					// 							// printf("Level %d:\n",depth);
					// 							for (j = 1; j <= session->stark->size[depth2]; j++) {
					// 								starknode* node2 = session->stark->level[depth2] + j;
					// 
					// 								if (node2->depth){
					// 										ubigraph_set_edge_attribute(tuple(node2->ubigraph.id, node->ubigraph.id), "visible", "false");
					// 								}
					// 
					// 							}
					// 						}
					// 					}
						
				
					
				}
				// vertex_id_t vid;
				// for (vid = 1; vid <= session->stark->ubicount; vid++) {
				// 	ubigraph_set_edge_attribute(cantor(vid, node->ubigraph.id), "visible", visible ? "true" : "false");
				// 	ubigraph_set_edge_attribute(cantor(node->ubigraph.id, vid), "visible", visible ? "true" : "false");
				// }
					
				
			}
		}
		
		for (depth = 1; session->stark->size[depth]; depth++) {
			for (i = 1; i <= session->stark->size[depth]; i++) {
				if (session->stark->level[depth][i].ubigraph.visible) {
					for (depth2 = depth-1; depth2 <= depth +1; depth2++) {
						for (j = 1; j <= session->stark->size[depth2]; j++) {
							if (session->stark->level[depth2][j].ubigraph.visible)
								ubigraph_set_edge_attribute(tuple(session->stark->level[depth][i].ubigraph.id, session->stark->level[depth2][j].ubigraph.id), "visible", "true");
						}
					}
				}
			}
		}
		
	} else {
		int level = atoi(data);
		stark_showhide_level(session, show, level);
	}
	return EXIT_SUCCESS;
}
#endif


#if defined(__INTEL_COMPILER)
//
// Disable ICC's remark #869: parameter was never referenced
// This is legal ANSI C code so we disable the remark that is turned on with /W4
//
#pragma warning ( disable : 869 )

#endif

// commands
int help(struct terminal_session* session, int argc, char** argv);
int quit(struct terminal_session* session, int argc, char** argv) {
	exit(0);
	return EXIT_SUCCESS;
}
int terminal_insert_sequence(struct terminal_session* session, int argc, char** argv) {
	if (argc) {
		// printf("Inserting %s\n",*argv);
		size_t length = strlen(*argv);
		insert_sequence(session->stark, *argv, length, length);
	}
	else
		printf("Cannot insert empty sequence.\n");
	return EXIT_SUCCESS;
}
int terminal_ubigraph_clean(struct terminal_session* session, int argc, char** argv) {
	int minK = 0;
	if (argc) {
		minK = atoi(*argv);
	}
	stark_unambiguos_clean(session->stark, minK);
	return EXIT_SUCCESS;
}
int terminal_stark_reset(struct terminal_session* session, int argc, char** argv) {
	stark_free(session->stark);
	stark_init(session->stark);
	return EXIT_SUCCESS;
}
int terminal_stark_print(struct terminal_session* session, int argc, char** argv) {
	if (argc)
		stark_print_depth(session->stark, atoi(*argv));
	else
		stark_print(session->stark,HARD_MAXDEPTH);
	return EXIT_SUCCESS;
}
int terminal_reverse(struct terminal_session* session, int argc, char** argv) {
	if (argc) {
		char* seq = *argv;
		reverse_complement_onspot(seq, strlen(seq));
		printf("%s\n", seq);
	}
	return EXIT_SUCCESS;
}
int terminal_stark_extract(struct terminal_session* session, int argc, char** argv) {
	if (argc) {
		char seq[1024];
		depth_t depth = atoi(argv[0]);
		offset_t offset = atoi(argv[1]);
		stark_extract_sequence(seq, session->stark, depth, offset);
		printf("%s\n",seq);
	} else {
		stark_extract_all(session->stark);
	}
	
	return EXIT_SUCCESS;
}
#ifdef ENCHAIN_SEQUENCES
int terminal_stark_extract_chain(struct terminal_session* session, int argc, char** argv) {
	if (argc) {
		char seq[1024];
		depth_t depth = atoi(argv[0]);
		offset_t offset = atoi(argv[1]);
		stark_start_chain(seq, session->stark, depth, offset);
		printf("%s\n",seq);
	} else {
		stark_extract_all(session->stark);
	}
	
	return EXIT_SUCCESS;
}
#endif
int terminal_stark_compact(struct terminal_session* session, int argc, char** argv) {
	stark_compact(session->stark);
	return EXIT_SUCCESS;
}
int terminal_stark_find(struct terminal_session* session, int argc, char** argv) {
	if (argc) {
		depth_t depth;
		offset_t offset;
		int result = stark_find(&depth, &offset, session->stark, *argv);
		
		if (result) {
			printf("Depth = %d\nOffset = %" PRI_OFFSET "\n", depth, offset);
		} else
			printf("Node does not exist.\n");
		
	} else
		printf("You have to provide a sequence to look for.\n");
	return EXIT_SUCCESS;
}
int terminal_stark_save(struct terminal_session* session, int argc, char** argv) {
	if (argc) {
		
		stark_serialize_file(*argv, session->stark);
		
	} else {
		printf("You have to specify a filename.\n");
	}
	
	return EXIT_SUCCESS;
}
int terminal_stark_load(struct terminal_session* session, int argc, char** argv) {
	if (argc) {
		
		stark_unserialize_file(*argv, session->stark);
		
		// if (stark) {
		// 	stark_free(session->stark);
		// 	session->stark = stark;
		// 	
		// }
		// else
		// 	printf("Failed to load stark.\n");
		
	} else {
		printf("You have to specify a filename.\n");
	}
	
	return EXIT_SUCCESS;
}
#ifdef UBIGRAPHAPI_H
int terminal_ubigraph_rebuild(struct terminal_session* session, int argc, char** argv) {
	stark_ubigraph_rebuild(session->stark);
	return EXIT_SUCCESS;
}
int terminal_stark_hide(struct terminal_session* session, int argc, char** argv) {
	for (; argc; argc--) {
		terminal_stark_showhide_wrap(session, *argv, 0);
		argv++;
	}
	return EXIT_SUCCESS;
}
int terminal_stark_show(struct terminal_session* session, int argc, char** argv) {
	for (; argc; argc--) {
		terminal_stark_showhide_wrap(session, *argv, 1);
		argv++;
	}
	return EXIT_SUCCESS;
}
int terminal_stark_mark(struct terminal_session* session, int argc, char** argv) {
	
	depth_t depth = atoi(argv[0]);
	offset_t offset = atoi(argv[1]);
	
	vertex_id_t vid = session->stark->level[depth][offset].ubigraph.id;
	//ubigraph_set_vertex_attribute(vid, "size", "2.0");
	if (argc > 2)
		ubigraph_set_vertex_attribute(vid, "color", argv[2]);
	
	return EXIT_SUCCESS;
}
int terminal_stark_mark_sequence(struct terminal_session* session, int argc, char** argv) {
	
	depth_t depth;
	offset_t offset;
	int found;
	
	while (argc) {
		
		found = stark_find(&depth, &offset, session->stark, *argv);
		
		if (found) {
			vertex_id_t vid = session->stark->level[depth][offset].ubigraph.id;
			ubigraph_set_vertex_attribute(vid, "size", "1.5");
			ubigraph_set_vertex_attribute(vid, "color", "#ffffff");
		}
			
		argc--;
		argv++;
	}
	
	
	
	
	
	return EXIT_SUCCESS;
}
int terminal_stark_mark_chain(struct terminal_session* session, int argc, char** argv) {
	
	depth_t depth;
	offset_t offset;
	vertex_id_t prev_vid = 0;
	int found;
	
	while (argc) {
		
		found = stark_find(&depth, &offset, session->stark, *argv);
		
		if (found) {
			vertex_id_t vid = session->stark->level[depth][offset].ubigraph.id;
			ubigraph_set_vertex_attribute(vid, "size", "1.5");
			ubigraph_set_vertex_attribute(vid, "color", "#ffffff");
			// ubigraph_set_vertex_attribute(vid, "visible", "true");
			
			if (prev_vid) {
				edge_id_t edge = ubigraph_new_edge(prev_vid, vid);
				ubigraph_change_edge_style(edge, -2);
				
			}
			
			prev_vid = vid;
		} else {
			
			prev_vid = 0;
		}
		
		
		argc--;
		argv++;
	}
	
	
	
	
	
	return EXIT_SUCCESS;
}


#endif

int terminal_example_next_line(struct terminal_session* session, int argc, char** argv) {
	
	if (!session->example_fp) {
		printf("There is currently no example loaded. Please use 'example' to load an example.\n");
		return -1;
	}
	
	char* readsuccess = NULL;
	
	if (!session->example_next_line) {
		session->example_next_line = malloc(2048);
		readsuccess = fgets(session->example_next_line, 2048, session->example_fp);
	}
	
	do {
		process_command_line(session, session->example_next_line);
		readsuccess = fgets(session->example_next_line, 2048, session->example_fp);
	} while (readsuccess && *session->example_next_line != '\n');
	
	// if (session->example_next_line) {
	// 	process_command_line(session, session->example_next_line);
	// } else {
	// 	session->example_next_line = malloc(2048);
	// }	
	// 
	// 
	// 
	// for (;bytesread = fgets(session->example_next_line, 2048, session->example_fp);) {
	// 	if (*session->example_next_line != '\n')
	// 		break;
	// } 
	
	while (readsuccess && *session->example_next_line == '\n')
		readsuccess = fgets(session->example_next_line, 2048, session->example_fp);
	
	if (!readsuccess) {
		printf("End of example.\n");
		fclose(session->example_fp);
		session->example_fp = NULL;
		free(session->example_next_line);
		session->example_next_line = NULL;
		
	} else {
		printf("next example command: %s", session->example_next_line);
	}
	
	
	return EXIT_SUCCESS;
}

int terminal_load_example(struct terminal_session* session, int argc, char** argv) {
	
	if (!argc) {
		printf("You have to specify a filename\n");
		return -1;
	} 
	
	if (session->example_fp)
		fclose(session->example_fp);
	if (session->example_next_line)
		free(session->example_next_line);
	
	FILE* fp = fopen(argv[0],"r");
	
	if (!fp) {
		printf("Cannot open file %s\n",argv[0]);
		return -2;
	} 
	
	session->example_fp = fp;
	
	terminal_example_next_line(session, 0, NULL);
	
	
	return EXIT_SUCCESS;
}

int terminal_stark_extract2(struct terminal_session* session, int argc, char** argv) {
	if (argc >= 3) {
//		char* sequence;
		depth_t depth = atoi(argv[1]);
		offset_t offset = atoi(argv[2]);
		int minK = atoi(argv[0]);
		
		stark_assemble_single(session->stark, minK, depth, offset);
	} else {
		// stark_extract_all(session->stark);
	}
	
	return EXIT_SUCCESS;
}

int terminal_stark_statistics(struct terminal_session* session, int argc, char** argv) {
	
	struct coverage_histogram_s* coverage_histogram = stark_statistics(session->stark, HARD_MAXDEPTH);
	
	printf("Nodes total: %lu\n", coverage_histogram->nodes_counted);
	printf("Sum coverage: %lu\n", coverage_histogram->sum_coverage);
	printf("Maximum Coverage: %u\n", coverage_histogram->max_coverage);
	printf("Average Coverage: %u\n", coverage_histogram->avg_coverage);
	
	if (argc > 0) {
		FILE* fp = fopen(*argv,"w");
		if (fp) {
			uint32_t coverage;
			for (coverage = 0; coverage <= coverage_histogram->max_coverage; coverage++)
				fprintf(fp,"%u\t%u\n", coverage, coverage_histogram->hist[coverage]);
				
			fclose(fp);
		}
	} 
	
	//free(coverage_histogram);
	
	return EXIT_SUCCESS;
}

int terminal_presentation_osascript_send_key_code(struct terminal_session* session, int argc, char** argv) {
	if (argc > 1) {
		char buffer[1024] = "osascript";
		size_t length = sizeof("osascript");
		char* fork_argv[2+4+(2*(argc-1))];
		int i, idx = 0;
		fork_argv[idx++] = buffer;
		fork_argv[idx++] = buffer + length;
		length += 1 + sprintf(buffer + length, "-e");
		fork_argv[idx++] = buffer + length;
		length += 1 + sprintf(buffer + length, "tell application \"%s\" to activate", argv[0]);
		fork_argv[idx++] = buffer + length;
		for (i = 1; i < argc; i++) {
			length += 1 + sprintf(buffer + length, "-e");
			fork_argv[idx++] = buffer + length;
			length += 1 + sprintf(buffer + length, "tell application \"System Events\" to key code %d", atoi(argv[i]));
			fork_argv[idx++] = buffer + length;
		}
		length += 1 + sprintf(buffer + length, "-e");
		fork_argv[idx++] = buffer + length;
		length += 1 + sprintf(buffer + length, "tell application \"Terminal\" to activate");
		fork_argv[idx++] = NULL;
		if (!fork()) {
			execvp(*fork_argv, fork_argv);
			perror("execvp");
			exit(-1);
		} else
			wait(0);
	} else {
		printf("You need to specify a process and a key code\n");
	}
	
	return EXIT_SUCCESS;
}

int terminal_skim_go_to_page(struct terminal_session* session, int argc, char** argv) {
	if (argc >= 2) {
		char buffer[256];
		snprintf(buffer, 256, "tell document \"%s\" of application \"Skim\" to go to page %d", argv[0], atoi(argv[1]));
		pid_t pid = fork();
		if (!pid) {
			char * args[] = {"osascript", "-e", buffer, NULL};
			execvp(*args, args);
			perror("execvp");
			exit(-1);
		} else if (pid > 0)
			wait(0);
		else
			perror("fork");
	} else {
		printf("You need to specify an open document and a page number\n");
	}
	
	return EXIT_SUCCESS;
}

int terminal_stark_get_neighbour(struct terminal_session* session, int argc, char** argv) {
	if (argc >= 3) {
		offset_t offset;
		depth_t depth;
		
		int result = stark_find(&depth, &offset, session->stark, *argv);
		
		if (!result) 
			printf("Node does not exist.\n");
		
		
		
		char direction = argv[1][0] == 'r' ? 1 : 0;
		char ch = foursymbol_index[argv[2][0]];
		
		offset_t neighbour_offset = stark_get_neighbour(session->stark, offset, depth, direction, ch);
		
		printf("Your neigbour is at offset %d\n", neighbour_offset);
		
	} else {
		printf("Usage: neighbour NODE f|r nucleotide\n");
	}
	
	return 0;
}

#include "stark_navigate.h"

void stark_substract_coverage_single(stark_t* const stark, depth_t depth, offset_t offset, coverage_t amount);

int terminal_stark_test_uplink(struct terminal_session* session, int argc, char** argv) {
	if (argc >= 2) {
		offset_t offset;
		depth_t depth;
		
		int result = stark_find(&depth, &offset, session->stark, *argv);
		
		if (!result) {
			printf("Node does not exist.\n");
			return -1;
		}
		
		stark_substract_coverage_single(session->stark, depth, offset, atoi(argv[1]));
		
	} else {
		printf("Usage: neighbour NODE f|r nucleotide\n");
	}
	return 0;
}

int terminal_stark_test_neighbour(struct terminal_session* session, int argc, char** argv) {
	if (argc >= 2) {
		offset_t offset;
		depth_t depth;
		
		int result = stark_find(&depth, &offset, session->stark, *argv);
		
		if (!result) {
			printf("Node does not exist.\n");
			return -1;
		}
		
		
		struct starknode_navigator_s nav = {.stark = session->stark, .depth = depth, .offset = offset, .direction = (argv[1][0] == 'r' ? 1 : 0), .node = session->stark->level[depth] + offset};
		
		// char direction = argv[1][0] == 'r' ? 1 : 0;
		char ch = foursymbol_index[argv[2][0]];
		
		offset_t * offsets[2];
		
		if(stark_navigate_shared_link_pair(offsets, &nav, ch)) {
			printf("Navigate failed\n");
		}
		
		printf("Your the two offsets are [%p, %p] = [%d, %d]\n", offsets[0], offsets[1], offsets[0][0], offsets[1][0]);
		
		
		/*
		offset_t * parent_offset_to_me = stark_navigate_parent_offset_to_self(session->stark, depth, offset, argv[1][0] == 'r' ? 1 : 0);
		
		printf("*(int *)%p = %d\n", parent_offset_to_me, *parent_offset_to_me);
		*/
		
	} else {
		printf("Usage: neighbour NODE f|r nucleotide\n");
	}
	return 0;
}

void stark_substract_coverage(stark_t* const stark, offset_t offset, depth_t depth, coverage_t amount, int_fast8_t share_forward, int_fast8_t share_back);

int terminal_stark_substract_coverage(struct terminal_session* session, int argc, char** argv) {
	if (argc >= 2) {
		offset_t offset;
		depth_t depth;
		
		int result = stark_find(&depth, &offset, session->stark, *argv);
		
		if (!result) {
			printf("Node does not exist.\n");
			return -1;
		}
		
		
		stark_substract_coverage(session->stark, offset, depth, atoi(argv[1]), -1, -1);
		
	} else {
		printf("Usage: substract NODE amount\n");
		return -1;
	}
	return 0;
}

// int terminal_stark_hierarchical_assembler_get_sequences(struct terminal_session* session, int argc, char** argv) {
// 	if (argc >= 1) {
// 		stark_hierarchical_assembler_get_sequences(session->stark, atoi(*argv));
// 		return 0;
// 	} else {
// 		printf("You have to specify minK\n");
// 		return -1;
// 	}
// 	return 0;
// }

// int terminal_stark_test_dfs(struct terminal_session* session, int argc, char** argv) {
// 	stark_test_dfs(session->stark);
// 	
// 	return 0;
// }

int terminal_stark_test_ha_grouping(struct terminal_session* session, int argc, char** argv) {
	if (argc >= 1) {
		printf("Currently unavailable\n");
		// stark_hierarchical_assembler_test(session->stark, atoi(*argv));
		return 0;
	} else {
		printf("You have to specify minK\n");
		return -1;
	}
}
int terminal_stark_hierarchical_assembler_init(struct terminal_session* session, int argc, char** argv) {
	if (argc >= 1) {
		// printf("Currently unavailable\n");
		session->hierarchical_assembler = malloc(sizeof(session->hierarchical_assembler[0]));
		stark_hierarchical_assembler_init (session->stark, session->hierarchical_assembler, atoi(*argv));
		return 0;
	} else {
		printf("You have to specify minK\n");
		return -1;
	}
}

int terminal_stark_hierarchical_assembler_test_and_merge_groups_alias(struct terminal_session* session, int argc, char** argv) {
	if (!session->hierarchical_assembler) {
		puts("Hirarchial assembler has not been initialised. Use hainit to initalise.");
		return -1;
	}
	int max_rounds = 0x10000;
	if (argc)
		max_rounds = atoi(*argv);
	
	size_t rounds = stark_hierarchical_assembler_test_and_merge_groups (session->hierarchical_assembler, max_rounds);
	printf("Merged groups in %zu rounds\n", rounds);
	return 0;
}

int terminal_stark_hierarchical_assembler_test_and_merge_groups(struct terminal_session* session, int argc, char** argv) {
	if (argc >= 1 && !strcmp(*argv, "groups")) {
		// printf("Currently unavailable\n");
		return terminal_stark_hierarchical_assembler_test_and_merge_groups_alias(session, argc - 1, argv+1);
	} else {
		printf("Unknown command\nUse help merge\n");
		return -1;
	}
}

int terminal_stark_hierarchical_assembler_print_all_group_contigs(struct terminal_session* session, int argc, char** argv) {
	if (!session->hierarchical_assembler) {
		puts("Hirarchial assembler has not been initialised. Use hainit to initalise.");
		return -1;
	}
	
	size_t bounds[2] = {0, INTPTR_MAX};
	if (argc > 0) 
		bounds[0] = atoi(*argv);
	if (argc > 1) 
		bounds[1] = atoi(argv[1]);
	
	stark_hierarchical_assembler_print_all_group_contigs (stdout, session->hierarchical_assembler, bounds);
	
	return 0;
}

const struct command commands[] = {
	{"help", "?", help, "Displays this help."},
	{"exit", "q", quit, "Exit."},
	{"insert", "i", terminal_insert_sequence, "Insert a genome sequence into the graph."},
	{"clean", "", terminal_ubigraph_clean, "Remove unambiguous nodes from the graph."},
	#ifdef UBIGRAPHAPI_H
	{"rebuild", "", terminal_ubigraph_rebuild, "Rebuild the Graph."},
	#endif
	{"reset", "", terminal_stark_reset, "Remove all nodes from the graph and reset it."},
	{"print", "", terminal_stark_print, "Print the Graph."},
	{"extract", "e", terminal_stark_extract, "Extract sequence."},
	#ifdef ENCHAIN_SEQUENCES
	{"echain", "ec", terminal_stark_extract_chain, "Extract chain sequence."},
	#endif
	// {"test", "", terminal_stark_test_dfs, "Test function, changes."},
	{"substract", "sub", terminal_stark_substract_coverage, "Substract coverage."},
	{"neighbour", "ne", terminal_stark_get_neighbour, "Print offset or a node's neigthbour."},
	{"reverse", "rev", terminal_reverse, "Print reverse complement of input sequence."},
	{"compact", "", terminal_stark_compact, "Compact the graph in Memory."},
	{"save", "", terminal_stark_save, "Save Graph to file. Arguments: filename"},
	{"load", "", terminal_stark_load, "Load Graph from file. Arguments: filename"},
	{"find", "", terminal_stark_find, "Find a sequence in the graph."},
	{"extract2", "e2", terminal_stark_extract2, "Extract-Assemble. Arguments mikK, startdepth, startoffset"},
	{"stats", "", terminal_stark_statistics, "Generate statistics. Arguments: histogram output file"},
	{"example", "ex", terminal_load_example, "Load command file with example commands."},
	{"next", "n", terminal_example_next_line, "Next example command."},
	{"osascript", "osa", terminal_presentation_osascript_send_key_code, "Send key code (second argument) to process (first argument) via applescript."},
	{"skimpage", "sp", terminal_skim_go_to_page, "tell document \"[arg1]\" of application \"Skim\" to go to page [arg2]"},
	#ifdef UBIGRAPHAPI_H
	{"hide", "h", terminal_stark_hide, "Hide Nodes. Arguments: all or level number"},
	{"show", "s", terminal_stark_show, "Show Nodes. Arguments: all, cleaned or level number"},
	{"mark", "m", terminal_stark_mark, "Marks a node making it big. Optionally specify an RGB encoded color preceeded by #."},
	{"marks", "ms", terminal_stark_mark_sequence, "Marks sequences making them big."},
	{"markc", "mc", terminal_stark_mark_chain, "Marks sequence chain making nodes big and creating edges."},
	#endif
	#ifdef HIERARCHIAL_ASSEMBLY
	// {"hags", "", terminal_stark_hierarchical_assembler_get_sequences, "stark_hierarchical_assembler_get_sequences"},
	{"hainit", "hai", terminal_stark_hierarchical_assembler_init, "Initialise Hirarchial assembler"},
	{"merge", "", terminal_stark_hierarchical_assembler_test_and_merge_groups, "Hirarchial assembler merge"},
	{"mg", "", terminal_stark_hierarchical_assembler_test_and_merge_groups_alias, "Hirarchial assembler merge groups"},
	{"hacontigs", "hac", terminal_stark_hierarchical_assembler_print_all_group_contigs, "Hirarchial assembler merge groups"},
	#endif
};

int help(struct terminal_session* session, int argc, char** argv) {
	int i;
	printf("%-10s %-5s DESCRIPTION\n","COMMAND", "ALIAS");
	for (i=0; i < sizeof commands/sizeof *commands;i++) {
		printf("%-10s %-5s %s\n",commands[i].name,commands[i].alias,commands[i].description);
	}
		
	return 0;
}

terminal_func_t select_str(char *s) {
	int i;
	for (i=0; i < sizeof(commands)/sizeof(struct command);i++) {
		if (!strcmp(s, commands[i].name) || !strcmp(s, commands[i].alias))
			return commands[i].function;
	}
		
	return NULL;
}









