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
#include <sys/time.h>
#include <sys/resource.h>
#include <omp.h>
#include <getopt.h>

#include "main.h"
#include "stark.h"
#include "distributor.h"
#include "stark_status.h"
#include "stark_numa.h"
#include "fastx.h"
#include "stark_R.h"
#include "stark_assemble_hierarchical.h"

#ifdef UV_SYSTEM
#include "uv_system.h"
#endif

enum {
	FASTA, FASTQ
	};

// void help(char* name) {
// 	fprintf(stderr,"Usage: %s [-N] [-m maxK] [-k minK] [-d dump_filename] [-n num_threads] [filename]\n",name);
// 	fprintf(stderr,"\nDefault is to read from - (stdin).\n");
// }

void return_distributor_page(struct stark_insert_threaded_s* arg) {
	stark_set_status_action("returning distributor page");
	distributor_put(arg->distributor, arg);
}

static inline int parse_int_or_char(const char * str) {
	char * endptr;
	int result = strtol(str, &endptr, 10);
	if (str == endptr)
		return *str;
	else
		return result;
}

int stark_count_children_cb(struct starknode_navigator_s * nav, void * _children) {
	size_t (*children)[9] = _children;
	int i;
	int num_children = 0;
	
	for (i = 0; i < 8; i++) {
		if (OFFSET_VALID(nav->node->child[0][i])) {
			num_children++;
		}
	}
	
	children[nav->depth][num_children]++;
	return 0;
}


static char *dump_filename = 0;
static int minK = 30;
static int maxK = HARD_MAXDEPTH;
static int trimN = 1;
static int trimQ = 0;
static int num_threads = 0;
static int statusport = 0;
static char * load_stark_filename = NULL;
static int in_files = 0;
static int max_rounds = 0x10000000;
static int min_contig_length = 200;
static FILE * in[2] = {NULL, NULL};
static int memory_conserve = 0;
static size_t max_memory = 0x100000000000000;



static struct option options_array[] = {
	  {"save-stark", 1, NULL, 'd'}
	, {"min-k", 1, NULL, 'k'}
	, {"max-k", 1, NULL, 'm'}
	, {"no-trimN", 0, &trimN, 0}
	, {"trimQ", 1, NULL, 'T'}
	, {"num-threads", 1, NULL, 'n'}
	, {"status-port", 1, NULL, 'p'}
	, {"load-stark", 1, NULL, 'L'}
	, {"max-rounds", 1, NULL, 'r'}
	, {"min-contig-length", 1, NULL, 'l'}
	, {"help", 0, NULL, 'h'}
	, {"conserve-memory", 0, &memory_conserve, 1}
	, {"maximum-memory", 1, NULL, 'M'}
	, {NULL, 0, NULL, 0}
};

static struct {
	char * description;
	int * default_value;
	char * default_text;
	char * argument_string;
} help_string[] = {
	  {"Set the filename for the stark dump after phase 1.", NULL, "none", "<string>"}
	, {"set the minimum k-mer length to use for cleaning and assembly", &minK, NULL, "<int>"}
	, {"Set the maximum k-mer length to use for import", &maxK, NULL, "<int>"}
	, {"Disable N trimming, all N's will be replaced with A's", NULL, NULL, NULL}
	, {"Set the FASTQ Quality cutoff value", NULL, "disabled", "<int/char>"}
	, {"Set the number of threads", NULL, "auto", "<int>"}
	, {"Set the TCP port for HTTP status socket", NULL, "disabled", "<int>"}
	, {"Load stark from a phase 1 dump previously created with -d/--save-stark", NULL, "none", "<name>"}
	, {"Set the maximum number of assembly rounds", &max_rounds, NULL, "<int>"}
	, {"Set the minimum contig length to output", &min_contig_length, NULL, "<int>"}
	, {"Conserve Memory (experimental)", NULL, "disabled", NULL}
	, {"Maximum memory", NULL, "64 petabytes", "<int>"}
	, {NULL, NULL, NULL}
};

static void help() {
	fprintf(stderr, "Usage: stark [options] [filename1 [filename2]]\n\n\tReads input from FASTA or FASTQ from stdin or up to two files\n\tDuplicate filtering is performed when reading from two files.\n\n");
	
	int i;
	
	for (i = 0; help_string[i].description && options_array[i].name; i++) {
		if (!options_array[i].flag)
			fprintf(stderr, "\t-%c %s\n", options_array[i].val, help_string[i].argument_string ? help_string[i].argument_string : "");
		
		// fprintf(stderr, "\t--%s %s", options_array[i].name, help_string[i].argument_string ? help_string[i].argument_string : "");
		fprintf(stderr, "\t--%s%c%s\n", options_array[i].name, help_string[i].argument_string ? '=' : ' ', help_string[i].argument_string ? help_string[i].argument_string : "");
		fprintf(stderr, "\t\t%s\n", help_string[i].description);
		if (help_string[i].default_text)
			fprintf(stderr, "\t\tDefault: %s\n", help_string[i].default_text);
		
		if (help_string[i].default_value)
			fprintf(stderr, "\t\tDefault: %d\n", *help_string[i].default_value);
		
		fprintf(stderr, "\n");
		
	}
	
	// exit(0);
}


int main(int argc, char** argv) {
	
	in[0] = stdin;
	
	// FILE* fp = stdin;
	// 
	// if (argc > 1) {
	// 	fp = fopen(argv[1],"r");
	// }
	
	// char *filename = 0;

	
	// -Wunused-variables when compiling without threads
	// #ifdef THREADED
	// #endif
	
	
	// FILE * in[2] = {stdin, NULL};
	// FILE * load_stark = NULL;
	
	for (;;) {
		int c;
		int option_index = 0;
		
		c = getopt_long(argc, argv, "d:k:m:T:n:p:L:r:l:hcM:", options_array, &option_index);
		
		if (c == -1)
			break;
		
		switch (c) {
			case 'k': minK = atoi(optarg) -2; break;
			case 'm': maxK = atoi(optarg); break;
			// case 'N': trimN = 0; break;
			case 'T': trimQ = parse_int_or_char(optarg); break;
			case 'p': statusport = atoi(optarg); break;
			case 'n':
				#ifdef THREADED
				num_threads = atoi(optarg);
				#else
				fprintf(stderr,"Threads are not supported in this build. Please build with -DTHREADED or another threaded option to enable threads.\n");
				#endif
			break;
			case 'd': dump_filename = optarg; break;
			case 'h': help(); exit(0); break;
			case 'L': load_stark_filename = optarg; break;
			case 'l': min_contig_length = atol(optarg); break;
			case 'r': max_rounds = atol(optarg); break;
			case 'c': memory_conserve = 1; break;
			case 'M': max_memory = atol(optarg); break;
		}
		
	}
	
	for (;optind < argc; optind++) {
		if (!strcmp(argv[optind], "-")) 
			in[in_files++] = stdin;
		else if (in_files < 2) {
			if (!(in[in_files++] = fopen(argv[optind], "r"))) {
				perror("fopen");
				fprintf(stderr, "Cannot open file %s\n", argv[optind]);
				exit(EXIT_FAILURE);
			}
		}
		else {
			fprintf(stderr,"Unknown argument %s\n",argv[optind]);
			help();
			exit(EXIT_FAILURE);
		}
	}
	
	struct usage {
		struct rusage start, io, rmdup, insert, clean, dump, stats, assembly;
	} usage;
	
	
	getrusage(RUSAGE_SELF, &usage.start);
		
	#if defined(STARK_STATUS_H)
	stark_status_init(statusport);
	stark_status_register_thread();
	#endif
	
	char timebuffer[128];
	
	#ifdef THREADED_IMPORT
	if (!getenv("__DPLACE_")) {
		int num_cpus = stark_shed_init(num_threads);
		if (num_cpus && !num_threads)
			omp_set_num_threads(num_cpus);
		else
			omp_set_num_threads(num_threads);
		stark_numa_lock_master_thread();
	}
	#endif
	
	stark_t* stark = calloc(sizeof(stark_t), 1);
	
	struct coverage_log_statistics_s* log_statistics = NULL;
	
	if (load_stark_filename) {
		#if defined(STARK_STATUS_H)
		stark_set_status_action("Unserializing stark.");
		#endif
		
		stark_unserialize_file(load_stark_filename, stark);
		
		log_statistics = stark_statistics_log(stark);
		
		getrusage(RUSAGE_SELF, &usage.dump);
		
	} else {
		#ifdef COMPRESSED_DATASTRUCTURE
		if (memory_conserve) {
			stark_phase1_init(stark->phase1 = malloc(sizeof(*stark->phase1)), maxK, max_memory);
		} else
		#endif
			stark_init(stark);

		// NEW OpenMP Code

		stark_set_status_action("Reading reads. I/O");
		struct fastx_read_list_s * fastx_read_list = fastx_read_list_create(in[0], in[1]);

		if (in_files < 2) {
			fastx_read_list_scan_se(fastx_read_list, in[0], FASTX_NONAME | FASTQ_NOQUAL);
			getrusage(RUSAGE_SELF, &usage.io);

			print_timediff(timebuffer, &usage.io, &usage.start);
			fprintf(stderr,"Read I/O on %zu reads took %s\n", fastx_read_list->numreads, timebuffer);
		} else {
			fastx_read_list_scan_pe(fastx_read_list, in, FASTX_NONAME);
			getrusage(RUSAGE_SELF, &usage.io);

			print_timediff(timebuffer, &usage.io, &usage.start);
			fprintf(stderr,"Read I/O on %zu reads took %s\n", fastx_read_list->numreads, timebuffer);

			stark_set_status_action("Sorting, removing duplicate reads.");
			fastx_pair_unique(fastx_read_list, -1);
		}
		getrusage(RUSAGE_SELF, &usage.rmdup);
		

		if (in_files == 2) {
			print_timediff(timebuffer, &usage.rmdup, &usage.io);
			fprintf(stderr,"Remove read duplicates took %s\n", timebuffer);
		}

		// volatile size_t read_num_global = 0;

		#ifdef THREADED_IMPORT
		stark_coverage_counter_cache_init();
		#endif

		size_t discarded_reads = 0;
		
		#if defined(STARK_STATUS_H)
		
		// #pragma omp master
		{
			stark_status.reads_total = fastx_read_list->numreads;
		}
		
		#endif
		
		// omp_set_dynamic(0);
		struct stark_read_coverage_statistics_s stark_read_coverage_statistics[maxK];
		memset(stark_read_coverage_statistics, 0, sizeof(struct stark_read_coverage_statistics_s) * maxK);

		#pragma omp parallel	default(none) \
								firstprivate(fastx_read_list, stark, trimN, trimQ, maxK, minK) \
								reduction(+ : discarded_reads) \
								shared(stark_read_coverage_statistics, usage, timebuffer, stderr, dump_filename, log_statistics, memory_conserve)
		{
			
			stark_autopin_to_cpu(omp_get_thread_num());
			

			#if defined(STARK_STATUS_H)
			stark_status_register_thread();			
			#endif
			
			
			#ifdef STARK_COVERAGE_CACHE_TOKEN
			stark_coverage_token_init();
			#endif

			ssize_t i;
			#ifdef COMPRESSED_DATASTRUCTURE
			if (memory_conserve) {
				#pragma omp for schedule(dynamic,0x1000) nowait
				for (i = 0; i < fastx_read_list->numreads; i++) {
					struct fastx_read_s * fastx_read = fastx_read_list->pointers[i];
					struct fastx_read_s trimmed_read;
					stark_set_status_action("Trimming read");
					fastx_read_trim(&trimmed_read, fastx_read, trimN, trimQ);

					if (trimmed_read.length.seq <= minK) {
						discarded_reads += 1;
					} else {
						stark_set_status_action("Inserting read");
						// insert_sequence(stark, fastx_read->data + trimmed_read.offset.seq, trimmed_read.length.seq, maxK);
						// stark_phase1_insert_one_sequence
						
						process_sequence_incomplete(stark, fastx_read->data + trimmed_read.offset.seq, trimmed_read.length.seq, maxK, stark_phase1_insert_one_sequence, stark->phase1);
						
						// #ifdef STARK_COVERAGE_CACHE_TOKEN
						// stark_increment_coverage_cached_tokenized_tryflush();
						// #endif
					}
				}
				
				#pragma omp master
				stark_phase1_compartment_fill_status(stark->phase1);
				
				#ifdef STARK_COVERAGE_CACHE_TOKEN
				stark_increment_coverage_cached_tokenized_waitflush();
				// stark_increment_coverage_cached_tokenized_free();
				stark_coverage_token_release();
				#endif

				stark_coverage_counter_cache_flush_free(stark);

				#if defined(UV_SYSTEM)
				uv_free_atomic();
				#endif

				#if defined(STARK_STATUS_H)
				stark_set_status_action("Sleeping 1");
				#endif

				#pragma omp barrier

				#pragma omp master
				{
					getrusage(RUSAGE_SELF, &usage.insert);

					fprintf(stderr, "%zu reads discarded due to trimming.\n", discarded_reads);

					size_t memory = stark_memory_usage(stark);

					size_t numreads = fastx_read_list->numreads;

					fprintf(stderr,"Memory usage after inserting %zu reads: %ldMB\n",numreads , memory / 1048576);
					print_timediff(timebuffer, &usage.insert, &usage.rmdup);
					fprintf(stderr,"Inserting %zu reads took %s\n",numreads , timebuffer);
					// #ifdef DEBUG
					// fprintf(stderr,"Each call to insert_one_sequence takes on average %zu nanoseconds\n",(total_insert_mis * 1000) / total_insertions);
					// #endif


				}
				
				exit(0);
			} else
			#endif
			 {
				#pragma omp for schedule(dynamic,0x1000) nowait
				for (i = 0; i < fastx_read_list->numreads; i++) {
					struct fastx_read_s * fastx_read = fastx_read_list->pointers[i];
					struct fastx_read_s trimmed_read;
					stark_set_status_action("Trimming read");
					fastx_read_trim(&trimmed_read, fastx_read, trimN, trimQ);

					if (trimmed_read.length.seq <= minK) {
						discarded_reads += 1;
					} else {
						stark_set_status_action("Inserting read");
						insert_sequence(stark, fastx_read->data + trimmed_read.offset.seq, trimmed_read.length.seq, maxK);

						#ifdef STARK_COVERAGE_CACHE_TOKEN
						stark_increment_coverage_cached_tokenized_tryflush();
						#endif
					}
				}
				
				#ifdef STARK_COVERAGE_CACHE_TOKEN
				stark_increment_coverage_cached_tokenized_waitflush();
				// stark_increment_coverage_cached_tokenized_free();
				stark_coverage_token_release();
				#endif

				stark_coverage_counter_cache_flush_free(stark);

				#if defined(UV_SYSTEM)
				uv_free_atomic();
				#endif

				#if defined(STARK_STATUS_H)
				stark_set_status_action("Sleeping 1");
				#endif

				#pragma omp barrier

				#pragma omp master
				{
					getrusage(RUSAGE_SELF, &usage.insert);

					fprintf(stderr, "%zu reads discarded due to trimming.\n", discarded_reads);

					size_t memory = stark_memory_usage(stark);

					size_t numreads = fastx_read_list->numreads;

					fprintf(stderr,"Memory usage after inserting %zu reads: %ldMB\n",numreads , memory / 1048576);
					print_timediff(timebuffer, &usage.insert, &usage.rmdup);
					fprintf(stderr,"Inserting %zu reads took %s\n",numreads , timebuffer);
					// #ifdef DEBUG
					// fprintf(stderr,"Each call to insert_one_sequence takes on average %zu nanoseconds\n",(total_insert_mis * 1000) / total_insertions);
					// #endif


				}
			}

			

			#if defined(STARK_STATUS_H)
			stark_set_status_action("Cleaning stark");
			#endif
			stark_unambiguos_clean(stark, minK);

			#pragma omp barrier

			#pragma omp master
			{


				#if defined(STARK_STATUS_H)
				stark_set_status_action("Compacting stark");
				#endif

				stark_compact(stark);

				getrusage(RUSAGE_SELF, &usage.clean);
				print_timediff(timebuffer, &usage.clean, &usage.insert);
				fprintf(stderr,"Cleaning and compacting took %s\n", timebuffer);

				size_t memory = stark_memory_usage(stark);

				fprintf(stderr,"Memory usage after clean & compact: %ldMB\n",memory / 1048576);
			}


			#if defined(STARK_STATUS_H)
			stark_set_status_action("Sleeping 1");
			#endif

			#pragma omp barrier

			#pragma omp single
			{
				#if defined(STARK_STATUS_H)
				stark_set_status_action("dumping stark to file");
				#endif

				if (dump_filename)
					stark_serialize_file(dump_filename, stark);
				getrusage(RUSAGE_SELF, &usage.dump);
				print_timediff(timebuffer, &usage.dump, &usage.clean);
				if (dump_filename)
					fprintf(stderr,"Dumping stark to file took %s\n", timebuffer);
			}

			#if defined(STARK_STATUS_H)
			stark_set_status_action("Sleeping 2");
			#endif

			#pragma omp single
			{
				#if defined(STARK_STATUS_H)
				stark_set_status_action("Building coverage histogram");
				#endif

				// struct coverage_histogram_s* coverage_histogram = stark_statistics(stark, HARD_MAXDEPTH);
				log_statistics = stark_statistics_log(stark);
				
				fprintf(stderr,
					"Nodes total: %lu\n"
					"Sum coverage: %lu\n"
					"Maximum Coverage: %u\n"
					"Average Coverage: %u\n"
					, log_statistics->nodes_counted
					, log_statistics->sum_coverage
					, log_statistics->max_coverage
					, log_statistics->avg_coverage
					);


				FILE* stats_fp = fopen("log_statistics.tab","w");
				if (stats_fp) {

					stark_print_log_statistics(stats_fp, log_statistics);

					fclose(stats_fp);
				}

				#if defined(STARK_STATUS_H)
				stark_set_status_action("Writing R file 1");
				#endif

				init_R();
				stark_generate_coverage_contour_plot_R(log_statistics, getenv("PBS_JOBNAME") ? getenv("PBS_JOBNAME") : "coverage");
			}

			#if defined(STARK_STATUS_H)
			stark_set_status_action("Sleeping 3");
			#endif

			{
				#if defined(STARK_STATUS_H)
				stark_set_status_action("Building read link stats");
				#endif

				struct stark_read_coverage_statistics_s stark_read_coverage_statistics_private[maxK];
				memset(stark_read_coverage_statistics_private, 0, sizeof(struct stark_read_coverage_statistics_s) * maxK);

				#pragma omp for schedule(dynamic,0x1000) nowait
				for (i = 0; i < fastx_read_list->numreads; i++) {
					struct fastx_read_s * read = fastx_read_list->pointers[i];
					struct fastx_read_s trimmed_read;
					stark_set_status_action("Trimming read");
					fastx_read_trim(&trimmed_read, read, trimN, trimQ);

					if (trimmed_read.length.seq <= minK) {
						// discarded_reads += 1;
					} else {
						stark_set_status_action("Building read statistics");
						process_sequence_incomplete(stark, read->data + trimmed_read.offset.seq, trimmed_read.length.seq, maxK, stark_read_link_statistics_subf, stark_read_coverage_statistics_private);
					}
				}

				#pragma omp critical  (merge_stats) 
				{
					for (i = 0; i < maxK; i++) {
						stark_read_coverage_statistics[i].count += stark_read_coverage_statistics_private[i].count;
						stark_read_coverage_statistics[i].square_sum += stark_read_coverage_statistics_private[i].square_sum;
						stark_read_coverage_statistics[i].sum += stark_read_coverage_statistics_private[i].sum;

						size_t j;
						for (j = 0; j < ( sizeof(stark_read_coverage_statistics[0].hist) / sizeof(stark_read_coverage_statistics[0].hist[0])); j++) {
							stark_read_coverage_statistics[i].hist[j] += stark_read_coverage_statistics_private[i].hist[j];
						}
					}
				}

			}

			#if defined(STARK_STATUS_H)
			stark_set_status_action("Sleeping");
			#endif

			#pragma omp barrier
			#pragma omp master
			{

				fastx_read_list_free(fastx_read_list);

				FILE* stats_fp = fopen("read_connection_statistics.tab","w");
				if (stats_fp) {

					fprintf(stats_fp, "k\tcount\tmean\tvariance\t0-var\n");
					for (i = 2; i < maxK && stark_read_coverage_statistics[i].count; i++) {

						fprintf(stats_fp, "%zu\t%zu\t%f\t%f\t%f"
							, i
							, stark_read_coverage_statistics[i].count
							// , stark_read_coverage_statistics[i].sum
							// , stark_read_coverage_statistics[i].square_sum
							, stark_read_coverage_statistics[i].sum / stark_read_coverage_statistics[i].count
							, (stark_read_coverage_statistics[i].square_sum - ((stark_read_coverage_statistics[i].sum * stark_read_coverage_statistics[i].sum) / stark_read_coverage_statistics[i].count)) / stark_read_coverage_statistics[i].count
							, stark_read_coverage_statistics[i].square_sum / stark_read_coverage_statistics[i].count
						);

						size_t j;
						for (j = 0; j < ( sizeof(stark_read_coverage_statistics[0].hist) / sizeof(stark_read_coverage_statistics[0].hist[0])); j++) {
							fprintf(stats_fp, "\t%zu", stark_read_coverage_statistics[i].hist[j]);
						}
						fprintf(stats_fp, "\n");
					}

					fclose(stats_fp);
				}
				
				stats_fp = fopen("kmer_child_statistics.tab","w");
				if (stats_fp) {
					fprintf(stats_fp, "k");
					int i;
					for (i = 0; i <= 8; i++)
						fprintf(stats_fp, "\t%d", i);
					fprintf(stats_fp, "\tavg\tpp2\n");
					
					size_t children[HARD_MAXDEPTH][9];
					memset(children, 0, sizeof(children));
					
					struct starknode_navigator_s iterator = {
						  .stark = stark
						, .depth = 1
						, .offset = 1
						, .direction = 1
					};

					stark_navigate_loop_forall(&iterator, 1, stark_count_children_cb, children);
					
					depth_t depth;
					for (depth = 1; depth < HARD_MAXDEPTH; depth ++) {
						if (stark->size[depth]) {
							fprintf(stats_fp, "%d", depth);
							size_t sum = 0;
							size_t accu = 0;
							size_t twoorless = 0;
							for (i = 0; i <= 8; i++) {
								fprintf(stats_fp, "\t%zu", children[depth][i]);
								sum += children[depth][i];
								accu += i * children[depth][i];
								if (i <= 2)
									twoorless += children[depth][i];
							}
								
								
							fprintf(stats_fp, "\t%.2lf\t%.2lf\n", (double)accu / sum, (double)twoorless / sum);
						} else
							break;
					}
					
					fclose(stats_fp);
				}
				
				getrusage(RUSAGE_SELF, &usage.stats);
				print_timediff(timebuffer, &usage.stats, &usage.dump);
				fprintf(stderr,"Making statistics took %s\n", timebuffer);
			}

			stark_set_status_action("Sleeping or dead.");
		}
		
		// #if defined(STARK_STATUS_H)
		// stark_set_status_action("dumping stark to file");
		// #endif
		// 
		// if (dump_filename)
		// 	stark_serialize_file(dump_filename, stark);
		// getrusage(RUSAGE_SELF, &usage.dump);
		// print_timediff(timebuffer, &usage.dump, &usage.clean);
		// if (dump_filename)
		// 	fprintf(stderr,"Dumping stark to file took %s\n", timebuffer);


		
	}
	
	
	struct stark_hierarchical_assembler_s hierarchical_assembler;
	getrusage(RUSAGE_SELF, &usage.stats);
	
	#pragma omp parallel	default(none) \
							firstprivate(stark, maxK, minK, max_rounds, min_contig_length) \
							shared(hierarchical_assembler, usage, timebuffer, stderr, stdout, stark_status)
	{
	
		#if defined(STARK_STATUS_H)
		stark_status_register_thread();
		stark_set_status_action("Sleeping");
		#endif
	
		stark_hierarchical_assembler_init(stark, &hierarchical_assembler, minK);
		
	
		#pragma omp master
		{
			
			#if defined(STARK_STATUS_H)
			stark_status.hierarchical_assembler = &hierarchical_assembler;
			#endif
			
			fprintf(stderr,"Initialised hirarchial assembler with %zu groups.\n", hierarchical_assembler.groups.size);
		}
		
		size_t rounds = stark_hierarchical_assembler_test_and_merge_groups_openmp_v2(&hierarchical_assembler, max_rounds);
		
		
		#pragma omp master
		{
			fprintf(stderr,"Merged down to %zu groups in %zu rounds.\n", hierarchical_assembler.groups.size, rounds);
		}
		
		#if defined(STARK_STATUS_H)
		stark_set_status_action("Sleeping");
		#endif
		
		#pragma omp master
		{
		getrusage(RUSAGE_SELF, &usage.assembly);
		print_timediff(timebuffer, &usage.assembly, &usage.stats);
		fprintf(stderr,"Assembling sequences took %s\n", timebuffer);
		}
		
		#pragma omp single
		{
			#if defined(STARK_STATUS_H)
			stark_set_status_action("Printing contigs");
			#endif
			
			size_t bounds[2] = {
				min_contig_length
				, SIZE_MAX >> 1
			};
			stark_hierarchical_assembler_print_all_group_contigs(stdout, &hierarchical_assembler, bounds);
		
		}
		
		
	}
	
	return 0;
	
}

