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

#include "main.h"
#include "stark.h"
#include "distributor.h"
#include "stark_status.h"
#include "stark_numa.h"
#include "fastx.h"
#include "stark_R.h"
#include "stark_assemble.h"
#include "list.h"

#ifdef UV_SYSTEM
#include "uv_system.h"
#endif

enum {
	FASTA, FASTQ
	};

void help(char* name) {
	fprintf(stderr,"Usage: %s [-N] [-m maxK] [-k minK] [-d dump_filename] [-n num_threads] [filename]\n",name);
	fprintf(stderr,"\nDefault is to read from - (stdin).\n");
}

static inline int parse_int_or_char(const char * str) {
	char * endptr;
	int result = strtol(str, &endptr, 10);
	if (str == endptr)
		return *str;
	else
		return result;
}

static inline size_t print_timediff(char * buffer, struct rusage *a, struct rusage *b) {
	struct timeval t;
	timeradd(&a->ru_utime, &a->ru_stime, &t);
	timersub(&t, &b->ru_utime, &t);
	timersub(&t, &b->ru_stime, &t);
	
	size_t sec = t.tv_sec;
	int h, min;
	if (sec) {
		if ((h = (sec / (60 * 60)))) {
			min = (sec % (60 * 60)) / 60;
			sec %= 60;
			return sprintf(buffer, "%d:%d:%zu", h, min, sec);
		} else if ((min = sec / 60)) {
			sec %= 60;
			return sprintf(buffer, "%d:%zu", min, sec);
		} else 
			return sprintf(buffer, "%zusec", sec);
	} else {
		return sprintf(buffer, "%lums", (unsigned long)t.tv_usec / 1000);
	}
}


struct stark_kmer_reference_analysis_s {
	struct {
		size_t agree;
		size_t disagree;
		} coverage;
	struct {
		size_t agree;
		size_t disagree;
		} nodecount;
};

int stark_refcount_cb(struct starknode_navigator_s * nav, struct stark_kmer_reference_analysis_s * stats) {
	if (nav->node->flags & STARK_FLAG_EXTRACTED) {
		stats[nav->depth].coverage.agree += nav->node->coverage;
		stats[nav->depth].nodecount.agree ++;
		nav->node->flags &= ~STARK_FLAG_EXTRACTED;
	} else {
		stats[nav->depth].coverage.disagree += nav->node->coverage;
		stats[nav->depth].nodecount.disagree ++;
	}
	
	return 0;
}


int stark_refmark_subf(void * aux, stark_t* const stark, struct uplevelcache* const ___RESTRICT uplevel, const depth_t depth) {
	
	starknode_t * backparent, * forwardparent;
	// if (!uplevel[1].offset)
	// 	return 1;
	// 
	// if (!uplevel->offset)
	// 	return 0;
	
	starknode_t* const tmp = stark->level[depth-1];
	backparent = tmp + uplevel->offset;
	forwardparent = tmp + uplevel[1].offset;
	
	const int_fast8_t forwardchar_index = uplevel[depth-1].symbol_index;
	// const int_fast8_t backchar_complement_index = FOURSYMBOL_COMPLEMENT(uplevel->symbol_index);
	const int_fast8_t uplevel0reverse = uplevel->reverse;
	// const int_fast8_t uplevel1reverse = uplevel[1].reverse;
	
	offset_t offset = backparent->child[uplevel0reverse][forwardchar_index];
	
	if (!OFFSET_VALID(offset))
		return 1;
	
	unsigned char rindex; //, palindrome;
	
	rindex = 0;
	// palindrome = 0;
	{
		size_t i,j = depth-1;
		char d;
		for (i = 0; i < depth-1; i++) {
			d = uplevel[i].symbol_index - FOURSYMBOL_COMPLEMENT(uplevel[j].symbol_index);
			if (d < 0) {
				break;
			}
			else if (d > 0) {
				rindex = 1; // reverse
				break;
			}
			j--;
		}

		// palindrome = 1;
		// sequence is a palindrome
		
		
	}
	
	
	starknode_t * node = stark->level[depth] + offset;
	
	node->flags |= STARK_FLAG_EXTRACTED;
	
	uplevel->offset = offset;
	uplevel->reverse = rindex;
	// uplevel->palindrome = palindrome;
	
	return 0;
}



int main(int argc, char** argv) {
	
	FILE* fp = stdin;
	
	if (argc > 1) {
		fp = fopen(argv[1],"r");
	}
	
	// char *filename = 0;
	char *dump_filename = 0;
	int minK = 30;
	int maxK = HARD_MAXDEPTH;
	// bool makedump = 0;
	int num_threads = 0;
	int trimN = 1;
	int trimQ = 0;
	
	uint16_t statusport = 0;
	
	FILE * in[2] = {stdin, NULL};
	FILE * load_stark = NULL;
	char * load_stark_filename = NULL;
	int in_files = 0;
	
	double link_strength = 0.1;
	
	list_t(FILE *) refs;
	list_init(&refs);
	
	int reading_refs = 0;
	
	for (++argv;*argv;argv++) {
		if (**argv == '-') {
			switch((*argv)[1]) {
				case 'k': minK = atoi(*(++argv)) -2; break;
				case 'm': maxK = atoi(*(++argv)); break;
				// case 'N': trimN = 0; break;
				case 'T': trimQ = parse_int_or_char(*(++argv)); break;
				case 'p': statusport = atoi(*(++argv)); break;
				case 'n':
					#ifdef THREADED
					num_threads = atoi(*(++argv));
					#else
					fprintf(stderr,"Threads are not supported in this build. Please build with -DTHREADED or another threaded option to enable threads.\n");
					#endif
				break;
				case 'd': dump_filename = *(++argv); break;
				case 'h': help("stark"); exit(0); break;
				case 'r': reading_refs = 1; break;
				case 'L': load_stark_filename = *(++argv); break;
				case 'c': link_strength = atof(*(++argv)); break;
				case '-': if (!memcmp((*argv) + 2,"help",4)) {
					help("stark");
					exit(0);
					} break;
				case 0: {
					if (reading_refs) {
						list_insert(&refs, stdin);
					} else {
						if (in_files < 2)
							in[in_files++] = stdin;
						else {
							help("stark");
							exit(EXIT_FAILURE);
						}
					}
				} break;
			}
		} else {
			if (reading_refs) {
				FILE * fp = fopen(*argv, "r");
				if (!fp) {
					perror("fopen");
					fprintf(stderr, "Cannot open file %s\n", *argv);
					exit(EXIT_FAILURE);
				}
				list_insert(&refs, fp);
			}
			else if (in_files < 2) {
				if (!(in[in_files++] = fopen(*argv, "r"))) {
					perror("fopen");
					fprintf(stderr, "Cannot open file %s\n", *argv);
					exit(EXIT_FAILURE);
				}
			}
			else
				{ fprintf(stderr,"Unknown argument %s\n",*argv);
				help("stark");
				exit(EXIT_FAILURE); }
		}
	}
	
	/*
	FILE *in = 0;
	
	if (filename && *filename != '-') {
		in = fopen(filename,"r");
		if (!in) {
			perror("fopen"); exit(EXIT_FAILURE);
		}
	} else {
		in = stdin;
		filename = "-";
	}
	*/
	
	DEBUG_MSG("sizeof(starknode_t) = %d", sizeof(starknode_t));
	
	struct usage {
		struct rusage start, io, rmdup, insert, clean, dump, assembly;
	} usage;
	
	
	getrusage(RUSAGE_SELF, &usage.start);
		
	#if defined(STARK_STATUS_H)
	stark_status_init(statusport);
	stark_status_register_thread();
	#endif
	
	char timebuffer[128];
	
	#ifdef THREADED_IMPORT
	if (!getenv("__DPLACE_")) {
		int num_cpus = stark_shed_init();
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
		stark_init(stark);

		// NEW OpenMP Code

		stark_set_status_action("Reading reads. I/O");
		// struct fastx_read_list_s * fastx_read_list = fastx_read_list_create(NULL, NULL);

		// if (in_files < 2) {
		// 	fastx_read_list_scan_se(fastx_read_list, in[0], FASTX_NONAME | FASTX_NOQUAL);
		// 	getrusage(RUSAGE_SELF, &usage.io);
		// 
		// 	print_timediff(timebuffer, &usage.io, &usage.start);
		// 	fprintf(stderr,"Read I/O on %zu reads took %s\n", fastx_read_list->numreads, timebuffer);
		// } else {
		// 	fastx_read_list_scan_pe(fastx_read_list, in, FASTX_NONAME);
		// 	getrusage(RUSAGE_SELF, &usage.io);
		// 
		// 	print_timediff(timebuffer, &usage.io, &usage.start);
		// 	fprintf(stderr,"Read I/O on %zu reads took %s\n", fastx_read_list->numreads, timebuffer);
		// 
		// 	stark_set_status_action("Sorting, removing duplicate reads.");
		// 	fastx_pair_unique(fastx_read_list, -1);
		// }
		// getrusage(RUSAGE_SELF, &usage.rmdup);


		if (in_files >= 2) {
			print_timediff(timebuffer, &usage.rmdup, &usage.io);
			fprintf(stderr,"Remove read duplicates took %s\n", timebuffer);
		}

		// volatile size_t read_num_global = 0;

		#ifdef THREADED_IMPORT
		stark_coverage_counter_cache_init();
		#endif

		size_t discarded_reads = 0;

		// omp_set_dynamic(0);
		struct stark_read_coverage_statistics_s stark_read_coverage_statistics[maxK];
		memset(stark_read_coverage_statistics, 0, sizeof(stark_read_coverage_statistics));
		size_t numreads = 0;
		struct fastx_read_list_s * fastx_read_list;

		#pragma omp parallel	default(none) \
								firstprivate(fastx_read_list, stark, trimN, trimQ, maxK, minK) \
								reduction(+ : discarded_reads) \
								shared(stark_read_coverage_statistics, usage, timebuffer, stderr, dump_filename, log_statistics, numreads, in_files, in, num_threads)
		{
			#ifdef STARK_COVERAGE_CACHE_TOKEN
			stark_coverage_token_init();
			#endif

			#if defined(STARK_STATUS_H)
			stark_status_register_thread();
			#endif
			
			size_t read_chunk_size = num_threads ? 0x1000 * num_threads : 0x1000;

			// single read read
			if (in_files < 2) {
				getrusage(RUSAGE_SELF, &usage.io);
				
				fastx_read_list = fastx_read_list_create(NULL, NULL);
				size_t this_round_numreads;
				do {
					this_round_numreads = 0;
					stark_set_status_action("Waiting to do I/O");
					#pragma omp critical
					{
						stark_set_status_action("Reading reads");
						while (fastx_read_list_scan_read(fastx_read_list, in[0], FASTX_NONAME | FASTQ_NOQUAL) > 0 && this_round_numreads < read_chunk_size)
							this_round_numreads++;

						numreads += this_round_numreads;
					}
					
					ssize_t i;
					for (i = 0; i < fastx_read_list->numreads; i++) {
						struct fastx_read_s * read = fastx_read_list->pointers[i];
						struct fastx_read_s trimmed_read;
						stark_set_status_action("Trimming read");
						fastx_read_trim(&trimmed_read, read, trimN, trimQ);

						if (trimmed_read.length.seq <= minK) {
							// discarded_reads += 1;
						} else {
							stark_set_status_action("Inserting read");
							insert_sequence(stark, read->data + trimmed_read.offset.seq, trimmed_read.length.seq, maxK);

							#ifdef STARK_COVERAGE_CACHE_TOKEN
							stark_increment_coverage_cached_tokenized_tryflush();
							#endif
						}
					}
					
					fastx_read_list->numreads = 0;
					fastx_read_list->memory_used = 0;
				} while (this_round_numreads == read_chunk_size);
				fastx_read_list_free(fastx_read_list);

				print_timediff(timebuffer, &usage.io, &usage.start);
				// fprintf(stderr,"Read I/O on %zu reads took %s\n", fastx_read_list->numreads, timebuffer);
			} else {
				
				#pragma omp master
				{
					fastx_read_list = fastx_read_list_create(NULL, NULL);

					fastx_read_list_scan_pe(fastx_read_list, in, FASTX_NONAME);
					getrusage(RUSAGE_SELF, &usage.io);

					print_timediff(timebuffer, &usage.io, &usage.start);
					fprintf(stderr,"Read I/O on %zu reads took %s\n", fastx_read_list->numreads, timebuffer);

					stark_set_status_action("Sorting, removing duplicate reads.");
					fastx_pair_unique(fastx_read_list, -1);

					numreads = fastx_read_list->numreads;
				}
				
				ssize_t i;
				#pragma omp for schedule(dynamic,0x1000) nowait
				for (i = 0; i < fastx_read_list->numreads; i++) {
					struct fastx_read_s * read = fastx_read_list->pointers[i];
					struct fastx_read_s trimmed_read;
					stark_set_status_action("Trimming read");
					fastx_read_trim(&trimmed_read, read, trimN, trimQ);

					if (trimmed_read.length.seq <= minK) {
						discarded_reads += 1;
					} else {
						stark_set_status_action("Inserting read");
						insert_sequence(stark, read->data + trimmed_read.offset.seq, trimmed_read.length.seq, maxK);

						#ifdef STARK_COVERAGE_CACHE_TOKEN
						stark_increment_coverage_cached_tokenized_tryflush();
						#endif
					}
				}
				
				#pragma omp master
				fastx_read_list_free(fastx_read_list);
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

				

				fprintf(stderr,"Memory usage after inserting %zu reads: %ldMiB\n",numreads , memory >> 20);
				print_timediff(timebuffer, &usage.insert, &usage.io);
				fprintf(stderr,"Inserting %zu reads took %s\n",numreads , timebuffer);
				// #ifdef DEBUG
				// fprintf(stderr,"Each call to insert_one_sequence takes on average %zu nanoseconds\n",(total_insert_mis * 1000) / total_insertions);
				// #endif
				

			}

			// -----------------------------------
			// no clean & compact in this build
			// -----------------------------------
			
			/*
			#if defined(STARK_STATUS_H)
			stark_set_status_action("Cleaning stark");
			#endif
			stark_unambiguos_clean(stark);

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
			*/


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
			
			// NO LINK STATS HERE
			
			/*
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
						process_sequence(stark, read->data + trimmed_read.offset.seq, trimmed_read.length.seq, maxK, stark_read_link_statistics_subf, stark_read_coverage_statistics_private);
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
			}
			*/

			stark_set_status_action("Sleeping or dead.");
		}
		
	}
	
	
	{
		if (refs.size > 0) {
			printf("REFID\tK\tNODECOUNTAGREE\tNODECOUNTDISAGREE\t%%\tMASSAGREE\tMASSDISAGREE\t%%\n");
		}
		
		int i;
		for (i = 0; i < refs.size; i++) {
			#if defined(STARK_STATUS_H)
			char status[128];
			sprintf(status, "Reading in Reference %d", i);
			stark_set_status_action(status);
			#endif
			
			struct fastx_read_list_s * fastx_reference = fastx_read_list_create(refs.list[i], NULL);

			fastx_read_list_scan_se(fastx_reference, refs.list[i], FASTQ_NOQUAL | FASTX_NONAME);
			
			
			#if defined(STARK_STATUS_H)
			stark_set_status_action("Marking k-mers");
			#endif
			
			int j;
			for (j = 0; j < fastx_reference->numreads; j++) {
				struct fastx_read_s * reference = fastx_reference->pointers[j];

				process_sequence_incomplete(stark, reference->data + reference->offset.seq, reference->length.seq, maxK, stark_refmark_subf, NULL);
			}
			
			fastx_read_list_free(fastx_reference);
			
			struct starknode_navigator_s iterator = {
				  .stark = stark
				, .depth = 2
				, .offset = 1
				, .direction = 1
			};
			
			struct stark_kmer_reference_analysis_s refstats[HARD_MAXDEPTH];
			memset(refstats, 0, sizeof(refstats));
			
			stark_navigate_loop_forall(&iterator, 1, stark_refcount_cb, refstats);
			
			for (j = 2; j < HARD_MAXDEPTH && j <= maxK && stark->size[j]; j++) {
				printf("%d\t%d\t%zu\t%zu\t%.1lf%%\t%zu\t%zu\t%.1lf%%\n"
					, i
					, j
					, refstats[j].nodecount.agree
					, refstats[j].nodecount.disagree
					, (double)(100.0 * refstats[j].nodecount.agree) / (double)(refstats[j].nodecount.disagree + refstats[j].nodecount.agree)
					, refstats[j].coverage.agree
					, refstats[j].coverage.disagree
					, (double)(100.0 * refstats[j].coverage.agree) / (double)(refstats[j].coverage.disagree + refstats[j].coverage.agree)
				);
			}
			
		}
	}
	
	list_free(refs);
	
	stark_free(stark);
	free(stark);
	if (log_statistics)
		free(log_statistics);
	
	return 0;
	
	#if defined(STARK_STATUS_H)
	stark_set_status_action("Assembling");
	#endif
	
}

