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
#include <stdint.h>
#include <string.h>
#include "fastx.h"

#define USAGE_STR \
"Usage: " "\n" \
"\n" \
"\t" "stark_fasta_rmdup [options]" "\n" \
"\n" \
"\t\t" "-h\t" "This help. You found it!.\n" \
"\t\t" "-i\t" "Next arguments that are not preceeded with a - (dash) are input filenames. A - (dash) stands for standard in. Up to two filenames are allowed. If omitted standard in will be used.\n" \
"\t\t" "-o\t" "Next arguments that are not preceeded with a - (dash) are output filenames. A - (dash) stands for standard out. Up to two filenames are allowed. If omitted standard out will be used.\n" \
"\t\t" "-Q CHAR/INT\t" "Join quality scores of duplicate reads using the provided score as a minimum offset. B for Illumina.\n" \
"\n" \
"Examples:\n"


// "\t\t" "-m INT\t" "Minimum read length for trimmed output. Nop if not used with -N -Q or -T. Default: half length.\n" 
// "\t\t" "-N\t" "Trim away preceeding or terminating N's. The next argument that is not preceeded with a - (dash) is an output filename for reads that do not pass the trimming criteria. May be omitted in which case those reads are not output. Shares output file with -T\n" 
// "\t\t" "-T CHAR/INT\t" "Trim away heads or tails of nucleotides with quality scores lower then pecified. The Quality score may be specified eigther as an integer or ASCII character encoding. The next argument that is not preceeded with a - (dash) is an output filename for reads that do not pass the trimming criteria. May be omitted in which case those reads are not output. Shares output file with -N\n"
// "stark_fasta_rmdup -i sequence_1.txt sequence_2.txt -o out_1.fq out_2.fq -N -T D discarded_reads.fq -Q B\n" 
// "stark_fasta_rmdup -o out_1.fq out_2.fq -N -T 58 -m 30\n (reads from standard in assuming paired end reads come one after another)\n"


#define PRINT_USAGE() fprintf(stderr, USAGE_STR)

static inline int parse_int_or_char(const char * str) {
	char * endptr;
	int result = strtol(str, &endptr, 10);
	if (str == endptr)
		return *str;
	else
		return result;
}

int main(int argc, char ** argv) {
	
	int trimN = 0;
	int trimQ = 0;
	int joinQ = -1;
	int32_t min_output_length = INT32_MAX;
	
	int in_files = 0;
	int out_files = 0;
	
	FILE * fq[2] = {stdin, NULL};
	FILE * out[2] = {stdout, NULL};
	FILE * trimfp = NULL;
	
	enum {IN, OUT, TRIM} reading_args = IN;
	
	for (++argv;*argv;argv++) {
		if (**argv == '-') {
			switch((*argv)[1]) {
				case 'i': reading_args = IN; break;
				case 'o': reading_args = OUT; break;
				// case 'N': trimN = 1; reading_args = TRIM; break;
				// case 'T': trimQ = parse_int_or_char(*(++argv)); reading_args = TRIM; break;
				case 'Q': joinQ = parse_int_or_char(*(++argv)); break;
				// case 'm': min_output_length = atoi(*(++argv)); break;
				case 'h': PRINT_USAGE(); exit(1); break;
				case '-': if (!memcmp((*argv) + 2,"help",4)) {
					PRINT_USAGE(); exit(1);
					} break;
				case '\0': {
					switch (reading_args) {
						case IN:
						if (in_files < 2)
							fq[in_files++] = stdin;
						else
							{ PRINT_USAGE(); exit(2); }
						break;
						case OUT:
						if (out_files < 2)
							out[out_files++] = stdout;
						else
							{ PRINT_USAGE(); exit(3); }
						break;
						case TRIM:
						if (!trimfp) {
							trimfp = stdout;
							fprintf(stderr, "Warning: outputting reads that don't pass the trim filter to standard out.\n");
						}
						else
							{ PRINT_USAGE(); exit(3); }
						break;
						default:
						PRINT_USAGE(); exit(4);
						break;
					}
				} break;
			}
		} else {
			switch (reading_args) {
				case IN:
				if (in_files < 2) {
					if (!(fq[in_files++] = fopen(*argv, "r"))) {
						perror("fopen");
						fprintf(stderr, "Cannot open file %s\n", *argv);
					}
				}
				else
					{ PRINT_USAGE(); exit(2); }
				break;
				case OUT:
				if (out_files < 2) {
					if (!(out[out_files++] = fopen(*argv, "w"))) {
						perror("fopen");
						fprintf(stderr, "Cannot open file %s\n", *argv);
					}
				}
				else
					{ PRINT_USAGE(); exit(3); }
				break;
				case TRIM:
				if (!trimfp) {
					if (!(trimfp = fopen(*argv, "w"))) {
						perror("fopen");
						fprintf(stderr, "Cannot open file %s\n", *argv);
					}
				}
				else
					{ PRINT_USAGE(); exit(3); }
				break;
				default:
				PRINT_USAGE(); exit(4);
				break;
			}
			// if (reading_args == IN) {
			// 	if (in_files < 2) {
			// 		if (!(fq[in_files++] = fopen(*argv, "r"))) {
			// 			perror(__FILE__ ":" TOSTRING(__LINE__) ": " "fopen");
			// 			fprintf(stderr, "Cannot open file %s\n", *argv);
			// 		}
			// 	}
			// 	else { PRINT_USAGE(); exit(5); }
			// } else 	if (reading_args == OUT) {
			// 	if (out_files < 2) {
			// 		if (!(out[out_files++] = fopen(*argv, "w"))) {
			// 			perror(__FILE__ ":" TOSTRING(__LINE__) ": " "fopen");
			// 			fprintf(stderr, "Cannot open file %s\n", *argv);
			// 		}
			// 	}
			// 	else {
			// 		PRINT_USAGE(); exit(6);
			// 	}
			// }
		}
	}
	
	if (!fq[1])
		fq[1] = fq[0];
	
	if (!out[1])
		out[1] = out[0];
	
	
	
	struct fastx_read_list_s * fastx_read_list = fastx_read_list_create();
	fastx_read_list_scan_pe(fastx_read_list, fq);
	
	fclose(fq[0]);
	if (fq[0] != fq[1])
		fclose(fq[1]);
	
	fprintf(stderr, "read in %zu reads, %zu pairs\n", fastx_read_list->numreads, fastx_read_list->numreads >> 1);
	fprintf(stderr, "Using quality joining on minimum %c ~ %d\n", joinQ, joinQ);
	
	fastx_pair_unique(fastx_read_list, joinQ);
	
	fprintf(stderr, "%zu unique pairs left\n", fastx_read_list->numreads >> 1);
	
	fastx_read_list_print_pe(fastx_read_list, out);
	
	/*
	fprintf(stderr, "%zu unique pairs left\n", stark_read_list->numpairs);
	
	stark_print_reads(stark_read_list, out[0], out[1]);
	*/
	fclose(out[0]);
	if (out[0] != out[1])
		fclose(out[1]);
	
	return 0;
}

