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


#include "main.h"
#include "stark.h"
#include <fcntl.h>
#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>
#include <errno.h>
#include <arpa/inet.h>

int check_endianess () {
	if (htonl(42) == 42)
		return 0; // BIG ENDIAN
	else
		return 1; // LITTLE ENDIAN
}

uint64_t (*htonll)(uint64_t);
uint64_t (*ntohll)(uint64_t);

uint64_t big_endian_htonll (uint64_t value) {
	return value;
}

uint64_t little_endian_htonll (uint64_t value) {
	const uint32_t high_part = htonl(value >> 32);
	const uint32_t low_part = htonl(value & 0xFFFFFFFFLL);

	return (((uint64_t)low_part << 32) | high_part);
}

uint64_t little_endian_ntohll (uint64_t value) {
	const uint32_t high_part = ntohl(value >> 32);
	const uint32_t low_part = ntohl(value & 0xFFFFFFFFLL);

	return (((uint64_t)low_part << 32) | high_part);
}

// __builtin_bswap64



const char CONTAINER_DEF[] = "STARK\n"; // DO NOT CHANGE THIS

struct stark_block_s {
	char type; // 
	/*
		0 for final block     // introduced in stark 0.3.0, allows detection of stark serializer end for concatenated input
		1 for option
		2 for data identifier
		3 for metadata
	*/
	char identifier[15];
	
	union {
		char padding[48];
		char metadata[48];
		
		struct {
			int32_t size; // option size in bytes
			int32_t offset; // option offset in block

		} option;
		
		struct {
			int64_t length; // data length in bytes
			int64_t blocksize; // block size in bytes
			int32_t size; // amount of elements in block
			int32_t padding;
			
			union {
				struct {
					int32_t depth;
				} starknode;
			} type;
		} data_identifier;
	} payload;
	
};

/*
struct stark_block_s {
	char type;
	char padding[63];
};

struct stark_metadata_s {
	char type; // 3 for metadata, 1 for option, 2 for data identifier
	char identifier[15]; // null terminated string identifier
	char data[48];
};

struct stark_option_s {
	char type; // 3 for metadata, 1 for option, 2 for data identifier
	char identifier[15]; // null terminated string identifier
	int32_t size; // option size in bytes
	int32_t offset; // option offset in block
	//char padding[40]; // for future format expansions
};

struct stark_data_identifier_s {
	char type; // 1 for option, 2 for data identifier
	char identifier[15]; // null terminated string identifier
	int32_t size; // amount of elements in block
	int64_t length; // data length in bytes
	int64_t blocksize; // block size in bytes
	//char padding[28]; // for future format expansions
};

struct stark_data_identifier_starknode_s {
	char type; // 1 for option, 2 for data identifier
	char identifier[15]; // null terminated string identifier
	int32_t size; // amount of elements in block
	int64_t length; // data length in bytes
	int64_t blocksize; // block size in bytes
	int32_t depth;
	//char padding[28]; // for future format expansions
};
*/

struct starknode_serialized_s {
	offset_t uplink[2];
	coverage_t coverage;
	unsigned char edges;
	unsigned char flags;
};


struct stark_block_s supported_options[] = {
	{1,		"uplink[0]",	{.option.size = sizeof(offset_t)}		},
	{1,		"uplink[1]",	{.option.size = sizeof(offset_t)}		},
	{1,		"coverage",		{.option.size = sizeof(coverage_t)}		},
	{1,		"edges",		{.option.size = sizeof(unsigned char)}	},
	{1,		"flags",		{.option.size = sizeof(unsigned char)}	},
};

static __inline__ void init_serializer() {
	if (htonl(42) == 42){
		htonll = big_endian_htonll;
		ntohll = big_endian_htonll;
	}
	else {
		htonll = little_endian_htonll;
		ntohll = little_endian_ntohll;
	}
	
	
	
	int i;
	int offset = 0;
	for (i = 0; i < (sizeof(supported_options) / sizeof(struct stark_block_s)); i++) {
		supported_options[i].payload.option.offset = offset;
		offset += supported_options[i].payload.option.size;
		//memset(&(supported_options[i].padding), 0, sizeof(supported_options[i].padding))
	}
	
}

int stark_serialize_fp(FILE* fp, stark_t* stark) {
	if (!check_endianess ()) {
		fprintf(stderr,"Serialization is only supported on little-endian processor architectures.\n");
		return 1;
	}
	
	init_serializer();
/*	
	int fd;
	
	fd = open(filename, O_WRONLY | O_CREAT | O_TRUNC);
	
	if (fd < 0) {
		fprintf(stderr,"Unable to open file %s\n",filename);
		return -1;
	}
	*/
	
	
	// int i, r;
	
	//r = write(fd, CONTAINER_DEF, writebytes = sizeof(CONTAINER_DEF) -1);
	
	if (!fwrite(CONTAINER_DEF, sizeof(CONTAINER_DEF) -1, 1, fp)) {
		fprintf(stderr,"Unable to write to stream\n");
		return -1;
	}
	
	
	struct stark_block_s stark_block;
	
	// struct stark_metadata_s stark_metadata;
	memset(&stark_block, 0, sizeof(stark_block));
	stark_block.type = 3;
	sprintf(stark_block.identifier, "%s", "version");
	sprintf(stark_block.payload.metadata, "%s", STARK_VERSION_STRING);
	//memcpy(&stark_block, &stark_metadata, sizeof(stark_metadata));
	
	// r = write(fd, &stark_block, writebytes = sizeof(stark_block));
	
	if (!fwrite(&stark_block, sizeof(stark_block), 1, fp)) {
		fprintf(stderr,"Unable to write to stream\n");
		return -1;
	}
	
	// write contained options
	/*
	for (i = 0; i < (sizeof(supported_options) / sizeof(struct stark_block_s)); i++) {
		//memset(&stark_block, 0, sizeof(stark_block));
		//memcpy(&stark_block, supported_options + i, sizeof(struct stark_option_s));
		
		
		r = write(fd, supported_options + i, sizeof(struct stark_block_s));

		if (r != writebytes) {
			fprintf(stderr,"Unable to write %d bytes to file %s\n",writebytes,filename);
			return -1;
		}
	}
	*/
	if (fwrite(supported_options, sizeof(stark_block), sizeof(supported_options) / sizeof(struct stark_block_s), fp) != sizeof(supported_options) / sizeof(struct stark_block_s)) {
		fprintf(stderr,"Unable to write to stream\n");
		return -1;
	}
	
	
	memset(&stark_block, 0, sizeof(stark_block));
	stark_block.type = 2;
	sprintf(stark_block.identifier, "%s", "starknodes");
	
	// struct stark_data_identifier_starknode_s stark_data_identifier_starknode;
	// stark_data_identifier_starknode.type = 2;
	// memcpy(stark_data_identifier_starknode.identifier, "starknodes", sizeof("starknodes"));
	// 
	
	
	struct starknode_serialized_s serialized;
	
	starknode_t starknode;
	depth_t depth;
	offset_t offset;
	
	for (depth = 1; stark->size[depth]; depth++) {
		// memset(&stark_block, 0, sizeof(stark_block));
		stark_block.payload.data_identifier.size = stark->size[depth];
		stark_block.payload.data_identifier.blocksize = sizeof(struct starknode_serialized_s);
		stark_block.payload.data_identifier.length = stark_block.payload.data_identifier.size * stark_block.payload.data_identifier.blocksize;
		stark_block.payload.data_identifier.type.starknode.depth = depth;
		
		//memcpy(&stark_block, &stark_data_identifier_starknode, sizeof(stark_data_identifier_starknode));
		// r = write(fd, &stark_block, writebytes = sizeof(struct stark_block_s));

		if (!fwrite(&stark_block, sizeof(stark_block), 1, fp)) {
			fprintf(stderr,"Unable to write to stream\n");
			return -1;
		}
		
		for (offset = 1; offset <= stark->size[depth]; offset++) {
			
			starknode = stark->level[depth][offset];
			
			serialized.uplink[0] = starknode.uplink[0]; //htonoffset_t(starknode.uplink[0]);
			serialized.uplink[1] = starknode.uplink[1]; //htonoffset_t(starknode.uplink[1]);
			
			serialized.coverage = starknode.coverage; //htoncoverage_t(starknode.coverage);
			serialized.edges = starknode.edges;
			serialized.flags = starknode.flags;
			
			// r = write(fd, &serialized, writebytes = sizeof(serialized));
			// 
			// if (r != writebytes) {
			// 	fprintf(stderr,"Unable to write %d bytes to file %s\n",writebytes,filename);
			// 	return -1;
			// }
			
			if (!fwrite(&serialized, sizeof(serialized), 1, fp)) {
				fprintf(stderr,"Unable to write to stream\n");
				return -1;
			}
			
		}
	}
	
	memset(&stark_block, 0, sizeof(stark_block));
	if (!fwrite(&stark_block, sizeof(stark_block), 1, fp)) {
		fprintf(stderr,"Unable to write to stream\n");
		return -1;
	}
	
	return 0;
}

/*
int readb(int fildes, void *buf, size_t nbyte) {
	int r, bytesleft = nbyte;
	
	while (bytesleft > 0) {
		r = read(fildes, buf, bytesleft);
		if (r <= 0)
			return r;
		
		bytesleft -= r;
		buf += r;
	}
	return nbyte;
}
*/

uint64_t intmask[] = {0, 0xFF, 0xFFFF, 0xFFFFFF, 0xFFFFFFFF, 0xFFFFFFFFFF, 0xFFFFFFFFFFFF, 0xFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF};

int stark_unserialize_fp(FILE* fp, stark_t* stark) {
	
	if (!check_endianess ()) {
		fprintf(stderr,"Serialization is only supported on little-endian processor architectures.\n");
		return 1;
	}
	
	// fprintf(stderr, "unserializing STARK on fd = %d\n",fd);
	
	int i;
	depth_t depth;
	offset_t offset;
			
	int j, blocks_read;
			
	struct stark_block_s reading_options[sizeof(supported_options) / sizeof(struct stark_block_s)];
	
	memcpy(reading_options, supported_options, sizeof(supported_options));
	
	for (i = 0; i < (sizeof(supported_options) / sizeof(struct stark_block_s)); i++) {
		reading_options[i].payload.option.size = 0;
		reading_options[i].payload.option.offset = 0;
	}
	// read options and metadata
	struct stark_block_s stark_block;
	
	starknode_t starknode;
	// struct starknode_serialized_s serialized;
		
	stark_free(stark);
	stark_init(stark);
	
	while (fread ( &stark_block, sizeof(stark_block), 1, fp )) {

		// if (b = 0; b < blocks_read; b++) {
			switch(stark_block.type) {
				case 0:
				DEBUG_MSG("Detected END-of-StarK, returning");
				return 0;
				break;
				case 1: // option

				for (i = 0; i < (sizeof(supported_options) / sizeof(struct stark_block_s)); i++) {
					if (!memcmp(reading_options[i].identifier, stark_block.identifier, strlen(reading_options[i].identifier))) {
						memcpy(reading_options + i, &stark_block, sizeof(*reading_options));
						if (reading_options[i].payload.option.size > 8)
							reading_options[i].payload.option.size = 8;
						if (reading_options[i].payload.option.size < 0)
							reading_options[i].payload.option.size = 0;
					}
				}

				break;
				case 2: // data block

					if (!memcmp("starknodes", stark_block.identifier, 10)) {
						/*
						struct stark_data_identifier_starknode_s {
							char type; // 1 for option, 2 for data identifier
							char identifier[15]; // null terminated string identifier
							int32_t size; // amount of elements in block
							int64_t length; // data length in bytes
							int64_t blocksize; // block size in bytes
							int32_t depth;
							//char padding[28]; // for future format expansions
						}
						*/


						depth = stark_block.payload.data_identifier.type.starknode.depth;

						// DEBUG_MSG("Loading Depth %d",depth);

						if (depth > 1) {
							if(!(stark->level[depth] = calloc(sizeof(starknode_t), stark_block.payload.data_identifier.size + 1))) {
								CRITICAL_ERROR("Out of Memory!");
							}
						}
						stark->size[depth] = stark_block.payload.data_identifier.size;
						stark->maxsize[depth] = stark_block.payload.data_identifier.size;


						char buffer[1024][stark_block.payload.data_identifier.blocksize];
						
						offset = 1;
						while (stark_block.payload.data_identifier.size) {
							blocks_read = fread( buffer, stark_block.payload.data_identifier.blocksize, 1023 < stark_block.payload.data_identifier.size ? 1023 : stark_block.payload.data_identifier.size, fp );
							if (!blocks_read){
								fprintf(stderr, "Error unserializing StarK, stream is corrupt of invalid format.");
								return -2;
								break;
							}
							
							stark_block.payload.data_identifier.size -= blocks_read;
							
							for (i = 0; i < blocks_read; i++) {

								memset(&starknode, 0, sizeof(starknode));
								j = 0;
								starknode.uplink[0] = *(uint64_t*)(buffer[i] + reading_options[j].payload.option.offset) & intmask[reading_options[j].payload.option.size]; j++;
								starknode.uplink[1] = *(uint64_t*)(buffer[i] + reading_options[j].payload.option.offset) & intmask[reading_options[j].payload.option.size]; j++;
								starknode.coverage = *(uint64_t*)(buffer[i] + reading_options[j].payload.option.offset) & intmask[reading_options[j].payload.option.size]; j++;
								starknode.edges = *(uint64_t*)(buffer[i] + reading_options[j].payload.option.offset) & intmask[reading_options[j].payload.option.size]; j++;
								starknode.flags = *(uint64_t*)(buffer[i] + reading_options[j].payload.option.offset) & intmask[reading_options[j].payload.option.size]; j++;
								starknode.depth = depth;

								// char buffer[1024];
								// stark_print_node_info(buffer, &starknode);
								// //buffer[1023] = 0;
								// printf("[%d," OFFSET_F "] : %s\n", depth, offset, buffer);

								stark->level[depth][offset] = starknode;

								int childnum = (starknode.flags & STARK_PARENT1_LINK_MASK) >> 4;
								if(STARK_NODE_IS_PALINDROME(stark->level[depth-1][starknode.uplink[0]])) {
									stark->level[depth-1][starknode.uplink[0]].child[0][childnum] = offset;
									stark->level[depth-1][starknode.uplink[0]].child[1][childnum] = offset;
								}
								else
									stark->level[depth-1][starknode.uplink[0]].child[starknode.flags & STARK_FLAG_PARENT1_REVERSE ? 1 : 0][childnum] = offset;

								childnum = (starknode.flags & STARK_PARENT2_LINK_MASK);
								if(STARK_NODE_IS_PALINDROME(stark->level[depth-1][starknode.uplink[1]])) {
									stark->level[depth-1][starknode.uplink[1]].child[0][childnum] = offset;
									stark->level[depth-1][starknode.uplink[1]].child[1][childnum] = offset;
								}
								else
									stark->level[depth-1][starknode.uplink[1]].child[starknode.flags & STARK_FLAG_PARENT2_REVERSE ? 0 : 1][childnum] = offset;

								offset++;
							}
							
						}
						
						

					}

				break;
				case 3: // metadata

				if (!memcmp("version", stark_block.identifier, 7)) {
					fprintf(stderr, "Unserializing StarK version %s\n", stark_block.payload.metadata);
				}

				break;
				//default: // unknown block
			}
		// }
	}
	
	DEBUG_MSG("stark unserialize reaches EOF");
	
	return 0;
	
}

int stark_serialize_file(char* filename, stark_t* stark) {
	FILE* fp = fopen(filename, "w");
	
	if (!fp) {
		fprintf(stderr,"Unable to open file %s\n",filename);
		perror("fopen");
		return -1;
	}
	
	stark_serialize_fp(fp, stark);
	fclose(fp);
	
	return 0;
	
}

int stark_unserialize_file(char* filename, stark_t* stark) {
	FILE* fp = fopen(filename, "r");
	
	if (!fp) {
		fprintf(stderr,"Unable to open file %s\n",filename);
		perror("fopen");
		return -1;
	}
		
	char buffer[16];
	fgets(buffer, 16, fp);
	
	if (memcmp(buffer, CONTAINER_DEF, sizeof(CONTAINER_DEF) -1)) {
		fprintf(stderr,"Input format is invalid\n");
	}
	
	stark_unserialize_fp(fp, stark);
	fclose(fp);
	
	return 0;
	
}


