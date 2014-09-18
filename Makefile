    # Copyright (C) 2014  Sergey Lamzin, https://github.com/sergeylamzin/stark

    # This file is part of the StarK genome assembler.

    # StarK is free software: you can redistribute it and/or modify it
    # under the terms of the GNU General Public License as published by
    # the Free Software Foundation, either version 3 of the License, or
    # (at your option) any later version.

    # StarK is distributed in the hope that it will be useful,
    # but WITHOUT ANY WARRANTY; without even the implied warranty of
    # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    # GNU General Public License for more details.

    # You should have received a copy of the GNU General Public License
    # along with Foobar.  If not, see <http://www.gnu.org/licenses/>.


# the compiler to use.
# CC=gcc
COMPILER=$(shell which cc)
CC=$(shell readlink $(COMPILER))
OFFSETL=32
CFLAGS += -std=gnu99 -Iinclude -I. -DOFFSET_L=$(OFFSETL) -D_GNU_SOURCE $(FLAGS)
DEBUGFLAGS=-g -O0 -Wall
LDFLAGS= -Wall

#-DDEBUG_SEQUENCES=32 
# UBILIBS=libs/resource.o -Llibs -L../libs -lm -lcurl  -lxmlrpc_client -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok -lxmlrpc -lubigraphclient
UBILIBS=-Llibs -L../libs -lm -lcurl -lxmlrpc_client -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok -lxmlrpc -lubigraphclient
UVLIBS=-lgru_alloc -lgru
# LIBS=$(GCCLIBS) $(ICCLIBS) $(PBS_LIBS)
VPATH = src:libs:include:sqlite3

ifneq (,$(findstring clang,$(CC)))
CFLAGS += -Wall -Wno-char-subscripts -Winline -finline-functions -funroll-loops
LDLIBS += -lpthread -ldl -lm
CFLAGS += -fopenmp
endif

ifneq (,$(findstring gcc,$(CC)))
CFLAGS += -Wall -Wno-char-subscripts -Winline -finline-functions -funroll-loops
LDLIBS += -lpthread -ldl -lm -lgomp
CFLAGS += -fopenmp
endif

ifneq (,$(findstring icc,$(CC)))
CFLAGS += -Wall -D_FORTIFY_SOURCE=0 -ipo -xHost -axAVX -openmp 
#ICCFLAGS+= -vec-report=5
#ICCFLAGS+= -fast
LDLIBS += -pthread -openmp -openmp-link static
# OPENMPSWITCH = -openmp 
ifeq ($(shell uname -s),Linux)
  CFLAGS += -openmp-threadprivate=compat
endif
endif

# unless otherwise determined assume Linux
SYSOBJS=getusage.o stark_numa.o 
ifeq ($(shell uname -s),Darwin)
  # For compiling on install machine
  SYSOBJS=getusage_darwin.o
else
	LDLIBS += -lrt
endif

#PBS status updates
-include /etc/pbs.conf
ifeq ($(PBS_EXEC),) 
	PBS_INCLUDE=
	PBS_LIBS=
else
	#include PBS status code
	CFLAGS += -I$(PBS_EXEC)/include -DPBS_STATUS
	LDLIBS += -L$(PBS_EXEC)/lib -lpbs
endif

#NUMA awaraness
ifeq ($(wildcard /usr/include/numa.h),) 
else 
    CFLAGS += -DUSE_LIBNUMA
	LDLIBS += -lnuma
endif


UBIOBJS=starkubi.o terminalubi.o test.o
STARKOBJS=stark.o stark_t.o stark_serialize_binary.o threadpool.o stark_status.o stark_increment_cache.o fastx.o stark_R.o stark_insert.o stark_assemble_hierarchical.o stark_phase1.2.o $(SYSOBJS)
TERMINALOBJS=terminal.o terminalubi.o test.o
BINARIES=interactive stark

REFCOMPAREOBJS=stark.o stark_t.o stark_serialize_binary.o threadpool.o stark_status.o stark_increment_cache.o fastx.o stark_R.o stark_insert.o $(SYSOBJS)

ALLOBJS=$(UBIOBJS) $(STARKOBJS) $(TERMINALOBJS) main.o 

all: interactive

stark: CFLAGS += -O3 -DTHREADED_IMPORT -DAMO_CACHE -DSTARK_COVERAGE_CACHE_TOKEN -DHIERARCHIAL_ASSEMBLY
stark: main.o $(STARKOBJS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LOADLIBES) $(LDLIBS)

profiling: PBS_INCLUDE=
profiling: PBS_LIBS=
profiling: CFLAGS += -O2 -g -DTHREADED_IMPORT -DAMO_CACHE -DDEBUG -DSTARK_COVERAGE_CACHE_TOKEN -DHIERARCHIAL_ASSEMBLY -fno-inline-functions -fno-inline -DNOCHECKS
profiling: main.o $(STARKOBJS)
		$(CC) $(LDFLAGS) -o stark $^ $(LDLIBS)

debug: CFLAGS += -O0 -g -DTHREADED_IMPORT -DAMO_CACHE -DDEBUG -DSTARK_COVERAGE_CACHE_TOKEN -DHIERARCHIAL_ASSEMBLY
debug: LDFLAGS += -g
debug: main.o $(STARKOBJS)
	$(CC) $(LDFLAGS) -o stark $^ $(LDLIBS)

# stark_uv: DEBUGFLAGS = -O3 -Wall -DTHREADED_IMPORT -DAMO_CACHE -DUV_SYSTEM -DSTARK_COVERAGE_CACHE_TOKEN
# stark_uv: main.o $(STARKOBJS) uv_system.o
# 	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS)

interactive: CFLAGS += -O0 -g -DINTERACTIVE -DUBIGRAPH -DDEBUG -DSTARK_ALLOC_MALLOC -DHIERARCHIAL_ASSEMBLY
interactive: LDLIBS += -lreadline $(UBILIBS)
interactive: LDFLAGS += -g
interactive: $(STARKOBJS) test.o terminal.o threadpool.o UbigraphThreaded.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS) 
	dsymutil interactive
	
interactive_noubi: UBIGRAPH =-DINTERACTIVE -DDEBUG -DSTARK_ALLOC_MALLOC -DHIERARCHIAL_ASSEMBLY
interactive_noubi: LDFLAGS += -g
interactive_noubi: $(STARKOBJS) test.o terminal.o threadpool.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS) 

stark_refcompare: CFLAGS += -O3 -DTHREADED_IMPORT -DAMO_CACHE -DSTARK_COVERAGE_CACHE_TOKEN
stark_refcompare: stark_refcompare.o $(REFCOMPAREOBJS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS) 
	
stark_refcompare_debug: CFLAGS += -O0 -g -DTHREADED_IMPORT
stark_refcompare_debug: LDFLAGS += -g
stark_refcompare_debug: stark_refcompare.o $(REFCOMPAREOBJS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LDLIBS) 

stark_commander: UBIGRAPH =-DINTERACTIVE
stark_commander: $(STARKOBJS) test.o terminal.o 
		$(CC) $(DEBUGFLAGS) -g -o $@ test.o $(STARKOBJS) terminal.o $(LIBS)
		dsymutil stark_commander

# starkubi.o: stark.c
# 	$(CC) $(CFLAGS) -c -o $@ $<
# terminalubi.o: terminal.c
# 	$(CC) $(CFLAGS) -c -o $@ $<

#starkubi.o: starkubi.c
#	$(CC) $(CFLAGS) -c $@ $^ -fnested-functions

# libs/resource.o: xmlrpc/resource.c
# 	echo disabled

libs/md5.o: md5.c
	mkdir -p libs
	$(CC) $(CFLAGS) -c -o $@ $<

threadpool.o: threadpool.h
$(UBIOBJS): UbigraphThreaded.h 
$(UBIOBJS) $(STARKOBJS) main.o: stark.h stark_t.h list.h uv_system.h stark_increment_cache.h stark_alloc_mmap.h stark_alloc_malloc.h threadpool.h stark_R.h stark_status.h stark_numa.h stark_navigate.h stark_assemble.h stark_assemble_hierarchical.h

$(TERMINALOBJS): terminal.h

$(ALLBOJS): Makefile main.h
# $(BINARIES): Makefile

sqlite3/sqlite3.o: sqlite3.c sqlite3.h
	$(CC) -O3 -Wall -c -o $@ $<

clean:
	rm -rf *.o libs/md5.o sqlite3/sqlite3.o stark_status_xsl.h stark_contour_log_stats_plot.h


distributor_test.o: distributor.h

dt: distributor_test.o
	$(CC) $(CFLAGS) -o distributor_test $< -lpthread
	dsymutil distributor_test
	
stark_status_xsl.h: stark_status.xsl
	xxd -i $< $@

stark_contour_log_stats_plot.h: stark_contour_log_stats_plot.R
	xxd -i $< $@

stark_status.o: stark_status_xsl.h
stark_R.o: stark_contour_log_stats_plot.h

main.o stark.o stark_refcompare.o stark_assemble_hierarchical.o: CFLAGS += $(OPENMPSWITCH)

fastx_rmdup: CFLAGS += -O3 -D_GNU_SOURCE
fastx_rmdup: fastx_rmdup.o fastx.o
	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $^ $(LDLIBS) 
