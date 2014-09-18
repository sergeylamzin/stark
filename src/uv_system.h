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


#ifdef UV_SYSTEM
#ifndef UV_SYSTEM_H
#define UV_SYSTEM_H

#include <stdint.h>
#include <uv/gru/gru.h>
#include <uv/gru/gru_instructions.h>


#ifdef UV_REV_1_WARS
#error UV REV 1 WARS unsupported

#endif
#define uv_gru_opdword(opcode, exopc, xtype, iaa0, iaa1, idef2, ima) \
 	   (1 << GRU_CB_ICMD_SHFT) | \
	   ((unsigned long)CBS_ACTIVE << GRU_ISTATUS_SHFT) | \
	   (idef2<< GRU_IDEF2_SHFT) | \
	   (iaa0 << GRU_CB_IAA0_SHFT) | \
	   (iaa1 << GRU_CB_IAA1_SHFT) | \
	   (ima << GRU_CB_IMA_SHFT) | \
	   (IMA_MAPPED << GRU_CB_IMA_SHFT) | \
	   (xtype << GRU_CB_XTYPE_SHFT) | \
	   (opcode << GRU_CB_OPC_SHFT) | \
	   (exopc << GRU_CB_EXOPC_SHFT)

#define GAMER_AMO_ADD_32_OPWORD __opdword(OP_GAMER, EOP_ER_ADD, XTYPE_W, IAA_RAM, 0, 0, 0)
#define GAMER_AMO_ADD_64_OPWORD __opdword(OP_GAMER, EOP_ER_ADD, XTYPE_DW, IAA_RAM, 0, 0, 0)

#ifndef UV_GRU_MAX_CBS_PER_THREAD
#define UV_GRU_MAX_CBS_PER_THREAD 32
#endif

// extern __thread gru_cookie_t gru_cookie[];
extern __thread gru_control_block_t *gru_control_block[];

#define gru_cb_is_idle(cb) (!(0x3000000 & ((volatile struct gru_instruction*)(cb))->tri0))
#define gru_cb_is_busy(cb) (0x3000000 & ((volatile struct gru_instruction*)(cb))->tri0)

#define rdtsc() ({uint64_t eax,edx; __asm__ ("rdtsc": "=a" (eax), "=d" (edx)); (eax | (edx << 32));})

extern __thread size_t times_called;
extern __thread size_t total_cycles_between;
extern __thread size_t times_waited;
extern __thread size_t total_cycles_waited;
extern __thread size_t last_call;
extern size_t gru_cbs_per_cpu;
extern __thread size_t uv_gru_num_context;

static inline void uv_gru_atomic_async(void * const addr, const uint64_t operand, const uint64_t opword) {
	times_called++;
	int num_context = 0;
	// while (gru_cb_is_busy(gru_control_block[num_context]))
	// 	num_context = (num_context + 1) % UV_GRU_CONTEXTS_PER_THREAD;
	
	uint64_t cycles_now = rdtsc();
	total_cycles_between += cycles_now - last_call;
	
	if (gru_cb_is_busy(gru_control_block[num_context])) {
		
		gru_wait_abort(gru_control_block[num_context]);
		uint64_t cycles_waited = rdtsc() - cycles_now;
		times_waited++;
		total_cycles_waited += cycles_waited;
	}
		
	
	((struct gru_instruction *)gru_control_block[num_context])->baddr0 = (long)addr;
	((struct gru_instruction *)gru_control_block[num_context])->op1_stride = operand;
	#ifdef UV_REV_1_WARS
		((struct gru_instruction *)gru_control_block[num_context])->nelem = 1;			// GRU 1.0 WAR
		extern void gru_implicit_abort_detected(struct gru_instruction *cb);
		struct gru_instruction_bits *bits = (struct gru_instruction_bits *)gru_control_block[num_context];
		if (bits->istatus >= CBS_ACTIVE)
			gru_implicit_abort_detected(gru_control_block[num_context]);
	#endif
	gru_start_instruction((struct gru_instruction *)gru_control_block[num_context], opword);
	/*
	__asm__ volatile (
		// "movq %1, %0\n\t" // this is so that any subsequent reads will have the most up to date version
		// the reason for this is that the processor will ensure proper ordering of temporal (default) read/writes. If you intend to read the status flag (which is part of the opword) before the non-temporal move (2 instructions down) leaves the write buffer then this will provide the correct answer.
		"sfence\n\t" // make sure all previous writes have been flushed out of the write buffer
		"movntiq %1, %0" // non-temporal write, this will invalidate the cacheline it is connected to as soon as it leaves the write buffer which will eventually happen.
	 	: "=m" (*(uint64_t*)(gru_control_block[0]))
		: "r" ((uint64_t)(opword)));
	*/
	last_call = rdtsc();
}

static inline void stark_uv_gru_increment_coverage_w_fallback(uint32_t * const addr, const uint32_t operand) {
	size_t i;
	for (i = 0; (i < uv_gru_num_context); i++) {
		if (gru_cb_is_idle(gru_control_block[i])) {
			times_called++;
			((struct gru_instruction *)gru_control_block[i])->baddr0 = (long)addr;
			((struct gru_instruction *)gru_control_block[i])->op1_stride = operand;
			#ifdef UV_REV_1_WARS
				((struct gru_instruction *)gru_control_block[i])->nelem = 1;			// GRU 1.0 WAR
				extern void gru_implicit_abort_detected(struct gru_instruction *cb);
				struct gru_instruction_bits *bits = (struct gru_instruction_bits *)gru_control_block[i];
				if (bits->istatus >= CBS_ACTIVE)
					gru_implicit_abort_detected(gru_control_block[i]);
			#endif
			gru_start_instruction((struct gru_instruction *)gru_control_block[i], GAMER_AMO_ADD_32_OPWORD);
			return;
		}
	}
	
	__sync_fetch_and_add(addr, operand);
}

// #define UV_SYSTEM_FLAG_DESTROY_GRU_CONTEXT_ON_PTHREAD_EXIT 1

int uv_init_atomic();
void uv_free_atomic();
int uv_init_global();

#endif
#endif
