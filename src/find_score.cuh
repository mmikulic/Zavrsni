#ifndef FIND_SCORE_BLOCK_CUH
#define FIND_SCORE_BLOCK_CUH

#include "cuda.h"
#include "structs.cuh"
#include "cub/cub/cub.cuh"

#define INIT_VAL 0
#define NEG_INF -(1 << 15)
#define THREADS 512

__global__ void find_score (gap *penalty, 
					   	   char *horizontal, 
					   	   char vertical, 
					   	   int row_len,  
					   	   data *input, 
					       data *output,
					       configuration config, 
					       int *total_max,
						   int *seq_last_idx,
						   int seq_count,
						   seg_val *aux) {
   	int id = blockIdx.x * config.block_size + threadIdx.x;
	int start = id * config.thread_chunk;
	int limit = min(start + config.thread_chunk, row_len);
	int seq_idx = placement(seq_last_idx, seq_count, start);
	int it = start;

	seg_val max_scr;
	max_scr.val = INIT_VAL;
	max_scr.seg_idx = 0;
	
	seg_val h_max;
	h_max.val = INIT_VAL;
	h_max.seg_idx = 0;
	
	seg_val identity;
	identity.val = 0;
	identity.seg_idx = 0;
	
	if (it == 0) {
		(*(output + it)).N = INIT_VAL;
		(*(output + it)).H = INIT_VAL;
		(*(output + it)).V = INIT_VAL;
		++it;
	}
	
	for (; it < limit; ++it) {
		if (*(horizontal + it - 1) == config.reset) {
			(*(output + it)).N = INIT_VAL;
			(*(output + it)).H = INIT_VAL;
			(*(output + it)).V = INIT_VAL;
			++seq_idx;
		} else {
			(*(output + it)).N = max(0, (*(horizontal + it - 1) == vertical) + 
										max((*(input + it - 1)).N,
											max((*(input + it - 1)).H, 
												(*(input + it - 1)).V
											   )
									   	   )
								   );
			(*(output + it)).V = max(max((*(input + it)).N - penalty->open, 0),
									max((*(input + it)).H - penalty->open, 
										(*(input + it)).V - penalty->extension
									   )
								   );
		}
		if ((*(output + it)).N > max_scr.val)
			max_scr.val = (*(output + it)).N;
	}
	
	//reduce faza za N i V
	__syncthreads();
	typedef cub::BlockReduce<seg_val, THREADS> BlockReduce;
	__shared__ typename BlockReduce::SmemStorage reduce_storage;
	max_scr = BlockReduce::Reduce(reduce_storage, max_scr, maxop<seg_val>());
	        //device reduce faza za N i V
	//__syncthreads();
	//max_scr = BlockReduce::Reduce(reduce_storage, max_scr, maxop<seg_val>());
	//int curr_size = config.grid_size;
	//while (curr_size > config.block_size) {
	//	if (threadIdx.x == 0)
	//			*(aux + blockIdx.x) = max_scr;
	//			curr_size = (curr_size + config.block_size - 1) / config.block_size;
	//	__syncthreads();
	//	if (id - threadIdx.x < curr_size) {
	//			max_scr = (id < curr_size ? *(aux + id) : identity);
	//			max_scr = BlockReduce::Reduce(reduce_storage, max_scr, maxop<seg_val>());
	//	}
	//}

	if (id == 0 && *total_max < max_scr.val)
		*total_max = max_scr.val;
	
	__syncthreads();
	for (it = start + (it == 0 ? 1 : 0); it < limit; ++it) {
		(*(output + it)).H = max(0, 
								max((*(output + it - 1)).N, 
									(*(output + it - 1)).V
								   ) 
								- penalty->open + it * penalty->extension);
		if ((*(output + it)).H > h_max.val) {
			h_max.val = (*(output + it)).H;
			h_max.seg_idx = seq_idx;
		}
	}
	
	//devicescan faza za H
		//scan unutar bloka
		//reduce unutar bloka
		//scan izmedu blokova
		//popunjavanje
	typedef cub::BlockScan<seg_val, THREADS> BlockScan;
	__shared__ typename BlockScan::SmemStorage scan_storage;
	
	//seg_val block_max = BlockReduce::Reduce(reduce_storage, h_max, maxop<seg_val>());
	//if (threadIdx.x == 0)
	//	*(aux + blockIdx.x) = block_max;
	//__syncthreads();
	//if (blockIdx.x == 0) {
	//	int proc_chunk = (config.grid_size + config.block_size - 1) / config.block_size;
	//	int proc_start = threadIdx.x * proc_chunk;
	//	int proc_limit = min(proc_start + proc_chunk, config.grid_size);
	//	seg_val proc_max;
	//	proc_max.val = 0;
	//	proc_max.seg_idx = 0;
	//	for (int i = proc_start; i < proc_limit; ++i) {
	//		if (proc_max < *(aux + i))
	//			proc_max = *(aux + i);
	//	}
	//	BlockScan::ExclusiveScan(scan_storage, proc_max, proc_max, identity, maxop<seg_val>());
	//	
	//	if (*(aux + proc_start) < proc_max)
	//		*(aux + proc_start) = proc_max;
	//	for (int i = proc_start + 1; i < proc_limit; ++i)
	//		*(aux + i) = mymax<seg_val>(*(aux + i), *(aux + i - 1));
	//}
	//
	//__syncthreads();
	//if (threadIdx.x == 0)
	//	h_max = mymax<seg_val>(*(aux + blockIdx.x), h_max);
	BlockScan::ExclusiveScan(scan_storage, h_max, h_max, identity, maxop<seg_val>());

	(*(output + start)).H = max(0, h_max.val - it * penalty->extension);
	for (it = start + 1; it < limit; ++it) {
		if (*(horizontal + it - 1) == config.reset)
			(*(output + it)).H = INIT_VAL;
		else
			(*(output + it)).H = max(0, (*(output + it - 1)).H - it * penalty->extension);
	}
}

#endif
