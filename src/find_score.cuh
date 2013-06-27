#ifndef FIND_SCORE_BLOCK_CUH
#define FIND_SCORE_BLOCK_CUH

#include "cuda.h"
#include "structs.h"
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

	//int threads = config.block_size;
	seg_val max_scr;
	seg_val h_max;
	max_scr.val = INIT_VAL;
	max_scr.idx = -1;
	h_max.val = INIT_VAL;
	h_max.seq = -1;
	
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
	max_scr.val = BlockReduce::Reduce(reduce_storage, max_scr, maxop<seg_val>());
	//dodat reduce izmedu blokova
	if (id == 0 && *total_max < max_scr.val)
		*total_max = max_scr.val;
	
	__syncthreads();	
	//H
	for (it = start + (it == 0 ? 1 : 0); it < limit; ++it) {
		(*(output + it)).H = max(0, 
								max((*(output + it - 1)).N, 
									(*(output + it - 1)).V
								   ) 
								- penalty->open + it * penalty->extension
							   );//check
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
	seg_val identity;
	identity.val = 0;
	identity.seg_idx = 0;
	typedef cub::BlockScan<seg_val, THREADS> BlockScan;
	__shared__ typename BlockScan::SmemStorage scan_storage;
	BlockScan::ExclusiveScan(scan_storage, h_max, h_max, identity, maxop<seg_val>());
	
	
	typedef cub::BlockScan<seg_val, THREADS> BlockScan;
	__shared__ typename BlockScan::SmemStorage scan_storage;
	BlockScan::ExclusiveScan(scan_storage, h_max, h_max, identity, maxop<seg_val>());

	(*(output + start)).H = max(0, max_scr.val - it * penalty->extension);
	for (it = start + 1; it < limit; ++it) {
		if (*(horizontal + it - 1) == config.reset)
			(*(output + it)).H = INIT_VAL;
		else
			(*(output + it)).H = max(0, (*(output + it - 1)).H - it * penalty->extension);
	}
}

#endif
