#ifndef FIND_SCORE_BLOCK_CUH
#define FIND_SCORE_BLOCK_CUH

#include "cuda.h"
#include "structs.h"
#include "cub/cub/cub.cuh"

#define INIT_VAL 0
#define NEG_INF -(1 << 15)

__global__ void find_score (gap *penalty, 
					   	   char *horizontal, 
					   	   char vertical, 
					   	   int row_len,  
					   	   data *input, 
					       data *output,
					       configuration config, 
					       int *total_max) {
   	int id = blockIdx.x * config.block_size + threadIdx.x;
	int start = id * config.thread_chunk;
	int limit = min(start + config.thread_chunk, row_len);
	
	int it = start;
	int max_scr = INIT_VAL;
	int h_max = INIT_VAL;
	//check
	if (it == 0) {
		*(output + it)->N = INIT_VAL;
		*(output + it)->H = INIT_VAL;
		*(output + it)->V = INIT_VAL;
		++it;
	}
	
	while (it < limit) {
		if (*(horizontal + it - 1) == config.reset) {
			*(output + it)->N = INIT_VAL;
			*(output + it)->H = INIT_VAL;
			*(output + it)->V = INIT_VAL;
		} else {
			*(output + it)->N = max(0, (*(horizontal + it - 1) == vertical) + 
										max(*(input + it - 1)->N,
											max(*(input + it - 1)->H, 
												*(input + it - 1)->V
											   )
									   	   )
								   );
			*(output + it)->H = max(0, 
									max(*(output + it - 1)->N, 
										*(output + it - 1)->V
									   ) 
									- penalty->open + it * penalty->extension
								   );//check
			*(output + it)->V = max(max(*(input + it)->N - penalty->open, 0),
									max(*(input + it)->H - penalty->open, 
										*(input + it)->V - penalty->extension
									   )
								   );
		}
		if (*(output + it)->N > max_scr)
			max_scr = *(output + it)->N;
		if (*(output + it)->H > h_max)
			h_max = *(output + it)->H;
		++it;
	}
	
	//device reduce faza za N i V
	__syncthreads();
	typedef cub::BlockReduce<int, config.thread_chunk> BlockReduce;
	__shared__ typename BlockReduce::SmemStorage smem_storage;
	max_scr = BlockReduce::Reduce(smem_storage, max_scr, maxop<int>());
	if (id == 0 && *total_max < max_scr)
		*total_max = max_scr;
	
	//devicescan faza za H
	__syncthreads();
	int identity = 0;
	typedef cub::BlockScan<int, config.thread_chunk> BlockScan;
	BlockScan::ExclusiveScan(smem_storage, h_max, h_max, identity, maxop<int>());

	(*(output + start))->H = max(0, max_scr - it * penalty->extension);
	for (it = start + 1; it < limit; ++it) {
		if (*(horizontal + it - 1) == config.reset)
			(*(output + it))->H = INIT_VAL;
		else
			(*(output + it))->H = max(0, *(output + it - 1)->H - it * penalty->extension);
	}
}

#endif
