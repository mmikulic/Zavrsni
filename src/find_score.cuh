#ifndef FIND_SCORE_CUH
#define FIND_SCORE_CUH

#include "cuda.h"
#include "structs.h"
#include "cub/cub.cuh"

#define INIT_VAL 0
#define NEG_INF -(1 << 15)

__global__ find_score (gap *penalty, 
					   	   char *horizontal, 
					   	   char vertical, 
					   	   int h_length,  
					   	   data *input, 
					       data *output,
					       int *aux,
					       int *aux_size,
					       configuration config, 
					       int *total_max) {
   	int id = blockIdx.x * config.block_size + threadIdx.x;
	int start = id * config.thread_chunk;
	int limit = min(start + config.thread_chunk, h_length);
	
	int it = start;
	int max_scr = INIT_VAL;
	//check
	if (it == 0) {
		*(output->N + it) = INIT_VAL;
		*(output->H + it) = INIT_VAL;//this is a bit out of place, but it should be faster this way
		*(output->V + it) = INIT_VAL;
		++it;
	}
	
	while (it < limit) {
		if (*(horizontal + it - 1) == config.reset) {
			*(output->N + it) = INIT_VAL;
			*(output->H + it) = INIT_VAL; //this is a bit out of place, but it should be faster this way
			*(output->V + it) = INIT_VAL;
		} else {
			*(output->N + it) = max(0, (*(horizontal + it - 1) == vertical) + 
										max(*(input->N + it - 1),
											max(*(input->H + it - 1), 
												*(input->V + it - 1)
											   )
									   	   )
								   );
			*(output->V + it) = max(max(*(input->N + it) - penalty->open, 0),
									max(*(input->H + it) - penalty->open, 
										*(input->V + it) - penalty->extension
									   )
								   );
		}
		if (*(output->N + it) > max_scr) {
			max_scr = *(output->N + it);
		}
		++it;
	}
	
	//device reduce faza za N i V
	__syncthreads();
	typedef cub::BlockReduce<int, config.thread_chunk> BlockReduce;
	__shared__ typename BlockReduce::SmemStorage smem_storage;
	max_scr = BlockReduce::Reduce(smem_storage, max_scr, MaxOperator<int>());
	int curr_size = *aux_size;
	while (curr_size > config.block_size) {
		if (threadIdx.x == 0)
			*(aux + blockIdx.x) = max_scr;
			curr_size = (curr_size + config.block_size - 1) / config.block_size;
		__syncthreads();
		if (id - threadIdx.x < curr_size) {
			max_scr = (id < curr_size ? *(aux + id), 0);
			max_scr = BlockReduce::Reduce(smem_storage, max_scr, maxop<int>());
		}
	}
	__syncthreads();
	if (id == 0 && *total_max < max_scr)
		*total_max = max_scr;	
	
	//H
	max_scr = INIT_VAL;
	it = start + (start == 0 ? 1 : 0);
	while (it < limit) {
		if (*(horizontal + it - 1) != config.reset) {
			*(output->H + it) = max(0, 
									max(*(output->N + it - 1), 
										*(output->V + it - 1)
									   ) 
									- penalty->open + it * penalty->extension
								   );
			if (*(output->N + it) > max_scr)
				max_scr = *(output->N + it);
		}
		++it;
	}
	
	//devicescan faza za H
	__syncthreads();
	int identity = 0;
	typedef cub::BlockScan<int, config.thread_chunk> BlockScan;
	BlockScan::ExclusiveScan(smem_storage, max_scr, max_scr, identity, maxop<int>());
	
	*(aux + id) = max_scr;
	curr_size = *aux_size;
	while(curr_size)
	
	

	max_scr = *(aux + id);
	*(output->H + start) = max(0, max_scr);
	for (it = start + 1; it < limit; ++it) {
		if (*(horizontal + it - 1) == config.reset)
			*(output->H + it) = INIT_VAL;
		else
			*(output->H + it) = max(0, *(output->H + it - 1) - it * penalty->extension);
	}
}

#endif
