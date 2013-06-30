#ifndef FIND_SCORE_BLOCK_CUH
#define FIND_SCORE_BLOCK_CUH

#include "cuda.h"
#include "structs.cuh"
#include "cub/cub/cub.cuh"

#define INIT_VAL 0
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
						   seg_val *aux,
						   int aux_size) {
   	int id = blockIdx.x * config.block_size + threadIdx.x;
	int start = id * config.thread_chunk;
	int limit = min(start + config.thread_chunk, row_len);
	int init_seq_idx = placement(seq_last_idx, seq_count, start);
	int it = start;
	int seq_idx = init_seq_idx;

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
		(output + it)->N = INIT_VAL;
		(output + it)->H = INIT_VAL;
		(output + it)->V = INIT_VAL;
		++it;
	}
	
	for (; it < limit; ++it) {
		if ((*(horizontal + it - 1)) == config.reset) {
			(output + it)->N = INIT_VAL;
			(output + it)->H = INIT_VAL;
			(output + it)->V = INIT_VAL;
		} else {
			(output + it)->N = max(0, (*(horizontal + it - 1) == vertical) + 
										max((input + it - 1)->N,
											max((input + it - 1)->H, 
												(input + it - 1)->V
											   )
									   	   )
								   );
			(output + it)->V = max(max((input + it)->N - penalty->open, 0),
									max((input + it)->H - penalty->open, 
										(input + it)->V - penalty->extension
									   )
								   );
		}
		if ((output + it)->N > max_scr.val)
			max_scr.val = (output + it)->N;
	}
	
	//block reduce faza za N i V
	__syncthreads();
	typedef cub::BlockReduce<seg_val, THREADS> BlockReduce;
	__shared__ typename BlockReduce::SmemStorage reduce_storage;
	max_scr = BlockReduce::Reduce(reduce_storage, max_scr, maxop<seg_val>());

	//------------------device reduce faza za N i V------------------
	
	__syncthreads();
	int proc_blocks = config.grid_size;
	int block_id = blockIdx.x * config.block_size;
	while (proc_blocks > 1) {
		if (blockIdx.x < proc_blocks && threadIdx.x == 0) {
			*(aux + blockIdx.x) = max_scr;
		}
		__syncthreads();
		proc_blocks = ceildiv(proc_blocks, config.block_size);
		if (block_id < proc_blocks) {
			max_scr = (id < proc_blocks ? *(aux + id) : identity);
			max_scr = BlockReduce::Reduce(reduce_storage, max_scr, maxop<seg_val>());
		}
	}
	
	//---------------------------------------------------------------
	
	if (id == 0 && *total_max < max_scr.val)
		*total_max = max_scr.val;

	//H
	__syncthreads();
	seq_idx = init_seq_idx;
	for (it = start + (start == 0 ? 1 : 0); it < limit; ++it) {
		if ((*(horizontal + it - 1)) == config.reset) {
			++seq_idx;
			continue;
		}
		(output + it)->H = max(0, 
								max((output + it - 1)->N, 
									(output + it - 1)->V
								   ) 
								- penalty->open + it * penalty->extension);
		if ((output + it)->H > h_max.val) {
			h_max.val = (output + it)->H;
			h_max.seg_idx = seq_idx;
		}
	}

	typedef cub::BlockScan<seg_val, THREADS> BlockScan;
	__shared__ typename BlockScan::SmemStorage scan_storage;
	
	//------------------device scan faza za H------------------	
	//if (config.grid_size > 1) {
	//	seg_val block_max = BlockReduce::Reduce(reduce_storage, h_max, maxop<seg_val>());
	//	
	//	if (threadIdx.x == 0) {
	//		(aux + blockIdx.x)->val = block_max.val;
	//		(aux + blockIdx.x)->seg_idx = block_max.seg_idx;
	//	}
	//
	//	__syncthreads();
	//	if (blockIdx.x == 0) {
	//		int proc_chunk = ceildiv(config.grid_size, config.block_size);
	//		int proc_start = threadIdx.x * proc_chunk;
	//		int proc_limit = min(proc_start + proc_chunk, config.grid_size);
	//		seg_val thread_max = identity;
	//		
	//		for (int i = proc_start; i < proc_limit; ++i)
	//			thread_max = mymax<seg_val>(thread_max, *(aux + i));
	//		BlockScan::ExclusiveScan(scan_storage, thread_max, thread_max, identity, maxop<seg_val>());
	//		for (int i = proc_limit - 1; i > proc_start; --i)
	//			*(aux + i) = *(aux + i - 1);
	//		*(aux + proc_start) = thread_max;
	//		if (threadIdx.x == 0)
	//			*flag = 1;
	//	}
	//	
	//	while(!(*flag));
	//	if (threadIdx.x == 0)
	//		h_max = *(aux + blockIdx.x);
	//}
	
	//---------------------------------------------------------
	
	//------------------device scan faza za H v2------------------
	seg_val block_max = BlockReduce::Reduce(reduce_storage, h_max, maxop<seg_val>());
	
	int step = 1;
	proc_blocks = config.grid_size;
	__syncthreads();
	while (proc_blocks > 1) {		
		if (blockIdx.x < proc_blocks && threadIdx.x == 0)
			(*(aux + (1 + blockIdx.x) * step - 1)) = block_max;
		proc_blocks = ceildiv(proc_blocks, config.block_size);
		__syncthreads();
		if (blockIdx.x < proc_blocks) {
			if ((1 + id) * step - 1 < config.grid_size)
				block_max = (*(aux + (1 + id) * step - 1));
			else
				block_max = identity;
			block_max = BlockReduce::Reduce(reduce_storage, block_max, maxop<seg_val>());
		}
		step *= config.block_size;
	}
	//odjebat ovo sta je vece od max grid
	//spustit scan
	step /= config.block_size;
	if (id == 0)
		(*(aux + (1 + blockIdx.x) * step - 1)) = block_max;
	
	while (step > 1) {
		if (blockIdx.x < proc_blocks) {
			if ((1 + id) * step - 1 < config.grid_size)
				block_max = (*(aux + (1 + id) * step - 1));
			else
				block_max = identity;
			BlockScan::ExclusiveScan(scan_storage, block_max, block_max, identity, maxop<seg_val>());
			if ((1 + id) * step - 1 < config.grid_size)
				(*(aux + (1 + id) * step - 1)) = block_max;
		}
	}
	//------------------device scan faza za H--------------------
	
	__syncthreads();
	//if (config.grid_size == 1) {
		BlockScan::ExclusiveScan(scan_storage, h_max, h_max, identity, maxop<seg_val>());
		(output + start)->H = max(0, h_max.val - it * penalty->extension);
		for (it = start + 1; it < limit; ++it) {
			if (*(horizontal + it - 1) == config.reset)
				(output + it)->H = INIT_VAL;
			else
				(output + it)->H = max(0, (output + it - 1)->H - it * penalty->extension);
		}
	//}
}

/*
__global__ void interblock_scan(seg_val *aux,
						   configuration config) {
	int proc_chunk = ceildiv(config.grid_size, config.block_size);
	int proc_start = threadIdx.x * proc_chunk;
	int proc_limit = min(proc_start + proc_chunk, config.grid_size);
	
	seg_val identity;
	identity.val = 0;
	identity.seg_idx = 0;
	
	seg_val thread_max = identity;
	
	for (int i = proc_start; i < proc_limit; ++i)
		thread_max = mymax<seg_val>(thread_max, *(aux + i));
	typedef cub::BlockScan<seg_val, THREADS> BlockScan;
	__shared__ typename BlockScan::SmemStorage scan_storage;
	BlockScan::ExclusiveScan(scan_storage, thread_max, thread_max, identity, maxop<seg_val>());
	for (int i = proc_limit - 1; i > proc_start; --i)
		*(aux + i) = *(aux + i - 1);
	*(aux + proc_start) = thread_max;
}

__global__ void update_H(data *output,
						 seg_val *aux,
						 gap *penalty,
						 int row_len,
						 configuration config,
						 char *horizontal,
						 int seq_count,
						 int *seq_last_idx) {
	int id = blockIdx.x * config.block_size + threadIdx.x;
	int start = id * config.thread_chunk;
	int limit = min(start + config.thread_chunk, row_len);
	seg_val h_max;
	h_max.val = INIT_VAL;
	h_max.seg_idx = 0;
	
	seg_val identity;
	identity.val = INIT_VAL;
	identity.seg_idx = 0;	
	
	int seq_idx = placement(seq_last_idx, seq_count, start);
	
	for (int it = start + (start == 0 ? 1 : 0); it < limit; ++it) {
		if ((*(horizontal + it - 1)) == config.reset) {
			++seq_idx;
			continue;
		}
		(output + it)->H = max(0, 
								max((output + it - 1)->N, 
									(output + it - 1)->V
								   ) 
								- penalty->open + it * penalty->extension);
		if ((output + it)->H > h_max.val) {
			h_max.val = (output + it)->H;
			h_max.seg_idx = seq_idx;
		}
	}
	
	typedef cub::BlockScan<seg_val, THREADS> BlockScan;
	__shared__ typename BlockScan::SmemStorage scan_storage;
	BlockScan::ExclusiveScan(scan_storage, h_max, h_max, identity, maxop<seg_val>());
	
	if (threadIdx.x == 0)
		h_max = *(aux + blockIdx.x);
	__syncthreads();
	BlockScan::ExclusiveScan(scan_storage, h_max, h_max, identity, maxop<seg_val>());
	
	(output + start)->H = max(0, h_max.val - start * penalty->extension);
	for (int it = start + 1; it < limit; ++it) {
		if (*(horizontal + it - 1) == config.reset)
			(output + it)->H = INIT_VAL;
		else
			(output + it)->H = max(0, (output + it - 1)->H - it * penalty->extension);
	}
}
*/
#endif
