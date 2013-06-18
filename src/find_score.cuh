#ifndef FIND_SCORE_CUH
#define FIND_SCORE_CUH

#include "cuda.h"
#include "structs.h"

#define NEG_INF -(1 << 15)

__global__ find_row_score (gap *penalty, 
					   	   char *horizontal, 
					   	   char vertical, 
					   	   int h_length,  
					   	   data *input, 
					       data *output,
					       configuration config, 
					       int *total_max) {
	int start = (blockIdx.x * config.block_size + threadIdx.x) * config.thread_chunk;
	int limit = min(start + config.thread_chunk, h_length);
	
	int it = start;
	int max_scr = NEG_INF;
	//check
	if (start == 0) {
		*(output->N + start) = NEG_INF;
		*(output->H + start) = NEG_INF;//this is a bit out of place, but it should be faster this way
		*(output->V + start) = NEG_INF;
	}
	
	for (it = start + 1; it < limit; ++it) {
		if (*(horizontal + it - 1) == config.reset) {
			*(output->N + it) = NEG_INF;
			*(output->H + it) = NEG_INF; //this is a bit out of place, but it should be faster this way
			*(output->V + it) = NEG_INF;
		} else {
			*(output->N + it) = (*(horizontal + it - 1) == vertical) + 
				max(*(input->N + it - 1), max(*(input->H + it - 1), *(input->V + it - 1)));
				//check if allowed
			*(output->V + it) = max(*(input->N + it) - penalty->open, 
				max(*(input->H + it) - penalty->open, *(input->V + it) - penalty->extension));
		}
		max_scr = max(max_scr, *(output->N + it));
	}
	
	//device reduce faza za N i V
	
	for (it = start + 1; it < limit; ++it) {
		if (*(horizontal + it - 1) == config.reset)
			continue;
		else {
			*(output->H + it) = max(*(output->N + it - 1), *(output->V + it - 1)) - penalty->open + it * penalty->extension;
			max_scr = max(max_scr, *(output->H + it));
		}
	}
	
	//devicescan faza za H
	
	*(output->H + start) = max_scr;
	for (it = start + 1; it < limit; ++it) {
		if (*(horizontal + it - 1) == config.reset)
			*(output->H + it) = NEG_INF;
		else
			//check
	}
}

#endif
