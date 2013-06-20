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
					       configuration config, 
					       int *total_max) {
	int start = (blockIdx.x * config.block_size + threadIdx.x) * config.thread_chunk;
	int limit = min(start + config.thread_chunk, h_length);
	
	int it = start;
	vector_element max_scr;
	max_scr.val = INIT_VAL;
	max_scr.idx = -1;
	//check
	if (start == 0) {
		*(output->N + start) = INIT_VAL;
		*(output->H + start) = INIT_VAL;//this is a bit out of place, but it should be faster this way
		*(output->V + start) = INIT_VAL;
	}
	
	for (it = start + 1; it < limit; ++it) {
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
				//check if allowed
			*(output->V + it) = max(max(*(input->N + it) - penalty->open, 0),
									max(*(input->H + it) - penalty->open, 
										*(input->V + it) - penalty->extension
									   )
								   );
		}
		if (*(output->N + it) > max_scr.val) {
			max_scr.val = *(output->N + it);
			max_scr.idx = it;
		}
	}
	
	//device reduce faza za N i V
	
	max_scr.val = INIT_VAL;
	max_scr.idx = -1;	
	for (it = start + 1; it < limit; ++it) {
		if (*(horizontal + it - 1) == config.reset)
			continue;
		else {
			*(output->H + it) = max(0, 
									max(*(output->N + it - 1), 
										*(output->V + it - 1)
									   ) 
									- penalty->open + it * penalty->extension);
			if (*(output->N + it) > max_scr.val) {
				max_scr.val = *(output->N + it);
				max_scr.idx = it;
			}
		}
	}
	
	//devicescan faza za H

	*(output->H + start) = max(0, max_scr.val);
	for (it = start + 1; it < limit; ++it) {
		if (*(horizontal + it - 1) == config.reset)
			*(output->H + it) = INIT_VAL;
		else
			*(output->H + it) = max(0, *(output->H + it - 1) - it * penalty->extension);
	}
}

#endif
