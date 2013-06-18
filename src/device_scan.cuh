#ifndef DEVICE_SCAN_CUH
#define DEVICE_SCAN_CUH

#include "structs.h"
#include "utils.cuh"
#include "cub/cub.cuh"

__device__ __host__

__device__ void deviceScan(int *in, int *out, int size) {
	int it;
	//pozvati devicereduce
	cub::DeviceReduce::Reduce(in, out, size, MaxOperator<int>());
	
	for (it = log2(size) - 1; it >= 0; --it) {
		
	}
	
}



#endif
