#ifndef DEVICE_SCAN_CUH
#define DEVICE_SCAN_CUH

#include "structs.h"
#include "utils.cuh"
#include "cub/cub.cuh"

__device__ void deviceReduce(int *array, int thoffset, int thlimit, int arraysize) {
	int limit = log2(arraysize);	
	offset = closestPow2(offset);
	
	for (d = 1; d < limit; ++d) {
		
	}
	
}

__device__ void prescan(int *array, int limit, int offset, int pow) {
	for ()
}

__device__ void deviceScan(int *array, int limit, int offset, int pow) {
	reduce(array, limit, offset, pow);
	prescan(array, limit, offset, pow);
}



#endif
