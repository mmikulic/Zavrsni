#ifndef DEVICE_SCAN_CUH
#define DEVICE_SCAN_CUH

#include "structs.h"
#include "utils.cuh"
#include "cub/cub.cuh"

__device__ void reduce(int *array, int limit, int offset, int pow) {
	int step = 1 << pow;
	for (int i = closestPow2(offset) - 1; i < limit; i += step) {
		*(array + i) = max(*(array + i - step), *(array + i));
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
