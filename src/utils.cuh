#ifndef UTILS_CUH
#define UTILS_CUH

template< typename T > struct maxop {
	__device__ __host__ maxop() {}
	__device__ __host__ T operator() (const T &a, const T &b) const {
	    return a < b ? b : a;
	}
};
/*
__device__ __host__ int log2(int n) {
	int pow = 0;
	while((1 << pow) < n)
		++pow;
	return pow;
}
__device__ __host__ int closestPow2(int n) {
	int num = 1;
	while(num <= n)
		num = num << 1;
	return num;
}
*/
__device__ __host__ bool operator< (const seg_val &a, const seg_val &b) const {
	if (a.seg_idx != b.seg_idx) return a.seg_idx > b.seg_idx;
	return a.val < b.val;
}

void exitWithMsg(const char *msg, int exitCode) {
	printf("ERROR\n");
	printf("%s\n\n", msg);
	exit(exitCode);
}

void safeAPIcall(cudaError_t err, int line) {
	if(err != cudaSuccess) {
		printf("Error in line %d\n", line);
		exitWithMsg(cudaGetErrorString(err), -2);
	}
}

void cudaSetAndCopyToDevice(void **dest, void *src, int size, int line) {
	safeAPIcall(cudaMalloc(dest, size), line);
	//printf("\n%d: malloc done\n", line);
	safeAPIcall(cudaMemcpy(*dest, src, size, cudaMemcpyHostToDevice), line);
}

void cudaCopyToHostAndFree(void *dest, void *src, int size, int line) {
	safeAPIcall(cudaMemcpy(dest, src, size, cudaMemcpyDeviceToHost), line);
	safeAPIcall(cudaFree(src), line);
}

__device__ int placement(int *array, int size, int element) {
	int middle, bottom, top;
	if (size == 1)
		return 0;
	
	bottom = 0;
	top = size - 1;
	while(bottom < top) {
		middle = (top + bottom) >> 1;
		if (*(array + middle) < element)
			bottom = middle + 1;
		if (*(array + middle) >= element)
			top = middle;
	}
	return top;
}
#endif
