#ifndef UTILS_CUH
#define UTILS_CUH

template< typename T > struct maxop {
	__device__ __host__ maxop() {}
	__device__ __host__ T operator() (const T &a, const T &b) const {
	    return a < b ? b : a;
	}
};

//check
__device__ __host__ bool operator< (const vector_element &a, const vector_element &b) const {
	if (index != other.index) return index > other.index;
	return value < other.value;
}

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

#endif
