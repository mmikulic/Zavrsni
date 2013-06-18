#ifndef UTILS_CUH
#define UTISL_CUH

template< typename T > struct MaxOperator {
	__device__ __host__ MaxOperator() {}
	__device__ __host__ T operator() (const T &a, const T &b) const {
	    return a < b ? b : a;
	}
};

__device__ __host__ int log2(int n) {
	int num = 1;
	int pow = 0;
	while(num < n) {
		num *= 2;
		++pow;
	}
	return pow;
}

#endif
