#ifndef MY_STRUCTS_CUH
#define MY_STRUCTS_CUH

typedef struct {
    int open;
    int extension;
} gap;

typedef struct {
	int N;
	int H;
	int V;
} data;

typedef struct {
	char reset;
	int block_size;
	int thread_chunk;
	int grid_size;
} configuration;

struct seg_val {
	int val;
	int seg_idx;
	__device__ __host__ bool operator< (const seg_val &b) const {
		if (seg_idx != b.seg_idx) return seg_idx < b.seg_idx;
		return val < b.val;
	}
};


#endif
