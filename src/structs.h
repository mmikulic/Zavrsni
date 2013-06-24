#ifndef MY_STRUCTS_H
#define MY_STRUCTS_H

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

typedef struct {
	int val;
	int idx;
} vector_element;

typedef struct {
	vector_element el;
	int row;
} matrix_element;


#endif
