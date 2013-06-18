#ifndef MY_STRUCTS_H
#define MY_STRUCTS_H

typedef struct {
    int open;
    int extension;
} gap;

typedef struct {
	int *N;
	int *H;
	int *V;
} data;

typedef struct {
	char reset;
	int block_size;
	int thread_chunk;
} configuration;


#endif
