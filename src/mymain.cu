#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "cuda.h"
#include "structs.h"
#include "device_scan.cuh"
#include "find_score.cuh"

char *get_protein(char *filename, int *s, char reset) {
	FILE *f = fopen(filename, "r");
	char c;
	char *protein;
	int ignore;
	int it, size;
	
	rewind(f);
	size = 0;
	while(fscanf(f, "%c", &c) != EOF) {
		printf("size = %d, reading c = '%c'\n", size, c);
		if (c == '>') {
			ignore = 1;
			if (size > 0)
				++size;
		}
		else if (ignore && c == '\n')
			ignore = 0;
		else if (!ignore && c >= 'A' && c <= 'Z')
			++size;
	}
	*s += size;
	printf("total size: %d\n\n", *s);
	fflush(stdout);
	
	protein = (char *)malloc((size + 1) * sizeof(char));
	it = 0;
	rewind(f);
	while(fscanf(f, "%c", &c) != EOF) {
		if (c == '>') {
			ignore = 1;
			if (it > 0) {
				*(protein + it) = reset;
				++it;
			}
		}
		else if (ignore && c == '\n')
			ignore = 0;
		else if (!ignore && c >= 'A' && c <= 'Z') {
			if (it >= size) {
				printf("char * size not enough!\n");
				fflush(stdout);
				return NULL;
			}
			*(protein + it) = c;
			++it;
		}
	}
	*(protein + size) = '\0';
	fclose(f);
	
	return protein;
}

void init(data *mat, int val, int size) {
	data->N = (int *)malloc(size * sizeof(int));
	data->H = (int *)malloc(size * sizeof(int));
	data->V = (int *)malloc(size * sizeof(int));
	
	for (int i = 0; i < size; ++i) {
		*(mat->N + i) = 0;
		*(mat->H + i) = 0;
		*(mat->V + i) = 0;
	}
}

int main(int argc, char **argv) {
	//variableinit
	
	printf("Welcome!\n");
	printf("This is a protein alignment software. The expected input is \n");
	printf("a list of files that contain protein representations.\n");
	
	if (argc < 3) {
		printf("I expect to receive at least two filenames as arguments.\n");
		exit(-1);
	}
	
	
//	printf("reset character?\n> ");
//	scanf(" %c", &(config.reset));
	
	int v_len = 0;	
	char *vertical = get_protein(argv[1], &v_len, config.reset);
	
	int h_len = 0;
	char *horizontal = get_protein(argv[2], &h_len, config.reset);
	for (int i = 3; i < argc; ++i) {
		strcat(horizontal, &(config.reset));
		++h_len;
		strcat(horizontal, get_protein(argv[i], &h_len, config.reset));
	}
	
	int mat_len = h_len + 1;
	configuration config;
	config.reset = '#';
	config.thread_chunk = 256;
	config.block_size = min(512, (mat_len + config.thread_chunk - 1) / 
														config.thread_chunk);
	config.grid_size = (mat_len + config.block_size * config.thread_chunk - 1) / 
									(config.block_size * config.thread_chunk);
	
	data matRow[2];//TODO: set up
	data devMatRow[2];//TODO: set up
	
	int *auxiliary = (int *)malloc(config.grid_size);
	int devAux;//set up
	
	return 0;
}
