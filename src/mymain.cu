#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "cuda.h"
#include "utils.cuh"
#include "structs.h"
#include "find_score.cuh"

char *get_protein(char *filename, int *s, char reset, int *seqcount) {
	FILE *f = fopen(filename, "r");
	char c;
	char *protein;
	int ignore = 0;
	int it, size;
	
	rewind(f);
	size = 0;
	while(fscanf(f, "%c", &c) != EOF) {
		//printf("size = %d, reading c = '%c'\n", size, c);
		if (c == '>') {
			ignore = 1;
			++(*seqcount);
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

void init(data **mat, int val, int size) {
	*mat = (data *)malloc(size * sizeof(data));
	
	for (int i = 0; i < size; ++i) {
		(*(*mat + i)).N = val;
		(*(*mat + i)).H = val;
		(*(*mat + i)).V = val;
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
	configuration config;
	config.reset = '#';
	
	int v_len = 0;
	int v_seq_count = 0;
	char *vertical = get_protein(argv[1], &v_len, config.reset, &v_seq_count);
	if (v_seq_count > 1)
		exitWithMsg("More than one vertical sequence!", -1);
	
	int h_len = 0;
	int seq_count = 0;
	char *horizontal = get_protein(argv[2], &h_len, config.reset, seq_count);
	for (int i = 3; i < argc; ++i) {
		strcat(horizontal, &(config.reset));
		++h_len;
		strcat(horizontal, get_protein(argv[i], &h_len, config.reset, &seq_count));
	}
	printf("Vertical: %s\n", vertical);
	printf("Horizontal: %s\n", horizontal);
	
	int *seq_last_idx = (int *)malloc(seq_count * sizeof(int));
	int seq_it = 0;
	for (int i = 0; i < h_len; ++i) {
		if (*(horizontal + i) == '#') {			
			*(seq_last_idx + seq_it) = i - 1;
			++seq_it;
		}
	}
	
	int row_len = h_len + 1;
	config.grid_size = 1;
	config.block_size = min(512, (row_len + 9) / 10);
	config.thread_chunk = (row_len + config.block_size - 1) / config.block_size;
	
	data *matRow[2];
	init(&matRow[0], 0, row_len);
	init(&matRow[1], 0, row_len);
	
	int total_max = -1;
	gap penalty;
	penalty.open = 12;
	penalty.extension = 2;
	
	//dev variables
	gap *dev_penalty;
	int *dev_total_max;
	int *dev_sec_last_idx;
	char *dev_h
	data *devMatRow[2];
	seg_val *dev_auxiliary;
	
	safeAPIcall(cudaMalloc((void **)&dev_auxiliary, config.block_size * sizeof(seg_val)), __LINE__);
	cudaSetAndCopyToDevice((void **)&devMatRow[0], matRow[0], row_len * sizeof(data), __LINE__);
	cudaSetAndCopyToDevice((void **)&devMatRow[1], matRow[1], row_len * sizeof(data), __LINE__);
	cudaSetAndCopyToDevice((void **)&dev_sec_last_idx, seq_last_idx, seq_count * sizeof(int), __LINE__);
	cudaSetAndCopyToDevice((void **)&dev_h, horizontal, h_len * sizeof(char), __LINE__);
	cudaSetAndCopyToDevice((void **)&dev_penalty, &penalty, sizeof(gap), __LINE__);
	cudaSetAndCopyToDevice((void **)&dev_total_max, &total_max, sizeof(int), __LINE__);
	
	int curr = 0;	
	clock_t start_time = clock();
	for (int i = 0; i < v_len; ++i) {
		find_score<<<config.grid_size, config.block_size>>>(
			dev_penalty,
			dev_h,
			vertical[i],
			row_len,
			devMatRow[curr ^ 1],
			devMatRow[curr],
			config,
			dev_total_max,
			dev_seq_last_idx,
			seq_count,
			dev_auxiliary
		);
		curr ^= 1;
	}
	clock_t end_time = clock();
	printf("Time: %d ticks\n", end_time-start_time);
	cudaCopyToHostAndFree(&total_max, dev_total_max, sizeof(int), __LINE__);
	printf("max local alignment: %d\n", total_max);
	
	safeAPIcall(cudaFree(devMatRow[0]), __LINE__);
	safeAPIcall(cudaFree(devMatRow[1]), __LINE__);
	safeAPIcall(cudaFree(dev_penalty), __LINE__);
	safeAPIcall(cudaFree(dev_h), __LINE__);
	safeAPIcall(cudaFree(dev_seq_last_idx), __LINE__);
		
	return 0;
}
