#ifndef CONFIG_CUH
#define CONFIG_CUH

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

CUDAcard findBestDevice() {
	int numOfDevices, bestDeviceNumber;

	cudaDeviceProp bestDeviceProps;
	
	safeAPIcall(cudaGetDeviceCount(&numOfDevices), __LINE__);
	
	int maxCores = -1;

	for (int i = 0; i < numOfDevices; ++i) {
		cudaDeviceProp currentDeviceProps;
		safeAPIcall(cudaGetDeviceProperties(&currentDeviceProps, i), __LINE__);
			
		int deviceCores = _ConvertSMVer2Cores(currentDeviceProps.major,
				currentDeviceProps.minor) * currentDeviceProps.multiProcessorCount;

		if (maxCores < deviceCores) {
			maxCores = deviceCores;
			bestDeviceNumber = i;
			bestDeviceProps = currentDeviceProps;
		}
	}

	if(maxCores < 0 || numOfDevices < 1)
		exitWithMsg("No CUDA capable card detected.", -2);
	
	CUDAcard gpu;
	gpu.cardNumber = bestDeviceNumber;
	gpu.major = bestDeviceProps.major;
	gpu.minor = bestDeviceProps.minor;
	gpu.cardsInSystem = numOfDevices;
	gpu.maxThreadsPerBlock = bestDeviceProps.maxThreadsDim[0];
	gpu.SMs = bestDeviceProps.multiProcessorCount;
	gpu.cudaCores = maxCores;
	gpu.globalMem = bestDeviceProps.totalGlobalMem;
	strcpy(gpu.name, bestDeviceProps.name);

	return gpu;
}

#endif
