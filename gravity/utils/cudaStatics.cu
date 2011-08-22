#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#ifdef UNIX
#include <stdint.h>
#include <unistd.h>
#endif
#include "mex.h"

// CUDA
#include "cuda.h"
#include "cuda_runtime.h"
#include "cublas.h"
#include "GPUmat.hh"

// static paramaters
static int init = 0;
static GPUmat *gm;

#include "cudaCommon.h"

__global__ void cukern_applySpecial_fade(double *arr, double *linAddrs, double *consts, double *fadeCoeff, int nSpecials);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (init == 0) {
    // Initialize function
    // mexLock();
    // load GPUmat
    gm = gmGetGPUmat();
    init = 1;
  }

  int numelarray, numelfixies;

  if( (nlhs != 0) || (nrhs != 5)) { mexErrMsgTxt("cudaStatics operator is cudaStatics(array, linearIndices, constants, fadeCoeffs, blockdim)"); }
  double **array = getGPUSourcePointers(prhs, 1, &numelarray,  0, gm);
  double **fixies= getGPUSourcePointers(prhs, 3, &numelfixies, 1, gm);

  int blockdim = (int)*mxGetPr(prhs[4]);

//printf("narray elements: %i; nfixies: %i; bsize: %i\n", numelarray, numelfixies, blockdim); fflush(stdout);

  dim3 griddim; griddim.x = numelfixies / blockdim + 1;
  if(griddim.x > 32768) {
    griddim.x = 32768;
    griddim.y = numelfixies/(blockdim*griddim.x) + 1;
    }
  cukern_applySpecial_fade<<<griddim, blockdim>>>(array[0], fixies[0], fixies[1], fixies[2], numelfixies);
}

__global__ void cukern_applySpecial_fade(double *arr, double *linAddrs, double *consts, double *fadeCoeff, int nSpecials)
{
int myAddr = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x*blockIdx.y);
if(myAddr >= nSpecials) return;

double f = fadeCoeff[myAddr];
int xaddr = (int)linAddrs[myAddr];

arr[xaddr] = f*consts[myAddr] + (1.0-f)*arr[xaddr];

}


