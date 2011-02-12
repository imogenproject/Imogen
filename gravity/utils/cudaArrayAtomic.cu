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

__global__ void cukern_ArraySetMin(double *array, double min,    int n);
__global__ void cukern_ArraySetMax(double *array, double max,    int n);
__global__ void cukern_ArrayFixNaN(double *array, double fixval, int n);

#define BLOCKDIM 128

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // At least 2 arguments expected
  // Input and result
  if (nrhs!=3)
     mexErrMsgTxt("Wrong number of arguments");
  if (init == 0) {
    // Initialize function
    // mexLock();
    // load GPUmat
    gm = gmGetGPUmat();
    init = 1;
  }

  // Get GPU array pointers
  GPUtype srcArray = gm->gputype.getGPUtype(prhs[0]);

  double val       = *mxGetPr(prhs[1]);
  double operation = *mxGetPr(prhs[2]);

  int numel;
  double **atomArray = getGPUSourcePointers(prhs, 1, &numel, 0, gm);

  dim3 blocksize; blocksize.x = BLOCKDIM; blocksize.y = blocksize.z = 1;
  dim3 gridsize; gridsize.y = gridsize.z = 1;

  gridsize.x = numel / BLOCKDIM;
  if(gridsize.x * BLOCKDIM < numel) gridsize.x++;

  switch((int)operation) {
    case 1: cukern_ArraySetMin<<<gridsize, blocksize>>>(atomArray[0], val, numel); break;
    case 2: cukern_ArraySetMax<<<gridsize, blocksize>>>(atomArray[0], val, numel); break;
    case 3: cukern_ArrayFixNaN<<<gridsize, blocksize>>>(atomArray[0], val, numel); break;
  }

}

__global__ void cukern_ArraySetMin(double *array, double min, int n)
{
int x = threadIdx.x + blockDim.x * blockIdx.x;
if(x >= n) return;

if(array[x] < min) array[x] = min;
}

__global__ void cukern_ArraySetMax(double *array, double max, int n)
{
int x = threadIdx.x + blockDim.x * blockIdx.x;
if(x >= n) return;

if(array[x] > max) array[x] = max;
}

__global__ void cukern_ArrayFixNaN(double *array, double fixval, int n)
{
int x = threadIdx.x + blockDim.x * blockIdx.x;
if(x >= n) return;

if (isnan( array[x] )) { array[x] = fixval; }

}

