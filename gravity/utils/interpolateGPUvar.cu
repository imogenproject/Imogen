#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#ifdef UoutputDimensionsIX
#include <stdint.h>
#include <unistd.h>
#endif

#include "mex.h"
#include "matrix.h"

// CUDA
#include "cuda.h"
#include "cuda_runtime.h"
#include "cublas.h"

// GPUmat
#include "GPUmat.hh"

// static paramaters
static int init = 0;
static GPUmat *gm;

#include "cudaKernels.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // At least 2 arguments expected
  // Input and result

if (init == 0) {
  gm = gmGetGPUmat();
  init = 1;
}

GPUtype srcArray;
GPUtype dstArray;
int direct;
int scaleFactor;
int outputDimensions[3];
const int *inputDimensions;
int *launchdims;

if (nrhs == 3) {
  // Get GPU array pointers if both are provided
  srcArray    = gm->gputype.getGPUtype(prhs[0]);
  dstArray    = gm->gputype.getGPUtype(prhs[2]);
//printf("Writing to 3rd argument\n"); fflush(stdout);
  } else if ((nlhs == 1) && (nrhs == 2)) {
  srcArray    = gm->gputype.getGPUtype(prhs[0]);
  } else mexErrMsgTxt("GPU interpolate error: either 3 RHS or 1LHS + 2RHS arguments required\n");

// Get scaling scaleFactoror
scaleFactor = (int)*mxGetPr(prhs[1]);
inputDimensions = gm->gputype.getSize(srcArray);

if (scaleFactor < 0) {
  // If scaling down, divide output size outputDimensions by scaleFactoror; We launch one thread per output cell
  scaleFactor = -scaleFactor;

  for(direct = 0; direct < 3; direct++) {
    outputDimensions[direct] = inputDimensions[direct] / scaleFactor;
    if(outputDimensions[direct]*scaleFactor < inputDimensions[direct]) outputDimensions[direct]++;
    }

  direct = -1;
  launchdims = outputDimensions;
  } else {
  // If scaling up, multiply output size outputDimensions by scaleFactoror; We launch one thread per input cell
//printf("Scaling up.\n");
  direct = 1;  
  outputDimensions[0] = inputDimensions[0] * scaleFactor;
  outputDimensions[1] = inputDimensions[1] * scaleFactor;
  outputDimensions[2] = inputDimensions[2] * scaleFactor;
  launchdims = (int *)inputDimensions;
  }

// Creating output array, it will match correct dimensions.
// If dest array is given, check for dimensional correctness.
if (nlhs == 1) {
  dstArray = gm->gputype.create(gpuDOUBLE, 3, outputDimensions, NULL);
  plhs[0] = gm->gputype.createMxArray(dstArray);
  } else {
  const int *dstdims = gm->gputype.getSize(dstArray);
  int d;
  for(d = 0; d < 3; d++) { if(dstdims[d] != outputDimensions[d]) mexErrMsgTxt("GPU interpolate error: destination array is wrong size.\n"); }
  }

dim3 gridsize;
gridsize.x = launchdims[0]/8;
gridsize.y = launchdims[1]/8;
gridsize.z = 1;

if(gridsize.x * 8 < launchdims[0]) gridsize.x++;
if(gridsize.y * 8 < launchdims[1]) gridsize.y++;

dim3 blocksize; blocksize.x = blocksize.y = 8;
blocksize.z = 1;

int nx = launchdims[0];
int ny = launchdims[1];
int nz = launchdims[2];
  if(nz == 0) nz = 1;

  if(direct > 0)
      upsampleKernel<<<gridsize, blocksize>>>((double *)gm->gputype.getGPUptr(srcArray), (double *)gm->gputype.getGPUptr(dstArray), scaleFactor, nx, ny, nz);
  else
    downsampleKernel<<<gridsize, blocksize>>>((double *)gm->gputype.getGPUptr(srcArray), (double *)gm->gputype.getGPUptr(dstArray), scaleFactor, nx, ny, nz);
  

}


