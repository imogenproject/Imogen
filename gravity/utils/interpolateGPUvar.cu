#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#ifdef UNIX
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
int fact;
int N[3];
const int *dims;
int *launchdims;

if (nrhs == 3) {
  // Get GPU array pointers if both are provided
  srcArray    = gm->gputype.getGPUtype(prhs[0]);
  dstArray    = gm->gputype.getGPUtype(prhs[2]);
printf("Writing to 3rd argument\n"); fflush(stdout);
} else if ((nlhs == 1) && (nrhs == 2)) {
  srcArray    = gm->gputype.getGPUtype(prhs[0]);
  } else mexErrMsgTxt("GPU interpolate up error: either 3 RHS or 1LHS + 2RHS arguments required\n");

// Get scaling factor
fact = (int)*mxGetPr(prhs[1]);
dims = gm->gputype.getSize(srcArray);

if (fact < 0) {
  // If scaling down, divide output size N by factor; We launch one thread per output cell
  direct = -1; 
  N[0] = -dims[0] / fact;
  N[1] = -dims[1] / fact;
  N[2] = -dims[2] / fact;
  launchdims = N;
}  else {
  // If scaling up, multiply output size N by factor; We launch one thread per input cell
  direct = 1;  
  N[0] = dims[0] * fact;
  N[1] = dims[1] * fact;
  N[2] = dims[2] * fact;
  launchdims = (int *)dims;
}
fact = abs(fact);

printf("direct: %i factor: %i\n", direct, fact); fflush(stdout);

if ((nlhs == 1) && (nrhs == 2)) {
  dstArray = gm->gputype.create(gpuDOUBLE, 3, N, NULL);
  plhs[0] = gm->gputype.createMxArray(dstArray);
}
  /* mex parameters are:
   RHS 0  - input array
   RHS 1  - scale factor
     RHS 2 or LHS 0 - increased resolution output array
  */

printf("Launch dimensions: %i %i %i\n", launchdims[0], launchdims[1], launchdims[2]); fflush(stdout);

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
    upResolutionKernel<<<gridsize, blocksize>>>((double *)gm->gputype.getGPUptr(srcArray), (double *)gm->gputype.getGPUptr(dstArray), fact, nx, ny, nz);
  else
    downResolutionKernel<<<gridsize, blocksize>>>((double *)gm->gputype.getGPUptr(srcArray), (double *)gm->gputype.getGPUptr(dstArray), fact, nx, ny, nz);
  

}


