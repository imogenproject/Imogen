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

__global__ void cukern_VanLeerLimiter(double *d1, double *d2, double *out, int n);

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
  int numel;
  double **sources = getGPUSourcePointers(prhs, 2, &numel, 0, gm);
  GPUtype srcRef = gm->gputype.getGPUtype(prhs[0]);

  double **dest    = makeGPUDestinationArrays(srcRef, plhs, 1, gm);

  int operation = (int)*mxGetPr(prhs[2]);

  dim3 gridsize; gridsize.x = numel / BLOCKDIM;
  if(gridsize.x *BLOCKDIM < numel) gridsize.x++;
  gridsize.y = gridsize.z = 1;

  dim3 blocksize;
  blocksize.x = BLOCKDIM;
  blocksize.y = blocksize.z = 1;

  switch((int)operation) {
    case 1: cukern_VanLeerLimiter<<<gridsize, blocksize>>>(sources[0], sources[1], dest[0], numel); break;
  }

}

__global__ void cukern_VanLeerLimiter(double *d1, double *d2, double *out, int n)
{

int x = threadIdx.x + blockIdx.x * blockDim.x;

if(x >= n) return;

double r;

r = d1[x] * d2[x];
if(r < 0.0) r = 0.0;

r = r / (d1[x] + d2[x]);

if (isnan(r)) { r = 0.0; }
  
out[x] = r;

}
