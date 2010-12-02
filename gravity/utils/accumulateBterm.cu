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
  if (nrhs!=4)
     mexErrMsgTxt("Wrong number of arguments");
  if (init == 0) {
    // Initialize function
    // mexLock();
    // load GPUmat
    gm = gmGetGPUmat();
    init = 1;
  }
  /* mex parameters are:
   0 Source term (that this applies the B operator to)
   1 Destination term (that this stores the result in)
   2 Precondition coefficient (one double)
   3 Precondition term (Accumulates successive scaled B operations
  */

  // Get GPU array pointers
  GPUtype srcArray    = gm->gputype.getGPUtype(prhs[0]);
  GPUtype dstArray    = gm->gputype.getGPUtype(prhs[1]);
  GPUtype accArray    = gm->gputype.getGPUtype(prhs[3]);

  // Get some control variables sorted out
  const int *dims    = gm->gputype.getSize(srcArray);

  dim3 gridsize;
  gridsize.x = dims[0]/EDGEDIM_BOP;
  gridsize.y = dims[1]/EDGEDIM_BOP;
  gridsize.z = 1;

  if(gridsize.x * EDGEDIM_BOP < dims[0]) gridsize.x++;
  if(gridsize.y * EDGEDIM_BOP < dims[1]) gridsize.y++;

  dim3 blocksize; blocksize.x = blocksize.y = EDGEDIM_BOP+2;
  blocksize.z = 1;

  int nx = dims[0];
  int ny = dims[1];
  int nz = dims[2];

  Laplacian_B_OperatorKernel<<<gridsize, blocksize>>>((double*)gm->gputype.getGPUptr(srcArray), (double*)gm->gputype.getGPUptr(dstArray), *mxGetPr(prhs[2]), (double*)gm->gputype.getGPUptr(accArray), nx, ny, nz, 4);

}
