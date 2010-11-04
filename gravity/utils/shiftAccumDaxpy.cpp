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
static CUfunction drvfuns[4];
static int init = 0;
static GPUmat *gm;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // At least 2 arguments expected
  // Input and result
  if (nrhs!=6)
     mexErrMsgTxt("Wrong number of arguments");
  if (init == 0) {
    // Initialize function
    // mexLock();
    // load GPUmat
    gm = gmGetGPUmat();
    init = 1;
  }
  /* mex parameters are:
   0 Source array [X or one of the y-edge-zeroed copies of it]
   1 Destination array [Accumulator, basically]
   2 Stack array 1 [We swap XY or YZ planes with this before copying to assure clean shift] - Must be all zeroes and of size max(NxNy, NyNz, NxNz)
   3 Stack array 2
   4 Shift directions
   5 Coefficient on shift */

  // Get GPU array pointers
  GPUtype srcArray    = gm->gputype.getGPUtype(prhs[0]);
  GPUtype dstArray    = gm->gputype.getGPUtype(prhs[1]);
  GPUtype stackArrayX = gm->gputype.getGPUtype(prhs[2]);
//GPUtype stackArrayY = gm->gputype.getGPUtype(prhs[3]);
  GPUtype stackArrayZ = gm->gputype.getGPUtype(prhs[3]);

  // Get some control variables sorted out
  double *shiftdirs  = mxGetPr(prhs[4]);
  const int *dims    = gm->gputype.getSize(srcArray);

  double alpha       = *mxGetPr(prhs[5]);

  int shifts[3];
  shifts[0] = (int)shiftdirs[0];
  shifts[1] = (int)shiftdirs[1];
  shifts[2] = (int)shiftdirs[2];

  double *cubSrc = (double*)gm->gputype.getGPUptr(srcArray);

  // Remove appropriate YZ plane if any
  double *cubDst = (double*)gm->gputype.getGPUptr(stackArrayX);

  if(shifts[0] == -1) cublasDswap(dims[1]*dims[2], cubSrc,             dims[0], cubDst, 1);
  if(shifts[0] ==  1) cublasDswap(dims[1]*dims[2], cubSrc + dims[0]-1, dims[0], cubDst, 1);

  // Remove appropriate XZ plane if any
  //stackSwapXZplane(cubSrc, (double*)gm->gputype.getGPUptr(stackArrayY), (int *)dims, shifts);

  // Remove appropriate XY plane if any
  cubDst = (double*)gm->gputype.getGPUptr(stackArrayZ);

  if(shifts[2] == -1) cublasDswap(dims[0]*dims[1], cubSrc,                               1, cubDst, 1);
  if(shifts[2] ==  1) cublasDswap(dims[0]*dims[1], cubSrc + dims[0]*dims[1]*(dims[2]-1), 1, cubDst, 1);

  // Decide the amount of offset to acheive desired shift
  int theta = shifts[0] + dims[0]*shifts[1] + dims[0]*dims[1]*shifts[2];
  int Ntot  = dims[0] * dims[1] * dims[2];

  cubDst = (double*)gm->gputype.getGPUptr(dstArray);
  if(theta >= 0) {
    cublasDaxpy(Ntot-theta, alpha, cubSrc,         1, cubDst + theta, 1);
  } else {
    cublasDaxpy(Ntot+theta, alpha, cubSrc - theta, 1, cubDst,         1);
  }

  // Replace the XY plane if it was removed
  cubDst = (double*)gm->gputype.getGPUptr(stackArrayZ);
  if(shifts[2] == -1) cublasDswap(dims[0]*dims[1], cubSrc,                               1, cubDst, 1);
  if(shifts[2] ==  1) cublasDswap(dims[0]*dims[1], cubSrc + dims[0]*dims[1]*(dims[2]-1), 1, cubDst, 1);

  // replace the XZ plane if it was removed
  //stackSwapXZplane(cubSrc, (double*)gm->gputype.getGPUptr(stackArrayY), (int *)dims, shifts);

  // Replace the YZ plane if it was removed
  cubDst = (double*)gm->gputype.getGPUptr(stackArrayX);
  if(shifts[0] == -1) cublasDswap(dims[1]*dims[2], cubSrc,             dims[0], cubDst, 1);
  if(shifts[0] ==  1) cublasDswap(dims[1]*dims[2], cubSrc + dims[0]-1, dims[0], cubDst, 1);
}

