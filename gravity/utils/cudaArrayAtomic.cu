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

__global__ void cukern_ArraySetMin(double *array, double min,    int nu, int nv, int nw);
__global__ void cukern_ArraySetMax(double *array, double max,    int nu, int nv, int nw);
__global__ void cukern_ArrayFixNaN(double *array, double fixval, int nu, int nv, int nw);

#define BLOCKDIM 8

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

  int numDims        = gm->gputype.getNdims(srcArray);
  const int *dims    = gm->gputype.getSize(srcArray);

  dim3 gridsize, blocksize;
  int nzToPass;

  switch(numDims) {
    case 2:
      gridsize.x = dims[1]/BLOCKDIM; if(gridsize.x * BLOCKDIM < dims[1]) gridsize.x++;
      gridsize.y = 1;
      gridsize.z = 1;
      blocksize.x = BLOCKDIM;
      blocksize.y = 1;
      blocksize.z = 1;
      nzToPass = 1;
      break;
    case 3:
      gridsize.x = dims[1]/BLOCKDIM; if(gridsize.x * BLOCKDIM < dims[1]) gridsize.x++;
      gridsize.y = dims[2]/BLOCKDIM; if(gridsize.y * BLOCKDIM < dims[2]) gridsize.y++;
      gridsize.z = 1;
      blocksize.x = BLOCKDIM;
      blocksize.y = BLOCKDIM;
      blocksize.z = 1;
      nzToPass = dims[2];
      break;
  }

//printf("%i %i %i %i %i %i\n", gridsize.x, gridsize.y, gridsize.z, blocksize.x, blocksize.y, blocksize.z);

  switch((int)operation) {
    case 1: cukern_ArraySetMin<<<gridsize, blocksize>>>((double*)gm->gputype.getGPUptr(srcArray), val, dims[0], dims[1], nzToPass); break;
    case 2: cukern_ArraySetMax<<<gridsize, blocksize>>>((double*)gm->gputype.getGPUptr(srcArray), val, dims[0], dims[1], nzToPass); break;
    case 3: cukern_ArrayFixNaN<<<gridsize, blocksize>>>((double*)gm->gputype.getGPUptr(srcArray), val, dims[0], dims[1], nzToPass); break;
  }

}

__global__ void cukern_ArraySetMin(double *array, double min, int nu, int nv, int nw)
{

int myY = threadIdx.x + blockDim.x*blockIdx.x;
int myZ = threadIdx.y + blockDim.y*blockIdx.y;

int myBaseaddr = nu*(myY + nv*myZ);

if((myY >= nv) || (myZ >= nw)) return;

int xcount;
for(xcount = 0; xcount < nu; xcount++) {
  if (array[myBaseaddr] < min) { array[myBaseaddr] = min; }
//array[myBaseaddr] = 5;
  myBaseaddr++;
  }

}

__global__ void cukern_ArraySetMax(double *array, double max, int nu, int nv, int nw)
{

int myY = threadIdx.x + blockDim.x*blockIdx.x;
int myZ = threadIdx.y + blockDim.y*blockIdx.y;

int myBaseaddr = nu*(myY + nv*myZ);

if((myY >= nv) || (myZ >= nw)) return;

int xcount;
for(xcount = 0; xcount < nu; xcount++) {
  if (array[myBaseaddr] > max) { array[myBaseaddr] = max; }

  myBaseaddr++;
  }

}


__global__ void cukern_ArrayFixNaN(double *array, double fixval, int nu, int nv, int nw)
{

int myY = threadIdx.x + blockDim.x*blockIdx.x;
int myZ = threadIdx.y + blockDim.y*blockIdx.y;

int myBaseaddr = nu*(myY + nv*myZ);

if((myY >= nv) || (myZ >= nw)) return;

int xcount;
for(xcount = 0; xcount < nu; xcount++) {
  if (isnan( array[myBaseaddr] )) { array[myBaseaddr] = fixval; }

  myBaseaddr++;
  }     

}

