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

__global__ void cukern_VanLeerLimiter(double *d1, double *d2, double *out, int nu, int nv, int nw);

#define BLOCKDIM 12

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

  // Get GPU array pointers
  GPUtype srcA = gm->gputype.getGPUtype(prhs[0]);
  GPUtype srcB = gm->gputype.getGPUtype(prhs[1]);
  GPUtype dest = gm->gputype.getGPUtype(prhs[2]);

  double operation = *mxGetPr(prhs[3]);

  int numDims        = gm->gputype.getNdims(srcA);
  const int *dims    = gm->gputype.getSize(srcA);

  dim3 gridsize;
  gridsize.x = dims[1]/BLOCKDIM;
  numDims == 3 ? gridsize.y = dims[2]/BLOCKDIM : gridsize.y = 1;
  gridsize.z = 1;

  if(gridsize.x * BLOCKDIM < dims[1]) gridsize.x++;
  if(numDims == 3) if(gridsize.y * BLOCKDIM < dims[2]) gridsize.y++;
  
  dim3 blocksize;
  blocksize.x = blocksize.y = BLOCKDIM;
  blocksize.z = 1;

//printf("%i %i %i %i %i %i\n", gridsize.x, gridsize.y, gridsize.z, blocksize.x, blocksize.y, blocksize.z);

  switch((int)operation) {
    case 1: cukern_VanLeerLimiter<<<gridsize, blocksize>>>((double*)gm->gputype.getGPUptr(srcA), (double*)gm->gputype.getGPUptr(srcB), (double*)gm->gputype.getGPUptr(dest), dims[0], dims[1], numDims == 3 ? dims[2] : 1); break;
  }

}

__global__ void cukern_VanLeerLimiter(double *d1, double *d2, double *out, int nu, int nv, int nw)
{

int myY = threadIdx.x + blockDim.x*blockIdx.x;
int myZ = threadIdx.y + blockDim.y*blockIdx.y;

int myBaseaddr = nu*(myY + nv*myZ);

if((myY >= nv) || (myZ >= nw)) return;

int xcount;
double r;

for(xcount = 0; xcount < nu; xcount++) {
  r = d1[myBaseaddr] * d2[myBaseaddr];
  if(r < 0.0) r = 0.0;

  r = r / (d1[myBaseaddr] + d2[myBaseaddr]);

  if (isnan(r)) { r = 0.0; }
//  if ((d1[myBaseaddr] + d2[myBaseaddr]) == 0.0) { r = 0.0; }
  
  out[myBaseaddr] = r;
//array[myBaseaddr] = 5;
  myBaseaddr++;
  }

}
