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

__global__ void cukern_ArrayTranspose2D(double *src, double *dst, int nx, int ny);

#define BDIM 16

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (init == 0) {
    // Initialize function
    // mexLock();
    // load GPUmat
    gm = gmGetGPUmat();
    init = 1;
  }

  dim3 blocksize; blocksize.x = blocksize.y = BDIM; blocksize.z = 1;
  int numel; dim3 gridsize;

  if((nlhs != 1) || (nrhs != 2)) { mexErrMsgTxt("cudaArrayRotate operator is rotated = cudaArrayRotate(array, dir)\n"); }
  double *rotdir = mxGetPr(prhs[1]);
  double **srcs = getGPUSourcePointers(prhs, 1, &numel, 0, gm);

  GPUtype src = gm->gputype.getGPUtype(prhs[0]);
  int ndims = gm->gputype.getNdims(src);
  const int *srcsize = gm->gputype.getSize(src);
  dim3 arrsize;
  double *destPtr;


  if((ndims == 3) && (srcsize[2] == 1)) ndims = 2;

  switch(ndims) {
    case 3: break;
    case 2:
      gridsize.x = srcsize[0] / BDIM; if(gridsize.x*BDIM < srcsize[0]) gridsize.x++;
      gridsize.y = srcsize[1] / BDIM; if(gridsize.y*BDIM < srcsize[1]) gridsize.y++;

      blocksize.x = blocksize.y = BDIM; blocksize.z = 1;

      int newsize[3]; newsize[0] = srcsize[1]; newsize[1] = srcsize[0]; newsize[2] = 1;
      GPUtype ra = gm->gputype.create(gpuDOUBLE, 2, newsize, NULL);
      plhs[0] = gm->gputype.createMxArray(ra);
      destPtr = (double *)gm->gputype.getGPUptr(ra);

//printf("%i %i\n%i %i\n%i %i\n", gridsize.x, gridsize.y, blocksize.x, blocksize.y, srcsize[0], srcsize[1]);

      cukern_ArrayTranspose2D<<<gridsize, blocksize>>>(srcs[0], destPtr, srcsize[0], srcsize[1]);
      break;      
  }


}

__global__ void cukern_ArrayTranspose2D(double *src, double *dst, int nx, int ny)
{
__shared__ double tmp[BDIM][BDIM];

int myx = threadIdx.x + BDIM*blockIdx.x;
int myy = threadIdx.y + BDIM*((blockIdx.y + blockIdx.x) % gridDim.y);

if((myx >= nx) || (myy >= ny)) return;

int myAddr = myx + nx*myy;

tmp[threadIdx.y][threadIdx.x] = src[myAddr];

__syncthreads();

myx = threadIdx.x + BDIM*((blockIdx.y + blockIdx.x) % gridDim.y);
myy = threadIdx.y + BDIM*blockIdx.x;
int myAddr2 = myx + ny*myy;

dst[myAddr2] =  tmp[threadIdx.x][threadIdx.y];

}

__global__ void cukern_ArrayExchangeY3D(double *src, double *dst, int nx, int ny, int nz)
{
__shared__ double tmp[BDIM][BDIM];
int saddr, daddr, da, zind;

saddr = threadIdx.x + BDIM*blockIdx.x;
daddr = threadIdx.y + BDIM*((blockIdx.y + blockIdx.x) % gridDim.y);

if((saddr >= nx) || (daddr > ny)) return;

da = nx*ny;

saddr = threadIdx.x + BDIM*blockIdx.x + nx*(threadIdx.y + BDIM*((blockIdx.y + blockIdx.x)%gridDim.y));
daddr = (threadIdx.x + BDIM*(blockIdx.y + blockIdx.x) % gridDim.y) + ny*(threadIdx.y + BDIM*blockIdx.x);

for(zind = 0; zind < nz; zind++) {
    tmp[threadIdx.y][threadIdx.x] = src[saddr];
    __syncthreads();
    dst[daddr] = tmp[threadIdx.x][threadIdx.y];

    saddr += da;
    daddr += da;
    }

}

__global__ void cukern_ArrayExchangeZ3D(double*src, double *dst, int nx, int ny, int nz)
{
__shared__ double tmp[BDIM][BDIM];
int saddr, daddr, zind;

saddr = threadIdx.x + BDIM*blockIdx.x;
daddr = threadIdx.y + BDIM*((blockIdx.y + blockIdx.x) % gridDim.y);

if((saddr >= nx) || (daddr > nz)) return;

saddr = threadIdx.x + BDIM*blockIdx.x + nx*ny*(threadIdx.y + BDIM*((blockIdx.y + blockIdx.x)%gridDim.y));
daddr = (threadIdx.x + BDIM*(blockIdx.y + blockIdx.x) % gridDim.y) + nx*ny*(threadIdx.y + BDIM*blockIdx.x);

for(zind = 0; zind < nz; zind++) {
    tmp[threadIdx.y][threadIdx.x] = src[saddr];
    __syncthreads();
    dst[daddr] = tmp[threadIdx.x][threadIdx.y];

    saddr += ny;
    daddr += ny;
    }


}
