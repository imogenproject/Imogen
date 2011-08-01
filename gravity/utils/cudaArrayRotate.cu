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
__global__ void cukern_ArrayExchangeY(double *src, double *dst, int nx, int ny, int nz);
__global__ void cukern_ArrayExchangeZ(double *src, double *dst, int nx, int ny, int nz);

#define BDIM 16

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
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
  int newsize[3];

  if((ndims == 3) && (srcsize[2] == 1)) ndims = 2;

  GPUtype ra;

  int indExchange = (int)*mxGetPr(prhs[1]);

  switch(ndims) {
    case 3:
      if(indExchange == 2) {
        gridsize.x = srcsize[0] / BDIM; if(gridsize.x*BDIM < srcsize[0]) gridsize.x++;
        gridsize.y = srcsize[1] / BDIM; if(gridsize.y*BDIM < srcsize[1]) gridsize.y++;

        blocksize.x = blocksize.y = BDIM; blocksize.z = 1;

        newsize[0] = srcsize[1]; newsize[1] = srcsize[0]; newsize[2] = srcsize[2];
        ra = gm->gputype.create(gpuDOUBLE, 3, newsize, NULL);
        plhs[0] = gm->gputype.createMxArray(ra);
        destPtr = (double *)gm->gputype.getGPUptr(ra);

        cukern_ArrayExchangeY<<<gridsize, blocksize>>>(srcs[0], destPtr, srcsize[0], srcsize[1], srcsize[2]);
        }
      if(indExchange == 3) {
        gridsize.x = srcsize[0] / BDIM; if(gridsize.x*BDIM < srcsize[0]) gridsize.x++;
        gridsize.y = srcsize[2] / BDIM; if(gridsize.y*BDIM < srcsize[2]) gridsize.y++;

        blocksize.x = blocksize.y = BDIM; blocksize.z = 1;

        newsize[0] = srcsize[2]; newsize[1] = srcsize[1]; newsize[2] = srcsize[0];
        ra = gm->gputype.create(gpuDOUBLE, 3, newsize, NULL);
        plhs[0] = gm->gputype.createMxArray(ra);
        destPtr = (double *)gm->gputype.getGPUptr(ra);

        cukern_ArrayExchangeZ<<<gridsize, blocksize>>>(srcs[0], destPtr, srcsize[0], srcsize[1], srcsize[2]);
        }
      break;
    case 2:
      gridsize.x = srcsize[0] / BDIM; if(gridsize.x*BDIM < srcsize[0]) gridsize.x++;
      gridsize.y = srcsize[1] / BDIM; if(gridsize.y*BDIM < srcsize[1]) gridsize.y++;

      blocksize.x = blocksize.y = BDIM; blocksize.z = 1;

      newsize[0] = srcsize[1]; newsize[1] = srcsize[0]; newsize[2] = 1;
      ra = gm->gputype.create(gpuDOUBLE, 2, newsize, NULL);
      plhs[0] = gm->gputype.createMxArray(ra);
      destPtr = (double *)gm->gputype.getGPUptr(ra);
      cukern_ArrayTranspose2D<<<gridsize, blocksize>>>(srcs[0], destPtr, srcsize[0], srcsize[1]);
      break;      
    }

}

__global__ void cukern_ArrayTranspose2D(double *src, double *dst, int nx, int ny)
{
__shared__ double tmp[BDIM][BDIM];

int myx = threadIdx.x + BDIM*blockIdx.x;
int myy = threadIdx.y + BDIM*((blockIdx.y + blockIdx.x) % gridDim.y);
int myAddr = myx + nx*myy;

if((myx < nx) && (myy < ny)) tmp[threadIdx.y][threadIdx.x] = src[myAddr];

__syncthreads();

myx = threadIdx.x + BDIM*((blockIdx.y + blockIdx.x) % gridDim.y);
myy = threadIdx.y + BDIM*blockIdx.x;
myAddr = myx + ny*myy;

if((myx < ny) && (myy < nx)) dst[myAddr] = tmp[threadIdx.x][threadIdx.y];

}

__global__ void cukern_ArrayExchangeY(double *src, double *dst, int nx, int ny, int nz)
{

__shared__ double tmp[BDIM][BDIM];

int myx = threadIdx.x + BDIM*blockIdx.x;
int myy = threadIdx.y + BDIM*((blockIdx.y + blockIdx.x) % gridDim.y);
int mySrcAddr = myx + nx*myy;
bool doRead = 0;
bool doWrite = 0;

if((myx < nx) && (myy < ny)) doRead = 1; 

myx = threadIdx.x + BDIM*((blockIdx.y + blockIdx.x) % gridDim.y);
myy = threadIdx.y + BDIM*blockIdx.x;
int myDstAddr = myx + ny*myy;

if((myx < ny) && (myy < nx)) doWrite = 1;

for(myx = 0; myx < nz; myx++) {
    if(doRead) tmp[threadIdx.y][threadIdx.x] = src[mySrcAddr];
    mySrcAddr += nx*ny;
    __syncthreads();

    if(doWrite) dst[myDstAddr] = tmp[threadIdx.x][threadIdx.y];
    myDstAddr += nx*ny;
    __syncthreads();
    }

}

__global__ void cukern_ArrayExchangeZ(double*src, double *dst, int nx, int ny, int nz)
{
__shared__ double tmp[BDIM][BDIM];

int myx = threadIdx.x + BDIM*blockIdx.x;
int myz = threadIdx.y + BDIM*((blockIdx.y + blockIdx.x) % gridDim.y);
int mySrcAddr = myx + nx*ny*myz;
bool doRead = 0;
bool doWrite = 0;

if((myx < nx) && (myz < nz)) doRead = 1;

myx = threadIdx.x + BDIM*((blockIdx.y + blockIdx.x) % gridDim.y);
myz = threadIdx.y + BDIM*blockIdx.x;
int myDstAddr = myx + nz*ny*myz;

if((myx < nz) && (myz < nx)) doWrite = 1;

for(myx = 0; myx < ny; myx++) {
    if(doRead) tmp[threadIdx.y][threadIdx.x] = src[mySrcAddr];
    mySrcAddr += nx;
    __syncthreads();

    if(doWrite) dst[myDstAddr] = tmp[threadIdx.x][threadIdx.y];
    myDstAddr += nz;
    __syncthreads();
    }


}

