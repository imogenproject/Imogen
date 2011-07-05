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

__global__ void cukern_blitzero(double *fuck, int nmax);
__global__ void cukern_ArrayTranspose2D(double *src, double *dst, int nx, int ny);
__global__ void cukern_ArrayExchangeY3D(double *src, double *dst, int nx, int ny, int nz);
__global__ void cukern_ArrayExchangeZ3D(double *src, double *dst, int nx, int ny, int nz);

#define BDIM 3

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
  int newsize[3];

//printf("srcsize: %i %i %i. nd=%i\n", srcsize[0], srcsize[1], srcsize[2], ndims);

  if((ndims == 3) && (srcsize[2] == 1)) ndims = 2;

  GPUtype ra;

  switch(ndims) {
    case 3: {
      gridsize.x = srcsize[0] / BDIM; if(gridsize.x*BDIM < srcsize[0]) gridsize.x++;
      int indExchange = (int)*mxGetPr(prhs[1]);
      if (indExchange == 2) {
        gridsize.y = srcsize[1] / BDIM; if(gridsize.y*BDIM < srcsize[1]) gridsize.y++;
        newsize[0] = srcsize[1]; newsize[1] = srcsize[0]; newsize[2] = srcsize[2];
        }
      if (indExchange == 3) {
        gridsize.y = srcsize[2] / BDIM; if(gridsize.y*BDIM < srcsize[2]) gridsize.y++;
        newsize[0] = srcsize[2]; newsize[1] = srcsize[1]; newsize[2] = srcsize[0];
        }
      blocksize.x = blocksize.y = BDIM; blocksize.z = 1;

      ra = gm->gputype.create(gpuDOUBLE, 3, newsize, NULL);
      plhs[0] = gm->gputype.createMxArray(ra);
      destPtr = (double *)gm->gputype.getGPUptr(ra);

//cukern_blitzero<<<gridsize.x*gridsize.y, BDIM>>>(destPtr, srcsize[0]*srcsize[1]*srcsize[2]);

      if (indExchange == 2) {
        cukern_ArrayExchangeY3D<<<gridsize, blocksize>>>(srcs[0], destPtr, srcsize[0], srcsize[1], srcsize[2]);
        }
      if (indExchange == 3) {
        cukern_ArrayExchangeZ3D<<<gridsize, blocksize>>>(srcs[0], destPtr, srcsize[0], srcsize[1], srcsize[2]); 
        }
      } break;
    case 2: {
      gridsize.x = srcsize[0] / BDIM; if(gridsize.x*BDIM < srcsize[0]) gridsize.x++;
      gridsize.y = srcsize[1] / BDIM; if(gridsize.y*BDIM < srcsize[1]) gridsize.y++;

      blocksize.x = blocksize.y = BDIM; blocksize.z = 1;

      newsize[0] = srcsize[1]; newsize[1] = srcsize[0]; newsize[2] = 1;
      ra = gm->gputype.create(gpuDOUBLE, 2, newsize, NULL);
      plhs[0] = gm->gputype.createMxArray(ra);
      destPtr = (double *)gm->gputype.getGPUptr(ra);

//printf("%i %i\n%i %i\n%i %i\n", gridsize.x, gridsize.y, blocksize.x, blocksize.y, srcsize[0], srcsize[1]);

      cukern_ArrayTranspose2D<<<gridsize, blocksize>>>(srcs[0], destPtr, srcsize[0], srcsize[1]);
      } break;      
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

int myx = threadIdx.x + BDIM*blockIdx.x;
int myy = threadIdx.y + BDIM*((blockIdx.y + blockIdx.x) % gridDim.y);

if((myx >= nx) || (myy >= ny)) return;

int myAddr = myx + nx*myy;

myx = threadIdx.x + BDIM*((blockIdx.y + blockIdx.x) % gridDim.y);
myy = threadIdx.y + BDIM*blockIdx.x;
int myAddr2 = myx + ny*myy;

int zct;
for(zct = 0; zct < nz; zct++) {
    tmp[threadIdx.y][threadIdx.x] = src[myAddr];
    __syncthreads();
    dst[myAddr2] =  myx;//tmp[threadIdx.x][threadIdx.y];
    __syncthreads();

    myAddr += nx*ny;
    myAddr2 += nx*ny;
    }

}

__global__ void cukern_blitzero(double *fuck, int nmax)
{
int n = BDIM*blockIdx.x + threadIdx.x;
fuck[n] = 999.0;
}

__global__ void cukern_ArrayExchangeZ3D(double*src, double *dst, int nx, int ny, int nz)
{
__shared__ double tmp[BDIM][BDIM];
int saddr, daddr, yind;

saddr = threadIdx.x + BDIM*blockIdx.x;
daddr = threadIdx.y + BDIM*((blockIdx.y + blockIdx.x) % gridDim.y);

if((saddr >= nx) || (daddr > nz)) return;

//saddr = threadIdx.x + BDIM*blockIdx.x + nx*ny*(threadIdx.y + BDIM*((blockIdx.y + blockIdx.x)%gridDim.y));
saddr = threadIdx.x + BDIM*blockIdx.x + nx*ny*(threadIdx.y + BDIM*blockIdx.y);
//daddr = threadIdx.x + (BDIM*(blockIdx.y + blockIdx.x) % gridDim.y) + nx*ny*(threadIdx.y + BDIM*blockIdx.x);
daddr = threadIdx.x + BDIM*blockIdx.y + nz*ny*(threadIdx.y + BDIM*blockIdx.x);

for(yind = 0; yind < ny; yind++) {
    tmp[threadIdx.y][threadIdx.x] = src[saddr];
    __syncthreads();
    dst[daddr] = tmp[threadIdx.x][threadIdx.y];

    saddr += nx;
    daddr += nz;
    __syncthreads(); 
    }

}
