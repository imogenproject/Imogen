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

__global__ void cukern_circshift3D(double *in, double *out, dim3 dimension, dim3 shift);
__global__ void cukern_circshift2D(double *in, double *out, dim3 dimension, dim3 shift);
__global__ void cukern_circshift1D(double *in, double *out, int dimension, int shift);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (init == 0) {
    // Initialize function
    // mexLock();
    // load GPUmat
    gm = gmGetGPUmat();
    init = 1;
  }

  dim3 blocksize; blocksize.x = blocksize.y = blocksize.z = 8;
  int numel; dim3 gridsize;

  if( (nlhs != 1) || (nrhs != 2)) { mexErrMsgTxt("circshift operator is shifted = cudaShift(orig, [nx ny nz])"); }
  double *shiftamt = mxGetPr(prhs[0]);
  double **srcs = getGPUSourcePointers(prhs, 1, &numel, 1, gm);

  dim3 shift;
  shift.x = (int)shiftamt[0];
  shift.y = (int)shiftamt[1];
  shift.z = (int)shiftamt[2];

  GPUtype src = gm->gputype.getGPUtype(prhs[1]);
  int ndims = gm->gputype.getNdims(src);
  const int *srcsize = gm->gputype.getSize(src);
  dim3 arrsize;
  double **destPtr;

  switch(ndims) {
    case 3:
    gridsize.x = srcsize[0] / 8; if(gridsize.x * 8 < srcsize[0]) gridsize.x++;
    gridsize.y = srcsize[1] / 8; if(gridsize.x * 8 < srcsize[1]) gridsize.y++;
    gridsize.z = srcsize[2] / 8; if(gridsize.x * 8 < srcsize[2]) gridsize.z++;
    arrsize.x = srcsize[0]; arrsize.y = srcsize[1]; arrsize.z = srcsize[2];

    destPtr = makeGPUDestinationArrays(gm->gputype.getGPUtype(prhs[1]), plhs, 1, gm);
    cukern_circshift3D<<<gridsize, blocksize>>>(srcs[0], destPtr[0], arrsize, shift);
    free(destPtr);
    break;
    case 2:
    gridsize.x = srcsize[0] / 8; if(gridsize.x * 8 < srcsize[0]) gridsize.x++;
    gridsize.y = srcsize[1] / 8; if(gridsize.x * 8 < srcsize[1]) gridsize.y++;
    gridsize.z = 1; blocksize.z = 1;
    arrsize.x = srcsize[0]; arrsize.y = srcsize[1]; arrsize.z = 1;

//intf("%i %i %i %i %i %i\n", gridsize.x, gridsize.y, arrsize.x, arrsize.y, shift.x, shift.y);

    destPtr = makeGPUDestinationArrays(gm->gputype.getGPUtype(prhs[1]), plhs, 1, gm);
    cukern_circshift2D<<<gridsize, blocksize>>>(srcs[0], destPtr[0], arrsize, shift);
    free(destPtr);
    break;
    case 1:
    gridsize.x = srcsize[0] / 8; if(gridsize.x * 8 < srcsize[0]) gridsize.x++;
    gridsize.y = 1; blocksize.y = 1;
    gridsize.z = 1; blocksize.z = 1;
    arrsize.x = srcsize[0]; arrsize.y = 1; arrsize.z = 1;

    destPtr = makeGPUDestinationArrays(gm->gputype.getGPUtype(prhs[1]), plhs, 1, gm);
    cukern_circshift1D<<<gridsize, blocksize>>>(srcs[0], destPtr[0], arrsize.x, shift.x);
    free(destPtr);
    break;
  }


}

//#define KERNEL_PREAMBLE int x = THREADLOOPS*(threadIdx.x + blockDim.x*blockIdx.x); if (x >= n) {return;} int imax; ((x+THREADLOOPS) > n) ? imax = n : imax = x + THREADLOOPS; for(; x < imax; x++)
//#define KERNEL_PREAMBLE int x = threadIdx.x + blockDim.x*blockIdx.x; if (x >= n) { return; }

// THIS KERNEL CALCULATES SOUNDSPEED 
__global__ void cukern_circshift3D(double *in, double *out, dim3 dimension, dim3 shift)
{
int idxX = threadIdx.x + 8*blockIdx.x;
int idxY = threadIdx.y + 8*blockIdx.y;
int idxZ = threadIdx.z + 8*blockIdx.z;

if((idxX >= dimension.x) || (idxY >= dimension.y) || (idxZ >= dimension.z)) return;

idxX = (idxX + shift.x) % dimension.x;
idxY = (idxY + shift.y) % dimension.y;
idxZ = (idxZ + shift.z) % dimension.z;

__shared__ double lblock[8][8][8];

lblock[idxX][idxY][idxZ] = in[idxX + dimension.x*(idxY + dimension.y*idxZ)];

__syncthreads();

idxX = threadIdx.x + 8*blockIdx.x;
idxY = threadIdx.y + 8*blockIdx.y;
idxZ = threadIdx.z + 8*blockIdx.z;

out[idxX + dimension.x*(idxY + dimension.y*idxZ)] = lblock[idxX][idxY][idxZ];

}

__global__ void cukern_circshift2D(double *in, double *out, dim3 dimension, dim3 shift)
{
int idxX = threadIdx.x + 8*blockIdx.x;
int idxY = threadIdx.y + 8*blockIdx.y;

if((idxX >= dimension.x) || (idxY >= dimension.y)) return;

idxX = (idxX + shift.x) % dimension.x;
idxY = (idxY + shift.y) % dimension.y;

__shared__ double lblock[8][8];

lblock[threadIdx.x][threadIdx.y] = in[idxX + dimension.x*idxY];

__syncthreads();

idxX = threadIdx.x + 8*blockIdx.x;
idxY = threadIdx.y + 8*blockIdx.y;

out[idxX + dimension.x*idxY] = lblock[threadIdx.x][threadIdx.y];

}

__global__ void cukern_circshift1D(double *in, double *out, int dimension, int shift)
{
int idxX = threadIdx.x + 8*blockIdx.x;

if((idxX >= dimension)) return;

idxX = (idxX + shift) % dimension;

__shared__ double lblock[8];

lblock[idxX] = in[idxX];

__syncthreads();

idxX = threadIdx.x + 8*blockIdx.x;

out[idxX] = lblock[idxX];

}
