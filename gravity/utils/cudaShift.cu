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

#define BLOCK_DIMENSION 8

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (init == 0) {
    // Initialize function
    // mexLock();
    // load GPUmat
    gm = gmGetGPUmat();
    init = 1;
  }

  dim3 blocksize; blocksize.x = blocksize.y = BLOCK_DIMENSION; blocksize.z = 1;
  int numel; dim3 gridsize;

  if( (nlhs != 1) || (nrhs != 2)) { mexErrMsgTxt("circshift operator is shifted = cudaShift([nx ny nz], orig)"); }
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
    gridsize.x = srcsize[0] / BLOCK_DIMENSION; if(gridsize.x * BLOCK_DIMENSION < srcsize[0]) gridsize.x++;
    gridsize.y = srcsize[1] / BLOCK_DIMENSION; if(gridsize.y * BLOCK_DIMENSION < srcsize[1]) gridsize.y++;
    gridsize.z = 1;
    arrsize.x = srcsize[0]; arrsize.y = srcsize[1]; arrsize.z = srcsize[2];

    destPtr = makeGPUDestinationArrays(gm->gputype.getGPUtype(prhs[1]), plhs, 1, gm);
    cukern_circshift3D<<<gridsize, blocksize>>>(srcs[0], destPtr[0], arrsize, shift);
    free(destPtr);
    break;
    case 2:
    gridsize.x = srcsize[0] / BLOCK_DIMENSION; if(gridsize.x * BLOCK_DIMENSION < srcsize[0]) gridsize.x++;
    gridsize.y = srcsize[1] / BLOCK_DIMENSION; if(gridsize.y * BLOCK_DIMENSION < srcsize[1]) gridsize.y++;
    gridsize.z = 1; blocksize.z = 1;
    arrsize.x = srcsize[0]; arrsize.y = srcsize[1]; arrsize.z = 1;

    destPtr = makeGPUDestinationArrays(gm->gputype.getGPUtype(prhs[1]), plhs, 1, gm);
    cukern_circshift2D<<<gridsize, blocksize>>>(srcs[0], destPtr[0], arrsize, shift);
    free(destPtr);
    break;
    case 1:
    gridsize.x = srcsize[0] / BLOCK_DIMENSION; if(gridsize.x * BLOCK_DIMENSION < srcsize[0]) gridsize.x++;
    gridsize.y = 1; blocksize.y = 1;
    gridsize.z = 1; blocksize.z = 1;
    arrsize.x = srcsize[0]; arrsize.y = 1; arrsize.z = 1;

    destPtr = makeGPUDestinationArrays(gm->gputype.getGPUtype(prhs[1]), plhs, 1, gm);
    cukern_circshift1D<<<gridsize, blocksize>>>(srcs[0], destPtr[0], arrsize.x, shift.x);
    free(destPtr);
    break;
  }


}

__global__ void cukern_circshift3D(double *in, double *out, dim3 dimension, dim3 shift)
{
__shared__ double lBlock[BLOCK_DIMENSION][BLOCK_DIMENSION];

int ctrZ;

int idxX = threadIdx.x + BLOCK_DIMENSION*blockIdx.x;
int idxY = threadIdx.y + BLOCK_DIMENSION*blockIdx.y;
int idxZ = 0;

if((idxX >= dimension.x) || (idxY >= dimension.y)) return;

int idxWrite = idxX + dimension.x * idxY;

idxX = (idxX + shift.x); idxX += (idxX < 0)*dimension.x;
idxY = (idxY + shift.y); idxY += (idxY < 0)*dimension.y;

idxX = idxX % dimension.x;
idxY = idxY % dimension.y;

idxZ = shift.z; idxZ += (idxZ < 0)*dimension.z;
idxZ = idxZ % dimension.z;

int idxRead = idxX + dimension.x * (idxY + dimension.y * idxZ);

for(ctrZ = 0; ctrZ < dimension.z; ctrZ++) {
    lBlock[threadIdx.x][threadIdx.y] = in[idxRead];
    __syncthreads();

    out[idxWrite] = lBlock[threadIdx.x][threadIdx.y];

    idxWrite += dimension.x*dimension.y;
    idxRead  += dimension.x*dimension.y;
    idxRead = idxRead % (dimension.x * dimension.y * dimension.z);    

    __syncthreads();

    }

}

__global__ void cukern_circshift2D(double *in, double *out, dim3 dimension, dim3 shift)
{
__shared__ double lBlock[BLOCK_DIMENSION][BLOCK_DIMENSION];

int idxX = threadIdx.x + BLOCK_DIMENSION*blockIdx.x;
int idxY = threadIdx.y + BLOCK_DIMENSION*blockIdx.y;

if((idxX >= dimension.x) || (idxY >= dimension.y)) return;

int idxWrite = idxX + dimension.x * idxY;

idxX = (idxX + shift.x); idxX += (idxX < 0)*dimension.x;
idxY = (idxY + shift.y); idxY += (idxY < 0)*dimension.y;

idxX = idxX % dimension.x;
idxY = idxY % dimension.y;

int idxRead = idxX + dimension.x * idxY;

lBlock[threadIdx.x][threadIdx.y] = in[idxRead];
__syncthreads();
out[idxWrite] = lBlock[threadIdx.x][threadIdx.y];

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

