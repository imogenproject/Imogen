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

  dim3 blocksize; blocksize.x = blocksize.y = 8; blocksize.z = 1;
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
    gridsize.x = srcsize[0] / 8; if(gridsize.x * 8 < srcsize[0]) gridsize.x++;
    gridsize.y = srcsize[1] / 8; if(gridsize.x * 8 < srcsize[1]) gridsize.y++;
    gridsize.z = 1;
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

__global__ void cukern_circshift3D(double *in, double *out, dim3 dimension, dim3 shift)
{
int idxXSrc = threadIdx.x + 8*blockIdx.x;
int idxYSrc = threadIdx.y + 8*blockIdx.y;
int idxZ = 0;

if((idxXSrc >= dimension.x) || (idxYSrc >= dimension.y)) return;

int idxXDest = (idxXSrc + shift.x) % dimension.x;
int idxYDest = (idxYSrc + shift.y) % dimension.y;

if (idxXDest < 0) idxXDest = idxXDest + dimension.x;
if (idxYDest < 0) idxYDest = idxYDest + dimension.y;

__shared__ double lblock[8][8];

int idxZsrc = shift.z; if(idxZsrc < 0) idxZsrc += dimension.z;

for (idxZ = 0; idxZ < dimension.z; idxZ++){
    
    lblock[threadIdx.x][threadIdx.y] = in[idxXSrc + dimension.x*(idxYSrc + dimension.y*idxZsrc)];
    __syncthreads();

    out[idxXDest + dimension.x*(idxYDest + dimension.y*idxZ)] = lblock[threadIdx.x][threadIdx.y];
    idxZsrc = ++idxZsrc % dimension.z;
    __syncthreads();
    
    }
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

out[idxX + dimension.x*idxY] = lblock[idxX][idxY];

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

