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

__global__ void cukern_DirectionalMax(double *d1, double *d2, double *out, int direct, int nx, int ny, int nz);
__global__ void cukern_GlobalMax(double *din, int n, double *dout);

#define BLOCKDIM 8

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // At least 2 arguments expected
  // Input and result
  if ( ( (nrhs!=3) && (nrhs!=1)) || (nlhs != 1) )
     mexErrMsgTxt("Wrong number of arguments (require 3 rhs, 1 lhs or 1 rhs, 1 lhs)");
  if (init == 0) {
    // Initialize function
    // mexLock();
    // load GPUmat
    gm = gmGetGPUmat();
    init = 1;
  }

if (nrhs == 3) {

  // Get GPU array pointers
  GPUtype srcA = gm->gputype.getGPUtype(prhs[0]);
  GPUtype srcB = gm->gputype.getGPUtype(prhs[1]);

  double operation = *mxGetPr(prhs[2]);

  int numDims        = gm->gputype.getNdims(srcA);
  const int *dimsA   = gm->gputype.getSize(srcA);
  const int *dimsB   = gm->gputype.getSize(srcB);

  GPUtype dest = gm->gputype.create(gpuDOUBLE, numDims, dimsA, NULL);
  plhs[0] = gm->gputype.createMxArray(dest);

  dim3 blocksize, gridsize, dims;

  dims.x = dimsA[0];
  numDims > 1 ? dims.y = dimsA[1] : dims.y = 1;
  numDims > 2 ? dims.z = dimsA[2] : dims.z = 1;

  blocksize.x = blocksize.y = BLOCKDIM; blocksize.z =1;

  switch((int)*mxGetPr(prhs[2])) {
    case 1:
	gridsize.x = dims.y / BLOCKDIM; if (gridsize.x * BLOCKDIM < dims.y) gridsize.x++;
	gridsize.y = dims.z / BLOCKDIM; if (gridsize.y * BLOCKDIM < dims.z) gridsize.y++;
	break;
    case 2:
        gridsize.x = dims.x / BLOCKDIM; if (gridsize.x * BLOCKDIM < dims.x) gridsize.x++;
        gridsize.y = dims.z / BLOCKDIM; if (gridsize.y * BLOCKDIM < dims.z) gridsize.y++;
	break;
    case 3:
	gridsize.x = dims.x / BLOCKDIM; if (gridsize.x * BLOCKDIM < dims.x) gridsize.x++;
        gridsize.y = dims.y / BLOCKDIM; if (gridsize.y * BLOCKDIM < dims.y) gridsize.y++;
	break;
    default: mexErrMsgTxt("Direction passed to directionalMaxFinder is not in { 1,2,3 }");
  }
  
  //printf("%i %i %i %i %i %i\n", gridsize.x, gridsize.y, gridsize.z, blocksize.x, blocksize.y, blocksize.z);
  cukern_DirectionalMax<<<gridsize, blocksize>>>((double*)gm->gputype.getGPUptr(srcA), (double*)gm->gputype.getGPUptr(srcB), (double*)gm->gputype.getGPUptr(dest), (int)*mxGetPr(prhs[2]), dims.x, dims.y, dims.z);
  } else {
    GPUtype srcA       = gm->gputype.getGPUtype(prhs[0]);
    int numel          = gm->gputype.getNumel(srcA);

    dim3 blocksize, gridsize;
    blocksize.x = 256; blocksize.y = blocksize.z = 1;
//printf("%i elements to start\n", numel);

    gridsize.x = numel / 256; if(gridsize.x * 256 < numel) gridsize.x++;
    gridsize.y = gridsize.z =1;
//printf("%i gridsize to start\n", gridsize.x);
    double *blkA; double *blkB = NULL;
    cudaMalloc(&blkA, gridsize.x * sizeof(double));

    cukern_GlobalMax<<<gridsize, blocksize>>>((double *)gm->gputype.getGPUptr(srcA), numel, blkA);

    while(gridsize.x > 256) { // Perform successive factor-of-256 fanins until we're left with <= 256 elements to search by CPU
      numel = gridsize.x;
      gridsize.x = gridsize.x / 256;
      if(gridsize.x * 256 < numel) gridsize.x++;
//double maxI[numel];
//cudaMemcpy(&maxI[0], blkA, sizeof(double)*numel, cudaMemcpyDeviceToHost);
//printf("Intermediate maxima: ");
//int q; for(q = 0; q < numel; q++) { printf("%.4lF ", maxI[q]); } printf("\n");


      if(blkB != NULL) cudaFree(blkB);
      blkB = blkA;
      cudaMalloc(&blkA, gridsize.x * sizeof(double));
      cukern_GlobalMax<<<gridsize, blocksize>>>(blkB, numel, blkA);

    }

  double maxes[gridsize.x];
  cudaMemcpy(&maxes[0], blkA, sizeof(double)*gridsize.x, cudaMemcpyDeviceToHost);

  mwSize dims[2];
  dims[0] = 1;
  dims[1] = 1;
  plhs[0] = mxCreateNumericArray (2, dims, mxDOUBLE_CLASS, mxREAL);

  double *d = mxGetPr(plhs[0]);
  d[0] = maxes[0];
//printf("final searching: ");
  for(numel = 1; numel < gridsize.x; numel++) { if(maxes[numel] > d[0]) d[0] = maxes[numel];  }
//printf("\n");
  }

}

__global__ void cukern_DirectionalMax(double *d1, double *d2, double *out, int direct, int nx, int ny, int nz)
{

int myU = threadIdx.x + blockDim.x*blockIdx.x;
int myV = threadIdx.y + blockDim.y*blockIdx.y;

double maxSoFar = -1e37;
int addrMax, myBaseaddr;

switch(direct) {
  case 1: { // Seek maxima in the X direction. U=y, V=z
    if ((myU >= ny) || (myV >= nz)) return;

    myBaseaddr = nx*(myU + ny*myV);
    addrMax = myBaseaddr + nx;

    for(; myBaseaddr < addrMax ; myBaseaddr++) {
      if ( abs(d1[myBaseaddr]) + d2[myBaseaddr] > maxSoFar) maxSoFar = abs(d1[myBaseaddr]) + d2[myBaseaddr];
      }

    myBaseaddr = nx*(myU + ny*myV);
    for(; myBaseaddr < addrMax ; myBaseaddr++) { out[myBaseaddr] = maxSoFar; }

  } break;
  case 2: { // Seek maxima in the Y direction. U=x, V=z
    if ((myU >= nx) || (myV >= nz)) return;

    myBaseaddr = myU + nx*ny*myV;
    addrMax = myBaseaddr + ny*nx;

    for(; myBaseaddr < addrMax ; myBaseaddr += nx) {
      if ( abs(d1[myBaseaddr]) + d2[myBaseaddr] > maxSoFar) maxSoFar = abs(d1[myBaseaddr]) + d2[myBaseaddr];
      }

    myBaseaddr = myU + nx*ny*myV;
    for(; myBaseaddr < addrMax ; myBaseaddr += nx) { out[myBaseaddr] = maxSoFar; }

  } break;
  case 3: { // Seeek maxima in the Z direction; U=x, V=y
  if ((myU >= nx) || (myV >= ny)) return;

    myBaseaddr = myU + nx*myV;
    addrMax = myBaseaddr + nx*ny*nz;

    for(; myBaseaddr < addrMax ; myBaseaddr += nx*ny) {
      if ( abs(d1[myBaseaddr]) + d2[myBaseaddr] > maxSoFar) maxSoFar = abs(d1[myBaseaddr]) + d2[myBaseaddr];
      }

    myBaseaddr = myU + nx*myV;
    for(; myBaseaddr < addrMax ; myBaseaddr += nx) { out[myBaseaddr] = maxSoFar; }

  } break;
}

}

// This must be invoked with 2^n threads, n >= 8 for efficiency
__global__ void cukern_GlobalMax(double *din, int n, double *dout)
{
int addr = threadIdx.x + blockDim.x*blockIdx.x;
__shared__ double loc[256];

if (addr >= n) { loc[threadIdx.x] = -1e37; return; }

loc[threadIdx.x] = din[addr];

__syncthreads();

// 256 threads here <-
if(( threadIdx.x % 2) > 0 ) return;
if (loc[threadIdx.x+1] > loc[threadIdx.x]) loc[threadIdx.x] = loc[threadIdx.x+1];
__syncthreads();
// 128 threads here <=
if(threadIdx.x % 4) return;
if (loc[threadIdx.x+2] > loc[threadIdx.x]) loc[threadIdx.x] = loc[threadIdx.x+2];
__syncthreads();
// 64 threads here <=
if(threadIdx.x % 8) return;
if (loc[threadIdx.x+4] > loc[threadIdx.x]) loc[threadIdx.x] = loc[threadIdx.x+4];
__syncthreads();
// 32 threads here <= (last full warp at n=8)
if(threadIdx.x % 16) return;
if (loc[threadIdx.x+8] > loc[threadIdx.x]) loc[threadIdx.x] = loc[threadIdx.x+8];
__syncthreads();
// 16 threads here <=
if(threadIdx.x % 32) return;
if (loc[threadIdx.x+16] > loc[threadIdx.x]) loc[threadIdx.x] = loc[threadIdx.x+16];
__syncthreads();
// 8 threads here <=

// Continuing with a nigh-empty warp is not profitable
// One serial loop over the last 8 values to identify the max.

if(threadIdx.x > 0) return;
for(addr = 32; addr < blockDim.x; addr += 32) { if (loc[addr] > loc[0]) loc[0] = loc[addr]; }

dout[blockIdx.x] = loc[0];
}

/*
__global__ void cukern_GlobalMax(double *din, int n, double *dout)
{

int x = blockIdx.x * blockDim.x + threadIdx.x;
__shared__ double locBloc[256];

double CsMax = -1e37;
locBloc[threadIdx.x] = -1e37;

if(x >= n) return; // If we get a very low resolution, save time & space on wasted threads

// Every block 
while(x < n) {
  if(din[x] > CsMax) CsMax = din[x];
  x += blockDim.x * gridDim.x;
  }

locBloc[threadIdx.x] = CsMax;

// Now we need the max of the stored shared array to write back to the global array
__syncthreads();

x = 2;
while(x < BLOCKDIM) {
  if(threadIdx.x % x != 0) break;

  if(locBloc[threadIdx.x + x/2] > locBloc[threadIdx.x]) locBloc[threadIdx.x] = locBloc[threadIdx.x + x/2];

  x *= 2;
  }

__syncthreads();

if(threadIdx.x == 0) dout[blockIdx.x] = locBloc[0];


}
*/
