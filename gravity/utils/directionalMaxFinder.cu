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

__global__ void cukern_DirectionalMax(double *d1, double *d2, double *out, int direct, int nx, int ny, int nz);
__global__ void cukern_GlobalMax(double *din, int n, double *dout);
__global__ void cukern_GlobalMax_forCFL(double *rho, double *cs, double *px, double *py, double *pz, int n, double *dout, int *dirOut);

#define BLOCKDIM 8
#define GLOBAL_BLKDIM 128

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // At least 2 arguments expected
  // Input and result
  if((nlhs == 0) || (nlhs > 2))
     mexErrMsgTxt("Either 1 return argument for simple & directional max or 2 for CFL max");

  if((nlhs == 2) && (nrhs != 5))
     mexErrMsgTxt("For CFL max require [max dir] = directionalMaxFinder(rho, soundspeed, px, py, pz)");

  if((nlhs == 1) && ((nrhs != 3) && (nrhs != 1)))
     mexErrMsgTxt("Either 1 or 3 arguments for one rturn argument");

  if (init == 0) {
    // Initialize function
    // mexLock();
    // load GPUmat
    gm = gmGetGPUmat();
    init = 1;
  }

switch(nrhs) {
  case 3: {

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

  } break;
  case 1: {
    GPUtype srcA       = gm->gputype.getGPUtype(prhs[0]);
    int numel          = gm->gputype.getNumel(srcA);

    dim3 blocksize, gridsize;
    blocksize.x = 256; blocksize.y = blocksize.z = 1;

//    gridsize.x = numel / 256; if(gridsize.x * 256 < numel) gridsize.x++;
    gridsize.x = 64;
    gridsize.y = gridsize.z =1;

    double *blkA;
    cudaMalloc(&blkA, gridsize.x * sizeof(double));

    cukern_GlobalMax<<<gridsize, blocksize>>>((double *)gm->gputype.getGPUptr(srcA), numel, blkA);

    double maxes[gridsize.x];
    cudaMemcpy(&maxes[0], blkA, sizeof(double)*gridsize.x, cudaMemcpyDeviceToHost);
    cudaFree(blkA);

    mwSize dims[2];
    dims[0] = 1;
    dims[1] = 1;
    plhs[0] = mxCreateNumericArray (2, dims, mxDOUBLE_CLASS, mxREAL);

    double *d = mxGetPr(plhs[0]);
    d[0] = maxes[0];

    for(numel = 1; numel < gridsize.x; numel++) { if(maxes[numel] > d[0]) d[0] = maxes[numel];  }
  } break;
  case 5: {
    int numel;
    double **arraysIn = getGPUSourcePointers(prhs, 5, &numel, 0, gm);

    dim3 blocksize, gridsize;
    blocksize.x = GLOBAL_BLKDIM; blocksize.y = blocksize.z = 1;

//    gridsize.x = numel / 256; if(gridsize.x * 256 < numel) gridsize.x++;
    gridsize.x = 64;
    gridsize.y = gridsize.z =1;

    double *blkA; int *blkB;
    cudaMalloc(&blkA, gridsize.x * sizeof(double));
    cudaMalloc(&blkB, gridsize.x * sizeof(int));

    cukern_GlobalMax_forCFL<<<gridsize, blocksize>>>(arraysIn[0], arraysIn[1], arraysIn[2], arraysIn[3], arraysIn[4], numel, blkA, blkB);

    double maxes[gridsize.x]; int maxIndices[gridsize.x];
    cudaMemcpy(&maxes[0], blkA, sizeof(double)*gridsize.x, cudaMemcpyDeviceToHost);
    cudaMemcpy(&maxIndices[0], blkB, sizeof(int)*gridsize.x, cudaMemcpyDeviceToHost);
    cudaFree(blkA);
    cudaFree(blkB);

    mwSize dims[2];
    dims[0] = 1;
    dims[1] = 1;
    plhs[0] = mxCreateNumericArray (2, dims, mxDOUBLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericArray (2, dims, mxDOUBLE_CLASS, mxREAL);

    double *maxout = mxGetPr(plhs[0]);
    double *dirout = mxGetPr(plhs[1]);
    maxout[0] = maxes[0];
    dirout[0] = maxIndices[0];
    
    for(numel = 1; numel < gridsize.x; numel++) { if(maxes[numel] > maxout[0]) { maxout[0] = maxes[numel]; dirout[0] = maxIndices[0]; }  }

  } break;
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

__global__ void cukern_GlobalMax(double *din, int n, double *dout)
{

int x = blockIdx.x * blockDim.x + threadIdx.x;
__shared__ double locBloc[256];

double CsMax = -1e37;
locBloc[threadIdx.x] = -1e37;
if(threadIdx.x == 0) dout[blockIdx.x] = locBloc[0]; // As a safety measure incase we return below

if(x >= n) return; // If we're fed a very small array, this will be easy

// Threads step through memory with a stride of (total # of threads), finding the max in this space
while(x < n) {
  if(din[x] > CsMax) CsMax = din[x];
  x += blockDim.x * gridDim.x;
  }
locBloc[threadIdx.x] = CsMax;

// Synchronize, then logarithmically fan in to identify each block's maximum
__syncthreads();

x = 2;
while(x < 256) {
  if(threadIdx.x % x != 0) break;

  if(locBloc[threadIdx.x + x/2] > locBloc[threadIdx.x]) locBloc[threadIdx.x] = locBloc[threadIdx.x + x/2];

  x *= 2;
  }

__syncthreads();

// Make sure the max is written and visible; each block writes one value. We test these 30 or so in CPU.
if(threadIdx.x == 0) dout[blockIdx.x] = locBloc[0];

}

// This is specifically for finding globalmax( max(abs(p_i))/rho + c_s) for the CFL constraint
__global__ void cukern_GlobalMax_forCFL(double *rho, double *cs, double *px, double *py, double *pz, int n, double *dout, int *dirOut)
{

int x = blockIdx.x * blockDim.x + threadIdx.x;
// In the end threads must share their maxima and fold them in logarithmically
__shared__ double locBloc[GLOBAL_BLKDIM];
__shared__ int locDir[GLOBAL_BLKDIM];

// Threadwise: The largest sound speed and it's directional index yet seen; Local comparision direction.
// temporary float values used to evaluate each cell
double CsMax = -1e37; int IndMax, locImax;
double tmpA, tmpB;

// Set all maxima to ~-infinity and index to invalid.
locBloc[threadIdx.x] = -1e37;
locDir[threadIdx.x] = 0;

// Have thread 0 write such to the globally shared values (overwrite garbage before we possibly get killed next line)
if(threadIdx.x == 0) { dout[blockIdx.x] = -1e37; dirOut[blockIdx.x] = 0; }

if(x >= n) return; // If we get a very low resolution, save time & space on wasted threads

// Jumping through memory,
while(x < n) {
  // Find the maximum |momentum| first; Convert it to velocity and add to soundspeed, then compare with this thread's previous max.
  tmpA = abs(px[x]);
  tmpB = abs(py[x]);

  if(tmpB > tmpA) { tmpA = tmpB; locImax = 2; } else { locImax = 1; }
  tmpB = abs(pz[x]);
  if(tmpB > tmpA) { tmpA = tmpB; locImax = 3; }  

  tmpA = tmpA / rho[x] + cs[x];

  if(tmpA > CsMax) { CsMax = tmpA; IndMax = locImax; }

  // Jump to next address to compare
  x += blockDim.x * gridDim.x;
  }

// Between them threads have surveyed entire array
// Flush threadwise maxima to shared memory
locBloc[threadIdx.x] = CsMax;
locDir[threadIdx.x] = IndMax;

// Now we need the max of the stored shared array to write back to the global array
__syncthreads();

x = 2;
while(x < GLOBAL_BLKDIM) {
  if(threadIdx.x % x != 0) return;

  if(locBloc[threadIdx.x + x/2] > locBloc[threadIdx.x]) {
	locBloc[threadIdx.x] = locBloc[threadIdx.x + x/2];
	locDir[threadIdx.x] = locDir[threadIdx.x + x/2];
	}

  x *= 2;
  }

__syncthreads();

if(threadIdx.x == 0) {
	dout[blockIdx.x] = locBloc[0];
	dirOut[blockIdx.x] = locDir[0];
	}


}


