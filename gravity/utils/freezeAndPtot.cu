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

__global__ void cukern_FreezeSpeed(double *rho, double *E, double *px, double *py, double *pz, double *bx, double *by, double *bz, double gam, double *freeze, double *ptot, int direct, int nx, int ny, int nz);

double **getSourcePointers(const mxArray *prhs[], int num, int *retNumel);
double **makeDestinationArrays(GPUtype src, mxArray *retArray[], int howmany);

#define BLOCKDIM 16

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  // At least 2 arguments expected
  // Input and result
  if ( (nrhs!=10) && (nrhs!=2))
     mexErrMsgTxt("Wrong number of arguments (require rhs, 1 lhs or 1 rhs, 1 lhs)");
  if (init == 0) {
    // Initialize function
    // mexLock();
    // load GPUmat
    gm = gmGetGPUmat();
    init = 1;
  }

  // Get GPU array pointers
  double direction = *mxGetPr(prhs[9]);

  GPUtype srcA       = gm->gputype.getGPUtype(prhs[0]);
  int numDims        = gm->gputype.getNdims(srcA);
  const int *dimsA   = gm->gputype.getSize(srcA);

  dim3 blocksize, gridsize, dims;

  dims.x = dimsA[0];
  numDims > 1 ? dims.y = dimsA[1] : dims.y = 1;
  numDims > 2 ? dims.z = dimsA[2] : dims.z = 1;

  blocksize.x = blocksize.y = BLOCKDIM; blocksize.z =1;

  switch((int)direction) {
    case 1:
	if(dims.z > 1) {
            gridsize.x = dims.y / BLOCKDIM; if (gridsize.x * BLOCKDIM < dims.y) gridsize.x++;
            gridsize.y = dims.z / BLOCKDIM; if (gridsize.y * BLOCKDIM < dims.z) gridsize.y++;
	} else {
            blocksize.x = 128; blocksize.y = 1;
            gridsize.x = dims.y / blocksize.x; if(gridsize.x * blocksize.x < dims.y) gridsize.x++;
            gridsize.y = 1;
	}
	break;
    case 2:
	if(dims.z > 1) {
            gridsize.x = dims.x / BLOCKDIM; if (gridsize.x * BLOCKDIM < dims.x) gridsize.x++;
            gridsize.y = dims.z / BLOCKDIM; if (gridsize.y * BLOCKDIM < dims.z) gridsize.y++;
	} else {
            blocksize.x = 128; blocksize.y = 1;
            gridsize.x = dims.x / blocksize.x; if(gridsize.x * blocksize.x < dims.x) gridsize.x++;
            gridsize.y = 1;
	}
	break;
    case 3:
	gridsize.x = dims.x / BLOCKDIM; if (gridsize.x * BLOCKDIM < dims.x) gridsize.x++;
        gridsize.y = dims.y / BLOCKDIM; if (gridsize.y * BLOCKDIM < dims.y) gridsize.y++;
	break;
    default: mexErrMsgTxt("Direction passed to directionalMaxFinder is not in { 1,2,3 }");
  }

  int numel;
  double **args = getSourcePointers(prhs, 8, &numel);
  double **ret = makeDestinationArrays(gm->gputype.getGPUtype(prhs[0]), plhs, 2);

  cukern_FreezeSpeed<<<gridsize, blocksize>>>(args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], *mxGetPr(prhs[8]), ret[0], ret[1], (int)direction, dims.x, dims.y, dims.z);
//   b                                         (double *rho, double *E, double *px, double *py, double *pz, double *bx, double *by, double *bz, double gam, double *freeze, double *ptot, int direct, int nx, int ny, int nz)
  free(ret);
  free(args);

}

double **getSourcePointers(const mxArray *prhs[], int num, int *retNumel)
{
  GPUtype src;
  double **gpuPointers = (double **)malloc(num * sizeof(double *));
  int iter;
  int numel = gm->gputype.getNumel(gm->gputype.getGPUtype(prhs[0]));
  for(iter = 0; iter < num; iter++) {
    src = gm->gputype.getGPUtype(prhs[iter]);
    if (gm->gputype.getNumel(src) != numel) { free(gpuPointers); mexErrMsgTxt("Fatal: Arrays contain nonequal number of elements."); }
    gpuPointers[iter] = (double *)gm->gputype.getGPUptr(src);
  }

retNumel[0] = numel;
return gpuPointers;
}

// Creates destination array that the kernels write to; Returns the GPU memory pointer, and assigns the LHS it's passed
double **makeDestinationArrays(GPUtype src, mxArray *retArray[], int howmany)
{
int d = gm->gputype.getNdims(src);
const int *ssize = gm->gputype.getSize(src);
int x;
int newsize[3];
for(x = 0; x < 3; x++) (x < d) ? newsize[x] = ssize[x] : newsize[x] = 1;

double **rvals = (double **)malloc(howmany*sizeof(double *));
int i;
for(i = 0; i < howmany; i++) {
  GPUtype ra = gm->gputype.create(gpuDOUBLE, d, newsize, NULL);
  retArray[i] = gm->gputype.createMxArray(ra);
  rvals[i] = (double *)gm->gputype.getGPUptr(ra);
  }

return rvals;

}

__global__ void cukern_FreezeSpeed(double *rho, double *E, double *px, double *py, double *pz, double *bx, double *by, double *bz, double gam, double *freeze, double *ptot, int direct, int nx, int ny, int nz)
{

int myU = threadIdx.x + blockDim.x*blockIdx.x;
int myV = threadIdx.y + blockDim.y*blockIdx.y;

double maxSoFar = -1e37;
double Cs;
double psqhf, bsqhf;
double gg1 = gam*(gam-1.0);
int addrMax, x;

switch(direct) {
  case 1: { // Seek maxima in the X direction. U=y, V=z
    if ((myU >= ny) || (myV >= nz)) return;

    x = nx*(myU + ny*myV);
    addrMax = x + nx;

    for(; x < addrMax ; x++) {
      psqhf = .5*(px[x]*px[x]+py[x]*py[x]+pz[x]*pz[x]);
      bsqhf = .5*(bx[x]*bx[x]+by[x]*by[x]+bz[x]*bz[x]);
      Cs    = sqrt(abs( (gg1*(E[x] - psqhf/rho[x]) + (4.0 - gg1)*bsqhf)/rho[x] )) + abs(px[x]/rho[x]);
      ptot[x] = (gam-1.0)*abs(E[x] - psqhf/rho[x]) + (2.0-gam)*bsqhf;

      if (Cs > maxSoFar) maxSoFar = Cs;
      }

    x = nx*(myU + ny*myV);
    __syncthreads();
    for(; x < addrMax ; x++) { freeze[x] = maxSoFar; }

  } break;
  case 2: { // Seek maxima in the Y direction. U=x, V=z
    if ((myU >= nx) || (myV >= nz)) return;

    x = myU + nx*ny*myV;
    addrMax = x + ny*nx;

    for(; x < addrMax ; x += nx) {
      psqhf = .5*(px[x]*px[x]+py[x]*py[x]+pz[x]*pz[x]);
      bsqhf = .5*(bx[x]*bx[x]+by[x]*by[x]+bz[x]*bz[x]);
      Cs    = sqrt(abs( (gg1*(E[x] - psqhf/rho[x]) + (4.0 - gg1)*bsqhf)/rho[x] )) + abs(py[x]/rho[x]);
      ptot[x] = (gam-1.0)*abs(E[x] - psqhf/rho[x]) + (2.0-gam)*bsqhf;

      if (Cs > maxSoFar) maxSoFar = Cs;
      }

    x = myU + nx*ny*myV;
    __syncthreads();
    for(; x < addrMax ; x += nx) { freeze[x] = maxSoFar; }

  } break;
  case 3: { // Seeek maxima in the Z direction; U=x, V=y
  if ((myU >= nx) || (myV >= ny)) return;

    x = myU + nx*myV;
    addrMax = x + nx*ny*nz;

    for(; x < addrMax ; x += nx*ny) {
      psqhf = .5*(px[x]*px[x]+py[x]*py[x]+pz[x]*pz[x]);
      bsqhf = .5*(bx[x]*bx[x]+by[x]*by[x]+bz[x]*bz[x]);
      Cs    = sqrt(abs( (gg1*(E[x] - psqhf/rho[x]) + (4.0 - gg1)*bsqhf)/rho[x] )) + abs(pz[x]/rho[x]); // actual sound speed plus |advection|
      ptot[x] = (gam-1.0)*abs(E[x] - psqhf/rho[x]) + (2.0-gam)*bsqhf;

      if (Cs > maxSoFar) maxSoFar = Cs;
      }

    x = myU + nx*myV;
    __syncthreads();
    for(; x < addrMax ; x += nx*ny) { freeze[x] = maxSoFar; }

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
