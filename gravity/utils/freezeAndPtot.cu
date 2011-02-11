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

__global__ void cukern_FreezeSpeed(double *rho, double *E, double *px, double *py, double *pz, double *bx, double *by, double *bz, double gam, double *freeze, double *ptot, int direct, int nu, int hu, int hv, int hw);

#define BLOCKDIM 64

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
  int direction = (int)*mxGetPr(prhs[9]);

  GPUtype srcReference = gm->gputype.getGPUtype(prhs[0]);
  int numDims          = gm->gputype.getNdims(srcReference);
  const int *dims      = gm->gputype.getSize(srcReference);

  dim3 arraySize;
  arraySize.x = dims[0];
  numDims > 1 ? arraySize.y = dims[1] : arraySize.y = 1;
  numDims > 2 ? arraySize.z = dims[2] : arraySize.z = 1;

  dim3 blocksize, gridsize;
  int hu, hv, hw, nu;

  blocksize.x = BLOCKDIM; blocksize.y = blocksize.z = 1;
  switch(direction) {
    case 1: // X direction flux: u = x, v = y, w = z;
      gridsize.x = arraySize.y;
      gridsize.y = arraySize.z;
      hu = 1; hv = arraySize.x; hw = arraySize.x * arraySize.x;
      nu = arraySize.x; break;
    case 2: // Y direction flux: u = y, v = x, w = z
      gridsize.x = arraySize.x;
      gridsize.y = arraySize.z;
      hu = arraySize.x; hv = 1; hw = arraySize.x * arraySize.y;
      nu = arraySize.y; break;
    case 3: // Z direction flux: u = z, v = x, w = y;
      gridsize.x = arraySize.x;
      gridsize.y = arraySize.y;
      hu = arraySize.x * arraySize.y; hv = 1; hw = arraySize.x;
      nu = arraySize.z; break;
    default: mexErrMsgTxt("Direction passed to directionalMaxFinder is not in { 1,2,3 }");
    }

  int numel;
  double **args = getGPUSourcePointers(prhs, 8, &numel, 0, gm);
  double **ret = makeGPUDestinationArrays(gm->gputype.getGPUtype(prhs[0]), plhs, 2, gm);

  cukern_FreezeSpeed<<<gridsize, blocksize>>>(args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], *mxGetPr(prhs[8]), ret[0], ret[1], direction, nu, hu, hv, hw);
//   b                                         (double *rho, double *E, double *px, double *py, double *pz, double *bx, double *by, double *bz, double gam, double *freeze, double *ptot, int direct, int nx, int ny, int nz)
  free(ret);
  free(args);

}

__global__ void cukern_FreezeSpeed(double *rho, double *E, double *px, double *py, double *pz, double *bx, double *by, double *bz, double gam, double *freeze, double *ptot, int direct, int nu, int hu, int hv, int hw)
{

int x = blockIdx.x * hv + blockIdx.y * hw + threadIdx.x * hu;
int addrMax = blockIdx.x * hv + blockIdx.y * hw + nu*hu;

double Cs, CsMax;
double psqhf, bsqhf;
double gg1 = gam*(gam-1.0);

__shared__ double locBloc[BLOCKDIM];

CsMax = 0.0;
locBloc[threadIdx.x] = 0.0;

if(x >= addrMax) return; // If we get a very low resolution, save time & space on wasted threads

while(x < addrMax) {
  psqhf = .5*(px[x]*px[x]+py[x]*py[x]+pz[x]*pz[x]);
  bsqhf = .5*(bx[x]*bx[x]+by[x]*by[x]+bz[x]*bz[x]);
  Cs    = sqrt(abs( (gg1*(E[x] - psqhf/rho[x]) + (4.0 - gg1)*bsqhf)/rho[x] )) + abs(px[x]/rho[x]);
  ptot[x] = (gam-1.0)*abs(E[x] - psqhf/rho[x]) + (2.0-gam)*bsqhf;

  if(Cs > CsMax) CsMax = Cs;

  x += blockDim.x * hu;
  }

locBloc[threadIdx.x] = CsMax;

// Now we need the max of the stored shared array to write back to the global array
__syncthreads();

x = 2;
while((x < BLOCKDIM) && (x < 2*nu)) {
  if(threadIdx.x % x != 0) break;

  if(locBloc[threadIdx.x + x/2] > locBloc[threadIdx.x]) locBloc[threadIdx.x] = locBloc[threadIdx.x + x/2];

  x *= 2;
  }

__syncthreads();

CsMax = locBloc[0];

x = blockIdx.x * hv + blockIdx.y * hw + threadIdx.x * hu;
while(x < addrMax) {
  freeze[x] = CsMax;
  x += blockDim.x * hu;
  }


}

