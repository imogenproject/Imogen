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

__global__ void cukern_FreezeSpeed(double *rho, double *E, double *px, double *py, double *pz, double *bx, double *by, double *bz, double gam, double *freeze, double *ptot, int nx);
__global__ void cukern_FreezeSpeed_hydro(double *rho, double *E, double *px, double *py, double *pz, double gam, double *freeze, double *ptot, int nx);

#define BLOCKDIM 64
#define MAXPOW   5

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  // At least 2 arguments expected
  // Input and result
  if ( (nrhs!=10) && (nrhs!=2))
     mexErrMsgTxt("Wrong number of arguments. Call using [ptot freeze] = FreezeAndPtot(mass, ener, momx, momy, momz, bz, by, bz, gamma, 1)");
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

  blocksize.x = BLOCKDIM; blocksize.y = blocksize.z = 1;
  gridsize.x = arraySize.y;
  gridsize.y = arraySize.z;

  int numel;
  double **args = getGPUSourcePointers(prhs, 8, &numel, 0, gm);
  double **ret = makeGPUDestinationArrays(gm->gputype.getGPUtype(prhs[0]), plhs, 1, gm); // ptotal array

  double *freezea;
  int newsize[2];
  newsize[0] = arraySize.y;
  newsize[1] = arraySize.z;
  GPUtype ra = gm->gputype.create(gpuDOUBLE, 2, newsize, NULL);
  plhs[1] = gm->gputype.createMxArray(ra);
  freezea = (double *)gm->gputype.getGPUptr(ra);

  int ispurehydro = (int)*mxGetPr(prhs[9]);

  if(ispurehydro) {
    cukern_FreezeSpeed_hydro<<<gridsize, blocksize>>>(args[0], args[1], args[2], args[3], args[4],  *mxGetPr(prhs[8]), freezea, ret[0], arraySize.x);
//                                                   (*rho,    *E,      *px,     *py,     *pz,      gam,              *freeze,  *ptot,  nx)
  } else {
    cukern_FreezeSpeed<<<gridsize, blocksize>>>(args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], *mxGetPr(prhs[8]), freezea, ret[0], arraySize.x);
//                                             (*rho,    *E,      *px,     *py,     *pz,     *bx,     *by,     *bz,     gam,              *freeze,  *ptot,  nx)
  }
  free(ret);
  free(args);

}

__global__ void cukern_FreezeSpeed(double *rho, double *E, double *px, double *py, double *pz, double *bx, double *by, double *bz, double gam, double *freeze, double *ptot, int nx)
{
/* gridDim = [ny nz], nx = nx */
int x = threadIdx.x + nx*(blockIdx.x + gridDim.x*blockIdx.y);
int addrMax = nx + nx*(blockIdx.x + gridDim.x*blockIdx.y);

double Cs, CsMax;
double psqhf, bsqhf;
double gg1 = gam*(gam-1.0);

__shared__ double locBloc[BLOCKDIM];

CsMax = 0.0;
locBloc[threadIdx.x] = 0.0;

if(x >= addrMax) return; // If we get a very low resolution

while(x < addrMax) {
  psqhf = .5*(px[x]*px[x]+py[x]*py[x]+pz[x]*pz[x]);
  bsqhf = .5*(bx[x]*bx[x]+by[x]*by[x]+bz[x]*bz[x]);
  Cs    = sqrt(abs( (gg1*(E[x] - psqhf/rho[x]) + (4.0 - gg1)*bsqhf)/rho[x] )) + abs(px[x]/rho[x]);
  ptot[x] = (gam-1.0)*abs(E[x] - psqhf/rho[x]) + (2.0-gam)*bsqhf;

  if(Cs > CsMax) CsMax = Cs;

  x += BLOCKDIM;
  }

locBloc[threadIdx.x] = CsMax;

// Now we need the max of the shared array to write back to the global C_f array
__syncthreads();

x = 2;
addrMax = 1;

while(x <= BLOCKDIM) {
  if(threadIdx.x % x != 0) return;

  if(locBloc[threadIdx.x + addrMax] > locBloc[threadIdx.x]) locBloc[threadIdx.x] = locBloc[threadIdx.x + addrMax];

  addrMax = x;
  x *= 2;
  }
//if (threadIdx.x > 0) return;
//CsMax = 0;
//for(x = 0; x < 64; x++) { if(locBloc[x] > CsMax) CsMax = locBloc[x]; }

freeze[blockIdx.x + gridDim.x*blockIdx.y] = locBloc[0];

}

__global__ void cukern_FreezeSpeed_hydro(double *rho, double *E, double *px, double *py, double *pz, double gam, double *freeze, double *ptot, int nx)
{
int x = threadIdx.x + nx*(blockIdx.x + gridDim.x*blockIdx.y);
int addrMax = nx + nx*(blockIdx.x + gridDim.x*blockIdx.y);

double Cs, CsMax;
double psqhf;
double gg1 = gam*(gam-1.0);

__shared__ double locBloc[BLOCKDIM];

CsMax = 0.0;
locBloc[threadIdx.x] = 0.0;

if(x >= addrMax) return; // If we get a very low resolution

while(x < addrMax) {
  psqhf = .5*(px[x]*px[x]+py[x]*py[x]+pz[x]*pz[x]);
  Cs    = sqrt(abs( (gg1*(E[x] - psqhf/rho[x]) )/rho[x] )) + abs(px[x]/rho[x]);
  ptot[x] = (gam-1.0)*abs(E[x] - psqhf/rho[x]);

  if(Cs > CsMax) CsMax = Cs;

  x += BLOCKDIM;
  }

locBloc[threadIdx.x] = CsMax;

// Now we need the max of the shared array to write back to the global C_f array
__syncthreads();

x = 2;
addrMax = 1;

while(x <= BLOCKDIM) {
  if(threadIdx.x % x != 0) return;

  if(locBloc[threadIdx.x + addrMax] > locBloc[threadIdx.x]) locBloc[threadIdx.x] = locBloc[threadIdx.x + addrMax];

  addrMax = x;
  x *= 2;
  }
//if (threadIdx.x > 0) return;
//CsMax = 0;
//for(x = 0; x < 64; x++) { if(locBloc[x] > CsMax) CsMax = locBloc[x]; }

freeze[blockIdx.x + gridDim.x*blockIdx.y] = locBloc[0];


}
