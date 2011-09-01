#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#ifdef UNIX
#include <stdint.h>
#include <unistd.h>
#endif
#include "mex.h"

// CUDA^M
#include "cuda.h"
#include "cuda_runtime.h"
#include "cublas.h"
#include "GPUmat.hh"

// static paramaters
static int init = 0;
static GPUmat *gm;

#include "cudaCommon.h"

#define BLOCKLEN 48
#define BLOCKLENP2 50
#define BLOCKLENP4 52

__global__ void cukern_magnetFullStep_uniform(double *velGrid, double *mag, double *bW, double *velFlow, double lambda, int nx);
__device__ void cukern_FluxLimiter_VanLeer(double deriv[2][BLOCKLENP4], double flux[2][BLOCKLENP4], int who);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // At least 2 arguments expected
  // Input and result
  if ((nrhs!=3) || (nlhs != 2)) mexErrMsgTxt("Wrong number of arguments: need [magW,velFlow] = cudaMagFullflux(mag, velgrid, lambda)\n");

  if (init == 0) {
    gm = gmGetGPUmat();
    init = 1;
  }

  // Get source array info and create destination arrays
  int numel;
  GPUtype srcReference = gm->gputype.getGPUtype(prhs[0]);

  double **srcs = getGPUSourcePointers(prhs, 2, &numel, 0, gm);
  double **dest = makeGPUDestinationArrays(srcReference,  plhs, 2, gm); 

  double lambda = *mxGetPr(prhs[2]);
  int numDims     = gm->gputype.getNdims(srcReference);
  const int *dims = gm->gputype.getSize(srcReference);

  dim3 arraySize;
  arraySize.x = dims[0];
  numDims > 1 ? arraySize.y = dims[1] : arraySize.y = 1;
  numDims > 2 ? arraySize.z = dims[2] : arraySize.z = 1;

  dim3 blocksize, gridsize;

  blocksize.x = BLOCKLEN+4;
  blocksize.y = blocksize.z = 1;

  gridsize.x = arraySize.y;
  gridsize.y = arraySize.z;

  if(arraySize.x > 1) {
    cukern_magnetFullStep_uniform<<<gridsize, blocksize>>>(srcs[0], srcs[1], dest[0], dest[1], lambda, arraySize.x);
  }


}


__global__ void cukern_magnetFullStep_uniform(double *velGrid, double *mag, double *bFull, double *velFlow, double lambda, int nx)
{
double v;
double b;
double bv;
double bflux;
double locVelFlow;
__shared__ double fluxLR[2][BLOCKLENP4];
__shared__ double derivLR[2][BLOCKLENP4];


/* Step 0 - obligatory annoying setup stuff (ASS) */
int I0 = nx*(blockIdx.x + gridDim.x * blockIdx.y);
int Xindex = (threadIdx.x-2);
int Xtrack = Xindex;
Xindex += nx*(threadIdx.x < 2);

int x;
bool doIflux = (threadIdx.x > 1) && (threadIdx.x < BLOCKLEN+2);

while(Xtrack < nx+2) {
    x = I0 + (Xindex % nx);

    b = mag[x];
    v = velGrid[x];

    //Step 1 - Calculate velocity flow
    fluxLR[0][threadIdx.x] = v;
    locVelFlow = ((fluxLR[0][threadIdx.x] + fluxLR[0][(threadIdx.x + 1)%BLOCKLENP4]) < 0); 

    //Step 2 - Calculate flux
    bv = b * v;
    fluxLR[0][threadIdx.x] = bv;
    __syncthreads();

    bflux = bv * (1 - locVelFlow) + fluxLR[0][(threadIdx.x + 1)%BLOCKLENP4]*locVelFlow;
    fluxLR[1][threadIdx.x] = bflux;
    __syncthreads();

        
    /* Step 3 - Differentiate fluxes & call limiter */
            /* left flux */
        derivLR[0][threadIdx.x] = fluxLR[0][(threadIdx.x-1)%BLOCKLENP4] - fluxLR[0][threadIdx.x]; /* left derivative */
        derivLR[1][threadIdx.x] = fluxLR[0][threadIdx.x] - fluxLR[0][(threadIdx.x+1)%BLOCKLENP4]; /* right derivative */
        __syncthreads();
        cukern_FluxLimiter_VanLeer(derivLR, fluxLR, 0);
        __syncthreads();

            /* Right flux */
        derivLR[0][threadIdx.x] = fluxLR[1][threadIdx.x] - fluxLR[1][(threadIdx.x-1)%BLOCKLENP4]; /* left derivative */
        derivLR[1][threadIdx.x] = fluxLR[1][(threadIdx.x+1)%BLOCKLENP4] - fluxLR[1][threadIdx.x]; /* right derivative */
        __syncthreads();
        cukern_FluxLimiter_VanLeer(derivLR, fluxLR, 1);

    if( doIflux && (Xindex < nx) ) {
            bFull[x] = b - lambda * ( bflux - fluxLR[1][threadIdx.x - 1] ); 
            velFlow[x] = locVelFlow;
        }

        Xindex += BLOCKLEN;
        Xtrack += BLOCKLEN;
        __syncthreads();
    }


    //bFull = mag(I).array - lambda * ( mag(I).flux(X).array - mag(I).flux(X).shift(X,-1) );



    //derivLR[0][threadIdx.x] = (bflux - flux[(threadIdx.x - 1)%BLOCKLENP4]) * (1 - locVelFlow) + (flux[threadIdx.x] - bflux) * locVelFlow;
    //__syncthreads();

    //derivLR[1][threadIdx.x] = (flux[threadIdx.x + 1)%BLOCKLENP4] - bflux) * (1 - locVelFlow) + (bflux - flux[(threadIdx.x + 2)%BLOCKLENP4]) * locVelFlow;
    //__syncthreads();
        
}

__device__ void cukern_FluxLimiter_VanLeer(double deriv[2][BLOCKLENP4], double flux[2][BLOCKLENP4], int who)
{

double r;

r = deriv[0][threadIdx.x] * deriv[1][threadIdx.x];
if(r < 0.0) r = 0.0;

r = r / ( deriv[0][threadIdx.x] + deriv[1][threadIdx.x]);
if (isnan(r)) { r = 0.0; }

flux[who][threadIdx.x] += r;

}
