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

__global__ void cukern_doCorrectorStep_uniform(double *rho, double *E, double *px, double *py, double *pz, double *bx, double *by, double *bz, double *P, double *Cfreeze, double *rhoW, double *enerW, double *pxW, double *pyW, double *pzW, double lambda, int nx);
__device__ void cukern_FluxLimiter_VanLeer(double deriv[2][36], double flux[2][36], int who);


#define BLOCKLEN 32

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // At least 2 arguments expected
  // Input and result
  if ((nrhs!=17) || (nlhs != 0)) mexErrMsgTxt("Wrong number of arguments: need cudaCorrectorFlux(rho, E, px, py, pz, bx, by, bz, P, rho_out, E_out, px_out, py_out, pz_out, C_freeze, lambda, purehydro?)\n");

  if (init == 0) {
    gm = gmGetGPUmat();
    init = 1;
  }

  // Get source array info and create destination arrays
  int numel;
  GPUtype srcReference = gm->gputype.getGPUtype(prhs[0]);

  double **srcs   = getGPUSourcePointers(prhs, 14, &numel, 0, gm);

  GPUtype cfreeze  = gm->gputype.getGPUtype(prhs[14]);
  double *gpu_cf = (double *)gm->gputype.getGPUptr(cfreeze);

  double lambda   = *mxGetPr(prhs[15]);
  int isPureHydro = (int)*mxGetPr(prhs[16]);

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
    cukern_doCorrectorStep_uniform<<<gridsize, blocksize>>>(srcs[0], srcs[1], srcs[2], srcs[3], srcs[4], srcs[5], srcs[6], srcs[7], srcs[8], gpu_cf, srcs[9], srcs[10], srcs[11], srcs[12], srcs[13], lambda, arraySize.x);
  }


}

/* blockidx.{xy} is our index in {yz}, and gridDim.{xy} gives the {yz} size */
/* Expect invocation with n+4 threads */
__global__ void cukern_doCorrectorStep_uniform(double *rho, double *E, double *px, double *py, double *pz, double *bx, double *by, double *bz, double *P, double *Cfreeze, double *rhoW, double *enerW, double *pxW, double *pyW, double *pzW, double lambda, int nx)
{
double Cinv, rhoinv;
double q_i[5];
double w_i;
__shared__ double fluxLR[2][36];
__shared__ double derivLR[2][36];
double *fluxdest;

/* Step 0 - obligatory annoying setup stuff (ASS) */
int I0 = nx*(blockIdx.x + gridDim.x * blockIdx.y);
int Xindex = (threadIdx.x-2);
int x; /* = Xindex % nx; */
int i;

/* Step 1 - calculate W values */
Cinv = 1.0/Cfreeze[I0];

while(Xindex < nx+2) {

    x = Xindex % nx;

    rhoinv = 1.0/rho[I0+x];
    q_i[0] = rho[I0+x];
    q_i[1] = E[I0+x];
    q_i[2] = px[I0+x];
    q_i[3] = py[I0+x];
    q_i[4] = pz[I0+x];

    /* rho, E, px, py, pz going down */
    /* Iterate over variables to flux */
    for(i = 0; i < 5; i++) {
        switch(i) {
            case 0: w_i = px[I0+x] * Cinv; break;
            case 1: w_i = (px[I0+x] * (E[I0+x] + P[I0+x]) - bx[I0+x]*(px[I0+x]*bx[I0+x]+py[I0+x]*by[I0+x]+pz[I0+x]*bz[I0+x]) ) * (rhoinv*Cinv); break;
            case 2: w_i = (px[I0+x]*px[I0+x]*rhoinv + P[I0+x] - bx[I0+x]*bx[I0+x])*Cinv; break;
            case 3: w_i = (px[I0+x]*py[I0+x]*rhoinv        - bx[I0+x]*by[I0+x])*Cinv; break;
            case 4: w_i = (px[I0+x]*pz[I0+x]*rhoinv        - bx[I0+x]*bz[I0+x])*Cinv; break;
            }

        /* Step 2 - decouple to L/R flux */
        fluxLR[0][threadIdx.x] = 0.5*(q_i[i] - w_i); /* Left  going */
        fluxLR[1][threadIdx.x] = 0.5*(q_i[i] + w_i); /* Right going */
        __syncthreads();

        /* Step 3 - Differentiate fluxes & call limiter */
            /* left flux */
        derivLR[0][threadIdx.x] = fluxLR[0][(threadIdx.x-1)%36] - fluxLR[0][threadIdx.x]; /* left derivative */
        derivLR[1][threadIdx.x] = fluxLR[0][threadIdx.x] - fluxLR[0][(threadIdx.x+1)%36]; /* right derivative */
        __syncthreads();
        cukern_FluxLimiter_VanLeer(derivLR, fluxLR, 0);
        __syncthreads();

            /* Right flux */
        derivLR[0][threadIdx.x] = fluxLR[1][threadIdx.x] - fluxLR[1][(threadIdx.x-1)%36]; /* left derivative */
        derivLR[1][threadIdx.x] = fluxLR[1][(threadIdx.x+1)%36] - fluxLR[1][threadIdx.x]; /* right derivative */
        __syncthreads();
        cukern_FluxLimiter_VanLeer(derivLR, fluxLR, 1);

        /* Step 4 - Perform flux and write to output array */
        __syncthreads();
        if( (threadIdx.x > 1) && (threadIdx.x < 34) && (Xindex < nx) ) {
            switch(i) {
                case 0: fluxdest = rhoW; break;
                case 1: fluxdest = enerW; break;
                case 2: fluxdest = pxW; break;
                case 3: fluxdest = pyW; break;
                case 4: fluxdest = pzW; break;
                }

            fluxdest[I0+x] -= lambda * ( fluxLR[0][threadIdx.x] - fluxLR[0][threadIdx.x+1] + \
                                         fluxLR[1][threadIdx.x] - fluxLR[1][threadIdx.x-1]  ) / Cinv;
            }
        __syncthreads();
        }

    Xindex += 32;
    }

}

__device__ void cukern_FluxLimiter_VanLeer(double deriv[2][36], double flux[2][36], int who)
{

double r;

r = deriv[0][threadIdx.x] * deriv[1][threadIdx.x];
if(r < 0.0) r = 0.0;

r = r / ( deriv[0][threadIdx.x] + deriv[1][threadIdx.x] );
if (isnan(r)) { r = 0.0; }

flux[who][threadIdx.x] += r;

}
