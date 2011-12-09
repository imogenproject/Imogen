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

__global__ void cukern_TVDStep_mhd_uniform(double *rho, double *E, double *px, double *py, double *pz, double *bx, double *by, double *bz, double *P, double *Cfreeze, double *rhoW, double *enerW, double *pxW, double *pyW, double *pzW, double lambda, int nx);
__global__ void cukern_TVDStep_hydro_uniform(double *rho, double *E, double *px, double *py, double *pz, double *P, double *Cfreeze, double *rhoW, double *enerW, double *pxW, double *pyW, double *pzW, double lambda, int nx);
__device__ void cukern_FluxLimiter_VanLeer(double deriv[2][BLOCKLENP4], double flux[2][BLOCKLENP4], int who);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // At least 2 arguments expected
  // Input and result
  if ((nrhs!=17) || (nlhs != 0)) mexErrMsgTxt("Wrong number of arguments: need cudaTVDStep(rho, E, px, py, pz, bx, by, bz, P, rho_out, E_out, px_out, py_out, pz_out, C_freeze, lambda, purehydro?)\n");

  if (init == 0) {
    gm = gmGetGPUmat();
    init = 1;
  }

  // Get source array info and create destination arrays
  int numel;
  GPUtype srcReference = gm->gputype.getGPUtype(prhs[0]);

  // Get the source gpu arrays as enumerated in the error message
  double **srcs   = getGPUSourcePointers(prhs, 14, &numel, 0, gm);

  // Get the GPU freeze speed. Differs in that it is not the same size
  GPUtype cfreeze  = gm->gputype.getGPUtype(prhs[14]);
  double *gpu_cf = (double *)gm->gputype.getGPUptr(cfreeze);

  // Get flux factor (dt / d) amd determine if we are doing the hydro case or the MHD case
  double lambda   = *mxGetPr(prhs[15]);
  int isPureHydro = (int)*mxGetPr(prhs[16]);

  // Do the usual rigamarole determining array sizes and GPU launch dimensions
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

  // Invoke the kernel
  if(arraySize.x > 1) {
    if(isPureHydro) {
      cukern_TVDStep_hydro_uniform<<<gridsize, blocksize>>>(srcs[0], srcs[1], srcs[2], srcs[3], srcs[4], srcs[8], gpu_cf, srcs[9], srcs[10], srcs[11], srcs[12], srcs[13], lambda, arraySize.x);
    } else {
      cukern_TVDStep_mhd_uniform<<<gridsize, blocksize>>>(srcs[0], srcs[1], srcs[2], srcs[3], srcs[4], srcs[5], srcs[6], srcs[7], srcs[8], gpu_cf, srcs[9], srcs[10], srcs[11], srcs[12], srcs[13], lambda, arraySize.x);
    }
  }


}

/* blockidx.{xy} is our index in {yz}, and gridDim.{xy} gives the {yz} size */
/* Expect invocation with n+4 threads */
__global__ void cukern_TVDStep_mhd_uniform(double *rho, double *E, double *px, double *py, double *pz, double *bx, double *by, double *bz, double *P, double *Cfreeze, double *rhoW, double *enerW, double *pxW, double *pyW, double *pzW, double lambda, int nx)
{
double Cinv, rhoinv;
double q_i[5];
double b_i[3];
double w_i;
__shared__ double fluxLR[2][BLOCKLENP4];
__shared__ double derivLR[2][BLOCKLENP4];
double *fluxdest;

/* Step 0 - obligatory annoying setup stuff (ASS) */
int I0 = nx*(blockIdx.x + gridDim.x * blockIdx.y);
int Xindex = (threadIdx.x-2);
int Xtrack = Xindex;
Xindex += nx*(threadIdx.x < 2);

int x; /* = Xindex % nx; */
int i;
bool doIflux = (threadIdx.x > 1) && (threadIdx.x < BLOCKLEN+2);

/* Step 1 - calculate W values */
Cinv = 1.0/Cfreeze[blockIdx.x + gridDim.x * blockIdx.y];

while(Xtrack < nx+2) {
    x = I0 + (Xindex % nx);

    rhoinv = 1.0/rho[x]; /* Preload all these out here */
    q_i[0] = rho[x];
    q_i[1] = E[x];       /* So we avoid multiple loops */
    q_i[2] = px[x];      /* over them inside the flux loop */
    q_i[3] = py[x];
    q_i[4] = pz[x];
    b_i[0] = bx[x];
    b_i[1] = by[x];
    b_i[2] = bz[x];

    /* rho, E, px, py, pz going down */
    /* Iterate over variables to flux */
    for(i = 0; i < 5; i++) {
        /* Step 1 - Calculate raw fluxes */
        switch(i) {
            case 0: w_i = q_i[2] * Cinv; break;
            case 1: w_i = (q_i[2] * (q_i[1] + P[x]) - b_i[0]*(q_i[2]*b_i[0]+q_i[3]*b_i[1]+q_i[4]*b_i[2]) ) * (rhoinv*Cinv); break;
            case 2: w_i = (q_i[2]*q_i[2]*rhoinv + P[x] - b_i[0]*b_i[0])*Cinv; break;
            case 3: w_i = (q_i[2]*q_i[3]*rhoinv        - b_i[0]*b_i[1])*Cinv; break;
            case 4: w_i = (q_i[2]*q_i[4]*rhoinv        - b_i[0]*b_i[2])*Cinv; break;
            }

        /* Step 2 - Decouple to L/R flux */
        fluxLR[0][threadIdx.x] = 0.5*(q_i[i] - w_i); /* Left  going flux */
        fluxLR[1][threadIdx.x] = 0.5*(q_i[i] + w_i); /* Right going flux */
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
        __syncthreads();

        /* Step 4 - Perform flux and write to output array */
       if( doIflux && (Xindex < nx) ) {
            switch(i) {
                case 0: fluxdest = rhoW; break;
                case 1: fluxdest = enerW; break;
                case 2: fluxdest = pxW; break;
                case 3: fluxdest = pyW; break;
                case 4: fluxdest = pzW; break;
                }

            fluxdest[x] -= lambda * ( fluxLR[0][threadIdx.x] - fluxLR[0][threadIdx.x+1] + \
                                      fluxLR[1][threadIdx.x] - fluxLR[1][threadIdx.x-1]  ) / Cinv; 

            }

        __syncthreads();
        }

    Xindex += BLOCKLEN;
    Xtrack += BLOCKLEN;
    }

}


__global__ void cukern_TVDStep_hydro_uniform(double *rho, double *E, double *px, double *py, double *pz, double *P, double *Cfreeze, double *rhoW, double *enerW, double *pxW, double *pyW, double *pzW, double lambda, int nx)
{
double Cinv, rhoinv;
double q_i[5];
double w_i;
__shared__ double fluxLR[2][BLOCKLENP4];
__shared__ double derivLR[2][BLOCKLENP4];
double *fluxdest;

/* Step 0 - obligatory annoying setup stuff (ASS) */
int I0 = nx*(blockIdx.x + gridDim.x * blockIdx.y);
int Xindex = (threadIdx.x-2);
int Xtrack = Xindex;
Xindex += nx*(threadIdx.x < 2);

int x; /* = Xindex % nx; */
int i;
bool doIflux = (threadIdx.x > 1) && (threadIdx.x < BLOCKLEN+2);

/* Step 1 - calculate W values */
Cinv = 1.0/Cfreeze[blockIdx.x + gridDim.x * blockIdx.y];

while(Xtrack < nx+2) {
    x = I0 + (Xindex % nx);

    rhoinv = 1.0/rho[x]; /* Preload all these out here */
    q_i[0] = rho[x];
    q_i[1] = E[x];       /* So we avoid multiple loops */
    q_i[2] = px[x];      /* over them inside the flux loop */
    q_i[3] = py[x];
    q_i[4] = pz[x];

    /* rho, E, px, py, pz going down */
    /* Iterate over variables to flux */
    for(i = 0; i < 5; i++) {
        /* Step 1 - Calculate raw fluxes */
        switch(i) {
            case 0: w_i = q_i[2] * Cinv; break;
            case 1: w_i = (q_i[2] * (q_i[1] + P[x]) - (q_i[2]+q_i[3]+q_i[4]) ) * (rhoinv * Cinv); break;
            case 2: w_i = (q_i[2]*q_i[2]*rhoinv + P[x]) * Cinv; break;
            case 3: w_i = (q_i[2]*q_i[3]*rhoinv) * Cinv; break;
            case 4: w_i = (q_i[2]*q_i[4]*rhoinv) * Cinv; break;
            }

        /* Step 2 - Decouple to L/R flux */
        fluxLR[0][threadIdx.x] = 0.5*(q_i[i] - w_i); /* Left  going flux */
        fluxLR[1][threadIdx.x] = 0.5*(q_i[i] + w_i); /* Right going flux */
        __syncthreads();

        /* Step 3 - Differentiate fluxes & call limiter */
            /* left flux */
        derivLR[0][threadIdx.x] = fluxLR[0][(threadIdx.x-1)%BLOCKLENP4] - fluxLR[0][threadIdx.x]; /* left derivative */
        derivLR[1][threadIdx.x] = fluxLR[0][threadIdx.x] - fluxLR[0][(threadIdx.x+1)%BLOCKLENP4]; /* right derivative */
        __syncthreads();
//        cukern_FluxLimiter_VanLeer(derivLR, fluxLR, 0);
        __syncthreads();

            /* Right flux */
        derivLR[0][threadIdx.x] = fluxLR[1][threadIdx.x] - fluxLR[1][(threadIdx.x-1)%BLOCKLENP4]; /* left derivative */
        derivLR[1][threadIdx.x] = fluxLR[1][(threadIdx.x+1)%BLOCKLENP4] - fluxLR[1][threadIdx.x]; /* right derivative */
        __syncthreads();
//        cukern_FluxLimiter_VanLeer(derivLR, fluxLR, 1);
        __syncthreads();

        /* Step 4 - Perform flux and write to output array */
       if( doIflux && (Xindex < nx) ) {
            switch(i) {
                case 0: fluxdest = rhoW; break;
                case 1: fluxdest = enerW; break;
                case 2: fluxdest = pxW; break;
                case 3: fluxdest = pyW; break;
                case 4: fluxdest = pzW; break;
                }

            fluxdest[x] -= lambda * ( fluxLR[0][threadIdx.x] - fluxLR[0][threadIdx.x+1] + \
                                      fluxLR[1][threadIdx.x] - fluxLR[1][threadIdx.x-1]  ) / Cinv;
            }

        __syncthreads();
        }

    Xindex += BLOCKLEN;
    Xtrack += BLOCKLEN;
    }

}



// NOTE: This function divides the fluxes by two, as this is NOT done in the left/right decoupling phase
// in order to avoid any more multiplies by 1/2 than necessary.
__device__ void cukern_FluxLimiter_VanLeer(double deriv[2][BLOCKLENP4], double flux[2][BLOCKLENP4], int who)
{

double r;

r = deriv[0][threadIdx.x] * deriv[1][threadIdx.x];
if(r < 0.0) r = 0.0;

r = r / ( deriv[0][threadIdx.x] + deriv[1][threadIdx.x]);
if (isnan(r)) { r = 0.0; }

flux[who][threadIdx.x] += r;

}
