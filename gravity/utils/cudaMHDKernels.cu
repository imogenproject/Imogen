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

//double **getGPUSourcePointers(const mxArray *prhs[], int num, int *retNumel);
//double **makeGPUDestinationArrays(GPUtype src, mxArray *retArray[], int howmany);
//double *makeDestinationArray(GPUtype src, mxArray *retArray[]);

#define OP_SOUNDSPEED 1
#define OP_GASPRESSURE 2
#define OP_TOTALPRESSURE 3
#define OP_MAGPRESSURE 4
#define OP_TOTALANDSND 5
#define OP_WARRAYS 6
#define OP_RELAXINGFLUX 7
#define OP_SEPERATELRFLUX 8
__global__ void cukern_Soundspeed(double *rho, double *E, double *px, double *py, double *pz, double *bx, double *by, double *bz, double *dout, double gam, int n);
__global__ void cukern_GasPressure(double *rho, double *E, double *px, double *py, double *pz, double *bx, double *by, double *bz, double *dout, double gam, int n);
__global__ void cukern_TotalPressure(double *rho, double *E, double *px, double *py, double *pz, double *bx, double *by, double *bz, double *dout, double gam, int n);
__global__ void cukern_MagneticPressure(double *bx, double *by, double *bz, double *dout, int n);
__global__ void cukern_TotalAndSound(double *rho, double *E, double *px, double *py, double *pz, double *bx, double *by, double *bz, double *total, double *sound, double gam, int n);
__global__ void cukern_CalcWArrays(double *rho, double *E, double *px, double *py, double *pz, double *bx, double *by, double *bz, double *P, double *Cfreeze, double *rhoW, double *enerW, double *pxW, double *pyW, double *pzW, int dir, int n);

__global__ void cukern_SeperateLRFlux(double *arr, double *wArr, double *left, double *right, int n);
__global__ void cukern_PerformFlux(double *array0, double *Cfreeze, double *fluxRa, double *fluxRb, double *fluxLa, double *fluxLb, double *out, double lambda, int n);

#define BLOCKWIDTH 256
#define THREADLOOPS 1


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (init == 0) {
    // Initialize function
    // mexLock();
    // load GPUmat
    gm = gmGetGPUmat();
    init = 1;
  }

  // Determine appropriate number of arguments for RHS
  if (nrhs < 2) mexErrMsgTxt("Require at least (computation type, input argument)");
  int operation = (int)*mxGetPr(prhs[0]);

  dim3 blocksize; blocksize.x = BLOCKWIDTH; blocksize.y = blocksize.z = 1;
  int numel; dim3 gridsize;

  // Select the appropriate kernel to invoke
  if((operation == OP_SOUNDSPEED) || (operation == OP_GASPRESSURE) || (operation == OP_TOTALPRESSURE)) {
    if( (nlhs != 1) || (nrhs != 10)) { mexErrMsgTxt("Soundspeed operator is Cs = cudaMHDKernels(1, rho, E, px, py, pz, bx, by, bz, gamma)"); }
    double gam = *mxGetPr(prhs[9]);
    double **srcs = getGPUSourcePointers(prhs, 8, &numel, 1, gm);
    gridsize.x = numel / (BLOCKWIDTH*THREADLOOPS); if(gridsize.x * (BLOCKWIDTH*THREADLOOPS) < numel) gridsize.x++;
    gridsize.y = gridsize.z =1;
    double **destPtr = makeGPUDestinationArrays(gm->gputype.getGPUtype(prhs[1]), plhs, 1, gm);
//printf("%i %i %i %i %i %i\n", blocksize.x, blocksize.y, blocksize.z, gridsize.x, gridsize.y, gridsize.z);
    switch(operation) {
      case OP_SOUNDSPEED:       cukern_Soundspeed<<<gridsize, blocksize>>>(srcs[0], srcs[1], srcs[2], srcs[3], srcs[4], srcs[5], srcs[6], srcs[7], destPtr[0], gam, numel); break;
      case OP_GASPRESSURE:     cukern_GasPressure<<<gridsize, blocksize>>>(srcs[0], srcs[1], srcs[2], srcs[3], srcs[4], srcs[5], srcs[6], srcs[7], destPtr[0], gam, numel); break;
      case OP_TOTALPRESSURE: cukern_TotalPressure<<<gridsize, blocksize>>>(srcs[0], srcs[1], srcs[2], srcs[3], srcs[4], srcs[5], srcs[6], srcs[7], destPtr[0], gam, numel); break;
    }
    free(destPtr);

  } else if((operation == OP_MAGPRESSURE)) {
    if( (nlhs != 1) || (nrhs != 4)) { mexErrMsgTxt("Magnetic pressure operator is Pm = cudaMHDKernels(4, bx, by, bz)"); }

    double **srcs = getGPUSourcePointers(prhs, 3, &numel, 1, gm);
    gridsize.x = numel / (BLOCKWIDTH*THREADLOOPS); if(gridsize.x * (BLOCKWIDTH*THREADLOOPS) < numel) gridsize.x++;
    gridsize.y = gridsize.z =1;
    double **destPtr = makeGPUDestinationArrays(gm->gputype.getGPUtype(prhs[1]), plhs, 1, gm);

    cukern_MagneticPressure<<<gridsize, blocksize>>>(srcs[0], srcs[1], srcs[2], destPtr[0], numel);
    free(destPtr); free(srcs);

  } else if((operation == OP_TOTALANDSND)) {
    if( (nlhs != 2) || (nrhs != 10)) { mexErrMsgTxt("Soundspeed operator is [Ptot Cs] = cudaMHDKernels(5, rho, E, px, py, pz, bx, by, bz, gamma)"); }
    double gam = *mxGetPr(prhs[9]);
    double **srcs = getGPUSourcePointers(prhs, 8, &numel, 1, gm);
    gridsize.x = numel / (BLOCKWIDTH*THREADLOOPS); if(gridsize.x * (BLOCKWIDTH*THREADLOOPS) < numel) gridsize.x++;
    gridsize.y = gridsize.z =1;
    double **destPtr = makeGPUDestinationArrays(gm->gputype.getGPUtype(prhs[1]), plhs, 2, gm);

    cukern_TotalAndSound<<<gridsize, blocksize>>>(srcs[0], srcs[1], srcs[2], srcs[3], srcs[4], srcs[5], srcs[6], srcs[7], destPtr[0], destPtr[1], gam, numel);
    free(destPtr); free(srcs);
  } else if ((operation == OP_WARRAYS)) {
    if( (nlhs != 5) || (nrhs != 12)) { mexErrMsgTxt("solving W operator is [rhoW enerW pxW pyW pzW] = cudaMHDKernels(6, rho, E, px, py, pz, bx, by, bz, P, cFreeze, direction)"); }
    int dir = (int)*mxGetPr(prhs[11]);
    double **srcs = getGPUSourcePointers(prhs, 10, &numel, 1, gm);
    gridsize.x = numel / (BLOCKWIDTH*THREADLOOPS); if(gridsize.x * (BLOCKWIDTH*THREADLOOPS) < numel) gridsize.x++;
    gridsize.y = gridsize.z =1;
    double **destPtr = makeGPUDestinationArrays(gm->gputype.getGPUtype(prhs[1]), plhs, 5, gm);

    cukern_CalcWArrays<<<gridsize, blocksize>>>(srcs[0], srcs[1], srcs[2], srcs[3], srcs[4], srcs[5], srcs[6], srcs[7], srcs[8], srcs[9], destPtr[0], destPtr[1], destPtr[2], destPtr[3], destPtr[4], dir, numel);
    free(destPtr); free(srcs);
  } else if ((operation == OP_RELAXINGFLUX)) {
    if( (nlhs != 1) || (nrhs != 8)) { mexErrMsgTxt("relaxing flux operator is fluxed = cudaMHDKernels(7, old, tempfreeze, right, right_shifted, left, left_shifted, lambda)"); }
    double lambda = *mxGetPr(prhs[7]);
    double **srcs = getGPUSourcePointers(prhs, 6, &numel, 1, gm);
    gridsize.x = numel / (BLOCKWIDTH*THREADLOOPS); if(gridsize.x * (BLOCKWIDTH*THREADLOOPS) < numel) gridsize.x++;
    gridsize.y = gridsize.z =1;
    double **destPtr = makeGPUDestinationArrays(gm->gputype.getGPUtype(prhs[1]), plhs, 1, gm);

    cukern_PerformFlux<<<gridsize, blocksize>>>(srcs[0], srcs[1], srcs[2], srcs[3], srcs[4], srcs[5], destPtr[0], lambda, numel);
    free(destPtr); free(srcs);
  } else if ((operation == OP_SEPERATELRFLUX)) {
    if ((nlhs != 2) || (nrhs != 3)) { mexErrMsgTxt("flux seperation operator is [Fl Fr] = cudaMHDKernels(8, array, wArray)"); }
    double **srcs = getGPUSourcePointers(prhs, 2, &numel, 1, gm);
    gridsize.x = numel / (BLOCKWIDTH*THREADLOOPS); if(gridsize.x * (BLOCKWIDTH*THREADLOOPS) < numel) gridsize.x++;
    gridsize.y = gridsize.z =1;
    double **destPtr = makeGPUDestinationArrays(gm->gputype.getGPUtype(prhs[1]), plhs, 2, gm);

    cukern_SeperateLRFlux<<<gridsize, blocksize>>>(srcs[0], srcs[1], destPtr[0], destPtr[1], numel);
    free(destPtr); free(srcs);
  }

}

//#define KERNEL_PREAMBLE int x = THREADLOOPS*(threadIdx.x + blockDim.x*blockIdx.x); if (x >= n) {return;} int imax; ((x+THREADLOOPS) > n) ? imax = n : imax = x + THREADLOOPS; for(; x < imax; x++)
#define KERNEL_PREAMBLE int x = threadIdx.x + blockDim.x*blockIdx.x; if (x >= n) { return; }

// THIS KERNEL CALCULATES SOUNDSPEED 
__global__ void cukern_Soundspeed(double *rho, double *E, double *px, double *py, double *pz, double *bx, double *by, double *bz, double *dout, double gam, int n)
{
double gg1 = gam*(gam-1.0);

KERNEL_PREAMBLE
dout[x] = sqrt(abs( (gg1*(E[x] - .5*(px[x]*px[x] + py[x]*py[x] + pz[x]*pz[x])/rho[x]) + (2.0 -.5*gg1)*(bx[x]*bx[x] + by[x]*by[x] + bz[x]*bz[x]))/rho[x] ));
}

// THIS KERNEL CALCULATES GAS PRESSURE
__global__ void cukern_GasPressure(double *rho, double *E, double *px, double *py, double *pz, double *bx, double *by, double *bz, double *dout, double gam, int n)
{
KERNEL_PREAMBLE
dout[x] = (gam-1.0)*abs(E[x] - .5*((px[x]*px[x]+py[x]*py[x]+pz[x]*pz[x])/rho[x] + bx[x]*bx[x]+by[x]*by[x]+bz[x]*bz[x]));
}

// THIS KERNEL CALCULATES TOTAL PRESSURE
__global__ void cukern_TotalPressure(double *rho, double *E, double *px, double *py, double *pz, double *bx, double *by, double *bz, double *dout, double gam, int n)
{
KERNEL_PREAMBLE
dout[x] = (gam-1.0)*abs(E[x] - .5*((px[x]*px[x]+py[x]*py[x]+pz[x]*pz[x])/rho[x])) + .5*(2.0-gam)*(bx[x]*bx[x]+by[x]*by[x]+bz[x]*bz[x]);
}

// THIS KERNEL CALCULATES MAGNETIC PRESSURE
__global__ void cukern_MagneticPressure(double *bx, double *by, double *bz, double *dout, int n)
{
KERNEL_PREAMBLE
dout[x] = .5*(bx[x]*bx[x]+by[x]*by[x]+bz[x]*bz[x]);
}

__global__ void cukern_TotalAndSound(double *rho, double *E, double *px, double *py, double *pz, double *bx, double *by, double *bz, double *total, double *sound, double gam, int n)
{
double gg1 = gam*(gam-1.0);
double psqhf, bsqhf;

KERNEL_PREAMBLE {
	psqhf = .5*(px[x]*px[x]+py[x]*py[x]+pz[x]*pz[x]);
	bsqhf = .5*(bx[x]*bx[x]+by[x]*by[x]+bz[x]*bz[x]);
	
	total[x] = (gam-1.0)*abs(E[x] - psqhf/rho[x]) + (2.0-gam)*bsqhf;
	sound[x]   = sqrt(abs( (gg1*(E[x] - psqhf/rho[x]) + (4.0 - gg1)*bsqhf)/rho[x] ));
	}
}

__global__ void cukern_CalcWArrays(double *rho, double *E, double *px, double *py, double *pz, double *bx, double *by, double *bz, double *P, double *Cfreeze, double *rhoW, double *enerW, double *pxW, double *pyW, double *pzW, int dir, int n)
{
double Cinv, rhoinv;

KERNEL_PREAMBLE {

Cinv = 1.0/Cfreeze[x];
rhoinv = 1.0/rho[x];

switch(dir) {
  case 1:
    rhoW[x]  = px[x] * Cinv;
    enerW[x] = (px[x] * (E[x] + P[x]) - bx[x]*(px[x]*bx[x]+py[x]*by[x]+pz[x]*bz[x]) ) * (rhoinv*Cinv);
    pxW[x]   = (px[x]*px[x]*rhoinv + P[x] - bx[x]*bx[x])*Cinv;
    pyW[x]   = (px[x]*py[x]*rhoinv        - bx[x]*by[x])*Cinv;
    pzW[x]   = (px[x]*pz[x]*rhoinv        - bx[x]*bz[x])*Cinv;
    break;
  case 2:
    rhoW[x]  = py[x] * Cinv;
    enerW[x] = (py[x] * (E[x] + P[x]) - by[x]*(px[x]*bx[x]+py[x]*by[x]+pz[x]*bz[x]) ) * (rhoinv*Cinv);
    pxW[x]   = (py[x]*px[x]*rhoinv        - by[x]*bx[x])*Cinv;
    pyW[x]   = (py[x]*py[x]*rhoinv + P[x] - by[x]*by[x])*Cinv;
    pzW[x]   = (py[x]*pz[x]*rhoinv        - by[x]*bz[x])*Cinv;
    break;
  case 3:
    rhoW[x]  = pz[x] * Cinv;
    enerW[x] = (pz[x] * (E[x] + P[x]) - bz[x]*(px[x]*bx[x]+py[x]*by[x]+pz[x]*bz[x]) ) * (rhoinv*Cinv);
    pxW[x]   = (pz[x]*px[x]*rhoinv        - bz[x]*bx[x])*Cinv;
    pyW[x]   = (pz[x]*py[x]*rhoinv        - bz[x]*by[x])*Cinv;
    pzW[x]   = (pz[x]*pz[x]*rhoinv + P[x] - bz[x]*bz[x])*Cinv;
    break;
  }

}
/*  mass.wArray    = mom(X).array ./ freezeSpd.array;

    %--- ENERGY DENSITY ---%
    ener.wArray    = velocity .* (ener.array + press) - mag(X).cellMag.array .* ...
                        ( mag(1).cellMag.array .* mom(1).array ...
                        + mag(2).cellMag.array .* mom(2).array ...
                        + mag(3).cellMag.array .* mom(3).array) ./ mass.array;
    ener.wArray    = ener.wArray ./ freezeSpd.array;

    %--- MOMENTUM DENSITY ---%
    for i=1:3
        mom(i).wArray    = (velocity .* mom(i).array + press*dirVec(i)...
                             - mag(X).cellMag.array .* mag(i).cellMag.array) ./ freezeSpd.array;
    end*/

}

__global__ void cukern_PerformFlux(double *array0, double *Cfreeze, double *fluxRa, double *fluxRb, double *fluxLa, double *fluxLb, double *out, double lambda, int n)
{
KERNEL_PREAMBLE 
out[x] = array0[x] - lambda*Cfreeze[x]*(fluxRa[x] - fluxRb[x] + fluxLa[x] - fluxLb[x]);

//v(i).store.array = v(i).array - 0.5*fluxFactor .* tempFreeze .* ...
//                        ( v(i).store.fluxR.array - v(i).store.fluxR.shift(X,-1) ...
//                        + v(i).store.fluxL.array - v(i).store.fluxL.shift(X,1) );
}

__global__ void cukern_SeperateLRFlux(double *arr, double *wArr, double *left, double *right, int n)
{
KERNEL_PREAMBLE {
	left[x]  = .5*(arr[x] - wArr[x]);
	right[x] = .5*(arr[x] + wArr[x]);
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
__shared__ derivLR[2][36];
double *fluxdest;

/* Step 0 - obligatory annoying setup stuff (ASS) */
int I0 = nx*(blockIdx.x + gridDim.x * blockIdx.y);
int Xindex = (threadIdx.x-2);
int x; /* = Xindex % nx; */
int i;

/* Step 1 - calculate W values */
Cinv = 1.0/Cfreeze[I0];

    while(Xindex < nx+2) {

    x = Xindex % nx

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
        fluxLR[0][threadIdx.x] = 0.5*(q_i[i] - w_i);
        fluxLR[1][threadIdx.x] = 0.5*(q_i[i] + w_i); 

        /* Step 3 - Differentiate fluxes & call limiter */
        derivLR[0][threadIdx.x] = fluxLR[0][threadIdx.x]        - fluxLR[0][(threadIdx.x-1)%36];
        derivLR[1][threadIdx.x] = fluxLR[0][(threadIdx.x+1)%36] - fluxLR[0][threadIdx.x];
        __syncthreads();
        fluxLimiter(derivLR, fluxLR, 0);

        derivLR[0][threadIdx.x] = fluxLR[1][threadIdx.x]        - fluxLR[1][(threadIdx.x-1)%36];
        derivLR[1][threadIdx.x] = fluxLR[1][(threadIdx.x+1)%36] - fluxLR[1][threadIdx.x];
        __syncthreads();
        fluxLimiter(derivLR, fluxLR, 1);

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

            fluxdest[I0+x] -= lambda * ( fluxLR[0][threadIdx.x] - fluxLR[0][threadIdx.x+1] + fluxLR[1][threadIdx.x] - fluxLR[1][threadIdx.x-1]  ) / Cinv;
            }
        }

    Xindex += 32;
    }

}
