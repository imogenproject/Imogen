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
#include "GPUmat.hh" // GPUmat drop-in toolkit

// static paramaters
static int init = 0;
static GPUmat *gm;

#include "cudaCommon.h" // This defines the getGPUSourcePointers and makeGPUDestinationArrays utility functions

__global__ void cukern_Wstep_mhd_uniform(double *rho, double *E, double *px, double *py, double *pz, double *bx, double *by, double *bz, double *P, double *Cfreeze, double *rhoW, double *enerW, double *pxW, double *pyW, double *pzW, double lambda, int nx);
__global__ void cukern_Wstep_hydro_uniform(double *rho, double *E, double *px, double *py, double *pz, double *P, double *Cfreeze, double *rhoW, double *enerW, double *pxW, double *pyW, double *pzW, double lambda, int nx);
__global__ void nullStep(fluidVarPtrs fluid, int numel);

#define BLOCKLEN 48
#define BLOCKLENP2 50
#define BLOCKLENP4 52

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // Check input and output argument counts

  // Input and result
  if ((nrhs!=12) || (nlhs != 5)) mexErrMsgTxt("Wrong number of arguments: need [5] = cudaWflux(rho, E, px, py, pz, bx, by, bz, Ptot, c_f, lambda, purehydro?)\n");

  if (init == 0) {
    gm = gmGetGPUmat();
    init = 1;
  }

  // Get source array info and create destination arrays
  int numel;
  GPUtype srcReference = gm->gputype.getGPUtype(prhs[0]);

  double **srcs = getGPUSourcePointers(prhs, 10, &numel, 0, gm);
  double **dest = makeGPUDestinationArrays(srcReference,  plhs, 5, gm);

  // Establish launch dimensions & a few other parameters
  int fluxDirection = 1;
  double lambda     = *mxGetPr(prhs[10]);
  int numDims        = gm->gputype.getNdims(srcReference);
  const int *dims    = gm->gputype.getSize(srcReference);

  dim3 arraySize;
  arraySize.x = dims[0];
  numDims > 1 ? arraySize.y = dims[1] : arraySize.y = 1;
  numDims > 2 ? arraySize.z = dims[2] : arraySize.z = 1;

  dim3 blocksize, gridsize;
  int hu, hv, hw, nu;

  // This bit is actually redundant now since arrays are always rotated so the fluid step is finite-differenced in the X direction
  blocksize.x = BLOCKLEN+4; blocksize.y = blocksize.z = 1;
  switch(fluxDirection) {
    case 1: // X direction flux: u = x, v = y, w = z;
      gridsize.x = arraySize.y;
      gridsize.y = arraySize.z;
      hu = 1; hv = arraySize.x; hw = arraySize.x * arraySize.x;
      nu = arraySize.x; break;
    case 2: // Y direction flux: u = y, v = x, w = z
      gridsize.x = arraySize.x;
      gridsize.y = arraySize.z;
      //hu = arraySize.x; hv = 1; hw = arraySize.x * arraySize.y;
      hv = arraySize.x; hu = 1; hw = arraySize.x * arraySize.y;
      nu = arraySize.y; break;
    case 3: // Z direction flux: u = z, v = x, w = y;
      gridsize.x = arraySize.x;
      gridsize.y = arraySize.y;
      hu = arraySize.x * arraySize.y; hv = 1; hw = arraySize.x;
      nu = arraySize.z; break;
    }

// It appears this is only used in the null step. It was used in a previous W step but that kernel was irreperably broken.
fluidVarPtrs fluid;
int i;
for(i = 0; i < 5; i++) { fluid.fluidIn[i] = srcs[i]; fluid.fluidOut[i] = dest[i]; }
fluid.B[0] = srcs[5];
fluid.B[1] = srcs[6];
fluid.B[2] = srcs[7];

fluid.Ptotal = srcs[8];
fluid.cFreeze = srcs[9];

// If the dimension has finite extent, performs actual step; If not, blits input arrays to output arrays
// NOTE: this situation should not occur, since the flux routine itself skips singleton dimensions for 1- and 2-d sims.
if(nu > 1) {
  int hydroOnly = (int)*mxGetPr(prhs[11]);
  
  if(hydroOnly == 1) {
    cukern_Wstep_hydro_uniform<<<gridsize, blocksize>>>(srcs[0], srcs[1], srcs[2], srcs[3], srcs[4], srcs[8], srcs[9], dest[0], dest[1], dest[2], dest[3], dest[4], lambda, arraySize.x);
    } else {
    cukern_Wstep_mhd_uniform<<<gridsize, blocksize>>>(srcs[0], srcs[1], srcs[2], srcs[3], srcs[4], srcs[5], srcs[6], srcs[7], srcs[8], srcs[9], dest[0], dest[1], dest[2], dest[3], dest[4], lambda, arraySize.x);
    }
  } else {
  nullStep<<<32, 128>>>(fluid, numel);
  }

}

__global__ void cukern_Wstep_mhd_uniform(double *rho, double *E, double *px, double *py, double *pz, double *bx, double *by, double *bz, double *P, double *Cfreeze, double *rhoW, double *enerW, double *pxW, double *pyW, double *pzW, double lambda, int nx)
{
double Cinv, rhoinv;
double q_i[5];
double b_i[3];
double w_i;
__shared__ double fluxLR[2][BLOCKLENP4];
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
        switch(i) {
            case 0: w_i = q_i[2] * Cinv; break;
            case 1: w_i = (q_i[2] * (q_i[1] + P[x]) - b_i[0]*(q_i[2]*b_i[0]+q_i[3]*b_i[1]+q_i[4]*b_i[2]) ) * (rhoinv*Cinv); break;
            case 2: w_i = (q_i[2]*q_i[2]*rhoinv + P[x] - b_i[0]*b_i[0])*Cinv; break;
            case 3: w_i = (q_i[2]*q_i[3]*rhoinv        - b_i[0]*b_i[1])*Cinv; break;
            case 4: w_i = (q_i[2]*q_i[4]*rhoinv        - b_i[0]*b_i[2])*Cinv; break;
            }

        /* Step 2 - decouple to L/R flux */
        fluxLR[0][threadIdx.x] = 0.5*(q_i[i] - w_i); /* Left  going flux */
        fluxLR[1][threadIdx.x] = 0.5*(q_i[i] + w_i); /* Right going flux */
        __syncthreads();

        /* Step 4 - Perform flux and write to output array */
        __syncthreads();
       if( doIflux && (Xindex < nx) ) {
            switch(i) {
                case 0: fluxdest = rhoW; break;
                case 1: fluxdest = enerW; break;
                case 2: fluxdest = pxW; break;
                case 3: fluxdest = pyW; break;
                case 4: fluxdest = pzW; break;
                }

            fluxdest[x] = q_i[i] - 0.5 * lambda * ( fluxLR[0][threadIdx.x] - fluxLR[0][threadIdx.x+1] + \
                                      fluxLR[1][threadIdx.x] - fluxLR[1][threadIdx.x-1]  ) / Cinv; 

            }

        __syncthreads();
        }

    Xindex += BLOCKLEN;
    Xtrack += BLOCKLEN;
    __syncthreads();
    }

}


__global__ void cukern_Wstep_hydro_uniform(double *rho, double *E, double *px, double *py, double *pz, double *P, double *Cfreeze, double *rhoW, double *enerW, double *pxW, double *pyW, double *pzW, double lambda, int nx)
{
double Cinv, rhoinv;
double q_i[5];
double w_i;
__shared__ double fluxLR[2][BLOCKLENP4];
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
        switch(i) {
            case 0: w_i = q_i[2] * Cinv; break;
            case 1: w_i = (q_i[2] * (q_i[1] + P[x])) * (rhoinv*Cinv); break;
            case 2: w_i = (q_i[2]*q_i[2]*rhoinv + P[x])*Cinv; break;
            case 3: w_i = (q_i[2]*q_i[3]*rhoinv       )*Cinv; break;
            case 4: w_i = (q_i[2]*q_i[4]*rhoinv       )*Cinv; break;
            }

        /* Step 2 - decouple to L/R flux */
        fluxLR[0][threadIdx.x] = 0.5*(q_i[i] - w_i); /* Left  going flux */
        fluxLR[1][threadIdx.x] = 0.5*(q_i[i] + w_i); /* Right going flux */
        __syncthreads();

        /* Step 4 - Perform flux and write to output array */
        __syncthreads();
       if( doIflux && (Xindex < nx) ) {
            switch(i) {
                case 0: fluxdest = rhoW; break;
                case 1: fluxdest = enerW; break;
                case 2: fluxdest = pxW; break;
                case 3: fluxdest = pyW; break;
                case 4: fluxdest = pzW; break;
                }

            fluxdest[x] = q_i[i] - 0.5 * lambda * ( fluxLR[0][threadIdx.x] - fluxLR[0][threadIdx.x+1] + \
                                      fluxLR[1][threadIdx.x] - fluxLR[1][threadIdx.x-1]  ) / Cinv; 

            }

        __syncthreads();
        }

    Xindex += BLOCKLEN;
    Xtrack += BLOCKLEN;
    __syncthreads();
    }

}

// Function simply blits fluid input variables to output, since with only 1 plane there's no derivative possible.
__global__ void nullStep(fluidVarPtrs fluid, int numel)
{
int idx0 = threadIdx.x + blockIdx.x*blockDim.x;
int didx = blockDim.x * gridDim.x;
int i;

while(idx0 < numel) {
  for(i = 0; i < 5; i++) { fluid.fluidOut[i][idx0] = fluid.fluidIn[i][idx0]; }

  idx0 += didx;
  }

}

