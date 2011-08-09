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

__global__ void cukern_Wstep_mhd_uniform(double *rho, double *E, double *px, double *py, double *pz, double *bx, double *by, double *bz, double *P, double *Cfreeze, double *rhoW, double *enerW, double *pxW, double *pyW, double *pzW, double lambda, int nx);
__global__ void cukern_Wstep_hydro_uniform(double *rho, double *E, double *px, double *py, double *pz, double *P, double *Cfreeze, double *rhoW, double *enerW, double *pxW, double *pyW, double *pzW, double lambda, int nx);
__global__ void nullStep(fluidVarPtrs fluid, int numel);


#define BLOCKLEN 48
#define BLOCKLENP2 50
#define BLOCKLENP4 52

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // At least 2 arguments expected
  // Input and result
  if ((nrhs!=3) || (nlhs != 1)) mexErrMsgTxt("Wrong number of arguments: need [magW] = cudaMagWflux(mag, velgrid, lambda)\n");

  if (init == 0) {
    gm = gmGetGPUmat();
    init = 1;
  }

  // Get source array info and create destination arrays
  int numel;
  GPUtype srcReference = gm->gputype.getGPUtype(prhs[0]);

  double **srcs = getGPUSourcePointers(prhs, 2, &numel, 0, gm);
  double **dest = makeGPUDestinationArrays(srcReference,  plhs, 1, gm);

  // Establish launch dimensions & a few other parameters
  int fluxDirection = 1;
  double lambda     = *mxGetPr(prhs[2]);
  int numDims        = gm->gputype.getNdims(srcReference);
  const int *dims    = gm->gputype.getSize(srcReference);

  dim3 arraySize;
  arraySize.x = dims[0];
  numDims > 1 ? arraySize.y = dims[1] : arraySize.y = 1;
  numDims > 2 ? arraySize.z = dims[2] : arraySize.z = 1;

  dim3 blocksize, gridsize;
  int hu, hv, hw, nu;

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
//      hu = arraySize.x; hv = 1; hw = arraySize.x * arraySize.y;
      hv = arraySize.x; hu = 1; hw = arraySize.x * arraySize.y;
      nu = arraySize.y; break;
    case 3: // Z direction flux: u = z, v = x, w = y;
      gridsize.x = arraySize.x;
      gridsize.y = arraySize.y;
      hu = arraySize.x * arraySize.y; hv = 1; hw = arraySize.x;
      nu = arraySize.z; break;
    }

if(nu > 1) {
    cukern_magnetWstep_unifirm<<<gridsize, blocksize>>>(srcs[0], srcs[1], dest[0], lambda, arraySize.x);
  } else {
  nullStep<<<32, 128>>>(fluid, numel);
  }

}

__global__ void cukern_magnetWstep_uniform(double *velGrid, double *mag, double *bW, double lambda, int nx)
{
double v_i;
double b_i;
__shared__ double  vel[BLOCKLENP4];
__shared__ double flux[BLOCKLENP4];

/* Step 0 - obligatory annoying setup stuff (ASS) */
int I0 = nx*(blockIdx.x + gridDim.x * blockIdx.y);
int Xindex = (threadIdx.x-2);
int Xtrack = Xindex;
Xindex += nx*(threadIdx.x < 2);

int x; /* = Xindex % nx; */
int i;
bool doIflux = (threadIdx.x > 1) && (threadIdx.x < BLOCKLEN+2);

/* Step 1 - calculate W values */

while(Xtrack < nx+2) {
    x = I0 + (Xindex % nx);

/*

    %-----------------------------------------------------------------------------------------------
    % Initialization
    %---------------
    fluxFactor = 0.5*run.time.dTime ./ run.DGRID{X};
    if isa(velGrid.array, 'GPUdouble')
        velocityFlow = double(velGrid.array + velGrid.shift(X,1));
        velocityFlow = GPUdouble( velocityFlow < 0.0 );
    else
        velocityFlow = ( (velGrid.array + velGrid.shift(X,1)) < 0.0 );
    end
    
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    %Half-Timestep predictor step (first-order upwind,not TVD)
    %+++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
    mag(I).store(X).fluxR.array = mag(I).array .* velGrid.array;
    mag(I).store(X).fluxR.array = mag(I).store(X).fluxR.array .* (1-velocityFlow) ...
							  + mag(I).store(X).fluxR.shift(X,1) .* velocityFlow;
    
    mag(I).store(X).array = mag(I).array ...
				- fluxFactor .* (mag(I).store(X).fluxR.array - mag(I).store(X).fluxR.shift(X,-1));
*/

    vel[x] = 
    rhoinv = 1.0/rho[x]; /* Preload all these out here */
    v_i[0] = px[x]*rhoinv;
    v_i[1] = py[x]*rhoinv;       /* So we avoid multiple loops */
    v_i[2] = px[x]*rhoinv;      /* over them inside the flux loop */
    b_i[0] = bx[x];
    b_i[1] = by[x];
    b_i[2] = bz[x];

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

