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

__global__ void cukern_magnetWstep_uniform(double *velGrid, double *mag, double *bW, double *velFlow, double lambda, int nx);

#define BLOCKLEN 48
#define BLOCKLENP2 50
#define BLOCKLENP4 52

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // At least 2 arguments expected
    // Input and result
    if ((nrhs!=3) || (nlhs != 2)) mexErrMsgTxt("Wrong number of arguments: need [magW,velFlow] = cudaMagWflux(mag, velgrid, lambda)\n");

    if (init == 0) {
        gm = gmGetGPUmat();
        init = 1;
    }

    // Get source array info and create destination arrays
    int numel;
    GPUtype srcReference = gm->gputype.getGPUtype(prhs[0]);

    double **srcs = getGPUSourcePointers(prhs, 2, &numel, 0, gm);
    double **dest = makeGPUDestinationArrays(srcReference,  plhs, 2, gm);

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
            hv = arraySize.x; hu = 1; hw = arraySize.x * arraySize.y;
            nu = arraySize.y; break;
        case 3: // Z direction flux: u = z, v = x, w = y;
            gridsize.x = arraySize.x;
            gridsize.y = arraySize.y;
            hu = arraySize.x * arraySize.y; hv = 1; hw = arraySize.x;
            nu = arraySize.z; break;
    }
    cukern_magnetWstep_uniform<<<gridsize , blocksize>>>(srcs[0], srcs[1], dest[0], dest[1], lambda, arraySize.x); 
}

__global__ void cukern_magnetWstep_uniform(double *velGrid, double *mag, double *bW, double *velFlow, double lambda, int nx)
{
    double v;
    double b;
    double bv;
    double locVelFlow;
    __shared__ double flux[BLOCKLENP4];

    /* Step 0 - obligatory annoying setup stuff (ASS) */
    int I0 = nx * (blockIdx.x + gridDim.x * blockIdx.y);
    int Xindex = (threadIdx.x-2);
    int Xtrack = Xindex;
    Xindex += nx * (threadIdx.x < 2);

    int x;
    bool doIflux = (threadIdx.x > 1) && (threadIdx.x < BLOCKLEN+2);

    while(Xtrack < nx + 2) {
        x = I0 + (Xindex % nx) ;

        v = velGrid[x];
        b = mag[x];

        // First step - calculate velocityflow
        flux[threadIdx.x] = v;
        locVelFlow = ((flux[threadIdx.x] + flux[(threadIdx.x + 1)%BLOCKLENP4]) < 0);
    
        // Second step - calculate flux
        bv = b * v;
        flux[threadIdx.x] = bv;
        __syncthreads();
        
        bv = bv * (1 - locVelFlow) + flux[(threadIdx.x + 1)%BLOCKLENP4] * locVelFlow;
        __syncthreads();
        
        flux[threadIdx.x] = bv;
        __syncthreads();

        // Third step - Perform flux and write to output array
        if( doIflux && (Xindex < nx) ) {
            bW[x] = b - lambda * ( bv - flux[threadIdx.x - 1] ); 
            velFlow[x] = locVelFlow;
        }

        Xindex += BLOCKLEN;
        Xtrack += BLOCKLEN;
        __syncthreads();
    }
}
