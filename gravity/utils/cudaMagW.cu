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

__global__ void cukern_magnetWstep_uniformX(double *mag, double *velGrid, double *bW, double *velFlow, double lambda, int nx);
__global__ void cukern_magnetWstep_uniformY(double *mag, double *velGrid, double *bW, double *velFlow, double lambda, int3 dims);
__global__ void cukern_magnetWstep_uniformZ(double *mag, double *velGrid, double *bW, double *velFlow, double lambda, int3 dims);

#define BLOCKDIMA 18
#define BLOCKDIMAM2 16
#define BLOCKDIMB 8

#define BLOCKLEN 128
#define BLOCKLENP4 132

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // At least 2 arguments expected
    // Input and result
    if ((nrhs!=4) || (nlhs != 2)) mexErrMsgTxt("Wrong number of arguments: need [magW,velFlow] = cudaMagWflux(mag, velgrid, lambda, dir)\n");

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
    int fluxDirection = (int)*mxGetPr(prhs[3]);
    double lambda     = *mxGetPr(prhs[2]);

    int numDims        = gm->gputype.getNdims(srcReference);
    const int *dims    = gm->gputype.getSize(srcReference);

    int3 arraySize;
    arraySize.x = dims[0];
    numDims > 1 ? arraySize.y = dims[1] : arraySize.y = 1;
    numDims > 2 ? arraySize.z = dims[2] : arraySize.z = 1;

    dim3 blocksize, gridsize;
    switch(fluxDirection) {
        case 1: // X direction flux. This is "priveleged" in that the shift and natural memory load directions align
            blocksize.x = BLOCKLEN+4; blocksize.y = blocksize.z = 1;
            gridsize.x = arraySize.y;
            gridsize.y = arraySize.z;
            cukern_magnetWstep_uniformX<<<gridsize , blocksize>>>(srcs[0], srcs[1], dest[0], dest[1], lambda, arraySize.x);
            break;
        case 2: // Y direction flux: u = y, v = x, w = z
            blocksize.x = BLOCKDIMB; blocksize.y = BLOCKDIMAM2;

            gridsize.x = arraySize.x / blocksize.x; gridsize.x += 1*(blocksize.x*gridsize.x < arraySize.x);
            gridsize.y = arraySize.y / blocksize.y; gridsize.y += 1*(blocksize.y*gridsize.y < arraySize.y);

            blocksize.y = BLOCKDIMA;

            cukern_magnetWstep_uniformY<<<gridsize , blocksize>>>(srcs[0], srcs[1], dest[0], dest[1], lambda, arraySize);
            break;
        case 3: // Z direction flux: u = z, v = x, w = y;
            blocksize.x = BLOCKDIMB; blocksize.y = BLOCKDIMAM2;

            gridsize.x = arraySize.x / blocksize.x; gridsize.x += 1*(blocksize.x*gridsize.x < arraySize.x);
            gridsize.y = arraySize.z / blocksize.y; gridsize.y += 1*(blocksize.y*gridsize.y < arraySize.z);

            blocksize.y = BLOCKDIMA;

            cukern_magnetWstep_uniformZ<<<gridsize , blocksize>>>(srcs[0], srcs[1], dest[0], dest[1], lambda, arraySize);

            break;
    }

}

__global__ void cukern_magnetWstep_uniformX(double *mag, double *velGrid, double *bW, double *velFlow, double lambda, int nx)
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
    __syncthreads();

    locVelFlow = (flux[threadIdx.x] + flux[(threadIdx.x + 1)%BLOCKLENP4]);
    if(locVelFlow < 0.0) { locVelFlow = 1.0; } else { locVelFlow = 0.0; }
    __syncthreads();

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

__global__ void cukern_magnetWstep_uniformY(double *mag, double *velGrid, double *bW, double *velFlow, double lambda, int3 dims)
{
double v, b, locVelFlow;

__shared__ double tile[BLOCKDIMB][BLOCKDIMA];
__shared__ double flux[BLOCKDIMB][BLOCKDIMA];

// Dimensions into the array
int myx = blockIdx.x*BLOCKDIMB + threadIdx.x;
int myy = blockIdx.y*BLOCKDIMAM2 + threadIdx.y - 1;

if((myx >= dims.x) || (myy > dims.y)) return; // we keep an extra Y thread for the finite diff.

bool IWrite = (threadIdx.y > 0) && (threadIdx.y <= BLOCKDIMAM2) && (myy < dims.y) && (myy >= 0);
// Exclude threads at the boundary of the fluxing direction from writing back

if(myy < 0) myy += dims.y; // wrap left edge back to right edge
myy = myy % dims.y; // wrap right edge back to left

int x = myx + dims.x*myy;
int z;

for(z = 0; z < dims.z; z++) {
    v = velGrid[x];
    b = mag[x];

    // first calculate velocityFlow
    tile[threadIdx.x][threadIdx.y] = v;
    flux[threadIdx.x][threadIdx.y] = b*v;
    __syncthreads();

    locVelFlow = (tile[threadIdx.x][threadIdx.y] + tile[threadIdx.x][(threadIdx.y+1) % BLOCKDIMA]);
    if(locVelFlow < 0.0) { locVelFlow = 1.0; } else { locVelFlow = 0.0; }

    __syncthreads();

    // Second step - calculate flux
    if(locVelFlow == 1) { tile[threadIdx.x][threadIdx.y] = flux[threadIdx.x][(threadIdx.y + 1)%BLOCKDIMA]; } else 
                        { tile[threadIdx.x][threadIdx.y] = flux[threadIdx.x][threadIdx.y]; }
   
    __syncthreads();

    // Third step - Perform flux and write to output array
    if( IWrite ) {
            bW[x] = b - lambda * ( tile[threadIdx.x][threadIdx.y] - tile[threadIdx.x][threadIdx.y-1]);
            velFlow[x] = locVelFlow;
        }

    x += dims.x*dims.y;
    __syncthreads(); 
    }

}

__global__ void cukern_magnetWstep_uniformZ(double *mag, double *velGrid, double *bW, double *velFlow, double lambda, int3 dims)
{
double v, b, locVelFlow;

__shared__ double tile[BLOCKDIMB][BLOCKDIMA];
__shared__ double flux[BLOCKDIMB][BLOCKDIMA];

int myx = blockIdx.x*BLOCKDIMB + threadIdx.x;
int myz = blockIdx.y*BLOCKDIMAM2 + threadIdx.y - 1;

if((myx >= dims.x) || (myz > dims.z)) return; // we keep an extra Y thread for the finite diff.

bool IWrite = (threadIdx.y > 0) && (threadIdx.y <= BLOCKDIMAM2) && (myz < dims.y) && (myz >= 0);
// Exclude threads at the boundary of the fluxing direction from writing back

if(myz < 0) myz += dims.z; // wrap left edge back to right edge
myz = myz % dims.z; // wrap right edge back to left

int x = myx + dims.x*dims.y*myz;
int y;

for(y = 0; y < dims.y; y++) {
    v = velGrid[x];
    b = mag[x];

    // first calculate velocityFlow
    tile[threadIdx.x][threadIdx.y] = v;
    flux[threadIdx.x][threadIdx.y] = b*v;
    __syncthreads();

    locVelFlow = (tile[threadIdx.x][threadIdx.y] + tile[threadIdx.x][(threadIdx.y+1) % BLOCKDIMA]);
    if(locVelFlow < 0.0) { locVelFlow = 1.0; } else { locVelFlow = 0.0; }

    __syncthreads();

    // Second step - calculate flux
    if(locVelFlow == 1) { tile[threadIdx.x][threadIdx.y] = flux[threadIdx.x][(threadIdx.y + 1)%BLOCKDIMA]; } else 
                        { tile[threadIdx.x][threadIdx.y] = flux[threadIdx.x][threadIdx.y]; }
   
    __syncthreads();

    // Third step - Perform flux and write to output array
    if( IWrite ) {
            bW[x] = b - lambda * ( tile[threadIdx.x][threadIdx.y] - tile[threadIdx.x][threadIdx.y-1]);
            velFlow[x] = locVelFlow;
        }

    x += dims.x;
    __syncthreads(); 
    }

}

