// Include standard headers
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#ifdef UNIX
#include <stdint.h>
#include <unistd.h>
#endif

// Standard Matlab MEX file includes
#include "mex.h"
#include "matrix.h"

// CUDA header includes
#include "cuda.h"
#include "cuda_runtime.h"
#include "cublas.h"

// GPUmat include
#include "GPUmat.hh"

// Headers for Imogen kernels
#include "cudaKernels.h"

// The B operator is used to calculate power series approximate inverses of the Laplace matrix operator
// as the meat of the polynomial preconditioner.
// It consists of the off-diagonal matrix elements of the operator with the diagonal normalized to 1
__global__ void Laplacian_B_OperatorKernel(double *src, double *dst, double accumCoefficient, double *accum, int nx, int ny, int nz, int order)
{

//if (blockIdx.x > 0) return;
//if (blockIdx.y > 1) return;

// Allocate local storage for 3 planes of data
// Since each plane depends only on itself and the planes above/below it, the strategy is to
// store only those 3 planes seen by the current plane and roll along, minimizing waste on loading edges.
__shared__ double locSrc[EDGEDIM_BOP+2][EDGEDIM_BOP+2][3];

// Decides if this thread computes B or is just a data loader for the edge values
bool IAmInterior = false;
if ((threadIdx.x > 0) && (threadIdx.x <= EDGEDIM_BOP) && (threadIdx.y > 0) && (threadIdx.y <= EDGEDIM_BOP)) { IAmInterior = true; }

// Calculate some global indices and offsets
int myx = threadIdx.x + EDGEDIM_BOP*blockIdx.x - 1;
int myy = threadIdx.y + EDGEDIM_BOP*blockIdx.y - 1;
int currentZ = 0;

int myAddr = myx + nx*myy;
// Determine if this thread would index outside the array bounds.
bool IAmInRange  = (myx >= 0) && (myy >= 0) && (myx < nx) && (myy < ny);

// Load up the first planes of data, deliberately one plane too high
locSrc[threadIdx.x][threadIdx.y][1] = 0.0;

// Setup the indirection that lets us quickly flip data planes
if(IAmInRange) {
        locSrc[threadIdx.x][threadIdx.y][2] = src[myAddr];
        } else { locSrc[threadIdx.x][threadIdx.y][2] = 0.0; }
int plidx[3];
plidx[0]=0;
plidx[1]=1;
plidx[2]=2;

__syncthreads();

while(currentZ < nz) {
        // Roll all 3 planes downwards
        plidx[0] = ++plidx[0] % 3;
        plidx[1] = ++plidx[1] % 3;
        plidx[2] = ++plidx[2] % 3;

        // Reload the top one
        if((currentZ < (nz-1)) && IAmInRange)
                locSrc[threadIdx.x][threadIdx.y][plidx[2]] = src[myAddr + nx*ny];
                else locSrc[threadIdx.x][threadIdx.y][plidx[2]] = 0.0;

        __syncthreads();

        // Compute B operator
        if(IAmInterior && IAmInRange) {
                double res;
                if(order == 2) {
                res   = ( locSrc[threadIdx.x+1][threadIdx.y  ][plidx[1]]
                        + locSrc[threadIdx.x-1][threadIdx.y  ][plidx[1]]
                        + locSrc[threadIdx.x  ][threadIdx.y+1][plidx[1]]
                        + locSrc[threadIdx.x  ][threadIdx.y-1][plidx[1]]
                        + locSrc[threadIdx.x  ][threadIdx.y  ][plidx[2]]
                        + locSrc[threadIdx.x  ][threadIdx.y  ][plidx[0]] ) / 6.0;
                } else {
                res   = ( locSrc[threadIdx.x-1][threadIdx.y  ][plidx[0]] // lower plane
                        + locSrc[threadIdx.x  ][threadIdx.y-1][plidx[0]]
                        + locSrc[threadIdx.x  ][threadIdx.y+1][plidx[0]]
                        + locSrc[threadIdx.x+1][threadIdx.y  ][plidx[0]]
                        + 2*locSrc[threadIdx.x  ][threadIdx.y][plidx[0]]

                        + locSrc[threadIdx.x-1][threadIdx.y-1][plidx[1]] // Center plane
                        + locSrc[threadIdx.x+1][threadIdx.y-1][plidx[1]]
                        + locSrc[threadIdx.x-1][threadIdx.y+1][plidx[1]]
                        + locSrc[threadIdx.x+1][threadIdx.y+1][plidx[1]]
                        +2*locSrc[threadIdx.x+1][threadIdx.y ][plidx[1]]
                        +2*locSrc[threadIdx.x-1][threadIdx.y ][plidx[1]]
                        +2*locSrc[threadIdx.x ][threadIdx.y+1][plidx[1]]
                        +2*locSrc[threadIdx.x ][threadIdx.y-1][plidx[1]]

                        + locSrc[threadIdx.x-1][threadIdx.y  ][plidx[2]] // upper plane
                        + locSrc[threadIdx.x  ][threadIdx.y-1][plidx[2]]
                        + locSrc[threadIdx.x  ][threadIdx.y+1][plidx[2]]
                        + locSrc[threadIdx.x+1][threadIdx.y  ][plidx[2]]
                        + 2*locSrc[threadIdx.x][threadIdx.y  ][plidx[2]] ) / 24.0;
                }
                

                dst[myAddr] = res;
//                dst[myAddr] = locSrc[threadIdx.x][threadIdx.y][plidx[1]];
//                accum[myAddr] = locSrc[threadIdx.x][threadIdx.y][plidx[1]];
                accum[myAddr] += accumCoefficient * res;
                }

	__syncthreads();

        // Increment z and address pointer
        currentZ++;
        myAddr += nx*ny;
        if(currentZ >= nz) break;
        }


}

// This is the core Multigrid algorithm for evaluating the integral poisson equation in logarithmic time
// WARNING: I cannot figure out why but THIS DOES NOT WORK ON THE FERMI ARCHITECTURE
#define LP_STATE_SIZE 4
#define THREADS_PER_BLOCK EDGEDIM_MGBC * EDGEDIM_MGBC
// Multilevel BC finding kernel
__global__ void mgbc_kernelplane(double **mass,
                          double **position,
                          short int *arrayDims,
                          int nLevels,
                          double *phi,
                          double *bvec,
                          double qConst, double selfPotRad)
{  

// Keep a local copy of the potentials this block is finding to reduce global IO
// Allocate enough loop stack space for all threads on shared mem
__shared__ double locPhiStore[EDGEDIM_MGBC][EDGEDIM_MGBC];
__shared__ short int loopstack[THREADS_PER_BLOCK * LP_STATE_SIZE * 10];

// If we are beyond the array to be written, return.
if(blockIdx.x*EDGEDIM_MGBC+threadIdx.x >= (int)bvec[9] ) return;
if(blockIdx.y*EDGEDIM_MGBC+threadIdx.y >= (int)bvec[10]) return;

locPhiStore[threadIdx.x][threadIdx.y] = 0.0;


double myX = (double)(blockIdx.x*EDGEDIM_MGBC+threadIdx.x)*bvec[3] + (double)(blockIdx.y*EDGEDIM_MGBC+threadIdx.y)*bvec[6] + bvec[0];
double myY = (double)(blockIdx.x*EDGEDIM_MGBC+threadIdx.x)*bvec[4] + (double)(blockIdx.y*EDGEDIM_MGBC+threadIdx.y)*bvec[7] + bvec[1];
double myZ = (double)(blockIdx.x*EDGEDIM_MGBC+threadIdx.x)*bvec[5] + (double)(blockIdx.y*EDGEDIM_MGBC+threadIdx.y)*bvec[8] + bvec[2];

short int CL = 0; // current level
double rad;
double deltaphi;

// Loop state. If you can't make function calls, you can still make a stack and use goto!
short int *lpst = &loopstack[(threadIdx.x + EDGEDIM_MGBC*threadIdx.y) * LP_STATE_SIZE * 10];
short int blockIter[3];

//Loop over coarsest level
// loopstack values: [index x y z, who_to_return_to]
for(blockIter[0] = 0; blockIter[0] < arrayDims[0]; blockIter[0]++) {
for(blockIter[1] = 0; blockIter[1] < arrayDims[1]; blockIter[1]++) {
for(blockIter[2] = 0; blockIter[2] < arrayDims[2]; blockIter[2]++) {
        lpst[0] = blockIter[0];
        lpst[1] = blockIter[1];
        lpst[2] = blockIter[2];
        lpst[3] = 0;
        
        goto MasterBlockIter;
        MasterBlockReturn:
        } } }

// Make one pushback to global memory
phi[(blockIdx.x*EDGEDIM_MGBC+threadIdx.x) + (blockIdx.y*EDGEDIM_MGBC+threadIdx.y)*(int)bvec[9] ] = locPhiStore[threadIdx.x][threadIdx.y];

return; // Function return

MasterBlockIter:
    
double *ptpos = &position[CL][3*(lpst[0] + arrayDims[3*CL+0]*(lpst[1] + arrayDims[3*CL+1]*lpst[2]))];    
rad = sqrt( (ptpos[0] - myX)*(ptpos[0] - myX) +
            (ptpos[1] - myY)*(ptpos[1] - myY) +
            (ptpos[2] - myZ)*(ptpos[2] - myZ));

if (rad < 1e-5) rad = selfPotRad;

deltaphi = mass[CL][lpst[0] + arrayDims[3*CL+0]*(lpst[1] + arrayDims[3*CL+1]*lpst[2])] / rad;
//deltaphi = mass[CL][lpst[0] + arrayDims[3*CL+0]*(lpst[1] + arrayDims[3*CL+1]*lpst[2])]*(1 + ptpos[0]*ptpos[0]/(rad*rad))/(rad*rad*rad);

//        if we are far enough away or can't go deeper, sum up this reduced block
//if( (rad > qConst) || (CL >= (nLevels-1)) ) {
if( (rad > qConst) || (CL >= 0) ) {
        locPhiStore[threadIdx.x][threadIdx.y] -= deltaphi;
        } else {
        CL++;
        qConst /= 2.0;
        // Otherwise, go up a level and recurse over the 8 subcells        
        for(lpst[4] = 2*lpst[0]; lpst[4] < 2*lpst[0]+2; lpst[4]++)
        for(lpst[5] = 2*lpst[1]; lpst[5] < 2*lpst[1]+2; lpst[5]++)
        for(lpst[6] = 2*lpst[2]; lpst[6] < 2*lpst[2]+2; lpst[6]++) {
                lpst[7] = 1;
//                CL++;
                lpst = &lpst[LP_STATE_SIZE]; // move up stack
                goto MasterBlockIter;

                iterReturn:
                lpst = &lpst[-LP_STATE_SIZE];
//                CL--;
                }
        CL--;
        qConst *= 2.0;
        }

// We're finished with this block, current level will be 0. Otherwise, jump back.
//if(CL > 0) { goto iterReturn; } else { goto MasterBlockReturn; }
if(lpst[3] == 1) { goto iterReturn; } else { goto MasterBlockReturn; }

}


// Every thread is responsible for calculating B at one of 512 locations in the block.
__global__ void SymmetricOperatorKernel(double *src, double *dst, int nx, int ny, int nz, double c1, double c2, double c3, double c4)
{

// 8K of shared memory, gone in an instant...
__shared__ double localSrc[10][10][10];

int trueBlkZ = blockIdx.x / (nx/8);
int trueBlkX = blockIdx.x - trueBlkZ*(nx/8);

// We are initially invoked with 10x10 threads in the XY plane
int myX = threadIdx.x + 8*trueBlkX - 1;
int myY = threadIdx.y + 8*blockIdx.y - 1;
int myZ = 8*trueBlkZ;

int myAddrBase = myX + nx*(myY + ny*myZ);
int counter;

if(trueBlkZ > 0) {
        trueBlkX = 0;
        myAddrBase -= nx*ny;
} else {
        trueBlkX = 1;
        localSrc[threadIdx.x][threadIdx.y][0] = 0.0;
}

if(trueBlkZ < nz) {
        trueBlkZ = 10;
} else {
        trueBlkZ = 9;
        localSrc[threadIdx.x][threadIdx.y][9] = 0.0;

}

//for(counter = 0; counter < 10; counter++) { localSrc[threadIdx.x][threadIdx.y][counter] = 0.0; }

// Load
for(counter = trueBlkX; counter < trueBlkZ; counter++) {

        if((myX >= 0) && (myX < nx) && (myY >= 0) && (myY < ny))
                localSrc[threadIdx.x][threadIdx.y][counter] = src[myAddrBase];
        else
                localSrc[threadIdx.x][threadIdx.y][counter] = 0.0;

        myAddrBase += nx*ny;
        }

if((threadIdx.x == 0) || (threadIdx.x == 9) || (threadIdx.y == 0) || (threadIdx.y == 9)) return;

// Make sure the local src block is ready.
__syncthreads();

myAddrBase = myX + nx*(myY + ny*myZ);
double res;

for(counter = 1; counter < 9; counter++) {
        res = c1 * localSrc[threadIdx.x][threadIdx.y][counter];
        if(c2 != 0.0) res += c2*( localSrc[threadIdx.x-1][threadIdx.y][counter]
                                + localSrc[threadIdx.x+1][threadIdx.y][counter]
                                + localSrc[threadIdx.x][threadIdx.y-1][counter]
                                + localSrc[threadIdx.x][threadIdx.y+1][counter]
                                + localSrc[threadIdx.x][threadIdx.y][counter-1]
                                + localSrc[threadIdx.x][threadIdx.y][counter+1] );

        if(c3 != 0.0) res += c3*( localSrc[threadIdx.x-1][threadIdx.y][counter-1]
                                + localSrc[threadIdx.x+1][threadIdx.y][counter-1]
                                + localSrc[threadIdx.x][threadIdx.y+1][counter-1]
                                + localSrc[threadIdx.x][threadIdx.y-1][counter-1]
                                + localSrc[threadIdx.x-1][threadIdx.y-1][counter]
                                + localSrc[threadIdx.x+1][threadIdx.y-1][counter]
                                + localSrc[threadIdx.x-1][threadIdx.y+1][counter]
                                + localSrc[threadIdx.x+1][threadIdx.y+1][counter]
                                + localSrc[threadIdx.x-1][threadIdx.y][counter+1]
                                + localSrc[threadIdx.x+1][threadIdx.y][counter+1]
                                + localSrc[threadIdx.x][threadIdx.y-1][counter+1]
                                + localSrc[threadIdx.x][threadIdx.y+1][counter+1] );

        if(c4 != 0.0) res += c4*( localSrc[threadIdx.x-1][threadIdx.y-1][counter-1]
                                + localSrc[threadIdx.x+1][threadIdx.y-1][counter-1]
                                + localSrc[threadIdx.x-1][threadIdx.y+1][counter-1]
                                + localSrc[threadIdx.x+1][threadIdx.y+1][counter-1]
                                + localSrc[threadIdx.x-1][threadIdx.y-1][counter+1]
                                + localSrc[threadIdx.x+1][threadIdx.y-1][counter+1]
                                + localSrc[threadIdx.x-1][threadIdx.y+1][counter+1]
                                + localSrc[threadIdx.x+1][threadIdx.y+1][counter+1] );

        dst[myAddrBase] = res;
        myAddrBase += nx*ny;
        }


}

// nx/ny/nz are size of source array, for which each cell has a thread
__global__ void upsampleKernel(double *src, double *dst, int factor, int nx, int ny, int nz)
{

int myx = threadIdx.x + blockDim.x * blockIdx.x;
if (myx >= nx) return;

int myy = threadIdx.y + blockDim.y * blockIdx.y;
if (myy >= ny) return;

int z;
double myval;

int lx = nx * factor;
int ly = ny * factor;

int chix, chiy, chiz;

for (z = 0; z < nz; z++) {
        myval = src[myx + nx*myy + nx*ny*z];

        for (chix = 0; chix < factor; chix++)
        for (chiy = 0; chiy < factor; chiy++)
        for (chiz = 0; chiz < factor; chiz++) {
//src[myx + nx*myy + nx*ny*z] = (double)(factor*myx + chix + lx*(factor*myy + chiy   + ly*(factor*z + chiz) ));
                dst[factor*myx + chix + lx*(factor*myy + chiy   + ly*(factor*z + chiz) )] = myval;
                }
        }

}

// nx/ny/nz are size of destination array for which each cell has a thread.
__global__ void downsampleKernel(double *src, double *dst, int factor, int nx, int ny, int nz)
{

int myx = threadIdx.x + blockDim.x * blockIdx.x;
if (myx >= nx) return;

int myy = threadIdx.y + blockDim.y * blockIdx.y;
if (myy >= ny) return;

int z;
double myval;

int lx = nx * factor;
int ly = ny * factor;

int chix, chiy, chiz;

for (z = 0; z < nz; z++) {
        myval = 0.0;

        for (chix = 0; chix < factor; chix++)
        for (chiy = 0; chiy < factor; chiy++)
        for (chiz = 0; chiz < factor; chiz++) {
                myval += src[factor*myx + chix + lx*(factor*myy + chiy   + ly*(factor*z + chiz) )];
                }
        dst[myx + nx*(myy + ny*z)] = myval;// / (double)(factor*factor*factor);
        }

}

// nx/ny/nz are size of destination array for which each cell has a thread.
__global__ void downsampleKernel_lowdim(double *src, double *dst, int factor, int nx, int ny, int nz)
{

int myx = threadIdx.x + blockDim.x * blockIdx.x;
if (myx >= nx) return;

int myy = threadIdx.y + blockDim.y * blockIdx.y;
if (myy >= ny) return;

int z;
double myval;

int lx = nx * factor;
int ly = ny * factor;

int chix, chiy, chiz;

for (z = 0; z < nz; z++) {
        myval = 0.0;

        for (chix = 0; chix < factor; chix++)
        for (chiy = 0; chiy < factor; chiy++)
        for (chiz = 0; chiz < factor; chiz++) {
                myval += src[factor*myx + chix + lx*(factor*myy + chiy   + ly*(factor*z + chiz) )];
                }
        dst[myx + nx*(myy + ny*z)] = myval;// / (double)(factor*factor*factor);
        }

}

