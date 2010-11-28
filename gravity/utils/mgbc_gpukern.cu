//=================================================================================================
//                                                                               I N C L U D E S 

// Basic includes
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

// Needed by unix because ???
#ifdef UNIX
#include <stdint.h>
#endif

// Matlab includes. Required to compile as a Mex function.
#include "mex.h"
#include "matrix.h"

// CUDA includes
#include "cuda.h"
#include "cuda_runtime.h"

// GPUmat includes
#include "GPUmat.hh"

//=================================================================================================
//                                                                         D E F I N I T I O N S

#define THREADS_PER_BLOCK 16
#define LP_STATE_SIZE 4

//=================================================================================================
//                                                               S T A T I C   V A R I A B L E S

/** Specifies whether or not the GPU has been initialized for use already. */
static int GPUInitialized = 0;

/** Reference to the GPUmat instance. */
static GPUmat *gpumat;

//=================================================================================================
//                                                                                 H E A D E R S 

// WTF is mgbc_kern suppose to mean. I get it after a while but what's wrong with naming it in
// English so other developers have a chance to figure this code out?
__global__ void mgbc_kern(double **mass, double **position, short int *arrayDims, int numLevels, 
                          double *points, double *phi,  int nPoints, double qConst, 
                          double selfPotRad);

//=================================================================================================
//                                                                             F U N C T I O N S

//_________________________________________________________________________________________________ mexFunction
/** Main entry point for the compiled mex function.
 * The prhs mxArray consists of the following items:
 * 0. rhos   - quantized mass density
 * 1. poss   - quantized positions
 * 2. points - the set of points to find phi at
 * 3. phi    - Blank array to write potential into
 * 4. cconst - Coarsening constant (lone CPU double)
 * 5. h      - the grid spacing */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // What's this for?
    //CUresult cudastatus = CUDA_SUCCESS;

    if (nrhs != 6) 
        mexErrMsgTxt("Wrong number of arguments");

    // Initialize static GPU instance if not already initialized by a previous call.
    if (GPUInitialized == 0) {
        gpumat           = gmGetGPUmat(); // Access GPUmat instance.
        GPUInitialized   = 1;
    }

    // Don't get what's going on here. Why is numLevels the same as numElements?
    int numLevels        = mxGetNumberOfElements(prhs[0]);
    short int *arraydims = (short int *)malloc(3*numLevels * sizeof(short int));
    
    // Allocate memory for the mass and position arrays
    double **rhoptr = (double **)malloc(numLevels * sizeof(double*));
    double **posptr = (double **)malloc(numLevels * sizeof(double*));
    
    // Build arrays of the size of every level and pointers to them 
    int j;
    for(j = 0; j < numLevels; j++) {
        mxArray *rhos = mxGetCell(prhs[0], j);
        mxArray *poss = mxGetCell(prhs[1], j);
    
        GPUtype rhoin = gpumat->gputype.getGPUtype(rhos);
        GPUtype posin = gpumat->gputype.getGPUtype(poss);
    
        int d  = gpumat->gputype.getNdims(rhoin);
        int *s = (int *)gpumat->gputype.getSize(rhoin);
    
        if(d == 2) {
            /* 2-D */
            arraydims[3*j + 0] = s[0];
            arraydims[3*j + 1] = s[1];
            arraydims[3*j + 2] = 1;
            } else {
            /* 3-D, it better not be higher-dimensional */
            arraydims[3*j+0] = s[0];
            arraydims[3*j+1] = s[1];
            arraydims[3*j+2] = s[2];
            }
    
        rhoptr[j] = (double *)gpumat->gputype.getGPUptr(rhoin); 
        posptr[j] = (double *)gpumat->gputype.getGPUptr(posin);
    }
    
    // Get rhoptr, posptr and arraydims to the GPU
    double **GPU_rhoptr;
      cudaMalloc(&GPU_rhoptr, sizeof(double *)*numLevels);
      cudaMemcpy(GPU_rhoptr, rhoptr, sizeof(double *)*numLevels, cudaMemcpyHostToDevice);
      
    double **GPU_posptr;
      cudaMalloc(&GPU_posptr, sizeof(double *)*numLevels);
      cudaMemcpy(GPU_posptr, posptr, sizeof(double *)*numLevels, cudaMemcpyHostToDevice);
      
    short int *GPU_arraydims;
      cudaMalloc(&GPU_arraydims, sizeof(short int)*3*numLevels);
      cudaMemcpy(GPU_arraydims, arraydims, sizeof(short int)*3*numLevels, 
                 cudaMemcpyHostToDevice);
    
    // Pass the following parameters:
    //   mass**, position**, quantized array bounds
    //   Array of X/Y/Z positions & # of points
    //   Quantization constant
    
    GPUtype pointsin       = gpumat->gputype.getGPUtype(prhs[2]);
      double *GPU_pointSet = (double *)gpumat->gputype.getGPUptr(pointsin);
      int numPoints        = gpumat->gputype.getNumel(pointsin) / 3;
    
    GPUtype phiout         = gpumat->gputype.getGPUtype(prhs[3]);
      double *GPU_phiret   = (double *)gpumat->gputype.getGPUptr(phiout);
    
    double coarseningConst = *mxGetPr(prhs[4]);
    
    // FIXME: If patches of 16 nearby cells are arranged serially, warp divergence will be reduced and speedup will be significant.
    
    // Determine number of threads/block & number of blocks
    // Assume the point to be arranced in a grid pattern
    int blockSize = THREADS_PER_BLOCK;
    int gridSize  = numPoints/blockSize;
    if ((gridSize*blockSize) < numPoints) 
        gridSize++;
    
    // Run the GPU kernel.
    mgbc_kern<<<gridSize, blockSize>>>(GPU_rhoptr, GPU_posptr, GPU_arraydims, numLevels,
                                       GPU_pointSet, GPU_phiret, numPoints, coarseningConst*pow(2,numLevels),
                                       *mxGetPr(prhs[5])/3.836451287);
        
    // Cleanup both GPU and CPU memory
    cudaFree(GPU_rhoptr);
    cudaFree(GPU_posptr);
    cudaFree(GPU_arraydims);
    free(rhoptr);
    free(posptr);
    free(arraydims);
    
    return;
}

//_________________________________________________________________________________________________ mgbc_kern
/** Multi-level boundary condition finding GPU kernel. */
__global__ void mgbc_kern(double **mass, double **position, short int *arrayDims, int numLevels,
                          double *points, double *phi, int nPoints, double qConst, 
                          double selfPotRad) {  

    // Create variables to be used during recursion.
    double rad;
    double deltaphi;
                          
    // Keep a local copy of the potentials this block is finding to reduce global IO
    // Allocate enough loop stack space for all threads on shared mem
    __shared__ double locPhiStore[THREADS_PER_BLOCK];
    __shared__ short int loopstack[THREADS_PER_BLOCK * LP_STATE_SIZE * 10];

    // Identify my thread/block/grid number and determine what point I'm finding phi at
    int myPointIdx = threadIdx.x + blockIdx.x*blockDim.x;
    if (myPointIdx > nPoints) return;

    locPhiStore[threadIdx.x] = 0.0;
    
    double myX = points[3*myPointIdx+0];
    double myY = points[3*myPointIdx+1];
    double myZ = points[3*myPointIdx+2];

    short int CL = 0; // current level

    // Loop state. If you can't make function calls, you can still make a stack and use goto!
    short int *loopState = &loopstack[threadIdx.x * LP_STATE_SIZE * 10];
    short int blockIter[3];

    //Loop over coarsest level
    // loopstack values: [index x y z, who_to_return_to]
    for(blockIter[0] = 0; blockIter[0] < arrayDims[0]; blockIter[0]++) {
        for(blockIter[1] = 0; blockIter[1] < arrayDims[1]; blockIter[1]++) {
            for(blockIter[2] = 0; blockIter[2] < arrayDims[2]; blockIter[2]++) {
                loopState[0] = blockIter[0];
                loopState[1] = blockIter[1];
                loopState[2] = blockIter[2];
                loopState[3] = 0;
                
                goto MasterBlockIter;
                MasterBlockReturn:
            } 
        }
    }

    // Make one pushback to global memory
    phi[myPointIdx] = locPhiStore[threadIdx.x];

    return;

    //---------------------------------------------------------------------------------------------
    // Pseudo-recursion portion of the algorithm exists below the function return and is accessible
    // by goto statements only. Unfortunately, this is necessary given the recursion limitations of
    // the GPU.
    //---------------------------------------------------------------------------------------------
    
    MasterBlockIter:
    
    double *ptpos = &position[CL][3*(loopState[0] + arrayDims[3*CL+0]*(loopState[1] 
                    + arrayDims[3*CL+1]*loopState[2]))];
                    
    rad           = sqrt( (ptpos[0] - myX)*(ptpos[0] - myX) +
                          (ptpos[1] - myY)*(ptpos[1] - myY) +
                          (ptpos[2] - myZ)*(ptpos[2] - myZ));

    if (rad < 1e-5) 
        rad = selfPotRad;
    
    deltaphi = mass[CL][loopState[0] + arrayDims[3*CL+0]*(loopState[1] + arrayDims[3*CL+1]*loopState[2])]/rad;
    
    //	if we are far enough away or can't go deeper, sum up this reduced block
    if( (rad > qConst) || (CL >= (numLevels-1)) ) //DEPRECATED METHOD: if( (rad  > qConst * pow(.5,CL)) || (CL >= (numLevels-1)) ) {
        locPhiStore[threadIdx.x] -= deltaphi;
        
    // Otherwise, go up a level and recurse over the 8 subcells	
    else {
        // Adjust to higher level for following recursion call.
        CL++;
        qConst /= 2.0;
        
        for(loopState[4] = 2*loopState[0]; loopState[4] < 2*loopState[0]+2; loopState[4]++) {
            for(loopState[5] = 2*loopState[1]; loopState[5] < 2*loopState[1]+2; loopState[5]++) {
                for(loopState[6] = 2*loopState[2]; loopState[6] < 2*loopState[2]+2; loopState[6]++) {
                    loopState[7] = 1;
                    loopState = &loopState[LP_STATE_SIZE]; // move up stack
                    goto MasterBlockIter;
    
                    iterReturn:
                    loopState = &loopState[-LP_STATE_SIZE];
                }
            }
        }
        
        // Return to current recursion state.
        CL--;
        qConst *= 2.0;
    }

    // If we're finished with this block, the current level will be 0. Otherwise, jump back for 
    // another level of recursive summing.
    if(loopState[3] == 1)
        goto iterReturn; 
    else 
        goto MasterBlockReturn;
}

