// Basic includes
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#ifdef UNIX
#include <stdint.h>
#endif

// Matlab
#include "mex.h"
#include "matrix.h"

// CUDA
#include "cuda.h"
#include "cuda_runtime.h"

// GPUmat
#include "GPUmat.hh"

// static paramaters
static GPUmat *gm;
static int init = 0;
#define THREADS_PER_BLOCK 16

// Header
__global__ void mgbc_kern(double **mass, double **position, short int *arrayDims, int nLevels, double *points, double *phi,  int nPoints, double qConst, double selfPotRad);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//CUresult cudastatus = CUDA_SUCCESS;

if (nrhs != 6) mexErrMsgTxt("Wrong number of arguments");

if (init == 0) {
// Initialize function
// load GPUmat
gm = gmGetGPUmat();

init = 1;
}
// Righthand side arguments:
// 0. rhos   - quantized density
// 1. poss   - quantized positions
// 2. points - the set of points to find phi at
// 3. phi    - Blank array to write potential into
// 4. cconst - Coarsening constant (lone CPU double)
// 5. h      - the grid spacing

// Pull the density/position cell arrays apart and create double **s
/* Get number of argument and allocate memory */
int nLevels = mxGetNumberOfElements(prhs[0]);

short int *arraydims = (short int *)malloc(3*nLevels * sizeof(short int));
double **rhoptr = (double **)malloc(nLevels * sizeof(double*));
double **posptr = (double **)malloc(nLevels * sizeof(double*));

int j;
/* Build arrays of the size of every level and *s to them */
// Pull the density/position cell arrays apart and create double **s
for(j = 0; j < nLevels; j++) {
	mxArray *rhos = mxGetCell(prhs[0], j);
	mxArray *poss = mxGetCell(prhs[1], j);

	GPUtype rhoin = gm->gputype.getGPUtype(rhos);
	GPUtype posin = gm->gputype.getGPUtype(poss);

	int d  = gm->gputype.getNdims(rhoin);
	int *s = (int *)gm->gputype.getSize(rhoin);

	if(d == 2) {
		/* 2-D */
		arraydims[3*j+0] = s[0];
		arraydims[3*j+1] = s[1];
		arraydims[3*j+2] = 1;
		} else {
		/* 3-D, it better not be higher-dimensional */
		arraydims[3*j+0] = s[0];
		arraydims[3*j+1] = s[1];
		arraydims[3*j+2] = s[2];
		}

	rhoptr[j] = (double *)gm->gputype.getGPUptr(rhoin); 
	posptr[j] = (double *)gm->gputype.getGPUptr(posin);
	}

// Get rhoptr, posptr and arraydims to the GPU
double **GPU_rhoptr;
  cudaMalloc(&GPU_rhoptr, sizeof(double *) * nLevels);
  cudaMemcpy(GPU_rhoptr, rhoptr, sizeof(double *) * nLevels, cudaMemcpyHostToDevice);
double **GPU_posptr;
  cudaMalloc(&GPU_posptr, sizeof(double *) * nLevels);
  cudaMemcpy(GPU_posptr, posptr, sizeof(double *) * nLevels, cudaMemcpyHostToDevice);
short int *GPU_arraydims;
  cudaMalloc(&GPU_arraydims, sizeof(short int) * 3 * nLevels);
  cudaMemcpy(GPU_arraydims, arraydims, sizeof(short int) * 3 * nLevels, cudaMemcpyHostToDevice);

// Pass the following parameters:
//   mass**, position**, quantized array bounds
//   Array of X/Y/Z positions & # of points
//   Quantization constant

GPUtype pointsin       = gm->gputype.getGPUtype(prhs[2]);
  double *GPU_pointSet = (double *)gm->gputype.getGPUptr(pointsin);
  int numPoints        = gm->gputype.getNumel(pointsin) / 3;

GPUtype phiout         = gm->gputype.getGPUtype(prhs[3]);
  double *GPU_phiret   = (double *)gm->gputype.getGPUptr(phiout);

double coarseningConst = *mxGetPr(prhs[4]);

// FIXME: If patches of 16 nearby cells are arranged serially, warp divergence will be reduced and speedup will be significant.

// Determine number of threads/block & number of blocks
// Assume the point to be arranced in a grid pattern
int blockSize = 16; //THREADS_PER_BLOCK;
int gridSize = numPoints/blockSize;
if ((gridSize*blockSize) < numPoints) gridSize++;

mgbc_kern<<<gridSize, blockSize>>>(GPU_rhoptr, GPU_posptr, GPU_arraydims, nLevels,
                                   GPU_pointSet, GPU_phiret, numPoints, coarseningConst*pow(2,nLevels),
				   *mxGetPr(prhs[5])/3.836451287);


cudaFree(GPU_rhoptr);
cudaFree(GPU_posptr);
cudaFree(GPU_arraydims);
free(rhoptr);
free(posptr);
free(arraydims);

return;
}

#define LP_STATE_SIZE 4
// Multilevel BC finding kernel
__global__ void mgbc_kern(double **mass,
                          double **position,
                          short int *arrayDims,
                          int nLevels,
                          double *points,
                          double *phi,
                          int nPoints,
                          double qConst, double selfPotRad)
{  

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
short int *lpst = &loopstack[threadIdx.x * LP_STATE_SIZE * 10];
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
phi[myPointIdx] = locPhiStore[threadIdx.x];

return; // Function return

double rad;
double deltaphi;
MasterBlockIter:
    
double *ptpos = &position[CL][3*(lpst[0] + arrayDims[3*CL+0]*(lpst[1] + arrayDims[3*CL+1]*lpst[2]))];    
rad = sqrt( (ptpos[0] - myX)*(ptpos[0] - myX) +
            (ptpos[1] - myY)*(ptpos[1] - myY) +
            (ptpos[2] - myZ)*(ptpos[2] - myZ));

if (rad < 1e-5) rad = selfPotRad;

deltaphi = mass[CL][lpst[0] + arrayDims[3*CL+0]*(lpst[1] + arrayDims[3*CL+1]*lpst[2])] /rad;

//	if we are far enough away or can't go deeper, sum up this reduced block
//if( (rad  > qConst * pow(.5,CL)) || (CL >= (nLevels-1)) ) {
if( (rad > qConst) || (CL >= (nLevels-1)) ) {
	locPhiStore[threadIdx.x] -= deltaphi;
	} else {
	CL++;
	qConst /= 2.0;
	// Otherwise, go up a level and recurse over the 8 subcells	
	for(lpst[4] = 2*lpst[0]; lpst[4] < 2*lpst[0]+2; lpst[4]++)
	for(lpst[5] = 2*lpst[1]; lpst[5] < 2*lpst[1]+2; lpst[5]++)
	for(lpst[6] = 2*lpst[2]; lpst[6] < 2*lpst[2]+2; lpst[6]++) {
		lpst[7] = 1;
//		CL++;
		lpst = &lpst[LP_STATE_SIZE]; // move up stack
		goto MasterBlockIter;

	        iterReturn:
		lpst = &lpst[-LP_STATE_SIZE];
//		CL--;
		}
	CL--;
	qConst *= 2.0;
	}

// We're finished with this block, current level will be 0. Otherwise, jump back.
//if(CL > 0) { goto iterReturn; } else { goto MasterBlockReturn; }
if(lpst[3] == 1) { goto iterReturn; } else { goto MasterBlockReturn; }

}

