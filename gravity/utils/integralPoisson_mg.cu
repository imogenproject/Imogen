// Basic includes
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

// Something at some point said include this
#ifdef UNIX
#include <stdint.h>
#endif

// Matlab includes - necessary for compilation as a mex & for matrix access
#include "mex.h"
#include "matrix.h"

// CUDA headers
#include "cuda.h"
#include "cuda_runtime.h"

// GPUmat headers
#include "GPUmat.hh"

// static paramaters from GPUmat
// Pointer to the GPUmat structure & have-been-initialized variable
static GPUmat *gm;
static int init = 0;

// Access to Imogen GPU kernels
#include "cudaKernels.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
if (nrhs != 6) mexErrMsgTxt("Wrong number of arguments");

// Initialize function if not initialized
if (init == 0) {
    gm = gmGetGPUmat();
    init = 1;
    }
// Righthand side arguments:
// 0. rhos   - quantized density; CELL ARRAY
// 1. poss   - quantized positions; CELL ARRAY
// 2. phi    - Array to write potential into; double[nu nv]
// 3. bvec   - [x0 delta1x delta1y delta1z delta2x delta2y delta2z nu nv], double type
// 4. cconst - Coarsening constant (one CPU double)
// 5. h      - the grid spacing (one CPU double)

// Pull the density/position cell arrays apart and create double **s
/* Get number of argumenta and allocate memory */
int numLevels = mxGetNumberOfElements(prhs[0]);
short int *arraydims = (short int *)malloc(3*numLevels * sizeof(short int));

// Allocate mmeory for pointers to mass and position mipmaps
double **rhoptr = (double **)malloc(numLevels * sizeof(double*));
double **posptr = (double **)malloc(numLevels * sizeof(double*));

int j;
// Read the GPU variables, get the raw memory pointers,
// fill out the arrays above
for(j = 0; j < numLevels; j++) {
	mxArray *rhos = mxGetCell(prhs[0], j);
	mxArray *poss = mxGetCell(prhs[1], j);

	GPUtype rhoin = gm->gputype.getGPUtype(rhos);
	GPUtype posin = gm->gputype.getGPUtype(poss);

	int d  = gm->gputype.getNdims(rhoin);
	int *s = (int *)gm->gputype.getSize(rhoin);

	if(d == 2) {
		/* 2-D, implictly set Nz = 1 */
		arraydims[3*j+0] = s[0];
		arraydims[3*j+1] = s[1];
		arraydims[3*j+2] = 1;
		} else {
		/* 3-D */
		arraydims[3*j+0] = s[0];
		arraydims[3*j+1] = s[1];
		arraydims[3*j+2] = s[2];
		}

	rhoptr[j] = (double *)gm->gputype.getGPUptr(rhoin); 
	posptr[j] = (double *)gm->gputype.getGPUptr(posin);
	}

// Having created all these arrays above, push them to the GPU
double **GPU_rhoptr;
  cudaMalloc(&GPU_rhoptr, sizeof(double *) * numLevels);
  cudaMemcpy(GPU_rhoptr, rhoptr, sizeof(double *) * numLevels, cudaMemcpyHostToDevice);
double **GPU_posptr;
  cudaMalloc(&GPU_posptr, sizeof(double *) * numLevels);
  cudaMemcpy(GPU_posptr, posptr, sizeof(double *) * numLevels, cudaMemcpyHostToDevice);
short int *GPU_arraydims;
  cudaMalloc(&GPU_arraydims, sizeof(short int) * 3 * numLevels);
  cudaMemcpy(GPU_arraydims, arraydims, sizeof(short int) * 3 * numLevels, cudaMemcpyHostToDevice);
double *GPU_bvec;
  cudaMalloc(&GPU_bvec, 11*sizeof(double));
  cudaMemcpy(GPU_bvec, mxGetPr(prhs[3]), 11*sizeof(double), cudaMemcpyHostToDevice);

// Pass the following parameters:
//   mass**, position**, quantized array bounds
//   Array of X/Y/Z positions & # of points
//   Quantization constant
GPUtype phiout         = gm->gputype.getGPUtype(prhs[2]);
double *GPU_phiret   = (double *)gm->gputype.getGPUptr(phiout);

double coarseningConst = *mxGetPr(prhs[4]);
double *bvec = mxGetPr(prhs[3]);

// FIXME: If patches of 16 nearby cells are arranged serially, warp divergence will be reduced and speedup will be significant.

// Determine number of threads/block & blocks/grid
// Values are set in cudaKernels.h at compile-time
dim3 blockSize;
dim3 gridSize;
blockSize.x = EDGEDIM_MGBC;
blockSize.y = EDGEDIM_MGBC;
blockSize.z = 1;
gridSize.x = (int)bvec[9]  / blockSize.x; if(gridSize.x * blockSize.x < (int)bvec[9] ) gridSize.x++;
gridSize.y = (int)bvec[10] / blockSize.y; if(gridSize.y * blockSize.y < (int)bvec[10]) gridSize.y++;
gridSize.z = 1;

// Execute kernel, then free all the stuff we allocated
// Do not free any input variables, Matlab/GPUmat handle those
mgbc_kernelplane<<<gridSize, blockSize>>>(GPU_rhoptr, GPU_posptr, GPU_arraydims, numLevels, GPU_phiret, GPU_bvec, coarseningConst*pow(2,numLevels), *mxGetPr(prhs[5])/3.836451287);// <- that number is wrong

cudaFree(GPU_bvec);
cudaFree(GPU_rhoptr);
cudaFree(GPU_posptr);
cudaFree(GPU_arraydims);
free(rhoptr);
free(posptr);
free(arraydims);

return;
}


