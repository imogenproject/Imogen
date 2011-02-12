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

typedef struct {
	double *fluidIn[5];
	double *fluidOut[5];

	double *B[3];

	double *Ptotal;
	double *cFreeze;
	} fluidVarPtrs;

#include "cudaCommon.h"

//double **getGPUSourcePointers(const mxArray *prhs[], int num, int *retNumel, int startat);
//double **makeGPUDestinationArrays(GPUtype src, mxArray *retArray[], int howmany);

__global__ void cudaWFluxKernel(fluidVarPtrs fluid, double lambda, int nu, int hu, int hv, int hw, int direction);
__global__ void nullStep(fluidVarPtrs fluid, int numel);

#define BLOCKLEN 126

/* Given the RHS and how many cuda arrays we expect, extracts a set of pointers to GPU memory for us
 Also conveniently checked for equal array extent and returns it for us */
/*double **getGPUSourcePointers(const mxArray *prhs[], int num, int *retNumel, int startat)
{
  GPUtype src;
  double **gpuPointers = (double **)malloc(num * sizeof(double *));
  int iter;
  int numel = gm->gputype.getNumel(gm->gputype.getGPUtype(prhs[startat]));

  for(iter = 0; iter < num; iter++) {
    src = gm->gputype.getGPUtype(prhs[startat + iter]);
    if (gm->gputype.getNumel(src) != numel) { free(gpuPointers); mexErrMsgTxt("Fatal: Arrays contain nonequal number of elements."); }
    gpuPointers[iter] = (double *)gm->gputype.getGPUptr(src);
  }

retNumel[0] = numel;
return gpuPointers;
}
*/
/* Creates destination array that the kernels write to; Returns the GPU memory pointer, and assigns the LHS it's passed */
/*
double **makeGPUDestinationArrays(GPUtype src, mxArray *retArray[], int howmany)
{
int d = gm->gputype.getNdims(src);
const int *ssize = gm->gputype.getSize(src);
int x;
int newsize[3];
for(x = 0; x < 3; x++) (x < d) ? newsize[x] = ssize[x] : newsize[x] = 1;

double **rvals = (double **)malloc(howmany*sizeof(double *));
int i;
for(i = 0; i < howmany; i++) {
  GPUtype ra = gm->gputype.create(gpuDOUBLE, d, newsize, NULL);
  retArray[i] = gm->gputype.createMxArray(ra);
  rvals[i] = (double *)gm->gputype.getGPUptr(ra);
  }

return rvals;

}*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // At least 2 arguments expected
  // Input and result
  if ((nrhs!=12) || (nlhs != 5)) mexErrMsgTxt("Wrong number of arguments: need [5] = cudaWflux(12)\n");

  if (init == 0) {
    gm = gmGetGPUmat();
    init = 1;
  }


//double **getGPUSourcePointers(const mxArray *prhs[], int num, int *retNumel, int startat);
//double **makeGPUDestinationArrays(GPUtype src, mxArray *retArray[], int howmany);

  // Get source array info and create destination arrays
  int numel;
  GPUtype srcReference = gm->gputype.getGPUtype(prhs[0]);

  double **srcs = getGPUSourcePointers(prhs, 10, &numel, 0, gm);
  double **dest = makeGPUDestinationArrays(srcReference,  plhs, 5, gm);

  // Establish launch dimensions & a few other parameters
  int fluxDirection = (int)*mxGetPr(prhs[11]);
  double lambda     = *mxGetPr(prhs[10]);
  int numDims        = gm->gputype.getNdims(srcReference);
  const int *dims    = gm->gputype.getSize(srcReference);

  dim3 arraySize;
  arraySize.x = dims[0];
  numDims > 1 ? arraySize.y = dims[1] : arraySize.y = 1;
  numDims > 2 ? arraySize.z = dims[2] : arraySize.z = 1;

  dim3 blocksize, gridsize;
  int hu, hv, hw, nu;

  blocksize.x = BLOCKLEN+2; blocksize.y = blocksize.z = 1;
  switch(fluxDirection) {
    case 1: // X direction flux: u = x, v = y, w = z;
      gridsize.x = arraySize.y;
      gridsize.y = arraySize.z;
      hu = 1; hv = arraySize.x; hw = arraySize.x * arraySize.x;
      nu = arraySize.x; break;
    case 2: // Y direction flux: u = y, v = x, w = z
      gridsize.x = arraySize.x;
      gridsize.y = arraySize.z;
      hu = arraySize.x; hv = 1; hw = arraySize.x * arraySize.y;
      nu = arraySize.y; break;
    case 3: // Z direction flux: u = z, v = x, w = y;
      gridsize.x = arraySize.x;
      gridsize.y = arraySize.y;
      hu = arraySize.x * arraySize.y; hv = 1; hw = arraySize.x;
      nu = arraySize.z; break;
    }

fluidVarPtrs fluid;
int i;
for(i = 0; i < 5; i++) { fluid.fluidIn[i] = srcs[i]; fluid.fluidOut[i] = dest[i]; }
fluid.B[0] = srcs[5];
fluid.B[1] = srcs[6];
fluid.B[2] = srcs[7];

fluid.Ptotal = srcs[8];
fluid.cFreeze = srcs[9];

if(nu > 1) {
  cudaWFluxKernel<<<gridsize, blocksize>>>(fluid, lambda, nu, hu, hv, hw, fluxDirection);
  } else {
  nullStep<<<32, 128>>>(fluid, numel);
  }

}


__global__ void cudaWFluxKernel(fluidVarPtrs fluid, double lambda, int nu, int hu, int hv, int hw, int direction)
{

// This is the shared flux array.
// It is long enough to hold both the left and right fluxes for a variable, with extra spaces that
// strategically leave the next flux variables for the left edge behind when our fLeft pointer is decremented
__shared__ double fluxArrayL[BLOCKLEN+6];
__shared__ double fluxArrayR[BLOCKLEN+6];

// Local copies of the hyperbolic state variables: The conserved quantities Qi (rho, energy,
// momentum), a W variable, and local copies of the pressure, inverse freezing speed and B
double Qi[5], Wi, lPress, Cinv, lBx, lBy, lBz;

// IndexBase identifies where we go into the array, i is our generic counter variable
int indexBase;// = blockIdx.x * hv + blockIdx.y * hw + ((threadIdx.x-1)%nu)* hu;
//int endAddress = blockIdx.x * hv + blockIdx.y * hw + nu*hu;
int i;

// Pointer into the shared flux array that's strategically decremented to leave the variables we need next behind
double *fLeft; double *fRight;
double rhoinv;

// Loop index - controls the advancement of the thread block through the line; Doflux tells threads to flux or not flux
int lpIdx;

for(lpIdx = threadIdx.x-1; lpIdx <= nu; lpIdx += BLOCKLEN) {
	// Our index in the U direction; If it exceeds the array size+1, return; If it equals it, don't flux that cell

	// Convert i to an index; Circular boundary conditions
	indexBase = blockIdx.x * hv + blockIdx.y * hw + (lpIdx % nu)* hu;
if(lpIdx < 0) indexBase = blockIdx.x * hv + blockIdx.y * hw + (nu-1)* hu;

	// Load local copies of hyperbolic state variables
	for(i = 0; i < 5; i++) { Qi[i] = fluid.fluidIn[i][indexBase]; }
	lPress = fluid.Ptotal[indexBase];
	Cinv   = 1.0 / fluid.cFreeze[indexBase];
	rhoinv = 1.0 / Qi[0];
	lBx = fluid.B[0][indexBase];
	lBy = fluid.B[1][indexBase];
	lBz = fluid.B[2][indexBase];

	// Start outselves off with 4 chances to move left
	fLeft  = fluxArrayL + 4;
	fRight = fluxArrayR + 4;

	// For each conserved quantity,
	for(i = 0; i < 5; i++) {
		// Calculate the W flux
		switch(5*direction + i) {
			case 5: Wi = Qi[2] * Cinv; break;
			case 6: Wi = (Qi[2] * (Qi[1] + lPress) - lBx*(Qi[2]*lBx+Qi[3]*lBy+Qi[4]*lBz) ) * Cinv * rhoinv; break;
			case 7: Wi = (Qi[2]*Qi[2]*rhoinv + lPress - lBx*lBx)*Cinv; break;
			case 8: Wi = (Qi[2]*Qi[3]*rhoinv          - lBx*lBy)*Cinv; break;
			case 9: Wi = (Qi[2]*Qi[4]*rhoinv          - lBx*lBz)*Cinv; break;

			case 10: Wi = Qi[3] * Cinv; break;
			case 11: Wi = (Qi[3] * (Qi[1] + lPress) - lBy*(Qi[2]*lBx+Qi[3]*lBy+Qi[4]*lBz) ) * Cinv * rhoinv; break;
			case 12: Wi = (Qi[3]*Qi[2]*rhoinv          - lBy*lBx)*Cinv; break;
			case 13: Wi = (Qi[3]*Qi[3]*rhoinv + lPress - lBy*lBy)*Cinv; break;
			case 14: Wi = (Qi[3]*Qi[4]*rhoinv          - lBy*lBz)*Cinv; break;
  
			case 15: Wi = Qi[4] * Cinv; break;
			case 16: Wi = (Qi[4] * (Qi[1] + lPress) - lBz*(Qi[2]*lBx+Qi[3]*lBy+Qi[4]*lBz) ) * Cinv * rhoinv; break;
			case 17: Wi = (Qi[4]*Qi[2]*rhoinv          - lBz*lBx)*Cinv; break;
			case 18: Wi = (Qi[4]*Qi[3]*rhoinv          - lBz*lBy)*Cinv; break;
			case 19: Wi = (Qi[4]*Qi[4]*rhoinv + lPress - lBz*lBz)*Cinv; break;
			}

		// Decouple into left & right going fluxes
		fLeft[threadIdx.x]  = .25*(Qi[i] - Wi);
		fRight[threadIdx.x] = .25*(Qi[i] + Wi);
//              fLeft[threadIdx.x]  = 5.0;
//              fRight[threadIdx.x] = 4;


		// Stop until all threads finish filling flux array
		__syncthreads();

		// Write updated quantities to global memory
//		if((threadIdx.x > 0) && (threadIdx.x <= BLOCKLEN)) fluid.fluidOut[i][indexBase] = Qi[i] - lambda*( fLeft[threadIdx.x] - fLeft[threadIdx.x+1] + fRight[threadIdx.x] - fRight[threadIdx.x-1] )/Cinv;
		if((lpIdx >= 0) && (lpIdx < nu)) fluid.fluidOut[i][indexBase] = Qi[i] - lambda*( fLeft[threadIdx.x] - fLeft[threadIdx.x+1] + fRight[threadIdx.x] - fRight[threadIdx.x-1] )/Cinv;

		// Make sure everyone's done reading shared memory before we futz with it again; Move left
		__syncthreads();
		fLeft--;
		fRight--;
		}

;
	// Thread zero is only to fill the left edge on the very first iteration - break
	if(threadIdx.x == 0) break;

	if(threadIdx.x < 6) {
		fluxArrayL[threadIdx.x-1] = fluxArrayL[threadIdx.x + BLOCKLEN];
		fluxArrayR[threadIdx.x-1] = fluxArrayR[threadIdx.x + BLOCKLEN];
		}

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

