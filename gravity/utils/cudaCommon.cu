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

#include "cudaCommon.h"

/* Given the RHS and how many cuda arrays we expect, extracts a set of pointers to GPU memory for us
 Also conveniently checked for equal array extent and returns it for us */
double **getGPUSourcePointers(const mxArray *prhs[], int num, int *retNumel, int startat, GPUmat *gm)
{
  GPUtype src;
  double **gpuPointers = (double **)malloc(num * sizeof(double *));
  int iter;
  int numel = gm->gputype.getNumel(gm->gputype.getGPUtype(prhs[startat]));

  for(iter = 0; iter < num; iter++) {
    src = gm->gputype.getGPUtype(prhs[startat + iter]);
/*    if (gm->gputype.getNumel(src) != numel) { free(gpuPointers); mexErrMsgTxt("Fatal: Arrays contain nonequal number of elements."); } */
    gpuPointers[iter] = (double *)gm->gputype.getGPUptr(src);
  }

retNumel[0] = numel;
return gpuPointers;
}

/* Creates destination array that the kernels write to; Returns the GPU memory pointer, and assigns the LHS it's passed */
double **makeGPUDestinationArrays(GPUtype src, mxArray *retArray[], int howmany, GPUmat *gm)
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

}

