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
static CUfunction drvfuns[4];
static int init = 0;
static GPUmat *gm;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // At least 2 arguments expected
  // Input and result
  if ((nrhs!=2) || (nlhs != 1))
     mexErrMsgTxt("Wrong number of arguments");
  if (init == 0) {
    gm = gmGetGPUmat();
    init = 1;
  }
  /* mex parameters are:
   0 array 1
   1 array 2
   This function is a wrapper for the CUBLAS double-prec dot product function pending newer matcuda
   */
  GPUtype arrayA = gm->gputype.getGPUtype(prhs[0]);
  GPUtype arrayB = gm->gputype.getGPUtype(prhs[1]);

  int numElements = gm->gputype.getNumel(arrayA);
  if (numElements != gm->gputype.getNumel(arrayB)) mexErrMsgTxt("Arrays contain different numbers of elements.\n");

  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL); //mxReal is our data-type
  double *AdotB = mxGetPr(plhs[0]);

  void *pointerA = (void *)gm->gputype.getGPUptr(arrayA);
  void *pointerB = (void *)gm->gputype.getGPUptr(arrayB);

  double *u = (double*)pointerA;
  double *v = (double *)pointerB; 

  AdotB[0] = cublasDdot(numElements, u, 1, v, 1);

}
