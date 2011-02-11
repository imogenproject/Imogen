typedef struct {
        double *fluidIn[5];
        double *fluidOut[5];

        double *B[3];

        double *Ptotal;
        double *cFreeze;
        } fluidVarPtrs;

double **getGPUSourcePointers(const mxArray *prhs[], int num, int *retNumel, int startat, GPUmat *gm);
double **makeGPUDestinationArrays(GPUtype src, mxArray *retArray[], int howmany, GPUmat *gm);

