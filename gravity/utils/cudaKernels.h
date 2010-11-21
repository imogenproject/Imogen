#define EDGEDIM_BOP 8
#define EDGEDIM_MGBC 8

__global__ void Laplacian_B_OperatorKernel(double *src, double *dst, double accumCoefficient, double *accum, int nx, int ny, int nz, int level);
__global__ void mgbc_kernelplane(double **mass, double **position, short int *arrayDims, int nLevels, double *phi,  double *bvec, double qConst, double selfPotRad);
__global__ void SymmetricOperatorKernel(double *src, double *dst, int nx, int ny, int nz, double c1, double c2, double c3, double c4);
__global__ void   upsampleKernel(double *src, double *dst, int factor, int nx, int ny, int nz);
__global__ void downsampleKernel(double *src, double *dst, int factor, int nx, int ny, int nz);

