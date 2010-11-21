# This script compiles the handful of binary modules needed for Imogen to support GPU gravitation.
# These must be rebuilt EVERY TIME GPUMAT IS UPDATED or Matlab will segfault when they are used.

# Change these to suit your system. CUDA_LDIR will be either lib or lib64
GPUMAT_DIR=/home/erik/GPUmat
MATLAB_DIR=/opt/Matlab/matlab2010b
CUDA_DIR = /usr/local/cuda
CUDA_LDIR = lib64
CUDA_ARCH = sm_13

# Don't mess with these
GPUMAT_INCL = -I$(GPUMAT_DIR)/modules/include
GPUMAT_CPP = $(GPUMAT_DIR)/modules/common/GPUmat.cpp

NVCC_OPTS = -arch $(CUDA_ARCH)
NVCC_INCL = -I$(MATLAB_DIR)/extern/include $(GPUMAT_INCL)

MEX_OPTS  = -DUNIX
MEX_INCL  = -I$(CUDA_DIR)/include $(GPUMAT_INCL) -I$(MATLAB_DIR)/extern/include
MEX_LDIR  = -L$(CUDA_DIR)/$(CUDA_LDIR)
MEX_LIBS  = -lcuda -lcudart -lcublas

###############################################################################
all: accumulateBterm symmetricLinearOperator wrap_cublasDdot integralPoisson_mg interpolate


kernels: cudaKernels.cu
	nvcc $(NVCC_OPTS) $(NVCC_INCL) -cuda cudaKernels.cu -o cudaKernels.cpp

# Build preconditioning operator
accumulateBterm: accumulateBterm.cu kernels
	nvcc $(NVCC_OPTS) $(NVCC_INCL) -cuda accumulateBterm.cu -o accumulateBterm.cpp
	mex $(MEX_OPTS) $(MEX_INCL) accumulateBterm.cpp $(GPUMAT_CPP) cudaKernels.cpp -o accumulateBterm $(MEX_LDIR) $(MEX_LIBS)
	rm -f accumulateBterm.cpp

# Build general symmetric linear operator function
symmetricLinearOperator: symmetricLinearOperator.cu kernels
	nvcc $(NVCC_OPTS) $(NVCC_INCL) -cuda symmetricLinearOperator.cu -o symmetricLinearOperator.cpp
	mex $(MEX_OPTS) $(MEX_INCL) symmetricLinearOperator.cpp $(GPUMAT_CPP) cudaKernels.cpp -o symmetricLinearOperator $(MEX_LDIR) $(MEX_LIBS)
	rm -f symmetricLinearOperator.cpp

interpolate: interpolateGPUvar.cu kernels
	nvcc $(NVCC_OPTS) $(NVCC_INCL) -cuda interpolateGPUvar.cu -o interpolateGPUvar.cpp
	mex $(MEX_OPTS) $(MEX_INCL) interpolateGPUvar.cpp $(GPUMAT_CPP) cudaKernels.cpp -o interpolateGPUvar $(MEX_LDIR) $(MEX_LIBS)
	rm -f interpolateGPUvar.cpp

# Build dot product operator
wrap_cublasDdot: wrap_cublasDdot.cpp
	mex $(MEX_OPTS) $(MEX_INCL) wrap_cublasDdot.cpp $(GPUMAT_CPP) $(MEX_LDIR) $(MEX_LIBS)

# Build boundary condition routine
mgbc_genkern: mgbc_genkern.cu
	nvcc $(NVCC_OPTS) $(NVCC_INCL) -cuda mgbc_genkern.cu -o mgbc_genkern.cpp
	mex $(MEX_OPTS) $(MEX_INCL) mgbc_genkern.cpp $(GPUMAT_CPP) -o mgbc_genkern $(MEX_LDIR) $(MEX_LIBS)
	rm -f mgbc_genkern.cpp

integralPoisson_mg: integralPoisson_mg.cu kernels
	nvcc $(NVCC_OPTS) $(NVCC_INCL) -cuda integralPoisson_mg.cu -o integralPoisson_mg.cpp
	mex $(MEX_OPTS) $(MEX_INCL) integralPoisson_mg.cpp $(GPUMAT_CPP) cudaKernels.cpp -o integralPoisson_mg $(MEX_LDIR) $(MEX_LIBS)
	rm -f integralPoisson_mg.cpp

clean:
	rm -f *mexa64 cudaKernels.cpp
