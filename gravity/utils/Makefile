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
all: gravity fluid

gravity: accumulateBterm symmetricLinearOperator wrap_cublasDdot integralPoisson_mg interpolate
fluid: cudaArrayAtomic fluxLimiterKernels directionalMaxFinder cudaMHDKernels cudaArrayRotate cudaWstep cudaTVDStep cudaMagW cudaMagFull freezeAndPtot cudaShift

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

cudaArrayAtomic: cudaArrayAtomic.cu cudaCommon
	nvcc $(NVCC_OPTS) $(NVCC_INCL) -cuda cudaArrayAtomic.cu -o cudaArrayAtomic.cpp
	mex $(MEX_OPTS) $(MEX_INCL) cudaArrayAtomic.cpp cudaCommon.cpp $(GPUMAT_CPP) -o cudaArrayAtomic $(MEX_LDIR) $(MEX_LIBS)
	rm -f cudaArrayAtomic.cpp

fluxLimiterKernels: fluxLimiterKernels.cu cudaCommon
	nvcc $(NVCC_OPTS) $(NVCC_INCL) -cuda fluxLimiterKernels.cu -o fluxLimiterKernels.cpp
	mex $(MEX_OPTS) $(MEX_INCL) fluxLimiterKernels.cpp cudaCommon.cpp $(GPUMAT_CPP) -o fluxLimiterKernels $(MEX_LDIR) $(MEX_LIBS)
	rm -f fluxLimiterKernels.cpp

directionalMaxFinder: directionalMaxFinder.cu cudaCommon
	nvcc $(NVCC_OPTS) $(NVCC_INCL) -cuda directionalMaxFinder.cu -o directionalMaxFinder.cpp
	mex $(MEX_OPTS) $(MEX_INCL) directionalMaxFinder.cpp cudaCommon.cpp $(GPUMAT_CPP) -o directionalMaxFinder $(MEX_LDIR) $(MEX_LIBS)
	rm -f directionalMaxFinder.cpp

cudaArrayRotate: cudaArrayRotate.cu cudaCommon
	nvcc $(NVCC_OPTS) $(NVCC_INCL) -cuda cudaArrayRotate.cu -o cudaArrayRotate.cpp
	mex $(MEX_OPTS) $(MEX_INCL) cudaArrayRotate.cpp cudaCommon.cpp $(GPUMAT_CPP) -o cudaArrayRotate $(MEX_LDIR) $(MEX_LIBS)
	rm -f cudaArrayRotate.cpp

cudaShift: cudaShift.cu cudaCommon
	nvcc $(NVCC_OPTS) $(NVCC_INCL) -cuda cudaShift.cu -o cudaShift.cpp
	mex $(MEX_OPTS) $(MEX_INCL) cudaShift.cpp cudaCommon.cpp $(GPUMAT_CPP) -o cudaShift $(MEX_LDIR) $(MEX_LIBS)
	rm -f cudaShift.cpp

cudaStatics: cudaStatics.cu cudaCommon
	nvcc $(NVCC_OPTS) $(NVCC_INCL) -cuda cudaStatics.cu -o cudaStatics.cpp
	mex $(MEX_OPTS) $(MEX_INCL) cudaStatics.cpp cudaCommon.cpp $(GPUMAT_CPP) -o cudaStatics $(MEX_LDIR) $(MEX_LIBS)
	rm -f cudaStatics.cpp

cudaMHDKernels: cudaMHDKernels.cu cudaCommon
	nvcc $(NVCC_OPTS) $(NVCC_INCL) -cuda cudaMHDKernels.cu -o cudaMHDKernels.cpp
	mex $(MEX_OPTS) $(MEX_INCL) cudaMHDKernels.cpp cudaCommon.cpp $(GPUMAT_CPP) -o cudaMHDKernels $(MEX_LDIR) $(MEX_LIBS)
	rm -f cudaMHDKernels.cpp

cudaLinearSubsagn: cudaLinearSubsagn.cu
	nvcc $(NVCC_OPTS) $(NVCC_INCL) -cuda cudaLinearSubsagn.cu -o cudaLinearSubsagn.cpp
	mex $(MEX_OPTS) $(MEX_INCL) cudaLinearSubsagn.cpp $(GPUMAT_CPP) -o cudaLinearSubsagn $(MEX_LDIR) $(MEX_LIBS)
	rm -f cudaLinearSubsagn.cpp

freezeAndPtot: freezeAndPtot.cu cudaCommon.cu
	nvcc $(NVCC_OPTS) $(NVCC_INCL) -cuda freezeAndPtot.cu -o freezeAndPtot.cpp
	mex $(MEX_OPTS) $(MEX_INCL) freezeAndPtot.cpp cudaCommon.cpp $(GPUMAT_CPP) -o freezeAndPtot $(MEX_LDIR) $(MEX_LIBS)
	rm -f freezeAndPtot.cpp

cudaWstep: cudaWstep.cu cudaCommon
	nvcc $(NVCC_OPTS) $(NVCC_INCL) -cuda cudaWstep.cu -o cudaWstep.cpp
	mex $(MEX_OPTS) $(MEX_INCL) cudaWstep.cpp cudaCommon.cpp $(GPUMAT_CPP) -o cudaWstep $(MEX_LDIR) $(MEX_LIBS)
	rm -f cudaWstep.cpp

cudaTVDStep: cudaTVDStep.cu cudaCommon
	nvcc $(NVCC_OPTS) $(NVCC_INCL) -cuda cudaTVDStep.cu -o cudaTVDStep.cpp
	mex $(MEX_OPTS) $(MEX_INCL) cudaTVDStep.cpp cudaCommon.cpp $(GPUMAT_CPP) -o cudaTVDStep $(MEX_LDIR) $(MEX_LIBS)
	rm -f cudaTVDStep.cpp

cudaMagW: cudaMagW.cu cudaCommon
	nvcc $(NVCC_OPTS) $(NVCC_INCL) -cuda cudaMagW.cu -o cudaMagW.cpp
	mex $(MEX_OPTS) $(MEX_INCL) cudaMagW.cpp cudaCommon.cpp $(GPUMAT_CPP) -o cudaMagW $(MEX_LDIR) $(MEX_LIBS)
	rm -f cudaMagW.cpp

cudaMagFull: cudaMagFull.cu cudaCommon
	nvcc $(NVCC_OPTS) $(NVCC_INCL) -cuda cudaMagFull.cu -o cudaMagFull.cpp
	mex $(MEX_OPTS) $(MEX_INCL) cudaMagFull.cpp cudaCommon.cpp $(GPUMAT_CPP) -o cudaMagFull $(MEX_LDIR) $(MEX_LIBS)
	rm -f cudaMagFull.cpp

cmaxForTimestep: cmaxForTimestep.cu cudaCommon
	nvcc $(NVCC_OPTS) $(NVCC_INCL) -cuda cmaxForTimestep.cu -o cmaxForTimestep.cpp
	mex $(MEX_OPTS) $(MEX_INCL) cmaxForTimestep.cpp cudaCommon.cpp $(GPUMAT_CPP) -o cmaxForTimestep $(MEX_LDIR) $(MEX_LIBS)
	rm -f cmaxForTimestep.cpp

cudaCommon: cudaCommon.h cudaCommon.cu
	nvcc $(NVCC_OPTS) $(NVCC_INCL) -cuda cudaCommon.cu -o cudaCommon.cpp

clean:
	rm -f *mexa64 cudaKernels.cpp
