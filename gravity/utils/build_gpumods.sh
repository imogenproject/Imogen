# This script compiles the handful of binary modules needed for Imogen to support GPU gravitation.
# These must be rebuilt EVERY TIME GPUMAT IS UPDATED or Matlab will segfault when they are used.

# Change these to suit your system. CUDA_LDIR will be probably be lib (or lib64 on mixed 32/64 bit systems)
GPUMAT_DIR="$HOME/GPUmat"
CUDA_DIR="/usr/local/cuda"
CUDA_LDIR="lib64"

GPUMAT_INCL="$GPUMAT_DIR/modules/include"
GPUMAT_CPP="$GPUMAT_DIR/modules/common/GPUmat.cpp"

# Build CUBLAS wrappers for linear solver
mex -DUNIX -I$CUDA_DIR/include -I$GPUMAT_INCL -L$CUDA_DIR/$CUDA_LDIR -lcuda -lcublas shiftAccumDaxpy.cpp $GPUMAT_CPP
mex -DUNIX -I$CUDA_DIR/include -I$GPUMAT_INCL -L$CUDA_DIR/$CUDA_LDIR -lcuda -lcublas wrap_cublasDdot.cpp $GPUMAT_CPP

# Build boundary condition routine
nvcc -arch sm_13 -I/opt/Matlab/Matlab2010a/extern/include -I$GPUMAT_INCL -cuda mgbc_gpukern.cu -o mgbc_gpukern.cpp
mex -DUNIX -I$CUDA_DIR/include -I$GPUMAT_INCL -L$CUDA_DIR/$CUDA_LDIR mgbc_gpukern.cpp $GPUMAT_CPP -o mgbc_gpukern -lcuda -lcudart

rm -f mgbc_gpukern.cpp
