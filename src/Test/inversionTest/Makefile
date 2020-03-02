OBJ = block_inverse_fortran.o block_inverse.o zblock_lu_cpp.o zblock_lu_CPU.o
OBJ_CUDA = inversionTest_cuda.o


ifeq ($(ARCH),macos)
  LIBS = -lgfortran -framework Accelerate
  CXX = g++-9 -DARCH_MACOS
  F77 = gfortran
else
ifeq ($(ARCH),ubuntu-cuda)
  LIBS = -L/usr/local/cuda/lib64 -lcudart -lcublas -L/usr/local/lib -lblas -llapack -lgfortran
  OBJ += $(OBJ_CUDA)
  CXX = g++ -DARCH_CUDA
  F77 = gfortran
  NVCC = nvcc -DARCH_CUDA
else
ifeq ($(ARCH),summit)
  LIBS = -L$(OLCF_CUDA_ROOT)/lib64 -lcudart -lcublas -L${OLCF_ESSL_ROOT}/lib64 -lessl -lgfortran
  OBJ += $(OBJ_CUDA)
  CXX = g++ -DARCH_IBM -DARCH_CUDA
  F77 = gfortran
  NVCC = nvcc -DARCH_IBM -DARCH_CUDA
else
  LIBS = -L/usr/local/lib -lblas -llapack -lgfortran
  CXX = g++ -DARCH_GENERIC
  F77 = gfortran
endif
endif
endif

# CXX = g++-9
# F77 = gfortran
# NVCC = nvcc

all: inversionTest

clean:
	rm -rf *.o inversionTest

inversionTest: inversionTest.cpp $(OBJ)
	$(CXX) -o inversionTest inversionTest.cpp $(OBJ) $(LIBS)

%.o : %.cpp
	$(CXX) $(INC_PATH) -c -o $@ $<

%.o : %.f
	$(F77) -c -o $@ $<

%.o : %.cu
	$(NVCC) -c -o $@ $<

