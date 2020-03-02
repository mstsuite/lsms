OBJ = block_inverse_fortran.o block_inverse.o zblock_lu_cpp.o zblock_lu_CPU.o
OBJ_CUDA = inversionTest_cuda.o


ifeq ($(ARCH),macos)
  LIBS = -lgfortran -framework Accelerate
  CXX = g++-9
  F77 = gfortran
else
ifeq ($(ARCH),ubuntu-cuda)
  LIBS = -L/usr/local/cuda/lib64 -lcudart -lcublas -L/usr/local/lib -lblas -llapack -lgfortran
  OBJ += $(OBJ_CUDA)
  CXX = g++
  F77 = gfortran
  NVCC = nvcc
else
  LIBS = -L/usr/local/lib -lblas -llapack -lgfortran
  CXX = g++
  F77 = gfortran
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

