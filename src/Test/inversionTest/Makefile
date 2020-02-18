OBJ = block_inverse_fortran.o block_inverse.o zblock_lu_cpp.o zblock_lu_CPU.o

LIBS = -lgfortran -framework Accelerate

CXX = g++-9
F77 = gfortran

all: inversionTest

clean:
	rm -rf *.o inversionTest

inversionTest: inversionTest.cpp $(OBJ)
	$(CXX) -o inversionTest inversionTest.cpp $(OBJ) $(LIBS)

%.o : %.cpp
	$(CXX) $(INC_PATH) -c -o $@ $<

%.o : %.f
	$(F77) -c -o $@ $<

