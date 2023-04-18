all: xcTest

XCLIB_INC = ../../../opt/include
XCLIB_LIB = ../../../opt/lib

xcTest: xcTest.cpp
	g++ -o xcTest -I$(XCLIB_INC) xcTest.cpp -L$(XCLIB_LIB) -lxc -lm

clean:
	rm -f xcTest
