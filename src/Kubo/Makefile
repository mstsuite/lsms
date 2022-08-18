test: CurrentMatrix.o testCurrentMatrix.o
	clang++ ../Misc/ClebschGordan.o SingleSiteScattering.o CurrentMatrix.o testCurrentMatrix.o -o test

CurrentMatrix.o: CurrentMatrix.cpp
	clang++ -I ../ -I ../../include/ -std=c++14 -o CurrentMatrix.o -c CurrentMatrix.cpp

testCurrentMatrix.o: testCurrentMatrix.cpp
	clang++ -I ../ -I ../../include/ -std=c++14 -o testCurrentMatrix.o -c testCurrentMatrix.cpp
