*****
Installing LSMS
*****

Here is an example of installing `LSMS` with GCC-10, OpenMPI and MKL. LibXC and Lua will be installed during the build. 

Clone the branch with the CMake build system

.. parsed-literal::

   git clone https://github.com/mstsuite/lsms.git


Create a build folder in a directory of your choice.

.. parsed-literal::

   mkdir build_lsms
   cd build_lsms


Run CMake. CMake will try to find all libaries automatically.

.. parsed-literal::

   cmake <path/to/lsms/root> \
   -DCMAKE_BUILD_TYPE=Release \
   -DCMAKE_CXX_COMPILER=mpic++ \
   -DCMAKE_C_COMPILER=gcc \
   -DCMAKE_Fortran_COMPILER=gfortran \
   -DCMAKE_CXX_FLAGS="-O3 -mtune=native" \
   -DCMAKE_Fortran_FLAGS="-O3 -mtune=native -fbacktrace -cpp -fallow-argument-mismatch" \
   -DBLA_VENDOR=Intel10_64lp


For a number of systems toolchain files are provided in the `toolchain` directory that provide the definitions to run CMake. Here is what the toolchain file (in toolchain/frontier-cray-hip.cmake) looks like for building on the Frontier supercomputer at OLCF using a combination of Cray and AMD compilers.

.. parsed-literal::
   Toolchain for building LSMS with HIP on OLCF Frontier with amd compiler

   set(CMAKE_BUILD_TYPE Release)
   set(BUILD_TESTING OFF)

   set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)

   set(SEARCH_LAPACK OFF)
   set(SEARCH_BLAS OFF)
   set(BLAS_LIBRARIES "$ENV{CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_cray.a")
   set(LAPACK_LIBRARIES "$ENV{CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_cray.a")
   set(USE_ACCELERATOR_HIP ON)
   set(MST_LINEAR_SOLVER_DEFAULT 0x0020)
   #set(MST_LINEAR_SOLVER_DEFAULT 0x0005)
   set(MST_BUILD_KKR_MATRIX_DEFAULT 0x3000)
   #set(MST_BUILD_KKR_MATRIX_DEFAULT 0x1000)
   set(AMDGPU_TARGETS "gfx90a")
   set(GPU_TARGETS "gfx90a")
   set(CMAKE_CXX_FLAGS "-std=c++17 --gcc-toolchain=/opt/cray/pe/gcc/12.2.0/snos -L/opt/cray/pe/gcc/12.2.0/snos/lib64 -mtune=native")

   # Loaded Modules:
   # module load PrgEnv-cray-amd
   # module load cray-hdf5
   # module load cmake libtool
   # module load rocm

The toolchain file makes it easy to keep track of all the build options. To use a toolchain file, you need to include it in the CMake command

.. parsed-literal::
   cmake -DCMAKE_TOOLCHAIN_FILE=<path/to/lsms/root>/toolchain/<system-name>.cmake \
   <path/to/lsms/root>

This is just an example for creating the build-system. `CMAKE_BUILD_TYPE` can be either Release or Debug. This has only an
effect if no compiler flags are specified. In this case cmake will preset some flags depending on the chosen 
build type. It is necessary to specify a C++,C and a Fortran Compiler. The C++ compiler should be preferably be specified
as the compiler wrapper, because it contains all the necessary flags and path already and it will be easier for CMake to choose the right
libraries. For the  `GNU` toolchain for example, one should prefer the compiler wrapper mpic++ over g++.
`CMAKE_CXX_FLAGS` and `CMAKE_Fortran_FLAGS` are the compiler flags that will be used to build the object files. 
OpenMP will be search automatically by CMake and corresponding flags are added automatically.
`LSMS` also needs LAPACK and BLAS. One can specify the desired type of LAPACK with the `BLA_VENDOR` definition. 
The corresponding documentation can be found here ([CMake](https://cmake.org/cmake/help/v3.18/module/FindLAPACK.html)).

Sometimes it is necessary to explicitly define the path to math libaries. In this case the automatic search has
to be turned of by setting `SEARCH_LAPACK` or `SEARCH_BLAS` to `OFF`. The libaries has to be specified explicitly, if 
the search is turned off.

.. parsed-literal::
   cmake <path/to/lsms/root> \
   -DCMAKE_BUILD_TYPE=Release \
   -DCMAKE_CXX_COMPILER=mpic++ \
   -DCMAKE_C_COMPILER=gcc \
   -DCMAKE_Fortran_COMPILER=gfortran \
   -DCMAKE_CXX_FLAGS="-O3 -mtune=native" \
   -DCMAKE_Fortran_FLAGS="-O3 -mtune=native -fbacktrace -cpp -fallow-argument-mismatch" \
   -DSEARCH_LAPACK=OFF \
   -DSEARCH_BLAS=OFF \
   -DLAPACK_LIBARIES=liblapack.so \
   -DBLAS_LIBARIES=libblas.so


The path to the LibXC and Lua libraries can be specified explicitly. In this case the libraries will not be installed
automatically.

.. parsed-literal::
   cmake <path/to/lsms/root> \
   -DCMAKE_BUILD_TYPE=Release \
   -DCMAKE_CXX_COMPILER=mpic++ \
   -DCMAKE_C_COMPILER=gcc \
   -DCMAKE_Fortran_COMPILER=gfortran \
   -DCMAKE_CXX_FLAGS="-O3 -mtune=native" \
   -DCMAKE_Fortran_FLAGS="-O3 -mtune=native -fbacktrace -cpp -fallow-argument-mismatch" \
   -Dlibxc_LIBRARIES=<path/libxc>/lib/libxc.a
   -Dlibxc_INCLUDE_DIR=<path/libxc>/include

Build all target.

.. parsed-literal::
   cmake --build . --parallel

The `LSMS` binary can then be found in the build directory under the bin folder.
