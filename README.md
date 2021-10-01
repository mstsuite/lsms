# LSMS
LSMS is a code for scalable first principles calculations of materials using multiple scattering theory.

## Citing LSMS
If you publish results obtained using LSMS we ask that you cite the following publications:

* Y. Wang, G. M. Stocks, W. A. Shelton, D. M. C. Nicholson, W. M. Temmerman, and Z. Szotek. _Order-n multiple scattering approach to electronic structure calculations_. Phys. Rev. Lett. **75**, 2867 (1995).

and if the GPU accelerated version was used, please cite additionally:

* M. Eisenbach, J. Larkin, J. Lutjens, S. Rennich, and J. H. Rogers. _GPU acceleration of the locally selfconsistent multiple scattering code for first principles calculation of the ground state and statistical physics of materials_. Computer Physics Communications **211**, 2 (2017).

and for calculations using Monte-Carlo simulations:

* M. Eisenbach, C.-G. Zhou, D. M. C. Nicholson, G. Brown, J. Larkin, and T. C. Schulthess. _A Scalable Method for Ab Initio Computation of Free Energies in Nanoscale Systems_. Proceedings of the Conference on High Performance Computing Networking, Storage and Analysis, ACM, New York, 64 (2009)

## Installation

[![CMake](https://github.com/francoMU/lsms/actions/workflows/cmake.yml/badge.svg)](https://github.com/francoMU/lsms/actions/workflows/cmake.yml)

### CMake

<div style="padding: 15px; border: 1px solid transparent; border-color: transparent; margin-bottom: 20px; border-radius: 4px; color: #a94442; background-color: #f2dede; border-color: #ebccd1;">
Warning

The current CMake can only build LSMS and not WL-LSMS

</div>

#### Prerequisites

- CMake: >=3.18
- make
- Autotools: (Optional) This is needed if LibXC is desired to be compiled automatically

The following libraries needs to be preinstalled

- HDF5
- LAPACK/BLAS or MKL
- MPI
- OpenMP

The following libaries can be preinstalled, but can also be installed automatically

- LibXC
- Lua


#### Examples

Here is an example of installing `LSMS` with GCC-10, OpenMPI and MKL. LibXC and Lua will be installed during the build. 

Clone the branch with the CMake build system

```bash
git clone -b feature/build_system git@github.com:francoMU/lsms.git
```

Create a build folder in a directory of your choice.

```bash
mkdir build_lsms
cd build_lsms
```

Run Cmake. CMake fill try to find all libaries automatically.

```bash
cmake <path/to/lsms/root> \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_COMPILER=mpic++ \
      -DCMAKE_C_COMPILER=gcc \
      -DCMAKE_CXX_COMPILER=gfortran \
      -DCMAKE_CXX_FLAGS="-O3 -mtune=native" \
      -DCMAKE_Fortran_FLAGS="-O3 -mtune=native -fbacktrace -cpp -fallow-argument-mismatch" \
      -DBLA_VENDOR=Intel10_64lp
```

This is just an example for creating the build-system. `CMAKE_BUILD_TYPE` can be either Release or Debug. This has only an
effect if no compiler flags are specified. In this case cmake will preset some flags depending on the chosen 
build type. It is necessary to specify a C++,C and a Fortran Compiler. Preferably, one specifies as a C++ the compiler 
wrapper, because it contains all the necessary flags and path already and it will be easier for CMake to choose the right
libraries. In case of the `gcc` toolchain, one should prefer the compiler wrapper mpic++ over g++. 
`CMAKE_CXX_FLAGS` and `CMAKE_Fortran_FLAGS` are the compiler flags that will be used to build the object files. 
OpenMP will be search automatically by CMake and corresponding flags are added automatically.
`LSMS` also needs LAPACK and BLAS. One can specify the desired type of LAPACK with the `BLA_VENDOR` definition. 
The corresponding documentation can be found here ([CMake](https://cmake.org/cmake/help/v3.18/module/FindLAPACK.html)).

Build all target.

```bash
cmake --build . --parallel
```

The `LSMS` binary can then be found in the build directory under ```bin```.
