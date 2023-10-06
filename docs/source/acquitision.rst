************
Getting LSMS
************

The code can be obtained from the github repo (https://github.com/mstsuite/lsms)

.. parsed-literal::
   git clone https://github.com/mstsuite/lsms.git

LSMS has the following dependencies

1. Fortran and C++
2. Cmake
3. HDF5
4. BLAS and LAPACK
5. MPI
6. (Optional) CUDA/ROCm for GPU acceleration

Lua and LibXC are included - if you would like to use a pre-existing lua or LibXC installation it can be linked during the build process
