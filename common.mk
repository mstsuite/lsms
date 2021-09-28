
# C compiler
CC := gcc
# C++ compiler
CXX := mpic++
# Fortran compiler
F90 := gfortran
# linker
LD := g++
# C flags
CFLAGS := -std=c11
# Fortran flags
FFLAGS := -g -fbacktrace
# Fortran preprocessor flags
FPPFLAGS := -cpp
# C++ flags
CXXFLAGS := -std=c++14
# C/C++ flags
CPPFLAGS := -g
# linker flags
LDFLAGS :=
# linker flags: libraries to link (e.g. -lfoo)
LDLIBS :=
# flags required for dependency generation; passed to compilers
DEPFLAGS = -MT $@ -MD -MP -MF $(DEPDIR)/$*.Td

# OpenMP
FFLAGS += -fopenmp
CPPFLAGS += -fopenmp

# HDF5
HDF5_INCLUDE_DIR=/usr/include/hdf5/serial
HDF5_LIBRARY_DIR=/usr/lib/x86_64-linux-gnu/hdf5
HDF5_LIBRARY=hdf5

# LibXC
LIBXC_INCLUDE_DIR=/home/fmoitzi/Libraries/libxc-5.1.5/include
LIBXC_LIBRARY_DIR=/home/fmoitzi/Libraries/libxc-5.1.5/bin
LIBXC_LIBRARY=xc

# Add definitions
CPPFLAGS += $(addprefix -D, USE_LIBXC)

# Lua
LUA_INCLUDE_DIR=/home/fmoitzi/Libraries/lua/include
LUA_LIBRARY_DIR=/home/fmoitzi/Libraries/lua/lib
LUA_LIBRARY=lua

# Lapack/Blas
LAPACK_INCLUDE_DIR=/usr/include/mkl
LAPACK_LIBRARY_DIR=/usr/lib/x86_64-linux-gnu/mkl
LAPACK_LIBRARY=lapack

# Add all include directories to the CPPFLAGS
CPPFLAGS += $(addprefix -I,$(HDF5_INCLUDE_DIR))
CPPFLAGS += $(addprefix -I,$(LIBXC_INCLUDE_DIR))
CPPFLAGS += $(addprefix -I,$(LUA_INCLUDE_DIR))

# Linker flags
LDLIBS += $(addprefix -L,$(HDF5_LIBRARY_DIR))
LDLIBS += $(addprefix -L,$(LIBXC_LIBRARY_DIR))
LDLIBS += $(addprefix -L,$(LUA_INCLUDE_DIR))

# Libaries
LDLIBS += $(addprefix -l,$(HDF5_LIBRARY))
LDLIBS += $(addprefix -l,$(LIBXC_LIBRARY))
LDLIBS += $(addprefix -l,$(LUA_LIBRARY))
