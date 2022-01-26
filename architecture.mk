
# C compiler
CC := gcc
# C++ compiler
CXX := mpic++
# Fortran compiler
F90 := gfortran
# linker
LD := mpic++
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
# linker flags: libraries to link
LDLIBS :=
# flags required for dependency generation; passed to compilers
DEPFLAGS = -MT $@ -MD -MP -MF $(DEPDIR)/$*.Td

USE_OPENMP=1

ifdef USE_OPENMP
	# OpenMP
	FFLAGS += -fopenmp
	CPPFLAGS += -fopenmp
endif

# HDF5
HDF5_INCLUDE_DIR=/usr/include/hdf5/serial
HDF5_LIBRARY_DIR=/usr/lib/x86_64-linux-gnu/hdf5/serial
HDF5_LIBRARY=hdf5

# LibXC
LIBXC_INCLUDE_DIR=/home/fmoitzi/Libraries/libxc-5.1.5/include
LIBXC_LIBRARY_DIR=/home/fmoitzi/Libraries/libxc-5.1.5/lib
LIBXC_LIBRARY=xc

# Add definitions
CPPFLAGS += $(addprefix -D, USE_LIBXC)

# Lua
LUA_INCLUDE_DIR=/home/fmoitzi/Libraries/lua/include
LUA_LIBRARY_DIR=/home/fmoitzi/Libraries/lua/lib
LUA_LIBRARY=lua

# Lapack/Blas
MKL_INCLUDE_DIR=/usr/include/mkl
MKL_LIBRARY_DIR=/usr/lib/x86_64-linux-gnu
MKL_LIBRARY=mkl_rt

# Add all include directories to the CPPFLAGS
CPPFLAGS += $(addprefix -I,$(HDF5_INCLUDE_DIR))
CPPFLAGS += $(addprefix -I,$(LIBXC_INCLUDE_DIR))
CPPFLAGS += $(addprefix -I,$(LUA_INCLUDE_DIR))
CPPFLAGS += $(addprefix -I,$(MKL_INCLUDE_DIR))

# Linker flags
LDLIBS += $(addprefix -L,/usr/lib/x86_64-linux-gnu)

LDLIBS += $(addprefix -L,$(HDF5_LIBRARY_DIR))
LDLIBS += $(addprefix -L,$(LIBXC_LIBRARY_DIR))
LDLIBS += $(addprefix -L,$(LUA_LIBRARY_DIR))
LDLIBS += $(addprefix -L,$(MKL_LIBRARY_DIR))

# Libaries
LDLIBS += $(addprefix -l,$(HDF5_LIBRARY))
LDLIBS += $(addprefix -l,$(LIBXC_LIBRARY))
LDLIBS += $(addprefix -l,$(LUA_LIBRARY))
LDLIBS += $(addprefix -l,$(MKL_LIBRARY))

# Gfortran
LDLIBS += $(addprefix -l,gfortran)

#
LDLIBS += $(addprefix -l,dl)

ifdef USE_OPENMP
	LDLIBS += $(addprefix -l,gomp)
endif