
local_dir := src/MultipoleMadelung

local_src := $(addprefix $(local_dir)/,\
	clm.f90 \
	gaunt_factor.cpp \
	gaunt_factor.f90 \
	gaussq.f90 \
	integer_factors.cpp \
	integer_factors.f90 \
	lattice_utils.cpp \
	legendre.f90 \
	madelung.cpp \
	MultipoleMadelung.cpp \
	spherical_harmonics.f90)

sources += $(local_src)

