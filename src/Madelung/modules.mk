
local_dir := src/Madelung

local_src := $(addprefix $(local_dir)/,\
		 cal_madelung_matrix.f \
         getkncut.f \
         getrscut.f \
         getstruc.f \
         interf.f \
         interfsmr.f \
         lattice.f \
         madewd.f \
         madewdj.f \
         madsum.f \
         ord3v.f \
         pqintg_c.f \
         lmfacts.f \
         bessj.f \
         Madelung.cpp)

sources += $(local_src)