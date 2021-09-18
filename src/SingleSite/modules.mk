
local_dir := src/SingleSite

local_src := $(addprefix $(local_dir)/,\
        			 SingleSiteScattering.cpp \
                     readSingleAtomData_hdf5.cpp \
                     readSingleAtomData_bigcell.cpp \
                     F_readSingleAtomData_bigcell.f90 \
                     single_scatterer_nonrel.f \
                     semrel_c.f \
                     scalar_m.f \
                     single_scatterer_rel.f \
                     spzwafu.f \
                     csbf.f	\
                     matops.f \
                     gjinv.f \
                     dirmag1-op.f \
                     dirmag2-op.f \
                     brmat.f \
                     writeSingleAtomData_hdf5.cpp \
                     writeSingleAtomData_bigcell.cpp \
                     F_writeSingleAtomData_bigcell.f90 \
                     checkAntiFerromagneticStatus.cpp)

sources += $(local_src)