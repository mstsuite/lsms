
local_dir := src/Potential

local_src := $(addprefix $(local_dir)/,\
				alpha2_c.f \
                calculateASAro3.f \
                calculateChargesPotential.cpp \
                epcorr.f \
                getqm_mt.f \
                getvmt.cpp \
                getXCName.cpp \
                interpolatePotential.cpp \
                libxcInterface.cpp \
                newexchg.f \
                newpot_c.f \
                PotentialShifter.cpp \
                newFunctionalInterface.cpp \
                lsf_functional.cpp \
                rs.f)

sources += $(local_src)

