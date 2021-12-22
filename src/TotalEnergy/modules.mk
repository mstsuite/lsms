
local_dir := src/TotalEnergy

local_src := $(addprefix $(local_dir)/,\
        	calculateTotalEnergy.cpp \
        	localTotalEnergy.cpp \
        	zeropt_c.f \
        	janake_c.f)

sources += $(local_src)