
local_dir := src/Core

local_src := $(addprefix $(local_dir)/,\
        calculateCoreStates.cpp \
        coreSolver.cpp \
        corslv_c.f \
        deepst_c.f \
        getcor_c.f \
        invals_c.f \
        inwhnk_c.f \
        inws_c.f \
        outws_c.f \
        richnk_c.f \
        semcst_c.f)

sources += $(local_src)