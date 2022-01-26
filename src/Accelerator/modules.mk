
local_dir := src/Accelerator

local_src := $(addprefix $(local_dir)/,\
        accelerator_initialize.F \
        accelerator_finalize.F \
        Accelerator.cpp \
        fortran.c)

sources += $(local_src)



