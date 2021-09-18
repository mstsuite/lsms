
local_dir := src/LuaInterface

local_src := $(addprefix $(local_dir)/,\
		 LuaInterface.cpp \
         LuaSupport.cpp)

sources += $(local_src)