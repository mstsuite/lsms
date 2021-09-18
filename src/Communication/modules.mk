
local_dir := src/Communication

local_src := $(addprefix $(local_dir)/,\
        			 distributeAtoms.cpp \
         			 LSMSCommunication.cpp \
                 	 REWLCommunication.cpp)

sources += $(local_src)