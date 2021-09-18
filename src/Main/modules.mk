
local_dir := src/Main

local_src := $(addprefix $(local_dir)/,\
		lsms.cpp \
        SystemParameters.cpp \
        read_input.cpp \
        PotentialIO.cpp \
        buildLIZandCommLists.cpp \
        energyContourIntegration.cpp \
        solveSingleScatterers.cpp \
        calculateDensities.cpp \
        calculateChemPot.cpp \
        checkConsistency.cpp \
        calculateEvec.cpp \
        initializeAtom.cpp \
        mixing.cpp \
        ReplicaExchangeWL.cpp \
        AlloyBankIO.cpp \
        rotateToGlobal.cpp \
        write_restart.cpp)

sources += $(local_src)