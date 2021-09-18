
local_dir := src/Misc

local_src := $(addprefix $(local_dir)/,\
		associatedLegendreFunction.cpp \
        bulirsch_stoer.f \
        bulirschStoerIntegrator.cpp \
        calculateGauntCoeficients.cpp \
        cgaunt_c.f \
        clock_time.c \
        cinterp.f \
        clebsch.f \
        cmtruni.f \
        cnewint.f \
        Coeficients.cpp \
        congauss_c.f \
        constraint.f \
        dfv.f \
        dfv_new.f \
        essl_workaround.f \
        fit.f \
        fitpot.f \
        fnpi.f \
        fstop.f \
        gaunt.f \
        getclm.f \
        ifacts_c.f \
        initwave.f \
        interp.f \
        matr.f \
        matrot1.f \
        mbeqa.f \
        mod_midpoint.f \
        newder.f \
        newint.f \
        plglmax.f \
        quadrature.cpp \
        readLastLine.cpp \
        ricbes.f \
        rotmat.f \
        rsimp.f \
        rwave.f \
        rzextr.f \
        spin_trafo.f \
        trltog.f \
        u_sigma_u.f \
        v_plus_minus.f \
        wrtmtx.f \
        ylag.f \
        zeroout.f \
        zsphbes.f \
        zsphbesj.f \
        zsphbesjh.f)

sources += $(local_src)