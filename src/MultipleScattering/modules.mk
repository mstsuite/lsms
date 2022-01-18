
local_dir := src/MultipleScattering

local_src := $(addprefix $(local_dir)/,\
        calculateTauMatrix.cpp \
        buildKKRMatrix_CPU.cpp \
        linearSolvers_CPU.cpp \
        makegij_c.f \
        setgij.f \
        block_inverse_fortran.f \
        zblock_lu.F \
        wasinv.f \
        zmar1.f \
        wasinv_p.f \
        zuqmx.f \
        zutfx.f \
        zucpx.f \
        zaxpby.f \
        zrandn.f \
        tau_inv_postproc.f \
        trgtol.f \
        green_function.f90 \
        gf_local.f90 \
        int_zz_zj.f90 \
        mdosms_c.f \
        mgreen_c.f \
        green_function_rel.f \
        write_kkrmat.f \
        relmtrx.f \
        gfill.f \
        gafill.f \
        magnet.f \
        magnetic_dens.f \
        new_dens.f \
        block_inverse.cpp \
        zblock_lu_cpp.cpp \
        zblock_lu_cublas.cpp)

sources += $(local_src)