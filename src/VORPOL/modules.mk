
local_dir := src/VORPOL

local_src := $(addprefix $(local_dir)/,\
					 calsig.f \
            		 caltnode.f\
                     celbnd.f \
                     chkbnd.f \
                     chkedge.f \
                     filter_edge.f \
                     inter_dip.f \
                     inter.f \
                     interstitial.f \
                     intphi.f \
                     intpl0.f \
                     invm3.f \
                     polyhedron.f \
                     rcritpts.f90 \
                     setup_boundary.f \
                     setup_vorpol.f \
                     sigma.f \
                     sort.f \
                     sortidx.f \
                     stepyll.f \
                     volvor.f \
                     setupVorpol.cpp)

sources += $(local_src)