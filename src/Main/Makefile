all: lsms wl-lsms rewl-lsms

ifdef USE_GPTL
  gptl_def=-DUSE_GPTL
else
  gptl_def=
endif

LINKER ?= $(CXX)

LSMS_OBJ = SystemParameters.o read_input.o PotentialIO.o buildLIZandCommLists.o \
           energyContourIntegration.o solveSingleScatterers.o \
           calculateDensities.o calculateChemPot.o checkConsistency.o \
           lsmsClass.o calculateEvec.o initializeAtom.o mixing.o \
           ReplicaExchangeWL.o AlloyBankIO.o rotateToGlobal.o \
           write_restart.o

clean:
	rm -f *.o *.a lsms $(TOP_DIR)/bin/lsms \
                      lsms_inter $(TOP_DIR)/bin/lsms_inter \
                      wl-lsms $(TOP_DIR)/bin/wl-lsms \
                      rewl-lsms $(TOP_DIR)/bin/rewl-lsms\
		      zblock_lu_driver

zblock_lu_driver: zblock_lu_cuda_driver.o libLSMS.a
	$(LINKER) -fopenmp -UBUILDKKRMATRIX_GPU $(INC_PATH) -DNO_HDF5 -o zblock_lu_driver zblock_lu_cuda_driver.o libLSMS.a $(LIB_DIR) -L../Accelerator -lAccelerator $(ADD_LIBS)

lsms: lsms.o libLSMS.a
	$(LINKER) $(INC_PATH) $(gptl_def) -o lsms lsms.o -L. -lLSMS $(LIBS) $(ADD_LIBS)
	cp lsms $(TOP_DIR)/bin

lsms_inter: lsms_inter.o
	$(LINKER) $(INC_PATH) -o lsms_inter lsms_inter.o $(LIBS)
	cp lsms_inter $(TOP_DIR)/bin

wl-lsms: wl_lsms.o libLSMS.a
	$(LINKER) $(INC_PATH) -o wl-lsms wl_lsms.o -L. -lLSMS $(LIBS) $(ADD_LIBS)
	cp wl-lsms $(TOP_DIR)/bin

rewl-lsms: rewl_lsms.o libLSMS.a
	$(LINKER) $(INC_PATH) -o rewl-lsms rewl_lsms.o -L. -lLSMS $(LIBS) $(ADD_LIBS)
	cp rewl-lsms $(TOP_DIR)/bin

%.o : %.cpp
	$(CXX) $(INC_PATH) -c -o $@ $<

libLSMS.a: $(LSMS_OBJ)
	$(AR) -r libLSMS.a $(LSMS_OBJ) 

