
export TOP_DIR = $(shell pwd)
export INC_PATH =
export LIBS := -L$(TOP_DIR)/lua/lib -llua $(TOP_DIR)/mjson/mjson.a

include architecture.h 

ifdef USE_LIBXC
  ADDITIONAL_TARGETS += libxc
  ADD_LIBS += -L$(TOP_DIR)/opt/lib/ -lxc
  INC_PATH += -I$(TOP_DIR)/opt/include/
  export OPT_DEFINES += -DUSE_LIBXC
endif

ifdef HAS_BACKTRACE
  export OPT_DEFINES += -DHAS_BACKTRACE
endif

all: liblua $(ADDITIONAL_TARGETS) libmjson LSMS
# all: liblua libjson $(ADDITIONAL_TARGETS) libmjson LSMS
# all: liblua LSMS Documentation

.PHONY: libxc clean LSMS Documentation liblua libjson libmjson Tools

clean:
	cd lua && $(MAKE) clean
	cd mjson && $(MAKE) clean
	cd src && $(MAKE) clean
	cd lib && $(MAKE) clean
	cd doc && $(MAKE) clean
	cd Tools && $(MAKE) clean
#	cd libjson && $(MAKE) clean

LSMS: liblua libmjson $(ADDITIONAL_TARGETS)
	cd src && $(MAKE)

Tools: liblua 
	cd Tools && $(MAKE)

Documentation:
	cd doc && $(MAKE)

liblua:
	cd lua; $(MAKE); $(MAKE) local

libmjson:
	cd mjson && $(MAKE)

# libjson:
#	cd libjson && $(MAKE)

zblock_lu_driver: liblua 
	cd src && $(MAKE) zblock_lu_driver


test: liblua $(ADDITIONAL_TARGETS) libmjson
	cd src && $(MAKE) test
