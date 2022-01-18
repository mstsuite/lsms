
# Defines common flags for compiling the program
include architecture.mk

# Collect information from each module in these four variables.
# Initialize them here as simple variables.
#programs :=
sources :=
#libraries :=
#extra_clean :=

# name of output binary
BINNAME := lsms
BINDIR := $(realpath .)/bin
BIN := $(BINDIR)/$(BINNAME)

include_dirs := src include
CPPFLAGS += $(addprefix -I ,$(include_dirs))
vpath %.h $(include_dirs)

include src/Accelerator/modules.mk
include src/Communication/modules.mk
include src/Core/modules.mk
include src/LuaInterface/modules.mk
include src/Madelung/modules.mk
include src/Main/modules.mk
include src/Misc/modules.mk
include src/MultipleScattering/modules.mk
include src/MultipoleMadelung/modules.mk
include src/Potential/modules.mk
include src/RadialGrid/modules.mk
include src/SingleSite/modules.mk
include src/TotalEnergy/modules.mk
include src/VORPOL/modules.mk

SRCS = $(sources)

# intermediate directory for generated object files
OBJDIR := .o
# intermediate directory for generated dependency files
DEPDIR := .d

# object files, auto generated from source files
OBJS := $(patsubst %,$(OBJDIR)/%.o,$(basename $(SRCS)))
# dependency files, auto generated from source files
DEPS := $(patsubst %,$(DEPDIR)/%.d,$(basename $(SRCS)))

# compilers (at least gcc and clang) don't create the subdirectories automatically
$(shell mkdir -p $(dir $(OBJS)) >/dev/null)
$(shell mkdir -p $(dir $(DEPS)) >/dev/null)
$(shell mkdir -p $(dir $(BIN)) >/dev/null)

# compile C source files
COMPILE.c = $(CC) $(DEPFLAGS) $(CFLAGS) $(CPPFLAGS) -c -o $@
# compile C++ source files
COMPILE.cc = $(CXX) $(DEPFLAGS) $(CXXFLAGS) $(CPPFLAGS) -c -o $@
# compile Fortran source files
COMPILE.f90 = $(F90) $(DEPFLAGS) $(FFLAGS) $(FPPFLAGS) -c -o $@
# link object files to binary
LINK.o = $(LD) $(LDFLAGS) -o $@
# precompile step
PRECOMPILE =
# postcompile step
POSTCOMPILE = mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d

all: $(BIN)

.PHONY: clean
clean:
	$(RM) -r $(OBJDIR) $(DEPDIR) $(BIN) $(BINDIR)

.PHONY: install
install:
	@echo no install tasks configured

.PHONY: convert
convert:
	@echo
	find . -type f -print0 | xargs -0 dos2unix

$(BIN): $(OBJS)
	$(LINK.o) $^ $(LDLIBS)

$(OBJDIR)/%.o: %.f
$(OBJDIR)/%.o: %.f $(DEPDIR)/%.d
	$(PRECOMPILE)
	$(COMPILE.f90) $<
	$(POSTCOMPILE)

$(OBJDIR)/%.o: %.F
$(OBJDIR)/%.o: %.F $(DEPDIR)/%.d
	$(PRECOMPILE)
	$(COMPILE.f90) $<
	$(POSTCOMPILE)

$(OBJDIR)/%.o: %.f90
$(OBJDIR)/%.o: %.f90 $(DEPDIR)/%.d
	$(PRECOMPILE)
	$(COMPILE.f90) $<
	$(POSTCOMPILE)

$(OBJDIR)/%.o: %.c
$(OBJDIR)/%.o: %.c $(DEPDIR)/%.d
	$(PRECOMPILE)
	$(COMPILE.c) $<
	$(POSTCOMPILE)

$(OBJDIR)/%.o: %.cpp
$(OBJDIR)/%.o: %.cpp $(DEPDIR)/%.d
	$(PRECOMPILE)
	$(COMPILE.cc) $<
	$(POSTCOMPILE)

$(OBJDIR)/%.o: %.cc
$(OBJDIR)/%.o: %.cc $(DEPDIR)/%.d
	$(PRECOMPILE)
	$(COMPILE.cc) $<
	$(POSTCOMPILE)

$(OBJDIR)/%.o: %.cxx
$(OBJDIR)/%.o: %.cxx $(DEPDIR)/%.d
	$(PRECOMPILE)
	$(COMPILE.cc) $<
	$(POSTCOMPILE)

.PRECIOUS: $(DEPDIR)/%.d
$(DEPDIR)/%.d: ;

include $(DEPS)