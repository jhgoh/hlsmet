## Makefile to build and run the met testbed programs
#APP = vivado_hls
APP = vitis_hls
SRCS = test_met.cpp
SRCS += $(wildcard src/*.cpp) $(wildcard src/*.h)
SRCS += $(wildcard bench/*.cpp) $(wildcard bench.*.h)
PRJPATH = proj0/solution1

all: synth

CSIMEXE := $(PRJPATH)/csim/build/csim.exe
$(CSIMEXE): tcl/do_met.tcl tcl/do_corrmet.tcl $(SRCS)
	@echo $(CSIMEXE) $(wildcard $(CSIMEXE))
ifeq ($(wildcard $(CSIMEXE)),)
	@echo "Cannot find $(CSIMEXE). Please build the csim first then continue"
	@echo "Available targets:"
	@echo "- mettest"
	@echo "- corrmettest"
	exit 1
endif

mettest: tcl/do_met.tcl $(SRCS)
	$(APP) -f tcl/do_met.tcl

corrmettest: tcl/do_corrmet.tcl $(SRCS)
	$(APP) -f tcl/do_corrmet.tcl

synth: tcl/do_synth.tcl | $(CSIMEXE)
	@echo do_synth
	$(APP) -f tcl/do_synth.tcl

clean:
	rm -rf proj*
