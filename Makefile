## Makefile to build and run the met testbed programs
#APP = vivado_hls
APP = vitis_hls
SRCS = test_met.cpp
SRCS += $(wildcard src/*.cpp) $(wildcard src/*.h)
SRCS += $(wildcard bench/*.cpp) $(wildcard bench.*.h)
PRJPATH = proj0/solution1

all: design

CSIMEXE := $(PRJPATH)/csim/build/csim.exe
csim: | $(CSIMEXE)

$(CSIMEXE): $(SRCS)
	$(APP) -f tcl/do_met.tcl

SYNTHOUT := $(PRJPATH)/syn
synth: $(CSIMEXE) | $(SYNTHOUT)

$(SYNTHOUT):
	$(APP) -f tcl/do_csynth.tcl

COSIMOUT := $(PRJPATH)/sim
cosim: $(SYNTHOUT) | $(COSIMOUT)

$(COSIMOUT):
	$(APP) -f tcl/do_cosim.tcl

design: $(COSIMOUT)
	$(APP) -f tcl/do_design.tcl

clean:
	rm -rf proj*
