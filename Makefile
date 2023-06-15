## Makefile to build and run the met testbed programs
#APP = vivado_hls
APP = vitis_hls

mettest:
	$(APP) -f do_met.tcl

all: mettest

clean:
	rm -rf proj*
