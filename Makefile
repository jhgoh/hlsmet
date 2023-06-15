APP = vivado_hls

mettest:
	$(APP) -f do_met.tcl

all: mettest

clean:
	rm -rf proj*
