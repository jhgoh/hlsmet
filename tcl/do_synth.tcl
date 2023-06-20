# open the project
open_project proj0
open_solution "solution1"

# do stuff
csynth_design
cosim_design -trace_level all
export_design -format ip_catalog  -vendor "cern-cms"

# exit Vivado HLS
exit
