#!/bin/csh

#starver dev
starver new
cd /global/homes/b/bsmckinz/STAR/Analysis/HFT/D0_Analysis/

root4star -b -q -x 'D0_Analysis_macro.cc(0,0,90000000,"Split_D0tree_200GeV_run14_preview2_31_40")' > log/run_dev2Split_D0tree_200GeV_run14_preview2_31_40.log
