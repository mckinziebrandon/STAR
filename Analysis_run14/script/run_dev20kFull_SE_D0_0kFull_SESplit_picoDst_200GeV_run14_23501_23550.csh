#!/bin/csh

#starver SL11d
#starver dev
starver new
cd /global/homes/a/aschmah/STAR/Analysis_run14/

root4star -b -q -x 'Analysis_Snurf_Macro.cc(0,9000000,0,"Split_picoDst_200GeV_run14_23501_23550","D0_0kFull_SE",400,8)' > log/run_dev20kFull_SE_D0_0kFull_SE_Split_picoDst_200GeV_run14_23501_23550.log
