#!/bin/csh

#starver SL11d
#starver dev
starver new
cd /global/homes/a/aschmah/STAR/Analysis_run14/

root4star -b -q -x 'Analysis_Snurf_Macro.cc(250001,2500000,1,"Split_picoDst_200GeV_run14_100901_100950","D0_250kFull_ME",400,8)' > log/run_dev2250kFull_ME_D0_250kFull_ME_Split_picoDst_200GeV_run14_100901_100950.log
