#!/bin/csh

#starver SL11d
#starver dev
starver new
cd /global/homes/a/aschmah/STAR/Analysis_run14/

root4star -b -q -x 'Analysis_Snurf_Macro.cc(0,100000,1,"Split_picoDst_200GeV_run14_105951_106000","D0_0k100k_ME",400,8)' > log/run_dev20k100k_ME_D0_0k100k_ME_Split_picoDst_200GeV_run14_105951_106000.log
