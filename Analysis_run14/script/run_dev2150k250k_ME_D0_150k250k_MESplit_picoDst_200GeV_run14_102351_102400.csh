#!/bin/csh

#starver SL11d
#starver dev
starver new
cd /global/homes/a/aschmah/STAR/Analysis_run14/

root4star -b -q -x 'Analysis_Snurf_Macro.cc(150001,250000,1,"Split_picoDst_200GeV_run14_102351_102400","D0_150k250k_ME",400,8)' > log/run_dev2150k250k_ME_D0_150k250k_ME_Split_picoDst_200GeV_run14_102351_102400.log
