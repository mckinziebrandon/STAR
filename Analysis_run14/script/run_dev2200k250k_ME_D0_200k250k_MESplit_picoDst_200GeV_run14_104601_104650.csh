#!/bin/csh

#starver SL11d
#starver dev
starver new
cd /global/homes/a/aschmah/STAR/Analysis_run14/

root4star -b -q -x 'Analysis_Snurf_Macro.cc(200001,250000,1,"Split_picoDst_200GeV_run14_104601_104650","D0_200k250k_ME",400,8)' > log/run_dev2200k250k_ME_D0_200k250k_ME_Split_picoDst_200GeV_run14_104601_104650.log
