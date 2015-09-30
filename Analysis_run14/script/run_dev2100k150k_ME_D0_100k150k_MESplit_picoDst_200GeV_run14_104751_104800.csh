#!/bin/csh

#starver SL11d
#starver dev
starver new
cd /global/homes/a/aschmah/STAR/Analysis_run14/

root4star -b -q -x 'Analysis_Snurf_Macro.cc(100001,150000,1,"Split_picoDst_200GeV_run14_104751_104800","D0_100k150k_ME",400,8)' > log/run_dev2100k150k_ME_D0_100k150k_ME_Split_picoDst_200GeV_run14_104751_104800.log
