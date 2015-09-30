#!/bin/csh

#starver SL11d
starver pro
#which root4star
cd /global/homes/a/aschmah/STAR/Analysis/Jet/Simulation/

 root4star -b -q 'JetSimAnalysis_Macro.cc(0,300000,0.3,"JetSim_list2_2")' >! log/run_devJetSim_file_13.log
