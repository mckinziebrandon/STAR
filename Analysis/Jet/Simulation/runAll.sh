#!/bin/bash
date

# qdel $(qselect -u aschmah)
#. ./runAll.sh Batch_JetSim_list
    
if [ $# -eq 1 ] 
then
  
  StartEvent=0;
  NEvents=300000;
  Radius=0.3;
  
  counter_file=0;
  for FILE in `cat $1`
  do
     echo "PYTHIA analysis started, file counter: " $counter_file 
     
      suffix_run="JetSim_file_"$counter_file;
      cp ./run.csh ./run_dev$suffix_run.csh
      echo -n " " >>./run_dev$suffix_run.csh
      echo -n "root4star -b -q 'JetSimAnalysis_Macro.cc(">>run_dev$suffix_run.csh
      echo -n $StartEvent","$NEvents","$Radius","'"'$FILE'"' >>run_dev$suffix_run.csh
      echo -n ')' >>run_dev$suffix_run.csh
      echo -n "' >! log/run_dev$suffix_run">>run_dev$suffix_run.csh
      #echo -n $suffix_run>>run_dev$suffix_run.csh
      echo ".log">>run_dev$suffix_run.csh
      
      qsub  -hard -l projectio=1,scratchfree=500,h_cpu=01:00:00,h_vmem=2.0G -o log/job_dev$suffix_run.log -e log/job_dev$suffix_run.err ./run_dev$suffix_run.csh
  
      mv run_dev$suffix_run.csh script/
      let "counter_file=counter_file+1"
  
  done

else
echo "**************************************************************************************"
echo "wrong number of parameters"
echo "**************************************************************************************" 
fi
