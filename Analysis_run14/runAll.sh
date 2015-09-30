#!/bin/bash
date
# . ./runAll.sh Batch_List_39GeV_d100 0k20k_SE 0
#               name of the list
#                                     suffix
#                                           same event or mixed event flag
#ifile=$1

    # Analysis_num:
    # 1   = phi
    # 2   = Lambda
    # 3   = K0s
    # 4   = Xi
    # 5   = Omega
    # 6   = D0 + KStar
    # 7   = KStar +/-
    # 8   = hadron v2
    # 88  = hadron spectra
    # 9   = vertex + event display
    # 10  = phi correction
    # 11  = event plane analysis
    # 111 = event plane analysis with different harmonics + eta gap values
    # 112 = BBC analysis
    # 12  = phi same event K+K+ and K-K-
    # 13  = phi V0
    # 14  = anti-Lambda
    # 15  = rho
    # 16  = anti-Xi (Xi+)
    # 17  = Xi(1530) -> Xi- + pi+
    # 18  = Xi(1530) -> Xi+ + pi-
    # 19  = anti-Omega (Omega+)
    # 20  = Merged Omega + Xi analysis + anti-particles
    # 200 = Omega(2250) analysis
    # 201 = d4s analysis
    # 21  = Sigma(1385) analysis
    # 22  = RefMult QA analysis
    # 23  = K->pi pi pi analysis
    # 24  = Theta+ -> K0S + p analysis
    # 25  = J/Psi -> e- + e+
    # 26  = Patrick's e+ e- analysis
    # 30  = proton-proton correlation analysis
    # 124 = Patrick's single tree analysis
    # 125 = Event analysis
    # 130 = m2 vs nSigmaP analysis
    # 300 = Jet analysis
    # 301 = Jet analysis - with randomized tracks
    # 302 = Fill Jet tree
    # 35  = LambdaC+ analysis
    # 351 = LambdaC+ V0 analysis
    # 360 = Ach, mean pt analysis
    # 400 = D0 HFT analysis
    
    analysis=400;

    # BeamTime_num:
    # 0 = 7.7
    # 1 = 11.5
    # 2 = 39
    # 3 = 62.4
    # 4 = 19.6
    # 5 = 27
    # 6 = 200
    # 7 = 14.5
    # 8 = 200 GeV run 14
    
    beamtime=8;
    
    # 0 = same event
    # 1 = mixed event

suffix="D0_"$2
if [ $# -eq 3 ] 
then                
for FILE in `cat $1`
do
    echo $FILE
    cp ./run.csh ./run_dev2$2_$suffix$FILE.csh
    
    echo -n "root4star -b -q -x 'Analysis_Snurf_Macro.cc(">>run_dev2$2_$suffix$FILE.csh
    echo -n '200001,250000,'>>run_dev2$2_$suffix$FILE.csh
    echo -n $3>>run_dev2$2_$suffix$FILE.csh 
    echo -n ',"'>>run_dev2$2_$suffix$FILE.csh
    echo -n $FILE>>run_dev2$2_$suffix$FILE.csh
    echo -n '","' >>run_dev2$2_$suffix$FILE.csh
    echo -n $suffix>>run_dev2$2_$suffix$FILE.csh
    echo -n '",' >>run_dev2$2_$suffix$FILE.csh 
    echo -n $analysis>>run_dev2$2_$suffix$FILE.csh
    echo -n ',' >>run_dev2$2_$suffix$FILE.csh
    echo -n $beamtime>>run_dev2$2_$suffix$FILE.csh
    echo -n ')' >>run_dev2$2_$suffix$FILE.csh
    echo -n "' > log/run_dev2">>run_dev2$2_$suffix$FILE.csh
    echo -n $2>>run_dev2$2_$suffix$FILE.csh
    echo -n "_">>run_dev2$2_$suffix$FILE.csh
    echo -n $suffix>>run_dev2$2_$suffix$FILE.csh
    echo -n "_">>run_dev2$2_$suffix$FILE.csh
    echo -n $FILE>>run_dev2$2_$suffix$FILE.csh
    echo ".log">>run_dev2$2_$suffix$FILE.csh

    qsub  -hard -l projectio=1,scratchfree=500,h_cpu=24:00:00,h_vmem=1.5G -o log/job_dev2$2_$FILE.log -e log/job_dev2$2_$suffix_$FILE.err ./run_dev2$2_$suffix$FILE.csh
   
    mv run_dev2$2_$suffix$FILE.csh script/

#    let "ifile+=1";
done
else
	echo "**************************************************************************************"
	echo "one file list as parameter and a suffix"
        echo "**************************************************************************************" 
fi