#!/bin/bash

ProcessNamesArray_mumu_2016_NominalRuns=("tZq_mumu_2016_nominal" "tZq_scaleup_mumu_2016_nominal" "tZq_scaledown_mumu_2016_nominal")

for i in ${!ProcessNamesArray_mumu_2016_NominalRuns[@]}; do
	if [ ${ProcessNamesArray_mumu_2016_NominalRuns[$i]} = "tZq_mumu_2016_nominal" ]
		then tmux new -d -s tZq_mumu_2016_nominal 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 0 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 30'
	elif [ ${ProcessNamesArray_mumu_2016_NominalRuns[$i]} = "tZq_scaleup_mumu_2016_nominal" ]
                then tmux new -d -s tZq_scaleup_mumu_2016_nominal 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 1 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 30'
	elif [ ${ProcessNamesArray_mumu_2016_NominalRuns[$i]} = "tZq_scaledown_mumu_2016_nominal" ]
                then tmux new -d -s tZq_scaledown_mumu_2016_nominal 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 2 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 30'
	else 
		echo "ERROR: Check ProcessNamesArray"
		break 
	fi 
done

echo "The tmux sessions for the SBR (mumu channel 2016) are running"




