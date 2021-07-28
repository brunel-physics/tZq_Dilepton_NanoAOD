#!/bin/bash

ProcessNamesArray=("tZq_scaleup" "tZq_scaledown") 

ChannelArray=("_ee" "_mumu")

YearArray=("2016")

for i in ${!ChannelArray[@]}; do

	tmux_string_ScaleUp="${ProcessNamesArray[0]}${ChannelArray[i]}_${YearArray[0]}"

	echo './bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 2 --npl 0 --sr 1 --sbr 0 --zjcr 0 --ttcr 0 --sys 0 --channel '$i' --dcc 0'

	tmux new -d -s $tmux_string_ScaleUp 'source /cvmfs/sft.cern.ch/lcg/views/LCG_96/x86_64-slc6-gcc8-opt/setup.sh; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 2 --npl 0 --sr 1 --sbr 0 --zjcr 0 --ttcr 0 --sys 0 --channel '$i' --dcc 0; sleep 86400'

	sleep 180

	tmux_string_ScaleDown="${ProcessNamesArray[1]}${ChannelArray[i]}_${YearArray[0]}"

        echo './bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 3 --npl 0 --sr 1 --sbr 0 --zjcr 0 --ttcr 0 --sys 0 --channel '$i' --dcc 0'

        tmux new -d -s $tmux_string_ScaleDown 'source /cvmfs/sft.cern.ch/lcg/views/LCG_96/x86_64-slc6-gcc8-opt/setup.sh; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 3 --npl 0 --sr 1 --sbr 0 --zjcr 0 --ttcr 0 --sys 0 --channel '$i' --dcc 0; sleep 86400'

        sleep 180

done



echo "The tmux sessions for tZq scale up and scale down in the SR (ee and mumu channels for 2016) are running."


