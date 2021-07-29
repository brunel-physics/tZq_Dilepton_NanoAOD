#!/bin/bash


SystematicNamesArray=("_Nominal"                       "_PU_ScaleUp"             "_PU_ScaleDown"           "_BTag_ScaleUp"              "_BTag_ScaleDown" 
		      "_JetSmearing_ScaleUp"           "_JetSmearing_ScaleDown"  "_JetResolution_ScaleUp"  "_JetResolution_ScaleDown"   "_LeptonEfficiencies_ScaleUp"
		      "_LeptonEfficiencies_ScaleDown"  "_PDF_ScaleUp"  		 "_PDF_ScaleDown"	   "_ME_Up"			"_ME_Down"
		      "_MET_Up"			       "_MET_Down"		 "_isr_up"		   "_isr_down"			"_fsr_up"
		      "_fsr_down")
	


ChannelArray=("_ee" "_mumu" "_emu")

YearArray=("2016" "2017" "2018")

for i in ${!SystematicNamesArray[@]}; do

	tmux_string="Data${SystematicNamesArray[i]}${ChannelArray[1]}_${YearArray[2]}"

	echo './bin/fulleventselectionMain.exe --mc 0 -y 2018 -p 159 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys '$i' --channel 2 --dcc 0'

	tmux new -d -s $tmux_string 'source /cvmfs/sft.cern.ch/lcg/views/LCG_96/x86_64-slc6-gcc8-opt/setup.sh; make clean; make; ./bin/fulleventselectionMain.exe --mc 0 -y 2018 -p 160 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys '$i' --channel 2 --dcc 0; sleep 86400'

	sleep 180

done



echo "The tmux sessions for the side band region (mumu channel for 2018 data) are running for all systematics."


