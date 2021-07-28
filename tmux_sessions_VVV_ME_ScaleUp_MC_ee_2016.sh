#!/bin/bash

ProcessNamesArray=("VVV_WWWTo4F"			"VVV_WWWTo4F_ext"			"VVV_WWZTo4F"
		   "VVV_WWZTo4F_ext"			"VVV_WZZ"				"VVV_WZZ_ext"
		   "VVV_ZZZ"				"VVV_ZZZ_ext")				
 			        



SystematicNamesArray=("_Nominal"                       "_PU_ScaleUp"             "_PU_ScaleDown"           "_BTag_ScaleUp"              "_BTag_ScaleDown" 
		      "_JetSmearing_ScaleUp"           "_JetSmearing_ScaleDown"  "_JetResolution_ScaleUp"  "_JetResolution_ScaleDown"   "_LeptonEfficiencies_ScaleUp"
		      "_LeptonEfficiencies_ScaleDown"  "_PDF_ScaleUp"  		 "_PDF_ScaleDown"	   "_ME_Up"			"_ME_Down"
		      "_MET_Up"			       "_MET_Down"		 "_isr_up"		   "_isr_down"			"_fsr_up"
		      "_fsr_down")
	


ChannelArray=("_ee" "_mumu" "_emu")

YearArray=("2016" "2017" "2018")

for i in ${!ProcessNamesArray[@]}; do

	tmux_string="${ProcessNamesArray[i]}${SystematicNamesArray[13]}${ChannelArray[0]}_${YearArray[0]}"

	j=$(($i + 96))

	echo './bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$j' --npl 0 --sr 1 --sbr 0 --zjcr 0 --ttcr 0 --sys 13 --channel 1 --dcc 0'

	tmux new -d -s $tmux_string 'source /cvmfs/sft.cern.ch/lcg/views/LCG_96/x86_64-slc6-gcc8-opt/setup.sh; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$j' --npl 0 --sr 1 --sbr 0 --zjcr 0 --ttcr 0 --sys 13 --channel 1 --dcc 0; sleep 86400'

	sleep 180

done



echo "The tmux sessions for VVV in the SR (ee channel for 2016, ME scale up) are running."


