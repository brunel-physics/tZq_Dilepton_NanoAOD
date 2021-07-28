#!/bin/bash

ProcessNamesArray=("ttbar_madgraph_ext"                 "ttbar_TTToHadronic"                    "ttbar_TTToHadronic_ext"
                   "ttbar_TTToSemileptonic"             "ttbar_TTToSemileptonic_ext"            "ttbar_atMCaNLO"
                   "ttbar_atMCaNLO_ext"                 "ttbar_inc")

ProcessNamesArray2=("ttbar_hdampUP"			"ttbar_hdampUP_ext"			"ttbar_hdampDOWN"
	            "ttbar_hdampDOWN_ext"		"TT_2l2nu_hdampUP_ext1"			"TT_2l2nu_hdampUP_ext2"
		    "TT_2l2nu_hdampDOWN"	        "TT_2l2nu_hdampDOWN_ext1"		"TT_2l2nu_hdampDOWN_ext2"
		    "TTToHadronic_hdampUP"		"TTToHadronic_hdampDOWN"		"TTToSemileptonic_hdampUP"
		    "TTToSemileptonic_hdampDOWN"        "TTToSemileptonic_hdampDOWN_ext")
		   
ProcessNamesArray3=("ttbar_isr_UP"			"ttbar_isr_DOWN" 			"ttbar_isr_DOwN_ext"
		    "ttbar_fsr_UP"			"ttbar_fsr_UP_ext"			"ttbar_fsr_DOWN"
		    "ttbar_fsr_DOWN_ext")		        
 			        



SystematicNamesArray=("_Nominal"                       "_PU_ScaleUp"             "_PU_ScaleDown"           "_BTag_ScaleUp"              "_BTag_ScaleDown" 
		      "_JetSmearing_ScaleUp"           "_JetSmearing_ScaleDown"  "_JetResolution_ScaleUp"  "_JetResolution_ScaleDown"   "_LeptonEfficiencies_ScaleUp"
		      "_LeptonEfficiencies_ScaleDown"  "_PDF_ScaleUp"  		 "_PDF_ScaleDown"	   "_ME_Up"			"_ME_Down"
		      "_MET_Up"			       "_MET_Down"		 "_isr_up"		   "_isr_down"			"_fsr_up"
		      "_fsr_down")
	


ChannelArray=("_ee" "_mumu" "_emu")

YearArray=("2016" "2017" "2018")

for i in ${!ProcessNamesArray[@]}; do

	tmux_string="${ProcessNamesArray[i]}${SystematicNamesArray[1]}${ChannelArray[0]}_${YearArray[0]}"

	j=$(($i + 30))

	echo './bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$j' --npl 0 --sr 1 --sbr 0 --zjcr 0 --ttcr 0 --sys 1 --channel 1 --dcc 0'

	tmux new -d -s $tmux_string 'source /cvmfs/sft.cern.ch/lcg/views/LCG_96/x86_64-slc6-gcc8-opt/setup.sh; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$j' --npl 0 --sr 1 --sbr 0 --zjcr 0 --ttcr 0 --sys 1 --channel 1 --dcc 0; sleep 86400'

	sleep 180

done


sleep 300

for i in ${!ProcessNamesArray2[@]}; do

        tmux_string2="${ProcessNamesArray2[i]}${SystematicNamesArray[1]}${ChannelArray[0]}_${YearArray[0]}"

        j=$(($i + 45))

        echo './bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$j' --npl 0 --sr 1 --sbr 0 --zjcr 0 --ttcr 0 --sys 1 --channel 1 --dcc 0'

        tmux new -d -s $tmux_string2 'source /cvmfs/sft.cern.ch/lcg/views/LCG_96/x86_64-slc6-gcc8-opt/setup.sh; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$j' --npl 0 --sr 1 --sbr 0 --zjcr 0 --ttcr 0 --sys 1 --channel 1 --dcc 0; sleep 86400'

        sleep 180

done

for i in ${!ProcessNamesArray3[@]}; do

        tmux_string3="${ProcessNamesArray3[i]}${SystematicNamesArray[1]}${ChannelArray[0]}_${YearArray[0]}"

        j=$(($i + 64))

        echo './bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$j' --npl 0 --sr 1 --sbr 0 --zjcr 0 --ttcr 0 --sys 1 --channel 1 --dcc 0'

        tmux new -d -s $tmux_string3 'source /cvmfs/sft.cern.ch/lcg/views/LCG_96/x86_64-slc6-gcc8-opt/setup.sh; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$j' --npl 0 --sr 1 --sbr 0 --zjcr 0 --ttcr 0 --sys 1 --channel 1 --dcc 0; sleep 86400'

        sleep 180

done

sleep 300


echo "The tmux sessions for ttbar in the SR (ee channel for 2016, PU scale up) are running."


