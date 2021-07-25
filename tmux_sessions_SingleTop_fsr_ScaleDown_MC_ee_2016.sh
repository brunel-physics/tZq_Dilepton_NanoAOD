#!/bin/bash

ProcessNamesArray=("SingleTop_tchannel_top"             "SingleTop_tchannel_top_ScaleUp"        "SingleTop_tchannel_top_ScaleDown"      
                   "SingleTop_tchannel_antitop"         "SingleTop_tchannel_antitop_ScaleUp"    "SingleTop_tchannel_antitop_ScaleDown"  
                   "SingleTop_schannel")
		   
ProcessNamesArray2=("SingleTop_tchannel_top_hdampUP"        "SingleTop_tchannel_top_hdampDOWN"	"SingleTop_tchannel_antitop_hdampUP"    
		    "SingleTop_tchannel_antitop_hdampDOWN")  


ProcessNamesArray3=("SingleTop_tW"	           "SingleTop_tW_ScaleUp" 	 "SingleTop_tW_ScaleDown"		    
		    "SingleTop_tbarW_ScaleUp"	   "SingleTop_tbarW_ScaleDown"   "SingleTop_tHq"			    
		    "SingleTop_tZq_W_lept_Z_had"   "SingleTop_tWZ_tWll")
 			        

SystematicNamesArray=("_Nominal"                       "_PU_ScaleUp"             "_PU_ScaleDown"           "_BTag_ScaleUp"              "_BTag_ScaleDown" 
		      "_JetSmearing_ScaleUp"           "_JetSmearing_ScaleDown"  "_JetResolution_ScaleUp"  "_JetResolution_ScaleDown"   "_LeptonEfficiencies_ScaleUp"
		      "_LeptonEfficiencies_ScaleDown"  "_PDF_ScaleUp"  		 "_PDF_ScaleDown"	   "_ME_Up"			"_ME_Down"
		      "_MET_Up"			       "_MET_Down"		 "_isr_up"		   "_isr_down"			"_fsr_up"
		      "_fsr_down")
	


ChannelArray=("_ee" "_mumu" "_emu")

YearArray=("2016" "2017" "2018")

for i in ${!ProcessNamesArray[@]}; do

	tmux_string="${ProcessNamesArray[i]}${SystematicNamesArray[20]}${ChannelArray[0]}_${YearArray[0]}"

	j=$(($i + 38))

	echo './bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$j' --npl 0 --sr 1 --sbr 0 --zjcr 0 --ttcr 0 --sys 20 --channel 1 --dcc 0'

	tmux new -d -s $tmux_string 'source /cvmfs/sft.cern.ch/lcg/views/LCG_96/x86_64-slc6-gcc8-opt/setup.sh; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$j' --npl 0 --sr 1 --sbr 0 --zjcr 0 --ttcr 0 --sys 20 --channel 1 --dcc 0; sleep 86400'

	sleep 240

done

sleep 180

for i in ${!ProcessNamesArray2[@]}; do

        tmux_string="${ProcessNamesArray2[i]}${SystematicNamesArray[20]}${ChannelArray[0]}_${YearArray[0]}"

        j=$(($i + 60))

        echo './bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$j' --npl 0 --sr 1 --sbr 0 --zjcr 0 --ttcr 0 --sys 20 --channel 1 --dcc 0'

        tmux new -d -s $tmux_string 'source /cvmfs/sft.cern.ch/lcg/views/LCG_96/x86_64-slc6-gcc8-opt/setup.sh; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$j' --npl 0 --sr 1 --sbr 0 --zjcr 0 --ttcr 0 --sys 20 --channel 1 --dcc 0; sleep 86400'

        sleep 240

done

sleep 180

for i in ${!ProcessNamesArray3[@]}; do

        tmux_string="${ProcessNamesArray3[i]}${SystematicNamesArray[20]}${ChannelArray[0]}_${YearArray[0]}"

        j=$(($i + 71))

        echo './bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$j' --npl 0 --sr 1 --sbr 0 --zjcr 0 --ttcr 0 --sys 20 --channel 1 --dcc 0'

        tmux new -d -s $tmux_string 'source /cvmfs/sft.cern.ch/lcg/views/LCG_96/x86_64-slc6-gcc8-opt/setup.sh; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$j' --npl 0 --sr 1 --sbr 0 --zjcr 0 --ttcr 0 --sys 20 --channel 1 --dcc 0; sleep 86400'

        sleep 240

done

echo "The tmux sessions for single top processes in the SR (ee channel for 2016, fsr scale down) are running."


