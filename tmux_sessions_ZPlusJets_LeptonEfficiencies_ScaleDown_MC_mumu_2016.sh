#!/bin/bash

ProcessNamesArray=("ZPlusJets_M50_aMCatNLO"             "ZPlusJets_M50_aMCatNLO_ext"            "ZPlusJets_M50_Madgraph" 
                   "ZPlusJets_M50_Madgraph_ext"         "ZPlusJets_M10To50_aMCatNLO"            "ZPlusJets_M10To50_aMCatNLO_ext"
                   "ZPlusJets_M10To50_Madgraph"         "ZPlusJets_M10To50_Madgraph_ext"        "ZPlusJets_PtBinned_0To50"
                   "ZPlusJets_PtBinned_50To100"         "ZPlusJets_PtBinned_50To100_ext"        "ZPlusJets_PtBinned_100To250"
                   "ZPlusJets_PtBinned_100To250_ext1"   "ZPlusJets_PtBinned_100To250_ext2"      "ZPlusJets_PtBinned_100To250_ext5"
                   "ZPlusJets_PtBinned_250To400"        "ZPlusJets_PtBinned_250To400_ext1"      "ZPlusJets_PtBinned_250To400_ext2"
                   "ZPlusJets_PtBinned_250To400_ext5"   "ZPlusJets_PtBinned_400To650"           "ZPlusJets_PtBinned_400To650_ext1"
                   "ZPlusJets_PtBinned_400To650_ext2"   "ZPlusJets_PtBinned_650ToInf"           "ZPlusJets_PtBinned_650ToInf_ext1"
                   "ZPlusJets_PtBinned_650ToInf_ext2")  



SystematicNamesArray=("_Nominal"                       "_PU_ScaleUp"             "_PU_ScaleDown"           "_BTag_ScaleUp"              "_BTag_ScaleDown" 
		      "_JetSmearing_ScaleUp"           "_JetSmearing_ScaleDown"  "_JetResolution_ScaleUp"  "_JetResolution_ScaleDown"   "_LeptonEfficiencies_ScaleUp"
		      "_LeptonEfficiencies_ScaleDown"  "_PDF_ScaleUp"  		 "_PDF_ScaleDown"	   "_ME_Up"			"_ME_Down"
		      "_MET_Up"			       "_MET_Down"		 "_isr_up"		   "_isr_down"			"_fsr_up"
		      "_fsr_down")
	


ChannelArray=("_ee" "_mumu" "_emu")

YearArray=("2016" "2017" "2018")

for i in ${!ProcessNamesArray[@]}; do

	tmux_string="${ProcessNamesArray[i]}${SystematicNamesArray[10]}${ChannelArray[1]}_${YearArray[0]}"

	j=$(($i + 3))

	echo './bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$j' --npl 0 --sr 1 --sbr 0 --zjcr 0 --ttcr 0 --sys 10 --channel 2 --dcc 0'

	tmux new -d -s $tmux_string 'source /cvmfs/sft.cern.ch/lcg/views/LCG_96/x86_64-slc6-gcc8-opt/setup.sh; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$j' --npl 0 --sr 1 --sbr 0 --zjcr 0 --ttcr 0 --sys 10 --channel 2 --dcc 0; sleep 86400'

	sleep 180

done



echo "The tmux sessions for z+jets in the SR (mumu channel for 2016, lepton efficiencies scale down) are running."


