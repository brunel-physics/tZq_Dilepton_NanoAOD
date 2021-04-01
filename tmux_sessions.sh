#!/bin/bash

    ############################################################# 2016 (MC, mumu channel) #################################################################################

#Nominal runs

ProcessNamesArray=("tZq_mumu_2016_nominal"          "tZq_scaleup_mumu_2016_nominal"   "tZq_scaledown_mumu_2016_nominal"
		   "ZPlusJets_M50_aMCatNLO"         "ZPlusJets_M50_aMCatNLO_ext"      "ZPlusJets_M50_Madgraph"
		   "ZPlusJets_M50_Madgraph_ext"     "ZPlusJets_M10To50_aMCatNLO"      "ZPlusJets_M10To50_aMCatNLO_ext" 
		   "ZPlusJets_M10To50_Madgraph"     "ZPlusJets_M10To50_Madgraph_ext"  "ZPlusJets_PtBinned_0To50"       
                   "ZPlusJets_PtBinned_50To100"     "ZPlusJets_PtBinned_50To100_ext"  "ZPlusJets_PtBinned_100To250")


for i in ${!ProcessNamesArray[@]}; do

	if [ ${ProcessNamesArray[$i]} = "tZq_mumu_2016_nominal" ]

		then tmux new -d -s tZq_mumu_2016_nominal 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 0 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 30'

	elif [ ${ProcessNamesArray[$i]} = "tZq_scaleup_mumu_2016_nominal" ]

                then tmux new -d -s tZq_scaleup_mumu_2016_nominal 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 1 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 30'

	elif [ ${ProcessNamesArray[$i]} = "tZq_scaledown_mumu_2016_nominal" ]

                then tmux new -d -s tZq_scaledown_mumu_2016_nominal 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 2 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 30'

	elif [ ${ProcessNamesArray[$i]} = "ZPlusJets_M50_aMCatNLO" ]

                then tmux new -d -s ZPlusJets_M50_aMCatNLO_mumu_2016_nominal 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 3 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 30'

	elif [ ${ProcessNamesArray[$i]} = "ZPlusJets_M50_Madgraph" ]

                then tmux new -d -s ZPlusJets_M50_Madgraph_nominal 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 5 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 30'

	elif [ ${ProcessNamesArray[$i]} = "ZPlusJets_M50_Madgraph_ext" ]

                then tmux new -d -s ZPlusJets_M50_Madgraph_ext_nominal 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 6 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 30'

	else 
		echo "ERROR: Check ProcessNamesArray for " ${ProcessNamesArray[$i]}
		break 

	fi 
done

echo "The tmux sessions for the SBR (mumu channel 2016) are running"


   ############################################################## 2017 (MC, mumu channel) #################################################################################

   ############################################################## 2018 (MC, mumu channel) #################################################################################
