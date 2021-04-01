#!/bin/bash

ProcessNamesArray=("tZq_mumu_2016_nominal"              "tZq_scaleup_mumu_2016_nominal"         "tZq_scaledown_mumu_2016_nominal"
                   "ZPlusJets_M50_aMCatNLO"             "ZPlusJets_M50_aMCatNLO_ext"            "ZPlusJets_M50_Madgraph"
                   "ZPlusJets_M50_Madgraph_ext"         "ZPlusJets_M10To50_aMCatNLO"            "ZPlusJets_M10To50_aMCatNLO_ext"
                   "ZPlusJets_M10To50_Madgraph"         "ZPlusJets_M10To50_Madgraph_ext"        "ZPlusJets_PtBinned_0To50"
                   "ZPlusJets_PtBinned_50To100"         "ZPlusJets_PtBinned_50To100_ext"        "ZPlusJets_PtBinned_100To250"
                   "ZPlusJets_PtBinned_100To250_ext1"   "ZPlusJets_PtBinned_100To250_ext2"      "ZPlusJets_PtBinned_100To250_ext5"
                   "ZPlusJets_PtBinned_250To400"        "ZPlusJets_PtBinned_250To400_ext1"      "ZPlusJets_PtBinned_250To400_ext2"
                   "ZPlusJets_PtBinned_250To400_ext5"   "ZPlusJets_PtBinned_400To650"           "ZPlusJets_PtBinned_400To650_ext1"
                   "ZPlusJets_PtBinned_400To650_ext2"   "ZPlusJets_PtBinned_650ToInf"           "ZPlusJets_PtBinned_650ToInf_ext1"
                   "ZPlusJets_PtBinned_650ToInf_ext2"   "ttbar_2l2nu"                           "ttbar_madgraph"
                   "ttbar_madgraph_ext"                 "ttbar_TTToHadronic"                    "ttbar_TTToHadronic_ext"
                   "ttbar_TTToSemileptonic"             "ttbar_TTToSemileptonic_ext"            "ttbar_atMCaNLO"
                   "ttbar_atMCaNLO_ext"                 "ttbar_inc"                             "SingleTop_tchannel_top"
                   "SingleTop_tchannel_top_ScaleUp"     "SingleTop_tchannel_top_ScaleDown"      "SingleTop_tchannel_antitop"
                   "SingleTop_tchannel_antitop_ScaleUp" "SingleTop_tchannel_antitop_ScaleDown"  "SingleTop_schannel"
		   "ttbar_hdampUP"			"ttbar_hdampUP_ext"			"ttbar_hdampDOWN"
	           "ttbar_hdampDOWN_ext"		"TT_2l2nu_hdampUP_ext1"			"TT_2l2nu_hdampUP_ext2"
		   "TT_2l2nu_hdampDOWN"			"TT_2l2nu_hdampDOWN_ext1"		"TT_2l2nu_hdampDOWN_ext2"
		   "TTToHadronic_hdampUP"		"TTToHadronic_hdampDOWN"		"TTToSemileptonic_hdampUP"
		   "TTToSemileptonic_hdampDOWN"		"TTToSemileptonic_hdampDOWN_ext"	"SingleTop_tchannel_top_hdampUP"
		   "SingleTop_tchannel_top_hdampDOWN"	"SingleTop_tchannel_antitop_hdampUP"    "SingleTop_tchannel_antitop_hdampDOWN"
		   "ttbar_isr_UP"			"ttbar_isr_DOWN" 			"ttbar_isr_DOwN_ext"
		   "ttbar_fsr_UP"			"ttbar_fsr_UP_ext"			"ttbar_fsr_DOWN"
		   "ttbar_fsr_DOWN_ext"		        "SingleTop_tW"				"SingleTop_tW_ScaleUp"
		   "SingleTop_tW_ScaleDown"		"SingleTop_tbarW_ScaleUp"		"SingleTop_tbarW_ScaleDown"
		   "SingleTop_tHq"			"SingleTop_tZq_W_lept_Z_had"		"SingleTop_tWZ_tWll"
		   "VV_ZZTo2l2nu"			"VV_ZZTo2l2nu_ext1"			"VV_ZZTo2l2nu_ext2"
		   "VV_ZZTo2l2Q"			"VV_ZZTo4L"				"VV_ZZTo4L_ext"
		   "VV_WZTo2l2Q"			"VV_WZTo3lNu"				"VV_WZTo3lNu_ext"
		   "VV_WZTo1l1Nu2Q"			"VV_WWTo2l2Nu"				"VV_WWToLNuQQ"
		   "VV_WWToLNuQQ_ext"			"VV_WGToLNuG"				"VV_ZGToLLG_ext"
		   "VVV_WWWTo4F"			"VVV_WWWTo4F_ext"			"VVV_WWZTo4F"
		   "VVV_WWZTo4F_ext"			"VVV_WZZ"				"VVV_WZZ_ext"
		   "VVV_ZZZ"				"VVV_ZZZ_ext"				"WPlusJets"
		   "WPlusJets_ext"			"ttbarV_ttWJetsToLNu"			"ttbarV_ttWJetsToLNu_ext"
		   "ttbarV_ttWJetsToQQ"			"ttbarV_ttZToLL"			"ttbarV_ttZToLL_ext"
		   "ttbarV_ttZToLLNuNu"			"ttbarV_ttZToLLNuNu_ext"		"ttbarV_ttZToLLNuNu_ext2"
		   "ttbarV_ttZToQQ"			"ttbarV_ttZToQQ_ext"			"ttbarV_ttHTobb"
		   "ttbarV_ttHTobb_ext"			"ttbarV_ttHToNonbb")



    ############################################################# 2016 (MC, mumu channel) #################################################################################

#Nominal runs

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
	
	elif [ ${ProcessNamesArray[$i]} = "ZPlusJets_M10To50_aMCatNLO" ]

                then tmux new -d -s ZPlusJets_M10To50_aMCatNLO_nominal 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 7 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 30'


	elif [ ${ProcessNamesArray[$i]} = "ZPlusJets_M10To50_aMCatNLO_ext" ]

                then tmux new -d -s ZPlusJets_M10To50_aMCatNLO_ext_nominal 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 8 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 30'

	elif [ ${ProcessNamesArray[$i]} = "ZPlusJets_M10To50_Madgraph" ]

                then tmux new -d -s ZPlusJets_M10To50_Madgraph_nominal 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 9 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 30'


	elif [ ${ProcessNamesArray[$i]} = "ZPlusJets_PtBinned_0To50" ]

                then tmux new -d -s ZPlusJets_PtBinned_0To50_nominal 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 11 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 30'

	elif [ ${ProcessNamesArray[$i]} = "ZPlusJets_PtBinned_50To100" ]

                then tmux new -d -s ZPlusJets_PtBinned_50To100_nominal 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 12 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 30'

	elif [ ${ProcessNamesArray[$i]} = "ZPlusJets_PtBinned_50To100_ext" ]

                then tmux new -d -s ZPlusJets_PtBinned_50To100_ext_nominal 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 13 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 30'

	elif [ ${ProcessNamesArray[$i]} = "ZPlusJets_PtBinned_100To250" ]

                then tmux new -d -s ZPlusJets_PtBinned_100To250_nominal 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 14 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 30'

	elif [ ${ProcessNamesArray[$i]} = "ZPlusJets_PtBinned_100To250_ext1" ]

                then tmux new -d -s ZPlusJets_PtBinned_100To250_ext1_nominal 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 15 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 30'

	elif [ ${ProcessNamesArray[$i]} = "ZPlusJets_PtBinned_100To250_ext2" ]
        
                then tmux new -d -s ZPlusJets_PtBinned_100To250_ext2_nominal 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 16 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 30'

	elif [ ${ProcessNamesArray[$i]} = "ZPlusJets_PtBinned_100To250_ext5" ]
        
                then tmux new -d -s ZPlusJets_PtBinned_100To250_ext5_nominal 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 17 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 30'

	elif [ ${ProcessNamesArray[$i]} = "ZPlusJets_PtBinned_250To400" ]
          
                then tmux new -d -s ZPlusJets_PtBinned_250To400_nominal 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 18 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 30'

	elif [ ${ProcessNamesArray[$i]} = "ZPlusJets_PtBinned_250To400_ext1" ]
          
                then tmux new -d -s ZPlusJets_PtBinned_250To400_ext1_nominal 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 19 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 30'

	elif [ ${ProcessNamesArray[$i]} = "ZPlusJets_PtBinned_250To400_ext2" ]

                then tmux new -d -s ZPlusJets_PtBinned_250To400_ext2_nominal 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 20 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 30'

	elif [ ${ProcessNamesArray[$i]} = "ZPlusJets_PtBinned_250To400_ext5" ]

                then tmux new -d -s ZPlusJets_PtBinned_250To400_ext5_nominal 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 21 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 30'	

	elif [ ${ProcessNamesArray[$i]} = "ZPlusJets_PtBinned_400To650" ]

                then tmux new -d -s ZPlusJets_PtBinned_400To650_nominal 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 22 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 30'

	elif [ ${ProcessNamesArray[$i]} = "ZPlusJets_PtBinned_400To650_ext1" ]

                then tmux new -d -s ZPlusJets_PtBinned_400To650_ext1_nominal 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 23 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 30'

	elif [ ${ProcessNamesArray[$i]} = "ZPlusJets_PtBinned_400To650_ext2" ]

                then tmux new -d -s ZPlusJets_PtBinned_400To650_ext2_nominal 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p 24 --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 30'

	else 
		echo "ERROR: Check ProcessNamesArray for " ${ProcessNamesArray[$i]}
		break 

	fi 
done

echo "The tmux sessions for the SBR (mumu channel 2016) are running"


   ############################################################## 2017 (MC, mumu channel) #################################################################################

   ############################################################## 2018 (MC, mumu channel) #################################################################################
