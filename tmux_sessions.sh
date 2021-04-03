#!/bin/bash

ProcessNamesArray=("tZq"				"tZq_scaleup"                           "tZq_scaledown"
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


SystematicNamesArray=("_Nominal"                       "_PU_ScaleUp"             "_PU_ScaleDown"           "_BTag_ScaleUp"              "_BTag_ScaleDown" 
		      "_JetSmearing_ScaleUp"           "_JetSmearing_ScaleDown"  "_JetResolution_ScaleUp"  "_JetResolution_ScaleDown"   "_LeptonEfficiencies_ScaleUp"
		      "_LeptonEfficiencies_ScaleDown"  "_PDF_ScaleUp"  		 "_PDF_ScaleDown"	   "_ME_Up"			"_ME_Down"
		      "_MET_Up"			       "_MET_Down"		 "_isr_up"		   "_isr_down"			"_fsr_up"
		      "_fsr_down")
	

ChannelArray=("_ee" "_mumu" "_emu")

YearArray=("2016" "2017" "2018")

#Nominal runs

#2016, mumu

for i in ${!ProcessNamesArray[@]}; do

	#nominal
	tmux_string_nominal_ee_2016="${ProcessNamesArray[i]}${SystematicNamesArray[0]}${ChannelArray[0]}_${YearArray[0]}"
	tmux_string_nominal_mumu_2016="${ProcessNamesArray[i]}${SystematicNamesArray[0]}${ChannelArray[1]}_${YearArray[0]}"
	tmux_string_nominal_ee_2017="${ProcessNamesArray[i]}${SystematicNamesArray[0]}${ChannelArray[0]}_${YearArray[1]}"
        tmux_string_nominal_mumu_2017="${ProcessNamesArray[i]}${SystematicNamesArray[0]}${ChannelArray[1]}_${YearArray[1]}"	
	tmux_string_nominal_ee_2018="${ProcessNamesArray[i]}${SystematicNamesArray[0]}${ChannelArray[0]}_${YearArray[2]}"
        tmux_string_nominal_mumu_2018="${ProcessNamesArray[i]}${SystematicNamesArray[0]}${ChannelArray[1]}_${YearArray[2]}"

	tmux new -d -s $tmux_string_nominal_ee_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 1 --dcc 0; sleep 60'

	tmux new -d -s $tmux_string_nominal_mumu_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 60'

	tmux new -d -s $tmux_string_nominal_ee_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_nominal_mumu_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 60'

	tmux new -d -s $tmux_string_nominal_ee_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_nominal_mumu_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 0 --channel 2 --dcc 0; sleep 60'

	#PU scale up
	tmux_string_PU_ScaleUp_ee_2016="${ProcessNamesArray[i]}${SystematicNamesArray[1]}${ChannelArray[0]}_${YearArray[0]}"
        tmux_string_PU_ScaleUp_mumu_2016="${ProcessNamesArray[i]}${SystematicNamesArray[1]}${ChannelArray[1]}_${YearArray[0]}"
        tmux_string_PU_ScaleUp_ee_2017="${ProcessNamesArray[i]}${SystematicNamesArray[1]}${ChannelArray[0]}_${YearArray[1]}"
        tmux_string_PU_ScaleUp_mumu_2017="${ProcessNamesArray[i]}${SystematicNamesArray[1]}${ChannelArray[1]}_${YearArray[1]}"
        tmux_string_PU_ScaleUp_ee_2018="${ProcessNamesArray[i]}${SystematicNamesArray[1]}${ChannelArray[0]}_${YearArray[2]}"
        tmux_string_PU_ScaleUp_mumu_2018="${ProcessNamesArray[i]}${SystematicNamesArray[1]}${ChannelArray[1]}_${YearArray[2]}"

	tmux new -d -s $tmux_string_PU_ScaleUp_ee_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 1 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_PU_ScaleUp_mumu_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 1 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_PU_ScaleUp_ee_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 1 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_PU_ScaleUp_mumu_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 1 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_PU_ScaleUp_ee_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 1 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_PU_ScaleUp_mumu_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 1 --channel 2 --dcc 0; sleep 60'

	#PU scale down
        tmux_string_PU_ScaleDown_ee_2016="${ProcessNamesArray[i]}${SystematicNamesArray[2]}${ChannelArray[0]}_${YearArray[0]}"
        tmux_string_PU_ScaleDown_mumu_2016="${ProcessNamesArray[i]}${SystematicNamesArray[2]}${ChannelArray[1]}_${YearArray[0]}"
        tmux_string_PU_ScaleDown_ee_2017="${ProcessNamesArray[i]}${SystematicNamesArray[2]}${ChannelArray[0]}_${YearArray[1]}"
        tmux_string_PU_ScaleDown_mumu_2017="${ProcessNamesArray[i]}${SystematicNamesArray[2]}${ChannelArray[1]}_${YearArray[1]}"
        tmux_string_PU_ScaleDown_ee_2018="${ProcessNamesArray[i]}${SystematicNamesArray[2]}${ChannelArray[0]}_${YearArray[2]}"
        tmux_string_PU_ScaleDown_mumu_2018="${ProcessNamesArray[i]}${SystematicNamesArray[2]}${ChannelArray[1]}_${YearArray[2]}"

	tmux new -d -s $tmux_string_PU_ScaleDown_ee_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 2 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_PU_ScaleDown_mumu_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 2 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_PU_ScaleDown_ee_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 2 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_PU_ScaleDown_mumu_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 2 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_PU_ScaleDown_ee_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 2 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_PU_ScaleDown_mumu_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 2 --channel 2 --dcc 0; sleep 60'

	#BTag scale up
	tmux_string_BTag_ScaleUp_ee_2016="${ProcessNamesArray[i]}${SystematicNamesArray[3]}${ChannelArray[0]}_${YearArray[0]}"
        tmux_string_BTag_ScaleUp_mumu_2016="${ProcessNamesArray[i]}${SystematicNamesArray[3]}${ChannelArray[1]}_${YearArray[0]}"
        tmux_string_BTag_ScaleUp_ee_2017="${ProcessNamesArray[i]}${SystematicNamesArray[3]}${ChannelArray[0]}_${YearArray[1]}"
        tmux_string_BTag_ScaleUp_mumu_2017="${ProcessNamesArray[i]}${SystematicNamesArray[3]}${ChannelArray[1]}_${YearArray[1]}"
        tmux_string_BTag_ScaleUp_ee_2018="${ProcessNamesArray[i]}${SystematicNamesArray[3]}${ChannelArray[0]}_${YearArray[2]}"
        tmux_string_BTag_ScaleUp_mumu_2018="${ProcessNamesArray[i]}${SystematicNamesArray[3]}${ChannelArray[1]}_${YearArray[2]}"

	tmux new -d -s $tmux_string_BTag_ScaleUp_ee_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 3 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_BTag_ScaleUp_mumu_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 3 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_BTag_ScaleUp_ee_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 3 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_BTag_ScaleUp_mumu_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 3 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_BTag_ScaleUp_ee_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 3 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_BTag_ScaleUp_mumu_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 3 --channel 2 --dcc 0; sleep 60'

	#BTag scale down
        tmux_string_BTag_ScaleDown_ee_2016="${ProcessNamesArray[i]}${SystematicNamesArray[4]}${ChannelArray[0]}_${YearArray[0]}"
        tmux_string_BTag_ScaleDown_mumu_2016="${ProcessNamesArray[i]}${SystematicNamesArray[4]}${ChannelArray[1]}_${YearArray[0]}"
        tmux_string_BTag_ScaleDown_ee_2017="${ProcessNamesArray[i]}${SystematicNamesArray[4]}${ChannelArray[0]}_${YearArray[1]}"
        tmux_string_BTag_ScaleDown_mumu_2017="${ProcessNamesArray[i]}${SystematicNamesArray[4]}${ChannelArray[1]}_${YearArray[1]}"
        tmux_string_BTag_ScaleDown_ee_2018="${ProcessNamesArray[i]}${SystematicNamesArray[4]}${ChannelArray[0]}_${YearArray[2]}"
        tmux_string_BTag_ScaleDown_mumu_2018="${ProcessNamesArray[i]}${SystematicNamesArray[4]}${ChannelArray[1]}_${YearArray[2]}"

	tmux new -d -s $tmux_string_BTag_ScaleDown_ee_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 4 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_BTag_ScaleDown_mumu_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 4 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_BTag_ScaleDown_ee_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 4 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_BTag_ScaleDown_mumu_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 4 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_BTag_ScaleDown_ee_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 4 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_BTag_ScaleDown_mumu_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 4 --channel 2 --dcc 0; sleep 60'

	#Jet smearing scale up
	tmux_string_JetSmearing_ScaleUp_ee_2016="${ProcessNamesArray[i]}${SystematicNamesArray[5]}${ChannelArray[0]}_${YearArray[0]}"
        tmux_string_JetSmearing_ScaleUp_mumu_2016="${ProcessNamesArray[i]}${SystematicNamesArray[5]}${ChannelArray[1]}_${YearArray[0]}"
        tmux_string_JetSmearing_ScaleUp_ee_2017="${ProcessNamesArray[i]}${SystematicNamesArray[5]}${ChannelArray[0]}_${YearArray[1]}"
        tmux_string_JetSmearing_ScaleUp_mumu_2017="${ProcessNamesArray[i]}${SystematicNamesArray[5]}${ChannelArray[1]}_${YearArray[1]}"
        tmux_string_JetSmearing_ScaleUp_ee_2018="${ProcessNamesArray[i]}${SystematicNamesArray[5]}${ChannelArray[0]}_${YearArray[2]}"
        tmux_string_JetSmearing_ScaleUp_mumu_2018="${ProcessNamesArray[i]}${SystematicNamesArray[5]}${ChannelArray[1]}_${YearArray[2]}"

	tmux new -d -s $tmux_string_JetSmearing_ScaleUp_ee_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 5 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_JetSmearing_ScaleUp_mumu_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 5 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_JetSmearing_ScaleUp_ee_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 5 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_JetSmearing_ScaleUp_mumu_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 5 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_JetSmearing_ScaleUp_ee_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 5 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_JetSmearing_ScaleUp_mumu_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 5 --channel 2 --dcc 0; sleep 60'

	#Jet smearing scale down
        tmux_string_JetSmearing_ScaleDown_ee_2016="${ProcessNamesArray[i]}${SystematicNamesArray[6]}${ChannelArray[0]}_${YearArray[0]}"
        tmux_string_JetSmearing_ScaleDown_mumu_2016="${ProcessNamesArray[i]}${SystematicNamesArray[6]}${ChannelArray[1]}_${YearArray[0]}"
        tmux_string_JetSmearing_ScaleDown_ee_2017="${ProcessNamesArray[i]}${SystematicNamesArray[6]}${ChannelArray[0]}_${YearArray[1]}"
        tmux_string_JetSmearing_ScaleDown_mumu_2017="${ProcessNamesArray[i]}${SystematicNamesArray[6]}${ChannelArray[1]}_${YearArray[1]}"
        tmux_string_JetSmearing_ScaleDown_ee_2018="${ProcessNamesArray[i]}${SystematicNamesArray[6]}${ChannelArray[0]}_${YearArray[2]}"
        tmux_string_JetSmearing_ScaleDown_mumu_2018="${ProcessNamesArray[i]}${SystematicNamesArray[6]}${ChannelArray[1]}_${YearArray[2]}"

	tmux new -d -s $tmux_string_JetSmearing_ScaleDown_ee_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 6 --channel 1 --dcc 0; sleep 60'
        
        tmux new -d -s $tmux_string_JetSmearing_ScaleDown_mumu_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 6 --channel 2 --dcc 0; sleep 60'
        
        tmux new -d -s $tmux_string_JetSmearing_ScaleDown_ee_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 6 --channel 1 --dcc 0; sleep 60'
        
        tmux new -d -s $tmux_string_JetSmearing_ScaleDown_mumu_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 6 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_JetSmearing_ScaleDown_ee_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 6 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_JetSmearing_ScaleDown_mumu_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 6 --channel 2 --dcc 0; sleep 60'

	#Jet resolution scale up
	tmux_string_JetResolution_ScaleUp_ee_2016="${ProcessNamesArray[i]}${SystematicNamesArray[7]}${ChannelArray[0]}_${YearArray[0]}"
        tmux_string_JetResolution_ScaleUp_mumu_2016="${ProcessNamesArray[i]}${SystematicNamesArray[7]}${ChannelArray[1]}_${YearArray[0]}"
        tmux_string_JetResolution_ScaleUp_ee_2017="${ProcessNamesArray[i]}${SystematicNamesArray[7]}${ChannelArray[0]}_${YearArray[1]}"
        tmux_string_JetResolution_ScaleUp_mumu_2017="${ProcessNamesArray[i]}${SystematicNamesArray[7]}${ChannelArray[1]}_${YearArray[1]}"
        tmux_string_JetResolution_ScaleUp_ee_2018="${ProcessNamesArray[i]}${SystematicNamesArray[7]}${ChannelArray[0]}_${YearArray[2]}"
        tmux_string_JetResolution_ScaleUp_mumu_2018="${ProcessNamesArray[i]}${SystematicNamesArray[7]}${ChannelArray[1]}_${YearArray[2]}"

	tmux new -d -s $tmux_string_JetResolution_ScaleUp_ee_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 7 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_JetResolution_ScaleUp_mumu_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 7 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_JetResolution_ScaleUp_ee_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 7 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_JetResolution_ScaleUp_mumu_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 7 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_JetResolution_ScaleUp_ee_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 7 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_JetResolution_ScaleUp_mumu_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 7 --channel 2 --dcc 0; sleep 60'

	#Jet resolution scale down
        tmux_string_JetResolution_ScaleDown_ee_2016="${ProcessNamesArray[i]}${SystematicNamesArray[8]}${ChannelArray[0]}_${YearArray[0]}"
        tmux_string_JetResolution_ScaleDown_mumu_2016="${ProcessNamesArray[i]}${SystematicNamesArray[8]}${ChannelArray[1]}_${YearArray[0]}"
        tmux_string_JetResolution_ScaleDown_ee_2017="${ProcessNamesArray[i]}${SystematicNamesArray[8]}${ChannelArray[0]}_${YearArray[1]}"
        tmux_string_JetResolution_ScaleDown_mumu_2017="${ProcessNamesArray[i]}${SystematicNamesArray[8]}${ChannelArray[1]}_${YearArray[1]}"
        tmux_string_JetResolution_ScaleDown_ee_2018="${ProcessNamesArray[i]}${SystematicNamesArray[8]}${ChannelArray[0]}_${YearArray[2]}"
        tmux_string_JetResolution_ScaleDown_mumu_2018="${ProcessNamesArray[i]}${SystematicNamesArray[8]}${ChannelArray[1]}_${YearArray[2]}"

	tmux new -d -s $tmux_string_JetResolution_ScaleDown_ee_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 8 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_JetResolution_ScaleDown_mumu_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 8 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_JetResolution_ScaleDown_ee_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 8 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_JetResolution_ScaleDown_mumu_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 8 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_JetResolution_ScaleDown_ee_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 8 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_JetResolution_ScaleDown_mumu_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 8 --channel 2 --dcc 0; sleep 60'

	#Lepton efficiencies scale up
	tmux_string_LeptonEfficiencies_ScaleUp_ee_2016="${ProcessNamesArray[i]}${SystematicNamesArray[9]}${ChannelArray[0]}_${YearArray[0]}"
        tmux_string_LeptonEfficiencies_ScaleUp_mumu_2016="${ProcessNamesArray[i]}${SystematicNamesArray[9]}${ChannelArray[1]}_${YearArray[0]}"
        tmux_string_LeptonEfficiencies_ScaleUp_ee_2017="${ProcessNamesArray[i]}${SystematicNamesArray[9]}${ChannelArray[0]}_${YearArray[1]}"
        tmux_string_LeptonEfficiencies_ScaleUp_mumu_2017="${ProcessNamesArray[i]}${SystematicNamesArray[9]}${ChannelArray[1]}_${YearArray[1]}"
        tmux_string_LeptonEfficiencies_ScaleUp_ee_2018="${ProcessNamesArray[i]}${SystematicNamesArray[9]}${ChannelArray[0]}_${YearArray[2]}"
        tmux_string_LeptonEfficiencies_ScaleUp_mumu_2018="${ProcessNamesArray[i]}${SystematicNamesArray[9]}${ChannelArray[1]}_${YearArray[2]}"

	tmux new -d -s $tmux_string_LeptonEfficiencies_ScaleUp_ee_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 9 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_LeptonEfficiencies_ScaleUp_mumu_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 9 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_LeptonEfficiencies_ScaleUp_ee_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 9 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_LeptonEfficiencies_ScaleUp_mumu_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 9 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_LeptonEfficiencies_ScaleUp_ee_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 9 --channel 1 --dcc 0; sleep 60'
        
        tmux new -d -s $tmux_string_LeptonEfficiencies_ScaleUp_mumu_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 9 --channel 2 --dcc 0; sleep 60'

	#Lepton efficiencies scale down
        tmux_string_LeptonEfficiencies_ScaleDown_ee_2016="${ProcessNamesArray[i]}${SystematicNamesArray[10]}${ChannelArray[0]}_${YearArray[0]}"
        tmux_string_LeptonEfficiencies_ScaleDown_mumu_2016="${ProcessNamesArray[i]}${SystematicNamesArray[10]}${ChannelArray[1]}_${YearArray[0]}"
        tmux_string_LeptonEfficiencies_ScaleDown_ee_2017="${ProcessNamesArray[i]}${SystematicNamesArray[10]}${ChannelArray[0]}_${YearArray[1]}"
        tmux_string_LeptonEfficiencies_ScaleDown_mumu_2017="${ProcessNamesArray[i]}${SystematicNamesArray[10]}${ChannelArray[1]}_${YearArray[1]}"
        tmux_string_LeptonEfficiencies_ScaleDown_ee_2018="${ProcessNamesArray[i]}${SystematicNamesArray[10]}${ChannelArray[0]}_${YearArray[2]}"
        tmux_string_LeptonEfficiencies_ScaleDown_mumu_2018="${ProcessNamesArray[i]}${SystematicNamesArray[10]}${ChannelArray[1]}_${YearArray[2]}"

	tmux new -d -s $tmux_string_LeptonEfficiencies_ScaleDown_ee_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 10 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_LeptonEfficiencies_ScaleDown_mumu_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 10 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_LeptonEfficiencies_ScaleDown_ee_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 10 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_LeptonEfficiencies_ScaleDown_mumu_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 10 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_LeptonEfficiencies_ScaleDown_ee_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 10 --channel 1 --dcc 0; sleep 60'
        
        tmux new -d -s $tmux_string_LeptonEfficiencies_ScaleDown_mumu_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 10 --channel 2 --dcc 0; sleep 60'

	#PDF_ScaleUp
	tmux_string_PDF_ScaleUp_ee_2016="${ProcessNamesArray[i]}${SystematicNamesArray[11]}${ChannelArray[0]}_${YearArray[0]}"
        tmux_string_PDF_ScaleUp_mumu_2016="${ProcessNamesArray[i]}${SystematicNamesArray[11]}${ChannelArray[1]}_${YearArray[0]}"
        tmux_string_PDF_ScaleUp_ee_2017="${ProcessNamesArray[i]}${SystematicNamesArray[11]}${ChannelArray[0]}_${YearArray[1]}"
        tmux_string_PDF_ScaleUp_mumu_2017="${ProcessNamesArray[i]}${SystematicNamesArray[11]}${ChannelArray[1]}_${YearArray[1]}"
        tmux_string_PDF_ScaleUp_ee_2018="${ProcessNamesArray[i]}${SystematicNamesArray[11]}${ChannelArray[0]}_${YearArray[2]}"
        tmux_string_PDF_ScaleUp_mumu_2018="${ProcessNamesArray[i]}${SystematicNamesArray[11]}${ChannelArray[1]}_${YearArray[2]}"

	tmux new -d -s $tmux_string_PDF_ScaleUp_ee_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 11 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_PDF_ScaleUp_mumu_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 11 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_PDF_ScaleUp_ee_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 11 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_PDF_ScaleUp_mumu_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 11 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_PDF_ScaleUp_ee_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 11 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_PDF_ScaleUp_mumu_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 11 --channel 2 --dcc 0; sleep 60'

	#PDF_ScaleDown
        tmux_string_PDF_ScaleDown_ee_2016="${ProcessNamesArray[i]}${SystematicNamesArray[12]}${ChannelArray[0]}_${YearArray[0]}"
        tmux_string_PDF_ScaleDown_mumu_2016="${ProcessNamesArray[i]}${SystematicNamesArray[12]}${ChannelArray[1]}_${YearArray[0]}"
        tmux_string_PDF_ScaleDown_ee_2017="${ProcessNamesArray[i]}${SystematicNamesArray[12]}${ChannelArray[0]}_${YearArray[1]}"
        tmux_string_PDF_ScaleDown_mumu_2017="${ProcessNamesArray[i]}${SystematicNamesArray[12]}${ChannelArray[1]}_${YearArray[1]}"
        tmux_string_PDF_ScaleDown_ee_2018="${ProcessNamesArray[i]}${SystematicNamesArray[12]}${ChannelArray[0]}_${YearArray[2]}"
        tmux_string_PDF_ScaleDown_mumu_2018="${ProcessNamesArray[i]}${SystematicNamesArray[12]}${ChannelArray[1]}_${YearArray[2]}"

	tmux new -d -s $tmux_string_PDF_ScaleDown_ee_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 12 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_PDF_ScaleDown_mumu_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 12 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_PDF_ScaleDown_ee_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 12 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_PDF_ScaleDown_mumu_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 12 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_PDF_ScaleDown_ee_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 12 --channel 1 --dcc 0; sleep 60'
        
        tmux new -d -s $tmux_string_PDF_ScaleDown_mumu_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 12 --channel 2 --dcc 0; sleep 60'

	#ME_Up
	tmux_string_ME_Up_ee_2016="${ProcessNamesArray[i]}${SystematicNamesArray[13]}${ChannelArray[0]}_${YearArray[0]}"
        tmux_string_ME_Up_mumu_2016="${ProcessNamesArray[i]}${SystematicNamesArray[13]}${ChannelArray[1]}_${YearArray[0]}"
        tmux_string_ME_Up_ee_2017="${ProcessNamesArray[i]}${SystematicNamesArray[13]}${ChannelArray[0]}_${YearArray[1]}"
        tmux_string_ME_Up_mumu_2017="${ProcessNamesArray[i]}${SystematicNamesArray[13]}${ChannelArray[1]}_${YearArray[1]}"
        tmux_string_ME_Up_ee_2018="${ProcessNamesArray[i]}${SystematicNamesArray[13]}${ChannelArray[0]}_${YearArray[2]}"
        tmux_string_ME_Up_mumu_2018="${ProcessNamesArray[i]}${SystematicNamesArray[13]}${ChannelArray[1]}_${YearArray[2]}"

	tmux new -d -s $tmux_string_ME_Up_ee_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 13 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_ME_Up_mumu_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 13 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_ME_Up_ee_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 13 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_ME_Up_mumu_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 13 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_ME_Up_ee_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 13 --channel 1 --dcc 0; sleep 60'
        
        tmux new -d -s $tmux_string_ME_Up_mumu_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 13 --channel 2 --dcc 0; sleep 60'

	#ME_Down
	tmux_string_ME_Down_ee_2016="${ProcessNamesArray[i]}${SystematicNamesArray[14]}${ChannelArray[0]}_${YearArray[0]}"
        tmux_string_ME_Down_mumu_2016="${ProcessNamesArray[i]}${SystematicNamesArray[14]}${ChannelArray[1]}_${YearArray[0]}"
        tmux_string_ME_Down_ee_2017="${ProcessNamesArray[i]}${SystematicNamesArray[14]}${ChannelArray[0]}_${YearArray[1]}"
        tmux_string_ME_Down_mumu_2017="${ProcessNamesArray[i]}${SystematicNamesArray[14]}${ChannelArray[1]}_${YearArray[1]}"
        tmux_string_ME_Down_ee_2018="${ProcessNamesArray[i]}${SystematicNamesArray[14]}${ChannelArray[0]}_${YearArray[2]}"
        tmux_string_ME_Down_mumu_2018="${ProcessNamesArray[i]}${SystematicNamesArray[14]}${ChannelArray[1]}_${YearArray[2]}"

	tmux new -d -s $tmux_string_ME_Down_ee_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 14 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_ME_Down_mumu_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 14 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_ME_Down_ee_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 14 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_ME_Down_mumu_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 14 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_ME_Down_ee_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 14 --channel 1 --dcc 0; sleep 60'
        
        tmux new -d -s $tmux_string_ME_Down_mumu_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 14 --channel 2 --dcc 0; sleep 60'

	#MET_Up
	tmux_string_MET_Up_ee_2016="${ProcessNamesArray[i]}${SystematicNamesArray[15]}${ChannelArray[0]}_${YearArray[0]}"
        tmux_string_MET_Up_mumu_2016="${ProcessNamesArray[i]}${SystematicNamesArray[15]}${ChannelArray[1]}_${YearArray[0]}"
        tmux_string_MET_Up_ee_2017="${ProcessNamesArray[i]}${SystematicNamesArray[15]}${ChannelArray[0]}_${YearArray[1]}"
        tmux_string_MET_Up_mumu_2017="${ProcessNamesArray[i]}${SystematicNamesArray[15]}${ChannelArray[1]}_${YearArray[1]}"
        tmux_string_MET_Up_ee_2018="${ProcessNamesArray[i]}${SystematicNamesArray[15]}${ChannelArray[0]}_${YearArray[2]}"
        tmux_string_MET_Up_mumu_2018="${ProcessNamesArray[i]}${SystematicNamesArray[15]}${ChannelArray[1]}_${YearArray[2]}"

	tmux new -d -s $tmux_string_MET_Up_ee_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 15 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_MET_Up_mumu_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 15 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_MET_Up_ee_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 15 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_MET_Up_mumu_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 15 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_MET_Up_ee_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 15 --channel 1 --dcc 0; sleep 60'
        
        tmux new -d -s $tmux_string_MET_Up_mumu_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 15 --channel 2 --dcc 0; sleep 60'

        #MET_Down
	tmux_string_MET_Down_ee_2016="${ProcessNamesArray[i]}${SystematicNamesArray[16]}${ChannelArray[0]}_${YearArray[0]}"
        tmux_string_MET_Down_mumu_2016="${ProcessNamesArray[i]}${SystematicNamesArray[16]}${ChannelArray[1]}_${YearArray[0]}"
        tmux_string_MET_Down_ee_2017="${ProcessNamesArray[i]}${SystematicNamesArray[16]}${ChannelArray[0]}_${YearArray[1]}"
        tmux_string_MET_Down_mumu_2017="${ProcessNamesArray[i]}${SystematicNamesArray[16]}${ChannelArray[1]}_${YearArray[1]}"
        tmux_string_MET_Down_ee_2018="${ProcessNamesArray[i]}${SystematicNamesArray[16]}${ChannelArray[0]}_${YearArray[2]}"
        tmux_string_MET_Down_mumu_2018="${ProcessNamesArray[i]}${SystematicNamesArray[16]}${ChannelArray[1]}_${YearArray[2]}"

	tmux new -d -s $tmux_string_MET_Down_ee_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 16 --channel 1 --dcc 0; sleep 60'
        
        tmux new -d -s $tmux_string_MET_Down_mumu_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 16 --channel 2 --dcc 0; sleep 60'
        
        tmux new -d -s $tmux_string_MET_Down_ee_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 16 --channel 1 --dcc 0; sleep 60'
        
        tmux new -d -s $tmux_string_MET_Down_mumu_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 16 --channel 2 --dcc 0; sleep 60'
        
        tmux new -d -s $tmux_string_MET_Down_ee_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 16 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_MET_Down_mumu_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 16 --channel 2 --dcc 0; sleep 60'

	#isr up
	tmux_string_isr_up_ee_2016="${ProcessNamesArray[i]}${SystematicNamesArray[17]}${ChannelArray[0]}_${YearArray[0]}"
        tmux_string_isr_up_mumu_2016="${ProcessNamesArray[i]}${SystematicNamesArray[17]}${ChannelArray[1]}_${YearArray[0]}"
        tmux_string_isr_up_ee_2017="${ProcessNamesArray[i]}${SystematicNamesArray[17]}${ChannelArray[0]}_${YearArray[1]}"
        tmux_string_isr_up_mumu_2017="${ProcessNamesArray[i]}${SystematicNamesArray[17]}${ChannelArray[1]}_${YearArray[1]}"
        tmux_string_isr_up_ee_2018="${ProcessNamesArray[i]}${SystematicNamesArray[17]}${ChannelArray[0]}_${YearArray[2]}"
        tmux_string_isr_up_mumu_2018="${ProcessNamesArray[i]}${SystematicNamesArray[17]}${ChannelArray[1]}_${YearArray[2]}"

	tmux new -d -s $tmux_string_isr_up_ee_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 17 --channel 1 --dcc 0; sleep 60'
        
        tmux new -d -s $tmux_string_isr_up_mumu_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 17 --channel 2 --dcc 0; sleep 60'
        
        tmux new -d -s $tmux_string_isr_up_ee_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 17 --channel 1 --dcc 0; sleep 60'
        
        tmux new -d -s $tmux_string_isr_up_mumu_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 17 --channel 2 --dcc 0; sleep 60'
        
        tmux new -d -s $tmux_string_isr_up_ee_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 17 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_isr_up_mumu_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 17 --channel 2 --dcc 0; sleep 60'

	#isr down
	tmux_string_isr_down_ee_2016="${ProcessNamesArray[i]}${SystematicNamesArray[18]}${ChannelArray[0]}_${YearArray[0]}"
        tmux_string_isr_down_mumu_2016="${ProcessNamesArray[i]}${SystematicNamesArray[18]}${ChannelArray[1]}_${YearArray[0]}"
        tmux_string_isr_down_ee_2017="${ProcessNamesArray[i]}${SystematicNamesArray[18]}${ChannelArray[0]}_${YearArray[1]}"
        tmux_string_isr_down_mumu_2017="${ProcessNamesArray[i]}${SystematicNamesArray[18]}${ChannelArray[1]}_${YearArray[1]}"
        tmux_string_isr_down_ee_2018="${ProcessNamesArray[i]}${SystematicNamesArray[18]}${ChannelArray[0]}_${YearArray[2]}"
        tmux_string_isr_down_mumu_2018="${ProcessNamesArray[i]}${SystematicNamesArray[18]}${ChannelArray[1]}_${YearArray[2]}"

	tmux new -d -s $tmux_string_isr_down_ee_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 18 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_isr_down_mumu_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 18 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_isr_down_ee_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 18 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_isr_down_mumu_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 18 --channel 2 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_isr_down_ee_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 18 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_isr_down_mumu_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 18 --channel 2 --dcc 0; sleep 60'

	#fsr up
	tmux_string_fsr_up_ee_2016="${ProcessNamesArray[i]}${SystematicNamesArray[19]}${ChannelArray[0]}_${YearArray[0]}"
        tmux_string_fsr_up_mumu_2016="${ProcessNamesArray[i]}${SystematicNamesArray[19]}${ChannelArray[1]}_${YearArray[0]}"
        tmux_string_fsr_up_ee_2017="${ProcessNamesArray[i]}${SystematicNamesArray[19]}${ChannelArray[0]}_${YearArray[1]}"
        tmux_string_fsr_up_mumu_2017="${ProcessNamesArray[i]}${SystematicNamesArray[19]}${ChannelArray[1]}_${YearArray[1]}"
        tmux_string_fsr_up_ee_2018="${ProcessNamesArray[i]}${SystematicNamesArray[19]}${ChannelArray[0]}_${YearArray[2]}"
        tmux_string_fsr_up_mumu_2018="${ProcessNamesArray[i]}${SystematicNamesArray[19]}${ChannelArray[1]}_${YearArray[2]}"

	tmux new -d -s $tmux_string_fsr_up_ee_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 19 --channel 1 --dcc 0; sleep 60'
        
        tmux new -d -s $tmux_string_fsr_up_mumu_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 19 --channel 2 --dcc 0; sleep 60'
        
        tmux new -d -s $tmux_string_fsr_up_ee_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 19 --channel 1 --dcc 0; sleep 60'
        
        tmux new -d -s $tmux_string_fsr_up_mumu_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 19 --channel 2 --dcc 0; sleep 60'
        
        tmux new -d -s $tmux_string_fsr_up_ee_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 19 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_fsr_up_mumu_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 19 --channel 2 --dcc 0; sleep 60'

        #fsr down
	tmux_string_fsr_down_ee_2016="${ProcessNamesArray[i]}${SystematicNamesArray[20]}${ChannelArray[0]}_${YearArray[0]}"
        tmux_string_fsr_down_mumu_2016="${ProcessNamesArray[i]}${SystematicNamesArray[20]}${ChannelArray[1]}_${YearArray[0]}"
        tmux_string_fsr_down_ee_2017="${ProcessNamesArray[i]}${SystematicNamesArray[20]}${ChannelArray[0]}_${YearArray[1]}"
        tmux_string_fsr_down_mumu_2017="${ProcessNamesArray[i]}${SystematicNamesArray[20]}${ChannelArray[1]}_${YearArray[1]}"
        tmux_string_fsr_down_ee_2018="${ProcessNamesArray[i]}${SystematicNamesArray[20]}${ChannelArray[0]}_${YearArray[2]}"
        tmux_string_fsr_down_mumu_2018="${ProcessNamesArray[i]}${SystematicNamesArray[20]}${ChannelArray[1]}_${YearArray[2]}"
		
	tmux new -d -s $tmux_string_fsr_down_ee_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 20 --channel 1 --dcc 0; sleep 60'
        
        tmux new -d -s $tmux_string_fsr_down_mumu_2016 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2016 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 20 --channel 2 --dcc 0; sleep 60'
        
        tmux new -d -s $tmux_string_fsr_down_ee_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 20 --channel 1 --dcc 0; sleep 60'
        
        tmux new -d -s $tmux_string_fsr_down_mumu_2017 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2017 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 20 --channel 2 --dcc 0; sleep 60'
    
        tmux new -d -s $tmux_string_fsr_down_ee_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 20 --channel 1 --dcc 0; sleep 60'

        tmux new -d -s $tmux_string_fsr_down_mumu_2018 'source ~/.bashrc; make clean; make; ./bin/fulleventselectionMain.exe --mc 1 -y 2018 -p '$i' --npl 0 --sr 1 --sbr 1 --zjcr 0 --ttcr 0 --sys 20 --channel 2 --dcc 0; sleep 60'

done



echo "The tmux sessions for the SBR (ee and mumu channels for 2016, 2017 and 2018) are running."


