#include <iostream>
#include <fstream>

using namespace std;

auto Normalisation_Calculation(const string& process, const string& year){

  float int_lumi, cross_section, NumberOfSimEvents;

  if(year == "2016"){

	int_lumi = 35883;

	if(process == "tZq"){NumberOfSimEvents = 13656784; cross_section = 0.0758;}
	else if(process == "ZPlusJets_M50_aMCatNLO"){NumberOfSimEvents = 120777245; cross_section = 5941.0;}
	else if(process == "ZPlusJets_M50_aMCatNLO_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_M50_Madgraph"){NumberOfSimEvents = 49748967; cross_section = 4963.0;}
        else if(process == "ZPlusJets_M50_Madgraph_ext"){NumberOfSimEvents = 96531428; cross_section = 4963.0;}
	else if(process == "ZPlusJets_M10To50_aMCatNLO"){NumberOfSimEvents = 67942840; cross_section = 18810.0;}
	else if(process == "ZPlusJets_M10To50_aMCatNLO_ext"){NumberOfSimEvents = 40154170; cross_section = 18810.0;}
	else if(process == "ZPlusJets_M10To50_Madgraph"){NumberOfSimEvents = 35114961; cross_section = 16270.0;}
        else if(process == "ZPlusJets_M10To50_Madgraph_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_0To50"){NumberOfSimEvents = 37458375; cross_section = 5352.57924;}
	else if(process == "ZPlusJets_PtBinned_50To100"){NumberOfSimEvents = 21847075; cross_section = 36381428;}
	else if(process == "ZPlusJets_PtBinned_50To100_ext"){NumberOfSimEvents = 108670239; cross_section = 363.81428;}
	else if(process == "ZPlusJets_PtBinned_100To250"){NumberOfSimEvents = 2046961; cross_section = 84.014804;}
	else if(process == "ZPlusJets_PtBinned_100To250_ext1"){NumberOfSimEvents = 2805972; cross_section = 84.014804;}
	else if(process == "ZPlusJets_PtBinned_100To250_ext2"){NumberOfSimEvents = 2991815; cross_section = 84.014804;}
	else if(process == "ZPlusJets_PtBinned_100To250_ext5"){NumberOfSimEvents = 76440229; cross_section = 84.014804;}
	else if(process == "ZPlusJets_PtBinned_250To400"){NumberOfSimEvents = 423976; cross_section = 3.228256512;}
	else if(process == "ZPlusJets_PtBinned_250To400_ext1"){NumberOfSimEvents = 590806; cross_section = 3.228256512;}
	else if(process == "ZPlusJets_PtBinned_250To400_ext2"){NumberOfSimEvents = 594317; cross_section = 3.228256512;}
	else if(process == "ZPlusJets_PtBinned_250To400_ext5"){NumberOfSimEvents = 19567800; cross_section = 3.228256512;}
	else if(process == "ZPlusJets_PtBinned_400To650"){NumberOfSimEvents = 432056; cross_section = 0.436041144;}
	else if(process == "ZPlusJets_PtBinned_400To650_ext1"){NumberOfSimEvents = 589842; cross_section = 0.436041144;}
	else if(process == "ZPlusJets_PtBinned_400To650_ext2"){NumberOfSimEvents = 604038; cross_section = 0.436041144;}
	else if(process == "ZPlusJets_PtBinned_650ToInf"){NumberOfSimEvents = 430691; cross_section = 0.040981055;}
        else if(process == "ZPlusJets_PtBinned_650ToInf_ext1"){NumberOfSimEvents = 599665; cross_section = 0.040981055;}
        else if(process == "ZPlusJets_PtBinned_650ToInf_ext2"){NumberOfSimEvents = 597526; cross_section = 0.040981055;}
	else if(process == "ttbar_2l2nu"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbar_madgraph"){NumberOfSimEvents = 6068369; cross_section = 56.86;}	
	else if(process == "ttbar_madgraph_ext"){NumberOfSimEvents = 24767666; cross_section = 56.86;}
	else if(process == "ttbar_TTToHadronic"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbar_TTToSemileptonic"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbar_atMCaNLO"){NumberOfSimEvents = 43768838; cross_section = 722.8;}
	else if(process == "SingleTop_tchannel_top"){NumberOfSimEvents = 67105876; cross_section = 136.02;}
        else if(process == "SingleTop_tchannel_top_ScaleUp"){NumberOfSimEvents = 5992440; cross_section = 136.02;}
	else if(process == "SingleTop_tchannel_top_ScaleDown"){NumberOfSimEvents = 64352832; cross_section = 136.02;} 	
	else if(process == "SingleTop_tchannel_antitop"){NumberOfSimEvents = 38811017; cross_section = 80.95;}
	else if(process == "SingleTop_tchannel_antitop_ScaleUp"){NumberOfSimEvents = 3970546; cross_section = 80.95;}
	else if(process == "SingleTop_tchannel_antitop_ScaleDown"){NumberOfSimEvents = 37359247; cross_section = 80.95;}
	else if(process == "SingleTop_schannel"){NumberOfSimEvents = 2989199; cross_section = 10.12;}
	else if(process == "ttbar_hdampUP"){NumberOfSimEvents = 29833668; cross_section = 730.6;}
	else if(process == "ttbar_hdampUP_ext"){NumberOfSimEvents = 28855428; cross_section = 730.6;}
	else if(process == "ttbar_hdampDOWN"){NumberOfSimEvents = 29047858; cross_section = 730.7;}
	else if(process == "ttbar_hdampDOWN_ext"){NumberOfSimEvents = 29229088; cross_section = 730.7;}
	else if(process == "SingleTop_tchannel_top_hdampUP"){NumberOfSimEvents = 0; cross_section = 0;} //need to check
	else if(process == "SingleTop_tchannel_top_hdampDOWN"){NumberOfSimEvents = 0; cross_section = 0;} //need to check
	else if(process == "ttbar_isr_UP"){NumberOfSimEvents = 58977100; cross_section = 730.6;}
	else if(process == "ttbar_isr_DOWN"){NumberOfSimEvents = 28409782; cross_section = 730.6;}
	else if(process == "ttbar_isr_DOWN_ext"){NumberOfSimEvents = 29915551; cross_section = 730.6;}
	else if(process == "ttbar_fsr_UP"){NumberOfSimEvents = 29632372; cross_section = 730.6;}
	else if(process == "ttbar_fsr_UP_ext"){NumberOfSimEvents = 29501065; cross_section = 730.6;}
	else if(process == "ttbar_fsr_DOWN"){NumberOfSimEvents = 29571600; cross_section = 730.6;}
        else if(process == "ttbar_fsr_DOWN_ext"){NumberOfSimEvents = 29571600; cross_section = 730.6;}
        else if(process == "SingleTop_tW"){NumberOfSimEvents = 6952830; cross_section = 38.09;}
	else if(process == "SingleTop_tW_ScaleUp"){NumberOfSimEvents = 997880; cross_section = 38.09;}
	else if(process == "SingleTop_tW_ScaleDown"){NumberOfSimEvents = 993640; cross_section = 38.09;}
	else if(process == "SingleTop_tbarW"){NumberOfSimEvents = 6933094; cross_section = 38.06;}
	else if(process == "SingleTop_tbarW_ScaleUp"){NumberOfSimEvents = 1000000; cross_section = 38.06;}
	else if(process == "SingleTop_tbarW_ScaleDown"){NumberOfSimEvents = 999068; cross_section = 38.06;}
	else if(process == "SingleTop_tHq"){NumberOfSimEvents = 3495799; cross_section = 0.2609;}
	else if(process == "SingleTop_tZq_W_lept_Z_had"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "SingleTop_tWZ_tWll"){NumberOfSimEvents = 50000; cross_section = 0.01104;}
	else if(process == "VV_ZZTo2l2nu"){NumberOfSimEvents = 8931750; cross_section = 0.5644;}
	else if(process == "VV_ZZTo2l2nu_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "VV_ZZTo2l2Q"){NumberOfSimEvents = 15462693; cross_section = 3.222;}
	else if(process == "VV_ZZTo4L"){NumberOfSimEvents = 10711278; cross_section = 1.204;}
	else if(process == "VV_WW1nuqq"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "VV_WZTo2l2Q"){NumberOfSimEvents = 26517272; cross_section = 5.606;}
	else if(process == "VV_WZTo3lNu"){NumberOfSimEvents = 1959179; cross_section = 4.688;}
        else if(process == "VV_WZTo1l2Nu2Q"){NumberOfSimEvents = 24311445; cross_section = 10.73;}
        else if(process == "VV_WWTo2l2Nu"){NumberOfSimEvents = 1999000; cross_section = 10.48;}
	else if(process == "VV_WWToLNuQQ"){NumberOfSimEvents = 1999200; cross_section = 43.53;}
	else if(process == "VV_WWToLNuQQ_ext"){NumberOfSimEvents = 6655400; cross_section = 43.53;}
        else if(process == "VV_WGToLNuG"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "VV_ZGToLLG"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "VVV_WWWTo4F"){NumberOfSimEvents = 240000; cross_section = 0.2086;}
	else if(process == "VVV_WWZTo4F"){NumberOfSimEvents = 250000; cross_section = 0.1651;}
	else if(process == "VVV_WZZ"){NumberOfSimEvents = 246800; cross_section = 0.05565;}
	else if(process == "VVV_ZZZ"){NumberOfSimEvents = 249237; cross_section = 0.01398;}
        else if(process == "WPlusJets"){NumberOfSimEvents = 22533326; cross_section = 60430.0;}
        else if(process == "WPlusJets_ext"){NumberOfSimEvents = 237263153; cross_section = 60430.0;}
        else if(process == "ttbarV_ttWJetsToLNu"){NumberOfSimEvents = 2160168; cross_section = 0.2001;}
	else if(process == "ttbarV_ttWJetsToLNu_ext"){NumberOfSimEvents = 3120397; cross_section = 0.2001;}
	else if(process == "ttbarV_ttWJetsToQQ"){NumberOfSimEvents = 833298; cross_section = 0.405;}
	else if(process == "ttbarV_ttZToLL"){NumberOfSimEvents = 1992438; cross_section = 0.2529;}
	else if(process == "ttbarV_ttZToLL_ext2"){NumberOfSimEvents = 5837781; cross_section = 0.2529;}
	else if(process == "ttbarV_ttZToLL_ext3"){NumberOfSimEvents = 5934228; cross_section = 0.2529;}
	else if(process == "ttbarV_ttZToQQ"){NumberOfSimEvents = 749400; cross_section = 0.5297;}
	else if(process == "ttbarV_ttZToQQ_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbarV_ttgamma"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbarV_ttgamma_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbarV_ttHTobb"){NumberOfSimEvents = 3872944; cross_section = 0.5638;}
	else if(process == "ttbarV_ttHToNonbb"){NumberOfSimEvents = 3981250; cross_section = 0.5638;}
	else{cout << "Double check the input process you have entered" << endl;}

  }
  else if(year == "2017"){

	int_lumi = 41528;

	if(process == "tZq"){NumberOfSimEvents = 13276146; cross_section = 0.07358;}
	else if(process == "ZPlusJets_M50_aMCatNLO"){NumberOfSimEvents = 27529915; cross_section = 6529.0;}
	else if(process == "ZPlusJets_M50_aMCatNLO_ext"){NumberOfSimEvents = 182104014; cross_section = 6529.0;}
	else if(process == "ZPlusJets_M50_Madgraph"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ZPlusJets_M50_Madgraph_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_M10To50_aMCatNLO"){NumberOfSimEvents = 316134; cross_section = 15810;}
	else if(process == "ZPlusJets_M10To50_aMCatNLO_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_M10To50_Madgraph"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ZPlusJets_M10To50_Madgraph_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_0To50"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_50To100"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_50To100_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_100To250"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_100To250_ext1"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_100To250_ext2"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_100To250_ext5"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_250To400"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_250To400_ext1"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_250To400_ext2"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_250To400_ext5"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_400To650"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_400To650_ext1"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_400To650_ext2"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_650ToInf"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ZPlusJets_PtBinned_650ToInf_ext1"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ZPlusJets_PtBinned_650ToInf_ext2"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbar_2l2nu"){NumberOfSimEvents = 69098644; cross_section = 88.29;}
	else if(process == "ttbar_madgraph"){NumberOfSimEvents = 6094476; cross_section = 56.86;}	
	else if(process == "ttbar_madgraph_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbar_TTToHadronic"){NumberOfSimEvents = 130091218; cross_section = 377.96;}
	else if(process == "ttbar_TTToSemileptonic"){NumberOfSimEvents = 110014744; cross_section = 365.34;}
	else if(process == "ttbar_atMCaNLO"){NumberOfSimEvents = 154280331; cross_section = 722.8;}
	else if(process == "SingleTop_tchannel_top"){NumberOfSimEvents = 122630600; cross_section = 136.02;}
        else if(process == "SingleTop_tchannel_top_ScaleUp"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "SingleTop_tchannel_top_ScaleDown"){NumberOfSimEvents = 0; cross_section = 0;} 	
	else if(process == "SingleTop_tchannel_antitop"){NumberOfSimEvents = 63620800; cross_section = 80.95;}
	else if(process == "SingleTop_tchannel_antitop_ScaleUp"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "SingleTop_tchannel_antitop_ScaleDown"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "SingleTop_schannel"){NumberOfSimEvents = 9883805; cross_section = 3.74;}
	else if(process == "ttbar_hdampUP"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbar_hdampUP_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbar_hdampDOWN"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbar_hdampDOWN_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "SingleTop_tchannel_top_hdampUP"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "SingleTop_tchannel_top_hdampDOWN"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbar_isr_UP"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbar_isr_DOWN"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbar_isr_DOWN_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbar_fsr_UP"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbar_fsr_UP_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbar_fsr_DOWN"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbar_fsr_DOWN_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "SingleTop_tW"){NumberOfSimEvents = 7945242; cross_section = 34.91;}
	else if(process == "SingleTop_tW_ScaleUp"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "SingleTop_tW_ScaleDown"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "SingleTop_tbarW"){NumberOfSimEvents = 7745276; cross_section = 34.91;}
	else if(process == "SingleTop_tbarW_ScaleUp"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "SingleTop_tbarW_ScaleDown"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "SingleTop_tHq"){NumberOfSimEvents = 3381548; cross_section = 0.3184;}
	else if(process == "SingleTop_tZq_W_lept_Z_had"){NumberOfSimEvents = 1000000; cross_section = 0.1573;}
	else if(process == "SingleTop_tWZ_tWll"){NumberOfSimEvents = 986000; cross_section = 0.01103;}
	else if(process == "VV_ZZTo2l2nu"){NumberOfSimEvents = 8744768; cross_section = 0.5644;}
	else if(process == "VV_ZZTo2l2nu_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "VV_ZZTo2l2Q"){NumberOfSimEvents = 27611672; cross_section = 3.222;}
	else if(process == "VV_ZZTo4L"){NumberOfSimEvents = 6964071; cross_section = 1.256;}
	else if(process == "VV_WW1nuqq"){NumberOfSimEvents = 8785360; cross_section = 45.99;}
	else if(process == "VV_WZTo2l2Q"){NumberOfSimEvents = 27582164; cross_section = 5.606;}
	else if(process == "VV_WZTo3lNu"){NumberOfSimEvents = 10987679; cross_section = 5.052;}
        else if(process == "VV_WZTo1l2Nu2Q"){NumberOfSimEvents = 4997672; cross_section = 45.68;}
        else if(process == "VV_WWTo2l2Nu"){NumberOfSimEvents = 2000000; cross_section = 11.08;}
	else if(process == "VV_WWToLNuQQ"){NumberOfSimEvents = 8785360; cross_section = 45.99;}
	else if(process == "VV_WWToLNuQQ_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "VV_WGToLNuG"){NumberOfSimEvents = 6283083; cross_section = 405.27;}
	else if(process == "VV_ZGToLLG"){NumberOfSimEvents = 30490034; cross_section = 51.50;}
        else if(process == "VVV_WWWTo4F"){NumberOfSimEvents = 232300; cross_section = 0.2086;}
	else if(process == "VVV_WWZTo4F"){NumberOfSimEvents = 250000; cross_section = 0.1651;}
	else if(process == "VVV_WZZ"){NumberOfSimEvents = 250000; cross_section = 0.05565;}
	else if(process == "VVV_ZZZ"){NumberOfSimEvents = 250000; cross_section = 0.01398;}
        else if(process == "WPlusJets"){NumberOfSimEvents = 30008250; cross_section = 52940.0;}
        else if(process == "WPlusJets_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbarV_ttWJetsToLNu"){NumberOfSimEvents = 4908905; cross_section = 0.2198;}
	else if(process == "ttbarV_ttWJetsToLNu_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbarV_ttWJetsToQQ"){NumberOfSimEvents = 811306; cross_section = 0.4316;}
	else if(process == "ttbarV_ttZToLL"){NumberOfSimEvents = 250000; cross_section = 0.05324;}
	else if(process == "ttbarV_ttZToLL_ext2"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbarV_ttZToLL_ext3"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbarV_ttZToQQ"){NumberOfSimEvents = 750000; cross_section = 0.5104;}
        else if(process == "ttbarV_ttZToQQ_ext"){NumberOfSimEvents = 8940000; cross_section = 0.5104;}
	else if(process == "ttbarV_ttgamma"){NumberOfSimEvents = 4642344; cross_section = 0.5804;}
	else if(process == "ttbarV_ttgamma_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbarV_ttHTobb"){NumberOfSimEvents = 8000000; cross_section = 0.5269;}
	else if(process == "ttbarV_ttHToNonbb"){NumberOfSimEvents = 7966779; cross_section = 0.5638;}
	else{cout << "Double check the input process you have entered" << endl;}	

  }
  else if(year == "2018"){

	int_lumi = 59688;

	if(process == "tZq"){NumberOfSimEvents = 13736000; cross_section = 0.07358;}
	else if(process == "ZPlusJets_M50_aMCatNLO"){NumberOfSimEvents = 997561; cross_section = 6529.0;}
	else if(process == "ZPlusJets_M50_aMCatNLO_ext"){NumberOfSimEvents = 193094040; cross_section = 6529.0;}
	else if(process == "ZPlusJets_M50_Madgraph"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ZPlusJets_M50_Madgraph_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_M10To50_aMCatNLO"){NumberOfSimEvents = 39392062; cross_section = 15810;}
	else if(process == "ZPlusJets_M10To50_aMCatNLO_ext"){NumberOfSimEvents = 46976952; cross_section = 15810;}
	else if(process == "ZPlusJets_M10To50_Madgraph"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ZPlusJets_M10To50_Madgraph_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_0To50"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_50To100"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_50To100_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_100To250"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_100To250_ext1"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_100To250_ext2"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_100To250_ext5"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_250To400"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_250To400_ext1"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_250To400_ext2"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_250To400_ext5"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_400To650"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_400To650_ext1"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_400To650_ext2"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_650ToInf"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ZPlusJets_PtBinned_650ToInf_ext1"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ZPlusJets_PtBinned_650ToInf_ext2"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbar_2l2nu"){NumberOfSimEvents = 64310000; cross_section = 88.29;}
	else if(process == "ttbar_madgraph"){NumberOfSimEvents = 28701360; cross_section = 54.23;}	
	else if(process == "ttbar_madgraph_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbar_TTToHadronic"){NumberOfSimEvents = 133664000; cross_section = 377.96;}
	else if(process == "ttbar_TTToSemileptonic"){NumberOfSimEvents = 101550000; cross_section = 365.34;}
	else if(process == "ttbar_atMCaNLO"){NumberOfSimEvents = 142155064; cross_section = 831.76;}
	else if(process == "SingleTop_tchannel_top"){NumberOfSimEvents = 154307600; cross_section = 136.02;}
        else if(process == "SingleTop_tchannel_top_ScaleUp"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "SingleTop_tchannel_top_ScaleDown"){NumberOfSimEvents = 0; cross_section = 0;} 	
	else if(process == "SingleTop_tchannel_antitop"){NumberOfSimEvents = 79090800; cross_section = 80.95;}
	else if(process == "SingleTop_tchannel_antitop_ScaleUp"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "SingleTop_tchannel_antitop_ScaleDown"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "SingleTop_schannel"){NumberOfSimEvents = 19965000; cross_section = 3.74;}
	else if(process == "ttbar_hdampUP"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbar_hdampUP_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbar_hdampDOWN"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbar_hdampDOWN_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "SingleTop_tchannel_top_hdampUP"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "SingleTop_tchannel_top_hdampDOWN"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbar_isr_UP"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbar_isr_DOWN"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbar_isr_DOWN_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbar_fsr_UP"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbar_fsr_UP_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbar_fsr_DOWN"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbar_fsr_DOWN_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "SingleTop_tW"){NumberOfSimEvents = 9598000; cross_section = 34.91;}
	else if(process == "SingleTop_tW_ScaleUp"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "SingleTop_tW_ScaleDown"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "SingleTop_tbarW"){NumberOfSimEvents = 7623000; cross_section = 34.97;}
	else if(process == "SingleTop_tbarW_ScaleUp"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "SingleTop_tbarW_ScaleDown"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "SingleTop_tHq"){NumberOfSimEvents = 3375995; cross_section = 0.3184;}
	else if(process == "SingleTop_tZq_W_lept_Z_had"){NumberOfSimEvents = 4977000; cross_section = 0.1518;}
	else if(process == "SingleTop_tWZ_tWll"){NumberOfSimEvents = 248600; cross_section = 0.01103;}
	else if(process == "VV_ZZTo2l2nu"){NumberOfSimEvents = 8382600; cross_section = 0.5644;}
	else if(process == "VV_ZZTo2l2nu_ext"){NumberOfSimEvents = 48046000; cross_section = 0.5644;}
        else if(process == "VV_ZZTo2l2Q"){NumberOfSimEvents = 27900469; cross_section = 3.222;}
	else if(process == "VV_ZZTo4L"){NumberOfSimEvents = 99009000; cross_section = 1.256;}
	else if(process == "VV_WW1nuqq"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "VV_WZTo2l2Q"){NumberOfSimEvents = 28193648; cross_section = 5.606;}
	else if(process == "VV_WZTo3lNu"){NumberOfSimEvents = 10749269; cross_section = 5.052;}
        else if(process == "VV_WZTo1l2Nu2Q"){NumberOfSimEvents = 18901469; cross_section = 10.73;}
        else if(process == "VV_WWTo2l2Nu"){NumberOfSimEvents = 7758900; cross_section = 11.08;}
	else if(process == "VV_WWToLNuQQ"){NumberOfSimEvents = 19199100; cross_section = 45.68;}
	else if(process == "VV_WWToLNuQQ_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "VV_WGToLNuG"){NumberOfSimEvents = 6108186; cross_section = 405.27;}
	else if(process == "VV_ZGToLLG"){NumberOfSimEvents = 13946364; cross_section = 51.50;}
        else if(process == "VVV_WWWTo4F"){NumberOfSimEvents = 240000; cross_section = 0.2086;}
	else if(process == "VVV_WWZTo4F"){NumberOfSimEvents = 250000; cross_section = 0.1651;}
	else if(process == "VVV_WZZ"){NumberOfSimEvents = 250000; cross_section = 0.05565;}
	else if(process == "VVV_ZZZ"){NumberOfSimEvents = 250000; cross_section = 0.01368;}
        else if(process == "WPlusJets"){NumberOfSimEvents = 71026861; cross_section = 52940.0;}
        else if(process == "WPlusJets_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbarV_ttWJetsToLNu"){NumberOfSimEvents = 4911941; cross_section = 0.2149;}
	else if(process == "ttbarV_ttWJetsToLNu_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbarV_ttWJetsToQQ"){NumberOfSimEvents = 835296; cross_section = 0.4316;}
	else if(process == "ttbarV_ttZToLL"){NumberOfSimEvents = 250000; cross_section = 0.05324;}
	else if(process == "ttbarV_ttZToLL_ext2"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbarV_ttZToLL_ext3"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ttbarV_ttZToQQ"){NumberOfSimEvents = 750000; cross_section = 0.5104;}
        else if(process == "ttbarV_ttZToQQ_ext"){NumberOfSimEvents = 8891000; cross_section = 0.5104;}
	else if(process == "ttbarV_ttgamma"){NumberOfSimEvents = 5968000; cross_section = 0.5804;}
	else if(process == "ttbarV_ttgamma_ext"){NumberOfSimEvents = 4940000; cross_section = 0.5804;}
	else if(process == "ttbarV_ttHTobb"){NumberOfSimEvents = 9580000; cross_section = 0.5269;}
	else if(process == "ttbarV_ttHToNonbb"){NumberOfSimEvents = 7525991; cross_section = 0.5638;}
	else{cout << "Double check the input process you have entered" << endl;}	
	

  }
  else{cout << "Choose a year out of 2016, 2017 or 2018" << endl;}



  float norm_factor = (int_lumi * cross_section) / NumberOfSimEvents;
  return norm_factor;


}




void NormalisationFactorsCalculator2(const string& year){

  ofstream NormalisationFactors;
  string normfilename;

  normfilename = "NormalisationFactors_" + year + ".txt";
  NormalisationFactors.open(normfilename.c_str());

  vector<string> Samples = {"tZq", "ZPlusJets_M50_aMCatNLO", "ZPlusJets_M50_aMCatNLO_ext", "ZPlusJets_M50_Madgraph", "ZPlusJets_M50_Madgraph_ext", "ZPlusJets_M10To50_aMCatNLO", 
			    "ZPlusJets_M10To50_aMCatNLO_ext", "ZPlusJets_M10To50_Madgraph", "ZPlusJets_M10To50_Madgraph_ext", "ZPlusJets_PtBinned_0To50", 
			    "ZPlusJets_PtBinned_50To100", "ZPlusJets_PtBinned_50To100_ext", "ZPlusJets_PtBinned_100To250", "ZPlusJets_PtBinned_100To250_ext1", 
			    "ZPlusJets_PtBinned_100To250_ext2", "ZPlusJets_PtBinned_100To250_ext5", "ZPlusJets_PtBinned_250To400", "ZPlusJets_PtBinned_250To400_ext1", 
			    "ZPlusJets_PtBinned_250To400_ext2", "ZPlusJets_PtBinned_250To400_ext5", "ZPlusJets_PtBinned_400To650", "ZPlusJets_PtBinned_400To650_ext1", 
			    "ZPlusJets_PtBinned_400To650_ext2", "ZPlusJets_PtBinned_650ToInf", "ZPlusJets_PtBinned_650ToInf_ext1", "ZPlusJets_PtBinned_650ToInf_ext2", 
			    "ttbar_2l2nu", "ttbar_madgraph", "ttbar_madgraph_ext", "ttbar_TTToHadronic", "ttbar_TTToSemileptonic", "ttbar_atMCaNLO", "SingleTop_tchannel_top", 
			    "SingleTop_tchannel_top_ScaleUp", "SingleTop_tchannel_top_ScaleDown", "SingleTop_tchannel_antitop", "SingleTop_tchannel_antitop_ScaleUp", 
			    "SingleTop_tchannel_antitop_ScaleUp", "SingleTop_schannel", "ttbar_hdampUP", "ttbar_hdampUP_ext", "ttbar_hdampDOWN", "ttbar_hdampDOWN_ext", 
			    "SingleTop_tchannel_top_hdampUP", "SingleTop_tchannel_top_hdampDOWN", "ttbar_isr_UP", "ttbar_isr_DOWN", "ttbar_isr_DOWN_ext", "ttbar_fsr_UP", 
			    "ttbar_fsr_UP_ext", "ttbar_fsr_DOWN", "ttbar_fsr_DOWN_ext", "SingleTop_tW", "SingleTop_tW_ScaleUp", "SingleTop_tW_ScaleDown", "SingleTop_tbarW", 
			    "SingleTop_tbarW_ScaleUp", "SingleTop_tbarW_ScaleDown", "SingleTop_tHq", "SingleTop_tZq_W_lept_Z_had", "SingleTop_tWZ_tWll", "VV_ZZTo2l2nu", 
			    "VV_ZZTo2l2nu_ext", "VV_ZZTo2l2Q", "VV_ZZTo4L", "VV_WW1nuqq", "VV_WZTo2l2Q", "VV_WZTo3lNu", "VV_WZTo1l2Nu2Q", "VV_WWTo2l2Nu", "VV_WWToLNuQQ", 
		            "VV_WWToLNuQQ_ext", "VV_WGToLNuG", "VV_ZGToLLG", "VVV_WWWTo4F", "VVV_WWZTo4F", "VVV_WZZ", "VVV_ZZZ", "WPlusJets","WPlusJets_ext",
			    "ttbarV_ttWJetsToLNu", "ttbarV_ttWJetsToLNu_ext", "ttbarV_ttWJetsToQQ", "ttbarV_ttZToLL", "ttbarV_ttZToLL_ext2", "ttbarV_ttZToLL_ext3", 
			    "ttbarV_ttgamma", "ttbarV_ttgamma_ext", "ttbarV_ttHTobb", "ttbarV_ttHToNonbb"};


  for(int i = 0; i < Samples.size(); i++){

 	NormalisationFactors << Normalisation_Calculation(Samples.at(i), year) << endl;

  }

 cout << "The output file NormalisationFactors_" << year << ".txt has been written." << endl;


}



void NormalisationFactorsCalculator(){

 NormalisationFactorsCalculator2("2016");
 NormalisationFactorsCalculator2("2017");
 NormalisationFactorsCalculator2("2018");

}
