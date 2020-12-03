#include <iostream>
#include <fstream>

using namespace std;

auto Normalisation_Calculation(const string& process, const string& year){

  float int_lumi, cross_section, NumberOfSimEvents;

  if(year == "2016"){

	int_lumi = 35883;

	if(process == "tZq"){NumberOfSimEvents = 13932600; cross_section = 0.07358;}
	else if(process == "tZq_scaleup"){NumberOfSimEvents = 6940360; cross_section = 0.0758;}
	else if(process == "tZq_scaledown"){NumberOfSimEvents = 6982676; cross_section = 0.0758;}
	else if(process == "ZPlusJets_M50_aMCatNLO"){NumberOfSimEvents = 120777245; cross_section = 5941;}
	else if(process == "ZPlusJets_M50_aMCatNLO_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_M50_Madgraph"){NumberOfSimEvents = 49748967; cross_section = 4963;}
	else if(process == "ZPlusJets_M50_Madgraph_ext"){NumberOfSimEvents = 93007332; cross_section = 4963;}
	else if(process == "ZPlusJets_M10To50_aMCatNLO"){NumberOfSimEvents = 67981236; cross_section = 18810;}
	else if(process == "ZPlusJets_M10To50_aMCatNLO_ext"){NumberOfSimEvents = 40364234; cross_section = 18810;}
	else if(process == "ZPlusJets_M10To50_Madgraph"){NumberOfSimEvents = 34909242; cross_section = 16290;}
	else if(process == "ZPlusJets_M10To50_Madgraph_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_PtBinned_0To50"){NumberOfSimEvents = 37458375; cross_section = 5375;}
	else if(process == "ZPlusJets_PtBinned_50To100"){NumberOfSimEvents = 21847075; cross_section = 354.8;}
	else if(process == "ZPlusJets_PtBinned_50To100_ext"){NumberOfSimEvents = 108670239; cross_section = 354.8;}
	else if(process == "ZPlusJets_PtBinned_100To250"){NumberOfSimEvents = 2046961; cross_section = 81.22;}
	else if(process == "ZPlusJets_PtBinned_100To250_ext1"){NumberOfSimEvents = 2805972; cross_section = 81.22;}
	else if(process == "ZPlusJets_PtBinned_100To250_ext2"){NumberOfSimEvents = 2991815; cross_section = 81.22;}
	else if(process == "ZPlusJets_PtBinned_100To250_ext5"){NumberOfSimEvents = 76440229; cross_section = 81.22;}
	else if(process == "ZPlusJets_PtBinned_250To400"){NumberOfSimEvents = 423976; cross_section = 2.991;}
	else if(process == "ZPlusJets_PtBinned_250To400_ext1"){NumberOfSimEvents = 590806; cross_section = 2.991;}
	else if(process == "ZPlusJets_PtBinned_250To400_ext2"){NumberOfSimEvents = 594317; cross_section = 2.991;}
        else if(process == "ZPlusJets_PtBinned_250To400_ext5"){NumberOfSimEvents = 19567800; cross_section = 2.991;}
        else if(process == "ZPlusJets_PtBinned_400To650"){NumberOfSimEvents = 432056; cross_section = 0.3882;}
        else if(process == "ZPlusJets_PtBinned_400To650_ext1"){NumberOfSimEvents = 589842; cross_section = 0.3882;}
        else if(process == "ZPlusJets_PtBinned_400To650_ext2"){NumberOfSimEvents = 604038; cross_section = 0.3882;}
        else if(process == "ZPlusJets_PtBinned_650ToInf"){NumberOfSimEvents = 430691; cross_section = 0.03737;}
        else if(process == "ZPlusJets_PtBinned_650ToInf_ext1"){NumberOfSimEvents = 599665; cross_section = 0.03737;}
        else if(process == "ZPlusJets_PtBinned_650ToInf_ext2"){NumberOfSimEvents = 597526; cross_section = 0.03737;}
        else if(process == "ttbar_2l2nu"){NumberOfSimEvents = 67860400; cross_section = 687.1;}
        else if(process == "ttbar_madgraph"){NumberOfSimEvents = 6068369; cross_section = 56.86;}
        else if(process == "ttbar_madgraph_ext"){NumberOfSimEvents = 24767666; cross_section = 56.86;}
        else if(process == "ttbar_TTToHadronic"){NumberOfSimEvents = 68302800; cross_section = 687.1;}
        else if(process == "ttbar_TTToHadronic_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbar_TTToSemileptonic"){NumberOfSimEvents = 107604800; cross_section = 687.1;}
        else if(process == "ttbar_TTToSemileptonic_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbar_atMCaNLO"){NumberOfSimEvents = 43556919; cross_section = 722.8;}
        else if(process == "ttbar_atMCaNLO_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbar_inc"){NumberOfSimEvents = 76367863; cross_section = 730.6;}
        else if(process == "SingleTop_tchannel_top"){NumberOfSimEvents = 58403420; cross_section = 136.02;}
	else if(process == "SingleTop_tchannel_top_ScaleUp"){NumberOfSimEvents = 5992440; cross_section = 136.02;}
        else if(process == "SingleTop_tchannel_top_ScaleDown"){NumberOfSimEvents = 65028288; cross_section = 136.02;}
        else if(process == "SingleTop_tchannel_antitop"){NumberOfSimEvents = 38811017; cross_section = 80.95;}
        else if(process == "SingleTop_tchannel_antitop_ScaleUp"){NumberOfSimEvents = 3970546; cross_section = 80.95;}
        else if(process == "SingleTop_tchannel_antitop_ScaleDown"){NumberOfSimEvents = 35658214; cross_section = 80.95;}
        else if(process == "SingleTop_schannel"){NumberOfSimEvents = 2917199; cross_section = 10.12;}
        else if(process == "ttbar_hdampUP"){NumberOfSimEvents = 29833668; cross_section = 730.6;}
        else if(process == "ttbar_hdampUP_ext"){NumberOfSimEvents = 28855428; cross_section = 730.6;}
        else if(process == "ttbar_hdampDOWN"){NumberOfSimEvents = 29124629; cross_section = 730.7;}
        else if(process == "ttbar_hdampDOWN_ext"){NumberOfSimEvents = 29229088; cross_section = 730.7;}
        else if(process == "TT_2l2nu_hdampUP"){NumberOfSimEvents = 9923800; cross_section = 687;}
        else if(process == "TT_2l2nu_hdampUP_ext1"){NumberOfSimEvents = 4965300; cross_section = 687;}
        else if(process == "TT_2l2nu_hdampUP_ext2"){NumberOfSimEvents = 29860800; cross_section = 687;}
        else if(process == "TT_2l2nu_hdampDOWN"){NumberOfSimEvents = 9963900; cross_section = 687.1;}
        else if(process == "TT_2l2nu_hdampDOWN_ext1"){NumberOfSimEvents = 4944800; cross_section = 687.1;}
        else if(process == "TT_2l2nu_hdampDOWN_ext2"){NumberOfSimEvents = 29190500; cross_section = 687.1;}
        else if(process == "TTToHadronic_hdampUP"){NumberOfSimEvents = 28695100; cross_section = 687;}
        else if(process == "TTToHadronic_hdampDOWN"){NumberOfSimEvents = 28900700; cross_section = 687.1;}
        else if(process == "TTToSemileptonic_hdampUP"){NumberOfSimEvents = 29671200; cross_section = 687;}
	else if(process == "TTToSemileptonic_hdampDOWN"){NumberOfSimEvents = 29818400; cross_section = 687.1;}
        else if(process == "TTToSemileptonic_hdampDOWN_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "SingleTop_tchannel_top_hdampUP"){NumberOfSimEvents = 5998200; cross_section = 113.3;}
        else if(process == "SingleTop_tchannel_top_hdampDOWN"){NumberOfSimEvents = 5922400; cross_section = 113.3;}
        else if(process == "SingleTop_tchannel_antitop_hdampUP"){NumberOfSimEvents = 3999400; cross_section = 67.91;}
        else if(process == "SingleTop_tchannel_antitop_hdampDOWN"){NumberOfSimEvents = 3999400; cross_section = 67.91;}
        else if(process == "ttbar_isr_UP"){NumberOfSimEvents = 58977100; cross_section = 630.6;}
        else if(process == "ttbar_isr_DOWN"){NumberOfSimEvents = 28504600; cross_section = 730.6;}
        else if(process == "ttbar_isr_DOWN_ext"){NumberOfSimEvents = 29915551; cross_section = 730.6;}
        else if(process == "ttbar_fsr_UP"){NumberOfSimEvents = 29632372; cross_section = 730.6;}
        else if(process == "ttbar_fsr_UP_ext"){NumberOfSimEvents = 29501065; cross_section = 730.6;}
        else if(process == "ttbar_fsr_DOWN"){NumberOfSimEvents = 29571600; cross_section = 730.6;}
        else if(process == "ttbar_fsr_DOWN_ext"){NumberOfSimEvents = 29602576; cross_section = 730.6;}
        else if(process == "SingleTop_tW"){NumberOfSimEvents = 4983500; cross_section = 34.91;}
        else if(process == "SingleTop_tW_ScaleUp"){NumberOfSimEvents = 997880; cross_section = 38.09;}
        else if(process == "SingleTop_tW_ScaleDown"){NumberOfSimEvents = 993640; cross_section = 38.09;}
        else if(process == "SingleTop_tbarW"){NumberOfSimEvents = 6933094; cross_section = 38.06;}
        else if(process == "SingleTop_tbarW_ScaleUp"){NumberOfSimEvents = 1000000; cross_section = 38.06;}
        else if(process == "SingleTop_tbarW_ScaleDown"){NumberOfSimEvents = 999068; cross_section = 38.06;}
	else if(process == "SingleTop_tHq"){NumberOfSimEvents = 3495799; cross_section = 0.2609;}
        else if(process == "SingleTop_tZq_W_lept_Z_had"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "SingleTop_tWZ_tWll"){NumberOfSimEvents = 50000; cross_section = 0.01104;}
        else if(process == "VV_ZZTo2l2nu"){NumberOfSimEvents = 8931750; cross_section = 0.5644;}
        else if(process == "VV_ZZTo2l2nu_ext1"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "VV_ZZTo2l2nu_ext2"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "VV_ZZTo2l2Q"){NumberOfSimEvents = 15462693; cross_section = 3.688;}
        else if(process == "VV_ZZTo4L"){NumberOfSimEvents = 10711278; cross_section = 1.204;}
        else if(process == "VV_ZZTo4L_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "VV_WZTo2l2Q"){NumberOfSimEvents = 26445785; cross_section = 5.606;}
        else if(process == "VV_WZTo3lNu"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "VV_WZTo3lNu_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "VV_WZTo1l1Nu2Q"){NumberOfSimEvents = 24311445; cross_section = 11.66;}
        else if(process == "VV_WWTo2l2Nu"){NumberOfSimEvents = 1999000; cross_section = 10.48;}
        else if(process == "VV_WWToLNuQQ"){NumberOfSimEvents = 1999200; cross_section = 43.53;}
        else if(process == "VV_WWToLNuQQ_ext"){NumberOfSimEvents = 6655400; cross_section = 43.53;}
        else if(process == "VV_WGToLNuG"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "VV_ZGToLLG"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "VV_ZGToLLG_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "VVV_WWWTo4F"){NumberOfSimEvents = 240000; cross_section = 0.2086;}
        else if(process == "VVV_WWWTo4F_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "VVV_WWZTo4F"){NumberOfSimEvents = 250000; cross_section = 0.1651;}
        else if(process == "VVV_WWZTo4F_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "VVV_WZZ"){NumberOfSimEvents = 246800; cross_section = 0.05565;}
        else if(process == "VVV_WZZ_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "VVV_ZZZ"){NumberOfSimEvents = 249237; cross_section = 0.01398;}
        else if(process == "VVV_ZZZ_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "WPlusJets"){NumberOfSimEvents = 22864149; cross_section = 60430;}
        else if(process == "WPlusJets_ext"){NumberOfSimEvents = 237263153; cross_section = 60430;}
        else if(process == "ttbarV_ttWJetsToLNu"){NumberOfSimEvents = 2160168; cross_section = 0.2001;}
        else if(process == "ttbarV_ttWJetsToLNu_ext"){NumberOfSimEvents = 3120397; cross_section = 0.2001;}
        else if(process == "ttbarV_ttWJetsToQQ"){NumberOfSimEvents = 833298; cross_section = 0.405;}
        else if(process == "ttbarV_ttZToLL"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbarV_ttZToLL_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbarV_ttZToLLNuNu"){NumberOfSimEvents = 1992438; cross_section = 0.2529;}
        else if(process == "ttbarV_ttZToLLNuNu_ext"){NumberOfSimEvents = 5837781; cross_section = 0.2529;}
        else if(process == "ttbarV_ttZToLLNuNu_ext2"){NumberOfSimEvents = 5862409; cross_section = 0.2529;}
        else if(process == "ttbarV_ttZToQQ"){NumberOfSimEvents = 749400; cross_section = 0.5297;}
	else if(process == "ttbarV_ttZToQQ_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbarV_ttHTobb"){NumberOfSimEvents = 3872944; cross_section = 0.5638;}
        else if(process == "ttbarV_ttHTobb_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbarV_ttHToNonbb"){NumberOfSimEvents = 3981250; cross_section = 0.5638;}
	else{throw std::logic_error("ERROR: Process name not recognised");}

  }
  else if(year == "2017"){

	int_lumi = 41528;

	if(process == "tZq"){NumberOfSimEvents = 13276146; cross_section = 0.07358;}
	else if(process == "tZq_scaleup"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "tZq_scaledown"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_M50_aMCatNLO"){NumberOfSimEvents = 27529915; cross_section = 6529;}
	else if(process == "ZPlusJets_M50_aMCatNLO_ext"){NumberOfSimEvents = 180273927; cross_section = 6529;}
	else if(process == "ZPlusJets_M50_Madgraph"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_M50_Madgraph_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_M10To50_aMCatNLO"){NumberOfSimEvents = 0; cross_section = 0;}
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
        else if(process == "ttbar_2l2nu"){NumberOfSimEvents = 68812194; cross_section = 88.29;}
        else if(process == "ttbar_madgraph"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbar_madgraph_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbar_TTToHadronic"){NumberOfSimEvents = 128829950; cross_section = 377.96;}
        else if(process == "ttbar_TTToHadronic_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbar_TTToSemileptonic"){NumberOfSimEvents = 109663792; cross_section = 365.34;}
        else if(process == "ttbar_TTToSemileptonic_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbar_atMCaNLO"){NumberOfSimEvents = 154280331; cross_section = 722.8;}
        else if(process == "ttbar_atMCaNLO_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbar_inc"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "SingleTop_tchannel_top"){NumberOfSimEvents = 122688200; cross_section = 113.3;}
	else if(process == "SingleTop_tchannel_top_ScaleUp"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "SingleTop_tchannel_top_ScaleDown"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "SingleTop_tchannel_antitop"){NumberOfSimEvents = 64761200; cross_section = 67.91;}
        else if(process == "SingleTop_tchannel_antitop_ScaleUp"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "SingleTop_tchannel_antitop_ScaleDown"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "SingleTop_schannel"){NumberOfSimEvents = 9883805; cross_section = 3.74;}
        else if(process == "ttbar_hdampUP"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbar_hdampUP_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbar_hdampDOWN"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbar_hdampDOWN_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "TT_2l2nu_hdampUP"){NumberOfSimEvents = 3288128; cross_section = 88.29;}
        else if(process == "TT_2l2nu_hdampUP_ext1"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "TT_2l2nu_hdampUP_ext2"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "TT_2l2nu_hdampDOWN"){NumberOfSimEvents = 5476459; cross_section = 88.29;}
        else if(process == "TT_2l2nu_hdampDOWN_ext1"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "TT_2l2nu_hdampDOWN_ext2"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "TTToHadronic_hdampUP"){NumberOfSimEvents = 27203920; cross_section = 377.96;}
        else if(process == "TTToHadronic_hdampDOWN"){NumberOfSimEvents = 27117982; cross_section = 377.96;}
        else if(process == "TTToSemileptonic_hdampUP"){NumberOfSimEvents = 23977012; cross_section = 365.34;}
	else if(process == "TTToSemileptonic_hdampDOWN"){NumberOfSimEvents = 26367765; cross_section = 365.34;}
        else if(process == "TTToSemileptonic_hdampDOWN_ext"){NumberOfSimEvents = 26849863; cross_section = 365.34;}
        else if(process == "SingleTop_tchannel_top_hdampUP"){NumberOfSimEvents = 1990000; cross_section = 113.3;}
        else if(process == "SingleTop_tchannel_top_hdampDOWN"){NumberOfSimEvents = 5965095; cross_section = 113.3;}
        else if(process == "SingleTop_tchannel_antitop_hdampUP"){NumberOfSimEvents = 976576; cross_section = 67.91;}
        else if(process == "SingleTop_tchannel_antitop_hdampDOWN"){NumberOfSimEvents = 992608; cross_section = 67.91;}
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
        else if(process == "SingleTop_tbarW"){NumberOfSimEvents = 7374736; cross_section = 34.97;}
        else if(process == "SingleTop_tbarW_ScaleUp"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "SingleTop_tbarW_ScaleDown"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "SingleTop_tHq"){NumberOfSimEvents = 9918994; cross_section = 0.729;}
        else if(process == "SingleTop_tZq_W_lept_Z_had"){NumberOfSimEvents = 1000000; cross_section = 0.1573;}
        else if(process == "SingleTop_tWZ_tWll"){NumberOfSimEvents = 986000; cross_section = 0.01103;}
        else if(process == "VV_ZZTo2l2nu"){NumberOfSimEvents = 8744768; cross_section = 0.5644;}
        else if(process == "VV_ZZTo2l2nu_ext1"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "VV_ZZTo2l2nu_ext2"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "VV_ZZTo2l2Q"){NumberOfSimEvents = 27441872; cross_section = 3.222;}
        else if(process == "VV_ZZTo4L"){NumberOfSimEvents = 6964071; cross_section = 1.325;}
        else if(process == "VV_ZZTo4L_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "VV_WZTo2l2Q"){NumberOfSimEvents = 27582164; cross_section = 5.606;}
        else if(process == "VV_WZTo3lNu"){NumberOfSimEvents = 10987679; cross_section = 5.052;}
        else if(process == "VV_WZTo3lNu_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "VV_WZTo1l1Nu2Q"){NumberOfSimEvents = 19086373; cross_section = 10.73;}
        else if(process == "VV_WWTo2l2Nu"){NumberOfSimEvents = 2000000; cross_section = 11.08;}
        else if(process == "VV_WWToLNuQQ"){NumberOfSimEvents = 8785360; cross_section = 45.99;}
        else if(process == "VV_WWToLNuQQ_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "VV_WGToLNuG"){NumberOfSimEvents = 6283083; cross_section = 463.9;}
        else if(process == "VV_ZGToLLG"){NumberOfSimEvents = 30490034; cross_section = 55.48;}
        else if(process == "VV_ZGToLLG_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "VVV_WWWTo4F"){NumberOfSimEvents = 240000; cross_section = 0.2086;}
        else if(process == "VVV_WWWTo4F_ext"){NumberOfSimEvents = 232300; cross_section = 0.2086;}
        else if(process == "VVV_WWZTo4F"){NumberOfSimEvents = 250000; cross_section = 0.1651;}
        else if(process == "VVV_WWZTo4F_ext"){NumberOfSimEvents = 250000; cross_section = 0.1651;}
        else if(process == "VVV_WZZ"){NumberOfSimEvents = 250000; cross_section = 0.05565;}
        else if(process == "VVV_WZZ_ext"){NumberOfSimEvents = 250000; cross_section = 0.05565;}
        else if(process == "VVV_ZZZ"){NumberOfSimEvents = 250000; cross_section = 0.01398;}
        else if(process == "VVV_ZZZ_ext"){NumberOfSimEvents = 250000; cross_section = 0.01398;}
        else if(process == "WPlusJets"){NumberOfSimEvents = 33073306; cross_section = 52940;}
        else if(process == "WPlusJets_ext"){NumberOfSimEvents = 30008250; cross_section = 52940;}
        else if(process == "ttbarV_ttWJetsToLNu"){NumberOfSimEvents = 4908905; cross_section = 0.2198;}
        else if(process == "ttbarV_ttWJetsToLNu_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbarV_ttWJetsToQQ"){NumberOfSimEvents = 811306; cross_section = 0.4316;}
        else if(process == "ttbarV_ttZToLL"){NumberOfSimEvents = 250000; cross_section = 0.05324;}
        else if(process == "ttbarV_ttZToLL_ext"){NumberOfSimEvents = 245978; cross_section = 0.05324;}
        else if(process == "ttbarV_ttZToLLNuNu"){NumberOfSimEvents = 7563490; cross_section = 0.2432;}
        else if(process == "ttbarV_ttZToLLNuNu_ext"){NumberOfSimEvents = 7932650; cross_section = 0.2432;}
        else if(process == "ttbarV_ttZToLLNuNu_ext2"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbarV_ttZToQQ"){NumberOfSimEvents = 750000; cross_section = 0.5104;}
	else if(process == "ttbarV_ttZToQQ_ext"){NumberOfSimEvents = 8940000; cross_section = 0.5104;}
        else if(process == "ttbarV_ttHTobb"){NumberOfSimEvents = 8000000; cross_section = 0.5269;}
        else if(process == "ttbarV_ttHTobb_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbarV_ttHToNonbb"){NumberOfSimEvents = 5499293; cross_section = 0.5638;}
	else{throw std::logic_error("ERROR: Process name not recognised");}	
	

  }
  else if(year == "2018"){

	int_lumi = 59688;

	if(process == "tZq"){NumberOfSimEvents = 13736000; cross_section = 0.07358;}
	else if(process == "tZq_scaleup"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "tZq_scaledown"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_M50_aMCatNLO"){NumberOfSimEvents = 997561; cross_section = 6529;}
	else if(process == "ZPlusJets_M50_aMCatNLO_ext"){NumberOfSimEvents = 193119590; cross_section = 6529;}
	else if(process == "ZPlusJets_M50_Madgraph"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_M50_Madgraph_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_M10To50_aMCatNLO"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_M10To50_aMCatNLO_ext"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "ZPlusJets_M10To50_Madgraph"){NumberOfSimEvents = 39392062; cross_section = 15810;}
	else if(process == "ZPlusJets_M10To50_Madgraph_ext"){NumberOfSimEvents = 46976952; cross_section = 15810;}
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
        else if(process == "ttbar_madgraph"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbar_madgraph_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbar_TTToHadronic"){NumberOfSimEvents = 128640000; cross_section = 377.96;}
        else if(process == "ttbar_TTToHadronic_ext"){NumberOfSimEvents = 198756000; cross_section = 377.96;}
        else if(process == "ttbar_TTToSemileptonic"){NumberOfSimEvents = 100790000; cross_section = 365.34;}
        else if(process == "ttbar_TTToSemileptonic_ext"){NumberOfSimEvents = 199829998; cross_section = 365.34;}
        else if(process == "ttbar_atMCaNLO"){NumberOfSimEvents = 142155064; cross_section = 831.76;}
        else if(process == "ttbar_atMCaNLO_ext"){NumberOfSimEvents = 162376405; cross_section = 831.76;}
        else if(process == "ttbar_inc"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "SingleTop_tchannel_top"){NumberOfSimEvents = 154307600; cross_section = 115.3;}
	else if(process == "SingleTop_tchannel_top_ScaleUp"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "SingleTop_tchannel_top_ScaleDown"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "SingleTop_tchannel_antitop"){NumberOfSimEvents = 79090800; cross_section = 69.09;}
        else if(process == "SingleTop_tchannel_antitop_ScaleUp"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "SingleTop_tchannel_antitop_ScaleDown"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "SingleTop_schannel"){NumberOfSimEvents = 19890000; cross_section = 3.74;}
        else if(process == "ttbar_hdampUP"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbar_hdampUP_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbar_hdampDOWN"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbar_hdampDOWN_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "TT_2l2nu_hdampUP"){NumberOfSimEvents = 5260000; cross_section = 687;}
        else if(process == "TT_2l2nu_hdampUP_ext1"){NumberOfSimEvents = 9903000; cross_section = 687;}
        else if(process == "TT_2l2nu_hdampUP_ext2"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "TT_2l2nu_hdampDOWN"){NumberOfSimEvents = 5458000; cross_section = 687.1;}
        else if(process == "TT_2l2nu_hdampDOWN_ext1"){NumberOfSimEvents = 9952000; cross_section = 687.1;}
        else if(process == "TT_2l2nu_hdampDOWN_ext2"){NumberOfSimEvents = 29312000; cross_section = 687.1;}
        else if(process == "TTToHadronic_hdampUP"){NumberOfSimEvents = 24965000; cross_section = 687;}
        else if(process == "TTToHadronic_hdampDOWN"){NumberOfSimEvents = 26425000; cross_section = 687.1;}
        else if(process == "TTToSemileptonic_hdampUP"){NumberOfSimEvents = 26964000; cross_section = 687;}
	else if(process == "TTToSemileptonic_hdampDOWN"){NumberOfSimEvents = 25904000; cross_section = 687.1;}
        else if(process == "TTToSemileptonic_hdampDOWN_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "SingleTop_tchannel_top_hdampUP"){NumberOfSimEvents = 1913789; cross_section = 115.3;}
        else if(process == "SingleTop_tchannel_top_hdampDOWN"){NumberOfSimEvents = 5855752; cross_section = 115.3;}
        else if(process == "SingleTop_tchannel_antitop_hdampUP"){NumberOfSimEvents = 997000; cross_section = 69.09;}
        else if(process == "SingleTop_tchannel_antitop_hdampDOWN"){NumberOfSimEvents = 981000; cross_section = 69.09;}
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
        else if(process == "SingleTop_tbarW"){NumberOfSimEvents = 7527000; cross_section = 34.97;}
        else if(process == "SingleTop_tbarW_ScaleUp"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "SingleTop_tbarW_ScaleDown"){NumberOfSimEvents = 0; cross_section = 0;}
	else if(process == "SingleTop_tHq"){NumberOfSimEvents = 29683974; cross_section = 0.729;}
        else if(process == "SingleTop_tZq_W_lept_Z_had"){NumberOfSimEvents = 4977000; cross_section = 0.1518;}
        else if(process == "SingleTop_tWZ_tWll"){NumberOfSimEvents = 248600; cross_section = 0.01103;}
        else if(process == "VV_ZZTo2l2nu"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "VV_ZZTo2l2nu_ext1"){NumberOfSimEvents = 8382600; cross_section = 0.6008;}
        else if(process == "VV_ZZTo2l2nu_ext2"){NumberOfSimEvents = 48046000; cross_section = 0.6008;}
        else if(process == "VV_ZZTo2l2Q"){NumberOfSimEvents = 27900469; cross_section = 3.222;}
        else if(process == "VV_ZZTo4L"){NumberOfSimEvents = 6689900; cross_section = 1.325;}
        else if(process == "VV_ZZTo4L_ext"){NumberOfSimEvents = 98870000; cross_section = 1.325;}
        else if(process == "VV_WZTo2l2Q"){NumberOfSimEvents = 28193648; cross_section = 5.606;}
        else if(process == "VV_WZTo3lNu"){NumberOfSimEvents = 10749269; cross_section = 5.052;}
        else if(process == "VV_WZTo3lNu_ext"){NumberOfSimEvents = 11248318; cross_section = 5.052;}
        else if(process == "VV_WZTo1l1Nu2Q"){NumberOfSimEvents = 18273533; cross_section = 10.73;}
        else if(process == "VV_WWTo2l2Nu"){NumberOfSimEvents = 7758900; cross_section = 11.08;}
        else if(process == "VV_WWToLNuQQ"){NumberOfSimEvents = 19199100; cross_section = 45.99;}
        else if(process == "VV_WWToLNuQQ_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "VV_WGToLNuG"){NumberOfSimEvents = 6108186; cross_section = 463.9;}
        else if(process == "VV_ZGToLLG"){NumberOfSimEvents = 13946364; cross_section = 55.48;}
        else if(process == "VV_ZGToLLG_ext"){NumberOfSimEvents = 19767732; cross_section = 55.48;}
	else if(process == "VVV_WWWTo4F"){NumberOfSimEvents = 240000; cross_section = 0.2086;}
        else if(process == "VVV_WWWTo4F_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "VVV_WWZTo4F"){NumberOfSimEvents = 250000; cross_section = 0.1676;}
        else if(process == "VVV_WWZTo4F_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "VVV_WZZ"){NumberOfSimEvents = 250000; cross_section = 0.05565;}
        else if(process == "VVV_WZZ_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "VVV_ZZZ"){NumberOfSimEvents = 250000; cross_section = 0.01398;}
        else if(process == "VVV_ZZZ_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "WPlusJets"){NumberOfSimEvents = 71026861; cross_section = 52940;}
        else if(process == "WPlusJets_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbarV_ttWJetsToLNu"){NumberOfSimEvents = 4911941; cross_section = 0.2149;}
        else if(process == "ttbarV_ttWJetsToLNu_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbarV_ttWJetsToQQ"){NumberOfSimEvents = 835296; cross_section = 0.4316;}
        else if(process == "ttbarV_ttZToLL"){NumberOfSimEvents = 250000; cross_section = 0.05324;}
        else if(process == "ttbarV_ttZToLL_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbarV_ttZToLLNuNu"){NumberOfSimEvents = 13280000; cross_section = 0.2432;}
        else if(process == "ttbarV_ttZToLLNuNu_ext"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbarV_ttZToLLNuNu_ext2"){NumberOfSimEvents = 0; cross_section = 0;}
        else if(process == "ttbarV_ttZToQQ"){NumberOfSimEvents = 750000; cross_section = 0.5104;}
	else if(process == "ttbarV_ttZToQQ_ext"){NumberOfSimEvents = 8822000; cross_section = 0.5104;}
        else if(process == "ttbarV_ttHTobb"){NumberOfSimEvents = 11835999; cross_section = 0.5269;}
        else if(process == "ttbarV_ttHTobb_ext"){NumberOfSimEvents = 9484000; cross_section = 0.5269;}
        else if(process == "ttbarV_ttHToNonbb"){NumberOfSimEvents = 7525991; cross_section = 0.5638;}
	else{throw std::logic_error("ERROR: Process name not recognised");}	

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

  vector<string> Samples = {"tZq", 			            "tZq_scaleup", 		              "tZq_scaledown",
			    "ZPlusJets_M50_aMCatNLO",               "ZPlusJets_M50_aMCatNLO_ext",      	      "ZPlusJets_M50_Madgraph",
			    "ZPlusJets_M50_Madgraph_ext",           "ZPlusJets_M10To50_aMCatNLO",             "ZPlusJets_M10To50_aMCatNLO_ext",
			    "ZPlusJets_M10To50_Madgraph",           "ZPlusJets_M10To50_Madgraph_ext",         "ZPlusJets_PtBinned_0To50",
			    "ZPlusJets_PtBinned_50To100", 	    "ZPlusJets_PtBinned_50To100_ext",         "ZPlusJets_PtBinned_100To250",
			    "ZPlusJets_PtBinned_100To250_ext1",     "ZPlusJets_PtBinned_100To250_ext2",       "ZPlusJets_PtBinned_100To250_ext5",
		            "ZPlusJets_PtBinned_250To400",          "ZPlusJets_PtBinned_250To400_ext1",       "ZPlusJets_PtBinned_250To400_ext2",
			    "ZPlusJets_PtBinned_250To400_ext5",     "ZPlusJets_PtBinned_400To650",            "ZPlusJets_PtBinned_400To650_ext1",
			    "ZPlusJets_PtBinned_400To650_ext2",     "ZPlusJets_PtBinned_650ToInf",            "ZPlusJets_PtBinned_650ToInf_ext1",
			    "ZPlusJets_PtBinned_650ToInf_ext2",     "ttbar_2l2nu",			      "ttbar_madgraph",
			    "ttbar_madgraph_ext",		    "ttbar_TTToHadronic",		      "ttbar_TTToHadronic_ext",
			    "ttbar_TTToSemileptonic",		    "ttbar_TTToSemileptonic_ext",	      "ttbar_atMCaNLO",
			    "ttbar_atMCaNLO_ext",		    "ttbar_inc",			      "SingleTop_tchannel_top",
		            "SingleTop_tchannel_top_ScaleUp",       "SingleTop_tchannel_top_ScaleDown",       "SingleTop_tchannel_antitop",
			    "SingleTop_tchannel_antitop_ScaleUp",   "SingleTop_tchannel_antitop_ScaleDown",   "SingleTop_schannel",
			    "ttbar_hdampUP",			    "ttbar_hdampUP_ext",		      "ttbar_hdampDOWN",
			    "ttbar_hdampDOWN_ext",		    "TT_2l2nu_hdampUP",			      "TT_2l2nu_hdampUP_ext1",
			    "TT_2l2nu_hdampUP_ext2",		    "TT_2l2nu_hdampDOWN",		      "TT_2l2nu_hdampDOWN_ext1",
			    "TT_2l2nu_hdampDOWN_ext2",		    "TTToHadronic_hdampUP",	              "TTToHadronic_hdampDOWN",
			    "TTToSemileptonic_hdampUP",		    "TTToSemileptonic_hdampDOWN",	      "TTToSemileptonic_hdampDOWN_ext",
	                    "SingleTop_tchannel_top_hdampUP",	    "SingleTop_tchannel_top_hdampDOWN",	      "SingleTop_tchannel_antitop_hdampUP",
			    "SingleTop_tchannel_antitop_hdampDOWN", "ttbar_isr_UP",			      "ttbar_isr_DOWN",
			    "ttbar_isr_DOWN_ext",		    "ttbar_fsr_UP",			      "ttbar_fsr_UP_ext",
			    "ttbar_fsr_DOWN",			    "ttbar_fsr_DOWN_ext",		      "SingleTop_tW",
			    "SingleTop_tW_ScaleUp",		    "SingleTop_tW_ScaleDown",		      "SingleTop_tbarW",
			    "SingleTop_tbarW_ScaleUp",		    "SingleTop_tbarW_ScaleDown",	      "SingleTop_tHq",
		            "SingleTop_tZq_W_lept_Z_had",	    "SingleTop_tWZ_tWll",		      "VV_ZZTo2l2nu",
			    "VV_ZZTo2l2nu_ext1",		    "VV_ZZTo2l2nu_ext2",		      "VV_ZZTo2l2Q",
		            "VV_ZZTo4L",			    "VV_ZZTo4L_ext",			      "VV_WZTo2l2Q",
			    "VV_WZTo3lNu",			    "VV_WZTo3lNu_ext",			      "VV_WZTo1l1Nu2Q",
			    "VV_WWTo2l2Nu",			    "VV_WWToLNuQQ",			      "VV_WWToLNuQQ_ext",
			    "VV_WGToLNuG",			    "VV_ZGToLLG",			      "VV_ZGToLLG_ext",
			    "VVV_WWWTo4F",			    "VVV_WWWTo4F_ext",			      "VVV_WWZTo4F",
			    "VVV_WWZTo4F_ext",			    "VVV_WZZ",				      "VVV_WZZ_ext",
			    "VVV_ZZZ",				    "VVV_ZZZ_ext",			      "WPlusJets",
			    "WPlusJets_ext",			    "ttbarV_ttWJetsToLNu",		      "ttbarV_ttWJetsToLNu_ext",
			    "ttbarV_ttWJetsToQQ",		    "ttbarV_ttZToLL",			      "ttbarV_ttZToLL_ext",
			    "ttbarV_ttZToLLNuNu",		    "ttbarV_ttZToLLNuNu_ext",		      "ttbarV_ttZToLLNuNu_ext2",
			    "ttbarV_ttZToQQ",			    "ttbarV_ttZToQQ_ext",		      "ttbarV_ttHTobb",
			    "ttbarV_ttHTobb_ext",		    "ttbarV_ttHToNonbb"}; 


  for(int i = 0; i < Samples.size(); i++){ NormalisationFactors << Normalisation_Calculation(Samples.at(i), year) << endl; }

  cout << "The output file NormalisationFactors_" << year << ".txt has been written." << endl; 


}



void NormalisationFactorsCalculator(){

 NormalisationFactorsCalculator2("2016");
 NormalisationFactorsCalculator2("2017");
 NormalisationFactorsCalculator2("2018");

}
