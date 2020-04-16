#include <iostream>
#include <fstream>

////This script won't compile but is just to copy and paste into main script (for RDataFrame column inputs)

auto EGammaTest2(const string& year, const string& type, const floats& pt, const floats& SuperclusterEta){

  floats OutputVector{};

  TFile* inputfile;


  for(int i = 0; i < pt.size(); i++){


  	if(year == "2016"){

  		if(type == "egammaEff" || type == "egammaEffSys"){

                	inputfile = new TFile("/home/eepgkkc/ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2016/egammaEffi_Tight_80X.txt_EGM2D.root", "READ");

        	}
        	else if(type == "egammaEffReco" || type == "egammaEffRecoSys"){

                	inputfile = new TFile("/home/eepgkkc/ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2016/egammaRecoEffi.txt_EGM2D.root", "READ");

        	}
		else{cout << "No EGamma SF input file found (2016)"<< endl;}

  	}
  	else if(year == "2017"){
 

  		if( (type == "egammaEffReco" || type == "egammaEffRecoSys") && pt.at(i) > 20){

    			inputfile = new TFile("/home/eepgkkc/ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2017/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root", "READ");

  		}
  		else if( (type == "egammaEffReco" || type == "egammaEffRecoSys") && pt.at(i) < 20){

    			inputfile = new TFile("/home/eepgkkc/ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2017/egammaEffi.txt_EGM2D_runBCDEF_passingRECO_lowEt.root", "READ");
  
  		}
  		else if(type == "egammaEff" || type == "egammaEffSys"){

    			inputfile = new TFile("/home/eepgkkc/ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2017/egammaEffi.txt_EGM2D_runBCDEF_passingTight94X.root", "READ");
  	
		}
		else{cout << "No EGamma SF input file found (2017)" << endl;}


  	}


  	TH2* histo = (TH2*)inputfile->GetObjectChecked("EGamma_SF2D", "TH2");



	if(abs(SuperClusterEta.at(i)) <= 2.5){

   		int SuperClusterEtaBin = histo->GetXaxis()->FindBin(SuperClusterEta.at(i));
        	int PtBin = histo->GetYaxis()->FindBin(pt.at(i));

		float EGammaSF;

		if(type == "egammaEffSys" || "egammaEffRecoSys"){
        		EGammaSF = histo->GetBinError(SuperClusterEtaBin, PtBin);
		}
		else{EGammaSF = histo->GetBinContent(SuperClusterEtaBin, PtBin);}

		OutputVector.push_back(EGammaSF); 
	}
	else{OutputVector.push_back(1.0);}


  }


  return OutputVector;

}





auto EGammaSF_egammaEff{[&year](){

 return EGammaTest2(year, "egammaEff");

}};



auto EGammaSF_egammaEffReco{[&year](){

 return EGammaTest2(year, "egammaEffReco"); 
 
}};


auto EGammaSF_egammaEff_Sys{[&year](){

 return EGammaTest2(year, "egammaEffSys");

}};



auto EGammaSF_egammaEffReco_Sys{[&year](){

 return EGammaTest2(year, "egammaEffRecoSys");

}};
