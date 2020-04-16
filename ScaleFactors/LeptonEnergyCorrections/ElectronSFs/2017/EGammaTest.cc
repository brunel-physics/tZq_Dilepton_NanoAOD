#include <iostream>
#include <fstream>



auto EGammaTest2(const string& year, const string& type){

  vector<float> OutputVector{};
  vector<float> SuperClusterEtaTest{0.125, -2.1, 1.9, 0.001, 0.5};
  vector<float> PtTest{10.1, 13.2, 19.1, 15.3, 250.0};

  TFile* inputfile;


  for(int i = 0; i < PtTest.size(); i++){


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
 

  		if( (type == "egammaEffReco" || type == "egammaEffRecoSys") && PtTest.at(i) > 20){

    			inputfile = new TFile("/home/eepgkkc/ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2017/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root", "READ");

  		}
  		else if( (type == "egammaEffReco" || type == "egammaEffRecoSys") && PtTest.at(i) < 20){

    			inputfile = new TFile("/home/eepgkkc/ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2017/egammaEffi.txt_EGM2D_runBCDEF_passingRECO_lowEt.root", "READ");
  
  		}
  		else if(type == "egammaEff" || type == "egammaEffSys"){

    			inputfile = new TFile("/home/eepgkkc/ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2017/egammaEffi.txt_EGM2D_runBCDEF_passingTight94X.root", "READ");
  	
		}
		else{cout << "No EGamma SF input file found (2017)" << endl;}


  	}


  	TH2* histo = (TH2*)inputfile->GetObjectChecked("EGamma_SF2D", "TH2");



	if(abs(SuperClusterEtaTest.at(i)) <= 2.5){

   		int SuperClusterEtaBin = histo->GetXaxis()->FindBin(SuperClusterEtaTest.at(i));
        	int PtBin = histo->GetYaxis()->FindBin(PtTest.at(i));
		
		float EGammaSF;

		if(type == "egammaEffSys" || type == "egammaEffRecoSys"){
        	 	EGammaSF = histo->GetBinError(SuperClusterEtaBin, PtBin);
		}
		else{EGammaSF = histo->GetBinContent(SuperClusterEtaBin, PtBin);}

		OutputVector.push_back(EGammaSF); 


	}
	else{OutputVector.push_back(1.0);}


  }


  return OutputVector;

}





auto EGammaTest(){

 //return EGammaTest2("2016", "egammaEff");
 //return EGammaTest2("2016", "egammaEffReco"); 
 //return EGammaTest2("2017", "egammaEff"); 
 //return EGammaTest2("2017", "egammaEffReco");
 

 //return EGammaTest2("2016", "egammaEffSys");
 //return EGammaTest2("2016", "egammaEffRecoSys"); 
 //return EGammaTest2("2017", "egammaEffSys"); 
 return EGammaTest2("2017", "egammaEffRecoSys"); 
 

}
