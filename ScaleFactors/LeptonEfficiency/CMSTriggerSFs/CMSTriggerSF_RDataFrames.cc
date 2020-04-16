#include <iostream>
#include <fstream>


auto CMSTriggerSF_RDataFrames2{[](const string& year, const string& UpOrDown, const floats& pts, const floats& etas){

  floats AbsEta = abs(etas);

  float lumiRunBCDEF_triggerSF = 19648.534;
  float lumiRunGH_triggerSF = 16144.444;

  TFile* inputfile_RunsBCDEF;
  TFile* inputfile_RunsGH;
  TH2* histo_RunsBCDEF;
  TH2* histo_RunsGH;


  if(year == "2016"){

  	//HLT Mu24 file (runs BCDEF)
	inputfile_RunsBCDEF = new TFile("/home/eepgkkc/ScaleFactors/LeptonEfficiency/CMSTriggerSFs/2016/HLT_Mu24_EfficienciesAndSF_RunBtoF.root", "READ");
	histo_RunsBCDEF = (TH2*)inputfile_RunsBCDEF->GetObjectChecked("IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio", "TH2");

        //HLT Mu24 ID file (runs GH)
        inputfile_RunsGH = new TFile("/home/eepgkkc/ScaleFactors/LeptonEfficiency/CMSTriggerSFs/2016/HLT_Mu24_EfficienciesAndSF_RunGtoH.root", "READ");
        histo_RunsGH = (TH2*)inputfile_RunsGH->GetObjectChecked("IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio", "TH2");

  }
  else if(year == "2017"){

    inputfile_RunsBCDEF = new TFile("/home/eepgkkc/ScaleFactors/LeptonEfficiency/CMSTriggerSFs/2017/HLT_Mu24_EfficienciesAndSF_RunBtoF.root", "READ");
    histo_RunsBCDEF = (TH2*)inputfile_RunsBCDEF->GetObjectChecked("IsoMu27_PtEtaBins/pt_abseta_ratio", "TH2");  

  }
  else{cout << "need to add code for 2018" << endl;}






  vector<float> CMSTriggerSFOutput{};
 
  for(int i = 0; i < pts.size(); i++){


  		int PtBin_RunsBCDEF = histo_RunsBCDEF->GetXaxis()->FindBin(pts.at(i));
        	int AbsEtaBin_RunsBCDEF = histo_RunsBCDEF->GetYaxis()->FindBin(AbsEta.at(i));

        	float CMSTriggerSF_RunsBCDEF = histo_RunsBCDEF->GetBinContent(PtBin_RunsBCDEF, AbsEtaBin_RunsBCDEF);
		float Error_RunsBCDEF = histo_RunsBCDEF->GetBinError(PtBin_RunsBCDEF, AbsEtaBin_RunsBCDEF);
	
		float Error_RunsBCDEFGH, CMSTriggerSF_RunsBCDEFGH, Error_RunsGH, CMSTriggerSF_RunsGH;

		int PtBin_RunsGH, AbsEtaBin_RunsGH;

		if(year == "2016"){

			PtBin_RunsBCDEF = histo_RunsBCDEF->GetXaxis()->FindBin(pts.at(i));
                        AbsEtaBin_RunsBCDEF = histo_RunsBCDEF->GetYaxis()->FindBin(AbsEta.at(i));
			PtBin_RunsGH = histo_RunsGH->GetXaxis()->FindBin(pts.at(i));
                	AbsEtaBin_RunsGH = histo_RunsGH->GetYaxis()->FindBin(AbsEta.at(i));

			CMSTriggerSF_RunsBCDEF = histo_RunsBCDEF->GetBinContent(PtBin_RunsBCDEF, AbsEtaBin_RunsBCDEF);
                        Error_RunsBCDEF = histo_RunsBCDEF->GetBinError(PtBin_RunsBCDEF, AbsEtaBin_RunsBCDEF);
                	CMSTriggerSF_RunsGH = histo_RunsGH->GetBinContent(PtBin_RunsGH, AbsEtaBin_RunsGH);
                	Error_RunsGH = histo_RunsGH->GetBinError(PtBin_RunsGH, AbsEtaBin_RunsGH);

			CMSTriggerSF_RunsBCDEFGH = ( (CMSTriggerSF_RunsBCDEF * lumiRunBCDEF_triggerSF) + (CMSTriggerSF_RunsGH * lumiRunGH_triggerSF) ) / (lumiRunBCDEF_triggerSF * lumiRunGH_triggerSF + 1.0e-06);
	
			Error_RunsBCDEFGH = ( (Error_RunsBCDEF * lumiRunBCDEF_triggerSF) + (Error_RunsGH * lumiRunGH_triggerSF) ) / (lumiRunBCDEF_triggerSF * lumiRunGH_triggerSF + 1.0e-06);

               
		        if(UpOrDown == "Up"){

				CMSTriggerSF_RunsBCDEFGH += Error_RunsBCDEFGH;
				CMSTriggerSFOutput.push_back(CMSTriggerSF_RunsBCDEFGH);

			}
			else if(UpOrDown == "Down"){

                                 CMSTriggerSF_RunsBCDEFGH -= Error_RunsBCDEFGH;
                                 CMSTriggerSFOutput.push_back(CMSTriggerSF_RunsBCDEFGH);
                                
                        }
                        else{CMSTriggerSFOutput.push_back(CMSTriggerSF_RunsBCDEFGH);}



		}
		else if(year == "2017"){

			if(UpOrDown == "Up"){
			
				CMSTriggerSF_RunsBCDEF += Error_RunsBCDEF;
                                CMSTriggerSFOutput.push_back(CMSTriggerSF_RunsBCDEF);

			}
			else if(UpOrDown == "Down"){
                        
                                CMSTriggerSF_RunsBCDEF -= Error_RunsBCDEF;
                                CMSTriggerSFOutput.push_back(CMSTriggerSF_RunsBCDEF);
                        
                        }
			else{CMSTriggerSFOutput.push_back(CMSTriggerSF_RunsBCDEF);}

        	}
		else{cout << "only 2016 and 2017 have so far been included" << endl;}


  }


  return CMSTriggerSFOutput; 

}};






auto CMSTriggerSF_RDataFrames_Central{[&year, &CMSTriggerSF_RDataFrames2](const floats& pts, const floats& etas){

  return CMSTriggerSF_RDataFrames2(year, " ", pts, etas);

}};


auto CMSTriggerSF_RDataFrames_Up{[&year, &CMSTriggerSF_RDataFrames2](const floats& pts, const floats& etas){

  return CMSTriggerSF_RDataFrames2(year, "Up", pts, etas);

}};


auto CMSTriggerSF_RDataFrames_Down{[&year, &CMSTriggerSF_RDataFrames2](const floats& pts, const floats& etas){

  return CMSTriggerSF_RDataFrames2(year, "Down", pts, etas);

}};
