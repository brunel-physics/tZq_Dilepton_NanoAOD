#include <iostream>
#include <fstream>


auto MuonSF{[](const string& type, const string& year, const string& UpOrDown, const floats& pt, consts floats& eta){


  AbsEta = abs(eta);

  float lumiRunBCDEF = 19713.888;
  float lumiRunGH = 16146.178;

  TFile* inputfile_RunsBCDEF;
  TFile* inputfile_RunsGH;
  TH2* histo_RunsBCDEF;
  TH2* histo_RunsGH;

  if(year == "2016"){

  	if(type == "ID" || type == "ID sys"){

		//Muon ID file (runs BCDEF)
		inputfile_RunsBCDEF = new TFile("/home/eepgkkc/ScaleFactors/LeptonEfficiency/MuonSFs/2016/MuonID_EfficienciesAndSF_BCDEF.root", "READ");
		histo_RunsBCDEF = (TH2*)inputfile_RunsBCDEF->GetObjectChecked("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio", "TH2");


                //Muon ID file (runs GH)
                inputfile_RunsGH = new TFile("/home/eepgkkc/ScaleFactors/LeptonEfficiency/MuonSFs/2016/MuonID_EfficienciesAndSF_GH.root", "READ");
                histo_RunsGH = (TH2*)inputfile_RunsGH->GetObjectChecked("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio", "TH2");

        } 
	else if(type == "Iso" || type == "Iso sys"){
        
                //Muon ISO file (runs BCDEF)
                inputfile_RunsBCDEF = new TFile("/home/eepgkkc/ScaleFactors/LeptonEfficiency/MuonSFs/2016/MuonISO_EfficienciesAndSF_BCDEF.root", "READ");
                histo_RunsBCDEF = (TH2*)inputfile_RunsBCDEF->GetObjectChecked("TightISO_TightID_pt_eta/pt_abseta_ratio", "TH2");
                

                //Muon ISO file (runs GH)
                inputfile_RunsGH = new TFile("/home/eepgkkc/ScaleFactors/LeptonEfficiency/MuonSFs/2016/MuonISO_EfficienciesAndSF_GH.root", "READ");
                histo_RunsGH = (TH2*)inputfile_RunsGH->GetObjectChecked("TightISO_TightID_pt_eta/pt_abseta_ratio", "TH2");

        }
	else{cout << "Please choose either ID or ISO for the type" << endl;}	

  }
  else if(year == "2017"){

  	if(type == "ID"){

  		//Muon ID file
  		inputfile_RunsBCDEF = new TFile("/home/eepgkkc/ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ID.root", "READ");
  		histo_RunsBCDEF = (TH2*)inputfile_RunsBCDEF->GetObjectChecked("NUM_TightID_DEN_genTracks_pt_abseta", "TH2");  

  	}
 	 else if(type == "ID sys"){

  		//Muon ID sys file
  		inputfile_RunsBCDEF = new TFile("/home/eepgkkc/ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ID_syst.root", "READ");
  		histo_RunsBCDEF = (TH2*)inputfile_RunsBCDEF->GetObjectChecked("NUM_TightID_DEN_genTracks_pt_abseta", "TH2");

  	}
 	 else if(type == "ID sys (stat)"){

        	//Muon ID sys (stat)
        	inputfile_RunsBCDEF = new TFile("/home/eepgkkc/ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ID_syst.root", "READ");
  		histo_RunsBCDEF = (TH2*)inputfile_RunsBCDEF->GetObjectChecked("NUM_TightID_DEN_genTracks_pt_abseta_stat", "TH2");

  	}
  	else if(type == "ID sys (syst)"){
 
        	//Muon ID sys (syst)
        	inputfile_RunsBCDEF = new TFile("/home/eepgkkc/ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ID_syst.root", "READ");
  		histo_RunsBCDEF = (TH2*)inputfile_RunsBCDEF->GetObjectChecked("NUM_TightID_DEN_genTracks_pt_abseta_syst", "TH2");

  	}
  	else if(type == "Iso"){

  		//Muon Iso file
  		inputfile_RunsBCDEF = new TFile("/home/eepgkkc/ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ISO.root", "READ");
  		histo_RunsBCDEF = (TH2*)inputfile_RunsBCDEF->GetObjectChecked("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta", "TH2"); 

  	}
  	else if(type == "Iso sys"){
  
  		//Muon Iso sys file
  		inputfile_RunsBCDEF = new TFile("/home/eepgkkc/ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ISO_syst.root", "READ");
  		histo_RunsBCDEF = (TH2*)inputfile_RunsBCDEF->GetObjectChecked("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta", "TH2");

  	}
  	else if(type == "Iso sys (stat)"){

   		inputfile_RunsBCDEF = new TFile("/home/eepgkkc/ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ISO_syst.root", "READ");
  		histo_RunsBCDEF = (TH2*)inputfile_RunsBCDEF->GetObjectChecked("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta_stat", "TH2");
  
  	}
  	else if(type == "Iso sys (syst)"){

  		inputfile_RunsBCDEF = new TFile("/home/eepgkkc/ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ISO_syst.root", "READ");
		histo_RunsBCDEF = (TH2*)inputfile_RunsBCDEF->GetObjectChecked("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta_syst", "TH2"); 
  
  	}
  	else{cout << "Choose a function input out of: ID, ID sys, ID sys (stat), ID sys (syst), Iso sys, Iso sys (stat), Iso sys (syst)." << endl;}

  }
  else{cout << "need to add code for 2018" << endl;}






  floats MuonSFOutput{};
 
  for(int i = 0; i < pt.size(); i++){

	if(pt.at(i) >= 20 && pt.at(i) <= 120 && AbsEta.at(i) <= 2.4){

  		int PtBin_RunsBCDEF = histo_RunsBCDEF->GetXaxis()->FindBin(pt.at(i));
        	int AbsEtaBin_RunsBCDEF = histo_RunsBCDEF->GetYaxis()->FindBin(AbsEta.at(i));

        	float MuonSF_RunsBCDEF = histo_RunsBCDEF->GetBinContent(PtBin_RunsBCDEF, AbsEtaBin_RunsBCDEF);
		float Error_RunsBCDEF = histo_RunsBCDEF->GetBinError(PtBin_RunsBCDEF, AbsEtaBin_RunsBCDEF);
	
		float Error_RunsBCDEFGH, MuonSF_RunsBCDEFGH, Error_RunsGH, MuonSF_RunsGH;

		int PtBin_RunsGH, AbsEtaBin_RunsGH;

		if(year == "2016"){

			PtBin_RunsGH = histo_RunsGH->GetXaxis()->FindBin(pt.at(i));
                	AbsEtaBin_RunsGH = histo_RunsGH->GetYaxis()->FindBin(AbsEta.at(i));

                	MuonSF_RunsGH = histo_RunsGH->GetBinContent(PtBin_RunsGH, AbsEtaBin_RunsGH);
                	Error_RunsGH = histo_RunsGH->GetBinError(PtBin_RunsGH, AbsEtaBin_RunsGH);

			MuonSF_RunsBCDEFGH = ( (MuonSF_RunsBCDEF * lumiRunBCDEF) + (MuonSF_RunsGH * lumiRunGH) ) / (lumiRunBCDEF * lumiRunGH + 1.0e-06);
	
			Error_RunsBCDEFGH = ( (Error_RunsBCDEF * lumiRunBCDEF) + (Error_RunsGH * lumiRunGH) ) / (lumiRunBCDEF * lumiRunGH + 1.0e-06);



			if(type == "ID sys"){
               
				if(UpOrDown == "Up"){

					MuonSF_RunsBCDEFGH += Error_RunsBCDEFGH + 0.01;
					MuonSFOutput.push_back(MuonSF_RunsBCDEFGH);

				}
				else if(UpOrDown == "Down"){

                                        MuonSF_RunsBCDEFGH -= Error_RunsBCDEFGH - 0.01;
                                        MuonSFOutput.push_back(MuonSF_RunsBCDEFGH);
                                
                                }
				else{cout << "Select an up or down uncertainty" << endl;}

                        }
			else if(type == "Iso sys"){
		
				if(UpOrDown == "Up"){

                                        MuonSF_RunsBCDEFGH += Error_RunsBCDEFGH + 0.005;
                                        MuonSFOutput.push_back(MuonSF_RunsBCDEFGH);
                                
                                }
                                else if(UpOrDown == "Down"){

                                        MuonSF_RunsBCDEFGH -= Error_RunsBCDEFGH - 0.005;
                                        MuonSFOutput.push_back(MuonSF_RunsBCDEFGH);

                                }
				else{cout << "Select an up or down uncertainty" << endl;}


			}
                        else{MuonSFOutput.push_back(MuonSF_RunsBCDEFGH);}



		}
		else if(year == "2017"){

			if(type == "ID sys (stat)" ||
		   	   type == "ID sys (syst)" ||
		   	   type == "Iso sys (stat)" || 
		   	   type == "Iso sys (syst)"){
		
				MuonSFOutput.push_back(Error_RunsBCDEF);
			}
			else{MuonSFOutput.push_back(MuonSF_RunsBCDEF);}

        	}
		else{cout << "only 2016 and 2017 have so far been included" << endl;}


	}
	else{MuonSFOutput.push_back(1.0);}



  }


  return MuonSFOutput; 

}};










auto MuonSFTest_ID{[&MuonSF](const floats& pt, consts floats& eta){

  return MuonSF("ID", year, " ", pt, eta);
  
}};



auto MuonSFTest_Iso{[&MuonSF](const floats& pt, consts floats& eta){

  return MuonSF("ID", year, " ", pt, eta);
    
}};




auto MuonSFTest_IDsys_syst{[&MuonSF](const floats& pt, consts floats& eta){
  
  if(year == "2016"){
        return MuonSF("ID sys", year, "Up", pt, eta);
  }
  else if(year == "2017"){
        return MuonSF("ID sys (syst)", year, " ", pt, eta);
  }  
  else{cout << "Need to add 2018" << endl;}

}};




auto MuonSFTest_IDsys_stat{[&MuonSF](const floats& pt, consts floats& eta){

  if(year == "2016"){
	return MuonSF("ID sys", year, "Down", pt, eta);
  }
  else if(year == "2017"){
  	return MuonSF("ID sys (stat)", year, " ", pt, eta);
  }
  else{cout << "Need to add 2018" << endl;}  

}};




auto MuonSFTest_Isosys_syst{[&MuonSF](const floats& pt, consts floats& eta){

  if(year == "2016"){
        return MuonSF("Iso sys", year, "Up", pt, eta);
  }
  else if(year == "2017"){
        return MuonSF("Iso sys (syst)", year, " ", pt, eta);
  }
  else{cout << "Need to add 2018" << endl;}

}};




auto MuonSFTest_Isosys_stat{[&MuonSF](const floats& pt, consts floats& eta){

  if(year == "2016"){
        return MuonSF("Iso sys", year, "Down", pt, eta);
  }
  else if(year == "2017"){
        return MuonSF("Iso sys (stat)", year, " ", pt, eta);
  }
  else{cout << "Need to add 2018" << endl;}

}};
