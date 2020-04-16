#include <iostream>
#include <fstream>

//script won't compile but made so that code can be copied for RDataFrame column inputs
//
auto MuonSFTest2(const string& type, const string& year, const string& run, const floats& MuonPt, const floats& MuonEta){


  floats AbsEta = abs(MuonEta);

  TFile* inputfile;
  TH2* histo;

  if(year == "2016"){

  	if((type == "ID" || type == "ID sys") && run == "BCDEF"){

		//Muon ID file (runs BCDEF)
		inputfile = new TFile("/home/eepgkkc/ScaleFactors/LeptonEfficiency/MuonSFs/2016/MuonID_EfficienciesAndSF_BCDEF.root", "READ");
		histo = (TH2*)inputfile->GetObjectChecked("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio", "TH2");

  	}
	else if((type == "ID" || type == "ID sys") && run == "GH"){

                //Muon ID file (runs GH)
                inputfile = new TFile("/home/eepgkkc/ScaleFactors/LeptonEfficiency/MuonSFs/2016/MuonID_EfficienciesAndSF_GH.root", "READ");
                histo = (TH2*)inputfile->GetObjectChecked("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio", "TH2");

        } 
	else if((type == "Iso" || type == "Iso sys") && run == "BCDEF"){
        
                //Muon ISO file (runs BCDEF)
                inputfile = new TFile("/home/eepgkkc/ScaleFactors/LeptonEfficiency/MuonSFs/2016/MuonISO_EfficienciesAndSF_BCDEF.root", "READ");
                histo = (TH2*)inputfile->GetObjectChecked("TightISO_TightID_pt_eta/pt_abseta_ratio", "TH2");
                
        }
	else if((type == "Iso" || type == "Iso sys") && run == "GH"){
                
                //Muon ISO file (runs GH)
                inputfile = new TFile("/home/eepgkkc/ScaleFactors/LeptonEfficiency/MuonSFs/2016/MuonISO_EfficienciesAndSF_GH.root", "READ");
                histo = (TH2*)inputfile->GetObjectChecked("TightISO_TightID_pt_eta/pt_abseta_ratio", "TH2");

        }
	else{cout << "Please choose either ID or ISO for the type, and BCDEF or GH for the run" << endl;}	

  }
  else if(year == "2017"){

  	if(type == "ID"){

  		//Muon ID file
  		inputfile = new TFile("/home/eepgkkc/ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ID.root", "READ");
  		histo = (TH2*)inputfile->GetObjectChecked("NUM_TightID_DEN_genTracks_pt_abseta", "TH2");  

  	}
 	 else if(type == "ID sys"){

  		//Muon ID sys file
  		inputfile = new TFile("/home/eepgkkc/ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ID_syst.root", "READ");
  		histo = (TH2*)inputfile->GetObjectChecked("NUM_TightID_DEN_genTracks_pt_abseta", "TH2");

  	}
 	 else if(type == "ID sys (stat)"){

        	//Muon ID sys (stat)
        	inputfile = new TFile("/home/eepgkkc/ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ID_syst.root", "READ");
  		histo = (TH2*)inputfile->GetObjectChecked("NUM_TightID_DEN_genTracks_pt_abseta_stat", "TH2");

  	}
  	else if(type == "ID sys (syst)"){
 
        	//Muon ID sys (syst)
        	inputfile = new TFile("/home/eepgkkc/ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ID_syst.root", "READ");
  		histo = (TH2*)inputfile->GetObjectChecked("NUM_TightID_DEN_genTracks_pt_abseta_syst", "TH2");

  	}
  	else if(type == "Iso"){

  		//Muon Iso file
  		inputfile = new TFile("/home/eepgkkc/ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ISO.root", "READ");
  		histo = (TH2*)inputfile->GetObjectChecked("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta", "TH2"); 

  	}
  	else if(type == "Iso sys"){
  
  		//Muon Iso sys file
  		inputfile = new TFile("/home/eepgkkc/ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ISO_syst.root", "READ");
  		histo = (TH2*)inputfile->GetObjectChecked("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta", "TH2");

  	}
  	else if(type == "Iso sys (stat)"){

   		inputfile = new TFile("/home/eepgkkc/ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ISO_syst.root", "READ");
  		histo = (TH2*)inputfile->GetObjectChecked("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta_stat", "TH2");
  
  	}
  	else if(type == "Iso sys (syst)"){

  		inputfile = new TFile("/home/eepgkkc/ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ISO_syst.root", "READ");
		histo = (TH2*)inputfile->GetObjectChecked("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta_syst", "TH2"); 
  
  	}
  	else{cout << "Choose a function input out of: ID, ID sys, ID sys (stat), ID sys (syst), Iso sys, Iso sys (stat), Iso sys (syst)." << endl;}

  }
  else{cout << "need to add code for 2016 and 2018" << endl;}



  vector<float> MuonSFOutput{};
 
  for(int i = 0; i < MuonPt.size(); i++){

	if(MuonPt >= 20 && MuonPt.at(i) <= 120 && AbsEta.at(i) <= 2.4){

  		int PtBin = histo->GetXaxis()->FindBin(MuonPt.at(i));
        	int AbsEtaBin = histo->GetYaxis()->FindBin(AbsEta.at(i));

        	float MuonSF = histo->GetBinContent(PtBin, AbsEtaBin);
		float Error = histo->GetBinError(PtBin, AbsEtaBin);

		if(type == "ID sys (stat)" ||
		   type == "ID sys (syst)" ||
		   type == "Iso sys (stat)" || 
		   type == "Iso sys (syst)"){
		
			MuonSFOutput.push_back(Error);
		}
		else{MuonSFOutput.push_back(MuonSF);}

        }
	else{MuonSFOutput.push_back(1.0);}


  }



  return MuonSFOutput;

}






auto MuonSFTest_ID{[&year, &run](const floats& MuonPt, const floats& MuonEta){


  if(year == "2016" && run == "BCDEF"){
  	return MuonSFTest2("ID", year, "BCDEF", MuonPt, MuonEta);
  }
  else if(year == "2016" && run == "GH"){
	return MuonSFTest2("ID", year, "GH", MuonPt, MuonEta);
  }
  else if(year == "2017"){
 	return MuonSFTest2("ID", year, "", MuonPt, MuonEta);
  }
  else{cout << "Muon ID SF only coded for 2016 and 2017" << endl;}


}};


auto MuonSFTest_Iso{[&year, &run](const floats& MuonPt, const floats& MuonEta){


  if(year == "2016" && run == "BCDEF"){
        return MuonSFTest2("Iso", year, "BCDEF", MuonPt, MuonEta);
  }     
  else if(year == "2016" && run == "GH"){
        return MuonSFTest2("Iso", year, "GH", MuonPt, MuonEta);
  }
  else if(year == "2017"){
        return MuonSFTest2("Iso", year, "", MuonPt, MuonEta);
  }
  else{cout << "Muon Iso SF only coded for 2016 and 2017" << endl;}


}};


auto MuonSFTest_ID_stat{[&year, &run](const floats& MuonPt, const floats& MuonEta){

  if(year == "2016" && run == "BCDEF"){
        return MuonSFTest2("ID sys", year, "BCDEF", MuonPt, MuonEta);
  }     
  else if(year == "2016" && run == "GH"){
        return MuonSFTest2("ID sys", year, "GH", MuonPt, MuonEta);
  }
  else if(year == "2017"){
        return MuonSFTest2("ID sys (stat)", year, "", MuonPt, MuonEta);
  }
  else{cout << "Muon ID SF stat only coded for 2016 and 2017" << endl;}


}};


auto MuonSFTest_ID_syst{[&year, &run](const floats& MuonPt, const floats& MuonEta){

  if(year == "2016" && run == "BCDEF"){
        return MuonSFTest2("ID sys", year, "BCDEF", MuonPt, MuonEta);
  }
  else if(year == "2016" && run == "GH"){
        return MuonSFTest2("ID sys", year, "GH", MuonPt, MuonEta);
  }
  else if(year == "2017"){
        return MuonSFTest2("ID sys (syst)", year, "", MuonPt, MuonEta);
  }
  else{cout << "Muon ID SF syst only coded for 2016 and 2017" << endl;}



}};



auto MuonSFTest_Iso_stat{[&year, &run](const floats& MuonPt, const floats& MuonEta){

  if(year == "2016" && run == "BCDEF"){
        return MuonSFTest2("Iso sys", year, "BCDEF", MuonPt, MuonEta);
  }
  else if(year == "2016" && run == "GH"){
        return MuonSFTest2("Iso sys", year, "GH", MuonPt, MuonEta);
  }
  else if(year == "2017"){
        return MuonSFTest2("Iso sys (stat)", year, "", MuonPt, MuonEta);
  }
  else{cout << "Muon Iso SF stat only coded for 2016 and 2017" << endl;}


}};




auto MuonSFTest_Iso_syst{[&year, &run](const floats& MuonPt, const floats& MuonEta){

  if(year == "2016" && run == "BCDEF"){
        return MuonSFTest2("Iso sys", year, "BCDEF", MuonPt, MuonEta);
  }
  else if(year == "2016" && run == "GH"){
        return MuonSFTest2("Iso sys", year, "GH", MuonPt, MuonEta);
  }
  else if(year == "2017"){
        return MuonSFTest2("Iso sys (syst)", year, "", MuonPt, MuonEta);
  }
  else{cout << "Muon Iso SF syst only coded for 2016 and 2017" << endl;}


}};
