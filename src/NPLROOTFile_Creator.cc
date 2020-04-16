#include<fstream>
#include<iostream>
#include<TFile.h>

using namespace std;


auto Hadd(const string& year, const bool& blinding){

  //tZq (signal), tHq, ttW, ttZ, WZ are used for the ratio

  if(blinding == true){

 	if(year == "2016"){

                gSystem->Exec("hadd Results_MCRatio_2016_ee_NPL_Blinded.root Results_t*q_2016_ee_NPL_Blinded.root Results_ttbarV_ttbarW_*_2016_ee_NPL_Blinded.root Results_ttbarV_ttbarZ_*_2016_ee_NPL_Blinded.root Results_Diboson_WZ_*_2016_ee_NPL_Blinded.root");
	
		gSystem->Exec("hadd Results_MCRatio_2016_mumu_NPL_Blinded.root Results_t*q_2016_mumu_NPL_Blinded.root Results_ttbarV_ttbarW_*_2016_mumu_NPL_Blinded.root Results_ttbarV_ttbarZ_*_2016_mumu_NPL_Blinded.root Results_Diboson_WZ_*_2016_mumu_NPL_Blinded.root");

		gSystem->Exec("hadd Results_AllMC_2016_ee_NPL_Blinded.root Results_MC_*_2016_ee_NPL_Blinded.root");

		gSystem->Exec("hadd Results_AllMC_2016_mumu_NPL_Blinded.root Results_MC_*_2016_mumu_NPL_Blinded.root");

		gSystem->Exec("hadd Results_AllData_2016_ee_NPL_Blinded.root Results_data_*_2016_ee_NPL_Blinded.root");

                gSystem->Exec("hadd Results_AllData_2016_mumu_NPL_Blinded.root Results_data_*_2016_mumu_NPL_Blinded.root");

        }
        else if(year == "2017"){
         
	        gSystem->Exec("hadd Results_MCRatio_2017_ee_NPL_Blinded.root Results_t*q_2017_ee_NPL_Blinded.root Results_ttbarV_ttbarW_*_2017_ee_NPL_Blinded.root Results_ttbarV_ttbarZ_*_2017_ee_NPL_Blinded.root Results_Diboson_WZ_*_2017_ee_NPL_Blinded.root");

		gSystem->Exec("hadd Results_MCRatio_2017_mumu_NPL_Blinded.root Results_t*q_2017_mumu_NPL_Blinded.root Results_ttbarV_ttbarW_*_2017_mumu_NPL_Blinded.root Results_ttbarV_ttbarZ_*_2017_mumu_NPL_Blinded.root Results_Diboson_WZ_*_2017_mumu_NPL_Blinded.root");        

		gSystem->Exec("hadd Results_AllMC_2017_ee_NPL_Blinded.root Results_MC_*_2017_ee_NPL_Blinded.root");

                gSystem->Exec("hadd Results_AllMC_2017_mumu_NPL_Blinded.root Results_MC_*_2017_mumu_NPL_Blinded.root");

		gSystem->Exec("hadd Results_AllData_2017_ee_NPL_Blinded.root Results_data_*_2017_ee_NPL_Blinded.root");

                gSystem->Exec("hadd Results_AllData_2017_mumu_NPL_Blinded.root Results_data_*_2017_mumu_NPL_Blinded.root");

	}
        else if(year == "2018"){
        
	        gSystem->Exec("hadd Results_MCRatio_2018_ee_NPL_Blinded.root Results_t*q_2018_ee_NPL_Blinded.root Results_ttbarV_ttbarW_*_2018_ee_NPL_Blinded.root Results_ttbarV_ttbarZ_*_2018_ee_NPL_Blinded.root Results_Diboson_WZ_*_2018_ee_NPL_Blinded.root");

		gSystem->Exec("hadd Results_MCRatio_2018_mumu_NPL_Blinded.root Results_t*q_2018_mumu_NPL_Blinded.root Results_ttbarV_ttbarW_*_2018_mumu_NPL_Blinded.root Results_ttbarV_ttbarZ_*_2018_mumu_NPL_Blinded.root Results_Diboson_WZ_*_2018_mumu_NPL_Blinded.root");

		gSystem->Exec("hadd Results_AllMC_2018_ee_NPL_Blinded.root Results_MC_*_2018_ee_NPL_Blinded.root");

                gSystem->Exec("hadd Results_AllMC_2018_mumu_NPL_Blinded.root Results_MC_*_2018_mumu_NPL_Blinded.root");

		gSystem->Exec("hadd Results_AllData_2018_ee_NPL_Blinded.root Results_data_*_2018_ee_NPL_Blinded.root");

                gSystem->Exec("hadd Results_AllData_2018_mumu_NPL_Blinded.root Results_data_*_2018_mumu_NPL_Blinded.root");

        }
        else{cout << "Choose a year out of 2016, 2017 or 2018" << endl;} 


  }
  else{

  	if(year == "2016"){
 
                gSystem->Exec("hadd Results_MCRatio_2016_ee_NPL.root Results_t*q_2016_ee_NPL.root Results_ttbarV_ttbarW_*_2016_ee_NPL.root Results_ttbarV_ttbarZ_*_2016_ee_NPL.root Results_Diboson_WZ_*_2016_ee_NPL.root");
 
		gSystem->Exec("hadd Results_MCRatio_2016_mumu_NPL.root Results_t*q_2016_mumu_NPL.root Results_ttbarV_ttbarW_*_2016_mumu_NPL.root Results_ttbarV_ttbarZ_*_2016_mumu_NPL.root Results_Diboson_WZ_*_2016_mumu_NPL.root"); 

		gSystem->Exec("hadd Results_AllMC_2016_ee_NPL.root Results_MC_*_2016_ee_NPL.root");

                gSystem->Exec("hadd Results_AllMC_2016_mumu_NPL.root Results_MC_*_2016_mumu_NPL.root");

		gSystem->Exec("hadd Results_AllData_2016_ee_NPL.root Results_data_*_2016_ee_NPL.root");

                gSystem->Exec("hadd Results_AllData_2016_mumu_NPL.root Results_data_*_2016_mumu_NPL.root");

         }
        else if(year == "2017"){
        
	        gSystem->Exec("hadd Results_MCRatio_2017_ee_NPL.root Results_t*q_2017_ee_NPL.root Results_ttbarV_ttbarW_*_2017_ee_NPL.root Results_ttbarV_ttbarZ_*_2017_ee_NPL.root Results_Diboson_WZ_*_2017_ee_NPL.root");
     
		gSystem->Exec("hadd Results_MCRatio_2017_mumu_NPL.root Results_t*q_2017_mumu_NPL.root Results_ttbarV_ttbarW_*_2017_mumu_NPL.root Results_ttbarV_ttbarZ_*_2017_mumu_NPL.root Results_Diboson_WZ_*_2017_mumu_NPL.root");

		gSystem->Exec("hadd Results_AllMC_2017_ee_NPL.root Results_MC_*_2017_ee_NPL.root");

                gSystem->Exec("hadd Results_AllMC_2017_mumu_NPL.root Results_MC_*_2017_mumu_NPL.root");

		gSystem->Exec("hadd Results_AllData_2017_ee_NPL.root Results_data_*_2017_ee_NPL.root");

                gSystem->Exec("hadd Results_AllData_2017_mumu_NPL.root Results_data_*_2017_mumu_NPL.root")

        }
        else if(year == "2018"){

                gSystem->Exec("hadd Results_MCRatio_2018_ee_NPL.root Results_t*q_2018_ee_NPL.root Results_ttbarV_ttbarW_*_2018_ee_NPL.root Results_ttbarV_ttbarZ_*_2018_ee_NPL.root Results_Diboson_WZ_*_2018_ee_NPL.root");

		gSystem->Exec("hadd Results_MCRatio_2018_mumu_NPL.root Results_t*q_2018_mumu_NPL.root Results_ttbarV_ttbarW_*_2018_mumu_NPL.root Results_ttbarV_ttbarZ_*_2018_mumu_NPL.root Results_Diboson_WZ_*_2018_mumu_NPL.root");

		gSystem->Exec("hadd Results_AllMC_2018_ee_NPL.root Results_MC_*_2018_ee_NPL.root");

                gSystem->Exec("hadd Results_AllMC_2018_mumu_NPL.root Results_MC_*_2018_mumu_NPL.root");

		gSystem->Exec("hadd Results_AllData_2018_ee_NPL.root Results_data_*_2018_ee_NPL.root");

                gSystem->Exec("hadd Results_AllData_2018_mumu_NPL.root Results_data_*_2018_mumu_NPL.root")


        }
        else{cout << "Choose a year out of 2016, 2017 or 2018" << endl;}


  }


}





auto NPLROOTFile_Creator2(const string& year, const bool& blinding){

 Hadd(year, blinding);

 TFile* AllMC_ee, AllMC_mumu, AllData_ee, AllData_mumu, MCRatio_ee, MCRatio_mumu;

 if(blinding == true){

        AllMC_ee = TFile::Open("Results_AllMC_" + year + "_ee_NPL_Blinded.root", "READ");
        AllMC_mumu = TFile::Open("Results_AllMC_" + year + "_mumu_NPL_Blinded.root", "READ");
	AllData_ee = TFile::Open("Results_AllData_" + year + "_ee_NPL_Blinded.root", "READ");
        AllData_mumu = TFile::Open("Results_AllData_" + year + "_mumu_NPL_Blinded.root", "READ");
	MCRatio_ee = TFile::Open("Results_MCRatio_" + year + "_ee_NPL_Blinded.root", "READ");
        MCRatio_mumu = TFile::Open("Results_MCRatio_" + year + "_mumu_NPL_Blinded.root", "READ");

 }
 else{
	AllMC_ee = TFile::Open("Results_AllMC_" + year + "_ee_NPL.root", "READ");
	AllMC_mumu = TFile::Open("Results_AllMC_" + year + "_mumu_NPL.root", "READ");
	AllData_ee = TFile::Open("Results_AllData_" + year + "_ee_NPL.root", "READ");
        AllData_mumu = TFile::Open("Results_AllData_" + year + "_mumu_NPL.root", "READ");
	MCRatio_ee = TFile::Open("Results_MCRatio_" + year + "_ee_NPL.root", "READ");
        MCRatio_mumu = TFile::Open("Results_MCRatio_" + year + "_mumu_NPL.root", "READ");

 }
 

 //For the subtraction between data and MC
 TH1* h_NData_SS_ee = (TH1*)AllData_ee->GetObjectChecked("SameSign", "TH1"); 
 TH1* h_NData_SS_mumu = (TH1*)AllData_mumu->GetObjectChecked("SameSign", "TH1");
 TH1* h_NMC_SS_ee = (TH1*)AllMC_ee->GetObjectChecked("SameSign", "TH1");
 TH1* h_NMC_SS_mumu = (TH1*)AllMC_mumu->GetObjectChecked("SameSign", "TH1");

 //For the ratio
 TH1* h_NMC_OS_NonPrompt_ee = (TH1*)MCRatio_ee->GetObjectChecked("OppositeSign", "TH1");
 TH1* h_NMC_OS_NonPrompt_mumu = (TH1*)MCRatio_mumu->GetObjectChecked("OppositeSign", "TH1");
 TH1* h_NMC_SS_NonPrompt_ee = (TH1*)MCRatio_ee->GetObjectChecked("SameSign", "TH1");
 TH1* h_NMC_SS_NonPrompt_mumu = (TH1*)MCRatio_mumu->GetObjectChecked("SameSign", "TH1");

 int NMC_OS_NonPrompt_ee = h_NMC_OS_NonPrompt_ee->GetNEntries();
 int NMC_OS_NonPrompt_mumu = h_NMC_OS_NonPrompt_mumu->GetNEntries();
 int h_NMC_SS_NonPrompt_ee = h_h_NMC_SS_NonPrompt_ee->GetNEntries();
 int h_NMC_SS_NonPrompt_mumu = h_h_NMC_SS_NonPrompt_mumu->GetNEntries();

 int ratio_ee = NMC_OS_NonPrompt_ee / NMC_SS_NonPrompt_ee;
 int ratio_mumu = NMC_OS_NonPrompt_mumu / NMC_SS_NonPrompt_mumu;

 //Calculating the number of opposite-sign non prompt data entries
 TH1* h_NData_OS_NonPrompt_ee; 
 TH1* h_NData_OS_NonPrompt_mumu;

 int nbins = h_NData_SS_ee->GetNbinsX();

 for(int i = 0; i < nbins; i++){

	h_NData_OS_NonPrompt_ee->SetBinContent.at(i) = ( ( h_NData_SS_ee->GetBinContent.at(i) - h_NMC_SS_ee->GetBinContent.at(i) ) * ratio_ee);
   	h_NData_OS_NonPrompt_mumu->SetBinContent.at(i) = ( ( h_NData_SS_mumu->GetBinContent.at(i) - h_NMC_SS_mumu->GetBinContent.at(i) ) * ratio_mumu);  

 }


 //Saving the histograms to an output file
 string NPL_output_file;

 if(blinding == true){NPL_output_file = "NPL_ee_output_" + year + "_" + "Blinded.root";}
 else{NPL_output_file = "NPL_ee_output_" + year + ".root";}

 TFile * NPL_output = new TFile(NPL_output_file.c_str());

 h_NData_OS_NonPrompt_ee->Write();
 h_NData_OS_NonPrompt_mumu->Write();
 h_NData_SS_ee->Write();
 h_NData_SS_mumu->Write();
 h_NMC_SS_ee->Write();
 h_NMC_SS_mumu->Write();
 h_NMC_OS_NonPrompt_ee->Write();
 h_NMC_OS_NonPrompt_mumu->Write();
 h_NMC_SS_NonPrompt_ee->Write();
 h_NMC_SS_NonPrompt_mumu->Write();


 NPL_output->Close();

}



auto NPLROOTFile_Creator(){

 bool blinding = true;
 vector<string> year = {"2016", "2017", "2018"};

 for(int i = 0; i < year.size(); i++){NPLROOTFile_Creator2(year.at(i), blinding);}

}
