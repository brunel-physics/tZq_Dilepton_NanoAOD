#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TH1D.h"
#include <boost/numeric/conversion/cast.hpp>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <iterator>
#include <string>
#include <algorithm>
#include <TFile.h>
#include "TEfficiency.h"
#include "TLorentzVector.h"
#include <iostream>
#include <fstream>
#include <random>
#include <TRandom3.h>
#include "../ScaleFactors/LeptonEnergyCorrections/RochesterCorrections/roccor.Run2.v3/RoccoR.cc"
#include "../ScaleFactors/LeptonEnergyCorrections/RochesterCorrections/roccor.Run2.v3/RoccoR.h"
#include "fulleventselectionAlgo.hpp"
#include "TCanvas.h"

fulleventselectionAlgo::fulleventselectionAlgo(){


}


fulleventselectionAlgo::~fulleventselectionAlgo()
{
}


using namespace ROOT; // RDataFrame's namespace
using ROOT::RDF::RNode;

using floats = ROOT::VecOps::RVec<float>;
using ints = ROOT::VecOps::RVec<int>;
using bools = ROOT::VecOps::RVec<bool>;
using chars = ROOT::VecOps::RVec<UChar_t>;
using doubles = ROOT::VecOps::RVec<double>;


std::fstream& GotoLine(std::fstream& file, unsigned int num){
		
    file.seekg(std::ios::beg);

    for(int i=0; i < num - 1; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }

    return file;

}




//For b tagging reweighting
int line_number = 0;
int Index = 0;
int output_index;
std::vector<int> line_number_vec{};

class CSVReader
{
	std::string fileName;
	std::string delimeter;
 
public:
	CSVReader(std::string filename, std::string delm = ",") :
			fileName(filename), delimeter(delm)
	{ }
 
	std::vector<std::vector<std::string> > getData();

};
 


std::vector<std::vector<std::string> > CSVReader::getData()
{
	std::ifstream file(fileName);
	std::vector<std::vector<std::string> > dataList;
	std::string line = "";
	int line_number = 0;
	
	while (getline(file, line))
	{
		line_number++;
		line_number_vec.push_back(line_number);
		std::vector<std::string> vec;
		boost::algorithm::split(vec, line, boost::is_any_of(delimeter));
		dataList.push_back(vec);
	}
	
	file.close();
	
	return dataList;

}
//b tagging reweighting ends here

//For the Rochester corrections
auto RochesterCorrections_testscript2(

const std::string& year, 
const std::string& process, 
const ints& MuonCharge, 
const floats& MuonPt,
const floats& MuonEta,
const floats& MuonPhi,
const ints& Muon_genPartIdx,
const ints& Muon_nTrackerLayers){


	std::string RoccoTextFile;

	if(year == "2016"){RoccoTextFile = "./ScaleFactors/LeptonEnergyCorrections/RochesterCorrections/roccor.Run2.v3/RoccoR2016.txt";}
	else if(year == "2017"){RoccoTextFile = "./ScaleFactors/LeptonEnergyCorrections/RochesterCorrections/roccor.Run2.v3/RoccoR2017.txt";}
	else if(year == "2018"){RoccoTextFile = "./ScaleFactors/LeptonEnergyCorrections/RochesterCorrections/roccor.Run2.v3/RoccoR2018.txt";}
	else{std::cout << "Error for rochester corrections: choose a year out of 2016, 2017 or 2018." << std::endl;}


	RoccoR rc{RoccoTextFile};
	doubles RochCorrVec{};

	ints s(MuonPt.size(), 0); //s is error set (default is 0)
	ints m(MuonPt.size(), 0); //m is error member (default is 0, ranges from 0 to nmembers-1)
	doubles u{gRandom->Rndm(), gRandom->Rndm()}; //u is a random number distributed uniformly between 0 and 1 (gRandom->Rndm());


	for(int i = 0; i < MuonPt.size(); i++){

		//scale factors for momentum of each muon:
		double RochCorrSF;
		double mcSF;


		if(process != "data_DoubleEGRunB" &&
		   process != "data_DoubleEGRunC" &&
		   process != "data_DoubleEGRunD" &&
		   process != "data_DoubleEGRunE" &&
		   process != "data_DoubleEGRunF" &&
		   process != "data_DoubleEGRunG" &&
		   process != "data_DoubleEGRunH" &&
		   process != "data_EGRunB" &&
                   process != "data_EGRunC" &&
                   process != "data_EGRunD" &&
		   process != "data_SingleElectronRunB" &&
		   process != "data_SingleElectronRunC" &&
		   process != "data_SingleElectronRunD" &&
		   process != "data_SingleElectronRunE" &&
		   process != "data_SingleElectronRunF" &&
		   process != "data_SingleElectronRunG" &&
		   process != "data_SingleElectronRunH" &&
		   process != "data_DoubleMuonRunB" &&
		   process != "data_DoubleMuonRunC" &&
		   process != "data_DoubleMuonRunD" &&
		   process != "data_DoubleMuonRunE" &&
		   process != "data_DoubleMuonRunF" &&
		   process != "data_DoubleMuonRunG" &&
		   process != "data_DoubleMuonRunH" &&
		   process != "data_SingleMuonRunB" &&
		   process != "data_SingleMuonRunC" &&
		   process != "data_SingleMuonRunD" &&
		   process != "data_SingleMuonRunE" &&
		   process != "data_SingleMuonRunF" &&
		   process != "data_SingleMuonRunG" &&
		   process != "data_SingleMuonRunH"){



			if(mcSF > 0){

				RochCorrSF = rc.kSpreadMC(MuonCharge.at(i), MuonPt.at(i), MuonEta.at(i), MuonPhi.at(i), Muon_genPartIdx.at(i), s.at(i), m.at(i)); //(recommended), MC scale and resolution correction when matched gen muon is available
			}
			else{
				RochCorrSF = rc.kSmearMC(MuonCharge.at(i), MuonPt.at(i), MuonEta.at(i), MuonPhi.at(i), Muon_nTrackerLayers.at(i), u.at(i), s.at(i), m.at(i)); //MC scale and extra smearing when matched gen muon is not available

			}

		}
		else{RochCorrSF = rc.kScaleDT(MuonCharge.at(i), MuonPt.at(i), MuonEta.at(i), MuonPhi.at(i), s.at(i), m.at(i));} //data

	
		RochCorrVec.push_back(RochCorrSF);


	}

	
	return RochCorrVec;


}


//Functions for reading the trigger efficiency and SF text files
double linereader_TriggerSF(const int& LineNumber, const std::string& InputTriggerSF_File, const std::string& year, const bool& blinding){

   std::string TriggerSF_TextFiles;

   if(InputTriggerSF_File == "Data_Central"){TriggerSF_TextFiles = "TriggerSF_Efficiency_Data_MET_" + year + ".txt";}
   else if(InputTriggerSF_File == "MC_Central"){TriggerSF_TextFiles = "TriggerSF_Efficiency_MC_ttbar_" + year + ".txt";}
   else if(InputTriggerSF_File == "Data_Uncert"){TriggerSF_TextFiles = "TriggerSF_EfficiencyUncerts_Data_MET_" + year + ".txt";}
   else if(InputTriggerSF_File == "MC_Uncert"){TriggerSF_TextFiles = "TriggerSF_EfficiencyUncerts_MC_ttbar_" + year + ".txt";}
   else if(InputTriggerSF_File == "SF_Central" || InputTriggerSF_File == "SF_Uncert"){TriggerSF_TextFiles = "TriggerSF_ScaleFactors_" + year + ".txt";}
   else{std::cout << "please choose an appropriate input text file for trigger SFs" << std::endl;}


   using namespace std;

   fstream file(TriggerSF_TextFiles.c_str());
   GotoLine(file, LineNumber);

   std::string line;
   file >> line;

   double Value = atof(line.c_str());
   return Value;

}




int linecounter_TriggerSF(const std::string& InputTriggerSF_File, const std::string& year, const bool& blinding){

   std::string TriggerSF_TextFiles;

 
  if(InputTriggerSF_File == "Data_Central"){TriggerSF_TextFiles = "TriggerSF_Efficiency_Data_MET_" + year + ".txt";}
   else if(InputTriggerSF_File == "MC_Central"){TriggerSF_TextFiles = "TriggerSF_Efficiency_MC_ttbar_" + year + ".txt";}
   else if(InputTriggerSF_File == "Data_Uncert"){TriggerSF_TextFiles = "TriggerSF_EfficiencyUncerts_Data_MET_" + year + ".txt";}
   else if(InputTriggerSF_File == "MC_Uncert"){TriggerSF_TextFiles = "TriggerSF_EfficiencyUncerts_MC_ttbar_" + year + ".txt";}
   else if(InputTriggerSF_File == "SF_Central" || InputTriggerSF_File == "SF_Uncert"){TriggerSF_TextFiles = "TriggerSF_ScaleFactors_" + year + ".txt";}
   else{std::cout << "please choose an appropriate input text file for trigger SFs" << std::endl;}



   int number_of_lines = 0;
   std::string line;
   std::ifstream myfile(TriggerSF_TextFiles.c_str());

   while (getline(myfile, line))
        ++number_of_lines;
        return number_of_lines;

}



auto textfilereader2_TriggerSF(const std::string& InputTriggerSF_File, const std::string& year, const bool& blinding){


   int NumberOfLines = linecounter_TriggerSF(InputTriggerSF_File, year, blinding);
   std::vector<double> Value;


   for(int i = 1; i < NumberOfLines+1; i++){
        Value.push_back(linereader_TriggerSF(i, InputTriggerSF_File, year, blinding));
   }

   return Value;

}




//For creating the output directories
//Creating a directory for all results
auto DirectoryCreator(const std::string& year, const bool& blinding, const bool& NPL){

	if(year == "2016"){

		if(blinding == true){

			if(NPL == false){

                        	system("mkdir BlindedResults_2016");
                        	system("mv *2016*.root BlindedResults_2016");
				system("mv *2016*.txt BlindedResults_2016");
				system("mkdir BTagEffPlots");
				system("mkdir Chi2Range");
				system("mkdir CutFlowReport");
				system("mkdir Resolution");
				system("mkdir Results");
				system("mkdir TriggerSF");
				system("mkdir GaussianFit");
				system("mkdir mW_versus_mTop");
				system("mv BTagEffPlots BlindedResults_2016/");
                                system("mv Chi2Range BlindedResults_2016/");
                                system("mv CutFlowReport BlindedResults_2016/");
                                system("mv Resolution BlindedResults_2016/");
                                system("mv Results BlindedResults_2016/");
                                system("mv TriggerSF BlindedResults_2016/");
                                system("mv GaussianFit BlindedResults_2016/");
                                system("mv mW_versus_mTop BlindedResults_2016/");
				system("mv BlindedResults_2016/BTagEffPlots_*.root BlindedResults_2016/BTagEffPlots");
				system("mv BlindedResults_2016/Chi2Range_*.txt BlindedResults_2016/Chi2Range");
				system("mv BlindedResults_2016/CutFlowReport_*.txt BlindedResults_2016/CutFlowReport");
				system("mv BlindedResults_2016/Resolution_*.txt BlindedResults_2016/Resolution");
				system("mv BlindedResults_2016/Results_*.root BlindedResults_2016/Results");
				system("mv BlindedResults_2016/TurnOnCurves_*.root BlindedResults_2016/TriggerSF");
				system("mv BlindedResults_2016/TriggerSF_*.txt BlindedResults_2016/TriggerSF");
				system("mv BlindedResults_2016/*_GaussianFit_*.root BlindedResults_2016/GaussianFit");
				system("mv BlindedResults_2016/*_mW_mTop_*.root BlindedResults_2016/mW_versus_mTop");

			}
			else{

				system("mkdir BlindedResults_NPL_2016");
                        	system("mv *NPL*.root BlindedResults_NPL_2016");
                        	system("mv *NPL*.root BlindedResults_NPL_2016");
				system("mkdir BTagEffPlots");
                                system("mkdir Chi2Range");
                                system("mkdir CutFlowReport");
                                system("mkdir Resolution");
                                system("mkdir Results");
                                system("mkdir TriggerSF");
                                system("mkdir GaussianFit");
                                system("mkdir mW_versus_mTop");
                                system("mv BTagEffPlots BlindedResults_NPL_2016/");
                                system("mv Chi2Range BlindedResults_NPL_2016/");
                                system("mv CutFlowReport BlindedResults_NPL_2016/");
                                system("mv Resolution BlindedResults_NPL_2016/");
                                system("mv Results BlindedResults_NPL_2016/");
                                system("mv TriggerSF BlindedResults_NPL_2016/");
                                system("mv GaussianFit BlindedResults_NPL_2016/");
                                system("mv mW_versus_mTop BlindedResults_NPL_2016/");
				system("mv UnblindedResults_NPL_2016/BTagEffPlots_*.root BlindedResults_NPL_2016/BTagEffPlots");
                                system("mv BlindedResults_NPL_2016/Chi2Range_*.txt BlindedResults_NPL_2016/Chi2Range");
                                system("mv BlindedResults_NPL_2016/CutFlowReport_*.txt BlindedResults_NPL_2016/CutFlowReport");
                                system("mv BlindedResults_NPL_2016/Resolution_*.txt BlindedResults_NPL_2016/Resolution");
                                system("mv BlindedResults_NPL_2016/Results_*.root BlindedResults_NPL_2016/Results");
                                system("mv BlindedResults_NPL_2016/TurnOnCurves_*.root BlindedResults_NPL_2016/TriggerSF");
                                system("mv BlindedResults_NPL_2016/TriggerSF_*.txt BlindedResults_NPL_2016/TriggerSF");
				system("mv BlindedResults_NPL_2016/*_GaussianFit_*.root BlindedResults_NPL_2016/GaussianFit");
                                system("mv BlindedResults_NPL_2016/*_mW_mTop_*.root BlindedResults_NPL_2016/mW_versus_mTop");

			}
	
                }
                else{

			if(NPL == false){

                        	system("mkdir UnblindedResults_2016");
                        	system("mv *2016*.root UnblindedResults_2016");
				system("mv *2016*.txt UnblindedResults_2016");
				system("mkdir BTagEffPlots");
                                system("mkdir Chi2Range");
                                system("mkdir CutFlowReport");
                                system("mkdir Resolution");
                                system("mkdir Results");
                                system("mkdir TriggerSF");
                                system("mkdir GaussianFit");
                                system("mkdir mW_versus_mTop");
                                system("mv BTagEffPlots UnblindedResults_2016/");
                                system("mv Chi2Range UnblindedResults_2016/");
                                system("mv CutFlowReport UnblindedResults_2016/");
                                system("mv Resolution UnblindedResults_2016/");
                                system("mv Results UnblindedResults_2016/");
                                system("mv TriggerSF UnblindedResults_2016/");
                                system("mv GaussianFit UnblindedResults_2016/");
                                system("mv mW_versus_mTop UnblindedResults_2016/");
                                system("mv UnblindedResults_2016/BTagEffPlots_*.root UnblindedResults_2016/BTagEffPlots");
                                system("mv UnblindedResults_2016/Chi2Range_*.txt UnblindedResults_2016/Chi2Range");
                                system("mv UnblindedResults_2016/CutFlowReport_*.txt UnblindedResults_2016/CutFlowReport");
                                system("mv UnblindedResults_2016/Resolution_*.txt UnblindedResults_2016/Resolution");
                                system("mv UnblindedResults_2016/Results_*.root UnblindedResults_2016/Results");
                                system("mv UnblindedResults_2016/TurnOnCurves_*.root UnblindedResults_2016/TriggerSF");
                                system("mv UnblindedResults_2016/TriggerSF_*.txt UnblindedResults_2016/TriggerSF");
				system("mv UnblindedResults_2016/*_GaussianFit_*.root UnblindedResults_2016/GaussianFit");
                                system("mv UnblindedResults_2016/*_mW_mTop_*.root UnblindedResults_2016/mW_versus_mTop");


			}
			else{
			
				system("mkdir UnblindedResults_NPL_2016");
                        	system("mv *NPL*.root UnblindedResults_NPL_2016");
                        	system("mv *NPL*.txt UnblindedResults_NPL_2016");
				system("mkdir BTagEffPlots");
                                system("mkdir Chi2Range");
                                system("mkdir CutFlowReport");
                                system("mkdir Resolution");
                                system("mkdir Results");
                                system("mkdir TriggerSF");
                                system("mkdir GaussianFit");
                                system("mkdir mW_versus_mTop");
                                system("mv BTagEffPlots UnblindedResults_NPL_2016/");
                                system("mv Chi2Range UnblindedResults_NPL_2016/");
                                system("mv CutFlowReport UnblindedResults_NPL_2016/");
                                system("mv Resolution UnblindedResults_NPL_2016/");
                                system("mv Results UnblindedResults_NPL_2016/");
                                system("mv TriggerSF UnblindedResults_NPL_2016/");
                                system("mv GaussianFit UnblindedResults_NPL_2016/");
                                system("mv mW_versus_mTop UnblindedResults_NPL_2016/");
                                system("mv UnblindedResults_NPL_2016/BTagEffPlots_*.root UnblindedResults_NPL_2016/BTagEffPlots");
                                system("mv UnblindedResults_NPL_2016/Chi2Range_*.txt UnblindedResults_NPL_2016/Chi2Range");
                                system("mv UnblindedResults_NPL_2016/CutFlowReport_*.txt UnblindedResults_NPL_2016/CutFlowReport");
                                system("mv UnblindedResults_NPL_2016/Resolution_*.txt UnblindedResults_NPL_2016/Resolution");
                                system("mv UnblindedResults_NPL_2016/Results_*.root UnblindedResults_NPL_2016/Results");
                                system("mv UnblindedResults_NPL_2016/TurnOnCurves_*.root UnblindedResults_NPL_2016/TriggerSF");
                                system("mv UnblindedResults_NPL_2016/TriggerSF_*.txt UnblindedResults_NPL_2016/TriggerSF");
				system("mv UnblindedResults_NPL_2016/*_GaussianFit_*.root UnblindedResults_NPL_2016/GaussianFit");
                                system("mv UnblindedResults_NPL_2016/*_mW_mTop_*.root UnblindedResults_NPL_2016/mW_versus_mTop");
				

			}

                }


	
	}
	else if(year == "2017"){

		if(blinding == true){

			if(NPL == false){

				system("mkdir BlindedResults_2017");
                                system("mv *2017*.root BlindedResults_2017");
                                system("mv *2017*.txt BlindedResults_2017");
                                system("mkdir BTagEffPlots");
                                system("mkdir Chi2Range");
                                system("mkdir CutFlowReport");
                                system("mkdir Resolution");
                                system("mkdir Results");
                                system("mkdir TriggerSF");
                                system("mkdir GaussianFit");
                                system("mkdir mW_versus_mTop");
                                system("mv BTagEffPlots BlindedResults_2017/");
                                system("mv Chi2Range BlindedResults_2017/");
                                system("mv CutFlowReport BlindedResults_2017/");
                                system("mv Resolution BlindedResults_2017/");
                                system("mv Results BlindedResults_2017/");
                                system("mv TriggerSF BlindedResults_2017/");
                                system("mv GaussianFit BlindedResults_2017/");
                                system("mv mW_versus_mTop BlindedResults_2017/");
                                system("mv BlindedResults_2017/BTagEffPlots_*.root BlindedResults_2017/BTagEffPlots");
                                system("mv BlindedResults_2017/Chi2Range_*.txt BlindedResults_2017/Chi2Range");
                                system("mv BlindedResults_2017/CutFlowReport_*.txt BlindedResults_2017/CutFlowReport");
                                system("mv BlindedResults_2017/Resolution_*.txt BlindedResults_2017/Resolution");
                                system("mv BlindedResults_2017/Results_*.root BlindedResults_2017/Results");
                                system("mv BlindedResults_2017/TurnOnCurves_*.root BlindedResults_2017/TriggerSF");
                                system("mv BlindedResults_2017/TriggerSF_*.txt BlindedResults_2017/TriggerSF");
				system("mv BlindedResults_2017/*_GaussianFit_*.root BlindedResults_2017/GaussianFit");
                                system("mv BlindedResults_2017/*_mW_mTop_*.root BlindedResults_2017/mW_versus_mTop");
			

			}
			else{
			
				system("mkdir BlindedResults_NPL_2017");
                                system("mv *NPL*.root BlindedResults_NPL_2017");
                                system("mv *NPL*.root BlindedResults_NPL_2017");
                                system("mkdir BTagEffPlots");
                                system("mkdir Chi2Range");
                                system("mkdir CutFlowReport");
                                system("mkdir Resolution");
                                system("mkdir Results");
                                system("mkdir TriggerSF");
                                system("mkdir GaussianFit");
                                system("mkdir mW_versus_mTop");
                                system("mv BTagEffPlots BlindedResults_NPL_2017/");
                                system("mv Chi2Range BlindedResults_NPL_2017/");
                                system("mv CutFlowReport BlindedResults_NPL_2017/");
                                system("mv Resolution BlindedResults_NPL_2017/");
                                system("mv Results BlindedResults_NPL_2017/");
                                system("mv TriggerSF BlindedResults_NPL_2017/");
                                system("mv GaussianFit BlindedResults_NPL_2017/");
                                system("mv mW_versus_mTop BlindedResults_NPL_2017/");
                                system("mv UnblindedResults_NPL_2017/BTagEffPlots_*.root BlindedResults_NPL_2017/BTagEffPlots");
                                system("mv BlindedResults_NPL_2017/Chi2Range_*.txt BlindedResults_NPL_2017/Chi2Range");
                                system("mv BlindedResults_NPL_2017/CutFlowReport_*.txt BlindedResults_NPL_2017/CutFlowReport");
                                system("mv BlindedResults_NPL_2017/Resolution_*.txt BlindedResults_NPL_2017/Resolution");
                                system("mv BlindedResults_NPL_2017/Results_*.root BlindedResults_NPL_2017/Results");
                                system("mv BlindedResults_NPL_2017/TurnOnCurves_*.root BlindedResults_NPL_2017/TriggerSF");
                                system("mv BlindedResults_NPL_2017/TriggerSF_*.txt BlindedResults_NPL_2017/TriggerSF");
				system("mv BlindedResults_NPL_2017/*_GaussianFit_*.root BlindedResults_NPL_2017/GaussianFit");
                                system("mv BlindedResults_NPL_2017/*_mW_mTop_*.root BlindedResults_NPL_2017/mW_versus_mTop");
			

			}

                }
                else{

			if(NPL == false){

				system("mkdir UnblindedResults_2017");
                                system("mv *2017*.root UnblindedResults_2017");
                                system("mv *2017*.txt UnblindedResults_2017");
                                system("mkdir BTagEffPlots");
                                system("mkdir Chi2Range");
                                system("mkdir CutFlowReport");
                                system("mkdir Resolution");
                                system("mkdir Results");
                                system("mkdir TriggerSF");
                                system("mkdir GaussianFit");
                                system("mkdir mW_versus_mTop");
                                system("mv BTagEffPlots UnblindedResults_2017/");
                                system("mv Chi2Range UnblindedResults_2017/");
                                system("mv CutFlowReport UnblindedResults_2017/");
                                system("mv Resolution UnblindedResults_2017/");
                                system("mv Results UnblindedResults_2017/");
                                system("mv TriggerSF UnblindedResults_2017/");
                                system("mv GaussianFit UnblindedResults_2017/");
                                system("mv mW_versus_mTop UnblindedResults_2017/");
                                system("mv UnblindedResults_2017/BTagEffPlots_*.root UnblindedResults_2017/BTagEffPlots");
                                system("mv UnblindedResults_2017/Chi2Range_*.txt UnblindedResults_2017/Chi2Range");
                                system("mv UnblindedResults_2017/CutFlowReport_*.txt UnblindedResults_2017/CutFlowReport");
                                system("mv UnblindedResults_2017/Resolution_*.txt UnblindedResults_2017/Resolution");
                                system("mv UnblindedResults_2017/Results_*.root UnblindedResults_2017/Results");
                                system("mv UnblindedResults_2017/TurnOnCurves_*.root UnblindedResults_2017/TriggerSF");
                                system("mv UnblindedResults_2017/TriggerSF_*.txt UnblindedResults_2017/TriggerSF");
				system("mv UnblindedResults_2017/*_GaussianFit_*.root UnblindedResults_2017/GaussianFit");
                                system("mv UnblindedResults_2017/*_mW_mTop_*.root UnblindedResults_2017/mW_versus_mTop");


			}
			else{

				system("mkdir UnblindedResults_NPL_2017");
                                system("mv *NPL*.root UnblindedResults_NPL_2017");
                                system("mv *NPL*.txt UnblindedResults_NPL_2017");
                                system("mkdir BTagEffPlots");
                                system("mkdir Chi2Range");
                                system("mkdir CutFlowReport");
                                system("mkdir Resolution");
                                system("mkdir Results");
                                system("mkdir TriggerSF");
                                system("mkdir GaussianFit");
                                system("mkdir mW_versus_mTop");
                                system("mv BTagEffPlots UnblindedResults_NPL_2017/");
                                system("mv Chi2Range UnblindedResults_NPL_2017/");
                                system("mv CutFlowReport UnblindedResults_NPL_2017/");
                                system("mv Resolution UnblindedResults_NPL_2017/");
                                system("mv Results UnblindedResults_NPL_2017/");
                                system("mv TriggerSF UnblindedResults_NPL_2017/");
                                system("mv GaussianFit UnblindedResults_NPL_2017/");
                                system("mv mW_versus_mTop UnblindedResults_NPL_2017/");
                                system("mv UnblindedResults_NPL_2017/BTagEffPlots_*.root UnblindedResults_NPL_2017/BTagEffPlots");
                                system("mv UnblindedResults_NPL_2017/Chi2Range_*.txt UnblindedResults_NPL_2017/Chi2Range");
                                system("mv UnblindedResults_NPL_2017/CutFlowReport_*.txt UnblindedResults_NPL_2017/CutFlowReport");
                                system("mv UnblindedResults_NPL_2017/Resolution_*.txt UnblindedResults_NPL_2017/Resolution");
                                system("mv UnblindedResults_NPL_2017/Results_*.root UnblindedResults_NPL_2017/Results");
			        system("mv UnblindedResults_NPL_2017/TurnOnCurves_*.root UnblindedResults_NPL_2017/TriggerSF");
                                system("mv UnblindedResults_NPL_2017/TriggerSF_*.txt UnblindedResults_NPL_2017/TriggerSF");                        
				system("mv UnblindedResults_NPL_2017/*_GaussianFit_*.root UnblindedResults_NPL_2017/GaussianFit");
                                system("mv UnblindedResults_NPL_2017/*_mW_mTop_*.root UnblindedResults_NPL_2017/mW_versus_mTop");
	
	
			}


                }

	
	}
	else if(year == "2018"){

		 if(blinding == true){
			
			if(NPL == false){

				system("mkdir BlindedResults_2018");
                                system("mv *2018*.root BlindedResults_2018");
                                system("mv *2018*.txt BlindedResults_2018");
                                system("mkdir BTagEffPlots");
                                system("mkdir Chi2Range");
                                system("mkdir CutFlowReport");
                                system("mkdir Resolution");
                                system("mkdir Results");
                                system("mkdir TriggerSF");
                                system("mkdir GaussianFit");
                                system("mkdir mW_versus_mTop");
                                system("mv BTagEffPlots BlindedResults_2018/");
                                system("mv Chi2Range BlindedResults_2018/");
                                system("mv CutFlowReport BlindedResults_2018/");
                                system("mv Resolution BlindedResults_2018/");
                                system("mv Results BlindedResults_2018/");
                                system("mv TriggerSF BlindedResults_2018/");
                                system("mv GaussianFit BlindedResults_2018/");
                                system("mv mW_versus_mTop BlindedResults_2018/");
                                system("mv BlindedResults_2018/BTagEffPlots_*.root BlindedResults_2018/BTagEffPlots");
                                system("mv BlindedResults_2018/Chi2Range_*.txt BlindedResults_2018/Chi2Range");
                                system("mv BlindedResults_2018/CutFlowReport_*.txt BlindedResults_2018/CutFlowReport");
                                system("mv BlindedResults_2018/Resolution_*.txt BlindedResults_2018/Resolution");
                                system("mv BlindedResults_2018/Results_*.root BlindedResults_2018/Results");
                                system("mv BlindedResults_2018/TurnOnCurves_*.root BlindedResults_2018/TriggerSF");
                                system("mv BlindedResults_2018/TriggerSF_*.txt BlindedResults_2018/TriggerSF");
				system("mv BlindedResults_2018/*_GaussianFit_*.root BlindedResults_2018/GaussianFit");
                                system("mv BlindedResults_2018/*_mW_mTop_*.root BlindedResults_2018/mW_versus_mTop");		

	
			}
			else{

				system("mkdir BlindedResults_NPL_2018");
                                system("mv *NPL*.root BlindedResults_NPL_2018");
                                system("mv *NPL*.root BlindedResults_NPL_2018");
                                system("mkdir BTagEffPlots");
                                system("mkdir Chi2Range");
                                system("mkdir CutFlowReport");
                                system("mkdir Resolution");
                                system("mkdir Results");
                                system("mkdir TriggerSF");
                                system("mkdir GaussianFit");
                                system("mkdir mW_versus_mTop");
                                system("mv BTagEffPlots BlindedResults_NPL_2018/");
                                system("mv Chi2Range BlindedResults_NPL_2018/");
                                system("mv CutFlowReport BlindedResults_NPL_2018/");
                                system("mv Resolution BlindedResults_NPL_2018/");
                                system("mv Results BlindedResults_NPL_2018/");
                                system("mv TriggerSF BlindedResults_NPL_2018/");
                                system("mv GaussianFit BlindedResults_NPL_2018/");
                                system("mv mW_versus_mTop BlindedResults_NPL_2018/");
                                system("mv UnblindedResults_NPL_2018/BTagEffPlots_*.root BlindedResults_NPL_2018/BTagEffPlots");
                                system("mv BlindedResults_NPL_2018/Chi2Range_*.txt BlindedResults_NPL_2018/Chi2Range");
                                system("mv BlindedResults_NPL_2018/CutFlowReport_*.txt BlindedResults_NPL_2018/CutFlowReport");
                                system("mv BlindedResults_NPL_2018/Resolution_*.txt BlindedResults_NPL_2018/Resolution");
                                system("mv BlindedResults_NPL_2018/Results_*.root BlindedResults_NPL_2018/Results");
                                system("mv BlindedResults_NPL_2018/TurnOnCurves_*.root BlindedResults_NPL_2018/TriggerSF");
                                system("mv BlindedResults_NPL_2018/TriggerSF_*.txt BlindedResults_NPL_2018/TriggerSF");
				system("mv BlindedResults_NPL_2018/*_GaussianFit_*.root BlindedResults_NPL_2018/GaussianFit");
                                system("mv BlindedResults_NPL_2018/*_mW_mTop_*.root BlindedResults_NPL_2018/mW_versus_mTop");
		

			}


                }
                else{

			if(NPL == false){

				system("mkdir UnblindedResults_2018");
                                system("mv *2018*.root UnblindedResults_2018");
                                system("mv *2018*.txt UnblindedResults_2018");
                                system("mkdir BTagEffPlots");
                                system("mkdir Chi2Range");
                                system("mkdir CutFlowReport");
                                system("mkdir Resolution");
                                system("mkdir Results");
                                system("mkdir TriggerSF");
                                system("mkdir GaussianFit");
                                system("mkdir mW_versus_mTop");
                                system("mv BTagEffPlots UnblindedResults_2018/");
                                system("mv Chi2Range UnblindedResults_2018/");
                                system("mv CutFlowReport UnblindedResults_2018/");
                                system("mv Resolution UnblindedResults_2018/");
                                system("mv Results UnblindedResults_2018/");
                                system("mv TriggerSF UnblindedResults_2018/");
                                system("mv GaussianFit UnblindedResults_2018/");
                                system("mv mW_versus_mTop UnblindedResults_2018/");
                                system("mv UnblindedResults_2018/BTagEffPlots_*.root UnblindedResults_2018/BTagEffPlots");
                                system("mv UnblindedResults_2018/Chi2Range_*.txt UnblindedResults_2018/Chi2Range");
                                system("mv UnblindedResults_2018/CutFlowReport_*.txt UnblindedResults_2018/CutFlowReport");
                                system("mv UnblindedResults_2018/Resolution_*.txt UnblindedResults_2018/Resolution");
                                system("mv UnblindedResults_2018/Results_*.root UnblindedResults_2018/Results");
                                system("mv UnblindedResults_2018/TurnOnCurves_*.root UnblindedResults_2018/TriggerSF");
                                system("mv UnblindedResults_2018/TriggerSF_*.txt UnblindedResults_2018/TriggerSF");
				system("mv UnblindedResults_2018/*_GaussianFit_*.root UnblindedResults_2018/GaussianFit");
                                system("mv UnblindedResults_2018/*_mW_mTop_*.root UnblindedResults_2018/mW_versus_mTop");



			}
			else{

				system("mkdir UnblindedResults_NPL_2018");
                                system("mv *NPL*.root UnblindedResults_NPL_2018");
                                system("mv *NPL*.txt UnblindedResults_NPL_2018");
                                system("mkdir BTagEffPlots");
                                system("mkdir Chi2Range");
                                system("mkdir CutFlowReport");
                                system("mkdir Resolution");
                                system("mkdir Results");
                                system("mkdir TriggerSF");
                                system("mkdir GaussianFit");
                                system("mkdir mW_versus_mTop");
                                system("mv BTagEffPlots UnblindedResults_NPL_2018/");
                                system("mv Chi2Range UnblindedResults_NPL_2018/");
                                system("mv CutFlowReport UnblindedResults_NPL_2018/");
                                system("mv Resolution UnblindedResults_NPL_2018/");
                                system("mv Results UnblindedResults_NPL_2018/");
                                system("mv TriggerSF UnblindedResults_NPL_2018/");
                                system("mv GaussianFit UnblindedResults_NPL_2018/");
                                system("mv mW_versus_mTop UnblindedResults_NPL_2018/");
                                system("mv UnblindedResults_NPL_2018/BTagEffPlots_*.root UnblindedResults_NPL_2018/BTagEffPlots");
                                system("mv UnblindedResults_NPL_2018/Chi2Range_*.txt UnblindedResults_NPL_2018/Chi2Range");
                                system("mv UnblindedResults_NPL_2018/CutFlowReport_*.txt UnblindedResults_NPL_2018/CutFlowReport");
                                system("mv UnblindedResults_NPL_2018/Resolution_*.txt UnblindedResults_NPL_2018/Resolution");
                                system("mv UnblindedResults_NPL_2018/Results_*.root UnblindedResults_NPL_2018/Results");
                                system("mv UnblindedResults_NPL_2018/TurnOnCurves_*.root UnblindedResults_NPL_2018/TriggerSF");
                                system("mv UnblindedResults_NPL_2018/TriggerSF_*.txt UnblindedResults_NPL_2018/TriggerSF");
				system("mv UnblindedResults_NPL_2018/*_GaussianFit_*.root UnblindedResults_NPL_2018/GaussianFit");
                                system("mv UnblindedResults_NPL_2018/*_mW_mTop_*.root UnblindedResults_NPL_2018/mW_versus_mTop");


			}

                }

	}
	else{std::cout << "Year must be 2016, 2017 or 2018" << std::endl;}


}



namespace{

constexpr double EndcapMinEta = 1.566;
constexpr double BarrelMaxEta = 1.4442;
double MaxTrackerEta;
double MinElectronPt;
double MinMuonPt;
double MinElectronPtEmu;
double MinMuonPtEmu;
double MaxElectronPt;
double MaxMuonPt;
int NumberOfSimulatedEvents_ee;
int NumberOfSimulatedEvents_mumu;


float W_stddev_ee;
float Top_stddev_ee;
float W_stddev_mumu;
float Top_stddev_mumu;

float Chi2_SR_ee;
float Chi2_SBR_ee;
float Chi2_SR_mumu;
float Chi2_SBR_mumu;


template<typename T>
[[gnu::const]] T select(const T& a, const ints& mask)
{
    return a[mask];
}

[[gnu::const]] auto delta_phi(const float phi1, const float phi2)
{
    return vdt::fast_atan2f(vdt::fast_sinf(phi1 - phi2), vdt::fast_cosf(phi1 - phi2));
}

[[gnu::const]] auto deltaR(const float eta1, const float phi1, const float eta2, const float phi2)
{
    return std::sqrt(std::pow(eta1 - eta2, 2) + std::pow(delta_phi(phi1, phi2), 2));
}




}

template<typename T, typename U>
[[gnu::const]] bool all_equal(const T& t, const U& u)
{
    return t == u;
}

template<typename T, typename U, typename... Types>
[[gnu::const]] bool all_equal(const T& t, const U& u, Types const&... args)
{
    return t == u && all_equal(u, args...);
}


[[gnu::const]] auto inv_mass(const floats& pts, const floats& etas, const floats& phis, const floats& ms)
{



    if (!all_equal(pts.size(), etas.size(), phis.size(), ms.size()))
    {
        throw std::logic_error("Collections must be the same size");
    }
    else if (pts.empty())
    {
        throw std::logic_error("Collections must not be empty");
    }

    TLorentzVector vec{};
    for (size_t i{0}; i < pts.size(); i++)
    {
        TLorentzVector p{};
        p.SetPtEtaPhiM(pts[i], etas[i], phis[i], ms[i]);
        vec += p;
    }
    return boost::numeric_cast<float>(vec.M());
}


[[gnu::const]] auto inv_mass_doubles(const doubles& pts, const doubles& etas, const doubles& phis, const doubles& ms)
{


    if (!all_equal(pts.size(), etas.size(), phis.size(), ms.size()))
    {
        throw std::logic_error("Collections must be the same size");
    }
    else if (pts.empty())
    {
        throw std::logic_error("Collections must not be empty");
    }

    TLorentzVector vec{};
    for (size_t i{0}; i < pts.size(); i++)
    {
        TLorentzVector p{};
        p.SetPtEtaPhiM(pts[i], etas[i], phis[i], ms[i]);
        vec += p;
    }
    return boost::numeric_cast<float>(vec.M());
}



//Reading values from the normalisation text file

double linereader(const int& LineNumber, const std::string& year){


   using namespace std;

   std::string NormFileString = "Normalisation/NormalisationFactors_" + year + ".txt"; 

   std::fstream file(NormFileString.c_str());
   GotoLine(file, LineNumber);

   std::string line;
   file >> line;

   double Value = atof(line.c_str());
   return Value;

}
//end of functions for normalisation










void fulleventselection_calculator(const std::string& process, const bool& blinding, const bool& NPL, const bool& SR, const bool& SBR, const bool& ZPlusJetsCR, const bool& ttbarCR, const std::string& year, const bool& PU_ScaleUp, const bool& PU_ScaleDown, const bool& BTag_ScaleUp, const bool& BTag_ScaleDown, const bool& JetSmearing_ScaleUp, const bool& JetSmearing_ScaleDown, const bool& JetResolution_ScaleUp, const bool& JetResolution_ScaleDown, const bool& LeptonEfficiencies_ScaleUp, const bool& LeptonEfficiencies_ScaleDown, const bool& PDF_ScaleUp, const bool& PDF_ScaleDown, const bool& ME_Up, const bool& ME_Down, const bool& MET_Up, const bool& MET_Down, const bool& isr_up, const bool& isr_down, const bool& fsr_up, const bool& fsr_down){


std::cout << '\n' << std::endl;
std::cout << "The process is: " << process << std::endl;
std::cout << "Blinding is: " << blinding << std::endl;
std::cout << "NPL is: " << NPL << std::endl;
std::cout << "SR is: " << SR << std::endl;
std::cout << "SBR is: " << SBR << std::endl;
std::cout << "ZPlusJetsCR is: " << ZPlusJetsCR << std::endl;
std::cout << "ttbarCR is: " << ttbarCR << std::endl;
std::cout << "The year is: " << year << std::endl;
std::cout << "PU_ScaleUp is: " << PU_ScaleUp << std::endl;
std::cout << "PU_ScaleDown is: " << PU_ScaleDown << std::endl;
std::cout << "BTag_ScaleUp is: " << BTag_ScaleUp << std::endl;
std::cout << "BTag_ScaleDown is: " << BTag_ScaleDown << std::endl;
std::cout << "JetSmearing_ScaleUp is: " << JetSmearing_ScaleUp << std::endl;
std::cout << "JetSmearing_ScaleDown is: " << JetSmearing_ScaleDown << std::endl;
std::cout << "JetResolution_ScaleUp is: " << JetResolution_ScaleUp << std::endl;
std::cout << "JetResolution_ScaleDown is: " << JetResolution_ScaleDown << std::endl;
std::cout << "LeptonEfficiencies_ScaleUp is: " << LeptonEfficiencies_ScaleUp << std::endl;
std::cout << "LeptonEfficiencies_ScaleDown is: " << LeptonEfficiencies_ScaleDown << std::endl;
std::cout << "PDF_ScaleUp is: " << PDF_ScaleUp << std::endl;
std::cout << "PDF_ScaleDown is: " << PDF_ScaleDown << std::endl;
std::cout << "ME_Up is: " << ME_Up << std::endl;
std::cout << "ME_Down is: " << ME_Down << std::endl;
std::cout << "MET_Up is: " << MET_Up << std::endl;
std::cout << "MET_Down is: " << MET_Down << std::endl;
std::cout << "isr_up is: " << isr_up << std::endl;
std::cout << "isr_down is: " << isr_down << std::endl;
std::cout << "fsr_up is: " << fsr_up << std::endl;
std::cout << "fsr_down is: " << fsr_down << std::endl;
std::cout << '\n' << std::endl;


//EnableImplicitMT();

std::vector<std::string> input_files;
std::ofstream CutFlowReport;
std::string cutflowstring;


if(process != "MC_triggerSF_ttbar" && process != "MC_triggerSF_ZPlusJets" && process != "Data_triggerSF"){

	if(blinding == false){

		if(PU_ScaleUp == true){

			if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_PU_ScaleUp_" + year + "_NPL.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_PU_ScaleUp_" + year + "_ZPlusJetsCR.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_PU_ScaleUp_" + year + "_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_PU_ScaleUp_" + year + "_NPL_ZPlusJetsCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_PU_ScaleUp_" + year + "_NPL_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_PU_ScaleUp_" + year + ".txt";}	



		}
		else if(PU_ScaleDown == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_PU_ScaleDown_" + year + "_NPL.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_PU_ScaleDown_" + year + "_ZPlusJetsCR.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_PU_ScaleDown_" + year + "_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_PU_ScaleDown_" + year + "_NPL_ZPlusJetsCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_PU_ScaleDown_" + year + "_NPL_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_PU_ScaleDown_" + year + ".txt";}



                }
		else if(BTag_ScaleUp == true){
                
                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_BTag_ScaleUp_" + year + "_NPL.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_BTag_ScaleUp_" + year + "_ZPlusJetsCR.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_BTag_ScaleUp_" + year + "_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_BTag_ScaleUp_" + year + "_NPL_ZPlusJetsCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_BTag_ScaleUp_" + year + "_NPL_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_BTag_ScaleUp_" + year + ".txt";}


                }
		else if(BTag_ScaleDown == true){
                
                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_BTag_ScaleDown_" + year + "_NPL.txt";
                        }       
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_BTag_ScaleDown_" + year + "_ZPlusJetsCR.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_BTag_ScaleDown_" + year + "_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_BTag_ScaleDown_" + year + "_NPL_ZPlusJetsCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_BTag_ScaleDown_" + year + "_NPL_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_BTag_ScaleDown_" + year + ".txt";}


                }
		else if(JetSmearing_ScaleUp == true){
                
                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_JetSmearing_ScaleUp_" + year + "_NPL.txt";
                        }       
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_JetSmearing_ScaleUp_" + year + "_ZPlusJetsCR.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_JetSmearing_ScaleUp_" + year + "_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_JetSmearing_ScaleUp_" + year + "_NPL_ZPlusJetsCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_JetSmearing_ScaleUp_" + year + "_NPL_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_JetSmearing_ScaleUp_" + year + ".txt";}


                }
		else if(JetSmearing_ScaleDown == true){
                
                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_JetSmearing_ScaleDown_" + year + "_NPL.txt";
                        }       
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_JetSmearing_ScaleDown_" + year + "_ZPlusJetsCR.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_JetSmearing_ScaleDown_" + year + "_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_JetSmearing_ScaleDown_" + year + "_NPL_ZPlusJetsCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_JetSmearing_ScaleDown_" + year + "_NPL_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_JetSmearing_ScaleDown_" + year + ".txt";}


                }
		else if(JetResolution_ScaleUp == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_JetResolution_ScaleUp_" + year + "_NPL.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_JetResolution_ScaleUp_" + year + "_ZPlusJetsCR.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_JetResolution_ScaleUp_" + year + "_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_JetResolution_ScaleUp_" + year + "_NPL_ZPlusJetsCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_JetResolution_ScaleUp_" + year + "_NPL_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_JetResolution_ScaleUp_" + year + ".txt";}


                }
                else if(JetResolution_ScaleDown == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_JetResolution_ScaleDown_" + year + "_NPL.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_JetResolution_ScaleDown_" + year + "_ZPlusJetsCR.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_JetResolution_ScaleDown_" + year + "_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_JetResolution_ScaleDown_" + year + "_NPL_ZPlusJetsCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_JetResolution_ScaleDown_" + year + "_NPL_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_JetResolution_ScaleDown_" + year + ".txt";}


                }
		else if(LeptonEfficiencies_ScaleUp == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_LeptonEfficiencies_ScaleUp_" + year + "_NPL.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_LeptonEfficiencies_ScaleUp_" + year + "_ZPlusJetsCR.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_LeptonEfficiencies_ScaleUp_" + year + "_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_LeptonEfficiencies_ScaleUp_" + year + "_NPL_ZPlusJetsCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_LeptonEfficiencies_ScaleUp_" + year + "_NPL_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_LeptonEfficiencies_ScaleUp_" + year + ".txt";}


                }
                else if(LeptonEfficiencies_ScaleDown == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_LeptonEfficiencies_ScaleDown_" + year + "_NPL.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_LeptonEfficiencies_ScaleDown_" + year + "_ZPlusJetsCR.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_LeptonEfficiencies_ScaleDown_" + year + "_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_LeptonEfficiencies_ScaleDown_" + year + "_NPL_ZPlusJetsCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_LeptonEfficiencies_ScaleDown_" + year + "_NPL_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_LeptonEfficiencies_ScaleDown_" + year + ".txt";}


                }
		else if(PDF_ScaleUp == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_PDF_ScaleUp_" + year + "_NPL.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_PDF_ScaleUp_" + year + "_ZPlusJetsCR.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_PDF_ScaleUp_" + year + "_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_PDF_ScaleUp_" + year + "_NPL_ZPlusJetsCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_PDF_ScaleUp_" + year + "_NPL_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_PDF_ScaleUp_" + year + ".txt";}


                }
                else if(PDF_ScaleDown == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_PDF_ScaleDown_" + year + "_NPL.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_PDF_ScaleDown_" + year + "_ZPlusJetsCR.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_PDF_ScaleDown_" + year + "_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_PDF_ScaleDown_" + year + "_NPL_ZPlusJetsCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_PDF_ScaleDown_" + year + "_NPL_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_PDF_ScaleDown_" + year + ".txt";}


                }	
		else if(ME_Up == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_ME_Up_" + year + "_NPL.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_ME_Up_" + year + "_ZPlusJetsCR.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_ME_Up_" + year + "_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_ME_Up_" + year + "_NPL_ZPlusJetsCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_ME_Up_" + year + "_NPL_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_ME_Up_" + year + ".txt";}


                }
                else if(ME_Down == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_ME_Down_" + year + "_NPL.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_ME_Down_" + year + "_ZPlusJetsCR.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_ME_Down_" + year + "_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_ME_Down_" + year + "_NPL_ZPlusJetsCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_ME_Down_" + year + "_NPL_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_ME_Down_" + year + ".txt";}


                }
		else if(MET_Up == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_MET_Up_" + year + "_NPL.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_MET_Up_" + year + "_ZPlusJetsCR.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_MET_Up_" + year + "_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_MET_Up_" + year + "_NPL_ZPlusJetsCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_MET_Up_" + year + "_NPL_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_MET_Up_" + year + ".txt";}


                }
                else if(MET_Down == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_MET_Down_" + year + "_NPL.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_MET_Down_" + year + "_ZPlusJetsCR.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_MET_Down_" + year + "_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_MET_Down_" + year + "_NPL_ZPlusJetsCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_MET_Down_" + year + "_NPL_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_MET_Down_" + year + ".txt";}


                }
		else if(isr_up == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_isr_up_" + year + "_NPL.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_isr_up_" + year + "_ZPlusJetsCR.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_isr_up_" + year + "_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_isr_up_" + year + "_NPL_ZPlusJetsCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_isr_up_" + year + "_NPL_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_isr_up_" + year + ".txt";}


                }
                else if(isr_down == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_isr_down_" + year + "_NPL.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_isr_down_" + year + "_ZPlusJetsCR.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_isr_down_" + year + "_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_isr_down_" + year + "_NPL_ZPlusJetsCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_isr_down_" + year + "_NPL_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_isr_down_" + year + ".txt";}


                }
		else if(fsr_up == true){
                
                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_fsr_up_" + year + "_NPL.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_fsr_up_" + year + "_ZPlusJetsCR.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_fsr_up_" + year + "_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_fsr_up_" + year + "_NPL_ZPlusJetsCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_fsr_up_" + year + "_NPL_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_fsr_up_" + year + ".txt";}


                }
                else if(fsr_down == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_fsr_down_" + year + "_NPL.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_fsr_down_" + year + "_ZPlusJetsCR.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_fsr_down_" + year + "_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_fsr_down_" + year + "_NPL_ZPlusJetsCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_fsr_down_" + year + "_NPL_ttbarCR.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_fsr_down_" + year + ".txt";}


                }	
		else{

			if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
        			cutflowstring = "CutFlowReport_" + process + "_" + year + "_NPL.txt";
			}
			else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
				cutflowstring = "CutFlowReport_" + process + "_" + year + "_ZPlusJetsCR.txt";
			}
			else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
				cutflowstring = "CutFlowReport_" + process + "_" + year + "_ttbarCR.txt";
			}
			else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
				cutflowstring = "CutFlowReport_" + process + "_" + year + "_NPL_ZPlusJetsCR.txt";
			}
			else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
				cutflowstring = "CutFlowReport_" + process + "_" + year + "_NPL_ttbarCR.txt";
			}
			else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
			else{cutflowstring = "CutFlowReport_" + process + "_" + year + ".txt";}

		}

	}
	else{


		if(PU_ScaleUp == true){

			if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_PU_ScaleUp_" + year + "_NPL_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_PU_ScaleUp_" + year + "_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_PU_ScaleUp_" + year + "_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_PU_ScaleUp_" + year + "_NPL_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_PU_ScaleUp_" + year + "_NPL_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_PU_ScaleUp_" + year + "_Blinded.txt";}	



		}
		else if(PU_ScaleDown == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_PU_ScaleDown_" + year + "_NPL_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_PU_ScaleDown_" + year + "_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_PU_ScaleDown_" + year + "_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_PU_ScaleDown_" + year + "_NPL_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_PU_ScaleDown_" + year + "_NPL_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_PU_ScaleDown_" + year + "_Blinded.txt";}



                }
		else if(BTag_ScaleUp == true){
                
                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_BTag_ScaleUp_" + year + "_NPL_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_BTag_ScaleUp_" + year + "_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_BTag_ScaleUp_" + year + "_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_BTag_ScaleUp_" + year + "_NPL_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_BTag_ScaleUp_" + year + "_NPL_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_BTag_ScaleUp_" + year + "_Blinded.txt";}


                }
		else if(BTag_ScaleDown == true){
                
                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_BTag_ScaleDown_" + year + "_NPL_Blinded.txt";
                        }       
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_BTag_ScaleDown_" + year + "_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_BTag_ScaleDown_" + year + "_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_BTag_ScaleDown_" + year + "_NPL_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_BTag_ScaleDown_" + year + "_NPL_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_BTag_ScaleDown_" + year + "_Blinded.txt";}


                }
		else if(JetSmearing_ScaleUp == true){
                
                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_JetSmearing_ScaleUp_" + year + "_NPL_Blinded.txt";
                        }       
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_JetSmearing_ScaleUp_" + year + "_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_JetSmearing_ScaleUp_" + year + "_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_JetSmearing_ScaleUp_" + year + "_NPL_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_JetSmearing_ScaleUp_" + year + "_NPL_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_JetSmearing_ScaleUp_" + year + "_Blinded.txt";}


                }
		else if(JetSmearing_ScaleDown == true){
                
                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_JetSmearing_ScaleDown_" + year + "_NPL_Blinded.txt";
                        }       
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_JetSmearing_ScaleDown_" + year + "_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_JetSmearing_ScaleDown_" + year + "_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_JetSmearing_ScaleDown_" + year + "_NPL_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_JetSmearing_ScaleDown_" + year + "_NPL_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_JetSmearing_ScaleDown_" + year + "_Blinded.txt";}


                }
		else if(JetResolution_ScaleUp == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_JetResolution_ScaleUp_" + year + "_NPL_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_JetResolution_ScaleUp_" + year + "_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_JetResolution_ScaleUp_" + year + "_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_JetResolution_ScaleUp_" + year + "_NPL_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_JetResolution_ScaleUp_" + year + "_NPL_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_JetResolution_ScaleUp_" + year + "_Blinded.txt";}


                }
                else if(JetResolution_ScaleDown == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_JetResolution_ScaleDown_" + year + "_NPL_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_JetResolution_ScaleDown_" + year + "_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_JetResolution_ScaleDown_" + year + "_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_JetResolution_ScaleDown_" + year + "_NPL_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_JetResolution_ScaleDown_" + year + "_NPL_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_JetResolution_ScaleDown_" + year + "_Blinded.txt";}


                }
		else if(LeptonEfficiencies_ScaleUp == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_LeptonEfficiencies_ScaleUp_" + year + "_NPL_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_LeptonEfficiencies_ScaleUp_" + year + "_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_LeptonEfficiencies_ScaleUp_" + year + "_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_LeptonEfficiencies_ScaleUp_" + year + "_NPL_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_LeptonEfficiencies_ScaleUp_" + year + "_NPL_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_LeptonEfficiencies_ScaleUp_" + year + "_Blinded.txt";}


                }
                else if(LeptonEfficiencies_ScaleDown == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_LeptonEfficiencies_ScaleDown_" + year + "_NPL_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_LeptonEfficiencies_ScaleDown_" + year + "_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_LeptonEfficiencies_ScaleDown_" + year + "_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_LeptonEfficiencies_ScaleDown_" + year + "_NPL_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_LeptonEfficiencies_ScaleDown_" + year + "_NPL_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_LeptonEfficiencies_ScaleDown_" + year + "_Blinded.txt";}


                }
		else if(PDF_ScaleUp == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_PDF_ScaleUp_" + year + "_NPL_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_PDF_ScaleUp_" + year + "_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_PDF_ScaleUp_" + year + "_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_PDF_ScaleUp_" + year + "_NPL_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_PDF_ScaleUp_" + year + "_NPL_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_PDF_ScaleUp_" + year + "_Blinded.txt";}


                }
                else if(PDF_ScaleDown == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_PDF_ScaleDown_" + year + "_NPL_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_PDF_ScaleDown_" + year + "_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_PDF_ScaleDown_" + year + "_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_PDF_ScaleDown_" + year + "_NPL_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_PDF_ScaleDown_" + year + "_NPL_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_PDF_ScaleDown_" + year + "_Blinded.txt";}


                }	
		else if(ME_Up == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_ME_Up_" + year + "_NPL_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_ME_Up_" + year + "_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_ME_Up_" + year + "_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_ME_Up_" + year + "_NPL_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_ME_Up_" + year + "_NPL_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_ME_Up_" + year + "Blinded.txt";}


                }
                else if(ME_Down == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_ME_Down_" + year + "_NPL_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_ME_Down_" + year + "_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_ME_Down_" + year + "_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_ME_Down_" + year + "_NPL_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_ME_Down_" + year + "_NPL_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_ME_Down_" + year + "_Blinded.txt";}


                }
		else if(MET_Up == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_MET_Up_" + year + "_NPL_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_MET_Up_" + year + "_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_MET_Up_" + year + "_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_MET_Up_" + year + "_NPL_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_MET_Up_" + year + "_NPL_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_MET_Up_" + year + "Blinded.txt";}


                }
                else if(MET_Down == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_MET_Down_" + year + "_NPL_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_MET_Down_" + year + "_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_MET_Down_" + year + "_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_MET_Down_" + year + "_NPL_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_MET_Down_" + year + "_NPL_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_MET_Down_" + year + "_Blinded.txt";}


                }
		else if(isr_up == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_isr_up_" + year + "_NPL_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_isr_up_" + year + "_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_isr_up_" + year + "_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_isr_up_" + year + "_NPL_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_isr_up_" + year + "_NPL_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_isr_up_" + year + "_Blinded.txt";}


                }
                else if(isr_down == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_isr_down_" + year + "_NPL_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_isr_down_" + year + "_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_isr_down_" + year + "_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_isr_down_" + year + "_NPL_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_isr_down_" + year + "_NPL_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_isr_down_" + year + "_Blinded.txt";}


                }
		else if(fsr_up == true){
                
                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_fsr_up_" + year + "_NPL_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_fsr_up_" + year + "_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_fsr_up_" + year + "_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_fsr_up_" + year + "_NPL_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_fsr_up_" + year + "_NPL_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_fsr_up_" + year + "_Blinded.txt";}


                }
                else if(fsr_down == true){

                        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_fsr_down_" + year + "_NPL_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_fsr_down_" + year + "_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_fsr_down_" + year + "_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                                cutflowstring = "CutFlowReport_" + process + "_fsr_down_" + year + "_NPL_ZPlusJetsCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                                cutflowstring = "CutFlowReport_" + process + "_fsr_down_" + year + "_NPL_ttbarCR_Blinded.txt";
                        }
                        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                        else{cutflowstring = "CutFlowReport_" + process + "_fsr_down_" + year + "_Blinded.txt";}


                }		
		else{

			if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                		cutflowstring = "CutFlowReport_" + process + "_" + year + "_NPL_Blinded.txt";
        		}
        		else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                		cutflowstring = "CutFlowReport_" + process + "_" + year + "_ZPlusJetsCR_Blinded.txt";
        		}
        		else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                		cutflowstring = "CutFlowReport_" + process + "_" + year + "_ttbarCR_Blinded.txt";
        		}
        		else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                		cutflowstring = "CutFlowReport_" + process + "_" + year + "_NPL_ZPlusJetsCR_Blinded.txt";
        		}
        		else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                		cutflowstring = "CutFlowReport_" + process + "_" + year + "_NPL_ttbarCR_Blinded.txt";
        		}
        		else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
        		else{cutflowstring = "CutFlowReport_" + process + "_" + year + "_Blinded.txt";}
		
		}


	}


	CutFlowReport.open(cutflowstring.c_str());


}


if(year == "2016"){

	if(process == "tZq"){input_files = {"/data/disk2/nanoAOD_2016/tZq_ll/*"};}
	else if(process == "ZPlusJets_M50_aMCatNLO"){input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_M50_aMCatNLO/*"};}
	else if(process == "ZPlusJets_M10To50_aMCatNLO"){input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_M10ToM50_aMCatNLO/*"};}
	else if(process == "ZPlusJets_M10To50_aMCatNLO_ext"){input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_M10ToM50_ext_aMCatNLO/*"};}
	else if(process == "ZPlusJets_M50_Madgraph"){input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_M50_Madgraph/*"};}
	else if(process == "ZPlusJets_M50_Madgraph_ext"){input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_M50_Madgraph_ext/*"};}
	else if(process == "ZPlusJets_M10To50_Madgraph"){input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_M10To50_Madgraph/*"};}
	else if(process == "ZPlusJets_PtBinned_0To50"){input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_0To50/*"};}
	else if(process == "ZPlusJets_PtBinned_50To100"){input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_50To100/*"};}
	else if(process == "ZPlusJets_PtBinned_50To100_ext"){input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_50To100_ext/*"};}
	else if(process == "ZPlusJets_PtBinned_100To250"){input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_100To250/*"};}
	else if(process == "ZPlusJets_PtBinned_100To250_ext1"){input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_100To250_ext1/*"};}
	else if(process == "ZPlusJets_PtBinned_100To250_ext2"){input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_100To250_ext2/*"};}
	else if(process == "ZPlusJets_PtBinned_100To250_ext5"){input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_100To250_ext5/*"};}
	else if(process == "ZPlusJets_PtBinned_250To400"){input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_250To400/*"};}
        else if(process == "ZPlusJets_PtBinned_250To400_ext1"){input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_250To400_ext1/*"};}
        else if(process == "ZPlusJets_PtBinned_250To400_ext2"){input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_250To400_ext2/*"};}
        else if(process == "ZPlusJets_PtBinned_400To650"){input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_400To650/*"};}
        else if(process == "ZPlusJets_PtBinned_400To650_ext1"){input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_400To650_ext1/*"};}
        else if(process == "ZPlusJets_PtBinned_400To650_ext2"){input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_400To650_ext2/*"};}
	else if(process == "ZPlusJets_PtBinned_650ToInf"){input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_650ToInf/*"};}
        else if(process == "ZPlusJets_PtBinned_650ToInf_ext1"){input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_650ToInf_ext1/*"};}
        else if(process == "ZPlusJets_PtBinned_650ToInf_ext2"){input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_650ToInf_ext2/*"};}
	else if(process == "SingleTop_tchannel_top"){input_files = {"/data/disk2/nanoAOD_2016/ST_tchannel_top/*"};}
	else if(process == "SingleTop_tchannel_tbar"){input_files = {"/data/disk2/nanoAOD_2016/ST_tchannel_antitop/*"};}
	else if(process == "SingleTop_schannel"){input_files = {"/data/disk2/nanoAOD_2016/ST_schannel/*"};}
	else if(process == "SingleTop_tW"){input_files = {"/data/disk2/nanoAOD_2016/tW_top/*"};}
	else if(process == "SingleTop_tbarW"){input_files = {"/data/disk2/nanoAOD_2016/tW_tbar/*"};}
	else if(process == "SingleTop_tHq"){input_files = {"/data/disk2/nanoAOD_2016/tHq/*"};}
	else if(process == "SingleTop_tWZ_tWll"){input_files = {"/data/disk2/nanoAOD_2016/tWZ_tWLL/*"};} 
	else if(process == "ttbar_inc"){input_files = {"/data/disk2/nanoAOD_2016/ttbar_inc/*"};}
	else if(process == "ttbar_madgraph_NanoAODv5"){input_files = {"/data/disk2/nanoAOD_2016/ttbar_madgraph/*"};}
	else if(process == "ttbar_aMCatNLO"){input_files = {"/data/disk2/nanoAOD_2016/ttbar_aMCatNLO/*"};}
	else if(process == "VV_ZZTo2L2Nu"){input_files = {"/data/disk2/nanoAOD_2016/ZZTo2L2Nu/*"};}
	else if(process == "VV_ZZTo2L2Q"){input_files = {"/data/disk2/nanoAOD_2016/ZZTo2L2Q/*"};}
	else if(process == "VV_ZZTo4L"){input_files = {"/data/disk2/nanoAOD_2016/ZZTo4L/*"};}
	else if(process == "VV_WZTo1L1Nu2Q"){input_files = {"/data/disk2/nanoAOD_2016/WZTo1L1Nu2Q/*"};}
	else if(process == "VV_WZTo2L2Q"){input_files = {"/data/disk2/nanoAOD_2016/WZTo2L2Q/*"};}
	else if(process == "VV_WWTo1L1Nu2Q"){input_files = {"/data/disk2/nanoAOD_2016/WWTo1L1Nu2Q/*"};}
	else if(process == "VV_WWTo2L2Nu"){input_files = {"/data/disk2/nanoAOD_2016/WWTo2L2Nu/*"};}
	else if(process == "VV_WWToLNuQQ"){input_files = {"/data/disk2/nanoAOD_2016/WWToLNuQQ/*"};}
	else if(process == "VV_WGToLNuG"){input_files = {"/data/disk2/nanoAOD_2016/WGToLNuG/*"};}
	else if(process == "VV_ZGToLLG"){input_files = {"/data/disk2/nanoAOD_2016/ZGToLLG/*"};}
	else if(process == "VVV_WWWTo4F"){input_files = {"/data/disk2/nanoAOD_2016/WWWTo4F/*"};}
	else if(process == "VVV_WWZ"){input_files = {"/data/disk2/nanoAOD_2016/WWZ/*"};}
	else if(process == "VVV_WZZ"){input_files = {"/data/disk2/nanoAOD_2016/WZZ/*"};}
	else if(process == "VVV_ZZZ"){input_files = {"/data/disk2/nanoAOD_2016/ZZZ/*"};}
	else if(process == "WPlusJets_WJetsToLNu"){input_files = {"/data/disk2/nanoAOD_2016/WJetsToLNu/*"};}
	else if(process == "ttbarV_ttWJetsToLNu"){input_files = {"/data/disk2/nanoAOD_2016/ttWJetsToLNu/*"};}
	else if(process == "ttbarV_ttWJetsToQQ"){input_files = {"/data/disk2/nanoAOD_2016/ttWJetsToQQ/*"};}
	else if(process == "ttbarV_ttgamma"){input_files = {"/data/disk2/nanoAOD_2016/ttgamma/*"};}
	else if(process == "ttbarV_ttZToQQ"){input_files = {"/data/disk2/nanoAOD_2016/ttZToQQ/*"};}
	else if(process == "ttbarV_ttHTobb"){input_files = {"/data/disk2/nanoAOD_2016/ttHTobb/*"};}
	else if(process == "ttbarV_ttHToNonbb"){input_files = {"/data/disk2/nanoAOD_2016/ttHToNonbb/*"};}
	else if(process == "ttbarV_ttZToLLNuNu"){input_files = {"/data/disk2/nanoAOD_2016/ttZToLLNuNu/*"};}
 	else if(process == "ttbarV_ttZToLLNuNu_ext2"){input_files = {"/data/disk2/nanoAOD_2016/ttZToLLNuNu_ext2/*"};}
	else if(process == "ttbarV_ttZToLLNuNu_ext3"){input_files = {"/data/disk2/nanoAOD_2016/ttZToLLNuNu_ext3/*"};}
	else if(process == "ttbarV_ttZToQQ"){input_files = {"/data/disk2/nanoAOD_2016/ttZToQQ/*"};}
	else if(process == "ttbarV_ttZToQQ_ext"){input_files = {"/data/disk2/nanoAOD_2016/ttZToQQ_ext/*"};}
	else if(process == "TT_hdampUP"){input_files = {"/data/disk2/nanoAOD_2016/TT_hdampUP/*"};}
	else if(process == "TT_hdampUP_ext"){input_files = {"/data/disk2/nanoAOD_2016/TT_hdampUP_ext/*"};}
	else if(process == "TT_hdampDOWN"){input_files = {"/data/disk2/nanoAOD_2016/TT_hdampDOWN/*"};}
        else if(process == "TT_hdampDOWN_ext"){input_files = {"/data/disk2/nanoAOD_2016/TT_hdampDOWN_ext/*"};}
	else if(process == "ST_tchannel_top_hdampup"){input_files = {"/data/disk2/nanoAOD_2016/ST_tchannel_top_hdampup/*"};}
	else if(process == "ST_tchannel_top_hdampdown"){input_files = {"/data/disk2/nanoAOD_2016/ST_tchannel_top_hdampdown/*"};}
	else if(process == "ST_tchannel_top_ScaleUp"){input_files = {"/data/disk2/nanoAOD_2016/ST_tchannel_top_ScaleUp/*"};}
	else if(process == "ST_tchannel_top_ScaleDown"){input_files = {"/data/disk2/nanoAOD_2016/ST_tchannel_top_ScaleDown_NanoAODv6/*"};}
	else if(process == "tW_tbar_ScaleUp"){input_files = {"/data/disk2/nanoAOD_2016/tW_tbar_ScaleUp/*"};}
	else if(process == "tW_tbar_ScaleDown"){input_files = {"/data/disk2/nanoAOD_2016/tW_tbar_ScaleDown/*"};}
	else if(process == "tW_top_ScaleUp"){input_files = {"/data/disk2/nanoAOD_2016/tW_top_ScaleUp/*"};}
        else if(process == "tW_top_ScaleDown"){input_files = {"/data/disk2/nanoAOD_2016/tW_top_ScaleDown/*"};}
	else if(process == "TT_isr_UP"){input_files = {"/data/disk2/nanoAOD_2016/TT_isr_UP/*"};}
        else if(process == "TT_isr_DOWN"){input_files = {"/data/disk2/nanoAOD_2016/TT_isr_DOWN/*"};}
        else if(process == "TT_isr_DOWN_ext"){input_files = {"/data/disk2/nanoAOD_2016/TT_isr_DOWN_ext/*"};}
	else if(process == "TT_fsr_UP"){input_files = {"/data/disk2/nanoAOD_2016/TT_fsr_UP/*"};}
	else if(process == "TT_fsr_UP_ext"){input_files = {"/data/disk2/nanoAOD_2016/TT_fsr_UP_ext/*"};}
	else if(process == "TT_fsr_DOWN"){input_files = {"/data/disk2/nanoAOD_2016/TT_fsr_DOWN/*"};}
        else if(process == "TT_fsr_DOWN_ext"){input_files = {"/data/disk2/nanoAOD_2016/TT_fsr_DOWN_ext/*"};}
	else if(process == "tW_tbar_ScaleUp"){input_files = {"/data/disk2/nanoAOD_2016/tW_tbar_ScaleUp/*"};}
	else if(process == "tW_tbar_ScaleDown"){input_files = {"/data/disk2/nanoAOD_2016/tW_tbar_ScaleDown/*"};}
	else if(process == "tW_top_ScaleUp"){input_files = {"/data/disk2/nanoAOD_2016/tW_top_ScaleUp/*"};}
        else if(process == "tW_top_ScaleDown"){input_files = {"/data/disk2/nanoAOD_2016/tW_top_ScaleDown/*"};}
	else if(process == "data_DoubleEGRunB"){input_files = {"/data/disk2/nanoAOD_2016/DoubleEGRun2016B/*"};}  
	else if(process == "data_DoubleEGRunC"){input_files = {"/data/disk2/nanoAOD_2016/DoubleEGRun2016C/*"};}
	else if(process == "data_DoubleEGRunD"){input_files = {"/data/disk2/nanoAOD_2016/DoubleEGRun2016D/*"};}
	else if(process == "data_DoubleEGRunE"){input_files = {"/data/disk2/nanoAOD_2016/DoubleEGRun2016E/*"};}
	else if(process == "data_DoubleEGRunF"){input_files = {"/data/disk2/nanoAOD_2016/DoubleEGRun2016F/*"};}
	else if(process == "data_DoubleEGRunG"){input_files = {"/data/disk2/nanoAOD_2016/DoubleEGRun2016G/*"};}
        else if(process == "data_DoubleEGRunH"){input_files = {"/data/disk2/nanoAOD_2016/DoubleEGRun2016H/*"};}
	else if(process == "data_DoubleMuonRunB"){input_files = {"/data/disk2/nanoAOD_2016/DoubleMuonRun2016B/*"};}
	else if(process == "data_DoubleMuonRunC"){input_files = {"/data/disk2/nanoAOD_2016/DoubleMuonRun2016C/*"};}
	else if(process == "data_DoubleMuonRunD"){input_files = {"/data/disk2/nanoAOD_2016/DoubleMuonRun2016D/*"};}
	else if(process == "data_DoubleMuonRunE"){input_files = {"/data/disk2/nanoAOD_2016/DoubleMuonRun2016E/*"};}
	else if(process == "data_DoubleMuonRunF"){input_files = {"/data/disk2/nanoAOD_2016/DoubleMuonRun2016F/*"};}
        else if(process == "data_DoubleMuonRunG"){input_files = {"/data/disk2/nanoAOD_2016/DoubleMuonRun2016G/*"};}
        else if(process == "data_DoubleMuonRunH"){input_files = {"/data/disk2/nanoAOD_2016/DoubleMuonRun2016H/*"};}
	else if(process == "data_SingleMuonRunB"){input_files = {"/data/disk3/nanoAOD_2016/SingleMuon_NanoAOD25Oct2019_RunB/*"};}
        else if(process == "data_SingleMuonRunC"){input_files = {"/data/disk3/nanoAOD_2016/SingleMuon_NanoAOD25Oct2019_RunC/*"};}
        else if(process == "data_SingleMuonRunD"){input_files = {"/data/disk3/nanoAOD_2016/SingleMuon_NanoAOD25Oct2019_RunD/*"};}
        else if(process == "data_SingleMuonRunE"){input_files = {"/data/disk3/nanoAOD_2016/SingleMuon_NanoAOD25Oct2019_RunE/*"};}
        else if(process == "data_SingleMuonRunF"){input_files = {"/data/disk3/nanoAOD_2016/SingleMuon_NanoAOD25Oct2019_RunF/*"};}
	else if(process == "data_SingleMuonRunG"){input_files = {"/data/disk3/nanoAOD_2016/SingleMuon_NanoAOD25Oct2019_RunG/*"};}
        else if(process == "data_SingleMuonRunH"){input_files = {"/data/disk3/nanoAOD_2016/SingleMuon_NanoAOD25Oct2019_RunH/*"};}
        else if(process == "data_SingleElectronRunB"){input_files = {"/data/disk3/nanoAOD_2016/SingleElectron_NanoAOD25Oct2019_RunB/*"};}
        else if(process == "data_SingleElectronRunC"){input_files = {"/data/disk3/nanoAOD_2016/SingleElectron_NanoAOD25Oct2019_RunC/*"};}
        else if(process == "data_SingleElectronRunD"){input_files = {"/data/disk3/nanoAOD_2016/SingleElectron_NanoAOD25Oct2019_RunD/*"};}
        else if(process == "data_SingleElectronRunE"){input_files = {"/data/disk3/nanoAOD_2016/SingleElectron_NanoAOD25Oct2019_RunE/*"};}
        else if(process == "data_SingleElectronRunF"){input_files = {"/data/disk3/nanoAOD_2016/SingleElectron_NanoAOD25Oct2019_RunF/*"};}
	else if(process == "data_SingleElectronRunG"){input_files = {"/data/disk3/nanoAOD_2016/SingleElectron_NanoAOD25Oct2019_RunG/*"};}
        else if(process == "data_SingleElectronRunH"){input_files = {"/data/disk3/nanoAOD_2016/SingleElectron_NanoAOD25Oct2019_RunH/*"};}
	else if(process == "Data_triggerSF"){input_files = {"/data/disk2/nanoAOD_2016/METRun2016B/*.root", "/data/disk2/nanoAOD_2016/METRun2016C/*.root", "/data/disk2/nanoAOD_2016/METRun2016D/*.root", "/data/disk2/nanoAOD_2016/METRun2016E/*.root", "/data/disk2/nanoAOD_2016/METRun2016F/*.root", "/data/disk2/nanoAOD_2016/METRun2016G/*.root", "/data/disk2/nanoAOD_2016/METRun2016H/*.root"};}
	else if(process == "MC_triggerSF_ttbar"){input_files = {"/data/disk2/nanoAOD_2016/ttbar_inc/*.root"};}
	else if(process == "MC_triggerSF_ZPlusJets"){input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets*/*"};}
	else if(process == "NPL_File_ee_Blinded"){input_files = {"NPL_ee_output_2016_Blinded.root"};}
	else if(process == "NPL_File_mumu_Blinded"){input_files = {"NPL_mumu_output_2016_Blinded.root"};}
	else if(process == "NPL_File_ee_Unblinded"){input_files = {"NPL_ee_output_2016.root"};}
        else if(process == "NPL_File_mumu_Unblinded"){input_files = {"NPL_mumu_output_2016.root"};}
	else{std::cout << "You inputted the process: " << process << " for the year " << year << ". Please input an MC signal, background of dataset name." << std::endl;}


}
else if(year == "2017"){

	if(process == "tZq"){input_files = {"/data/disk0/nanoAOD_2017/tZq_ll/*"};}
	else if(process == "ZPlusJets_M50_aMCatNLO"){input_files = {"/data/disk0/nanoAOD_2017/DYJetsToLL_NanoAODv5/*"};}
	else if(process == "ZPlusJets_M50_aMCatNLO_ext"){input_files = {"/data/disk0/nanoAOD_2017/DYJetsToLL_ext_NanoAODv5/*"};}
	else if(process == "ZPlusJets_M10To50_Madgraph"){input_files = {"/data/disk0/nanoAOD_2017/DYJetsToLL_M10to50/*"};}
	else if(process == "SingleTop_tchannel_top"){input_files = {"/data/disk0/nanoAOD_2017/ST_tchannel_top/*"};}
	else if(process == "SingleTop_tchannel_tbar"){input_files = {"/data/disk0/nanoAOD_2017/ST_tchannel_tbar/*"};}
	else if(process == "SingleTop_schannel"){input_files = {"/data/disk0/nanoAOD_2017/ST_schannel/*"};}
	else if(process == "SingleTop_tW"){input_files = {"/data/disk0/nanoAOD_2017/ST_tW/*"};}
	else if(process == "SingleTop_tbarW"){input_files = {"/data/disk0/nanoAOD_2017/ST_tbarW/*"};}
	else if(process == "SingleTop_tHq"){input_files = {"/data/disk0/nanoAOD_2017/tHq/*"};}
	else if(process == "SingleTop_tZq_W_lept_Z_had"){input_files = {"/data/disk0/nanoAOD_2017/tZq_W_lept_Z_had/*"};} 
	else if(process == "SingleTop_tWZ_tWll"){input_files = {"/data/disk0/nanoAOD_2017/tWZ_tWLL_NanoAODv5/*"};} 
	else if(process == "ttbar_2l2nu"){input_files = {"/data/disk0/nanoAOD_2017/ttbar_2l2nu_NanoAODv5/*"};}
	else if(process == "ttbar_madgraph_NanoAODv5"){input_files = {"/data/disk0/nanoAOD_2017/ttbar_madgraph_NanoAODv5/*"};}
	else if(process == "ttbar_TTToHadronic"){input_files = {"/data/disk0/nanoAOD_2017/TTToHadronic/*"};}
	else if(process == "ttbar_TTToSemileptonic"){input_files = {"/data/disk0/nanoAOD_2017/TTToSemileptonic/*"};}
	else if(process == "ttbar_aMCatNLO"){input_files = {"/data/disk0/nanoAOD_2017/ttbar_aMCatNLO/*"};}
	else if(process == "VV_ZZTo2Q2Nu"){input_files = {"/data/disk0/nanoAOD_2017/ZZTo2Q2Nu/*"};}
	else if(process == "VV_ZZTo2L2Nu"){input_files = {"/data/disk0/nanoAOD_2017/ZZTo2L2Nu/*"};}
	else if(process == "VV_ZZTo2L2Q"){input_files = {"/data/disk0/nanoAOD_2017/ZZTo2L2Q/*"};}
	else if(process == "VV_ZZTo4L"){input_files = {"/data/disk0/nanoAOD_2017/ZZTo4L/*"};}
	else if(process == "VV_WZTo1L1Nu2Q"){input_files = {"/data/disk0/nanoAOD_2017/WZTo1L1Nu2Q/*"};}
	else if(process == "VV_WZTo2L2Q"){input_files = {"/data/disk0/nanoAOD_2017/WZTo2L2Q/*"};}
	else if(process == "VV_WZTo3LNu"){input_files = {"/data/disk0/nanoAOD_2017/WZTo3LNu/*"};}
	else if(process == "VV_WWTo1L1Nu2Q"){input_files = {"/data/disk0/nanoAOD_2017/WWTo1L1Nu2Q/*"};}
	else if(process == "VV_WWTo2L2Nu"){input_files = {"/data/disk0/nanoAOD_2017/WWTo2L2Nu/*"};}
	else if(process == "VV_WWToLNuQQ"){input_files = {"/data/disk0/nanoAOD_2017/WWToLNuQQ/*"};}
	else if(process == "VV_WWTolnuqq"){input_files = {"/data/disk0/nanoAOD_2017/WWTolnuqq/*"};}
	else if(process == "VV_WGToLNuG"){input_files = {"/data/disk0/nanoAOD_2017/WGToLNuG/*"};}
	else if(process == "VV_ZGToLLG"){input_files = {"/data/disk0/nanoAOD_2017/ZGToLLG/*"};}
	else if(process == "VVV_WWWTo4F"){input_files = {"/data/disk0/nanoAOD_2017/WWWTo4F/*"};}
	else if(process == "VVV_WWZTo4F"){input_files = {"/data/disk0/nanoAOD_2017/WWZTo4F/*"};}
	else if(process == "VVV_WZZ"){input_files = {"/data/disk0/nanoAOD_2017/WZZ/*"};}
	else if(process == "VVV_ZZZ"){input_files = {"/data/disk0/nanoAOD_2017/ZZZ/*"};}
	else if(process == "WPlusJets_WJetsToLNu"){input_files = {"/data/disk0/nanoAOD_2017/WJetsToLNu/*"};}
	else if(process == "ttbarV_ttWJetsToLNu"){input_files = {"/data/disk0/nanoAOD_2017/ttWJetsToLNu/*"};}
	else if(process == "ttbarV_ttWJetsToQQ"){input_files = {"/data/disk0/nanoAOD_2017/ttWJetsToQQ/*"};}
	else if(process == "ttbarV_ttgamma"){input_files = {"/data/disk0/nanoAOD_2017/ttgamma/*"};}
	else if(process == "ttbarV_ttZToLL"){input_files = {"/data/disk0/nanoAOD_2017/ttZToLL/*"};}
	else if(process == "ttbarV_ttHTobb"){input_files = {"/data/disk0/nanoAOD_2017/ttHTobb/*"};}
	else if(process == "ttbarV_ttHToNonbb"){input_files = {"/data/disk0/nanoAOD_2017/ttHToNonbb/*"};}
	else if(process == "ttbarV_ttZToLLNuNu"){input_files = {"/data/disk0/nanoAOD_2017/ttZToLLNuNu/*"};}
	else if(process == "ttbarV_ttZToQQ"){input_files = {"/data/disk0/nanoAOD_2017/ttZToQQ/*"};}
	else if(process == "ttbarV_ttZToQQ_ext"){input_files = {"/data/disk0/nanoAOD_2017/ttZToQQ_ext/*"};}
	else if(process == "data_DoubleEGRunB"){input_files = {"/data/disk0/nanoAOD_2017/DoubleEGRun2017B/*"};}
	else if(process == "data_DoubleEGRunC"){input_files = {"/data/disk0/nanoAOD_2017/DoubleEGRun2017C/*"};}
	else if(process == "data_DoubleEGRunD"){input_files = {"/data/disk0/nanoAOD_2017/DoubleEGRun2017D/*"};}
	else if(process == "data_DoubleEGRunE"){input_files = {"/data/disk0/nanoAOD_2017/DoubleEGRun2017E/*"};}
	else if(process == "data_DoubleEGRunF"){input_files = {"/data/disk0/nanoAOD_2017/DoubleEGRun2017F/*"};}
	else if(process == "data_DoubleMuonRunB"){input_files = {"/data/disk0/nanoAOD_2017/DoubleMuonRun2017B/*"};}
	else if(process == "data_DoubleMuonRunC"){input_files = {"/data/disk0/nanoAOD_2017/DoubleMuonRun2017C/*"};}
	else if(process == "data_DoubleMuonRunD"){input_files = {"/data/disk0/nanoAOD_2017/DoubleMuonRun2017D/*"};}
	else if(process == "data_DoubleMuonRunE"){input_files = {"/data/disk0/nanoAOD_2017/DoubleMuonRun2017E/*"};}
	else if(process == "data_DoubleMuonRunF"){input_files = {"/data/disk0/nanoAOD_2017/DoubleMuonRun2017F/*"};}
	else if(process == "data_SingleMuonRunB"){input_files = {"/data/disk3/nanoAOD_2017/SingleMuon_NanoAOD25Oct2019_RunB/*"};}
	else if(process == "data_SingleMuonRunC"){input_files = {"/data/disk3/nanoAOD_2017/SingleMuon_NanoAOD25Oct2019_RunC/*"};}
	else if(process == "data_SingleMuonRunD"){input_files = {"/data/disk3/nanoAOD_2017/SingleMuon_NanoAOD25Oct2019_RunD/*"};}
	else if(process == "data_SingleMuonRunE"){input_files = {"/data/disk3/nanoAOD_2017/SingleMuon_NanoAOD25Oct2019_RunE/*"};}
	else if(process == "data_SingleMuonRunF"){input_files = {"/data/disk3/nanoAOD_2017/SingleMuon_NanoAOD25Oct2019_RunF/*"};}
	else if(process == "data_SingleElectronRunB"){input_files = {"/data/disk3/nanoAOD_2017/SingleElectron_NanoAOD25Oct2019_RunB/*"};}
        else if(process == "data_SingleElectronRunC"){input_files = {"/data/disk3/nanoAOD_2017/SingleElectron_NanoAOD25Oct2019_RunC/*"};}
        else if(process == "data_SingleElectronRunD"){input_files = {"/data/disk3/nanoAOD_2017/SingleElectron_NanoAOD25Oct2019_RunD/*"};}
        else if(process == "data_SingleElectronRunE"){input_files = {"/data/disk3/nanoAOD_2017/SingleElectron_NanoAOD25Oct2019_RunE/*"};}
        else if(process == "data_SingleElectronRunF"){input_files = {"/data/disk3/nanoAOD_2017/SingleElectron_NanoAOD25Oct2019_RunF/*"};}
	else if(process == "Data_triggerSF"){input_files = {"/data/disk0/nanoAOD_2017/METRun2017B/*.root", "/data/disk0/nanoAOD_2017/METRun2017C/*.root", "/data/disk0/nanoAOD_2017/METRun2017D/*.root", "/data/disk0/nanoAOD_2017/METRun2017E/*.root", "/data/disk0/nanoAOD_2017/METRun2017F/*"};}
	else if(process == "MC_triggerSF_ttbar"){input_files = {"/data/disk0/nanoAOD_2017/ttbar_2l2nu/*.root", "/data/disk0/nanoAOD_2017/ttbar_madgraph_NanoAODv5/*", "/data/disk0/nanoAOD_2017/TTToHadronic/*", "/data/disk0/nanoAOD_2017/TTToSemileptonic/*", "/data/disk0/nanoAOD_2017/ttbar_aMCatNLO/*"};}
	else if(process == "MC_triggerSF_ZPlusJets"){input_files = {"/data/disk0/nanoAOD_2017/DYJetsToLL_NanoAODv5/*", "/data/disk0/nanoAOD_2017/DYJetsToLL_ext_NanoAODv5/*", "/data/disk0/nanoAOD_2017/DYJetsToLL_M10to50/*"};}
	else if(process == "NPL_File_ee_Blinded"){input_files = {"NPL_ee_output_2017_Blinded.root"};}
        else if(process == "NPL_File_mumu_Blinded"){input_files = {"NPL_mumu_output_2017_Blinded.root"};}
        else if(process == "NPL_File_ee_Unblinded"){input_files = {"NPL_ee_output_2017.root"};}
        else if(process == "NPL_File_mumu_Unblinded"){input_files = {"NPL_mumu_output_2017.root"};}
	else{std::cout << "You inputted the process: " << process << " for the year " << year << ". Please input an MC signal, background of dataset name." << std::endl;}

}
else if(year == "2018"){

	if(process == "tZq"){input_files = {"/data/disk1/nanoAOD_2018/tZq_ll/*"};}
	else if(process == "ZPlusJets_M50_aMCatNLO"){input_files = {"/data/disk1/nanoAOD_2018/DYJetsToLL_NanoAODv5/*"};}
	else if(process == "ZPlusJets_M50_aMCatNLO_ext"){input_files = {"/data/disk1/nanoAOD_2018/DYJetsToLL_ext_NanoAODv5/*"};}
	else if(process == "ZPlusJets_M10To50_Madgraph"){input_files = {"/data/disk1/nanoAOD_2018/DYJetsToLL_M10to50/*"};}
	else if(process == "SingleTop_tchannel_top"){input_files = {"/data/disk1/nanoAOD_2018/ST_tchannel_top/*"};}
	else if(process == "SingleTop_tchannel_tbar"){input_files = {"/data/disk1/nanoAOD_2018/ST_tchannel_tbar/*"};}
	else if(process == "SingleTop_schannel"){input_files = {"/data/disk1/nanoAOD_2018/ST_schannel/*"};}
	else if(process == "SingleTop_tW"){input_files = {"/data/disk1/nanoAOD_2018/ST_tW/*"};}
	else if(process == "SingleTop_tbarW"){input_files = {"/data/disk1/nanoAOD_2018/ST_tbarW/*"};}
	else if(process == "SingleTop_tHq"){input_files = {"/data/disk1/nanoAOD_2018/tHq/*"};}
	else if(process == "SingleTop_tZq_W_lept_Z_had"){input_files = {"/data/disk1/nanoAOD_2018/tZq_W_lept_Z_had/*"};} 
	else if(process == "SingleTop_tWZ_tWll"){input_files = {"/data/disk1/nanoAOD_2018/tWZ_tWLL_NanoAODv5/*"};} 
	else if(process == "ttbar_2l2nu"){input_files = {"/data/disk1/nanoAOD_2018/ttbar_2l2nu/*"};}
	else if(process == "ttbar_madgraph_NanoAODv5"){input_files = {"/data/disk1/nanoAOD_2018/ttbar_madgraph/*"};}
	else if(process == "ttbar_TTToHadronic"){input_files = {"/data/disk1/nanoAOD_2018/TTToHadronic/*"};}
	else if(process == "ttbar_TTToSemileptonic"){input_files = {"/data/disk1/nanoAOD_2018/TTToSemileptonic/*"};}
	else if(process == "ttbar_aMCatNLO"){input_files = {"/data/disk1/nanoAOD_2018/ttbar_aMCatNLO/*"};}
	else if(process == "VV_ZZTo2Q2Nu"){input_files = {"/data/disk1/nanoAOD_2018/ZZTo2Q2Nu/*"};}
	else if(process == "VV_ZZTo2L2Nu"){input_files = {"/data/disk1/nanoAOD_2018/ZZTo2L2Nu/*"};}
	else if(process == "VV_ZZTo2L2Q"){input_files = {"/data/disk1/nanoAOD_2018/ZZTo2L2Q/*"};}
	else if(process == "VV_ZZTo4L"){input_files = {"/data/disk1/nanoAOD_2018/ZZTo4L/*"};}
	else if(process == "VV_WZTo1L1Nu2Q"){input_files = {"/data/disk1/nanoAOD_2018/WZTo1L1Nu2Q/*"};}
	else if(process == "VV_WZTo2L2Q"){input_files = {"/data/disk1/nanoAOD_2018/WZTo2L2Q/*"};}
	else if(process == "VV_WZTo3LNu"){input_files = {"/data/disk1/nanoAOD_2018/WZTo3LNu/*"};}
	else if(process == "VV_WWTo1L1Nu2Q"){input_files = {"/data/disk1/nanoAOD_2018/WWTo1L1Nu2Q/*"};}
	else if(process == "VV_WWTo2L2Nu"){input_files = {"/data/disk1/nanoAOD_2018/WWTo2L2Nu/*"};}
	else if(process == "VV_WWToLNuQQ"){input_files = {"/data/disk1/nanoAOD_2018/WWToLNuQQ/*"};}
	else if(process == "VV_WGToLNuG"){input_files = {"/data/disk1/nanoAOD_2018/WGToLNuG/*"};}
	else if(process == "VV_ZGToLLG"){input_files = {"/data/disk1/nanoAOD_2018/ZGToLLG/*"};}
	else if(process == "VVV_WWWTo4F"){input_files = {"/data/disk1/nanoAOD_2018/WWWTo4F/*"};}
	else if(process == "VVV_WWZTo4F"){input_files = {"/data/disk1/nanoAOD_2018/WWZTo4F/*"};}
	else if(process == "VVV_WZZ"){input_files = {"/data/disk1/nanoAOD_2018/WZZ/*"};}
	else if(process == "VVV_ZZZ"){input_files = {"/data/disk1/nanoAOD_2018/ZZZ/*"};}
	else if(process == "WPlusJets_WJetsToLNu"){input_files = {"/data/disk1/nanoAOD_2018/WJetsToLNu/*"};}
	else if(process == "ttbarV_ttWJetsToLNu"){input_files = {"/data/disk1/nanoAOD_2018/ttWJetsToLNu/*"};}
	else if(process == "ttbarV_ttWJetsToQQ"){input_files = {"/data/disk1/nanoAOD_2018/ttWJetsToQQ/*"};}
	else if(process == "ttbarV_ttgamma"){input_files = {"/data/disk1/nanoAOD_2018/ttgamma/*"};}
	else if(process == "ttbarV_ttZToLL"){input_files = {"/data/disk1/nanoAOD_2018/ttZToLL/*"};}
	else if(process == "ttbarV_ttHTobb"){input_files = {"/data/disk1/nanoAOD_2018/ttHTobb/*"};}
	else if(process == "ttbarV_ttHToNonbb"){input_files = {"/data/disk1/nanoAOD_2018/ttHToNonbb/*"};}
	else if(process == "ttbarV_ttZToLLNuNu"){input_files = {"/data/disk1/nanoAOD_2018/ttZToLLNuNu/*"};}
	else if(process == "ttbarV_ttZToQQ"){input_files = {"/data/disk1/nanoAOD_2018/ttZToQQ/*"};}
	else if(process == "ttbarV_ttZToQQ_ext"){input_files = {"/data/disk1/nanoAOD_2018/ttZToQQ_ext/*"};}
	else if(process == "data_EGRunB"){input_files = {"/data/disk3/nanoAOD_2018/EGammaRunB/*"};}
	else if(process == "data_EGRunC"){input_files = {"/data/disk3/nanoAOD_2018/EGammaRunC/*"};}
	else if(process == "data_EGRunD"){input_files = {"/data/disk3/nanoAOD_2018/EGammaRunD/*"};}
	else if(process == "data_DoubleMuonRunB"){input_files = {"/data/disk1/nanoAOD_2018/DoubleMuonRun2018B/*"};}
	else if(process == "data_DoubleMuonRunC"){input_files = {"/data/disk1/nanoAOD_2018/DoubleMuonRun2018C/*"};}
	else if(process == "data_DoubleMuonRunD"){input_files = {"/data/disk1/nanoAOD_2018/DoubleMuonRun2018D/*"};}
	else if(process == "data_SingleMuonRunB"){input_files = {"/data/disk3/nanoAOD_2018/SingleMuon_NanoAOD25Oct_2019_RunB/*"};}
        else if(process == "data_SingleMuonRunC"){input_files = {"/data/disk3/nanoAOD_2018/SingleMuon_NanoAOD25Oct_2019_RunC/*"};}
        else if(process == "data_SingleMuonRunD"){input_files = {"/data/disk3/nanoAOD_2018/SingleMuon_NanoAOD25Oct_2019_RunD/*"};}
	else if(process == "Data_triggerSF"){input_files = {"/data/disk1/nanoAOD_2018/METRun2018B/*.root", "/data/disk1/nanoAOD_2018/METRun2018C/*.root", "/data/disk1/nanoAOD_2018/METRun2018D/*.root"};}
	else if(process == "MC_triggerSF_ttbar"){input_files = {"/data/disk0/nanoAOD_2017/ttbar_2l2nu/*.root", "/data/disk0/nanoAOD_2017/ttbar_madgraph_NanoAODv5/*", "/data/disk0/nanoAOD_2017/TTToHadronic/*", "/data/disk0/nanoAOD_2017/TTToSemileptonic/*", "/data/disk0/nanoAOD_2017/ttbar_aMCatNLO/*"};}
	else if(process == "MC_triggerSF_ZPlusJets"){input_files = {"/data/disk0/nanoAOD_2018/DYJetsToLL_NanoAODv5/*", "/data/disk0/nanoAOD_2018/DYJetsToLL_ext_NanoAODv5/*", "/data/disk0/nanoAOD_2018/DYJetsToLL_M10to50/*"};}
	else if(process == "NPL_File_ee_Blinded"){input_files = {"NPL_ee_output_2018_Blinded.root"};}
        else if(process == "NPL_File_mumu_Blinded"){input_files = {"NPL_mumu_output_2018_Blinded.root"};}
        else if(process == "NPL_File_ee_Unblinded"){input_files = {"NPL_ee_output_2018.root"};}
        else if(process == "NPL_File_mumu_Unblinded"){input_files = {"NPL_mumu_output_2018.root"};}
	else{std::cout << "You inputted the process: " << process << " for the year " << year << ". Please input an MC signal, background of dataset name." << std::endl;}


}
else{std::cout << "Script only for 2016, 2017 or 2018 samples" << std::endl;}





RDataFrame d("Events", input_files);
auto d_dataframe = d.Range(0, 1000);

//RDataFrame d_dataframe("Events", input_files);





//Functions for events that pass the ee, mumu, or emu triggers
//To prevent double counting single and double lepton datasets
auto SingleElectron{[&year](

const bool& HLT_Ele32_WPTight_Gsf_L1DoubleEG,
const bool& HLT_Ele32_eta2p1_WPTight_Gsf,
const bool& HLT_Ele35_WPTight_Gsf,
const bool& HLT_Ele25_eta2p1_WPTight_Gsf, 
const bool& HLT_Ele27_WPTight_Gsf)-> bool{


  //for 2016 see: https://twiki.cern.ch/twiki/bin/view/CMS/TopTriggerYear2016
  //for 2017 see: https://twiki.cern.ch/twiki/bin/view/CMS/TopTriggerYear2017
  //for 2018 see: https://twiki.cern.ch/twiki/bin/view/CMS/TopTriggerYear2018

  if(year == "2016"){

	return  HLT_Ele25_eta2p1_WPTight_Gsf > 0 ||
		HLT_Ele27_WPTight_Gsf > 0 ||
		HLT_Ele32_eta2p1_WPTight_Gsf > 0;

  }
  else if(year == "2017"){

  	return //HLT_Ele32_WPTight_Gsf_L1DoubleEG > 0 || 
   	       HLT_Ele35_WPTight_Gsf > 0; 

  }
  else if(year == "2018"){

	return HLT_Ele32_WPTight_Gsf_L1DoubleEG > 0;

  }
  else{std::cout << "Choose a year out of 2016, 2017 or 2018 for the trigger paths" << std::endl;}


}};


auto DoubleElectron{[&year](

const bool& HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL,
const bool& HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ

)->bool{


  //for 2016 see: https://twiki.cern.ch/twiki/bin/view/CMS/TopTriggerYear2016
  //for 2017 see: https://twiki.cern.ch/twiki/bin/view/CMS/TopTriggerYear2017
  //for 2018 see: https://twiki.cern.ch/twiki/bin/view/CMS/TopTriggerYear2018

  if(year == "2016"){

	return HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ > 0;

  }
  else if(year == "2017"){
  
	return HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL > 0 ||
   	       HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ > 0; 

  }
  else if(year == "2018"){
	
	return HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL > 0 || 
               HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ > 0;

  }
  else{std::cout << "Choose a year out of 2016, 2017 or 2018 for the trigger paths" << std::endl;}


}};



auto SingleMuon{[&year](

const bool& HLT_IsoMu24,
const bool& HLT_IsoMu27,
const bool& HLT_IsoMu24_eta2p1

)->bool{

  
  //for 2016 see: https://twiki.cern.ch/twiki/bin/view/CMS/TopTriggerYear2016
  //for 2017 see: https://twiki.cern.ch/twiki/bin/view/CMS/TopTriggerYear2017
  //for 2018 see: https://twiki.cern.ch/twiki/bin/view/CMS/TopTriggerYear2018

  if(year == "2016"){

  	return HLT_IsoMu24  > 0 ||
	       HLT_IsoMu24_eta2p1  > 0;

  }
  else if(year == "2017"){

  	return HLT_IsoMu24 > 0 || 
   	       HLT_IsoMu27 > 0; 

  }
  else if(year == "2018"){

	return HLT_IsoMu24 > 0;

  }
  else{std::cout << "Please choose a year out of 2016, 2017 or 2018 for the trigger paths" << std::endl;}

}};



auto DoubleMuon{[&year](

const bool& HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, 
const bool& HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8, 
const bool& HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8)->bool{


  //for 2016 see: https://twiki.cern.ch/twiki/bin/view/CMS/TopTriggerYear2016
  //for 2017 see: https://twiki.cern.ch/twiki/bin/view/CMS/TopTriggerYear2017
  //for 2018 see: https://twiki.cern.ch/twiki/bin/view/CMS/TopTriggerYear2018


  if(year == "2016"){return HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ > 0;}
  else if(year == "2017"){

  	return HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ > 0 || 
   	       HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 > 0; 

  }
  else if(year == "2018"){

	return HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 > 0 ||
	       HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 > 0;

  }
  else{std::cout << "Choose a year out of 2016, 2017 or 2018 for the trigger paths" << std::endl;}


}};




auto MuonElectron{[&year](

const bool& HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, 
const bool& HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, 
const bool& HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, 
const bool& HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,
const bool& HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,
const bool& HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL)->bool{



 //for 2016 see: https://twiki.cern.ch/twiki/bin/view/CMS/TopTriggerYear2016
 //for 2017 see: https://twiki.cern.ch/twiki/bin/view/CMS/TopTriggerYear2017
 //for 2018 see: https://twiki.cern.ch/twiki/bin/view/CMS/TopTriggerYear2018


 if(year == "2016"){

	return //HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ > 0 ||  (branch not present in MET)
               //HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ > 0 ||  (branch not present in MET)
               //HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ > 0 ||   (branch not present in MET)
	       HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL > 0 ||
               //HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL > 0 || (branch not present in MET)
               HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL > 0;


 }
 else if(year == "2017"){

 	return HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ > 0 || 
   	       HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ > 0 || 
  	       HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ > 0; 

 }
 else if(year == "2018"){

	return HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ > 0 ||
               HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ > 0 ||
               HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ > 0;


 }
 else{std::cout << "Choose a year out of 2016, 2017 or 2018 for the trigger paths" << std::endl;}


}};











if(year == "2016"){
        if(ttbarCR == false){MinElectronPt = 15; MaxElectronPt = 35; MinMuonPt = 20; MaxMuonPt = 26; MaxTrackerEta = 2.4;}
        else{MinElectronPt = 25; MinMuonPt = 25;}
}
else if(year == "2017" || year == "2018"){
	if(ttbarCR == false){MinElectronPt = 15; MaxElectronPt = 38; MinMuonPt = 20; MaxMuonPt = 29; MaxTrackerEta = 2.5;}
        else{MinElectronPt = 25; MinMuonPt = 25;}
}
else{std::cout << "Choose the year out of 2016, 2017 or 2018, and choose ttbarCR as either true or false";}



std::string JetMassInput, JetPtInput, JetEtaInput, JetPhiInput;

if(

process != "data_DoubleEGRunB" &&
process != "data_DoubleEGRunC" &&
process != "data_DoubleEGRunD" &&
process != "data_DoubleEGRunE" &&
process != "data_DoubleEGRunF" &&
process != "data_DoubleEGRunG" &&
process != "data_DoubleEGRunH" &&
process != "data_EGRunB" &&
process != "data_EGRunC" &&
process != "data_EGRunD" &&
process != "data_SingleElectronRunB" &&
process != "data_SingleElectronRunC" &&
process != "data_SingleElectronRunD" &&
process != "data_SingleElectronRunE" &&
process != "data_SingleElectronRunF" &&
process != "data_SingleElectronRunG" &&
process != "data_SingleElectronRunH" &&
process != "data_DoubleMuonRunB" &&
process != "data_DoubleMuonRunC" &&
process != "data_DoubleMuonRunD" &&
process != "data_DoubleMuonRunE" &&
process != "data_DoubleMuonRunF" &&
process != "data_DoubleMuonRunG" &&
process != "data_DoubleMuonRunH" &&
process != "data_SingleMuonRunB" &&
process != "data_SingleMuonRunC" &&
process != "data_SingleMuonRunD" &&
process != "data_SingleMuonRunE" &&
process != "data_SingleMuonRunF" &&
process != "data_SingleMuonRunG" &&
process != "data_SingleMuonRunH" 

){


  JetMassInput = "SmearedJetMass";
  JetPtInput = "SmearedJetPt";
  JetEtaInput = "SmearedJetEta";
  JetPhiInput = "SmearedJetPhi";


}
else if(

process == "data_DoubleEGRunB" ||
process == "data_DoubleEGRunC" ||
process == "data_DoubleEGRunD" ||
process == "data_DoubleEGRunE" ||
process == "data_DoubleEGRunF" ||
process == "data_DoubleEGRunG" ||
process == "data_DoubleEGRunH" ||
process == "data_EGRunB" ||
process == "data_EGRunC" ||
process == "data_EGRunD" ||
process == "data_SingleElectronRunB" ||
process == "data_SingleElectronRunC" ||
process == "data_SingleElectronRunD" ||
process == "data_SingleElectronRunE" ||
process == "data_SingleElectronRunF" ||
process == "data_SingleElectronRunG" ||
process == "data_SingleElectronRunH" ||
process == "data_DoubleMuonRunB" ||
process == "data_DoubleMuonRunC" ||
process == "data_DoubleMuonRunD" ||
process == "data_DoubleMuonRunE" ||
process == "data_DoubleMuonRunF" ||
process == "data_DoubleMuonRunG" ||
process == "data_DoubleMuonRunH" ||
process == "data_SingleMuonRunB" ||
process == "data_SingleMuonRunC" ||
process == "data_SingleMuonRunD" ||
process == "data_SingleMuonRunE" ||
process == "data_SingleMuonRunF" ||
process == "data_SingleMuonRunG" ||
process == "data_SingleMuonRunH" 

){


  JetMassInput = "Jet_mass";
  JetPtInput = "Jet_pt";
  JetEtaInput = "Jet_eta";
  JetPhiInput = "Jet_phi";

}
else{std::cout << "Please choose MC or data as the input" << std::endl;}



//Events that only pass the lepton selection criteria


auto LeadingElectron_dz_ECALBarrel_function{[](

const float& LeadingElectron_pT,
const floats& Electron_eta_Selection,
const floats& Electron_dz){


  floats OutVec{};

  for(int i = 0; i < Electron_eta_Selection.size(); i++){

	if(LeadingElectron_pT > MaxElectronPt && abs(Electron_eta_Selection.at(i)) < 1.479 && Electron_dz.at(i) < 0.1){OutVec.push_back(Electron_dz.at(i));}
	else{float zero = 0.0; OutVec.push_back(zero);}

  }

  return OutVec;  


}};


auto SubleadingElectron_dz_ECALBarrel_function{[](

const float& SubleadingElectron_pT,
const floats& Electron_eta_Selection,
const floats& Electron_dz){


  floats OutVec{};

  for(int i = 0; i < Electron_eta_Selection.size(); i++){

        if(SubleadingElectron_pT > MinElectronPt && abs(Electron_eta_Selection.at(i)) < 1.479 && Electron_dz.at(i) < 0.1){OutVec.push_back(Electron_dz.at(i));}
	else{float zero = 0.0; OutVec.push_back(zero);}

  }

  return OutVec;


}};



auto LeadingElectron_dz_ECALEndcaps_function{[](

const float& LeadingElectron_pT,
const floats& Electron_eta_Selection,
const floats& Electron_dz){


  floats OutVec{};

  for(int i = 0; i < Electron_eta_Selection.size(); i++){

        if(LeadingElectron_pT > MaxElectronPt && abs(Electron_eta_Selection.at(i)) > 1.479 && abs(Electron_eta_Selection.at(i)) < 3.0 && Electron_dz.at(i) < 0.2){OutVec.push_back(Electron_dz.at(i));}
	else{float zero = 0.0; OutVec.push_back(zero);}

  }

  return OutVec;

}};



auto SubleadingElectron_dz_ECALEndcaps_function{[](

const float& SubleadingElectron_pT,
const floats& Electron_eta_Selection,
const floats& Electron_dz){


  floats OutVec{};

  for(int i = 0; i < Electron_eta_Selection.size(); i++){

        if(SubleadingElectron_pT > MinElectronPt && abs(Electron_eta_Selection.at(i)) > 1.479 && abs(Electron_eta_Selection.at(i)) < 3.0 && Electron_dz.at(i) <  0.2){OutVec.push_back(Electron_dz.at(i));}
	else{float zero = 0.0; OutVec.push_back(zero);}

  }

  return OutVec;


}};



auto LeadingElectron_dxy_ECALBarrel_function{[](

const float& LeadingElectron_pT,
const floats& Electron_eta_Selection,
const floats& Electron_dxy){


  floats OutVec{};

  for(int i = 0; i < Electron_eta_Selection.size(); i++){

        if(LeadingElectron_pT > MaxElectronPt && abs(Electron_eta_Selection.at(i)) < 1.479  && Electron_dxy.at(i) < 0.05){OutVec.push_back(Electron_dxy.at(i));}
	else{float zero = 0.0; OutVec.push_back(zero);}

  }

  return OutVec;

}};



auto SubleadingElectron_dxy_ECALBarrel_function{[](

const float& SubleadingElectron_pT,
const floats& Electron_eta_Selection,
const floats& Electron_dxy){


  floats OutVec{};

  for(int i = 0; i < Electron_eta_Selection.size(); i++){

        if(SubleadingElectron_pT > MinElectronPt && abs(Electron_eta_Selection.at(i)) < 1.479 && Electron_dxy.at(i) < 0.05){OutVec.push_back(Electron_dxy.at(i));}
	else{float zero = 0.0; OutVec.push_back(zero);}

  }

  return OutVec;

}};

auto LeadingElectron_dxy_ECALEndcaps_function{[](

const float& LeadingElectron_pT,
const floats& Electron_eta_Selection,
const floats& Electron_dxy){


  floats OutVec{};

  for(int i = 0; i < Electron_eta_Selection.size(); i++){

        if(LeadingElectron_pT > MaxElectronPt && abs(Electron_eta_Selection.at(i)) > 1.479 && abs(Electron_eta_Selection.at(i)) < 3.0  && Electron_dxy.at(i) < 0.1){OutVec.push_back(Electron_dxy.at(i));}
	else{float zero = 0.0; OutVec.push_back(zero);}

  }

  return OutVec;


}};

auto SubleadingElectron_dxy_ECALEndcaps_function{[](

const float& SubleadingElectron_pT,
const floats& Electron_eta_Selection,
const floats& Electron_dxy){


  floats OutVec{};

  for(int i = 0; i < Electron_eta_Selection.size(); i++){

        if(SubleadingElectron_pT > MinElectronPt && abs(Electron_eta_Selection.at(i)) > 1.479 && abs(Electron_eta_Selection.at(i)) < 3.0 && Electron_dxy.at(i) < 0.1){OutVec.push_back(Electron_dxy.at(i));}
	else{float zero = 0.0; OutVec.push_back(zero);}

  }

  return OutVec;


}};







auto OppositeSign{[](const ints& charges){

  return charges.size() == 2 ? signbit(charges.at(0)) != signbit(charges.at(1)) : false;

}};

auto OppositeSign_emu{[](const ints& charges1, const ints& charges2){

  return (charges1.size() == 1 && charges2.size() == 1) ? signbit(charges1.at(0)) != signbit(charges2.at(0)) : false;

}};


auto SameSign{[](const ints& charges){

  return charges.size() == 2 ? signbit(charges.at(0)) == signbit(charges.at(1)) : false;

}};


auto OppositeSignNonPrompt{[](const ints& charges, const chars& Lepton_genPartFlav){

  bool OppositeSignChargeCheck = charges.size() == 2 ? signbit(charges.at(0)) != signbit(charges.at(1)) : false;
  bool LeptonNonPromptCheck = all_of(Lepton_genPartFlav.begin(), Lepton_genPartFlav.end(), [](int i){return i != 1;});

  return OppositeSignChargeCheck && (LeptonNonPromptCheck == 1);

}};


auto SameSignNonPrompt{[](const ints& charges, const chars& Lepton_genPartFlav){
  

  bool SameSignChargeCheck = charges.size() == 2 ? signbit(charges.at(0)) == signbit(charges.at(1)) : false;  
  bool LeptonNonPromptCheck = all_of(Lepton_genPartFlav.begin(), Lepton_genPartFlav.end(), [](int i){return i != 1;});

  return SameSignChargeCheck && (LeptonNonPromptCheck == 1);

}};



auto OppositeSignPrompt{[](const ints& charges, const chars& Lepton_genPartFlav){

  
  bool OppositeSignChargeCheck = charges.size() == 2 ? signbit(charges.at(0)) != signbit(charges.at(1)) : false;
  bool LeptonPromptCheck = all_of(Lepton_genPartFlav.begin(), Lepton_genPartFlav.end(), [](int i){return i == 1;});

  return OppositeSignChargeCheck && (LeptonPromptCheck == 1);

}};




auto SameSignPrompt{[](const ints& charges, const chars& Lepton_genPartFlav){

  
  bool SameSignChargeCheck = charges.size() == 2 ? signbit(charges.at(0)) == signbit(charges.at(1)) : false;
  bool LeptonPromptCheck = all_of(Lepton_genPartFlav.begin(), Lepton_genPartFlav.end(), [](int i){return i == 1;});  

  return SameSignChargeCheck && (LeptonPromptCheck == 1);

}};




auto ElectronsFunction{[](

const int targetID,
const floats& Electron_pt,
const floats& Electron_eta,
const ints& Electron_cutBased,
const bools& Electron_isPFcand

){
 
  return (Electron_pt > MinElectronPt && (abs(Electron_eta) < MaxTrackerEta && (abs(Electron_eta) < 1.442 || abs(Electron_eta) > 1.566) ) && Electron_cutBased >= targetID && Electron_isPFcand);

}};


auto ElectronsFunctionEmu{[](

const int targetID,
const floats& Electron_pt,
const floats& Electron_eta,
const ints& Electron_cutBased,
const bools& Electron_isPFcand

){

  return (Electron_pt > MinElectronPtEmu && (abs(Electron_eta) < MaxTrackerEta && (abs(Electron_eta) < 1.442 || abs(Electron_eta) > 1.566) ) && Electron_cutBased >= targetID && Electron_isPFcand);

}};



auto TightElectronsFunction{[&ElectronsFunction](

const floats& Electron_pt,
const floats& Electron_eta,
const ints& Electron_cutBased,
const bools& Electron_isPFcand


){

  return ElectronsFunction(4, Electron_pt, Electron_eta, Electron_cutBased, Electron_isPFcand);

}};


auto TightElectronsFunctionEmu{[&ElectronsFunctionEmu](

const floats& Electron_pt,
const floats& Electron_eta,
const ints& Electron_cutBased,
const bools& Electron_isPFcand


){

  return ElectronsFunctionEmu(4, Electron_pt, Electron_eta, Electron_cutBased, Electron_isPFcand);

}};


auto LooseElectronsFunction{[&ElectronsFunction](

const floats& Electron_pt,
const floats& Electron_eta,
const ints& Electron_cutBased,
const bools& Electron_isPFcand


){

  return ElectronsFunction(1, Electron_pt, Electron_eta, Electron_cutBased, Electron_isPFcand);

}};


auto LooseElectronsFunctionEmu{[&ElectronsFunctionEmu](

const floats& Electron_pt,
const floats& Electron_eta,
const ints& Electron_cutBased,
const bools& Electron_isPFcand


){

  return ElectronsFunctionEmu(1, Electron_pt, Electron_eta, Electron_cutBased, Electron_isPFcand);

}};


auto MuonsFunction{[](

const float target_iso, 
const bools& isPFs, 
const floats& Muon_pt, 
const floats& Muon_eta, 
const bools& ids, 
const floats& isos

){

  return (isPFs && Muon_pt > MinMuonPt && abs(Muon_eta) < MaxTrackerEta && ids && isos <= target_iso);

}};


auto MuonsFunctionEmu{[](

const float target_iso,
const bools& isPFs,
const floats& Muon_pt,
const floats& Muon_eta,
const bools& ids,
const floats& isos

){

  return (isPFs && Muon_pt > MinMuonPtEmu && abs(Muon_eta) < MaxTrackerEta && ids && isos <= target_iso);

}};


auto TightMuonsFunction{[&MuonsFunction](const bools& isPFs, const floats& pts, const floats& etas, const bools& ids, const floats& isos) {

  return MuonsFunction(0.25, isPFs, pts, etas, ids, isos);

}};


auto TightMuonsFunctionEmu{[&MuonsFunctionEmu](const bools& isPFs, const floats& pts, const floats& etas, const bools& ids, const floats& isos) {

  return MuonsFunctionEmu(0.25, isPFs, pts, etas, ids, isos);

}};


auto LooseMuonsFunction{[&MuonsFunction](const bools& isPFs, const floats& pts, const floats& etas, const bools& ids, const floats& isos) {

  return MuonsFunction(0.15, isPFs, pts, etas, ids, isos);

}};


auto LooseMuonsFunctionEmu{[&MuonsFunctionEmu](const bools& isPFs, const floats& pts, const floats& etas, const bools& ids, const floats& isos) {

  return MuonsFunctionEmu(0.15, isPFs, pts, etas, ids, isos);

}};


auto lep_cut_ee{[](

const floats& tight_ele_pts, 
const floats& loose_ele_pts, 
const bool os,
const unsigned int& nElectron,
const floats& LeadingElectron_dz_ECALBarrel,
const floats& LeadingElectron_dxy_ECALBarrel,
const floats& LeadingElectron_dz_ECALEndcaps,
const floats& LeadingElectron_dxy_ECALEndcaps,
const floats& SubleadingElectron_dz_ECALBarrel,
const floats& SubleadingElectron_dxy_ECALBarrel,
const floats& SubleadingElectron_dz_ECALEndcaps,
const floats& SubleadingElectron_dxy_ECALEndcaps
){


  const bool ele_cut{tight_ele_pts.size() == 2 && tight_ele_pts.size() == loose_ele_pts.size()};
  bool lead_pt_cut{false};

  lead_pt_cut = tight_ele_pts.empty() ? false : *max_element(tight_ele_pts.begin(), tight_ele_pts.end()) > MaxElectronPt;


  return

  os &&
  lead_pt_cut &&
  ele_cut &&
  nElectron == 2 &&
  LeadingElectron_dz_ECALBarrel.at(0) < 0.1 &&
  LeadingElectron_dxy_ECALBarrel.at(0) < 0.05 &&
  LeadingElectron_dz_ECALEndcaps.at(0) < 0.2 &&
  LeadingElectron_dxy_ECALEndcaps.at(0) < 0.1 &&
  SubleadingElectron_dz_ECALBarrel.at(0) < 0.1 &&
  SubleadingElectron_dxy_ECALBarrel.at(0) < 0.05 &&
  SubleadingElectron_dz_ECALEndcaps.at(0) < 0.2 &&
  SubleadingElectron_dxy_ECALEndcaps.at(0) < 0.1;

}};



auto lep_cut_mumu{[](const floats& tight_mu_pts, const floats& loose_mu_pts, const bool os, const unsigned int nMuon) {

        
  const bool mu_cut{tight_mu_pts.size() == 2 && tight_mu_pts.size() == loose_mu_pts.size()};
  bool lead_pt_cut{false};

  lead_pt_cut = tight_mu_pts.empty() ? false : *std::max_element(tight_mu_pts.begin(), tight_mu_pts.end()) > MaxMuonPt;

  return 

  os && 
  lead_pt_cut && 
  mu_cut &&
  nMuon == 2;

  }};



//emu only for ttbar control region and trigger SF calculations
auto lep_cut_emu{[](

const floats& tight_ele_pts,
const floats& loose_ele_pts,
const bool os,
const unsigned int& nElectron,
const floats& tight_mu_pts, 
const floats& loose_mu_pts, 
const unsigned int nMuon

){


  const bool emu_cut{tight_ele_pts.size() == 1 && tight_mu_pts.size() == 1 && (tight_mu_pts.size() == loose_mu_pts.size()) && (tight_mu_pts.size() == loose_mu_pts.size())};

  return 

  os &&
  emu_cut &&
  nElectron == 1 &&
  nMuon == 1;

}};





auto LeadingVariable{[](const floats& variable){


  if(variable.size() > 0){

  float first_largest_value = variable.at(0);

        for(int i = 1; i < variable.size(); i++){

                if(variable.at(i) > first_largest_value){
                        first_largest_value = variable.at(i);

                }

        }

  return first_largest_value;

  }
  else{float zero = 0.0; return zero;}

}};


auto LeadingVariableEmu{[](const floats& variable1, const floats& variable2){
  
  float first_largest_value_electron, first_largest_value_muon; 

  if(variable1.size() > 0){

  first_largest_value_electron = variable1.at(0);

        for(int i = 1; i < variable1.size(); i++){

                if(variable1.at(i) > first_largest_value_electron){
                        first_largest_value_electron = variable1.at(i);

                }

        }

  }

  if(variable2.size() > 0){
  
  first_largest_value_muon = variable2.at(0);
        
        for(int i = 1; i < variable2.size(); i++){
                
                if(variable2.at(i) > first_largest_value_muon){
                        first_largest_value_muon = variable2.at(i);
                
                }
        
        }

  }


  if(first_largest_value_electron > first_largest_value_muon){return first_largest_value_electron;}
  else{return first_largest_value_muon;}


}};



auto SubleadingVariable{[](const floats& variable){


  if(variable.size() == 0){float zero = 0.0; return zero;}
  else{

  	if(variable.size() > 1){

  	float first_largest_value = variable.at(0);

	  	for(int i = 1; i < variable.size(); i++){

			if(variable.at(i) > first_largest_value){
				first_largest_value = variable.at(i);

			}

	  	}

  	float second_largest_value = INT_MIN;

	  	for(int i = 0; i < variable.size(); i++){

			if( (variable.at(i) > second_largest_value) && (variable.at(i) < first_largest_value) ){
				second_largest_value = variable.at(i);
			}
	
	  	}


  	return second_largest_value;


  	}
  	else{return variable.at(0);}

  }


}};



auto SubleadingVariableEmu{[](const floats& variable1, const floats& variable2){

 
  float first_largest_value_electron, first_largest_value_muon;

  if(variable1.size() > 0){

  first_largest_value_electron = variable1.at(0);

        for(int i = 1; i < variable1.size(); i++){

                if(variable1.at(i) > first_largest_value_electron){
                        first_largest_value_electron = variable1.at(i);

                }

        }

  }

  if(variable2.size() > 0){

  first_largest_value_muon = variable2.at(0);

        for(int i = 1; i < variable2.size(); i++){

                if(variable2.at(i) > first_largest_value_muon){
                        first_largest_value_muon = variable2.at(i);

                }

        }

  }


  if(first_largest_value_electron < first_largest_value_muon){return first_largest_value_electron;}
  else{return first_largest_value_muon;}


}};


auto ThirdLeadingVariable{[](const floats& variable){


  if(variable.size() > 2){

  float first_largest_value = variable.at(0);

        for(int i = 1; i < variable.size(); i++){

                if(variable.at(i) > first_largest_value){
                        first_largest_value = variable.at(i);

                }

        }

  float second_largest_value = INT_MIN;

        for(int i = 0; i < variable.size(); i++){

                if( (variable.at(i) > second_largest_value) && (variable.at(i) < first_largest_value) ){
                        second_largest_value = variable.at(i);
                }
        
        }


  float third_largest_value = INT_MIN;

	for(int i = 0; i < variable.size(); i++){
	
		if( (variable.at(i) > third_largest_value) && (variable.at(i) < second_largest_value) ){
		
			third_largest_value = variable.at(i);			

		}
	
	}


  return third_largest_value;


  }

}};



auto FourthLeadingVariable{[](const floats& variable){


  if(variable.size() > 3){

  float first_largest_value = variable.at(0);

        for(int i = 1; i < variable.size(); i++){

                if(variable.at(i) > first_largest_value){
                        first_largest_value = variable.at(i);

                }

        }

  float second_largest_value = INT_MIN;

        for(int i = 0; i < variable.size(); i++){

                if( (variable.at(i) > second_largest_value) && (variable.at(i) < first_largest_value) ){
                        second_largest_value = variable.at(i);
                }

        }


  float third_largest_value = INT_MIN;

        for(int i = 0; i < variable.size(); i++){

                if( (variable.at(i) > third_largest_value) && (variable.at(i) < second_largest_value) ){

                        third_largest_value = variable.at(i);

                }

        }


  float fourth_largest_value = INT_MIN;

	for(int i = 0; i < variable.size(); i++){

                if( (variable.at(i) > fourth_largest_value) && (variable.at(i) < third_largest_value) ){

                        fourth_largest_value = variable.at(i);

                }

        }


  return fourth_largest_value;


  }


}};




auto MET_function{[](const floats& MET_input){

  return MET_input;

}};




auto deltaRcheck_float{[](

const float& Object1_eta,
const float& Object1_phi,
const float& Object2_eta,
const float& Object2_phi
){

  float dR = sqrt(pow(Object1_eta - Object2_eta, 2) + pow(Object1_phi - Object2_phi, 2));
  return dR;

}};

auto deltaRcheck_4floats{[](

const floats& Object1_eta,
const floats& Object1_phi,
const floats& Object2_eta,
const floats& Object2_phi

){

  floats dR = sqrt(pow(Object1_eta - Object2_eta, 2) + pow(Object1_phi - Object2_phi, 2));
  return dR;

}};


auto deltaRcheck_floats{[](const floats& Object1_eta, const floats& Object1_phi, const floats& Object2_eta, const floats& Object2_phi) {

  floats min_dRs{};

  if(Object2_phi.size() > 1){

  	transform(Object1_eta.begin(), Object1_eta.end(), Object1_phi.begin(), std::back_inserter(min_dRs), [&](float Object1_eta, float Object1_phi) { return std::min(deltaR(Object1_eta, Object1_phi, Object2_eta.at(0), Object2_phi.at(0)), deltaR(Object1_eta, Object1_phi, Object2_eta.at(1), Object2_phi.at(1))); });

  }
  else{

	transform(Object1_eta.begin(), Object1_eta.end(), Object1_phi.begin(), std::back_inserter(min_dRs), [&](float Object1_eta, float Object1_phi) { return deltaR(Object1_eta, Object1_phi, Object2_eta.at(0), Object2_phi.at(0)); });

  }

  return min_dRs;
 
}};



auto deltaRcheck_Top_function{[](

const doubles& Object1_phi_Selection,
const doubles& Object1_eta_Selection,
const float& Object2_eta_Selection,
const float& Object2_phi_Selection
){

  doubles dR = sqrt(pow(Object1_eta_Selection - Object2_eta_Selection, 2) + pow(Object1_phi_Selection - Object2_phi_Selection, 2));
  return dR;

}};



auto deltaRcheck_WTop_function{[](

const floats& Object1_phi_Selection,
const floats& Object1_eta_Selection,
const doubles& Object2_eta_Selection,
const doubles& Object2_phi_Selection
){

  doubles dR_vec{};

  for(int i = 0; i < Object1_phi_Selection.size(); i++){

  	double dR = sqrt(pow(Object1_eta_Selection.at(i) - Object2_eta_Selection.at(0), 2) + pow(Object1_phi_Selection.at(i) - Object2_phi_Selection.at(0), 2));
  	dR_vec.push_back(dR);

  }

  return dR_vec;


}};


auto deltaRcheck_W_function{[](

const doubles& Object1_phi_Selection,
const doubles& Object1_eta_Selection,
const doubles& Object2_eta_Selection,
const doubles& Object2_phi_Selection
){

  doubles dR = sqrt(pow(Object1_eta_Selection - Object2_eta_Selection, 2) + pow(Object1_phi_Selection - Object2_phi_Selection, 2));
  return dR;

}};



auto deltaRcheck_W_function2{[](

const doubles& Object1_phi_Selection,
const doubles& Object1_eta_Selection,
const float& Object2_eta_Selection,
const float& Object2_phi_Selection
){
 
  doubles dR = sqrt(pow(Object1_eta_Selection - Object2_eta_Selection, 2) + pow(Object1_phi_Selection - Object2_phi_Selection, 2));
  return dR;

}};



auto DeltaPhi_function{[](

const floats& Object1_phi_Selection,
const floats& Object2_phi_Selection

){

  floats dPhi = abs(Object1_phi_Selection - Object2_phi_Selection);
  return dPhi;

}};




auto DeltaPhi_function2{[](

const doubles& Object1_phi_Selection,
const doubles& Object2_phi_Selection

){

  doubles dPhi = abs(Object1_phi_Selection - Object2_phi_Selection);
  return dPhi;


}};



auto DeltaPhi_function3{[](

const doubles& Object1_phi_Selection,
const floats& Object2_phi_Selection

){

  doubles dPhi = abs(Object1_phi_Selection - Object2_phi_Selection);
  return dPhi;


}};



auto DeltaPhi_function4{[](

const floats& Object1_phi,
const doubles& Object2_phi

){

 doubles dPhi_vec{};

 for(int i = 0; i < Object1_phi.size(); i++){

 	double dPhi = Object1_phi.at(i) - Object2_phi.at(0);
	dPhi_vec.push_back(dPhi);

 }

 return dPhi_vec;

}};



auto DeltaPhi_doublesandfloat{[](

const doubles& Object1_phi,
const float& Object2_phi

){

  doubles dPhi = abs(Object1_phi - Object2_phi);
  return dPhi;

}};

auto DeltaPhi_floatandfloat{[](

const float& Object1_phi,
const float& Object2_phi

){

  double dPhi = abs(Object1_phi - Object2_phi);
  return dPhi;

}};




auto tight_jets_function{[&year](

const floats& Jet_pt_Selection,
const floats& Jet_eta_Selection,
const ints& Jet_jetId_Selection,
const floats& dRJet_lep){

  int JetId;

  if(year == "2016"){JetId = 1;} //1 is loose 
  else if(year == "2017" || year == "2018"){JetId = 2;} //2 is tight
  else{std::cout << "Choose a year out of 2016, 2017 or 2018" << std::endl;}

  return

  Jet_pt_Selection > 30 &&
  Jet_eta_Selection < 4.7 &&
  Jet_jetId_Selection >= JetId &&
  dRJet_lep > 0.4;


}};


auto jet_selection_function{[](const ints& tight_jets) {

  auto njet{count_if(tight_jets.begin(), tight_jets.end(), [](int i) { return i; })};
  return njet >= 4 && njet <= 6;

}};


auto SumSquared2LeadingJets_pT{[](

const float& LeadingJetPt,
const float& SubleadingJetPt

){

  double SumSquaredPt = pow(LeadingJetPt + SubleadingJetPt, 2);
  return SumSquaredPt;


}};


auto JetPtSum{[](

const float& LeadingJetPt,
const float& SubleadingJetPt,
const float& ThirdJetPt,
const float& FourthJetPt

){

  float JetPtSum = LeadingJetPt + SubleadingJetPt + ThirdJetPt + FourthJetPt;
  return JetPtSum;

}};


auto JetEtaSum{[](

const float& LeadingJetEta,
const float& SubleadingJetEta,
const float& ThirdJetEta,
const float& FourthJetEta

){

  float JetEtaSum = LeadingJetEta + SubleadingJetEta + ThirdJetEta + FourthJetEta;
  return JetEtaSum;

}};


auto JetPhiSum{[](

const float& LeadingJetPhi,
const float& SubleadingJetPhi,
const float& ThirdJetPhi,
const float& FourthJetPhi

){

  float JetPhiSum = LeadingJetPhi + SubleadingJetPhi + ThirdJetPhi + FourthJetPhi;
  return JetPhiSum;

}};



auto LepPtSum{[](

const float& LeadingLepPt,
const float& SubleadingLepPt

){

  float LepPtSum = LeadingLepPt + SubleadingLepPt;
  return LepPtSum;

}};



auto LepEtaSum{[](

const float& LeadingLepEta,
const float& SubleadingLepEta

){

  float LepEtaSum = LeadingLepEta + SubleadingLepEta;
  return LepEtaSum;

}};



auto LepPhiSum{[](

const float& LeadingLepPhi,
const float& SubleadingLepPhi

){

  float LepPhiSum = LeadingLepPhi + SubleadingLepPhi;
  return LepPhiSum;

}};



auto HT{[](const float& Pt){
 
  float HT = abs(Pt);
  return HT;

}};


auto HT_double{[](const doubles& Pt){

  doubles HT = abs(Pt);
  return HT;

}};


auto HT_floats{[](const floats& Pt){

  floats HT = abs(Pt);
  return HT;

}};


auto TotJetHT{[](

const float& LeadingJetHT,
const float& SubleadingJetHT,
const float& ThirdJetHT,
const float& FourthJetHT

){
  
  float TotJetHT = LeadingJetHT + SubleadingJetHT + ThirdJetHT + FourthJetHT;
  return TotJetHT;

}};



auto TotLepHT{[](

const float& LeadingLeptonHT,
const float& SubleadingLeptonHT

){

  float TotLepHT = LeadingLeptonHT + SubleadingLeptonHT;
  return TotLepHT;

}};

auto TotHTOverTotpT{[](const float& TotHT, const float& TotpT){

  float TotHTOverTotpT = TotHT / TotpT;
  return TotHTOverTotpT;


}};

auto TotHTOverTotpT_floats{[](const floats& TotHT, const floats& TotpT){

  floats TotHTOverTotpT = TotHT / TotpT;
  return TotHTOverTotpT;


}};

auto InvMass_AllJets{[](

const float& LeadingJetPt,
const float& SubleadingJetPt,
const float& ThirdJetPt,
const float& FourthJetPt,
const float& LeadingJetEta,
const float& SubleadingJetEta,
const float& ThirdJetEta,
const float& FourthJetEta,
const float& LeadingJetPhi,
const float& SubleadingJetPhi,
const float& ThirdJetPhi,
const float& FourthJetPhi,
const float& LeadingJetMass,
const float& SubleadingJetMass,
const float& ThirdJetMass,
const float& FourthJetMass,
const unsigned int& nJet 

){

  TLorentzVector Jet1 = {};
  TLorentzVector Jet2 = {};
  TLorentzVector Jet3 = {};
  TLorentzVector Jet4 = {};


  Jet1.SetPtEtaPhiM(LeadingJetPt, LeadingJetEta, LeadingJetPhi, LeadingJetMass);
  Jet2.SetPtEtaPhiM(SubleadingJetPt, SubleadingJetEta, SubleadingJetPhi, SubleadingJetMass);
  Jet3.SetPtEtaPhiM(ThirdJetPt, ThirdJetEta, ThirdJetPhi, ThirdJetMass);
  Jet4.SetPtEtaPhiM(FourthJetPt, FourthJetEta, FourthJetPhi, FourthJetMass);

  float InvMassAllJets = (Jet1 + Jet2 + Jet3 + Jet4).M();

  return InvMassAllJets;

}};

auto InvMass_3Jets{[](

const float& LeadingJetPt,
const float& SubleadingJetPt,
const float& ThirdJetPt,
const float& LeadingJetEta,
const float& SubleadingJetEta,
const float& ThirdJetEta,
const float& LeadingJetPhi,
const float& SubleadingJetPhi,
const float& ThirdJetPhi,
const float& LeadingJetMass,
const float& SubleadingJetMass,
const float& ThirdJetMass,
const unsigned int nJet

){

  TLorentzVector Jet1 = {};
  TLorentzVector Jet2 = {};
  TLorentzVector Jet3 = {};

  Jet1.SetPtEtaPhiM(LeadingJetPt, LeadingJetEta, LeadingJetPhi, LeadingJetMass);
  Jet2.SetPtEtaPhiM(SubleadingJetPt, SubleadingJetEta, SubleadingJetPhi, SubleadingJetMass);
  Jet3.SetPtEtaPhiM(ThirdJetPt, ThirdJetEta, ThirdJetPhi, ThirdJetMass);

  float InvMass3Jets = (Jet1 + Jet2 + Jet3).M();

  return InvMass3Jets;

}};


auto bjet_id{[](const ints& tight_jets, const floats& btags, const floats& etas) {
     
	return tight_jets && (btags > 0.8838f) && (etas < MaxTrackerEta);
}};


auto nonbjet_id{[](const ints& tight_jets, const floats& btags, const floats& etas) {

  return tight_jets && (btags == 0) && (etas < MaxTrackerEta);


}};


auto bjet_cut{[](const ints& bjets) {

        const auto nbjet{std::count_if(bjets.begin(), bjets.end(), [](int i) { return i; })};
        return nbjet >= 1 && nbjet <= 2;

}};

//Lambda functions between lines 1211 and 1242 are only for calculating b-tagging efficiency
//For the numerators
auto BTAGEFF_bjet_id_WP{[](const ints& tight_jets, const floats& btags, const floats& etas, const ints& Jet_partonFlavour) {

	return abs(Jet_partonFlavour) == 5 && btags > 0.8838f && abs(etas) < MaxTrackerEta;
	
}};


auto BTAGEFF_charm_id_WP{[](const ints& tight_jets, const floats& btags, const floats& etas, const ints& Jet_partonFlavour) {

       return abs(Jet_partonFlavour) == 4 && btags > 0.8838f && abs(etas) < MaxTrackerEta;

}};



auto BTAGEFF_lightjets_id_WP{[](const ints& tight_jets, const floats& btags, const floats& etas, const ints& Jet_partonFlavour) {
        
      return abs(Jet_partonFlavour) > 0 && abs(Jet_partonFlavour) < 4 && btags > 0.8838f && abs(etas) < MaxTrackerEta;

}};



auto BTAGEFF_gluon_id_WP{[](const ints& tight_jets, const floats& btags, const floats& etas, const ints& Jet_partonFlavour) {
                
      return abs(Jet_partonFlavour) == 21 && btags > 0.8838f && abs(etas) < MaxTrackerEta;

}};


auto BTAGEFF_nonbjet_id_WP{[](const ints& tight_jets, const floats& btags, const floats& etas, const ints& Jet_partonFlavour){

    return abs(Jet_partonFlavour) != 5 && btags > 0.8838f && abs(etas) < MaxTrackerEta;

}};



//For the denominators
auto BTAGEFF_bjet_id{[](const ints& tight_jets, const floats& etas, const ints& Jet_partonFlavour) {

	return abs(Jet_partonFlavour) == 5 && abs(etas) < MaxTrackerEta;

}};



auto BTAGEFF_charm_id{[](const ints& tight_jets, const floats& etas, const ints& Jet_partonFlavour) {

	return abs(Jet_partonFlavour) == 4 && abs(etas) < MaxTrackerEta;

}};




auto BTAGEFF_lightjets_id{[](const ints& tight_jets, const floats& etas, const ints& Jet_partonFlavour) {

	return abs(Jet_partonFlavour) > 0 && abs(Jet_partonFlavour) < 4 && abs(etas) < MaxTrackerEta;

}};



auto BTAGEFF_gluon_id{[](const ints& tight_jets, const floats& etas, const ints& Jet_partonFlavour) {

        return abs(Jet_partonFlavour) == 21 && abs(etas) < MaxTrackerEta;

}};



auto BTAGEFF_nonbjet_id{[](const ints& tight_jets, const floats& etas, const ints& Jet_partonFlavour){

	 return abs(Jet_partonFlavour) != 5 && abs(etas) < MaxTrackerEta;

}};





auto numberofbjets{[](const ints& bjets) {

        const auto nbjet{std::count_if(bjets.begin(), bjets.end(), [](int i) { return i; })};
        return nbjet;

}};



auto bjet_variable{[](

const floats& Jet_variable,
const unsigned int& nJet,
const ints& lead_bjet

){


  floats vec{};

  for(int i = 0; i < nJet; i++){
        if(lead_bjet.at(i) == 1){ 
		vec.push_back(Jet_variable.at(i));
	}

  }

  return vec;

}};



auto BLorentzVector{[](

const floats& bjet_pt,
const floats& bjet_eta,
const floats& bjet_phi,
const floats& bjet_mass

){

  auto BJets = TLorentzVector{};

  for(int i = 0; i < bjet_pt.size(); i++){

	auto Vec = TLorentzVector{};
	Vec.SetPtEtaPhiM(bjet_pt.at(i), bjet_eta.at(i), bjet_phi.at(i), bjet_mass.at(i));
	BJets += Vec;

  }



  return BJets;

}};


auto LeadingBJetOutputDiscriminant{[](

const float& LeadingJetpT,
const floats& Jet_btagCSVV2,
const ints& tight_jets,
const floats& Jet_eta_Selection

){

  return LeadingJetpT && (Jet_btagCSVV2  > 0.8838) && tight_jets && (abs(Jet_eta_Selection) < MaxTrackerEta);

}};

auto SubleadingBJetOutputDiscriminant{[](

const float& SubleadingJetpT,
const floats& Jet_btagCSVV2,
const ints& tight_jets,
const floats& Jet_eta_Selection

){

  return SubleadingJetpT && (Jet_btagCSVV2  > 0.8838) && tight_jets && (abs(Jet_eta_Selection) < MaxTrackerEta);

}};

auto ThirdBJetOutputDiscriminant{[](

const float& ThirdJetpT,
const floats& Jet_btagCSVV2,
const ints& tight_jets,
const floats& Jet_eta_Selection

){

  return ThirdJetpT && (Jet_btagCSVV2  > 0.8838) && tight_jets && (abs(Jet_eta_Selection) < MaxTrackerEta);

}};

auto FourthBJetOutputDiscriminant{[](

const float& FourthJetpT,
const floats& Jet_btagCSVV2,
const ints& tight_jets,
const floats& Jet_eta_Selection

){

  return FourthJetpT && (Jet_btagCSVV2  > 0.8838) && tight_jets && (abs(Jet_eta_Selection) < MaxTrackerEta);

}};



auto BJetOutputDiscriminant{[](

const ints& BJetBTags,
const floats& Jet_btagCSVV2 
){

  floats btagoutput{};

  for(int i = 0; i < BJetBTags.size(); i++){
	if(BJetBTags.at(i) != 0){
		btagoutput.push_back(Jet_btagCSVV2.at(i));
	}
  }

  return btagoutput;

}};



// W mass cut

constexpr float W_MASS = 80.385f;
constexpr float W_MASS_CUT = 20.f;

auto find_lead_mask{[](const ints& mask, const floats& vals) {
  

  const auto masked_vals{mask * vals};
  const auto max_idx{boost::numeric_cast<size_t>(std::distance(masked_vals.begin(), max_element(masked_vals.begin(), masked_vals.end())))};
  ints lead_mask(masked_vals.size(), 0); // must be ()
  lead_mask.at(max_idx) = 1;
  return lead_mask;


}};


auto find_w_pair{[](const floats& pts, const floats& etas, const floats& phis, const floats& ms, const ints& tight_jets, const ints& lead_bjet) {

double w_reco_mass{std::numeric_limits<double>::infinity()};
size_t jet_index_1{std::numeric_limits<size_t>::max()};
size_t jet_index_2{std::numeric_limits<size_t>::max()};
const size_t njets{pts.size()};


for (size_t i{0}; i < njets; ++i){
	for (size_t j{i + 1}; j < njets; ++j)
            {
                if (tight_jets[i] != 0 && tight_jets[j] != 0
                    && lead_bjet[i] != 1 && lead_bjet[j] != 1)
                {
                    continue;
                }

                auto jet1{TLorentzVector{}};
                auto jet2{TLorentzVector{}};
                jet1.SetPtEtaPhiM(pts.at(i), etas.at(i), phis.at(i), ms.at(i));
                jet2.SetPtEtaPhiM(pts.at(j), etas.at(j), phis.at(j), ms.at(j));

                if (const double reco_mass{(jet1 + jet2).M()}; std::abs(W_MASS - reco_mass) < std::abs(W_MASS - w_reco_mass))
                {
                    w_reco_mass = reco_mass;
                    jet_index_1 = i;
                    jet_index_2 = j;
                }
            }
        }

        ints w_pair(njets, 0);
        w_pair.at(jet_index_1) = 1;
        w_pair.at(jet_index_2) = 1;
        return w_pair;
}};



auto w_mass_cut{[](const float& w_mass) {
		
  return ( abs(w_mass - W_MASS) < W_MASS_CUT );

}};



auto w_mass_cut_ZPlusJetsCR{[](const float& w_mass, const float& MET_sumEt) {
	
  return abs(w_mass - W_MASS) > W_MASS_CUT && (MET_sumEt < 50);

}};



auto WPairJet1{[](const floats& pts, const floats& etas, const floats& phis, const floats& ms, const ints& tight_jets, const ints& lead_bjet) {


double w_reco_mass{std::numeric_limits<double>::infinity()};
size_t jet_index_1{std::numeric_limits<size_t>::max()};
size_t jet_index_2{std::numeric_limits<size_t>::max()};
const size_t njets{pts.size()};

auto jet1{TLorentzVector{}};
auto jet2{TLorentzVector{}};


for (size_t i{0}; i < njets; ++i){
        for (size_t j{i + 1}; j < njets; ++j)
            {
                if (tight_jets[i] != 0 && tight_jets[j] != 0
                    && lead_bjet[i] != 1 && lead_bjet[j] != 1)
                {
                    continue;
                }

                jet1.SetPtEtaPhiM(pts.at(i), etas.at(i), phis.at(i), ms.at(i));
                jet2.SetPtEtaPhiM(pts.at(j), etas.at(j), phis.at(j), ms.at(j));

                if (const double reco_mass{(jet1 + jet2).M()}; std::abs(W_MASS - reco_mass) < std::abs(W_MASS - w_reco_mass))
                {
                    w_reco_mass = reco_mass;
                    jet_index_1 = i;
                    jet_index_2 = j;
                }
            }
        }

	jet1.SetPtEtaPhiM(pts.at(jet_index_1), etas.at(jet_index_1), phis.at(jet_index_1), ms.at(jet_index_1));
	jet2.SetPtEtaPhiM(pts.at(jet_index_2), etas.at(jet_index_2), phis.at(jet_index_2), ms.at(jet_index_2));
        return jet1;
    
}};


auto WPairJet2{[](const floats& pts, const floats& etas, const floats& phis, const floats& ms, const ints& tight_jets, const ints& lead_bjet) {


double w_reco_mass{std::numeric_limits<double>::infinity()};
size_t jet_index_1{std::numeric_limits<size_t>::max()};
size_t jet_index_2{std::numeric_limits<size_t>::max()};
const size_t njets{pts.size()};

auto jet1{TLorentzVector{}};
auto jet2{TLorentzVector{}};


for (size_t i{0}; i < njets; ++i){
        for (size_t j{i + 1}; j < njets; ++j)
            {
                if (tight_jets[i] != 0 && tight_jets[j] != 0
                    && lead_bjet[i] != 1 && lead_bjet[j] != 1)
                {
                    continue;
                }

                jet1.SetPtEtaPhiM(pts.at(i), etas.at(i), phis.at(i), ms.at(i));
                jet2.SetPtEtaPhiM(pts.at(j), etas.at(j), phis.at(j), ms.at(j));

                if (const double reco_mass{(jet1 + jet2).M()}; std::abs(W_MASS - reco_mass) < std::abs(W_MASS - w_reco_mass))
                {
                    w_reco_mass = reco_mass;
                    jet_index_1 = i;
                    jet_index_2 = j;
                }
            }
        }

	jet1.SetPtEtaPhiM(pts.at(jet_index_1), etas.at(jet_index_1), phis.at(jet_index_1), ms.at(jet_index_1));
        jet2.SetPtEtaPhiM(pts.at(jet_index_2), etas.at(jet_index_2), phis.at(jet_index_2), ms.at(jet_index_2));
        return jet2;
    
}};



auto TLorentzVectorMass{[](const TLorentzVector& object){

  doubles vec{};
  vec.push_back(object.M());
  return vec;

}};


auto TLorentzVectorMass_float{[](const TLorentzVector& object){
  
  floats vec{};
  vec.push_back(object.M());
  return vec;

}};


auto TLorentzVectorPt{[](const TLorentzVector& object){

  doubles vec{};
  vec.push_back(object.Pt());
  return vec;

}};



auto TLorentzVectorPt_float{[](const TLorentzVector& object){
  
  floats vec{};
  vec.push_back(object.Pt());
  return vec;

}};


auto TLorentzVectorPhi{[](const TLorentzVector& object){

  doubles vec{};
  vec.push_back(object.Phi());
  return vec;

}};


auto TLorentzVectorPhi_float{[](const TLorentzVector& object){

  floats vec{};
  vec.push_back(object.Phi());
  return vec;

}};


auto TLorentzVectorEta{[](const TLorentzVector& object){

  doubles vec{};
  vec.push_back(object.Eta());
  return vec;

}};



auto TLorentzVectorEta_float{[](const TLorentzVector& object){

  floats vec{};
  vec.push_back(object.Eta());
  return vec;

}};







constexpr float Z_MASS{91.1876f};
constexpr float Z_MASS_CUT{20.f};

auto z_mass_cut{[](const float& z_mass) {

  return abs(z_mass - Z_MASS) < Z_MASS_CUT;

}};


auto RecoZ{[](

const unsigned int& nLepton,
const float& LeadingleptonPt,
const float& LeadingleptonEta,
const float& LeadingleptonPhi,
const float& LeadingleptonMass,
const float& SubleadingleptonPt,
const float& SubleadingleptonEta,
const float& SubleadingleptonPhi,
const float& SubleadingleptonMass

){

  TLorentzVector ZBoson = {};
  TLorentzVector LeadingLepton = {};
  TLorentzVector SubleadingLepton = {};

  LeadingLepton.SetPtEtaPhiM(LeadingleptonPt, LeadingleptonEta, LeadingleptonPhi, LeadingleptonMass);
  SubleadingLepton.SetPtEtaPhiM(SubleadingleptonPt, SubleadingleptonEta, SubleadingleptonPhi, SubleadingleptonMass);

  ZBoson = LeadingLepton + SubleadingLepton;

  return ZBoson;

}};


auto RecoZHT{[](const doubles& RecoZPt){

  doubles RecoZHT = abs(RecoZPt);
  return RecoZHT;

}};


auto RecoWHT{[](const floats& RecoWPt){

  floats RecoWHT = abs(RecoWPt);
  return RecoWHT;

}};

auto WLorentzVector{[](

const floats& w_pair_pt,
const floats& w_pair_eta, 
const floats& w_pair_phi, 
const float& w_mass, 
const ints& w_reco_jets

){

  const auto nRecoWBosons{std::count_if(w_reco_jets.begin(), w_reco_jets.end(), [](int i) { return i; })};

  auto RecoW = TLorentzVector{};
  
  for(int i = 0; i < nRecoWBosons; i++){

	  auto Vec = TLorentzVector{};
	  Vec.SetPtEtaPhiM(w_pair_pt.at(i), w_pair_eta.at(i), w_pair_phi.at(i), w_mass);
	  RecoW += Vec;

  }


  return RecoW;

}};

constexpr float TOP_MASS = 173.3;

auto top_reconstruction_function{[](

const floats& bjets_pt,
const floats& bjets_eta,
const floats& bjets_phi,
const floats& bjets_mass,
const unsigned int& nJet,
const floats& w_pair_pt,
const floats& w_pair_eta,
const floats& w_pair_phi,
const float& w_mass, 
const ints& w_reco_jets

){

  auto reco_top = TLorentzVector{}; 
  auto BJets = TLorentzVector{};
  auto RecoW = TLorentzVector{};

  double top_reco_mass = std::numeric_limits<double>::infinity();
  size_t index_1{std::numeric_limits<size_t>::max()};
  const size_t num{w_pair_pt.size()};

  for(int i = 0; i < num; i++){


  BJets.SetPtEtaPhiM(bjets_pt.at(0), bjets_eta.at(0), bjets_phi.at(0), bjets_mass.at(0));
  RecoW.SetPtEtaPhiM(w_pair_pt.at(i), w_pair_eta.at(i), w_pair_phi.at(i), w_mass);
		
  const double reco_mass = (RecoW + BJets).M(); 

  if(abs(TOP_MASS - reco_mass) < abs(TOP_MASS - top_reco_mass)){

	  top_reco_mass = reco_mass;
	  index_1 = i;

  }


  }


  BJets.SetPtEtaPhiM(bjets_pt.at(0), bjets_eta.at(0), bjets_phi.at(0), bjets_mass.at(0));
  RecoW.SetPtEtaPhiM(w_pair_pt.at(index_1), w_pair_eta.at(index_1), w_pair_phi.at(index_1), w_mass);
  reco_top = RecoW + BJets;	


  return reco_top;

}};



auto TotalHT_System{[](

const doubles& RecoZHT,
const floats& RecoWHT,
const doubles& Top_HT,
const float& TotLepHT,
const float& TotJetHT

){

  floats TotalHTSystem = RecoZHT + RecoWHT.at(0) + Top_HT + TotLepHT + TotJetHT;
  return TotalHTSystem;

}};

auto TotalPt_System{[](

const doubles& RecoZPt,
const floats& RecoWPt,
const doubles& Top_Pt,
const float& TotLepPt,
const float& TotJetPt

){

  floats TotalPtSystem = RecoZPt + RecoWPt.at(0) + Top_Pt + TotLepPt + TotJetPt;
  return TotalPtSystem;

}};


auto TotalEta_System{[](

const doubles& RecoZEta,
const floats& RecoWEta,
const doubles& Top_Eta,
const float& TotLepEta,
const float& TotJetEta

){

  doubles TotalEtaSystem = RecoZEta + RecoWEta.at(0) + Top_Eta + TotLepEta + TotJetEta;
  return TotalEtaSystem;

}};


auto TotalPhi_System{[](

const doubles& RecoZPhi,
const floats& RecoWPhi,
const doubles& Top_Phi,
const float& TotLepPhi,
const float& TotJetPhi

){

  doubles TotalPhiSystem = RecoZPhi + RecoWPhi.at(0) + Top_Phi + TotLepPhi + TotJetPhi;
  return TotalPhiSystem;


}};


//Minimum delta R between the Z boson candidate and any jet
auto MinDeltaR{[](

const unsigned int& nJet,
const doubles& RecoZPhi,
const doubles& RecoZEta,
const floats& Jet_Phi_Selection,
const floats& Jet_eta_Selection
){

    doubles output_vec;
  
    for(int i; i < nJet; i++){

    	double DeltaR = sqrt(pow(RecoZPhi.at(i) - Jet_Phi_Selection.at(i), 2) + pow(RecoZEta.at(i) - Jet_eta_Selection.at(i), 2));
    	double DeltaR2 = sqrt(pow(RecoZPhi.at(i+1) - Jet_Phi_Selection.at(i+1), 2) + pow(RecoZEta.at(i+1) - Jet_eta_Selection.at(i+1), 2));

    	double Output = (DeltaR2 < DeltaR) ? DeltaR2 : DeltaR;  
    	output_vec.push_back(Output);

    }

    return output_vec;

}};


//Minimum delta phi between the Z boson candidate and any jet
auto MinDeltaPhi{[](

const unsigned int& nJet,
const doubles& RecoZPhi,
const floats& Jet_Phi_Selection
){


  double output;
  doubles output_vec{};

  for(int i; i < nJet; i++){

    double DeltaPhi = RecoZPhi.at(i) - Jet_Phi_Selection.at(i);
    double DeltaPhi2 = RecoZPhi.at(i+1) - Jet_Phi_Selection.at(i+1);

    output = (DeltaPhi2 < DeltaPhi) ? DeltaPhi2 : DeltaPhi;
    output_vec.push_back(output);

  }

  return output_vec;

}};


auto dR_LeadingLepton_LeadingBJet{[](

const floats& bjeteta,
const float& LeadingLeptonEta,
const floats& bjetphi,
const float& LeadingLeptonPhi

){

  doubles DeltaR = sqrt(pow(LeadingLeptonPhi - bjetphi, 2) + pow(LeadingLeptonEta - bjeteta, 2));
  return DeltaR;

}};



auto dR_SubleadingLepton_LeadingBJet{[](

const floats& bjeteta,
const float& SubleadingLeptonEta,
const floats& bjetphi,
const float& SubleadingLeptonPhi

){

  doubles DeltaR = sqrt(pow(SubleadingLeptonPhi - bjetphi, 2) + pow(SubleadingLeptonEta - bjeteta, 2));
  return DeltaR;

}};


auto DeltaPhi_Lepton_BJet{[](

const floats& Jet_phi_Selection,
const float& LeptonPhi

){

  doubles DeltaPhi = LeptonPhi - Jet_phi_Selection;
  return DeltaPhi;

}};



auto TransverseWMass{[](

const double& dPhi_j1j2,
const doubles& WPairJet1Pt,
const doubles& WPairJet2Pt

){

  doubles mtW = sqrt(2 * WPairJet1Pt * WPairJet2Pt * (1 - cos(dPhi_j1j2)) );
  return mtW;

}};


auto filter_function{[](

const bool& Flag_goodVertices_Selection, 
const bool& Flag_globalSuperTightHalo2016Filter_Selection, 
const bool& Flag_HBHENoiseFilter_Selection, 
const bool& Flag_HBHENoiseIsoFilter_Selection, 
const bool& Flag_EcalDeadCellTriggerPrimitiveFilter_Selection, 
const bool& Flag_BadPFMuonFilter_Selection, 
const bool& Flag_BadChargedCandidateFilter_Selection, 
const bool& Flag_ecalBadCalibFilter_Selection, 
const bool& Flag_eeBadScFilter_Selection
)-> bool{


return  

Flag_goodVertices_Selection > 0 || 
Flag_globalSuperTightHalo2016Filter_Selection > 0 || 
Flag_HBHENoiseFilter_Selection > 0 || 
Flag_HBHENoiseIsoFilter_Selection > 0 || 
Flag_EcalDeadCellTriggerPrimitiveFilter_Selection > 0 || 
Flag_BadPFMuonFilter_Selection > 0 || 
Flag_BadChargedCandidateFilter_Selection > 0 || 
Flag_ecalBadCalibFilter_Selection > 0 || 
Flag_eeBadScFilter_Selection > 0;

}};



//Lambda function for jet smearing simulation corrections (scaling method)
//Twiki link: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures

//Reading the JER and SF text files
std::string FileNameJetSmear;


auto RowReader2{[&FileNameJetSmear, &year](

const int& LineSpecified, 
const bool& sigmaJER, 
const bool& SF, 
const bool& up, 
const bool& down,
const floats& Jet_eta,
const floats& Jet_rho,
const floats& Jet_pt) { 


  float Col1, Col2, Col3, Col4, Col5, Col6, Col7, Col8, Col9, Col10, Col11;
  
  if(year == "2016"){

  	if(sigmaJER == true && SF == false && up == false && down == false){FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2016/Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt";}
  	else if(sigmaJER == false && SF == true && up == false && down == false){FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2016/Summer16_25nsV1_MC_SF_AK4PFchs.txt";}
  	else if(sigmaJER == false && SF == false && up == true && down == false){FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2016/Summer16_25nsV1_MC_SF_AK4PFchs.txt";}
  	else if(sigmaJER == false && SF == false && up == false && down == true){FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2016/Summer16_25nsV1_MC_SF_AK4PFchs.txt";}
  	else{std::cout << "Please enter an appropriate file name" << std::endl;}

  }
  else if(year == "2017"){

	if(sigmaJER == true && SF == false && up == false && down == false){FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2017/Fall17_V3_MC_PtResolution_AK4PFchs.txt";}
        else if(sigmaJER == false && SF == true && up == false && down == false){FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2017/Fall17_V3_MC_SF_AK4PFchs.txt";}
        else if(sigmaJER == false && SF == false && up == true && down == false){FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2017/Fall17_V3_MC_SF_AK4PFchs.txt";}
        else if(sigmaJER == false && SF == false && up == false && down == true){FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2017/Fall17_V3_MC_SF_AK4PFchs.txt";}
        else{std::cout << "Please enter an appropriate file name" << std::endl;}

  }
  else if(year == "2018"){

  	if(sigmaJER == true && SF == false && up == false && down == false){FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2018/Autumn18_V1_MC_PtResolution_AK4PFchs.txt";}
  	else if(sigmaJER == false && SF == true && up == false && down == false){FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2018/Autumn18_V1_MC_SF_AK4PFchs.txt";}
  	else if(sigmaJER == false && SF == false && up == true && down == false){FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2018/Autumn18_V1_MC_SF_AK4PFchs.txt";}
  	else if(sigmaJER == false && SF == false && up == false && down == true){FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2018/Autumn18_V1_MC_SF_AK4PFchs.txt";}
  	else{std::cout << "Please enter an appropriate file name" << std::endl;}

  }
  else{std::cout << "The year can only be 2016, 2017 or 2018" << std::endl;}


  std::ifstream file;
  file.open(FileNameJetSmear);

  if (file.good())
  {
    std::string str = "";

    int line_number = 0;
	
	 while(getline(file, str) && line_number != LineSpecified){
		++line_number;
	}
	if(line_number == LineSpecified){
		if(sigmaJER == true && SF == false && up == false && down == false){
			file >> Col1;
                	file >> Col2;
                	file >> Col3;
                	file >> Col4;
                	file >> Col5;
                	file >> Col6;
                	file >> Col7;
                	file >> Col8;
                	file >> Col9;
                	file >> Col10;
                	file >> Col11;
		}
		else if(sigmaJER == false && SF == true && up == false && down == false){
			file >> Col1;
                        file >> Col2;
                        file >> Col3;
                        file >> Col4;
                        file >> Col5;
                        file >> Col6;
		}
		else if(sigmaJER == false && SF == false && up == true && down == false){
                        file >> Col1;
                        file >> Col2;
                        file >> Col3;
                        file >> Col4;
                        file >> Col5;
                        file >> Col6;
                }
		else if(sigmaJER == false && SF == false && up == false && down == true){
                        file >> Col1;
                        file >> Col2;
                        file >> Col3;
                        file >> Col4;
                        file >> Col5;
                        file >> Col6;
                }
		else{std::cout << "Please enter an appropriate file name" << std::endl;}

	}

  }
 
  file.close(); 
  
  bool Jet_eta_check = all_of(Jet_eta.begin(), Jet_eta.end(), [&Col1, &Col2](int i){return ( i > abs(Col1) && i < abs(Col2) );});
  bool Jet_rho_check = all_of(Jet_rho.begin(), Jet_rho.end(), [&Col3, &Col4](int i){return ( i > abs(Col3) && i < abs(Col4) );});
  bool Jet_pt_check = all_of(Jet_pt.begin(), Jet_pt.end(), [&Col6, &Col7](int i){return ( i > abs(Col6) && i < abs(Col7) );} );

 
  if(sigmaJER == true && SF == false && up == false && down == false){

        floats answer = sqrt( Col8*abs(Col8) / (Jet_pt*Jet_pt)+Col9*Col9*pow(Jet_pt,Col11)+Col10*Col10 );

  	if( Jet_eta_check == true && Jet_rho_check == true && Jet_pt_check == true ){
		return answer;
  	}
  	else{
		floats Vec1(answer.size(), 0.0);
		return Vec1;}

  }
  else if(sigmaJER == false && SF == true && up == false && down == false){
        
	floats Col4Vec = {};
        Col4Vec.push_back(Col4);

 	if(Jet_eta_check == true){
		return Col4Vec;

	}
	else{	floats Vec2(Col4Vec.size(), 0.0);
		return Vec2;}
  }
  else if(sigmaJER == false && SF == false && up == true && down == false){

        floats UpValueVec = {};
        float UpValue = Col6 - Col4;
        UpValueVec.push_back(UpValue);

	if(Jet_eta_check == true){
                return UpValueVec;

	}
        else{	floats Vec3(UpValueVec.size(), 0.0);
		return Vec3;}
  }
  else if(sigmaJER == false && SF == false && up == false && down == true){
	
        floats DownValueVec = {};
        float DownValue = Col4 - Col5;
        DownValueVec.push_back(DownValue);

        if(Jet_eta_check == true){
                return DownValueVec;
	}
        else{	floats Vec4(DownValueVec.size(), 0.0);
		return Vec4;}
  }
  else{std::cout << "bools cannot be all true or all false" << std::endl; std::cout << "sigmaJER = " << sigmaJER << std::endl; std::cout << "SF = " << SF << std::endl; std::cout << "up = " << up << std::endl; std::cout << "down = " << down << std::endl;} 



}}; 


auto linecounter{[&FileNameJetSmear, &year](const bool& sigmaJER, const bool& SF, const bool& up, const bool& down){ 

   int number_of_lines = 0;
   std::string line;


   if(year == "2016"){

        if(sigmaJER == true && SF == false && up == false && down == false){FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2016/Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt";}
        else if(sigmaJER == false && SF == true && up == false && down == false){FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2016/Summer16_25nsV1_MC_SF_AK4PFchs.txt";}
        else if(sigmaJER == false && SF == false && up == true && down == false){FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2016/Summer16_25nsV1_MC_SF_AK4PFchs.txt";}
        else if(sigmaJER == false && SF == false && up == false && down == true){FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2016/Summer16_25nsV1_MC_SF_AK4PFchs.txt";}
        else{std::cout << "Please enter an appropriate file name" << std::endl;}

  }
  else if(year == "2017"){

        if(sigmaJER == true && SF == false && up == false && down == false){FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2017/Fall17_V3_MC_PtResolution_AK4PFchs.txt";}
        else if(sigmaJER == false && SF == true && up == false && down == false){FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2017/Fall17_V3_MC_SF_AK4PFchs.txt";}
        else if(sigmaJER == false && SF == false && up == true && down == false){FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2017/Fall17_V3_MC_SF_AK4PFchs.txt";}
        else if(sigmaJER == false && SF == false && up == false && down == true){FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2017/Fall17_V3_MC_SF_AK4PFchs.txt";}
        else{std::cout << "Please enter an appropriate file name" << std::endl;}

  }
  else if(year == "2018"){

        if(sigmaJER == true && SF == false && up == false && down == false){FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2018/Autumn18_V1_MC_PtResolution_AK4PFchs.txt";}
        else if(sigmaJER == false && SF == true && up == false && down == false){FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2018/Autumn18_V1_MC_SF_AK4PFchs.txt";}
        else if(sigmaJER == false && SF == false && up == true && down == false){FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2018/Autumn18_V1_MC_SF_AK4PFchs.txt";}
        else if(sigmaJER == false && SF == false && up == false && down == true){FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2018/Autumn18_V1_MC_SF_AK4PFchs.txt";}
        else{std::cout << "Please enter an appropriate file name" << std::endl;}

  }
  else{std::cout << "The year can only be 2016, 2017 or 2018" << std::endl;} 


   std::ifstream myfile(FileNameJetSmear);

   while (getline(myfile, line))
        ++number_of_lines;
    	return number_of_lines;

}};




auto RowReader3{[&RowReader2, &linecounter](

const bool& SigmaJER, 
const bool& JetSmearScaleFactor, 
const bool& Up, 
const bool& Down,
const floats& Jet_eta, 
const floats& Jet_rho, 
const floats& Jet_pt
){

  for(int i = 0; i < linecounter(SigmaJER, JetSmearScaleFactor, Up, Down) + 1; i++){

		std::string quantity; 

		if(SigmaJER == true && JetSmearScaleFactor == false && Up == false && Down == false){quantity = "sigma JER";}
   		else if(SigmaJER == false && JetSmearScaleFactor == true && Up == false && Down == false){quantity = "SF";}
   		else if(SigmaJER == false && JetSmearScaleFactor == false && Up == true && Down == false){quantity = "SF (up variation)";}
   		else if(SigmaJER == false && JetSmearScaleFactor == false && Up == false && Down == true){quantity = "SF (down variation)";}
   		else{std::cout << "Please enter an appropriate file name" << std::endl;}

		return RowReader2(i, SigmaJER, JetSmearScaleFactor, Up, Down, Jet_eta, Jet_rho, Jet_pt) > 0.0;

  }


}};



//sigma_JER reads this file for 2017: https://github.com/cms-jet/JRDatabase/blob/master/textFiles/Fall17_V3_MC/Fall17_V3_MC_PtResolution_AK4PF.txt
auto sigma_JER{[&RowReader3](const floats& Jet_eta, const floats& Jet_rho,const floats& Jet_pt){

  bool SigmaJER = true;
  bool JetSmearScaleFactor = false;
  bool Up = false;
  bool Down = false;

  return RowReader3(SigmaJER, JetSmearScaleFactor, Up, Down, Jet_eta, Jet_rho, Jet_pt);

}};


auto sigma_JER_up{[&RowReader3](const floats& Jet_eta, const floats& Jet_rho,const floats& Jet_pt){

  bool SigmaJER = false;
  bool JetSmearScaleFactor = false;
  bool Up = true;
  bool Down = false;

  return RowReader3(SigmaJER, JetSmearScaleFactor, Up, Down, Jet_eta, Jet_rho, Jet_pt);

}};


auto sigma_JER_down{[&RowReader3](const floats& Jet_eta, const floats& Jet_rho,const floats& Jet_pt){

  bool SigmaJER = false;
  bool JetSmearScaleFactor = false;
  bool Up = false;
  bool Down = true;
  
  return RowReader3(SigmaJER, JetSmearScaleFactor, Up, Down, Jet_eta, Jet_rho, Jet_pt);

}};

//SJER reads this file for 2017: https://github.com/cms-jet/JRDatabase/blob/master/textFiles/Fall17_V3_MC/Fall17_V3_MC_SF_AK4PF.txt 
auto SJER_nominal{[&RowReader3](const floats& Jet_eta, const floats& Jet_rho, const floats& Jet_pt){

  bool SigmaJER = false;
  bool JetSmearScaleFactor = true;
  bool Up = false;
  bool Down = false;
 
  return RowReader3(SigmaJER, JetSmearScaleFactor, Up, Down, Jet_eta, Jet_rho, Jet_pt);


}};

auto SJER_up{[&RowReader3](const floats& Jet_eta, const floats& Jet_rho, const floats& Jet_pt){

  bool SigmaJER = false;
  bool JetSmearScaleFactor = false;
  bool Up = true;
  bool Down = false;
  
  return RowReader3(SigmaJER, JetSmearScaleFactor, Up, Down, Jet_eta, Jet_rho, Jet_pt);


}};

auto SJER_down{[&RowReader3](const floats& Jet_eta, const floats& Jet_rho, const floats& Jet_pt){

  bool SigmaJER = false;
  bool JetSmearScaleFactor = false;
  bool Up = false;
  bool Down = true;
  
  return RowReader3(SigmaJER, JetSmearScaleFactor, Up, Down, Jet_eta, Jet_rho, Jet_pt);


}};


auto dRCone_condition{[](

const floats& eta,
const floats& phi,
const floats& eta_ptcl,
const floats& phi_ptcl

){


   bool PhiSize = (phi.size() == phi_ptcl.size()) ? true : false;
   bool EtaSize = (eta.size() == eta_ptcl.size()) ? true : false;

   if(PhiSize == true && EtaSize == true){

        doubles dphi = phi - phi_ptcl;
        doubles deta = eta - eta_ptcl;
        doubles deltaR = sqrt( pow(dphi, 2) + pow(deta, 2) );
        const double RCone = 0.4;

	return deltaR < RCone / 2;
   }
   else{
   	ints OutputVec(eta.size(), 0);
	return OutputVec;
   } 


}};


auto pT_ptcl_condition{[](const floats& pT, const floats& pT_ptcl, const ints& sigma_JER){

 bool PtSize = (pT.size() == pT_ptcl.size()) ? true : false;

 if(PtSize == true){
 	return abs(pT - pT_ptcl) < 3 * sigma_JER.at(0) * pT;
 }
 else{
	if(pT.size() < pT_ptcl.size()){ints OutputVec(pT.size(), 0); return OutputVec;}
	else{ints OutputVec(pT_ptcl.size(), 0); return OutputVec;}

 }


}};


auto GreaterThanZero{[](const floats& sJER_nominal){

  return (sJER_nominal*sJER_nominal - 1);

}};


//Calculating the jet smearing correction factor using the hybrid method
auto MaxComparison{[&GreaterThanZero](const ints& sJER_nominal){

 floats MaximumFloats = sqrt(sJER_nominal*sJER_nominal - 1);
 int size = MaximumFloats.size();

 if(sJER_nominal.at(0) > 0){
	floats MaximumFloats = sqrt(sJER_nominal*sJER_nominal - 1);
 	return MaximumFloats;
 }
 else{
 	floats ZeroVec(size, 0.0);
	return ZeroVec;
 }

}};


auto JetSmearingFunction_HybridMethod{[&dRCone_condition, &pT_ptcl_condition, &MaxComparison](

const floats& pT,
const floats& eta,      
const floats& phi, 
const floats& pT_ptcl, 
const floats& eta_ptcl, 
const floats& phi_ptcl, 
const ints& sJER_nominal, 
const ints& sigma_JER){


  floats cJER;
 
  bool dRCone_check = all_of(dRCone_condition(eta, phi, eta_ptcl, phi_ptcl).begin(), dRCone_condition(eta, phi, eta_ptcl, phi_ptcl).end(), [](int i){return i > 0;});

  bool pT_ptcl_check = all_of(pT_ptcl_condition(pT, pT_ptcl, sigma_JER).begin(), pT_ptcl_condition(pT, pT_ptcl, sigma_JER).end(), [](int i){ return i != 0;});

  bool PtSize = (pT.size() == pT_ptcl.size()) ? true : false;


  if(dRCone_check == true && pT_ptcl_check == true && PtSize == true){

	floats cJER = 1 + ( (sJER_nominal.at(0) - 1) * ( (pT - pT_ptcl) / pT ) );

	floats cJER_vec{};

	for(int i = 0; i < cJER.size(); i++){

		if(cJER.at(i) > 0){ cJER_vec.push_back(cJER.at(i)); }
		else{cJER_vec.push_back(0.);}		

	}


	return cJER_vec;

  }
  else{
  	  float N = gRandom->Gaus(0, sJER_nominal.at(0));
          cJER = 1.0 + (N * MaxComparison(sJER_nominal));     
	  return cJER; 
  }
  


}};


auto ApplyCJER{[](

const floats& JetPt, 
const floats& JetEta,
const floats& JetPhi,
const floats& JetMass, 
const floats& cJER, 
const unsigned int& nJet

){


  std::vector<TLorentzVector> OutputVec{};

  for(int i = 0; i < nJet; i++){

    TLorentzVector JetFourMomentum_New{};
    float JetPt_new = JetPt.at(i) * cJER.at(0);
    float JetEta_new = JetEta.at(i) * cJER.at(0);
    float JetPhi_new = JetPhi.at(i) * cJER.at(0);
    float JetMass_new = JetMass.at(i) * cJER.at(0);


    JetFourMomentum_New.SetPtEtaPhiM(JetPt_new, JetEta_new, JetPhi_new, JetMass_new);

    OutputVec.push_back(JetFourMomentum_New);

  }

 return OutputVec; 

}};


//To get the properties of the smeared jet
auto GetSmearedJetPt{[](std::vector<TLorentzVector> SmearedJet4Momentum, const floats& JetPt){

 floats NewPtVec = {};

 for(int i = 0; i < JetPt.size(); i++){

        float NewPt = (SmearedJet4Momentum.at(i)).Pt();
 	NewPtVec.push_back(NewPt);

 }

 return NewPtVec;

}};


auto GetSmearedJetPhi{[](std::vector<TLorentzVector> SmearedJet4Momentum, const floats& JetPhi){

 floats NewPhiVec{};


 for(int i = 0; i < JetPhi.size(); i++){

	float NewPhi = (SmearedJet4Momentum.at(i)).Phi();
        NewPhiVec.push_back(NewPhi);

 }
 
 return NewPhiVec;

}};

 
auto GetSmearedJetEta{[](std::vector<TLorentzVector> SmearedJet4Momentum, const floats& JetEta){


 floats NewEtaVec = {};

 for(int i = 0; i < JetEta.size(); i++){

        float NewEta = (SmearedJet4Momentum.at(i)).Eta();
        NewEtaVec.push_back(NewEta);

 }

 return NewEtaVec;

}};


auto GetSmearedJetMass{[](std::vector<TLorentzVector> SmearedJet4Momentum, const floats& JetMass){


 floats NewMassVec = {};

 for(int i = 0; i < JetMass.size(); i++){

	float NewPhi = (SmearedJet4Momentum.at(i)).Phi();
        float NewPt = (SmearedJet4Momentum.at(i)).Pt();
        float NewEta = (SmearedJet4Momentum.at(i)).Eta();
        float NewMass = (SmearedJet4Momentum.at(i)).M();
		
	NewMassVec.push_back(NewMass);

 }

 return NewMassVec;

}};


//For the Rochester corrections
auto MuonFourMomentum{[](

const floats& Muon_pt,
const floats& Muon_eta,
const floats& Muon_phi,
const floats& Muon_mass

){

  TLorentzVector Muon4Mo{};
  
  for(int i = 0; i < Muon_pt.size(); i++){
	
  	TLorentzVector vec{};
  	vec.SetPtEtaPhiM(Muon_pt.at(i), Muon_eta.at(i), Muon_phi.at(i), Muon_mass.at(i));
  	Muon4Mo += vec;

  }
  return Muon4Mo;

}};



auto RochCorrVec_Function{[&process, &year](

const ints& MuonCharge, 
const floats& MuonPt, 
const floats& MuonEta, 
const floats& MuonPhi, 
const ints& Muon_genPartIdx, 
const ints& Muon_nTrackerLayers

){

  floats CorrectionFactor = RochesterCorrections_testscript2(year, process, MuonCharge, MuonPt, MuonEta, MuonPhi, Muon_genPartIdx, Muon_nTrackerLayers);
  return CorrectionFactor;

}};



auto RochCorrVec_Function_data{[&process, &year](

const ints& MuonCharge,
const floats& MuonPt,
const floats& MuonEta,
const floats& MuonPhi,
const ints& DummyColumnInts,
const ints& Muon_nTrackerLayers

){

  floats CorrectionFactor = RochesterCorrections_testscript2(year, process, MuonCharge, MuonPt, MuonEta, MuonPhi, DummyColumnInts, Muon_nTrackerLayers);
  return CorrectionFactor;

}};




auto RochCorrMuon4Mo{[](const TLorentzVector& Muon4Mo, const floats& RochCorrVec){

  TLorentzVector NewVec{};

  double NewVecMass = Muon4Mo.M() * RochCorrVec.at(0);
  double NewVecPt = Muon4Mo.Pt() * RochCorrVec.at(0);
  double NewVecPhi = Muon4Mo.Phi() * RochCorrVec.at(0);
  double NewVecEta = Muon4Mo.Eta() * RochCorrVec.at(0);

  NewVec.SetPtEtaPhiM(NewVecPt, NewVecEta, NewVecPhi, NewVecMass);
  return NewVec;

}};



//For the normalisation factors
auto NormalisationFactorFunction{[&process, &year](){

  std::vector<std::string> ProcessStrings = {" ", "tZq", "ZPlusJets_M50_aMCatNLO", "ZPlusJets_M50_aMCatNLO_ext", "ZPlusJets_M50_Madgraph", "ZPlusJets_M50_Madgraph_ext",
				    "ZPlusJets_M10To50_aMCatNLO", "ZPlusJets_M10To50_aMCatNLO_ext", "ZPlusJets_M10To50_Madgraph", "ZPlusJets_M10To50_Madgraph_ext",
			            "ZPlusJets_PtBinned_0To50", "ZPlusJets_PtBinned_50To100", "ZPlusJets_PtBinned_50To100_ext", "ZPlusJets_PtBinned_100To250",
				    "ZPlusJets_PtBinned_100To250_ext1", "ZPlusJets_PtBinned_100To250_ext2", "ZPlusJets_PtBinned_100To250_ext5",
				    "ZPlusJets_PtBinned_250To400", "ZPlusJets_PtBinned_250To400_ext1", "ZPlusJets_PtBinned_250To400_ext2",
				    "ZPlusJets_PtBinned_250To400_ext5", "ZPlusJets_PtBinned_400To650", "ZPlusJets_PtBinned_400To650_ext1", 
				    "ZPlusJets_PtBinned_400To650_ext2", "ZPlusJets_PtBinned_650ToInf", "ZPlusJets_PtBinned_650ToInf_ext1", 
				    "ZPlusJets_PtBinned_650ToInf_ext2", "ttbar_2l2nu", "ttbar_madgraph", "ttbar_madgraph_ext", "ttbar_TTToHadronic", 
				    "ttbar_TTToSemileptonic", "ttbar_atMCaNLO", "ttbar_inc", "SingleTop_tchannel_top", "SingleTop_tchannel_top_ScaleUp", 
				    "SingleTop_tchannel_top_ScaleDown", "SingleTop_tchannel_antitop", "SingleTop_schannel", "ttbar_hdampUP", 
				    "ttbar_hdampUP_ext", "ttbar_hdampDOWN", "ttbar_hdampDOWN_ext", "SingleTop_tchannel_top_hdampUP", 
				    "SingleTop_tchannel_top_hdampDOWN", "ttbar_isr_UP", "ttbar_isr_DOWN", "ttbar_isr_DOWN_ext",
				    "ttbar_fsr_UP", "ttbar_fsr_UP_ext", "ttbar_fsr_DOWN", "ttbar_fsr_DOWN_ext", 
				    "SingleTop_tW", "SingleTop_tW_ScaleUp", "SingleTop_tW_ScaleDown", "SingleTop_tbarW", 
				    "SingleTop_tbarW_ScaleUp", "SingleTop_tbarW_ScaleDown", "SingleTop_tHq", "SingleTop_tZq_W_lept_Z_had", 
				    "SingleTop_tWZ_tWll", "VV_ZZTo2l2nu", "VV_ZZTo2l2nu_ext", "VV_ZZTo2l2Q", "VV_ZZTo4L", "VV_WW1nuqq", "VV_WZTo2l2Q", 
				    "VV_WZTo3lNu", "VV_WZTo1l2Nu2Q", "VV_WWTo2l2Nu", "VV_WWToLNuQQ", "VV_WWToLNuQQ_ext", "VV_WGToLNuG", "VV_ZGToLLG", 
				    "VVV_WWWTo4F", "VVV_WWZTo4F", "VVV_WZZ", "VVV_ZZZ", "WPlusJets", "WPlusJets_ext", "ttbarV_ttWJetsToLNu", 
				    "ttbarV_ttWJetsToLNu_ext", "ttbarV_ttZToLLNuNu", "ttbarV_ttWJetsToQQ", "ttbarV_ttZToLL", "ttbarV_ttZToLL_ext2", "ttbarV_ttZToLL_ext3", 
				    "ttbarV_ttgamma", "ttbarV_ttgamma_ext", "ttbarV_ttHTobb", "ttbarV_ttHToNonbb"};


  for(int i = 1; i < ProcessStrings.size(); i++){

	std::cout << "process = " << process << std::endl;
	std::cout << "ProcessStrings.at(i) = " << ProcessStrings.at(i) << std::endl;

	if(process == ProcessStrings.at(i)){return linereader(i, year);}
	else{continue;}

  }


}};









//Electron selection and reconstruction SFs
//2016
TFile* EGammaEff_inputfile_2016 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2016/egammaEffi_Tight_80X.txt_EGM2D.root", "READ");
TFile* EGammaEffSys_inputfile_2016 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2016/egammaEffi_Tight_80X.txt_EGM2D.root", "READ");
TFile* EGammaEffReco_inputfile_2016 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2016/egammaRecoEffi.txt_EGM2D.root", "READ");
TFile* EGammaEffRecoSys_inputfile_2016 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2016/egammaRecoEffi.txt_EGM2D.root", "READ");

//2017
TFile* EGammaEffReco_HigherPt_inputfile_2017 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2017/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root", "READ");
TFile* EGammaEffRecoSys_HigherPt_inputfile_2017 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2017/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root", "READ"); 
TFile* EGammaEffReco_LowPt_inputfile_2017 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2017/egammaEffi.txt_EGM2D_runBCDEF_passingRECO_lowEt.root", "READ");
TFile* EGammaEffRecoSys_LowPt_inputfile_2017 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2017/egammaEffi.txt_EGM2D_runBCDEF_passingRECO_lowEt.root", "READ");
TFile* EGammaEff_inputfile_2017 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2017/egammaEffi.txt_EGM2D_runBCDEF_passingTight94X.root", "READ");
TFile* EGammaEffSys_inputfile_2017 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2017/egammaEffi.txt_EGM2D_runBCDEF_passingTight94X.root", "READ");

//2018
TFile* EGammaEff_inputfile_2018 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2018/2018_ElectronTight.root", "READ"); //need to double check if this is the right file
TFile* EGammaEffSys_inputfile_2018 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2018/2018_ElectronTight.root", "READ");
TFile* EGammaEffReco_inputfile_2018 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2018/egammaEffi.txt_EGM2D_updatedAll.root", "READ");
TFile* EGammaEffRecoSys_inputfile_2018 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2018/egammaEffi.txt_EGM2D_updatedAll.root", "READ");

//Histograms
//2016
TH2* EGammaEff2016_histo = (TH2*)(EGammaEff_inputfile_2016->Get("EGamma_SF2D")->Clone());
EGammaEff2016_histo->SetDirectory(nullptr);
TH2* EGammaEffSys2016_histo = (TH2*)(EGammaEffSys_inputfile_2016->Get("EGamma_SF2D")->Clone());
EGammaEffSys2016_histo->SetDirectory(nullptr);
TH2* EGammaEffReco2016_histo = (TH2*)(EGammaEffReco_inputfile_2016->Get("EGamma_SF2D")->Clone());
EGammaEffReco2016_histo->SetDirectory(nullptr);
TH2* EGammaEffRecoSys2016_histo = (TH2*)(EGammaEffRecoSys_inputfile_2016->Get("EGamma_SF2D")->Clone());
EGammaEffRecoSys2016_histo->SetDirectory(nullptr);

//2017
TH2* EGammaEff2017_histo = (TH2*)(EGammaEff_inputfile_2017->Get("EGamma_SF2D")->Clone());
EGammaEff2017_histo->SetDirectory(nullptr);
TH2* EGammaEffSys2017_histo = (TH2*)(EGammaEffSys_inputfile_2017->Get("EGamma_SF2D")->Clone());
EGammaEffSys2017_histo->SetDirectory(nullptr);
TH2* EGammaEffReco_LowPt_2017_histo = (TH2*)(EGammaEffReco_LowPt_inputfile_2017->Get("EGamma_SF2D")->Clone());
EGammaEffReco_LowPt_2017_histo->SetDirectory(nullptr); 
TH2* EGammaEffRecoSys_LowPt_2017_histo = (TH2*)(EGammaEffRecoSys_LowPt_inputfile_2017->Get("EGamma_SF2D")->Clone());
EGammaEffRecoSys_LowPt_2017_histo->SetDirectory(nullptr);
TH2* EGammaEffReco_HigherPt_2017_histo = (TH2*)(EGammaEffReco_HigherPt_inputfile_2017->Get("EGamma_SF2D")->Clone());
EGammaEffReco_HigherPt_2017_histo->SetDirectory(nullptr);
TH2* EGammaEffRecoSys_HigherPt_2017_histo = (TH2*)(EGammaEffRecoSys_HigherPt_inputfile_2017->Get("EGamma_SF2D")->Clone());
EGammaEffRecoSys_HigherPt_2017_histo->SetDirectory(nullptr);

//2018
TH2* EGammaEff2018_histo = (TH2*)(EGammaEff_inputfile_2018->Get("EGamma_SF2D")->Clone());
EGammaEff2018_histo->SetDirectory(nullptr);
TH2* EGammaEffSys2018_histo = (TH2*)(EGammaEffSys_inputfile_2018->Get("EGamma_SF2D")->Clone());
EGammaEffSys2018_histo->SetDirectory(nullptr);
TH2* EGammaEffReco2018_histo = (TH2*)(EGammaEffReco_inputfile_2018->Get("EGamma_SF2D")->Clone());
EGammaEffReco2018_histo->SetDirectory(nullptr);
TH2* EGammaEffRecoSys2018_histo = (TH2*)(EGammaEffRecoSys_inputfile_2018->Get("EGamma_SF2D")->Clone());
EGammaEffRecoSys2018_histo->SetDirectory(nullptr);




//EGamma SF functions
auto EGammaFunction{[&EGammaEff2016_histo,     	     	     &EGammaEffSys2016_histo,
		     &EGammaEffReco2016_histo, 	     	     &EGammaEffRecoSys2016_histo,
		     &EGammaEff2017_histo,                   &EGammaEffSys2017_histo, 
		     &EGammaEffReco_LowPt_2017_histo,        &EGammaEffRecoSys_LowPt_2017_histo,
		     &EGammaEffReco_HigherPt_2017_histo,     &EGammaEffRecoSys_HigherPt_2017_histo,
		     &EGammaEff2018_histo,	             &EGammaEffSys2018_histo,
		     &EGammaEffReco2018_histo,	             &EGammaEffRecoSys2018_histo
		     ](const std::string& year, const std::string& type, const floats& pt, const floats& SuperClusterEta){


   std::cout << "inside EGammaFunction" << std::endl;

   floats OutputVector{};
   floats OutputVectorFinal{};

   for(int i = 0; i < pt.size(); i++){

  	if( abs(SuperClusterEta.at(i)) < 2.5 ){

		//2016
		int Bin_EGammaEff2016 = EGammaEff2016_histo->FindBin( SuperClusterEta.at(i), pt.at(i) );
		int Bin_EGammaEffSys2016 = EGammaEffSys2016_histo->FindBin( SuperClusterEta.at(i), pt.at(i) );
		int Bin_EGammaEffReco2016 = EGammaEffReco2016_histo->FindBin( SuperClusterEta.at(i), pt.at(i) );
                int Bin_EGammaEffRecoSys2016 = EGammaEffSys2016_histo->FindBin( SuperClusterEta.at(i), pt.at(i) );

		//2017
		int Bin_EGammaEff2017 = EGammaEff2017_histo->FindBin(SuperClusterEta.at(i), pt.at(i));
		int Bin_EGammaEffSys2017 = EGammaEffSys2017_histo->FindBin( SuperClusterEta.at(i), pt.at(i) );
		int Bin_EGammaEffReco_LowPt_2017 = EGammaEffReco_LowPt_2017_histo->FindBin( SuperClusterEta.at(i), pt.at(i) );
		int Bin_EGammaEffRecoSys_LowPt_2017 = EGammaEffRecoSys_LowPt_2017_histo->FindBin( SuperClusterEta.at(i), pt.at(i) );
		int Bin_EGammaEffReco_HigherPt_2017 = EGammaEffReco_HigherPt_2017_histo->FindBin( SuperClusterEta.at(i), pt.at(i) );
                int Bin_EGammaEffRecoSys_HigherPt_2017 = EGammaEffRecoSys_HigherPt_2017_histo->FindBin( SuperClusterEta.at(i), pt.at(i) );


		//2018
		int Bin_EGammaEff2018 = EGammaEff2018_histo->FindBin( SuperClusterEta.at(i), pt.at(i) );
                int Bin_EGammaEffSys2018 = EGammaEffSys2018_histo->FindBin( SuperClusterEta.at(i), pt.at(i) );
                int Bin_EGammaEffReco2018 = EGammaEffReco2018_histo->FindBin( SuperClusterEta.at(i), pt.at(i) );
                int Bin_EGammaEffRecoSys2018 = EGammaEffSys2018_histo->FindBin( SuperClusterEta.at(i), pt.at(i) );

		float EGammaSF;

		if(year == "2016"){
			if(type == "EGammaEffSys"){EGammaSF = EGammaEffSys2016_histo->GetBinError(Bin_EGammaEffSys2016);}
			else if(type == "EGammaEffRecoSys"){EGammaSF = EGammaEffRecoSys2016_histo->GetBinError(Bin_EGammaEffRecoSys2016);}
			else if(type == "EGammaEff"){EGammaSF = EGammaEff2016_histo->GetBinContent(Bin_EGammaEff2016);}
			else if(type == "EGammaEffReco"){EGammaSF = EGammaEffReco2016_histo->GetBinContent(Bin_EGammaEffReco2016);}
			else{std::cout << "Choose a type out of EGammaEffSys, EGammaEffRecoSys, EGammaEff or EGammaEffReco for 2016" << std::endl;}
		}
		else if(year == "2017"){
			if(type == "EGammaEffSys"){EGammaSF = EGammaEffSys2017_histo->GetBinError(Bin_EGammaEffSys2017);}
                        else if(type == "EGammaEffRecoSys" && pt.at(i) <= 20){EGammaSF = EGammaEffRecoSys_LowPt_2017_histo->GetBinError(Bin_EGammaEffRecoSys_LowPt_2017);}
			else if(type == "EGammaEffRecoSys" && pt.at(i) > 20){EGammaSF = EGammaEffRecoSys_HigherPt_2017_histo->GetBinError(Bin_EGammaEffRecoSys_HigherPt_2017);}
                        else if(type == "EGammaEff"){std::cout << "here 2" << std::endl; EGammaSF = EGammaEff2017_histo->GetBinContent(Bin_EGammaEff2017);}
			else if(type == "EGammaEffReco" && pt.at(i) <= 20){EGammaSF = EGammaEffReco_LowPt_2017_histo->GetBinContent(Bin_EGammaEffReco_LowPt_2017);}
                        else if(type == "EGammaEffReco" && pt.at(i) > 20){EGammaSF = EGammaEffReco_HigherPt_2017_histo->GetBinContent(Bin_EGammaEffReco_HigherPt_2017);}
                        else{std::cout << "Choose a type out of EGammaEffSys, EGammaEffRecoSys, EGammaEff or EGammaEffReco for 2017" << std::endl;}

		}
		else if(year == "2018"){
			if(type == "EGammaEffSys"){EGammaSF = EGammaEffSys2018_histo->GetBinError(Bin_EGammaEffSys2018);}
                        else if(type == "EGammaEffRecoSys"){EGammaSF = EGammaEffRecoSys2018_histo->GetBinError(Bin_EGammaEffRecoSys2018);}
                        else if(type == "EGammaEff"){EGammaSF = EGammaEff2018_histo->GetBinContent(Bin_EGammaEff2018);}
                        else if(type == "EGammaEffReco"){EGammaSF = EGammaEffReco2018_histo->GetBinContent(Bin_EGammaEffReco2018);}
                        else{std::cout << "Choose a type out of EGammaEffSys, EGammaEffRecoSys, EGammaEff or EGammaEffReco for 2018" << std::endl;}

		}
	

		OutputVector.push_back(EGammaSF);
		
  	}
  	else{OutputVector.push_back(1.0);}



  } //end of for loop


  for(int i = 0; i < OutputVector.size(); i++){

	if(OutputVector.at(i) == 0){OutputVectorFinal.push_back(1.0);}
	else{OutputVectorFinal.push_back( OutputVector.at(i) );}

  }

  return OutputVectorFinal.at(0);

}};


EGammaEff_inputfile_2016->Close();
EGammaEffSys_inputfile_2016->Close();
EGammaEffReco_inputfile_2016->Close();
EGammaEffRecoSys_inputfile_2016->Close(); 
EGammaEffReco_HigherPt_inputfile_2017->Close();
EGammaEffRecoSys_HigherPt_inputfile_2017->Close();
EGammaEffReco_LowPt_inputfile_2017->Close();
EGammaEffRecoSys_LowPt_inputfile_2017->Close();
EGammaEff_inputfile_2017->Close();
EGammaEffSys_inputfile_2017->Close();
EGammaEff_inputfile_2018->Close();
EGammaEffSys_inputfile_2018->Close();
EGammaEffReco_inputfile_2018->Close();
EGammaEffRecoSys_inputfile_2018->Close();



auto EGammaSF_egammaEff{[&year, &EGammaFunction](const floats& Electron_pt_Selection, const floats& SuperClusterEta){

  return EGammaFunction(year, "EGammaEff", Electron_pt_Selection, SuperClusterEta);

}};





auto EGammaSF_egammaEffReco{[&year, &EGammaFunction](const floats& Electron_pt_Selection, const floats& SuperClusterEta){

  return EGammaFunction(year, "EGammaEffReco", Electron_pt_Selection, SuperClusterEta);

}};





auto EGammaSF_egammaEff_Sys{[&year, &EGammaFunction](const floats& Electron_pt_Selection, const floats& SuperClusterEta){

  return EGammaFunction(year, "EGammaEffSys", Electron_pt_Selection, SuperClusterEta);

}};





auto EGammaSF_egammaEffReco_Sys{[&year, &EGammaFunction](const floats& Electron_pt_Selection, const floats& SuperClusterEta){

  return EGammaFunction(year, "EGammaEffRecoSys", Electron_pt_Selection, SuperClusterEta);

}};








//Lambda functions for lepton efficiencies

auto MuonSF{[&year](const std::string& type, const std::string& year, const std::string& UpOrDown, const floats& pt, const floats& eta){

  std::cout << "inside MuonSF" << std::endl;

  floats AbsEta = abs(eta);

  float lumiRunBCDEF = 19713.888;
  float lumiRunGH = 16146.178;

  TFile* inputfile_RunsBCDEF;
  TFile* inputfile_RunsGH;
  TH2* histo_RunsBCDEF;
  TH2* histo_RunsGH;
  TFile* inputfile_RunsABCD;
  TH2* histo_RunsABCD;


  if(year == "2016"){

  	if(type == "ID" || type == "ID sys"){

		//Muon ID file (runs BCDEF)
		inputfile_RunsBCDEF = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2016/MuonID_EfficienciesAndSF_BCDEF.root", "READ");
		histo_RunsBCDEF = (TH2*)(inputfile_RunsBCDEF->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio")->Clone());
		histo_RunsBCDEF->SetDirectory(nullptr);


                //Muon ID file (runs GH)
                inputfile_RunsGH = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2016/MuonID_EfficienciesAndSF_GH.root", "READ");
                histo_RunsGH = (TH2*)(inputfile_RunsGH->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio")->Clone());
		histo_RunsGH->SetDirectory(nullptr);

        } 
	else if(type == "Iso" || type == "Iso sys"){
        
                //Muon ISO file (runs BCDEF)
                inputfile_RunsBCDEF = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2016/MuonISO_EfficienciesAndSF_BCDEF.root", "READ");
                histo_RunsBCDEF = (TH2*)(inputfile_RunsBCDEF->Get("TightISO_TightID_pt_eta/pt_abseta_ratio")->Clone());
		histo_RunsBCDEF->SetDirectory(nullptr);                


                //Muon ISO file (runs GH)
                inputfile_RunsGH = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2016/MuonISO_EfficienciesAndSF_GH.root", "READ");
                histo_RunsGH = (TH2*)(inputfile_RunsGH->Get("TightISO_TightID_pt_eta/pt_abseta_ratio")->Clone());
		histo_RunsGH->SetDirectory(nullptr);

        }
	else{std::cout << "Please choose either ID or ISO for the type" << std::endl;}	

  }
  else if(year == "2017"){

  	if(type == "ID"){

  		//Muon ID file
  		inputfile_RunsBCDEF = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ID.root", "READ");
  		histo_RunsBCDEF = (TH2*)(inputfile_RunsBCDEF->Get("NUM_TightID_DEN_genTracks_pt_abseta")->Clone());  
		histo_RunsBCDEF->SetDirectory(nullptr);

  	}
 	 else if(type == "ID sys"){

  		//Muon ID sys file
  		inputfile_RunsBCDEF = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ID_syst.root", "READ");
  		histo_RunsBCDEF = (TH2*)(inputfile_RunsBCDEF->Get("NUM_TightID_DEN_genTracks_pt_abseta")->Clone());
		histo_RunsBCDEF->SetDirectory(nullptr);

  	}
 	 else if(type == "ID sys (stat)"){

        	//Muon ID sys (stat)
        	inputfile_RunsBCDEF = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ID_syst.root", "READ");
  		histo_RunsBCDEF = (TH2*)(inputfile_RunsBCDEF->Get("NUM_TightID_DEN_genTracks_pt_abseta_stat")->Clone());
		histo_RunsBCDEF->SetDirectory(nullptr);

  	}
  	else if(type == "ID sys (syst)"){
 
        	//Muon ID sys (syst)
        	inputfile_RunsBCDEF = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ID_syst.root", "READ");
  		histo_RunsBCDEF = (TH2*)(inputfile_RunsBCDEF->Get("NUM_TightID_DEN_genTracks_pt_abseta_syst")->Clone());
		histo_RunsBCDEF->SetDirectory(nullptr);

  	}
  	else if(type == "Iso"){

  		//Muon Iso file
  		inputfile_RunsBCDEF = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ISO.root", "READ");
  		histo_RunsBCDEF = (TH2*)(inputfile_RunsBCDEF->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta")->Clone());
		histo_RunsBCDEF->SetDirectory(nullptr); 

  	}
  	else if(type == "Iso sys"){
  
  		//Muon Iso sys file
  		inputfile_RunsBCDEF = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ISO_syst.root", "READ");
  		histo_RunsBCDEF = (TH2*)(inputfile_RunsBCDEF->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta")->Clone());
		histo_RunsBCDEF->SetDirectory(nullptr);

  	}
  	else if(type == "Iso sys (stat)"){

   		inputfile_RunsBCDEF = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ISO_syst.root", "READ");
  		histo_RunsBCDEF = (TH2*)(inputfile_RunsBCDEF->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta_stat")->Clone());
		histo_RunsBCDEF->SetDirectory(nullptr);  

  	}
  	else if(type == "Iso sys (syst)"){

  		inputfile_RunsBCDEF = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ISO_syst.root", "READ");
		histo_RunsBCDEF = (TH2*)(inputfile_RunsBCDEF->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta_syst")->Clone()); 
		histo_RunsBCDEF->SetDirectory(nullptr);  

  	}
  	else{std::cout << "Choose a function input out of: ID, ID sys, ID sys (stat), ID sys (syst), Iso sys, Iso sys (stat), Iso sys (syst)." << std::endl;}

  }
  else if(year == "2018"){

	if(type == "ID" || type == "ID sys"){

                //Muon ID file (runs ABCD)
		inputfile_RunsABCD = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2018/RunABCD_SF_ID.root", "READ"); //need to double check if root file is correct
                histo_RunsABCD = (TH2*)(inputfile_RunsABCD->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio")->Clone());
		histo_RunsABCD->SetDirectory(nullptr);


        }
        else if(type == "Iso" || type == "Iso sys"){

                //Muon ISO file (runs ABCD)
                inputfile_RunsABCD = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2018/RunABCD_SF_ISO.root", "READ"); //need to double check if root file is correct
                histo_RunsABCD = (TH2*)(inputfile_RunsABCD->Get("TightISO_TightID_pt_eta/pt_abseta_ratio")->Clone());
		histo_RunsABCD->SetDirectory(nullptr);


        }
        else{std::cout << "Please choose either ID or ISO for the type" << std::endl;}	

  }
  else{std::cout << "Select a year out of 2016, 2017 or 2018." << std::endl;}

 

  floats MuonSFOutput{};
 
  for(int i = 0; i < pt.size(); i++){

	if(pt.at(i) >= 20 && pt.at(i) <= 120 && AbsEta.at(i) <= MaxTrackerEta){

  		int PtBin_RunsBCDEF = histo_RunsBCDEF->GetXaxis()->FindBin(pt.at(i));
        	int AbsEtaBin_RunsBCDEF = histo_RunsBCDEF->GetYaxis()->FindBin(AbsEta.at(i));

		//for 2016 and 2017
        	float MuonSF_RunsBCDEF = histo_RunsBCDEF->GetBinContent(PtBin_RunsBCDEF, AbsEtaBin_RunsBCDEF);
		float Error_RunsBCDEF = histo_RunsBCDEF->GetBinError(PtBin_RunsBCDEF, AbsEtaBin_RunsBCDEF);
	
		//for 2018
		float MuonSF_RunsABCD = histo_RunsBCDEF->GetBinContent(PtBin_RunsBCDEF, AbsEtaBin_RunsBCDEF);
                float Error_RunsABCD = histo_RunsBCDEF->GetBinError(PtBin_RunsBCDEF, AbsEtaBin_RunsBCDEF);

		//for 2016
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
				else{std::cout << "Select an up or down uncertainty" << std::endl;}

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
				else{std::cout << "Select an up or down uncertainty" << std::endl;}


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
		else if(year == "2018"){


                        if(type == "ID sys (stat)" ||
                           type == "ID sys (syst)" ||
                           type == "Iso sys (stat)" ||
                           type == "Iso sys (syst)"){

                                MuonSFOutput.push_back(Error_RunsABCD);
                        }
                        else{MuonSFOutput.push_back(MuonSF_RunsABCD);}

                }
		else{std::cout << "Code only for 2016, 2017 or 2018." << std::endl;}


	}
	else{std::cout << "inside else statement for pushing back muon SF output" << std::endl; auto One = 1.0; MuonSFOutput.push_back(One);}



  }


  inputfile_RunsBCDEF->Close();
  inputfile_RunsGH->Close();
  inputfile_RunsABCD->Close();

  return MuonSFOutput.at(0); 

}};





auto MuonSFTest_ID{[&MuonSF, &year](const floats& pt, const floats& eta){

  return MuonSF("ID", year, " ", pt, eta);
  
}};



auto MuonSFTest_Iso{[&MuonSF, &year](const floats& pt, const floats& eta){

  return MuonSF("Iso", year, " ", pt, eta);
    
}};




auto MuonSFTest_ID_sys_syst{[&MuonSF, &year](const floats& pt, const floats& eta){
  
  if(year == "2016"){
        return MuonSF("ID sys", year, "Up", pt, eta);
  }
  else if(year == "2017"){
        return MuonSF("ID sys (syst)", year, " ", pt, eta);
  }  
  else{std::cout << "Need to add 2018" << std::endl;}

}};




auto MuonSFTest_ID_sys_stat{[&MuonSF, &year](const floats& pt, const floats& eta){

  if(year == "2016"){
	return MuonSF("ID sys", year, "Down", pt, eta);
  }
  else if(year == "2017"){
  	return MuonSF("ID sys (stat)", year, " ", pt, eta);
  }
  else{std::cout << "Need to add 2018" << std::endl;}  

}};




auto MuonSFTest_Iso_sys_syst{[&MuonSF, &year](const floats& pt, const floats& eta){

  if(year == "2016"){
        return MuonSF("Iso sys", year, "Up", pt, eta);
  }
  else if(year == "2017"){
        return MuonSF("Iso sys (syst)", year, " ", pt, eta);
  }
  else{std::cout << "Need to add 2018" << std::endl;}

}};




auto MuonSFTest_Iso_sys_stat{[&MuonSF, &year](const floats& pt, const floats& eta){

  if(year == "2016"){
        return MuonSF("Iso sys", year, "Down", pt, eta);
  }
  else if(year == "2017"){
        return MuonSF("Iso sys (stat)", year, " ", pt, eta);
  }
  else{std::cout << "Need to add 2018" << std::endl;}

}};









//dummy column
auto DummyColumnFunction{[](const floats& pts){

  return pts;

}};

auto DummyColumnFunctionInts{[](const ints& charges){

  return charges;

}};


//PS weight lambda function
auto PSWeight{[&year, &process](floats& PSWeight, const floats& pts){

  floats Ones(4, 1.0);

  if(year == "2017" || year == "2018"){

  	if(process == "tZq" ||
	   process == "SingleTop_tbarW" ||
	   process == "SingleTop_schannel" ||
	   process == "SingleTop_tchannel_top" ||
	   process == "SingleTop_tchannel_tbar" ||
	   process == "ttbarV_ttgamma" ||
	   process == "ttbar_TTToHadronic" ||
	   process == "ttbar_TTToSemileptonic"){return PSWeight;}
	else{return Ones;}

 
 }
 else{return Ones;}


}};



//Functions for trigger SFs
///Lambda function for the MET triggers


auto MET_triggers_function{[&year](

const bool& HLT_MET200,
const bool& HLT_MET250,
const bool& HLT_PFMET120_PFMHT120_IDTight,
const bool& HLT_PFMET170_HBHECleaned,
const bool& HLT_PFHT300_PFMET100,
//const bool& HLT_MET105_IsoTrk50, //(not in run B)
//const bool& HLT_MET120_IsoTrk50, //(not in run B)
//const bool& HLT_HT430_DisplacedDijet40_DisplacedTrack, //(not in run B)
//const bool& HLT_HT650_DisplacedDijet60_Inclusive, //(not in run B)
//const bool& HLT_HT750_DisplacedDijet80_Inclusive, //(not in run B)
//const bool& HLT_PFMET120_PFMHT120_IDTight_HFCleaned, //(not in any)
//const bool& HLT_PFMET120_PFMHT120_IDTight_L1ETMnoHF, //(not in any)
//const bool& HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned, //(not in any)
//const bool& HLT_PFMET120_PFMHT120_IDTight_PFHT60, //(not in run B)
const bool& HLT_PFMET130_PFMHT130_IDTight,
const bool& HLT_PFMET140_PFMHT140_IDTight,
//const bool& HLT_PFMET200_HBHE_BeamHaloCleaned, //(not in run B)
//const bool& HLT_PFMET250_HBHECleaned, //(not in run B)
//const bool& HLT_PFMET300_HBHECleaned, //(not in run B)
//const bool& HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned, (not in any)
//const bool& HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_L1ETMnoHF, (not in any)
//const bool& HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60, //(not in run B)
const bool& HLT_PFMETNoMu120_PFMHTNoMu120_IDTight,
//const bool& HLT_PFMETNoMu130_PFMHTNoMu130_IDTight, //(not in run B)
//const bool& HLT_PFMETNoMu140_PFMHTNoMu140_IDTight, //(not in run B)
const bool& HLT_PFHT1050,
const bool& HLT_PFHT180,
const bool& HLT_PFHT500_PFMET100_PFMHT100_IDTight,
const bool& HLT_PFHT500_PFMET110_PFMHT110_IDTight,
const bool& HLT_PFHT700_PFMET85_PFMHT85_IDTight,
const bool& HLT_PFHT700_PFMET95_PFMHT95_IDTight,
const bool& HLT_PFHT800_PFMET75_PFMHT75_IDTight,
const bool& HLT_PFHT800_PFMET85_PFMHT85_IDTight)->bool{

  if(year == "2016"){
	
	return HLT_MET200 > 0 ||
        HLT_MET250 > 0 ||
        HLT_PFMET120_PFMHT120_IDTight > 0 ||
        HLT_PFMET170_HBHECleaned > 0 ||
        HLT_PFHT300_PFMET100 > 0;

  }
  else if(year == "2017" || year == "2018"){

  return 

  //HLT_MET105_IsoTrk50 > 0 || //(not in run B)
  //HLT_MET120_IsoTrk50 > 0 || //(not in run B)
  //HLT_HT430_DisplacedDijet40_DisplacedTrack > 0 || //(not in run B)
  //HLT_HT650_DisplacedDijet60_Inclusive > 0 || //(not in run B)
  //HLT_HT750_DisplacedDijet80_Inclusive > 0 || //(not in run B)
  //HLT_PFMET120_PFMHT120_IDTight_HFCleaned > 0 || //(not in any)
  //HLT_PFMET120_PFMHT120_IDTight_L1ETMnoHF > 0 || //(not in any)
  //HLT_PFMET120_PFMHT120_IDTight_PFHT60_HFCleaned > 0 || //(not in any)
  //HLT_PFMET120_PFMHT120_IDTight_PFHT60 > 0 || //(not in run B)
  HLT_PFMET120_PFMHT120_IDTight > 0 ||
  HLT_PFMET130_PFMHT130_IDTight > 0 ||
  HLT_PFMET140_PFMHT140_IDTight > 0 ||
  //HLT_PFMET200_HBHE_BeamHaloCleaned > 0 || //(not in run B)
  //HLT_PFMET250_HBHECleaned > 0 || //(not in run B)
  //HLT_PFMET300_HBHECleaned > 0 || //(not in run B)
  //HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_HFCleaned > 0 || (not in any)
  //HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_L1ETMnoHF > 0 || (not in any)
  //HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 > 0 || //(not in run B)
  HLT_PFMETNoMu120_PFMHTNoMu120_IDTight > 0 ||
  //HLT_PFMETNoMu130_PFMHTNoMu130_IDTight > 0 || //(not in run B)
  //HLT_PFMETNoMu140_PFMHTNoMu140_IDTight > 0 || //(not in run B)
  HLT_PFHT1050 > 0 ||
  HLT_PFHT180 > 0 ||
  HLT_PFHT500_PFMET100_PFMHT100_IDTight > 0 ||
  HLT_PFHT500_PFMET110_PFMHT110_IDTight > 0 ||
  HLT_PFHT700_PFMET85_PFMHT85_IDTight > 0 ||
  HLT_PFHT700_PFMET95_PFMHT95_IDTight > 0 ||
  HLT_PFHT800_PFMET75_PFMHT75_IDTight > 0 ||
  HLT_PFHT800_PFMET85_PFMHT85_IDTight > 0;

  }
  else{std::cout << "Choose a year out of 2016, 2017 or 2018" << std::endl;}

}};






//Events that pass the selection criteria and the lepton triggers

auto ee_selection_LL_Trig_function{[&year](

const bool& HLT_Ele32_WPTight_Gsf_L1DoubleEG, 
const bool& HLT_Ele35_WPTight_Gsf, 
const bool& HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL, 
const bool& HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, 
const bool& HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, 
const bool& HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8, 
const bool& HLT_IsoMu27, 
const bool& HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, 
const bool& HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, 
const bool& HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, 
const bool& HLT_Ele25_eta2p1_WPTight_Gsf, 
const bool& HLT_Ele27_WPTight_Gsf, 
const bool& HLT_Ele32_eta2p1_WPTight_Gsf,
const bool& HLT_IsoMu24,
const bool& HLT_IsoMu24_eta2p1,
const bool& HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,
const bool& HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,
const bool& HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,
const bool& HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8
)-> bool{

  if(year == "2016"){

	return //single or double electron and not any of the others
	
	(HLT_Ele25_eta2p1_WPTight_Gsf > 0 || //single electron 
	 HLT_Ele27_WPTight_Gsf > 0 || //single electron
	 HLT_Ele32_eta2p1_WPTight_Gsf > 0 || //single electron
	 HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ > 0) //double electron 
	 && 
	(HLT_IsoMu24 <= 0 || //single muon
	 HLT_IsoMu24_eta2p1 <= 0 || //single muon
	 HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ <= 0 || //double muon
	 HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 <= 0 || //double muon
	 HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || //muon electron
	 HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || //muon electron
	 HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || //muon electron
	 HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL <= 0 || //muon electron
	 HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL <= 0 || //muon electron
	 HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL <= 0); //muon electron


  }
  else if(year == "2017"){

  return 

  (HLT_Ele32_WPTight_Gsf_L1DoubleEG > 0 || //single electron
   HLT_Ele35_WPTight_Gsf > 0 || //single electron
   HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL > 0 || //double electron
   HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ > 0) && // double electron
  (//HLT_IsoMu24 <= 0 || //single muon
   HLT_IsoMu27 <= 0 || //single muon
   HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ <= 0 || //double muon
   HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 <= 0 || //double muon
   //HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 <= 0 || //double muon (not in MET run B)
   //HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL <= 0 || //muon+electron (not in MET run B)
   HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || //muon+electron
   //HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL <= 0 || //muon+electron (not in MET run B)
   HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || //muon+electron
   //HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL <= 0 || //muon+electron (not in MET run B)
   HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0); //muon+electron

  }
  else if(year == "2018"){

	return (HLT_Ele32_WPTight_Gsf_L1DoubleEG > 0 || //single electron
		HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL > 0 || //single electron
		HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ > 0) //double electron

		&&

		(HLT_IsoMu24 <= 0 || //single muon
		 HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 <= 0 || //double muon
		 HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 <= 0 || //double muon
		 HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || //muon+electron
		 HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || //muon+electron
		 HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0); //muon+electron


  }
  else{std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl;}

}};


auto mumu_selection_LL_Trig_function{[&year](

const bool& HLT_Ele32_WPTight_Gsf_L1DoubleEG, 
const bool& HLT_Ele35_WPTight_Gsf, 
const bool& HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL, 
const bool& HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, 
const bool& HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, 
const bool& HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8, 
//const bool& HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8, (not in MET Run B)
const bool& HLT_IsoMu27, 
//const bool& HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, (not in MET run B)
const bool& HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, 
//const bool& HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, (not in MET run B)
const bool& HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, 
//const bool& HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, (not in MET run B)
const bool& HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,
const bool& HLT_Ele25_eta2p1_WPTight_Gsf,
const bool& HLT_Ele27_WPTight_Gsf, 
const bool& HLT_Ele32_eta2p1_WPTight_Gsf,
const bool& HLT_IsoMu24,
const bool& HLT_IsoMu24_eta2p1,
const bool& HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,
const bool& HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,
const bool& HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,
const bool& HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8)->bool{


  if(year == "2016"){

        return //single or double muon and not any of the others

	(HLT_IsoMu24 > 0 || //single muon
         HLT_IsoMu24_eta2p1 > 0 || //single muon
         HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ > 0 || //double muon
         HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 > 0) //double muon

	&&

        (HLT_Ele25_eta2p1_WPTight_Gsf <= 0 || //single electron 
         HLT_Ele27_WPTight_Gsf <= 0 || //single electron
         HLT_Ele32_eta2p1_WPTight_Gsf <= 0 || //single electron
         HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || //double electron 
         HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || //muon electron
         HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || //muon electron
         HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || //muon electron
         HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL <= 0 || //muon electron
         HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL <= 0 || //muon electron
         HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL <= 0); //muon electron


  }
  else if(year == "2017"){

  return
 
  (HLT_IsoMu27 > 0  || //single muon
  HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ > 0 || //double muon
  HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 > 0) //|| //double muon
  //HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 > 0) && //double muon (not in MET Run B)
  && (//HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL <= 0 || //muon+electron (not in MET Run B)
  HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || //muon+electron
  //HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL <= 0 || //muon+electron (not in MET Run B)
  HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || //muon+electron
  //HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL <= 0 || //muon+electron (not in MET Run B)
  HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || //muon+electron
  HLT_Ele32_WPTight_Gsf_L1DoubleEG <= 0 || //single electron
  HLT_Ele35_WPTight_Gsf <= 0 || //single electron
  HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL <= 0 || //double electron
  HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ <= 0);  //double electron

  }
  else if(year == "2018"){

	return //single or double muon and not any of the others
	
	(HLT_IsoMu24 > 0 || //single muon
         HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 > 0 || //double muon
         HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 > 0) //double muon

	&&

	(HLT_Ele32_WPTight_Gsf_L1DoubleEG <= 0 || //single electron
         HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL <= 0 || //single electron
         HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || //double electron
	 HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || //muon+electron
         HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || //muon+electron
         HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0); //muon+electron

  }
  else{std::cout << "Choose 2016, 2017 or 2018" << std::endl;}
  

}};


auto emu_selection_LL_Trig_function{[&year](

const bool& HLT_Ele32_WPTight_Gsf_L1DoubleEG, 
const bool& HLT_Ele35_WPTight_Gsf, 
const bool& HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL, 
const bool& HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, 
const bool& HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, 
const bool& HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8, 
//const bool& HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8, (not in MET Run B)
const bool& HLT_IsoMu27, 
//const bool& HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, (not in MET Run B)
const bool& HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, 
//const bool& HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, (not in MET Run B) 
const bool& HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, 
//const bool& HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, (not in MET Run B)
const bool& HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,
const bool& HLT_Ele25_eta2p1_WPTight_Gsf,
const bool& HLT_Ele27_WPTight_Gsf, 
const bool& HLT_Ele32_eta2p1_WPTight_Gsf,
const bool& HLT_IsoMu24,
const bool& HLT_IsoMu24_eta2p1,
const bool& HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,
const bool& HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,
const bool& HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,
const bool& HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8)->bool{


  if(year == "2016"){

	return //single lepton or muon+electron and not any of the others

        (HLT_IsoMu24 > 0 || //single muon
         HLT_IsoMu24_eta2p1 > 0 || //single muon
	 HLT_Ele25_eta2p1_WPTight_Gsf > 0 || //single electron 
         HLT_Ele27_WPTight_Gsf > 0 || //single electron
         HLT_Ele32_eta2p1_WPTight_Gsf > 0 || //single electron
	 HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ > 0 || //muon electron
         HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ > 0 || //muon electron
         HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ > 0 || //muon electron
         HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL > 0 || //muon electron
         HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL > 0 || //muon electron
         HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL > 0) //muon electron

         &&

         (HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ <= 0 || //double muon
          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 <= 0 || //double muon
          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ <= 0); //double electron 



  }
  else if(year == "2017"){

  	return 

  	(HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL <= 0 || // double electron
  	HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || //double electron
  	HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ <= 0 || //double muon
  	HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 <= 0 //|| //double muon
  	//HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 <= 0) && //double muon (not in MET Run B)
  	) &&

  	(HLT_Ele32_WPTight_Gsf_L1DoubleEG > 0 || //single electron
  	HLT_Ele35_WPTight_Gsf > 0 || //single electron
  	HLT_IsoMu27 > 0 || //single muon
  	//HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL > 0 || //muon+electron (not in MET Run B)
  	HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ > 0 || //muon+electron
  	//HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL > 0 || //muon+electron (not in MET Run B)
  	HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ > 0 || //muon+electron 
  	//HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL > 0 || //muon+electron (not in MET Run B)
  	HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ > 0); //muon+electron

  }
  else if(year == "2018"){

	 return //single lepton or muon+electron and not any of the others

        (HLT_IsoMu24 > 0 || //single muon
	 HLT_Ele32_WPTight_Gsf_L1DoubleEG > 0 || //single electron
         HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL > 0 || //single electron
         HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ > 0 || //muon+electron
         HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ > 0 || //muon+electron
         HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ > 0) //muon+electron

	 &&

         (HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 > 0 || //double muon
          HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 > 0 || //double muon
          HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ <= 0); //double electron
	
  }
  else{std::cout << "Input 2016, 2017 or 2018" << std::endl;}

}};








//MET triggers
std::vector<std::string> MET_triggers;

if(year == "2016"){ 

	MET_triggers = {"HLT_MET200",
		        "HLT_MET250",
		        "HLT_PFMET120_PFMHT120_IDTight",
			"HLT_PFMET170_HBHECleaned",
			"HLT_PFHT300_PFMET100",
			"DummyBool",
			"DummyBool",
			"DummyBool",
			"DummyBool",
                        "DummyBool",
                        "DummyBool",
                        "DummyBool",
			"DummyBool",
                        "DummyBool",
                        "DummyBool",
                        "DummyBool"};


}
else if(year == "2017" || year == "2018"){

	MET_triggers = {"DummyBool",
			"DummyBool",
			"DummyBool",
                        "DummyBool",
			"DummyBool",
			"HLT_PFMET130_PFMHT130_IDTight",
			"HLT_PFMET140_PFMHT140_IDTight",
			"HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
			"HLT_PFHT1050",
			"HLT_PFHT180",
			"HLT_PFHT500_PFMET100_PFMHT100_IDTight",
			"HLT_PFHT500_PFMET110_PFMHT110_IDTight",
			"HLT_PFHT700_PFMET85_PFMHT85_IDTight",
			"HLT_PFHT700_PFMET95_PFMHT95_IDTight",
			"HLT_PFHT800_PFMET75_PFMHT75_IDTight",
			"HLT_PFHT800_PFMET85_PFMHT85_IDTight"};

}
else{std::cout << "Code is only for the years 2016, 2017 and 2018" << std::endl;}


//dummy lambda function
auto DummyBool{[](const bool& dummyinput){

  return dummyinput > 0;

}};



std::vector<std::string> leptontriggers_strings;

if(year == "2016"){

	leptontriggers_strings = {

	"DummyBool", //single electron (for 2018 only)
	"DummyBool", //single electron (for 2017 only)
	"DummyBool", //double electron (for 2017 and 2018)
	"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", //double electron (for all years)
	"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", //double muon (for 2016 and 2017)
	"DummyBool", //double muon (for 2017 and 2018)
	"DummyBool", //single muon (for 2017 only)
//	"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", //muon+electron (for all years) (branch not present in MET 2016)
//	"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", //muon+electron (for all years) (branch not present in MET 2016)
//	"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", //muon+electron (for all years) (branch not present in MET 2016)
	"DummyBool",
	"DummyBool",
	"DummyBool",
	"HLT_Ele25_eta2p1_WPTight_Gsf", //single electron (for 2016 only)
	"HLT_Ele27_WPTight_Gsf", //single electron (for 2016 only)
	"HLT_Ele32_eta2p1_WPTight_Gsf", //single electron (for 2016 only)
	"HLT_IsoMu24", //single muon (all years)
//	"HLT_IsoMu24_eta2p1", //single muon (for 2016 only) (branch not present in MET 2016)
	"DummyBool",
	"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", //muon+electron (for 2016 only)
//	"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", //muon+electron (for 2016 only) (branch not present in MET 2016)
	"DummyBool",
	"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", //muon+electron (for 2016 only)
	"DummyBool" //double muon (for 2018 only)

};

}
else if(year == "2017"){

	leptontriggers_strings = {

        "DummyBool", //single electron (for 2018 only)
        "HLT_Ele35_WPTight_Gsf", //single electron (for 2017 only)
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", //double electron (for 2017 and 2018)
        "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", //double electron (for all years)
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", //double muon (for 2016 and 2017)
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", //double muon (for 2017 and 2018)
        "HLT_IsoMu27", //single muon (for 2017 only)
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", //muon+electron (for all years)
        "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", //muon+electron (for all years)
        "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", //muon+electron (for all years)
        "DummyBool", //single electron (for 2016 only)
        "DummyBool", //single electron (for 2016 only)
        "DummyBool", //single electron (for 2016 only)
        "HLT_IsoMu24", //single muon (all years)
        "DummyBool", //single muon (for 2016 only)
        "DummyBool", //muon+electron (for 2016 only)
        "DummyBool", //muon+electron (for 2016 only)
        "DummyBool", //muon+electron (for 2016 only)
        "DummyBool" //double muon (for 2018 only)	

};

}
else if(year == "2018"){

	leptontriggers_strings = {

	"HLT_Ele32_WPTight_Gsf_L1DoubleEG", //single electron (for 2018 only)
	"DummyBool", //single electron (for 2017 only)
	"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", //double electron (for 2017 and 2018)
	"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", //double electron (for all years)
	"DummyBool", //double muon (for 2016 and 2017)
	"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",  //double muon (for 2017 and 2018)
	"DummyBool", //single muon (for 2017 only)
	"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", //muon+electron (for all years)
	"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", //muon+electron (for all years)
	"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", //muon+electron (for all years)
	"DummyBool", //single electron (for 2016 only)
	"DummyBool", //single electron (for 2016 only)
	"DummyBool", //single electron (for 2016 only)
	"HLT_IsoMu24", //single muon (all years)
	"DummyBool", //single muon (for 2016 only)	
	"DummyBool", //muon+electron (for 2016 only)
	"DummyBool", //muon+electron (for 2016 only)
	"DummyBool", //muon+electron (for 2016 only)
	"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8" //double muon (for 2018 only)

};

} 
else{std::cout << "choose 2016, 2017 or 2018" << std::endl;}






std::vector<std::string> sJER_sigmaJER_strings = {"Jet_eta", "fixedGridRhoFastjetAll", "Jet_pt"};

std::string SJER;
std::string SIGMAJER;

if(JetSmearing_ScaleUp == true){SJER = "sJER_up";} 
else if(JetSmearing_ScaleDown == true){SJER = "sJER_down";}
else{SJER = "sJER_Nominal";}

if(JetResolution_ScaleUp == true){SIGMAJER = "sigma_JER_up";}
else if(JetResolution_ScaleDown == true){SIGMAJER = "sigma_JER_down";}
else{SIGMAJER = "sigma_JER";}


std::vector<std::string> JetSmearingStrings = {"Jet_pt", "Jet_eta", "Jet_phi", "GenJet_pt", "GenJet_eta", "GenJet_phi", SJER, SIGMAJER};
std::vector<std::string> ApplyCJER_strings = {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "cJER", "nJet"};


std::vector<std::string> flags_strings = {
"Flag_goodVertices", 
"Flag_globalSuperTightHalo2016Filter", 
"Flag_HBHENoiseFilter", 
"Flag_HBHENoiseIsoFilter", 
"Flag_EcalDeadCellTriggerPrimitiveFilter", 
"Flag_BadPFMuonFilter", 
"Flag_BadChargedCandidateFilter", 
"Flag_ecalBadCalibFilter", 
"Flag_eeBadScFilter"
};


std::vector<std::string> mumu_strings = {
"nMuon", 
"Muon_tightId_Selection", 
"Muon_pfIsoId_Selection", 
"Muon_eta", 
"Muon_charge_Selection", 
"Muon_pt_Selection", 
"Muon_isPFcand_Selection"};



std::vector<std::string> deltaR_JetE_strings = {
JetEtaInput,
JetPhiInput,
"z_lep_eta",
"z_lep_phi"
};

std::vector<std::string> deltaR_JetMu_strings = {
JetEtaInput,
JetPhiInput,
"z_lep_eta",
"z_lep_phi"
};

std::vector<std::string> deltaR_j1j2_strings = {
"LeadingJetEta",
"LeadingJetPhi",
"SubleadingJetEta",
"SubleadingJetPhi",
};

std::vector<std::string> InvMass_AllJets_strings = {
"LeadingJetpT",
"SubleadingJetpT",
"ThirdJetpT",
"FourthJetpT",
"LeadingJetEta",
"SubleadingJetEta",
"ThirdJetEta",
"FourthJetEta",
"LeadingJetPhi",
"SubleadingJetPhi",
"ThirdJetPhi",
"FourthJetPhi",
"LeadingJetMass",
"SubleadingJetMass",
"ThirdJetMass",
"FourthJetMass",
"nJet"
};


std::vector<std::string> InvMass_3Jets_strings = {
"LeadingJetpT",
"SubleadingJetpT",
"ThirdJetpT",
"LeadingJetEta",
"SubleadingJetEta",
"ThirdJetEta",
"LeadingJetPhi",
"SubleadingJetPhi",
"ThirdJetPhi",
"LeadingJetMass",
"SubleadingJetMass",
"ThirdJetMass",
"nJet"
};

std::vector<std::string> bjet_mass_strings = {JetMassInput, "nJet", "lead_bjet"};
std::vector<std::string> bjet_eta_strings = {JetEtaInput, "nJet", "lead_bjet"};
std::vector<std::string> bjet_pt_strings = {JetPtInput, "nJet", "lead_bjet"};
std::vector<std::string> bjet_phi_strings = {JetPhiInput, "nJet", "lead_bjet"};

std::vector<std::string> nonbjet_mass_strings = {JetMassInput, "nJet", "nonbjets"};
std::vector<std::string> nonbjet_eta_strings = {JetEtaInput, "nJet", "nonbjets"};
std::vector<std::string> nonbjet_pt_strings = {JetPtInput, "nJet", "nonbjets"};
std::vector<std::string> nonbjet_phi_strings = {JetPhiInput, "nJet", "nonbjets"};


std::vector<std::string> w_pair_strings = {
JetPtInput,
JetEtaInput,
JetPhiInput,
JetMassInput,
"Jet_jetId",
"lead_bjet"};


std::vector<std::string> top_strings = {
"bjetpt",
"bjeteta",
"bjetphi",
"bjetmass",
"nJet",
"w_pair_pt", 
"w_pair_eta", 
"w_pair_phi", 
"w_mass",
"w_reco_jets"};

std::vector<std::string> deltaR_Top_LeadingElectron_strings = {
"Top_Phi",
"Top_Eta",
"LeadingElectronEta",
"LeadingElectronPhi"
};

std::vector<std::string> deltaR_Top_LeadingMuon_strings = {
"Top_Phi",
"Top_Eta",
"LeadingMuonEta",
"LeadingMuonPhi"
};

std::vector<std::string> deltaR_Top_SubleadingElectron_strings = {
"Top_Phi",
"Top_Eta",
"SubleadingElectronEta",
"SubleadingElectronPhi"
};

std::vector<std::string> deltaR_Top_SubleadingMuon_strings = {
"Top_Phi",
"Top_Eta",
"SubleadingMuonEta",
"SubleadingMuonPhi"
};

std::vector<std::string> deltaR_Top_LeadingJet_strings = {
"Top_Phi",
"Top_Eta",
"SubleadingJetEta",
"SubleadingJetPhi"
};

std::vector<std::string> deltaR_Top_SubleadingJet_strings = {
"Top_Phi",
"Top_Eta",
"SubleadingJetEta",
"SubleadingJetPhi"
};

std::vector<std::string> deltaR_Top_ThirdJet_strings = {
"Top_Phi",
"Top_Eta",
"ThirdJetEta",
"ThirdJetPhi"
};

std::vector<std::string> deltaR_Top_FourthJet_strings = {
"Top_Phi",
"Top_Eta",
"FourthJetEta",
"FourthJetPhi"
};


std::vector<std::string> deltaR_WJet1_WJet2_strings = {
"WPairJet1Phi",
"WPairJet1Eta",
"WPairJet2Phi",
"WPairJet2Eta",
};

std::vector<std::string> deltaR_Z_Top_strings = {
"RecoZPhi",
"RecoZEta",
"Top_Phi",
"Top_Eta",
};

std::vector<std::string> deltaR_WJet1_LeadingElectron_strings = {
"WPairJet1Phi",
"WPairJet1Eta",
"LeadingElectronPhi",
"LeadingElectronEta"
};

std::vector<std::string> deltaR_WJet1_SubleadingElectron_strings = {
"WPairJet1Phi",
"WPairJet1Eta",
"SubleadingElectronPhi",
"SubleadingElectronEta"
};

std::vector<std::string> deltaR_WJet1_LeadingMuon_strings = {
"WPairJet1Phi",
"WPairJet1Eta",
"LeadingMuonPhi",
"LeadingMuonEta"
};

std::vector<std::string> deltaR_WJet1_SubleadingMuon_strings = {
"WPairJet1Phi",
"WPairJet1Eta",
"SubleadingMuonPhi",
"SubleadingMuonEta"
};


std::vector<std::string> deltaR_WJet1_LeadingJet_strings = {
"WPairJet1Phi",
"WPairJet1Eta",
"LeadingJetPhi",
"LeadingJetEta"
};

std::vector<std::string> deltaR_WJet1_SubleadingJet_strings = {
"WPairJet1Phi",
"WPairJet1Eta",
"SubleadingJetPhi",
"SubleadingJetEta"
};

std::vector<std::string> deltaR_WJet1_ThirdJet_strings = {
"WPairJet1Phi",
"WPairJet1Eta",
"ThirdJetPhi",
"ThirdJetEta"
};

std::vector<std::string> deltaR_WJet1_FourthJet_strings = {
"WPairJet1Phi",
"WPairJet1Eta",
"FourthJetPhi",
"FourthJetEta"
};

std::vector<std::string> deltaR_WJet2_LeadingElectron_strings = {
"WPairJet2Phi",
"WPairJet2Eta",
"LeadingElectronPhi",
"LeadingElectronEta"
};

std::vector<std::string> deltaR_WJet2_SubleadingElectron_strings = {
"WPairJet2Phi",
"WPairJet2Eta",
"SubleadingElectronPhi",
"SubleadingElectronEta"
};


std::vector<std::string> deltaR_WJet2_LeadingMuon_strings = {
"WPairJet2Phi",
"WPairJet2Eta",
"LeadingMuonPhi",
"LeadingMuonEta"
};

std::vector<std::string> deltaR_WJet2_SubleadingMuon_strings = {
"WPairJet2Phi",
"WPairJet2Eta",
"SubleadingMuonPhi",
"SubleadingMuonEta"
};


std::vector<std::string> deltaR_WJet2_LeadingJet_strings = {
"WPairJet2Phi",
"WPairJet2Eta",
"LeadingJetPhi",
"LeadingJetEta"
};

std::vector<std::string> deltaR_WJet2_SubleadingJet_strings = {
"WPairJet2Phi",
"WPairJet2Eta",
"SubleadingJetPhi",
"SubleadingJetEta"
};

std::vector<std::string> deltaR_WJet2_ThirdJet_strings = {
"WPairJet2Phi",
"WPairJet2Eta",
"ThirdJetPhi",
"ThirdJetEta"
};

std::vector<std::string> deltaR_WJet2_FourthJet_strings = {
"WPairJet2Phi",
"WPairJet2Eta",
"FourthJetPhi",
"FourthJetEta"
};


std::vector<std::string> deltaR_Z_LeadingJet_strings = {
"RecoZPhi",
"RecoZEta",
"LeadingJetPhi",
"LeadingJetEta"
};

std::vector<std::string> deltaR_Z_SubleadingJet_strings = {
"RecoZPhi",
"RecoZEta",
"SubleadingJetPhi",
"SubleadingJetEta"
};

std::vector<std::string> deltaR_Z_ThirdJet_strings = {
"RecoZPhi",
"RecoZEta",
"ThirdJetPhi",
"ThirdJetEta"
};

std::vector<std::string> deltaR_Z_FourthJet_strings = {
"RecoZPhi",
"RecoZEta",
"FourthJetPhi",
"FourthJetEta"
};

std::vector<std::string> deltaR_Top_W_strings = {
"w_pair_eta",
"w_pair_phi",
"Top_Eta",
"Top_Phi"
};

std::vector<std::string> deltaR_W_Z_strings = {
"w_pair_phi",
"w_pair_eta",
"RecoZEta",
"RecoZPhi"
};

std::vector<std::string> RecoZstrings_ee{
"nElectron",
"LeadingElectron_pT",
"LeadingElectronEta",
"LeadingElectronPhi",
"LeadingElectronMass",
"SubleadingElectron_pT",
"SubleadingElectronEta",
"SubleadingElectronPhi",
"SubleadingElectronMass"
};

std::vector<std::string> RecoZstrings_mumu = {
"nMuon",
"LeadingMuon_pT",
"LeadingMuonEta",
"LeadingMuonPhi",
"LeadingMuonMass",
"SubleadingMuon_pT",
"SubleadingMuonEta",
"SubleadingMuonPhi",
"SubleadingMuonMass"
};


std::vector<std::string> TotalHT_System_strings = {
"RecoZHT",
"RecoWHT",
"Top_HT",
"TotLepHT",
"TotJetHT"
};

std::vector<std::string> TotalPt_System_strings = {
"RecoZPt",
"w_pair_pt",
"ReweightedTopPt",
"LepPtSum",
"JetPtSum"
};

std::vector<std::string> TotalPt_Systemstrings_data = {
"RecoZPt",
"w_pair_pt",
"Top_Pt",
"LepPtSum",
"JetPtSum"
};


std::vector<std::string> TotalEta_System_strings = {
"RecoZEta",
"w_pair_eta",
"Top_Eta",
"LepEtaSum",
"JetEtaSum"
};


std::vector<std::string> TotalPhi_System_strings = {
"RecoZEta",
"w_pair_eta",
"Top_Eta",
"LepPhiSum",
"JetPhiSum"
};


std::vector<std::string> dR_Z_WPairJet1_strings = {
"RecoZPhi",
"RecoZEta",
"WPairJet1Eta",
"WPairJet1Phi"
};

std::vector<std::string> dR_Z_WPairJet2_strings = {
"RecoZPhi",
"RecoZEta",
"WPairJet2Eta",
"WPairJet2Phi"
};


std::vector<std::string> MinDeltaR_strings = {
"nJet",
"RecoZPhi",
"RecoZEta",
JetPhiInput,
JetEtaInput
};


std::vector<std::string> MinDeltaPhi_strings = {
"nJet",
"RecoZPhi",
JetPhiInput
};

std::vector<std::string> dR_LeadingLepton_LeadingBJetstrings_ee = {
"bjeteta",
"LeadingElectronEta",
"bjetphi",
"LeadingElectronPhi"
};

std::vector<std::string> dR_SubleadingLepton_LeadingBJetstrings_ee = {
"bjeteta",
"SubleadingElectronEta",
"bjetphi",
"SubleadingElectronPhi"
};

std::vector<std::string> dR_LeadingLepton_LeadingBJetstrings_mumu = {
"bjeteta",
"LeadingMuonEta",
"bjetphi",
"LeadingMuonPhi"
};

std::vector<std::string> dR_SubleadingLepton_LeadingBJetstrings_mumu = {
"bjeteta",
"SubleadingMuonEta",
"bjetphi",
"SubleadingMuonPhi"
};


std::vector<std::string> DeltaPhi_Leadinglepton_BJet_string_ee = {
JetPhiInput,
"LeadingElectronPhi"
};

std::vector<std::string> DeltaPhi_Leadinglepton_BJet_string_mumu = {
JetPhiInput,
"LeadingMuonPhi"
};


std::vector<std::string> DeltaPhi_Subleadinglepton_BJet_string_ee = {
JetPhiInput,
"SubleadingElectronPhi"
};

std::vector<std::string> DeltaPhi_Subleadinglepton_BJet_string_mumu = {
JetPhiInput,
"SubleadingMuonPhi"
};



std::vector<std::string> BJetOutputDiscriminantStrings_Leading = {
"LeadingJetpT",
"Jet_btagCSVV2",
"tight_jets",
JetEtaInput
}; 

std::vector<std::string> BJetOutputDiscriminantStrings_Subleading = {
"SubleadingJetpT",
"Jet_btagCSVV2",
"tight_jets",
JetEtaInput
};

std::vector<std::string> BJetOutputDiscriminantStrings_Third = {
"ThirdJetpT",
"Jet_btagCSVV2",
"tight_jets",
JetEtaInput
};

std::vector<std::string> BJetOutputDiscriminantStrings_Fourth = {
"FourthJetpT",
"Jet_btagCSVV2",
"tight_jets",
JetEtaInput
};

std::vector<std::string> lep_cut_ee_strings;
std::vector<std::string> lep_cut_mumu_strings;

if(NPL == true){

	lep_cut_ee_strings = {
	"Electron_pt_Selection", 
	"LooseElectron_pt_Selection", 
	"SameSign",
	"nElectron",
	"LeadingElectron_dz_ECALBarrel",
	"LeadingElectron_dxy_ECALBarrel",
	"LeadingElectron_dz_ECALEndcaps",
	"LeadingElectron_dxy_ECALEndcaps",
	"SubleadingElectron_dz_ECALBarrel",
	"SubleadingElectron_dxy_ECALBarrel",
	"SubleadingElectron_dz_ECALEndcaps",
	"SubleadingElectron_dxy_ECALEndcaps"
	
	};


	lep_cut_mumu_strings = {
	"Muon_pt_Selection",
	"LooseMuon_pt_Selection", 
	"SameSign", 
	"nMuon"
	};

}
else{

	lep_cut_ee_strings = {
	"Electron_pt_Selection",
	"LooseElectron_pt_Selection",
	"OppositeSign",
	"nElectron",
	"LeadingElectron_dz_ECALBarrel",
	"LeadingElectron_dxy_ECALBarrel",
	"LeadingElectron_dz_ECALEndcaps",
	"LeadingElectron_dxy_ECALEndcaps",
	"SubleadingElectron_dz_ECALBarrel",
	"SubleadingElectron_dxy_ECALBarrel",
	"SubleadingElectron_dz_ECALEndcaps",
	"SubleadingElectron_dxy_ECALEndcaps"
	
	};

	lep_cut_mumu_strings = {
	"Muon_pt_Selection",  
	"LooseMuon_pt_Selection", 
	"OppositeSign", 
	"nMuon"
	};

}


std::vector<std::string> lep_cut_emu_strings = {
"Electron_pt_SelectionEmu",
"LooseElectron_pt_SelectionEmu",
"OppositeSign_emu",
"nElectron",
"Muon_pt_SelectionEmu", 
"LooseMuon_pt_SelectionEmu",
"nMuon"

};


//Event cleaning
auto d_EventCleaning = d_dataframe.Filter(filter_function, flags_strings, "Event cleaning filter");

//Pile up modelling
//2016
TFile *dataPileupFile_2016 = new TFile("./ScaleFactors/PileUp/2016/truePileupTest.root", "READ");
TH1F *dataPU_2016 = (TH1F*)(dataPileupFile_2016->Get("pileup")->Clone());
TFile *mcPileupFile_2016 = new TFile("./ScaleFactors/PileUp/2016/pileupMC.root", "READ");
TH1F* mcPU_2016 = (TH1F*)(mcPileupFile_2016->Get("pileup")->Clone());

//2016 part 1
TFile *dataPileupFile_2016_part1 = new TFile("./ScaleFactors/PileUp/2016/truePileupTest_part1.root", "READ");
TH1F *dataPU_2016_part1 = (TH1F*)(dataPileupFile_2016_part1->Get("pileup")->Clone());
TFile *mcPileupFile_2016_part1 = new TFile("./ScaleFactors/PileUp/2016/pileupMC.root", "READ");
TH1F* mcPU_2016_part1 = (TH1F*)(mcPileupFile_2016_part1->Get("pileup")->Clone());

//2016 part 2
TFile *dataPileupFile_2016_part2 = new TFile("./ScaleFactors/PileUp/2016/truePileupTest_part2.root", "READ");
TH1F *dataPU_2016_part2 = (TH1F*)(dataPileupFile_2016_part2->Get("pileup")->Clone());
TFile *mcPileupFile_2016_part2 = new TFile("./ScaleFactors/PileUp/2016/pileupMC.root", "READ");
TH1F* mcPU_2016_part2 = (TH1F*)(mcPileupFile_2016_part2->Get("pileup")->Clone());

//2017
TFile *dataPileupFile_2017 = new TFile("./ScaleFactors/PileUp/2017/truePileupTest.root", "READ");
TH1F *dataPU_2017 = (TH1F*)(dataPileupFile_2017->Get("pileup")->Clone());
TFile *mcPileupFile_2017 = new TFile("./ScaleFactors/PileUp/2017/pileupMC.root", "READ");
TH1F* mcPU_2017 = (TH1F*)(mcPileupFile_2017->Get("pileup")->Clone());


//2018
TFile *dataPileupFile_2018 = new TFile("./ScaleFactors/PileUp/2018/MyDataPileupHistogram2018.root", "READ");
TH1F *dataPU_2018 = (TH1F*)(dataPileupFile_2018->Get("pileup")->Clone());
TFile *mcPileupFile_2018 = new TFile("./ScaleFactors/PileUp/2018/pileupMC2018.root", "READ");
TH1F* mcPU_2018 = (TH1F*)(mcPileupFile_2018->Get("pileup")->Clone());




//Systematic files
//2016
TFile *systUpFile_2016 = new TFile("./ScaleFactors/PileUp/2016/truePileupUp.root", "READ");
TH1F *pileupUpHist_2016 = (TH1F*)(systUpFile_2016->Get("pileup")->Clone());
TFile *systDownFile_2016 = new TFile("./ScaleFactors/PileUp/2016/truePileupDown.root", "READ");
TH1F *pileupDownHist_2016 = (TH1F*)(systDownFile_2016->Get("pileup")->Clone());
//part 1
TFile *systUpFile_2016_part1 = new TFile("./ScaleFactors/PileUp/2016/truePileupUp_part1.root", "READ");
TH1F *pileupUpHist_2016_part1 = (TH1F*)(systUpFile_2016_part1->Get("pileup")->Clone());
TFile *systDownFile_2016_part1 = new TFile("./ScaleFactors/PileUp/2016/truePileupDown_part1.root", "READ");
TH1F *pileupDownHist_2016_part1 = (TH1F*)(systDownFile_2016_part1->Get("pileup")->Clone());

//part 2
TFile *systUpFile_2016_part2 = new TFile("./ScaleFactors/PileUp/2016/truePileupUp_part2.root", "READ");
TH1F *pileupUpHist_2016_part2 = (TH1F*)(systUpFile_2016_part2->Get("pileup")->Clone());
TFile *systDownFile_2016_part2 = new TFile("./ScaleFactors/PileUp/2016/truePileupDown_part2.root", "READ");
TH1F *pileupDownHist_2016_part2 = (TH1F*)(systDownFile_2016_part2->Get("pileup")->Clone());


TH1F *puReweight_2016 = (TH1F*)(dataPU_2016->Clone());
puReweight_2016->Scale(1.0 / puReweight_2016->Integral());
mcPU_2016->Scale(1.0 / mcPU_2016->Integral());
puReweight_2016->Divide(mcPU_2016);
puReweight_2016->SetDirectory(nullptr);

TH1F *puReweight_2016_part1 = (TH1F*)(dataPU_2016_part1->Clone());
puReweight_2016_part1->Scale(1.0 / puReweight_2016_part1->Integral());
mcPU_2016_part1->Scale(1.0 / mcPU_2016_part1->Integral());
puReweight_2016_part1->Divide(mcPU_2016_part1);
puReweight_2016_part1->SetDirectory(nullptr);

TH1F *puReweight_2016_part2 = (TH1F*)(dataPU_2016_part2->Clone());
puReweight_2016_part2->Scale(1.0 / puReweight_2016_part2->Integral());
mcPU_2016_part2->Scale(1.0 / mcPU_2016_part2->Integral());
puReweight_2016_part2->Divide(mcPU_2016_part2);
puReweight_2016_part2->SetDirectory(nullptr);


// 2017
TFile *systUpFile_2017 = new TFile("./ScaleFactors/PileUp/2017/truePileupUp.root", "READ");
TH1F *pileupUpHist_2017 = (TH1F*)(systUpFile_2017->Get("pileup")->Clone());
TFile *systDownFile_2017 = new TFile("./ScaleFactors/PileUp/2017/truePileupDown.root", "READ");
TH1F *pileupDownHist_2017 = (TH1F*)(systDownFile_2017->Get("pileup")->Clone());


TH1F *puReweight_2017 = (TH1F*)(dataPU_2017->Clone());
puReweight_2017->Scale(1.0 / puReweight_2017->Integral());
mcPU_2017->Scale(1.0 / mcPU_2017->Integral());
puReweight_2017->Divide(mcPU_2017);
puReweight_2017->SetDirectory(nullptr);

//2018
TFile *systUpFile_2018 = new TFile("./ScaleFactors/PileUp/2018/MyDataPileupHistogramScaleUp2018.root", "READ");
TH1F *pileupUpHist_2018 = (TH1F*)(systUpFile_2018->Get("pileup")->Clone());
TFile *systDownFile_2018 = new TFile("./ScaleFactors/PileUp/2018/MyDataPileupHistogramScaleDown2018.root", "READ");
TH1F *pileupDownHist_2018 = (TH1F*)(systDownFile_2018->Get("pileup")->Clone());


TH1F *puReweight_2018 = (TH1F*)(dataPU_2018->Clone());
puReweight_2018->Scale(1.0 / puReweight_2018->Integral());
mcPU_2018->Scale(1.0 / mcPU_2018->Integral());
puReweight_2018->Divide(mcPU_2018);
puReweight_2018->SetDirectory(nullptr);



///Systematic sample
//2016
TH1F *puSystUp_2016 = (TH1F*)(pileupUpHist_2016->Clone());
puSystUp_2016->Scale(1.0 / puSystUp_2016->Integral());
puSystUp_2016->Divide(mcPU_2016);
puSystUp_2016->SetDirectory(nullptr);
TH1F *puSystDown_2016 = (TH1F*)(pileupDownHist_2016->Clone());
puSystDown_2016->Scale(1.0 / puSystDown_2016->Integral());
puSystDown_2016->Divide(mcPU_2016);
puSystDown_2016->SetDirectory(nullptr);

//2016 part 1
TH1F *puSystUp_2016_part1 = (TH1F*)(pileupUpHist_2016_part1->Clone());
puSystUp_2016_part1->Scale(1.0 / puSystUp_2016_part1->Integral());
puSystUp_2016_part1->Divide(mcPU_2016_part1);
puSystUp_2016_part1->SetDirectory(nullptr);
TH1F *puSystDown_2016_part1 = (TH1F*)(pileupDownHist_2016_part1->Clone());
puSystDown_2016_part1->Scale(1.0 / puSystDown_2016_part1->Integral());
puSystDown_2016_part1->Divide(mcPU_2016_part1);
puSystDown_2016_part1->SetDirectory(nullptr);

//2016 part 2
TH1F *puSystUp_2016_part2 = (TH1F*)(pileupUpHist_2016_part2->Clone());
puSystUp_2016_part2->Scale(1.0 / puSystUp_2016_part2->Integral());
puSystUp_2016_part2->Divide(mcPU_2016_part2);
puSystUp_2016_part2->SetDirectory(nullptr);
TH1F *puSystDown_2016_part2 = (TH1F*)(pileupDownHist_2016_part2->Clone());
puSystDown_2016_part2->Scale(1.0 / puSystDown_2016_part2->Integral());
puSystDown_2016_part2->Divide(mcPU_2016_part2);
puSystDown_2016_part2->SetDirectory(nullptr);

//2017
TH1F *puSystUp_2017 = (TH1F*)(pileupUpHist_2017->Clone());
puSystUp_2017->Scale(1.0 / puSystUp_2017->Integral());
puSystUp_2017->Divide(mcPU_2017);
puSystUp_2017->SetDirectory(nullptr);
TH1F *puSystDown_2017 = (TH1F*)(pileupDownHist_2017->Clone());
puSystDown_2017->Scale(1.0 / puSystDown_2017->Integral());
puSystDown_2017->Divide(mcPU_2017);
puSystDown_2017->SetDirectory(nullptr);

//2018
TH1F *puSystUp_2018 = (TH1F*)(pileupUpHist_2018->Clone());
puSystUp_2018->Scale(1.0 / puSystUp_2018->Integral());
puSystUp_2018->Divide(mcPU_2018);
puSystUp_2018->SetDirectory(nullptr);
TH1F *puSystDown_2018 = (TH1F*)(pileupDownHist_2018->Clone());
puSystDown_2018->Scale(1.0 / puSystDown_2018->Integral());
puSystDown_2018->Divide(mcPU_2018);
puSystDown_2018->SetDirectory(nullptr);


//Closing the pile up files
dataPileupFile_2016->Close();
mcPileupFile_2016->Close();
systUpFile_2016->Close();
systDownFile_2016->Close();

dataPileupFile_2016_part1->Close();
mcPileupFile_2016_part1->Close();
systUpFile_2016_part1->Close();
systDownFile_2016_part1->Close();

dataPileupFile_2016_part2->Close();
mcPileupFile_2016_part2->Close();
systUpFile_2016_part2->Close();
systDownFile_2016_part2->Close();

dataPileupFile_2017->Close();
mcPileupFile_2017->Close();
systUpFile_2017->Close();
systDownFile_2017->Close();

dataPileupFile_2018->Close();
mcPileupFile_2018->Close();
systUpFile_2018->Close();
systDownFile_2018->Close();


//Implementing the PU modelling


auto PU_function{[&puReweight_2016, &puReweight_2016_part1, &puReweight_2016_part2, &puReweight_2017, &puReweight_2018, &year](int PV_npvs_input){

  float PU_Weight_input;

  if(year == "2016"){
        PU_Weight_input = puReweight_2016->GetBinContent(puReweight_2016->GetXaxis()->FindBin(PV_npvs_input));
  }
  else if(year == "2017"){
  	PU_Weight_input = puReweight_2017->GetBinContent(puReweight_2017->GetXaxis()->FindBin(PV_npvs_input));
  }
  else if(year == "2018"){
        PU_Weight_input = puReweight_2018->GetBinContent(puReweight_2018->GetXaxis()->FindBin(PV_npvs_input));
  }
  else{std::cout << "Choose a year out of 2016, 2017 or 2018 for the PU function" << std::endl;}


  return PU_Weight_input;

}};



///MC starts here

if(

process != "data_DoubleEGRunB" &&
process != "data_DoubleEGRunC" &&
process != "data_DoubleEGRunD" &&
process != "data_DoubleEGRunE" &&
process != "data_DoubleEGRunF" &&
process != "data_DoubleEGRunG" &&
process != "data_DoubleEGRunH" &&
process != "data_SingleElectronRunB" &&
process != "data_SingleElectronRunC" &&
process != "data_SingleElectronRunD" &&
process != "data_SingleElectronRunE" &&
process != "data_SingleElectronRunF" &&
process != "data_SingleElectronRunG" &&
process != "data_SingleElectronRunH" &&
process != "data_DoubleMuonRunB" &&
process != "data_DoubleMuonRunC" &&
process != "data_DoubleMuonRunD" &&
process != "data_DoubleMuonRunE" &&
process != "data_DoubleMuonRunF" &&
process != "data_DoubleMuonRunG" &&
process != "data_DoubleMuonRunH" &&
process != "data_SingleMuonRunB" &&
process != "data_SingleMuonRunC" &&
process != "data_SingleMuonRunD" &&
process != "data_SingleMuonRunE" &&
process != "data_SingleMuonRunF" &&
process != "data_SingleMuonRunG" &&
process != "data_SingleMuonRunH" 

){
  
  //Filtering events that pass the ee selection criteria
  auto d_ee_selection_defines = d_EventCleaning.Define("DummyBool", DummyBool, {"HLT_PFHT250"})
					       .Define("PU", PU_function, {"PV_npvs"})
                                           .Define("TightElectrons", TightElectronsFunction, {"Electron_pt", "Electron_eta", "Electron_cutBased", "Electron_isPFcand"})
                                          .Define("Electron_pt_Selection", select<floats>, {"Electron_pt", "TightElectrons"})
					  .Define("Electron_phi_Selection", select<floats>, {"Electron_phi", "TightElectrons"})
					  .Define("Electron_eta_Selection", select<floats>, {"Electron_eta", "TightElectrons"})
                                          .Define("Electron_charge_Selection", select<ints>, {"Electron_charge", "TightElectrons"})
                                          .Define("LooseElectrons", LooseElectronsFunction, {"Electron_pt", "Electron_eta", "Electron_cutBased", "Electron_isPFcand"})
                                          .Define("LooseElectron_pt_Selection", select<floats>, {"Electron_pt", "LooseElectrons"})
                                          .Define("LooseElectron_charge_Selection", select<ints>, {"Electron_charge", "LooseElectrons"})
                                          .Define("OppositeSign", OppositeSign, {"Electron_charge_Selection"})
                                          .Define("SameSign", SameSign, {"Electron_charge_Selection"})
                                          .Define("LeadingElectron_pT", LeadingVariable, {"Electron_pt_Selection"})
                                          .Define("SubleadingElectron_pT", SubleadingVariable, {"Electron_pt_Selection"})
                                          .Define("LeadingElectronPhi", LeadingVariable, {"Electron_phi_Selection"})
                                          .Define("SubleadingElectronPhi", SubleadingVariable, {"Electron_phi_Selection"})
					  .Define("LeadingElectronEta", LeadingVariable, {"Electron_eta_Selection"})
                                          .Define("SubleadingElectronEta", SubleadingVariable, {"Electron_eta_Selection"})
                                          .Define("LeadingElectronMass", LeadingVariable, {"Electron_mass"})
                                          .Define("SubleadingElectronMass", SubleadingVariable, {"Electron_mass"})
                                          .Define("LeadingElectron_RelIso_Selection", LeadingVariable, {"Electron_jetRelIso"})
                                          .Define("SubleadingElectron_RelIso_Selection", SubleadingVariable, {"Electron_jetRelIso"})
                                          .Define("LeadingElectron_dz_ECALBarrel", LeadingElectron_dz_ECALBarrel_function, {"LeadingElectron_pT", "Electron_eta_Selection", "Electron_dz"})
                                          .Define("LeadingElectron_dxy_ECALBarrel", LeadingElectron_dxy_ECALBarrel_function, {"LeadingElectron_pT", "Electron_eta_Selection", "Electron_dxy"})
                                          .Define("LeadingElectron_dz_ECALEndcaps", LeadingElectron_dz_ECALEndcaps_function, {"LeadingElectron_pT", "Electron_eta_Selection", "Electron_dz"})
                                          .Define("LeadingElectron_dxy_ECALEndcaps", LeadingElectron_dxy_ECALEndcaps_function, {"LeadingElectron_pT", "Electron_eta_Selection", "Electron_dxy"})
                                          .Define("SubleadingElectron_dz_ECALBarrel", SubleadingElectron_dz_ECALBarrel_function, {"SubleadingElectron_pT", "Electron_eta_Selection", "Electron_dz"})
                                          .Define("SubleadingElectron_dxy_ECALBarrel", SubleadingElectron_dxy_ECALBarrel_function, {"SubleadingElectron_pT", "Electron_eta_Selection", "Electron_dxy"})
                                          .Define("SubleadingElectron_dz_ECALEndcaps", SubleadingElectron_dz_ECALEndcaps_function, {"SubleadingElectron_pT", "Electron_eta_Selection", "Electron_dz"})
                                          .Define("SubleadingElectron_dxy_ECALEndcaps", SubleadingElectron_dxy_ECALEndcaps_function, {"SubleadingElectron_pT", "Electron_eta_Selection", "Electron_dxy"});




auto d_mumu_selection_defines = d_EventCleaning.Define("DummyBool", DummyBool, {"HLT_PFHT250"})
					       .Define("PU", PU_function, {"PV_npvs"})
                                               .Define("TightMuons", TightMuonsFunction, {"Muon_isPFcand", "Muon_pt", "Muon_eta", "Muon_tightId", "Muon_pfRelIso04_all"})
                                               .Define("Muon_pt_Selection", select<floats>, {"Muon_pt", "TightMuons"})
                                               .Define("Muon_eta_Selection", select<floats>, {"Muon_eta", "TightMuons"})
                                               .Define("Muon_phi_Selection", select<floats>, {"Muon_phi", "TightMuons"})
                                               .Define("Muon_mass_Selection", select<floats>, {"Muon_mass", "TightMuons"})
                                               .Define("Muon_charge_Selection", select<ints>, {"Muon_charge", "TightMuons"})
                                               .Define("LooseMuons", LooseMuonsFunction, {"Muon_isPFcand", "Muon_pt", "Muon_eta", "Muon_softId", "Muon_pfRelIso04_all"})
                                               .Define("LooseMuon_pt_Selection", select<floats>, {"Muon_pt", "LooseMuons"})
                                               .Define("LooseMuon_charge_Selection", select<ints>, {"Muon_charge", "LooseMuons"})
					       .Define("OppositeSign", OppositeSign, {"Muon_charge_Selection"})
                                               .Define("SameSign", SameSign, {"Muon_charge_Selection"})
                                               .Define("LeadingMuon_pT", LeadingVariable, {"Muon_pt_Selection"})
                                               .Define("SubleadingMuon_pT", SubleadingVariable, {"Muon_pt_Selection"})
                                               .Define("LeadingMuonPhi", LeadingVariable, {"Muon_phi_Selection"})
                                               .Define("SubleadingMuonPhi", SubleadingVariable, {"Muon_phi_Selection"})
                                               .Define("LeadingMuonMass", LeadingVariable, {"Muon_mass_Selection"})
                                               .Define("SubleadingMuonMass", SubleadingVariable, {"Muon_mass_Selection"})
                                               .Define("LeadingMuon_RelIso_Selection", LeadingVariable, {"Muon_jetRelIso"})
                                               .Define("SubleadingMuon_RelIso_Selection", SubleadingVariable, {"Muon_jetRelIso"})
                                               .Define("LeadingMuonEta", LeadingVariable, {"Muon_eta_Selection"})
                                               .Define("SubleadingMuonEta", SubleadingVariable, {"Muon_eta_Selection"});
        

     auto d_emu_selection_defines = d_EventCleaning.Define("DummyBool", DummyBool, {"HLT_PFHT250"})
						   .Define("PU", PU_function, {"PV_npvs"})
                                           	   .Define("TightElectronsEmu", TightElectronsFunctionEmu, {"Electron_pt", "Electron_eta", "Electron_cutBased", "Electron_isPFcand"})
						   .Define("LooseMuonsEmu", LooseMuonsFunctionEmu, {"Muon_isPFcand", "Muon_pt", "Muon_eta", "Muon_softId", "Muon_pfRelIso04_all"})
                                          	   .Define("Electron_pt_SelectionEmu", select<floats>, {"Electron_pt", "TightElectronsEmu"})
						   .Define("Electron_eta_SelectionEmu", select<floats>, {"Electron_eta", "TightElectronsEmu"})
                                          	   .Define("Electron_charge_SelectionEmu", select<ints>, {"Electron_charge", "TightElectronsEmu"})
                                                   .Define("LooseElectronsEmu", LooseElectronsFunctionEmu, {"Electron_pt", "Electron_eta", "Electron_cutBased", "Electron_isPFcand"})
                                          	   .Define("LooseElectron_pt_SelectionEmu", select<floats>, {"Electron_pt", "LooseElectronsEmu"})
                                          	   .Define("LooseElectron_charge_SelectionEmu", select<ints>, {"Electron_charge", "LooseElectronsEmu"})
						   .Define("Electron_phi_SelectionEmu", select<floats>, {"Electron_phi", "TightElectronsEmu"})
						   .Define("TightMuonsEmu", TightMuonsFunctionEmu, {"Muon_isPFcand", "Muon_pt", "Muon_eta", "Muon_tightId", "Muon_pfRelIso04_all"})
                                                   .Define("Muon_pt_SelectionEmu", select<floats>, {"Muon_pt", "TightMuonsEmu"})
                                            	   .Define("Muon_eta_SelectionEmu", select<floats>, {"Muon_eta", "TightMuonsEmu"})
                                             	   .Define("Muon_phi_SelectionEmu", select<floats>, {"Muon_phi", "TightMuonsEmu"})
                                                   .Define("Muon_mass_SelectionEmu", select<floats>, {"Muon_mass", "TightMuonsEmu"})
                                             	   .Define("Muon_charge_SelectionEmu", select<ints>, {"Muon_charge", "TightMuonsEmu"})
                                             	   .Define("LooseMuon_pt_SelectionEmu", select<floats>, {"Muon_pt", "LooseMuonsEmu"})
                                             	   .Define("LooseMuon_charge_SelectionEmu", select<ints>, {"Muon_charge", "LooseMuonsEmu"})
                                                   .Define("OppositeSign_emu", OppositeSign_emu, {"Electron_charge_SelectionEmu", "Muon_charge_SelectionEmu"})
                                             	   .Define("LeadingElectronMuon_pT", LeadingVariableEmu, {"Electron_pt_SelectionEmu", "Muon_pt_SelectionEmu"})
                                            	   .Define("SubleadingElectronMuon_pT", SubleadingVariableEmu, {"Electron_pt_SelectionEmu", "Muon_pt_SelectionEmu"})
                                             	   .Define("LeadingElectronMuonPhi", LeadingVariableEmu, {"Electron_phi_SelectionEmu", "Muon_phi_SelectionEmu"})
                                                   .Define("SubleadingElectronMuonPhi", SubleadingVariableEmu, {"Electron_phi_SelectionEmu", "Muon_phi_SelectionEmu"})
                                             	   .Define("LeadingElectronMuonEta", LeadingVariableEmu, {"Electron_eta_SelectionEmu", "Muon_eta_SelectionEmu"})
                                             	   .Define("SubleadingElectronMuonEta", SubleadingVariableEmu, {"Electron_eta_SelectionEmu", "Muon_eta_SelectionEmu"});



auto d_ee_selection = d_ee_selection_defines.Filter(lep_cut_ee, lep_cut_ee_strings, "lepton cut (ee)");
auto d_mumu_selection = d_mumu_selection_defines.Filter(lep_cut_mumu, lep_cut_mumu_strings, "lepton cut (mumu)");


/*

if(ttbarCR == false && process != "Data_triggerSF" && process != "MC_triggerSF"){

	d_ee_selection = d_ee_selection_defines.Filter(lep_cut_ee, lep_cut_ee_strings, "lepton cut (ee)");
	d_mumu_selection = d_mumu_selection_defines.Filter(lep_cut_mumu, lep_cut_mumu_strings, "lepton cut (mumu)");

}
else if(ttbarCR == false && (process == "Data_triggerSF" || process == "MC_triggerSF")){
	
        d_ee_selection = d_ee_selection_defines.Filter(lep_cut_ee, lep_cut_ee_strings, "lepton cut (ee)");
        d_mumu_selection = d_mumu_selection_defines.Filter(lep_cut_mumu, lep_cut_mumu_strings, "lepton cut (mumu)");

}
else if(ttbarCR == true){

       auto d_ee_selection = d_emu_selection_defines.Filter(lep_cut_emu, lep_cut_emu_strings, "lepton cut (emu)");
       auto d_mumu_selection = d_emu_selection_defines.Filter(lep_cut_emu, lep_cut_emu_strings, "lepton cut (emu)");

}
else{std::cout << "please choose ttbarCR as true or false" << std::endl;}
*/



  if( (process == "MC_triggerSF_ttbar" || process == "MC_triggerSF_ZPlusJets" || process == "Data_triggerSF") ){

	auto d_emu_selection = d_emu_selection_defines.Filter(lep_cut_emu, lep_cut_emu_strings, "lepton cut (emu)");

  	//Filtering events that pass the MET filters and the selection criteria (for trigger SF calculation)
  	auto d_MET_And_Selection_ee = d_ee_selection.Filter(MET_triggers_function, MET_triggers);
  	auto d_MET_And_Selection_mumu = d_mumu_selection.Filter(MET_triggers_function, MET_triggers);
  	auto d_MET_And_Selection_emu = d_emu_selection.Filter(MET_triggers_function, MET_triggers);


  	//Filtering events that pass the dilepton triggers and the selection criteria (for trigger SF calculation)
  	auto d_DileptonTriggers_And_Selection_ee = d_ee_selection.Filter(ee_selection_LL_Trig_function, leptontriggers_strings);
  	auto d_DileptonTriggers_And_Selection_mumu = d_mumu_selection.Filter(mumu_selection_LL_Trig_function, leptontriggers_strings);
  	auto d_DileptonTriggers_And_Selection_emu = d_emu_selection.Filter(emu_selection_LL_Trig_function, leptontriggers_strings);

  	//Filtering events that pass the selection criteria, MET filters and dilepton triggers (for trigger SF calculation)
  	auto d_MET_DileptonTriggers_And_Selection_ee = d_DileptonTriggers_And_Selection_ee.Filter(MET_triggers_function, MET_triggers);
  	auto d_MET_DileptonTriggers_And_Selection_mumu = d_DileptonTriggers_And_Selection_mumu.Filter(MET_triggers_function, MET_triggers);
  	auto d_MET_DileptonTriggers_And_Selection_emu = d_DileptonTriggers_And_Selection_emu.Filter(MET_triggers_function, MET_triggers);

        
  	//counting numbers of entries that pass the selection criteria (for trigger SF calculation)
  	float TriggerSF_EventWeight_SelectionCriteria_ee, TriggerSF_EventWeight_DileptonAndSelectionCriteria_ee, TriggerSF_EventWeight_METAndSelectionCriteria_ee, TriggerSF_EventWeight_METDileptonAndSelectionCriteria_ee;

	float TriggerSF_EventWeight_SelectionCriteria_mumu, TriggerSF_EventWeight_DileptonAndSelectionCriteria_mumu, TriggerSF_EventWeight_METAndSelectionCriteria_mumu, TriggerSF_EventWeight_METDileptonAndSelectionCriteria_mumu;
	
	float TriggerSF_EventWeight_SelectionCriteria_emu, TriggerSF_EventWeight_DileptonAndSelectionCriteria_emu, TriggerSF_EventWeight_METAndSelectionCriteria_emu, TriggerSF_EventWeight_METDileptonAndSelectionCriteria_emu;

 	if(process == "MC_triggerSF_ttbar" || "MC_triggerSF_ZPlusJets"){

		TriggerSF_EventWeight_SelectionCriteria_ee = *d_ee_selection.Sum("PU");
		TriggerSF_EventWeight_DileptonAndSelectionCriteria_ee = *d_DileptonTriggers_And_Selection_ee.Sum("PU");
		TriggerSF_EventWeight_METAndSelectionCriteria_ee = *d_MET_And_Selection_ee.Sum("PU");
		TriggerSF_EventWeight_METDileptonAndSelectionCriteria_ee = *d_MET_DileptonTriggers_And_Selection_ee.Sum("PU");

		TriggerSF_EventWeight_SelectionCriteria_mumu = *d_mumu_selection.Sum("PU");
                TriggerSF_EventWeight_DileptonAndSelectionCriteria_mumu = *d_DileptonTriggers_And_Selection_mumu.Sum("PU");
                TriggerSF_EventWeight_METAndSelectionCriteria_mumu = *d_MET_And_Selection_mumu.Sum("PU");
                TriggerSF_EventWeight_METDileptonAndSelectionCriteria_mumu = *d_MET_DileptonTriggers_And_Selection_mumu.Sum("PU");


		TriggerSF_EventWeight_SelectionCriteria_emu = *d_emu_selection.Sum("PU");
                TriggerSF_EventWeight_DileptonAndSelectionCriteria_emu = *d_DileptonTriggers_And_Selection_emu.Sum("PU");
                TriggerSF_EventWeight_METAndSelectionCriteria_emu = *d_MET_And_Selection_emu.Sum("PU");
                TriggerSF_EventWeight_METDileptonAndSelectionCriteria_emu = *d_MET_DileptonTriggers_And_Selection_emu.Sum("PU");

	}
	else{

		TriggerSF_EventWeight_SelectionCriteria_ee = 1.0;
                TriggerSF_EventWeight_DileptonAndSelectionCriteria_ee = 1.0;
                TriggerSF_EventWeight_METAndSelectionCriteria_ee = 1.0;
                TriggerSF_EventWeight_METDileptonAndSelectionCriteria_ee = 1.0;


		TriggerSF_EventWeight_SelectionCriteria_mumu = 1.0;
                TriggerSF_EventWeight_DileptonAndSelectionCriteria_mumu = 1.0;
                TriggerSF_EventWeight_METAndSelectionCriteria_mumu = 1.0;
                TriggerSF_EventWeight_METDileptonAndSelectionCriteria_mumu = 1.0;

		TriggerSF_EventWeight_SelectionCriteria_emu = 1.0;
                TriggerSF_EventWeight_DileptonAndSelectionCriteria_emu = 1.0;
                TriggerSF_EventWeight_METAndSelectionCriteria_emu = 1.0;
                TriggerSF_EventWeight_METDileptonAndSelectionCriteria_emu = 1.0;


        }
  	

  	float N_SelectionCriteria_ee = ( *d_ee_selection.Count() ) * TriggerSF_EventWeight_SelectionCriteria_ee;
  	float N_SelectionCriteria_mumu = ( *d_mumu_selection.Count() ) * TriggerSF_EventWeight_SelectionCriteria_mumu;
  	float N_SelectionCriteria_emu = ( *d_emu_selection.Count() ) * TriggerSF_EventWeight_SelectionCriteria_emu;


  	//counting the numbers of entries that pass the selection criteria and dilepton triggers (for trigger SF calculation)
  	float N_Dilepton_And_SelectionCriteria_ee = ( *d_DileptonTriggers_And_Selection_ee.Count() ) * TriggerSF_EventWeight_DileptonAndSelectionCriteria_ee;
  	float N_Dilepton_And_SelectionCriteria_mumu = ( *d_DileptonTriggers_And_Selection_mumu.Count() ) * TriggerSF_EventWeight_DileptonAndSelectionCriteria_mumu;
  	float N_Dilepton_And_SelectionCriteria_emu = ( *d_DileptonTriggers_And_Selection_emu.Count() ) * TriggerSF_EventWeight_DileptonAndSelectionCriteria_emu;


  	//counting the numbers of entries that pass the selection criteria and the MET filters (for trigger SF calculation)
	float  N_MET_And_SelectionCriteria_ee = ( *d_MET_And_Selection_ee.Count() ) * TriggerSF_EventWeight_METAndSelectionCriteria_ee;
	float N_MET_And_SelectionCriteria_mumu = ( *d_MET_And_Selection_mumu.Count() ) * TriggerSF_EventWeight_METAndSelectionCriteria_mumu;
        float N_MET_And_SelectionCriteria_emu = ( *d_MET_And_Selection_emu.Count() ) * TriggerSF_EventWeight_METAndSelectionCriteria_emu;

  	//counting the numbers of entries that pass the selection criteria and the dilepton triggers (for trigger SF calculation)
	float N_MET_DileptonTriggers_And_SelectionCriteria_ee = ( *d_MET_DileptonTriggers_And_Selection_ee.Count() ) * TriggerSF_EventWeight_METDileptonAndSelectionCriteria_ee; 
	float N_MET_DileptonTriggers_And_SelectionCriteria_mumu = ( *d_MET_DileptonTriggers_And_Selection_mumu.Count() ) * TriggerSF_EventWeight_METDileptonAndSelectionCriteria_mumu; 
	float N_MET_DileptonTriggers_And_SelectionCriteria_emu = ( *d_MET_DileptonTriggers_And_Selection_emu.Count() ) * TriggerSF_EventWeight_METDileptonAndSelectionCriteria_emu;  

     

	//Calculating the efficiencies
	float Eff_ee = (N_MET_DileptonTriggers_And_SelectionCriteria_ee) / (N_MET_And_SelectionCriteria_ee + 1.0e-06); //1.0e-06 to prevent nan values
	float Eff_mumu = (N_MET_DileptonTriggers_And_SelectionCriteria_mumu) / (N_MET_And_SelectionCriteria_mumu + 1.0e-06); //1.0e-06 to prevent nan values
	float Eff_emu = (N_MET_DileptonTriggers_And_SelectionCriteria_emu) / (N_MET_And_SelectionCriteria_emu + 1.0e-06); //1.0e-06 to prevent nan values

	
        std::string TriggerSF_Efficiency;

        if(process == "Data_triggerSF"){ TriggerSF_Efficiency = "TriggerSF_Efficiency_Data_MET_" + year + ".txt";}
	else if(process == "MC_triggerSF_ttbar"){ TriggerSF_Efficiency = "TriggerSF_Efficiency_MC_ttbar_" + year + ".txt";}
	else if(process == "MC_triggerSF_ZPlusJets"){ TriggerSF_Efficiency = "TriggerSF_Efficiency_MC_ZPlusJets_" + year + ".txt";}

	std::ofstream TriggerSF_Efficiency_File;
	TriggerSF_Efficiency_File.open(TriggerSF_Efficiency.c_str());	

	TriggerSF_Efficiency_File << Eff_ee << '\n'
				  << Eff_mumu << '\n'
				  << Eff_emu << '\n' << std::endl;


	//alpha values
        float Eff_MET_Dilepton_Selection_ee = (N_MET_DileptonTriggers_And_SelectionCriteria_ee) / (N_SelectionCriteria_ee + 1.0e-06); 
	float Eff_MET_Dilepton_Selection_mumu = (N_MET_DileptonTriggers_And_SelectionCriteria_mumu) / (N_SelectionCriteria_mumu + 1.0e-06);
	float Eff_MET_Dilepton_Selection_emu = (N_MET_DileptonTriggers_And_SelectionCriteria_emu) / (N_SelectionCriteria_emu + 1.0e-06);

	float Eff_Dilepton_Selection_ee = (N_Dilepton_And_SelectionCriteria_ee) / (N_SelectionCriteria_ee + 1.0e-06);
        float Eff_Dilepton_Selection_mumu = (N_Dilepton_And_SelectionCriteria_mumu) / (N_SelectionCriteria_mumu + 1.0e-06);
        float Eff_Dilepton_Selection_emu = (N_Dilepton_And_SelectionCriteria_emu) / (N_SelectionCriteria_emu + 1.0e-06);

	float Eff_MET_Selection_ee = (N_MET_And_SelectionCriteria_ee) / (N_SelectionCriteria_ee + 1.0e-06);
        float Eff_MET_Selection_mumu = (N_MET_And_SelectionCriteria_mumu) / (N_SelectionCriteria_mumu + 1.0e-06);
        float Eff_MET_Selection_emu = (N_MET_And_SelectionCriteria_emu) / (N_SelectionCriteria_emu + 1.0e-06);

	float Alpha_ee = (Eff_Dilepton_Selection_ee * Eff_MET_Selection_ee) / (Eff_MET_Dilepton_Selection_ee + 1.0e-06);
	float Alpha_mumu = (Eff_Dilepton_Selection_mumu * Eff_MET_Selection_mumu) / (Eff_MET_Dilepton_Selection_mumu + 1.0e-06);
	float Alpha_emu = (Eff_Dilepton_Selection_emu * Eff_MET_Selection_emu) / (Eff_MET_Dilepton_Selection_emu + 1.0e-06);


	std::string TriggerSF_Alpha;

        if(process == "Data_triggerSF"){ TriggerSF_Alpha = "TriggerSF_Alpha_Data_MET_" + year + ".txt";}
        else if(process == "MC_triggerSF_ttbar"){ TriggerSF_Alpha = "TriggerSF_Alpha_MC_ttbar_" + year + ".txt";}
	else if(process == "MC_triggerSF_ZPlusJets"){ TriggerSF_Alpha = "TriggerSF_Alpha_MC_ZPlusJets_" + year + ".txt";}
        else{std::cout << "process must be either Data_triggerSF, MC_triggerSF_ttbar or MC_triggerSF_ZPlusJets" << std::endl;}

       std::ofstream TriggerSF_Alpha_File;
        TriggerSF_Alpha_File.open(TriggerSF_Alpha.c_str());

        TriggerSF_Alpha_File << Alpha_ee << '\n'
                             << Alpha_mumu << '\n'
                             << Alpha_emu << '\n' << std::endl;



	//Uncertainties for the efficiencies
	double level = 0.60;

	float Eff_ee_UpperUncert = Eff_ee - TEfficiency::ClopperPearson(N_MET_And_SelectionCriteria_ee, N_MET_DileptonTriggers_And_SelectionCriteria_ee, level, true);
	float Eff_ee_LowerUncert = Eff_ee - TEfficiency::ClopperPearson(N_MET_And_SelectionCriteria_ee, N_MET_DileptonTriggers_And_SelectionCriteria_ee, level, false);
	float Eff_mumu_UpperUncert = Eff_mumu - TEfficiency::ClopperPearson(N_MET_And_SelectionCriteria_mumu, N_MET_DileptonTriggers_And_SelectionCriteria_mumu, level, true);
        float Eff_mumu_LowerUncert = Eff_mumu - TEfficiency::ClopperPearson(N_MET_And_SelectionCriteria_mumu, N_MET_DileptonTriggers_And_SelectionCriteria_mumu,level, false);
	float Eff_emu_UpperUncert = Eff_emu - TEfficiency::ClopperPearson(N_MET_And_SelectionCriteria_emu, N_MET_DileptonTriggers_And_SelectionCriteria_emu,level, true);
        float Eff_emu_LowerUncert = Eff_emu - TEfficiency::ClopperPearson(N_MET_And_SelectionCriteria_emu, N_MET_DileptonTriggers_And_SelectionCriteria_emu,level, false);



	std::string TriggerSF_EfficiencyUncerts;

        if(process == "Data_triggerSF"){ TriggerSF_EfficiencyUncerts = "TriggerSF_EfficiencyUncerts_Data_MET_" + year + ".txt";}
        else if(process == "MC_triggerSF_ttbar"){ TriggerSF_EfficiencyUncerts = "TriggerSF_EfficiencyUncerts_MC_ttbar_" + year + ".txt";}
	else if(process == "MC_triggerSF_ZPlusJets"){ TriggerSF_EfficiencyUncerts = "TriggerSF_EfficiencyUncerts_MC_ZPlusJets_" + year + ".txt";}
        else{std::cout << "process must be either Data_triggerSF, MC_triggerSF_ttbar, MC_triggerSF_ZPlusJets" << std::endl;}

       std::ofstream TriggerSF_EfficiencyUncerts_File;
        TriggerSF_EfficiencyUncerts_File.open(TriggerSF_EfficiencyUncerts.c_str());

        TriggerSF_EfficiencyUncerts_File << Eff_ee_UpperUncert << '\n'
                             		 << Eff_ee_LowerUncert << '\n'
					 << Eff_mumu_UpperUncert << '\n'
                                         << Eff_mumu_LowerUncert << '\n'
					 << Eff_emu_UpperUncert << '\n'
                                         << Eff_emu_LowerUncert << '\n' << std::endl;


	//Turn on curves
	int NumBins = 40;

	//Weights for the turn on curves
	float Weight_ee = N_MET_DileptonTriggers_And_SelectionCriteria_ee / N_SelectionCriteria_ee;
	float Weight_mumu = N_MET_DileptonTriggers_And_SelectionCriteria_mumu / N_SelectionCriteria_mumu;
	float Weight_emu = N_MET_DileptonTriggers_And_SelectionCriteria_emu / N_SelectionCriteria_emu;

	auto TurnOnCurveWeight_ee{[&Weight_ee](const float& PU){float weight = PU * Weight_ee; return weight;}};
	auto TurnOnCurveWeight_mumu{[&Weight_mumu](const float& PU){float weight = PU * Weight_mumu; return weight;}};
	auto TurnOnCurveWeight_emu{[&Weight_emu](const float& PU){float weight = PU * Weight_emu; return weight;}};

	auto Weighted_ee_dataframe = d_ee_selection.Define("ee_weight", TurnOnCurveWeight_ee, {"PU"});
	auto Weighted_mumu_dataframe = d_mumu_selection.Define("mumu_weight", TurnOnCurveWeight_mumu, {"PU"});	
	auto Weighted_emu_dataframe = d_emu_selection.Define("emu_weight", TurnOnCurveWeight_emu, {"PU"});

        //ee
	auto h_LeadingElectronPt_ee = Weighted_ee_dataframe.Profile1D({"h_LeadingElectronPt_ee", "Leading electron p_{T} (ee)", NumBins, 0, 300}, "LeadingElectron_pT", "ee_weight");
	auto h_SubleadingElectronPt_ee = Weighted_ee_dataframe.Profile1D({"h_SubleadingElectronPt_ee", "Subleading electron p_{T} (ee)", NumBins, 0, 300}, "SubleadingElectron_pT", "ee_weight");
	auto h_LeadingElectronEta_ee = Weighted_ee_dataframe.Profile1D({"h_LeadingElectronEta_ee", "Leading electron #eta (ee)", NumBins, -5, 5}, "LeadingElectronEta", "ee_weight");
        auto h_SubleadingElectronEta_ee = Weighted_ee_dataframe.Profile1D({"h_SubleadingElectronEta_ee", "Subleading electron #eta (ee)", NumBins, -5, 5}, "SubleadingElectronEta", "ee_weight");

	//mumu
	auto h_LeadingMuonPt_mumu = Weighted_mumu_dataframe.Profile1D({"h_LeadingMuonPt_mumu", "Leading muon p_{T} (mumu)", NumBins, 0, 300}, "LeadingMuon_pT", "mumu_weight");
        auto h_SubleadingMuonPt_mumu = Weighted_mumu_dataframe.Profile1D({"h_SubleadingMuonPt_ee", "Subleading muon p_{T} (mumu)", NumBins, 0, 300}, "SubleadingMuon_pT", "mumu_weight");
        auto h_LeadingMuonEta_mumu = Weighted_mumu_dataframe.Profile1D({"h_LeadingMuonEta_ee", "Leading muon #eta (ee)", NumBins, -5, 5}, "LeadingMuonEta", "mumu_weight");
        auto h_SubleadingMuonEta_mumu = Weighted_mumu_dataframe.Profile1D({"h_SubleadingMuonEta_ee", "Subleading muon #eta (ee)", NumBins, -5, 5}, "SubleadingMuonEta", "mumu_weight");


	//emu
	auto h_LeadingElectronMuonPt_emu = Weighted_emu_dataframe.Profile1D({"h_LeadingElectronMuonPt_emu", "Leading electron-muon p_{T} (emu)", NumBins, 0, 300}, "LeadingElectronMuon_pT", "emu_weight");
        auto h_SubleadingElectronMuonPt_emu = Weighted_emu_dataframe.Profile1D({"h_SubleadingElectronMuonPt_emu", "Subleading electon-muon p_{T} (emu)", NumBins, 0, 300}, "SubleadingElectronMuon_pT", "emu_weight");
        auto h_LeadingElectronMuonEta_emu = Weighted_emu_dataframe.Profile1D({"h_LeadingElectronMuonEta_emu", "Leading electron-muon #eta (emu)", NumBins, -5, 5}, "LeadingElectronMuonEta", "emu_weight");
        auto h_SubleadingElectronMuonEta_emu = Weighted_emu_dataframe.Profile1D({"h_SubleadingElectronMuonEta_emu", "Subleading electron-muon #eta (emu)", NumBins, -5, 5}, "SubleadingElectronMuonEta", "emu_weight");


	//2D histograms
	//ee
	auto h_LeadingVsSubleading_ElectronPt_ee = Weighted_ee_dataframe.Profile2D({"h_LeadingVsSubleading_ElectronPt_ee", "Leading electron p_{T} vs subleading electron p_{T} (ee)", NumBins, 0, 300, NumBins, 0, 300}, "SubleadingElectron_pT", "LeadingElectron_pT", "ee_weight");

	auto h_LeadingVsSubleading_ElectronEta_ee = Weighted_ee_dataframe.Profile2D({"h_LeadingVsSubleading_ElectronEta_ee", "Leading electron #eta vs subleading electron #eta (ee)", NumBins, -5, 5, NumBins, -5, 5}, "SubleadingElectronEta", "LeadingElectronEta", "ee_weight");
	
	//mumu
	auto h_LeadingVsSubleading_MuonPt_mumu = Weighted_mumu_dataframe.Profile2D({"h_LeadingVsSubleading_MuonPt_mumu", "Leading muon p_{T} vs subleading muon p_{T} (mumu)", NumBins, 0, 300, NumBins, 0, 300}, "SubleadingMuon_pT", "LeadingMuon_pT", "mumu_weight");

        auto h_LeadingVsSubleading_MuonEta_mumu = Weighted_mumu_dataframe.Profile2D({"h_LeadingVsSubleading_MuonEta_mumu", "Leading muon #eta vs subleading muon #eta (mumu)", NumBins, -5, 5, NumBins, -5, 5}, "SubleadingMuonEta", "LeadingMuonEta", "mumu_weight");
	
	//emu	
	auto h_LeadingVsSubleading_ElectronMuonPt_emu = Weighted_emu_dataframe.Profile2D({"h_LeadingVsSubleading_ElectronMuonPt_emu", "Leading electron-muon p_{T} vs subleading electron-muon p_{T} (mumu)", NumBins, 0, 300, NumBins, 0, 300}, "SubleadingElectronMuon_pT", "LeadingElectronMuon_pT", "emu_weight");
        
        auto h_LeadingVsSubleading_ElectronMuonEta_emu = Weighted_emu_dataframe.Profile2D({"h_LeadingVsSubleading_ElectronMuonEta_emu", "Leading electron-muon #eta vs subleading electron-muon #eta (emu)", NumBins, -5, 5, NumBins, -5, 5}, "SubleadingElectronMuonEta", "LeadingElectronMuonEta", "emu_weight");


	//Errors for the turn on curves
	for (Int_t bin = 1; bin != NumBins + 1; bin++){

		double errUp, errDown, error;

		//ee
		//Leading lepton pt
		errUp = (h_LeadingElectronPt_ee->GetBinContent(bin) - TEfficiency::ClopperPearson(h_LeadingElectronPt_ee->GetBinEntries(bin), h_LeadingElectronPt_ee->GetBinEntries(bin) * h_LeadingElectronPt_ee->GetBinContent(bin), level, true));

		errDown = (h_LeadingElectronPt_ee->GetBinContent(bin) - TEfficiency::ClopperPearson(h_LeadingElectronPt_ee->GetBinEntries(bin), h_LeadingElectronPt_ee->GetBinEntries(bin) * h_LeadingElectronPt_ee->GetBinContent(bin), level, false));

		error = errUp > errDown ? errUp : errDown;

		h_LeadingElectronPt_ee->SetBinError(bin, error);
	


		//Subleading lepton pt
                errUp = (h_SubleadingElectronPt_ee->GetBinContent(bin) - TEfficiency::ClopperPearson(h_SubleadingElectronPt_ee->GetBinEntries(bin), h_SubleadingElectronPt_ee->GetBinEntries(bin) * h_SubleadingElectronPt_ee->GetBinContent(bin), level, true));

                errDown = (h_SubleadingElectronPt_ee->GetBinContent(bin) - TEfficiency::ClopperPearson(h_SubleadingElectronPt_ee->GetBinEntries(bin), h_SubleadingElectronPt_ee->GetBinEntries(bin) * h_SubleadingElectronPt_ee->GetBinContent(bin), level, false));

                error = errUp > errDown ? errUp : errDown;

                h_SubleadingElectronPt_ee->SetBinError(bin, error);


		//Leading lepton eta
                errUp = (h_LeadingElectronEta_ee->GetBinContent(bin) - TEfficiency::ClopperPearson(h_LeadingElectronEta_ee->GetBinEntries(bin), h_LeadingElectronEta_ee->GetBinEntries(bin) * h_LeadingElectronEta_ee->GetBinContent(bin), level, true));

                errDown = (h_LeadingElectronEta_ee->GetBinContent(bin) - TEfficiency::ClopperPearson(h_LeadingElectronEta_ee->GetBinEntries(bin), h_LeadingElectronEta_ee->GetBinEntries(bin) * h_LeadingElectronEta_ee->GetBinContent(bin), level, false));

                error = errUp > errDown ? errUp : errDown;

                h_LeadingElectronEta_ee->SetBinError(bin, error);



                //Subleading lepton eta
                errUp = (h_SubleadingElectronEta_ee->GetBinContent(bin) - TEfficiency::ClopperPearson(h_SubleadingElectronEta_ee->GetBinEntries(bin), h_SubleadingElectronEta_ee->GetBinEntries(bin) * h_SubleadingElectronEta_ee->GetBinContent(bin), level, true));

                errDown = (h_SubleadingElectronEta_ee->GetBinContent(bin) - TEfficiency::ClopperPearson(h_SubleadingElectronEta_ee->GetBinEntries(bin), h_SubleadingElectronEta_ee->GetBinEntries(bin) * h_SubleadingElectronEta_ee->GetBinContent(bin), level, false));

                error = errUp > errDown ? errUp : errDown;

                h_SubleadingElectronEta_ee->SetBinError(bin, error);	




		//mumu
		//Leading lepton pt
                errUp = (h_LeadingMuonPt_mumu->GetBinContent(bin) - TEfficiency::ClopperPearson(h_LeadingMuonPt_mumu->GetBinEntries(bin), h_LeadingMuonPt_mumu->GetBinEntries(bin) * h_LeadingMuonPt_mumu->GetBinContent(bin), level, true));

                errDown = (h_LeadingMuonPt_mumu->GetBinContent(bin) - TEfficiency::ClopperPearson(h_LeadingMuonPt_mumu->GetBinEntries(bin), h_LeadingMuonPt_mumu->GetBinEntries(bin) * h_LeadingMuonPt_mumu->GetBinContent(bin), level, false));

                error = errUp > errDown ? errUp : errDown;

                h_LeadingMuonPt_mumu->SetBinError(bin, error);



                //Subleading lepton pt
                errUp = (h_SubleadingMuonPt_mumu->GetBinContent(bin) - TEfficiency::ClopperPearson(h_SubleadingMuonPt_mumu->GetBinEntries(bin), h_SubleadingMuonPt_mumu->GetBinEntries(bin) * h_SubleadingMuonPt_mumu->GetBinContent(bin), level, true));

                errDown = (h_SubleadingMuonPt_mumu->GetBinContent(bin) - TEfficiency::ClopperPearson(h_SubleadingMuonPt_mumu->GetBinEntries(bin), h_SubleadingMuonPt_mumu->GetBinEntries(bin) * h_SubleadingMuonPt_mumu->GetBinContent(bin), level, false));

                error = errUp > errDown ? errUp : errDown;

                h_SubleadingMuonPt_mumu->SetBinError(bin, error);


                //Leading lepton eta
                errUp = (h_LeadingMuonEta_mumu->GetBinContent(bin) - TEfficiency::ClopperPearson(h_LeadingMuonEta_mumu->GetBinEntries(bin), h_LeadingMuonEta_mumu->GetBinEntries(bin) * h_LeadingMuonEta_mumu->GetBinContent(bin), level, true));

                errDown = (h_LeadingMuonEta_mumu->GetBinContent(bin) - TEfficiency::ClopperPearson(h_LeadingMuonEta_mumu->GetBinEntries(bin), h_LeadingMuonEta_mumu->GetBinEntries(bin) * h_LeadingMuonEta_mumu->GetBinContent(bin), level, false));

                error = errUp > errDown ? errUp : errDown;

                h_LeadingMuonEta_mumu->SetBinError(bin, error);



                //Subleading lepton eta
                errUp = (h_SubleadingMuonEta_mumu->GetBinContent(bin) - TEfficiency::ClopperPearson(h_SubleadingMuonEta_mumu->GetBinEntries(bin), h_SubleadingMuonEta_mumu->GetBinEntries(bin) * h_SubleadingMuonEta_mumu->GetBinContent(bin), level, true));

                errDown = (h_SubleadingMuonEta_mumu->GetBinContent(bin) - TEfficiency::ClopperPearson(h_SubleadingMuonEta_mumu->GetBinEntries(bin), h_SubleadingMuonEta_mumu->GetBinEntries(bin) * h_SubleadingMuonEta_mumu->GetBinContent(bin), level, false));

                error = errUp > errDown ? errUp : errDown;

                h_SubleadingMuonEta_mumu->SetBinError(bin, error);		


		//emu
		//Leading lepton pt
                errUp = (h_LeadingElectronMuonPt_emu->GetBinContent(bin) - TEfficiency::ClopperPearson(h_LeadingElectronMuonPt_emu->GetBinEntries(bin), h_LeadingElectronMuonPt_emu->GetBinEntries(bin) * h_LeadingElectronMuonPt_emu->GetBinContent(bin), level, true));

                errDown = (h_LeadingElectronMuonPt_emu->GetBinContent(bin) - TEfficiency::ClopperPearson(h_LeadingElectronMuonPt_emu->GetBinEntries(bin), h_LeadingElectronMuonPt_emu->GetBinEntries(bin) * h_LeadingElectronMuonPt_emu->GetBinContent(bin), level, false));

                error = errUp > errDown ? errUp : errDown;

                h_LeadingElectronMuonPt_emu->SetBinError(bin, error);



                //Subleading lepton pt
                errUp = (h_SubleadingElectronMuonPt_emu->GetBinContent(bin) - TEfficiency::ClopperPearson(h_SubleadingElectronMuonPt_emu->GetBinEntries(bin), h_SubleadingElectronMuonPt_emu->GetBinEntries(bin) * h_SubleadingElectronMuonPt_emu->GetBinContent(bin), level, true));

                errDown = (h_SubleadingElectronMuonPt_emu->GetBinContent(bin) - TEfficiency::ClopperPearson(h_SubleadingElectronMuonPt_emu->GetBinEntries(bin), h_SubleadingElectronMuonPt_emu->GetBinEntries(bin) * h_SubleadingElectronMuonPt_emu->GetBinContent(bin), level, false));

                error = errUp > errDown ? errUp : errDown;

                h_SubleadingElectronMuonPt_emu->SetBinError(bin, error);


                //Leading lepton eta
                errUp = (h_LeadingElectronMuonEta_emu->GetBinContent(bin) - TEfficiency::ClopperPearson(h_LeadingElectronMuonEta_emu->GetBinEntries(bin), h_LeadingElectronMuonEta_emu->GetBinEntries(bin) * h_LeadingElectronMuonEta_emu->GetBinContent(bin), level, true));

                errDown = (h_LeadingElectronMuonEta_emu->GetBinContent(bin) - TEfficiency::ClopperPearson(h_LeadingElectronMuonEta_emu->GetBinEntries(bin), h_LeadingElectronMuonEta_emu->GetBinEntries(bin) * h_LeadingElectronMuonEta_emu->GetBinContent(bin), level, false));

                error = errUp > errDown ? errUp : errDown;

                h_LeadingElectronMuonEta_emu->SetBinError(bin, error);



                //Subleading lepton eta
                errUp = (h_SubleadingElectronMuonEta_emu->GetBinContent(bin) - TEfficiency::ClopperPearson(h_SubleadingElectronMuonEta_emu->GetBinEntries(bin), h_SubleadingElectronMuonEta_emu->GetBinEntries(bin) * h_SubleadingElectronMuonEta_emu->GetBinContent(bin), level, true));

                errDown = (h_SubleadingElectronMuonEta_emu->GetBinContent(bin) - TEfficiency::ClopperPearson(h_SubleadingElectronMuonEta_emu->GetBinEntries(bin), h_SubleadingElectronMuonEta_emu->GetBinEntries(bin) * h_SubleadingElectronMuonEta_emu->GetBinContent(bin), level, false));

                error = errUp > errDown ? errUp : errDown;

                h_SubleadingElectronMuonEta_emu->SetBinError(bin, error);		
	

	}

	//Saving the turn on curves to an output file
	std::string TurnOnCurvesOutput;

	if(blinding == false){

        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                TurnOnCurvesOutput = "TurnOnCurves_" + process + "_" + year + "_NPL.root";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                TurnOnCurvesOutput = "TurnOnCurves_" + process + "_" + year + "_ZPlusJetsCR.root";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                TurnOnCurvesOutput = "TurnOnCurves_" + process + "_" + year + "_ttbarCR.root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                TurnOnCurvesOutput = "TurnOnCurves_" + process + "_" + year + "_NPL_ZPlusJetsCR.root";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                TurnOnCurvesOutput = "TurnOnCurves_" + process + "_" + year + "_NPL_ttbarCR.root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
        else{TurnOnCurvesOutput = "TurnOnCurves_" + process + "_" + year + ".root";}

	}
	else{


        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                TurnOnCurvesOutput = "TurnOnCurves_" + process + "_" + year + "_NPL_Blinded.root";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                TurnOnCurvesOutput = "TurnOnCurves_" + process + "_" + year + "_ZPlusJetsCR_Blinded.root";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                TurnOnCurvesOutput = "TurnOnCurves_" + process + "_" + year + "_ttbarCR_Blinded.root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                TurnOnCurvesOutput = "TurnOnCurves_" + process + "_" + year + "_NPL_ZPlusJetsCR_Blinded.root";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                TurnOnCurvesOutput = "TurnOnCurves_" + process + "_" + year + "_NPL_ttbarCR_Blinded.root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
        else{TurnOnCurvesOutput = "TurnOnCurves_" + process + "_" + year + "_Blinded.root";}


	}		


	TFile* TurnOnCurvesFile = new TFile(TurnOnCurvesOutput.c_str(), "RECREATE");
	
	//1D histograms
	h_LeadingElectronPt_ee->GetYaxis()->SetTitle("Number of events");
	h_SubleadingElectronPt_ee->GetYaxis()->SetTitle("Number of events");
	h_LeadingElectronEta_ee->GetYaxis()->SetTitle("Number of events");
	h_SubleadingElectronEta_ee->GetYaxis()->SetTitle("Number of events");
	h_LeadingMuonPt_mumu->GetYaxis()->SetTitle("Number of events");
	h_SubleadingMuonPt_mumu->GetYaxis()->SetTitle("Number of events");
	h_LeadingMuonEta_mumu->GetYaxis()->SetTitle("Number of events");
	h_SubleadingMuonEta_mumu->GetYaxis()->SetTitle("Number of events");
	h_LeadingElectronMuonPt_emu->GetYaxis()->SetTitle("Number of events");
        h_SubleadingElectronMuonPt_emu->GetYaxis()->SetTitle("Number of events");
        h_LeadingElectronMuonEta_emu->GetYaxis()->SetTitle("Number of events");
        h_SubleadingElectronMuonEta_emu->GetYaxis()->SetTitle("Number of events");

	h_LeadingElectronPt_ee->GetXaxis()->SetTitle("p_{T}");
        h_SubleadingElectronPt_ee->GetXaxis()->SetTitle("p_{T}");
        h_LeadingElectronEta_ee->GetXaxis()->SetTitle("#eta");
        h_SubleadingElectronEta_ee->GetXaxis()->SetTitle("#eta");
        h_LeadingMuonPt_mumu->GetXaxis()->SetTitle("p_{T}");
        h_SubleadingMuonPt_mumu->GetXaxis()->SetTitle("p_{T}");
        h_LeadingMuonEta_mumu->GetXaxis()->SetTitle("#eta");
        h_SubleadingMuonEta_mumu->GetXaxis()->SetTitle("#eta");
        h_LeadingElectronMuonPt_emu->GetXaxis()->SetTitle("p_{T}");
        h_SubleadingElectronMuonPt_emu->GetXaxis()->SetTitle("p_{T}");
        h_LeadingElectronMuonEta_emu->GetXaxis()->SetTitle("#eta");
        h_SubleadingElectronMuonEta_emu->GetXaxis()->SetTitle("#eta");

	h_LeadingElectronPt_ee->Write();
        h_SubleadingElectronPt_ee->Write();
        h_LeadingElectronEta_ee->Write();
        h_SubleadingElectronEta_ee->Write();
        h_LeadingMuonPt_mumu->Write();
        h_SubleadingMuonPt_mumu->Write();
        h_LeadingMuonEta_mumu->Write();
        h_SubleadingMuonEta_mumu->Write();
        h_LeadingElectronMuonPt_emu->Write();
        h_SubleadingElectronMuonPt_emu->Write();
        h_LeadingElectronMuonEta_emu->Write();
        h_SubleadingElectronMuonEta_emu->Write();

	//2D histograms
	h_LeadingVsSubleading_ElectronPt_ee->GetYaxis()->SetTitle("Leading lepton p_{T}");
	h_LeadingVsSubleading_ElectronEta_ee->GetYaxis()->SetTitle("Leading lepton #eta");
	h_LeadingVsSubleading_MuonPt_mumu->GetYaxis()->SetTitle("Leading lepton p_{T}");
	h_LeadingVsSubleading_MuonEta_mumu->GetYaxis()->SetTitle("Leading lepton #eta");
	h_LeadingVsSubleading_ElectronMuonPt_emu->GetYaxis()->SetTitle("Leading lepton p_{T}");
	h_LeadingVsSubleading_ElectronMuonEta_emu->GetYaxis()->SetTitle("Leading lepton #eta");

	h_LeadingVsSubleading_ElectronPt_ee->GetXaxis()->SetTitle("Subleading lepton p_{T}");
	h_LeadingVsSubleading_ElectronEta_ee->GetXaxis()->SetTitle("Subleading lepton #eta");
	h_LeadingVsSubleading_MuonPt_mumu->GetXaxis()->SetTitle("Subleading lepton p_{T}");
	h_LeadingVsSubleading_MuonEta_mumu->GetXaxis()->SetTitle("Subleading lepton #eta");
	h_LeadingVsSubleading_ElectronMuonPt_emu->GetXaxis()->SetTitle("Subleading lepton p_{T}");
	h_LeadingVsSubleading_ElectronMuonEta_emu->GetXaxis()->SetTitle("Subleading lepton #eta");

	h_LeadingVsSubleading_ElectronPt_ee->Write();
	h_LeadingVsSubleading_ElectronEta_ee->Write();
	h_LeadingVsSubleading_MuonPt_mumu->Write();
	h_LeadingVsSubleading_MuonEta_mumu->Write();
	h_LeadingVsSubleading_ElectronMuonPt_emu->Write();
	h_LeadingVsSubleading_ElectronMuonEta_emu->Write();

	TurnOnCurvesFile->Close();

	return;

  } //end of if statement for trigger SF processes


  //Reading the efficiency text files to calculate the trigger SFs
  float MC_Efficiency_Central_ee = ( textfilereader2_TriggerSF("MC_Central", year, blinding) ).at(0);
  float MC_Efficiency_Central_mumu = ( textfilereader2_TriggerSF("MC_Central", year, blinding) ).at(1);
  float MC_Efficiency_Central_emu = ( textfilereader2_TriggerSF("MC_Central", year, blinding) ).at(2);


  float Data_Efficiency_Central_ee = ( textfilereader2_TriggerSF("Data_Central", year, blinding) ).at(0);
  float Data_Efficiency_Central_mumu = ( textfilereader2_TriggerSF("Data_Central", year, blinding) ).at(1);
  float Data_Efficiency_Central_emu = ( textfilereader2_TriggerSF("Data_Central", year, blinding) ).at(2);


  float SF_ee = Data_Efficiency_Central_ee / (MC_Efficiency_Central_ee + 1.0e-06);
  float SF_mumu = Data_Efficiency_Central_mumu / (MC_Efficiency_Central_mumu + 1.0e-06);
  float SF_emu = Data_Efficiency_Central_emu / (MC_Efficiency_Central_emu + 1.0e-06);

  float MC_Efficiency_UpperUncert_ee = ( textfilereader2_TriggerSF("MC_Uncert", year, blinding) ).at(0);
  float MC_Efficiency_LowerUncert_ee = ( textfilereader2_TriggerSF("MC_Uncert", year, blinding) ).at(1);
  float MC_Efficiency_UpperUncert_mumu = ( textfilereader2_TriggerSF("MC_Uncert", year, blinding) ).at(2);
  float MC_Efficiency_LowerUncert_mumu = ( textfilereader2_TriggerSF("MC_Uncert", year, blinding) ).at(3);
  float MC_Efficiency_UpperUncert_emu = ( textfilereader2_TriggerSF("MC_Uncert", year, blinding) ).at(4);
  float MC_Efficiency_LowerUncert_emu = ( textfilereader2_TriggerSF("MC_Uncert", year, blinding) ).at(5);

  float Data_Efficiency_UpperUncert_ee = ( textfilereader2_TriggerSF("Data_Uncert", year, blinding) ).at(0);
  float Data_Efficiency_LowerUncert_ee = ( textfilereader2_TriggerSF("Data_Uncert", year, blinding) ).at(1);
  float Data_Efficiency_UpperUncert_mumu = ( textfilereader2_TriggerSF("Data_Uncert", year, blinding) ).at(2);
  float Data_Efficiency_LowerUncert_mumu = ( textfilereader2_TriggerSF("Data_Uncert", year, blinding) ).at(3);
  float Data_Efficiency_UpperUncert_emu = ( textfilereader2_TriggerSF("Data_Uncert", year, blinding) ).at(4);
  float Data_Efficiency_LowerUncert_emu = ( textfilereader2_TriggerSF("Data_Uncert", year, blinding) ).at(5);


  //Uncertainties for SFs
  //ee
  double SF_UpperUncert_ee = ((Data_Efficiency_Central_ee + Data_Efficiency_UpperUncert_ee) / (MC_Efficiency_Central_ee - MC_Efficiency_LowerUncert_ee + 1.0e-06)) - SF_ee;
  double SF_LowerUncert_ee = ((Data_Efficiency_Central_ee + Data_Efficiency_LowerUncert_ee)/ (MC_Efficiency_Central_ee - MC_Efficiency_UpperUncert_ee + 1.0e-06)) - SF_ee;


  double SF_Uncert_ee = 0.0;
  if (SF_UpperUncert_ee > SF_LowerUncert_ee){SF_Uncert_ee = SF_UpperUncert_ee;}
  else{SF_Uncert_ee = SF_LowerUncert_ee;}

  //mumu
  double SF_UpperUncert_mumu = ((Data_Efficiency_Central_mumu + Data_Efficiency_UpperUncert_mumu) / (MC_Efficiency_Central_mumu - MC_Efficiency_LowerUncert_mumu + 1.0e-06)) - SF_mumu;
  double SF_LowerUncert_mumu = ((Data_Efficiency_Central_mumu + Data_Efficiency_LowerUncert_mumu)/ (MC_Efficiency_Central_mumu - MC_Efficiency_UpperUncert_mumu + 1.0e-06)) - SF_mumu;


  double SF_Uncert_mumu = 0.0;
  if (SF_UpperUncert_mumu > SF_LowerUncert_mumu){SF_Uncert_mumu = SF_UpperUncert_mumu;}
  else{SF_Uncert_mumu = SF_LowerUncert_mumu;}


  //emu
  double SF_UpperUncert_emu = ((Data_Efficiency_Central_emu + Data_Efficiency_UpperUncert_emu) / (MC_Efficiency_Central_emu - MC_Efficiency_LowerUncert_emu + 1.0e-06)) - SF_emu;
  double SF_LowerUncert_emu = ((Data_Efficiency_Central_emu + Data_Efficiency_LowerUncert_emu)/ (MC_Efficiency_Central_emu - MC_Efficiency_UpperUncert_emu + 1.0e-06)) - SF_emu;


  double SF_Uncert_emu = 0.0;
  if (SF_UpperUncert_emu > SF_LowerUncert_emu){SF_Uncert_emu = SF_UpperUncert_emu;}
  else{SF_Uncert_emu = SF_LowerUncert_emu;}


  if(PU_ScaleUp == false && PU_ScaleDown == false && BTag_ScaleUp == false && BTag_ScaleDown == false && JetSmearing_ScaleUp == false && JetSmearing_ScaleDown == false && JetResolution_ScaleUp == false && JetResolution_ScaleDown == false && LeptonEfficiencies_ScaleUp == false && LeptonEfficiencies_ScaleDown == false && PDF_ScaleUp == false && PDF_ScaleDown == false && ME_Up == false && ME_Down == false && MET_Up == false && MET_Down == false && isr_up == false && isr_down == false && fsr_up == false && fsr_down == false){

  	std::string TriggerSF_ScaleFactors = "TriggerSF_ScaleFactors_" + year + ".txt";
  	std::ofstream TriggerSF_ScaleFactors_File;
  	TriggerSF_ScaleFactors_File.open(TriggerSF_ScaleFactors.c_str());

  	TriggerSF_ScaleFactors_File << SF_ee << '\n'
			      	    << SF_mumu << '\n'
			      	    << SF_emu << '\n' 
			      	    << SF_Uncert_ee << '\n'
			     	    << SF_Uncert_mumu << '\n'
			      	    << SF_Uncert_emu << '\n' << std::endl;


  }


  //Filtering events with a reconstructed Z boson
  auto d_ee_recoZ_selection = d_ee_selection.Define("OppositeSignNonPrompt", OppositeSignNonPrompt, {"Electron_charge_Selection", "Electron_genPartFlav"})
                                          .Define("OppositeSignPrompt", OppositeSignPrompt, {"Electron_charge_Selection", "Electron_genPartFlav"})
                                          .Define("SameSignNonPrompt", SameSignNonPrompt, {"Electron_charge_Selection", "Electron_genPartFlav"})
                                          .Define("SameSignPrompt", SameSignPrompt, {"Electron_charge_Selection", "Electron_genPartFlav"})
					  .Define("z_lep_eta", "Electron_eta[TightElectrons]")
                                          .Define("z_lep_phi", "Electron_phi[TightElectrons]")
                                          .Define("z_lep_mass", "Electron_mass[TightElectrons]")
                                          .Define("z_lep_pt", "Electron_pt[TightElectrons]")
                                          .Define("z_mass", inv_mass, {"z_lep_pt", "z_lep_eta", "z_lep_phi", "z_lep_mass"})
                                          .Define("RecoZ", RecoZ, RecoZstrings_ee)
                                          .Define("RecoZPt", TLorentzVectorPt, {"RecoZ"})
                                          .Define("RecoZPhi", TLorentzVectorPhi, {"RecoZ"})
                                          .Define("RecoZEta", TLorentzVectorEta, {"RecoZ"})
                                          .Define("dR_ll", deltaRcheck_float, {"LeadingElectronEta", "LeadingElectronPhi", "SubleadingElectronEta", "SubleadingElectronPhi"})
                                          .Define("dPhi_ll", DeltaPhi_floatandfloat, {"LeadingElectronPhi", "SubleadingElectronPhi"})
                                          .Filter(z_mass_cut, {"z_mass"}, "Z mass cut (ee channel)");


  auto d_mumu_recoZ_selection = d_mumu_selection.Define("DummyColumn", DummyColumnFunction, {"Muon_pt"})
						.Define("Muon_genPartIdx_Selection", select<ints>, {"Muon_genPartIdx", "TightMuons"})
                                                .Define("Muon_nTrackerLayers_Selection", select<ints>, {"Muon_nTrackerLayers", "TightMuons"})
                                                .Define("MuonFourMomentum", MuonFourMomentum, {"Muon_pt_Selection", "Muon_eta_Selection", "Muon_phi_Selection", "Muon_mass_Selection"})
                                                .Define("RochCorrVec", RochCorrVec_Function, {"Muon_charge_Selection", "Muon_pt_Selection", "Muon_eta_Selection", "Muon_phi_Selection", "Muon_genPartIdx_Selection", "Muon_nTrackerLayers_Selection"})
                                                .Define("MuonFourMomentum_RochCorr", RochCorrMuon4Mo, {"MuonFourMomentum", "RochCorrVec"})
                                                .Define("MuonPt_RochCorr", TLorentzVectorPt_float, {"MuonFourMomentum_RochCorr"})
                                                .Define("MuonEta_RochCorr", TLorentzVectorEta_float, {"MuonFourMomentum_RochCorr"})
                                                .Define("MuonPhi_RochCorr", TLorentzVectorPhi_float, {"MuonFourMomentum_RochCorr"})
                                                .Define("MuonMass_RochCorr", TLorentzVectorMass_float, {"MuonFourMomentum_RochCorr"})
					        .Define("OppositeSignNonPrompt", OppositeSignNonPrompt, {"Muon_charge_Selection", "Muon_genPartFlav"})
                                                .Define("OppositeSignPrompt", OppositeSignPrompt, {"Muon_charge_Selection", "Muon_genPartFlav"})
                                                .Define("SameSignNonPrompt", SameSignNonPrompt, {"Muon_charge_Selection", "Muon_genPartFlav"})
                                                .Define("SameSignPrompt", SameSignPrompt, {"Muon_charge_Selection", "Muon_genPartFlav"})
					        .Define("z_lep_eta", "MuonEta_RochCorr")
                                                .Define("z_lep_phi", "MuonPhi_RochCorr")
                                                .Define("z_lep_mass", "MuonMass_RochCorr")
                                                .Define("z_lep_pt", "MuonPt_RochCorr")
                                                .Define("z_mass", inv_mass, {"z_lep_pt", "z_lep_eta", "z_lep_phi", "z_lep_mass"})
                                                .Define("RecoZ", RecoZ, RecoZstrings_mumu)
                                                .Define("RecoZPt", TLorentzVectorPt, {"RecoZ"})
                                                .Define("RecoZPhi", TLorentzVectorPhi, {"RecoZ"})
                                                .Define("RecoZEta", TLorentzVectorPt, {"RecoZ"})
                                                .Define("dR_ll", deltaRcheck_float, {"LeadingMuonPhi", "LeadingMuonEta", "SubleadingMuonEta", "SubleadingMuonPhi"})
                                                .Define("dPhi_ll", DeltaPhi_floatandfloat, {"LeadingMuonPhi", "SubleadingMuonPhi"})
                                                .Filter(z_mass_cut, {"z_mass"}, "Z mass cut (mumu channel)");

 

  auto d_ee_recoZ_jets_selection = d_ee_recoZ_selection.Define("sJER_Nominal", SJER_nominal, sJER_sigmaJER_strings)
                      				       .Define("sJER_up", SJER_up, sJER_sigmaJER_strings)
                      				       .Define("sJER_down", SJER_down, sJER_sigmaJER_strings)
                      				       .Define("sigma_JER", sigma_JER, sJER_sigmaJER_strings)
						       .Define("sigma_JER_up", sigma_JER_up, sJER_sigmaJER_strings)
						       .Define("sigma_JER_down", sigma_JER_down, sJER_sigmaJER_strings)
                      				       .Define("cJER", JetSmearingFunction_HybridMethod, JetSmearingStrings)
                      				       .Define("SmearedJet4Momentum", ApplyCJER, ApplyCJER_strings)
                      				       .Define("SmearedJetPt", GetSmearedJetPt, {"SmearedJet4Momentum", "Jet_pt"})
                      				       .Define("SmearedJetPhi", GetSmearedJetPhi, {"SmearedJet4Momentum", "Jet_phi"})
                      				       .Define("SmearedJetEta", GetSmearedJetEta, {"SmearedJet4Momentum", "Jet_eta"})
                      				       .Define("SmearedJetMass", GetSmearedJetMass, {"SmearedJet4Momentum", "Jet_mass"})
		      				       .Define("LeadingJetMass", LeadingVariable, {JetMassInput})
                                                       .Define("SubleadingJetMass", SubleadingVariable, {JetMassInput})
                                                       .Define("ThirdJetMass", ThirdLeadingVariable, {JetMassInput})
                                                       .Define("FourthJetMass", FourthLeadingVariable, {JetMassInput})
                                                       .Define("LeadingJetpT", LeadingVariable, {JetPtInput})
                                                       .Define("SubleadingJetpT", SubleadingVariable, {JetPtInput})
                                                       .Define("ThirdJetpT", ThirdLeadingVariable, {JetPtInput})
                                                       .Define("FourthJetpT", FourthLeadingVariable, {JetPtInput})
                                                       .Define("SumSquaredPt", SumSquared2LeadingJets_pT, {"LeadingJetpT", "SubleadingJetpT"})
                                                       .Define("JetPtSum", JetPtSum, {"LeadingJetpT", "SubleadingJetpT", "ThirdJetpT", "FourthJetpT"})
                                                       .Define("LeadingJetEta", LeadingVariable, {JetEtaInput})
                                                       .Define("SubleadingJetEta", SubleadingVariable, {JetEtaInput})
                                                       .Define("ThirdJetEta", ThirdLeadingVariable, {JetEtaInput})
                                                       .Define("FourthJetEta", FourthLeadingVariable, {JetEtaInput})
                                                       .Define("LeadingJetPhi", LeadingVariable, {JetPhiInput})
                                                       .Define("SubleadingJetPhi", SubleadingVariable, {JetPhiInput})
                                                       .Define("ThirdJetPhi", ThirdLeadingVariable, {JetPhiInput})
                                                       .Define("FourthJetPhi", FourthLeadingVariable, {JetPhiInput})
                                                       .Define("dRJet_e", deltaRcheck_floats, deltaR_JetE_strings)
                                                       .Define("dR_j1j2", deltaRcheck_float, deltaR_j1j2_strings)
                                                       .Define("dPhi_j1j2", DeltaPhi_floatandfloat, {"LeadingJetPhi", "SubleadingJetPhi"})
                                                       .Define("LeadingJetHT", HT, {"LeadingJetpT"})
                                                       .Define("SubleadingJetHT", HT, {"SubleadingJetpT"})
                                                       .Define("ThirdJetHT", HT, {"ThirdJetpT"})
                                                       .Define("FourthJetHT", HT, {"FourthJetpT"})
                                                       .Define("TotJetHT", TotJetHT, {"LeadingJetHT", "SubleadingJetHT", "ThirdJetHT", "FourthJetHT"})
                                                       .Define("LeadingElectronHT", HT, {"LeadingElectron_pT"})
                                                       .Define("SubleadingElectronHT", HT, {"SubleadingElectron_pT"})
                                                       .Define("TotLepHT", TotLepHT, {"LeadingElectronHT", "SubleadingElectronHT"})
                                                       .Define("TotHTOverTotpT_Jets", TotHTOverTotpT, {"TotJetHT", "JetPtSum"})
                                                       .Define("LepPtSum", LepPtSum, {"LeadingElectron_pT", "SubleadingElectron_pT"})
                                                       .Define("LepEtaSum", LepEtaSum, {"LeadingElectron_pT", "SubleadingElectron_pT"})
                                                       .Define("LepPhiSum", LepPhiSum, {"LeadingElectron_pT", "SubleadingElectron_pT"})
                                                       .Define("TotHTOverTotpT_Leptons", TotHTOverTotpT, {"TotLepHT", "LepPtSum"})
                                                       .Define("InvMassAllJets", InvMass_AllJets, InvMass_AllJets_strings)
                                                       .Define("InvMass3Jets", InvMass_3Jets, InvMass_3Jets_strings)
                                                       .Define("JetEtaSum", JetEtaSum, {"LeadingJetEta", "SubleadingJetEta", "ThirdJetEta", "FourthJetEta"})
                                                       .Define("JetPhiSum", JetPhiSum, {"LeadingJetPhi", "SubleadingJetPhi", "ThirdJetPhi", "FourthJetPhi"})
                                                       .Define("tight_jets", tight_jets_function, {JetPtInput, JetEtaInput, "Jet_jetId", "dRJet_e"})
                                                       .Filter(jet_selection_function, {"tight_jets"}, "jet cut (ee channel)");


  auto d_ee_recoZ_jets_bjets_selection = d_ee_recoZ_jets_selection.Define("bjets", bjet_id, {"tight_jets", "Jet_btagCSVV2", JetEtaInput})
                                                                  .Define("nbjets", numberofbjets, {"bjets"})
                                                                  .Define("BTAGEFF_bjet_id_WP", BTAGEFF_bjet_id_WP, {"tight_jets", "Jet_btagCSVV2", JetEtaInput, "Jet_partonFlavour"})
								  .Define("BTAGEFF_nonbjet_id_WP", BTAGEFF_nonbjet_id_WP, {"tight_jets", "Jet_btagCSVV2", JetEtaInput, "Jet_partonFlavour"})
                                                                  .Define("BTAGEFF_charm_id_WP", BTAGEFF_charm_id_WP, {"tight_jets", "Jet_btagCSVV2", JetEtaInput, "Jet_partonFlavour"})
                                                                  .Define("BTAGEFF_lightjets_id_WP", BTAGEFF_lightjets_id_WP, {"tight_jets", "Jet_btagCSVV2", JetEtaInput, "Jet_partonFlavour"})
                                                                  .Define("BTAGEFF_gluon_id_WP", BTAGEFF_gluon_id_WP, {"tight_jets", "Jet_btagCSVV2", JetEtaInput, "Jet_partonFlavour"})
                                                                  .Define("BTAGEFF_bjet_id", BTAGEFF_bjet_id, {"tight_jets", JetEtaInput, "Jet_partonFlavour"})
								  .Define("BTAGEFF_nonbjet_id", BTAGEFF_nonbjet_id, {"tight_jets", JetEtaInput, "Jet_partonFlavour"})
                                                                  .Define("BTAGEFF_charm_id", BTAGEFF_charm_id, {"tight_jets", JetEtaInput, "Jet_partonFlavour"})
                                                                  .Define("BTAGEFF_lightjets_id", BTAGEFF_lightjets_id, {"tight_jets", JetEtaInput, "Jet_partonFlavour"})
                                                                  .Define("BTAGEFF_gluon_id", BTAGEFF_gluon_id, {"tight_jets", JetEtaInput, "Jet_partonFlavour"})
                                                                  .Define("BTAGEFF_bjet_pt_num", select<floats>, {JetPtInput, "BTAGEFF_bjet_id_WP"})
                                                                  .Define("BTAGEFF_bjet_eta_num", select<floats>, {JetEtaInput, "BTAGEFF_bjet_id_WP"})
								  .Define("BTAGEFF_nonbjet_pt_num", select<floats>, {JetPtInput, "BTAGEFF_nonbjet_id_WP"})
                                                                  .Define("BTAGEFF_nonbjet_eta_num", select<floats>, {JetEtaInput, "BTAGEFF_nonbjet_id_WP"})
                                                                  .Define("BTAGEFF_charm_pt_num", select<floats>, {JetPtInput, "BTAGEFF_charm_id_WP"})
                                                                  .Define("BTAGEFF_charm_eta_num", select<floats>, {JetEtaInput, "BTAGEFF_charm_id_WP"})
                                                                  .Define("BTAGEFF_lightjets_pt_num", select<floats>, {JetPtInput, "BTAGEFF_lightjets_id_WP"})
                                                                  .Define("BTAGEFF_lightjets_eta_num", select<floats>, {JetEtaInput, "BTAGEFF_lightjets_id_WP"})
                                                                  .Define("BTAGEFF_gluon_pt_num", select<floats>, {JetPtInput, "BTAGEFF_gluon_id_WP"})
                                                                  .Define("BTAGEFF_gluon_eta_num", select<floats>, {JetEtaInput, "BTAGEFF_gluon_id_WP"})
                                                                  .Define("BTAGEFF_bjet_pt_denom", select<floats>, {JetPtInput, "BTAGEFF_bjet_id"})
                                                                  .Define("BTAGEFF_bjet_eta_denom", select<floats>, {JetEtaInput, "BTAGEFF_bjet_id"})
								  .Define("BTAGEFF_nonbjet_pt_denom", select<floats>, {JetPtInput, "BTAGEFF_nonbjet_id"})
                                                                  .Define("BTAGEFF_nonbjet_eta_denom", select<floats>, {JetEtaInput, "BTAGEFF_nonbjet_id"})
                                                                  .Define("BTAGEFF_charm_pt_denom", select<floats>, {JetPtInput, "BTAGEFF_charm_id"})
                                                                  .Define("BTAGEFF_charm_eta_denom", select<floats>, {JetEtaInput, "BTAGEFF_charm_id"})
                                                                  .Define("BTAGEFF_lightjets_pt_denom", select<floats>, {JetPtInput, "BTAGEFF_lightjets_id"})
                                                                  .Define("BTAGEFF_lightjets_eta_denom", select<floats>, {JetEtaInput, "BTAGEFF_lightjets_id"})
                                                                  .Define("BTAGEFF_gluon_pt_denom", select<floats>, {JetPtInput, "BTAGEFF_gluon_id"})
                                                                  .Define("BTAGEFF_gluon_eta_denom", select<floats>, {JetEtaInput, "BTAGEFF_gluon_id"})
								  .Filter(bjet_cut, {"bjets"}, "b jet cut (ee channel)");
								
  auto d_mumu_recoZ_jets_selection = d_mumu_recoZ_selection.Define("sJER_Nominal", SJER_nominal, sJER_sigmaJER_strings)
                      					   .Define("sJER_up", SJER_up, sJER_sigmaJER_strings)
                      					   .Define("sJER_down", SJER_down, sJER_sigmaJER_strings)
                      					   .Define("sigma_JER", sigma_JER, sJER_sigmaJER_strings)
							   .Define("sigma_JER_up", sigma_JER_up, sJER_sigmaJER_strings)
                                                           .Define("sigma_JER_down", sigma_JER_down, sJER_sigmaJER_strings)
                      					   .Define("cJER", JetSmearingFunction_HybridMethod, JetSmearingStrings)
                      					   .Define("SmearedJet4Momentum", ApplyCJER, ApplyCJER_strings)
                      					   .Define("SmearedJetPt", GetSmearedJetPt, {"SmearedJet4Momentum", "Jet_pt"})
                      					   .Define("SmearedJetPhi", GetSmearedJetPhi, {"SmearedJet4Momentum", "Jet_phi"})
                      					   .Define("SmearedJetEta", GetSmearedJetEta, {"SmearedJet4Momentum", "Jet_eta"})
                      					   .Define("SmearedJetMass", GetSmearedJetMass, {"SmearedJet4Momentum", "Jet_mass"})
		      					   .Define("LeadingJetMass", LeadingVariable, {JetMassInput})
                                                           .Define("SubleadingJetMass", SubleadingVariable, {JetMassInput})
                                                           .Define("ThirdJetMass", ThirdLeadingVariable, {JetMassInput})
                                                           .Define("FourthJetMass", FourthLeadingVariable, {JetMassInput})
                                                           .Define("LeadingJetpT", LeadingVariable, {JetPtInput})
                                                           .Define("SubleadingJetpT", SubleadingVariable, {JetPtInput})
                                                           .Define("ThirdJetpT", ThirdLeadingVariable, {JetPtInput})
                                                           .Define("FourthJetpT", FourthLeadingVariable, {JetPtInput})
                                                           .Define("SumSquaredPt", SumSquared2LeadingJets_pT, {"LeadingJetpT", "SubleadingJetpT"})
                                                           .Define("JetPtSum", JetPtSum, {"LeadingJetpT", "SubleadingJetpT", "ThirdJetpT", "FourthJetpT"})
                                                           .Define("LeadingJetEta", LeadingVariable, {JetEtaInput})
                                                           .Define("SubleadingJetEta", SubleadingVariable, {JetEtaInput})
                                                           .Define("ThirdJetEta", ThirdLeadingVariable, {JetEtaInput})
                                                           .Define("FourthJetEta", FourthLeadingVariable, {JetEtaInput})
                                                           .Define("LeadingJetPhi", LeadingVariable, {JetPhiInput})
                                                           .Define("SubleadingJetPhi", SubleadingVariable, {JetPhiInput})
                                                           .Define("ThirdJetPhi", ThirdLeadingVariable, {JetPhiInput})
                                                           .Define("FourthJetPhi", FourthLeadingVariable, {JetPhiInput})
                                                           .Define("dRJet_mu", deltaRcheck_floats, deltaR_JetMu_strings)
                                                           .Define("dR_j1j2", deltaRcheck_float, deltaR_j1j2_strings)
                                                           .Define("dPhi_j1j2", DeltaPhi_floatandfloat, {"LeadingJetPhi", "SubleadingJetPhi"})
                                                           .Define("LeadingJetHT", HT, {"LeadingJetpT"})
                                                           .Define("SubleadingJetHT", HT, {"SubleadingJetpT"})
                                                           .Define("ThirdJetHT", HT, {"ThirdJetpT"})
                                                           .Define("FourthJetHT", HT, {"FourthJetpT"})
                                                           .Define("TotJetHT", TotJetHT, {"LeadingJetHT", "SubleadingJetHT", "ThirdJetHT", "FourthJetHT"})
                                                           .Define("LeadingMuonHT", HT, {"LeadingMuon_pT"})
                                                           .Define("SubleadingMuonHT", HT, {"SubleadingMuon_pT"})
                                                           .Define("TotLepHT", TotLepHT, {"LeadingMuonHT", "SubleadingMuonHT"})
                                                           .Define("TotHTOverTotpT_Jets", TotHTOverTotpT, {"TotJetHT", "JetPtSum"})
                                                           .Define("LepPtSum", LepPtSum, {"LeadingMuon_pT", "SubleadingMuon_pT"})
                                                           .Define("LepEtaSum", LepEtaSum, {"LeadingMuonEta", "SubleadingMuonEta"})
                                                           .Define("LepPhiSum", LepPhiSum, {"LeadingMuonPhi", "SubleadingMuonPhi"})
                                                           .Define("TotHTOverTotpT_Leptons", TotHTOverTotpT, {"TotLepHT", "LepPtSum"})
                                                           .Define("InvMassAllJets", InvMass_AllJets, InvMass_AllJets_strings)
                                                           .Define("InvMass3Jets", InvMass_3Jets, InvMass_3Jets_strings)
                                                           .Define("JetEtaSum", JetEtaSum, {"LeadingJetEta", "SubleadingJetEta", "ThirdJetEta", "FourthJetEta"})
                                                           .Define("JetPhiSum", JetPhiSum, {"LeadingJetPhi", "SubleadingJetPhi", "ThirdJetPhi", "FourthJetPhi"})
                                                           .Define("tight_jets", tight_jets_function, {JetPtInput, JetEtaInput, "Jet_jetId", "dRJet_mu"})
                                                           .Filter(jet_selection_function, {"tight_jets"}, "jet cut (mumu channel)");

   auto d_mumu_recoZ_jets_bjets_selection = d_mumu_recoZ_jets_selection.Define("bjets", bjet_id, {"tight_jets", "Jet_btagCSVV2", JetEtaInput})
                                                                      .Define("nbjets", numberofbjets, {"bjets"})
                                                                      .Define("BTAGEFF_bjet_id_WP", BTAGEFF_bjet_id_WP, {"tight_jets", "Jet_btagCSVV2", JetEtaInput, "Jet_partonFlavour"})
                                                                      .Define("BTAGEFF_charm_id_WP", BTAGEFF_charm_id_WP, {"tight_jets", "Jet_btagCSVV2", JetEtaInput, "Jet_partonFlavour"})
                                                                      .Define("BTAGEFF_lightjets_id_WP", BTAGEFF_lightjets_id_WP, {"tight_jets", "Jet_btagCSVV2", JetEtaInput, "Jet_partonFlavour"})
                                                                      .Define("BTAGEFF_gluon_id_WP", BTAGEFF_gluon_id_WP, {"tight_jets", "Jet_btagCSVV2", JetEtaInput, "Jet_partonFlavour"})
								      .Define("BTAGEFF_nonbjet_id_WP", BTAGEFF_nonbjet_id_WP, {"tight_jets", "Jet_btagCSVV2", JetEtaInput, "Jet_partonFlavour"})
                                                                      .Define("BTAGEFF_bjet_id", BTAGEFF_bjet_id, {"tight_jets", JetEtaInput, "Jet_partonFlavour"})
                                                                      .Define("BTAGEFF_charm_id", BTAGEFF_charm_id, {"tight_jets", JetEtaInput, "Jet_partonFlavour"})
                                                                      .Define("BTAGEFF_lightjets_id", BTAGEFF_lightjets_id, {"tight_jets", JetEtaInput, "Jet_partonFlavour"})
                                                                      .Define("BTAGEFF_gluon_id", BTAGEFF_gluon_id, {"tight_jets", JetEtaInput, "Jet_partonFlavour"})
								      .Define("BTAGEFF_nonbjet_id", BTAGEFF_nonbjet_id, {"tight_jets", JetEtaInput, "Jet_partonFlavour"})
                                                                      .Define("BTAGEFF_bjet_pt_num", select<floats>, {JetPtInput, "BTAGEFF_bjet_id_WP"})
                                                                      .Define("BTAGEFF_bjet_eta_num", select<floats>, {JetEtaInput, "BTAGEFF_bjet_id_WP"})
                                                                      .Define("BTAGEFF_charm_pt_num", select<floats>, {JetPtInput, "BTAGEFF_charm_id_WP"})
                                                                      .Define("BTAGEFF_charm_eta_num", select<floats>, {JetEtaInput, "BTAGEFF_charm_id_WP"})
                                                                      .Define("BTAGEFF_lightjets_pt_num", select<floats>, {JetPtInput, "BTAGEFF_lightjets_id_WP"})
                                                                      .Define("BTAGEFF_lightjets_eta_num", select<floats>, {JetEtaInput, "BTAGEFF_lightjets_id_WP"})
                                                                      .Define("BTAGEFF_gluon_pt_num", select<floats>, {JetPtInput, "BTAGEFF_gluon_id_WP"})
                                                                      .Define("BTAGEFF_gluon_eta_num", select<floats>, {JetEtaInput, "BTAGEFF_gluon_id_WP"})
								      .Define("BTAGEFF_nonbjet_pt_num", select<floats>, {JetPtInput, "BTAGEFF_nonbjet_id_WP"})
                                                                      .Define("BTAGEFF_nonbjet_eta_num", select<floats>, {JetEtaInput, "BTAGEFF_nonbjet_id_WP"})
                                                                      .Define("BTAGEFF_bjet_pt_denom", select<floats>, {JetPtInput, "BTAGEFF_bjet_id"})
                                                                      .Define("BTAGEFF_bjet_eta_denom", select<floats>, {JetEtaInput, "BTAGEFF_bjet_id"})
                                                                      .Define("BTAGEFF_charm_pt_denom", select<floats>, {JetPtInput, "BTAGEFF_charm_id"})
                                                                      .Define("BTAGEFF_charm_eta_denom", select<floats>, {JetEtaInput, "BTAGEFF_charm_id"})
                                                                      .Define("BTAGEFF_lightjets_pt_denom", select<floats>, {JetPtInput, "BTAGEFF_lightjets_id"})
                                                                      .Define("BTAGEFF_lightjets_eta_denom", select<floats>, {JetEtaInput, "BTAGEFF_lightjets_id"})
                                                                      .Define("BTAGEFF_gluon_pt_denom", select<floats>, {JetPtInput, "BTAGEFF_gluon_id"})
                                                                      .Define("BTAGEFF_gluon_eta_denom", select<floats>, {JetEtaInput, "BTAGEFF_gluon_id"})
								      .Define("BTAGEFF_nonbjet_pt_denom", select<floats>, {JetPtInput, "BTAGEFF_nonbjet_id"})
                                                                      .Define("BTAGEFF_nonbjet_eta_denom", select<floats>, {JetEtaInput, "BTAGEFF_nonbjet_id"})
                                                                      .Filter(bjet_cut, {"bjets"}, "b jet cut (ee channel)");


std::string BTagEffOutput = process + "_BTagEffPlots.root";


if(blinding == false){

	if(PU_ScaleUp == true){

		if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PU_ScaleUp_" + year + "_NPL.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PU_ScaleUp_" + year + "_ZPlusJetsCR.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PU_ScaleUp_" + year + "_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PU_ScaleUp_" + year + "_NPL_ZPlusJetsCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PU_ScaleUp_" + year + "_NPL_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_PU_ScaleUp_" + year + ".root";}

	}
	else if(PU_ScaleDown == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PU_ScaleDown_" + year + "_NPL.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PU_ScaleDown_" + year + "_ZPlusJetsCR.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PU_ScaleDown_" + year + "_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PU_ScaleDown_" + year + "_NPL_ZPlusJetsCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PU_ScaleDown_" + year + "_NPL_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_PU_ScaleDown_" + year + ".root";}

        }
	else if(BTag_ScaleUp == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_BTag_ScaleUp_" + year + "_NPL.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_BTag_ScaleUp_" + year + "_ZPlusJetsCR.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_BTag_ScaleUp_" + year + "_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_BTag_ScaleUp_" + year + "_NPL_ZPlusJetsCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_BTag_ScaleUp_" + year + "_NPL_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_BTag_ScaleUp_" + year + ".root";}

        }
        else if(BTag_ScaleDown == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_BTag_ScaleDown_" + year + "_NPL.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_BTag_ScaleDown_" + year + "_ZPlusJetsCR.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_BTag_ScaleDown_" + year + "_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_BTag_ScaleDown_" + year + "_NPL_ZPlusJetsCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_BTag_ScaleDown_" + year + "_NPL_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_BTag_ScaleDown_" + year + ".root";}

        }
	else if(JetSmearing_ScaleUp == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetSmearing_ScaleUp_" + year + "_NPL.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetSmearing_ScaleUp_" + year + "_ZPlusJetsCR.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetSmearing_ScaleUp_" + year + "_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetSmearing_ScaleUp_" + year + "_NPL_ZPlusJetsCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetSmearing_ScaleUp_" + year + "_NPL_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_JetSmearing_ScaleUp_" + year + ".root";}

        }
        else if(JetSmearing_ScaleDown == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetSmearing_ScaleDown_" + year + "_NPL.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetSmearing_ScaleDown_" + year + "_ZPlusJetsCR.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetSmearing_ScaleDown_" + year + "_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetSmearing_ScaleDown_" + year + "_NPL_ZPlusJetsCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetSmearing_ScaleDown_" + year + "_NPL_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_JetSmearing_ScaleDown_" + year + ".root";}

        }
	else if(JetResolution_ScaleUp == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetResolution_ScaleUp_" + year + "_NPL.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetResolution_ScaleUp_" + year + "_ZPlusJetsCR.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetResolution_ScaleUp_" + year + "_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetResolution_ScaleUp_" + year + "_NPL_ZPlusJetsCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetResolution_ScaleUp_" + year + "_NPL_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_JetResolution_ScaleUp_" + year + ".root";}

        }
        else if(JetResolution_ScaleDown == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetResolution_ScaleDown_" + year + "_NPL.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetResolution_ScaleDown_" + year + "_ZPlusJetsCR.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetResolution_ScaleDown_" + year + "_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetResolution_ScaleDown_" + year + "_NPL_ZPlusJetsCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetResolution_ScaleDown_" + year + "_NPL_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_JetResolution_ScaleDown_" + year + ".root";}

        }
	else if(LeptonEfficiencies_ScaleUp == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_LeptonEfficiencies_ScaleUp_" + year + "_NPL.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_LeptonEfficiencies_ScaleUp_" + year + "_ZPlusJetsCR.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_LeptonEfficiencies_ScaleUp_" + year + "_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_LeptonEfficiencies_ScaleUp_" + year + "_NPL_ZPlusJetsCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_LeptonEfficiencies_ScaleUp_" + year + "_NPL_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_LeptonEfficiencies_ScaleUp_" + year + ".root";}

        }
        else if(LeptonEfficiencies_ScaleDown == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_LeptonEfficiencies_ScaleDown_" + year + "_NPL.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_LeptonEfficiencies_ScaleDown_" + year + "_ZPlusJetsCR.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_LeptonEfficiencies_ScaleDown_" + year + "_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_LeptonEfficiencies_ScaleDown_" + year + "_NPL_ZPlusJetsCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_LeptonEfficiencies_ScaleDown_" + year + "_NPL_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_LeptonEfficiencies_ScaleDown_" + year + ".root";}

        }
	else if(PDF_ScaleUp == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PDF_ScaleUp_" + year + "_NPL.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PDF_ScaleUp_" + year + "_ZPlusJetsCR.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PDF_ScaleUp_" + year + "_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PDF_ScaleUp_" + year + "_NPL_ZPlusJetsCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PDF_ScaleUp_" + year + "_NPL_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_PDF_ScaleUp_" + year + ".root";}

        }
        else if(PDF_ScaleDown == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PDF_ScaleDown_" + year + "_NPL.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PDF_ScaleDown_" + year + "_ZPlusJetsCR.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PDF_ScaleDown_" + year + "_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PDF_ScaleDown_" + year + "_NPL_ZPlusJetsCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PDF_ScaleDown_" + year + "_NPL_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_PDF_ScaleDown_" + year + ".root";}

        }
	else if(ME_Up == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_ME_Up_" + year + "_NPL.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_ME_Up_" + year + "_ZPlusJetsCR.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_ME_Up_" + year + "_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_ME_Up_" + year + "_NPL_ZPlusJetsCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_ME_Up_" + year + "_NPL_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_ME_Up_" + year + ".root";}

        }
        else if(ME_Down == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_ME_Down_" + year + "_NPL.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_ME_Down_" + year + "_ZPlusJetsCR.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_ME_Down_" + year + "_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_ME_Down_" + year + "_NPL_ZPlusJetsCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_ME_Down_" + year + "_NPL_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_ME_Down_" + year + ".root";}

        }
	else if(MET_Up == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_MET_Up_" + year + "_NPL.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_MET_Up_" + year + "_ZPlusJetsCR.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_MET_Up_" + year + "_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_MET_Up_" + year + "_NPL_ZPlusJetsCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_MET_Up_" + year + "_NPL_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_MET_Up_" + year + ".root";}

        }
        else if(MET_Down == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_MET_Down_" + year + "_NPL.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_MET_Down_" + year + "_ZPlusJetsCR.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_MET_Down_" + year + "_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_MET_Down_" + year + "_NPL_ZPlusJetsCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_MET_Down_" + year + "_NPL_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_MET_Down_" + year + ".root";}

        }
	else if(isr_up == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_isr_up_" + year + "_NPL.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_isr_up_" + year + "_ZPlusJetsCR.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_isr_up_" + year + "_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_isr_up_" + year + "_NPL_ZPlusJetsCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_isr_up_" + year + "_NPL_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_isr_up_" + year + ".root";}

        }
        else if(isr_down == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_isr_down_" + year + "_NPL.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_isr_down_" + year + "_ZPlusJetsCR.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_isr_down_" + year + "_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_isr_down_" + year + "_NPL_ZPlusJetsCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_isr_down_" + year + "_NPL_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_isr_down_" + year + ".root";}

        }
	else if(fsr_up == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_fsr_up_" + year + "_NPL.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_fsr_up_" + year + "_ZPlusJetsCR.root";
                }       
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_fsr_up_" + year + "_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_fsr_up_" + year + "_NPL_ZPlusJetsCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_fsr_up_" + year + "_NPL_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_fsr_up_" + year + ".root";}

        }
        else if(fsr_down == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_fsr_down_" + year + "_NPL.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_fsr_down_" + year + "_ZPlusJetsCR.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_fsr_down_" + year + "_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_fsr_down_" + year + "_NPL_ZPlusJetsCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_fsr_down_" + year + "_NPL_ttbarCR.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_fsr_down_" + year + ".root";}

        }
	else{

        	if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                	BTagEffOutput = "BTagEffPlots_" + process + "_" + year + "_NPL.root";
        	}
        	else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                	BTagEffOutput = "BTagEffPlots_" + process + "_" + year + "_ZPlusJetsCR.root";
        	}
        	else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                	BTagEffOutput = "BTagEffPlots_" + process + "_" + year + "_ttbarCR.root";
        	}
        	else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                	BTagEffOutput = "BTagEffPlots_" + process + "_" + year + "_NPL_ZPlusJetsCR.root";
        	}
        	else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                	BTagEffOutput = "BTagEffPlots_" + process + "_" + year + "_NPL_ttbarCR.root";
        	}
        	else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
        	else{BTagEffOutput = "BTagEffPlots_" + process + "_" + year + ".root";}

	}

}
else{


	if(PU_ScaleUp == true){

		if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PU_ScaleUp_" + year + "_NPL_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PU_ScaleUp_" + year + "_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PU_ScaleUp_" + year + "_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PU_ScaleUp_" + year + "_NPL_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PU_ScaleUp_" + year + "_NPL_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_PU_ScaleUp_" + year + "_Blinded.root";}

	}
	else if(PU_ScaleDown == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PU_ScaleDown_" + year + "_NPL_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PU_ScaleDown_" + year + "_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PU_ScaleDown_" + year + "_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PU_ScaleDown_" + year + "_NPL_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PU_ScaleDown_" + year + "_NPL_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_PU_ScaleDown_" + year + "_Blinded.root";}

        }
	else if(BTag_ScaleUp == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_BTag_ScaleUp_" + year + "_NPL_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_BTag_ScaleUp_" + year + "_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_BTag_ScaleUp_" + year + "_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_BTag_ScaleUp_" + year + "_NPL_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_BTag_ScaleUp_" + year + "_NPL_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_BTag_ScaleUp_" + year + "_Blinded.root";}

        }
        else if(BTag_ScaleDown == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_BTag_ScaleDown_" + year + "_NPL_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_BTag_ScaleDown_" + year + "_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_BTag_ScaleDown_" + year + "_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_BTag_ScaleDown_" + year + "_NPL_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_BTag_ScaleDown_" + year + "_NPL_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_BTag_ScaleDown_" + year + "_Blinded.root";}

        }
	else if(JetSmearing_ScaleUp == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetSmearing_ScaleUp_" + year + "_NPL_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetSmearing_ScaleUp_" + year + "_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetSmearing_ScaleUp_" + year + "_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetSmearing_ScaleUp_" + year + "_NPL_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetSmearing_ScaleUp_" + year + "_NPL_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_JetSmearing_ScaleUp_" + year + "_Blinded.root";}

        }
        else if(JetSmearing_ScaleDown == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetSmearing_ScaleDown_" + year + "_NPL_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetSmearing_ScaleDown_" + year + "_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetSmearing_ScaleDown_" + year + "_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetSmearing_ScaleDown_" + year + "_NPL_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetSmearing_ScaleDown_" + year + "_NPL_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_JetSmearing_ScaleDown_" + year + "_Blinded.root";}

        }
	else if(JetResolution_ScaleUp == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetResolution_ScaleUp_" + year + "_NPL_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetResolution_ScaleUp_" + year + "_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetResolution_ScaleUp_" + year + "_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetResolution_ScaleUp_" + year + "_NPL_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetResolution_ScaleUp_" + year + "_NPL_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_JetResolution_ScaleUp_" + year + "_Blinded.root";}

        }
        else if(JetResolution_ScaleDown == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetResolution_ScaleDown_" + year + "_NPL_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetResolution_ScaleDown_" + year + "_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetResolution_ScaleDown_" + year + "_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetResolution_ScaleDown_" + year + "_NPL_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_JetResolution_ScaleDown_" + year + "_NPL_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_JetResolution_ScaleDown_" + year + "_Blinded.root";}

        }
	else if(LeptonEfficiencies_ScaleUp == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_LeptonEfficiencies_ScaleUp_" + year + "_NPL_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_LeptonEfficiencies_ScaleUp_" + year + "_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_LeptonEfficiencies_ScaleUp_" + year + "_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_LeptonEfficiencies_ScaleUp_" + year + "_NPL_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_LeptonEfficiencies_ScaleUp_" + year + "_NPL_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_LeptonEfficiencies_ScaleUp_" + year + "_Blinded.root";}

        }
        else if(LeptonEfficiencies_ScaleDown == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_LeptonEfficiencies_ScaleDown_" + year + "_NPL_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_LeptonEfficiencies_ScaleDown_" + year + "_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_LeptonEfficiencies_ScaleDown_" + year + "_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_LeptonEfficiencies_ScaleDown_" + year + "_NPL_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_LeptonEfficiencies_ScaleDown_" + year + "_NPL_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_LeptonEfficiencies_ScaleDown_" + year + "_Blinded.root";}

        }
	else if(PDF_ScaleUp == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PDF_ScaleUp_" + year + "_NPL_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PDF_ScaleUp_" + year + "_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PDF_ScaleUp_" + year + "_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PDF_ScaleUp_" + year + "_NPL_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PDF_ScaleUp_" + year + "_NPL_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_PDF_ScaleUp_" + year + "_Blinded.root";}

        }
        else if(PDF_ScaleDown == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PDF_ScaleDown_" + year + "_NPL_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PDF_ScaleDown_" + year + "_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PDF_ScaleDown_" + year + "_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PDF_ScaleDown_" + year + "_NPL_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_PDF_ScaleDown_" + year + "_NPL_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_PDF_ScaleDown_" + year + "_Blinded.root";}

        }
	else if(ME_Up == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_ME_Up_" + year + "_NPL_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_ME_Up_" + year + "_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_ME_Up_" + year + "_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_ME_Up_" + year + "_NPL_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_ME_Up_" + year + "_NPL_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_ME_Up_" + year + "_Blinded.root";}

        }
        else if(ME_Down == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_ME_Down_" + year + "_Blinded_NPL.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_ME_Down_" + year + "_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_ME_Down_" + year + "_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_ME_Down_" + year + "_NPL_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_ME_Down_" + year + "_NPL_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_ME_Down_" + year + "_Blinded.root";}

        }
	else if(MET_Up == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_MET_Up_" + year + "_NPL_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_MET_Up_" + year + "_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_MET_Up_" + year + "_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_MET_Up_" + year + "_NPL_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_MET_Up_" + year + "_NPL_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_MET_Up_" + year + "_Blinded.root";}

        }
        else if(MET_Down == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_MET_Down_" + year + "_Blinded_NPL.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_MET_Down_" + year + "_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_MET_Down_" + year + "_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_MET_Down_" + year + "_NPL_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_MET_Down_" + year + "_NPL_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_MET_Down_" + year + "_Blinded.root";}

        }
	else if(isr_up == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_isr_up_" + year + "_NPL_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_isr_up_" + year + "_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_isr_up_" + year + "_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_isr_up_" + year + "_NPL_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_isr_up_" + year + "_NPL_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_isr_up_" + year + "_Blinded.root";}

        }
        else if(isr_down == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_isr_down_" + year + "_NPL_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_isr_down_" + year + "_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_isr_down_" + year + "_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_isr_down_" + year + "_NPL_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_isr_down_" + year + "_NPL_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_isr_down_" + year + "_Blinded.root";}

        }
	else if(fsr_up == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_fsr_up_" + year + "_NPL_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_fsr_up_" + year + "_ZPlusJetsCR_Blinded.root";
                }       
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_fsr_up_" + year + "_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_fsr_up_" + year + "_NPL_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_fsr_up_" + year + "_NPL_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_fsr_up_" + year + "_Blinded.root";}

        }
        else if(fsr_down == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_fsr_down_" + year + "_NPL_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_fsr_down_" + year + "_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_fsr_down_" + year + "_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        BTagEffOutput = "BTagEffPlots_" + process + "_fsr_down_" + year + "_NPL_ZPlusJetsCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        BTagEffOutput = "BTagEffPlots_" + process + "_fsr_down_" + year + "_NPL_ttbarCR_Blinded.root";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{BTagEffOutput = "BTagEffPlots_" + process + "_fsr_down_" + year + "_Blinded.root";}

        }
	else{

        	if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                	BTagEffOutput = "BTagEffPlots_" + process + "_" + year + "_NPL_Blinded.root";
        	}
        	else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                	BTagEffOutput = "BTagEffPlots_" + process + "_" + year + "_ZPlusJetsCR_Blinded.root";
       	 	}
        	else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                	BTagEffOutput = "BTagEffPlots_" + process + "_" + year + "_ttbarCR_Blinded.root";
        	}
        	else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                	BTagEffOutput = "BTagEffPlots_" + process + "_" + year + "_NPL_ZPlusJetsCR_Blinded.root";
        	}
        	else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                	BTagEffOutput = "BTagEffPlots_" + process + "_" + year + "_NPL_ttbarCR_Blinded.root";
        	}
        	else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
        	else{BTagEffOutput = "BTagEffPlots_" + process + "_" + year + "_Blinded.root";}

	}


}	


	TFile* BTagEffPlots = new TFile(BTagEffOutput.c_str(), "RECREATE");
	int NBins = 40;
	double minpt = 0;
	double maxpt = 500;
	double mineta = -3;
	double maxeta = 3;

        auto h_bjet_ee_num = d_ee_recoZ_jets_bjets_selection.Histo2D({"h_bjet_ee_num", "h_bjet_ee_num", NBins, mineta, maxeta, NBins, minpt, maxpt}, {"BTAGEFF_bjet_eta_num"}, {"BTAGEFF_bjet_pt_num"});
	auto h_nonbjet_ee_num = d_ee_recoZ_jets_bjets_selection.Histo2D({"h_nonbjet_ee_num", "h_nonbjet_ee_num", NBins, mineta, maxeta, NBins, minpt, maxpt}, {"BTAGEFF_nonbjet_eta_num"}, {"BTAGEFF_nonbjet_pt_num"});

        auto h_charm_ee_num = d_ee_recoZ_jets_bjets_selection.Histo2D({"h_charm_ee_num", "h_charm_ee_num", NBins, mineta, maxpt, NBins, minpt, maxpt}, {"BTAGEFF_charm_eta_num"}, {"BTAGEFF_charm_pt_num"});
        auto h_lightjets_ee_num = d_ee_recoZ_jets_bjets_selection.Histo2D({"h_lightjets_ee_num", "h_lightjets_ee_num", NBins, mineta, maxeta, NBins, minpt, maxpt}, {"BTAGEFF_lightjets_eta_num"}, {"BTAGEFF_lightjets_pt_num"});
        auto h_gluon_ee_num = d_ee_recoZ_jets_bjets_selection.Histo2D({"h_gluon_ee_num", "h_gluon_ee_num", NBins, mineta, maxeta, NBins, minpt, maxpt}, {"BTAGEFF_gluon_eta_num"}, {"BTAGEFF_gluon_pt_num"});

	auto h_bjet_mumu_num = d_mumu_recoZ_jets_bjets_selection.Histo2D({"h_bjet_mumu_num", "h_bjet_mumu_num", NBins, mineta, maxeta, NBins, minpt, maxpt}, {"BTAGEFF_bjet_eta_num"}, {"BTAGEFF_bjet_pt_num"});
	auto h_nonbjet_mumu_num = d_mumu_recoZ_jets_bjets_selection.Histo2D({"h_nonbjet_mumu_num", "h_nonbjet_mumu_num", NBins, mineta, maxeta, NBins, minpt, maxpt}, {"BTAGEFF_nonbjet_eta_num"}, {"BTAGEFF_nonbjet_pt_num"});

        auto h_charm_mumu_num = d_mumu_recoZ_jets_bjets_selection.Histo2D({"h_charm_mumu_num", "h_charm_mumu_num", NBins, mineta, maxeta, NBins, minpt, maxpt}, {"BTAGEFF_charm_eta_num"}, {"BTAGEFF_charm_pt_num"});
        auto h_lightjets_mumu_num = d_mumu_recoZ_jets_bjets_selection.Histo2D({"h_lightjets_mumu_num", "h_lightjets_mumu_num", NBins, mineta, maxeta, NBins, minpt, maxpt}, {"BTAGEFF_lightjets_eta_num"}, {"BTAGEFF_lightjets_pt_num"});
        auto h_gluon_mumu_num = d_mumu_recoZ_jets_bjets_selection.Histo2D({"h_gluon_mumu_num", "h_gluon_mumu_num", NBins, mineta, maxeta, NBins, minpt, maxpt}, {"BTAGEFF_gluon_eta_num"}, {"BTAGEFF_gluon_pt_num"});


        auto h_bjet_ee_denom = d_ee_recoZ_jets_bjets_selection.Histo2D({"h_bjet_ee_denom", "h_bjet_ee_denom", NBins, mineta, maxeta, NBins, minpt, maxpt}, {"BTAGEFF_bjet_eta_denom"}, {"BTAGEFF_bjet_pt_denom"});
	auto h_nonbjet_ee_denom = d_ee_recoZ_jets_bjets_selection.Histo2D({"h_nonbjet_ee_denom", "h_nonbjet_ee_denom", NBins, mineta, maxeta, NBins, minpt, maxpt}, {"BTAGEFF_nonbjet_eta_denom"}, {"BTAGEFF_nonbjet_pt_denom"});

        auto h_charm_ee_denom = d_ee_recoZ_jets_bjets_selection.Histo2D({"h_charm_ee_denom", "h_charm_ee_denom", NBins, mineta, maxeta, NBins, minpt, maxpt}, {"BTAGEFF_charm_eta_denom"}, {"BTAGEFF_charm_pt_denom"});
        auto h_lightjets_ee_denom = d_ee_recoZ_jets_bjets_selection.Histo2D({"h_lightjets_ee_denom", "h_lightjets_ee_denom", NBins, mineta, maxeta, NBins, minpt, maxpt}, {"BTAGEFF_lightjets_eta_denom"}, {"BTAGEFF_lightjets_pt_denom"});
        auto h_gluon_ee_denom = d_ee_recoZ_jets_bjets_selection.Histo2D({"h_gluon_ee_denom", "h_gluon_ee_denom", NBins, mineta, maxeta, NBins, minpt, maxpt}, {"BTAGEFF_gluon_eta_denom"}, {"BTAGEFF_gluon_pt_denom"});
 
       auto h_bjet_mumu_denom = d_mumu_recoZ_jets_bjets_selection.Histo2D({"h_bjet_mumu_denom", "h_bjet_mumu_denom", NBins, mineta, maxeta, NBins, minpt, maxpt}, {"BTAGEFF_bjet_eta_denom"}, {"BTAGEFF_bjet_pt_denom"});
	auto h_nonbjet_mumu_denom = d_mumu_recoZ_jets_bjets_selection.Histo2D({"h_bjet_mumu_denom", "h_nonbjet_mumu_denom", NBins, mineta, maxeta, NBins, minpt, maxpt}, {"BTAGEFF_nonbjet_eta_denom"}, {"BTAGEFF_nonbjet_pt_denom"});	

        auto h_charm_mumu_denom = d_mumu_recoZ_jets_bjets_selection.Histo2D({"h_charm_mumu_denom", "h_charm_mumu_denom", NBins, mineta, maxeta, NBins, minpt, maxpt}, {"BTAGEFF_charm_eta_denom"}, {"BTAGEFF_charm_pt_denom"});
        auto h_lightjets_mumu_denom = d_mumu_recoZ_jets_bjets_selection.Histo2D({"h_lightjets_mumu_denom", "h_lightjets_mumu_denom", NBins, mineta, maxeta, NBins, minpt, maxpt}, {"BTAGEFF_lightjets_eta_denom"}, {"BTAGEFF_lightjets_pt_denom"});
        auto h_gluon_mumu_denom = d_mumu_recoZ_jets_bjets_selection.Histo2D({"h_gluon_mumu_denom", "h_gluon_mumu_denom", NBins, mineta, maxeta, NBins, minpt, maxpt}, {"BTAGEFF_gluon_eta_denom"}, {"BTAGEFF_gluon_pt_denom"});


	h_bjet_ee_num->GetXaxis()->SetTitle("#eta");
	h_nonbjet_ee_num->GetXaxis()->SetTitle("#eta");
	h_charm_ee_num->GetXaxis()->SetTitle("#eta");
	h_lightjets_ee_num->GetXaxis()->SetTitle("#eta");
	h_gluon_ee_num->GetXaxis()->SetTitle("#eta");
	h_bjet_mumu_num->GetXaxis()->SetTitle("#eta");
	h_nonbjet_mumu_num->GetXaxis()->SetTitle("#eta");
	h_charm_mumu_num->GetXaxis()->SetTitle("#eta");
	h_lightjets_mumu_num->GetXaxis()->SetTitle("#eta");
	h_gluon_mumu_num->GetXaxis()->SetTitle("#eta");


	h_bjet_ee_denom->GetXaxis()->SetTitle("#eta");
	h_nonbjet_ee_denom->GetXaxis()->SetTitle("#eta");
	h_charm_ee_denom->GetXaxis()->SetTitle("#eta");
	h_lightjets_ee_denom->GetXaxis()->SetTitle("#eta");
	h_gluon_ee_denom->GetXaxis()->SetTitle("#eta");
	h_bjet_mumu_denom->GetXaxis()->SetTitle("#eta");
	h_nonbjet_mumu_denom->GetXaxis()->SetTitle("#eta");
	h_charm_mumu_denom->GetXaxis()->SetTitle("#eta");
	h_lightjets_mumu_denom->GetXaxis()->SetTitle("#eta");
	h_gluon_mumu_denom->GetXaxis()->SetTitle("#eta");


	h_bjet_ee_num->GetYaxis()->SetTitle("p_{T}");
	h_nonbjet_ee_num->GetYaxis()->SetTitle("p_{T}");
	h_charm_ee_num->GetYaxis()->SetTitle("p_{T}");
	h_lightjets_ee_num->GetYaxis()->SetTitle("p_{T}");
	h_gluon_ee_num->GetYaxis()->SetTitle("p_{T}");
	h_bjet_mumu_num->GetYaxis()->SetTitle("p_{T}");
	h_nonbjet_mumu_num->GetYaxis()->SetTitle("p_{T}");
	h_charm_mumu_num->GetYaxis()->SetTitle("p_{T}");
	h_lightjets_mumu_num->GetYaxis()->SetTitle("p_{T}");
	h_gluon_mumu_num->GetYaxis()->SetTitle("p_{T}");

	h_bjet_ee_denom->GetYaxis()->SetTitle("p_{T}");
	h_nonbjet_ee_denom->GetYaxis()->SetTitle("p_{T}");
	h_charm_ee_denom->GetYaxis()->SetTitle("p_{T}");
	h_lightjets_ee_denom->GetYaxis()->SetTitle("p_{T}");
	h_gluon_ee_denom->GetYaxis()->SetTitle("p_{T}");
	h_bjet_mumu_denom->GetYaxis()->SetTitle("p_{T}");
	h_nonbjet_mumu_denom->GetYaxis()->SetTitle("p_{T}");
	h_charm_mumu_denom->GetYaxis()->SetTitle("p_{T}");
	h_lightjets_mumu_denom->GetYaxis()->SetTitle("p_{T}");
	h_gluon_mumu_denom->GetYaxis()->SetTitle("p_{T}");

	h_bjet_ee_num->Write();
	h_nonbjet_ee_num->Write();
	h_charm_ee_num->Write();
	h_lightjets_ee_num->Write();
	h_gluon_ee_num->Write();
	h_bjet_mumu_num->Write();
	h_nonbjet_mumu_num->Write();
	h_charm_mumu_num->Write();
	h_lightjets_mumu_num->Write();
	h_gluon_mumu_num->Write();

	h_bjet_ee_denom->Write();
	h_nonbjet_ee_denom->Write();
	h_charm_ee_denom->Write();
	h_lightjets_ee_denom->Write();
	h_gluon_ee_denom->Write();
	h_bjet_mumu_denom->Write();
	h_nonbjet_mumu_denom->Write();
	h_charm_mumu_denom->Write();
	h_lightjets_mumu_denom->Write();
	h_gluon_mumu_denom->Write();


	TH2F* h_bjet_ee = new TH2F("h_bjet_ee", "h_bjet_ee", NBins, mineta, maxeta, NBins, minpt, maxpt);
	TH2F* h_nonbjet_ee = new TH2F("h_nonbjet_ee", "h_nonbjet_ee", NBins, mineta, maxeta, NBins, minpt, maxpt);
	TH2F* h_charm_ee = new TH2F("h_charm_ee", "h_charm_ee", NBins, mineta, maxeta, NBins, minpt, maxpt);
	TH2F* h_lightjets_ee = new TH2F("h_lightjets_ee", "h_lightjets_ee", NBins, mineta, maxeta, NBins, minpt, maxpt);
	TH2F* h_gluon_ee = new TH2F("h_gluon_ee", "h_gluon_ee", NBins, mineta, maxeta, NBins, minpt, maxpt);

	TH2F* h_bjet_mumu = new TH2F("h_bjet_mumu", "h_bjet_mumu", NBins, mineta, maxeta, NBins, minpt, maxpt);
	TH2F* h_nonbjet_mumu = new TH2F("h_bjet_mumu", "h_nonbjet_mumu", NBins, mineta, maxeta, NBins, minpt, maxpt);
        TH2F* h_charm_mumu = new TH2F("h_charm_mumu", "h_charm_mumu", NBins, mineta, maxeta, NBins, minpt, maxpt);
        TH2F* h_lightjets_mumu = new TH2F("h_lightjets_mumu", "h_lightjets_mumu", NBins, mineta, maxeta, NBins, minpt, maxpt);
        TH2F* h_gluon_mumu = new TH2F("h_gluon_mumu", "h_gluon_mumu", NBins, mineta, maxeta, NBins, minpt, maxpt);	



	h_bjet_ee = (TH2F*)h_bjet_ee_num->Clone();
	h_bjet_ee->Divide(h_bjet_ee_denom.GetPtr());
	h_bjet_mumu = (TH2F*)h_bjet_mumu_num->Clone();
	h_bjet_mumu->Divide(h_bjet_mumu_denom.GetPtr());
	h_charm_ee = (TH2F*)h_charm_ee_num->Clone();
	h_charm_ee->Divide(h_charm_ee_denom.GetPtr());
	h_charm_mumu = (TH2F*)h_charm_mumu_num->Clone();
	h_charm_mumu->Divide(h_charm_mumu_denom.GetPtr());
	h_lightjets_ee = (TH2F*)h_lightjets_ee_num->Clone();
	h_lightjets_ee->Divide(h_lightjets_ee_denom.GetPtr());
	h_lightjets_mumu = (TH2F*)h_lightjets_mumu_num->Clone();
	h_lightjets_mumu->Divide(h_lightjets_mumu_denom.GetPtr());
	h_gluon_ee = (TH2F*)h_gluon_ee_num->Clone();
	h_gluon_ee->Divide(h_gluon_ee_denom.GetPtr());
	h_gluon_mumu = (TH2F*)h_gluon_mumu_num->Clone();
	h_gluon_mumu->Divide(h_gluon_mumu_denom.GetPtr());
	h_nonbjet_ee = (TH2F*)h_nonbjet_ee_num->Clone();
        h_nonbjet_ee->Divide(h_nonbjet_ee_denom.GetPtr());
        h_nonbjet_mumu = (TH2F*)h_nonbjet_mumu_num->Clone();
        h_nonbjet_mumu->Divide(h_nonbjet_mumu_denom.GetPtr());


	h_bjet_ee->SetTitle("h_bjet_ee");
	h_nonbjet_ee->SetTitle("h_nonbjet_ee");
        h_charm_ee->SetTitle("h_charm_ee");
        h_lightjets_ee->SetTitle("h_lightjets_ee");
        h_gluon_ee->SetTitle("h_gluon_ee");
        h_bjet_mumu->SetTitle("h_bjet_mumu");
	h_nonbjet_mumu->SetTitle("h_nonbjet_mumu");
        h_charm_mumu->SetTitle("h_charm_mumu");
        h_lightjets_mumu->SetTitle("h_lightjets_mumu");
        h_gluon_mumu->SetTitle("h_gluon_mumu");
	

	h_bjet_ee->GetXaxis()->SetTitle("#eta");
	h_nonbjet_ee->GetXaxis()->SetTitle("#eta");
	h_charm_ee->GetXaxis()->SetTitle("#eta");
	h_lightjets_ee->GetXaxis()->SetTitle("#eta");
	h_gluon_ee->GetXaxis()->SetTitle("#eta");
	h_bjet_mumu->GetXaxis()->SetTitle("#eta");
	h_nonbjet_mumu->GetXaxis()->SetTitle("#eta");
	h_charm_mumu->GetXaxis()->SetTitle("#eta");
	h_lightjets_mumu->GetXaxis()->SetTitle("#eta");
	h_gluon_mumu->GetXaxis()->SetTitle("#eta");


	h_bjet_ee->GetYaxis()->SetTitle("p_{T}");
	h_nonbjet_ee->GetYaxis()->SetTitle("p_{T}");
	h_charm_ee->GetYaxis()->SetTitle("p_{T}");
	h_lightjets_ee->GetYaxis()->SetTitle("p_{T}");
	h_gluon_ee->GetYaxis()->SetTitle("p_{T}");
	h_bjet_mumu->GetYaxis()->SetTitle("p_{T}");
	h_nonbjet_mumu->GetYaxis()->SetTitle("p_{T}");
	h_charm_mumu->GetYaxis()->SetTitle("p_{T}");
	h_lightjets_mumu->GetYaxis()->SetTitle("p_{T}");
	h_gluon_mumu->GetYaxis()->SetTitle("p_{T}");


	h_bjet_ee->Write();
	h_nonbjet_ee->Write();
	h_charm_ee->Write();
	h_lightjets_ee->Write();
	h_gluon_ee->Write();
	h_bjet_mumu->Write();
	h_nonbjet_mumu->Write();
	h_charm_mumu->Write();
	h_lightjets_mumu->Write();
	h_gluon_mumu->Write();


	BTagEffPlots->Close();




//continue b tag here

auto EffBTaggedFunction_ee{[&h_bjet_ee_num, &h_bjet_ee_denom, &NBins](const floats& DummyColumn, const floats& pts, const floats& etas){

  floats BTaggedEff{};

  for(int i = 0; i < pts.size(); i++){

	int PtNum = h_bjet_ee_num->GetXaxis()->FindBin(pts.at(i));
	int EtaNum = h_bjet_ee_num->GetYaxis()->FindBin(etas.at(i));

	int PtDenom = h_bjet_ee_denom->GetXaxis()->FindBin(pts.at(i));
	int EtaDenom = h_bjet_ee_denom->GetYaxis()->FindBin(etas.at(i));

	float Numerator = h_bjet_ee_num->GetBinContent(PtNum, EtaNum);
	float Denominator = h_bjet_ee_denom->GetBinContent(PtDenom, EtaDenom);

	float eff = Numerator / Denominator;
	
	if(!isnan(eff) && !isinf(eff) && eff > 0){BTaggedEff.push_back(eff);}
        else{BTaggedEff.push_back(1.);}


  }

  return BTaggedEff;


}};




auto EffBTaggedFunction_mumu{[&h_bjet_mumu_num, &h_bjet_mumu_denom, &NBins](const floats& DummyColumn, const floats& pts, const floats& etas){

  floats BTaggedEff{};

  for(int i = 0; i < pts.size(); i++){

        int PtNum = h_bjet_mumu_num->GetXaxis()->FindBin(pts.at(i));
        int EtaNum = h_bjet_mumu_num->GetYaxis()->FindBin(etas.at(i));

        int PtDenom = h_bjet_mumu_denom->GetXaxis()->FindBin(pts.at(i));
        int EtaDenom = h_bjet_mumu_denom->GetYaxis()->FindBin(etas.at(i));

        float Numerator = h_bjet_mumu_num->GetBinContent(PtNum, EtaNum);
        float Denominator = h_bjet_mumu_denom->GetBinContent(PtDenom, EtaDenom);

        float eff = Numerator / Denominator;

        if(!isnan(eff) && !isinf(eff) && eff > 0){BTaggedEff.push_back(eff);}
	else{BTaggedEff.push_back(1.);}


  }

  return BTaggedEff;


}};


auto EffNonBTaggedFunction_ee{[&h_nonbjet_ee_num, &h_nonbjet_ee_denom, &NBins](const floats& DummyColumn, const floats& pts, const floats& etas){

  floats NonBTaggedEff{};

  for(int i = 0; i < pts.size(); i++){

        int PtNum = h_nonbjet_ee_num->GetXaxis()->FindBin(pts.at(i));
        int EtaNum = h_nonbjet_ee_num->GetYaxis()->FindBin(etas.at(i));
    
        int PtDenom = h_nonbjet_ee_denom->GetXaxis()->FindBin(pts.at(i));
        int EtaDenom = h_nonbjet_ee_denom->GetYaxis()->FindBin(etas.at(i));

        float Numerator = h_nonbjet_ee_num->GetBinContent(PtNum, EtaNum);
        float Denominator = h_nonbjet_ee_denom->GetBinContent(PtDenom, EtaDenom);

        float eff = Numerator / Denominator;

        if(!isnan(eff) && !isinf(eff) && eff > 0){NonBTaggedEff.push_back(eff);}
        else{NonBTaggedEff.push_back(1.);}


  }

  return NonBTaggedEff;


}};


auto EffNonBTaggedFunction_mumu{[&h_nonbjet_mumu_num, &h_nonbjet_mumu_denom, &NBins](const floats& DummyColumn, const floats& pts, const floats& etas){

  floats NonBTaggedEff{};

  for(int i = 0; i < pts.size(); i++){

        int PtNum = h_nonbjet_mumu_num->GetXaxis()->FindBin(pts.at(i));
        int EtaNum = h_nonbjet_mumu_num->GetYaxis()->FindBin(etas.at(i));

        int PtDenom = h_nonbjet_mumu_denom->GetXaxis()->FindBin(pts.at(i));
        int EtaDenom = h_nonbjet_mumu_denom->GetYaxis()->FindBin(etas.at(i));

        float Numerator = h_nonbjet_mumu_num->GetBinContent(PtNum, EtaNum);
        float Denominator = h_nonbjet_mumu_denom->GetBinContent(PtDenom, EtaDenom);

        float eff = Numerator / Denominator;

        if(!isnan(eff) && !isinf(eff) && eff > 0){NonBTaggedEff.push_back(eff);}
        else{NonBTaggedEff.push_back(1.);}


  }

  return NonBTaggedEff;


}};



auto EffBTaggedProduct{[](const floats& EffBTagged){

  float initial = 1;

  for(int i = 0; i < EffBTagged.size(); i++ ){

    initial = EffBTagged.at(i) * initial;

  }

  return initial;


}};





auto EffNonBTaggedProduct{[](const floats& EffNonBTagged){

  float initial = 1;

  for(int i = 0; i < EffNonBTagged.size(); i++ ){

  	initial = (1 - EffNonBTagged.at(i)) * initial;

  }

  return initial;


}};



auto ProbBTagMCFunction{[](const float& EffBTaggedProduct, const float& EffNonBTaggedProduct){

  float MCProb = EffBTaggedProduct * EffNonBTaggedProduct; 
  return MCProb;

}};



//reading the csv file to obtain the b tagging scale factor for each event
auto CMSBTagSF_Function{[&BTag_ScaleUp, &BTag_ScaleDown](const floats& pts, const floats etas, const floats CSVv2Discr, bool BTagOrNot, const ints& Jet_partonFlavour){

  floats ResultVector{};

  for(int j = 0; j < Jet_partonFlavour.size(); j++){


	CSVReader reader("./ScaleFactors/BTaggingEfficiency/CSVv2_94XSF_V2_B_F.csv");
	std::vector<std::vector<std::string> > dataList = reader.getData();
        
	std::vector<std::string> OutputVec{}; 
	std::vector<std::string> outputstringvec{};

	std::string number;

        if(BTagOrNot == true){number = "1";}
	else{number = "0";}

	std::vector<std::string> CSVv2OperatingPointTest(pts.size(), number); 

	std::string MeasurementTypeString;

	if(abs(Jet_partonFlavour.at(j)) == 4 || abs(Jet_partonFlavour.at(j)) == 5){MeasurementTypeString = "mujets";}
	else if( (abs(Jet_partonFlavour.at(j)) >= 0 && abs(Jet_partonFlavour.at(j)) < 4) || ( abs(Jet_partonFlavour.at(j)) == 21 ) ){MeasurementTypeString = "incl";}
	else{std::cout << "Not charm, bjet, gluon or light jet. Flavour = " << Jet_partonFlavour.at(j) << std::endl;}

        std::vector<std::string> MeasurementTypeTest(pts.size(), MeasurementTypeString); 

	std::string systematic_type_string;

	if(BTag_ScaleUp == true){systematic_type_string = "up";}
        else if(BTag_ScaleDown == true){systematic_type_string = "down";}
 	else{systematic_type_string = "central";}

        std::vector<std::string> SysTypeTest(pts.size(), "central");
        std::vector<std::string> JetFlavourTest(pts.size(), "0"); 

        std::vector<std::string> EtaTest{};


	for(int i = 0; i < etas.size(); i++){

		std::stringstream ss;
		ss << etas.at(i);
		std::string EtaString(ss.str());

		EtaTest.push_back(EtaString);

	}

	std::vector<std::string> PtTest{};


        for(int i = 0; i < pts.size(); i++){

                std::stringstream ss;
                ss << pts.at(i);
                std::string PtString(ss.str());

                PtTest.push_back(PtString);

        }


	std::vector<std::string> DiscrTest{};

        for(int i = 0; i < pts.size(); i++){


                std::stringstream ss;
                ss << CSVv2Discr.at(i);
                std::string CSVv2DiscrString(ss.str());

                DiscrTest.push_back(CSVv2DiscrString);

        }



	std::vector<std::string> OutVec{};
	std::vector<std::string> FinalOutVec{};



for(int i = 0; i < CSVv2OperatingPointTest.size(); i++){
	
	for(std::vector<std::string> vec : dataList)
	{
		for(std::string data : vec)
		{	
			OutputVec.push_back(data);

		}
	}


		for(std::vector<std::string> vec : dataList)
        	{
                	for(std::string data : vec)
                	{
				
				std::string VecAt1String = vec.at(1);
				std::string VecAt2String = vec.at(2);
				std::string VecAt3String = vec.at(3);
			        std::string VecAt4String = vec.at(4);
				std::string VecAt5String = vec.at(5);
				std::string VecAt6String = vec.at(6);
				std::string VecAt7String = vec.at(7);
				std::string VecAt8String = vec.at(8);
				std::string VecAt9String = vec.at(9);
				

				VecAt1String.erase(0, 0);
				VecAt2String.erase(0, 0);
                                VecAt3String.erase(0, 0);
                                VecAt4String.erase(0, 0);
				VecAt1String.erase(remove(VecAt1String.begin(), VecAt1String.end(), ' '), VecAt1String.end());		
				VecAt2String.erase(remove(VecAt2String.begin(), VecAt2String.end(), ' '), VecAt2String.end());
				VecAt3String.erase(remove(VecAt3String.begin(), VecAt3String.end(), ' '), VecAt3String.end());
                                VecAt4String.erase(remove(VecAt4String.begin(), VecAt4String.end(), ' '), VecAt4String.end());	


				float VecAt4Float = stof(VecAt4String);
                        	float VecAt5Float = stof(VecAt5String);
				float VecAt6Float = stof(VecAt6String);
				float VecAt7Float = stof(VecAt7String);
				float VecAt8Float = stof(VecAt8String);
                        	float VecAt9Float = stof(VecAt9String);


				float PtTestFloat = stof(PtTest.at(i));
				float EtaTestFloat = stof(EtaTest.at(i));
				float DiscrTestFloat = stof(DiscrTest.at(i));



				if( (vec.at(0) == CSVv2OperatingPointTest.at(i)) 
    				&& (VecAt1String == MeasurementTypeTest.at(i))
    				&& (VecAt2String == SysTypeTest.at(i))
    	  			&& (VecAt3String == JetFlavourTest.at(i))
    				&& (VecAt4Float < EtaTestFloat)
    				&& (VecAt5Float > EtaTestFloat)
    	  			&& (VecAt6Float < PtTestFloat)     
    				&& (VecAt7Float > PtTestFloat)
    				&& (VecAt8Float < DiscrTestFloat) 
    				&& (VecAt9Float > DiscrTestFloat)
  				){
					OutVec.push_back(vec.at(10));					
				}
				else if( (vec.at(0) != CSVv2OperatingPointTest.at(i))
                        	|| (vec.at(1) != MeasurementTypeTest.at(i))
                        	|| (vec.at(2) != SysTypeTest.at(i))
                        	|| (vec.at(3) != JetFlavourTest.at(i))
                        	|| (VecAt4Float > EtaTestFloat)
                        	|| (VecAt5Float < EtaTestFloat)
                        	|| (VecAt6Float > PtTestFloat)     
                        	|| (VecAt7Float < PtTestFloat)
                        	|| (VecAt8Float > DiscrTestFloat) 
                        	|| (VecAt9Float < DiscrTestFloat)){OutVec.push_back("0");}
				else{std::cout << "double check criteria" << std::endl;}
			
                	}

        	}



	std::vector<std::string> NewOutVec{};
	std::vector<std::string> Zeroes{}; 
	Zeroes.push_back("0");
	Zeroes.push_back("0");

	bool check = all_of(OutVec.begin(), OutVec.end(), [](std::string s){return s == "0";});


	if(OutVec.size() != 0 && check == false){
		for(int i = 0; i < OutVec.size(); i++){

			if(OutVec.at(i) != "0"){NewOutVec.push_back(OutVec.at(i));}
	
		}
	
		std::string outputString;

	
		if(NewOutVec.size() > 11){

			outputString = NewOutVec.at( ((i+1)*11)-1 );
		}
		else{
			outputString = NewOutVec.at(0); 
		}

		outputString.erase(outputString.begin()+1);
                outputString.erase(outputString.begin());
                outputString.erase(outputString.end()-2);
                outputString.erase(outputString.end()-1);
		
                std::string::size_type pos = 0;
 
                while ((pos = outputString.find('x', pos)) != std::string::npos)
                {
                        outputString.replace(pos, 1, PtTest.at(i));
                        pos += 2;
                }
                
		outputstringvec.push_back(outputString);
		FinalOutVec.push_back(outputstringvec.at(i));

	
	}
	else{FinalOutVec.push_back(Zeroes.at(0));}


}//end of for loop



//Evaluating the mathematical expression in the std::string
std::string ConcatenatedString, ConcatenatedString2, ConcatenatedString3, ConcatenatedString4;
std::string ConcatenatedString5, ConcatenatedString6, ConcatenatedString7, ConcatenatedString8;
std::string ConcatenatedString9, ConcatenatedString10, ConcatenatedString11, ConcatenatedString12;
std::string ConcatenatedString13, ConcatenatedString14, ConcatenatedString15, ConcatenatedString16;
std::vector<char> VecForConcString{};
std::vector<char> VecForConcString2{};
std::vector<char> VecForConcString3{};
std::vector<char> VecForConcString4{};
std::vector<char> VecForConcString5{};
std::vector<char> VecForConcString6{};
std::vector<char> VecForConcString7{};
std::vector<char> VecForConcString8{};
std::vector<char> VecForConcString9{};
std::vector<char> VecForConcString10{};
std::vector<char> VecForConcString11{};
std::vector<char> VecForConcString12{};
std::vector<char> VecForConcString13{};
std::vector<char> VecForConcString14{};
std::vector<char> VecForConcString15{};
std::vector<char> VecForConcString16{};
int index, index2, index3, index4, index5, index6, index7, index8, index9, index10, index11, index12, index13, index14, index15, index16;
float result;

for(int i = 0; i < FinalOutVec.size(); i++){

	std::string FirstElement = FinalOutVec.at(i);

	if(FirstElement.at(0) != '('){

		int LastIndex;

		//first
		for(int i = 0; i < FirstElement.length(); i++){

			if(FirstElement.at(i) != ')'&& 
	   	   	   FirstElement.at(i) != '(' &&
	   	   	   FirstElement.at(i) != '*' &&
	   	   	   FirstElement.at(i) != '/' &&
	   	   	   FirstElement.at(i) != '+' &&
	  	   	   FirstElement.at(i) != '-'){VecForConcString.push_back(FirstElement.at(i)); LastIndex = i;}
			else if(i == 0 && FirstElement.at(i) == '('){continue;}
			else{index = i; break;}
		}


		for(int i = 0; i < VecForConcString.size(); i++){
			if(i == 0){ConcatenatedString = VecForConcString.at(i);}
			else{ConcatenatedString += VecForConcString.at(i);}

		}	


		float ConcatenatedStringToFloat = stof(ConcatenatedString);

		if(LastIndex == FirstElement.length()-1){ResultVector.push_back(ConcatenatedStringToFloat);}
		else{
	
			int Min1;

			if(FirstElement.at(index) == '+' && 
			   FirstElement.at(index+1) == '(' && 
			   FirstElement.at(index+2) == '-' && 
			   FirstElement.at(index+3) == '('){
			
				Min1 = index+4;
			}
			else if(FirstElement.at(index) == '+' &&
                           	FirstElement.at(index+1) == '(' && 
                          	FirstElement.at(index+2) == '(' && 
                           	FirstElement.at(index+3) == '-' &&
				FirstElement.at(index+2) == '('){
                        
                                	Min1 = index+5;
                        }
			else if(FirstElement.at(index) == '*' &&
				FirstElement.at(index+1) == '(' &&
				FirstElement.at(index+2) == '('){
			
					Min1 = index+3;

			}
			else if(FirstElement.at(index) == '+' &&
				FirstElement.at(index+1) == '-'){
			
					Min1 = index+2;

			}
			else if(FirstElement.at(index) == '+' &&
				FirstElement.at(index+1) == '(' &&
				FirstElement.at(index+2) == '(' &&
				FirstElement.at(index+3) == '-' &&
				FirstElement.at(index+4) == '('){

				Min1 = index+5;
			

			}
			else{Min1 = index+1;}

			//second
			for(int i = Min1; i < FirstElement.length(); i++){

        			if(FirstElement.at(i) == 'e' && FirstElement.at(i+1) == '-'){
					
					VecForConcString2.push_back(FirstElement.at(i));
					VecForConcString2.push_back(FirstElement.at(i+1));
					VecForConcString2.push_back(FirstElement.at(i+2));
					VecForConcString2.push_back(FirstElement.at(i+3));
					
					index2 = i+4;
					break;

				}
				else if(FirstElement.at(i) != ')'&&
           	  		   	FirstElement.at(i) != '(' &&
           	  		   	FirstElement.at(i) != '*' &&
           	  		   	FirstElement.at(i) != '/' &&
           	  		   	FirstElement.at(i) != '+' &&
                  		   	FirstElement.at(i) != '-'){VecForConcString2.push_back(FirstElement.at(i));}
               			else{index2 = i; break;}
			
			}



			for(int i = 0; i < VecForConcString2.size(); i++){
        			if(i == 0){ConcatenatedString2 = VecForConcString2.at(i);}
        			else{ConcatenatedString2 += VecForConcString2.at(i);}

			}

			float ConcatenatedStringToFloat2 = stof(ConcatenatedString2);
			int Min2;

			if(FirstElement.at(index2) == '*' && 
			   FirstElement.at(index2+1) == '(' && 
			   FirstElement.at(index2+2) == 'l' && 
			   FirstElement.at(index2+3) == 'o' &&
			   FirstElement.at(index2+4) == 'g' &&
			   FirstElement.at(index2+5) == '('){Min2 = index2+6;}
                        else{Min2 = index2+2;}

			//third
			for(int i = Min2; i < FirstElement.length(); i++){

        			if(FirstElement.at(i) != ')'&&
           	   		   FirstElement.at(i) != '(' &&
           	   		   FirstElement.at(i) != '*' &&
           	   		   FirstElement.at(i) != '/' &&
           	   		   FirstElement.at(i) != '+' &&
           	   		   FirstElement.at(i) != '-'){VecForConcString3.push_back(FirstElement.at(i));}
        			else{index4 = i; break;}
			}


			for(int i = 0; i < VecForConcString3.size(); i++){
        			if(i == 0){ConcatenatedString3 = VecForConcString3.at(i);}
        			else{ConcatenatedString3 += VecForConcString3.at(i);}

			}


			float ConcatenatedStringToFloat3 = stof(ConcatenatedString3);
			int Min3;

			if(FirstElement.at(index4) == '+'){Min3 = index4+1;}
			else{Min3 = index4+1;}


			//fourth
			for(int i = Min3; i < FirstElement.length(); i++){

        			if(FirstElement.at(i) != ')'&&
           	   		   FirstElement.at(i) != '(' &&
           	   		   FirstElement.at(i) != '*' &&
           	   		   FirstElement.at(i) != '/' &&
           	   	           FirstElement.at(i) != '+' &&
           	   		   FirstElement.at(i) != '-'){VecForConcString4.push_back(FirstElement.at(i));}
        			else{index5 = i; break;}
			}


			for(int i = 0; i < VecForConcString4.size(); i++){
        			if(i == 0){ConcatenatedString4 = VecForConcString4.at(i);}
        			else{ConcatenatedString4 += VecForConcString4.at(i);}

			}
			

			float ConcatenatedStringToFloat4 = stof(ConcatenatedString4);


			int Min4;

			if(FirstElement.at(index5) == ')' && 
			   FirstElement.at(index5+1) == '*' && 
			   FirstElement.at(index5+2) == '(' && 
			   FirstElement.at(index5+3) == 'l' && 
			   FirstElement.at(index5+4) == 'o' && 
			   FirstElement.at(index5+5) == 'g' && 
			   FirstElement.at(index5+6) == '('){Min4 = index5+7;}
			else if(FirstElement.at(index5) == ')' &&
				FirstElement.at(index5+1) == ')' &&
				FirstElement.at(index5+2) == '/' &&
				FirstElement.at(index5+3) == '('){Min4 = index5+4;}
                        else{Min4 = index5+2;} 


			//fifth
			for(int i = Min4; i < FirstElement.length(); i++){

        			if(FirstElement.at(i) != ')'&&
           	   		   FirstElement.at(i) != '(' &&
           	   		   FirstElement.at(i) != '*' &&
           	   		   FirstElement.at(i) != '/' &&
           	   		   FirstElement.at(i) != '+' &&
           	   		   FirstElement.at(i) != '-'){VecForConcString5.push_back(FirstElement.at(i));}
        			else{index6 = i; break;}
			}


			for(int i = 0; i < VecForConcString5.size(); i++){
        			if(i == 0){ConcatenatedString5 = VecForConcString5.at(i);}
        			else{ConcatenatedString5 += VecForConcString5.at(i);}

			}

	
			float ConcatenatedStringToFloat5 = stof(ConcatenatedString5);

			int Min5;

			if(FirstElement.at(index6) == '+' && FirstElement.at(index6+1) != '('){Min5 = index6+1;}
			else if(FirstElement.at(index6) == '+' && FirstElement.at(index6+1) == '('){Min5 = index6+2;}
			else{Min5 = index6+1;}


			//sixth
			for(int i = Min5; i < FirstElement.length(); i++){

        			if(FirstElement.at(i) != ')'&&
           	   		   FirstElement.at(i) != '(' &&
           	   		   FirstElement.at(i) != '*' &&
           	   		   FirstElement.at(i) != '/' &&
           	   		   FirstElement.at(i) != '+' &&
           	   		   FirstElement.at(i) != '-'){VecForConcString6.push_back(FirstElement.at(i));}
        			else{index7 = i; break;}
			}


			for(int i = 0; i < VecForConcString6.size(); i++){
        			if(i == 0){ConcatenatedString6 = VecForConcString6.at(i);}
        			else{ConcatenatedString6 += VecForConcString6.at(i);}

			}
	

			float ConcatenatedStringToFloat6 = stof(ConcatenatedString6);

			int Min6;

                        if(FirstElement.at(index7) == ')' &&
                           FirstElement.at(index7+1) == '*' &&
                           FirstElement.at(index7+2) == '('){Min6 = index7+3;}
			else if(FirstElement.at(index7) == '*'){Min6 = index7+1;}
                        else{Min6 = index7+3;}
			

			//seventh
			for(int i = Min6; i < FirstElement.length(); i++){

        			if(FirstElement.at(i) != ')'&&
           	   		   FirstElement.at(i) != '(' &&
           	   		   FirstElement.at(i) != '*' &&
           	   		   FirstElement.at(i) != '/' &&
           	   		   FirstElement.at(i) != '+' &&
           	   		   FirstElement.at(i) != '-'){VecForConcString7.push_back(FirstElement.at(i));}
        			else{index8 = i; break;}
			}


			for(int i = 0; i < VecForConcString7.size(); i++){
        			if(i == 0){ConcatenatedString7 = VecForConcString7.at(i);}
        			else{ConcatenatedString7 += VecForConcString7.at(i);}

			}


			float ConcatenatedStringToFloat7 = stof(ConcatenatedString7);
			int Min7;
                        
			if(FirstElement.at(index8) == '-' &&
                           FirstElement.at(index8+1) == '(' &&
                           FirstElement.at(index8+2) == '-' &&
			   FirstElement.at(index8+3) == '('){Min7 = index8+4;}
                        else{Min7 = index8+1;}
	


			//for an output containing 7 floats
			if(FirstElement.at(index6) == '+' && 
			   FirstElement.at(index6+1) == '(' &&
			   FirstElement.at(index8) == ')' &&
                           FirstElement.at(index8+1) == ')' &&
                           FirstElement.at(index8+2) == ')'){

				
				result = ConcatenatedStringToFloat*((ConcatenatedStringToFloat2+(ConcatenatedStringToFloat3*ConcatenatedStringToFloat4))/(ConcatenatedStringToFloat5+(ConcatenatedStringToFloat6*ConcatenatedStringToFloat7)));

				ResultVector.push_back(result);


			}
			else{


				//eighth
				for(int i = Min7; i < FirstElement.length(); i++){

        				if(FirstElement.at(i) != ')'&&
           	   		   	   FirstElement.at(i) != '(' &&
           	   		   	   FirstElement.at(i) != '*' &&
           	   		   	   FirstElement.at(i) != '/' &&
           	   		   	   FirstElement.at(i) != '+' &&
           	   		   	   FirstElement.at(i) != '-'){VecForConcString8.push_back(FirstElement.at(i));}
        				else{index9 = i; break;}
				}     




				for(int i = 0; i < VecForConcString8.size(); i++){
        				if(i == 0){ConcatenatedString8 = VecForConcString8.at(i);}
        				else{ConcatenatedString8 += VecForConcString8.at(i);}
    
				}   


				float ConcatenatedStringToFloat8 = stof(ConcatenatedString8);
				int Min8;

				if(FirstElement.at(index9) == '*' &&
			   	   FirstElement.at(index9+1) == 'l' &&
                           	   FirstElement.at(index9+2) == 'o' &&
                           	   FirstElement.at(index9+3) == 'g' &&
                           	   FirstElement.at(index9+4) == '('){Min8 = index9+5;}
				else if(FirstElement.at(index9) == ')' &&
					FirstElement.at(index9+1) == ')'){break;}
				else{Min8 = index9+2;}



				//ninth
				for(int i = Min8; i < FirstElement.length(); i++){

        				if(FirstElement.at(i) != ')'&&
           	   		   	   FirstElement.at(i) != '(' &&
           	   		   	   FirstElement.at(i) != '*' &&
           	   		   	   FirstElement.at(i) != '/' &&
           	   		  	   FirstElement.at(i) != '+' &&
           	   		   	   FirstElement.at(i) != '-'){VecForConcString9.push_back(FirstElement.at(i));}
        				else{index10 = i; break;}
				}


				for(int i = 0; i < VecForConcString9.size(); i++){
        				if(i == 0){ConcatenatedString9 = VecForConcString9.at(i);}
        				else{ConcatenatedString9 += VecForConcString9.at(i);}

				}


				float ConcatenatedStringToFloat9 = stof(ConcatenatedString9);
				int Min9;

				if(FirstElement.at(index10) == '+'){Min9 = index10+1;}
				else{Min9 = index10 + 1;}


				int LastIndex2;

				//tenth
				for(int i = Min9; i < FirstElement.length(); i++){

        				if(FirstElement.at(i) != ')'&&
           	   		   	   FirstElement.at(i) != '(' &&
           	   		   	   FirstElement.at(i) != '*' &&
           	   		   	   FirstElement.at(i) != '/' &&
           	   		   	   FirstElement.at(i) != '+' &&
           	   		   	   FirstElement.at(i) != '-'){VecForConcString10.push_back(FirstElement.at(i)); LastIndex2 = i;}
        				else{index11 = i; break;}
				}


				for(int i = 0; i < VecForConcString10.size(); i++){
        				if(i == 0){ConcatenatedString10 = VecForConcString10.at(i);}
        				else{ConcatenatedString10 += VecForConcString10.at(i);}

				}
		

				float ConcatenatedStringToFloat10 = stof(ConcatenatedString10);

				//result for an equation containing 10 floats and ending in "))))))))' 
			 	if(LastIndex2 == FirstElement.length()-9 &&
                            	   FirstElement.at(index) == '+' &&
                            	   FirstElement.at(index+1) == '(' &&
                            	   FirstElement.at(index+2) == '-' &&
                           	   FirstElement.at(index+3) == '('){
				

			    		result = ConcatenatedStringToFloat+(-(ConcatenatedStringToFloat2*(log(ConcatenatedStringToFloat3+ConcatenatedStringToFloat4)*(log(ConcatenatedStringToFloat5+ConcatenatedStringToFloat6)*(ConcatenatedStringToFloat7-(-(ConcatenatedStringToFloat8*log(ConcatenatedStringToFloat9+ConcatenatedStringToFloat10))))))));

                                	ResultVector.push_back(result);

				}
                        	else{


					int Min10;

					if(FirstElement.at(index11+8) == '-'){Min10 = index11+9;}
					else{Min10 = index11+3;}


					int LastIndex3;

					//eleventh
					for(int i = Min10; i < FirstElement.length(); i++){

        					if(FirstElement.at(i) != ')'&&
           	   		   	   	   FirstElement.at(i) != '(' &&
           	   		  	   	   FirstElement.at(i) != '*' &&
           	   		   	   	   FirstElement.at(i) != '/' &&
           	   		   	   	   FirstElement.at(i) != '+' &&
           	   		   	   	   FirstElement.at(i) != '-'){VecForConcString11.push_back(FirstElement.at(i)); LastIndex3 = i;}
        					else{index12 = i; break;}
					}


					for(int i = 0; i < VecForConcString11.size(); i++){
        					if(i == 0){ConcatenatedString11 = VecForConcString11.at(i);}
        					else{ConcatenatedString11 += VecForConcString11.at(i);}

					}

					float ConcatenatedStringToFloat11 = stof(ConcatenatedString11);


					//for an equation with 11 floats
					if(LastIndex == FirstElement.length()-1 && FirstElement.at(index11+8) == '-'){

						result = ConcatenatedStringToFloat+((-(ConcatenatedStringToFloat2*exp(-5)*(log(ConcatenatedStringToFloat3+ConcatenatedStringToFloat4)*(log(ConcatenatedStringToFloat5+ConcatenatedStringToFloat6)*(ConcatenatedStringToFloat7-(-(ConcatenatedStringToFloat8*log(ConcatenatedStringToFloat9+ConcatenatedStringToFloat10))))))))-ConcatenatedStringToFloat11);

						ResultVector.push_back(result);

					}
					else if(FirstElement.at(index) == '+' &&
                                		FirstElement.at(index+1) == '(' &&
                                		FirstElement.at(index+2) == '(' &&
                                		FirstElement.at(index+3) == '-' &&
                                		FirstElement.at(index+4) == '(' &&
						FirstElement.at(index11+8) == '-'){

                                		result = ConcatenatedStringToFloat+((-(ConcatenatedStringToFloat2*(log(ConcatenatedStringToFloat3+ConcatenatedStringToFloat4)*(log(ConcatenatedStringToFloat5+ConcatenatedStringToFloat6)*(ConcatenatedStringToFloat7-(-(ConcatenatedStringToFloat8*log(ConcatenatedStringToFloat9+ConcatenatedStringToFloat10))))))))-ConcatenatedStringToFloat11);

						ResultVector.push_back(result);
                        
                        		}
					else{
				
						//twelfth
						for(int i = index12+1; i < FirstElement.length(); i++){
	
        						if(FirstElement.at(i) != ')'&&
           	   			   	   	   FirstElement.at(i) != '(' &&
           	   			   	   	   FirstElement.at(i) != '*' &&
           	   		   	   		   FirstElement.at(i) != '/' &&
           	   		   		   	   FirstElement.at(i) != '+' &&
           	   		   	   		   FirstElement.at(i) != '-'){VecForConcString12.push_back(FirstElement.at(i));}
        						else{index13 = i; break;}
						}


						for(int i = 0; i < VecForConcString12.size(); i++){
        						if(i == 0){ConcatenatedString12 = VecForConcString12.at(i);}
        						else{ConcatenatedString12 += VecForConcString12.at(i);}

						}

						float ConcatenatedStringToFloat12 = stof(ConcatenatedString12);

						//thirteenth
						for(int i = index13+1; i < FirstElement.length(); i++){

        						if(FirstElement.at(i) != ')'&&
           	   		   	   	   	   FirstElement.at(i) != '(' &&
           	   		   	   	   	   FirstElement.at(i) != '*' &&
           	   		   	   	   	   FirstElement.at(i) != '/' &&
           	   		   	   	  	   FirstElement.at(i) != '+' &&
           	   		   	   	   	   FirstElement.at(i) != '-'){VecForConcString13.push_back(FirstElement.at(i));}
        						else{index14 = i; break;}
						}


						for(int i = 0; i < VecForConcString13.size(); i++){
        						if(i == 0){ConcatenatedString13 = VecForConcString13.at(i);}
        						else{ConcatenatedString13 += VecForConcString13.at(i);}

						}

						float ConcatenatedStringToFloat13 = stof(ConcatenatedString13);

						//Calculating the result for an equation containing 13 floats

						if(FirstElement.at(index) == '+'){result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2;}
						else if(FirstElement.at(index) == '-'){result = ConcatenatedStringToFloat - ConcatenatedStringToFloat2;}
						else if(FirstElement.at(index) == '*'){result = ConcatenatedStringToFloat * ConcatenatedStringToFloat2;}
						else{std::cout << "FirstElement.at(index) is " << FirstElement.at(index) << ". This is not a +, - or *. Output has been set to zero" << std::endl; result = 0;}


						if(FirstElement.at(index2) == '*' && FirstElement.at(index2+1) == '(' && FirstElement.at(index2+2) == '-'){
							result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3);
						}
						else{std::cout << "error" << std::endl;}


						if(FirstElement.at(index4) == '+'){
							result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + ConcatenatedStringToFloat4);
						}
						else{std::cout << "error message 2" << std::endl;}


						if(FirstElement.at(index5) == '*' && FirstElement.at(index5+1) == '('){
							result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + ConcatenatedStringToFloat4*(ConcatenatedStringToFloat5));
						}
						else{std::cout << "error message 3" << std::endl;}

						if(FirstElement.at(index6) == '+'){
        						result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + ConcatenatedStringToFloat4*(ConcatenatedStringToFloat5 + ConcatenatedStringToFloat6));
						}
						else{std::cout << "error message 4" << std::endl;}


						if(FirstElement.at(index7) == '*'){
        						result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + ConcatenatedStringToFloat4*(ConcatenatedStringToFloat5 + ConcatenatedStringToFloat6*(-ConcatenatedStringToFloat7)));
						}	   
						else{std::cout << "error message 5" << std::endl;}


						if(FirstElement.at(index8) == '+'){
 	       						result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + ConcatenatedStringToFloat4*(ConcatenatedStringToFloat5 + ConcatenatedStringToFloat6*(-ConcatenatedStringToFloat7+ConcatenatedStringToFloat8)));
						}
						else{std::cout << "error message 5" << std::endl;}

						if(FirstElement.at(index9) == '*' && FirstElement.at(index9+1) == '('){
							result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + ConcatenatedStringToFloat4*(ConcatenatedStringToFloat5 + ConcatenatedStringToFloat6*(-ConcatenatedStringToFloat7+ConcatenatedStringToFloat8*(ConcatenatedStringToFloat9))));
						}
						else{std::cout << "error message 6" << std::endl;}

						if(FirstElement.at(index10) == '+'){
        						result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + ConcatenatedStringToFloat4*(ConcatenatedStringToFloat5 + ConcatenatedStringToFloat6*(-ConcatenatedStringToFloat7+ConcatenatedStringToFloat8*(ConcatenatedStringToFloat9+ConcatenatedStringToFloat10))));
						}
						else{std::cout << "error message 7" << std::endl;}

						if(FirstElement.at(index11) == '*' && FirstElement.at(index11+1) == '(' && FirstElement.at(index11+2) == '-'){
        						result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + ConcatenatedStringToFloat4*(ConcatenatedStringToFloat5 + ConcatenatedStringToFloat6*(-ConcatenatedStringToFloat7+ConcatenatedStringToFloat8*(ConcatenatedStringToFloat9+ConcatenatedStringToFloat10*(-ConcatenatedStringToFloat11)))));
						}
						else{std::cout << "error message 8" << std::endl;}


						if(FirstElement.at(index12) == '+'){
        						result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + ConcatenatedStringToFloat4*(ConcatenatedStringToFloat5 + ConcatenatedStringToFloat6*(-ConcatenatedStringToFloat7+ConcatenatedStringToFloat8*(ConcatenatedStringToFloat9+ConcatenatedStringToFloat10*(-ConcatenatedStringToFloat11+ConcatenatedStringToFloat12)))));
						}	
						else{std::cout << "error message 9" << std::endl;}

						if(FirstElement.at(index13) == '*'){
        						result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + ConcatenatedStringToFloat4*(ConcatenatedStringToFloat5 + ConcatenatedStringToFloat6*(-ConcatenatedStringToFloat7+ConcatenatedStringToFloat8*(ConcatenatedStringToFloat9+ConcatenatedStringToFloat10*(-ConcatenatedStringToFloat11+ConcatenatedStringToFloat12*ConcatenatedStringToFloat13)))));
						}
						else{std::cout << "error message 10" << std::endl;}

						ResultVector.push_back(result);

						}

				}

			}

	 	}//extra bracket added here

	}
	else{
		//first
        	for(int i = 1; i < FirstElement.length(); i++){

                	if(FirstElement.at(i) != ')'&&
                   	   FirstElement.at(i) != '(' &&
                   	   FirstElement.at(i) != '*' &&
                   	   FirstElement.at(i) != '/' &&
                   	   FirstElement.at(i) != '+' &&
                   	   FirstElement.at(i) != '-'){VecForConcString.push_back(FirstElement.at(i));}
                	else{index = i; break;}
        	}


        	for(int i = 0; i < VecForConcString.size(); i++){
                	if(i == 0){ConcatenatedString = VecForConcString.at(i);}
                	else{ConcatenatedString += VecForConcString.at(i);}

        	}

        	float ConcatenatedStringToFloat = stof(ConcatenatedString);

		int Minimum;

		if(FirstElement.at(index) == '+' &&
		   FirstElement.at(index+1) == '(' &&
		   FirstElement.at(index+2) == '-' &&
		   FirstElement.at(index+3) == '('){

			Minimum = index+4;

		}
		else{Minimum = index+3;}


		//second	
		for(int i = Minimum; i < FirstElement.length(); i++){

                	if(FirstElement.at(i) != ')'&&
                   	   FirstElement.at(i) != '(' &&
                   	   FirstElement.at(i) != '*' &&
                   	   FirstElement.at(i) != '/' &&
                   	   FirstElement.at(i) != '+' &&
                   	   FirstElement.at(i) != '-'){VecForConcString2.push_back(FirstElement.at(i));}
                	else{index2 = i; break;}
        	}


        	for(int i = 0; i < VecForConcString2.size(); i++){
                	if(i == 0){ConcatenatedString2 = VecForConcString2.at(i);}
                	else{ConcatenatedString2 += VecForConcString2.at(i);}

        	}
        	
		float ConcatenatedStringToFloat2 = stof(ConcatenatedString2);

		int Minimum2;

		if(FirstElement.at(index2) == '*' && FirstElement.at(index2+1) == '('){Minimum2 = index2+ 6;}
		else{Minimum2 = index2+2;}


		//third
		for(int i = Minimum2; i < FirstElement.length(); i++){

                	if(FirstElement.at(i) != ')'&&
                   	   FirstElement.at(i) != '(' &&
                   	   FirstElement.at(i) != '*' &&
                   	   FirstElement.at(i) != '/' &&
                   	   FirstElement.at(i) != '+' &&
                   	   FirstElement.at(i) != '-'){VecForConcString3.push_back(FirstElement.at(i));}
                	else{index3 = i; break;}
        	}


        	for(int i = 0; i < VecForConcString3.size(); i++){
                	if(i == 0){ConcatenatedString3 = VecForConcString3.at(i);}
                	else{ConcatenatedString3 += VecForConcString3.at(i);}

        	}

        	float ConcatenatedStringToFloat3 = stof(ConcatenatedString3);

		//fourth
		for(int i = index3+1; i < FirstElement.length(); i++){

                	if(FirstElement.at(i) != ')'&&
                   	   FirstElement.at(i) != '(' &&
                   	   FirstElement.at(i) != '*' &&
                   	   FirstElement.at(i) != '/' &&
                   	   FirstElement.at(i) != '+' &&
                   	   FirstElement.at(i) != '-'){VecForConcString4.push_back(FirstElement.at(i));}
                	else{index4 = i; break;}
        	}


        	for(int i = 0; i < VecForConcString4.size(); i++){
                	if(i == 0){ConcatenatedString4 = VecForConcString4.at(i);}
                	else{ConcatenatedString4 += VecForConcString4.at(i);}

        	}

        	float ConcatenatedStringToFloat4 = stof(ConcatenatedString4);

		int Minimum3;


		if(FirstElement.at(index4) == ')' && FirstElement.at(index4+1) == '*' && FirstElement.at(index4+2) == '('){Minimum3 = index4+7;}
		else{Minimum3 = index4+4;}


		//fifth
		for(int i = Minimum3; i < FirstElement.length(); i++){

                	if(FirstElement.at(i) != ')'&&
                   	   FirstElement.at(i) != '(' &&
                   	   FirstElement.at(i) != '*' &&
                   	   FirstElement.at(i) != '/' &&
                   	   FirstElement.at(i) != '+' &&
                   	   FirstElement.at(i) != '-'){VecForConcString5.push_back(FirstElement.at(i));}
               	 	else{index5 = i; break;}
        	}


        	for(int i = 0; i < VecForConcString5.size(); i++){
                	if(i == 0){ConcatenatedString5 = VecForConcString5.at(i);}
                	else{ConcatenatedString5 += VecForConcString5.at(i);}

        	}

        	float ConcatenatedStringToFloat5 = stof(ConcatenatedString5);

		int Minimum4;
		
		if(FirstElement.at(index4) == ')' && FirstElement.at(index4+1) == '*' && FirstElement.at(index4+2) == '('){Minimum4 = index5+1;}
		else{Minimum4 = index5+2;}	

		
		//sixth
		for(int i = Minimum4; i < FirstElement.length(); i++){

                	if(FirstElement.at(i) != ')'&&
                   	   FirstElement.at(i) != '(' &&
                   	   FirstElement.at(i) != '*' &&
                   	   FirstElement.at(i) != '/' &&
                   	   FirstElement.at(i) != '+' &&
                   	   FirstElement.at(i) != '-'){VecForConcString6.push_back(FirstElement.at(i));}
                	else{index6 = i; break;}
        	}


        	for(int i = 0; i < VecForConcString6.size(); i++){
                	if(i == 0){ConcatenatedString6 = VecForConcString6.at(i);}
                	else{ConcatenatedString6 += VecForConcString6.at(i);}

        	}

        	float ConcatenatedStringToFloat6 = stof(ConcatenatedString6);
	
		//seventh
		for(int i = index6+1; i < FirstElement.length(); i++){

                	if(FirstElement.at(i) != ')'&&
                   	   FirstElement.at(i) != '(' &&
                   	   FirstElement.at(i) != '*' &&
                   	   FirstElement.at(i) != '/' &&
                   	   FirstElement.at(i) != '+' &&
                   	   FirstElement.at(i) != '-'){VecForConcString7.push_back(FirstElement.at(i));}
                	else{index7 = i; break;}
        	}


        	for(int i = 0; i < VecForConcString7.size(); i++){
                	if(i == 0){ConcatenatedString7 = VecForConcString7.at(i);}
                	else{ConcatenatedString7 += VecForConcString7.at(i);}

        	}

        	float ConcatenatedStringToFloat7 = stof(ConcatenatedString7);


		int MinValue;


		if(FirstElement.at(index) == '+' &&
		   FirstElement.at(index+1) == '(' &&
		   FirstElement.at(index+2) == '-' &&
		   FirstElement.at(index+3) == '('){
		
			MinValue = index7+7;

		}
		else{MinValue = index7+5;}

		//seventh
        	for(int i = MinValue; i < FirstElement.length(); i++){

                	if(FirstElement.at(i) != ')'&&
                   	   FirstElement.at(i) != '(' &&
                  	   FirstElement.at(i) != '*' &&
                   	   FirstElement.at(i) != '/' &&
                   	   FirstElement.at(i) != '+' &&
                   	   FirstElement.at(i) != '-'){VecForConcString8.push_back(FirstElement.at(i));}
                	else{index8 = i; break;}
        	}


        	for(int i = 0; i < VecForConcString8.size(); i++){
                	if(i == 0){ConcatenatedString8 = VecForConcString8.at(i);}
                	else{ConcatenatedString8 += VecForConcString8.at(i);}

        	}

        	float ConcatenatedStringToFloat8 = stof(ConcatenatedString8);


		//result

		if(FirstElement.at(index) == '*' && 
		   FirstElement.at(index+1) == '(' && 
		   FirstElement.at(index+2) == '(' &&
		   FirstElement.at(index2) == '+' && 
		   FirstElement.at(index2+1) == '(' &&
		   FirstElement.at(index3) == '*' &&
		   FirstElement.at(index4) == ')' && 
		   FirstElement.at(index4+1) == ')' && 
		   FirstElement.at(index4+2) == '/' && 
		   FirstElement.at(index4+3) == '(' &&
		   FirstElement.at(index5) == '+' && 
		   FirstElement.at(index5+1) == '(' &&
		   FirstElement.at(index6) == '*' &&
		   FirstElement.at(index7) == ')' && 
		   FirstElement.at(index7+1) == ')' && 
		   FirstElement.at(index7+2) == ')' && 
		   FirstElement.at(index7+3) == ')' && 
		   FirstElement.at(index7+4) == '-'){

			result = (ConcatenatedStringToFloat*((ConcatenatedStringToFloat2+(ConcatenatedStringToFloat3*ConcatenatedStringToFloat4))/(ConcatenatedStringToFloat5+(ConcatenatedStringToFloat6*ConcatenatedStringToFloat7))))-ConcatenatedStringToFloat8;

        	}
		else if(FirstElement.at(index7+4) == '+'){

                        result = (ConcatenatedStringToFloat*((ConcatenatedStringToFloat2+(ConcatenatedStringToFloat3*ConcatenatedStringToFloat4))/(ConcatenatedStringToFloat5+(ConcatenatedStringToFloat6*ConcatenatedStringToFloat7))))+ConcatenatedStringToFloat8;

                }
        	else if(FirstElement.at(index) == '+' &&
			FirstElement.at(index+1) == '(' &&
			FirstElement.at(index+2) == '-' &&
			FirstElement.at(index+3) == '('){
	

			//eighth
                	for(int i = index8+5; i < FirstElement.length(); i++){

                        	if(FirstElement.at(i) != ')'&&
                           	   FirstElement.at(i) != '(' &&
                           	   FirstElement.at(i) != '*' &&
                           	   FirstElement.at(i) != '/' &&
                           	   FirstElement.at(i) != '+' &&
                           	   FirstElement.at(i) != '-'){VecForConcString9.push_back(FirstElement.at(i));}
                        	else{index9 = i; break;}
                	}


                	for(int i = 0; i < VecForConcString9.size(); i++){
                        	if(i == 0){ConcatenatedString9 = VecForConcString9.at(i);}
                        	else{ConcatenatedString9 += VecForConcString9.at(i);}

                	}
                	
			float ConcatenatedStringToFloat9 = stof(ConcatenatedString9);


			//ninth
                        for(int i = index9+4; i < FirstElement.length(); i++){

                                if(FirstElement.at(i) != ')'&&
                                   FirstElement.at(i) != '(' &&
                                   FirstElement.at(i) != '*' &&
                                   FirstElement.at(i) != '/' &&
                                   FirstElement.at(i) != '+' &&
                                   FirstElement.at(i) != '-'){VecForConcString10.push_back(FirstElement.at(i));}
                                else{index10 = i; break;}
                        }


                        for(int i = 0; i < VecForConcString10.size(); i++){
                                if(i == 0){ConcatenatedString10 = VecForConcString10.at(i);}
                                else{ConcatenatedString10 += VecForConcString10.at(i);}

                        }

                        float ConcatenatedStringToFloat10 = stof(ConcatenatedString10);

			//tenth
                        for(int i = index10+9; i < FirstElement.length(); i++){
          
                                if(FirstElement.at(i) != ')'&&
                                   FirstElement.at(i) != '(' &&
                                   FirstElement.at(i) != '*' &&
                                   FirstElement.at(i) != '/' &&
                                   FirstElement.at(i) != '+' &&
                                   FirstElement.at(i) != '-'){VecForConcString11.push_back(FirstElement.at(i));}
                                else{index11 = i; break;}
                        } 
    

                        for(int i = 0; i < VecForConcString11.size(); i++){
                                if(i == 0){ConcatenatedString11 = VecForConcString11.at(i);}
                                else{ConcatenatedString11 += VecForConcString11.at(i);}
    
                        }
			
                        float ConcatenatedStringToFloat11 = stof(ConcatenatedString11);


		

			result = (ConcatenatedStringToFloat+(-(ConcatenatedStringToFloat2*(log(ConcatenatedStringToFloat3+ConcatenatedStringToFloat4)*(log(ConcatenatedStringToFloat5+ConcatenatedStringToFloat6)*(ConcatenatedStringToFloat7-(-(ConcatenatedStringToFloat8*log(ConcatenatedStringToFloat9+ConcatenatedStringToFloat10)))))))))+ConcatenatedStringToFloat11;
			

		}
		else{std::cout << "ERROR" << std::endl;}

		ResultVector.push_back(result);

	}	

	VecForConcString.clear();
	VecForConcString2.clear();
	VecForConcString3.clear();
	VecForConcString4.clear();
	VecForConcString5.clear();
	VecForConcString6.clear();
	VecForConcString7.clear();
	VecForConcString8.clear();
	VecForConcString9.clear();
	VecForConcString10.clear();
	VecForConcString11.clear();
	VecForConcString12.clear();
	VecForConcString13.clear();


}//end of for loop


}//end of first for loop

return ResultVector;

}};

std::cout << "before CMSBTagSF" << std::endl;

auto CMSBTagSF{[&CMSBTagSF_Function](const floats& pts, const floats etas, const floats CSVv2Discr, const ints& Jet_partonFlavour){

 return CMSBTagSF_Function(pts, etas, CSVv2Discr, true, Jet_partonFlavour);

}};


auto CMSNonBTagSF{[&CMSBTagSF_Function](const floats& pts, const floats etas, const floats CSVv2Discr, const ints& Jet_partonFlavour){

 return CMSBTagSF_Function(pts, etas, CSVv2Discr, false, Jet_partonFlavour);

}};



auto EffBTaggedProductData{[](const floats& EffBTagged, const floats& CMSBTagSF){

  float initial = 1;
  float output;

  for(int i = 0; i < EffBTagged.size(); i++){

        output = (CMSBTagSF.at(0)*EffBTagged.at(i)) * initial;

  }

  return output;

}};




auto EffNonBTaggedProductData{[](const floats& EffNonBTagged, const floats& CMSNonBTagSF){

  float initial = 1;

  int size = (CMSNonBTagSF.size() < EffNonBTagged.size()) ? CMSNonBTagSF.size() : EffNonBTagged.size();

  for(int i = 0; i < size; i++){

  		initial = (1 - (CMSNonBTagSF.at(i)*EffNonBTagged.at(i)) ) * initial;

  }

  return initial;

}};







auto ProbBTagDataFunction{[](const float& EffBTaggedProductData, const float& EffNonBTaggedProductData){

  float DataProb = EffBTaggedProductData * EffNonBTaggedProductData;
  return DataProb;

}};




auto BTagWeightFunction{[](const float& ProbBTagMC, const float& ProbBTagData){

	float BTagWeight = (ProbBTagData) / (ProbBTagMC);
	
        if( !isnan(BTagWeight) && !isinf(BTagWeight) ){return BTagWeight;}
	else{float One = 1.0; return One;}


}};




//Filtering events with a reconstructed W boson
auto d_ee_recoZ_jets_bjets_recoW_selection_defines = d_ee_recoZ_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "SmearedJetPt"})
                 						.Define("w_reco_jets", find_w_pair, {"SmearedJetPt", "SmearedJetPhi", "SmearedJetEta", "SmearedJetMass", "Jet_jetId", "lead_bjet"})
                 						.Define("w_pair_pt", select<floats>, {"SmearedJetPt", "w_reco_jets"})
                 						.Define("w_pair_eta", select<floats>, {"SmearedJetEta", "w_reco_jets"})
                 						.Define("w_pair_phi", select<floats>, {"SmearedJetPhi", "w_reco_jets"})
                 						.Define("w_pair_mass", select<floats>, {"SmearedJetMass", "w_reco_jets"})
                 						.Define("w_mass", inv_mass, {"w_pair_pt", "w_pair_eta", "w_pair_phi", "w_pair_mass"})
								.Define("WPairJet1", WPairJet1, {"SmearedJetPt", "SmearedJetPhi", "SmearedJetEta", "SmearedJetMass", "Jet_jetId", "lead_bjet"})
                 						.Define("WPairJet2", WPairJet2, {"SmearedJetPt", "SmearedJetPhi", "SmearedJetEta", "SmearedJetMass", "Jet_jetId", "lead_bjet"})
								.Define("WPairJet1Pt", TLorentzVectorPt, {"WPairJet1"})
								.Define("WPairJet1Eta", TLorentzVectorEta, {"WPairJet1"})
								.Define("WPairJet1Phi", TLorentzVectorPhi, {"WPairJet1"})
								.Define("WPairJet1Mass", TLorentzVectorMass, {"WPairJet1"})
								.Define("WPairJet2Pt", TLorentzVectorPt, {"WPairJet2"})
                                                                .Define("WPairJet2Eta", TLorentzVectorEta, {"WPairJet2"})
                                                                .Define("WPairJet2Phi", TLorentzVectorPhi, {"WPairJet2"})
                                                                .Define("WPairJet2Mass", TLorentzVectorMass, {"WPairJet2"})
								.Define("dR_WJet1_WJet2", deltaRcheck_W_function, deltaR_WJet1_WJet2_strings)
								.Define("dWj1j2", DeltaPhi_function2, {"WPairJet1Phi", "WPairJet2Phi"})
								.Define("dR_WJet1_LeadingE", deltaRcheck_W_function2, deltaR_WJet1_LeadingElectron_strings)
								.Define("dR_WJet1_SubleadingE", deltaRcheck_W_function2, deltaR_WJet1_SubleadingElectron_strings)
								.Define("dR_WJet2_LeadingE", deltaRcheck_W_function2, deltaR_WJet2_LeadingElectron_strings)
                                                                .Define("dR_WJet2_SubleadingE", deltaRcheck_W_function2, deltaR_WJet2_SubleadingElectron_strings)
								.Define("dR_WJet1_LeadingJet", deltaRcheck_W_function2, deltaR_WJet1_LeadingJet_strings)
                                                                .Define("dR_WJet1_SubleadingJet", deltaRcheck_W_function2, deltaR_WJet1_SubleadingJet_strings)
                                                                .Define("dR_WJet1_ThirdJet", deltaRcheck_W_function2, deltaR_WJet1_ThirdJet_strings)
                                                                .Define("dR_WJet1_FourthJet", deltaRcheck_W_function2, deltaR_WJet1_FourthJet_strings)
								.Define("dR_WJet2_LeadingJet", deltaRcheck_W_function2, deltaR_WJet2_LeadingJet_strings)
                                                                .Define("dR_WJet2_SubleadingJet", deltaRcheck_W_function2, deltaR_WJet2_SubleadingJet_strings)
                                                                .Define("dR_WJet2_ThirdJet", deltaRcheck_W_function2, deltaR_WJet2_ThirdJet_strings)
                                                                .Define("dR_WJet2_FourthJet", deltaRcheck_W_function2, deltaR_WJet2_FourthJet_strings)
								.Define("dPhi_WJet1_LeadingE", DeltaPhi_doublesandfloat, {"WPairJet1Phi", "LeadingElectronPhi"})
                                                                .Define("dPhi_WJet1_SubleadingE", DeltaPhi_doublesandfloat, {"WPairJet1Phi", "SubleadingElectronPhi"})
                                                                .Define("dPhi_WJet2_LeadingE", DeltaPhi_doublesandfloat, {"WPairJet2Phi", "LeadingElectronPhi"})
                                                                .Define("dPhi_WJet2_SubleadingE", DeltaPhi_doublesandfloat, {"WPairJet2Phi", "SubleadingElectronPhi"})
                                                                .Define("dPhi_WJet1_LeadingJet", DeltaPhi_doublesandfloat, {"WPairJet1Phi", "LeadingJetPhi"})
                                                                .Define("dPhi_WJet1_SubleadingJet", DeltaPhi_doublesandfloat, {"WPairJet1Phi", "SubleadingJetPhi"})
                                                                .Define("dPhi_WJet1_ThirdJet", DeltaPhi_doublesandfloat, {"WPairJet1Phi", "ThirdJetPhi"})
                                                                .Define("dPhi_WJet1_FourthJet", DeltaPhi_doublesandfloat, {"WPairJet1Phi", "FourthJetPhi"})
                                                                .Define("dPhi_WJet2_LeadingJet", DeltaPhi_doublesandfloat, {"WPairJet2Phi", "LeadingJetPhi"})
                                                                .Define("dPhi_WJet2_SubleadingJet", DeltaPhi_doublesandfloat, {"WPairJet2Phi", "SubleadingJetPhi"})
                                                                .Define("dPhi_WJet2_ThirdJet", DeltaPhi_doublesandfloat, {"WPairJet2Phi", "ThirdJetPhi"})
                                                                .Define("dPhi_WJet2_FourthJet", DeltaPhi_doublesandfloat, {"WPairJet2Phi", "FourthJetPhi"})
								.Define("WJet1HT",  HT_double, {"WPairJet1Pt"})
                                                                .Define("WJet2HT",  HT_double, {"WPairJet2Pt"})
								.Define("RecoWHT", RecoWHT, {"w_pair_pt"})
								.Define("mtW", TransverseWMass, {"dPhi_j1j2", "WPairJet1Pt", "WPairJet2Pt"});

auto d_ee_recoZ_jets_bjets_recoW_selection = d_ee_recoZ_jets_bjets_recoW_selection_defines.Filter(w_mass_cut, {"w_mass"}, "W mass cut (ee channel)");

if(ZPlusJetsCR == true){auto d_ee_recoZ_jets_bjets_recoW_selection = d_ee_recoZ_jets_bjets_recoW_selection_defines.Filter(w_mass_cut_ZPlusJetsCR, {"w_mass", "MET_sumEt"}, "W mass cut (ee channel)");}
else{auto d_ee_recoZ_jets_bjets_recoW_selection = d_ee_recoZ_jets_bjets_recoW_selection_defines.Filter(w_mass_cut, {"w_mass"}, "W mass cut (ee channel)");}


auto d_mumu_recoZ_jets_bjets_recoW_selection_defines = d_mumu_recoZ_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", "SmearedJetPt"})
							            .Define("w_reco_jets", find_w_pair, {"SmearedJetPt", "SmearedJetPhi", "SmearedJetEta", "SmearedJetMass", "Jet_jetId", "lead_bjet"})
                                                                    .Define("w_pair_pt", select<floats>, {"SmearedJetPt", "w_reco_jets"})
                                                                    .Define("w_pair_eta", select<floats>, {"SmearedJetEta", "w_reco_jets"})
                                                                    .Define("w_pair_phi", select<floats>, {"SmearedJetPhi", "w_reco_jets"})
                                                                    .Define("w_pair_mass", select<floats>, {"SmearedJetMass", "w_reco_jets"})
                                                                    .Define("w_mass", inv_mass, {"w_pair_pt", "w_pair_eta", "w_pair_phi", "w_pair_mass"})
								    .Define("WPairJet1", WPairJet1, {"SmearedJetPt", "SmearedJetPhi", "SmearedJetEta", "SmearedJetMass", "Jet_jetId", "lead_bjet"})
                                                                    .Define("WPairJet2", WPairJet2, {"SmearedJetPt", "SmearedJetPhi", "SmearedJetEta", "SmearedJetMass", "Jet_jetId", "lead_bjet"})
								    .Define("WPairJet1Pt", TLorentzVectorPt, {"WPairJet1"})
                                                                    .Define("WPairJet1Eta", TLorentzVectorEta, {"WPairJet1"})
                                                                    .Define("WPairJet1Phi", TLorentzVectorPhi, {"WPairJet1"})
                                                                    .Define("WPairJet1Mass", TLorentzVectorMass, {"WPairJet1"})
                                                                    .Define("WPairJet2Pt", TLorentzVectorPt, {"WPairJet2"})
                                                                    .Define("WPairJet2Eta", TLorentzVectorEta, {"WPairJet2"})
                                                                    .Define("WPairJet2Phi", TLorentzVectorPhi, {"WPairJet2"})
                                                                    .Define("WPairJet2Mass", TLorentzVectorMass, {"WPairJet2"})
								    .Define("dR_WJet1_WJet2", deltaRcheck_W_function, deltaR_WJet1_WJet2_strings)
                                                                    .Define("dWj1j2", DeltaPhi_function2, {"WPairJet1Phi", "WPairJet2Phi"})
                                                                    .Define("dR_WJet1_LeadingMu", deltaRcheck_W_function2, deltaR_WJet1_LeadingMuon_strings)
                                                                    .Define("dR_WJet1_SubleadingMu", deltaRcheck_W_function2, deltaR_WJet1_SubleadingMuon_strings)
                                                                    .Define("dR_WJet2_LeadingMu", deltaRcheck_W_function2, deltaR_WJet2_LeadingMuon_strings)
                                                                    .Define("dR_WJet2_SubleadingMu", deltaRcheck_W_function2, deltaR_WJet2_SubleadingMuon_strings)
                                                                    .Define("dR_WJet1_LeadingJet", deltaRcheck_W_function2, deltaR_WJet1_LeadingJet_strings)
                                                                    .Define("dR_WJet1_SubleadingJet", deltaRcheck_W_function2, deltaR_WJet1_SubleadingJet_strings)
                                                                    .Define("dR_WJet1_ThirdJet", deltaRcheck_W_function2, deltaR_WJet1_ThirdJet_strings)
                                                                    .Define("dR_WJet1_FourthJet", deltaRcheck_W_function2, deltaR_WJet1_FourthJet_strings)
                                                                    .Define("dR_WJet2_LeadingJet", deltaRcheck_W_function2, deltaR_WJet2_LeadingJet_strings)
                                                                    .Define("dR_WJet2_SubleadingJet", deltaRcheck_W_function2, deltaR_WJet2_SubleadingJet_strings)
                                                                    .Define("dR_WJet2_ThirdJet", deltaRcheck_W_function2, deltaR_WJet2_ThirdJet_strings)
                                                                    .Define("dR_WJet2_FourthJet", deltaRcheck_W_function2, deltaR_WJet2_FourthJet_strings)
                                                                    .Define("dPhi_WJet1_LeadingMu", DeltaPhi_doublesandfloat, {"WPairJet1Phi", "LeadingMuonPhi"})
                                                                    .Define("dPhi_WJet1_SubleadingMu", DeltaPhi_doublesandfloat, {"WPairJet1Phi", "SubleadingMuonPhi"})
                                                                    .Define("dPhi_WJet2_LeadingMu", DeltaPhi_doublesandfloat, {"WPairJet2Phi", "LeadingMuonPhi"})
                                                                    .Define("dPhi_WJet2_SubleadingMu", DeltaPhi_doublesandfloat, {"WPairJet2Phi", "SubleadingMuonPhi"})
                                                                    .Define("dPhi_WJet1_LeadingJet", DeltaPhi_doublesandfloat, {"WPairJet1Phi", "LeadingJetPhi"})
                                                                    .Define("dPhi_WJet1_SubleadingJet", DeltaPhi_doublesandfloat, {"WPairJet1Phi", "SubleadingJetPhi"})
                                                                    .Define("dPhi_WJet1_ThirdJet", DeltaPhi_doublesandfloat, {"WPairJet1Phi", "ThirdJetPhi"})
                                                                    .Define("dPhi_WJet1_FourthJet", DeltaPhi_doublesandfloat, {"WPairJet1Phi", "FourthJetPhi"})
                                                                    .Define("dPhi_WJet2_LeadingJet", DeltaPhi_doublesandfloat, {"WPairJet2Phi", "LeadingJetPhi"})
                                                                    .Define("dPhi_WJet2_SubleadingJet", DeltaPhi_doublesandfloat, {"WPairJet2Phi", "SubleadingJetPhi"})
                                                                    .Define("dPhi_WJet2_ThirdJet", DeltaPhi_doublesandfloat, {"WPairJet2Phi", "ThirdJetPhi"})
                                                                    .Define("dPhi_WJet2_FourthJet", DeltaPhi_doublesandfloat, {"WPairJet2Phi", "FourthJetPhi"})
								    .Define("WJet1HT",  HT_double, {"WPairJet1Pt"})
								    .Define("WJet2HT",  HT_double, {"WPairJet2Pt"})
								    .Define("RecoWHT", RecoWHT, {"w_pair_pt"})
								    .Define("mtW", TransverseWMass, {"dPhi_j1j2", "WPairJet1Pt", "WPairJet2Pt"});



auto d_mumu_recoZ_jets_bjets_recoW_selection = d_mumu_recoZ_jets_bjets_recoW_selection_defines.Filter(w_mass_cut, {"w_mass"}, "W mass cut (mumu channel)");

if(ZPlusJetsCR == true){auto d_mumu_recoZ_jets_bjets_recoW_selection = d_mumu_recoZ_jets_bjets_recoW_selection_defines.Filter(w_mass_cut_ZPlusJetsCR, {"w_mass", "MET_sumEt"}, "W mass cut (mumu channel)");}
else{auto d_mumu_recoZ_jets_bjets_recoW_selection = d_mumu_recoZ_jets_bjets_recoW_selection_defines.Filter(w_mass_cut, {"w_mass"}, "W mass cut (mumu channel)");}




//Print cut report
std::cout << "before print cut flow report" << std::endl;

auto allCutsReport = d.Report();

std::cout << "after allCutsReport. Need to change dataframe input when not running on a range." << std::endl;


for (auto&& cutInfo: allCutsReport){

        CutFlowReport << cutInfo.GetName() << '\t' << cutInfo.GetAll() << '\t' << cutInfo.GetPass() << '\t' << cutInfo.GetEff() << " %" << std::endl;

}


std::cout << "after the for loop for cut flow report" << std::endl;




//Filtering events with a reconstructed top quark
auto d_ee_recoZ_jets_bjets_recoW_recoT_selection = d_ee_recoZ_jets_bjets_recoW_selection.Define("RecoW", WLorentzVector, {"w_pair_pt", "w_pair_eta", "w_pair_phi", "w_mass", "w_reco_jets"})
											.Define("bjetmass", bjet_variable, bjet_mass_strings)
											.Define("bjetpt", bjet_variable, bjet_pt_strings)
											.Define("bjeteta", bjet_variable, bjet_eta_strings)
											.Define("bjetphi", bjet_variable, bjet_phi_strings)
											.Define("BJets", BLorentzVector, {"bjetpt", "bjeteta", "bjetphi", "bjetmass"})
											.Define("RecoTop", top_reconstruction_function, top_strings)
											.Define("Top_Pt", TLorentzVectorPt, {"RecoTop"})
											.Define("Top_Eta", TLorentzVectorEta, {"RecoTop"})
											.Define("Top_Phi", TLorentzVectorPhi, {"RecoTop"})
											.Define("Top_Mass", TLorentzVectorMass, {"RecoTop"})
											.Define("Top_HT", HT_double, {"Top_Pt"})
											.Define("dR_Top_LeadingElectron", deltaRcheck_Top_function, deltaR_Top_LeadingElectron_strings)
											.Define("dR_Top_SubleadingElectron", deltaRcheck_Top_function, deltaR_Top_SubleadingElectron_strings)
											.Define("dR_Top_LeadingJet", deltaRcheck_Top_function, deltaR_Top_LeadingJet_strings)
											.Define("dR_Top_SubleadingJet", deltaRcheck_Top_function, deltaR_Top_SubleadingJet_strings)
											.Define("dR_Top_ThirdJet", deltaRcheck_Top_function, deltaR_Top_ThirdJet_strings)
											.Define("dR_Top_FourthJet", deltaRcheck_Top_function, deltaR_Top_FourthJet_strings)
										        .Define("dR_Top_W", deltaRcheck_WTop_function, deltaR_Top_W_strings)
											.Define("dPhi_Wj1_Top", DeltaPhi_function2, {"WPairJet1Phi", "Top_Phi"})
											.Define("dPhi_Wj2_Top", DeltaPhi_function2, {"WPairJet2Phi", "Top_Phi"})
										        .Define("dR_Z_Top", deltaRcheck_W_function, deltaR_Z_Top_strings)
											.Define("dPhi_Z_Top", DeltaPhi_function2, {"Top_Phi", "RecoZPhi"})
											.Define("dR_Z_WPairJet1", deltaRcheck_W_function, dR_Z_WPairJet1_strings)
											.Define("dR_Z_WPairJet2", deltaRcheck_W_function, dR_Z_WPairJet2_strings)
											.Define("dPhi_Z_WPairJet1", DeltaPhi_function2, {"RecoZPhi", "WPairJet1Phi"})
                                                                                        .Define("dPhi_Z_WPairJet2", DeltaPhi_function2, {"RecoZPhi", "WPairJet2Phi"})
											.Define("MinDeltaR", MinDeltaR, MinDeltaR_strings)
											.Define("MinDeltaPhi", MinDeltaPhi, MinDeltaPhi_strings)
											.Define("dR_LeadingLepton_LeadingBJet", dR_LeadingLepton_LeadingBJet, dR_LeadingLepton_LeadingBJetstrings_ee)
											.Define("dR_SubleadingLepton_LeadingBJet", dR_SubleadingLepton_LeadingBJet, dR_SubleadingLepton_LeadingBJetstrings_ee)
											.Define("DeltaPhi_Leadinglepton_BJet", DeltaPhi_Lepton_BJet, DeltaPhi_Leadinglepton_BJet_string_ee)
                                                                                        .Define("DeltaPhi_Subleadinglepton_BJet", DeltaPhi_Lepton_BJet, DeltaPhi_Subleadinglepton_BJet_string_ee)		
										    	.Define("MET", MET_function, {"MET_sumEt"})
										        .Define("LeadingBJetOutputDiscriminant", LeadingBJetOutputDiscriminant, BJetOutputDiscriminantStrings_Leading)
                                                        			        .Define("SubleadingBJetOutputDiscriminant", SubleadingBJetOutputDiscriminant, BJetOutputDiscriminantStrings_Subleading)
                                                        				.Define("ThirdBJetOutputDiscriminant", ThirdBJetOutputDiscriminant, BJetOutputDiscriminantStrings_Third)
                                                        				.Define("FourthBJetOutputDiscriminant", FourthBJetOutputDiscriminant, BJetOutputDiscriminantStrings_Fourth)
                                                        				.Define("LeadingBJetOutputDiscriminant2", BJetOutputDiscriminant, {"LeadingBJetOutputDiscriminant", "Jet_btagCSVV2"})
                                                        				.Define("SubleadingBJetOutputDiscriminant2", BJetOutputDiscriminant, {"SubleadingBJetOutputDiscriminant", "Jet_btagCSVV2"})
                                                        				.Define("ThirdBJetOutputDiscriminant2", BJetOutputDiscriminant, {"ThirdBJetOutputDiscriminant", "Jet_btagCSVV2"})
                                                        				.Define("FourthBJetOutputDiscriminant2", BJetOutputDiscriminant, {"FourthBJetOutputDiscriminant", "Jet_btagCSVV2"})
                                                                                        .Define("dPhi_W_Top", DeltaPhi_function4, {"w_pair_phi", "Top_Phi"})
										        .Define("dR_Z_LeadingJet", deltaRcheck_W_function2, deltaR_Z_LeadingJet_strings)
                                                                            	        .Define("dR_Z_SubleadingJet", deltaRcheck_W_function2, deltaR_Z_SubleadingJet_strings)
                                                                            	        .Define("dR_Z_ThirdJet", deltaRcheck_W_function2, deltaR_Z_ThirdJet_strings)
                                                                            	        .Define("dR_Z_FourthJet", deltaRcheck_W_function2, deltaR_Z_FourthJet_strings)
                                                                            	        .Define("dPhi_LeadingJet_Z", DeltaPhi_doublesandfloat, {"RecoZPhi", "LeadingJetPhi"})
                                                                            		.Define("dPhi_SubleadingJet_Z", DeltaPhi_doublesandfloat, {"RecoZPhi", "SubleadingJetPhi"})
                                                                            	        .Define("dPhi_ThirdJet_Z", DeltaPhi_doublesandfloat, {"RecoZPhi", "ThirdJetPhi"})
                                                                            	        .Define("dPhi_FourthJet_Z", DeltaPhi_doublesandfloat, {"RecoZPhi", "FourthJetPhi"})
                                                                            	        .Define("dR_W_Z", deltaRcheck_WTop_function, deltaR_W_Z_strings)
                                                                            	        .Define("RecoZHT", HT_double, {"RecoZPt"})
                                                                            		.Define("dPhi_W_Z", DeltaPhi_function4, {"w_pair_phi", "RecoZPhi"})
                                                                                        .Define("TotalEta_System", TotalEta_System, TotalEta_System_strings)
                                                                                        .Define("TotalPhi_System", TotalPhi_System, TotalPhi_System_strings)
											.Define("InvTopMass", inv_mass_doubles, {"Top_Pt", "Top_Eta", "Top_Phi", "Top_Mass"});


std::cout << "before d_mumu_recoZ_jets_bjets_recoW_recoT_selection" << std::endl;


auto d_mumu_recoZ_jets_bjets_recoW_recoT_selection = d_mumu_recoZ_jets_bjets_recoW_selection.Define("RecoW", WLorentzVector, {"w_pair_pt", "w_pair_eta", "w_pair_phi", "w_mass", "w_reco_jets"})
                                                                                            .Define("bjetmass", bjet_variable, bjet_mass_strings)
                                                                                            .Define("bjetpt", bjet_variable, bjet_pt_strings)
                                                                                            .Define("bjeteta", bjet_variable, bjet_eta_strings)
                                                                                            .Define("bjetphi", bjet_variable, bjet_phi_strings)
                                                                                            .Define("BJets", BLorentzVector, {"bjetpt", "bjeteta", "bjetphi", "bjetmass"})
											    .Define("RecoTop", top_reconstruction_function, top_strings)
											    .Define("Top_Pt", TLorentzVectorPt, {"RecoTop"})
											    .Define("Top_Eta", TLorentzVectorEta, {"RecoTop"})
                                                                                            .Define("Top_Phi", TLorentzVectorPhi, {"RecoTop"})
											    .Define("Top_Mass", TLorentzVectorMass, {"RecoTop"})
											    .Define("Top_HT", HT_double, {"Top_Pt"})
											    .Define("dR_Top_LeadingMuon", deltaRcheck_Top_function, deltaR_Top_LeadingMuon_strings)
                                                                                            .Define("dR_Top_SubleadingMuon", deltaRcheck_Top_function, deltaR_Top_SubleadingMuon_strings)
                                                                                            .Define("dR_Top_LeadingJet", deltaRcheck_Top_function, deltaR_Top_LeadingJet_strings)
                                                                                            .Define("dR_Top_SubleadingJet", deltaRcheck_Top_function, deltaR_Top_SubleadingJet_strings)
                                                                                            .Define("dR_Top_ThirdJet", deltaRcheck_Top_function, deltaR_Top_ThirdJet_strings)
                                                                                            .Define("dR_Top_FourthJet", deltaRcheck_Top_function, deltaR_Top_FourthJet_strings)
											    .Define("dR_Top_W", deltaRcheck_WTop_function, deltaR_Top_W_strings)
											    .Define("dPhi_Wj1_Top", DeltaPhi_function2, {"WPairJet1Phi", "Top_Phi"})
											    .Define("dPhi_Wj2_Top", DeltaPhi_function2, {"WPairJet2Phi", "Top_Phi"})
											    .Define("dR_Z_Top", deltaRcheck_W_function, deltaR_Z_Top_strings)
											    .Define("dPhi_Z_Top", DeltaPhi_function2, {"Top_Phi", "RecoZPhi"})
											    .Define("MinDeltaR", MinDeltaR, MinDeltaR_strings)
										            .Define("MinDeltaPhi", MinDeltaPhi, MinDeltaPhi_strings)
											    .Define("dR_LeadingLepton_LeadingBJet", dR_LeadingLepton_LeadingBJet, dR_LeadingLepton_LeadingBJetstrings_mumu)
                                                                                            .Define("dR_SubleadingLepton_LeadingBJet", dR_SubleadingLepton_LeadingBJet, dR_SubleadingLepton_LeadingBJetstrings_mumu)
											    .Define("DeltaPhi_Leadinglepton_BJet", DeltaPhi_Lepton_BJet, DeltaPhi_Leadinglepton_BJet_string_mumu)
											    .Define("DeltaPhi_Subleadinglepton_BJet", DeltaPhi_Lepton_BJet, DeltaPhi_Subleadinglepton_BJet_string_mumu)
											    .Define("MET", MET_function, {"MET_sumEt"})
											    .Define("LeadingBJetOutputDiscriminant", LeadingBJetOutputDiscriminant, BJetOutputDiscriminantStrings_Leading)
                                                        				    .Define("SubleadingBJetOutputDiscriminant", SubleadingBJetOutputDiscriminant, BJetOutputDiscriminantStrings_Subleading)
                                                        				    .Define("ThirdBJetOutputDiscriminant", ThirdBJetOutputDiscriminant, BJetOutputDiscriminantStrings_Third)
                                                        				    .Define("FourthBJetOutputDiscriminant", FourthBJetOutputDiscriminant, BJetOutputDiscriminantStrings_Fourth)
                                                        				    .Define("LeadingBJetOutputDiscriminant2", BJetOutputDiscriminant, {"LeadingBJetOutputDiscriminant", "Jet_btagCSVV2"})
                                                        				    .Define("SubleadingBJetOutputDiscriminant2", BJetOutputDiscriminant, {"SubleadingBJetOutputDiscriminant", "Jet_btagCSVV2"})
                                                        				    .Define("ThirdBJetOutputDiscriminant2", BJetOutputDiscriminant, {"ThirdBJetOutputDiscriminant", "Jet_btagCSVV2"})
                                                        				    .Define("FourthBJetOutputDiscriminant2", BJetOutputDiscriminant, {"FourthBJetOutputDiscriminant", "Jet_btagCSVV2"})
											    .Define("dPhi_W_Top", DeltaPhi_function4, {"w_pair_phi", "Top_Phi"})
											    .Define("dR_Z_LeadingJet", deltaRcheck_W_function2, deltaR_Z_LeadingJet_strings)
                                                                                	    .Define("dR_Z_SubleadingJet", deltaRcheck_W_function2, deltaR_Z_SubleadingJet_strings)
                                                                                	    .Define("dR_Z_ThirdJet", deltaRcheck_W_function2, deltaR_Z_ThirdJet_strings)
                                                                                	    .Define("dR_Z_FourthJet", deltaRcheck_W_function2, deltaR_Z_FourthJet_strings)
                                                                                	    .Define("dPhi_LeadingJet_Z", DeltaPhi_doublesandfloat, {"RecoZPhi", "LeadingJetPhi"})
                                                                                	    .Define("dPhi_SubleadingJet_Z", DeltaPhi_doublesandfloat, {"RecoZPhi", "SubleadingJetPhi"})
                                                                                	    .Define("dPhi_ThirdJet_Z", DeltaPhi_doublesandfloat, {"RecoZPhi", "ThirdJetPhi"})
                                                                                	    .Define("dPhi_FourthJet_Z", DeltaPhi_doublesandfloat, {"RecoZPhi", "FourthJetPhi"})
                                                                                	    .Define("dR_W_Z", deltaRcheck_WTop_function, deltaR_W_Z_strings)
                                                                                	    .Define("RecoZHT", HT_double, {"RecoZPt"})
                                                                                	    .Define("dR_Z_WPairJet1", deltaRcheck_W_function, dR_Z_WPairJet1_strings)
                                                                                	    .Define("dR_Z_WPairJet2", deltaRcheck_W_function, dR_Z_WPairJet2_strings)
                                                                                	    .Define("dPhi_Z_WPairJet1", DeltaPhi_function2, {"RecoZPhi", "WPairJet1Phi"})
                                                                                	    .Define("dPhi_Z_WPairJet2", DeltaPhi_function2, {"RecoZPhi", "WPairJet2Phi"})
                                                                                	    .Define("dPhi_W_Z", DeltaPhi_function4, {"w_pair_phi", "RecoZPhi"})
                                                                                            .Define("TotalEta_System", TotalEta_System, TotalEta_System_strings)
                                                                                            .Define("TotalPhi_System", TotalPhi_System, TotalPhi_System_strings)
											    .Define("InvTopMass", inv_mass_doubles, {"Top_Pt", "Top_Eta", "Top_Phi", "Top_Mass"});



//lambda functions for top quark pT reweighting
auto UnweightedTopPt{[](const doubles& pts){

        return pts;

}};

std::cout << "before d_TopReweighted_ee" << std::endl;

auto d_TopReweighted_ee = d_ee_recoZ_jets_bjets_recoW_recoT_selection.Define("UnweightedTopPt", UnweightedTopPt, {"Top_Pt"});
auto d_TopReweighted_mumu = d_mumu_recoZ_jets_bjets_recoW_recoT_selection.Define("UnweightedTopPt", UnweightedTopPt, {"Top_Pt"});

std::cout << "after d_TopReweighted_mumu" << std::endl;

auto h_WeightedTop_ee = d_TopReweighted_ee.Histo1D({"h_WeightedTop_ee", "h_WeightedTop_ee", NBins, 0, 200}, "Top_Pt");
auto h_WeightedTop_mumu = d_TopReweighted_mumu.Histo1D({"h_WeightedTop_mumu", "h_WeightedTop_mumu", NBins, 0, 200}, "Top_Pt");

std::cout << "after h_WeightedTop_mumu" << std::endl;


if(process == "ttbar_2l2nu" ||
    process == "ttbar_madgraph_NanoAODv5" ||
    process == "ttbar_TTToHadronic" ||
    process == "ttbar_TTToSemileptonic"){

	auto TopReweighting_topquark{[](

		const int& GenPart_pdgId,
		const int& GenPart_statusFlags,
		const doubles& Top_pt

	){

		return GenPart_pdgId == 6 && GenPart_statusFlags == 13 && Top_pt > 0;

	}};

	auto TopReweighting_antitopquark{[](

		const int& GenPart_pdgId,
		const int& GenPart_statusFlags,
		const doubles& Top_pt

	){

		return GenPart_pdgId == -6 && GenPart_statusFlags == 13 && Top_pt > 0;

	}};



	auto TopReweighting_weight{[](

		const ints& TopReweighting_topquark,
		const ints& TopReweighting_antitopquark

	){

		doubles SF_top = exp(-0.0615-(0.00005* TopReweighting_topquark) );
		doubles SF_antitop = exp(-0.0615-(0.00005* TopReweighting_antitopquark) );
		doubles weight = sqrt( SF_top * SF_antitop);
		
		return weight;

	}};



	std::cout << "before d_TopReweighted_ee" << std::endl;

	d_TopReweighted_ee = d_ee_recoZ_jets_bjets_recoW_recoT_selection.Define("TopReweighting_topquark", TopReweighting_topquark, {"GenPart_pdgId", "GenPart_statusFlags", "Top_Pt"})
									.Define("TopReweighting_antitopquark", TopReweighting_antitopquark, {"GenPart_pdgId", "GenPart_statusFlags", "Top_Pt"})
									.Define("TopWeight", TopReweighting_weight, {"TopReweighting_topquark", "TopReweighting_antitopquark"});



	d_TopReweighted_mumu = d_mumu_recoZ_jets_bjets_recoW_recoT_selection.Define("TopReweighting_topquark", TopReweighting_topquark, {"GenPart_pdgId", "GenPart_statusFlags", "Top_Pt"})
                                                                            .Define("TopReweighting_antitopquark", TopReweighting_antitopquark, {"GenPart_pdgId", "GenPart_statusFlags", "Top_Pt"})
                                                                            .Define("TopWeight", TopReweighting_weight, {"TopReweighting_topquark", "TopReweighting_antitopquark"});

	
	std::cout << "after d_TopReweighted_mumu" << std::endl;

	h_WeightedTop_ee = d_TopReweighted_ee.Histo1D({"h_WeightedTop_ee", "h_WeightedTop_ee", NBins, 0, 200}, "Top_Pt", "TopWeight");

	std::cout << "h_WeightedTop_ee" << std::endl;	

	h_WeightedTop_mumu = d_TopReweighted_mumu.Histo1D({"h_WeightedTop_mumu", "h_WeightedTop_mumu", NBins, 0, 200}, "Top_Pt", "TopWeight");

	std::cout << "h_WeightedTop_mumu" << std::endl;	

	auto WeightedTop_ee_Function{[&h_WeightedTop_ee](){
	
		const auto NumOfBins = h_WeightedTop_ee->GetNbinsX(); 
		doubles Vec_ee(NumOfBins);

		for (int i = 0; i < NumOfBins; i++){
			auto content = h_WeightedTop_ee->GetBinContent(i);
			Vec_ee.push_back(content);
		}

		return Vec_ee;

	}};

	auto WeightedTop_mumu_Function{[&h_WeightedTop_mumu](){

                const auto NumOfBins = h_WeightedTop_mumu->GetNbinsX();
                doubles Vec_mumu(NumOfBins);

                for (int i = 0; i < NumOfBins; i++){
                        auto content = h_WeightedTop_mumu->GetBinContent(i);
                        Vec_mumu.push_back(content);
                }
                
                return Vec_mumu;

        }};

	d_TopReweighted_ee = d_ee_recoZ_jets_bjets_recoW_recoT_selection.Define("ReweightedTopPt", WeightedTop_ee_Function, {});
	d_TopReweighted_mumu = d_mumu_recoZ_jets_bjets_recoW_recoT_selection.Define("ReweightedTopPt", WeightedTop_mumu_Function, {});

}
else{

	std::cout << "in else, before d_TopReweighted_ee" << std::endl;

	d_TopReweighted_ee = d_ee_recoZ_jets_bjets_recoW_recoT_selection.Define("ReweightedTopPt", UnweightedTopPt, {"Top_Pt"});
	d_TopReweighted_mumu = d_mumu_recoZ_jets_bjets_recoW_recoT_selection.Define("ReweightedTopPt", UnweightedTopPt, {"Top_Pt"});

	std::cout << "in else, after d_TopReweighted_ee" << std::endl;

}




//For ME_Up and ME_Down
ints SummedWeights(14, 0);

auto NominalWeight{[&PDF_ScaleUp, &PDF_ScaleDown](const floats& LHEPdfWeight, const floats& LHEWeight_originalXWGTUP){

  std::cout << "inside NominalWeight" << std::endl;

  float PdfMin = 1.0;
  float PdfMax = 1.0;

  //For the min and max Pdf weights
  for(int i = 0; i < LHEPdfWeight.size(); i++){


        float LHEDivision = LHEPdfWeight.at(i) / LHEWeight_originalXWGTUP.at(0); //the size of LHEWeight_originalXWGTUP is always 1

        if(LHEDivision > PdfMax){PdfMax = LHEDivision;}
        else{continue;}

        if(LHEDivision < PdfMin){PdfMin = LHEDivision;}
        else{continue;}

  }


  if(PDF_ScaleUp == true){return PdfMax;}
  else if(PDF_ScaleDown == true){return PdfMin;}
  else{float One = 1.0; return One;}


}};




auto ME_uncert_function{[&SummedWeights](const floats& LHEPdfWeight, const floats& LHEWeight_originalXWGTUP, const floats& ReturnedPSWeight){

  std::cout << "inside ME_uncert_function" << std::endl;
  std::cout << "LHEPdfWeight.size() = " << LHEPdfWeight.size() << std::endl;
  std::cout << "LHEWeight_originalXWGTUP.size() = " << LHEWeight_originalXWGTUP.size() << std::endl;

  floats pdf = LHEPdfWeight / LHEWeight_originalXWGTUP.at(0);

  std::cout << "first pdf.size() = " << pdf.size() << std::endl;

  if(pdf.size() > 0){
  	for(int i = 0; i < pdf.size(); i++){ std::cout << "pdf.size() = " << pdf.size() << std::endl; pdf.at(i) >= 0.0 ? SummedWeights[0]++ : SummedWeights[1]++;} //pdf weight
  }

  std::cout << "ReturnedPSWeight.size() = " << ReturnedPSWeight.size() << std::endl;

  ReturnedPSWeight.at(1) >= 0.0 ? SummedWeights[2]++ : SummedWeights[3]++; //fsr down
  ReturnedPSWeight.at(0) >= 0.0 ? SummedWeights[4]++ : SummedWeights[5]++; //isr down
  (ReturnedPSWeight.at(1) * ReturnedPSWeight.at(0)) >= 0.0 ? SummedWeights[6]++ : SummedWeights[7]++; //both isr and fsr down
  ReturnedPSWeight.at(3) >= 0.0 ? SummedWeights[8]++ : SummedWeights[9]++; //fsr up
  ReturnedPSWeight.at(2) >= 0.0 ? SummedWeights[10]++ : SummedWeights[11]++; //isr up
  (ReturnedPSWeight.at(3) * ReturnedPSWeight.at(2)) >= 0.0 ? SummedWeights[12]++ : SummedWeights[13]++; //both isr and fsr up

  std::cout << "after ReturnedPSWeight" << std::endl;

  //SF is:  (total num of +ively-weighted events - total num of -ively-weighted events) / (total num of +ively-weighted events - total num of -ively-weighted events)

  int TotalNumPositive = SummedWeights[0] + SummedWeights[2] + SummedWeights[4] + SummedWeights[6] + SummedWeights[8] + SummedWeights[10] + SummedWeights[12]; 
  int TotalNumNegative = SummedWeights[1] + SummedWeights[3] + SummedWeights[5] + SummedWeights[7] + SummedWeights[9] + SummedWeights[11] + SummedWeights[13]; 


  float ME_SF = (TotalNumPositive + TotalNumNegative) / (TotalNumPositive - TotalNumNegative);

  std::cout << "after float ME_SF" << std::endl;

  return ME_SF;

}};



//Histogram for ME uncertainties
auto ME_histo_function{[&SummedWeights](){

  std::cout << "inside ME_histo_function" << std::endl;

  ints numerators;

  for(int i; i < SummedWeights.size(); i+=2){int output = SummedWeights[i] + SummedWeights[i+1]; std::cout << "output = " << output << std::endl; numerators.push_back(output);}

  return numerators;

}};


//SFs for ME up and down
auto GeneratorWeight{[&SummedWeights, &ME_Up, &ME_Down](const ints& ME_numerator_histo, const float& CalculatedNominalWeight, const floats& ReturnedPSWeight){

	std::cout << "inside GeneratorWeight" << std::endl;

 	int TotalNumPositive = SummedWeights[0] + SummedWeights[2] + SummedWeights[4] + SummedWeights[6] + SummedWeights[8] + SummedWeights[10] + SummedWeights[12];


	if(ME_Up == true){

		float generatorWeight_ScaleUp = (TotalNumPositive / ME_numerator_histo.at(7)) * ( (ReturnedPSWeight.at(3) * ReturnedPSWeight.at(2)) / abs(CalculatedNominalWeight) ); 
                return generatorWeight_ScaleUp; 
	
	}
	else if(ME_Down == true){	

		float generatorWeight_ScaleDown =  (TotalNumPositive / ME_numerator_histo.at(1)) * ( (ReturnedPSWeight.at(1) * ReturnedPSWeight.at(0)) / abs(CalculatedNominalWeight) ); 
		return generatorWeight_ScaleDown; 

	}
	else{	

		float generatorWeight = (TotalNumPositive / ME_numerator_histo.at(4)) * ( CalculatedNominalWeight / abs(CalculatedNominalWeight) );
		return generatorWeight; 

	}

}};





auto DummyColumnFunction{[](const floats& pts){

	return pts;

}};


std::string PSWeightString_ee;
std::string PSWeightString_mumu;

if( (year == "2017" || year == "2018") &&
     (process == "tZq" ||
     process == "SingleTop_tbarW" ||
     process == "SingleTop_schannel" ||
     process == "SingleTop_tchannel_top" ||
     process == "SingleTop_tchannel_tbar" ||
     process == "ttbarV_ttgamma" ||
     process == "ttbar_TTToHadronic" ||
     process == "ttbar_TTToSemileptonic") ){

	PSWeightString_ee = "PSWeight";
	PSWeightString_mumu = "PSWeight";

}
else{PSWeightString_ee = "Electron_pt_Selection"; PSWeightString_mumu = "Muon_pt_Selection";}



auto d_WeightedEvents_ee = d_TopReweighted_ee.Define("TotalHT_System", TotalHT_System, TotalHT_System_strings)
                                             .Define("TotalPt_System", TotalPt_System, TotalPt_System_strings)
					     .Define("TotHTOverTotpT_System", TotHTOverTotpT_floats, {"TotalHT_System", "TotalPt_System"})
					     .Define("DummyColumnBJet", DummyColumnFunction, {"bjetpt"})
					     .Define("CMSBTagSF", CMSBTagSF, {"bjetpt", "bjeteta", "Jet_btagCSVV2", "Jet_partonFlavour"})
					     .Define("nonbjets", nonbjet_id, {"tight_jets", "Jet_btagCSVV2", JetEtaInput})
                                             .Define("notbjetpt", bjet_variable, nonbjet_pt_strings)
                                             .Define("notbjeteta", bjet_variable, nonbjet_eta_strings)
  					     .Define("CMSNonBTagSF", CMSNonBTagSF, {"notbjetpt", "notbjeteta", "Jet_btagCSVV2", "Jet_partonFlavour"})
					     .Define("EffBTagged", EffBTaggedFunction_ee, {"DummyColumnBJet", JetPtInput, JetEtaInput})
			    		     .Define("EffNonBTagged", EffNonBTaggedFunction_ee, {"DummyColumnBJet", JetPtInput, JetEtaInput})
					     .Define("EffBTaggedProduct", EffBTaggedProduct, {"EffBTagged"})
					     .Define("EffNonBTaggedProduct", EffNonBTaggedProduct, {"EffNonBTagged"})
					     .Define("EffBTaggedProductData", EffBTaggedProductData, {"EffBTagged", "CMSBTagSF"})
                                             .Define("EffNonBTaggedProductData", EffNonBTaggedProductData, {"EffNonBTagged", "CMSNonBTagSF"})
					     .Define("ProbBTagMC", ProbBTagMCFunction, {"EffBTaggedProduct", "EffNonBTaggedProduct"})
 					     .Define("ProbBTagData", ProbBTagDataFunction, {"EffBTaggedProductData", "EffNonBTaggedProductData"})
					     .Define("BTagWeight", BTagWeightFunction, {"ProbBTagMC", "ProbBTagData"})
					     .Define("EGammaSF_egammaEff", EGammaSF_egammaEff, {"Electron_pt_Selection", "Electron_eta_Selection"})
					     .Define("EGammaSF_egammaEffSys", EGammaSF_egammaEff_Sys, {"Electron_pt_Selection", "Electron_eta_Selection"})
					     .Define("EGammaSF_egammaEffReco", EGammaSF_egammaEffReco, {"Electron_pt_Selection", "Electron_eta_Selection"})
					     .Define("EGammaSF_egammaEffRecoSys", EGammaSF_egammaEffReco_Sys, {"Electron_pt_Selection", "Electron_eta_Selection"})
					     .Define("ReturnedPSWeight", PSWeight, {PSWeightString_ee, "Electron_pt_Selection"})
					     .Define("CalculatedNominalWeight", NominalWeight, {"LHEPdfWeight", "LHEWeight_originalXWGTUP"})
					     .Define("ME_SF", ME_uncert_function, {"LHEPdfWeight", "LHEWeight_originalXWGTUP", "ReturnedPSWeight"})
					     .Define("ME_numerator_histo", ME_histo_function, {})
					     .Define("CalculatedGeneratorWeight", GeneratorWeight, {"ME_numerator_histo", "CalculatedNominalWeight", "ReturnedPSWeight"});
								      


auto d_WeightedEvents_mumu = d_TopReweighted_mumu.Define("TotalHT_System", TotalHT_System, TotalHT_System_strings)
                                                 .Define("TotalPt_System", TotalPt_System, TotalPt_System_strings)
				 		 .Define("TotHTOverTotpT_System", TotHTOverTotpT_floats, {"TotalHT_System", "TotalPt_System"})
						 .Define("DummyColumnBJet", DummyColumnFunction, {"bjetpt"})
						 .Define("CMSBTagSF", CMSNonBTagSF, {"bjetpt", "bjeteta", "Jet_btagCSVV2", "Jet_partonFlavour"})
                                                 .Define("nonbjets", nonbjet_id, {"tight_jets", "Jet_btagCSVV2", JetEtaInput})
                                                 .Define("notbjetpt", bjet_variable, nonbjet_pt_strings)
                                                 .Define("notbjeteta", bjet_variable, nonbjet_eta_strings)
                                                 .Define("CMSNonBTagSF", CMSNonBTagSF, {"notbjetpt", "notbjeteta", "Jet_btagCSVV2", "Jet_partonFlavour"})
                                                 .Define("EffBTagged", EffBTaggedFunction_mumu, {"DummyColumnBJet", JetPtInput, JetEtaInput})
                                                 .Define("EffNonBTagged", EffNonBTaggedFunction_mumu, {"DummyColumnBJet", JetPtInput, JetEtaInput})
                                                 .Define("EffBTaggedProduct", EffBTaggedProduct, {"EffBTagged"})
                                                 .Define("EffNonBTaggedProduct", EffNonBTaggedProduct, {"EffNonBTagged"})
                                                 .Define("EffBTaggedProductData", EffBTaggedProductData, {"EffBTagged", "CMSBTagSF"})
                                                 .Define("EffNonBTaggedProductData", EffNonBTaggedProductData, {"EffNonBTagged", "CMSNonBTagSF"})
                                                 .Define("ProbBTagMC", ProbBTagMCFunction, {"EffBTaggedProduct", "EffNonBTaggedProduct"})
                                                 .Define("ProbBTagData", ProbBTagDataFunction, {"EffBTaggedProductData", "EffNonBTaggedProductData"})
                                                 .Define("BTagWeight", BTagWeightFunction, {"ProbBTagMC", "ProbBTagData"})
						 .Define("MuonSFTest_ID", MuonSFTest_ID, {"MuonPt_RochCorr", "MuonEta_RochCorr"})
						 .Define("MuonSFTest_Iso", MuonSFTest_Iso, {"MuonPt_RochCorr", "MuonEta_RochCorr"})
                                                 .Define("MuonSFTest_ID_sys_syst", MuonSFTest_ID_sys_syst, {"MuonPt_RochCorr", "MuonEta_RochCorr"})
						 .Define("MuonSFTest_ID_sys_stat", MuonSFTest_ID_sys_stat, {"MuonPt_RochCorr", "MuonEta_RochCorr"})
						 .Define("MuonSFTest_Iso_sys_syst", MuonSFTest_Iso_sys_syst, {"MuonPt_RochCorr", "MuonEta_RochCorr"})
                                                 .Define("MuonSFTest_Iso_sys_stat", MuonSFTest_Iso_sys_stat, {"MuonPt_RochCorr", "MuonEta_RochCorr"})
						 .Define("ReturnedPSWeight", PSWeight, {PSWeightString_mumu, "MuonPt_RochCorr"})
					   	 .Define("CalculatedNominalWeight", NominalWeight, {"LHEPdfWeight", "LHEWeight_originalXWGTUP"})
						 .Define("ME_SF", ME_uncert_function, {"LHEPdfWeight", "LHEWeight_originalXWGTUP", "ReturnedPSWeight"})
						 .Define("ME_numerator_histo", ME_histo_function, {})
						 .Define("CalculatedGeneratorWeight", GeneratorWeight, {"ME_numerator_histo", "CalculatedNominalWeight", "ReturnedPSWeight"});
				


//lambda function for implementing the MET uncertainties
auto METUncertFunction{[&MET_Up, &MET_Down](

const floats& MET_MetUnclustEnUpDeltaX, 
const floats& MET_MetUnclustEnUpDeltaY, 
const floats& MET_phi,
const floats& MET_sumEt,
std::vector<TLorentzVector> SmearedJet4Momentum,
const floats& Jet_pt, 
const floats& Jet_eta, 
const floats& Jet_phi, 
const floats& Jet_mass){

  std::vector<TLorentzVector> metVecOriginal{};
  floats metVecOriginal_px;
  floats metVecOriginal_py;

  std::vector<TLorentzVector> metVec{};
  std::vector<TLorentzVector> UnsmearedJet{};
  floats SmearedJetPxVec;
  floats SmearedJetPyVec;
  floats UnsmearedJetPx;
  floats UnsmearedJetPy;

  std::cout << "inside METUncertFunction" << std::endl;

  //TLorentzVector for unsmeared jets
  for(int i = 0; i < Jet_pt.size(); i++){ ( UnsmearedJet.at(i) ).SetPtEtaPhiM(Jet_pt.at(i), Jet_eta.at(i), Jet_phi.at(i), Jet_mass.at(i)); }

  //Obtaining the px and py of unsmeared jets
  for(int i = 0; i < UnsmearedJet.size(); i++){ UnsmearedJetPx.push_back( (UnsmearedJet.at(i)).Px() ); }
  for(int i = 0; i < UnsmearedJet.size(); i++){ UnsmearedJetPy.push_back( (UnsmearedJet.at(i)).Py() ); }
   

  //Obtaining the px and py of smeared jets
  for(int i = 0; i < SmearedJet4Momentum.size(); i++){
 
  	float SmearedJetPx = ( SmearedJet4Momentum.at(i) ).Px();
	float SmearedJetPy = ( SmearedJet4Momentum.at(i) ).Py();
 	SmearedJetPxVec.push_back(SmearedJetPx);
	SmearedJetPyVec.push_back(SmearedJetPy);

  }

  //Original MET vector
  for(int i = 0; i < MET_phi.size(); i++){ 

	(metVecOriginal.at(i)).SetPtEtaPhiE(MET_sumEt.at(i), 0, MET_phi.at(i), MET_sumEt.at(i)); 
	metVecOriginal_px.push_back( (metVecOriginal.at(i)).Px() );
	metVecOriginal_py.push_back( (metVecOriginal.at(i)).Py() );

  }

  floats MET_px_up =  metVecOriginal_px + MET_MetUnclustEnUpDeltaX;
  floats MET_py_up =  metVecOriginal_py + MET_MetUnclustEnUpDeltaY;
  floats MET_px_down =  metVecOriginal_px - MET_MetUnclustEnUpDeltaX;
  floats MET_py_down =  metVecOriginal_py - MET_MetUnclustEnUpDeltaY;  

  //For the nominal MET and MET uncertainties
  
  floats UnclusteredEnergyUp = sqrt( pow(MET_px_up, 2) + pow(MET_py_up, 2) );
  floats UnclusteredEnergyDown = sqrt( pow(MET_px_down, 2) + pow(MET_py_down, 2) );

  for(int i = 0; i < MET_phi.size(); i++){

  	if(MET_Up == true){ (metVec.at(i)).SetPtEtaPhiE(UnclusteredEnergyUp.at(i), 0, MET_phi.at(i), UnclusteredEnergyUp.at(i));}
  	else if(MET_Down == true){ (metVec.at(i)).SetPtEtaPhiE(UnclusteredEnergyDown.at(i), 0, MET_phi.at(i), UnclusteredEnergyDown.at(i));}
  	else{ (metVec.at(i)).SetPtEtaPhiE(MET_sumEt.at(i), 0, MET_phi.at(i), MET_sumEt.at(i));}

 }

 //Propagating the jet smearing to the MET
 
 for(int i = 0; i < SmearedJetPxVec.size(); i++){
 
 	( metVec.at(i) ).SetPx( (metVec.at(i)).Px() + UnsmearedJetPx.at(i));
        ( metVec.at(i) ).SetPy( (metVec.at(i)).Py() + UnsmearedJetPy.at(i));
 	( metVec.at(i) ).SetPx( (metVec.at(i)).Px() - SmearedJetPxVec.at(i));
        ( metVec.at(i) ).SetPy( (metVec.at(i)).Py() - SmearedJetPyVec.at(i));
 
 }


  return metVec;

}};




std::vector<std::string> MET_uncert_strings = {

"MET_MetUnclustEnUpDeltaX",
"MET_MetUnclustEnUpDeltaY",
"MET_phi",
"MET_sumEt",
"SmearedJet4Momentum",
"Jet_pt",
"Jet_eta",
"Jet_phi",
"Jet_mass"

};


//Defining the new MET column
auto d_WeightedEvents_withMET_ee = d_WeightedEvents_ee.Define("newMET", METUncertFunction, MET_uncert_strings);
auto d_WeightedEvents_withMET_mumu = d_WeightedEvents_mumu.Define("newMET", METUncertFunction, MET_uncert_strings);


int nbins = 40;


//Chi^2 calculation using MC samples

std::string Filename;


if(process == "tZq"){


  if(PU_ScaleUp == true){

	if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_PU_ScaleUp_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_ZPlusJetsCR_PU_ScaleUp_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_ttbarCR_PU_ScaleUp_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ZPlusJetsCR_PU_ScaleUp_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ttbarCR_PU_ScaleUp_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
        else{Filename = process + "_AfterFullSelection_GaussianFit_PU_ScaleUp_" + year + ".root";}


  }
  else if(PU_ScaleDown == true){

        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_PU_ScaleDown_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_ZPlusJetsCR_PU_ScaleDown_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_ttbarCR_PU_ScaleDown_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ZPlusJetsCR_PU_ScaleDown_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ttbarCR_PU_ScaleDown_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
        else{Filename = process + "_AfterFullSelection_GaussianFit_PU_ScaleDown_" + year + ".root";}


  }
  else if(BTag_ScaleUp == true){

        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_BTag_ScaleUp_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_ZPlusJetsCR_BTag_ScaleUp_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_ttbarCR_BTag_ScaleUp_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ZPlusJetsCR_BTag_ScaleUp_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ttbarCR_BTag_ScaleUp_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
        else{Filename = process + "_AfterFullSelection_GaussianFit_BTag_ScaleUp_" + year + ".root";}


  }
  else if(BTag_ScaleDown == true){

        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_BTag_ScaleDown_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_ZPlusJetsCR_BTag_ScaleDown_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_ttbarCR_BTag_ScaleDown_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ZPlusJetsCR_BTag_ScaleDown_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ttbarCR_BTag_ScaleDown_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
        else{Filename = process + "_AfterFullSelection_GaussianFit_BTag_ScaleDown_" + year + ".root";}


  }
  else if(JetSmearing_ScaleUp == true){

        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_JetSmearing_ScaleUp_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_ZPlusJetsCR_JetSmearing_ScaleUp_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_ttbarCR_JetSmearing_ScaleUp_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ZPlusJetsCR_JetSmearing_ScaleUp_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ttbarCR_JetSmearing_ScaleUp_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
        else{Filename = process + "_AfterFullSelection_GaussianFit_JetSmearing_ScaleUp_" + year + ".root";}


  }
  else if(JetSmearing_ScaleDown == true){

        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_JetSmearing_ScaleDown_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_ZPlusJetsCR_JetSmearing_ScaleDown_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_ttbarCR_JetSmearing_ScaleDown_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ZPlusJetsCR_JetSmearing_ScaleDown_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ttbarCR_JetSmearing_ScaleDown_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
        else{Filename = process + "_AfterFullSelection_GaussianFit_JetSmearing_ScaleDown_" + year + ".root";}


  }
  else if(JetResolution_ScaleUp == true){
  
        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_JetResolution_ScaleUp_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_ZPlusJetsCR_JetResolution_ScaleUp_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_ttbarCR_JetResolution_ScaleUp_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ZPlusJetsCR_JetResolution_ScaleUp_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ttbarCR_JetResolution_ScaleUp_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
        else{Filename = process + "_AfterFullSelection_GaussianFit_JetResolution_ScaleUp_" + year + ".root";}


  }
  else if(JetResolution_ScaleDown == true){

        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_JetResolution_ScaleDown_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_ZPlusJetsCR_JetResolution_ScaleDown_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_ttbarCR_JetResolution_ScaleDown_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ZPlusJetsCR_JetResolution_ScaleDown_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ttbarCR_JetResolution_ScaleDown_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
        else{Filename = process + "_AfterFullSelection_GaussianFit_JetResolution_ScaleDown_" + year + ".root";}


  }
  else if(LeptonEfficiencies_ScaleUp == true){

        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_LeptonEfficiencies_ScaleUp_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_ZPlusJetsCR_LeptonEfficiencies_ScaleUp_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_ttbarCR_LeptonEfficiencies_ScaleUp_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ZPlusJetsCR_LeptonEfficiencies_ScaleUp_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ttbarCR_LeptonEfficiencies_ScaleUp_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
        else{Filename = process + "_AfterFullSelection_GaussianFit_LeptonEfficiencies_ScaleUp_" + year + ".root";}


  }
  else if(LeptonEfficiencies_ScaleDown == true){

        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_LeptonEfficiencies_ScaleDown_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_ZPlusJetsCR_LeptonEfficiencies_ScaleDown_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_ttbarCR_LeptonEfficiencies_ScaleDown_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ZPlusJetsCR_LeptonEfficiencies_ScaleDown_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
               Filename = process + "_AfterFullSelection_GaussianFit_NPL_ttbarCR_LeptonEfficiencies_ScaleDown_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
        else{Filename = process + "_AfterFullSelection_GaussianFit_LeptonEfficiencies_ScaleDown_" + year + ".root";}


  }
  else if(PDF_ScaleUp == true){
  
        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_PDF_ScaleUp_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_ZPlusJetsCR_PDF_ScaleUp_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_ttbarCR_PDF_ScaleUp_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ZPlusJetsCR_PDF_ScaleUp_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ttbarCR_PDF_ScaleUp_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
        else{Filename = process + "_AfterFullSelection_GaussianFit_PDF_ScaleUp_" + year + ".root";}


  }
  else if(PDF_ScaleDown == true){

        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_PDF_ScaleDown_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_ZPlusJetsCR_PDF_ScaleDown_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_ttbarCR_PDF_ScaleDown_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ZPlusJetsCR_PDF_ScaleDown_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ttbarCR_PDF_ScaleDown_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
        else{Filename = process + "_AfterFullSelection_GaussianFit_PDF_ScaleDown_" + year + ".root";}


  }
  else if(ME_Up == true){
  
        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ME_Up_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_ZPlusJetsCR_ME_Up_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_ttbarCR_ME_Up_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ZPlusJetsCR_ME_Up_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ttbarCR_PDF_ME_Up_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
        else{Filename = process + "_AfterFullSelection_GaussianFit_ME_Up_" + year + ".root";}


  }
  else if(ME_Down == true){

        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ME_Down_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_ZPlusJetsCR_ME_Down_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_ttbarCR_ME_Down_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ZPlusJetsCR_ME_Down_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ttbarCR_ME_Down_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
        else{Filename = process + "_AfterFullSelection_GaussianFit_ME_Down_" + year + ".root";}


  }
  else if(MET_Up == true){
  
        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_MET_Up_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_ZPlusJetsCR_MET_Up_" + year + ".root";

        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_ttbarCR_PDF_MET_Up_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ZPlusJetsCR_MET_Up_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ttbarCR_MET_Up_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
        else{Filename = process + "_AfterFullSelection_GaussianFit_MET_Up_" + year + ".root";}


  }
  else if(MET_Down == true){

        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_MET_Down_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_ZPlusJetsCR_MET_Down_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_ttbarCR_MET_Down_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ZPlusJetsCR_MET_Down_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ttbarCR_MET_Down_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
        else{Filename = process + "_AfterFullSelection_GaussianFit_MET_Down_" + year + ".root";}


  }
  else if(isr_up == true){
  
        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_isr_up_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_ZPlusJetsCR_isr_up_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_ttbarCR_isr_up_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ZPlusJetsCR_isr_up_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ttbarCR_isr_up_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
        else{Filename = process + "_AfterFullSelection_GaussianFit_isr_up_" + year + ".root";}


  }
  else if(isr_down == true){

        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_isr_down_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_ZPlusJetsCR_isr_down_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_ttbarCR_isr_down_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ZPlusJetsCR_isr_down_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ttbarCR_isr_down_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
        else{Filename = process + "_AfterFullSelection_GaussianFit_isr_down_" + year + ".root";}


  }
  else if(fsr_up == true){

        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_fsr_up_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_ZPlusJetsCR_fsr_up_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_ttbarCR_fsr_up_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ZPlusJetsCR_fsr_up_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ttbarCR_fsr_up_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
        else{Filename = process + "_AfterFullSelection_GaussianFit_fsr_up_" + year + ".root";}


  }
  else if(fsr_down == true){

        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_fsr_down_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_ZPlusJetsCR_fsr_down_" + year + ".root";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_ttbarCR_fsr_down_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ZPlusJetsCR_fsr_down_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                Filename = process + "_AfterFullSelection_GaussianFit_NPL_ttbarCR_fsr_down_" + year + ".root";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
        else{Filename = process + "_AfterFullSelection_GaussianFit_fsr_down_" + year + ".root";}


  }
  else{

	std::cout << "inside the else statement for nominal" << std::endl;

	if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
        	Filename = process + "_AfterFullSelection_GaussianFit_NPL_" + year + ".root";
	}
	else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
		Filename = process + "_AfterFullSelection_GaussianFit_ZPlusJetsCR_" + year + ".root";
	}
	else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
		Filename = process + "_AfterFullSelection_GaussianFit_ttbarCR_" + year + ".root";
	}
	else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
		Filename = process + "_AfterFullSelection_GaussianFit_NPL_ZPlusJetsCR_" + year + ".root";
	}
	else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
		Filename = process + "_AfterFullSelection_GaussianFit_NPL_ttbarCR_" + year + ".root";
	}
	else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
	else{Filename = process + "_AfterFullSelection_GaussianFit_" + year + ".root";}

  }



	TFile* FittedHistosOutput = new TFile(Filename.c_str(), "RECREATE");
        FittedHistosOutput->cd();


	auto h_WMass_ee = d_WeightedEvents_withMET_ee.Histo1D({"h_WMass_ee", "Mass distribution of the W mass candidate (ee channel)", nbins, 0, 150}, "w_mass");
	auto h_InvTopMass_ee = d_WeightedEvents_withMET_ee.Histo1D({"h_InvTopMass_ee", "Mass distribution of the top candidate (ee channel)", nbins, 0, 500}, "InvTopMass");
	auto h_WMass_mumu = d_WeightedEvents_withMET_mumu.Histo1D({"h_WMass_mumu", "Mass distribution of the W mass candidate (mumu channel)", nbins, 0, 150}, "w_mass");
	auto h_InvTopMass_mumu = d_WeightedEvents_withMET_mumu.Histo1D({"h_InvTopMass_mumu", "Mass distribution of the top candidate (mumu channel)", nbins, 0, 500}, "InvTopMass");

	h_WMass_ee->Fit("gaus");
	h_InvTopMass_ee->Fit("gaus");
	h_WMass_ee->Write();
	h_InvTopMass_ee->Write();

	h_WMass_mumu->Fit("gaus");
        h_InvTopMass_mumu->Fit("gaus");
        h_WMass_mumu->Write();
        h_InvTopMass_mumu->Write();

	W_stddev_ee = h_WMass_ee->GetStdDev(); //Finding the resolution (same as the standard deviation)
	Top_stddev_ee = h_InvTopMass_ee->GetStdDev(); //Finding the resolution (same as the standard deviation)
	W_stddev_mumu = h_WMass_mumu->GetStdDev(); //Finding the resolution (same as the standard deviation)
        Top_stddev_mumu = h_InvTopMass_mumu->GetStdDev(); //Finding the resolution (same as the standard deviation)

	FittedHistosOutput->Close();


	std::string tZq_WAndTop_Filename; 

	
	if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
		tZq_WAndTop_Filename = "tZq_AfterFullSelection_GaussianFit_NPL_" + year + ".root";
	}
	else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
		tZq_WAndTop_Filename = "tZq_AfterFullSelection_GaussianFit_ZPlusJetsCR_" + year + ".root";
	}
	else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
		tZq_WAndTop_Filename = "tZq_AfterFullSelection_GaussianFit_ttbarCR_" + year + ".root";
	}
	else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){ 
	 	tZq_WAndTop_Filename = "tZq_AfterFullSelection_GaussianFit_NPL_ZPlusJetsCR_" + year + ".root";
	}
	else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){ 
	 	tZq_WAndTop_Filename = "tZq_AfterFullSelection_GaussianFit_NPL_ttbarCR_" + year + ".root";
	}
	else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
	else{tZq_WAndTop_Filename = "tZq_AfterFullSelection_GaussianFit_" + year + ".root";}





	TFile* tZq_WAndTop_File = new TFile{tZq_WAndTop_Filename.c_str(), "READ"};

	TH1* WHist_ee = (TH1*)tZq_WAndTop_File->GetObjectChecked("h_WMass_ee", "TH1");
	TH1* TopHist_ee = (TH1*)tZq_WAndTop_File->GetObjectChecked("h_InvTopMass_ee", "TH1");
	TH1* WHist_mumu = (TH1*)tZq_WAndTop_File->GetObjectChecked("h_WMass_mumu", "TH1");
	TH1* TopHist_mumu = (TH1*)tZq_WAndTop_File->GetObjectChecked("h_InvTopMass_mumu", "TH1");

	W_stddev_ee = WHist_ee->GetStdDev();
	Top_stddev_ee = TopHist_ee->GetStdDev();
	W_stddev_mumu = WHist_mumu->GetStdDev();
	Top_stddev_mumu = TopHist_mumu->GetStdDev();

	tZq_WAndTop_File->Close();



	//Write the nominal mass and resolution values to a text file

	std::string filenamestring_variable;

	if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
		filenamestring_variable = "Resolution_" + process + "_" + year + "_NPL.txt";
	}
	else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
		filenamestring_variable = "Resolution_" + process + "_" + year + "_ZPlusJetsCR.txt";
	}
	else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
		filenamestring_variable = "Resolution_" + process + "_" + year + "_ttbarCR.txt";
	}
	else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
		filenamestring_variable = "Resolution_" + process + "_" + year + "_NPL_ZPlusJetsCR.txt";
	}
	else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
		filenamestring_variable = "Resolution_" + process + "_" + year + "_NPL_ttbarCR.txt";	
	}
	else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
	else{filenamestring_variable = "Resolution_" + process + "_" + year + ".txt";}



	std::ofstream Resolution;
	Resolution.open(filenamestring_variable.c_str());

	Resolution << "W_stddev_ee: " << W_stddev_ee << '\n'
	   	   << "Top_stddev_ee: " << Top_stddev_ee << '\n' 
	   	   << "W_stddev_mumu: " << W_stddev_mumu << '\n'
           	   << "Top_stddev_mumu: " << Top_stddev_mumu << '\n' << std::endl;



} //end of if statement for tZq

std::vector<float> CutRanges_ee = {};

//Lambda function for chi squared calculation (calculated using MC but applied to both MC and data)
auto chi2_ee{[&process, &CutRanges_ee, &SBR](const float& w_mass, const float& Top_Mass){
	

  //5 sigma is 99.99994%
  float FiveSigma = 0.9999994;
  float FiveSigmaOverTwo = FiveSigma/2;

  float LowerBound = W_MASS - (W_MASS * (1 - FiveSigmaOverTwo));
  float UpperBound = W_MASS + (W_MASS * (1 - FiveSigmaOverTwo));

  //calculating chi2
  float chi2 = pow(( (w_mass - W_MASS) / W_stddev_ee), 2) + pow(( (Top_Mass - TOP_MASS) / Top_stddev_ee), 2);


  if(process == "tZq" && SBR == true){

  	//returning chi2 values only for when w_mass is within 5 sigma of the known W mass 
  	if(w_mass > LowerBound && w_mass < UpperBound){
		CutRanges_ee.push_back(chi2);
		return chi2;
	}	
	else{
		std::cout << "w_mass is not within 5 sigma of the mean W mass value (ee)" << std::endl;
		float out = 999.0;
                return out;
	}

  }
  else{return chi2;}


}};


std::vector<float> CutRanges_mumu = {};

auto chi2_mumu{[&process, &CutRanges_mumu, &SBR](const float& w_mass, const float& Top_Mass){


  //5 sigma is 99.99994%
  float FiveSigma = 0.9999994;
  float FiveSigmaOverTwo = FiveSigma/2;

  float LowerBound = W_MASS - (W_MASS * (1 - FiveSigmaOverTwo));
  float UpperBound = W_MASS + (W_MASS * (1 - FiveSigmaOverTwo));


  //calculating chi2
  float chi2 = pow(( (w_mass - W_MASS) / W_stddev_mumu), 2) + pow(( (Top_Mass - TOP_MASS) / Top_stddev_mumu), 2);


  if(process == "tZq" && SBR == true){
  
	//returning chi2 values only for when w_mass is within 5 sigma of the known W mass 
        if(w_mass > LowerBound && w_mass < UpperBound){
		CutRanges_mumu.push_back(chi2);
		return chi2;
	}
        else{std::cout << "w_mass is not within 5 sigma of the mean W mass value (mumu)" << std::endl;
	     float out = 999.0;
	     return out;
	}


  }
  else{return chi2;}

 
}};





std::string Chi2Range_string;


//Section for experimental blinding
if(blinding == true && (SBR == true || SR == true)){


	auto Blinding_ee =  d_WeightedEvents_withMET_ee.Define("chi2", chi2_ee, {"w_mass", "InvTopMass"});
	auto Blinding_mumu =  d_WeightedEvents_withMET_mumu.Define("chi2", chi2_mumu, {"w_mass", "InvTopMass"});


	if(process == "tZq"){
		NumberOfSimulatedEvents_ee = *( Blinding_ee.Filter("chi2").Count() );
		NumberOfSimulatedEvents_mumu = *( Blinding_mumu.Filter("chi2").Count() );	
	}

	int OneSigmaOfNumEvents_ee = NumberOfSimulatedEvents_ee * 0.68;
	int OneSigmaOfNumEvents_mumu = NumberOfSimulatedEvents_mumu * 0.68;

	auto histo_chi2_ee = Blinding_ee.Histo1D("chi2");
	auto histo_chi2_mumu = Blinding_mumu.Histo1D("chi2");

	TAxis *xaxis_ee = histo_chi2_ee->GetXaxis();
	TAxis *xaxis_mumu = histo_chi2_mumu->GetXaxis();

	double MaxBin_ee = xaxis_ee->GetBinCenter( histo_chi2_ee->FindLastBinAbove() );
	double MinBin_ee = xaxis_ee->GetBinCenter( histo_chi2_ee->FindFirstBinAbove() ) ;
	double MaxBin_mumu = xaxis_mumu->GetBinCenter( histo_chi2_mumu->FindLastBinAbove() );
        double MinBin_mumu = xaxis_mumu->GetBinCenter( histo_chi2_mumu->FindFirstBinAbove() );

	int NumBins_ee = MaxBin_ee - MinBin_ee;
	int NumBins_mumu = MaxBin_mumu - MinBin_mumu;


	auto histo_chi2_ee_rebinned = Blinding_ee.Histo1D({"histo_chi2_ee_rebinned", "histo_chi2_ee_rebinned", 2*NumBins_ee, MinBin_ee, MaxBin_ee}, {"chi2"});
	auto histo_chi2_mumu_rebinned = Blinding_mumu.Histo1D({"histo_chi2_mumu_rebinned", "histo_chi2_mumu_rebinned", 2*NumBins_mumu, MinBin_mumu, MaxBin_mumu}, {"chi2"});

	TAxis * histo_chi2_ee_rebinned_x = histo_chi2_ee_rebinned->GetXaxis();
	TAxis * histo_chi2_mumu_rebinned_x = histo_chi2_mumu_rebinned->GetXaxis();

	int total_ee = 0;
	int total_mumu = 0;

	for(int i = 0; i < histo_chi2_ee_rebinned->GetEntries(); i++){

		auto NumberOfEvents_ee = histo_chi2_ee_rebinned->GetBinContent(i);
		total_ee += NumberOfEvents_ee;

		if(total_ee >= OneSigmaOfNumEvents_ee){Chi2_SR_ee = histo_chi2_ee_rebinned_x->GetBinCenter(i); break;}
		else{continue;}

	}

	for(int i = 0; i < histo_chi2_mumu_rebinned->GetEntries(); i++){

                auto NumberOfEvents_mumu = histo_chi2_mumu_rebinned->GetBinContent(i);
                total_mumu += NumberOfEvents_mumu;
		
                if(total_mumu >= OneSigmaOfNumEvents_mumu){Chi2_SR_mumu = histo_chi2_mumu_rebinned_x->GetBinCenter(i); break;}
                else{continue;}

        }


	//chi2 range for ee

	auto Chi2Cut_ee{[&SBR, &SR](const float& Chi2){	

	  if(SBR == true){return Chi2_SR_ee < Chi2 && Chi2 < Chi2_SBR_ee;}
	  else if(SR == true){return Chi2 < Chi2_SR_ee;}
	  else{std::cout << "SB and SR cannot both be false" << std::endl;}

	}};

	auto AfterChi2Cut_ee = Blinding_ee.Define("AfterChi2Cut_ee", Chi2Cut_ee, {"chi2"}).Filter(Chi2Cut_ee, {"chi2"});

	//chi2 range for mumu

        auto Chi2Cut_mumu{[&SBR, &SR](const float& Chi2){

          if(SBR == true){return Chi2_SR_mumu < Chi2 && Chi2 < Chi2_SBR_mumu;}
          else if(SR == true){return Chi2 < Chi2_SR_mumu;}
          else{std::cout << "SB and SR cannot both be false" << std::endl;}

        }};

        auto AfterChi2Cut_mumu = Blinding_mumu.Define("AfterChi2Cut_ee", Chi2Cut_mumu, {"chi2"}).Filter(Chi2Cut_mumu, {"chi2"});	
		

	



	//Naming the chi2 text file
	if(PU_ScaleUp == true){

		if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_PU_ScaleUp_" + year + "_NPL.txt";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_PU_ScaleUp_" + year + "_ZPlusJetsCR.txt";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_PU_ScaleUp_" + year + "_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_PU_ScaleUp_" + year + "_NPL_ZPlusJetsCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){	
                        Chi2Range_string = "Chi2Range_" + process + "_PU_ScaleUp_" + year + "_NPL_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{Chi2Range_string = "Chi2Range_" + process + "_PU_ScaleUp_" + year + ".txt";}


	}
	else if(PU_ScaleDown == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_PU_ScaleDown_" + year + "_NPL.txt";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_PU_ScaleDown_" + year + "_ZPlusJetsCR.txt";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_PU_ScaleDown_" + year + "_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_PU_ScaleDown_" + year + "_NPL_ZPlusJetsCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_PU_ScaleDown_" + year + "_NPL_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{Chi2Range_string = "Chi2Range_" + process + "_PU_ScaleDown_" + year + ".txt";}


        }
	else if(BTag_ScaleUp == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_BTag_ScaleUp_" + year + "_NPL.txt";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_BTag_ScaleUp_" + year + "_ZPlusJetsCR.txt";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_BTag_ScaleUp_" + year + "_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_BTag_ScaleUp_" + year + "_NPL_ZPlusJetsCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_BTag_ScaleUp_" + year + "_NPL_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{Chi2Range_string = "Chi2Range_" + process + "_BTag_ScaleUp_" + year + ".txt";}


        }
        else if(BTag_ScaleDown == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_BTag_ScaleDown_" + year + "_NPL.txt";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_BTag_ScaleDown_" + year + "_ZPlusJetsCR.txt";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_BTag_ScaleDown_" + year + "_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_BTag_ScaleDown_" + year + "_NPL_ZPlusJetsCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_BTag_ScaleDown_" + year + "_NPL_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{Chi2Range_string = "Chi2Range_" + process + "_BTag_ScaleDown_" + year + ".txt";}


        }
	else if(JetSmearing_ScaleUp == true){
        
                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_JetSmearing_ScaleUp_" + year + "_NPL.txt";
                }       
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_JetSmearing_ScaleUp_" + year + "_ZPlusJetsCR.txt";
                }       
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_JetSmearing_ScaleUp_" + year + "_ttbarCR.txt";
                }       
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_JetSmearing_ScaleUp_" + year + "_NPL_ZPlusJetsCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_JetSmearing_ScaleUp_" + year + "_NPL_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{Chi2Range_string = "Chi2Range_" + process + "_JetSmearing_ScaleUp_" + year + ".txt";}


        }
        else if(JetSmearing_ScaleDown == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_JetSmearing_ScaleDown_" + year + "_NPL.txt";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_JetSmearing_ScaleDown_" + year + "_ZPlusJetsCR.txt";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_JetSmearing_ScaleDown_" + year + "_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_JetSmearing_ScaleDown_" + year + "_NPL_ZPlusJetsCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_JetSmearing_ScaleDown_" + year + "_NPL_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{Chi2Range_string = "Chi2Range_" + process + "_JetSmearing_ScaleDown_" + year + ".txt";}


        }
	else if(JetResolution_ScaleUp == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_JetResolution_ScaleUp_" + year + "_NPL.txt";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_JetResolution_ScaleUp_" + year + "_ZPlusJetsCR.txt";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_JetResolution_ScaleUp_" + year + "_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_JetResolution_ScaleUp_" + year + "_NPL_ZPlusJetsCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_JetResolution_ScaleUp_" + year + "_NPL_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{Chi2Range_string = "Chi2Range_" + process + "_JetResolution_ScaleUp_" + year + ".txt";}


        }
        else if(JetResolution_ScaleDown == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_JetResolution_ScaleDown_" + year + "_NPL.txt";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_JetResolution_ScaleDown_" + year + "_ZPlusJetsCR.txt";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_JetResolution_ScaleDown_" + year + "_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_JetResolution_ScaleDown_" + year + "_NPL_ZPlusJetsCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_JetResolution_ScaleDown_" + year + "_NPL_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{Chi2Range_string = "Chi2Range_" + process + "_JetResolution_ScaleDown_" + year + ".txt";}


        }
	else if(LeptonEfficiencies_ScaleUp == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_LeptonEfficiencies_ScaleUp_" + year + "_NPL.txt";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_LeptonEfficiencies_ScaleUp_" + year + "_ZPlusJetsCR.txt";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_LeptonEfficiencies_ScaleUp_" + year + "_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_LeptonEfficiencies_ScaleUp_" + year + "_NPL_ZPlusJetsCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_LeptonEfficiencies_ScaleUp_" + year + "_NPL_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{Chi2Range_string = "Chi2Range_" + process + "_LeptonEfficiencies_ScaleUp_" + year + ".txt";}


        }
        else if(LeptonEfficiencies_ScaleDown == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_LeptonEfficiencies_ScaleDown_" + year + "_NPL.txt";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_LeptonEfficiencies_ScaleDown_" + year + "_ZPlusJetsCR.txt";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_LeptonEfficiencies_ScaleDown_" + year + "_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_LeptonEfficiencies_ScaleDown_" + year + "_NPL_ZPlusJetsCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_LeptonEfficiencies_ScaleDown_" + year + "_NPL_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{Chi2Range_string = "Chi2Range_" + process + "_JetResolution_ScaleDown_" + year + ".txt";}


        }
	else if(PDF_ScaleUp == true){
                
                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_PDF_ScaleUp_" + year + "_NPL.txt";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_PDF_ScaleUp_" + year + "_ZPlusJetsCR.txt";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_PDF_ScaleUp_" + year + "_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_PDF_ScaleUp_" + year + "_NPL_ZPlusJetsCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_PDF_ScaleUp_" + year + "_NPL_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{Chi2Range_string = "Chi2Range_" + process + "_PDF_ScaleUp_" + year + ".txt";}

        
        }
        else if(PDF_ScaleDown == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_PDF_ScaleDown_" + year + "_NPL.txt";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_PDF_ScaleDown_" + year + "_ZPlusJetsCR.txt";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_PDF_ScaleDown_" + year + "_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_PDF_ScaleDown_" + year + "_NPL_ZPlusJetsCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_PDF_ScaleDown_" + year + "_NPL_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{Chi2Range_string = "Chi2Range_" + process + "_PDF_ScaleDown_" + year + ".txt";}


        }
	else if(ME_Up == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_ME_Up_" + year + "_NPL.txt";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_ME_Up_" + year + "_ZPlusJetsCR.txt";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_ME_Up_" + year + "_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_ME_Up_" + year + "_NPL_ZPlusJetsCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_ME_Up_" + year + "_NPL_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{Chi2Range_string = "Chi2Range_" + process + "_ME_Up_" + year + ".txt";}
                
                
        }
        else if(ME_Down == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_ME_Down_" + year + "_NPL.txt";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_ME_Down_" + year + "_ZPlusJetsCR.txt";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_ME_Down_" + year + "_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_ME_Down_" + year + "_NPL_ZPlusJetsCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_ME_Down_" + year + "_NPL_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{Chi2Range_string = "Chi2Range_" + process + "_ME_Down_" + year + ".txt";}


        }
	else if(MET_Up == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_MET_Up_" + year + "_NPL.txt";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_MET_Up_" + year + "_ZPlusJetsCR.txt";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_MET_Up_" + year + "_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_MET_Up_" + year + "_NPL_ZPlusJetsCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_MET_Up_" + year + "_NPL_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{Chi2Range_string = "Chi2Range_" + process + "_MET_Up_" + year + ".txt";}


        }
        else if(MET_Down == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_MET_Down_" + year + "_NPL.txt";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_MET_Down_" + year + "_ZPlusJetsCR.txt";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_MET_Down_" + year + "_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_MET_Down_" + year + "_NPL_ZPlusJetsCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_MET_Down_" + year + "_NPL_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{Chi2Range_string = "Chi2Range_" + process + "_MET_Down_" + year + ".txt";}


        }
	else if(isr_up == true){
        
                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_isr_up_" + year + "_NPL.txt";
                }       
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_isr_up_" + year + "_ZPlusJetsCR.txt";
                }       
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_isr_up_" + year + "_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_isr_up_" + year + "_NPL_ZPlusJetsCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_isr_up_" + year + "_NPL_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{Chi2Range_string = "Chi2Range_" + process + "_isr_up_" + year + ".txt";}


        }
        else if(isr_down == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_isr_down_" + year + "_NPL.txt";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_isr_down_" + year + "_ZPlusJetsCR.txt";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_isr_down_" + year + "_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_isr_down_" + year + "_NPL_ZPlusJetsCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_isr_down_" + year + "_NPL_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{Chi2Range_string = "Chi2Range_" + process + "_isr_down_" + year + ".txt";}


        }
	else if(fsr_up == true){

                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_fsr_up_" + year + "_NPL.txt";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_fsr_up_" + year + "_ZPlusJetsCR.txt";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_fsr_up_" + year + "_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_fsr_up_" + year + "_NPL_ZPlusJetsCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_fsr_up_" + year + "_NPL_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{Chi2Range_string = "Chi2Range_" + process + "_fsr_up_" + year + ".txt";}
                        
                
        }       
        else if(fsr_down == true){
                
                if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_fsr_down_" + year + "_NPL.txt";
                }
                else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_fsr_down_" + year + "_ZPlusJetsCR.txt";
                }
                else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_fsr_down_" + year + "_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                        Chi2Range_string = "Chi2Range_" + process + "_fsr_down_" + year + "_NPL_ZPlusJetsCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                        Chi2Range_string = "Chi2Range_" + process + "_fsr_down_" + year + "_NPL_ttbarCR.txt";
                }
                else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
                else{Chi2Range_string = "Chi2Range_" + process + "_fsr_down_" + year + ".txt";}


        }
	else{

		if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
      			Chi2Range_string = "Chi2Range_" + process + "_" + year + "_NPL.txt";
		}
		else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
      			Chi2Range_string = "Chi2Range_" + process + "_" + year + "_ZPlusJetsCR.txt";
		}
		else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
			Chi2Range_string = "Chi2Range_" + process + "_" + year + "_ttbarCR.txt";
		}
		else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
			Chi2Range_string = "Chi2Range_" + process + "_" + year + "_NPL_ZPlusJetsCR.txt";
		}
		else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
			Chi2Range_string = "Chi2Range_" + process + "_" + year + "_NPL_ttbarCR.txt";
		}
		else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
		else{Chi2Range_string = "Chi2Range_" + process + "_" + year + ".txt";}

	}

	

	std::ofstream Chi2Range;
	Chi2Range.open(Chi2Range_string.c_str());

	Chi2Range << "Chi2_SR_ee: " << Chi2_SR_ee << '\n'
		  << "Chi2_SBR_ee: " << Chi2_SBR_ee << '\n'
		  << "Chi2_SR_mumu: " << Chi2_SR_mumu << '\n'
                  << "Chi2_SBR_mumu: " << Chi2_SBR_mumu << std::endl;



	//Saving the histograms to output root files
	std::string OutRootFile_ee;
	std::string OutRootFile_mumu;  

        if(NPL == false){
  
      		OutRootFile_ee  = "Results_MC_" + process + "_" + year + "_ee_Blinded.root";
        	OutRootFile_mumu = "Results_MC_" + process + "_" + year + "_mumu_Blinded.root";
	
	}
 	else{

		OutRootFile_ee  = "Results_MC_" + process + "_" + year + "_ee_NPL_Blinded.root";
                OutRootFile_mumu = "Results_MC_" + process + "_" + year + "_mumu_NPL_Blinded.root";

	}


	auto colNames_ee = AfterChi2Cut_ee.GetColumnNames();
	auto colNames_mumu = AfterChi2Cut_mumu.GetColumnNames();

	auto N_Columns_ee = colNames_ee.size();
	auto N_Columns_mumu = colNames_mumu.size();

	
	TFile * output_ee = new TFile(OutRootFile_ee.c_str(), "RECREATE");

	auto ZMassType = AfterChi2Cut_ee.GetColumnType("z_mass");
	auto EventWeightType = AfterChi2Cut_ee.GetColumnType("EventWeight");

	auto ZMassEntries = AfterChi2Cut_ee.Filter("z_mass").Count();
	auto EventWeightEntries = AfterChi2Cut_ee.Filter("EventWeight").Count();

	std::cout << "ZMassType = " << ZMassType << ". EventWeightType = " << EventWeightType << std::endl;
	std::cout << "ZMassEntries = " << *ZMassEntries << ". EventWeightEntries = " << *EventWeightEntries << std::endl;

	//auto h_z_mass_weighted = AfterChi2Cut_ee.Histo1D({"h_z_mass_weighted", "z_mass", 40, 0, 150}, "z_mass", "EventWeight");
	
	//auto h_EventWeight = AfterChi2Cut_ee.Histo1D({"h_EventWeight", "EventWeight", 40, 0, 1}, "EventWeight");
	auto h_ZMass = AfterChi2Cut_ee.Histo1D({"h_ZMass", "ZMass", 40, 0, 150}, "z_mass");
	auto h_BTagWeight = AfterChi2Cut_ee.Histo1D({"h_BTagWeight", "BTagWeight", 40, 0, 1}, "BTagWeight");

	//h_EventWeight->Write();
	h_ZMass->Write();
	h_BTagWeight->Write();
	output_ee->Close();

	std::cout << "after output_ee->Close()" << std::endl;

}
else{

	std::cout << "inside the else statment for when blinding is false" << std::endl;

	std::string OutRootFile_ee_unblinded;
	std::string OutRootFile_mumu_unblinded;


	if(NPL == false){
		OutRootFile_ee_unblinded = "Results_MC_" + process + "_" + year + "_ee.root";
		OutRootFile_mumu_unblinded = "Results_MC_" + process + "_" + year + "_mumu.root";
	}
	else{
		OutRootFile_ee_unblinded = "Results_MC_" + process + "_" + year + "_ee_NPL.root";
        	OutRootFile_mumu_unblinded = "Results_MC_" + process + "_" + year + "_mumu_NPL.root";
	}



	//ee output distributions
	TFile * output_ee = new TFile(OutRootFile_ee_unblinded.c_str(), "RECREATE");
	output_ee->cd();

	auto ZWeight_ee{[&NormalisationFactorFunction, &SF_ee, 			     &SF_Uncert_ee,
		         &LeptonEfficiencies_ScaleUp,  &LeptonEfficiencies_ScaleDown,
		         &PDF_ScaleUp,                 &PDF_ScaleDown,
		         &isr_up, 		       &isr_down,
		         &fsr_up,		       &fsr_down
		        ](const float& PU, const float& BTagWeight, const floats& ReturnedPSWeight, const float& CalculatedNominalWeight, const float& EGammaSF_egammaEff, const float& EGammaSF_egammaEffReco, const float& EGammaSF_egammaEffSys, const float& EGammaSF_egammaEffRecoSys, const float& CalculatedGeneratorWeight, const float& ME_SF){


  			if(LeptonEfficiencies_ScaleUp == true){return PU * NormalisationFactorFunction() * BTagWeight * (SF_ee += SF_Uncert_ee) * CalculatedNominalWeight * EGammaSF_egammaEffSys * EGammaSF_egammaEffRecoSys * CalculatedGeneratorWeight * ME_SF;}
  			else if(LeptonEfficiencies_ScaleDown == true){return PU * NormalisationFactorFunction() * (SF_ee -= SF_Uncert_ee) * CalculatedNominalWeight * EGammaSF_egammaEffSys * EGammaSF_egammaEffRecoSys * CalculatedGeneratorWeight * ME_SF;}
  			else if(PDF_ScaleUp == true){return PU * NormalisationFactorFunction() * BTagWeight * SF_ee * CalculatedNominalWeight * EGammaSF_egammaEff * EGammaSF_egammaEffReco * CalculatedGeneratorWeight * ME_SF;}
  			else if(PDF_ScaleDown == true){return PU * NormalisationFactorFunction() * BTagWeight * SF_ee * CalculatedNominalWeight * EGammaSF_egammaEff * EGammaSF_egammaEffReco * CalculatedGeneratorWeight * ME_SF;}
  			else if(isr_up == true){return PU * NormalisationFactorFunction() * BTagWeight * SF_ee * ReturnedPSWeight.at(2) * CalculatedNominalWeight * EGammaSF_egammaEff * EGammaSF_egammaEffReco * CalculatedGeneratorWeight * ME_SF;}
  			else if(isr_down == true){return PU * NormalisationFactorFunction() * BTagWeight * SF_ee * ReturnedPSWeight.at(0) * CalculatedNominalWeight * EGammaSF_egammaEff * EGammaSF_egammaEffReco * CalculatedGeneratorWeight * ME_SF;}
  			else if(fsr_up == true){return PU * NormalisationFactorFunction() * BTagWeight * SF_ee * ReturnedPSWeight.at(3) * CalculatedNominalWeight * EGammaSF_egammaEff * EGammaSF_egammaEffReco * CalculatedGeneratorWeight * ME_SF;}
  			else if(fsr_down == true){return PU * NormalisationFactorFunction() * BTagWeight * SF_ee * ReturnedPSWeight.at(1) * CalculatedNominalWeight * EGammaSF_egammaEff * EGammaSF_egammaEffReco * CalculatedGeneratorWeight * ME_SF;}
  			else{return PU * NormalisationFactorFunction() * BTagWeight * SF_ee * CalculatedNominalWeight * EGammaSF_egammaEff * EGammaSF_egammaEffReco * CalculatedGeneratorWeight * ME_SF;}
  
	}};

/*
"MuonSFTest_ID"
"MuonSFTest_Iso"
"MuonSFTest_ID_sys_syst"
"MuonSFTest_ID_sys_stat"
"MuonSFTest_Iso_sys_syst"
"MuonSFTest_Iso_sys_stat"
*/


        auto h_ZMass = d_WeightedEvents_withMET_ee.Define("ZWeight", ZWeight_ee, {"PU", "BTagWeight", "ReturnedPSWeight", "CalculatedNominalWeight", "EGammaSF_egammaEff", "EGammaSF_egammaEffReco", "EGammaSF_egammaEffSys", "EGammaSF_egammaEffRecoSys", "CalculatedGeneratorWeight", "ME_SF"}).Histo1D({"h_ZMass", "ZMass", 40, 0, 150}, "z_mass", "ZWeight");


	h_ZMass->Write();
	output_ee->Close();
	


}



std::cout << "before NPL estimation" << std::endl;

//NPL estimation
if(NPL == true){

std::cout << "inside NPL = true" << std::endl;

	//For the SF (Number of opposite-sign non prompt MC events / number of same-sign non-prompt MC events)
	//ee
	int Num_OppositeSignNonPrompt_MC_ee = *( d_WeightedEvents_ee.Filter("OppositeSignNonPrompt").Count() );
	int Num_SameSignNonPrompt_MC_ee = *( d_WeightedEvents_ee.Filter("SameSignNonPrompt").Count() );

	int NPL_SF_ee;

	if(Num_OppositeSignNonPrompt_MC_ee != 0 && Num_SameSignNonPrompt_MC_ee != 0){
		 NPL_SF_ee = Num_OppositeSignNonPrompt_MC_ee / Num_SameSignNonPrompt_MC_ee;
	}
	else{NPL_SF_ee = 99999;}

	//mumu
	int Num_OppositeSignNonPrompt_MC_mumu = *( d_WeightedEvents_mumu.Filter("OppositeSignNonPrompt").Count() );
        int Num_SameSignNonPrompt_MC_mumu = *( d_WeightedEvents_mumu.Filter("SameSignNonPrompt").Count() );

	int NPL_SF_mumu; 
 
        if(Num_OppositeSignNonPrompt_MC_mumu != 0 && Num_SameSignNonPrompt_MC_mumu != 0){ 
                 NPL_SF_mumu = Num_OppositeSignNonPrompt_MC_mumu / Num_SameSignNonPrompt_MC_mumu;
        }
        else{NPL_SF_mumu = 99999;}

	std::ofstream NPL_numbers_MC;
	std::string NPLNumbersStringMC = "NPL_SF_MC_" + process + "_" + year + ".txt";
	NPL_numbers_MC.open(NPLNumbersStringMC.c_str()); 

        NPL_numbers_MC << "Num_OppositeSignNonPrompt_MC_ee = " << Num_OppositeSignNonPrompt_MC_ee << '\n'
        	       << "Num_SameSignNonPrompt_MC_ee = " << Num_SameSignNonPrompt_MC_ee << '\n'
        	       << "Num_OppositeSignNonPrompt_MC_mumu = " << Num_OppositeSignNonPrompt_MC_mumu << '\n'
        	       << "Num_SameSignNonPrompt_MC_mumu = " << Num_SameSignNonPrompt_MC_mumu << '\n'
		       << "NPL_SF_ee = " << NPL_SF_ee << '\n'
		       << "NPL_SF_mumu = " << NPL_SF_mumu << std::endl;


	//For  (number of same sign data events)- (the expected number of real same sign events with charge misidentification)
	//For this, I have already saved a histogram of the number of same sign MC events to the ROOT output file (branch = "SameSign")
	//Will subtract this from the number of same sign events in data and then multiply by the scale factor saved in the text file

}




												

}
else if(

//data starts here
process == "data_DoubleEGRunB" ||
process == "data_DoubleEGRunC" ||
process == "data_DoubleEGRunD" ||
process == "data_DoubleEGRunE" ||
process == "data_DoubleEGRunF" ||
process == "data_DoubleEGRunG" ||
process == "data_DoubleEGRunH" ||
process == "data_EGRunB" ||
process == "data_EGRunC" ||
process == "data_EGRunD" ||
process == "data_SingleElectronRunB" ||
process == "data_SingleElectronRunC" ||
process == "data_SingleElectronRunD" ||
process == "data_SingleElectronRunE" ||
process == "data_SingleElectronRunF" ||
process == "data_SingleElectronRunG" ||
process == "data_SingleElectronRunH" ||
process == "data_DoubleMuonRunB" ||
process == "data_DoubleMuonRunC" ||
process == "data_DoubleMuonRunD" ||
process == "data_DoubleMuonRunE" ||
process == "data_DoubleMuonRunF" ||
process == "data_DoubleMuonRunG" ||
process == "data_DoubleMuonRunH" ||
process == "data_SingleMuonRunB" ||
process == "data_SingleMuonRunC" ||
process == "data_SingleMuonRunD" ||
process == "data_SingleMuonRunE" ||
process == "data_SingleMuonRunF" ||
process == "data_SingleMuonRunG" ||
process == "data_SingleMuonRunH" 
){


	if(process == "data_DoubleEGRunB" ||
   process == "data_DoubleEGRunC" ||
   process == "data_DoubleEGRunD" ||
   process == "data_DoubleEGRunE" ||
   process == "data_DoubleEGRunF" ||
   process == "data_DoubleEGRunG" ||
   process == "data_DoubleEGRunH" ||
   process == "data_DoubleMuonRunB" ||
   process == "data_DoubleMuonRunC" ||
   process == "data_DoubleMuonRunD" ||
   process == "data_DoubleMuonRunE" ||
   process == "data_DoubleMuonRunF" ||
   process == "data_DoubleMuonRunG" ||
   process == "data_DoubleMuonRunH"){




	//Filtering events that pass the double electron triggers 
	auto EventsPassingDoubleElectronTrigger = d_dataframe.Define("EventsPassingDoubleElectronTrigger", DoubleElectron, {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"})
							     .Filter(DoubleElectron, {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"});



	//Filtering events that pass the double muon triggers
	std::string Mass3p8Filter;
	std::string mumuFilter2;


 	if(year == "2018"){Mass3p8Filter = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", mumuFilter2 = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ";}
	else if(year == "2017"){mumuFilter2 = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", Mass3p8Filter = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ";}
	else if(year == "2016"){Mass3p8Filter = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ"; mumuFilter2 = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ";}
	else{std::cout << "Choose the year out of 2016, 2017 or 2018" << std::endl;}


	auto EventsPassingDoubleMuonTrigger = d_dataframe.Define("EventsPassingDoubleMuonTrigger", DoubleMuon, {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", mumuFilter2, Mass3p8Filter})
							 .Filter(DoubleMuon, {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", mumuFilter2, Mass3p8Filter});


	//Filtering events that pass the muon electron triggers
	std::string EmuString1;
        std::string EmuString2;
        std::string EmuString3;
	std::string EmuString4;
	std::string EmuString5;
	std::string EmuString6;

	if(year == "2016"){

		EmuString1 = "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ";
                EmuString2 = "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ";
                EmuString3 = "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ";
		EmuString4 = "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL"; 
		EmuString5 = "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL"; 
		EmuString6 = "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL";

	}
	else if(year == "2017" || year == "2018"){

		EmuString1 = "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ";
                EmuString2 = "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ";
                EmuString3 = "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ";
                EmuString4 = "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ";
                EmuString5 = "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ";
                EmuString6 = "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ";

        }
	else{std::cout << "Choose a year out of 2016, 2017 or 2018" << std::endl;}
 	

	
	auto EventsPassingElectronMuonTrigger = d_dataframe.Define("EventsPassingElectronMuonTrigger", MuonElectron, {EmuString1, EmuString2, EmuString3, EmuString4, EmuString5, EmuString6})
							   .Filter(MuonElectron, {EmuString1, EmuString2, EmuString3, EmuString4, EmuString5, EmuString6});




	//std::strings for output file names
	std::string Name_ee = "EventsPassingDoubleElectronTrigger_" + process + "_" + year + ".root";
	std::string Name_mumu = "EventsPassingDoubleMuonTrigger_" + process + "_" + year + ".root";
	std::string Name_emu = "EventsPassingElectronMuonTrigger_" + process + "_" + year + ".root";


	//Saving events that pass the double electron triggers and their event numbers
	auto Snapshot_EventsPassingDoubleElectronTrigger = EventsPassingDoubleElectronTrigger.Snapshot("Events", Name_ee.c_str(), "");

	//Saving events that pass the double muon triggers and their event numbers
        auto Snapshot_EventsPassingDoubleMuonTrigger = EventsPassingDoubleMuonTrigger.Snapshot("Events", Name_mumu.c_str(), "");	

	//Saving events that pass the double muon triggers and their event numbers
        auto Snapshot_EventsPassingElectronMuonTrigger = EventsPassingElectronMuonTrigger.Snapshot("Events", Name_emu.c_str(), ""); 

	 


}
else if(

process == "data_SingleElectronRunB" ||
process == "data_SingleElectronRunC" ||
process == "data_SingleElectronRunD" ||
process == "data_SingleElectronRunE" ||
process == "data_SingleElectronRunF" ||
process == "data_SingleElectronRunG" ||
process == "data_SingleElectronRunH" ||
process == "data_SingleMuonRunB" ||
process == "data_SingleMuonRunC" ||
process == "data_SingleMuonRunD" ||
process == "data_SingleMuonRunE" ||
process == "data_SingleMuonRunF" ||
process == "data_SingleMuonRunG" ||
process == "data_SingleMuonRunH"


){

	std::cout << "before Filtering events that pass the single electron triggers" << std::endl;

        //Filtering events that pass the single electron triggers 
	std::string SingleElectronTrigger1;
	std::string SingleElectronTrigger2;
	std::string SingleElectronTrigger3;
	std::string SingleElectronTrigger4;
        std::string SingleElectronTrigger5;

	if(year == "2016"){

		SingleElectronTrigger1 = "HLT_Ele32_eta2p1_WPTight_Gsf";
		SingleElectronTrigger2 = "HLT_Ele32_eta2p1_WPTight_Gsf";
		SingleElectronTrigger3 = "HLT_Ele25_eta2p1_WPTight_Gsf";
		SingleElectronTrigger4 = "HLT_Ele25_eta2p1_WPTight_Gsf";
                SingleElectronTrigger5 = "HLT_Ele27_WPTight_Gsf";

	}
	else if(year == "2017"){

		SingleElectronTrigger1 = "HLT_Ele35_WPTight_Gsf";
                SingleElectronTrigger2 = "HLT_Ele35_WPTight_Gsf";
                SingleElectronTrigger3 = "HLT_Ele35_WPTight_Gsf";
		SingleElectronTrigger4 = "HLT_Ele35_WPTight_Gsf";
		SingleElectronTrigger5 = "HLT_Ele35_WPTight_Gsf";
	}
	else if(year == "2018"){

		SingleElectronTrigger1 = "HLT_Ele32_WPTight_Gsf_L1DoubleEG";
                SingleElectronTrigger2 = "HLT_Ele32_WPTight_Gsf_L1DoubleEG";
		SingleElectronTrigger3 = "HLT_Ele32_WPTight_Gsf_L1DoubleEG";
		SingleElectronTrigger4 = "HLT_Ele32_WPTight_Gsf_L1DoubleEG";
                SingleElectronTrigger5 = "HLT_Ele32_WPTight_Gsf_L1DoubleEG";

	}
	else{std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl;}


/* 2nd 4th 5th for 2016
 * 3rd for 2017
 * 1st for 2018
const bool& HLT_Ele32_WPTight_Gsf_L1DoubleEG,
const bool& HLT_Ele32_eta2p1_WPTight_Gsf,
const bool& HLT_Ele35_WPTight_Gsf,
const bool& HLT_Ele25_eta2p1_WPTight_Gsf,
const bool& HLT_Ele27_WPTight_Gsf
*/

	std::cout << "before EventsPassingSingleElectronTrigger" << std::endl;

        auto EventsPassingSingleElectronTrigger = d_dataframe.Define("EventsPassingSingleElectronTrigger", SingleElectron, {SingleElectronTrigger1, SingleElectronTrigger2, SingleElectronTrigger3, SingleElectronTrigger4, SingleElectronTrigger5})
                                                             .Filter(SingleElectron, {SingleElectronTrigger1, SingleElectronTrigger2, SingleElectronTrigger3, SingleElectronTrigger4, SingleElectronTrigger5});

	std::cout << "before //Filtering events that pass the single muon triggers" << std::endl;


        //Filtering events that pass the single muon triggers
        auto EventsPassingSingleMuonTrigger = d_dataframe.Define("EventsPassingSingleMuonTrigger", SingleMuon, {"HLT_IsoMu24", "HLT_IsoMu27", "HLT_IsoMu24_eta2p1"})
                                                         .Filter(SingleMuon, {"HLT_IsoMu24", "HLT_IsoMu27", "HLT_IsoMu24_eta2p1"});



        //std::strings for output file names
        std::string Name_electron = "EventsPassingSingleElectronTrigger_" + process + "_" + year + ".root";
        std::string Name_muon = "EventsPassingSingleMuonTrigger_" + process + "_" + year + ".root";

	std::cout << "before Saving events that pass the single electron triggers and their event numbers" << std::endl;
	
        //Saving events that pass the single electron triggers and their event numbers
        auto Snapshot_EventsPassingSingleElectronTrigger = EventsPassingSingleElectronTrigger.Snapshot("Events", Name_electron.c_str(), "");

	std::cout << "before Saving events that pass the single muon triggers and their event numbers" << std::endl;

        //Saving events that pass the single muon triggers and their event numbers
        auto Snapshot_EventsPassingSingleMuonTrigger = EventsPassingSingleMuonTrigger.Snapshot("Events", Name_muon.c_str(), "");

	 


}


std::cout << "before d_DoubleCountCheck" << std::endl;

auto d_DoubleCountCheck = d_EventCleaning.Define("DummyColumn3", DummyColumnFunction, {"Electron_pt"});


if(process != "data_EGRunB" &&
   process != "data_EGRunC" &&
   process != "data_EGRunD"){


	std::cout << "before ee_file" << std::endl;
	
	//Preventing double counting of single and double lepton datasets 
	std::string ee_file = "EventsPassingDoubleElectronTrigger_" + process + "_" + year + ".root";
	std::string mumu_file = "EventsPassingDoubleMuonTrigger_" + process + "_" + year + ".root";
	std::string emu_file = "EventsPassingElectronMuonTrigger_" + process + "_" + year + ".root";
	std::string e_file = "EventsPassingSingleElectronTrigger_" + process + "_" + year + ".root";
	std::string mu_file = "EventsPassingSingleMuonTrigger_" + process + "_" + year + ".root";

	std::cout << "after mu_file" << std::endl;

	TFile * ee_file_Open = new TFile(ee_file.c_str(), "READ");
	TFile * mumu_file_Open = new TFile(mumu_file.c_str(), "READ");
	TFile * emu_file_Open = new TFile(emu_file.c_str(), "READ");
	TFile * e_file_Open = new TFile(e_file.c_str(), "READ");
	TFile * mu_file_Open = new TFile(mu_file.c_str(), "READ");

	std::cout << "before h_EventNumber_ee" << std::endl;

	TH1* h_EventNumber_ee = (TH1*)ee_file_Open->GetObjectChecked("event", "TH1");
	
	std::cout << "before h_EventNumber_mumu" << std::endl;

	TH1* h_EventNumber_mumu = (TH1*)mumu_file_Open->GetObjectChecked("event", "TH1");
	TH1* h_EventNumber_emu = (TH1*)emu_file_Open->GetObjectChecked("event", "TH1");
	TH1* h_EventNumber_e = (TH1*)ee_file_Open->GetObjectChecked("event", "TH1");
	TH1* h_EventNumber_mu = (TH1*)mu_file_Open->GetObjectChecked("event", "TH1");

	std::cout << "after h_EventNumber_mu" << std::endl;

	//Returning events that have not been double-counted
	//For ee and single electron 
	auto DoubleCountCheck_ee_channel{[&h_EventNumber_ee, &h_EventNumber_e](){

	  std::cout << "inside oubleCountCheck_ee_channel" << std::endl;

	  floats EventNumbers{};

	  for(int i = 0; i < h_EventNumber_ee->GetNbinsX(); i++){

		std::cout << "h_EventNumber_ee->GetNbinsX() = " << h_EventNumber_ee->GetNbinsX() << std::endl;

		float DoubleElectronEventNumber = h_EventNumber_ee->GetBinContent(i);

		for(int j = 0; j < h_EventNumber_e->GetNbinsX(); j++){

			if( DoubleElectronEventNumber  != h_EventNumber_e->GetBinContent(j) ){ EventNumbers.push_back(DoubleElectronEventNumber); }
			else{continue;}
		
		}

	  }


	 return EventNumbers.at(0);


	}};

	//For mumu and single muon
	auto DoubleCountCheck_mumu_channel{[&h_EventNumber_mumu, &h_EventNumber_mu](){

	  std::cout << "inside DoubleCountCheck_mumu_channel" << std::endl;

          floats EventNumbers{};

          for(int i = 0; i < h_EventNumber_mumu->GetNbinsX(); i++){
        
                float DoubleMuonEventNumber = h_EventNumber_mumu->GetBinContent(i);

                for(int j = 0; j < h_EventNumber_mu->GetNbinsX(); j++){

                        if( DoubleMuonEventNumber  != h_EventNumber_mu->GetBinContent(j) ){ EventNumbers.push_back(DoubleMuonEventNumber); }
                        else{continue;}

                }

          }


         return EventNumbers.at(0);


        }};


	//For emu and single electron
	auto DoubleCountCheck1_emu_channel{[&h_EventNumber_emu, &h_EventNumber_e](){

	  std::cout << "inside DoubleCountCheck1_emu_channel" << std::endl;

          floats EventNumbers{};

          for(int i = 0; i < h_EventNumber_emu->GetNbinsX(); i++){

                float MuonElectronEventNumber = h_EventNumber_emu->GetBinContent(i);

                for(int j = 0; j < h_EventNumber_e->GetNbinsX(); j++){

                        if( MuonElectronEventNumber  != h_EventNumber_e->GetBinContent(j) ){ EventNumbers.push_back(MuonElectronEventNumber); }
                        else{continue;}

                }

          }


        return EventNumbers.at(0);


        }};


	//For emu and single muon
        auto DoubleCountCheck2_emu_channel{[&h_EventNumber_emu, &h_EventNumber_mu](){

	 std::cout << "inside DoubleCountCheck2_emu_channel" << std::endl;

          floats EventNumbers{};

          for(int i = 0; i < h_EventNumber_emu->GetNbinsX(); i++){

                float MuonElectronEventNumber = h_EventNumber_emu->GetBinContent(i);

                for(int j = 0; j < h_EventNumber_mu->GetNbinsX(); j++){

                        if( MuonElectronEventNumber  != h_EventNumber_mu->GetBinContent(j) ){ EventNumbers.push_back(MuonElectronEventNumber); }
                        else{continue;}

                }

          }


         return EventNumbers.at(0);


        }};	


	std::cout << "before PreventDoubleCount_FilterFunction_Not2018" << std::endl;

	auto PreventDoubleCount_FilterFunction_Not2018{[&DoubleCountCheck_ee_channel, &DoubleCountCheck_mumu_channel, &DoubleCountCheck1_emu_channel, DoubleCountCheck2_emu_channel](const unsigned int& event){

		std::cout << "inside PreventDoubleCount_FilterFunction_Not2018" << std::endl;

		return (event == DoubleCountCheck_ee_channel()) && (event == DoubleCountCheck_mumu_channel()) && (event == DoubleCountCheck1_emu_channel()) && (event == DoubleCountCheck2_emu_channel());

	}};

	std::cout << "before PreventDoubleCount_FilterFunction_2018" << std::endl;

	//The double and single electron datasets for 2018 are merged. Only need to check for double counting for muons. 
	auto PreventDoubleCount_FilterFunction_2018{[&DoubleCountCheck_mumu_channel, &DoubleCountCheck1_emu_channel, DoubleCountCheck2_emu_channel](const unsigned int& event){

		std::cout << "inside PreventDoubleCount_FilterFunction_2018" << std::endl;

                return (event == DoubleCountCheck_mumu_channel()) && (event == DoubleCountCheck1_emu_channel()) && (event == DoubleCountCheck2_emu_channel());

        }};

  std::cout << "before d_DoubleCountCheck" << std::endl;

  auto d_DoubleCountCheck = d_EventCleaning.Define("DummyColumn1", DummyColumnFunction, {"Electron_pt"});

  std::cout << "before year != 2018" << std::endl;

  if(year != "2018"){ auto d_DoubleCountCheck = d_EventCleaning.Filter(PreventDoubleCount_FilterFunction_Not2018, {"event"}, "Double counting filter");}
  else{ auto d_DoubleCountCheck = d_EventCleaning.Filter(PreventDoubleCount_FilterFunction_2018, {"event"}, "Double counting filter");}


  std::cout << "after else{ auto d_DoubleCountCheck" << std::endl;

}

std::cout << "before GoldenJsonReader" << std::endl;

//Using the golden ison file to filter events
auto GoldenJsonReader{[&year](){

 std::cout << "GoldenJsonReader" << std::endl;

 std::string GoldenJsonFileName;

 if(year == "2016"){GoldenJsonFileName = "./ScaleFactors/GoldenJSON/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt";}
 else if(year == "2017"){GoldenJsonFileName = "./ScaleFactors/GoldenJSON/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt";}
 else if(year == "2018"){GoldenJsonFileName = "./ScaleFactors/GoldenJSON/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt";}
 else{std::cout << "Choose the year out of 2016, 2017 or 2018" << std::endl;}


 std::ifstream myReadFile;
 myReadFile.open(GoldenJsonFileName);
 static char output[100];
 std::vector<std::string> GoldenJsonOutput{};
  

 if (myReadFile.is_open()) {
 while (!myReadFile.eof()) {


    myReadFile >> output;
    GoldenJsonOutput.push_back(output);

 }
}

myReadFile.close();
return GoldenJsonOutput;

}};


std::cout << "before GoldenJson_SplitChars" << std::endl;


auto GoldenJson_SplitChars{[&year, &GoldenJsonReader](){

  std::cout << "inside GoldenJson_SplitChars" << std::endl;

  std::vector<char> out{};

  for(int i = 0; i < (GoldenJsonReader()).size(); i++){

  	std::string element = GoldenJsonReader().at(i);

  	for(int j = 0; j < element.size(); j++){

  		out.push_back(element.at(j));  	

  	}

  }

  return out;

}};


std::cout << "before RunNumberCheck" << std::endl;

auto RunNumberCheck{[&year, &GoldenJson_SplitChars](const unsigned int& InputRunNumber){

  std::cout << "inside RunNumberCheck" << std::endl;
  
 std::vector<char> EventsVector{}; 

 for(int i = 0; i < (GoldenJson_SplitChars()).size(); i++){

	int RunNumBeingRead;

 	if(  GoldenJson_SplitChars().at(i+1) == '"' &&
	    (GoldenJson_SplitChars().at(i+2) == '2' ||
	     GoldenJson_SplitChars().at(i+2) == '3')  ){ 

		int digit1 = GoldenJson_SplitChars().at(i+2) - '0'; 
		int digit2 = GoldenJson_SplitChars().at(i+3) - '0'; 
		int digit3 = GoldenJson_SplitChars().at(i+4) - '0'; 
		int digit4 = GoldenJson_SplitChars().at(i+5) - '0'; 
		int digit5 = GoldenJson_SplitChars().at(i+6) - '0'; 
		int digit6 = GoldenJson_SplitChars().at(i+7) - '0';


		int run = (digit1*100000) + (digit2*10000) + (digit3*1000) + (digit4*100) + (digit5*10) + digit6;
	
		RunNumBeingRead = run;

		if(run == InputRunNumber){


			for(int j = 2; j < (GoldenJson_SplitChars()).size(); j++){

				if(GoldenJson_SplitChars().at(i+10) == '[' && 
				   GoldenJson_SplitChars().at(i+11) == '['){	

					if( GoldenJson_SplitChars().at( (i+10)+j ) == ']' &&
					    GoldenJson_SplitChars().at( (i+10)+(j+1) ) == ']'){

					
						for(int k = (i+10); k < ((i+10)+(j+2)); k++){

							EventsVector.push_back(GoldenJson_SplitChars().at(k));
						
						}

						return EventsVector;	
											

					}
					else{std::cout << "GoldenJson_SplitChars().at((i+11)+j) and at (i+11)+(j+1) are not ']' " << std::endl; continue;}


				}
				else{std::cout << "GoldenJson_SplitChars().at(i+10) is not '[' " << std::endl; continue;}


			}	

		}
		else{continue;}

	}
	else{std::cout << "The run number of " << RunNumBeingRead << " does not match the input run number of " << InputRunNumber << std::endl;}


 }


}};


std::cout << "before ReturnRunNumAndEventRanges" << std::endl;

auto ReturnRunNumAndEventRanges{[&year, &RunNumberCheck](const unsigned int& InputRunNumber){

 std::cout << "inside ReturnRunNumAndEventRanges" << std::endl;
 
 std::vector<int> RunNumAndEvents{};

 RunNumAndEvents.push_back(InputRunNumber);

 std::vector<char> Runs = RunNumberCheck(InputRunNumber);


 for(int i = 0; i < Runs.size(); i++){


 	if(Runs.at(i) == ']' && Runs.at(i+1) == ']'){break;}
	else if( isdigit(Runs.at(i)) || Runs.at(i) == ',' || Runs.at(i) == ' ' || Runs.at(i) == ']'){continue;}
 	else if( (Runs.at(i) == '[' && Runs.at(i+1) == '[') ||
		 (Runs.at(i) == '[' && isdigit(Runs.at(i+1))) ){


		if(  isdigit( Runs.at(i+1) ) && //For the min value being a 4 digit number
      	     	     isdigit( Runs.at(i+2) ) &&
      	     	     isdigit( Runs.at(i+3) ) && 
		     isdigit( Runs.at(i+4) ) &&
      	     	     Runs.at(i+5) == ','){


			int Min_Digit1 = Runs.at(i+1) - '0';
			int Min_Digit2 = Runs.at(i+2) - '0';
			int Min_Digit3 = Runs.at(i+3) - '0';
			int Min_Digit4 = Runs.at(i+4) - '0';

			int EventNumber_Min = (Min_Digit1*1000) + (Min_Digit2*100) + (Min_Digit3*10) + Min_Digit4;
			RunNumAndEvents.push_back(EventNumber_Min);

			if(isdigit( Runs.at(i+6) ) &&
                           isdigit( Runs.at(i+7) ) &&
                           isdigit( Runs.at(i+8) ) &&
                           isdigit( Runs.at(i+9) ) &&
                           Runs.at(i+10) == ']' ){ //For the min value being a 4 digit number and the max value being a 4 digit number
                                   

                                        int Max_Digit1 = Runs.at(i+6) - '0';
                                        int Max_Digit2 = Runs.at(i+7) - '0';
                                        int Max_Digit3 = Runs.at(i+8) - '0';
                                        int Max_Digit4 = Runs.at(i+9) - '0';

                                        int EventNumber_Max = (Max_Digit1*1000) + (Max_Digit2*100) + (Max_Digit3*10) + Max_Digit4;
                                        RunNumAndEvents.push_back(EventNumber_Max);


                                }
			else if(isdigit( Runs.at(i+6) ) && //For the min value being a 4 digit number and the max value being a 3 digit number
      	   	   	        isdigit( Runs.at(i+7) ) &&
      	   	   	        isdigit( Runs.at(i+8) ) &&
      	   	   	        Runs.at(i+9) == ']'){


				int Max_Digit1 = Runs.at(i+6) - '0';
        			int Max_Digit2 = Runs.at(i+7) - '0';
        			int Max_Digit3 = Runs.at(i+8) - '0';
	
				int EventNumber_Max = (Max_Digit1*100) + (Max_Digit2*10) + Max_Digit3;
        			RunNumAndEvents.push_back(EventNumber_Max);
		

			}
			else if(isdigit( Runs.at(i+6) ) && //For the min value being a 4 digit number and the max value being a 2 digit number 
         			isdigit( Runs.at(i+7) ) &&
         			Runs.at(i+8) == ']' ){

					int Max_Digit1 = Runs.at(i+6) - '0';
                			int Max_Digit2 = Runs.at(i+7) - '0';

                			int EventNumber_Max = (Max_Digit1*10) + Max_Digit2;
                			RunNumAndEvents.push_back(EventNumber_Max);


			}
			else if(isdigit( Runs.at(i+6) ) &&
                                Runs.at(i+7) == ']' ){ //For the min value being a 4 digit number and the max value being a 1 digit number


                                        int Max_Digit1 = Runs.at(i+6) - '0';

                                        int EventNumber_Max = Max_Digit1;
                                        RunNumAndEvents.push_back(EventNumber_Max);


                        }
                       else{std::cout << "error" << std::endl;}

 		}	
 		else if(  isdigit( Runs.at(i+1) ) && //For the min value being a 3 digit number
      	     	     isdigit( Runs.at(i+2) ) &&
      	     	     isdigit( Runs.at(i+3) ) && 
      	     	     Runs.at(i+4) == ','){


			int Min_Digit1 = Runs.at(i+1) - '0';
			int Min_Digit2 = Runs.at(i+2) - '0';
			int Min_Digit3 = Runs.at(i+3) - '0';

			int EventNumber_Min = (Min_Digit1*100) + (Min_Digit2*10) + Min_Digit3;
			RunNumAndEvents.push_back(EventNumber_Min);

			if(isdigit( Runs.at(i+5) ) &&
                           isdigit( Runs.at(i+6) ) &&
                           isdigit( Runs.at(i+7) ) &&
                           isdigit( Runs.at(i+8) ) &&
                           Runs.at(i+9) == ']' ){ //For the min value being a 3 digit number and the max value being a 4 digit number
                                   

                                        int Max_Digit1 = Runs.at(i+5) - '0';
                                        int Max_Digit2 = Runs.at(i+6) - '0';
                                        int Max_Digit3 = Runs.at(i+7) - '0';
                                        int Max_Digit4 = Runs.at(i+8) - '0';


                                        int EventNumber_Max = (Max_Digit1*1000) + (Max_Digit2*100) + (Max_Digit3*10) + Max_Digit4;

                                        RunNumAndEvents.push_back(EventNumber_Max);


                                }
			else if(isdigit( Runs.at(i+5) ) && //For the min value being a 3 digit number and the max value being a 3 digit number
      	   	   	        isdigit( Runs.at(i+6) ) &&
      	   	   	        isdigit( Runs.at(i+7) ) &&
      	   	   	        Runs.at(i+8) == ']'){


				int Max_Digit1 = Runs.at(i+5) - '0';
        			int Max_Digit2 = Runs.at(i+6) - '0';
        			int Max_Digit3 = Runs.at(i+7) - '0';
	
				int EventNumber_Max = (Max_Digit1*100) + (Max_Digit2*10) + Max_Digit3;
        			RunNumAndEvents.push_back(EventNumber_Max);
		

			}
			else if(isdigit( Runs.at(i+5) ) && //For the min value being a 3 digit number and the max value being a 2 digit number 
         			isdigit( Runs.at(i+6) ) &&
         			Runs.at(i+7) == ']' ){

					int Max_Digit1 = Runs.at(i+5) - '0';
                			int Max_Digit2 = Runs.at(i+6) - '0';

                			int EventNumber_Max = (Max_Digit1*10) + Max_Digit2;
                			RunNumAndEvents.push_back(EventNumber_Max);


			}
			else if(isdigit( Runs.at(i+5) ) &&
                                Runs.at(i+6) == ']' ){ //For the min value being a 3 digit number and the max value being a 1 digit number


                                        int Max_Digit1 = Runs.at(i+5) - '0';

                                        int EventNumber_Max = Max_Digit1;
                                        RunNumAndEvents.push_back(EventNumber_Max);


                        }
                       else{std::cout << "error" << std::endl;}

 		}
 		else if(isdigit( Runs.at(i+1) ) && //For the min value being a 2 digit number
         		isdigit( Runs.at(i+2) ) &&
         		Runs.at(i+3) == ',' ){


				int Min_Digit1 = Runs.at(i+1) - '0';
        			int Min_Digit2 = Runs.at(i+2) - '0';

        			int EventNumber_Min = (Min_Digit1*10) + Min_Digit2;
        			RunNumAndEvents.push_back(EventNumber_Min);

				if(isdigit( Runs.at(i+4) ) &&
                                   isdigit( Runs.at(i+5) ) &&
				   isdigit( Runs.at(i+6) ) &&
                                   isdigit( Runs.at(i+7) ) &&
                                   Runs.at(i+8) == ']' ){ //For the min value being a 2 digit number and the max value being a 4 digit number


                                        int Max_Digit1 = Runs.at(i+4) - '0';
                                        int Max_Digit2 = Runs.at(i+5) - '0';
					int Max_Digit3 = Runs.at(i+6) - '0';
                                        int Max_Digit4 = Runs.at(i+7) - '0';

                                        int EventNumber_Max = (Max_Digit1*1000) + (Max_Digit2*100) + (Max_Digit3*10) * Max_Digit4;
                                        RunNumAndEvents.push_back(EventNumber_Max);


                                }
				else if(isdigit( Runs.at(i+4) ) && //For the min value being a 2 digit number and the max value being a 3 digit nummber
           	   	   	   isdigit( Runs.at(i+5) ) &&
           	   	   	   isdigit( Runs.at(i+6) ) &&
           	   	   	   Runs.at(i+7) == ']'){


                			int Max_Digit1 = Runs.at(i+4) - '0';
                			int Max_Digit2 = Runs.at(i+5) - '0';
                			int Max_Digit3 = Runs.at(i+6) - '0';

                			int EventNumber_Max = (Max_Digit1*100) + (Max_Digit2*10) + Max_Digit3;
                			RunNumAndEvents.push_back(EventNumber_Max);


        			}
        			else if(isdigit( Runs.at(i+4) ) &&
                			isdigit( Runs.at(i+5) ) &&
                			Runs.at(i+6) == ']' ){ //For the min value being a 2 digit number and the max value being a 2 digit number


                        		int Max_Digit1 = Runs.at(i+4) - '0';
                        		int Max_Digit2 = Runs.at(i+5) - '0';

                        		int EventNumber_Max = (Max_Digit1*10) + Max_Digit2;
                        		RunNumAndEvents.push_back(EventNumber_Max);


        			}
				else if(isdigit( Runs.at(i+4) ) &&
                                        Runs.at(i+5) == ']' ){ //For the min value being a 2 digit number and the max value being a 1 digit number


                                        int Max_Digit1 = Runs.at(i+4) - '0';

                                        int EventNumber_Max = Max_Digit1;
                                        RunNumAndEvents.push_back(EventNumber_Max);


                                }
				else{std::cout << "error" <<  "Runs.at(i+1) = " << Runs.at(i+1) << '\n' << "Runs.at(i+2) = " << Runs.at(i+2) << std::endl;}

   	}
	else if(  isdigit( Runs.at(i+1) ) && //For the min value being a 1 digit number
      	     	  Runs.at(i+2) == ',' &&
		  isdigit(Runs.at(i+3)) ){


			int Min_Digit1 = Runs.at(i+1) - '0';	

			int EventNumber_Min = Min_Digit1;

			RunNumAndEvents.push_back(EventNumber_Min);

			if(isdigit( Runs.at(i+3) ) &&
                           isdigit( Runs.at(i+4) ) &&
                           isdigit( Runs.at(i+5) ) &&
                           isdigit( Runs.at(i+6) ) &&
                           Runs.at(i+7) == ']' ){ //For the min value being a 1 digit number and the max value being a 4 digit number
                                   

                                        int Max_Digit1 = Runs.at(i+3) - '0';
                                        int Max_Digit2 = Runs.at(i+4) - '0';
                                        int Max_Digit3 = Runs.at(i+5) - '0';
                                        int Max_Digit4 = Runs.at(i+6) - '0';

                                        int EventNumber_Max = (Max_Digit1*1000) + (Max_Digit2*100) + (Max_Digit3*10) * Max_Digit4;
                                        RunNumAndEvents.push_back(EventNumber_Max);


                                }
			else if(isdigit( Runs.at(i+3) ) && //For the min value being a 1 digit number and the max value being a 3 digit number
      	   	   	        isdigit( Runs.at(i+4) ) &&
      	   	   	        isdigit( Runs.at(i+5) ) &&
      	   	   	        Runs.at(i+6) == ']'){


				int Max_Digit1 = Runs.at(i+3) - '0';
        			int Max_Digit2 = Runs.at(i+4) - '0';
        			int Max_Digit3 = Runs.at(i+5) - '0';
	
				int EventNumber_Max = (Max_Digit1*100) + (Max_Digit2*10) + Max_Digit3;
        			RunNumAndEvents.push_back(EventNumber_Max);
		

			}
			else if(isdigit( Runs.at(i+3) ) && //For the min value being a 1 digit number and the max value being a 2 digit number 
         			isdigit( Runs.at(i+4) ) &&
         			Runs.at(i+5) == ']' ){

					int Max_Digit1 = Runs.at(i+3) - '0';
                			int Max_Digit2 = Runs.at(i+4) - '0';

                			int EventNumber_Max = (Max_Digit1*10) + Max_Digit2;
                			RunNumAndEvents.push_back(EventNumber_Max);


			}
			else if(isdigit( Runs.at(i+3) ) &&
                                Runs.at(i+4) == ']' ){ //For the min value being a 1 digit number and the max value being a 1 digit number


                                        int Max_Digit1 = Runs.at(i+3) - '0';

                                        int EventNumber_Max = Max_Digit1;
                                        RunNumAndEvents.push_back(EventNumber_Max);


                        }
                       else{std::cout << "error" << std::endl;}

 		}	
		else{std::cout << "INSIDE THE ELSE STATEMENT" << '\n' << "Runs.at(i) = " << Runs.at(i) << '\n' << "Runs.at(i+1) = " << Runs.at(i+1) << '\n' << "Runs.at(i+2) = " << Runs.at(i+2) << std::endl;}


	}	 


  }


  return RunNumAndEvents;


}};


std::cout << "before RunAndLumiFilterFunction" << std::endl;

auto RunAndLumiFilterFunction{[&ReturnRunNumAndEventRanges](const unsigned int& InputRunNumber, const unsigned int& luminosityBlock){

  std::cout << "inside RunAndLumiFilterFunction" << std::endl;

  if( InputRunNumber == ReturnRunNumAndEventRanges(InputRunNumber).at(0) ){

	for(int i = 1; i < ReturnRunNumAndEventRanges(InputRunNumber).size(); i+=2){

		int MinLumi = ReturnRunNumAndEventRanges(InputRunNumber).at(i);
		int MaxLumi = ReturnRunNumAndEventRanges(InputRunNumber).at(i+1);

		if(luminosityBlock > MinLumi && luminosityBlock < MaxLumi){return InputRunNumber && luminosityBlock;}
		else{continue;}

	}

  }
  else{return false;}


}};

std::cout << "before d_GoldenJsonFilteredEvents" << std::endl;


//Filtering events with the run numbers and lumi ranges from the golden json file
auto d_GoldenJsonFilteredEvents = d_EventCleaning.Define("DummyColumn2", DummyColumnFunction, {"Electron_pt"});


if(process != "data_EGRunB" &&
   process != "data_EGRunC" &&
   process != "data_EGRunD"){

	auto d_GoldenJsonFilteredEvents = d_DoubleCountCheck.Filter(RunAndLumiFilterFunction, {"run", "luminosityBlock"}, "GoldenJson filter");

}
else{auto d_GoldenJsonFilteredEvents = d_EventCleaning.Filter(RunAndLumiFilterFunction, {"run", "luminosityBlock"}, "GoldenJson filter");}




 //Filtering events that pass the ee selection criteria
auto d_ee_selection_defines = d_GoldenJsonFilteredEvents.Define("PU", PU_function, {"PV_npvs"})
                                           .Define("TightElectrons", TightElectronsFunction, {"Electron_pt", "Electron_eta", "Electron_cutBased", "Electron_isPFcand"})
                                          .Define("Electron_pt_Selection", select<floats>, {"Electron_pt", "TightElectrons"})
                                          .Define("Electron_charge_Selection", select<ints>, {"Electron_charge", "TightElectrons"})
                                          .Define("LooseElectrons", LooseElectronsFunction, {"Electron_pt", "Electron_eta", "Electron_cutBased", "Electron_isPFcand"})
                                          .Define("LooseElectron_pt_Selection", select<floats>, {"Electron_pt", "LooseElectrons"})
                                          .Define("LooseElectron_charge_Selection", select<ints>, {"Electron_charge", "LooseElectrons"})
                                          .Define("OppositeSign", OppositeSign, {"Electron_charge_Selection"})
                                          .Define("SameSign", SameSign, {"Electron_charge_Selection"})
                                          .Define("Electron_eta_Selection", {"Electron_eta"})
                                          .Define("LeadingElectron_pT", LeadingVariable, {"Electron_pt_Selection"})
                                          .Define("SubleadingElectron_pT", SubleadingVariable, {"Electron_pt_Selection"})
                                          .Define("LeadingElectronPhi", LeadingVariable, {"Electron_phi"})
                                          .Define("SubleadingElectronPhi", SubleadingVariable, {"Electron_phi"})
                                          .Define("LeadingElectronMass", LeadingVariable, {"Electron_mass"})
                                          .Define("SubleadingElectronMass", SubleadingVariable, {"Electron_mass"})
                                          .Define("LeadingElectron_RelIso_Selection", LeadingVariable, {"Electron_jetRelIso"})
                                          .Define("SubleadingElectron_RelIso_Selection", SubleadingVariable, {"Electron_jetRelIso"})
                                          .Define("LeadingElectronEta", LeadingVariable, {"Electron_eta_Selection"})
                                          .Define("SubleadingElectronEta", SubleadingVariable, {"Electron_eta_Selection"})
                                          .Define("LeadingElectron_dz_ECALBarrel", LeadingElectron_dz_ECALBarrel_function, {"LeadingElectron_pT", "Electron_eta_Selection", "Electron_dz"})
                                          .Define("LeadingElectron_dxy_ECALBarrel", LeadingElectron_dxy_ECALBarrel_function, {"LeadingElectron_pT", "Electron_eta_Selection", "Electron_dxy"})
                                          .Define("LeadingElectron_dz_ECALEndcaps", LeadingElectron_dz_ECALEndcaps_function, {"LeadingElectron_pT", "Electron_eta_Selection", "Electron_dz"})
                                          .Define("LeadingElectron_dxy_ECALEndcaps", LeadingElectron_dxy_ECALEndcaps_function, {"LeadingElectron_pT", "Electron_eta_Selection", "Electron_dxy"})
                                          .Define("SubleadingElectron_dz_ECALBarrel", SubleadingElectron_dz_ECALBarrel_function, {"SubleadingElectron_pT", "Electron_eta_Selection", "Electron_dz"})
                                          .Define("SubleadingElectron_dxy_ECALBarrel", SubleadingElectron_dxy_ECALBarrel_function, {"SubleadingElectron_pT", "Electron_eta_Selection", "Electron_dxy"})
                                          .Define("SubleadingElectron_dz_ECALEndcaps", SubleadingElectron_dz_ECALEndcaps_function, {"SubleadingElectron_pT", "Electron_eta_Selection", "Electron_dz"})
                                          .Define("SubleadingElectron_dxy_ECALEndcaps", SubleadingElectron_dxy_ECALEndcaps_function, {"SubleadingElectron_pT", "Electron_eta_Selection", "Electron_dxy"});


 std::cout << "before  d_ee_selection" << std::endl;

 auto d_ee_selection = d_ee_selection_defines.Filter(lep_cut_ee, lep_cut_ee_strings, "lepton cut (ee)");

 std::cout << "after d_ee_selection" << std::endl;

if(ttbarCR == false){auto d_ee_selection = d_ee_selection_defines.Filter(lep_cut_ee, lep_cut_ee_strings, "lepton cut (ee)");}
else{auto d_ee_selection = d_ee_selection_defines.Filter(lep_cut_emu, lep_cut_emu_strings, "lepton cut (emu)");}


std::cout << "before d_mumu_selection_defines" << std::endl;


//Filtering events that pass the mumu selection criteria
auto d_mumu_selection_defines = d_GoldenJsonFilteredEvents.Define("PU", PU_function, {"PV_npvs"})
                                             .Define("TightMuons", TightMuonsFunction, {"Muon_isPFcand", "Muon_pt", "Muon_eta", "Muon_tightId", "Muon_pfRelIso04_all"})
                                             .Define("Muon_pt_Selection", select<floats>, {"Muon_pt", "TightMuons"})
                                             .Define("Muon_charge_Selection", select<ints>, {"Muon_charge", "TightMuons"})
                                             .Define("LooseMuons", LooseMuonsFunction, {"Muon_isPFcand", "Muon_pt", "Muon_eta", "Muon_softId", "Muon_pfRelIso04_all"})
                                             .Define("LooseMuon_pt_Selection", select<floats>, {"Muon_pt", "LooseMuons"})
                                             .Define("LooseMuon_charge_Selection", select<ints>, {"Muon_charge", "LooseMuons"})
                                             .Define("OppositeSign", OppositeSign, {"Muon_charge_Selection"})
                                             .Define("SameSign", SameSign, {"Muon_charge_Selection"})
                                             .Define("Muon_eta_Selection", {"Muon_eta"})
                                             .Define("LeadingMuon_pT", LeadingVariable, {"Muon_pt_Selection"})
                                             .Define("SubleadingMuon_pT", SubleadingVariable, {"Muon_pt_Selection"})
                                             .Define("LeadingMuonPhi", LeadingVariable, {"Muon_phi"})
                                             .Define("SubleadingMuonPhi", SubleadingVariable, {"Muon_phi"})
                                             .Define("LeadingMuonMass", LeadingVariable, {"Muon_mass"})
                                             .Define("SubleadingMuonMass", SubleadingVariable, {"Muon_mass"})
                                             .Define("LeadingMuon_RelIso_Selection", LeadingVariable, {"Muon_jetRelIso"})
                                             .Define("SubleadingMuon_RelIso_Selection", SubleadingVariable, {"Muon_jetRelIso"})
                                             .Define("LeadingMuonEta", LeadingVariable, {"Muon_eta_Selection"})
                                             .Define("SubleadingMuonEta", SubleadingVariable, {"Muon_eta_Selection"})
					     .Define("DummyColumnInts", DummyColumnFunctionInts, {"Muon_charge"})
                                             .Define("MuonFourMomentum", MuonFourMomentum, {"Muon_pt_Selection", "Muon_eta", "Muon_phi", "Muon_mass"})
                                             .Define("RochCorrVec", RochCorrVec_Function_data, {"Muon_charge_Selection", "Muon_pt_Selection", "Muon_eta", "Muon_phi", "DummyColumnInts", "Muon_nTrackerLayers"})
                                             .Define("MuonFourMomentum_RochCorr", RochCorrMuon4Mo, {"MuonFourMomentum", "RochCorrVec"})
					     .Define("MuonPt_RochCorr", TLorentzVectorPt, {"MuonFourMomentum_RochCorr"})
                                             .Define("MuonEta_RochCorr", TLorentzVectorEta, {"MuonFourMomentum_RochCorr"})
                                             .Define("MuonPhi_RochCorr", TLorentzVectorPhi, {"MuonFourMomentum_RochCorr"})
                                             .Define("MuonMass_RochCorr", TLorentzVectorMass, {"MuonFourMomentum_RochCorr"});

 std::cout << "before d_mumu_selection" << std::endl;

 auto d_mumu_selection = d_mumu_selection_defines.Filter(lep_cut_mumu, lep_cut_mumu_strings, "lepton cut (mumu)");

 std::cout << "before if(ttbarCR == false){" << std::endl;

 if(ttbarCR == false){auto d_mumu_selection = d_mumu_selection_defines.Filter(lep_cut_mumu, lep_cut_mumu_strings, "lepton cut (mumu)");}
 else{auto d_mumu_selection = d_mumu_selection_defines.Filter(lep_cut_emu, lep_cut_emu_strings, "lepton cut (emu)");}
 

 std::cout << "before d_ee_recoZ_selection" << std::endl;


  //Filtering events with a reconstructed Z boson
 auto d_ee_recoZ_selection = d_ee_selection.Define("z_lep_eta", "Electron_eta[TightElectrons]")
                                          .Define("z_lep_phi", "Electron_phi[TightElectrons]")
                                          .Define("z_lep_mass", "Electron_mass[TightElectrons]")
                                          .Define("z_lep_pt", "Electron_pt[TightElectrons]")
                                          .Define("z_mass", inv_mass, {"z_lep_pt", "z_lep_eta", "z_lep_phi", "z_lep_mass"})
                                          .Define("RecoZ", RecoZ, RecoZstrings_ee)
                                          .Define("RecoZPt", TLorentzVectorPt, {"RecoZ"})
                                          .Define("RecoZPhi", TLorentzVectorPhi, {"RecoZ"})
                                          .Define("RecoZEta", TLorentzVectorEta, {"RecoZ"})
                                          .Define("dR_ll", deltaRcheck_float, {"LeadingElectronEta", "LeadingElectronPhi", "SubleadingElectronEta", "SubleadingElectronPhi"})
                                          .Define("dPhi_ll", DeltaPhi_floatandfloat, {"LeadingElectronPhi", "SubleadingElectronPhi"})
                                          .Filter(z_mass_cut, {"z_mass"}, "Z mass cut (ee channel)");

 std::cout << "before d_mumu_recoZ_selection" << std::endl;

 auto d_mumu_recoZ_selection = d_mumu_selection.Define("z_lep_eta", "Muon_eta[TightMuons]")
                                              .Define("z_lep_phi", "Muon_phi[TightMuons]")
                                              .Define("z_lep_mass", "Muon_mass[TightMuons]")
                                              .Define("z_lep_pt", "Muon_pt[TightMuons]")
                                              .Define("z_mass", inv_mass, {"z_lep_pt", "z_lep_eta", "z_lep_phi", "z_lep_mass"})
                                              .Define("RecoZ", RecoZ, RecoZstrings_mumu)
                                              .Define("RecoZPt", TLorentzVectorPt, {"RecoZ"})
                                              .Define("RecoZPhi", TLorentzVectorPhi, {"RecoZ"})
                                              .Define("RecoZEta", TLorentzVectorPt, {"RecoZ"})
                                              .Define("dR_ll", deltaRcheck_float, {"LeadingMuonPhi", "LeadingMuonEta", "SubleadingMuonEta", "SubleadingMuonPhi"})
                                              .Define("dPhi_ll", DeltaPhi_floatandfloat, {"LeadingMuonPhi", "SubleadingMuonPhi"})
                                              .Filter(z_mass_cut, {"z_mass"}, "Z mass cut (mumu channel)");


  std::cout << "before d_ee_recoZ_jets_selection" << std::endl;

  auto d_ee_recoZ_jets_selection = d_ee_recoZ_selection.Define("LeadingJetMass", LeadingVariable, {JetMassInput})
                                         	     .Define("SubleadingJetMass", SubleadingVariable, {JetMassInput})
                                         	     .Define("ThirdJetMass", ThirdLeadingVariable, {JetMassInput})
                                         	     .Define("FourthJetMass", FourthLeadingVariable, {JetMassInput})
					 	     .Define("LeadingJetpT", LeadingVariable, {JetPtInput})
					 	     .Define("SubleadingJetpT", SubleadingVariable, {JetPtInput})
					 	     .Define("ThirdJetpT", ThirdLeadingVariable, {JetPtInput})
                                         	     .Define("FourthJetpT", FourthLeadingVariable, {JetPtInput})
					 	     .Define("SumSquaredPt", SumSquared2LeadingJets_pT, {"LeadingJetpT", "SubleadingJetpT"})
					 	     .Define("JetPtSum", JetPtSum, {"LeadingJetpT", "SubleadingJetpT", "ThirdJetpT", "FourthJetpT"})
					 	     .Define("LeadingJetEta", LeadingVariable, {JetEtaInput})
                                         	     .Define("SubleadingJetEta", SubleadingVariable, {JetEtaInput})
                                                     .Define("ThirdJetEta", ThirdLeadingVariable, {JetEtaInput})
                                         	     .Define("FourthJetEta", FourthLeadingVariable, {JetEtaInput})
					 	     .Define("LeadingJetPhi", LeadingVariable, {JetPhiInput})
                                         	     .Define("SubleadingJetPhi", SubleadingVariable, {JetPhiInput})
                                         	     .Define("ThirdJetPhi", ThirdLeadingVariable, {JetPhiInput})
                                         	     .Define("FourthJetPhi", FourthLeadingVariable, {JetPhiInput})
					 	     .Define("dRJet_e", deltaRcheck_floats, deltaR_JetE_strings)
				         	     .Define("dR_j1j2", deltaRcheck_float, deltaR_j1j2_strings)
					 	     .Define("dPhi_j1j2", DeltaPhi_floatandfloat, {"LeadingJetPhi", "SubleadingJetPhi"})
					 	     .Define("LeadingJetHT", HT, {"LeadingJetpT"})
                                         	     .Define("SubleadingJetHT", HT, {"SubleadingJetpT"})
                                         	     .Define("ThirdJetHT", HT, {"ThirdJetpT"})
                                         	     .Define("FourthJetHT", HT, {"FourthJetpT"})
                                         	     .Define("TotJetHT", TotJetHT, {"LeadingJetHT", "SubleadingJetHT", "ThirdJetHT", "FourthJetHT"})
					 	     .Define("LeadingElectronHT", HT, {"LeadingElectron_pT"})
                                         	     .Define("SubleadingElectronHT", HT, {"SubleadingElectron_pT"})
					 	     .Define("TotLepHT", TotLepHT, {"LeadingElectronHT", "SubleadingElectronHT"})
					 	     .Define("TotHTOverTotpT_Jets", TotHTOverTotpT, {"TotJetHT", "JetPtSum"})
					 	     .Define("LepPtSum", LepPtSum, {"LeadingElectron_pT", "SubleadingElectron_pT"})
					 	     .Define("LepEtaSum", LepEtaSum, {"LeadingElectron_pT", "SubleadingElectron_pT"})
					 	     .Define("LepPhiSum", LepPhiSum, {"LeadingElectron_pT", "SubleadingElectron_pT"})
					 	     .Define("TotHTOverTotpT_Leptons", TotHTOverTotpT, {"TotLepHT", "LepPtSum"})
					 	     .Define("InvMassAllJets", InvMass_AllJets, InvMass_AllJets_strings)
					 	     .Define("InvMass3Jets", InvMass_3Jets, InvMass_3Jets_strings)
					 	     .Define("JetEtaSum", JetEtaSum, {"LeadingJetEta", "SubleadingJetEta", "ThirdJetEta", "FourthJetEta"})
                                         	     .Define("JetPhiSum", JetPhiSum, {"LeadingJetPhi", "SubleadingJetPhi", "ThirdJetPhi", "FourthJetPhi"})
					 	     .Define("tight_jets", tight_jets_function, {JetPtInput, JetEtaInput, "Jet_jetId", "dRJet_e"})
					 	     .Filter(jet_selection_function, {"tight_jets"}, "jet cut (ee channel)");
  

  std::cout << "before d_ee_recoZ_jets_bjets_selection" << std::endl;

  auto d_ee_recoZ_jets_bjets_selection = d_ee_recoZ_jets_selection.Define("bjets", bjet_id, {"tight_jets", "Jet_btagCSVV2", JetEtaInput})
                                                                  .Define("nbjets", numberofbjets, {"bjets"})
                                                                  .Filter(bjet_cut, {"bjets"}, "b jet cut (ee channel)");


  std::cout << "before d_mumu_recoZ_jets_selection" << std::endl; 

  auto d_mumu_recoZ_jets_selection = d_mumu_recoZ_selection.Define("LeadingJetMass", LeadingVariable, {JetMassInput})
                                                         .Define("SubleadingJetMass", SubleadingVariable, {JetMassInput})
                                                         .Define("ThirdJetMass", ThirdLeadingVariable, {JetMassInput})
                                                         .Define("FourthJetMass", FourthLeadingVariable, {JetMassInput})
					                 .Define("LeadingJetpT", LeadingVariable, {JetPtInput})
                                                         .Define("SubleadingJetpT", SubleadingVariable, {JetPtInput})
                                                         .Define("ThirdJetpT", ThirdLeadingVariable, {JetPtInput})
                                                         .Define("FourthJetpT", FourthLeadingVariable, {JetPtInput})
					                 .Define("SumSquaredPt", SumSquared2LeadingJets_pT, {"LeadingJetpT", "SubleadingJetpT"})
					                 .Define("JetPtSum", JetPtSum, {"LeadingJetpT", "SubleadingJetpT", "ThirdJetpT", "FourthJetpT"})
                                                         .Define("LeadingJetEta", LeadingVariable, {JetEtaInput})
                                                         .Define("SubleadingJetEta", SubleadingVariable, {JetEtaInput})
                                                         .Define("ThirdJetEta", ThirdLeadingVariable, {JetEtaInput})
                                                         .Define("FourthJetEta", FourthLeadingVariable, {JetEtaInput})
					                 .Define("LeadingJetPhi", LeadingVariable, {JetPhiInput})
                                                         .Define("SubleadingJetPhi", SubleadingVariable, {JetPhiInput})
                                                         .Define("ThirdJetPhi", ThirdLeadingVariable, {JetPhiInput})
                                                         .Define("FourthJetPhi", FourthLeadingVariable, {JetPhiInput})
                                                         .Define("dRJet_mu", deltaRcheck_floats, deltaR_JetMu_strings)
					                 .Define("dR_j1j2", deltaRcheck_float, deltaR_j1j2_strings)
					                 .Define("dPhi_j1j2", DeltaPhi_floatandfloat, {"LeadingJetPhi", "SubleadingJetPhi"})
					                 .Define("LeadingJetHT", HT, {"LeadingJetpT"})
					                 .Define("SubleadingJetHT", HT, {"SubleadingJetpT"})
					                 .Define("ThirdJetHT", HT, {"ThirdJetpT"})
					                 .Define("FourthJetHT", HT, {"FourthJetpT"})
					                 .Define("TotJetHT", TotJetHT, {"LeadingJetHT", "SubleadingJetHT", "ThirdJetHT", "FourthJetHT"})
                                                         .Define("LeadingMuonHT", HT, {"LeadingMuon_pT"})
                                                         .Define("SubleadingMuonHT", HT, {"SubleadingMuon_pT"})
                                                         .Define("TotLepHT", TotLepHT, {"LeadingMuonHT", "SubleadingMuonHT"})
					                 .Define("TotHTOverTotpT_Jets", TotHTOverTotpT, {"TotJetHT", "JetPtSum"})
                                                         .Define("LepPtSum", LepPtSum, {"LeadingMuon_pT", "SubleadingMuon_pT"})
					                 .Define("LepEtaSum", LepEtaSum, {"LeadingMuonEta", "SubleadingMuonEta"})
					                 .Define("LepPhiSum", LepPhiSum, {"LeadingMuonPhi", "SubleadingMuonPhi"})
					                 .Define("TotHTOverTotpT_Leptons", TotHTOverTotpT, {"TotLepHT", "LepPtSum"})
					                 .Define("InvMassAllJets", InvMass_AllJets, InvMass_AllJets_strings)
					                 .Define("InvMass3Jets", InvMass_3Jets, InvMass_3Jets_strings)
					                 .Define("JetEtaSum", JetEtaSum, {"LeadingJetEta", "SubleadingJetEta", "ThirdJetEta", "FourthJetEta"})
                                                         .Define("JetPhiSum", JetPhiSum, {"LeadingJetPhi", "SubleadingJetPhi", "ThirdJetPhi", "FourthJetPhi"})
				                         .Define("tight_jets", tight_jets_function, {JetPtInput, JetEtaInput, "Jet_jetId", "dRJet_mu"})
				                         .Filter(jet_selection_function, {"tight_jets"}, "jet cut (mumu channel)");

  std::cout << "before d_mumu_recoZ_jets_bjets_selection" << std::endl;

  auto d_mumu_recoZ_jets_bjets_selection = d_mumu_recoZ_jets_selection.Define("bjets", bjet_id, {"tight_jets", "Jet_btagCSVV2", JetEtaInput})
                                                                      .Define("nbjets", numberofbjets, {"bjets"})
                                                                      .Filter(bjet_cut, {"bjets"}, "b jet cut (mumu channel)"); 


  std::cout << "before Filtering events with a reconstructed W boson" << std::endl;

 //Filtering events with a reconstructed W boson
 auto d_ee_recoZ_jets_bjets_recoW_selection_defines = d_ee_recoZ_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", JetPtInput})
                 						.Define("w_reco_jets", find_w_pair, {JetPtInput, JetPhiInput, JetEtaInput, JetMassInput, "Jet_jetId", "lead_bjet"})
                 						.Define("w_pair_pt", select<floats>, {JetPtInput, "w_reco_jets"})
                 						.Define("w_pair_eta", select<floats>, {JetEtaInput, "w_reco_jets"})
                 						.Define("w_pair_phi", select<floats>, {JetPhiInput, "w_reco_jets"})
                 						.Define("w_pair_mass", select<floats>, {JetMassInput, "w_reco_jets"})
                 						.Define("w_mass", inv_mass, {"w_pair_pt", "w_pair_eta", "w_pair_phi", "w_pair_mass"})
								.Define("WPairJet1", WPairJet1, {JetPtInput, JetPhiInput, JetEtaInput, JetMassInput, "Jet_jetId", "lead_bjet"})
                 						.Define("WPairJet2", WPairJet2, {JetPtInput, JetPhiInput, JetEtaInput, JetMassInput, "Jet_jetId", "lead_bjet"})
								.Define("WPairJet1Pt", TLorentzVectorPt, {"WPairJet1"})
								.Define("WPairJet1Eta", TLorentzVectorEta, {"WPairJet1"})
								.Define("WPairJet1Phi", TLorentzVectorPhi, {"WPairJet1"})
								.Define("WPairJet1Mass", TLorentzVectorMass, {"WPairJet1"})
								.Define("WPairJet2Pt", TLorentzVectorPt, {"WPairJet2"})
                                                                .Define("WPairJet2Eta", TLorentzVectorEta, {"WPairJet2"})
                                                                .Define("WPairJet2Phi", TLorentzVectorPhi, {"WPairJet2"})
                                                                .Define("WPairJet2Mass", TLorentzVectorMass, {"WPairJet2"})
								.Define("dR_WJet1_WJet2", deltaRcheck_W_function, deltaR_WJet1_WJet2_strings)
								.Define("dWj1j2", DeltaPhi_function2, {"WPairJet1Phi", "WPairJet2Phi"})
								.Define("dR_WJet1_LeadingE", deltaRcheck_W_function2, deltaR_WJet1_LeadingElectron_strings)
								.Define("dR_WJet1_SubleadingE", deltaRcheck_W_function2, deltaR_WJet1_SubleadingElectron_strings)
								.Define("dR_WJet2_LeadingE", deltaRcheck_W_function2, deltaR_WJet2_LeadingElectron_strings)
                                                                .Define("dR_WJet2_SubleadingE", deltaRcheck_W_function2, deltaR_WJet2_SubleadingElectron_strings)
								.Define("dR_WJet1_LeadingJet", deltaRcheck_W_function2, deltaR_WJet1_LeadingJet_strings)
                                                                .Define("dR_WJet1_SubleadingJet", deltaRcheck_W_function2, deltaR_WJet1_SubleadingJet_strings)
                                                                .Define("dR_WJet1_ThirdJet", deltaRcheck_W_function2, deltaR_WJet1_ThirdJet_strings)
                                                                .Define("dR_WJet1_FourthJet", deltaRcheck_W_function2, deltaR_WJet1_FourthJet_strings)
								.Define("dR_WJet2_LeadingJet", deltaRcheck_W_function2, deltaR_WJet2_LeadingJet_strings)
                                                                .Define("dR_WJet2_SubleadingJet", deltaRcheck_W_function2, deltaR_WJet2_SubleadingJet_strings)
                                                                .Define("dR_WJet2_ThirdJet", deltaRcheck_W_function2, deltaR_WJet2_ThirdJet_strings)
                                                                .Define("dR_WJet2_FourthJet", deltaRcheck_W_function2, deltaR_WJet2_FourthJet_strings)
								.Define("dPhi_WJet1_LeadingE", DeltaPhi_doublesandfloat, {"WPairJet1Phi", "LeadingElectronPhi"})
                                                                .Define("dPhi_WJet1_SubleadingE", DeltaPhi_doublesandfloat, {"WPairJet1Phi", "SubleadingElectronPhi"})
                                                                .Define("dPhi_WJet2_LeadingE", DeltaPhi_doublesandfloat, {"WPairJet2Phi", "LeadingElectronPhi"})
                                                                .Define("dPhi_WJet2_SubleadingE", DeltaPhi_doublesandfloat, {"WPairJet2Phi", "SubleadingElectronPhi"})
                                                                .Define("dPhi_WJet1_LeadingJet", DeltaPhi_doublesandfloat, {"WPairJet1Phi", "LeadingJetPhi"})
                                                                .Define("dPhi_WJet1_SubleadingJet", DeltaPhi_doublesandfloat, {"WPairJet1Phi", "SubleadingJetPhi"})
                                                                .Define("dPhi_WJet1_ThirdJet", DeltaPhi_doublesandfloat, {"WPairJet1Phi", "ThirdJetPhi"})
                                                                .Define("dPhi_WJet1_FourthJet", DeltaPhi_doublesandfloat, {"WPairJet1Phi", "FourthJetPhi"})
                                                                .Define("dPhi_WJet2_LeadingJet", DeltaPhi_doublesandfloat, {"WPairJet2Phi", "LeadingJetPhi"})
                                                                .Define("dPhi_WJet2_SubleadingJet", DeltaPhi_doublesandfloat, {"WPairJet2Phi", "SubleadingJetPhi"})
                                                                .Define("dPhi_WJet2_ThirdJet", DeltaPhi_doublesandfloat, {"WPairJet2Phi", "ThirdJetPhi"})
                                                                .Define("dPhi_WJet2_FourthJet", DeltaPhi_doublesandfloat, {"WPairJet2Phi", "FourthJetPhi"})
								.Define("WJet1HT",  HT_double, {"WPairJet1Pt"})
                                                                .Define("WJet2HT",  HT_double, {"WPairJet2Pt"})
								.Define("RecoWHT", RecoWHT, {"w_pair_pt"})
								.Define("mtW", TransverseWMass, {"dPhi_j1j2", "WPairJet1Pt", "WPairJet2Pt"});


 std::cout << "before d_ee_recoZ_jets_bjets_recoW_selection" << std::endl;

auto d_ee_recoZ_jets_bjets_recoW_selection = d_ee_recoZ_jets_bjets_recoW_selection_defines.Filter(w_mass_cut, {"w_mass"}, "W mass cut (ee channel)");

 std::cout << "before ZPlusJetsCR == true" << std::endl;

if(ZPlusJetsCR == true){auto d_ee_recoZ_jets_bjets_recoW_selection = d_ee_recoZ_jets_bjets_recoW_selection_defines.Filter(w_mass_cut_ZPlusJetsCR, {"w_mass", "MET_sumEt"}, "W mass cut (ee channel)");}
else{auto d_ee_recoZ_jets_bjets_recoW_selection = d_ee_recoZ_jets_bjets_recoW_selection_defines.Filter(w_mass_cut, {"w_mass"}, "W mass cut (ee channel)");}


std::cout << "before d_mumu_recoZ_jets_bjets_recoW_selection_defines" << std::endl;


auto d_mumu_recoZ_jets_bjets_recoW_selection_defines = d_mumu_recoZ_jets_bjets_selection.Define("lead_bjet", find_lead_mask, {"bjets", JetPtInput})
							            .Define("w_reco_jets", find_w_pair, {JetPtInput, JetPhiInput, JetEtaInput, JetMassInput, "Jet_jetId", "lead_bjet"})
                                                                    .Define("w_pair_pt", select<floats>, {JetPtInput, "w_reco_jets"})
                                                                    .Define("w_pair_eta", select<floats>, {JetEtaInput, "w_reco_jets"})
                                                                    .Define("w_pair_phi", select<floats>, {JetPhiInput, "w_reco_jets"})
                                                                    .Define("w_pair_mass", select<floats>, {JetMassInput, "w_reco_jets"})
                                                                    .Define("w_mass", inv_mass, {"w_pair_pt", "w_pair_eta", "w_pair_phi", "w_pair_mass"})
								    .Define("WPairJet1", WPairJet1, {JetPtInput, JetPhiInput, JetEtaInput, JetMassInput, "Jet_jetId", "lead_bjet"})
                                                                    .Define("WPairJet2", WPairJet2, {JetPtInput, JetPhiInput, JetEtaInput, JetMassInput, "Jet_jetId", "lead_bjet"})
								    .Define("WPairJet1Pt", TLorentzVectorPt, {"WPairJet1"})
                                                                    .Define("WPairJet1Eta", TLorentzVectorEta, {"WPairJet1"})
                                                                    .Define("WPairJet1Phi", TLorentzVectorPhi, {"WPairJet1"})
                                                                    .Define("WPairJet1Mass", TLorentzVectorMass, {"WPairJet1"})
                                                                    .Define("WPairJet2Pt", TLorentzVectorPt, {"WPairJet2"})
                                                                    .Define("WPairJet2Eta", TLorentzVectorEta, {"WPairJet2"})
                                                                    .Define("WPairJet2Phi", TLorentzVectorPhi, {"WPairJet2"})
                                                                    .Define("WPairJet2Mass", TLorentzVectorMass, {"WPairJet2"})
								    .Define("dR_WJet1_WJet2", deltaRcheck_W_function, deltaR_WJet1_WJet2_strings)
                                                                    .Define("dWj1j2", DeltaPhi_function2, {"WPairJet1Phi", "WPairJet2Phi"})
                                                                    .Define("dR_WJet1_LeadingMu", deltaRcheck_W_function2, deltaR_WJet1_LeadingMuon_strings)
                                                                    .Define("dR_WJet1_SubleadingMu", deltaRcheck_W_function2, deltaR_WJet1_SubleadingMuon_strings)
                                                                    .Define("dR_WJet2_LeadingMu", deltaRcheck_W_function2, deltaR_WJet2_LeadingMuon_strings)
                                                                    .Define("dR_WJet2_SubleadingMu", deltaRcheck_W_function2, deltaR_WJet2_SubleadingMuon_strings)
                                                                    .Define("dR_WJet1_LeadingJet", deltaRcheck_W_function2, deltaR_WJet1_LeadingJet_strings)
                                                                    .Define("dR_WJet1_SubleadingJet", deltaRcheck_W_function2, deltaR_WJet1_SubleadingJet_strings)
                                                                    .Define("dR_WJet1_ThirdJet", deltaRcheck_W_function2, deltaR_WJet1_ThirdJet_strings)
                                                                    .Define("dR_WJet1_FourthJet", deltaRcheck_W_function2, deltaR_WJet1_FourthJet_strings)
                                                                    .Define("dR_WJet2_LeadingJet", deltaRcheck_W_function2, deltaR_WJet2_LeadingJet_strings)
                                                                    .Define("dR_WJet2_SubleadingJet", deltaRcheck_W_function2, deltaR_WJet2_SubleadingJet_strings)
                                                                    .Define("dR_WJet2_ThirdJet", deltaRcheck_W_function2, deltaR_WJet2_ThirdJet_strings)
                                                                    .Define("dR_WJet2_FourthJet", deltaRcheck_W_function2, deltaR_WJet2_FourthJet_strings)
                                                                    .Define("dPhi_WJet1_LeadingMu", DeltaPhi_doublesandfloat, {"WPairJet1Phi", "LeadingMuonPhi"})
                                                                    .Define("dPhi_WJet1_SubleadingMu", DeltaPhi_doublesandfloat, {"WPairJet1Phi", "SubleadingMuonPhi"})
                                                                    .Define("dPhi_WJet2_LeadingMu", DeltaPhi_doublesandfloat, {"WPairJet2Phi", "LeadingMuonPhi"})
                                                                    .Define("dPhi_WJet2_SubleadingMu", DeltaPhi_doublesandfloat, {"WPairJet2Phi", "SubleadingMuonPhi"})
                                                                    .Define("dPhi_WJet1_LeadingJet", DeltaPhi_doublesandfloat, {"WPairJet1Phi", "LeadingJetPhi"})
                                                                    .Define("dPhi_WJet1_SubleadingJet", DeltaPhi_doublesandfloat, {"WPairJet1Phi", "SubleadingJetPhi"})
                                                                    .Define("dPhi_WJet1_ThirdJet", DeltaPhi_doublesandfloat, {"WPairJet1Phi", "ThirdJetPhi"})
                                                                    .Define("dPhi_WJet1_FourthJet", DeltaPhi_doublesandfloat, {"WPairJet1Phi", "FourthJetPhi"})
                                                                    .Define("dPhi_WJet2_LeadingJet", DeltaPhi_doublesandfloat, {"WPairJet2Phi", "LeadingJetPhi"})
                                                                    .Define("dPhi_WJet2_SubleadingJet", DeltaPhi_doublesandfloat, {"WPairJet2Phi", "SubleadingJetPhi"})
                                                                    .Define("dPhi_WJet2_ThirdJet", DeltaPhi_doublesandfloat, {"WPairJet2Phi", "ThirdJetPhi"})
                                                                    .Define("dPhi_WJet2_FourthJet", DeltaPhi_doublesandfloat, {"WPairJet2Phi", "FourthJetPhi"})
								    .Define("WJet1HT",  HT_double, {"WPairJet1Pt"})
								    .Define("WJet2HT",  HT_double, {"WPairJet2Pt"})
								    .Define("RecoWHT", RecoWHT, {"w_pair_pt"})
								    .Define("mtW", TransverseWMass, {"dPhi_j1j2", "WPairJet1Pt", "WPairJet2Pt"});


std::cout << "before d_mumu_recoZ_jets_bjets_recoW_selection" << std::endl;


auto d_mumu_recoZ_jets_bjets_recoW_selection = d_mumu_recoZ_jets_bjets_recoW_selection_defines.Filter(w_mass_cut, {"w_mass"}, "W mass cut (mumu channel)");

if(ZPlusJetsCR == true){auto d_mumu_recoZ_jets_bjets_recoW_selection = d_mumu_recoZ_jets_bjets_recoW_selection_defines.Filter(w_mass_cut_ZPlusJetsCR, {"w_mass", "MET_sumEt"}, "W mass cut (mumu channel)");}
else{auto d_mumu_recoZ_jets_bjets_recoW_selection = d_mumu_recoZ_jets_bjets_recoW_selection_defines.Filter(w_mass_cut, {"w_mass"}, "W mass cut (mumu channel)");}


std::cout << "before filtering events with a reconstructed top quark" << std::endl;

//Filtering events with a reconstructed top quark
auto d_ee_recoZ_jets_bjets_recoW_recoT_selection = d_ee_recoZ_jets_bjets_recoW_selection.Define("RecoW", WLorentzVector, {"w_pair_pt", "w_pair_eta", "w_pair_phi", "w_mass", "w_reco_jets"})
											.Define("bjetmass", bjet_variable, bjet_mass_strings)
											.Define("bjetpt", bjet_variable, bjet_pt_strings)
											.Define("bjeteta", bjet_variable, bjet_eta_strings)
											.Define("bjetphi", bjet_variable, bjet_phi_strings)
											.Define("BJets", BLorentzVector, {"bjetpt", "bjeteta", "bjetphi", "bjetmass"})
											.Define("RecoTop", top_reconstruction_function, top_strings)
											.Define("Top_Pt", TLorentzVectorPt, {"RecoTop"})
											.Define("Top_Eta", TLorentzVectorEta, {"RecoTop"})
											.Define("Top_Phi", TLorentzVectorPhi, {"RecoTop"})
											.Define("Top_Mass", TLorentzVectorMass, {"RecoTop"})
											.Define("Top_HT", HT_double, {"Top_Pt"})
											.Define("dR_Top_LeadingElectron", deltaRcheck_Top_function, deltaR_Top_LeadingElectron_strings)
											.Define("dR_Top_SubleadingElectron", deltaRcheck_Top_function, deltaR_Top_SubleadingElectron_strings)
											.Define("dR_Top_LeadingJet", deltaRcheck_Top_function, deltaR_Top_LeadingJet_strings)
											.Define("dR_Top_SubleadingJet", deltaRcheck_Top_function, deltaR_Top_SubleadingJet_strings)
											.Define("dR_Top_ThirdJet", deltaRcheck_Top_function, deltaR_Top_ThirdJet_strings)
											.Define("dR_Top_FourthJet", deltaRcheck_Top_function, deltaR_Top_FourthJet_strings)
										        .Define("dR_Top_W", deltaRcheck_WTop_function, deltaR_Top_W_strings)
											.Define("dPhi_Wj1_Top", DeltaPhi_function2, {"WPairJet1Phi", "Top_Phi"})
											.Define("dPhi_Wj2_Top", DeltaPhi_function2, {"WPairJet2Phi", "Top_Phi"})
										        .Define("dR_Z_Top", deltaRcheck_W_function, deltaR_Z_Top_strings)
											.Define("dPhi_Z_Top", DeltaPhi_function2, {"Top_Phi", "RecoZPhi"})
											.Define("dR_Z_WPairJet1", deltaRcheck_W_function, dR_Z_WPairJet1_strings)
											.Define("dR_Z_WPairJet2", deltaRcheck_W_function, dR_Z_WPairJet2_strings)
											.Define("dPhi_Z_WPairJet1", DeltaPhi_function2, {"RecoZPhi", "WPairJet1Phi"})
                                                                                        .Define("dPhi_Z_WPairJet2", DeltaPhi_function2, {"RecoZPhi", "WPairJet2Phi"})
											.Define("MinDeltaR", MinDeltaR, MinDeltaR_strings)
											.Define("MinDeltaPhi", MinDeltaPhi, MinDeltaPhi_strings)
											.Define("dR_LeadingLepton_LeadingBJet", dR_LeadingLepton_LeadingBJet, dR_LeadingLepton_LeadingBJetstrings_ee)
											.Define("dR_SubleadingLepton_LeadingBJet", dR_SubleadingLepton_LeadingBJet, dR_SubleadingLepton_LeadingBJetstrings_ee)
											.Define("DeltaPhi_Leadinglepton_BJet", DeltaPhi_Lepton_BJet, DeltaPhi_Leadinglepton_BJet_string_ee)
                                                                                        .Define("DeltaPhi_Subleadinglepton_BJet", DeltaPhi_Lepton_BJet, DeltaPhi_Subleadinglepton_BJet_string_ee)		
										    	.Define("MET", MET_function, {"MET_sumEt"})
										        .Define("LeadingBJetOutputDiscriminant", LeadingBJetOutputDiscriminant, BJetOutputDiscriminantStrings_Leading)
                                                        			        .Define("SubleadingBJetOutputDiscriminant", SubleadingBJetOutputDiscriminant, BJetOutputDiscriminantStrings_Subleading)
                                                        				.Define("ThirdBJetOutputDiscriminant", ThirdBJetOutputDiscriminant, BJetOutputDiscriminantStrings_Third)
                                                        				.Define("FourthBJetOutputDiscriminant", FourthBJetOutputDiscriminant, BJetOutputDiscriminantStrings_Fourth)
                                                        				.Define("LeadingBJetOutputDiscriminant2", BJetOutputDiscriminant, {"LeadingBJetOutputDiscriminant", "Jet_btagCSVV2"})
                                                        				.Define("SubleadingBJetOutputDiscriminant2", BJetOutputDiscriminant, {"SubleadingBJetOutputDiscriminant", "Jet_btagCSVV2"})
                                                        				.Define("ThirdBJetOutputDiscriminant2", BJetOutputDiscriminant, {"ThirdBJetOutputDiscriminant", "Jet_btagCSVV2"})
                                                        				.Define("FourthBJetOutputDiscriminant2", BJetOutputDiscriminant, {"FourthBJetOutputDiscriminant", "Jet_btagCSVV2"})
                                                                                        .Define("dPhi_W_Top", DeltaPhi_function4, {"w_pair_phi", "Top_Phi"})
										        .Define("dR_Z_LeadingJet", deltaRcheck_W_function2, deltaR_Z_LeadingJet_strings)
                                                                            	        .Define("dR_Z_SubleadingJet", deltaRcheck_W_function2, deltaR_Z_SubleadingJet_strings)
                                                                            	        .Define("dR_Z_ThirdJet", deltaRcheck_W_function2, deltaR_Z_ThirdJet_strings)
                                                                            	        .Define("dR_Z_FourthJet", deltaRcheck_W_function2, deltaR_Z_FourthJet_strings)
                                                                            	        .Define("dPhi_LeadingJet_Z", DeltaPhi_doublesandfloat, {"RecoZPhi", "LeadingJetPhi"})
                                                                            		.Define("dPhi_SubleadingJet_Z", DeltaPhi_doublesandfloat, {"RecoZPhi", "SubleadingJetPhi"})
                                                                            	        .Define("dPhi_ThirdJet_Z", DeltaPhi_doublesandfloat, {"RecoZPhi", "ThirdJetPhi"})
                                                                            	        .Define("dPhi_FourthJet_Z", DeltaPhi_doublesandfloat, {"RecoZPhi", "FourthJetPhi"})
                                                                            	        .Define("dR_W_Z", deltaRcheck_WTop_function, deltaR_W_Z_strings)
                                                                            	        .Define("RecoZHT", HT_double, {"RecoZPt"})
                                                                            		.Define("dPhi_W_Z", DeltaPhi_function4, {"w_pair_phi", "RecoZPhi"})
											.Define("TotalHT_System", TotalHT_System, TotalHT_System_strings)
                                                                                        .Define("TotalPt_System", TotalPt_System, TotalPt_Systemstrings_data)
                                                                                        .Define("TotalEta_System", TotalEta_System, TotalEta_System_strings)
                                                                                        .Define("TotalPhi_System", TotalPhi_System, TotalPhi_System_strings)
                                                                                        .Define("TotHTOverTotpT_System", TotHTOverTotpT_floats, {"TotalHT_System", "TotalPt_System"})
											.Define("InvTopMass", inv_mass_doubles, {"Top_Pt", "Top_Eta", "Top_Phi", "Top_Mass"});



std::cout << "before d_mumu_recoZ_jets_bjets_recoW_recoT_selection" << std::endl;

auto d_mumu_recoZ_jets_bjets_recoW_recoT_selection = d_mumu_recoZ_jets_bjets_recoW_selection.Define("RecoW", WLorentzVector, {"w_pair_pt", "w_pair_eta", "w_pair_phi", "w_mass", "w_reco_jets"})
                                                                                            .Define("bjetmass", bjet_variable, bjet_mass_strings)
                                                                                            .Define("bjetpt", bjet_variable, bjet_pt_strings)
                                                                                            .Define("bjeteta", bjet_variable, bjet_eta_strings)
                                                                                            .Define("bjetphi", bjet_variable, bjet_phi_strings)
                                                                                            .Define("BJets", BLorentzVector, {"bjetpt", "bjeteta", "bjetphi", "bjetmass"})
											    .Define("RecoTop", top_reconstruction_function, top_strings)
											    .Define("Top_Pt", TLorentzVectorPt, {"RecoTop"})
											    .Define("Top_Eta", TLorentzVectorEta, {"RecoTop"})
                                                                                            .Define("Top_Phi", TLorentzVectorPhi, {"RecoTop"})
											    .Define("Top_Mass", TLorentzVectorMass, {"RecoTop"})
											    .Define("Top_HT", HT_double, {"Top_Pt"})
											    .Define("dR_Top_LeadingMuon", deltaRcheck_Top_function, deltaR_Top_LeadingMuon_strings)
                                                                                            .Define("dR_Top_SubleadingMuon", deltaRcheck_Top_function, deltaR_Top_SubleadingMuon_strings)
                                                                                            .Define("dR_Top_LeadingJet", deltaRcheck_Top_function, deltaR_Top_LeadingJet_strings)
                                                                                            .Define("dR_Top_SubleadingJet", deltaRcheck_Top_function, deltaR_Top_SubleadingJet_strings)
                                                                                            .Define("dR_Top_ThirdJet", deltaRcheck_Top_function, deltaR_Top_ThirdJet_strings)
                                                                                            .Define("dR_Top_FourthJet", deltaRcheck_Top_function, deltaR_Top_FourthJet_strings)
											    .Define("dR_Top_W", deltaRcheck_WTop_function, deltaR_Top_W_strings)
											    .Define("dPhi_Wj1_Top", DeltaPhi_function2, {"WPairJet1Phi", "Top_Phi"})
											    .Define("dPhi_Wj2_Top", DeltaPhi_function2, {"WPairJet2Phi", "Top_Phi"})
											    .Define("dR_Z_Top", deltaRcheck_W_function, deltaR_Z_Top_strings)
											    .Define("dPhi_Z_Top", DeltaPhi_function2, {"Top_Phi", "RecoZPhi"})
											    .Define("MinDeltaR", MinDeltaR, MinDeltaR_strings)
										            .Define("MinDeltaPhi", MinDeltaPhi, MinDeltaPhi_strings)
											    .Define("dR_LeadingLepton_LeadingBJet", dR_LeadingLepton_LeadingBJet, dR_LeadingLepton_LeadingBJetstrings_mumu)
                                                                                            .Define("dR_SubleadingLepton_LeadingBJet", dR_SubleadingLepton_LeadingBJet, dR_SubleadingLepton_LeadingBJetstrings_mumu)
											    .Define("DeltaPhi_Leadinglepton_BJet", DeltaPhi_Lepton_BJet, DeltaPhi_Leadinglepton_BJet_string_mumu)
											    .Define("DeltaPhi_Subleadinglepton_BJet", DeltaPhi_Lepton_BJet, DeltaPhi_Subleadinglepton_BJet_string_mumu)
											    .Define("MET", MET_function, {"MET_sumEt"})
											    .Define("LeadingBJetOutputDiscriminant", LeadingBJetOutputDiscriminant, BJetOutputDiscriminantStrings_Leading)
                                                        				    .Define("SubleadingBJetOutputDiscriminant", SubleadingBJetOutputDiscriminant, BJetOutputDiscriminantStrings_Subleading)
                                                        				    .Define("ThirdBJetOutputDiscriminant", ThirdBJetOutputDiscriminant, BJetOutputDiscriminantStrings_Third)
                                                        				    .Define("FourthBJetOutputDiscriminant", FourthBJetOutputDiscriminant, BJetOutputDiscriminantStrings_Fourth)
                                                        				    .Define("LeadingBJetOutputDiscriminant2", BJetOutputDiscriminant, {"LeadingBJetOutputDiscriminant", "Jet_btagCSVV2"})
                                                        				    .Define("SubleadingBJetOutputDiscriminant2", BJetOutputDiscriminant, {"SubleadingBJetOutputDiscriminant", "Jet_btagCSVV2"})
                                                        				    .Define("ThirdBJetOutputDiscriminant2", BJetOutputDiscriminant, {"ThirdBJetOutputDiscriminant", "Jet_btagCSVV2"})
                                                        				    .Define("FourthBJetOutputDiscriminant2", BJetOutputDiscriminant, {"FourthBJetOutputDiscriminant", "Jet_btagCSVV2"})
											    .Define("dPhi_W_Top", DeltaPhi_function4, {"w_pair_phi", "Top_Phi"})
											    .Define("dR_Z_LeadingJet", deltaRcheck_W_function2, deltaR_Z_LeadingJet_strings)
                                                                                	    .Define("dR_Z_SubleadingJet", deltaRcheck_W_function2, deltaR_Z_SubleadingJet_strings)
                                                                                	    .Define("dR_Z_ThirdJet", deltaRcheck_W_function2, deltaR_Z_ThirdJet_strings)
                                                                                	    .Define("dR_Z_FourthJet", deltaRcheck_W_function2, deltaR_Z_FourthJet_strings)
                                                                                	    .Define("dPhi_LeadingJet_Z", DeltaPhi_doublesandfloat, {"RecoZPhi", "LeadingJetPhi"})
                                                                                	    .Define("dPhi_SubleadingJet_Z", DeltaPhi_doublesandfloat, {"RecoZPhi", "SubleadingJetPhi"})
                                                                                	    .Define("dPhi_ThirdJet_Z", DeltaPhi_doublesandfloat, {"RecoZPhi", "ThirdJetPhi"})
                                                                                	    .Define("dPhi_FourthJet_Z", DeltaPhi_doublesandfloat, {"RecoZPhi", "FourthJetPhi"})
                                                                                	    .Define("dR_W_Z", deltaRcheck_WTop_function, deltaR_W_Z_strings)
                                                                                	    .Define("RecoZHT", HT_double, {"RecoZPt"})
                                                                                	    .Define("dR_Z_WPairJet1", deltaRcheck_W_function, dR_Z_WPairJet1_strings)
                                                                                	    .Define("dR_Z_WPairJet2", deltaRcheck_W_function, dR_Z_WPairJet2_strings)
                                                                                	    .Define("dPhi_Z_WPairJet1", DeltaPhi_function2, {"RecoZPhi", "WPairJet1Phi"})
                                                                                	    .Define("dPhi_Z_WPairJet2", DeltaPhi_function2, {"RecoZPhi", "WPairJet2Phi"})
                                                                                	    .Define("dPhi_W_Z", DeltaPhi_function4, {"w_pair_phi", "RecoZPhi"})
											    .Define("TotalHT_System", TotalHT_System, TotalHT_System_strings)
                                                                                            .Define("TotalPt_System", TotalPt_System, TotalPt_Systemstrings_data)
                                                                                            .Define("TotalEta_System", TotalEta_System, TotalEta_System_strings)
                                                                                            .Define("TotalPhi_System", TotalPhi_System, TotalPhi_System_strings)
                                                                                            .Define("TotHTOverTotpT_System", TotHTOverTotpT_floats, {"TotalHT_System", "TotalPt_System"})
											    .Define("InvTopMass", inv_mass_doubles, {"Top_Pt", "Top_Eta", "Top_Phi", "Top_Mass"});
												
												 

int nbins = 40;

std::cout << "before Unblinded data histograms start here" << std::endl;



//Unblinded data histograms start here
 auto h_WMass_ee = d_ee_recoZ_jets_bjets_recoW_recoT_selection.Histo1D({"h_WMass_ee", "Mass distribution of the W mass candidate (ee channel)", nbins, 0, 150}, "w_mass");
  auto h_InvTopMass_ee = d_ee_recoZ_jets_bjets_recoW_recoT_selection.Histo1D({"h_InvTopMass_ee", "Mass distribution of the top candidate (ee channel)", nbins, 0, 500}, "InvTopMass");
  auto h_WMass_mumu = d_mumu_recoZ_jets_bjets_recoW_recoT_selection.Histo1D({"h_WMass_mumu", "Mass distribution of the W mass candidate (mumu channel)", nbins, 0, 150}, "w_mass");
  auto h_InvTopMass_mumu = d_mumu_recoZ_jets_bjets_recoW_recoT_selection.Histo1D({"h_InvTopMass_mumu", "Mass distribution of the top candidate (mumu channel)", nbins, 0, 500}, "InvTopMass");

std::cout << "before chi2 for data" << std::endl;


//Chi^2 calculation using MC samples

double W_stddev_ee, Top_stddev_ee, W_stddev_mumu, Top_stddev_mumu;

std::string Filename;
std::string tZq_WAndTop_Filename; 





if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
        tZq_WAndTop_Filename = "tZq_AfterFullSelection_GaussianFit_NPL_" + year + ".root";
}
else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
        tZq_WAndTop_Filename = "tZq_AfterFullSelection_GaussianFit_ZPlusJetsCR_" + year + ".root";
}
else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
        tZq_WAndTop_Filename = "tZq_AfterFullSelection_GaussianFit_ttbarCR_" + year + ".root";
}
else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
         tZq_WAndTop_Filename = "tZq_AfterFullSelection_GaussianFit_NPL_ZPlusJetsCR_" + year + ".root";
}
else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
         tZq_WAndTop_Filename = "tZq_AfterFullSelection_GaussianFit_NPL_ttbarCR_" + year + ".root";
}
else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){std::cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << std::endl;}
else{tZq_WAndTop_Filename = "tZq_AfterFullSelection_GaussianFit_" + year + ".root";}

std::cout << "before tZq_WAndTop_File" << std::endl;

TFile* tZq_WAndTop_File = TFile::Open(tZq_WAndTop_Filename.c_str(), "READ");

TH1* WHist_ee = (TH1*)tZq_WAndTop_File->GetObjectChecked("h_WMass_ee", "TH1");
TH1* TopHist_ee = (TH1*)tZq_WAndTop_File->GetObjectChecked("h_InvTopMass_ee", "TH1");
TH1* WHist_mumu = (TH1*)tZq_WAndTop_File->GetObjectChecked("h_WMass_mumu", "TH1");
TH1* TopHist_mumu = (TH1*)tZq_WAndTop_File->GetObjectChecked("h_InvTopMass_mumu", "TH1");

W_stddev_ee  = WHist_ee->GetStdDev();
Top_stddev_ee = TopHist_ee->GetStdDev();
W_stddev_mumu  = WHist_mumu->GetStdDev();
Top_stddev_mumu = TopHist_mumu->GetStdDev();

tZq_WAndTop_File->Close();



std::cout << "before chi2_ee" << std::endl;


//Lambda function for chi squared calculation (calculated using MC but applied to both MC and data)
auto chi2_ee{[&W_stddev_ee, &Top_stddev_ee](const float& w_mass, const float& Top_Mass){

  std::cout << "inside chi2_ee for data" << std::endl;

  //calculating chi2
  float chi2 = pow(( (w_mass - W_MASS) / W_stddev_ee), 2) + pow(( (Top_Mass - TOP_MASS) / Top_stddev_ee), 2);
  return chi2;

}};


std::cout << "before chi2_mumu" << std::endl;


auto chi2_mumu{[&W_stddev_mumu, &Top_stddev_mumu](const float& w_mass, const float& Top_Mass){

  std::cout << "inside chi2_mumu for data" << std::endl;

  //calculating chi2
  float chi2 = pow(( (w_mass - W_MASS) / W_stddev_mumu), 2) + pow(( (Top_Mass - TOP_MASS) / Top_stddev_mumu), 2);
  return chi2;	

}};



std::cout << "before Chi2Range_string for data" << std::endl;

std::string Chi2Range_string, BlindedHistosRootFile;

if(blinding == true && (SBR == true || SR == true)){

    	std::cout << "before Blinding_ee for data" << std::endl;

	auto Blinding_ee =  d_ee_recoZ_jets_bjets_recoW_recoT_selection.Define("chi2", chi2_ee, {"w_mass", "InvTopMass"});
	auto Blinding_mumu =  d_mumu_recoZ_jets_bjets_recoW_recoT_selection.Define("chi2", chi2_mumu, {"w_mass", "InvTopMass"});


	std::cout << "before chi2_filter_ee for data" << std::endl;

	auto chi2_filter_ee{[&SBR, &SR](const float& chi2_ee){

		if(SBR == true){return chi2_ee < Chi2_SBR_ee && chi2_ee > Chi2_SR_ee;}
		else if(SR == true){return chi2_ee < Chi2_SR_ee;}
		else{std::cout << "SR and SBR cannot both be false" << std::endl;}

	}};

	
	std::cout << "before chi2_filter_mumu for data" << std::endl;

	auto chi2_filter_mumu{[&SBR, &SR](const float& chi2_mumu){
                
                if(SBR == true){return chi2_mumu < Chi2_SBR_mumu && chi2_mumu > Chi2_SR_mumu;}
                else if(SR == true){return chi2_mumu < Chi2_SR_mumu;} 
                else{std::cout << "SR and SBR cannot both be false" << std::endl;}
        
        }};


	std::cout << "before Blinding_ee_filtered for data" << std::endl; 

	auto Blinding_ee_filtered = Blinding_ee.Filter(chi2_filter_ee, {"chi2"}, "chi squared filter (ee)");
	auto Blinding_mumu_filtered = Blinding_mumu.Filter(chi2_filter_mumu, {"chi2"}, "chi squared filter (mumu)");


	std::cout << "before data snapshots" << std::endl;

	 //snapshots to save the histograms to output root files
        std::string OutRootFile_ee = "Results_" + process + "_" + year + "_ee_Blinded.root";
        std::string OutRootFile_mumu = "Results_" + process + "_" + year + "_mumu_Blinded.root";

	auto snapshot_ee = Blinding_ee_filtered.Snapshot("Events", OutRootFile_ee.c_str(), ".*");
	auto snapshot_mumu = Blinding_mumu_filtered.Snapshot("Events", OutRootFile_mumu.c_str(), ".*");
	

}
else{

 std::cout << "before OutRootFile_ee_unblinded for data" << std::endl;

 std::string OutRootFile_ee_unblinded = "Results_" + process + "_" + year + "_ee.root";
 std::string OutRootFile_mumu_unblinded = "Results_" + process + "_" + year + "_mumu.root";

 std::cout << "before snapshot_ee_unblinded" << std::endl;

 auto snapshot_ee_unblinded = d_ee_recoZ_jets_bjets_recoW_recoT_selection.Snapshot("Events", OutRootFile_ee_unblinded.c_str(), ".*");

 std::cout << "after snapshot_ee_unblinded" << std::endl;

 //auto snapshot_mumu_unblinded = d_mumu_recoZ_jets_bjets_recoW_recoT_selection.Snapshot("Events", OutRootFile_mumu_unblinded.c_str(), ".*");


 std::cout << "after snapshot_mumu_unblinded" << std::endl;

}




std::cout << "before Print cut report for data" << std::endl;

//Print cut report
auto allCutsReport{d_dataframe.Report()};


for (auto&& cutInfo: allCutsReport)

{CutFlowReport << cutInfo.GetName() << '\t' << cutInfo.GetAll() << '\t' << cutInfo.GetPass() << '\t' << cutInfo.GetEff() << " %" << std::endl;}



}
else{std::cout << "Please select data or MC" << std::endl;}





}










auto fulleventselection2(const bool& blinding, const bool& NPL, const bool& SR, const bool& SBR, const bool& ZPlusJetsCR, const bool& ttbarCR, const std::string& year, const bool& PU_ScaleUp, const bool& PU_ScaleDown, const bool& BTag_ScaleUp, const bool& BTag_ScaleDown, const bool& JetSmearing_ScaleUp, const bool& JetSmearing_ScaleDown, const bool& JetResolution_ScaleUp, const bool& JetResolution_ScaleDown, const bool& LeptonEfficiencies_ScaleUp, const bool& LeptonEfficiencies_ScaleDown, const bool& PDF_ScaleUp, const bool& PDF_ScaleDown, const bool& ME_Up, const bool& ME_Down, const bool& MET_Up, const bool& MET_Down, const bool& isr_up, const bool& isr_down, const bool& fsr_up, const bool& fsr_down){


  std::vector<std::string> Processes;

  if(year == "2016"){

        Processes = {/*"MC_triggerSF_ttbar", "Data_triggerSF", "tZq", "ZPlusJets_M50_aMCatNLO", "ZPlusJets_M10To50_aMCatNLO", "ZPlusJets_M10To50_aMCatNLO_ext", 
		     "ZPlusJets_M50_Madgraph", "ZPlusJets_M50_Madgraph_ext", "ZPlusJets_M10To50_Madgraph",
                     "ttbar_madgraph_NanoAODv5", "ttbar_aMCatNLO", "ttbar_inc", "SingleTop_schannel",
                     "SingleTop_tchannel_top", "SingleTop_tchannel_tbar", "SingleTop_tHq", "SingleTop_tW", "SingleTop_tbarW",
                     "SingleTop_tWZ_tWll","VV_ZZTo2L2Nu", "VV_ZZTo2L2Q", "VV_ZZTo4L", "VV_WZTo1L1Nu2Q",
                     "VV_WZTo2L2Q", "VV_WWTo2L2Nu", "VV_WWToLNuQQ", "VVV_WWWTo4F",
                     "VVV_WWZ", "VVV_WZZ", "VVV_ZZZ", "WPlusJets_WJetsToLNu", "ttbarV_ttWJetsToLNu", "ttbarV_ttWJetsToQQ", 
                     "ttbarV_ttHTobb", "ttbarV_ttHToNonbb", "ttbarV_ttZToLLNuNu", "ttbarV_ttZToLLNuNu_ext2", "ttbarV_ttZToLLNuNu_ext3", "ttbarV_ttZToQQ",
                     "TT_hdampUP", "TT_hdampUP_ext", "TT_hdampDOWN", "TT_hdampDOWN_ext", */ /*"ST_tchannel_top_hdampup", "ST_tchannel_top_hdampdown",
                     "ST_tchannel_top_ScaleUp", "ST_tchannel_top_ScaleDown", "tW_tbar_ScaleUp", "tW_tbar_ScaleDown", "tW_top_ScaleUp", "tW_top_ScaleDown",
                     "TT_isr_UP", "TT_isr_DOWN", "TT_isr_DOWN_ext", "TT_fsr_UP", "TT_fsr_UP_ext", "TT_fsr_DOWN", "TT_fsr_DOWN_ext",*/ "data_DoubleEGRunB",
                     "data_DoubleEGRunC", "data_DoubleEGRunD", "data_DoubleEGRunE", "data_DoubleEGRunF", "data_DoubleEGRunG", "data_DoubleEGRunH",
                     "data_SingleElectronRunB", "data_SingleElectronRunC", "data_SingleElectronRunD", "data_SingleElectronRunE", "data_SingleElectronRunF",
                     "data_SingleElectronRunG", "data_SingleElectronRunH", "data_DoubleMuonRunB", "data_DoubleMuonRunC", "data_DoubleMuonRunD",
                     "data_DoubleMuonRunE", "data_SingleElectronRunF", "data_SingleElectronRunG", "data_DoubleMuonRunH", "data_SingleMuonRunB",
                     "data_SingleMuonRunC", "data_SingleMuonRunD", "data_SingleMuonRunE", "data_SingleMuonRunF", "data_SingleMuonRunG",
                     "data_SingleMuonRunH"};

 }
  else if(year == "2017"){

  	Processes = {/*"MC_triggerSF_ttbar", "Data_triggerSF", "tZq", "ZPlusJets_M50_aMCatNLO", "ZPlusJets_M50_aMCatNLO_ext", "ZPlusJets_M10To50_Madgraph", "ttbar_2l2nu",
		     "ttbar_madgraph_NanoAODv5", "ttbar_TTToHadronic", "ttbar_TTToSemileptonic", "ttbar_aMCatNLO", "SingleTop_schannel",
	      	     "SingleTop_tchannel_top", "SingleTop_tchannel_tbar", "SingleTop_tHq", "SingleTop_tW", "SingleTop_tbarW",
	             "SingleTop_tZq_W_lept_Z_had", "SingleTop_tWZ_tWll", "VV_ZZTo2Q2Nu", "VV_ZZTo2L2Nu", "VV_ZZTo2L2Q", "VV_ZZTo4L", "VV_WZTo1L1Nu2Q", 
	             "VV_WZTo2L2Q", "VV_WZTo3LNu", "VV_WWTo1L1Nu2Q", "VV_WWTo2L2Nu", "VV_WWToLNuQQ", "VV_WGToLNuG", "VV_ZGToLLG", "VVV_WWWTo4F", 
	             "VVV_WWZTo4F", "VVV_WZZ", "VVV_ZZZ", "WPlusJets_WJetsToLNu", "ttbarV_ttWJetsToLNu", "ttbarV_ttWJetsToQQ", "ttbarV_ttgamma",  
	             "ttbarV_ttZToLL", "ttbarV_ttHTobb", "ttbarV_ttHToNonbb", "ttbarV_ttZToLLNuNu", "ttbarV_ttZToQQ", "ttbarV_ttZToQQ_ext",
		     "data_DoubleEGRunB", "data_DoubleEGRunC", "data_DoubleEGRunD", "data_DoubleEGRunE", "data_DoubleEGRunF", "data_SingleElectronRunB", 
		     "data_SingleElectronRunC", "data_SingleElectronRunD", "data_SingleElectronRunE", "data_SingleElectronRunF", "data_DoubleMuonRunB", 
		     "data_DoubleMuonRunC", "data_DoubleMuonRunD", "data_DoubleMuonRunE", "data_DoubleMuonRunF", "data_SingleMuonRunB",*/ 
		     "data_SingleMuonRunC", "data_SingleMuonRunD", "data_SingleMuonRunE", "data_SingleMuonRunF"};

 }
 else if(year == "2018"){

	Processes = {"MC_triggerSF_ttbar", "Data_triggerSF", "tZq", "ZPlusJets_M50_aMCatNLO", "ZPlusJets_M50_aMCatNLO_ext", "ZPlusJets_M10To50_Madgraph", "ttbar_2l2nu",
                     "ttbar_madgraph_NanoAODv5", "ttbar_TTToHadronic", "ttbar_TTToSemileptonic", "ttbar_aMCatNLO", "SingleTop_schannel",
                     "SingleTop_tchannel_top", "SingleTop_tchannel_tbar", "SingleTop_tHq", "SingleTop_tW", "SingleTop_tbarW",
                     "SingleTop_tZq_W_lept_Z_had", "SingleTop_tWZ_tWll", "VV_ZZTo2Q2Nu", "VV_ZZTo2L2Nu", "VV_ZZTo2L2Q", "VV_ZZTo4L", "VV_WZTo1L1Nu2Q",
                     "VV_WZTo2L2Q", "VV_WZTo3LNu", "VV_WWTo1L1Nu2Q", "VV_WWTo2L2Nu", "VV_WWToLNuQQ", "VV_WGToLNuG", "VV_ZGToLLG", "VVV_WWWTo4F",
                     "VVV_WWZTo4F", "VVV_WZZ", "VVV_ZZZ", "WPlusJets_WJetsToLNu", "ttbarV_ttWJetsToLNu", "ttbarV_ttWJetsToQQ", "ttbarV_ttgamma",
                     "ttbarV_ttZToLL", "ttbarV_ttHTobb", "ttbarV_ttHToNonbb", "ttbarV_ttZToLLNuNu", "ttbarV_ttZToQQ", "ttbarV_ttZToQQ_ext", 
		     "data_EGRunB", "data_EGRunC", "data_EGRunD", "data_DoubleMuonRunB", "data_DoubleMuonRunC", "data_DoubleMuonRunD", "data_SingleMuonRunB", 
		     "data_SingleMuonRunC", "data_SingleMuonRunD"};


 }
 else{std::cout << "Error: Choose a year out of 2016, 2017 or 2018" << std::endl;}



//looping over the process names
  for(int i = 0; i < Processes.size(); i++){

  	fulleventselection_calculator(Processes.at(i), blinding, NPL, SR, SBR, ZPlusJetsCR, ttbarCR, year, PU_ScaleUp, PU_ScaleDown, BTag_ScaleUp, BTag_ScaleDown, JetSmearing_ScaleUp, JetSmearing_ScaleDown, JetResolution_ScaleUp, JetResolution_ScaleDown, LeptonEfficiencies_ScaleUp, LeptonEfficiencies_ScaleDown, PDF_ScaleUp, PDF_ScaleDown, ME_Up, ME_Down, MET_Up, MET_Down, isr_up, isr_down, fsr_up, fsr_down);

  }




}


//NPL background estimation
auto Hadd(const std::string& year, const bool& blinding){

  //tZq (signal), tHq, ttW, ttZ, WZ are used for the ratio

  if(blinding == true){

 	if(year == "2016"){

                system("hadd Results_MCRatio_2016_ee_NPL_Blinded.root Results_t*q_2016_ee_NPL_Blinded.root Results_ttbarV_ttbarW_*_2016_ee_NPL_Blinded.root Results_ttbarV_ttbarZ_*_2016_ee_NPL_Blinded.root Results_Diboson_WZ_*_2016_ee_NPL_Blinded.root");
	
		system("hadd Results_MCRatio_2016_mumu_NPL_Blinded.root Results_t*q_2016_mumu_NPL_Blinded.root Results_ttbarV_ttbarW_*_2016_mumu_NPL_Blinded.root Results_ttbarV_ttbarZ_*_2016_mumu_NPL_Blinded.root Results_Diboson_WZ_*_2016_mumu_NPL_Blinded.root");

		system("hadd Results_AllMC_2016_ee_NPL_Blinded.root Results_MC_*_2016_ee_NPL_Blinded.root");

		system("hadd Results_AllMC_2016_mumu_NPL_Blinded.root Results_MC_*_2016_mumu_NPL_Blinded.root");

		system("hadd Results_AllData_2016_ee_NPL_Blinded.root Results_data_*_2016_ee_NPL_Blinded.root");

                system("hadd Results_AllData_2016_mumu_NPL_Blinded.root Results_data_*_2016_mumu_NPL_Blinded.root");

        }
        else if(year == "2017"){
         
	        system("hadd Results_MCRatio_2017_ee_NPL_Blinded.root Results_t*q_2017_ee_NPL_Blinded.root Results_ttbarV_ttbarW_*_2017_ee_NPL_Blinded.root Results_ttbarV_ttbarZ_*_2017_ee_NPL_Blinded.root Results_Diboson_WZ_*_2017_ee_NPL_Blinded.root");

		system("hadd Results_MCRatio_2017_mumu_NPL_Blinded.root Results_t*q_2017_mumu_NPL_Blinded.root Results_ttbarV_ttbarW_*_2017_mumu_NPL_Blinded.root Results_ttbarV_ttbarZ_*_2017_mumu_NPL_Blinded.root Results_Diboson_WZ_*_2017_mumu_NPL_Blinded.root");        

		system("hadd Results_AllMC_2017_ee_NPL_Blinded.root Results_MC_*_2017_ee_NPL_Blinded.root");

                system("hadd Results_AllMC_2017_mumu_NPL_Blinded.root Results_MC_*_2017_mumu_NPL_Blinded.root");

		system("hadd Results_AllData_2017_ee_NPL_Blinded.root Results_data_*_2017_ee_NPL_Blinded.root");

                system("hadd Results_AllData_2017_mumu_NPL_Blinded.root Results_data_*_2017_mumu_NPL_Blinded.root");

	}
        else if(year == "2018"){
        
	        system("hadd Results_MCRatio_2018_ee_NPL_Blinded.root Results_t*q_2018_ee_NPL_Blinded.root Results_ttbarV_ttbarW_*_2018_ee_NPL_Blinded.root Results_ttbarV_ttbarZ_*_2018_ee_NPL_Blinded.root Results_Diboson_WZ_*_2018_ee_NPL_Blinded.root");

		system("hadd Results_MCRatio_2018_mumu_NPL_Blinded.root Results_t*q_2018_mumu_NPL_Blinded.root Results_ttbarV_ttbarW_*_2018_mumu_NPL_Blinded.root Results_ttbarV_ttbarZ_*_2018_mumu_NPL_Blinded.root Results_Diboson_WZ_*_2018_mumu_NPL_Blinded.root");

		system("hadd Results_AllMC_2018_ee_NPL_Blinded.root Results_MC_*_2018_ee_NPL_Blinded.root");

                system("hadd Results_AllMC_2018_mumu_NPL_Blinded.root Results_MC_*_2018_mumu_NPL_Blinded.root");

		system("hadd Results_AllData_2018_ee_NPL_Blinded.root Results_data_*_2018_ee_NPL_Blinded.root");

                system("hadd Results_AllData_2018_mumu_NPL_Blinded.root Results_data_*_2018_mumu_NPL_Blinded.root");

        }
        else{std::cout << "Choose a year out of 2016, 2017 or 2018" << std::endl;} 


  }
  else{

  	if(year == "2016"){
 
                system("hadd Results_MCRatio_2016_ee_NPL.root Results_t*q_2016_ee_NPL.root Results_ttbarV_ttbarW_*_2016_ee_NPL.root Results_ttbarV_ttbarZ_*_2016_ee_NPL.root Results_Diboson_WZ_*_2016_ee_NPL.root");
 
		system("hadd Results_MCRatio_2016_mumu_NPL.root Results_t*q_2016_mumu_NPL.root Results_ttbarV_ttbarW_*_2016_mumu_NPL.root Results_ttbarV_ttbarZ_*_2016_mumu_NPL.root Results_Diboson_WZ_*_2016_mumu_NPL.root"); 

		system("hadd Results_AllMC_2016_ee_NPL.root Results_MC_*_2016_ee_NPL.root");

                system("hadd Results_AllMC_2016_mumu_NPL.root Results_MC_*_2016_mumu_NPL.root");

		system("hadd Results_AllData_2016_ee_NPL.root Results_data_*_2016_ee_NPL.root");

                system("hadd Results_AllData_2016_mumu_NPL.root Results_data_*_2016_mumu_NPL.root");

         }
        else if(year == "2017"){
        
	        system("hadd Results_MCRatio_2017_ee_NPL.root Results_t*q_2017_ee_NPL.root Results_ttbarV_ttbarW_*_2017_ee_NPL.root Results_ttbarV_ttbarZ_*_2017_ee_NPL.root Results_Diboson_WZ_*_2017_ee_NPL.root");
     
		system("hadd Results_MCRatio_2017_mumu_NPL.root Results_t*q_2017_mumu_NPL.root Results_ttbarV_ttbarW_*_2017_mumu_NPL.root Results_ttbarV_ttbarZ_*_2017_mumu_NPL.root Results_Diboson_WZ_*_2017_mumu_NPL.root");

		system("hadd Results_AllMC_2017_ee_NPL.root Results_MC_*_2017_ee_NPL.root");

                system("hadd Results_AllMC_2017_mumu_NPL.root Results_MC_*_2017_mumu_NPL.root");

		system("hadd Results_AllData_2017_ee_NPL.root Results_data_*_2017_ee_NPL.root");

                system("hadd Results_AllData_2017_mumu_NPL.root Results_data_*_2017_mumu_NPL.root");

        }
        else if(year == "2018"){

                system("hadd Results_MCRatio_2018_ee_NPL.root Results_t*q_2018_ee_NPL.root Results_ttbarV_ttbarW_*_2018_ee_NPL.root Results_ttbarV_ttbarZ_*_2018_ee_NPL.root Results_Diboson_WZ_*_2018_ee_NPL.root");

		system("hadd Results_MCRatio_2018_mumu_NPL.root Results_t*q_2018_mumu_NPL.root Results_ttbarV_ttbarW_*_2018_mumu_NPL.root Results_ttbarV_ttbarZ_*_2018_mumu_NPL.root Results_Diboson_WZ_*_2018_mumu_NPL.root");

		system("hadd Results_AllMC_2018_ee_NPL.root Results_MC_*_2018_ee_NPL.root");

                system("hadd Results_AllMC_2018_mumu_NPL.root Results_MC_*_2018_mumu_NPL.root");

		system("hadd Results_AllData_2018_ee_NPL.root Results_data_*_2018_ee_NPL.root");

                system("hadd Results_AllData_2018_mumu_NPL.root Results_data_*_2018_mumu_NPL.root");


        }
        else{std::cout << "Choose a year out of 2016, 2017 or 2018" << std::endl;}


  }


}





auto NPLROOTFile_Creator2(const std::string& year, const bool& blinding){

 Hadd(year, blinding);

 TFile* AllMC_ee;
 TFile* AllMC_mumu;
 TFile* AllData_ee;
 TFile* AllData_mumu;
 TFile* MCRatio_ee;
 TFile* MCRatio_mumu;

 std::string AllMC_ee_File;
 std::string AllMC_mumu_File;
 std::string AllData_ee_File;
 std::string AllData_mumu_File;
 std::string MCRatio_ee_File;
 std::string MCRatio_mumu_File;

 if(blinding == true){

 	AllMC_ee_File = "Results_AllMC_" + year + "_ee_NPL_Blinded.root";
	AllMC_mumu_File = "Results_AllMC_" + year + "_mumu_NPL_Blinded.root";
	AllData_ee_File = "Results_AllData_" + year + "_ee_NPL_Blinded.root";
	AllData_mumu_File = "Results_AllData_" + year + "_mumu_NPL_Blinded.root";
	MCRatio_ee_File = "Results_MCRatio_" + year + "_ee_NPL_Blinded.root";
	MCRatio_mumu_File = "Results_MCRatio_" + year + "_mumu_NPL_Blinded.root";
 }
 else{

	AllMC_ee_File = "Results_AllMC_" + year + "_ee_NPL.root";
        AllMC_mumu_File = "Results_AllMC_" + year + "_mumu_NPL.root";
        AllData_ee_File = "Results_AllData_" + year + "_ee_NPL.root";
        AllData_mumu_File = "Results_AllData_" + year + "_mumu_NPL.root";
        MCRatio_ee_File = "Results_MCRatio_" + year + "_ee_NPL.root";
        MCRatio_mumu_File = "Results_MCRatio_" + year + "_mumu_NPL.root";

 }
 
 AllMC_ee = new TFile(AllMC_ee_File.c_str(), "READ");
 AllMC_mumu = new TFile(AllMC_mumu_File.c_str(), "READ");
 AllData_ee = new TFile(AllData_ee_File.c_str(), "READ");
 AllData_mumu = new TFile(AllData_mumu_File.c_str(), "READ");
 MCRatio_ee = new TFile(MCRatio_ee_File.c_str(), "READ");
 MCRatio_mumu = new TFile(MCRatio_mumu_File.c_str(), "READ");


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

 int NMC_OS_NonPrompt_ee = h_NMC_OS_NonPrompt_ee->GetEntries();
 int NMC_OS_NonPrompt_mumu = h_NMC_OS_NonPrompt_mumu->GetEntries();
 int NMC_SS_NonPrompt_ee = h_NMC_SS_NonPrompt_ee->GetEntries();
 int NMC_SS_NonPrompt_mumu = h_NMC_SS_NonPrompt_mumu->GetEntries();

 int ratio_ee = NMC_OS_NonPrompt_ee / NMC_SS_NonPrompt_ee;
 int ratio_mumu = NMC_OS_NonPrompt_mumu / NMC_SS_NonPrompt_mumu;

 //Calculating the number of opposite-sign non prompt data entries
 TH1* h_NData_OS_NonPrompt_ee; 
 TH1* h_NData_OS_NonPrompt_mumu;

 int nbins = h_NData_SS_ee->GetNbinsX();

 for(int i = 0; i < nbins; i++){

	double ee_content = ( ( h_NData_SS_ee->GetBinContent(i) - h_NMC_SS_ee->GetBinContent(i) ) * ratio_ee);
	double mumu_content = ( ( h_NData_SS_mumu->GetBinContent(i) - h_NMC_SS_mumu->GetBinContent(i) ) * ratio_mumu);

	h_NData_OS_NonPrompt_ee->SetBinContent(i,ee_content);
   	h_NData_OS_NonPrompt_mumu->SetBinContent(i, mumu_content);

 }


 //Saving the histograms to an output file
 std::string NPL_output_file;

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






//Main function
void fulleventselectionAlgo::fulleventselection(){

  time_t now = time(0);
  tm* localtm = localtime(&now);
  std::cout << "The script started running:" << " " << asctime(localtm) << std::endl;


//  fulleventselection2(blinding, NPL, SR, SBR, ZPlusJetsCR, ttbarCR, year, PU_ScaleUp, PU_ScaleDown, BTag_ScaleUp, BTag_ScaleDown, JetSmearing_ScaleUp, JetSmearing_ScaleDown, JetResolution_ScaleUp, JetResolution_ScaleDown, LeptonEfficiencies_ScaleUp, LeptonEfficiencies_ScaleDown, PDF_ScaleUp, PDF_ScaleDown, ME_Up, ME_Down, MET_Up, MET_Down, isr_up, isr_down, fsr_up, fsr_down);


  bool blinding = false;
  std::vector<bool> NPL = {false/*, true*/};
//  std::vector<bool> ZPlusJetsCR = {false, true}; 
//  std::vector<bool> ttbarCR = {false, true};

  bool SR = true;
  bool SBR = false;
  bool ZPlusJetsCR = false;
  bool ttbarCR = false;

  std::vector<std::string> year = {/*"2016", */"2017"/*, "2018"*/};
  

  for(int i = 0; i < year.size(); i++){

	for(int j = 0; j < NPL.size(); j++){

  		//Nominal
  		fulleventselection2(blinding, NPL.at(j), SR, SBR, ZPlusJetsCR, ttbarCR, year.at(i), false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false);

 			if(NPL.at(i) == false){

				//PU_ScaleUp
				fulleventselection2(blinding, NPL.at(j), SR, SBR, ZPlusJetsCR, ttbarCR, year.at(i), true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false);

				//PU_ScaleDown
				fulleventselection2(blinding, NPL.at(j), SR, SBR, ZPlusJetsCR, ttbarCR, year.at(i), false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false);

				//BTag_ScaleUp
				fulleventselection2(blinding, NPL.at(j), SR, SBR, ZPlusJetsCR, ttbarCR, year.at(i), false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false);

				//BTag_ScaleDown
				fulleventselection2(blinding, NPL.at(j), SR, SBR, ZPlusJetsCR, ttbarCR, year.at(i), false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false);

				//JetSmearing_ScaleUp
				fulleventselection2(blinding, NPL.at(j), SR, SBR, ZPlusJetsCR, ttbarCR, year.at(i), false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false);

				//JetSmearing_ScaleDown
				fulleventselection2(blinding, NPL.at(j), SR, SBR, ZPlusJetsCR, ttbarCR, year.at(i), false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false, false);

				//JetResolution_ScaleUp
				fulleventselection2(blinding, NPL.at(j), SR, SBR, ZPlusJetsCR, ttbarCR, year.at(i), false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false, false);

				//JetResolution_ScaleDown
				fulleventselection2(blinding, NPL.at(j), SR, SBR, ZPlusJetsCR, ttbarCR, year.at(i), false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false, false);

				//LeptonEfficiencies_ScaleUp
				fulleventselection2(blinding, NPL.at(j), SR, SBR, ZPlusJetsCR, ttbarCR, year.at(i), false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false, false);

				//LeptonEfficiencies_ScaleDown
				fulleventselection2(blinding, NPL.at(j), SR, SBR, ZPlusJetsCR, ttbarCR, year.at(i), false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false, false);

				//PDF_ScaleUp 
        			fulleventselection2(blinding, NPL.at(j), SR, SBR, ZPlusJetsCR, ttbarCR, year.at(i), false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false, false);

        			//PDF_ScaleDown 
        			fulleventselection2(blinding, NPL.at(j), SR, SBR, ZPlusJetsCR, ttbarCR, year.at(i), false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false, false);

				//ME_Up
				fulleventselection2(blinding, NPL.at(j), SR, SBR, ZPlusJetsCR, ttbarCR, year.at(i), false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false, false);

				//ME_Down
				fulleventselection2(blinding, NPL.at(j), SR, SBR, ZPlusJetsCR, ttbarCR, year.at(i), false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false);

				//MET_Up
                        	fulleventselection2(blinding, NPL.at(j), SR, SBR, ZPlusJetsCR, ttbarCR, year.at(i), false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false, false);

                        	//MET_Down
                        	fulleventselection2(blinding, NPL.at(j), SR, SBR, ZPlusJetsCR, ttbarCR, year.at(i), false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false, false);	


   				if(year.at(i) == "2017" || year.at(i) == "2018"){

					//isr_up
					fulleventselection2(blinding, NPL.at(j), SR, SBR, ZPlusJetsCR, ttbarCR, year.at(i), false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false, false);

					//isr_down
					fulleventselection2(blinding, NPL.at(j), SR, SBR, ZPlusJetsCR, ttbarCR, year.at(i), false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false, false);

					//fsr_up
					fulleventselection2(blinding, NPL.at(j), SR, SBR, ZPlusJetsCR, ttbarCR, year.at(i), false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true, false);

					//fsr_down  
					fulleventselection2(blinding, NPL.at(j), SR, SBR, ZPlusJetsCR, ttbarCR, year.at(i), false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, true);

				}



  			}
  			else{


				//Creating the NPL root file
				NPLROOTFile_Creator2(year.at(i), blinding);


				//Running over the NPL root file
				//Nominal
				if(blinding == true){

					fulleventselection_calculator("NPL_File_ee_Blinded", blinding, NPL.at(j), SR, SBR, ZPlusJetsCR, ttbarCR, year.at(i), false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false);

					fulleventselection_calculator("NPL_File_mumu_Blinded", blinding, NPL.at(j), SR, SBR, ZPlusJetsCR, ttbarCR, year.at(i), false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false);

				}
				else{

					fulleventselection_calculator("NPL_File_ee_Unblinded", blinding, NPL.at(j), SR, SBR, ZPlusJetsCR, ttbarCR, year.at(i), false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false);

        				fulleventselection_calculator("NPL_File_mumu_Unblinded", blinding, NPL.at(j), SR, SBR, ZPlusJetsCR, ttbarCR, year.at(i), false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false);

				}



			}



  	}//end of for loop for NPL

  }//end of for loop for year



  //Saving the outputs to a directory
  for(int i = 0; i < year.size(); i++){

	for(int j = 0; j < NPL.size(); j++){

		DirectoryCreator(year.at(i), blinding, NPL.at(j));

	}

  }

 

  //Printing out the time the script finished running
  time_t now2 = time(0);
  tm* localtm2 = localtime(&now2);
  std::cout << "The script finished running at: " << asctime(localtm2) << std::endl;



}


