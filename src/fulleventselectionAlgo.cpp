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

constexpr double EndcapMinEta = 1.566;
constexpr double BarrelMaxEta = 1.4442;
double MaxTrackerEta;
double MinElectronPt;
double MinMuonPt;
double MinElectronPtEmu;
double MinMuonPtEmu;
double MaxElectronPt;
double MaxMuonPt;
constexpr float Z_MASS{91.1876f};
constexpr float Z_MASS_CUT{20.f}; 
constexpr float W_MASS = 80.385f;
constexpr float W_MASS_CUT = 20.f;
constexpr float TOP_MASS = 173.3;
float Chi2_SR;
float Chi2_SBR;


template<typename T>
[[gnu::const]] T select(const T& a, const ints& mask)
{
  return a[mask];
}

template<typename T, typename U> //for the all equal function
[[gnu::const]] bool all_equal(const T& t, const U& u)
{
  return t == u;
}

template<typename T, typename U, typename... Types>
[[gnu::const]] bool all_equal(const T& t, const U& u, Types const&... args)
{
    return t == u && all_equal(u, args...);
}



std::fstream& GotoLine(std::fstream& file, unsigned int num){
		
  file.seekg(std::ios::beg);

  for(unsigned int i =0; i < num - 1; ++i){
     file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  }

  return file;

}




//For b tagging reweighting
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



float N_SelectionCriteria_ee_MC, N_SelectionCriteria_mumu_MC, N_SelectionCriteria_emu_MC;
float N_MET_And_LeptonSelection_ee_MC, N_MET_And_LeptonSelection_mumu_MC, N_MET_And_LeptonSelection_emu_MC;
float N_LeptonTriggersAndSelectionCriteria_ee_MC, N_LeptonTriggersAndSelectionCriteria_mumu_MC, N_LeptonTriggersAndSelectionCriteria_emu_MC;
float N_MET_LeptonTriggers_SelectionCriteria_ee_MC, N_MET_LeptonTriggers_SelectionCriteria_mumu_MC, N_MET_LeptonTriggers_SelectionCriteria_emu_MC;
float Eff_ee_MC, Eff_mumu_MC, Eff_emu_MC;
float Eff_MET_LeptonTriggers_SelectionCriteria_ee_MC, Eff_MET_LeptonTriggers_SelectionCriteria_mumu_MC, Eff_MET_LeptonTriggers_SelectionCriteria_emu_MC;
float Eff_LeptonTriggers_SelectionCriteria_ee_MC, Eff_LeptonTriggers_SelectionCriteria_mumu_MC, Eff_LeptonTriggers_SelectionCriteria_emu_MC;
float Eff_MET_SelectionCriteria_ee_MC, Eff_MET_SelectionCriteria_mumu_MC, Eff_MET_SelectionCriteria_emu_MC;
float Alpha_ee_MC, Alpha_mumu_MC, Alpha_emu_MC;
float Eff_UpperUncert_ee_MC, Eff_UpperUncert_mumu_MC, Eff_UpperUncert_emu_MC;
float Eff_LowerUncert_ee_MC, Eff_LowerUncert_mumu_MC, Eff_LowerUncert_emu_MC;

float N_SelectionCriteria_ee_DATA, N_SelectionCriteria_mumu_DATA, N_SelectionCriteria_emu_DATA;
float N_MET_And_LeptonSelection_ee_DATA, N_MET_And_LeptonSelection_mumu_DATA, N_MET_And_LeptonSelection_emu_DATA;
float N_LeptonTriggersAndSelectionCriteria_ee_DATA, N_LeptonTriggersAndSelectionCriteria_mumu_DATA, N_LeptonTriggersAndSelectionCriteria_emu_DATA;
float N_MET_LeptonTriggers_SelectionCriteria_ee_DATA, N_MET_LeptonTriggers_SelectionCriteria_mumu_DATA, N_MET_LeptonTriggers_SelectionCriteria_emu_DATA;
float Eff_ee_DATA, Eff_mumu_DATA, Eff_emu_DATA;
float Eff_MET_LeptonTriggers_SelectionCriteria_ee_DATA, Eff_MET_LeptonTriggers_SelectionCriteria_mumu_DATA, Eff_MET_LeptonTriggers_SelectionCriteria_emu_DATA;
float Eff_LeptonTriggers_SelectionCriteria_ee_DATA, Eff_LeptonTriggers_SelectionCriteria_mumu_DATA, Eff_LeptonTriggers_SelectionCriteria_emu_DATA;
float Eff_MET_SelectionCriteria_ee_DATA, Eff_MET_SelectionCriteria_mumu_DATA, Eff_MET_SelectionCriteria_emu_DATA;
float Alpha_ee_DATA, Alpha_mumu_DATA, Alpha_emu_DATA;
float Eff_UpperUncert_ee_DATA, Eff_UpperUncert_mumu_DATA, Eff_UpperUncert_emu_DATA;
float Eff_LowerUncert_ee_DATA, Eff_LowerUncert_mumu_DATA, Eff_LowerUncert_emu_DATA;
float SF_ee, SF_mumu, SF_emu;
float SF_Uncert_ee, SF_Uncert_mumu, SF_Uncert_emu;
 

void tZq_NanoAOD_Output(const int& MCInt,  	    const int& ProcessInt,  const int& NPLInt,     const int& SRInt,          const int& SBRInt, 
		 	const int& ZPlusJetsCRInt,  const int& ttbarCRInt,  const int& YearInt,    const int& SystematicInt,  const int& ChannelInt, 
			const int& DoubleCountCheckInt){

  std::cout << "YearInt = " << YearInt << std::endl; 
  
  ROOT::RDF::RResultPtr<TH2D> h_bjet_num;
  ROOT::RDF::RResultPtr<TH2D> h_bjet_denom;
  ROOT::RDF::RResultPtr<TH2D> h_nonbjet_num;
  ROOT::RDF::RResultPtr<TH2D> h_nonbjet_denom;
  ROOT::RDF::RResultPtr<TH2D> h_charm_num;
  ROOT::RDF::RResultPtr<TH2D> h_charm_denom;
  ROOT::RDF::RResultPtr<TH2D> h_lightjets_num;
  ROOT::RDF::RResultPtr<TH2D> h_lightjets_denom;
  ROOT::RDF::RResultPtr<TH2D> h_gluon_num;
  ROOT::RDF::RResultPtr<TH2D> h_gluon_denom;

  std::vector<std::string> input_files;
 
  std::string SJER;
  std::string SIGMAJER;

  std::string FileNameJetSmear;
  std::string NPLNumbersString; 

  std::string SampleType;
  std::string Channel;
  std::string Process;
  std::string NonPromptLepton;
  std::string SignalRegion;
  std::string SideBandRegion; 
  std::string ZPlusJetsControlRegion;
  std::string ttbarControlRegion;
  std::string Year;
  std::string Systematic;

  std::string JetMassInput;
  std::string JetPtInput;
  std::string JetEtaInput; 
  std::string JetPhiInput;
  
  std::string PSWeightString;
  std::string DoubleCountString;

  TFile* EGammaEff_inputfile_2016 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2016/egammaEffi_Tight_80X.txt_EGM2D.root", "READ");
  TFile* EGammaEffSys_inputfile_2016 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2016/egammaEffi_Tight_80X.txt_EGM2D.root", "READ");
  TFile* EGammaEffReco_inputfile_2016 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2016/egammaRecoEffi.txt_EGM2D.root", "READ");
  TFile* EGammaEffRecoSys_inputfile_2016 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2016/egammaRecoEffi.txt_EGM2D.root", "READ");

  TFile* EGammaEffReco_HigherPt_inputfile_2017 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2017/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root", "READ");
  TFile* EGammaEffRecoSys_HigherPt_inputfile_2017 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2017/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root", "READ"); 
  TFile* EGammaEffReco_LowPt_inputfile_2017 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2017/egammaEffi.txt_EGM2D_runBCDEF_passingRECO_lowEt.root", "READ");
  TFile* EGammaEffRecoSys_LowPt_inputfile_2017 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2017/egammaEffi.txt_EGM2D_runBCDEF_passingRECO_lowEt.root", "READ");
  TFile* EGammaEff_inputfile_2017 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2017/egammaEffi.txt_EGM2D_runBCDEF_passingTight94X.root", "READ");
  TFile* EGammaEffSys_inputfile_2017 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2017/egammaEffi.txt_EGM2D_runBCDEF_passingTight94X.root", "READ");

  TFile* EGammaEff_inputfile_2018 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2018/2018_ElectronTight.root", "READ"); 
  TFile* EGammaEffSys_inputfile_2018 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2018/2018_ElectronTight.root", "READ");
  TFile* EGammaEffReco_inputfile_2018 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2018/egammaEffi.txt_EGM2D_updatedAll.root", "READ");
  TFile* EGammaEffRecoSys_inputfile_2018 = new TFile("./ScaleFactors/LeptonEnergyCorrections/ElectronSFs/2018/egammaEffi.txt_EGM2D_updatedAll.root", "READ");

  TH2* EGammaEff2016_histo = dynamic_cast<TH2*>(EGammaEff_inputfile_2016->Get("EGamma_SF2D")->Clone());
  EGammaEff2016_histo->SetDirectory(nullptr);
  TH2* EGammaEffSys2016_histo = dynamic_cast<TH2*>(EGammaEffSys_inputfile_2016->Get("EGamma_SF2D")->Clone());
  EGammaEffSys2016_histo->SetDirectory(nullptr);
  TH2* EGammaEffReco2016_histo = dynamic_cast<TH2*>(EGammaEffReco_inputfile_2016->Get("EGamma_SF2D")->Clone());
  EGammaEffReco2016_histo->SetDirectory(nullptr);
  TH2* EGammaEffRecoSys2016_histo = dynamic_cast<TH2*>(EGammaEffRecoSys_inputfile_2016->Get("EGamma_SF2D")->Clone());
  EGammaEffRecoSys2016_histo->SetDirectory(nullptr);

  TH2* EGammaEff2017_histo = dynamic_cast<TH2*>(EGammaEff_inputfile_2017->Get("EGamma_SF2D")->Clone());
  EGammaEff2017_histo->SetDirectory(nullptr);
  TH2* EGammaEffSys2017_histo = dynamic_cast<TH2*>(EGammaEffSys_inputfile_2017->Get("EGamma_SF2D")->Clone());
  EGammaEffSys2017_histo->SetDirectory(nullptr);
  TH2* EGammaEffReco_LowPt_2017_histo = dynamic_cast<TH2*>(EGammaEffReco_LowPt_inputfile_2017->Get("EGamma_SF2D")->Clone());
  EGammaEffReco_LowPt_2017_histo->SetDirectory(nullptr); 
  TH2* EGammaEffRecoSys_LowPt_2017_histo = dynamic_cast<TH2*>(EGammaEffRecoSys_LowPt_inputfile_2017->Get("EGamma_SF2D")->Clone());
  EGammaEffRecoSys_LowPt_2017_histo->SetDirectory(nullptr);
  TH2* EGammaEffReco_HigherPt_2017_histo = dynamic_cast<TH2*>(EGammaEffReco_HigherPt_inputfile_2017->Get("EGamma_SF2D")->Clone());
  EGammaEffReco_HigherPt_2017_histo->SetDirectory(nullptr);
  TH2* EGammaEffRecoSys_HigherPt_2017_histo = dynamic_cast<TH2*>(EGammaEffRecoSys_HigherPt_inputfile_2017->Get("EGamma_SF2D")->Clone());
  EGammaEffRecoSys_HigherPt_2017_histo->SetDirectory(nullptr);

  TH2* EGammaEff2018_histo = dynamic_cast<TH2*>(EGammaEff_inputfile_2018->Get("EGamma_SF2D")->Clone());
  EGammaEff2018_histo->SetDirectory(nullptr);
  TH2* EGammaEffSys2018_histo = dynamic_cast<TH2*>(EGammaEffSys_inputfile_2018->Get("EGamma_SF2D")->Clone());
  EGammaEffSys2018_histo->SetDirectory(nullptr);
  TH2* EGammaEffReco2018_histo = dynamic_cast<TH2*>(EGammaEffReco_inputfile_2018->Get("EGamma_SF2D")->Clone());
  EGammaEffReco2018_histo->SetDirectory(nullptr);
  TH2* EGammaEffRecoSys2018_histo = dynamic_cast<TH2*>(EGammaEffRecoSys_inputfile_2018->Get("EGamma_SF2D")->Clone());
  EGammaEffRecoSys2018_histo->SetDirectory(nullptr);

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


  TFile * inputfile_RunsBCDEF_ID_2016 = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2016/MuonID_EfficienciesAndSF_BCDEF.root", "READ");
  TH2* histo_RunsBCDEF_ID_2016 = dynamic_cast<TH2*>(inputfile_RunsBCDEF_ID_2016->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio")->Clone());
  histo_RunsBCDEF_ID_2016->SetDirectory(nullptr);

  TFile* inputfile_RunsGH_ID_2016 = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2016/MuonID_EfficienciesAndSF_GH.root", "READ");
  TH2* histo_RunsGH_ID_2016 = dynamic_cast<TH2*>(inputfile_RunsGH_ID_2016->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/pt_abseta_ratio")->Clone());
  histo_RunsGH_ID_2016->SetDirectory(nullptr);

  TFile* inputfile_RunsBCDEF_ISO_2016 = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2016/MuonISO_EfficienciesAndSF_BCDEF.root", "READ");
  TH2* histo_RunsBCDEF_ISO_2016 = dynamic_cast<TH2*>(inputfile_RunsBCDEF_ISO_2016->Get("TightISO_TightID_pt_eta/pt_abseta_ratio")->Clone());
  histo_RunsBCDEF_ISO_2016->SetDirectory(nullptr);

  TFile* inputfile_RunsGH_ISO_2016 = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2016/MuonISO_EfficienciesAndSF_GH.root", "READ");
  TH2* histo_RunsGH_ISO_2016 = dynamic_cast<TH2*>(inputfile_RunsGH_ISO_2016->Get("TightISO_TightID_pt_eta/pt_abseta_ratio")->Clone());
  histo_RunsGH_ISO_2016->SetDirectory(nullptr);

  TFile* inputfile_RunsBCDEF_ID_2017 = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ID.root", "READ");
  TH2* histo_RunsBCDEF_ID_2017 = dynamic_cast<TH2*>(inputfile_RunsBCDEF_ID_2017->Get("NUM_TightID_DEN_genTracks_pt_abseta")->Clone());
  histo_RunsBCDEF_ID_2017->SetDirectory(nullptr);

  TFile* inputfile_RunsBCDEF_ID_Sys_2017 = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ID_syst.root", "READ");
  TH2* histo_RunsBCDEF_ID_Sys_2017 = dynamic_cast<TH2*>(inputfile_RunsBCDEF_ID_Sys_2017->Get("NUM_TightID_DEN_genTracks_pt_abseta")->Clone());
  histo_RunsBCDEF_ID_Sys_2017->SetDirectory(nullptr);

  TFile* inputfile_RunsBCDEF_ID_Sys_Stat_2017 = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ID_syst.root", "READ");
  TH2* histo_RunsBCDEF_ID_Sys_Stat_2017 = dynamic_cast<TH2*>(inputfile_RunsBCDEF_ID_Sys_Stat_2017->Get("NUM_TightID_DEN_genTracks_pt_abseta_stat")->Clone());
  histo_RunsBCDEF_ID_Sys_Stat_2017->SetDirectory(nullptr);

  TFile* inputfile_RunsBCDEF_ID_Sys_Syst_2017 = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ID_syst.root", "READ");
  TH2* histo_RunsBCDEF_ID_Sys_Syst_2017 = dynamic_cast<TH2*>(inputfile_RunsBCDEF_ID_Sys_Syst_2017->Get("NUM_TightID_DEN_genTracks_pt_abseta_syst")->Clone());
  histo_RunsBCDEF_ID_Sys_Syst_2017->SetDirectory(nullptr);

  TFile* inputfile_RunsBCDEF_ISO_2017 = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ISO.root", "READ");
  TH2* histo_RunsBCDEF_ISO_2017 = dynamic_cast<TH2*>(inputfile_RunsBCDEF_ISO_2017->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta")->Clone());
  histo_RunsBCDEF_ISO_2017->SetDirectory(nullptr);


  TFile* inputfile_RunsBCDEF_ISO_Sys_2017 = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ISO_syst.root", "READ");
  TH2* histo_RunsBCDEF_ISO_Sys_2017 = dynamic_cast<TH2*>(inputfile_RunsBCDEF_ISO_Sys_2017->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta")->Clone());
  histo_RunsBCDEF_ISO_Sys_2017->SetDirectory(nullptr);

  TFile* inputfile_RunsBCDEF_ISO_Sys_Stat_2017 = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ISO_syst.root", "READ");
  TH2* histo_RunsBCDEF_ISO_Sys_Stat_2017 = dynamic_cast<TH2*>(inputfile_RunsBCDEF_ISO_Sys_Stat_2017->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta_stat")->Clone());
  histo_RunsBCDEF_ISO_Sys_Stat_2017->SetDirectory(nullptr);

  TFile* inputfile_RunsBCDEF_ISO_Sys_Syst_2017 = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2017/Muon_RunBCDEF_SF_ISO_syst.root", "READ");
  TH2* histo_RunsBCDEF_ISO_Sys_Syst_2017 = dynamic_cast<TH2*>(inputfile_RunsBCDEF_ISO_Sys_Syst_2017->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta_syst")->Clone());
  histo_RunsBCDEF_ISO_Sys_Syst_2017->SetDirectory(nullptr);

  TFile* inputfile_RunsABCD_ID_2018 = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2018/RunABCD_SF_ID.root", "READ"); //need to double check if root file is correct
  TH2* histo_RunsABCD_ID_2018 = dynamic_cast<TH2*>(inputfile_RunsABCD_ID_2018->Get("NUM_TightID_DEN_TrackerMuons_pt_abseta")->Clone());
  histo_RunsABCD_ID_2018->SetDirectory(nullptr);

  TFile* inputfile_RunsABCD_ISO_2018 = new TFile("./ScaleFactors/LeptonEfficiency/MuonSFs/2018/RunABCD_SF_ISO.root", "READ"); //need to double check if root file is correct
  TH2* histo_RunsABCD_ISO_2018 = dynamic_cast<TH2*>(inputfile_RunsABCD_ISO_2018->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta")->Clone());
  histo_RunsABCD_ISO_2018->SetDirectory(nullptr);

  
  inputfile_RunsBCDEF_ID_2016->Close(); 
  inputfile_RunsGH_ID_2016->Close();
  inputfile_RunsBCDEF_ISO_2016->Close();
  inputfile_RunsGH_ISO_2016->Close();
  inputfile_RunsBCDEF_ID_2017->Close();
  inputfile_RunsBCDEF_ID_Sys_2017->Close();
  inputfile_RunsBCDEF_ID_Sys_Stat_2017->Close();
  inputfile_RunsBCDEF_ID_Sys_Syst_2017->Close();
  inputfile_RunsBCDEF_ISO_2017->Close();
  inputfile_RunsBCDEF_ISO_Sys_2017->Close();
  inputfile_RunsBCDEF_ISO_Sys_Stat_2017->Close();
  inputfile_RunsBCDEF_ISO_Sys_Syst_2017->Close();
  inputfile_RunsABCD_ID_2018->Close();
  inputfile_RunsABCD_ISO_2018->Close();

  //Pile up modelling
  //2016
  TFile *dataPileupFile_2016 = new TFile("./ScaleFactors/PileUp/2016/truePileupTest.root", "READ");
  TH1D *dataPU_2016 = dynamic_cast<TH1D*>(dataPileupFile_2016->Get("pileup")->Clone());

  TFile *mcPileupFile_2016 = new TFile("./ScaleFactors/PileUp/2016/pileupMC.root", "READ");
  TH1D* mcPU_2016 = dynamic_cast<TH1D*>(mcPileupFile_2016->Get("pileup")->Clone());

  //2016 part 1
  TFile *dataPileupFile_2016_part1 = new TFile("./ScaleFactors/PileUp/2016/truePileupTest_part1.root", "READ");
  TH1D *dataPU_2016_part1 = dynamic_cast<TH1D*>(dataPileupFile_2016_part1->Get("pileup")->Clone());
  TFile *mcPileupFile_2016_part1 = new TFile("./ScaleFactors/PileUp/2016/pileupMC.root", "READ");
  TH1D* mcPU_2016_part1 = dynamic_cast<TH1D*>(mcPileupFile_2016_part1->Get("pileup")->Clone());

  //2016 part 2
  TFile *dataPileupFile_2016_part2 = new TFile("./ScaleFactors/PileUp/2016/truePileupTest_part2.root", "READ");
  TH1D *dataPU_2016_part2 = dynamic_cast<TH1D*>(dataPileupFile_2016_part2->Get("pileup")->Clone());
  TFile *mcPileupFile_2016_part2 = new TFile("./ScaleFactors/PileUp/2016/pileupMC.root", "READ");
  TH1D* mcPU_2016_part2 = dynamic_cast<TH1D*>(mcPileupFile_2016_part2->Get("pileup")->Clone());

  //2017
  TFile *dataPileupFile_2017 = new TFile("./ScaleFactors/PileUp/2017/truePileupTest.root", "READ");
  TH1D *dataPU_2017 = dynamic_cast<TH1D*>(dataPileupFile_2017->Get("pileup")->Clone());
  TFile *mcPileupFile_2017 = new TFile("./ScaleFactors/PileUp/2017/pileupMC.root", "READ");
  TH1D* mcPU_2017 = dynamic_cast<TH1D*>(mcPileupFile_2017->Get("pileup")->Clone());

  //2018
  TFile *dataPileupFile_2018 = new TFile("./ScaleFactors/PileUp/2018/MyDataPileupHistogram2018.root", "READ");
  TH1D *dataPU_2018 = dynamic_cast<TH1D*>(dataPileupFile_2018->Get("pileup")->Clone());
  TFile *mcPileupFile_2018 = new TFile("./ScaleFactors/PileUp/2018/pileupMC2018.root", "READ");
  TH1D* mcPU_2018 = dynamic_cast<TH1D*>(mcPileupFile_2018->Get("pileup")->Clone());

  //Systematic files
  //2016
  TFile *systUpFile_2016 = new TFile("./ScaleFactors/PileUp/2016/truePileupUp.root", "READ");
  TH1D *pileupUpHist_2016 = dynamic_cast<TH1D*>(systUpFile_2016->Get("pileup")->Clone());
  TFile *systDownFile_2016 = new TFile("./ScaleFactors/PileUp/2016/truePileupDown.root", "READ");
  TH1D *pileupDownHist_2016 = dynamic_cast<TH1D*>(systDownFile_2016->Get("pileup")->Clone());

  //part 1
  TFile *systUpFile_2016_part1 = new TFile("./ScaleFactors/PileUp/2016/truePileupUp_part1.root", "READ");
  TH1D *pileupUpHist_2016_part1 = dynamic_cast<TH1D*>(systUpFile_2016_part1->Get("pileup")->Clone());
  TFile *systDownFile_2016_part1 = new TFile("./ScaleFactors/PileUp/2016/truePileupDown_part1.root", "READ");
  TH1D *pileupDownHist_2016_part1 = dynamic_cast<TH1D*>(systDownFile_2016_part1->Get("pileup")->Clone());

  //part 2
  TFile *systUpFile_2016_part2 = new TFile("./ScaleFactors/PileUp/2016/truePileupUp_part2.root", "READ");
  TH1D *pileupUpHist_2016_part2 = dynamic_cast<TH1D*>(systUpFile_2016_part2->Get("pileup")->Clone());
  TFile *systDownFile_2016_part2 = new TFile("./ScaleFactors/PileUp/2016/truePileupDown_part2.root", "READ");
  TH1D *pileupDownHist_2016_part2 = dynamic_cast<TH1D*>(systDownFile_2016_part2->Get("pileup")->Clone());

  TH1D *puReweight_2016 = dynamic_cast<TH1D*>(dataPU_2016->Clone());
  puReweight_2016->Scale(1.0 / puReweight_2016->Integral());
  mcPU_2016->Scale(1.0 / mcPU_2016->Integral());
  puReweight_2016->Divide(mcPU_2016);
  puReweight_2016->SetDirectory(nullptr);

  TH1D *puReweight_2016_part1 = dynamic_cast<TH1D*>(dataPU_2016_part1->Clone());
  puReweight_2016_part1->Scale(1.0 / puReweight_2016_part1->Integral());
  mcPU_2016_part1->Scale(1.0 / mcPU_2016_part1->Integral());
  puReweight_2016_part1->Divide(mcPU_2016_part1);
  puReweight_2016_part1->SetDirectory(nullptr);

  TH1D *puReweight_2016_part2 = dynamic_cast<TH1D*>(dataPU_2016_part2->Clone());
  puReweight_2016_part2->Scale(1.0 / puReweight_2016_part2->Integral());
  mcPU_2016_part2->Scale(1.0 / mcPU_2016_part2->Integral());
  puReweight_2016_part2->Divide(mcPU_2016_part2);
  puReweight_2016_part2->SetDirectory(nullptr);

  // 2017
  TFile *systUpFile_2017 = new TFile("./ScaleFactors/PileUp/2017/truePileupUp.root", "READ");
  TH1D *pileupUpHist_2017 = dynamic_cast<TH1D*>(systUpFile_2017->Get("pileup")->Clone());
  TFile *systDownFile_2017 = new TFile("./ScaleFactors/PileUp/2017/truePileupDown.root", "READ");
  TH1D *pileupDownHist_2017 = dynamic_cast<TH1D*>(systDownFile_2017->Get("pileup")->Clone());

  TH1D *puReweight_2017 = dynamic_cast<TH1D*>(dataPU_2017->Clone());
  puReweight_2017->Scale(1.0 / puReweight_2017->Integral());
  mcPU_2017->Scale(1.0 / mcPU_2017->Integral());
  puReweight_2017->Divide(mcPU_2017);
  puReweight_2017->SetDirectory(nullptr);

  //2018
  TFile *systUpFile_2018 = new TFile("./ScaleFactors/PileUp/2018/MyDataPileupHistogramScaleUp2018.root", "READ");
  TH1D *pileupUpHist_2018 = dynamic_cast<TH1D*>(systUpFile_2018->Get("pileup")->Clone());
  TFile *systDownFile_2018 = new TFile("./ScaleFactors/PileUp/2018/MyDataPileupHistogramScaleDown2018.root", "READ");
  TH1D *pileupDownHist_2018 = dynamic_cast<TH1D*>(systDownFile_2018->Get("pileup")->Clone());

  TH1D *puReweight_2018 = dynamic_cast<TH1D*>(dataPU_2018->Clone());
  puReweight_2018->Scale(1.0 / puReweight_2018->Integral());
  mcPU_2018->Scale(1.0 / mcPU_2018->Integral());
  puReweight_2018->Divide(mcPU_2018);
  puReweight_2018->SetDirectory(nullptr);

  ///Systematic sample
  //2016
  TH1D *puSystUp_2016 = dynamic_cast<TH1D*>(pileupUpHist_2016->Clone());
  puSystUp_2016->Scale(1.0 / puSystUp_2016->Integral());
  puSystUp_2016->Divide(mcPU_2016);
  puSystUp_2016->SetDirectory(nullptr);
  TH1D *puSystDown_2016 = dynamic_cast<TH1D*>(pileupDownHist_2016->Clone());
  puSystDown_2016->Scale(1.0 / puSystDown_2016->Integral());
  puSystDown_2016->Divide(mcPU_2016);
  puSystDown_2016->SetDirectory(nullptr);

  //2016 part 1
  TH1D *puSystUp_2016_part1 = dynamic_cast<TH1D*>(pileupUpHist_2016_part1->Clone());
  puSystUp_2016_part1->Scale(1.0 / puSystUp_2016_part1->Integral());
  puSystUp_2016_part1->Divide(mcPU_2016_part1);
  puSystUp_2016_part1->SetDirectory(nullptr);
  TH1D *puSystDown_2016_part1 = dynamic_cast<TH1D*>(pileupDownHist_2016_part1->Clone());
  puSystDown_2016_part1->Scale(1.0 / puSystDown_2016_part1->Integral());
  puSystDown_2016_part1->Divide(mcPU_2016_part1);
  puSystDown_2016_part1->SetDirectory(nullptr);

  //2016 part 2
  TH1D *puSystUp_2016_part2 = dynamic_cast<TH1D*>(pileupUpHist_2016_part2->Clone());
  puSystUp_2016_part2->Scale(1.0 / puSystUp_2016_part2->Integral());
  puSystUp_2016_part2->Divide(mcPU_2016_part2);
  puSystUp_2016_part2->SetDirectory(nullptr);
  TH1D *puSystDown_2016_part2 = dynamic_cast<TH1D*>(pileupDownHist_2016_part2->Clone());
  puSystDown_2016_part2->Scale(1.0 / puSystDown_2016_part2->Integral());
  puSystDown_2016_part2->Divide(mcPU_2016_part2);
  puSystDown_2016_part2->SetDirectory(nullptr);

  //2017
  TH1D *puSystUp_2017 = dynamic_cast<TH1D*>(pileupUpHist_2017->Clone());
  puSystUp_2017->Scale(1.0 / puSystUp_2017->Integral());
  puSystUp_2017->Divide(mcPU_2017);
  puSystUp_2017->SetDirectory(nullptr);
  TH1D *puSystDown_2017 = dynamic_cast<TH1D*>(pileupDownHist_2017->Clone());
  puSystDown_2017->Scale(1.0 / puSystDown_2017->Integral());
  puSystDown_2017->Divide(mcPU_2017);
  puSystDown_2017->SetDirectory(nullptr);

  //2018
  TH1D *puSystUp_2018 = dynamic_cast<TH1D*>(pileupUpHist_2018->Clone());
  puSystUp_2018->Scale(1.0 / puSystUp_2018->Integral());
  puSystUp_2018->Divide(mcPU_2018);
  puSystUp_2018->SetDirectory(nullptr);
  TH1D *puSystDown_2018 = dynamic_cast<TH1D*>(pileupDownHist_2018->Clone());
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



  //Setting the SampleType, Channel, Process, NPL, SR, SBR, ZPlusJetsCR, ttbarCR, Year and Systematic strings for the output file names
  switch(MCInt){

	case 0: SampleType = "data";
		JetMassInput = "Jet_mass"; JetPtInput = "Jet_pt"; JetEtaInput = "Jet_eta"; JetPhiInput = "Jet_phi";
		break;

	case 1: SampleType = "MC";
		JetMassInput = "SmearedJetMass"; JetPtInput = "SmearedJetPt"; JetEtaInput = "SmearedJetEta"; JetPhiInput = "SmearedJetPhi";
		break;

  }

  switch(ProcessInt){

	case 0: Process = "tZq"; 

		switch(YearInt){
			case 2016: input_files = {"/data/disk2/nanoAOD_2016/tZq_ll/*"}; break;
			case 2017: input_files = {"/data/disk0/nanoAOD_2017/tZq_ll/*"}; break;
			case 2018: input_files = {"/data/disk1/nanoAOD_2018/tZq_ll/*"}; break;
			default: std::cout << "Inside the tZq switch statement. Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
		}

		break;

	case 1: Process = "ZPlusJets_M50_aMCatNLO";
		
                switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_M50_aMCatNLO/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/DYJetsToLL_NanoAODv5/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/DYJetsToLL_NanoAODv5/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                }
                
                break;

	case 2: Process = "ZPlusJets_M50_aMCatNLO_ext";
	
		switch(YearInt){
                        case 2016: break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/DYJetsToLL_ext_NanoAODv5/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/DYJetsToLL_ext_NanoAODv5/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                }

                break;

	case 3: Process = "ZPlusJets_M50_Madgraph";

		switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_M50_Madgraph/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                }

                break;

	case 4: Process = "ZPlusJets_M50_Madgraph_ext";

		switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_M50_Madgraph_ext/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                }

                break;

	case 5: Process = "ZPlusJets_M10To50_aMCatNLO";

		switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_M10ToM50_aMCatNLO/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                }

                break;

        case 6: Process = "ZPlusJets_M10To50_aMCatNLO_ext";

		switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_M10ToM50_ext_aMCatNLO/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                }

                break;

	case 7: Process = "ZPlusJets_M10To50_Madgraph";

		switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_M10To50_Madgraph/*"}; break;
                        case 2017: break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/DYJetsToLL_M10to50/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                }

                break;

	case 8: Process = "ZPlusJets_PtBinned_0To50";

		switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_0To50/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                }

		break;

	case 9: Process = "ZPlusJets_PtBinned_50To100";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_50To100/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

	case 10: Process = "ZPlusJets_PtBinned_50To100_ext";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_50To100_ext/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

	case 11: Process = "ZPlusJets_PtBinned_100To250";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_100To250/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

	case 12: Process = "ZPlusJets_PtBinned_100To250_ext1";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_100To250_ext1/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

	case 13: Process = "ZPlusJets_PtBinned_100To250_ext2";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_100To250_ext2/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

	case 14: Process = "ZPlusJets_PtBinned_100To250_ext5";
		
		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_100To250_ext5/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 15: Process = "ZPlusJets_PtBinned_250To400";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_250To400/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 16: Process = "ZPlusJets_PtBinned_250To400_ext1";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_250To400_ext1/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 17: Process = "ZPlusJets_PtBinned_250To400_ext2";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_250To400_ext2/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 18: Process = "ZPlusJets_PtBinned_250To400_ext5";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_250To400_ext5/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 19: Process = "ZPlusJets_PtBinned_400To650";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_400To650/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

	case 20: Process = "ZPlusJets_PtBinned_400To650_ext1";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_400To650_ext1/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 21: Process = "ZPlusJets_PtBinned_400To650_ext2";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_400To650_ext2/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 22: Process = "ZPlusJets_PtBinned_650ToInf";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_650ToInf/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 23: Process = "ZPlusJets_PtBinned_650ToInf_ext1";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_650ToInf_ext1/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 25: Process = "ZPlusJets_PtBinned_650ToInf_ext2";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZPlusJets_PtBinned_650ToInf_ext2/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 26: Process = "ttbar_2l2nu";

		 switch(YearInt){
                        case 2016: break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/ttbar_2l2nu_NanoAODv5/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/ttbar_2l2nu/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

	case 27: Process = "ttbar_madgraph";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ttbar_madgraph/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/TTToHadronic/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/ttbar_madgraph/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 28: Process = "ttbar_madgraph_ext";

		 switch(YearInt){
                        case 2016: break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 29: Process = "ttbar_TTToHadronic";
	
		 switch(YearInt){
                        case 2016: break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/TTToHadronic/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/TTToHadronic/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 30: Process = "ttbar_TTToSemileptonic";

		 switch(YearInt){
                        case 2016: break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/TTToSemileptonic/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/TTToSemileptonic/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 31: Process = "ttbar_atMCaNLO";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ttbar_aMCatNLO/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/ttbar_aMCatNLO/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/ttbar_aMCatNLO/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 32: Process = "ttbar_inc";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ttbar_inc/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

	case 33: Process = "SingleTop_tchannel_top";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ST_tchannel_top/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/ST_tchannel_top/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/ST_tchannel_top/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 34: Process = "SingleTop_tchannel_top_ScaleUp";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ST_tchannel_top_ScaleUp/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 35: Process = "SingleTop_tchannel_top_ScaleDown";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ST_tchannel_top_ScaleDown_NanoAODv6/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 36: Process = "SingleTop_tchannel_antitop";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ST_tchannel_antitop/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/ST_tchannel_tbar/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/ST_tchannel_tbar/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 37: Process = "SingleTop_tchannel_antitop_ScaleUp";

		 switch(YearInt){
                        case 2016: break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 38: Process = "SingleTop_tchannel_antitop_ScaleDown";

		 switch(YearInt){
                        case 2016: break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

	case 39: Process = "SingleTop_schannel";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ST_schannel/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/ST_schannel/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/ST_schannel/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 40: Process = "ttbar_hdampUP";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/TT_hdampUP/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 41: Process = "ttbar_hdampUP_ext";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/TT_hdampUP_ext/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 42: Process = "ttbar_hdampDOWN";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/TT_hdampDOWN/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 43: Process = "ttbar_hdampDOWN_ext";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/TT_hdampDOWN_ext/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 44: Process = "SingleTop_tchannel_top_hdampUP";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ST_tchannel_top_hdampup/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

	case 45: Process = "SingleTop_tchannel_top_hdampDOWN";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ST_tchannel_top_hdampdown/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 46: Process = "ttbar_isr_UP";

		 switch(YearInt){
                        case 2016: input_files =  {"/data/disk2/nanoAOD_2016/TT_isr_UP/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 47: Process = "ttbar_isr_DOWN";
		
		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/TT_isr_DOWN/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 48: Process = "ttbar_isr_DOWN_ext";
	
		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/TT_isr_DOWN_ext/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 49: Process = "ttbar_fsr_UP";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/TT_fsr_UP/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 50: Process = "ttbar_fsr_UP_ext";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/TT_fsr_UP_ext/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

	case 51: Process = "ttbar_fsr_DOWN";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/TT_fsr_DOWN/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 52: Process = "ttbar_fsr_DOWN_ext";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/TT_fsr_DOWN_ext/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 53: Process = "SingleTop_tW";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/tW_top/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/ST_tW/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/ST_tW/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 54: Process = "SingleTop_tW_ScaleUp";
	
		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/tW_top_ScaleUp/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 55: Process = "SingleTop_tW_ScaleDown";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/tW_top_ScaleDown/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 56: Process = "SingleTop_tbarW";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/tW_tbar/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/ST_tbarW/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/ST_tbarW/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

	case 57: Process = "SingleTop_tbarW_ScaleUp";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/tW_tbar_ScaleUp/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 58: Process = "SingleTop_tbarW_ScaleDown";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/tW_tbar_ScaleDown/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 59: Process = "SingleTop_tHq";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/tHq/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/tHq/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/tHq/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 60: Process = "SingleTop_tZq_W_lept_Z_had";

		 switch(YearInt){
                        case 2016: break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/tZq_W_lept_Z_had/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/tZq_W_lept_Z_had/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }
	
                 break;

        case 61: Process = "SingleTop_tWZ_tWll";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/tWZ_tWLL/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/tWZ_tWLL_NanoAODv5/*"}; break;
                        case 2018: input_files =  {"/data/disk1/nanoAOD_2018/tWZ_tWLL_NanoAODv5/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 62: Process = "VV_ZZTo2l2nu";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZZTo2L2Nu/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/ZZTo2L2Nu/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/ZZTo2L2Nu/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

	case 63: Process = "VV_ZZTo2l2nu_ext1";

		 switch(YearInt){
                        case 2016: break;
                        case 2017: break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/ZZTo2L2Nu_ext/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

	case 64: Process = "VV_ZZTo2l2nu_ext2";

                 switch(YearInt){
                        case 2016: break;
                        case 2017: break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/ZZTo2L2Nu_ext2/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 65: Process = "VV_ZZTo2l2Q";
	
		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZZTo2L2Q/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/ZZTo2L2Q/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/ZZTo2L2Q/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 66: Process = "VV_ZZTo4L";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZZTo4L/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/ZZTo4L/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/ZZTo4L/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 67: Process = "VV_WW1nuqq";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/WWTo1L1Nu2Q/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/WWTolnuqq/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/WWToLNuQQ/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 68: Process = "VV_WZTo2l2Q";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/WZTo2L2Q/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/WZTo2L2Q/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/WZTo2L2Q/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 69: Process = "VV_WZTo3lNu";

		 switch(YearInt){
                        case 2016: break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/WZTo3LNu/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/WZTo3LNu/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 70: Process = "VV_WZTo1l1Nu2Q";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/WZTo1L1Nu2Q/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/WZTo1L1Nu2Q/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/WZTo1L1Nu2Q/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 71: Process = "VV_WWTo2l2Nu";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/WWTo2L2Nu/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/WWTo2L2Nu/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/WWTo2L2Nu/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 72: Process = "VV_WWToLNuQQ";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/WWToLNuQQ/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/WWToLNuQQ/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/WWToLNuQQ/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 73: Process = "VV_WWToLNuQQ_ext";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/WWToQQ_ext/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 74: Process = "VV_WGToLNuG";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/WGToLNuG/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/WGToLNuG/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/WGToLNuG/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 75: Process = "VV_ZGToLLG";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZGToLLG/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/ZGToLLG/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/ZGToLLG/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

	case 76: Process = "VVV_WWWTo4F";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/WWWTo4F/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/WWWTo4F/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/WWWTo4F/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 77: Process = "VVV_WWZTo4F";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/WWZ/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/WWZTo4F/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/WWZTo4F/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 78: Process = "VVV_WZZ";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/WZZ/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/WZZ/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/WZZ/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 79: Process = "VVV_ZZZ";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ZZZ/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/ZZZ/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/ZZZ/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

	case 80: Process = "WPlusJets";
		
		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/WJetsToLNu/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/WJetsToLNu/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/WJetsToLNu/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 81: Process = "WPlusJets_ext";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/WJetsToLNu_ext2/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 82: Process = "ttbarV_ttWJetsToLNu";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ttWJetsToLNu/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/ttWJetsToLNu/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/ttWJetsToLNu/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 83: Process = "ttbarV_ttWJetsToLNu_ext";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ttWJetsToLNu/ttWJetsToLNu_ext2/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

	case 84: Process = "ttbarV_ttWJetsToQQ";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ttWJetsToQQ/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/ttWJetsToQQ/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/ttWJetsToQQ/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 85: Process = "ttbarV_ttZToLL";

		 switch(YearInt){
                        case 2016: break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/ttZToLL/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/ttZToLL/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

	case 86: Process = "ttbarV_ttZToLLNuNu";

                 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ttZToLLNuNu/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/ttZToLLNuNu/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/ttZToLLNuNu/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

	case 87: Process = "ttbarV_ttZToLLNuNu_ext2";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ttZToLLNuNu_ext2/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 88: Process = "ttbarV_ttZToLLNuNu_ext3";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ttZToLLNuNu_ext3/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;
	
	case 89: Process = "ttbarV_ttZToQQ";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ttZToQQ/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/ttZToQQ/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/ttZToQQ/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 90: Process = "ttbarV_ttZToQQ_ext";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ttZToQQ_ext/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/ttZToQQ_ext/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/ttZToQQ_ext/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

        case 91: Process = "ttbarV_ttgamma";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ttgamma/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/ttgamma/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/ttgamma/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;
   
        case 92: Process = "ttbarV_ttgamma_ext";

		 switch(YearInt){
                        case 2016: break;
                        case 2017: break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/ttgamma_ext/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

	case 93: Process = "ttbarV_ttHTobb";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ttHTobb/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/ttHTobb/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/ttHTobb/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

	case 94: Process = "ttbarV_ttHToNonbb";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ttHToNonbb/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/ttHToNonbb/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/ttHToNonbb/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

	case 95: Process = "TriggerSF_MC";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/ttbar_inc/*.root"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/ttbar_2l2nu/*.root"}; break;
                        case 2018: input_files = {"/data/disk0/nanoAOD_2017/ttbar_2l2nu/*.root"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

	case 96: Process = "TriggerSF_DATA";

                 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/METRun2016B/*.root", "/data/disk2/nanoAOD_2016/METRun2016C/*.root", "/data/disk2/nanoAOD_2016/METRun2016D/*.root", "/data/disk2/nanoAOD_2016/METRun2016E/*.root", "/data/disk2/nanoAOD_2016/METRun2016F/*.root", "/data/disk2/nanoAOD_2016/METRun2016G/*.root", "/data/disk2/nanoAOD_2016/METRun2016H/*.root"}; break;

                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/METRun2017B/*.root", "/data/disk0/nanoAOD_2017/METRun2017C/*.root", "/data/disk0/nanoAOD_2017/METRun2017D/*.root", "/data/disk0/nanoAOD_2017/METRun2017E/*.root", "/data/disk0/nanoAOD_2017/METRun2017F/*"}; break;

                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/METRun2018B/*.root", "/data/disk1/nanoAOD_2018/METRun2018C/*.root", "/data/disk1/nanoAOD_2018/METRun2018D/*.root"}; break;

                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

                 break;

	case 97: Process = "Data_DoubleEGRunB";

		 switch(YearInt){
			case 2016: input_files = {"/data/disk2/nanoAOD_2016/DoubleEGRun2016B/*"}; break;
			case 2017: input_files = {"/data/disk0/nanoAOD_2017/DoubleEGRun2017B/*"}; break;
			case 2018: break;
			default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
		 }

	case 98: Process = "Data_DoubleEGRunC";

		switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/DoubleEGRun2016C/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/DoubleEGRun2017C/*"}; break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

        case 99: Process = "Data_DoubleEGRunD";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/DoubleEGRun2016D/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/DoubleEGRun2017D/*"}; break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

	case 100: Process = "Data_DoubleEGRunE";

		switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/DoubleEGRun2016E/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/DoubleEGRun2017E/*"}; break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }


	case 101: Process = "Data_DoubleEGRunF";

		 switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/DoubleEGRun2016F/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/DoubleEGRun2017F/*"}; break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

	case 102: Process = "Data_DoubleEGRunG";

		switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/DoubleEGRun2016G/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }


	case 103: Process = "Data_DoubleEGRunH";

		  switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/DoubleEGRun2016H/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                 }

	case 104: Process = "Data_DoubleMuonRunB";

		  switch(YearInt){
			case 2016: input_files = {"/data/disk2/nanoAOD_2016/DoubleMuonRun2016B/*"}; break;
			case 2017: input_files = {"/data/disk0/nanoAOD_2017/DoubleMuonRun2017B/*"}; break;
			case 2018: input_files = {"/data/disk1/nanoAOD_2018/DoubleMuonRun2018B/*"}; break;
			default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
		  }
	case 105: Process = "Data_DoubleMuonRunC";

                  switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/DoubleMuonRun2016C/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/DoubleMuonRun2017C/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/DoubleMuonRun2018C/*"}; break;
			default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 106: Process = "Data_DoubleMuonRunD";

                  switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/DoubleMuonRun2016D/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/DoubleMuonRun2017D/*"}; break;
                        case 2018: input_files = {"/data/disk1/nanoAOD_2018/DoubleMuonRun2018D/*"}; break;
			default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 107: Process = "Data_DoubleMuonRunE";

                  switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/DoubleMuonRun2016E/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/DoubleMuonRun2017E/*"}; break;
                        case 2018: break;
			default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 108: Process = "Data_DoubleMuonRunF";

                  switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/DoubleMuonRun2016F/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/DoubleMuonRun2017F/*"}; break;
                        case 2018: break;
			default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 109: Process = "Data_DoubleMuonRunG";

                  switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/DoubleMuonRun2016G/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 110: Process = "Data_DoubleMuonRunH";

                  switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/DoubleMuonRun2016H/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 111: Process = "Data_MuonEGRunB";

                  switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/MuonEGRun2016B/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/MuonEGRun2017B/*"}; break;
                        case 2018: input_files = {"/data/disk3/nanoAOD_2018/MuonEGRun2018B/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 112: Process = "Data_MuonEGRunC";

                  switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/MuonEGRun2016C/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/MuonEGRun2017C/*"}; break;
                        case 2018: input_files = {"/data/disk3/nanoAOD_2018/MuonEGRun2018C/*"}; break;
			default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }
	
	case 113: Process = "Data_MuonEGRunD";

                  switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/MuonEGRun2016D/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/MuonEGRun2017D/*"}; break;
                        case 2018: input_files = {"/data/disk3/nanoAOD_2018/MuonEGRun2018D/*"}; break;
			default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

        case 114: Process = "Data_MuonEGRunE";

                  switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/MuonEGRun2016E/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/MuonEGRun2017E/*"}; break;
                        case 2018: break;
			default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 115: Process = "Data_MuonEGRunF";

                  switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/MuonEGRun2016F/*"}; break;
                        case 2017: input_files = {"/data/disk0/nanoAOD_2017/MuonEGRun2017F/*"}; break;
                        case 2018: break;
			default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 116: Process = "Data_MuonEGRunG";

                  switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/MuonEGRun2016G/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 117: Process = "Data_MuonEGRunH";

                  switch(YearInt){
                        case 2016: input_files = {"/data/disk2/nanoAOD_2016/MuonEGRun2016H/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 118: Process = "Data_SingleMuonRunB";

                  switch(YearInt){
                        case 2016: input_files = {"/data/disk3/nanoAOD_2016/SingleMuon_NanoAOD25Oct2019_RunB/*"}; break;
                        case 2017: input_files = {"/data/disk3/nanoAOD_2017/SingleMuon_NanoAOD25Oct2019_RunB/*"}; break;
                        case 2018: input_files = {"/data/disk3/nanoAOD_2018/SingleMuon_NanoAOD25Oct_2019_RunB/*"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 119: Process = "Data_SingleMuonRunC";

                  switch(YearInt){
                        case 2016: input_files = {"/data/disk3/nanoAOD_2016/SingleMuon_NanoAOD25Oct2019_RunC/*"}; break;
                        case 2017: input_files = {"/data/disk3/nanoAOD_2017/SingleMuon_NanoAOD25Oct2019_RunC/*"}; break;
                        case 2018: input_files = {"/data/disk3/nanoAOD_2018/SingleMuon_NanoAOD25Oct_2019_RunC/*"}; break;
			default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 120: Process = "Data_SingleMuonRunD";

                  switch(YearInt){
                        case 2016: input_files = {"/data/disk3/nanoAOD_2016/SingleMuon_NanoAOD25Oct2019_RunD/*"}; break;
                        case 2017: input_files = {"/data/disk3/nanoAOD_2017/SingleMuon_NanoAOD25Oct2019_RunD/*"}; break;
                        case 2018: input_files = {"/data/disk3/nanoAOD_2018/SingleMuon_NanoAOD25Oct_2019_RunD/*"}; break;
			default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 121: Process = "Data_SingleMuonRunE";

                  switch(YearInt){
                        case 2016: input_files = {"/data/disk3/nanoAOD_2016/SingleMuon_NanoAOD25Oct2019_RunE/*"}; break;
                        case 2017: input_files = {"/data/disk3/nanoAOD_2017/SingleMuon_NanoAOD25Oct2019_RunE/*"}; break;
			case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 122: Process = "Data_SingleMuonRunF";

                  switch(YearInt){
			case 2016: input_files = {"/data/disk3/nanoAOD_2016/SingleMuon_NanoAOD25Oct2019_RunF/*"}; break;
                        case 2017: input_files = {"/data/disk3/nanoAOD_2017/SingleMuon_NanoAOD25Oct2019_RunF/*"}; break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 123: Process = "Data_SingleMuonRunG";

                  switch(YearInt){
                        case 2016: input_files = {"/data/disk3/nanoAOD_2016/SingleMuon_NanoAOD25Oct2019_RunG/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 124: Process = "Data_SingleMuonRunH";

                  switch(YearInt){
                        case 2016: input_files = {"/data/disk3/nanoAOD_2016/SingleMuon_NanoAOD25Oct2019_RunH/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }
 
	case 125: Process = "Data_SingleElectronRunB";

                  switch(YearInt){
                        case 2016: input_files = {"/data/disk3/nanoAOD_2016/SingleElectron_NanoAOD25Oct2019_RunB/*"}; break;
                        case 2017: input_files = {"/data/disk3/nanoAOD_2017/SingleElectron_NanoAOD25Oct2019_RunB/*"}; break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 126: Process = "Data_SingleElectronRunC";

                  switch(YearInt){
                        case 2016: input_files = {"/data/disk3/nanoAOD_2016/SingleElectron_NanoAOD25Oct2019_RunC/*"}; break;
                        case 2017: input_files = {"/data/disk3/nanoAOD_2017/SingleElectron_NanoAOD25Oct2019_RunC/*"}; break;
                        case 2018: break;
			default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 127: Process = "Data_SingleElectronRunD";

                  switch(YearInt){
			case 2016: input_files = {"/data/disk3/nanoAOD_2016/SingleElectron_NanoAOD25Oct2019_RunD/*"}; break;
                        case 2017: input_files = {"/data/disk3/nanoAOD_2017/SingleElectron_NanoAOD25Oct2019_RunD/*"}; break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 128: Process = "Data_SingleElectronRunE";

                  switch(YearInt){
                        case 2016: input_files = {"/data/disk3/nanoAOD_2016/SingleElectron_NanoAOD25Oct2019_RunE/*"}; break;
                        case 2017: input_files = {"/data/disk3/nanoAOD_2017/SingleElectron_NanoAOD25Oct2019_RunE/*"}; break;
                        case 2018: break;
			default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 129: Process = "Data_SingleElectronRunF";

                  switch(YearInt){
                        case 2016: input_files = {"/data/disk3/nanoAOD_2016/SingleElectron_NanoAOD25Oct2019_RunF/*"}; break;
                        case 2017: input_files = {"/data/disk3/nanoAOD_2017/SingleElectron_NanoAOD25Oct2019_RunF/*"}; break;
                        case 2018: break;
			default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 130: Process = "Data_SingleElectronRunG";

                  switch(YearInt){
                        case 2016: input_files = {"/data/disk3/nanoAOD_2016/SingleElectron_NanoAOD25Oct2019_RunG/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 131: Process = "Data_SingleElectronRunH";

                  switch(YearInt){
                        case 2016: input_files = {"/data/disk3/nanoAOD_2016/SingleElectron_NanoAOD25Oct2019_RunH/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 132: Process = "Data_EGRunB";

		  switch(YearInt){ 
                        case 2016: input_files = {"/data/disk3/nanoAOD_2018/EGammaRunB/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 133: Process = "Data_EGRunC";

                  switch(YearInt){
                        case 2016: input_files = {"/data/disk3/nanoAOD_2018/EGammaRunC/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 134: Process = "Data_EGRunD";

                  switch(YearInt){
                        case 2016: input_files = {"/data/disk3/nanoAOD_2018/EGammaRunD/*"}; break;
                        case 2017: break;
                        case 2018: break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 135: Process = "NPL_File_ee_Blinded";

                  switch(YearInt){
                        case 2016: input_files = {"NPL_ee_output_2016_Blinded.root"}; break;
                        case 2017: input_files = {"NPL_ee_output_2017_Blinded.root"}; break;
                        case 2018: input_files = {"NPL_ee_output_2018_Blinded.root"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 136: Process = "NPL_File_ee_Unblinded";

                  switch(YearInt){
                        case 2016: input_files = {"NPL_ee_output_2016_Unblinded.root"}; break;
                        case 2017: input_files = {"NPL_ee_output_2017_Unblinded.root"}; break;
                        case 2018: input_files = {"NPL_ee_output_2018_Unblinded.root"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	case 137: Process = "NPL_File_mumu_Blinded";

                  switch(YearInt){
                        case 2016: input_files = {"NPL_mumu_output_2016_Blinded.root"}; break;
                        case 2017: input_files = {"NPL_mumu_output_2017_Blinded.root"}; break;
                        case 2018: input_files = {"NPL_mumu_output_2018_Blinded.root"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

        case 138: Process = "NPL_File_mumu_Unblinded";

                  switch(YearInt){
                        case 2016: input_files = {"NPL_mumu_output_2016_Unblinded.root"}; break;
                        case 2017: input_files = {"NPL_mumu_output_2017_Unblinded.root"}; break;
                        case 2018: input_files = {"NPL_mumu_output_2018_Unblinded.root"}; break;
                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018" << std::endl; break;
                  }

	default: std::cout << "ERROR: Double check the value of ProcessInt." << std::endl;
		 break;

  }

  std::cout << '\n' << std::endl;
  std::cout << '\n' << std::endl;
  std::cout << '\n' << std::endl;
  std::cout << '\n' << std::endl;
  std::cout << '\n' << std::endl;
  std::cout << "Process = " << Process << std::endl;
  std::cout << '\n' << std::endl;
  std::cout << '\n' << std::endl;
  std::cout << '\n' << std::endl;
  std::cout << '\n' << std::endl;
  std::cout << '\n' << std::endl;

  switch(ChannelInt){
        
        case 1: Channel = "ee";
        break;
        
        case 2: Channel = "mumu";
        break;

	case 3: Channel = "emu";
        break;

	default: std::cout << "ERROR: Please choose a channel out of 1 (ee), 2 (mumu) or 3 (emu)." << std::endl; break;

  }

  switch(NPLInt){

	case 1: NonPromptLepton = "NPL";
        break;

	default: NonPromptLepton = "";
        break;
        
  }

  switch(SRInt){

	case 1: SignalRegion = "SR";
	break;

	default: SignalRegion = "";
	break;

  }

  switch(SBRInt){

        case 1: SideBandRegion = "SBR";
        break;

        default: SideBandRegion = "";
        break;
  
  }

  switch(SBRInt){

        case 1: SideBandRegion = "SBR";
        break;

        default: SideBandRegion = "";
        break;

  }

  switch(ZPlusJetsCRInt){

	case 1: ZPlusJetsControlRegion = "ZPlusJetsCR";
        break;

        default: ZPlusJetsControlRegion = "";
        break;

  }

  switch(ttbarCRInt){

        case 1: ttbarControlRegion = "ttbarCR";

        	switch(YearInt){
			case 2016: MinElectronPt = 25; MaxElectronPt = 35; MinMuonPt = 25; MaxMuonPt = 26; MaxTrackerEta = 2.4; break;
			case 2017: MinElectronPt = 25; MaxElectronPt = 38; MinMuonPt = 25; MaxMuonPt = 29; MaxTrackerEta = 2.5; break;
			case 2018: MinElectronPt = 25; MaxElectronPt = 38; MinMuonPt = 25; MaxMuonPt = 29; MaxTrackerEta = 2.5; break;
			default: std::cout << "ERROR: Year must be 2016, 2017 or 2018." << std::endl; break;

		}

		break;

        default: ttbarControlRegion = "";

		 switch(YearInt){
                        case 2016: MinElectronPt = 15; MaxElectronPt = 35; MinMuonPt = 20; MaxMuonPt = 26; MaxTrackerEta = 2.4; 
                        case 2017: MinElectronPt = 15; MaxElectronPt = 38; MinMuonPt = 20; MaxMuonPt = 29; MaxTrackerEta = 2.5; break;
                        case 2018: MinElectronPt = 15; MaxElectronPt = 38; MinMuonPt = 20; MaxMuonPt = 29; MaxTrackerEta = 2.5; break;
                        default: std::cout << "ERROR: Year must be 2016, 2017 or 2018." << std::endl; break;
                 }

       		 break;

  }


  switch(YearInt){

	case 2016: Year = "2016"; break;
	case 2017: Year = "2017"; break;
	case 2018: Year = "2018"; break;

	default: std::cout << "ERROR: Year must be 2016, 2017 or 2018." << std::endl;	
		 break;

  }

  switch(SystematicInt){

	case 0: Systematic = "Nominal"; SJER = "sJER_Nominal"; SIGMAJER = "sigma_JER"; break;
	case 1: Systematic = "PU_ScaleUp"; SJER = "sJER_Nominal"; SIGMAJER = "sigma_JER"; break;
	case 2: Systematic = "PU_ScaleDown"; SJER = "sJER_Nominal"; SIGMAJER = "sigma_JER"; break;
	case 3: Systematic = "BTag_ScaleUp"; SJER = "sJER_Nominal"; SIGMAJER = "sigma_JER"; break;
	case 4: Systematic = "BTag_ScaleDown"; SJER = "sJER_Nominal"; SIGMAJER = "sigma_JER"; break;
        case 5: Systematic = "JetSmearing_ScaleUp"; SJER = "sJER_up"; SIGMAJER = "sigma_JER"; break;
        case 6: Systematic = "JetSmearing_ScaleDown"; SJER = "sJER_down"; SIGMAJER = "sigma_JER"; break;
        case 7: Systematic = "JetResolution_ScaleUp"; SJER = "sJER_Nominal"; SIGMAJER = "sigma_JER_up"; break;
        case 8: Systematic = "JetResolution_ScaleDown"; SJER = "sJER_Nominal"; SIGMAJER = "sigma_JER_down"; break;
        case 9: Systematic = "LeptonEfficiencies_ScaleUp"; SJER = "sJER_Nominal"; SIGMAJER = "sigma_JER"; break;
        case 10: Systematic = "LeptonEfficiencies_ScaleDown"; SJER = "sJER_Nominal"; SIGMAJER = "sigma_JER"; break;
        case 11: Systematic = "PDF_ScaleUp"; SJER = "sJER_Nominal"; SIGMAJER = "sigma_JER"; break;
        case 12: Systematic = "PDF_ScaleDown"; SJER = "sJER_Nominal"; SIGMAJER = "sigma_JER"; break;
        case 13: Systematic = "ME_Up"; SJER = "sJER_Nominal"; SIGMAJER = "sigma_JER"; break;
        case 14: Systematic = "ME_Down"; SJER = "sJER_Nominal"; SIGMAJER = "sigma_JER"; break;
        case 15: Systematic = "MET_Up"; SJER = "sJER_Nominal"; SIGMAJER = "sigma_JER"; break;
        case 16: Systematic = "MET_Down"; SJER = "sJER_Nominal"; SIGMAJER = "sigma_JER"; break;
        case 17: Systematic = "isr_up"; SJER = "sJER_Nominal"; SIGMAJER = "sigma_JER"; break;
        case 18: Systematic = "isr_down"; SJER = "sJER_Nominal"; SIGMAJER = "sigma_JER"; break;
        case 19: Systematic = "fsr_up"; SJER = "sJER_Nominal"; SIGMAJER = "sigma_JER"; break;
        case 20: Systematic = "fsr_down"; SJER = "sJER_Nominal"; SIGMAJER = "sigma_JER"; break;
	default: std::cout << "ERROR: SystematicInt must be between 1 and 20." << std::endl; break;

  }



  switch(YearInt){
                
  	case 2016: PSWeightString = "LeptonPt"; 

		   break;

        case 2017: switch(ProcessInt){ 
                   	case 0: PSWeightString = "PSWeight"; break;
                        case 29: PSWeightString = "PSWeight"; break;
                        case 30: PSWeightString = "PSWeight"; break;
                        case 33: PSWeightString = "PSWeight"; break;
                        case 36: PSWeightString = "PSWeight"; break;
                        case 39: PSWeightString = "PSWeight"; break;
                        case 56: PSWeightString = "PSWeight"; break;
                        case 89: PSWeightString = "PSWeight"; break;
                        default: PSWeightString = "LeptonPt"; break;
                    }
	
		   
		   break;

        case 2018: switch(ProcessInt){ 
			case 0: PSWeightString = "PSWeight"; break;
                        case 29: PSWeightString = "PSWeight"; break;
                        case 30: PSWeightString = "PSWeight"; break;
                        case 33: PSWeightString = "PSWeight"; break;
                        case 36: PSWeightString = "PSWeight"; break;
                        case 39: PSWeightString = "PSWeight"; break;
                        case 56: PSWeightString = "PSWeight"; break;
                        case 89: PSWeightString = "PSWeight"; break;
                        default: PSWeightString = "LeptonPt"; break;                  
		    }

		 
		 break;

	default: std::cout << "ERROR: Inside the switch statement for PSWeightString. Please choose the year as 2016, 2017 or 2018." << std::endl; break;
        
  }



  //Lambda functions start here
  //Lambda function for the event cleaning
  auto filter_function{[](const bool& Flag_goodVertices_Selection,    			 const bool& Flag_globalSuperTightHalo2016Filter_Selection, 
			  const bool& Flag_HBHENoiseFilter_Selection, 			 const bool& Flag_HBHENoiseIsoFilter_Selection, 
			  const bool& Flag_EcalDeadCellTriggerPrimitiveFilter_Selection, const bool& Flag_BadPFMuonFilter_Selection, 
			  const bool& Flag_BadChargedCandidateFilter_Selection, 	 const bool& Flag_ecalBadCalibFilter_Selection, 
			  const bool& Flag_eeBadScFilter_Selection)-> bool{


  	//std::cout << "print 1" << std::endl;

	return  Flag_goodVertices_Selection > 0       	     || Flag_globalSuperTightHalo2016Filter_Selection > 0     || Flag_HBHENoiseFilter_Selection > 0 || 
		Flag_HBHENoiseIsoFilter_Selection > 0 	     || Flag_EcalDeadCellTriggerPrimitiveFilter_Selection > 0 || Flag_BadPFMuonFilter_Selection > 0 || 
		Flag_BadChargedCandidateFilter_Selection > 0 || Flag_ecalBadCalibFilter_Selection > 0 		      || Flag_eeBadScFilter_Selection > 0;

  }};


  //Lambda functions for filtering data events using the golden json file. This is to remove events with bad lumisections.
  auto GoldenJsonReader{[&YearInt](){

	//std::cout << "print 2" << std::endl;

  	std::string GoldenJsonFileName;

   	switch(YearInt){

		case 2016: GoldenJsonFileName = "./ScaleFactors/GoldenJSON/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt"; break;
        	case 2017: GoldenJsonFileName = "./ScaleFactors/GoldenJSON/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt"; break;
        	case 2018: GoldenJsonFileName = "./ScaleFactors/GoldenJSON/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt"; break;
        	default: std::cout << "Choose the year out of 2016, 2017 or 2018" << std::endl; break;

   	}

   	std::ifstream myReadFile;
   	myReadFile.open(GoldenJsonFileName);
   	static char output[100];
   	std::vector<std::string> GoldenJsonOutput{};
  
   	if (myReadFile.is_open()) { while (!myReadFile.eof()) {myReadFile >> output; GoldenJsonOutput.push_back(output);}}

  	myReadFile.close();
  	return GoldenJsonOutput;

  }};


  auto GoldenJson_SplitChars{[&YearInt, &GoldenJsonReader](){

	//std::cout << "print 3" << std::endl;

  	std::vector<char> out{};

    	for(long unsigned int i = 0; i < (GoldenJsonReader()).size(); i++){

  		std::string element = GoldenJsonReader().at(i);
  	  	for(long unsigned int j = 0; j < element.size(); j++){out.push_back(element.at(j));}

    	}

    	return out;

  }};
 

  auto RunNumberCheck{[&YearInt, &GoldenJson_SplitChars](const unsigned int& InputRunNumber){

	//std::cout << "print 4" << std::endl;

  	std::vector<char> EventsVector{}; 

  	for(long unsigned int i = 0; i < (GoldenJson_SplitChars()).size(); i++){

		unsigned int RunNumBeingRead;

 		if(  GoldenJson_SplitChars().at(i+1) == '"' && (GoldenJson_SplitChars().at(i+2) == '2' || GoldenJson_SplitChars().at(i+2) == '3')  ){ 

			int digit1 = GoldenJson_SplitChars().at(i+2) - '0'; 
			int digit2 = GoldenJson_SplitChars().at(i+3) - '0'; 
			int digit3 = GoldenJson_SplitChars().at(i+4) - '0'; 
			int digit4 = GoldenJson_SplitChars().at(i+5) - '0'; 
			int digit5 = GoldenJson_SplitChars().at(i+6) - '0'; 
			int digit6 = GoldenJson_SplitChars().at(i+7) - '0';

			int run = (digit1*100000) + (digit2*10000) + (digit3*1000) + (digit4*100) + (digit5*10) + digit6;
			RunNumBeingRead = run;

			if(run == InputRunNumber){

				for(long unsigned int j = 2; j < (GoldenJson_SplitChars()).size(); j++){

					if(GoldenJson_SplitChars().at(i+10) == '[' && GoldenJson_SplitChars().at(i+11) == '['){	

						if( GoldenJson_SplitChars().at( (i+10)+j ) == ']' && GoldenJson_SplitChars().at( (i+10)+(j+1) ) == ']'){

							for(int k = (i+10); k < ((i+10)+(j+2)); k++){EventsVector.push_back(GoldenJson_SplitChars().at(k));}
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


  auto ReturnRunNumAndEventRanges{[&YearInt, &RunNumberCheck](const unsigned int& InputRunNumber){

	//std::cout << "print 5" << std::endl;

   	std::vector<int> RunNumAndEvents{};
  	RunNumAndEvents.push_back(InputRunNumber);
   	std::vector<char> Runs = RunNumberCheck(InputRunNumber);

   	for(long unsigned int i = 0; i < Runs.size(); i++){

 		if(Runs.at(i) == ']' && Runs.at(i+1) == ']'){break;}
		else if( isdigit(Runs.at(i)) || Runs.at(i) == ',' || Runs.at(i) == ' ' || Runs.at(i) == ']'){continue;}
 		else if( (Runs.at(i) == '[' && Runs.at(i+1) == '[') || (Runs.at(i) == '[' && isdigit(Runs.at(i+1))) ){

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


 
  auto RunAndLumiFilterFunction{[&ReturnRunNumAndEventRanges, &MCInt](const unsigned int& InputRunNumber, const unsigned int& luminosityBlock){

      //std::cout << "print 6" << std::endl;	

      switch(MCInt){
      
          case 0: if( InputRunNumber == ReturnRunNumAndEventRanges(InputRunNumber).at(0) ){
                  	for(long unsigned int i = 1; i < ReturnRunNumAndEventRanges(InputRunNumber).size(); i+=2){

                        	int MinLumi = ReturnRunNumAndEventRanges(InputRunNumber).at(i);
                        	int MaxLumi = ReturnRunNumAndEventRanges(InputRunNumber).at(i+1);

                        	if(luminosityBlock > MinLumi && luminosityBlock < MaxLumi){return InputRunNumber && luminosityBlock;}
                        	else{continue;}

                    	}

                 }
                 else{return false;}
              
          case 1: return  InputRunNumber && luminosityBlock;
              
      }


  }};  

  //To prevent double counting datasets
  auto DoubleCountCheckLeptonTriggers{[&YearInt](const int& TriggerType,
					 	 const bool& HLT_Ele25_eta2p1_WPTight_Gsf, 		         const bool& HLT_Ele27_WPTight_Gsf, 
					 	 const bool& HLT_Ele32_eta2p1_WPTight_Gsf, 		         const bool& HLT_Ele32_WPTight_Gsf_L1DoubleEG,
					 	 const bool& HLT_Ele35_WPTight_Gsf, 	   		         const bool& HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL,
					 	 const bool& HLT_IsoMu24, 		   		         const bool& HLT_IsoMu24_eta2p1, 
					 	 const bool& HLT_IsoMu27, 		   		         const bool& HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ,
					 	 const bool& HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8,          const bool& HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8,
					 	 const bool& HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,          const bool& HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,
					 	 const bool& HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, const bool& HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,
					 	 const bool& HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,    const bool& HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,
					 	 const bool& HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL){


	//std::cout << "print 7" << std::endl;

	switch(TriggerType){

		case 0: switch(YearInt){//single electron
						case 2016: return (HLT_Ele25_eta2p1_WPTight_Gsf > 0 || HLT_Ele27_WPTight_Gsf > 0 || HLT_Ele32_eta2p1_WPTight_Gsf > 0);
						case 2017: return (HLT_Ele32_WPTight_Gsf_L1DoubleEG > 0 || HLT_Ele35_WPTight_Gsf > 0);
						case 2018: return (HLT_Ele32_WPTight_Gsf_L1DoubleEG > 0 || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL > 0);
						default: std::cout << "Please choose a year out of 2016, 2017 or 2018." << std::endl;
					}

		case 1: switch(YearInt){//single muon
                                   	case 2016: return (HLT_IsoMu24 <= 0 || HLT_IsoMu24_eta2p1 <= 0);
                                        case 2017: return (HLT_IsoMu27 <= 0);
                                        case 2018: return HLT_IsoMu24 <= 0;
                                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018." << std::endl;
                                   }

		case 2: switch(YearInt){//double muon
                                        case 2016: return (HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ <= 0 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 <= 0);
                                        case 2017: return (HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ <= 0 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 <= 0);
                                        case 2018: return (HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 <= 0 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 <= 0);
                                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018." << std::endl;
                                   }

		case 3: switch(YearInt){//double electron
                                        case 2016: return HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ > 0;
                                        case 2017: return (HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL > 0 || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ > 0);
                                        case 2018: return HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ > 0;
                                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018." << std::endl;
                                   }

		case 4: switch(YearInt){//muon EG
                                        case 2016: return (HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || 
                                            	   	  HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || 
                                            		  HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || 
                                            		  HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL <= 0 || 
                                            		  HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL <= 0 || 
                                            		  HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL <= 0);

                                        case 2017: return (HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || 
                                            		   HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0 ||
                                            		   HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0); 
                                        
					case 2018: return (HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || 
                                           		   HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || 
                                           		   HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0); 

                                        default: std::cout << "Please choose a year out of 2016, 2017 or 2018." << std::endl;

                                   }

		default: std::cout << "Choose a case out of SingleElectron, SingleMuon, DoubleElectron, DoubleMuon or MuonEG." << std::endl; break;

	}

  }};

  auto DoubleCountCheck_EventFunction{[&Process, &DoubleCountCheckLeptonTriggers](const bool& HLT_Ele25_eta2p1_WPTight_Gsf,
									          const bool& HLT_Ele27_WPTight_Gsf,
                                   	         				  const bool& HLT_Ele32_eta2p1_WPTight_Gsf,                       
										  const bool& HLT_Ele32_WPTight_Gsf_L1DoubleEG,
                                   	         			          const bool& HLT_Ele35_WPTight_Gsf,                             
										  const bool& HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL,
                                   	        				  const bool& HLT_IsoMu24,       
										  const bool& HLT_IsoMu24_eta2p1,
                                   	         				  const bool& HLT_IsoMu27,                                        
										  const bool& HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ,
                                   	         				  const bool& HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8,          
										  const bool& HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8,
                                   	         				  const bool& HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,          
										  const bool& HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,
                                   	         				  const bool& HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, 
										  const bool& HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,
                                   	         			          const bool& HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,    
										  const bool& HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,
                                   	         				  const bool& HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, 
										  const ULong64_t& event){


	//std::cout << "print 8" << std::endl;

	bool Single_E = DoubleCountCheckLeptonTriggers(0, 
				       		      HLT_Ele25_eta2p1_WPTight_Gsf, 			    HLT_Ele27_WPTight_Gsf, 
				       		      HLT_Ele32_eta2p1_WPTight_Gsf, 			    HLT_Ele32_WPTight_Gsf_L1DoubleEG,
                                       		      HLT_Ele35_WPTight_Gsf,        			    HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL,
                                       		      HLT_IsoMu24,          			  	    HLT_IsoMu24_eta2p1,
                                       		      HLT_IsoMu27,                                         HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ,
                                       		      HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8,           HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8,
                                       		      HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,           HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,
                                       		      HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,  HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,
                                       		      HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,     HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,
                                       		      HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);

	bool Single_Mu = DoubleCountCheckLeptonTriggers(1,
				       		      HLT_Ele25_eta2p1_WPTight_Gsf,                        HLT_Ele27_WPTight_Gsf,
                                       		      HLT_Ele32_eta2p1_WPTight_Gsf,                        HLT_Ele32_WPTight_Gsf_L1DoubleEG,
                                       		      HLT_Ele35_WPTight_Gsf,                               HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL,
                                       		      HLT_IsoMu24,                                         HLT_IsoMu24_eta2p1,
                                       		      HLT_IsoMu27,                                         HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ,
                                       		      HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8,           HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8,
                                       		      HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,           HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,
                                       		      HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,  HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,
                                       		      HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,     HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,
                                       		      HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);

	bool Double_E = DoubleCountCheckLeptonTriggers(2,
				       		      HLT_Ele25_eta2p1_WPTight_Gsf,                        HLT_Ele27_WPTight_Gsf,
                                       		      HLT_Ele32_eta2p1_WPTight_Gsf,                        HLT_Ele32_WPTight_Gsf_L1DoubleEG,
                                       		      HLT_Ele35_WPTight_Gsf,                               HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL,
                                       		      HLT_IsoMu24,                                         HLT_IsoMu24_eta2p1,
                                       		      HLT_IsoMu27,                                         HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ,
                                       		      HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8,           HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8,
                                      		      HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,           HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,
                                       		      HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,  HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,
                                       		      HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,     HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,
                                       		      HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);

	bool Double_Mu = DoubleCountCheckLeptonTriggers(3,
				       		      HLT_Ele25_eta2p1_WPTight_Gsf,                        HLT_Ele27_WPTight_Gsf,
                                       		      HLT_Ele32_eta2p1_WPTight_Gsf,                        HLT_Ele32_WPTight_Gsf_L1DoubleEG,
                                       		      HLT_Ele35_WPTight_Gsf,                               HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL,
                                       		      HLT_IsoMu24,                                         HLT_IsoMu24_eta2p1,
                                       		      HLT_IsoMu27,                                         HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ,
                                       		      HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8,           HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8,
                                       		      HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,           HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,
                                       		      HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,  HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,
                                       		      HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,     HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,
                                       		      HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);

	bool ElecMu = DoubleCountCheckLeptonTriggers(4,
				       		     HLT_Ele25_eta2p1_WPTight_Gsf,                        HLT_Ele27_WPTight_Gsf,
                                       		     HLT_Ele32_eta2p1_WPTight_Gsf,                        HLT_Ele32_WPTight_Gsf_L1DoubleEG,
                                       	             HLT_Ele35_WPTight_Gsf,                               HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL,
                                       		     HLT_IsoMu24,                                         HLT_IsoMu24_eta2p1,
                                       		     HLT_IsoMu27,                                         HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ,
                                       		     HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8,           HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8,
                                       		     HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,           HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,
                                       		     HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,  HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ,
                                       		     HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,     HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,
                                       		     HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);

	if(Process == "Data_DoubleEGRunB"){return Double_E == true;}
	else if(Process == "Data_DoubleEGRunC"){return Double_E == true;}
	else if(Process == "Data_DoubleEGRunD"){return Double_E == true;}
	else if(Process == "Data_DoubleEGRunE"){return Double_E == true;}
	else if(Process == "Data_DoubleEGRunF"){return Double_E == true;}
	else if(Process == "Data_DoubleEGRunG"){return Double_E == true;}
	else if(Process == "Data_DoubleEGRunH"){return Double_E == true;}
	else if(Process == "Data_SingleElectronRunB"){return Single_E == true;}
	else if(Process == "Data_SingleElectronRunC"){return Single_E == true;}
	else if(Process == "Data_SingleElectronRunD"){return Single_E == true;}
	else if(Process == "Data_SingleElectronRunE"){return Single_E == true;}
	else if(Process == "Data_SingleElectronRunF"){return Single_E == true;}
	else if(Process == "Data_SingleElectronRunG"){return Single_E == true;}
	else if(Process == "Data_SingleElectronRunH"){return Single_E == true;}
	else if(Process == "Data_DoubleMuonRunB"){return Double_Mu == true;}
	else if(Process == "Data_DoubleMuonRunC"){return Double_Mu == true;}
	else if(Process == "Data_DoubleMuonRunD"){return Double_Mu == true;}
	else if(Process == "Data_DoubleMuonRunE"){return Double_Mu == true;}
	else if(Process == "Data_DoubleMuonRunF"){return Double_Mu == true;}
	else if(Process == "Data_DoubleMuonRunG"){return Double_Mu == true;}
	else if(Process == "Data_DoubleMuonRunH"){return Double_Mu == true;}
        else if(Process == "Data_SingleMuonRunB"){return Single_Mu == true;}
	else if(Process == "Data_SingleMuonRunC"){return Single_Mu == true;}
	else if(Process == "Data_SingleMuonRunD"){return Single_Mu == true;}
	else if(Process == "Data_SingleMuonRunE"){return Single_Mu == true;}
	else if(Process == "Data_SingleMuonRunF"){return Single_Mu == true;}
	else if(Process == "Data_SingleMuonRunG"){return Single_Mu == true;}
	else if(Process == "Data_SingleMuonRunH"){return Single_Mu == true;}
	else if(Process == "Data_MuonEGRunB"){return ElecMu == true;}
	else if(Process == "Data_MuonEGRunC"){return ElecMu == true;}
	else if(Process == "Data_MuonEGRunD"){return ElecMu == true;}
	else if(Process == "Data_MuonEGRunE"){return ElecMu == true;}
	else if(Process == "Data_MuonEGRunF"){return ElecMu == true;}
	else if(Process == "Data_MuonEGRunG"){return ElecMu == true;}
	else if(Process == "Data_MuonEGRunH"){return ElecMu == true;}
	else{return event > 0;}


  }};


  //Lambda function for the pile up modelling
  auto PU_function{[&puReweight_2016, &puReweight_2016_part1, &puReweight_2016_part2, &puReweight_2017, &puReweight_2018, &YearInt](int PV_npvs_input){

  	////std::cout << "print 9" << std::endl;

      	float PU_Weight_input;

      	switch(YearInt){

        	case 2016: PU_Weight_input = puReweight_2016->GetBinContent(puReweight_2016->GetXaxis()->FindBin(PV_npvs_input)); break;
          	case 2017: PU_Weight_input = puReweight_2017->GetBinContent(puReweight_2017->GetXaxis()->FindBin(PV_npvs_input)); break;
          	case 2018: PU_Weight_input = puReweight_2018->GetBinContent(puReweight_2018->GetXaxis()->FindBin(PV_npvs_input)); break;
          	default: std::cout << "Choose a year out of 2016, 2017 or 2018 for the PU function" << std::endl; break;

      	}

      	return PU_Weight_input;

  }};

  //Lambda functions for the electron selection
  auto ElectronsFunction{[](const int targetID, const floats& Electron_pt, const floats& Electron_eta, const ints& Electron_cutBased, const bools& Electron_isPFcand){
 
  	//std::cout << "print 9" << std::endl;
  	return (Electron_pt > MinElectronPt && (abs(Electron_eta) < MaxTrackerEta && (abs(Electron_eta) < 1.442 || abs(Electron_eta) > 1.566) ) && 
		Electron_cutBased >= targetID && Electron_isPFcand);

  }};

  auto MuonsFunction{[](const float target_iso, const bools& isPFs, const floats& Muon_pt, const floats& Muon_eta, const bools& ids, const floats& isos){

  	//std::cout << "print 10" << std::endl;
  	return (isPFs && Muon_pt > MinMuonPt && abs(Muon_eta) < MaxTrackerEta && ids && isos <= target_iso);

  }};



  auto TightLeptonsFunction{[&ChannelInt, &ElectronsFunction, &MuonsFunction](const floats& Electron_pt,      const floats& Electron_eta, const ints& Electron_cutBased, 
									      const bools& Electron_isPFcand, const bools& isPFs,         const floats& pts, 
									      const floats& etas,             const bools& ids,           const floats& isos){

  	//std::cout << "print 11" << std::endl;

	switch(ChannelInt){

  		case 1: return ElectronsFunction(4, Electron_pt, Electron_eta, Electron_cutBased, Electron_isPFcand);
		case 2: return MuonsFunction(0.25, isPFs, pts, etas, ids, isos);
		case 3: return ElectronsFunction(4, Electron_pt, Electron_eta, Electron_cutBased, Electron_isPFcand) || MuonsFunction(0.25, isPFs, pts, etas, ids, isos);
		default: std::cout << "Choose the ee, mumu or emu channel" << std::endl; break; 

	}

  }};



  auto LeptonVariableFunctionFloats{[&ChannelInt](const floats& Electron_input, const floats& Muon_input){

	//std::cout << "print 12" << std::endl; 

	floats Emu_vector_floats{};

	switch(ChannelInt){

		case 1: return Electron_input; 
		case 2: return Muon_input;

		case 3: for(int i = 0; i < Electron_input.size(); i++){Emu_vector_floats.push_back(Electron_input.at(i));}
			for(int i = 0; i < Muon_input.size(); i++){Emu_vector_floats.push_back(Muon_input.at(i));}
			return Emu_vector_floats;
	
		default: std::cout << "ERROR: Choose ChannelInt to be 1 (for ee), 2 (for mumu) or 3 (for emu)." << std::endl; break;

	}

  }};

  auto LeptonVariableFunctionInts{[&ChannelInt](const ints& Electron_input, const ints& Muon_input){

        //std::cout << "print 13" << std::endl;

	ints Emu_vector_ints{};

        switch(ChannelInt){

                case 1: return Electron_input;
                case 2: return Muon_input;

                case 3: for(int i = 0; i < Electron_input.size(); i++){Emu_vector_ints.push_back(Electron_input.at(i));}
                        for(int i = 0; i < Muon_input.size(); i++){Emu_vector_ints.push_back(Muon_input.at(i));}
                        return Emu_vector_ints;

                default: std::cout << "ERROR: Choose ChannelInt to be 1 (for ee), 2 (for mumu) or 3 (for emu)." << std::endl; break;

        }

  }};


  auto LeptonVariableFunctionChars{[&ChannelInt](const chars& Electron_input, const chars& Muon_input){

        //std::cout << "print 14" << std::endl;

	chars Emu_vector_chars{};

        switch(ChannelInt){

                case 1: return Electron_input;
                case 2: return Muon_input;

                case 3: for(int i = 0; i < Electron_input.size(); i++){Emu_vector_chars.push_back(Electron_input.at(i));}
                        for(int i = 0; i < Muon_input.size(); i++){Emu_vector_chars.push_back(Muon_input.at(i));}
                        return Emu_vector_chars;

                default: std::cout << "ERROR: Choose ChannelInt to be 1 (for ee), 2 (for mumu) or 3 (for emu)." << std::endl; break;

        }

  }};

  auto LooseLeptonsFunction{[&ChannelInt, &ElectronsFunction, &MuonsFunction](const floats& Electron_pt,      const floats& Electron_eta, const ints& Electron_cutBased, 
									      const bools& Electron_isPFcand, const bools& isPFs,         const floats& pts, 
									      const floats& etas,             const bools& ids,           const floats& isos){

  	//std::cout << "print 15" << std::endl;

	switch(ChannelInt){

  		case 1: return ElectronsFunction(1, Electron_pt, Electron_eta, Electron_cutBased, Electron_isPFcand);
		case 2: return MuonsFunction(0.15, isPFs, pts, etas, ids, isos);
		case 3: return ElectronsFunction(1, Electron_pt, Electron_eta, Electron_cutBased, Electron_isPFcand) || MuonsFunction(0.15, isPFs, pts, etas, ids, isos);
		default: std::cout << "Choose the ee, mumu or emu channel" << std::endl; break;

	}

  }};

  auto OppositeSign{[&ChannelInt](const ints& charges){

  	//std::cout << "print 16" << std::endl;

 	return charges.size() == 2 ? signbit(charges.at(0)) != signbit(charges.at(1)) : false;

  }};

  auto SameSign{[&ChannelInt](const ints& charges){

  	//std::cout << "print 14" << std::endl;
	return charges.size() == 2 ? signbit(charges.at(0)) == signbit(charges.at(1)) : false;

  }};

  auto LeadingVariable{[&ChannelInt](const floats& variable){

  	//std::cout << "print 17" << std::endl;

	if(variable.size() > 0){

  		float first_largest_value = variable.at(0);

        	for(long unsigned int i = 1; i < variable.size(); i++){

                	if(variable.at(i) > first_largest_value){
                        	first_largest_value = variable.at(i);
                	}

        	}

  		return first_largest_value;

  	}
  	
	else{float zero = 0.0; return zero;}


  }};


  auto SubleadingVariable{[&ChannelInt](const floats& variable){

  	//std::cout << "print 18" << std::endl;

  	if(variable.size() == 0){float zero = 0.0; return zero;}
  	else{

  		if(variable.size() > 1){

  			float first_largest_value = variable.at(0);

	  		for(long unsigned int i = 1; i < variable.size(); i++){
				if(variable.at(i) > first_largest_value){
					first_largest_value = variable.at(i);
				}
	  		}

  			float second_largest_value = INT_MIN;

	  		for(long unsigned int i = 0; i < variable.size(); i++){
				if( (variable.at(i) > second_largest_value) && (variable.at(i) < first_largest_value) ){
					second_largest_value = variable.at(i);
				}
	
	  		}

  			return second_largest_value;

  		}
  		else{return variable.at(0);}

  	}



  }};

  auto ThirdLeadingVariable{[](const floats& variable){

  	//std::cout << "print 19" << std::endl;

  	if(variable.size() > 2){

  		float first_largest_value = variable.at(0);

        	for(long unsigned int i = 1; i < variable.size(); i++){

                	if(variable.at(i) > first_largest_value){
                        	first_largest_value = variable.at(i);
                	}

        	}

  		float second_largest_value = INT_MIN;

        	for(long unsigned int i = 0; i < variable.size(); i++){

                	if( (variable.at(i) > second_largest_value) && (variable.at(i) < first_largest_value) ){
                        	second_largest_value = variable.at(i);
                	}
        
        	}


  		float third_largest_value = INT_MIN;

		for(long unsigned int i = 0; i < variable.size(); i++){
	
			if( (variable.at(i) > third_largest_value) && (variable.at(i) < second_largest_value) ){
				third_largest_value = variable.at(i);			
			}
	
		}


  		return third_largest_value;

  	}

  }};



  auto FourthLeadingVariable{[](const floats& variable){

  	//std::cout << "print 20" << std::endl;

  	if(variable.size() > 3){

  		float first_largest_value = variable.at(0);

        	for(long unsigned int i = 1; i < variable.size(); i++){

                	if(variable.at(i) > first_largest_value){
                        	first_largest_value = variable.at(i);
                	}

        	}

  		float second_largest_value = INT_MIN;

        	for(long unsigned int i = 0; i < variable.size(); i++){

                	if( (variable.at(i) > second_largest_value) && (variable.at(i) < first_largest_value) ){
                        	second_largest_value = variable.at(i);
                	}

        	}


  		float third_largest_value = INT_MIN;

        	for(long unsigned int i = 0; i < variable.size(); i++){

                	if( (variable.at(i) > third_largest_value) && (variable.at(i) < second_largest_value) ){
                        	third_largest_value = variable.at(i);
                	}

        	}


  		float fourth_largest_value = INT_MIN;

		for(long unsigned int i = 0; i < variable.size(); i++){

                	if( (variable.at(i) > fourth_largest_value) && (variable.at(i) < third_largest_value) ){
                        	fourth_largest_value = variable.at(i);
                	}

        	}


  		return fourth_largest_value;

  	}


  }};



  auto MET_Triggers_Function{[&YearInt, &ProcessInt](const bool& HLT_MET200, 
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
						     const bool& HLT_PFHT800_PFMET85_PFMHT85_IDTight, 
						     const ULong64_t& event){


	//std::cout << "print 29" << std::endl;

	switch(ProcessInt){

		case 95: switch(YearInt){

				case 2016: return HLT_MET200 > 0 || HLT_MET250 > 0 || 
					          HLT_PFMET120_PFMHT120_IDTight > 0 || HLT_PFMET170_HBHECleaned > 0 || HLT_PFHT300_PFMET100 > 0;

  				default: return //HLT_MET105_IsoTrk50 > 0 || //(not in run B)
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

		case 96: switch(YearInt){

                                case 2016: return HLT_MET200 > 0 || HLT_MET250 > 0 ||
                                                  HLT_PFMET120_PFMHT120_IDTight > 0 || HLT_PFMET170_HBHECleaned > 0 || HLT_PFHT300_PFMET100 > 0;

                                default: return //HLT_MET105_IsoTrk50 > 0 || //(not in run B)
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

	      default: return event > 0;
  	
	}


  }};


  auto Lepton_Triggers_Function{[&ChannelInt, &YearInt, &ProcessInt](const bool& HLT_Ele32_WPTight_Gsf_L1DoubleEG, 
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
								     const bool& HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8,
								     const ULong64_t& event

  ){

	//std::cout << "print 30" << std::endl;

	if(ProcessInt == 95 || 96){

	  switch(ChannelInt){

		case 1: switch(YearInt){
				case 2016: return //single or double electron and not any of the others
	
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
				
				case 2017: return
					
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

				case 2018: return 

					   (HLT_Ele32_WPTight_Gsf_L1DoubleEG > 0 || //single electron
					    HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL > 0 || //single electron
					    HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ > 0) //double electron

					   &&

					  (HLT_IsoMu24 <= 0 || //single muon
		 			   HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 <= 0 || //double muon
		 			   HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 <= 0 || //double muon
		 			   HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || //muon+electron
		 		           HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || //muon+electron
		 			   HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0); //muon+electron


				default: std::cout << "The year must be 2016, 2017 or 2018." << std::endl; break;
			}

		case 2: switch(YearInt){
                                case 2016: return //single or double muon and not any of the others

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

                                case 2017: return
 
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

                                case 2018: return //single or double muon and not any of the others
	
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

                                default: std::cout << "The year must be 2016, 2017 or 2018." << std::endl; break;
                        }

		case 3: switch(YearInt){
                                case 2016: return 

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

                                case 2017: return 

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

                                case 2018: return //single lepton or muon+electron and not any of the others

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

                                default: std::cout << "The year must be 2016, 2017 or 2018." << std::endl; break;
                        }

		default: std::cout << "ChannelInt must be 1 (for ee), 2 (for mumu) or 3 (for emu)." << std::endl; break;


	  }

	}
	else{return event > 0;}

  }};


  auto NumberOfLeptonsFunction{[&ChannelInt](const unsigned int& nElectron, const unsigned int& nMuon){

	switch(ChannelInt){
		case 1: return nElectron == 2;
		case 2: return nMuon == 2;
		case 3: return nElectron == 1 && nMuon == 1;
		default: std::cout << "ChannelInt must be 1 (for ee), 2 (for mumu), or 3 (for emu)" << std::endl;
	}

  }};

  auto LeptonCut{[&ChannelInt](const floats& tight_lepton_pts,   const floats& loose_lepton_pts,   const bool& os, 
			       const unsigned int& nElectron,    const unsigned int& nMuon,        const floats& Electron_dz,   
			       const floats& Electron_dxy,       const float& LeadingLeptonPt,     const float& SubleadingLeptonPt, 
			       const floats& LeptonEta){

  	//std::cout << "print 31" << std::endl;

  	const bool lepton_cut{tight_lepton_pts.size() == 2 && tight_lepton_pts.size() == loose_lepton_pts.size()};
  	bool lead_pt_cut{false};

	float MaxLeptonPt;

	switch(ChannelInt){
		case 1: MaxLeptonPt = MaxElectronPt; break;
		case 2: MaxLeptonPt = MaxMuonPt; break;
		case 3: MaxLeptonPt = MaxElectronPt; break; //both MaxMuonPt and MaxElectronPt are 25 for the emu channel
		default: std::cout << "ERROR: Channel must be ee, mumu or emu." << std::endl; break;
	}

	lead_pt_cut = tight_lepton_pts.empty() ? false : *max_element(tight_lepton_pts.begin(), tight_lepton_pts.end()) > MaxLeptonPt;


	switch(ChannelInt){

		case 1: return 
			os 				         &&  lead_pt_cut 			 && 
			lepton_cut 				 &&  nElectron == 2;/*			 &&
			abs(LeptonEta) < 1.442			 &&  abs(LeptonEta) > 1.566		 && 
  			LeadingLeptonPt > MaxElectronPt          &&  SubleadingLeptonPt > MinElectronPt  && 
			(abs(LeptonEta) < 1.442     		 &&  Electron_dz < 0.1 			 && Electron_dxy < 0.05)  || //barrel region
                        (abs(LeptonEta) > 1.566                  &&  abs(LeptonEta) < 3.0   		 && Electron_dz < 0.2 && Electron_dxy < 0.1); //endcaps
*/

		case 2: return os && lead_pt_cut && lepton_cut && nMuon == 2;

		case 3: return os && lepton_cut && nElectron == 1 && nMuon == 1;

		default: std::cout << "ChannelInt must be 1 (for ee), 2 (for mumu) or 3 (for emu)." << std::endl; break;

	}


  }};
 
  auto OppositeSignNonPrompt{[](const ints& charges, const chars& Lepton_genPartFlav){

  	//std::cout << "print 32" << std::endl;
  	bool OppositeSignChargeCheck = charges.size() == 2 ? signbit(charges.at(0)) != signbit(charges.at(1)) : false;
  	bool LeptonNonPromptCheck = all_of(Lepton_genPartFlav.begin(), Lepton_genPartFlav.end(), [](int i){return i != 1;});

  	return OppositeSignChargeCheck && (LeptonNonPromptCheck == 1);

  }};

  auto OppositeSignPrompt{[](const ints& charges, const chars& Lepton_genPartFlav){

  	//std::cout << "print 33" << std::endl;
  	bool OppositeSignChargeCheck = charges.size() == 2 ? signbit(charges.at(0)) != signbit(charges.at(1)) : false;
  	bool LeptonPromptCheck = all_of(Lepton_genPartFlav.begin(), Lepton_genPartFlav.end(), [](int i){return i == 1;});

  	return OppositeSignChargeCheck && (LeptonPromptCheck == 1);

  }};


  auto SameSignNonPrompt{[](const ints& charges, const chars& Lepton_genPartFlav){
  
  	//std::cout << "print 34" << std::endl;

  	bool SameSignChargeCheck = charges.size() == 2 ? signbit(charges.at(0)) == signbit(charges.at(1)) : false;  
  	bool LeptonNonPromptCheck = all_of(Lepton_genPartFlav.begin(), Lepton_genPartFlav.end(), [](int i){return i != 1;});

  	return SameSignChargeCheck && (LeptonNonPromptCheck == 1);

  }};

  auto SameSignPrompt{[](const ints& charges, const chars& Lepton_genPartFlav){

  	//std::cout << "print 35" << std::endl;  

  	bool SameSignChargeCheck = charges.size() == 2 ? signbit(charges.at(0)) == signbit(charges.at(1)) : false;
  	bool LeptonPromptCheck = all_of(Lepton_genPartFlav.begin(), Lepton_genPartFlav.end(), [](int i){return i == 1;});  

  	return SameSignChargeCheck && (LeptonPromptCheck == 1);

  }};


  auto inv_mass{[](const floats& pts, const floats& etas, const floats& phis, const floats& ms)
   {

	//std::cout << "print 36" << std::endl;

    	if (!all_equal(pts.size(), etas.size(), phis.size(), ms.size())){throw std::logic_error("Collections must be the same size");}
    	else if(pts.empty()){throw std::logic_error("Collections must not be empty");}

    	TLorentzVector vec{};

    	for (size_t i{0}; i < pts.size(); i++){
        	TLorentzVector p{};
        	p.SetPtEtaPhiM(pts[i], etas[i], phis[i], ms[i]);
        	vec += p;
    	}
	
    	return boost::numeric_cast<float>(vec.M());

  }};

  auto RecoZ{[](const float& LeadingleptonPt,    const float& LeadingleptonEta,    const float& LeadingleptonPhi,    const float& LeadingleptonMass,
		const float& SubleadingleptonPt, const float& SubleadingleptonEta, const float& SubleadingleptonPhi, const float& SubleadingleptonMass){

  	//std::cout << "print 37" << std::endl;

  	TLorentzVector ZBoson = {};
  	TLorentzVector LeadingLepton = {};
  	TLorentzVector SubleadingLepton = {};

  	LeadingLepton.SetPtEtaPhiM(LeadingleptonPt, LeadingleptonEta, LeadingleptonPhi, LeadingleptonMass);
  	SubleadingLepton.SetPtEtaPhiM(SubleadingleptonPt, SubleadingleptonEta, SubleadingleptonPhi, SubleadingleptonMass);

  	ZBoson = LeadingLepton + SubleadingLepton;

  	return ZBoson;

  }};

  
  auto TLorentzVectorVariable{[](const int& VariableChoice, const TLorentzVector& object){

  	//std::cout << "print 38" << std::endl;

  	doubles vec{};
	
	switch(VariableChoice){
		case 1: vec.push_back(object.Pt()); break;
		case 2: vec.push_back(object.Phi()); break;
		case 3: vec.push_back(object.Eta()); break;
		case 4: vec.push_back(object.M()); break;
		default: break;
	}

  	return vec;

  }};

  auto TLorentzVectorVariablePt{[&TLorentzVectorVariable](const TLorentzVector& object){
	//std::cout << "print 158" << std::endl; 
	return TLorentzVectorVariable(1, object);
  }};
  
  auto TLorentzVectorVariablePhi{[&TLorentzVectorVariable](const TLorentzVector& object){
	//std::cout << "print 159" << std::endl; 
	return TLorentzVectorVariable(2, object);
  }};

  auto TLorentzVectorVariableEta{[&TLorentzVectorVariable](const TLorentzVector& object){
	//std::cout << "print 160" << std::endl;
	return TLorentzVectorVariable(3, object);
  }};
  
  auto TLorentzVectorVariableMass{[&TLorentzVectorVariable](const TLorentzVector& object){
	//std::cout << "print 161" << std::endl; 
	return TLorentzVectorVariable(4, object);
  }};


  auto deltaRcheck_float{[](const float& Object1_eta, const float& Object1_phi, const float& Object2_eta, const float& Object2_phi){

  	//std::cout << "print 41" << std::endl;

  	float dR = sqrt(pow(Object1_eta - Object2_eta, 2) + pow(Object1_phi - Object2_phi, 2));
  	return dR;

  }};


  auto DeltaPhi_floatandfloat{[](const float& Object1_phi, const float& Object2_phi){

  	//std::cout << "print 42" << std::endl;

  	double dPhi = abs(Object1_phi - Object2_phi);
  	return dPhi;

  }};

  auto LeptonFourMomentumFunction{[](const floats& Muon_pt, const floats& Muon_eta, const floats& Muon_phi, const floats& Muon_mass){
  
  	//std::cout << "print 43" << std::endl;

  	TLorentzVector Muon4Mo{};
  
  	for(long unsigned int i = 0; i < Muon_pt.size(); i++){
  		TLorentzVector vec{};
  		vec.SetPtEtaPhiM(Muon_pt.at(i), Muon_eta.at(i), Muon_phi.at(i), Muon_mass.at(i));
  		Muon4Mo += vec;
  	}

  	return Muon4Mo;

  }};

  auto RochesterCorrections_testscript2{[](const int& YearInteger, const int& MonteCarloInt, const ints& MuonCharge, const floats& MuonPt,
					   const floats& MuonEta, const floats& MuonPhi, const ints& Muon_genPartIdx, const ints& Muon_nTrackerLayers){

	//std::cout << "print 44" << std::endl;

	std::string RoccoTextFile;

	switch(YearInteger){

		case 2016: RoccoTextFile = "./ScaleFactors/LeptonEnergyCorrections/RochesterCorrections/roccor.Run2.v3/RoccoR2016.txt"; break;
		case 2017: RoccoTextFile = "./ScaleFactors/LeptonEnergyCorrections/RochesterCorrections/roccor.Run2.v3/RoccoR2017.txt"; break;
		case 2018: RoccoTextFile = "./ScaleFactors/LeptonEnergyCorrections/RochesterCorrections/roccor.Run2.v3/RoccoR2018.txt"; break;
		default: std::cout << "Error for rochester corrections: choose a year out of 2016, 2017 or 2018." << std::endl; break;
	}

	RoccoR rc{RoccoTextFile};
	doubles RochCorrVec{};

	ints s(MuonPt.size(), 0); //s is error set (default is 0)
	ints m(MuonPt.size(), 0); //m is error member (default is 0, ranges from 0 to nmembers-1)
	doubles u{gRandom->Rndm(), gRandom->Rndm()}; //u is a random number distributed uniformly between 0 and 1 (gRandom->Rndm());


	for(unsigned int i = 0; i < MuonPt.size(); i++){

		//scale factors for momentum of each muon:
		double RochCorrSF;
		double mcSF;


		if(MonteCarloInt == 1){ 

			if(Muon_genPartIdx.size() > 0 && Muon_nTrackerLayers.size() > 0){

				if(mcSF > 0){
	
					RochCorrSF = rc.kSpreadMC(MuonCharge.at(i), MuonPt.at(i), MuonEta.at(i), MuonPhi.at(i), Muon_genPartIdx.at(0), s.at(i), m.at(i)); //(recommended), MC scale and resolution correction when matched gen muon is available
				}
				else{
					RochCorrSF = rc.kSmearMC(MuonCharge.at(i), MuonPt.at(i), MuonEta.at(i), MuonPhi.at(i), Muon_nTrackerLayers.at(0), u.at(i), s.at(i), m.at(i)); //MC scale and extra smearing when matched gen muon is not available

				}
			}
			else{RochCorrSF = 1.0;}

		}
		else{RochCorrSF = rc.kScaleDT(MuonCharge.at(i), MuonPt.at(i), MuonEta.at(i), MuonPhi.at(i), s.at(i), m.at(i));} //data

	
		RochCorrVec.push_back(RochCorrSF);


	}

	
	return RochCorrVec;


  }};


  auto RochCorrVec_Function{[&MCInt, &YearInt, &RochesterCorrections_testscript2](const ints& MuonCharge,      const floats& MuonPt,           const floats& MuonEta, 
									          const floats& MuonPhi,       const ints& Muon_genPartIdx,    const ints& Muon_nTrackerLayers){

  	//std::cout << "print 45" << std::endl;

  	floats CorrectionFactor = RochesterCorrections_testscript2(YearInt, MCInt, MuonCharge, MuonPt, MuonEta, MuonPhi, Muon_genPartIdx, Muon_nTrackerLayers);
  	return CorrectionFactor;

  }};


  auto RochCorrVec_Function_data{[&YearInt, &RochesterCorrections_testscript2, &MCInt](const ints& MuonCharge, const floats& MuonPt,        const floats& MuonEta,
						   			               const floats& MuonPhi,  const ints& DummyColumnInts, const ints& Muon_nTrackerLayers){

  	//std::cout << "print 46" << std::endl;

  	floats CorrectionFactor = RochesterCorrections_testscript2(YearInt, MCInt, MuonCharge, MuonPt, MuonEta, MuonPhi, DummyColumnInts, Muon_nTrackerLayers);
  	return CorrectionFactor;

  }};

  auto RochCorrMuon4Mo{[&ChannelInt](const TLorentzVector& Muon4Mo, const floats& RochCorrVec){

  	//std::cout << "print 47" << std::endl;

  	TLorentzVector NewVec{};

	double NewVecMass = Muon4Mo.M() * RochCorrVec.at(0);
        double NewVecPt = Muon4Mo.Pt() * RochCorrVec.at(0);
  	double NewVecPhi = Muon4Mo.Phi() * RochCorrVec.at(0);
  	double NewVecEta = Muon4Mo.Eta() * RochCorrVec.at(0);

        switch(ChannelInt){
        
                case 2: NewVec.SetPtEtaPhiM(NewVecPt, NewVecEta, NewVecPhi, NewVecMass); break;
		default: NewVec.SetPtEtaPhiM(Muon4Mo.Pt(), Muon4Mo.Eta(), Muon4Mo.Phi(), Muon4Mo.M()); break;

	}

  	return NewVec;

  }};

  auto TLorentzVector_float{[](const int& VariableOption, const TLorentzVector& object){
  
  	//std::cout << "print 48" << std::endl;

  	floats vec{};

	switch(VariableOption){
		case 1: vec.push_back(object.Pt()); break;
		case 2: vec.push_back(object.Phi()); break;
		case 3: vec.push_back(object.Eta()); break;
		case 4: vec.push_back(object.M()); break;
		default: break;
	}

  	return vec;

  }};

  auto TLorentzVector_float_pt{[&TLorentzVector_float](const TLorentzVector& object){
	//std::cout << "print 154" << std::endl; 
	return TLorentzVector_float(1, object);
  }};

  auto TLorentzVector_float_phi{[&TLorentzVector_float](const TLorentzVector& object){
	//std::cout << "print 155" << std::endl; 
	return TLorentzVector_float(2, object);
  }};

  auto TLorentzVector_float_eta{[&TLorentzVector_float](const TLorentzVector& object){
	//std::cout << "print 156" << std::endl; 
	return TLorentzVector_float(3, object);
  }};
  
  auto TLorentzVector_float_mass{[&TLorentzVector_float](const TLorentzVector& object){
	//std::cout << "print 157" << std::endl; 
	return TLorentzVector_float(4, object);
  }};

  auto z_mass_cut{[](const float& z_mass) {

  	//std::cout << "print 49" << std::endl;

  	return abs(z_mass - Z_MASS) < Z_MASS_CUT;

  }};


  auto RowReader2{[&FileNameJetSmear, &YearInt](const int& LineSpecified, const bool& sigmaJER, const bool& SF, const bool& up, 
					        const bool& down, const floats& Jet_eta, const floats& Jet_rho, const floats& Jet_pt) { 

  	//std::cout << "print 50" << std::endl;

  	float Col1, Col2, Col3, Col4, Col5, Col6, Col7, Col8, Col9, Col10, Col11;
 
	switch(YearInt){

		case 2016: if(sigmaJER == true && SF == false && up == false && down == false){
				FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2016/Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt";
			   }
  			   else if(sigmaJER == false && SF == true && up == false && down == false){
				FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2016/Summer16_25nsV1_MC_SF_AK4PFchs.txt";
			   }
  			   else if(sigmaJER == false && SF == false && up == true && down == false){
				FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2016/Summer16_25nsV1_MC_SF_AK4PFchs.txt";
			   }
  			   else if(sigmaJER == false && SF == false && up == false && down == true){
				FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2016/Summer16_25nsV1_MC_SF_AK4PFchs.txt";
			   }
  			   else{std::cout << "Please enter an appropriate file name" << std::endl;}

			   break;
 
  		case 2017: if(sigmaJER == true && SF == false && up == false && down == false){
				FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2017/Fall17_V3_MC_PtResolution_AK4PFchs.txt";
			   }
        		   else if(sigmaJER == false && SF == true && up == false && down == false){
				FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2017/Fall17_V3_MC_SF_AK4PFchs.txt";
			   }
        		   else if(sigmaJER == false && SF == false && up == true && down == false){
				FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2017/Fall17_V3_MC_SF_AK4PFchs.txt";
			   }
        		   else if(sigmaJER == false && SF == false && up == false && down == true){
				FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2017/Fall17_V3_MC_SF_AK4PFchs.txt";
			   }
        		   else{std::cout << "Please enter an appropriate file name" << std::endl;}

			   break;

		case 2018: if(sigmaJER == true && SF == false && up == false && down == false){
				FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2018/Autumn18_V1_MC_PtResolution_AK4PFchs.txt";
			   }
  			   else if(sigmaJER == false && SF == true && up == false && down == false){
				FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2018/Autumn18_V1_MC_SF_AK4PFchs.txt";
			   }
  			   else if(sigmaJER == false && SF == false && up == true && down == false){
				FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2018/Autumn18_V1_MC_SF_AK4PFchs.txt";
			   }
  			   else if(sigmaJER == false && SF == false && up == false && down == true){
				FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2018/Autumn18_V1_MC_SF_AK4PFchs.txt";
			   }
  			   else{std::cout << "Please enter an appropriate file name" << std::endl;}

			   break;

  		default: std::cout << "The year can only be 2016, 2017 or 2018" << std::endl; break;

	}

  	std::ifstream file;
  	file.open(FileNameJetSmear);

  	if (file.good())
  	{
    		std::string str = "";

    		int line_num = 0;
	
	 	while(getline(file, str) && line_num != LineSpecified){
			++line_num;
		}
		if(line_num == LineSpecified){
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
 

  	floats AnswerVec{};
 
  	for(unsigned int i = 0; i < Jet_pt.size(); i++){

		if(  (Jet_eta.at(i) > abs(Col1) && Jet_eta.at(i) < abs(Col2)) && 
	     	     (Jet_rho.at(0) > abs(Col3) && Jet_rho.at(0) < abs(Col4)) &&
	             (Jet_pt.at(i) > abs(Col6) && Jet_pt.at(i) < abs(Col7) ) ){

  			if(sigmaJER == true && SF == false && up == false && down == false){

        			float answer = sqrt( Col8*abs(Col8) / (Jet_pt.at(i)*Jet_pt.at(i))+Col9*Col9*pow(Jet_pt.at(i),Col11)+Col10*Col10 );
				AnswerVec.push_back(answer);

			}
 	 		else if(sigmaJER == false && SF == true && up == false && down == false){
        
        			AnswerVec.push_back(Col4);
  			
			}
  			else if(sigmaJER == false && SF == false && up == true && down == false){

        			float UpValue = Col6 - Col4;
        			AnswerVec.push_back(UpValue);

			}
  			else if(sigmaJER == false && SF == false && up == false && down == true){
	
        			float DownValue = Col4 - Col5;
        			AnswerVec.push_back(DownValue);

  			}
  			else{std::cout << "bools cannot be all true or all false" << std::endl; 
			     std::cout << "sigmaJER = " << sigmaJER << std::endl; 
			     std::cout << "SF = " << SF << std::endl; 
			     std::cout << "up = " << up << std::endl; 
			     std::cout << "down = " << down << std::endl;} 

		}
		else{float zero = 0.0; AnswerVec.push_back(zero);}

  
   	} //end of for loop


   	return AnswerVec;


  }}; 



  auto linecounter{[&FileNameJetSmear, &YearInt](const bool& sigmaJER, const bool& SF, const bool& up, const bool& down){ 

  	//std::cout << "print 51" << std::endl;

   	int number_of_lines = 0;
   	std::string line;

	switch(YearInt){

		case 2016: if(sigmaJER == true && SF == false && up == false && down == false){
				FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2016/Summer16_25nsV1_MC_PtResolution_AK4PFchs.txt";
			   }
        		   else if(sigmaJER == false && SF == true && up == false && down == false){
				FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2016/Summer16_25nsV1_MC_SF_AK4PFchs.txt";
			   }
        		   else if(sigmaJER == false && SF == false && up == true && down == false){
				FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2016/Summer16_25nsV1_MC_SF_AK4PFchs.txt";
			   }
        		   else if(sigmaJER == false && SF == false && up == false && down == true){
				FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2016/Summer16_25nsV1_MC_SF_AK4PFchs.txt";
			   }
        		   else{std::cout << "Please enter an appropriate file name" << std::endl;}

			   break;

		case 2017: if(sigmaJER == true && SF == false && up == false && down == false){
				FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2017/Fall17_V3_MC_PtResolution_AK4PFchs.txt";
			   }
        		   else if(sigmaJER == false && SF == true && up == false && down == false){
				FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2017/Fall17_V3_MC_SF_AK4PFchs.txt";
			   }
        		   else if(sigmaJER == false && SF == false && up == true && down == false){
				FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2017/Fall17_V3_MC_SF_AK4PFchs.txt";
			   }
        		   else if(sigmaJER == false && SF == false && up == false && down == true){
				FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2017/Fall17_V3_MC_SF_AK4PFchs.txt";
			   }
        		   else{std::cout << "Please enter an appropriate file name" << std::endl;}

			   break;

		case 2018: if(sigmaJER == true && SF == false && up == false && down == false){
				FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2018/Autumn18_V1_MC_PtResolution_AK4PFchs.txt";
			   }
        		   else if(sigmaJER == false && SF == true && up == false && down == false){
				FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2018/Autumn18_V1_MC_SF_AK4PFchs.txt";
			   }
        		   else if(sigmaJER == false && SF == false && up == true && down == false){
				FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2018/Autumn18_V1_MC_SF_AK4PFchs.txt";
			   }
        		  else if(sigmaJER == false && SF == false && up == false && down == true){
				FileNameJetSmear = "./ScaleFactors/JECs/JetSmearing/2018/Autumn18_V1_MC_SF_AK4PFchs.txt";
			   }
        		  else{std::cout << "Please enter an appropriate file name" << std::endl;}

			  break;

  		default: std::cout << "The year can only be 2016, 2017 or 2018" << std::endl; break;

	} 


   	std::ifstream myfile(FileNameJetSmear);

   	while (getline(myfile, line))
        	++number_of_lines;
    		return number_of_lines;

  }};




  auto RowReader3{[&RowReader2, &linecounter](const bool& SigmaJER,  const bool& JetSmearScaleFactor, const bool& Up, const bool& Down,
					      const floats& Jet_eta, const floats& Jet_rho, const floats& Jet_pt){

  	//std::cout << "print 52" << std::endl;

  	int k;

  	for(int i = 0; i < linecounter(SigmaJER, JetSmearScaleFactor, Up, Down) + 1; i++){

		std::string quantity; 

		if(SigmaJER == true && JetSmearScaleFactor == false && Up == false && Down == false){quantity = "sigma JER";}
   		else if(SigmaJER == false && JetSmearScaleFactor == true && Up == false && Down == false){quantity = "SF";}
   		else if(SigmaJER == false && JetSmearScaleFactor == false && Up == true && Down == false){quantity = "SF (up variation)";}
   		else if(SigmaJER == false && JetSmearScaleFactor == false && Up == false && Down == true){quantity = "SF (down variation)";}
   		else{std::cout << "Please enter an appropriate file name" << std::endl;}



		bool check = any_of(RowReader2(i, SigmaJER, JetSmearScaleFactor, Up, Down, Jet_eta, Jet_rho, Jet_pt).begin(),
				    RowReader2(i, SigmaJER, JetSmearScaleFactor, Up, Down, Jet_eta, Jet_rho, Jet_pt).end(),
				    [](float j){return j != 0;});


		if(check == 1){k = i; break;}
		else{continue;}


  	}

  	float factor;

  	for(long unsigned int i = 0; i < RowReader2(k, SigmaJER, JetSmearScaleFactor, Up, Down, Jet_eta, Jet_rho, Jet_pt).size(); i++){

		if(RowReader2(k, SigmaJER, JetSmearScaleFactor, Up, Down, Jet_eta, Jet_rho, Jet_pt).at(i) != 0){factor = RowReader2(k, SigmaJER, JetSmearScaleFactor, Up, Down, Jet_eta, Jet_rho, Jet_pt).at(i);}
		else{continue;}

  	}

  	return factor;

  }};

 auto sigma_JER{[&RowReader3](const floats& Jet_eta, const floats& Jet_rho,const floats& Jet_pt){

  //std::cout << "print 53" << std::endl;

  bool SigmaJER = true;
  bool JetSmearScaleFactor = false;
  bool Up = false;
  bool Down = false;

  return RowReader3(SigmaJER, JetSmearScaleFactor, Up, Down, Jet_eta, Jet_rho, Jet_pt);

}};


auto sigma_JER_up{[&RowReader3](const floats& Jet_eta, const floats& Jet_rho,const floats& Jet_pt){

  //std::cout << "print 54" << std::endl;

  bool SigmaJER = false;
  bool JetSmearScaleFactor = false;
  bool Up = true;
  bool Down = false;

  return RowReader3(SigmaJER, JetSmearScaleFactor, Up, Down, Jet_eta, Jet_rho, Jet_pt);

}};


auto sigma_JER_down{[&RowReader3](const floats& Jet_eta, const floats& Jet_rho,const floats& Jet_pt){

  //std::cout << "print 55" << std::endl;

  bool SigmaJER = false;
  bool JetSmearScaleFactor = false;
  bool Up = false;
  bool Down = true;
  
  return RowReader3(SigmaJER, JetSmearScaleFactor, Up, Down, Jet_eta, Jet_rho, Jet_pt);

}};
 


  auto SJER_Nominal_Function{[&RowReader3](const floats& Jet_eta, const floats& Jet_rho, const floats& Jet_pt){

  	//std::cout << "print 53" << std::endl;

  	bool SigmaJER = false;
  	bool JetSmearScaleFactor = true;
  	bool Up = false;
  	bool Down = false;
 
  	return RowReader3(SigmaJER, JetSmearScaleFactor, Up, Down, Jet_eta, Jet_rho, Jet_pt);

  }};

  auto SJER_Up_Function{[&RowReader3](const floats& Jet_eta, const floats& Jet_rho, const floats& Jet_pt){

  	//std::cout << "print 54" << std::endl;

  	bool SigmaJER = false;
  	bool JetSmearScaleFactor = false;
  	bool Up = true;
  	bool Down = false;
  
  	return RowReader3(SigmaJER, JetSmearScaleFactor, Up, Down, Jet_eta, Jet_rho, Jet_pt);

  }};

  auto SJER_Down_Function{[&RowReader3](const floats& Jet_eta, const floats& Jet_rho, const floats& Jet_pt){

  	//std::cout << "print 55" << std::endl;

  	bool SigmaJER = false;
  	bool JetSmearScaleFactor = false;
  	bool Up = false;
  	bool Down = true;
  
  	return RowReader3(SigmaJER, JetSmearScaleFactor, Up, Down, Jet_eta, Jet_rho, Jet_pt);

  }};


  auto MaxComparison{[](const float& sJER_nominal){

  	//std::cout << "print 56" << std::endl;

 	float MaximumFloats = sqrt(sJER_nominal*sJER_nominal - 1);

 	if(MaximumFloats > 0){return MaximumFloats;}
 	else{float zero = 0.0; return zero;}

  }};


  auto JetSmearingFunction_HybridMethod{[&MaxComparison](const floats& pT, const floats& eta, const floats& phi, const floats& pT_ptcl, const floats& eta_ptcl, 
							 const floats& phi_ptcl, const float& sJER_nominal, const float& sigma_JER_input, const ints& Jet_genJetIdx){

  	//std::cout << "print 57" << std::endl;

  	floats cJER_vec{};

  	for(long unsigned int i = 0; i < pT.size(); i++){

		float cJER_Scaling;
		float N = gRandom->Gaus(0, sigma_JER_input);
        	float cJER_Stochastic = 1.0 + ( N * MaxComparison(sJER_nominal) );

  		if(Jet_genJetIdx.at(i) != -1){

			long unsigned int j = Jet_genJetIdx.at(i);

				if( j < pT_ptcl.size() ){

					double dphi = phi.at(i) - phi_ptcl.at(j);
        				double deta = eta.at(i) - eta_ptcl.at(j);
        				double deltaR = sqrt( pow(dphi, 2) + pow(deta, 2) );
        				const double RCone = 0.4;

 					if( (abs(pT.at(i) - pT_ptcl.at(j)) < 3 * sigma_JER_input * pT.at(i)) && (deltaR == RCone / 2) ){

						cJER_Scaling = 1 + ( (sJER_nominal - 1) * ( (pT.at(i) - pT_ptcl.at(j)) / pT.at(i) ) );
						cJER_vec.push_back(cJER_Scaling);
		
					}
					else{cJER_vec.push_back(cJER_Stochastic);}

				}
				else{cJER_vec.push_back(cJER_Stochastic);}

  		}
  		else{cJER_vec.push_back(cJER_Stochastic);}


  	}

  	return cJER_vec;
  
  }};



  auto ApplyCJER{[](const floats& JetPt, const floats& JetEta, const floats& JetPhi, const floats& JetMass, const floats& cJER, const unsigned int& nJet){

  	//std::cout << "print 58" << std::endl;

  	std::vector<TLorentzVector> OutputVec{};

  	for(unsigned int i = 0; i < nJet; i++){

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


  auto GetSmearedJetPt{[](std::vector<TLorentzVector> SmearedJet4Momentum, const floats& JetPt){

  	//std::cout << "print 59" << std::endl;

  	floats NewPtVec = {};

 	for(long unsigned int i = 0; i < JetPt.size(); i++){

        	float NewPt = (SmearedJet4Momentum.at(i)).Pt();
 		NewPtVec.push_back(NewPt);

	}

 	return NewPtVec;

  }};


  auto GetSmearedJetPhi{[](std::vector<TLorentzVector> SmearedJet4Momentum, const floats& JetPhi){

  	//std::cout << "print 60" << std::endl;

 	floats NewPhiVec{};

 	for(long unsigned int i = 0; i < JetPhi.size(); i++){

		float NewPhi = (SmearedJet4Momentum.at(i)).Phi();
        	NewPhiVec.push_back(NewPhi);

 	}
 
	return NewPhiVec;

  }};

 
  auto GetSmearedJetEta{[](std::vector<TLorentzVector> SmearedJet4Momentum, const floats& JetEta){

  	//std::cout << "print 61" << std::endl;

 	floats NewEtaVec = {};

 	for(long unsigned int i = 0; i < JetEta.size(); i++){

        	float NewEta = (SmearedJet4Momentum.at(i)).Eta();
        	NewEtaVec.push_back(NewEta);

 	}

 	return NewEtaVec;

  }};


  auto GetSmearedJetMass{[](std::vector<TLorentzVector> SmearedJet4Momentum, const floats& JetMass){

  	//std::cout << "print 62" << std::endl;

 	floats NewMassVec = {};

 	for(long unsigned int i = 0; i < JetMass.size(); i++){

        	float NewMass = (SmearedJet4Momentum.at(i)).M();	
		NewMassVec.push_back(NewMass);

 	}

 	return NewMassVec;

  }};

  
  auto SumSquared2LeadingJets_pT{[](const float& LeadingJetPt, const float& SubleadingJetPt){

  	//std::cout << "print 63" << std::endl;

  	double SumSquaredPt = pow(LeadingJetPt + SubleadingJetPt, 2);
  	return SumSquaredPt;

  }};


  auto JetSum{[](const float& LeadingJet, const float& SubleadingJet, const float& ThirdJet, const float& FourthJet){

  	//std::cout << "print 64" << std::endl;

  	float JetSumOutput = LeadingJet + SubleadingJet + ThirdJet + FourthJet;
  	return JetSumOutput;

  }};

  auto delta_phi{[](const float phi1, const float phi2){
  	return vdt::fast_atan2f(vdt::fast_sinf(phi1 - phi2), vdt::fast_cosf(phi1 - phi2));
  }};

  auto deltaR{[&delta_phi](const float eta1, const float phi1, const float eta2, const float phi2){
  	return std::sqrt(std::pow(eta1 - eta2, 2) + std::pow(delta_phi(phi1, phi2), 2));
  }};

  auto deltaRcheck_floats{[&deltaR](const floats& Object1_eta, const floats& Object1_phi, const floats& Object2_eta, const floats& Object2_phi) {

  	//std::cout << "print 65" << std::endl;

  	floats min_dRs{};

	 if(Object2_phi.size() > 1){

  		transform(Object1_eta.begin(), Object1_eta.end(), Object1_phi.begin(), std::back_inserter(min_dRs), [&](float Object1_eta, float Object1_phi) { return std::min(deltaR(Object1_eta, Object1_phi, Object2_eta.at(0), Object2_phi.at(0)), deltaR(Object1_eta, Object1_phi, Object2_eta.at(1), Object2_phi.at(1))); });

  	}
  	else{

		transform(Object1_eta.begin(), Object1_eta.end(), Object1_phi.begin(), std::back_inserter(min_dRs), [&](float Object1_eta, float Object1_phi) { return deltaR(Object1_eta, Object1_phi, Object2_eta.at(0), Object2_phi.at(0)); });

  	}	

  	return min_dRs;
 
  }};


  auto HT{[](const float& Pt){
 
  	//std::cout << "print 66" << std::endl;

  	float HTOutput = abs(Pt);
  	return HTOutput;

  }};

  auto TotJetHT{[](const float& LeadingJetHT, const float& SubleadingJetHT, const float& ThirdJetHT, const float& FourthJetHT){
  
  	//std::cout << "print 67" << std::endl;

  	float TotJetHTOutput = LeadingJetHT + SubleadingJetHT + ThirdJetHT + FourthJetHT;
  	return TotJetHTOutput;

  }};

  auto TotLepHT{[](const float& LeadingLeptonHT, const float& SubleadingLeptonHT){

  	//std::cout << "print 68" << std::endl;

  	float TotLepHTOutput = LeadingLeptonHT + SubleadingLeptonHT;
  	return TotLepHTOutput;

  }};

  
  auto TotHTOverTotpT{[](const float& TotHT, const float& TotpT){

  	//std::cout << "print 69" << std::endl;

  	float TotHTOverTotpTOutput = TotHT / TotpT;
  	return TotHTOverTotpTOutput;

  }};


  auto LepSum{[](const float& LeadingLep, const float& SubleadingLep){

  	//std::cout << "print 70" << std::endl;

  	float LepSumOutput = LeadingLep + SubleadingLep;
  	return LepSumOutput;

  }};


  auto InvMass_AllJets{[](const float& LeadingJetPt,   const float& SubleadingJetPt,   const float& ThirdJetPt,   const float& FourthJetPt,
			  const float& LeadingJetEta,  const float& SubleadingJetEta,  const float& ThirdJetEta,  const float& FourthJetEta,
		    	  const float& LeadingJetPhi,  const float& SubleadingJetPhi,  const float& ThirdJetPhi,  const float& FourthJetPhi,
			  const float& LeadingJetMass, const float& SubleadingJetMass, const float& ThirdJetMass, const float& FourthJetMass){

  	//std::cout << "print 71" << std::endl;

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


  auto InvMass_3Jets{[](const float& LeadingJetPt,   const float& SubleadingJetPt,   const float& ThirdJetPt,
			const float& LeadingJetEta,  const float& SubleadingJetEta,  const float& ThirdJetEta,
			const float& LeadingJetPhi,  const float& SubleadingJetPhi,  const float& ThirdJetPhi,
			const float& LeadingJetMass, const float& SubleadingJetMass, const float& ThirdJetMass){
  
  	//std::cout << "print 72" << std::endl;

  	TLorentzVector Jet1 = {};
  	TLorentzVector Jet2 = {};
  	TLorentzVector Jet3 = {};

  	Jet1.SetPtEtaPhiM(LeadingJetPt, LeadingJetEta, LeadingJetPhi, LeadingJetMass);
  	Jet2.SetPtEtaPhiM(SubleadingJetPt, SubleadingJetEta, SubleadingJetPhi, SubleadingJetMass);
  	Jet3.SetPtEtaPhiM(ThirdJetPt, ThirdJetEta, ThirdJetPhi, ThirdJetMass);

  	float InvMass3Jets = (Jet1 + Jet2 + Jet3).M();
  	return InvMass3Jets;

  }};


  auto tight_jets_function{[&YearInt](const floats& Jet_pt_Selection, const floats& Jet_eta_Selection, const ints& Jet_jetId_Selection, const floats& dRJet_lep){

  	//std::cout << "print 73" << std::endl;

  	int JetId;

	switch(YearInt){
  		case 2016: JetId = 1; break; //1 is loose 
  		case 2017: JetId = 2; break; //2 is tight
		case 2018: JetId = 2; break; //2 is tight
  		default: std::cout << "Choose a year out of 2016, 2017 or 2018" << std::endl; break;
	}

  	return Jet_pt_Selection > 30 && Jet_eta_Selection < 4.7 && Jet_jetId_Selection >= JetId && dRJet_lep > 0.4;

  }};


  auto jet_selection_function{[](const ints& tight_jets) {

  	//std::cout << "print 74" << std::endl;

  	auto njet{count_if(tight_jets.begin(), tight_jets.end(), [](int i) { return i; })};
  	return njet >= 4 && njet <= 6;

  }};

  auto bjet_id{[](const ints& tight_jets, const floats& btags, const floats& etas) {
     
        //std::cout << "print 75" << std::endl;
	return tight_jets && (btags > 0.8838f) && (etas < MaxTrackerEta);
  
  }};


  auto numberofbjets{[](const ints& bjets) {

	//std::cout << "print 76" << std::endl;
        const auto nbjet{std::count_if(bjets.begin(), bjets.end(), [](int i) { return i; })};
        return nbjet;

  }};

  auto BTAGEFF_bjet_id_WP{[](const ints& tight_jets, const floats& btags, const floats& etas, const ints& Jet_partonFlavour) {

	//std::cout << "print 77" << std::endl;
	return abs(Jet_partonFlavour) == 5 && btags > 0.8838f && abs(etas) < MaxTrackerEta;
	
  }};


  auto BTAGEFF_charm_id_WP{[](const ints& tight_jets, const floats& btags, const floats& etas, const ints& Jet_partonFlavour) {

	//std::cout << "print 78" << std::endl;
        return abs(Jet_partonFlavour) == 4 && btags > 0.8838f && abs(etas) < MaxTrackerEta;

  }};



  auto BTAGEFF_lightjets_id_WP{[](const ints& tight_jets, const floats& btags, const floats& etas, const ints& Jet_partonFlavour) {
        
      	//std::cout << "print 79" << std::endl;
      	return abs(Jet_partonFlavour) > 0 && abs(Jet_partonFlavour) < 4 && btags > 0.8838f && abs(etas) < MaxTrackerEta;

  }};



  auto BTAGEFF_gluon_id_WP{[](const ints& tight_jets, const floats& btags, const floats& etas, const ints& Jet_partonFlavour) {
                
      	//std::cout << "print 80" << std::endl;
      	return abs(Jet_partonFlavour) == 21 && btags > 0.8838f && abs(etas) < MaxTrackerEta;

  }};


  auto BTAGEFF_nonbjet_id_WP{[](const ints& tight_jets, const floats& btags, const floats& etas, const ints& Jet_partonFlavour){

  	//std::cout << "print 81" << std::endl;
    	return abs(Jet_partonFlavour) != 5 && btags > 0.8838f && abs(etas) < MaxTrackerEta;

  }};


  auto BTAGEFF_bjet_id{[](const ints& tight_jets, const floats& etas, const ints& Jet_partonFlavour) {

  	//std::cout << "print 82" << std::endl;
	return abs(Jet_partonFlavour) == 5 && abs(etas) < MaxTrackerEta;

  }};


  auto BTAGEFF_charm_id{[](const ints& tight_jets, const floats& etas, const ints& Jet_partonFlavour) {

	//std::cout << "print 83" << std::endl;
	return abs(Jet_partonFlavour) == 4 && abs(etas) < MaxTrackerEta;

  }};


  auto BTAGEFF_lightjets_id{[](const ints& tight_jets, const floats& etas, const ints& Jet_partonFlavour) {

	//std::cout << "print 84" << std::endl;
	return abs(Jet_partonFlavour) > 0 && abs(Jet_partonFlavour) < 4 && abs(etas) < MaxTrackerEta;

  }};


  auto BTAGEFF_gluon_id{[](const ints& tight_jets, const floats& etas, const ints& Jet_partonFlavour) {

	//std::cout << "print 85" << std::endl;
        return abs(Jet_partonFlavour) == 21 && abs(etas) < MaxTrackerEta;

  }};


 auto BTAGEFF_nonbjet_id{[](const ints& tight_jets, const floats& etas, const ints& Jet_partonFlavour){
	
	//std::cout << "print 86" << std::endl;
	return abs(Jet_partonFlavour) != 5 && abs(etas) < MaxTrackerEta;

  }};

  auto bjet_cut{[](const ints& bjets) {

        //std::cout << "print 87" << std::endl;

        const auto nbjet{std::count_if(bjets.begin(), bjets.end(), [](int i) { return i; })};
        return nbjet >= 1 && nbjet <= 2;

  }};

  auto find_lead_mask{[](const ints& mask, const floats& vals) {
  
  	//std::cout << "print 88" << std::endl;

  	const auto masked_vals{mask * vals};
  	const auto max_idx{boost::numeric_cast<size_t>(std::distance(masked_vals.begin(), max_element(masked_vals.begin(), masked_vals.end())))};
  	ints lead_mask(masked_vals.size(), 0); 
  	lead_mask.at(max_idx) = 1;
  	return lead_mask;

  }};

  auto WPair{[](const floats& pts, const floats& etas, const floats& phis, const floats& ms, const ints& tight_jets, const ints& lead_bjet) {

        //std::cout << "print 89" << std::endl;

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

        auto jet1{TLorentzVector{}};
        auto jet2{TLorentzVector{}};

        jet1.SetPtEtaPhiM(pts.at(jet_index_1), etas.at(jet_index_1), phis.at(jet_index_1), ms.at(jet_index_1));
        jet2.SetPtEtaPhiM(pts.at(jet_index_2), etas.at(jet_index_2), phis.at(jet_index_2), ms.at(jet_index_2));

	return w_pair;

  }};


  auto WPairJet1{[](const floats& pts, const floats& etas, const floats& phis, const floats& ms, const ints& tight_jets, const ints& lead_bjet) {

	//std::cout << "print 90" << std::endl;

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

        auto jet1{TLorentzVector{}};
        auto jet2{TLorentzVector{}};

        jet1.SetPtEtaPhiM(pts.at(jet_index_1), etas.at(jet_index_1), phis.at(jet_index_1), ms.at(jet_index_1));
        jet2.SetPtEtaPhiM(pts.at(jet_index_2), etas.at(jet_index_2), phis.at(jet_index_2), ms.at(jet_index_2));

        return jet1;
 
  }};


  auto WPairJet2{[](const floats& pts, const floats& etas, const floats& phis, const floats& ms, const ints& tight_jets, const ints& lead_bjet) {

	//std::cout << "print 91" << std::endl;

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

        auto jet1{TLorentzVector{}};
        auto jet2{TLorentzVector{}};

        jet1.SetPtEtaPhiM(pts.at(jet_index_1), etas.at(jet_index_1), phis.at(jet_index_1), ms.at(jet_index_1));
        jet2.SetPtEtaPhiM(pts.at(jet_index_2), etas.at(jet_index_2), phis.at(jet_index_2), ms.at(jet_index_2));

        return jet2;   
 
  }};


  auto deltaRcheck_W_function{[](const doubles& Object1_phi_Selection, const doubles& Object1_eta_Selection,
				 const doubles& Object2_eta_Selection, const doubles& Object2_phi_Selection){

  	//std::cout << "print 93" << std::endl;

  	doubles dR = sqrt(pow(Object1_eta_Selection - Object2_eta_Selection, 2) + pow(Object1_phi_Selection - Object2_phi_Selection, 2));
  	return dR;

  }};


  auto DeltaPhi_function2{[](const doubles& Object1_phi_Selection, const doubles& Object2_phi_Selection){

  	//std::cout << "print 94" << std::endl;

  	doubles dPhi = abs(Object1_phi_Selection - Object2_phi_Selection);
  	return dPhi;

  }};


  auto deltaRcheck_W_function2{[](const doubles& Object1_phi_Selection, const doubles& Object1_eta_Selection,
				  const float& Object2_eta_Selection,   const float& Object2_phi_Selection){
 
  	//std::cout << "print 95" << std::endl;

  	doubles dR = sqrt(pow(Object1_eta_Selection - Object2_eta_Selection, 2) + pow(Object1_phi_Selection - Object2_phi_Selection, 2));
  	return dR;

  }};

  auto DeltaPhi_doublesandfloat{[](const doubles& Object1_phi, const float& Object2_phi){

  	//std::cout << "print 96" << std::endl;

  	doubles dPhi = abs(Object1_phi - Object2_phi);
  	return dPhi;

  }};
  
  auto HT_double{[](const doubles& Pt){

  	//std::cout << "print 97" << std::endl;

  	doubles HT_Output = abs(Pt);
  	return HT_Output;

  }};

  auto RecoWHT{[](const floats& RecoWPt){

  	//std::cout << "print 98" << std::endl;

  	floats RecoWHTOutput = abs(RecoWPt);
  	return RecoWHTOutput;

  }};

  auto TransverseWMass{[](const double& dPhi_j1j2, const doubles& WPairJet1Pt, const doubles& WPairJet2Pt){

  	//std::cout << "print 99" << std::endl;

  	doubles mtW = sqrt(2 * WPairJet1Pt * WPairJet2Pt * (1 - cos(dPhi_j1j2)) );
  	return mtW;

  }};


  auto w_mass_cut{[&ZPlusJetsCRInt](const float& w_mass, const float& MET_sumEt) {
	
  	//std::cout << "print 100" << std::endl;
      
  	switch(ZPlusJetsCRInt){
      
    		case 0: return ( abs(w_mass - W_MASS) < W_MASS_CUT ); break;
    		case 1: return abs(w_mass - W_MASS) > W_MASS_CUT && (MET_sumEt < 50); break;
      
        }

  }};


  auto WLorentzVector{[](const floats& w_pair_pt, const floats& w_pair_eta, const floats& w_pair_phi, const float& w_mass, const ints& w_reco_jets){

  	//std::cout << "print 101" << std::endl;

  	const auto nRecoWBosons{std::count_if(w_reco_jets.begin(), w_reco_jets.end(), [](int i) { return i; })};
  	auto RecoW = TLorentzVector{};
  
  	for(int i = 0; i < nRecoWBosons; i++){

	  	auto Vec = TLorentzVector{};
	  	Vec.SetPtEtaPhiM(w_pair_pt.at(i), w_pair_eta.at(i), w_pair_phi.at(i), w_mass);
	  	RecoW += Vec;

  	}

  	return RecoW;

  }};


  auto bjet_variable{[](const floats& Jet_variable, const unsigned int& nJet, const ints& lead_bjet){

  	//std::cout << "print 102" << std::endl;

  	floats vec{};

  	for(unsigned int i = 0; i < nJet; i++){
        	if(lead_bjet.at(i) == 1){ 
			vec.push_back(Jet_variable.at(i));
		}

  	}

  	return vec;

  }};


  auto BLorentzVector{[](const floats& bjet_pt, const floats& bjet_eta, const floats& bjet_phi, const floats& bjet_mass){

  	//std::cout << "print 103" << std::endl;

  	auto BJets = TLorentzVector{};

  	for(long unsigned int i = 0; i < bjet_pt.size(); i++){

		auto Vec = TLorentzVector{};
		Vec.SetPtEtaPhiM(bjet_pt.at(i), bjet_eta.at(i), bjet_phi.at(i), bjet_mass.at(i));
		BJets += Vec;

  	}

  	return BJets;

  }};


  auto top_reconstruction_function{[](const floats& bjets_pt,  const floats& bjets_eta,  const floats& bjets_phi,  const floats& bjets_mass,
				      const floats& w_pair_pt, const floats& w_pair_eta, const floats& w_pair_phi, const float& w_mass ){

  	//std::cout << "print 104" << std::endl;

  	auto reco_top = TLorentzVector{}; 
  	auto BJets = TLorentzVector{};
  	auto RecoW = TLorentzVector{};

  	double top_reco_mass = std::numeric_limits<double>::infinity();
  	size_t index_1{std::numeric_limits<size_t>::max()};
  	const size_t num{w_pair_pt.size()};

  	for(unsigned int i = 0; i < num; i++){

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

  
  auto deltaRcheck_Top_function{[](const doubles& Object1_phi_Selection, const doubles& Object1_eta_Selection,
				   const float& Object2_eta_Selection,   const float& Object2_phi_Selection){

  	//std::cout << "print 106" << std::endl;

  	doubles dR = sqrt(pow(Object1_eta_Selection - Object2_eta_Selection, 2) + pow(Object1_phi_Selection - Object2_phi_Selection, 2));
  	return dR;

  }};

  auto deltaRcheck_WTop_function{[](const floats& Object1_phi_Selection,  const floats& Object1_eta_Selection,
				    const doubles& Object2_eta_Selection, const doubles& Object2_phi_Selection){

  	//std::cout << "print 107" << std::endl;

  	doubles dR_vec{};

  	for(long unsigned int i = 0; i < Object1_phi_Selection.size(); i++){

  		double dR = sqrt(pow(Object1_eta_Selection.at(i) - Object2_eta_Selection.at(0), 2) + pow(Object1_phi_Selection.at(i) - Object2_phi_Selection.at(0), 2));
  		dR_vec.push_back(dR);

  	}

  	return dR_vec;

  }};

  
  auto MinDeltaR{[](const unsigned int& nJet, const doubles& RecoZPhi, const doubles& RecoZEta, const floats& Jet_Phi_Selection, const floats& Jet_eta_Selection){

  	//std::cout << "print 108" << std::endl;

    	doubles output_vec;
	double Output;  

    	for(unsigned int i = 0; i < nJet; i++){

		if(RecoZEta.size() > 1){
    			double DeltaR = sqrt(pow(RecoZPhi.at(i) - Jet_Phi_Selection.at(i), 2) + pow(RecoZEta.at(i) - Jet_eta_Selection.at(i), 2));
    			double DeltaR2 = sqrt(pow(RecoZPhi.at(i+1) - Jet_Phi_Selection.at(i+1), 2) + pow(RecoZEta.at(i+1) - Jet_eta_Selection.at(i+1), 2));

    			Output = (DeltaR2 < DeltaR) ? DeltaR2 : DeltaR;  
		
		}
		else{Output = sqrt(pow(RecoZPhi.at(0) - Jet_Phi_Selection.at(i), 2) + pow(RecoZEta.at(0) - Jet_eta_Selection.at(i), 2));}

		output_vec.push_back(Output);

    	}

    	return output_vec;

  }};


  auto MinDeltaPhi{[](const unsigned int& nJet, const doubles& RecoZPhi, const floats& Jet_Phi_Selection){

  	//std::cout << "print 109" << std::endl;

  	double output;
  	doubles output_vec{};

  	for(unsigned int i = 0; i < nJet; i++){

		if(RecoZPhi.size() > 1){
    			
			double DeltaPhi = std::abs(RecoZPhi.at(i) - Jet_Phi_Selection.at(i));
    			double DeltaPhi2 = std::abs(RecoZPhi.at(i+1) - Jet_Phi_Selection.at(i+1));

    			output = (DeltaPhi2 < DeltaPhi) ? DeltaPhi2 : DeltaPhi;

		}
		else{output = std::abs(RecoZPhi.at(0) - Jet_Phi_Selection.at(i));}

    		output_vec.push_back(output);

	}

  	return output_vec;

  }};


  auto dR_Lepton_LeadingBJet_Function{[](const floats& bjeteta, const float& LeptonEta, const floats& bjetphi, const float& LeptonPhi){

  	//std::cout << "print 110" << std::endl;

  	doubles DeltaR = sqrt(pow(LeptonPhi - bjetphi, 2) + pow(LeptonEta - bjeteta, 2));
  	return DeltaR;

  }};


  auto DeltaPhi_Lepton_BJet{[](const floats& Jet_phi_Selection, const float& LeptonPhi){

  	//std::cout << "print 111" << std::endl;

  	doubles DeltaPhi = abs(LeptonPhi - Jet_phi_Selection);
  	return DeltaPhi;

  }};


  auto MET_function{[](const floats& MET_input){

  	//std::cout << "print 112" << std::endl;
  	return MET_input;
  
  }};

  
  auto BJetOutputDiscriminantFunction{[](const float& JetPt, const floats& Jet_btagCSVV2, const ints& tight_jets, const floats& Jet_eta_Selection){

  	//std::cout << "print 112" << std::endl;
  	return JetPt && (Jet_btagCSVV2  > 0.8838) && tight_jets && (abs(Jet_eta_Selection) < MaxTrackerEta);

  }};


  auto DeltaPhi_function4{[](const floats& Object1_phi, const doubles& Object2_phi){

  	//std::cout << "print 113" << std::endl;

 	doubles dPhi_vec{};

 	for(long unsigned int i = 0; i < Object1_phi.size(); i++){

 		double dPhi = Object1_phi.at(i) - Object2_phi.at(0);
		dPhi_vec.push_back(dPhi);

 	}

 	return dPhi_vec;

  }};


  auto TotalVariable_System{[](const doubles& RecoZInput, const floats& RecoWInput, const doubles& TopInput, const float& TotLepInput, const float& TotJetInput){

  	//std::cout << "print 114" << std::endl;

  	doubles TotalSystemOutput = RecoZInput + RecoWInput.at(0) + TopInput + TotLepInput + TotJetInput;
  	return TotalSystemOutput;

  }};


  auto inv_mass_doubles{[](const doubles& pts, const doubles& etas, const doubles& phis, const doubles& ms){

	//std::cout << "print 115" << std::endl;

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

  }};

  auto UnweightedTopPt{[](const doubles& pts){

	//std::cout << "print 116" << std::endl;
        return pts;

  }};


  auto TopReweighting_topquark{[](const ints& GenPart_pdgId, const ints& GenPart_statusFlags, const floats& GenPart_pt){

  	//std::cout << "print 117" << std::endl;
	return GenPart_pdgId == 6 && GenPart_statusFlags == 13 && GenPart_pt > 0; 

  }};

  auto TopReweighting_antitopquark{[](const ints& GenPart_pdgId, const ints& GenPart_statusFlags, const floats& GenPart_pt){
		
	//std::cout << "print 118" << std::endl;
	return GenPart_pdgId == -6 && GenPart_statusFlags == 13 && GenPart_pt > 0; 

  }};


  auto TopReweighting_weight{[&ProcessInt](const ints& TopReweighting_topquark_input, const ints& TopReweighting_antitopquark_input){

	//std::cout << "print 119" << std::endl;

	doubles SF_top = exp(-0.0615-(0.00005* TopReweighting_topquark_input) );
	doubles SF_antitop = exp(-0.0615-(0.00005* TopReweighting_antitopquark_input) );

	doubles weight;

	switch(ProcessInt){
		case 26: weight = sqrt( SF_top * SF_antitop); break;
		case 27: weight = sqrt( SF_top * SF_antitop); break;
		case 28: weight = sqrt( SF_top * SF_antitop); break;
		case 29: weight = sqrt( SF_top * SF_antitop); break;
		case 30: weight = sqrt( SF_top * SF_antitop); break;
		default: doubles weight_vec(SF_top.size(), 1.0); weight = weight_vec; break;
	}
 
	return weight;

  }};

  
  auto TotHTOverTotpT_doubles{[](const doubles& TotHT, const doubles& TotpT){

  	//std::cout << "print 120" << std::endl;

  	floats TotHTOverTotpTOutput = TotHT / TotpT;
  	return TotHTOverTotpTOutput;

  }};

  
  auto CMSBTagSF_Function{[&SystematicInt](const floats& pts, const floats etas, const floats CSVv2Discr, bool BTagOrNot, const ints& Jet_partonFlavour){

  	//std::cout << "print 121" << std::endl;

  	floats ResultVector{};

  	for(long unsigned int j = 0; j < Jet_partonFlavour.size(); j++){

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

		switch(SystematicInt){
			case 3: systematic_type_string = "up"; break;
        		case 4: systematic_type_string = "down"; break;
 			default: systematic_type_string = "central"; break;
		}

        	std::vector<std::string> SysTypeTest(pts.size(), "central");
        	std::vector<std::string> JetFlavourTest(pts.size(), "0"); 
        	std::vector<std::string> EtaTest{};


		for(long unsigned int i = 0; i < etas.size(); i++){
			std::stringstream ss;
			ss << etas.at(i);
			std::string EtaString(ss.str());
			EtaTest.push_back(EtaString);
		}

		std::vector<std::string> PtTest{};

        	for(long unsigned int i = 0; i < pts.size(); i++){
                	std::stringstream ss;
                	ss << pts.at(i);
                	std::string PtString(ss.str());
                	PtTest.push_back(PtString);
        	}

		std::vector<std::string> DiscrTest{};

        	for(long unsigned int i = 0; i < pts.size(); i++){
                	std::stringstream ss;
                	ss << CSVv2Discr.at(i);
                	std::string CSVv2DiscrString(ss.str());
                	DiscrTest.push_back(CSVv2DiscrString);
        	}

		std::vector<std::string> OutVec{};
		std::vector<std::string> FinalOutVec{};

		for(long unsigned int i = 0; i < CSVv2OperatingPointTest.size(); i++){
	
			for(std::vector<std::string> vec : dataList){
				for(std::string data : vec){OutputVec.push_back(data);}
			}

		for(std::vector<std::string> vec : dataList){
                	for(std::string data : vec){
				
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
			for(long unsigned int k = 0; k < OutVec.size(); k++){if(OutVec.at(k) != "0"){NewOutVec.push_back(OutVec.at(k));}}
	
			std::string outputString;

			if(NewOutVec.size() > 11){outputString = NewOutVec.at( ((i+1)*11)-1 );}
			else{outputString = NewOutVec.at(0); }

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


	//Evaluating the mathematical expression in the string
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
	int index, index2, index3, index4, index5, index6, index7, index8, index9, index10, index11, index12, index13, index14;
	float result;

	for(long unsigned int i = 0; i < FinalOutVec.size(); i++){

		std::string FirstElement = FinalOutVec.at(i);

		if(FirstElement.at(0) != '('){

			long unsigned int LastIndex;

			//first
			for(long unsigned int k = 0; k < FirstElement.length(); k++){

				if(FirstElement.at(k) != ')'&& 
	   	   	   	   FirstElement.at(k) != '(' &&
	   	   	   	   FirstElement.at(k) != '*' &&
	   	   	   	   FirstElement.at(k) != '/' &&
	   	   	   	   FirstElement.at(k) != '+' &&
	  	   	   	   FirstElement.at(k) != '-'){VecForConcString.push_back(FirstElement.at(k)); LastIndex = k;}
				else if(k == 0 && FirstElement.at(k) == '('){continue;}
				else{index = k; break;}
			}

			for(long unsigned int k = 0; k < VecForConcString.size(); k++){
				if(k == 0){ConcatenatedString = VecForConcString.at(k);}
				else{ConcatenatedString += VecForConcString.at(k);}
			}	

			float ConcatenatedStringToFloat = stof(ConcatenatedString);

			if(LastIndex == FirstElement.length()-1){ResultVector.push_back(ConcatenatedStringToFloat);}
			else{
	
				int Min1;

				if(FirstElement.at(index) == '+' && FirstElement.at(index+1) == '(' && FirstElement.at(index+2) == '-' && FirstElement.at(index+3) == '('){
		
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
				else if(FirstElement.at(index) == '+' && FirstElement.at(index+1) == '-'){
			
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
				for(long unsigned int k = Min1; k < FirstElement.length(); k++){

        				if(FirstElement.at(k) == 'e' && FirstElement.at(k+1) == '-'){
					
						VecForConcString2.push_back(FirstElement.at(k));
						VecForConcString2.push_back(FirstElement.at(k+1));
						VecForConcString2.push_back(FirstElement.at(k+2));
						VecForConcString2.push_back(FirstElement.at(k+3));
					
						index2 = k+4;
						break;

					}
					else if(FirstElement.at(k) != ')'&&
           	  		   		FirstElement.at(k) != '(' &&
           	  		   		FirstElement.at(k) != '*' &&
           	  		   		FirstElement.at(k) != '/' &&
           	  		   		FirstElement.at(k) != '+' &&
                  		   		FirstElement.at(k) != '-'){VecForConcString2.push_back(FirstElement.at(k));}
               				else{index2 = k; break;}
			
				}

				for(long unsigned int k = 0; k < VecForConcString2.size(); k++){
        				if(k == 0){ConcatenatedString2 = VecForConcString2.at(k);}
        				else{ConcatenatedString2 += VecForConcString2.at(k);}
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
				for(long unsigned int k = Min2; k < FirstElement.length(); k++){

        				if(FirstElement.at(k) != ')'&&
           	   		   	   FirstElement.at(k) != '(' &&
           	   		   	   FirstElement.at(k) != '*' &&
           	   		   	   FirstElement.at(k) != '/' &&
           	   		   	   FirstElement.at(k) != '+' &&
           	   		   	   FirstElement.at(k) != '-'){VecForConcString3.push_back(FirstElement.at(k));}
        				else{index4 = k; break;}
				}

				for(long unsigned int k = 0; k < VecForConcString3.size(); k++){
        				if(k == 0){ConcatenatedString3 = VecForConcString3.at(k);}
        				else{ConcatenatedString3 += VecForConcString3.at(k);}
				}


				float ConcatenatedStringToFloat3 = stof(ConcatenatedString3);
				int Min3;

				if(FirstElement.at(index4) == '+'){Min3 = index4+1;}
				else{Min3 = index4+1;}

				//fourth
				for(long unsigned int k = Min3; k < FirstElement.length(); k++){

        				if(FirstElement.at(k) != ')'&&
           	   		   	   FirstElement.at(k) != '(' &&
           	   		   	   FirstElement.at(k) != '*' &&
           	   		   	   FirstElement.at(k) != '/' &&
           	   	           	   FirstElement.at(k) != '+' &&
           	   		  	   FirstElement.at(k) != '-'){VecForConcString4.push_back(FirstElement.at(k));}
        				else{index5 = k; break;}
				}


				for(long unsigned int k = 0; k < VecForConcString4.size(); k++){
        				if(k == 0){ConcatenatedString4 = VecForConcString4.at(k);}
        				else{ConcatenatedString4 += VecForConcString4.at(k);}
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
				for(long unsigned int k = Min4; k < FirstElement.length(); k++){

        				if(FirstElement.at(k) != ')'&&
           	   		   	   FirstElement.at(k) != '(' &&
           	   		   	   FirstElement.at(k) != '*' &&
           	   		  	   FirstElement.at(k) != '/' &&
           	   		   	   FirstElement.at(k) != '+' &&
           	   		   	   FirstElement.at(k) != '-'){VecForConcString5.push_back(FirstElement.at(k));}
        				else{index6 = k; break;}
				}

				for(long unsigned int k = 0; k < VecForConcString5.size(); k++){
        				if(k == 0){ConcatenatedString5 = VecForConcString5.at(k);}
        				else{ConcatenatedString5 += VecForConcString5.at(k);}
				}

	
				float ConcatenatedStringToFloat5 = stof(ConcatenatedString5);

				int Min5;

				if(FirstElement.at(index6) == '+' && FirstElement.at(index6+1) != '('){Min5 = index6+1;}
				else if(FirstElement.at(index6) == '+' && FirstElement.at(index6+1) == '('){Min5 = index6+2;}
				else{Min5 = index6+1;}


				//sixth
				for(long unsigned int k = Min5; k < FirstElement.length(); k++){

        				if(FirstElement.at(k) != ')'&&
           	   		   	   FirstElement.at(k) != '(' &&
           	   		   	   FirstElement.at(k) != '*' &&
           	   		   	   FirstElement.at(k) != '/' &&
           	   		   	   FirstElement.at(k) != '+' &&
           	   		   	   FirstElement.at(k) != '-'){VecForConcString6.push_back(FirstElement.at(k));}
        				else{index7 = k; break;}
				}

				for(long unsigned int k = 0; k < VecForConcString6.size(); k++){
        				if(k == 0){ConcatenatedString6 = VecForConcString6.at(k);}
        				else{ConcatenatedString6 += VecForConcString6.at(k);}
				}
	
				float ConcatenatedStringToFloat6 = stof(ConcatenatedString6);

				int Min6;

                       	 	if(FirstElement.at(index7) == ')' && FirstElement.at(index7+1) == '*' && FirstElement.at(index7+2) == '('){Min6 = index7+3;}
				else if(FirstElement.at(index7) == '*'){Min6 = index7+1;}
                        	else{Min6 = index7+3;}
			
				//seventh
				for(long unsigned int k = Min6; k < FirstElement.length(); k++){

        				if(FirstElement.at(k) != ')'&&
           	   		   	   FirstElement.at(k) != '(' &&
           	   		   	   FirstElement.at(k) != '*' &&
           	   		   	   FirstElement.at(k) != '/' &&
           	   		   	   FirstElement.at(k) != '+' &&
           	   		   	   FirstElement.at(k) != '-'){VecForConcString7.push_back(FirstElement.at(k));}
        				else{index8 = k; break;}
				}

				for(long unsigned int k = 0; k < VecForConcString7.size(); k++){
        				if(k == 0){ConcatenatedString7 = VecForConcString7.at(k);}
        				else{ConcatenatedString7 += VecForConcString7.at(k);}
				}


				float ConcatenatedStringToFloat7 = stof(ConcatenatedString7);
				int Min7;
                        
				if(FirstElement.at(index8) == '-' && FirstElement.at(index8+1) == '(' && FirstElement.at(index8+2) == '-' && FirstElement.at(index8+3) == '('){
					Min7 = index8+4;
				}
                        	else{Min7 = index8+1;}
	

				//for an output containing 7 floats
				if(FirstElement.at(index6) == '+' && FirstElement.at(index6+1) == '(' && FirstElement.at(index8) == ')' && FirstElement.at(index8+1) == ')' &&
                           	   FirstElement.at(index8+2) == ')'){

					result = ConcatenatedStringToFloat*((ConcatenatedStringToFloat2+(ConcatenatedStringToFloat3*ConcatenatedStringToFloat4))/
						 (ConcatenatedStringToFloat5+(ConcatenatedStringToFloat6*ConcatenatedStringToFloat7)));

					ResultVector.push_back(result);

				}
				else{

					//eighth
					for(long unsigned int k = Min7; k < FirstElement.length(); k++){

        					if(FirstElement.at(k) != ')'&&
           	   		   	   	   FirstElement.at(k) != '(' &&
           	   		   	   	   FirstElement.at(k) != '*' &&
           	   		   	   	   FirstElement.at(k) != '/' &&
           	   		   	   	   FirstElement.at(k) != '+' &&
           	   		   	   	   FirstElement.at(k) != '-'){VecForConcString8.push_back(FirstElement.at(k));}
        					else{index9 = k; break;}
					}     

					for(long unsigned int k = 0; k < VecForConcString8.size(); k++){
        					if(k == 0){ConcatenatedString8 = VecForConcString8.at(k);}
        					else{ConcatenatedString8 += VecForConcString8.at(k);}
   					}   

					float ConcatenatedStringToFloat8 = stof(ConcatenatedString8);
					int Min8;

					if(FirstElement.at(index9) == '*' &&
			   	   	   FirstElement.at(index9+1) == 'l' &&
                           	   	   FirstElement.at(index9+2) == 'o' &&
                           	   	   FirstElement.at(index9+3) == 'g' &&
                           	   	   FirstElement.at(index9+4) == '('){Min8 = index9+5;}
					else if(FirstElement.at(index9) == ')' && FirstElement.at(index9+1) == ')'){break;}
					else{Min8 = index9+2;}

					//ninth
					for(long unsigned int k = Min8; k < FirstElement.length(); k++){

        					if(FirstElement.at(k) != ')'&&
           	   		   	   	   FirstElement.at(k) != '(' &&
           	   		   	  	   FirstElement.at(k) != '*' &&
           	   		   	   	   FirstElement.at(k) != '/' &&
           	   		  	   	   FirstElement.at(k) != '+' &&
           	   		   	   	   FirstElement.at(k) != '-'){VecForConcString9.push_back(FirstElement.at(k));}
        					else{index10 = k; break;}
					}

					for(long unsigned int k = 0; k < VecForConcString9.size(); k++){
        					if(k == 0){ConcatenatedString9 = VecForConcString9.at(k);}
        					else{ConcatenatedString9 += VecForConcString9.at(k);}
					}

					float ConcatenatedStringToFloat9 = stof(ConcatenatedString9);
					int Min9;

					if(FirstElement.at(index10) == '+'){Min9 = index10+1;}
					else{Min9 = index10 + 1;}

					int LastIndex2;

					//tenth
					for(long unsigned int k = Min9; k < FirstElement.length(); k++){

        					if(FirstElement.at(k) != ')'&&
           	   		   	   	   FirstElement.at(k) != '(' &&
           	   		   	   	   FirstElement.at(k) != '*' &&
           	   		   	   	   FirstElement.at(k) != '/' &&
           	   		   	   	   FirstElement.at(k) != '+' &&
           	   		   	   	   FirstElement.at(k) != '-'){VecForConcString10.push_back(FirstElement.at(k)); LastIndex2 = k;}
        					else{index11 = k; break;}
					}

					for(long unsigned int k = 0; k < VecForConcString10.size(); k++){
        					if(k == 0){ConcatenatedString10 = VecForConcString10.at(k);}
        					else{ConcatenatedString10 += VecForConcString10.at(k);}
					}
		

					float ConcatenatedStringToFloat10 = stof(ConcatenatedString10);

					//result for an equation containing 10 floats and ending in "))))))))' 
			 		if(LastIndex2 == FirstElement.length()-9 &&
                            	   	   FirstElement.at(index) == '+' &&
                            	   	   FirstElement.at(index+1) == '(' &&
                            	  	   FirstElement.at(index+2) == '-' &&
                           	   	   FirstElement.at(index+3) == '('){
				

			    			result = ConcatenatedStringToFloat+(-(ConcatenatedStringToFloat2*(log(ConcatenatedStringToFloat3+
							 ConcatenatedStringToFloat4)*(log(ConcatenatedStringToFloat5+ConcatenatedStringToFloat6)*(ConcatenatedStringToFloat7-
							 (-(ConcatenatedStringToFloat8*log(ConcatenatedStringToFloat9+ConcatenatedStringToFloat10))))))));

                                		ResultVector.push_back(result);

					}
                        		else{

						int Min10;

						if(FirstElement.at(index11+8) == '-'){Min10 = index11+9;}
						else{Min10 = index11+3;}

						int LastIndex3;

						//eleventh
						for(long unsigned int k = Min10; k < FirstElement.length(); k++){

        						if(FirstElement.at(k) != ')'&&
           	   		   	   	   	   FirstElement.at(k) != '(' &&
           	   		  	   	   	   FirstElement.at(k) != '*' &&
           	   		   	   	   	   FirstElement.at(k) != '/' &&
           	   		   	   	   	   FirstElement.at(k) != '+' &&
           	   		   	   	   	   FirstElement.at(k) != '-'){VecForConcString11.push_back(FirstElement.at(k)); LastIndex3 = k;}
        						else{index12 = k; break;}
						
						}

						for(long unsigned int k = 0; k < VecForConcString11.size(); k++){
        						if(k == 0){ConcatenatedString11 = VecForConcString11.at(k);}
        						else{ConcatenatedString11 += VecForConcString11.at(k);}
						}

						float ConcatenatedStringToFloat11 = stof(ConcatenatedString11);

						//for an equation with 11 floats
						if(LastIndex == FirstElement.length()-1 && FirstElement.at(index11+8) == '-'){

							result = ConcatenatedStringToFloat+((-(ConcatenatedStringToFloat2*exp(-5)*(log(ConcatenatedStringToFloat3+
								 ConcatenatedStringToFloat4)*(log(ConcatenatedStringToFloat5+
								 ConcatenatedStringToFloat6)*(ConcatenatedStringToFloat7-(-
								 (ConcatenatedStringToFloat8*log(ConcatenatedStringToFloat9+ConcatenatedStringToFloat10))))))))-
								 ConcatenatedStringToFloat11);

							ResultVector.push_back(result);
					
						}
						else if(FirstElement.at(index) == '+' &&
                                			FirstElement.at(index+1) == '(' &&
                                			FirstElement.at(index+2) == '(' &&
                                			FirstElement.at(index+3) == '-' &&
                                			FirstElement.at(index+4) == '(' &&
							FirstElement.at(index11+8) == '-'){

                                			result = ConcatenatedStringToFloat+((-(ConcatenatedStringToFloat2*(log(ConcatenatedStringToFloat3
								 +ConcatenatedStringToFloat4)*(log(ConcatenatedStringToFloat5
								 +ConcatenatedStringToFloat6)*(ConcatenatedStringToFloat7-(-
								 (ConcatenatedStringToFloat8*log(ConcatenatedStringToFloat9+ConcatenatedStringToFloat10))))))))
								 -ConcatenatedStringToFloat11);

							ResultVector.push_back(result);
                        
                        			}
						else{
				
							//twelfth
							for(long unsigned int k = index12+1; k < FirstElement.length(); k++){
	
        							if(FirstElement.at(k) != ')'&&
           	   			   	   	   	   FirstElement.at(k) != '(' &&
           	   			   	   	   	   FirstElement.at(k) != '*' &&
           	   		   	   		   	   FirstElement.at(k) != '/' &&
           	   		   		   	   	   FirstElement.at(k) != '+' &&
           	   		   	   		   	   FirstElement.at(k) != '-'){VecForConcString12.push_back(FirstElement.at(k));}
        							else{index13 = k; break;}
							}

							for(long unsigned  int k = 0; k < VecForConcString12.size(); k++){
        							if(k == 0){ConcatenatedString12 = VecForConcString12.at(k);}
        							else{ConcatenatedString12 += VecForConcString12.at(k);}
							}

							float ConcatenatedStringToFloat12 = stof(ConcatenatedString12);

							//thirteenth
							for(long unsigned int k = index13+1; k < FirstElement.length(); k++){

        							if(FirstElement.at(k) != ')'&&
           	   		   	   	   	   	   FirstElement.at(k) != '(' &&
           	   		   	   	   	   	   FirstElement.at(k) != '*' &&
           	   		   	   	   	   	   FirstElement.at(k) != '/' &&
           	   		   	   	  	  	   FirstElement.at(k) != '+' &&
           	   		   	   	   	   	   FirstElement.at(k) != '-'){VecForConcString13.push_back(FirstElement.at(k));}
        							else{index14 = k; break;}
							}

							for(long unsigned int k = 0; k < VecForConcString13.size(); k++){
        							if(k == 0){ConcatenatedString13 = VecForConcString13.at(k);}
        							else{ConcatenatedString13 += VecForConcString13.at(k);}
							}

							float ConcatenatedStringToFloat13 = stof(ConcatenatedString13);

							//Calculating the result for an equation containing 13 floats

							if(FirstElement.at(index) == '+'){result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2;}
							else if(FirstElement.at(index) == '-'){result = ConcatenatedStringToFloat - ConcatenatedStringToFloat2;}
							else if(FirstElement.at(index) == '*'){result = ConcatenatedStringToFloat * ConcatenatedStringToFloat2;}
							else{std::cout << "FirstElement.at(index) is " << FirstElement.at(index) 
								       << ". This is not a +, - or *. Output has been set to zero" << std::endl; result = 0;}

							if(FirstElement.at(index2) == '*' && FirstElement.at(index2+1) == '(' && FirstElement.at(index2+2) == '-'){
								result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3);
							}
							else{std::cout << "error" << std::endl;}


							if(FirstElement.at(index4) == '+'){
								result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + 
									 ConcatenatedStringToFloat4);
							}
							else{std::cout << "error message 2" << std::endl;}

							if(FirstElement.at(index5) == '*' && FirstElement.at(index5+1) == '('){
								result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 
									 + ConcatenatedStringToFloat4*(ConcatenatedStringToFloat5));
							}
							else{std::cout << "error message 3" << std::endl;}

							if(FirstElement.at(index6) == '+'){
        							result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + 
									 ConcatenatedStringToFloat4*(ConcatenatedStringToFloat5 + ConcatenatedStringToFloat6));
							}
							else{std::cout << "error message 4" << std::endl;}


							if(FirstElement.at(index7) == '*'){
        							result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + 
									 ConcatenatedStringToFloat4*(ConcatenatedStringToFloat5 + 
									 ConcatenatedStringToFloat6*(-ConcatenatedStringToFloat7)));
							}	   
							else{std::cout << "error message 5" << std::endl;}


							if(FirstElement.at(index8) == '+'){
 	       							result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + 
									 ConcatenatedStringToFloat4*(ConcatenatedStringToFloat5 + 
									 ConcatenatedStringToFloat6*(-ConcatenatedStringToFloat7+ConcatenatedStringToFloat8)));
							}
							else{std::cout << "error message 5" << std::endl;}

							if(FirstElement.at(index9) == '*' && FirstElement.at(index9+1) == '('){
								result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + 
									 ConcatenatedStringToFloat4*(ConcatenatedStringToFloat5 + 
									 ConcatenatedStringToFloat6*(-ConcatenatedStringToFloat7+
									 ConcatenatedStringToFloat8*(ConcatenatedStringToFloat9))));
							}
							else{std::cout << "error message 6" << std::endl;}

							if(FirstElement.at(index10) == '+'){
        							result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + 
									 ConcatenatedStringToFloat4*(ConcatenatedStringToFloat5 + 
									 ConcatenatedStringToFloat6*(-ConcatenatedStringToFloat7+
									 ConcatenatedStringToFloat8*(ConcatenatedStringToFloat9+ConcatenatedStringToFloat10))));
							}
							else{std::cout << "error message 7" << std::endl;}

							if(FirstElement.at(index11) == '*' && FirstElement.at(index11+1) == '(' && FirstElement.at(index11+2) == '-'){
        							result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + 
									 ConcatenatedStringToFloat4*(ConcatenatedStringToFloat5 + 
									 ConcatenatedStringToFloat6*(-ConcatenatedStringToFloat7+
									 ConcatenatedStringToFloat8*(ConcatenatedStringToFloat9+
									 ConcatenatedStringToFloat10*(-ConcatenatedStringToFloat11)))));
							}
							else{std::cout << "error message 8" << std::endl;}


							if(FirstElement.at(index12) == '+'){
        							result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + 
									 ConcatenatedStringToFloat4*(ConcatenatedStringToFloat5 + 
									 ConcatenatedStringToFloat6*(-ConcatenatedStringToFloat7 +
									 ConcatenatedStringToFloat8*(ConcatenatedStringToFloat9 + 
									 ConcatenatedStringToFloat10*(-ConcatenatedStringToFloat11+ConcatenatedStringToFloat12)))));
							}	
							else{std::cout << "error message 9" << std::endl;}

							if(FirstElement.at(index13) == '*'){
        							result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + 
									 ConcatenatedStringToFloat4*(ConcatenatedStringToFloat5 + 
									 ConcatenatedStringToFloat6*(-ConcatenatedStringToFloat7+
									 ConcatenatedStringToFloat8*(ConcatenatedStringToFloat9+
									 ConcatenatedStringToFloat10*(-ConcatenatedStringToFloat11+
									 ConcatenatedStringToFloat12*ConcatenatedStringToFloat13)))));
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
        		for(long unsigned int k = 1; k < FirstElement.length(); k++){

                		if(FirstElement.at(k) != ')'&&
                   	   	   FirstElement.at(k) != '(' &&
                   	   	   FirstElement.at(k) != '*' &&
                   	   	   FirstElement.at(k) != '/' &&
                   	   	   FirstElement.at(k) != '+' &&
                   	   	   FirstElement.at(k) != '-'){VecForConcString.push_back(FirstElement.at(k));}
                		else{index = k; break;}
        		}

        		for(long unsigned int k = 0; k < VecForConcString.size(); k++){
                		if(k == 0){ConcatenatedString = VecForConcString.at(k);}
                		else{ConcatenatedString += VecForConcString.at(k);}
        		}

        		float ConcatenatedStringToFloat = stof(ConcatenatedString);

			int Minimum;

			if(FirstElement.at(index) == '+' && FirstElement.at(index+1) == '(' && FirstElement.at(index+2) == '-' && FirstElement.at(index+3) == '('){
				Minimum = index+4;
			}
			else{Minimum = index+3;}


			//second	
			for(long unsigned int k = Minimum; k < FirstElement.length(); k++){

                		if(FirstElement.at(k) != ')'&&
                   	   	   FirstElement.at(k) != '(' &&
                   	   	   FirstElement.at(k) != '*' &&
                   	   	   FirstElement.at(k) != '/' &&
                   	   	   FirstElement.at(k) != '+' &&
                   	   	   FirstElement.at(k) != '-'){VecForConcString2.push_back(FirstElement.at(k));}
                		else{index2 = k; break;}
        		}

        		for(long unsigned int k = 0; k < VecForConcString2.size(); k++){
                		if(k == 0){ConcatenatedString2 = VecForConcString2.at(k);}
                		else{ConcatenatedString2 += VecForConcString2.at(k);}
        		}
        	
			float ConcatenatedStringToFloat2 = stof(ConcatenatedString2);

			int Minimum2;

			if(FirstElement.at(index2) == '*' && FirstElement.at(index2+1) == '('){Minimum2 = index2+ 6;}
			else{Minimum2 = index2+2;}

			//third
			for(long unsigned int k = Minimum2; k < FirstElement.length(); k++){

                		if(FirstElement.at(k) != ')'&&
                   	   	   FirstElement.at(k) != '(' &&
                   	   	   FirstElement.at(k) != '*' &&
                   	  	   FirstElement.at(k) != '/' &&
                   	   	   FirstElement.at(k) != '+' &&
                   	   	   FirstElement.at(k) != '-'){VecForConcString3.push_back(FirstElement.at(k));}
                		else{index3 = k; break;}
        		}

        		for(long unsigned int k = 0; k < VecForConcString3.size(); k++){
                		if(k == 0){ConcatenatedString3 = VecForConcString3.at(k);}
                		else{ConcatenatedString3 += VecForConcString3.at(k);}
        		}

        		float ConcatenatedStringToFloat3 = stof(ConcatenatedString3);

			//fourth
			for(long unsigned int k = index3+1; k < FirstElement.length(); k++){

                		if(FirstElement.at(k) != ')'&&
                   	   	   FirstElement.at(k) != '(' &&
                   	   	   FirstElement.at(k) != '*' &&
                   	   	   FirstElement.at(k) != '/' &&
                   	   	   FirstElement.at(k) != '+' &&
                   	   	   FirstElement.at(k) != '-'){VecForConcString4.push_back(FirstElement.at(k));}
                		else{index4 = k; break;}
        		}

        		for(long unsigned int k = 0; k < VecForConcString4.size(); k++){
                		if(k == 0){ConcatenatedString4 = VecForConcString4.at(k);}
                		else{ConcatenatedString4 += VecForConcString4.at(k);}
        		}

        		float ConcatenatedStringToFloat4 = stof(ConcatenatedString4);

			int Minimum3;

			if(FirstElement.at(index4) == ')' && FirstElement.at(index4+1) == '*' && FirstElement.at(index4+2) == '('){Minimum3 = index4+7;}
			else{Minimum3 = index4+4;}

			//fifth
			for(long unsigned int k = Minimum3; k < FirstElement.length(); k++){

                		if(FirstElement.at(k) != ')'&&
                   	   	   FirstElement.at(k) != '(' &&
                   	   	   FirstElement.at(k) != '*' &&
                   	   	   FirstElement.at(k) != '/' &&
                   	   	   FirstElement.at(k) != '+' &&
                   	   	   FirstElement.at(k) != '-'){VecForConcString5.push_back(FirstElement.at(k));}
               	 		else{index5 = k; break;}
        		}

        		for(long unsigned int k = 0; k < VecForConcString5.size(); k++){
                		if(k == 0){ConcatenatedString5 = VecForConcString5.at(k);}
                		else{ConcatenatedString5 += VecForConcString5.at(k);}
        		}

        		float ConcatenatedStringToFloat5 = stof(ConcatenatedString5);

			int Minimum4;
		
			if(FirstElement.at(index4) == ')' && FirstElement.at(index4+1) == '*' && FirstElement.at(index4+2) == '('){Minimum4 = index5+1;}
			else{Minimum4 = index5+2;}	
		
			//sixth
			for(long unsigned int k = Minimum4; k < FirstElement.length(); k++){

                		if(FirstElement.at(k) != ')'&&
                   	   	   FirstElement.at(k) != '(' &&
                   	   	   FirstElement.at(k) != '*' &&
                   	   	   FirstElement.at(k) != '/' &&
                   	   	   FirstElement.at(k) != '+' &&
                   	   	   FirstElement.at(k) != '-'){VecForConcString6.push_back(FirstElement.at(k));}
                		else{index6 = k; break;}
        		}

        		for(long unsigned int k = 0; k < VecForConcString6.size(); k++){
                		if(k == 0){ConcatenatedString6 = VecForConcString6.at(k);}
                		else{ConcatenatedString6 += VecForConcString6.at(k);}
        		}

        		float ConcatenatedStringToFloat6 = stof(ConcatenatedString6);
	
			//seventh
			for(long unsigned int k = index6+1; k < FirstElement.length(); k++){

                		if(FirstElement.at(k) != ')'&&
                   	   	   FirstElement.at(k) != '(' &&
                   	   	   FirstElement.at(k) != '*' &&
                   	   	   FirstElement.at(k) != '/' &&
                   	   	   FirstElement.at(k) != '+' &&
                   	  	   FirstElement.at(k) != '-'){VecForConcString7.push_back(FirstElement.at(k));}
                		else{index7 = k; break;}
        		}

        		for(long unsigned int k = 0; k < VecForConcString7.size(); k++){
                		if(k == 0){ConcatenatedString7 = VecForConcString7.at(k);}
                		else{ConcatenatedString7 += VecForConcString7.at(k);}
        		}

        		float ConcatenatedStringToFloat7 = stof(ConcatenatedString7);

			int MinValue;

			if(FirstElement.at(index) == '+' && FirstElement.at(index+1) == '(' && FirstElement.at(index+2) == '-' && FirstElement.at(index+3) == '('){
				MinValue = index7+7;
			}
			else{MinValue = index7+5;}

			//seventh
        		for(long unsigned int k = MinValue; k < FirstElement.length(); k++){

                		if(FirstElement.at(k) != ')'&&
                   	   	   FirstElement.at(k) != '(' &&
                  	   	   FirstElement.at(k) != '*' &&
                   	   	   FirstElement.at(k) != '/' &&
                   	  	   FirstElement.at(k) != '+' &&
                   	   	   FirstElement.at(k) != '-'){VecForConcString8.push_back(FirstElement.at(k));}
                		else{index8 = k; break;}
        		}

        		for(long unsigned int k = 0; k < VecForConcString8.size(); k++){
                		if(k == 0){ConcatenatedString8 = VecForConcString8.at(k);}
                		else{ConcatenatedString8 += VecForConcString8.at(k);}
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

				result = (ConcatenatedStringToFloat*((ConcatenatedStringToFloat2+(ConcatenatedStringToFloat3*ConcatenatedStringToFloat4))/
					 (ConcatenatedStringToFloat5+(ConcatenatedStringToFloat6*ConcatenatedStringToFloat7))))-ConcatenatedStringToFloat8;

        		}
			else if(FirstElement.at(index7+4) == '+'){

                        	result = (ConcatenatedStringToFloat*((ConcatenatedStringToFloat2+(ConcatenatedStringToFloat3*ConcatenatedStringToFloat4))/
					 (ConcatenatedStringToFloat5+(ConcatenatedStringToFloat6*ConcatenatedStringToFloat7))))+ConcatenatedStringToFloat8;

                	}
        		else if(FirstElement.at(index) == '+' &&
				FirstElement.at(index+1) == '(' &&
				FirstElement.at(index+2) == '-' &&
				FirstElement.at(index+3) == '('){
	
					//eighth
                			for(long unsigned int p = index8+5; p < FirstElement.length(); p++){

                        			if(FirstElement.at(p) != ')'&&
                           	   		   FirstElement.at(p) != '(' &&
                           	   	           FirstElement.at(p) != '*' &&
                           	   	           FirstElement.at(p) != '/' &&
                           	   	           FirstElement.at(p) != '+' &&
                           	   	           FirstElement.at(p) != '-'){VecForConcString9.push_back(FirstElement.at(p));}
                        			else{index9 = p; break;}
                			}

                			for(long unsigned int p = 0; p < VecForConcString9.size(); p++){
                        			if(p == 0){ConcatenatedString9 = VecForConcString9.at(p);}
                        			else{ConcatenatedString9 += VecForConcString9.at(p);}
                			}
                	
					float ConcatenatedStringToFloat9 = stof(ConcatenatedString9);

					//ninth
                        		for(long unsigned int p = index9+4; p < FirstElement.length(); p++){

                                		if(FirstElement.at(p) != ')'&&
                                   		   FirstElement.at(p) != '(' &&
                                   		   FirstElement.at(p) != '*' &&
                                   		   FirstElement.at(p) != '/' &&
                                   		   FirstElement.at(p) != '+' &&
                                   		   FirstElement.at(p) != '-'){VecForConcString10.push_back(FirstElement.at(p));}
                                		else{index10 = p; break;}
                        		}

                        		for(long unsigned int p = 0; p < VecForConcString10.size(); p++){
                                		if(p == 0){ConcatenatedString10 = VecForConcString10.at(p);}
                                		else{ConcatenatedString10 += VecForConcString10.at(p);}
                        		}

                        		float ConcatenatedStringToFloat10 = stof(ConcatenatedString10);

					//tenth
                        		for(long unsigned int p = index10+9; p < FirstElement.length(); p++){
          
                                		if(FirstElement.at(p) != ')'&&
                                   		   FirstElement.at(p) != '(' &&
                                   		   FirstElement.at(p) != '*' &&
                                   		   FirstElement.at(p) != '/' &&
                                   		   FirstElement.at(p) != '+' &&
                                   		   FirstElement.at(p) != '-'){VecForConcString11.push_back(FirstElement.at(p));}
                                		else{index11 = p; break;}
                        		}	 
    
                        		for(long unsigned int p = 0; p < VecForConcString11.size(); p++){
                                		if(p == 0){ConcatenatedString11 = VecForConcString11.at(p);}
                                		else{ConcatenatedString11 += VecForConcString11.at(p);}
                      			}
			
                        		float ConcatenatedStringToFloat11 = stof(ConcatenatedString11);

					result = (ConcatenatedStringToFloat+(-(ConcatenatedStringToFloat2*(log(ConcatenatedStringToFloat3+
						  ConcatenatedStringToFloat4)*(log(ConcatenatedStringToFloat5+ConcatenatedStringToFloat6)*(ConcatenatedStringToFloat7-
						 (-(ConcatenatedStringToFloat8*log(ConcatenatedStringToFloat9+ConcatenatedStringToFloat10)))))))))+ConcatenatedStringToFloat11;
			
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

  auto CMSBTagSF{[&CMSBTagSF_Function](const floats& pts, const floats etas, const floats CSVv2Discr, const ints& Jet_partonFlavour){

 	//std::cout << "print 122" << std::endl;
 	return CMSBTagSF_Function(pts, etas, CSVv2Discr, true, Jet_partonFlavour);

  }};

  auto nonbjet_id{[](const ints& tight_jets, const floats& btags, const floats& etas) {

  	//std::cout << "print 123" << std::endl;
  	return tight_jets && (btags == 0) && (etas < MaxTrackerEta);

  }};

  auto CMSNonBTagSF{[&CMSBTagSF_Function](const floats& pts, const floats etas, const floats CSVv2Discr, const ints& Jet_partonFlavour){

 	//std::cout << "print 124" << std::endl;
 	return CMSBTagSF_Function(pts, etas, CSVv2Discr, false, Jet_partonFlavour);

  }};
  
  auto EffBTaggedFunction{[&h_bjet_num, &h_nonbjet_num, &h_charm_num, &h_lightjets_num, &h_gluon_num, 
			   &h_bjet_denom, &h_nonbjet_denom, &h_charm_denom, &h_lightjets_denom, &h_gluon_denom](
			   const int& HistOption, const floats& pts, const floats& etas){

  	//std::cout << "print 125" << std::endl;

  	floats BTaggedEff{};

	ROOT::RDF::RResultPtr<TH2D> EffHistogramNumerator;
	ROOT::RDF::RResultPtr<TH2D> EffHistogramDenominator;

	switch(HistOption){
		case 0: EffHistogramNumerator = h_bjet_num; EffHistogramDenominator = h_bjet_denom; break;
		case 1: EffHistogramNumerator = h_nonbjet_num; EffHistogramDenominator = h_nonbjet_denom; break;
		case 2: EffHistogramNumerator = h_charm_num; EffHistogramDenominator = h_charm_denom; break;
		case 3: EffHistogramNumerator = h_lightjets_num; EffHistogramDenominator = h_lightjets_denom; break;
		case 4: EffHistogramNumerator = h_gluon_num; EffHistogramDenominator = h_gluon_denom; break;
        }

  	for(long unsigned int i = 0; i < pts.size(); i++){

		int PtNum = EffHistogramNumerator->GetXaxis()->FindBin(pts.at(i));
		int EtaNum = EffHistogramNumerator->GetYaxis()->FindBin(etas.at(i));

		int PtDenom = EffHistogramDenominator->GetXaxis()->FindBin(pts.at(i));
		int EtaDenom = EffHistogramDenominator->GetYaxis()->FindBin(etas.at(i));

		float Numerator = EffHistogramNumerator->GetBinContent(PtNum, EtaNum);
		float Denominator = EffHistogramDenominator->GetBinContent(PtDenom, EtaDenom);

		float eff = Numerator / Denominator;
	
		if(!isnan(eff) && !isinf(eff) && eff > 0){BTaggedEff.push_back(eff);}
        	else{BTaggedEff.push_back(1.);}

  	}

 	return BTaggedEff;

  }};

  auto EffBTagged_Function{[&EffBTaggedFunction](const floats& pts, const floats& etas){

  	return EffBTaggedFunction(0, pts, etas);

  }};

  auto EffNonBTagged_Function{[&EffBTaggedFunction](const floats& pts, const floats& etas){

        return EffBTaggedFunction(1, pts, etas);

  }};

  auto EffBTaggedProduct{[](const floats& EffBTagged){
  
  	//std::cout << "print 126" << std::endl;

  	float initial = 1;

  	for(long unsigned int i = 0; i < EffBTagged.size(); i++ ){initial = EffBTagged.at(i) * initial;}

  	return initial;

  }};

  auto EffNonBTaggedProduct{[](const floats& EffNonBTagged){

  	//std::cout << "print 127" << std::endl;
 
  	float initial = 1;

  	for(long unsigned int i = 0; i < EffNonBTagged.size(); i++ ){initial = (1 - EffNonBTagged.at(i)) * initial;}

  	return initial;

  }};

  auto EffBTaggedProductData{[](const floats& EffBTagged, const floats& CMSBTagSFInput){

  	//std::cout << "print 128" << std::endl;

  	float initial = 1;
  	float output;

  	for(long unsigned int i = 0; i < EffBTagged.size(); i++){output = (CMSBTagSFInput.at(0)*EffBTagged.at(i)) * initial;}

  	return output;

  }};



  auto EffNonBTaggedProductData{[](const floats& EffNonBTagged, const floats& CMSNonBTagSFInput){

  	//std::cout << "print 129" << std::endl;

  	float initial = 1;

  	int size = (CMSNonBTagSFInput.size() < EffNonBTagged.size()) ? CMSNonBTagSFInput.size() : EffNonBTagged.size();

  	for(int i = 0; i < size; i++){initial = (1 - (CMSNonBTagSFInput.at(i)*EffNonBTagged.at(i)) ) * initial;}

  	return initial;

  }};


  auto ProbBTagMCFunction{[](const float& EffBTaggedProductInput, const float& EffNonBTaggedProductInput){

  	//std::cout << "print 130" << std::endl;

  	float MCProb = EffBTaggedProductInput * EffNonBTaggedProductInput; 
  	return MCProb;

  }};


  auto ProbBTagDataFunction{[](const float& EffBTaggedProductDataInput, const float& EffNonBTaggedProductDataInput){

  	//std::cout << "print 131" << std::endl;
  
  	float DataProb = EffBTaggedProductDataInput * EffNonBTaggedProductDataInput;
  	return DataProb;
  
  }};


  auto BTagWeightFunction{[](const float& ProbBTagMC, const float& ProbBTagData){

  	//std::cout << "print 131" << std::endl;

	float BTagWeight = (ProbBTagData) / (ProbBTagMC);
	
        if( !isnan(BTagWeight) && !isinf(BTagWeight) ){return BTagWeight;}
	else{float One = 1.0; return One;}

  }};

  auto EGammaFunction{[&EGammaEff2016_histo,     	     &EGammaEffSys2016_histo,
		       &EGammaEffReco2016_histo, 	     &EGammaEffRecoSys2016_histo,
		       &EGammaEff2017_histo,                 &EGammaEffSys2017_histo, 
		       &EGammaEffReco_LowPt_2017_histo,      &EGammaEffRecoSys_LowPt_2017_histo,
		       &EGammaEffReco_HigherPt_2017_histo,   &EGammaEffRecoSys_HigherPt_2017_histo,
		       &EGammaEff2018_histo,	             &EGammaEffSys2018_histo,
		       &EGammaEffReco2018_histo,	     &EGammaEffRecoSys2018_histo
		     ](const int& YearInput, const std::string& type, const floats& pt, const floats& SuperClusterEta){


  	//std::cout << "print 132" << std::endl;

   	floats OutputVector{};
   	floats OutputVectorFinal{};

   	for(long unsigned int i = 0; i < pt.size(); i++){

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

			if(YearInput == 2016){
				if(type == "EGammaEffSys"){EGammaSF = EGammaEffSys2016_histo->GetBinError(Bin_EGammaEffSys2016);}
				else if(type == "EGammaEffRecoSys"){EGammaSF = EGammaEffRecoSys2016_histo->GetBinError(Bin_EGammaEffRecoSys2016);}
				else if(type == "EGammaEff"){EGammaSF = EGammaEff2016_histo->GetBinContent(Bin_EGammaEff2016);}
				else if(type == "EGammaEffReco"){EGammaSF = EGammaEffReco2016_histo->GetBinContent(Bin_EGammaEffReco2016);}
				else{std::cout << "Choose a type out of EGammaEffSys, EGammaEffRecoSys, EGammaEff or EGammaEffReco for 2016" << std::endl;}
			}
			else if(YearInput == 2017){
				if(type == "EGammaEffSys"){EGammaSF = EGammaEffSys2017_histo->GetBinError(Bin_EGammaEffSys2017);}
                        	else if(type == "EGammaEffRecoSys" && pt.at(i) <= 20){EGammaSF = EGammaEffRecoSys_LowPt_2017_histo->GetBinError(Bin_EGammaEffRecoSys_LowPt_2017);}
				else if(type == "EGammaEffRecoSys" && pt.at(i) > 20){EGammaSF = EGammaEffRecoSys_HigherPt_2017_histo->GetBinError(Bin_EGammaEffRecoSys_HigherPt_2017);}
                        	else if(type == "EGammaEff"){EGammaSF = EGammaEff2017_histo->GetBinContent(Bin_EGammaEff2017);}
				else if(type == "EGammaEffReco" && pt.at(i) <= 20){EGammaSF = EGammaEffReco_LowPt_2017_histo->GetBinContent(Bin_EGammaEffReco_LowPt_2017);}
                        	else if(type == "EGammaEffReco" && pt.at(i) > 20){EGammaSF = EGammaEffReco_HigherPt_2017_histo->GetBinContent(Bin_EGammaEffReco_HigherPt_2017);}
                        	else{std::cout << "Choose a type out of EGammaEffSys, EGammaEffRecoSys, EGammaEff or EGammaEffReco for 2017" << std::endl;}
			}
			else if(YearInput == 2018){
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


  	for(long unsigned int i = 0; i < OutputVector.size(); i++){
		if(OutputVector.at(i) == 0){OutputVectorFinal.push_back(1.0);}
		else{OutputVectorFinal.push_back( OutputVector.at(i) );}
  	}

  	return OutputVectorFinal.at(0);

  }};


  auto EGammaSF_egammaEff{[&YearInt, &EGammaFunction](const floats& Electron_pt_Selection, const floats& SuperClusterEta){

  	//std::cout << "print 133" << std::endl;
  	return EGammaFunction(YearInt, "EGammaEff", Electron_pt_Selection, SuperClusterEta);

  }};

  auto EGammaSF_egammaEff_Sys{[&YearInt, &EGammaFunction](const floats& Electron_pt_Selection, const floats& SuperClusterEta){

  	//std::cout << "print 134" << std::endl;
  	return EGammaFunction(YearInt, "EGammaEffSys", Electron_pt_Selection, SuperClusterEta);

  }};

  auto EGammaSF_egammaEffReco{[&YearInt, &EGammaFunction](const floats& Electron_pt_Selection, const floats& SuperClusterEta){

  	//std::cout << "print 135" << std::endl;
  	return EGammaFunction(YearInt, "EGammaEffReco", Electron_pt_Selection, SuperClusterEta);

  }};


  auto EGammaSF_egammaEffReco_Sys{[&YearInt, &EGammaFunction](const floats& Electron_pt_Selection, const floats& SuperClusterEta){

  	//std::cout << "print 136" << std::endl;
  	return EGammaFunction(YearInt, "EGammaEffRecoSys", Electron_pt_Selection, SuperClusterEta);

  }};


  auto MuonSF{[&Year, 				   &histo_RunsBCDEF_ID_2016,   &histo_RunsGH_ID_2016, 	      &histo_RunsBCDEF_ISO_2016, 	     
	       &histo_RunsGH_ISO_2016,  	   &histo_RunsBCDEF_ID_2017,   &histo_RunsBCDEF_ID_Sys_2017,  &histo_RunsBCDEF_ID_Sys_Stat_2017, 
	       &histo_RunsBCDEF_ID_Sys_Syst_2017,  &histo_RunsBCDEF_ISO_2017,  &histo_RunsBCDEF_ISO_Sys_2017, &histo_RunsBCDEF_ISO_Sys_Stat_2017,
	       &histo_RunsBCDEF_ISO_Sys_Syst_2017, &histo_RunsABCD_ID_2018,    &histo_RunsABCD_ISO_2018
              ](const std::string& type, const int& YearInt, const std::string& UpOrDown, const floats& pt, const floats& eta){

  	//std::cout << "print 137" << std::endl;

  	floats AbsEta = abs(eta);
  	float lumiRunBCDEF = 19713.888;
  	float lumiRunGH = 16146.178;

  	floats MuonSFOutput{};

  	for(long unsigned int i = 0; i < pt.size(); i++){

  		if(pt.at(i) >= 20 && pt.at(i) <= 120 && AbsEta.at(i) <= MaxTrackerEta){ 

			if(YearInt == 2016){

				float MuonSF_RunsBCDEF_ID_2016 = histo_RunsBCDEF_ID_2016->GetBinContent( histo_RunsBCDEF_ID_2016->FindBin(pt.at(i), AbsEta.at(i)) );
				float MuonSF_RunsGH_ID_2016 = histo_RunsGH_ID_2016->GetBinContent( histo_RunsGH_ID_2016->FindBin(pt.at(i), AbsEta.at(i)) );
				float MuonSF_RunsBCDEF_ISO_2016 = histo_RunsBCDEF_ISO_2016->GetBinContent( histo_RunsBCDEF_ISO_2016->FindBin(pt.at(i), AbsEta.at(i)) );
				float MuonSF_RunsGH_ISO_2016 = histo_RunsGH_ISO_2016->GetBinContent( histo_RunsGH_ISO_2016->FindBin(pt.at(i), AbsEta.at(i)) );
				float Error_RunsBCDEF_ID_2016 = histo_RunsBCDEF_ID_2016->GetBinError( histo_RunsBCDEF_ID_2016->FindBin(pt.at(i), AbsEta.at(i)) );
				float Error_RunsGH_ID_2016 = histo_RunsGH_ID_2016->GetBinError( histo_RunsGH_ID_2016->FindBin(pt.at(i), AbsEta.at(i)) );
				float Error_RunsBCDEF_ISO_2016 = histo_RunsBCDEF_ISO_2016->GetBinError( histo_RunsBCDEF_ISO_2016->FindBin(pt.at(i), AbsEta.at(i)) );
				float Error_RunsGH_ISO_2016 = histo_RunsGH_ISO_2016->GetBinError( histo_RunsGH_ISO_2016->FindBin(pt.at(i), AbsEta.at(i)) );

				float Error_RunsBCDEFGH, MuonSF_RunsBCDEFGH, Error_RunsBCDEF, MuonSF_RunsBCDEF, Error_RunsGH, MuonSF_RunsGH;

				if(type == "ID" || type == "ID sys"){

					MuonSF_RunsBCDEF = MuonSF_RunsBCDEF_ID_2016; MuonSF_RunsGH = MuonSF_RunsGH_ID_2016;
					Error_RunsBCDEF = Error_RunsBCDEF_ID_2016; Error_RunsGH = Error_RunsGH_ID_2016;
			
				}		
				else if(type == "Iso" || type == "Iso sys"){

					MuonSF_RunsBCDEF = MuonSF_RunsBCDEF_ISO_2016; MuonSF_RunsGH = MuonSF_RunsGH_ISO_2016;
					Error_RunsBCDEF = Error_RunsBCDEF_ISO_2016; Error_RunsGH = Error_RunsGH_ISO_2016;

				}
				else{std::cout << "For 2016, type must be ID, ID sys, Iso or Iso sys" << std::endl;}


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
			else if(YearInt == 2017){

				float MuonSF_RunsBCDEF_ID_2017 = histo_RunsBCDEF_ID_2017->GetBinContent( histo_RunsBCDEF_ID_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				float MuonSF_RunsBCDEF_ID_Sys_2017 = histo_RunsBCDEF_ID_Sys_2017->GetBinContent( histo_RunsBCDEF_ID_Sys_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				float MuonSF_RunsBCDEF_ID_Sys_Stat_2017 = histo_RunsBCDEF_ID_Sys_Stat_2017->GetBinContent( histo_RunsBCDEF_ID_Sys_Stat_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				float MuonSF_RunsBCDEF_ID_Sys_Syst_2017 = histo_RunsBCDEF_ID_Sys_Syst_2017->GetBinContent( histo_RunsBCDEF_ID_Sys_Syst_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				float MuonSF_RunsBCDEF_ISO_2017 = histo_RunsBCDEF_ISO_2017->GetBinContent( histo_RunsBCDEF_ISO_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				float MuonSF_RunsBCDEF_ISO_Sys_2017 = histo_RunsBCDEF_ISO_Sys_2017->GetBinContent( histo_RunsBCDEF_ISO_Sys_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				float MuonSF_RunsBCDEF_ISO_Sys_Stat_2017 = histo_RunsBCDEF_ISO_Sys_Stat_2017->GetBinContent( histo_RunsBCDEF_ISO_Sys_Stat_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				float MuonSF_RunsBCDEF_ISO_Sys_Syst_2017 = histo_RunsBCDEF_ISO_Sys_Syst_2017->GetBinContent( histo_RunsBCDEF_ISO_Sys_Syst_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				float Error_RunsBCDEF_ID_2017 = histo_RunsBCDEF_ID_2017->GetBinError( histo_RunsBCDEF_ID_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				float Error_RunsBCDEF_ID_Sys_2017 = histo_RunsBCDEF_ID_Sys_2017->GetBinError( histo_RunsBCDEF_ID_Sys_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				float Error_RunsBCDEF_ID_Sys_Stat_2017 = histo_RunsBCDEF_ID_Sys_Stat_2017->GetBinError( histo_RunsBCDEF_ID_Sys_Stat_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				float Error_RunsBCDEF_ID_Sys_Syst_2017 = histo_RunsBCDEF_ID_Sys_Syst_2017->GetBinError( histo_RunsBCDEF_ID_Sys_Syst_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				float Error_RunsBCDEF_ISO_2017 = histo_RunsBCDEF_ISO_2017->GetBinError( histo_RunsBCDEF_ISO_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				float Error_RunsBCDEF_ISO_Sys_2017 = histo_RunsBCDEF_ISO_Sys_2017->GetBinError( histo_RunsBCDEF_ISO_Sys_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				float Error_RunsBCDEF_ISO_Sys_Stat_2017 = histo_RunsBCDEF_ISO_Sys_Stat_2017->GetBinError( histo_RunsBCDEF_ISO_Sys_Stat_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				float Error_RunsBCDEF_ISO_Sys_Syst_2017 = histo_RunsBCDEF_ISO_Sys_Syst_2017->GetBinError( histo_RunsBCDEF_ISO_Sys_Syst_2017->FindBin(pt.at(i), AbsEta.at(i)) );

				if(type == "ID"){MuonSFOutput.push_back(MuonSF_RunsBCDEF_ID_2017);}
				else if(type == "Iso"){MuonSFOutput.push_back(MuonSF_RunsBCDEF_ISO_2017);}
				else if(type == "ID sys"){MuonSFOutput.push_back(Error_RunsBCDEF_ID_Sys_2017);}
				else if(type == "ID sys (stat)"){MuonSFOutput.push_back(Error_RunsBCDEF_ID_Sys_Stat_2017);} 
				else if(type == "ID sys (syst)"){MuonSFOutput.push_back(Error_RunsBCDEF_ID_Sys_Syst_2017);}
				else if(type == "Iso sys"){MuonSFOutput.push_back(Error_RunsBCDEF_ISO_Sys_2017);} 
				else if(type == "Iso sys (stat)"){MuonSFOutput.push_back(Error_RunsBCDEF_ISO_Sys_Stat_2017);} 
				else if(type == "Iso sys (syst)"){MuonSFOutput.push_back(Error_RunsBCDEF_ISO_Sys_Syst_2017);}
				else{std::cout << "Incorrect type for 2017 muon SF function" << std::endl;}

			}
			else if(YearInt == 2018){

				float MuonSF_RunsABCD_ID_2018 = histo_RunsABCD_ID_2018->GetBinContent( histo_RunsABCD_ID_2018->FindBin(pt.at(i), AbsEta.at(i)) );
				float MuonSF_RunsABCD_ISO_2018 = histo_RunsABCD_ISO_2018->GetBinContent( histo_RunsABCD_ISO_2018->FindBin(pt.at(i), AbsEta.at(i)) );
				float Error_RunsABCD_ID_2018 = histo_RunsABCD_ID_2018->GetBinError( histo_RunsABCD_ID_2018->FindBin(pt.at(i), AbsEta.at(i)) );
				float Error_RunsABCD_ISO_2018 = histo_RunsABCD_ISO_2018->GetBinError( histo_RunsABCD_ISO_2018->FindBin(pt.at(i), AbsEta.at(i)) );

				if(type == "ID"){MuonSFOutput.push_back(MuonSF_RunsABCD_ID_2018);}
				else if(type == "ID sys"){MuonSFOutput.push_back(Error_RunsABCD_ID_2018);}
				else if(type == "Iso"){MuonSFOutput.push_back(MuonSF_RunsABCD_ISO_2018);}
				else if(type == "Iso sys"){MuonSFOutput.push_back(Error_RunsABCD_ISO_2018);}
				else{std::cout << "Error with Muon SF type (2018)" << std::endl;}

			}
			else{std::cout << "Code only for 2016, 2017 or 2018." << std::endl;}

		}
		else{std::cout << "inside else statement for pushing back muon SF output" << std::endl; float One = 1.0; MuonSFOutput.push_back(One);}


  	}

  	return MuonSFOutput.at(0); 

  }};

  auto MuonSFTest_ID{[&MuonSF, &YearInt](const floats& pt, const floats& eta){

  	//std::cout << "print 138" << std::endl;
  	return MuonSF("ID", YearInt, " ", pt, eta);
  
  }};

  auto MuonSFTest_Iso{[&MuonSF, &YearInt](const floats& pt, const floats& eta){

  	//std::cout << "print 139" << std::endl;
  	return MuonSF("Iso", YearInt, " ", pt, eta);
    
  }};

  auto MuonSFTest_ID_sys_syst{[&MuonSF, &YearInt](const floats& pt, const floats& eta){
  
  	//std::cout << "print 140" << std::endl;

	switch(YearInt){
		case 2016: return MuonSF("ID sys", YearInt, "Up", pt, eta);
        	case 2017: return MuonSF("ID sys (syst)", YearInt, " ", pt, eta);
  		case 2018: std::cout << "Need to add 2018" << std::endl; 
		default: std::cout << "ERROR: Please choose the year out of 2016, 2017 or 2018" << std::endl; break;
	}
		
  }};


  auto MuonSFTest_ID_sys_stat{[&MuonSF, &YearInt](const floats& pt, const floats& eta){

  	//std::cout << "print 141" << std::endl;

	switch(YearInt){
		case 2016: return MuonSF("ID sys", YearInt, "Down", pt, eta);
  		case 2017: return MuonSF("ID sys (stat)", YearInt, " ", pt, eta);
  		case 2018: std::cout << "Need to add 2018" << std::endl;
		default: std::cout << "ERROR: Please choose the year out of 2016, 2017 or 2018" << std::endl; break; 
	}

  }};


  auto MuonSFTest_Iso_sys_syst{[&MuonSF, &YearInt](const floats& pt, const floats& eta){

  	//std::cout << "print 142" << std::endl;

	switch(YearInt){
		case 2016: return MuonSF("Iso sys", YearInt, "Up", pt, eta);
        	case 2017: return MuonSF("Iso sys (syst)", YearInt, " ", pt, eta);
  		case 2018: std::cout << "Need to add 2018" << std::endl; 
		default: std::cout << "ERROR: Please choose the year out of 2016, 2017 or 2018" << std::endl; break;
	}

  }};


  auto MuonSFTest_Iso_sys_stat{[&MuonSF, &YearInt](const floats& pt, const floats& eta){

  	//std::cout << "print 143" << std::endl;

	switch(YearInt){
		case 2016: return MuonSF("Iso sys", YearInt, "Down", pt, eta);
        	case 2017: return MuonSF("Iso sys (stat)", YearInt, " ", pt, eta);
  		case 2018: std::cout << "Need to add 2018" << std::endl;
		default: std::cout << "ERROR: Please choose the year out of 2016, 2017 or 2018" << std::endl; break;
	}

  }};




  auto PSWeightFunction{[&YearInt, &ProcessInt](const floats& PSWeightInput){

  	//std::cout << "print 144" << std::endl;

  	floats Ones(4, 1.0);

	switch(YearInt){

		case 2016: return Ones;
		case 2017: switch(ProcessInt){
				case 0: return PSWeightInput;
				case 29: return PSWeightInput;
				case 30: return PSWeightInput;
				case 33: return PSWeightInput;
				case 36: return PSWeightInput;
				case 39: return PSWeightInput;
				case 56: return PSWeightInput;
				case 89: return PSWeightInput;
				default: return Ones;
			    }
		case 2018: switch(ProcessInt){
                            	case 0: return PSWeightInput;
                                case 29: return PSWeightInput;
                                case 30: return PSWeightInput;
                                case 33: return PSWeightInput;
                                case 36: return PSWeightInput;
                                case 39: return PSWeightInput;
                                case 56: return PSWeightInput;
                                case 89: return PSWeightInput;
                                default: return Ones;
			    }

	}


  }};

  
  auto NominalWeight{[&SystematicInt](const floats& LHEPdfWeight, const floats& LHEWeight_originalXWGTUP){

  	//std::cout << "print 145" << std::endl;

  	float PdfMin = 1.0;
  	float PdfMax = 1.0;

  	//For the min and max Pdf weights
  	for(unsigned int i = 0; i < LHEPdfWeight.size(); i++){

        	float LHEDivision = LHEPdfWeight.at(i) / LHEWeight_originalXWGTUP.at(0); //the size of LHEWeight_originalXWGTUP is always 1

        	if(LHEDivision > PdfMax){PdfMax = LHEDivision;}
        	else{continue;}

        	if(LHEDivision < PdfMin){PdfMin = LHEDivision;}
        	else{continue;}

  	}

	switch(SystematicInt){
		case 11: return PdfMax;
		case 12: return PdfMin;
		default: float One = 1.0; return One;
	}

  }};

  ints SummedWeights(14, 0);

  auto ME_uncert_function{[&SummedWeights](const floats& LHEPdfWeight, const floats& LHEWeight_originalXWGTUP, const floats& ReturnedPSWeight){

  	//std::cout << "print 146" << std::endl;

  	floats pdf = LHEPdfWeight / LHEWeight_originalXWGTUP.at(0);

  	for(long unsigned int i = 0; i < pdf.size(); i++){pdf.at(i) >= 0.0 ? SummedWeights[0]++ : SummedWeights[1]++;} //pdf weight

  		ReturnedPSWeight.at(1) >= 0.0 ? SummedWeights[2]++ : SummedWeights[3]++; //fsr down
  		ReturnedPSWeight.at(0) >= 0.0 ? SummedWeights[4]++ : SummedWeights[5]++; //isr down
  		(ReturnedPSWeight.at(1) * ReturnedPSWeight.at(0)) >= 0.0 ? SummedWeights[6]++ : SummedWeights[7]++; //both isr and fsr down
  		ReturnedPSWeight.at(3) >= 0.0 ? SummedWeights[8]++ : SummedWeights[9]++; //fsr up
  		ReturnedPSWeight.at(2) >= 0.0 ? SummedWeights[10]++ : SummedWeights[11]++; //isr up
  		(ReturnedPSWeight.at(3) * ReturnedPSWeight.at(2)) >= 0.0 ? SummedWeights[12]++ : SummedWeights[13]++; //both isr and fsr up


  		int TotalNumPositive = SummedWeights[0] + SummedWeights[2] + SummedWeights[4] + SummedWeights[6] + SummedWeights[8] + SummedWeights[10] + SummedWeights[12]; 
  		int TotalNumNegative = SummedWeights[1] + SummedWeights[3] + SummedWeights[5] + SummedWeights[7] + SummedWeights[9] + SummedWeights[11] + SummedWeights[13]; 

  		float ME_SF = (TotalNumPositive + TotalNumNegative) / (TotalNumPositive - TotalNumNegative);
 		return ME_SF;

  }};




  auto GeneratorWeight{[&SummedWeights, &SystematicInt](const float& CalculatedNominalWeight, const floats& ReturnedPSWeight){

	//std::cout << "print 147" << std::endl;


	int TotalNumPositive = SummedWeights[0] + SummedWeights[2] + SummedWeights[4] + SummedWeights[6] + SummedWeights[8] + SummedWeights[10] + SummedWeights[12];
        int TotalNumNegative = SummedWeights[1] + SummedWeights[3] + SummedWeights[5] + SummedWeights[7] + SummedWeights[9] + SummedWeights[11] + SummedWeights[13];

	float genweight;

	switch(SystematicInt){
		case 13: genweight = ((TotalNumPositive + TotalNumNegative)/(SummedWeights[12] - SummedWeights[13])) * ((ReturnedPSWeight.at(3) * ReturnedPSWeight.at(2))/abs(CalculatedNominalWeight)); 
			 break;

		case 14: genweight =  ((TotalNumPositive + TotalNumNegative)/(SummedWeights[6] - SummedWeights[7])) * ((ReturnedPSWeight.at(1) * ReturnedPSWeight.at(0))/abs(CalculatedNominalWeight)); 
			 break;
		default: genweight = ((TotalNumPositive + TotalNumNegative) / (TotalNumPositive - TotalNumNegative)) * ( CalculatedNominalWeight / abs(CalculatedNominalWeight) );
			 break;
	}

	return genweight;

  }};


  auto OriginalMetFunction{[&SystematicInt](const floats& MET_sumEt, const floats& MET_phi){

	//std::cout << "print 163" << std::endl; 

	std::vector<TLorentzVector> OriginalMET{};
	TLorentzVector OriginalMET_Element{};

	for(int i = 0; i < MET_phi.size(); i++){

                OriginalMET_Element.SetPtEtaPhiE(MET_sumEt.at(i), 0, MET_phi.at(i), MET_sumEt.at(i));
                OriginalMET.push_back(OriginalMET_Element);

	}

        return OriginalMET;

  }};


  auto ScaledMetFunction{[&SystematicInt](std::vector<TLorentzVector> OriginalMET, const floats& MET_sumEt, const floats& MET_phi, const floats& MET_MetUnclustEnUpDeltaX,  const floats& MET_MetUnclustEnUpDeltaY){

	//std::cout << "print 164" << std::endl;

	std::vector<TLorentzVector> ScaledMET{};
	TLorentzVector ScaledMET_Element{};
	floats metVecOriginal_px;
	floats metVecOriginal_py;

	for(int i = 0; i < OriginalMET.size(); i++){
	
		metVecOriginal_px.push_back( (OriginalMET.at(i)).Px() );
		metVecOriginal_py.push_back( (OriginalMET.at(i)).Py() );
	
	}


	floats MET_px_up =  metVecOriginal_px + MET_MetUnclustEnUpDeltaX;
        floats MET_py_up =  metVecOriginal_py + MET_MetUnclustEnUpDeltaY;
        floats MET_px_down =  metVecOriginal_px - MET_MetUnclustEnUpDeltaX;
        floats MET_py_down =  metVecOriginal_py - MET_MetUnclustEnUpDeltaY;

        floats UnclusteredEnergyUp = sqrt( pow(MET_px_up, 2) + pow(MET_py_up, 2) );
        floats UnclusteredEnergyDown = sqrt( pow(MET_px_down, 2) + pow(MET_py_down, 2) );

        for(long unsigned int i = 0; i < MET_phi.size(); i++){

		if(SystematicInt == 15){
			ScaledMET_Element.SetPtEtaPhiE(UnclusteredEnergyUp.at(i), 0, MET_phi.at(i), UnclusteredEnergyUp.at(i));
		}
                else if(SystematicInt == 16){
			ScaledMET_Element.SetPtEtaPhiE(UnclusteredEnergyDown.at(i), 0, MET_phi.at(i), UnclusteredEnergyDown.at(i));
		}
                else{ScaledMET_Element.SetPtEtaPhiE(MET_sumEt.at(i), 0, MET_phi.at(i), MET_sumEt.at(i));}

                ScaledMET.push_back(ScaledMET_Element);

        }	

	return ScaledMET;

  }};

  auto UnsmearedJetTLorentzVectorFunction{[](const floats& Jet_pt, const floats& Jet_phi, const floats& Jet_eta, const floats& Jet_mass){


	//std::cout << "print 165" << std::endl;

  	std::vector<TLorentzVector> UnsmearedJetVector{};
	TLorentzVector UnsmearedJetVector_Element{};	

	for(int i = 0; i < Jet_pt.size(); i++){
		UnsmearedJetVector_Element.SetPtEtaPhiM(Jet_pt.at(i), Jet_phi.at(i), Jet_eta.at(i), Jet_mass.at(i));
		UnsmearedJetVector.push_back(UnsmearedJetVector_Element);
	}

	return UnsmearedJetVector;
  }};


  auto METUncertFunction{[&SystematicInt](std::vector<TLorentzVector> OriginalMET,           std::vector<TLorentzVector> SmearedJet4Momentum, 
				          std::vector<TLorentzVector> UnsmearedJet4Momentum){

  	//std::cout << "print 148" << std::endl;

	std::vector<TLorentzVector> NewMetVector{};

	for(int i = 0; i < SmearedJet4Momentum.size(); i++){NewMetVector.push_back(OriginalMET.at(0) + SmearedJet4Momentum.at(i) - UnsmearedJet4Momentum.at(i));}
	
	return NewMetVector;

  }};




  auto linereader{[&Year](const int& LineNumber, const std::string YearChoice){
        
        //std::cout << "print 150" << std::endl;
        using namespace std;
        
        std::string NormFileString = "src/Normalisation/NormalisationFactors_" + YearChoice + ".txt";
  
        std::fstream file(NormFileString.c_str());
        GotoLine(file, LineNumber);
        
        std::string line;
        file >> line;
        
        double Value = atof(line.c_str());
        return Value;
 
 }};

  auto NormalisationFactorFunction{[&Year, &ProcessInt, &MCInt, &linereader](){
        
        //std::cout << "print 151" << std::endl;

        switch(MCInt){
                case 1: return linereader(ProcessInt+1, Year);
                default: double one = 1.0; return one;
        }
  }};


  auto EventWeight{[&NormalisationFactorFunction, &ChannelInt, &SystematicInt, &ttbarCRInt]
		    (const float& PUInput, 		        const float& BTagWeightInput, 		     const floats& ReturnedPSWeightInput, 
		     const float& EGammaSF_egammaEffInput,      const float& EGammaSF_egammaEffRecoInput, 
		     const float& EGammaSF_egammaEffSysInput,   const float& EGammaSF_egammaEffRecoSysInput, const float& CalculatedGeneratorWeightInput, 
		     const float& ME_SFInput, 			const doubles& TopWeightInput, 		     const float& CalculatedNominalWeightInput, 
		     const float& MuonSFTest_IDInput, 		const float& MuonSFTest_IsoInput, 	     const float& MuonSFTest_ID_sys_systInput, 
		     const float& MuonSFTest_ID_sys_statInput,  const float& MuonSFTest_Iso_sys_systInput,   const float& MuonSFTest_Iso_sys_statInput){


			//std::cout << "print 149" << std::endl;

			float TrigSF_ee, TrigSFUncert_ee, TrigSF_mumu, TrigSFUncert_mumu;

			switch(ttbarCRInt){
				case 0: TrigSF_ee = SF_ee; TrigSF_mumu = SF_mumu; TrigSFUncert_ee = SF_Uncert_ee; TrigSFUncert_mumu = SF_Uncert_mumu; break;
				default: TrigSF_ee = SF_emu; TrigSF_mumu = SF_emu; TrigSFUncert_ee = SF_Uncert_emu; TrigSFUncert_mumu = SF_Uncert_emu; break;

			}

			float EventWeightOutput;

			std::cout << "channel int = " << ChannelInt << " SystematicInt = " << SystematicInt << std::endl;

			switch(ChannelInt){
				case 1: 
					switch(SystematicInt){ 
						case 9: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * 
							(TrigSF_ee += TrigSFUncert_ee) * CalculatedNominalWeightInput * EGammaSF_egammaEffSysInput * 
							EGammaSF_egammaEffRecoSysInput * CalculatedGeneratorWeightInput * TopWeightInput.at(0);

							break;

                        			case 10: EventWeightOutput = PUInput * NormalisationFactorFunction() * (TrigSF_ee -= TrigSFUncert_ee) * 
					  		 CalculatedNominalWeightInput * 
						         EGammaSF_egammaEffSysInput * EGammaSF_egammaEffRecoSysInput * CalculatedGeneratorWeightInput * 
							 TopWeightInput.at(0);

							 break;
				
						case 11: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF_ee * 
							 CalculatedNominalWeightInput * 
							 EGammaSF_egammaEffInput * EGammaSF_egammaEffRecoInput * CalculatedGeneratorWeightInput * TopWeightInput.at(0);
                        
							 break;

						case 12: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF_ee * 
						         		     CalculatedNominalWeightInput * EGammaSF_egammaEffInput * EGammaSF_egammaEffRecoInput * 
									     CalculatedGeneratorWeightInput * TopWeightInput.at(0);
				
							 break;

						case 17: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF_ee * ReturnedPSWeightInput.at(2) * 
							 CalculatedNominalWeightInput * EGammaSF_egammaEffInput * EGammaSF_egammaEffRecoInput * CalculatedGeneratorWeightInput * 
							 TopWeightInput.at(0);

							 break;

                        			case 18: EventWeightOutput =  PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF_ee * 
							 ReturnedPSWeightInput.at(0) * 
							 CalculatedNominalWeightInput * EGammaSF_egammaEffInput * EGammaSF_egammaEffRecoInput * CalculatedGeneratorWeightInput * TopWeightInput.at(0);

							 break;

                        			case 19: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF_ee * ReturnedPSWeightInput.at(3) * 
							 CalculatedNominalWeightInput * EGammaSF_egammaEffInput * EGammaSF_egammaEffRecoInput * CalculatedGeneratorWeightInput  
							 * TopWeightInput.at(0);

							 break;

                        			case 20: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF_ee * ReturnedPSWeightInput.at(1) * 
							 CalculatedNominalWeightInput * EGammaSF_egammaEffInput * EGammaSF_egammaEffRecoInput * CalculatedGeneratorWeightInput  
							 * TopWeightInput.at(0);

							 break;

						default: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF_ee * 
							 CalculatedNominalWeightInput * 
							 EGammaSF_egammaEffInput * EGammaSF_egammaEffRecoInput * CalculatedGeneratorWeightInput * TopWeightInput.at(0);
							
							 break;

					}
  
				case 2: switch(SystematicInt){
						case 9: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * (TrigSF_mumu += TrigSFUncert_mumu) * 
							CalculatedNominalWeightInput * MuonSFTest_ID_sys_systInput * MuonSFTest_Iso_sys_systInput * 
							CalculatedGeneratorWeightInput * TopWeightInput.at(0);

							break;

                        			case 10: EventWeightOutput = PUInput * NormalisationFactorFunction() * (TrigSF_mumu -= TrigSFUncert_mumu) * 
							 CalculatedNominalWeightInput * 
							 MuonSFTest_ID_sys_statInput * MuonSFTest_Iso_sys_statInput * CalculatedGeneratorWeightInput * 
							 TopWeightInput.at(0);
				
							 break;

						case 11: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF_mumu * 
							 CalculatedNominalWeightInput * 
							 MuonSFTest_IDInput * MuonSFTest_IsoInput * CalculatedGeneratorWeightInput * TopWeightInput.at(0);

							 break;

                        			case 12: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF_mumu * 
						         CalculatedNominalWeightInput * 
							 MuonSFTest_IDInput * MuonSFTest_IsoInput * CalculatedGeneratorWeightInput * TopWeightInput.at(0);

							 break;
	
                        			case 17: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF_mumu * 
							 ReturnedPSWeightInput.at(2) * 
							 CalculatedNominalWeightInput * MuonSFTest_IDInput * MuonSFTest_IsoInput * CalculatedGeneratorWeightInput * 
							 TopWeightInput.at(0);

							 break;

						case 18: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF_mumu *
							 ReturnedPSWeightInput.at(0) * 
							 CalculatedNominalWeightInput * MuonSFTest_IDInput * MuonSFTest_IsoInput * CalculatedGeneratorWeightInput * 
						         TopWeightInput.at(0);

							 break; 

                        			case 19: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF_mumu * 
							 ReturnedPSWeightInput.at(3) * 
							 CalculatedNominalWeightInput * MuonSFTest_IDInput * MuonSFTest_IsoInput * CalculatedGeneratorWeightInput * 
							 TopWeightInput.at(0);

							 break;

                        			case 20: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF_mumu * 
							 ReturnedPSWeightInput.at(1) * 
							 CalculatedNominalWeightInput * MuonSFTest_IDInput * MuonSFTest_IsoInput * CalculatedGeneratorWeightInput * 
							 TopWeightInput.at(0);

							 break;

                        			default: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF_mumu * 
							 CalculatedNominalWeightInput * 
							 MuonSFTest_IDInput * MuonSFTest_IsoInput * CalculatedGeneratorWeightInput * TopWeightInput.at(0);

							 break;

					}

				default: std::cout << "ERROR: Channel for event weights must be ee or mumu. ChannelInt is: " << ChannelInt << std::endl; break;

		}



  	float FinalEventWeight = (EventWeightOutput/abs(EventWeightOutput)) * ME_SFInput; 

	std::cout << '\n' << std::endl;
	std::cout << '\n' << std::endl;
	std::cout << "PUInput = " << PUInput << std::endl;
        std::cout << "NormalisationFactorFunction() = " << NormalisationFactorFunction() << std::endl;
        std::cout << "BTagWeightInput = " << BTagWeightInput << std::endl;
        std::cout << "TrigSF_ee = " << TrigSF_ee << std::endl;
        std::cout << "CalculatedNominalWeightInput = " <<  CalculatedNominalWeightInput << std::endl;
        std::cout << "EGammaSF_egammaEffInput = " << EGammaSF_egammaEffInput << std::endl;
        std::cout << "EGammaSF_egammaEffRecoInput = " << EGammaSF_egammaEffRecoInput << std::endl;
        std::cout << "CalculatedGeneratorWeightInput = " << CalculatedGeneratorWeightInput << std::endl;
        std::cout << "TopWeightInput.at(0) = " << TopWeightInput.at(0) << std::endl;
	std::cout << "ME_SFInput = " << ME_SFInput << std::endl;
	std::cout << "FinalEventWeight = " << FinalEventWeight << std::endl;
	std::cout << '\n' << std::endl;
        std::cout << '\n' << std::endl;	

	if(!isnan(FinalEventWeight) && !isinf(FinalEventWeight) && (FinalEventWeight > 0)){return FinalEventWeight;}
	else{std::cout << "Final event weight is either a nan, inf or 0." << std::endl; float One = 1.0; return One;}     

  }};



  std::vector<float> CutRanges = {};
  float W_stddev;
  float Top_stddev;

  auto Chi2Function{[&ProcessInt, &CutRanges, &SBRInt, &SystematicInt, &W_stddev, &Top_stddev](const float& w_mass, const float& Top_Mass){

  	//std::cout << "print 152" << std::endl;
	
  	float FiveSigmaW = 5*W_stddev;

  	//calculating chi2
  	float chi2 = pow(( (w_mass - W_MASS) / W_stddev), 2) + pow(( (Top_Mass - TOP_MASS) / Top_stddev), 2);

  	float LowerBound = W_MASS - FiveSigmaW;
  	float UpperBound = W_MASS + FiveSigmaW;

	switch(ProcessInt){
		case 0: switch(SystematicInt){
				case 0: switch(SBRInt){ 
						case 1: //returning chi2 values only for when w_mass is within 5 sigma of the known W mass 
  							if(w_mass > LowerBound && w_mass < UpperBound){ CutRanges.push_back(chi2); return chi2;}	
							else{
								std::cout << "w_mass is not within 5 sigma of the mean W mass value (ee)" << std::endl;
								float out = 999.0;
                						return out;
							}
						default: return chi2;

  					}

				default: return chi2;

			}

		default: return chi2;
	}

  }};

  auto Chi2Cut{[&SBRInt, &SRInt](const float& Chi2){	

  	//std::cout << "print 153" << std::endl;

	  switch(SBRInt){
		case 1: return Chi2_SR < Chi2 && Chi2 < Chi2_SBR; break;
		default: break;
	  }

	  switch(SRInt){
		case 1: return Chi2 < Chi2_SR; break;
		default: break;

	  }

  }};


  //Input file selection
  //EnableImplicitMT(); //to enable multithreading
  RDataFrame d("Events", input_files); //accessing the events TTree of the input file
  
  auto d_Range = d.Range(0, 5000);

  //Event cleaning
  auto d_EventCleaning = d_Range.Filter(filter_function, {"Flag_goodVertices",              "Flag_globalSuperTightHalo2016Filter",     "Flag_HBHENoiseFilter", 
						    "Flag_HBHENoiseIsoFilter",        "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_BadPFMuonFilter",
                                                    "Flag_BadChargedCandidateFilter", "Flag_ecalBadCalibFilter",                 "Flag_eeBadScFilter"}, "Event cleaning filter");

 //Filtering events using the golden json file (for data not MC)
 auto d_GoldenJson = d_EventCleaning.Filter(RunAndLumiFilterFunction, {"run", "luminosityBlock"}, "GoldenJson filter");

 std::vector<std::string> DoubleCountCheckStrings;

 if(Process == "Data_DoubleEGRunB"       || Process == "Data_DoubleEGRunC"       || Process == "Data_DoubleEGRunD"       || 
    Process == "Data_DoubleEGRunE"       || Process == "Data_DoubleEGRunF"       || Process == "Data_DoubleEGRunG"       || 
    Process == "Data_DoubleEGRunH"       || Process == "Data_DoubleMuonRunB"     || Process == "Data_DoubleMuonRunC"     || 
    Process == "Data_DoubleMuonRunD"     || Process == "Data_DoubleMuonRunE"     || Process == "Data_DoubleMuonRunF"     ||
    Process == "Data_DoubleMuonRunG"     || Process == "Data_DoubleMuonRunH"     || Process == "Data_MuonEGRunB"         ||
    Process == "Data_MuonEGRunC"         || Process == "Data_MuonEGRunD"	 || Process == "Data_MuonEGRunE"         ||
    Process == "Data_MuonEGRunF"         || Process == "Data_MuonEGRunG"         || Process == "Data_MuonEGRunH"         ||
    Process == "Data_SingleMuonRunB"     || Process == "Data_SingleMuonRunC"     || Process == "Data_SingleMuonRunD"     ||
    Process == "Data_SingleMuonRunE"     || Process == "Data_SingleMuonRunF"     || Process == "Data_SingleMuonRunG"     ||
    Process == "Data_SingleMuonRunH"     || Process == "Data_SingleElectronRunB" || Process == "Data_SingleElectronRunC" ||
    Process == "Data_SingleElectronRunD" || Process == "Data_SingleElectronRunE" || Process == "Data_SingleElectronRunF" ||
    Process == "Data_SingleElectronRunG" || Process == "Data_SingleElectronRunH"){

	DoubleCountCheckStrings = {"HLT_Ele25_eta2p1_WPTight_Gsf",              	 "HLT_Ele27_WPTight_Gsf",                     
			           "HLT_Ele32_eta2p1_WPTight_Gsf",	 		 "HLT_Ele32_WPTight_Gsf_L1DoubleEG",          
				   "HLT_Ele35_WPTight_Gsf",                     	 "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
                                   "HLT_IsoMu24",		                	 "HLT_IsoMu24_eta2p1",                        
				   "HLT_IsoMu27",			        	 "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",       
				   "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", 	 "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
                                   "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 	 "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
				   "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                                   "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",    "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
                                   "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",     "event"};

  }
 else{

	DoubleCountCheckStrings = {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
				   "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
				   "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
				   "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
				   "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
				   "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
				   "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
				   "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
				   "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
				   "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "event"};


 }

 //Preventing the double-counting of events in the single and double lepton datasets
 auto d_DoubleCountCheck = d_GoldenJson.Filter(DoubleCountCheck_EventFunction, DoubleCountCheckStrings); 

 if(DoubleCountCheckInt == 1){

    std::string DoubleCountCheckFile = "DoubleCountCheck_" + Process + "_" + Systematic + "_" + NonPromptLepton + "_" +
                                       SignalRegion + "_" + SideBandRegion + "_" + ZPlusJetsControlRegion + "_" + ttbarControlRegion + "_" + Year + ".root";

    auto DoubleCountCheckSnapshot = d_DoubleCountCheck.Snapshot("Events", DoubleCountCheckFile.c_str());

    return;

 }



 //Lepton selection
 auto d_LeptonSelection = d_GoldenJson.Define("PU", PU_function, {"PV_npvs"})
                                      .Define("TightLeptons", TightLeptonsFunction, {"Electron_pt", "Electron_eta", "Electron_cutBased", "Electron_isPFcand",
                                                   				     "Muon_isPFcand", "Muon_pt", "Muon_eta", "Muon_tightId", "Muon_pfRelIso04_all"})
				      .Define("LeptonPt", LeptonVariableFunctionFloats, {"Electron_pt", "Muon_pt"})
				      .Define("LeptonPhi", LeptonVariableFunctionFloats, {"Electron_phi", "Muon_phi"})
				      .Define("LeptonEta", LeptonVariableFunctionFloats, {"Electron_eta", "Muon_eta"})
				      .Define("LeptonCharge", LeptonVariableFunctionInts, {"Electron_charge", "Muon_charge"})
				      .Define("LeptonMass", LeptonVariableFunctionFloats, {"Electron_mass", "Muon_mass"})
				      .Define("LeptonJetRelIso", LeptonVariableFunctionFloats, {"Electron_jetRelIso", "Muon_jetRelIso"})
                                      .Define("TightLeptonsPt", select<floats>, {"LeptonPt", "TightLeptons"})
                                      .Define("TightLeptonsPhi", select<floats>, {"LeptonPhi", "TightLeptons"})
                                      .Define("TightLeptonsEta", select<floats>, {"LeptonEta", "TightLeptons"})
				      .Define("TightLeptonsMass", select<floats>, {"LeptonMass", "TightLeptons"})
                                      .Define("TightLeptonsCharge", select<ints>, {"LeptonCharge", "TightLeptons"})
				      .Define("TightLeptonsJetRelIso", select<floats>, {"LeptonJetRelIso", "TightLeptons"})
                                      .Define("LooseLeptons", LooseLeptonsFunction, {"Electron_pt", "Electron_eta", "Electron_cutBased", "Electron_isPFcand",
                                                   				       "Muon_isPFcand", "Muon_pt", "Muon_eta", "Muon_softId", "Muon_pfRelIso04_all"})
                                      .Define("LooseLeptonsPt", select<floats>, {"LeptonPt", "LooseLeptons"})
				      .Define("LooseLeptonsPhi", select<floats>, {"LeptonPhi", "LooseLeptons"})
				      .Define("LooseLeptonsEta", select<floats>, {"LeptonEta", "LooseLeptons"})
                                      .Define("LooseLeptonsCharge", select<ints>, {"LeptonCharge", "LooseLeptons"})
				      .Define("LooseLeptonsMass", select<floats>, {"LeptonMass", "LooseLeptons"})
                                      .Define("OppositeSign", OppositeSign, {"LeptonCharge"})
                                      .Define("SameSign", SameSign, {"LeptonCharge"})
                                      .Define("LeadingLeptonPt", LeadingVariable, {"TightLeptonsPt"})
                                      .Define("SubleadingLeptonPt", SubleadingVariable, {"TightLeptonsPt"})
                                      .Define("LeadingLeptonPhi", LeadingVariable, {"TightLeptonsPhi"})
                                      .Define("SubleadingLeptonPhi", SubleadingVariable, {"TightLeptonsPhi"})
                                      .Define("LeadingLeptonEta", LeadingVariable, {"TightLeptonsEta"})
                                      .Define("SubleadingLeptonEta", SubleadingVariable, {"TightLeptonsEta"})
                                      .Define("LeadingLeptonMass", LeadingVariable, {"TightLeptonsMass"})
                                      .Define("SubleadingLeptonMass", SubleadingVariable, {"TightLeptonsMass"})
                                      .Define("LeadingLepton_RelIso_Selection", LeadingVariable, {"TightLeptonsJetRelIso"})
                                      .Define("SubleadingLepton_RelIso_Selection", SubleadingVariable, {"TightLeptonsJetRelIso"})
				      .Filter(LeptonCut, {"TightLeptonsPt", "LooseLeptonsPt", "OppositeSign", "nElectron", "nMuon", "Electron_dz", "Electron_dxy",
							  "LeadingLeptonPt", "SubleadingLeptonPt", "LeptonEta"}, "lepton cut");

  std::string LeptonSelectionFile = "LeptonSelection_" + Process + "_" + Systematic + "_" + Channel + "_" + NonPromptLepton + "_" +
                                     SignalRegion + "_" + SideBandRegion + "_" + ZPlusJetsControlRegion + "_" + ttbarControlRegion + "_" + Year + ".root";

  auto Snapshot_LeptonSelection = d_LeptonSelection.Snapshot("Events", LeptonSelectionFile.c_str());


  //Calculating the trigger scale factors
  std::vector<std::string> MET_Triggers_Strings;
  std::vector<std::string> Lepton_Triggers_Strings;


   switch(ProcessInt){

	case 95: switch(YearInt){

			case 2016: MET_Triggers_Strings = {"HLT_MET200", 		    		"HLT_MET250",
		        				   "HLT_PFMET120_PFMHT120_IDTight", 		"HLT_PFMET170_HBHECleaned",
							   "HLT_PFHT300_PFMET100",	    		"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
							   "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
							   "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
							   "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
							   "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
							   "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "event"};

				   Lepton_Triggers_Strings = {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 	    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
							      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 	    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
							      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", 		    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
							      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",	    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
							      "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
							      "HLT_Ele25_eta2p1_WPTight_Gsf", 	     		    "HLT_Ele27_WPTight_Gsf",
							      "HLT_Ele32_eta2p1_WPTight_Gsf",	     		    "HLT_IsoMu24",
							      "HLT_IsoMu24_eta2p1",		     		    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
							      "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
							      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 	    "event"}; 


				   break;

			case 2017: MET_Triggers_Strings = {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
							   "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
							   "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_PFMET130_PFMHT130_IDTight", 
							   "HLT_PFMET140_PFMHT140_IDTight",             "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
							   "HLT_PFHT1050",                              "HLT_PFHT180",
							   "HLT_PFHT500_PFMET100_PFMHT100_IDTight",     "HLT_PFHT500_PFMET110_PFMHT110_IDTight",
							   "HLT_PFHT700_PFMET85_PFMHT85_IDTight",       "HLT_PFHT700_PFMET95_PFMHT95_IDTight",
							   "HLT_PFHT800_PFMET75_PFMHT75_IDTight",       "HLT_PFHT800_PFMET85_PFMHT85_IDTight", "event"};


				  Lepton_Triggers_Strings = {"HLT_Ele32_WPTight_Gsf_L1DoubleEG",                   "HLT_Ele35_WPTight_Gsf",
                                                              "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",             "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                                              "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
                                                              "HLT_IsoMu27",                                        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                                              "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                                                              "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",          "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                                              "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",          "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                                              "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",          "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                                              "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",          "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                                              "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",          "event"};


				   break;

			case 2018: MET_Triggers_Strings = {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                                           "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                                           "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_PFMET130_PFMHT130_IDTight",
                                                           "HLT_PFMET140_PFMHT140_IDTight",             "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
                                                           "HLT_PFHT1050",                              "HLT_PFHT180",
                                                           "HLT_PFHT500_PFMET100_PFMHT100_IDTight",     "HLT_PFHT500_PFMET110_PFMHT110_IDTight",
                                                           "HLT_PFHT700_PFMET85_PFMHT85_IDTight",       "HLT_PFHT700_PFMET95_PFMHT95_IDTight",
                                                           "HLT_PFHT800_PFMET75_PFMHT75_IDTight",       "HLT_PFHT800_PFMET85_PFMHT85_IDTight", "event"}; 


				   Lepton_Triggers_Strings = {"HLT_Ele32_WPTight_Gsf_L1DoubleEG", 	   	    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
							      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",     	    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
							      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 	    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", 
							      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 	    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
							      "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
							      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 	    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
							      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 	    "HLT_IsoMu24",
							      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 	    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
							      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 	    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
							      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8", 	    "event"};


				   break;
	
		 };


		 break;

	case 96: switch(YearInt){

                        case 2016: MET_Triggers_Strings = {"HLT_MET200",                    		"HLT_MET250",
                                                           "HLT_PFMET120_PFMHT120_IDTight", 		"HLT_PFMET170_HBHECleaned",
                                                           "HLT_PFHT300_PFMET100",          		"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                                           "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
							   "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
							   "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
							   "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
							   "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "event"};

				   Lepton_Triggers_Strings = {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",          "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                                              "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",          "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                                              "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
                                                              "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",          "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                                              "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                                                              "HLT_Ele25_eta2p1_WPTight_Gsf",                       "HLT_Ele27_WPTight_Gsf",
                                                              "HLT_Ele32_eta2p1_WPTight_Gsf",                       "HLT_IsoMu24",
                                                              "HLT_IsoMu24_eta2p1",                                 "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
                                                              "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
                                                              "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",          "event"};


                                   break;

                        case 2017: MET_Triggers_Strings = {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                                           "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                                           "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_PFMET130_PFMHT130_IDTight",
                                                           "HLT_PFMET140_PFMHT140_IDTight",             "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
                                                           "HLT_PFHT1050",                              "HLT_PFHT180",
                                                           "HLT_PFHT500_PFMET100_PFMHT100_IDTight",     "HLT_PFHT500_PFMET110_PFMHT110_IDTight",
                                                           "HLT_PFHT700_PFMET85_PFMHT85_IDTight",       "HLT_PFHT700_PFMET95_PFMHT95_IDTight",
                                                           "HLT_PFHT800_PFMET75_PFMHT75_IDTight",       "HLT_PFHT800_PFMET85_PFMHT85_IDTight", "event"};


				   Lepton_Triggers_Strings = {"HLT_Ele32_WPTight_Gsf_L1DoubleEG", 		    "HLT_Ele35_WPTight_Gsf",
							      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", 	    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
							      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", 		    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", 
							      "HLT_IsoMu27",			     		    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
							      "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
							      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 	    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
							      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 	    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
							      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 	    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
							      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 	    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
							      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 	    "event"};


                                   break;

                        case 2018: MET_Triggers_Strings = {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                                           "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                                           "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_PFMET130_PFMHT130_IDTight",
                                                           "HLT_PFMET140_PFMHT140_IDTight",             "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
                                                           "HLT_PFHT1050",                              "HLT_PFHT180",
                                                           "HLT_PFHT500_PFMET100_PFMHT100_IDTight",     "HLT_PFHT500_PFMET110_PFMHT110_IDTight",
                                                           "HLT_PFHT700_PFMET85_PFMHT85_IDTight",       "HLT_PFHT700_PFMET95_PFMHT95_IDTight",
                                                           "HLT_PFHT800_PFMET75_PFMHT75_IDTight",       "HLT_PFHT800_PFMET85_PFMHT85_IDTight", "event"};


				   Lepton_Triggers_Strings = {"HLT_Ele32_WPTight_Gsf_L1DoubleEG",                   "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                                              "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",             "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                                              "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",          "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
                                                              "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",          "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                                              "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                                                              "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",          "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                                              "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",          "HLT_IsoMu24",
                                                              "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",          "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                                              "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",          "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                                              "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",        "event"};



                                   break;

                 }; 


		 break;

	default: MET_Triggers_Strings = {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
					 "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
					 "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
					 "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
					 "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
                                         "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
                                         "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
                                         "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "event"};

		Lepton_Triggers_Strings = {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                           "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                           "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                           "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                           "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                           "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                           "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                           "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
					   "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 
					   "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", "event"};
		
		break;

  }


  auto d_MET_And_LeptonSelection = d_LeptonSelection.Filter(MET_Triggers_Function, MET_Triggers_Strings);
  auto d_LeptonTriggers_And_LeptonSelection = d_LeptonSelection.Filter(Lepton_Triggers_Function, Lepton_Triggers_Strings);
  auto d_MET_LeptonTriggers_LeptonSelection = d_LeptonTriggers_And_LeptonSelection.Filter(MET_Triggers_Function, MET_Triggers_Strings);

  float LeptonSelection_EventWeight, MET_And_LeptonSelection_EventWeight, LeptonTriggersAndSelectionCriteria_EventWeight, MET_LeptonTriggers_SelectionCriteria_EventWeight;

  switch(ProcessInt){
	case 95: LeptonSelection_EventWeight = *d_LeptonSelection.Sum("PU");
	         MET_And_LeptonSelection_EventWeight = *d_MET_And_LeptonSelection.Sum("PU");
	         LeptonTriggersAndSelectionCriteria_EventWeight = *d_LeptonTriggers_And_LeptonSelection.Sum("PU");
                 MET_LeptonTriggers_SelectionCriteria_EventWeight = *d_MET_LeptonTriggers_LeptonSelection.Sum("PU");
	         break;

	default: LeptonSelection_EventWeight = 1.0;
		 MET_And_LeptonSelection_EventWeight = 1.0;
		 LeptonTriggersAndSelectionCriteria_EventWeight = 1.0;
		 MET_LeptonTriggers_SelectionCriteria_EventWeight = 1.0;
		 break;
  }

  float N_SelectionCriteria = *(d_LeptonSelection.Count()) * LeptonSelection_EventWeight;
  float N_MET_And_SelectionCriteria = *(d_MET_And_LeptonSelection.Count()) * MET_And_LeptonSelection_EventWeight;
  float N_LeptonTriggersAndSelectionCriteria = *(d_LeptonTriggers_And_LeptonSelection.Count()) * LeptonTriggersAndSelectionCriteria_EventWeight;
  float N_MET_LeptonTriggers_SelectionCriteria = *(d_MET_LeptonTriggers_LeptonSelection.Count()) * MET_LeptonTriggers_SelectionCriteria_EventWeight;

  //Calculating the efficiency
  float Eff = (N_MET_LeptonTriggers_SelectionCriteria) / (N_SelectionCriteria + 1.0e-06);
  
  //Calculating alpha
  float Eff_MET_LeptonTriggers_SelectionCriteria = (N_MET_LeptonTriggers_SelectionCriteria) / (N_SelectionCriteria + 1.0e-06);
  float Eff_LeptonTriggers_SelectionCriteria = (N_LeptonTriggersAndSelectionCriteria) / (N_SelectionCriteria + 1.0e-06);
  float Eff_MET_SelectionCriteria = (N_MET_And_SelectionCriteria) / (N_SelectionCriteria + 1.0e-06);

  float Alpha = (Eff_LeptonTriggers_SelectionCriteria * Eff_MET_SelectionCriteria) / (Eff_MET_LeptonTriggers_SelectionCriteria + 1.0e-06);

  //Calculating the uncertainties for the efficiencies
  double level = 0.60;
  float Eff_UpperUncert = Eff - TEfficiency::ClopperPearson(N_MET_And_SelectionCriteria, N_MET_LeptonTriggers_SelectionCriteria, level, true);
  float Eff_LowerUncert = Eff - TEfficiency::ClopperPearson(N_MET_And_SelectionCriteria, N_MET_LeptonTriggers_SelectionCriteria, level, false);

  //Turn on curves
  int NumBins = 40;
  float Weight = N_MET_LeptonTriggers_SelectionCriteria / N_SelectionCriteria;
  auto TurnOnCurveWeight{[&Weight](const float& PU){float weight = PU * Weight; return weight;}};

  auto Weighted_dataframe = d_LeptonSelection.Define("weight", TurnOnCurveWeight, {"PU"});

  auto h_LeadingLeptonPt = Weighted_dataframe.Profile1D({"h_LeadingLeptonPt", "Leading lepton p_{T}", NumBins, 0, 300}, "LeadingLeptonPt", "weight");
  auto h_SubleadingLeptonPt = Weighted_dataframe.Profile1D({"h_SubleadingLeptonPt", "Subleading lepton p_{T}", NumBins, 0, 300}, "SubleadingLeptonPt", "weight");
  auto h_LeadingLeptonEta = Weighted_dataframe.Profile1D({"h_LeadingLeptonEta", "Leading lepton #eta", NumBins, -5, 5}, "LeadingLeptonEta", "weight");
  auto h_SubleadingLeptonEta = Weighted_dataframe.Profile1D({"h_SubleadingLeptonEta", "Subleading lepton #eta", NumBins, -5, 5}, "SubleadingLeptonEta", "weight");

  auto h_LeadingVsSubleading_LeptonPt = Weighted_dataframe.Profile2D({"h_LeadingVsSubleading_LeptonPt", "Leading lepton p_{T} vs subleading lepton p_{T}", NumBins, 0, 300, NumBins, 0, 300}, "SubleadingLeptonPt", "LeadingLeptonPt", "weight");

  auto h_LeadingVsSubleading_LeptonEta = Weighted_dataframe.Profile2D({"h_LeadingVsSubleading_LeptonEta", "Leading lepton #eta vs subleading lepton #eta", NumBins, -5, 5, NumBins, -5, 5}, "SubleadingLeptonEta", "LeadingLeptonEta", "weight");

  //Errors for the turn on curves
  for (Int_t bin = 1; bin != NumBins + 1; bin++){

  	double errUp, errDown, error;

	//Leading lepton pt
	errUp = (h_LeadingLeptonPt->GetBinContent(bin) - TEfficiency::ClopperPearson(h_LeadingLeptonPt->GetBinEntries(bin), h_LeadingLeptonPt->GetBinEntries(bin) * h_LeadingLeptonPt->GetBinContent(bin), level, true));

	errDown = (h_LeadingLeptonPt->GetBinContent(bin) - TEfficiency::ClopperPearson(h_LeadingLeptonPt->GetBinEntries(bin), h_LeadingLeptonPt->GetBinEntries(bin) * h_LeadingLeptonPt->GetBinContent(bin), level, false));

	error = errUp > errDown ? errUp : errDown;

	h_LeadingLeptonPt->SetBinError(bin, error);
	
	//Subleading lepton pt
        errUp = (h_SubleadingLeptonPt->GetBinContent(bin) - TEfficiency::ClopperPearson(h_SubleadingLeptonPt->GetBinEntries(bin), h_SubleadingLeptonPt->GetBinEntries(bin) * h_SubleadingLeptonPt->GetBinContent(bin), level, true));

        errDown = (h_SubleadingLeptonPt->GetBinContent(bin) - TEfficiency::ClopperPearson(h_SubleadingLeptonPt->GetBinEntries(bin), h_SubleadingLeptonPt->GetBinEntries(bin) * h_SubleadingLeptonPt->GetBinContent(bin), level, false));

        error = errUp > errDown ? errUp : errDown;

        h_SubleadingLeptonPt->SetBinError(bin, error);

	//Leading lepton eta
        errUp = (h_LeadingLeptonEta->GetBinContent(bin) - TEfficiency::ClopperPearson(h_LeadingLeptonEta->GetBinEntries(bin), h_LeadingLeptonEta->GetBinEntries(bin) * h_LeadingLeptonEta->GetBinContent(bin), level, true));

        errDown = (h_LeadingLeptonEta->GetBinContent(bin) - TEfficiency::ClopperPearson(h_LeadingLeptonEta->GetBinEntries(bin), h_LeadingLeptonEta->GetBinEntries(bin) * h_LeadingLeptonEta->GetBinContent(bin), level, false));

        error = errUp > errDown ? errUp : errDown;

        h_LeadingLeptonEta->SetBinError(bin, error);

        //Subleading lepton eta
        errUp = (h_SubleadingLeptonEta->GetBinContent(bin) - TEfficiency::ClopperPearson(h_SubleadingLeptonEta->GetBinEntries(bin), h_SubleadingLeptonEta->GetBinEntries(bin) * h_SubleadingLeptonEta->GetBinContent(bin), level, true));

        errDown = (h_SubleadingLeptonEta->GetBinContent(bin) - TEfficiency::ClopperPearson(h_SubleadingLeptonEta->GetBinEntries(bin), h_SubleadingLeptonEta->GetBinEntries(bin) * h_SubleadingLeptonEta->GetBinContent(bin), level, false));

        error = errUp > errDown ? errUp : errDown;

        h_SubleadingLeptonEta->SetBinError(bin, error);

  }


  if((ProcessInt == 95 || ProcessInt == 96) && SystematicInt == 0){

  	std::string TurnOnCurvesOutput = "TurnOnCurves_" + Process + "_" + Systematic + "_" + Channel + "_" + NonPromptLepton + "_" +
                                    SignalRegion + "_" + SideBandRegion + "_" + ZPlusJetsControlRegion + "_" + ttbarControlRegion + "_" + Year + ".root";

  	TFile* TurnOnCurvesFile = new TFile(TurnOnCurvesOutput.c_str(), "RECREATE");

  	h_LeadingLeptonPt->GetYaxis()->SetTitle("Number of events");
  	h_SubleadingLeptonPt->GetYaxis()->SetTitle("Number of events");
  	h_LeadingLeptonEta->GetYaxis()->SetTitle("Number of events");
  	h_SubleadingLeptonEta->GetYaxis()->SetTitle("Number of events");
  	h_LeadingLeptonPt->GetXaxis()->SetTitle("p_{T}");
  	h_SubleadingLeptonPt->GetXaxis()->SetTitle("p_{T}");
  	h_LeadingLeptonEta->GetXaxis()->SetTitle("#eta");
  	h_SubleadingLeptonEta->GetXaxis()->SetTitle("#eta");
  	h_LeadingVsSubleading_LeptonPt->GetYaxis()->SetTitle("Leading lepton p_{T}");
  	h_LeadingVsSubleading_LeptonEta->GetYaxis()->SetTitle("Leading lepton #eta");
  	h_LeadingVsSubleading_LeptonPt->GetXaxis()->SetTitle("Subleading lepton p_{T}");
  	h_LeadingVsSubleading_LeptonEta->GetXaxis()->SetTitle("Subleading lepton #eta");
  	h_LeadingLeptonPt->Write();
  	h_SubleadingLeptonPt->Write();
  	h_LeadingLeptonEta->Write();
 	h_SubleadingLeptonEta->Write();
  	h_LeadingVsSubleading_LeptonPt->Write();
  	h_LeadingVsSubleading_LeptonEta->Write();

  	TurnOnCurvesFile->Close();

  }

  switch(NPLInt){
  	case 0: switch(ZPlusJetsCRInt){
			case 0: switch(ttbarCRInt){
					case 0: switch(SystematicInt){ 
  							case 0: switch(ProcessInt){
  									case 95: switch(ChannelInt){
										 	case 1: 
								
											N_SelectionCriteria_ee_MC = N_SelectionCriteria;
											N_MET_And_LeptonSelection_ee_MC = N_MET_And_SelectionCriteria;
										        N_LeptonTriggersAndSelectionCriteria_ee_MC = N_LeptonTriggersAndSelectionCriteria;
											N_MET_LeptonTriggers_SelectionCriteria_ee_MC = N_MET_LeptonTriggers_SelectionCriteria;
											Eff_ee_MC = Eff;
											Eff_MET_LeptonTriggers_SelectionCriteria_ee_MC = Eff_MET_LeptonTriggers_SelectionCriteria;
										        Eff_LeptonTriggers_SelectionCriteria_ee_MC = Eff_LeptonTriggers_SelectionCriteria;
											Eff_MET_SelectionCriteria_ee_MC = Eff_MET_SelectionCriteria;
											Alpha_ee_MC = Alpha;
											Eff_UpperUncert_ee_MC = Eff_UpperUncert;
											Eff_LowerUncert_ee_MC = Eff_LowerUncert;
											break;
		
											case 2: 

											N_SelectionCriteria_mumu_MC = N_SelectionCriteria;
                                							N_MET_And_LeptonSelection_mumu_MC = N_MET_And_SelectionCriteria;
                                							N_LeptonTriggersAndSelectionCriteria_mumu_MC = N_LeptonTriggersAndSelectionCriteria;
                               								N_MET_LeptonTriggers_SelectionCriteria_mumu_MC = N_MET_LeptonTriggers_SelectionCriteria;
                                							Eff_mumu_MC = Eff;
                                							Eff_MET_LeptonTriggers_SelectionCriteria_mumu_MC = Eff_MET_LeptonTriggers_SelectionCriteria;
                                							Eff_LeptonTriggers_SelectionCriteria_mumu_MC = Eff_LeptonTriggers_SelectionCriteria;
                                							Eff_MET_SelectionCriteria_mumu_MC = Eff_MET_SelectionCriteria;
                                							Alpha_mumu_MC = Alpha;
											Eff_UpperUncert_mumu_MC = Eff_UpperUncert;
                                							Eff_LowerUncert_mumu_MC = Eff_LowerUncert;
											break;

											case 3: 
	
											N_SelectionCriteria_emu_MC = N_SelectionCriteria;
                                							N_MET_And_LeptonSelection_emu_MC = N_MET_And_SelectionCriteria;
                                							N_LeptonTriggersAndSelectionCriteria_emu_MC = N_LeptonTriggersAndSelectionCriteria;
                                							N_MET_LeptonTriggers_SelectionCriteria_emu_MC = N_MET_LeptonTriggers_SelectionCriteria;
                                							Eff_emu_MC = Eff;
                                							Eff_MET_LeptonTriggers_SelectionCriteria_emu_MC = Eff_MET_LeptonTriggers_SelectionCriteria;
                                							Eff_LeptonTriggers_SelectionCriteria_emu_MC = Eff_LeptonTriggers_SelectionCriteria;
                                							Eff_MET_SelectionCriteria_emu_MC = Eff_MET_SelectionCriteria;
                                							Alpha_emu_MC = Alpha;
											Eff_UpperUncert_emu_MC = Eff_UpperUncert;
                                							Eff_LowerUncert_emu_MC = Eff_LowerUncert;
											break;
			
											default: std::cout << "ERROR: ChannelInt must be 1 (for ee), 2 (for mumu) or 3 (for emu)." << std::endl; break;

									}


									case 96: switch(ChannelInt){
                        								case 1: 

											N_SelectionCriteria_ee_DATA = N_SelectionCriteria;
                                							N_MET_And_LeptonSelection_ee_DATA = N_MET_And_SelectionCriteria;
                                							N_LeptonTriggersAndSelectionCriteria_ee_DATA = N_LeptonTriggersAndSelectionCriteria;
                                							N_MET_LeptonTriggers_SelectionCriteria_ee_DATA = N_MET_LeptonTriggers_SelectionCriteria;
                                							Eff_ee_DATA = Eff;
                                							Eff_MET_LeptonTriggers_SelectionCriteria_ee_DATA = Eff_MET_LeptonTriggers_SelectionCriteria;
                                							Eff_LeptonTriggers_SelectionCriteria_ee_DATA = Eff_LeptonTriggers_SelectionCriteria;
                                							Eff_MET_SelectionCriteria_ee_DATA = Eff_MET_SelectionCriteria;
                                							Alpha_ee_DATA = Alpha;
											Eff_UpperUncert_ee_DATA = Eff_UpperUncert;
                                							Eff_LowerUncert_ee_DATA = Eff_LowerUncert;
                                							break;

                        								case 2: 

											N_SelectionCriteria_mumu_DATA = N_SelectionCriteria;
                                							N_MET_And_LeptonSelection_mumu_DATA = N_MET_And_SelectionCriteria;
                                							N_LeptonTriggersAndSelectionCriteria_mumu_DATA = N_LeptonTriggersAndSelectionCriteria;
                                							N_MET_LeptonTriggers_SelectionCriteria_mumu_DATA = N_MET_LeptonTriggers_SelectionCriteria;
                                							Eff_mumu_DATA = Eff;
                                						Eff_MET_LeptonTriggers_SelectionCriteria_mumu_DATA = Eff_MET_LeptonTriggers_SelectionCriteria;
                                							Eff_LeptonTriggers_SelectionCriteria_mumu_DATA = Eff_LeptonTriggers_SelectionCriteria;
                                							Eff_MET_SelectionCriteria_mumu_DATA = Eff_MET_SelectionCriteria;
                               								Alpha_mumu_DATA = Alpha;
											Eff_UpperUncert_mumu_DATA = Eff_UpperUncert;
                                							Eff_LowerUncert_mumu_DATA = Eff_LowerUncert;
                                							break;

                        								case 3: 

											N_SelectionCriteria_emu_DATA = N_SelectionCriteria;
                                							N_MET_And_LeptonSelection_emu_DATA = N_MET_And_SelectionCriteria;
                                							N_LeptonTriggersAndSelectionCriteria_emu_DATA = N_LeptonTriggersAndSelectionCriteria;
                                							N_MET_LeptonTriggers_SelectionCriteria_emu_DATA = N_MET_LeptonTriggers_SelectionCriteria;
                                							Eff_emu_DATA = Eff;
                                							Eff_MET_LeptonTriggers_SelectionCriteria_emu_DATA = Eff_MET_LeptonTriggers_SelectionCriteria;
                                							Eff_LeptonTriggers_SelectionCriteria_emu_DATA = Eff_LeptonTriggers_SelectionCriteria;
                                							Eff_MET_SelectionCriteria_emu_DATA = Eff_MET_SelectionCriteria;
                                							Alpha_emu_DATA = Alpha;
											Eff_UpperUncert_emu_DATA = Eff_UpperUncert;
                                							Eff_LowerUncert_emu_DATA = Eff_LowerUncert;
                                							break;

                        								default: std::cout << "ERROR: ChannelInt must be 1 (for ee), 2 (for mumu) or 3 (for emu)." << std::endl; break;

                								}

									default: break;

 							 }
					
						default:break;
					
					}

				default: break;


			}

		default: break;


	}

    default: break;

  }

  //Calculating the trigger scale factors
  SF_ee = Eff_ee_DATA/(Eff_ee_MC + 1.0e-06);
  SF_mumu = Eff_ee_DATA/(Eff_mumu_MC + 1.0e-06);
  SF_emu = Eff_ee_DATA/(Eff_emu_MC + 1.0e-06);

  //Uncertainties in the trigger scale factors
  float SF_UpperUncert_ee = ((Eff_ee_DATA + Eff_UpperUncert_ee_DATA) / (Eff_ee_MC - Eff_LowerUncert_ee_MC + 1.0e-06)) - SF_ee;
  float SF_LowerUncert_ee = ((Eff_ee_DATA + Eff_LowerUncert_ee_DATA)/ (Eff_ee_MC - Eff_UpperUncert_ee_MC + 1.0e-06)) - SF_ee;
  float SF_UpperUncert_mumu = ((Eff_mumu_DATA + Eff_UpperUncert_mumu_DATA) / (Eff_mumu_MC - Eff_LowerUncert_mumu_MC + 1.0e-06)) - SF_mumu;
  float SF_LowerUncert_mumu = ((Eff_mumu_DATA + Eff_LowerUncert_mumu_DATA)/ (Eff_mumu_MC - Eff_UpperUncert_mumu_MC + 1.0e-06)) - SF_mumu;
  float SF_UpperUncert_emu = ((Eff_emu_DATA + Eff_UpperUncert_emu_DATA) / (Eff_emu_MC - Eff_LowerUncert_emu_MC + 1.0e-06)) - SF_emu;
  float SF_LowerUncert_emu = ((Eff_emu_DATA + Eff_LowerUncert_emu_DATA)/ (Eff_emu_MC - Eff_UpperUncert_emu_MC + 1.0e-06)) - SF_emu;


  SF_Uncert_ee = 0.0;
  if (SF_UpperUncert_ee > SF_LowerUncert_ee){SF_Uncert_ee = SF_UpperUncert_ee;}
  else{SF_Uncert_ee = SF_LowerUncert_ee;}

  SF_Uncert_mumu = 0.0;
  if (SF_UpperUncert_mumu > SF_LowerUncert_mumu){SF_Uncert_mumu = SF_UpperUncert_mumu;}
  else{SF_Uncert_mumu = SF_LowerUncert_mumu;}

  SF_Uncert_emu = 0.0;
  if (SF_UpperUncert_emu > SF_LowerUncert_emu){SF_Uncert_emu = SF_UpperUncert_emu;}
  else{SF_Uncert_emu = SF_LowerUncert_emu;}

  std::string TriggerSFValuesFile = "TriggerSFValues_" + Process + "_" + Systematic + "_" + Channel + "_" + NonPromptLepton + "_" +
                                     SignalRegion + "_" + SideBandRegion + "_" + ZPlusJetsControlRegion + "_" + ttbarControlRegion + "_" + Year + ".txt";

  std::ofstream TriggerSFValues;
  TriggerSFValues.open(TriggerSFValuesFile.c_str());

  TriggerSFValues << "N_SelectionCriteria_ee_MC = " << N_SelectionCriteria_ee_MC << '\n'
		  << "N_SelectionCriteria_mumu_MC = " << N_SelectionCriteria_mumu_MC << '\n'
		  << "N_SelectionCriteria_emu_MC = " << N_SelectionCriteria_emu_MC << '\n'
		  << "N_SelectionCriteria_ee_DATA = " << N_SelectionCriteria_ee_DATA << '\n'
                  << "N_SelectionCriteria_mumu_DATA = " << N_SelectionCriteria_mumu_DATA << '\n'
                  << "N_SelectionCriteria_emu_DATA = " << N_SelectionCriteria_emu_DATA << '\n'
		  << "N_MET_And_LeptonSelection_ee_MC = " << N_MET_And_LeptonSelection_ee_MC << '\n'
		  << "N_MET_And_LeptonSelection_mumu_MC = " << N_MET_And_LeptonSelection_mumu_MC << '\n'
		  << "N_MET_And_LeptonSelection_emu_MC = " << N_MET_And_LeptonSelection_emu_MC << '\n'
		  << "N_MET_And_LeptonSelection_ee_DATA = " << N_MET_And_LeptonSelection_ee_DATA << '\n'
                  << "N_MET_And_LeptonSelection_mumu_DATA = " << N_MET_And_LeptonSelection_mumu_DATA << '\n'
                  << "N_MET_And_LeptonSelection_emu_DATA = " << N_MET_And_LeptonSelection_emu_DATA << '\n'
		  << "N_MET_LeptonTriggers_SelectionCriteria_ee_MC = " << N_MET_LeptonTriggers_SelectionCriteria_ee_MC << '\n'
		  << "N_MET_LeptonTriggers_SelectionCriteria_mumu_MC = " << N_MET_LeptonTriggers_SelectionCriteria_mumu_MC << '\n'
		  << "N_MET_LeptonTriggers_SelectionCriteria_emu_MC = " << N_MET_LeptonTriggers_SelectionCriteria_emu_MC << '\n'
		  << "N_MET_LeptonTriggers_SelectionCriteria_ee_DATA = " << N_MET_LeptonTriggers_SelectionCriteria_ee_DATA << '\n'
                  << "N_MET_LeptonTriggers_SelectionCriteria_mumu_DATA = " << N_MET_LeptonTriggers_SelectionCriteria_mumu_DATA << '\n'
                  << "N_MET_LeptonTriggers_SelectionCriteria_emu_DATA = " << N_MET_LeptonTriggers_SelectionCriteria_emu_DATA << '\n'
		  << "Eff_ee_MC = " << Eff_ee_MC << '\n'
		  << "Eff_mumu_MC = " << Eff_mumu_MC << '\n'
		  << "Eff_emu_MC = " << Eff_emu_MC << '\n'
		  << "Eff_ee_DATA = " << Eff_ee_DATA << '\n'
                  << "Eff_mumu_DATA = " << Eff_mumu_DATA << '\n'
                  << "Eff_emu_DATA = " << Eff_emu_DATA << '\n'
                  << "Alpha_ee_MC = " << Alpha_ee_MC << '\n'
		  << "Alpha_mumu_MC = " << Alpha_mumu_MC << '\n'
		  << "Alpha_emu_MC = " << Alpha_emu_MC << '\n'
		  << "Alpha_ee_DATA = " << Alpha_ee_DATA << '\n'
                  << "Alpha_mumu_DATA = " << Alpha_mumu_DATA << '\n'
                  << "Alpha_emu_DATA = " << Alpha_emu_DATA << '\n'                               
                  << "Eff_UpperUncert_ee_MC = " << Eff_UpperUncert_ee_MC << '\n'
		  << "Eff_UpperUncert_mumu_MC = " << Eff_UpperUncert_mumu_MC << '\n'
		  << "Eff_UpperUncert_emu_MC = " << Eff_UpperUncert_emu_MC << '\n'
		  << "Eff_LowerUncert_ee_MC = " << Eff_LowerUncert_ee_MC << '\n'
                  << "Eff_LowerUncert_mumu_MC = " << Eff_LowerUncert_mumu_MC << '\n'
                  << "Eff_LowerUncert_emu_MC = " << Eff_LowerUncert_emu_MC << '\n'
		  << "Eff_UpperUncert_ee_DATA = " << Eff_UpperUncert_ee_DATA << '\n'
                  << "Eff_UpperUncert_mumu_DATA = " << Eff_UpperUncert_mumu_DATA << '\n'
                  << "Eff_UpperUncert_emu_DATA = " << Eff_UpperUncert_emu_DATA << '\n'
                  << "Eff_LowerUncert_ee_DATA = " << Eff_LowerUncert_ee_DATA << '\n'
                  << "Eff_LowerUncert_mumu_DATA = " << Eff_LowerUncert_mumu_DATA << '\n'
                  << "Eff_LowerUncert_emu_DATA = " << Eff_LowerUncert_emu_DATA << '\n' 
		  << "SF_ee = " << SF_ee << '\n'
		  << "SF_mumu = " << SF_mumu << '\n'
		  << "SF_emu = " << SF_emu << '\n'
		  << "SF_Uncert_ee = " << SF_Uncert_ee << '\n'
                  << "SF_Uncert_mumu = " << SF_Uncert_mumu << '\n'
                  << "SF_Uncert_emu = " << SF_Uncert_emu << '\n'
   		  << std::endl;

  if(ProcessInt == 95 || ProcessInt == 96){return;}


  //Z boson candidate reconstruction
  auto d_ZCandidateReco = d_LeptonSelection.Define("LeptonGenPartFlav", LeptonVariableFunctionChars, {"Electron_genPartFlav", "Muon_genPartFlav"})
					   .Define("TightLeptonsGenPartFlav", select<chars>, {"LeptonGenPartFlav", "TightLeptons"})
				           .Define("OppositeSignNonPrompt", OppositeSignNonPrompt, {"TightLeptonsCharge", "TightLeptonsGenPartFlav"})
                                           .Define("OppositeSignPrompt", OppositeSignPrompt, {"TightLeptonsCharge", "TightLeptonsGenPartFlav"})
                                           .Define("SameSignNonPrompt", SameSignNonPrompt, {"TightLeptonsCharge", "TightLeptonsGenPartFlav"})
                                           .Define("SameSignPrompt", SameSignPrompt, {"TightLeptonsCharge", "TightLeptonsGenPartFlav"})
                                           .Define("RecoZ", RecoZ, {"LeadingLeptonPt", "LeadingLeptonEta", "LeadingLeptonPhi", "LeadingLeptonMass",
								    "SubleadingLeptonPt", "SubleadingLeptonEta", "SubleadingLeptonPhi", "SubleadingLeptonMass"})
                                           .Define("RecoZPt", TLorentzVectorVariablePt, {"RecoZ"})
                                           .Define("RecoZPhi", TLorentzVectorVariablePhi, {"RecoZ"})
                                           .Define("RecoZEta", TLorentzVectorVariableEta, {"RecoZ"})
					   .Define("RecoZMass", TLorentzVectorVariableMass, {"RecoZ"})
                                           .Define("dR_ll", deltaRcheck_float, {"LeadingLeptonEta", "LeadingLeptonPhi", "SubleadingLeptonEta", "SubleadingLeptonPhi"})
                                           .Define("dPhi_ll", DeltaPhi_floatandfloat, {"LeadingLeptonPhi", "SubleadingLeptonPhi"})
                                           .Define("LeptonFourMomentum", LeptonFourMomentumFunction, {"TightLeptonsPt", "TightLeptonsEta", "TightLeptonsPhi", "TightLeptonsMass"})
                                           .Define("RochCorrVec", RochCorrVec_Function, {"TightLeptonsCharge", "TightLeptonsPt", "TightLeptonsEta", "TightLeptonsPhi", 
										         "Muon_genPartIdx", "Muon_nTrackerLayers"})
                                           .Define("MuonFourMomentum_RochCorr", RochCorrMuon4Mo, {"LeptonFourMomentum", "RochCorrVec"})
                                           .Define("LeptonPt_RochCorr", TLorentzVector_float_pt, {"MuonFourMomentum_RochCorr"})
                                           .Define("LeptonEta_RochCorr", TLorentzVector_float_eta, {"MuonFourMomentum_RochCorr"})
                                           .Define("LeptonPhi_RochCorr", TLorentzVector_float_phi, {"MuonFourMomentum_RochCorr"})
                                           .Define("LeptonMass_RochCorr", TLorentzVector_float_mass, {"MuonFourMomentum_RochCorr"})
                                           .Define("z_mass", inv_mass, {"LeptonPt_RochCorr", "LeptonEta_RochCorr", "LeptonPhi_RochCorr", "LeptonMass_RochCorr"})
                                           .Filter(z_mass_cut, {"z_mass"}, "Z mass cut");


  std::string ZCandidateRecoFile = "ZCandidateReco_" + Process + "_" + Systematic + "_" + Channel + "_" + NonPromptLepton + "_" +
                                    SignalRegion + "_" + SideBandRegion + "_" + ZPlusJetsControlRegion + "_" + ttbarControlRegion + "_" + Year + ".root";
  
  auto Snapshot_ZCandidateReco = d_ZCandidateReco.Snapshot("Events", ZCandidateRecoFile.c_str());

  //Jet selection
  auto d_JetSelection = d_ZCandidateReco.Define("sJER_Nominal", SJER_Nominal_Function, {"Jet_eta", "fixedGridRhoFastjetAll", "Jet_pt"})
                                        .Define("sJER_up", SJER_Up_Function, {"Jet_eta", "fixedGridRhoFastjetAll", "Jet_pt"})
                                        .Define("sJER_down", SJER_Down_Function, {"Jet_eta", "fixedGridRhoFastjetAll", "Jet_pt"})
                      		        .Define("sigma_JER", sigma_JER, {"Jet_eta", "fixedGridRhoFastjetAll", "Jet_pt"})
			                .Define("sigma_JER_up", sigma_JER_up, {"Jet_eta", "fixedGridRhoFastjetAll", "Jet_pt"})
					.Define("sigma_JER_down", sigma_JER_down, {"Jet_eta", "fixedGridRhoFastjetAll", "Jet_pt"})
                      			.Define("cJER", JetSmearingFunction_HybridMethod, {{"Jet_pt", "Jet_eta", "Jet_phi", "GenJet_pt", "GenJet_eta", "GenJet_phi", SJER, SIGMAJER, "Jet_genJetIdx"}})
                      		        .Define("SmearedJet4Momentum", ApplyCJER, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "cJER", "nJet"})
                      			.Define("SmearedJetPt", GetSmearedJetPt, {"SmearedJet4Momentum", "Jet_pt"})
                      			.Define("SmearedJetPhi", GetSmearedJetPhi, {"SmearedJet4Momentum", "Jet_phi"})
                      			.Define("SmearedJetEta", GetSmearedJetEta, {"SmearedJet4Momentum", "Jet_eta"})
                      			.Define("SmearedJetMass", GetSmearedJetMass, {"SmearedJet4Momentum", "Jet_mass"})
		      		        .Define("LeadingJetMass", LeadingVariable, {JetMassInput})
                                        .Define("SubleadingJetMass", SubleadingVariable, {JetMassInput})
                                        .Define("ThirdJetMass", ThirdLeadingVariable, {JetMassInput})
                                        .Define("FourthJetMass", FourthLeadingVariable, {JetMassInput})
                                        .Define("LeadingJetPt", LeadingVariable, {JetPtInput})
                                        .Define("SubleadingJetPt", SubleadingVariable, {JetPtInput})
                                        .Define("ThirdJetPt", ThirdLeadingVariable, {JetPtInput})
                                        .Define("FourthJetPt", FourthLeadingVariable, {JetPtInput})
                                        .Define("SumSquaredPt", SumSquared2LeadingJets_pT, {"LeadingJetPt", "SubleadingJetPt"})
                                        .Define("JetPtSum", JetSum, {"LeadingJetPt", "SubleadingJetPt", "ThirdJetPt", "FourthJetPt"})
                                        .Define("LeadingJetEta", LeadingVariable, {JetEtaInput})
                                        .Define("SubleadingJetEta", SubleadingVariable, {JetEtaInput})
                                        .Define("ThirdJetEta", ThirdLeadingVariable, {JetEtaInput})
                                        .Define("FourthJetEta", FourthLeadingVariable, {JetEtaInput})
                                        .Define("LeadingJetPhi", LeadingVariable, {JetPhiInput})
                                        .Define("SubleadingJetPhi", SubleadingVariable, {JetPhiInput})
                                        .Define("ThirdJetPhi", ThirdLeadingVariable, {JetPhiInput})
                                        .Define("FourthJetPhi", FourthLeadingVariable, {JetPhiInput})
                                        .Define("dRJet_Lepton", deltaRcheck_floats, {JetEtaInput, JetPhiInput, "LeptonEta_RochCorr", "LeptonPhi_RochCorr"})
                                        .Define("dR_j1j2", deltaRcheck_float, {"LeadingJetEta", "LeadingJetPhi", "SubleadingJetEta", "SubleadingJetPhi"})
                                        .Define("dPhi_j1j2", DeltaPhi_floatandfloat, {"LeadingJetPhi", "SubleadingJetPhi"})
                                        .Define("LeadingJetHT", HT, {"LeadingJetPt"})
                                        .Define("SubleadingJetHT", HT, {"SubleadingJetPt"})
                                        .Define("ThirdJetHT", HT, {"ThirdJetPt"})
                                        .Define("FourthJetHT", HT, {"FourthJetPt"})
                                        .Define("TotJetHT", TotJetHT, {"LeadingJetHT", "SubleadingJetHT", "ThirdJetHT", "FourthJetHT"})
                                        .Define("LeadingLeptonHT", HT, {"LeadingLeptonPt"})
                                        .Define("SubleadingLeptonHT", HT, {"SubleadingLeptonPt"})
                                        .Define("TotLepHT", TotLepHT, {"LeadingLeptonHT", "SubleadingLeptonHT"})
                                        .Define("TotHTOverTotpT_Jets", TotHTOverTotpT, {"TotJetHT", "JetPtSum"})
                                        .Define("LepPtSum", LepSum, {"LeadingLeptonPt", "SubleadingLeptonPt"})
                                        .Define("LepEtaSum", LepSum, {"LeadingLeptonPt", "SubleadingLeptonPt"})
                                        .Define("LepPhiSum", LepSum, {"LeadingLeptonPt", "SubleadingLeptonPt"})
                                        .Define("TotHTOverTotpT_Leptons", TotHTOverTotpT, {"TotLepHT", "LepPtSum"})
                                        .Define("InvMassAllJets", InvMass_AllJets, {"LeadingJetPt",   "SubleadingJetPt",   "ThirdJetPt",   "FourthJetPt",
										    "LeadingJetEta",  "SubleadingJetEta",  "ThirdJetEta",  "FourthJetEta",
										    "LeadingJetPhi",  "SubleadingJetPhi",  "ThirdJetPhi",  "FourthJetPhi",
										    "LeadingJetMass", "SubleadingJetMass", "ThirdJetMass", "FourthJetMass"})
                                        .Define("InvMass3Jets", InvMass_3Jets, {"LeadingJetPt",   "SubleadingJetPt",   "ThirdJetPt",
										"LeadingJetEta",  "SubleadingJetEta",  "ThirdJetEta",
									        "LeadingJetPhi",  "SubleadingJetPhi",  "ThirdJetPhi",
									        "LeadingJetMass", "SubleadingJetMass", "ThirdJetMass"})
                                        .Define("JetEtaSum", JetSum, {"LeadingJetEta", "SubleadingJetEta", "ThirdJetEta", "FourthJetEta"})
                                        .Define("JetPhiSum", JetSum, {"LeadingJetPhi", "SubleadingJetPhi", "ThirdJetPhi", "FourthJetPhi"})
                                        .Define("tight_jets", tight_jets_function, {JetPtInput, JetEtaInput, "Jet_jetId", "dRJet_Lepton"})
                                        .Filter(jet_selection_function, {"tight_jets"}, "jet cut");

  std::string JetSelectionFile = "JetSelection_" + Process + "_" + Systematic + "_" + Channel + "_" + NonPromptLepton + "_" +
                                 SignalRegion + "_" + SideBandRegion + "_" + ZPlusJetsControlRegion + "_" + ttbarControlRegion + "_" + Year + ".root";
  
  auto Snapshot_JetSelection = d_JetSelection.Snapshot("Events", JetSelectionFile.c_str());


  //B jet selection
  auto d_BJetSelection = d_JetSelection.Define("bjets", bjet_id, {"tight_jets", "Jet_btagCSVV2", JetEtaInput})
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
			 

  std::string BJetSelectionFile = "BJetSelection_" + Process + "_" + Systematic + "_" + Channel + "_" + NonPromptLepton + "_" +
                                 SignalRegion + "_" + SideBandRegion + "_" + ZPlusJetsControlRegion + "_" + ttbarControlRegion + "_" + Year + ".root";

  auto Snapshot_BJetSelection = d_BJetSelection.Snapshot("Events", BJetSelectionFile.c_str());


  //For the b-tagging efficiencies
  std::string BTagString = "BTagEffPlots_" + Process + "_" + Systematic + "_" + Channel + "_" + NonPromptLepton + "_" +
                            SignalRegion + "_" + SideBandRegion + "_" + ZPlusJetsControlRegion + "_" + ttbarControlRegion + "_" + Year + ".root";

  TFile* BTagEffPlots = new TFile(BTagString.c_str(), "RECREATE");
  double minpt = 0;
  double maxpt = 500;
  double mineta = -3;
  double maxeta = 3;

  h_bjet_num = d_BJetSelection.Histo2D({"h_bjet_ee_num", "h_bjet_ee_num", NumBins, mineta, maxeta, NumBins, minpt, maxpt}, {"BTAGEFF_bjet_eta_num"}, {"BTAGEFF_bjet_pt_num"});
  
  h_nonbjet_num = d_BJetSelection.Histo2D({"h_nonbjet_ee_num", "h_nonbjet_ee_num", NumBins, mineta, maxeta, NumBins, minpt, maxpt}, {"BTAGEFF_nonbjet_eta_num"}, {"BTAGEFF_nonbjet_pt_num"});

  h_charm_num = d_BJetSelection.Histo2D({"h_charm_ee_num", "h_charm_ee_num", NumBins, mineta, maxpt, NumBins, minpt, maxpt}, {"BTAGEFF_charm_eta_num"}, {"BTAGEFF_charm_pt_num"});
  
  h_lightjets_num = d_BJetSelection.Histo2D({"h_lightjets_ee_num", "h_lightjets_ee_num", NumBins, mineta, maxeta, NumBins, minpt, maxpt}, {"BTAGEFF_lightjets_eta_num"}, {"BTAGEFF_lightjets_pt_num"});
  
  h_gluon_num = d_BJetSelection.Histo2D({"h_gluon_ee_num", "h_gluon_ee_num", NumBins, mineta, maxeta, NumBins, minpt, maxpt}, {"BTAGEFF_gluon_eta_num"}, {"BTAGEFF_gluon_pt_num"});


  h_bjet_denom = d_BJetSelection.Histo2D({"h_bjet_ee_denom", "h_bjet_ee_denom", NumBins, mineta, maxeta, NumBins, minpt, maxpt}, {"BTAGEFF_bjet_eta_denom"}, {"BTAGEFF_bjet_pt_denom"});
   
  h_nonbjet_denom = d_BJetSelection.Histo2D({"h_nonbjet_ee_denom", "h_nonbjet_ee_denom", NumBins, mineta, maxeta, NumBins, minpt, maxpt}, {"BTAGEFF_nonbjet_eta_denom"}, {"BTAGEFF_nonbjet_pt_denom"});

  h_charm_denom = d_BJetSelection.Histo2D({"h_charm_ee_denom", "h_charm_ee_denom", NumBins, mineta, maxeta, NumBins, minpt, maxpt}, {"BTAGEFF_charm_eta_denom"}, {"BTAGEFF_charm_pt_denom"});
   
  h_lightjets_denom = d_BJetSelection.Histo2D({"h_lightjets_ee_denom", "h_lightjets_ee_denom", NumBins, mineta, maxeta, NumBins, minpt, maxpt}, {"BTAGEFF_lightjets_eta_denom"}, {"BTAGEFF_lightjets_pt_denom"});
   
  h_gluon_denom = d_BJetSelection.Histo2D({"h_gluon_ee_denom", "h_gluon_ee_denom", NumBins, mineta, maxeta, NumBins, minpt, maxpt}, {"BTAGEFF_gluon_eta_denom"}, {"BTAGEFF_gluon_pt_denom"});
 
  h_bjet_num->GetXaxis()->SetTitle("#eta");
  h_nonbjet_num->GetXaxis()->SetTitle("#eta");
  h_charm_num->GetXaxis()->SetTitle("#eta");
  h_lightjets_num->GetXaxis()->SetTitle("#eta");
  h_gluon_num->GetXaxis()->SetTitle("#eta");
  h_bjet_denom->GetXaxis()->SetTitle("#eta");
  h_nonbjet_denom->GetXaxis()->SetTitle("#eta");
  h_charm_denom->GetXaxis()->SetTitle("#eta");
  h_lightjets_denom->GetXaxis()->SetTitle("#eta");
  h_gluon_denom->GetXaxis()->SetTitle("#eta");

  h_bjet_num->GetYaxis()->SetTitle("p_{T}");
  h_nonbjet_num->GetYaxis()->SetTitle("p_{T}");
  h_charm_num->GetYaxis()->SetTitle("p_{T}");
  h_lightjets_num->GetYaxis()->SetTitle("p_{T}");
  h_gluon_num->GetYaxis()->SetTitle("p_{T}");

  h_bjet_denom->GetYaxis()->SetTitle("p_{T}");
  h_nonbjet_denom->GetYaxis()->SetTitle("p_{T}");
  h_charm_denom->GetYaxis()->SetTitle("p_{T}");
  h_lightjets_denom->GetYaxis()->SetTitle("p_{T}");
  h_gluon_denom->GetYaxis()->SetTitle("p_{T}");

  h_bjet_num->Write();
  h_nonbjet_num->Write();
  h_charm_num->Write();
  h_lightjets_num->Write();
  h_gluon_num->Write();
  h_bjet_denom->Write();
  h_nonbjet_denom->Write();
  h_charm_denom->Write();
  h_lightjets_denom->Write();
  h_gluon_denom->Write();

  TH2D* h_bjet = new TH2D("h_bjet", "h_bjet", NumBins, mineta, maxeta, NumBins, minpt, maxpt);
  TH2D* h_nonbjet = new TH2D("h_nonbjet", "h_nonbjet", NumBins, mineta, maxeta, NumBins, minpt, maxpt);
  TH2D* h_charm = new TH2D("h_charm", "h_charm", NumBins, mineta, maxeta, NumBins, minpt, maxpt);
  TH2D* h_lightjets = new TH2D("h_lightjets", "h_lightjets", NumBins, mineta, maxeta, NumBins, minpt, maxpt);
  TH2D* h_gluon = new TH2D("h_gluon", "h_gluon", NumBins, mineta, maxeta, NumBins, minpt, maxpt);

  h_bjet = dynamic_cast<TH2D*>(h_bjet_num->Clone());
  h_bjet->Divide(h_bjet_denom.GetPtr());
  h_charm = dynamic_cast<TH2D*>(h_charm_num->Clone());
  h_charm->Divide(h_charm_denom.GetPtr());
  h_lightjets = dynamic_cast<TH2D*>(h_lightjets_num->Clone());
  h_lightjets->Divide(h_lightjets_denom.GetPtr());
  h_gluon = dynamic_cast<TH2D*>(h_gluon_num->Clone());
  h_gluon->Divide(h_gluon_denom.GetPtr());
  h_nonbjet = dynamic_cast<TH2D*>(h_nonbjet_num->Clone());
  h_nonbjet->Divide(h_nonbjet_denom.GetPtr());

  h_bjet->SetTitle("h_bjet");
  h_nonbjet->SetTitle("h_nonbjet");
  h_charm->SetTitle("h_charm");
  h_lightjets->SetTitle("h_lightjets");
  h_gluon->SetTitle("h_gluon");
  h_bjet->GetXaxis()->SetTitle("#eta");
  h_nonbjet->GetXaxis()->SetTitle("#eta");
  h_charm->GetXaxis()->SetTitle("#eta");
  h_lightjets->GetXaxis()->SetTitle("#eta");
  h_gluon->GetXaxis()->SetTitle("#eta");
  h_bjet->GetYaxis()->SetTitle("p_{T}");
  h_nonbjet->GetYaxis()->SetTitle("p_{T}");
  h_charm->GetYaxis()->SetTitle("p_{T}");
  h_lightjets->GetYaxis()->SetTitle("p_{T}");
  h_gluon->GetYaxis()->SetTitle("p_{T}");

  h_bjet->Write();
  h_nonbjet->Write();
  h_charm->Write();
  h_lightjets->Write();
  h_gluon->Write();

  BTagEffPlots->Close();


  //Reconstrucing the W boson candidate
  auto d_WCandReco = d_BJetSelection.Define("lead_bjet", find_lead_mask, {"bjets", "SmearedJetPt"})
                 		    .Define("w_reco_jets", WPair, {"SmearedJetPt", "SmearedJetPhi", "SmearedJetEta", "SmearedJetMass", "Jet_jetId", "lead_bjet"})
                 		    .Define("w_pair_pt", select<floats>, {"SmearedJetPt", "w_reco_jets"})
                 		    .Define("w_pair_eta", select<floats>, {"SmearedJetEta", "w_reco_jets"})
                 		    .Define("w_pair_phi", select<floats>, {"SmearedJetPhi", "w_reco_jets"})
                 		    .Define("w_pair_mass", select<floats>, {"SmearedJetMass", "w_reco_jets"})
                 		    .Define("w_mass", inv_mass, {"w_pair_pt", "w_pair_eta", "w_pair_phi", "w_pair_mass"})
				    .Define("WPairJet1", WPairJet1, {"SmearedJetPt", "SmearedJetPhi", "SmearedJetEta", "SmearedJetMass", "Jet_jetId", "lead_bjet"})
                 		    .Define("WPairJet2", WPairJet2, {"SmearedJetPt", "SmearedJetPhi", "SmearedJetEta", "SmearedJetMass", "Jet_jetId", "lead_bjet"})
				    .Define("WPairJet1Pt", TLorentzVectorVariablePt, {"WPairJet1"})
				    .Define("WPairJet1Eta", TLorentzVectorVariableEta, {"WPairJet1"})
				    .Define("WPairJet1Phi", TLorentzVectorVariablePhi, {"WPairJet1"})
				    .Define("WPairJet1Mass", TLorentzVectorVariableMass, {"WPairJet1"})
				    .Define("WPairJet2Pt", TLorentzVectorVariablePt, {"WPairJet2"})
                                    .Define("WPairJet2Eta", TLorentzVectorVariableEta, {"WPairJet2"})
                                    .Define("WPairJet2Phi", TLorentzVectorVariablePhi, {"WPairJet2"})
                                    .Define("WPairJet2Mass", TLorentzVectorVariableMass, {"WPairJet2"})
				    .Define("dR_WJet1_WJet2", deltaRcheck_W_function, {"WPairJet1Phi", "WPairJet1Eta", "WPairJet2Phi", "WPairJet2Eta"})
				    .Define("dWj1j2", DeltaPhi_function2, {"WPairJet1Phi", "WPairJet2Phi"})
				    .Define("dR_WJet1_LeadingLepton", deltaRcheck_W_function2, {"WPairJet1Phi", "WPairJet1Eta", "LeadingLeptonPhi", "LeadingLeptonEta"})
				    .Define("dR_WJet1_SubleadingLepton", deltaRcheck_W_function2, {"WPairJet1Phi", "WPairJet1Eta", "SubleadingLeptonPhi", "SubleadingLeptonEta"})
				    .Define("dR_WJet2_LeadingLepton", deltaRcheck_W_function2, {"WPairJet2Phi", "WPairJet2Eta", "LeadingLeptonPhi", "LeadingLeptonEta"})
                                    .Define("dR_WJet2_SubleadingLepton", deltaRcheck_W_function2, {"WPairJet2Phi", "WPairJet2Eta", "SubleadingLeptonPhi", "SubleadingLeptonEta"})
				    .Define("dR_WJet1_LeadingJet", deltaRcheck_W_function2, {"WPairJet1Phi", "WPairJet1Eta", "LeadingJetPhi", "LeadingJetEta"})
                                    .Define("dR_WJet1_SubleadingJet", deltaRcheck_W_function2, {"WPairJet1Phi", "WPairJet1Eta", "SubleadingJetPhi", "SubleadingJetEta"})
                                    .Define("dR_WJet1_ThirdJet", deltaRcheck_W_function2, {"WPairJet1Phi", "WPairJet1Eta", "ThirdJetPhi", "ThirdJetEta"})
                                    .Define("dR_WJet1_FourthJet", deltaRcheck_W_function2, {"WPairJet1Phi", "WPairJet1Eta", "FourthJetPhi", "FourthJetEta"})
			            .Define("dR_WJet2_LeadingJet", deltaRcheck_W_function2, {"WPairJet2Phi", "WPairJet2Eta", "LeadingJetPhi", "LeadingJetEta"})
                                    .Define("dR_WJet2_SubleadingJet", deltaRcheck_W_function2, {"WPairJet2Phi", "WPairJet2Eta", "SubleadingJetPhi", "SubleadingJetEta"})
                                    .Define("dR_WJet2_ThirdJet", deltaRcheck_W_function2, {"WPairJet2Phi", "WPairJet2Eta", "ThirdJetPhi", "ThirdJetEta"})
                                    .Define("dR_WJet2_FourthJet", deltaRcheck_W_function2, {"WPairJet2Phi", "WPairJet2Eta", "FourthJetPhi", "FourthJetEta"})
				    .Define("dPhi_WJet1_LeadingLepton", DeltaPhi_doublesandfloat, {"WPairJet1Phi", "LeadingLeptonPhi"})
                                    .Define("dPhi_WJet1_SubleadingLepton", DeltaPhi_doublesandfloat, {"WPairJet1Phi", "SubleadingLeptonPhi"})
                                    .Define("dPhi_WJet2_LeadingLepton", DeltaPhi_doublesandfloat, {"WPairJet2Phi", "LeadingLeptonPhi"})
                                    .Define("dPhi_WJet2_SubleadingLepton", DeltaPhi_doublesandfloat, {"WPairJet2Phi", "SubleadingLeptonPhi"})
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
				    .Define("mtW", TransverseWMass, {"dPhi_j1j2", "WPairJet1Pt", "WPairJet2Pt"})
				    .Filter(w_mass_cut, {"w_mass", "MET_sumEt"}, "W mass cut");

   
  std::string WCandRecoFile = "WCandReco_" + Process + "_" + Systematic + "_" + Channel + "_" + NonPromptLepton + "_" +
                              SignalRegion + "_" + SideBandRegion + "_" + ZPlusJetsControlRegion + "_" + ttbarControlRegion + "_" + Year + ".root";

  auto Snapshot_WCandReco = d_WCandReco.Snapshot("Events", WCandRecoFile.c_str());


  //Reconstructing the top quark candidate
  auto d_TopCandReco = d_WCandReco.Define("RecoW", WLorentzVector, {"w_pair_pt", "w_pair_eta", "w_pair_phi", "w_mass", "w_reco_jets"})
				  .Define("bjetmass", bjet_variable, {JetMassInput, "nJet", "lead_bjet"})
				  .Define("bjetpt", bjet_variable, {JetPtInput, "nJet", "lead_bjet"})
			          .Define("bjeteta", bjet_variable, {JetEtaInput, "nJet", "lead_bjet"})
				  .Define("bjetphi", bjet_variable, {JetPhiInput, "nJet", "lead_bjet"})
				  .Define("BJets", BLorentzVector, {"bjetpt", "bjeteta", "bjetphi", "bjetmass"})
				  .Define("RecoTop", top_reconstruction_function, {"bjetpt", "bjeteta", "bjetphi", "bjetmass", 	
										   "w_pair_pt", "w_pair_eta", "w_pair_phi", "w_mass"})
			          .Define("Top_Pt", TLorentzVectorVariablePt, {"RecoTop"})
			          .Define("Top_Eta", TLorentzVectorVariableEta, {"RecoTop"})
			          .Define("Top_Phi", TLorentzVectorVariablePhi, {"RecoTop"})
			          .Define("Top_Mass", TLorentzVectorVariableMass, {"RecoTop"})
				  .Define("Top_HT", HT_double, {"Top_Pt"})
				  .Define("dR_Top_LeadingElectron", deltaRcheck_Top_function, {"Top_Phi", "Top_Eta", "LeadingLeptonEta", "LeadingLeptonPhi"})
			          .Define("dR_Top_SubleadingElectron", deltaRcheck_Top_function, {"Top_Phi", "Top_Eta", "SubleadingLeptonEta", "SubleadingLeptonPhi"})
				  .Define("dR_Top_LeadingJet", deltaRcheck_Top_function, {"Top_Phi", "Top_Eta", "LeadingJetEta", "LeadingJetPhi"})
			          .Define("dR_Top_SubleadingJet", deltaRcheck_Top_function, {"Top_Phi", "Top_Eta", "SubleadingJetEta", "SubleadingJetPhi"})
				  .Define("dR_Top_ThirdJet", deltaRcheck_Top_function, {"Top_Phi", "Top_Eta", "ThirdJetEta", "ThirdJetPhi"})
				  .Define("dR_Top_FourthJet", deltaRcheck_Top_function, {"Top_Phi", "Top_Eta", "FourthJetEta", "FourthJetPhi"})
			          .Define("dR_Top_W", deltaRcheck_WTop_function, {"w_pair_eta", "w_pair_phi", "Top_Eta", "Top_Phi"})
				  .Define("dPhi_Wj1_Top", DeltaPhi_function2, {"WPairJet1Phi", "Top_Phi"})
				  .Define("dPhi_Wj2_Top", DeltaPhi_function2, {"WPairJet2Phi", "Top_Phi"})
			          .Define("dR_Z_Top", deltaRcheck_W_function, {"RecoZPhi", "RecoZEta", "Top_Phi", "Top_Eta"})
				  .Define("dPhi_Z_Top", DeltaPhi_function2, {"Top_Phi", "RecoZPhi"})
				  .Define("dR_Z_WPairJet1", deltaRcheck_W_function, {"RecoZPhi", "RecoZEta", "WPairJet1Eta", "WPairJet1Phi"})
				  .Define("dR_Z_WPairJet2", deltaRcheck_W_function, {"RecoZPhi", "RecoZEta", "WPairJet2Eta", "WPairJet2Phi"})
			          .Define("dPhi_Z_WPairJet1", DeltaPhi_function2, {"RecoZPhi", "WPairJet1Phi"})
                                  .Define("dPhi_Z_WPairJet2", DeltaPhi_function2, {"RecoZPhi", "WPairJet2Phi"})
				  .Define("MinDeltaR", MinDeltaR, {"nJet", "RecoZPhi", "RecoZEta", JetPhiInput, JetEtaInput})
				  .Define("MinDeltaPhi", MinDeltaPhi, {"nJet", "RecoZPhi", JetPhiInput})
				  .Define("dR_LeadingLepton_LeadingBJet", dR_Lepton_LeadingBJet_Function, {"bjeteta", "LeadingLeptonEta", "bjetphi", "LeadingLeptonPhi"})
			          .Define("dR_SubleadingLepton_LeadingBJet", dR_Lepton_LeadingBJet_Function, {"bjeteta", "SubleadingLeptonEta", "bjetphi", "SubleadingLeptonPhi"})
				  .Define("DeltaPhi_Leadinglepton_BJet", DeltaPhi_Lepton_BJet, {JetPhiInput, "LeadingLeptonPhi"})
                                  .Define("DeltaPhi_Subleadinglepton_BJet", DeltaPhi_Lepton_BJet, {JetPhiInput, "SubleadingLeptonPhi"})		
				  .Define("MET", MET_function, {"MET_sumEt"})
			          .Define("LeadingBJetOutputDiscriminant", BJetOutputDiscriminantFunction, {"LeadingJetPt", "Jet_btagCSVV2", "tight_jets", JetEtaInput})
                                  .Define("SubleadingBJetOutputDiscriminant", BJetOutputDiscriminantFunction, {"SubleadingJetPt", "Jet_btagCSVV2", "tight_jets", JetEtaInput})
                                  .Define("ThirdBJetOutputDiscriminant", BJetOutputDiscriminantFunction, {"ThirdJetPt", "Jet_btagCSVV2", "tight_jets", JetEtaInput})
                                  .Define("FourthBJetOutputDiscriminant", BJetOutputDiscriminantFunction, {"FourthJetPt", "Jet_btagCSVV2", "tight_jets", JetEtaInput})
                                  .Define("dPhi_W_Top", DeltaPhi_function4, {"w_pair_phi", "Top_Phi"})
				  .Define("dR_Z_LeadingJet", deltaRcheck_W_function2, {"RecoZPhi", "RecoZEta", "LeadingJetPhi", "LeadingJetEta"})
                                  .Define("dR_Z_SubleadingJet", deltaRcheck_W_function2, {"RecoZPhi", "RecoZEta", "SubleadingJetPhi", "SubleadingJetEta"})
                                  .Define("dR_Z_ThirdJet", deltaRcheck_W_function2, {"RecoZPhi", "RecoZEta", "ThirdJetPhi", "ThirdJetEta"})
                                  .Define("dR_Z_FourthJet", deltaRcheck_W_function2, {"RecoZPhi", "RecoZEta", "FourthJetPhi", "FourthJetEta"})
                                  .Define("dPhi_LeadingJet_Z", DeltaPhi_doublesandfloat, {"RecoZPhi", "LeadingJetPhi"})
                                  .Define("dPhi_SubleadingJet_Z", DeltaPhi_doublesandfloat, {"RecoZPhi", "SubleadingJetPhi"})
                                  .Define("dPhi_ThirdJet_Z", DeltaPhi_doublesandfloat, {"RecoZPhi", "ThirdJetPhi"})
                                  .Define("dPhi_FourthJet_Z", DeltaPhi_doublesandfloat, {"RecoZPhi", "FourthJetPhi"})
                                  .Define("dR_W_Z", deltaRcheck_WTop_function, {"w_pair_phi", "w_pair_eta", "RecoZEta", "RecoZPhi"})
                                  .Define("RecoZHT", HT_double, {"RecoZPt"})
                                  .Define("dPhi_W_Z", DeltaPhi_function4, {"w_pair_phi", "RecoZPhi"})
                                  .Define("TotalEta_System", TotalVariable_System, {"RecoZEta", "w_pair_eta", "Top_Eta", "LepEtaSum", "JetEtaSum"})
                                  .Define("TotalPhi_System", TotalVariable_System, {"RecoZPhi", "w_pair_phi", "Top_Phi", "LepPhiSum", "JetPhiSum"})
				  .Define("InvTopMass", inv_mass_doubles, {"Top_Pt", "Top_Eta", "Top_Phi", "Top_Mass"})
				  .Define("UnweightedTopPt", UnweightedTopPt, {"Top_Pt"})
				  .Define("TopReweighting_topquark", TopReweighting_topquark, {"GenPart_pdgId", "GenPart_statusFlags", "GenPart_pt"})
                                  .Define("TopReweighting_antitopquark", TopReweighting_antitopquark, {"GenPart_pdgId", "GenPart_statusFlags", "GenPart_pt"})
                                  .Define("TopWeight", TopReweighting_weight, {"TopReweighting_topquark", "TopReweighting_antitopquark"});


  std::string TopCandRecoFile = "TopCandReco_" + Process + "_" + Systematic + "_" + Channel + "_" + NonPromptLepton + "_" +
                                SignalRegion + "_" + SideBandRegion + "_" + ZPlusJetsControlRegion + "_" + ttbarControlRegion + "_" + Year + ".root";

  auto Snapshot_TopCandReco = d_TopCandReco.Snapshot("Events", TopCandRecoFile.c_str());


  auto d_EventWeightDefines = d_TopCandReco.Define("TotalHT_System", TotalVariable_System, {"RecoZHT", "RecoWHT", "Top_HT", "TotLepHT", "TotJetHT"})
                                           .Define("TotalPt_System", TotalVariable_System, {"RecoZPt", "w_pair_pt", "Top_Pt", "LepPtSum", "JetPtSum"})
					   .Define("TotHTOverTotpT_System", TotHTOverTotpT_doubles, {"TotalHT_System", "TotalPt_System"})
					   .Define("CMSBTagSF", CMSBTagSF, {"bjetpt", "bjeteta", "Jet_btagCSVV2", "Jet_partonFlavour"})
					   .Define("nonbjets", nonbjet_id, {"tight_jets", "Jet_btagCSVV2", JetEtaInput})
                                           .Define("notbjetpt", bjet_variable, {JetPtInput, "nJet", "nonbjets"})
                                           .Define("notbjeteta", bjet_variable, {JetEtaInput, "nJet", "nonbjets"})
  					   .Define("CMSNonBTagSF", CMSNonBTagSF, {"notbjetpt", "notbjeteta", "Jet_btagCSVV2", "Jet_partonFlavour"})
					   .Define("EffBTagged", EffBTagged_Function, {JetPtInput, JetEtaInput})
					   .Define("EffNonBTagged", EffNonBTagged_Function, {JetPtInput, JetEtaInput})
					   .Define("EffBTaggedProduct", EffBTaggedProduct, {"EffBTagged"})
					   .Define("EffNonBTaggedProduct", EffNonBTaggedProduct, {"EffNonBTagged"})
					   .Define("EffBTaggedProductData", EffBTaggedProductData, {"EffBTagged", "CMSBTagSF"})
                                           .Define("EffNonBTaggedProductData", EffNonBTaggedProductData, {"EffNonBTagged", "CMSNonBTagSF"})
					   .Define("ProbBTagMC", ProbBTagMCFunction, {"EffBTaggedProduct", "EffNonBTaggedProduct"})
 					   .Define("ProbBTagData", ProbBTagDataFunction, {"EffBTaggedProductData", "EffNonBTaggedProductData"})
					   .Define("BTagWeight", BTagWeightFunction, {"ProbBTagMC", "ProbBTagData"})
					   .Define("EGammaSF_egammaEff", EGammaSF_egammaEff, {"TightLeptonsPt", "TightLeptonsEta"})
					   .Define("EGammaSF_egammaEffSys", EGammaSF_egammaEff_Sys, {"TightLeptonsPt", "LeptonEta"})
					   .Define("EGammaSF_egammaEffReco", EGammaSF_egammaEffReco, {"TightLeptonsPt", "TightLeptonsEta"})
					   .Define("EGammaSF_egammaEffRecoSys", EGammaSF_egammaEffReco_Sys, {"TightLeptonsPt", "TightLeptonsEta"})
					   .Define("MuonSFTest_ID", MuonSFTest_ID, {"LeptonPt_RochCorr", "LeptonEta_RochCorr"})
					   .Define("MuonSFTest_Iso", MuonSFTest_Iso, {"LeptonPt_RochCorr", "LeptonEta_RochCorr"})
                                           .Define("MuonSFTest_ID_sys_syst", MuonSFTest_ID_sys_syst, {"LeptonPt_RochCorr", "LeptonEta_RochCorr"})
				           .Define("MuonSFTest_ID_sys_stat", MuonSFTest_ID_sys_stat, {"LeptonPt_RochCorr", "LeptonEta_RochCorr"})
					   .Define("MuonSFTest_Iso_sys_syst", MuonSFTest_Iso_sys_syst, {"LeptonPt_RochCorr", "LeptonEta_RochCorr"})
                                           .Define("MuonSFTest_Iso_sys_stat", MuonSFTest_Iso_sys_stat, {"LeptonPt_RochCorr", "LeptonEta_RochCorr"})
					   .Define("ReturnedPSWeight", PSWeightFunction, {PSWeightString})
					   .Define("CalculatedNominalWeight", NominalWeight, {"LHEPdfWeight", "LHEWeight_originalXWGTUP"})
					   .Define("ME_SF", ME_uncert_function, {"LHEPdfWeight", "LHEWeight_originalXWGTUP", "ReturnedPSWeight"})
					   .Define("CalculatedGeneratorWeight", GeneratorWeight, {"CalculatedNominalWeight", "ReturnedPSWeight"})
					   .Define("OriginalMET", OriginalMetFunction, {"MET_sumEt", "MET_phi"})
					   .Define("ScaledMET", ScaledMetFunction, {"OriginalMET", "MET_sumEt", "MET_phi", "MET_MetUnclustEnUpDeltaX", "MET_MetUnclustEnUpDeltaY"})
					   .Define("UnsmearedJet4Momentum", UnsmearedJetTLorentzVectorFunction, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass"})
					   .Define("newMET", METUncertFunction, {"ScaledMET", "SmearedJet4Momentum", "UnsmearedJet4Momentum"})
					   .Define("EventWeight", EventWeight, {"PU", "BTagWeight", "ReturnedPSWeight", "EGammaSF_egammaEff", 
										"EGammaSF_egammaEffReco", "EGammaSF_egammaEffSys", "EGammaSF_egammaEffRecoSys", 
										"CalculatedGeneratorWeight", "ME_SF", "TopWeight", "CalculatedNominalWeight", "MuonSFTest_ID", "MuonSFTest_Iso", 
										"MuonSFTest_ID_sys_syst", "MuonSFTest_ID_sys_stat", "MuonSFTest_Iso_sys_syst", 
										"MuonSFTest_Iso_sys_stat"});
								      


  std::string EventWeightFile = "EventWeight_" + Process + "_" + Systematic + "_" + Channel + "_" + NonPromptLepton + "_" +
                                SignalRegion + "_" + SideBandRegion + "_" + ZPlusJetsControlRegion + "_" + ttbarControlRegion + "_" + Year + ".root";

  auto Snapshot_EventWeight = d_EventWeightDefines.Snapshot("Events", EventWeightFile.c_str());

  if(ProcessInt == 0 && SystematicInt == 0){

  	auto h_WMass = d_EventWeightDefines.Histo1D("w_mass", "EventWeight");
  	auto h_InvTopMass = d_EventWeightDefines.Histo1D("InvTopMass", "EventWeight"); 
	auto h_WMass_Unweighted = d_EventWeightDefines.Histo1D("w_mass");
        auto h_InvTopMass_Unweighted = d_EventWeightDefines.Histo1D("InvTopMass");
	auto h_EventWeight = d_EventWeightDefines.Histo1D("EventWeight");	

	h_WMass->Fit("gaus"); 
        h_InvTopMass->Fit("gaus");
  	h_WMass_Unweighted->Fit("gaus"); 
  	h_InvTopMass_Unweighted->Fit("gaus");

  	W_stddev = h_WMass->GetStdDev(); 
  	Top_stddev = h_InvTopMass->GetStdDev();

	std::cout << "W_stddev = " << std::endl;
	std::cout << "Top_stddev = " << std::endl;

  	std::string GaussianFitsFileString = "GaussianFits_" + Process + "_" + Systematic + "_" + Channel + "_" + NonPromptLepton + "_" +
                                             SignalRegion + "_" + SideBandRegion + "_" + ZPlusJetsControlRegion + "_" + ttbarControlRegion + "_" + Year + ".root";

        TFile * GaussianFitsFile = new TFile(GaussianFitsFileString.c_str(), "RECREATE");

	h_WMass->Write();
	h_InvTopMass->Write();
	h_WMass_Unweighted->Write();
        h_InvTopMass_Unweighted->Write();
	h_EventWeight->Write();

	GaussianFitsFile->Close();

  

  	//Write the nominal mass and resolution values to a text file 
  	std::ofstream Resolution;
  	std::string ResolutionString = "Resolution_" + Process + "_" + Systematic + "_" + Channel + "_" + NonPromptLepton + "_" +
                                       SignalRegion + "_" + SideBandRegion + "_" + ZPlusJetsControlRegion + "_" + ttbarControlRegion + "_" + Year + ".txt";

  	Resolution.open(ResolutionString.c_str());

  	Resolution << "W_stddev: " << W_stddev << '\n'
	           << "Top_stddev: " << Top_stddev << std::endl;

  }

  auto d_Blinding =  d_EventWeightDefines.Define("chi2", Chi2Function, {"w_mass", "InvTopMass"});

  std::cout << "ProcessInt = " << ProcessInt << std::endl;
  std::cout << "SystematicInt = " << SystematicInt << std::endl;


  if(ProcessInt == 0 && SystematicInt == 0){

  	int NumberOfSimulatedEvents = *( d_Blinding.Filter("chi2 != 999.0").Count() );	
	int OneSigmaOfNumEvents = NumberOfSimulatedEvents * 0.68;
			
	auto histo_chi2 = d_Blinding.Histo1D("chi2");
	TAxis *xaxis = histo_chi2->GetXaxis();
	double MaxBin = xaxis->GetBinCenter( histo_chi2->FindLastBinAbove() );
	double MinBin = xaxis->GetBinCenter( histo_chi2->FindFirstBinAbove() ) ;

	int NumBins = MaxBin - MinBin;

	auto histo_chi2_rebinned = d_Blinding.Histo1D({"histo_chi2_rebinned", "histo_chi2_rebinned", 2*NumBins, MinBin, MaxBin}, {"chi2"});

	TAxis * histo_chi2_rebinned_x = histo_chi2_rebinned->GetXaxis();

	int total = 0;

	for(int i = 0; i < histo_chi2_rebinned->GetEntries(); i++){
		
		auto NumberOfEvents = histo_chi2_rebinned->GetBinContent(i);
		total += NumberOfEvents;

		if(total >= OneSigmaOfNumEvents){Chi2_SR = histo_chi2_rebinned_x->GetBinCenter(i); break;}
		else{continue;}

	}

	std::cout << "before MaxElement" << std::endl;

	//for Chi2_SBR
	auto MaxElement = std::max_element(CutRanges.begin(), CutRanges.end());
	auto DistToMax = std::distance(CutRanges.begin(), MaxElement);
	Chi2_SBR = CutRanges.at(DistToMax);

	std::ofstream Chi2Range;      
	std::string Chi2Range_string = "Chi2Range_" + Process + "_" + Systematic + "_" + Channel + "_" + NonPromptLepton + "_" +
                                       SignalRegion + "_" + SideBandRegion + "_" + ZPlusJetsControlRegion + "_" + ttbarControlRegion + "_" + Year + ".txt";
 
        Chi2Range.open(Chi2Range_string.c_str());

        Chi2Range << "Chi2_SR: " << Chi2_SR << '\n'
                  << "Chi2_SBR: " << Chi2_SBR << std::endl;

  }

  auto d_Blinded = d_Blinding.Define("AfterChi2Cut", Chi2Cut, {"chi2"}).Filter(Chi2Cut, {"chi2"});
  auto colNames = d_Blinding.GetDefinedColumnNames();
  const auto N_Columns = colNames.size();

  std::string OutRootFileStart;

  switch(MCInt){case 1: OutRootFileStart = "Results_MC_"; break;
		default: OutRootFileStart = "Results_Data_"; break;}

  std::string OutRootFile = OutRootFileStart + Process + "_" + Systematic + "_" + Channel + "_" + NonPromptLepton + "_" +
                            SignalRegion + "_" + SideBandRegion + "_" + ZPlusJetsControlRegion + "_" + ttbarControlRegion + "_" + Year + ".txt";

  TFile * output = new TFile(OutRootFile.c_str(), "RECREATE");
  output->cd();

  ROOT::RDF::RResultPtr<TH1D> histo[N_Columns] = {};

  for(long unsigned int i = 0; i < N_Columns; i++){

  	auto ColName = colNames.at(i);

	if(ColName != "PU"                      && ColName != "BTagWeight"                && ColName != "ReturnedPSWeight"              &&
           ColName != "CalculatedNominalWeight" && ColName != "EGammaSF_egammaEff"        && ColName != "EGammaSF_egammaEffReco"        &&
           ColName != "EGammaSF_egammaEffSys"   && ColName != "EGammaSF_egammaEffRecoSys" && ColName != "CalculatedGeneratorWeight"     &&
           ColName != "ME_SF"                   &&
           ColName != "RecoZ"                   && ColName != "SmearedJet4Momentum"       && ColName != "WPairJet1"                     && 
           ColName != "WPairJet2"               && ColName != "RecoW"                     && ColName != "BJets"                         && 
           ColName != "RecoTop"                 && ColName != "MinDeltaR"                     &&
           ColName != "MinDeltaPhi"             && ColName != "newMET"                    && ColName != "EventWeight" ){

           	std::cout << "ColName = " << ColName << std::endl;

                histo[i] = d_Blinded.Histo1D(ColName.c_str(), "EventWeight");
                histo[i]->Write();
                        
        }
 	else if(ColName  == "PU"                      || ColName == "BTagWeight"                || ColName == "ReturnedPSWeight"          ||
                ColName  == "CalculatedNominalWeight" || ColName == "EGammaSF_egammaEff"        || ColName == "EGammaSF_egammaEffReco"    ||
                ColName  == "EGammaSF_egammaEffSys"   || ColName == "EGammaSF_egammaEffRecoSys" || ColName == "CalculatedGeneratorWeight" ||
                ColName  == "ME_SF"                   || ColName == "EventWeight" ){

                histo[i] = d_Blinded.Histo1D(ColName.c_str());
                histo[i]->Write();

        }		
        else{continue;}

  }

  output->Close();	


  //Cut flow report
  auto allCutsReport = d.Report();
  std::ofstream CutFlowReport;
  std::string CutFlowReportString = "CutFlowReport_" + Process + "_" + Systematic + "_" + Channel + "_" + NonPromptLepton + "_" +
                                    SignalRegion + "_" + SideBandRegion + "_" + ZPlusJetsControlRegion + "_" + ttbarControlRegion + "_" + Year + ".txt"; 

  CutFlowReport.open(CutFlowReportString.c_str());

  for(auto&& cutInfo: allCutsReport){
  	CutFlowReport << cutInfo.GetName() << '\t' << cutInfo.GetAll() << '\t' << cutInfo.GetPass() << '\t' << cutInfo.GetEff() << " %" << std::endl;
  }


}


//Main function is here
void fulleventselectionAlgo::fulleventselection(){

  int MC_Selection = 1;
  std::vector<int> Process_Selection = {95, 96, 0}; //95 for trigger SF MC, 96 for trigger SF data, 0 for tZq
  int NPL_Selection = 0;
  int SR_Selection = 1;
  int SBR_Selection = 1;
  int ZPlusJetsCR_Selection = 0;
  int ttbarCR_Selection = 0;
  int Year_Selection = 2017;
  int Systematic_Selection = 0;
  int Channel_Selection = 1;
  int DoubleCountCheck_Selection = 0; //set this to 1 when running over double electron, double muon, single electron, single muon or MuonEG samples


  for(int i = 0; Process_Selection.size(); i++){

  	tZq_NanoAOD_Output(MC_Selection, 	    Process_Selection.at(i), NPL_Selection,        SR_Selection,         SBR_Selection,             ZPlusJetsCR_Selection, 
		           ttbarCR_Selection,       Year_Selection,          Systematic_Selection, Channel_Selection,    DoubleCountCheck_Selection);

  }


}



