#include "fulleventselectionAlgo.hpp"
#include "RoccoR.hpp"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TH1D.h"
#include <TFile.h>
#include "TEfficiency.h"
#include "TLorentzVector.h"
//#include "LorentzVector.h"
#include "Math/Vector4D.h"
#include <TRandom3.h>
#include "TCanvas.h"

#include <vector>
#include <iterator>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <random>


#include <boost/numeric/conversion/cast.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>


fulleventselectionAlgo::fulleventselectionAlgo()
    : MC_Selection_{},
    Year_Selection_{}, // No default value set as it currently is required to be set by the user
    Process_Selection_{},
    NPL_Selection_{},
    SR_Selection_{},
    SBR_Selection_{},
    ZPlusJetsCR_Selection_{},
    ttbarCR_Selection_{},
    Systematic_Selection_{},
    Channel_Selection_{},
    DoubleCountCheck_Selection_{} 
//    , anotherVariable_{defaultValue}

{}

fulleventselectionAlgo::~fulleventselectionAlgo() {}

void fulleventselectionAlgo::parseCommandLineArguements(int argc, char* argv[]){

    gErrorIgnoreLevel = kInfo;
    // Set up environment a little.
    std::cerr << std::setprecision(6) << std::fixed;
    std::cout << std::setprecision(6) << std::fixed;

    namespace po = boost::program_options;
    po::options_description desc("Options");
  
   desc.add_options()("help,h", "Print this message.")(
        "mc, m",
        po::value<int>(&MC_Selection_)->required(),
        "MC selection (e.g. 0 is data, 1 is MC)")(
        "year,y",
        po::value<int>(&Year_Selection_)->required(), // No default value, has to be explicitly set by the user
        "Year selection. Either 2016, 2017, or 2018.")(
        "process,p",
	po::value<int>(&Process_Selection_)->required(),
	"Process selection (e.g. 0 is tZq)")(
        "npl, n",
        po::value<int>(&NPL_Selection_)->required(),
        "NPL selection (e.g. 0 is for not the NPL run, 1 is for the NPL run)")(
        "sr, S",
        po::value<int>(&SR_Selection_)->required(), // No default value, has to be explicitly set by the user
        "Signal region. Set SR_Selection to 1 for results in the signal region.")(
        "sbr, s",
        po::value<int>(&SBR_Selection_)->required(),
        "Side band region selection. Set both SR and SBR equal to 1 for results in the side band region.")(
        "zjcr, z",
        po::value<int>(&ZPlusJetsCR_Selection_)->required(),
        "Z+jets control region selection. Set to 1 for results in the z+jets control region, otherwise set to 0.")(
        "ttcr, t",
        po::value<int>(&ttbarCR_Selection_)->required(),
        "ttbar control region selection. Set to 1 for results in the ttbar control region, otherwise set to 0.")(
        "sys, u",
        po::value<int>(&Systematic_Selection_)->required(),
        "Systematic selection. 0 = nominal run")(
        "channel, c",
        po::value<int>(&Channel_Selection_)->required(),
        "Channel selection. ee channel = 1, mumu channel = 2, emu channel = 3")(
        "dcc, d",
        po::value<int>(&DoubleCountCheck_Selection_)->required(),
        "Set to 1 when running over the single or double lepton datasets.");

    po::variables_map vm;

    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);

        if (vm.count("help")) {
            std::cout << desc;
            std::exit(0);
        }

        po::notify(vm);
    }
    catch (const std::logic_error& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cerr << "Use -h or --help for help." << std::endl;
        std::exit(1);
    }
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
constexpr float TOP_MASS = 172.6;

float Chi2_SR; 
float Chi2_SBR;


bool batch = true;

//method to convert the text file to vector of ROOT files 
std::vector<std::string>
TxtFileListToVecStr(const std::string f)
{
    std::ifstream infile(f, std::ifstream::in);

    int i = 0;

    std::vector<std::string> v;
    std::string line("");

    while (std::getline(infile, line))
    {

	i++;

        line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        line.erase(std::remove(line.begin(), line.end(), '\t'), line.end());
        if (line.find_first_not_of(' ') == std::string::npos)
            continue;
        if (line[line.find_first_not_of(' ')] == '#')
            continue;

        if (batch)
        {
            v.push_back("root://cms-xrd-global.cern.ch//" + line);
        }
        else
        {
            v.push_back("root://" + line);
        }
        
    }

    infile.close();
    return v;
}



template<typename T>
[[gnu::const]] T select(const T& a, const ints& mask)
{
  std::cout << "print 200" << std::endl;
  return a[mask];
}

template<typename T>
[[gnu::const]] T select_floats(const T& a, const floats& mask)
{
  std::cout << "print 201" << std::endl;
  return a[mask];
}

template<typename T, typename U> //for the all equal function
[[gnu::const]] bool all_equal(const T& t, const U& u)
{
  std::cout << "print 202" << std::endl;
  return t == u;
}

template<typename T, typename U, typename... Types>
[[gnu::const]] bool all_equal(const T& t, const U& u, Types const&... args)
{
    std::cout << "print 203" << std::endl;
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



double N_SelectionCriteria_MC;
double N_MET_And_LeptonSelection_MC; 
double N_LeptonTriggersAndSelectionCriteria_MC;
double N_MET_LeptonTriggers_SelectionCriteria_MC;
double Eff_MC;
double Eff_MET_LeptonTriggers_SelectionCriteria_MC;
double Eff_LeptonTriggers_SelectionCriteria_MC;
double Eff_MET_SelectionCriteria_MC;
double Alpha_MC;
double Eff_UpperUncert_MC;
double Eff_LowerUncert_MC;

double N_SelectionCriteria_DATA; 
double N_MET_And_LeptonSelection_DATA; 
double N_LeptonTriggersAndSelectionCriteria_DATA; 
double N_MET_LeptonTriggers_SelectionCriteria_DATA; 
double Eff_DATA; 
double Eff_MET_LeptonTriggers_SelectionCriteria_DATA;
double Eff_LeptonTriggers_SelectionCriteria_DATA; 
double Eff_MET_SelectionCriteria_DATA; 
double Alpha_DATA; 
double Eff_UpperUncert_DATA; 
double Eff_LowerUncert_DATA; 

double TrigSF;
double TrigSF_Uncert;


void tZq_NanoAOD_Output(const int& MCInt,  	    const int& ProcessInt,  const int& NPLInt,     const int& SRInt,          const int& SBRInt, 
		 	const int& ZPlusJetsCRInt,  const int& ttbarCRInt,  const int& YearInt,    const int& SystematicInt,  const int& ChannelInt, 
			const int& DoubleCountCheckInt){

  std::cout << "MCInt = " << MCInt << std::endl;
  std::cout << "YearInt = " << YearInt << std::endl; 
  std::cout << "ProcessInt = " << ProcessInt << std::endl;
  std::cout << "NPLInt = " << NPLInt << std::endl;
  std::cout << "SRInt = " << SRInt << std::endl;
  std::cout << "SBRInt = " << SBRInt << std::endl;
  std::cout << "ZPlusJetsCRInt = " << ZPlusJetsCRInt << std::endl;
  std::cout << "ttbarCRInt = " << ttbarCRInt << std::endl;
  std::cout << "SystematicInt = " << SystematicInt << std::endl;
  std::cout << "ChannelInt = " << ChannelInt << std::endl;
  std::cout << "DoubleCountCheckInt = " << DoubleCountCheckInt << std::endl;


  ROOT::RDF::RResultPtr<TH2D> h_bjet_num;
  ROOT::RDF::RResultPtr<TH2D> h_bjet_denom;
  ROOT::RDF::RResultPtr<TH2D> h_nonbjet_num;
  ROOT::RDF::RResultPtr<TH2D> h_nonbjet_denom;

  TH2D * h_bjet;
  TH2D * h_nonbjet;
 
  std::string SJER;

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
  std::string GeneratorWeightString; 
  std::vector<std::string> LeptonGenPartFlavStrings;
  std::vector<std::string> RochCorrVecStrings;
  std::vector<std::string> JetSmearingStrings;
  std::vector<std::string> TightSmearedJetsHadronFlavourStrings;
  std::vector<std::string> TopReweighting_topquark_strings;
  std::vector<std::string> TopReweighting_antitopquark_strings;
  std::vector<std::string> PDFWeightStrings;

  std::string JetMassInput;
  std::string JetPtInput;
  std::string JetEtaInput; 
  std::string JetPhiInput;
  
  std::string PSWeightString;
  std::string DoubleCountString;
	
  std::string HessianOrMC;
  std::string Tune;

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


  TH2* histo_RunsABCD_ID_2018_stat = dynamic_cast<TH2*>(inputfile_RunsABCD_ID_2018->Get("NUM_TightID_DEN_TrackerMuons_pt_abseta_stat")->Clone());
  histo_RunsABCD_ID_2018_stat->SetDirectory(nullptr);


  TH2* histo_RunsABCD_ID_2018_syst = dynamic_cast<TH2*>(inputfile_RunsABCD_ID_2018->Get("NUM_TightID_DEN_TrackerMuons_pt_abseta_syst")->Clone());
  histo_RunsABCD_ID_2018_syst->SetDirectory(nullptr);


  TH2* histo_RunsABCD_ISO_2018_stat = dynamic_cast<TH2*>(inputfile_RunsABCD_ISO_2018->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta_stat")->Clone());
  histo_RunsABCD_ISO_2018_stat->SetDirectory(nullptr);


  TH2* histo_RunsABCD_ISO_2018_syst = dynamic_cast<TH2*>(inputfile_RunsABCD_ISO_2018->Get("NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta_syst")->Clone());
  histo_RunsABCD_ISO_2018_syst->SetDirectory(nullptr);


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



  std::cout << "MCInt = " << MCInt << std::endl;

  //Setting the SampleType, Channel, Process, NPL, SR, SBR, ZPlusJetsCR, ttbarCR, Year and Systematic strings for the output file names
  switch(SystematicInt){

        case 0: Systematic = "Nominal"; SJER = "sJER_Nominal"; break;
        case 1: Systematic = "PU_ScaleUp"; SJER = "sJER_Nominal"; break;
        case 2: Systematic = "PU_ScaleDown"; SJER = "sJER_Nominal"; break;
        case 3: Systematic = "BTag_ScaleUp"; SJER = "sJER_Nominal"; break;
        case 4: Systematic = "BTag_ScaleDown"; SJER = "sJER_Nominal"; break;
        case 5: Systematic = "JetSmearing_ScaleUp"; SJER = "sJER_up"; break;
        case 6: Systematic = "JetSmearing_ScaleDown"; SJER = "sJER_down"; break;
        case 7: Systematic = "JEC_ScaleUp"; SJER = "sJER_Nominal"; break;
        case 8: Systematic = "JEC_ScaleDown"; SJER = "sJER_Nominal"; break;
        case 9: Systematic = "LeptonEfficiencies_ScaleUp"; SJER = "sJER_Nominal"; break;
        case 10: Systematic = "LeptonEfficiencies_ScaleDown"; SJER = "sJER_Nominal"; break;
        case 11: Systematic = "PDF_ScaleUp"; SJER = "sJER_Nominal"; break;
        case 12: Systematic = "PDF_ScaleDown"; SJER = "sJER_Nominal"; break;
        case 13: Systematic = "ME_Up"; SJER = "sJER_Nominal"; break;
        case 14: Systematic = "ME_Down"; SJER = "sJER_Nominal"; break;
        case 15: Systematic = "MET_Up"; SJER = "sJER_Nominal"; break;
        case 16: Systematic = "MET_Down"; SJER = "sJER_Nominal"; break;
        case 17: Systematic = "isr_up"; SJER = "sJER_Nominal"; break;
        case 18: Systematic = "isr_down"; SJER = "sJER_Nominal"; break;
        case 19: Systematic = "fsr_up"; SJER = "sJER_Nominal"; break;
        case 20: Systematic = "fsr_down"; SJER = "sJER_Nominal"; break;
        default: std::cout << "ERROR: SystematicInt must be between 1 and 20." << std::endl; break;

  } 
  

  switch(MCInt){

	case 0: SampleType = "data";
		JetMassInput = "Jet_mass"; JetPtInput = "Jet_pt"; JetEtaInput = "Jet_eta"; JetPhiInput = "Jet_phi"; 
		GeneratorWeightString = "CaloMET_pt"; //input is not used for the function for data, so have put anything here
		LeptonGenPartFlavStrings = {"Electron_cleanmask", "Electron_cleanmask"}; //just a random char RDataFrame to use for data, since it won't be used 
		RochCorrVecStrings = {"TightLeptonsCharge", "TightLeptonsPt", "TightLeptonsEta", "TightLeptonsPhi", "Muon_charge", "Muon_nTrackerLayers"};
		JetSmearingStrings = {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_pt", "Jet_eta", "Jet_phi", SJER, "sigma_JER", "nJet"}; //repeated variables aren't used anyway since a correction factor of 1 is applied for data
		TightSmearedJetsHadronFlavourStrings = {"Jet_jetId", "tight_jets"};
		TopReweighting_topquark_strings = {"nJet", "nJet", "Jet_pt"}; //using any input string here as a weight of 1 is returned anyway for data
                TopReweighting_antitopquark_strings = {"nJet", "nJet", "Jet_pt"}; //using any input string here as a weight of 1 is returned anyway for data
		PDFWeightStrings = {"PV_score", "nMuon"}; //any float and UInt here for data, as a weight of 1 is applied anyway
		break;

	case 1: SampleType = "MC";
		JetMassInput = "SmearedJetMass"; JetPtInput = "SmearedJetPt"; JetEtaInput = "SmearedJetEta"; JetPhiInput = "SmearedJetPhi";
		GeneratorWeightString = "genWeight";
		LeptonGenPartFlavStrings = {"Electron_genPartFlav", "Muon_genPartFlav"};
		RochCorrVecStrings = {"TightLeptonsCharge", "TightLeptonsPt", "TightLeptonsEta", "TightLeptonsPhi", "Muon_genPartIdx", "Muon_nTrackerLayers"};
		JetSmearingStrings = {"Jet_pt", "Jet_eta", "Jet_phi", "GenJet_pt", "GenJet_eta", "GenJet_phi", SJER, "sigma_JER", "Jet_genJetIdx"};
		TightSmearedJetsHadronFlavourStrings = {"Jet_hadronFlavour", "tight_jets"};
		TopReweighting_topquark_strings = {"GenPart_pdgId", "GenPart_statusFlags", "GenPart_pt"};
		TopReweighting_antitopquark_strings = {"GenPart_pdgId", "GenPart_statusFlags", "GenPart_pt"};
		PDFWeightStrings = {"LHEPdfWeight", "nLHEPdfWeight"};
		break;

  }

  switch(ProcessInt){

	case 0: Process = "tZq"; 

		switch(YearInt){
			case 2016: HessianOrMC = "Hessian"; break;
			case 2017: HessianOrMC = "Hessian"; break;
			case 2018: HessianOrMC = "Hessian"; break;
			default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
		}

		break;

	case 1: Process = "tZq_scaleup"; 

		switch(YearInt){
			case 2016: HessianOrMC = "MC"; break;
			case 2017: break;
			case 2018: break;
			default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
		}

		break;

	case 2: Process = "tZq_scaledown"; 

                switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                	default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
		}
		
		break;


	case 3: Process = "ZPlusJets_M50_aMCatNLO";
		
                switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;

                        case 2017: HessianOrMC = "Hessian"; break;

                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                }
                
                break;

	case 4: Process = "ZPlusJets_M50_aMCatNLO_ext";
	
		switch(YearInt){
                        case 2016: break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                }

                break;

	case 5: Process = "ZPlusJets_M50_Madgraph";

		switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                }

                break;

	case 6: Process = "ZPlusJets_M50_Madgraph_ext";

		switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                }

                break;

	case 7: Process = "ZPlusJets_M10To50_aMCatNLO";

		switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                }

                break;

        case 8: Process = "ZPlusJets_M10To50_aMCatNLO_ext";

		switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                }

                break;

	case 9: Process = "ZPlusJets_M10To50_Madgraph";

		switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                }

                break;

	case 10: Process = "ZPlusJets_M10To50_Madgraph_ext";

                switch(YearInt){
                        case 2016: break;
                        case 2017: break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                }

                break;

	case 11: Process = "ZPlusJets_PtBinned_0To50";

		switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                }

		break;

	case 12: Process = "ZPlusJets_PtBinned_50To100";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 13: Process = "ZPlusJets_PtBinned_50To100_ext";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 14: Process = "ZPlusJets_PtBinned_100To250";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 15: Process = "ZPlusJets_PtBinned_100To250_ext1";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 16: Process = "ZPlusJets_PtBinned_100To250_ext2";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 17: Process = "ZPlusJets_PtBinned_100To250_ext5";
		
		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 18: Process = "ZPlusJets_PtBinned_250To400";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 19: Process = "ZPlusJets_PtBinned_250To400_ext1";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 20: Process = "ZPlusJets_PtBinned_250To400_ext2";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 21: Process = "ZPlusJets_PtBinned_250To400_ext5";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 22: Process = "ZPlusJets_PtBinned_400To650";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 23: Process = "ZPlusJets_PtBinned_400To650_ext1";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 24: Process = "ZPlusJets_PtBinned_400To650_ext2";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 25: Process = "ZPlusJets_PtBinned_650ToInf";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 26: Process = "ZPlusJets_PtBinned_650ToInf_ext1";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 27: Process = "ZPlusJets_PtBinned_650ToInf_ext2";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 28: Process = "ttbar_2l2nu";

		 switch(YearInt){
                        case 2016: HessianOrMC = "Hessian"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 29: Process = "ttbar_madgraph";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break; 
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 30: Process = "ttbar_madgraph_ext";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 31: Process = "ttbar_TTToHadronic";
	
		 switch(YearInt){
                        case 2016: HessianOrMC = "Hessian"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 32: Process = "ttbar_TTToHadronic_ext";

                 switch(YearInt){
                        case 2016: break;
                        case 2017: break; 
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 33: Process = "ttbar_TTToSemileptonic";

		 switch(YearInt){
                        case 2016: HessianOrMC = "Hessian"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 34: Process = "ttbar_TTToSemileptonic_ext";

                 switch(YearInt){
                        case 2016: break;
                        case 2017: break; 
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 35: Process = "ttbar_atMCaNLO";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;
	
	case 36: Process = "ttbar_atMCaNLO_ext";

                 switch(YearInt){
                        case 2016: break; 
                        case 2017: break;
			case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;


        case 37: Process = "ttbar_inc";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 38: Process = "SingleTop_tchannel_top";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 39: Process = "SingleTop_tchannel_top_ScaleUp";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 40: Process = "SingleTop_tchannel_top_ScaleDown";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 41: Process = "SingleTop_tchannel_antitop";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 42: Process = "SingleTop_tchannel_antitop_ScaleUp";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 43: Process = "SingleTop_tchannel_antitop_ScaleDown";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 44: Process = "SingleTop_schannel";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 45: Process = "ttbar_hdampUP";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 46: Process = "ttbar_hdampUP_ext";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 47: Process = "ttbar_hdampDOWN";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 48: Process = "ttbar_hdampDOWN_ext";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 49: Process = "TT_2l2nu_hdampUP";

                 switch(YearInt){
                        case 2016: HessianOrMC = "Hessian"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 50: Process = "TT_2l2nu_hdampUP_ext1";

                 switch(YearInt){
                        case 2016: HessianOrMC = "Hessian"; break;
                        case 2017: break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }
         
                 break;

	case 51: Process = "TT_2l2nu_hdampUP_ext2";

                 switch(YearInt){
                        case 2016: HessianOrMC = "Hessian"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }
         
                 break;

	case 52: Process = "TT_2l2nu_hdampDOWN";

                 switch(YearInt){
                        case 2016: HessianOrMC = "Hessian"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }
         
                 break;
	
	case 53: Process = "TT_2l2nu_hdampDOWN_ext1";

                 switch(YearInt){
                        case 2016: HessianOrMC = "Hessian"; break;
                        case 2017: break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;
	
	case 54: Process = "TT_2l2nu_hdampDOWN_ext2";

                 switch(YearInt){
                        case 2016: HessianOrMC = "Hessian"; break;
                        case 2017: break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;


	case 55: Process = "TTToHadronic_hdampUP";

                 switch(YearInt){
                        case 2016: HessianOrMC = "Hessian"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 56: Process = "TTToHadronic_hdampDOWN";

                 switch(YearInt){
                        case 2016: HessianOrMC = "Hessian"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;
	
	case 57: Process = "TTToSemileptonic_hdampUP";

                 switch(YearInt){
                        case 2016: HessianOrMC = "Hessian"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 58: Process = "TTToSemileptonic_hdampDOWN";

                 switch(YearInt){
                        case 2016: HessianOrMC = "Hessian"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 59: Process = "TTToSemileptonic_hdampDOWN_ext";

                 switch(YearInt){
                        case 2016: break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;


        case 60: Process = "SingleTop_tchannel_top_hdampUP";

		 switch(YearInt){
                        case 2016: HessianOrMC = "Hessian"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 61: Process = "SingleTop_tchannel_top_hdampDOWN";

		 switch(YearInt){
                        case 2016: HessianOrMC = "Hessian"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
			default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

    	case 62: Process = "SingleTop_tchannel_antitop_hdampUP";

                 switch(YearInt){
                        case 2016: HessianOrMC = "Hessian"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 63: Process = "SingleTop_tchannel_antitop_hdampDOWN";

                 switch(YearInt){
                        case 2016: HessianOrMC = "Hessian"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 64: Process = "ttbar_isr_UP";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 65: Process = "ttbar_isr_DOWN";
		
		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 66: Process = "ttbar_isr_DOWN_ext";
	
		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 67: Process = "ttbar_fsr_UP";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 68: Process = "ttbar_fsr_UP_ext";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 69: Process = "ttbar_fsr_DOWN";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 70: Process = "ttbar_fsr_DOWN_ext";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 71: Process = "SingleTop_tW";

		 switch(YearInt){
                        case 2016: HessianOrMC = "Hessian"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 72: Process = "SingleTop_tW_ScaleUp";
	
		 switch(YearInt){
                        case 2016: HessianOrMC = "Hessian"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 73: Process = "SingleTop_tW_ScaleDown";

		 switch(YearInt){
                        case 2016: HessianOrMC = "Hessian"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 74: Process = "SingleTop_tbarW";

		 switch(YearInt){
                        case 2016: HessianOrMC = "Hessian"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 75: Process = "SingleTop_tbarW_ScaleUp";

		 switch(YearInt){
                        case 2016: HessianOrMC = "Hessian"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 76: Process = "SingleTop_tbarW_ScaleDown";

		 switch(YearInt){
                        case 2016: HessianOrMC = "Hessian"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 77: Process = "SingleTop_tHq";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 78: Process = "SingleTop_tZq_W_lept_Z_had";

		 switch(YearInt){
                        case 2016: break;
                        case 2017: HessianOrMC = "MC"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }
	
                 break;

        case 79: Process = "SingleTop_tWZ_tWll";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: HessianOrMC = "MC"; break;
                        case 2018: HessianOrMC = "MC"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 80: Process = "VV_ZZTo2l2nu";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 81: Process = "VV_ZZTo2l2nu_ext1";

		 switch(YearInt){
                        case 2016: break;
                        case 2017: break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 82: Process = "VV_ZZTo2l2nu_ext2";

                 switch(YearInt){
                        case 2016: break;
                        case 2017: break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 83: Process = "VV_ZZTo2l2Q";
	
		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 84: Process = "VV_ZZTo4L";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 85: Process = "VV_ZZTo4L_ext";

                 switch(YearInt){
                        case 2016: break;
                        case 2017: break; 
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;
	

        case 86: Process = "VV_WZTo2l2Q";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 87: Process = "VV_WZTo3lNu";

		 switch(YearInt){
                        case 2016: break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 88: Process = "VV_WZTo3lNu_ext";

                 switch(YearInt){
                        case 2016: break;
                        case 2017: break; 
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;


        case 89: Process = "VV_WZTo1l1Nu2Q";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 90: Process = "VV_WWTo2l2Nu";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 91: Process = "VV_WWToLNuQQ";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 92: Process = "VV_WWToLNuQQ_ext";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 93: Process = "VV_WGToLNuG";

		 switch(YearInt){
                        case 2016: break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 94: Process = "VV_ZGToLLG";

		 switch(YearInt){
                        case 2016: break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 95: Process = "VV_ZGToLLG_ext";

                 switch(YearInt){
                        case 2016: break;
                        case 2017: break; 
			case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 96: Process = "VVV_WWWTo4F";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: HessianOrMC = "MC"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 97: Process = "VVV_WWWTo4F_ext";

                 switch(YearInt){
                        case 2016: break; 
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: break;
			default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 98: Process = "VVV_WWZTo4F";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: HessianOrMC = "MC"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 99: Process = "VVV_WWZTo4F_ext";

                 switch(YearInt){
                        case 2016: break; 
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: break;
			default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 100: Process = "VVV_WZZ";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: HessianOrMC = "MC"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 101: Process = "VVV_WZZ_ext";

                 switch(YearInt){
                        case 2016: break;
			case 2017: HessianOrMC = "MC"; break;
                        case 2018: break;
			default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 102: Process = "VVV_ZZZ";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: HessianOrMC = "MC"; break;
                        case 2018: HessianOrMC = "MC"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 103: Process = "VVV_ZZZ_ext";

                 switch(YearInt){
                        case 2016: break;
			case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: break;
			default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 104: Process = "WPlusJets";
		
		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 105: Process = "WPlusJets_ext";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 106: Process = "ttbarV_ttWJetsToLNu";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 107: Process = "ttbarV_ttWJetsToLNu_ext";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 108: Process = "ttbarV_ttWJetsToQQ";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 109: Process = "ttbarV_ttZToLL";

		 switch(YearInt){
                        case 2016: break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 110: Process = "ttbarV_ttZToLL_ext";

                 switch(YearInt){
                        case 2016: break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: break;
			default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 111: Process = "ttbarV_ttZToLLNuNu";

                 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 112: Process = "ttbarV_ttZToLLNuNu_ext";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 113: Process = "ttbarV_ttZToLLNuNu_ext2";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: break;
                        case 2018: break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;
	
	case 114: Process = "ttbarV_ttZToQQ";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

        case 115: Process = "ttbarV_ttZToQQ_ext";

		 switch(YearInt){
                        case 2016: break; 
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 116: Process = "ttbarV_ttHTobb";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: HessianOrMC = "Hessian"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 117: Process = "ttbarV_ttHTobb_ext";

                 switch(YearInt){
                        case 2016: break;
                        case 2017: break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 118: Process = "ttbarV_ttHToNonbb";

		 switch(YearInt){
                        case 2016: HessianOrMC = "MC"; break;
                        case 2017: HessianOrMC = "MC"; break;
                        case 2018: HessianOrMC = "Hessian"; break;
                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018"); break;
                 }

                 break;

	case 119: Process = "TriggerSF_MC"; break;

	case 120: Process = "TriggerSF_DATA"; break;

	case 121: Process = "Data_DoubleEGRunB"; break;

	case 122: Process = "Data_DoubleEGRunC"; break;

        case 123: Process = "Data_DoubleEGRunD"; break;

	case 124: Process = "Data_DoubleEGRunE"; break;


	case 125: Process = "Data_DoubleEGRunF"; break;

	case 126: Process = "Data_DoubleEGRunG"; break;

	case 127: Process = "Data_DoubleEGRunH"; break;

	case 128: Process = "Data_DoubleMuonRunB"; break;		

	case 129: Process = "Data_DoubleMuonRunC"; break;

	case 130: Process = "Data_DoubleMuonRunD"; break;

	case 131: Process = "Data_DoubleMuonRunE"; break;

	case 132: Process = "Data_DoubleMuonRunF"; break;

	case 133: Process = "Data_DoubleMuonRunG"; break;

	case 134: Process = "Data_DoubleMuonRunH"; break;

	case 135: Process = "Data_MuonEGRunB"; break;

	case 136: Process = "Data_MuonEGRunC"; break;
	
	case 137: Process = "Data_MuonEGRunD"; break;

        case 138: Process = "Data_MuonEGRunE"; break;

	case 139: Process = "Data_MuonEGRunF"; break;

	case 140: Process = "Data_MuonEGRunG"; break;

	case 141: Process = "Data_MuonEGRunH"; break;

	case 142: Process = "Data_SingleMuonRunB"; break;

	case 143: Process = "Data_SingleMuonRunC"; break;

	case 144: Process = "Data_SingleMuonRunD"; break;

	case 145: Process = "Data_SingleMuonRunE"; break;

	case 146: Process = "Data_SingleMuonRunF"; break;

	case 147: Process = "Data_SingleMuonRunG"; break;

	case 148: Process = "Data_SingleMuonRunH"; break;
 
	case 149: Process = "Data_SingleElectronRunB"; break;

	case 150: Process = "Data_SingleElectronRunC"; break;

	case 151: Process = "Data_SingleElectronRunD"; break;

	case 152: Process = "Data_SingleElectronRunE"; break;

	case 153: Process = "Data_SingleElectronRunF"; break;

	case 154: Process = "Data_SingleElectronRunG"; break;

	case 155: Process = "Data_SingleElectronRunH"; break;

	case 156: Process = "Data_EGRunB"; break;

	case 157: Process = "Data_EGRunC"; break;

	case 158: Process = "Data_EGRunD"; break;

	case 159: Process = "Data_Combined"; break;

	default: throw std::logic_error("Double check the value of ProcessInt."); break;

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

	default: throw std::logic_error("ChannelInt must be 1 (for ee), 2 (for mumu) or 3 (for emu)."); break;

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



  switch(YearInt){
                
  	case 2016: PSWeightString = "LeptonPt"; 

		   break;

        case 2017: switch(ProcessInt){ 
                   	case 0: PSWeightString = "PSWeight"; break; //tZq
			case 28: PSWeightString = "PSWeight"; break; //ttbar (2l2nu)
			case 31: PSWeightString = "PSWeight"; break; //ttbar (to hadronic)
			case 32: PSWeightString = "PSWeight"; break; //ttbar (to hadronic ext)
			case 33: PSWeightString = "PSWeight"; break; //ttbar (to semileptonic)
			case 34: PSWeightString = "PSWeight"; break; //ttbar (to semileptonic ext)
                        case 38: PSWeightString = "PSWeight"; break; //single top t-channel (top)
			case 41: PSWeightString = "PSWeight"; break; //single top t-channel (antitop)
			case 44: PSWeightString = "PSWeight"; break; //single top s-channel
			case 71: PSWeightString = "PSWeight"; break; //single top (tW)
			case 74: PSWeightString = "PSWeight"; break; //single top (tbarW)
			case 90: PSWeightString = "PSWeight"; break; //WWTo2L2Nu
			case 106: PSWeightString = "PSWeight"; break; //ttWToLNu
                        case 107: PSWeightString = "PSWeight"; break; //ttWToLNu_ext
			case 114: PSWeightString = "PSWeight"; break; //ttZToQQ
                        case 115: PSWeightString = "PSWeight"; break; //ttZToQQ_ext
                        default: PSWeightString = "LeptonPt"; break;
                    }
	
		   
		   break;

        case 2018: switch(ProcessInt){ 
			case 0: PSWeightString = "PSWeight"; break; //tZq
			case 28: PSWeightString = "PSWeight"; break; //ttbar (2l2nu)
                        case 31: PSWeightString = "PSWeight"; break; //ttbar (to hadronic)
                        case 32: PSWeightString = "PSWeight"; break; //ttbar (to hadronic ext)
                        case 33: PSWeightString = "PSWeight"; break; //ttbar (to semileptonic)
                        case 34: PSWeightString = "PSWeight"; break; //ttbar (to semileptonic ext)
                        case 38: PSWeightString = "PSWeight"; break; //single top t-channel (top)
                        case 41: PSWeightString = "PSWeight"; break; //single top t-channel (antitop)
                        case 44: PSWeightString = "PSWeight"; break; //single top s-channel
                        case 71: PSWeightString = "PSWeight"; break; //single top (tW)
                        case 74: PSWeightString = "PSWeight"; break; //single top (tbarW)
			case 90: PSWeightString = "PSWeight"; break; //WWTo2L2Nu
			case 106: PSWeightString = "PSWeight"; break; //ttWToLNu
                        case 107: PSWeightString = "PSWeight"; break; //ttWToLNu_ext
			case 114: PSWeightString = "PSWeight"; break; //ttZToQQ
                        case 115: PSWeightString = "PSWeight"; break; //ttZToQQ_ext
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


  	std::cout << "print 1" << std::endl;

	return  Flag_goodVertices_Selection > 0       	     || Flag_globalSuperTightHalo2016Filter_Selection > 0     || Flag_HBHENoiseFilter_Selection > 0 || 
		Flag_HBHENoiseIsoFilter_Selection > 0 	     || Flag_EcalDeadCellTriggerPrimitiveFilter_Selection > 0 || Flag_BadPFMuonFilter_Selection > 0 || 
		Flag_BadChargedCandidateFilter_Selection > 0 || Flag_ecalBadCalibFilter_Selection > 0 		      || Flag_eeBadScFilter_Selection > 0;

  }};



  //Lambda functions for filtering data events using the golden json file. This is to remove events with bad lumisections.
  std::vector<std::string> GoldenJsonOutput{};
  std::vector<char> GoldenJson_SplitChars_Output{};
  std::vector<char> EventsVector{};
  std::vector<int> RunNumAndEvents{};
  std::vector<int> BadRuns{};
  std::vector<int> GoodRuns{};

  //Reading the golden json file and storing its contents into the vector GoldenJsonOutput
  auto GoldenJsonReader{[&YearInt, &GoldenJsonOutput](){

	std::cout << "print 2" << std::endl;

	if(GoldenJsonOutput.size() == 0){

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
  
   		if (myReadFile.is_open()) { while (!myReadFile.eof()) {myReadFile >> output; GoldenJsonOutput.push_back(output);}}

  		myReadFile.close();
	
	}

	return 0;

  }};


  //Iterating over the contents of the GoldenJsonOutput vector to save each individual character as its own element in the vector GoldenJson_SplitChars_Output
  auto GoldenJson_SplitChars{[&YearInt, &GoldenJsonReader, &GoldenJsonOutput, &GoldenJson_SplitChars_Output](){

	std::cout << "print 3" << std::endl;

	if(GoldenJson_SplitChars_Output.size() == 0){

		GoldenJsonReader();

    		for(long unsigned int i = 0; i < GoldenJsonOutput.size(); i++){

  			std::string element = GoldenJsonOutput.at(i);
  	  		for(long unsigned int j = 0; j < element.size(); j++){GoldenJson_SplitChars_Output.push_back(element.at(j));}

    		}

	}

	return 0;

  }};
 

  auto RunNumberCheck{[&GoodRuns, &YearInt, &GoldenJsonOutput, &GoldenJson_SplitChars, &GoldenJson_SplitChars_Output, &EventsVector, &BadRuns](const unsigned int& InputRunNumber){

	std::cout << "print 4" << std::endl;

	GoldenJson_SplitChars();

	bool GoodRunsCheck = any_of(GoodRuns.begin(), GoodRuns.end(), [&InputRunNumber](int i){return i == InputRunNumber;}); //check to see if the input run number is already in the GoodRuns vector

	if(GoodRunsCheck == true){return 0;}
	

  	for(long unsigned int i = 0; i < GoldenJson_SplitChars_Output.size(); i++){

		std::cout << "i = " << i << " . GoldenJson_SplitChars_Output.size() = " <<  GoldenJson_SplitChars_Output.size() << std::endl;
	
		 if(i == (GoldenJson_SplitChars_Output.size() - 1)){

                        std::cout << "i is equal to GoldenJson_SplitChars_Output.size() - 1. i =  " << i  
				  << " GoldenJson_SplitChars_Output.size() = " << GoldenJson_SplitChars_Output.size() << std::endl;

                        bool BadRunsCheck = any_of(BadRuns.begin(), BadRuns.end(), [&InputRunNumber](int i){return i == InputRunNumber;}); //check to see if the input run number is already in the BadRuns vector
          
                        if(BadRunsCheck == true){std::cout << "BadRunsCheck is true. BadRuns.size() = " << BadRuns.size() << std::endl; return 0;}
                        else{
                                std::cout << "No matches. Adding " << InputRunNumber << " to the BadRuns vector." << std::endl;
                                BadRuns.push_back(InputRunNumber);
                                return 0;

                        }
    
                }

	
		unsigned int RunNumBeingRead;

 		if(  GoldenJson_SplitChars_Output.at(i+1) == '"' && (GoldenJson_SplitChars_Output.at(i+2) == '2' || GoldenJson_SplitChars_Output.at(i+2) == '3')  ){ 

			int digit1 = GoldenJson_SplitChars_Output.at(i+2) - '0';	
			int digit2 = GoldenJson_SplitChars_Output.at(i+3) - '0';
			int digit3 = GoldenJson_SplitChars_Output.at(i+4) - '0'; 
			int digit4 = GoldenJson_SplitChars_Output.at(i+5) - '0';
			int digit5 = GoldenJson_SplitChars_Output.at(i+6) - '0';
			int digit6 = GoldenJson_SplitChars_Output.at(i+7) - '0';

			int run = (digit1*100000) + (digit2*10000) + (digit3*1000) + (digit4*100) + (digit5*10) + digit6;
			RunNumBeingRead = run;

			if(run == InputRunNumber){

				std::cout << '\n' << std::endl;
				std::cout << '\n' << std::endl;
				std::cout << "run " << RunNumBeingRead << " is equal to the InputRunNumber of " << InputRunNumber << std::endl;
				std::cout << '\n' << std::endl;
				std::cout << '\n' << std::endl;

				for(long unsigned int j = 2; j < GoldenJson_SplitChars_Output.size(); j++){

					if(GoldenJson_SplitChars_Output.at(i+10) == '[' && GoldenJson_SplitChars_Output.at(i+11) == '['){	

						if( GoldenJson_SplitChars_Output.at( (i+10)+j ) == ']' && GoldenJson_SplitChars_Output.at( (i+10)+(j+1) ) == ']'){

							for(int k = (i+10); k < ((i+10)+(j+2)); k++){
								EventsVector.push_back(GoldenJson_SplitChars_Output.at(k));
							}

							return 0;	
											
						}
						else{continue;}

					}
					else{continue;}

				}	

			}
			else{continue;}

		}
		else{std::cout << "The run number of " << RunNumBeingRead << " does not match the input run number of " << InputRunNumber << std::endl;
		     continue;}


   	}


  }};


  //Function to save the event numbers that match the run to the RunNumAndEvents vector
  auto ReturnRunNumAndEventRanges{[&GoodRuns, &YearInt, &RunNumberCheck, &RunNumAndEvents, &GoldenJson_SplitChars_Output, &GoldenJsonOutput, &EventsVector, &BadRuns](const unsigned int& InputRunNumber){

	std::cout << "print 5" << std::endl;

	RunNumberCheck(InputRunNumber);

	bool BadRunsCheck2 = any_of(BadRuns.begin(), BadRuns.end(), [&InputRunNumber](int i){return i == InputRunNumber;}); //check to see if the input run number is already in the BadRuns vector

	bool GoodRunsCheck2 = any_of(GoodRuns.begin(), GoodRuns.end(), [&InputRunNumber](int i){return i == InputRunNumber;}); //check to see if the input run number is already in the GoodRuns vector

	std::cout << '\n' << std::endl;
	std::cout << '\n' << std::endl;
	std::cout << "inside ReturnRunNumAndEventRanges" << std::endl;
	std::cout << "BadRunsCheck2 = " << BadRunsCheck2  << std::endl;
	std::cout << "GoodRunsCheck2 = " << GoodRunsCheck2 << std::endl;
	std::cout << '\n' << std::endl;
	std::cout << '\n' << std::endl;	


        if(BadRunsCheck2 == true){std::cout << "BadRunsCheck is true. BadRuns.size() = " << BadRuns.size() << std::endl; return 0;}
	else if(GoodRunsCheck2 == true){return 0;}
	else{
 
		std::cout << "inside the else for RunNumCheck" << std::endl;

		GoodRuns.push_back(InputRunNumber);

  		RunNumAndEvents.push_back(InputRunNumber);
   		std::vector<char> Runs = EventsVector;

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
				else{std::cout << "INSIDE THE ELSE STATEMENT." << std::endl;}


		}	 


    	}

   }


    EventsVector.clear();
    return 0;

  }};


 
  auto RunAndLumiFilterFunction{[&ReturnRunNumAndEventRanges, &MCInt, &RunNumAndEvents, &BadRuns, &GoodRuns](const unsigned int& InputRunNumber, const unsigned int& luminosityBlock){

     std::cout << "print 6" << std::endl;	

     bool BadRunsCheck3 = any_of(BadRuns.begin(), BadRuns.end(), [&InputRunNumber](int i){return i == InputRunNumber;}); //check to see if the input run number is already in the BadRuns vector

     bool GoodRunsCheck3 = any_of(GoodRuns.begin(), GoodRuns.end(), [&InputRunNumber](int i){return i == InputRunNumber;}); //check to see if the input run number is already in the GoodRuns vector

     bool RunNumAndEventsCheck = any_of(RunNumAndEvents.begin(), RunNumAndEvents.end(), [&InputRunNumber](int i){return i == InputRunNumber;});    

      switch(MCInt){
      
          case 0: 
     		  if(BadRunsCheck3 == true){std::cout << "InputRunNumber = " << InputRunNumber << " BadRunsCheck is true. BadRuns.size() = " << BadRuns.size() << std::endl; return false;}

     		  if(GoodRunsCheck3 == false){ReturnRunNumAndEventRanges(InputRunNumber);}

		  if(RunNumAndEvents.size() == 0){return false;}
		  else if(RunNumAndEventsCheck == true){

//				for(int i = 0; i < RunNumAndEvents.size(); i++){std::cout << "RunNumAndEvents.at(i) = " << RunNumAndEvents.at(i) << std::endl;}    
			
				std::cout << "InputRunNumber = " << InputRunNumber << std::endl;
	
				std::vector<int>::iterator it_RunNumAndEvents = std::find(RunNumAndEvents.begin(), RunNumAndEvents.end(), InputRunNumber);
				int index_RunNumAndEvents = std::distance(RunNumAndEvents.begin(), it_RunNumAndEvents);
				int DistToEnd = std::distance(it_RunNumAndEvents, RunNumAndEvents.end());

				ints FinalRunNumAndEventsVec{};

				std::cout << "index_RunNumAndEvents = " << index_RunNumAndEvents << std::endl;
				std::cout << "DistToEnd = " << DistToEnd << std::endl;

				for(int i = index_RunNumAndEvents; i < (index_RunNumAndEvents + DistToEnd); i++){FinalRunNumAndEventsVec.push_back(RunNumAndEvents.at(i));}

//				std::cout << "FinalRunNumAndEventsVec.size() = " << FinalRunNumAndEventsVec.size() << std::endl;

//				for(int i = 0; i < FinalRunNumAndEventsVec.size(); i++){std::cout << "FinalRunNumAndEventsVec.at(i) = " << FinalRunNumAndEventsVec.at(i) << std::endl;}
					
				for(int i = 0; i < FinalRunNumAndEventsVec.size()-2; i+=2){

					if(luminosityBlock > FinalRunNumAndEventsVec.at(i+1) && luminosityBlock < FinalRunNumAndEventsVec.at(i+2)){

						std::cout << "luminosityBlock of " << luminosityBlock << " is between " << FinalRunNumAndEventsVec.at(i+1) << " and " 
							  << FinalRunNumAndEventsVec.at(i+2) << ". Returning true." << std::endl;

						return true;
					}
					else{std::cout << "luminosityBlock of " << luminosityBlock << " is not between " << FinalRunNumAndEventsVec.at(i+1) << " and " 
                                                       << FinalRunNumAndEventsVec.at(i+2) << ". Continuing." << std::endl; continue;}

				}
			
				//throw std::logic_error("break");
				std::cout << "returning false" << std::endl;
				return false;	

                 }
                 else{return false;}


		 break;
              
          case 1: return  InputRunNumber && luminosityBlock; break;
              
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


	std::cout << "print 7" << std::endl;

	switch(TriggerType){

		case 0: switch(YearInt){//single electron
						case 2016: return (HLT_Ele25_eta2p1_WPTight_Gsf > 0 || HLT_Ele27_WPTight_Gsf > 0 || HLT_Ele32_eta2p1_WPTight_Gsf > 0);
						case 2017: return (HLT_Ele32_WPTight_Gsf_L1DoubleEG > 0 || HLT_Ele35_WPTight_Gsf > 0);
						case 2018: return (HLT_Ele32_WPTight_Gsf_L1DoubleEG > 0 || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL > 0);
						default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018."); break;
					}

		case 1: switch(YearInt){//single muon
                                   	case 2016: return HLT_IsoMu24 <= 0;/* || HLT_IsoMu24_eta2p1 <= 0);*/ //the commented out one is not in 2016
                                        case 2017: return (HLT_IsoMu27 <= 0);
                                        case 2018: return HLT_IsoMu24 <= 0;
                                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018."); break; 
                                   }

		case 2: switch(YearInt){//double muon
                                        case 2016: return HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ <= 0;
                                        case 2017: return (HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ <= 0 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 <= 0);
                                        case 2018: return (HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 <= 0 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8 <= 0);
					default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018."); break;
                                   }

		case 3: switch(YearInt){//double electron
                                        case 2016: return HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ > 0;
                                        case 2017: return (HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL > 0 || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ > 0);
                                        case 2018: return HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ > 0;
					default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018."); break;
                                   }

		case 4: switch(YearInt){//muon EG
                                        case 2016: return //HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || not in 2016 
                                            	   	  //HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || not in all 2016 runs 
                                            		  //HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || not in 2016 
                                            		  HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL <= 0; //|| 
                                            		  //HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL <= 0 || not in 2016 
                                            		  //HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL <= 0); //not in 2016 run H

                                        case 2017: return (HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || 
                                            		   HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0 ||
                                            		   HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0); 
                                        
					case 2018: return (HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || 
                                           		   HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || 
                                           		   HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0); 

                                        default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018."); break;

                                   }

		default: throw std::logic_error("Choose a case out of SingleElectron, SingleMuon, DoubleElectron, DoubleMuon or MuonEG."); break; 

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


	std::cout << "print 8" << std::endl;

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
  auto PU_function{[&puReweight_2016, &puReweight_2016_part1, &puReweight_2016_part2, &puReweight_2017, &puReweight_2018, 
		    &puSystUp_2016,  &puSystDown_2016, &puSystUp_2017,  &puSystDown_2017, &puSystUp_2018,  &puSystDown_2018, 
	            &YearInt, &SystematicInt](int PV_npvs_input){

  	std::cout << "print 9 PU_function" << std::endl;

      	double PU_Weight_input;

      	switch(YearInt){

        	case 2016: switch(SystematicInt){
				case 1: PU_Weight_input = puSystUp_2016->GetBinContent(puSystUp_2016->GetXaxis()->FindBin(PV_npvs_input)); break;
			        case 2: PU_Weight_input = puSystDown_2016->GetBinContent(puSystDown_2016->GetXaxis()->FindBin(PV_npvs_input)); break;
				default: PU_Weight_input = puReweight_2016->GetBinContent(puReweight_2016->GetXaxis()->FindBin(PV_npvs_input)); break;
			   }

			   break;

          	case 2017: switch(SystematicInt){
				case 1: PU_Weight_input = puSystUp_2017->GetBinContent(puSystUp_2017->GetXaxis()->FindBin(PV_npvs_input)); break;
                                case 2: PU_Weight_input = puSystDown_2017->GetBinContent(puSystDown_2017->GetXaxis()->FindBin(PV_npvs_input)); break;
				default: PU_Weight_input = puReweight_2017->GetBinContent(puReweight_2017->GetXaxis()->FindBin(PV_npvs_input)); break;

			   }

			   break;

          	case 2018: switch(SystematicInt){
				case 1: PU_Weight_input = puSystUp_2018->GetBinContent(puSystUp_2018->GetXaxis()->FindBin(PV_npvs_input)); break;
                                case 2: PU_Weight_input = puSystDown_2018->GetBinContent(puSystDown_2018->GetXaxis()->FindBin(PV_npvs_input)); break;
			   	default: PU_Weight_input = puReweight_2018->GetBinContent(puReweight_2018->GetXaxis()->FindBin(PV_npvs_input)); break;
			   }

			   break;

		default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018."); break;

      	}

      	return PU_Weight_input;

  }};

  //Lambda functions for the electron selection
  auto ElectronsFunction{[](const int targetID, const floats& Electron_pt, const floats& Electron_eta, const ints& Electron_cutBased, const bools& Electron_isPFcand){
 
  	std::cout << "print 9" << std::endl;
  	return (Electron_pt > MinElectronPt && (abs(Electron_eta) < MaxTrackerEta && (abs(Electron_eta) < 1.442 || abs(Electron_eta) > 1.566) ) && 
		Electron_cutBased >= targetID && Electron_isPFcand);

  }};

  auto MuonsFunction{[](const float target_iso, const bools& isPFs, const floats& Muon_pt, const floats& Muon_eta, const bools& ids, const floats& isos){

  	std::cout << "print 10" << std::endl;
  	return (isPFs && Muon_pt > MinMuonPt && abs(Muon_eta) < MaxTrackerEta && ids && isos <= target_iso);

  }};

  auto EmuFunction{[&ElectronsFunction, &MuonsFunction](const int targetID,     const floats& Electron_pt, const floats& Electron_eta, const ints& Electron_cutBased, const bools& Electron_isPFcand,
		      const float target_iso, const bools& isPFs,        const floats& Muon_pt,      const floats& Muon_eta,        const bools& ids, const floats& isos){

	std::cout << "print 10 part 2" << std::endl;

	return ElectronsFunction(targetID, Electron_pt, Electron_eta, Electron_cutBased, Electron_isPFcand).size() == 1 && MuonsFunction(target_iso, isPFs, Muon_pt, Muon_eta, ids, isos).size() == 1;

  }};

  auto EmuFunctionCount{[&EmuFunction](const int targetID,     const floats& Electron_pt, const floats& Electron_eta, const ints& Electron_cutBased, const bools& Electron_isPFcand,
                      const float target_iso, const bools& isPFs,        const floats& Muon_pt,      const floats& Muon_eta,        const bools& ids, const floats& isos){

	int n_emu = EmuFunction(targetID, Electron_pt, Electron_eta, Electron_cutBased, Electron_isPFcand, target_iso, isPFs, Muon_pt, Muon_eta, ids, isos);
	
	int EmuVecSize = Electron_pt.size() + Muon_pt.size();

	ints emu_out;
	
  	for(int i = 0; i < EmuVecSize; i++){emu_out.push_back(n_emu);}

	return emu_out;

  }};

  auto TightLeptonsFunction{[&ChannelInt, &ElectronsFunction, &MuonsFunction, &EmuFunctionCount](const floats& Electron_pt,      const floats& Electron_eta, const ints& Electron_cutBased, 
                                                                                            const bools& Electron_isPFcand, const bools& isPFs,         const floats& pts, 
                                                                                            const floats& etas,             const bools& ids,           const floats& isos){

        std::cout << "print 11" << std::endl;

        switch(ChannelInt){

                case 1: return ElectronsFunction(4, Electron_pt, Electron_eta, Electron_cutBased, Electron_isPFcand);

                case 2: return MuonsFunction(0.25, isPFs, pts, etas, ids, isos);

                case 3: return EmuFunctionCount(4, Electron_pt, Electron_eta, Electron_cutBased, Electron_isPFcand, 0.25, isPFs, pts, etas, ids, isos); 

                default: throw std::logic_error("ChannelInt must be 1 (for ee), 2 (for mumu) or 3 (for emu)."); break;

        }

  }};

  auto LeptonVariableFunctionFloats{[&ChannelInt](const floats& Electron_input, const floats& Muon_input){

	std::cout << "print 12" << std::endl; 

	floats Emu_vector_floats{};

	switch(ChannelInt){

		case 1: return Electron_input; 
		case 2: return Muon_input;

		case 3: if(Electron_input.size() > 0){

				for(int i = 0; i < Electron_input.size(); i++){
					Emu_vector_floats.push_back(Electron_input.at(i));
				}	

			}

			if(Muon_input.size() > 0){
			
				for(int i = 0; i < Muon_input.size(); i++){
					Emu_vector_floats.push_back(Muon_input.at(i));
				}

			}

			return Emu_vector_floats;

		default: throw std::logic_error("ChannelInt must be 1 (for ee), 2 (for mumu) or 3 (for emu)."); break;

	}

  }};

  auto LeptonVariableFunctionInts{[&ChannelInt](const ints& Electron_input, const ints& Muon_input){

        std::cout << "print 13" << std::endl;

	ints Emu_vector_ints{};

        switch(ChannelInt){

                case 1: return Electron_input;
                case 2: return Muon_input;

                case 3: if(Electron_input.size() > 0){

                                for(int i = 0; i < Electron_input.size(); i++){
                                        Emu_vector_ints.push_back(Electron_input.at(i));
                                }

                        }

                        if(Muon_input.size() > 0){

                                for(int i = 0; i < Muon_input.size(); i++){
                                        Emu_vector_ints.push_back(Muon_input.at(i));
                                }

                        } 

			return Emu_vector_ints;

			default: throw std::logic_error("ChannelInt must be 1 (for ee), 2 (for mumu) or 3 (for emu)."); break;

        }

  }};


  auto LeptonVariableFunctionChars{[&ChannelInt, &MCInt](const chars& Electron_input, const chars& Muon_input){

        std::cout << "print 14" << std::endl;

	chars LeptonVariableFunctionChars_Zeroes(2, 0);
	
	switch(MCInt){
		case 0: LeptonVariableFunctionChars_Zeroes; break;
		default: break;
	};

	chars Emu_vector_chars{};

        switch(ChannelInt){

                case 1: return Electron_input; break;
                case 2: return Muon_input; break;

                case 3: if(Electron_input.size() > 0){

                                for(int i = 0; i < Electron_input.size(); i++){
                                        Emu_vector_chars.push_back(Electron_input.at(i));
                                }

                        }

                        if(Muon_input.size() > 0){

                                for(int i = 0; i < Muon_input.size(); i++){
                                        Emu_vector_chars.push_back(Muon_input.at(i));
                                }

                        } 

                        return Emu_vector_chars;

			break;

		default: throw std::logic_error("ChannelInt must be 1 (for ee), 2 (for mumu) or 3 (for emu)."); break;

        }

  }};


  auto LooseLeptonsFunction{[&ChannelInt, &ElectronsFunction, &MuonsFunction, &EmuFunctionCount](const floats& Electron_pt,      const floats& Electron_eta, const ints& Electron_cutBased, 
									                    const bools& Electron_isPFcand, const bools& isPFs,         const floats& pts, 
									                    const floats& etas,             const bools& ids,           const floats& isos){

  	std::cout << "print 15" << std::endl;

	switch(ChannelInt){

  		case 1: return ElectronsFunction(1, Electron_pt, Electron_eta, Electron_cutBased, Electron_isPFcand); 

		case 2: return MuonsFunction(0.15, isPFs, pts, etas, ids, isos);  

		case 3: return EmuFunctionCount(1, Electron_pt, Electron_eta, Electron_cutBased, Electron_isPFcand, 0.15, isPFs, pts, etas, ids, isos); 

		default: throw std::logic_error("ChannelInt must be 1 (for ee), 2 (for mumu) or 3 (for emu)."); break;

	}

  }};

  auto OppositeSign{[&ChannelInt](const ints& charges){

  	std::cout << "print 16" << std::endl;

	return charges.size() == 2 ? signbit(charges.at(0)) != signbit(charges.at(1)) : false;

  }};

  auto SameSign{[&ChannelInt](const ints& charges){

  	std::cout << "print 14" << std::endl;
	return charges.size() == 2 ? signbit(charges.at(0)) == signbit(charges.at(1)) : false;

  }};

  auto LeadingVariable{[&ChannelInt](const floats& variable){

  	std::cout << "print 17" << std::endl;
	
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

  	std::cout << "print 18" << std::endl;

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

  	std::cout << "print 19" << std::endl;

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

  	std::cout << "print 20" << std::endl;

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


	std::cout << "print 29" << std::endl;

	switch(ProcessInt){

		case 112: switch(YearInt){

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

		case 113: switch(YearInt){

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

	std::cout << "print 30" << std::endl;

	if(ProcessInt == 119 || 120){

	  switch(ChannelInt){

		case 1: switch(YearInt){
				case 2016: return //single or double electron and not any of the others
	
					   (HLT_Ele25_eta2p1_WPTight_Gsf > 0 || //single electron 
	 			 	    HLT_Ele27_WPTight_Gsf > 0 || //single electron
	 			 	    HLT_Ele32_eta2p1_WPTight_Gsf > 0 || //single electron
	 			 	    HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ > 0) //double electron 
	 				    && 
					   (HLT_IsoMu24 <= 0 || //single muon
				 	    //HLT_IsoMu24_eta2p1 <= 0 || //single muon //not in 2016
	 			 	    HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ <= 0 || //double muon
	 			 	    //HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 <= 0 || //double muon //not in 2016
	 			 	    //HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || //muon electron (not in 2016)
	 			 	    //HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || //muon electron (not in all 2016 runs)
	 			 	    //HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || //muon electron (not in 2016)
	 			 	    HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL <= 0 || //muon electron
	 			 	    //HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL <= 0 || //muon electron (not in 2016)
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
         				    //HLT_IsoMu24_eta2p1 > 0 || //single muon //not in 2016
         				    HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ > 0 //double muon
         				    //HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 > 0) //double muon ..not in 2016
					    )

					    &&

        				  (HLT_Ele25_eta2p1_WPTight_Gsf <= 0 || //single electron 
         				   HLT_Ele27_WPTight_Gsf <= 0 || //single electron
         				   HLT_Ele32_eta2p1_WPTight_Gsf <= 0 || //single electron
         				   HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || //double electron 
         				   //HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || //muon electron (not in 2016)
         				   //HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || //muon electron (not in all 2016 runs)
         				   //HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ <= 0 || //muon electron (not in 2016)
         			           HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL <= 0 || //muon electron
         				   //HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL <= 0 || //muon electron (not in 2016)
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
         				    //HLT_IsoMu24_eta2p1 > 0 || //single muon //not in 2016
	 				    HLT_Ele25_eta2p1_WPTight_Gsf > 0 || //single electron 
         				    HLT_Ele27_WPTight_Gsf > 0 || //single electron
         				    HLT_Ele32_eta2p1_WPTight_Gsf > 0 || //single electron
	 				    //HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ > 0 || //muon electron (not in 2016)
         				    //HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ > 0 || //muon electron (not in all 2016 runs)
         				    //HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ > 0 || //muon electron (not in 2016)
         				    HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL > 0 || //muon electron
         				    //HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL > 0 || //muon electron (not in 2016)
         				    HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL > 0) //muon electron

         				    &&

         				   (HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ <= 0 || //double muon
          				    //HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8 <= 0 || //double muon //not in 2016
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

                                default: throw std::logic_error("The year must be 2016, 2017 or 2018."); break;
                        }

		default: throw std::logic_error("ChannelInt must be 1 (for ee), 2 (for mumu) or 3 (for emu)."); break;


	  }

	}
	else{return event > 0;}

  }};


  auto NumberOfLeptonsFunction{[&ChannelInt](const unsigned int& nElectron, const unsigned int& nMuon){

	std::cout << "print NumberOfLeptonsFunction" << std::endl;

	switch(ChannelInt){
		case 1: return nElectron == 2;
		case 2: return nMuon == 2;
		case 3: return nElectron == 1 && nMuon == 1;
		default: throw std::logic_error("ChannelInt must be 1 (for ee), 2 (for mumu) or 3 (for emu)."); break;
	}

  }};

  auto Electron_dxy_dz_Function{[&ChannelInt](const floats& Electron_dz, const floats& Electron_dxy, const float& LeadingLeptonPt, const float& SubleadingLeptonPt, const floats& LeptonEta){


	std::cout << "print Electron_dxy_dz_Function" << std::endl;

	if(ChannelInt == 1){

		 return Electron_dz.size() > 0                   && Electron_dz.size() < 3               &&
                        abs(LeptonEta) < 1.442                   &&  abs(LeptonEta) > 1.566              &&
                        LeadingLeptonPt > MaxElectronPt          &&  SubleadingLeptonPt > MinElectronPt  &&
                        (abs(LeptonEta) < 1.442                  &&  Electron_dz < 0.1                   && Electron_dxy < 0.05)  || //barrel region
                        (abs(LeptonEta) > 1.566                  &&  abs(LeptonEta) < 3.0                && Electron_dz < 0.2 && Electron_dxy < 0.1); //endcaps

	}
	else{return abs(LeptonEta) < 3.0;} 
		

  }}; 



  auto LeptonCut{[&ChannelInt](const bool& os,                   const unsigned int& nElectron,    const unsigned int& nMuon,        
			       const floats& Electron_dz,        const floats& Electron_dxy,       const float& LeadingLeptonPt,     
			       const float& SubleadingLeptonPt,  const floats& LeptonEta,          const ints& Electron_dxy_dz,
			       const floats& tight_lepton_pts,   const floats& loose_lepton_pts){

  	std::cout << "print 31" << std::endl;

	const bool lepton_cut{tight_lepton_pts.size() == 2 && tight_lepton_pts.size() == loose_lepton_pts.size()};
  	bool lead_pt_cut{false};

	float MaxLeptonPt;

	switch(ChannelInt){
		case 1: MaxLeptonPt = MaxElectronPt; break;
		case 2: MaxLeptonPt = MaxMuonPt; break;
		case 3: MaxLeptonPt = MaxElectronPt; break; //both MaxMuonPt and MaxElectronPt are 25 for the emu channel
		default: throw std::logic_error("ChannelInt must be 1 (for ee), 2 (for mumu) or 3 (for emu)."); break;	
	}

	lead_pt_cut = tight_lepton_pts.empty() ? false : *max_element(tight_lepton_pts.begin(), tight_lepton_pts.end()) > MaxLeptonPt;

	bool Electron_dxy_dz_bool = any_of(Electron_dxy_dz.begin(), Electron_dxy_dz.end(), [&Electron_dxy_dz](int i = 0){return i > 0;});

	switch(ChannelInt){

		case 1: return os && lead_pt_cut && lepton_cut && nElectron == 2 && Electron_dxy_dz_bool;
		
		case 2: return os && lead_pt_cut && lepton_cut && nMuon == 2 && Electron_dxy_dz_bool;

		case 3: return os && lead_pt_cut && lepton_cut && nElectron == 1 && nMuon == 1 && Electron_dxy_dz_bool;

		default: throw std::logic_error("ChannelInt must be 1 (for ee), 2 (for mumu) or 3 (for emu)."); break;

	}


  }};
 
  auto OppositeSignNonPrompt{[](const ints& charges, const chars& Lepton_genPartFlav){

  	std::cout << "print 32" << std::endl;
	bool OppositeSignChargeCheck = charges.size() == 2 ? signbit(charges.at(0)) != signbit(charges.at(1)) : false;
  	bool LeptonNonPromptCheck = all_of(Lepton_genPartFlav.begin(), Lepton_genPartFlav.end(), [](int i){return i != 1;});

  	return OppositeSignChargeCheck && (LeptonNonPromptCheck == 1);

  }};

  auto OppositeSignPrompt{[](const ints& charges, const chars& Lepton_genPartFlav){

  	std::cout << "print 33" << std::endl;

	bool OppositeSignChargeCheck = charges.size() == 2 ? signbit(charges.at(0)) != signbit(charges.at(1)) : false;
  	bool LeptonPromptCheck = all_of(Lepton_genPartFlav.begin(), Lepton_genPartFlav.end(), [](int i){return i == 1;});

  	return OppositeSignChargeCheck && (LeptonPromptCheck == 1);

  }};


  auto SameSignNonPrompt{[](const ints& charges, const chars& Lepton_genPartFlav){
  
  	std::cout << "print 34" << std::endl;

	bool SameSignChargeCheck = charges.size() == 2 ? signbit(charges.at(0)) == signbit(charges.at(1)) : false;  
  	bool LeptonNonPromptCheck = all_of(Lepton_genPartFlav.begin(), Lepton_genPartFlav.end(), [](int i){return i != 1;});

  	return SameSignChargeCheck && (LeptonNonPromptCheck == 1);

  }};

  auto SameSignPrompt{[](const ints& charges, const chars& Lepton_genPartFlav){

  	std::cout << "print 35" << std::endl;  

	bool SameSignChargeCheck = charges.size() == 2 ? signbit(charges.at(0)) == signbit(charges.at(1)) : false;
  	bool LeptonPromptCheck = all_of(Lepton_genPartFlav.begin(), Lepton_genPartFlav.end(), [](int i){return i == 1;});  

  	return SameSignChargeCheck && (LeptonPromptCheck == 1);

  }};


  auto inv_mass{[](const floats& pts, const floats& etas, const floats& phis, const floats& ms)
   {

	std::cout << "print 36" << std::endl;

/*	std::cout << "pts.size() = " << pts.size() << std::endl;
	std::cout << "etas.size() = " << etas.size() << std::endl;
	std::cout << "phis.size() = " << phis.size() << std::endl;
	std::cout << "ms.size() = " << ms.size() << std::endl;
*/

    	//if (!all_equal(pts.size(), etas.size(), phis.size(), ms.size())){throw std::logic_error("Collections must be the same size");}
    	//else if(pts.empty()){throw std::logic_error("Collections must not be empty");}

    	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> vec{};

    	for (size_t i{0}; i < pts.size(); i++){
        	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> p{};
		p.SetPt(pts[i]);
        	p.SetEta(etas[i]);
		p.SetPhi(phis[i]);
		p.SetM(ms[i]);
        	vec += p;
    	}
	
    	return boost::numeric_cast<float>(vec.M());

  }};

  auto RecoZ{[](const float& LeadingleptonPt,    const float& LeadingleptonEta,    const float& LeadingleptonPhi,    const float& LeadingleptonMass,
		const float& SubleadingleptonPt, const float& SubleadingleptonEta, const float& SubleadingleptonPhi, const float& SubleadingleptonMass){

  	std::cout << "print 37" << std::endl;

  	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> ZBoson = {};
  	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> LeadingLepton = {};
  	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> SubleadingLepton = {};

  	LeadingLepton.SetPt(LeadingleptonPt);
	LeadingLepton.SetEta(LeadingleptonEta);
	LeadingLepton.SetPhi(LeadingleptonPhi);
	LeadingLepton.SetM(LeadingleptonMass);

  	SubleadingLepton.SetPt(SubleadingleptonPt);
	SubleadingLepton.SetEta(SubleadingleptonEta);
	SubleadingLepton.SetPhi(SubleadingleptonPhi);
	SubleadingLepton.SetM(SubleadingleptonMass);

  	ZBoson = LeadingLepton + SubleadingLepton;

  	return ZBoson;

  }};

  
  auto TLorentzVectorVariable{[](const int& VariableChoice, const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>& object){

  	std::cout << "print 38" << std::endl;

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

  auto TLorentzVectorVariablePt{[&TLorentzVectorVariable](const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>& object){
	std::cout << "print 158" << std::endl; 
	return TLorentzVectorVariable(1, object);
  }};
  
  auto TLorentzVectorVariablePhi{[&TLorentzVectorVariable](const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>& object){
	std::cout << "print 159" << std::endl; 
	return TLorentzVectorVariable(2, object);
  }};

  auto TLorentzVectorVariableEta{[&TLorentzVectorVariable](const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>& object){
	std::cout << "print 160" << std::endl;
	return TLorentzVectorVariable(3, object);
  }};
  
  auto TLorentzVectorVariableMass{[&TLorentzVectorVariable](const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>& object){
	std::cout << "print 161" << std::endl; 
	return TLorentzVectorVariable(4, object);
  }};


  auto deltaRcheck_float{[](const float& Object1_eta, const float& Object1_phi, const float& Object2_eta, const float& Object2_phi){

  	std::cout << "print 41" << std::endl;

  	float dR = sqrt(pow(Object1_eta - Object2_eta, 2) + pow(Object1_phi - Object2_phi, 2));
  	return dR;

  }};


  auto DeltaPhi_floatandfloat{[](const float& Object1_phi, const float& Object2_phi){

  	std::cout << "print 42" << std::endl;

  	double dPhi = abs(Object1_phi - Object2_phi);
  	return dPhi;

  }};

  auto LeptonFourMomentumFunction{[](const floats& Muon_pt, const floats& Muon_eta, const floats& Muon_phi, const floats& Muon_mass){
  
  	std::cout << "print 43" << std::endl;

  	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> Muon4Mo{};
  
  	for(long unsigned int i = 0; i < Muon_pt.size(); i++){

  		ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> vec{};

  		vec.SetPt(Muon_pt.at(i));
		vec.SetEta(Muon_eta.at(i));
		vec.SetPhi(Muon_phi.at(i));
		vec.SetM(Muon_mass.at(i));

  		Muon4Mo += vec;
  	}

  	return Muon4Mo;

  }};

  auto RochesterCorrections_testscript2{[&MCInt](const int& YearInteger, const int& MonteCarloInt, const ints& MuonCharge, const floats& MuonPt,
					         const floats& MuonEta, const floats& MuonPhi, const ints& Muon_genPartIdx, const ints& Muon_nTrackerLayers){

	std::cout << "print 44" << std::endl;


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


		switch(MCInt){

			case 1: 

				if(Muon_genPartIdx.size() > 0 && Muon_nTrackerLayers.size() > 0){

					if(mcSF > 0){
	
						RochCorrSF = rc.kSpreadMC(MuonCharge.at(i), MuonPt.at(i), MuonEta.at(i), MuonPhi.at(i), Muon_genPartIdx.at(0), s.at(i), m.at(i)); //(recommended), MC scale and resolution correction when matched gen muon is available
					}
					else{
						RochCorrSF = rc.kSmearMC(MuonCharge.at(i), MuonPt.at(i), MuonEta.at(i), MuonPhi.at(i), Muon_nTrackerLayers.at(0), u.at(i), s.at(i), m.at(i)); //MC scale and extra smearing when matched gen muon is not available

					}
				}
				else{RochCorrSF = 1.0;}

				break;
		
			default: RochCorrSF = rc.kScaleDT(MuonCharge.at(i), MuonPt.at(i), MuonEta.at(i), MuonPhi.at(i), s.at(i), m.at(i)); //data
				 break;

		}
	
		RochCorrVec.push_back(RochCorrSF);


	}

	
	return RochCorrVec;


  }};


  auto RochCorrVec_Function{[&MCInt, &YearInt, &RochesterCorrections_testscript2](const ints& MuonCharge,      const floats& MuonPt,           const floats& MuonEta, 
									          const floats& MuonPhi,       const ints& Muon_genPartIdx,    const ints& Muon_nTrackerLayers){

  	std::cout << "print 45" << std::endl;

	floats CorrectionFactor = RochesterCorrections_testscript2(YearInt, MCInt, MuonCharge, MuonPt, MuonEta, MuonPhi, Muon_genPartIdx, Muon_nTrackerLayers);

  	return CorrectionFactor;

  }};



  auto RochCorrMuon4Mo{[&ChannelInt](const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>& Muon4Mo, const floats& RochCorrVec){

  	std::cout << "print 47" << std::endl;

  	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> NewVec{};

	double NewVecMass = Muon4Mo.M() * RochCorrVec.at(0);
        double NewVecPt = Muon4Mo.Pt() * RochCorrVec.at(0);
  	double NewVecPhi = Muon4Mo.Phi() * RochCorrVec.at(0);
  	double NewVecEta = Muon4Mo.Eta() * RochCorrVec.at(0);

        switch(ChannelInt){
        
                case 2: NewVec.SetPt(NewVecPt); 
			NewVec.SetEta(NewVecEta);
			NewVec.SetPhi(NewVecPhi);
			NewVec.SetM(NewVecMass);

			break;

		default: NewVec.SetPt(Muon4Mo.Pt()); 
			 NewVec.SetEta(Muon4Mo.Eta());
			 NewVec.SetPhi(Muon4Mo.Phi());
			 NewVec.SetM(Muon4Mo.M());			

			 break;

	}

  	return NewVec;

  }};

  auto TLorentzVector_float{[](const int& VariableOption, const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>& object){
  
  	std::cout << "print 48" << std::endl;

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

  auto TLorentzVector_float_pt{[&TLorentzVector_float](const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>& object){
	return TLorentzVector_float(1, object);
  }};

  auto TLorentzVector_float_phi{[&TLorentzVector_float](const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>& object){
	return TLorentzVector_float(2, object);
  }};

  auto TLorentzVector_float_eta{[&TLorentzVector_float](const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>& object){
	return TLorentzVector_float(3, object);
  }};
  
  auto TLorentzVector_float_mass{[&TLorentzVector_float](const ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>& object){
	return TLorentzVector_float(4, object);
  }};

  auto z_mass_cut{[](const float& z_mass) {

  	std::cout << "print 49" << std::endl;

  	return !isinf(z_mass);//abs(z_mass - Z_MASS) < Z_MASS_CUT;

  }};


  auto RowReader2{[&FileNameJetSmear, &YearInt](const int& LineSpecified, const bool& sigmaJER, const bool& SF, const bool& up, 
					        const bool& down, const floats& Jet_eta, const floats& Jet_rho, const floats& Jet_pt) { 

  	std::cout << "print 50" << std::endl;

  	double Col1, Col2, Col3, Col4, Col5, Col6, Col7, Col8, Col9, Col10, Col11;
 
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

	std::cout << '\n' << std::endl;
	std::cout << '\n' << std::endl;
	std::cout << '\n' << std::endl;
	std::cout << "FileNameJetSmear = " << FileNameJetSmear << std::endl;
	std::cout << '\n' << std::endl;
        std::cout << '\n' << std::endl;
        std::cout << '\n' << std::endl;

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

  		if(sigmaJER == true && SF == false && up == false && down == false){

			if(  (Jet_eta.at(i) > abs(Col1) && Jet_eta.at(i) < abs(Col2)) &&
                             (Jet_rho.at(0) > abs(Col3) && Jet_rho.at(0) < abs(Col4)) &&
                             (Jet_pt.at(i) > abs(Col6) && Jet_pt.at(i) < abs(Col7) ) ){

        			float answer = sqrt( Col8*abs(Col8) / (Jet_pt.at(i)*Jet_pt.at(i))+Col9*Col9*pow(Jet_pt.at(i),Col11)+Col10*Col10 );
				AnswerVec.push_back(answer);

			}

		}
 	 	else if(sigmaJER == false && SF == true && up == false && down == false){
        		
			if(Jet_eta.at(i) > abs(Col1) && Jet_eta.at(i) < abs(Col2)){

        			AnswerVec.push_back(Col4);
  			
			}
		}
  		else if(sigmaJER == false && SF == false && up == true && down == false){

			if(Jet_eta.at(i) > abs(Col1) && Jet_eta.at(i) < abs(Col2)){

        			float UpValue = Col6;
        			AnswerVec.push_back(UpValue);

			}

		}
  		else if(sigmaJER == false && SF == false && up == false && down == true){
	
			if(Jet_eta.at(i) > abs(Col1) && Jet_eta.at(i) < abs(Col2)){

        			float DownValue = Col5;
        			AnswerVec.push_back(DownValue);

			}
		}
  		else{
			std::cout << "sigmaJER = " << sigmaJER << std::endl;
                        std::cout << "SF = " << SF << std::endl;
                        std::cout << "up = " << up << std::endl;
                        std::cout << "down = " << down << std::endl;

			throw std::logic_error("bools cannot be all true or all false"); 

		}

  
   	} //end of for loop


	std::cout << '\n' << std::endl;
        std::cout << "AnswerVec = " << AnswerVec << std::endl;
        std::cout << '\n' << std::endl;

	if(AnswerVec.size() == 0){
		std::cout << "Jet_eta = " << Jet_eta << std::endl; 
		std::cout << "Col1 = " << Col1 << std::endl;
		std::cout << "Col2 = " << Col2 << std::endl;
		std::cout << "Jet_rho = " << Jet_rho << std::endl;
		std::cout << "Jet_pt = " << Jet_pt << std::endl;
	}

   	return AnswerVec;


  }}; 



  auto linecounter{[&FileNameJetSmear, &YearInt](const bool& sigmaJER, const bool& SF, const bool& up, const bool& down){ 

  	std::cout << "print 51" << std::endl;

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

  	std::cout << "print 52" << std::endl;

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

	std::cout << "k = " << k << std::endl;

  	float factor;

  	for(long unsigned int i = 0; i < RowReader2(k, SigmaJER, JetSmearScaleFactor, Up, Down, Jet_eta, Jet_rho, Jet_pt).size(); i++){

		if(RowReader2(k, SigmaJER, JetSmearScaleFactor, Up, Down, Jet_eta, Jet_rho, Jet_pt).at(i) != 0){factor = RowReader2(k, SigmaJER, JetSmearScaleFactor, Up, Down, Jet_eta, Jet_rho, Jet_pt).at(i);}
		else{continue;}

  	}

  	return factor;

  }};

 auto sigma_JER{[&RowReader3](const floats& Jet_eta, const floats& Jet_rho,const floats& Jet_pt){

  std::cout << "print 53" << std::endl;

  bool SigmaJER = true;
  bool JetSmearScaleFactor = false;
  bool Up = false;
  bool Down = false;

  std::cout << '\n' << std::endl;
  std::cout << '\n' << std::endl;
  std::cout << "RowReader3(SigmaJER, JetSmearScaleFactor, Up, Down, Jet_eta, Jet_rho, Jet_pt) = " << RowReader3(SigmaJER, JetSmearScaleFactor, Up, Down, Jet_eta, Jet_rho, Jet_pt) << std::endl;
  std::cout << '\n' << std::endl;
  std::cout << '\n' << std::endl;

  return RowReader3(SigmaJER, JetSmearScaleFactor, Up, Down, Jet_eta, Jet_rho, Jet_pt);

}};


  auto SJER_Nominal_Function{[&RowReader3](const floats& Jet_eta, const floats& Jet_rho, const floats& Jet_pt){

  	std::cout << "print 53" << std::endl;

  	bool SigmaJER = false;
  	bool JetSmearScaleFactor = true;
  	bool Up = false;
  	bool Down = false;
 
  	return RowReader3(SigmaJER, JetSmearScaleFactor, Up, Down, Jet_eta, Jet_rho, Jet_pt);

  }};

  auto SJER_Up_Function{[&RowReader3](const floats& Jet_eta, const floats& Jet_rho, const floats& Jet_pt){

  	std::cout << "print 54" << std::endl;

  	bool SigmaJER = false;
  	bool JetSmearScaleFactor = false;
  	bool Up = true;
  	bool Down = false;
  
  	return RowReader3(SigmaJER, JetSmearScaleFactor, Up, Down, Jet_eta, Jet_rho, Jet_pt);

  }};

  auto SJER_Down_Function{[&RowReader3](const floats& Jet_eta, const floats& Jet_rho, const floats& Jet_pt){

  	std::cout << "print 55" << std::endl;

  	bool SigmaJER = false;
  	bool JetSmearScaleFactor = false;
  	bool Up = false;
  	bool Down = true;
  
  	return RowReader3(SigmaJER, JetSmearScaleFactor, Up, Down, Jet_eta, Jet_rho, Jet_pt);

  }};


  auto MaxComparison{[](const float& sJER_nominal){

  	std::cout << "print 56" << std::endl;

 	float MaximumFloats = sqrt(sJER_nominal*sJER_nominal - 1);

 	if(MaximumFloats > 0){return MaximumFloats;}
 	else{float zero = 0.0; return zero;}

  }};


  auto JetSmearingFunction_HybridMethod{[&MaxComparison, &MCInt](const floats& pT, const floats& eta, const floats& phi, const floats& pT_ptcl, const floats& eta_ptcl, 
							         const floats& phi_ptcl, const float& sJER_nominal, const float& sigma_JER_input, const ints& Jet_genJetIdx){

  	std::cout << "print 57" << std::endl;

  	floats cJER_vec{};

  	for(long unsigned int i = 0; i < pT.size(); i++){

		switch(MCInt){

			case 0: cJER_vec.push_back(1.0); break;
			default:

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

				break;

			}

  	}

  	return cJER_vec;
  
  }};



  auto ApplyCJER{[](const floats& JetPt, const floats& JetEta, const floats& JetPhi, const floats& JetMass, const floats& cJER, const unsigned int& nJet){

  	std::cout << "print 58" << std::endl;

  	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> OutputVec{};

  	for(unsigned int i = 0; i < nJet; i++){

    		ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> JetFourMomentum_New{};
    		float JetPt_new = JetPt.at(i) * cJER.at(0);
    		float JetEta_new = JetEta.at(i) * cJER.at(0);
    		float JetPhi_new = JetPhi.at(i) * cJER.at(0);
    		float JetMass_new = JetMass.at(i) * cJER.at(0);
		
    		JetFourMomentum_New.SetPt(JetPt_new);
		JetFourMomentum_New.SetEta(JetEta_new);
		JetFourMomentum_New.SetPhi(JetPhi_new);
		JetFourMomentum_New.SetM(JetMass_new);
    		OutputVec += JetFourMomentum_New;

  	}

 	return OutputVec; 

  }};


 auto GetSmearedJetPt{[](const floats& JetPt, const floats& JetEta, const floats& JetPhi, const floats& JetMass, const floats& cJER, const unsigned int& nJet){

	std::cout << "print 59" << std::endl;

	floats new_pts{};

        for(unsigned int i = 0; i < nJet; i++){

                float JetPt_new = JetPt.at(i) * cJER.at(0);
                new_pts.push_back(JetPt_new);

        }

        return new_pts;

  }};


  auto GetSmearedJetPhi{[](const floats& JetPt, const floats& JetEta, const floats& JetPhi, const floats& JetMass, const floats& cJER, const unsigned int& nJet){
        
        std::cout << "print 60" << std::endl;

        floats new_phis{};

        for(unsigned int i = 0; i < nJet; i++){
  
                float JetPhi_new = JetPhi.at(i) * cJER.at(0); 
                new_phis.push_back(JetPhi_new);

        }

        return new_phis;

  }};

  
  auto GetSmearedJetEta{[](const floats& JetPt, const floats& JetEta, const floats& JetPhi, const floats& JetMass, const floats& cJER, const unsigned int& nJet){
 
        std::cout << "print 61" << std::endl;

        floats new_etas{};

        for(unsigned int i = 0; i < nJet; i++){
        
                float JetEta_new = JetEta.at(i) * cJER.at(0);
                new_etas.push_back(JetEta_new);
  
        }
        
        return new_etas;
        
  }}; 


  auto GetSmearedJetMass{[](const floats& JetPt, const floats& JetEta, const floats& JetPhi, const floats& JetMass, const floats& cJER, const unsigned int& nJet){
 
        std::cout << "print 62" << std::endl;

        floats new_mass{};

        for(unsigned int i = 0; i < nJet; i++){
        
                float JetMass_new = JetMass.at(i) * cJER.at(0);
                new_mass.push_back(JetMass_new);
  
        }
        
        return new_mass;
        
  }}; 

  auto SumSquared2LeadingJets_pT{[](const float& LeadingJetPt, const float& SubleadingJetPt){

  	std::cout << "print 63" << std::endl;

  	double SumSquaredPt = pow(LeadingJetPt + SubleadingJetPt, 2);
  	return SumSquaredPt;

  }};


  auto JetSum{[](const float& LeadingJet, const float& SubleadingJet, const float& ThirdJet, const float& FourthJet){

  	std::cout << "print 64" << std::endl;

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

  	std::cout << "print 65" << std::endl;

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
 
  	std::cout << "print 66" << std::endl;

  	float HTOutput = abs(Pt);
  	return HTOutput;

  }};

  auto TotJetHT{[](const float& LeadingJetHT, const float& SubleadingJetHT, const float& ThirdJetHT, const float& FourthJetHT){
  
  	std::cout << "print 67" << std::endl;

  	float TotJetHTOutput = LeadingJetHT + SubleadingJetHT + ThirdJetHT + FourthJetHT;
  	return TotJetHTOutput;

  }};

  auto TotLepHT{[](const float& LeadingLeptonHT, const float& SubleadingLeptonHT){

  	std::cout << "print 68" << std::endl;

  	float TotLepHTOutput = LeadingLeptonHT + SubleadingLeptonHT;
  	return TotLepHTOutput;

  }};

  
  auto TotHTOverTotpT{[](const float& TotHT, const float& TotpT){

  	std::cout << "print 69" << std::endl;

  	float TotHTOverTotpTOutput = TotHT / TotpT;
  	return TotHTOverTotpTOutput;

  }};


  auto LepSum{[](const float& LeadingLep, const float& SubleadingLep){

  	std::cout << "print 70" << std::endl;

  	float LepSumOutput = LeadingLep + SubleadingLep;
  	return LepSumOutput;

  }};


  auto InvMass_AllJets{[](const float& LeadingJetPt,   const float& SubleadingJetPt,   const float& ThirdJetPt,   const float& FourthJetPt,
			  const float& LeadingJetEta,  const float& SubleadingJetEta,  const float& ThirdJetEta,  const float& FourthJetEta,
		    	  const float& LeadingJetPhi,  const float& SubleadingJetPhi,  const float& ThirdJetPhi,  const float& FourthJetPhi,
			  const float& LeadingJetMass, const float& SubleadingJetMass, const float& ThirdJetMass, const float& FourthJetMass){

  	std::cout << "print 71" << std::endl;

  	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> Jet1 = {};
  	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> Jet2 = {};
 	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> Jet3 = {};
  	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> Jet4 = {};

  	Jet1.SetPt(LeadingJetPt);
	Jet1.SetEta(LeadingJetEta);
	Jet1.SetPhi(LeadingJetPhi);
	Jet1.SetM(LeadingJetMass);

  	Jet2.SetPt(SubleadingJetPt);
	Jet2.SetEta(SubleadingJetEta);
	Jet2.SetPhi(SubleadingJetPhi);
	Jet2.SetM(SubleadingJetMass);

	Jet3.SetPt(ThirdJetPt);
        Jet3.SetEta(ThirdJetEta);
        Jet3.SetPhi(ThirdJetPhi);
        Jet3.SetM(ThirdJetMass);

	Jet4.SetPt(FourthJetPt);
        Jet4.SetEta(FourthJetEta);
        Jet4.SetPhi(FourthJetPhi);
        Jet4.SetM(FourthJetMass);
	
  	float InvMassAllJets = (Jet1 + Jet2 + Jet3 + Jet4).M();
  	return InvMassAllJets;

  }};


  auto InvMass_3Jets{[](const float& LeadingJetPt,   const float& SubleadingJetPt,   const float& ThirdJetPt,
			const float& LeadingJetEta,  const float& SubleadingJetEta,  const float& ThirdJetEta,
			const float& LeadingJetPhi,  const float& SubleadingJetPhi,  const float& ThirdJetPhi,
			const float& LeadingJetMass, const float& SubleadingJetMass, const float& ThirdJetMass){
  
  	std::cout << "print 72" << std::endl;

  	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> Jet1 = {};
  	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> Jet2 = {};
  	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> Jet3 = {};

  	Jet1.SetPt(LeadingJetPt);
	Jet1.SetEta(LeadingJetEta);
	Jet1.SetPhi(LeadingJetPhi);
	Jet1.SetM(LeadingJetMass); 

	Jet2.SetPt(SubleadingJetPt);
        Jet2.SetEta(SubleadingJetEta);
        Jet2.SetPhi(SubleadingJetPhi);
        Jet2.SetM(SubleadingJetMass);

	Jet3.SetPt(ThirdJetPt);
        Jet3.SetEta(ThirdJetEta);
        Jet3.SetPhi(ThirdJetPhi);
        Jet3.SetM(ThirdJetMass);


  	float InvMass3Jets = (Jet1 + Jet2 + Jet3).M();
  	return InvMass3Jets;

  }};


  auto tight_jets_function{[&YearInt](const floats& Jet_pt_Selection, const floats& Jet_eta_Selection, const ints& Jet_jetId_Selection, const floats& dRJet_lep){

  	std::cout << "print 73" << std::endl;

	std::cout << "Jet_pt_Selection.size() = " << Jet_pt_Selection.size() << std::endl;
	std::cout << "Jet_eta_Selection.size() = " << Jet_eta_Selection.size() << std::endl;
	std::cout << "Jet_jetId_Selection.size() = " << Jet_jetId_Selection.size() << std::endl;
	std::cout << "dRJet_lep.size() = " << dRJet_lep.size() << std::endl;

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

  	std::cout << "print 74" << std::endl;

  	auto njet{count_if(tight_jets.begin(), tight_jets.end(), [](int i) { return i; })};
  	return njet >= 4 && njet <= 6;

  }};

  auto bjet_id{[](const ints& tight_jets, const floats& btags, const floats& etas) {
     
        std::cout << "print 75" << std::endl;
	return /*tight_jets &&*/ (btags > 0.8838f) && (etas < MaxTrackerEta);
  
  }};


  auto numberofbjets{[](const ints& bjets) {

	std::cout << "print 76" << std::endl;
        const auto nbjet{std::count_if(bjets.begin(), bjets.end(), [](int i) { return i; })};
        return nbjet;

  }};

  auto BTAGEFF_bjet_id_WP{[](const ints& tight_jets, const floats& btags, const floats& etas, const ints& Jet_hadronFlavour) {

	std::cout << "print 77" << std::endl;
	return abs(Jet_hadronFlavour) == 0 && btags > 0.8838f && abs(etas) < MaxTrackerEta;
	
  }};

  auto BTAGEFF_nonbjet_id_WP{[](const ints& tight_jets, const floats& btags, const floats& etas, const ints& Jet_hadronFlavour){

  	std::cout << "print 81" << std::endl;
    	return (abs(Jet_hadronFlavour) == 1 || abs(Jet_hadronFlavour) == 2) && btags > 0 && abs(etas) < MaxTrackerEta;

  }};


  auto BTAGEFF_bjet_id{[](const ints& tight_jets, const floats& etas, const ints& Jet_hadronFlavour) {

  	std::cout << "print 82" << std::endl;
	return abs(Jet_hadronFlavour) == 0 && abs(etas) < MaxTrackerEta;

  }};


 auto BTAGEFF_nonbjet_id{[](const ints& tight_jets, const floats& etas, const ints& Jet_hadronFlavour){
	
	std::cout << "print 86" << std::endl;
	return (abs(Jet_hadronFlavour) == 1 || abs(Jet_hadronFlavour) == 2) && abs(etas) < MaxTrackerEta;

  }};

  auto bjet_cut{[](const ints& bjets) {

        std::cout << "print 87" << std::endl;

        const auto nbjet{std::count_if(bjets.begin(), bjets.end(), [](int i) { return i; })};
        return nbjet >= 1 && nbjet <= 2;

  }};

  auto find_lead_mask{[](const ints& mask, const floats& vals) {
  
  	std::cout << "print 88" << std::endl;

  	const auto masked_vals{mask * vals};
  	const auto max_idx{boost::numeric_cast<size_t>(std::distance(masked_vals.begin(), max_element(masked_vals.begin(), masked_vals.end())))};
  	ints lead_mask(masked_vals.size(), 0); 
  	lead_mask.at(max_idx) = 1;
  	return lead_mask;

  }};

  auto WPair{[](const floats& pts, const floats& etas, const floats& phis, const floats& ms, const ints& tight_jets, const ints& lead_bjet) {

        std::cout << "print 89" << std::endl;

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

                        //auto jet1{TLorentzVector{}};
                        //auto jet2{TLorentzVector{}};
                        //jet1.SetPtEtaPhiM(pts.at(i), etas.at(i), phis.at(i), ms.at(i));
                        //jet2.SetPtEtaPhiM(pts.at(j), etas.at(j), phis.at(j), ms.at(j));

			ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> jet1{};
			ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> jet2{};
	
			jet1.SetPt(pts.at(i));
			jet1.SetEta(etas.at(i));
			jet1.SetPhi(phis.at(i));
			jet1.SetM(ms.at(i));

			jet2.SetPt(pts.at(j));
                        jet2.SetEta(etas.at(j));
                        jet2.SetPhi(phis.at(j));
                        jet2.SetM(ms.at(j));


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

        auto jet1{ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>{}};
        auto jet2{ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>{}};

        jet1.SetPt(pts.at(jet_index_1));
	jet1.SetEta(etas.at(jet_index_1));
	jet1.SetPhi(phis.at(jet_index_1));
	jet1.SetM(ms.at(jet_index_1));

        jet2.SetPt(pts.at(jet_index_2));
	jet2.SetEta(etas.at(jet_index_2));
	jet2.SetPhi(phis.at(jet_index_2));
	jet2.SetM(ms.at(jet_index_2));

	return w_pair;

  }};


  auto WPairJet1{[](const floats& pts, const floats& etas, const floats& phis, const floats& ms, const ints& tight_jets, const ints& lead_bjet) {

	std::cout << "print 90" << std::endl;

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

                        auto jet1{ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>{}};
                        auto jet2{ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>{}};

                        jet1.SetPt(pts.at(i));
			jet1.SetEta(etas.at(i));
			jet1.SetPhi(phis.at(i));
			jet1.SetM(ms.at(i));

                        jet2.SetPt(pts.at(j));
			jet2.SetEta(etas.at(j));
			jet2.SetPhi(phis.at(j));
			jet2.SetM(ms.at(j));

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

        auto jet1{ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>{}};
        auto jet2{ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>{}};

        jet1.SetPt(pts.at(jet_index_1));
	jet1.SetEta(etas.at(jet_index_1));
	jet1.SetPhi(phis.at(jet_index_1));
	jet1.SetM(ms.at(jet_index_1));
     
        jet2.SetPt(pts.at(jet_index_2));
	jet2.SetEta(etas.at(jet_index_2));
	jet2.SetPhi(phis.at(jet_index_2));
	jet2.SetM(ms.at(jet_index_2));

        return jet1;
 
  }};


  auto WPairJet2{[](const floats& pts, const floats& etas, const floats& phis, const floats& ms, const ints& tight_jets, const ints& lead_bjet) {

	std::cout << "print 91" << std::endl;

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

                        auto jet1{ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>{}};
                        auto jet2{ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>{}};
              
	                jet1.SetPt(pts.at(i));
			jet1.SetEta(etas.at(i));
			jet1.SetPhi(phis.at(i));
			jet1.SetM(ms.at(i));

                        jet2.SetPt(pts.at(j));
			jet2.SetEta(etas.at(j));
			jet2.SetPhi(phis.at(j));
			jet2.SetM(ms.at(j));

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

        auto jet1{ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>{}};
        auto jet2{ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>{}};

        jet1.SetPt(pts.at(jet_index_1));
  	jet1.SetEta(etas.at(jet_index_1));
	jet1.SetPhi(phis.at(jet_index_1));
	jet1.SetM(ms.at(jet_index_1));

        jet2.SetPt(pts.at(jet_index_2));
	jet2.SetEta(etas.at(jet_index_2));
	jet2.SetPhi(phis.at(jet_index_2));
	jet2.SetM(ms.at(jet_index_2));

        return jet2;   
 
  }};


  auto deltaRcheck_W_function{[](const doubles& Object1_phi_Selection, const doubles& Object1_eta_Selection,
				 const doubles& Object2_eta_Selection, const doubles& Object2_phi_Selection){

  	std::cout << "print 93" << std::endl;

  	doubles dR = sqrt(pow(Object1_eta_Selection - Object2_eta_Selection, 2) + pow(Object1_phi_Selection - Object2_phi_Selection, 2));
  	return dR;

  }};


  auto deltaRcheck_DoublesAndFloats{[](const doubles& Object1_phi_Selection, const doubles& Object1_eta_Selection,
                                       const floats& Object2_eta_Selection, const floats& Object2_phi_Selection){

        std::cout << "print 93 part 2" << std::endl;

        std::cout << "Object1_phi_Selection.size() = " << Object1_phi_Selection.size() << std::endl;
	std::cout << "Object1_phi_Selection.size() = " << Object1_phi_Selection.size() << std::endl;
        std::cout << "Object2_phi_Selection.size() = " << Object2_phi_Selection.size() << std::endl;
        std::cout << "Object2_phi_Selection.size() = " << Object2_phi_Selection.size() << std::endl;

	doubles dR;

	if( Object2_phi_Selection.size() > 0 ){
        	dR = sqrt(pow(Object1_eta_Selection - Object2_eta_Selection, 2) + pow(Object1_phi_Selection - Object2_phi_Selection, 2));
	}
	else{dR = sqrt(pow(Object1_eta_Selection, 2) + pow(Object1_phi_Selection, 2));}

        return dR;

  }}; 

  auto DeltaPhi_function2{[](const doubles& Object1_phi_Selection, const doubles& Object2_phi_Selection){

  	std::cout << "print 94" << std::endl;

  	doubles dPhi = abs(Object1_phi_Selection - Object2_phi_Selection);
  	return dPhi;

  }};


  auto deltaRcheck_W_function2{[](const doubles& Object1_phi_Selection, const doubles& Object1_eta_Selection,
				  const float& Object2_eta_Selection,   const float& Object2_phi_Selection){
 
  	std::cout << "print 95" << std::endl;

  	doubles dR = sqrt(pow(Object1_eta_Selection - Object2_eta_Selection, 2) + pow(Object1_phi_Selection - Object2_phi_Selection, 2));
  	return dR;

  }};

  auto DeltaPhi_doublesandfloat{[](const doubles& Object1_phi, const float& Object2_phi){

  	std::cout << "print 96" << std::endl;

  	doubles dPhi = abs(Object1_phi - Object2_phi);
  	return dPhi;

  }};
 

  auto DeltaPhi_doublesandfloats{[](const doubles& Object1_phi, const floats& Object2_phi){

        std::cout << "print 96" << std::endl;

	std::cout << "Object1_phi.size() = " << Object1_phi.size() << std::endl;
	std::cout << "Object2_phi.size() = " << Object2_phi.size() << std::endl;

	doubles dPhi;

	if( Object2_phi.size() > 0 ){
        	dPhi = abs(Object1_phi - Object2_phi);
	}
	else{dPhi = abs(Object1_phi);}

        return dPhi;

  }};
 
 
  auto HT_double{[](const doubles& Pt){

  	std::cout << "print 97" << std::endl;

  	doubles HT_Output = abs(Pt);
  	return HT_Output;

  }};

  auto RecoWHT{[](const floats& RecoWPt){

  	std::cout << "print 98" << std::endl;

  	floats RecoWHTOutput = abs(RecoWPt);
  	return RecoWHTOutput;

  }};

  auto TransverseWMass{[](const double& dPhi_j1j2, const doubles& WPairJet1Pt, const doubles& WPairJet2Pt){

  	std::cout << "print 99" << std::endl;

  	doubles mtW = sqrt(2 * WPairJet1Pt * WPairJet2Pt * (1 - cos(dPhi_j1j2)) );
  	return mtW;

  }};


  auto w_mass_cut{[&ZPlusJetsCRInt](const float& w_mass, const float& MET_sumEt) {
	
  	std::cout << "print 100" << std::endl;
	
	switch(ZPlusJetsCRInt){
      
    		case 0: return !isinf(w_mass); break; //( abs(w_mass - W_MASS) < W_MASS_CUT ); break;
    		case 1: return !isinf(w_mass) && (MET_sumEt < 50); break;
      
        }

  }};


  auto WLorentzVector{[](const floats& w_pair_pt, const floats& w_pair_eta, const floats& w_pair_phi, const float& w_mass, const ints& w_reco_jets){

  	std::cout << "print 101" << std::endl;

  	const auto nRecoWBosons{std::count_if(w_reco_jets.begin(), w_reco_jets.end(), [](int i) { return i; })};
  	auto RecoW = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>{};
  
  	for(int i = 0; i < nRecoWBosons; i++){

	  	auto Vec = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>{};

	  	Vec.SetPt(w_pair_pt.at(i));
		Vec.SetEta(w_pair_eta.at(i));
		Vec.SetPhi(w_pair_phi.at(i));
		Vec.SetM(w_mass); 

	 	RecoW += Vec;

  	}

  	return RecoW;

  }};

  auto NumberOfSmearedTightJetsFunction{[](const ints& tight_jets){

	ints ntightjets_vec;
	const auto ntightjets{std::count_if(tight_jets.begin(), tight_jets.end(), [](int i) { return i; })};

	ntightjets_vec.push_back(ntightjets);
        return ntightjets_vec;

  }};


  auto bjet_variable{[](const floats& Jet_variable, const ints& nJet, const ints& lead_bjet){

  	std::cout << "print 102" << std::endl;

  	floats vec{};

  	for(int i = 0; i < nJet.size(); i++){
        	if(lead_bjet.at(i) == 1){ 
			vec.push_back(Jet_variable.at(i));
		}

  	}

  	return vec;

  }};


  auto BLorentzVector{[](const floats& bjet_pt, const floats& bjet_eta, const floats& bjet_phi, const floats& bjet_mass){

  	std::cout << "print 103" << std::endl;

  	auto BJets = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>{};

  	for(long unsigned int i = 0; i < bjet_pt.size(); i++){

		auto Vec = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>{};

		Vec.SetPt(bjet_pt.at(i));
		Vec.SetEta(bjet_eta.at(i));
		Vec.SetPhi(bjet_phi.at(i));
		Vec.SetM(bjet_mass.at(i));

		BJets += Vec;

  	}

  	return BJets;

  }};


  auto top_reconstruction_function{[](const floats& bjets_pt, const floats& bjets_eta, const floats& bjets_phi, const floats& bjets_mass,
				      const floats& w_pair_pt, const floats& w_pair_eta, const floats& w_pair_phi, const float& w_mass ){

  	std::cout << "print 104" << std::endl;

	float leadingbjetpt; float leadingbjeteta; float leadingbjetphi; float leadingbjetmass;

	if(bjets_pt.size() > 1){
		
		leadingbjetpt = (bjets_pt.at(0) > bjets_pt.at(1)) ? bjets_pt.at(0) : bjets_pt.at(1);
		leadingbjeteta = (bjets_eta.at(0) > bjets_eta.at(1)) ? bjets_eta.at(0) : bjets_eta.at(1);
		leadingbjetphi = (bjets_phi.at(0) > bjets_phi.at(1)) ? bjets_phi.at(0) : bjets_phi.at(1);
		leadingbjetmass = (bjets_mass.at(0) > bjets_mass.at(1)) ? bjets_mass.at(0) : bjets_mass.at(1);
	}
	else{leadingbjetpt = bjets_pt.at(0); leadingbjeteta = bjets_eta.at(0); leadingbjetphi = bjets_phi.at(0); leadingbjetmass = bjets_mass.at(0);}

  	auto reco_top = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>{}; 
  	auto BJets = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>{};
  	auto RecoW = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>{};

  	double top_reco_mass = std::numeric_limits<double>::infinity();
  	size_t index_1{std::numeric_limits<size_t>::max()};
  	const size_t num{w_pair_pt.size()};

  	for(unsigned int i = 0; i < num; i++){

  		BJets.SetPt(leadingbjetpt);
  		BJets.SetEta(leadingbjeteta);
		BJets.SetPhi(leadingbjetphi);
		BJets.SetM(leadingbjetmass);

		RecoW.SetPt(w_pair_pt.at(i));
		RecoW.SetEta(w_pair_eta.at(i));
		RecoW.SetPhi(w_pair_phi.at(i));
		RecoW.SetM(w_mass);		

  		const double reco_mass = (RecoW + BJets).M(); 

  		if(abs(TOP_MASS - reco_mass) < abs(TOP_MASS - top_reco_mass)){

	  		top_reco_mass = reco_mass;
	  		index_1 = i;

  		}

  	}


  	BJets.SetPt(leadingbjetpt);
	BJets.SetEta(leadingbjeteta);
	BJets.SetPhi(leadingbjetphi);
	BJets.SetM(leadingbjetmass);

  	RecoW.SetPt(w_pair_pt.at(index_1));
	RecoW.SetEta(w_pair_eta.at(index_1));
	RecoW.SetPhi(w_pair_phi.at(index_1));
	RecoW.SetM(w_mass);

  	reco_top = RecoW + BJets;	


  	return reco_top;

  }};

  
  auto deltaRcheck_Top_function{[](const doubles& Object1_phi_Selection, const doubles& Object1_eta_Selection,
				   const float& Object2_eta_Selection,   const float& Object2_phi_Selection){

  	std::cout << "print 106" << std::endl;

  	doubles dR = sqrt(pow(Object1_eta_Selection - Object2_eta_Selection, 2) + pow(Object1_phi_Selection - Object2_phi_Selection, 2));
  	return dR;

  }};

  auto deltaRcheck_WTop_function{[](const floats& Object1_phi_Selection,  const floats& Object1_eta_Selection,
				    const doubles& Object2_eta_Selection, const doubles& Object2_phi_Selection){

  	std::cout << "print 107" << std::endl;

  	doubles dR_vec{};

  	for(long unsigned int i = 0; i < Object1_phi_Selection.size(); i++){

  		double dR = sqrt(pow(Object1_eta_Selection.at(i) - Object2_eta_Selection.at(0), 2) + pow(Object1_phi_Selection.at(i) - Object2_phi_Selection.at(0), 2));
  		dR_vec.push_back(dR);

  	}

  	return dR_vec;

  }};

  
  auto MinDeltaR{[](const ints& nJet, const doubles& RecoZPhi, const doubles& RecoZEta, const floats& Jet_Phi_Selection, const floats& Jet_eta_Selection){

  	std::cout << "print 108" << std::endl;

    	doubles output_vec;
	double Output;  

    	for(int i = 0; i < nJet.size(); i++){

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


  auto MaxDeltaR{[](const ints& nJet, const doubles& RecoZPhi, const doubles& RecoZEta, const floats& Jet_Phi_Selection, const floats& Jet_eta_Selection){

        std::cout << "print 108" << std::endl;

        doubles output_vec;
        double Output;

        for(int i = 0; i < nJet.size(); i++){

                if(RecoZEta.size() > 1){
                        double DeltaR = sqrt(pow(RecoZPhi.at(i) - Jet_Phi_Selection.at(i), 2) + pow(RecoZEta.at(i) - Jet_eta_Selection.at(i), 2));
                        double DeltaR2 = sqrt(pow(RecoZPhi.at(i+1) - Jet_Phi_Selection.at(i+1), 2) + pow(RecoZEta.at(i+1) - Jet_eta_Selection.at(i+1), 2));

                        Output = (DeltaR2 > DeltaR) ? DeltaR2 : DeltaR;

                }
                else{Output = sqrt(pow(RecoZPhi.at(0) - Jet_Phi_Selection.at(i), 2) + pow(RecoZEta.at(0) - Jet_eta_Selection.at(i), 2));}

                output_vec.push_back(Output);

        }

        return output_vec;

  }}; 


  auto MinDeltaPhi{[](const ints& nJet, const doubles& RecoZPhi, const floats& Jet_Phi_Selection){

  	std::cout << "print 109" << std::endl;

  	double output;
  	doubles output_vec{};

  	for(int i = 0; i < nJet.size(); i++){

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


  auto MaxDeltaPhi{[](const ints& nJet, const doubles& RecoZPhi, const floats& Jet_Phi_Selection){
        
        std::cout << "print 109" << std::endl;
        
        double output;
        doubles output_vec{};

        for(int i = 0; i < nJet.size(); i++){

                if(RecoZPhi.size() > 1){

                        double DeltaPhi = std::abs(RecoZPhi.at(i) - Jet_Phi_Selection.at(i));
                        double DeltaPhi2 = std::abs(RecoZPhi.at(i+1) - Jet_Phi_Selection.at(i+1));

                        output = (DeltaPhi2 > DeltaPhi) ? DeltaPhi2 : DeltaPhi;

                }
                else{output = std::abs(RecoZPhi.at(0) - Jet_Phi_Selection.at(i));}

                output_vec.push_back(output);

        }

        return output_vec;

  }};

  auto dR_Lepton_LeadingBJet_Function{[](const floats& bjeteta, const float& LeptonEta, const floats& bjetphi, const float& LeptonPhi){

  	std::cout << "print 110" << std::endl;

  	doubles DeltaR = sqrt(pow(LeptonPhi - bjetphi, 2) + pow(LeptonEta - bjeteta, 2));
  	return DeltaR;

  }};


  auto DeltaPhi_Lepton_BJet{[](const floats& Jet_phi_Selection, const float& LeptonPhi){

  	std::cout << "print 111" << std::endl;

  	doubles DeltaPhi = abs(LeptonPhi - Jet_phi_Selection);
  	return DeltaPhi;

  }};


  auto MET_function{[](const floats& MET_input){

  	std::cout << "print 112" << std::endl;
  	return MET_input;
  
  }};

  
  auto BJetOutputDiscriminantFunction{[](const float& JetPt, const floats& Jet_btagCSVV2, const ints& tight_jets, const floats& Jet_eta_Selection){

  	std::cout << "print 112" << std::endl;
  	return JetPt && (Jet_btagCSVV2  > 0.8838) /*&& tight_jets */&& (abs(Jet_eta_Selection) < MaxTrackerEta);

  }};


  auto DeltaPhi_function4{[](const floats& Object1_phi, const doubles& Object2_phi){

  	std::cout << "print 113" << std::endl;

 	doubles dPhi_vec{};

 	for(long unsigned int i = 0; i < Object1_phi.size(); i++){

 		double dPhi = Object1_phi.at(i) - Object2_phi.at(0);
		dPhi_vec.push_back(dPhi);

 	}

 	return dPhi_vec;

  }};


  auto TotalVariable_System{[](const doubles& RecoZInput, const floats& RecoWInput, const doubles& TopInput, const float& TotLepInput, const float& TotJetInput){

  	std::cout << "print 114" << std::endl;

  	doubles TotalSystemOutput = RecoZInput + RecoWInput.at(0) + TopInput + TotLepInput + TotJetInput;
  	return TotalSystemOutput;

  }};


  auto inv_mass_doubles{[](const doubles& pts, const doubles& etas, const doubles& phis, const doubles& ms){

	std::cout << "print 115" << std::endl;
/*
	std::cout << "pts.size() = " << pts.size() << std::endl;
        std::cout << "etas.size() = " << etas.size() << std::endl;
        std::cout << "phis.size() = " << phis.size() << std::endl;
        std::cout << "ms.size() = " << ms.size() << std::endl;
*/
    	if (!all_equal(pts.size(), etas.size(), phis.size(), ms.size()))
    	{
        	throw std::logic_error("Collections must be the same size");
    	}
    	else if (pts.empty())
   	{
        	throw std::logic_error("Collections must not be empty");
    	}

    	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> vec{};
    	for (size_t i{0}; i < pts.size(); i++)
    	{
        	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> p{};

        	p.SetPt(pts[i]);
		p.SetEta(etas[i]);
		p.SetPhi(phis[i]);
		p.SetM(ms[i]);

        	vec += p;
    	}

    	return boost::numeric_cast<float>(vec.M());

  }};

  auto UnweightedTopPt{[](const doubles& pts){

	std::cout << "print 116" << std::endl;
        return pts;

  }};


  auto TopReweighting_topquark{[&MCInt](const ints& GenPart_pdgId, const ints& GenPart_statusFlags, const floats& GenPart_pt){

  	std::cout << "print 117" << std::endl;
	std::cout << "GenPart_pdgId.size() = " << GenPart_pdgId.size() << std::endl;
	std::cout << "GenPart_statusFlags.size() = " << GenPart_statusFlags.size() << std::endl;
	std::cout << "GenPart_pt.size() = " << GenPart_pt.size() << std::endl;

	switch(MCInt){

		case 1:	return GenPart_pdgId == 6 && GenPart_statusFlags == 13 && GenPart_pt > 0; break;
		default: ints out_vec_top(GenPart_pdgId.size(), 1.0); return out_vec_top; break; //putting anything here since a weight of 1 is applied for data

	} 

  }};

  auto TopReweighting_antitopquark{[&MCInt](const ints& GenPart_pdgId, const ints& GenPart_statusFlags, const floats& GenPart_pt){
		
	std::cout << "print 118" << std::endl;

	switch(MCInt){

		case 1: return GenPart_pdgId == -6 && GenPart_statusFlags == 13 && GenPart_pt > 0; break;
		default: ints out_vec_antitop(GenPart_pdgId.size(), 1.0); return out_vec_antitop; break; //putting anything here since a weight of 1 is applied for data

	} 

  }};


  auto TopReweighting_weight{[&ProcessInt](const ints& TopReweighting_topquark_input, const ints& TopReweighting_antitopquark_input){

	std::cout << "print 119" << std::endl;

	doubles SF_top = exp(-0.0615-(0.00005* TopReweighting_topquark_input) );
	doubles SF_antitop = exp(-0.0615-(0.00005* TopReweighting_antitopquark_input) );

	doubles weight;

	switch(ProcessInt){
		case 28: weight = sqrt( SF_top * SF_antitop); break; //ttbar_2l2nu
		case 29: weight = sqrt( SF_top * SF_antitop); break; //ttbar_madgraph
		case 30: weight = sqrt( SF_top * SF_antitop); break; //ttbar_madgraph_ext
		case 31: weight = sqrt( SF_top * SF_antitop); break; //ttbar_TTToHadronic
		case 32: weight = sqrt( SF_top * SF_antitop); break; //ttbar_TTToHadronic_ext
		case 33: weight = sqrt( SF_top * SF_antitop); break; //ttbar_TTToSemileptonic
		case 34: weight = sqrt( SF_top * SF_antitop); break; //ttbar_TTToSemileptonic_ext
		case 35: weight = sqrt( SF_top * SF_antitop); break; //ttbar_atMCaNLO
		case 36: weight = sqrt( SF_top * SF_antitop); break; //ttbar_atMCaNLO_ext
		case 37: weight = sqrt( SF_top * SF_antitop); break; //ttbar_inc
		case 45: weight = sqrt( SF_top * SF_antitop); break; //ttbar_hdampUP
		case 46: weight = sqrt( SF_top * SF_antitop); break; //ttbar_hdampUP_ext
		case 47: weight = sqrt( SF_top * SF_antitop); break; //ttbar_hdampDOWN
		case 48: weight = sqrt( SF_top * SF_antitop); break; //ttbar_hdampDOWN_ext
		case 49: weight = sqrt( SF_top * SF_antitop); break; //TT_2l2nu_hdampUP
		case 50: weight = sqrt( SF_top * SF_antitop); break; //TT_2l2nu_hdampUP_ext1
		case 51: weight = sqrt( SF_top * SF_antitop); break; //TT_2l2nu_hdampUP_ext2
		case 52: weight = sqrt( SF_top * SF_antitop); break; //TT_2l2nu_hdampDOWN
                case 53: weight = sqrt( SF_top * SF_antitop); break; //TT_2l2nu_hdampDOWN_ext1
                case 54: weight = sqrt( SF_top * SF_antitop); break; //TT_2l2nu_hdampDOWN_ext2
		case 55: weight = sqrt( SF_top * SF_antitop); break; //ttbar (to hadronic, hdamp up)
		case 56: weight = sqrt( SF_top * SF_antitop); break; //ttbar (to hadronic, hdamp down)
		case 57: weight = sqrt( SF_top * SF_antitop); break; //ttbar (to semileptonic, hdamp up)
                case 58: weight = sqrt( SF_top * SF_antitop); break; //ttbar (to semileptonic, hdamp down)
		case 59: weight = sqrt( SF_top * SF_antitop); break; //ttbar (to semileptonic, hdamp down (ext))
		case 64: weight = sqrt( SF_top * SF_antitop); break; //ttbar_isr_UP
		case 65: weight = sqrt( SF_top * SF_antitop); break; //ttbar_isr_DOWN
		case 66: weight = sqrt( SF_top * SF_antitop); break; //ttbar_isr_DOWN_ext
		case 67: weight = sqrt( SF_top * SF_antitop); break; //ttbar_fsr_UP
		case 68: weight = sqrt( SF_top * SF_antitop); break; //ttbar_fsr_UP_ext
		case 69: weight = sqrt( SF_top * SF_antitop); break; //ttbar_fsr_DOWN
		case 70: weight = sqrt( SF_top * SF_antitop); break; //ttbar_fsr_DOWN_ext
		default: doubles weight_vec(SF_top.size(), 1.0); weight = weight_vec; break;
	}
 
	return weight;

  }};

  
  auto TotHTOverTotpT_doubles{[](const doubles& TotHT, const doubles& TotpT){

  	std::cout << "print 120" << std::endl;

  	floats TotHTOverTotpTOutput = TotHT / TotpT;
  	return TotHTOverTotpTOutput;

  }};

  
  auto CMSBTagSF_Function{[&SystematicInt](const floats& pts, const floats etas, const floats CSVv2Discr, bool BTagOrNot, const ints& Jet_hadronFlavour){

  	std::cout << "print 121" << std::endl;

  	doubles ResultVector{};

  	for(long unsigned int j = 0; j < Jet_hadronFlavour.size(); j++){

		CSVReader reader("./ScaleFactors/BTaggingEfficiency/CSVv2_94XSF_V2_B_F.csv");
		std::vector<std::vector<std::string> > dataList = reader.getData();
		std::vector<std::string> OutputVec{}; 
		std::vector<std::string> outputstringvec{};

		std::string number;

        	if(BTagOrNot == true){number = "1";}
		else{number = "0";}

		std::vector<std::string> CSVv2OperatingPointTest(pts.size(), number); 
		std::string MeasurementTypeString;

		if(abs(Jet_hadronFlavour.at(j)) == 0 || abs(Jet_hadronFlavour.at(j)) == 1){MeasurementTypeString = "mujets";}
		else if(abs(Jet_hadronFlavour.at(j)) == 2){MeasurementTypeString = "incl";}
		else{std::cout << "Not charm, bjet, gluon or light jet. Flavour = " << Jet_hadronFlavour.at(j) << std::endl; MeasurementTypeString = "0";}

        	std::vector<std::string> MeasurementTypeTest(pts.size(), MeasurementTypeString); 
		std::string systematic_type_string;

		switch(SystematicInt){
			case 3: systematic_type_string = "up"; break;
        		case 4: systematic_type_string = "down"; break;
 			default: systematic_type_string = "central"; break;
		}

        	std::vector<std::string> SysTypeTest(pts.size(), "central");
        	std::vector<std::string> JetFlavourTest{}; 
        	std::vector<std::string> EtaTest{};
		std::vector<std::string> PtTest{};
		std::vector<std::string> DiscrTest{};

		for(long unsigned int i = 0; i < Jet_hadronFlavour.size(); i++){

			std::stringstream ss;
                        ss << abs(Jet_hadronFlavour.at(i));
                        std::string Jet_hadronFlavourString(ss.str());
                        JetFlavourTest.push_back(Jet_hadronFlavourString);

		}


		for(long unsigned int i = 0; i < etas.size(); i++){
			std::stringstream ss;
			ss << etas.at(i);
			std::string EtaString(ss.str());
			EtaTest.push_back(EtaString);
		}


        	for(long unsigned int i = 0; i < pts.size(); i++){
                	std::stringstream ss;
                	ss << pts.at(i);
                	std::string PtString(ss.str());
                	PtTest.push_back(PtString);
        	}


        	for(long unsigned int i = 0; i < CSVv2Discr.size(); i++){

                	std::stringstream ss;
                	ss << CSVv2Discr.at(i);
                	std::string CSVv2DiscrString(ss.str());
                	DiscrTest.push_back(CSVv2DiscrString);
        	}


		std::vector<std::string> FinalOutVec{};


		//for loop to run over input events to see if the values of the variables match
		//if they do, the string in column 11 of the csv file (the equation string) is pushed back to the FinalOutVec vector
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
						std::string outputString = vec.at(10);

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
 
						FinalOutVec.push_back(outputString);
						break;
					
					}
					else{continue;}
			
                		}

        		}

		
		//if no matching row is found, a string equal to "1.0" is pushed back to FinalOutVec
		//for this case, a SF of 1.0 will be returned 
		if(i == 0 && FinalOutVec.empty()){FinalOutVec.push_back("1.0");}
		else if( FinalOutVec.size() != i+1 ){FinalOutVec.push_back("1.0");}
		else{continue;}


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


	doubles FinalVector{};

	for(long unsigned int i = 0; i < pts.size(); i++){
		FinalVector.push_back(ResultVector.at(i));
	}

	return FinalVector;
	
	//return ResultVector;

  }};

  auto CMSBTagSF{[&CMSBTagSF_Function](const floats& pts, const floats etas, const floats CSVv2Discr, const ints& Jet_hadronFlavour){

 	std::cout << "print 122" << std::endl;

/* 	std::cout << '\n' << std::endl;
        std::cout << '\n' << std::endl;
        std::cout << '\n' << std::endl;
        std::cout << "inside CMSBTagSF" << std::endl;
        std::cout << "pts.size() = " << pts.size() << std::endl;
        std::cout << "etas.size() = " << etas.size() << std::endl;
        std::cout << "CSVv2Discr.size() = " << CSVv2Discr.size() << std::endl;
        std::cout << "Jet_hadronFlavour.size() = " << Jet_hadronFlavour.size() << std::endl;
	std::cout << "pts = " << pts << std::endl;
        std::cout << "etas = " << etas << std::endl;
        std::cout << "CSVv2Discr = " << CSVv2Discr << std::endl;
        std::cout << "Jet_hadronFlavour = " << Jet_hadronFlavour << std::endl;
	std::cout << "CMSBTagSF_Function(pts, etas, CSVv2Discr, true, Jet_hadronFlavour) = " << CMSBTagSF_Function(pts, etas, CSVv2Discr, true, Jet_hadronFlavour) << std::endl;
        std::cout << '\n' << std::endl;
        std::cout << '\n' << std::endl;
        std::cout << '\n' << std::endl;
*/

 	return CMSBTagSF_Function(pts, etas, CSVv2Discr, true, Jet_hadronFlavour);

  }};

  auto nonbjet_id{[](const ints& tight_jets, const floats& btags, const floats& etas) {

  	std::cout << "print 123" << std::endl;
  	return /*tight_jets && */(btags >= 0) && (etas < MaxTrackerEta);

  }};

  auto CMSNonBTagSF{[&CMSBTagSF_Function](const floats& pts, const floats etas, const floats CSVv2Discr, const ints& Jet_hadronFlavour){

 	std::cout << "print 124" << std::endl;
/*
 	std::cout << '\n' << std::endl;
	std::cout << '\n' << std::endl;
	std::cout << '\n' << std::endl; 
 	std::cout << "inside CMSNonBTagSF" << std::endl;
	std::cout << "pts.size() = " << pts.size() << std::endl;
	std::cout << "etas.size() = " << etas.size() << std::endl;
	std::cout << "CSVv2Discr.size() = " << CSVv2Discr.size() << std::endl;
	std::cout << "Jet_hadronFlavour.size() = " << Jet_hadronFlavour.size() << std::endl;
	std::cout << "pts = " << pts << std::endl;
        std::cout << "etas = " << etas << std::endl;
        std::cout << "CSVv2Discr = " << CSVv2Discr << std::endl;
        std::cout << "Jet_hadronFlavour = " << Jet_hadronFlavour << std::endl;
	std::cout << "CMSBTagSF_Function(pts, etas, CSVv2Discr, false, Jet_hadronFlavour) = " << CMSBTagSF_Function(pts, etas, CSVv2Discr, false, Jet_hadronFlavour) << std::endl;
	std::cout << '\n' << std::endl;
	std::cout << '\n' << std::endl;
	std::cout << '\n' << std::endl;
*/	

 	return CMSBTagSF_Function(pts, etas, CSVv2Discr, false, Jet_hadronFlavour);

  }};
  
  auto EffBTaggedFunction{[&h_bjet, &h_nonbjet](
			   const int& HistOption, const floats& pts, const floats& etas, const floats& CSVv2Discr, const ints& JetFlav){

  	std::cout << "print 125" << std::endl;

  	doubles BTaggedEff{};

	double eff;
	int PtBin;
	int EtaBin;

	for(long unsigned int i = 0; i < pts.size(); i++){

			switch(HistOption){
				case 0: PtBin = h_bjet->GetYaxis()->FindBin(pts.at(i));
					EtaBin = h_bjet->GetXaxis()->FindBin(etas.at(i));
				
					eff = h_bjet->GetBinContent(EtaBin, PtBin);
					if(eff == 0){eff = 1.0;}

					break;

				case 1: PtBin = h_nonbjet->GetYaxis()->FindBin(pts.at(i));
                        		EtaBin = h_nonbjet->GetXaxis()->FindBin(etas.at(i));
                       	 	
                                	eff = h_nonbjet->GetBinContent(EtaBin, PtBin);
                                	if(eff == 0){eff = 1.0;}

                        		break;

				default: throw std::logic_error("HistOption must be 0 or 1"); 
        		}

			BTaggedEff.push_back(eff);
		
	}

 	return BTaggedEff;

  }};



  auto EffBTagged_Function{[&EffBTaggedFunction](const floats& pts, const floats& etas, const floats& CSVv2Discr, const ints& JetFlav){

  	return EffBTaggedFunction(0, pts, etas, CSVv2Discr, JetFlav);

  }};

  auto EffNonBTagged_Function{[&EffBTaggedFunction](const floats& pts, const floats& etas, const floats& CSVv2Discr, const ints& JetFlav){

        return EffBTaggedFunction(1, pts, etas, CSVv2Discr, JetFlav);

  }};

  auto ProductOperator_E_i_Function{[](const doubles& EffBTagged){
  
  	std::cout << "print 126" << std::endl;
/*
	std::cout << '\n' << std::endl;
        std::cout << '\n' << std::endl;
        std::cout << '\n' << std::endl;
	std::cout << "inside EffBTaggedProduct" << std::endl;
	std::cout << "EffBTagged = " << EffBTagged << std::endl;
	std::cout << '\n' << std::endl;
        std::cout << '\n' << std::endl;
        std::cout << '\n' << std::endl;
*/
  	double initial = 1;

	if(EffBTagged.size() > 0){
  		for(long unsigned int i = 0; i < EffBTagged.size(); i++ ){initial = EffBTagged.at(i) * initial;}
	}
	else{initial = 1 * initial;}
/*
	std::cout << '\n' << std::endl;
	std::cout << '\n' << std::endl;
	std::cout << "initial in EffBTaggedProduct = " << initial << std::endl;
	std::cout << '\n' << std::endl;
	std::cout << '\n' << std::endl;
*/
  	return initial;

  }};

  auto ProductOperator_1_Minus_E_j_Function{[](const doubles& EffNonBTagged){

  	std::cout << "print 127" << std::endl;

//	std::cout << '\n' << std::endl;
//        std::cout << '\n' << std::endl;
//        std::cout << '\n' << std::endl;
//	std::cout << "inside EffNonBTaggedProduct" << std::endl;
//	std::cout << "EffNonBTagged" << EffNonBTagged << std::endl;
//	std::cout << '\n' << std::endl;
//       std::cout << '\n' << std::endl;
//        std::cout << '\n' << std::endl; 

  	double initial = 1;

  	for(long unsigned int i = 0; i < EffNonBTagged.size(); i++ ){
		if(EffNonBTagged.at(i) == 1){initial = 1 * initial;}
		else{initial = (1 - EffNonBTagged.at(i)) * initial;}
	}

//	std::cout << '\n' << std::endl; 
// 	std::cout << '\n' << std::endl;
//      std::cout << "initial in NonEffBTaggedProduct = " << initial << std::endl;
//      std::cout << '\n' << std::endl;
//      std::cout << '\n' << std::endl;

  	return initial;

  }};

  auto ProductOperator_SFi_Times_Ei_Function{[](const doubles& EffBTagged, const doubles& CMSBTagSFInput){

  	std::cout << "print 128" << std::endl;

	if(EffBTagged.size() != CMSBTagSFInput.size()){
                std::cout << "EffBTagged = " << EffBTagged << std::endl;
                std::cout << "CMSBTagSFInput = " << CMSBTagSFInput << std::endl;
                throw std::logic_error("Eff and CMS SF vectors are not the same size");
        }

  	double initial = 1;
/*
	std::cout << '\n' << std::endl;
        std::cout << '\n' << std::endl;
        std::cout << '\n' << std::endl;
  	std::cout << "inside EffBTaggedProductData" << std::endl;
	std::cout << "EffBTagged = " << EffBTagged << std::endl;
	std::cout << "CMSBTagSFInput = " << CMSBTagSFInput << std::endl;
*/
	if(CMSBTagSFInput.size() > 0 && EffBTagged.size() > 0){
  		for(long unsigned int i = 0; i < EffBTagged.size(); i++){initial = (CMSBTagSFInput.at(i)*EffBTagged.at(i)) * initial;}
	}
	//else{throw std::logic_error("Size of btag SF vector is zero");}

/*
	std::cout << "inital = " << initial << std::endl;
	std::cout << '\n' << std::endl;
        std::cout << '\n' << std::endl;
        std::cout << '\n' << std::endl;
*/

  	return initial;

  }};



  auto ProductOperator_1_Minus_SFj_Times_Ej_Function{[](const doubles& EffNonBTagged, const doubles& CMSNonBTagSFInput){

  	std::cout << "print 129" << std::endl;

  	double initial = 1;

	if(EffNonBTagged.size() != CMSNonBTagSFInput.size()){
		std::cout << "EffNonBTagged = " << EffNonBTagged << std::endl;
        	std::cout << "CMSNonBTagSFInput = " << CMSNonBTagSFInput << std::endl;
		throw std::logic_error("Eff and CMS SF vectors are not the same size");
	}

/*
	std::cout << '\n' << std::endl;
	std::cout << '\n' << std::endl;
	std::cout << '\n' << std::endl;
	std::cout << "inside EffNonBTaggedProductData" << std::endl;
	std::cout << "EffNonBTagged = " << EffNonBTagged << std::endl;
	std::cout << "CMSNonBTagSFInput = " << CMSNonBTagSFInput << std::endl;
*/
  	for(int i = 0; i < EffNonBTagged.size(); i++){
		if(CMSNonBTagSFInput.size() > 0){
			if(CMSNonBTagSFInput.at(i) != 1 || EffNonBTagged.at(i) != 1){
				initial = (1 - (CMSNonBTagSFInput.at(i)*EffNonBTagged.at(i)) ) * initial;
			}
			else{initial = 1 * initial;}
		}
	}

/*
	std::cout << '\n' << std::endl;
        std::cout << '\n' << std::endl;
        std::cout << '\n' << std::endl;
	std::cout << "INITIAL = " << initial << std::endl;
	std::cout << '\n' << std::endl;
*/
  	return initial;

  }};


  auto ProbBTagMCFunction{[](const double& ProductOperator_E_i_Input, const double& ProductOperator_1_Minus_E_j_Input){

  	std::cout << "print 130" << std::endl;

	//std::cout << "EffBTaggedProductInput = " << EffBTaggedProductInput << std::endl;
	//std::cout << "EffNonBTaggedProductInput = " << EffNonBTaggedProductInput << std::endl;

	double MCProb = ProductOperator_E_i_Input * ProductOperator_1_Minus_E_j_Input;
  	return MCProb;

  }};


  auto ProbBTagDataFunction{[](const double& ProductOperator_SFi_Times_Ei_Input, const double& ProductOperator_1_Minus_SFj_Times_Ej_Input){

  	std::cout << "print 131" << std::endl;
 
	//std::cout << "EffBTaggedProductDataInput = " << EffBTaggedProductDataInput << std::endl;
	//std::cout << "EffNonBTaggedProductDataInput = " << EffNonBTaggedProductDataInput << std::endl;
 
  	double DataProb = ProductOperator_SFi_Times_Ei_Input * ProductOperator_1_Minus_SFj_Times_Ej_Input;
  	return DataProb;
  
  }};


  auto BTagWeightFunction{[](const double& ProbBTagMC, const double& ProbBTagData){

  	std::cout << "print 131" << std::endl;

	double BTagWeight = (ProbBTagData) / (ProbBTagMC);

        if( !isnan(BTagWeight) && !isinf(BTagWeight) && BTagWeight != 0){return BTagWeight;}
	else{throw std::logic_error("BTagWeight is either nan, infinity or zero"); /*double One = 1.0; return One;*/}

  }};

  auto EGammaFunction{[&EGammaEff2016_histo,     	     &EGammaEffSys2016_histo,
		       &EGammaEffReco2016_histo, 	     &EGammaEffRecoSys2016_histo,
		       &EGammaEff2017_histo,                 &EGammaEffSys2017_histo, 
		       &EGammaEffReco_LowPt_2017_histo,      &EGammaEffRecoSys_LowPt_2017_histo,
		       &EGammaEffReco_HigherPt_2017_histo,   &EGammaEffRecoSys_HigherPt_2017_histo,
		       &EGammaEff2018_histo,	             &EGammaEffSys2018_histo,
		       &EGammaEffReco2018_histo,	     &EGammaEffRecoSys2018_histo
		     ](const int& YearInput, const std::string& type, const floats& pt, const floats& SuperClusterEta){


  	std::cout << "print 132" << std::endl;

   	doubles OutputVector{};
   	doubles OutputVectorFinal{};

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

			double EGammaSF;

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

  	std::cout << "print 133" << std::endl;
  	return EGammaFunction(YearInt, "EGammaEff", Electron_pt_Selection, SuperClusterEta);

  }};

  auto EGammaSF_egammaEff_Sys{[&YearInt, &EGammaFunction](const floats& Electron_pt_Selection, const floats& SuperClusterEta){

  	std::cout << "print 134" << std::endl;
  	return EGammaFunction(YearInt, "EGammaEffSys", Electron_pt_Selection, SuperClusterEta);

  }};

  auto EGammaSF_egammaEffReco{[&YearInt, &EGammaFunction](const floats& Electron_pt_Selection, const floats& SuperClusterEta){

  	std::cout << "print 135" << std::endl;
  	return EGammaFunction(YearInt, "EGammaEffReco", Electron_pt_Selection, SuperClusterEta);

  }};


  auto EGammaSF_egammaEffReco_Sys{[&YearInt, &EGammaFunction](const floats& Electron_pt_Selection, const floats& SuperClusterEta){

  	std::cout << "print 136" << std::endl;
  	return EGammaFunction(YearInt, "EGammaEffRecoSys", Electron_pt_Selection, SuperClusterEta);

  }};


  auto MuonSF{[&Year, 				   &histo_RunsBCDEF_ID_2016,     &histo_RunsGH_ID_2016, 	      &histo_RunsBCDEF_ISO_2016, 	     
	       &histo_RunsGH_ISO_2016,  	   &histo_RunsBCDEF_ID_2017,     &histo_RunsBCDEF_ID_Sys_2017,  &histo_RunsBCDEF_ID_Sys_Stat_2017, 
	       &histo_RunsBCDEF_ID_Sys_Syst_2017,  &histo_RunsBCDEF_ISO_2017,    &histo_RunsBCDEF_ISO_Sys_2017, &histo_RunsBCDEF_ISO_Sys_Stat_2017,
	       &histo_RunsBCDEF_ISO_Sys_Syst_2017, &histo_RunsABCD_ID_2018,      &histo_RunsABCD_ISO_2018,      &histo_RunsABCD_ID_2018_stat,
	       &histo_RunsABCD_ID_2018_syst,       &histo_RunsABCD_ISO_2018_stat, &histo_RunsABCD_ISO_2018_syst
              ](const std::string& type, const int& YearInt, const std::string& UpOrDown, const floats& pt, const floats& eta){

  	std::cout << "print 137" << std::endl;

  	floats AbsEta = abs(eta);
  	float lumiRunBCDEF = 19713.888;
  	float lumiRunGH = 16146.178;

  	doubles MuonSFOutput{};

  	for(long unsigned int i = 0; i < pt.size(); i++){

  		if(pt.at(i) >= 20 && pt.at(i) <= 120 && AbsEta.at(i) <= MaxTrackerEta){ 

			if(YearInt == 2016){

				double MuonSF_RunsBCDEF_ID_2016 = histo_RunsBCDEF_ID_2016->GetBinContent( histo_RunsBCDEF_ID_2016->FindBin(pt.at(i), AbsEta.at(i)) );
				double MuonSF_RunsGH_ID_2016 = histo_RunsGH_ID_2016->GetBinContent( histo_RunsGH_ID_2016->FindBin(pt.at(i), AbsEta.at(i)) );
				double MuonSF_RunsBCDEF_ISO_2016 = histo_RunsBCDEF_ISO_2016->GetBinContent( histo_RunsBCDEF_ISO_2016->FindBin(pt.at(i), AbsEta.at(i)) );
				double MuonSF_RunsGH_ISO_2016 = histo_RunsGH_ISO_2016->GetBinContent( histo_RunsGH_ISO_2016->FindBin(pt.at(i), AbsEta.at(i)) );
				double Error_RunsBCDEF_ID_2016 = histo_RunsBCDEF_ID_2016->GetBinError( histo_RunsBCDEF_ID_2016->FindBin(pt.at(i), AbsEta.at(i)) );
				double Error_RunsGH_ID_2016 = histo_RunsGH_ID_2016->GetBinError( histo_RunsGH_ID_2016->FindBin(pt.at(i), AbsEta.at(i)) );
				double Error_RunsBCDEF_ISO_2016 = histo_RunsBCDEF_ISO_2016->GetBinError( histo_RunsBCDEF_ISO_2016->FindBin(pt.at(i), AbsEta.at(i)) );
				double Error_RunsGH_ISO_2016 = histo_RunsGH_ISO_2016->GetBinError( histo_RunsGH_ISO_2016->FindBin(pt.at(i), AbsEta.at(i)) );

				double Error_RunsBCDEFGH, MuonSF_RunsBCDEFGH, Error_RunsBCDEF, MuonSF_RunsBCDEF, Error_RunsGH, MuonSF_RunsGH;

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

				double MuonSF_RunsBCDEF_ID_2017 = histo_RunsBCDEF_ID_2017->GetBinContent( histo_RunsBCDEF_ID_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				double MuonSF_RunsBCDEF_ID_Sys_2017 = histo_RunsBCDEF_ID_Sys_2017->GetBinContent( histo_RunsBCDEF_ID_Sys_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				double MuonSF_RunsBCDEF_ID_Sys_Stat_2017 = histo_RunsBCDEF_ID_Sys_Stat_2017->GetBinContent( histo_RunsBCDEF_ID_Sys_Stat_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				double MuonSF_RunsBCDEF_ID_Sys_Syst_2017 = histo_RunsBCDEF_ID_Sys_Syst_2017->GetBinContent( histo_RunsBCDEF_ID_Sys_Syst_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				double MuonSF_RunsBCDEF_ISO_2017 = histo_RunsBCDEF_ISO_2017->GetBinContent( histo_RunsBCDEF_ISO_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				double MuonSF_RunsBCDEF_ISO_Sys_2017 = histo_RunsBCDEF_ISO_Sys_2017->GetBinContent( histo_RunsBCDEF_ISO_Sys_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				double MuonSF_RunsBCDEF_ISO_Sys_Stat_2017 = histo_RunsBCDEF_ISO_Sys_Stat_2017->GetBinContent( histo_RunsBCDEF_ISO_Sys_Stat_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				double MuonSF_RunsBCDEF_ISO_Sys_Syst_2017 = histo_RunsBCDEF_ISO_Sys_Syst_2017->GetBinContent( histo_RunsBCDEF_ISO_Sys_Syst_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				double Error_RunsBCDEF_ID_2017 = histo_RunsBCDEF_ID_2017->GetBinError( histo_RunsBCDEF_ID_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				double Error_RunsBCDEF_ID_Sys_2017 = histo_RunsBCDEF_ID_Sys_2017->GetBinError( histo_RunsBCDEF_ID_Sys_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				double Error_RunsBCDEF_ID_Sys_Stat_2017 = histo_RunsBCDEF_ID_Sys_Stat_2017->GetBinError( histo_RunsBCDEF_ID_Sys_Stat_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				double Error_RunsBCDEF_ID_Sys_Syst_2017 = histo_RunsBCDEF_ID_Sys_Syst_2017->GetBinError( histo_RunsBCDEF_ID_Sys_Syst_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				double Error_RunsBCDEF_ISO_2017 = histo_RunsBCDEF_ISO_2017->GetBinError( histo_RunsBCDEF_ISO_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				double Error_RunsBCDEF_ISO_Sys_2017 = histo_RunsBCDEF_ISO_Sys_2017->GetBinError( histo_RunsBCDEF_ISO_Sys_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				double Error_RunsBCDEF_ISO_Sys_Stat_2017 = histo_RunsBCDEF_ISO_Sys_Stat_2017->GetBinError( histo_RunsBCDEF_ISO_Sys_Stat_2017->FindBin(pt.at(i), AbsEta.at(i)) );
				double Error_RunsBCDEF_ISO_Sys_Syst_2017 = histo_RunsBCDEF_ISO_Sys_Syst_2017->GetBinError( histo_RunsBCDEF_ISO_Sys_Syst_2017->FindBin(pt.at(i), AbsEta.at(i)) );

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


				double MuonSF_RunsABCD_ID_2018 = histo_RunsABCD_ID_2018->GetBinContent( histo_RunsABCD_ID_2018->FindBin(pt.at(i), AbsEta.at(i)) );
				double MuonSF_RunsABCD_ISO_2018 = histo_RunsABCD_ISO_2018->GetBinContent( histo_RunsABCD_ISO_2018->FindBin(pt.at(i), AbsEta.at(i)) );
				double Error_RunsABCD_ID_2018 = histo_RunsABCD_ID_2018->GetBinError( histo_RunsABCD_ID_2018->FindBin(pt.at(i), AbsEta.at(i)) );
                                double Error_RunsABCD_ISO_2018 = histo_RunsABCD_ISO_2018->GetBinError( histo_RunsABCD_ISO_2018->FindBin(pt.at(i), AbsEta.at(i)) );

				double Error_RunsABCD_ID_2018_stat = histo_RunsABCD_ID_2018_stat->GetBinError( histo_RunsABCD_ID_2018_stat->FindBin(pt.at(i), AbsEta.at(i)) );
				double Error_RunsABCD_ISO_2018_stat = histo_RunsABCD_ISO_2018_stat->GetBinError( histo_RunsABCD_ISO_2018_stat->FindBin(pt.at(i), AbsEta.at(i)) );

				double Error_RunsABCD_ID_2018_syst = histo_RunsABCD_ID_2018_syst->GetBinError( histo_RunsABCD_ID_2018_syst->FindBin(pt.at(i), AbsEta.at(i)) );
                                double Error_RunsABCD_ISO_2018_syst = histo_RunsABCD_ISO_2018_syst->GetBinError( histo_RunsABCD_ISO_2018_syst->FindBin(pt.at(i), AbsEta.at(i)) );

				if(type == "ID"){MuonSFOutput.push_back(MuonSF_RunsABCD_ID_2018);}
				if(type == "ID sys"){MuonSFOutput.push_back(Error_RunsABCD_ID_2018);}
				else if(type == "ID sys (syst)"){MuonSFOutput.push_back(Error_RunsABCD_ID_2018_syst);}
				else if(type == "ID sys (stat)"){MuonSFOutput.push_back(Error_RunsABCD_ID_2018_stat);}
				else if(type == "Iso"){MuonSFOutput.push_back(MuonSF_RunsABCD_ISO_2018);}
				else if(type == "Iso sys"){MuonSFOutput.push_back(Error_RunsABCD_ISO_2018);}
				else if(type == "Iso sys (syst)"){MuonSFOutput.push_back(Error_RunsABCD_ISO_2018_syst);}
				else if(type == "Iso sys (stat)"){MuonSFOutput.push_back(Error_RunsABCD_ISO_2018_stat);}
				else{std::cout << "Error with Muon SF type (2018)" << std::endl;}

			}
			else{throw std::logic_error("Year must be 2016, 2017 or 2018");}

		}
		else{double One = 1.0; MuonSFOutput.push_back(One);}


  	}

  	return MuonSFOutput.at(0); 

  }};

  auto MuonSFTest_ID{[&MuonSF, &YearInt](const floats& pt, const floats& eta){

  	std::cout << "print 138" << std::endl;
  	return MuonSF("ID", YearInt, " ", pt, eta);
  
  }};

  auto MuonSFTest_Iso{[&MuonSF, &YearInt](const floats& pt, const floats& eta){

  	std::cout << "print 139" << std::endl;
  	return MuonSF("Iso", YearInt, " ", pt, eta);
    
  }};

  auto MuonSFTest_ID_sys_syst{[&MuonSF, &YearInt](const floats& pt, const floats& eta){
  
  	std::cout << "print 140" << std::endl;

	switch(YearInt){
		case 2016: return MuonSF("ID sys", YearInt, "Up", pt, eta);
        	case 2017: return MuonSF("ID sys (syst)", YearInt, " ", pt, eta);
  		case 2018: return MuonSF("ID sys (syst)", YearInt, " ", pt, eta); 
		default: throw std::logic_error("Year must be 2016, 2017 or 2018"); break;
	}
		
  }};


  auto MuonSFTest_ID_sys_stat{[&MuonSF, &YearInt](const floats& pt, const floats& eta){

  	std::cout << "print 141" << std::endl;

	switch(YearInt){
		case 2016: return MuonSF("ID sys", YearInt, "Down", pt, eta);
  		case 2017: return MuonSF("ID sys (stat)", YearInt, " ", pt, eta);
  		case 2018: return MuonSF("ID sys (stat)", YearInt, " ", pt, eta);
		default: throw std::logic_error("Year must be 2016, 2017 or 2018"); break; 
	}

  }};


  auto MuonSFTest_Iso_sys_syst{[&MuonSF, &YearInt](const floats& pt, const floats& eta){

  	std::cout << "print 142" << std::endl;

	switch(YearInt){
		case 2016: return MuonSF("Iso sys", YearInt, "Up", pt, eta);
        	case 2017: return MuonSF("Iso sys (syst)", YearInt, " ", pt, eta);
  		case 2018: return MuonSF("Iso sys (syst)", YearInt, " ", pt, eta); 
		default: throw std::logic_error("Year must be 2016, 2017 or 2018"); break;
	}

  }};


  auto MuonSFTest_Iso_sys_stat{[&MuonSF, &YearInt](const floats& pt, const floats& eta){

  	std::cout << "print 143" << std::endl;

	switch(YearInt){
		case 2016: return MuonSF("Iso sys", YearInt, "Down", pt, eta);
        	case 2017: return MuonSF("Iso sys (stat)", YearInt, " ", pt, eta);
  		case 2018: return MuonSF("Iso sys (stat)", YearInt, " ", pt, eta);
		default: throw std::logic_error("Year must be 2016, 2017 or 2018"); break;
	}

  }};




  auto PSWeightFunction{[&YearInt, &ProcessInt](const floats& PSWeightInput){

  	std::cout << "print 144" << std::endl;

  	doubles Ones(4, 1.0);
	doubles PSWeightInput_Doubles;

	for(int i = 0; i < PSWeightInput.size(); i++){PSWeightInput_Doubles.push_back(PSWeightInput.at(i));}

	switch(YearInt){

		case 2016: return Ones;
		case 2017: switch(ProcessInt){
                        	case 0: return PSWeightInput_Doubles; break; //tZq
				case 28: return PSWeightInput_Doubles; break; //ttbar (2l2nu)
                        	case 31: return PSWeightInput_Doubles; break; //ttbar (to hadronic)
                        	case 32: return PSWeightInput_Doubles; break; //ttbar (to hadronic ext)
                        	case 33: return PSWeightInput_Doubles; break; //ttbar (to semileptonic)
                        	case 34: return PSWeightInput_Doubles; break; //ttbar (to semileptonic ext)
                        	case 38: return PSWeightInput_Doubles; break; //single top t-channel (top)
                        	case 41: return PSWeightInput_Doubles; break; //single top t-channel (antitop)
                        	case 44: return PSWeightInput_Doubles; break; //single top s-channel
                        	case 71: return PSWeightInput_Doubles; break; //single top (tW)
                        	case 74: return PSWeightInput_Doubles; break; //single top (tbarW)
				case 90: return PSWeightInput_Doubles; break; //WWTo2L2Nu
				case 106: return PSWeightInput_Doubles; break; //ttWToLNu
                                case 107: return PSWeightInput_Doubles; break; //ttWToLNu_ext
				case 114: return PSWeightInput_Doubles; break; //ttZToQQ
				case 115: return PSWeightInput_Doubles; break; //ttZToQQ_ext
				default: return Ones;			    
			}
		case 2018: switch(ProcessInt){
				case 0: return PSWeightInput_Doubles; break; //tZq
				case 28: return PSWeightInput_Doubles; break; //ttbar (2l2nu)
                                case 31: return PSWeightInput_Doubles; break; //ttbar (to hadronic)
                                case 32: return PSWeightInput_Doubles; break; //ttbar (to hadronic ext)
                                case 33: return PSWeightInput_Doubles; break; //ttbar (to semileptonic)
                                case 34: return PSWeightInput_Doubles; break; //ttbar (to semileptonic ext)
                                case 38: return PSWeightInput_Doubles; break; //single top t-channel (top)
                                case 41: return PSWeightInput_Doubles; break; //single top t-channel (antitop)
                                case 44: return PSWeightInput_Doubles; break; //single top s-channel
                                case 71: return PSWeightInput_Doubles; break; //single top (tW)
                                case 74: return PSWeightInput_Doubles; break; //single top (tbarW)
				case 90: return PSWeightInput_Doubles; break; //WWTo2L2Nu
				case 106: return PSWeightInput_Doubles; break; //ttWToLNu
                                case 107: return PSWeightInput_Doubles; break; //ttWToLNu_ext
				case 114: return PSWeightInput_Doubles; break; //ttZToQQ
                                case 115: return PSWeightInput_Doubles; break; //ttZToQQ_ext
                                default: return Ones;   
 
			}

	}


  }};

  
  auto PDFWeight{[&MCInt, &SystematicInt, &HessianOrMC](const floats& LHEPdfWeight, const unsigned int& nLHEPdfWeight){

  	std::cout << "print 145" << std::endl;

	double PdfUncert;
	double NominalPdfWeight = LHEPdfWeight.at(0);

  	//For the up and down PDF uncertainties
  	
	switch(MCInt){

		case 0: double One; return One; break;

		case 1:

			if(HessianOrMC == "Hessian"){
  				//For Hessian PDF sets
  				double PdfUncert_Hessian_Squared;
  	
  				for(unsigned int k = 1; k < nLHEPdfWeight; k++){PdfUncert_Hessian_Squared += pow((LHEPdfWeight.at(k) - LHEPdfWeight.at(0)), 2);}

				double PdfUncert_Hessian = sqrt(PdfUncert_Hessian_Squared);
				PdfUncert = PdfUncert_Hessian;
			}
			else if(HessianOrMC == "MC"){
				//For MC PDF sets
				double SumOfLHEPdfWeights;
	
				for(unsigned int k = 1; k < nLHEPdfWeight; k++){SumOfLHEPdfWeights += LHEPdfWeight.at(k);}

				double MeanPdfWeight = (1/nLHEPdfWeight) * SumOfLHEPdfWeights;	
				double PartOf_PdfUncert_MC;

				for(unsigned int k = 1; k < nLHEPdfWeight; k++){PartOf_PdfUncert_MC += pow((LHEPdfWeight.at(k) - MeanPdfWeight), 2);}

				double PdfUncert_MC = sqrt( ( 1/(nLHEPdfWeight-1) ) * PartOf_PdfUncert_MC ); 
				PdfUncert = PdfUncert_MC;
			}
			else{throw std::logic_error("Please set the HessianOrMC string to be Hessian or MC.");}

			switch(SystematicInt){
				case 11: return NominalPdfWeight + PdfUncert; break;
				case 12: return NominalPdfWeight - PdfUncert; break;
				default: return NominalPdfWeight; break;
			}

			break;

		default: throw std::logic_error("MCInt must be 0 for data or 1 for MC"); break;

  	}

  }};


  double SumOfSigns_GenWeight;
  double SumOfGenWeights;
  double gen_weightSF;

  auto GeneratorWeight{[&MCInt, &SystematicInt, &gen_weightSF](const float& InputGenWeight){

	std::cout << "print 147" << std::endl;

	double CalculatedGenWeight;

        switch(MCInt){

		case 0: CalculatedGenWeight = 1.0; break;
		
		case 1: CalculatedGenWeight = gen_weightSF;

			break;

		default: throw std::logic_error("MCInt must be 0 or 1"); break;

	}

	return CalculatedGenWeight;

  }};

  auto GeneratorWeightFilterFunction{[&YearInt, &ProcessInt](const float& InputGenWeight){

	if(  (YearInt == 2016 && (ProcessInt == 5 || ProcessInt == 6 || ProcessInt == 9 || ProcessInt == 29 || ProcessInt == 30 || ProcessInt == 77 || ProcessInt == 79) ) ||
       (YearInt == 2017 && (ProcessInt == 77 || ProcessInt == 79 || ProcessInt == 93 || ProcessInt == 104 || ProcessInt == 105) ) ||
       (YearInt == 2018 && (ProcessInt == 0 || ProcessInt == 9 || ProcessInt == 10 || ProcessInt == 77 || ProcessInt == 79 || ProcessInt == 93 || ProcessInt == 104))){


		if(InputGenWeight > 0){return true;}
		else{return false;}

	}
	else{return true;}

  }};


  auto OriginalMetFunction{[&SystematicInt](const floats& MET_sumEt, const floats& MET_phi){

	std::cout << "print 163" << std::endl; 

	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> OriginalMET{};
	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> OriginalMET_Element{};
	
	for(int i = 0; i < MET_phi.size(); i++){

        	OriginalMET_Element.SetPt(MET_sumEt.at(i));
		OriginalMET_Element.SetEta(0);
		OriginalMET_Element.SetPhi(MET_phi.at(i));
		OriginalMET_Element.SetM(MET_sumEt.at(i));

        	OriginalMET += OriginalMET_Element;
	
	}

        return OriginalMET;

  }};


  auto ScaledMetFunction{[&SystematicInt](ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> OriginalMET, const floats& MET_sumEt, const floats& MET_phi, const floats& MET_MetUnclustEnUpDeltaX,  const floats& MET_MetUnclustEnUpDeltaY){

	std::cout << "print 164" << std::endl;

	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> ScaledMET{};
	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> ScaledMET_Element{};
	floats metVecOriginal_px;
	floats metVecOriginal_py;

	metVecOriginal_px.push_back( OriginalMET.Px() );
	metVecOriginal_py.push_back( OriginalMET.Py() );
	
	floats MET_px_up =  metVecOriginal_px + MET_MetUnclustEnUpDeltaX;
        floats MET_py_up =  metVecOriginal_py + MET_MetUnclustEnUpDeltaY;
        floats MET_px_down =  metVecOriginal_px - MET_MetUnclustEnUpDeltaX;
        floats MET_py_down =  metVecOriginal_py - MET_MetUnclustEnUpDeltaY;

        floats UnclusteredEnergyUp = sqrt( pow(MET_px_up, 2) + pow(MET_py_up, 2) );
        floats UnclusteredEnergyDown = sqrt( pow(MET_px_down, 2) + pow(MET_py_down, 2) );

        for(long unsigned int i = 0; i < MET_phi.size(); i++){

		if(SystematicInt == 15){

			ScaledMET_Element.SetPt(UnclusteredEnergyUp.at(i));
			ScaledMET_Element.SetEta(0);
			ScaledMET_Element.SetPhi(MET_phi.at(i));
			ScaledMET_Element.SetM(UnclusteredEnergyUp.at(i));

		}
                else if(SystematicInt == 16){

			ScaledMET_Element.SetPt(UnclusteredEnergyDown.at(i));
			ScaledMET_Element.SetEta(0);
			ScaledMET_Element.SetPhi(MET_phi.at(i));
			ScaledMET_Element.SetM(UnclusteredEnergyDown.at(i));

		}
                else{
			ScaledMET_Element.SetPt(MET_sumEt.at(i));
			ScaledMET_Element.SetEta(0);
			ScaledMET_Element.SetPhi(MET_phi.at(i));
			ScaledMET_Element.SetM(MET_sumEt.at(i));
		}

                ScaledMET += ScaledMET_Element;

        }	

	return ScaledMET;

  }};

  auto UnsmearedJetTLorentzVectorFunction{[](const floats& Jet_pt, const floats& Jet_phi, const floats& Jet_eta, const floats& Jet_mass){


	std::cout << "print 165" << std::endl;

  	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> UnsmearedJetVector{};
	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> UnsmearedJetVector_Element{};	

	for(int i = 0; i < Jet_pt.size(); i++){

		UnsmearedJetVector_Element.SetPt(Jet_pt.at(i));
		UnsmearedJetVector_Element.SetEta(Jet_eta.at(i));
		UnsmearedJetVector_Element.SetPhi(Jet_phi.at(i));
		UnsmearedJetVector_Element.SetM(Jet_mass.at(i));

		UnsmearedJetVector += UnsmearedJetVector_Element;
	}

	return UnsmearedJetVector;
  }};


  auto METUncertFunction{[&SystematicInt](ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> OriginalMET, ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> SmearedJet4Momentum, 
				          ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> UnsmearedJet4Momentum){

  	std::cout << "print 148" << std::endl;

	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> NewMetVector{};
	ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>> FinalVec{};

	double OriginalMET_Pt = OriginalMET.Pt();	
	double OriginalMET_Eta = OriginalMET.Eta();
	double OriginalMET_Phi = OriginalMET.Phi();
	double OriginalMET_M = OriginalMET.M();		

	double SmearedJet4Momentum_Pt = SmearedJet4Momentum.Pt();
	double SmearedJet4Momentum_Eta = SmearedJet4Momentum.Eta();
	double SmearedJet4Momentum_Phi = SmearedJet4Momentum.Phi();
	double SmearedJet4Momentum_M = SmearedJet4Momentum.M();

	double UnsmearedJet4Momentum_Pt = UnsmearedJet4Momentum.Pt();
        double UnsmearedJet4Momentum_Eta = UnsmearedJet4Momentum.Eta();
        double UnsmearedJet4Momentum_Phi = UnsmearedJet4Momentum.Phi();
        double UnsmearedJet4Momentum_M = UnsmearedJet4Momentum.M();

	FinalVec.SetM(OriginalMET_M + SmearedJet4Momentum_M - UnsmearedJet4Momentum_M);
	FinalVec.SetEta(OriginalMET_Eta + SmearedJet4Momentum_Eta - UnsmearedJet4Momentum_Eta);
	FinalVec.SetPhi(OriginalMET_Phi + SmearedJet4Momentum_Phi - UnsmearedJet4Momentum_Phi);
	FinalVec.SetM(OriginalMET_M + SmearedJet4Momentum_M - UnsmearedJet4Momentum_M);

	NewMetVector += FinalVec;

	return NewMetVector;

  }};




  auto linereader{[&Year](const int& LineNumber, const std::string YearChoice){
        
        std::cout << "print 150" << std::endl;
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
        
        std::cout << "print 151" << std::endl;

	double OutputNormFactor;

        switch(MCInt){
                case 1: OutputNormFactor = linereader(ProcessInt+1, Year); return OutputNormFactor;
                default: double one = 1.0; return one;
        }
  }};


  auto EventWeight{[&NormalisationFactorFunction, &ChannelInt, &SystematicInt, &ttbarCRInt]
		    (const double& PUInput, 		         const double& BTagWeightInput, 	       const doubles& ReturnedPSWeightInput, 
		     const double& EGammaSF_egammaEffInput,      const double& EGammaSF_egammaEffRecoInput, 
		     const double& EGammaSF_egammaEffSysInput,   const double& EGammaSF_egammaEffRecoSysInput, const double& CalculatedGeneratorWeightInput, 
		     const doubles& TopWeightInput, 	         const double& CalculatedPDFWeightInput, 
		     const double& MuonSFTest_IDInput, 		 const double& MuonSFTest_IsoInput, 	       const double& MuonSFTest_ID_sys_systInput, 
		     const double& MuonSFTest_ID_sys_statInput,  const double& MuonSFTest_Iso_sys_systInput,   const double& MuonSFTest_Iso_sys_statInput){


			std::cout << "print 149" << std::endl;

			double EventWeightOutput;

			switch(ChannelInt){
				case 1: switch(SystematicInt){ 
						case 9: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * 
							(TrigSF += TrigSF_Uncert) * CalculatedPDFWeightInput * EGammaSF_egammaEffSysInput * 
							EGammaSF_egammaEffRecoSysInput * CalculatedGeneratorWeightInput * TopWeightInput.at(0);

							break;

                        			case 10: EventWeightOutput = PUInput * NormalisationFactorFunction() * (TrigSF -= TrigSF_Uncert) * 
					  		 CalculatedPDFWeightInput * 
						         EGammaSF_egammaEffSysInput * EGammaSF_egammaEffRecoSysInput * CalculatedGeneratorWeightInput * 
							 TopWeightInput.at(0);

							 break;
				
						case 11: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF * 
							 CalculatedPDFWeightInput * 
							 EGammaSF_egammaEffInput * EGammaSF_egammaEffRecoInput * CalculatedGeneratorWeightInput * TopWeightInput.at(0);
                        
							 break;

						case 12: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF * 
						         		     CalculatedPDFWeightInput * EGammaSF_egammaEffInput * EGammaSF_egammaEffRecoInput * 
									     CalculatedGeneratorWeightInput * TopWeightInput.at(0);
				
							 break;

						case 13: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF *
                                                         ReturnedPSWeightInput.at(2) * ReturnedPSWeightInput.at(3) * CalculatedPDFWeightInput *
                                                         EGammaSF_egammaEffInput * EGammaSF_egammaEffRecoInput * CalculatedGeneratorWeightInput * TopWeightInput.at(0);

							 break;

						case 14: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF *
                                                         ReturnedPSWeightInput.at(1) * ReturnedPSWeightInput.at(4) * CalculatedPDFWeightInput *
                                                         EGammaSF_egammaEffInput * EGammaSF_egammaEffRecoInput * CalculatedGeneratorWeightInput * TopWeightInput.at(0);

                                                         break;

						case 17: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF * ReturnedPSWeightInput.at(2) * 
							 CalculatedPDFWeightInput * EGammaSF_egammaEffInput * EGammaSF_egammaEffRecoInput * CalculatedGeneratorWeightInput * 
							 TopWeightInput.at(0);

							 break;

                        			case 18: EventWeightOutput =  PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF * 
							 ReturnedPSWeightInput.at(0) * 
							 CalculatedPDFWeightInput * EGammaSF_egammaEffInput * EGammaSF_egammaEffRecoInput * CalculatedGeneratorWeightInput * TopWeightInput.at(0);

							 break;

                        			case 19: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF * ReturnedPSWeightInput.at(3) * 
							 CalculatedPDFWeightInput * EGammaSF_egammaEffInput * EGammaSF_egammaEffRecoInput * CalculatedGeneratorWeightInput  
							 * TopWeightInput.at(0);

							 break;

                        			case 20: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF * ReturnedPSWeightInput.at(1) * 
							 CalculatedPDFWeightInput * EGammaSF_egammaEffInput * EGammaSF_egammaEffRecoInput * CalculatedGeneratorWeightInput  
							 * TopWeightInput.at(0);

							 break;

						default: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF * 
							 CalculatedPDFWeightInput * 
							 EGammaSF_egammaEffInput * EGammaSF_egammaEffRecoInput * CalculatedGeneratorWeightInput * TopWeightInput.at(0);
					
							 break;		

					}
 
					break;
 
				case 2: switch(SystematicInt){
						case 9: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * (TrigSF += TrigSF_Uncert) * 
							CalculatedPDFWeightInput * MuonSFTest_IDInput * MuonSFTest_IsoInput * CalculatedGeneratorWeightInput * TopWeightInput.at(0);

							break;

                        			case 10: EventWeightOutput = PUInput * NormalisationFactorFunction() * (TrigSF -= TrigSF_Uncert) * 
							 CalculatedPDFWeightInput * MuonSFTest_IDInput * MuonSFTest_IsoInput * CalculatedGeneratorWeightInput * 
							 TopWeightInput.at(0);
				
							 break;

						case 11: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF * 
							 CalculatedPDFWeightInput * 
							 MuonSFTest_IDInput * MuonSFTest_IsoInput * CalculatedGeneratorWeightInput * TopWeightInput.at(0);

							 break;

                        			case 12: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF * 
						         CalculatedPDFWeightInput * 
							 MuonSFTest_IDInput * MuonSFTest_IsoInput * CalculatedGeneratorWeightInput * TopWeightInput.at(0);

							 break;
	

						case 13: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF *
                                                         CalculatedPDFWeightInput * ReturnedPSWeightInput.at(2) * ReturnedPSWeightInput.at(3) *
                                                         MuonSFTest_IDInput * MuonSFTest_IsoInput * CalculatedGeneratorWeightInput * TopWeightInput.at(0);

							 break;

						case 14: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF *
                                                         CalculatedPDFWeightInput * ReturnedPSWeightInput.at(1) * ReturnedPSWeightInput.at(4) *
                                                         MuonSFTest_IDInput * MuonSFTest_IsoInput * CalculatedGeneratorWeightInput * TopWeightInput.at(0);

                                                         break;

                        			case 17: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF * 
							 ReturnedPSWeightInput.at(2) * 
							 CalculatedPDFWeightInput * MuonSFTest_IDInput * MuonSFTest_IsoInput * CalculatedGeneratorWeightInput * 
							 TopWeightInput.at(0);

							 break;

						case 18: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF *
							 ReturnedPSWeightInput.at(0) * 
							 CalculatedPDFWeightInput * MuonSFTest_IDInput * MuonSFTest_IsoInput * CalculatedGeneratorWeightInput * 
						         TopWeightInput.at(0);

							 break; 

                        			case 19: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF * 
							 ReturnedPSWeightInput.at(3) * 
							 CalculatedPDFWeightInput * MuonSFTest_IDInput * MuonSFTest_IsoInput * CalculatedGeneratorWeightInput * 
							 TopWeightInput.at(0);

							 break;

                        			case 20: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF * 
							 ReturnedPSWeightInput.at(1) * 
							 CalculatedPDFWeightInput * MuonSFTest_IDInput * MuonSFTest_IsoInput * CalculatedGeneratorWeightInput * 
							 TopWeightInput.at(0);

							 break;

                        			default: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF * 
							 CalculatedPDFWeightInput * 
							 MuonSFTest_IDInput * MuonSFTest_IsoInput * CalculatedGeneratorWeightInput * TopWeightInput.at(0);

							 break;

					}

					break;

					case 3: switch(SystematicInt){ 
						
							case 9: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * 
								(TrigSF += TrigSF_Uncert) * CalculatedPDFWeightInput * EGammaSF_egammaEffSysInput * 
								EGammaSF_egammaEffRecoSysInput * MuonSFTest_ID_sys_systInput * MuonSFTest_Iso_sys_systInput * CalculatedGeneratorWeightInput * TopWeightInput.at(0);

							break;

                        				case 10: EventWeightOutput = PUInput * NormalisationFactorFunction() * (TrigSF -= TrigSF_Uncert) * 
					  		 	 CalculatedPDFWeightInput * MuonSFTest_ID_sys_systInput * MuonSFTest_Iso_sys_systInput *
						         	 EGammaSF_egammaEffSysInput * EGammaSF_egammaEffRecoSysInput * CalculatedGeneratorWeightInput * 
							 	 TopWeightInput.at(0);

							 break;
				
							case 11: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF * 
							  	 CalculatedPDFWeightInput * MuonSFTest_ID_sys_systInput * MuonSFTest_Iso_sys_systInput *
							 	 EGammaSF_egammaEffInput * EGammaSF_egammaEffRecoInput * CalculatedGeneratorWeightInput * TopWeightInput.at(0);
                        
							 break;

							case 12: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF * 
						         		     CalculatedPDFWeightInput * EGammaSF_egammaEffInput * EGammaSF_egammaEffRecoInput * 
									     MuonSFTest_ID_sys_systInput * MuonSFTest_Iso_sys_systInput *
									     CalculatedGeneratorWeightInput * TopWeightInput.at(0);
				
							 break;

							case 13: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF *
                                                          	 ReturnedPSWeightInput.at(2) * ReturnedPSWeightInput.at(3) * CalculatedPDFWeightInput *
                                                                 EGammaSF_egammaEffInput * EGammaSF_egammaEffRecoInput * MuonSFTest_ID_sys_systInput * MuonSFTest_Iso_sys_systInput *								      CalculatedGeneratorWeightInput * TopWeightInput.at(0);

							 break;

							case 14: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF *
                                                                 ReturnedPSWeightInput.at(1) * ReturnedPSWeightInput.at(4) * CalculatedPDFWeightInput *
                                                                 EGammaSF_egammaEffInput * EGammaSF_egammaEffRecoInput * MuonSFTest_ID_sys_systInput * MuonSFTest_Iso_sys_systInput * CalculatedGeneratorWeightInput * TopWeightInput.at(0);

                                                         break;

							case 17: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF * ReturnedPSWeightInput.at(2) * 
							         CalculatedPDFWeightInput * EGammaSF_egammaEffInput * EGammaSF_egammaEffRecoInput * MuonSFTest_ID_sys_systInput * MuonSFTest_Iso_sys_systInput * CalculatedGeneratorWeightInput * 
							         TopWeightInput.at(0);

							 break;

                        				case 18: EventWeightOutput =  PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF * 
							         ReturnedPSWeightInput.at(0) * 
							         CalculatedPDFWeightInput * EGammaSF_egammaEffInput * EGammaSF_egammaEffRecoInput * MuonSFTest_ID_sys_systInput * MuonSFTest_Iso_sys_systInput * CalculatedGeneratorWeightInput * TopWeightInput.at(0);

							 break;

                        				case 19: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF * ReturnedPSWeightInput.at(3) * 
							 	 CalculatedPDFWeightInput * EGammaSF_egammaEffInput * EGammaSF_egammaEffRecoInput * MuonSFTest_ID_sys_systInput * MuonSFTest_Iso_sys_systInput * CalculatedGeneratorWeightInput  
							 	 * TopWeightInput.at(0);

							 break;

                        				case 20: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF * ReturnedPSWeightInput.at(1) * 
							 	 CalculatedPDFWeightInput * EGammaSF_egammaEffInput * EGammaSF_egammaEffRecoInput * MuonSFTest_ID_sys_systInput * MuonSFTest_Iso_sys_systInput * CalculatedGeneratorWeightInput  
							         * TopWeightInput.at(0);

							 break;

							default: EventWeightOutput = PUInput * NormalisationFactorFunction() * BTagWeightInput * TrigSF * 
							 	 CalculatedPDFWeightInput * 
							         EGammaSF_egammaEffInput * EGammaSF_egammaEffRecoInput * MuonSFTest_ID_sys_systInput * MuonSFTest_Iso_sys_systInput * CalculatedGeneratorWeightInput * TopWeightInput.at(0);
					
							 break;		
  
					}

					break; 

					default: throw std::logic_error("ChannelInt must be 1 (for ee), 2 (for mumu) or 3 (for emu)."); break;

		}


	std::cout << '\n' << std::endl;
	std::cout << '\n' << std::endl;
	std::cout << "PUInput = " << PUInput << std::endl;
        std::cout << "NormalisationFactorFunction() = " << NormalisationFactorFunction() << std::endl;
        std::cout << "BTagWeightInput = " << BTagWeightInput << std::endl;
        std::cout << "TrigSF = " << TrigSF << std::endl;
        std::cout << "CalculatedPDFWeightInput = " <<  CalculatedPDFWeightInput << std::endl;
        std::cout << "EGammaSF_egammaEffInput = " << EGammaSF_egammaEffInput << std::endl;
        std::cout << "EGammaSF_egammaEffRecoInput = " << EGammaSF_egammaEffRecoInput << std::endl;
        std::cout << "CalculatedGeneratorWeightInput = " << CalculatedGeneratorWeightInput << std::endl;
        std::cout << "TopWeightInput.at(0) = " << TopWeightInput.at(0) << std::endl;
	std::cout << "EventWeightOutput = " << EventWeightOutput << std::endl;
	std::cout << '\n' << std::endl;
        std::cout << '\n' << std::endl;	

	if(!isnan(EventWeightOutput) && !isinf(EventWeightOutput) && (EventWeightOutput > 0)){return EventWeightOutput;}
	else{std::cout << "Final event weight is either a nan, inf or 0." << std::endl; double One = 1.0; return One;}     

  }};



  auto OverallNormalisationFunction{[&NormalisationFactorFunction, &ChannelInt, &SystematicInt, &ttbarCRInt]
		    		    (const double& PUInput, 		         const double& BTagWeightInput, 	       const doubles& ReturnedPSWeightInput, 
		     		     const double& EGammaSF_egammaEffInput,      const double& EGammaSF_egammaEffRecoInput, 
		     		     const double& EGammaSF_egammaEffSysInput,   const double& EGammaSF_egammaEffRecoSysInput, const double& CalculatedGeneratorWeightInput, 
		     		     const doubles& TopWeightInput, 	         const double& CalculatedPDFWeightInput, 
		                     const double& MuonSFTest_IDInput, 		 const double& MuonSFTest_IsoInput, 	       const double& MuonSFTest_ID_sys_systInput, 
		                     const double& MuonSFTest_ID_sys_statInput,  const double& MuonSFTest_Iso_sys_systInput,   const double& MuonSFTest_Iso_sys_statInput){


			std::cout << "print 149" << std::endl;

			double OverallNormalisationOutput;

			switch(ChannelInt){
				case 1: switch(SystematicInt){ 
						case 9: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + 
							(TrigSF += TrigSF_Uncert) + CalculatedPDFWeightInput + EGammaSF_egammaEffSysInput + 
							EGammaSF_egammaEffRecoSysInput + CalculatedGeneratorWeightInput + TopWeightInput.at(0);

							break;

                        			case 10: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + (TrigSF -= TrigSF_Uncert) + 
					  		 CalculatedPDFWeightInput + 
						         EGammaSF_egammaEffSysInput + EGammaSF_egammaEffRecoSysInput + CalculatedGeneratorWeightInput + 
							 TopWeightInput.at(0);

							 break;
				
						case 11: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + TrigSF + 
							 CalculatedPDFWeightInput + 
							 EGammaSF_egammaEffInput + EGammaSF_egammaEffRecoInput + CalculatedGeneratorWeightInput + TopWeightInput.at(0);
                        
							 break;

						case 12: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + TrigSF + 
						         		     CalculatedPDFWeightInput + EGammaSF_egammaEffInput + EGammaSF_egammaEffRecoInput + 
									     CalculatedGeneratorWeightInput + TopWeightInput.at(0);
				
							 break;

						case 13: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + TrigSF +
                                                         (ReturnedPSWeightInput.at(2) + ReturnedPSWeightInput.at(3)) + CalculatedPDFWeightInput +
                                                         EGammaSF_egammaEffInput + EGammaSF_egammaEffRecoInput + CalculatedGeneratorWeightInput + TopWeightInput.at(0);

							 break;

						case 14: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + TrigSF +
                                                         (ReturnedPSWeightInput.at(1) * ReturnedPSWeightInput.at(4)) + CalculatedPDFWeightInput +
                                                         EGammaSF_egammaEffInput + EGammaSF_egammaEffRecoInput + CalculatedGeneratorWeightInput + TopWeightInput.at(0);

                                                         break;

						case 17: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + TrigSF + ReturnedPSWeightInput.at(2) + 
							 CalculatedPDFWeightInput + EGammaSF_egammaEffInput + EGammaSF_egammaEffRecoInput + CalculatedGeneratorWeightInput + 
							 TopWeightInput.at(0);

							 break;

                        			case 18: OverallNormalisationOutput =  PUInput + NormalisationFactorFunction() + BTagWeightInput + TrigSF + 
							 ReturnedPSWeightInput.at(0) + 
							 CalculatedPDFWeightInput + EGammaSF_egammaEffInput + EGammaSF_egammaEffRecoInput + CalculatedGeneratorWeightInput + TopWeightInput.at(0);

							 break;

                        			case 19: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + TrigSF + ReturnedPSWeightInput.at(3) + 
							 CalculatedPDFWeightInput + EGammaSF_egammaEffInput + EGammaSF_egammaEffRecoInput + CalculatedGeneratorWeightInput  
							 + TopWeightInput.at(0);

							 break;

                        			case 20: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + TrigSF + ReturnedPSWeightInput.at(1) + 
							 CalculatedPDFWeightInput + EGammaSF_egammaEffInput + EGammaSF_egammaEffRecoInput + CalculatedGeneratorWeightInput  
							 + TopWeightInput.at(0);

							 break;

						default: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + TrigSF + 
							 CalculatedPDFWeightInput + 
							 EGammaSF_egammaEffInput + EGammaSF_egammaEffRecoInput + CalculatedGeneratorWeightInput + TopWeightInput.at(0);
					
							 break;		

					}
 
					break;
 
				case 2: switch(SystematicInt){
						case 9: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + (TrigSF += TrigSF_Uncert) + 
							CalculatedPDFWeightInput + MuonSFTest_ID_sys_systInput + MuonSFTest_Iso_sys_systInput + 
							CalculatedGeneratorWeightInput + TopWeightInput.at(0);

							break;

                        			case 10: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + (TrigSF -= TrigSF_Uncert) + 
							 CalculatedPDFWeightInput + 
							 MuonSFTest_ID_sys_statInput + MuonSFTest_Iso_sys_statInput + CalculatedGeneratorWeightInput + 
							 TopWeightInput.at(0);
				
							 break;

						case 11: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + TrigSF + 
							 CalculatedPDFWeightInput + 
							 MuonSFTest_IDInput + MuonSFTest_IsoInput + CalculatedGeneratorWeightInput + TopWeightInput.at(0);

							 break;

                        			case 12: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + TrigSF + 
						         CalculatedPDFWeightInput + 
							 MuonSFTest_IDInput + MuonSFTest_IsoInput + CalculatedGeneratorWeightInput + TopWeightInput.at(0);

							 break;
	

						case 13: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + TrigSF +
                                                         CalculatedPDFWeightInput + (ReturnedPSWeightInput.at(2) * ReturnedPSWeightInput.at(3)) +
                                                         MuonSFTest_IDInput + MuonSFTest_IsoInput + CalculatedGeneratorWeightInput + TopWeightInput.at(0);

							 break;

						case 14: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + TrigSF +
                                                         CalculatedPDFWeightInput + (ReturnedPSWeightInput.at(1) * ReturnedPSWeightInput.at(4)) +
                                                         MuonSFTest_IDInput + MuonSFTest_IsoInput + CalculatedGeneratorWeightInput + TopWeightInput.at(0);

                                                         break;

                        			case 17: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + TrigSF + 
							 ReturnedPSWeightInput.at(2) + 
							 CalculatedPDFWeightInput + MuonSFTest_IDInput + MuonSFTest_IsoInput + CalculatedGeneratorWeightInput + 
							 TopWeightInput.at(0);

							 break;

						case 18: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + TrigSF +
							 ReturnedPSWeightInput.at(0) + 
							 CalculatedPDFWeightInput + MuonSFTest_IDInput + MuonSFTest_IsoInput + CalculatedGeneratorWeightInput + 
						         TopWeightInput.at(0);

							 break; 

                        			case 19: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + TrigSF + 
							 ReturnedPSWeightInput.at(3) + 
							 CalculatedPDFWeightInput + MuonSFTest_IDInput + MuonSFTest_IsoInput + CalculatedGeneratorWeightInput + 
							 TopWeightInput.at(0);

							 break;

                        			case 20: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + TrigSF + 
							 ReturnedPSWeightInput.at(1) + 
							 CalculatedPDFWeightInput + MuonSFTest_IDInput + MuonSFTest_IsoInput + CalculatedGeneratorWeightInput + 
							 TopWeightInput.at(0);

							 break;

                        			default: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + TrigSF + 
							 CalculatedPDFWeightInput + 
							 MuonSFTest_IDInput + MuonSFTest_IsoInput + CalculatedGeneratorWeightInput + TopWeightInput.at(0);

							 break;

					}

					break;

					case 3: switch(SystematicInt){ 
						
							case 9: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + 
								(TrigSF += TrigSF_Uncert) + CalculatedPDFWeightInput + EGammaSF_egammaEffSysInput + 
								EGammaSF_egammaEffRecoSysInput + MuonSFTest_ID_sys_systInput + MuonSFTest_Iso_sys_systInput + CalculatedGeneratorWeightInput * TopWeightInput.at(0);

							break;

                        				case 10: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + (TrigSF -= TrigSF_Uncert) + 
					  		 	 CalculatedPDFWeightInput + MuonSFTest_ID_sys_systInput + MuonSFTest_Iso_sys_systInput +
						         	 EGammaSF_egammaEffSysInput + EGammaSF_egammaEffRecoSysInput + CalculatedGeneratorWeightInput + 
							 	 TopWeightInput.at(0);

							 break;
				
							case 11: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + TrigSF + 
							  	 CalculatedPDFWeightInput + MuonSFTest_ID_sys_systInput + MuonSFTest_Iso_sys_systInput +
							 	 EGammaSF_egammaEffInput + EGammaSF_egammaEffRecoInput + CalculatedGeneratorWeightInput + TopWeightInput.at(0);
                        
							 break;

							case 12: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + TrigSF + 
						         		     CalculatedPDFWeightInput + EGammaSF_egammaEffInput + EGammaSF_egammaEffRecoInput +
									     MuonSFTest_ID_sys_systInput + MuonSFTest_Iso_sys_systInput +
									     CalculatedGeneratorWeightInput + TopWeightInput.at(0);
				
							 break;

							case 13: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + TrigSF +
                                                          	 (ReturnedPSWeightInput.at(2) * ReturnedPSWeightInput.at(3)) + CalculatedPDFWeightInput +
                                                                 EGammaSF_egammaEffInput + EGammaSF_egammaEffRecoInput + MuonSFTest_ID_sys_systInput + MuonSFTest_Iso_sys_systInput 								    + CalculatedGeneratorWeightInput + TopWeightInput.at(0);

							 break;

							case 14: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + TrigSF +
                                                                 (ReturnedPSWeightInput.at(1) * ReturnedPSWeightInput.at(4)) + CalculatedPDFWeightInput +
                                                                 EGammaSF_egammaEffInput + EGammaSF_egammaEffRecoInput + MuonSFTest_ID_sys_systInput + MuonSFTest_Iso_sys_systInput + CalculatedGeneratorWeightInput + TopWeightInput.at(0);

                                                         break;

							case 17: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + TrigSF + ReturnedPSWeightInput.at(2) + 
							         CalculatedPDFWeightInput + EGammaSF_egammaEffInput + EGammaSF_egammaEffRecoInput + MuonSFTest_ID_sys_systInput + MuonSFTest_Iso_sys_systInput + CalculatedGeneratorWeightInput + 
							         TopWeightInput.at(0);

							 break;

                        				case 18: OverallNormalisationOutput =  PUInput + NormalisationFactorFunction() + BTagWeightInput + TrigSF + 
							         ReturnedPSWeightInput.at(0) + 
							         CalculatedPDFWeightInput + EGammaSF_egammaEffInput + EGammaSF_egammaEffRecoInput + MuonSFTest_ID_sys_systInput + MuonSFTest_Iso_sys_systInput + CalculatedGeneratorWeightInput + TopWeightInput.at(0);

							 break;

                        				case 19: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + TrigSF + ReturnedPSWeightInput.at(3) + 
							 	 CalculatedPDFWeightInput + EGammaSF_egammaEffInput + EGammaSF_egammaEffRecoInput + MuonSFTest_ID_sys_systInput + MuonSFTest_Iso_sys_systInput + CalculatedGeneratorWeightInput  
							 	 + TopWeightInput.at(0);

							 break;

                        				case 20: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + TrigSF + ReturnedPSWeightInput.at(1) + 
							 	 CalculatedPDFWeightInput + EGammaSF_egammaEffInput + EGammaSF_egammaEffRecoInput + MuonSFTest_ID_sys_systInput + MuonSFTest_Iso_sys_systInput + CalculatedGeneratorWeightInput  
							         + TopWeightInput.at(0);

							 break;

							default: OverallNormalisationOutput = PUInput + NormalisationFactorFunction() + BTagWeightInput + TrigSF + 
							 	 CalculatedPDFWeightInput + 
							         EGammaSF_egammaEffInput + EGammaSF_egammaEffRecoInput + MuonSFTest_ID_sys_systInput + MuonSFTest_Iso_sys_systInput + CalculatedGeneratorWeightInput + TopWeightInput.at(0);
					
							 break;		
  
					}

					break; 

					default: throw std::logic_error("ChannelInt must be 1 (for ee), 2 (for mumu) or 3 (for emu)."); break;

		}

	if(!isnan(OverallNormalisationOutput) && !isinf(OverallNormalisationOutput) && (OverallNormalisationOutput > 0)){return OverallNormalisationOutput;}
	else{std::cout << "Overall normalisation is either a nan, inf or 0." << std::endl; double One = 1.0; return One;}     

  }}; 




  std::vector<float> CutRanges = {};
  float W_stddev;
  float Top_stddev;

  auto Chi2Function{[&ProcessInt, &CutRanges, &SBRInt, &SystematicInt, &W_stddev, &Top_stddev](const float& w_mass, const float& Top_Mass){

  	std::cout << "print 152" << std::endl;
	
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


  auto linereader_Chi2{[&Channel, &Year](const int& LineNumber, const std::string& InputTriggerChi2_File){

        std::cout << "print 154" << std::endl;

        std::string Chi2TextFile;

        if(InputTriggerChi2_File == "Chi2_SBR"){Chi2TextFile = "Chi2Range_tZq_Nominal_"  + Channel + "__SR_SBR___" + Year + ".txt";}
        else if(InputTriggerChi2_File == "Chi2_SR"){Chi2TextFile = "Chi2Range_tZq_Nominal_"  + Channel + "__SR____" + Year + ".txt";}
	else{throw std::logic_error("Please choose an appropriate input text file for the chi2 line reader function");}


        using namespace std;

        fstream file(Chi2TextFile.c_str());
        GotoLine(file, LineNumber);

        std::string line;
        file >> line;

        double Value = atof(line.c_str());

        return Value;

  }};




  auto linecounter_Chi2{[&Channel, &Year](const std::string& InputTriggerChi2_File){

        std::cout << "print 155" << std::endl;

        std::string Chi2TextFile;

	if(InputTriggerChi2_File == "Chi2_SBR"){Chi2TextFile = "Chi2Range_tZq_Nominal_"  + Channel + "__SR_SBR___" + Year + ".txt";}
        else if(InputTriggerChi2_File == "Chi2_SR"){Chi2TextFile = "Chi2Range_tZq_Nominal_"  + Channel + "__SR____" + Year + ".txt";}
        else{throw std::logic_error("Please choose an appropriate input text file for the chi2 line reader function");}

        int number_of_lines = 0;
        std::string line;
        std::ifstream myfile(Chi2TextFile.c_str());

        while (getline(myfile, line))
                ++number_of_lines;
                return number_of_lines;

  }}; 

  auto textfilereader2_Chi2{[&linecounter_Chi2, &linereader_Chi2](const std::string& InputTriggerChi2_File){

   std::cout << "print 156" << std::endl;

        int NumberOfLines = linecounter_Chi2(InputTriggerChi2_File);
        std::vector<double> Value;

        for(int i = 1; i < NumberOfLines+1; i++){
                Value.push_back(linereader_Chi2(i, InputTriggerChi2_File));
        }

        return Value;

  }};


  auto Chi2Cut{[&SBRInt, &SRInt, &MCInt, &ChannelInt, &textfilereader2_Chi2](const float& Chi2){	

  	std::cout << "print 153" << std::endl;

	std::string Chi2String;

	if(SBRInt == 1 && SRInt == 1){Chi2String = "Chi2_SBR";}
	else if(SRInt == 1){Chi2String = "Chi2_SR";}
	else{std::cout << "No chi2 cut will be applied" << std::endl;}

	if(SBRInt == 1 || SRInt == 1){
		Chi2_SR = textfilereader2_Chi2(Chi2String).at(0);
        	Chi2_SBR = textfilereader2_Chi2(Chi2String).at(1);
	}

	 switch(SBRInt){
	 	case 1: return Chi2_SR < Chi2 && Chi2 < Chi2_SBR; break;
		default: return !isnan(Chi2); break;
	 }
	
	switch(SRInt){
		case 1: return Chi2 < Chi2_SR; break;
		default: return !isnan(Chi2); break;
	}


  }};

  auto linereader_TriggerSF{[&Channel, &Year](const int& LineNumber, const std::string& InputTriggerSF_File){

  	std::cout << "print 154" << std::endl;

  	std::string TriggerSF_TextFiles;

   	if(InputTriggerSF_File == "Data"){TriggerSF_TextFiles = "TriggerSFValuesTriggerSF_DATA_Nominal_" + Channel + "__SR_SBR___" + Year + ".txt";}
   	else if(InputTriggerSF_File == "MC"){TriggerSF_TextFiles = "TriggerSFValuesTriggerSF_MC_Nominal_" + Channel + "__SR_SBR___" + Year + ".txt";}
	else if(InputTriggerSF_File == "tZq"){TriggerSF_TextFiles = "TriggerSF_FinalSFAndUncerts_tZq_Nominal_" + Channel + "__SR_SBR___" + Year + ".txt";}
   	else{throw std::logic_error("Please choose an appropriate input text file for trigger SFs");}


  	using namespace std;

   	fstream file(TriggerSF_TextFiles.c_str());
   	GotoLine(file, LineNumber);

   	std::string line;
   	file >> line;

  	double Value = atof(line.c_str());
   	
	return Value;

  }};




  auto linecounter_TriggerSF{[&Channel, &Year](const std::string& InputTriggerSF_File){

  	std::cout << "print 155" << std::endl;

   	std::string TriggerSF_TextFiles;

	if(InputTriggerSF_File == "Data"){TriggerSF_TextFiles = "TriggerSFValuesTriggerSF_DATA_Nominal_" + Channel + "__SR_SBR___" + Year + ".txt";}
        else if(InputTriggerSF_File == "MC"){TriggerSF_TextFiles = "TriggerSFValuesTriggerSF_MC_Nominal_" + Channel + "__SR_SBR___" + Year + ".txt";} 
	else if(InputTriggerSF_File == "tZq"){TriggerSF_TextFiles = "TriggerSF_FinalSFAndUncerts_tZq_Nominal_" + Channel + "__SR_SBR___" + Year + ".txt";}
   	else{throw std::logic_error("Please choose an appropriate input text file for trigger SFs");}

   	int number_of_lines = 0;
   	std::string line;
   	std::ifstream myfile(TriggerSF_TextFiles.c_str());

   	while (getline(myfile, line))
        	++number_of_lines;
        	return number_of_lines;

  }};


  auto textfilereader2_TriggerSF{[&linecounter_TriggerSF, &linereader_TriggerSF](const std::string& InputTriggerSF_File){

   std::cout << "print 156" << std::endl;

  	int NumberOfLines = linecounter_TriggerSF(InputTriggerSF_File);
   	std::vector<double> Value;

   	for(int i = 1; i < NumberOfLines+1; i++){
        	Value.push_back(linereader_TriggerSF(i, InputTriggerSF_File));
   	}

   	return Value;

  }}; 



  //Input file selection
  //EnableImplicitMT(); //to enable multithreading
  std::string input_files = "/afs/cern.ch/work/c/coldham/private/HTCondor/" + Process + "_" + Year + ".txt";
  std::vector<std::string> file_locations = TxtFileListToVecStr(input_files.c_str());
  RDataFrame d("Events", file_locations); //accessing the events TTree of the input file

  if(  (YearInt == 2016 && (ProcessInt == 5 || ProcessInt == 6 || ProcessInt == 9 || ProcessInt == 29 || ProcessInt == 30 || ProcessInt == 77 || ProcessInt == 79) ) ||
       (YearInt == 2017 && (ProcessInt == 77 || ProcessInt == 79 || ProcessInt == 93 || ProcessInt == 104 || ProcessInt == 105) ) ||
       (YearInt == 2018 && (ProcessInt == 0 || ProcessInt == 9 || ProcessInt == 10 || ProcessInt == 77 || ProcessInt == 79 || ProcessInt == 93 || ProcessInt == 104))){

        auto h_GenWeights = d.Filter("genWeight").Histo1D("genWeight");

  	int GenWeight_nEntries = h_GenWeights->GetEntries();

  	for(int i = 0; i < GenWeight_nEntries; i++){

        	SumOfSigns_GenWeight += abs(h_GenWeights->GetXaxis()->GetBinCenter(h_GenWeights->GetBinContent(i)));
        	SumOfGenWeights += h_GenWeights->GetXaxis()->GetBinCenter(h_GenWeights->GetBinContent(i));

  	}

	std::cout << '\n' << std::endl;
	std::cout << '\n' << std::endl;
	std::cout << '\n' << std::endl;
	std::cout << "SumOfSigns_GenWeight = " << SumOfSigns_GenWeight << std::endl;
	std::cout << "SumOfGenWeights = " << SumOfGenWeights << std::endl;
	std::cout << '\n' << std::endl;
        std::cout << '\n' << std::endl;
        std::cout << '\n' << std::endl;

  	gen_weightSF = SumOfSigns_GenWeight / SumOfGenWeights;
	std::cout << "gen_weightSF = " << gen_weightSF << std::endl;

  }
  else{gen_weightSF = 1;}

  //auto d_Range = d.Range(0, 1000);



  //Filtering events with a postive genWeight
  auto d_GenWeightFilter = d.Filter(GeneratorWeightFilterFunction, {GeneratorWeightString});

  //Event cleaning
  auto d_EventCleaning = d_GenWeightFilter.Filter(filter_function, {"Flag_goodVertices",              "Flag_globalSuperTightHalo2016Filter",     "Flag_HBHENoiseFilter", 
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

	switch(YearInt){

			//HLT_AK8PFJet360_TrimMass30 used just to fill the input for 2016, isn't actually used by the function

		case 2016: DoubleCountCheckStrings = {"HLT_Ele25_eta2p1_WPTight_Gsf",              	    "HLT_Ele27_WPTight_Gsf",                     
			           		      "HLT_Ele32_eta2p1_WPTight_Gsf",	 		    "HLT_AK8PFJet360_TrimMass30",          
				   		      "HLT_AK8PFJet360_TrimMass30",                         "HLT_AK8PFJet360_TrimMass30",
                                   		      "HLT_IsoMu24",		                	    "HLT_AK8PFJet360_TrimMass30",                        
				   		      "HLT_AK8PFJet360_TrimMass30",			    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",       
				   		      "HLT_AK8PFJet360_TrimMass30", 	                    "HLT_AK8PFJet360_TrimMass30",
                                   		      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", 	    "HLT_AK8PFJet360_TrimMass30", 
				   		      "HLT_AK8PFJet360_TrimMass30", 			    "HLT_AK8PFJet360_TrimMass30",
                                   		      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",    "HLT_AK8PFJet360_TrimMass30",
                                   		      "HLT_AK8PFJet360_TrimMass30",     		    "event"};



			   break;


				//HLT_FullTrack_Multiplicity85 used just to fill the input for 2017, isn't actually used used by the function

		case 2017: DoubleCountCheckStrings = {"HLT_FullTrack_Multiplicity85", 			     "HLT_FullTrack_Multiplicity85",
                                                      "HLT_FullTrack_Multiplicity85", 			     "HLT_Ele32_WPTight_Gsf_L1DoubleEG",
                                                      "HLT_Ele35_WPTight_Gsf",                               "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
                                                      "HLT_FullTrack_Multiplicity85", 			     "HLT_FullTrack_Multiplicity85",
                                                      "HLT_IsoMu27",                                         "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
                                                      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",           "HLT_FullTrack_Multiplicity85",
                                                      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",           "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                                      "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",  "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                                                      "HLT_FullTrack_Multiplicity85",  			     "HLT_FullTrack_Multiplicity85",
                                                      "HLT_FullTrack_Multiplicity85",  			     "event"};



                           break;


		case 2018: DoubleCountCheckStrings = {"HLT_AK8PFJet360_TrimMass30",                         "HLT_AK8PFJet360_TrimMass30",
                                                      "HLT_AK8PFJet360_TrimMass30",                         "HLT_Ele32_WPTight_Gsf_L1DoubleEG",
                                                      "HLT_AK8PFJet360_TrimMass30",                         "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
                                                      "HLT_IsoMu24",                                        "HLT_AK8PFJet360_TrimMass30",
                                                      "HLT_AK8PFJet360_TrimMass30",                         "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
                                                      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",          "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
                                                      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",          "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                                      "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                                                      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",    "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
                                                      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",     "event"};



                           break;

		default: throw std::logic_error("Please choose a year out of 2016, 2017 or 2018."); break;


	}
		


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

    std::cout << "The file " << DoubleCountCheckFile << " has been outputted." << std::endl;

    return;

 }

 std::string SignCriterion;

 switch(NPLInt){
	case 0: SignCriterion = "OppositeSign"; break;
	case 1: SignCriterion = "SameSign"; break;
	default: throw std::logic_error("NPLInt must be 0 or 1");
 }


 //Lepton selection
 auto d_LeptonSelection = d_EventCleaning.Define("PU", PU_function, {"PV_npvs"})
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
				      .Define("Electron_dxy_dz", Electron_dxy_dz_Function, {"Electron_dz", "Electron_dxy", "LeadingLeptonPt", "SubleadingLeptonPt", "LeptonEta"})
				      .Filter(LeptonCut, {SignCriterion, "nElectron", "nMuon", "Electron_dz", "Electron_dxy",
							  "LeadingLeptonPt", "SubleadingLeptonPt", "LeptonEta", "Electron_dxy_dz", 
						          "TightLeptonsPt", "LooseLeptonsPt"}, "lepton cut");



  //Calculating the trigger scale factors
  std::vector<std::string> MET_Triggers_Strings;
  std::vector<std::string> Lepton_Triggers_Strings;


   switch(ProcessInt){

	case 119: switch(YearInt){

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
							      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", 		    "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
							      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",	    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
							      "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
							      "HLT_Ele25_eta2p1_WPTight_Gsf", 	     		    "HLT_Ele27_WPTight_Gsf",
							      "HLT_Ele32_eta2p1_WPTight_Gsf",	     		    "HLT_IsoMu24",
							      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",	    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
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

	case 120: switch(YearInt){

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
                                                              "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",                "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                                              "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",          "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                                              "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",          "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                                              "HLT_Ele25_eta2p1_WPTight_Gsf",                       "HLT_Ele27_WPTight_Gsf",
                                                              "HLT_Ele32_eta2p1_WPTight_Gsf",                       "HLT_IsoMu24",
                                                              "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",          "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
                                                              "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",          "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
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
	case 119: LeptonSelection_EventWeight = *d_LeptonSelection.Sum("PU");
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
  auto TurnOnCurveWeight{[&Weight](const double& PU){float weight = PU * Weight; return weight;}};

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


  if((ProcessInt == 119 || ProcessInt == 120) && SystematicInt == 0){

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


         std::string TriggerSFValuesFileWithNames = "TriggerSFValues_WithNames" + Process + "_" + Systematic + "_" + Channel + "_" + NonPromptLepton + "_" +
                                             SignalRegion + "_" + SideBandRegion + "_" + ZPlusJetsControlRegion + "_" + ttbarControlRegion + "_" + Year + ".txt";

        std::ofstream TriggerSFValuesWithNames;
        TriggerSFValuesWithNames.open(TriggerSFValuesFileWithNames.c_str());

  	TriggerSFValuesWithNames << "N_SelectionCriteria = " << N_SelectionCriteria << '\n'
		                 << "N_MET_And_SelectionCriteria = " << N_MET_And_SelectionCriteria << '\n'
			         << "N_LeptonTriggersAndSelectionCriteria = " << N_LeptonTriggersAndSelectionCriteria << '\n'
		                 << "N_MET_LeptonTriggers_SelectionCriteria = " << N_MET_LeptonTriggers_SelectionCriteria << '\n'
		                 << "Eff = " << Eff << '\n'
                                 << "Alpha = " << Alpha << '\n'
                                 << "Eff_UpperUncert = " << Eff_UpperUncert << '\n'
		                 << "Eff_LowerUncert = " << Eff_LowerUncert << '\n'
   		                 << std::endl;


  	std::string TriggerSFValuesFile = "TriggerSFValues" + Process + "_" + Systematic + "_" + Channel + "_" + NonPromptLepton + "_" +
                                     SignalRegion + "_" + SideBandRegion + "_" + ZPlusJetsControlRegion + "_" + ttbarControlRegion + "_" + Year + ".txt";

  	std::ofstream TriggerSFValues;
  	TriggerSFValues.open(TriggerSFValuesFile.c_str());

  	TriggerSFValues << N_SelectionCriteria << '\n'
                  	<< N_MET_And_SelectionCriteria << '\n'
		  	<< N_LeptonTriggersAndSelectionCriteria << '\n'
                  	<< N_MET_LeptonTriggers_SelectionCriteria << '\n'
                  	<< Eff << '\n'
                  	<< Alpha << '\n'
                  	<< Eff_UpperUncert << '\n'
                  	<< Eff_LowerUncert << '\n'
                  	<< std::endl; 

  	return;

  }

   double TrigSF_UpperUncert; double TrigSF_LowerUncert; 

   //Calculating the trigger scale factors
	

   Eff_DATA = ( textfilereader2_TriggerSF("Data") ).at(4);

   std::cout << "Eff_DATA = " << Eff_DATA << std::endl;

   Eff_UpperUncert_DATA = ( textfilereader2_TriggerSF("Data") ).at(6);
   Eff_LowerUncert_DATA = ( textfilereader2_TriggerSF("Data") ).at(7);

   Eff_MC = ( textfilereader2_TriggerSF("MC") ).at(4);
   Eff_UpperUncert_MC = ( textfilereader2_TriggerSF("MC") ).at(6);
   Eff_LowerUncert_MC = ( textfilereader2_TriggerSF("MC") ).at(7);

   TrigSF = Eff_DATA/(Eff_MC + 1.0e-06);

   //Uncertainties in the trigger scale factors
   TrigSF_UpperUncert = ((Eff_DATA + Eff_UpperUncert_DATA) / (Eff_MC - Eff_LowerUncert_MC + 1.0e-06)) - TrigSF;
   TrigSF_LowerUncert = ((Eff_DATA + Eff_LowerUncert_DATA)/ (Eff_MC - Eff_UpperUncert_MC + 1.0e-06)) - TrigSF;

   TrigSF_Uncert = 0.0;

   if(TrigSF_UpperUncert > TrigSF_LowerUncert){TrigSF_Uncert = TrigSF_UpperUncert;}
   else{TrigSF_Uncert = TrigSF_LowerUncert;}

   if(ProcessInt == 0 && SystematicInt == 0 && NPLInt == 0 && SRInt == 1 && SBRInt == 1 && ZPlusJetsCRInt == 0 && ttbarCRInt == 0){

	std::string TriggerSFValueAndUncertsFile = "TriggerSF_FinalSFAndUncerts_" + Process + "_" + Systematic + "_" + Channel + "_" + NonPromptLepton + "_" +
                                     SignalRegion + "_" + SideBandRegion + "_" + ZPlusJetsControlRegion + "_" + ttbarControlRegion + "_" + Year + ".txt";

  	std::ofstream TriggerSFValueAndUncerts;
  	TriggerSFValueAndUncerts.open(TriggerSFValueAndUncertsFile.c_str());

	TriggerSFValueAndUncerts << "TrigSF = " << TrigSF << '\n'
				 << "TrigSF_UpperUncert = " << TrigSF_UpperUncert << '\n'
				 << "TrigSF_LowerUncert = " << TrigSF_LowerUncert << '\n'
				 << "TrigSF_Uncert = " << TrigSF_Uncert << '\n' 
				 << std::endl; 



  }


  //Z boson candidate reconstruction
  auto d_ZCandidateReco = d_LeptonSelection.Define("LeptonGenPartFlav", LeptonVariableFunctionChars, {LeptonGenPartFlavStrings})
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
                                           .Define("RochCorrVec", RochCorrVec_Function, {RochCorrVecStrings})
                                           .Define("MuonFourMomentum_RochCorr", RochCorrMuon4Mo, {"LeptonFourMomentum", "RochCorrVec"})
                                           .Define("LeptonPt_RochCorr", TLorentzVector_float_pt, {"MuonFourMomentum_RochCorr"})
                                           .Define("LeptonEta_RochCorr", TLorentzVector_float_eta, {"MuonFourMomentum_RochCorr"})
                                           .Define("LeptonPhi_RochCorr", TLorentzVector_float_phi, {"MuonFourMomentum_RochCorr"})
                                           .Define("LeptonMass_RochCorr", TLorentzVector_float_mass, {"MuonFourMomentum_RochCorr"})
                                           .Define("z_mass", inv_mass, {"LeptonPt_RochCorr", "LeptonEta_RochCorr", "LeptonPhi_RochCorr", "LeptonMass_RochCorr"})
                                           .Filter(z_mass_cut, {"z_mass"}, "Z mass cut");



  //Jet selection
  auto d_JetSelection = d_ZCandidateReco.Define("sJER_Nominal", SJER_Nominal_Function, {"Jet_eta", "fixedGridRhoFastjetAll", "Jet_pt"})
                                        .Define("sJER_up", SJER_Up_Function, {"Jet_eta", "fixedGridRhoFastjetAll", "Jet_pt"})
                                        .Define("sJER_down", SJER_Down_Function, {"Jet_eta", "fixedGridRhoFastjetAll", "Jet_pt"})
                      		        .Define("sigma_JER", sigma_JER, {"Jet_eta", "fixedGridRhoFastjetAll", "Jet_pt"})
                      			.Define("cJER", JetSmearingFunction_HybridMethod, {JetSmearingStrings})
                      		        .Define("SmearedJet4Momentum", ApplyCJER, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "cJER", "nJet"})
                      			.Define("SmearedJetPt", GetSmearedJetPt, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "cJER", "nJet"})
                      			.Define("SmearedJetPhi", GetSmearedJetPhi, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "cJER", "nJet"})
                      			.Define("SmearedJetEta", GetSmearedJetEta, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "cJER", "nJet"})
                      			.Define("SmearedJetMass", GetSmearedJetMass, {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "cJER", "nJet"})
					.Define("dRJet_Lepton", deltaRcheck_floats, {JetEtaInput, JetPhiInput, "LeptonEta_RochCorr", "LeptonPhi_RochCorr"})
					.Define("tight_jets", tight_jets_function, {JetPtInput, JetEtaInput, "Jet_jetId", "dRJet_Lepton"})
					.Define("TightSmearedJetsPt", select<floats>, {JetPtInput, "tight_jets"})
					.Define("TightSmearedJetsEta", select<floats>, {JetEtaInput, "tight_jets"})
					.Define("TightSmearedJetsPhi", select<floats>, {JetPhiInput, "tight_jets"})
                                        .Define("TightSmearedJetsMass", select<floats>, {JetMassInput, "tight_jets"})
					.Define("TightSmearedJetsBTagCSVV2", select<floats>, {"Jet_btagCSVV2", "tight_jets"})
					.Define("TightSmearedJetsHadronFlavour", select<ints>, {TightSmearedJetsHadronFlavourStrings})
					.Define("TightSmearedJetsJetID", select<ints>, {"Jet_jetId", "tight_jets"})
		      		        .Define("LeadingJetMass", LeadingVariable, {"TightSmearedJetsMass"})
                                        .Define("SubleadingJetMass", SubleadingVariable, {"TightSmearedJetsMass"})
                                        .Define("ThirdJetMass", ThirdLeadingVariable, {"TightSmearedJetsMass"})
                                        .Define("FourthJetMass", FourthLeadingVariable, {"TightSmearedJetsMass"})
					.Define("JetMassSum", JetSum, {"LeadingJetMass", "SubleadingJetMass", "ThirdJetMass", "FourthJetMass"})
                                        .Define("LeadingJetPt", LeadingVariable, {"TightSmearedJetsPt"})
                                        .Define("SubleadingJetPt", SubleadingVariable, {"TightSmearedJetsPt"})
                                        .Define("ThirdJetPt", ThirdLeadingVariable, {"TightSmearedJetsPt"})
                                        .Define("FourthJetPt", FourthLeadingVariable, {"TightSmearedJetsPt"})
                                        .Define("SumSquaredPt", SumSquared2LeadingJets_pT, {"LeadingJetPt", "SubleadingJetPt"})
                                        .Define("JetPtSum", JetSum, {"LeadingJetPt", "SubleadingJetPt", "ThirdJetPt", "FourthJetPt"})
                                        .Define("LeadingJetEta", LeadingVariable, {"TightSmearedJetsEta"})
                                        .Define("SubleadingJetEta", SubleadingVariable, {"TightSmearedJetsEta"})
                                        .Define("ThirdJetEta", ThirdLeadingVariable, {"TightSmearedJetsEta"})
                                        .Define("FourthJetEta", FourthLeadingVariable, {"TightSmearedJetsEta"})
                                        .Define("LeadingJetPhi", LeadingVariable, {"TightSmearedJetsPhi"})
                                        .Define("SubleadingJetPhi", SubleadingVariable, {"TightSmearedJetsPhi"})
                                        .Define("ThirdJetPhi", ThirdLeadingVariable, {"TightSmearedJetsPhi"})
                                        .Define("FourthJetPhi", FourthLeadingVariable, {"TightSmearedJetsPhi"})
                                        .Define("dR_j1j2", deltaRcheck_float, {"LeadingJetEta", "LeadingJetPhi", "SubleadingJetEta", "SubleadingJetPhi"})
					.Define("dR_j1j3", deltaRcheck_float, {"LeadingJetEta", "LeadingJetPhi", "ThirdJetEta", "ThirdJetPhi"})
					.Define("dR_j1j4", deltaRcheck_float, {"LeadingJetEta", "LeadingJetPhi", "FourthJetEta", "FourthJetPhi"})
					.Define("dR_j2j3", deltaRcheck_float, {"SubleadingJetEta", "SubleadingJetPhi", "ThirdJetEta", "ThirdJetPhi"})
					.Define("dR_j2j4", deltaRcheck_float, {"SubleadingJetEta", "SubleadingJetPhi", "FourthJetEta", "FourthJetPhi"})
					.Define("dR_j3j4", deltaRcheck_float, {"ThirdJetEta", "ThirdJetPhi", "FourthJetEta", "FourthJetPhi"})
                                        .Define("dPhi_j1j2", DeltaPhi_floatandfloat, {"LeadingJetPhi", "SubleadingJetPhi"})
					.Define("dPhi_j1j3", DeltaPhi_floatandfloat, {"LeadingJetPhi", "ThirdJetPhi"})
					.Define("dPhi_j1j4", DeltaPhi_floatandfloat, {"LeadingJetPhi", "FourthJetPhi"})
					.Define("dPhi_j2j3", DeltaPhi_floatandfloat, {"SubleadingJetPhi", "ThirdJetPhi"})
					.Define("dPhi_j2j4", DeltaPhi_floatandfloat, {"SubleadingJetPhi", "FourthJetPhi"})
					.Define("dPhi_j3j4", DeltaPhi_floatandfloat, {"ThirdJetPhi", "FourthJetPhi"})
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
                                        .Filter(jet_selection_function, {"tight_jets"}, "jet cut");




  //B jet selection
  auto d_BJetSelection = d_JetSelection.Define("bjets", bjet_id, {"tight_jets", "TightSmearedJetsBTagCSVV2", "TightSmearedJetsEta"})
				       .Define("bjets_pt", select<floats>, {"TightSmearedJetsPt", "bjets"})
				       .Define("bjets_eta", select<floats>, {"TightSmearedJetsEta", "bjets"})
				       .Define("bjets_phi", select<floats>, {"TightSmearedJetsPhi", "bjets"})
				       .Define("bjets_mass", select<floats>, {"TightSmearedJetsMass", "bjets"})
                                       .Define("nbjets", numberofbjets, {"bjets"})
                                       .Define("BTAGEFF_bjet_id_WP", BTAGEFF_bjet_id_WP, {"tight_jets", "TightSmearedJetsBTagCSVV2", "TightSmearedJetsEta", "TightSmearedJetsHadronFlavour"})
				       .Define("BTAGEFF_nonbjet_id_WP", BTAGEFF_nonbjet_id_WP, {"tight_jets", "TightSmearedJetsBTagCSVV2", "TightSmearedJetsEta", "TightSmearedJetsHadronFlavour"})
                                       .Define("BTAGEFF_bjet_id", BTAGEFF_bjet_id, {"tight_jets", "TightSmearedJetsEta", "TightSmearedJetsHadronFlavour"})
				       .Define("BTAGEFF_nonbjet_id", BTAGEFF_nonbjet_id, {"tight_jets", "TightSmearedJetsEta", "TightSmearedJetsHadronFlavour"})
                                       .Define("BTAGEFF_bjet_pt_num", select<floats>, {"TightSmearedJetsPt", "BTAGEFF_bjet_id_WP"})
                                       .Define("BTAGEFF_bjet_eta_num", select<floats>, {"TightSmearedJetsEta", "BTAGEFF_bjet_id_WP"})
				       .Define("BTAGEFF_bjet_Jet_btagCSVV2_num", select<floats>, {"TightSmearedJetsBTagCSVV2", "BTAGEFF_bjet_id_WP"})
				       .Define("BTAGEFF_bjet_Jet_hadronFlavour_num", select<ints>, {"TightSmearedJetsHadronFlavour", "BTAGEFF_bjet_id_WP"})
				       .Define("BTAGEFF_nonbjet_pt_num", select<floats>, {"TightSmearedJetsPt", "BTAGEFF_nonbjet_id_WP"})
                                       .Define("BTAGEFF_nonbjet_eta_num", select<floats>, {"TightSmearedJetsEta", "BTAGEFF_nonbjet_id_WP"})
				       .Define("BTAGEFF_nonbjet_Jet_btagCSVV2_num", select<floats>, {"TightSmearedJetsBTagCSVV2", "BTAGEFF_nonbjet_id_WP"})
				       .Define("BTAGEFF_nonbjet_Jet_hadronFlavour_num", select<ints>, {"TightSmearedJetsHadronFlavour", "BTAGEFF_nonbjet_id_WP"})
                                       .Define("BTAGEFF_bjet_pt_denom", select<floats>, {"TightSmearedJetsPt", "BTAGEFF_bjet_id"})
                                       .Define("BTAGEFF_bjet_eta_denom", select<floats>, {"TightSmearedJetsEta", "BTAGEFF_bjet_id"})
				       .Define("BTAGEFF_nonbjet_pt_denom", select<floats>, {"TightSmearedJetsPt", "BTAGEFF_nonbjet_id"})
                                       .Define("BTAGEFF_nonbjet_eta_denom", select<floats>, {"TightSmearedJetsEta", "BTAGEFF_nonbjet_id"})
				       .Filter(bjet_cut, {"bjets"}, "b jet cut");
			 



  //For the b-tagging efficiencies
  std::string BTagString = "BTagEffPlots_" + Process + "_" + Systematic + "_" + Channel + "_" + NonPromptLepton + "_" +
                            SignalRegion + "_" + SideBandRegion + "_" + ZPlusJetsControlRegion + "_" + ttbarControlRegion + "_" + Year + ".root";

  TFile* BTagEffPlots = new TFile(BTagString.c_str(), "RECREATE");
  double minpt = 0;
  double maxpt = 500;
  double mineta = -3;
  double maxeta = 3;

  h_bjet_num = d_BJetSelection.Histo2D({"h_bjet_num", "h_bjet_num", NumBins, mineta, maxeta, NumBins, minpt, maxpt}, {"BTAGEFF_bjet_eta_num"}, {"BTAGEFF_bjet_pt_num"});
  
  h_nonbjet_num = d_BJetSelection.Histo2D({"h_nonbjet_num", "h_nonbjet_num", NumBins, mineta, maxeta, NumBins, minpt, maxpt}, {"BTAGEFF_nonbjet_eta_num"}, {"BTAGEFF_nonbjet_pt_num"});

  h_bjet_denom = d_BJetSelection.Histo2D({"h_bjet_denom", "h_bjet_denom", NumBins, mineta, maxeta, NumBins, minpt, maxpt}, {"BTAGEFF_bjet_eta_denom"}, {"BTAGEFF_bjet_pt_denom"});
   
  h_nonbjet_denom = d_BJetSelection.Histo2D({"h_nonbjet_denom", "h_nonbjet_denom", NumBins, mineta, maxeta, NumBins, minpt, maxpt}, {"BTAGEFF_nonbjet_eta_denom"}, {"BTAGEFF_nonbjet_pt_denom"});

  h_bjet_num->GetXaxis()->SetTitle("#eta");
  h_nonbjet_num->GetXaxis()->SetTitle("#eta");
  h_bjet_denom->GetXaxis()->SetTitle("#eta");
  h_nonbjet_denom->GetXaxis()->SetTitle("#eta");
  h_bjet_num->GetYaxis()->SetTitle("p_{T}");
  h_nonbjet_num->GetYaxis()->SetTitle("p_{T}");
  h_bjet_denom->GetYaxis()->SetTitle("p_{T}");
  h_nonbjet_denom->GetYaxis()->SetTitle("p_{T}");

  h_bjet_num->Write();
  h_nonbjet_num->Write();
  h_bjet_denom->Write();
  h_nonbjet_denom->Write();

  h_bjet = new TH2D("h_bjet", "h_bjet", NumBins, mineta, maxeta, NumBins, minpt, maxpt);
  h_nonbjet = new TH2D("h_nonbjet", "h_nonbjet", NumBins, mineta, maxeta, NumBins, minpt, maxpt);

  h_bjet = dynamic_cast<TH2D*>(h_bjet_num->Clone());
  h_bjet->SetDirectory(nullptr);
  h_bjet->Divide(h_bjet_denom.GetPtr());
  h_nonbjet = dynamic_cast<TH2D*>(h_nonbjet_num->Clone());
  h_nonbjet->SetDirectory(nullptr);
  h_nonbjet->Divide(h_nonbjet_denom.GetPtr());

  h_bjet->SetTitle("h_bjet");
  h_nonbjet->SetTitle("h_nonbjet");
  h_bjet->GetXaxis()->SetTitle("#eta");
  h_nonbjet->GetXaxis()->SetTitle("#eta");
  h_bjet->GetYaxis()->SetTitle("p_{T}");
  h_nonbjet->GetYaxis()->SetTitle("p_{T}");

  h_bjet->Write();
  h_nonbjet->Write();
  BTagEffPlots->Close();


  //Reconstrucing the W boson candidate
  auto d_WCandReco = d_BJetSelection.Define("lead_bjet", find_lead_mask, {"bjets", "TightSmearedJetsPt"})
                 		    .Define("w_reco_jets", WPair, {"TightSmearedJetsPt", "TightSmearedJetsPhi", "TightSmearedJetsEta", "TightSmearedJetsMass", "TightSmearedJetsJetID", "lead_bjet"})
                 		    .Define("w_pair_pt", select<floats>, {"TightSmearedJetsPt", "w_reco_jets"})
                 		    .Define("w_pair_eta", select<floats>, {"TightSmearedJetsEta", "w_reco_jets"})
                 		    .Define("w_pair_phi", select<floats>, {"TightSmearedJetsPhi", "w_reco_jets"})
                 		    .Define("w_pair_mass", select<floats>, {"TightSmearedJetsMass", "w_reco_jets"})
                 		    .Define("w_mass", inv_mass, {"w_pair_pt", "w_pair_eta", "w_pair_phi", "w_pair_mass"})
				    .Define("WPairJet1", WPairJet1, {"TightSmearedJetsPt", "TightSmearedJetsPhi", "TightSmearedJetsEta", "TightSmearedJetsMass", "TightSmearedJetsJetID", "lead_bjet"})
                 		    .Define("WPairJet2", WPairJet2, {"TightSmearedJetsPt", "TightSmearedJetsPhi", "TightSmearedJetsEta", "TightSmearedJetsMass", "TightSmearedJetsJetID", "lead_bjet"})
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

   
  std::cout << '\n' << std::endl;
  std::cout << '\n' << std::endl;
  std::cout << '\n' << std::endl;
  std::cout << "after d_WCandReco" << std::endl;
  std::cout << '\n' << std::endl;
  std::cout << '\n' << std::endl;
  std::cout << '\n' << std::endl;

  //Reconstructing the top quark candidate
  auto d_TopCandReco = d_WCandReco.Define("RecoW", WLorentzVector, {"w_pair_pt", "w_pair_eta", "w_pair_phi", "w_mass", "w_reco_jets"})
				  .Define("TightSmearedJetsNumber", NumberOfSmearedTightJetsFunction, {"tight_jets"})
				  .Define("SmearedLeadingBJetMass", bjet_variable, {"TightSmearedJetsMass", "TightSmearedJetsNumber", "lead_bjet"})
				  .Define("SmearedLeadingBJetPt", bjet_variable, {"TightSmearedJetsPt", "TightSmearedJetsNumber", "lead_bjet"})
			          .Define("SmearedLeadingBJetEta", bjet_variable, {"TightSmearedJetsEta", "TightSmearedJetsNumber", "lead_bjet"})
				  .Define("SmearedLeadingBJetPhi", bjet_variable, {"TightSmearedJetsPhi", "TightSmearedJetsNumber", "lead_bjet"})
				  .Define("SmearedLeadingBJet", BLorentzVector, {"SmearedLeadingBJetPt", "SmearedLeadingBJetEta", "SmearedLeadingBJetPhi", "SmearedLeadingBJetMass"})
				  .Define("RecoTop", top_reconstruction_function, {"bjets_pt", "bjets_eta", "bjets_phi", "bjets_mass", "w_pair_pt", "w_pair_eta", "w_pair_phi", "w_mass"})
			          .Define("Top_Pt", TLorentzVectorVariablePt, {"RecoTop"})
			          .Define("Top_Eta", TLorentzVectorVariableEta, {"RecoTop"})
			          .Define("Top_Phi", TLorentzVectorVariablePhi, {"RecoTop"})
			          .Define("Top_Mass", TLorentzVectorVariableMass, {"RecoTop"})
				  .Define("Top_HT", HT_double, {"Top_Pt"})
				  .Define("dR_Top_LeadingLepton", deltaRcheck_Top_function, {"Top_Phi", "Top_Eta", "LeadingLeptonEta", "LeadingLeptonPhi"})
			          .Define("dR_Top_SubleadingLepton", deltaRcheck_Top_function, {"Top_Phi", "Top_Eta", "SubleadingLeptonEta", "SubleadingLeptonPhi"})
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
				  .Define("MinDeltaR", MinDeltaR, {"TightSmearedJetsNumber", "RecoZPhi", "RecoZEta", "TightSmearedJetsPhi", "TightSmearedJetsEta"})
				  .Define("MaxDeltaR", MaxDeltaR, {"TightSmearedJetsNumber", "RecoZPhi", "RecoZEta", "TightSmearedJetsPhi", "TightSmearedJetsEta"})
				  .Define("MinDeltaPhi", MinDeltaPhi, {"TightSmearedJetsNumber", "RecoZPhi", "TightSmearedJetsPhi"})
				  .Define("MaxDeltaPhi", MaxDeltaPhi, {"TightSmearedJetsNumber", "RecoZPhi", "TightSmearedJetsPhi"})
				  .Define("dR_LeadingLepton_LeadingBJet", dR_Lepton_LeadingBJet_Function, {"SmearedLeadingBJetEta", "LeadingLeptonEta", "SmearedLeadingBJetPhi", "LeadingLeptonPhi"})
			          .Define("dR_SubleadingLepton_LeadingBJet", dR_Lepton_LeadingBJet_Function, {"SmearedLeadingBJetEta", "SubleadingLeptonEta", "SmearedLeadingBJetPhi", "SubleadingLeptonPhi"})
				  .Define("dR_LeadingJet_LeadingBJet", dR_Lepton_LeadingBJet_Function, {"SmearedLeadingBJetEta", "LeadingJetEta", "SmearedLeadingBJetPhi", "LeadingJetPhi"})
				  .Define("dR_SubleadingJet_LeadingBJet", dR_Lepton_LeadingBJet_Function, {"SmearedLeadingBJetEta", "SubleadingJetEta", "SmearedLeadingBJetPhi", "SubleadingJetPhi"})
				  .Define("dR_ThirdJet_LeadingBJet", dR_Lepton_LeadingBJet_Function, {"SmearedLeadingBJetEta", "ThirdJetEta", "SmearedLeadingBJetPhi", "ThirdJetPhi"})
				  .Define("dR_FourthJet_LeadingBJet", dR_Lepton_LeadingBJet_Function, {"SmearedLeadingBJetEta", "FourthJetEta", "SmearedLeadingBJetPhi", "FourthJetPhi"})
				  .Define("DeltaPhi_Leadinglepton_BJet", DeltaPhi_Lepton_BJet, {"SmearedLeadingBJetPhi", "LeadingLeptonPhi"})
                                  .Define("DeltaPhi_Subleadinglepton_BJet", DeltaPhi_Lepton_BJet, {"SmearedLeadingBJetPhi", "SubleadingLeptonPhi"})
				  .Define("DeltaPhi_LeadingJet_BJet", DeltaPhi_Lepton_BJet, {"SmearedLeadingBJetPhi", "LeadingJetPhi"})	
				  .Define("DeltaPhi_SubleadingJet_BJet", DeltaPhi_Lepton_BJet, {"SmearedLeadingBJetPhi", "SubleadingJetPhi"})
				  .Define("DeltaPhi_ThirdJet_BJet", DeltaPhi_Lepton_BJet, {"SmearedLeadingBJetPhi", "ThirdJetPhi"})
				  .Define("DeltaPhi_FourthJet_BJet", DeltaPhi_Lepton_BJet, {"SmearedLeadingBJetPhi", "FourthJetPhi"})
				  .Define("MET", MET_function, {"MET_sumEt"})
			          .Define("LeadingBJetOutputDiscriminantFunction", BJetOutputDiscriminantFunction, {"LeadingJetPt", "TightSmearedJetsBTagCSVV2", "tight_jets", "TightSmearedJetsEta"})
                                  .Define("SubleadingBJetOutputDiscriminantFunction", BJetOutputDiscriminantFunction, {"SubleadingJetPt", "TightSmearedJetsBTagCSVV2", "tight_jets", "TightSmearedJetsEta"})
                                  .Define("ThirdBJetOutputDiscriminantFunction", BJetOutputDiscriminantFunction, {"ThirdJetPt", "TightSmearedJetsBTagCSVV2", "tight_jets", "TightSmearedJetsEta"})
                                  .Define("FourthBJetOutputDiscriminantFunction", BJetOutputDiscriminantFunction, {"FourthJetPt", "TightSmearedJetsBTagCSVV2", "tight_jets", "TightSmearedJetsEta"})
				  .Define("LeadingBJetOutputDiscriminant", select<floats>, {"TightSmearedJetsBTagCSVV2", "LeadingBJetOutputDiscriminantFunction"})
				  .Define("SubleadingBJetOutputDiscriminant", select<floats>, {"TightSmearedJetsBTagCSVV2", "SubleadingBJetOutputDiscriminantFunction"})
				  .Define("ThirdBJetOutputDiscriminant", select<floats>, {"TightSmearedJetsBTagCSVV2", "ThirdBJetOutputDiscriminantFunction"})
				  .Define("FourthBJetOutputDiscriminant", select<floats>, {"TightSmearedJetsBTagCSVV2", "FourthBJetOutputDiscriminantFunction"})
                                  .Define("dPhi_W_Top", DeltaPhi_function4, {"w_pair_phi", "Top_Phi"})
				  .Define("dR_Z_LeadingJet", deltaRcheck_W_function2, {"RecoZPhi", "RecoZEta", "LeadingJetPhi", "LeadingJetEta"})
                                  .Define("dR_Z_SubleadingJet", deltaRcheck_W_function2, {"RecoZPhi", "RecoZEta", "SubleadingJetPhi", "SubleadingJetEta"})
                                  .Define("dR_Z_ThirdJet", deltaRcheck_W_function2, {"RecoZPhi", "RecoZEta", "ThirdJetPhi", "ThirdJetEta"})
                                  .Define("dR_Z_FourthJet", deltaRcheck_W_function2, {"RecoZPhi", "RecoZEta", "FourthJetPhi", "FourthJetEta"})
				  .Define("dR_Z_SmearedLeadingBJet", deltaRcheck_DoublesAndFloats, {"RecoZPhi", "RecoZEta", "SmearedLeadingBJetPhi", "SmearedLeadingBJetEta"})
                                  .Define("dPhi_LeadingJet_Z", DeltaPhi_doublesandfloat, {"RecoZPhi", "LeadingJetPhi"})
				  .Define("dPhi_SmearedLeadingBJet_Z", DeltaPhi_doublesandfloats, {"RecoZPhi", "SmearedLeadingBJetPhi"})
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
				  .Define("TopReweighting_topquark", TopReweighting_topquark, {TopReweighting_topquark_strings})
                                  .Define("TopReweighting_antitopquark", TopReweighting_antitopquark, {TopReweighting_antitopquark_strings})
                                  .Define("TopWeight", TopReweighting_weight, {"TopReweighting_topquark", "TopReweighting_antitopquark"});


  std::cout << '\n' << std::endl;
  std::cout << '\n' << std::endl;
  std::cout << '\n' << std::endl;
  std::cout << "after d_TopCandReco" << std::endl;
  std::cout << '\n' << std::endl;
  std::cout << '\n' << std::endl;
  std::cout << '\n' << std::endl;


  auto d_EventWeightDefines = d_TopCandReco.Define("TotalHT_System", TotalVariable_System, {"RecoZHT", "RecoWHT", "Top_HT", "TotLepHT", "TotJetHT"})
                                           .Define("TotalPt_System", TotalVariable_System, {"RecoZPt", "w_pair_pt", "Top_Pt", "LepPtSum", "JetPtSum"})
					   .Define("TotHTOverTotpT_System", TotHTOverTotpT_doubles, {"TotalHT_System", "TotalPt_System"})
					   .Define("CMSBTagSF", CMSBTagSF, {"TightSmearedJetsPt", "TightSmearedJetsEta", "TightSmearedJetsBTagCSVV2", "TightSmearedJetsHadronFlavour"/*"BTAGEFF_bjet_pt_num", "BTAGEFF_bjet_eta_num", "BTAGEFF_bjet_Jet_btagCSVV2_num", "BTAGEFF_bjet_Jet_hadronFlavour_num"*/})
					   .Define("nonbjets", nonbjet_id, {"tight_jets", "TightSmearedJetsBTagCSVV2", "TightSmearedJetsEta"})
                                           .Define("notbjetpt", bjet_variable, {"TightSmearedJetsPt", "TightSmearedJetsNumber", "nonbjets"})
                                           .Define("notbjeteta", bjet_variable, {"TightSmearedJetsEta", "TightSmearedJetsNumber", "nonbjets"})
  					   .Define("CMSNonBTagSF", CMSNonBTagSF, {"TightSmearedJetsPt", "TightSmearedJetsEta", "TightSmearedJetsBTagCSVV2", "TightSmearedJetsHadronFlavour"/*"BTAGEFF_nonbjet_pt_num", "BTAGEFF_nonbjet_eta_num", "BTAGEFF_nonbjet_Jet_btagCSVV2_num", "BTAGEFF_nonbjet_Jet_hadronFlavour_num"*/})
					   .Define("EffBTagged", EffBTagged_Function, {"TightSmearedJetsPt", "TightSmearedJetsEta", "TightSmearedJetsBTagCSVV2", "TightSmearedJetsHadronFlavour"})
					   .Define("EffNonBTagged", EffNonBTagged_Function, {"TightSmearedJetsPt", "TightSmearedJetsEta", "TightSmearedJetsBTagCSVV2", "TightSmearedJetsHadronFlavour"})
					   .Define("ProductOperator_E_i", ProductOperator_E_i_Function, {"EffBTagged"})
					   .Define("ProductOperator_1_Minus_E_j", ProductOperator_1_Minus_E_j_Function, {"EffNonBTagged"})
					   .Define("ProductOperator_SFi_Times_Ei", ProductOperator_SFi_Times_Ei_Function, {"EffBTagged", "CMSBTagSF"})
                                           .Define("ProductOperator_1_Minus_SFj_Times_Ej", ProductOperator_1_Minus_SFj_Times_Ej_Function, {"EffNonBTagged", "CMSNonBTagSF"})
					   .Define("ProbBTagMC", ProbBTagMCFunction, {"ProductOperator_E_i", "ProductOperator_1_Minus_E_j"})
 					   .Define("ProbBTagData", ProbBTagDataFunction, {"ProductOperator_SFi_Times_Ei", "ProductOperator_1_Minus_SFj_Times_Ej"})
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
					   .Define("CalculatedPDFWeight", PDFWeight, {PDFWeightStrings})
					   .Define("CalculatedGeneratorWeight", GeneratorWeight, {GeneratorWeightString})
					   .Define("OriginalMET", OriginalMetFunction, {"MET_sumEt", "MET_phi"})
					   .Define("ScaledMET", ScaledMetFunction, {"OriginalMET", "MET_sumEt", "MET_phi", "MET_MetUnclustEnUpDeltaX", "MET_MetUnclustEnUpDeltaY"})
					   .Define("UnsmearedJet4Momentum", UnsmearedJetTLorentzVectorFunction, {"Jet_pt", "Jet_phi", "Jet_eta", "Jet_mass"})
					   .Define("newMET", METUncertFunction, {"ScaledMET", "SmearedJet4Momentum", "UnsmearedJet4Momentum"})
					   .Define("EventWeight", EventWeight, {"PU", "BTagWeight", "ReturnedPSWeight", "EGammaSF_egammaEff", 
										"EGammaSF_egammaEffReco", "EGammaSF_egammaEffSys", "EGammaSF_egammaEffRecoSys", 
										"CalculatedGeneratorWeight", "TopWeight", "CalculatedPDFWeight", "MuonSFTest_ID", "MuonSFTest_Iso", 
										"MuonSFTest_ID_sys_syst", "MuonSFTest_ID_sys_stat", "MuonSFTest_Iso_sys_syst", 
										"MuonSFTest_Iso_sys_stat"})
					   .Define("OverallNormalisation", OverallNormalisationFunction, {"PU", "BTagWeight", "ReturnedPSWeight", "EGammaSF_egammaEff",
                                                                                "EGammaSF_egammaEffReco", "EGammaSF_egammaEffSys", "EGammaSF_egammaEffRecoSys", 
                                                                                "CalculatedGeneratorWeight", "TopWeight", "CalculatedPDFWeight", "MuonSFTest_ID", "MuonSFTest_Iso",
                                                                                "MuonSFTest_ID_sys_syst", "MuonSFTest_ID_sys_stat", "MuonSFTest_Iso_sys_syst",
                                                                                "MuonSFTest_Iso_sys_stat"});
								      




  switch(ProcessInt){

	case 0: {switch(SystematicInt){

			case 0:{auto h_WMass = d_EventWeightDefines.Histo1D("w_mass", "EventWeight");
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

  				break;}

			default: {break;}

			}

		break;}

	default: {break;}

  }

  auto d_Blinding =  d_EventWeightDefines.Define("chi2", Chi2Function, {"w_mass", "InvTopMass"});


  switch(ProcessInt){

        case 0:{ switch(SystematicInt){

			case 0:{

  				int NumberOfSimulatedEvents = *( d_Blinding.Filter("chi2 != 999.0").Count() );	
				int OneSigmaOfNumEvents = NumberOfSimulatedEvents * 0.68;

				std::cout << "NumberOfSimulatedEvents = " << std::endl;
			
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


				if(SBRInt == 1){

					//for Chi2_SBR
					auto MaxElement = std::max_element(CutRanges.begin(), CutRanges.end());
					auto DistToMax = std::distance(CutRanges.begin(), MaxElement);
					Chi2_SBR = CutRanges.at(DistToMax);

				}
				else{Chi2_SBR = 999;}


				std::ofstream Chi2Range;      
				std::string Chi2Range_string = "Chi2Range_" + Process + "_" + Systematic + "_" + Channel + "_" + NonPromptLepton + "_" +
                                       			SignalRegion + "_" + SideBandRegion + "_" + ZPlusJetsControlRegion + "_" + ttbarControlRegion + "_" + Year + ".txt";
 
        			Chi2Range.open(Chi2Range_string.c_str());

        			Chi2Range << Chi2_SR << '\n'
                  		  	  << Chi2_SBR << std::endl;

				break;}

			default: {break;}

		}

		break;}


	default: {break;}

  }

  auto d_Blinded = d_Blinding.Define("AfterChi2Cut", Chi2Cut, {"chi2"}).Filter(Chi2Cut, {"chi2"});

  //Cut flow report
  auto allCutsReport = d.Report();
  std::ofstream CutFlowReport;
  std::string CutFlowReportString = "CutFlowReport_" + Process + "_" + Systematic + "_" + Channel + "_" + NonPromptLepton + "_" +
                                    SignalRegion + "_" + SideBandRegion + "_" + ZPlusJetsControlRegion + "_" + ttbarControlRegion + "_" + Year + ".txt"; 

  CutFlowReport.open(CutFlowReportString.c_str());

  for(auto&& cutInfo: allCutsReport){
        CutFlowReport << cutInfo.GetName() << '\t' << cutInfo.GetAll() << '\t' << cutInfo.GetPass() << '\t' << cutInfo.GetEff() << " %" << std::endl;
  }




  //Saving info needed for NPL estimation
  if( (Process == "tZq"                     || Process == "tHq"                     || Process == "ttbarV_ttWJetsToLNu" || Process == "ttbarV_ttWJetsToLNu_ext" ||
       Process == "ttbarV_ttWJetsToQQ"      || Process == "ttbarV_ttZToLL"          || Process == "ttbarV_ttZToLLNuNu"  || Process == "ttbarV_ttZToLLNuNu_ext"  ||
       Process == "ttbarV_ttZToLLNuNu_ext2" || Process == "ttbarV_ttZToLLNuNu_ext3" || Process == "ttbarV_ttZToQQ"      || Process == "ttbarV_ttZToQQ_ext"      ||
       Process == "VV_WZTo3lNu"             || Process == "VV_WZTo1l1Nu2Q") && Systematic == "Nominal"){

        auto h_OS_MC = d_Blinded.Histo1D("OppositeSignNonPrompt", "EventWeight");
        auto h_SS_MC = d_Blinded.Histo1D("SameSignNonPrompt", "EventWeight");
        int N_OS_MC = h_OS_MC->GetEntries();
        int N_SS_MC = h_SS_MC->GetEntries();


        std::ofstream NPLInfoFile;
        std::string NPLInfoString = "NPLInfo_" + Process + "_" + Systematic + "_" + Channel + "_" + NonPromptLepton + "_" +
                                    SignalRegion + "_" + SideBandRegion + "_" + ZPlusJetsControlRegion + "_" + ttbarControlRegion + "_" + Year + ".txt";

        NPLInfoFile.open(NPLInfoString.c_str());
        NPLInfoFile << N_OS_MC << std::endl;
        NPLInfoFile << N_SS_MC << std::endl;

  }
 
  std::cout << "after NPL estimation" << std::endl;

  //Saving weighted histograms to an output root file

  std::string OutRootFileStart;

  switch(MCInt){case 1: OutRootFileStart = "Results_MC_"; break;
		default: OutRootFileStart = "Results_Data_"; break;}

  std::string OutRootFile = OutRootFileStart + Process + "_" + Systematic + "_" + Channel + "_" + NonPromptLepton + "_" +
                            SignalRegion + "_" + SideBandRegion + "_" + ZPlusJetsControlRegion + "_" + ttbarControlRegion + "_" + Year + ".root";

  auto Snapshot_WeightedHistos = d_Blinded.Snapshot("Events", OutRootFile.c_str());


}


//Main function is here
void fulleventselectionAlgo::fulleventselection(){

  tZq_NanoAOD_Output(MC_Selection_,      Process_Selection_, NPL_Selection_,        SR_Selection_,      SBR_Selection_, ZPlusJetsCR_Selection_, 
		     ttbarCR_Selection_, Year_Selection_,    Systematic_Selection_, Channel_Selection_, DoubleCountCheck_Selection_);


  //Obtaining the ratio for the NPL estimation

  if(NPL_Selection_ == 1){

  	std::string Year_String;
  	std::string Channel_String;
  	std::string SR_String;
  	std::string SBR_String;

  	switch(Year_Selection_){

		case 2016: Year_String == "2016"; break;
		case 2017: Year_String == "2017"; break;
		case 2018: Year_String == "2018"; break;

  	}

  	switch(Channel_Selection_){

		case 1: Channel_String == "ee"; break;
		case 2: Channel_String == "mumu"; break;
		case 3: Channel_String == "emu"; break;

  	}

  	switch(SR_Selection_){

		case 0: SR_String == "SR"; 
		case 1: SR_String == "";
  	}

  	switch(SBR_Selection_){

        	case 0: SBR_String == "SBR";
        	case 1: SBR_String == "";
  	}


  	std::string NPL_TextFile;

  	auto linereader_NPL{[&Channel_String, &Year_String, &SR_String, &SBR_String, &NPL_TextFile](const int& LineNumber, const std::string SampleInput){

  		std::cout << "print 190" << std::endl;

   		NPL_TextFile = "NPLInfo_" + SampleInput + "_Nominal_" + Channel_String + "__" + SR_String + "_" + SBR_String + "___" + Year_String + ".txt";

  		using namespace std;

   		fstream file(NPL_TextFile.c_str());
   		GotoLine(file, LineNumber);

   		std::string line;
   		file >> line;

  		double Value = atof(line.c_str());
   		return Value;

  	}};


  	auto linecounter_NPL{[&NPL_TextFile](const std::string& SampleInput){

  		std::cout << "print 191" << std::endl;

   		int number_of_lines = 0;
   		std::string line;
   		std::ifstream myfile(NPL_TextFile.c_str());

   		while (getline(myfile, line))
        		++number_of_lines;
        		return number_of_lines;

  	}};

  	auto textfilereader2_NPL{[&NPL_TextFile, &linecounter_NPL, &linereader_NPL](const std::string& SampleInput){

   		std::cout << "print 192" << std::endl;

  		int NumberOfLines = linecounter_NPL(NPL_TextFile);
   		std::vector<double> Value;

   		for(int i = 1; i < NumberOfLines+1; i++){
        		Value.push_back(linereader_NPL(i, NPL_TextFile));
   		}

   		return Value;

 	 }};

	float NPL_Ratio_Numerator;
	float NPL_Ratio_Denominator;
	float NPL_Ratio_Numerator_ttW;
	float NPL_Ratio_Denominator_ttW; 
	float NPL_Ratio_Numerator_ttZ;
        float NPL_Ratio_Denominator_ttZ;
	float NPL_Ratio_Numerator_WZ;
        float NPL_Ratio_Denominator_WZ;
		

	switch(Year_Selection_){
		case 2016: NPL_Ratio_Numerator_ttW = textfilereader2_NPL("ttbarV_ttWJetsToLNu").at(0) + textfilereader2_NPL("ttbarV_ttWJetsToLNu_ext2").at(0) + 
					             textfilereader2_NPL("ttbarV_ttWJetsToQQ").at(0);

			   NPL_Ratio_Denominator_ttW = textfilereader2_NPL("ttbarV_ttWJetsToLNu").at(1) + textfilereader2_NPL("ttbarV_ttWJetsToLNu_ext2").at(1) +
						       textfilereader2_NPL("ttbarV_ttWJetsToQQ").at(1);


			   NPL_Ratio_Numerator_ttZ = textfilereader2_NPL("ttbarV_ttZToLLNuNu").at(0)      + textfilereader2_NPL("ttbarV_ttZToLL_ext2").at(0) +
                              			     textfilereader2_NPL("ttbarV_ttZToLLNuNu_ext3").at(0) + textfilereader2_NPL("ttbarV_ttZToQQ").at(0);

			   NPL_Ratio_Denominator_ttZ = textfilereader2_NPL("ttbarV_ttZToLLNuNu").at(1)      + textfilereader2_NPL("ttbarV_ttZToLL_ext2").at(1) +
                                                       textfilereader2_NPL("ttbarV_ttZToLLNuNu_ext3").at(1) + textfilereader2_NPL("ttbarV_ttZToQQ").at(1);

			   NPL_Ratio_Numerator_WZ = textfilereader2_NPL("VV_WZTo2l2Q").at(0) + textfilereader2_NPL("VV_WZTo1l1Nu2Q").at(0);
			   NPL_Ratio_Denominator_WZ = textfilereader2_NPL("VV_WZTo2l2Q").at(1) + textfilereader2_NPL("VV_WZTo1l1Nu2Q").at(1); 

			   break;

		case 2017: NPL_Ratio_Numerator_ttW = textfilereader2_NPL("ttbarV_ttWJetsToLNu").at(0) + textfilereader2_NPL("ttbarV_ttWJetsToQQ").at(0);
			   NPL_Ratio_Denominator_ttW = textfilereader2_NPL("ttbarV_ttWJetsToLNu").at(1) + textfilereader2_NPL("ttbarV_ttWJetsToQQ").at(1);

			   NPL_Ratio_Numerator_ttZ = textfilereader2_NPL("ttbarV_ttZToLL").at(0) + textfilereader2_NPL("ttbarV_ttZToLLNuNu").at(0) +
						     textfilereader2_NPL("ttbarV_ttZToQQ").at(0) + textfilereader2_NPL("ttbarV_ttZToQQ_ext").at(0);
	
			   NPL_Ratio_Denominator_ttZ = textfilereader2_NPL("ttbarV_ttZToLL").at(1) + textfilereader2_NPL("ttbarV_ttZToLLNuNu").at(1) +
                                                       textfilereader2_NPL("ttbarV_ttZToQQ").at(1) + textfilereader2_NPL("ttbarV_ttZToQQ_ext").at(1);
 

			   NPL_Ratio_Numerator_WZ = textfilereader2_NPL("VV_WZTo2l2Q").at(0) + textfilereader2_NPL("VV_WZTo1l1Nu2Q").at(0) + 
						    textfilereader2_NPL("VV_WZTo3lNu").at(0);

			   NPL_Ratio_Denominator_WZ = textfilereader2_NPL("VV_WZTo2l2Q").at(1) + textfilereader2_NPL("VV_WZTo1l1Nu2Q").at(1) + 
                                                      textfilereader2_NPL("VV_WZTo3lNu").at(1);

			   break;

		case 2018: NPL_Ratio_Numerator_ttW = textfilereader2_NPL("ttbarV_ttWJetsToLNu").at(0) + textfilereader2_NPL("ttbarV_ttWJetsToQQ").at(0);
			   NPL_Ratio_Denominator_ttW = textfilereader2_NPL("ttbarV_ttWJetsToLNu").at(1) + textfilereader2_NPL("ttbarV_ttWJetsToQQ").at(1);

			   NPL_Ratio_Numerator_ttZ = textfilereader2_NPL("ttbarV_ttZToLL").at(0) + textfilereader2_NPL("ttbarV_ttZToLLNuNu").at(0) +
                                                     textfilereader2_NPL("ttbarV_ttZToQQ").at(0) + textfilereader2_NPL("ttbarV_ttZToQQ_ext").at(0);

			   NPL_Ratio_Denominator_ttZ = textfilereader2_NPL("ttbarV_ttZToLL").at(1) + textfilereader2_NPL("ttbarV_ttZToLLNuNu").at(1) +
                                                       textfilereader2_NPL("ttbarV_ttZToQQ").at(1) + textfilereader2_NPL("ttbarV_ttZToQQ_ext").at(1);


			   NPL_Ratio_Numerator_WZ = textfilereader2_NPL("VV_WZTo2l2Q").at(0) + textfilereader2_NPL("VV_WZTo1l1Nu2Q").at(0) + 
                                                    textfilereader2_NPL("VV_WZTo3lNu").at(0) + textfilereader2_NPL("VV_WZTo3lNu_ext").at(0);

                           NPL_Ratio_Denominator_WZ = textfilereader2_NPL("VV_WZTo2l2Q").at(1) + textfilereader2_NPL("VV_WZTo1l1Nu2Q").at(1) +
                                                      textfilereader2_NPL("VV_WZTo3lNu").at(1) + textfilereader2_NPL("VV_WZTo3lNu_ext").at(1);

			   break;

		default: std::cout << "Error: The year must be 2016, 2017 or 2018" << std::endl; break;
        
        }

	NPL_Ratio_Numerator = textfilereader2_NPL("tZq").at(0) + textfilereader2_NPL("SingleTop_tHq").at(0) + 
			      NPL_Ratio_Numerator_ttW 	       + NPL_Ratio_Numerator_ttZ		    +
			      NPL_Ratio_Numerator_WZ;

	NPL_Ratio_Denominator = textfilereader2_NPL("tZq").at(1) + textfilereader2_NPL("SingleTop_tHq").at(1) +               
			        NPL_Ratio_Denominator_ttW 	 + NPL_Ratio_Denominator_ttZ		      + 
				NPL_Ratio_Denominator_WZ;
  
        float NPL_Weight = NPL_Ratio_Numerator / NPL_Ratio_Denominator;

	std::ofstream NPLWeightFile;
        std::string NPLWeightString = "NPLWeight_" + Year_String + ".txt";

        NPLWeightFile.open(NPLWeightString.c_str());
        NPLWeightFile << NPL_Weight << std::endl;

  }



}




