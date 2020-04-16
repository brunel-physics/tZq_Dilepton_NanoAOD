#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TH1D.h"
#include <TFile.h>
#include <iostream>
#include <fstream>

using namespace ROOT; // RDataFrame's namespace
using namespace std;

void SampleCheck(){


EnableImplicitMT();

TFile* output = new TFile("ElectronPt.root", "RECREATE");

vector<string> input_files = {"/data/disk0/nanoAOD_2017/DYJetsToLL_M10to50/*.root",
			      "/data/disk0/nanoAOD_2017/ttbar_2l2nu/*.root",
		              "/data/disk0/nanoAOD_2017/TTbar_inc_nominal/*.root",
		              "/data/disk0/nanoAOD_2017/ttbar_madgraph/*.root",
		              "/data/disk0/nanoAOD_2017/ttbar_madgraph_ext/*.root",
		              "/data/disk0/nanoAOD_2017/TTToHadronic/*.root",
		              "/data/disk0/nanoAOD_2017/TTToSemileptonic/*.root",
	                      "/data/disk0/nanoAOD_2017/ttbar_aMCatNLO/*.root",
		              "/data/disk0/nanoAOD_2017/ttbar_inc/*.root",
			      "/data/disk0/nanoAOD_2017/ST_tchannel_top/*.root",
		              "/data/disk0/nanoAOD_2017/ST_tchannel_tbar/*.root",
		              "/data/disk0/nanoAOD_2017/ST_schannel/*.root",
		              "/data/disk0/nanoAOD_2017/ST_tW/*.root",
		              "/data/disk0/nanoAOD_2017/ST_tbarW/*.root",
		              "/data/disk0/nanoAOD_2017/tHq/*.root", 
		              "/data/disk0/nanoAOD_2017/tZq_W_lept_Z_had/*root",
			      "/data/disk0/nanoAOD_2017/ZZTo2Q2Nu/*.root",
		       	      "/data/disk0/nanoAOD_2017/ZZTo2L2Nu/*.root",
		              "/data/disk0/nanoAOD_2017/ZZTo2L2Q/*.root",
		              "/data/disk0/nanoAOD_2017/ZZTo4L/*.root",
		              "/data/disk0/nanoAOD_2017/WZTo1L1Nu2Q/*.root",
                              "/data/disk0/nanoAOD_2017/WZTo2L2Q/*.root",
                              "/data/disk0/nanoAOD_2017/WZTo3LNu/*.root",
                              "/data/disk0/nanoAOD_2017/WWTo1L1Nu2Q/*.root",
		              "/data/disk0/nanoAOD_2017/WWTo2L2Nu/*.root",
                              "/data/disk0/nanoAOD_2017/WWToLNuQQ/*.root"
			      "/data/disk0/nanoAOD_2017/WWWTo4F/*.root",
		              "/data/disk0/nanoAOD_2017/WWZTo4F/*.root",
		              "/data/disk0/nanoAOD_2017/WZZ/*.root",
		              "/data/disk0/nanoAOD_2017/ZZZ/*.root",
			      "/data/disk0/nanoAOD_2017/WJetsToLNu/*.root",
			      "/data/disk0/nanoAOD_2017/ttWJetsToLNu/*.root",
		              "/data/disk0/nanoAOD_2017/ttWJetsToQQ/*.root",
		              "/data/disk0/nanoAOD_2017/ttgamma/*.root",
		              "/data/disk0/nanoAOD_2017/ttZToLL/*.root",
	                      "/data/disk0/nanoAOD_2017/ttHTobb/*.root",
	                      "/data/disk0/nanoAOD_2017/ttHToNonbb/*.root",
	                      "/data/disk0/nanoAOD_2017/ttZToLLNuNu/*.root",
	                      "/data/disk0/nanoAOD_2017/ttZToQQ/*.root",
	                      "/data/disk0/nanoAOD_2017/ttZToQQ_ext/*.root",
		              "/data/disk0/nanoAOD_2017/tWZ_tWll/*.root",
			      "/data/disk0/nanoAOD_2017/DoubleEGRun2017B/*.root",
		       	      "/data/disk0/nanoAOD_2017/DoubleEGRun2017C/*.root",
		              "/data/disk0/nanoAOD_2017/DoubleEGRun2017D/*.root",
		              "/data/disk0/nanoAOD_2017/DoubleEGRun2017E/*.root",
		              "/data/disk0/nanoAOD_2017/DoubleEGRun2017F/*.root",
		              "/data/disk0/nanoAOD_2017/DoubleMuonRun2017B/*.root",
		              "/data/disk0/nanoAOD_2017/DoubleMuonRun2017C/*.root",
		              "/data/disk0/nanoAOD_2017/DoubleMuonRun2017D/*.root",
		              "/data/disk0/nanoAOD_2017/DoubleMuonRun2017E/*.root",
		              "/data/disk0/nanoAOD_2017/DoubleMuonRun2017F/*.root"};

RDataFrame d_dataframe("Events", input_files);

auto ElectronPt = d_dataframe.Histo1D("Electron_pt");
ElectronPt->Draw();
ElectronPt->Write();
output->Close();



}
