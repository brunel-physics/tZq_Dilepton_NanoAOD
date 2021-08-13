#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TH1D.h"
#include <boost/numeric/conversion/cast.hpp>
#include <algorithm>
#include <TFile.h>
#include <iostream>
#include <fstream>

using namespace ROOT; // RDataFrame's namespace
using namespace std;


void process_comparison_plotter(const string& year, const string& variable_name, const string& channel, const string& xaxis_name, const string& systematic, const string& region, const string& GraphTitle){

  std::string tZq_file = "tZq_" + systematic + "_" + channel + "_" + region + "_" + year + "_Combined.root";
  std::string ZPlusJets_file = "ZPlusJets_" + systematic + "_" + channel + "_" + region + "_" + year + "_Combined.root";
  std::string ttbar_file = "ttbar_" + systematic + "_" + channel + "_" + region + "_" + year + "_Combined.root";
  std::string SingleTop_file = "SingleTop_" + systematic + "_" + channel + "_" + region + "_" + year + "_Combined.root";
  std::string VV_file = "VV_" + systematic + "_" + channel + "_" + region + "_" + year + "_Combined.root";
  std::string VVV_file = "VVV_" + systematic + "_" + channel + "_" + region + "_" + year + "_Combined.root";
  std::string WPlusJets_file = "WPlusJets_" + systematic + "_" + channel + "_" + region + "_" + year + "_Combined.root";
  std::string ttbarV_file = "ttbarV_" + systematic + "_" + channel + "_" + region + "_" + year + "_Combined.root";
  std::string data_file = "data_" + channel + "_" + region + "_" + year + "_Combined.root";

  RDataFrame d_tZq("Events", tZq_file.c_str());
  RDataFrame d_ZPlusJets("Events", ZPlusJets_file.c_str());
  RDataFrame d_ttbar("Events", ttbar_file.c_str());
  RDataFrame d_SingleTop("Events", SingleTop_file.c_str());
  //RDataFrame d_VV("Events", VV_file.c_str());
  RDataFrame d_VVV("Events", VVV_file.c_str());
  RDataFrame d_WPlusJets("Events", WPlusJets_file.c_str());
  RDataFrame d_ttbarV("Events", ttbarV_file.c_str());
  RDataFrame d_data("Events", data_file.c_str());

  double min; double max; 

  int nbins = 25;

  if(xaxis_name == "Mass [GeV]"){min = 0; max = 1000;}
  else if(xaxis_name == "p_{T} [GeV]"){min = 0; max = 500;}
  else if(xaxis_name == "#eta"){min = -4; max = 4;}
  else if(xaxis_name == "#phi"){min = -5; max = 5;}
  else{throw std::logic_error("ERROR: Check axis_name");}

  auto Zeroes{[](){

  	ROOT::VecOps::RVec<double> Out(1000, 0.);
	return Out;
  }};

  auto IncreaseEventWeight_1000{[](const double& EventWeightInput){

	return EventWeightInput * 1000;

  }};



  //Applying the event weight to the histograms
  std::cout << "print 1" << std::endl;
  auto h_tZq = d_tZq.Histo1D({"h_tZq", variable_name.c_str(), nbins, min, max}, variable_name.c_str(), "EventWeight");
  std::cout << "print 2" << std::endl;
  auto h_ZPlusJets = d_ZPlusJets.Histo1D({"h_ZPlusJets", variable_name.c_str(), nbins, min, max}, variable_name.c_str(), "EventWeight");
  std::cout << "print 3" << std::endl;
  auto h_ttbar = d_ttbar.Histo1D({"h_ttbar", variable_name.c_str(), nbins, min, max}, variable_name.c_str(), "EventWeight");
  std::cout << "print 4" << std::endl;
  auto h_SingleTop = d_SingleTop.Histo1D({"h_SingleTop", variable_name.c_str(), nbins, min, max}, variable_name.c_str(), "EventWeight");
  std::cout << "print 5" << std::endl;
  //auto h_VV = d_VV.Histo1D({"h_VV", variable_name.c_str(), nbins, min, max}, variable_name.c_str(), "EventWeight");
  std::cout << "print 6" << std::endl;
  auto h_VVV = d_VVV.Histo1D({"h_VVV", variable_name.c_str(), nbins, min, max}, variable_name.c_str(), "EventWeight");
  std::cout << "print 7" << std::endl;
  auto h_WPlusJets = d_WPlusJets.Histo1D({"h_WPlusJets", variable_name.c_str(), nbins, min, max}, variable_name.c_str(), "EventWeight");
  std::cout << "print 8" << std::endl;
  auto h_ttbarV = d_ttbarV.Histo1D({"h_ttbarV", variable_name.c_str(), nbins, min, max}, variable_name.c_str(), "EventWeight");
  std::cout << "print 9" << std::endl;
  auto h_data = d_ttbarV.Define("VecOfZeroes", Zeroes, {}).Histo1D({"h_data", variable_name.c_str(), nbins, min, max}, "VecOfZeroes");
  //auto h_data = d_data.Histo1D({"h_data", variable_name.c_str(), nbins, min, max}, variable_name.c_str(), "EventWeight");
  std::cout << "print 10" << std::endl;

  //For the canvas
  int W = 800;
  int H = 600;
  int H_ref = 600;
  int W_ref = 800;

  float T = 0.08 * H_ref;
  float B = 0.12 * H_ref;
  float L = 0.12 * W_ref;
  float R = 0.04 * W_ref;

  TCanvas* c1 =new TCanvas("c1"," ",50,50,W,H);
  c1->cd();
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetFrameFillStyle(0);
  c1->SetFrameBorderMode(0);
  c1->SetLeftMargin(L / W);
  c1->SetRightMargin(R / W);
  c1->SetTopMargin(T / H);
  c1->SetBottomMargin(B / H);
  c1->SetTickx(0);
  c1->SetTicky(0);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  h_tZq->SetFillColor(600); //kBlue
  h_ZPlusJets->SetFillColor(632); //kRed
  h_ttbar->SetFillColor(418); //kGreen+2
  h_SingleTop->SetFillColor(618); //kMagenta+2
  //h_VV->SetFillColor(802); //kOrange+2
  h_VVV->SetFillColor(400); //kYellow
  h_WPlusJets->SetFillColor(425); //kCyan-7
  h_ttbarV->SetFillColor(920); //kGray
  h_tZq->SetLineColor(600); //kBlue
  h_ZPlusJets->SetLineColor(632); //kRed
  h_ttbar->SetLineColor(418); //kGreen+2
  h_SingleTop->SetLineColor(618); //kMagenta+2
  //h_VV->SetLineColor(802); //kOrange+2
  h_VVV->SetLineColor(400); //kYellow
  h_WPlusJets->SetLineColor(425); //kCyan-7
  h_ttbarV->SetLineColor(920); //kGray
  h_data->SetMarkerColor(1); //kBlack
  h_data->SetMarkerStyle(20);
  h_data->SetMarkerSize(1.0);

  h_tZq->SetTitle("tZq");
  h_ZPlusJets->SetTitle("Z+jets");
  h_ttbar->SetTitle("t#bar{t}");
  h_SingleTop->SetTitle("Single top");
  //h_VV->SetTitle("VV");
  h_VVV->SetTitle("VVV");
  h_WPlusJets->SetTitle("W+jets");
  h_ttbarV->SetTitle("t#bar{t}V");
  h_data->SetTitle("data");

  h_tZq->Scale(1./h_tZq->Integral());
  h_ZPlusJets->Scale(1./h_ZPlusJets->Integral());
  h_ttbar->Scale(1./h_ttbar->Integral());
  h_SingleTop->Scale(1./h_SingleTop->Integral());
  //h_VV->Scale(1./h_VV->Integral());
  h_VVV->Scale(1./h_VVV->Integral());
  h_WPlusJets->Scale(1./h_WPlusJets->Integral());
  h_ttbarV->Scale(1./h_ttbarV->Integral());
  h_data->Scale(1./h_data->Integral());

  TPad *pad = new TPad("pad","pad",0.01, 0.315, 0.99, 0.99);
  pad->SetTopMargin(0);
  pad->SetFillColor(0);
  pad->SetBorderMode(0);
  pad->SetFrameFillStyle(0);
  pad->SetFrameBorderMode(0);
  pad->SetLeftMargin(L / W);
  pad->SetRightMargin(R / W);
  pad->SetTopMargin(T / H);
  pad->SetBottomMargin(B / H * 0.3);
  pad->SetLogy();
  pad->SetTickx(0);
  pad->SetTicky(0);
  pad->Draw();
  pad->cd();

  THStack *MC_Stack = new THStack("MC_Stack",GraphTitle.c_str());
  MC_Stack->Add(h_VVV.GetPtr());
  //MC_Stack->Add(h_VV.GetPtr());
  MC_Stack->Add(h_ttbarV.GetPtr());
  MC_Stack->Add(h_ttbar.GetPtr());
  MC_Stack->Add(h_ZPlusJets.GetPtr());
  MC_Stack->Add(h_WPlusJets.GetPtr());
  MC_Stack->Add(h_SingleTop.GetPtr());

  //normalising the stack: https://root-forum.cern.ch/t/normalize-thsatck/21807/6  
  TH1F *MC_Stack_Normalised = new TH1F(*((TH1F *)(MC_Stack->GetStack()->Last())));
  //MC_Stack_Normalised->SetNameTitle("MC_Stack_Normalised", "my normalized stack");
  MC_Stack_Normalised->Scale(1./MC_Stack_Normalised->Integral("width"));

  //MC_Stack->Scale(1./MC_Stack->Integral());

  //MC_Stack->Draw("HIST");
  MC_Stack_Normalised->Draw("HIST");
  h_tZq->Draw("HIST SAME");
  h_data->Draw("HIST PSAME");

  gPad->BuildLegend();

/*
  MC_Stack->GetHistogram()->GetXaxis()->SetTickLength(0);
  MC_Stack->GetHistogram()->GetXaxis()->SetLabelOffset(999);
  MC_Stack->GetHistogram()->GetYaxis()->SetTitle("Events");
  MC_Stack->GetHistogram()->GetYaxis()->SetTitleOffset(0.7);
  MC_Stack->GetHistogram()->GetYaxis()->SetTitleSize(0.05);
*/

  MC_Stack_Normalised->GetXaxis()->SetTickLength(0);
  MC_Stack_Normalised->GetXaxis()->SetLabelOffset(999);
  MC_Stack_Normalised->GetYaxis()->SetTitle("Events");
  MC_Stack_Normalised->GetYaxis()->SetTitleOffset(0.7);
  MC_Stack_Normalised->GetYaxis()->SetTitleSize(0.05);

  TPaveText* ptext1 = new TPaveText(0.1, 1.0, 0.6, 0.94, "NDCCBR");
  TText *t1=ptext1->AddText(GraphTitle.c_str());
  ptext1->SetFillStyle(0);
  ptext1->SetBorderSize(0);
  ptext1->Draw();

  TPaveText* ptext3 = new TPaveText(0.7, 1.0, 0.95, 0.7, "NDCCBR");
  TText *t3=ptext3->AddText("CMS Work in progress");
  t3->SetTextFont(52);
  ptext3->SetFillStyle(0);
  ptext3->SetBorderSize(0);
  ptext3->Draw();

  std::string lumi;

  if(year == "2016"){lumi = "35.88 fb^{-1} (13 TeV)";}
  else if(year == "2017"){lumi = "41.86 fb^{-1} (13 TeV)";}
  else if(year == "2018"){lumi = "59.699 fb^{-1} (13 TeV)";}
  else{throw std::logic_error("ERROR: Year must be 2016, 2017 or 2018.");}

  TPaveText* ptext4 = new TPaveText(0.7, 1.0, 1.0, 0.94, "NDCCBR");
  TText *t4=ptext4->AddText(lumi.c_str());
  ptext4->SetFillStyle(0);
  ptext4->SetBorderSize(0);
  ptext4->Draw();

  c1->cd();

  TPad *pad2 = new TPad("pad2", "pad2", 0.01, 0.01, 0.99, 0.3275);
  pad2->SetTopMargin(0);
  pad2->SetFillColor(0);
  pad2->SetBorderMode(0);
  pad2->SetFrameFillStyle(0);
  pad2->SetFrameBorderMode(0);
  pad2->SetLeftMargin(L / W);
  pad2->SetRightMargin(R / W);
  pad2->SetTopMargin(T / H);
  pad2->SetBottomMargin(B / H * 2.1);
  pad2->SetTickx(0);
  pad2->SetTicky(0);
  pad2->SetGridy(1);
  pad2->Draw();
  pad2->cd();

  //for the ratio plot
  TH1D * rp = (TH1D*)(h_data->Clone());
  rp->Divide((TH1D*)MC_Stack->GetStack()->Last());
  rp->SetStats(false);
  rp->GetXaxis()->SetNdivisions(6, 5, 0);
  rp->GetYaxis()->SetNdivisions(6, 5, 0);
  rp->GetXaxis()->SetLabelSize(0.08);
  rp->GetXaxis()->SetLabelOffset(0.030);
  rp->GetXaxis()->SetTitleSize(0.11);
  rp->GetYaxis()->SetLabelSize(0.08);
  rp->GetYaxis()->SetTitle("Data/MC");
  rp->GetYaxis()->SetTitleSize(0.12);
  rp->GetYaxis()->SetTitleOffset(L / W * 3.);
  rp->GetYaxis()->CenterTitle();
  rp->GetXaxis()->SetTitle(xaxis_name.c_str());
  rp->SetMinimum(0.5);
  rp->SetMaximum(1.5);
  rp->Draw();

  //Updating the canvas and saving the output pdf file
  std::string OutputPlot = variable_name + "_" + channel + "_" + systematic + "_" + region + "_" + year + ".pdf";
  
  c1->Update();
  c1->SaveAs(OutputPlot.c_str());

}



void process_comparison_2016(){

  gSystem->Exec("mkdir Plots");

  process_comparison_plotter("2016", "InvTopMass", "ee", "Mass [GeV]", "Nominal", "NoChi2Cut", "Invariant mass of the top quark candidate (ee channel)");
/*  process_comparison_plotter("2016", "w_mass", "ee", "Mass [GeV]", "Nominal", "NoChi2Cut", "Reconstructed mass of the W quark candidate (ee channel)");
  process_comparison_plotter("2016", "z_mass", "ee", "Mass [GeV]", "Nominal", "NoChi2Cut", "Reconstructed mass of the Z boson candidate (ee channel)");
  process_comparison_plotter("2016", "LeadingLeptonPt", "ee", "p_{T} [GeV]", "Nominal", "NoChi2Cut", "Transverse momentum of the leading lepton (ee channel)");
  process_comparison_plotter("2016", "SubleadingLeptonPt", "ee", "p_{T} [GeV]", "Nominal", "NoChi2Cut", "Transverse momentum of the subleading lepton (ee channel)");
  process_comparison_plotter("2016", "LeadingLeptonEta", "ee", "#eta", "Nominal", "NoChi2Cut", "#eta of the leading lepton (ee channel)");
  process_comparison_plotter("2016", "SubleadingLeptonEta", "ee", "#eta", "Nominal", "NoChi2Cut", "#eta of the subleading lepton (ee channel)");
  process_comparison_plotter("2016", "LeadingLeptonPhi", "ee", "#phi", "Nominal", "NoChi2Cut", "#phi of the leading lepton (ee channel)");
  process_comparison_plotter("2016", "SubleadingLeptonPhi", "ee", "#phi", "Nominal", "NoChi2Cut", "#phi of the subleading lepton (ee channel)");
  process_comparison_plotter("2016", "LeadingLeptonMass", "ee", "Mass [GeV]", "Nominal", "NoChi2Cut", "Mass of the leading lepton (ee channel)");
  process_comparison_plotter("2016", "SubleadingLeptonMass", "ee", "Mass [GeV]", "Nominal", "NoChi2Cut", "Mass of the subleading lepton (ee channel)");
  process_comparison_plotter("2016", "SmearedJetPt", "ee", "p_{T} [GeV]", "Nominal", "NoChi2Cut", "Transverse momenta of smeared jets (ee channel)");
  process_comparison_plotter("2016", "SmearedJetPhi", "ee", "#phi", "Nominal", "NoChi2Cut", "#phi of smeared jets (ee channel)");
  process_comparison_plotter("2016", "SmearedJetEta", "ee", "#eta", "Nominal", "NoChi2Cut", "#eta of smeared jets (ee channel)");
  process_comparison_plotter("2016", "SmearedJetMass", "ee", "Mass [GeV]", "Nominal", "NoChi2Cut", "Mass of smeared jets (ee channel)");
  process_comparison_plotter("2016", "TightSmearedJetsPt", "ee", "p_{T} [GeV]", "Nominal", "NoChi2Cut", "Transverse momenta of tight smeared jets (ee channel)");
  process_comparison_plotter("2016", "TightSmearedJetsPhi", "ee", "#phi", "Nominal", "NoChi2Cut", "#phi of tight smeared jets (ee channel)");
  process_comparison_plotter("2016", "TightSmearedJetsEta", "ee", "#eta", "Nominal", "NoChi2Cut", "#eta of tight smeared jets (ee channel)");
  process_comparison_plotter("2016", "TightSmearedJetsMass", "ee", "Mass [GeV]", "Nominal", "NoChi2Cut", "Mass of tight smeared jets (ee channel)");
  process_comparison_plotter("2016", "LeadingJetPt", "ee", "p_{T} [GeV]", "Nominal", "NoChi2Cut", "Transverse momenta of the leading tight smeared jet (ee channel)");
  process_comparison_plotter("2016", "LeadingJetPhi", "ee", "#phi", "Nominal", "NoChi2Cut", "#phi of the leading tight smeared jet (ee channel)");
  process_comparison_plotter("2016", "LeadingJetEta", "ee", "#eta", "Nominal", "NoChi2Cut", "#eta of the leading tight smeared jet (ee channel)");
  process_comparison_plotter("2016", "LeadingJetMass", "ee", "Mass [GeV]", "Nominal", "NoChi2Cut", "Mass of the leading tight smeared jet (ee channel)");
  process_comparison_plotter("2016", "SubleadingJetPt", "ee", "p_{T} [GeV]", "Nominal", "NoChi2Cut", "Transverse momenta of the subleading tight smeared jet (ee channel)");
  process_comparison_plotter("2016", "SubleadingJetPhi", "ee", "#phi", "Nominal", "NoChi2Cut", "#phi of the subleading tight smeared jet (ee channel)");
  process_comparison_plotter("2016", "SubleadingJetEta", "ee", "#eta", "Nominal", "NoChi2Cut", "#eta of the subleading tight smeared jet (ee channel)");
  process_comparison_plotter("2016", "SubleadingJetMass", "ee", "Mass [GeV]", "Nominal", "NoChi2Cut", "Mass of the subleading tight smeared jet (ee channel)");
  process_comparison_plotter("2016", "ThirdJetPt", "ee", "p_{T} [GeV]", "Nominal", "NoChi2Cut", "Transverse momenta of the third tight smeared jet (ee channel)");
  process_comparison_plotter("2016", "ThirdJetPhi", "ee", "#phi", "Nominal", "NoChi2Cut", "#phi of the third tight smeared jet (ee channel)");
  process_comparison_plotter("2016", "ThirdJetEta", "ee", "#eta", "Nominal", "NoChi2Cut", "#eta of the third tight smeared jet (ee channel)");
  process_comparison_plotter("2016", "ThirdJetMass", "ee", "Mass [GeV]", "Nominal", "NoChi2Cut", "Mass of the third tight smeared jet (ee channel)"); 
  process_comparison_plotter("2016", "FourthJetPt", "ee", "p_{T} [GeV]", "Nominal", "NoChi2Cut", "Transverse momenta of the fourth tight smeared jet (ee channel)");
  process_comparison_plotter("2016", "FourthJetPhi", "ee", "#phi", "Nominal", "NoChi2Cut", "#phi of the fourth tight smeared jet (ee channel)");
  process_comparison_plotter("2016", "FourthJetEta", "ee", "#eta", "Nominal", "NoChi2Cut", "#eta of the fourth tight smeared jet (ee channel)");
  process_comparison_plotter("2016", "FourthJetMass", "ee", "Mass [GeV]", "Nominal", "NoChi2Cut", "Mass of the fourth tight smeared jet (ee channel)");
*/
  gSystem->Exec("mv *.pdf Plots/");

}
