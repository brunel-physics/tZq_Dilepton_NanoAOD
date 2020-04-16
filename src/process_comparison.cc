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


void process_comparison_plotter(const TString& variable_name, const TString& xaxis_name, const TString& GraphTitle){


//The root files for normalised histograms were combined for each process using the linux command "hadd".
//This was achieved using the bash script hadd_NormalisedHistos.sh
//To run that bash script, type "./hadd_NormalisedHistos.sh" into your terminal


//Reading the root files that contain the combined histograms
TFile *tZq_ll = TFile::Open("tZq_Combined_AfterFullSelection.root", "READ");
TFile *ZPlusJets = TFile::Open("ZPlusJets_Combined_AfterFullSelection.root", "READ");
TFile *ttbar = TFile::Open("ttbar_Combined_AfterFullSelection.root", "READ");
TFile *SingleTop = TFile::Open("SingleTop_Combined_AfterFullSelection.root", "READ");
TFile *VV = TFile::Open("VV_Combined_AfterFullSelection.root", "READ");
TFile *VVV = TFile::Open("VVV_Combined_AfterFullSelection.root", "READ");
TFile *WPlusJets = TFile::Open("WPlusJets_Combined_AfterFullSelection.root", "READ");
TFile *ttbarV = TFile::Open("ttbarV_Combined_AfterFullSelection.root", "READ");
TFile *data = TFile::Open("data_Combined_AfterFullSelection.root", "READ");


//Histograms for events that pass all cuts
TH1* h_tZq_Combined_AllCuts = (TH1*)tZq_ll->GetObjectChecked("h_" + variable_name + "_AllCuts_Blinded", "TH1");
TH1* h_ZPlusJets_Combined_AllCuts = (TH1*)ZPlusJets->GetObjectChecked("h_" + variable_name + "_AllCuts_Blinded", "TH1");
TH1* h_ttbar_Combined_AllCuts = (TH1*)ttbar->GetObjectChecked("h_" + variable_name + "_AllCuts_Blinded", "TH1");	
TH1* h_SingleTop_Combined_AllCuts = (TH1*)SingleTop->GetObjectChecked("h_" + variable_name + "_AllCuts_Blinded", "TH1");
TH1* h_VV_Combined_AllCuts = (TH1*)VV->GetObjectChecked("h_" + variable_name + "_AllCuts_Blinded", "TH1");
TH1* h_VVV_Combined_AllCuts = (TH1*)VVV->GetObjectChecked("h_" + variable_name + "_AllCuts_Blinded", "TH1");
TH1* h_WPlusJets_Combined_AllCuts = (TH1*)WPlusJets->GetObjectChecked("h_" + variable_name + "_AllCuts_Blinded", "TH1");
TH1* h_ttbarV_Combined_AllCuts = (TH1*)ttbarV->GetObjectChecked("h_" + variable_name + "_AllCuts_Blinded", "TH1");
TH1* h_data_Combined_AllCuts = (TH1*)data->GetObjectChecked("h_" + variable_name + "_AllCuts_Blinded", "TH1");


h_tZq_Combined_AllCuts->Sumw2();
h_ZPlusJets_Combined_AllCuts->Sumw2();
h_ttbar_Combined_AllCuts->Sumw2();
h_SingleTop_Combined_AllCuts->Sumw2();
h_VV_Combined_AllCuts->Sumw2();
h_VVV_Combined_AllCuts->Sumw2();
h_WPlusJets_Combined_AllCuts->Sumw2();
h_ttbarV_Combined_AllCuts->Sumw2();
h_data_Combined_AllCuts->Sumw2();

int W = 800;
int H = 600;
int H_ref = 600;
int W_ref = 800;

// references for T, B, L, R
float T = 0.08 * H_ref;
float B = 0.12 * H_ref;
float L = 0.12 * W_ref;
float R = 0.04 * W_ref;

//TCanvas* c1 =new TCanvas("c1"," ",200,10,700,500);
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

h_tZq_Combined_AllCuts->SetFillColor(600); //kBlue
h_ZPlusJets_Combined_AllCuts->SetFillColor(632); //kRed
h_ttbar_Combined_AllCuts->SetFillColor(418); //kGreen+2
h_SingleTop_Combined_AllCuts->SetFillColor(618); //kMagenta+2
h_VV_Combined_AllCuts->SetFillColor(802); //kOrange+2
h_VVV_Combined_AllCuts->SetFillColor(400); //kYellow
h_WPlusJets_Combined_AllCuts->SetFillColor(425); //kCyan-7
h_ttbarV_Combined_AllCuts->SetFillColor(920); //kGray
h_tZq_Combined_AllCuts->SetLineColor(600); //kBlue
h_ZPlusJets_Combined_AllCuts->SetLineColor(632); //kRed
h_ttbar_Combined_AllCuts->SetLineColor(418); //kGreen+2
h_SingleTop_Combined_AllCuts->SetLineColor(618); //kMagenta+2
h_VV_Combined_AllCuts->SetLineColor(802); //kOrange+2
h_VVV_Combined_AllCuts->SetLineColor(400); //kYellow
h_WPlusJets_Combined_AllCuts->SetLineColor(425); //kCyan-7
h_ttbarV_Combined_AllCuts->SetLineColor(920); //kGray
h_data_Combined_AllCuts->SetMarkerColor(1); //kBlack
h_data_Combined_AllCuts->SetMarkerStyle(20);
h_data_Combined_AllCuts->SetMarkerSize(1.0);

h_tZq_Combined_AllCuts->SetTitle("tZq");
h_ZPlusJets_Combined_AllCuts->SetTitle("Z+jets");
h_ttbar_Combined_AllCuts->SetTitle("t#bar{t}");
h_SingleTop_Combined_AllCuts->SetTitle("Single top");
h_VV_Combined_AllCuts->SetTitle("VV");
h_VVV_Combined_AllCuts->SetTitle("VVV");
h_WPlusJets_Combined_AllCuts->SetTitle("W+jets");
h_ttbarV_Combined_AllCuts->SetTitle("t#bar{t}V");
h_data_Combined_AllCuts->SetTitle("data");


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
pad->SetTickx(0);
pad->SetTicky(0);
pad->Draw();
pad->cd();



THStack *MC_Stack = new THStack("MC_Stack",GraphTitle);
MC_Stack->Add(h_ZPlusJets_Combined_AllCuts);
MC_Stack->Add(h_ttbar_Combined_AllCuts);
MC_Stack->Add(h_VV_Combined_AllCuts);
MC_Stack->Add(h_WPlusJets_Combined_AllCuts);
MC_Stack->Add(h_SingleTop_Combined_AllCuts);
MC_Stack->Add(h_ttbarV_Combined_AllCuts);
MC_Stack->Add(h_VVV_Combined_AllCuts);
MC_Stack->Add(h_tZq_Combined_AllCuts);
MC_Stack->Draw("HIST");
h_data_Combined_AllCuts->Draw("HIST PSAME");

gPad->BuildLegend();

MC_Stack->GetHistogram()->GetXaxis()->SetTickLength(0);
MC_Stack->GetHistogram()->GetXaxis()->SetLabelOffset(999);
MC_Stack->GetHistogram()->GetYaxis()->SetTitle("Events");
MC_Stack->GetHistogram()->GetYaxis()->SetTitleOffset(0.7);
MC_Stack->GetHistogram()->GetYaxis()->SetTitleSize(0.05);

TPaveText* ptext1 = new TPaveText(0.1, 1.0, 0.6, 0.94, "NDCCBR");
TText *t1=ptext1->AddText(GraphTitle);
ptext1->SetFillStyle(0);
ptext1->SetBorderSize(0);
ptext1->Draw();

//TPaveText* ptext3 = new TPaveText(0.7, 1.0, 0.9, 0.7, "NDCCBR");
TPaveText* ptext3 = new TPaveText(0.7, 1.0, 0.95, 0.7, "NDCCBR");
TText *t3=ptext3->AddText("CMS Work in progress");
t3->SetTextFont(52);
ptext3->SetFillStyle(0);
ptext3->SetBorderSize(0);
ptext3->Draw();

TPaveText* ptext4 = new TPaveText(0.7, 1.0, 1.0, 0.94, "NDCCBR");
TText *t4=ptext4->AddText("41.86 fb^{-1} (13 TeV)");
ptext4->SetFillStyle(0);
ptext4->SetBorderSize(0);
ptext4->Draw();


c1->cd();

//TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
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


TH1D * rp = (TH1D*)(h_data_Combined_AllCuts->Clone());
rp->Divide((TH1D*)MC_Stack->GetStack()->Last());

cout << "MC entries in the stack: " << ((TH1D*)MC_Stack->GetStack()->Last())->GetEntries() << endl;
cout << "Data entries: " << h_data_Combined_AllCuts->GetEntries() << endl;

for(int i = 0; i < rp->GetEntries(); i++){
        cout << "Ratio for " << variable_name << " is " << rp->GetBinContent(i) << endl;
}

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
rp->GetXaxis()->SetTitle(xaxis_name);
rp->SetMinimum(0.5);
rp->SetMaximum(1.5);
rp->Draw();
c1->Update();
c1->SaveAs(variable_name + "_AllCuts.pdf");


}



void process_comparison(){

//ee
process_comparison_plotter("Top_Mass_ee", "Mass [GeV]", "Mass of the top quark candidate (ee channel)");
process_comparison_plotter("Top_Pt_ee", "p_{T} [GeV]", "Transverse momentum of the top quark candidate (ee channel)");
//process_comparison_plotter("Top_Eta_ee", "#eta", "#eta of the top quark candidate (ee channel)");
//process_comparison_plotter("Top_Phi_ee", "#phi", "#phi of the top quark candidate (ee channel)");
process_comparison_plotter("w_mass_ee", "Mass [GeV]", "Mass of the W boson candidate (ee channel)");
//process_comparison_plotter("w_pair_pt_ee", "p_{T} [GeV]", "Transverse momentum of the W boson candidate (ee channel)");
//process_comparison_plotter("w_pair_eta_ee", "#eta", "#eta of the W boson candidate (ee channel)");
//process_comparison_plotter("w_pair_phi_ee", "#phi", "#phi of the W boson candidate (ee channel)");
process_comparison_plotter("z_mass_ee", "Mass [GeV]", "Mass of the Z boson candidate (ee channel)");
//process_comparison_plotter("RecoZPt_ee", "p_{T} [GeV]", "Transverse momentum of the Z boson candidate (ee channel)");
//process_comparison_plotter("RecoZEta_ee", "#eta", "#eta of the Z boson candidate (ee channel)");
//process_comparison_plotter("RecoZPhi_ee", "#phi", "#phi of the Z boson candidate (ee channel)");

//mumu
process_comparison_plotter("Top_Mass_mumu", "Mass [GeV]", "Mass of the top quark candidate (mumu channel)");
process_comparison_plotter("Top_Pt_mumu", "p_{T} [GeV]", "Transverse momentum of the top quark candidate (mumu channel)");
//process_comparison_plotter("Top_Eta_mumu", "#eta", "#eta of the top quark candidate (mumu channel)");
//process_comparison_plotter("Top_Phi_mumu", "#phi", "#phi of the top quark candidate (mumu channel)");
process_comparison_plotter("w_mass_mumu", "Mass [GeV]", "Mass of the W boson candidate (mumu channel)");
//process_comparison_plotter("w_pair_pt_mumu", "p_{T} [GeV]", "Transverse momentum of the W boson candidate (mumu channel)");
//process_comparison_plotter("w_pair_eta_mumu", "#eta", "#eta of the W boson candidate (mumu channel)");
//process_comparison_plotter("w_pair_phi_mumu", "#phi", "#phi of the W boson candidate (mumu channel)");
process_comparison_plotter("z_mass_mumu", "Mass [GeV]", "Mass of the Z boson candidate (mumu channel)");
//process_comparison_plotter("RecoZPt_mumu", "p_{T} [GeV]", "Transverse momentum of the Z boson candidate (mumu channel)");
//process_comparison_plotter("RecoZEta_mumu", "#eta", "#eta of the Z boson candidate (mumu channel)");
//process_comparison_plotter("RecoZPhi_mumu", "#phi", "#phi of the Z boson candidate (mumu channel)");


}
