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


void process_comparison_plotter(const TString& variable_name, const TString& channel, const TString& xaxis_name, const TString& systematic, const TString& region, const TString& GraphTitle){

  RDataFrame d_tZq("Events", {"tZq_Combined.root"});
  RDataFrame d_ZPlusJets("Events", {"ZPlusJets_Combined.root"});
  RDataFrame d_ttbar("Events", {"ttbar_Combined.root"});
  RDataFrame d_SingleTop("Events", {"SingleTop_Combined.root"});
  RDataFrame d_VV("Events", {"VV_Combined.root"});
  RDataFrame d_VVV("Events", {"VVV_Combined.root"});
  RDataFrame d_WPlusJets("Events", {"WPlusJets_Combined.root"});
  RDataFrame d_ttbarV("Events", {"ttbarV_Combined.root"});
  RDataFrame d_data("Events", {"data_Combined.root"});

  //Applying the event weight to the histograms
  auto h_tZq = d_tZq.Filter(variable_name).Histo1D({variable_name, "EventWeight"});
  auto h_ZPlusJets = d_ZPlusJets.Filter(variable_name).Histo1D({variable_name, "EventWeight"});
  auto h_SingleTop = d_SingleTop.Filter(variable_name).Histo1D({variable_name, "EventWeight"});
  auto h_VV = d_VV.Filter(variable_name).Histo1D({variable_name, "EventWeight"});
  auto h_VVV = d_VVV.Filter(variable_name).Histo1D({variable_name, "EventWeight"});
  auto h_WPlusJets = d_WPlusJets.Filter(variable_name).Histo1D({variable_name, "EventWeight"});
  auto h_ttbarV = d_ttbarV.Filter(variable_name).Histo1D({variable_name, "EventWeight"});
  auto h_data = d_data.Filter(variable_name).Histo1D({variable_name, "EventWeight"});

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
  h_VV->SetFillColor(802); //kOrange+2
  h_VVV->SetFillColor(400); //kYellow
  h_WPlusJets->SetFillColor(425); //kCyan-7
  h_ttbarV->SetFillColor(920); //kGray
  h_tZq->SetLineColor(600); //kBlue
  h_ZPlusJets->SetLineColor(632); //kRed
  h_ttbar->SetLineColor(418); //kGreen+2
  h_SingleTop->SetLineColor(618); //kMagenta+2
  h_VV->SetLineColor(802); //kOrange+2
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
  h_VV->SetTitle("VV");
  h_VVV->SetTitle("VVV");
  h_WPlusJets->SetTitle("W+jets");
  h_ttbarV->SetTitle("t#bar{t}V");
  h_data->SetTitle("data");

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
  MC_Stack->Add(h_ZPlusJets);
  MC_Stack->Add(h_ttbar);
  MC_Stack->Add(h_VV);
  MC_Stack->Add(h_WPlusJets);
  MC_Stack->Add(h_SingleTop);
  MC_Stack->Add(h_ttbarV);
  MC_Stack->Add(h_VVV);
  MC_Stack->Add(h_tZq);
  MC_Stack->Draw("HIST");
  h_data->Draw("HIST PSAME");

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
  rp->GetXaxis()->SetTitle(xaxis_name);
  rp->SetMinimum(0.5);
  rp->SetMaximum(1.5);
  rp->Draw();

  //Updating the canvasand saving the output pdf file
  c1->Update();
  c1->SaveAs(variable_name + "_" + channel+ "_" + systematic + "_" + region + ".pdf");


}



void process_comparison(){

  process_comparison_plotter("InvTopMass", "ee", "Mass [GeV]", "Nominal", "SBR", "Invariant mass of the top quark candidate (ee channel)");

}
