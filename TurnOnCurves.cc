#include <TFile.h>
#include <iostream>
#include <fstream>
#include "TH1D.h"


void TurnOnCurvesProducer(const int& YearInt, const int& ChannelInt, const string& VariableName, const string& WhichLepton){

  std::string ChannelString;
  std::string ChannelStringSymbol;
  std::string xaxis_name;
  std::string GraphTitle;
  std::string MC_File_String;
  std::string Data_File_String;
  std::string YearString;  

  std::string Histo_String = "h_" + WhichLepton + "Lepton" + VariableName;

  switch(ChannelInt){

	case 1: ChannelString = "ee"; ChannelStringSymbol = "ee"; break;
	case 2: ChannelString = "mumu"; ChannelStringSymbol = "#mu#mu"; break;
	case 3: ChannelString = "emu"; ChannelStringSymbol = "e#mu"; break;

	default: throw std::logic_error("ERROR: ChannelInt must be 1 (for ee), 2 (for mumu) or 3 (for emu)"); break;

  }

  if(VariableName != "Pt" && VariableName != "Eta"){throw std::logic_error("ERROR: VariableName must be Pt or Eta");}

  if(VariableName == "Pt"){xaxis_name = "p_{T}";}
  else{xaxis_name = "#eta";}

  switch(YearInt){

	case 2016: MC_File_String = "TurnOnCurves_TriggerSF_MC_Nominal_" + ChannelString + "__SR_SBR___2016.root";
		   Data_File_String = "TurnOnCurves_TriggerSF_DATA_Nominal_" + ChannelString + "__SR_SBR___2016.root";
		   YearString = "2016";
	           
                   break;

	case 2017: MC_File_String = "TurnOnCurves_TriggerSF_MC_Nominal_" + ChannelString + "__SR_SBR___2017.root";
                   Data_File_String = "TurnOnCurves_TriggerSF_DATA_Nominal_" + ChannelString + "__SR_SBR___2017.root";
		   YearString = "2017";
			
		   break;

	case 2018: MC_File_String = "TurnOnCurves_TriggerSF_MC_Nominal_" + ChannelString + "__SR_SBR___2018.root";
                   Data_File_String = "TurnOnCurves_TriggerSF_DATA_Nominal_" + ChannelString + "__SR_SBR___2018.root";
		   YearString = "2018";

		   break;

	default: throw std::logic_error("ERROR: Year must be 2016, 2017 or 2018"); break;
	
  }

  if(WhichLepton != "Leading" && WhichLepton != "Subleading"){throw std::logic_error("ERROR: WhichLepton must be Leading or Subleading");}

  if(WhichLepton == "Leading"){GraphTitle = "Leading lepton " + xaxis_name + " (" + ChannelStringSymbol + ", " + YearString + ")";}
  else{GraphTitle = "Subleading lepton " + xaxis_name + " (" + ChannelStringSymbol + ", " + YearString + ")";}


  TFile *MCFile = TFile::Open(MC_File_String.c_str(), "READ");
  TFile *DataFile = TFile::Open(Data_File_String.c_str(), "READ");

  TH1* h_MC = (TH1*)MCFile->GetObjectChecked(Histo_String.c_str(), "TH1");
  TH1* h_Data = (TH1*)DataFile->GetObjectChecked(Histo_String.c_str(), "TH1");

  int W = 800;
  int H = 600;
  int H_ref = 600;
  int W_ref = 800;

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
  gPad->BuildLegend();

  h_MC->SetMarkerColor(2); //kRed
  h_MC->SetMarkerStyle(20);
  h_MC->SetMarkerSize(1.0);
  h_MC->GetYaxis()->SetTitle("Efficiency");
  h_Data->SetMarkerColor(1); //kBlack
  h_Data->SetMarkerStyle(20);
  h_Data->SetMarkerSize(1.0);

  h_MC->Draw("HIST");
  h_Data->Draw("HIST PSAME");


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

  TH1D * rp = (TH1D*)(h_Data->Clone());
  rp->Divide((TH1D*)(h_MC->Clone()));

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
  c1->Update();

  std::string OutputFileName = "TurnOnCurves_" + WhichLepton + VariableName + "_" + ChannelString + "_" + YearString + ".pdf";
  c1->SaveAs(OutputFileName.c_str());

}




void TurnOnCurves(){

  //2016
  TurnOnCurvesProducer(2016, 1, "Pt", "Leading");
  TurnOnCurvesProducer(2016, 2, "Pt", "Leading");
  //TurnOnCurvesProducer(2016, 3, "Pt", "Leading");

  TurnOnCurvesProducer(2016, 1, "Pt", "Subleading");
  TurnOnCurvesProducer(2016, 2, "Pt", "Subleading");
  //TurnOnCurvesProducer(2016, 3, "Pt", "Subleading");

  TurnOnCurvesProducer(2016, 1, "Eta", "Leading");
  TurnOnCurvesProducer(2016, 2, "Eta", "Leading");
  //TurnOnCurvesProducer(2016, 3, "Eta", "Leading");

  TurnOnCurvesProducer(2016, 1, "Eta", "Subleading");
  TurnOnCurvesProducer(2016, 2, "Eta", "Subleading");
  //TurnOnCurvesProducer(2016, 3, "Eta", "Subleading");

  //2017
  TurnOnCurvesProducer(2017, 1, "Pt", "Leading");
  TurnOnCurvesProducer(2017, 2, "Pt", "Leading");
  //TurnOnCurvesProducer(2017, 3, "Pt", "Leading");

  TurnOnCurvesProducer(2017, 1, "Pt", "Subleading");
  TurnOnCurvesProducer(2017, 2, "Pt", "Subleading");
  //TurnOnCurvesProducer(2017, 3, "Pt", "Subleading");

  TurnOnCurvesProducer(2017, 1, "Eta", "Leading");
  TurnOnCurvesProducer(2017, 2, "Eta", "Leading");
  //TurnOnCurvesProducer(2017, 3, "Eta", "Leading");

  TurnOnCurvesProducer(2017, 1, "Eta", "Subleading");
  TurnOnCurvesProducer(2017, 2, "Eta", "Subleading");
  //TurnOnCurvesProducer(2017, 3, "Eta", "Subleading");

  //2018
  TurnOnCurvesProducer(2018, 1, "Pt", "Leading");
  TurnOnCurvesProducer(2018, 2, "Pt", "Leading");
  //TurnOnCurvesProducer(2018, 3, "Pt", "Leading");

  TurnOnCurvesProducer(2018, 1, "Pt", "Subleading");
  TurnOnCurvesProducer(2018, 2, "Pt", "Subleading");
  //TurnOnCurvesProducer(2018, 3, "Pt", "Subleading");

  TurnOnCurvesProducer(2018, 1, "Eta", "Leading");
  TurnOnCurvesProducer(2018, 2, "Eta", "Leading");
  //TurnOnCurvesProducer(2018, 3, "Eta", "Leading");

  TurnOnCurvesProducer(2018, 1, "Eta", "Subleading");
  TurnOnCurvesProducer(2018, 2, "Eta", "Subleading");
  //TurnOnCurvesProducer(2018, 3, "Eta", "Subleading");

}
