#include "TGraph.h"
#include "TROOT.h"

using namespace std;


float NanoAOD_MCEff_ee, NanoAOD_MCEff_mumu, NanoAOD_MCEff_emu;
float NanoAOD_MCEff_UpUncert_ee, NanoAOD_MCEff_UpUncert_mumu, NanoAOD_MCEff_UpUncert_emu;
float NanoAOD_MCEff_DownUncert_ee, NanoAOD_MCEff_DownUncert_mumu, NanoAOD_MCEff_DownUncert_emu;

float NanoAOD_DataEff_ee, NanoAOD_DataEff_mumu, NanoAOD_DataEff_emu;
float NanoAOD_DataEff_UpUncert_ee, NanoAOD_DataEff_UpUncert_mumu, NanoAOD_DataEff_UpUncert_emu;
float NanoAOD_DataEff_DownUncert_ee, NanoAOD_DataEff_DownUncert_mumu, NanoAOD_DataEff_DownUncert_emu;

float NanoAOD_Alpha_ee, NanoAOD_Alpha_mumu, NanoAOD_Alpha_emu;

float NanoAOD_TriggerSF_ee, NanoAOD_TriggerSF_mumu, NanoAOD_TriggerSF_emu;
float NanoAOD_TriggerSF_UpUncert_ee, NanoAOD_TriggerSF_UpUncert_mumu, NanoAOD_TriggerSF_UpUncert_emu;
float NanoAOD_TriggerSF_DownUncert_ee, NanoAOD_TriggerSF_DownUncert_mumu, NanoAOD_TriggerSF_DownUncert_emu;


float MiniAOD_MCEff_ee, MiniAOD_MCEff_mumu, MiniAOD_MCEff_emu;
float MiniAOD_MCEff_UpUncert_ee, MiniAOD_MCEff_UpUncert_mumu, MiniAOD_MCEff_UpUncert_emu;
float MiniAOD_MCEff_DownUncert_ee, MiniAOD_MCEff_DownUncert_mumu, MiniAOD_MCEff_DownUncert_emu;

float MiniAOD_DataEff_ee, MiniAOD_DataEff_mumu, MiniAOD_DataEff_emu;
float MiniAOD_DataEff_UpUncert_ee, MiniAOD_DataEff_UpUncert_mumu, MiniAOD_DataEff_UpUncert_emu;
float MiniAOD_DataEff_DownUncert_ee, MiniAOD_DataEff_DownUncert_mumu, MiniAOD_DataEff_DownUncert_emu;

float MiniAOD_Alpha_ee, MiniAOD_Alpha_mumu, MiniAOD_Alpha_emu;

float MiniAOD_TriggerSF_ee, MiniAOD_TriggerSF_mumu, MiniAOD_TriggerSF_emu;
float MiniAOD_TriggerSF_UpUncert_ee, MiniAOD_TriggerSF_UpUncert_mumu, MiniAOD_TriggerSF_UpUncert_emu;
float MiniAOD_TriggerSF_DownUncert_ee, MiniAOD_TriggerSF_DownUncert_mumu, MiniAOD_TriggerSF_DownUncert_emu;



void TriggerSFComparison_MiniAndNanoAOD2(const std::string& year){

  const int n = 3;
  double Channels[] = {1, 2, 3};

  if(year == "2016"){

	NanoAOD_MCEff_ee = 0.970543; NanoAOD_MCEff_mumu = 0.976452; NanoAOD_MCEff_emu = 0;
        NanoAOD_MCEff_UpUncert_ee = 0.008; NanoAOD_MCEff_UpUncert_mumu = 0.004; NanoAOD_MCEff_UpUncert_emu = 0;
        NanoAOD_MCEff_DownUncert_ee = -0.008; NanoAOD_MCEff_DownUncert_mumu = -0.004; NanoAOD_MCEff_DownUncert_emu = 0;

        NanoAOD_DataEff_ee = 0.972363; NanoAOD_DataEff_mumu = 0.959516; NanoAOD_DataEff_emu = 0;
        NanoAOD_DataEff_UpUncert_ee = 0.003; NanoAOD_DataEff_UpUncert_mumu = 0.007; NanoAOD_DataEff_UpUncert_emu = 0;
        NanoAOD_DataEff_DownUncert_ee = -0.002; NanoAOD_DataEff_DownUncert_mumu = -0.001; NanoAOD_DataEff_DownUncert_emu = 0;

        NanoAOD_Alpha_ee = 0.999999; NanoAOD_Alpha_mumu = 0.999999; NanoAOD_Alpha_emu = 0;

        NanoAOD_TriggerSF_ee = 1.00187; NanoAOD_TriggerSF_mumu = 0.982655; NanoAOD_TriggerSF_emu = 0;
        NanoAOD_TriggerSF_UpUncert_ee = 0.003; NanoAOD_TriggerSF_UpUncert_mumu = 0.001; NanoAOD_TriggerSF_UpUncert_emu = 0;
        NanoAOD_TriggerSF_DownUncert_ee = -0.0010; NanoAOD_TriggerSF_DownUncert_mumu = -0.001; NanoAOD_TriggerSF_DownUncert_emu = 0;


        MiniAOD_MCEff_ee = 0.988; MiniAOD_MCEff_mumu = 0.992; MiniAOD_MCEff_emu = 0.977;
        MiniAOD_MCEff_UpUncert_ee = 0.001; MiniAOD_MCEff_UpUncert_mumu = 0.001; MiniAOD_MCEff_UpUncert_emu = 0.005;
        MiniAOD_MCEff_DownUncert_ee = -0.001; MiniAOD_MCEff_DownUncert_mumu = -0.001; MiniAOD_MCEff_DownUncert_emu = -0.005;

        MiniAOD_DataEff_ee = 0.976; MiniAOD_DataEff_mumu = 0.985; MiniAOD_DataEff_emu = 0.964;
        MiniAOD_DataEff_UpUncert_ee = 0.001; MiniAOD_DataEff_UpUncert_mumu = 0.001; MiniAOD_DataEff_UpUncert_emu = 0.011;
        MiniAOD_DataEff_DownUncert_ee = -0.001; MiniAOD_DataEff_DownUncert_mumu = -0.001; MiniAOD_DataEff_DownUncert_emu = -0.011;

        MiniAOD_Alpha_ee = 0.99890; MiniAOD_Alpha_mumu = 1.00151; MiniAOD_Alpha_emu = 0.98775;

        MiniAOD_TriggerSF_ee = 0.987; MiniAOD_TriggerSF_mumu = 0.993; MiniAOD_TriggerSF_emu = 0.987;
        MiniAOD_TriggerSF_UpUncert_ee = 0.001; MiniAOD_TriggerSF_UpUncert_mumu = 0.000; MiniAOD_TriggerSF_UpUncert_emu = 0.007;
        MiniAOD_TriggerSF_DownUncert_ee = -0.001; MiniAOD_TriggerSF_DownUncert_mumu = -0.000; MiniAOD_TriggerSF_DownUncert_emu = -0.007;


  }
  else if(year == "2017"){

	NanoAOD_MCEff_ee = 0.974412; NanoAOD_MCEff_mumu = 0.973075; NanoAOD_MCEff_emu = 0;
        NanoAOD_MCEff_UpUncert_ee = -0.002; NanoAOD_MCEff_UpUncert_mumu = -0.002; NanoAOD_MCEff_UpUncert_emu = 0;
        NanoAOD_MCEff_DownUncert_ee = 0.003; NanoAOD_MCEff_DownUncert_mumu = 0.001; NanoAOD_MCEff_DownUncert_emu = 0;

        NanoAOD_DataEff_ee = 0.965362; NanoAOD_DataEff_mumu = 0.956032; NanoAOD_DataEff_emu = 0;
        NanoAOD_DataEff_UpUncert_ee = -0.001; NanoAOD_DataEff_UpUncert_mumu = 0.001; NanoAOD_DataEff_UpUncert_emu = 0;
        NanoAOD_DataEff_DownUncert_ee = 0.005; NanoAOD_DataEff_DownUncert_mumu = 0.001; NanoAOD_DataEff_DownUncert_emu = 0;

        NanoAOD_Alpha_ee = 0.9999; NanoAOD_Alpha_mumu = 0.9999; NanoAOD_Alpha_emu = 0;

        NanoAOD_TriggerSF_ee = 0.990711; NanoAOD_TriggerSF_mumu = 0.982484; NanoAOD_TriggerSF_emu = 0;
        NanoAOD_TriggerSF_UpUncert_ee = -0.001; NanoAOD_TriggerSF_UpUncert_mumu = 0.001; NanoAOD_TriggerSF_UpUncert_emu = 0;
        NanoAOD_TriggerSF_DownUncert_ee = 0.005; NanoAOD_TriggerSF_DownUncert_mumu = 0.001; NanoAOD_TriggerSF_DownUncert_emu = 0;


        MiniAOD_MCEff_ee = 0.92849; MiniAOD_MCEff_mumu = 0.99509; MiniAOD_MCEff_emu = 0.91272;
        MiniAOD_MCEff_UpUncert_ee = -0.00041; MiniAOD_MCEff_UpUncert_mumu = -0.00008; MiniAOD_MCEff_UpUncert_emu = -0.00244;
        MiniAOD_MCEff_DownUncert_ee = 0.00041; MiniAOD_MCEff_DownUncert_mumu = 0.00008; MiniAOD_MCEff_DownUncert_emu = 0.00251;

        MiniAOD_DataEff_ee = 0.86448; MiniAOD_DataEff_mumu = 0.96692; MiniAOD_DataEff_emu = 0.87028;
        MiniAOD_DataEff_UpUncert_ee = -0.00277; MiniAOD_DataEff_UpUncert_mumu = -0.00035; MiniAOD_DataEff_UpUncert_emu = -0.00709;
        MiniAOD_DataEff_DownUncert_ee = 0.00282; MiniAOD_DataEff_DownUncert_mumu = 0.00036; MiniAOD_DataEff_DownUncert_emu = 0.00747;

        MiniAOD_Alpha_ee = 0.99695; MiniAOD_Alpha_mumu = 0.99982; MiniAOD_Alpha_emu = 0.99144;

        MiniAOD_TriggerSF_ee = 0.93106; MiniAOD_TriggerSF_mumu = 0.97170; MiniAOD_TriggerSF_emu = 0.95350;
        MiniAOD_TriggerSF_UpUncert_ee = 0.001; MiniAOD_TriggerSF_UpUncert_mumu = 0.001; MiniAOD_TriggerSF_UpUncert_emu = 0.001;
        MiniAOD_TriggerSF_DownUncert_ee = 0.001; MiniAOD_TriggerSF_DownUncert_mumu = 0.001; MiniAOD_TriggerSF_DownUncert_emu = 0.001;


  }
  else{std::cout << "ERROR: year must be 2016 or 2017" << std::endl;}


  double NanoAOD_YValues_MCEff[] = {NanoAOD_MCEff_ee, NanoAOD_MCEff_mumu, NanoAOD_MCEff_emu};
  double NanoAOD_YValues_ErrorsDown_MCEff[] = {NanoAOD_MCEff_DownUncert_ee, NanoAOD_MCEff_DownUncert_mumu, NanoAOD_MCEff_DownUncert_emu};
  double NanoAOD_YValues_ErrorsUp_MCEff[] = {NanoAOD_MCEff_UpUncert_ee, NanoAOD_MCEff_UpUncert_mumu, NanoAOD_MCEff_UpUncert_emu};

  double MiniAOD_YValues_MCEff[] = {MiniAOD_MCEff_ee, MiniAOD_MCEff_mumu, MiniAOD_MCEff_emu};
  double MiniAOD_YValues_ErrorsDown_MCEff[] = {MiniAOD_MCEff_DownUncert_ee, MiniAOD_MCEff_DownUncert_mumu, MiniAOD_MCEff_DownUncert_emu};
  double MiniAOD_YValues_ErrorsUp_MCEff[] = {MiniAOD_MCEff_UpUncert_ee, MiniAOD_MCEff_UpUncert_mumu, MiniAOD_MCEff_UpUncert_emu};

  double NanoAOD_YValues_DataEff[] = {NanoAOD_DataEff_ee, NanoAOD_DataEff_mumu, NanoAOD_DataEff_emu};
  double NanoAOD_YValues_ErrorsDown_DataEff[] = {NanoAOD_DataEff_DownUncert_ee, NanoAOD_DataEff_DownUncert_mumu, NanoAOD_DataEff_DownUncert_emu};
  double NanoAOD_YValues_ErrorsUp_DataEff[] = {NanoAOD_DataEff_UpUncert_ee, NanoAOD_DataEff_UpUncert_mumu, NanoAOD_DataEff_UpUncert_emu};

  double MiniAOD_YValues_DataEff[] = {MiniAOD_DataEff_ee, MiniAOD_DataEff_mumu, MiniAOD_DataEff_emu};
  double MiniAOD_YValues_ErrorsDown_DataEff[] = {MiniAOD_DataEff_DownUncert_ee, MiniAOD_DataEff_DownUncert_mumu, MiniAOD_DataEff_DownUncert_emu};
  double MiniAOD_YValues_ErrorsUp_DataEff[] = {MiniAOD_DataEff_UpUncert_ee, MiniAOD_DataEff_UpUncert_mumu, MiniAOD_DataEff_UpUncert_emu};
  
  double NanoAOD_YValues_Sf[] = {NanoAOD_TriggerSF_ee, NanoAOD_TriggerSF_mumu, NanoAOD_TriggerSF_emu};
  double NanoAOD_YValues_ErrorsDown_Sf[] = {NanoAOD_TriggerSF_DownUncert_ee, NanoAOD_TriggerSF_DownUncert_mumu, NanoAOD_TriggerSF_DownUncert_emu};
  double NanoAOD_YValues_ErrorsUp_Sf[] = {NanoAOD_TriggerSF_UpUncert_ee, NanoAOD_TriggerSF_UpUncert_mumu, NanoAOD_TriggerSF_UpUncert_emu};

  double MiniAOD_YValues_Sf[] = {MiniAOD_TriggerSF_ee, MiniAOD_TriggerSF_mumu, MiniAOD_TriggerSF_emu};
  double MiniAOD_YValues_ErrorsDown_Sf[] = {MiniAOD_TriggerSF_DownUncert_ee, MiniAOD_TriggerSF_DownUncert_mumu, MiniAOD_TriggerSF_DownUncert_emu};
  double MiniAOD_YValues_ErrorsUp_Sf[] = {MiniAOD_TriggerSF_UpUncert_ee, MiniAOD_TriggerSF_UpUncert_mumu, MiniAOD_TriggerSF_UpUncert_emu};

  double NanoAOD_YValues_Alpha[] = {NanoAOD_Alpha_ee, NanoAOD_Alpha_mumu, NanoAOD_Alpha_emu};
  double NanoAOD_YValues_ErrorsDown_Alpha[] = {0, 0, 0};
  double NanoAOD_YValues_ErrorsUp_Alpha[] = {0, 0, 0};

  double MiniAOD_YValues_Alpha[] = {MiniAOD_Alpha_ee, MiniAOD_Alpha_mumu, MiniAOD_Alpha_emu};
  double MiniAOD_YValues_ErrorsDown_Alpha[] = {0, 0, 0};
  double MiniAOD_YValues_ErrorsUp_Alpha[] = {0, 0, 0};


  TGraphAsymmErrors* h_NanoAOD_MCEff = new TGraphAsymmErrors(n, Channels, NanoAOD_YValues_MCEff, 0, 0, NanoAOD_YValues_ErrorsDown_MCEff, NanoAOD_YValues_ErrorsUp_MCEff);
  TGraphAsymmErrors* h_MiniAOD_MCEff = new TGraphAsymmErrors(n, Channels, MiniAOD_YValues_MCEff, 0, 0, MiniAOD_YValues_ErrorsDown_MCEff, MiniAOD_YValues_ErrorsUp_MCEff);

  TCanvas* c1 =new TCanvas("c1"," ",200,10,700,500);
  c1->cd();

  TMultiGraph *gr = new TMultiGraph();
  std::string MCEffTitlesString = "MC Efficiency Values " + year;
  std::string MCEffValuesPdfString = "MCEffValues_" + year + ".pdf"; 

  gr->SetTitle(MCEffTitlesString.c_str());

  h_NanoAOD_MCEff->SetMarkerColor(4);
  h_NanoAOD_MCEff->SetMarkerStyle(8);
  h_MiniAOD_MCEff->SetMarkerColor(1);
  h_MiniAOD_MCEff->SetMarkerStyle(8);
  h_NanoAOD_MCEff->SetTitle("Nano AOD");
  h_MiniAOD_MCEff->SetTitle("Mini AOD");

  gr->Add(h_NanoAOD_MCEff);
  gr->Add(h_MiniAOD_MCEff);
  gr->GetYaxis()->SetTitle("Efficiency");
  gr->GetXaxis()->SetTitle("Channel");
  gr->GetXaxis()->SetTitleOffset(1.5);
  gr->GetXaxis()->SetBinLabel(1, "ee");
  gr->GetXaxis()->SetBinLabel(50, "mumu");
  gr->GetXaxis()->SetBinLabel(100, "emu");

  gr->Draw("AP");
  c1->BuildLegend();
  c1->SetGridx(1);
  c1->SetGridy(1);
  c1->SaveAs(MCEffValuesPdfString.c_str());



  TGraphAsymmErrors* h_NanoAOD_DataEff = new TGraphAsymmErrors(n, Channels, NanoAOD_YValues_DataEff, 0, 0, NanoAOD_YValues_ErrorsDown_DataEff, NanoAOD_YValues_ErrorsUp_DataEff);
  TGraphAsymmErrors* h_MiniAOD_DataEff = new TGraphAsymmErrors(n, Channels, MiniAOD_YValues_DataEff, 0, 0, MiniAOD_YValues_ErrorsDown_DataEff, MiniAOD_YValues_ErrorsUp_DataEff);

  TCanvas* c2 =new TCanvas("c2"," ",200,10,700,500);
  c2->cd();

  TMultiGraph *gr2 = new TMultiGraph();
  std::string DataEffTitlesString = "Data Efficiency Values " + year;
  std::string DataEffValuesPdfString = "DataEffValues_" + year + ".pdf";

  gr2->SetTitle(DataEffTitlesString.c_str());

  h_NanoAOD_DataEff->SetMarkerColor(4);
  h_NanoAOD_DataEff->SetMarkerStyle(8);
  h_MiniAOD_DataEff->SetMarkerColor(1);
  h_MiniAOD_DataEff->SetMarkerStyle(8);
  h_NanoAOD_DataEff->SetTitle("Nano AOD");
  h_MiniAOD_DataEff->SetTitle("Mini AOD");

  gr2->Add(h_NanoAOD_DataEff);
  gr2->Add(h_MiniAOD_DataEff);
  gr2->GetYaxis()->SetTitle("Efficiency");
  gr2->GetXaxis()->SetTitle("Channel");
  gr2->GetXaxis()->SetTitleOffset(1.5);
  gr2->GetXaxis()->SetBinLabel(1, "ee");
  gr2->GetXaxis()->SetBinLabel(50, "mumu");
  gr2->GetXaxis()->SetBinLabel(100, "emu");
  gr2->Draw("AP");

  c2->BuildLegend(0.2,0.2,0.5,0.4);
  c2->SetGridx(1);
  c2->SetGridy(1);
  c2->SaveAs(DataEffValuesPdfString.c_str());


  TGraphAsymmErrors* h_NanoAOD_Sf = new TGraphAsymmErrors(n, Channels, NanoAOD_YValues_Sf, 0, 0, NanoAOD_YValues_ErrorsDown_Sf, NanoAOD_YValues_ErrorsUp_Sf);
  TGraphAsymmErrors* h_MiniAOD_Sf = new TGraphAsymmErrors(n, Channels, MiniAOD_YValues_Sf, 0, 0, MiniAOD_YValues_ErrorsDown_Sf, MiniAOD_YValues_ErrorsUp_Sf);

  TCanvas* c3 =new TCanvas("c3"," ",200,10,700,500);
  c3->cd();

  TMultiGraph *gr3 = new TMultiGraph();
  std::string SFValuesTitle = "Scale Factor Values " + year;
  std::string SFValuesPdf = "SFValues_" + year + ".pdf";

  gr3->SetTitle("Scale Factor Values");

  h_NanoAOD_Sf->SetMarkerColor(4);
  h_NanoAOD_Sf->SetMarkerStyle(8);
  h_MiniAOD_Sf->SetMarkerColor(1);
  h_MiniAOD_Sf->SetMarkerStyle(8);
  h_NanoAOD_Sf->SetTitle("Nano AOD");
  h_MiniAOD_Sf->SetTitle("Mini AOD");

  gr3->Add(h_NanoAOD_Sf);
  gr3->Add(h_MiniAOD_Sf);
  gr3->GetYaxis()->SetTitle("Scale Factor");
  gr3->GetXaxis()->SetTitle("Channel");
  gr3->GetXaxis()->SetTitleOffset(1.5);
  gr3->GetXaxis()->SetBinLabel(1, "ee");
  gr3->GetXaxis()->SetBinLabel(50, "mumu");
  gr3->GetXaxis()->SetBinLabel(100, "emu");

  gr3->Draw("AP");
  c3->BuildLegend(0.6,0.1, 0.9,0.3);
  c3->SetGridx(1);
  c3->SetGridy(1);
  c3->SaveAs(SFValuesPdf.c_str());



  TGraphAsymmErrors* h_NanoAOD_Alpha = new TGraphAsymmErrors(n, Channels, NanoAOD_YValues_Alpha, 0, 0, NanoAOD_YValues_ErrorsDown_Alpha, NanoAOD_YValues_ErrorsUp_Alpha);
  TGraphAsymmErrors* h_MiniAOD_Alpha = new TGraphAsymmErrors(n, Channels, MiniAOD_YValues_Alpha, 0, 0, MiniAOD_YValues_ErrorsDown_Alpha, MiniAOD_YValues_ErrorsUp_Alpha);

  TCanvas* c4 =new TCanvas("c4"," ",200,10,700,500);
  c4->cd();

  TMultiGraph *gr4 = new TMultiGraph();
  std::string AlphaValuesTitle = "Alpha Values " + year;
  std::string AlphaValuesPdf = "AlphaValues_" + year + ".pdf";

  gr4->SetTitle(AlphaValuesTitle.c_str());

  h_NanoAOD_Alpha->SetMarkerColor(4);
  h_NanoAOD_Alpha->SetMarkerStyle(8);
  h_MiniAOD_Alpha->SetMarkerColor(1);
  h_MiniAOD_Alpha->SetMarkerStyle(8);
  h_NanoAOD_Alpha->SetTitle("Nano AOD");
  h_MiniAOD_Alpha->SetTitle("Mini AOD");

  gr4->Add(h_NanoAOD_Alpha);
  gr4->Add(h_MiniAOD_Alpha);
  gr4->GetYaxis()->SetTitle("Alpha");
  gr4->GetXaxis()->SetTitle("Channel");
  gr4->GetXaxis()->SetTitleOffset(1.5);
  gr4->GetXaxis()->SetBinLabel(1, "ee");
  gr4->GetXaxis()->SetBinLabel(50, "mumu");
  gr4->GetXaxis()->SetBinLabel(100, "emu");

  gr4->Draw("AP");
  c4->BuildLegend();
  c4->SetGridx(1);
  c4->SetGridy(1);
  c4->SaveAs(AlphaValuesPdf.c_str());


}





void TriggerSFComparison_MiniAndNanoAOD(){

    TriggerSFComparison_MiniAndNanoAOD2("2016");
    TriggerSFComparison_MiniAndNanoAOD2("2017");
  
}
 
