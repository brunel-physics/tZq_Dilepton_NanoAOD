//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//
// Script to plot background vs signal simulations for tZq analysis
//
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <algorithm>
#include <TFile.h>
#include "TH1D.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"

using namespace std;
using namespace ROOT; // RDataFrame's namespace

//////////////////////////////////////////////////////////////////////////
// Class for variables under consideration
//////////////////////////////////////////////////////////////////////////

class variable
{
    public:
        string name;
        string title;
        string graph_title;
        double xlow, xhigh;
        int nbins;
        bool log_scale;

        void set_range(double, double);
};

void variable::set_range(double _xlow, double _xhigh)
{
    xlow = _xlow;
    xhigh = _xhigh;
}

//////////////////////////////////////////////////////////////////////////
// Define the variables
//////////////////////////////////////////////////////////////////////////

vector <variable> define_variables()
{
    vector <variable> variables;
    variable var;
    
    // Generic
    var.nbins = 25;
    var.log_scale = true;
    
    var.set_range(0, 1000);
    
    // Invariant top mass
    var.name = "InvTopMass";
    var.title = "m_{t} [GeV]";
    var.graph_title = "Invariant mass of the top quark candidate";
    variables.push_back(var);
    
    // W mass
    var.name = "w_mass";
    var.title = "m_{W} [GeV]";
    var.graph_title = "Reconstructed mass of the W boson candidate";
    variables.push_back(var);
    
    // Z mass
    var.name = "z_mass";
    var.title = "m_{Z} [GeV]";
    var.graph_title = "Reconstructed mass of the Z boson candidate";
    variables.push_back(var);
    
    var.set_range(0, 20);
    
    // Leading lepton mass
    var.name = "LeadingLeptonMass";
    var.title = "m_{lep1} [GeV]";
    var.graph_title = "Mass of the leading lepton";
    variables.push_back(var);
    
    // Subleading lepton mass
    var.name = "SubleadingLeptonMass";
    var.title = "m_{lep2} [GeV]";
    var.graph_title = "Mass of the subleading lepton";
    variables.push_back(var);
    
    var.set_range(0, 200);
    
    // Smeared jet mass
    var.name = "SmearedJetMass";
    var.title = "m_{smeared jets} [GeV]";
    var.graph_title = "Mass of smeared jets";
    variables.push_back(var);
    
    // Tight smeared jets mass
    var.name = "TightSmearedJetsMass";
    var.title = "m_{tight smeared jets} [GeV]";
    var.graph_title = "Mass of tight smeared jets";
    variables.push_back(var);
    
    // Leading jet mass
    var.name = "LeadingJetMass";
    var.title = "m_{jet1} [GeV]";
    var.graph_title = "Mass of the leading tight smeared jet";
    variables.push_back(var);
    
    // Subleading jet mass
    var.name = "SubleadingJetMass";
    var.title = "m_{jet2} [GeV]";
    var.graph_title = "Mass of the subleading tight smeared jet";
    variables.push_back(var);
    
    // Third jet mass
    var.name = "ThirdJetMass";
    var.title = "m_{jet3} [GeV]";
    var.graph_title = "Mass of the third tight smeared jet";
    variables.push_back(var);
    
    // Fourth jet mass
    var.name = "FourthJetMass";
    var.title = "m_{jet4} [GeV]";
    var.graph_title = "Mass of the fourth tight smeared jet";
    variables.push_back(var);
    
    return variables;
}

//////////////////////////////////////////////////////////////////////////
// Class for processes under consideration
//////////////////////////////////////////////////////////////////////////

class process
{
    public:
        string name;
        string title;
        string draw_mode;
        
        // Marker properties
        int marker_style;
        int marker_size;
        int marker_colour;
        
        // Line properties
        int line_width;
        int line_colour;
        
        // Fill properties
        int fill_colour;
};

//////////////////////////////////////////////////////////////////////////
// Define the processes
//////////////////////////////////////////////////////////////////////////

vector <process> define_processes(bool include_signal, bool include_data)
{
    vector <process> processes;
    process proc;
    
    // Generic
    proc.draw_mode = "f";
    proc.line_width = 1;
    proc.line_colour = 1; // kBlack
    
    // Z+jets
    proc.name = "ZPlusJets";
    proc.title = "Z+jets";
    proc.fill_colour = 632; // kRed
    processes.push_back(proc);
    
    // ttbar
    proc.name = "ttbar";
    proc.title = "t#bar{t}";
    proc.fill_colour = 418; // kGreen+2
    processes.push_back(proc);
    
    // Single top
    proc.name = "SingleTop";
    proc.title = "Single top";
    proc.fill_colour = 618; // kMagenta+2
    processes.push_back(proc);
    
    // Diboson
    proc.name = "VV";
    proc.title = "Diboson";
    proc.fill_colour = 802; // kOrange+2
    //processes.push_back(proc);
    
    // Triboson
    proc.name = "VVV";
    proc.title = "Triboson";
    proc.fill_colour = 400; // kYellow
    processes.push_back(proc);
    
    // W+jets
    proc.name = "WPlusJets";
    proc.title = "W+jets";
    proc.fill_colour = 425; // kCyan-7
    processes.push_back(proc);
    
    // ttbar+V
    proc.name = "ttbarV";
    proc.title = "t#bar{t}+V";
    proc.fill_colour = 920; // kGray
    processes.push_back(proc);
    
    // tZq
    if(include_signal)
    {
        proc.name = "tZq";
        proc.title = "tZq";
        proc.draw_mode = "l";
        proc.line_width = 2;
        proc.line_colour = 600; // kBlue
        proc.fill_colour = 0; // Hollow
        processes.push_back(proc);
    }

    // Data
    if(include_data)
    {
        proc.name = "data";
        proc.title = "Data";
        proc.draw_mode = "pe";
        proc.marker_style = 20;
        proc.marker_size = 2;
        proc.marker_colour = 1; // kBlack
        processes.push_back(proc);
    }

    return processes;
}

//////////////////////////////////////////////////////////////////////////
// Plotter
//////////////////////////////////////////////////////////////////////////

void plotter(variable var, string year, string channel, string systematic, string region, bool include_signal, bool include_data, bool normalise)
{

    // Define the included processes
    vector <process> processes = define_processes(include_signal, include_data);
    
    // Form histograms for each process
    map <string, ROOT::RDF::RResultPtr<TH1D>> hist;
    for(const auto & proc: processes)
    {
        string file = proc.name + "_" + systematic + "_" + channel + "_" + region + "_" + year + "_Combined.root";
        RDataFrame data_frame("Events", file.c_str());
        hist[proc.name] = data_frame.Histo1D({proc.name.c_str(), var.name.c_str(), var.nbins, var.xlow, var.xhigh}, var.name.c_str(), "EventWeight");
    }
    
    // Set up the canvas
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TCanvas * c1 = new TCanvas("c1", "canvas", 960, 800);
    double pad1_ylow;
    pad1_ylow = include_data ? 0.25 : 0;
    TPad * pad1 = new TPad("pad1", "pad1", 0, pad1_ylow, 1, 1);
    if(var.log_scale) pad1->SetLogy();
    c1->cd();
    pad1->Draw();
    pad1->cd();
    
    // Format the processes, and add the backgrounds to a stack
    THStack * bg_stack = new THStack("bg_stack", var.graph_title.c_str());
    double total_integral = 0.0;
    for(auto proc = processes.rbegin(); proc != processes.rend(); ++proc)
    {
        // Format
        hist[proc->name]->SetMarkerStyle(proc->marker_style);
        hist[proc->name]->SetMarkerSize(proc->marker_size);
        hist[proc->name]->SetMarkerColor(proc->marker_colour);
        hist[proc->name]->SetLineWidth(proc->line_width);
        hist[proc->name]->SetLineColor(proc->line_colour);
        hist[proc->name]->SetFillColor(proc->fill_colour);

        // Add backgrounds to stack
        if((proc->name != "tZq") && (proc->name != "data"))
        {
            bg_stack->Add(hist[proc->name].GetPtr());
            total_integral += hist[proc->name]->Integral("width");
        }
    }
    
    // Normalise if requested
    if(normalise)
    {
        hist["tZq"]->Scale(1/hist["tZq"]->Integral("width"));
        for(const auto & proc: processes)
            if((proc.name != "tZq") && (proc.name != "data"))
                hist[proc.name]->Scale(1/total_integral);
    }
    
    // Draw the components
    bg_stack->Draw("hist");
    if(include_signal) hist["tZq"]->Draw("hist same");
    if(include_data) hist["data"]->Draw("pe same");
    
    // Axis titles
    //bg_stack->SetMinimum(1e-1);
    bg_stack->GetXaxis()->SetTitle(var.title.c_str());
    bg_stack->GetYaxis()->SetTitle("Events per bin");
    
    // Legend
    TLegend * leg = new TLegend(0.65, 0.6, 0.89, 0.89);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    if(include_data) leg->AddEntry(hist["data"].GetPtr(), "Data", "pe");
    if(include_signal) leg->AddEntry(hist["tZq"].GetPtr(), "tZq (signal)", "l");
    for(const auto & proc: processes)
        if((proc.name != "tZq") && (proc.name != "data"))
            leg->AddEntry(hist[proc.name].GetPtr(), proc.title.c_str(), "f");
    leg->Draw();

    // CMS text
    TPaveText * CMS = new TPaveText(0.1, 0.88, 0.2, 0.98, "blNDC");
    TText * CMS_text = CMS->AddText("CMS");
    CMS_text->SetTextFont(61);
    CMS->SetFillStyle(0);
    CMS->SetBorderSize(0);
    CMS->Draw();
    
    // Work in progress text
    TPaveText * WIP = new TPaveText(0.21, 0.88, 0.41, 0.95, "blNDC");
    TText * WIP_text = WIP->AddText("Work in progress");
    WIP_text->SetTextFont(52);
    WIP->SetFillStyle(0);
    WIP->SetBorderSize(0);
    WIP->Draw();
    
    // Lumi text
    string lumi_str;
    switch(stoi(year))
    {
        case 2016:
            lumi_str = "35.9 fb^{-1} (13 TeV)"; break;
        case 2017:
            lumi_str = "41.9 fb^{-1} (13 TeV)"; break;
        case 2018:
            lumi_str = "59.7 fb^{-1} (13 TeV)"; break;
    }
    TPaveText * lumi = new TPaveText(0.65, 0.875, 0.9, 0.975, "blNDC");
    TText * lumi_text = lumi->AddText(lumi_str.c_str());
    lumi->SetFillStyle(0);
    lumi->SetBorderSize(0);
    lumi->Draw();
    
    // Channel text
    TPaveText * ch = new TPaveText(0.4, 0.84, 0.57, 0.9, "blNDC");
    if(channel == "ee")
        ch->AddText("ee channel");
    else if(channel == "mumu")
        ch->AddText("#mu#mu channel");
    else
        ch->AddText("e#mu channel");
    ch->SetFillStyle(0);
    ch->SetBorderSize(0);
    ch->Draw();

    // Output
    string norm_str = normalise ? "_norm" : "";
    string outfile = "Plots/" + var.name + norm_str + "_" + channel + "_" + systematic + "_" + region + "_" + year + ".pdf";
    c1->Update();
    c1->SaveAs(outfile.c_str());
    c1->Delete();

}

//////////////////////////////////////////////////////////////////////////
// Main
//////////////////////////////////////////////////////////////////////////

void process_comparison(){

    // Make the Plots directory if it does not exist
    gSystem->Exec("mkdir -p Plots");

    // Define the variables
    vector <variable> variables = define_variables();
    
    // Define years, channels, regions
    vector <string> years{"2016"};
    vector <string> channels{"ee"};
    vector <string> systematics{"Nominal"};
    vector <string> regions{"NoChi2Cut"};
    
    // Make the plots
    for(const auto & year: years)
        for(const auto & channel: channels)
            for(const auto & systematic: systematics)
                for(const auto & region: regions)
                    for(const auto & var: variables)
                        plotter(var, year, channel, systematic, region, true, false, false);
}
