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
    var.log_scale = false;
   /* 
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
   
    var.set_range(0, 250);
 
    // Z mass
    var.name = "z_mass";
    var.title = "m_{Z} [GeV]";
    var.graph_title = "Reconstructed mass of the Z boson candidate";
    variables.push_back(var);
    
    var.set_range(0, 100);
    
    // Leading lepton mass
    var.name = "LeadingLeptonMass";
    var.title = "m_{lep1} [GeV]";
    var.graph_title = "Mass of the leading lepton";
    variables.push_back(var);    

    var.set_range(0, 0.15);

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
   
    var.set_range(0, 100);
 
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
   
    var.set_range(0, 20);
 
    // Third jet mass
    var.name = "ThirdJetMass";
    var.title = "m_{jet3} [GeV]";
    var.graph_title = "Mass of the third tight smeared jet";
    variables.push_back(var);
   
    var.set_range(0, 20);
 
    // Fourth jet mass
    var.name = "FourthJetMass";
    var.title = "m_{jet4} [GeV]";
    var.graph_title = "Mass of the fourth tight smeared jet";
    variables.push_back(var);

    //Lepton pt
    var.name = "LeptonPt";
    var.title = "p_{T} [GeV]";
    var.graph_title = "Transverse momenta of all leptons";
    variables.push_back(var);

    var.set_range(-4, 4);

    //Lepton phi
    var.name = "LeptonPhi";
    var.title = "#phi";
    var.graph_title = "Pseudorapidity of all leptons";
    variables.push_back(var);

    var.set_range(-5, 5);   

    //Lepton eta  
    var.name = "LeptonEta";
    var.title = "#eta";
    var.graph_title = "#eta of all leptons";
    variables.push_back(var);  

    var.set_range(-2, 2);

    //Lepton charge
    var.name = "LeptonCharge";
    var.title = "charge";
    var.graph_title = "Charge of all leptons";
    variables.push_back(var);

    var.set_range(0, 50);

    //Lepton mass
    var.name = "LeptonMass";
    var.title = "m_{l} [GeV]";
    var.graph_title = "Mass of all leptons";
    variables.push_back(var);
   
    var.set_range(0, 1);

    //Lepton jet relative isolation
    var.name = "LeptonJetRelIso";
    var.title = "Lepton-jet relative isolation";
    var.graph_title = "Lepton-jet relative isolation";
    variables.push_back(var);

    var.set_range(0, 350);

    //Tight lepton pt
    var.name = "TightLeptonsPt";
    var.title = "p_{T} [GeV]";
    var.graph_title = "Transverse momenta of all tight leptons";
    variables.push_back(var);
  
    var.set_range(-4, 4);

    //Tight lepton phi
    var.name = "TightLeptonsPhi";
    var.title = "#phi";
    var.graph_title = "Pseudorapidity of all tight leptons";
    variables.push_back(var);

    var.set_range(-5, 5);

    //Tight lepton eta
    var.name = "TightLeptonsEta";
    var.title = "#eta";
    var.graph_title = "#eta of all tight leptons";
    variables.push_back(var);

    var.set_range(-2, 2);

    //Tight lepton charge
    var.name = "TightLeptonsCharge";
    var.title = "charge";
    var.graph_title = "Charge of all tight leptons";
    variables.push_back(var);

    var.set_range(0, 0.5);

    //Tight lepton mass
    var.name = "TightLeptonsMass";
    var.title = "m_{l} [GeV]";
    var.graph_title = "Mass of all tight leptons";
    variables.push_back(var);

    var.set_range(0, 1);

    //Tight lepton jet relative isolation
    var.name = "TightLeptonsJetRelIso";
    var.title = "Tight lepton-jet relative isolation";
    var.graph_title = "Tight lepton-jet relative isolation";
    variables.push_back(var);

    var.set_range(0, 350);
   
    //Loose leptons pt
    var.name = "LooseLeptonsPt";
    var.title = "p_{T} [GeV]";
    var.graph_title = "Transverse momenta of all loose leptons";
    variables.push_back(var);    

    var.set_range(-4, 4);

    //Loose leptons phi
    var.name = "LooseLeptonsPhi";
    var.title = "#phi";
    var.graph_title = "Pseudorapidity of all loose leptons";
    variables.push_back(var);    

    var.set_range(-5, 5);
 
    //Loose leptons eta
    var.name = "LooseLeptonsEta";
    var.title = "#eta";
    var.graph_title = "#eta of all loose leptons";
    variables.push_back(var);

    var.set_range(-2, 2);

    //Loose leptons charge
    var.name = "LooseLeptonsCharge";
    var.title = "charge";
    var.graph_title = "Charge of all loose leptons";
    variables.push_back(var);

    var.set_range(0, 0.5);

    //Loose leptons mass
    var.name = "LooseLeptonsMass";
    var.title = "m_{l} [GeV]";
    var.graph_title = "Mass of all loose leptons";
    variables.push_back(var);

    var.set_range(0, 5);

    //Opposite sign leptons
    var.name = "OppositeSign";
    var.title = "Opposite sign leptons";
    var.graph_title = "Opposite sign leptons";
    variables.push_back(var);

    //Same sign leptons
    var.name = "SameSign";
    var.title = "Same sign leptons";
    var.graph_title = "Same sign leptons";
    variables.push_back(var);

    var.set_range(0, 350);

    //Lepton pt
    var.name = "LeadingLeptonPt";
    var.title = "p_{T} [GeV]";
    var.graph_title = "Transverse momenta of the leading lepton";
    variables.push_back(var);

    //Subleading lepton pt
    var.name = "SubleadingLeptonPt";
    var.title = "p_{T} [GeV]";
    var.graph_title = "Transverse momenta of the subleading lepton";
    variables.push_back(var);  

    var.set_range(-4, 4);

    //Leading lepton phi
    var.name = "LeadingLeptonPhi";
    var.title = "#phi";
    var.graph_title = "#phi of the leading lepton";
    variables.push_back(var);

    //Subleading lepton phi
    var.name = "SubleadingLeptonPhi";
    var.title = "#phi";
    var.graph_title = "#phi of the subleading lepton";
    variables.push_back(var);

    var.set_range(-5, 5);

    //Leading lepton eta
    var.name = "LeadingLeptonEta";
    var.title = "#eta";
    var.graph_title = "#eta of the leading lepton";
    variables.push_back(var);

    //Subleading lepton eta
    var.name = "SubleadingLeptonEta";
    var.title = "#eta";
    var.graph_title = "#eta of the subleading lepton";
    variables.push_back(var);

    var.set_range(0, 0.5);

    //Leading lepton mass
    var.name = "LeadingLeptonMass";
    var.title = "m_{l} [GeV]";
    var.graph_title = "Mass of the leading lepton";
    variables.push_back(var);

    //Subleading lepton mass
    var.name = "SubleadingLeptonMass";
    var.title = "m_{l} [GeV]";
    var.graph_title = "Mass of the subleading lepton";
    variables.push_back(var);

    var.set_range(0, 1);
    
    //Leading lepton relative isolation
    var.name = "LeadingLepton_RelIso_Selection";
    var.title = "relative isolation";
    var.graph_title = "Relative isolation of the subleading lepton";
    variables.push_back(var);

    var.set_range(-1, 2);

    //Lepton parton flavour (generator-level)
    var.name = "LeptonGenPartFlav";
    var.title = "Parton flavour of generator-level leptons";
    var.graph_title = "parton flavour";   
    variables.push_back(var);

    //Tight lepton parton flavour (generator-level)
    var.name = "TightLeptonsGenPartFlav";
    var.title = "Parton flavour of tight generator-level leptons";
    var.graph_title = "parton flavour";
    variables.push_back(var);

    var.set_range(0, 5);

    //Number of opposite sign non prompt leptons
    var.name = "OppositeSignNonPrompt";
    var.title = "Opposite sign non prompt leptons";
    var.graph_title = "Opposite sign non prompt leptons";
    variables.push_back(var);   

    //Number of opposite sign prompt leptons
    var.name = "OppositeSignPrompt";
    var.title = "Opposite sign prompt leptons";
    var.graph_title = "Opposite sign prompt leptons";
    variables.push_back(var);

    //Number of same sign non prompt leptons
    var.name = "SameSignNonPrompt";
    var.title = "Same sign non prompt leptons";
    var.graph_title = "Same sign non prompt leptons";
    variables.push_back(var);    

    //Number of same sign prompt leptons
    var.name = "SameSignPrompt";
    var.title = "Same sign prompt leptons";
    var.graph_title = "Same sign prompt leptons";
    variables.push_back(var);

    var.set_range(0, 400);

    //Pt of the reco Z
    var.name = "RecoZPt";
    var.title = "p_{T} [GeV]";
    var.graph_title = "Transverse momentum of the reconstructed Z quark candidate";
    variables.push_back(var); 

    var.set_range(-4, 4);
 
    //Phi of the reco Z
    var.name = "RecoZPhi";
    var.title = "#phi";
    var.graph_title = "Pseudorapidity of the reconstructed Z quark candidate";
    variables.push_back(var);	
  
    var.set_range(-5, 5);

    //Eta of the reco Z
    var.name = "RecoZEta";
    var.title = "#eta";
    var.graph_title = "#eta of the reconstructed Z quark candidate";
    variables.push_back(var);

    var.set_range(0, 250);

    //Mass of the reco Z
    var.name = "RecoZMass";
    var.title = "m_{Z} [GeV]";
    var.graph_title = "#eta of the reconstructed Z quark candidate";
    variables.push_back(var);

    var.set_range(0, 6);

    //Delta R between the final state leptons
    var.name = "dR_ll";
    var.title = "#Delta(R)_{ll}";
    var.graph_title = "#Delta(R) between the lepton pair";
    variables.push_back(var);

    var.set_range(0, 6);

    //Delta phi between the final state leptons
    var.name = "dPhi_ll";
    var.title = "#Delta(#phi)_{ll}";
    var.graph_title = "#Delta(#phi) between the lepton pair";
    variables.push_back(var);

    var.set_range(0, 750);

    //Pt of rochester-corrected leptons
    var.name = "LeptonPt_RochCorr";
    var.title = "p_{T} [GeV]";
    var.graph_title = "Transverse momenta of the rochester-corrected leptons";   
    variables.push_back(var);

    var.set_range(-5, 5);

    //Eta of rochester-corrected leptons
    var.name = "LeptonEta_RochCorr";
    var.title = "#eta";
    var.graph_title = "#eta of the rochester-corrected leptons";
    variables.push_back(var);

    var.set_range(-4, 4);

    //Phi of rochester-corrected leptons
    var.name = "LeptonPhi_RochCorr";
    var.title = "#phi";
    var.graph_title = "#phi of the rochester-corrected leptons";
    variables.push_back(var);

    var.set_range(0, 250);

    //Mass of rochester-corrected leptons
    var.name = "LeptonMass_RochCorr";
    var.title = "m_{l}";
    var.graph_title = "Mass of the rochester-corrected leptons";
    variables.push_back(var);
*/
    var.set_range(0, 350);

    //Reco Z mass
    var.name = "z_mass";
    var.title = "m_{ll}";
    var.graph_title = "Reconstructed mass of the Z boson candidate";
    variables.push_back(var);
/*
    var.set_range(0, 3.5);

    //sJER scale factor (nominal)
    var.name = "sJER_Nominal";
    var.title = "s_{JER}";
    var.graph_title = "Data-to-simulation core resolution factor (nominal)";
    variables.push_back(var);   

    //sJER scale factor (up)
    var.name = "sJER_up";
    var.title = "s_{JER}";
    var.graph_title = "Data-to-simulation core resolution factor (up variation)";
    variables.push_back(var);

    //sJER scale factor (down)
    var.name = "sJER_down";
    var.title = "s_{JER}";
    var.graph_title = "Data-to-simulation core resolution factor (down variation)";
    variables.push_back(var);

    var.set_range(0, 1);

    //sigma JER (nominal)
    var.name = "sigma_JER";
    var.title = "#sigma_{JER}";
    var.graph_title = "Transverse momenta resolution factor (nominal)";
    variables.push_back(var);

    var.set_range(0, 3);

    //sigma JER (up)
    var.name = "sigma_JER_up";
    var.title = "#sigma_{JER}";
    var.graph_title = "Transverse momenta resolution factor (up variation)";
    variables.push_back(var);


    //sigma JER (down)
    var.name = "sigma_JER_down";
    var.title = "#sigma_{JER}";
    var.graph_title = "Transverse momenta resolution factor (down variation)";
    variables.push_back(var);

    var.set_range(0, 2.5);

    //cJER
    var.name = "cJER";
    var.title = "Jet energy smearing correction factor";
    var.graph_title = "Jet energy smearing correction factor";
    variables.push_back(var);

    var.set_range(0, 400);

    //Smeeared jet pt
    var.name = "SmearedJetPt";
    var.title = "p_{T} [GeV]";
    var.graph_title = "Transverse momenta of smeared jets";
    variables.push_back(var);

    var.set_range(-4, 4);

    //Smeared jet phi
    var.name = "SmearedJetPhi";
    var.title = "#phi";
    var.graph_title = "#phi of smeared jets";
    variables.push_back(var);

    var.set_range(-6, 6);
   
    //Smeared jet eta
    var.name = "SmearedJetEta";
    var.title = "#eta";
    var.graph_title = "Pseudorapidity of smeared jets";
    variables.push_back(var);

    var.set_range(0, 50);

    //Smeared jet mass
    var.name = "SmearedJetMass";
    var.title = "m_{j} [GeV]";
    var.graph_title = "Mass of smeared jets";
    variables.push_back(var);

    var.set_range(0, 6);

    //Delta R between jets and leptons 
    var.name = "dRJet_Lepton";
    var.title = "#Delta(R)_jl";
    var.graph_title = "#Delta(R) between jets and leptons";
    variables.push_back(var);

    var.set_range(0, 450);

    //Pt of tight smeared jets
    var.name = "TightSmearedJetsPt";
    var.title = "p_{T} [GeV]";
    var.graph_title = "Transverse momenta of tight smeared jets";
    variables.push_back(var);

    var.set_range(-4, 4);

    //Phi of tight smeared jets
    var.name = "TightSmearedJetsPhi";
    var.title = "#phi";
    var.graph_title = "#phi of tight smeared jets";
    variables.push_back(var);

    var.set_range(-6, 6);

    //Eta of tight smeared jets
    var.name = "TightSmearedJetsEta";
    var.title = "#eta";
    var.graph_title = "Pseudorapidity of tight smeared jets";
    variables.push_back(var);
   
    var.set_range(0, 100);

    //Mass of tight smeared jets
    var.name = "TightSmearedJetsMass";
    var.title = "m_{j} [GeV]";
    var.graph_title = "Mass of tight smeared jets";
    variables.push_back(var);

    var.set_range(0, 1);

    //CSVv2 discrimminant of tight smeared jets
    var.name = "TightSmearedJetsBTagCSVV2";
    var.title = "CSVv2 discriminant";
    var.graph_title = "CSVv2 discriminant of tight smeared jets";
    variables.push_back(var);   

    var.set_range(-10, 10);

    //Hadron flavour of tight smeared jets
    var.name = "TightSmearedJetsHadronFlavour";
    var.title = "Hadron flavour";
    var.graph_title = "Hadron flavour of tight smeared jets";
    variables.push_back(var);

    //Jet ID of tight smeared jets
    var.name = "TightSmearedJetsJetID";
    var.title = "Jet ID";
    var.graph_title = "Jet ID of tight smeared jets";
    variables.push_back(var);

    var.set_range(0, 1000);

    //Leading jet mass
    var.name = "LeadingJetMass";
    var.title = "m_{j} [GeV]";
    var.graph_title = "Mass of the leading jet";
    variables.push_back(var);

    var.set_range(0, 70);

    //Subleading jet mass
    var.name = "SubleadingJetMass";
    var.title = "m_{j} [GeV]";
    var.graph_title = "Mass of the subleading jet";
    variables.push_back(var);

    //Third jet mass
    var.name = "ThirdJetMass";
    var.title = "m_{j} [GeV]";
    var.graph_title = "Mass of the third jet";
    variables.push_back(var);

    //Fourth jet mass
    var.name = "FourthJetMass";
    var.title = "m_{j} [GeV]";
    var.graph_title = "Mass of the fourth jet";
    variables.push_back(var);

    var.set_range(0, 500);

    //Leading jet pt
    var.name = "LeadingJetPt";
    var.title = "p_{T} [GeV]";
    var.graph_title = "Transverse momentum of the leading jet";
    variables.push_back(var);

    //Subleading jet pt
    var.name = "SubleadingJetPt";
    var.title = "p_{T} [GeV]";
    var.graph_title = "Transverse momentum of the subleading jet";
    variables.push_back(var);

    var.set_range(0, 300);

    //Third jet pt
    var.name = "ThirdJetPt";
    var.title = "p_{T} [GeV]";
    var.graph_title = "Transverse momentum of the third jet";
    variables.push_back(var);


    var.set_range(0, 200);

    //Fourth jet pt
    var.name = "FourthJetPt";
    var.title = "p_{T} [GeV]";
    var.graph_title = "Transverse momentum of the fourth jet";
    variables.push_back(var);


    var.set_range(0, 200000);

    //Sum of pt squared of all jets
    var.name = "SumSquaredPt";
    var.title = "#sum p_{T}^{2} [GeV]";
    var.graph_title = "Summation of the transverse momenta of all jets squared";
    variables.push_back(var);

    var.set_range(0, 2000);

    //Sum of jet pt
    var.name = "JetPtSum";
    var.title = "#sum p_{T} [GeV]";
    var.graph_title = "Summation of the transverse momenta of all jets";
    variables.push_back(var);

    var.set_range(-6, 6);

    //Leading jet eta
    var.name = "LeadingJetEta";
    var.title = "#eta";
    var.graph_title = "Pseudorapidity of the leading jet";
    variables.push_back(var);

    //Subleading jet eta
    var.name = "SubleadingJetEta";
    var.title = "#eta";
    var.graph_title = "Pseudorapidity of the subleading jet";
    variables.push_back(var);

    //Third jet eta
    var.name = "ThirdJetEta";
    var.title = "#eta";
    var.graph_title = "Pseudorapidity of the third jet";
    variables.push_back(var);

    //Fourth jet eta
    var.name = "FourthJetEta";
    var.title = "#eta";
    var.graph_title = "Pseudorapidity of the fourth jet";
    variables.push_back(var);

    var.set_range(-4, 4);

    //Leading jet phi
    var.name = "LeadingJetPhi";
    var.title = "#phi";
    var.graph_title = "#phi of the leading jet";
    variables.push_back(var);

    //Subleading jet phi
    var.name = "SubleadingJetPhi";
    var.title = "#phi";
    var.graph_title = "#phi of the subleading jet";
    variables.push_back(var);

    //Third jet phi
    var.name = "ThirdJetPhi";
    var.title = "#phi";
    var.graph_title = "#phi of the third jet";
    variables.push_back(var);

    //Fourth jet phi
    var.name = "FourthJetPhi";
    var.title = "#phi";
    var.graph_title = "#phi of the fourth jet";
    variables.push_back(var);

    var.set_range(0, 6);

    //Delta R between the leading and subleading jets
    var.name = "dR_j1j2";
    var.title = "#Delta(R)_{j1j2}";
    var.graph_title = "#Delta(R) between the leading and subleading jets";
    variables.push_back(var);

    //Delta phi between leading and subleading jets
    var.name = "dPhi_j1j2";
    var.title = "#Delta(#phi)_{j1j2}";
    var.graph_title = "#Delta(#phi) between the leading and subleading jets";
    variables.push_back(var);

    var.set_range(0, 1000);

    //Leading jet HT
    var.name = "LeadingJetHT";
    var.title = "H_{T}";
    var.graph_title = "H_{T} of the leading jet";
    variables.push_back(var);

    var.set_range(0, 400);

    //Subleading jet HT
    var.name = "SubleadingJetHT";
    var.title = "H_{T}";
    var.graph_title = "H_{T} of the subleading jet";
    variables.push_back(var);

    var.set_range(0, 300);

    //Third jet HT
    var.name = "ThirdJetHT";
    var.title = "H_{T}";
    var.graph_title = "H_{T} of the third jet";
    variables.push_back(var);

    var.set_range(0, 200);

    //Fourth jet HT
    var.name = "FourthJetHT";
    var.title = "H_{T}";
    var.graph_title = "H_{T} of the fourth jet";
    variables.push_back(var);

    var.set_range(0, 2000);

    //Total HT of all jets
    var.name = "TotJetHT";
    var.title = "H_{T}";
    var.graph_title = "Total H_{T} of all jets";
    variables.push_back(var);

    var.set_range(0, 300);

    //Leading lepton HT
    var.name = "LeadingLeptonHT";
    var.title = "H_{T}";
    var.graph_title = "H_{T} of the leading lepton";
    variables.push_back(var);    

    //Subleading lepton HT
    var.name = "SubleadingLeptonHT";
    var.title = "H_{T}";
    var.graph_title = "H_{T} of the subleading lepton";
    variables.push_back(var);

    var.set_range(0, 500);

    //Total HT of all leptons
    var.name = "TotLepHT";
    var.title = "H_{T}";
    var.graph_title = "Total H_{T} of all leptons";
    variables.push_back(var);

    var.set_range(0, 2);

    //Total HT over total pt of all jets
    var.name = "TotHTOverTotpT_Jets";
    var.title = "#sum H_{T}/p_{T}";
    var.graph_title = "H_{T}/p_{T} for all jets";
    variables.push_back(var);

    var.set_range(0, 1000);

    //Lepton pt sum
    var.name = "LepPtSum";
    var.title = "#sum p_{T}_{l} [GeV]";
    var.graph_title = "Sum of the transverse momenta of all leptons";
    variables.push_back(var);

    var.set_range(0, 400);

    //Lepton eta sum
    var.name = "LepEtaSum";
    var.title = "#sum #eta_{l}";
    var.graph_title = "Sum of the pseudorapidities of all leptons";
    variables.push_back(var);

    var.set_range(0, 400);

    //Lepton phi sum
    var.name = "LepPhiSum";
    var.title = "#sum#phi_{l}";
    var.graph_title = "Sum of the azimuthal angles of all leptons";
    variables.push_back(var);

    var.set_range(0, 2);

    //Total HT over total pt of all leptons
    var.name = "TotHTOverTotpT_Leptons";
    var.title = "H_{T}/p_{T}";
    var.graph_title = "H_{T}/p_{T} for all leptons";
    variables.push_back(var);

    var.set_range(0, 2000);

    //Invariant mass of all jets 
    var.name = "InvMassAllJets";
    var.title = "m_{j} [GeV]";
    var.graph_title = "Invariant mass of all jets";
    variables.push_back(var);


    //Invariant mass of the three leading jets
    var.name = "InvMass3Jets";
    var.title = "m_{j} [GeV]";
    var.graph_title = "Invariant mass of the three most leading jets";
    variables.push_back(var);

    var.set_range(-20, 20);

    //Sum of jet etas
    var.name = "JetEtaSum";
    var.title = "#sum#eta_{j}";
    var.graph_title = "Sum of the pseudorapidities of all jets";
    variables.push_back(var);    

    var.set_range(-15, 15);

    //Sum of jet phis
    var.name = "JetPhiSum";
    var.title = "#sum#phi_{j}";
    var.graph_title = "Sum of the azimuthal angles of all jets";
    variables.push_back(var);

    var.set_range(0, 500);

    //pt of bjets
    var.name = "bjets_pt";
    var.title = "p_{T} [GeV]";
    var.graph_title = "Transverse momentum of bjet candidates";
    variables.push_back(var);

    var.set_range(-5, 5);

    //eta of bjets
    var.name = "bjets_eta";
    var.title = "#eta";
    var.graph_title = "Pseudorapidities of bjet candidates";
    variables.push_back(var);

    var.set_range(-4, 4);

    //phi of bjets
    var.name = "bjets_phi";
    var.title = "#phi";
    var.graph_title = "Azimuthal angles of bjet candidates";
    variables.push_back(var);

    var.set_range(0, 100);

    //mass of bjets
    var.name = "bjets_mass";
    var.title = "m_{b} [GeV]";
    var.graph_title = "Masses of bjet candidates";
    variables.push_back(var);

    //var.name = "nbjets";
    //var.title = "Number of bjets";
    //var.graph_title = "Number of bjet candidates";
    //variables.push_back(var);
    
    var.set_range(0, 1);

    var.name = "BTAGEFF_bjet_id_WP";
    var.title = "BTAGEFF_bjet_id_WP";
    var.graph_title = "BTAGEFF_bjet_id_WP";
    variables.push_back(var);

    var.name = "BTAGEFF_nonbjet_id_WP";
    var.title = "BTAGEFF_nonbjet_id_WP";
    var.graph_title = "BTAGEFF_nonbjet_id_WP";
    variables.push_back(var);

    var.name = "BTAGEFF_bjet_id";
    var.title = "BTAGEFF_bjet_id";
    var.graph_title = "BTAGEFF_bjet_id";
    variables.push_back(var);

    var.name = "BTAGEFF_nonbjet_id";
    var.title = "BTAGEFF_nonbjet_id";
    var.graph_title = "BTAGEFF_nonbjet_id";
    variables.push_back(var);

    var.set_range(0, 300);

    var.name = "BTAGEFF_bjet_pt_num";
    var.title = "BTAGEFF_bjet_pt_num";
    var.graph_title = "BTAGEFF_bjet_pt_num";
    variables.push_back(var);
  
    var.set_range(-5, 5); 

    var.name = "BTAGEFF_bjet_eta_num";
    var.title = "BTAGEFF_bjet_eta_num";
    var.graph_title = "BTAGEFF_bjet_eta_num";
    variables.push_back(var);
 
    var.set_range(0.85, 1.05);

    var.name = "BTAGEFF_bjet_Jet_btagCSVV2_num";
    var.title = "BTAGEFF_bjet_Jet_btagCSVV2_num";
    var.graph_title = "BTAGEFF_bjet_Jet_btagCSVV2_num";
    variables.push_back(var);

    var.set_range(-1, 1);

    var.name = "BTAGEFF_bjet_Jet_hadronFlavour_num";
    var.title = "BTAGEFF_bjet_Jet_hadronFlavour_num";
    var.graph_title = "BTAGEFF_bjet_Jet_hadronFlavour_num";
    variables.push_back(var);

    var.set_range(0, 300);

    var.name = "BTAGEFF_nonbjet_pt_num";
    var.title = "BTAGEFF_nonbjet_pt_num";
    var.graph_title = "BTAGEFF_nonbjet_pt_num";
    variables.push_back(var);

    var.set_range(-5, 5);

    var.name = "BTAGEFF_nonbjet_eta_num";
    var.title = "BTAGEFF_nonbjet_eta_num";
    var.graph_title = "BTAGEFF_nonbjet_eta_num";
    variables.push_back(var);

    var.set_range(0, 1);

    var.name = "BTAGEFF_nonbjet_Jet_btagCSVV2_num";
    var.title = "BTAGEFF_nonbjet_Jet_btagCSVV2_num";
    var.graph_title = "BTAGEFF_nonbjet_Jet_btagCSVV2_num";
    variables.push_back(var);
 
    var.set_range(-10, 10);

    var.name = "BTAGEFF_nonbjet_Jet_hadronFlavour_num";
    var.title = "BTAGEFF_nonbjet_Jet_hadronFlavour_num";
    var.graph_title = "BTAGEFF_nonbjet_Jet_hadronFlavour_num";
    variables.push_back(var);

    var.set_range(0, 300);

    var.name = "BTAGEFF_bjet_pt_denom";
    var.title = "BTAGEFF_bjet_pt_denom";
    var.graph_title = "BTAGEFF_bjet_pt_denom";
    variables.push_back(var);

    var.set_range(-5, 5);

    var.name = "BTAGEFF_bjet_eta_denom";
    var.title = "BTAGEFF_bjet_eta_denom";
    var.graph_title = "BTAGEFF_bjet_eta_denom";
    variables.push_back(var);

    var.set_range(0, 300);

    var.name = "BTAGEFF_nonbjet_pt_denom";
    var.title = "BTAGEFF_nonbjet_pt_denom";
    var.graph_title = "BTAGEFF_nonbjet_pt_denom";
    variables.push_back(var);

    var.name = "BTAGEFF_nonbjet_eta_denom";
    var.title = "BTAGEFF_nonbjet_eta_denom";
    var.graph_title = "BTAGEFF_nonbjet_eta_denom";
    variables.push_back(var);

    var.name = "BTAGEFF_nonbjet_eta_denom";
    var.title = "BTAGEFF_nonbjet_eta_denom";
    var.graph_title = "BTAGEFF_nonbjet_eta_denom";
    variables.push_back(var);

    var.set_range(0, 500);

    //pt of the reco W
    var.name = "w_pair_pt";
    var.title = "p_{T} [GeV]";
    var.graph_title = "Transverse momentum of the reconstructed W boson candidate";
    variables.push_back(var);

    var.set_range(-5, 5);

    //eta of the reco W
    var.name = "w_pair_eta";
    var.title = "#eta"; 
    var.graph_title = "Pseudorapidity of the reconstructed W boson candidate";   
    variables.push_back(var);

    var.set_range(-4, 4);

    //phi of the reco W
    var.name = "w_pair_phi";
    var.title = "#phi";       
    var.graph_title = "#phi of the reconstructed W boson candidate";
    variables.push_back(var);

    var.set_range(0, 100);

    //mass of the reco W
    var.name = "w_pair_mass";
    var.title = "m_{qq} [GeV]";
    var.graph_title = "Mass of the reconstructed W boson candidate";
    variables.push_back(var);

    var.set_range(0, 1000);

    //invariant mass of the reco W
    var.name = "w_mass";
    var.title = "m_{qq} [GeV]";
    var.graph_title = "Invariant mass of the reconstructed W boson candidate";
    variables.push_back(var);
    
    var.set_range(0, 500);

    //leading jet pt from the reco W 
    var.name = "WPairJet1Pt";
    var.title = "p_{T} [GeV]";
    var.graph_title = "Transverse momenta of the leading jet in the W boson candidate reconstruction";
    variables.push_back(var);

    var.set_range(0, 200);

    //subleading jet pt from the reco W 
    var.name = "WPairJet2Pt";
    var.title = "p_{T} [GeV]";
    var.graph_title = "Transverse momenta of the subleading jet in the W boson candidate reconstruction";
    variables.push_back(var);

    var.set_range(-5, 5);

    //leading jet eta from the reco W 
    var.name = "WPairJet1Eta";
    var.title = "#eta";
    var.graph_title = "Pseudorapidity of the leading jet in the W boson candidate reconstruction";
    variables.push_back(var);

    //subleading jet eta from the reco W
    var.name = "WPairJet2Eta";
    var.title = "#eta";
    var.graph_title = "Pseudorapidity of the subleading jet in the W boson candidate reconstruction";
    variables.push_back(var);

    var.set_range(-4, 4);

    //leading jet phi from the reco W
    var.name = "WPairJet1Phi";
    var.title = "#phi";
    var.graph_title = "Azimuthal angle of the leading jet in the W boson candidate reconstruction";
    variables.push_back(var);

    //subleading jet phi from the reco W
    var.name = "WPairJet2Phi";
    var.title = "#phi";
    var.graph_title = "Azimuthal angle of the subleading jet in the W boson candidate reconstruction";
    variables.push_back(var);

    var.set_range(0, 100);

    //leading jet mass from the reco W
    var.name = "WPairJet1Mass";
    var.title = "m_{j_{1}} [GeV]";
    var.graph_title = "Mass of the leading jet in the W boson candidate reconstruction";
    variables.push_back(var);

    var.set_range(0, 200);

    //subleading jet mass from the reco W
    var.name = "WPairJet2Pt";
    var.title = "m_{j_{2}} [GeV]";
    var.graph_title = "Mass of the subleading jet in the W boson candidate reconstruction";
    variables.push_back(var);

    var.set_range(0, 6);

    //Delta R between the leading and subleading jets in the W pair
    var.name = "dR_WJet1_WJet2";
    var.title = "#Delta(R)_{j_{1}j_{2}}";
    var.graph_title = "#Delta(R) between the leading and subleading jets in the W boson candidate reconstruction";
    variables.push_back(var);

    var.set_range(0, 6);

    //Delta phi between the leading and subleading jets in the W pair
    var.name = "dWj1j2";
    var.title = "#Delta(#phi)_{j_{1}j_{2}}";
    var.graph_title = "#Delta(#phi) between the leading and subleading jets in the W boson candidate reconstruction";
    variables.push_back(var);

    var.set_range(0, 6);
 
    //Delta R between the leading jet in the W boson candidate reconstruction and the leading lepton
    var.name = "dR_WJet1_LeadingLepton";
    var.title = "#Delta(R)_{j_{1}l_{1}}";
    var.graph_title = "#Delta(R) between the leading jet in the W boson candidate reconstruction and the leading lepton";
    variables.push_back(var);

    //Delta R between the leading jet in the W boson candidate reconstruction and the subleading lepton
    var.name = "dR_WJet1_SubleadingLepton";
    var.title = "#Delta(R)_{j_{1}l_{2}}";
    var.graph_title = "#Delta(R) between the leading jet in the W boson candidate reconstruction and the subleading lepton";
    variables.push_back(var);

    //Delta R between the subleading jet in the W boson candidate reconstruction and the leading lepton
    var.name = "dR_WJet2_LeadingLepton";
    var.title = "#Delta(R)_{j_{2}l_{1}}";
    var.graph_title = "#Delta(R) between the subleading jet in the W boson candidate reconstruction and the leading lepton";
    variables.push_back(var);

    //Delta R between the subleading jet in the W boson candidate reconstruction and the subleading lepton
    var.name = "dR_WJet2_SubleadingLepton";
    var.title = "#Delta(R)_{j_{2}l_{2}}";
    var.graph_title = "#Delta(R) between the subleading jet in the W boson candidate reconstruction and the subleading lepton";
    variables.push_back(var);

    //Delta R between the leading jet in the W boson candidate reconstruction and the leading jet
    var.name = "dR_WJet1_LeadingJet";
    var.title = "#Delta(R)_{W_{j_{1}}j_{1}}";
    var.graph_title = "#Delta(R) between the leading jet in the W boson candidate reconstruction and the leading jet of the system";
    variables.push_back(var);

    //Delta R between the leading jet in the W boson candidate reconstruction and the subleading jet
    var.name = "dR_WJet1_SubleadingJet";
    var.title = "#Delta(R)_{W_{j_{1}}j_{2}}";
    var.graph_title = "#Delta(R) between the leading jet in the W boson candidate reconstruction and the subleading jet of the system";
    variables.push_back(var);

    //Delta R between the leading jet in the W boson candidate reconstruction and the third jet
    var.name = "dR_WJet1_ThirdJet";
    var.title = "#Delta(R)_{W_{j_{1}}j_{3}}";
    var.graph_title = "#Delta(R) between the leading jet in the W boson candidate reconstruction and the third jet of the system";
    variables.push_back(var);

    //Delta R between the leading jet in the W boson candidate reconstruction and the fourth jet 
    var.name = "dR_WJet1_FourthJet";
    var.title = "#Delta(R)_{W_{j_{1}}j_{4}}";
    var.graph_title = "#Delta(R) between the leading jet in the W boson candidate reconstruction and the fourth jet of the system";
    variables.push_back(var);

    //Delta R between the subleading jet in the W boson candidate reconstruction and the leading jet
    var.name = "dR_WJet2_LeadingJet";
    var.title = "#Delta(R)_{W_{j_{2}}j_{1}}";
    var.graph_title = "#Delta(R) between the subleading jet in the W boson candidate reconstruction and the leading jet of the system";
    variables.push_back(var);

    //Delta R between the subleading jet in the W boson candidate reconstruction and the subleading jet
    var.name = "dR_WJet2_SubleadingJet";
    var.title = "#Delta(R)_{W_{j_{2}}j_{2}}";
    var.graph_title = "#Delta(R) between the subleading jet in the W boson candidate reconstruction and the subleading jet of the system";
    variables.push_back(var);

    //Delta R between the subleading jet in the W boson candidate reconstruction and the third jet
    var.name = "dR_WJet2_ThirdJet";
    var.title = "#Delta(R)_{W_{j_{2}}j_{3}}";
    var.graph_title = "#Delta(R) between the subleading jet in the W boson candidate reconstruction and the third jet of the system";
    variables.push_back(var);
 
    //Delta R between the subleading jet in the W boson candidate reconstruction and the fourth jet
    var.name = "dR_WJet2_FourthJet";
    var.title = "#Delta(R)_{W_{j_{2}}j_{4}}";
    var.graph_title = "#Delta(R) between the subleading jet in the W boson candidate reconstruction and the fourth jet of the system";
    variables.push_back(var);

    var.set_range(0, 500);

    //HT of the leading jet in the W pair
    var.name = "WJet1HT";
    var.title = "H_{T} [GeV]";
    var.graph_title = "H_{T} of the leading jet in the W boson candidate reconstruction";
    variables.push_back(var);

    var.set_range(0, 200);

    //HT of the subleading jet in the W pair
    var.name = "WJet2HT";
    var.title = "H_{T} [GeV]";
    var.graph_title = "H_{T} of the subleading jet in the W boson candidate reconstruction";
    variables.push_back(var);

    var.set_range(0, 500);

    //HT of the reco W
    var.name = "RecoWHT";
    var.title = "H_{T} [GeV]";
    var.graph_title = "H_{T} of the reconstructed W boson candidate";
    variables.push_back(var); 

    var.set_range(0, 400);

    //Transverse W mass
    var.name = "mtW";
    var.title = "m_{T}_{jj} [GeV]";
    var.graph_title = "Transverse mass of the reconstructed W boson candidate";
    variables.push_back(var);

    var.set_range(0, 10);

    //Number of tight smeared jets
    var.name = "TightSmearedJetsNumber";
    var.title = "number";
    var.graph_title = "Number of tight smeared jets";
    variables.push_back(var);

    var.set_range(0, 100);

    //Mass of the leading smeared bjet
    var.name = "SmearedLeadingBJetMass";
    var.title = "m_{b} [GeV]";
    var.graph_title = "Mass of the smeared leading b jet";
    variables.push_back(var);

    var.set_range(0, 1000);

    //Pt of the leading smeared bjet
    var.name = "SmearedLeadingBJetPt";
    var.title = "p_{T} [GeV]";
    var.graph_title = "p_{T} of the smeared leading b jet";
    variables.push_back(var);

    var.set_range(-5, 5);

    //Eta of the leading smeared bjet
    var.name = "SmearedLeadingBJetEta";
    var.title = "#eta";
    var.graph_title = "Pseudorapidity of the smeared leading b jet";
    variables.push_back(var);

    var.set_range(-4, 4);

    //Phi of the leading smeared bjet
    var.name = "SmearedLeadingBJetPhi";
    var.title = "#phi";
    var.graph_title = "Azimuthal angle of the smeared leading b jet";
    variables.push_back(var);  

    var.set_range(0, 1000);

    //Top pt
    var.name = "Top_Pt";
    var.title = "p_{T}";
    var.graph_title = "Transverse momentum of the reconstructed top quark candidate";
    variables.push_back(var);

    var.set_range(-5, 5);

    //Top eta
    var.name = "Top_Eta";
    var.title = "#eta";
    var.graph_title = "Pseudorapidity of the reconstructed top quark candidate";
    variables.push_back(var);

    var.set_range(-4, 4);

    //Top phi
    var.name = "Top_Phi";
    var.title = "#phi";
    var.graph_title = "Azimuthal angle of the reconstructed top quark candidate";
    variables.push_back(var);

    var.set_range(0, 1000);

    //Top mass
    var.name = "Top_Mass";
    var.title = "m_{t} [GeV]";
    var.graph_title = "Mass of the reconstructed top quark candidate";
    variables.push_back(var);

    //Top HT
    var.name = "Top_HT";
    var.title = "H_{T}";
    var.graph_title = "H_{T} of the reconstructed top quark candidate";
    variables.push_back(var);

    var.set_range(0, 6);

    //Delta R between the reco top and leading lepton
    var.name = "dR_Top_LeadingElectron";
    var.title = "#Delta(R)_{l_{1}t}";
    var.graph_title = "#Delta(R) between the reconstructed top quark candidate and the leading lepton";
    variables.push_back(var);


    //Delta R between the reco top and subleading lepton
    var.name = "dR_Top_SubleadingElectron";
    var.title = "#Delta(R)_{l_{2}t}";
    var.graph_title = "#Delta(R) between the reconstructed top quark candidate and the subleading lepton";
    variables.push_back(var);

    //Delta R between the reco top and leading jet
    var.name = "dR_Top_LeadingJet";
    var.title = "#Delta(R)_{j_{1}t}";
    var.graph_title = "#Delta(R) between the reconstructed top quark candidate and the leading jet";
    variables.push_back(var);

    //Delta R between the reco top and subleading jet
    var.name = "dR_Top_SubleadingJet";
    var.title = "#Delta(R)_{j_{2}t}";
    var.graph_title = "#Delta(R) between the reconstructed top quark candidate and the subleading jet";
    variables.push_back(var);

    //Delta R between the reco top and third jet
    var.name = "dR_Top_ThirdJet";
    var.title = "#Delta(R)_{j_{3}t}";
    var.graph_title = "#Delta(R) between the reconstructed top quark candidate and the third jet";
    variables.push_back(var);

    //Delta R between the reco top and fourth jet
    var.name = "dR_Top_FourthJet";
    var.title = "#Delta(R)_{j_{4}t}";
    var.graph_title = "#Delta(R) between the reconstructed top quark candidate and the fourth jet";
    variables.push_back(var);

    //Delta R between the reco top and reco W
    var.name = "dR_Top_W";
    var.title = "#Delta(R)_{tW}";
    var.graph_title = "#Delta(R) between the reconstructed top quark candidate and the reconstructed W boson candidate";
    variables.push_back(var);
 
    //Delta phi between the reconstructed top quark candidate and the leading jet from the W reconstruction
    var.name = "dPhi_Wj1_Top";
    var.title = "#Delta(#phi)_{W_{j1}t}";
    var.graph_title = "#Delta(#phi) between the reconstructed top quark candidate and the leading jet from the W reconstruction";
    variables.push_back(var);

    //Delta phi between the reconstructed top quark candidate and the subleading jet from the W reconstruction
    var.name = "dPhi_Wj2_Top";
    var.title = "#Delta(#phi)_{W_{j2}t}";
    var.graph_title = "#Delta(#phi) between the reconstructed top quark candidate and the subleading jet from the W reconstruction";
    variables.push_back(var);

    //Delta R between the reconstructed Z boson and top quark candidates
    var.name = "dR_Z_Top";
    var.title = "#Delta(R)_{Zt}";
    var.graph_title = "#Delta(R) between the reconstructed Z boson and top quark candidates";
    variables.push_back(var);

    //Delta phi between the reconstructed Z boson and top quark candidate
    var.name = "dPhi_Z_Top";
    var.title = "#Delta(#phi)_{Zt}";
    var.graph_title = "#Delta(#phi) between the reconstructed Z boson and top quark candidates";
    variables.push_back(var);

    //Delta R between the reconstructed Z boson and the leading jet in the W boson candidate reconstruction
    var.name = "dR_Z_WPairJet1";
    var.title = "#Delta(R)_{W_{j1}Z}";
    var.graph_title = "#Delta(R) between the reconstructed Z boson and the leading jet in the W boson candidate reconstruction";
    variables.push_back(var);

    //Delta R between between the reconstructed Z boson and the subleading jet in the W boson candidate reconstruction
    var.name = "dR_Z_WPairJet2";
    var.title = "#Delta(R)_{W_{j2}Z}";
    var.graph_title = "#Delta(R) between the reconstructed Z boson and the subleading jet in the W boson candidate reconstruction";
    variables.push_back(var);

    //Delta phi between the reconstructed Z boson and the leading jet in the W boson candidate reconstruction
    var.name = "dPhi_Z_WPairJet1";
    var.title = "#Delta(#phi)_{W_{j1}Z}";
    var.graph_title = "#Delta(#phi) between the reconstructed Z boson and the leading jet in the W boson candidate reconstruction";
    variables.push_back(var);

    //Delta phi between the reconstructed Z boson and the subleading jet in the W boson candidate reconstruction
    var.name = "dPhi_Z_WPairJet2";
    var.title = "#Delta(#phi)_{W_{j2}Z}";
    var.graph_title = "#Delta(#phi) between the reconstructed Z boson and the subleading jet in the W boson candidate reconstruction";
    variables.push_back(var);

    //Min delta R
    var.name = "MinDeltaR";
    var.title = "#Delta(R)_{min}";
    var.graph_title = "Minimum #Delta(R) of the system";
    variables.push_back(var);

    //Min delta phi
    var.name = "MinDeltaPhi";
    var.title = "#Delta(#phi)_{min}";
    var.graph_title = "Minimum #Delta(#phi) of the system";
    variables.push_back(var); 

    //Delta R between the leading lepton and leading bjet
    var.name = "dR_LeadingLepton_LeadingBJet";
    var.title = "#Delta(R)_{l_{1}b_{1}}";
    var.graph_title = "#Delta(R) between the leading lepton and leading bjet";
    variables.push_back(var);
 
    //Delta R between the subleading lepton and leading bjet
    var.name = "dR_SubleadingLepton_LeadingBJet";
    var.title = "#Delta(R)_{l_{2}b_{1}}";
    var.graph_title = "#Delta(R) between the subleading lepton and leading bjet";
    variables.push_back(var);

    //Delta phi between the leading lepton and leading bjet 
    var.name = "DeltaPhi_Leadinglepton_BJet";
    var.title = "#Delta(#phi)_{l_{1}b_{1}}";
    var.graph_title = "#Delta(#phi) between the leading lepton and leading bjet";
    variables.push_back(var);

    //Delta phi between the subleading lepton and leading bjet
    var.name = "DeltaPhi_Subleadinglepton_BJet";
    var.title = "#Delta(#phi)_{l_{2}b_{1}}";
    var.graph_title = "#Delta(#phi) between the subleading lepton and leading bjet";
    variables.push_back(var);

    var.set_range(0, 3000);

    //MET
    var.name = "MET";
    var.title = "Uncorrected MET";
    var.graph_title = "Missing E_{T}";
    variables.push_back(var);

    var.set_range(0, 1);

    //CSVv2 discrimminant of the leading bjet
    var.name = "LeadingBJetOutputDiscriminant";
    var.title = "CSVv2 discriminant";
    var.graph_title = "CSVv2 discriminant of the leading bjet";
    variables.push_back(var);

    //CSVv2 discrimminant of the subleading bjet
    var.name = "SubleadingBJetOutputDiscriminant";
    var.title = "CSVv2 discriminant";
    var.graph_title = "CSVv2 discriminant of the subleading bjet";
    variables.push_back(var);

    //CSVv2 discrimminant of the third bjet
    var.name = "ThirdBJetOutputDiscriminant";
    var.title = "CSVv2 discriminant";
    var.graph_title = "CSVv2 discriminant of the third bjet";
    variables.push_back(var);

    //CSVv2 discrimminant of the fourth bjet
    var.name = "FourthBJetOutputDiscriminant";
    var.title = "CSVv2 discriminant";
    var.graph_title = "CSVv2 discriminant of the fourth bjet";
    variables.push_back(var);
   
    var.set_range(0, 6);

    //Delta phi between the reconstructed top quark and W boson candidates
    var.name = "dPhi_W_Top";
    var.title = "#Delta(#phi)_{tW}";
    var.graph_title = "#Delta(#phi) between the reconstructed top quark and W boson candidates";
    variables.push_back(var);   

    //Delta R between the reconstructed Z boson candidate and the leading jet
    var.name = "dR_Z_LeadingJet";
    var.title = "#Delta(R)_{j_{1}Z}";
    var.graph_title = "#Delta(R) between the reconstructed Z boson candidate and the leading jet";
    variables.push_back(var);

    //Delta R between the reconstructed Z boson candidate and the subleading jet
    var.name = "dR_Z_SubleadingJet";
    var.title = "#Delta(R)_{j_{2}Z}";
    var.graph_title = "#Delta(R) between the reconstructed Z boson candidate and the subleading jet";
    variables.push_back(var);

    //Delta R between the reconstructed Z boson candidate and the third jet
    var.name = "dR_Z_ThirdJet";
    var.title = "#Delta(R)_{j_{3}Z}";
    var.graph_title = "#Delta(R) between the reconstructed Z boson candidate and the third jet";
    variables.push_back(var);

    //Delta R between the reconstructed Z boson candidate and the fourth jet
    var.name = "dR_Z_FourthJet";
    var.title = "#Delta(R)_{j_{4}Z}";
    var.graph_title = "#Delta(R) between the reconstructed Z boson candidate and the fourth jet";
    variables.push_back(var);

    //Delta phi between the reconstructed Z boson candidate and the leading bjet
    var.name = "dPhi_LeadingJet_Z";
    var.title = "#Delta(#phi)_{b_{1}Z}";
    var.graph_title = "#Delta{#phi between the reconstructed Z boson candidate and the leading bjet";
    variables.push_back(var);
  
    //Delta phi between the reconstructed Z boson candidate and the subleading bjet
    var.name = "dPhi_SubleadingJet_Z";
    var.title = "#Delta(#phi)_{b_{2}Z}";
    var.graph_title = "#Delta(#phi) between the reconstructed Z boson candidate and the subleading bjet";
    variables.push_back(var);

    //Delta phi between the reconstructed Z boson candidate and the third bjet
    var.name = "dPhi_ThirdJet_Z";
    var.title = "#Delta(#phi)_{b_{3}Z}";
    var.graph_title = "#Delta(#phi) between the reconstructed Z boson candidate and the third bjet";
    variables.push_back(var);

    var.set_range(0, 6);

    //Delta phi between the reconstructed Z boson candidate and the fourth bjet
    var.name = "dPhi_FourthJet_Z";
    var.title = "#Delta(#phi)_{b_{4}Z}";
    var.graph_title = "#Delta(#phi) between the reconstructed Z boson candidate and the fourth bjet";
    variables.push_back(var);

    //Delta R between the reconstructed W and Z boson candidates
    var.name = "dR_W_Z";
    var.title = "#Delta(R)_{WZ}";
    var.graph_title = "#Delta(R) between the reconstructed W and Z boson candidates";
    variables.push_back(var);

    var.set_range(0, 500);

    //HT of the reco Z
    var.name = "RecoZHT";
    var.title = "H_{T} [GeV]";
    var.graph_title = "H_{T} of the reconstructed Z boson candidate";
    variables.push_back(var);

    var.set_range(0, 6);

    //Delta phi between the reconstructed W and Z boson candidates
    var.name = "dPhi_W_Z";
    var.title = "#Delta(#phi)_{WZ}";
    var.graph_title = "#Delta(#phi) between the reconstructed W and Z boson candidates";
    variables.push_back(var);

    var.set_range(0, 500);
   
    //Total eta of the system
    var.name = "TotalEta_System";
    var.title = "#sum #eta";
    var.graph_title = "Total #eta of all objects in the system";
    variables.push_back(var);

    //Total phi of the system
    var.name = "TotalPhi_System";
    var.title = "#sum #phi";
    var.graph_title = "Total #phi of the system";
    variables.push_back(var);

    var.set_range(0, 1000);

    //Invariant top mass
    var.name = "InvTopMass";
    var.title = "m_{t} [GeV]";
    var.graph_title = "Invariant mass of the reconstructed top quark candidate";
    variables.push_back(var); 

    var.set_range(0, 5000);

    //Total HT of the system
    var.name = "TotalHT_System";
    var.title = "H_{T} [GeV]";
    var.graph_title = "Total H_{T} of the system";
    variables.push_back(var);
  

    //Total pT of the system
    var.name = "TotalPt_System";
    var.title = "p_{T} [GeV]";
    var.graph_title = "Total p_{T} of the system";
    variables.push_back(var);

    var.set_range(0, 2);

    //Total HT over pT over the system
    var.name = "TotHTOverTotpT_System";
    var.title = "H_{T}/p_{T}";
    var.graph_title = "Total H_{T}/p_{T} of the system";
    variables.push_back(var);

    var.set_range(0, 1000);

    //pt of non b-tagged jets
    var.name = "notbjetpt";
    var.title = "p_{T} [GeV]";
    var.graph_title = "p_{T} of non b-tagged jets";
    variables.push_back(var);

    var.set_range(-5, 5);

    //eta of non b-tagged jets
    var.name = "notbjeteta";
    var.title = "#eta";
    var.graph_title = "Pseudorapidity of non b-tagged jets";
    variables.push_back(var);

    var.set_range(0.8, 1.2);

    //CMS btagging efficiency SF
    var.name = "CMSBTagSF";
    var.title = "CMSBTagSF";
    var.graph_title = "CMS-provided btagging efficiency scale factor (b-tagged jets)";
    variables.push_back(var);

    //CMS btagging efficiency SF for non btagged jets
    var.name = "CMSNonBTagSF";
    var.title = "CMSNonBTagSF";
    var.graph_title = "CMS-provided btagging efficiency scale factor (non b-tagged jets)";
    variables.push_back(var);

    //B-tagged jets efficiency
    var.name = "EffBTagged";
    var.title = "#epsilon_{btagged}";
    var.graph_title = "Efficiency of identifying b-tagged jets";  
    variables.push_back(var);

    //Non-tagged jets efficiency
    var.name = "EffNonBTagged";
    var.title = "#epsilon_{nonbtagged}";
    var.graph_title = "Efficiency of identifying non b-tagged jets";
    variables.push_back(var);

    var.set_range(0, 1.2);

    //Products of the efficiencies for btagged jets
    var.name = "ProductOperator_E_i";
    var.title = "#prod #epsilon_{i}";
    var.graph_title = "#prod #epsilon_{i}";
    variables.push_back(var);

    var.set_range(0.99, 1.01);

    //Products of  1-efficiencies for non-btagged jets
    var.name = "ProductOperator_1_Minus_E_j";
    var.title = "#prod 1-#epsilon_{j}";
    var.graph_title = "#prod 1-#epsilon_{j}";
    variables.push_back(var);

    var.set_range(0, 1.0);

    //Product of the CMS btagging SF times btagging efficiency
    var.name = "ProductOperator_SFi_Times_Ei";
    var.title = "#prod SF_{i}#epsilon_{i}";
    var.graph_title = "#prod Sf_{i}#epsilon_{i}";
    variables.push_back(var);

    //Product of the 1 - non btag SF times non btag efficiency
    var.name = "ProductOperator_1_Minus_SFj_Times_Ej";
    var.title = "#prod 1-SF_{j}#epsilon_{j}";
    var.graph_title = "#prod 1-SF_{j}#epsilon_{j}";
    variables.push_back(var);
  
    var.set_range(0, 1.2);
 
    //MC btagging efficiency
    var.name = "ProbBTagMC";
    var.title = "P(MC)";
    var.graph_title = "P(MC) – btagging efficiency";
    variables.push_back(var);
   
    //Data btagging efficiency
    var.name = "ProbBTagData";
    var.title = "P(Data)";
    var.graph_title = "P(Data) – btagging efficiency";
    variables.push_back(var);

    var.set_range(0, 3);
    
    //Pile up weight
    var.name = "PU";
    var.title = "pile up weight";
    var.graph_title = "Pile up weights";
    variables.push_back(var);


    //B tagging efficiency weight
    var.name = "BTagWeight";
    var.title = "#omega";
    var.graph_title = "B-tagging efficiency weight";
    variables.push_back(var);

    var.set_range(0.8, 1.2);

    //egamma SF
    var.name = "EGammaSF_egammaEff";
    var.title = "E#gamma SF";
    var.graph_title = "E#gamma SF";
    variables.push_back(var);
    
    var.set_range(0, 0.2);

    //egamma systematics
    var.name = "EGammaSF_egammaEffSys";
    var.title = "E#gamma SF sys";
    var.graph_title = "E#gamma SF sys";
    variables.push_back(var);

    var.set_range(0.8, 1.2);

    //egamma reco SF
    var.name = "EGammaSF_egammaEffReco";
    var.title = "E#gamma SF reco";
    var.graph_title = "E#gamma SF reco";
    variables.push_back(var);

    var.set_range(0, 0.2);

    //egamma reco systematics
    var.name = "EGammaSF_egammaEffRecoSys";
    var.title = "E#gamma SF reco sys";
    var.graph_title = "E#gamma SF reco sys";
    variables.push_back(var);

    var.set_range(0, 1.1);

    //Muon ID Sf
    var.name = "MuonSFTest_ID";
    var.title = "Muon ID SF";
    var.graph_title = "Muon ID SF";
    variables.push_back(var);

    //Muon ID SF stat
    var.name = "MuonSFTest_ID_sys_stat";
    var.title = "Muon ID SF (stat.)";
    var.graph_title = "Muon ID SF (stat.)";
    variables.push_back(var);
 
    //Muon ID SF syst
    var.name = "MuonSFTest_ID_sys_syst";
    var.title = "Muon ID SF (syst.)";
    var.graph_title = "Muon ID SF (syst.)";
    variables.push_back(var);

    //Muon Iso SF
    var.name = "MuonSFTest_Iso";
    var.title = "Muon Iso SF";
    var.graph_title = "Muon Iso SF";
    variables.push_back(var);

    //Muon Iso SF stat
    var.name = "MuonSFTest_Iso_sys_stat";
    var.title = "Muon Iso SF (stat.)";
    var.graph_title = "Muon Iso SF (stat.)";
    variables.push_back(var);

    //Muon Iso SF syst
    var.name = "MuonSFTest_Iso_sys_syst";
    var.title = "Muon Iso SF (syst.)";
    var.graph_title = "Muon Iso SF (syst.)";
    variables.push_back(var);

    var.set_range(0.99, 1.01);

    //parton shower weights
    var.name = "ReturnedPSWeight";
    var.title = "PS weight";
    var.graph_title = "Parton shower weight";
    variables.push_back(var);

    var.set_range(0.8, 1.6);

    //PDF weights
    var.name = "CalculatedPDFWeight";
    var.title = "Parton distribution function weight";
    var.graph_title = "Parton distribution function weight";
    variables.push_back(var);

    var.set_range(0.99, 1.01);

    //Generator weight
    var.name = "CalculatedGeneratorWeight";
    var.title = "Generator weight";
    var.graph_title = "Generator weight";
    variables.push_back(var);

    var.set_range(-1, 5);

    //Event weight
    var.name = "EventWeight";
    var.title = "Event weight";
    var.graph_title = "Event weight";
    variables.push_back(var);

    var.set_range(0, 1000);   
 
    //Missing transverse energy (uncorrected) 
    //var.name = "OriginalMET";
    //var.title = "Original MET";
    //var.graph_title = "Original MET";
    //variables.push_back(var);

    //var.name = "ScaledMET";
    //var.title = "Scaled MET";
    //var.graph_title = "Scaled MET";
    //variables.push_back(var);
  
    //var.name = "newMET";
    //var.title = "Corrected MET";
    //var.graph_title = "Corrected MET";
    //variables.push_back(var);  

    //var.name = "UnsmearedJet4Momentum";
    //var.title = "Four momenta of unsmeared jets";
    //var.graph_title = "Four momenta";
    //variables.push_back(var);

    var.set_range(0, 30);

    //Overall normalisation
    var.name = "OverallNormalisation";
    var.title = "Overall normalisation";
    var.graph_title = "Overall normalisation";
    variables.push_back(var);

    //chi2 variable for experimental blinding (before cut)
    //var.name = "chi2";
    //var.title = "#chi^{2}";
    //var.graph_title = "#chi^{2} (before filter)";
    //variables.push_back(var);

    //chi2 variable for experimental blinding (after cut)
    //var.name = "AfterChi2Cut";
    //var.title = "#chi^{2}";
    //var.graph_title = "#chi^{2} (after filter)";
    //variables.push_back(var);
*/
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
                        plotter(var, year, channel, systematic, region, true, false, true);
}
