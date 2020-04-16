#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include <vector>

void JEC(){

// Instantiate uncertainty sources
const int nsrc = 33;
const char* srcnames[nsrc] =
  {"Absolute", "HighPtExtra",  "SinglePionECAL", "SinglePionHCAL",
   "FlavorQCD", "Time",
   "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF",
   "RelativePtBB","RelativePtEC1", "RelativePtEC2", "RelativePtHF", "RelativeFSR",
   "RelativeStatEC2", "RelativeStatHF",
   "PileUpDataMC", 
   "PileUpPtBB", "PileUpPtEC", "PileUpPtHF","PileUpBias",
   "SubTotalPileUp","SubTotalRelative","SubTotalPt","SubTotalMC",
   "Total","TotalNoFlavor",
   "FlavorZJet","FlavorPhotonJet","FlavorPureGluon","FlavorPureQuark","FlavorPureCharm","FlavorPureBottom"};
std::vector<JetCorrectionUncertainty*> vsrc(nsrc);

for (int isrc = 0; isrc < nsrc; isrc++) {

   const char *name = srcnames[isrc];
   JetCorrectorParameters *p = new JetCorrectorParameters("2017.txt", name);
   JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
   vsrc[isrc] = unc;
} // for isrc

// Total uncertainty for reference
JetCorrectionUncertainty *total = new JetCorrectionUncertainty(*(new JetCorrectorParameters("2017.txt", "Total")));

// Calculate uncertainty per source and as a total
double jetpt(50.);
double jeteta(2.4);
double sum2_up(0), sum2_dw(0);

for (int isrc = 0; isrc < nsrc; isrc++) {

      JetCorrectionUncertainty *unc = vsrc[isrc];
      unc->setJetPt(jetpt);
      unc->setJetEta(jeteta);
      double sup = unc->getUncertainty(true); // up variation
      unc->setJetPt(jetpt);
      unc->setJetEta(jeteta);
      double sdw = unc->getUncertainty(false); // down variation

      sum2_up += pow(max(sup,sdw),2);
      sum2_dw += pow(min(sup,sdw),2);
} // for isrc

total->setJetPt(jetpt);
total->setJetEta(jeteta);
double uncert = total->getUncertainty(true);

// Check that quadratic sum of sources equals total uncertainty
assert(fabs(uncert - sqrt(sum2_up)) < 1e-3);


}
