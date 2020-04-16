#include <fstream>
#include <iostream>
#include <sstream>
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TH1D.h"
#include <algorithm>
#include <TFile.h>
#include "TLorentzVector.h"
#include <iostream>
#include <fstream>

using namespace ROOT; // RDataFrame's namespace
using namespace std;

using floats = ROOT::VecOps::RVec<float>;
using ints = ROOT::VecOps::RVec<int>;
using bools = ROOT::VecOps::RVec<bool>;
using chars = ROOT::VecOps::RVec<UChar_t>;
using doubles = ROOT::VecOps::RVec<double>;


float JetSmearingTest(){

  vector<string>input_files = {"/data/disk0/nanoAOD_2017/tZq_ll/1283E24E-F279-6F42-B94E-9F73F91F47A3.root"};
  RDataFrame d_dataframe("Events", input_files);

  auto sigmaJER{[](){


  }};




}

