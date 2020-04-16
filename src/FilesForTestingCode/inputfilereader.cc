#include "ROOT/RDataFrame.hxx"

using namespace ROOT; // RDataFrame's namespace
using namespace std;

void inputfilereader2(const string& process, const string& variable){

//string input_file = process + "_AfterFullSelection.root";

string input_file = process + "_" + variable + "_" + "ee_AfterFullSelection.root";

RDataFrame dataframe(variable.c_str(), input_file.c_str());

auto histo = dataframe.Histo1D(process.c_str());

}

void inputfilereader(){

inputfilereader2("tZq", "Electron_pt_Selection");

}
