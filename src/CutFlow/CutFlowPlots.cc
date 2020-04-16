#include <iostream>
#include <fstream>
#include <limits>
#include <stdlib.h> //for atof
#include <vector> 

using namespace std;


std::fstream& GotoLine(std::fstream& file, unsigned int num){
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}


double linereader(const int& LineNumber, const string& process, const string& channel){
    
   string InputFile = "CutFlowReport_" + process + "_" + channel + "_Blinded.txt";
   fstream file(InputFile.c_str());
   GotoLine(file, LineNumber);

   string line;

   double Value = atof(line.c_str());
   return Value;

}

int linecounter(const string& process, const string& channel){ 
    
   int number_of_lines = 0;
   string line;
   string InputFile = "CutFlowReport_" + process + "_" + channel + "_Blinded.txt";
   ifstream myfile(InputFile.c_str());

   while (getline(myfile, line))
        ++number_of_lines;
    	return number_of_lines;

}


auto textfilereader2(const string& process, const string& channel){

   int NumberOfLines = linecounter(process, channel);
   vector<double> Value;

   for(int i = 1; i < NumberOfLines+1; i++){
   	Value.push_back(linereader(i, process, channel));
   }

   return Value;

}


auto textfilereader3{[](const string& process, const string& channel){

  int max = (textfilereader2(process, channel)).size();
  auto vector1 = textfilereader2(process, channel);
  vector<float> OutputVector;  

  cout << vector1.at(1) << endl;
/*
  cout << vector1.at(2) << endl; //1 (ee)
  cout << vector1.at(5) << endl;
  cout << vector1.at(6) << endl; //2 (mumu)
  cout << vector1.at(9) << endl;
  cout << vector1.at(10) << endl; //3 (ee)
  cout << vector1.at(13) << endl;
  cout << vector1.at(14) << endl; //4 (mumu)
  cout << vector1.at(17) << endl;
  cout << vector1.at(18) << endl; //5 (ee)
  cout << vector1.at(21) << endl;
  cout << vector1.at(22) << endl; //6 (mumu)
  cout << vector1.at(25) << endl;
  cout << vector1.at(26) << endl; //7 (ee)
  cout << vector1.at(29) << endl;
  cout << vector1.at(30) << endl; //8 (mumu)
  cout << vector1.at(33) << endl;
  cout << vector1.at(34) << endl; //9 (ee)
  cout << vector1.at(37) << endl;
  cout << vector1.at(38) << endl; //10 (mumu)
  cout << vector1.at(41) << endl;
  cout << vector1.at(42) << endl; //11 (ee)
*/
  OutputVector.push_back( vector1.at(1) );
/*  OutputVector.push_back( vector1.at(2) );
  OutputVector.push_back( vector1.at(5) );
  OutputVector.push_back( vector1.at(6) );
  OutputVector.push_back( vector1.at(9) );
  OutputVector.push_back( vector1.at(10) );
*/

  return OutputVector;

}};


auto CutFlowPlots(){

	return textfilereader3("SingleTop_schannel", "ee");

} 

