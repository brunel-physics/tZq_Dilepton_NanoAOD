#include <fstream>
#include <limits>
#include <stdlib.h> //for atof
#include <vector> 
#include <iostream>

std::fstream& GotoLine(std::fstream& file, unsigned int num){
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}


double linereader(const int& LineNumber){
    
   using namespace std;

   fstream file("NormalisationFactorsTest_JustValues.txt");
   GotoLine(file, LineNumber);

   string line;
   file >> line;
   cin.get();

   double Value = atof(line.c_str());
   return Value;

}

int linecounter(){ 
    
   int number_of_lines = 0;
   string line;
   ifstream myfile("NormalisationFactorsTest_JustValues.txt");

   while (getline(myfile, line))
        ++number_of_lines;
    	return number_of_lines;

}


auto textfilereader2(){
//auto textfilereader2(){

   int NumberOfLines = linecounter();
   vector<double> Value;

   for(int i = 1; i < /*NumberOfLines+1*/2; i++){
   	Value.push_back(linereader(i));
   }

   return Value;

}


auto textfilereader{[](){

  int max = (textfilereader2()).size();
  auto vector1 = textfilereader2();
  vector<double> vector2 (max, 1);
  vector<double> OutputVector;  

  cout << vector1.at(0) << endl;
  cout << vector2.at(0) << endl;

  for(int i = 0; i < max; i++){
	OutputVector.push_back( vector1.at(i) + vector2.at(i) );
  }

  return OutputVector;

}};

