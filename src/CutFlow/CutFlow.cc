#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;


float Col1, Col2, Col3, Col4;
string FileName;


vector<float> RowReader2(const int& LineSpecified, const string& process, const string& channel) { 


  string FileName = "CutFlowReport_" + process + "_" + channel + "_Blinded.txt";

  ifstream file;
  file.open(FileName.c_str());

  if (file.good())
  {
    string str = "";

    int line_number = 0;
	
	 while(getline(file, str) && line_number != LineSpecified){
		++line_number;
	}
	if(line_number == LineSpecified){
			file >> Col1;
                	file >> Col2;
                	file >> Col3;
                	file >> Col4;
	}

  }
 
  file.close(); 
	
  
  vector<float> OutputVector{};
  OutputVector.push_back(Col3);
  OutputVector.push_back(Col4);

  return OutputVector;

} 





int linecounter(const string& process, const string& channel){ 
    
   int number_of_lines = 0;
   string line;

   string FileName = "CutFlowReport_" + process + "_" + channel + "_Blinded.txt";

   ifstream myfile(FileName);

   while (getline(myfile, line))
        ++number_of_lines;
    	return number_of_lines;

}



auto RowReader3(const string& process, const string& channel){

//  vector<vector<float>> Out{};

  for(int i = 0; i < linecounter(process, channel) + 1; i++){

//	Out.push_back(RowReader2(i, process, channel));
	return RowReader2(i, process, channel);
  }


}



auto CutFlow(){

  return RowReader3("SingleTop_schannel", "ee");

}

