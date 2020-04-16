#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

float Col1, Col4, Col5, Col6, Col7, Col8, Col9, Col10;
string FileName, Col2, Col3, Col11;


string BTagFileReader2(const int& LineSpecified, const string& year){

  if(year == "2016"){FileName = " ";}
  else if(year == "2017"){FileName = "/ScaleFactors/BTaggingEfficiency/CSVv2_94XSF_V2_B_F.csv";}
  else if(year == "2018"){FileName = " ";}
  else{cout << "Please choose a year out of 2016, 2017 or 2018." << endl;}


  cout << "before file open" << endl;

  ifstream file;
  file.open(FileName);
  cout << "after file open" << endl;
  if (file.good())
  {
    string str = "";

    int line_number = 0;
cout << "before while" << endl;	
	 while(getline(file, str) && line_number != LineSpecified){
		++line_number;
	}
	if(line_number == LineSpecified){
           cout << "reading columns" << endl;
		file >> Col1;
                file >> Col2;
                file >> Col3;
                file >> Col4;
                file >> Col5;
                file >> Col6;
                file >> Col7;
                file >> Col8;
                file >> Col9;
                file >> Col10;
                file >> Col11;
         cout << "after reading columns" << endl;
       } 

   cout << "before file close" << endl;
   file.close(); 

   int Test_CSVv2OperatingPoint = 1;
   string Test_MeasurementType = "comb";
   string Test_SysType = "central";
   int Test_JetFlavour = 1;
   
   float Test_Eta = 2.9;
   float Test_Pt = 150;
   float Test_Discrim = 0.8;
   string equation = Col11;

   cout << "before if statement" << endl;
   if( Test_CSVv2OperatingPoint == Col1 && Test_MeasurementType == Col2 && Test_SysType == Col3 && Test_JetFlavour == Col4 &&
       (Test_Eta > Col5 && Test_Eta < Col6) && (Test_Eta > Col5 && Test_Eta < Col6) && (Test_Pt > Col7 && Test_Pt < Col8) &&
       (Test_Discrim > Col9 && Test_Discrim < Col10) ){

   	return equation;

   }
   else{return "0";}



 }
 else{cout << "Error reading the input file" << endl;}

}


int linecounter(const string& year){ 
   cout << "line counter" << endl;   
   int number_of_lines = 0;
   string line;

   if(year == "2016"){FileName = " ";}
   else if(year == "2017"){FileName = "/ScaleFactors/BTaggingEfficiency/CSVv2_94XSF_V2_B_F.csv";}
   else if(year == "2018"){FileName = " ";}
   else{cout << "Please input a year out of 2016, 2017 or 2018." << endl;}

   ifstream myfile(FileName);

   while (getline(myfile, line))
        ++number_of_lines;
    	return number_of_lines;

}



string BTagFileReader3(const string& year){
  cout << "before for loop" << endl;
  for(int i = 0; i < linecounter(year) + 1; i++){

	if(BTagFileReader2(i, year) != "0"){

		cout << "The values on line number " << i+2 << " were used to calculate the b-tagging SF." << endl;
		return BTagFileReader2(i, year);
	}
	else{return "0";}

  }


}


string BTagFileReader(){

  return BTagFileReader3("2017");

}
