#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;


float Col1, Col2, Col5, Col7, Col8, Col9, Col10, Col11, Col12, Col13, Col14, Col15, Col16, Col17, Col18, Col19, Col20;
string Col3, Col4, Col6, Col21, FileName;
string ZeroOutput = "0";

string RowReader2(const int& LineSpecified) { 


  FileName = "/home/eepgkkc/ScaleFactors/BTaggingEfficiency/CSVv2_94XSF_V2_B_F.csv";
  ifstream file;
  file.open(FileName);

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
                	file >> Col5;
			file >> Col6;
                	file >> Col7;
			file >> Col8;
                	file >> Col9;
			file >> Col10;
                	file >> Col11;
			file >> Col12;
                	file >> Col13;
			file >> Col14;
                	file >> Col15;
			file >> Col16;
                	file >> Col17;
			file >> Col18;
                	file >> Col19;
			file >> Col20;
                	file >> Col21;
		}
	else{cout << "Please enter an appropriate file name" << endl;}

	}

  //}
 
  file.close(); 

  int CSVWorkingPoint = 0;
  string MeasurementType = "comb";	
  string SysType = "central";
  int JetFlavour = 1;
  float testeta = 2.9;
  float testpt = 150;
  float discr = 0.88;
  

  cout << "Col1 = " << Col1 << endl;
  cout << "Col2 = " << Col2 << endl;
  cout << "Col3 = " << Col3 << endl;
  cout << "Col4 = " << Col4 << endl;
  cout << "Col5 = " << Col5 << endl;
  cout << "Col6 = " << Col6 << endl;
  cout << "Col7 = " << Col7 << endl;
  cout << "Col8 = " << Col8 << endl;
  cout << "Col9 = " << Col9 << endl;
  cout << "Col10 = " << Col10 << endl;
  cout << "Col11 = " << Col11 << endl;
  cout << "Col12 = " << Col12 << endl;
  cout << "Col13 = " << Col13 << endl;
  cout << "Col14 = " << Col14 << endl;
  cout << "Col15 = " << Col15 << endl;
  cout << "Col16 = " << Col16 << endl;
  cout << "Col17 = " << Col17 << endl;
  cout << "Col18 = " << Col18 << endl;
  cout << "Col19 = " << Col19 << endl;
  cout << "Col21 = " << Col21 << endl;


  if( (CSVWorkingPoint == Col1) && 
      (MeasurementType == Col3) &&
      (SysType == Col5) &&
      (JetFlavour == Col7) &&
      (testeta > Col9 && testeta < Col11) &&
      (testpt > Col13 && testpt < Col15) &&
      (discr > Col17 && discr < Col19)
    ){
  	return Col21;
     }
  else{return ZeroOutput;}

} 




int linecounter(){ 
    
   int number_of_lines = 0;
   string line;

   FileName = "/home/eepgkkc/ScaleFactors/BTaggingEfficiency/CSVv2_94XSF_V2_B_F.csv";
   ifstream myfile(FileName);

   while (getline(myfile, line))
        ++number_of_lines;
    	return number_of_lines;

}



string RowReader3(){

  for(int i = 0; i < linecounter() + 1; i++){

	if(RowReader2(i) != "0"){
		cout << "RowReader2(i) = " << RowReader2(i) << endl;
		return RowReader2(i);
	}
	else{cout << "Row reader 2 is equal to 0" << endl;}

  }


}



string RowReader_CSV(){

  return RowReader3();

}

