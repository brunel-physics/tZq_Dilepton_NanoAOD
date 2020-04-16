#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;


float Col1, Col2, Col3, Col4, Col5, Col6, Col7, Col8, Col9, Col10, Col11;
string FileName;


float RowReader2(const int& LineSpecified, const bool& sigmaJER, const bool& SF, const bool& up, const bool& down) { 

  if(sigmaJER == true && SF == false && up == false && down == false){FileName = "Fall17_V3_MC_PtResolution_AK4PF.txt";}
  else if(sigmaJER == false && SF == true && up == false && down == false){FileName = "Fall17_V3_MC_SF_AK4PF.txt";}
  else if(sigmaJER == false && SF == false && up == true && down == false){FileName = "Fall17_V3_MC_SF_AK4PF.txt";}
  else if(sigmaJER == false && SF == false && up == false && down == true){FileName = "Fall17_V3_MC_SF_AK4PF.txt";}
  else{cout << "Please enter an appropriate file name" << endl;}

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
		if(sigmaJER == true && SF == false && up == false && down == false){
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
		}
		else if(sigmaJER == false && SF == true && up == false && down == false){
			file >> Col1;
                        file >> Col2;
                        file >> Col3;
                        file >> Col4;
                        file >> Col5;
                        file >> Col6;
		}
		else if(sigmaJER == false && SF == false && up == true && down == false){
                        file >> Col1;
                        file >> Col2;
                        file >> Col3;
                        file >> Col4;
                        file >> Col5;
                        file >> Col6;
                }
		else if(sigmaJER == false && SF == false && up == false && down == true){
                        file >> Col1;
                        file >> Col2;
                        file >> Col3;
                        file >> Col4;
                        file >> Col5;
                        file >> Col6;
                }
		else{cout << "Please enter an appropriate file name" << endl;}

	}

  }
 
  file.close(); 
	
  float testeta = 2.9;
  float testrho = 6.50;
  float testpt = 150;
  
  if(sigmaJER == true && SF == false && up == false && down == false){

        float answer = sqrt( Col8*abs(Col8) / (testpt*testpt)+Col9*Col9*pow(testpt,Col11)+Col10*Col10 );

  	if( (testeta > abs(Col1) && testeta < abs(Col2)) && (testrho > Col3 && testrho < Col4) && (testpt > Col6 && testpt < Col7) ){
		return answer;
  	}
  	else{return 0;}

  }
  else if(sigmaJER == false && SF == true && up == false && down == false){
 	if(testeta > Col1 && testeta < Col2){
		return Col4;
	}
	else{return 0;}
  }
  else if(sigmaJER == false && SF == false && up == true && down == false){
	if(testeta > Col1 && testeta < Col2){
        	float UpValue = Col6 - Col4;
                return UpValue;
	}
        else{return 0;}
  }
  else if(sigmaJER == false && SF == false && up == false && down == true){
        if(testeta > Col1 && testeta < Col2){
        	float DownValue = Col4 - Col5;
                return DownValue;
	}
        else{return 0;}
  }
  else{cout << "bools cannot be all true or all false" << endl;} 



} 


int linecounter(const bool& sigmaJER, const bool& SF, const bool& up, const bool& down){ 
    
   int number_of_lines = 0;
   string line;

   if(sigmaJER == true && SF == false && up == false && down == false){FileName = "Fall17_V3_MC_PtResolution_AK4PF.txt";}
   else if(sigmaJER == false && SF == true && up == false && down == false){FileName = "Fall17_V3_MC_SF_AK4PF.txt";}
   else if(sigmaJER == false && SF == false && up == true && down == false){FileName = "Fall17_V3_MC_SF_AK4PF.txt";}
   else if(sigmaJER == false && SF == false && up == false && down == true){FileName = "Fall17_V3_MC_SF_AK4PF.txt";}
   else{cout << "Please enter an appropriate file name" << endl;}

   ifstream myfile(FileName);

   while (getline(myfile, line))
        ++number_of_lines;
    	return number_of_lines;

}



float RowReader3(const bool& SigmaJER, const bool& JetSmearScaleFactor, const bool& Up, const bool& Down){

  for(int i = 0; i < linecounter(SigmaJER, JetSmearScaleFactor, Up, Down) + 1; i++){

	if(RowReader2(i, SigmaJER, JetSmearScaleFactor, Up, Down) != 0){
		string quantity; 

		if(SigmaJER == true && JetSmearScaleFactor == false && Up == false && Down == false){quantity = "sigma JER";}
   		else if(SigmaJER == false && JetSmearScaleFactor == true && Up == false && Down == false){quantity = "SF";}
   		else if(SigmaJER == false && JetSmearScaleFactor == false && Up == true && Down == false){quantity = "SF (up variation)";}
   		else if(SigmaJER == false && JetSmearScaleFactor == false && Up == false && Down == true){quantity = "SF (down variation)";}
   		else{cout << "Please enter an appropriate file name" << endl;}

		cout << "The values on line number " << i+2 << " were used to calculate the " << quantity << endl;
		return RowReader2(i, SigmaJER, JetSmearScaleFactor, Up, Down);
	}
	else{continue;}

  }


}



float RowReader(){

  bool SigmaJER = true;
  bool JetSmearScaleFactor = false;
  bool Up = false;
  bool Down = false;

  return RowReader3(SigmaJER, JetSmearScaleFactor, Up, Down);

}

