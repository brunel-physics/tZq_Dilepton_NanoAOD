#include <fstream>
#include <limits>
#include <stdlib.h> //for atof
#include <vector> 
#include <iostream>

string process = "test";
string channel = "ee";

std::fstream& GotoLine(std::fstream& file, unsigned int num){
    file.seekg(std::ios::beg);
    for(int i=0; i < num - 1; ++i){
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }
    return file;
}


double linereader_BTagMC(const int& LineNumber, const bool& blinding, const bool& NPL, const bool& ZPlusJetsCR, const bool& ttbarCR){
    

   string BTagInputFileMC;

   if(blinding == false){
	
	if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
        	BTagInputFileMC = "BTagProbabilityValues_MC_" + process + "_" + channel + "_NPL.txt";
	}
	else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
		BTagInputFileMC = "BTagProbabilityValues_MC_" + process + "_" + channel + "_ZPlusJetsCR.txt";
	}
	else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
		BTagInputFileMC = "BTagProbabilityValues_MC_" + process + "_" + channel + "_ttbarCR.txt";
	}
	else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
		BTagInputFileMC = "BTagProbabilityValues_MC_" + process + "_" + channel + "_NPL_ZPlusJetsCR.txt";
	}
	else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
		BTagInputFileMC = "BTagProbabilityValues_MC_" + process + "_" + channel + "_NPL_ttbarCR.txt";
	}
	else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << endl;}
	else{BTagInputFileMC = "BTagProbabilityValues_MC_" + process + "_" + channel + ".txt";}

  }
  else{

	if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                BTagInputFileMC = "BTagProbabilityValues_MC_" + process + "_" + channel + "_NPL_Blinded.txt";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                BTagInputFileMC = "BTagProbabilityValues_MC_" + process + "_" + channel + "_ZPlusJetsCR_Blinded.txt";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                BTagInputFileMC = "BTagProbabilityValues_MC_" + process + "_" + channel + "_ttbarCR_Blinded.txt";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                BTagInputFileMC = "BTagProbabilityValues_MC_" + process + "_" + channel + "_NPL_ZPlusJetsCR_Blinded.txt";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                BTagInputFileMC = "BTagProbabilityValues_MC_" + process + "_" + channel + "_NPL_ttbarCR_Blinded.txt";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << endl;}
        else{BTagInputFileMC = "BTagProbabilityValues_MC_" + process + "_" + channel + "_Blinded.txt";}


 } 


   fstream file(BTagInputFileMC.c_str());
   GotoLine(file, LineNumber);

   string line;
   file >> line;

   double Value = atof(line.c_str());
   return Value;

}



double linereader_BTagData(const int& LineNumber, const bool& blinding, const bool& NPL, const bool& ZPlusJetsCR, const bool& ttbarCR){

   string BTagInputFileData;

   if(blinding == false){

        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                BTagInputFileData = "BTagProbabilityValues_Data_" + process + "_" + channel + "_NPL.txt";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                BTagInputFileData = "BTagProbabilityValues_Data_" + process + "_" + channel + "_ZPlusJetsCR.txt";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                BTagInputFileData = "BTagProbabilityValues_Data_" + process + "_" + channel + "_ttbarCR.txt";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                BTagInputFileData = "BTagProbabilityValues_Data_" + process + "_" + channel + "_NPL_ZPlusJetsCR.txt";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                BTagInputFileData = "BTagProbabilityValues_Data_" + process + "_" + channel + "_NPL_ttbarCR.txt";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << endl;}
        else{BTagInputFileData = "BTagProbabilityValues_Data_" + process + "_" + channel + ".txt";}

  }
  else{

        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                BTagInputFileData = "BTagProbabilityValues_Data_" + process + "_" + channel + "_NPL_Blinded.txt";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                BTagInputFileData = "BTagProbabilityValues_Data_" + process + "_" + channel + "_ZPlusJetsCR_Blinded.txt";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                BTagInputFileData = "BTagProbabilityValues_Data_" + process + "_" + channel + "_ttbarCR_Blinded.txt";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                BTagInputFileData = "BTagProbabilityValues_Data_" + process + "_" + channel + "_NPL_ZPlusJetsCR_Blinded.txt";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                BTagInputFileData = "BTagProbabilityValues_Data_" + process + "_" + channel + "_NPL_ttbarCR_Blinded.txt";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << endl;}
        else{BTagInputFileData = "BTagProbabilityValues_Data_" + process + "_" + channel + "_Blinded.txt";}


 } 


   fstream file(BTagInputFileData.c_str());
   GotoLine(file, LineNumber);

   string line;
   file >> line;

   double Value = atof(line.c_str());
   return Value;

}



int linecounter_BTagMC(const bool& blinding, const bool& NPL, const bool& ZPlusJetsCR, const bool& ttbarCR){ 
    
   int number_of_lines = 0;
   string line;
   string BTagInputFileMC; 

  if(blinding == false){
        
        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                BTagInputFileMC = "BTagProbabilityValues_MC_" + process + "_" + channel + "_NPL.txt";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                BTagInputFileMC = "BTagProbabilityValues_MC_" + process + "_" + channel + "_ZPlusJetsCR.txt";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                BTagInputFileMC = "BTagProbabilityValues_MC_" + process + "_" + channel + "_ttbarCR.txt";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){ 
                BTagInputFileMC = "BTagProbabilityValues_MC_" + process + "_" + channel + "_NPL_ZPlusJetsCR.txt";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){ 
                BTagInputFileMC = "BTagProbabilityValues_MC_" + process + "_" + channel + "_NPL_ttbarCR.txt";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << endl;}
        else{BTagInputFileMC = "BTagProbabilityValues_MC_" + process + "_" + channel + ".txt";}
  
  }
  else{
        
        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                BTagInputFileMC = "BTagProbabilityValues_MC_" + process + "_" + channel + "_NPL_Blinded.txt";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                BTagInputFileMC = "BTagProbabilityValues_MC_" + process + "_" + channel + "_ZPlusJetsCR_Blinded.txt";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                BTagInputFileMC = "BTagProbabilityValues_MC_" + process + "_" + channel + "_ttbarCR_Blinded.txt";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){ 
                BTagInputFileMC = "BTagProbabilityValues_MC_" + process + "_" + channel + "_NPL_ZPlusJetsCR_Blinded.txt";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){ 
                BTagInputFileMC = "BTagProbabilityValues_MC_" + process + "_" + channel + "_NPL_ttbarCR_Blinded.txt";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << endl;}
        else{BTagInputFileMC = "BTagProbabilityValues_MC_" + process + "_" + channel + "_Blinded.txt";}

 
 } 
   cout << "BTagInputFileMC = " << BTagInputFileMC << endl;

   ifstream myfile(BTagInputFileMC.c_str());

   while (getline(myfile, line))
        ++number_of_lines;
    	return number_of_lines;

}

int linecounter_BTagData(const bool& blinding, const bool& NPL, const bool& ZPlusJetsCR, const bool& ttbarCR){

   int number_of_lines = 0;
   string line;

   string BTagInputFileData;

   if(blinding == false){

        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                BTagInputFileData = "BTagProbabilityValues_Data_" + process + "_" + channel + "_NPL.txt";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                BTagInputFileData = "BTagProbabilityValues_Data_" + process + "_" + channel + "_ZPlusJetsCR.txt";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                BTagInputFileData = "BTagProbabilityValues_Data_" + process + "_" + channel + "_ttbarCR.txt";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                BTagInputFileData = "BTagProbabilityValues_Data_" + process + "_" + channel + "_NPL_ZPlusJetsCR.txt";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                BTagInputFileData = "BTagProbabilityValues_Data_" + process + "_" + channel + "_NPL_ttbarCR.txt";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << endl;}
        else{BTagInputFileData = "BTagProbabilityValues_Data_" + process + "_" + channel + ".txt";}

  }
  else{

        if(NPL == true && ZPlusJetsCR == false & ttbarCR == false){
                BTagInputFileData = "BTagProbabilityValues_Data_" + process + "_" + channel + "_NPL_Blinded.txt";
        }
        else if(NPL == false && ZPlusJetsCR == true & ttbarCR == false){
                BTagInputFileData = "BTagProbabilityValues_Data_" + process + "_" + channel + "_ZPlusJetsCR_Blinded.txt";
        }
        else if(NPL == false && ZPlusJetsCR == false & ttbarCR == true){
                BTagInputFileData = "BTagProbabilityValues_Data_" + process + "_" + channel + "_ttbarCR_Blinded.txt";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == false){
                BTagInputFileData = "BTagProbabilityValues_Data_" + process + "_" + channel + "_NPL_ZPlusJetsCR_Blinded.txt";
        }
        else if(NPL == true && ZPlusJetsCR == false & ttbarCR == true){
                BTagInputFileData = "BTagProbabilityValues_Data_" + process + "_" + channel + "_NPL_ttbarCR_Blinded.txt";
        }
        else if(NPL == true && ZPlusJetsCR == true & ttbarCR == true){cout << "Error: NPL, ZPlusJetsCR and ttbarCR cannot all be true." << endl;}
        else{BTagInputFileData = "BTagProbabilityValues_Data_" + process + "_" + channel + "_Blinded.txt";}


 }


   ifstream myfile(BTagInputFileData.c_str());

   while (getline(myfile, line))
        ++number_of_lines;
        return number_of_lines;

}



auto BTag_textfilereader2(const bool& blinding, const bool& NPL, const bool& ZPlusJetsCR, const bool& ttbarCR){

   int NumberOfLines_MCFile = linecounter_BTagMC(blinding, NPL, ZPlusJetsCR, ttbarCR);
   int NumberOfLines_DataFile = linecounter_BTagData(blinding, NPL, ZPlusJetsCR, ttbarCR);
   vector<double> BTag_WeightValues;


   cout << "NumberOfLines_MCFile = " << NumberOfLines_MCFile << endl;
   cout << "NumberOfLines_DataFile = " << NumberOfLines_DataFile << endl;

   cout << '\n' << endl;

   for(int i = 1; i < NumberOfLines_MCFile+1; i++){

	cout << "linereader_Data(i) = " << linereader_BTagData(i, blinding, NPL, ZPlusJetsCR, ttbarCR) << endl;
	cout << "linereader_MC(i) = " << linereader_BTagMC(i, blinding, NPL, ZPlusJetsCR, ttbarCR) << endl;

	double BTagWeight = linereader_BTagData(i, blinding, NPL, ZPlusJetsCR, ttbarCR) / linereader_BTagMC(i, blinding, NPL, ZPlusJetsCR, ttbarCR); 
   	
	BTag_WeightValues.push_back(BTagWeight);

   }

   cout << '\n' << endl;

   for(int i = 0; i < BTag_WeightValues.size(); i++){

	cout << "BTag_WeightValues.at(i) = " << BTag_WeightValues.at(i) << endl;

   }

   return BTag_WeightValues;

}




auto BTag_textfilereader(){

	return BTag_textfilereader2(false, false, false, false);

}
