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


string TriggerSF_TextFiles;


double linereader_TriggerSF(const int& LineNumber, const string& InputTriggerSF_File){
    
   if(InputTriggerSF_File == "Data_Central"){TriggerSF_TextFiles = "TriggerSF_Efficiency_Data.txt";}
   else if(InputTriggerSF_File == "MC_Central"){TriggerSF_TextFiles = "TriggerSF_Efficiency_MC.txt";}
   else if(InputTriggerSF_File == "Data_Uncert"){TriggerSF_TextFiles = "TriggerSF_EfficiencyUncerts_Data.txt";}
   else if(InputTriggerSF_File == "MC_Uncert"){TriggerSF_TextFiles = "TriggerSF_EfficiencyUncerts_MC.txt";}
   else if(InputTriggerSF_File == "SF_Central"){TriggerSF_TextFiles = "TriggerSF_ScaleFactors.txt";}
   else if(InputTriggerSF_File == "SF_Uncert"){TriggerSF_TextFiles = "TriggerSF_ScaleFactors_Uncerts.txt";}
   else{cout << "please choose an appropriate input text file for trigger SFs" << endl;}


   using namespace std;

   fstream file(TriggerSF_TextFiles);
   GotoLine(file, LineNumber);

   string line;
   file >> line;

   double Value = atof(line.c_str());
   return Value;

}




int linecounter_TriggerSF(const string& InputTriggerSF_File){ 

   if(InputTriggerSF_File == "Data_Central"){TriggerSF_TextFiles = "TriggerSF_Efficiency_Data.txt";}
   else if(InputTriggerSF_File == "MC_Central"){TriggerSF_TextFiles = "TriggerSF_Efficiency_MC.txt";}
   else if(InputTriggerSF_File == "Data_Uncert"){TriggerSF_TextFiles = "TriggerSF_EfficiencyUncerts_Data.txt";}
   else if(InputTriggerSF_File == "MC_Uncert"){TriggerSF_TextFiles = "TriggerSF_EfficiencyUncerts_MC.txt";}
   else if(InputTriggerSF_File == "SF_Central"){TriggerSF_TextFiles = "TriggerSF_ScaleFactors.txt";}
   else if(InputTriggerSF_File == "SF_Uncert"){TriggerSF_TextFiles = "TriggerSF_ScaleFactors_Uncerts.txt";}
   else{cout << "please choose an appropriate input text file for trigger SFs" << endl;} 

 
   int number_of_lines = 0;
   string line;
   ifstream myfile(TriggerSF_TextFiles);

   while (getline(myfile, line))
        ++number_of_lines;
    	return number_of_lines;

}


auto textfilereader2_TriggerSF(const string& InputTriggerSF_File){

   int NumberOfLines = linecounter_TriggerSF(InputTriggerSF_File);
   vector<double> Value;

   for(int i = 1; i < NumberOfLines+1; i++){
   	Value.push_back(linereader_TriggerSF(i, InputTriggerSF_File));
   }

   return Value;

}


auto textfilereader_TriggerSF{[](){

  //auto vector1 = textfilereader2_TriggerSF("MC_Central");
  //auto vector1 = textfilereader2_TriggerSF("Data_Central");
  //auto vector1 = textfilereader2_TriggerSF("MC_Uncert");
  auto vector1 = textfilereader2_TriggerSF("Data_Uncert");
  //auto vector1 = textfilereader2_TriggerSF("SF_Central");
  //auto vector1 = textfilereader2_TriggerSF("SF_Uncert"); 
  //auto vector1 = textfilereader2_TriggerSF("SF_Uncert");

  //return vector1.at(0); //ee up or ee central
  //return vector1.at(1); //ee down or mumu central
  //return vector1.at(2); //mumu up or emu central
  //return vector1.at(3); //mumu down
  //return vector1.at(4); //emu up
  //return vector1.at(5); //emu down

}};

