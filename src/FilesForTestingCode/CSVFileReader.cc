#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <string>
#include <algorithm>
#include <boost/algorithm/string.hpp>

int line_number = 0;
int Index = 0;
int output_index;
vector<int> line_number_vec{};

class CSVReader
{
	std::string fileName;
	std::string delimeter;
 
public:
	CSVReader(std::string filename, std::string delm = ",") :
			fileName(filename), delimeter(delm)
	{ }
 
	std::vector<std::vector<std::string> > getData();

};
 


std::vector<std::vector<std::string> > CSVReader::getData()
{
	std::ifstream file(fileName);
	std::vector<std::vector<std::string> > dataList;
	std::string line = "";
	int line_number = 0;
	
	while (getline(file, line))
	{
		line_number++;
		line_number_vec.push_back(line_number);
		std::vector<std::string> vec;
		boost::algorithm::split(vec, line, boost::is_any_of(delimeter));
		dataList.push_back(vec);
	}
	
	file.close();
	
	return dataList;

}


auto CSVFileReader()
{
	CSVReader reader("/home/eepgkkc/ScaleFactors/BTaggingEfficiency/CSVv2_94XSF_V2_B_F.csv");
	std::vector<std::vector<std::string> > dataList = reader.getData();
        
	vector<string> OutputVec{}; 
	vector<string> outputstringvec{};

	vector<string> CSVv2OperatingPointTest{"3", "1", "2", "0", "2", "0", "0", "1"}; //no spaces after start
  	vector<string> MeasurementTypeTest{"iterativefit", "comb", "mujets", "comb", "mujets", "comb", "mujets", "comb"}; //one space after start
 	vector<string> SysTypeTest{"central", "down", "up", "central", "down", "down", "up", "down"}; // one space after start

	vector<string> JetFlavourTest{"2", "0", "1", "1", "0", "1", "0", "0"}; //one space after start
        vector<string> EtaTest{"1.7", "1.7", "2.0", "0", "-2.4", "-1.3", "-0.04", "1.7"}; 
        vector<string> PtTest{"150", "25", "601", "880", "22", "250", "29.3", "25"}; 
        vector<string> DiscrTest{"-14", "0.4", "0.95", "0.5", "0.8", "0.722", "0.33", "0.4"}; 

	vector<string> OutVec{};
	vector<string> FinalOutVec{};

for(int i = 0; i < CSVv2OperatingPointTest.size(); i++){
	
	for(std::vector<std::string> vec : dataList)
	{
		for(std::string data : vec)
		{	
			OutputVec.push_back(data);

		}
	}


		for(std::vector<std::string> vec : dataList)
        	{
                	for(std::string data : vec)
                	{
				
				string VecAt1String = vec.at(1);
				string VecAt2String = vec.at(2);
				string VecAt3String = vec.at(3);
			        string VecAt4String = vec.at(4);
				string VecAt5String = vec.at(5);
				string VecAt6String = vec.at(6);
				string VecAt7String = vec.at(7);
				string VecAt8String = vec.at(8);
				string VecAt9String = vec.at(9);
				

				VecAt1String.erase(0, 0);
				VecAt2String.erase(0, 0);
                                VecAt3String.erase(0, 0);
                                VecAt4String.erase(0, 0);
				VecAt1String.erase(remove(VecAt1String.begin(), VecAt1String.end(), ' '), VecAt1String.end());		
				VecAt2String.erase(remove(VecAt2String.begin(), VecAt2String.end(), ' '), VecAt2String.end());
				VecAt3String.erase(remove(VecAt3String.begin(), VecAt3String.end(), ' '), VecAt3String.end());
                                VecAt4String.erase(remove(VecAt4String.begin(), VecAt4String.end(), ' '), VecAt4String.end());	



				float VecAt4Float = stof(VecAt4String);
                        	float VecAt5Float = stof(VecAt5String);
				float VecAt6Float = stof(VecAt6String);
				float VecAt7Float = stof(VecAt7String);
				float VecAt8Float = stof(VecAt8String);
                        	float VecAt9Float = stof(VecAt9String);

				float PtTestFloat = stof(PtTest.at(i));
				float EtaTestFloat = stof(EtaTest.at(i));
				float DiscrTestFloat = stof(DiscrTest.at(i));
			

				if( (vec.at(0) == CSVv2OperatingPointTest.at(i)) 
    				&& (VecAt1String == MeasurementTypeTest.at(i))
    				&& (VecAt2String == SysTypeTest.at(i))
    	  			&& (VecAt3String == JetFlavourTest.at(i))
    				&& (VecAt4Float < EtaTestFloat)
    				&& (VecAt5Float > EtaTestFloat)
    	  			&& (VecAt6Float < PtTestFloat)     
    				&& (VecAt7Float > PtTestFloat)
    				&& (VecAt8Float < DiscrTestFloat) 
    				&& (VecAt9Float > DiscrTestFloat)
  				){
					OutVec.push_back(vec.at(10));					
				}
				else if( (vec.at(0) != CSVv2OperatingPointTest.at(i))
                        	|| (vec.at(1) != MeasurementTypeTest.at(i))
                        	|| (vec.at(2) != SysTypeTest.at(i))
                        	|| (vec.at(3) != JetFlavourTest.at(i))
                        	|| (VecAt4Float > EtaTestFloat)
                        	|| (VecAt5Float < EtaTestFloat)
                        	|| (VecAt6Float > PtTestFloat)     
                        	|| (VecAt7Float < PtTestFloat)
                        	|| (VecAt8Float > DiscrTestFloat) 
                        	|| (VecAt9Float < DiscrTestFloat)){OutVec.push_back("0");}
				else{cout << "double check criteria" << endl;}
			
                	}

        	}




	vector<string> NewOutVec{};
	vector<string> Zeroes{}; 
	Zeroes.push_back("0");
	Zeroes.push_back("0");

	bool check = all_of(OutVec.begin(), OutVec.end(), [](string s){return s == "0";});


	if(OutVec.size() != 0 && check == false){
		for(int i = 0; i < OutVec.size(); i++){

			if(OutVec.at(i) != "0"){NewOutVec.push_back(OutVec.at(i));}
	
		}
	
		string outputstring;
	
		if(NewOutVec.size() > 11){
			outputstring = NewOutVec.at( ((i+1)*11)-1 );
		}
		else{
			outputstring = NewOutVec.at(0); 
		}

		outputstring.erase(outputstring.begin()+1);
                outputstring.erase(outputstring.begin());
                outputstring.erase(outputstring.end()-2);
                outputstring.erase(outputstring.end()-1);
		
                string::size_type pos = 0;
 
                while ((pos = outputstring.find('x', pos)) != string::npos)
                {
                        outputstring.replace(pos, 1, PtTest.at(i));
                        pos += 2;
                }
                
		outputstringvec.push_back(outputstring);
		FinalOutVec.push_back(outputstringvec.at(i));

	
	}
	else{FinalOutVec.push_back(Zeroes.at(0));}



}//end of for loop



//Evaluating the mathematical expression in the string
string ConcatenatedString, ConcatenatedString2, ConcatenatedString3, ConcatenatedString4;
string ConcatenatedString5, ConcatenatedString6, ConcatenatedString7, ConcatenatedString8;
string ConcatenatedString9, ConcatenatedString10, ConcatenatedString11, ConcatenatedString12;
string ConcatenatedString13, ConcatenatedString14, ConcatenatedString15, ConcatenatedString16;
vector<char> VecForConcString{};
vector<char> VecForConcString2{};
vector<char> VecForConcString3{};
vector<char> VecForConcString4{};
vector<char> VecForConcString5{};
vector<char> VecForConcString6{};
vector<char> VecForConcString7{};
vector<char> VecForConcString8{};
vector<char> VecForConcString9{};
vector<char> VecForConcString10{};
vector<char> VecForConcString11{};
vector<char> VecForConcString12{};
vector<char> VecForConcString13{};
vector<char> VecForConcString14{};
vector<char> VecForConcString15{};
vector<char> VecForConcString16{};
vector<float> ResultVector{};
int index, index2, index3, index4, index5, index6, index7, index8, index9, index10, index11, index12, index13, index14, index15, index16;
float result;

for(int i = 0; i < FinalOutVec.size(); i++){

	string FirstElement = FinalOutVec.at(i);

	if(FirstElement.at(0) != '('){

		int LastIndex;

		//first
		for(int i = 0; i < FirstElement.length(); i++){

			if(FirstElement.at(i) != ')'&& 
	   	   	   FirstElement.at(i) != '(' &&
	   	   	   FirstElement.at(i) != '*' &&
	   	   	   FirstElement.at(i) != '/' &&
	   	   	   FirstElement.at(i) != '+' &&
	  	   	   FirstElement.at(i) != '-'){VecForConcString.push_back(FirstElement.at(i)); LastIndex = i;}
			else if(i == 0 && FirstElement.at(i) == '('){continue;}
			else{index = i; break;}
		}


		for(int i = 0; i < VecForConcString.size(); i++){
			if(i == 0){ConcatenatedString = VecForConcString.at(i);}
			else{ConcatenatedString += VecForConcString.at(i);}

		}	
		
		float ConcatenatedStringToFloat = stof(ConcatenatedString);
	
		if(LastIndex == FirstElement.length()-1){ResultVector.push_back(ConcatenatedStringToFloat);}
		else{
	
			int Min1;

			if(FirstElement.at(index) == '+' && 
			   FirstElement.at(index+1) == '(' && 
			   FirstElement.at(index+2) == '-' && 
			   FirstElement.at(index+3) == '('){
			
				Min1 = index+4;
			}
			else if(FirstElement.at(index) == '+' &&
                           	FirstElement.at(index+1) == '(' && 
                          	FirstElement.at(index+2) == '(' && 
                           	FirstElement.at(index+3) == '-' &&
				FirstElement.at(index+2) == '('){
                        
                                	Min1 = index+5;
                        }
			else if(FirstElement.at(index) == '*' &&
				FirstElement.at(index+1) == '(' &&
				FirstElement.at(index+2) == '('){
			
					Min1 = index+3;

			}
			else if(FirstElement.at(index) == '+' &&
				FirstElement.at(index+1) == '-'){
			
					Min1 = index+2;

			}
			else if(FirstElement.at(index) == '+' &&
				FirstElement.at(index+1) == '(' &&
				FirstElement.at(index+2) == '(' &&
				FirstElement.at(index+3) == '-' &&
				FirstElement.at(index+4) == '('){

				Min1 = index+5;
			

			}
			else{Min1 = index+1;}

			//second
			for(int i = Min1; i < FirstElement.length(); i++){

        			if(FirstElement.at(i) == 'e' && FirstElement.at(i+1) == '-'){
					
					VecForConcString2.push_back(FirstElement.at(i));
					VecForConcString2.push_back(FirstElement.at(i+1));
					VecForConcString2.push_back(FirstElement.at(i+2));
					VecForConcString2.push_back(FirstElement.at(i+3));
					
					index2 = i+4;
					break;

				}
				else if(FirstElement.at(i) != ')'&&
           	  		   	FirstElement.at(i) != '(' &&
           	  		   	FirstElement.at(i) != '*' &&
           	  		   	FirstElement.at(i) != '/' &&
           	  		   	FirstElement.at(i) != '+' &&
                  		   	FirstElement.at(i) != '-'){VecForConcString2.push_back(FirstElement.at(i));}
               			else{index2 = i; break;}
			
			}



			for(int i = 0; i < VecForConcString2.size(); i++){
        			if(i == 0){ConcatenatedString2 = VecForConcString2.at(i);}
        			else{ConcatenatedString2 += VecForConcString2.at(i);}

			}

			float ConcatenatedStringToFloat2 = stof(ConcatenatedString2);

			int Min2;

			if(FirstElement.at(index2) == '*' && 
			   FirstElement.at(index2+1) == '(' && 
			   FirstElement.at(index2+2) == 'l' && 
			   FirstElement.at(index2+3) == 'o' &&
			   FirstElement.at(index2+4) == 'g' &&
			   FirstElement.at(index2+5) == '('){Min2 = index2+6;}
                        else{Min2 = index2+2;}

			//third
			for(int i = Min2; i < FirstElement.length(); i++){

        			if(FirstElement.at(i) != ')'&&
           	   		   FirstElement.at(i) != '(' &&
           	   		   FirstElement.at(i) != '*' &&
           	   		   FirstElement.at(i) != '/' &&
           	   		   FirstElement.at(i) != '+' &&
           	   		   FirstElement.at(i) != '-'){VecForConcString3.push_back(FirstElement.at(i));}
        			else{index4 = i; break;}
			}


			for(int i = 0; i < VecForConcString3.size(); i++){
        			if(i == 0){ConcatenatedString3 = VecForConcString3.at(i);}
        			else{ConcatenatedString3 += VecForConcString3.at(i);}

			}

			float ConcatenatedStringToFloat3 = stof(ConcatenatedString3);

			int Min3;

			if(FirstElement.at(index4) == '+'){Min3 = index4+1;}
			else{Min3 = index4+1;}


			//fourth
			for(int i = Min3; i < FirstElement.length(); i++){

        			if(FirstElement.at(i) != ')'&&
           	   		   FirstElement.at(i) != '(' &&
           	   		   FirstElement.at(i) != '*' &&
           	   		   FirstElement.at(i) != '/' &&
           	   	           FirstElement.at(i) != '+' &&
           	   		   FirstElement.at(i) != '-'){VecForConcString4.push_back(FirstElement.at(i));}
        			else{index5 = i; break;}
			}


			for(int i = 0; i < VecForConcString4.size(); i++){
        			if(i == 0){ConcatenatedString4 = VecForConcString4.at(i);}
        			else{ConcatenatedString4 += VecForConcString4.at(i);}

			}

			float ConcatenatedStringToFloat4 = stof(ConcatenatedString4);

			int Min4;

			if(FirstElement.at(index5) == ')' && 
			   FirstElement.at(index5+1) == '*' && 
			   FirstElement.at(index5+2) == '(' && 
			   FirstElement.at(index5+3) == 'l' && 
			   FirstElement.at(index5+4) == 'o' && 
			   FirstElement.at(index5+5) == 'g' && 
			   FirstElement.at(index5+6) == '('){Min4 = index5+7;}
                        else{Min4 = index5+2;}


			//fifth
			for(int i = Min4; i < FirstElement.length(); i++){

        			if(FirstElement.at(i) != ')'&&
           	   		   FirstElement.at(i) != '(' &&
           	   		   FirstElement.at(i) != '*' &&
           	   		   FirstElement.at(i) != '/' &&
           	   		   FirstElement.at(i) != '+' &&
           	   		   FirstElement.at(i) != '-'){VecForConcString5.push_back(FirstElement.at(i));}
        			else{index6 = i; break;}
			}


			for(int i = 0; i < VecForConcString5.size(); i++){
        			if(i == 0){ConcatenatedString5 = VecForConcString5.at(i);}
        			else{ConcatenatedString5 += VecForConcString5.at(i);}

			}

			float ConcatenatedStringToFloat5 = stof(ConcatenatedString5);

			int Min5;

			if(FirstElement.at(index6) == '+'){Min5 = index6+1;}
			else{Min5 = index6+1;}

			//sixth
			for(int i = Min5; i < FirstElement.length(); i++){

        			if(FirstElement.at(i) != ')'&&
           	   		   FirstElement.at(i) != '(' &&
           	   		   FirstElement.at(i) != '*' &&
           	   		   FirstElement.at(i) != '/' &&
           	   		   FirstElement.at(i) != '+' &&
           	   		   FirstElement.at(i) != '-'){VecForConcString6.push_back(FirstElement.at(i));}
        			else{index7 = i; break;}
			}


			for(int i = 0; i < VecForConcString6.size(); i++){
        			if(i == 0){ConcatenatedString6 = VecForConcString6.at(i);}
        			else{ConcatenatedString6 += VecForConcString6.at(i);}

			}
	
			float ConcatenatedStringToFloat6 = stof(ConcatenatedString6);

			int Min6;

                        if(FirstElement.at(index8) == ')' &&
                           FirstElement.at(index8+1) == '*' &&
                           FirstElement.at(index8+2) == '('){Min6 = index8+3;}
                        else{Min6 = index7+3;}
			

			//seventh
			for(int i = Min6; i < FirstElement.length(); i++){

        			if(FirstElement.at(i) != ')'&&
           	   		   FirstElement.at(i) != '(' &&
           	   		   FirstElement.at(i) != '*' &&
           	   		   FirstElement.at(i) != '/' &&
           	   		   FirstElement.at(i) != '+' &&
           	   		   FirstElement.at(i) != '-'){VecForConcString7.push_back(FirstElement.at(i));}
        			else{index8 = i; break;}
			}


			for(int i = 0; i < VecForConcString7.size(); i++){
        			if(i == 0){ConcatenatedString7 = VecForConcString7.at(i);}
        			else{ConcatenatedString7 += VecForConcString7.at(i);}

			}

			float ConcatenatedStringToFloat7 = stof(ConcatenatedString7);


			int Min7;
                        
			if(FirstElement.at(index8) == '-' &&
                           FirstElement.at(index8+1) == '(' &&
                           FirstElement.at(index8+2) == '-' &&
			   FirstElement.at(index8+3) == '('){Min7 = index8+4;}
                        else{Min7 = index8+1;}

			//eighth
			for(int i = Min7; i < FirstElement.length(); i++){

        			if(FirstElement.at(i) != ')'&&
           	   		   FirstElement.at(i) != '(' &&
           	   		   FirstElement.at(i) != '*' &&
           	   		   FirstElement.at(i) != '/' &&
           	   		   FirstElement.at(i) != '+' &&
           	   		   FirstElement.at(i) != '-'){VecForConcString8.push_back(FirstElement.at(i));}
        			else{index9 = i; break;}
			}     


			for(int i = 0; i < VecForConcString8.size(); i++){
        			if(i == 0){ConcatenatedString8 = VecForConcString8.at(i);}
        			else{ConcatenatedString8 += VecForConcString8.at(i);}
    
			}   

			float ConcatenatedStringToFloat8 = stof(ConcatenatedString8);

			int Min8;

			if(FirstElement.at(index9) == '*' &&
			   FirstElement.at(index9+1) == 'l' &&
                           FirstElement.at(index9+2) == 'o' &&
                           FirstElement.at(index9+3) == 'g' &&
                           FirstElement.at(index9+4) == '('){Min8 = index9+5;}
			else{Min8 = index9+2;}

			//ninth
			for(int i = Min8; i < FirstElement.length(); i++){

        			if(FirstElement.at(i) != ')'&&
           	   		   FirstElement.at(i) != '(' &&
           	   		   FirstElement.at(i) != '*' &&
           	   		   FirstElement.at(i) != '/' &&
           	   		   FirstElement.at(i) != '+' &&
           	   		   FirstElement.at(i) != '-'){VecForConcString9.push_back(FirstElement.at(i));}
        			else{index10 = i; break;}
			}


			for(int i = 0; i < VecForConcString9.size(); i++){
        			if(i == 0){ConcatenatedString9 = VecForConcString9.at(i);}
        			else{ConcatenatedString9 += VecForConcString9.at(i);}

			}

			float ConcatenatedStringToFloat9 = stof(ConcatenatedString9);

			int Min9;

			if(FirstElement.at(index10) == '+'){Min9 = index10+1;}
			else{Min9 = index10 + 1;}


			int LastIndex2;

			//tenth
			for(int i = Min9; i < FirstElement.length(); i++){

        			if(FirstElement.at(i) != ')'&&
           	   		   FirstElement.at(i) != '(' &&
           	   		   FirstElement.at(i) != '*' &&
           	   		   FirstElement.at(i) != '/' &&
           	   		   FirstElement.at(i) != '+' &&
           	   		   FirstElement.at(i) != '-'){VecForConcString10.push_back(FirstElement.at(i)); LastIndex2 = i;}
        			else{index11 = i; break;}
			}


			for(int i = 0; i < VecForConcString10.size(); i++){
        			if(i == 0){ConcatenatedString10 = VecForConcString10.at(i);}
        			else{ConcatenatedString10 += VecForConcString10.at(i);}

			}
		
			float ConcatenatedStringToFloat10 = stof(ConcatenatedString10);


			//result for an equation containing 10 floats and ending in "))))))))' 
			 if(LastIndex2 == FirstElement.length()-9 &&
                            FirstElement.at(index) == '+' &&
                            FirstElement.at(index+1) == '(' &&
                            FirstElement.at(index+2) == '-' &&
                            FirstElement.at(index+3) == '('){
				

			    	result = ConcatenatedStringToFloat+(-(ConcatenatedStringToFloat2*(log(ConcatenatedStringToFloat3+ConcatenatedStringToFloat4)*(log(ConcatenatedStringToFloat5+ConcatenatedStringToFloat6)*(ConcatenatedStringToFloat7-(-(ConcatenatedStringToFloat8*log(ConcatenatedStringToFloat9+ConcatenatedStringToFloat10))))))));

                                ResultVector.push_back(result);

			}
                        else{


				int Min10;

				if(FirstElement.at(index11+8) == '-'){Min10 = index11+9;}
				else{Min10 = index11+3;}


				int LastIndex3;

				//eleventh
				for(int i = Min10; i < FirstElement.length(); i++){

        				if(FirstElement.at(i) != ')'&&
           	   		   	   FirstElement.at(i) != '(' &&
           	   		  	   FirstElement.at(i) != '*' &&
           	   		   	   FirstElement.at(i) != '/' &&
           	   		   	   FirstElement.at(i) != '+' &&
           	   		   	   FirstElement.at(i) != '-'){VecForConcString11.push_back(FirstElement.at(i)); LastIndex3 = i;}
        				else{index12 = i; break;}
				}


				for(int i = 0; i < VecForConcString11.size(); i++){
        				if(i == 0){ConcatenatedString11 = VecForConcString11.at(i);}
        				else{ConcatenatedString11 += VecForConcString11.at(i);}

				}

				float ConcatenatedStringToFloat11 = stof(ConcatenatedString11);


				//for an equation with 11 floats
				if(LastIndex == FirstElement.length()-1 && FirstElement.at(index11+8) == '-'){

					result = ConcatenatedStringToFloat+((-(ConcatenatedStringToFloat2*exp(-5)*(log(ConcatenatedStringToFloat3+ConcatenatedStringToFloat4)*(log(ConcatenatedStringToFloat5+ConcatenatedStringToFloat6)*(ConcatenatedStringToFloat7-(-(ConcatenatedStringToFloat8*log(ConcatenatedStringToFloat9+ConcatenatedStringToFloat10))))))))-ConcatenatedStringToFloat11);

					ResultVector.push_back(result);

				}
				else if(FirstElement.at(index) == '+' &&
                                	FirstElement.at(index+1) == '(' &&
                                	FirstElement.at(index+2) == '(' &&
                                	FirstElement.at(index+3) == '-' &&
                                	FirstElement.at(index+4) == '(' &&
					FirstElement.at(index11+8) == '-'){

                                	result = ConcatenatedStringToFloat+((-(ConcatenatedStringToFloat2*(log(ConcatenatedStringToFloat3+ConcatenatedStringToFloat4)*(log(ConcatenatedStringToFloat5+ConcatenatedStringToFloat6)*(ConcatenatedStringToFloat7-(-(ConcatenatedStringToFloat8*log(ConcatenatedStringToFloat9+ConcatenatedStringToFloat10))))))))-ConcatenatedStringToFloat11);

					ResultVector.push_back(result);
                        
                        	}
				else{
				
					//twelfth
					for(int i = index12+1; i < FirstElement.length(); i++){

        					if(FirstElement.at(i) != ')'&&
           	   		   	   	   FirstElement.at(i) != '(' &&
           	   		   	   	   FirstElement.at(i) != '*' &&
           	   		   	   	   FirstElement.at(i) != '/' &&
           	   		   	   	   FirstElement.at(i) != '+' &&
           	   		   	   	   FirstElement.at(i) != '-'){VecForConcString12.push_back(FirstElement.at(i));}
        					else{index13 = i; break;}
					}


					for(int i = 0; i < VecForConcString12.size(); i++){
        					if(i == 0){ConcatenatedString12 = VecForConcString12.at(i);}
        					else{ConcatenatedString12 += VecForConcString12.at(i);}

					}

					float ConcatenatedStringToFloat12 = stof(ConcatenatedString12);

					//thirteenth
					for(int i = index13+1; i < FirstElement.length(); i++){

        					if(FirstElement.at(i) != ')'&&
           	   		   	   	   FirstElement.at(i) != '(' &&
           	   		   	   	   FirstElement.at(i) != '*' &&
           	   		   	   	   FirstElement.at(i) != '/' &&
           	   		   	   	   FirstElement.at(i) != '+' &&
           	   		   	   	   FirstElement.at(i) != '-'){VecForConcString13.push_back(FirstElement.at(i));}
        					else{index14 = i; break;}
					}


					for(int i = 0; i < VecForConcString13.size(); i++){
        					if(i == 0){ConcatenatedString13 = VecForConcString13.at(i);}
        					else{ConcatenatedString13 += VecForConcString13.at(i);}

					}

					float ConcatenatedStringToFloat13 = stof(ConcatenatedString13);

					//Calculating the result for an equation containing 13 floats

					if(FirstElement.at(index) == '+'){result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2;}
					else if(FirstElement.at(index) == '-'){result = ConcatenatedStringToFloat - ConcatenatedStringToFloat2;}
					else if(FirstElement.at(index) == '*'){result = ConcatenatedStringToFloat * ConcatenatedStringToFloat2;}
					else{cout << "FirstElement.at(index) is " << FirstElement.at(index) << ". This is not a +, - or *. Output has been set to zero" << endl; result = 0;}


					if(FirstElement.at(index2) == '*' && FirstElement.at(index2+1) == '(' && FirstElement.at(index2+2) == '-'){
						result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3);
					}
					else{cout << "error" << endl;}


					if(FirstElement.at(index4) == '+'){
						result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + ConcatenatedStringToFloat4);
					}
					else{cout << "error message 2" << endl;}


					if(FirstElement.at(index5) == '*' && FirstElement.at(index5+1) == '('){
						result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + ConcatenatedStringToFloat4*(ConcatenatedStringToFloat5));
					}
					else{cout << "error message 3" << endl;}

					if(FirstElement.at(index6) == '+'){
        					result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + ConcatenatedStringToFloat4*(ConcatenatedStringToFloat5 + ConcatenatedStringToFloat6));
					}
					else{cout << "error message 4" << endl;}


					if(FirstElement.at(index7) == '*'){
        					result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + ConcatenatedStringToFloat4*(ConcatenatedStringToFloat5 + ConcatenatedStringToFloat6*(-ConcatenatedStringToFloat7)));
					}   
					else{cout << "error message 4" << endl;}


					if(FirstElement.at(index8) == '+'){
 	       					result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + ConcatenatedStringToFloat4*(ConcatenatedStringToFloat5 + ConcatenatedStringToFloat6*(-ConcatenatedStringToFloat7+ConcatenatedStringToFloat8)));
					}
					else{cout << "error message 5" << endl;}

					if(FirstElement.at(index9) == '*' && FirstElement.at(index9+1) == '('){
						result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + ConcatenatedStringToFloat4*(ConcatenatedStringToFloat5 + ConcatenatedStringToFloat6*(-ConcatenatedStringToFloat7+ConcatenatedStringToFloat8*(ConcatenatedStringToFloat9))));
					}
					else{cout << "error message 6" << endl;}

					if(FirstElement.at(index10) == '+'){
        					result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + ConcatenatedStringToFloat4*(ConcatenatedStringToFloat5 + ConcatenatedStringToFloat6*(-ConcatenatedStringToFloat7+ConcatenatedStringToFloat8*(ConcatenatedStringToFloat9+ConcatenatedStringToFloat10))));
					}
					else{cout << "error message 7" << endl;}

					if(FirstElement.at(index11) == '*' && FirstElement.at(index11+1) == '(' && FirstElement.at(index11+2) == '-'){
        					result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + ConcatenatedStringToFloat4*(ConcatenatedStringToFloat5 + ConcatenatedStringToFloat6*(-ConcatenatedStringToFloat7+ConcatenatedStringToFloat8*(ConcatenatedStringToFloat9+ConcatenatedStringToFloat10*(-ConcatenatedStringToFloat11)))));
					}
					else{cout << "error message 8" << endl;}


					if(FirstElement.at(index12) == '+'){
        					result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + ConcatenatedStringToFloat4*(ConcatenatedStringToFloat5 + ConcatenatedStringToFloat6*(-ConcatenatedStringToFloat7+ConcatenatedStringToFloat8*(ConcatenatedStringToFloat9+ConcatenatedStringToFloat10*(-ConcatenatedStringToFloat11+ConcatenatedStringToFloat12)))));
					}	
					else{cout << "error message 9" << endl;}

					if(FirstElement.at(index13) == '*'){
        					result = ConcatenatedStringToFloat + ConcatenatedStringToFloat2*(-ConcatenatedStringToFloat3 + ConcatenatedStringToFloat4*(ConcatenatedStringToFloat5 + ConcatenatedStringToFloat6*(-ConcatenatedStringToFloat7+ConcatenatedStringToFloat8*(ConcatenatedStringToFloat9+ConcatenatedStringToFloat10*(-ConcatenatedStringToFloat11+ConcatenatedStringToFloat12*ConcatenatedStringToFloat13)))));
					}
					else{cout << "error message 10" << endl;}

					ResultVector.push_back(result);

				}

			}

		}

	}
	else{
		//first
        	for(int i = 1; i < FirstElement.length(); i++){

                	if(FirstElement.at(i) != ')'&&
                   	   FirstElement.at(i) != '(' &&
                   	   FirstElement.at(i) != '*' &&
                   	   FirstElement.at(i) != '/' &&
                   	   FirstElement.at(i) != '+' &&
                   	   FirstElement.at(i) != '-'){VecForConcString.push_back(FirstElement.at(i));}
                	else{index = i; break;}
        	}


        	for(int i = 0; i < VecForConcString.size(); i++){
                	if(i == 0){ConcatenatedString = VecForConcString.at(i);}
                	else{ConcatenatedString += VecForConcString.at(i);}

        	}

        	float ConcatenatedStringToFloat = stof(ConcatenatedString);

		int Minimum;

		if(FirstElement.at(index) == '+' &&
		   FirstElement.at(index+1) == '(' &&
		   FirstElement.at(index+2) == '-' &&
		   FirstElement.at(index+3) == '('){

			Minimum = index+4;

		}
		else{Minimum = index+3;}


		//second	
		for(int i = Minimum; i < FirstElement.length(); i++){

                	if(FirstElement.at(i) != ')'&&
                   	   FirstElement.at(i) != '(' &&
                   	   FirstElement.at(i) != '*' &&
                   	   FirstElement.at(i) != '/' &&
                   	   FirstElement.at(i) != '+' &&
                   	   FirstElement.at(i) != '-'){VecForConcString2.push_back(FirstElement.at(i));}
                	else{index2 = i; break;}
        	}


        	for(int i = 0; i < VecForConcString2.size(); i++){
                	if(i == 0){ConcatenatedString2 = VecForConcString2.at(i);}
                	else{ConcatenatedString2 += VecForConcString2.at(i);}

        	}
        	
		float ConcatenatedStringToFloat2 = stof(ConcatenatedString2);

		int Minimum2;

		if(FirstElement.at(index2) == '*' && FirstElement.at(index2+1) == '('){Minimum2 = index2+ 6;}
		else{Minimum2 = index2+2;}


		//third
		for(int i = Minimum2; i < FirstElement.length(); i++){

                	if(FirstElement.at(i) != ')'&&
                   	   FirstElement.at(i) != '(' &&
                   	   FirstElement.at(i) != '*' &&
                   	   FirstElement.at(i) != '/' &&
                   	   FirstElement.at(i) != '+' &&
                   	   FirstElement.at(i) != '-'){VecForConcString3.push_back(FirstElement.at(i));}
                	else{index3 = i; break;}
        	}


        	for(int i = 0; i < VecForConcString3.size(); i++){
                	if(i == 0){ConcatenatedString3 = VecForConcString3.at(i);}
                	else{ConcatenatedString3 += VecForConcString3.at(i);}

        	}

        	float ConcatenatedStringToFloat3 = stof(ConcatenatedString3);

		//fourth
		for(int i = index3+1; i < FirstElement.length(); i++){

                	if(FirstElement.at(i) != ')'&&
                   	   FirstElement.at(i) != '(' &&
                   	   FirstElement.at(i) != '*' &&
                   	   FirstElement.at(i) != '/' &&
                   	   FirstElement.at(i) != '+' &&
                   	   FirstElement.at(i) != '-'){VecForConcString4.push_back(FirstElement.at(i));}
                	else{index4 = i; break;}
        	}


        	for(int i = 0; i < VecForConcString4.size(); i++){
                	if(i == 0){ConcatenatedString4 = VecForConcString4.at(i);}
                	else{ConcatenatedString4 += VecForConcString4.at(i);}

        	}

        	float ConcatenatedStringToFloat4 = stof(ConcatenatedString4);

		int Minimum3;


		if(FirstElement.at(index4) == ')' && FirstElement.at(index4+1) == '*' && FirstElement.at(index4+2) == '('){Minimum3 = index4+7;}
		else{Minimum3 = index4+4;}


		//fifth
		for(int i = Minimum3; i < FirstElement.length(); i++){

                	if(FirstElement.at(i) != ')'&&
                   	   FirstElement.at(i) != '(' &&
                   	   FirstElement.at(i) != '*' &&
                   	   FirstElement.at(i) != '/' &&
                   	   FirstElement.at(i) != '+' &&
                   	   FirstElement.at(i) != '-'){VecForConcString5.push_back(FirstElement.at(i));}
               	 	else{index5 = i; break;}
        	}


        	for(int i = 0; i < VecForConcString5.size(); i++){
                	if(i == 0){ConcatenatedString5 = VecForConcString5.at(i);}
                	else{ConcatenatedString5 += VecForConcString5.at(i);}

        	}

        	float ConcatenatedStringToFloat5 = stof(ConcatenatedString5);

		int Minimum4;
		
		if(FirstElement.at(index4) == ')' && FirstElement.at(index4+1) == '*' && FirstElement.at(index4+2) == '('){Minimum4 = index5+1;}
		else{Minimum4 = index5+2;}	

		
		//sixth
		for(int i = Minimum4; i < FirstElement.length(); i++){

                	if(FirstElement.at(i) != ')'&&
                   	   FirstElement.at(i) != '(' &&
                   	   FirstElement.at(i) != '*' &&
                   	   FirstElement.at(i) != '/' &&
                   	   FirstElement.at(i) != '+' &&
                   	   FirstElement.at(i) != '-'){VecForConcString6.push_back(FirstElement.at(i));}
                	else{index6 = i; break;}
        	}


        	for(int i = 0; i < VecForConcString6.size(); i++){
                	if(i == 0){ConcatenatedString6 = VecForConcString6.at(i);}
                	else{ConcatenatedString6 += VecForConcString6.at(i);}

        	}

        	float ConcatenatedStringToFloat6 = stof(ConcatenatedString6);
	
		//seventh
		for(int i = index6+1; i < FirstElement.length(); i++){

                	if(FirstElement.at(i) != ')'&&
                   	   FirstElement.at(i) != '(' &&
                   	   FirstElement.at(i) != '*' &&
                   	   FirstElement.at(i) != '/' &&
                   	   FirstElement.at(i) != '+' &&
                   	   FirstElement.at(i) != '-'){VecForConcString7.push_back(FirstElement.at(i));}
                	else{index7 = i; break;}
        	}


        	for(int i = 0; i < VecForConcString7.size(); i++){
                	if(i == 0){ConcatenatedString7 = VecForConcString7.at(i);}
                	else{ConcatenatedString7 += VecForConcString7.at(i);}

        	}

        	float ConcatenatedStringToFloat7 = stof(ConcatenatedString7);


		int MinValue;


		if(FirstElement.at(index) == '+' &&
		   FirstElement.at(index+1) == '(' &&
		   FirstElement.at(index+2) == '-' &&
		   FirstElement.at(index+3) == '('){
		
			MinValue = index7+7;

		}
		else{MinValue = index7+5;}

		//seventh
        	for(int i = MinValue; i < FirstElement.length(); i++){

                	if(FirstElement.at(i) != ')'&&
                   	   FirstElement.at(i) != '(' &&
                  	   FirstElement.at(i) != '*' &&
                   	   FirstElement.at(i) != '/' &&
                   	   FirstElement.at(i) != '+' &&
                   	   FirstElement.at(i) != '-'){VecForConcString8.push_back(FirstElement.at(i));}
                	else{index8 = i; break;}
        	}


        	for(int i = 0; i < VecForConcString8.size(); i++){
                	if(i == 0){ConcatenatedString8 = VecForConcString8.at(i);}
                	else{ConcatenatedString8 += VecForConcString8.at(i);}

        	}

        	float ConcatenatedStringToFloat8 = stof(ConcatenatedString8);


		//result

		if(FirstElement.at(index) == '*' && 
		   FirstElement.at(index+1) == '(' && 
		   FirstElement.at(index+2) == '(' &&
		   FirstElement.at(index2) == '+' && 
		   FirstElement.at(index2+1) == '(' &&
		   FirstElement.at(index3) == '*' &&
		   FirstElement.at(index4) == ')' && 
		   FirstElement.at(index4+1) == ')' && 
		   FirstElement.at(index4+2) == '/' && 
		   FirstElement.at(index4+3) == '(' &&
		   FirstElement.at(index5) == '+' && 
		   FirstElement.at(index5+1) == '(' &&
		   FirstElement.at(index6) == '*' &&
		   FirstElement.at(index7) == ')' && 
		   FirstElement.at(index7+1) == ')' && 
		   FirstElement.at(index7+2) == ')' && 
		   FirstElement.at(index7+3) == ')' && 
		   FirstElement.at(index7+4) == '-'){

			result = (ConcatenatedStringToFloat*((ConcatenatedStringToFloat2+(ConcatenatedStringToFloat3*ConcatenatedStringToFloat4))/(ConcatenatedStringToFloat5+(ConcatenatedStringToFloat6*ConcatenatedStringToFloat7))))-ConcatenatedStringToFloat8;

        	}
		else if(FirstElement.at(index7+4) == '+'){

                        result = (ConcatenatedStringToFloat*((ConcatenatedStringToFloat2+(ConcatenatedStringToFloat3*ConcatenatedStringToFloat4))/(ConcatenatedStringToFloat5+(ConcatenatedStringToFloat6*ConcatenatedStringToFloat7))))+ConcatenatedStringToFloat8;

                }
        	else if(FirstElement.at(index) == '+' &&
			FirstElement.at(index+1) == '(' &&
			FirstElement.at(index+2) == '-' &&
			FirstElement.at(index+3) == '('){
	

			//eighth
                	for(int i = index8+5; i < FirstElement.length(); i++){

                        	if(FirstElement.at(i) != ')'&&
                           	   FirstElement.at(i) != '(' &&
                           	   FirstElement.at(i) != '*' &&
                           	   FirstElement.at(i) != '/' &&
                           	   FirstElement.at(i) != '+' &&
                           	   FirstElement.at(i) != '-'){VecForConcString9.push_back(FirstElement.at(i));}
                        	else{index9 = i; break;}
                	}


                	for(int i = 0; i < VecForConcString9.size(); i++){
                        	if(i == 0){ConcatenatedString9 = VecForConcString9.at(i);}
                        	else{ConcatenatedString9 += VecForConcString9.at(i);}

                	}
                	
			float ConcatenatedStringToFloat9 = stof(ConcatenatedString9);


			//ninth
                        for(int i = index9+4; i < FirstElement.length(); i++){

                                if(FirstElement.at(i) != ')'&&
                                   FirstElement.at(i) != '(' &&
                                   FirstElement.at(i) != '*' &&
                                   FirstElement.at(i) != '/' &&
                                   FirstElement.at(i) != '+' &&
                                   FirstElement.at(i) != '-'){VecForConcString10.push_back(FirstElement.at(i));}
                                else{index10 = i; break;}
                        }


                        for(int i = 0; i < VecForConcString10.size(); i++){
                                if(i == 0){ConcatenatedString10 = VecForConcString10.at(i);}
                                else{ConcatenatedString10 += VecForConcString10.at(i);}

                        }

                        float ConcatenatedStringToFloat10 = stof(ConcatenatedString10);

			//tenth
                        for(int i = index10+9; i < FirstElement.length(); i++){
          
                                if(FirstElement.at(i) != ')'&&
                                   FirstElement.at(i) != '(' &&
                                   FirstElement.at(i) != '*' &&
                                   FirstElement.at(i) != '/' &&
                                   FirstElement.at(i) != '+' &&
                                   FirstElement.at(i) != '-'){VecForConcString11.push_back(FirstElement.at(i));}
                                else{index11 = i; break;}
                        } 
    

                        for(int i = 0; i < VecForConcString11.size(); i++){
                                if(i == 0){ConcatenatedString11 = VecForConcString11.at(i);}
                                else{ConcatenatedString11 += VecForConcString11.at(i);}
    
                        }
			
                        float ConcatenatedStringToFloat11 = stof(ConcatenatedString11);


		

			result = (ConcatenatedStringToFloat+(-(ConcatenatedStringToFloat2*(log(ConcatenatedStringToFloat3+ConcatenatedStringToFloat4)*(log(ConcatenatedStringToFloat5+ConcatenatedStringToFloat6)*(ConcatenatedStringToFloat7-(-(ConcatenatedStringToFloat8*log(ConcatenatedStringToFloat9+ConcatenatedStringToFloat10)))))))))+ConcatenatedStringToFloat11;
			

		}
		else{cout << "ERROR" << endl;}

		ResultVector.push_back(result);

	}	

	VecForConcString.clear();
	VecForConcString2.clear();
	VecForConcString3.clear();
	VecForConcString4.clear();
	VecForConcString5.clear();
	VecForConcString6.clear();
	VecForConcString7.clear();
	VecForConcString8.clear();
	VecForConcString9.clear();
	VecForConcString10.clear();
	VecForConcString11.clear();
	VecForConcString12.clear();
	VecForConcString13.clear();


}//end of for loop


return ResultVector;

}


