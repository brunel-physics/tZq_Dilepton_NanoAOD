#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;


vector<float> OutputFloats{};


auto CutFlowReader2(const int& LineSpecified, const string& process) { 

  float One, Two, Three;

  string FileName = "CutFlowReport_" + process + "_Blinded.txt";
  ifstream file;
  file.open(FileName.c_str());

  int line_number = 0; 

  if (file.good())
  {
    string str = "";

 //   int line_number = 0;
	
	 while(getline(file, str) && line_number != LineSpecified){
		++line_number;
	}

 }

	if(line_number == LineSpecified){
			
 			if(LineSpecified == 0 || LineSpecified == 1){
			
				string Col1, Col2, Col3; 
				float Col4, Col5, Col6;
 				file >> Col1;
                		file >> Col2;
                		file >> Col3;
                		file >> Col4;
    				file >> Col5;
				file >> Col6;

				One = Col4;
				Two = Col5;
				Three = Col6;

			}
			else if(LineSpecified == 2 || 
				LineSpecified == 3 || 
				LineSpecified == 5 || 
				LineSpecified == 7 ||
				LineSpecified == 8 ||
				LineSpecified == 9 ||
				LineSpecified == 10){
				
				string Col1, Col2, Col3, Col4, Col5;
                                float Col6, Col7, Col8;
                                file >> Col1;
                                file >> Col2;
                                file >> Col3;
                                file >> Col4;
                                file >> Col5;
                                file >> Col6;
				file >> Col7;
				file >> Col8;
				
				One = Col6;
                                Two = Col7;
                                Three = Col8;


			}
			else if(LineSpecified == 4 || LineSpecified == 6){

                                string Col1, Col2, Col3, Col4;
				float Col5, Col6, Col7; 

                                file >> Col1;
                                file >> Col2;
                                file >> Col3;
                                file >> Col4;
                                file >> Col5;
                                file >> Col6;
				file >> Col7;


				One = Col5;
                                Two = Col6;
                                Three = Col7;

                        }


	}

//  }
 

  file.close();
  OutputFloats.push_back(One);
  OutputFloats.push_back(Two);
  OutputFloats.push_back(Three);

  cout << "OutputFloats.size() = " << OutputFloats.size() << endl;

  return OutputFloats;


} 






auto CutFlowReader(){


//lines:
//1 = lepton (ee), 2 = lepton (mumu), 3 = Z mass (ee), 4 = Z mass (mumu), 
//5 = jet (ee), 6 = bjet (ee), 7 = jet (mumu), 8 = bjet (mumu), 9 = W mass (ee), 10 = W mass (mumu)

// string process = "SingleTop_tchannel_top";
 string process = "tZq";

 float One = CutFlowReader2(0, process).at(1);
 float Two = CutFlowReader2(1, process).at(4);
 float Three = CutFlowReader2(2, process).at(7);
 float Four = CutFlowReader2(3, process).at(10);
 float Five = CutFlowReader2(4, process).at(13);
 float Six = CutFlowReader2(5, process).at(16);
 float Seven = CutFlowReader2(6, process).at(19);
 float Eight = CutFlowReader2(7, process).at(22);
 float Nine = CutFlowReader2(8, process).at(25);
 float Ten = CutFlowReader2(9, process).at(28);

 cout << "One is = " << One << endl;
 cout << "Two is = " << Two << endl;
 cout << "Three is = " << Three << endl;
 cout << "Four is = " << Four << endl;
 cout << "Five is = " << Five << endl;
 cout << "Six is = " << Six << endl;
 cout << "Seven is = " << Seven << endl;
 cout << "Eight is = " << Eight << endl;
 cout << "Nine is = " << Nine << endl;
 cout << "Ten is = " << Ten << endl;



}


