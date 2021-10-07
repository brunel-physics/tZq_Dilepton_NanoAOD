#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;



auto CutFlowReader(const int& LineSpecified, const string& process, const string& channel, const string& year) { 

  vector<float> OutputFloats{};

  float One, Two, Three;

  string FileName;

  if(process == "Data"){FileName = "CutFlowReport_Data_Combined_Nominal_" + channel + "__SR_SBR___" + year + ".txt";}
  else{FileName = "CutFlowReport_" + process + "_Nominal_" + channel + "__SR____" + year + ".txt";}

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
			
 			if(LineSpecified == 3 || LineSpecified == 4){
			
				string Col1, Col2, Col3; 
				float Col4, Col5, Col6;
 				file >> Col1;
                		file >> Col2;
                		file >> Col3;
                		file >> Col4;
    				file >> Col5;
				file >> Col6;


				if(process == "tZq"){
					std::cout << "LineSpecified = " << LineSpecified << std::endl;
					std::cout << "Col1 = " << Col1 << std::endl;	
					std::cout << "Col2 = " << Col2 << std::endl;
					std::cout << "Col3 = " << Col3 << std::endl;
					std::cout << "Col4 = " << Col4 << std::endl;
					std::cout << "Col5 = " << Col5 << std::endl;
					std::cout << "Col6 = " << Col6 << std::endl;
				}

				One = Col4;
				Two = Col5;
				Three = Col6;

			}
			else if(LineSpecified == 0 || LineSpecified == 1 || LineSpecified == 2){ 
				
				string Col1, Col2;
                                float Col3, Col4, Col5;
                                file >> Col1;
                                file >> Col2;
                                file >> Col3;
                                file >> Col4;
                                file >> Col5;
			
				if(process == "tZq"){
					std::cout << "LineSpecified = " << LineSpecified << std::endl;
                                	std::cout << "Col1 = " << Col1 << std::endl;
                                	std::cout << "Col2 = " << Col2 << std::endl;
                                	std::cout << "Col3 = " << Col3 << std::endl;
                                	std::cout << "Col4 = " << Col4 << std::endl;
                                	std::cout << "Col5 = " << Col5 << std::endl;
				}
	
				One = Col3;
                                Two = Col4;
                                Three = Col5;


			}
			else{throw std::logic_error("ERROR: Check the line specified.");} 


	}

  }
 

  file.close();
  OutputFloats.push_back(One);
  OutputFloats.push_back(Two);
  OutputFloats.push_back(Three);


  return OutputFloats;


} 








auto CutFlowFunction(const string& ch, const string& y){

   const int num = 5;
   const char *Cuts[num] = {"Lepton cut", "Z mass cut", "Jet cut", "b jet cut", "W mass cut"};

   //tZq started with 13276146 event
   const float tZq[num] = {CutFlowReader(0, "tZq", ch, y).at(1),
			    CutFlowReader(1, "tZq", ch, y).at(1),
			    CutFlowReader(2, "tZq", ch, y).at(1), 
			    CutFlowReader(3, "tZq", ch, y).at(1),	
			    CutFlowReader(4, "tZq", ch, y).at(1)};


   //Z+jets M10To50
   const float ZPlusJets_M10To50[num] = {CutFlowReader(0, "ZPlusJets_M10To50", ch, y).at(1),
                            		  CutFlowReader(1, "ZPlusJets_M10To50", ch, y).at(1),
                            		  CutFlowReader(2, "ZPlusJets_M10To50", ch, y).at(1),         
                            		  CutFlowReader(3, "ZPlusJets_M10To50", ch, y).at(1),         
                            		  CutFlowReader(4, "ZPlusJets_M10To50", ch, y).at(1)};


   //Z+jets M50
   const float ZPlusJets_M50[num] = {CutFlowReader(0, "ZPlusJets_M50_1", ch, y).at(1),
                                      CutFlowReader(1, "ZPlusJets_M50_1", ch, y).at(1),
                                      CutFlowReader(2, "ZPlusJets_M50_1", ch, y).at(1),
                                      CutFlowReader(3, "ZPlusJets_M50_1", ch, y).at(1),
                                      CutFlowReader(4, "ZPlusJets_M50_1", ch, y).at(1) };


   //Z+jets M50 (ext)
   const float ZPlusJets_M50_ext[num] = {CutFlowReader(0, "ZPlusJets_M50_ext", ch, y).at(1),
                                          CutFlowReader(1, "ZPlusJets_M50_ext", ch, y).at(1),
                                          CutFlowReader(2, "ZPlusJets_M50_ext", ch, y).at(1),   
                                          CutFlowReader(3, "ZPlusJets_M50_ext", ch, y).at(1),   
                                          CutFlowReader(4, "ZPlusJets_M50_ext", ch, y).at(1) };


   //ttbar (2l2nu)
   const float ttbar_2l2nu[num] = {CutFlowReader(0, "ttbar_2l2nu", ch, y).at(1),
                                    CutFlowReader(1, "ttbar_2l2nu", ch, y).at(1),
                                    CutFlowReader(2, "ttbar_2l2nu", ch, y).at(1),   
                                    CutFlowReader(3, "ttbar_2l2nu", ch, y).at(1),   
                                    CutFlowReader(4, "ttbar_2l2nu", ch, y).at(1) };


   //ttbar (madgraph)
   const float ttbar_madgraph[num] = {CutFlowReader(0, "ttbar_madgraph", ch, y).at(1),
                                       CutFlowReader(1, "ttbar_madgraph", ch, y).at(1),
                                       CutFlowReader(2, "ttbar_madgraph", ch, y).at(1),   
                                       CutFlowReader(3, "ttbar_madgraph", ch, y).at(1),   
                                       CutFlowReader(4, "ttbar_madgraph", ch, y).at(1) };


   //ttbar (TTToHadronic)
   const float ttbar_TTToHadronic[num] = {CutFlowReader(0, "ttbar_TTToHadronic", ch, y).at(1),
                                          CutFlowReader(1, "ttbar_TTToHadronic", ch, y).at(1),
                                          CutFlowReader(2, "ttbar_TTToHadronic", ch, y).at(1),   
                                          CutFlowReader(3, "ttbar_TTToHadronic", ch, y).at(1),   
                                          CutFlowReader(4, "ttbar_TTToHadronic", ch, y).at(1) };


   //ttbar (TTToSemileptonic)
   const float ttbar_TTToSemileptonic[num] = {CutFlowReader(0, "ttbar_TTToSemileptonic", ch, y).at(1),
                                               CutFlowReader(1, "ttbar_TTToSemileptonic", ch, y).at(1),
                                               CutFlowReader(2, "ttbar_TTToSemileptonic", ch, y).at(1),   
                                               CutFlowReader(3, "ttbar_TTToSemileptonic", ch, y).at(1),   
                                               CutFlowReader(4, "ttbar_TTToSemileptonic", ch, y).at(1) };

   //ttbar (aMCatNLO)
   const float ttbar_aMCatNLO[num] = {CutFlowReader(0, "ttbar_aMCatNLO", ch, y).at(1),
                                       CutFlowReader(1, "ttbar_aMCatNLO", ch, y).at(1),
                                       CutFlowReader(2, "ttbar_aMCatNLO", ch, y).at(1),   
                                       CutFlowReader(3, "ttbar_aMCatNLO", ch, y).at(1),   
                                       CutFlowReader(4, "ttbar_aMCatNLO", ch, y).at(1) };

   
   //single top (s-channel)
   const float singletop_schannel[num] = {CutFlowReader(0, "Singletop_schannel", ch, y).at(1),
                                          CutFlowReader(1, "SingleTop_schannel", ch, y).at(1),
                                          CutFlowReader(2, "SingleTop_schannel", ch, y).at(1),   
                                          CutFlowReader(3, "SingleTop_schannel", ch, y).at(1),   
                                          CutFlowReader(4, "SingleTop_schannel", ch, y).at(1) };



   //single top (t-channel top)
   const float singletop_tchanneltop[num] = {CutFlowReader(0, "SingleTop_tchannel_top", ch, y).at(1),
                                              CutFlowReader(1, "SingleTop_tchannel_top", ch, y).at(1),
                                              CutFlowReader(2, "SingleTop_tchannel_top", ch, y).at(1),   
                                              CutFlowReader(3, "SingleTop_tchannel_top", ch, y).at(1),   
                                              CutFlowReader(4, "SingleTop_tchannel_top", ch, y).at(1) };


   //single top (t-channel antitop)
   const float singletop_tchanneltbar[num] = {CutFlowReader(0, "SingleTop_tchannel_tbar", ch, y).at(1),
                                               CutFlowReader(1, "SingleTop_tchannel_tbar", ch, y).at(1),
                                               CutFlowReader(2, "SingleTop_tchannel_tbar", ch, y).at(1),   
                                               CutFlowReader(3, "SingleTop_tchannel_tbar", ch, y).at(1),   
                                               CutFlowReader(4, "SingleTop_tchannel_tbar", ch, y).at(1) };



   //single top (tHq) 
   const float singletop_tHq[num] = {CutFlowReader(0, "SingleTop_tHq", ch, y).at(1),
                                      CutFlowReader(1, "SingleTop_tHq", ch, y).at(1),
                                      CutFlowReader(2, "SingleTop_tHq", ch, y).at(1),   
                                      CutFlowReader(3, "SingleTop_tHq", ch, y).at(1),   
                                      CutFlowReader(4, "SingleTop_tHq", ch, y).at(1) };


   //single top (tW) 
   const float singletop_tW[num] = {CutFlowReader(0, "SingleTop_tW", ch, y).at(1),
                                     CutFlowReader(1, "SingleTop_tW", ch, y).at(1),
                                     CutFlowReader(2, "SingleTop_tW", ch, y).at(1),   
                                     CutFlowReader(3, "SingleTop_tW", ch, y).at(1),   
                                     CutFlowReader(4, "SingleTop_tW", ch, y).at(1) };



   //single top (tbarW) 
   const float singletop_tbarW[num] = {CutFlowReader(0, "SingleTop_tbarW", ch, y).at(1),
                                        CutFlowReader(1, "SingleTop_tbarW", ch, y).at(1),
                                        CutFlowReader(2, "SingleTop_tbarW", ch, y).at(1),   
                                        CutFlowReader(3, "SingleTop_tbarW", ch, y).at(1),   
                                        CutFlowReader(4, "SingleTop_tbarW", ch, y).at(1) };


   //single top (tZq_W_lept_Z_had) 
   const float singletop_tZq_W_lept_Z_had[num] = {CutFlowReader(0, "SingleTop_tZq_W_lept_Z_had", ch, y).at(1),
                                          	   CutFlowReader(1, "SingleTop_tZq_W_lept_Z_had", ch, y).at(1),
                                          	   CutFlowReader(2, "SingleTop_tZq_W_lept_Z_had", ch, y).at(1),   
                                          	   CutFlowReader(3, "SingleTop_tZq_W_lept_Z_had", ch, y).at(1),   
                                          	   CutFlowReader(4, "SingleTop_tZq_W_lept_Z_had", ch, y).at(1) };


   //single top (tWZ_tWll)
   const float singletop_tWZ_tWll[num] = {CutFlowReader(0, "SingleTop_tWZ_tWll", ch, y).at(1),
                                          CutFlowReader(1, "SingleTop_tWZ_tWll", ch, y).at(1),
                                          CutFlowReader(2, "SingleTop_tWZ_tWll", ch, y).at(1),   
                                          CutFlowReader(3, "SingleTop_tWZ_tWll", ch, y).at(1),   
                                          CutFlowReader(4, "SingleTop_tWZ_tWll", ch, y).at(1) };


   //Diboson (ZZTo2Q2Nu)
   const float VV_ZZTo2Q2Nu[num] = {CutFlowReader(0, "VV_ZZTo2Q2Nu", ch, y).at(1),
                                     CutFlowReader(1, "VV_ZZTo2Q2Nu", ch, y).at(1),
                                     CutFlowReader(2, "VV_ZZTo2Q2Nu", ch, y).at(1),   
                                     CutFlowReader(3, "VV_ZZTo2Q2Nu", ch, y).at(1),   
                                     CutFlowReader(4, "VV_ZZTo2Q2Nu", ch, y).at(1) };


   //Diboson (ZZTo2L2Nu)
   const float VV_ZZTo2L2Nu[num] = {CutFlowReader(0, "VV_ZZTo2L2Nu", ch, y).at(1),
                                     CutFlowReader(1, "VV_ZZTo2L2Nu", ch, y).at(1),
                                     CutFlowReader(2, "VV_ZZTo2L2Nu", ch, y).at(1),   
                                     CutFlowReader(3, "VV_ZZTo2L2Nu", ch, y).at(1),   
                                     CutFlowReader(4, "VV_ZZTo2L2Nu", ch, y).at(1) };


   //Diboson (ZZTo2L2Q)
   const float VV_ZZTo2L2Q[num] = {CutFlowReader(0, "VV_ZZTo2L2Q", ch, y).at(1),
                                    CutFlowReader(1, "VV_ZZTo2L2Q", ch, y).at(1),
                                    CutFlowReader(2, "VV_ZZTo2L2Q", ch, y).at(1),   
                                    CutFlowReader(3, "VV_ZZTo2L2Q", ch, y).at(1),   
                                    CutFlowReader(4, "VV_ZZTo2L2Q", ch, y).at(1) };


   //Diboson (ZZTo4L)
   const float VV_ZZTo4L[num] = {CutFlowReader(0, "VV_ZZTo4L", ch, y).at(1),
                                  CutFlowReader(1, "VV_ZZTo4L", ch, y).at(1),
                                  CutFlowReader(2, "VV_ZZTo4L", ch, y).at(1),   
                                  CutFlowReader(3, "VV_ZZTo4L", ch, y).at(1),   
                                  CutFlowReader(4, "VV_ZZTo4L", ch, y).at(1) };


   //Diboson (WZTo1L1Nu2Q)
   const float VV_WZTo1L1Nu2Q[num] = {CutFlowReader(0, "VV_WZTo1L1Nu2Q", ch, y).at(1),
                                       CutFlowReader(1, "VV_WZTo1L1Nu2Q", ch, y).at(1),
                                       CutFlowReader(2, "VV_WZTo1L1Nu2Q", ch, y).at(1),   
                                       CutFlowReader(3, "VV_WZTo1L1Nu2Q", ch, y).at(1),   
                                       CutFlowReader(4, "VV_WZTo1L1Nu2Q", ch, y).at(1) };


   //Diboson (WZTo2L2Q)
   const float VV_WZTo2L2Q[num] = {CutFlowReader(0, "VV_WZTo2L2Q", ch, y).at(1),
                                    CutFlowReader(1, "VV_WZTo2L2Q", ch, y).at(1),
                                    CutFlowReader(2, "VV_WZTo2L2Q", ch, y).at(1),   
                                    CutFlowReader(3, "VV_WZTo2L2Q", ch, y).at(1),   
                                    CutFlowReader(4, "VV_WZTo2L2Q", ch, y).at(1)};


   //Diboson (WZTo3LNu)
   const float VV_WZTo3LNu[num] = {CutFlowReader(0, "VV_WZTo3LNu", ch, y).at(1),
                                    CutFlowReader(1, "VV_WZTo3LNu", ch, y).at(1),
                                    CutFlowReader(2, "VV_WZTo3LNu", ch, y).at(1),   
                                    CutFlowReader(3, "VV_WZTo3LNu", ch, y).at(1),   
                                    CutFlowReader(4, "VV_WZTo3LNu", ch, y).at(1)};


   //Diboson (WWTo1L1Nu2Q)
   const float VV_WWTo1L1Nu2Q[num] = {CutFlowReader(0, "VV_WWTo1L1Nu2Q", ch, y).at(1),
                                       CutFlowReader(1, "VV_WWTo1L1Nu2Q", ch, y).at(1),
                                       CutFlowReader(2, "VV_WWTo1L1Nu2Q", ch, y).at(1),   
                                       CutFlowReader(3, "VV_WWTo1L1Nu2Q", ch, y).at(1),   
                                       CutFlowReader(4, "VV_WWTo1L1Nu2Q", ch, y).at(1) };


   //Diboson (WWTo2L2Nu)
   const float VV_WWTo2L2Nu[num] = {CutFlowReader(0, "VV_WWTo2L2Nu", ch, y).at(1),
                                     CutFlowReader(1, "VV_WWTo2L2Nu", ch, y).at(1),
                                     CutFlowReader(2, "VV_WWTo2L2Nu", ch, y).at(1),   
                                     CutFlowReader(3, "VV_WWTo2L2Nu", ch, y).at(1),   
                                     CutFlowReader(4, "VV_WWTo2L2Nu", ch, y).at(1)};


   //Diboson (WWToLNuQQ)
   const float VV_WWToLNuQQ[num] = {CutFlowReader(0, "VV_WWToLNuQQ", ch, y).at(1),
                                     CutFlowReader(1, "VV_WWToLNuQQ", ch, y).at(1),
                                     CutFlowReader(2, "VV_WWToLNuQQ", ch, y).at(1),   
                                     CutFlowReader(3, "VV_WWToLNuQQ", ch, y).at(1),   
                                     CutFlowReader(4, "VV_WWToLNuQQ", ch, y).at(1)};

   
   //Diboson (WGToLNuG)
   const float VV_WGToLNuG[num] = {CutFlowReader(0, "VV_WGToLNuG", ch, y).at(1),
                                    CutFlowReader(1, "VV_WGToLNuG", ch, y).at(1),
                                    CutFlowReader(2, "VV_WGToLNuG", ch, y).at(1),   
                                    CutFlowReader(3, "VV_WGToLNuG", ch, y).at(1),   
                                    CutFlowReader(4, "VV_WGToLNuG", ch, y).at(1) };


   //Diboson (ZGToLLG)
   const float VV_ZGToLLG[num] = {CutFlowReader(0, "VV_ZGToLLG", ch, y).at(1),
                                   CutFlowReader(1, "VV_ZGToLLG", ch, y).at(1),
                                   CutFlowReader(2, "VV_ZGToLLG", ch, y).at(1),   
                                   CutFlowReader(3, "VV_ZGToLLG", ch, y).at(1),   
                                   CutFlowReader(4, "VV_ZGToLLG", ch, y).at(1) };


   //Triboson (WWWTo4F)
   const float VVV_WWWTo4F[num] = {CutFlowReader(0, "VVV_WWWTo4F", ch, y).at(1),
                                    CutFlowReader(1, "VVV_WWWTo4F", ch, y).at(1),
                                    CutFlowReader(2, "VVV_WWWTo4F", ch, y).at(1),   
                                    CutFlowReader(3, "VVV_WWWTo4F", ch, y).at(1),   
                                    CutFlowReader(4, "VVV_WWWTo4F", ch, y).at(1) };


   //Triboson (WWZTo4F)
   const float VVV_WWZTo4F[num] = {CutFlowReader(0, "VVV_WWZTo4F", ch, y).at(1),
                                    CutFlowReader(1, "VVV_WWZTo4F", ch, y).at(1),
                                    CutFlowReader(2, "VVV_WWZTo4F", ch, y).at(1),   
                                    CutFlowReader(3, "VVV_WWZTo4F", ch, y).at(1),   
                                    CutFlowReader(4, "VVV_WWZTo4F", ch, y).at(1) };


   //Triboson (WZZTo4F)
   const float VVV_WZZTo4F[num] = {CutFlowReader(0, "VVV_WZZTo4F", ch, y).at(1),
                                    CutFlowReader(1, "VVV_WZZTo4F", ch, y).at(1),
                                    CutFlowReader(2, "VVV_WZZTo4F", ch, y).at(1),   
                                    CutFlowReader(3, "VVV_WZZTo4F", ch, y).at(1),   
                                    CutFlowReader(4, "VVV_WZZTo4F", ch, y).at(1) };


   //Triboson (ZZZTo4F)
   const float VVV_ZZZTo4F[num] = {CutFlowReader(0, "VVV_ZZZTo4F", ch, y).at(1),
                                    CutFlowReader(1, "VVV_ZZZTo4F", ch, y).at(1),
                                    CutFlowReader(2, "VVV_ZZZTo4F", ch, y).at(1),   
                                    CutFlowReader(3, "VVV_ZZZTo4F", ch, y).at(1),   
                                    CutFlowReader(4, "VVV_ZZZTo4F", ch, y).at(1)};


   //w+jets
   const float WPlusJets_WJetsToLNu[num] = {CutFlowReader(0, "WPlusJets_WJetsToLNu", ch, y).at(1),
                                             CutFlowReader(1, "WPlusJets_WJetsToLNu", ch, y).at(1),
                                             CutFlowReader(2, "WPlusJets_WJetsToLNu", ch, y).at(1),   
                                             CutFlowReader(3, "WPlusJets_WJetsToLNu", ch, y).at(1),   
                                             CutFlowReader(4, "WPlusJets_WJetsToLNu", ch, y).at(1) };


   //ttbarV (ttWJetsToLNu) 
   const float ttbarV_ttWJetsToLNu[num] = {CutFlowReader(0, "ttbarV_ttWJetsToLNu", ch, y).at(1),
                                            CutFlowReader(1, "ttbarV_ttWJetsToLNu", ch, y).at(1),
                                            CutFlowReader(2, "ttbarV_ttWJetsToLNu", ch, y).at(1),   
                                            CutFlowReader(3, "ttbarV_ttWJetsToLNu", ch, y).at(1),   
                                            CutFlowReader(4, "ttbarV_ttWJetsToLNu", ch, y).at(1) };


   //ttbarV (ttWJetsToQQ) 
   const float ttbarV_ttWJetsToQQ[num] = {CutFlowReader(0, "ttbarV_ttWJetsToQQ", ch, y).at(1),
                                          CutFlowReader(1, "ttbarV_ttWJetsToQQ", ch, y).at(1),
                                          CutFlowReader(2, "ttbarV_ttWJetsToQQ", ch, y).at(1),   
                                          CutFlowReader(3, "ttbarV_ttWJetsToQQ", ch, y).at(1),   
                                          CutFlowReader(4, "ttbarV_ttWJetsToQQ", ch, y).at(1)};


   //ttbarV (ttgamma) 
   const float ttbarV_ttgamma[num] = {CutFlowReader(0, "ttbarV_ttgamma", ch, y).at(1),
                                       CutFlowReader(1, "ttbarV_ttgamma", ch, y).at(1),
                                       CutFlowReader(2, "ttbarV_ttgamma", ch, y).at(1),   
                                       CutFlowReader(3, "ttbarV_ttgamma", ch, y).at(1),   
                                       CutFlowReader(4, "ttbarV_ttgamma", ch, y).at(1)};


   //ttbarV (ttZToLL) 
   const float ttbarV_ttZToLL[num] = {CutFlowReader(0, "ttbarV_ttZToLL", ch, y).at(1),
                                       CutFlowReader(1, "ttbarV_ttZToLL", ch, y).at(1),
                                       CutFlowReader(2, "ttbarV_ttZToLL", ch, y).at(1),   
                                       CutFlowReader(3, "ttbarV_ttZToLL", ch, y).at(1),   
                                       CutFlowReader(4, "ttbarV_ttZToLL", ch, y).at(1)};


   //ttbarV (ttHTobb) 
   const float ttbarV_ttHTobb[num] = {CutFlowReader(0, "ttbarV_ttHTobb", ch, y).at(1),
                                       CutFlowReader(1, "ttbarV_ttHTobb", ch, y).at(1),
                                       CutFlowReader(2, "ttbarV_ttHTobb", ch, y).at(1),   
                                       CutFlowReader(3, "ttbarV_ttHTobb", ch, y).at(1),   
                                       CutFlowReader(4, "ttbarV_ttHTobb", ch, y).at(1) };


   //ttbarV (ttHToNonbb) 
   const float ttbarV_ttHToNonbb[num] = {CutFlowReader(0, "ttbarV_ttHToNonbb", ch, y).at(1),
                                          CutFlowReader(1, "ttbarV_ttHToNonbb", ch, y).at(1),
                                          CutFlowReader(2, "ttbarV_ttHToNonbb", ch, y).at(1),   
                                          CutFlowReader(3, "ttbarV_ttHToNonbb", ch, y).at(1),   
                                          CutFlowReader(4, "ttbarV_ttHToNonbb", ch, y).at(1)};


   //ttbarV (ttZToLL) 
   const float ttbarV_ttZToLLNuNu[num] = {CutFlowReader(0, "ttbarV_ttZToLLNuNu", ch, y).at(1),
                                           CutFlowReader(1, "ttbarV_ttZToLLNuNu", ch, y).at(1),
                                           CutFlowReader(2, "ttbarV_ttZToLLNuNu", ch, y).at(1),   
                                           CutFlowReader(3, "ttbarV_ttZToLLNuNu", ch, y).at(1),   
                                           CutFlowReader(4, "ttbarV_ttZToLLNuNu", ch, y).at(1)};


   //ttbarV (ttZToQQ) 
   const float ttbarV_ttZToQQ[num] = {CutFlowReader(0, "ttbarV_ttZToQQ_1", ch, y).at(1),
                                       CutFlowReader(1, "ttbarV_ttZToQQ_1", ch, y).at(1),
                                       CutFlowReader(2, "ttbarV_ttZToQQ_1", ch, y).at(1),   
                                       CutFlowReader(3, "ttbarV_ttZToQQ_1", ch, y).at(1),   
                                       CutFlowReader(4, "ttbarV_ttZToQQ_1", ch, y).at(1) };



   //ttbarV (ttZToQQ ext) 
   const float ttbarV_ttZToQQ_ext[num] = {CutFlowReader(0, "ttbarV_ttZToQQ_ext", ch, y).at(1),
                                           CutFlowReader(1, "ttbarV_ttZToQQ_ext", ch, y).at(1),
                                           CutFlowReader(2, "ttbarV_ttZToQQ_ext", ch, y).at(1),   
                                           CutFlowReader(3, "ttbarV_ttZToQQ_ext", ch, y).at(1),   
                                           CutFlowReader(4, "ttbarV_ttZToQQ_ext", ch, y).at(1)};


   //Data
   const float data[num] = {CutFlowReader(0, "Data", ch, y).at(1),
                               CutFlowReader(1, "Data", ch, y).at(1),
                               CutFlowReader(2, "Data", ch, y).at(1),   
                               CutFlowReader(3, "Data", ch, y).at(1),   
                               CutFlowReader(4, "Data", ch, y).at(1)};



   //Initialising arrays for totals
   float ZPlusJets[num] = {};
   float ttbar[num] = {};
   float singletop[num] = {};
   float VV[num] = {};
   float VVV[num] = {};
   float ttbarV[num] = {};
   float Data[num] = {};

   for(int i = 0; i < num; i++){

  	//Arrays for each type of process (arrays for processes of the same type are added)
  	//Z+jets
  	ZPlusJets[i] = ZPlusJets_M10To50[i] + ZPlusJets_M50[i] + ZPlusJets_M50_ext[i];

  	//ttbar
  	ttbar[i] = ttbar_2l2nu[i] + ttbar_madgraph[i] + ttbar_TTToHadronic[i] + ttbar_TTToSemileptonic[i] + ttbar_aMCatNLO[i];

  	//single top
  	singletop[i] = singletop_schannel[i] + singletop_tchanneltop[i] + singletop_tchanneltbar[i] + singletop_tHq[i] + singletop_tW[i] + singletop_tbarW[i] + singletop_tZq_W_lept_Z_had[i] + singletop_tWZ_tWll[i];

  	//diboson
  	VV[i] = VV_ZZTo2Q2Nu[i] + VV_ZZTo2L2Nu[i] + VV_ZZTo2L2Q[i] + VV_ZZTo4L[i] + VV_WZTo1L1Nu2Q[i] + VV_WZTo2L2Q[i] + VV_WZTo3LNu[i] + VV_WWTo1L1Nu2Q[i] + VV_WWTo2L2Nu[i] + VV_WWToLNuQQ[i] + VV_WGToLNuG[i] + VV_ZGToLLG[i];


  	//triboson  
  	VVV[i] = VVV_WWWTo4F[i] + VVV_WWZTo4F[i] + VVV_WZZTo4F[i] + VVV_ZZZTo4F[i];

  	//ttbarV
  	ttbarV[i] = ttbarV_ttWJetsToLNu[i] + ttbarV_ttWJetsToQQ[i] + ttbarV_ttgamma[i] + ttbarV_ttZToLL[i] + ttbarV_ttHTobb[i] + ttbarV_ttHToNonbb[i] + ttbarV_ttZToLLNuNu[i] + ttbarV_ttZToLLNuNu[i] + ttbarV_ttZToQQ[i] + ttbarV_ttZToQQ_ext[i];

        //data 
        Data[i] = data[i];

  }


   //For the canvas
   int W = 800;
   int H = 600;
   int H_ref = 600;
   int W_ref = 800;
   
   float T = 0.08 * H_ref;
   float B = 0.12 * H_ref;
   float L = 0.12 * W_ref;
   float R = 0.04 * W_ref;

   TCanvas *c1 = new TCanvas("c1","cut flow canvas", 50, 50, W, H);
   c1->SetGrid();
   c1->cd();
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetFrameFillStyle(0);
   c1->SetFrameBorderMode(0);
   c1->SetLeftMargin(L / W);
   c1->SetRightMargin(R / W);
   c1->SetTopMargin(T / H);
   c1->SetBottomMargin(B / H);
   c1->SetTickx(0);
   c1->SetTicky(0); 

   TH1F *h_tZq = new TH1F("h_tZq", "tZq", num, 0, num);
   TH1F *h_ZPlusJets = new TH1F("h_ZPlusJets", "Z+jets", num, 0, num);
   TH1F *h_ttbar = new TH1F("h_ttbar", "t#bar{t}", num, 0, num);
   TH1F *h_singletop = new TH1F("h_singletop", "single top", num, 0, num);
   TH1F *h_VV = new TH1F("h_VV", "VV", num, 0, num);
   TH1F *h_VVV = new TH1F("h_VVV", "VVV", num, 0, num);
   TH1F *h_WPlusJets = new TH1F("h_WPlusJets", "W+jets", num, 0, num);
   TH1F *h_ttbarV = new TH1F("h_ttbarV", "t#bar{t}V", num, 0, num);
   TH1F *h_data = new TH1F("h_data", "data", num, 0, num);

   for(int i = 0; i < num; i++){

	h_tZq->Fill(i, tZq[i]);
	h_ZPlusJets->Fill(i, ZPlusJets[i]);
	h_ttbar->Fill(i, ttbar[i]);
	h_singletop->Fill(i, singletop[i]);
	h_VV->Fill(i, VV[i]);
	h_VVV->Fill(i, VVV[i]);
        h_WPlusJets->Fill(i, WPlusJets_WJetsToLNu[i]);
	h_ttbarV->Fill(i, ttbarV[i]);
	h_data->Fill(i, Data[i]);

   }

  

   for (int i=1; i <= num; i++){

	h_tZq->GetXaxis()->SetBinLabel(i,Cuts[i-1]);
	h_ZPlusJets->GetXaxis()->SetBinLabel(i,Cuts[i-1]);
        h_ttbar->GetXaxis()->SetBinLabel(i,Cuts[i-1]);
        h_singletop->GetXaxis()->SetBinLabel(i,Cuts[i-1]);
        h_VV->GetXaxis()->SetBinLabel(i,Cuts[i-1]);
        h_VVV->GetXaxis()->SetBinLabel(i,Cuts[i-1]);
        h_WPlusJets->GetXaxis()->SetBinLabel(i,Cuts[i-1]);
        h_ttbarV->GetXaxis()->SetBinLabel(i,Cuts[i-1]);
        h_data->GetXaxis()->SetBinLabel(i,Cuts[i-1]);

   }


   TPad *pad = new TPad("pad","pad",0.01, 0.315, 0.99, 0.99);
   pad->SetTopMargin(0);
   pad->SetFillColor(0);
   pad->SetBorderMode(0);
   pad->SetFrameFillStyle(0);
   pad->SetFrameBorderMode(0);
   pad->SetLeftMargin(L / W);
   pad->SetRightMargin(R / W);
   pad->SetTopMargin(T / H);
   pad->SetBottomMargin(B / H * 0.3);
   pad->SetTickx(0);
   pad->SetTicky(0);
   pad->Draw();
   pad->cd(); 
  
 

   //VV (ee)
   h_VV->GetXaxis()->SetTitle("Cut");
   h_VV->GetXaxis()->SetTitleOffset(1.5);
   h_VV->GetYaxis()->SetTitle("Events");
   h_VV->GetYaxis()->SetTitleOffset(1.2);
   h_VV->GetYaxis()->SetTitleSize(0.05);
   h_VV->SetFillColor(602); //kOrange+2
   h_VV->SetLineColor(602);
   h_VV->SetStats(0);
   h_VV->Draw("HIST");

   //ZPlusJets (ee)
   h_ZPlusJets->SetFillColor(632); //kRed
   h_ZPlusJets->SetLineColor(632);
   h_ZPlusJets->SetStats(0);
   h_ZPlusJets->Draw("HIST SAME"); 

   //single top (ee)
   h_singletop->SetFillColor(618); //kMagenta+2
   h_singletop->SetLineColor(618);
   h_singletop->SetStats(0);
   h_singletop->Draw("HIST SAME"); 
 
   //tZq (ee)
   h_tZq->SetFillColor(600); //kBlue
   h_tZq->SetLineColor(600);
   h_tZq->SetStats(0);
   h_tZq->Draw("HIST SAME"); 

   //ttbar (ee)
   h_ttbar->SetFillColor(418); //kGreen+2
   h_ttbar->SetLineColor(418);
   h_ttbar->SetStats(0);
   h_ttbar->Draw("HIST SAME"); 

   //ttbarV (ee)
   h_ttbarV->SetFillColor(920); //kGray 
   h_ttbarV->SetLineColor(920);
   h_ttbarV->SetStats(0);
   h_ttbarV->Draw("HIST SAME"); 

   //VVV (ee)
   h_VVV->SetFillColor(400); //kYellow
   h_VVV->SetLineColor(400);
   h_VVV->SetStats(0);
   h_VVV->Draw("HIST SAME");

   //W+jets (ee)
   h_WPlusJets->SetFillColor(425); //kCyan-7 
   h_WPlusJets->SetLineColor(425);
   h_WPlusJets->SetStats(0);
   h_WPlusJets->Draw("HIST SAME");


   //data (ee) 
   h_data->SetMarkerColor(1); //kBlack
   h_data->SetMarkerStyle(20);
   h_data->SetMarkerSize(1);
   h_data->SetStats(0);
   h_data->Draw("P SAME");

   //Legend
   //TLegend * leg = new TLegend(0.65, 0.6, 0.89, 0.89);
   TLegend * leg = new TLegend(0.9, 0.5, 0.8, 0.9);
   leg->SetFillColor(0);
   leg->SetBorderSize(0);
   leg->AddEntry(h_tZq, "tZq", "f");
   leg->AddEntry(h_ZPlusJets, "Z+jets", "f");
   leg->AddEntry(h_ttbar, "t#bar{t}", "f");
   leg->AddEntry(h_ttbarV, "t#bar{t}V", "f");
   leg->AddEntry(h_singletop, "ST", "f");
   leg->AddEntry(h_WPlusJets, "W+jets", "f");
   leg->AddEntry(h_VVV, "VVV", "f");
   leg->AddEntry(h_VV, "VV", "f");
   leg->AddEntry(h_data, "Data", "p");
   leg->Draw();

   // CMS text
   //TPaveText * CMS = new TPaveText(0.1, 0.88, 0.2, 0.98, "blNDC");
   TPaveText * CMS = new TPaveText(0.1, 0.92, 0.2, 0.99, "blNDC");
   TText * CMS_text = CMS->AddText("CMS");
   CMS_text->SetTextFont(61);
   CMS->SetFillStyle(0);
   CMS->SetBorderSize(0);
   CMS->Draw();

   // Work in progress text
   TPaveText * WIP = new TPaveText(0.18, 0.91, 0.41, 0.99, "blNDC");
   TText * WIP_text = WIP->AddText("Work in progress");
   WIP_text->SetTextFont(52);
   WIP->SetFillStyle(0);
   WIP->SetBorderSize(0);
   WIP->Draw();

   // Lumi text
   string lumi_str;

   switch(stoi(y)){
       case 2016: lumi_str = "35.9 fb^{-1} (13 TeV)"; break;
       case 2017: lumi_str = "41.9 fb^{-1} (13 TeV)"; break;
       case 2018: lumi_str = "59.7 fb^{-1} (13 TeV)"; break;
    }

    TPaveText * lumi = new TPaveText(0.65, 0.92, 0.99, 0.99, "blNDC");
    TText * lumi_text = lumi->AddText(lumi_str.c_str());
    lumi->SetFillStyle(0);
    lumi->SetBorderSize(0);
    lumi->Draw();

    // Channel text
    TPaveText * chan = new TPaveText(0.5, 0.92, 0.67, 0.99, "blNDC");

    if(ch == "ee"){chan->AddText("ee channel");}
    else if(ch == "mumu"){chan->AddText("#mu#mu channel");}
    else{chan->AddText("e#mu channel");}

    chan->SetFillStyle(0);
    chan->SetBorderSize(0);
    chan->Draw();
 

   //Remove title
   gStyle->SetOptTitle(0);

   h_VV->GetXaxis()->SetLabelOffset(999);
   h_VV->GetXaxis()->SetLabelSize(0);

   c1->cd();

   TPad *pad2 = new TPad("pad2", "pad2", 0.01, 0.01, 0.99, 0.3275);
   pad2->SetTopMargin(0);
   pad2->SetFillColor(0);
   pad2->SetBorderMode(0);
   pad2->SetFrameFillStyle(0);
   pad2->SetFrameBorderMode(0);
   pad2->SetLeftMargin(L / W);
   pad2->SetRightMargin(R / W);
   pad2->SetTopMargin(T / H);
   pad2->SetBottomMargin(B / H * 2.1);
   pad2->SetTickx(0);
   pad2->SetTicky(0);
   pad2->SetGridy(1);
   pad2->Draw();
   pad2->cd(); 

   //for the ratio plot
   TH1F *h_MC = new TH1F("h_MC", "MC", num, 0, num);
  
   h_MC->Add(h_tZq);
   h_MC->Add(h_ZPlusJets);
   h_MC->Add(h_ttbar);
   h_MC->Add(h_singletop);
   h_MC->Add(h_VV);
   h_MC->Add(h_VVV);
   h_MC->Add(h_WPlusJets);
   h_MC->Add(h_ttbarV);


   TH1D * rp = (TH1D*)(h_data->Clone());
   rp->Divide(h_MC);
   rp->SetStats(false);
   rp->GetXaxis()->SetNdivisions(6, 5, 0);
   rp->GetYaxis()->SetNdivisions(6, 5, 0);
   rp->GetXaxis()->SetLabelSize(0.15);
   rp->GetXaxis()->SetLabelOffset(0.030);
   rp->GetXaxis()->SetTitleSize(0.11);
   rp->GetYaxis()->SetLabelSize(0.08);
   rp->GetYaxis()->SetTitle("Data/MC");
   rp->GetYaxis()->SetTitleOffset(0.6);
   rp->GetYaxis()->SetTitleSize(0.1);
   rp->GetYaxis()->CenterTitle();
   rp->SetMinimum(0.5);
   rp->SetMaximum(1.5);
   rp->Draw(); 


   //Saving the histograms in pdf files
   std::string output = "CutFlow_" + ch + "_" + y + ".pdf";
   h_VV->SaveAs(output.c_str());

}


auto CutFlow(){

  CutFlowFunction("ee", "2016");

}
