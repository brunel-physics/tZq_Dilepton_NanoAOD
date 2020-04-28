#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;



auto CutFlowReader(const int& LineSpecified, const string& process) { 

  vector<float> OutputFloats{};

  float One, Two, Three;

  string FileName = "CutFlowReport_" + process + "_Blinded.txt";
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

  }
 

  file.close();
  OutputFloats.push_back(One);
  OutputFloats.push_back(Two);
  OutputFloats.push_back(Three);


  return OutputFloats;


} 








auto CutFlow(){

   const int num = 5;
   const char *Cuts[num] = {"lepton", "z mass", "jets", "bjets", "w mass"};

   cout << "CutFlowReader(0, tZq).at(1) = " << CutFlowReader(0, "WPlusJets").at(1) << endl;
   cout << "CutFlowReader(1, tZq).at(1) = " << CutFlowReader(1, "WPlusJets").at(1) << endl;
   cout << "CutFlowReader(2, tZq).at(1) = " << CutFlowReader(2, "WPlusJets").at(1) << endl;
   cout << "CutFlowReader(3, tZq).at(1) = " << CutFlowReader(3, "WPlusJets").at(1) << endl;
   cout << "CutFlowReader(4, tZq).at(1) = " << CutFlowReader(4, "WPlusJets").at(1) << endl;
   cout << "CutFlowReader(5, tZq).at(1) = " << CutFlowReader(5, "WPlusJets").at(1) << endl;
   cout << "CutFlowReader(6, tZq).at(1) = " << CutFlowReader(6, "WPlusJets").at(1) << endl;
   cout << "CutFlowReader(7, tZq).at(1) = " << CutFlowReader(7, "WPlusJets").at(1) << endl;
   cout << "CutFlowReader(8, tZq).at(1) = " << CutFlowReader(8, "WPlusJets").at(1) << endl;
   cout << "CutFlowReader(9, tZq).at(1) = " << CutFlowReader(9, "WPlusJets").at(1) << endl;

   //tZq started with 13276146 event
   const float tZq_ee[num] = {CutFlowReader(0, "tZq").at(1),
			    CutFlowReader(2, "tZq").at(1),
			    CutFlowReader(4, "tZq").at(1), 
			    CutFlowReader(5, "tZq").at(1),	
			    CutFlowReader(8, "tZq").at(1)};


   const float tZq_mumu[num] = {CutFlowReader(1, "tZq").at(1),
			     CutFlowReader(3, "tZq").at(1),
			     CutFlowReader(6, "tZq").at(1),
			     CutFlowReader(7, "tZq").at(1),
			     CutFlowReader(9, "tZq").at(1)};
 


   //Z+jets M10To50
   const float ZPlusJets_M10To50_ee[num] = {CutFlowReader(0, "ZPlusJets_M10To50").at(1),
                            		  CutFlowReader(2, "ZPlusJets_M10To50").at(1),
                            		  CutFlowReader(4, "ZPlusJets_M10To50").at(1),         
                            		  CutFlowReader(5, "ZPlusJets_M10To50").at(1),         
                            		  CutFlowReader(8, "ZPlusJets_M10To50").at(1)};

   const float ZPlusJets_M10To50_mumu[num] = {CutFlowReader(1, "ZPlusJets_M10To50").at(1),
                             	            CutFlowReader(3, "ZPlusJets_M10To50").at(1),
                             		    CutFlowReader(6, "ZPlusJets_M10To50").at(1),
                             		    CutFlowReader(7, "ZPlusJets_M10To50").at(1),
                             		    CutFlowReader(9, "ZPlusJets_M10To50").at(1)};

   //Z+jets M50
   const float ZPlusJets_M50_ee[num] = {CutFlowReader(0, "ZPlusJets_M50_1").at(1),
                                      CutFlowReader(2, "ZPlusJets_M50_1").at(1),
                                      CutFlowReader(4, "ZPlusJets_M50_1").at(1),
                                      CutFlowReader(5, "ZPlusJets_M50_1").at(1),
                                      CutFlowReader(8, "ZPlusJets_M50_1").at(1) };

   const float ZPlusJets_M50_mumu[num] = {CutFlowReader(1, "ZPlusJets_M50_1").at(1),
                                        CutFlowReader(3, "ZPlusJets_M50_1").at(1),
                                        CutFlowReader(6, "ZPlusJets_M50_1").at(1),
                                        CutFlowReader(7, "ZPlusJets_M50_1").at(1),
                                        CutFlowReader(9, "ZPlusJets_M50_1").at(1) };

   //Z+jets M50 (ext)
   const float ZPlusJets_M50_ext_ee[num] = {CutFlowReader(0, "ZPlusJets_M50_ext").at(1),
                                          CutFlowReader(2, "ZPlusJets_M50_ext").at(1),
                                          CutFlowReader(4, "ZPlusJets_M50_ext").at(1),   
                                          CutFlowReader(5, "ZPlusJets_M50_ext").at(1),   
                                          CutFlowReader(8, "ZPlusJets_M50_ext").at(1) };

   const float ZPlusJets_M50_ext_mumu[num] = {CutFlowReader(1, "ZPlusJets_M50_ext").at(1),
                                            CutFlowReader(3, "ZPlusJets_M50_ext").at(1),
                                            CutFlowReader(6, "ZPlusJets_M50_ext").at(1),
                                            CutFlowReader(7, "ZPlusJets_M50_ext").at(1),
                                            CutFlowReader(9, "ZPlusJets_M50_ext").at(1) };

   //ttbar (2l2nu)
   const float ttbar_2l2nu_ee[num] = {CutFlowReader(0, "ttbar_2l2nu").at(1),
                                    CutFlowReader(2, "ttbar_2l2nu").at(1),
                                    CutFlowReader(4, "ttbar_2l2nu").at(1),   
                                    CutFlowReader(5, "ttbar_2l2nu").at(1),   
                                    CutFlowReader(8, "ttbar_2l2nu").at(1) };

   const float ttbar_2l2nu_mumu[num] = {CutFlowReader(1, "ttbar_2l2nu").at(1),
                                      CutFlowReader(3, "ttbar_2l2nu").at(1),
                                      CutFlowReader(6, "ttbar_2l2nu").at(1),
                                      CutFlowReader(7, "ttbar_2l2nu").at(1),
                                      CutFlowReader(9, "ttbar_2l2nu").at(1) };

   //ttbar (madgraph)
   const float ttbar_madgraph_ee[num] = {CutFlowReader(0, "ttbar_madgraph").at(1),
                                       CutFlowReader(2, "ttbar_madgraph").at(1),
                                       CutFlowReader(4, "ttbar_madgraph").at(1),   
                                       CutFlowReader(5, "ttbar_madgraph").at(1),   
                                       CutFlowReader(8, "ttbar_madgraph").at(1) };

   const float ttbar_madgraph_mumu[num] = {CutFlowReader(1, "ttbar_madgraph").at(1),
                                         CutFlowReader(3, "ttbar_madgraph").at(1),
                                         CutFlowReader(6, "ttbar_madgraph").at(1),
                                         CutFlowReader(7, "ttbar_madgraph").at(1),
                                         CutFlowReader(9, "ttbar_madgraph").at(1) }; 

   //ttbar (TTToHadronic)
   const float ttbar_TTToHadronic_ee[num] = {CutFlowReader(0, "ttbar_TTToHadronic").at(1),
                                          CutFlowReader(2, "ttbar_TTToHadronic").at(1),
                                          CutFlowReader(4, "ttbar_TTToHadronic").at(1),   
                                          CutFlowReader(5, "ttbar_TTToHadronic").at(1),   
                                          CutFlowReader(8, "ttbar_TTToHadronic").at(1) };

   const float ttbar_TTToHadronic_mumu[num] = {CutFlowReader(1, "ttbar_TTToHadronic").at(1),
                                             CutFlowReader(3, "ttbar_TTToHadronic").at(1),
                                             CutFlowReader(6, "ttbar_TTToHadronic").at(1),
                                             CutFlowReader(7, "ttbar_TTToHadronic").at(1),
                                             CutFlowReader(9, "ttbar_TTToHadronic").at(1) };

   //ttbar (TTToSemileptonic)
   const float ttbar_TTToSemileptonic_ee[num] = {CutFlowReader(0, "ttbar_TTToSemileptonic").at(1),
                                               CutFlowReader(2, "ttbar_TTToSemileptonic").at(1),
                                               CutFlowReader(4, "ttbar_TTToSemileptonic").at(1),   
                                               CutFlowReader(5, "ttbar_TTToSemileptonic").at(1),   
                                               CutFlowReader(8, "ttbar_TTToSemileptonic").at(1) };

   const float ttbar_TTToSemileptonic_mumu[num] = {CutFlowReader(1, "ttbar_TTToSemileptonic").at(1),
                                      	         CutFlowReader(3, "ttbar_TTToSemileptonic").at(1),
                                      	         CutFlowReader(6, "ttbar_TTToSemileptonic").at(1),
                                      	         CutFlowReader(7, "ttbar_TTToSemileptonic").at(1),
                                      		 CutFlowReader(9, "ttbar_TTToSemileptonic").at(1) };

   //ttbar (aMCatNLO)
   const float ttbar_aMCatNLO_ee[num] = {CutFlowReader(0, "ttbar_aMCatNLO").at(1),
                                       CutFlowReader(2, "ttbar_aMCatNLO").at(1),
                                       CutFlowReader(4, "ttbar_aMCatNLO").at(1),   
                                       CutFlowReader(5, "ttbar_aMCatNLO").at(1),   
                                       CutFlowReader(8, "ttbar_aMCatNLO").at(1) };

   const float ttbar_aMCatNLO_mumu[num] = {CutFlowReader(1, "ttbar_aMCatNLO").at(1),
                                         CutFlowReader(3, "ttbar_aMCatNLO").at(1),
                                         CutFlowReader(6, "ttbar_aMCatNLO").at(1),
                                         CutFlowReader(7, "ttbar_aMCatNLO").at(1),
                                         CutFlowReader(9, "ttbar_aMCatNLO").at(1) };
   
   //single top (s-channel)
   const float singletop_schannel_ee[num] = {CutFlowReader(0, "Singletop_schannel").at(1),
                                          CutFlowReader(2, "SingleTop_schannel").at(1),
                                          CutFlowReader(4, "SingleTop_schannel").at(1),   
                                          CutFlowReader(5, "SingleTop_schannel").at(1),   
                                          CutFlowReader(8, "SingleTop_schannel").at(1) };

   const float singletop_schannel_mumu[num] = {CutFlowReader(1, "SingleTop_schannel").at(1),
                                             CutFlowReader(3, "SingleTop_schannel").at(1),
                                             CutFlowReader(6, "SingleTop_schannel").at(1),
                                             CutFlowReader(7, "SingleTop_schannel").at(1),
                                             CutFlowReader(9, "SingleTop_schannel").at(1) };


   //single top (t-channel top)
   const float singletop_tchanneltop_ee[num] = {CutFlowReader(0, "SingleTop_tchannel_top").at(1),
                                              CutFlowReader(2, "SingleTop_tchannel_top").at(1),
                                              CutFlowReader(4, "SingleTop_tchannel_top").at(1),   
                                              CutFlowReader(5, "SingleTop_tchannel_top").at(1),   
                                              CutFlowReader(8, "SingleTop_tchannel_top").at(1) };

   const float singletop_tchanneltop_mumu[num] = {CutFlowReader(1, "SingleTop_tchannel_top").at(1),
                                      		CutFlowReader(3, "SingleTop_tchannel_top").at(1),
                                      		CutFlowReader(6, "SingleTop_tchannel_top").at(1),
                                      		CutFlowReader(7, "SingleTop_tchannel_top").at(1),
                                      		CutFlowReader(9, "SingleTop_tchannel_top").at(1) };

   //single top (t-channel antitop)
   const float singletop_tchanneltbar_ee[num] = {CutFlowReader(0, "SingleTop_tchannel_tbar").at(1),
                                               CutFlowReader(2, "SingleTop_tchannel_tbar").at(1),
                                               CutFlowReader(4, "SingleTop_tchannel_tbar").at(1),   
                                               CutFlowReader(5, "SingleTop_tchannel_tbar").at(1),   
                                               CutFlowReader(8, "SingleTop_tchannel_tbar").at(1) };

   const float singletop_tchanneltbar_mumu[num] = {CutFlowReader(1, "SingleTop_tchannel_tbar").at(1),
                                      		 CutFlowReader(3, "SingleTop_tchannel_tbar").at(1),
                                      		 CutFlowReader(6, "SingleTop_tchannel_tbar").at(1),
                                      		 CutFlowReader(7, "SingleTop_tchannel_tbar").at(1),
                                      		 CutFlowReader(9, "SingleTop_tchannel_tbar").at(1) };


   //single top (tHq) 
   const float singletop_tHq_ee[num] = {CutFlowReader(0, "SingleTop_tHq").at(1),
                                      CutFlowReader(2, "SingleTop_tHq").at(1),
                                      CutFlowReader(4, "SingleTop_tHq").at(1),   
                                      CutFlowReader(5, "SingleTop_tHq").at(1),   
                                      CutFlowReader(8, "SingleTop_tHq").at(1) };

   const float singletop_tHq_mumu[num] = {CutFlowReader(1, "SingleTop_tHq").at(1),
                                        CutFlowReader(3, "SingleTop_tHq").at(1),
                                        CutFlowReader(6, "SingleTop_tHq").at(1),
                                        CutFlowReader(7, "SingleTop_tHq").at(1),
                                        CutFlowReader(9, "SingleTop_tHq").at(1) };

   //single top (tW) 
   const float singletop_tW_ee[num] = {CutFlowReader(0, "SingleTop_tW").at(1),
                                     CutFlowReader(2, "SingleTop_tW").at(1),
                                     CutFlowReader(4, "SingleTop_tW").at(1),   
                                     CutFlowReader(5, "SingleTop_tW").at(1),   
                                     CutFlowReader(8, "SingleTop_tW").at(1) };

   const float singletop_tW_mumu[num] = {CutFlowReader(1, "SingleTop_tW").at(1),
                                       CutFlowReader(3, "SingleTop_tW").at(1),
                                       CutFlowReader(6, "SingleTop_tW").at(1),
                                       CutFlowReader(7, "SingleTop_tW").at(1),
                                       CutFlowReader(9, "SingleTop_tW").at(1) };



   //single top (tbarW) 
   const float singletop_tbarW_ee[num] = {CutFlowReader(0, "SingleTop_tbarW").at(1),
                                        CutFlowReader(2, "SingleTop_tbarW").at(1),
                                        CutFlowReader(4, "SingleTop_tbarW").at(1),   
                                        CutFlowReader(5, "SingleTop_tbarW").at(1),   
                                        CutFlowReader(8, "SingleTop_tbarW").at(1) };

   const float singletop_tbarW_mumu[num] = {CutFlowReader(1, "SingleTop_tbarW").at(1),
                                          CutFlowReader(3, "SingleTop_tbarW").at(1),
                                          CutFlowReader(6, "SingleTop_tbarW").at(1),
                                          CutFlowReader(7, "SingleTop_tbarW").at(1),
                                          CutFlowReader(9, "SingleTop_tbarW").at(1) };

   //single top (tZq_W_lept_Z_had) 
   const float singletop_tZq_W_lept_Z_had_ee[num] = {CutFlowReader(0, "SingleTop_tZq_W_lept_Z_had").at(1),
                                          	   CutFlowReader(2, "SingleTop_tZq_W_lept_Z_had").at(1),
                                          	   CutFlowReader(4, "SingleTop_tZq_W_lept_Z_had").at(1),   
                                          	   CutFlowReader(5, "SingleTop_tZq_W_lept_Z_had").at(1),   
                                          	   CutFlowReader(8, "SingleTop_tZq_W_lept_Z_had").at(1) };

   const float singletop_tZq_W_lept_Z_had_mumu[num] = {CutFlowReader(1, "SingleTop_tZq_W_lept_Z_had").at(1),
                                      		     CutFlowReader(3, "SingleTop_tZq_W_lept_Z_had").at(1),
                                      		     CutFlowReader(6, "SingleTop_tZq_W_lept_Z_had").at(1),
                                      		     CutFlowReader(7, "SingleTop_tZq_W_lept_Z_had").at(1),
                                      		     CutFlowReader(9, "SingleTop_tZq_W_lept_Z_had").at(1) };

   //single top (tWZ_tWll)
   const float singletop_tWZ_tWll_ee[num] = {CutFlowReader(0, "SingleTop_tWZ_tWll").at(1),
                                          CutFlowReader(2, "SingleTop_tWZ_tWll").at(1),
                                          CutFlowReader(4, "SingleTop_tWZ_tWll").at(1),   
                                          CutFlowReader(5, "SingleTop_tWZ_tWll").at(1),   
                                          CutFlowReader(8, "SingleTop_tWZ_tWll").at(1) };

   const float singletop_tWZ_tWll_mumu[num] = {CutFlowReader(1, "SingleTop_tWZ_tWll").at(1),
                                             CutFlowReader(3, "SingleTop_tWZ_tWll").at(1),
                                             CutFlowReader(6, "SingleTop_tWZ_tWll").at(1),
                                             CutFlowReader(7, "SingleTop_tWZ_tWll").at(1),
                                             CutFlowReader(9, "SingleTop_tWZ_tWll").at(1) };

   //Diboson (ZZTo2Q2Nu)
   const float VV_ZZTo2Q2Nu_ee[num] = {CutFlowReader(0, "VV_ZZTo2Q2Nu").at(1),
                                     CutFlowReader(2, "VV_ZZTo2Q2Nu").at(1),
                                     CutFlowReader(4, "VV_ZZTo2Q2Nu").at(1),   
                                     CutFlowReader(5, "VV_ZZTo2Q2Nu").at(1),   
                                     CutFlowReader(8, "VV_ZZTo2Q2Nu").at(1) };

   const float VV_ZZTo2Q2Nu_mumu[num] = {CutFlowReader(1, "VV_ZZTo2Q2Nu").at(1),
                                      CutFlowReader(3, "VV_ZZTo2Q2Nu").at(1),
                                      CutFlowReader(6, "VV_ZZTo2Q2Nu").at(1),
                                      CutFlowReader(7, "VV_ZZTo2Q2Nu").at(1),
                                      CutFlowReader(9, "VV_ZZTo2Q2Nu").at(1) };

   //Diboson (ZZTo2L2Nu)
   const float VV_ZZTo2L2Nu_ee[num] = {CutFlowReader(0, "VV_ZZTo2L2Nu").at(1),
                                     CutFlowReader(2, "VV_ZZTo2L2Nu").at(1),
                                     CutFlowReader(4, "VV_ZZTo2L2Nu").at(1),   
                                     CutFlowReader(5, "VV_ZZTo2L2Nu").at(1),   
                                     CutFlowReader(8, "VV_ZZTo2L2Nu").at(1) };

   const float VV_ZZTo2L2Nu_mumu[num] = {CutFlowReader(1, "VV_ZZTo2L2Nu").at(1),
                                      CutFlowReader(3, "VV_ZZTo2L2Nu").at(1),
                                      CutFlowReader(6, "VV_ZZTo2L2Nu").at(1),
                                      CutFlowReader(7, "VV_ZZTo2L2Nu").at(1),
                                      CutFlowReader(9, "VV_ZZTo2L2Nu").at(1) };

   //Diboson (ZZTo2L2Q)
   const float VV_ZZTo2L2Q_ee[num] = {CutFlowReader(0, "VV_ZZTo2L2Q").at(1),
                                    CutFlowReader(2, "VV_ZZTo2L2Q").at(1),
                                    CutFlowReader(4, "VV_ZZTo2L2Q").at(1),   
                                    CutFlowReader(5, "VV_ZZTo2L2Q").at(1),   
                                    CutFlowReader(8, "VV_ZZTo2L2Q").at(1) };

   const float VV_ZZTo2L2Q_mumu[num] = {CutFlowReader(1, "VV_ZZTo2L2Q").at(1),
                                      CutFlowReader(3, "VV_ZZTo2L2Q").at(1),
                                      CutFlowReader(6, "VV_ZZTo2L2Q").at(1),
                                      CutFlowReader(7, "VV_ZZTo2L2Q").at(1),
                                      CutFlowReader(9, "VV_ZZTo2L2Q").at(1) };

   //Diboson (ZZTo4L)
   const float VV_ZZTo4L_ee[num] = {CutFlowReader(0, "VV_ZZTo4L").at(1),
                                  CutFlowReader(2, "VV_ZZTo4L").at(1),
                                  CutFlowReader(4, "VV_ZZTo4L").at(1),   
                                  CutFlowReader(5, "VV_ZZTo4L").at(1),   
                                  CutFlowReader(8, "VV_ZZTo4L").at(1) };

   const float VV_ZZTo4L_mumu[num] = {CutFlowReader(1, "VV_ZZTo4L").at(1),
                                    CutFlowReader(3, "VV_ZZTo4L").at(1),
                                    CutFlowReader(6, "VV_ZZTo4L").at(1),
                                    CutFlowReader(7, "VV_ZZTo4L").at(1),
                                    CutFlowReader(9, "VV_ZZTo4L").at(1)};

   //Diboson (WZTo1L1Nu2Q)
   const float VV_WZTo1L1Nu2Q_ee[num] = {CutFlowReader(0, "VV_WZTo1L1Nu2Q").at(1),
                                       CutFlowReader(2, "VV_WZTo1L1Nu2Q").at(1),
                                       CutFlowReader(4, "VV_WZTo1L1Nu2Q").at(1),   
                                       CutFlowReader(5, "VV_WZTo1L1Nu2Q").at(1),   
                                       CutFlowReader(8, "VV_WZTo1L1Nu2Q").at(1) };

   const float VV_WZTo1L1Nu2Q_mumu[num] = {CutFlowReader(1, "VV_WZTo1L1Nu2Q").at(1),
                                         CutFlowReader(3, "VV_WZTo1L1Nu2Q").at(1),
                                         CutFlowReader(6, "VV_WZTo1L1Nu2Q").at(1),
                                         CutFlowReader(7, "VV_WZTo1L1Nu2Q").at(1),
                                         CutFlowReader(9, "VV_WZTo1L1Nu2Q").at(1)};

   //Diboson (WZTo2L2Q)
   const float VV_WZTo2L2Q_ee[num] = {CutFlowReader(0, "VV_WZTo2L2Q").at(1),
                                    CutFlowReader(2, "VV_WZTo2L2Q").at(1),
                                    CutFlowReader(4, "VV_WZTo2L2Q").at(1),   
                                    CutFlowReader(5, "VV_WZTo2L2Q").at(1),   
                                    CutFlowReader(8, "VV_WZTo2L2Q").at(1)};

   const float VV_WZTo2L2Q_mumu[num] = {CutFlowReader(1, "VV_WZTo2L2Q").at(1),
                                      CutFlowReader(3, "VV_WZTo2L2Q").at(1),
                                      CutFlowReader(6, "VV_WZTo2L2Q").at(1),
                                      CutFlowReader(7, "VV_WZTo2L2Q").at(1),
                                      CutFlowReader(9, "VV_WZTo2L2Q").at(1)};

   //Diboson (WZTo3LNu)
   const float VV_WZTo3LNu_ee[num] = {CutFlowReader(0, "VV_WZTo3LNu").at(1),
                                    CutFlowReader(2, "VV_WZTo3LNu").at(1),
                                    CutFlowReader(4, "VV_WZTo3LNu").at(1),   
                                    CutFlowReader(5, "VV_WZTo3LNu").at(1),   
                                    CutFlowReader(8, "VV_WZTo3LNu").at(1)};

   const float VV_WZTo3LNu_mumu[num] = {CutFlowReader(1, "VV_WZTo3LNu").at(1),
                                      CutFlowReader(3, "VV_WZTo3LNu").at(1),
                                      CutFlowReader(6, "VV_WZTo3LNu").at(1),
                                      CutFlowReader(7, "VV_WZTo3LNu").at(1),
                                      CutFlowReader(9, "VV_WZTo3LNu").at(1)};

   //Diboson (WWTo1L1Nu2Q)
   const float VV_WWTo1L1Nu2Q_ee[num] = {CutFlowReader(0, "VV_WWTo1L1Nu2Q").at(1),
                                       CutFlowReader(2, "VV_WWTo1L1Nu2Q").at(1),
                                       CutFlowReader(4, "VV_WWTo1L1Nu2Q").at(1),   
                                       CutFlowReader(5, "VV_WWTo1L1Nu2Q").at(1),   
                                       CutFlowReader(8, "VV_WWTo1L1Nu2Q").at(1) };

   const float VV_WWTo1L1Nu2Q_mumu[num] = {CutFlowReader(1, "VV_WWTo1L1Nu2Q").at(1),
                                         CutFlowReader(3, "VV_WWTo1L1Nu2Q").at(1),
                                         CutFlowReader(6, "VV_WWTo1L1Nu2Q").at(1),
                                         CutFlowReader(7, "VV_WWTo1L1Nu2Q").at(1),
                                         CutFlowReader(9, "VV_WWTo1L1Nu2Q").at(1)};

   //Diboson (WWTo2L2Nu)
   const float VV_WWTo2L2Nu_ee[num] = {CutFlowReader(0, "VV_WWTo2L2Nu").at(1),
                                     CutFlowReader(2, "VV_WWTo2L2Nu").at(1),
                                     CutFlowReader(4, "VV_WWTo2L2Nu").at(1),   
                                     CutFlowReader(5, "VV_WWTo2L2Nu").at(1),   
                                     CutFlowReader(8, "VV_WWTo2L2Nu").at(1)};

   const float VV_WWTo2L2Nu_mumu[num] = {CutFlowReader(1, "VV_WWTo2L2Nu").at(1),
                                      CutFlowReader(3, "VV_WWTo2L2Nu").at(1),
                                      CutFlowReader(6, "VV_WWTo2L2Nu").at(1),
                                      CutFlowReader(7, "VV_WWTo2L2Nu").at(1),
                                      CutFlowReader(9, "VV_WWTo2L2Nu").at(1)}; 

   //Diboson (WWToLNuQQ)
   const float VV_WWToLNuQQ_ee[num] = {CutFlowReader(0, "VV_WWToLNuQQ").at(1),
                                     CutFlowReader(2, "VV_WWToLNuQQ").at(1),
                                     CutFlowReader(4, "VV_WWToLNuQQ").at(1),   
                                     CutFlowReader(5, "VV_WWToLNuQQ").at(1),   
                                     CutFlowReader(8, "VV_WWToLNuQQ").at(1)};

   const float VV_WWToLNuQQ_mumu[num] = {CutFlowReader(1, "VV_WWToLNuQQ").at(1),
                                      CutFlowReader(3, "VV_WWToLNuQQ").at(1),
                                      CutFlowReader(6, "VV_WWToLNuQQ").at(1),
                                      CutFlowReader(7, "VV_WWToLNuQQ").at(1),
                                      CutFlowReader(9, "VV_WWToLNuQQ").at(1)};
   
   //Diboson (WGToLNuG)
   const float VV_WGToLNuG_ee[num] = {CutFlowReader(0, "VV_WGToLNuG").at(1),
                                    CutFlowReader(2, "VV_WGToLNuG").at(1),
                                    CutFlowReader(4, "VV_WGToLNuG").at(1),   
                                    CutFlowReader(5, "VV_WGToLNuG").at(1),   
                                    CutFlowReader(8, "VV_WGToLNuG").at(1) };

   const float VV_WGToLNuG_mumu[num] = {CutFlowReader(1, "VV_WGToLNuG").at(1),
                                      CutFlowReader(3, "VV_WGToLNuG").at(1),
                                      CutFlowReader(6, "VV_WGToLNuG").at(1),
                                      CutFlowReader(7, "VV_WGToLNuG").at(1),
                                      CutFlowReader(9, "VV_WGToLNuG").at(1)};

   //Diboson (ZGToLLG)
   const float VV_ZGToLLG_ee[num] = {CutFlowReader(0, "VV_ZGToLLG").at(1),
                                   CutFlowReader(2, "VV_ZGToLLG").at(1),
                                   CutFlowReader(4, "VV_ZGToLLG").at(1),   
                                   CutFlowReader(5, "VV_ZGToLLG").at(1),   
                                   CutFlowReader(8, "VV_ZGToLLG").at(1) };

   const float VV_ZGToLLG_mumu[num] = {CutFlowReader(1, "VV_ZGToLLG").at(1),
                                     CutFlowReader(3, "VV_ZGToLLG").at(1),
                                     CutFlowReader(6, "VV_ZGToLLG").at(1),
                                     CutFlowReader(7, "VV_ZGToLLG").at(1),
                                     CutFlowReader(9, "VV_ZGToLLG").at(1)};

   //Triboson (WWWTo4F)
   const float VVV_WWWTo4F_ee[num] = {CutFlowReader(0, "VVV_WWWTo4F").at(1),
                                    CutFlowReader(2, "VVV_WWWTo4F").at(1),
                                    CutFlowReader(4, "VVV_WWWTo4F").at(1),   
                                    CutFlowReader(5, "VVV_WWWTo4F").at(1),   
                                    CutFlowReader(8, "VVV_WWWTo4F").at(1) };

   const float VVV_WWWTo4F_mumu[num] = {CutFlowReader(1, "VVV_WWWTo4F").at(1),
                                      CutFlowReader(3, "VVV_WWWTo4F").at(1),
                                      CutFlowReader(6, "VVV_WWWTo4F").at(1),
                                      CutFlowReader(7, "VVV_WWWTo4F").at(1),
                                      CutFlowReader(9, "VVV_WWWTo4F").at(1)};

   //Triboson (WWZTo4F)
   const float VVV_WWZTo4F_ee[num] = {CutFlowReader(0, "VVV_WWZTo4F").at(1),
                                    CutFlowReader(2, "VVV_WWZTo4F").at(1),
                                    CutFlowReader(4, "VVV_WWZTo4F").at(1),   
                                    CutFlowReader(5, "VVV_WWZTo4F").at(1),   
                                    CutFlowReader(8, "VVV_WWZTo4F").at(1) };

   const float VVV_WWZTo4F_mumu[num] = {CutFlowReader(1, "VVV_WWZTo4F").at(1),
                                      CutFlowReader(3, "VVV_WWZTo4F").at(1),
                                      CutFlowReader(6, "VVV_WWZTo4F").at(1),
                                      CutFlowReader(7, "VVV_WWZTo4F").at(1),
                                      CutFlowReader(9, "VVV_WWZTo4F").at(1)};

   //Triboson (WZZTo4F)
   const float VVV_WZZTo4F_ee[num] = {CutFlowReader(0, "VVV_WZZTo4F").at(1),
                                    CutFlowReader(2, "VVV_WZZTo4F").at(1),
                                    CutFlowReader(4, "VVV_WZZTo4F").at(1),   
                                    CutFlowReader(5, "VVV_WZZTo4F").at(1),   
                                    CutFlowReader(8, "VVV_WZZTo4F").at(1) };

   const float VVV_WZZTo4F_mumu[num] = {CutFlowReader(1, "VVV_WZZTo4F").at(1),
                                      CutFlowReader(3, "VVV_WZZTo4F").at(1),
                                      CutFlowReader(6, "VVV_WZZTo4F").at(1),
                                      CutFlowReader(7, "VVV_WZZTo4F").at(1),
                                      CutFlowReader(9, "VVV_WZZTo4F").at(1)};

   //Triboson (ZZZTo4F)
   const float VVV_ZZZTo4F_ee[num] = {CutFlowReader(0, "VVV_ZZZTo4F").at(1),
                                    CutFlowReader(2, "VVV_ZZZTo4F").at(1),
                                    CutFlowReader(4, "VVV_ZZZTo4F").at(1),   
                                    CutFlowReader(5, "VVV_ZZZTo4F").at(1),   
                                    CutFlowReader(8, "VVV_ZZZTo4F").at(1)};

   const float VVV_ZZZTo4F_mumu[num] = {CutFlowReader(1, "VVV_ZZZTo4F").at(1),
                                      CutFlowReader(3, "VVV_ZZZTo4F").at(1),
                                      CutFlowReader(6, "VVV_ZZZTo4F").at(1),
                                      CutFlowReader(7, "VVV_ZZZTo4F").at(1),
                                      CutFlowReader(9, "VVV_ZZZTo4F").at(1)};

   //w+jets
   const float WPlusJets_WJetsToLNu_ee[num] = {CutFlowReader(0, "WPlusJets_WJetsToLNu").at(1),
                                             CutFlowReader(2, "WPlusJets_WJetsToLNu").at(1),
                                             CutFlowReader(4, "WPlusJets_WJetsToLNu").at(1),   
                                             CutFlowReader(5, "WPlusJets_WJetsToLNu").at(1),   
                                             CutFlowReader(8, "WPlusJets_WJetsToLNu").at(1) };

   const float WPlusJets_WJetsToLNu_mumu[num] = {CutFlowReader(1, "WPlusJets_WJetsToLNu").at(1),
                                               CutFlowReader(3, "WPlusJets_WJetsToLNu").at(1),
                                               CutFlowReader(6, "WPlusJets_WJetsToLNu").at(1),
                                               CutFlowReader(7, "WPlusJets_WJetsToLNu").at(1),
                                               CutFlowReader(9, "WPlusJets_WJetsToLNu").at(1)};

   //ttbarV (ttWJetsToLNu) 
   const float ttbarV_ttWJetsToLNu_ee[num] = {CutFlowReader(0, "ttbarV_ttWJetsToLNu").at(1),
                                            CutFlowReader(2, "ttbarV_ttWJetsToLNu").at(1),
                                            CutFlowReader(4, "ttbarV_ttWJetsToLNu").at(1),   
                                            CutFlowReader(5, "ttbarV_ttWJetsToLNu").at(1),   
                                            CutFlowReader(8, "ttbarV_ttWJetsToLNu").at(1) };

   const float ttbarV_ttWJetsToLNu_mumu[num] = {CutFlowReader(1, "ttbarV_ttWJetsToLNu").at(1),
                                      	      CutFlowReader(3, "ttbarV_ttWJetsToLNu").at(1),
                                      	      CutFlowReader(6, "ttbarV_ttWJetsToLNu").at(1),
                                      	      CutFlowReader(7, "ttbarV_ttWJetsToLNu").at(1),
                                      	      CutFlowReader(9, "ttbarV_ttWJetsToLNu").at(1)};

   //ttbarV (ttWJetsToQQ) 
   const float ttbarV_ttWJetsToQQ_ee[num] = {CutFlowReader(0, "ttbarV_ttWJetsToQQ").at(1),
                                          CutFlowReader(2, "ttbarV_ttWJetsToQQ").at(1),
                                          CutFlowReader(4, "ttbarV_ttWJetsToQQ").at(1),   
                                          CutFlowReader(5, "ttbarV_ttWJetsToQQ").at(1),   
                                          CutFlowReader(8, "ttbarV_ttWJetsToQQ").at(1)};

   const float ttbarV_ttWJetsToQQ_mumu[num] = {CutFlowReader(1, "ttbarV_ttWJetsToQQ").at(1),
                                             CutFlowReader(3, "ttbarV_ttWJetsToQQ").at(1),
                                             CutFlowReader(6, "ttbarV_ttWJetsToQQ").at(1),
                                             CutFlowReader(7, "ttbarV_ttWJetsToQQ").at(1),
                                             CutFlowReader(9, "ttbarV_ttWJetsToQQ").at(1) };

   //ttbarV (ttgamma) 
   const float ttbarV_ttgamma_ee[num] = {CutFlowReader(0, "ttbarV_ttgamma").at(1),
                                       CutFlowReader(2, "ttbarV_ttgamma").at(1),
                                       CutFlowReader(4, "ttbarV_ttgamma").at(1),   
                                       CutFlowReader(5, "ttbarV_ttgamma").at(1),   
                                       CutFlowReader(8, "ttbarV_ttgamma").at(1)};

   const float ttbarV_ttgamma_mumu[num] = {CutFlowReader(1, "ttbarV_ttgamma").at(1),
                                         CutFlowReader(3, "ttbarV_ttgamma").at(1),
                                         CutFlowReader(6, "ttbarV_ttgamma").at(1),
                                         CutFlowReader(7, "ttbarV_ttgamma").at(1),
                                         CutFlowReader(9, "ttbarV_ttgamma").at(1) };

   //ttbarV (ttZToLL) 
   const float ttbarV_ttZToLL_ee[num] = {CutFlowReader(0, "ttbarV_ttZToLL").at(1),
                                       CutFlowReader(2, "ttbarV_ttZToLL").at(1),
                                       CutFlowReader(4, "ttbarV_ttZToLL").at(1),   
                                       CutFlowReader(5, "ttbarV_ttZToLL").at(1),   
                                       CutFlowReader(8, "ttbarV_ttZToLL").at(1)};

   const float ttbarV_ttZToLL_mumu[num] = {CutFlowReader(1, "ttbarV_ttZToLL").at(1),
                                         CutFlowReader(3, "ttbarV_ttZToLL").at(1),
                                         CutFlowReader(6, "ttbarV_ttZToLL").at(1),
                                         CutFlowReader(7, "ttbarV_ttZToLL").at(1),
                                         CutFlowReader(9, "ttbarV_ttZToLL").at(1) };

   //ttbarV (ttHTobb) 
   const float ttbarV_ttHTobb_ee[num] = {CutFlowReader(0, "ttbarV_ttHTobb").at(1),
                                       CutFlowReader(2, "ttbarV_ttHTobb").at(1),
                                       CutFlowReader(4, "ttbarV_ttHTobb").at(1),   
                                       CutFlowReader(5, "ttbarV_ttHTobb").at(1),   
                                       CutFlowReader(8, "ttbarV_ttHTobb").at(1) };

   const float ttbarV_ttHTobb_mumu[num] = {CutFlowReader(1, "ttbarV_ttHTobb").at(1),
                                         CutFlowReader(3, "ttbarV_ttHTobb").at(1),
                                         CutFlowReader(6, "ttbarV_ttHTobb").at(1),
                                         CutFlowReader(7, "ttbarV_ttHTobb").at(1),
                                         CutFlowReader(9, "ttbarV_ttHTobb").at(1) };

   //ttbarV (ttHToNonbb) 
   const float ttbarV_ttHToNonbb_ee[num] = {CutFlowReader(0, "ttbarV_ttHToNonbb").at(1),
                                          CutFlowReader(2, "ttbarV_ttHToNonbb").at(1),
                                          CutFlowReader(4, "ttbarV_ttHToNonbb").at(1),   
                                          CutFlowReader(5, "ttbarV_ttHToNonbb").at(1),   
                                          CutFlowReader(8, "ttbarV_ttHToNonbb").at(1)};

   const float ttbarV_ttHToNonbb_mumu[num] = {CutFlowReader(1, "ttbarV_ttHToNonbb").at(1),
                                            CutFlowReader(3, "ttbarV_ttHToNonbb").at(1),
                                            CutFlowReader(6, "ttbarV_ttHToNonbb").at(1),
                                            CutFlowReader(7, "ttbarV_ttHToNonbb").at(1),
                                            CutFlowReader(9, "ttbarV_ttHToNonbb").at(1) };

   //ttbarV (ttZToLL) 
   const float ttbarV_ttZToLLNuNu_ee[num] = {CutFlowReader(0, "ttbarV_ttZToLLNuNu").at(1),
                                           CutFlowReader(2, "ttbarV_ttZToLLNuNu").at(1),
                                           CutFlowReader(4, "ttbarV_ttZToLLNuNu").at(1),   
                                           CutFlowReader(5, "ttbarV_ttZToLLNuNu").at(1),   
                                           CutFlowReader(8, "ttbarV_ttZToLLNuNu").at(1)};

   const float ttbarV_ttZToLLNuNu_mumu[num] = {CutFlowReader(1, "ttbarV_ttZToLLNuNu").at(1),
                                             CutFlowReader(3, "ttbarV_ttZToLLNuNu").at(1),
                                             CutFlowReader(6, "ttbarV_ttZToLLNuNu").at(1),
                                             CutFlowReader(7, "ttbarV_ttZToLLNuNu").at(1),
                                             CutFlowReader(9, "ttbarV_ttZToLLNuNu").at(1)};


   //ttbarV (ttZToQQ) 
   const float ttbarV_ttZToQQ_ee[num] = {CutFlowReader(0, "ttbarV_ttZToQQ_1").at(1),
                                       CutFlowReader(2, "ttbarV_ttZToQQ_1").at(1),
                                       CutFlowReader(4, "ttbarV_ttZToQQ_1").at(1),   
                                       CutFlowReader(5, "ttbarV_ttZToQQ_1").at(1),   
                                       CutFlowReader(8, "ttbarV_ttZToQQ_1").at(1) };

   const float ttbarV_ttZToQQ_mumu[num] = {CutFlowReader(1, "ttbarV_ttZToQQ_1").at(1),
                                         CutFlowReader(3, "ttbarV_ttZToQQ_1").at(1),
                                         CutFlowReader(6, "ttbarV_ttZToQQ_1").at(1),
                                         CutFlowReader(7, "ttbarV_ttZToQQ_1").at(1),
                                         CutFlowReader(9, "ttbarV_ttZToQQ_1").at(1) };


   //ttbarV (ttZToQQ ext) 
   const float ttbarV_ttZToQQ_ext_ee[num] = {CutFlowReader(0, "ttbarV_ttZToQQ_ext").at(1),
                                           CutFlowReader(2, "ttbarV_ttZToQQ_ext").at(1),
                                           CutFlowReader(4, "ttbarV_ttZToQQ_ext").at(1),   
                                           CutFlowReader(5, "ttbarV_ttZToQQ_ext").at(1),   
                                           CutFlowReader(8, "ttbarV_ttZToQQ_ext").at(1)};

   const float ttbarV_ttZToQQ_ext_mumu[num] = {CutFlowReader(1, "ttbarV_ttZToQQ_ext").at(1),
                                             CutFlowReader(3, "ttbarV_ttZToQQ_ext").at(1),
                                             CutFlowReader(6, "ttbarV_ttZToQQ_ext").at(1),
                                             CutFlowReader(7, "ttbarV_ttZToQQ_ext").at(1),
                                             CutFlowReader(9, "ttbarV_ttZToQQ_ext").at(1)};


   //data (DoubleEG)
   const float data_DoubleEGRunB_ee[num] = {CutFlowReader(0, "data_DoubleEGRunB").at(1),
                                          CutFlowReader(2, "data_DoubleEGRunB").at(1),
                                          CutFlowReader(4, "data_DoubleEGRunB").at(1),   
                                          CutFlowReader(5, "data_DoubleEGRunB").at(1),   
                                          CutFlowReader(8, "data_DoubleEGRunB").at(1)};

   const float data_DoubleEGRunC_ee[num] = {CutFlowReader(0, "data_DoubleEGRunC").at(1),
                                          CutFlowReader(2, "data_DoubleEGRunC").at(1),
                                          CutFlowReader(4, "data_DoubleEGRunC").at(1),   
                                          CutFlowReader(5, "data_DoubleEGRunC").at(1),   
                                          CutFlowReader(8, "data_DoubleEGRunC").at(1) };

   const float data_DoubleEGRunD_ee[num] = {CutFlowReader(0, "data_DoubleEGRunD").at(1),
                                          CutFlowReader(2, "data_DoubleEGRunD").at(1),
                                          CutFlowReader(4, "data_DoubleEGRunD").at(1),   
                                          CutFlowReader(5, "data_DoubleEGRunD").at(1),   
                                          CutFlowReader(8, "data_DoubleEGRunD").at(1) };

   const float data_DoubleEGRunE_ee[num] = {CutFlowReader(0, "data_DoubleEGRunE").at(1),
                                          CutFlowReader(2, "data_DoubleEGRunE").at(1),
                                          CutFlowReader(4, "data_DoubleEGRunE").at(1),   
                                          CutFlowReader(5, "data_DoubleEGRunE").at(1),   
                                          CutFlowReader(8, "data_DoubleEGRunE").at(1) };

   const float data_DoubleEGRunF_ee[num] = {CutFlowReader(0, "data_DoubleEGRunF").at(1),
                                          CutFlowReader(2, "data_DoubleEGRunF").at(1),
                                          CutFlowReader(4, "data_DoubleEGRunF").at(1),   
                                          CutFlowReader(5, "data_DoubleEGRunF").at(1),   
                                          CutFlowReader(8, "data_DoubleEGRunF").at(1) };

   const float data_DoubleEGRunB_mumu[num] = {CutFlowReader(1, "data_DoubleEGRunB").at(1),
                                            CutFlowReader(3, "data_DoubleEGRunB").at(1),
                                            CutFlowReader(6, "data_DoubleEGRunB").at(1),
                                            CutFlowReader(7, "data_DoubleEGRunB").at(1),
                                            CutFlowReader(9, "data_DoubleEGRunB").at(1) };

   const float data_DoubleEGRunC_mumu[num] = {CutFlowReader(1, "data_DoubleEGRunC").at(1),
                                            CutFlowReader(3, "data_DoubleEGRunC").at(1),
                                            CutFlowReader(6, "data_DoubleEGRunC").at(1),
                                            CutFlowReader(7, "data_DoubleEGRunC").at(1),
                                            CutFlowReader(9, "data_DoubleEGRunC").at(1) };

   const float data_DoubleEGRunD_mumu[num] = {CutFlowReader(1, "data_DoubleEGRunD").at(1),
                                            CutFlowReader(3, "data_DoubleEGRunD").at(1),
                                            CutFlowReader(6, "data_DoubleEGRunD").at(1),
                                            CutFlowReader(7, "data_DoubleEGRunD").at(1),
                                            CutFlowReader(9, "data_DoubleEGRunD").at(1) };

   const float data_DoubleEGRunE_mumu[num] = {CutFlowReader(1, "data_DoubleEGRunE").at(1),
                                            CutFlowReader(3, "data_DoubleEGRunE").at(1),
                                            CutFlowReader(6, "data_DoubleEGRunE").at(1),
                                            CutFlowReader(7, "data_DoubleEGRunE").at(1),
                                            CutFlowReader(9, "data_DoubleEGRunE").at(1) };

   const float data_DoubleEGRunF_mumu[num] = {CutFlowReader(1, "data_DoubleEGRunF").at(1),
                                            CutFlowReader(3, "data_DoubleEGRunF").at(1),
                                            CutFlowReader(6, "data_DoubleEGRunF").at(1),
                                            CutFlowReader(7, "data_DoubleEGRunF").at(1),
                                            CutFlowReader(9, "data_DoubleEGRunF").at(1) };


   //data (DoubleMuon)
   const float data_DoubleMuonRunB_ee[num] = {CutFlowReader(0, "data_DoubleMuonRunB").at(1),
                                            CutFlowReader(2, "data_DoubleMuonRunB").at(1),
                                            CutFlowReader(4, "data_DoubleMuonRunB").at(1),   
                                            CutFlowReader(5, "data_DoubleMuonRunB").at(1),   
                                            CutFlowReader(8, "data_DoubleMuonRunB").at(1) };

   const float data_DoubleMuonRunC_ee[num] = {CutFlowReader(0, "data_DoubleMuonRunC").at(1),
                                            CutFlowReader(2, "data_DoubleMuonRunC").at(1),
                                            CutFlowReader(4, "data_DoubleMuonRunC").at(1),   
                                            CutFlowReader(5, "data_DoubleMuonRunC").at(1),   
                                            CutFlowReader(8, "data_DoubleMuonRunC").at(1) };

   const float data_DoubleMuonRunD_ee[num] = {CutFlowReader(0, "data_DoubleMuonRunD").at(1),
                                            CutFlowReader(2, "data_DoubleMuonRunD").at(1),
                                            CutFlowReader(4, "data_DoubleMuonRunD").at(1),   
                                            CutFlowReader(5, "data_DoubleMuonRunD").at(1),   
                                            CutFlowReader(8, "data_DoubleMuonRunD").at(1) };

   const float data_DoubleMuonRunE_ee[num] = {CutFlowReader(0, "data_DoubleMuonRunE").at(1),
                                            CutFlowReader(2, "data_DoubleMuonRunE").at(1),
                                            CutFlowReader(4, "data_DoubleMuonRunE").at(1),   
                                            CutFlowReader(5, "data_DoubleMuonRunE").at(1),   
                                            CutFlowReader(8, "data_DoubleMuonRunE").at(1) };

   const float data_DoubleMuonRunF_ee[num] = {CutFlowReader(0, "data_DoubleMuonRunF").at(1),
                                            CutFlowReader(2, "data_DoubleMuonRunF").at(1),
                                            CutFlowReader(4, "data_DoubleMuonRunF").at(1),   
                                            CutFlowReader(5, "data_DoubleMuonRunF").at(1),   
                                            CutFlowReader(8, "data_DoubleMuonRunF").at(1) };
 

   const float data_DoubleMuonRunB_mumu[num] = {CutFlowReader(1, "data_DoubleMuonRunB").at(1),
                                              CutFlowReader(3, "data_DoubleMuonRunB").at(1),
                                              CutFlowReader(6, "data_DoubleMuonRunB").at(1),
                                              CutFlowReader(7, "data_DoubleMuonRunB").at(1),
                                              CutFlowReader(9, "data_DoubleMuonRunB").at(1) };

   const float data_DoubleMuonRunC_mumu[num] = {CutFlowReader(1, "data_DoubleMuonRunC").at(1),
                                              CutFlowReader(3, "data_DoubleMuonRunC").at(1),
                                              CutFlowReader(6, "data_DoubleMuonRunC").at(1),
                                              CutFlowReader(7, "data_DoubleMuonRunC").at(1),
                                              CutFlowReader(9, "data_DoubleMuonRunC").at(1) };

   const float data_DoubleMuonRunD_mumu[num] = {CutFlowReader(1, "data_DoubleMuonRunD").at(1),
                                              CutFlowReader(3, "data_DoubleMuonRunD").at(1),
                                              CutFlowReader(6, "data_DoubleMuonRunD").at(1),
                                              CutFlowReader(7, "data_DoubleMuonRunD").at(1),
                                              CutFlowReader(9, "data_DoubleMuonRunD").at(1)};

   const float data_DoubleMuonRunE_mumu[num] = {CutFlowReader(1, "data_DoubleMuonRunE").at(1),
                                              CutFlowReader(3, "data_DoubleMuonRunE").at(1),
                                              CutFlowReader(6, "data_DoubleMuonRunE").at(1),
                                              CutFlowReader(7, "data_DoubleMuonRunE").at(1),
                                              CutFlowReader(9, "data_DoubleMuonRunE").at(1)};

   const float data_DoubleMuonRunF_mumu[num] = {CutFlowReader(1, "data_DoubleMuonRunF").at(1),
                                              CutFlowReader(3, "data_DoubleMuonRunF").at(1),
                                              CutFlowReader(6, "data_DoubleMuonRunF").at(1),
                                              CutFlowReader(7, "data_DoubleMuonRunF").at(1),
                                              CutFlowReader(9, "data_DoubleMuonRunF").at(1)};



   //Initialising arrays for totals
   float ZPlusJets_ee[num] = {};
   float ZPlusJets_mumu[num] = {};
   float ttbar_ee[num] = {};
   float ttbar_mumu[num] = {};
   float singletop_ee[num] = {};
   float singletop_mumu[num] = {};
   float VV_ee[num] = {};
   float VV_mumu[num] = {};
   float VVV_ee[num] = {};
   float VVV_mumu[num] = {};
   float ttbarV_ee[num] = {};
   float ttbarV_mumu[num] = {};
   float data_ee[num] = {};
   float data_mumu[num] = {};

   for(int i = 0; i < num; i++){

  	//Arrays for each type of process (arrays for processes of the same type are added)
  	//Z+jets
  	ZPlusJets_ee[i] = ZPlusJets_M10To50_ee[i] + ZPlusJets_M50_ee[i] + ZPlusJets_M50_ext_ee[i];
  	ZPlusJets_mumu[i] = ZPlusJets_M10To50_mumu[i] + ZPlusJets_M50_mumu[i] + ZPlusJets_M50_ext_mumu[i];

  	//ttbar
  	ttbar_ee[i] = ttbar_2l2nu_ee[i] + ttbar_madgraph_ee[i] + ttbar_TTToHadronic_ee[i] + ttbar_TTToSemileptonic_ee[i] + ttbar_aMCatNLO_ee[i];
  	ttbar_mumu[i] = ttbar_2l2nu_mumu[i] + ttbar_madgraph_mumu[i] + ttbar_TTToHadronic_mumu[i] + ttbar_TTToSemileptonic_mumu[i] + ttbar_aMCatNLO_mumu[i];

  	//single top
  	singletop_ee[i] = singletop_schannel_ee[i] + singletop_tchanneltop_ee[i] + singletop_tchanneltbar_ee[i] + singletop_tHq_ee[i] + singletop_tW_ee[i] + singletop_tbarW_ee[i] + singletop_tZq_W_lept_Z_had_ee[i] + singletop_tWZ_tWll_ee[i];

  	singletop_mumu[i] = singletop_schannel_mumu[i] + singletop_tchanneltop_mumu[i] + singletop_tchanneltbar_mumu[i] + singletop_tHq_mumu[i] + singletop_tW_mumu[i] + singletop_tbarW_mumu[i] + singletop_tZq_W_lept_Z_had_mumu[i] + singletop_tWZ_tWll_mumu[i];

  	//diboson
  	VV_ee[i] = VV_ZZTo2Q2Nu_ee[i] + VV_ZZTo2L2Nu_ee[i] + VV_ZZTo2L2Q_ee[i] + VV_ZZTo4L_ee[i] + VV_WZTo1L1Nu2Q_ee[i] + VV_WZTo2L2Q_ee[i] + VV_WZTo3LNu_ee[i] + VV_WWTo1L1Nu2Q_ee[i] + VV_WWTo2L2Nu_ee[i] + VV_WWToLNuQQ_ee[i] + VV_WGToLNuG_ee[i] + VV_ZGToLLG_ee[i];

  	VV_mumu[i] = VV_ZZTo2Q2Nu_mumu[i] + VV_ZZTo2L2Nu_mumu[i] + VV_ZZTo2L2Q_mumu[i] + VV_ZZTo4L_mumu[i] + VV_WZTo1L1Nu2Q_mumu[i] + VV_WZTo2L2Q_mumu[i] + VV_WZTo3LNu_mumu[i] + VV_WWTo1L1Nu2Q_mumu[i] + VV_WWTo2L2Nu_mumu[i] + VV_WWToLNuQQ_mumu[i] + VV_WGToLNuG_mumu[i] + VV_ZGToLLG_mumu[i];

  	//triboson  
  	VVV_ee[i] = VVV_WWWTo4F_ee[i] + VVV_WWZTo4F_ee[i] + VVV_WZZTo4F_ee[i] + VVV_ZZZTo4F_ee[i];
  	VVV_mumu[i] = VVV_WWWTo4F_mumu[i] + VVV_WWZTo4F_mumu[i] + VVV_WZZTo4F_mumu[i] + VVV_ZZZTo4F_mumu[i];

  	//ttbarV
  	ttbarV_ee[i] = ttbarV_ttWJetsToLNu_ee[i] + ttbarV_ttWJetsToQQ_ee[i] + ttbarV_ttgamma_ee[i] + ttbarV_ttZToLL_ee[i] + ttbarV_ttHTobb_ee[i] + ttbarV_ttHToNonbb_ee[i] + ttbarV_ttZToLLNuNu_ee[i] + ttbarV_ttZToLLNuNu_ee[i] + ttbarV_ttZToQQ_ee[i] + ttbarV_ttZToQQ_ext_ee[i];

  	ttbarV_mumu[i] = ttbarV_ttWJetsToLNu_mumu[i] + ttbarV_ttWJetsToQQ_mumu[i] + ttbarV_ttgamma_mumu[i] + ttbarV_ttZToLL_mumu[i] + ttbarV_ttHTobb_mumu[i] + ttbarV_ttHToNonbb_mumu[i] + ttbarV_ttZToLLNuNu_mumu[i] + ttbarV_ttZToLLNuNu_mumu[i] + ttbarV_ttZToQQ_mumu[i] + ttbarV_ttZToQQ_ext_mumu[i];

       //data 
       data_ee[i] = data_DoubleEGRunB_ee[i] + data_DoubleEGRunC_ee[i] + data_DoubleEGRunD_ee[i] + data_DoubleEGRunE_ee[i] + data_DoubleEGRunF_ee[i] + data_DoubleMuonRunB_ee[i] + data_DoubleMuonRunC_ee[i] + data_DoubleMuonRunD_ee[i] + data_DoubleMuonRunE_ee[i] + data_DoubleMuonRunF_ee[i];

       data_mumu[i] = data_DoubleMuonRunB_mumu[i] + data_DoubleMuonRunC_mumu[i] + data_DoubleMuonRunD_mumu[i] + data_DoubleMuonRunE_mumu[i] + data_DoubleMuonRunF_mumu[i] + data_DoubleEGRunB_mumu[i] + data_DoubleEGRunC_mumu[i] + data_DoubleEGRunD_mumu[i] + data_DoubleEGRunE_mumu[i] + data_DoubleEGRunF_mumu[i];
  	

  }


   TCanvas *c1 = new TCanvas("c1","cut flow canvas",10,10,900,500);
   c1->SetGrid();
   c1->SetBottomMargin(0.15);

   TH1F *h_tZq_ee = new TH1F("h_tZq_ee", "tZq", num, 0, num);
   TH1F *h_tZq_mumu = new TH1F("h_tZq_mumu", "tZq", num, 0, num);
   TH1F *h_ZPlusJets_ee = new TH1F("h_ZPlusJets_ee", "Z+jets", num, 0, num);
   TH1F *h_ZPlusJets_mumu = new TH1F("h_ZPlusJets_mumu", "Z+jets", num, 0, num);
   TH1F *h_ttbar_ee = new TH1F("h_ttbar_ee", "t#bar{t}", num, 0, num);
   TH1F *h_ttbar_mumu = new TH1F("h_ttbar_mumu", "t#bar{t}", num, 0, num);	
   TH1F *h_singletop_ee = new TH1F("h_singletop_ee", "single top", num, 0, num);
   TH1F *h_singletop_mumu = new TH1F("h_singletop_mumu", "single top", num, 0, num);
   TH1F *h_VV_ee = new TH1F("h_VV_ee", "VV", num, 0, num);
   TH1F *h_VV_mumu = new TH1F("h_VV_mumu", "VV", num, 0, num);
   TH1F *h_VVV_ee = new TH1F("h_VVV_ee", "VVV", num, 0, num);
   TH1F *h_VVV_mumu = new TH1F("h_VVV_mumu", "VVV", num, 0, num);
   TH1F *h_WPlusJets_ee = new TH1F("h_WPlusJets_ee", "W+jets", num, 0, num);
   TH1F *h_WPlusJets_mumu = new TH1F("h_WPlusJets_mumu", "W+jets", num, 0, num);
   TH1F *h_ttbarV_ee = new TH1F("h_ttbarV_ee", "t#bar{t}V", num, 0, num);
   TH1F *h_ttbarV_mumu = new TH1F("h_ttbarV_mumu", "t#bar{t}V", num, 0, num);
   TH1F *h_data_ee = new TH1F("h_data_ee", "data", num, 0, num);
   TH1F *h_data_mumu = new TH1F("h_data_mumu", "data", num, 0, num);


   for(int i = 0; i < num; i++){

	cout << "tZq_ee[i] = " << tZq_ee[i] << endl;

	h_tZq_ee->Fill(i, tZq_ee[i]);
	h_tZq_mumu->Fill(i, tZq_mumu[i]);
	h_ZPlusJets_ee->Fill(i, ZPlusJets_ee[i]);
	h_ZPlusJets_mumu->Fill(i, ZPlusJets_mumu[i]);
	h_ttbar_ee->Fill(i, ttbar_ee[i]);
	h_ttbar_mumu->Fill(i, ttbar_mumu[i]);
	h_singletop_ee->Fill(i, singletop_ee[i]);
	h_singletop_mumu->Fill(i, singletop_mumu[i]);
	h_VV_ee->Fill(i, VV_ee[i]);
	h_VV_mumu->Fill(i, VV_mumu[i]);
	h_VVV_ee->Fill(i, VVV_ee[i]);
	h_VVV_mumu->Fill(i, VVV_mumu[i]);
        h_WPlusJets_ee->Fill(i, WPlusJets_WJetsToLNu_ee[i]);
	h_WPlusJets_mumu->Fill(i, WPlusJets_WJetsToLNu_mumu[i]);
	h_ttbarV_ee->Fill(i, ttbarV_ee[i]);
	h_ttbarV_mumu->Fill(i, ttbarV_mumu[i]);
	h_data_ee->Fill(i, data_ee[i]);
	h_data_mumu->Fill(i, data_mumu[i]);

	cout << "tZq_ee[i] = " << tZq_ee[i] << endl;

   }


   for (int i=1; i <= num; i++){

	h_tZq_ee->GetXaxis()->SetBinLabel(i,Cuts[i-1]);
	h_tZq_mumu->GetXaxis()->SetBinLabel(i,Cuts[i-1]);
	h_ZPlusJets_ee->GetXaxis()->SetBinLabel(i,Cuts[i-1]);
        h_ZPlusJets_mumu->GetXaxis()->SetBinLabel(i,Cuts[i-1]);
        h_ttbar_ee->GetXaxis()->SetBinLabel(i,Cuts[i-1]);
        h_ttbar_mumu->GetXaxis()->SetBinLabel(i,Cuts[i-1]);
        h_singletop_ee->GetXaxis()->SetBinLabel(i,Cuts[i-1]);
        h_singletop_mumu->GetXaxis()->SetBinLabel(i,Cuts[i-1]);
        h_VV_ee->GetXaxis()->SetBinLabel(i,Cuts[i-1]);
        h_VV_mumu->GetXaxis()->SetBinLabel(i,Cuts[i-1]);
        h_VVV_ee->GetXaxis()->SetBinLabel(i,Cuts[i-1]);
        h_VVV_mumu->GetXaxis()->SetBinLabel(i,Cuts[i-1]);
        h_WPlusJets_ee->GetXaxis()->SetBinLabel(i,Cuts[i-1]);
        h_WPlusJets_mumu->GetXaxis()->SetBinLabel(i,Cuts[i-1]);
        h_ttbarV_ee->GetXaxis()->SetBinLabel(i,Cuts[i-1]);
        h_ttbarV_mumu->GetXaxis()->SetBinLabel(i,Cuts[i-1]);
        h_data_ee->GetXaxis()->SetBinLabel(i,Cuts[i-1]);
        h_data_mumu->GetXaxis()->SetBinLabel(i,Cuts[i-1]);


  }



   //Drawing the histograms
   //tZq (ee)
   h_tZq_ee->GetXaxis()->SetTitle("Cut");
   h_tZq_ee->GetXaxis()->SetTitleOffset(1.5);
   h_tZq_ee->GetYaxis()->SetTitle("Events");
   h_tZq_ee->GetYaxis()->SetTitleOffset(1.0);
   h_tZq_ee->SetFillColor(600); //kBlue
   h_tZq_ee->SetLineColor(600); //kBlue
   h_tZq_ee->SetStats(0);
   h_tZq_ee->Draw("HIST");

   //ZPlusJets (ee)
   h_ZPlusJets_ee->GetXaxis()->SetTitle("Cut");
   h_ZPlusJets_ee->GetXaxis()->SetTitleOffset(1.5);
   h_ZPlusJets_ee->GetYaxis()->SetTitle("Events");
   h_ZPlusJets_ee->GetYaxis()->SetTitleOffset(1.0);
   h_ZPlusJets_ee->SetFillColor(632); //kRed
   h_ZPlusJets_ee->SetLineColor(632); //KRed
   h_ZPlusJets_ee->SetStats(0);
   h_ZPlusJets_ee->Draw("HIST SAME");

   //ttbar (ee)
   h_ttbar_ee->GetXaxis()->SetTitle("Cut");
   h_ttbar_ee->GetXaxis()->SetTitleOffset(1.5);
   h_ttbar_ee->GetYaxis()->SetTitle("Events");
   h_ttbar_ee->GetYaxis()->SetTitleOffset(1.0);
   h_ttbar_ee->SetFillColor(418); //kGreen+2
   h_ttbar_ee->SetLineColor(418); //kGreen+2
   h_ttbar_ee->SetStats(0);
   h_ttbar_ee->Draw("HIST SAME");

   //single top (ee)
   h_singletop_ee->GetXaxis()->SetTitle("Cut");
   h_singletop_ee->GetXaxis()->SetTitleOffset(1.5);
   h_singletop_ee->GetYaxis()->SetTitle("Events");
   h_singletop_ee->GetYaxis()->SetTitleOffset(1.0);
   h_singletop_ee->SetFillColor(618); //kMagenta+2
   h_singletop_ee->SetLineColor(618); //kMagenta+2
   h_singletop_ee->SetStats(0);
   h_singletop_ee->Draw("HIST SAME");

   //VV (ee)
   h_VV_ee->GetXaxis()->SetTitle("Cut");
   h_VV_ee->GetXaxis()->SetTitleOffset(1.5);
   h_VV_ee->GetYaxis()->SetTitle("Events");
   h_VV_ee->GetYaxis()->SetTitleOffset(1.0);
   h_VV_ee->SetFillColor(602); //kOrange+2
   h_VV_ee->SetLineColor(602); //kOrange+2
   h_VV_ee->SetStats(0);
   h_VV_ee->Draw("HIST SAME");

   //VVV (ee)
   h_VVV_ee->GetXaxis()->SetTitle("Cut");
   h_VVV_ee->GetXaxis()->SetTitleOffset(1.5);
   h_VVV_ee->GetYaxis()->SetTitle("Events");
   h_VVV_ee->GetYaxis()->SetTitleOffset(1.0);
   h_VVV_ee->SetFillColor(400); //kYellow
   h_VVV_ee->SetLineColor(400); //kYellow
   h_VVV_ee->SetStats(0);
   h_VVV_ee->Draw("HIST SAME");

   //W+jets (ee)
   h_WPlusJets_ee->GetXaxis()->SetTitle("Cut");
   h_WPlusJets_ee->GetXaxis()->SetTitleOffset(1.5);
   h_WPlusJets_ee->GetYaxis()->SetTitle("Events");
   h_WPlusJets_ee->GetYaxis()->SetTitleOffset(1.0);
   h_WPlusJets_ee->SetFillColor(425); //kCyan-7 
   h_WPlusJets_ee->SetLineColor(425); //kCyan-7
   h_WPlusJets_ee->SetStats(0);
   h_WPlusJets_ee->Draw("HIST SAME");

   //ttbarV (ee)
   h_ttbarV_ee->GetXaxis()->SetTitle("Cut");
   h_ttbarV_ee->GetXaxis()->SetTitleOffset(1.5);
   h_ttbarV_ee->GetYaxis()->SetTitle("Events");
   h_ttbarV_ee->GetYaxis()->SetTitleOffset(1.0);
   h_ttbarV_ee->SetFillColor(920); //kGray 
   h_ttbarV_ee->SetLineColor(920); //kGray
   h_ttbarV_ee->SetStats(0);
   h_ttbarV_ee->Draw("HIST SAME");

   //data (ee) 
   h_data_ee->GetXaxis()->SetTitle("Cut");
   h_data_ee->GetXaxis()->SetTitleOffset(1.5);
   h_data_ee->GetYaxis()->SetTitle("Events");
   h_data_ee->GetYaxis()->SetTitleOffset(1.0);
   h_data_ee->SetMarkerColor(1); //kBlack
   h_data_ee->SetMarkerStyle(20);
   h_data_ee->SetMarkerSize(1.0);
   h_data_ee->SetStats(0);
   h_data_ee->Draw("SAME");

   //tZq (mumu)
   h_tZq_mumu->GetXaxis()->SetTitle("Cut");
   h_tZq_mumu->GetXaxis()->SetTitleOffset(1.5);
   h_tZq_mumu->GetYaxis()->SetTitle("Events");
   h_tZq_mumu->GetYaxis()->SetTitleOffset(1.0);
   h_tZq_mumu->SetFillColor(600); //kBlue
   h_tZq_mumu->SetLineColor(600); //kBlue
   h_tZq_mumu->SetStats(0);
   h_tZq_mumu->Draw("HIST");


   //ttbar (mumu)
   h_ttbar_mumu->GetXaxis()->SetTitle("Cut");
   h_ttbar_mumu->GetXaxis()->SetTitleOffset(1.5);
   h_ttbar_mumu->GetYaxis()->SetTitle("Events");
   h_ttbar_mumu->GetYaxis()->SetTitleOffset(1.0);
   h_ttbar_mumu->SetFillColor(418); //kGreen+2
   h_ttbar_mumu->SetLineColor(418); //kGreen+2
   h_ttbar_mumu->SetStats(0);
   h_ttbar_mumu->Draw("HIST SAME");

   //ZPlusJets (mumu)
   h_ZPlusJets_mumu->GetXaxis()->SetTitle("Cut");
   h_ZPlusJets_mumu->GetXaxis()->SetTitleOffset(1.5);
   h_ZPlusJets_mumu->GetYaxis()->SetTitle("Events");
   h_ZPlusJets_mumu->GetYaxis()->SetTitleOffset(1.0);
   h_ZPlusJets_mumu->SetFillColor(632); //kRed
   h_ZPlusJets_mumu->SetLineColor(632); //kRed
   h_ZPlusJets_mumu->SetStats(0);
   h_ZPlusJets_mumu->Draw("HIST SAME"); 

   //single top (mumu)
   h_singletop_mumu->GetXaxis()->SetTitle("Cut");
   h_singletop_mumu->GetXaxis()->SetTitleOffset(1.5);
   h_singletop_mumu->GetYaxis()->SetTitle("Events");
   h_singletop_mumu->GetYaxis()->SetTitleOffset(1.0);
   h_singletop_mumu->SetFillColor(618); //kMagenta+2
   h_singletop_mumu->SetLineColor(618); //kMagenta+2
   h_singletop_mumu->SetStats(0);
   h_singletop_mumu->Draw("HIST SAME");


   //VV (mumu)
   h_VV_mumu->GetXaxis()->SetTitle("Cut");
   h_VV_mumu->GetXaxis()->SetTitleOffset(1.5);
   h_VV_mumu->GetYaxis()->SetTitle("Events");
   h_VV_mumu->GetYaxis()->SetTitleOffset(1.0);
   h_VV_mumu->SetFillColor(802); //kOrange+2
   h_VV_mumu->SetLineColor(802); //kOrange+2
   h_VV_mumu->SetStats(0);
   h_VV_mumu->Draw("HIST SAME");
 
   //VVV (mumu)
   h_VVV_mumu->GetXaxis()->SetTitle("Cut");
   h_VVV_mumu->GetXaxis()->SetTitleOffset(1.5);
   h_VVV_mumu->GetYaxis()->SetTitle("Events");
   h_VVV_mumu->GetYaxis()->SetTitleOffset(1.0);
   h_VVV_mumu->SetFillColor(400); //kYellow
   h_VVV_mumu->SetLineColor(400); //kYellow
   h_VVV_mumu->SetStats(0);
   h_VVV_mumu->Draw("HIST SAME");

   //W+jets (mumu)
   h_WPlusJets_mumu->GetXaxis()->SetTitle("Cut");
   h_WPlusJets_mumu->GetXaxis()->SetTitleOffset(1.5);
   h_WPlusJets_mumu->GetYaxis()->SetTitle("Events");
   h_WPlusJets_mumu->GetYaxis()->SetTitleOffset(1.0);
   h_WPlusJets_mumu->SetFillColor(425); //kCyan-7 
   h_WPlusJets_mumu->SetLineColor(425); //kCyan-7
   h_WPlusJets_mumu->SetStats(0);
   h_WPlusJets_mumu->Draw("HIST SAME");

   //ttbarV (mumu)
   h_ttbarV_mumu->GetXaxis()->SetTitle("Cut");
   h_ttbarV_mumu->GetXaxis()->SetTitleOffset(1.5);
   h_ttbarV_mumu->GetYaxis()->SetTitle("Events");
   h_ttbarV_mumu->GetYaxis()->SetTitleOffset(1.0);
   h_ttbarV_mumu->SetFillColor(920); //kGray
   h_ttbarV_mumu->SetLineColor(920); //kGray
   h_ttbarV_mumu->SetStats(0);
   h_ttbarV_mumu->Draw("HIST SAME");

   //data (mumu) 
   h_data_mumu->GetXaxis()->SetTitle("Cut");
   h_data_mumu->GetXaxis()->SetTitleOffset(1.5);
   h_data_mumu->GetYaxis()->SetTitle("Events");
   h_data_mumu->GetYaxis()->SetTitleOffset(1.0);
   h_data_mumu->SetMarkerColor(1); //kBlack
   h_data_mumu->SetMarkerStyle(20);
   h_data_mumu->SetMarkerSize(1.0);
   h_data_mumu->SetStats(0);
   h_data_mumu->Draw("SAME");


   //Legend
   gPad->BuildLegend();
   //Remove title
   gStyle->SetOptTitle(0);

   //Saving the histograms in pdf files
   h_tZq_ee->SaveAs("CutFlow_ee.pdf");
   h_tZq_mumu->SaveAs("CutFlow_mumu.pdf");

}

