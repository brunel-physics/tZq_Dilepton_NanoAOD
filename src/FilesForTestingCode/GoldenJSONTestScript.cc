#include <fstream>
#include <string>
#include<iostream>
#include<fstream>

using namespace std;



auto GoldenJsonReader(const string& year) {

 string GoldenJsonFileName;

 if(year == "2016"){GoldenJsonFileName = "/home/eepgkkc/ScaleFactors/GoldenJSON/Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt";}
 else if(year == "2017"){GoldenJsonFileName = "/home/eepgkkc/ScaleFactors/GoldenJSON/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt";}
 else if(year == "2018"){GoldenJsonFileName = "/home/eepgkkc/ScaleFactors/GoldenJSON/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt";}
 else{cout << "Choose the year out of 2016, 2017 or 2018" << endl;}


 ifstream myReadFile;
 myReadFile.open(GoldenJsonFileName);
 static char output[100];
 vector<string> GoldenJsonOutput{};
  

 if (myReadFile.is_open()) {
 while (!myReadFile.eof()) {


    myReadFile >> output;
    GoldenJsonOutput.push_back(output);

 }
}

myReadFile.close();
return GoldenJsonOutput;

}




auto GoldenJson_SplitChars(const string& year){

  vector<char> out{};

  for(int i = 0; i < (GoldenJsonReader(year)).size(); i++){

  	string element = GoldenJsonReader(year).at(i);

  	for(int j = 0; j < element.size(); j++){

  		out.push_back(element.at(j));  	

  	}

  }

  return out;

}

int InputRunNumber = 297101; 

auto RunNumberCheck(const string& year){

  
 vector<char> EventsVector{}; 

 for(int i = 0; i < (GoldenJson_SplitChars(year)).size(); i++){

 	if(  GoldenJson_SplitChars(year).at(i+1) == '"' &&
	    (GoldenJson_SplitChars(year).at(i+2) == '2' ||
	     GoldenJson_SplitChars(year).at(i+2) == '3')  ){ 

		int digit1 = GoldenJson_SplitChars(year).at(i+2) - '0'; 
		int digit2 = GoldenJson_SplitChars(year).at(i+3) - '0'; 
		int digit3 = GoldenJson_SplitChars(year).at(i+4) - '0'; 
		int digit4 = GoldenJson_SplitChars(year).at(i+5) - '0'; 
		int digit5 = GoldenJson_SplitChars(year).at(i+6) - '0'; 
		int digit6 = GoldenJson_SplitChars(year).at(i+7) - '0';

		cout << "digit1 = " << digit1 << endl;
		cout << "digit2 = " << digit2 << endl;
		cout << "digit3 = " << digit3 << endl;
		cout << "digit4 = " << digit4 << endl;
		cout << "digit5 = " << digit5 << endl;
		cout << "digit6 = " << digit6 << endl;

		int run = (digit1*100000) + (digit2*10000) + (digit3*1000) + (digit4*100) + (digit5*10) + digit6;

		cout << "The run number being read is: " << run << endl;

		if(run == InputRunNumber){

			cout << "The Golden Json run number and the input run number are the same" << endl;

			for(int j = 2; j < (GoldenJson_SplitChars(year)).size(); j++){

				if(GoldenJson_SplitChars(year).at(i+10) == '[' && 
				   GoldenJson_SplitChars(year).at(i+11) == '['){	

					if( GoldenJson_SplitChars(year).at( (i+10)+j ) == ']' &&
					    GoldenJson_SplitChars(year).at( (i+10)+(j+1) ) == ']'){

					
						for(int k = (i+10); k < ((i+10)+(j+2)); k++){

							EventsVector.push_back(GoldenJson_SplitChars(year).at(k));
						
						}

						return EventsVector;	
											

					}
					else{cout << "GoldenJson_SplitChars(year).at((i+11)+j) and at (i+11)+(j+1) are not ']' " << endl; continue;}


				}
				else{cout << "GoldenJson_SplitChars(year).at(i+10) is not '[' " << endl; continue;}


			}	

		}
		else{continue;}

	}
	else{cout << "The run number does not match the input run number" << endl;}


 }


}


auto ReturnRunNumAndEventRanges(const string& year){

 vector<int> RunNumAndEvents{};

 RunNumAndEvents.push_back(InputRunNumber);

 vector<char> Runs = RunNumberCheck(year);


 for(int i = 0; i < Runs.size(); i++){


 	if(Runs.at(i) == ']' && Runs.at(i+1) == ']'){break;}
	else if( isdigit(Runs.at(i)) || Runs.at(i) == ',' || Runs.at(i) == ' ' || Runs.at(i) == ']'){continue;}
 	else if( (Runs.at(i) == '[' && Runs.at(i+1) == '[') ||
		 (Runs.at(i) == '[' && isdigit(Runs.at(i+1))) ){

		cout << "inside the third else if" << endl;

		if(  isdigit( Runs.at(i+1) ) && //For the min value being a 4 digit number
      	     	     isdigit( Runs.at(i+2) ) &&
      	     	     isdigit( Runs.at(i+3) ) && 
		     isdigit( Runs.at(i+4) ) &&
      	     	     Runs.at(i+5) == ','){

			cout << "inside if for the min value being a 4 digit number" << endl;

			int Min_Digit1 = Runs.at(i+1) - '0';
			int Min_Digit2 = Runs.at(i+2) - '0';
			int Min_Digit3 = Runs.at(i+3) - '0';
			int Min_Digit4 = Runs.at(i+4) - '0';

			int EventNumber_Min = (Min_Digit1*1000) + (Min_Digit2*100) + (Min_Digit3*10) + Min_Digit4;
			RunNumAndEvents.push_back(EventNumber_Min);

			if(isdigit( Runs.at(i+6) ) &&
                           isdigit( Runs.at(i+7) ) &&
                           isdigit( Runs.at(i+8) ) &&
                           isdigit( Runs.at(i+9) ) &&
                           Runs.at(i+10) == ']' ){ //For the min value being a 4 digit number and the max value being a 4 digit number
                                   
                                        cout << "Runs.at(i+6) = " << Runs.at(i+6) << endl;
                                        cout << "Runs.at(i+7) = " << Runs.at(i+7) << endl;
                                        cout << "Runs.at(i+8) = " << Runs.at(i+8) << endl;
                                        cout << "Runs.at(i+9) = " << Runs.at(i+9) << endl;
                                        cout << "Runs.at(i+10) = " << Runs.at(i+10) << endl;

                                        int Max_Digit1 = Runs.at(i+6) - '0';
                                        int Max_Digit2 = Runs.at(i+7) - '0';
                                        int Max_Digit3 = Runs.at(i+8) - '0';
                                        int Max_Digit4 = Runs.at(i+9) - '0';

                                        int EventNumber_Max = (Max_Digit1*1000) + (Max_Digit2*100) + (Max_Digit3*10) + Max_Digit4;
                                        RunNumAndEvents.push_back(EventNumber_Max);


                                }
			else if(isdigit( Runs.at(i+6) ) && //For the min value being a 4 digit number and the max value being a 3 digit number
      	   	   	        isdigit( Runs.at(i+7) ) &&
      	   	   	        isdigit( Runs.at(i+8) ) &&
      	   	   	        Runs.at(i+9) == ']'){

				cout << "Runs.at(i+6) = " << Runs.at(i+6) << endl;
				cout << "Runs.at(i+7) = " << Runs.at(i+7) << endl;
				cout << "Runs.at(i+8) = " << Runs.at(i+8) << endl;

				int Max_Digit1 = Runs.at(i+6) - '0';
        			int Max_Digit2 = Runs.at(i+7) - '0';
        			int Max_Digit3 = Runs.at(i+8) - '0';
	
				int EventNumber_Max = (Max_Digit1*100) + (Max_Digit2*10) + Max_Digit3;
        			RunNumAndEvents.push_back(EventNumber_Max);
		

			}
			else if(isdigit( Runs.at(i+6) ) && //For the min value being a 4 digit number and the max value being a 2 digit number 
         			isdigit( Runs.at(i+7) ) &&
         			Runs.at(i+8) == ']' ){

					int Max_Digit1 = Runs.at(i+6) - '0';
                			int Max_Digit2 = Runs.at(i+7) - '0';

                			int EventNumber_Max = (Max_Digit1*10) + Max_Digit2;
                			RunNumAndEvents.push_back(EventNumber_Max);


			}
			else if(isdigit( Runs.at(i+6) ) &&
                                Runs.at(i+7) == ']' ){ //For the min value being a 4 digit number and the max value being a 1 digit number

                                        cout << "Runs.at(i+6) = " << Runs.at(i+6) << endl;
                                        cout << "Runs.at(i+7) = " << Runs.at(i+7) << endl;

                                        int Max_Digit1 = Runs.at(i+6) - '0';

                                        int EventNumber_Max = Max_Digit1;
                                        RunNumAndEvents.push_back(EventNumber_Max);


                        }
                       else{cout << "error" << endl;}

 		}	
 		else if(  isdigit( Runs.at(i+1) ) && //For the min value being a 3 digit number
      	     	     isdigit( Runs.at(i+2) ) &&
      	     	     isdigit( Runs.at(i+3) ) && 
      	     	     Runs.at(i+4) == ','){

			cout << "inside else if for the min value being a 3 digit number" << endl;

			int Min_Digit1 = Runs.at(i+1) - '0';
			int Min_Digit2 = Runs.at(i+2) - '0';
			int Min_Digit3 = Runs.at(i+3) - '0';

			int EventNumber_Min = (Min_Digit1*100) + (Min_Digit2*10) + Min_Digit3;
			RunNumAndEvents.push_back(EventNumber_Min);

			if(isdigit( Runs.at(i+5) ) &&
                           isdigit( Runs.at(i+6) ) &&
                           isdigit( Runs.at(i+7) ) &&
                           isdigit( Runs.at(i+8) ) &&
                           Runs.at(i+9) == ']' ){ //For the min value being a 3 digit number and the max value being a 4 digit number
                                   
                                        cout << "Runs.at(i+5) = " << Runs.at(i+5) << endl;
                                        cout << "Runs.at(i+6) = " << Runs.at(i+6) << endl;
                                        cout << "Runs.at(i+7) = " << Runs.at(i+7) << endl;
                                        cout << "Runs.at(i+8) = " << Runs.at(i+8) << endl;
                                        cout << "Runs.at(i+9) = " << Runs.at(i+9) << endl;

                                        int Max_Digit1 = Runs.at(i+5) - '0';
                                        int Max_Digit2 = Runs.at(i+6) - '0';
                                        int Max_Digit3 = Runs.at(i+7) - '0';
                                        int Max_Digit4 = Runs.at(i+8) - '0';

					cout << "Max_Digit4 = " << Max_Digit4 << endl;

                                        int EventNumber_Max = (Max_Digit1*1000) + (Max_Digit2*100) + (Max_Digit3*10) + Max_Digit4;

					cout << "EventNumber_Max = " << EventNumber_Max << endl;
                                        RunNumAndEvents.push_back(EventNumber_Max);


                                }
			else if(isdigit( Runs.at(i+5) ) && //For the min value being a 3 digit number and the max value being a 3 digit number
      	   	   	        isdigit( Runs.at(i+6) ) &&
      	   	   	        isdigit( Runs.at(i+7) ) &&
      	   	   	        Runs.at(i+8) == ']'){

				cout << "Runs.at(i+5) = " << Runs.at(i+5) << endl;
				cout << "Runs.at(i+6) = " << Runs.at(i+6) << endl;
				cout << "Runs.at(i+7) = " << Runs.at(i+7) << endl;

				int Max_Digit1 = Runs.at(i+5) - '0';
        			int Max_Digit2 = Runs.at(i+6) - '0';
        			int Max_Digit3 = Runs.at(i+7) - '0';
	
				int EventNumber_Max = (Max_Digit1*100) + (Max_Digit2*10) + Max_Digit3;
        			RunNumAndEvents.push_back(EventNumber_Max);
		

			}
			else if(isdigit( Runs.at(i+5) ) && //For the min value being a 3 digit number and the max value being a 2 digit number 
         			isdigit( Runs.at(i+6) ) &&
         			Runs.at(i+7) == ']' ){

					int Max_Digit1 = Runs.at(i+5) - '0';
                			int Max_Digit2 = Runs.at(i+6) - '0';

                			int EventNumber_Max = (Max_Digit1*10) + Max_Digit2;
                			RunNumAndEvents.push_back(EventNumber_Max);


			}
			else if(isdigit( Runs.at(i+5) ) &&
                                Runs.at(i+6) == ']' ){ //For the min value being a 3 digit number and the max value being a 1 digit number

                                        cout << "Runs.at(i+5) = " << Runs.at(i+5) << endl;
                                        cout << "Runs.at(i+6) = " << Runs.at(i+6) << endl;

                                        int Max_Digit1 = Runs.at(i+5) - '0';

                                        int EventNumber_Max = Max_Digit1;
                                        RunNumAndEvents.push_back(EventNumber_Max);


                        }
                       else{cout << "error" << endl;}

 		}
 		else if(isdigit( Runs.at(i+1) ) && //For the min value being a 2 digit number
         		isdigit( Runs.at(i+2) ) &&
         		Runs.at(i+3) == ',' ){

				cout << "inside else if for the min value being a 2 digit number" << endl;

				int Min_Digit1 = Runs.at(i+1) - '0';
        			int Min_Digit2 = Runs.at(i+2) - '0';

        			int EventNumber_Min = (Min_Digit1*10) + Min_Digit2;
        			RunNumAndEvents.push_back(EventNumber_Min);

				if(isdigit( Runs.at(i+4) ) &&
                                   isdigit( Runs.at(i+5) ) &&
				   isdigit( Runs.at(i+6) ) &&
                                   isdigit( Runs.at(i+7) ) &&
                                   Runs.at(i+8) == ']' ){ //For the min value being a 2 digit number and the max value being a 4 digit number

                                        cout << "Runs.at(i+4) = " << Runs.at(i+4) << endl;
                                        cout << "Runs.at(i+5) = " << Runs.at(i+5) << endl;
                                        cout << "Runs.at(i+6) = " << Runs.at(i+6) << endl;
					cout << "Runs.at(i+7) = " << Runs.at(i+7) << endl;
                                        cout << "Runs.at(i+8) = " << Runs.at(i+8) << endl;

                                        int Max_Digit1 = Runs.at(i+4) - '0';
                                        int Max_Digit2 = Runs.at(i+5) - '0';
					int Max_Digit3 = Runs.at(i+6) - '0';
                                        int Max_Digit4 = Runs.at(i+7) - '0';

                                        int EventNumber_Max = (Max_Digit1*1000) + (Max_Digit2*100) + (Max_Digit3*10) * Max_Digit4;
                                        RunNumAndEvents.push_back(EventNumber_Max);


                                }
				else if(isdigit( Runs.at(i+4) ) && //For the min value being a 2 digit number and the max value being a 3 digit nummber
           	   	   	   isdigit( Runs.at(i+5) ) &&
           	   	   	   isdigit( Runs.at(i+6) ) &&
           	   	   	   Runs.at(i+7) == ']'){

					cout << "Runs.at(i+4) = " << Runs.at(i+4) << endl;
                			cout << "Runs.at(i+5) = " << Runs.at(i+5) << endl;
                			cout << "Runs.at(i+6) = " << Runs.at(i+6) << endl;

                			int Max_Digit1 = Runs.at(i+4) - '0';
                			int Max_Digit2 = Runs.at(i+5) - '0';
                			int Max_Digit3 = Runs.at(i+6) - '0';

                			int EventNumber_Max = (Max_Digit1*100) + (Max_Digit2*10) + Max_Digit3;
                			RunNumAndEvents.push_back(EventNumber_Max);


        			}
        			else if(isdigit( Runs.at(i+4) ) &&
                			isdigit( Runs.at(i+5) ) &&
                			Runs.at(i+6) == ']' ){ //For the min value being a 2 digit number and the max value being a 2 digit number

					cout << "Runs.at(i+4) = " << Runs.at(i+4) << endl;
					cout << "Runs.at(i+5) = " << Runs.at(i+5) << endl;
                        		cout << "Runs.at(i+6) = " << Runs.at(i+6) << endl;

                        		int Max_Digit1 = Runs.at(i+4) - '0';
                        		int Max_Digit2 = Runs.at(i+5) - '0';

                        		int EventNumber_Max = (Max_Digit1*10) + Max_Digit2;
                        		RunNumAndEvents.push_back(EventNumber_Max);


        			}
				else if(isdigit( Runs.at(i+4) ) &&
                                        Runs.at(i+5) == ']' ){ //For the min value being a 2 digit number and the max value being a 1 digit number

                                        cout << "Runs.at(i+4) = " << Runs.at(i+4) << endl;
                                        cout << "Runs.at(i+5) = " << Runs.at(i+5) << endl;

                                        int Max_Digit1 = Runs.at(i+4) - '0';

                                        int EventNumber_Max = Max_Digit1;
                                        RunNumAndEvents.push_back(EventNumber_Max);


                                }
				else{cout << "error" <<  "Runs.at(i+1) = " << Runs.at(i+1) << '\n' << "Runs.at(i+2) = " << Runs.at(i+2) << endl;}

   	}
	else if(  isdigit( Runs.at(i+1) ) && //For the min value being a 1 digit number
      	     	  Runs.at(i+2) == ',' &&
		  isdigit(Runs.at(i+3)) ){

   			cout << "inside else if for the min value being a 1 digit number" << endl;
			cout << "Runs.at(i+1) = " << Runs.at(i+1) << endl;

			int Min_Digit1 = Runs.at(i+1) - '0';	

			int EventNumber_Min = Min_Digit1;

			RunNumAndEvents.push_back(EventNumber_Min);

			if(isdigit( Runs.at(i+3) ) &&
                           isdigit( Runs.at(i+4) ) &&
                           isdigit( Runs.at(i+5) ) &&
                           isdigit( Runs.at(i+6) ) &&
                           Runs.at(i+7) == ']' ){ //For the min value being a 1 digit number and the max value being a 4 digit number
                                   
                                        cout << "Runs.at(i+3) = " << Runs.at(i+3) << endl;
                                        cout << "Runs.at(i+4) = " << Runs.at(i+4) << endl;
                                        cout << "Runs.at(i+5) = " << Runs.at(i+5) << endl;
                                        cout << "Runs.at(i+6) = " << Runs.at(i+6) << endl;
                                        cout << "Runs.at(i+7) = " << Runs.at(i+7) << endl;

                                        int Max_Digit1 = Runs.at(i+3) - '0';
                                        int Max_Digit2 = Runs.at(i+4) - '0';
                                        int Max_Digit3 = Runs.at(i+5) - '0';
                                        int Max_Digit4 = Runs.at(i+6) - '0';

                                        int EventNumber_Max = (Max_Digit1*1000) + (Max_Digit2*100) + (Max_Digit3*10) * Max_Digit4;
                                        RunNumAndEvents.push_back(EventNumber_Max);


                                }
			else if(isdigit( Runs.at(i+3) ) && //For the min value being a 1 digit number and the max value being a 3 digit number
      	   	   	        isdigit( Runs.at(i+4) ) &&
      	   	   	        isdigit( Runs.at(i+5) ) &&
      	   	   	        Runs.at(i+6) == ']'){

				cout << "Runs.at(i+3) = " << Runs.at(i+3) << endl;
				cout << "Runs.at(i+4) = " << Runs.at(i+4) << endl;
				cout << "Runs.at(i+5) = " << Runs.at(i+5) << endl;

				int Max_Digit1 = Runs.at(i+3) - '0';
        			int Max_Digit2 = Runs.at(i+4) - '0';
        			int Max_Digit3 = Runs.at(i+5) - '0';
	
				int EventNumber_Max = (Max_Digit1*100) + (Max_Digit2*10) + Max_Digit3;
        			RunNumAndEvents.push_back(EventNumber_Max);
		

			}
			else if(isdigit( Runs.at(i+3) ) && //For the min value being a 1 digit number and the max value being a 2 digit number 
         			isdigit( Runs.at(i+4) ) &&
         			Runs.at(i+5) == ']' ){

					int Max_Digit1 = Runs.at(i+3) - '0';
                			int Max_Digit2 = Runs.at(i+4) - '0';

                			int EventNumber_Max = (Max_Digit1*10) + Max_Digit2;
                			RunNumAndEvents.push_back(EventNumber_Max);


			}
			else if(isdigit( Runs.at(i+3) ) &&
                                Runs.at(i+4) == ']' ){ //For the min value being a 1 digit number and the max value being a 1 digit number

                                        cout << "Runs.at(i+3) = " << Runs.at(i+3) << endl;
                                        cout << "Runs.at(i+4) = " << Runs.at(i+4) << endl;

                                        int Max_Digit1 = Runs.at(i+3) - '0';

                                        int EventNumber_Max = Max_Digit1;
                                        RunNumAndEvents.push_back(EventNumber_Max);


                        }
                       else{cout << "error" << endl;}

 		}	
		else{cout << "INSIDE THE ELSE STATEMENT" << '\n' << "Runs.at(i) = " << Runs.at(i) << '\n' << "Runs.at(i+1) = " << Runs.at(i+1) << '\n' << "Runs.at(i+2) = " << Runs.at(i+2) << endl;}


	}	 


  }


  return RunNumAndEvents;


}




auto GoldenJSONTestScript(){

 return ReturnRunNumAndEventRanges("2017");

}
