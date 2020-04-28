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
	vector<int> indices{};

	string CSVv2OperatingPointTest = "3"; //no spaces after start
  	string MeasurementTypeTest = " iterativefit"; //one space after start
 	string SysTypeTest = " central"; // one space after start
	string JetFlavourTest = " 2"; //one space after start
        string EtaTest = " 1.7"; //one space after start
        string PtTest = " 61"; //one space after start
	string PtTestMax = "  61"; //two spaces after start
        string DiscrTest = " 0.8"; //one space after start


	for(std::vector<std::string> vec : dataList)
	{
		for(std::string data : vec)
		{	
			OutputVec.push_back(data);

		}
	}


	for(int i = 1; i < 1918; i++){

		int n;
		if(i == 1){n = 10*i; indices.push_back(n);}
		else{n = indices.at(i-2) + 11; indices.push_back(n);}

		string outputstring = OutputVec.at(n);
		outputstring.erase(outputstring.begin()+1);
		outputstring.erase(outputstring.begin());
		outputstring.erase(outputstring.end()-2);
		outputstring.erase(outputstring.end()-1);

		string::size_type pos = 0;
 
    		while ((pos = outputstring.find('x', pos)) != string::npos)
    		{
        		outputstring.replace(pos, 1, PtTest);
        		pos += 2;
    		}
	
		outputstringvec.push_back(outputstring);

	}

vector<string> OutVec{};

	for(std::vector<std::string> vec : dataList)
        {
                for(std::string data : vec)
                {

		        cout << "vec.at(0).length() = " << vec.at(0).length() << endl;
                        cout << "CSVv2OperatingPointTest.length() = " << CSVv2OperatingPointTest.length() << endl;
			cout << "vec.at(1).length() = " << vec.at(1).length() << endl;
			cout << "MeasurementTypeTest.length() = " << MeasurementTypeTest.length() << endl;
			cout << "vec.at(2).length() = " << vec.at(2).length() << endl;
                        cout << "SysTypeTest.length() = " << SysTypeTest.length() << endl;
			cout << "vec.at(3).length() = " << vec.at(3).length() << endl;
			cout << "JetFlavourTest.length() = " << JetFlavourTest.length() << endl;
			cout << "vec.at(4).length() = " << vec.at(4).length() << endl;
                        cout << "EtaTest.length() = " << EtaTest.length() << endl;
			cout << "vec.at(5).length() = " << vec.at(5).length() << endl;
                        cout << "EtaTest.length() = " << EtaTest.length() << endl;
			cout << "vec.at(6).length() = " << vec.at(6).length() << endl;
			cout << "vec.at(6) = " << vec.at(6) << endl;
                        cout << "PtTest.length() = " << PtTest.length() << endl;
			cout << "PtTest = " << PtTest << endl;
		  	cout << "vec.at(7).length() = " << vec.at(7).length() << endl;
                        cout << "vec.at(7) = " << vec.at(7) << endl;
                        cout << "PtTestMax.length() = " << PtTestMax.length() << endl;
                        cout << "PtTestMax = " << PtTestMax << endl;
			cout << "vec.at(8).length() = " << vec.at(8).length() << endl;
                        cout << "vec.at(8) = " << vec.at(8) << endl;
                        cout << "DiscrTest.length() = " << DiscrTest.length() << endl;
                        cout << "DiscrTest = " << DiscrTest << endl;
			cout << "vec.at(9).length() = " << vec.at(9).length() << endl;
                        cout << "vec.at(9) = " << vec.at(9) << endl;
                        cout << "DiscrTest.length() = " << DiscrTest.length() << endl;
                        cout << "DiscrTest = " << DiscrTest << endl;
		
	
			if( (vec.at(0) == CSVv2OperatingPointTest) 
    			&& (vec.at(1) == MeasurementTypeTest)
    			&& (vec.at(2) == SysTypeTest)
    	  		&& (vec.at(3) == JetFlavourTest)
    			&& (vec.at(4) < EtaTest)
    			&& (vec.at(5) > EtaTest)
    	  		&& (vec.at(6) < PtTest)     
    			&& (vec.at(7) > PtTestMax)
    			&& (vec.at(8) < DiscrTest) 
    			&& (vec.at(9) > DiscrTest)
  			){
				OutVec.push_back(vec.at(10)); 						
			}
			else if( (vec.at(0) != CSVv2OperatingPointTest)
                        || (vec.at(1) != MeasurementTypeTest)
                        || (vec.at(2) != SysTypeTest)
                        || (vec.at(3) != JetFlavourTest)
                        || (vec.at(4) > EtaTest)
                        || (vec.at(5) < EtaTest)
                        || (vec.at(6) > PtTest)     
                        || (vec.at(7) < PtTestMax)
                        || (vec.at(8) > DiscrTest) 
                        || (vec.at(9) < DiscrTest)){OutVec.push_back("0");}
			else{cout << "double check criteria" << endl;}


                }
        }



//bool any_of(OutVec.begin(), OutVec.end(), [](string i){return i != "0";});
//cout << " any_of = " << any_of << endl;


return OutVec;

}


