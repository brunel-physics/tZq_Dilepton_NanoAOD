#include "/home/eepgkkc/ScaleFactors/RochesterCorrections/roccor.Run2.v3/RoccoR.cc"
#include "/home/eepgkkc/ScaleFactors/RochesterCorrections/roccor.Run2.v3/RoccoR.h"

auto RochesterCorrections_testscript2(const string& year, const string& process){

	string RoccoTextFile;

	if(year == "2016"){RoccoTextFile = "/home/eepgkkc/ScaleFactors/RochesterCorrections/roccor.Run2.v3/RoccoR2016.txt";}
	else if(year == "2017"){RoccoTextFile = "/home/eepgkkc/ScaleFactors/RochesterCorrections/roccor.Run2.v3/RoccoR2017.txt";}
	else if(year == "2018"){RoccoTextFile = "/home/eepgkkc/ScaleFactors/RochesterCorrections/roccor.Run2.v3/RoccoR2018.txt";}
	else{cout << "Error for rochester corrections: choose a year out of 2016, 2017 or 2018." << endl;}


	RoccoR rc{RoccoTextFile};
	vector<double> RochCorrVec{};

	vector<int> MuonCharge{1, -1};
	vector<float> MuonPt{150.0, 95.6};
	vector<float> MuonEta{2.0, -0.01};
	vector<float> MuonPhi{1.7, -2.6};
	vector<int> s{0, 0}; //s is error set (default is 0)
	vector<int> m{0, 0}; //m is error member (default is 0, ranges from 0 to nmembers-1)
	vector<int> MuonGen{1, 0};
	vector<double> u{gRandom->Rndm(), gRandom->Rndm()}; //u is a random number distributed uniformly between 0 and 1 (gRandom->Rndm());
	vector<int> Muon_nTrackerLayers{4, 6}; //trackerLayersWithMeasurement (nl in the README)


	for(int i = 0; i < MuonPt.size(); i++){

		//scale factors for momentum of each muon:
		double RochCorrSF;
		double mcSF;


		if(process != "data_DoubleEGRun2017B" &&
		   process != "data_DoubleEGRun2017C" &&
		   process != "data_DoubleEGRun2017D" &&
		   process != "data_DoubleEGRun2017E" &&
		   process != "data_DoubleEGRun2017F" &&
		   process != "data_DoubleMuonRun2017B" &&
		   process != "data_DoubleMuonRun2017C" &&
		   process != "data_DoubleMuonRun2017D" &&
		   process != "data_DoubleMuonRun2017E" &&
		   process != "data_DoubleMuonRun2017F"){

			if(mcSF > 0){

				RochCorrSF = rc.kSpreadMC(MuonCharge.at(i), MuonPt.at(i), MuonEta.at(i), MuonPhi.at(i), MuonGen.at(i), s.at(i), m.at(i)); //(recommended), MC scale and resolution correction when matched gen muon is available
			}
			else{
				RochCorrSF = rc.kSmearMC(MuonCharge.at(i), MuonPt.at(i), MuonEta.at(i), MuonPhi.at(i), Muon_nTrackerLayers.at(i), u.at(i), s.at(i), m.at(i)); //MC scale and extra smearing when matched gen muon is not available

			}

		}
		else{RochCorrSF = rc.kScaleDT(MuonCharge.at(i), MuonPt.at(i), MuonEta.at(i), MuonPhi.at(i), s.at(i), m.at(i));} //data

	
		RochCorrVec.push_back(RochCorrSF);


	}

	
	return RochCorrVec;


}



auto RochesterCorrections_testscript(){

	return RochesterCorrections_testscript2("2017", "data_DoubleEGRun2017B");
	//return RochesterCorrections_testscript2("2017", "tZq");

}
