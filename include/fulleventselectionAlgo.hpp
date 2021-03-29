#ifndef _fulleventselectionAlgo_h_
#define _fulleventselectionAlgo_h_


class fulleventselectionAlgo/* : public edm::EDAnalyzer*/
{

  public:
  //Constructor
  fulleventselectionAlgo();
  ~fulleventselectionAlgo();

  void parseCommandLineArguements(int argc, char* argv[]); 
  void fulleventselection();

  //private:
  //edm::EDGetTokenT<LHEEventProduct> externalLHEToken_;

  int MC_Selection_;
  int Year_Selection_;
  int Process_Selection_;
  int NPL_Selection_;
  int SR_Selection_;
  int SBR_Selection_;
  int ZPlusJetsCR_Selection_;
  int ttbarCR_Selection_;
  int Systematic_Selection_;
  int Channel_Selection_;
  int DoubleCountCheck_Selection_;


};


#endif
