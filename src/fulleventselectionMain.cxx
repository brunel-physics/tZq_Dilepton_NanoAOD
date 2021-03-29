#include "fulleventselectionAlgo.hpp"

#include <iostream>

int main(int argc, char* argv[]){

  fulleventselectionAlgo fulleventselectionMain;

  fulleventselectionMain.parseCommandLineArguements(argc, argv);
  fulleventselectionMain.fulleventselection();


}
