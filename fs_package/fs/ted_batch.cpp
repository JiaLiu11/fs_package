#include "EOS.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

using namespace std;

int main(void)
{
  EOS eos;
  eos.loadEOSFromFile("s95p-PCE/EOS_converted.dat", "s95p-PCE/coeff.dat");  //load S95 EOS table

  ifstream inFile("tdec_table.dat");
  ofstream of;
  of.open("edec_converted.dat", std::ios_base::out);

  double t_temp, ed_temp;
  while (!inFile.eof())
  {
    inFile >> t_temp;  //safty procedure, reach the end of the file
    t_temp = t_temp/1000.0; // convert to MeV
    ed_temp= eos.edFromT(t_temp);
    of << setprecision(12) << ed_temp << endl;
  }
  inFile.close();
  of.close();
}
