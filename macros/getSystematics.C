#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

#include "TCanvas.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLine.h"


double mesonSyst(const double & mass){

  double syst = 0;
  //syst on BR meson->gamma decay modes
  if (mass < 0.135) syst += pow(0.0003,2)+pow(0.0005,2);
  else if (mass < 0.548) syst += pow(0.005,2)+pow(0.011,2);
  else if (mass<0.648) syst += pow(0.034,2)+pow(0.011,2);
  else syst += pow(0.036,2)+pow(0.038,2);

  //QCD+PDF on sigma non-diff - placeholder
  syst += pow(0.05,2);

  return sqrt(syst);
  
}

double syst(const std::string & process, const double & mass, int sign=1){
  if (process.find("meson") != process.npos)
    return mesonSyst(mass);
  else if (process.find("qcd") != process.npos){
    if (sign>0) return 0.05;
    else return 0.15;
  }
  else return 0.2;
}
