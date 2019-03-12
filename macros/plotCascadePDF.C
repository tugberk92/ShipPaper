#include <iostream>
#include <fstream> 
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
#include <iomanip>


#include "TFile.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPaletteAxis.h"

#include "TDRStyle.h"

int plotCascadePDF(){//main
  
  SetTdrStyle();

  TFile *fileIn = TFile::Open("CascadeParticles_500evts.root");

  
  const unsigned nP = 5;
  std::string particle[nP] = {"Pi0","Eta","Omega","EtaPrime","Proton"};

  TCanvas *myc = new TCanvas("myc","myc",1);
  TCanvas *mycP = new TCanvas("mycP","mycP",1);

  TH2F *hPDF[nP];
  TH2F *hPDFprimary[nP];

  gStyle->SetOptStat(0);
  
  fileIn->cd();
  TTree *mesonTree = (TTree*)gDirectory->Get("Mesons");
  TTree *protonTree = (TTree*)gDirectory->Get("Protons");
  
  
  for (unsigned iP(0); iP<nP; ++iP){
    std::ostringstream lVar,lCut;
    lVar << "theta" << particle[iP] << ":p" << particle[iP]
	 << ">>h" << iP << "(100,0,100,100,-1,1)";
    lCut << "isPrimary" << particle[iP] << " == 1";
