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

int plotPbremPDF(){//main
  
  SetTdrStyle();

  const unsigned nF = 2;
  TFile *file[nF];
  //file[0] = TFile::Open("ParaPhoton_eps5e-06_m0.5.root");
  //file[1] = TFile::Open("ParaPhoton_eps1e-08_m1.1.root");
  file[0] = TFile::Open("ParaPhoton_eps1e-07_m0.3.root");
  file[1] = TFile::Open("ParaPhoton_eps1e-07_m2.0.root");
 
  //std::string masseps[nF] = {"eps5e-06_m0.5","eps1e-08_m1.1"};
  std::string masseps[nF] = {"eps1e-07_m0.3","eps1e-07_m2.0"};
  
  TCanvas *myc = new TCanvas("myc","myc",1000,500);
  myc->Divide(2,1);

  TH2F *hPDF[nF];

  gStyle->SetOptStat(0);
  
  for (unsigned iF(0); iF<nF; ++iF){
    file[iF]->cd();
    hPDF[iF] = (TH2F*)gDirectory->Get(("hPDF_"+masseps[iF]).c_str());
    myc->cd(iF+1);
    gPad->SetRightMargin(0.15);
    gPad->SetLeftMargin(0.15);
    gPad->SetTopMargin(0.1);

    
    hPDF[iF]->Rebin2D(10,1);
    hPDF[iF]->Draw("Surf3z");
    hPDF[iF]->SetTitle(";p (GeV); #theta (rad); f(p,#theta)");
    hPDF[iF]->GetYaxis()->SetRangeUser(-0.1,0.1);

    hPDF[iF]->GetXaxis()->SetTitleOffset(1.5);
    hPDF[iF]->GetYaxis()->SetTitleOffset(2);
    hPDF[iF]->GetZaxis()->SetTitleOffset(1.5);

    TLatex lat;
    if (iF==0) lat.DrawLatexNDC(0.2,0.93,"m_{#gamma^{D}} = 0.3 GeV");
    else lat.DrawLatexNDC(0.2,0.93,"m_{#gamma^{D}} = 2.0 GeV");
    
    gPad->Update();
    TPaletteAxis *palette = (TPaletteAxis*)hPDF[iF]->GetListOfFunctions()->FindObject("palette");

    palette->SetX1NDC(0.88);
    palette->SetX2NDC(0.93);
    gPad->Modified();
    gPad->Update();
    

  }
  
  myc->Update();
  myc->Print("figures/pbrem_PDF.pdf");
  
  return 0;


};//main
