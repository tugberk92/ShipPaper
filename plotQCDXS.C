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
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPaletteAxis.h"
#include "TGraphErrors.h"

#include "TDRStyle.h"

double xsParam(double *mass, double*p){
  if (mass[0] > 3) return exp(-5.673-0.8869*mass[0]);
  else if (mass[0] > 1.4) return exp(-3.802-1.532*mass[0]);
  else return 0.0586-0.09037*mass[0] + 0.0360743*mass[0]*mass[0];
};

int plotQCDXS(){//main
  
  SetTdrStyle();

  TFile *file = TFile::Open("crossSectionsQCDbis.root");
  file->cd();
  
  const unsigned nE = 6;
  double epsilon[nE] = {1e-7,8e-8,6e-8,4e-8,2e-8,1e-8};

  int color[nE] = {1,2,3,4,6,7};
  
  TCanvas *myc = new TCanvas("myc","myc",1);
  myc->cd();
  gPad->SetLogy(1);
  gStyle->SetOptStat(0);

  TGraphErrors *grE[nE];

  TLegend *leg = new TLegend(0.66,0.5,0.92,0.92);
  leg->SetFillColor(10);
  
  for (unsigned iE(0); iE<nE; ++iE){//loop on epsilon
    std::ostringstream lStr;
    lStr<<"grNorm_"<<iE;
    grE[iE] = (TGraphErrors*)gDirectory->Get(lStr.str().c_str());

    grE[iE]->SetLineColor(color[iE]);
    grE[iE]->SetMarkerColor(color[iE]);
    grE[iE]->SetMarkerStyle(20+iE);
    grE[iE]->SetTitle(";M_{#gamma^{D}} (GeV); #sigma_{QCD}/#epsilon^{2} (mb)");
    grE[iE]->Draw(iE==0?"APE":"PEsame");
    grE[iE]->GetXaxis()->SetRangeUser(1,6.2);
    grE[iE]->GetYaxis()->SetRangeUser(0.00001,0.01);

    lStr.str("");
    lStr << "#epsilon = " << epsilon[iE];
    leg->AddEntry(grE[iE],lStr.str().c_str(),"P");
  }

  TF1 *myParam = new TF1("myParam",xsParam,1,10,0);
  myParam->SetLineWidth(2);
  myParam->SetLineColor(9);
  myParam->Draw("same");

  leg->AddEntry(myParam,"Fit","L");
  leg->Draw("same");

  myc->Update();  
  myc->Print("figures/qcdXSnorm_massFitBis.pdf");
  
  return 0;


};//main
