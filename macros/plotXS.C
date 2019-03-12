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
#include "TH1F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPaletteAxis.h"
#include "TGraphErrors.h"

#include "TDRStyle.h"

int plotXS(){//main
  
  SetTdrStyle();

  TFile *file = TFile::Open("histXS.root");
  file->cd();
  
  TCanvas *myc = new TCanvas("myc","myc",1);
  myc->cd();
  gPad->SetGridy(1);
  gStyle->SetOptStat(0);

  TLegend *leg = new TLegend(0.55,0.42,0.90,0.73);
  leg->SetFillColor(10);
  
  TH1F *hM = (TH1F*)gDirectory->Get("DPxsMeson");
  TH1F *hP = (TH1F*)gDirectory->Get("DPxsPbrem");
  TH1F *hQ = (TH1F*)gDirectory->Get("DPxsQCD");

  hM->SetLineColor(2);
  hP->SetLineColor(8);
  hQ->SetLineColor(4);
  hM->SetLineWidth(2);
  hP->SetLineWidth(2);
  hQ->SetLineWidth(2);

  hQ->GetXaxis()->SetRangeUser(0.1,4.9);
  hQ->Draw("hist");

  hQ->GetXaxis()->SetTitle("m_{#gamma^{D}} (GeV.c^{-2})");
  
  hM->Draw("histsame");
  hP->Draw("histsame");
  
  leg->AddEntry(hM,"Meson","L");
  leg->AddEntry(hP,"Pbrem","L");
  leg->AddEntry(hQ,"QCD","L");
 
  leg->Draw("same");

  myc->Update();  
  myc->Print("figures/XSvsMass.pdf");
  
  return 0;


};//main
