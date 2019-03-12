#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLegend.h"

#include <string>
#include <vector>
#include <iostream>
#include <sstream>

#include "TDRStyle.h"

int plotBR() {//main
  
  SetTdrStyle();

  TFile *fin = TFile::Open("BRnew.root");
  fin->cd();

  TCanvas *myc = new TCanvas("myc","myc",1);
  myc->cd();

  const unsigned nH = 5;

  TH1F *hBR[nH];
  hBR[0] = (TH1F*)gDirectory->Get("hBR_e");
  hBR[1] = (TH1F*)gDirectory->Get("hBR_mu");
  hBR[2] = (TH1F*)gDirectory->Get("hBR_tau");
  hBR[3] = (TH1F*)gDirectory->Get("hBR_had");
  hBR[4] = (TH1F*)gDirectory->Get("hBR_tot");

  gStyle->SetOptStat(0);
  
  TLegend *leg = new TLegend(0.64,0.7,0.94,0.94);
  leg->SetFillColor(10);

  hBR[4]->SetTitle("");
  hBR[4]->Draw("hist");
  hBR[4]->GetXaxis()->SetTitle("m_{#gamma^{D}} (GeV.c^{-2})");
  hBR[4]->GetXaxis()->SetRangeUser(0,4.5);
  //leg->AddEntry(hBR[4],"Total","L");

  for (unsigned iH(0); iH<nH-1; ++iH){
    hBR[iH]->SetLineColor(iH<3?iH+2:iH+3);
    hBR[iH]->Draw("histsame");
  }
  leg->AddEntry(hBR[0],"#gamma^{D}#rightarrow e^{+}e^{-}","L");
  leg->AddEntry(hBR[1],"#gamma^{D}#rightarrow #mu^{+}#mu^{-}","L");
  leg->AddEntry(hBR[2],"#gamma^{D}#rightarrow #tau^{+}#tau^{-}","L");
  leg->AddEntry(hBR[3],"#gamma^{D}#rightarrow hadrons","L");
  leg->Draw("same");

  myc->Update();
  myc->Print("figures/BRvsmass.pdf");

  return 1;
}//main
