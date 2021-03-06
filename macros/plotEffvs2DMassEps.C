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
#include "TGraphErrors.h"
#include "TLegend.h"


#include "plotUtilities.C"
#include "TDRStyle.h"


int plotEffvs2DMassEps(){//main

  SetTdrStyle();

  //gStyle->SetPadTopMargin(0.05);
  //gStyle->SetPadBottomMargin(0.10);
  //gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetOptStat(0);
  
  TLatex lat;
  char lbuf[500];

  const std::string prod="/home/amagnan/SOFTWARE/SHIP/ShipDPAnalysis/data/190216/";

  const unsigned nP = 3;
  std::string proc[nP] ={"meson","pbrem","qcd"};
  
  const unsigned nDM = 4;//7;
  std::string decayMode[nDM] = {"e","mu","tau","hadron"};//"pi","ka","mix","sum"};
  std::string decayModeStr[nDM] = {"ee","#mu#mu","#tau#tau","hh+X"};//"#pi#pi","KK","hh+X","pp+X"};

  int color[9] = {1,2,3,4,6,7,8,9,5};
  
  TCanvas *mycVessel[nP];
  TCanvas *mycReco[nP];
  for (unsigned iP(0); iP<nP; ++iP){
    std::ostringstream lname;
    lname << "mycVessel" << iP ;
    mycVessel[iP] = new TCanvas(lname.str().c_str(),"",1000,1000);
    mycVessel[iP]->Divide(nDM/2,2);
    lname.str("");
    lname << "mycReco" << iP ;
    mycReco[iP] = new TCanvas(lname.str().c_str(),"",1000,1000);
    mycReco[iP]->Divide(nDM/2,2);
  }

  gStyle->SetOptStat(0);
  
  std::vector<Proba>lProb[nDM][nP];
  TH2F *h2DVessel[nDM][nP];
  TH2F *h2DReco[nDM][nP];

  double maxV = 0;
  double maxR = 0;
  
  for (unsigned iDM(0); iDM<nDM; ++iDM){
    
    std::cout << " - Processing " << decayMode[iDM] << std::endl;
    
    for (unsigned iP(0); iP<nP; ++iP){
      
      std::cout << " - Processing " << proc[iP] << std::endl;
      std::ostringstream lname;
      lname << prod << proc[iP] << "_Ana_" << decayMode[iDM] << ".dat";
      readProba(lname.str(),lProb[iDM][iP],0);
    
      std::sort(lProb[iDM][iP].begin(), lProb[iDM][iP].end(), customSort<Proba>);

      std::ostringstream label;
      label << "h2DVessel_" << decayMode[iDM] << "_" << proc[iP];
      h2DVessel[iDM][iP] = new TH2F(label.str().c_str(),";log(m_{#gamma^{D}}) (log(GeV)); log(#varepsilon); P_{vessel} (%)",iP==2?11:30,iP==2?0:-2,log10(11.),60,-9,-3);
      label.str("");
      label << "h2DReco_" << decayMode[iDM] << "_" << proc[iP];
      h2DReco[iDM][iP] = new TH2F(label.str().c_str(),";log(m_{#gamma^{D}}) (log(GeV)); log(#varepsilon); P_{reco} (%)",iP==2?11:30,iP==2?0:-2,log10(11.),60,-9,-3);

      for (unsigned iele(0); iele<lProb[iDM][iP].size(); ++iele){
	int bin = h2DVessel[iDM][iP]->FindBin(log10(lProb[iDM][iP][iele].mass),log10(lProb[iDM][iP][iele].eps));
	if (h2DVessel[iDM][iP]->GetBinContent(bin)==0 && lProb[iDM][iP][iele].Pvessel<=1. && lProb[iDM][iP][iele].Pvessel==lProb[iDM][iP][iele].Pvessel)
	  {
	    h2DVessel[iDM][iP]->Fill(log10(lProb[iDM][iP][iele].mass),log10(lProb[iDM][iP][iele].eps),lProb[iDM][iP][iele].Pvessel*100.);
	    if (maxV<lProb[iDM][iP][iele].Pvessel) maxV = lProb[iDM][iP][iele].Pvessel;
	  }
	if (h2DReco[iDM][iP]->GetBinContent(bin)==0 && lProb[iDM][iP][iele].Preco<=1 && lProb[iDM][iP][iele].Preco==lProb[iDM][iP][iele].Preco)
	  {
	    h2DReco[iDM][iP]->Fill(log10(lProb[iDM][iP][iele].mass),log10(lProb[iDM][iP][iele].eps),lProb[iDM][iP][iele].Preco*100.);
	    if (maxR<lProb[iDM][iP][iele].Preco) maxR = lProb[iDM][iP][iele].Preco;
	  }
      }

    }
  }
  
  std::cout << " Histos filled. MaxV = " << maxV << " maxR = " << maxR << std::endl;

  for (unsigned iP(0); iP<nP; ++iP){
    for (unsigned iDM(0); iDM<nDM; ++iDM){
      
      mycVessel[iP]->cd(iDM+1);
      h2DVessel[iDM][iP]->GetZaxis()->SetRangeUser(0,std::min(100.,maxV*100.));
      h2DVessel[iDM][iP]->Draw("colz");
      lat.DrawLatexNDC(0.2,0.96,(proc[iP]+" #gamma^{D}#rightarrow "+decayModeStr[iDM]).c_str());
      
      mycReco[iP]->cd(iDM+1);
      h2DReco[iDM][iP]->GetZaxis()->SetRangeUser(0,100.);
      h2DReco[iDM][iP]->Draw("colz");
      lat.DrawLatexNDC(0.2,0.96,(proc[iP]+" #gamma^{D}#rightarrow "+decayModeStr[iDM]).c_str());
      
    }//loop over Decay Mode
    
    mycVessel[iP]->Update();
    mycVessel[iP]->Print(("../figures/EffVesselvs2DMassEps_"+proc[iP]+".pdf").c_str());

    mycReco[iP]->Update();
    mycReco[iP]->Print(("../figures/EffRecovs2DMassEps_"+proc[iP]+".pdf").c_str());
    
    }//loop over prod

  return 1;
  
  const unsigned neps = 5;
  TH1F *h1DVessel[nDM][nP][neps];
  TH1F *h1DReco[nDM][nP][neps];

  for (unsigned iDM(0); iDM<nDM; ++iDM){
    for (unsigned iP(0); iP<nP; ++iP){
 
      unsigned ieps = 0;
      for (unsigned iY(0); iY<h2DVessel[iDM][iP]->GetNbinsY()+2 && ieps < neps; ++iY){

	if (iY==h2DVessel[iDM][iP]->GetYaxis()->FindBin(-8+0.5*ieps+0.01)){
	  std::ostringstream label;
	  label << "h1DVessel_" << decayMode[iDM] << "_" << proc[iP] << "_" << ieps;
	  h1DVessel[iDM][iP][ieps] = new TH1F(label.str().c_str(),";log(m_{#gamma^{D}}) (log(GeV));P_{vessel} (%)",iP==2?10:30,iP==2?0:-2,1);
	  label.str("");
	  label << "h1DReco_" << decayMode[iDM] << "_" << proc[iP] << "_" << ieps;
	  h1DReco[iDM][iP][ieps] = new TH1F(label.str().c_str(),";log(m_{#gamma^{D}}) (log(GeV));P_{reco} (%)",iP==2?10:30,iP==2?0:-2,1);
	  for (unsigned iX(0); iX<h2DVessel[iDM][iP]->GetNbinsX()+2; ++iX){
	    h1DVessel[iDM][iP][ieps]->Fill(h2DVessel[iDM][iP]->GetXaxis()->GetBinCenter(iX),h2DVessel[iDM][iP]->GetBinContent(iX,iY));
	    h1DReco[iDM][iP][ieps]->Fill(h2DReco[iDM][iP]->GetXaxis()->GetBinCenter(iX),h2DReco[iDM][iP]->GetBinContent(iX,iY));
	  }
	  ieps++;
	}

      }

      mycVessel[iDM]->cd(iP+1);
      for (ieps=0; ieps<neps; ++ieps){
	h1DVessel[iDM][iP][ieps]->GetYaxis()->SetRangeUser(0,std::min(100.,maxV*100.));
	h1DVessel[iDM][iP][ieps]->SetMarkerColor(color[ieps]);
	h1DVessel[iDM][iP][ieps]->SetLineColor(color[ieps]);
	h1DVessel[iDM][iP][ieps]->SetMarkerStyle(20+ieps);
	h1DVessel[iDM][iP][ieps]->Draw(ieps==0?"PE":"PEsame");
      }
      lat.DrawLatexNDC(0.2,0.96,(proc[iP]+" #gamma^{D}#rightarrow "+decayModeStr[iDM]).c_str());
      
      mycReco[iDM]->cd(iP+1);
      for (ieps=0; ieps<neps; ++ieps){
	h1DReco[iDM][iP][ieps]->GetYaxis()->SetRangeUser(0,100.);
	h1DReco[iDM][iP][ieps]->SetMarkerColor(color[ieps]);
	h1DReco[iDM][iP][ieps]->SetLineColor(color[ieps]);
	h1DReco[iDM][iP][ieps]->SetMarkerStyle(20+ieps);
	h1DReco[iDM][iP][ieps]->Draw(ieps==0?"PE":"PEsame");
      }
      lat.DrawLatexNDC(0.2,0.96,(proc[iP]+" #gamma^{D}#rightarrow "+decayModeStr[iDM]).c_str());
      
    }//loop over prod

    mycVessel[iDM]->Update();
    mycVessel[iDM]->Print(("../figures/EffVesselvs1DMass_"+decayMode[iDM]+".pdf").c_str());

    mycReco[iDM]->Update();
    mycReco[iDM]->Print(("../figures/EffRecovs1DMass_"+decayMode[iDM]+".pdf").c_str());
    
  }
  
  return 0;
}//
