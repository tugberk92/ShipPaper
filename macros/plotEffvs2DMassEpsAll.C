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


int plotEffvs2DMassEpsAll(){//main

  SetTdrStyle();

  //gStyle->SetPadTopMargin(0.05);
  //gStyle->SetPadBottomMargin(0.10);
  //gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetOptStat(0);
  
  TLatex lat;
  char lbuf[500];

  const std::string prod="/home/amagnan/SOFTWARE/SHIP/ShipDPAnalysis/data/190216/";

  const unsigned nP =3;
  std::string proc[nP] ={"meson","pbrem","qcd"};
  
  const unsigned nDM = 1;//7;
  std::string decayMode[nDM] = {"all"};//"e","mu","tau","hadron"};//"pi","ka","mix","sum"};
  std::string decayModeStr[nDM] = {"pp+X"};//"ee","#mu#mu","#tau#tau","hh+X"};//"#pi#pi","KK","hh+X","pp+X"};

  int color[9] = {1,2,3,4,6,7,8,9,5};
  
  TCanvas *mycGood[nDM];
  TCanvas *mycVessel[nDM];
  TCanvas *mycReco[nDM];
  for (unsigned iDM(0); iDM<nDM; ++iDM){
    std::ostringstream lname;
    lname << "mycVessel" << iDM ;
    mycVessel[iDM] = new TCanvas(lname.str().c_str(),"",1500,500);
    mycVessel[iDM]->Divide(nP,1);
    lname.str("");
    lname << "mycReco" << iDM ;
    mycReco[iDM] = new TCanvas(lname.str().c_str(),"",1500,500);
    mycReco[iDM]->Divide(nP,1);
    lname.str("");
    lname << "mycGood" << iDM ;
    mycGood[iDM] = new TCanvas(lname.str().c_str(),"",1500,500);
    mycGood[iDM]->Divide(nP,1);
 }

  gStyle->SetOptStat(0);
  
  std::vector<Proba2>lProb[nDM][nP];
  TH2F *h2DVessel[nDM][nP];
  TH2F *h2DReco[nDM][nP];
  TH2F *h2DGood[nDM][nP];

  double maxV = 0;
  double maxR = 0;
  double maxG = 0;
  double maxDPm = 0;
  
  for (unsigned iDM(0); iDM<nDM; ++iDM){
    
    std::cout << " - Processing " << decayMode[iDM] << std::endl;
    
    for (unsigned iP(0); iP<nP; ++iP){
      
      std::cout << " - Processing " << proc[iP] << std::endl;
      std::ostringstream lname;
      lname << prod << proc[iP] << "_Ana_" << decayMode[iDM] << ".dat";
      readProba2(lname.str(),lProb[iDM][iP],0);
    
      std::sort(lProb[iDM][iP].begin(), lProb[iDM][iP].end(), customSort<Proba2>);

      std::ostringstream label;
      label << "h2DVessel_" << decayMode[iDM] << "_" << proc[iP];
      h2DVessel[iDM][iP] = new TH2F(label.str().c_str(),";log(m_{#gamma^{D}} (GeV)); log(#varepsilon); P_{vessel} (%)",iP==2?11:30,iP==2?0:-2,log10(11.),60,-9,-3);
      label.str("");
      label << "h2DReco_" << decayMode[iDM] << "_" << proc[iP];
      h2DReco[iDM][iP] = new TH2F(label.str().c_str(),";log(m_{#gamma^{D}} (GeV)); log(#varepsilon); P_{reco} (%)",iP==2?11:30,iP==2?0:-2,log10(11.),60,-9,-3);
      label.str("");
      label << "h2DGood_" << decayMode[iDM] << "_" << proc[iP];
      h2DGood[iDM][iP] = new TH2F(label.str().c_str(),";log(m_{#gamma^{D}} (GeV)); log(#varepsilon); P_{good} (%)",iP==2?11:30,iP==2?0:-2,log10(11.),60,-9,-3);

      for (unsigned iele(0); iele<lProb[iDM][iP].size(); ++iele){
	int bin = h2DVessel[iDM][iP]->FindBin(log10(lProb[iDM][iP][iele].mass),log10(lProb[iDM][iP][iele].eps));
	if (h2DVessel[iDM][iP]->GetBinContent(bin)==0 && lProb[iDM][iP][iele].Pvessel<=1. && lProb[iDM][iP][iele].Pvessel==lProb[iDM][iP][iele].Pvessel)
	  {
	    h2DVessel[iDM][iP]->Fill(log10(lProb[iDM][iP][iele].mass),log10(lProb[iDM][iP][iele].eps),lProb[iDM][iP][iele].Pvessel*100.);
	    if (maxV<lProb[iDM][iP][iele].Pvessel) {
	      maxV = lProb[iDM][iP][iele].Pvessel;
	      maxDPm = lProb[iDM][iP][iele].mass;
	    }
	  }
	if (h2DReco[iDM][iP]->GetBinContent(bin)==0 && lProb[iDM][iP][iele].Preco<=1 && lProb[iDM][iP][iele].Preco==lProb[iDM][iP][iele].Preco)
	  {
	    h2DReco[iDM][iP]->Fill(log10(lProb[iDM][iP][iele].mass),log10(lProb[iDM][iP][iele].eps),lProb[iDM][iP][iele].Preco*100.);
	    if (maxR<lProb[iDM][iP][iele].Preco) maxR = lProb[iDM][iP][iele].Preco;
	  }
	if (h2DGood[iDM][iP]->GetBinContent(bin)==0 && lProb[iDM][iP][iele].Pgood<=1 && lProb[iDM][iP][iele].Pgood==lProb[iDM][iP][iele].Pgood)
	  {
	    h2DGood[iDM][iP]->Fill(log10(lProb[iDM][iP][iele].mass),log10(lProb[iDM][iP][iele].eps),lProb[iDM][iP][iele].Pgood*100.);
	    if (maxG<lProb[iDM][iP][iele].Pgood) maxG = lProb[iDM][iP][iele].Pgood;
	  }
      }

    }
  }
  
  std::cout << " Histos filled. MaxV = " << maxV << " for mass " << maxDPm << " maxR = " << maxR << " maxG = " << maxG << std::endl;

  for (unsigned iDM(0); iDM<nDM; ++iDM){
    for (unsigned iP(0); iP<nP; ++iP){
      
      mycVessel[iDM]->cd(iP+1);
      h2DVessel[iDM][iP]->GetZaxis()->SetRangeUser(0,std::min(100.,maxV*100.));
      h2DVessel[iDM][iP]->Draw("colz");
      lat.DrawLatexNDC(0.2,0.96,(proc[iP]+" #gamma^{D}#rightarrow "+decayModeStr[iDM]).c_str());
      
      mycReco[iDM]->cd(iP+1);
      h2DReco[iDM][iP]->GetZaxis()->SetRangeUser(0,100.);
      h2DReco[iDM][iP]->Draw("colz");
      lat.DrawLatexNDC(0.2,0.96,(proc[iP]+" #gamma^{D}#rightarrow "+decayModeStr[iDM]).c_str());
      
      mycGood[iDM]->cd(iP+1);
      h2DGood[iDM][iP]->GetZaxis()->SetRangeUser(0,100.);
      h2DGood[iDM][iP]->Draw("colz");
      lat.DrawLatexNDC(0.2,0.96,(proc[iP]+" #gamma^{D}#rightarrow "+decayModeStr[iDM]).c_str());
      
    }//loop over prod

    mycVessel[iDM]->Update();
    mycVessel[iDM]->Print(("../figures/EffVesselvs2DMassEps_"+decayMode[iDM]+".pdf").c_str());

    mycReco[iDM]->Update();
    mycReco[iDM]->Print(("../figures/EffRecovs2DMassEps_"+decayMode[iDM]+".pdf").c_str());
    
    mycGood[iDM]->Update();
    mycGood[iDM]->Print(("figures/EffGoodvs2DMassEps_"+decayMode[iDM]+".pdf").c_str());
    
  }//loop over Decay Mode

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
