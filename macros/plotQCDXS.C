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

TPad* plot_ratio(TCanvas *canv, bool up){
  canv->SetFillColor      (0);
  canv->SetBorderMode     (0);
  canv->SetBorderSize     (10);
  // Set margins to reasonable defaults
  canv->SetLeftMargin     (0.17);
  canv->SetRightMargin    (0.05);
  canv->SetTopMargin      (0.05);
  canv->SetBottomMargin   (0.18);
  // Setup a frame which makes sense
  canv->SetFrameFillStyle (0);
  canv->SetFrameLineStyle (0);
  canv->SetFrameBorderMode(0);
  canv->SetFrameBorderSize(10);
  canv->SetFrameFillStyle (0);
  canv->SetFrameLineStyle (0);
  canv->SetFrameBorderMode(0);
  canv->SetFrameBorderSize(10);      

  canv->cd();
  TPad *pad = 0;
  if (up){
    pad = new TPad("upper","pad",0, 0.26 ,1 ,1);
    pad->SetBottomMargin(0.01);
    pad->SetTopMargin(0.05);
    pad->Draw();
    pad->cd();
    return pad;
  }
  else {
    pad = new TPad("lower","pad",0, 0   ,1 ,0.26);  
    pad->SetTopMargin(0.05);
    pad->SetBottomMargin(0.24);
    pad->Draw();
    return pad;
  }

};

double xsParam(double *mass, double*p){
  //if (mass[0] > 3) return exp(-5.673-0.8869*mass[0]);
  //else if (mass[0] > 1.4) return exp(-3.802-1.532*mass[0]);
  //else return 0.0586-0.09037*mass[0] + 0.0360743*mass[0]*mass[0];
  if (mass[0] > 3) return exp(-5.928-0.8669*mass[0]);
  else return exp(-4.1477-1.4745*mass[0]);
};

int plotQCDXS(){//main
  
  SetTdrStyle();
  const unsigned nVar = 10;
  const std::string suffix[13] = {"central","muFlo","muFhi","muRlo","muRhi","alphaSpdf","MSTW08","CT09MC1","CTEQ6L1","NNPDF31","NNPDF31alphaS"};

  const std::string legLabel[13] = {"NNPDF23",
				      "0.5#times #mu_{F}","2#times#mu_{F}",
				      "0.5#times #mu_{R}","2#times#mu_{R}",
				    //"#alpha_{S}=0.11","#alpha_{S}=0.15",
				      "#alpha_{S}=0.119",
				      "MSTW08","CT09MC1",
				      "CTEQ6L1",
				      "NNPDF31","#alpha_{S}=0.118"
  };
  
  const unsigned nE = 1;
  double epsilon[6] = {1e-7,8e-8,6e-8,4e-8,2e-8,1e-8};

  int color[13] = {1,2,2,3,3,1,6,7,8,9,9};
  
  TCanvas *myc = new TCanvas("myc","myc",1);
  gStyle->SetOptStat(0);
  TPad *upper = plot_ratio(myc, true);
  TPad *lower = plot_ratio(myc, false);

  TGraphErrors *grE[nE][nVar];
  TGraphErrors *grRatio[nE][nVar];
  TGraphErrors *grRatioTotP[nE];
  TGraphErrors *grRatioTotM[nE];
  TLegend *leg = new TLegend(0.83,0.27,0.94,0.94);
  leg->SetFillColor(10);
  TLegend *leg1 = new TLegend(0.7,0.6,0.94,0.94);
  leg1->SetFillColor(10);

  TLatex lat;
  char lbuf[500];
 
  for (unsigned iV(0); iV<nVar; ++iV){

    std::cout << " - Processing var " << suffix[iV] << std::endl;
    
    TFile *file = TFile::Open(("/home/amagnan/SOFTWARE/pythia8230/examples/crossSectionsQCD_"+suffix[iV]+".root").c_str());
    if (!file) continue;
    file->cd();
    
    for (unsigned iE(0); iE<nE; ++iE){//loop on epsilon

      std::cout << " --- Processing eps " << epsilon[iE] << std::endl;
      
      std::ostringstream lStr;
      lStr<<"gr_"<<iE<<"_0";
      TGraphErrors* grtmp = (TGraphErrors*)gDirectory->Get(lStr.str().c_str());
      if (!grtmp) {
	std::cout << " -- histo " << lStr.str() << " not found" << std::endl;
	continue;
      }
      grE[iE][iV] = (TGraphErrors*)grtmp->Clone((suffix[iV]+lStr.str()).c_str());

      upper->cd();
      gPad->SetLogy(1);
      
      grE[iE][iV]->SetLineColor(color[iV]);
      grE[iE][iV]->SetMarkerColor(color[iV]);
      grE[iE][iV]->SetMarkerStyle(20+iV);
      grE[iE][iV]->SetTitle(";; #sigma_{QCD}/#epsilon^{2} (mb)");
      for (unsigned iP(0); iP<grE[iE][iV]->GetN(); ++iP){
	double x0=0,y0=0;
	grE[iE][iV]->GetPoint(iP,x0,y0);
	grE[iE][iV]->SetPoint(iP,x0,y0/pow(epsilon[iE],2));
      }
      if (iV==0 || iV>5) {
	grE[iE][iV]->Draw(iV==0?"APLE":"PLEsame");
	grE[iE][iV]->GetXaxis()->SetRangeUser(1.4,6.2);
	grE[iE][iV]->GetYaxis()->SetRangeUser(0.00001,0.01);
	leg1->AddEntry(grE[iE][iV],legLabel[iV].c_str(),"P");
      }
      //lStr.str("");
      //lStr << "#epsilon = " << epsilon[iE];
      //leg->AddEntry(grE[iE][iV],lStr.str().c_str(),"P");

      if (iV>0 && iV<6) {
	//get ratio to central value
	grRatio[iE][iV] = (TGraphErrors*)grE[iE][iV]->Clone(("Ratio_"+suffix[iV]+lStr.str()).c_str());
	if (iV==1) {
	  grRatioTotP[iE] = (TGraphErrors*)grE[iE][iV]->Clone(("RatioTotPlus"+lStr.str()).c_str());
	  grRatioTotM[iE] = (TGraphErrors*)grE[iE][iV]->Clone(("RatioTotMinus"+lStr.str()).c_str());
	}
	for (unsigned iP(0); iP<grE[iE][iV]->GetN(); ++iP){
	  double x0=0,y0=0;
	  double x1=0,y1=0;
	  double x2p=0,y2p=0;
	  double x2m=0,y2m=0;
	  grE[iE][0]->GetPoint(iP,x0,y0);
	  grE[iE][iV]->GetPoint(iP,x1,y1);
	  if (iV>1) {
	    grRatioTotP[iE]->GetPoint(iP,x2p,y2p);
	    grRatioTotM[iE]->GetPoint(iP,x2m,y2m);
	  }
	  if (fabs(x1-x0)>0.001) continue;
	  grRatio[iE][iV]->SetPoint(iP,x0,y1/y0-1);
	  if (y1-y0>0) grRatioTotP[iE]->SetPoint(iP,x0,sqrt(pow((y1-y0)/y0,2)+pow(y2p,2)));
	  else grRatioTotM[iE]->SetPoint(iP,x0,-1.*sqrt(pow((y1-y0)/y0,2)+pow(y2m,2)));	  
	}
	
	lower->cd();
	gPad->SetGridy(1);
	grRatio[iE][iV]->SetTitle(";M_{#gamma^{D}} (GeV); #delta#sigma/#sigma");
	grRatio[iE][iV]->GetXaxis()->SetLabelSize(0.1);
	grRatio[iE][iV]->GetYaxis()->SetLabelSize(0.1);
	grRatio[iE][iV]->GetXaxis()->SetTitleSize(0.1);
	grRatio[iE][iV]->GetYaxis()->SetTitleSize(0.1);
	grRatio[iE][iV]->GetYaxis()->SetTitleOffset(0.5);
	grRatio[iE][iV]->Draw(iV==1?"AL":"Lsame");
	grRatio[iE][iV]->GetXaxis()->SetRangeUser(1.4,6.2);
	grRatio[iE][iV]->GetYaxis()->SetRangeUser(-0.25,0.25);
	leg->AddEntry(grRatio[iE][iV],legLabel[iV].c_str(),"L");
      }

    }
  }//loop on variations

  upper->cd();
  TF1 *myParam = new TF1("myParam",xsParam,1,10,0);
  myParam->SetLineWidth(2);
  myParam->SetLineColor(9);
  myParam->Draw("same");

  lat.DrawLatexNDC(0.17,0.85,"Pythia V8.2.30, pp#rightarrow#gamma^{D}");

  //lat.SetTextColor(9);
  //lat.DrawLatexNDC(0.72,0.34,"Fit function");
  //grE[0][0]->Fit("expo","R","",1.4,3);
  //grE[0][0]->Fit("expo","R","",3,10);
  
  leg1->AddEntry(myParam,"Fit","L");
  leg1->Draw("same");


  lower->cd();
  grRatioTotP[0]->SetLineColor(6);
  grRatioTotM[0]->SetLineColor(6);
  grRatioTotP[0]->SetLineWidth(2);
  grRatioTotM[0]->SetLineWidth(2);
 
  grRatioTotP[0]->Draw("Lsame");
  grRatioTotM[0]->Draw("Lsame");
  leg->AddEntry(grRatioTotP[0],"Total","L");
      
  leg->Draw("same");

  myc->Update();  
  myc->Print("figures/qcdXSnorm_massFit_WithSysts.pdf");
  
  return 0;

};//main
