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


int plotBRvsmass(){//main

  SetTdrStyle();

  //gStyle->SetPadTopMargin(0.05);
  //gStyle->SetPadBottomMargin(0.10);
  //gStyle->SetPadLeftMargin(0.1);
  //gStyle->SetPadRightMargin(0.05);
  gStyle->SetOptStat(0);
  
  TLatex lat;
  char lbuf[500];

  const std::string prod="/home/amagnan/SOFTWARE/SHIP/ShipDPAnalysis/data/190216/";

  const unsigned nP = 3;
  std::string proc[nP] ={"meson","pbrem","qcd"};

  const double massMin = 0.1;//0.001
  
  //const unsigned nDM = 7;
  //std::string decayMode[nDM] = {"e","mu","tau","pi","ka","mix","sum"};
  const unsigned nDM = 4;
  std::string decayMode[nDM] = {"e","mu","tau","hadron"};

  int color[9] = {1,2,3,4,6,7,8,9,5};

  const unsigned nC = 3;
  TCanvas *mycC[nC];
  for (unsigned iC(0); iC<nC; ++iC){
    std::ostringstream lname;
    lname << "myc" << iC ;
    mycC[iC] = new TCanvas(lname.str().c_str(),"",1);
  }


  TLegend *leg = new TLegend(0.7,0.65,0.94,0.94);
  leg->SetFillColor(10);
  
  std::vector<Proba>lProb[nDM];
  std::vector<Decay>lDecay;
  std::map<double,BRmean> lDecayMap[nDM];
  TGraphErrors *grDecay[nDM];

  std::map<double,std::vector<double> > lInfo;
  std::pair<std::map<double,std::vector<double> >::iterator,bool> lInsertInfo;

  
  //read number of events
  for (unsigned iP(0); iP<nP; ++iP){
      
    std::cout << " - Processing " << proc[iP] << std::endl;
    std::ostringstream lname;
    lname << prod << proc[iP] << "_Ana_sum.dat";
    readDecay(lname.str(),lDecay,massMin);
  }
  std::sort(lDecay.begin(), lDecay.end(), customSort<Decay>);
 
  
  for (unsigned iDM(0); iDM<nDM; ++iDM){
    
    std::cout << " - Processing " << decayMode[iDM] << std::endl;
    
    for (unsigned iP(0); iP<nP; ++iP){
      
      std::cout << " - Processing " << proc[iP] << std::endl;
      std::ostringstream lname;
      lname << prod << proc[iP] << "_Ana_" << decayMode[iDM] << ".dat";
      readProba(lname.str(),lProb[iDM],massMin);
    }
    
    std::sort(lProb[iDM].begin(), lProb[iDM].end(), customSort<Proba>);
    
    std::pair<std::map<double,BRmean>::iterator,bool> lInsert;

    double prevEps = 0;
    for (unsigned iele(0); iele<lProb[iDM].size(); ++iele){
      BRmean lBR;
      lBR.BR = lProb[iDM][iele].BR;
      lBR.BR2 = pow(lProb[iDM][iele].BR,2);
      lBR.nPts = 1;
      lInsert = lDecayMap[iDM].insert(std::pair<double,BRmean>(lProb[iDM][iele].mass,lBR));
      if (lInsert.second){
	prevEps = lProb[iDM][iele].eps;
      }
      else {
	if (fabs(lProb[iDM][iele].eps-prevEps)<1e-10) continue;
	(lInsert.first)->second.BR += lProb[iDM][iele].BR;
	(lInsert.first)->second.BR2 += pow(lProb[iDM][iele].BR,2);
	(lInsert.first)->second.nPts += 1;
	prevEps = lProb[iDM][iele].eps;
      }
    }
    std::ostringstream label;
    label << "hDM_" << decayMode[iDM];// << "_" << proc[iP];
    grDecay[iDM] = new TGraphErrors();
    grDecay[iDM]->SetName(label.str().c_str());
    grDecay[iDM]->SetTitle(";m_{#gamma^{D}} (GeV); BR");
    grDecay[iDM]->SetMinimum(0);
    grDecay[iDM]->SetMaximum(1);
    
    std::map<double,BRmean>::iterator lIter = lDecayMap[iDM].begin();
    unsigned iPt = 0;
    for ( ; lIter != lDecayMap[iDM].end(); ++lIter){
      grDecay[iDM]->SetPoint(iPt,lIter->first,lIter->second.calcMean());
      grDecay[iDM]->SetPointError(iPt,0,lIter->second.calcMeanError());
      iPt++;
      std::vector<double> lTmp;
      lTmp.push_back(lIter->second.calcMean());
      lInsertInfo = lInfo.insert(std::pair<double, std::vector<double> >(lIter->first,lTmp));
      if (!lInsertInfo.second){
	(lInsertInfo.first)->second.push_back(lIter->second.calcMean());
      }
    }
    
    mycC[0]->cd();
    grDecay[iDM]->SetMarkerColor(color[iDM]);
    grDecay[iDM]->SetLineColor(color[iDM]);
    grDecay[iDM]->SetMarkerStyle(20+iDM);
    grDecay[iDM]->Draw(iDM==0?"APL":"PLsame");
    
    //leg->AddEntry(grDecay[iDM],decayMode[iDM].c_str(),"P");
    
  }//loop over Decay Mode
  
  std::map<double,std::vector<double> >::iterator lIter = lInfo.begin();
  for ( ; lIter != lInfo.end(); ++lIter){
    std::cout <<lIter->first << " " ;
    double total = 0;
    for (unsigned iE(0); iE<lIter->second.size(); ++iE){
      std::cout << lIter->second[iE] << " ";
      total += lIter->second[iE];
    }
    std::cout << total << std::endl;
  }
  leg->AddEntry(grDecay[0],"#gamma^{D}#rightarrow ee","P");
  leg->AddEntry(grDecay[1],"#gamma^{D}#rightarrow #mu#mu","P");
  leg->AddEntry(grDecay[2],"#gamma^{D}#rightarrow #tau#tau","P");
  leg->AddEntry(grDecay[3],"#gamma^{D}#rightarrow hh+X ","P");

  
  leg->Draw("same");
  mycC[0]->Update();
  mycC[0]->Print("figures/DecayBRvsMass_allProd.pdf");

  mycC[1]->cd();
  //TLegend *leg2 = new TLegend(0.75,0.46,0.93,0.75);
  //leg2->SetFillColor(10);
  gPad->SetLogx(0);
  grDecay[0]->Draw("APL");
  //leg2->AddEntry(grDecay[0],"#gamma^{D}#rightarrow ee","P");
  grDecay[1]->Draw("PLsame");
  //leg2->AddEntry(grDecay[1],"#gamma^{D}#rightarrow #mu#mu","P");
  grDecay[2]->Draw("PLsame");
  //leg2->AddEntry(grDecay[2],"#gamma^{D}#rightarrow #tau#tau","P");
  grDecay[3]->Draw("PLsame");
  //leg2->AddEntry(grDecay[3],"#gamma^{D}#rightarrow hh+X ","P");

  TGraphErrors *grall = (TGraphErrors*)grDecay[0]->Clone("grDecay_all");

  int n1=0,n2=0,n3=0;
  for (unsigned iPt(0); iPt<grall->GetN(); ++iPt){
    double x1,y1,x2,y2,x3,y3,x0,y0;
    double yerr1=grDecay[1]->GetErrorY(n1);
    double yerr2=grDecay[2]->GetErrorY(n2);
    double yerr3=grDecay[3]->GetErrorY(n3);
    double yerr0=grDecay[0]->GetErrorY(iPt);
    grDecay[1]->GetPoint(n1,x1,y1);
    grDecay[2]->GetPoint(n2,x2,y2);
    grDecay[3]->GetPoint(n3,x3,y3);
    grDecay[0]->GetPoint(iPt,x0,y0);
    double sum = y0;
    double sum2 = pow(yerr0,2);
    if (fabs(x1-x0)<0.0001){
      sum += y1;
      sum2 += pow(yerr1,2);
      n1++;
    }
    if (fabs(x2-x0)<0.0001){
      sum += y2;
      sum2 += pow(yerr2,2);
      n2++;
    }
    if (fabs(x3-x0)<0.0001){
      sum += y3;
      sum2 += pow(yerr3,2);
      n3++;
    }
    grall->SetPoint(iPt,x0,sum);
    grall->SetPointError(iPt,0,sqrt(sum2));
  }
  grall->SetMarkerStyle(25);
  grall->SetMarkerColor(6);
  grall->SetLineColor(6);
  
  grall->Draw("PLsame");
  //leg->AddEntry(grall,"#gamma^{D}#rightarrow pp+X ","P");
	     
  leg->Draw("same");
  
  mycC[1]->Update();
  mycC[1]->Print("figures/DecayBRvsMass_final.pdf");

  mycC[2]->cd();
  TLegend *leg3 = new TLegend(0.15,0.13,0.37,0.45);
  leg3->SetFillColor(10);
  gPad->SetLogx(1);
  grDecay[0]->Draw("APL");
  grDecay[0]->GetXaxis()->SetRangeUser(0.15,10);
  grDecay[0]->GetYaxis()->SetRangeUser(0.0,1.0);
  grDecay[1]->Draw("PLsame");
  grDecay[2]->Draw("PLsame");
  grDecay[3]->Draw("PLsame");
  //grall->Draw("PLsame");
	     
  leg3->AddEntry(grDecay[0],"#gamma^{D}#rightarrow ee","P");
  leg3->AddEntry(grDecay[1],"#gamma^{D}#rightarrow #mu#mu","P");
  leg3->AddEntry(grDecay[2],"#gamma^{D}#rightarrow #tau#tau","P");
  leg3->AddEntry(grDecay[3],"#gamma^{D}#rightarrow hh+X ","P");
  //leg3->AddEntry(grall,"#gamma^{D}#rightarrow pp+X","P");

  leg->Draw("same");
  
  mycC[2]->Update();
  mycC[2]->Print("../figures/DecayBRvsMass_final_logX.pdf");
 
    


  return 0;
}//
