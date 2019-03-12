#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <functional>
#include <algorithm>

#include "TCanvas.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "TLine.h"

#include "plotUtilities.C"
#include "TDRStyle.h"

int plotLimit(){
  SetTdrStyle();
  
 
  const unsigned nP = 4;
  std::string proc[nP] ={"meson","pbrem","qcd","all"}; 

  std::string syst[3] = {"","_up","_down"};
  
  TLatex lat;
  char lbuf[500];
  
  TCanvas *mycC[nP+1];
  for (unsigned iP(0); iP<nP+1; ++iP){
    std::ostringstream lname;
    lname << "mycC" << iP ;
    mycC[iP] = new TCanvas(lname.str().c_str(),"",1);
  }
 
  
  std::vector<Limit>lLimit[nP][3];
  std::vector<Limit>lLimitU[nP][3];
  std::vector<Limit>lLimitD[nP][3];
  TGraph *gr[nP][3];

  double mMax[nP] = {0,0,0,0};
  
  for (unsigned iP(0); iP<nP; ++iP){

    std::cout << " - Processing " << proc[iP] << std::endl;
    std::ostringstream lname;
    for (unsigned iL(0); iL<3; ++iL) {
      lname.str("");
      lname << "ForAMM" << proc[iP] << syst[iL] << ".txt";
      readLimit(lname.str(),lLimit[iP][iL],iP==2?1.35:0.095);

      std::sort(lLimit[iP][iL].begin(), lLimit[iP][iL].end(), customSortEps<Limit>(0));
    
      lname.str("");
      lname << "gr" << proc[iP] << syst[iL];
      gr[iP][iL] = new TGraph();
      gr[iP][iL]->SetName(lname.str().c_str());
      gr[iP][iL]->SetTitle("; m_{#gamma^{DP}} (GeV); #varepsilon");

      mMax[iP] = 0;
      unsigned turnOverPt1 = 0;
      unsigned turnOverPt2 = 0;
      for (unsigned ipt(0); ipt<lLimit[iP][iL].size(); ++ipt){
	if (lLimit[iP][iL][ipt].mass >= mMax[iP]) {
	  if (lLimit[iP][iL][ipt].mass > mMax[iP]) turnOverPt1 = ipt;
	  mMax[iP] = lLimit[iP][iL][ipt].mass;
	  turnOverPt2 = ipt;
	}
      }
      std::cout << " -- Found mMax = " << mMax[iP] << " between points " << turnOverPt1 << " and " << turnOverPt2 << std::endl;

      
      for (unsigned ipt(0); ipt<lLimit[iP][iL].size(); ++ipt){
	if (ipt<turnOverPt2) lLimitU[iP][iL].push_back(lLimit[iP][iL][ipt]);
	else lLimitD[iP][iL].push_back(lLimit[iP][iL][ipt]);
      }
     std::sort(lLimitU[iP][iL].begin(), lLimitU[iP][iL].end(), customSort<Limit>);
     std::sort(lLimitD[iP][iL].begin(), lLimitD[iP][iL].end(), customSortInv<Limit>);

      
      unsigned counter = 0;

      //bool skip[lLimit[iP][iL].size()];
      //for (unsigned ipt(0); ipt<lLimit[iP][iL].size(); ++ipt){
      //skip[ipt] = false;
      //}
      for (unsigned ipt(0); ipt<lLimitU[iP][iL].size(); ++ipt){
	/*if (ipt<turnOverPt2){
	  for (unsigned ib(ipt); ib<turnOverPt2; ++ib){
	    if (lLimit[iP][iL][ib].mass < lLimit[iP][iL][ipt].mass) skip[ib] = true;
	  }
	}
	else if (ipt>turnOverPt2){
	  for (unsigned ib(ipt); ib<lLimit[iP][iL].size(); ++ib){
	    if (lLimit[iP][iL][ib].mass > lLimit[iP][iL][ipt].mass) skip[ib] = true;
	  }
	}
	//if (skip[ipt]) continue;
	*/
	gr[iP][iL]->SetPoint(counter,lLimitU[iP][iL][ipt].mass,lLimitU[iP][iL][ipt].eps);
	counter++;
	//std::cout << ipt << " " << lLimit[iP][ipt].mass << " " << lLimit[iP][ipt].eps << std::endl;
      }
      for (unsigned ipt(0); ipt<lLimitD[iP][iL].size(); ++ipt){
	gr[iP][iL]->SetPoint(counter,lLimitD[iP][iL][ipt].mass,lLimitD[iP][iL][ipt].eps);
	counter++;
      }      
      //gr[iP][iL]->SetPoint(lLimit[iP][iL].size(),5.,1.e-9);
      
      mycC[iP]->cd();
      gPad->SetLogy(1);
      gr[iP][iL]->SetMinimum(1e-8);
      gr[iP][iL]->SetMaximum(iP<2?1e-5 : 1e-6);
      gr[iP][iL]->SetLineColor(iL+1);
      gr[iP][iL]->SetLineWidth(2);
      gr[iP][iL]->Draw(iL==0?"AL":"Lsame");
      if (iL==0) lat.DrawLatexNDC(0.1,0.95,proc[iP].c_str());
      
    }

    mycC[iP]->Update();
    mycC[iP]->Print(("figures/Limitvsmass_"+proc[iP]+".pdf").c_str());
      

  }
  
  mycC[nP]->cd();
  gPad->SetLogy(1);
  for (unsigned iP(0); iP<nP; ++iP){
    gr[iP][0]->SetLineColor(iP+1);
    gr[iP][0]->SetLineWidth(2);
    gr[iP][0]->SetMinimum(1e-8);
    gr[iP][0]->SetMaximum(3e-5);
  }
  gr[3][0]->Draw("AL");
  gr[2][0]->Draw("Lsame");
  gr[1][0]->Draw("Lsame");
  gr[0][0]->Draw("Lsame");
    
  mycC[nP]->Update();
  mycC[nP]->Print("figures/LimitvsmassSuperposed.pdf");







  return 0;
}//main
