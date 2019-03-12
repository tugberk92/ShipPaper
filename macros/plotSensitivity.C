#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <map>

#include "TCanvas.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLine.h"
#include "TLegend.h"

#include "plotUtilities.C"
#include "getSystematics.C"
#include "TDRStyle.h"

void fillMap(const std::vector<Rate> & lRate,
	     std::map<std::pair<double,double>,std::pair<unsigned,double> > & lMap){
  
  std::pair<std::map<std::pair<double,double>,std::pair<unsigned,double> >::iterator,bool> lInsert;
  unsigned nTot = lRate.size();

  for (unsigned iN(0); iN<nTot; ++iN){
    lInsert = lMap.insert(std::pair<std::pair<double,double>,std::pair<unsigned,double> >(std::pair<double,double>(lRate[iN].mass,lRate[iN].eps),std::pair<unsigned,double>(1,lRate[iN].rate1)));
    if (!lInsert.second) {
      std::cout << " -- already found:  m " << lRate[iN].mass << " eps " << lRate[iN].eps << " adding " << lInsert.first->second.second << " and " << lRate[iN].rate1 << " nFound=" << lInsert.first->second.first << std::endl;
      lInsert.first->second.second += lRate[iN].rate1;
      lInsert.first->second.first += 1;
    }
  }
  std::map<std::pair<double,double>,std::pair<unsigned,double> >::iterator lIter = lMap.begin();
  
  for ( ; lIter != lMap.end(); ++lIter){
    lIter->second.second = lIter->second.second/lIter->second.first;
  }

};

void fillTotalMap(std::map<std::pair<double,double>,std::pair<unsigned,double> > & lMapOut, const std::map<std::pair<double,double>,std::pair<unsigned,double> > & lMapIn){

  std::pair<std::map<std::pair<double,double>,std::pair<unsigned,double> >::iterator,bool> lInsert;
  std::map<std::pair<double,double>,std::pair<unsigned,double> >::const_iterator lIter = lMapIn.begin();
  
  for ( ; lIter != lMapIn.end(); ++lIter){
    lInsert = lMapOut.insert(std::pair<std::pair<double,double>,std::pair<unsigned,double> >(std::pair<double,double>(lIter->first.first,lIter->first.second),std::pair<unsigned,double>(1,lIter->second.second)));
    if (!lInsert.second) {
      //std::cout << " -- already found:  m " << lIter->first.first << " eps " << lIter->first.second << " adding " << lInsert.first->second.second << " and " << lIter->second.second << " nFound=" << lInsert.first->second.first << std::endl;
      lInsert.first->second.second += lIter->second.second;
      lInsert.first->second.first += 1;
    }
    
  }

};

void fillRateFromMap(const std::map<std::pair<double,double>,std::pair<unsigned,double> > & lMap,
		     std::vector<Rate> & lRate){
  
  std::map<std::pair<double,double>,std::pair<unsigned,double> >::const_iterator lIter = lMap.begin();
  lRate.clear();
  for ( ; lIter != lMap.end(); ++lIter){
    Rate ltmp;
    ltmp.mass = lIter->first.first;
    ltmp.eps = lIter->first.second;
    ltmp.rate1 = lIter->second.second;
    bool oneProc = ((ltmp.mass>=0.9 && ltmp.mass < 1.4) || ltmp.mass > 2.2);
    unsigned nProc = lIter->second.first;
    if (oneProc ||
	(!oneProc && nProc == 2)
	){
      lRate.push_back(ltmp);
      std::cout << nProc << " m " << ltmp.mass << " eps " << ltmp.eps << " rate " << ltmp.rate1 << std::endl;
    }
  }

};

int plotSensitivity(){
  
  SetTdrStyle();
  
  //gStyle->SetPadTopMargin(0.05);
  //gStyle->SetPadBottomMargin(0.10);
  //gStyle->SetPadLeftMargin(0.11);
  //gStyle->SetPadRightMargin(0.15);
  gStyle->SetOptStat(0);

  const double clsval = 2.3;
  
  const std::string prod="/home/amagnan/SOFTWARE/SHIP/ShipDPAnalysis/data/190216/";

  TLatex lat;
  char lbuf[500];
  
  const unsigned nP = 4;
  std::string proc[nP] ={"meson","pbrem","qcd","all"}; 
  std::ofstream foutCLS[nP][3];
  
  std::vector<Rate>lRate[nP];
  TH2F *hSens[nP];

  TCanvas *mycC[nP];
  for (unsigned iP(0); iP<nP; ++iP){
    std::ostringstream lname;
    lname << "mycC" << iP ;
    mycC[iP] = new TCanvas(lname.str().c_str(),"",1);
  }

  TLegend *leg = new TLegend(0.6,0.6,0.94,0.94);
  leg->SetFillColor(10);

  
  
  std::map<std::pair<double,double>,std::pair<unsigned,double> > lMap[nP];

  for (unsigned iP(0); iP<nP; ++iP){

    std::cout << " - Processing " << proc[iP] << std::endl;
    
    std::ostringstream lname;
    lname.str("");
    lname << prod << proc[iP] << "_Ana_rate1.dat";
    if (iP<3) {
      readRate(lname.str(),lRate[iP],iP==2?1.35:0.001);
      fillMap(lRate[iP],lMap[iP]);
    }
    else{
      fillTotalMap(lMap[3],lMap[0]);
      fillTotalMap(lMap[3],lMap[1]);
      fillTotalMap(lMap[3],lMap[2]);
      fillRateFromMap(lMap[iP],lRate[iP]);
    }
    
    std::sort(lRate[iP].begin(), lRate[iP].end(), customSort<Rate>);

    lname.str("");
    lname << "ForAMM" << proc[iP];
    for (unsigned iL(0); iL<3; ++iL) {
      if (iL==0) foutCLS[iP][iL].open((lname.str()+".txt").c_str());
      else if (iL==1) foutCLS[iP][iL].open((lname.str()+"_up.txt").c_str());
      else foutCLS[iP][iL].open((lname.str()+"_down.txt").c_str());
    }
    
    //    for (unsigned iN(0); iN<lRate[iP].size(); ++iN){
    //std::cout << lRate[iP][iN].mass << " "
    //		<< lRate[iP][iN].eps << " "
    //		<< lRate[iP][iN].rate2 << std::endl;
    //}

    lname.str("");
    lname << "hSensitivity_" << proc[iP];
    hSens[iP] = new TH2F(lname.str().c_str(),";mass (GeV);#gamma' coupling to SM log_{10}(#varepsilon);rate for 2.10^{10} p.o.t.",120,0,9.2,100,-9,-3);
    if (!hSens[iP]) return 1;

    double mass = 0;
    unsigned nM = 0;
    unsigned nE = 0;
    TGraph *gr[3] = {0,0,0};
    TGraph *grTest[3] = {0,0,0};
    mycC[iP]->cd();
    gPad->SetRightMargin(0.05);
    TLine *line = 0;
    TLine *line1 = 0;
    TLine *line2 = 0;
    
    
    mycC[iP]->Print(("figures/CheckRate_"+proc[iP]+".pdf[").c_str());

    bool first = true;
    double avgSys = 0;
    for (unsigned iN(0); iN<lRate[iP].size(); ++iN){
      
      hSens[iP]->Fill(lRate[iP][iN].mass,log10(lRate[iP][iN].eps),lRate[iP][iN].rate1);
      //std::cout << lRate[iP][iN].mass << " "
      //	<< lRate[iP][iN].eps << " "
      //	<< lRate[iP][iN].rate1///(lRate[iP][iN].eps*lRate[iP][iN].eps) << " "
	//	<< (lRate[iP][iN].rate1-lRate[iP][iN].rate2)/lRate[iP][iN].rate1
	//	<< std::endl;

      double tmpRate = lRate[iP][iN].rate1;
      //if (tmpRate>50) tmpRate=50;
      //if (tmpRate<0.05) tmpRate=0.05;
      if (fabs(lRate[iP][iN].mass-mass)>0.0001 || iN==lRate[iP].size()-1) {

	if (iN==lRate[iP].size()-1){
	  gr[0]->SetPoint(nE,log10(lRate[iP][iN].eps),log10(tmpRate));
	  gr[1]->SetPoint(nE,log10(lRate[iP][iN].eps),log10(tmpRate*(1+syst(proc[iP],mass,1))));
	  gr[2]->SetPoint(nE,log10(lRate[iP][iN].eps),log10(tmpRate*(1-syst(proc[iP],mass,-1))));
	  nE += 1;
	}
	if (gr[0]) {
	  gPad->SetLogy(0);
	  gPad->SetLogx(0);
	  gPad->SetGridx(0);
	  gr[0]->SetMarkerStyle(20);
	  gr[0]->SetMarkerColor(1);
	  gr[0]->SetLineColor(1);
	  gr[0]->SetTitle(";#gamma' coupling to SM log(#varepsilon);log(rate) for 2.10^{10} p.o.t.");
	  gr[0]->Draw("AP");
	  gr[0]->GetYaxis()->SetRangeUser(-6,6);
	  sprintf(lbuf,"m_{#gamma^{D}} = %3.3f GeV",mass);
	  lat.DrawLatexNDC(0.2,0.8,lbuf);
	  line = new TLine(-9,log10(clsval),-5,log10(clsval));
	  line->SetLineColor(6);
	  line->Draw("same");

	  
	  std::cout << " doing interpolation" << std::endl;
	  //interpolate and get 2.3 event points.
	  std::vector<std::pair<double,double> > limitVec[3];
	  lname.str("");
	  lname << "hTest_" << proc[iP] << "_M" << static_cast<unsigned>(mass*1000);

	  for (unsigned iL(0); iL<3; ++iL) {
	    grTest[iL]=new TGraph();
	    if (iL==0) grTest[iL]->SetName(lname.str().c_str());
	    else if (iL==1) grTest[iL]->SetName((lname.str()+"_up").c_str());
	    else  grTest[iL]->SetName((lname.str()+"_down").c_str());
	  }
	  for (unsigned ieps(0); ieps<1000000;++ieps){
	    double tmpeps=-9.+ieps*6./1000000;
	    double limit[3];
	    double diff[3];
	    for (unsigned iL(0); iL<3; ++iL) {
	      limit[iL] = gr[iL]->Eval(tmpeps);
	      diff[iL] = fabs(limit[iL]-log10(clsval));
	      if (diff[iL]<0.01) limitVec[iL].push_back(std::pair<double,double>(diff[iL],tmpeps));
	      if (ieps%1000==0){
		grTest[iL]->SetPoint(ieps/1000,tmpeps,limit[iL]);
	      }
	    }
	    //if (ieps>700&&ieps<800&&iP==2&&mass==1.0) std::cout << "check " << mass << " " << tmpeps << " " << limit << std::endl;
	  }
	  for (unsigned iL(0); iL<3; ++iL) {
	    if (iL==0) grTest[iL]->SetLineWidth(2);
	    grTest[iL]->SetLineStyle(1+iL);
	    grTest[iL]->SetLineWidth(2);
	    grTest[iL]->SetMarkerColor(1);
	    grTest[iL]->SetLineColor(9);
	    grTest[iL]->Draw("Lsame");
	  }
	  unsigned inext[3] = {0,0,0};
	  for (unsigned iL(0); iL<3; ++iL) {
	    if (limitVec[iL].size()>0){
	      std::sort(limitVec[iL].begin(),limitVec[iL].end(),sortPair);
	      //for (unsigned iele(0); iele<limitVec.size();++iele){
	      //if (iP==1) std::cout << mass << " " << iele << " " << limitVec[iele].second << " " << limitVec[iele].first << std::endl;
	      //}
	      std::cout << "mass: " << mass << " eps " << pow(10,limitVec[iL][0].second) << " diff " << limitVec[iL][0].first << std::endl;
	      while (fabs(limitVec[iL][inext[iL]].second-limitVec[iL][0].second)<0.001 && inext[iL] < limitVec[iL].size()-1) {
		inext[iL]++;
		//	      if (inext>limitVec.size()-1) break;
	      }
	      std::cout << "Syst " << iL << " mass: " << mass << " inext " << inext[iL] << " eps " << pow(10,limitVec[iL][inext[iL]].second) << " diff " << limitVec[iL][inext[iL]].first << std::endl;
	      foutCLS[iP][iL] << mass << " " << pow(10,limitVec[iL][0].second) << std::endl;
	      foutCLS[iP][iL] << mass << " " << pow(10,limitVec[iL][inext[iL]].second) << std::endl;
	      line1 = new TLine(limitVec[iL][0].second,0,limitVec[iL][0].second,0.6);
	      line1->SetLineColor(iL==0?6:9);
	      line1->SetLineStyle(1+iL);
	      line1->Draw("same");
	      line2 = new TLine(limitVec[iL][inext[iL]].second,0,limitVec[iL][inext[iL]].second,0.6);
	      line2->SetLineColor(iL==0?6:9);
	      line2->SetLineStyle(1+iL);
	      line2->Draw("same");
	    } else {
	      std::cout << " -- did not find any point next to " << clsval << std::endl;
	    }
	  }

	  if (iP==0 && first){
	    leg->AddEntry(gr[0],"Central","P");
	    leg->AddEntry(grTest[1],"+syst","L");
	    leg->AddEntry(grTest[2],"-syst","L");
	    first = false;
	  }
	  leg->Draw("same");
	  mycC[iP]->Update();
	  mycC[iP]->Print(("figures/CheckRate_"+proc[iP]+".pdf").c_str());
	  sprintf(lbuf,"m%d",static_cast<unsigned>(mass*1000));
	  mycC[iP]->Print(("figures/CheckRate_"+proc[iP]+"_"+lbuf+".C").c_str());
	}
	nE=0;
	mass = lRate[iP][iN].mass;
	lname.str("");
	lname << "hCheckRate_" << proc[iP] << "_M" << static_cast<unsigned>(mass*1000);
	std::cout << lname.str() << std::endl;
	for (unsigned iL(0); iL<3; ++iL) {
	  gr[iL] = new TGraph();
	  //gr[0]->SetBit(TGraph::kIsSortedX);
	}
	gr[0]->SetName(lname.str().c_str());
	gr[1]->SetName((lname.str()+"_up").c_str());
	gr[2]->SetName((lname.str()+"_down").c_str());
	gr[0]->SetPoint(0,log10(lRate[iP][iN].eps),log10(tmpRate));
	gr[1]->SetPoint(0,log10(lRate[iP][iN].eps),log10(tmpRate*(1+syst(proc[iP],mass,1))));
	gr[2]->SetPoint(0,log10(lRate[iP][iN].eps),log10(tmpRate*(1-syst(proc[iP],mass,-1))));
	avgSys += syst(proc[iP],mass,1);

	nE += 1;
	nM += 1;
      }//if different mass
      else {
	gr[0]->SetPoint(nE,log10(lRate[iP][iN].eps),log10(tmpRate));
	gr[1]->SetPoint(nE,log10(lRate[iP][iN].eps),log10(tmpRate*(1+syst(proc[iP],mass,1))));
	gr[2]->SetPoint(nE,log10(lRate[iP][iN].eps),log10(tmpRate*(1-syst(proc[iP],mass,-1))));

	nE += 1;
      }
      
    }//loop on rates

    avgSys = avgSys/nM;

    std::cout << " --- Found " << nM << " mass points. Avg syst = " << avgSys << std::endl;


    mycC[iP]->Print(("figures/CheckRate_"+proc[iP]+".pdf]").c_str());
    for (unsigned iL(0); iL<3; ++iL) { 
      foutCLS[iP][iL].close();
    }
  }//loop on proc mode
  
  
  TCanvas *myc[nP];
  for (unsigned iP(0); iP<nP; ++iP){
    std::ostringstream lname;
    lname << "myc" << iP ;
    myc[iP] = new TCanvas(lname.str().c_str(),"",1);
    myc[iP]->SetRightMargin(0.1);
  }
  
  for (unsigned iP(0); iP<nP; ++iP){
    myc[iP]->cd();
    gPad->SetLogx(1);
    gPad->SetLogy(0);
    gPad->SetLogz(1);
    if (iP==0) hSens[iP]->GetXaxis()->SetRangeUser(0,9.1);
    else if (iP==1) hSens[iP]->GetXaxis()->SetRangeUser(0,9.1);
    else if (iP==2) hSens[iP]->GetXaxis()->SetRangeUser(0,9.1);
    hSens[iP]->GetZaxis()->SetRangeUser(0.01,100);
    hSens[iP]->Draw("colz");
    lat.DrawLatexNDC(0.1,1.0,proc[iP].c_str());
    myc[iP]->Update();
    myc[iP]->Print(("figures/Sensitivity_"+proc[iP]+".pdf").c_str());    
  }
  
  
    
  return 0;
}
