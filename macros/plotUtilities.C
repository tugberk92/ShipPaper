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

struct Rate{
  double mass;
  double eps;
  double rate1;
  //double rate2;
  //double vtxEff;
};

struct Limit{
  double mass;
  double eps;
  
};

//bool operator<(const Limit & left, const Limit & right){
//return left.mass < right.mass;
//};

struct Proba{
  double mass;
  double eps;
  double BR;
  double Pvessel;
  double Preco;
};

struct Proba2{
  double mass;
  double eps;
  double Pgood;
  double BR;
  double Pvessel;
  double Preco;
};

struct Decay{
  double mass;
  double eps;
  unsigned Ngen;
  unsigned NDP;
  unsigned Ngood;
  unsigned Nvessel;
  unsigned Nreco;
};

struct BRmean{
  double BR;
  double BR2;
  unsigned nPts;
  double calcMean(){
    if (nPts>0) return BR/nPts;
    else return 0;
  }
  double calcSigma(){
    if (nPts>1) return sqrt(BR2/(nPts-1)-pow(calcMean(),2));
    else return 0;
  }
  double calcMeanError(){
    if (nPts>0) return calcSigma()/sqrt(nPts);
    else return 0;
  }
};

template <class T>
static bool customSort(const T & a,
		       const T & b)
{
  if (fabs(a.mass-b.mass)<0.0001) return a.eps >= b.eps;
  else return a.mass <= b.mass;
};

template <class T>
static bool customSortInv(const T & a,
			  const T & b)
{
  if (fabs(a.mass-b.mass)<0.0001) return a.eps >= b.eps;
  else return a.mass >= b.mass;
};

template <class T>
class customSortEps{

public:
  customSortEps(const double & mMax=0) : mMax_(mMax){}
  bool operator()(const T & a,const T & b){
    //if (a.mass<=mMax_) return a.mass>b.mass;
    //else return a.mass < b.mass;
    
    return a.eps >= b.eps;// && ((a.mass<=mMax_ && a.mass>=b.mass) || (a.mass>mMax_ && a.mass<=b.mass));
  };
  
private:
  double mMax_;
  
};

static bool sortPair(const std::pair<double,double> & a,
		     const std::pair<double,double> & b)
  
{
  return a.first <= b.first;
};

void readRate(const std::string & prodfile,
	      std::vector<Rate> & ratevec,
	       const double & massMin = 0){

  std::ifstream lInput;
  lInput.open(prodfile.c_str());
  if(!lInput.is_open()) {
    std::cerr << "Unable to open file: " << prodfile << ". Exit..." << std::endl;
    exit(1);
  }
  
  while(1){
    Rate tmp;
    tmp.mass = 0;
    tmp.eps = 0;
    tmp.rate1 = 0;
    //tmp.rate2 = 0;
    //tmp.vtxEff = 1;
    lInput>>tmp.mass>>tmp.eps>>tmp.rate1;//>>tmp.rate2>>tmp.vtxEff;

    //protect against nan
    if (tmp.rate1!=tmp.rate1) tmp.rate1=0;
    //if (tmp.rate2!=tmp.rate2) tmp.rate2=0;

    //protect against 0 values
    if (tmp.mass>massMin){// && tmp.rate2>0.00000){
      ratevec.push_back(tmp);
    } else if (tmp.mass==0) {
      std::cout << " -- 0's for " << tmp.mass << " " << tmp.eps << " " << tmp.rate1
	//<< " " << tmp.rate2
		<< std::endl;
    }
    if(lInput.eof()){
      break;
    }
  }

  std::cout << " Found " << ratevec.size() << " points for file " << prodfile << std::endl;
  
  lInput.close();
};

void readLimit(const std::string & prodfile,
	       std::vector<Limit> & ratevec,
	       const double & massMin = 0){

  std::ifstream lInput;
  lInput.open(prodfile.c_str());
  if(!lInput.is_open()) {
    std::cerr << "Unable to open file: " << prodfile << ". Exit..." << std::endl;
    exit(1);
  }
  
  while(1){
    Limit tmp;
    tmp.mass = 0;
    tmp.eps = 0;
    lInput>>tmp.mass>>tmp.eps;

    if (tmp.mass>massMin){// && tmp.rate2>0.00000){
      ratevec.push_back(tmp);
    } else if (tmp.mass==0){
      std::cout << " -- 0's for " << tmp.mass << " " << tmp.eps
		<< std::endl;
    }
    if(lInput.eof()){
      break;
    }
  }

  std::cout << " Found " << ratevec.size() << " points for file " << prodfile << std::endl;
  
  lInput.close();
};

void readProba(const std::string & prodfile,
	       std::vector<Proba> & ratevec,
	       const double minMass = 0){

  std::ifstream lInput;
  lInput.open(prodfile.c_str());
  if(!lInput.is_open()) {
    std::cerr << "Unable to open file: " << prodfile << ". Exit..." << std::endl;
    exit(1);
  }
  
  while(1){
    Proba tmp;
    tmp.mass = 0;
    tmp.eps = 0;
    tmp.BR = 0;
    tmp.Pvessel = 0;
    tmp.Preco = 0;
    lInput>>tmp.mass>>tmp.eps>>tmp.BR>>tmp.Pvessel>>tmp.Preco;
    //protect against 0 values
    if (tmp.mass>minMass && tmp.BR>0.0){
      ratevec.push_back(tmp);
    } else if (tmp.mass==0){
      std::cout << " -- 0's for " << tmp.mass << " " << tmp.eps << " " << tmp.BR << std::endl;
    }
    if(lInput.eof()){
      break;
    }
  }

  std::cout << " Found " << ratevec.size() << " points for file " << prodfile << std::endl;
  
  lInput.close();
};


void readProba2(const std::string & prodfile,
	       std::vector<Proba2> & ratevec,
	       const double minMass = 0){

  std::ifstream lInput;
  lInput.open(prodfile.c_str());
  if(!lInput.is_open()) {
    std::cerr << "Unable to open file: " << prodfile << ". Exit..." << std::endl;
    exit(1);
  }
  
  while(1){
    Proba2 tmp;
    tmp.mass = 0;
    tmp.eps = 0;
    tmp.Pgood = 0;
    tmp.BR = 0;
    tmp.Pvessel = 0;
    tmp.Preco = 0;
    lInput>>tmp.mass>>tmp.eps>>tmp.Pgood>>tmp.BR>>tmp.Pvessel>>tmp.Preco;
    //protect against 0 values
    if (tmp.mass>minMass && tmp.BR>0.0){
      ratevec.push_back(tmp);
    } else if (tmp.mass==0){
      std::cout << " -- 0's for " << tmp.mass << " " << tmp.eps << " " << tmp.BR << std::endl;
    }
    if(lInput.eof()){
      break;
    }
  }

  std::cout << " Found " << ratevec.size() << " points for file " << prodfile << std::endl;
  
  lInput.close();
};


void readDecay(const std::string & prodfile,
	       std::vector<Decay> & ratevec,
	       const double minMass = 0){

  std::ifstream lInput;
  lInput.open(prodfile.c_str());
  if(!lInput.is_open()) {
    std::cerr << "Unable to open file: " << prodfile << ". Exit..." << std::endl;
    exit(1);
  }
  
  while(1){
    Decay tmp;
    tmp.mass = 0;
    tmp.eps = 0;
    tmp.Ngen = 0;
    tmp.NDP = 0;
    tmp.Ngood = 0;
    tmp.Nvessel = 0;
    tmp.Nreco = 0;

    lInput>>tmp.mass>>tmp.eps>>tmp.Ngen>>tmp.NDP>>tmp.Ngood>>tmp.Nvessel>>tmp.Nreco;
    //protect against 0 values
    if (tmp.mass>minMass && tmp.Ngood>0){
      ratevec.push_back(tmp);
    } else if (tmp.mass==0){
      std::cout << " -- 0's for " << tmp.mass << " " << tmp.eps << " " << tmp.Ngood << std::endl;
    }
    if(lInput.eof()){
      break;
    }
  }

  std::cout << " Found " << ratevec.size() << " points for file " << prodfile << std::endl;
  
  lInput.close();
};

