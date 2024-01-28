// prttools - useful functions 
// original author: Roman Dzhygadlo - GSI Darmstadt 

#ifndef prttools_h
#define prttools_h 1

#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TSpline.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TMath.h"
#include "TChain.h"
#include "TGaxis.h"
#include "TColor.h"
#include "TString.h"
#include "TArrayD.h"
//#include "TSpectrum.h"
//#include "TSpectrum2.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include "TPaveStats.h"
#include "TObjString.h"
#include "TApplication.h"
#include <TLegend.h>
#include <TAxis.h>
#include <TPaletteAxis.h>
#include <TRandom.h>
#include <TCutG.h>
#include <TKey.h>
#include "TPRegexp.h"
#include "TFitResult.h"

#include <iostream>
#include <fstream>
#include <sstream>

const int prt_nmcp = 24;
const int prt_npix = 16 * 16; 

const int prt_maxdircch = prt_nmcp*prt_npix;
const int prt_maxnametdc=10000;

TRandom  prt_rand;
TChain*  prt_ch = 0;
int    prt_entries(0),prt_particle(0),prt_geometry(2023),prt_beamx(0),prt_beamz(0);
int    prt_last_layoutId,prt_last_max,prt_last_maxz,prt_last_minz;
double prt_theta(0), prt_test1(0),prt_test2(0),prt_mom(0),prt_phi(0);
TString  prt_savepath(""),prt_info("");
TH2F*    prt_hdigi[prt_nmcp];
//TSpectrum *prt_spect = new TSpectrum(2);

const int prt_ntdcm = 80; //41
int prt_ntdc = prt_ntdcm;
int prt_maxch = prt_ntdc*48;
TString prt_tdcsid[prt_ntdcm];
const int prt_maxchm = prt_ntdcm*48;
int map_tdc[prt_maxnametdc];
int map_mpc[prt_maxchm/prt_npix][prt_npix];
int map_mcp[prt_maxchm];
int map_pix[prt_maxchm];
int map_row[prt_maxchm];
int map_col[prt_maxchm];

int prt_pid(0), prt_pdg[]={-11,13,211,321,2212};
double prt_mass[]={0.000511,0.1056584,0.139570,0.49368,0.9382723};
TString  prt_name[]={"e","muon","pion","kaon","proton"};
TString  prt_lname[]={"e","#mu","#pi","K","p"};
int    prt_color[]={kOrange+6,kCyan+1,kBlue+1,kRed+1,kRed+1};
double prt_particleArray[3000];

TF1 *prt_gaust;
/*TVector3 prt_fit(TH1 *h, double range = 3, double threshold=20, double limit=2, int peakSearch=1,int bkg = 0, TString opt="MQ"){
  int binmax = h->GetMaximumBin();
  double xmax = h->GetXaxis()->GetBinCenter(binmax);
  if(bkg==0) prt_gaust = new TF1("prt_gaust","[0]*exp(-0.5*((x-[1])/[2])^2)",xmax-range,xmax+range);
  else if(bkg==1) prt_gaust = new TF1("prt_gaust","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+x*[4]",xmax-range,xmax+range);
  prt_gaust->SetNpx(500);
  prt_gaust->SetParNames("const","mean","sigma");
  prt_gaust->SetLineColor(2);
  double integral = h->Integral(h->GetXaxis()->FindBin(xmax-range),h->GetXaxis()->FindBin(xmax+range));
  double xxmin, xxmax, sigma1(0), mean1(0), mean2(0);
  xxmax = xmax;
  xxmin = xxmax;
  int nfound(1);
  if(integral>threshold){
    
    if(peakSearch == 1){
      prt_gaust->SetParameter(1,xmax);
      prt_gaust->SetParameter(2,0.005);
      prt_gaust->SetParLimits(2,0.003,limit);
      h->Fit("prt_gaust",opt,"",xxmin-range, xxmax+range);
    }
    
    if(peakSearch>1){
      nfound = prt_spect->Search(h,4,"goff",0.1);
      std::cout<<"nfound  "<<nfound <<std::endl;
      if(nfound==1){
	prt_gaust =new TF1("prt_gaust","gaus(0)",xmax-range,xmax+range);
	prt_gaust->SetNpx(500);
	prt_gaust->SetParameter(1,prt_spect->GetPositionX()[0]);
      }else if(nfound>=2){
	double p1 = prt_spect->GetPositionX()[0];
	double p2 = prt_spect->GetPositionX()[1];
	if(p1>p2) {
	  xxmax = p1;
	  xxmin = p2;
	}else {
	  xxmax = p2;
	  xxmin = p1;
	}
	if(peakSearch==20){
	  xxmax=xxmin;
	  prt_gaust =new TF1("prt_gaust","gaus(0)",xxmin-range,xxmin+range);
	  prt_gaust->SetNpx(500);
	  prt_gaust->SetParameter(1,prt_spect->GetPositionX()[0]);
	}else{
	  prt_gaust =new TF1("prt_gaust","gaus(0)+gaus(3)",xmax-range,xmax+range);
	  prt_gaust->SetNpx(500);
	  prt_gaust->SetParameter(0,1000);
	  prt_gaust->SetParameter(3,1000);
	
	  prt_gaust->FixParameter(1,xxmin);
	  prt_gaust->FixParameter(4,xxmax);
	  prt_gaust->SetParameter(2,0.1);
	  prt_gaust->SetParameter(5,0.1);
	  h->Fit("prt_gaust",opt,"",xxmin-range, xxmax+range);
	  prt_gaust->ReleaseParameter(1);
	  prt_gaust->ReleaseParameter(4);
	}
      }
    
      prt_gaust->SetParameter(2,0.2);
      prt_gaust->SetParameter(5,0.2);
    }

    h->Fit("prt_gaust",opt,"",xxmin-range, xxmax+range);
    mean1 = prt_gaust->GetParameter(1);
    sigma1 = prt_gaust->GetParameter(2);
    if(sigma1>10) sigma1=10;
    
    if(peakSearch == 2){
      mean2 = (nfound==1) ? prt_gaust->GetParameter(1) : prt_gaust->GetParameter(4);
      // sigma2 = (nfound==1) ? prt_gaust->GetParameter(2) : prt_gaust->GetParameter(5);
    }
  }
  delete prt_gaust;
  return TVector3(mean1,sigma1,mean2);
  }

TGraph *prt_fitslices(TH2F *hh,double minrange=0, double maxrange=0, double fitrange=1,int rebin=1,int ret=0){
  TH2F *h =(TH2F*) hh->Clone("h");
  h->RebinY(rebin);
  int point(0);
  TGraph *gres = new TGraph();
  for (int i=1;i<h->GetNbinsY();i++){
    double x = h->GetYaxis()->GetBinCenter(i);
    TH1D* hp;
    if(minrange!=maxrange){
      TCutG *cut = new TCutG("prt_onepeakcut",5);
      cut->SetVarX("y");
      cut->SetVarY("x");
      cut->SetPoint(0,minrange,-1E6);
      cut->SetPoint(1,minrange, 1E6);
      cut->SetPoint(2,maxrange, 1E6);
      cut->SetPoint(3,maxrange,-1E6);
      cut->SetPoint(4,minrange,-1E6);
    
      hp = h->ProjectionX(Form("bin%d",i),i,i,"[prt_onepeakcut]");
    }else{
      hp = h->ProjectionX(Form("bin%d",i),i,i);
    }

    TVector3 res = prt_fit((TH1F*)hp,fitrange,100,2,1,1);
    double y=0;
    if(ret==0) y = res.X();
    if(ret==1) y = res.Y();
    if(ret==2) y = res.X() + 0.5*res.Y();
    if(ret==3) y = res.X() - 0.5*res.Y();
    if(y==0 || y<minrange || y>maxrange) continue;
    
    gres->SetPoint(point,y,x);
    gres->SetLineWidth(2);
    gres->SetLineColor(kRed);
    point++;
  }
  return gres;
  }*/

void prt_createMap(int setupid = 2019){
  prt_geometry = setupid;

  TGaxis::SetMaxDigits(4);
  for(int i=0; i<prt_maxnametdc; i++) map_tdc[i]=-1;
  for(int i=0; i<prt_ntdc; i++){
    int dec = TString::BaseConvert(prt_tdcsid[i],16,10).Atoi();
    map_tdc[dec]=i;
  }
  
  //  for(int ch=0; ch<prt_maxdircch; ch++){
  for(int ch=0; ch<prt_maxch; ch++){
    int mcp = ch/prt_npix;
    int pix = ch%prt_npix;
    int col = pix/2 - 8*(pix/16);
    int row = pix%2 + 2*(pix/16);
    // pix = col+sqrt(prt_npix)*row;

    map_mpc[mcp][pix]=ch;
    map_mcp[ch] = mcp;
    map_pix[ch] = pix;
    map_row[ch] = row;
    map_col[ch] = col;
  } 

  for(int i=0; i<5; i++){
    prt_particleArray[prt_pdg[i]]=i;
  }
  prt_particleArray[212]=2;
}

TString prt_randstr(int len = 10){
  TString str = ""; 
  static const char alphanum[] =
    "0123456789"
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz";

  for (int i = 0; i < len; ++i) {
    str += alphanum[rand() % (sizeof(alphanum) - 1)];
  }
  return str;
}

// layoutId == 5    - 5 row's design for the PANDA Barrel DIRC
// layoutId == 2015 - cern 2015
// layoutId == 2016 - cern 2016
// layoutId == 2017 - cern 2017
// layoutId == 2018 - cern 2018
// layoutId == 2021 - new 3.6 row's design for the PANDA Barrel DIRC
// layoutId == 2023 - new 2x4 layout for the PANDA Barrel DIRC
// layoutId == 2031 - EIC DIRC beam test
// layoutId == 2030 - EIC DIRC prism
// layoutId == 2032 - EIC DIRC focusing prism 
TH1 * prt_cdigi_th;
TCanvas *prt_drawDigi(int layoutId = 0, double maxz = 0, double minz = 0, TCanvas *cdigi = NULL){

  prt_last_layoutId=layoutId;
  prt_last_maxz=maxz;
  prt_last_minz=minz;

  if(prt_geometry==2021) layoutId=2021;
  // if(prt_geometry==2019) layoutId=2018;

  TString sid = prt_randstr(3);
  if(cdigi) cdigi->cd();
  else cdigi = new TCanvas("hp="+sid,"hp_"+sid,800,400);

  TPad* prt_hpads[prt_nmcp];
  TPad* prt_hpglobal;
  
  if(layoutId==2015 ||  layoutId==5) prt_hpglobal = new TPad(sid,"T",0.04,0.04,0.88,0.96);
  else if(layoutId==2021) prt_hpglobal = new TPad(sid,"T",0.12,0.02,0.78,0.98);
  else if(layoutId==2016) prt_hpglobal = new TPad(sid,"T",0.2,0.02,0.75,0.98);
  else if(layoutId==2017) prt_hpglobal = new TPad(sid,"T",0.15,0.02,0.80,0.98);
  else if(layoutId==2018) prt_hpglobal = new TPad(sid,"T",0.05,0.07,0.9,0.93);
  else if(layoutId==2023) prt_hpglobal = new TPad(sid,"T",0.073,0.02,0.877,0.98);
  else if(layoutId==2030) prt_hpglobal = new TPad(sid,"T",0.10,0.025,0.82,0.975);
  else if(layoutId==2032) prt_hpglobal = new TPad(sid,"T",0.04,0.025,0.91,0.975);
  else if(layoutId==2031) prt_hpglobal = new TPad(sid,"T",0.12,0.01,0.80,0.99);
  else  prt_hpglobal = new TPad(sid,"T",0.04,0.04,0.96,0.96);

  prt_hpglobal->SetFillStyle(0);
  prt_hpglobal->Draw();  
  prt_hpglobal->cd();
  
  int nrow = 3, ncol = 5;
  if(layoutId ==2016) ncol=3;
  if(layoutId ==2017) ncol=4;
  if(layoutId ==2018 || layoutId ==2023) {nrow=2; ncol=4;}
  if(layoutId ==2021) ncol=4;
  if(layoutId ==2030) {nrow=4; ncol=6;}
  if(layoutId ==2032) {nrow=4; ncol=7;}
  if(layoutId ==2031) {nrow=3; ncol=4;}
  
  if(layoutId > 1){
    float tbw(0.02), tbh(0.01), shift(0),shiftw(0.02),shifth(0),margin(0.01);
    int padi(0);
      
    for(int i=0; i<ncol; i++){
      for(int j=0; j<nrow; j++){
	if(j==1) shift = -0.028;
	else shift = 0;
	shifth=0;
	if(layoutId == 5) {shift =0; shiftw=0.001; tbw=0.001; tbh=0.001;}
	if(layoutId == 2021) {
	  if(i==0 && j == nrow-1) continue;
	  shift =0; shiftw=0.001; tbw=0.001; tbh=0.001;
	  if(i==0) shifth=0.167;
	}
	if(layoutId == 2016) {
	  shift = -0.01; shiftw=0.01; tbw=0.03; tbh=0.006;
	  if(j==1) shift += 0.015;
	}
	if(layoutId == 2017) {
	  margin= 0.1;
	  shift = 0; shiftw=0.01; tbw=0.005; tbh=0.006;
	  //if(j==1) shift += 0.015;
	}
	if(layoutId == 2018) {
	  margin= 0.1;
	  shift = 0; shiftw=0.01; tbw=0.005; tbh=0.006;
	}
	if(layoutId == 2023) {
	  margin= 0.1;
	  shift = 0; shiftw=0.01; tbw=0.0015; tbh=0.042;
	}
	if(layoutId == 2030) {
	  margin= 0.1;
	  shift = 0; shiftw=0.01; tbw=0.001; tbh=0.001;
	  padi=j*ncol+i;
	}
	if(layoutId == 2032) {
	  margin= 0.1;
	  shift = 0; shiftw=0.01; tbw=0.001; tbh=0.001;
	  padi=j*ncol+i;
	}
	if(layoutId == 2031) {
	  margin= 0.1;
	  shift = 0; shiftw=0.01; tbw=0.001; tbh=0.001;
	  padi=i*nrow+j;
	}

	prt_hpads[padi] =  new TPad(sid+Form("P%d",i*10+j),"T",
				    i/(ncol+2*margin)+tbw+shift+shiftw,
				    j/(double)nrow+tbh+shifth,
				    (i+1)/(ncol+2*margin)-tbw+shift+shiftw,
				    (1+j)/(double)nrow-tbh+shifth, 21);
	prt_hpads[padi]->SetFillColor(kCyan-8);
	prt_hpads[padi]->SetMargin(0.055,0.055,0.055,0.055);
	prt_hpads[padi]->Draw();
	padi++;
      }
    }
  }else{
    float tbw(0.02), tbh(0.01), shift(0),shiftw(-0.02);
    int padi(0);
    for(int ii=0; ii<ncol; ii++){
      for(int j=0; j<nrow; j++){
	if(j==1) shift = 0.04;
	else shift = 0;
	prt_hpads[padi] =  new TPad(Form("P%d",ii*10+j),"T", ii/(double)ncol+tbw+shift+shiftw , j/(double)nrow+tbh, (ii+1)/(double)ncol-tbw+shift+shiftw, (1+j)/(double)nrow-tbh, 21);
	prt_hpads[padi]->SetFillColor(kCyan-8);
	prt_hpads[padi]->SetMargin(0.04,0.04,0.04,0.04);
	prt_hpads[padi]->Draw(); 
	padi++;
      }
    }
  }

  int np;
  double max=0;

  {
    double tmax;
    if(maxz==0){
      for(int p=0; p<nrow*ncol;p++){
	tmax = prt_hdigi[p]->GetBinContent(prt_hdigi[p]->GetMaximumBin());
	if(max<tmax) max = tmax;
      }
    }else{
      max = maxz;
    }
  
    if(maxz==-2 || minz==-2){ // optimize range
      for(int p=0; p<nrow*ncol;p++){
	tmax = prt_hdigi[p]->GetMaximum();
	if(max<tmax) max = tmax;
      }
      int tbins = 2000;
      TH1F *h = new TH1F("","",tbins,0,max);
      for(int p=0; p<nrow*ncol;p++){
	for(int i=0; i<8; i++){
	  for(int j=0; j<8; j++){
	    double val = prt_hdigi[p]->GetBinContent(i+1,j+1);
	    if(val!=0) h->Fill(val);

	  }
	}
      }
      double integral;
      for(int i=0; i<tbins; i++){
	integral = h->Integral(0,i);
	if(integral>0) {
	  if(minz==-2) minz = h->GetBinCenter(i);
	  break;
	} 
      }

      for(int i=tbins; i>0; i--){
	integral = h->Integral(i,tbins);
	if(integral>10) {
	  if(maxz==-2) max = h->GetBinCenter(i);
	  break;
	} 
      }
    }
  }

  prt_last_max = max;
  int nnmax(0);
  for(int p=0; p<nrow*ncol;p++){
    if(layoutId == 1 || layoutId == 4)  np =p%nrow*ncol + p/3;
    else np = p;

    if(layoutId == 6 && p>10) continue;
    
    prt_hpads[p]->cd();
    prt_hdigi[np]->Draw("col");//"col+text"
    if(maxz==-1)  max = prt_hdigi[np]->GetBinContent(prt_hdigi[np]->GetMaximumBin());
    if(nnmax<prt_hdigi[np]->GetEntries()) nnmax=np;
    prt_hdigi[np]->SetMaximum(max);
    prt_hdigi[np]->SetMinimum(minz);
  }
  
  // prt_cdigi_palette = (TPaletteAxis*)prt_hdigi[nnmax]->GetListOfFunctions()->FindObject("palette");
  // prt_cdigi_palette->SetX1NDC(0.89);
  // prt_cdigi_palette->SetY1NDC(0.1);
  // prt_cdigi_palette->SetX2NDC(0.93);
  // prt_cdigi_palette->SetY2NDC(0.90);
  
  cdigi->cd();
  TPaletteAxis* prt_cdigi_palette;  
  if(layoutId==2018 || layoutId==2023)  prt_cdigi_palette = new TPaletteAxis(0.89,0.1,0.93,0.90,(TH1 *) prt_hdigi[nnmax]);
  else if(layoutId==2032) prt_cdigi_palette = new TPaletteAxis(0.91,0.1,0.94,0.90,(TH1 *) prt_hdigi[nnmax]);
  else prt_cdigi_palette = new TPaletteAxis(0.82,0.1,0.86,0.90,(TH1 *) prt_hdigi[nnmax]);
  prt_cdigi_palette->SetName("prt_palette");  
  prt_cdigi_palette->Draw();

  cdigi->Modified();
  cdigi->Update();
  
  return cdigi;
}


TString prt_getPixData(TString s="m,p,v\n", int layoutId=1){
  int nrow = 3, ncol = 5, np, nmax=0, npix=8;
  if(layoutId ==2016) ncol=3;
  if(layoutId ==2017) ncol=4;
  if(layoutId ==2018 || layoutId ==2023) {nrow=2; ncol=4;}
  if(layoutId ==2021) ncol=4;
  if(layoutId ==2030) {nrow=4; ncol=6; npix=16;}
  if(layoutId ==2032) {nrow=4; ncol=7; npix=16;}
  if(layoutId ==2031) {nrow=3; ncol=4; npix=16;}
  
  int nnmax(0);  
  for(int p=0; p<nrow*ncol;p++){
    if(layoutId == 1 || layoutId == 4)  np =p%nrow*ncol + p/3;
    else if(layoutId == 2030)  np =p%ncol*nrow + p/ncol;
    else np = p;

    if(prt_last_maxz==-1) prt_last_max = prt_hdigi[p]->GetBinContent(prt_hdigi[p]->GetMaximumBin());
    if(nnmax<prt_hdigi[p]->GetEntries()) nnmax=p;
    prt_hdigi[p]->SetMaximum(prt_last_max);
    prt_hdigi[p]->SetMinimum(prt_last_minz);
    for(int i=1; i<=npix; i++){
      for(int j=1; j<=npix; j++){
  	double weight = (double)(prt_hdigi[p]->GetBinContent(i,j))/(double)prt_last_max*255;
  	if(weight>255) weight=255;
  	if(weight > 0) s += Form("%d,%d,%d\n", np, (i-1)*npix+j-1, (int)weight);
      }
    }
  }
  return s;
}

void prt_initDigi(int type=1){  
  if(type == 1){
    for(int m=0; m<prt_nmcp;m++){
      if(prt_hdigi[m]) prt_hdigi[m]->Reset("M");
      else{
	prt_hdigi[m] = new TH2F( Form("mcp%d", m),Form("mcp%d", m),8,0.,8.,8,0.,8.);
	prt_hdigi[m]->SetStats(0);
	prt_hdigi[m]->SetTitle(0);
	prt_hdigi[m]->GetXaxis()->SetNdivisions(10);
	prt_hdigi[m]->GetYaxis()->SetNdivisions(10);
	prt_hdigi[m]->GetXaxis()->SetLabelOffset(100);
	prt_hdigi[m]->GetYaxis()->SetLabelOffset(100);
	prt_hdigi[m]->GetXaxis()->SetTickLength(1);
	prt_hdigi[m]->GetYaxis()->SetTickLength(1);
	prt_hdigi[m]->GetXaxis()->SetAxisColor(15);
	prt_hdigi[m]->GetYaxis()->SetAxisColor(15);
      }
    }
  }
  if(type == 2){ //eic
    for(int m=0; m<prt_nmcp;m++){
      if(prt_hdigi[m]) prt_hdigi[m]->Reset("M");
      else{
	prt_hdigi[m] = new TH2F( Form("mcp%d", m),Form("mcp%d", m),16,0,16,16,0,16);
	prt_hdigi[m]->SetStats(0);
	prt_hdigi[m]->SetTitle(0);
	prt_hdigi[m]->GetXaxis()->SetNdivisions(20);
	prt_hdigi[m]->GetYaxis()->SetNdivisions(20);
	prt_hdigi[m]->GetXaxis()->SetLabelOffset(100);
	prt_hdigi[m]->GetYaxis()->SetLabelOffset(100);
	prt_hdigi[m]->GetXaxis()->SetTickLength(1);
	prt_hdigi[m]->GetYaxis()->SetTickLength(1);
	prt_hdigi[m]->GetXaxis()->SetAxisColor(15);
	prt_hdigi[m]->GetYaxis()->SetAxisColor(15);
      }
    }
  }
}

void prt_resetDigi(){
  for(int m=0; m<prt_nmcp;m++){	
    prt_hdigi[m]->Reset("M");
  }
}

void prt_axisHits800x500(TH2 * hist){
  hist->SetStats(0);
  hist->SetTitle(Form("%d hits",(int)hist->GetEntries()));
  hist->GetXaxis()->SetTitle("z, [cm]");
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(0.8);
  hist->GetYaxis()->SetTitle("y, [cm]");
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleOffset(0.7);
}

void prt_axisAngle800x500(TH2 * hist){
  hist->SetStats(0);
  hist->SetTitle(Form("%d hits",(int)hist->GetEntries()));
  hist->GetXaxis()->SetTitle("#theta, [degree]");
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(0.8);
  hist->GetYaxis()->SetTitle("photons per track, [#]");
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleOffset(0.7);
}

void prt_axisAngle800x500(TH1 * hist){
  hist->SetStats(0);
  hist->SetTitle(Form("%d hits",(int)hist->GetEntries()));
  hist->GetXaxis()->SetTitle("#theta, [degree]");
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(0.8);
  hist->GetYaxis()->SetTitle("photons per track, [#]");
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleOffset(0.7);
}

void prt_axisTime800x500(TH2 * hist){
  hist->GetXaxis()->SetTitle("time, [ns]");
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(0.8);
  hist->GetYaxis()->SetTitle("entries, #");
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleOffset(0.7);
  hist->SetLineColor(1);
}

void prt_axisTime800x500(TH1 * hist, TString xtitle = "time [ns]"){
  TGaxis::SetMaxDigits(3);
  hist->GetXaxis()->SetTitle(xtitle);
  hist->GetXaxis()->SetTitleSize(0.06);
  hist->GetXaxis()->SetTitleOffset(0.8);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetTitle("entries [#]");
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleOffset(0.7);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->SetLineColor(1);
}

void prt_setPrettyStyle(){
  // Canvas printing details: white bg, no borders.
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(0);

  // Canvas frame printing details: white bg, no borders.
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(0);

  // Plot title details: centered, no bg, no border, nice font.
  gStyle->SetTitleX(0.1);
  gStyle->SetTitleW(0.8);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFillColor(0);

  // Font details for titles and labels.
  gStyle->SetTitleFont(42, "xyz");
  gStyle->SetTitleFont(42, "pad");
  gStyle->SetLabelFont(42, "xyz");
  gStyle->SetLabelFont(42, "pad");

  // Details for stat box.
  gStyle->SetStatColor(0);
  gStyle->SetStatFont(42);
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatX(0.975);
  gStyle->SetStatY(0.9);

  // gStyle->SetOptStat(0);
}

void prt_setGStyle(TGraph *g, int id){
  int prt_coll[]={kBlack,kRed+1,kGreen,  kBlue,  4,kCyan-6,kOrange,  7,8,9,10,1,1,1};
  int prt_colm[]={kBlack,kRed+1,kGreen+2,kBlue+1,4,kCyan-6,kOrange+1,7,8,9,10,1,1,1};

  int cl=(id<10)? prt_coll[id]:id;
  int cm=(id<10)? prt_colm[id]:id;
  g->SetLineColor(cl);
  g->SetMarkerColor(cm);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.8);
  g->SetName(Form("gr_%d",id));
}

void prt_setRootPalette(int pal = 0){

  // pal =  1: rainbow\n"
  // pal =  2: reverse-rainbow\n"
  // pal =  3: amber\n"
  // pal =  4: reverse-amber\n"
  // pal =  5: blue/white\n"
  // pal =  6: white/blue\n"
  // pal =  7: red temperature\n"
  // pal =  8: reverse-red temperature\n"
  // pal =  9: green/white\n"
  // pal = 10: white/green\n"
  // pal = 11: orange/blue\n"
  // pal = 12: blue/orange\n"
  // pal = 13: white/black\n"
  // pal = 14: black/white\n"

  const int NRGBs = 5;
  const int NCont = 255;
  gStyle->SetNumberContours(NCont);

  if (pal < 1 && pal> 15) return;
  else pal--;

  double stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  double red[15][NRGBs]   = {{ 0.00, 0.00, 0.87, 1.00, 0.51 },
			       { 0.51, 1.00, 0.87, 0.00, 0.00 },
			       { 0.17, 0.39, 0.62, 0.79, 1.00 },
			       { 1.00, 0.79, 0.62, 0.39, 0.17 },
			       { 0.00, 0.00, 0.00, 0.38, 1.00 },
			       { 1.00, 0.38, 0.00, 0.00, 0.00 },
			       { 0.00, 0.50, 0.89, 0.95, 1.00 },
			       { 1.00, 0.95, 0.89, 0.50, 0.00 },
			       { 0.00, 0.00, 0.38, 0.75, 1.00 },
			       { 0.00, 0.34, 0.61, 0.84, 1.00 },
			       { 0.75, 1.00, 0.24, 0.00, 0.00 },
			       { 0.00, 0.00, 0.24, 1.00, 0.75 },
			       { 0.00, 0.34, 0.61, 0.84, 1.00 },
			       { 1.00, 0.84, 0.61, 0.34, 0.00 },
			       { 0.00, 0.00, 0.80, 1.00, 0.80 }
  };
  double green[15][NRGBs] = {{ 0.00, 0.81, 1.00, 0.20, 0.00 },		    
			       { 0.00, 0.20, 1.00, 0.81, 0.00 },
			       { 0.01, 0.02, 0.39, 0.68, 1.00 },
			       { 1.00, 0.68, 0.39, 0.02, 0.01 },
			       { 0.00, 0.00, 0.38, 0.76, 1.00 },
			       { 1.00, 0.76, 0.38, 0.00, 0.00 },
			       { 0.00, 0.00, 0.27, 0.71, 1.00 },
			       { 1.00, 0.71, 0.27, 0.00, 0.00 },
			       { 0.00, 0.35, 0.62, 0.85, 1.00 },
			       { 1.00, 0.75, 0.38, 0.00, 0.00 },
			       { 0.24, 1.00, 0.75, 0.18, 0.00 },
			       { 0.00, 0.18, 0.75, 1.00, 0.24 },
			       { 0.00, 0.34, 0.61, 0.84, 1.00 },
			       { 1.00, 0.84, 0.61, 0.34, 0.00 },
			       { 0.00, 0.85, 1.00, 0.30, 0.00 }		
  };
  double blue[15][NRGBs]  = {{ 0.51, 1.00, 0.12, 0.00, 0.00 },
			       { 0.00, 0.00, 0.12, 1.00, 0.51 },
			       { 0.00, 0.09, 0.18, 0.09, 0.00 },
			       { 0.00, 0.09, 0.18, 0.09, 0.00 },
			       { 0.00, 0.47, 0.83, 1.00, 1.00 },
			       { 1.00, 1.00, 0.83, 0.47, 0.00 },
			       { 0.00, 0.00, 0.00, 0.40, 1.00 },
			       { 1.00, 0.40, 0.00, 0.00, 0.00 },
			       { 0.00, 0.00, 0.00, 0.47, 1.00 },
			       { 1.00, 0.47, 0.00, 0.00, 0.00 },
			       { 0.00, 0.62, 1.00, 0.68, 0.12 },
			       { 0.12, 0.68, 1.00, 0.62, 0.00 },
			       { 0.00, 0.34, 0.61, 0.84, 1.00 },
			       { 1.00, 0.84, 0.61, 0.34, 0.00 },
			       { 0.60, 1.00, 0.10, 0.00, 0.00 }
  };


  TColor::CreateGradientColorTable(NRGBs, stops, red[pal], green[pal], blue[pal], NCont);
 
}
int prt_getColorId(int ind, int style =0){
  int cid = 1;
  if(style==0) {
    cid=ind+1;
    if(cid==5) cid =8;
    if(cid==3) cid =kOrange+2;
  }
  if(style==1) cid=ind+300;
  return cid;
}

int prt_shiftHist(TH1 *hist, double double_shift){
  int bins=hist->GetXaxis()->GetNbins();
  double xmin=hist->GetXaxis()->GetBinLowEdge(1);
  double xmax=hist->GetXaxis()->GetBinUpEdge(bins);
  double_shift=double_shift*(bins/(xmax-xmin));
  int shift=0;
  if(double_shift<0) shift=TMath::FloorNint(double_shift);
  if(double_shift>0) shift=TMath::CeilNint(double_shift);
  if(shift==0) return 0;
  if(shift>0){
    for(int i=1; i<=bins; i++){
      if(i+shift<=bins) hist->SetBinContent(i,hist->GetBinContent(i+shift));
      if(i+shift>bins) hist->SetBinContent(i,0);
    }
    return 0;
  }
  if(shift<0){
    for(int i=bins; i>0; i--){
      if(i+shift>0) hist->SetBinContent(i,hist->GetBinContent(i+shift));
      if(i+shift<=0) hist->SetBinContent(i,0);
    }    
    return 0;
  }
  return 1;
} 

void prt_addInfo(TString str){
  prt_info += str+"\n";
}

void prt_writeInfo(TString filename){
  std::ofstream myfile;
  myfile.open (filename);
  myfile << prt_info+"\n";
  myfile.close();
}

void prt_writeString(TString filename, TString str){
  std::ofstream myfile;
  myfile.open (filename);
  myfile << str+"\n";
  myfile.close();
  std::cout<<"output: "<<filename<<std::endl;  
}

TString prt_createDir(TString inpath=""){
  if(inpath != "") prt_savepath = inpath;
  TString finalpath = prt_savepath;

  if(finalpath =="") return "";
  
  if(prt_savepath.EndsWith("auto")) {
    TString dir = prt_savepath.ReplaceAll("auto","data");
    gSystem->mkdir(dir);
    TDatime *time = new TDatime();
    TString path(""), stime = Form("%d.%d.%d", time->GetDay(),time->GetMonth(),time->GetYear()); 
    gSystem->mkdir(dir+"/"+stime);
    for(int i=0; i<1000; i++){
      path = stime+"/"+Form("arid-%d",i);
      if(gSystem->mkdir(dir+"/"+path)==0) break;
    }
    gSystem->Unlink(dir+"/last");
    gSystem->Symlink(path, dir+"/last");
    finalpath = dir+"/"+path;
    prt_savepath=finalpath;
  }else{
    gSystem->mkdir(prt_savepath,kTRUE);
  }
  
  prt_writeInfo(finalpath+"/readme");
  return finalpath;
}

void prt_canvasPrint(TPad *c, TString name="", TString path="", int what=0){
  c->Modified();
  c->Update();
  c->Print(path+"/"+name+".png");
  if(what>0) c->Print(path+"/"+name+".C");
  if(what>1) c->Print(path+"/"+name+".pdf");
  if(what>2) c->Print(path+"/"+name+".eps");    
}

void prt_set_style(TCanvas *c){
  // prt_setRootPalette(1);
  if(fabs(c->GetBottomMargin()-0.1)<0.001) c->SetBottomMargin(0.12);
  TIter next(c->GetListOfPrimitives());
  TObject *obj;
	
  while((obj = next())){
    if(obj->InheritsFrom("TH1")){
      TH1F *hh = (TH1F*)obj;
      hh->GetXaxis()->SetTitleSize(0.06);
      hh->GetYaxis()->SetTitleSize(0.06);

      hh->GetXaxis()->SetLabelSize(0.05);
      hh->GetYaxis()->SetLabelSize(0.05);

      hh->GetXaxis()->SetTitleOffset(0.85);
      hh->GetYaxis()->SetTitleOffset(0.76);

      if(c->GetWindowHeight()>700){
	hh->GetXaxis()->SetTitleSize(0.045);
	hh->GetYaxis()->SetTitleSize(0.045);

	hh->GetXaxis()->SetLabelSize(0.035);
	hh->GetYaxis()->SetLabelSize(0.035);

	hh->GetXaxis()->SetTitleOffset(0.85);
	hh->GetYaxis()->SetTitleOffset(0.98);
	c->SetRightMargin(0.12);
      }
					   
      if(c->GetWindowWidth()<=600){
      	hh->GetXaxis()->SetTitleSize(0.05);
      	hh->GetYaxis()->SetTitleSize(0.05);

      	hh->GetXaxis()->SetLabelSize(0.04);
      	hh->GetYaxis()->SetLabelSize(0.04);

      	hh->GetXaxis()->SetTitleOffset(0.85);
      	hh->GetYaxis()->SetTitleOffset(0.92);
      	c->SetRightMargin(0.12);
      }


      if(fabs(c->GetBottomMargin()-0.12)<0.001){
	TPaletteAxis *palette = (TPaletteAxis*)hh->GetListOfFunctions()->FindObject("palette");
	if(palette) {
	  palette->SetY1NDC(0.12);
	  c->Modified();
	}
      }
      c->Modified();
      c->Update();
    }
    
    if(obj->InheritsFrom("TGraph")){
      TGraph *gg = (TGraph*)obj;
      gg->GetXaxis()->SetLabelSize(0.05);
      gg->GetXaxis()->SetTitleSize(0.06);
      gg->GetXaxis()->SetTitleOffset(0.84);

      gg->GetYaxis()->SetLabelSize(0.05);
      gg->GetYaxis()->SetTitleSize(0.06);
      gg->GetYaxis()->SetTitleOffset(0.8);
    }
 
    if(obj->InheritsFrom("TMultiGraph")){
      TMultiGraph *gg = (TMultiGraph*)obj;
      gg->GetXaxis()->SetLabelSize(0.05);
      gg->GetXaxis()->SetTitleSize(0.06);
      gg->GetXaxis()->SetTitleOffset(0.84);

      gg->GetYaxis()->SetLabelSize(0.05);
      gg->GetYaxis()->SetTitleSize(0.06);
      gg->GetYaxis()->SetTitleOffset(0.8);
    }
    
    if(obj->InheritsFrom("TF1")){
      TF1 *f = (TF1*)obj;
      f->SetNpx(500);
    }
  }
}

void prt_save(TPad *c= NULL,TString path="", int what=0, int style=0){
  TString name = c->GetName();
  Bool_t batch = gROOT->IsBatch();
  gROOT->SetBatch(1);
  
  if(c && path != "") {
    int w = 800, h = 400;
    if(style != -1){
      if(style == 1) {w = 800; h = 500;}
      if(style == 2) {w = 800; h = 600;}
      if(style == 3) {w = 800; h = 400;}
      if(style == 5) {w = 800; h = 900;} 
      if(style == 0){ 
    	w = ((TCanvas*)c)->GetWindowWidth();
    	h = ((TCanvas*)c)->GetWindowHeight();
      }
      
      TCanvas *cc;
      if(TString(c->GetName()).Contains("hp") || TString(c->GetName()).Contains("cdigi")) {
	cc = prt_drawDigi(prt_last_layoutId,prt_last_maxz,prt_last_minz);
	cc->SetCanvasSize(800,400);
	if(name.Contains("=")) name =  name.Tokenize('=')->First()->GetName();	
      }else{
      	cc = new TCanvas(TString(c->GetName())+"exp","cExport",0,0,w,h);
      	cc = (TCanvas*) c->DrawClone();
	cc->SetCanvasSize(w,h);
      }
      
      if(style == 0) prt_set_style(cc);
      
      prt_canvasPrint(cc,name,path,what);
    }else{
      c->SetCanvasSize(w,h);
      prt_canvasPrint(c,name,path,what);
    }
  }
  
  gROOT->SetBatch(batch);
}

TString prt_createSubDir(TString dir="dir"){
  gSystem->mkdir(dir);
  return dir;
}

TList *prt_canvaslist;

TCanvas *prt_canvasGet(TString name="c"){
  TIter next(prt_canvaslist);
  TCanvas *c = 0;
  while((c = (TCanvas*) next())){
    if(c->GetName()==name || name=="*") break;
  }
  return c;
}

void prt_canvasAdd(TString name="c",int w=800, int h=400){
  if(!prt_canvasGet(name)){
    if(!prt_canvaslist) prt_canvaslist = new TList();
    TCanvas *c = new TCanvas(name,name,0,0,w,h); 
    prt_canvaslist->Add(c);
  }
}

void prt_canvasAdd(TCanvas *c){
  if(!prt_canvaslist) prt_canvaslist = new TList();
  c->cd();
  prt_canvaslist->Add(c);
}

void prt_canvasCd(){
  
}

void prt_canvasDel(TString name="c"){
  TIter next(prt_canvaslist);
  TCanvas *c=0;
  while((c = (TCanvas*) next())){
    if(c->GetName()==name || name=="*") prt_canvaslist->Remove(c);
    c->Delete();
  }
}

// style = 0 - for web blog
// style = 1 - for talk 
// what = 0 - save in png
// what = 1 - save in png, C
// what = 2 - save in png, C, pdf
// what = 3 - save in png, C, pdf, eps 
void prt_canvasSave(int what=1, int style=0, Bool_t rm=false){
  TIter next(prt_canvaslist);
  TCanvas *c=0;
  TString path = prt_createDir();
  while((c = (TCanvas*) next())){
    prt_set_style(c);
    prt_save(c, path, what,style);
    if(rm){
      prt_canvaslist->Remove(c);
      c->Close();
    }
  }
}

// path - folder for saving
void prt_canvasSave(TString path, int what=1, int style=0, Bool_t rm=false){
  prt_savepath = path;
  prt_canvasSave(what,style,rm);
}

void prt_set_style(){
  TIter next(prt_canvaslist);
  TCanvas *c=0;
  TString path = prt_createDir();
  while((c = (TCanvas*) next())){
    prt_set_style(c);
  }
}

void prt_waitPrimitive(TString name, TString prim=""){
  TIter next(prt_canvaslist);
  TCanvas *c=0;
  while((c = (TCanvas*) next())){
    if(TString(c->GetName())==name){
      c->Modified(); 
      c->Update(); 
      c->WaitPrimitive(prim);
    }
  }
}

double prt_integral(TH1F *h,double xmin, double xmax){
  TAxis *axis = h->GetXaxis();
  int bmin = axis->FindBin(xmin);
  int bmax = axis->FindBin(xmax);
  double integral = h->Integral(bmin,bmax);
  integral -= h->GetBinContent(bmin)*(xmin-axis->GetBinLowEdge(bmin))/axis->GetBinWidth(bmin);
  integral -= h->GetBinContent(bmax)*(axis->GetBinUpEdge(bmax)-xmax)/axis->GetBinWidth(bmax);
  return integral;
}

void prt_normalize(TH1F* hists[],int size){
  // for(int i=0; i<size; i++){
  //   hists[i]->Scale(1/hists[i]->Integral(), "width"); 
  // }
  
  double max = 0;
  double min = 0;
  for(int i=0; i<size; i++){
    double tmax =  hists[i]->GetBinContent(hists[i]->GetMaximumBin());
    double tmin = hists[i]->GetMinimum();
    if(tmax>max) max = tmax;
    if(tmin<min) min = tmin;
  }
  max += 0.05*max;
  for(int i=0; i<size; i++){
    hists[i]->GetYaxis()->SetRangeUser(min,max);
  }
}

void prt_normalizeto(TH1F* hists[],int size, double max=1){
  for(int i=0; i<size; i++){
    double tmax =  hists[i]->GetBinContent(hists[i]->GetMaximumBin());
    if(tmax>0)hists[i]->Scale(max/tmax);
  }
}

void prt_normalize(TH1F* h1,TH1F* h2){
  double max = (h1->GetMaximum()>h2->GetMaximum())? h1->GetMaximum() : h2->GetMaximum();
  max += max*0.1;
  h1->GetYaxis()->SetRangeUser(0,max);
  h2->GetYaxis()->SetRangeUser(0,max);
}

// just x for now
TGraph* prt_smooth(TGraph* g,int smoothness=1){
  double x, y;
  int n = g->GetN();
  TH1F *h = new TH1F("h","h",g->GetN(),0,n);
  TGraph *gr = new TGraph();
  gr->SetName(g->GetName());
  for(auto i=0; i<n; i++){
    g->GetPoint(i,x,y);
    h->Fill(i,x);
  }

  h->Smooth(smoothness);
  
  for(auto i=0;i<n;i++){
    g->GetPoint(i,x,y);
    gr->SetPoint(i,h->GetBinContent(i),y);
  }
  return gr;
}

int prt_get_pid(int pdg){
  int pid=0;
  if(pdg==-11)   pid=0; //e
  if(pdg==13)   pid=1; //mu
  if(pdg==211)  pid=2; //pi
  if(pdg==321)  pid=3; //K
  if(pdg==2212) pid=4; //p
  return pid;
}

double prt_get_momentum_from_tof(double dist,double dtof){
  double s = dtof*0.299792458/dist;
  double x = s*s;
  double a = prt_mass[2]*prt_mass[2]; //pi
  double b = prt_mass[4]*prt_mass[4]; //p
  double p = sqrt((a - 2*sqrt(a*a+a*b*x-2*a*b+b*b)/s + b)/(x - 4));

  return p;
}

// return TOF difference [ns] for mominum p [GeV] and flight path l [m]
double prt_get_tof_diff(int pid1=211, int pid2=321, double p=1, double l=2){
  double c = 299792458;
  double m1 = prt_mass[prt_get_pid(pid1)];
  double m2 = prt_mass[prt_get_pid(pid2)];
  double td = l*(sqrt(p*p+m1*m1)-sqrt(p*p+m2*m2))/(p*c)*1E9;

  // relativistic
  // td = l*(m1*m1-m2*m2)/(2*p*p*c)*1E9;
  
  return td;
}

bool prt_ispath(TString path){
  Long_t *id(0),*size(0),*flags(0),*modtime(0);
  return !gSystem->GetPathInfo(path,id,size,flags,modtime);
}

int prt_get3digit(TString str){
  TPRegexp e("[0-9][0-9][0-9]");
  return ((TString)str(e)).Atoi();
}

#endif
