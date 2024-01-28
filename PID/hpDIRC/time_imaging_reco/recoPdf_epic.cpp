// Author: Nilanga Wickramaarachchi
// Based on https://github.com/rdom/eicdirc - by Roman Dzhygadlo

#include "prttools.cpp"
#include <TVirtualFitter.h>
#include <TKey.h>
#include <TRandom.h>

const int nch(24*256);
TF1 *pdff[nch],*pdfs[nch];
TH1F *hpdff[nch],*hpdfs[nch];
double nhits_cut, nhits_cut_val;

void recoPdf_epic(TString in="eicrecon.root", TString pdf="eicrecon.pdf.root", int theta_ang = 30, double timeres=0.1, int pid =3, TString nameid="") // use pid = 0 for e+
{

  TGaxis::SetMaxDigits(4);
  
  TCanvas *cc = new TCanvas("cc","cc");

  TH1F *hl[5],*hll[5],*hnph[5];
  TH1F *Hist_leadtime = new TH1F("leadtime",";LE time [ns]; entries [#]", 2000, 0, 100);
  for(int i=0; i<5; i++){
    hl[i] = new TH1F(Form("hl_%d",i),";LE time [ns]; entries [#]", 2000,0,100);    
    //hll[i]= new TH1F(Form("ll_i%d",i),";ln L("+prt_lname[pid]+") - ln L(#pi); entries [#]",240,-80,80);
    //hll[i]= new TH1F(Form("ll_i%d",i),";ln L("+prt_lname[pid]+") - ln L(#pi); entries [#]",300,-300,300);
    hll[i]= new TH1F(Form("ll_i%d",i),";ln L("+prt_lname[pid]+") - ln L(#pi); entries [#]",400,-200,200);
    hnph[i] = new TH1F(Form("hnph_%d",i),";multiplicity [#]; entries [#]", 220,0,220);
    hnph[i]->SetLineColor(prt_color[i]);
    hll[i]->SetLineColor(prt_color[i]);
  }

  TFile f(pdf);
  TTree *tree_nph = (TTree*) f.Get("nph_pip");
  tree_nph->SetBranchAddress("nhits_cut",&nhits_cut);               
  tree_nph->GetEntry(0);
  nhits_cut_val = nhits_cut;
  //std::cout << "number of hits cut at " << nhits_cut_val << std::endl;
  
  int rebin=timeres/(100/2000.);
  std::cout<<"rebin "<<rebin <<std::endl;
  
  int integ1(0), integ2(0);
  for(int i=0; i < nch; i++){
    hpdff[i] = (TH1F*)f.Get(Form("hf_%d",i));
    hpdfs[i] = (TH1F*)f.Get(Form("hs_%d",i));
    if(rebin >0) hpdff[i]->Rebin(rebin);
    if(rebin >0) hpdfs[i]->Rebin(rebin);
    integ1+= hpdff[i]->Integral();
    integ2+= hpdfs[i]->Integral();
    //hpdff[i]->Smooth(1);
    //hpdfs[i]->Smooth(1);
  }

  prt_ch = new TChain("hpDIRCrawHits/dirctree");
  prt_ch->Add(in);
  int nEvents = prt_ch->GetEntries();
  std::cout << "Entries = " << nEvents << std::endl;
  
  int hit_size = 0;
  int Particle_id = 0;

  const int arr_size = 10000;

  int mcp_num[arr_size], pixel_id[arr_size];
  Double_t lead_time[arr_size];
  double track_momentum_at_bar[arr_size][3];
  int track_id[arr_size];

  prt_ch->SetBranchAddress("nhits", &hit_size);
  prt_ch->SetBranchAddress("mcp_id", &mcp_num);
  prt_ch->SetBranchAddress("pixel_id", &pixel_id);  
  prt_ch->SetBranchAddress("hit_time", &lead_time);

  int printstep=1000;
  double time = 0;
  
  double nph[5]={0};
  double sigma_nph[5]={0};
  double nph_err[5] = {0};
  double sigma_nph_err[5] = {0};
  double sep_err;

  int tnph(0),totalf(0),totals(0), ch(0);
  

  for (int ievent=0; ievent < nEvents; ievent++)
    {
      prt_ch->GetEntry(ievent);

      if(ievent < nEvents/2) Particle_id = 211;
      else Particle_id = 321;

      //if((Particle_id == 211 && ievent > 2499) || (Particle_id == 11 && ievent > 27499)) continue;
      if((Particle_id == 211 && ievent > 2499) || (Particle_id == 321 && ievent > 27499)) continue;

      int nHits = hit_size;
      	
      //if(nHits < nhits_cut_val) continue;
      if(ievent%printstep==0 && ievent!=0) cout<< "Event # "<< ievent<< " # hits "<< nHits <<endl;

      
      int id = prt_get_pid(Particle_id);    
      double aminf,amins, sum(0),sumf(0),sums(0);
      tnph = 0;    
      if(hll[id]->GetEntries()>4000) continue;
    
      for(int h=0; h < nHits; h++)
	{      
	  if(pixel_id[h] < 0) continue;
	  ch = 256 * mcp_num[h] + pixel_id[h];
	  time = lead_time[h] + gRandom->Gaus(0, timeres);
	  
	  aminf = hpdff[ch]->GetBinContent(hpdff[ch]->FindBin(time)); 
	  amins = hpdfs[ch]->GetBinContent(hpdfs[ch]->FindBin(time));
	  tnph++;
      
	  //double noise = 1e-6; 
	  double noise = 0.5e-6; 

	  sumf+=TMath::Log((aminf+noise));
	  sums+=TMath::Log((amins+noise));    
            
	  hl[id]->Fill(time);
	  if(id == 2) Hist_leadtime->Fill(time);
	}

      if(tnph>4) hnph[id]->Fill(tnph);
      sum = sumf-sums;
      if(fabs(sum)<0.1) continue;
        
      hll[id]->Fill(sum);      
    }
  
  gStyle->SetOptStat(0);

  prt_theta = theta_ang;

  //TString name = Form("%1.0f_%1.2f_pik_6GeV_epic",prt_theta,timeres);
  TString name = "pik_6GeV_epic";
  //TString name = "pie_1.2GeV_epic"; 

  prt_canvasAdd("nph_"+name,800,400);
  for(int i=0; i<5; i++)
    {
      if(hnph[i]->GetEntries()<20) continue;
      TFitResultPtr r = hnph[i]->Fit("gaus","SQ","",5,250);
      nph[i] = r->Parameter(1);
      nph_err[i] = r->ParError(1);
      sigma_nph[i] = r->Parameter(2);
      sigma_nph_err[i] = r->ParError(2);      
    }

  hnph[2]->Draw();
  hnph[pid]->Draw("sames hist");

  prt_canvasGet("nph_"+name)->Update();
  
  TPaveText *pt = new TPaveText();
  pt->AddText(Form("N(#pi^{+}) = %1.2f #pm %1.2f", nph[2], nph_err[2]));
  ((TText*)pt->GetListOfLines()->Last())->SetTextColor(prt_color[2]);

  pt->AddText(Form("N(K^{+}) = %1.2f #pm %1.2f", nph[3], nph_err[3]));
  //pt->AddText(Form("N(e^{+}) = %1.2f #pm %1.2f", nph[0], nph_err[0]));
  ((TText*)pt->GetListOfLines()->Last())->SetTextColor(prt_color[3]);
  //((TText*)pt->GetListOfLines()->Last())->SetTextColor(prt_color[0]);

  pt->SetX1NDC(0.12);
  pt->SetX2NDC(0.24);
  pt->SetY1NDC(0.65);
  pt->SetY2NDC(0.85);
  pt->SetShadowColor(0);
  pt->SetFillColor(0);
  pt->SetLineColor(0);
  pt->Draw();

  prt_canvasAdd("sep_"+name,800,500);  
  prt_normalize(hll,5);
  
  TF1 *ff;
  double m1(0),m2(0),s1(100),s2(100); 
  double dm1=0,dm2=0,ds1=0,ds2=0;
  if(hll[2]->GetEntries()>10){
    hll[2]->Fit("gaus","Sq");
    ff = hll[2]->GetFunction("gaus");
    ff->SetLineColor(1);
    m1=ff->GetParameter(1);
    s1=ff->GetParameter(2);
    dm1=ff->GetParError(1);
    ds1=ff->GetParError(2);

    if(pid==0){ //handle tails      
      hll[2]->Fit("gaus","S","",-300, m1+1.8*s1);
      ff = hll[2]->GetFunction("gaus");
      m1=ff->GetParameter(1);
      s1=ff->GetParameter(2);
      dm1=ff->GetParError(1);
      ds1=ff->GetParError(2);
    }
  }

  if(hll[pid]->GetEntries()>10){
    hll[pid]->Fit("gaus","Sq");
    ff = hll[pid]->GetFunction("gaus");
    ff->SetLineColor(1);
    m2=ff->GetParameter(1);
    s2=ff->GetParameter(2);
    dm2=ff->GetParError(1);
    ds2=ff->GetParError(2);

    if(pid==0){ ///handle tails
      hll[pid]->Fit("gaus","S","",m2-1.8*s2, 300);
      ff = hll[pid]->GetFunction("gaus");      
      m2=ff->GetParameter(1);
      s2=ff->GetParameter(2);
      dm2=ff->GetParError(1);
      ds2=ff->GetParError(2);
    }
  }

  double sep = (fabs(m1-m2))/(0.5*(s1+s2));
  
  double e1,e2,e3,e4;
  e1=2/(s1+s2)*dm1;
  e2=-2/(s1+s2)*dm2;
  e3=-((2*(m1 - m2))/((s1 + s2)*(s1 + s2)))*ds1;
  e4=-((2*(m1 - m2))/((s1 + s2)*(s1 + s2)))*ds2;
  sep_err=sqrt(e1*e1+e2*e2+e3*e3+e4*e4); 
  
  /*std::cout << "mean1 = " << m1 << " +/- " << dm1 << std::endl;
  std::cout << "mean2 = " << m2 << " +/- " << dm2 <<std::endl;
  std::cout << "sigma1 = " << s1 << " +/- " << ds1 <<std::endl;
  std::cout << "sigma2 = " << s2 << " +/- " << ds2 <<std::endl;
  */
  std::cout<<in<<" separation "<< sep << "+/-" << sep_err << std::endl;
  
  hll[2]->SetTitle(Form("#theta = %1.2f deg       separation = (%1.2f #pm %1.2f)#sigma",prt_theta, sep, sep_err));
  hll[2]->Draw();
  hll[pid]->Draw("same");

  for(int i=0; i<5; i++){
    hl[i]->Scale(1/hl[i]->GetMaximum());
  }
    
  prt_normalize(hl,3);
  prt_canvasAdd("tim_"+name,800,400);
  hl[2]->Draw();
  hl[pid]->SetLineColor(4);
  hl[2]->Draw("same");
  hl[pid]->SetLineColor(2);
  hl[pid]->Draw("same");

  prt_canvasAdd("time_pip_"+name,800,400);
  gStyle->SetOptStat(1);
  Hist_leadtime->Draw();

  prt_canvasSave("."+nameid,0);
  
  TFile fc(in+"_r.root","recreate");
  TTree *tc = new TTree("reco","reco");
  tc->Branch("theta",&prt_theta,"theta/D");
  tc->Branch("sep",&sep,"sep/D");
  tc->Branch("sep_err",&sep_err,"sep_err/D");
  tc->Branch("timeres",&timeres,"timeres/D");
  tc->Branch("nph",&nph,"nph[5]/D");
  tc->Branch("nph_err",&nph_err,"nph_err[5]/D");
  tc->Branch("sigma_nph",&sigma_nph,"sigma_nph[5]/D");
  tc->Branch("sigma_nph_err",&sigma_nph_err,"sigma_nph_err[5]/D");
  tc->Fill();
  tc->Write();
}
