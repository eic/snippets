void Perform_Check_Entries_ep(){
	
TFile *f = TFile::Open("../../D0_Jan25_Signal/HF_reco/helix/SignalD0.root");;
TTree *treeMLSig = (TTree*)f->Get("treeMLSig");;
printf("=================Signal Checks=================\n");;
printf("Signal Entries:------> y:(-3, -1) \t pt: (0.,1.0)-----------------\n");;
treeMLSig->Draw("mass_D0>>h","(mass_D0 > 1.6 && mass_D0 < 2.5) && d0xy_pi>0.02 && d0xy_k>0.02 && dca_12 < 0.07 && costheta > 0.95 && decay_length > 0.05 && dca_D0 < 0.1 && (y_D0>-3 && y_D0<-1) && (pt_D0>0 && pt_D0<1)");
TH1D *h = (TH1D*)gDirectory->Get("h"); 
printf("%1.1f\n",h->GetEntries());

printf("Signal Entries:------> y:(-3, -1) \t pt: (1.,2.0)-----------------\n");;
treeMLSig->Draw("mass_D0>>h","(mass_D0 > 1.6 && mass_D0 < 2.5) && d0xy_pi>0.02 && d0xy_k>0.02 && dca_12 < 0.07 && costheta > 0.95 && decay_length > 0.05 && dca_D0 < 0.1 && (y_D0>-3 && y_D0<-1) && (pt_D0>1 && pt_D0<2)");
h = (TH1D*)gDirectory->Get("h"); 
printf("%1.1f\n",h->GetEntries());

printf("Signal Entries:------> y:(-3, -1) \t pt: (2.0,5.0)-----------------\n");;
treeMLSig->Draw("mass_D0>>h","(mass_D0 > 1.6 && mass_D0 < 2.5) && d0xy_pi>0.02 && d0xy_k>0.02 && dca_12 < 0.07 && costheta > 0.95 && decay_length > 0.05 && dca_D0 < 0.1 && (y_D0>-3 && y_D0<-1) && (pt_D0>2 && pt_D0<5)");
h = (TH1D*)gDirectory->Get("h"); 
printf("%1.1f\n",h->GetEntries());

printf("Signal Entries:------> y:(-1, 1) \t pt: (0.,1.0)-----------------\n");;
treeMLSig->Draw("mass_D0>>h","(mass_D0 > 1.6 && mass_D0 < 2.5) && (d0xy_pi>0.02 && d0xy_pi<10.) && (d0xy_k>0.02 && d0xy_k<10.) && decay_length <100. && (y_D0>-1 && y_D0<1) && (pt_D0>0 && pt_D0<1)");
h = (TH1D*)gDirectory->Get("h"); 
printf("%1.1f\n",h->GetEntries());

printf("Signal Entries:------> y:(-1, 1) \t pt: (1.,2.0)-----------------\n");;
treeMLSig->Draw("mass_D0>>h","(mass_D0 > 1.6 && mass_D0 < 2.5) && (d0xy_pi>0.02 && d0xy_pi<10.) && (d0xy_k>0.02 && d0xy_k<10.) && decay_length <100. && (y_D0>-1 && y_D0<1) && (pt_D0>1 && pt_D0<2)");
h = (TH1D*)gDirectory->Get("h"); 
printf("%1.1f\n",h->GetEntries());

printf("Signal Entries:------> y:(-1, 1) \t pt: (2.0,5.0)-----------------\n");;
treeMLSig->Draw("mass_D0>>h","(mass_D0 > 1.6 && mass_D0 < 2.5) && (d0xy_pi>0.02 && d0xy_pi<10.) && (d0xy_k>0.02 && d0xy_k<10.) && decay_length <100. && (y_D0>-1 && y_D0<1) && (pt_D0>2 && pt_D0<5)");
h = (TH1D*)gDirectory->Get("h"); 
printf("%1.1f\n",h->GetEntries());

printf("Signal Entries:------> y:(1, 3) \t pt: (0.,1.0)-----------------\n");;
treeMLSig->Draw("mass_D0>>h","(mass_D0 > 1.6 && mass_D0 < 2.5) && (d0xy_pi>0.02 && d0xy_pi<10.) && (d0xy_k>0.02 && d0xy_k<10.) && decay_length <100. && (y_D0>1 && y_D0<3) && (pt_D0>0 && pt_D0<1)");
h = (TH1D*)gDirectory->Get("h"); 
printf("%1.1f\n",h->GetEntries());

printf("Signal Entries:------> y:(1, 3) \t pt: (1.,2.0)-----------------\n");;
treeMLSig->Draw("mass_D0>>h","(mass_D0 > 1.6 && mass_D0 < 2.5) && (d0xy_pi>0.02 && d0xy_pi<10.) && (d0xy_k>0.02 && d0xy_k<10.) && decay_length <100. && (y_D0>1 && y_D0<3) && (pt_D0>1 && pt_D0<2)");
h = (TH1D*)gDirectory->Get("h"); 
printf("%1.1f\n",h->GetEntries());

printf("Signal Entries:------> y:(1, 3) \t pt: (2.0,5.0)-----------------\n");;
treeMLSig->Draw("mass_D0>>h","(mass_D0 > 1.6 && mass_D0 < 2.5) && (d0xy_pi>0.02 && d0xy_pi<10.) && (d0xy_k>0.02 && d0xy_k<10.) && decay_length <100. && (y_D0>1 && y_D0<3) && (pt_D0>2 && pt_D0<5)");
h = (TH1D*)gDirectory->Get("h"); 
printf("%1.1f\n",h->GetEntries());

printf("=================Bkg Checks=================\n");;
f = TFile::Open("../../D0_Jan25_Bkg/HF_reco/helix/BkgD0.root");;
TTree *treeMLBkg = (TTree*)f->Get("treeMLBkg");;	

printf("Bkg Entries:------> y:(-3, -1) \t pt: (0.,1.0)-----------------\n");;
treeMLBkg->Draw("mass_D0>>h","(mass_D0 > 1.6 && mass_D0 < 2.5) && d0xy_pi>0.02 && d0xy_k>0.02 && dca_12 < 0.07 && costheta > 0.95 && decay_length > 0.05 && dca_D0 < 0.1 && (y_D0>-3 && y_D0<-1) && (pt_D0>0 && pt_D0<1)");
h = (TH1D*)gDirectory->Get("h"); 
printf("%1.1f\n",h->GetEntries());

printf("Bkg Entries:------> y:(-3, -1) \t pt: (1.,2.0)-----------------\n");;
treeMLBkg->Draw("mass_D0>>h","(mass_D0 > 1.6 && mass_D0 < 2.5) && d0xy_pi>0.02 && d0xy_k>0.02 && dca_12 < 0.07 && costheta > 0.95 && decay_length > 0.05 && dca_D0 < 0.1 && (y_D0>-3 && y_D0<-1) && (pt_D0>1 && pt_D0<2)");
h = (TH1D*)gDirectory->Get("h"); 
printf("%1.1f\n",h->GetEntries());

printf("Bkg Entries:------> y:(-3, -1) \t pt: (2.0,5.0)-----------------\n");;
treeMLBkg->Draw("mass_D0>>h","(mass_D0 > 1.6 && mass_D0 < 2.5) && d0xy_pi>0.02 && d0xy_k>0.02 && dca_12 < 0.07 && costheta > 0.95 && decay_length > 0.05 && dca_D0 < 0.1 && (y_D0>-3 && y_D0<-1) && (pt_D0>2 && pt_D0<5)");
h = (TH1D*)gDirectory->Get("h"); 
printf("%1.1f\n",h->GetEntries());

printf("Bkg Entries:------> y:(-1, 1) \t pt: (0.,1.0)-----------------\n");;
treeMLBkg->Draw("mass_D0>>h","(mass_D0 > 1.6 && mass_D0 < 2.5) && (d0xy_pi>0.02 && d0xy_pi<10.) && (d0xy_k>0.02 && d0xy_k<10.) && decay_length <100. && (y_D0>-1 && y_D0<1) && (pt_D0>0 && pt_D0<1)");
h = (TH1D*)gDirectory->Get("h"); 
printf("%1.1f\n",h->GetEntries());

printf("Bkg Entries:------> y:(-1, 1) \t pt: (1.,2.0)-----------------\n");;
treeMLBkg->Draw("mass_D0>>h","(mass_D0 > 1.6 && mass_D0 < 2.5) && (d0xy_pi>0.02 && d0xy_pi<10.) && (d0xy_k>0.02 && d0xy_k<10.) && decay_length <100. && (y_D0>-1 && y_D0<1) && (pt_D0>1 && pt_D0<2)");
h = (TH1D*)gDirectory->Get("h"); 
printf("%1.1f\n",h->GetEntries());

printf("Bkg Entries:------> y:(-1, 1) \t pt: (2.0,5.0)-----------------\n");;
treeMLBkg->Draw("mass_D0>>h","(mass_D0 > 1.6 && mass_D0 < 2.5) && (d0xy_pi>0.02 && d0xy_pi<10.) && (d0xy_k>0.02 && d0xy_k<10.) && decay_length <100. && (y_D0>-1 && y_D0<1) && (pt_D0>2 && pt_D0<5)");
h = (TH1D*)gDirectory->Get("h"); 
printf("%1.1f\n",h->GetEntries());

printf("Bkg Entries:------> y:(1, 3) \t pt: (0.,1.0)-----------------\n");;
treeMLBkg->Draw("mass_D0>>h","(mass_D0 > 1.6 && mass_D0 < 2.5) && (d0xy_pi>0.02 && d0xy_pi<10.) && (d0xy_k>0.02 && d0xy_k<10.) && decay_length <100. && (y_D0>1 && y_D0<3) && (pt_D0>0 && pt_D0<1)");
h = (TH1D*)gDirectory->Get("h"); 
printf("%1.1f\n",h->GetEntries());

printf("Bkg Entries:------> y:(1, 3) \t pt: (1.,2.0)-----------------\n");;
treeMLBkg->Draw("mass_D0>>h","(mass_D0 > 1.6 && mass_D0 < 2.5) && (d0xy_pi>0.02 && d0xy_pi<10.) && (d0xy_k>0.02 && d0xy_k<10.) && decay_length <100. && (y_D0>1 && y_D0<3) && (pt_D0>1 && pt_D0<2)");
h = (TH1D*)gDirectory->Get("h"); 
printf("%1.1f\n",h->GetEntries());

printf("Bkg Entries:------> y:(1, 3) \t pt: (2.0,5.0)-----------------\n");;
treeMLBkg->Draw("mass_D0>>h","(mass_D0 > 1.6 && mass_D0 < 2.5) && (d0xy_pi>0.02 && d0xy_pi<10.) && (d0xy_k>0.02 && d0xy_k<10.) && decay_length <100. && (y_D0>1 && y_D0<3) && (pt_D0>2 && pt_D0<5)");
h = (TH1D*)gDirectory->Get("h"); 
printf("%1.1f\n",h->GetEntries());	
}
