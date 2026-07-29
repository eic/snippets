// Code to check the orgin of D0 (prompt vs non-prompt) at the simulation level.
// Shyam Kumar; shyam.kumar@ba.infn.it; shyam055119@gmail.com

int GetD0OriginFlavor(int imc, const TTreeReaderArray<int> &mcPartPdg, const TTreeReaderArray<unsigned int> &mcPartParent_begin, const TTreeReaderArray<unsigned int> &mcPartParent_end, const TTreeReaderArray<int> &mcPartParent_index);
void Check_Origin_D0meson(TString filelist = "test.list"){
     gStyle->SetOptStat(0); 
    // Read all files
	int nfiles = 0;
     char filename[512];
	TChain *chain = new TChain("events");

	ifstream *inputstream = new ifstream;
     inputstream->open(filelist.Data());
     if(!inputstream) printf("[Error] Cannot open file list: %s\n", filelist.Data());
     
	while (inputstream->good())
	  {
	  inputstream->getline(filename, 512);
	  if (inputstream->good())
	  {
	    TFile *ftmp = TFile::Open(filename, "READ");
	    if (!ftmp || !ftmp->IsOpen() || !ftmp->GetNkeys())
	    {
		 printf("[e] Skipping bad file: %s\n", filename);
		 if (ftmp) { ftmp->Close(); delete ftmp; }
		 continue; 
	    }
	    cout << "[i] Add " << nfiles << "th file: " << filename << endl;
	    chain->Add(filename);
	    nfiles++;

	    ftmp->Close(); 
	    delete ftmp;
	  }
	}

    inputstream->close();
    printf("[i] Read in %d files with %lld events in total\n", nfiles, chain->GetEntries());
    TTreeReader myReader(chain);

     // MC Track Properties
     TTreeReaderArray<int> mcPartGenStatus = {myReader, "MCParticles.generatorStatus"};
     TTreeReaderArray<int> mcPartPdg = {myReader, "MCParticles.PDG"};
     TTreeReaderArray<double> mcPartMass = {myReader, "MCParticles.mass"};
	TTreeReaderArray<double> mcEndPointX = {myReader, "MCParticles.endpoint.x"};
	TTreeReaderArray<double> mcEndPointY = {myReader, "MCParticles.endpoint.y"};
	TTreeReaderArray<double> mcEndPointZ = {myReader, "MCParticles.endpoint.z"};
	TTreeReaderArray<double> mcMomPx = {myReader, "MCParticles.momentum.x"};
     TTreeReaderArray<double> mcMomPy = {myReader, "MCParticles.momentum.y"};
     TTreeReaderArray<double> mcMomPz = {myReader, "MCParticles.momentum.z"};
     TTreeReaderArray<double> mcPartVx = {myReader, "MCParticles.vertex.x"};
     TTreeReaderArray<double> mcPartVy = {myReader, "MCParticles.vertex.y"};
     TTreeReaderArray<double> mcPartVz = {myReader, "MCParticles.vertex.z"};
     TTreeReaderArray<int> mcPartParent_index = {myReader, "_MCParticles_parents.index"};
     TTreeReaderArray<unsigned int> mcPartDaughter_begin = {myReader, "MCParticles.daughters_begin"};
     TTreeReaderArray<unsigned int> mcPartParent_begin = {myReader, "MCParticles.parents_begin"};
     TTreeReaderArray<unsigned int> mcPartParent_end = {myReader, "MCParticles.parents_end"};
     TTreeReaderArray<unsigned int> mcPartDaughter_end = {myReader, "MCParticles.daughters_end"};
     TTreeReaderArray<int> mcPartDaughter_index = {myReader, "_MCParticles_daughters.index"};
     
     TCanvas *c = new TCanvas("c","c",2000,1200);
     c->SetMargin(0.10, 0.05 ,0.12,0.07);
     TH1F *hEventStat = new TH1F("hEventStat", ";;Entries (a.u.)", 6, 0, 6);
     hEventStat->GetXaxis()->SetBinLabel(1, "MC events");
     hEventStat->GetXaxis()->SetBinLabel(2, "D^{0}");
     hEventStat->GetXaxis()->SetBinLabel(3, "NonPrompt D^{0}");
     hEventStat->GetXaxis()->SetBinLabel(4, "Prompt D^{0}");
     hEventStat->GetXaxis()->SetBinLabel(5, "originHF==#pm 5 (Beauty)");
     hEventStat->GetXaxis()->SetBinLabel(6, "originHF==#pm 4 (Charm)");
     hEventStat->SetLineWidth(2);     
     hEventStat->SetMarkerSize(2); 
     hEventStat->SetTitleSize(0.04,"XY"); 
     hEventStat->SetLabelSize(0.04,"XY");   

  TFile *fout = new TFile("output_file.root","recreate");
  TH2D *hptD0vsdl = new TH2D("hptD0vsdl","hptD0vsdl;pD0 (GeV/c); dl (mm)",100.,0.,20.,100,0.,5.0); 
  int count = -1;
  const double eps = 1e-6; // 1 nm in mm
  TDatabasePDG *pdgDB = TDatabasePDG::Instance();
  while (myReader.Next()) 
  {
   count++;	
 //  if (count!=nevent) continue;
 //  printf("Event No. = %d\n",count);
   hEventStat->Fill(0.5);	
    // find MC primary vertex
    int nMCPart = mcPartMass.GetSize();

    TVector3 vertex_mc(-999., -999., -999.);
    TVector3 SV_mc(-999., -999., -999.);
    TVector3 pD0(0.,0.,0.);
   for(int imc=0; imc<nMCPart; imc++)
	  {
	  if(mcPartGenStatus[imc] == 4 && mcPartPdg[imc] == 11)
	    {
	      vertex_mc.SetXYZ(mcEndPointX[imc], mcEndPointY[imc], mcEndPointZ[imc]); // collision vertex
	      break;
	    }
	   }
	   
	for(int imc=0; imc<nMCPart; imc++)
	{
	  if(fabs(mcPartPdg[imc]) != 421) continue;
       int parent_index_begin = mcPartParent_index[mcPartParent_begin[imc]];
       int parent_index_end   = mcPartParent_index[mcPartParent_end[imc]];
	  hEventStat->Fill(1.5);
	  TVector3 prodVD0 (mcPartVx[imc],mcPartVy[imc],mcPartVz[imc]); 
	  double dist = (prodVD0-vertex_mc).Mag(); 
	  if (dist > eps) hEventStat->Fill(2.5);
	  else hEventStat->Fill(3.5);
	  //printf("MCParticle = %d, PDG = %d, Dist = %f \n",imc,mcPartPdg[imc],dist); 
	  // classify origin: 4 = charm, 5 = beauty
       int originHF = GetD0OriginFlavor(imc,
                                     mcPartPdg,
                                     mcPartParent_begin,
                                     mcPartParent_end,
                                     mcPartParent_index);
         // printf("MC D0 (index %d): PDG=%d → originHF = %d\n",imc, mcPartPdg[imc], originHF);
       if(fabs(originHF)==5) hEventStat->Fill(4.5);
       else if(fabs(originHF)==4) hEventStat->Fill(5.5);          
      }	   
	 }
      c->cd();
      hEventStat->Draw("hist-text");
      c->SaveAs("hist_origin_D0meson.png");
	 fout->cd();
	 hEventStat->Write();
      fout->Close();
}
// return 5 if any beauty ancestor is found, 4 if (at least) a charm ancestor but no beauty,
int GetD0OriginFlavor(int imc,
                      const TTreeReaderArray<int> &mcPartPdg,
                      const TTreeReaderArray<unsigned int> &mcPartParent_begin,
                      const TTreeReaderArray<unsigned int> &mcPartParent_end,
                      const TTreeReaderArray<int> &mcPartParent_index)
{
  int origin = 0;        // 0 = unknown, 4 = charm, 5 = beauty
  int current = imc;
  const int maxDepth = 20;  // safety

  for (int depth = 0; depth < maxDepth; ++depth) {

    int pdg = TMath::Abs(mcPartPdg[current]);

    // charm quark or charm hadron
    if (pdg == 4 || (pdg/100 == 4) || (pdg/1000 == 4)) {
      if (origin == 0) origin = 4;
    }

    // beauty quark or beauty hadron
    if (pdg == 5 || (pdg/100 == 5) || (pdg/1000 == 5)) {
      origin = 5;
      break;
    }

    int firstParentIdx = mcPartParent_begin[current];
    int lastParentIdx  = mcPartParent_end[current];

    if (firstParentIdx >= lastParentIdx) {
      break;
    }

    int parentParticleIndex = mcPartParent_index[firstParentIdx];
    if (parentParticleIndex < 0) break;

    current = parentParticleIndex;
  }

  return origin;
}
