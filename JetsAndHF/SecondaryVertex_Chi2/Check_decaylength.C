// Code to check the decaylength of D0 at the simulation level.
// Shyam Kumar; shyam.kumar@ba.infn.it; shyam055119@gmail.com

void Check_decaylength(TString filelist = "test.list"){

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

  TFile *fout = new TFile("output_file.root","recreate");
  TH2D *hptD0vsdl = new TH2D("hptD0vsdl","hptD0vsdl;pD0 (GeV/c); dl (mm)",100.,0.,20.,100,0.,5.0); 
  int count = -1;
  TDatabasePDG *pdgDB = TDatabasePDG::Instance();
  while (myReader.Next()) 
  {
   count++;	
 //  if (count!=nevent) continue;
   printf("Event No. = %d\n",count);	
    // find MC primary vertex
    int nMCPart = mcPartMass.GetSize();

    TVector3 vertex_mc(-999., -999., -999.);
    TVector3 SV_mc(-999., -999., -999.);
    TVector3 pD0(0.,0.,0.);
   for(int imc=0; imc<nMCPart; imc++)
	  {
	  if(mcPartGenStatus[imc] == 4 && mcPartPdg[imc] == 11)
	    {
	      vertex_mc.SetXYZ(mcEndPointX[imc], mcEndPointY[imc], mcEndPointZ[imc]);
	      break;
	    }
	   }
	   
	for(int imc=0; imc<nMCPart; imc++)
	{
	  if(fabs(mcPartPdg[imc]) != 421) continue;
	  pD0.SetXYZ(mcMomPx[imc],mcMomPy[imc],mcMomPz[imc]);
	  int nDaughters = mcPartDaughter_end[imc]-mcPartDaughter_begin[imc];
	  if(nDaughters!=2) continue;
	  // A "status 21" particle might represent the before-hadronization D⁰.
       //A "status 1" particle is the final, stable D⁰ after decays or transport.
	  int parent_index_begin = mcPartParent_index[mcPartParent_begin[imc]];
	  int parent_index_end = mcPartParent_index[mcPartParent_end[imc]];
	  TString particle1 = pdgDB->GetParticle(mcPartPdg[parent_index_begin])->GetName();
	  TString particle2 = pdgDB->GetParticle(mcPartPdg[parent_index_end])->GetName();
    printf("Parent Name  = (%s, %s) \n",particle1.Data(),particle2.Data());
	  // find D0 that decay into pi+K
	  bool is_pik_decay = false;	  
	  int daug_index_1 = mcPartDaughter_index[mcPartDaughter_begin[imc]];
	  int daug_index_2 = mcPartDaughter_index[mcPartDaughter_begin[imc]+1];
	  int daug_pdg_1 = mcPartPdg[daug_index_1];
	  int daug_pdg_2 = mcPartPdg[daug_index_2];
	  if( (fabs(daug_pdg_1)==321 && fabs(daug_pdg_2)==211) || (fabs(daug_pdg_1)==211 && fabs(daug_pdg_2)==321) )
	  {
	      is_pik_decay = true;
	      SV_mc.SetXYZ(mcPartVx[daug_index_1],mcPartVy[daug_index_1],mcPartVz[daug_index_1]);
	      SV_mc.Print();
	      hptD0vsdl->Fill(pD0.Mag(),(SV_mc-vertex_mc).Mag());
	  }
	  
    }	   
	 
	 }
	 fout->cd();
   hptD0vsdl->Write();
   fout->Close();

}
