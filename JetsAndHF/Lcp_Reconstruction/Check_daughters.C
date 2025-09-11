void Check_daughters(TString filelist = "test_Lc.list")
{
    
    
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
   TTreeReader treereader(chain);

  // MC  
  TTreeReaderArray<int> mcPartGenStatus = {treereader, "MCParticles.generatorStatus"};
  TTreeReaderArray<int> mcPartPdg = {treereader, "MCParticles.PDG"};
  TTreeReaderArray<float> mcPartCharge = {treereader, "MCParticles.charge"};
  TTreeReaderArray<unsigned int> mcPartParent_begin = {treereader, "MCParticles.parents_begin"};
  TTreeReaderArray<unsigned int> mcPartParent_end = {treereader, "MCParticles.parents_end"};
  TTreeReaderArray<int> mcPartParent_index = {treereader, "_MCParticles_parents.index"};
  TTreeReaderArray<unsigned int> mcPartDaughter_begin = {treereader, "MCParticles.daughters_begin"};
  TTreeReaderArray<unsigned int> mcPartDaughter_end = {treereader, "MCParticles.daughters_end"};
  TTreeReaderArray<int> mcPartDaughter_index = {treereader, "_MCParticles_daughters.index"};
  TTreeReaderArray<double> mcPartMass = {treereader, "MCParticles.mass"};
  TTreeReaderArray<double> mcPartVx = {treereader, "MCParticles.vertex.x"};
  TTreeReaderArray<double> mcPartVy = {treereader, "MCParticles.vertex.y"};
  TTreeReaderArray<double> mcPartVz = {treereader, "MCParticles.vertex.z"};
  TTreeReaderArray<double> mcMomPx = {treereader, "MCParticles.momentum.x"};
  TTreeReaderArray<double> mcMomPy = {treereader, "MCParticles.momentum.y"};
  TTreeReaderArray<double> mcMomPz = {treereader, "MCParticles.momentum.z"};
  TTreeReaderArray<double> mcEndPointX = {treereader, "MCParticles.endpoint.x"};
  TTreeReaderArray<double> mcEndPointY = {treereader, "MCParticles.endpoint.y"};
  TTreeReaderArray<double> mcEndPointZ = {treereader, "MCParticles.endpoint.z"};	
  
    TDatabasePDG *pdgDB = TDatabasePDG::Instance();
    Int_t iEvent = 0;
    while(treereader.Next())
    {
   //	printf("Event No. = %d \n", iEvent);
    	int nMCPart = mcPartMass.GetSize();
    	for(int imc=0; imc<nMCPart; imc++)
	{
	  if(fabs(mcPartPdg[imc]) != 4122) continue;
	  
	  TString particle_name = pdgDB->GetParticle(mcPartPdg[imc])->GetName();
 
	  int nDaughters = mcPartDaughter_end[imc]-mcPartDaughter_begin[imc];
	  if(nDaughters!=3) continue;
	  int daug_index_1 = mcPartDaughter_index[mcPartDaughter_begin[imc]];
	  int daug_index_2 = mcPartDaughter_index[mcPartDaughter_begin[imc]+1];
	  int daug_index_3 = mcPartDaughter_index[mcPartDaughter_begin[imc]+2];
	 // int daug_index_4 = mcPartDaughter_index[mcPartDaughter_begin[imc]+3];
	  int daug_pdg_1 = mcPartPdg[daug_index_1];
	  int daug_pdg_2 = mcPartPdg[daug_index_2];
	  int daug_pdg_3 = mcPartPdg[daug_index_3]; 
	 // int daug_pdg_4 = mcPartPdg[daug_index_4];
	  if (fabs(daug_pdg_1)!=13 && fabs(daug_pdg_2)!=13 && fabs(daug_pdg_3)!=13) continue;
	  printf("Event No. = %d \n", iEvent);
	  printf("Parent Particle = %s \n", particle_name.Data());
	  printf("PDG name = %d \t %d \t %d \t \n",daug_pdg_1, daug_pdg_2, daug_pdg_3);
	  printf("Daughters name = %s \t %s \t %s \t \n", pdgDB->GetParticle(daug_pdg_1)->GetName(),pdgDB->GetParticle(daug_pdg_2)->GetName(), pdgDB->GetParticle(daug_pdg_3)->GetName());

    	
    	
    	
    }
  
	iEvent++;
	
	}
	delete inputstream;
	delete chain;
}
