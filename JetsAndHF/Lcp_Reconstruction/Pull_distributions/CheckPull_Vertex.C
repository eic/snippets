// Code to draw the pull distribution of the primary vertex
// Shyam Kumar; INFN Bari; shyam.kumar@ba.infn.it; shyam055119@gmail.com

void CheckPull_Vertex(TString filelist = "test_D0.list"){
	
	gStyle->SetPalette(kRainBow);
	gStyle->SetTitleSize(0.045,"XY");	
	gStyle->SetTitleSize(0.04,"XY");	
	gStyle->SetLabelSize(0.04,"XY");	
	gStyle->SetTitleOffset(1.0,"XY");	
	gStyle->SetOptStat(1);
	gStyle->SetOptFit(1);
	gStyle->SetOptTitle(0);
	gStyle->SetGridColor(kBlack);     
	gStyle->SetGridWidth(2);        
	gStyle->SetGridStyle(2);	
	
	TH1F *hPullVtxX = new TH1F("hPullVtxX", "Pull x distribution of Primary vertex;(Vx_{rec}-Vx_{mc})/#sigma_{vx}; Entries (a.u.)", 1000, -10., 10.);
     TH1F *hPullVtxY = new TH1F("hPullVtxY", "Pull y distribution of Primary vertex;(Vy_{rec}-Vy_{mc})/#sigma_{vy}; Entries (a.u.)", 1000, -10., 10.);
     TH1F *hPullVtxZ = new TH1F("hPullVtxZ", "Pull z distribution of Primary vertex;(Vz_{rec}-Vz_{mc})/#sigma_{vz}; Entries (a.u.)", 2000, -20, 20);
       
     TH1F *hPullVtxX_1 = new TH1F("hPullVtxX_1", "Pull x distribution (Ntracks = 1) of Primary vertex;(Vx_{rec}-Vx_{mc})/#sigma_{vx}; Entries (a.u.)", 1000, -10., 10.);
     TH1F *hPullVtxY_1 = new TH1F("hPullVtxY_1", "Pull y distribution (Ntracks = 1) of Primary vertex;(Vy_{rec}-Vy_{mc})/#sigma_{vy}; Entries (a.u.)", 1000, -10., 10.);
     TH1F *hPullVtxZ_1 = new TH1F("hPullVtxZ_1", "Pull z distribution (Ntracks = 1) of Primary vertex;(Vz_{rec}-Vz_{mc})/#sigma_{vz}; Entries (a.u.)", 2000, -20, 20);
     
     TH1F *hPullVtxX_2 = new TH1F("hPullVtxX_2", "Pull x distribution (Ntracks > 1) of Primary vertex;(Vx_{rec}-Vx_{mc})/#sigma_{vx}; Entries (a.u.)", 1000, -10., 10.);
     TH1F *hPullVtxY_2 = new TH1F("hPullVtxY_2", "Pull y distribution (Ntracks > 1) of Primary vertex;(Vy_{rec}-Vy_{mc})/#sigma_{vy}; Entries (a.u.)", 1000, -10., 10.);
     TH1F *hPullVtxZ_2 = new TH1F("hPullVtxZ_2", "Pull z distribution (Ntracks > 1) of Primary vertex;(Vz_{rec}-Vz_{mc})/#sigma_{vz}; Entries (a.u.)", 2000, -20, 20);
     
     TH1F *hResVtxX = new TH1F("hResVtxX", "Residual x distribution (Ntracks > 1) of Primary vertex;(Vx_{rec}-Vx_{mc}) (mm); Entries (a.u.)", 1000, -10., 10.);
     TH1F *hResVtxY = new TH1F("hResVtxY", "Residual y distribution (Ntracks > 1) of Primary vertex;(Vy_{rec}-Vy_{mc}) (mm); Entries (a.u.)", 1000, -10., 10.);
     TH1F *hResVtxZ = new TH1F("hResVtxZ", "Residual z distribution (Ntracks > 1) of Primary vertex;(Vz_{rec}-Vz_{mc}) (mm); Entries (a.u.)", 2000, -5., 5.);
     TH1F *hChi2 = new TH1F("hChi2", "Chi2 distribution (Ntracks > 1) of Primary vertex; #chi^{2}; Entries (a.u.)", 2000, 0., 200.);
     
     TFile *fout = new TFile(Form("Pull_distribution_vertex_%s.root",filelist.Data()),"recreate");
     
	
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
  	TTreeReaderArray<int> mcPartGenStatus = {treereader, "MCParticles.generatorStatus"};
     TTreeReaderArray<int> mcPartPdg = {treereader, "MCParticles.PDG"};
     TTreeReaderArray<double> mcPartMass = {treereader, "MCParticles.mass"};
	TTreeReaderArray<double> mcEndPointX = {treereader, "MCParticles.endpoint.x"};
	TTreeReaderArray<double> mcEndPointY = {treereader, "MCParticles.endpoint.y"};
	TTreeReaderArray<double> mcEndPointZ = {treereader, "MCParticles.endpoint.z"};
	
	TTreeReaderArray<float> CTVx = {treereader, "CentralTrackVertices.position.x"};
	TTreeReaderArray<float> CTVy = {treereader, "CentralTrackVertices.position.y"};
	TTreeReaderArray<float> CTVz = {treereader, "CentralTrackVertices.position.z"};
	TTreeReaderArray<float> CTVerr_xx = {treereader, "CentralTrackVertices.positionError.xx"};
	TTreeReaderArray<float> CTVerr_yy = {treereader, "CentralTrackVertices.positionError.yy"};
	TTreeReaderArray<float> CTVerr_zz = {treereader, "CentralTrackVertices.positionError.zz"};
	TTreeReaderArray<float> CTVchi2 = {treereader, "CentralTrackVertices.chi2"};
	TTreeReaderArray<int> prim_vtx_index = {treereader, "PrimaryVertices_objIdx.index"};
     TTreeReaderArray<float> rcCharge = {treereader, "ReconstructedChargedParticles.charge"};
	
	int nevents = 0;
     while(treereader.Next())
    {
      if(nevents%1000==0) printf("\n[i] Events %d\n",nevents);
      // find MC primary vertex
      int nMCPart = mcPartMass.GetSize();

      TVector3 vertex_mc(-999., -999., -999.);
      for(int imc=0; imc<nMCPart; imc++)
	{
	  if(mcPartGenStatus[imc] == 4 && mcPartPdg[imc] == 11)
	    {
	      vertex_mc.SetXYZ(mcEndPointX[imc], mcEndPointY[imc], mcEndPointZ[imc]);
	      break;
	    }
	}
      // get RC primary vertex and it's error
     TVector3 vertex_rc(-999., -999., -999.);
     TVector3 err_vertex_rc(-999., -999., -999.);
     double chi2 = 0.;
     if(prim_vtx_index.GetSize()>0)
	{
	  int rc_vtx_index = prim_vtx_index[0];
	  vertex_rc.SetXYZ(CTVx[rc_vtx_index], CTVy[rc_vtx_index], CTVz[rc_vtx_index]);
	  err_vertex_rc.SetXYZ(sqrt(CTVerr_xx[rc_vtx_index]), sqrt(CTVerr_yy[rc_vtx_index]), sqrt(CTVerr_zz[rc_vtx_index]));
	  chi2 = CTVchi2[rc_vtx_index];
	}
	hPullVtxX->Fill((vertex_rc.x()-vertex_mc.x())/err_vertex_rc.x()); 
	hPullVtxY->Fill((vertex_rc.y()-vertex_mc.y())/err_vertex_rc.y()); 
     hPullVtxZ->Fill((vertex_rc.z()-vertex_mc.z())/err_vertex_rc.z());
     
     // Separate two cases one with one track while other more than one track  
     if (rcCharge.GetSize()==1){
	hPullVtxX_1->Fill((vertex_rc.x()-vertex_mc.x())/err_vertex_rc.x()); 
	hPullVtxY_1->Fill((vertex_rc.y()-vertex_mc.y())/err_vertex_rc.y()); 
     hPullVtxZ_1->Fill((vertex_rc.z()-vertex_mc.z())/err_vertex_rc.z());
     }
     else{
	hPullVtxX_2->Fill((vertex_rc.x()-vertex_mc.x())/err_vertex_rc.x()); 
	hPullVtxY_2->Fill((vertex_rc.y()-vertex_mc.y())/err_vertex_rc.y()); 
     hPullVtxZ_2->Fill((vertex_rc.z()-vertex_mc.z())/err_vertex_rc.z());
     
     hResVtxX->Fill((vertex_rc.x()-vertex_mc.x())); 
	hResVtxY->Fill((vertex_rc.y()-vertex_mc.y())); 
     hResVtxZ->Fill((vertex_rc.z()-vertex_mc.z()));
     hChi2->Fill(chi2);
     	
     }
           
     nevents++;
	}
	
	
	fout->cd();
	hPullVtxX->Write();
	hPullVtxY->Write();
	hPullVtxZ->Write();
	hPullVtxX_1->Write();
	hPullVtxY_1->Write();
	hPullVtxZ_1->Write();
	hPullVtxX_2->Write();
	hPullVtxY_2->Write();
	hPullVtxZ_2->Write();	
	hResVtxX->Write();
	hResVtxY->Write();
	hResVtxZ->Write();	
	hChi2->Write();	
	fout->Close();
	delete inputstream;
	delete chain;
}
