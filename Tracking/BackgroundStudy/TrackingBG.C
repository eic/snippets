void TrackingBG(TString type="Forced"){

    //Define Style
    gStyle->SetOptStat(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetLabelSize(0.035,"X");
    gStyle->SetLabelSize(0.035,"Y");
    //gStyle->SetLabelOffset(0.01,"X");
    //gStyle->SetLabelOffset(0.01,"Y");
    gStyle->SetTitleXSize(0.04);
    gStyle->SetTitleXOffset(0.9);
    gStyle->SetTitleYSize(0.04);
    gStyle->SetTitleYOffset(0.9);

    TFile *fout = new TFile("TrackingBG_10x275_Forced.root","RECREATE");

    //Define histograms
    fout->cd();
    //--
    //Eta distributions for generated particles
    TH1 * hMC1 = new TH1D("hMC1","",100,-15,10); //DIS
    hMC1->SetLineColor(kBlack); hMC1->SetLineWidth(2);
    hMC1->GetXaxis()->SetTitle("Particle #eta");hMC1->GetXaxis()->CenterTitle();

    TH1 * hMC2 = new TH1D("hMC2","",100,-15,10); //SR
    hMC2->SetLineColor(kBlue); hMC2->SetLineWidth(2);
    hMC2->GetXaxis()->SetTitle("Particle #eta");hMC2->GetXaxis()->CenterTitle();

    TH1 * hMC3 = new TH1D("hMC3","",100,-15,10); //Bremstrahlung
    hMC3->SetLineColor(kOrange); hMC3->SetLineWidth(2);
    hMC3->GetXaxis()->SetTitle("Particle #eta");hMC3->GetXaxis()->CenterTitle();

    TH1 * hMC4 = new TH1D("hMC4","",100,-15,10); //Coulomb
    hMC4->SetLineColor(kGreen); hMC4->SetLineWidth(2);
    hMC4->GetXaxis()->SetTitle("Particle #eta");hMC4->GetXaxis()->CenterTitle();

    TH1 * hMC5 = new TH1D("hMC5","",100,-15,10); //Touschek
    hMC5->SetLineColor(kRed); hMC5->SetLineWidth(2);
    hMC5->GetXaxis()->SetTitle("Particle #eta");hMC5->GetXaxis()->CenterTitle();

    TH1 * hMC6 = new TH1D("hMC6","",100,-15,10); //Proton beam gas
    hMC6->SetLineColor(kMagenta); hMC6->SetLineWidth(2);
    hMC6->GetXaxis()->SetTitle("Particle #eta");hMC6->GetXaxis()->CenterTitle();
    //--
    //Digitized hit rate on SVT Disks
    const int N_disks = 10;
    double disk_min[N_disks] = {-1040.0, -860.0, -660.0, -460.0, -260.0, 240.0, 440.0, 690.0, 990.0 , 1340.0};
    double disk_max[N_disks] = {-1000.0, -840.0, -640.0, -440.0, -240.0, 260.0, 460.0, 710.0, 1010.0, 1360.0};
    TString disk_name[N_disks] = {"E-Si Disk 4","E-Si Disk 3","E-Si Disk 2","E-Si Disk 1","E-Si Disk 0",
                                 "H-Si Disk 0","H-Si Disk 1","H-Si Disk 2","H-Si Disk 3","H-Si Disk 4"};

    TH2 *hSVTDiskRate[N_disks];
    
    //Loop over momentum points
    for(int idisk = 0; idisk < N_disks; idisk++){

        TString hname = TString::Format("hSVTDiskRate_%d", idisk);
        TString htitle = TString::Format("Digitized hit Rate per RSU per 1 ms: %s",disk_name[idisk].Data());
        
        hSVTDiskRate[idisk] = new TH2D(hname,htitle,50,-500,500,50,-500,500);
        hSVTDiskRate[idisk]->GetXaxis()->SetTitle("X [mm]");hSVTDiskRate[idisk]->GetXaxis()->CenterTitle();
        hSVTDiskRate[idisk]->GetYaxis()->SetTitle("Y [mm]");hSVTDiskRate[idisk]->GetYaxis()->CenterTitle();
    
    }
    //--
    //Digitized hit rate on SVT Vertexing layers
    const int N_vtx = 3;
    double vtx_min[N_vtx] = {30.,46.,115.};
    double vtx_max[N_vtx] = {42.,60.,130.};

    TString vtx_name[N_vtx] = {"SVT L0","SVT L1","SVT L2"};

    TH2 *hSVTVtxRate[N_vtx];
    
    //Loop over momentum points
    for(int ivtx = 0; ivtx < N_vtx; ivtx++){

        TString hname = TString::Format("hSVTVtxRate_%d", ivtx);
        TString htitle = TString::Format("Digitized hit Rate per RSU per 1 ms: %s",vtx_name[ivtx].Data());
        
        hSVTVtxRate[ivtx] = new TH2D(hname,htitle,20,-200,200,40,-400,400);
        hSVTVtxRate[ivtx]->GetXaxis()->SetTitle("Z [mm]");hSVTVtxRate[ivtx]->GetXaxis()->CenterTitle();
        hSVTVtxRate[ivtx]->GetYaxis()->SetTitle("r #phi [mm]");hSVTVtxRate[ivtx]->GetYaxis()->CenterTitle();
    
    }
    //--
    //Digitized hit rate on SVT Barrel layers
    const int N_barrel = 2;
    double barrel_min[N_barrel] = {260.,410.};
    double barrel_max[N_barrel] = {280.,450.};

    TString barrel_name[N_barrel] = {"SVT L3","SVT L4"};

    TH2 *hSVTbarrelRate[N_barrel];
    
    //Loop over momentum points
    for(int ibarrel = 0; ibarrel < N_barrel; ibarrel++){

        TString hname = TString::Format("hSVTbarrelRate_%d", ibarrel);
        TString htitle = TString::Format("Digitized hit Rate per RSU per 1 ms: %s", barrel_name[ibarrel].Data());
        
        hSVTbarrelRate[ibarrel] = new TH2D(hname,htitle,40,-400,400,150,-1500,1500);
        hSVTbarrelRate[ibarrel]->GetXaxis()->SetTitle("Z [mm]");hSVTbarrelRate[ibarrel]->GetXaxis()->CenterTitle();
        hSVTbarrelRate[ibarrel]->GetYaxis()->SetTitle("r #phi [mm]");hSVTbarrelRate[ibarrel]->GetYaxis()->CenterTitle();
    
    }
    //--

    //Open file
    TString input = "10x275_" + type + "/eicrecon_*.root";

    TChain *tree = new TChain("events");
    tree->Add(input.Data()); //Events with backgrounds present

    cout<<"Running analysis on "<<input<<"!"<<endl;
    cout<<"Analyzing "<<tree->GetEntries()<<" events!"<<endl;

    //Create Array Reader
    TTreeReader tr(tree);

    TTreeReaderArray<int>   gen_status(tr, "MCParticles.generatorStatus");
    TTreeReaderArray<int>   gen_pid(tr, "MCParticles.PDG");
    TTreeReaderArray<double> gen_px(tr, "MCParticles.momentum.x");
    TTreeReaderArray<double> gen_py(tr, "MCParticles.momentum.y");
    TTreeReaderArray<double> gen_pz(tr, "MCParticles.momentum.z");
    TTreeReaderArray<double> gen_mass(tr, "MCParticles.mass");
    TTreeReaderArray<float> gen_charge(tr, "MCParticles.charge");
    TTreeReaderArray<double> gen_vx(tr, "MCParticles.vertex.x");
    TTreeReaderArray<double> gen_vy(tr, "MCParticles.vertex.y");
    TTreeReaderArray<double> gen_vz(tr, "MCParticles.vertex.z");
    TTreeReaderArray<float> SiEC_RecHit_posx(tr,"SiEndcapTrackerRecHits.position.x");
    TTreeReaderArray<float> SiEC_RecHit_posy(tr,"SiEndcapTrackerRecHits.position.y");
    TTreeReaderArray<float> SiEC_RecHit_posz(tr,"SiEndcapTrackerRecHits.position.z");
    TTreeReaderArray<float> SiVtx_RecHit_posx(tr,"SiBarrelVertexRecHits.position.x");
    TTreeReaderArray<float> SiVtx_RecHit_posy(tr,"SiBarrelVertexRecHits.position.y");
    TTreeReaderArray<float> SiVtx_RecHit_posz(tr,"SiBarrelVertexRecHits.position.z");
    TTreeReaderArray<float> Sibarrel_RecHit_posx(tr,"SiBarrelTrackerRecHits.position.x");
    TTreeReaderArray<float> Sibarrel_RecHit_posy(tr,"SiBarrelTrackerRecHits.position.y");
    TTreeReaderArray<float> Sibarrel_RecHit_posz(tr,"SiBarrelTrackerRecHits.position.z");

    //Define other variables
    int counter(0);
    TLorentzVector gen_vec;

    //Loop over events
    while (tr.Next()) {

	    if(counter%100==0) cout<<"Analyzing event "<<counter<<endl;
	    counter++;

        //Loop over generated particles
        for(int igen=0;igen<gen_status.GetSize();igen++){
	        
            auto status = gen_status[igen];
            gen_vec.SetXYZM(gen_px[igen],gen_py[igen],gen_pz[igen],gen_mass[igen]);
            auto eta = gen_vec.Eta();

            if( status==6001 )
                hMC6->Fill(eta);
            else if( status==5001 )
                hMC5->Fill(eta);
            else if( status==4001 )
                hMC4->Fill(eta);
            else if( status==3001 )
                hMC3->Fill(eta);
            else if( status==2001 )
                hMC2->Fill(eta);
            else if( status==1 )
                hMC1->Fill(eta);

        } //End loop over generated particles

        //Loop over Digitized hits for SVT disks
        for(int ihit=0;ihit < SiEC_RecHit_posx.GetSize();ihit++){

            for(int idisk = 0; idisk < N_disks; idisk++){
                if( SiEC_RecHit_posz[ihit] > disk_min[idisk] && SiEC_RecHit_posz[ihit] < disk_max[idisk] )
                    hSVTDiskRate[idisk]->Fill(SiEC_RecHit_posx[ihit],SiEC_RecHit_posy[ihit]);
            }

        } // End Loop over Digitized hits for SVT disks

        //Loop over Digitized hits for SVT Vtx layers
        for(int ihit=0;ihit < SiVtx_RecHit_posx.GetSize();ihit++){

            for(int ivtx = 0; ivtx < N_vtx; ivtx++){

                auto hit_r = std::hypot(SiVtx_RecHit_posx[ihit],SiVtx_RecHit_posy[ihit]);
                auto hit_phi = std::atan2(SiVtx_RecHit_posy[ihit],SiVtx_RecHit_posx[ihit]);
                auto hit_z = SiVtx_RecHit_posz[ihit];

                if( hit_r > vtx_min[ivtx] && hit_r < vtx_max[ivtx] )
                    hSVTVtxRate[ivtx]->Fill(hit_z,hit_r*hit_phi);
            }

        } // End Loop over Digitized hits for SVT Vtx layers

        //Loop over Digitized hits for SVT Barrel layers
        for(int ihit=0;ihit < Sibarrel_RecHit_posx.GetSize();ihit++){

            for(int ibarrel = 0; ibarrel < N_barrel; ibarrel++){

                auto hit_r = std::hypot(Sibarrel_RecHit_posx[ihit],Sibarrel_RecHit_posy[ihit]);
                auto hit_phi = std::atan2(Sibarrel_RecHit_posy[ihit],Sibarrel_RecHit_posx[ihit]);
                auto hit_z = Sibarrel_RecHit_posz[ihit];

                if( hit_r > barrel_min[ibarrel] && hit_r < barrel_max[ibarrel] )
                    hSVTbarrelRate[ibarrel]->Fill(hit_z,hit_r*hit_phi);
            }

        } // End Loop over Digitized hits for SVT Barrel layers

    } //End loop over events

    //Check SR rate
    cout <<"------"<<endl;
    cout <<"Total number of SR photons generated: "<< hMC2->GetEntries() << endl;
    cout <<"SR photon rate [MHz]: " << 1.0e-6 * ( hMC2->GetEntries() / (2.0e-6 * tree->GetEntries()) ) << endl;
    cout <<"------"<<endl;

    //Scale hit rate histograms to 1 ms.
    //Each 'event' is 2us.
    auto SFac = 1.0e-3 / (2.0e-6 * tree->GetEntries());

    for(int idisk = 0; idisk < N_disks; idisk++){
        hSVTDiskRate[idisk]->Scale(SFac);
    }

    for(int ivtx = 0; ivtx < N_vtx; ivtx++){
        hSVTVtxRate[ivtx]->Scale(SFac);
    }

    for(int ibarrel = 0; ibarrel < N_barrel; ibarrel++){
        hSVTbarrelRate[ibarrel]->Scale(SFac);
    }

    //Write rate histograms to output ROOT file
    fout->cd();
    for(int idisk = 0; idisk < N_disks; idisk++){
        hSVTDiskRate[idisk]->Write();
    }

    for(int ivtx = 0; ivtx < N_vtx; ivtx++){
        hSVTVtxRate[ivtx]->Write();
    }

    for(int ibarrel = 0; ibarrel < N_barrel; ibarrel++){
        hSVTbarrelRate[ibarrel]->Write();
    }

    // Make plots
    TCanvas *c1 = new TCanvas("c1");
    
    auto frame1 = c1->DrawFrame(-15,0.1,10,1e8);
    frame1->SetTitle("10x275 GeV: Forced DIS Configuration");
    frame1->GetXaxis()->SetTitle("Particle #eta");frame1->GetXaxis()->CenterTitle();
    gPad->SetLogy();

    hMC1->Draw("same");
    hMC2->Draw("same");
    hMC3->Draw("same");
    hMC4->Draw("same");
    hMC5->Draw("same");
    hMC6->Draw("same");

    TLegend *leg1 = new TLegend(0.175, 0.65, 0.375, 0.85, "MCParticles.generatorStatus");
    leg1->AddEntry(hMC1, "1: DIS", "l");
    leg1->AddEntry(hMC2, "2001: SR", "l");
    leg1->AddEntry(hMC3, "3001: Bremstrahlung", "l");
    leg1->AddEntry(hMC4, "4001: Coulomb", "l");
    leg1->AddEntry(hMC5, "5001: Touschek", "l");
    leg1->AddEntry(hMC6, "6001: Proton beam gas", "l");
    leg1->Draw();

    c1->Print("plots/TrackingBG_10x275_Forced.pdf[");
    c1->Print("plots/TrackingBG_10x275_Forced.pdf");

    TCanvas *c2[N_disks];
    for(int idisk = 0; idisk < N_disks; idisk++){
        c2[idisk] = new TCanvas(Form("c2[%d]",idisk));
        hSVTDiskRate[idisk]->Draw("colz");
        c2[idisk]->Print("plots/TrackingBG_10x275_Forced.pdf");
        
    }

    TCanvas *c3[N_vtx];
    for(int ivtx = 0; ivtx < N_vtx; ivtx++){
        c3[ivtx] = new TCanvas(Form("c3[%d]",ivtx));
        hSVTVtxRate[ivtx]->Draw("colz");
        c3[ivtx]->Print("plots/TrackingBG_10x275_Forced.pdf");
        
    }

    TCanvas *c4[N_barrel];
    for(int ibarrel = 0; ibarrel < N_barrel; ibarrel++){
        c4[ibarrel] = new TCanvas(Form("c4[%d]",ibarrel));
        hSVTbarrelRate[ibarrel]->Draw("colz");
        c4[ibarrel]->Print("plots/TrackingBG_10x275_Forced.pdf");
        
    }

    c4[N_barrel-1]->Print("plots/TrackingBG_10x275_Forced.pdf]");

    // Close output ROOT file
    fout->Close();

}