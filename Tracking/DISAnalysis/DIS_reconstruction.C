//------------------
// energy_set 0 is 18x275 GeV
// energy_set 1 is 10x100 GeV
void DIS_reconstruction(int energy_set = 0, int Q2_set = 1, bool use_campaign = 0){

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

    //Define histograms
 
    //---Eta distributions---
    TH1 *h1a = new TH1D("h1a","Generated Charged Particles",100,-4,4);
    h1a->GetXaxis()->SetTitle("#eta_{gen.}");h1a->GetXaxis()->CenterTitle();
    h1a->SetLineColor(kTeal);h1a->SetLineWidth(2);

    TH1 *h1a1 = new TH1D("h1a1","Generated Charged Particles",100,-4,4);  //Minimum momentum cut of Pt > 200 MeV/c
    h1a1->GetXaxis()->SetTitle("#eta_{gen.}");h1a1->GetXaxis()->CenterTitle();
    h1a1->SetLineColor(kRed);h1a1->SetLineWidth(2);

    TH1 *h1a2 = new TH1D("h1a2","Generated Charged Particles",100,-4,4);  //Minimum momentum cut of Pt > 500 MeV/c
    h1a2->GetXaxis()->SetTitle("#eta_{gen.}");h1a2->GetXaxis()->CenterTitle();
    h1a2->SetLineColor(kBlack);h1a2->SetLineWidth(2);
    h1a2->SetFillColor(kBlack);h1a2->SetFillStyle(3244);

    TH1 *h1b = new TH1D("h1b","Reconstructed Real-seeded tracks",100,-4,4);
    h1b->GetXaxis()->SetTitle("#eta_{rec.}");h1b->GetXaxis()->CenterTitle();
    h1b->SetLineColor(kGreen);h1b->SetLineWidth(2);

    TH1 *h1b1 = new TH1D("h1b1","Reconstructed Real-seeded tracks",100,-4,4); //Minimum momentum cut of Pt > 200 MeV/c
    h1b1->GetXaxis()->SetTitle("#eta_{rec.}");h1b1->GetXaxis()->CenterTitle();
    h1b1->SetLineColor(kBlue);h1b1->SetLineWidth(2);

    TH1 *h1b2 = new TH1D("h1b2","Reconstructed Real-seeded tracks",100,-4,4); //Minimum momentum cut of Pt > 500 MeV/c
    h1b2->GetXaxis()->SetTitle("#eta_{rec.}");h1b2->GetXaxis()->CenterTitle();
    h1b2->SetLineColor(kRed);h1b2->SetLineWidth(2);
    h1b2->SetMarkerColor(kRed);h1b2->SetMarkerStyle(kFullCrossX);//h1b2->SetMarkerSize(0.75);

    TH1 *h1c = new TH1D("h1c","Reconstructed Truth-seeded tracks",100,-4,4);
    h1c->GetXaxis()->SetTitle("#eta_{rec.}");h1c->GetXaxis()->CenterTitle();
    h1c->SetLineColor(kRed);h1c->SetLineWidth(2);

    TH1 *h1c1 = new TH1D("h1c1","Reconstructed Truth-seeded tracks",100,-4,4); //Minimum momentum cut of Pt > 200 MeV/c
    h1c1->GetXaxis()->SetTitle("#eta_{rec.}");h1c1->GetXaxis()->CenterTitle();
    h1c1->SetLineColor(kOrange);h1c1->SetLineWidth(2);

    TH1 *h1c2 = new TH1D("h1c2","Reconstructed Truth-seeded tracks",100,-4,4); //Minimum momentum cut of Pt > 500 MeV/c
    h1c2->GetXaxis()->SetTitle("#eta_{rec.}");h1c2->GetXaxis()->CenterTitle();
    h1c2->SetLineColor(kMagenta);h1c2->SetLineWidth(2);
    h1c2->SetMarkerColor(kMagenta);h1c2->SetMarkerStyle(kFullCrossX);//h1c2->SetMarkerSize(0.75);

    TH1 *h1rb1 = new TH1D("h1rb1","",100,-4,4); //Real-seeded tracks (Pt > 200 MeV/c cut)
    TH1 *h1rc1 = new TH1D("h1rc1","",100,-4,4); //Truth-seeded tracks (Pt > 200 MeV/c cut)
    TH1 *h1rb2 = new TH1D("h1rb2","",100,-4,4); //Real-seeded tracks (Pt > 500 MeV/c cut)
    TH1 *h1rc2 = new TH1D("h1rc2","",100,-4,4); //Truth-seeded tracks (Pt > 500 MeV/c cut)

    //Transverse Momentum distributions
    //For generated charged particles and tracks, place |eta| < 3.5 cut
    TH1 *h2a = new TH1D("h2a","Generated Charged Particles",100,0,15);
    h2a->GetXaxis()->SetTitle("P_{T} [GeV/c]");h2a->GetXaxis()->CenterTitle();
    h2a->SetLineColor(kBlack);h2a->SetLineWidth(2);
    h2a->SetFillColor(kBlack);h2a->SetFillStyle(3244);

    TH1 *h2b = new TH1D("h2b","Reconstructed Real-seeded tracks",100,0,15);
    h2b->GetXaxis()->SetTitle("P_{T} [GeV/c]");h2b->GetXaxis()->CenterTitle();
    h2b->SetLineColor(kRed);h2b->SetLineWidth(2);
    h2b->SetMarkerColor(kRed);h2b->SetMarkerStyle(kFullCrossX);

    TH1 *h2c = new TH1D("h2c","Reconstructed Truth-seeded tracks",100,0,15);
    h2c->GetXaxis()->SetTitle("P_{T} [GeV/c]");h2c->GetXaxis()->CenterTitle();
    h2c->SetLineColor(kMagenta);h2c->SetLineWidth(2);
    h2c->SetMarkerColor(kMagenta);h2c->SetMarkerStyle(kFullCrossX);

    //Transverse Momentum distributions for generated charged particles
    //with matching track. Only consider particles w/ |eta| < 3.5
    TH1 *h3a = new TH1D("h3a","",100,0,15);
    h3a->GetXaxis()->SetTitle("P_{T} [GeV/c]");h3a->GetXaxis()->CenterTitle();
    h3a->SetLineColor(kRed);h3a->SetLineWidth(2);
    h3a->SetMarkerColor(kRed);h3a->SetMarkerStyle(kFullCircle);    

    TH1 *h3b = new TH1D("h3b","",100,0,15);
    h3b->GetXaxis()->SetTitle("P_{T} [GeV/c]");h3b->GetXaxis()->CenterTitle();
    h3b->SetLineColor(kMagenta);h3b->SetLineWidth(2);
    h3b->SetMarkerColor(kMagenta);h3b->SetMarkerStyle(kOpenCircle);

    TH1 *h3ra = new TH1D("h3ra","",100,0,15);
    TH1 *h3rb = new TH1D("h3rb","",100,0,15);

    //File information
    TString path;
    if(use_campaign) path = "./campaign_input/";
    else path = "./input/";

    TString beam_energies;
    if(energy_set == 0) beam_energies = "18x275";
    else if(energy_set == 1) beam_energies = "10x100";

    TString run_name;
    if(use_campaign) run_name.Form("list_%s_Q2_%d.txt",beam_energies.Data(),Q2_set);
    else run_name.Form("%s/Q2_%d/eicrecon_*root",beam_energies.Data(),Q2_set);
    
    //Open File
    TString input = path + run_name;
    TChain *tree = new TChain("events");

    if(use_campaign){
	std::ifstream in(input);
	std::string file("");
	while (in >> file)
		tree->Add(file.data());
    	in.close();
    }
    else{
    	tree->Add(input.Data());
   }

    cout<<"Running file "<<input<<"!"<<endl;
    cout<<"Analyzing "<<tree->GetEntries()<<" events!"<<endl;

    //Create Array Reader
    TTreeReader tr(tree);

    TTreeReaderArray<int>   gen_status(tr, "MCParticles.generatorStatus");
    TTreeReaderArray<int>   gen_pid(tr, "MCParticles.PDG");
    TTreeReaderArray<double> gen_px(tr, "MCParticles.momentum.x");
    TTreeReaderArray<double> gen_py(tr, "MCParticles.momentum.y");
    TTreeReaderArray<double> gen_pz(tr, "MCParticles.momentum.z");
    TTreeReaderArray<double> gen_mass(tr, "MCParticles.mass"); //Not important here
    TTreeReaderArray<float> gen_charge(tr, "MCParticles.charge");
    TTreeReaderArray<double> gen_vx(tr, "MCParticles.vertex.x");
    TTreeReaderArray<double> gen_vy(tr, "MCParticles.vertex.y");
    TTreeReaderArray<double> gen_vz(tr, "MCParticles.vertex.z");

    TTreeReaderArray<float> rec_px(tr, "ReconstructedChargedParticles.momentum.x");
    TTreeReaderArray<float> rec_py(tr, "ReconstructedChargedParticles.momentum.y");
    TTreeReaderArray<float> rec_pz(tr, "ReconstructedChargedParticles.momentum.z");
    TTreeReaderArray<float> rec_mass(tr, "ReconstructedChargedParticles.mass");
    TTreeReaderArray<int> rec_type(tr, "ReconstructedChargedParticles.type"); //Type 0: Match to generated particle
    TTreeReaderArray<int> rec_pdg(tr, "ReconstructedChargedParticles.PDG"); //Uses PID lookup tables

    TTreeReaderArray<float> rec_ts_px(tr, "ReconstructedTruthSeededChargedParticles.momentum.x");
    TTreeReaderArray<float> rec_ts_py(tr, "ReconstructedTruthSeededChargedParticles.momentum.y");
    TTreeReaderArray<float> rec_ts_pz(tr, "ReconstructedTruthSeededChargedParticles.momentum.z");
    TTreeReaderArray<float> rec_ts_mass(tr, "ReconstructedTruthSeededChargedParticles.mass");
    TTreeReaderArray<int> rec_ts_type(tr, "ReconstructedTruthSeededChargedParticles.type"); //Match to generated particle
    TTreeReaderArray<int> rec_ts_pdg(tr, "ReconstructedTruthSeededChargedParticles.PDG"); //Uses PID lookup tables

    TTreeReaderArray<unsigned int> assoc_rec_id(tr, "ReconstructedChargedParticleAssociations.recID");
    TTreeReaderArray<unsigned int> assoc_sim_id(tr, "ReconstructedChargedParticleAssociations.simID");
    TTreeReaderArray<float> assoc_weight(tr, "ReconstructedChargedParticleAssociations.weight");

    TTreeReaderArray<unsigned int> assoc_ts_rec_id(tr, "ReconstructedTruthSeededChargedParticleAssociations.recID");
    TTreeReaderArray<unsigned int> assoc_ts_sim_id(tr, "ReconstructedTruthSeededChargedParticleAssociations.simID");
    TTreeReaderArray<float> assoc_ts_weight(tr, "ReconstructedTruthSeededChargedParticleAssociations.weight");

    //Other variables
    TLorentzVector gen_vec;
    TVector3 gen_vertex;

    TLorentzVector rec_vec;
    TVector3 track_vec; //Reconstructed track momentum vector
    int counter(0);

    //Loop over events
    while (tr.Next()) {

	    if(counter%100==0) cout<<"Analyzing event "<<counter<<endl;
	    counter++;

        //Loop over generated particles
        for(int igen=0;igen<gen_status.GetSize();igen++){
	        
            auto charge = gen_charge[igen];
            auto status = gen_status[igen];

            //Require final-state, charged particle (no secondaries)
            if( status==1 && fabs(charge) > 0.01 ){
                gen_vec.SetXYZM(gen_px[igen],gen_py[igen],gen_pz[igen],gen_mass[igen]);
                gen_vertex.SetXYZ(gen_vx[igen],gen_vy[igen],gen_vz[igen]);
                
		auto eta = gen_vec.Eta();
		auto pt = gen_vec.Pt();

                //Fill eta histograms
                h1a->Fill(eta);
                if( pt > 0.2 ) h1a1->Fill(eta);
                if( pt > 0.5 ) h1a2->Fill(eta);

		//Fill Pt histogram
		if( fabs(eta) < 3.5 ) h2a->Fill(pt);		

		//For |eta| < 3.5, find if there is an associated reconstructed track with weight > 0.8
		
		//Real-seeded tracks
		for(int iassoc = 0; iassoc < assoc_weight.GetSize(); iassoc++){

			if(fabs(eta) < 3.5 && assoc_sim_id[iassoc] == igen && assoc_weight[iassoc]>0.8){
				//For a given track,, only one MC Particle can have a weight > 0.8
				//---
				//It is possible that a given MC Particle is associated to multiple
				//tracks with weight > 0.8, although this is unlikely due to shared hit
				//cut in ambiguity solver
				h3a->Fill(pt);	
			}
		}

		//Truth-seeded tracks
		for(int iassoc = 0; iassoc < assoc_ts_weight.GetSize(); iassoc++){

                        if(fabs(eta) < 3.5 && assoc_ts_sim_id[iassoc] == igen && assoc_ts_weight[iassoc]>0.8){
                                //For a given track,, only one MC Particle can have a weight > 0.8
                                //---
                                //It is possible that a given MC Particle is associated to multiple
                                //tracks with weight > 0.8, although this is unlikely due to shared hit
                                //cut in ambiguity solver
                                h3b->Fill(pt);
                        }
                }

            } //if generated final-state charged particle
        } //End loop over generated particles
        
        //Loop over reconstructed real-seeded charged particles (copy of tracks with PID info)
        int rec_mult = rec_type.GetSize();

        for(int irec=0;irec<rec_mult;irec++){

            rec_vec.SetXYZM(rec_px[irec],rec_py[irec],rec_pz[irec],rec_mass[irec]);

	    auto eta = rec_vec.Eta();
            auto pt = rec_vec.Pt();            

            //Fill eta histograms
            h1b->Fill(eta);
            if( pt > 0.2 ) h1b1->Fill(eta);
            if( pt > 0.5 ) h1b2->Fill(eta);
            //if( rec_vec.Pt() > 0.2 && rec_type[irec] == 0 ) h1b2->Fill(rec_vec.Eta());
	    
	    //Fill pt histogram
	    if( fabs(eta) < 3.5 ) h2b->Fill(pt);

        } //End loop over reconstructed particles

        //Loop over reconstructed truth-seeded charged particles (copy of tracks with PID info)
        int rec_ts_mult = rec_ts_type.GetSize();

        for(int irec=0;irec<rec_ts_mult;irec++){

            rec_vec.SetXYZM(rec_ts_px[irec],rec_ts_py[irec],rec_ts_pz[irec],rec_ts_mass[irec]);
            
	    auto eta = rec_vec.Eta();
            auto pt = rec_vec.Pt();
	
            //Fill eta histograms
            h1c->Fill(rec_vec.Eta());
            if( rec_vec.Pt() > 0.2 ) h1c1->Fill(rec_vec.Eta());
            if( rec_vec.Pt() > 0.5 ) h1c2->Fill(rec_vec.Eta());
            //if( rec_vec.Pt() > 0.2 && rec_ts_type[irec] == 0 ) h1c2->Fill(rec_vec.Eta());

	    //Fill pt histogram
	    if( fabs(eta) < 3.5 ) h2c->Fill(pt);

        } //End loop over reconstructed particles

    } //End loop over events

    //Make ratio histograms
    h1rb1 = (TH1*) h1b1->Clone("h1rb1");
    h1rb1->Divide(h1a1);
    h1rb1->SetTitle("Ratio of recontructed to generated particle counts: P_{T} > 200 MeV/c");
    h1rb1->GetXaxis()->SetTitle("#eta");h1rb1->GetXaxis()->CenterTitle();
    h1rb1->GetYaxis()->SetTitle("Ratio");h1rb1->GetYaxis()->CenterTitle();

    h1rc1 = (TH1*) h1c1->Clone("h1rc1");
    h1rc1->Divide(h1a1);
    h1rc1->SetTitle("Ratio of recontructed to generated particle counts: P_{T} > 200 MeV/c");
    h1rc1->GetXaxis()->SetTitle("#eta");h1rc1->GetXaxis()->CenterTitle();
    h1rc1->GetYaxis()->SetTitle("Ratio");h1rc1->GetYaxis()->CenterTitle();

    h1rb2 = (TH1*) h1b2->Clone("h1rb2");
    h1rb2->Divide(h1a2);
    h1rb2->SetTitle("Ratio of recontructed to generated particle counts: P_{T} > 500 MeV/c");
    h1rb2->GetXaxis()->SetTitle("#eta");h1rb2->GetXaxis()->CenterTitle();
    h1rb2->GetYaxis()->SetTitle("Ratio");h1rb2->GetYaxis()->CenterTitle();

    h1rc2 = (TH1*) h1c2->Clone("h1rc2");
    h1rc2->Divide(h1a2);
    h1rc2->SetTitle("Ratio of recontructed to generated particle counts: P_{T} > 500 MeV/c");
    h1rc2->GetXaxis()->SetTitle("#eta");h1rc2->GetXaxis()->CenterTitle();
    h1rc2->GetYaxis()->SetTitle("Ratio");h1rc2->GetYaxis()->CenterTitle();

    h3ra = (TH1*) h3a->Clone("h3ra");
    h3ra->Divide(h2a);

    h3rb = (TH1*) h3b->Clone("h3rb");
    h3rb->Divide(h2a);

    //Make plots
    //Generated charged particles
    TCanvas *c1a = new TCanvas("c1a");
    h1a->Draw();
    h1a1->Draw("same");
    h1a2->Draw("same");

    TLegend *leg1a = new TLegend(0.125,0.7,0.375,0.875);
    leg1a->SetBorderSize(0);leg1a->SetFillStyle(0);
    leg1a->SetHeader(Form("Pythia8: %s GeV, Q^{2} > %d GeV^{2}",beam_energies.Data(),Q2_set));
    leg1a->AddEntry("h1a","All generated charged particles","l");
    leg1a->AddEntry("h1a1","+ P_{T} > 200 MeV/c","l");
    leg1a->AddEntry("h1a2","+ P_{T} > 500 MeV/c","l");
    leg1a->Draw();

    //Real-seeded tracks
    TCanvas *c1b = new TCanvas("c1b");
    h1b->Draw();
    h1b1->Draw("same");
    h1b2->Draw("same");

    TLegend *leg1b = new TLegend(0.25,0.7,0.5,0.875);
    leg1b->SetBorderSize(0);leg1b->SetFillStyle(0);
    leg1b->SetHeader(Form("Pythia8: %s GeV, Q^{2} > %d GeV^{2}",beam_energies.Data(),Q2_set));
    leg1b->AddEntry("h1b","All real-seeded tracks","l");
    leg1b->AddEntry("h1b1","+ P_{T} > 200 MeV/c","l");
    leg1b->AddEntry("h1b2","+ P_{T} > 500 MeV/c","l");
    leg1b->Draw();

    //Truth-seeded tracks
    TCanvas *c1c = new TCanvas("c1c");
    h1c->Draw();
    h1c1->Draw("same");
    h1c2->Draw("same");

    TLegend *leg1c = new TLegend(0.25,0.7,0.5,0.875);
    leg1c->SetBorderSize(0);leg1c->SetFillStyle(0);
    leg1c->SetHeader(Form("Pythia8: %s GeV, Q^{2} > %d GeV^{2}",beam_energies.Data(),Q2_set));
    leg1c->AddEntry("h1c","All truth-seeded tracks","l");
    leg1c->AddEntry("h1c1","+ P_{T} > 200 MeV/c","l");
    leg1c->AddEntry("h1c2","+ P_{T} > 500 MeV/c","l");
    leg1c->Draw();

    //Comparison 1
    TCanvas *c1d = new TCanvas("c1d");
    auto frame_d1 = c1d->DrawFrame(-4,0,4,1.2*h1a1->GetMaximum());
    frame_d1->GetXaxis()->SetTitle("#eta_{gen} or #eta_{rec}");frame_d1->GetXaxis()->CenterTitle();
    h1a1->Draw("same");
    h1b1->Draw("same");
    h1c1->Draw("same");

    TLegend *leg1d = new TLegend(0.125,0.675,0.575,0.875);
    leg1d->SetBorderSize(0);leg1d->SetFillStyle(0);
    leg1d->SetHeader(Form("Pythia8: %s GeV, Q^{2} > %d GeV^{2}",beam_energies.Data(),Q2_set));
    leg1d->AddEntry("h1a1","Generated charged particles w/ P_{T} > 200 MeV/c","l");
    leg1d->AddEntry("h1b1","Real-seeded tracks w/ P_{T} > 200 MeV/c","l");
    leg1d->AddEntry("h1c1","Truth-seeded tracks w/ P_{T} > 200 MeV/c","l");
    leg1d->Draw();

    //Comparison 2a
    TCanvas *c1e = new TCanvas("c1e");
    auto frame_e1 = c1e->DrawFrame(-4,0,4,1.2*h1a1->GetMaximum());
    frame_e1->GetXaxis()->SetTitle("#eta_{gen} or #eta_{rec}");frame_e1->GetXaxis()->CenterTitle();
    h1a2->Draw("same");
    h1b2->Draw("P same");

    TLegend *leg1e = new TLegend(0.125,0.675,0.575,0.875);
    leg1e->SetBorderSize(0);leg1e->SetFillStyle(0);
    leg1e->SetHeader(Form("Pythia8: %s GeV, Q^{2} > %d GeV^{2}",beam_energies.Data(),Q2_set));
    leg1e->AddEntry("h1a2","Generated charged particles w/ P_{T} > 500 MeV/c","fl");
    leg1e->AddEntry("h1b2","Real-seeded tracks w/ P_{T} > 500 MeV/c","p");
    leg1e->Draw();

    //Comparison 2b
    TCanvas *c1e1 = new TCanvas("c1e1");
    frame_e1->Draw();
    h1a2->Draw("same");
    h1c2->Draw("P same");

    TLegend *leg1e1 = new TLegend(0.125,0.675,0.575,0.875);
    leg1e1->SetBorderSize(0);leg1e1->SetFillStyle(0);
    leg1e1->SetHeader(Form("Pythia8: %s GeV, Q^{2} > %d GeV^{2}",beam_energies.Data(),Q2_set));
    leg1e1->AddEntry("h1a2","Generated charged particles w/ P_{T} > 500 MeV/c","fl");
    leg1e1->AddEntry("h1c2","Truth-seeded tracks w/ P_{T} > 500 MeV/c","p");
    leg1e1->Draw();

    //Comparison 1 ratio
    TCanvas *c1f = new TCanvas("c1f");
    h1rb1->Draw();
    h1rc1->Draw("same");

    TLegend *leg1f = new TLegend(0.575,0.25,0.875,0.45);
    leg1f->SetBorderSize(0);leg1f->SetFillStyle(0);
    leg1f->SetHeader(Form("Pythia8: %s GeV, Q^{2} > %d GeV^{2}",beam_energies.Data(),Q2_set));
    leg1f->AddEntry("h1rb1","Real-seeded tracking","l");
    leg1f->AddEntry("h1rc1","Truth-seeded tracking","l");
    leg1f->Draw();

    //Comparison 2 ratio
    TCanvas *c1g = new TCanvas("c1g");
    h1rb2->Draw();
    h1rc2->Draw("same");

    TLegend *leg1g = new TLegend(0.575,0.25,0.875,0.45);
    leg1g->SetBorderSize(0);leg1g->SetFillStyle(0);
    leg1g->SetHeader(Form("Pythia8: %s GeV, Q^{2} > %d GeV^{2}",beam_energies.Data(),Q2_set));
    leg1g->AddEntry("h1rb2","Real-seeded tracking","l");
    leg1g->AddEntry("h1rc2","Truth-seeded tracking","l");
    leg1g->Draw();

    //PT spectra -- Generated charged particles and Reconstructed truth-seeded tracks
    TCanvas *c2a = new TCanvas("c2a");
    auto frame_2a = c2a->DrawFrame(0,0.1,15,1.2*h2a->GetMaximum());
    frame_2a->GetXaxis()->SetTitle("P_{T,gen} or P_{T,rec} [GeV/c]");frame_2a->GetXaxis()->CenterTitle();
    gPad->SetLogy();
    h2a->Draw("same");
    h2b->Draw("P same");

    TLegend *leg2a = new TLegend(0.375,0.675,0.875,0.875);
    leg2a->SetBorderSize(0);leg2a->SetFillStyle(0);
    leg2a->SetHeader(Form("Pythia8: %s GeV, Q^{2} > %d GeV^{2}",beam_energies.Data(),Q2_set));
    leg2a->AddEntry("h2a","Generated charged particles w/ -3.5 < #eta < 3.5","fl");
    leg2a->AddEntry("h2b","Real-seeded tracks w/ -3.5 < #eta < 3.5","p");
    leg2a->Draw();

    //PT spectra -- Generated charged particles and Reconstructed truth-seeded tracks
    TCanvas *c2b = new TCanvas("c2b");
    frame_2a->Draw();
    gPad->SetLogy();
    h2a->Draw("same");
    h2c->Draw("P same");

    TLegend *leg2b = new TLegend(0.375,0.675,0.875,0.875);
    leg2b->SetBorderSize(0);leg2b->SetFillStyle(0);
    leg2b->SetHeader(Form("Pythia8: %s GeV, Q^{2} > %d GeV^{2}",beam_energies.Data(),Q2_set));
    leg2b->AddEntry("h2a","Generated charged particles w/ -3.5 < #eta < 3.5","fl");
    leg2b->AddEntry("h2c","Truth-seeded tracks w/ -3.5 < #eta < 3.5","p");
    leg2b->Draw();

    //Efficiency as a function of PT
    TCanvas *c3a = new TCanvas("c3a");
    auto frame_3a = c3a->DrawFrame(0,0,15,1.1);
    frame_3a->GetXaxis()->SetTitle("P_{T} [GeV/c]");frame_3a->GetXaxis()->CenterTitle();
    frame_3a->GetYaxis()->SetTitle("Efficiency");frame_3a->GetYaxis()->CenterTitle();

    h3ra->Draw("P same");
    h3rb->Draw("P same");

    TLegend *leg3a = new TLegend(0.4,0.3,0.8,0.5);
    leg3a->SetBorderSize(0);leg3a->SetFillStyle(0);
    leg3a->SetHeader(Form("Pythia8: %s GeV, Q^{2} > %d GeV^{2}",beam_energies.Data(),Q2_set));
    leg3a->AddEntry("h3ra","Real-Seeded tracks","p");
    leg3a->AddEntry("h3rb","Truth-seeded tracks","p");
    leg3a->Draw();

    //Print plots to file
    c1a->Print(Form("plots/DIS_reconstruction_E_%s_Q2_%d.pdf[",beam_energies.Data(),Q2_set));
    c1a->Print(Form("plots/DIS_reconstruction_E_%s_Q2_%d.pdf",beam_energies.Data(),Q2_set));
    c1b->Print(Form("plots/DIS_reconstruction_E_%s_Q2_%d.pdf",beam_energies.Data(),Q2_set));
    c1c->Print(Form("plots/DIS_reconstruction_E_%s_Q2_%d.pdf",beam_energies.Data(),Q2_set));
    c1d->Print(Form("plots/DIS_reconstruction_E_%s_Q2_%d.pdf",beam_energies.Data(),Q2_set));
    c1e->Print(Form("plots/DIS_reconstruction_E_%s_Q2_%d.pdf",beam_energies.Data(),Q2_set));
    c1e1->Print(Form("plots/DIS_reconstruction_E_%s_Q2_%d.pdf",beam_energies.Data(),Q2_set));
    c1f->Print(Form("plots/DIS_reconstruction_E_%s_Q2_%d.pdf",beam_energies.Data(),Q2_set));
    c1g->Print(Form("plots/DIS_reconstruction_E_%s_Q2_%d.pdf",beam_energies.Data(),Q2_set));
    c2a->Print(Form("plots/DIS_reconstruction_E_%s_Q2_%d.pdf",beam_energies.Data(),Q2_set));
    c2b->Print(Form("plots/DIS_reconstruction_E_%s_Q2_%d.pdf",beam_energies.Data(),Q2_set));
    c3a->Print(Form("plots/DIS_reconstruction_E_%s_Q2_%d.pdf",beam_energies.Data(),Q2_set));
    c3a->Print(Form("plots/DIS_reconstruction_E_%s_Q2_%d.pdf]",beam_energies.Data(),Q2_set));

}
