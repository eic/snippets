//------------------
void track_eff_fixed_mom(std::string gen_part = "muminus"){

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

    //Define simulated momentum points
    const int N_mom = 6;
    double fixed_mom[N_mom] = {0.5,1.0,2.0,5.0,10.0,20.0}; //In GeV/c

    //Tracker Efficiency Histograms
    TH1 *h_eta_gen[N_mom]; //Generated particle eta distribution
    TH1 *h_eta_gen_ts_eff[N_mom]; //Generated particle eta distribution w/ at least 1 truth-seeded track found
    TH1 *h_eta_ts_ratio[N_mom]; //Efficiency ratio for truth-seeded tracking
    TH1 *h_eta_gen_rs_eff[N_mom]; //Generated particle eta distribution w/ at least 1 real-seeded track found
    TH1 *h_eta_rs_ratio[N_mom]; //Efficiency ratio for real-seeded tracking

    //Only keep events with -3.5<eta_true<3.5
    TH1 *h_phi_gen[N_mom]; //Generated particle phi distribution
    TH1 *h_phi_gen_ts_eff[N_mom]; //Generated particle phi distribution w/ at least 1 truth-seeded track found
    TH1 *h_phi_ts_ratio[N_mom]; //Efficiency ratio for truth-seeded tracking
    TH1 *h_phi_gen_rs_eff[N_mom]; //Generated particle phi distribution w/ at least 1 real-seeded track found
    TH1 *h_phi_rs_ratio[N_mom]; //Efficiency ratio for real-seeded tracking

    //Loop over momentum points
    for(int imom = 0; imom < N_mom; imom++){

        cout<<endl<<"Analyzing momentum = "<< fixed_mom[imom] <<"GeV/c point!"<< endl;

        //Define histograms
        //---Eta---
        h_eta_gen[imom] = new TH1D(Form("h_eta_gen[%d]",imom),Form("Tracker Efficiency vs. generated particle #eta: p_{gen} = %.2f GeV/c",fixed_mom[imom]),100,-4,4);
        h_eta_gen[imom]->GetXaxis()->SetTitle("#eta");h_eta_gen[imom]->GetXaxis()->CenterTitle();
        h_eta_gen[imom]->GetYaxis()->SetTitle("Efficiency");h_eta_gen[imom]->GetYaxis()->CenterTitle();
        h_eta_gen[imom]->SetLineColor(kBlue);h_eta_gen[imom]->SetLineWidth(2);

        h_eta_gen_ts_eff[imom] = new TH1D(Form("h_eta_gen_ts_eff[%d]",imom),Form("Tracker Efficiency vs. generated particle #eta: p_{gen} = %.2f GeV/c",fixed_mom[imom]),100,-4,4);
        h_eta_gen_ts_eff[imom]->GetXaxis()->SetTitle("#eta");h_eta_gen_ts_eff[imom]->GetXaxis()->CenterTitle();
        h_eta_gen_ts_eff[imom]->GetYaxis()->SetTitle("Efficiency");h_eta_gen_ts_eff[imom]->GetYaxis()->CenterTitle();
        h_eta_gen_ts_eff[imom]->SetLineColor(kBlue);h_eta_gen_ts_eff[imom]->SetLineWidth(2);

        h_eta_gen_rs_eff[imom] = new TH1D(Form("h_eta_gen_rs_eff[%d]",imom),Form("Tracker Efficiency vs. generated particle #eta: p_{gen} = %.2f GeV/c",fixed_mom[imom]),100,-4,4);
        h_eta_gen_rs_eff[imom]->GetXaxis()->SetTitle("#eta");h_eta_gen_rs_eff[imom]->GetXaxis()->CenterTitle();
        h_eta_gen_rs_eff[imom]->GetYaxis()->SetTitle("Efficiency");h_eta_gen_rs_eff[imom]->GetYaxis()->CenterTitle();
        h_eta_gen_rs_eff[imom]->SetLineColor(kGreen);h_eta_gen_rs_eff[imom]->SetLineWidth(2);

	//---Phi---
	 h_phi_gen[imom] = new TH1D(Form("h_phi_gen[%d]",imom),Form("Tracker Efficiency vs. generated particle #phi: p_{gen} = %.2f GeV/c",fixed_mom[imom]),100,-3.40,3.40);
        h_phi_gen[imom]->GetXaxis()->SetTitle("#phi");h_phi_gen[imom]->GetXaxis()->CenterTitle();
        h_phi_gen[imom]->GetYaxis()->SetTitle("Efficiency");h_phi_gen[imom]->GetYaxis()->CenterTitle();
        h_phi_gen[imom]->SetLineColor(kBlue);h_phi_gen[imom]->SetLineWidth(2);

        h_phi_gen_ts_eff[imom] = new TH1D(Form("h_phi_gen_ts_eff[%d]",imom),Form("Tracker Efficiency vs. generated particle #phi: p_{gen} = %.2f GeV/c",fixed_mom[imom]),100,-3.40,3.40);
        h_phi_gen_ts_eff[imom]->GetXaxis()->SetTitle("#phi");h_phi_gen_ts_eff[imom]->GetXaxis()->CenterTitle();
        h_phi_gen_ts_eff[imom]->GetYaxis()->SetTitle("Efficiency");h_phi_gen_ts_eff[imom]->GetYaxis()->CenterTitle();
        h_phi_gen_ts_eff[imom]->SetLineColor(kBlue);h_phi_gen_ts_eff[imom]->SetLineWidth(2);

        h_phi_gen_rs_eff[imom] = new TH1D(Form("h_phi_gen_rs_eff[%d]",imom),Form("Tracker Efficiency vs. generated particle #phi: p_{gen} = %.2f GeV/c",fixed_mom[imom]),100,-3.40,3.40);
        h_phi_gen_rs_eff[imom]->GetXaxis()->SetTitle("#phi");h_phi_gen_rs_eff[imom]->GetXaxis()->CenterTitle();
        h_phi_gen_rs_eff[imom]->GetYaxis()->SetTitle("Efficiency");h_phi_gen_rs_eff[imom]->GetYaxis()->CenterTitle();
        h_phi_gen_rs_eff[imom]->SetLineColor(kGreen);h_phi_gen_rs_eff[imom]->SetLineWidth(2);

        //File information
        float x_gen;
        float y_gen;
        float z_gen;
        TString run_name;
        TString path = "./output/";
        int pid_code;

        if(gen_part=="muminus"){

            pid_code = 13;

            if(imom==0){
                run_name = "eicrecon_out_500MeV.root";
                x_gen = 0;
                y_gen = 0;
                z_gen = 0;
            }
            else if(imom==1){
                run_name = "eicrecon_out_1GeV.root";
                x_gen = 0;
                y_gen = 0;
                z_gen = 0;
            }
            else if(imom==2){
                run_name = "eicrecon_out_2GeV.root";
                x_gen = 0;
                y_gen = 0;
                z_gen = 0;
            }
            else if(imom==3){
                run_name = "eicrecon_out_5GeV.root"; 
                x_gen = 0;
                y_gen = 0;
                z_gen = 0;
            }
	    else if(imom==4){
                run_name = "eicrecon_out_10GeV.root";
                x_gen = 0;
                y_gen = 0;
                z_gen = 0;
            }
  	    else if(imom==5){
                run_name = "eicrecon_out_20GeV.root";
                x_gen = 0;
                y_gen = 0;
                z_gen = 0;
            }
        }

        //Open File
        TString input = path + run_name;
        TFile *f = new TFile(input.Data());
        TTree *tree = (TTree*) f->Get("events");

        //Create Array Reader
        TTreeReader tr(tree);

        TTreeReaderArray<int>   gen_status(tr, "MCParticles.generatorStatus");
        TTreeReaderArray<int>   gen_pid(tr, "MCParticles.PDG");
        TTreeReaderArray<float> gen_px(tr, "MCParticles.momentum.x");
        TTreeReaderArray<float> gen_py(tr, "MCParticles.momentum.y");
        TTreeReaderArray<float> gen_pz(tr, "MCParticles.momentum.z");
        TTreeReaderArray<double> gen_mass(tr, "MCParticles.mass"); //Not important here
        TTreeReaderArray<float> gen_charge(tr, "MCParticles.charge");
        TTreeReaderArray<double> gen_vx(tr, "MCParticles.vertex.x");
        TTreeReaderArray<double> gen_vy(tr, "MCParticles.vertex.y");
        TTreeReaderArray<double> gen_vz(tr, "MCParticles.vertex.z");

        TTreeReaderArray<float> track_ts_qoverp(tr, "CentralCKFTruthSeededTrackParameters.qOverP");
        TTreeReaderArray<float> track_ts_theta(tr, "CentralCKFTruthSeededTrackParameters.theta");
        TTreeReaderArray<float> track_ts_phi(tr, "CentralCKFTruthSeededTrackParameters.phi");
        TTreeReaderArray<float> track_ts_loca(tr, "CentralCKFTruthSeededTrackParameters.loc.a");
        TTreeReaderArray<float> track_ts_locb(tr, "CentralCKFTruthSeededTrackParameters.loc.b");

        TTreeReaderArray<float> track_rs_qoverp(tr, "CentralCKFTrackParameters.qOverP");
        TTreeReaderArray<float> track_rs_theta(tr, "CentralCKFTrackParameters.theta");
        TTreeReaderArray<float> track_rs_phi(tr, "CentralCKFTrackParameters.phi");
        TTreeReaderArray<float> track_rs_loca(tr, "CentralCKFTrackParameters.loc.a");
        TTreeReaderArray<float> track_rs_locb(tr, "CentralCKFTrackParameters.loc.b");

        //Other variables
        TLorentzVector gen_vec;
        TVector3 gen_vertex;
        float charge;
        bool found_primary(false);
        int counter(0);

        //Loop over events
        while (tr.Next()) {

	    if(counter%100==0) cout<<"Analyzing event "<<counter<<endl;
	    counter++;

            //Reset variables
            found_primary = false;

            //Loop over generated particles, select primary particle (assuming single particle)
            for(int igen=0;igen<gen_status.GetSize();igen++){
	            //PID requirement so that background particles are not used, but may not be needed because of 'break'
                if(gen_status[igen]==1 && gen_pid[igen]==pid_code){
                    gen_vec.SetXYZM(gen_px[igen],gen_py[igen],gen_pz[igen],gen_mass[igen]);
                    gen_vertex.SetXYZ(gen_vx[igen],gen_vy[igen],gen_vz[igen]);
                    charge = gen_charge[igen];
                    found_primary = true;
                    break;
                }
            }

            //Require primary particle for all track results
            if(found_primary){
 
                h_eta_gen[imom]->Fill(gen_vec.Eta());
		if( std::abs(gen_vec.Eta())<3.0 ) h_phi_gen[imom]->Fill(gen_vec.Phi());

                //Truth-seeded tracking
                int track_ts_mult = track_ts_qoverp.GetSize();
                if(track_ts_mult>0) {
                    h_eta_gen_ts_eff[imom]->Fill(gen_vec.Eta());
		    if( std::abs(gen_vec.Eta())<3.0 ) h_phi_gen_ts_eff[imom]->Fill(gen_vec.Phi());
                }

                //Real seeded tracking
                int track_rs_mult = track_rs_qoverp.GetSize();
                if(track_rs_mult>0) {
                    h_eta_gen_rs_eff[imom]->Fill(gen_vec.Eta());
		    if( std::abs(gen_vec.Eta())<3.0 ) h_phi_gen_rs_eff[imom]->Fill(gen_vec.Phi());
                }
            }
        } // End loop over events

        //Divide histograms
        h_eta_gen_ts_eff[imom] = (TH1*) h_eta_gen_ts_eff[imom]->Clone(Form("h_eta_gen_ts_eff[%d]",imom)); 
        h_eta_gen_ts_eff[imom]->Divide(h_eta_gen[imom]);

        h_eta_gen_rs_eff[imom] = (TH1*) h_eta_gen_rs_eff[imom]->Clone(Form("h_eta_gen_rs_eff[%d]",imom)); 
        h_eta_gen_rs_eff[imom]->Divide(h_eta_gen[imom]);

	h_phi_gen_ts_eff[imom] = (TH1*) h_phi_gen_ts_eff[imom]->Clone(Form("h_phi_gen_ts_eff[%d]",imom));
        h_phi_gen_ts_eff[imom]->Divide(h_phi_gen[imom]);

        h_phi_gen_rs_eff[imom] = (TH1*) h_phi_gen_rs_eff[imom]->Clone(Form("h_phi_gen_rs_eff[%d]",imom));
        h_phi_gen_rs_eff[imom]->Divide(h_phi_gen[imom]);

    } // End loop over momentum points


    //Make plots
    TCanvas *c1_ts = new TCanvas("c1_ts");
    c1_ts->Divide(2,2);
    c1_ts->cd(1); h_eta_gen_ts_eff[0]->Draw();
    c1_ts->cd(2); h_eta_gen_ts_eff[1]->Draw();
    c1_ts->cd(3); h_eta_gen_ts_eff[2]->Draw();
    c1_ts->cd(4); h_eta_gen_ts_eff[3]->Draw();

    TCanvas *c1_rs = new TCanvas("c1_rs");
    c1_rs->Divide(2,2);
    c1_rs->cd(1); h_eta_gen_rs_eff[0]->Draw();
    c1_rs->cd(2); h_eta_gen_rs_eff[1]->Draw();
    c1_rs->cd(3); h_eta_gen_rs_eff[2]->Draw();
    c1_rs->cd(4); h_eta_gen_rs_eff[3]->Draw();
 
    TCanvas *c2_ts = new TCanvas("c2_ts");
    c2_ts->Divide(2,2);
    c2_ts->cd(1); h_eta_gen_ts_eff[4]->Draw();
    c2_ts->cd(2); h_eta_gen_ts_eff[5]->Draw();

    TCanvas *c2_rs = new TCanvas("c2_rs");
    c2_rs->Divide(2,2);
    c2_rs->cd(1); h_eta_gen_rs_eff[4]->Draw();
    c2_rs->cd(2); h_eta_gen_rs_eff[5]->Draw();

    TCanvas *c3_ts = new TCanvas("c3_ts");
    c3_ts->Divide(2,2);
    c3_ts->cd(1); h_phi_gen_ts_eff[0]->Draw();
    c3_ts->cd(2); h_phi_gen_ts_eff[1]->Draw();
    c3_ts->cd(3); h_phi_gen_ts_eff[2]->Draw();
    c3_ts->cd(4); h_phi_gen_ts_eff[3]->Draw();

    TCanvas *c3_rs = new TCanvas("c3_rs");
    c3_rs->Divide(2,2);
    c3_rs->cd(1); h_phi_gen_rs_eff[0]->Draw();
    c3_rs->cd(2); h_phi_gen_rs_eff[1]->Draw();
    c3_rs->cd(3); h_phi_gen_rs_eff[2]->Draw();
    c3_rs->cd(4); h_phi_gen_rs_eff[3]->Draw();

    TCanvas *c4_ts = new TCanvas("c4_ts");
    c4_ts->Divide(2,2);
    c4_ts->cd(1); h_phi_gen_ts_eff[4]->Draw();
    c4_ts->cd(2); h_phi_gen_ts_eff[5]->Draw();

    TCanvas *c4_rs = new TCanvas("c4_rs");
    c4_rs->Divide(2,2);
    c4_rs->cd(1); h_phi_gen_rs_eff[4]->Draw();
    c4_rs->cd(2); h_phi_gen_rs_eff[5]->Draw();

    //Print plots to file
    c1_ts->Print("plots/track_eff_fixed_mom.pdf[");
    c1_ts->Print("plots/track_eff_fixed_mom.pdf");
    c1_rs->Print("plots/track_eff_fixed_mom.pdf");
    c2_ts->Print("plots/track_eff_fixed_mom.pdf");
    c2_rs->Print("plots/track_eff_fixed_mom.pdf");
    c3_ts->Print("plots/track_eff_fixed_mom.pdf");
    c3_rs->Print("plots/track_eff_fixed_mom.pdf");
    c4_ts->Print("plots/track_eff_fixed_mom.pdf");
    c4_rs->Print("plots/track_eff_fixed_mom.pdf");
    c4_rs->Print("plots/track_eff_fixed_mom.pdf]");

}
