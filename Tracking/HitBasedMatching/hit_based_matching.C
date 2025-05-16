//------------------
// energy_set 0 is 18x275 GeV
// energy_set 1 is 10x100 GeV
void hit_based_matching(int energy_set = 0, int Q2_set = 1){

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
 
    //Track purity
    TH1 *h1a = new TH1D("h1a","Truth-seeded tracks",100,0,1.1);
    h1a->GetXaxis()->SetTitle("Fraction of track measurements from a given MC Particle");h1a->GetXaxis()->CenterTitle();
    h1a->GetYaxis()->SetTitle("Number of Tracks"); h1a->GetYaxis()->CenterTitle();
    h1a->SetLineColor(kBlack);h1a->SetLineWidth(2);

    TH1 *h1b = new TH1D("h1b","Real-seeded tracks",100,0,1.1);
    h1b->GetXaxis()->SetTitle("Fraction of track measurements from a given MC Particle");h1b->GetXaxis()->CenterTitle();
    h1b->GetYaxis()->SetTitle("Number of Tracks"); h1b->GetYaxis()->CenterTitle();
    h1b->SetLineColor(kBlack);h1b->SetLineWidth(2);

    //Momentum reconstruction
    TH2 *h2a = new TH2D("h2a","Truth-seeded tracks",100,0,50,100,0,50);
    h2a->GetXaxis()->SetTitle("Associated generated particle momentum [GeV/c]"); h2a->GetXaxis()->CenterTitle();
    h2a->GetYaxis()->SetTitle("Reconstructed track momentum [GeV/c]"); h2a->GetYaxis()->CenterTitle();

    TH2 *h2b = new TH2D("h2b","Real-seeded tracks",100,0,50,100,0,50);
    h2b->GetXaxis()->SetTitle("Associated generated particle momentum [GeV/c]"); h2b->GetXaxis()->CenterTitle();
    h2b->GetYaxis()->SetTitle("Reconstructed track momentum [GeV/c]"); h2b->GetYaxis()->CenterTitle();

    //Momentum reconstruction -- -3 < eta_gen < +3
    TH2 *h3a = new TH2D("h3a","Truth-seeded tracks",100,0,50,100,0,50);
    h3a->GetXaxis()->SetTitle("Associated generated particle momentum [GeV/c]"); h3a->GetXaxis()->CenterTitle();
    h3a->GetYaxis()->SetTitle("Reconstructed track momentum [GeV/c]"); h3a->GetYaxis()->CenterTitle();

    TH2 *h3b = new TH2D("h3b","Real-seeded tracks",100,0,50,100,0,50);
    h3b->GetXaxis()->SetTitle("Associated generated particle momentum [GeV/c]"); h3b->GetXaxis()->CenterTitle();
    h3b->GetYaxis()->SetTitle("Reconstructed track momentum [GeV/c]"); h3b->GetYaxis()->CenterTitle();

    //Status of associated particles for entirely pure tracks
    TH1 *h4a = new TH1D("h4a","Truth-seeded tracks: Entirely pure tracks",2,0,2);
    h4a->GetXaxis()->SetBinLabel(1,"Secondary Particle");
    h4a->GetXaxis()->SetBinLabel(2,"Primary Particle");
    h4a->GetXaxis()->SetTitle("Status of associated MC particle");h4a->GetXaxis()->CenterTitle();
    h4a->GetYaxis()->SetTitle("Number of Tracks"); h4a->GetYaxis()->CenterTitle();
    h4a->SetLineColor(kBlack);h4a->SetLineWidth(2);

    TH1 *h4b = new TH1D("h4b","Real-seeded tracks: Entirely pure tracks",2,0,2);
    h4b->GetXaxis()->SetBinLabel(1,"Secondary Particle");
    h4b->GetXaxis()->SetBinLabel(2,"Primary Particle");
    h4b->GetXaxis()->SetTitle("Status of associated MC particle");h4b->GetXaxis()->CenterTitle();
    h4b->GetYaxis()->SetTitle("Number of Tracks"); h4b->GetYaxis()->CenterTitle();
    h4b->SetLineColor(kBlack);h4b->SetLineWidth(2);

    //Status of associated particle for mixed tracks
    TH1 *h5a = new TH1D("h5a","Truth-seeded tracks: Status of most-associated particle for mixed tracks",2,0,2);
    h5a->GetXaxis()->SetBinLabel(1,"Secondary Particle");
    h5a->GetXaxis()->SetBinLabel(2,"Primary Particle");
    h5a->GetXaxis()->SetTitle("Status of associated MC particle");h5a->GetXaxis()->CenterTitle();
    h5a->GetYaxis()->SetTitle("Number of Tracks"); h5a->GetYaxis()->CenterTitle();
    h5a->SetLineColor(kBlack);h5a->SetLineWidth(2);

    TH1 *h6a = new TH1D("h6a","Truth-seeded tracks: Status of less-associated particles for mixed tracks",2,0,2);
    h6a->GetXaxis()->SetBinLabel(1,"Secondary Particle");
    h6a->GetXaxis()->SetBinLabel(2,"Primary Particle");
    h6a->GetXaxis()->SetTitle("Status of associated MC particle");h6a->GetXaxis()->CenterTitle();
    h6a->GetYaxis()->SetTitle("Number of Tracks"); h6a->GetYaxis()->CenterTitle();
    h6a->SetLineColor(kBlack);h6a->SetLineWidth(2);

    TH1 *h5b = new TH1D("h5b","Real-seeded tracks: Status of most-associated particle for mixed tracks",2,0,2);
    h5b->GetXaxis()->SetBinLabel(1,"Secondary Particle");
    h5b->GetXaxis()->SetBinLabel(2,"Primary Particle");
    h5b->GetXaxis()->SetTitle("Status of associated MC particle");h5b->GetXaxis()->CenterTitle();
    h5b->GetYaxis()->SetTitle("Number of Tracks"); h5b->GetYaxis()->CenterTitle();
    h5b->SetLineColor(kBlack);h5b->SetLineWidth(2);

    TH1 *h6b = new TH1D("h6b","Real-seeded tracks: Status of less-associated particles for mixed tracks",2,0,2);
    h6b->GetXaxis()->SetBinLabel(1,"Secondary Particle");
    h6b->GetXaxis()->SetBinLabel(2,"Primary Particle");
    h6b->GetXaxis()->SetTitle("Status of associated MC particle");h6b->GetXaxis()->CenterTitle();
    h6b->GetYaxis()->SetTitle("Number of Tracks"); h6b->GetYaxis()->CenterTitle();
    h6b->SetLineColor(kBlack);h6b->SetLineWidth(2);

    //File information
    TString path = "./input/";
    TString run_name;
    TString beam_energies;

    if(energy_set == 0 && Q2_set == 1){
        beam_energies = "18x275";
        run_name = "Q2_1/eicrecon_*.root"; //Pythia8 DIS events with 18x275 GeV and Q2>1 GeV2
    }
    else if(energy_set == 0 && Q2_set == 10){
        beam_energies = "18x275";
        run_name = "eicrecon_out_E_18_275_Q2_10.root"; //Pythia8 DIS events with 18x275 GeV and Q2>10 GeV2
    }
    else if(energy_set == 0 && Q2_set == 100){
        beam_energies = "18x275";
        run_name = "Q2_100/eicrecon_*.root"; //Pythia8 DIS events with 18x275 GeV and Q2>100 GeV2
    }
    else if(energy_set == 1 && Q2_set == 1){
        beam_energies = "10x100";
        run_name = "eicrecon_out_E_10_100_Q2_1.root"; //Pythia8 DIS events with 10x100 GeV and Q2>1 GeV2
    }
    else if(energy_set == 1 && Q2_set == 10){
        beam_energies = "10x100";
        run_name = "eicrecon_out_E_10_100_Q2_10.root"; //Pythia8 DIS events with 10x100 GeV and Q2>10 GeV2
    }
    else if(energy_set == 1 && Q2_set == 100){
        beam_energies = "10x100";
        run_name = "eicrecon_out_E_10_100_Q2_100.root"; //Pythia8 DIS events with 10x100 GeV and Q2>100 GeV2
    }

    //Open File
    TString input = path + run_name;
    TChain *tree = new TChain("events");
    tree->Add(input.Data());

    cout<<"Analyzing "<<tree->GetEntries()<<" events!"<<endl;

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

    TTreeReaderArray<float> hit_assoc_weight_ts(tr,"CentralCKFTruthSeededTrackAssociations.weight");
    TTreeReaderArray<int>   hit_assoc_mc_ts(tr, "_CentralCKFTruthSeededTrackAssociations_sim.index");
    TTreeReaderArray<int>   hit_assoc_track_ts(tr, "_CentralCKFTruthSeededTrackAssociations_rec.index");

    TTreeReaderArray<float> hit_assoc_weight_rs(tr,"CentralCKFTrackAssociations.weight");
    TTreeReaderArray<int>   hit_assoc_mc_rs(tr, "_CentralCKFTrackAssociations_sim.index");
    TTreeReaderArray<int>   hit_assoc_track_rs(tr, "_CentralCKFTrackAssociations_rec.index");

    TTreeReaderArray<float> track_qoverp_ts(tr, "CentralCKFTruthSeededTrackParameters.qOverP");
    TTreeReaderArray<float> track_qoverp_rs(tr, "CentralCKFTrackParameters.qOverP");
    
    //Other variables
    TLorentzVector gen_vec;
    TVector3 gen_vertex;

    //Counters for truth-seeded associations
    int pure_track_counter_ts(0);
    int mixed_track_counter_ts(0);

    //Counter for real-seeded associations
    int pure_track_counter_rs(0);
    int mixed_track_counter_rs(0);

    TLorentzVector rec_vec;
    int counter(0);

    //Loop over events
    while (tr.Next()) {

	    if(counter%100==0) cout<<"Analyzing event "<<counter<<endl;
	    counter++;

        //Define vector of tuples containing (track index, mc particle status, association weight)
        std::vector<std::tuple<int, int, float>> vec;

        //Loop over (truth-seeded) hit-based associations
        for(int iassoc=0;iassoc<hit_assoc_weight_ts.GetSize();iassoc++){

            auto assoc_weight = hit_assoc_weight_ts[iassoc];

            h1a->Fill(assoc_weight);

            //Count pure and mixed tracks (avoid double-counting of mixed tracks)
            if(assoc_weight > 0.99)
                pure_track_counter_ts++;
            else if(iassoc==0)
                mixed_track_counter_ts++;
            else if(iassoc>0 && hit_assoc_track_ts[iassoc] != hit_assoc_track_ts[iassoc-1]) //I assume the associated track index is sorted
                mixed_track_counter_ts++;

            //Get the associated track and MC particle
            //N.B. We use the fact that (right now) the track collection index is the same as the track parameter collection index
            int mc_index = hit_assoc_mc_ts[iassoc];
            int trk_index = hit_assoc_track_ts[iassoc];

            gen_vec.SetXYZM(gen_px[mc_index],gen_py[mc_index],gen_pz[mc_index],gen_mass[mc_index]);
            auto gen_mom = gen_vec.P();

            auto trk_mom = fabs(1./track_qoverp_ts[trk_index]);

            h2a->Fill(gen_mom,trk_mom);
            if( fabs(gen_vec.Eta()) < 3 )
                h3a->Fill(gen_mom,trk_mom);

            //Check status (i.e. primary or secondary) of associated particle
            int status = gen_status[mc_index];

            if(assoc_weight > 0.99){
                h4a->Fill(status);
            }
            else{ //Fill vector of tuples for mixed tracks
                vec.push_back(std::make_tuple(trk_index,status,assoc_weight));
            }

        }

        //Sort vector by first tuple index (not needed since should already be in ascending order)
        std::sort(vec.begin(), vec.end());

        //Group based on the first tuple index (i.e. track index)
        size_t ivec = 0;
        while (ivec < vec.size()) {

            int current_first = std::get<0>(vec[ivec]);
            float current_max_weight = std::get<2>(vec[ivec]);
            size_t start = ivec;
            size_t index_max_weight = ivec;

            //Find all consecutive elements with the same first value
            while (ivec < vec.size() && std::get<0>(vec[ivec]) == current_first) {
                float current_weight = std::get<2>(vec[ivec]);
                if( current_weight > current_max_weight ){ 
                    current_max_weight = current_weight; 
                    index_max_weight = ivec;
                }
                ++ivec;
            }

            //Fill histograms with MC particle status values
            if (ivec - start > 1) {
                for (size_t ii = start; ii < ivec; ++ii){
                    
                    if(ii == index_max_weight)
                        h5a->Fill( std::get<1>(vec[ii]) );
                    else
                        h6a->Fill( std::get<1>(vec[ii]) );
                }
            }
        } //End while loop over vec

        //Clear vec and reset ivec for real-seeded tracking use
        vec.clear();
        ivec = 0;
        //End of truth-seeded tracking portion

        //Loop over (real-seeded) hit-based associations
        for(int iassoc=0;iassoc<hit_assoc_weight_rs.GetSize();iassoc++){

            auto assoc_weight = hit_assoc_weight_rs[iassoc];

            h1b->Fill(assoc_weight);

            //Count pure and mixed tracks (avoid double-counting of mixed tracks)
            if(assoc_weight > 0.99)
                pure_track_counter_rs++;
            else if(iassoc==0)
                mixed_track_counter_rs++;
            else if(iassoc>0 && hit_assoc_track_rs[iassoc] != hit_assoc_track_rs[iassoc-1]) //I assume the associated track index is sorted
                mixed_track_counter_rs++;
 
            //Get the associated track and MC particle
            //N.B. We use the fact that (right now) the track collection index is the same as the track parameter collection index
            int mc_index = hit_assoc_mc_rs[iassoc];
            int trk_index = hit_assoc_track_rs[iassoc];

            gen_vec.SetXYZM(gen_px[mc_index],gen_py[mc_index],gen_pz[mc_index],gen_mass[mc_index]);
            auto gen_mom = gen_vec.P();

            auto trk_mom = fabs(1./track_qoverp_rs[trk_index]);

            h2b->Fill(gen_mom,trk_mom);

            if( fabs(gen_vec.Eta()) < 3 )
                h3b->Fill(gen_mom,trk_mom);

            //Check status (i.e. primary or secondary) of associated particle
            int status = gen_status[mc_index];

            if(assoc_weight > 0.99){
                h4b->Fill(status);
            }
            else{ //Fill vector of tuples for mixed tracks
                vec.push_back(std::make_tuple(trk_index,status,assoc_weight));
            }
 
        }

        //Sort vector by first tuple index (not needed since should already be in ascending order)
        std::sort(vec.begin(), vec.end());

        //Group based on the first tuple index (i.e. track index)
        while (ivec < vec.size()) {

            int current_first = std::get<0>(vec[ivec]);
            float current_max_weight = std::get<2>(vec[ivec]);
            size_t start = ivec;
            size_t index_max_weight = ivec;

            //Find all consecutive elements with the same first value
            while (ivec < vec.size() && std::get<0>(vec[ivec]) == current_first) {
                float current_weight = std::get<2>(vec[ivec]);
                if( current_weight > current_max_weight ){ 
                    current_max_weight = current_weight; 
                    index_max_weight = ivec;
                }
                ++ivec;
            }

            //Fill histograms with MC particle status values
            if (ivec - start > 1) {
                for (size_t ii = start; ii < ivec; ++ii){
                    
                    if(ii == index_max_weight)
                        h5b->Fill( std::get<1>(vec[ii]) );
                    else
                        h6b->Fill( std::get<1>(vec[ii]) );
                }
            }
        } //End while loop over vec
        //End of real-seeded tracking portion

    } //End loop over events

    //Calculate fraction of pure tracks
    float frac_pure_ts = ((float) pure_track_counter_ts) / ((float) (pure_track_counter_ts + mixed_track_counter_ts));
    float frac_pure_rs = ((float) pure_track_counter_rs) / ((float) (pure_track_counter_rs + mixed_track_counter_rs));

    //Make plots
    //----
    //Association weights -- truth-seeding
    TCanvas *c1a = new TCanvas("c1a");
    c1a->SetLogy();
    h1a->Draw();

    TLegend *leg1a = new TLegend(0.125,0.7,0.475,0.875);
    leg1a->SetBorderSize(0);leg1a->SetFillStyle(0);
    leg1a->SetHeader(Form("Pythia8: %s GeV, Q^{2} > %d GeV^{2}",beam_energies.Data(),Q2_set));
    leg1a->AddEntry(h1a, Form("Fraction of entirely pure tracks = %.3f",frac_pure_ts),"l");
    leg1a->Draw();

    //Association weights -- real-seeding
    TCanvas *c1b = new TCanvas("c1b");
    c1b->SetLogy();
    h1b->Draw();

    TLegend *leg1b = new TLegend(0.125,0.7,0.475,0.875);
    leg1b->SetBorderSize(0);leg1b->SetFillStyle(0);
    leg1b->SetHeader(Form("Pythia8: %s GeV, Q^{2} > %d GeV^{2}",beam_energies.Data(),Q2_set));
    leg1b->AddEntry(h1b, Form("Fraction of entirely pure tracks = %.3f",frac_pure_rs),"l");
    leg1b->Draw();

    //Track vs. Associated MC particle Momentum -- truth-seeding
    TCanvas *c2a = new TCanvas("c2a");
    c2a->SetLogz();
    h2a->Draw("colz");

    TLatex *tex2a = new TLatex(25,8,Form("Pythia8: %s GeV, Q^{2} > %d GeV^{2}",beam_energies.Data(),Q2_set));
    tex2a->SetTextSize(0.03);
    tex2a->Draw();

    //Track vs. Associated MC particle Momentum -- real-seeding
    TCanvas *c2b = new TCanvas("c2b");
    c2b->SetLogz();
    h2b->Draw("colz");

    TLatex *tex2b = new TLatex(25,8,Form("Pythia8: %s GeV, Q^{2} > %d GeV^{2}",beam_energies.Data(),Q2_set));
    tex2b->SetTextSize(0.03);
    tex2b->Draw();

    //Track vs. Associated MC particle Momentum -- truth-seeding, eta_gen restriction
    TCanvas *c3a = new TCanvas("c3a");
    c3a->SetLogz();
    h3a->Draw("colz");

    TLatex *tex3a = new TLatex();
    tex3a->SetTextSize(0.03);
    tex3a->DrawLatex(25,8,Form("Pythia8: %s GeV, Q^{2} > %d GeV^{2}",beam_energies.Data(),Q2_set));
    tex3a->DrawLatex(25,5," -3 < #eta_{gen} < +3");

    //Track vs. Associated MC particle Momentum -- real-seeding, eta_gen restriction
    TCanvas *c3b = new TCanvas("c3b");
    c3b->SetLogz();
    h3b->Draw("colz");

    TLatex *tex3b = new TLatex();
    tex3b->SetTextSize(0.03);
    tex3b->DrawLatex(25,8,Form("Pythia8: %s GeV, Q^{2} > %d GeV^{2}",beam_energies.Data(),Q2_set));
    tex3b->DrawLatex(25,5," -3 < #eta_{gen} < +3");

    //Status of associated particles for entirely pure tracks -- truth seeding
    TCanvas *c4a = new TCanvas("c4a");
    c4a->SetLogy();
    h4a->Draw();

    TLegend *leg4a = new TLegend(0.125,0.7,0.475,0.875);
    leg4a->SetBorderSize(0);leg4a->SetFillStyle(0);
    leg4a->SetHeader(Form("Pythia8: %s GeV, Q^{2} > %d GeV^{2}",beam_energies.Data(),Q2_set));
    leg4a->Draw();

    //Status of associated particles for entirely pure tracks -- real seeding
    TCanvas *c4b = new TCanvas("c4b");
    c4b->SetLogy();
    h4b->Draw();

    TLegend *leg4b = new TLegend(0.125,0.7,0.475,0.875);
    leg4b->SetBorderSize(0);leg4b->SetFillStyle(0);
    leg4b->SetHeader(Form("Pythia8: %s GeV, Q^{2} > %d GeV^{2}",beam_energies.Data(),Q2_set));
    leg4b->Draw();

    //Status of associated particles for mixed tracks -- truth seeding
    TCanvas *c5a = new TCanvas("c5a");
    c5a->SetLogy();
    h5a->Draw();

    TLegend *leg5a = new TLegend(0.125,0.7,0.475,0.875);
    leg5a->SetBorderSize(0);leg5a->SetFillStyle(0);
    leg5a->SetHeader(Form("Pythia8: %s GeV, Q^{2} > %d GeV^{2}",beam_energies.Data(),Q2_set));
    leg5a->Draw();

    TCanvas *c6a = new TCanvas("c6a");
    c6a->SetLogy();
    h6a->Draw();

    TLegend *leg6a = new TLegend(0.525,0.7,0.875,0.875);
    leg6a->SetBorderSize(0);leg6a->SetFillStyle(0);
    leg6a->SetHeader(Form("Pythia8: %s GeV, Q^{2} > %d GeV^{2}",beam_energies.Data(),Q2_set));
    leg6a->Draw();

    //Status of associated particles for mixed tracks -- real seeding
    TCanvas *c5b = new TCanvas("c5b");
    c5b->SetLogy();
    h5b->Draw();

    TLegend *leg5b = new TLegend(0.125,0.7,0.475,0.875);
    leg5b->SetBorderSize(0);leg5b->SetFillStyle(0);
    leg5b->SetHeader(Form("Pythia8: %s GeV, Q^{2} > %d GeV^{2}",beam_energies.Data(),Q2_set));
    leg5b->Draw();

    TCanvas *c6b = new TCanvas("c6b");
    c6b->SetLogy();
    h6b->Draw();

    TLegend *leg6b = new TLegend(0.525,0.7,0.875,0.875);
    leg6b->SetBorderSize(0);leg6b->SetFillStyle(0);
    leg6b->SetHeader(Form("Pythia8: %s GeV, Q^{2} > %d GeV^{2}",beam_energies.Data(),Q2_set));
    leg6b->Draw();

    //Print plots to file
    c1a->Print(Form("plots/hit_based_matching_%s_Q2_%d.pdf[",beam_energies.Data(),Q2_set));
    c1a->Print(Form("plots/hit_based_matching_%s_Q2_%d.pdf",beam_energies.Data(),Q2_set));
    c1b->Print(Form("plots/hit_based_matching_%s_Q2_%d.pdf",beam_energies.Data(),Q2_set));
    c2a->Print(Form("plots/hit_based_matching_%s_Q2_%d.pdf",beam_energies.Data(),Q2_set));
    c2b->Print(Form("plots/hit_based_matching_%s_Q2_%d.pdf",beam_energies.Data(),Q2_set));
    c3a->Print(Form("plots/hit_based_matching_%s_Q2_%d.pdf",beam_energies.Data(),Q2_set));
    c3b->Print(Form("plots/hit_based_matching_%s_Q2_%d.pdf",beam_energies.Data(),Q2_set));
    c4a->Print(Form("plots/hit_based_matching_%s_Q2_%d.pdf",beam_energies.Data(),Q2_set));
    c4b->Print(Form("plots/hit_based_matching_%s_Q2_%d.pdf",beam_energies.Data(),Q2_set));
    c5a->Print(Form("plots/hit_based_matching_%s_Q2_%d.pdf",beam_energies.Data(),Q2_set));
    c6a->Print(Form("plots/hit_based_matching_%s_Q2_%d.pdf",beam_energies.Data(),Q2_set));
    c5b->Print(Form("plots/hit_based_matching_%s_Q2_%d.pdf",beam_energies.Data(),Q2_set));
    c6b->Print(Form("plots/hit_based_matching_%s_Q2_%d.pdf",beam_energies.Data(),Q2_set));
    c6b->Print(Form("plots/hit_based_matching_%s_Q2_%d.pdf]",beam_energies.Data(),Q2_set));

}
