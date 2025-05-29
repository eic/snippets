#include <edm4hep/SimTrackerHitCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <edm4eic/Measurement2DCollection.h>
#include <edm4eic/TrackCollection.h>
#include <edm4eic/MCRecoTrackerHitAssociationCollection.h>
#include <podio/Frame.h>
#include <podio/ROOTReader.h>

//SimHit collections
std::vector<std::string> sim_coll_names{
				 	"VertexBarrelHits","SiBarrelHits","TrackerEndcapHits",
					"MPGDBarrelHits","BackwardMPGDEndcapHits","ForwardMPGDEndcapHits","OuterMPGDBarrelHits",
					"TOFBarrelHits","TOFEndcapHits"
					};

//TrackerHit RecHit (digitized) collections
std::vector<std::string> rec_coll_names{
					"SiBarrelVertexRecHits","SiBarrelTrackerRecHits","SiEndcapTrackerRecHits",
					"MPGDBarrelRecHits","BackwardMPGDEndcapRecHits","ForwardMPGDEndcapRecHits","OuterMPGDBarrelRecHits",
					"TOFBarrelRecHits","TOFEndcapRecHits"
					};

//Collection IDs for RecHit collections.
std::vector<unsigned int> rec_coll_ids;

//Simulation file to use
//std::string input_file = "input/eicrecon_etof_test_local.root";
//std::string input_file = "input/eicrecon_etof_test.root";
std::string input_file = "input/eicrecon_out_2GeV.root";


//------------------
//Template to access 'sign' of radius in (x,y) plane
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

//-----------------
//Function to set Collection IDs
vector<unsigned int> get_coll_ids(vector<string> coll_names){

	vector<unsigned int> coll_ids;
	
        TFile *f = new TFile(input_file.c_str());
        TTree *t = (TTree*) f->Get("podio_metadata");

        TTreeReader tr(t);
        TTreeReaderArray<unsigned int> c_ids(tr,"events___idTable.m_collectionIDs");
        TTreeReaderArray<string> c_names(tr,"events___idTable.m_names");	

	tr.Next();

	//Find the associated collection ID
	for(int iname = 0; iname < coll_names.size(); iname++) {
		for(int icol=0;icol<c_ids.GetSize();icol++){
			if( c_names[icol] == coll_names[iname] )
				coll_ids.push_back(c_ids[icol]);
		}
	}

	delete f;

	return coll_ids;
}

//-----------------
// Main function
void hit_matching(){

	//Track collection to use: real-seeded or truth-seeded
	std::string trk_coll = "CentralCKFTracks";

	//Print out information event-by-event
	bool print_evt_info = false;

	//Using new Barrel MPDG digitization
	bool new_mpgd_digi = false;

	//Plot x-y positions of last Forward TOF layer
	bool ftof_back_xy = false;

	//Set collection Ids
	rec_coll_ids = get_coll_ids(rec_coll_names);

	//Create map between RecHit collection IDs and names
    	std::unordered_map<unsigned int, std::string> RecHitCollMap;
    	for (int i = 0; i < rec_coll_ids.size(); i++) {
		RecHitCollMap[rec_coll_ids[i]] = rec_coll_names[i];
    	}

	//Print out Map information
	cout << endl << "Created Map between Collection IDs and Names:"<<endl;
	for (const auto& pair : RecHitCollMap) {
        	cout << "Collection ID: "<< pair.first << " | Collection Name: " << pair.second << endl;
	}

	//Create map between RecHit collection names and vector index
	std::unordered_map<std::string, int> index_map;
		for (size_t i = 0; i < rec_coll_names.size(); i++) {
        	index_map[rec_coll_names[i]] = i;
    	}

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
	TH1 *h1 = new TH1D("h1","Number of SimHits associated w/ primary particle",10,0,10);
	h1->SetLineWidth(3);h1->SetLineColor(kBlue);
	h1->GetXaxis()->SetTitle("Number of hits");h1->GetXaxis()->CenterTitle();
	h1->GetXaxis()->CenterLabels();
	h1->GetYaxis()->SetTitle("Number of tracks");h1->GetYaxis()->CenterTitle();

	TH1 *h2 = new TH1D("h2","Fraction of associated SimHits that survive 'digitization'",25,0,1.1);
	h2->SetLineWidth(3);h2->SetLineColor(kBlue);
	h2->GetXaxis()->SetTitle("Fraction");h2->GetXaxis()->CenterTitle();
	h2->GetYaxis()->SetTitle("Number of tracks");h2->GetYaxis()->CenterTitle();

	TH1 *h3 = new TH1D("h3","Fraction of track measurement hits that are associated w/ primary particle",25,0,1.1);
	h3->SetLineWidth(3);h3->SetLineColor(kBlue);
	h3->GetXaxis()->SetTitle("Fraction");h3->GetXaxis()->CenterTitle();
	h3->GetYaxis()->SetTitle("Number of tracks");h3->GetYaxis()->CenterTitle();

	TH1 *h4 = new TH1D("h4","Fraction of associated digitized hits classified as measurement hits",25,0,1.1);
	h4->SetLineWidth(3);h4->SetLineColor(kBlue);
	h4->GetXaxis()->SetTitle("Fraction");h4->GetXaxis()->CenterTitle();
	h4->GetYaxis()->SetTitle("Number of tracks");h4->GetYaxis()->CenterTitle();

	TH1 *h5 = new TH1D("h5","Fraction of track outlier hits that are associated w/ primary particle",25,0,1.1);
	h5->SetLineWidth(3);h5->SetLineColor(kBlue);
	h5->GetXaxis()->SetTitle("Fraction");h5->GetXaxis()->CenterTitle();
	h5->GetYaxis()->SetTitle("Number of tracks");h5->GetYaxis()->CenterTitle();

	TH1 *h6 = new TH1D("h6","Fraction of associated digitized hits classified as outlier hits",25,0,1.1);
	h6->SetLineWidth(3);h6->SetLineColor(kBlue);
	h6->GetXaxis()->SetTitle("Fraction");h6->GetXaxis()->CenterTitle();
	h6->GetYaxis()->SetTitle("Number of tracks");h6->GetYaxis()->CenterTitle();

	TH1 *h7 = new TH1D("h7","Fraction of associated digitized hits missing from track",25,0,1.1);
	h7->SetLineWidth(3);h7->SetLineColor(kBlue);
	h7->GetXaxis()->SetTitle("Fraction");h7->GetXaxis()->CenterTitle();
	h7->GetYaxis()->SetTitle("Number of tracks");h7->GetYaxis()->CenterTitle();

	//Number of RecHits per collection
	//1) All associated RecHits for events w/ a reconstructed track
	//2) All associated measurement RecHits for events w/ a reconstructed track
	//3) All associated outlier RecHits for events w/ a reconstructed track
	int nBins = rec_coll_names.size();

    	TH1 *h8a = new TH1D("h8a","Events with a reconstructed track", nBins, 0, nBins); //Associated digitized hits
	h8a->SetLineWidth(3);h8a->SetLineColor(kBlue);
	h8a->GetYaxis()->SetTitle("Counts");h8a->GetYaxis()->CenterTitle();

	TH1 *h8b = new TH1D("h8b","Events with a reconstructed track", nBins, 0, nBins); //Associated track measurement hits
        h8b->SetFillColor(kGreen+2);h8b->SetFillStyle(3004);h8b->SetLineWidth(0);
        h8b->GetYaxis()->SetTitle("Counts");h8b->GetYaxis()->CenterTitle();

	TH1 *h8c = new TH1D("h8c","Events with a reconstructed track", nBins, 0, nBins); //Associated outlier hits
        h8c->SetFillColor(kRed);h8c->SetFillStyle(3005);h8c->SetLineWidth(0);
        h8c->GetYaxis()->SetTitle("Counts");h8c->GetYaxis()->CenterTitle();

    	//Assign bin labels
    	for (size_t i = 0; i < rec_coll_names.size(); i++) {
        	h8a->GetXaxis()->SetBinLabel(i + 1, rec_coll_names[i].c_str());
		h8b->GetXaxis()->SetBinLabel(i + 1, rec_coll_names[i].c_str());
		h8c->GetXaxis()->SetBinLabel(i + 1, rec_coll_names[i].c_str());
    	}

	//Number of RecHits per collection
        //1) All associated RecHits
        //2) All associated Measurement2D Hits
        TH1 *h9a = new TH1D("h9a","", nBins, 0, nBins); //Associated RecHits
        h9a->SetLineWidth(3);h9a->SetLineColor(kBlue);
        h9a->GetYaxis()->SetTitle("Counts");h9a->GetYaxis()->CenterTitle();

        TH1 *h9b = new TH1D("h9b","", nBins, 0, nBins); //Associated Measurement2D
        h9b->SetFillColor(kMagenta+2);h9b->SetFillStyle(3004);h9b->SetLineWidth(0);
        h9b->GetYaxis()->SetTitle("Counts");h9b->GetYaxis()->CenterTitle();

	//Assign bin labels
	for (size_t i = 0; i < rec_coll_names.size(); i++) {
                h9a->GetXaxis()->SetBinLabel(i + 1, rec_coll_names[i].c_str());
		h9b->GetXaxis()->SetBinLabel(i + 1, rec_coll_names[i].c_str());
	}

	//R vs. Z: Events with a reconstructed track
	TH2 *hh1a = new TH2D("hh1a","All associated digitized hits",1000,-1800,2000,1000,-800,800);
	hh1a->GetXaxis()->SetTitle("z [mm]");hh1a->GetXaxis()->CenterTitle();
	hh1a->GetYaxis()->SetTitle("signed r [mm]");hh1a->GetYaxis()->CenterTitle();

	TH2 *hh1b = new TH2D("hh1b","Associated track measurement hits",1000,-1800,2000,1000,-800,800);
	hh1b->GetXaxis()->SetTitle("z [mm]");hh1b->GetXaxis()->CenterTitle();
	hh1b->GetYaxis()->SetTitle("signed r [mm]");hh1b->GetYaxis()->CenterTitle();

	TH2 *hh1c = new TH2D("hh1c","Associated track outlier hits",1000,-1800,2000,1000,-800,800);
	hh1c->GetXaxis()->SetTitle("z [mm]");hh1c->GetXaxis()->CenterTitle();
	hh1c->GetYaxis()->SetTitle("signed r [mm]");hh1c->GetYaxis()->CenterTitle();

	TH2 *hh1d = new TH2D("hh1d","Associated digitized hits missing from track",1000,-1800,2000,1000,-800,800);
	hh1d->GetXaxis()->SetTitle("z [mm]");hh1d->GetXaxis()->CenterTitle();
	hh1d->GetYaxis()->SetTitle("signed r [mm]");hh1d->GetYaxis()->CenterTitle();	

	//Forward TOF (last layer) y vs. x: Events with a reconstructed track
	TH2 *hh2a = new TH2D("hh2a","All associated digitized hits",1000,-600,600,1000,-600,600);
        hh2a->GetXaxis()->SetTitle("x [mm]");hh2a->GetXaxis()->CenterTitle();
        hh2a->GetYaxis()->SetTitle("y [mm]");hh2a->GetYaxis()->CenterTitle();

        TH2 *hh2b = new TH2D("hh2b","Associated track measurement hits",1000,-600,600,1000,-600,600);
        hh2b->GetXaxis()->SetTitle("x [mm]");hh2b->GetXaxis()->CenterTitle();
        hh2b->GetYaxis()->SetTitle("y [mm]");hh2b->GetYaxis()->CenterTitle();

        TH2 *hh2c = new TH2D("hh2c","Associated track outlier hits",1000,-600,600,1000,-600,600);
        hh2c->GetXaxis()->SetTitle("x [mm]");hh2c->GetXaxis()->CenterTitle();
        hh2c->GetYaxis()->SetTitle("y [mm]");hh2c->GetYaxis()->CenterTitle();

        TH2 *hh2d = new TH2D("hh2d","Associated digitized hits missing from track",1000,-600,600,1000,-600,600);
        hh2d->GetXaxis()->SetTitle("x [mm]");hh2d->GetXaxis()->CenterTitle();
        hh2d->GetYaxis()->SetTitle("y [mm]");hh2d->GetYaxis()->CenterTitle();	


	//Open file	
	podio::ROOTReader r;
	r.openFile(input_file);

	auto nevents = r.getEntries(podio::Category::Event);
	cout<<"---------------"<<endl;
	cout<<"Total number of events = "<<nevents<<"!"<<endl;
	cout<<"---------------"<<endl; 

	//Loop over events. Can I do this more easily?
	for( unsigned int ievent = 0; ievent < nevents; ievent++){

		if((ievent%1000)==0) cout<<"Processed event "<<ievent<<endl;

		//Get event
		auto f = podio::Frame(r.readNextEntry(podio::Category::Event));
		if(print_evt_info)
			cout<<"For event "<<ievent<<":"<<endl;

		//RawHit / SimHit association collection
		auto& raw_hit_assocs = f.get<edm4eic::MCRecoTrackerHitAssociationCollection>("CentralTrackingRawHitAssociations");

		//Total number of SimHits
		int tot_simhits(0);
		//Number of SimHits that are associated with primary particle
		int matched_simhits(0);

		//Loop over Simhit collections
		if(print_evt_info) 
			cout<<"Number of SimHit collections = "<<sim_coll_names.size()<<endl;
		for(int icoll = 0; icoll < sim_coll_names.size(); icoll++){
			//Simhits
			auto coll_name = sim_coll_names.at(icoll);
			auto& simhit_coll = f.get<edm4hep::SimTrackerHitCollection>(coll_name);
			auto num_simhits = simhit_coll.size();

			if(print_evt_info){
				printf("\nNumber of %s = %lu \n",coll_name.c_str(),num_simhits);
				cout<<"Hit Collection ID | Hit Index | Hit CellID | Hit Quality | MC Index | MC PDG ID | MC Status"<<endl;
			}
			for(int ihit = 0; ihit<num_simhits; ihit++){
				auto hit = simhit_coll.at(ihit);
				auto hit_collid = hit.getObjectID().collectionID;
				auto hit_index = hit.getObjectID().index;
				auto hit_cellid = hit.getCellID();
				auto quality = hit.getQuality();

				//MCParticle associated with simhit
				auto hit_mcpart = hit.getMCParticle();
				auto mc_index = hit_mcpart.getObjectID().index;
				auto mc_pdg = hit_mcpart.getPDG();
				auto mc_status = hit_mcpart.getGeneratorStatus();

				if(print_evt_info)
					printf("%u | %d | %lu | %d | %d| %d | %d \n",hit_collid, hit_index, hit_cellid, quality, mc_index, mc_pdg, mc_status);

				if(mc_status==1 && quality==0) matched_simhits++;
				tot_simhits++;

				//For new Barrel MPGDs digitization, a single SimHit usually is converted into 2 RecHits
				//So, we add an additional count for those to keep the numbers consistent
				if( new_mpgd_digi && (coll_name=="MPGDBarrelHits" || coll_name=="OuterMPGDBarrelHits") ){
					if(mc_status==1 && quality==0) matched_simhits++;
                                	tot_simhits++;
				}
			} //Loop over SimHits in a given collection
		} //Loop over SimHit collections

		//-------------
		//Total number of RecHits associated to primary particle
		int tot_rechits(0);
		
		//Loop over RecHit collections
		if(print_evt_info){
			cout<<endl<<endl;
			cout<<"Number of RecHit collections = "<<rec_coll_names.size()<<endl;
		}
		for(int icoll = 0; icoll < rec_coll_names.size(); icoll++){
			//RecHits
			auto coll_name = rec_coll_names.at(icoll);
			auto& rechit_coll = f.get<edm4eic::TrackerHitCollection>(coll_name);
			auto num_rechits = rechit_coll.size();

			if(print_evt_info){
				printf("\nNumber of %s = %lu \n",coll_name.c_str(),num_rechits);
				cout<<"Hit Collection ID | Hit Index | Hit CellID"<<endl;
			}
			for(int ihit = 0; ihit<num_rechits; ihit++){

				auto hit = rechit_coll.at(ihit);
				auto hit_collid = hit.getObjectID().collectionID;
				auto hit_index = hit.getObjectID().index;
				auto hit_cellid = hit.getCellID();

				//RecHit positions
				auto hit_x = hit.getPosition().x;
				auto hit_y = hit.getPosition().y;
				auto hit_z = hit.getPosition().z;
				double hit_radius = sgn<double>(hit_x) * sqrt(hit_x*hit_x + hit_y*hit_y);
				
				if(print_evt_info)
					printf("%u | %d | %lu \n",hit_collid, hit_index, hit_cellid);

				//Find corresponding SimHit and check if it is associated w/ primary particle
				auto raw_hit = hit.getRawHit();

				for (const auto raw_hit_assoc : raw_hit_assocs) {
					if (raw_hit_assoc.getRawHit() == raw_hit) {
						auto sim_hit = raw_hit_assoc.getSimHit();
                        			auto mc_particle = sim_hit.getMCParticle();
						auto mc_status = sim_hit.getMCParticle().getGeneratorStatus();
						auto quality = sim_hit.getQuality();

						if(mc_status==1 && quality==0){ 
							tot_rechits++;
						
							//If there is at least 1 track, fill histogram of collection RecHits
							auto& track_coll = f.get<edm4eic::TrackCollection>(trk_coll);
                					auto num_tracks = track_coll.size();

							if(num_tracks > 0 && (index_map.find(coll_name) != index_map.end()) ){
								h8a->Fill( index_map[coll_name] );
								hh1a->Fill(hit_z,hit_radius);

								if(coll_name=="TOFEndcapRecHits" && hit_z>1860.)
									hh2a->Fill(hit_x,hit_y);
							}
							
							//Without track requirement
							if( index_map.find(coll_name) != index_map.end() ){
                                                                h9a->Fill( index_map[coll_name] );
                                                        }
			
						} //If hit associated to primary particle
					}
				} //Loop over RawHit/SimHit associations 
			} //Loop over RecHits in a given collection
		} //Loop over RecHit collections

		//-------------
		//Loop over Measurement2D hits (independent of track reconstruction)
		auto& meas2d_coll = f.get<edm4eic::Measurement2DCollection>("CentralTrackerMeasurements");
		for(int imeas = 0; imeas<meas2d_coll.size(); imeas++){
			
			auto rec_hit = meas2d_coll.at(imeas).getHits().at(0); //One RecHit exists for a given Measurement2D
                        
			auto cellid = rec_hit.getCellID();
                        auto collid = rec_hit.getObjectID().collectionID;

			//Get collection name for RecHit
                        std::string coll_name;
                        if (RecHitCollMap.find(collid) != RecHitCollMap.end())
                        	coll_name = RecHitCollMap.at(collid);
                        else
                        	coll_name = "UnknownCollection";

			//Find corresponding SimHit and check if it is associated w/ primary particle
                        auto raw_hit = rec_hit.getRawHit();

                        for (const auto raw_hit_assoc : raw_hit_assocs) {
                        	if (raw_hit_assoc.getRawHit() == raw_hit) {
                                	auto sim_hit = raw_hit_assoc.getSimHit();
                                       	auto mc_particle = sim_hit.getMCParticle();
                                        auto mc_status = sim_hit.getMCParticle().getGeneratorStatus();
                                        auto quality = sim_hit.getQuality();
					
                                        if(mc_status==1 && quality==0){
                                        	if( index_map.find(coll_name) != index_map.end() ){
                                                	h9b->Fill( index_map[coll_name] );
                                                }
                                        } //If hit associated to primary particle
                        	}
			} //Loop over RawHit/SimHit associations
		} //Loop over Measurement2D points

		//-------------
		//RecHits that are measurement/outlier hits of a given track
		auto& track_coll = f.get<edm4eic::TrackCollection>(trk_coll);
		auto num_tracks = track_coll.size();

		//Will want to adjust this for events with multiple tracks (not single particle)
		int num_meas(0); //number of measurements
		int num_out(0); //number of outliers
		int num_meas_corr(0); //number of correctly-identified measurements
		int num_out_wrong(0); //number of wrongly-identified outliers

		//Print track and associated measurement RecHit info
		if(print_evt_info)
			cout<<endl<<endl<<"Number of tracks = "<<num_tracks<<endl;	

		for(int itrack = 0; itrack<num_tracks; itrack++){
			
			auto track = track_coll.at(itrack); //Track

			//-------------
			//Get measurements
			if(print_evt_info)
				cout<<"Track Index | Measurement number | Measurement CellID | Collection | Should be measurement?"<<endl;
			auto meas2d = track.getMeasurements(); //Associated Measurement2D array
			for( int imeas = 0; imeas < meas2d.size(); imeas++){
				
				auto rec_hit = meas2d.at(imeas).getHits().at(0); //One RecHit exists for a given Measurement2D 

				auto cellid = rec_hit.getCellID();
				auto collid = rec_hit.getObjectID().collectionID;
		
				//RecHit positions
                        	auto hit_x = rec_hit.getPosition().x;
                        	auto hit_y = rec_hit.getPosition().y;
                        	auto hit_z = rec_hit.getPosition().z;
                        	double hit_radius = sgn<double>(hit_x) * sqrt(hit_x*hit_x + hit_y*hit_y);
	
				//Get collection name for RecHit	
				std::string coll_name;
				if (RecHitCollMap.find(collid) != RecHitCollMap.end())
					coll_name = RecHitCollMap.at(collid);
				else
					coll_name = "UnknownCollection";

				//Find associated SimHit and see if hit should be used in track fit
				std::string should_be_meas = "N";
				auto raw_hit = rec_hit.getRawHit();

				for (const auto raw_hit_assoc : raw_hit_assocs) {
					if (raw_hit_assoc.getRawHit() == raw_hit) {
						auto sim_hit = raw_hit_assoc.getSimHit();
                        			auto mc_particle = sim_hit.getMCParticle();
						auto mc_status = sim_hit.getMCParticle().getGeneratorStatus();
						auto quality = sim_hit.getQuality();

						if(mc_status==1 && quality==0){
							should_be_meas = "Y";
							num_meas_corr++;
				
							//Fill associated measurement hit collection histogram
							if( (index_map.find(coll_name) != index_map.end()) ){
                                                                h8b->Fill( index_map[coll_name] );
								hh1b->Fill(hit_z,hit_radius);

								if(coll_name=="TOFEndcapRecHits" && hit_z>1860.)
                                                                        hh2b->Fill(hit_x,hit_y);
                                                        }

						} //If hit associated to primary particle
					}
				} //Loop over RawHit/SimHit associations 

				if(print_evt_info)
					printf("%d | %d | %lu | %s | %s \n",itrack, imeas, cellid, coll_name.c_str(), should_be_meas.c_str());
				num_meas++;
			}

			//-------------
			//Get outliers (currently stored in Trajectory)
			if(print_evt_info)
				cout<<"Track Index | Outlier number | Measurement CellID | Collection | Should be measurement?"<<endl;
			auto traj = track.getTrajectory();
			auto out2d = traj.getOutliers_deprecated();

			for( int iout = 0; iout < out2d.size(); iout++){
				
				auto out_hit = out2d.at(iout).getHits().at(0); //One RecHit exists for a given Measurement2D 

				auto cellid = out_hit.getCellID();
				auto collid = out_hit.getObjectID().collectionID;

				//RecHit positions
                                auto hit_x = out_hit.getPosition().x;
                                auto hit_y = out_hit.getPosition().y;
                                auto hit_z = out_hit.getPosition().z;
                                double hit_radius = sgn<double>(hit_x) * sqrt(hit_x*hit_x + hit_y*hit_y);
				
                                //Get collection name for RecHit
                                std::string coll_name;
                                if (RecHitCollMap.find(collid) != RecHitCollMap.end())
                                        coll_name = RecHitCollMap.at(collid);
                                else
                                        coll_name = "UnknownCollection";

				//Find associated SimHit and see if hit should be used in track fit
				std::string should_be_meas = "N";
				auto raw_hit = out_hit.getRawHit();

				for (const auto raw_hit_assoc : raw_hit_assocs) {
					if (raw_hit_assoc.getRawHit() == raw_hit) {
						auto sim_hit = raw_hit_assoc.getSimHit();
                        			auto mc_particle = sim_hit.getMCParticle();
						auto mc_status = sim_hit.getMCParticle().getGeneratorStatus();
						auto quality = sim_hit.getQuality();

						if(mc_status==1 && quality==0){
							should_be_meas = "Y";
							num_out_wrong++;

							//Fill associated outlier hit collection histogram
                                                        if( (index_map.find(coll_name) != index_map.end()) ){
                                                                h8c->Fill( index_map[coll_name] );
								hh1c->Fill(hit_z,hit_radius);

								if(coll_name=="TOFEndcapRecHits" && hit_z>1860.)
                                                                        hh2c->Fill(hit_x,hit_y);
                                                        }
						}
					}
				} //Loop over RawHit/SimHit associations

				if(print_evt_info)
					printf("%d | %d | %lu | %s | %s \n",itrack, iout, cellid, coll_name.c_str(), should_be_meas.c_str());
				num_out++;
			}
			
		} //End loop over tracks

		//Print some event statistics
		if(print_evt_info){
			cout<<endl<<"Event statistics:"<<endl;
			cout<<"SimHits (Barrel MPGD corrected) associated with primary particle = "<<matched_simhits<<" / "<<tot_simhits<<endl;
			cout<<"RecHits associated with primary particle = "<<tot_rechits<<endl;
			cout<<"Number of correctly-identified measurement hits = "<<num_meas_corr<<" / "<<num_meas<<endl;
			cout<<"Number of outlier hits that should be measurement hits = "<<num_out_wrong<<" / "<<num_out<<endl;
		}

		//Total number of hits (meas. + out.) in track that are associated to a primary-particle hit
		auto track_assoc_hit_tot = num_meas_corr+num_out_wrong;
		//Total number of digitized hits (RecHits) associated with primary particle completely missing from track
		auto simhits_missing = matched_simhits - track_assoc_hit_tot;
		auto rechits_missing = simhits_missing - ( matched_simhits - tot_rechits ); //Don't count any hits lost during digitization
		
		if(print_evt_info) 
			cout<<"Number of RecHits associated with primary particle completely missing from track = "<<rechits_missing<<" / "<<tot_rechits<<endl;

		//Add line between events
		if(print_evt_info)
			cout<<endl<<"---------------"<<endl;

		//Fill histograms
		if(num_tracks>0){ //Require at least 1 track
			h1->Fill(matched_simhits);
			if(matched_simhits>0) 
				h2->Fill( ((float) tot_rechits) / matched_simhits );
			if(num_meas>0) 
				h3->Fill( ((float) num_meas_corr) / num_meas );
			if(tot_rechits>0) 
				h4->Fill( ((float) num_meas_corr) / tot_rechits );
			if(num_out>0) 
				h5->Fill( ((float) num_out_wrong) / num_out );
			if(tot_rechits>0) 
				h6->Fill( ((float) num_out_wrong) / tot_rechits );
			if(tot_rechits>0) 
				h7->Fill( ((float) rechits_missing) / tot_rechits );
		}

	} //End loop over events

	//Scale and Subtract histograms
	TH1* h8bs = (TH1D*)h8b->Clone("h8bs");
	h8bs->GetYaxis()->SetTitle("Fraction of associated digitized hits");
	h8bs->Divide(h8a);

	TH1* h8cs = (TH1D*)h8c->Clone("h8cs");
	h8cs->GetYaxis()->SetTitle("Fraction of associated digitized hits");
        h8cs->Divide(h8a);

	hh1d->Add(hh1a,hh1b,1,-1); //Subtract track measurement hits from associated hits
	hh1d->Add(hh1c,-1); //Also subtract outlier hits

	hh2d->Add(hh2a,hh2b,1,-1); //Subtract track measurement hits from associated hits
        hh2d->Add(hh2c,-1); //Also subtract outlier hits

	//Histogram stacking
	THStack *stack1 = new THStack("stack1", "Stacked Histograms");
    	stack1->Add(h8b);
    	stack1->Add(h8c);

	THStack *stack2 = new THStack("stack2", "Events with a reconstructed track; ; Fraction of associated digitized hits");
        stack2->Add(h8bs);
        stack2->Add(h8cs);

	//Make plots
	TCanvas *c1 = new TCanvas("c1");
	h1->Draw();

	TLatex *tex1 = new TLatex();
	tex1->SetTextSize(0.03);
	//tex1->DrawLatexNDC(0.2,0.7,"Events with a reconstructed truth-seeded track"); //Uncomment if using truth-seeded tracks
	tex1->DrawLatexNDC(0.2,0.7,"Events with a reconstructed real-seeded track");

	TCanvas *c2 = new TCanvas("c2");
	h2->Draw();

	TCanvas *c3 = new TCanvas("c3");
	h3->Draw();

	TCanvas *c4 = new TCanvas("c4");
	h4->Draw();

	TCanvas *c5 = new TCanvas("c5");
	h5->Draw();

	TCanvas *c6 = new TCanvas("c6");
	h6->Draw();

	TCanvas *c7 = new TCanvas("c7");
	h7->Draw();

	TCanvas *c8 = new TCanvas("c8");
        h8a->Draw(); 
    	gPad->Update(); // Update the pad to set axis range correctly
	stack1->Draw("HIST same");

	TLegend *legend = new TLegend(0.4, 0.75, 0.9, 0.9); 
	legend->AddEntry(h8a,"Associated Digitized Hits (RecHits)","l");
	legend->AddEntry(h8b,"Associated Tracker Measurement Hits", "f");
   	legend->AddEntry(h8c,"Associated Tracker Outlier Hits", "f");
	legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.035);
	legend->Draw();

	TCanvas *c8s = new TCanvas("c8s");
	c8s->SetGrid();
	stack2->Draw("HIST");

	TCanvas *c9 = new TCanvas("c9");
	h9a->Draw();
	h9b->Draw("same");

	TLegend *legend2 = new TLegend(0.4, 0.75, 0.9, 0.9);
        legend2->AddEntry(h9a,"Associated Digitized Hits (RecHits)","l");
        legend2->AddEntry(h9b,"Associated Measurement2D Hits", "f");
        legend2->SetBorderSize(0);
        legend2->SetFillStyle(0);
        legend2->SetTextSize(0.035);
        legend2->Draw();

	TCanvas *c10a = new TCanvas("c10a");
	hh1a->Draw();

	TCanvas *c10b = new TCanvas("c10b");
        hh1b->Draw();

	TCanvas *c10c = new TCanvas("c10c");
        hh1c->Draw();

	TCanvas *c10d = new TCanvas("c10d");
        hh1d->Draw();

	TCanvas *c11a = new TCanvas("c11a");
        hh2a->Draw();

	TCanvas *c11b = new TCanvas("c11b");
        hh2b->Draw();

        TCanvas *c11c = new TCanvas("c11c");
        hh2c->Draw();

        TCanvas *c11d = new TCanvas("c11d");
        hh2d->Draw();

	//Print plots to file
	cout<<endl;
	c1->Print("plots/hit_matching.pdf[");
	c1->Print("plots/hit_matching.pdf");
	c2->Print("plots/hit_matching.pdf");
	c3->Print("plots/hit_matching.pdf");
	c4->Print("plots/hit_matching.pdf");
	c5->Print("plots/hit_matching.pdf");
	c6->Print("plots/hit_matching.pdf");
	c7->Print("plots/hit_matching.pdf");
	c8->Print("plots/hit_matching.pdf");
	c8s->Print("plots/hit_matching.pdf");
	c9->Print("plots/hit_matching.pdf");
	c10a->Print("plots/hit_matching.pdf");
	c10b->Print("plots/hit_matching.pdf");
	c10c->Print("plots/hit_matching.pdf");
	c10d->Print("plots/hit_matching.pdf");
	c11a->Print("plots/hit_matching.pdf");
	c11b->Print("plots/hit_matching.pdf");
	c11c->Print("plots/hit_matching.pdf");
	c11d->Print("plots/hit_matching.pdf");
	c11d->Print("plots/hit_matching.pdf]");
}
