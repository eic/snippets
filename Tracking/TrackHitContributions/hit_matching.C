// Histogram hits, inclusively and in association to tracks.
// - Input = podio output of EICrecon, symlinked to "hit_matching.input.root".
// - MPGDs may be 2DStrip instead of pixel.
//  - If indeed, they are to be of the Multiple Sensitive Volume type (see
//   "EICrecon/src/algorithms/digi/MPGDTrackerDigi.cc").
//  - Whether they are processed as 2DStrip, requires:
//   i) A TGeometry file in the $PWD, named "hit_matching.geometry.root". File
//    to be produced by "dd_web_display". E.g.:
//      dd_web_display -o hit_matching.geometry.root \
//        --export $DETECTOR_PATH/$DETECTOR_CONFIG.xml
//  ii) The parsing of the TGeometry file, on a per MPGD basis, determines that
//    we are dealing w/ a MSV.
//  - When (i) and (ii) are fulfilled...
//  ...Strips (named 'p' and 'n') are histo'd independently, as well as 'p&n'.
//  ...3D info is obtained by combining the meaningful piece of 'p' and 'n'
//    RecHits...
//  ...This, possibly after having transformed them from WorldRS to LocalRS (
//    case of Outer Barrel): a dedicated "Geometry" class is used. It's quite
//    general, providing for World<->Local for any MPGD (and more if extended).
//  - NOTA BENE: this 2DStrip processing is valid only if the geometry of the
//   TGeometry file matches that of the data. Which matching is the
//   responsibility of the user.

#include <podio/ObjectID.h>
#include <podio/ROOTReader.h>
#include <podio/Frame.h>
#include <edm4hep/Vector3f.h>
#include <edm4hep/SimTrackerHitCollection.h>
#include <edm4eic/TrackerHitCollection.h>
#include <edm4eic/MCRecoTrackerHitAssociationCollection.h>
#include <edm4eic/TrackCollection.h>
#include <edm4eic/Measurement2DCollection.h>

#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TLegend.h>
#include <THStack.h>
#include <TMath.h>
#include <TFile.h>

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
std::string input_file = "hit_matching.input.root";

//------------------
//Template to access 'sign' of radius in (x,y) plane
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

//-----------------
//Function to set Collection IDs
vector<unsigned int> get_coll_ids(vector<string> coll_names){

	vector<unsigned int> coll_ids;

	podio::ROOTReader reader;
	reader.openFile(input_file.c_str());

	unsigned nEntries = reader.getEntries(podio::Category::Event);

	//Find the associated collection IDs
	for (unsigned i = 0; i < nEntries; ++i) {

		if(i > 0) break;

    		auto frameData = reader.readEntry(podio::Category::Event, i);
    		auto frame     = podio::Frame(std::move(frameData));

    		auto collectionNames = frame.getAvailableCollections();

		// Find associated collection ID for each track collection
		for(int iname = 0; iname < (int)coll_names.size(); iname++) {

			for (const auto& name : collectionNames) {
				const auto* coll = frame.get(name);
      				if (coll) {
        				unsigned int collectionID = coll->getID();
					if( name == coll_names[iname] )
					coll_ids.push_back(collectionID);
				}
			}

		} //Loop over track collections
	} //Loop over events

	return coll_ids;
}

//-----------------
// Declaration of "index_map", to be defined infra
std::unordered_map<std::string, int> index_map;

//-----------------
//Geometry class
#include <TGeoManager.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
class Geometry {
public:
  Geometry(int verbose);
  bool is2DStripMPGD(std::string coll_name);
  bool getGeoMats(int mpgd, bool &isFullyMSV);
  bool WorldToLocal(int idet, unsigned long cellID,
		    double *gpos,  // In mm
		    double *lpos); // In mm
  bool LocalToWorld(int idet, unsigned long cellID,
		    double *lpos,  // In mm
		    double *gpos); // In mm
  // ***** SINGLETON CLASS
  static inline Geometry* address;
  static Geometry* Ptr()
  {
    if (address) {
      return(address);
    }
    else {
      printf("** Geometry::Ptr(): Inconsistency: The object of Geometry class is not yet created or already destructed\n");
      return 0;
    }
  };
  static constexpr int N_MPGDs = 4;
  unsigned int MPGD2DStrips;
private:
  vector<vector<TGeoHMatrix>> geoDetMats[N_MPGDs];
  int verboseLevel;
  // Tags allowing to find MPGD layers and modules
  static constexpr const char *tag1s[N_MPGDs] =
    {"Inner","Outer","Endcap",  "Endcap"};
  static constexpr const char *tag2s[N_MPGDs] =
    {0,      0,      "Backward","Forward"};
};

//-----------------
//MPGDs 2DStrip:
int is2DStripMPGD(std::string coll_name)
{
  int mpgd;
  if      (coll_name=="MPGDBarrelRecHits")         mpgd = 0;
  else if (coll_name=="OuterMPGDBarrelRecHits")    mpgd = 1;
  else if (coll_name=="BackwardMPGDEndcapRecHits") mpgd = 2;
  else if (coll_name=="ForwardMPGDEndcapRecHits")  mpgd = 3;
  else return -1;
  if (mpgd>=0) {
    Geometry *geometry = Geometry::Ptr();
    if ((geometry->MPGD2DStrips)&0x1<<mpgd) return mpgd;
  }
  return -1;
}
// Bit specifying strip type in cellID:
using CellID = std::uint64_t;
const CellID pStripBit = ((std::int64_t)0x1)<<28;
const CellID nStripBit = ((std::int64_t)0x2)<<28;
// Binning:
// - Based on detector's index (="index_map"),
// - 2DStrip MPGD case: subdivide into THREE sub-bins:
//  index_map +0/3 for pStrips ('p'),
//  index_map +1/3 for nStrips ('n'),
//  index_map +2/3 for the AND ('p&n').
float getIndex(std::string& coll_name, std::int64_t cellid, int &pn)
{
  // - Returns floating point index, encoding either a fraction (2DStrip case)
  //  or an integer (all other cases).
  //   Returns <0 index, if input "coll_name" unknown
  // - In the 2DStrip case, "pn" encodes upon return the strip coordinate type,
  //  based on the "strip" field of "cellid".
  //   If the latter is invalid, again a <0 index is returned.
  float index = -1;
  if (index_map.find(coll_name) != index_map.end()) {
    index = index_map[coll_name]; pn = 0;
    if (is2DStripMPGD(coll_name)>=0) { // 2DStrip MPGD case
      if      (  cellid&nStripBit) { index += 1/3.; pn = 1; }
      else if (!(cellid&pStripBit)) {
	static unsigned int warnings = 0; unsigned int bit = 0x1<<int(index);
	if (!(warnings&bit)) { // No strip bit => pixelized
	  // Warning: only once to not clutter stdout.y
	  printf("\"%s\" Hit 0x%lx has not strip bit\n",
		 coll_name.c_str(),cellid);
	  warnings |= bit;
	}
	pn = 3; index = -1;
      }
    }
  }
  return index + /* to make for rounding approx. */ .1;
}
//Book keeping of MPGD strip hits
using SimID2RecIDs = std::map<CellID,vector<CellID>>;
SimID2RecIDs::iterator pANDn(CellID simID, CellID recID, SimID2RecIDs &simID2RecIDs)
{
  // Workflow is such that cases where strip field is neither 'p' nor 'n' are excluded
  int pn = (recID&pStripBit) ? 0 : 1; unsigned int bpn = 0x1<<pn;
  SimID2RecIDs::iterator is = simID2RecIDs.find(simID);
  if (is == simID2RecIDs.end()) {
    vector<CellID> recIDs; recIDs.push_back(recID);
    simID2RecIDs[simID] = recIDs;
  } else {
    vector<CellID> &recIDs = is->second; recIDs.push_back(recID);
    if (recIDs.size()==2) {
      CellID prvID = recIDs[0];
      int prvpn = (prvID&pStripBit) ? 0 : 1; unsigned int prvB = 0x1<<prvpn;
      if ((prvB|bpn) == 0x3) // Both 'p' and 'n', successively
	return is;
    }
    //else
    // Size =1: incomplete set of RecHits
    // Size >2: ambiguous case
  }
  return simID2RecIDs.end();
}
bool get2DStripRZ(const string coll_name, const edm4eic::TrackerHitCollection& rechit_coll, SimID2RecIDs::iterator im,
		  double &hit_radius, float &hit_z);

//-----------------
// Main function
void hit_matching(int verbose = 0, int nEvtMx = 0){

	//Track collection to use: real-seeded or truth-seeded
	std::string trk_coll = "CentralCKFTracks";

	//Print out information event-by-event
	int print_evt_info = verbose;

	//Plot x-y positions of last Forward TOF layer
	bool ftof_back_xy = false;

	//Set collection Ids
	rec_coll_ids = get_coll_ids(rec_coll_names);

	//Create map between RecHit collection IDs and names
    	std::unordered_map<unsigned int, std::string> RecHitCollMap;
	for (int i = 0; i < (int)rec_coll_ids.size(); i++) {
		RecHitCollMap[rec_coll_ids[i]] = rec_coll_names[i];
    	}

	//Print out Map information
	cout << endl << "Created Map between Collection IDs and Names:"<<endl;
	for (const auto& pair : RecHitCollMap) {
        	cout << "Collection ID: "<< pair.first << " | Collection Name: " << pair.second << endl;
	}

	//Create map between RecHit collection names and vector index...
	//...making room for two indices for the 2DStrip MPGDs.
	size_t nRecColls = rec_coll_names.size();
	for (size_t i = 0; i < nRecColls; i++) {
		index_map[rec_coll_names[i]] = i;
	}

	//Geometry
	try {
		new Geometry(verbose);
	} catch (const std::runtime_error& error) {
		cerr << error.what() << endl;
		return;
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

	//Binning:
	int nBins = nRecColls;
	// 3 bins for each of the 2DStrip MPGDs: one per strip + one for p&n
	for (size_t i = 0; i < nRecColls; i++) {
          if (is2DStripMPGD(rec_coll_names[i])>=0) nBins += 2;
	}
	double *bins = new double[nBins+1];
	size_t idx, jdx; double bin; for (idx=jdx = 0, bin = 0; idx < nRecColls; idx++, jdx++) {
          string coll_name = rec_coll_names[idx];
          if (is2DStripMPGD(coll_name)>=0) {
            bins[jdx] = bin; bin += 1/3.;
            jdx++; bins[jdx] = bin; bin += 1/3.;
            jdx++; bins[jdx] = bin; bin += 1/3.;
          }
          else {
            bins[jdx] = bin; bin += 1;
          }
	}
	bins[jdx] = bin;

	TH1 *h8a = new TH1D("h8a","Events with a reconstructed track", nBins, bins); //Associated digitized hits
	h8a->SetLineWidth(3);h8a->SetLineColor(kBlue);
	h8a->GetXaxis()->SetTickLength(0);
	h8a->GetYaxis()->SetTitle("Counts");h8a->GetYaxis()->CenterTitle();
	h8a->GetYaxis()->SetMaxDigits(3);

	TH1 *h8b = new TH1D("h8b","Events with a reconstructed track", nBins, bins); //Associated track measurement hits
        h8b->SetFillColor(kGreen+2);h8b->SetFillStyle(3004);h8b->SetLineWidth(0);
        h8b->GetYaxis()->SetTitle("Counts");h8b->GetYaxis()->CenterTitle();

	TH1 *h8c = new TH1D("h8c","Events with a reconstructed track", nBins, bins); //Associated outlier hits
        h8c->SetFillColor(kRed);h8c->SetFillStyle(3005);h8c->SetLineWidth(0);
        h8c->GetYaxis()->SetTitle("Counts");h8c->GetYaxis()->CenterTitle();

	//Number of RecHits per collection
        //1) All associated RecHits
        //2) All associated Measurement2D Hits
        TH1 *h9a = new TH1D("h9a","", nBins, bins); //Associated RecHits
        h9a->SetLineWidth(3);h9a->SetLineColor(kBlue);
        h9a->GetXaxis()->SetTickLength(0);
        h9a->GetYaxis()->SetTitle("Counts");h9a->GetYaxis()->CenterTitle();
        h9a->GetYaxis()->SetMaxDigits(3);

        TH1 *h9b = new TH1D("h9b","", nBins, bins); //Associated Measurement2D
        h9b->SetFillColor(kMagenta+2);h9b->SetFillStyle(3004);h9b->SetLineWidth(0);
        h9b->GetYaxis()->SetTitle("Counts");h9b->GetYaxis()->CenterTitle();

	// Assign bin labels
	for (idx=jdx = 0; idx < nRecColls; idx++) {
          string coll_name = rec_coll_names[idx];
          int stripMPGD = is2DStripMPGD(coll_name);
          int pnMx = stripMPGD>=0 ? 3 : 1;
          for (int pn = 0; pn<pnMx; pn++, jdx++) {
            if (stripMPGD>=0) {
              switch (pn) {
              case 0: coll_name += string("-p");   break;
              case 1: coll_name =  string("-n");   break;
              case 2: coll_name =  string("-p AND -n"); break;
              }
            }
            h8a->GetXaxis()->SetBinLabel(jdx + 1, coll_name.c_str());
            h8b->GetXaxis()->SetBinLabel(jdx + 1, coll_name.c_str());
            h8c->GetXaxis()->SetBinLabel(jdx + 1, coll_name.c_str());
            h9a->GetXaxis()->SetBinLabel(jdx + 1, coll_name.c_str());
            h9b->GetXaxis()->SetBinLabel(jdx + 1, coll_name.c_str());
          }
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
	for( unsigned int ievent = 0; ievent < (nEvtMx ? nEvtMx : nevents); ievent++){
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

		//Loop over SimHit collections
		if(print_evt_info)
			cout<<"Number of SimHit collections = "<<sim_coll_names.size()<<endl;
		for(int icoll = 0; icoll < (int)sim_coll_names.size(); icoll++){
			//SimHits
			auto coll_name = sim_coll_names.at(icoll);
			auto& simhit_coll = f.get<edm4hep::SimTrackerHitCollection>(coll_name);
			auto num_simhits = simhit_coll.size();

			if(print_evt_info){
				printf("\nNumber of %s = %lu \n",coll_name.c_str(),num_simhits);
				cout<<"Hit Collection ID | Hit Index | Hit CellID | Hit Quality | MC Index | MC PDG ID | MC Status"<<endl;
			}
			for(int ihit = 0; ihit<(int)num_simhits; ihit++){
				auto hit = simhit_coll.at(ihit);
				auto hit_collid = hit.getObjectID().collectionID;
				auto hit_index = hit.getObjectID().index;
				auto hit_cellid = hit.getCellID();
				auto quality = hit.getQuality();

				//MCParticle associated with simhit
				auto hit_mcpart = hit.getParticle();
				auto mc_index = hit_mcpart.getObjectID().index;
				auto mc_pdg = hit_mcpart.getPDG();
				auto mc_status = hit_mcpart.getGeneratorStatus();

				if(print_evt_info)
					printf("%8x | %2d | %08lx,%08lx | %x | %d | %d | %d\n",hit_collid, hit_index, hit_cellid&0xffffffff, hit_cellid>>32, quality, mc_index, mc_pdg, mc_status);

				if(mc_status==1 && quality==0) matched_simhits++;
				tot_simhits++;

				//For 2DStrip MPGDs digitization, a single SimHit (or set thereof)
				// usually is converted into 2 RecHits.
				//=> Add an additional count to keep the numbers, APPROX., consistent
				if(is2DStripMPGD(coll_name)>=0) {
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
		for(int icoll = 0; icoll < (int)rec_coll_names.size(); icoll++){
			//RecHits
			auto coll_name = rec_coll_names.at(icoll);
			auto& rechit_coll = f.get<edm4eic::TrackerHitCollection>(coll_name);
			auto num_rechits = rechit_coll.size();

			int stripMPGD = is2DStripMPGD(coll_name);

			if(print_evt_info){
				printf("\nNumber of %s = %lu \n",coll_name.c_str(),num_rechits);
				cout<<"Hit Collection ID | Hit Index | Hit CellID"<<endl;
			}
			//Book keeping of MPGD strip hits
			SimID2RecIDs simID2RecIDs;
			for(int ihit = 0; ihit<(int)num_rechits; ihit++){

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

				//Find corresponding SimHit(s) and check if it is associated w/ primary particle
				//(Note: A given RecHit can be associated to several SimHits...
				// ...Because of the accumulation of several primary particles.
				// ...For 2DStrip MPGDs (in the MultipleSensitiveVolume implemantation), because a single particle can create several SimHits.)

				//First get RawHit
				//(Note: edm4eic provides for only one RecHit->RawHit, while in the
				//Clustered MPGD case, several RawHits (a cluster of them, that is) are,
				//most of the time, associated to the RecHit. Yet, since all such
				//RawHits are associated to the same SimHits, this does not affect the
				//overall RecHit->SimHits association.)
				auto raw_hit = hit.getRawHit();

				for (const auto raw_hit_assoc : raw_hit_assocs) {
					if (raw_hit_assoc.getRawHit() == raw_hit) {
						auto sim_hit = raw_hit_assoc.getSimHit();
                        			auto mc_particle = sim_hit.getParticle();
						auto mc_status = sim_hit.getParticle().getGeneratorStatus();
						auto quality = sim_hit.getQuality();

						if(mc_status==1 && quality==0){ 
							tot_rechits++;
						
							//If there is at least 1 track, fill histogram of collection RecHits
							auto& track_coll = f.get<edm4eic::TrackCollection>(trk_coll);
                					auto num_tracks = track_coll.size();

							//Get index, adding (1||2)/3 in 2DStrip MPGD case.
							int pn; // p|n of 2DStrip MPGD
							float index = getIndex(coll_name,hit_cellid,pn);
							if (pn == 3) continue; // It's a MPGD w/ no valid strip bit...
							//Fill histogram of collection RecHits
							if (index > 0) { 
                                                          // With and without track requirement
                                                          if(num_tracks > 0) h8a->Fill( index );
                                                          h9a->Fill( index );
                                                          if (stripMPGD<0) { //Histogram RecHit positions
                                                            hh1a->Fill(hit_z,hit_radius);

                                                            if(coll_name=="TOFEndcapRecHits" && hit_z>1860.)
                                                              hh2a->Fill(hit_x,hit_y);
                                                          }
							}

							if (stripMPGD>=0) { // 2DStrip MPGDs: do we have 'p&n'?
                                                          CellID simID = sim_hit.getCellID();
                                                          SimID2RecIDs:: iterator im =
                                                            pANDn(simID,hit_cellid,simID2RecIDs);
                                                          if (im != simID2RecIDs.end()) {
                                                            float pnIndex = int(index) + 2/3. + .1;
                                                            if(num_tracks > 0) {
                                                              h8a->Fill( pnIndex ); //Histogram 'p&n'
                                                              if (get2DStripRZ(coll_name,rechit_coll,im,hit_radius,hit_z))
											hh1a->Fill(hit_z,hit_radius); //Histogram RecHit positions
                                                            }
                                                            h9a->Fill( pnIndex );
                                                            if(print_evt_info)
                                                              printf("%8x | %2d | %08lx,%08lx | %5.0f | %4.0f\n",
                                                                     hit_collid, hit_index, hit_cellid&0xffffffff, hit_cellid>>32,
                                                                     hit_z,hit_radius);
                                                          }
							}
			
						} //If hit associated to primary particle
					}
				} //Loop over RawHit/SimHit associations 
			} //Loop over RecHits in a given collection
		} //Loop over RecHit collections

		//-------------
		//Loop over Measurement2D hits (independent of track reconstruction)
		auto& meas2d_coll = f.get<edm4eic::Measurement2DCollection>("CentralTrackerMeasurements");
		//Book keeping of MPGD strip hits.
		SimID2RecIDs simID2RecIDs;
		for(int imeas = 0; imeas<(int)meas2d_coll.size(); imeas++){

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
                                       	auto mc_particle = sim_hit.getParticle();
                                        auto mc_status = sim_hit.getParticle().getGeneratorStatus();
                                        auto quality = sim_hit.getQuality();
					
					if(mc_status==1 && quality==0){
						//Get index, adding 1||2/3 in 2DStrip MPGD case.
                                                int pn; // p|n of 2DStrip MPGD
                                                float index = getIndex(coll_name,cellid,pn);
                                                if (pn == 3) continue; // It's a MPGD w/ no valid strip bit...
                                                if(index > 0) h9b->Fill( index );

                                                //2DStrip MPGDs
                                                int stripMPGD = is2DStripMPGD(coll_name); if (stripMPGD>=0) {
                                                  // Do we have 'p&n'?
                                                  CellID simID = sim_hit.getCellID();
                                                  SimID2RecIDs:: iterator im =
                                                    pANDn(simID,cellid,simID2RecIDs);
                                                  if (im != simID2RecIDs.end()) {
                                                    float pnIndex = int(index) + 2/3. + .1;
                                                    h9b->Fill(pnIndex);
                                                  }
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

		for(int itrack = 0; itrack<(int)num_tracks; itrack++){
			
			auto track = track_coll.at(itrack); //Track

			//-------------
			//Get measurements
			if(print_evt_info)
				cout<<"Track Index | Measurement number | Measurement CellID | Collection | Should be measurement?"<<endl;
			auto meas2d = track.getMeasurements(); //Associated Measurement2D array
			//Book keeping of MPGD strip hits.
			SimID2RecIDs simID2RecIDs;
			for( int imeas = 0; imeas < (int)meas2d.size(); imeas++){

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

				int stripMPGD = is2DStripMPGD(coll_name);

				//Find associated SimHit and see if hit should be used in track fit
				std::string should_be_meas = "N";
				auto raw_hit = rec_hit.getRawHit();

				for (const auto raw_hit_assoc : raw_hit_assocs) {
					if (raw_hit_assoc.getRawHit() == raw_hit) {
						auto sim_hit = raw_hit_assoc.getSimHit();
                        			auto mc_particle = sim_hit.getParticle();
						auto mc_status = sim_hit.getParticle().getGeneratorStatus();
						auto quality = sim_hit.getQuality();

						if(mc_status==1 && quality==0){
							should_be_meas = "Y";
							num_meas_corr++;
				
							//Get index, adding 1||2/3 in 2DStrip MPGD case,
							int pn; // p|n of 2DStrip MPGD
							float index = getIndex(coll_name,cellid,pn);
							if (pn == 3) continue; // It's a MPGD w/ no valid strip bit...
							//Fill associated measurement hit collection histogram
							if(index > 0) {
                                                          h8b->Fill( index );
                                                          if (stripMPGD<0) { //Histogram RecHit positions
                                                            hh1b->Fill(hit_z,hit_radius);

                                                            if(coll_name=="TOFEndcapRecHits" && hit_z>1860.)
                                                              hh2b->Fill(hit_x,hit_y);
                                                          }
							}

							if (stripMPGD>=0) { // 2DStrip MPGDs: do we have 'p&n'?
                                                          CellID simID = sim_hit.getCellID();
                                                          SimID2RecIDs:: iterator im =
                                                            pANDn(simID,cellid,simID2RecIDs);
                                                          if (im != simID2RecIDs.end()) {
                                                            float pnIndex = int(index) + 2/3. + .1;
                                                            h8b->Fill( pnIndex ); //Histogram "p&n"
                                                            auto& rechit_coll = f.get<edm4eic::TrackerHitCollection>(coll_name);
                                                            if (get2DStripRZ(coll_name,rechit_coll,im,hit_radius,hit_z))
                                                              hh1b->Fill(hit_z,hit_radius); //Histogram RecHit positions
                                                          }
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

			for( int iout = 0; iout < (int)out2d.size(); iout++){
				
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

				int stripMPGD = is2DStripMPGD(coll_name);

				//Find associated SimHit and see if hit should be used in track fit
				std::string should_be_meas = "N";
				auto raw_hit = out_hit.getRawHit();

				for (const auto raw_hit_assoc : raw_hit_assocs) {
					if (raw_hit_assoc.getRawHit() == raw_hit) {
						auto sim_hit = raw_hit_assoc.getSimHit();
                        			auto mc_particle = sim_hit.getParticle();
						auto mc_status = sim_hit.getParticle().getGeneratorStatus();
						auto quality = sim_hit.getQuality();

						if(mc_status==1 && quality==0){
							should_be_meas = "Y";
							num_out_wrong++;

							//Get index, adding (1||2)/3 in 2DStrip MPGD case.
							int pn; // p|n of 2DStrip MPGD
							float index = getIndex(coll_name,cellid,pn);
							if (pn == 3) continue; // It's a MPGD w/ no valid strip bit...
							//Fill associated outlier hit collection histogram
							if(index > 0){
                                                          h8c->Fill( index );
                                                          if (stripMPGD<0) { //Histogram RecHit positions
                                                            hh1c->Fill(hit_z,hit_radius);

                                                            if(coll_name=="TOFEndcapRecHits" && hit_z>1860.)
                                                              hh2c->Fill(hit_x,hit_y);
                                                          }
							}

							if (stripMPGD>=0) { // 2DStrip MPGDs: do we have 'p&n'?
                                                          CellID simID = sim_hit.getCellID();
                                                          SimID2RecIDs:: iterator im =
                                                            pANDn(simID,cellid,simID2RecIDs);
                                                          if (im != simID2RecIDs.end()) {
                                                            float pnIndex = int(index) + 2/3. + .1;
                                                            h8c->Fill( pnIndex ); //Histogram 'p&n' outliers
                                                            auto& rechit_coll = f.get<edm4eic::TrackerHitCollection>(coll_name);
                                                            if (get2DStripRZ(coll_name,rechit_coll,im,hit_radius,hit_z))
                                                              hh1c->Fill(hit_z,hit_radius); //Histogram RecHit positions
                                                          }
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
	hh1d->SetMinimum(0);

	hh2d->Add(hh2a,hh2b,1,-1); //Subtract track measurement hits from associated hits
        hh2d->Add(hh2c,-1); //Also subtract outlier hits
	hh2d->SetMinimum(0);

	//Histogram stacking
	THStack *stack1 = new THStack("stack1", "Stacked Histograms");
    	stack1->Add(h8b);
    	stack1->Add(h8c);

	THStack *stack2 = new THStack("stack2", "Events with a reconstructed track; ; Fraction of associated digitized hits");
        stack2->Add(h8bs);
        stack2->Add(h8cs);

	//Make plots
	TCanvas *c1 = new TCanvas("c1","c1");
	h1->Draw();

	TLatex *tex1 = new TLatex();
	tex1->SetTextSize(0.03);
	//tex1->DrawLatexNDC(0.2,0.7,"Events with a reconstructed truth-seeded track"); //Uncomment if using truth-seeded tracks
	tex1->DrawLatexNDC(0.2,0.7,"Events with a reconstructed real-seeded track");

	TCanvas *c2 = new TCanvas("c2","c2");
	h2->Draw();

	TCanvas *c3 = new TCanvas("c3","c3");
	h3->Draw();

	TCanvas *c4 = new TCanvas("c4","c4");
	h4->Draw();

	TCanvas *c5 = new TCanvas("c5","c5");
	h5->Draw();

	TCanvas *c6 = new TCanvas("c6","c6");
	h6->Draw();

	TCanvas *c7 = new TCanvas("c7","c7");
	h7->Draw();

	TCanvas *c8 = new TCanvas("c8","c8"); c8->SetBottomMargin(.15);
        h8a->Draw(); 
	std::function<void(TH1*)> resetMaxAtMPGD = [](TH1 *h) {
          // Rescale if maximum falls on 2DStrip MPGDs, lest histo collide w/ legend.
          double yMx; int bjn, bMx; for (bjn=bMx = 1, yMx =0; bjn<=h->GetNbinsX();
                                         bjn++) {
            double y = h->GetBinContent(bjn); if (y>yMx) { yMx = y; bMx = bjn; }
          }
          const char *nameMx = h->GetXaxis()->GetBinLabel(bMx);
          if (strstr(nameMx,"MPGD") || !strcmp(nameMx,"-n")) h->SetMaximum(yMx*1.3);
	};
	resetMaxAtMPGD(h8a);
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

	TCanvas *c8s = new TCanvas("c8s","c8s"); c8s->SetBottomMargin(.15);
	c8s->SetGrid();
	stack2->Draw("HIST");

	TCanvas *c9 = new TCanvas("c9","c9"); c9->SetBottomMargin(.15);
	h9a->Draw();
	resetMaxAtMPGD(h9a);
	h9b->Draw("same");

	TLegend *legend2 = new TLegend(0.4, 0.75, 0.9, 0.9);
        legend2->AddEntry(h9a,"Associated Digitized Hits (RecHits)","l");
        legend2->AddEntry(h9b,"Associated Measurement2D Hits", "f");
        legend2->SetBorderSize(0);
        legend2->SetFillStyle(0);
        legend2->SetTextSize(0.035);
        legend2->Draw();

	TCanvas *c10a = new TCanvas("c10a","c10a");
	hh1a->Draw();
	TCanvas *c10b = new TCanvas("c10b","c10b");
	hh1b->Draw();
	TCanvas *c10c = new TCanvas("c10c","c10c");
	hh1c->Draw();
	TCanvas *c10d = new TCanvas("c10d","c10d");
	hh1d->Draw();

	TCanvas *c11a, *c11b, *c11c, *c11d;
	if(ftof_back_xy){
          c11a = new TCanvas("c11a","c11a");
          hh2a->Draw();
          c11b = new TCanvas("c11b","c11b");
          hh2b->Draw();
          c11c = new TCanvas("c11c","c11c");
          hh2c->Draw();
          c11d = new TCanvas("c11d","c11d");
          hh2d->Draw();
	}

	//Print plots to files
	TCanvas *cs[] = {c1,c2,c3,c4,c5,c6,c7,c8,c8s,c9,
          /* */    c10a,c10b,c10c,c10d};
	int nCanvas = sizeof(cs)/sizeof(TCanvas*);
	TCanvas *c11s[] = {c11a,c11b,c11c,c11d};
	int nC11s = sizeof(c11s)/sizeof(TCanvas*);
	TCanvas *lastC;
	if (ftof_back_xy){ lastC = c11s[nC11s-1]; nC11s -= 1; }
	else             { lastC = cs[nCanvas-1]; nCanvas -= 1; }
	//PDF file
	c1->Print("plots/hit_matching.pdf(");
	for (int ic = 1; ic<nCanvas; ic++) cs[ic]->Print("plots/hit_matching.pdf");
	if (ftof_back_xy){
          for (int ic = 1; ic<nC11s; ic++) c11s[ic]->Print("plots/hit_matching.pdf");
	}
	lastC->Print("plots/hit_matching.pdf)");
	//ROOT file
	TFile *fout = TFile::Open("plots/hit_matching.TCanvas.root","RECREATE");
	for (int ic = 1; ic<nCanvas; ic++) cs[ic]->Write();
	if (ftof_back_xy){
          for (int ic = 1; ic<nC11s; ic++) c11s[ic]->Write();
	}
	lastC->Write();
	fout->Close();
}
// **********************************************************************
// ********************     2DStrip MPGDs     ***************************
// **********************************************************************
bool get2DStripRZ(const string coll_name, const edm4eic::TrackerHitCollection& rechit_coll, SimID2RecIDs::iterator im,
		  double &hit_radius, float &hit_z)
{
  // Upon input arg. "im" pointing to two cellIDs (expected to be associated 'p'
  // and 'n' RecHits) return the signed radius and Z coordinate of their
  // combination.
  vector<CellID> &recIDs = im->second;
  //Get the (x,y,z) of the corresponding two hits
  double gposPN[2][3];
  auto num_rechits = rechit_coll.size();
  unsigned int ihit, match; for(ihit=match = 0; ihit<num_rechits; ihit++){
    auto hit = rechit_coll.at(ihit); CellID hit_cellid = hit.getCellID();
    for (int j = 0; j<2; j++) {
      if (hit_cellid == recIDs[j]) {
	match |= 0x1<<j;
	gposPN[j][0] = hit.getPosition().x;
	gposPN[j][1] = hit.getPosition().y;
	gposPN[j][2] = hit.getPosition().z;
      }
    }
  }
  if (match != 0x3) { // Inconsistency: error message and return false
    printf("** get2DStripRZ: Error retrieving RecHit of \"%s\" w/ cellID:/n",
	   coll_name.c_str());
    for (int j = 0; j < 2; j++) {
      if (match&0x1<<j) {
	CellID cID = recIDs[j]; printf(" (0x%08lx,0x%08lx)",cID&0xffffffff,cID>>32);
      }
    }
    return false;
  }
  //Get R and Z out (x,y,z)'s
  int stripMPGD = is2DStripMPGD(coll_name);
  if (stripMPGD == 1) {      // ***** Outer Barrel
    // Outer Barrel segmentation is expected to be CartesianGridUV. =>
    // RecHits coordinates expressed in World Reference System are not affected
    // by large uncertainties.
    // => We have to:
    // - Transform to LocalRS.
    // - There, combine U of 'p' strips w/ V of 'n' to obtain Local (X,Y,Z).
    // - Transform back to WorldRS.
    Geometry *geometry = Geometry::Ptr();
    if (!geometry) {
      printf("** get2DStripRZ: No access to geometry\n");
      return false;
    }
    double lpos[3], Ur=0, Vr=0;
    for (int j = 0; j < 2; j++) {
      CellID recID = recIDs[j]; int pn = (recID&pStripBit) ? 0 : 1;
      double *gpos = gposPN[j]; 
      if (geometry->WorldToLocal(1,recID,gpos,lpos)) {
	if (pn==0) Ur = lpos[1]-lpos[0];
	else       Vr = lpos[1]+lpos[0];
      }
      else return false;
    }
    // p&n Local -> World
    lpos[0] = (Vr-Ur)/2; lpos[1] = (Vr+Ur)/2;
    double gpos[3];  if (!geometry->LocalToWorld(1,recIDs[0],lpos,gpos))
		       return false;
    hit_radius = sqrt(gpos[0]*gpos[0]+gpos[1]*gpos[1])*sgn<double>(gpos[0]);
    hit_z = gpos[2];
  }
  else if (stripMPGD == 0) { // ***** CyMBaL
    for (int j = 0; j < 2; j++) {
      CellID recID = recIDs[j]; int pn = (recID&pStripBit) ? 0 : 1;
      double *gpos = gposPN[j];
      if (pn == 0) { // phi strip
	double &X = gpos[0], &Y = gpos[1];
	hit_radius = sgn<double>(X) * sqrt(X*X + Y*Y);
      }
      else           // z strip
	hit_z = gpos[2];
    }
  }
  else {                     // ***** ECT
    double X=0, Y=0;
    for (int j = 0; j < 2; j++) {
      CellID recID = recIDs[j]; int pn = (recID&pStripBit) ? 0 : 1;
      double *gpos = gposPN[j];
      if (pn == 0) { // x strip
	X = gpos[0];
	hit_z = gpos[2]; // Z can be taken indifferently from x or y strip
      }
      else           // z strip
	Y = gpos[1];
    }
    hit_radius = sgn<double>(X) * sqrt(X*X + Y*Y);
  }
  return true;
}
// ********** GEOMETRY
Geometry::Geometry(int verbose)
{
  // ***** CONSTRUCTOR
  // - Explore the TGeometry to:
  //  i) Determine which of the MPGDs are 2DStrip (in fact what's what's rather
  //    determined is whether they have MultipleSensitiveVolume, which as of
  //    2026/04, is the only way MPGDs can be 2DStrip).
  // => Stored in data member "MPGD2DStrips"
  // ii) Retrieve the World-To-Local transforms.
  // => Stored in "geoDetMatrix[detector][layer][module]".
  // - Assumptions (some effort is made to limit their number):
  //   - All MPGD modules can be reached via the TGeometry tree using a minimal
  //    set of tags (see "findLayersModules").
  //     Layer's node name includes "(L|l)ayer". Module's include "(M|m)odule".
  //     Sensitive surface is on a node named...
  //    ..."ReferenceThinGap" for 2DStrip (in fact MSV) MPGDs,
  //    ..."DriftGap" otherwise.
  //   - Encoding of layer# and module# numbering in IDDescriptor starts @ 0 (
  //    as opposed to @ 1).
  //   - Possibly others...
  
  // Init
  MPGD2DStrips = 0;
  if(address == 0) //!< Protection against multiple instances
    address = this;
  verboseLevel = verbose;

  // Get the TGeometry.
  // - Let's first check that the file exists in PWD, in order to avoid the
  //  possible warning message printed out by "TGeoManager::Import" when the
  //  file is missing.
  const char geoFN[] = "hit_matching.geometry.root";
  TSystemDirectory *sys = new TSystemDirectory(".",".");
  if (!sys->GetListOfFiles()->FindObject(geoFN)) {
    printf("No \"%s\" file in $PWD. => All MPGDs to be processed as pixels\n",
	   geoFN);
    return;
  }
  // - Import the TGeometry
  gGeoManager->SetVerboseLevel(verbose);
  if (!gGeoManager->Import(geoFN)) {
    printf("Error Importing \"%s\" into gGeoManager. => All MPGDs to be processed as pixels\n",
	   geoFN);
    return;
  }

  // Set "MPGD2DStrips" and "geoDetMats".
  for (int mpgd = 0; mpgd<N_MPGDs; mpgd++) {
    bool isMSV = false;
    if (getGeoMats(mpgd,isMSV)) {
      if (isMSV) MPGD2DStrips |= 0x1<<mpgd;
    }
  }
};
bool findLayersModules(int LM,       // 0: Layer, 1: Module 
		       string &path, // Initial path
		       const char *tag1, const char *tag2,
		       int verbose,
		       vector<string>& targetNodes);
bool findDriftGap(const char *modulePath, int verbose,
		  const TGeoHMatrix *&m, bool &isMSV);
bool Geometry::getGeoMats(int mpgd, bool &isFullyMSV)
{
  // 

  isFullyMSV = true;

  vector<vector<TGeoHMatrix>> &geoMats = geoDetMats[mpgd];

  const char *tag1 = tag1s[mpgd], *tag2 = tag2s[mpgd];
  vector<string> layerNodes;
  string path("/"); gGeoManager->CdTop();
  path += string(gGeoManager->GetCurrentNode()->GetName());
  bool ok = findLayersModules(0,path,tag1,tag2,verboseLevel,layerNodes);
  if (ok) {
    int nLayers = layerNodes.size();
    if (verboseLevel>1)
      printf("\"%s%c%s\": Found %d Layer(s): \n",
	     tag1,tag2?',':'\0',tag2?tag2:"\0",nLayers);
    for (int layer = 0; layer<nLayers; layer++) {
      const char *layerPath = layerNodes[layer].c_str();
      if (verboseLevel>1)
	printf("%2d/%2d: \"%s\"\n",layer,nLayers,layerPath);
      vector<string> moduleNodes;
      path = layerNodes[layer]; gGeoManager->cd(path.c_str());
      ok = findLayersModules(1,path,tag1,tag2,verboseLevel,moduleNodes);
      if (!ok) break;
      int nModules = moduleNodes.size();
      if (verboseLevel>1)
	printf("\"%s%c%s\", Layer %d: Found %d Modules: \n",
	       tag1,tag2?',':'\0',tag2?tag2:"\0",layer,nModules);
      vector<TGeoHMatrix> gMats;
      for (int module = 0; module<nModules; module++) {
	const char *modulePath = moduleNodes[module].c_str();
	if (verboseLevel>1)
	  printf("%2d/%2d: \"%s\"\n",module,nModules,modulePath);
	const TGeoHMatrix *m = 0; bool isMSV = false;
	if (findDriftGap(modulePath,verboseLevel,m,isMSV)) {
	  // cm -> mm
	  TGeoHMatrix mmm(*m);
	  const double *trcm = mmm.GetTranslation(); double trmm[3];
	  for (int xyz = 0; xyz<3; xyz++) trmm[xyz] = trcm[xyz]*10;
	  mmm.SetTranslation(trmm);
	  gMats.push_back(mmm);
	  isFullyMSV &= isMSV;
	}
	else {
	  if (verboseLevel) 
	    printf("** getGeoMats: Could not get TGeoHMatrix of Drift Gap in node \"%s\"\n",
		   modulePath);
	}
      }
      geoMats.push_back(gMats);
    }
  }
  if (!ok) {
    printf("** getGeoMats: No Module found for \"%s%c%s\"\n",
	   tag1,tag2?',':'\0',tag2?tag2:"\0");
    isFullyMSV = false;
  }
  if (verboseLevel) {
    int nLayers = geoMats.size();
    if (verboseLevel)
      printf("\"%s%c%s\": Found %d Layers: \n",
	     tag1,tag2?',':'\0',tag2?tag2:"\0",nLayers);
    for (int layer = 0; layer<nLayers; layer++) {
      vector<TGeoHMatrix>& gMats = geoMats[layer];
      int nModules = gMats.size();
      if (verboseLevel)
	printf("\"%s%c%s\", Layer %d: Found %d Modules: \n",
	       tag1,tag2?',':'\0',tag2?tag2:"\0",layer,nModules);
      for (int module = 0; module<nModules; module++) {
	TGeoHMatrix& geoMat = gMats[module];
	const double *tr = geoMat.GetTranslation();
	if (verboseLevel>1)
	  printf("%2d/%d: %8.4f,%8.4f,%8.4f cm\n",
		 module,nModules,tr[0],tr[1],tr[2]);
      }
    }
  }
  return true;
}
bool findLayersModules(int LM, string &path,
		       const char *tag1, const char *tag2,
		       int verbose,
		       vector<string>& targetNodes)
{
  // Find all TGeoNodes w/ name including:
  //  i) input tags ("tag2" can be =0, it's needed for Endcaps which share a
  //    common "assembly" node),
  //  AND if LM=0, the word "MPGD" (it's needed for the earlier levels, to
  //    raise ambiguities),
  // ii) LM=0: the word "Layer" (or "layer").
  //     LM=1: the word "Module" (or "module").
  // Collect their path in "targetNodes".
  // - Starting point is input "path"
  // - Then node fulfilling (i) is recursively searched...
  // - ...until nodes fulfilling (ii) are found (in the mean time, only one
  //  node is allowed, else we wouldn't know which path to then follow).
  // This strategy is presumed to work in all kinds of circumstance, possibly
  // differing from the situation as of 2026/04, where we have, e.g.:
  //  "/world_volume_1
  //    /OuterBarrelMPGDSubAssembly_10
  //      /MPGDOuterBarrel_0
  //        /MPGDOuterBarrel_layer0_0
  //          /MPGDOuterBarrelModule_0,2,...,23"
  // A maximum depth of 10 is still imposed in addition.
  const TGeoNode *node = gGeoManager->GetCurrentNode(), *found;
  for (int depth = 0; depth<10; depth++) {
    const TObjArray* subNodeArray = node->GetNodes();
    TObject *o = subNodeArray->First();
    found = 0; while (o) {
      if (!o->IsA()->InheritsFrom(TGeoNode::Class())) {
	printf("** findLayersModules: Object \"%s\" of Class \"%s\" among daughter of TGeoNode \"%s\"\n",
	       o->GetName(),o->ClassName(),node->GetName());
	return false;
      }
      const TGeoNode *subNode = (TGeoNode*)o;
      const char *name = subNode->GetName();
      if (verbose>2) {
	const TGeoMatrix *m = subNode->GetMatrix();
	const double *tr = m->GetTranslation();
	printf("  \"%s\": %8.4f,%8.4f,%8.4f cm\n",name,tr[0],tr[1],tr[2]);
      }
      // Apply requirement (i)
      if ((LM==1 || strstr(name,"MPGD")) &&
	  (strstr(name,tag1) || tag2 && strstr(name,tag2))) {
	bool isTarget = LM==0 ? strstr(name,"Layer")  || strstr(name,"layer")
	  /* */               : strstr(name,"Module") || strstr(name,"module");
	if (!isTarget && found) {
	  printf("** findLayersModules: Ambiguity: TGeoNode \"%s\"(=> not a %s) has more than one subNode: \"%s\" and \"%s\" ",
		 node->GetName(),LM==0?"Layer":"Module",found->GetName(),name);
	  return false;
	}
	if (isTarget) // Requirement (ii): Is Add path to the list
	  targetNodes.push_back(path+string("/")+string(subNode->GetName()));
	found = subNode;
      }
      o = subNodeArray->After(o);
    }
    if (!found) return false;
    if (targetNodes.empty()) { // Increment "path" and keep searching
      path += string("/") + string(found->GetName());
      node = found;
      if (verbose>1) {
	const TGeoMatrix *m = found->GetMatrix();
	const double *tr = m->GetTranslation();
	printf("\"%s\": %8.4f,%8.4f,%8.4f cm\n",path.c_str(),tr[0],tr[1],tr[2]);
      }
    }
    else
      break;
  }
  return !targetNodes.empty();
}
bool findDriftGap(const char *modulePath, int verbose,
		  const TGeoHMatrix *&m, bool &isMSV)
{
  // Return TGeoHMatrix of World-to-Local of Reference SubVolume, which name
  // is assumed to be "ReferenceThinGap".
  if (!gGeoManager->cd(modulePath)) {
    printf("** findDriftGap: Bad input path \"%s\"\n",
	   modulePath);
    return false;
  }
  TGeoNode *node = gGeoManager->GetCurrentNode();
  const TObjArray* subNodeArray = node->GetNodes();
  TObject *o = subNodeArray->First();
  const TGeoNode* found = 0; while (o) {
    if (!o->IsA()->InheritsFrom(TGeoNode::Class())) {
      printf("** findDriftGap: Object \"%s\" of Class \"%s\" among daughter of TGeoNode \"%s\"\n",
	     o->GetName(),o->ClassName(),node->GetName());
      return false;
    }
    const TGeoNode *subNode = (TGeoNode*)o;
    const char *name = subNode->GetName();
    if (verbose>2) {
      const TGeoMatrix *m = subNode->GetMatrix();
      const double *tr = m->GetTranslation();
      printf("  \"%s\": %8.4f,%8.4f,%8.4f cm\n",name,tr[0],tr[1],tr[2]);
    }
    isMSV = strstr(name,"ReferenceThinGap");
    if (isMSV || strstr(name,"DriftGap")) {
      found = subNode; break;
    }
    o = subNodeArray->After(o);
  }
  if (!found) return false;
  string nodePath = string(modulePath)+string("/")+string(found->GetName());
  if (!gGeoManager->cd(nodePath.c_str())) {
    printf("** findDriftGap: Inconsistency: gGeoManager cannot cd(\"%s\")\n",
	   nodePath.c_str());
    return false;
  }
  m = gGeoManager->GetCurrentMatrix();
  if (verbose>1) {
    const double *tr = m->GetTranslation();
    printf("\"%s\": %8.4f,%8.4f,%8.4f cm\n",found->GetName(),tr[0],tr[1],tr[2]);
  }
  return true;
}
void cellID2LayerModule(int idet, unsigned long cellID,
			unsigned int &layer, unsigned int &module)
{
  // Assuming IDDescriptor to be:
  if (idet<2) { // Barrel:
    // <id>system:8,layer:4,module:12,         strip:28:4,u:-16,v:-16</id>
    layer =  cellID>>8&0xf, module = cellID>>12&0xfff;
  }
  else {        // ECT
    // <id>system:8,layer:2,module:6,sensor:16,strip:28:4,x:32:-16,z:-16</id>
    layer =  cellID>>8&0x3, module = cellID>>10&0x3f;
  }
}
bool Geometry::WorldToLocal(int idet, unsigned long cellID,
			    double *gpos, // In mm
			    double *lpos) // In mm
{
  unsigned int layer, module; cellID2LayerModule(idet,cellID,layer,module);
  bool ok = false; 
  if (idet<N_MPGDs) {
    vector<vector<TGeoHMatrix>> &geoMats = geoDetMats[idet];
    size_t nLayers = geoMats.size(), nModules = 0; if (layer<nLayers) {
      if (module<geoMats[layer].size()) ok = true;
    }
    if (ok) {
      const TGeoHMatrix &geoMat = geoMats[layer][module];
      geoMat.MasterToLocal(gpos,lpos);
    }
  }
  if (!ok) {
    printf("** Geometry::WorldToLocal: No Geometry found for mpgd=%d, cellID = 0x%08lx => (layer,module) = (%u,%u)\n",idet,cellID&0xffffffff,layer,module);
    return false;
  }
  return true;
}
bool Geometry::LocalToWorld(int idet, unsigned long cellID,
			    double *lpos, // in mm
			    double *gpos) // in mm
{
  unsigned int layer, module; cellID2LayerModule(idet,cellID,layer,module);
  bool ok = false; 
  if (idet<N_MPGDs) {
    vector<vector<TGeoHMatrix>> &geoMats = geoDetMats[idet];
    size_t nLayers = geoMats.size(), nModules = 0; if (layer<nLayers) {
      if (module<geoMats[layer].size()) ok = true;
    }
    if (ok) {
      const TGeoHMatrix &geoMat = geoMats[layer][module];
      geoMat.LocalToMaster(lpos,gpos);
    }
  }
  if (!ok) {
    printf("** Geometry::LocalToWorld: No Geometry found for mpgd=%d, cellID = 0x%08lx => (layer,module) = (%u,%u)\n",idet,cellID&0xffffffff,layer,module);
    return false;
  }
  return true;
}
