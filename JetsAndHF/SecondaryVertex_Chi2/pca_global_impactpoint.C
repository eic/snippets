// Example code provide by Barak Schmookler
// for calculating DCA to a specific point
// Important first execute this command: source /opt/detector/epic-main/bin/thisepic.sh 
// then do root -l pca_global_impactpoint.C

#include "DD4hep/Detector.h"
#include "DD4hep/DetElement.h"
#include "ActsPlugins/DD4hep/ConvertDD4hepDetector.hpp"
#include "ActsPlugins/DD4hep/DD4hepFieldAdapter.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "ActsPlugins/DD4hep/ConvertDD4hepDetector.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"


void pca_global_impactpoint(TString filelist = "test.list"){

     // Read all files
   bool debug = false;
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

    TTreeReader tr(chain);


    // Load DD4Hep geometry
    dd4hep::Detector& detector = dd4hep::Detector::getInstance();
    detector.fromCompact("/opt/detector/epic-main/share/epic/epic_craterlake.xml");
    dd4hep::DetElement geometry = detector.world();

    // Convert DD4Hep geometry to tracking geometry
    Acts::GeometryContext trackingGeoCtx;
    auto logger = Acts::getDefaultLogger("DD4hepConversion", Acts::Logging::Level::INFO);
    Acts::BinningType bTypePhi = Acts::equidistant;
    Acts::BinningType bTypeR = Acts::equidistant;
    Acts::BinningType bTypeZ = Acts::equidistant;
    double layerEnvelopeR = Acts::UnitConstants::mm;
    double layerEnvelopeZ = Acts::UnitConstants::mm;
    double defaultLayerThickness = Acts::UnitConstants::fm;
    using ActsPlugins::sortDetElementsByID;

    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry{nullptr};
    trackingGeometry = ActsPlugins::convertDD4hepDetector(geometry,*logger,bTypePhi,bTypeR,bTypeZ,layerEnvelopeR,layerEnvelopeZ,defaultLayerThickness,sortDetElementsByID,trackingGeoCtx);

    // Define Perigee surface at which reconstructed track parameters are set
    auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(Acts::Vector3(0,0,0));

    // Get Magnetic field context
    Acts::MagneticFieldContext fieldctx;
    std::shared_ptr<const ActsPlugins::DD4hepFieldAdapter> field_provider = std::make_shared<const ActsPlugins::DD4hepFieldAdapter>(detector.field());
    Acts::MagneticFieldProvider::Cache field_cache = field_provider->makeCache(fieldctx);

	// Stepper and Propagator
    using Stepper    = Acts::EigenStepper<>;
    using Propagator = Acts::Propagator<Stepper>;

    Stepper stepper(field_provider);
    Propagator propagator(stepper);

	// Create Impact Point Estimator
    Acts::ImpactPointEstimator::Config ImPoEs_cfg(field_provider,std::make_shared<Propagator>(propagator));

    Acts::ImpactPointEstimator::State ImPoEs_state;
    ImPoEs_state.fieldCache = field_cache;

    Acts::ImpactPointEstimator ImPoEs(ImPoEs_cfg);

    // Create 'vertex' at particle's creation point -- which is (x,y,z) = (1,0,0) mm
    Acts::Vector3 vtx_pos(1.0 * Acts::UnitConstants::mm, 0, 0);

    // Create another 'vertex' at (2,0,0) mm and check the distance at the DCA to this point
    Acts::Vector3 vtx_pos2(2.0 * Acts::UnitConstants::mm, 0, 0);

    // Define Style
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

    // Define histograms
    TH2 *h1 = new TH2D("h1"," particle generated at (x,y,z) = (1,0,0) mm",100,-1.5,1.5,100,-1.5,1.5);
    h1->GetXaxis()->SetTitle("Track global x at beamline (z-axis) POCA [mm]");h1->GetXaxis()->CenterTitle();
    h1->GetYaxis()->SetTitle("Track global y at beamline (z-axis) POCA [mm]");h1->GetYaxis()->CenterTitle();

    TH1 *h2 = new TH1D("h2"," particle generated at (x,y,z) = (1,0,0) mm",100,-0.02,1);
    h2->GetXaxis()->SetTitle("Track distance to generation point at 3D DCA [mm]");h2->GetXaxis()->CenterTitle();
    h2->SetLineWidth(2);h2->SetLineColor(kBlue);

    TH1 *h3 = new TH1D("h3"," particle generated at (x,y,z) = (1,0,0) mm",100,0,2);
    h3->GetXaxis()->SetTitle("Track distance to (x,y,z) = (2,0,0) mm at 3D DCA [mm]");h3->GetXaxis()->CenterTitle();
    h3->SetLineWidth(2);h3->SetLineColor(kBlue);

    TH1 *h3a = new TH1D("h3a"," particle generated at (x,y,z) = (1,0,0) mm",100,0,2);
    h3a->GetXaxis()->SetTitle("Track xy distance to (x,y,z) = (2,0,0) mm at 3D DCA [mm]");h3a->GetXaxis()->CenterTitle();
    h3a->SetLineWidth(2);h3a->SetLineColor(kBlue);

    TH1 *h3b = new TH1D("h3b"," particle generated at (x,y,z) = (1,0,0) mm",100,-1,1);
    h3b->GetXaxis()->SetTitle("Track z signed distance to (x,y,z) = (2,0,0) mm at 3D DCA [mm]");h3b->GetXaxis()->CenterTitle();
    h3b->SetLineWidth(2);h3b->SetLineColor(kBlue);

    TH1 *h4a = new TH2D("h4a"," particle generated at (x,y,z) = (1,0,0) mm",100,-3.2,3.2,100,0,2);
    h4a->GetXaxis()->SetTitle("Generated particle #phi [Rad]");h4a->GetXaxis()->CenterTitle();
    h4a->GetYaxis()->SetTitle("Track xy distance to (x,y,z) = (2,0,0) mm at 3D DCA [mm]");h4a->GetYaxis()->CenterTitle();

    TH1 *h4b = new TH2D("h4b"," particle generated at (x,y,z) = (1,0,0) mm",100,0,3.2,100,-1,1);
    h4b->GetXaxis()->SetTitle("Generated particle #theta [Rad]");h4b->GetXaxis()->CenterTitle();
    h4b->GetYaxis()->SetTitle("Track z signed distance to (x,y,z) = (2,0,0) mm at 3D DCA [mm]");h4b->GetYaxis()->CenterTitle();

    // Load input ROOT file
    int pid_code = 13;
    std::string coll = "CentralCKFTrackParameters"; //Real-seeded track parameters
    TLorentzVector gen_vec; //Generated particle four momentum

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

    TTreeReaderArray<float> track_qoverp(tr, Form("%s.qOverP",coll.c_str()));
    TTreeReaderArray<float> track_theta(tr, Form("%s.theta",coll.c_str()));
    TTreeReaderArray<float> track_phi(tr, Form("%s.phi",coll.c_str()));
    TTreeReaderArray<float> track_loca(tr, Form("%s.loc.a",coll.c_str()));
    TTreeReaderArray<float> track_locb(tr, Form("%s.loc.b",coll.c_str()));

    //Loop over events
    int counter(0);
    while (tr.Next()) {

        cout<<"Analyzing event "<<counter<<endl;
        

        //Loop over generated particles, select primary particle (assuming  particle)
        for(size_t igen=0;igen<gen_status.GetSize();igen++){

            if(gen_status[igen]==1 && gen_pid[igen]==pid_code){ //PID code requirement not really needed...
                gen_vec.SetXYZM(gen_px[igen],gen_py[igen],gen_pz[igen],gen_mass[igen]);
                break;
            }
        } //End loop over generated particles

        //Loop over tracks
        size_t track_mult = track_qoverp.GetSize();
        for(size_t itrack=0;itrack<track_mult;itrack++){
         printf("Track No: %zu \n",itrack);

         auto loc_a = track_loca[itrack];
         auto loc_b = track_locb[itrack];
         auto phi = track_phi[itrack];
         auto theta = track_theta[itrack];
         auto qoverP = track_qoverp[itrack];

            // Create BoundTrackParamters
         Acts::BoundVector params;

         params(Acts::eBoundLoc0)   = loc_a;
         params(Acts::eBoundLoc1)   = loc_b;
         params(Acts::eBoundPhi)    = phi;
         params(Acts::eBoundTheta)  = theta;
         params(Acts::eBoundQOverP) = qoverP;
         params(Acts::eBoundTime)   = 0; 

            //FIXME: Set covariance matrix based on input ROOT file information
         Acts::BoundSquareMatrix cov = Acts::BoundSquareMatrix::Zero();

         Acts::BoundTrackParameters track_parameters(perigee,params,cov,Acts::ParticleHypothesis::pion());

            //---- Part 1: Convert from local coordinates to global coordinates at beamline (z-axis) POCA ----
         Acts::Vector2 localpos( loc_a, loc_b );
         Acts::Vector3 direction(sin(theta)*cos(phi), 
            sin(theta)*sin(phi), 
            cos(theta));
         if (debug) printf("Local Pos (%f,%f) \n",loc_a,loc_b);                         

            // Convert to global coordinates for PCA at perigee surface
         auto global_perigee = perigee->localToGlobal(trackingGeoCtx,localpos,direction);

            /*std::cout << "Tracking Geometry Context Type: " << typeid(trackingGeoCtx).name() << std::endl;
            // Correct way to check surface type perigee surface never distinguishes between plane or cylinder
            // Get surface name https://github.com/acts-project/acts/blob/bbd3fcc56f1c059c72d0d1caf040e5f5fa596edb/Core/include/Acts/Surfaces/Surface.hpp#L62C1-L75C5
            int perigee_type = perigee->type();
            std::string surfaceTypeStr = perigee->name(); 
            std::cout<<"Perigee Type: "<<perigee_type<<" Perigee Name:  "<<surfaceTypeStr<<std::endl;
            
            const auto& bounds = perigee->bounds();
            std::cout << "Surface Bounds Type ID: " << bounds.type() <<typeid(bounds).name()<< std::endl;*/
         Acts::Transform3 globalTransform = perigee->transform(trackingGeoCtx);
         std::cout << "Global Transformation Matrix:\n" << globalTransform.matrix() << std::endl;


         if( !global_perigee.isZero() ){
            h1->Fill(global_perigee.x(),global_perigee.y());
        }
        TVector3 pos(-1*loc_a* sin(phi), loc_a*cos(phi),loc_b);
        
        double eta = pos.Eta();
            if (debug && fabs(eta)>3){ // apply eta cut
                printf("Global Pos ACTS: (%f,%f,%f) \n",global_perigee.x(),global_perigee.y(),global_perigee.z());  
                printf("Global Pos Helix swimming: (%f,%f,%f) \n",pos.x(),pos.y(),pos.z());  
            }                       

            //---- Part 2: Get track parameters at 3D DCA to creation point at (x,y,z) = (1,0,0) mm ----
            auto result = ImPoEs.estimate3DImpactParameters(trackingGeoCtx,fieldctx,track_parameters,vtx_pos,ImPoEs_state);

            if(result.ok()){
              Acts::BoundTrackParameters trk_boundpar_vtx = result.value();
              const auto& trk_vtx_params  = trk_boundpar_vtx.parameters();

              h2->Fill(trk_vtx_params[Acts::eBoundLoc0]);
          }

            //---- Part 2a: Get track parameters at 3D DCA to (x,y,z) = (2,0,0) mm ----
          auto result2 = ImPoEs.estimate3DImpactParameters(trackingGeoCtx,fieldctx,track_parameters,vtx_pos2,ImPoEs_state);

          if(result2.ok()){
              Acts::BoundTrackParameters trk_boundpar_vtx2 = result2.value();
              const auto& trk_vtx_params2  = trk_boundpar_vtx2.parameters();

              h3->Fill(trk_vtx_params2[Acts::eBoundLoc0]);

                //Get global position at 3D DCA
              auto trk_vtx2_gbl_pos = trk_boundpar_vtx2.position(trackingGeoCtx);

              auto delta_x = trk_vtx2_gbl_pos.x() - 2.;
              auto delta_y = trk_vtx2_gbl_pos.y() - 0.;
              auto delta_z = trk_vtx2_gbl_pos.z() - 0.;

              auto xy_distance = std::hypot(delta_x,delta_y);

              h3a->Fill(xy_distance);
              h3b->Fill(delta_z);
              h4a->Fill(gen_vec.Phi(),xy_distance);
              h4b->Fill(gen_vec.Theta(),delta_z);
          }

        } // End loop over tracks
        counter++;
        
    } // End loop over events

    // Make plots
    TCanvas *c1 = new TCanvas("c1");
    h1->Draw();

    TCanvas *c2 = new TCanvas("c2");
    h2->Draw();

    TCanvas *c3 = new TCanvas("c3");
    h3->Draw();

    TCanvas *c3a = new TCanvas("c3a");
    h3a->Draw();

    TCanvas *c3b = new TCanvas("c3b");
    h3b->Draw();

    TCanvas *c4a = new TCanvas("c4a");
    h4a->Draw();

    TCanvas *c4b = new TCanvas("c4b");
    h4b->Draw();

    //Print plots to file
    c1->Print("pca_global_impactpoint.pdf[");
    c1->Print("pca_global_impactpoint.pdf");
    c2->Print("pca_global_impactpoint.pdf");
    c3->Print("pca_global_impactpoint.pdf");
    c3a->Print("pca_global_impactpoint.pdf");
    c3b->Print("pca_global_impactpoint.pdf");
    c4a->Print("pca_global_impactpoint.pdf");
    c4b->Print("pca_global_impactpoint.pdf");
    c3b->Print("pca_global_impactpoint.pdf]");
}
