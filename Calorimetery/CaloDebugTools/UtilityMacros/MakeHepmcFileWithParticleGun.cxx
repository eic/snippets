//////////////////////////////////////////////////////////////
// 'MakeHepmcFileWithParticleGun.cxx'
// Dhevan Gangadharan, Nico Schmidt
//
// Derived from emcal_barrel_particles_gen.cxx. Generates
// HepMC files for events of 1 or more particles.
//////////////////////////////////////////////////////////////
#include <cmath>
#include <string>
#include <iostream>
#include <math.h>
#include <random>
#include <cctype>
#include <fmt/core.h>

#include "TFile.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TF2.h"
#include "TH2D.h"
#include "TLorentzVector.h"

#include "HepMC3/GenEvent.h"
#include "HepMC3/Print.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"

using namespace std;
using namespace HepMC3;

std::tuple <int, double> extract_particle_parameters(std::string particle_name);

// converter center: -55609.5, full thickness = 1 mm
// end of spectrometer magnet: -56390
double Vz = 0; // Primary vertex location in mm

double Z = 1; // ion beam particle charge
double electronPz = -18;
double hadronPz = 275;

double prefactor = 2.3179; // 4 alpha r_e^2 (mb)
double protonMass = 0.938272;
double electronMass = 0.51099895e-3;
double photonMass = 0;
double muonMass = 0.1056583745;
double pionZeroMass = 0.1349768;
double pionMass = 0.13957039;

// output name [Derek, 02.22.2024]
const string sOutBase = "forBHCalMultiPartSim.e10h11pipXpim.d8m3y2024";

void MakeHepmcFileWithParticleGun(
  int n_events = 1e4,
  double E_start = 10.0,
  double E_end = 10.0,
  double eta_min = -1.1,
  double eta_max = 1.1,
  string out_fname= sOutBase + ".hepmc"
) {

  // make output root file name [Derek, 02.22.2024]
  const string sOutRootName = sOutBase + ".hepmc.root";
 
  TFile *fout = new TFile(sOutRootName.data(), "RECREATE");
  
  cout << "Generating " << n_events << " events" << endl;
  cout << "particle energy range: " << E_start << " - " << E_end << " GeV" << endl;
  cout << "particle eta range: " << eta_min << " - " << eta_max << endl;


  WriterAscii hepmc_output(out_fname);
  
  int events_parsed = 0;
  
  GenEvent evt(Units::GEV, Units::MM);

  // Random number generator
  TRandom* r1 = new TRandom3();
  
  // Create events
  for (events_parsed = 0; events_parsed < n_events; events_parsed++) {
    // FourVector( px, py, pz, e, pdgid, status )

    // Create a vertex for the beam particles
    GenVertexPtr v1 = std::make_shared<GenVertex>();
    
    // Generate beam particles
    // status 1 (final state particle), 
    // status 2 (decayed physical particle)
    // status 3 (Documentation line, can be used to indicate in/out particles of hard process)
    // status 4 (beam particle)
    GenParticlePtr p1 = std::make_shared<GenParticle>(
        FourVector(0.0, 0.0, electronPz, sqrt( pow(electronPz, 2) + pow(electronMass, 2) ) ), 11, 4);
    GenParticlePtr p2 = std::make_shared<GenParticle>(
        FourVector(0.0, 0.0, hadronPz, sqrt( pow(hadronPz, 2) + pow(protonMass, 2) ) ), 2212, 4);
    GenParticlePtr p2_2 = std::make_shared<GenParticle>(
        FourVector(0.0, 0.0, hadronPz, sqrt( pow(hadronPz, 2) + pow(protonMass, 2) ) ), 2212, 1);
    
    v1->add_particle_in( p1 );
    v1->add_particle_in( p2 );
    v1->add_particle_out( p2_2 ); // make the output particle a beam proton (not to be in analysis)

    /////////////////////////////////////////////////////////////////
    // Now add particles of interest at separate vertices
    // vertices must have at least 1 outgoing and 1 incoming particle

    double E = r1->Uniform(E_start, E_end);
    double p = sqrt( pow(E,2) - pow(pionMass,2) );
    double eta1 = r1->Uniform(eta_min, eta_max);
    double eta2 = r1->Uniform(eta_min, eta_max);
    double theta1 = 2*atan(exp(-eta1));
    double theta2 = 2*atan(exp(-eta2));
    // double theta1 = TMath::Pi();
    // double theta2 = TMath::Pi();
    double phi1 = r1->Uniform(0, 2*TMath::Pi());
    double phi2 = r1->Uniform(0, 2*TMath::Pi());
    auto [id, mass] = extract_particle_parameters( "piplus" );

    // pion 1 
    GenVertexPtr v2 = std::make_shared<GenVertex>();
    v2->set_position( FourVector(0.0, 0.0, 0.0, 0) ); // in mm
    // v2->set_position( FourVector(0.0, 120.0, -64000, 0) ); // in mm
    GenParticlePtr p3_in = std::make_shared<GenParticle>( 
        FourVector(p*sin(theta1)*cos(phi1), p*sin(theta1)*sin(phi1), p*cos(theta1), E), id, 3);
    GenParticlePtr p3_out = std::make_shared<GenParticle>( 
        FourVector(p*sin(theta1)*cos(phi1), p*sin(theta1)*sin(phi1), p*cos(theta1), E), id, 1);
 
    v2->add_particle_in( p3_in );
    v1->add_particle_out( p3_out );


    E = r1->Uniform(E_start, E_end);
    p = sqrt( pow(E,2) - pow(pionMass,2) );

    auto [id2, mass2] = extract_particle_parameters( "piminus" );
    // pion 2
    GenVertexPtr v3 = std::make_shared<GenVertex>();
    v3->set_position( FourVector(0.0, 0.0, 0.0, 0) ); // in mm
    GenParticlePtr p4_in = std::make_shared<GenParticle>( 
        FourVector(p*sin(theta2)*cos(phi2), p*sin(theta2)*sin(phi2), p*cos(theta2), E), id2, 3);
    GenParticlePtr p4_out = std::make_shared<GenParticle>( 
        FourVector(p*sin(theta2)*cos(phi2), p*sin(theta2)*sin(phi2), p*cos(theta2), E), id2, 1);

    v3->add_particle_in( p4_in );
    v1->add_particle_out( p4_out );

    evt.add_vertex( v1 );
    //evt.add_vertex( v2 );
    //evt.add_vertex( v3 );


    // Shift the whole event to specific point
    evt.shift_position_to( FourVector(0,0, Vz, 0) );

    if (events_parsed == 0) {
      std::cout << "First event: " << std::endl;
      Print::listing(evt);
    }

    hepmc_output.write_event(evt);
    
    if (events_parsed % 10000 == 0) {
      std::cout << "Event: " << events_parsed << std::endl;
    }
    
    evt.clear();
  }
  hepmc_output.close();
  std::cout << "Events parsed and written: " << events_parsed << std::endl;

  fout->Close();

}

//----------------------------------------------------------------------------
// Returns particle pdgID and mass in [GeV]
std::tuple <int, double> extract_particle_parameters(std::string particle_name) {
  if (particle_name == "electron") return std::make_tuple(11,    electronMass);
  if (particle_name == "photon")   return std::make_tuple(22,    photonMass);
  if (particle_name == "positron") return std::make_tuple(-11,   electronMass);
  if (particle_name == "proton")   return std::make_tuple(2212,  protonMass);
  if (particle_name == "muon")     return std::make_tuple(13,    muonMass);
  if (particle_name == "pi0")      return std::make_tuple(111,   pionZeroMass);
  if (particle_name == "piplus")   return std::make_tuple(211,   pionMass);
  if (particle_name == "piminus")  return std::make_tuple(-211,  pionMass);

  std::cout << "wrong particle name" << std::endl;
  abort();
}

