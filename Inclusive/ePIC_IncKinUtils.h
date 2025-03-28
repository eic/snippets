//---------------------------------------------------------------------------------------
//
// Utility functions; calculation of kinematic quantities for inclusive ePIC analyses
//
// Q2, x and y for all 5 methods used for the InclusiveKinematics branch of EICrecon
//     - Electron, JB, DA, Sigma and eSigma
//
// Author: O. Jevons, 27/02/25
//
//---------------------------------------------------------------------------------------

#include "TMath.h"

// Aliases for common 3/4-vector types
using P3EVector=ROOT::Math::PxPyPzEVector;
using P3MVector=ROOT::Math::PxPyPzMVector;
using MomVector=ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<Double_t>,ROOT::Math::DefaultCoordinateSystemTag>;

//-----------------------------------------------------------------------------------------------------------------------------
// FUNCTION DEFINITIONS
//
// NOTE: Templating applied for brevity
//       Functions are valid for any type which contains the operators X(), Y(), Z(), E(), Pt(), M2() and Dot()
//       This includes TLorentzVector (legacy) and ALL variants of ROOT::Math::LorentzVector class
//
//
// NOTE 2: Order of vectors (whichever are required) is ALWAYS as follows (for process ep -> e'p'X - change scattered baryon as necessary)
//
//         e -> p -> e' -> p' -> X
//
//-----------------------------------------------------------------------------------------------------------------------------

// Calculate DIS kinematics - electron method
// e + p -> e' + X
// Q2
template <typename V>
Double_t calcQ2_Elec(const V& e, const V& ep){
  return -(e-ep).M2();
}
// x
template <typename V>
Double_t calcX_Elec(const V& e, const V& p, const V& ep){
  P3EVector q(e.X()-ep.X(), e.Y()-ep.Y(), e.Z()-ep.Z(), e.E()-ep.E());
  Double_t q2 = -q.Dot(q);
  Double_t den = 2*q.Dot(p);

  return q2/den;
}
// y
template <typename V>
Double_t calcY_Elec(const V& e, const V& p, const V& ep){
  P3EVector q(e.X()-ep.X(), e.Y()-ep.Y(), e.Z()-ep.Z(), e.E()-ep.E());
  Double_t num = p.Dot(q);
  Double_t den = p.Dot(e);
  
  return num/den;
}
// Calculate all 3 kinematic variables (Q2, x, y) simultaneously and pass them to holding variables
template <typename V>
void calcKin_Elec(const V& e, const V& p, const V& ep, Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;

  // Calculate kinematics
  P3EVector q(e.X()-ep.X(), e.Y()-ep.Y(), e.Z()-ep.Z(), e.E()-ep.E());
  Float_t Q2_e = -q.Dot(q);
  
  Float_t q_dot_p = q.Dot(p);
  Float_t x_e = Q2_e/(2*q_dot_p);

  Float_t e_dot_p = e.Dot(p);
  Float_t y_e = q_dot_p/e_dot_p;

  // Export variables
  Q2 = Q2_e;
  x = x_e;
  y = y_e;
}

// Calculate inclusive kinematics (Q2, x, y) with the JB method
// SIDIS case ep -> e'p'X
// Q2
template <typename V>
Double_t calcQ2_JB(const V& e, const V& pp, const V& X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());
  Float_t Q2_jb = pT2_h / (1-y_jb);

  return Q2_jb;
}
// x
template <typename V>
Double_t calcX_JB(const V& e, const V& p, const V& pp, const V& X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());
  Float_t Q2_jb = pT2_h / (1-y_jb);
  Float_t x_jb = Q2_jb / (4*e.E()*p.E()*y_jb);
  
  return x_jb;
}
// y
template <typename V>
Double_t calcY_JB(const V& e, const V& pp, const V& X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();

  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());

  return y_jb;
}
// Calculate all 3 kinematic variables (Q2, x, y) simultaneously and pass them to holding variables 
template <typename V>
void calcKin_JB(const V& e, const V& p, const V& pp, const V& X, Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
  
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());
  Float_t Q2_jb = pT2_h / (1-y_jb);
  Float_t x_jb = Q2_jb / (4*e.E()*p.E()*y_jb);
  
  // Export kinematic variables
  Q2 = Q2_jb;
  x = x_jb;
  y = y_jb;
}

// DIS case ep -> e'X
// Q2
template <typename V>
Double_t calcQ2_JB(const V& e, const V& HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());
  Float_t Q2_jb = pT2_h / (1-y_jb);

  return Q2_jb;
}
// x
template <typename V>
Double_t calcX_JB(const V& e, const V& p, const V& HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());
  Float_t Q2_jb = pT2_h / (1-y_jb);
  Float_t x_jb = Q2_jb / (4*e.E()*p.E()*y_jb);
  
  return x_jb;
}
// y
template <typename V>
Double_t calcY_JB(const V& e, const V& HFS){
  // Intermediate variables
  Float_t sigma_h = HFS.E() - HFS.Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());

  return y_jb;
}
// Calculate all 3 kinematic variables (Q2, x, y) simultaneously and pass them to holding variables 
template <typename V>
void calcKin_JB(const V& e, const V& p, const V& HFS, Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
  
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  
  // Kinematic variables
  Float_t y_jb = sigma_h / (2*e.E());
  Float_t Q2_jb = pT2_h / (1-y_jb);
  Float_t x_jb = Q2_jb / (4*e.E()*p.E()*y_jb);
  
  // Export kinematic variables
  Q2 = Q2_jb;
  x = x_jb;
  y = y_jb;
}

// Calculate inclusive kinematics (Q2, x, y) with the DA method
// SIDIS case ep -> e'p'X
// Q2
template <typename V>
Double_t calcQ2_DA(const V& e, const V& ep, const V& pp, const V& X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t Q2_da = 4*e.E()*e.E()*(1./TMath::Tan(ep.Theta()/2))*(1./(TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2)));

  return Q2_da;
}
// x
template <typename V>
Double_t calcX_DA(const V& e, const V& p, const V& ep, const V& pp, const V& X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t y_da = TMath::Tan(gamma/2) / (TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2));
  Float_t Q2_da = 4*e.E()*e.E()*(1./TMath::Tan(ep.Theta()/2))*(1./(TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2)));
  Float_t x_da = Q2_da/(4*e.E()*p.E()*y_da);

  return x_da;
}
// y
template <typename V>
Double_t calcY_DA(const V& ep, const V& pp, const V& X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t y_da = TMath::Tan(gamma/2) / (TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2));

  return y_da;
}
// Calculate all 3 kinematic variables (Q2, x, y) simultaneously and pass them to holding variables 
template <typename V>
void calcKin_DA(const V& e, const V& p, const V& ep, const V& pp, const V& X, Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
 
  // Calculate intermediate quantities
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t y_da = TMath::Tan(gamma/2) / (TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2));
  Float_t Q2_da = 4*e.E()*e.E()*(1./TMath::Tan(ep.Theta()/2))*(1./(TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2)));
  Float_t x_da = Q2_da/(4*e.E()*p.E()*y_da);

  // Export kinematic variables
  Q2 = Q2_da;
  x = x_da;
  y = y_da;
}

// DIS case ep -> e'X
// Q2
template <typename V>
Double_t calcQ2_DA(const V& e, const V& ep, const V& HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t Q2_da = 4*e.E()*e.E()*(1./TMath::Tan(ep.Theta()/2))*(1./(TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2)));

  return Q2_da;
}
// x
template <typename V>
Double_t calcX_DA(const V& e, const V& p, const V& ep, const V& HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t y_da = TMath::Tan(gamma/2) / (TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2));
  Float_t Q2_da = 4*e.E()*e.E()*(1./TMath::Tan(ep.Theta()/2))*(1./(TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2)));
  Float_t x_da = Q2_da/(4*e.E()*p.E()*y_da);

  return x_da;
}
// y
template <typename V>
Double_t calcY_DA(const V& ep, const V& HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t y_da = TMath::Tan(gamma/2) / (TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2));

  return y_da;
}
// Calculate all 3 kinematic variables (Q2, x, y) simultaneously and pass them to holding variables 
template <typename V>
void calcKin_DA(const V& e, const V& p, const V& ep, const V& HFS, Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
 
  // Calculate intermediate quantities
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t S2_h = sigma_h*sigma_h;
  Float_t gamma = TMath::ACos((pT2_h-S2_h)/(pT2_h+S2_h));
  
  // Kinematic variables
  Float_t y_da = TMath::Tan(gamma/2) / (TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2));
  Float_t Q2_da = 4*e.E()*e.E()*(1./TMath::Tan(ep.Theta()/2))*(1./(TMath::Tan(ep.Theta()/2) +  TMath::Tan(gamma/2)));
  Float_t x_da = Q2_da/(4*e.E()*p.E()*y_da);

  // Export kinematic variables
  Q2 = Q2_da;
  x = x_da;
  y = y_da;
}

// Calculate inclusive kinematics (Q2, x, y) with the Sigma method
// SIDIS case ep -> e'p'X
// Q2
template <typename V>
Double_t calcQ2_Sigma(const V& ep, const V& pp, const V& X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  
  return Q2_sigma;
}
// x
template <typename V>
Double_t calcX_Sigma(const V& e, const V& p, const V& ep, const V& pp, const V& X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);
  
  return x_sigma;
}
// y
template <typename V>
Double_t calcY_Sigma(const V& ep, const V& pp, const V& X){
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;

  return y_sigma;
}
// Calculate all 3 kinematic variables (Q2, x, y) simultaneously and pass them to holding variables 
template <typename V>
void calcKin_Sigma(const V& e, const V& p, const V& ep, const V& pp, const V& X, Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
  
  // Intermediate variables
  Float_t pT2_h = TMath::Power((pp+X).X(),2) + TMath::Power((pp+X).Y(),2);
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);
  
  // Export kinematic variables
  Q2 = Q2_sigma;
  x = x_sigma;
  y = y_sigma;
}

// DIS case ep -> e'X
// Q2
template <typename V>
Double_t calcQ2_Sigma(const V& ep, const V& HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  
  return Q2_sigma;
}
// x
template <typename V>
Double_t calcX_Sigma(const V& e, const V& p, const V& ep, const V& HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);
  
  return x_sigma;
}
// y
template <typename V>
Double_t calcY_Sigma(const V& ep, const V& HFS){
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;

  return y_sigma;
}
// Calculate all 3 kinematic variables (Q2, x, y) simultaneously and pass them to holding variables 
template <typename V>
void calcKin_Sigma(const V& e, const V& p, const V& ep, const V& HFS, Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
  
  // Intermediate variables
  Float_t pT2_h = TMath::Power(HFS.X(),2) + TMath::Power(HFS.Y(),2);
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Kinematic variables
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);
  
  // Export kinematic variables
  Q2 = Q2_sigma;
  x = x_sigma;
  y = y_sigma;
}

// Calculate inclusive kinematics (Q2, x, y) with the e-Sigma method
// SIDIS case ep -> e'p'X
// Q2
template <typename V>
Double_t calcQ2_ESigma(const V& e, const V& ep){
  // Intermediate variables
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();

  // Electron-based kinematics
  Float_t y_e = 1 - sigma_e / (2*e.E());
  Float_t Q2_e = pT2_e / (1-y_e);

  return Q2_e;
}
// x
template <typename V>
Double_t calcX_ESigma(const V& e, const V& p, const V& ep, const V& pp, const V& X){
  // Intermediate variables
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Hadron-based kinematics
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);

  return x_sigma;
}
// y
template <typename V>
Double_t calcY_ESigma(const V& e, const V& p, const V& ep, const V& pp, const V& X){
  // Intermediate variables
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Electron-based kinematics
  Float_t y_e = 1 - sigma_e / (2*e.E());
  Float_t Q2_e = pT2_e / (1-y_e);

  // Hadron-based kinematics
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);

  // Mixed kinematics
  Float_t y_esigma = Q2_e / (4*e.E()*p.E()*x_sigma);

  return y_esigma;
}
// Calculate all 3 kinematic variables (Q2, x, y) simultaneously and pass them to holding variables 
template <typename V>
void calcKin_ESigma(const V& e, const V& p, const V& ep, const V& pp, const V& X, Float_t &Q2, Float_t &x, Float_t &y){
  // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
  
  // Intermediate variables
  Float_t sigma_h = (pp+X).E() - (pp+X).Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Electron-based kinematics
  Float_t y_e = 1 - sigma_e / (2*e.E());
  Float_t Q2_e = pT2_e / (1-y_e);

  // Hadron-based kinematics
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);

  // Mixed kinematics
  Float_t y_esigma = Q2_e / (4*e.E()*p.E()*x_sigma);

  // Export kinematics
  Q2 = Q2_e;
  x = x_sigma;
  y = y_esigma;
}

// DIS case ep -> e'X
// NO NEED TO REDEFINE Q2 - ONLY DEPENDS ON BEAM AND SCATTERED ELECTRON
// x
template <typename V>
Double_t calcX_ESigma(const V& e, const V& p, const V& ep, const V& HFS){
  // Intermediate variables
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Hadron-based kinematics
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);

  return x_sigma;
}
// y
template <typename V>
Double_t calcY_ESigma(const V& e, const V& p, const V& ep, const V& HFS){
  // Intermediate variables
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Electron-based kinematics
  Float_t y_e = 1 - sigma_e / (2*e.E());
  Float_t Q2_e = pT2_e / (1-y_e);

  // Hadron-based kinematics
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);

  // Mixed kinematics
  Float_t y_esigma = Q2_e / (4*e.E()*p.E()*x_sigma);

  return y_esigma;
}
// Calculate all 3 kinematic variables (Q2, x, y) simultaneously and pass them to holding variables 
template <typename V>
void calcKin_ESigma(const V& e, const V& p, const V& ep, const V& HFS, Float_t &Q2, Float_t &x, Float_t &y){
   // Reset kinematic variables
  Q2 = 0;
  x = 0;
  y = 0;
  
  // Intermediate variables
  Float_t sigma_h = HFS.E() - HFS.Z();
  Float_t pT2_e = ep.Pt() * ep.Pt();
  Float_t sigma_e = ep.E() - ep.Z();
  Float_t sigma_tot = sigma_e + sigma_h;

  // Electron-based kinematics
  Float_t y_e = 1 - sigma_e / (2*e.E());
  Float_t Q2_e = pT2_e / (1-y_e);

  // Hadron-based kinematics
  Float_t y_sigma = sigma_h/sigma_tot;
  Float_t Q2_sigma = pT2_e / (1-y_sigma);
  Float_t x_sigma = Q2_sigma / (4*e.E()*p.E()*y_sigma);

  // Mixed kinematics
  Float_t y_esigma = Q2_e / (4*e.E()*p.E()*x_sigma);

  // Export kinematics
  Q2 = Q2_e;
  x = x_sigma;
  y = y_esigma;
}
