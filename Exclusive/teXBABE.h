// Stephen JD Kay - University of York
// Stephen.Kay@york.ac.uk
// 01/07/25
// Header file for a generic method of calculating teXBABE as defined by the tRECO convention
// This offers an alternative to ePIC_ExcKinUtils where ONLY teXBABE is defined and calculated - https://github.com/eic/snippets/blob/main/Exclusive/ePIC_ExcKinUtils.h
// Method expects a three body final state where the expected baryon (BA) state is known and the scattered electron (eSc) and other hadronic state(s) X, are measured and well known

PxPyPzEVector Vec_BACorr;

// Explicit method to calculate teXBABE. Inputs are the detected scattered electron, detected X, detected baryon (BA), mass of the expected baryon, Electron beam and hadron beam
// All inputs OTHER than the BA mass are expected as PxPyPzEVectors
Double_t Calc_teXBABE(PxPyPzEVector Vec_eSc, PxPyPzEVector Vec_X, PxPyPzEVector Vec_BA, Double_t Mass_BA, PxPyPzEVector Vec_EBeam, PxPyPzEVector Vec_HBeam){
  PxPyPzEVector Vec_PMiss = (Vec_EBeam + Vec_HBeam) - (Vec_eSc + Vec_X);
  // "Correct" the BA vector - utilise knowledge of the final state and the POSITION resolution of the original BA track
  Vec_BACorr.SetPxPyPzE(Vec_PMiss.P()*sin(Vec_BA.Theta())*cos(Vec_BA.Phi()), Vec_PMiss.P()*sin(Vec_BA.Theta())*sin(Vec_BA.Phi()), Vec_PMiss.P()*cos(Vec_BA.Theta()), sqrt(pow(Vec_PMiss.P(),2)+(pow(Mass_BA,2))));
  Vec_teXBABE = ;
  Double_t t_eXBABE = -1*((Vec_HBeam - Vec_BACorr).mag2());
  return(t_eXBABE);
}

// Shortened method to calculate teXBABE where inputs are PMiss, BA vector, BA mass and HBeam only
Double_t Calc_teXBABE(PxPyPzEVector Vec_PMiss, PxPyPzEVector Vec_BA, Double_t Mass_BA, PxPyPzEVector Vec_HBeam){
  // "Correct" the BA vector - utilise knowledge of the final state and the POSITION resolution of the original BA track
  Vec_BACorr.SetPxPyPzE(Vec_PMiss.P()*sin(Vec_BA.Theta())*cos(Vec_BA.Phi()), Vec_PMiss.P()*sin(Vec_BA.Theta())*sin(Vec_BA.Phi()), Vec_PMiss.P()*cos(Vec_BA.Theta()), sqrt(pow(Vec_PMiss.P(),2)+(pow(Mass_BA,2))));
  Double_t t_eXBABE = -1*((Vec_HBeam - Vec_BACorr).mag2());
  return(t_eXBABE);
}

// Shortest method, inputs are corrected BA vector and HBeam only
Double_t Calc_teXBABE(PxPyPzEVector Vec_BACorr, PxPyPzEVector Vec_HBeam){
  Double_t t_eXBABE = -1*((Vec_HBeam - Vec_BACorr).mag2());
  return(t_eXBABE);
}
