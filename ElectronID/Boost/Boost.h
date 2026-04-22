// Modified from EICrecon Boost.h

// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2022 Wouter Deconinck, Barak Schmookler

#pragma once

#include <Math/Vector4D.h>
#include <Math/LorentzRotation.h>
#include <Math/LorentzVector.h>
#include <Math/RotationX.h>
#include <Math/RotationY.h>
#include <Math/Boost.h>

using ROOT::Math::PxPyPzEVector;

namespace eicrecon {

using ROOT::Math::LorentzRotation;

inline LorentzRotation determine_boost(PxPyPzEVector ei, PxPyPzEVector pi) {

  using ROOT::Math::Boost;
  using ROOT::Math::RotationX;
  using ROOT::Math::RotationY;

  // Step 1: Find the needed boosts and rotations from the incoming lepton and hadron beams
  // (note, this will give you a perfect boost, in principle you will not know the beam momenta exactly and should use an average)

  PxPyPzEVector eo = ei;
  PxPyPzEVector po = pi;

  // cout << "lab e and p: " << Form("(%f, %f, %f, %f) (%f, %f, %f, %f)", ei.Px(), ei.Py(), ei.Pz(), ei.E(), pi.Px(), pi.Py(), pi.Pz(), pi.E()) << endl;

  // Define the Boost to make beams back-to-back
  const auto cmBoost = (ei + pi).BoostToCM();
  const Boost boost_to_cm(cmBoost);

  // Boost to COM frame
  pi = boost_to_cm(pi);
  ei = boost_to_cm(ei);
  // cout << "CM boosted e and p: " << Form("(%f, %f, %f, %f) (%f, %f, %f, %f)", ei.Px(), ei.Py(), ei.Pz(), ei.E(), pi.Px(), pi.Py(), pi.Pz(), pi.E()) << endl;

  // This will boost beams from a center of momentum frame back to (nearly) their original energies
  PxPyPzEVector eh(0, 0, -1*sqrt(pow(eo.E(),2)-pow(eo.M(),2)), eo.E());
  PxPyPzEVector ph(0, 0,    sqrt(pow(po.E(),2)-pow(po.M(),2)), po.E());
  // cout << "headon e and p: " << Form("(%f, %f, %f, %f) (%f, %f, %f, %f)", eh.Px(), eh.Py(), eh.Pz(), eh.E(), ph.Px(), ph.Py(), ph.Pz(), ph.E()) << endl;

  const auto hoBoost = (eh + ph).BoostToCM();
  // const Boost headon_to_cm(hoBoost);
  const Boost boost_to_headon(-hoBoost);

  // PxPyPzEVector et = headon_to_cm(eh);
  // PxPyPzEVector pt = headon_to_cm(ph);
  // cout << "headon boosted to CM e and p: " << Form("(%f, %f, %f, %f) (%f, %f, %f, %f)", et.Px(), et.Py(), et.Pz(), et.E(), pt.Px(), pt.Py(), pt.Pz(), pt.E()) << endl;

  // PxPyPzEVector erb = boost_to_headon(et);
  // PxPyPzEVector prb = boost_to_headon(pt);
  // cout << "reversed boost e and p: " << Form("(%f, %f, %f, %f) (%f, %f, %f, %f)", erb.Px(), erb.Py(), erb.Pz(), erb.E(), prb.Px(), prb.Py(), prb.Pz(), prb.E()) << endl;
  
  // Boost and rotate the incoming beams to find the proper rotations TLorentzVector

  // Rotate to head-on
  RotationY rotAboutY(-1.0 * atan2(pi.Px(), pi.Pz())); // Rotate to remove x component of beams
  RotationX rotAboutX(+1.0 * atan2(pi.Py(), pi.Pz())); // Rotate to remove y component of beams

  // PxPyPzEVector er = rotAboutX(rotAboutY(ei));
  // PxPyPzEVector pr = rotAboutX(rotAboutY(pi));
  // cout << "rotation e and p: " << Form("(%f, %f, %f, %f) (%f, %f, %f, %f)", er.Px(), er.Py(), er.Pz(), er.E(), pr.Px(), pr.Py(), pr.Pz(), pr.E()) << endl;

  // PxPyPzEVector ef = boost_to_headon(er);
  // PxPyPzEVector pf = boost_to_headon(pr);
  // cout << "final boost e and p: " << Form("(%f, %f, %f, %f) (%f, %f, %f, %f)", ef.Px(), ef.Py(), ef.Pz(), ef.E(), pf.Px(), pf.Py(), pf.Pz(), pf.E()) << endl;

  // cout << "**" << endl;

  // final matrix: P' = [BtoH][RX][RY][BtoCM]P <-- Matrix multi. goes from R to L
  LorentzRotation tf;
  tf *= boost_to_headon;
  tf *= rotAboutX;
  tf *= rotAboutY;
  tf *= boost_to_cm;

  // PxPyPzEVector em = tf(eo);
  // PxPyPzEVector pm = tf(po);
  // cout << "boost matrix e and p: " << Form("(%f, %f, %f, %f) (%f, %f, %f, %f)", em.Px(), em.Py(), em.Pz(), em.E(), pm.Px(), pm.Py(), pm.Pz(), pm.E()) << endl;

  return tf;
}

inline PxPyPzEVector apply_boost(const LorentzRotation& tf, PxPyPzEVector part) {

  // Step 2: Apply boosts and rotations to any particle 4-vector
  // (here too, choices will have to be made as to what the 4-vector is for reconstructed particles)

  // Boost and rotate particle 4-momenta into the headon frame
  tf(part);
  return part;
}

} // namespace eicrecon
