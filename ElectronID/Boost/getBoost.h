#pragma once

#include "Beam.h"
#include "Boost.h"

#include <Math/LorentzVector.h>
using ROOT::Math::PxPyPzEVector;

#include <Math/LorentzRotation.h>
using ROOT::Math::LorentzRotation;

LorentzRotation getBoost(double eE, double eN, double mE, double mN) {

    TVector3 ve(0,0,-sqrt(eE*eE-mE*mE));
    TVector3 vn(0,0,sqrt(eN*eN-mN*mN));
    
    const PxPyPzEVector ei(
        eicrecon::round_beam_four_momentum(
            ve,
            mE,
            {-1*eE},
            0.0)
        );

    const PxPyPzEVector ni(
        eicrecon::round_beam_four_momentum(
            vn,
            mN,
            {eN},
            -0.025)
        );

    // std::cout << "initial electron beam: " << Form("(%f, %f, %f, %f)", ei.Px(), ei.Py(), ei.Pz(), ei.E()) << std::endl;
    // std::cout << "initial nucleon beam: " << Form("(%f, %f, %f, %f)", ni.Px(), ni.Py(), ni.Pz(), ni.E()) << std::endl;
    // std::cout << "Created vector, e, p (mrad) " << ei.Theta() * 1000 << " " << ni.Theta() * 1000 << std::endl;

    LorentzRotation boost = eicrecon::determine_boost(ei, ni);

    PxPyPzEVector boosted_e = boost(ei);
    PxPyPzEVector boosted_n = boost(ni);

    // std::cout << "Boosted electron beam: " << Form("(%f, %f, %f, %f)", boosted_e.Px(), boosted_e.Py(), boosted_e.Pz(), boosted_e.E()) << std::endl;
    // std::cout << "Boosted nucleon beam: " << Form("(%f, %f, %f, %f)", boosted_n.Px(), boosted_n.Py(), boosted_n.Pz(), boosted_n.E()) << std::endl;
    // std::cout << "Incoming headon, e, p (mrad) " << boosted_e.Theta() * 1000 << " " << boosted_n.Theta() * 1000 << std::endl;

    return boost;
}
