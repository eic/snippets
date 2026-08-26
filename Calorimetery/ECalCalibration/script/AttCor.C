#include <cmath>
#include <fstream>
#include <iostream>
#include "TMath.h"
#include "TDatabasePDG.h"
#include "TVector3.h"
#include <unordered_map>

void AttCor(const std::string& side,
            const std::string& inputTxt,
            const std::string& outputTxt,
            double thetaMin = 0,
            double thetaMax = 180) {
    std::ifstream input(inputTxt.c_str());
    std::ofstream output(outputTxt.c_str());

    // read input
    double x, y, temp, npe, z, momx, momy, momz;
    int status, mcId, pdg, nevent;

    // convert pdg to mass
    auto pdgDB = TDatabasePDG::Instance();

    // sum energies
    std::unordered_map<int, double> mcId2TotE, mcId2Pdg, mcId2Mom, mcIdX, mcIdY, mcIdZ;
    std::unordered_map<int, int> mcStatus;
    double prevMomx = 0, prevMomy = 0, prevMomz = 0;
    while((input >> x >> y >> temp >> npe >> z >> status >> momx >> momy >> momz >> pdg >> mcId >> temp >> nevent)) {
         if(temp < -999) {
              if(nevent % 100 == 0) std::cout << "Processing event " << nevent << std::endl;
              // end of event line
              // summarized the event
              for(const auto& key : mcId2TotE) {
                  TVector3 vec(mcIdX[key.first], mcIdY[key.first], mcIdZ[key.first]);
                  if(thetaMin <= vec.Theta()*TMath::RadToDeg() && vec.Theta()*TMath::RadToDeg() <= thetaMax) {
                      auto particle = pdgDB -> GetParticle(mcId2Pdg[key.first]);
                      double mass = 0;
                      if(particle) mass = particle -> Mass(); // in GeV/c^2
                      double totE = std::sqrt(mcId2Mom[key.first]*mcId2Mom[key.first] + mass*mass);
                      output << key.first << " " 
                             << mcId2Pdg[key.first] << " " 
                             << vec.Theta()*TMath::RadToDeg() << " " 
                             << key.second << " " 
                             << totE << " " 
                             << vec.Phi() * TMath::RadToDeg() << " "
                             << vec.Eta() << " " 
                             << vec.Phi() * TMath::RadToDeg() << " "
                             << vec.Eta() << " "
                             << mcStatus[key.first] << " " 
                             << nevent << std::endl;
                  }
              }
              output << "-9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 " << nevent << std::endl;
              mcId2TotE.clear();
              mcId2Pdg.clear();
              mcIdX.clear();
              mcIdY.clear();
              mcIdZ.clear();
          }
          else {
              mcIdX[mcId] += x;
              mcIdY[mcId] += y;
              mcIdZ[mcId] += z;
              mcId2Mom[mcId] = TVector3(momx, momy, momz).Mag();
              mcId2TotE[mcId] += npe;
              mcStatus[mcId] = status;
              mcId2Pdg[mcId] = pdg;
          }
     }

}
