#include <iostream>
#include "TMath.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TChain.h"
#include "TH2.h"
#include "TVector3.h"
#include <fstream>
#include <string>
#include <fstream>
#include <unordered_map>

void writeCALOROCToTxt(const std::string& input,
                       const std::string& outputTxt) {
    TChain chain("events");
    chain.Add(input.c_str());

    TTreeReader reader(&chain);

    TTreeReaderArray<int> genStatus(reader, "MCParticles.generatorStatus");
    TTreeReaderArray<double> momX(reader, "MCParticles.momentum.x");
    TTreeReaderArray<double> momY(reader, "MCParticles.momentum.y");
    TTreeReaderArray<double> momZ(reader, "MCParticles.momentum.z");
    TTreeReaderArray<int> PDG(reader, "MCParticles.PDG");
    TTreeReaderArray<double> mass(reader, "MCParticles.mass");

    TTreeReaderArray<unsigned long> CellID(reader, "EcalBarrelScFiRecHits.cellID");
    TTreeReaderArray<float> CALOROCE(reader, "EcalBarrelScFiRecHits.energy");
    TTreeReaderArray<float> CALOROCX(reader, "EcalBarrelScFiRecHits.position.x");
    TTreeReaderArray<float> CALOROCY(reader, "EcalBarrelScFiRecHits.position.y");
    TTreeReaderArray<float> CALOROCZ(reader, "EcalBarrelScFiRecHits.position.z");

    TTreeReaderArray<unsigned long> ScfiCellID(reader, "EcalBarrelScFiPNpeHits.cellID");
    TTreeReaderArray<unsigned int> ScfiBegin(reader, "EcalBarrelScFiPNpeHits.contributions_begin");
    TTreeReaderArray<unsigned int> ScfiEnd(reader, "EcalBarrelScFiPNpeHits.contributions_end");
    TTreeReaderArray<int> ScfiMCcId(reader, "_EcalBarrelScFiPNpeHits_contributions.index");
    TTreeReaderArray<int> ScfiMCId(reader, "_EcalBarrelScFiPAttenuatedHitContributions_particle.index");


    std::ofstream output(outputTxt.c_str());
    std::unordered_map<int, double> dp;

    int nevent = 0;
    while(reader.Next()) {
        if(nevent % 10 == 0) std::cout << "Processing event " << nevent << std::endl;
        bool emptyEvent = true;

        std::vector<int> Att2MC(ScfiCellID.GetSize(), -1);
        std::unordered_map<unsigned long, int> cellID2AttId;

        for(size_t i = 0; i < ScfiCellID.GetSize(); ++i){
            for(size_t j = ScfiBegin[i]; j < ScfiEnd[i]; ++j){
                int cid = ScfiMCcId[j];
                Att2MC[i] = ScfiMCId[cid]; // BIG ASSUMPATIONS: ScfiMCId[j] all returns the same id for ScfiBegin[i] <= j < ScfiEnd[i]
            }
            cellID2AttId[ScfiCellID[i]] = i;
        }


        for(size_t i = 0; i < CellID.GetSize(); ++i) {
            auto it = cellID2AttId.find(CellID[i]);
            if(it == cellID2AttId.end()) continue;

            int attId = it -> second;
            int mcId = Att2MC[attId];
            output << CALOROCX[i] << " " 
                   << CALOROCY[i] << " " 
                   << CALOROCZ[i] << " " 
                   << CALOROCE[i] << " " 
                   << CALOROCZ[i] << " " 
                   << genStatus[mcId] << " " 
                   << momX[mcId] << " " 
                   << momY[mcId] << " " 
                   << momZ[mcId] << " " 
                   << PDG[mcId] << " " 
                   << mcId << " " 
                   << CALOROCE[i] 
                   << " " << nevent << std::endl;
        }
        for(int i = 0; i < 12; ++i) output << "-9999 ";
        output << nevent << std::endl; // this serves as line break for next event
        ++nevent;
    }
}
