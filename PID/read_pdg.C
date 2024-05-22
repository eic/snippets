#include <edm4eic/ReconstructedParticleCollection.h>
#include <podio/Frame.h>
#include <podio/ROOTFrameReader.h>

void read_pdg(std::string input_file="simu/podio_output.root") {
    auto reader = podio::ROOTFrameReader();
    reader.openFile(input_file);
    for (size_t i = 0; i < reader.getEntries(podio::Category::Event); i++) {
        auto frame = podio::Frame(reader.readNextEntry(podio::Category::Event));
        std::cerr << "event " << i << std::endl;
        auto &part_coll = frame.get<edm4eic::ReconstructedParticleCollection>("ReconstructedChargedParticles");
        for (auto part : part_coll) {
            std::cerr << "particle energy " << part.getEnergy() << std::endl;;
            std::cerr << "PID PDG " << part.getPDG() << std::endl;;
	    for (auto partid : part.getParticleIDs()) {
                std::cerr << "\tdetector " << partid.getType()
                          << " PDG " << partid.getPDG()
			  << " probability " << partid.getLikelihood();
	    }
        }
    }
}
