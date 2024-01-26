#include <podio/Frame.h>
#include <podio/ROOTFrameReader.h>
#include <edm4eic/ReconstructedParticleCollection.h>

void test(std::string input_file="2024-01-10T21-11-08+00-00_0179803eb6f56c96c459a6af4f187d9d88205f56/rec_dis_18x275_minQ2=1000_craterlake.edm4eic.root") {
    auto reader = podio::ROOTFrameReader();
    reader.openFile(input_file);
    for (size_t i = 0; i < reader.getEntries(podio::Category::Event); i++) {
        auto frame = podio::Frame(reader.readNextEntry(podio::Category::Event));
        std::cerr << "event " << i << std::endl;
        auto &part_coll = frame.get<edm4eic::ReconstructedParticleCollection>("ReconstructedChargedParticles");
        for (auto part : part_coll) {
            std::cerr << "particle momentum magnitude " << part.getMomentum().x << std::endl;;
            std::cerr << "particle Loc-a " << part.getTracks()[0].getTrajectory().getTrackParameters()[0].getLoc()[0] << std::endl;;
            std::cerr << "particle Loc-b " << part.getTracks()[0].getTrajectory().getTrackParameters()[0].getLoc()[1] << std::endl;;
        }
    }
}
