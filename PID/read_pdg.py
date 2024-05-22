import os

import dd4hep
from podio.reading import get_reader

desc = dd4hep.Detector.getInstance()
desc.fromXML("{}/{}.xml".format(os.environ["DETECTOR_PATH"], os.environ["DETECTOR_CONFIG"]))
name_lookup = {
    desc.constantAsLong("BackwardRICH_ID") : "pfRICH",
    desc.constantAsLong("BarrelDIRC_ID") : "hpDIRC",
    desc.constantAsLong("BarrelTOF_ID") : "TOF",
    desc.constantAsLong("ForwardRICH_ID") : "DRICH",
}

def read_pdg(input_file="simu/podio_output.root"):
    reader = get_reader(input_file)
    for i, frame in enumerate(reader.get("events")):
        print(f"event {i}")
        for part in frame.get("ReconstructedChargedParticles"):
            print(f"particle energy {part.getEnergy()}")
            for partid in part.getParticleIDs():
                detector_name = name_lookup.get(partid.getType())
                print(f"\tdetector {detector_name} PDG {partid.getPDG()} probability {partid.getLikelihood()}")

read_pdg()
