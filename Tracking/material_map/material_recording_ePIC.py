#!/usr/bin/env python3

import os
import warnings
from pathlib import Path
import argparse

import acts
from acts.examples import (
    GaussianVertexGenerator,
    ParametricParticleGenerator,
    FixedMultiplicityGenerator,
    EventGenerator,
    RandomNumbers,
)

import acts.examples.dd4hep
import acts.examples.geant4
import acts.examples.geant4.dd4hep
from acts.examples.odd import getOpenDataDetector
from acts.examples.geant4 import GdmlDetectorConstructionFactory

u = acts.UnitConstants

import build_geometry_ePIC as geo

_material_recording_executed = False


def runMaterialRecording(
    detectorConstructionFactory,
    outputDir,
    outputName="geant4_material_tracks.root",
    tracksPerEvent=1000,
    s=None,
    etaRange=(-8, 8),
):
    global _material_recording_executed
    if _material_recording_executed:
        warnings.warn("Material recording already ran in this process. Expect crashes")
    _material_recording_executed = True

    rnd = RandomNumbers(seed=228)

    evGen = EventGenerator(
        level=acts.logging.INFO,
        generators=[
            EventGenerator.Generator(
                multiplicity=FixedMultiplicityGenerator(n=1),
                vertex=GaussianVertexGenerator(
                    stddev=acts.Vector4(0, 0, 0, 0),
                    mean=acts.Vector4(0, 0, 0, 0),
                ),
                particles=ParametricParticleGenerator(
                    pdg=acts.PdgParticle.eInvalid,
                    charge=0,
                    randomizeCharge=False,
                    mass=0,
                    p=(1 * u.GeV, 10 * u.GeV),
                    eta=etaRange,
                    numParticles=tracksPerEvent,
                    etaUniform=True,
                ),
            )
        ],
        outputParticles="particles_initial",
        randomNumbers=rnd,
    )

    s.addReader(evGen)

    g4Alg = acts.examples.geant4.Geant4MaterialRecording(
        level=acts.logging.INFO,
        detectorConstructionFactory=detectorConstructionFactory,
        randomNumbers=rnd,
        inputParticles=evGen.config.outputParticles,
        outputMaterialTracks="material-tracks",
    )

    s.addAlgorithm(g4Alg)

    outpath=os.path.join(outputDir, outputName)
    s.addWriter(
        acts.examples.RootMaterialTrackWriter(
            prePostStep=True,
            recalculateTotals=True,
            inputMaterialTracks="material-tracks",
            filePath=outpath,
            level=acts.logging.INFO,
        )
    )
    print("Done! Recorded steps in "+outpath)

    return s

# python material_recording_ePIC.py -i epic_craterlake_matmap.xml -n 1000 -t 1000
if "__main__" == __name__:

    p = argparse.ArgumentParser(description="Record the Geant4 materials with geantino")
    p.add_argument(
        "-n", "--nevents", type=int, default=1000, help="Number of events to generate"
    )
    p.add_argument(
        "-t", "--ntracks", type=int, default=1000, help="Particle tracks per event"
    )
    p.add_argument(
        "-i", "--xmlFile", type=str, default=os.environ.get("DETECTOR_PATH", "")+"epic_craterlake.xml", help="DD4hep input file"
    )
    p.add_argument(
        "-o", "--outputName", type=str, default="geant4_material_tracks.root", help="Name of the output rootfile"
    )
    p.add_argument(
        "--eta_min",
        type=float,
        default=-8.0,
        help="eta min (optional)",
    )
    p.add_argument(
        "--eta_max",
        type=float,
        default=8.0,
        help="eta max (optional)",
    )
    args = p.parse_args()

    detectorConstructionFactory = None

    detector, trackingGeometry, decorators = geo.buildePICGeometry(
        args.xmlFile)

    detectorConstructionFactory = (
            acts.examples.geant4.dd4hep.DDG4DetectorConstructionFactory(detector)
        )

    runMaterialRecording(
        detectorConstructionFactory=detectorConstructionFactory,
        tracksPerEvent=args.ntracks,
        outputDir=os.getcwd(),
        outputName=args.outputName,
        etaRange=(args.eta_min, args.eta_max),
        s=acts.examples.Sequencer(events=args.nevents, numThreads=1),
    ).run()

