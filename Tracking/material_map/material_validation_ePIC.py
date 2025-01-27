#!/usr/bin/env python3
import os
from typing import Optional
import argparse
from pathlib import Path

from acts.examples import (
    Sequencer,
    ProcessCode,
    RootMaterialTrackWriter,
    JsonFormat,
)

import acts
import acts.examples.dd4hep

from acts import (
    Vector4,
    UnitConstants as u,
    SurfaceMaterialMapper,
    VolumeMaterialMapper,
    Navigator,
    Propagator,
    StraightLineStepper,
    MaterialMapJsonConverter,
)


import json
import build_geometry_ePIC as geo


def runMaterialValidation(
    trackingGeometry,
    decorators,
    outputDir,
    outputName,
    nevents=1000,
    ntests=1000,
    dumpPropagationSteps=False,
    s=None,
):
    s = s or Sequencer(events=nevents, numThreads=-1)

    for decorator in decorators:
        s.addContextDecorator(decorator)
        # for decorator in decorators:
        # assert decorator.decorate(context) == ProcessCode.SUCCESS


    nav = acts.Navigator(
        trackingGeometry=trackingGeometry,
        resolveSensitive=True,
        resolveMaterial=True,
        resolvePassive=True,
    )

    stepper = StraightLineStepper()

    # prop    = Propagator(stepper, nav)
        # mapper = SurfaceMaterialMapper(level=acts.logging.INFO, propagator=propagator)
        # mmAlgCfg.materialSurfaceMapper = mapper

    prop = acts.examples.ConcretePropagator(acts.Propagator(stepper, nav))

    alg = acts.examples.PropagationAlgorithm(
        propagatorImpl=prop,
        level=acts.logging.INFO,
        randomNumberSvc=acts.examples.RandomNumbers(),
        ntests=ntests,
        sterileLogger=False,
        recordMaterialInteractions=True,
        d0Sigma=0,
        z0Sigma=0,
    )

    s.addAlgorithm(alg)

    outpath = os.path.join(outputDir, outputName)
    s.addWriter(
        RootMaterialTrackWriter(
            level=acts.logging.INFO,
            inputMaterialTracks=alg.config.propagationMaterialCollection,
            filePath=outpath,
            storeSurface=True,
            storeVolume=True,
        )
    )
    print("Done! Output at "+outpath)
    outpath = outputDir + "/propagation_steps.root"
    if dumpPropagationSteps:
        s.addWriter(
            acts.examples.RootPropagationStepsWriter(
                level=acts.logging.INFO,
                inputMaterialTracks=alg.config.propagationStepCollection,
                filePath=outpath,
            )
        )
    print("      Propogation steps dumped at "+outpath)

    return s

# python material_validation_ePIC.py --xmlFile epic_craterlake_matmap.xml
if "__main__" == __name__:

    p = argparse.ArgumentParser(
        description="Script to produce propogation validation for ePIC material mapping."
    )
    p.add_argument(
        "--xmlFile",
        default=os.environ.get("DETECTOR_PATH", "")+"epic_craterlake.xml",
        help="input xml file containing ePIC geometry",
    )
    p.add_argument(
        "--matFile",
        type=str,
        default="material-map.json",
        help="input material map file, can be either Json or Cbor",
    )
    p.add_argument(
        "--outputName",
        type=str,
        default="propagation-material.root",
        help="customized name of the output rootfile",
    )
    p.add_argument(
        "-n","--nevents",
        type=int,
        default=100,
        help="number of events to run",
    )
    p.add_argument(
        "-t","--ntests",
        type=int,
        default=100,
        help="number of tests per event",
    )

    args = p.parse_args()
    print(args)
    
    if 'json' in args.matFile:
        fmt = JsonFormat.Json

    elif 'cbor' in args.matFile:
        fmt = JsonFormat.Cbor
    else:
        print('ERROR(material_validation_ePIC.py): please provide a material map file in .json or .cbor format')
        exit()


    detector, trackingGeometry, decorators = geo.buildePICGeometry(args.xmlFile, args.matFile)

    runMaterialValidation(
        trackingGeometry, decorators, outputDir=os.getcwd(), outputName=args.outputName, nevents=args.nevents, ntests=args.ntests
        ).run()

