#!/usr/bin/env python3
import os
from typing import Optional
import argparse
from pathlib import Path

from acts.examples import (
    Sequencer,
    WhiteBoard,
    AlgorithmContext,
    ProcessCode,
    RootMaterialTrackReader,
    RootMaterialTrackWriter,
    CsvTrackingGeometryWriter,
    ObjTrackingGeometryWriter,
    MaterialMapping,
    JsonMaterialWriter,
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


def runMaterialMapping(
    trackingGeometry,
    decorators,
    outputDir,
    inputDir,
    inputFile,
    mapName="material-map",
    mapFormat=JsonFormat.Json,
    mapSurface=True,
    mapVolume=False,
    readCachedSurfaceInformation=False,
    mappingStep=1,
    s=None,
):
    s = s or Sequencer(numThreads=1)

    for decorator in decorators:
        s.addContextDecorator(decorator)

    wb = WhiteBoard(acts.logging.INFO)

    context = AlgorithmContext(0, 0, wb)

    for decorator in decorators:
        assert decorator.decorate(context) == ProcessCode.SUCCESS

    # Read material step information from a ROOT TTRee
    s.addReader(
        RootMaterialTrackReader(
            level=acts.logging.INFO,
            outputMaterialTracks="material-tracks",
            fileList=[
                os.path.join(
                    inputDir,
                    inputFile
                )
            ],
            readCachedSurfaceInformation=readCachedSurfaceInformation,
        )
    )

    stepper = StraightLineStepper()

    mmAlgCfg = MaterialMapping.Config(context.geoContext, context.magFieldContext)
    mmAlgCfg.trackingGeometry = trackingGeometry
    mmAlgCfg.inputMaterialTracks = "material-tracks"

    if mapSurface:
        navigator = Navigator(
            trackingGeometry=trackingGeometry,
            resolveSensitive=True, # not required for material mapping
            resolveMaterial=True,
            resolvePassive=True,
        )
        propagator = Propagator(stepper, navigator)
        mapper = SurfaceMaterialMapper(level=acts.logging.INFO, propagator=propagator)
        mmAlgCfg.materialSurfaceMapper = mapper

    if mapVolume:
        navigator = Navigator(
            trackingGeometry=trackingGeometry,
        )
        propagator = Propagator(stepper, navigator)
        mapper = VolumeMaterialMapper(
            level=acts.logging.INFO, propagator=propagator, mappingStep=mappingStep
        )
        mmAlgCfg.materialVolumeMapper = mapper

    jmConverterCfg = MaterialMapJsonConverter.Config(
        processSensitives=True,# not required for material mapping
        processApproaches=True,
        processRepresenting=True,# not required for material mapping
        processBoundaries=True,
        processVolumes=True,
        context=context.geoContext,
    )

    outpath = os.path.join(outputDir, mapName)
    jmw = JsonMaterialWriter(
        level=acts.logging.VERBOSE,
        converterCfg=jmConverterCfg,
        fileName=outpath,
        writeFormat=mapFormat,
    )
    print("Done! Generated material map at "+outpath)

    mmAlgCfg.materialWriters = [jmw]

    s.addAlgorithm(MaterialMapping(level=acts.logging.INFO, config=mmAlgCfg))

    outpath = os.path.join(outputDir,mapName + "_tracks.root")
            
    s.addWriter(
        RootMaterialTrackWriter(
            level=acts.logging.INFO,
            inputMaterialTracks=mmAlgCfg.mappingMaterialCollection,
            filePath=outpath,
            storeSurface=True,
            storeVolume=True,
        )
    )    
    print("      Propogated Geantino tracks at "+outpath)

    return s


if "__main__" == __name__:

    p = argparse.ArgumentParser(
        description="Script to generate material map for ePIC geometry"
    )
    p.add_argument(
        "--xmlFile",
        default=os.environ.get("DETECTOR_PATH", "")+"epic_craterlake.xml",
        help="input xml file containing ePIC geometry",
    )
    p.add_argument(
        "--stepFile",
        type=str,
        default="geant4_material_tracks.root",
        help="input rootfile containing material steps.",
    )
    p.add_argument(
        "--geoFile",
        type=str,
        default="geometry-map.json",
        help="input json file to define volumes and layers used in material mapping",
    )
    p.add_argument(
        "--matFile",
        type=str,
        default="material-map.json",
        help="output filename for the generated material map, can be json and cbor formats",
    )

    # p.add_argument(
    #     '-v',"--run_validation",
    #     action ="store_true",
    #     default=False,
    #     help="Run validation",
    # )
    args = p.parse_args()
    print(args)

    if 'json' in args.matFile:
        fmt = JsonFormat.Json

    elif 'cbor' in args.matFile:
        fmt = JsonFormat.Cbor
    else:
        print('ERROR(material_mapping_ePIC.py): please provide a material map file in .json or .cbor format')
        exit()

    mapName = args.matFile.split('.')[0]

    detector, trackingGeometry, decorators = geo.buildePICGeometry(
        args.xmlFile, args.geoFile)

    runMaterialMapping(
        trackingGeometry,
        decorators,
        outputDir = os.getcwd(),
        inputDir  = os.getcwd(),
        inputFile = args.stepFile,
        readCachedSurfaceInformation=False,
        mapVolume= False,
        mapName  = mapName,
        mapFormat=fmt
    ).run()

    