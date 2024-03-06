#!/usr/bin/env python3
import os
import argparse
from pathlib import Path

from acts.examples import (
    AlignedDetector,
    WhiteBoard,
    AlgorithmContext,
    ProcessCode,
    JsonSurfacesWriter,
    JsonMaterialWriter,
    JsonFormat,
)

import acts
import json

from acts import MaterialMapJsonConverter
import build_geometry_ePIC as geo


def runGeometry(
    trackingGeometry,
    # decorators,
    outputDir="",
    outputName="geometry-map",
    outputObj=False,
    outputCsv=False,
    outputJson=True,
    outputRoot=False,
):
    wb = WhiteBoard(acts.logging.INFO)
    context = AlgorithmContext(0, 0, wb)

    writer = JsonSurfacesWriter(
        level=acts.logging.INFO,
        trackingGeometry=trackingGeometry,
        outputDir=outputDir,
        writePerEvent=True,
        writeSensitive=True,
    )
    writer.write(context)

    jmConverterCfg = MaterialMapJsonConverter.Config(
        processSensitives=True,
        processApproaches=True,
        processRepresenting=True,
        processBoundaries=True,
        processVolumes=True,
        processNonMaterial=True,
        context=context.geoContext,
    )

    outpath = os.path.join(outputDir, outputName)
    jmw = JsonMaterialWriter(
        level=acts.logging.VERBOSE,
        converterCfg=jmConverterCfg,
        fileName=outpath,
        writeFormat=JsonFormat.Json,
    )
    print("Done! Output geometry file at "+outpath+".json")
    jmw.write(trackingGeometry)


if "__main__" == __name__:
    p = argparse.ArgumentParser(
        description="Script to generate geometry-map.json for ePIC geometry"
    )
    p.add_argument(
        "-i","--xmlFile",
        default=os.environ.get("DETECTOR_PATH", "")+"epic_craterlake.xml",
        help="Input xml file containing ePIC geometry",
    )
    p.add_argument(
        "-o","--geoFile",
        type=str,
        default="geometry-map.json",
        help="Output Json file to define volumes and layers used in material mapping",
    )

    args = p.parse_args()
    print(args)
    detector, trackingGeometry, decorators = geo.buildePICGeometry(
        args.xmlFile)

    runGeometry(trackingGeometry, outputDir=os.getcwd(), outputName=args.geoFile.split(".")[0])

    # Uncomment if you want to create the geometry id mapping for DD4hep
    # dd4hepIdGeoIdMap = acts.examples.dd4hep.createDD4hepIdGeoIdMap(trackingGeometry)
    # dd4hepIdGeoIdValueMap = {}
    # for key, value in dd4hepIdGeoIdMap.items():
    #     dd4hepIdGeoIdValueMap[key] = value.value()

    # with open('odd-dd4hep-geoid-mapping.json', 'w') as outfile:
    #    json.dump(dd4hepIdGeoIdValueMap, outfile)