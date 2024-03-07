## Stand alone function to build ePIC geometry with ACTS python bindings
## for material mapping
## Shujie Li, 03, 2024
#!/usr/bin/env python3
import argparse
from pathlib import Path

# from acts.examples import (
#     # Sequencer,
#     # WhiteBoard,
#     # AlgorithmContext,
#     # ProcessCode,
#     # RootMaterialTrackReader,
#     # RootMaterialTrackWriter,
#     # CsvTrackingGeometryWriter,
#     # ObjTrackingGeometryWriter,
#     MaterialMapping,
#     JsonMaterialWriter,
#     JsonSurfacesWriter,
#     JsonFormat,
# )

import acts
import acts.examples.dd4hep

from acts import (
    Vector4,
    # UnitConstants as u,
    # SurfaceMaterialMapper,
    # VolumeMaterialMapper,
    # Navigator,
    # Propagator,
    # StraightLineStepper,
    MaterialMapJsonConverter
)

import json

def buildePICGeometry(
    xmlFile,
    jsonFile="",
    logLevel=acts.logging.WARNING,
):
    customLogLevel = acts.examples.defaultLogging(logLevel=logLevel)
    logger = acts.logging.getLogger("BuildePICGeometry")

    matDeco = None
    if len(jsonFile)>0:
        file = Path(jsonFile)
        logger.info("Adding material from %s", file.absolute())
        matDeco = acts.IMaterialDecorator.fromFile(
            file,
            level=customLogLevel(maxLevel=acts.logging.INFO),
        )

    dd4hepConfig = acts.examples.dd4hep.DD4hepGeometryService.Config(
        xmlFileNames=[xmlFile],
        logLevel=logLevel,
        dd4hepLogLevel=customLogLevel(),
    )
    detector = acts.examples.dd4hep.DD4hepDetector()

    config = acts.MaterialMapJsonConverter.Config()

    trackingGeometry, deco = detector.finalize(dd4hepConfig, matDeco)

    return detector, trackingGeometry, deco