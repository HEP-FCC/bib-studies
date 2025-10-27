#!/usr/bin/env python

import json
import math
from optparse import OptionParser
import re

import dd4hep as dd4hepModule
from ROOT import dd4hep

from helpers import GeoFile


######################################
# option parser
parser = OptionParser()
parser.add_option('-d', '--detGeoFile',
                  type=str, default='$K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml',
                  help='compact XML detector geometry file')

(options, args) = parser.parse_args()

geo_file = GeoFile(options.detGeoFile)

######################################
# Functions

match geo_file.short:
    case "ALLEGRO":
        from ALLEGRO import get_cells_map
    case "IDEA":
        from IDEA import get_cells_map
    case "CLD":
        from CLD import get_cells_map
    case "ILD_FCCee":
        from ILD_FCCee import get_cells_map
    case _:
        raise NotImplementedError(f"get_cells_map not implemented for detector {geo_file.short}")


##################################################
# Geometry part

# useful doxy:
# https://dd4hep.web.cern.ch/dd4hep/reference/classdd4hep_1_1Detector.html
# https://dd4hep.web.cern.ch/dd4hep/reference/classdd4hep_1_1SensitiveDetector.html
# https://dd4hep.web.cern.ch/dd4hep/reference/classdd4hep_1_1DetElement.html
# https://dd4hep.web.cern.ch/dd4hep/reference/classdd4hep_1_1Volume.html
# https://dd4hep.web.cern.ch/dd4hep/reference/classdd4hep_1_1Box.html
# https://dd4hep.web.cern.ch/dd4hep/reference/DD4hep_2Detector_8h_source.html

# Load geometry
det = dd4hep.Detector.getInstance()
det.fromXML(geo_file.path)

# Get subdetectors
dict_sub = {}

for subdet_name, sens_det in det.sensitiveDetectors():
    print("==== name = ", subdet_name) 

    # Get the DetElement object
    subdet = det.detector(subdet_name)
    print("        ID = ", subdet.id())
    print("        type = ", subdet.type())
    print("        typeFlag = ", subdet.typeFlag())

    hitsCollection = str(sens_det.hitsCollection)
    print("        hits = ", hitsCollection)

    det_element_cells = get_cells_map(det, subdet, subdet_name)
    print("        det_element_cells = ", det_element_cells)

    # Get max dimensions for the subdetector in x,y,z (half-length form 0)
    volume = subdet.volume()
    box = volume.boundingBox()

    # dimensions seem off (ask Brieuc) - for the moment multiply the max by a factor 10
    max_z = 10*math.ceil(box.z()) #max z rounded up
    max_r = 10*math.ceil(math.sqrt(math.pow(box.x(),2)+math.pow(box.y(),2))) #max r rounded up
    print("        max_z = ", max_z)
    print("        max_r = ", max_r)
    if(box.x()!=box.y()): print ('WARNING: different X/Y in the detector bounding box X:',box.x(),'Y:',box.y())

    if subdet.id()>0.:
        dict_sub[str(subdet_name)]={
            'id': int(subdet.id()),
            'typeFlag': subdet.typeFlag(),
            'hitsCollection': hitsCollection,
            'det_element_cells': det_element_cells,
            'max_z': max_z, 
            'max_r': max_r,
        }

with open(geo_file.name+'_DetectorDimensions.json', 'w') as fp:
    json.dump(dict_sub, fp, indent=2)
