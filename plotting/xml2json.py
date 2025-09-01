#!/usr/bin/env python

from collections import defaultdict
import json
import math
from optparse import OptionParser
import re

import dd4hep as dd4hepModule
from ROOT import dd4hep


######################################
# option parser
parser = OptionParser()
parser.add_option('-d', '--detGeoFile',
                  type=str, default='$K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml',
                  help='compact XML detector geometry file')

(options, args) = parser.parse_args()

geo_file_name = options.detGeoFile


######################################
# Functions

skip_pattern = r"(supportTube)|(cryo)"
re_skip = re.compile(skip_pattern)

def get_cells(detector, n_cells = 0):
    sub_detectors = detector.children()
    if sub_detectors.size() == 0:
        # print("Counting", detector.GetName())
        # print("Counting", detector.id())
        n_cells += 1
    else:
        for d in sub_detectors:
            d_name = str(d[0])
            if re_skip.match(d_name):
                print("Skipping sub detector:", d_name)
                continue
            n_cells = get_cells(d[1], n_cells)
    return n_cells


def get_cells_map(detector, sub_det, name):

    cells_map = defaultdict(int)

    # N.B. this could be affected from naming scheme changes
    match name:
        case "EMEC_turbine":
            # Method 1, parse the constants from the geo file
            n_wheels = detector.constantAsLong("nWheels")

            for nw in range(n_wheels):
                # Uses "Calib" instead of ReadOut layers (need to check with a detector expert)
                r_layers = detector.constantAsLong(f"ECalEndcapNumCalibRhoLayersWheel{nw+1}")
                z_layers = detector.constantAsLong(f"ECalEndcapNumCalibZLayersWheel{nw+1}")
                n_units = detector.constantAsLong(f"nUnitCells{nw+1}")

                cells_map[f"wheel{nw+1}"] = r_layers * z_layers * n_units
                cells_map[f"wheel-{nw+1}"] = r_layers * z_layers * n_units

        case "HCalThreePartsEndcap":
            # Split only by positive and negative layers
            for de_name, de in sub_det.children():
                layer_name = str(de_name)
                if "layer-" in layer_name:
                    cells_map["layer_-1"] += get_cells(de)
                elif "layer" in layer_name:
                    cells_map["layer_1"] += get_cells(de)

        case "ECalBarrel":
            # # Method 2, loop only over the LAr "bath" elements
            # bath = sub_det.child("bath")
            # for n, d in bath.children():
            #     # naming scheme 'activeX_Y' where X is the unit number and Y the wheel number
            #     n = str(n)
            #     # if not ("active" in n):
            #     #     continue

            #     # # Get the cells for all the wheels
            #     # nw = int(str(n)[-1])
            #     #cells_map[f"wheel{nw}"] += get_cells(d)
            #     cells_map[n] += get_cells(d)

            n_layers = detector.constantAsLong(f"ECalBarrelNumLayers")
            cells_map =  {f"layer{i}": 1 for i in range(n_layers)}

        case _:
            # Loop over detector elements
            for de_name, de in sub_det.children():
                # layer_id = de.id()
                # print(f"         layer_name (id): {de_name} ({layer_id})")
                if re_skip.match(str(de_name)):
                        print("Skipping sub detector:", de_name)
                        continue
                cells_map[str(de_name)] = get_cells(de)
    return cells_map


##################################################
# Geometry part

# useful doxy:
# https://dd4hep.web.cern.ch/dd4hep/reference/classdd4hep_1_1Detector.html
# https://dd4hep.web.cern.ch/dd4hep/reference/classdd4hep_1_1SensitiveDetector.html
# https://dd4hep.web.cern.ch/dd4hep/reference/classdd4hep_1_1DetElement.html
# https://dd4hep.web.cern.ch/dd4hep/reference/classdd4hep_1_1Volume.html
# https://dd4hep.web.cern.ch/dd4hep/reference/classdd4hep_1_1Box.html
# https://dd4hep.web.cern.ch/dd4hep/reference/DD4hep_2Detector_8h_source.html

# Get detector name
detector = geo_file_name.split("/")[-1].strip(".xml")
# Load geometry
det = dd4hep.Detector.getInstance()
det.fromXML(geo_file_name)
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
    #mySolid=volume.solid()
    #print('dims',mySolid.dimensions())
    #print('IsCylType()',mySolid.IsCylType())
    #print("        dimensions = ", [box.x(), box.y(), box.z()])
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
    # Fill a dictionary to match via ID subdetectors (togheter with its max dimensions)to the SimHits collections
    #if subdet.id()>0.:
    #    binning_string = str(nbins_z)+','+str(-max_z)+'.,'+str(max_z)+'.,'+str(nbins_r)+','+str(-max_r)+'.,'+str(max_r)+'.'
    #    print("        binning = ", binning_string)
    #    dict_sub[str(subdet.id())] ={"binning":  binning_string}

with open(geo_file_name.split('/')[-1].strip('.xml')+'_DetectorDimensions.json', 'w') as fp:
    json.dump(dict_sub, fp, indent=2)
