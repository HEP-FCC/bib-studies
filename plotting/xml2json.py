#!/usr/bin/env python

import dd4hep as dd4hepModule
from ROOT import dd4hep
from optparse import OptionParser
import math
import os

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

def get_cells(detector, n_cells = 0):
    sub_detectors = detector.children() # std::map<string,det>
    if sub_detectors.size() == 0:
        #print("increase", n_cells)
        n_cells += 1
    else:
        for d in sub_detectors:
            n_cells = get_cells(d[1], n_cells)
    return n_cells

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
#subdet = det.detector("VertexBarrel")
#subdets = det.detectors();
subdets = det.sensitiveDetectors();
for subdet_name, sens_det in subdets:
    print("==== name = ", subdet_name) 

    # Get the DetElement object
    subdet = det.detector(subdet_name)
    print("        ID = ", subdet.id())
    print("        type = ", subdet.type())

    # Get the hits collection name
    hitsCollection = str(sens_det.hitsCollection)
    n_cells = get_cells(subdet)
    print("        hits = ", hitsCollection)
    print("        n_cells = ", n_cells)

    # for layer_name, layer in subdet.children():
    #     layer_id = layer.id()
    #     print(f"         layer_name (id): {layer_name} ({layer_id})")
    

    # Get max dimensions for the subdetector in x,y,z (half-length form 0)
    volume = subdet.volume()
    box = volume.boundingBox()
    #mySolid=volume.solid()
    #print('dims',mySolid.dimensions())
    #print('IsCylType()',mySolid.IsCylType())
    #print("        dimensions = ", [box.x(), box.y(), box.z()])
    # dimensions seem off (ask Brieuc) - for the moment multiply the max by a factor 10
    max_z = 10*math.ceil(box.z()) #max z rounded up
    max_r = max(10*box.x(),10*box.y())#10*math.ceil(math.sqrt(math.pow(box.x(),2)+math.pow(box.y(),2))) #max r rounded up
    if(box.x()!=box.y()): print ('WARNING: different X/Y in the detector bounding box X:',box.x(),'Y:',box.y())

    if subdet.id()>0.:
      dict_sub[str(subdet_name)]={
          'id': int(subdet.id()),
          'hitsCollection': hitsCollection,
          'n_cells': n_cells,
          'max_z': max_z, 
          'max_r': max_r,
          }
    # Fill a dictionary to match via ID subdetectors (togheter with its max dimensions)to the SimHits collections
    #if subdet.id()>0.:
    #    binning_string = str(nbins_z)+','+str(-max_z)+'.,'+str(max_z)+'.,'+str(nbins_r)+','+str(-max_r)+'.,'+str(max_r)+'.'
    #    print("        binning = ", binning_string)
    #    dict_sub[str(subdet.id())] ={"binning":  binning_string}

import json
with open(geo_file_name.split('/')[-1].strip('.xml')+'_DetectorDimensions.json', 'w') as fp:
    json.dump(dict_sub, fp, indent=2)
