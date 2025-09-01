#!/usr/bin/env python

import dd4hep as dd4hepModule
from ROOT import dd4hep
from podio import root_io
from optparse import OptionParser
import math
import subprocess
from glob import glob
import os

#TODO: # fix this script

######################################
# option parser

parser = OptionParser()
parser.add_option('-d', '--detGeoFile',
                  type=str, default='$K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml',
                  help='compact XML detector geometry file')
parser.add_option('-i', '--inputFile',
                  type=str, default='/eos/home-s/sfranche/FCC/BIB/data/aciarma_4IP_2024may29/Z/DDSim_output/bib_v1/ALLEGRO_o1_v03_1_r2025-05-29.root',
                  help='input file simulated with the same geo as detGeoFile')
parser.add_option('-n', '--numberOfFiles',
                  type=int, default=100,
                  help='number of input files to be considered in drawhits.py called by this script')
parser.add_option('--nBinsR',
                  type=int, default=200,
                  help='binning in radial direction')
parser.add_option('--nBinsZ',
                  type=int, default=200,
                  help='binning in z direction')
parser.add_option('-s', '--sample',
                  type=str, default='ipc',
                  help='sample name to save plots')

(options, args) = parser.parse_args()

geo_file_name = options.detGeoFile
input_file_path = options.inputFile
nfiles = options.numberOfFiles
nbins_r = options.nBinsR
nbins_z = options.nBinsZ
sample_name = options.sample

# Prepare input
input_path = glob(input_file_path)

input_files = []
if os.path.isdir(input_path[0]):
    for p in input_path:
        input_files += glob(p+"/*.root")
elif isinstance(input_path, str):
    print("Input is a single file: ", input_path)
    input_files = [input_path]
else:
    input_files = input_path


##################################################
# Geometry part

# useful doxy:
# https://dd4hep.web.cern.ch/dd4hep/reference/classdd4hep_1_1Detector.html
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
for subdet_name in subdets:
    print("==== name = ", subdet_name[0]) 
    subdet = det.detector(subdet_name[0])
    print("        ID = ", subdet.id())
    print("        type = ", subdet.type())

    # Get max dimensions for the subdetector in x,y,z (half-length form 0)
    volume = subdet.volume()
    box = volume.boundingBox()
    print("        dimensions = ", [box.x(), box.y(), box.z()])
    # dimensions seem off (ask Brieuc) - for the moment multiply the max by a factor 10
    max_z = 10*math.ceil(box.z()) #max z rounded up
    max_r = 10*math.ceil(math.sqrt(math.pow(box.x(),2)+math.pow(box.y(),2))) #max r rounded up

    # Fill a dictionary to match via ID subdetectors (togheter with its max dimensions)to the SimHits collections
    if subdet.id()>0.:
        binning_string = str(nbins_z)+','+str(-max_z)+'.,'+str(max_z)+'.,'+str(nbins_r)+','+str(-max_r)+'.,'+str(max_r)+'.'
        print("        binning = ", binning_string)
        dict_sub[str(subdet.id())] ={"binning":  binning_string}


##################################################
# Simulation part

print("input file = ", input_files[0])
podio_reader = root_io.Reader(input_files[0])
dict_collection = {}

# Get collections from sim file metadata
metadata = podio_reader.get("metadata")[0]
collections_list = metadata.parameters
for collection_encoding in collections_list:
    collection = collection_encoding.replace("__CellIDEncoding", "")
    #print ("-- collection_encoding = ", collection_encoding)
    #print ("-- collection = ", collection)
    cellid_encoding = metadata.get_parameter(collection_encoding)
    decoder = dd4hep.BitFieldCoder(cellid_encoding)
    dict_collection[collection] = {"decoder": decoder}

black_list = ["DRcaloSiPMreadout"]
for event in podio_reader.get("events"):
    for collection_name, collection_decoder in dict_collection.items():

        # Skip collections that are (temporary) blacklisted
        if collection_name in black_list:
            print(f"NOTE: {collection_name} is blacklisted, skipping...")
            continue

        # get first hit of collection to extract the detector ID and then stop
        for hit in event.get(collection_name):
            cellID = hit.getCellID()
            subdetID = collection_decoder["decoder"].get(cellID, "system")
            #print ("-- collection = ", collection_name)
            #print ("-- subdetID = ", subdetID)
            if str(subdetID) in dict_sub:
                dict_sub[str(subdetID)]["collection"] = collection_name
            break


print(dict_sub)
# run python script to create the hit maps for each subdetector
for subdet_id, value in dict_sub.items():
    if "collection" in value:
        # if (dict_sub[subdet_id]["collection"] == "VertexBarrelCollection"):
        print(" plots for subdet id = ", subdet_id)
        subprocess.run([
            "drawhits.py",
            "--infilePath="+input_file_path,
            "--collection="+dict_sub[subdet_id]["collection"],
            "--numberOfFiles="+str(nfiles),
            "--sample=ipc_"+detector,
            "-m",
            "--map_binning="+dict_sub[subdet_id]["binning"],
        ])
