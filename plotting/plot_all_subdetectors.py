#!/usr/bin/env python

import dd4hep as dd4hepModule
from ROOT import dd4hep
from podio import root_io
import argparse
import subprocess
from glob import glob
import os
import time

from helpers import DetFilePath, print_header

######################################
# option parser

def parse_args():
    parser = argparse.ArgumentParser(description='plot_all_subdetectors.py', 
        epilog='Example:\ndrawhits.py -s MuonTaggerBarrel -d $BIB_STUDIES/detectors_dicts/ALLEGRO_o1_v03_DetectorDimensions.json -e -1 -D 0',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-D', '--debugLevel',
                    type=int, default=1,
                    help='debug level (0:quiet, 1:default, 2:verbose)')
    parser.add_argument('-d', '--detGeoFile',
                      type=str, default='$K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml',
                      help='compact XML detector geometry file')
    parser.add_argument('-i', '--inputFile',
                      type=str, default='/eos/home-s/sfranche/FCC/BIB/data/aciarma_4IP_2024may29/Z/DDSim_output/bib_v1/ALLEGRO_o1_v03_1_r2025-05-29.root',
                      help='input file simulated with the **same** geo as detGeoFile')
    parser.add_argument('-n', '--numberOfFiles',
                      type=int, default=100,
                      help='number of input files to be considered in drawhits.py called by this script')
    parser.add_argument('-e', '--numberOfEvents',
                    type=int, default=1,
                    help='number of events per file to consider, put -1 to take all events')
    # parser.add_argument('--nBinsR',
    #                   type=int, default=200,
    #                   help='binning in radial direction')
    # parser.add_argument('--nBinsZ',
    #                   type=int, default=200,
    #                   help='binning in z direction')
    parser.add_argument('-s', '--sample',
                      type=str, default='ipc',
                      help='sample name to save plots')

    return (parser.parse_args())

options = parse_args()

debug_level = options.debugLevel
geo_file_name = options.detGeoFile
det_file =  DetFilePath(options.detGeoFile)
input_file_path = options.inputFile
nfiles = options.numberOfFiles
nevts_per_file = options.numberOfEvents
# nbins_r = options.nBinsR
# nbins_z = options.nBinsZ
sample_name = options.sample

detector = det_file.name #includes version
assumptions_file_name = "$BIB_STUDIES/detectors_dicts/"+detector+"_assumptions.json"
#read assumptions file
import json
with open(os.path.expandvars(assumptions_file_name)) as f:
    assumptions = json.load(f)


# # Prepare input
# input_path = glob(input_file_path)

# input_files = []
# if os.path.isdir(input_path[0]):
#     for p in input_path:
#         input_files += glob(p+"/*.root")
# elif isinstance(input_path, str):
#     print("Input is a single file: ", input_path)
#     input_files = [input_path]
# else:
#     input_files = input_path


##################################################
# Geometry part

# useful doxy:
# https://dd4hep.web.cern.ch/dd4hep/reference/classdd4hep_1_1Detector.html
# https://dd4hep.web.cern.ch/dd4hep/reference/classdd4hep_1_1DetElement.html
# https://dd4hep.web.cern.ch/dd4hep/reference/classdd4hep_1_1Volume.html
# https://dd4hep.web.cern.ch/dd4hep/reference/classdd4hep_1_1Box.html
# https://dd4hep.web.cern.ch/dd4hep/reference/DD4hep_2Detector_8h_source.html

# Get detector name
# detector = geo_file_name.split("/")[-1].strip(".xml")

# # Load geometry
# det = dd4hep.Detector.getInstance()
# det.fromXML(geo_file_name)

# # Get subdetectors
# dict_sub = {}
# #subdet = det.detector("VertexBarrel")
# #subdets = det.detectors();
# subdets = det.sensitiveDetectors();
# for subdet_name in subdets:
#     print("==== name = ", subdet_name[0]) 
#     subdet = det.detector(subdet_name[0])
#     print("        ID = ", subdet.id())
#     print("        type = ", subdet.type())

#     # Get max dimensions for the subdetector in x,y,z (half-length form 0)
#     volume = subdet.volume()
#     box = volume.boundingBox()
#     print("        dimensions = ", [box.x(), box.y(), box.z()])
#     # dimensions seem off (ask Brieuc) - for the moment multiply the max by a factor 10
#     max_z = 10*math.ceil(box.z()) #max z rounded up
#     max_r = 10*math.ceil(math.sqrt(math.pow(box.x(),2)+math.pow(box.y(),2))) #max r rounded up

#     # Fill a dictionary to match via ID subdetectors (togheter with its max dimensions)to the SimHits collections
#     if subdet.id()>0.:
#         binning_string = str(nbins_z)+','+str(-max_z)+'.,'+str(max_z)+'.,'+str(nbins_r)+','+str(-max_r)+'.,'+str(max_r)+'.'
#         print("        binning = ", binning_string)
#         dict_sub[str(subdet.id())] ={"binning":  binning_string}


##################################################
# Simulation part

# print("input file = ", input_files[0])
# podio_reader = root_io.Reader(input_files[0])
# dict_collection = {}

# # Get collections from sim file metadata
# metadata = podio_reader.get("metadata")[0]
# collections_list = metadata.parameters
# for collection_encoding in collections_list:
#     collection = collection_encoding.replace("__CellIDEncoding", "")
#     cellid_encoding = metadata.get_parameter(collection_encoding)
#     decoder = dd4hep.BitFieldCoder(cellid_encoding)
#     dict_collection[collection] = {"decoder": decoder}

# black_list = ["DRcaloSiPMreadout"]
# for event in podio_reader.get("events"):
#     for collection_name, collection_decoder in dict_collection.items():

#         # Skip collections that are (temporary) blacklisted
#         if collection_name in black_list:
#             print(f"NOTE: {collection_name} is blacklisted, skipping...")
#             continue

#         # get first hit of collection to extract the detector ID and then stop
#         for hit in event.get(collection_name):
#             cellID = hit.getCellID()
#             subdetID = collection_decoder["decoder"].get(cellID, "system")
#             #print ("-- collection = ", collection_name)
#             #print ("-- subdetID = ", subdetID)
#             if str(subdetID) in dict_sub:
#                 dict_sub[str(subdetID)]["collection"] = collection_name
#             break


#print(dict_sub)

# run python script to create the hit maps for each subdetector
for subdet_id,dummy_values in assumptions.items():
        # if (dict_sub[subdet_id]["collection"] == "VertexBarrelCollection"):
        print_header(f"Plotting subdetector ID: {subdet_id}")
        
        # Start timing for this subdetector
        start_time = time.time()
        
        subprocess.run([
            "drawhits.py",
            "--debugLevel="+str(debug_level),
            "--infilePath="+input_file_path,
            "--subDetector="+subdet_id,
            #"--collection="+dict_sub[subdet_id]["collection"],
            "--numberOfFiles="+str(nfiles),
            "--numberOfEvents="+str(nevts_per_file),
            "--sample="+sample_name+"_"+detector,
            #"-m",
            #"--map_binning="+dict_sub[subdet_id]["binning"],
        ])
        
        # Calculate and print elapsed time for this subdetector
        elapsed_time = time.time() - start_time
        print(f"Time for {subdet_id}: {elapsed_time:.2f} seconds")
