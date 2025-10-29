#!/usr/bin/env python

import argparse
from collections import defaultdict
import re

import ROOT

from helpers import load_json, simplify_dict
from constants import b_to_GB, MHz_to_Hz
from visualization import setup_root_style, draw_hist


######################################
# option parser

parser = argparse.ArgumentParser(description='hits2bw.py',
        epilog='Example:\nhits2bw.py -i <path_to_hits_histograms.root> -d $BIB_STUDIES/detectors_dicts/ALLEGRO_o1_v03_DetectorDimensions.json -a $BIB_STUDIES/detectors_dicts/ALLEGRO_o1_v03_assumptions.json',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--inputFile',
                  type=str, default='',
                  help='path to input file containing hits histograms.')
parser.add_argument('-o', '--outputFile',
                  type=str, default='bandwidths',
                  help='Name of output root file, the sample name will be added as prefix and the detector as suffix.')
parser.add_argument('-d', '--detDictFile',
                  type=str, default='$BIB_STUDIES/detectors_dicts/ALLEGRO_o1_v03_DetectorDimensions.json',
                  help='JSON dictionary with some key detector parameters.')
parser.add_argument('-a', '--assumptions',
                  type=str, default='$BIB_STUDIES/detectors_dicts/ALLEGRO_o1_v03_assumptions.json',
                  help='JSON dictionary with assumptions for bandwidth estimates.')
parser.add_argument('-r', '--rate',
                  type=float, default=52.,
                  help='hit rate in MHz.')

options = parser.parse_args()

input_file_path = options.inputFile
output_file_name = options.outputFile
detector_dict_path = options.detDictFile
assumptions_path = options.assumptions
rate = options.rate

######################################
# style

setup_root_style(stat_box=False)
#to avoid canvas showing up slowing us down
ROOT.gROOT.SetBatch(True)
#to avoid canvas->Print printouts
ROOT.gErrorIgnoreLevel = ROOT.kWarning


#######################################
# functions


#######################################
# parse the input path and/or files

input_file_name = input_file_path.split("/")[-1].strip(".root")
sub_detector = re.sub(".*[0-9]+(evt_)","",input_file_name) #.group(0)

print(f"Reading file '{input_file_name}' (sub detector: {sub_detector})")

input_file = ROOT.TFile(input_file_path, "READ")
detector_dict = load_json(detector_dict_path, sub_detector)
assumptions_dict = load_json(assumptions_path, sub_detector)

hits_collection = detector_dict["hitsCollection"]
strategy = assumptions_dict["strategy"]
hit_size = assumptions_dict["hit_size"]
multipliers = assumptions_dict["multipliers"]

# Update layer related dictionary to have identical keys
channels = simplify_dict(detector_dict["det_element_cells"])
if isinstance(hit_size, dict):
    hit_size_tmp = simplify_dict(hit_size)

    # Convert to defaultdict to avoid KeyErrors
    # for extra layers not defined in the assumptions
    # (e.g. the 0 or the N+1 layers)
    hit_size = defaultdict(int)
    hit_size.update(hit_size_tmp)

# Pick the right histogram depending on the defined strategy
h_name = ""
match strategy:
    case "hit_counts":
        h_name = f"h_avg_hits_x_layer_{hits_collection}"
    case "occupancy":
        h_name = f"h_avg_occ_x_layer_{hits_collection}"
    case _:
        raise AttributeError(f"Unknown strategy defined to compute the bandwidth ({strategy})")

h_bw = input_file.Get(h_name).Clone()


#######################################
# convert the counts/occupancy histogram to BW

for b in range(1, h_bw.GetNbinsX()+1):
    counts = h_bw.GetBinContent(b)
    error = h_bw.GetBinError(b)
    layer_n = h_bw.GetBinCenter(b)

    # convert hits / occupancy to GB / s
    scale_factor = 1
    try:
        scale_factor *= rate * MHz_to_Hz * hit_size * b_to_GB
    except TypeError:
        scale_factor *= rate * MHz_to_Hz * hit_size[layer_n] * b_to_GB

    if strategy == "occupancy":
        scale_factor *= channels[layer_n] * 0.01

    # Consider additional modifiers
    for m in multipliers.values():
        scale_factor *= m

    h_bw.SetBinContent(b, counts * scale_factor)
    h_bw.SetBinError(b, error * scale_factor)

#######################################
# Output the results

tot_bw = h_bw.Integral()
tot_bw_msg = f"Total bandwidth = {tot_bw:.2f} GB/s"
print(tot_bw_msg)

draw_hist(h_bw, "Layer", "Bandwidth [GB/s]", f"{input_file_name}_bw_per_layer", tot_bw_msg)

# Write the histograms to the output file
output_file_name = f"{input_file_name}_{output_file_name}.root"
with ROOT.TFile(output_file_name,"RECREATE") as f:
    h_bw.Write()

print("Histograms saved in:", output_file_name)