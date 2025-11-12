#!/usr/bin/env python

import argparse
from collections import defaultdict
import re

import ROOT

from helpers import load_json, simplify_dict, layer_number_from_string, is_endcap
from constants import b_to_GB, MHz_to_Hz, cm2_to_mm2
from visualization import setup_root_style, draw_hist


######################################
# option parser

parser = argparse.ArgumentParser(description='hits2highLevelEstimations.py',
        epilog='Example:\nhits2highLevelEstimations.py -i <path_to_hits_histograms.root> -d $BIB_STUDIES/detectors_dicts/ALLEGRO_o1_v03_DetectorDimensions.json -a $BIB_STUDIES/detectors_dicts/ALLEGRO_o1_v03_assumptions.json',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--inputFile',
                  type=str, default='',
                  help='path to input file containing hits histograms.')
parser.add_argument('-o', '--outputFile',
                  type=str, default='high-level_estimates',
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

detector_type = detector_dict["typeFlag"]

hits_collection = detector_dict["hitsCollection"]
strategy = assumptions_dict["strategy"]
hit_size = assumptions_dict["hit_size"]
multipliers = assumptions_dict["multipliers"]

pixel_size_u = simplify_dict(assumptions_dict["pixel_size_u"])
pixel_size_v = simplify_dict(assumptions_dict["pixel_size_v"])
print(f"Pixel size u: {pixel_size_u}, pixel size v: {pixel_size_v}")

# Update layer related dictionary to have identical keys
n_cells = simplify_dict(detector_dict["det_element_cells"])
print("Number of cells: ",n_cells)
sensor_size_map = simplify_dict(detector_dict["sensor_size_map"])
sensors_per_module_map = simplify_dict(detector_dict["sensors_per_module_map"])
print(f"Sensor size: {sensor_size_map}, sensors per module: {sensors_per_module_map}")

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
h_bw.SetNameTitle(f"{input_file_name}_bw_per_layer",f"{input_file_name}_bw_per_layer")
h_avg_hit_rate = input_file.Get(h_name).Clone()
h_max_hit_rate = input_file.Get(h_name).Clone() 
h_max_hit_rate.Reset()
h_avg_occ_cell_per_layer = input_file.Get(h_name).Clone()
h_max_cell_occ = input_file.Get(h_name).Clone() 
h_max_cell_occ.Reset()

layer_cells = {}
for ln, cells in detector_dict["det_element_cells"].items():
    layer_cells[ln] = cells

n_layers = len(layer_cells.keys())
layer_binning = [n_layers, -0.5, n_layers - 0.5]
if is_endcap(detector_type):
    print("Endcap detector detected, adjusting layer binning accordingly")
    max_l = int(n_layers / 2) + 0.5
    layer_binning = [n_layers + 1, -max_l, +max_l]


hist_n_cells = ROOT.TH1D("hist_n_cells", "Number of Cells per Layer;Layer;Number of Cells", *layer_binning) # This is the number of sensors in case of semiconductor detectors
hist_sensor_size = ROOT.TH1D("hist_sensor_size", "Sensor size per Layer;Layer;Sensor size [mm^{2}]", *layer_binning)
hist_sensors_per_module = ROOT.TH1D("hist_sensors_per_module", "Sensors per module;Layer;Sensors per module", *layer_binning)
hist_pixel_size_u = ROOT.TH1D("hist_pixel_size_u", "Pixel Size U per Layer;Layer;Pixel Size U [mm]", *layer_binning)
hist_pixel_size_v = ROOT.TH1D("hist_pixel_size_v", "Pixel Size V per Layer;Layer;Pixel Size V [mm]", *layer_binning)

# Fill the histograms with data
for layer, value in n_cells.items():
    print(int(layer+1+(len(n_cells)/2 if is_endcap(detector_type) else 0)), value)
    hist_n_cells.SetBinContent(int(layer+1+(len(n_cells)/2 if is_endcap(detector_type) else 0)), value)

for layer, value in sensor_size_map.items():
    hist_sensor_size.SetBinContent(int(layer+1+(len(n_cells)/2 if is_endcap(detector_type) else 0)), value)

for layer, value in sensors_per_module_map.items():
    hist_sensors_per_module.SetBinContent(int(layer+1+(len(n_cells)/2 if is_endcap(detector_type) else 0)), value)

# Currently not used but kept for future use
hist_module_size = hist_sensor_size.Clone()
hist_module_size.Multiply(hist_sensors_per_module)
hist_module_size.SetNameTitle("hist_module_size", "Module size per Layer;Layer;Module size [mm^{2}]")

for layer, value in pixel_size_u.items():
    hist_pixel_size_u.SetBinContent(int(layer+1+(len(n_cells)/2 if is_endcap(detector_type) else 0)), value)

for layer, value in pixel_size_v.items():
    hist_pixel_size_v.SetBinContent(int(layer+1+(len(n_cells)/2 if is_endcap(detector_type) else 0)), value)


#######################################
# convert the counts/occupancy histogram to estimate of high-level properties (bandwidth, hit rate, etc.)

# Bandwidth
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
        scale_factor *= cells[layer_n] * 0.01

    # Consider additional modifiers
    for m in multipliers.values():
        scale_factor *= m

    h_bw.SetBinContent(b, counts * scale_factor)
    h_bw.SetBinError(b, error * scale_factor)

scale_factor = 1
for m in multipliers.values():
    scale_factor *= m

# Average hit rate per layer
h_avg_hit_rate.Scale(rate*cm2_to_mm2*scale_factor)
h_avg_hit_rate.Divide(hist_n_cells*hist_sensor_size)
h_avg_hit_rate.SetNameTitle(f"{input_file_name}_avg_hit_rate_per_layer", f"{input_file_name}_avg_hit_rate_per_layer")
draw_hist(h_avg_hit_rate, "Layer", "Average hit rate [MHz/cm^{2}]", f"{input_file_name}_hit_rate_per_layer")

# Average cell occupancy per layer
h_avg_occ_cell_per_layer.Scale(scale_factor)
h_avg_occ_cell_per_layer.Divide(hist_n_cells*hist_sensor_size/hist_pixel_size_u/hist_pixel_size_v)
h_avg_occ_cell_per_layer.SetNameTitle(f"{input_file_name}_avg_occ_per_layer", f"{input_file_name}_avg_occ_per_layer;Layer;Average pixel occupancy per event")
draw_hist(h_avg_occ_cell_per_layer, "Layer", "Average pixel occupancy per event", f"{input_file_name}_avg_pixel_occupancy_per_layer")



h_avg_hit_rate_per_cell = {}
h_occ_per_cell = {}

for i, (ln, cells) in enumerate(detector_dict["det_element_cells"].items()):
    if is_endcap(detector_type):
        i_layer_bin = int(ln + len(detector_dict["det_element_cells"])/2) + 1 # to skip layer 0 in case of disk
    else:
        i_layer_bin = ln + 1

    # Hit rate per cell (i.e. module in semiconductor detectors)
    h_avg_hit_rate_per_cell[ln] = input_file.Get(f"h_avg_hits_x_layer{ln}_per_cell_{hits_collection}").Clone()
    h_avg_hit_rate_per_cell[ln].Scale(rate*cm2_to_mm2*scale_factor/hist_module_size.GetBinContent(i_layer_bin))
    h_avg_hit_rate_per_cell[ln].SetNameTitle(f"{input_file_name}_hitRate_layer{ln}_per_cell", f"{input_file_name}_hitRate_layer{ln}_per_cell;Module;Average hit rate per module [MHz/cm^{2}]" )
    draw_hist(h_avg_hit_rate_per_cell[ln], "Module", "Average hit rate [MHz/cm^{2}]", f"{input_file_name}_hitRate_layer{ln}_per_cell", )

    # Extract maximal hit rate per cell
    h_max_hit_rate.SetBinContent(i_layer_bin, h_avg_hit_rate_per_cell[ln].GetMaximum())
    h_max_hit_rate.SetBinError(i_layer_bin, h_avg_hit_rate_per_cell[ln].GetBinError(h_avg_hit_rate_per_cell[ln].GetMaximumBin()))

    # Occupancy per cell (i.e. pixel occupancy in semiconductor detector)
    h_occ_per_cell[ln] = input_file.Get(f"h_avg_hits_x_layer{ln}_per_cell_{hits_collection}").Clone()
    h_occ_per_cell[ln].Scale(scale_factor/(hist_module_size.GetBinContent(i_layer_bin)/hist_pixel_size_u.GetBinContent(i_layer_bin)/hist_pixel_size_v.GetBinContent(i_layer_bin)))
    h_occ_per_cell[ln].SetNameTitle(f"{input_file_name}_occ_per_cell_layer{ln}", f"{input_file_name}_occ_per_cell_layer{ln};Module;Pixel occupancy per event" )
    draw_hist(h_occ_per_cell[ln], "Module", "Pixel occupancy per event", f"{input_file_name}_occ_per_cell_layer{ln}")

    # Extract maximal occupancy per cell
    h_max_cell_occ.SetBinContent(i_layer_bin, h_occ_per_cell[ln].GetMaximum())
    h_max_cell_occ.SetBinError(i_layer_bin, h_occ_per_cell[ln].GetBinError(h_occ_per_cell[ln].GetMaximumBin()))

h_max_hit_rate.SetNameTitle(f"{input_file_name}_max_hit_rate_per_cell", f"{input_file_name}_max_hit_rate_per_cell;Layer;Maximal hit rate per module [MHz/cm^{2}]")
draw_hist(h_max_hit_rate, "Layer", "Maximal hit rate [MHz/cm^{2}]", f"{input_file_name}_max_hit_rate_per_layer")

h_max_cell_occ.SetNameTitle(f"{input_file_name}_max_cell_occupancy", f"{input_file_name}_max_cell_occupancy;Layer;Maximal pixel occupancy per event")
draw_hist(h_max_cell_occ, "Layer", "Maximal pixel occupancy per event", f"{input_file_name}_max_pixel_occupancy_per_layer")

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
    h_avg_hit_rate.Write()
    h_max_hit_rate.Write()
    h_avg_occ_cell_per_layer.Write()
    h_max_cell_occ.Write()
    hist_n_cells.Write()
    hist_sensor_size.Write()
    hist_sensors_per_module.Write()
    hist_module_size.Write()
    hist_pixel_size_u.Write()
    hist_pixel_size_v.Write()
    for ln, cells in detector_dict["det_element_cells"].items():
        h_avg_hit_rate_per_cell[ln].Write()
        h_occ_per_cell[ln].Write()
print("Histograms saved in:", output_file_name)