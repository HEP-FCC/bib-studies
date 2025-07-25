#!/usr/bin/env python

from collections import Counter
from glob import glob
import json
import math
import os
from optparse import OptionParser
import re

import dd4hep as dd4hepModule
from podio import root_io
import ROOT



######################################
# style
ROOT.gStyle.SetOptStat(0000)
ROOT.gStyle.SetOptFit(0000)

######################################
# option parser

parser = OptionParser()
parser.add_option('-i', '--infilePath',
                  type=str, default='/eos/home-s/sfranche/FCC/BIB/data/aciarma_4IP_2024may29/Z/DDSim_output/bib_v1/',
                  help='path to directory with input files')
parser.add_option('-n', '--numberOfFiles',
                  type=int, default=-1,
                  help='number of IPC files to consider (1 event per file), put -1 to take all files')
parser.add_option('-t', '--tree',
                  type=str, default='events',
                  help='name of the tree in the root file')
parser.add_option('--sample',
                  type=str, default='ipc',
                  help='sample name to save plots')
parser.add_option('-d', '--detDictFile',
                  type=str, default="",
                  help='JSON dictionary with some key detector parameters')
parser.add_option('-s', '--subDetector',
                  type=str, default='VertexBarrel',
                  help='variable against with to draw the plots, options are: eta, phi, pt or mu')
parser.add_option('-p', '--draw_profiles',
                  action="store_true",
                  help='activate drawing of profile plots')
parser.add_option('-m', '--draw_maps',
                  action="store_true",
                  help='activate drawing of number of maps plots')
parser.add_option('--map_binning',
                  type=str, default='200,-200., 200., 200, -50., 50.',
                  help='comma separated binning for map plots (2D plots): nx,xmin,xmax,ny,ymin,ymax')

(options, args) = parser.parse_args()

nfiles = options.numberOfFiles
path = options.infilePath
tree_name = options.tree
sample_name = options.sample
detector_dict_path = options.detDictFile
sub_detector = options.subDetector
draw_maps = options.draw_maps
draw_profiles = options.draw_profiles
binning = []
for i,x in enumerate(options.map_binning.split(',')):
    if i==0 or i==3:
        binning.append(int(x))
    else:
        binning.append(float(x))
#binning = [float(x) for x in options.map_binning.split(',')]


#######################################
# functions

def myText(x, y, color, text, size=0.04, angle=0.):
    l = ROOT.TLatex()
    l.SetNDC()
    l.SetTextColor(color)
    l.PaintLatex( x, y, angle, size, text)
    l.DrawLatex(x,y,text)


def draw_map(hist, titlex, titley, plot_title, message=''):

    cm1 = ROOT.TCanvas("cm1_"+message, "cm1_"+message, 800, 600)

    hist.GetXaxis().SetTitle(titlex)
    hist.GetYaxis().SetTitle(titley)
    hist.Draw("COLZ")
    entries = hist.GetEntries()

    myText(0.20, 0.9, ROOT.kBlack, message + '  , entries = ' + str(entries))

    cm1.Print(plot_title+".pdf")
    cm1.Print(plot_title+".png")


def draw_profile(hist, titlex, titley, plot_title, message=''):

    cm1 = ROOT.TCanvas("cm1_"+message, "cm1_"+message, 800, 600)

    hist.GetXaxis().SetTitle(titlex)
    hist.GetYaxis().SetTitle(titley)
    hist.Draw("HistE")
    entries = hist.GetEntries()

    myText(0.20, 0.9, ROOT.kBlack, message + '  , entries = ' + str(entries))

    cm1.Print(plot_title+".pdf")
    cm1.Print(plot_title+".png")


#######################################
# parse the input path

input_path = glob(path)

list_files = []
if os.path.isdir(input_path[0]):
    for p in input_path:
        list_files += glob(p+"/*.root")
elif isinstance(input_path, str):
    list_files = [input_path]
else:
    list_files = input_path

# list nfiles in the path directory ending in .root  to consider for plotting
list_input_files = []
for f in sorted(list_files): #use sorted list to make sure to always take the same files
    if f.endswith(".root"):
        list_input_files.append( os.path.join(path, f) )
    if nfiles > 0 and len(list_input_files) >= nfiles: #if nfiles=-1 run all files
        break

# Load the detector json dictionary
detector_dict = {}
with open(detector_dict_path,"r") as f:
    detector_dict = json.load(f)[sub_detector]


# Count  the number of cells per layer, merging different sides
# leveraging the current naming scheme of active layers
# TODO: this can use instead directly the layer IDs once they're better understood
layer_cells = Counter()
for de, cells in detector_dict["det_element_cells"].items():
    try:
        layer_name = re.search("layer(_)?(-)?[0-9]+",de).group(0)
    except AttributeError:
        print("Warning: Couldn't read layer number from", de)
        continue

    # Get the layer number, ignore the side (sign -)
    ln = int(layer_name.replace("layer","").replace("_","").replace("-",""))
    
    layer_cells[ln] += cells


#get_layer_cells(detector_dict["det_element_cells"])
n_layers = len(layer_cells.keys())

#######################################
# prepare histograms

collection = detector_dict["hitsCollection"]
z_range = int(detector_dict["max_z"]*1.1)
r_range = int(detector_dict["max_r"]*1.1)
n_bins_z = int(z_range * 2)  # 1mm binning  
n_bins_r = int(r_range * 2)  # 1mm binning  

# Profile histograms
h_hit_E = ROOT.TH1D("h_hit_E_"+collection, "h_hit_E_"+collection, 100, 0, 1)
h_particle_E = ROOT.TH1D("h_particle_E_"+collection, "h_particle_E_"+collection, 100, 0, 1)
h_hits_x_layer = ROOT.TH1D("h_hits_x_layer_"+collection, "h_hits_x_layer_"+collection, n_layers, -0.5, n_layers-0.5)
h_occ_x_layer = ROOT.TH1D("h_occ_x_layer_"+collection, "h_hits_x_layer_"+collection, n_layers, -0.5, n_layers-0.5)

# Hit maps definitions
hist_zr = ROOT.TH2D("hist_zr_"+collection, "hist_zr_"+collection, n_bins_z, -z_range, z_range, n_bins_r, -r_range, -r_range)
hist_xy = ROOT.TH2D("hist_xy_"+collection, "hist_xy_"+collection, n_bins_r, -r_range, -r_range, n_bins_r, -r_range, -r_range)
hist_zphi = ROOT.TH2D("hist_zphi_"+collection, "hist_zphi_"+collection, n_bins_z, -z_range, z_range, 500, -3.5, 3.5)

fill_weight = 1. / len(list_input_files)

#######################################
# run the event loop

podio_reader = root_io.Reader(list_input_files)

# Check if collection is valid and setup a cell ID decoder
metadata = podio_reader.get("metadata")[0]

id_encoding = metadata.get_parameter(collection+"__CellIDEncoding")
decoder = ROOT.dd4hep.BitFieldCoder(id_encoding)

for i,event in enumerate(podio_reader.get(tree_name)):
    
    if i % 100 == 0: # print a message every 100 processed files
        print('processing file number = ' + str(i) + ' / ' + str(nfiles) )

    fired_cells = []
    fired_cells_x_layer = Counter()

    # Loop over the hits
    for hit in event.get(collection):
        # Hits doxy:
        # https://edm4hep.web.cern.ch/classedm4hep_1_1_mutable_sim_tracker_hit.html
        # https://edm4hep.web.cern.ch/classedm4hep_1_1_mutable_sim_calorimeter_hit.html
        # https://edm4hep.web.cern.ch/classedm4hep_1_1_m_c_particle.html

        is_calo_hit = "CalorimeterHit" in str(hit)

        cell_id = hit.cellID()

        # Skip cells firing multiple times
        if cell_id in fired_cells:
            # print(cell_id,"has already fired in this event!")
            continue
        else:
            fired_cells.append(cell_id)

        # side = decoder.get(cell_id, "side")
        layer = decoder.get(cell_id, "layer")
        #slayer = decoder.get(cell_id, "superlayer")
        # print(f"side = {side}, layer = {layer}")
        fired_cells_x_layer[layer] += 1

        x = hit.getPosition().x
        y = hit.getPosition().y
        z = hit.getPosition().z
        r = math.sqrt(math.pow(x, 2) + math.pow(y, 2))
        phi = math.acos(x/r) * math.copysign(1, y)

        if x < 0.:
            r = -1. * r

        # Deposited energy, converted to MeV
        E_hit = 0
        if is_calo_hit:
            E_hit = hit.getEnergy() * 1e3  # For calorimeter hits
        else:
            E_hit = hit.getEDep() * 1e3  # For tracker hits

        # Fill the hits histograms
        h_hit_E.Fill(E_hit, fill_weight)
        h_hits_x_layer.Fill(layer, fill_weight)

        hist_zr.Fill(z, r, fill_weight)
        hist_xy.Fill(x, y, fill_weight)
        hist_zphi.Fill(z, phi, fill_weight)

        if not is_calo_hit:
            particle = hit.getParticle()
            E_particle = particle.getEnergy()

            # Fill MC particle histograms
            h_particle_E.Fill(E_particle, fill_weight)

    # <--- end of the hits loop

    print
    for l, fc in fired_cells_x_layer.items():
        #print(f" fc:{fc}, l:{l},  nc:{layer_cells[l]} ")
        layer_occupancy = fc / layer_cells[l] * 100
        #print(f" fc:{fc}, nc:{layer_cells[l]}  , occ:{layer_occupancy}")

        h_occ_x_layer.Fill(l, layer_occupancy * fill_weight)
# <--- end of the event loop


if draw_maps:
    draw_map(hist_zr, "z [mm]", "r [mm]", sample_name+"_map_zr_"+str(nfiles)+"evt_"+sub_detector, collection)
    draw_map(hist_xy, "x [mm]", "y [mm]", sample_name+"_map_xy_"+str(nfiles)+"evt_"+sub_detector, collection)
    draw_map(hist_zphi, "z [mm]", "#phi [rad]", sample_name+"_map_zphi_"+str(nfiles)+"evt_"+sub_detector, collection)

if draw_profiles: 
    draw_profile(h_hit_E, "Deposited energy [MeV]", "Hits / events",  sample_name+"_hit_E_"+str(nfiles)+"evt_"+sub_detector, collection)
    draw_profile(h_hits_x_layer, "Layer number", "Hits / events",  sample_name+"_hits_x_layer_"+str(nfiles)+"evt_"+sub_detector, collection)
    draw_profile(h_occ_x_layer, "Layer number", "Occupancy [%] / events",  sample_name+"_occ_x_layer_"+str(nfiles)+"evt_"+sub_detector, collection)
    
    draw_profile(h_particle_E, "MC particle energy [GeV]", "Particle / events",  sample_name+"_particle_E_"+str(nfiles)+"evt_"+sub_detector, collection)