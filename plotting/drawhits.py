#!/usr/bin/env python

from collections import Counter, defaultdict
import json
import math
import numpy as np
from optparse import OptionParser
import re

import dd4hep as dd4hepModule
from ROOT import dd4hep
from podio import root_io
import ROOT

from constants import C_MM_NS
from helpers import path_to_list, sorted_n_files
from visualization import setup_root_style, draw_hist, draw_map



######################################
# style

setup_root_style()
#to avoid canvas showing up slowing us down
ROOT.gROOT.SetBatch(True)
#to avoid canvas->Print printouts
ROOT.gErrorIgnoreLevel = ROOT.kWarning

######################################
# option parser

parser = OptionParser()
parser.add_option('-i', '--infilePath',
                  type=str, default='/eos/home-s/sfranche/FCC/samples/bib/ipc/aciarma_4IP_2024may29_Z/CADbp_ALLEGRO_o1_v03_r2025-05-29_3998.root',
                  help='path to input file or input directory')
parser.add_option('-o', '--outputFile',
                  type=str, default='hits',
                  help='Name of output root file, the sample name will be added as prefix and the detector as suffix')
parser.add_option('-n', '--numberOfFiles',
                  type=int, default=1,
                  help='number of files to consider (1 event per file), put -1 to take all files')
parser.add_option('-e', '--numberOfEvents',
                  type=int, default=1,
                  help='number of events per file to consider, put -1 to take all events')
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
parser.add_option('-p', '--draw_hists',
                  action="store_true",
                  help='activate drawing of profile plots')
parser.add_option('-m', '--draw_maps',
                  action="store_true",
                  help='activate drawing of number of maps plots')
parser.add_option('-r','--bin_width_r',
                  type=int, default=1,
                  help='bin width for radius r in mm (default 1 mm)')
parser.add_option('-z','--bin_width_z',
                  type=int, default=1,
                  help='bin width for z in mm (default 1 mm)')
parser.add_option('-a', '--bin_width_phi',
                  type=float, default=0.01,
                  help='bin width for azimuthal angle phi')

(options, args) = parser.parse_args()

n_files = options.numberOfFiles
events_per_file = options.numberOfEvents
path = options.infilePath
output_file_name = options.outputFile
tree_name = options.tree
sample_name = options.sample
detector_dict_path = options.detDictFile
sub_detector = options.subDetector
draw_maps = options.draw_maps
draw_hists = options.draw_hists
bw_r = options.bin_width_r
bw_z = options.bin_width_z
bw_phi = options.bin_width_phi


#######################################
# functions

# Detector types from:
# https://github.com/AIDASoft/DD4hep/blob/master/DDCore/include/DD4hep/DetType.h
is_calo = lambda x: (x & dd4hep.DetType.CALORIMETER) == dd4hep.DetType.CALORIMETER
is_endcap = lambda x: (x & dd4hep.DetType.ENDCAP) == dd4hep.DetType.ENDCAP


def get_layer(cell_id, decoder, detector, dtype):
    """
    Run the decoder differently for each detector
    """

    #print(f"layer={layer}, side={side}")

    match detector:
        case "HCalThreePartsEndcap":
            # Get only the side of the hit
            # type = 0,1,2 for positive side
            # type = 3,4,5 for negative side
            tp =  decoder.get(cell_id,"type") #FIXME: Appears to be always 0?! double check using nightly sim. events
            #print(f"layer={layer}, side={side}, type={tp}")
            if tp > 2:
                return -1
            return 1

        case "EMEC_turbine":
            side = decoder.get(cell_id, "side")
            wheel = decoder.get(cell_id, "wheel") + 1
            return  wheel * side

        case "ECalBarrel": 
            return  0

        case "DCH_v2":
            # Number of layers per super layer could be read from geo file
            nl_x_sl = 8
            layer = decoder.get(cell_id, "layer")
            super_layer = decoder.get(cell_id, "superlayer")
            return (super_layer * nl_x_sl) + layer + 1

        case _:
            # Default way: side * layer, where side should be +/- 1
            layer = decoder.get(cell_id, "layer")

            # Get the side, if available
            side = 0
            if is_endcap(dtype):
                print("endcap!")
                side = decoder.get(cell_id, "side")

            if side != 0:
                layer *= side

            return layer


#######################################
# parse the input path and/or files


# list n_files in the path directory ending in .root to consider for plotting
# use sorted list to make sure to always take the same files
list_files = path_to_list(path)
list_input_files = sorted_n_files(list_files, n_files)

# Load the detector json dictionary
detector_dict = {}
with open(detector_dict_path,"r") as f:
    _dict = json.load(f)
    try:
        detector_dict = _dict[sub_detector]
    except KeyError:
        raise KeyError(f"'{sub_detector}' is not available, valid sub-detector entries are: "
                       +" | ".join(_dict.keys()))

# Parse the layer name in to an integer number
print("Retrieve the number of cells per layer")
layer_cells = {}
n_tot_cells = 0
for de, cells in detector_dict["det_element_cells"].items():
    #for ECalBArrel take the cell from bath - TODO: improve code clarity if possible
    if de == 'bath':
        layer_cells[0] = cells
        n_tot_cells += cells
    #do detector with actual layers or wheels
    side = 0
    try:
        side_name = re.search("(side)(_)?(-)?[0-9]+",de).group(0)
        side = int(re.sub(r"[^0-9-]","",side_name))
        if abs(side) > 1:
            raise NotImplementedError(f"Found side = {side} for '{de}', not supported")
    except AttributeError:
        pass

    layer_name = "000"
    try:
        layer_name = re.search("(layer|wheel)(_)?(-)?[0-9]+",de).group(0) 
    except AttributeError:
        print("Warning: Couldn't read layer number from", de)
        continue

    
    # Get the layer number, sign based on the side
    ln = int(re.sub(r"[^0-9-]","",layer_name))
    print (" === layer_name = ", layer_name)
    print (" === ln = ", ln)
    if side != 0:
        ln *= side

    layer_cells[ln] = cells
    n_tot_cells += cells

print(">>>",layer_cells)

n_layers = len(layer_cells.keys())


#######################################
# prepare histograms

# get the hits collection name
collection = detector_dict["hitsCollection"]

# list of the histograms that will be saved in output ROOT file
histograms = []

# Set some of the binning
z_range = int(detector_dict["max_z"]*1.1)
r_range = int(detector_dict["max_r"]*1.1)
phi_range = 3.5
z_binning = [int(z_range * 2 / bw_z), -z_range, z_range]  # mm binning  
r_binning = [int(r_range * 2 / bw_r), -r_range, r_range]  # mm binning  
phi_binning = [int(phi_range * 2 / bw_phi), -phi_range, phi_range ] # rad binning


# Profile histograms
h_hit_E = ROOT.TH1D("h_hit_E_"+collection, "h_hit_E_"+collection, 100, 0, 5)
h_particle_E = ROOT.TH1D("h_particle_E_"+collection, "h_particle_E_"+collection, 500, 0, 50)
h_particle_pt = ROOT.TH1D("h_particle_pt_"+collection, "h_particle_pt_"+collection, 500, 0, 50)
h_particle_eta = ROOT.TH1D("h_particle_eta_"+collection, "h_particle_eta_"+collection, 100, -5, 5)
histograms += [h_hit_E, h_particle_E, h_particle_pt, h_particle_pt, h_particle_eta]

# ID histogram - use alphanumeric labels, form https://root.cern/doc/master/hist004__TH1__labels_8C.html
h_particle_ID = ROOT.TH1D("h_particle_ID_"+collection, "h_particle_ID_"+collection, 1, 0, 1)
h_particle_ID.SetCanExtend(ROOT.TH1.kAllAxes)   # Allow both axes to extend past the initial range
histograms += [h_particle_ID]

h_hits_x_layer = ROOT.TH1D("h_hits_x_layer_"+collection, "h_hits_x_layer_"+collection, n_layers, -0.5, n_layers-0.5)
h_avg_occ_x_layer = ROOT.TH1D("h_avg_occ_x_layer_"+collection, "h_avg_occ_x_layer_"+collection, n_layers, -0.5, n_layers-0.5)
histograms += [h_hits_x_layer, h_avg_occ_x_layer]

h_occ_x_layer = {}
log_bins = np.logspace(-2,2,50)
for l in layer_cells.keys():
    h_occ_x_layer[l] = ROOT.TH1D(f"h_occ_x_layer{l}_{collection}", f"h_occ_x_layer{l}_{collection}", len(log_bins)-1, log_bins)
    histograms += [h_occ_x_layer[l]]
h_occ = ROOT.TH1D(f"h_occ_tot_{collection}", f"h_occ_tot_{collection}", len(log_bins)-1, log_bins)
histograms += [h_occ]

# Hit maps definitions
hist_zr = ROOT.TH2D("hist_zr_"+collection, "hist_zr_"+collection+"; ; ; hits/(%d#times%d) mm^{2} per event"%(bw_z, bw_r), *z_binning, *r_binning)
hist_xy = ROOT.TH2D("hist_xy_"+collection, "hist_xy_"+collection+"; ; ; hits/(%d#times%d) mm^{2} per event"%(bw_r, bw_r), *r_binning, *r_binning)
hist_zphi = ROOT.TH2D("hist_zphi_"+collection, "hist_zphi_"+collection+"; ; ; hits/(%1.2f#times%d)rad#timesmm per event"%(bw_phi,bw_z), *z_binning, *phi_binning)
histograms += [hist_zr, hist_xy, hist_zphi]

# Timing histograms
h_hit_t = ROOT.TH1D("hist_hit_t_"+collection, "hist_hit_t_"+collection, 200, 0, 5)
h_hit_t_x_layer = ROOT.TH2D("hist_hit_t_map_"+collection, "hist_hit_t_"+collection+"; ; ; hits / events", 200, 0, 10, n_layers, -0.5, n_layers-0.5)
h_hit_t_corr = ROOT.TH1D("hist_hit_t_corr_"+collection, "hist_hit_t_corr_"+collection, 200, -5, 5)
histograms += [h_hit_t, h_hit_t_x_layer, h_hit_t_corr]

# Hit densities per layer
h_z_density_vs_layer_mm = {}
h_phi_density_vs_layer = {}
h_zphi_density_vs_layer = {}
for l in layer_cells.keys():
    h_z_density_vs_layer_mm[l] = ROOT.TH1D(f"hist_z_density_vs_layer{l}_mm_{collection}", f"hist_z_density_vs_layer{l}_mm_{collection}", *z_binning)
    h_phi_density_vs_layer[l] = ROOT.TH1D(f"hist_phi_density_vs_layer{l}_{collection}", f"hist_phi_density_vs_layer{l}_{collection}", *phi_binning)
    h_zphi_density_vs_layer[l] = ROOT.TH2D(f"hist_zphi_vs_layer{l}"+collection, f"hist_zphi_vs_layer{l}"+collection+"; ; ; hits/(%1.2f#times%d)rad#timesmm per event"%(bw_phi,bw_z), *z_binning, *phi_binning)
    histograms += [h_z_density_vs_layer_mm[l], h_phi_density_vs_layer[l], h_zphi_density_vs_layer[l]]

n_events = events_per_file * len(list_input_files)

fill_weight = 1. / (n_events)

#######################################
# run the event loop

print("Initializing the PODIO reader...")

podio_reader = root_io.Reader(list_input_files)

# Check if collection is valid and setup a cell ID decoder
metadata = podio_reader.get("metadata")[0]

id_encoding = metadata.get_parameter(collection+"__CellIDEncoding")
decoder = ROOT.dd4hep.BitFieldCoder(id_encoding)
#print("HERE",decoder.fieldDescription()) #get possible values

max_occ_per_layer = defaultdict(float)
max_occ_per_layer_evt = defaultdict(str)
processed_events = 0

for i,event in enumerate(podio_reader.get(tree_name)):

    if i >= n_events and events_per_file > 0:
        break

    if i % 100 == 0: # print a message every 100 processed files
        print('processing event number = ' + str(i) + ' / ' + str(n_events) )

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
        cell_fired = False

        # Skip cells firing multiple times
        if cell_id in fired_cells:
            # print(cell_id,"has already fired in this event!")
            cell_fired = True
        else:
            fired_cells.append(cell_id)

        #TODO: Handle better cell_id info
        layer = get_layer(cell_id, decoder, sub_detector, detector_dict["typeFlag"])

        x_mm = hit.getPosition().x
        y_mm = hit.getPosition().y
        z_mm = hit.getPosition().z
        r_mm = math.sqrt(math.pow(x_mm, 2) + math.pow(y_mm, 2))
        phi = math.acos(x_mm/r_mm) * math.copysign(1, y_mm)

        if is_calo_hit:
            t = -999  # Timing not available for MutableSimCalorimeterHit
        else:
            t = hit.getTime()

        hit_distance = (x_mm**2 + y_mm**2 + z_mm**2)**0.5

        if x_mm < 0.:
            r_mm = -1. * r_mm

        # Deposited energy, converted to MeV
        E_hit = 0
        if is_calo_hit:
            E_hit = hit.getEnergy() * 1e3  # For calorimeter hits
        else:
            E_hit = hit.getEDep() * 1e3  # For tracker hits

        # For layer percentage occupancy, count cells only once
        if not cell_fired:
            fired_cells_x_layer[layer] += 1

        # Fill the hits histograms
        h_hit_E.Fill(E_hit, fill_weight)
        h_hits_x_layer.Fill(layer, fill_weight)

        hist_zr.Fill(z_mm, r_mm, fill_weight)
        hist_xy.Fill(x_mm, y_mm, fill_weight)
        hist_zphi.Fill(z_mm, phi, fill_weight)

        h_hit_t.Fill(t, fill_weight)
        h_hit_t_x_layer.Fill(t, layer, fill_weight)

        h_hit_t_corr.Fill(t - (hit_distance / C_MM_NS), fill_weight)

        h_z_density_vs_layer_mm[layer].Fill(z_mm, fill_weight)
        h_phi_density_vs_layer[layer].Fill(phi, fill_weight)
        h_zphi_density_vs_layer[layer].Fill(z_mm, phi, fill_weight)

        if not is_calo_hit:
            particle = hit.getParticle()

            particle_p4 = ROOT.Math.LorentzVector('ROOT::Math::PxPyPzM4D<double>')(particle.getMomentum().x, particle.getMomentum().y, particle.getMomentum().z, particle.getMass())

            # Fill MC particle histograms
            h_particle_E.Fill(particle_p4.E(), fill_weight)
            h_particle_pt.Fill(particle_p4.pt(), fill_weight)
            h_particle_eta.Fill(particle_p4.eta(), fill_weight)
            h_particle_ID.Fill(str(particle.getPDG()), fill_weight)

    # <--- end of the hits loop
    
    tot_occupancy = 0
    for l, fc in fired_cells_x_layer.items():
        tot_occupancy += fc
        #print(f" fc:{fc}, l:{l},  nc:{layer_cells[l]} ")
        layer_occupancy = fc / layer_cells[l] * 100
        #print(f" fc:{fc}, nc:{layer_cells[l]}  , occ:{layer_occupancy}")

        h_avg_occ_x_layer.Fill(l, layer_occupancy * fill_weight)
        h_occ_x_layer[l].Fill(layer_occupancy, fill_weight)

        if layer_occupancy > max_occ_per_layer[l]:
            max_occ_per_layer[l] = layer_occupancy

            max_occ_file = int(i / events_per_file)
            max_occ_event = i -  max_occ_file * events_per_file
            max_occ_per_layer_evt[l] = f"event: {max_occ_event}, file: {list_input_files[max_occ_file]}"

    tot_occupancy = tot_occupancy / n_tot_cells * 100
    h_occ.Fill(tot_occupancy, fill_weight)

    processed_events +=1
# <--- end of the event loop




print("####################")
print("# Occupancy max values:")
for i in max_occ_per_layer.keys():
    print(f"max per layer {i}: {max_occ_per_layer[i]} ({max_occ_per_layer_evt[i]})")
print("####################")

if draw_maps:
    draw_map(hist_zr, "z [mm]", "r [mm]", sample_name+"_map_zr_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_map(hist_xy, "x [mm]", "y [mm]", sample_name+"_map_xy_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_map(hist_zphi, "z [mm]", "#phi [rad]", sample_name+"_map_zphi_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_map(h_hit_t_x_layer, "timing [ns]", "layer number", sample_name+"_map_timing_"+str(n_events)+"evt_"+sub_detector, collection)
    
    for l, h in h_zphi_density_vs_layer.items():
        draw_map(h, "z [mm]", "#phi [rad]", sample_name+f"_map_zphi_layer{l}_"+str(n_events)+"evt_"+sub_detector, collection)

if draw_hists: 
    draw_hist(h_hit_E, "Deposited energy [MeV]", "Hits / events",  sample_name+"_hit_E_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_hist(h_hit_t, "Timing [ns]", "Hits / events",  sample_name+"_hit_t_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_hist(h_hit_t_corr, "Timing - TOF [ns]", "Hits / events",  sample_name+"_hit_t_corr_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_hist(h_hits_x_layer, "Layer number", "Hits / events",  sample_name+"_hits_x_layer_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_hist(h_avg_occ_x_layer, "Layer number", "Occupancy [%] / events",  sample_name+"_avg_occ_x_layer_"+str(n_events)+"evt_"+sub_detector, collection, log_y=True)
    for l, h  in h_occ_x_layer.items():
        draw_hist(h, "Occupancy [%]", "Entries / events",  sample_name+f"_occ_x_layer{l}_"+str(n_events)+"evt_"+sub_detector, collection, log_x=True)
    draw_hist(h_occ, "Occupancy [%]", "Entries / events",  sample_name+f"_occ_tot_"+str(n_events)+"evt_"+sub_detector, collection, log_x=True)

    draw_hist(h_particle_E, "MC particle energy [GeV]", "Particle / events",  sample_name+"_particle_E_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_hist(h_particle_pt, "MC particle p_{T} [GeV]", "Particle / events",  sample_name+"_particle_pt_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_hist(h_particle_eta, "MC particle #eta", "Particle / events",  sample_name+"_particle_eta_"+str(n_events)+"evt_"+sub_detector, collection)

    h_particle_ID.GetXaxis().LabelsOption("v>")  # vertical labels, sorted by decreasing values
    draw_hist(h_particle_ID, "MC particle PDG ID", "Hits / events",  sample_name+"_particle_ID_"+str(n_events)+"evt_"+sub_detector, collection)

    for l in h_z_density_vs_layer_mm.keys():
        draw_hist(h_z_density_vs_layer_mm[l], "z [mm]","Hits/event",  sample_name+f"_zDensity_layer{l}_"+str(n_events)+"evt_"+sub_detector, collection)
        draw_hist(h_phi_density_vs_layer[l],  "phi","Hits/event",  sample_name+f"_phiDensity_layer{l}_"+str(n_events)+"evt_"+sub_detector, collection)

# Write the histograms to the output file
output_file_name = f"{sample_name}_{output_file_name}_{sub_detector}.root"
with ROOT.TFile(output_file_name,"RECREATE") as f:
    for h in histograms:
        h.Write()

print("Histograms saved in:", output_file_name)


if processed_events != n_events:
    print(f"ERROR: processed events ({processed_events}) != total expected events ({n_events}) !!!")
    print("       ===> Histograms normalization is WRONG!")
