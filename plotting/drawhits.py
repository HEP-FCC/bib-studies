#!/usr/bin/env python

from collections import defaultdict
import math
import numpy as np
import argparse

# import dd4hep as dd4hepModule
# from ROOT import dd4hep
from podio import root_io
import ROOT

from particle import Particle

from constants import C_MM_NS
from helpers import path_to_list, sorted_n_files, load_json, layer_number_from_string, simplify_dict, is_calo, is_endcap, DetFilePath
from visualization import setup_root_style, draw_hist, draw_map


######################################
# option parser

def parse_args():
    parser = argparse.ArgumentParser(description='drawhits.py', 
        epilog='Example:\ndrawhits.py -s MuonTaggerBarrel -d $BIB_STUDIES/detectors_dicts/ALLEGRO_o1_v03_DetectorDimensions.json -e -1 -D 0',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-D', '--debugLevel',
                    type=int, default=1,
                    help='debug level (0:quiet, 1:default, 2:verbose)')
    parser.add_argument('-i', '--infilePath',
                    type=str, default='/eos/home-s/sfranche/FCC/samples/bib/ipc/aciarma_4IP_2024may29_Z/CADbp_ALLEGRO_o1_v03_r2025-05-29_3998.root',
                    help='path to input file or input directory')
    #parser.add_argument('-o', '--outputFileTag',
    #                type=str, default='hits',
    #                help='prefix of output root file, the sample name will be added as prefix and the detector as suffix')
    parser.add_argument('-n', '--numberOfFiles',
                    type=int, default=1,
                    help='number of files to consider (1 event per file), put -1 to take all files')
    parser.add_argument('-e', '--numberOfEvents',
                    type=int, default=1,
                    help='number of events per file to consider, put -1 to take all events')
    parser.add_argument('-t', '--tree',
                    type=str, default='events',
                    help='name of the tree in the root file')
    parser.add_argument('--sample',
                    type=str, default='ipc',
                    help='Brief sample name for plot titles and output files (e.g. ipc, injbkg, )')
    parser.add_argument('-d', '--detDictFile',
                    type=str, default='$BIB_STUDIES/detectors_dicts/ALLEGRO_o1_v03_DetectorDimensions.json',
                    help='JSON dictionary with some key detector parameters')
    parser.add_argument('-a', '--assumptions',
                    type=str, default=None,
                    help='JSON dictionary with assumptions')
    parser.add_argument('-s', '--subDetector',
                    type=str, default='VertexBarrel',
                    help='Name of the sub detector system')
    parser.add_argument('-p', '--draw_hists',
                    action="store_true",
                    help='activate drawing of 1D histograms')
    parser.add_argument('-m', '--draw_maps',
                    action="store_true",
                    help='activate drawing of number of maps plots')
    parser.add_argument('-r','--bin_width_r',
                    type=int, default=None,
                    help='bin width for radius r in mm')
    parser.add_argument('-z','--bin_width_z',
                    type=int, default=None,
                    help='bin width for z in mm')
    parser.add_argument('--bin_width_phi',
                    type=float, default=None,
                    help='bin width for azimuthal angle phi')
    parser.add_argument('--stat_box',
                    action="store_true",
                    help='draw the ROOT stat box')
    parser.add_argument('--skip_layers',
                    action="store_true",
                    help='Skip plots per single layer, useful for sub-detectors with many layers (e.g. drift chamber)')
    parser.add_argument('--integration_time',
                    type=int, default=1,
                    help='Integration time in  number of events, to estimate hits pile-up.')
    parser.add_argument('--energy_per_layer',
                    action="store_true",
                    help='Save histograms of hit energy per layer.')
    parser.add_argument('--digi',
                    action='store_true',
                    help='run on digitized hits (requires to pass also the assumptions dict with info on the digitized collection).')
    parser.add_argument('--e_cut',
                    action='store_true',
                    help='apply energy cut (requires to pass also the assumptions dict with energy threshold).')

    return (parser.parse_args())

options = parse_args()

debug = options.debugLevel
n_files = options.numberOfFiles
events_per_file = options.numberOfEvents
input_path = options.infilePath
tree_name = options.tree
sample_name = options.sample
det_file =  DetFilePath(options.detDictFile)
assumptions_path = options.assumptions
sub_detector = options.subDetector
draw_maps = options.draw_maps
draw_hists = options.draw_hists
bw_r = options.bin_width_r
bw_z = options.bin_width_z
bw_phi = options.bin_width_phi
stat_box = options.stat_box
skip_plot_per_layer = options.skip_layers
integration_time = options.integration_time
energy_per_layer = options.energy_per_layer
digi = options.digi
e_cut = options.e_cut

######################################
# style

setup_root_style(stat_box)
#to avoid canvas showing up slowing us down
ROOT.gROOT.SetBatch(True)
#to avoid canvas->Print printouts
ROOT.gErrorIgnoreLevel = ROOT.kWarning


#######################################
# functions

print("Importing detector specific functions for:", det_file.short)
match det_file.short:
    case "ALLEGRO":
        from ALLEGRO import get_layer,get_module,get_sensor
    case "IDEA":
        from IDEA import get_layer,get_module,get_sensor
    case "CLD":
        from CLD import get_layer
    case "ILD_FCCee":
        from ILD_FCCee import get_layer
    case _:
        raise NotImplementedError(f"get_layer not implemented for detector {det_file.short}")

#todo: fix
def draw_eta_line(eta_value):
    """Draw a line of constant pseudorapidity eta on the current canvas."""
    eta = eta_value
    theta = 2.0 * math.atan(math.exp(-eta))
    tan_theta = math.tan(theta)

    zmax = 3000.0  # mm (adapt to your plot range)
    n = 1000

    graph = ROOT.TGraph(n)
    #for i in range(n):
    z = 0#i * zmax / (n - 1)
    r = abs(z) * tan_theta
    graph.SetPoint(0, z, r)
    z=300
    r = abs(z) * tan_theta
    graph.SetPoint(1, z, r)

    graph.SetLineColor(ROOT.kRed)
    #graph.SetLineStyle(ROOT.kDashed)
    graph.SetLineWidth(2)
    graph.Draw("*Lsame")

#######################################
# parse the input path and/or files


# list n_files in the path directory ending in .root to consider for plotting
# use sorted list to make sure to always take the same files
if(debug>0): print("Parsing input path:", input_path)
#make a suffic for eventual filename from input_path
#if input is a .root file, then use the filename without extension,
#if input is a directory, then use the directory name
suffix_from_input = ""
if(input_path.endswith(".root")):
    suffix_from_input = input_path.split("/")[-1].replace(".root","")
else:
    suffix_from_input = input_path.strip("/").split("/")[-1]
#last protection: replace any space or special char in suffix_from_input with underscore
suffix_from_input = "".join([c if c.isalnum() else "_" for c in suffix_from_input])

list_files = path_to_list(input_path)
list_input_files = sorted_n_files(list_files, n_files)

print("Initializing the PODIO reader...")
podio_reader = root_io.Reader(list_input_files)

# Check if collection is valid and setup a cell ID decoder
metadata = podio_reader.get("metadata")[0]

nEvents_fullFileList = len(podio_reader.get("events"))
print(f"Number of events in the full file list: {nEvents_fullFileList}")

# Read the detector and assumptions dictionaries
detector_dict = load_json(det_file.path, sub_detector)
detector_type = detector_dict["typeFlag"]
assumptions = load_json(assumptions_path, sub_detector) if assumptions_path else None  #TODO: set default path to the /templates folder

# get the hits collection name
collection =  detector_dict["hitsCollection"]
if digi:
    collection = assumptions["digitized_hits"]["collection"]
id_encoding = metadata.get_parameter(collection+"__CellIDEncoding")
decoder = ROOT.dd4hep.BitFieldCoder(id_encoding)
#print("HERE",decoder.fieldDescription()) #get possible values

# Parse the layer name in to an integer number
print("Retrieve the number of cells per layer")
layer_cells = {}
n_tot_cells = 0
for det_element, cells in detector_dict["det_element_cells"].items():
    ln = layer_number_from_string(det_element)
    layer_cells[ln] = cells
    n_tot_cells += cells

n_layers = len(layer_cells.keys())
print(">>>",dict(sorted(layer_cells.items())))

print("Retrieve the number of sensors per module/cell (assuming this is the same everywhere within the layer")
sensors_per_module_map = {}
if detector_dict.get("sensors_per_module_map") != {}:
    sensors_per_module_map = simplify_dict(detector_dict["sensors_per_module_map"])
else:
    sensors_per_module_map = {i: 1 for i in layer_cells.keys()}  # default to 1 sensor per module
    print(">>> 'sensors_per_module_map' not found in detector dictionary, assuming 1 sensor per module for all layers.")

print(">>>", sensors_per_module_map)

if not layer_cells:
    print("WARNING: map of cells per layer is empty.")
    n_tot_cells = 1

if (e_cut or digi) and assumptions == None:
    raise ValueError(f"Given the input settings (e_cut={e_cut}, digi={digi})," +
                     "an 'assumptions' JSON is needed as input (-a/--assumptions <assumptions.json>)")

# Set the energy thresholds, if required
E_thr_MeV = 0
if e_cut:
    if digi:
        E_thr_MeV = assumptions["digitized_hits"]["E_thr_MeV"]
    else:
        E_thr_MeV = assumptions["E_thr_MeV"]

    # if the thresholds are defined for each layer, simplify the dictionary keys
    if isinstance(E_thr_MeV,dict):
        E_thr_MeV = simplify_dict(E_thr_MeV) 

    print(">>> Cutting Hits with energy below:", E_thr_MeV)

#set of cellIDs that fired in the run (set=> unique entries)
fullRun_fired_cellIDs = set()

#######################################
# prepare histograms

# list of the histograms that will be saved in output ROOT file
histograms = []

# Set some of the binning
x_range = int(detector_dict["max_r"]*1.2)
y_range = int(detector_dict["max_r"]*1.2)
z_range = int(detector_dict["max_z"]*1.2)
r_range = int(detector_dict["max_r"]*1.2)
eta_range = 5
eta_bins = 100
eta_bin_size = (2*eta_range)/eta_bins
if debug>1: print(f"z_range: {z_range}, r_range: {r_range}")
phi_range = 3.5

if not bw_z:
    bw_z = (z_range * 2) /  100
bw_x = (x_range * 2) /  1000
bw_y = (y_range * 2) /  1000
if not bw_r:
    bw_r = (r_range * 2) /  100
if not bw_phi:
    bw_phi = (phi_range * 2) /  100

x_binning = [int(x_range * 2 / bw_x), -x_range, x_range]  # mm binning
y_binning = [int(y_range * 2 / bw_y), -y_range, y_range]  # mm binning
z_binning = [int(z_range * 2 / bw_z), -z_range, z_range]  # mm binning
r_binning = [int(r_range * 2 / bw_r), -r_range, r_range]  # mm binning
phi_binning = [int(phi_range * 2 / bw_phi), -phi_range, phi_range ] # rad binning
if debug>1: print(f"z_binning: {z_binning}, r_binning: {r_binning}, phi_binning: {phi_binning}")

layer_binning = [n_layers + 1, -0.5, n_layers + 0.5]
if is_endcap(detector_type):    # Binning from -negative to positive layer numbers, with a dummy 0 layer in the middle
    max_l = int(n_layers / 2) + 0.5
    layer_binning = [n_layers + 1, -max_l, +max_l]


# Profile histograms
histograms += [
    h_hit_x_mm               := ROOT.TH1D("h_hit_x_mm_"+collection               , "h_hit_x_mm_"+collection               +"; [mm] ; hits;", 10000, -x_range, x_range),
    h_hit_y_mm               := ROOT.TH1D("h_hit_y_mm_"+collection               , "h_hit_y_mm_"+collection               +"; [mm] ; hits;", 10000, -y_range, y_range),
    h_hit_z_mm               := ROOT.TH1D("h_hit_z_mm_"+collection               , "h_hit_z_mm_"+collection               +"; [mm] ; hits;", 10000, -z_range, z_range),
    h_hit_r_mm               := ROOT.TH1D("h_hit_r_mm_"+collection               , "h_hit_r_mm_"+collection               +"; [mm] ; hits;",10000, -r_range, r_range),
    h_hit_eta                := ROOT.TH1D("h_hit_eta_"+collection                , "h_hit_eta_"+collection                +"; [eta]; hits;", eta_bins*2, -eta_range, eta_range),
    h_hit_rateDensity_VS_eta := ROOT.TH1D("h_hit_rateDensity_VS_eta_"+collection , "h_hit_rateDensity_VS_eta_"+collection +"; [eta]; [MHz/cm^2];", eta_bins, 0, eta_range),
    h_hit_E_MeV    := ROOT.TH1D("h_hit_E_MeV_"+collection   , "h_hit_E_MeV_"+collection   +"; [MeV]; hits;", 500, 0, 50),
    h_hit_E_keV    := ROOT.TH1D("h_hit_E_keV_"+collection   , "h_hit_E_keV_"+collection   +"; [keV]; hits;", 500, 0, 500),
    h_hit_particle_E   := ROOT.TH1D("h_hit_particle_E_"+collection  , "h_hit_particle_E_"+collection  +"; [E] ; hits;", 500, 0, 50),
    h_hit_particle_pt  := ROOT.TH1D("h_hit_particle_pt_"+collection , "h_hit_particle_pt_"+collection +"; [E/c] ; hits;", 500, 0, 50),
    h_hit_particle_eta := ROOT.TH1D("h_hit_particle_eta_"+collection, "h_hit_particle_eta_"+collection+"; [eta] ; hits;", 100, -5, 5),
]

# ID histogram - use alphanumeric labels, form https://root.cern/doc/master/hist004__TH1__labels_8C.html
histograms += [ h_hit_particle_ID := ROOT.TH1D("h_hit_particle_ID_"+collection, "h_hit_particle_ID_"+collection, 1, 0, 1) ]
h_hit_particle_ID.SetCanExtend(ROOT.TH1.kAllAxes)   # Allow both axes to extend past the initial range
histograms += [ h_hit_particle_ID_E_MeV := ROOT.TH2D("h_hit_particle_ID_E_MeV_"+collection, "h_hit_particle_ID_E_MeV_"+collection, 1, 0, 1, 500, 0, 50) ]
h_hit_particle_ID_E_MeV.SetCanExtend(ROOT.TH1.kXaxis)   # Allow both axes to extend past the initial range

h_hit_E_MeV_x_layer = {}
if energy_per_layer:
    for l in layer_cells.keys():
        h_hit_E_MeV_x_layer[l] = ROOT.TH1D(f"h_hit_E_MeV_layer{l}_{collection}", f"h_hit_E_MeV_layer{l}_{collection}; hit energy [MeV]", 500, 0, 50)
        histograms += [h_hit_E_MeV_x_layer[l]]

# Hit maps definitions
histograms += [
    hit_zr := ROOT.TH2D("hit_zr_"+collection, "hit_zr_"+collection+";  z (bin=%dmm) ;r (bin=%dmm) ; hits/(%d#times%d) mm^{2} per event"%(bw_z, bw_r,bw_z, bw_r), *z_binning, *r_binning),
    hit_xy := ROOT.TH2D("hit_xy_"+collection, "hit_xy_"+collection+";  x (bin=%dmm) ;y (bin=%dmm) ; hits/(%d#times%d) mm^{2} per event"%(bw_r, bw_r,bw_r, bw_r), *r_binning, *r_binning),
    hit_zphi := ROOT.TH2D("hit_zphi_"+collection, "hit_zphi_"+collection+"; phi (bin=%1.2frad) ; z (bin=%dmm); hits/(%1.2f#times%d)rad#timesmm per event"%(bw_phi,bw_z,bw_phi,bw_z), *z_binning, *phi_binning),
]

# Timing histograms
histograms += [
    h_hit_t := ROOT.TH1D("hist_hit_t_"+collection, "hist_hit_t_"+collection, 200, 0, 20),
    h_hit_t_x_layer := ROOT.TH2D("hist_hit_t_map_"+collection, "hist_hit_t_"+collection+"; ; ; hits / events", 200, 0, 20, *layer_binning),
    h_hit_t_corr := ROOT.TH1D("hist_hit_t_corr_"+collection, "hist_hit_t_corr_"+collection, 200, -10, 10),
]

histograms += [
    h_avg_hits_x_layer := ROOT.TH1D("h_avg_hits_x_layer_"+collection, "h_avg_hits_x_layer_"+collection+";layer index;avg hits/evt", *layer_binning),
    h_avg_firing_cells_x_layer := ROOT.TH1D("h_avg_firing_cells_x_layer_"+collection, "h_avg_firing_cells_x_layer_"+collection+";layer index;avg cellsFired/evt", *layer_binning),
    h_avg_occ_x_layer := ROOT.TH1D("h_avg_occ_x_layer_"+collection, "h_avg_occ_x_layer_"+collection+";layer index;avg occupancy/evt [%]", *layer_binning),
]

# Occupancy histograms, defined as fraction of cells that fired
# Requires that the total cells are correctly computed and stored in detector_dict
h_occ_x_layer = {}
log_bins = np.hstack([0,np.logspace(-2,2,50)]) #add 0 value to log bins (to avoid underflow in some histos)
for l in layer_cells.keys():
    h_occ_x_layer[l] = ROOT.TH1D(f"h_occ_x_layer{l}_{collection}", f"h_occ_x_layer{l}_{collection} ; fired cells / total layer cells [%];;", len(log_bins)-1, log_bins)
    histograms += [h_occ_x_layer[l]]
histograms += [ h_occ := ROOT.TH1D(f"h_occ_tot_{collection}", f"h_occ_tot_{collection} ; total fired cells / total layer cells [%]", len(log_bins)-1, log_bins) ]


# Hit densities per layer
if not skip_plot_per_layer:
    h_z_density_vs_layer_mm = {}
    h_phi_density_vs_layer = {}
    h_zphi_density_vs_layer = {}
    h_xy_density_vs_layer = {}
    for l in layer_cells.keys():
        h_z_density_vs_layer_mm[l] = ROOT.TH1D(f"hist_z_density_vs_layer{l}_mm_{collection}", f"hist_z_density_vs_layer{l}_mm_{collection}", *z_binning)
        h_phi_density_vs_layer[l] = ROOT.TH1D(f"hist_phi_density_vs_layer{l}_{collection}", f"hist_phi_density_vs_layer{l}_{collection}", *phi_binning)
        h_zphi_density_vs_layer[l] = ROOT.TH2D(f"hist_zphi_vs_layer{l}_{collection}", f"hist_zphi_vs_layer{l}_"+collection+"; ; ; hits/(%1.2f#times%d)rad#timesmm per event"%(bw_phi,bw_z), *z_binning, *phi_binning)
        h_xy_density_vs_layer[l] = ROOT.TH2D(f"hist_xy_vs_layer{l}_{collection}", f"hist_xy_vs_layer{l}_{collection};  x (bin=%dmm) ;y (bin=%dmm) ; hits/(%d#times%d) mm^{2} per event"%(bw_r, bw_r,bw_r, bw_r), *r_binning, *r_binning)

        histograms += [h_z_density_vs_layer_mm[l], h_phi_density_vs_layer[l], h_zphi_density_vs_layer[l], h_xy_density_vs_layer[l]]

    # Per-module histograms
    h_avg_hits_x_layer_x_module = {}
    bin_start = 0
    cell_binning = {}
    for l in layer_cells.keys():
        cell_binning[l] = [int(layer_cells[l]/sensors_per_module_map[l]), -0.5 + bin_start, layer_cells[l]/sensors_per_module_map[l]-0.5  + bin_start]

        h_avg_hits_x_layer_x_module[l] = ROOT.TH1D(f"h_avg_hits_x_layer{l}_x_module_{collection}", f"h_avg_hits_x_layer{l}_x_module_{collection};Module ID;Hits / event", *cell_binning[l])
        histograms += [h_avg_hits_x_layer_x_module[l]]
        if(sub_detector=="VertexBarrel" or sub_detector=="SiWrB"):
            bin_start += layer_cells[l]/sensors_per_module_map[l]

    h_pu_x_layer = {}
    for l in layer_cells.keys():
        h_pu_x_layer[l] = ROOT.TH1D(f"h_pu_x_layer{l}_{collection}", f"h_pu_x_layer{l}_{collection}", integration_time+1, -0.5, integration_time+0.5)
        histograms += [h_pu_x_layer[l]]
        
# Pile-up histograms
histograms += [
    h_tot_pu := ROOT.TH1D(f"h_tot_pu_{collection}", f"h_tot_pu_{collection}", integration_time+1, -0.5, integration_time+0.5),
    h_avg_pu_x_layer := ROOT.TH1D(f"h_avg_pu_x_layer_{collection}", f"h_avg_pu_x_layer_{collection}", *layer_binning)
]


if events_per_file==-1: 
    events_per_file = nEvents_fullFileList
n_events = events_per_file * len(list_input_files)


fill_weight = 1. / (n_events)
if debug>1: print(f"fill weight: {fill_weight}")

#######################################
# run the event loop

# Counters for the event / hit loops
processed_events = 0
processed_hits = 0

integrating_channels = defaultdict(int)
pile_up_counter = defaultdict(int)

for i,event in enumerate(podio_reader.get(tree_name)):

    if i >= n_events and events_per_file > 0:
        break

    if i % 100 == 0: # print a message every 100 processed files
        print('processing event number = ' + str(i) + ' / ' + str(n_events) )

    fired_cells = []
    fired_cells_x_layer = {l: 0 for l in layer_cells.keys()}  # Initialize to 0 for all the layers to ensure that also events with 0 occupancy are correctly counted

    # Loop over the hits
    for hit in event.get(collection):
        if(debug>1): print("Processing hit:", hit)
        # Hits doxy:
        # https://edm4hep.web.cern.ch/classedm4hep_1_1_mutable_sim_tracker_hit.html
        # https://edm4hep.web.cern.ch/classedm4hep_1_1_mutable_sim_calorimeter_hit.html
        # https://edm4hep.web.cern.ch/classedm4hep_1_1_m_c_particle.html

        is_calo_hit = is_calo(detector_type)

        cell_id = hit.cellID()
        cell_fired = False

        # skip cells firing multiple times
        if cell_id in fired_cells:
            # print(cell_id,"has already fired in this event!")
            cell_fired = True
        else:
            fired_cells.append(cell_id)
            # count cellIDs that fired in the run
            fullRun_fired_cellIDs.add(cell_id)

        #TODO: Handle better cell_id info
        layer_n = get_layer(cell_id, decoder, sub_detector, detector_type)
        module_n = get_module(cell_id, decoder, sub_detector, detector_type)
        # sensor_n = get_sensor(cell_id, decoder, sub_detector, detector_type)  # not used for now

        E_hit_thr = 0
        if isinstance(E_thr_MeV, dict):
            E_hit_thr = E_thr_MeV[layer_n]

        if is_calo_hit:
            t = -999  # Timing not available for MutableSimCalorimeterHit
        else:
            E_hit_thr = E_thr_MeV

        # Deposited energy, converted to MeV
        E_hit = 0
        if is_calo_hit:
            E_hit = hit.getEnergy() * 1e3  # For calorimeter hits
        else:
            E_hit = hit.getEDep() * 1e3  # For tracker hits

        # Apply cut deposited energy
        if E_hit >= E_hit_thr or (not e_cut):

            x_mm = hit.getPosition().x
            y_mm = hit.getPosition().y
            z_mm = hit.getPosition().z
            r_mm = math.sqrt(math.pow(x_mm, 2) + math.pow(y_mm, 2))
            phi = math.acos(x_mm/r_mm) * math.copysign(1, y_mm)
            theta = math.atan2(r_mm, z_mm)
            eta  = -math.log(math.tan(theta / 2.0))
            #eta bin area in cm2
            area_cm2 = 2 * math.pi * (r_mm/10)**2 * math.cosh(eta) * eta_bin_size


            if is_calo_hit:
                t = -999  # Timing not available for MutableSimCalorimeterHit
            else:
                t = hit.getTime()

            hit_distance = (x_mm**2 + y_mm**2 + z_mm**2)**0.5

            if x_mm < 0.:
                r_mm = -1. * r_mm

            # For layer percentage occupancy, count cells only once
            if not cell_fired:
                fired_cells_x_layer[layer_n] += 1

            if(debug>1): print(" cell_id:", cell_id, " layer_n:", layer_n, " x_mm:", x_mm, " y_mm:", y_mm, " z_mm:", z_mm, " r_mm:", r_mm, " phi:", phi)

            # Fill the hits histograms
            h_hit_x_mm.Fill(x_mm, fill_weight)
            h_hit_y_mm.Fill(y_mm, fill_weight)
            h_hit_z_mm.Fill(z_mm, fill_weight)
            h_hit_r_mm.Fill(r_mm, fill_weight)
            h_hit_eta.Fill(eta, fill_weight)
            #foreach event, fill the eta bin, scaled by the bin area in cm2 => <hits>/evt/cm2, but 40MHz evt rate => multiply to getMHz/cm2
            h_hit_rateDensity_VS_eta.Fill(abs(eta), 40.0*1./area_cm2)  # hits/cm2 => X40MHz for MHz/cm2
            h_hit_E_MeV.Fill(E_hit, fill_weight)
            h_hit_E_keV.Fill(E_hit * 1e3, fill_weight)
            h_avg_hits_x_layer.Fill(layer_n, fill_weight)
            if not cell_fired:
                h_avg_firing_cells_x_layer.Fill(layer_n, fill_weight)

            if not skip_plot_per_layer:
                h_avg_hits_x_layer_x_module[layer_n].Fill(module_n, fill_weight)

            hit_zr.Fill(z_mm, r_mm, fill_weight)
            hit_xy.Fill(x_mm, y_mm, fill_weight)
            hit_zphi.Fill(z_mm, phi, fill_weight)

            h_hit_t.Fill(t, fill_weight)
            h_hit_t_x_layer.Fill(t, layer_n, fill_weight)

            h_hit_t_corr.Fill(t - (hit_distance / C_MM_NS), fill_weight)

            if layer_cells and not skip_plot_per_layer:
                h_z_density_vs_layer_mm[layer_n].Fill(z_mm, fill_weight)
                h_phi_density_vs_layer[layer_n].Fill(phi, fill_weight)
                h_zphi_density_vs_layer[layer_n].Fill(z_mm, phi, fill_weight)
                h_xy_density_vs_layer[layer_n].Fill(x_mm, y_mm, fill_weight)

                if energy_per_layer:
                    h_hit_E_MeV_x_layer[layer_n].Fill(E_hit, fill_weight)

            if not is_calo_hit:
                particle = hit.getParticle()

                particle_p4 = ROOT.Math.LorentzVector('ROOT::Math::PxPyPzM4D<double>')(particle.getMomentum().x, particle.getMomentum().y, particle.getMomentum().z, particle.getMass())

                # Fill MC particle histograms
                h_hit_particle_E.Fill(particle_p4.E(), fill_weight)
                h_hit_particle_pt.Fill(particle_p4.pt(), fill_weight)
                h_hit_particle_eta.Fill(particle_p4.eta(), fill_weight)
                pdg=particle.getPDG()
                h_hit_particle_ID.Fill(Particle.from_pdgid(pdg).name+"["+str(pdg)+"]", fill_weight)
                h_hit_particle_ID_E_MeV.Fill(Particle.from_pdgid(pdg).name, particle_p4.E(), fill_weight)

            # Monitor in channels that are integrating signal
            if integration_time > 1:
                if ((layer_n, cell_id) in integrating_channels) and (integrating_channels[layer_n, cell_id] > 1):

                    # count one hit per event during integration time as pileup
                    if not cell_fired:
                        pile_up_counter[layer_n, cell_id] += 1

                else:
                    # Start integration, reset pile-up counter
                    integrating_channels[layer_n, cell_id] = integration_time
                    pile_up_counter[layer_n, cell_id] = 0

        processed_hits += 1
    # <--- end of the hits loop

    # Compute the per-event occupancy
    tot_occupancy = 0

    for l, fc in fired_cells_x_layer.items():
        tot_occupancy += fc
        if layer_cells:
            layer_occupancy = fc / layer_cells[l] * 100
            h_avg_occ_x_layer.Fill(l, layer_occupancy * fill_weight)
            h_occ_x_layer[l].Fill(layer_occupancy, fill_weight)

    tot_occupancy = tot_occupancy / n_tot_cells * 100
    h_occ.Fill(tot_occupancy, fill_weight)

    if integration_time > 1:
        # Update the integration and pileup counters
        reset_channels = []
        for c in integrating_channels.keys():
            integrating_channels[c] -= 1
            if integrating_channels[c] < 1:
                reset_channels.append(c)

        for c in reset_channels:
            integrating_channels.pop(c)
            pu = pile_up_counter.pop(c)
            h_tot_pu.Fill(pu)
            h_pu_x_layer[c[0]].Fill(pu)

    processed_events +=1
# <--- end of the event loop

print("Loop done!")
print(" - Processed events:", processed_events)
print(" - Processed hits:", processed_hits, 
      f"(avg. per event: {processed_hits/processed_events})")
print(" - Number of unique cells fired:", len(fullRun_fired_cellIDs), 
      f"(avg. per event: {len(fullRun_fired_cellIDs)/processed_events})")

# Compute the average pile up hits
if integration_time > 1:
    for l, h in h_pu_x_layer.items():
        h_avg_pu_x_layer.Fill(l, h.GetMean())
        #h_avg_pu_x_layer.SetBinError(l, h.GetMeanError()) #TODO: set the right error

# # Draw TEST histos
# cm1 = ROOT.TCanvas("cm1_", "cm1_", 800, 600)
# h_hit_rateDensity_VS_eta.Draw("COLZ")
# #draw_eta_line(1)
# cm1.SetLogy()
# cm1.Print("test.pdf")
cm2 = ROOT.TCanvas("cm2_", "cm2_", 800, 600)
hit_zr.Draw("COLZ")
#draw_eta_line(1)
graph = ROOT.TGraph(4)
graph.SetPoint(0, -z_range, -z_range* 0.851)
graph.SetPoint(1,  z_range,  z_range* 0.851)
graph.SetPoint(2,  z_range, -z_range* 0.851)
graph.SetPoint(3, -z_range,  z_range* 0.851)

graph.SetLineColor(ROOT.kBlue)
graph.SetLineStyle(ROOT.kDashed)
graph.SetLineWidth(2)
graph.Draw("Lsame")
cm2.Print(f"testEta_{sample_name}_{n_events}evt_{sub_detector}_{suffix_from_input}.pdf")

if draw_maps:
    draw_map(hit_zr, "z [mm]", "r [mm]", sample_name+"_map_zr_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_map(hit_xy, "x [mm]", "y [mm]", sample_name+"_map_xy_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_map(hit_zphi, "z [mm]", "#phi [rad]", sample_name+"_map_zphi_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_map(h_hit_t_x_layer, "timing [ns]", "layer number", sample_name+"_map_timing_"+str(n_events)+"evt_"+sub_detector, collection)

    if not skip_plot_per_layer:
        for l, h in h_zphi_density_vs_layer.items():
            draw_map(h, "z [mm]", "#phi [rad]", sample_name+f"_map_zphi_layer{l}_"+str(n_events)+"evt_"+sub_detector, collection)

if draw_hists: 
    draw_hist(h_hit_E_MeV, "Deposited energy [MeV]", "Hits / events",  sample_name+"_hit_E_MeV_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_hist(h_hit_E_keV, "Deposited energy [keV]", "Hits / events",  sample_name+"_hit_E_keV_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_hist(h_hit_t, "Timing [ns]", "Hits / events",  sample_name+"_hit_t_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_hist(h_hit_t_corr, "Timing - TOF [ns]", "Hits / events",  sample_name+"_hit_t_corr_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_hist(h_avg_hits_x_layer, "Layer number", "Hits / events",  sample_name+"_avg_hits_x_layer_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_hist(h_avg_firing_cells_x_layer, "Layer number", "Firing cells / events",  sample_name+"_avg_firing_cells_x_layer_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_hist(h_occ, "Occupancy [%]", "Entries / events",  sample_name+f"_occ_tot_"+str(n_events)+"evt_"+sub_detector, collection, log_x=True)

    if energy_per_layer:
        for l, h in h_hit_E_MeV_x_layer.items():
            draw_hist(h, "Deposited energy [MeV]", "Hits per event",  sample_name+f"_hit_E_MeV_x_layer{l}_"+str(n_events)+"evt_"+sub_detector, collection)

    if not skip_plot_per_layer:
        for l, h in h_occ_x_layer.items():
            draw_hist(h, "Occupancy [%]", "Entries / events",  sample_name+f"_occ_x_layer{l}_"+str(n_events)+"evt_"+sub_detector, collection, log_x=True)

    draw_hist(h_hit_particle_E, "MC particle energy [GeV]", "Particle / events",  sample_name+"_particle_E_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_hist(h_hit_particle_pt, "MC particle p_{T} [GeV]", "Particle / events",  sample_name+"_particle_pt_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_hist(h_hit_particle_eta, "MC particle #eta", "Particle / events",  sample_name+"_particle_eta_"+str(n_events)+"evt_"+sub_detector, collection)

    h_hit_particle_ID.GetXaxis().LabelsOption("v>")  # vertical labels, sorted by decreasing values
    draw_hist(h_hit_particle_ID, "MC particle PDG ID", "Hits / events",  sample_name+"_particle_ID_"+str(n_events)+"evt_"+sub_detector, collection)

    if not skip_plot_per_layer:
        for l in h_z_density_vs_layer_mm.keys():
            draw_hist(h_z_density_vs_layer_mm[l], "z [mm]","Hits/event",  sample_name+f"_zDensity_layer{l}_"+str(n_events)+"evt_"+sub_detector, collection)
            draw_hist(h_phi_density_vs_layer[l],  "phi","Hits/event",  sample_name+f"_phiDensity_layer{l}_"+str(n_events)+"evt_"+sub_detector, collection)
            draw_hist(h_xy_density_vs_layer[l], "x [mm]", "y [mm]", sample_name+f"_xyDensity_layer{l}_"+str(n_events)+"evt_"+sub_detector, collection, draw_opt="colz")

    if integration_time > 1:
        draw_hist(h_tot_pu, "Number of pileup hits", "Entries",  sample_name+"_tot_pu_"+str(n_events)+"evt_"+sub_detector, collection)
        draw_hist(h_avg_pu_x_layer, "Layer", "Average pileup hits",  sample_name+"_avg_pu_x_layer_"+str(n_events)+"evt_"+sub_detector, collection)

print("Writing histograms...")
# Write the histograms to the output file
output_file_name = f"{sample_name}_{n_events}evt_{sub_detector}_{suffix_from_input}.root"
with ROOT.TFile(output_file_name,"RECREATE") as f:
    for h in histograms:
        if(debug>1): print("Writing histo:", h.GetName())
        h.Write()

print("Histograms saved in:", output_file_name)


if processed_events != n_events:
    print(f"ERROR: processed events ({processed_events}) != total expected events ({n_events}) !!!")
    print("       ===> Histograms normalization is WRONG!")
