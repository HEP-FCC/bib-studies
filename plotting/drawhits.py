#!/usr/bin/env python

from collections import Counter, defaultdict
import math
import numpy as np
import argparse

import dd4hep as dd4hepModule
from ROOT import dd4hep
from podio import root_io
import ROOT

from constants import C_MM_NS
from helpers import path_to_list, sorted_n_files, load_json, layer_number_from_string
from visualization import setup_root_style, draw_hist, draw_map


######################################
# option parser

def parse_args():
    parser = argparse.ArgumentParser(description='drawhits.py', 
        epilog='Example:\ndrawhits.py -s MuonTaggerBarrel -d ALLEGRO_o1_v03_DetectorDimensions.json -e -1 -D 0',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-D', '--debugLevel',
                    type=int, default=1,
                    help='debug level (0:quiet, 1:default, 2:verbose)')
    parser.add_argument('-i', '--infilePath',
                    type=str, default='/eos/home-s/sfranche/FCC/samples/bib/ipc/aciarma_4IP_2024may29_Z/CADbp_ALLEGRO_o1_v03_r2025-05-29_3998.root',
                    help='path to input file or input directory')
    parser.add_argument('-o', '--outputFile',
                    type=str, default='hits',
                    help='Name of output root file, the sample name will be added as prefix and the detector as suffix')
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
                    help='sample name to save plots')
    parser.add_argument('-d', '--detDictFile',
                    type=str, default='ALLEGRO_o1_v03_DetectorDimensions.json',
                    help='JSON dictionary with some key detector parameters')
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
    parser.add_argument('-a', '--bin_width_phi',
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

    return(parser.parse_args())

options = parse_args()

debug = options.debugLevel
n_files = options.numberOfFiles
events_per_file = options.numberOfEvents
input_path = options.infilePath
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
stat_box = options.stat_box
skip_plot_per_layer = options.skip_layers
integration_time = options.integration_time
energy_per_layer = options.energy_per_layer

######################################
# style

setup_root_style(stat_box)
#to avoid canvas showing up slowing us down
ROOT.gROOT.SetBatch(True)
#to avoid canvas->Print printouts
ROOT.gErrorIgnoreLevel = ROOT.kWarning


#######################################
# functions

# Detector types from:
# https://github.com/AIDASoft/DD4hep/blob/master/DDCore/include/DD4hep/DetType.h
is_calo = lambda x: (x & dd4hep.DetType.CALORIMETER) == dd4hep.DetType.CALORIMETER
is_endcap = lambda x: (x & dd4hep.DetType.ENDCAP) == dd4hep.DetType.ENDCAP #DetType_ENDCAP in xml


def get_layer(cell_id, decoder, detector, dtype):
    """
    Run the decoder differently for each detector
    """

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

        # case "ECalBarrel": 
        #     return  0

        case "DCH_v2":
            # Number of layers per super layer could be read from geo file
            nl_x_sl = 8
            layer = decoder.get(cell_id, "layer")
            super_layer = decoder.get(cell_id, "superlayer")
            return (super_layer * nl_x_sl) + layer + 1

        case "VertexDisks":
            # shift layer number by 1 to remove degeneracy of layer 0
            layer = decoder.get(cell_id, "layer") + 1
            side = decoder.get(cell_id, "side")

            return layer * side

        case "MuonTaggerEndcap":
            # Default way: side * layer, where side should be +/- 1
            layer = decoder.get(cell_id, "layer") + 1
            #probably no side available in decoder (see xml readout part)
            theta = decoder.get(cell_id, "theta")
            #print(f"layer={layer}, theta={theta}")

            side=-1
            if theta<168:
              side = 1

            return layer * side

        case _:           
            layer = decoder.get(cell_id, "layer")

            # Get the side, if available
            side = 0
            if is_endcap(dtype):
                side = decoder.get(cell_id, "side")

            if side != 0:
                layer *= side

            return layer


#######################################
# parse the input path and/or files


# list n_files in the path directory ending in .root to consider for plotting
# use sorted list to make sure to always take the same files
list_files = path_to_list(input_path)
list_input_files = sorted_n_files(list_files, n_files)

detector_dict = load_json(detector_dict_path, sub_detector)
detector_type = detector_dict["typeFlag"]

# Parse the layer name in to an integer number
print("Retrieve the number of cells per layer")
layer_cells = {}
n_tot_cells = 0
for de, cells in detector_dict["det_element_cells"].items():
    #for ECalBArrel take the cell from bath - TODO: improve code clarity if possible
    if de == 'bath':
        layer_cells[0] = cells
        n_tot_cells += cells

    ln = layer_number_from_string(de)

    layer_cells[ln] = cells
    n_tot_cells += cells

if not layer_cells:
    print("WARNING: map of cells per layer is empty.")
    n_tot_cells = 1

print(">>>",layer_cells)

n_layers = len(layer_cells.keys())


#######################################
# prepare histograms

# get the hits collection name
collection = detector_dict["hitsCollection"]

# list of the histograms that will be saved in output ROOT file
histograms = []

# Set some of the binning
z_range = int(detector_dict["max_z"]*1.2)
r_range = int(detector_dict["max_r"]*1.2)
if debug>1: print(f"z_range: {z_range}, r_range: {r_range}")
phi_range = 3.5

if not bw_z:
    bw_z = (z_range * 2) /  100
if not bw_r:
    bw_r = (r_range * 2) /  100
if not bw_phi:
    bw_phi = (phi_range * 2) /  100

z_binning = [int(z_range * 2 / bw_z), -z_range, z_range]  # mm binning
r_binning = [int(r_range * 2 / bw_r), -r_range, r_range]  # mm binning
phi_binning = [int(phi_range * 2 / bw_phi), -phi_range, phi_range ] # rad binning
if debug>1: print(f"z_binning: {z_binning}, r_binning: {r_binning}, phi_binning: {phi_binning}")

layer_binning = [n_layers + 1, -0.5, n_layers + 0.5]
if is_endcap(detector_type):
    max_l = int(n_layers / 2) + 0.5
    layer_binning = [n_layers + 1, -max_l, +max_l]


# Profile histograms
histograms += [
    h_hit_E_MeV := ROOT.TH1D("h_hit_E_MeV_"+collection, "h_hit_E_MeV_"+collection, 100, 0, 5),
    h_hit_E_keV := ROOT.TH1D("h_hit_E_keV_"+collection, "h_hit_E_MeV_"+collection, 500, 0, 500),
    h_particle_E := ROOT.TH1D("h_particle_E_"+collection, "h_particle_E_"+collection, 500, 0, 50),
    h_particle_pt := ROOT.TH1D("h_particle_pt_"+collection, "h_particle_pt_"+collection, 500, 0, 50),
    h_particle_eta := ROOT.TH1D("h_particle_eta_"+collection, "h_particle_eta_"+collection, 100, -5, 5),
]

h_hit_E_MeV_x_layer = {}
if energy_per_layer:
    for l in layer_cells.keys():
        h_hit_E_MeV_x_layer[l] = ROOT.TH1D(f"h_hit_E_MeV_layer{l}_{collection}", f"h_hit_E_MeV_layer{l}_{collection}", 100, 0, 5)
        histograms += [h_hit_E_MeV_x_layer[l]]

# ID histogram - use alphanumeric labels, form https://root.cern/doc/master/hist004__TH1__labels_8C.html
histograms += [ h_particle_ID := ROOT.TH1D("h_particle_ID_"+collection, "h_particle_ID_"+collection, 1, 0, 1) ]
h_particle_ID.SetCanExtend(ROOT.TH1.kAllAxes)   # Allow both axes to extend past the initial range

histograms += [
    h_avg_hits_x_layer := ROOT.TH1D("h_avg_hits_x_layer_"+collection, "h_avg_hits_x_layer_"+collection, *layer_binning),
    h_avg_occ_x_layer := ROOT.TH1D("h_avg_occ_x_layer_"+collection, "h_avg_occ_x_layer_"+collection, *layer_binning),
]

# Occupancy histograms, defined as fraction of cells that fired
# Requires that the total cells are correctly computed and stored in detector_dict
h_occ_x_layer = {}
log_bins = np.logspace(-2,2,50)
for l in layer_cells.keys():
    h_occ_x_layer[l] = ROOT.TH1D(f"h_occ_x_layer{l}_{collection}", f"h_occ_x_layer{l}_{collection}", len(log_bins)-1, log_bins)
    histograms += [h_occ_x_layer[l]]
histograms += [ h_occ := ROOT.TH1D(f"h_occ_tot_{collection}", f"h_occ_tot_{collection}", len(log_bins)-1, log_bins) ]

# Hit maps definitions
histograms += [
    hist_zr := ROOT.TH2D("hist_zr_"+collection, "hist_zr_"+collection+"; ; ; hits/(%d#times%d) mm^{2} per event"%(bw_z, bw_r), *z_binning, *r_binning),
    hist_xy := ROOT.TH2D("hist_xy_"+collection, "hist_xy_"+collection+"; ; ; hits/(%d#times%d) mm^{2} per event"%(bw_r, bw_r), *r_binning, *r_binning),
    hist_zphi := ROOT.TH2D("hist_zphi_"+collection, "hist_zphi_"+collection+"; ; ; hits/(%1.2f#times%d)rad#timesmm per event"%(bw_phi,bw_z), *z_binning, *phi_binning),
]

# Timing histograms
histograms += [
    h_hit_t := ROOT.TH1D("hist_hit_t_"+collection, "hist_hit_t_"+collection, 200, 0, 20),
    h_hit_t_x_layer := ROOT.TH2D("hist_hit_t_map_"+collection, "hist_hit_t_"+collection+"; ; ; hits / events", 200, 0, 20, *layer_binning),
    h_hit_t_corr := ROOT.TH1D("hist_hit_t_corr_"+collection, "hist_hit_t_corr_"+collection, 200, -10, 10),
]

# Hit densities per layer
h_z_density_vs_layer_mm = {}
h_phi_density_vs_layer = {}
h_zphi_density_vs_layer = {}
for l in layer_cells.keys():
    h_z_density_vs_layer_mm[l] = ROOT.TH1D(f"hist_z_density_vs_layer{l}_mm_{collection}", f"hist_z_density_vs_layer{l}_mm_{collection}", *z_binning)
    h_phi_density_vs_layer[l] = ROOT.TH1D(f"hist_phi_density_vs_layer{l}_{collection}", f"hist_phi_density_vs_layer{l}_{collection}", *phi_binning)
    h_zphi_density_vs_layer[l] = ROOT.TH2D(f"hist_zphi_vs_layer{l}"+collection, f"hist_zphi_vs_layer{l}"+collection+"; ; ; hits/(%1.2f#times%d)rad#timesmm per event"%(bw_phi,bw_z), *z_binning, *phi_binning)
    histograms += [h_z_density_vs_layer_mm[l], h_phi_density_vs_layer[l], h_zphi_density_vs_layer[l]]

# Pile-up histograms
histograms += [
    h_tot_pu := ROOT.TH1D(f"h_tot_pu_{collection}", f"h_tot_pu_{collection}", integration_time+1, -0.5, integration_time+0.5),
    h_avg_pu_x_layer := ROOT.TH1D(f"h_avg_pu_x_layer_{collection}", f"h_avg_pu_x_layer_{collection}", *layer_binning)
]

h_pu_x_layer = {}
for l in layer_cells.keys():
    h_pu_x_layer[l] = ROOT.TH1D(f"h_pu_x_layer{l}_{collection}", f"h_pu_x_layer{l}_{collection}", integration_time+1, -0.5, integration_time+0.5)
    histograms += [h_pu_x_layer[l]]

n_events = events_per_file * len(list_input_files)

fill_weight = 1. / (n_events)
if debug>1: print(f"fill weight: {fill_weight}")

#######################################
# run the event loop

print("Initializing the PODIO reader...")

podio_reader = root_io.Reader(list_input_files)

# Check if collection is valid and setup a cell ID decoder
metadata = podio_reader.get("metadata")[0]

id_encoding = metadata.get_parameter(collection+"__CellIDEncoding")
decoder = ROOT.dd4hep.BitFieldCoder(id_encoding)
#print("HERE",decoder.fieldDescription()) #get possible values

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
    fired_cells_x_layer = Counter()

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

        #TODO: Handle better cell_id info
        layer_n = get_layer(cell_id, decoder, sub_detector, detector_type)

        x_mm = hit.getPosition().x
        y_mm = hit.getPosition().y
        z_mm = hit.getPosition().z
        r_mm = math.sqrt(math.pow(x_mm, 2) + math.pow(y_mm, 2))
        phi = math.acos(x_mm/r_mm) * math.copysign(1, y_mm)

        if(debug>1): print(" cell_id:", cell_id, " layer_n:", layer_n, " x_mm:", x_mm, " y_mm:", y_mm, " z_mm:", z_mm, " r_mm:", r_mm, " phi:", phi)

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
            fired_cells_x_layer[layer_n] += 1

        # Fill the hits histograms
        h_hit_E_MeV.Fill(E_hit, fill_weight)
        h_hit_E_keV.Fill(E_hit * 1e3, fill_weight)
        h_avg_hits_x_layer.Fill(layer_n, fill_weight)

        hist_zr.Fill(z_mm, r_mm, fill_weight)
        hist_xy.Fill(x_mm, y_mm, fill_weight)
        hist_zphi.Fill(z_mm, phi, fill_weight)

        h_hit_t.Fill(t, fill_weight)
        h_hit_t_x_layer.Fill(t, layer_n, fill_weight)

        h_hit_t_corr.Fill(t - (hit_distance / C_MM_NS), fill_weight)

        if layer_cells:
            h_z_density_vs_layer_mm[layer_n].Fill(z_mm, fill_weight)
            h_phi_density_vs_layer[layer_n].Fill(phi, fill_weight)
            h_zphi_density_vs_layer[layer_n].Fill(z_mm, phi, fill_weight)

            if energy_per_layer:
                h_hit_E_MeV_x_layer[layer_n].Fill(E_hit, fill_weight)

        if not is_calo_hit:
            particle = hit.getParticle()

            particle_p4 = ROOT.Math.LorentzVector('ROOT::Math::PxPyPzM4D<double>')(particle.getMomentum().x, particle.getMomentum().y, particle.getMomentum().z, particle.getMass())

            # Fill MC particle histograms
            h_particle_E.Fill(particle_p4.E(), fill_weight)
            h_particle_pt.Fill(particle_p4.pt(), fill_weight)
            h_particle_eta.Fill(particle_p4.eta(), fill_weight)
            h_particle_ID.Fill(str(particle.getPDG()), fill_weight)

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
        #print(f" fc:{fc}, l:{l},  nc:{layer_cells[l]} ")

        if layer_cells:
            layer_occupancy = fc / layer_cells[l] * 100
            #print(f" fc:{fc}, nc:{layer_cells[l]}  , occ:{layer_occupancy}")

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

# Compute the average pile up hits
if integration_time > 1:
    for l, h in h_pu_x_layer.items():
        h_avg_pu_x_layer.Fill(l, h.GetMean())
        #h_avg_pu_x_layer.SetBinError(l, h.GetMeanError()) #TODO: set the right error

# Draw the histograms
if draw_maps:
    draw_map(hist_zr, "z [mm]", "r [mm]", sample_name+"_map_zr_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_map(hist_xy, "x [mm]", "y [mm]", sample_name+"_map_xy_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_map(hist_zphi, "z [mm]", "#phi [rad]", sample_name+"_map_zphi_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_map(h_hit_t_x_layer, "timing [ns]", "layer number", sample_name+"_map_timing_"+str(n_events)+"evt_"+sub_detector, collection)

    if not skip_plot_per_layer:
        for l, h in h_zphi_density_vs_layer.items():
            draw_map(h, "z [mm]", "#phi [rad]", sample_name+f"_map_zphi_layer{l}_"+str(n_events)+"evt_"+sub_detector, collection)

if draw_hists: 
    draw_hist(h_hit_E_MeV, "Deposited energy [MeV]", "Hits / events",  sample_name+"_hit_E_MeV_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_hist(h_hit_E_keV, "Deposited energy [keV]", "Hits / events",  sample_name+"_hit_E_keV_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_hist(h_hit_t, "Timing [ns]", "Hits / events",  sample_name+"_hit_t_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_hist(h_hit_t_corr, "Timing - TOF [ns]", "Hits / events",  sample_name+"_hit_t_corr_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_hist(h_avg_hits_x_layer, "Layer number", "Hits / events",  sample_name+"_hits_x_layer_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_hist(h_avg_occ_x_layer, "Layer number", "Average occupancy [%]",  sample_name+"_avg_occ_x_layer_"+str(n_events)+"evt_"+sub_detector, collection, log_y=False)
    draw_hist(h_occ, "Occupancy [%]", "Entries / events",  sample_name+f"_occ_tot_"+str(n_events)+"evt_"+sub_detector, collection, log_x=True)

    if energy_per_layer:
        for l, h in h_hit_E_MeV_x_layer.items():
            draw_hist(h, "Deposited energy [MeV]", "Hits per event",  sample_name+f"_hit_E_MeV_x_layer{l}_"+str(n_events)+"evt_"+sub_detector, collection)

    if not skip_plot_per_layer:
        for l, h in h_occ_x_layer.items():
            draw_hist(h, "Occupancy [%]", "Entries / events",  sample_name+f"_occ_x_layer{l}_"+str(n_events)+"evt_"+sub_detector, collection, log_x=True)

    draw_hist(h_particle_E, "MC particle energy [GeV]", "Particle / events",  sample_name+"_particle_E_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_hist(h_particle_pt, "MC particle p_{T} [GeV]", "Particle / events",  sample_name+"_particle_pt_"+str(n_events)+"evt_"+sub_detector, collection)
    draw_hist(h_particle_eta, "MC particle #eta", "Particle / events",  sample_name+"_particle_eta_"+str(n_events)+"evt_"+sub_detector, collection)

    h_particle_ID.GetXaxis().LabelsOption("v>")  # vertical labels, sorted by decreasing values
    draw_hist(h_particle_ID, "MC particle PDG ID", "Hits / events",  sample_name+"_particle_ID_"+str(n_events)+"evt_"+sub_detector, collection)

    if not skip_plot_per_layer:
        for l in h_z_density_vs_layer_mm.keys():
            draw_hist(h_z_density_vs_layer_mm[l], "z [mm]","Hits/event",  sample_name+f"_zDensity_layer{l}_"+str(n_events)+"evt_"+sub_detector, collection)
            draw_hist(h_phi_density_vs_layer[l],  "phi","Hits/event",  sample_name+f"_phiDensity_layer{l}_"+str(n_events)+"evt_"+sub_detector, collection)
    
    if integration_time > 1:
         draw_hist(h_tot_pu, "Number of pileup hits", "Entries",  sample_name+"_tot_pu_"+str(n_events)+"evt_"+sub_detector, collection)
         draw_hist(h_avg_pu_x_layer, "Layer", "Average pileup hits",  sample_name+"_avg_pu_x_layer_"+str(n_events)+"evt_"+sub_detector, collection)

print("Writing histograms...")
# Write the histograms to the output file
output_file_name = f"{sample_name}_{output_file_name}_{n_events}evt_{sub_detector}.root"
with ROOT.TFile(output_file_name,"RECREATE") as f:
    for h in histograms:
        if(debug>1): print("Writing histo:", h.GetName())
        h.Write()

print("Histograms saved in:", output_file_name)


if processed_events != n_events:
    print(f"ERROR: processed events ({processed_events}) != total expected events ({n_events}) !!!")
    print("       ===> Histograms normalization is WRONG!")
