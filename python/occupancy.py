#!/usr/bin/env python

import argparse
import math
from collections import Counter
from fractions import Fraction

# import dd4hep as dd4hepModule
# from ROOT import dd4hep
from podio import root_io
import ROOT

from helpers import path_to_list, sorted_n_files, load_json, layer_number_from_string, simplify_dict, is_calo, is_endcap, DetFilePath
from visualization import setup_root_style, draw_hist, myText


######################################
# option parser

def parse_args():
    parser = argparse.ArgumentParser(description='occupancy.py',
        epilog='Example:\noccupancy.py -s MuonTaggerBarrel -d $BIB_STUDIES/detectors_dicts/ALLEGRO_o1_v03_DetectorDimensions.json -n 1 -e -1 -D 0',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-D', '--debugLevel',
                    type=int, default=1,
                    help='debug level (0:quiet, 1:default, 2:verbose)')
    parser.add_argument('-i', '--infilePath',
                    type=str, default='/eos/experiment/fcc/users/j/jaeyserm/bdsim/LCC_v2short_halo/output_100093_edm4hep.root',
                    help='path to input file or input directory')
    parser.add_argument('-n', '--numberOfFiles',
                    type=int, default=1,
                    help='number of files to consider (1 event per file), put -1 to take all files')
    parser.add_argument('-e', '--numberOfEvents',
                    type=int, default=1,
                    help='number of events per file to consider, put -1 to take all events')
    parser.add_argument('--sample',
                    type=str, default='ipc',
                    help='Brief sample name for plot titles and output files (e.g. ipc, injbkg, )')
    parser.add_argument('-s', '--subDetector',
                    type=str, default='VertexBarrel',
                    help='Name of the sub detector system')
    parser.add_argument('--collection',
                    type=str, default=None,
                    help='hits collection to read, overriding the one in the detector/assumptions dict (e.g. OverlayDCHCollection)')
    parser.add_argument('-p', '--draw_hists',
                    action="store_true",
                    help='activate drawing of 1D histograms')
    parser.add_argument('--stat_box',
                    action="store_true",
                    help='draw the ROOT stat box')
    parser.add_argument('--digi',
                    action='store_true',
                    help='run on digitized hits (requires to pass also the assumptions dict with info on the digitized collection).')
    parser.add_argument('-d', '--detDictFile',
                    type=str, default='$BIB_STUDIES/detectors_dicts/ALLEGRO_o1_v03_DetectorDimensions.json',
                    help='JSON dictionary with some key detector parameters')
    parser.add_argument('-a', '--assumptions',
                    type=str, default=None,
                    help='JSON dictionary with assumptions')
    parser.add_argument('--e_cut',
                    action='store_true',
                    help='apply energy cut (requires to pass also the assumptions dict with energy threshold).')
    parser.add_argument('--per_event',
                    action='store_true',
                    help='each event is one readout window (e.g. overlaid + digitized windows): occupancy of every event, '
                         'then mean over the events with its standard error. Default: the whole input is one readout window.')
    parser.add_argument('--bx_per_event',
                    type=str, default=None,
                    help='bunch crossings represented by one input event, only to check the window length '
                         '(e.g. 1 for raw IPC, 1/2020 for SR halo, 31 for overlaid windows)')
    parser.add_argument('--window_bx',
                    type=float, default=31,
                    help='expected length of the readout window in BX (warning if the input does not match)')

    return (parser.parse_args())

options = parse_args()

debug = options.debugLevel
n_files = options.numberOfFiles
events_per_file = options.numberOfEvents
input_path = options.infilePath
sample_name = options.sample
det_file =  DetFilePath(options.detDictFile)
assumptions_path = options.assumptions
sub_detector = options.subDetector
draw_hists = options.draw_hists
stat_box = options.stat_box
digi = options.digi
e_cut = options.e_cut
per_event = options.per_event
bx_per_event = Fraction(options.bx_per_event) if options.bx_per_event else None
window_bx = options.window_bx

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
        from ALLEGRO import get_layer
    case "IDEA":
        from IDEA import get_layer
    case "CLD":
        from CLD import get_layer
    case "ILD_FCCee":
        from ILD_FCCee import get_layer
    case _:
        raise NotImplementedError(f"get_layer not implemented for detector {det_file.short}")


def spread_error(s1, s2, n):
    """
    sqrt( n/(n-1) * sum_i (x_i - mean)^2 ), from s1 = sum x_i and s2 = sum x_i^2 over the n events
    """
    if n < 2:
        return 0.
    return math.sqrt(max(0., n / (n - 1) * (s2 - s1 * s1 / n)))  # max(): rounding can make it slightly negative


def occupancy_jackknife(fired, s1, s2, n, cells):
    """
    Occupancy of the full sample (one readout window) [%], with the leave-one-event-out jackknife error.

    u_i = cells fired ONLY by event i (s1 = sum u_i, s2 = sum u_i^2). Leaving event i out removes u_i cells
    from the union, and the leave-one-out values are rescaled by n/(n-1) to the length of the full window,
    which gives sigma = sqrt( n/(n-1) * sum_i (u_i - <u>)^2 ) / cells.
    The event is the unit, so the clustering of the cells (one particle fires several cells) is taken into account.
    No fired cell: the error is the 95% CL upper limit, 3 / cells.
    """
    if fired == 0:
        return 0., 3 / cells * 100
    return fired / cells * 100, spread_error(s1, s2, n) / cells * 100


def occupancy_mean(s1, s2, n, cells):
    """
    Mean occupancy per event (= per readout window) [%], with its standard error std / sqrt(n),
    which is also the jackknife error of a mean. s1 = sum fired_i, s2 = sum fired_i^2 over the n events.
    No fired cell in any event: the error is the 95% CL upper limit on the mean, 3 / (n * cells).
    """
    if s1 == 0:
        return 0., 3 / (n * cells) * 100
    return s1 / n / cells * 100, spread_error(s1, s2, n) / n / cells * 100


def draw_occupancy_x_layer(hist, tot_occ, tot_err, plot_title, message=''):
    """
    Draw the occupancy per layer, with the total occupancy of the sub-detector as a dashed line
    """
    cm1 = ROOT.TCanvas("cm1_"+message, "cm1_"+message, 800, 600)

    hist.GetXaxis().SetTitle("Layer number")
    hist.GetYaxis().SetTitle("Occupancy [%]")
    y_max = max([hist.GetBinContent(b) + hist.GetBinError(b) for b in range(1, hist.GetNbinsX()+1)] + [tot_occ])
    hist.SetMinimum(0)
    hist.SetMaximum(1.2 * y_max if y_max > 0 else 1)
    hist.SetMarkerStyle(ROOT.kFullCircle)
    hist.Draw("E")

    x_min = hist.GetXaxis().GetXmin()
    x_max = hist.GetXaxis().GetXmax()
    tot_line = ROOT.TLine(x_min, tot_occ, x_max, tot_occ)
    tot_line.SetLineColor(ROOT.kRed)
    tot_line.SetLineStyle(ROOT.kDashed)
    tot_line.SetLineWidth(2)
    tot_line.Draw()

    myText(0.20, 0.9, ROOT.kBlack, message)
    myText(0.55, 0.9, ROOT.kRed, f"total = {tot_occ:.3g} #pm {tot_err:.2g} %")

    cm1.Print(plot_title+".pdf")
    #cm1.Print(plot_title+".png")

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

if (e_cut or digi) and assumptions == None:
    raise ValueError(f"Given the input settings (e_cut={e_cut}, digi={digi})," +
                     "an 'assumptions' JSON is needed as input (-a/--assumptions <assumptions.json>)")

# get the hits collection name
collection =  detector_dict["hitsCollection"]
if digi:
    collection = assumptions["digitized_hits"]["collection"]
if options.collection:
    collection = options.collection
if collection+"__CellIDEncoding" in metadata.parameters:
    id_encoding = metadata.get_parameter(collection+"__CellIDEncoding")
else:
    # e.g. DCH_DigiCollection: the digitizer does not store an encoding, but keeps the sim hit cellIDs
    print(f"WARNING: no CellIDEncoding stored for {collection}, using the one of {detector_dict['hitsCollection']}")
    id_encoding = metadata.get_parameter(detector_dict["hitsCollection"]+"__CellIDEncoding")
decoder = ROOT.dd4hep.BitFieldCoder(id_encoding)
#print("HERE",decoder.fieldDescription()) #get possible values

# Pixel detectors: the cellID identifies a pixel, while the JSON counts the cells per sensor
id_fields = [f.split(":")[0] for f in id_encoding.split(",")]
if "x" in id_fields and "y" in id_fields:
    print(f"WARNING: the {collection} cellIDs include x/y (pixel) fields, while the cells per layer in the",
          "detector JSON count sensors => occupancy = fired pixels / sensors!")

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

if not layer_cells:
    print("WARNING: map of cells per layer is empty.")
    n_tot_cells = 1

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

#######################################
# prepare histograms

# list of the histograms that will be saved in output ROOT file
histograms = []

layer_binning = [n_layers + 1, -0.5, n_layers + 0.5]
if is_endcap(detector_type):    # Binning from -negative to positive layer numbers, with a dummy 0 layer in the middle
    max_l = int(n_layers / 2) + 0.5
    layer_binning = [n_layers + 1, -max_l, +max_l]


# Occupancy, defined as fraction of cells that fired
# Requires that the total cells are correctly computed and stored in detector_dict
# Default: occupancy of the full sample (all events share the same readout window), one value per layer
# --per_event: mean over the events (each event is one readout window) of the occupancy of every event
fired_title = "avg. fired cells / window" if per_event else "fired cells"
histograms += [
    h_fired_cells_x_layer := ROOT.TH1D("h_fired_cells_x_layer_"+collection, "h_fired_cells_x_layer_"+collection+";layer index;"+fired_title, *layer_binning),
    h_occ_x_layer := ROOT.TH1D("h_occ_x_layer_"+collection, "h_occ_x_layer_"+collection+";layer index;fired cells / total layer cells [%]", *layer_binning),
    h_occ := ROOT.TH1D("h_occ_tot_"+collection, "h_occ_tot_"+collection+";;total fired cells / total cells [%]", 1, 0, 1),
]


n_events = nEvents_fullFileList  # already the total over all input files
if events_per_file != -1:
    n_events = min(events_per_file * len(list_input_files), nEvents_fullFileList)

# Length of the readout window: the whole input (default) or one event (--per_event)
if bx_per_event is None:
    print("WARNING: --bx_per_event not given, the length of the readout window in BX is not checked")
else:
    input_window_bx = float(bx_per_event) if per_event else float(n_events * bx_per_event)
    print(f">>> Readout window: {input_window_bx:.4g} BX ({'1 event' if per_event else f'{n_events} events'}"
          f" x {float(bx_per_event):.4g} BX/event), expected {window_bx:g} BX")
    if not math.isclose(input_window_bx, window_bx, rel_tol=1e-3):
        print(f"WARNING: the readout window is {input_window_bx:.4g} BX instead of {window_bx:g} BX:",
              "the occupancy is", "UNDER" if input_window_bx < window_bx else "OVER", "estimated",
              "(bias of the window length, not included in the error)")

#######################################
# run the event loop

# Counters for the event / hit loops
processed_events = 0
processed_hits = 0

# Place these data structures outside the loop to avoid double-counting across different events
cell_n_events = Counter()  # cellID -> number of events in which it fired (the keys are the fired cells of the run)
cell_layer = {}            # cellID -> layer
event_cells = []           # fired cells of every event (default mode, to find the cells fired by one event only)
# --per_event: sums over the events of the fired cells per layer / in total, and of their squares
sum_fired_x_layer = {l: 0 for l in layer_cells.keys()}
sum2_fired_x_layer = {l: 0 for l in layer_cells.keys()}
sum_fired = sum2_fired = 0
fired_per_event = []       # total fired cells of every event

for i, event in enumerate(podio_reader.get('events')):

    if i >= n_events and events_per_file > 0:
        break

    if i % 100 == 0: # print a message every 100 processed files
        print('processing event number = ' + str(i) + ' / ' + str(n_events) )

    fired_cells = set()  # cellIDs that fired in this event (set => unique entries)

    # Loop over the hits
    for hit in event.get(collection):
        if(debug>1): print("Processing hit:", hit)
        # Hits doxy:
        # https://edm4hep.web.cern.ch/classedm4hep_1_1_mutable_sim_tracker_hit.html
        # https://edm4hep.web.cern.ch/classedm4hep_1_1_mutable_sim_calorimeter_hit.html
        # https://edm4hep.web.cern.ch/classedm4hep_1_1_m_c_particle.html

        is_calo_hit = is_calo(detector_type)

        cell_id = hit.getCellID()

        #TODO: Handle better cell_id info
        layer_n = get_layer(cell_id, decoder, sub_detector, detector_type)

        # Energy threshold, per layer or the same for all the layers
        E_hit_thr = E_thr_MeV[layer_n] if isinstance(E_thr_MeV, dict) else E_thr_MeV

        # Deposited energy, converted to MeV
        E_hit = 0
        if is_calo_hit:
            E_hit = hit.getEnergy() * 1e3  # For calorimeter hits
        else:
            E_hit = hit.getEDep() * 1e3  # For tracker hits

        # Apply cut deposited energy, before flagging the cell as fired
        if E_hit >= E_hit_thr or (not e_cut):

            # For layer percentage occupancy, count cells only once per event
            fired_cells.add(cell_id)
            cell_layer[cell_id] = layer_n

        processed_hits += 1
    # <--- end of the hits loop

    if per_event:
        fired_x_layer = Counter(cell_layer[c] for c in fired_cells)
        for l in layer_cells.keys():
            sum_fired_x_layer[l] += fired_x_layer[l]
            sum2_fired_x_layer[l] += fired_x_layer[l]**2
        sum_fired += len(fired_cells)
        sum2_fired += len(fired_cells)**2
        fired_per_event.append(len(fired_cells))
    else:
        cell_n_events.update(fired_cells)
        event_cells.append(fired_cells)

    processed_events +=1
# <--- end of the event loop

# Compute the per-layer and total occupancy of the full sample
print("Loop done!")
print(" - Processed events:", processed_events)
print(" - Processed hits:", processed_hits,
      f"(avg. per event: {processed_hits/processed_events})")

n = processed_events
unit = "readout windows (events)" if per_event else "events"
if n < 2:
    print(f"WARNING: only {n} {unit}: the error cannot be computed and is set to 0")
elif n < 10:
    print(f"WARNING: only {n} {unit}: the error is unreliable below ~10 (jackknife / standard error)")

if per_event:
    # Mean over the readout windows, error = std / sqrt(n)
    # (s1, s2) = sums over the events of the fired cells and of their squares
    sums_x_layer = {l: (sum_fired_x_layer[l], sum2_fired_x_layer[l]) for l in layer_cells.keys()}
    sums_tot = (sum_fired, sum2_fired)
    layer_occupancy = lambda l: occupancy_mean(*sums_x_layer[l], n, layer_cells[l])
    tot_occupancy, tot_occupancy_err = occupancy_mean(*sums_tot, n, n_tot_cells)
    layer_fired = lambda l: sums_x_layer[l][0] / n
else:
    # Union of all the events (one readout window), error = leave-one-event-out jackknife
    # (s1, s2) = sums over the events of the cells fired ONLY by that event, and of their squares
    fired_x_layer = Counter(cell_layer[c] for c in cell_n_events)
    sum_unique_x_layer, sum2_unique_x_layer = Counter(), Counter()
    sum_unique = sum2_unique = 0
    for cells in event_cells:
        unique = [c for c in cells if cell_n_events[c] == 1]
        for l, u in Counter(cell_layer[c] for c in unique).items():
            sum_unique_x_layer[l] += u
            sum2_unique_x_layer[l] += u**2
        sum_unique += len(unique)
        sum2_unique += len(unique)**2
    layer_occupancy = lambda l: occupancy_jackknife(fired_x_layer[l], sum_unique_x_layer[l], sum2_unique_x_layer[l], n, layer_cells[l])
    tot_occupancy, tot_occupancy_err = occupancy_jackknife(len(cell_n_events), sum_unique, sum2_unique, n, n_tot_cells)
    layer_fired = lambda l: fired_x_layer[l]

if layer_cells:
    if debug>0: print(f"{'layer':>6s}  {'avg. fired / cells' if per_event else 'fired / cells':>18s}  occupancy [%]")
    for l in sorted(layer_cells.keys()):
        layer_occ, layer_err = layer_occupancy(l)
        fc = layer_fired(l)
        h_fired_cells_x_layer.SetBinContent(h_fired_cells_x_layer.FindBin(l), fc)
        h_fired_cells_x_layer.SetBinError(h_fired_cells_x_layer.FindBin(l), layer_err / 100 * layer_cells[l])
        h_occ_x_layer.SetBinContent(h_occ_x_layer.FindBin(l), layer_occ)
        h_occ_x_layer.SetBinError(h_occ_x_layer.FindBin(l), layer_err)
        if debug>0: print(f"{l:>6d}  {fc:>9.6g} / {layer_cells[l]:<8g}  {layer_occ:.4g} +- {layer_err:.2g}")
    h_fired_cells_x_layer.SetEntries(len(layer_cells))
    h_occ_x_layer.SetEntries(len(layer_cells))

h_occ.SetBinContent(1, tot_occupancy)
h_occ.SetBinError(1, tot_occupancy_err)
h_occ.SetEntries(1)

if per_event:
    print(f" - Mean fired cells per readout window: {sum_fired / n:.6g} / {n_tot_cells}",
          f"=> mean total occupancy: {tot_occupancy:.4g} +- {tot_occupancy_err:.2g} % (standard error over {n} windows)")
    # Distribution of the total occupancy of the readout windows
    occ_per_event = [f / n_tot_cells * 100 for f in fired_per_event]
    occ_max = max(occ_per_event) if occ_per_event and max(occ_per_event) > 0 else 1
    histograms += [ h_occ_windows := ROOT.TH1D("h_occ_windows_"+collection, "h_occ_windows_"+collection+";total fired cells / total cells [%];readout windows", 50, 0, 1.2 * occ_max) ]
    for o in occ_per_event:
        h_occ_windows.Fill(o)
else:
    print(" - Number of unique cells fired:", len(cell_n_events), "/", n_tot_cells,
          f"=> total occupancy: {tot_occupancy:.4g} +- {tot_occupancy_err:.2g} % (jackknife over {n} events)")

mode_tag = "_perEvent" if per_event else ""
if draw_hists:
    draw_hist(h_fired_cells_x_layer, "Layer number", fired_title.capitalize(),  sample_name+"_fired_cells_x_layer_"+str(n_events)+"evt"+mode_tag+"_"+sub_detector, collection)
    draw_occupancy_x_layer(h_occ_x_layer, tot_occupancy, tot_occupancy_err, sample_name+"_occ_x_layer_"+str(n_events)+"evt"+mode_tag+"_"+sub_detector, collection)
    if per_event:
        draw_hist(h_occ_windows, "Occupancy [%]", "Readout windows",  sample_name+"_occ_windows_"+str(n_events)+"evt_"+sub_detector, collection)

# Open ROOT file and create directory structure for organized output
output_file_name = f"{sample_name}_{n_events}evt{mode_tag}_{sub_detector}_{suffix_from_input}.root"
output_file = ROOT.TFile(output_file_name, "RECREATE")
dir_per_layer = output_file.mkdir("per_layer")


# Categorize and write histograms to appropriate directories
for h in histograms:
    if(debug>1): print("Writing histo:", h.GetName())

    hname = h.GetName()

    # Categorize and write to appropriate directory
    if ("_layer" in hname) or ("_x_layer" in hname):
        dir_per_layer.cd()
    else:
        output_file.cd()  # root directory for general histograms
    h.Write()

output_file.Close()
print("Histograms saved in:", output_file_name)


if processed_events != n_events:
    print(f"ERROR: processed events ({processed_events}) != total expected events ({n_events}) !!!")
    print("       ===> Not all the requested events were processed, the occupancy is underestimated!")




# Functionalities:
#
# Parse N input files with M events per file
# Parse subdetector collecion (DCH)
# Parse output file name
#
# Iterate over every event of every file, and collect the fired cells of every event (above the energy threshold, if any)
# Default: all the events are in the same readout window, so each cell counts once in the whole sample
#   occupancy = unique fired cells / total cells (no per-event normalization), total and per layer
#   error = leave-one-event-out jackknife, from the cells fired by one event only
# --per_event: every event is one readout window (e.g. overlaid + digitized windows)
#   occupancy = mean over the events of fired cells / total cells, error = std / sqrt(n_events)
# Optional (--bx_per_event): warn if the readout window is not --window_bx BX long (bias, not in the error)
#
# Optional: Plot the occupancy (and fired cells) vs layer number, with the total occupancy as a dashed line
#           (--per_event: also the distribution of the total occupancy of the windows)
#
#