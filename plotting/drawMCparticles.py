#!/usr/bin/env python

"""
Kinematics and production vertices of the MC particles of an event, independently of the hits they create
(for the particles that create hits, see drawhits.py).

Works on simulated files (MCParticles) and, with --mcCollection, also on generator-level edm4hep files
(e.g. --mcCollection Pairs for the incoherent pair creation samples), to inspect the BIB particles before ddsim.
"""

import argparse
import math
from collections import Counter, defaultdict

import numpy as np
from podio import root_io
import ROOT

from helpers import path_to_list, sorted_n_files, print_header
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

parser = argparse.ArgumentParser(description='drawmcparticles.py',
        epilog='Example:\ndrawmcparticles.py -i <sim_file.root> -e 100 --sample myparticles\n'
               'drawmcparticles.py -i <generator_level_file.root> --mcCollection Pairs -e 100',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--infilePath',
                  type=str, required=True,
                  help='path to input file or input directory')
parser.add_argument('-o', '--outputFile',
                  type=str, default='particles',
                  help='Name of output root file, the sample name will be added as prefix')
parser.add_argument('-n', '--numberOfFiles',
                  type=int, default=1,
                  help='number of files to consider (1 event per file), put -1 to take all files')
parser.add_argument('-e', '--numberOfEvents',
                  type=int, default=1,
                  help='number of events per file to consider, put -1 to take all events')
parser.add_argument('-t', '--tree',
                  type=str, default='events',
                  help='name of the tree in the root file')
parser.add_argument('--mcCollection',
                  type=str, default='MCParticles',
                  help='name of the MC particle collection, e.g. Pairs for the generator-level IPC files')
parser.add_argument('--sample',
                  type=str, default='ipc',
                  help='sample name to save plots')

options = parser.parse_args()

n_files = options.numberOfFiles
events_per_file = options.numberOfEvents
path = options.infilePath
output_file_name = options.outputFile
tree_name = options.tree
mc_collection = options.mcCollection
sample_name = options.sample


#######################################
# parse the input path and/or files

print_header("Importing and reading files", oneline=True)

# list n_files in the path directory ending in .root to consider for plotting
# use sorted list to make sure to always take the same files
list_files = path_to_list(path)
list_input_files = sorted_n_files(list_files, n_files)

print("# Initializing the PODIO reader...")
podio_reader = root_io.Reader(list_input_files)

nEvents_fullFileList = len(podio_reader.get(tree_name))
print(f"Number of events in the full file list: {nEvents_fullFileList}")

# number of events to process: the requested ones, capped to the available ones
n_events = nEvents_fullFileList
if events_per_file > 0:
    n_events = min(events_per_file * len(list_input_files), nEvents_fullFileList)

fill_weight = 1. / (n_events)


#######################################
# prepare histograms

# Binning estimated based on ALLEGRO size: ~12m long, ~ 5m radius
r_binning = [500, 0, 5000]       # radius
sr_binning = [500, -5000, 5000]  # signed radius
phi_binning = [70, -3.5, 3.5]
z_binning = [120, -6000, 6000]

# Particle categories, with a short description for the histogram subtitles
categories = {
    "all": "all MC particles",
    "gen": "generator-level MC particles (not created in the simulation)",
    "sim": "MC particles created in the simulation",
    "calorimeter": "MC particles decaying in the calorimeter",
    "tracker": "MC particles decaying in the tracker",
}

# Kinematic variables: name -> (axis title, scale factor w.r.t. the value in GeV, binning)
variables = {
    "E_GeV":  ("MC particle energy [GeV]", 1., [500, 0, 50]),
    "E_MeV":  ("MC particle energy [MeV]", 1e3, [1000, 0, 1000]),
    "E_keV":  ("MC particle energy [keV]", 1e6, [1000, 0, 1000]),
    "pt_GeV": ("MC particle p_{T} [GeV]", 1., [500, 0, 50]),
    "pt_MeV": ("MC particle p_{T} [MeV]", 1e3, [1000, 0, 1000]),
    "pt_keV": ("MC particle p_{T} [keV]", 1e6, [1000, 0, 1000]),
    "eta":    ("MC particle #eta", 1., [200, -10, 10]),
    "phi":    ("MC particle #phi", 1., phi_binning),
}

histograms = defaultdict(list)  # output folder -> histograms
h_particles = defaultdict(dict)

# Standard histograms for each particle category: kinematics + pdgID
for c, desc in categories.items():
    for var, (title, scale, binning) in variables.items():
        name = f"h_{c}_{var}"
        h_particles[c][var] = ROOT.TH1D(name, f"#splitline{{{name}}}{{Kinematics of the {desc}}}; {title}; particles / event", *binning)

    # ID histogram - use alphanumeric labels, form https://root.cern/doc/master/hist004__TH1__labels_8C.html
    name = f"h_{c}_pdgID"
    h_particles[c]["ID"] = ROOT.TH1D(name, f"#splitline{{{name}}}{{PDG ID of the {desc}}}; PDG ID; particles / event", 1, 0, 1)
    h_particles[c]["ID"].SetCanExtend(ROOT.TH1.kAllAxes)   # Allow both axes to extend past the initial range

    histograms[c] += list(h_particles[c].values())


# Production vertex of the particles created in the simulation
vtx_note = "Production vertex of the MC particles created in the simulation"
def vtx_hist1d(name, what, unit_title, binning):
    return ROOT.TH1D(name, f"#splitline{{{name}}}{{{vtx_note}, {what}}}; {unit_title}; particles / event", *binning)
def vtx_hist2d(name, what, titlex, titley, binning_x, binning_y):
    return ROOT.TH2D(name, f"#splitline{{{name}}}{{{vtx_note}, {what}}}; {titlex}; {titley}; particles / event", *binning_x, *binning_y)

histograms["sim_vertex"] += [
    h_sim_vertex_r_mm := vtx_hist1d("sim_vertex_r_mm", "radius", "origin r [mm]", r_binning),
    h_sim_vertex_r_mm_tracker := vtx_hist1d("sim_vertex_r_mm_tracker", "radius, tracker region", "origin r [mm]", [200, 0, 2000]),
    h_sim_vertex_r_mm_vtx_det := vtx_hist1d("sim_vertex_r_mm_vtx_det", "radius, vertex detector region", "origin r [mm]", [400, 0, 400]),

    h_sim_vertex_z_mm := vtx_hist1d("sim_vertex_z_mm", "z", "origin z [mm]", z_binning),
    h_sim_vertex_z_mm_tracker := vtx_hist1d("sim_vertex_z_mm_tracker", "z, tracker region", "origin z [mm]", [600, -3000, 3000]),

    h_sim_vertex_phi := vtx_hist1d("sim_vertex_phi", "azimuthal angle", "origin #phi", phi_binning),

    m_sim_vertex_xy_mm := vtx_hist2d("sim_vertex_xy", "x-y", "origin x [mm]", "origin y [mm]", sr_binning, sr_binning),
    m_sim_vertex_xy_mm_tracker := vtx_hist2d("sim_vertex_xy_tracker", "x-y, tracker region", "origin x [mm]", "origin y [mm]", [200, -2000, 2000], [200, -2000, 2000]),
    m_sim_vertex_xy_mm_vtx_det := vtx_hist2d("sim_vertex_xy_vtx_det", "x-y, vertex detector region", "origin x [mm]", "origin y [mm]", [400, -400, 400], [400, -400, 400]),

    m_sim_vertex_zr_mm := vtx_hist2d("sim_vertex_zr", "z-r", "origin z [mm]", "origin r [mm]", z_binning, r_binning),
    m_sim_vertex_zr_mm_tracker := vtx_hist2d("sim_vertex_zr_tracker", "z-r, tracker region", "origin z [mm]", "origin r [mm]", [600, -3000, 3000], [200, 0, 2000]),
    m_sim_vertex_zr_mm_vtx_det := vtx_hist2d("sim_vertex_zr_vtx_det", "z-r, vertex detector region", "origin z [mm]", "origin r [mm]", [500, -1000, 1000], [400, 0, 400]),

    m_sim_vertex_zphi_mm := vtx_hist2d("sim_vertex_zphi", "z-#phi", "origin z [mm]", "origin #phi", z_binning, phi_binning),
]


#######################################
# functions

def fill_arrays(h, values, weights):
    """Fill a 1D histogram with all the values of an event at once (much faster than one Fill call per value)"""
    h.FillN(len(values), values, weights)

def fill_arrays_2d(h, xs, ys, weights):
    h.FillN(len(xs), xs, ys, weights)

def as_array(values):
    return np.ascontiguousarray(values, dtype=np.float64)


#######################################
# run the event loop

print_header("Processing events", oneline=True)

processed_events = 0
processed_particles = 0
pdg_counts = defaultdict(Counter)  # category -> {PDG ID label: number of particles}
progress_step = max(1, n_events // 50) if n_events > 0 else 1  # update the progress every 2%

for i, event in enumerate(podio_reader.get(tree_name)):

    if i >= n_events and events_per_file > 0:
        break

    if (i % progress_step == 0) or (i + 1 == n_events):
        print(f"\rprocessing event {min(i + 1, n_events)}/{n_events} ({int(100 * min(i + 1, n_events) / n_events):3d}%)", end="", flush=True)

    try:
        particles = event.get(mc_collection)
    except KeyError:
        raise SystemExit(f"\nERROR: collection '{mc_collection}' not found in the event, use --mcCollection to choose another one")

    # values of this event per category: E, pt, eta, phi
    values = {c: ([], [], [], []) for c in categories}
    vertex = ([], [], [])  # x, y, z of the particles created in the simulation

    # Loop over the particles
    # Particles doxy:
    # https://edm4hep.web.cern.ch/classedm4hep_1_1_m_c_particle.html
    for particle in particles:
        # kinematics computed directly (equivalent to a PxPyPzM4D lorentz vector, but much faster in python)
        mom = particle.getMomentum()
        px, py, pz = mom.x, mom.y, mom.z
        pt = math.hypot(px, py)
        energy = math.sqrt(px * px + py * py + pz * pz + particle.getMass() ** 2)
        eta = math.asinh(pz / pt) if pt > 0 else math.copysign(math.inf, pz)
        phi = math.atan2(py, px)
        pdg_label = str(particle.getPDG())

        cats = ["all"]
        if particle.isCreatedInSimulation():
            cats.append("sim")
            vtx = particle.getVertex()
            vertex[0].append(vtx.x)
            vertex[1].append(vtx.y)
            vertex[2].append(vtx.z)
        else:
            cats.append("gen")
        if particle.isDecayedInCalorimeter():
            cats.append("calorimeter")
        if particle.isDecayedInTracker():
            cats.append("tracker")

        for c in cats:
            for lst, v in zip(values[c], (energy, pt, eta, phi)):
                lst.append(v)
            pdg_counts[c][pdg_label] += 1

        processed_particles += 1
    # <--- end of particles loop

    # Fill the histograms of this event
    for c, (energies, pts, etas, phis) in values.items():
        if not energies:
            continue
        w = as_array(np.full(len(energies), fill_weight))
        arrays = {"E": as_array(energies), "pt": as_array(pts), "eta": as_array(etas), "phi": as_array(phis)}
        for var, (title, scale, binning) in variables.items():
            quantity = var.split("_")[0]
            fill_arrays(h_particles[c][var], as_array(arrays[quantity] * scale) if scale != 1. else arrays[quantity], w)

    if vertex[0]:
        x, y, z = as_array(vertex[0]), as_array(vertex[1]), as_array(vertex[2])
        r = np.hypot(x, y)
        phi_v = np.arctan2(y, x)
        w = as_array(np.full(len(x), fill_weight))
        for h in (h_sim_vertex_r_mm, h_sim_vertex_r_mm_tracker, h_sim_vertex_r_mm_vtx_det):
            fill_arrays(h, r, w)
        for h in (h_sim_vertex_z_mm, h_sim_vertex_z_mm_tracker):
            fill_arrays(h, z, w)
        fill_arrays(h_sim_vertex_phi, phi_v, w)
        for h in (m_sim_vertex_xy_mm, m_sim_vertex_xy_mm_tracker, m_sim_vertex_xy_mm_vtx_det):
            fill_arrays_2d(h, x, y, w)
        for h in (m_sim_vertex_zr_mm, m_sim_vertex_zr_mm_tracker, m_sim_vertex_zr_mm_vtx_det):
            fill_arrays_2d(h, z, r, w)
        fill_arrays_2d(m_sim_vertex_zphi_mm, z, phi_v, w)

    processed_events += 1
# <--- end of event loop
if n_events > 0:
    print()

# PDG ID histograms have text labels, fill them once with the counts of the whole run
for c, counts in pdg_counts.items():
    for label, n in counts.items():
        h_particles[c]["ID"].Fill(label, n * fill_weight)

print_header("Done: statistics", oneline=True)
print(" - Processed events:", processed_events)
print(" - Processed particles:", processed_particles,
      f"(avg. per event: {processed_particles/max(processed_events, 1):.1f})")


#######################################
# Draw histograms

print_header("Drawing histograms", oneline=True)

for c, h_dict in h_particles.items():
    for var, (title, scale, binning) in variables.items():
        quantity = var.split("_")[0]
        draw_hist(h_dict[var], title, "Particles per event", f"{sample_name}_{c}_particle_{var}_{n_events}evt", c+"_particles")

    h_dict["ID"].GetXaxis().LabelsOption("v>")  # vertical labels, sorted by decreasing values
    draw_hist(h_dict["ID"], "MC particle PDG ID", "Particles per event", f"{sample_name}_{c}_particle_ID_{n_events}evt", c+"_particles")

# Draw the vertex origin position of sim particles
draw_hist(h_sim_vertex_r_mm,  "Origin r [mm]", "Particles per event", f"{sample_name}_sim_vertex_r_mm")
draw_hist(h_sim_vertex_r_mm_tracker,  "Origin r [mm]", "Particles per event", f"{sample_name}_sim_vertex_r_mm_tracker")
draw_hist(h_sim_vertex_r_mm_vtx_det,  "Origin r [mm]", "Particles per event", f"{sample_name}_sim_vertex_r_mm_vtx_det")
draw_hist(h_sim_vertex_z_mm, "Origin z [mm]", "Particles per event", f"{sample_name}_sim_vertex_z_mm")
draw_hist(h_sim_vertex_z_mm_tracker, "Origin z [mm]", "Particles per event", f"{sample_name}_sim_vertex_z_mm_tracker")
draw_hist(h_sim_vertex_phi, "Origin #phi", "Particles per event", f"{sample_name}_sim_vertex_phi")

draw_map(m_sim_vertex_xy_mm,"Origin x [mm]","Origin y [mm]",  f"{sample_name}_sim_vertex_xy_mm")
draw_map(m_sim_vertex_xy_mm_tracker,"Origin x [mm]","Origin y [mm]",  f"{sample_name}_sim_vertex_xy_mm_tracker")
draw_map(m_sim_vertex_xy_mm_vtx_det,"Origin x [mm]","Origin y [mm]",  f"{sample_name}_sim_vertex_xy_mm_vtx_det")
draw_map(m_sim_vertex_zr_mm, "Origin z [mm]","Origin r [mm]",  f"{sample_name}_sim_vertex_zr_mm")
draw_map(m_sim_vertex_zr_mm_tracker, "Origin z [mm]","Origin r [mm]",  f"{sample_name}_sim_vertex_zr_mm_tracker")
draw_map(m_sim_vertex_zr_mm_vtx_det, "Origin z [mm]","Origin r [mm]",  f"{sample_name}_sim_vertex_zr_mm_vtx_det")
draw_map(m_sim_vertex_zphi_mm, "Origin z [mm]","Origin #phi",  f"{sample_name}_sim_vertex_zphi_mm")


#######################################
# Write histograms to output file, one folder per particle category

output_file_name = f"{sample_name}_{output_file_name}.root"
with ROOT.TFile(output_file_name, "RECREATE") as f:
    for folder, hists in histograms.items():
        d = f.mkdir(folder)
        d.cd()
        for h in hists:
            h.Write()

print("Histograms saved in:", output_file_name)

if processed_events != n_events:
    print(f"ERROR: processed events ({processed_events}) != total expected events ({n_events}) !!!")
    print("       ===> Histograms normalization is WRONG!")
