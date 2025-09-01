#!/usr/bin/env python

from collections import defaultdict
import math
from optparse import OptionParser

from podio import root_io
import ROOT

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
                  type=str, default='particles',
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


(options, args) = parser.parse_args()

n_files = options.numberOfFiles
events_per_file = options.numberOfEvents
path = options.infilePath
output_file_name = options.outputFile
tree_name = options.tree
sample_name = options.sample


#######################################
# parse the input path and/or files

# list n_files in the path directory ending in .root to consider for plotting
# use sorted list to make sure to always take the same files
list_files = path_to_list(path)
list_input_files = sorted_n_files(list_files, n_files)


#######################################
# functions

def fill_kinematics(h_dict, p, p4, w):
    h_dict["E_GeV"].Fill(p4.E(), w)
    h_dict["E_MeV"].Fill(p4.E() * 1e3, w)
    h_dict["E_keV"].Fill(p4.E() * 1e6, w)
    h_dict["pt_GeV"].Fill(p4.pt(), w)
    h_dict["pt_MeV"].Fill(p4.pt() * 1e3, w)
    h_dict["pt_keV"].Fill(p4.pt() * 1e6, w)
    h_dict["eta"].Fill(p4.eta(), w)
    h_dict["phi"].Fill(p4.phi(), w)
    h_dict["ID"].Fill(str(p.getPDG()), w)


#######################################
# prepare histograms

histograms = []

# Binning estimated based on ALLEGRO size: ~12m long, ~ 5m radius
r_binning = [500, 0, 5000]       # radius
sr_binning = [500, -5000, 5000]  # signed radius
phi_binning = [70, -3.5, 3.5]
z_binning = [120, -6000, 6000]

# Particle collection name
collections = [
    "all",
    "gen",
    "sim",
    "calorimeter",
    "tracker",
]

h_particles = defaultdict(dict)

# Standard histograms for each particle collection: kinematics + pdgID
for c in collections:
    h_particles[c]["E_GeV"] = ROOT.TH1D(f"h_{c}_E_GeV", f"h_{c}_E_GeV", 500, 0, 50)
    h_particles[c]["E_MeV"] = ROOT.TH1D(f"h_{c}_E_MeV", f"h_{c}_E_MeV", 1000, 0, 1000)
    h_particles[c]["E_keV"] = ROOT.TH1D(f"h_{c}_E_keV", f"h_{c}_E_keV", 1000, 0, 1000)
    h_particles[c]["pt_GeV"] = ROOT.TH1D(f"h_{c}_pt_GeV", f"h_{c}_pt_GeV", 500, 0, 50)
    h_particles[c]["pt_MeV"] = ROOT.TH1D(f"h_{c}_pt_MeV", f"h_{c}_pt_MeV", 1000, 0, 1000)
    h_particles[c]["pt_keV"] = ROOT.TH1D(f"h_{c}_pt_keV", f"h_{c}_pt_keV", 1000, 0, 1000)
    h_particles[c]["eta"] = ROOT.TH1D(f"h_{c}_eta", f"h_{c}_eta", 100, -5, 5)
    h_particles[c]["phi"] = ROOT.TH1D(f"h_{c}_phi", f"h_{c}_phi", *phi_binning)

    # ID histogram - use alphanumeric labels, form https://root.cern/doc/master/hist004__TH1__labels_8C.html
    h_particles[c]["ID"] = ROOT.TH1D(f"h_{c}_pdgID_", f"h_{c}_pdgID_", 1, 0, 1)
    h_particles[c]["ID"].SetCanExtend(ROOT.TH1.kAllAxes)   # Allow both axes to extend past the initial range

    histograms += list(h_particles[c].values())




# Vertex origin of sim particles
h_sim_vertex_r_mm = ROOT.TH1D("sim_vertex_r_mm", "sim_vertex_r_mm", *r_binning)
h_sim_vertex_r_mm_tracker = ROOT.TH1D("sim_vertex_r_mm_tracker", "sim_vertex_r_mm_tracker", 200, 0, 2000)
h_sim_vertex_r_mm_vtx_det = ROOT.TH1D("sim_vertex_r_mm_vtx_det", "sim_vertex_r_mm_vtx_det", 400, 0, 400)
h_sim_vertex_z_mm = ROOT.TH1D("sim_vertex_z_mm", "sim_vertex_z_mm", *z_binning)
h_sim_vertex_z_mm_tracker = ROOT.TH1D("sim_vertex_z_mm_tracker", "sim_vertex_z_mm_tracker", 600, -3000, 3000)
h_sim_vertex_phi = ROOT.TH1D("sim_vertex_phi", "sim_vertex_phi", *phi_binning)
histograms += [
    h_sim_vertex_r_mm, h_sim_vertex_r_mm_tracker, h_sim_vertex_r_mm_vtx_det, 
    h_sim_vertex_phi, 
    h_sim_vertex_z_mm, h_sim_vertex_z_mm_tracker ]

m_sim_vertex_xy_mm = ROOT.TH2D("sim_vertex_xy", "sim_vertex_xy", *sr_binning, *sr_binning)
m_sim_vertex_xy_mm_tracker = ROOT.TH2D("sim_vertex_xy_tracker", "sim_vertex_xy", 200, -2000, 2000, 200, -2000, 2000)
m_sim_vertex_xy_mm_vtx_det = ROOT.TH2D("sim_vertex_xy_vtx_det", "sim_vertex_xy", 400, -400, 400, 400, -400, 400)
m_sim_vertex_zr_mm = ROOT.TH2D("sim_vertex_zr", "sim_vertex_zr", *z_binning, *r_binning)
m_sim_vertex_zr_mm_tracker = ROOT.TH2D("sim_vertex_zr_tracker", "sim_vertex_zr_tracker", 600, -3000, 3000, 200, 0, 2000)
m_sim_vertex_zr_mm_vtx_det = ROOT.TH2D("sim_vertex_zr_vtx_det", "sim_vertex_zr_vtx_det", 500, -1000, 1000, 400, 0, 400)
m_sim_vertex_zphi_mm = ROOT.TH2D("sim_vertex_zphi", "sim_vertex_zphi", *z_binning, *phi_binning)

histograms += [
    m_sim_vertex_xy_mm, m_sim_vertex_xy_mm_tracker, m_sim_vertex_xy_mm_vtx_det, 
    m_sim_vertex_zr_mm, m_sim_vertex_zr_mm_tracker, m_sim_vertex_zr_mm_vtx_det, 
    m_sim_vertex_zphi_mm]


n_events = events_per_file * len(list_input_files)

fill_weight = 1. / (n_events)

#######################################
# run the event loop

print("Initializing the PODIO reader...")

podio_reader = root_io.Reader(list_input_files)

processed_events = 0
processed_particles = 0

for i,event in enumerate(podio_reader.get(tree_name)):

    if i >= n_events and events_per_file > 0:
        break

    if i % 100 == 0: # print a message every 100 processed files
        print('processing event number = ' + str(i) + ' / ' + str(n_events) )

    # Loop over the particle
    for particle in event.get("MCParticles"):
        # Particles doxy:
        # https://edm4hep.web.cern.ch/classedm4hep_1_1_m_c_particle.html

        # retrieve the four momentum
        particle_p4 = ROOT.Math.LorentzVector('ROOT::Math::PxPyPzM4D<double>')(
            particle.getMomentum().x, 
            particle.getMomentum().y, 
            particle.getMomentum().z, 
            particle.getMass())

        # Fill MC particle histograms
        fill_kinematics(h_particles["all"], particle, particle_p4, fill_weight)


        #if not particle.isBackscatter(): #  exclude particles originating from a calorimeter shower

        if particle.isCreatedInSimulation():
            # Sim particles
            fill_kinematics(h_particles["sim"], particle, particle_p4, fill_weight)

            # Access the production vertex of the particle
            vtx = particle.getVertex()
            x_mm = particle.getVertex().x
            y_mm = particle.getVertex().y
            z_mm = particle.getVertex().z
            r_mm = math.sqrt(math.pow(x_mm, 2) + math.pow(y_mm, 2))
            phi = math.acos(x_mm/r_mm) * math.copysign(1, y_mm)

            h_sim_vertex_r_mm.Fill(r_mm, fill_weight)
            h_sim_vertex_r_mm_tracker.Fill(r_mm, fill_weight)
            h_sim_vertex_r_mm_vtx_det.Fill(r_mm, fill_weight)
            h_sim_vertex_phi.Fill(phi, fill_weight)
            h_sim_vertex_z_mm.Fill(z_mm, fill_weight)
            h_sim_vertex_z_mm_tracker.Fill(z_mm, fill_weight)
            m_sim_vertex_xy_mm.Fill(x_mm, y_mm ,fill_weight)
            m_sim_vertex_xy_mm_tracker.Fill(x_mm, y_mm ,fill_weight)
            m_sim_vertex_xy_mm_vtx_det.Fill(x_mm, y_mm ,fill_weight)
            m_sim_vertex_zr_mm.Fill(z_mm, r_mm, fill_weight)
            m_sim_vertex_zr_mm_tracker.Fill(z_mm, r_mm, fill_weight)
            m_sim_vertex_zr_mm_vtx_det.Fill(z_mm, r_mm, fill_weight)
            m_sim_vertex_zphi_mm.Fill(z_mm, phi, fill_weight)

        else:
            # Generator particles
            fill_kinematics(h_particles["gen"], particle, particle_p4, fill_weight)

        # Particles decaying in the calo
        if particle.isDecayedInCalorimeter():
            fill_kinematics(h_particles["calorimeter"], particle, particle_p4, fill_weight)

        # Particles decaying in the tacker
        if particle.isDecayedInTracker():
            fill_kinematics(h_particles["tracker"], particle, particle_p4, fill_weight)

        processed_particles += 1 
    # <--- end of particles loop

    processed_events +=1
# <--- end of event loop

print("Loop done!")
print(" - Processed events:", processed_events)
print(" - Processed particles:", processed_particles, 
      f"(avg. per event: {processed_particles/processed_events})")

#######################################
# Draw histograms

for c, h_dict in h_particles.items():
    draw_hist(h_dict["E_GeV"],  "MC particle energy [GeV]", "Particles per event", f"{sample_name}_{c}_particle_E_GeV_{n_events}evt",  c+"_particles")
    draw_hist(h_dict["E_MeV"],  "MC particle energy [MeV]", "Particles per event", f"{sample_name}_{c}_particle_E_MeV_{n_events}evt",  c+"_particles")
    draw_hist(h_dict["E_keV"],  "MC particle energy [keV]", "Particles per event", f"{sample_name}_{c}_particle_E_keV_{n_events}evt",  c+"_particles")
    draw_hist(h_dict["pt_GeV"], "MC particle p_{T} [GeV]", "Particles per event",  f"{sample_name}_{c}_particle_pt_GeV_{n_events}evt", c+"_particles")
    draw_hist(h_dict["pt_MeV"], "MC particle p_{T} [MeV]", "Particles per event",  f"{sample_name}_{c}_particle_pt_MeV_{n_events}evt", c+"_particles")
    draw_hist(h_dict["pt_keV"], "MC particle p_{T} [keV]", "Particles per event",  f"{sample_name}_{c}_particle_pt_keV_{n_events}evt", c+"_particles")
    draw_hist(h_dict["eta"],    "MC particle #eta", "Particles per event",         f"{sample_name}_{c}_particle_eta_{n_events}evt",    c+"_particles")
    draw_hist(h_dict["phi"],    "MC particle #phi", "Particles per event",         f"{sample_name}_{c}_particle_phi_{n_events}evt",    c+"_particles")

    h_dict["ID"].GetXaxis().LabelsOption("v>")  # vertical labels, sorted by decreasing values
    draw_hist(h_dict["ID"], "MC particle PDG ID", "Hits / events", f"{sample_name}_{c}_particle_ID_{n_events}evt", c+"_particles")

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
# Write histograms to output file

output_file_name = f"{sample_name}_{output_file_name}.root"
with ROOT.TFile(output_file_name,"RECREATE") as f:
    for h in histograms:
        h.Write()

print("Histograms saved in:", output_file_name)

if processed_events != n_events:
    print(f"ERROR: processed events ({processed_events}) != total expected events ({n_events}) !!!")
    print("       ===> Histograms normalization is WRONG!")


