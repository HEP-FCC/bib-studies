#!/usr/bin/env python

import ROOT
from podio import root_io
import os
import math
from optparse import OptionParser
from glob import glob

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
parser.add_option('-s', '--sample',
                  type=str, default='ipc',
                  help='sample name to save plots')
parser.add_option('-c', '--collection',
                  type=str, default='VertexBarrelCollection',
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
collection = options.collection
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

# Profile histograms
h_hit_E = ROOT.TH1D("h_hit_E_"+collection, "h_hit_E_"+collection, 50, 0, 50)

# Hit maps definitions
hist_zr = ROOT.TH2D("hist_zr_"+collection, "hist_zr_"+collection, binning[0], binning[1], binning[2], binning[3], binning[4], binning[5])
hist_xy = ROOT.TH2D("hist_xy_"+collection, "hist_xy_"+collection, binning[3], binning[4], binning[5], binning[3], binning[4], binning[5])
hist_zphi = ROOT.TH2D("hist_zphi_"+collection, "hist_zphi_"+collection, binning[0], binning[1], binning[2], 500, -3.5, 3.5)

fill_weight = 1. / len(list_input_files)

for i,input_file_path in enumerate(list_input_files):
    if i % 100 == 0: # print a message every 100 processed files
        print('processing file number = ' + str(i) + ' / ' + str(nfiles) )
    podio_reader = root_io.Reader(input_file_path)

    for event in podio_reader.get(tree_name):

        for hit in event.get(ROOT.TString(collection).Data()): #convert python string in a root TString
            # Hits doxy:
            # https://edm4hep.web.cern.ch/classedm4hep_1_1_mutable_sim_tracker_hit.html
            # https://edm4hep.web.cern.ch/classedm4hep_1_1_mutable_sim_calorimeter_hit.html

            x = hit.getPosition().x
            y = hit.getPosition().y
            z = hit.getPosition().z
            r = math.sqrt(math.pow(x, 2) + math.pow(y, 2))
            if x < 0.:
                r = -1. * r
            phi = math.acos(x/r) * math.copysign(1, y)

            # For Calo
            # E = hit.getEnergy() * 1e-3 # The doxy says it should be in GeV, but seems MeV?

            # For tracker
            E = hit.getEDep() * 1e-3 # The doxy says it should be in GeV, but seems MeV?

            # Fill the histograms
            h_hit_E.Fill(x, fill_weight)

            hist_zr.Fill(z, r, fill_weight)
            hist_xy.Fill(x, y, fill_weight)
            hist_zphi.Fill(z, phi, fill_weight)

if draw_maps:
    draw_map(hist_zr, "z [mm]", "r [mm]", sample_name+"_map_zr_"+str(nfiles)+"evt_"+collection, collection)
    draw_map(hist_xy, "x [mm]", "y [mm]", sample_name+"_map_xy_"+str(nfiles)+"evt_"+collection, collection)
    draw_map(hist_zphi, "z [mm]", "#phi [rad]", sample_name+"_map_zphi_"+str(nfiles)+"evt_"+collection, collection)

if draw_profiles: 
    draw_profile(h_hit_E, "Deposited Energy [GeV]", "Hits / event",  sample_name+"_hit_E_"+str(nfiles)+"evt_"+collection, collection)