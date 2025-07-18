import ROOT
from podio import root_io
import dd4hep as dd4hepModule
from ROOT import dd4hep
import os
import math
from optparse import OptionParser

######################################
# style
ROOT.gStyle.SetOptStat(0000)
ROOT.gStyle.SetOptFit(0000)

######################################
#option parser
                  
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
binning = []
for i,x in enumerate(options.map_binning.split(',')):
    if i==0 or i==3:
        binning.append(int(x))
    else:
        binning.append(float(x))
#binning = [float(x) for x in options.map_binning.split(',')]

#######################################
#functions

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
    
    myText(0.20, 0.9, ROOT.kBlack, message)
    
    cm1.Print(plot_title+".pdf")
    cm1.Print(plot_title+".png")

#######################################

# list nfiles in the path directory ending in .root  to consider for plotting
list_input_files = []
for f in sorted(os.listdir(path)): #use sorted list to make sure to always take the same files
    if f.endswith(".root"):
        list_input_files.append( os.path.join(path, f) )
    if nfiles > 0 and len(list_input_files) >= nfiles: #if nfiles=-1 run all files
        break

    
hist_zr = ROOT.TH2D( "hist_zr_"+collection, "hist_zr_"+collection, binning[0], binning[1], binning[2], binning[3], binning[4], binning[5])

for i,input_file_path in enumerate(list_input_files):
    if i % 100 == 0: # print a message every 100 processed files
        print('processing file number = ' + str(i) + ' / ' + str(nfiles) )
    podio_reader = root_io.Reader(input_file_path)

    for event in podio_reader.get(tree_name):
        
        for hit in event.get(ROOT.TString(collection).Data()): #convert python string in a root TString
            x = hit.getPosition().x
            y = hit.getPosition().y
            z = hit.getPosition().z
            r = math.sqrt(math.pow(x,2)+math.pow(y,2))
            if x<0.:
                r = -1.*r                
            hist_zr.Fill(z, r, 1./len(list_input_files))
            
if(draw_maps):
    draw_map(hist_zr, "z [mm]", "r [mm]", sample_name+"_map_zr_"+str(nfiles)+"evt_"+collection, collection)


