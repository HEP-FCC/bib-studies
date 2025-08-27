"""
    Plotting and visualization tools
"""
import ROOT


def setup_root_style(stat_box=True):
    """
    Setup the general style for ROOT plots
    """
    gs = ROOT.gStyle

    if stat_box:
        gs.SetOptStat(111111)
        gs.SetOptFit(111111)
    else:
        gs.SetOptStat(0000)
        gs.SetOptFit(0000)

    gs.SetPadTickX(1)
    gs.SetPadTickY(1)
    
    gs.SetPadTopMargin(0.15)
    gs.SetPadRightMargin(0.15)
    gs.SetPadBottomMargin(0.15)
    gs.SetPadLeftMargin(0.15)

    gs.SetTextFont(42)
    gs.SetTextSize(0.05)

    gs.SetTitleXOffset(1.4)
    gs.SetTitleYOffset(1.4)
    gs.SetTitleOffset(1.4,"z")

    gs.SetMarkerSize(1.2)
    gs.SetHistLineWidth(2)


def myText(x, y, color, text, size=0.04, angle=0.):
    l = ROOT.TLatex()
    l.SetNDC()
    l.SetTextColor(color)
    l.PaintLatex( x, y, angle, size, text)
    l.DrawLatex(x,y,text)


def draw_map(hist, titlex, titley, plot_title, message='', log_x = False, log_y=False):

    cm1 = ROOT.TCanvas("cm1_"+message, "cm1_"+message, 800, 600)

    hist.GetXaxis().SetTitle(titlex)
    hist.GetYaxis().SetTitle(titley)
    hist.Draw("COLZ")
    entries = hist.GetEntries()

    myText(0.20, 0.9, ROOT.kBlack, message + '  , entries = ' + str(entries))

    if log_x:
        cm1.SetLogx()
    if log_y:
        cm1.SetLogy()

    cm1.Print(plot_title+".pdf")
    #cm1.Print(plot_title+".png")


def draw_hist(hist, titlex, titley, plot_title, message='', log_x = False, log_y=False, draw_opt="HistE"):

    cm1 = ROOT.TCanvas("cm1_"+message, "cm1_"+message, 800, 600)

    hist.GetXaxis().SetTitle(titlex)
    hist.GetYaxis().SetTitle(titley)
    hist.Draw(draw_opt)
    entries = hist.GetEntries()

    myText(0.20, 0.9, ROOT.kBlack, message + '  , entries = ' + str(entries))

    if log_x:
        cm1.SetLogx()
    if log_y:
        cm1.SetLogy()

    cm1.Print(plot_title+".pdf")
    #cm1.Print(plot_title+".png")
