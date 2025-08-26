"""
    Collection of plotting and visualization tools
"""
import ROOT


def setup_root_style():
    """
    Setup the general style for ROOT plots
    """
    gs = ROOT.gStyle

    gs.SetOptStat(111111)
    gs.SetOptFit(111111)
    
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
