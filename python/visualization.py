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


def draw_overlay(hists, titlex, titley, plot_title, message='', log_x=False, log_y=False, legend_title='', save_pdf=True):
    """
    Overlay several 1D histograms (dict {legend label: histogram}) on one canvas, e.g. one per layer.
    The histograms are cloned, so the originals are left untouched. Empty histograms are skipped.
    The canvas is printed to <plot_title>.pdf if save_pdf, and returned (None if there is nothing to draw),
    so that it can be saved in a ROOT file.

    """
    colors = [ROOT.kBlack, ROOT.kRed+1, ROOT.kBlue+1, ROOT.kGreen+2, ROOT.kMagenta+1, ROOT.kOrange+7,
              ROOT.kCyan+1, ROOT.kViolet+1, ROOT.kPink+1, ROOT.kTeal+2, ROOT.kGray+2]
    hists = {label: h for label, h in hists.items() if h.GetEntries() > 0}
    if not hists:
        print(f"WARNING: no filled histogram to overlay for {plot_title}, plot skipped")
        return None

    cm1 = ROOT.TCanvas("c_"+plot_title, "c_"+plot_title, 800, 600)
    if log_x:
        cm1.SetLogx()
    if log_y:
        cm1.SetLogy()
    legend = ROOT.TLegend(0.35, 0.67, 0.84, 0.83)  # inside the frame
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.03)
    if legend_title:
        legend.SetHeader(legend_title)
    if len(hists) > 4:
        legend.SetNColumns(2 if len(hists) <= 12 else 3)

    clones = []
    for i, (label, h) in enumerate(hists.items()):
        c = h.Clone(f"{h.GetName()}_overlay")
        c.SetDirectory(0)
        c.SetTitle("")  # the title of the first histogram would be drawn, but it only describes that one
        c.SetLineColor(colors[i % len(colors)])
        c.SetLineStyle(1 + i // len(colors))  # change the style once the colors are used up
        c.SetLineWidth(2)
        ROOT.SetOwnership(c, False)  # the canvas is returned and may outlive this function, so the histograms must too
        clones.append(c)
        legend.AddEntry(c, label, "l")

    # common vertical range, so that no histogram is cut
    ymax = max(c.GetMaximum() for c in clones)
    positive = [c.GetBinContent(b) for c in clones for b in range(1, c.GetNbinsX() + 1) if c.GetBinContent(b) > 0]
    ymin = min(positive) if positive else 0
    # zoom the x axis on the range where at least one histogram is filled
    filled = [b for c in clones for b in range(1, c.GetNbinsX() + 1) if c.GetBinContent(b) > 0]
    if filled and not log_x:
        xlo = min(clones[0].GetXaxis().GetBinLowEdge(min(filled)), clones[0].GetXaxis().GetBinLowEdge(max(filled)))
        xhi = clones[0].GetXaxis().GetBinUpEdge(max(filled))
        margin = 0.05 * (xhi - xlo)
        clones[0].GetXaxis().SetRangeUser(max(xlo - margin, clones[0].GetXaxis().GetXmin()), min(xhi + margin, clones[0].GetXaxis().GetXmax()))

    y_lo = 0.5 * ymin if log_y else 0
    y_hi = 30 * ymax if log_y else 1.6 * ymax  # headroom for the legend
    for i, c in enumerate(clones):
        if i == 0:
            c.GetXaxis().SetTitle(titlex)
            c.GetYaxis().SetTitle(titley)
            c.SetMaximum(y_hi)
            c.SetMinimum(y_lo)
            c.Draw("HIST")
        else:
            c.Draw("HIST SAME")
    ROOT.SetOwnership(legend, False)
    legend.Draw()
    myText(0.20, 0.9, ROOT.kBlack, message)

    if save_pdf:
        cm1.Print(plot_title+".pdf")
    return cm1


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
