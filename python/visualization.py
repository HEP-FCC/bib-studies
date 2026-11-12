"""
    Plotting and visualization tools
"""
import math

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


def _axis_ticks(lo, hi, log):
    """Values of the major and minor ticks of an axis between lo and hi"""
    if log:
        major, minor = [], []
        for n in range(math.floor(math.log10(lo)), math.ceil(math.log10(hi)) + 1):
            if lo <= 10**n <= hi:
                major.append(10**n)
            minor += [m * 10**n for m in range(2, 10) if lo <= m * 10**n <= hi]
        return major, minor
    raw = (hi - lo) / 5
    magnitude = 10**math.floor(math.log10(raw))
    step = next(f * magnitude for f in (1, 2, 5, 10) if f * magnitude >= raw)
    k_lo, k_hi = math.ceil(lo / (step / 5) - 1e-9), math.floor(hi / (step / 5) + 1e-9)
    major = [k * step / 5 for k in range(k_lo, k_hi + 1) if k % 5 == 0]
    minor = [k * step / 5 for k in range(k_lo, k_hi + 1) if k % 5 != 0]
    return major, minor


def draw_axis_ndc(x, y_bottom, y_top, lo, hi, log, title, label_size=0.04, title_size=0.045):
    """
    Draw a vertical axis (ticks and labels on its right) at the NDC position x, between y_bottom and y_top, for the values lo to hi.
    It is made of plain lines and text in the pad (no TGaxis, no sub-pad), so that it does not depend on how a viewer (e.g. JSROOT) handles those.
    """
    def position(v):
        fraction = (math.log10(v) - math.log10(lo)) / (math.log10(hi) - math.log10(lo)) if log else (v - lo) / (hi - lo)
        return y_bottom + fraction * (y_top - y_bottom)

    def label(v):
        if not log:
            return f"{v:g}"
        n = round(math.log10(v))
        return "1" if n == 0 else "10" if n == 1 else f"10^{{{n}}}"

    line = ROOT.TLine()
    line.SetLineColor(ROOT.kBlack)
    line.DrawLineNDC(x, y_bottom, x, y_top)
    major, minor = _axis_ticks(lo, hi, log)
    for v in minor:
        line.DrawLineNDC(x, position(v), x + 0.006, position(v))
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(label_size)
    latex.SetTextAlign(12)  # left, vertically centered
    for v in major:
        line.DrawLineNDC(x, position(v), x + 0.012, position(v))
        latex.DrawLatex(x + 0.018, position(v), label(v))
    latex.SetTextSize(title_size)
    latex.SetTextAngle(90)
    latex.SetTextAlign(22)  # centered
    latex.DrawLatex(x + 0.105, (y_bottom + y_top) / 2, title)


def draw_overlay(hists, titlex, titley, plot_title, message='', log_x=False, log_y=False, legend_title='', save_pdf=True, right_axes=None, title=''):
    """
    Overlay several 1D histograms (dict {legend label: histogram}) on one canvas, e.g. one per layer.
    The histograms are cloned, so the originals are left untouched. Empty histograms are skipped.
    The canvas is printed to <plot_title>.pdf if save_pdf, and returned (None if there is nothing to draw),
    so that it can be saved in a ROOT file.

    right_axes: list of (scale, title) or (scale, title, note) to add axes on the right, which show the
    main y axis multiplied by 'scale' (e.g. to convert a density per event to a rate), the optional note
    (e.g. the assumptions used for the scale) is written under the message.

    title: title of the plot (and of the canvas), written on top, with the message (e.g. the collection) under it.
    """
    right_axes = right_axes or []
    colors = [ROOT.kBlack, ROOT.kRed+1, ROOT.kBlue+1, ROOT.kGreen+2, ROOT.kMagenta+1, ROOT.kOrange+7,
              ROOT.kCyan+1, ROOT.kViolet+1, ROOT.kPink+1, ROOT.kTeal+2, ROOT.kGray+2]
    hists = {label: h for label, h in hists.items() if h.GetEntries() > 0}
    if not hists:
        print(f"WARNING: no filled histogram to overlay for {plot_title}, plot skipped")
        return None

    cm1 = ROOT.TCanvas("c_"+plot_title, title or "c_"+plot_title, 800, 600)
    right_margin = 0.15 + 0.12 * len(right_axes)  # room for the axes on the right
    cm1.SetRightMargin(right_margin)
    if log_x:
        cm1.SetLogx()
    if log_y:
        cm1.SetLogy()
    legend = ROOT.TLegend(0.35, 0.67, 0.99 - right_margin, 0.83)  # inside the frame
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

    def text(x, y, string, size):  # TLatex with its own size (myText always draws with the default one)
        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextFont(42)
        latex.SetTextSize(size)
        latex.DrawLatex(x, y, string)

    if title:
        text(0.20, 0.95, title, 0.045)
        text(0.20, 0.905, message, 0.035)
        note_y = 0.865
    else:
        text(0.20, 0.9, message, 0.04)
        note_y = 0.86

    # axes on the right: the main y axis times a scale factor, with the same range and position as the frame of the main pad
    if right_axes and (y_lo > 0 or not log_y):
        y_bottom, y_top = cm1.GetBottomMargin(), 1 - cm1.GetTopMargin()
        for i, (scale, axis_title, *note) in enumerate(right_axes):
            draw_axis_ndc(1 - right_margin + i * 0.12, y_bottom, y_top, y_lo * scale, y_hi * scale, log_y, axis_title,
                          label_size=ROOT.gStyle.GetLabelSize("Y"), title_size=ROOT.gStyle.GetTitleSize("Y"))
            if note and note[0]:
                text(0.20, note_y - 0.04 * i, note[0], 0.03)

    if save_pdf:
        cm1.Print(plot_title+".pdf")
    return cm1


def draw_hist(hist, titlex, titley, plot_title, message='', log_x = False, log_y=False, draw_opt="HistE", top_axis=None):
    """top_axis = (scale, title) adds a second x axis on top whose values are scale * the bottom ones (e.g. c*t from t)."""

    cm1 = ROOT.TCanvas("cm1_"+message, "cm1_"+message, 800, 600)
    if top_axis:
        cm1.SetTopMargin(0.25)  # room for the labels and the title of the second axis, below the plot title

    hist.GetXaxis().SetTitle(titlex)
    hist.GetYaxis().SetTitle(titley)
    hist.Draw(draw_opt)
    entries = hist.GetEntries()

    myText(0.20, 0.9, ROOT.kBlack, message + '  , entries = ' + str(entries))

    if log_x:
        cm1.SetLogx()
    if log_y:
        cm1.SetLogy()

    if top_axis:
        scale, top_title = top_axis
        cm1.Update()
        x = hist.GetXaxis()
        y_top = 10 ** cm1.GetUymax() if log_y else cm1.GetUymax()
        axis = ROOT.TGaxis(x.GetXmin(), y_top, x.GetXmax(), y_top, scale * x.GetXmin(), scale * x.GetXmax(), 510, "-")
        axis.SetTitle(top_title)
        axis.SetLabelFont(x.GetLabelFont())
        axis.SetLabelSize(x.GetLabelSize())
        axis.SetTitleFont(x.GetTitleFont())
        axis.SetTitleSize(x.GetTitleSize())
        axis.Draw()

    cm1.Print(plot_title+".pdf")
    #cm1.Print(plot_title+".png")
