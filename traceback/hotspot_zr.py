#!/usr/bin/env python3
"""TH2D of the z-r position of the vertices of one case, in the window of the "zoom on the MDI region" maps.

Default case: hit_sibling_mdi (the red/orange points of the maps, "Ancestor kink in MDI"). The count per
bin shows whether the kinks are concentrated in a hotspot (e.g. at the crotch of the two beam pipes) or
spread along the beam pipe. The MDI outline of the z-r maps is drawn on top; the 10 fullest bins are printed.

Reads the trees_f*_n*.jsonl files of the unmerged traceback output, like plot_merged_maps.py.
Run in the Key4hep nightly (ROOT; hit_traceback.py, imported for the MDI dimensions, needs podio / uproot):
    python3 hotspot_zr.py TREES_DIR OUTPUT_DIR [--case hit_sibling_mdi] [--dz 10] [--dr 5] [--logz]
        [--zmin -1450 --zmax -1150 --rmax 50]   # zoom on the hotspot
Writes OUTPUT_DIR/hotspot_zr_<case>[_zoom][_logz].{png,root} (the TH2D is in the root file).
The window defaults to the one of the zoomed z-r map (|z| <= 2450 mm, r <= 250 mm).
"""

import argparse
import glob
import json
import math
import os

import numpy as np
import ROOT

from hit_traceback import HALF_CROSSING, CRYO_RMAX, CRYO_ZMIN, IP_RMAX, IP_ZMAX, LEG_RMAX, LEG_ZMAX

ZMAX_ZOOM, RMAX_ZOOM = LEG_ZMAX + 100, 250.0       # default window: the zoomed z-r map [mm]


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("trees_dir", help="directory with the trees_f*_n*.jsonl files")
    ap.add_argument("output_dir", help="the png and the root file go here")
    ap.add_argument("--case", default="hit_sibling_mdi", help="vertex case to histogram (default: hit_sibling_mdi)")
    ap.add_argument("--dz", type=float, default=10.0, help="bin width in z [mm] (default: 10)")
    ap.add_argument("--dr", type=float, default=5.0, help="bin width in r [mm] (default: 5)")
    ap.add_argument("--zmin", type=float, default=-ZMAX_ZOOM, help="lower edge in z [mm]")
    ap.add_argument("--zmax", type=float, default=ZMAX_ZOOM, help="upper edge in z [mm]")
    ap.add_argument("--rmax", type=float, default=RMAX_ZOOM, help="upper edge in r [mm] (lower edge is 0)")
    ap.add_argument("--logz", action="store_true", help="logarithmic colour scale")
    return ap.parse_args()


def outline():
    """MDI outline in (z, r), the same lines as in hit_traceback.plot_maps: list of TGraph"""
    def graph(z, r):
        g = ROOT.TGraph(len(z), np.asarray(z, dtype="float64"), np.asarray(r, dtype="float64"))
        g.SetLineColor(ROOT.kGray + 2)
        g.SetLineWidth(1)
        return g
    zz = np.linspace(-LEG_ZMAX, LEG_ZMAX, 201)
    graphs = [graph([-IP_ZMAX, -IP_ZMAX, IP_ZMAX, IP_ZMAX], [0, IP_RMAX, IP_RMAX, 0]),
              graph(zz, HALF_CROSSING * abs(zz) + LEG_RMAX),
              graph(zz, np.maximum(0, HALF_CROSSING * abs(zz) - LEG_RMAX))]
    for sign in (1, -1):
        graphs.append(graph(sign * np.array([CRYO_ZMIN, CRYO_ZMIN, LEG_ZMAX, LEG_ZMAX]), [0, CRYO_RMAX, CRYO_RMAX, 0]))
    return graphs


def main():
    args = parse_args()
    paths = sorted(glob.glob(os.path.join(args.trees_dir, "trees_f*_n*.jsonl")))
    if not paths:
        raise SystemExit(f"no trees_f*_n*.jsonl in {args.trees_dir}")
    os.makedirs(args.output_dir, exist_ok=True)

    nz, nr = round((args.zmax - args.zmin) / args.dz), round(args.rmax / args.dr)
    h = ROOT.TH2D(f"h_zr_{args.case}", f"{args.case} vertices, MDI zoom;z [mm];r [mm];vertices per bin",
                  nz, args.zmin, args.zmax, nr, 0, args.rmax)
    h.SetStats(False)
    n_case = 0
    for p in paths:
        with open(p) as fi:
            for line in fi:
                for n in json.loads(line)["nodes"].values():
                    if n["label"] in ("mdi", "other") and n["case"] == args.case:
                        n_case += 1
                        h.Fill(n["z"], math.hypot(n["x"], n["y"]))
    n_win = int(h.Integral())
    print(f"{args.case}: {n_case} vertices, {n_win} in the window z [{args.zmin:.0f}, {args.zmax:.0f}] r [0, {args.rmax:.0f}] "
          f"({nz} x {nr} bins of {args.dz} x {args.dr} mm)")

    top = sorted(((h.GetBinContent(i, j), i, j) for i in range(1, nz + 1) for j in range(1, nr + 1)), reverse=True)[:10]
    print("fullest bins:")
    for c, i, j in top:
        print(f"  {int(c):6d}  z [{h.GetXaxis().GetBinLowEdge(i):7.0f}, {h.GetXaxis().GetBinUpEdge(i):7.0f}]"
              f"  r [{h.GetYaxis().GetBinLowEdge(j):5.0f}, {h.GetYaxis().GetBinUpEdge(j):5.0f}]")

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetPalette(ROOT.kBird)
    c = ROOT.TCanvas("c", "", 1400, 560)
    c.SetRightMargin(0.14)
    c.SetLogz(args.logz)
    h.Draw("COLZ")
    graphs = outline()                              # keep the references alive until the canvas is saved
    for g in graphs:
        g.Draw("L SAME")
    zoomed = (args.zmin, args.zmax, args.rmax) != (-ZMAX_ZOOM, ZMAX_ZOOM, RMAX_ZOOM)
    tag = args.case + ("_zoom" if zoomed else "") + ("_logz" if args.logz else "")
    c.SaveAs(os.path.join(args.output_dir, f"hotspot_zr_{tag}.png"))
    fo = ROOT.TFile(os.path.join(args.output_dir, f"hotspot_zr_{tag}.root"), "RECREATE")
    h.Write()
    fo.Close()
    print("wrote", os.path.join(args.output_dir, f"hotspot_zr_{tag}.{{png,root}}"))


if __name__ == "__main__":
    main()
