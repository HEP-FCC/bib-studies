#!/usr/bin/env python3
"""x-y, x-z and z-r vertex maps from the traceback trees (trees_*.jsonl of hit_traceback.py).

All plotting of the traceback lives here; hit_traceback.py only traces, labels and skims, and
writes the trees (one tree per line, with the region, label and case of every vertex). So the
maps can be made from any selection of trees files: one condor chunk, a few chunks, all chunks
of a run, or one merged trees file.

Input: one or more trees files and/or directories (a directory = all trees_*.jsonl in it).
Trees files are JSON lines, so merging them is concatenation: --merged-trees FILE writes all
the trees that were read (every collection) into one file, which can be plotted again later.
(Alternative for a merged EDM4hep skim: hit_traceback.py <skim dir> <out> --no-skim rebuilds
one trees file from the skim; validated node-by-node against the raw run on 2026-10-24.)

One set of three maps per plotted collection (--collections, default: every VTX / DCH / SiWr
collection with trees in the input; calorimeter / muon trees of older all-collection runs are
not plotted), each with a full panel (autoscaled to the plotted vertices) and a zoom
on the MDI region. Plotted: the mdi and other vertices (cases hit_ancestor_mdi, hit_sibling_mdi,
hit_ancestor_other, hit_sibling_other), each once with equal weight; the own-subdetector,
injection and unres_dropped vertices are not plotted (an empty map says what it did not plot).
--primary-kinks-only: only the kinks of SR primaries before their own real hits in the
collection (vertices with primary_kink, trees written since 2026-11-06), e.g. the Compton
recoil in the MDI of a primary photon that then deposits energy in the DCH itself.

Output: vertex_map_{xy,xz,zr}_<collection>_<tag>[_primary_kinks][_first_kink].png. The tag is
f<first>_n<files> when every input is a chunk file trees_f<first>_n<files>.jsonl (contiguous
chunks expected, gaps are reported; the skipped, truncated input files are counted in
<files>, as in the chunks), the name of the file for a single other trees file, otherwise
--tag is needed.

The labels are taken from the trees, not recomputed. Trees written before 2026-11-06 carry the
DCH labels for every collection: a warning is printed when the labels of a non-DCH collection
do not match its own region (regenerate those trees with hit_traceback.py). Old trees without
a "coll" field are DCH trees.

Run in the Key4hep nightly (hit_traceback.py, imported for the MDI dimensions and the labels,
needs podio / uproot):
    python3 plot_merged_maps.py TREES [TREES ...] -o OUTPUT_DIR [--collections all|a,b]
                                [--primary-kinks-only] [--first-kink-only] [--tag T] [--merged-trees FILE]
e.g.
    python3 plot_merged_maps.py .../traceback_allcoll/trees -o .../traceback_allcoll/plots
    python3 plot_merged_maps.py .../trees/trees_f01000_n200.jsonl -o /tmp/plots --collections SiWrBCollection
    python3 plot_merged_maps.py .../trees -o /tmp/plots --collections DCHCollection --primary-kinks-only
"""

import argparse
import glob
import json
import math
import os
import re
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from hit_traceback import (CRYO_RMAX, CRYO_ZMIN, HALF_CROSSING, IP_RMAX, IP_ZMAX, LEG_RMAX, LEG_ZMAX,
                           coarse_label, hit_summary, own_region)

CHUNK = re.compile(r"trees_f(\d+)_n(\d+)\.jsonl$")

REGION_NAMES = {"vtx": "VTX", "dch": "DCH", "siwr": "SiWr"}   # legend text

# vertex labels that go on the maps (own region and injection are stored, not plotted), colour by case
PLOT_LABELS = ("mdi", "other")
CASE_COLORS = {"hit_ancestor_mdi": "#000000", "hit_sibling_mdi": "#eb6834",
               "hit_ancestor_other": "#0a2fff", "hit_sibling_other": "#1baf7a"}
# {det} = the own subdetector of the plotted collection (REGION_NAMES)
CASE_LEGEND = {"hit_ancestor_mdi": "Ancestor in MDI", "hit_sibling_mdi": "Ancestor kink in MDI",
               "hit_ancestor_other": "Ancestor outside MDI and {det}", "hit_sibling_other": "Ancestor kink outside MDI and {det}"}
# marker area (default 8) and drawing order (default 1, higher = on top) of single cases: the few
# ancestor vertices are larger than the kinks and drawn last, so that they are not hidden by the kinks
CASE_SIZE = {"hit_ancestor_mdi": 20, "hit_ancestor_other": 20}
CASE_ZORDER = {"hit_ancestor_mdi": 3, "hit_ancestor_other": 3}
# grid spacing [mm] of the maps, (horizontal, vertical) axis; the zoomed-out x-z map keeps the default grid
GRID_STEPS = {("xy", False): (100, 100), ("xy", True): (25, 25), ("xz", True): (100, 25), ("zr", True): (100, 25)}


# ---------------------------------------------------------------------------
# Plot the maps
# ---------------------------------------------------------------------------
def first_kinks_only(storage):
    """The trees as the maps need them with --first-kink-only: per tree (one root photon of one event)
    the kinks (older siblings) are cut to the first kink of any region and the first kink in the MDI
    (the same vertex if the first kink is in the MDI); the ancestors (chain vertices) stay as they are.
    "First" = smallest production time. If the first vertex in the MDI of a tree is an ancestor (a
    black one, born at or before the first kink in the MDI), the tree has no kink in the MDI to plot.
    Only a view for plotting: the trees of the storage are not changed.
    """
    out = []
    for t in storage:
        nodes = t["nodes"]
        kinks = [n for n in nodes.values() if n["role"] == "sibling"]
        keep = set()
        if kinks:
            keep.add(min(kinks, key=lambda n: n["t"])["idx"])           # first kink of any region
        mdi = [n for n in nodes.values() if n["case"] in ("hit_ancestor_mdi", "hit_sibling_mdi")]
        if mdi:
            first_mdi = min(mdi, key=lambda n: (n["t"], n["case"] != "hit_ancestor_mdi"))   # a tie goes to the ancestor
            if first_mdi["case"] == "hit_sibling_mdi":
                keep.add(first_mdi["idx"])                              # first kink in the MDI
            else:                                                       # first vertex in the MDI is an ancestor:
                keep = {i for i in keep if nodes[i]["case"] != "hit_sibling_mdi"}   # no orange kink for this tree
        out.append(dict(t, nodes={i: n for i, n in nodes.items() if n["role"] != "sibling" or i in keep}))
    return out


def primary_kinks_only(storage):
    """Only the trees with real hits of their SR primary, cut to the kinks of the primary before its hits
    (primary_kink). A view for plotting, the trees of the storage are not changed."""
    return [dict(t, nodes={i: n for i, n in t["nodes"].items() if n.get("primary_kink")})
            for t in storage if t.get("n_primary_hits")]


def plot_maps(storage, coll, output_dir, tag, first_kink_only=False, primary_only=False):
    title_tag = tag
    if primary_only:                                # before the first-kink cut: first kinks of the primary
        storage = primary_kinks_only(storage)
        tag += "_primary_kinks"
    if first_kink_only:
        storage = first_kinks_only(storage)
        tag += "_first_kink"

    own = own_region(coll)
    det = REGION_NAMES.get(own, coll)               # legend: "outside MDI and <det>"
    pts = {case: [] for case in CASE_COLORS}
    for t in storage:
        for n in t["nodes"].values():
            if n["label"] in PLOT_LABELS and n["case"] in pts:   # each vertex once, equal weight
                pts[n["case"]].append((n["x"], n["y"], n["z"], math.hypot(n["x"], n["y"])))   # (unres_dropped is not plotted)
    n_plotted = sum(len(v) for v in pts.values())
    n_not_plotted = Counter(n["case"] for t in storage for n in t["nodes"].values()
                            if not (n["label"] in PLOT_LABELS and n["case"] in pts))

    for view, (a, b) in {"xy": (0, 1), "xz": (2, 0), "zr": (2, 3)}.items():     # columns of pts: x, y, z, r
        fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
        for ax, zoom in zip(axes, (False, True)):
            for case, color in CASE_COLORS.items():
                v = np.array(pts[case]).reshape(-1, 4)
                ax.scatter(v[:, a], v[:, b], s=CASE_SIZE.get(case, 8), color=color, edgecolors="none", alpha=0.7,
                           zorder=CASE_ZORDER.get(case, 1), label=f"{CASE_LEGEND[case].format(det=det)} ({len(v)})")
            if not n_plotted:                       # say why the map is empty instead of leaving a blank frame
                why = ", ".join(f"{c} {k}" for c, k in sorted(n_not_plotted.items())) or "none"
                ax.text(0.5, 0.5, f"no plotted vertices\n{len(storage)} trees, vertices not plotted: {why}",
                        transform=ax.transAxes, ha="center", va="center", fontsize=9, color="#555555")
            if view == "xz":                        # MDI region outline
                kw = dict(color="#888888", lw=0.8)
                ax.plot([-IP_ZMAX, IP_ZMAX, IP_ZMAX, -IP_ZMAX, -IP_ZMAX], [IP_RMAX, IP_RMAX, -IP_RMAX, -IP_RMAX, IP_RMAX], **kw)
                zz = np.array([-LEG_ZMAX, LEG_ZMAX])
                for sign in (1, -1):
                    for off in (LEG_RMAX, -LEG_RMAX):
                        ax.plot(zz, sign * HALF_CROSSING * zz + off, **kw)
                    zc = sign * np.array([CRYO_ZMIN, CRYO_ZMIN, LEG_ZMAX, LEG_ZMAX, CRYO_ZMIN])
                    ax.plot(zc, [-CRYO_RMAX, CRYO_RMAX, CRYO_RMAX, -CRYO_RMAX, -CRYO_RMAX], **kw)
            if view == "zr":                        # MDI region outline in (z, r): the two legs coincide
                kw = dict(color="#888888", lw=0.8)
                ax.plot([-IP_ZMAX, -IP_ZMAX, IP_ZMAX, IP_ZMAX], [0, IP_RMAX, IP_RMAX, 0], **kw)
                zz = np.linspace(-LEG_ZMAX, LEG_ZMAX, 201)
                ax.plot(zz, HALF_CROSSING * abs(zz) + LEG_RMAX, **kw)
                ax.plot(zz, np.maximum(0, HALF_CROSSING * abs(zz) - LEG_RMAX), **kw)
                for sign in (1, -1):
                    zc = sign * np.array([CRYO_ZMIN, CRYO_ZMIN, LEG_ZMAX, LEG_ZMAX])
                    ax.plot(zc, [0, CRYO_RMAX, CRYO_RMAX, 0], **kw)
            if zoom:
                if view == "xy":
                    ax.set_xlim(-300, 300), ax.set_ylim(-300, 300)
                else:
                    ax.set_xlim(-LEG_ZMAX - 100, LEG_ZMAX + 100), ax.set_ylim(-250 if view == "xz" else 0, 250)
                ax.set_title("zoom on the MDI region")
            else:
                ax.set_title("all plotted vertices")
                ax.legend(loc="upper right", fontsize=9)
            ax.set_xlabel("xyzr"[a] + " [mm]")
            ax.set_ylabel("xyzr"[b] + " [mm]")
            ax.grid(color="#e0e0e0", lw=0.5)
            if (view, zoom) in GRID_STEPS:          # finer grid: a gridline every step, labels stay at the default ticks
                for axis, step in zip((ax.xaxis, ax.yaxis), GRID_STEPS[view, zoom]):
                    axis.set_minor_locator(MultipleLocator(step))
                ax.grid(which="both", color="#e0e0e0", lw=0.5)
                ax.tick_params(which="minor", length=0)
        fig.suptitle(f"vertices relevant to {coll} hits, {view}, {title_tag}"
                     + (", kinks of SR primaries before their own hits" if primary_only else "")
                     + (", first kinks only" if first_kink_only else ""))
        fig.tight_layout()
        out = os.path.join(output_dir, f"vertex_map_{view}_{coll}_{tag}.png")
        fig.savefig(out, dpi=130)
        plt.close(fig)
        print("wrote", out)


# ---------------------------------------------------------------------------
# Read the trees, main function
# ---------------------------------------------------------------------------
def parse_args():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("trees", nargs="+",
                    help="trees files (*.jsonl) and/or directories (all trees_*.jsonl in them)")
    ap.add_argument("-o", "--output-dir", required=True, help="the plots go here (keep it on EOS)")
    ap.add_argument("--collections", default="all",
                    help="collections to plot: 'all' (every VTX / DCH / SiWr collection with trees in the input, "
                         "default) or a comma-separated list, e.g. DCHCollection,SiWrBCollection")
    ap.add_argument("--primary-kinks-only", action="store_true",
                    help="plot only the kinks of SR primaries before their own real hits in the collection "
                         "(primary_kink vertices). Adds _primary_kinks to the output names")
    ap.add_argument("--first-kink-only", action="store_true",
                    help="per tree, plot only the first kink of any region and the first kink in the MDI (no kink in "
                         "the MDI if the first vertex in the MDI is an ancestor); all ancestors stay. "
                         "Adds _first_kink to the output names")
    ap.add_argument("--tag", default=None,
                    help="name tag of the plots instead of the automatic one (see the description)")
    ap.add_argument("--merged-trees", default=None,
                    help="also write all trees read (every collection) into this one .jsonl file")
    return ap.parse_args()


def find_trees_files(paths):
    files = []
    for p in paths:
        if os.path.isdir(p):
            found = sorted(glob.glob(os.path.join(p, "trees_*.jsonl")))
            if not found:
                raise SystemExit(f"no trees_*.jsonl in {p}")
            files += found
        elif os.path.isfile(p):
            files.append(p)
        else:
            raise SystemExit(f"no such file or directory: {p}")
    return sorted(set(files), key=files.index)


def make_tag(files):
    chunks = [CHUNK.search(os.path.basename(f)) for f in files]
    if all(chunks):
        chunks = sorted((int(m.group(1)), int(m.group(2))) for m in chunks)
        for (f0, n0), (f1, _) in zip(chunks, chunks[1:]):
            if f0 + n0 != f1:
                print(f"WARNING: chunk gap or overlap between file {f0 + n0} and {f1}")
        return f"f{chunks[0][0]:05d}_n{sum(n for _, n in chunks)}"
    if len(files) == 1:
        return re.sub(r"^trees_", "", os.path.basename(files[0]).removesuffix(".jsonl"))
    raise SystemExit("the inputs are not all trees_f<first>_n<files>.jsonl chunk files: give --tag")


def main():
    args = parse_args()
    files = find_trees_files(args.trees)
    merged = None
    if args.merged_trees:
        # a merged file from an earlier run in an input directory would be read twice: leave it out
        files = [f for f in files if os.path.abspath(f) != os.path.abspath(args.merged_trees)]
        if not files:
            raise SystemExit("no trees files left to read")
    tag = args.tag or make_tag(files)
    os.makedirs(args.output_dir, exist_ok=True)
    if args.merged_trees:
        os.makedirs(os.path.dirname(os.path.abspath(args.merged_trees)), exist_ok=True)
        merged = open(args.merged_trees, "w")

    by_coll = {}                                    # {collection: [trees]}
    for p in files:
        with open(p) as fi:
            for line in fi:
                if merged:
                    merged.write(line if line.endswith("\n") else line + "\n")
                t = json.loads(line)
                t["nodes"] = {int(i): n for i, n in t["nodes"].items()}   # json turns the int keys into str
                by_coll.setdefault(t.get("coll", "DCHCollection"), []).append(t)   # old trees files have no coll: DCH
    if merged:
        merged.close()
        print("wrote", args.merged_trees)

    if args.collections == "all":
        colls = sorted(c for c in by_coll if own_region(c) is not None)
        other = sorted(c for c in by_coll if own_region(c) is None)
        if other:
            print(f"not plotted (not VTX / DCH / SiWr, cannot be traced back with the sample): {other}")
    else:
        colls = [c for c in args.collections.split(",") if c]
        missing = [c for c in colls if c not in by_coll]
        if missing:
            print(f"WARNING: no trees of {missing} in the input (collections with trees: {sorted(by_coll)}); "
                  "their maps are empty")

    print(f"trees files {len(files)}, tag {tag}")
    for coll in colls:
        storage = by_coll.get(coll, [])
        own = own_region(coll)
        print(f"\n{coll} (own region {own}): {len(storage)} trees, {sum(len(t['nodes']) for t in storage)} vertices")
        print("  vertex cases:", dict(sum((Counter(t["cases"]) for t in storage), Counter())))
        print(f"  hits by the cases on their path:", dict(hit_summary(storage)))
        stale = sum(n["label"] != coarse_label(n["region"], own)
                    for t in storage for n in t["nodes"].values() if "region" in n)
        if stale:
            print(f"  WARNING: {stale} vertices have labels that do not match the own region {own} of {coll}: "
                  "trees written with the old DCH-only labels, regenerate them with hit_traceback.py")
        if storage and "n_primary_hits" not in storage[0]:
            print("  NOTE: trees written before 2026-11-06 have no primary-hit fields (primary hits count as unres_dropped)")
        plot_maps(storage, coll, args.output_dir, tag, first_kink_only=args.first_kink_only,
                  primary_only=args.primary_kinks_only)


if __name__ == "__main__":
    main()
