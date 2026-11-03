#!/usr/bin/env python3
"""x-y, x-z and z-r vertex maps of the whole traceback sample (all condor chunks together).

Draws the same figures as the per-chunk plots of hit_traceback.py (vertex_map_xy_*.png,
vertex_map_xz_*.png): the maps come from the same function, hit_traceback.plot_maps,
so the layout, colours, MDI outline, zoom panel and legend are identical. Only the input
differs: instead of the trees of one chunk, it reads all trees_f*_n*.jsonl files of the
unmerged output and plots them in one figure.

The maps need the case and the label of every vertex, which hit_traceback.py stores in
the trees_*.jsonl files (the skimmed EDM4hep files, merged or not, do not carry them),
so the trees are read, not the merged EDM4hep file.

Run in the Key4hep nightly (hit_traceback.py imports podio / uproot):
    python3 plot_merged_maps.py TREES_DIR OUTPUT_DIR
e.g.
    python3 plot_merged_maps.py .../traceback_skimmed/unmerged/trees .../traceback_skimmed/merged/plots

The output files are vertex_map_{xy,xz,zr}_f<first>_n<files>.png, with <first> and <files> the
first input file index and the number of input files of all chunks together, in the naming
of the chunk plots (the skipped, truncated input files are counted in <files>, as in the chunks).
"""

import argparse
import glob
import json
import os
import re
from collections import Counter

from hit_traceback import hit_summary, plot_maps

CHUNK = re.compile(r"trees_f(\d+)_n(\d+)\.jsonl$")


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("trees_dir", help="directory with the trees_f*_n*.jsonl files of the chunks")
    ap.add_argument("output_dir", help="the plots go here (keep it on EOS)")
    ap.add_argument("--collection", default="DCHCollection",
                    help="hit collection of the traceback, only for the plot title (default: DCHCollection)")
    ap.add_argument("--first-kink-only", action="store_true",
                    help="per tree, plot only the first kink of any region and the first kink in the MDI (no kink in "
                         "the MDI if the first vertex in the MDI is an ancestor); all ancestors stay. "
                         "Adds _first_kink to the output names")
    ap.add_argument("--tag", default=None,
                    help="name tag of the plots instead of f<first>_n<files> (e.g. for trees rebuilt from a merged skim)")
    return ap.parse_args()


def main():
    args = parse_args()
    paths = sorted(glob.glob(os.path.join(args.trees_dir, "trees_f*_n*.jsonl")))
    if not paths:
        raise SystemExit(f"no trees_f*_n*.jsonl in {args.trees_dir}")

    # the chunks have to be contiguous: [first, first + n) of the sorted list of input files
    chunks = sorted((int(m.group(1)), int(m.group(2))) for m in map(CHUNK.search, paths))
    for (f0, n0), (f1, _) in zip(chunks, chunks[1:]):
        if f0 + n0 != f1:
            print(f"WARNING: chunk gap or overlap between file {f0 + n0} and {f1}")
    first = chunks[0][0]
    n_files = sum(n for _, n in chunks)
    tag = args.tag or f"f{first:05d}_n{n_files}"

    storage = []                                    # all hit-family trees, as in hit_traceback.main
    for p in paths:
        with open(p) as fi:
            for line in fi:
                t = json.loads(line)
                t["nodes"] = {int(i): n for i, n in t["nodes"].items()}   # json turns the int keys into str
                storage.append(t)
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"chunks {len(paths)}, input files {n_files}")
    print(f"trees {len(storage)}, {sum(len(t['nodes']) for t in storage)} vertices")
    print("  vertex cases:", dict(sum((Counter(t["cases"]) for t in storage), Counter())))
    print(f"  {args.collection} hits by the cases on their path:", dict(hit_summary(storage)))

    plot_maps(storage, args.collection, args.output_dir, tag, first_kink_only=args.first_kink_only)


if __name__ == "__main__":
    main()
