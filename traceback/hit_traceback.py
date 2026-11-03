#!/usr/bin/env python3
"""Trace subdetector hits back to their MDI origin, and skim MCParticles.

Simple version of the traceback algorithm (pseudo-code revision of 2026-10-21):
  1. for every hit of --collection (DCH) in the event, grow a family tree from
     the hit particle up to the oldest parent (the root, the primary photon).
     For every photon ancestor, also note its daughters born before the chain
     daughter (older siblings = earlier scatter points of that photon);
  2. give every particle (vertex) of the trees a region label (injection / mdi /
     dch / other) and a case (hit_ancestor_mdi, hit_sibling_mdi, hit_ancestor_other,
     hit_sibling_other, unres_dropped, dch, injection);
  3. plot the x-y, x-z and z-r maps of all mdi and other vertices;
  4. skim MCParticles to the particles in the trees (storage) and write EDM4hep
     files with only MCParticles, --collection and EventHeader (one per input
     file, or one merged file). All other collections are dropped.

A truncated or corrupt input file (size != fEND, or any read error) is skipped as a
whole, with nothing of it in the trees or the skim; skipped_<tag>.txt lists the files
and reasons. The exit code is 0 unless no file at all could be processed.

Input: ddsim EDM4hep outputs of option F
(--part.keepAllParticles True --part.minimalKineticEnergy 0, user handler on).
The traceback reads with uproot (only the needed branches), the skim with podio
(it has to rewrite every collection with a particle link, uproot cannot write
EDM4hep).

Run in the Key4hep nightly that wrote the files:
    python3 hit_traceback.py INPUT_DIR OUTPUT_DIR [--n-files N] [--merge]
"""

# ---------------------------------------------------------------------------
# Headers, includes
# ---------------------------------------------------------------------------
import argparse
import glob
import json
import math
import os
from collections import Counter

import edm4hep                 # import before podio reads a frame: in the nightly, a second
import numpy as np             # dictionary load returns 1 and edm4hep raises
import podio
import uproot
from podio.root_io import Reader, Writer

BIT_CREATED_IN_SIM = 30       # MCParticle simulatorStatus bit
BIT_PRODUCED_BY_SECONDARY = 30  # SimTrackerHit quality bit

# MDI region: union of cylinders (ALLEGRO_o1_v03_shape_SRmasktest dimensions)
HALF_CROSSING = 0.015         # each beam axis is tilted by +-15 mrad in x-z [rad]
IP_RMAX = 35.0                # central + cone beam pipe, ConeBeamPipe_Rmax = 32.85 mm + wall [mm]
IP_ZMAX = 1190.0              # SeparatedBeamPipe_z [mm]
LEG_RMAX = 25.0               # separated pipes, SeparatedBeamPipe_rmax = 15 mm + wall + margin [mm]
LEG_ZMAX = 2350.0             # end of the legs (and of the cryostat cylinder) [mm]
CRYO_RMAX = 225.0             # cryostat envelope, cryo_out_r [mm]
CRYO_ZMIN = 1400.0            # cryostat nose start, cryo_nose_min_z [mm]
# DCH region: envelope of ALLEGRO_o1_v03 DectDimensions.xml (DCH_*_total)
DCH_RMIN, DCH_RMAX, DCH_ZMAX = 350.0, 2000.0, 3000.0

# vertex labels that go on the maps (dch and injection are stored, not plotted), colour by case
PLOT_LABELS = ("mdi", "other")
CASE_COLORS = {"hit_ancestor_mdi": "#000000", "hit_sibling_mdi": "#eb6834",
               "hit_ancestor_other": "#0a2fff", "hit_sibling_other": "#1baf7a"}
CASE_LEGEND = {"hit_ancestor_mdi": "Ancestor in MDI", "hit_sibling_mdi": "Ancestor kink in MDI",
               "hit_ancestor_other": "Ancestor outside MDI and DCH", "hit_sibling_other": "Ancestor kink outside MDI and DCH"}
HIT_CASES = tuple(CASE_COLORS)  # vertex cases on the path of a hit: chain (ancestor) or older sibling (kink), in the MDI or elsewhere
# marker area (default 8) and drawing order (default 1, higher = on top) of single cases: the few
# ancestor vertices are larger than the kinks and drawn last, so that they are not hidden by the kinks
CASE_SIZE = {"hit_ancestor_mdi": 20, "hit_ancestor_other": 20}
CASE_ZORDER = {"hit_ancestor_mdi": 3, "hit_ancestor_other": 3}
# grid spacing [mm] of the maps, (horizontal, vertical) axis; the zoomed-out x-z map keeps the default grid
GRID_STEPS = {("xy", False): (100, 100), ("xy", True): (25, 25), ("xz", True): (100, 25), ("zr", True): (100, 25)}

MC_FIELDS = {
    "pdg": "PDG", "gen": "generatorStatus", "sim": "simulatorStatus", "t": "time", "m": "mass",
    "x": "vertex.x", "y": "vertex.y", "z": "vertex.z",
    "px": "momentum.x", "py": "momentum.y", "pz": "momentum.z",
    "par_b": "parents_begin", "par_e": "parents_end", "dau_b": "daughters_begin", "dau_e": "daughters_end",
}


# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
def parse_args():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_dir", help="directory with the ddsim EDM4hep outputs (*.root)")
    ap.add_argument("output_dir", help="trees, plots and skimmed files go here (keep it on EOS)")
    ap.add_argument("--collection", default="DCHCollection",
                    help="hit collection to trace and keep in the skim (DCH settings hardcoded)")
    ap.add_argument("--n-files", type=int, default=None, help="number of files to process (default: all)")
    ap.add_argument("--first-file", type=int, default=0, help="index of the first file in the sorted list (condor chunks)")
    ap.add_argument("--n-events", type=int, default=None, help="events per file (default: all)")
    ap.add_argument("--merge", action="store_true", help="write one merged skimmed file instead of one per input")
    ap.add_argument("--no-skim", action="store_true", help="traceback and plots only")
    return ap.parse_args()


# ---------------------------------------------------------------------------
# Event access (uproot): one dict of numpy arrays per event
# ---------------------------------------------------------------------------
def check_complete(path):
    """raise if the ROOT file is truncated: size on disk != fEND in the header (EOS mv truncations)"""
    size = os.path.getsize(path)
    with open(path, "rb") as f:
        head = f.read(20)
    if len(head) < 20 or head[:4] != b"root":
        raise OSError(f"not a ROOT file ({size} bytes)")
    version = int.from_bytes(head[4:8], "big")
    fend = int.from_bytes(head[12:20] if version >= 1000000 else head[12:16], "big")
    if fend != size:
        raise OSError(f"truncated: {size} bytes on disk, fEND {fend}")


def read_events(path, n_events, coll):
    tree = uproot.open(path)["events"]
    branches = {key: f"MCParticles/MCParticles.{name}" for key, name in MC_FIELDS.items()}
    branches |= {
        "par_idx": "_MCParticles_parents/_MCParticles_parents.index",
        "dau_idx": "_MCParticles_daughters/_MCParticles_daughters.index",
    }
    # hits of coll (SimTrackerHit): particle index, energy, quality
    hit_branches = {
        "p": f"_{coll}_particle/_{coll}_particle.index",
        "e": f"{coll}/{coll}.eDep",
        "q": f"{coll}/{coll}.quality",
    }
    arrays = tree.arrays(list(branches.values()) + list(hit_branches.values()), entry_stop=n_events, library="np")
    for e in range(len(arrays[branches["pdg"]])):
        ev = {key: arrays[b][e] for key, b in branches.items()}
        ev["hits"] = {key: arrays[b][e] for key, b in hit_branches.items()}
        yield ev


def parents(ev, i):
    return ev["par_idx"][ev["par_b"][i]:ev["par_e"][i]]

def daughters(ev, i):
    return ev["dau_idx"][ev["dau_b"][i]:ev["dau_e"][i]]

def created_in_sim(ev, i):
    return bool((ev["sim"][i] >> BIT_CREATED_IN_SIM) & 1)

def ekin_kev(ev, i):
    m = ev["m"][i]
    e = math.sqrt(ev["px"][i] ** 2 + ev["py"][i] ** 2 + ev["pz"][i] ** 2 + m ** 2)
    return 1e6 * (e - m)


# ---------------------------------------------------------------------------
# Helper function: define detector regions
# ---------------------------------------------------------------------------
def in_mdi(x, y, z):
    """union of the IP cylinder, the two beam-pipe legs at +-15 mrad and the cryostat"""
    r = math.hypot(x, y)
    if r <= IP_RMAX and abs(z) <= IP_ZMAX:
        return True
    if abs(z) > LEG_ZMAX:
        return False
    if r <= CRYO_RMAX and abs(z) >= CRYO_ZMIN:
        return True
    for sign in (1, -1):
        ux, uz = sign * math.sin(HALF_CROSSING), math.cos(HALF_CROSSING)   # leg axis through the IP
        s = x * ux + z * uz                                                 # position along the axis
        if x * x + y * y + z * z - s * s <= LEG_RMAX ** 2:                  # distance to the axis
            return True
    return False

def in_dch(x, y, z):
    return DCH_RMIN <= math.hypot(x, y) <= DCH_RMAX and abs(z) <= DCH_ZMAX

def region_label(ev, i):
    if not created_in_sim(ev, i):
        return "injection"          # primary: its vertex is the SR source, not an interaction
    x, y, z = ev["x"][i], ev["y"][i], ev["z"][i]
    if in_mdi(x, y, z):
        return "mdi"
    if in_dch(x, y, z):
        return "dch"
    return "other"


# ---------------------------------------------------------------------------
# Helper function: get older siblings
# ---------------------------------------------------------------------------
def older_siblings(ev, i, q):
    """daughters of q born before particle i, sorted by time"""
    d = daughters(ev, q)
    d = d[ev["t"][d] < ev["t"][i]]
    return d[np.argsort(ev["t"][d], kind="stable")]


# ---------------------------------------------------------------------------
# Helper function: create trees
# ---------------------------------------------------------------------------
def make_node(ev, i, role):
    return dict(
        idx=int(i), pdg=int(ev["pdg"][i]), role=role, label=region_label(ev, i), case=None,
        x=float(ev["x"][i]), y=float(ev["y"][i]), z=float(ev["z"][i]), t=float(ev["t"][i]),
        ekin_kev=ekin_kev(ev, i),
        component=None,       # TODO: geometry component (MC vertices carry no volume name)
        parent=None, children=[], siblings_before=[], walked=False, root_status=None,
        n_hits=0, edep_kev=0.0, n_dropped=0,   # own hits (hit particles only)
    )

def grow_tree(ev, i, edep_kev, dropped, nodes):
    """Walk from hit particle i up to the root (the while-loop).

    nodes: {MCParticle index: node} of this event, shared by all hits. When the
    walk reaches a node that was already walked (an existing leg), it stops:
    the rest of the leg up to the root is already there.
    """
    if i not in nodes:
        nodes[i] = make_node(ev, i, "chain")
    leaf = nodes[i]
    leaf["n_hits"] += 1
    leaf["edep_kev"] += edep_kev
    leaf["n_dropped"] += dropped

    while True:
        node = nodes[i]
        if node["walked"]:                          # existing leg
            break

        node["walked"], node["role"] = True, "chain"
        par = parents(ev, i)
        if len(par) == 0:                           # no parent: this is the root
            complete = not created_in_sim(ev, i) and ev["gen"][i] == 1
            node["root_status"] = "complete" if complete else "broken"
            break
        q = int(par[0])                             # ddsim stores one parent
        node["parent"] = q
        if ev["pdg"][q] == 22:                      # photon parent: older daughters are earlier scatters
            for s in older_siblings(ev, i, q):
                s = int(s)
                if s not in nodes:
                    nodes[s] = make_node(ev, s, "sibling")
                    nodes[s]["parent"] = q
                node["siblings_before"].append(s)
        if q not in nodes:
            nodes[q] = make_node(ev, q, "chain")
        if i not in nodes[q]["children"]:
            nodes[q]["children"].append(i)
        i = q


def vertex_case(node):
    """classification of one vertex (particle) of a tree"""
    if node["root_status"] == "broken":
        return "unres_dropped"      # born in the simulation, but its parent was not stored
    if node["label"] in ("mdi", "other"):   # other = outside MDI and DCH (e.g. backscatter)
        return f"hit_{'ancestor' if node['role'] == 'chain' else 'sibling'}_{node['label']}"
    return node["label"]                    # dch or injection


def split_trees(nodes, file_name, event):
    """group the nodes of one event by the root of their tree (one tree per root photon)"""
    def root_of(i):
        while nodes[i]["parent"] is not None:
            i = nodes[i]["parent"]
        return i

    trees = {}
    for i, n in nodes.items():
        n["case"] = vertex_case(n)
        r = root_of(i)
        t = trees.setdefault(r, dict(file=file_name, event=event, root=r, root_status=nodes[r]["root_status"],
                                     nodes={}, cases=Counter()))
        t["nodes"][i] = n
        t["cases"][n["case"]] += 1
    return list(trees.values())


def hit_summary(trees):
    """per hit: which vertex cases are on its path (chain + older siblings)"""
    out = Counter()
    for t in trees:
        nodes = t["nodes"]
        for leaf in nodes.values():
            n_hits = leaf["n_hits"]
            if n_hits == 0:
                continue
            chain = [leaf["idx"]]
            while nodes[chain[-1]]["parent"] is not None:
                chain.append(nodes[chain[-1]]["parent"])
            cases = {nodes[j]["case"] for j in chain} | {nodes[s]["case"] for c in chain for s in nodes[c]["siblings_before"]}
            if leaf["n_dropped"] or len(chain) == 1 or "unres_dropped" in cases:
                key = "unres_dropped"
            elif not cases & set(HIT_CASES):
                key = "unres_nomdi"     # path has only dch / injection vertices
            else:
                key = "+".join(sorted(cases & set(HIT_CASES)))
            out[key] += n_hits
    return out


# ---------------------------------------------------------------------------
# Skim the files (podio): rebuild MCParticles and the traced hit collection
# ---------------------------------------------------------------------------
def copy_parameters(old, new):
    for key in old.parameters:                      # e.g. nPrimaries
        for typ in old.get_param_info(key):
            new.put_parameter(key, old.get_parameter(key, as_type=typ), as_type=typ)


def skim_frame(frame, keep, coll):
    new = podio.Frame()
    copy_parameters(frame, new)
    old_mc = frame.get("MCParticles")
    remap = {old: n for n, old in enumerate(keep)}
    mc = edm4hep.MCParticleCollection()
    muts = []
    for old in keep:
        m = old_mc[old].clone(False)
        mc.push_back(m)
        muts.append(m)
    for old in keep:
        p = old_mc[old]
        for q in p.getParents():
            if q.getObjectID().index in remap:
                muts[remap[old]].addToParents(muts[remap[q.getObjectID().index]])
        for q in p.getDaughters():
            if q.getObjectID().index in remap:
                muts[remap[old]].addToDaughters(muts[remap[q.getObjectID().index]])

    # hits of coll: every linked particle is a tree node, so in keep
    out = edm4hep.SimTrackerHitCollection()
    for h in frame.get(coll):
        m = h.clone(False)
        p = h.getParticle()
        if p.isAvailable():
            m.setParticle(muts[remap[p.getObjectID().index]])
        out.push_back(m)
    new.put(out, coll)

    if "EventHeader" in frame.getAvailableCollections():
        header = edm4hep.EventHeaderCollection()
        for h in frame.get("EventHeader"):
            header.push_back(h.clone())
        new.put(header, "EventHeader")
    new.put(mc, "MCParticles")                      # all other collections are dropped
    return new


def skim_file(path, keeps, coll, write_metadata):
    """keeps: list of kept MCParticle indices, one per event (same order as uproot).
    Returns (reader, [(frame, category)]); nothing is written, so a file that fails
    half way leaves no partial events in the output. Keep the reader alive until written."""
    reader = Reader(path)
    events = list(reader.get("events"))
    if len(events) != len(keeps):
        raise RuntimeError(f"podio reads {len(events)} events, uproot {len(keeps)}")
    frames = [(skim_frame(frame, keep, coll), "events") for frame, keep in zip(events, keeps)]
    for cat in ("runs", "metadata"):
        if cat in reader.categories and (cat == "runs" or write_metadata):
            frames += [(frame, cat) for frame in reader.get(cat)]
    return reader, frames


def write_frames(writer, frames):
    for frame, cat in frames:
        writer.write_frame(frame, cat)


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


def plot_maps(storage, coll, output_dir, tag, first_kink_only=False):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator

    title_tag = tag
    if first_kink_only:
        storage = first_kinks_only(storage)
        tag += "_first_kink"

    pts = {case: [] for case in CASE_COLORS}
    for t in storage:
        for n in t["nodes"].values():
            if n["label"] in PLOT_LABELS and n["case"] in pts:   # each vertex once, equal weight
                pts[n["case"]].append((n["x"], n["y"], n["z"], math.hypot(n["x"], n["y"])))   # (unres_dropped is not plotted)

    for view, (a, b) in {"xy": (0, 1), "xz": (2, 0), "zr": (2, 3)}.items():     # columns of pts: x, y, z, r
        fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
        for ax, zoom in zip(axes, (False, True)):
            for case, color in CASE_COLORS.items():
                v = np.array(pts[case]).reshape(-1, 4)
                ax.scatter(v[:, a], v[:, b], s=CASE_SIZE.get(case, 8), color=color, edgecolors="none", alpha=0.7,
                           zorder=CASE_ZORDER.get(case, 1), label=f"{CASE_LEGEND[case]} ({len(v)})")
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
        fig.suptitle(f"vertices relevant to {coll} hits, {view}, {title_tag}" + (", first kinks only" if first_kink_only else ""))
        fig.tight_layout()
        out = os.path.join(output_dir, f"vertex_map_{view}_{tag}.png")
        fig.savefig(out, dpi=130)
        plt.close(fig)
        print("wrote", out)


# ---------------------------------------------------------------------------
# Process one file
# ---------------------------------------------------------------------------
def trace_file(f, args):
    """Trees and kept MCParticle indices of one file; returns (trees, keeps, n_mc, n_no_particle)"""
    name = os.path.basename(f)
    check_complete(f)
    storage, keeps = [], []
    n_mc = n_no_particle = 0
    for e, ev in enumerate(read_events(f, args.n_events, args.collection)):
        # create the hit family tree for each hit
        nodes = {}
        d = ev["hits"]
        for h in range(len(d["p"])):
            i = int(d["p"][h])
            if i < 0:
                n_no_particle += 1
                continue
            dropped = bool((d["q"][h] >> BIT_PRODUCED_BY_SECONDARY) & 1)
            grow_tree(ev, i, 1e6 * float(d["e"][h]), dropped, nodes)
        trees = split_trees(nodes, name, e)
        storage += trees

        # every particle in the trees of this event -> retained MCParticles
        keeps.append(sorted(i for t in trees for i in t["nodes"]))
        n_mc += len(ev["pdg"])
    return storage, keeps, n_mc, n_no_particle


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    files = sorted(glob.glob(os.path.join(args.input_dir, "*.root")))
    last = None if args.n_files is None else args.first_file + args.n_files
    files = files[args.first_file:last]
    tag = f"f{args.first_file:05d}_n{len(files)}"
    os.makedirs(args.output_dir, exist_ok=True)

    writer = None
    if not args.no_skim:
        skim_dir = os.path.join(args.output_dir, "skim")
        os.makedirs(skim_dir, exist_ok=True)
        if args.merge:
            writer = Writer(os.path.join(skim_dir, f"skim_merged_{tag}.root"))

    storage = []                                    # all hit-family trees of the input_dir
    n_mc = n_keep = 0
    n_no_particle = 0
    n_done = 0
    skipped = []                                    # (file, reason): truncated or corrupt inputs
    for k, f in enumerate(files):
        name = os.path.basename(f)
        # a file that fails anywhere is skipped as a whole: nothing of it enters storage or the skim
        out_file = None
        try:
            trees, keeps, mc, no_part = trace_file(f, args)
            if not args.no_skim:
                reader, frames = skim_file(f, keeps, args.collection, write_metadata=(n_done == 0))
                if args.merge:
                    write_frames(writer, frames)
                else:
                    out_file = os.path.join(skim_dir, name.replace(".root", "_skim.root"))
                    w = Writer(out_file)
                    write_frames(w, frames)
                    w.finish()
                    out_file = None
        except Exception as exc:
            reason = f"{type(exc).__name__}: {exc}"
            print(f"[{k + 1}/{len(files)}] {name}: SKIPPED, {reason}", flush=True)
            skipped.append((name, reason))
            if out_file and os.path.exists(out_file):
                os.remove(out_file)                 # partial per-file skim
            continue
        storage += trees
        n_mc += mc
        n_keep += sum(len(keep) for keep in keeps)
        n_no_particle += no_part
        n_done += 1
        print(f"[{k + 1}/{len(files)}] {name}: {len(storage)} trees so far", flush=True)
    if writer is not None:
        writer.finish()

    if skipped:
        out = os.path.join(args.output_dir, f"skipped_{tag}.txt")
        with open(out, "w") as fo:
            for name, reason in skipped:
                fo.write(f"{name}\t{reason}\n")
        print(f"\nSKIPPED {len(skipped)} of {len(files)} files, listed in {out}")
    if files and not n_done:
        raise SystemExit("no input file could be processed")

    # output: one tree per line
    out = os.path.join(args.output_dir, f"trees_{tag}.jsonl")
    with open(out, "w") as fo:
        for t in storage:
            fo.write(json.dumps(t, default=sorted) + "\n")
    print("wrote", out)

    # summary
    print(f"\ntrees {len(storage)}, {sum(len(t['nodes']) for t in storage)} vertices")
    print(f"  vertex cases:", dict(sum((t["cases"] for t in storage), Counter())))
    print(f"  {args.collection} hits by the cases on their path:", dict(hit_summary(storage)))
    if n_no_particle:
        print("  hits without particle link:", n_no_particle)
    print(f"files processed {n_done} of {len(files)}")
    print(f"MCParticles kept {n_keep} of {n_mc} ({100 * n_keep / max(n_mc, 1):.2f}%)")

    plot_maps(storage, args.collection, args.output_dir, tag)


if __name__ == "__main__":
    main()
