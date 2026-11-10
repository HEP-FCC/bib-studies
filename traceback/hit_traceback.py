#!/usr/bin/env python3
"""Trace subdetector hits back to their MDI origin, and skim MCParticles.

Simple version of the traceback algorithm (pseudo-code revision of 2026-10-21):
  1. for every hit of every tracker hit collection of the VTX, DCH and SiWr in the
     event (found automatically, or the list of --hit-collections), grow a family tree
     from the hit particle up to the oldest parent (the root, the primary photon).
     For every photon ancestor, also note its daughters born before the chain
     daughter (older siblings = earlier scatter points of that photon).
     The trees are built per collection (every tree has a "coll" field, e.g.
     DCHCollection, SiWrBCollection). Only these subdetectors are traced: with the
     user particle handler of the sample, particles born outside the tracker region
     are not stored, so calorimeter / muon hits cannot be traced back;
     Primary hits: if the hit particle is the SR primary photon itself and the hit is
     real (producedBySecondary not set: the photon step deposited the energy, e.g. a
     Compton or photoelectric electron that was not produced), its own daughters born
     before the hit time are its kinks (scatter points before the hit, e.g. a Compton
     recoil in the MDI). They are stored like older siblings (role sibling, cases
     hit_sibling_mdi / hit_sibling_other), with primary_kink = true; the primary node
     gets primary_hits (number of its real hits) and kinks_before_hit, the tree
     n_primary_hits;
  2. give every particle (vertex) of the trees a region (injection / mdi / vtx / dch /
     siwr / other), a coarse label and a case. The label is relative to the
     subdetector of the hit collection (its "own" region, from the collection name,
     COLLECTION_REGIONS): injection / mdi / <own region> / other, so for DCHCollection
     it is injection / mdi / dch / other, for SiWrBCollection injection / mdi / siwr /
     other (a vertex in the DCH is then "other"). Cases:
     hit_ancestor_mdi, hit_sibling_mdi, hit_ancestor_other, hit_sibling_other,
     unres_dropped, <own region>, injection;
  3. write the trees (trees_<tag>.jsonl, one tree per line, with the regions, labels
     and cases) and a per-collection summary (collections_<tag>.json);
  4. skim MCParticles to the particles in the trees (storage) and write EDM4hep
     files with only MCParticles, EventHeader and the traced tracker hit collections,
     one per input file or one merged file. All other collections are dropped.
No plots: the vertex maps are made from the trees files by plot_merged_maps.py
(any chunk, set of chunks or merged trees file).

A truncated or corrupt input file (size != fEND, or any read error) is skipped as a
whole, with nothing of it in the trees or the skim; skipped_<tag>.txt lists the files
and reasons. The exit code is 0 unless no file at all could be processed.

Input: ddsim EDM4hep outputs:
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
# Subdetector regions: envelopes of ALLEGRO_o1_v03 (DectDimensions.xml), boxes (rmin, rmax, |z|min, |z|max) [mm]
# symmetric in z. A vertex belongs to the first subdetector of this ordered table whose box contains it (the MDI
# is checked before, in in_mdi). Barrel and endcaps are two boxes. To add a subdetector: one line.
SUBDETECTOR_BOXES = {
    "vtx":     [(13.0, 330.0, 0.0, 964.0)],                   # VTX_r_min, VTX_r_max, VTXD_z_max_total
    "dch":     [(350.0, 2000.0, 0.0, 3000.0)],                # DCH gas volume radii, DCH_half_length_total (as in the DCH-only version)
    "siwr":    [(2040.0, 2080.0, 0.0, 3150.0),                # barrel: SiWrB_*
                (350.0, 2040.0, 3020.0, 3090.0)],             # disks: SiWrD_*
}
# own region of a hit collection (the subdetector its hits are in), by collection-name prefix (ddsim names of
# ALLEGRO_o1_v03). The own region plays the part the DCH plays for DCHCollection: its vertices get the label/case
# <own region> (and plot_merged_maps.py does not draw them); every other subdetector region is "other".
# Only the tracker-region subdetectors whose hits can be traced back with the sample (VTX, DCH, SiWr); tracker
# hit collections of other names are not traced.
COLLECTION_REGIONS = (("DCH", "dch"), ("Vertex", "vtx"), ("SiWr", "siwr"))

# vertex cases on the path of a hit: chain (ancestor) or older sibling (kink), in the MDI or elsewhere
# (the plotted cases of plot_merged_maps.py)
HIT_CASES = ("hit_ancestor_mdi", "hit_sibling_mdi", "hit_ancestor_other", "hit_sibling_other")

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
                    help="hit collection of the detailed printed summary; always traced")
    ap.add_argument("--hit-collections", default="all",
                    help="tracker hit collections to trace and keep in the skim: 'all' (every SimTrackerHit "
                         "collection of the VTX, DCH and SiWr in the file, COLLECTION_REGIONS) or a comma-separated "
                         "list of names")
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


def find_hit_collections(tree, wanted="all"):
    """tracker hit collections of the events tree to trace (sorted list of names), from the branch types:
    SimTrackerHit collections with an own region (VTX, DCH, SiWr: COLLECTION_REGIONS).
    wanted: "all" or a comma-separated list of collection names (an unknown or untraceable name raises).
    """
    found = sorted(name for name, typ in tree.typenames().items()
                   if "/" not in name and "SimTrackerHitData" in typ and own_region(name) is not None)
    if wanted != "all":
        names = [w for w in wanted.split(",") if w]
        missing = [w for w in names if w not in found]
        if missing:
            raise RuntimeError(f"no such VTX / DCH / SiWr tracker hit collection {missing}, the file has {found}")
        found = names
    return found


def read_events(path, n_events, colls):
    """colls: list of tracker hit collection names from find_hit_collections"""
    tree = uproot.open(path)["events"]
    branches = {key: f"MCParticles/MCParticles.{name}" for key, name in MC_FIELDS.items()}
    branches |= {
        "par_idx": "_MCParticles_parents/_MCParticles_parents.index",
        "dau_idx": "_MCParticles_daughters/_MCParticles_daughters.index",
    }
    # hits of every collection: particle index, energy, quality, time (for the kinks of a primary before its hit)
    hit_branches = {coll: {
        "p": f"_{coll}_particle/_{coll}_particle.index",
        "e": f"{coll}/{coll}.eDep",
        "q": f"{coll}/{coll}.quality",
        "t": f"{coll}/{coll}.time",
    } for coll in colls}
    wanted = list(branches.values()) + [b for d in hit_branches.values() for b in d.values()]
    arrays = tree.arrays(wanted, entry_stop=n_events, library="np")
    for e in range(len(arrays[branches["pdg"]])):
        ev = {key: arrays[b][e] for key, b in branches.items()}
        ev["hits"] = {coll: {key: arrays[b][e] for key, b in d.items()} for coll, d in hit_branches.items()}
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

def in_boxes(x, y, z, boxes):
    r = math.hypot(x, y)
    return any(rmin <= r <= rmax and zmin <= abs(z) <= zmax for rmin, rmax, zmin, zmax in boxes)

def region_of(ev, i):
    """fine region of the vertex of particle i: injection, mdi, a key of SUBDETECTOR_BOXES (the first
    in the table whose box contains it), or other"""
    if not created_in_sim(ev, i):
        return "injection"          # primary: its vertex is the SR source, not an interaction
    x, y, z = ev["x"][i], ev["y"][i], ev["z"][i]
    if in_mdi(x, y, z):
        return "mdi"
    for name, boxes in SUBDETECTOR_BOXES.items():
        if in_boxes(x, y, z, boxes):
            return name
    return "other"

def own_region(coll):
    """subdetector region of the hits of collection coll (COLLECTION_REGIONS), None if the name is not known"""
    for prefix, region in COLLECTION_REGIONS:
        if coll.startswith(prefix):
            return region
    return None

def coarse_label(region, own):
    """label relative to the hit collection: injection / mdi / own (the region of the hit collection, e.g. dch) / other
    (every other subdetector, and outside all of them)"""
    return region if region in ("injection", "mdi") or region == own else "other"


# ---------------------------------------------------------------------------
# Helper function: get older siblings
# ---------------------------------------------------------------------------
def daughters_before(ev, q, t):
    """daughters of q born before time t [ns], sorted by time"""
    d = daughters(ev, q)
    d = d[ev["t"][d] < t]
    return d[np.argsort(ev["t"][d], kind="stable")]

def older_siblings(ev, i, q):
    """daughters of q born before particle i, sorted by time"""
    return daughters_before(ev, q, ev["t"][i])

def is_sr_primary(ev, i):
    """the SR primary photon (generator particle, not created in the simulation)"""
    return not created_in_sim(ev, i) and ev["gen"][i] == 1 and ev["pdg"][i] == 22


# ---------------------------------------------------------------------------
# Helper function: create trees
# ---------------------------------------------------------------------------
def make_node(ev, i, role, own):
    """own: region of the hit collection of this tree (own_region)"""
    region = region_of(ev, i)
    return dict(
        idx=int(i), pdg=int(ev["pdg"][i]), role=role, region=region, label=coarse_label(region, own), case=None,
        x=float(ev["x"][i]), y=float(ev["y"][i]), z=float(ev["z"][i]), t=float(ev["t"][i]),
        ekin_kev=ekin_kev(ev, i),
        component=None,       # TODO: geometry component (MC vertices carry no volume name)
        parent=None, children=[], siblings_before=[], walked=False, root_status=None,
        n_hits=0, edep_kev=0.0, n_dropped=0,   # own hits (hit particles only)
        primary_hits=0, kinks_before_hit=[],   # SR primary only: own real hits, own daughters born before them
        primary_kink=False,                    # a kink of an SR primary before one of its own real hits
    )

def grow_tree(ev, i, edep_kev, dropped, t_hit, nodes, own):
    """Walk from hit particle i up to the root (the while-loop).

    own: region of the hit collection (own_region), for the labels.
    t_hit: time of the hit [ns], for the kinks of a primary before its own hit.
    nodes: {MCParticle index: node} of this event and collection, shared by all its hits. When the
    walk reaches a node that was already walked (an existing leg), it stops:
    the rest of the leg up to the root is already there.
    """
    if i not in nodes:
        nodes[i] = make_node(ev, i, "chain", own)
    leaf = nodes[i]
    leaf["n_hits"] += 1
    leaf["edep_kev"] += edep_kev
    leaf["n_dropped"] += dropped

    # real hit of the SR primary itself: its own daughters born before the hit are its kinks on the way to the hit
    # (an attributed hit, producedBySecondary, was made by a dropped secondary: no kinks, it stays unres_dropped)
    if not dropped and is_sr_primary(ev, i):
        leaf["primary_hits"] += 1
        for s in daughters_before(ev, i, t_hit):
            s = int(s)
            if s not in nodes:
                nodes[s] = make_node(ev, s, "sibling", own)
                nodes[s]["parent"] = i
            nodes[s]["primary_kink"] = True
            if s not in leaf["kinks_before_hit"]:
                leaf["kinks_before_hit"].append(s)

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
                    nodes[s] = make_node(ev, s, "sibling", own)
                    nodes[s]["parent"] = q
                node["siblings_before"].append(s)
        if q not in nodes:
            nodes[q] = make_node(ev, q, "chain", own)
        if i not in nodes[q]["children"]:
            nodes[q]["children"].append(i)
        i = q


def vertex_case(node):
    """classification of one vertex (particle) of a tree"""
    if node["root_status"] == "broken":
        return "unres_dropped"      # born in the simulation, but its parent was not stored
    if node["label"] in ("mdi", "other"):   # other = outside MDI and the own subdetector (e.g. backscatter)
        return f"hit_{'ancestor' if node['role'] == 'chain' else 'sibling'}_{node['label']}"
    return node["label"]                    # own region (e.g. dch) or injection


def split_trees(nodes, file_name, event, coll):
    """group the nodes of one event (one hit collection) by the root of their tree (one tree per root photon)"""
    def root_of(i):
        while nodes[i]["parent"] is not None:
            i = nodes[i]["parent"]
        return i

    trees = {}
    for i, n in nodes.items():
        n["case"] = vertex_case(n)
        r = root_of(i)
        t = trees.setdefault(r, dict(file=file_name, event=event, coll=coll, root=r, root_status=nodes[r]["root_status"],
                                     n_primary_hits=nodes[r]["primary_hits"], nodes={}, cases=Counter()))
        t["nodes"][i] = n
        t["cases"][n["case"]] += 1
    return list(trees.values())


def hit_summary(trees):
    """per hit: which vertex cases are on its path (chain + older siblings).
    Real hits of the SR primary itself: primary_hit+<cases of its kinks before the hit>, or primary_hit_nokink
    (the kinks are those before any of its hits: for a primary with several hits, a kink between two
    of them counts for all of them)"""
    out = Counter()
    for t in trees:
        nodes = t["nodes"]
        for leaf in nodes.values():
            n_hits = leaf["n_hits"]
            if n_hits == 0:
                continue
            real = leaf.get("primary_hits", 0)      # (trees written before 2026-11-06 have no primary fields)
            if real:
                if n_hits > real:
                    out["unres_dropped"] += n_hits - real      # its attributed (producedBySecondary) hits
                cases = {nodes[s]["case"] for s in leaf["kinks_before_hit"]} & set(HIT_CASES)
                out["+".join(["primary_hit"] + sorted(cases)) if cases else "primary_hit_nokink"] += real
                continue
            chain = [leaf["idx"]]
            while nodes[chain[-1]]["parent"] is not None:
                chain.append(nodes[chain[-1]]["parent"])
            cases = {nodes[j]["case"] for j in chain} | {nodes[s]["case"] for c in chain for s in nodes[c]["siblings_before"]}
            if leaf["n_dropped"] or len(chain) == 1 or "unres_dropped" in cases:
                key = "unres_dropped"
            elif not cases & set(HIT_CASES):
                key = "unres_nomdi"     # path has only own-region (e.g. dch) / injection vertices
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


def skim_frame(frame, keep, colls):
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

    def relink(src, dst):
        """point dst (a new hit) to the new copy of the particle of src, if it has one"""
        p = src.getParticle()
        if p.isAvailable():                         # every linked particle is a tree node, so in keep
            dst.setParticle(muts[remap[p.getObjectID().index]])

    for coll in colls:
        out = edm4hep.SimTrackerHitCollection()
        for h in frame.get(coll):
            m = h.clone(False)
            relink(h, m)
            out.push_back(m)
        new.put(out, coll)

    if "EventHeader" in frame.getAvailableCollections():
        header = edm4hep.EventHeaderCollection()
        for h in frame.get("EventHeader"):
            header.push_back(h.clone())
        new.put(header, "EventHeader")
    new.put(mc, "MCParticles")                      # all other collections are dropped
    return new


def skim_file(path, keeps, colls, write_metadata):
    """keeps: list of kept MCParticle indices, one per event (same order as uproot).
    Returns (reader, [(frame, category)]); nothing is written, so a file that fails
    half way leaves no partial events in the output. Keep the reader alive until written."""
    reader = Reader(path)
    events = list(reader.get("events"))
    if len(events) != len(keeps):
        raise RuntimeError(f"podio reads {len(events)} events, uproot {len(keeps)}")
    frames = [(skim_frame(frame, keep, colls), "events") for frame, keep in zip(events, keeps)]
    for cat in ("runs", "metadata"):
        if cat in reader.categories and (cat == "runs" or write_metadata):
            frames += [(frame, cat) for frame in reader.get(cat)]
    return reader, frames


def write_frames(writer, frames):
    for frame, cat in frames:
        writer.write_frame(frame, cat)


# ---------------------------------------------------------------------------
# Process one file
# ---------------------------------------------------------------------------
def trace_file(f, args):
    """Trees and kept MCParticle indices of one file;
    returns (colls, trees, keeps, n_mc, n_no_particle), colls = the traced collections of this file
    and n_no_particle = Counter {collection: hits without particle link}"""
    name = os.path.basename(f)
    check_complete(f)
    colls = find_hit_collections(uproot.open(f)["events"], args.hit_collections)
    storage, keeps = [], []
    n_mc = 0
    n_no_particle = Counter()
    for e, ev in enumerate(read_events(f, args.n_events, colls)):
        # create the hit family trees of each collection (own nodes: a particle with hits in two collections is in both)
        trees = []
        for coll, d in ev["hits"].items():
            nodes = {}
            own = own_region(coll)                  # the labels are relative to the subdetector of the collection
            for h in range(len(d["p"])):
                i = int(d["p"][h])
                if i < 0:
                    n_no_particle[coll] += 1
                    continue
                dropped = bool((d["q"][h] >> BIT_PRODUCED_BY_SECONDARY) & 1)
                grow_tree(ev, i, 1e6 * float(d["e"][h]), dropped, float(d["t"][h]), nodes, own)
            trees += split_trees(nodes, name, e, coll)
        storage += trees

        # every particle in the trees of this event (all collections) -> retained MCParticles
        keeps.append(sorted({i for t in trees for i in t["nodes"]}))
        n_mc += len(ev["pdg"])
    return colls, storage, keeps, n_mc, n_no_particle


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    if args.hit_collections != "all":               # the plot collection is always traced
        names = [w for w in args.hit_collections.split(",") if w]
        if args.collection not in names:
            names.append(args.collection)
        args.hit_collections = ",".join(names)
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
    n_no_particle = Counter()                       # {collection: hits without particle link}
    all_colls = set()                               # traced collections over the processed files
    n_done = 0
    skipped = []                                    # (file, reason): truncated or corrupt inputs
    for k, f in enumerate(files):
        name = os.path.basename(f)
        # a file that fails anywhere is skipped as a whole: nothing of it enters storage or the skim
        out_file = None
        try:
            colls, trees, keeps, mc, no_part = trace_file(f, args)
            if not args.no_skim:
                reader, frames = skim_file(f, keeps, colls, write_metadata=(n_done == 0))
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
        all_colls |= set(colls)
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
    out = out_trees = os.path.join(args.output_dir, f"trees_{tag}.jsonl")
    with open(out, "w") as fo:
        for t in storage:
            fo.write(json.dumps(t, default=sorted) + "\n")
    print("wrote", out)

    # summary
    print(f"\ntrees {len(storage)}, {sum(len(t['nodes']) for t in storage)} vertices, {len(all_colls)} hit collections")
    summary = {}                                    # per collection, also written to collections_<tag>.json
    for coll in sorted(all_colls):
        ct = [t for t in storage if t["coll"] == coll]
        summary[coll] = dict(own_region=own_region(coll), trees=len(ct),
                             vertices=sum(len(t["nodes"]) for t in ct),
                             hits=sum(n["n_hits"] for t in ct for n in t["nodes"].values()),
                             hits_without_particle=n_no_particle[coll],
                             primary_hits=sum(t["n_primary_hits"] for t in ct),
                             trees_with_primary_hits=sum(t["n_primary_hits"] > 0 for t in ct),
                             primary_kinks=sum(n["primary_kink"] for t in ct for n in t["nodes"].values()),
                             regions=dict(Counter(n["region"] for t in ct for n in t["nodes"].values())),
                             cases=dict(sum((Counter(t["cases"]) for t in ct), Counter())),
                             hit_paths=dict(hit_summary(ct)))
        s = summary[coll]
        print(f"  {coll} (own region {s['own_region']}): {s['trees']} trees, {s['vertices']} vertices, {s['hits']} hits"
              + (f", {s['hits_without_particle']} hits without particle link" if s["hits_without_particle"] else "")
              + f", {s['primary_hits']} real hits of SR primaries ({s['trees_with_primary_hits']} trees, "
                f"{s['primary_kinks']} primary kinks)")
        if s["trees"]:
            print(f"    vertex cases: {s['cases']}\n    hits by the cases on their path: {s['hit_paths']}")
    out = os.path.join(args.output_dir, f"collections_{tag}.json")
    with open(out, "w") as fo:
        json.dump(summary, fo, indent=1)
    print("wrote", out)

    # --collection (DCH) in detail
    dch = [t for t in storage if t["coll"] == args.collection]
    print(f"\n{args.collection}: {len(dch)} trees, {sum(len(t['nodes']) for t in dch)} vertices")
    print(f"  vertex cases:", dict(sum((t["cases"] for t in dch), Counter())))
    print(f"  {args.collection} hits by the cases on their path:", dict(hit_summary(dch)))
    print(f"files processed {n_done} of {len(files)}")
    print(f"MCParticles kept {n_keep} of {n_mc} ({100 * n_keep / max(n_mc, 1):.2f}%)")
    print(f"maps: python3 plot_merged_maps.py {out_trees} <plot dir>")


if __name__ == "__main__":
    main()
