import os

import ROOT
ROOT.gROOT.SetBatch(True)

METADATA_CATEGORIES = ("metadata", "meta")
MARKER_KEY = "nPrimariesTotal"


###### Functions to compute weight factor depending on total edm4hep file submissions
def read_metadata(path, strict=False, category=None):
    
    if not os.path.exists(path):
        if strict:
            raise FileNotFoundError(path)
        return {}

    try:
        from podio import root_io
        reader = root_io.Reader(path)
        available = [c for c in reader.categories if c != "events"]
        if category is not None:
            candidates = [category]
        else:
            # Preferred names first, then anything else the file happens to
            # have, so an unexpected rename still resolves.
            candidates = ([c for c in METADATA_CATEGORIES if c in available] +
                        [c for c in available if c not in METADATA_CATEGORIES])
        for name in candidates:
            frames = reader.get(name)
            if len(frames) == 0:
                continue
            # Bind the Frame to a name and build the dict while it is alive
            frame = frames[0]
            keys = list(frame.parameters)
            if category is None and MARKER_KEY not in keys:
                continue        # ddsim's own 'metadata' frame, not ours
            return {key: frame.get_parameter(key) for key in keys}
        raise KeyError(f"{path}: no metadata frame with '{MARKER_KEY}' "
                    f"(categories: {available})")
    except Exception:               # noqa: BLE001 -- not podio, unreadable, or no metadata
        if strict:
            raise
        return {}

def sr_weight(path, nevents):
    """
    The SR weight for `nevents` files of the simulation sample that `path` belongs to.

    Multiply any photon count by this to get counts per bunch crossing, one beam.

    `path` is any one file of the sample -- only its constants are used, so it
    does not matter which. `nevents` is how many events the analysis actually
    read, across all files. The metadata counts primaries, not events, so the
    per-event primary count from the file converts between them; every job in a
    sample uses the same ngenerate, which is what makes one file enough.

        w = sr_weight("halo_00123.root", chain.GetEntries())
    """
    if nevents <= 0:
        raise ValueError(f"nevents must be positive, got {nevents}")

    meta = read_metadata(path, strict=True)
    primaries_per_event = float(meta["nPrimariesTotal"]) / float(meta["nEvents"])

    weight = float(meta["chargeFraction"]) * float(meta["bunchIntensity"]) / (nevents * primaries_per_event)
    return (weight)


def print_normalisation(root_files, events_file):
    """
    Report how much background is being submitted, in event frames and in bunch crossings.
    """
    n_files = len(root_files)
    total_frames = events_file * n_files
    meta = read_metadata(root_files[0], category="metadata")

    print("-" * 72)
    print("edm4hep background normalisation")
    print(f"  events per file : {events_file}")
    print(f"  files submitted : {n_files}")
    print(f"  total frames    : {total_frames}")

    weight_line = None
    try:
        # sr_weight is counts per BX per particle, so its reciprocal is how many
        # BX these frames add up to.
        w = sr_weight(root_files[0], total_frames)

    except Exception:       # not an SR sample, or no weight metadata
        # The guinea-pig samples carry no SR constants. Use the generator's own
        # over-production factor instead, defaulting to one BX per frame.
        bw_scale = float(meta.get("bw_scale", 1.0))
        bx = total_frames * bw_scale
        if bw_scale != 1:
            weight_line = (f"per-BX weight   : 1/{bx:g} = {1.0 / bx:.6g}   "
                           f"(bw_scale={bw_scale:g}, so one frame is {bw_scale:g} BX)")
    else:
        bx = 1.0 / w
        weight_line = (f"sr_weight       : {w:.6g}   "
                       f"(chargeFraction={meta.get('chargeFraction')}, "
                       f"bunchIntensity={float(meta.get('bunchIntensity')):g}, "
                       f"nPrimariesTotal={meta.get('nPrimariesTotal')})")

    print(f"  BX covered      : {bx:g}")
    if weight_line:
        print(f"  {weight_line}")

    # OverlayTiming consumes one background frame per BX, and neither sample is
    # anywhere near one frame per BX, so spell out the real multiplicity.
    print(f"  frames per BX   : {total_frames / bx:.4g}" "   <- ideal OverlayTiming's NumberBackground value (if no jobs fail).")
    print(" ")
    print("-" * 72)

def count_events(path):
    """
    Number of podio event frames in the edm4hep file at `path`.

    This is what `ddsim -N` counts: one frame, however many particles it holds.
    The `events` tree is the authority here rather than the metadata's `nEvents`
    -- it is what ddsim will actually read, and the guinea-pig IPC samples carry
    no `nEvents` key at all.

    `check_metadata` additionally spins up a podio Reader to cross-check the tree
    against the metadata's `nEvents`.
    """
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise IOError(f"Could not open {path}")
    try:
        tree = f.Get("events")
        if not tree:
            raise KeyError(f"{path}: no 'events' tree -- not an edm4hep/podio file?")
        n = int(tree.GetEntries())
    finally:
        f.Close()

    return n