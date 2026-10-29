#!/usr/bin/env python3
"""Merge many single-event edm4hep (podio ROOT) files into chunked output files.

Events are concatenated (one output event per input event). Collections matching
any --drop regex are not written, e.g. calorimeter / muon tagger hits.

The event data is copied at the ROOT level: the branches of the dropped collections
are disabled and the compressed baskets of the kept ones are copied as they are
(TTree::CopyEntries in "fast" mode), so dropped collections are never read or
decompressed. hadd / podio.root_merge.merge_files cannot drop collections, and the
podio_metadata (collection type table) must only list the kept collections, so it
is written by podio from a single reduced frame and copied into the output.

Run inside the Key4hep nightly environment first:
    source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh   # (alias sknight)

Examples:
    python3 merge_edm4hep_files.py --input-dir /path/to/files
    python3 merge_edm4hep_files.py --input-dir /path/to/files --max-files 5 --files-per-output 2 --output-dir /tmp/mergetest
"""

import argparse
import glob
import os
import re
import sys
import tempfile
import time


# Authenticate to EOS with the Kerberos ticket: with a grid certificate in ~/.globus,
# xrootd otherwise tries GSI first and asks for the PEM pass phrase
os.environ.setdefault("XrdSecPROTOCOL", "krb5,unix")

DEFAULT_DROP = []
META_TREES = ("metadata", "podio_metadata")  # written from podio, not copied from the inputs


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-dir", default=None, help="directory with the input files")
    parser.add_argument("--pattern", default="*_edm4hep_*.root", help="glob pattern for input files")
    parser.add_argument("--output-dir", default=None, help="output directory (default: <input-dir>/../merged)")
    parser.add_argument(
        "--prefix",
        default=None,
        help="output file name prefix (default: ALLEGRO_o1_v03_SR_halo_merged, with _noMCP appended for --no-mcparticles)",
    )
    parser.add_argument(
        "--files-per-output", type=int, default=0, help="input files per output file (0 = all into one file)"
    )
    parser.add_argument(
        "--drop",
        action="append",
        default=None,
        metavar="REGEX",
        help=f"drop event collections matching this regex (repeatable, default: {DEFAULT_DROP})",
    )

    parser.add_argument(
        "--no-mcparticles",
        action="store_true",
        help="also drop MCParticles (~98%% of the file size); hit->MCParticle links become unavailable",
    )
    parser.add_argument("--max-files", type=int, default=0, help="only use the first N input files (for testing)")

    args = parser.parse_args()
    if args.input_dir is None:
        parser.error("--input-dir is required")
    if args.drop is None:
        args.drop = list(DEFAULT_DROP)
    if args.no_mcparticles:
        args.drop.append(r"^MCParticles$")
    if args.prefix is None:
        args.prefix = "ALLEGRO_o1_v03_SR_halo_merged" + ("_noMCP" if args.no_mcparticles else "")
    if args.output_dir is None:
        args.output_dir = os.path.join(os.path.dirname(os.path.abspath(args.input_dir)), "merged")
    return args


def run_number(path):
    """Sort key: the run number in '..._<N>_edm4hep_...', falling back to the file name."""
    match = re.search(r"_(\d+)_edm4hep", os.path.basename(path))
    return (int(match.group(1)) if match else float("inf"), os.path.basename(path))


def select_collections(names, drop):
    drop_re = [re.compile(p) for p in drop]
    kept, dropped = [], []
    for name in sorted(names):
        if not any(r.search(name) for r in drop_re):
            kept.append(name)
        else:
            dropped.append(name)
    return kept, dropped


def read_url(path):
    """Read EOS user files through xrootd (about 2x faster per file than the fuse mount)."""
    if path.startswith("/eos/user/") or path.startswith("/eos/home-"):
        return "root://eosuser.cern.ch/" + path
    return path


def open_root_file(path):
    """Open an input file through xrootd, falling back to the plain path. Returns None if unreadable."""
    for url in dict.fromkeys((read_url(path), path)):
        try:
            tfile = ROOT.TFile.Open(url)
        except OSError:
            continue
        if tfile and not tfile.IsZombie():
            return tfile
    return None


def tree_names(tfile):
    return set(k.GetName() for k in tfile.GetListOfKeys() if k.GetClassName() == "TTree")


def branch_names(tree):
    return set(b.GetName() for b in tree.GetListOfBranches())


def collection_of_branch(branch, collections):
    """podio branch naming: collection X owns 'X', 'X_objIdx' (subset collections) and '_X_<member>'."""
    if branch in collections:
        return branch
    if branch.endswith("_objIdx") and branch[: -len("_objIdx")] in collections:
        return branch[: -len("_objIdx")]
    if branch.startswith("_"):
        owners = [c for c in collections if branch.startswith(f"_{c}_")]
        if owners:
            return max(owners, key=len)
    return None  # e.g. the frame parameter branches GP*


def read_layout(path, drop):
    """Reference layout from the first file: trees, event branches, kept/dropped collections and branches."""
    reader = root_io.Reader(path)
    events = list(reader.get("events"))
    if not events:
        sys.exit(f"First file {path} has no events")
    collections = set(str(c) for c in events[0].getAvailableCollections())
    kept, dropped = select_collections(collections, drop)
    del events, reader

    tfile = open_root_file(path)
    if tfile is None:
        sys.exit(f"Cannot open first file {path}")
    trees = tree_names(tfile)
    branches = branch_names(tfile.Get("events"))
    compression = tfile.GetCompressionSettings()
    tfile.Close()
    dropped_branches = sorted(b for b in branches if collection_of_branch(b, collections) in dropped)
    return {
        "trees": trees,
        "branches": branches,
        "kept": kept,
        "dropped": dropped,
        "dropped_branches": dropped_branches,
        "compression": compression,
    }


def merge_chunk(chunk, out_path, layout):
    """Merge the files of one chunk into out_path, removing the partial output on failure."""
    try:
        return _merge_chunk(chunk, out_path, layout)
    except BaseException:
        if os.path.exists(out_path):
            os.remove(out_path)
        raise


def _merge_chunk(chunk, out_path, layout):
    """Returns (used_files, skipped[(file, reason)], n_events)."""
    out_file = None
    out_trees = {}
    used, skipped = [], []
    n_events = 0
    copy_trees = sorted(layout["trees"] - set(META_TREES))

    for i_file, path in enumerate(chunk):
        if i_file and i_file % 100 == 0:
            print(f"    {i_file}/{len(chunk)} files read, {n_events} events", flush=True)
        tfile = open_root_file(path)
        if tfile is None:  # corrupt / unreadable file
            skipped.append((path, "cannot open"))
            continue
        try:
            trees = tree_names(tfile)
            if trees != layout["trees"]:
                skipped.append((path, f"trees {sorted(trees)} != {sorted(layout['trees'])}"))
                continue
            in_events = tfile.Get("events")
            if in_events.GetEntries() == 0:
                skipped.append((path, "no events"))
                continue
            if branch_names(in_events) != layout["branches"]:
                skipped.append((path, "event branch layout differs from the first file"))
                continue

            if out_file is None:
                out_file = ROOT.TFile(out_path, "RECREATE", "", layout["compression"])
            for name in copy_trees:
                in_tree = tfile.Get(name)
                if name == "events":
                    for b in layout["dropped_branches"]:
                        in_tree.SetBranchStatus(b, 0)
                        if in_tree.GetBranch(b).GetListOfBranches().GetEntries():
                            in_tree.SetBranchStatus(f"{b}.*", 0)  # split sub-branches are not disabled with the parent
                if name not in out_trees:
                    out_file.cd()
                    out_trees[name] = in_tree.CloneTree(0)  # empty tree with the active branches only
                n_in = in_tree.GetEntries()
                n_before = out_trees[name].GetEntries()
                out_trees[name].CopyEntries(in_tree, -1, "fast")  # returns bytes, not entries
                n_out = out_trees[name].GetEntries() - n_before
                if n_out != n_in:
                    raise RuntimeError(f"{path}: copied {n_out} of {n_in} '{name}' entries")
                if name == "events":
                    n_events += n_out
        finally:
            tfile.Close()
        used.append(path)

    if out_file is None:
        return used, skipped, 0

    out_file.cd()
    for tree in out_trees.values():
        tree.Write("", ROOT.TObject.kOverwrite)
    write_podio_metadata(out_file, used, layout)
    out_file.Close()
    return used, skipped, n_events


def write_podio_metadata(out_file, used, layout):
    """Write one reduced frame per category with podio to a temporary file, so that its
    podio_metadata lists only the kept collections, then copy podio_metadata (and the
    metadata frame of the first input, with MergedInputFiles) from there into the output."""
    first = root_io.Reader(read_url(used[0]))
    metadata_frames = list(first.get("metadata")) if "metadata" in layout["trees"] else []

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = os.path.join(tmp_dir, "podio_metadata.root")
        writer = root_io.Writer(tmp_path)
        for cat in sorted(layout["trees"] - set(META_TREES)):
            frames = first.get(cat, coll_names=layout["kept"]) if cat == "events" else first.get(cat)
            writer.write_frame(next(iter(frames)), cat)
        for frame in metadata_frames:
            frame.put_parameter("MergedInputFiles", [os.path.basename(p) for p in used])
            writer.write_frame(frame, "metadata")
        writer.finish()
        del writer, first

        tmp_file = ROOT.TFile.Open(tmp_path)
        for name in META_TREES:
            tree = tmp_file.Get(name)
            if not tree:
                continue  # inputs without a metadata tree
            out_file.cd()
            tree.CloneTree(-1, "fast").Write("", ROOT.TObject.kOverwrite)
        tmp_file.Close()


def main():
    args = parse_args()

    files = sorted(glob.glob(os.path.join(args.input_dir, args.pattern)), key=run_number)
    if args.max_files > 0:
        files = files[: args.max_files]
    if not files:
        sys.exit(f"No files matching {args.pattern} in {args.input_dir}")

    # Reference layout from the first file
    layout = read_layout(files[0], args.drop)
    if not layout["kept"]:
        sys.exit("All collections would be dropped, check --drop / --no-mcparticles")
    dropped = layout["dropped"]

    size = args.files_per_output if args.files_per_output > 0 else len(files)
    chunks = [files[i : i + size] for i in range(0, len(files), size)]
    width = max(3, len(str(len(chunks) - 1)))
    out_paths = [os.path.join(args.output_dir, f"{args.prefix}_part{i:0{width}d}.root") for i in range(len(chunks))]

    print(f"Input : {len(files)} files from {args.input_dir}")
    print(f"Output: {len(chunks)} file(s) of up to {size} inputs in {args.output_dir}")
    print(f"Trees : {sorted(layout['trees'])}")
    print(f"Dropping {len(dropped)} collections: {', '.join(dropped)}")
    print(f"Dropping {len(layout['dropped_branches'])} event branches: {', '.join(layout['dropped_branches'])}")

    os.makedirs(args.output_dir, exist_ok=True)
    # Replace any previous run with the same prefix, including parts beyond the new chunk count
    old_parts = sorted(glob.glob(os.path.join(args.output_dir, f"{glob.escape(args.prefix)}_part*.root")))
    for path in old_parts:
        os.remove(path)
    if old_parts:
        print(f"Removed {len(old_parts)} output file(s) from a previous run")
    all_skipped = []
    total_events = total_used = bytes_in = bytes_out = 0
    t_start = time.time()

    for i, (chunk, out) in enumerate(zip(chunks, out_paths)):
        t0 = time.time()
        used, skipped, n_events = merge_chunk(chunk, out, layout)
        for path, reason in skipped:
            print(f"  WARNING: skipped {os.path.basename(path)}: {reason}")
        all_skipped.extend(skipped)
        if not used:
            print(f"[{i + 1}/{len(chunks)}] {os.path.basename(out)}: no valid inputs, nothing written")
            continue
        total_events += n_events
        total_used += len(used)
        bytes_in += sum(os.path.getsize(p) for p in used)
        bytes_out += os.path.getsize(out)
        print(
            f"[{i + 1}/{len(chunks)}] {os.path.basename(out)}: {n_events} events from {len(used)} files, "
            f"{os.path.getsize(out) / 1e9:.2f} GB, {time.time() - t0:.0f} s"
        )

    print("\nSummary")
    print(f"  inputs merged : {total_used}")
    print(f"  inputs skipped: {len(all_skipped)}")
    for path, reason in all_skipped:
        print(f"    {path}: {reason}")
    print(f"  events written: {total_events}")
    if bytes_in:
        print(f"  size in/out   : {bytes_in / 1e9:.2f} GB -> {bytes_out / 1e9:.2f} GB ({100 * bytes_out / bytes_in:.1f}%)")
    print(f"  wall time     : {time.time() - t_start:.0f} s")


try:
    import ROOT
    from podio import root_io
except ImportError:
    sys.exit("Cannot import podio. Set up Key4hep first: source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh (sknight)")

if __name__ == "__main__":
    main()
