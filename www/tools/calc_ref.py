#!/usr/bin/env python3
"""Reference implementation of the rate model (the page re-implements the same maths in JS).

  Gbps/layer = occupancy * rate[MHz] * bits_per_unit * multiplier / 1000 * n_cells * SF
  fibers_required = Gbps / (line_rate * encoding_efficiency)
  fibers_allowed  = fibers_per_module * n_cells / cells_per_module

Usage:
  calc_ref.py [DETECTOR] [--sample ID] [--tech ID] [--sf X]   print a per-layer table
  calc_ref.py --check                                           compare with the sheet snapshot
"""
import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
TESTS = os.path.normpath(os.path.join(HERE, "..", "tests"))


def load(*path):
    with open(os.path.join(*path)) as f:
        return json.load(f)


def sample_gbps(sub, layer, occ_mean, rate_MHz, sf):
    r = sub["readout"]
    return occ_mean * rate_MHz * r["bits_per_unit"] * r["multiplier"] / 1000 * layer["n_cells"] * sf


def link_gbps(tech):
    return tech["line_rate_Gbps"] * tech["encoding_efficiency"]


def source_rate(src, settings):
    """Rate in MHz; null in sources.json (IPC) means f_bc of the selected working point."""
    if src["rate_MHz"] is not None:
        return src["rate_MHz"]
    return load(DATA, "working_points.json")[settings["working_point"]]["f_bc_MHz"]


def layer_rows(det, samples, sources, settings, techs, sample_id, tech_id=None, sf=None):
    """Yield (sub, layer, gbps, fibers_required, fibers_allowed); a source mixture (Z) sums its samples."""
    sf = settings["safety_factor"] if sf is None else sf
    tech = techs[tech_id or settings["default_technology"]]
    by_id = {s["id"]: s for s in samples}
    smp = by_id[sample_id]
    src = sources[smp["source"]]
    comps = [(sample_id, 1.0)]
    if "mixture" in src:
        comps = [(m["sample"], m["weight"]) for m in src["mixture"]]
    for sub in det["subdetectors"]:
        for layer in sub["layers"]:
            if layer["n_cells"] is None:
                continue
            m = 0.0
            seen = False
            for sid, w in comps:
                o = by_id[sid]["occupancy"].get(sub["id"], {}).get(str(layer["id"]))
                if o:
                    m += w * o["mean"]
                    seen = True
            if not seen:
                continue
            g = sample_gbps(sub, layer, m, source_rate(src, settings), sf)
            allowed = layer["fibers_per_module"] * layer["n_cells"] / sub["cells_per_module"]
            yield sub, layer, g, g / link_gbps(tech), allowed


def check():
    settings = load(DATA, "settings.json")
    sources = load(DATA, "sources.json")
    techs = load(DATA, "technologies.json")
    det = load(DATA, "ALLEGRO_o1_v03", "detector.json")
    samples = load(DATA, "ALLEGRO_o1_v03", "samples.json")["samples"]
    ref = load(TESTS, "sheet_reference_ALLEGRO.json")
    # IPC: the sheet formula points at the GHC_V25.3-4 column; Z: 0.7 Zbb + 0.1 Zmumu
    checks = [("IPC GHC_V25.3-4", "GHC_V25.3-4", ref["ipc_GHC_V25.3-4_Gbps_per_layer"], source_rate(sources["ipc"], settings)),
              ("Z (Zbb/Zmumu mixture)", "Zbb", ref["z_Gbps_per_layer"], source_rate(sources["z"], settings))]
    bad = 0
    for name, sid, expected, rate in checks:
        n = 0
        worst = 0.0
        for sub, layer, g, _, _ in layer_rows(det, samples, sources, settings, techs, sid):
            exp = expected.get(sub["id"], {}).get(str(layer["id"]))
            if exp is None:
                continue
            n += 1
            # occupancies are transcribed at 0.01 % precision: allow that rounding band + 1 %
            band = sample_gbps(sub, layer, 0.5e-4, rate, settings["safety_factor"])
            dev = abs(g - exp)
            worst = max(worst, dev / max(exp, 1e-9))
            if dev > band + 0.01 * exp + 0.005:
                bad += 1
                print("  MISMATCH %s %s/%s: got %.3f, sheet %.3f (band %.3f)" % (name, sub["id"], layer["id"], g, exp, band))
        print("%s: %d layers compared, worst relative deviation %.1f%%" % (name, n, 100 * worst))
    print("FAIL" if bad else "OK")
    return 1 if bad else 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("detector", nargs="?", default="ALLEGRO_o1_v03")
    ap.add_argument("--sample")
    ap.add_argument("--tech")
    ap.add_argument("--sf", type=float)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    if a.check:
        sys.exit(check())
    settings = load(DATA, "settings.json")
    sources = load(DATA, "sources.json")
    techs = load(DATA, "technologies.json")
    det = load(DATA, a.detector, "detector.json")
    samples = load(DATA, a.detector, "samples.json")["samples"]
    sid = a.sample or settings["default_samples"]["ipc"]
    print("%-22s %6s %12s %12s %12s" % ("subdetector", "layer", "Gbps", "fib req", "fib allowed"))
    tot = 0.0
    for sub, layer, g, fr, fa in layer_rows(det, samples, sources, settings, techs, sid, a.tech, a.sf):
        tot += g
        print("%-22s %6s %12.2f %12.1f %12.1f" % (sub["id"], layer["id"], g, fr, fa))
    print("total %.1f Gbps (%.2f Tbps)" % (tot, tot / 1000))


if __name__ == "__main__":
    main()
