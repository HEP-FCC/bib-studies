#!/usr/bin/env python3
"""One-off seed of data/ from the "FCCee TDAQ numbers" sheet (snapshot of 2026-01-12).

The sheet is NOT the long-term source: later the numbers come from the ROOT files
(script TBD) or are typed by hand. This only gives the page something realistic to
show. Occupancies were transcribed from the displayed cells (rounded to 0.01 %).

Usage: seed_from_sheet.py        (rewrites data/ and tests/sheet_reference_ALLEGRO.json)
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
TESTS = os.path.normpath(os.path.join(HERE, "..", "tests"))

COMPONENTS = [{"id": "core", "label": "Core"}, {"id": "halo", "label": "Halo"}]
SHEET_NOTE = "seeded from sheet snapshot 2026-01-12, occupancies rounded to 0.01%"


def dump(obj, *path):
    p = os.path.join(*path)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w") as f:
        json.dump(obj, f, indent=1)
        f.write("\n")


def tri(s):
    """'mean/stdev/p95' in percent -> dict of fractions."""
    if s is None:
        return None
    m, sd, p = (float(x) / 100 for x in s.split("/"))
    return {"mean": m, "stdev": sd, "p95": p}


def mean_only(s):
    return None if s is None else {"mean": float(s) / 100, "stdev": None, "p95": None}


# layer tuple: (id, cells, IPC GHC_V25.3-4, IPC LCC_V106, injection mean %, Zbb, Zmumu, sheet IPC Gbps/L, fibers/module)
# percent strings are "mean/stdev/p95"; None = no value in the sheet
def L(i, cells, ghc, lcc, inj, zbb, zmm, sheet, fpm=2):
    return dict(id=i, cells=cells, ghc=ghc, lcc=lcc, inj=inj, zbb=zbb, zmm=zmm, sheet=sheet, fpm=fpm)


PIXEL = "pixel_cluster"
CELL = "cell_hit"
WIRE = "wire_hit"

ALLEGRO = [
    dict(id="VertexBarrel", group="Vertex", model=PIXEL, bits=18, mult=3, sens=2, ch=1000, layers=[
        L(0, 90, "18.21/8.74/33.80", "26.53/13.93/53.43", "0.00", "13.61/5.47/23.62", "2.48/1.08/3.34", 230.06),
        L(1, 240, "2.57/2.51/7.79", "3.21/3.45/9.84", "0.42", "6.47/2.84/11.93", "0.90/0.58/1.18", 86.60),
        L(2, 540, "0.64/0.91/2.54", "0.88/1.35/3.45", "0.74", "3.70/1.74/7.12", "0.43/0.28/0.64", 48.14),
        L(3, 736, "0.02/0.07/0.13", "0.04/0.16/0.14", "0.54", "2.32/1.32/4.80", "0.19/0.16/0.39", 2.27),
        L(4, 3264, "0.01/0.03/0.06", "0.01/0.03/0.06", "0.18", "0.57/0.41/1.30", "0.04/0.03/0.07", 5.96),
    ]),
    dict(id="VertexDisks", group="Vertex", model=PIXEL, bits=18, mult=3, sens=2, ch=1000, layers=[
        L(-1, 976, "0.11/0.15/0.39", "0.18/0.20/0.52", "2.05", "0.68/0.63/1.93", "0.07/0.10/0.24", 14.53),
        L(-2, 992, "0.05/0.09/0.30", "0.08/0.12/0.34", "1.82", "0.33/0.42/1.26", "0.03/0.07/0.19", 7.38),
        L(-3, 864, "0.05/0.10/0.24", "0.09/0.13/0.33", "2.08", "0.23/0.33/0.98", "0.01/0.05/0.13", 6.19),
        L(1, 976, "0.10/0.14/0.37", "0.18/0.20/0.53", "1.74", "0.64/0.59/1.86", "0.07/0.10/0.24", 13.43),
        L(2, 992, "0.05/0.09/0.20", "0.09/0.12/0.34", "2.02", "0.31/0.39/1.20", "0.03/0.07/0.19", 6.96),
        L(3, 864, "0.05/0.09/0.24", "0.09/0.13/0.33", "2.20", "0.22/0.32/0.88", "0.02/0.06/0.13", 5.58),
    ]),
    dict(id="SiWrD", group="Silicon wrapper", model=PIXEL, bits=18, mult=3, sens=1, ch=1000, layers=[
        L(-1, 1284, "0.42/0.21/0.81", "0.32/0.18/0.61", "0.55", "0.48/0.39/1.20", "0.03/0.05/0.08", 74.81),
        L(-2, 1244, "0.43/0.21/0.80", "0.32/0.19/0.69", "0.32", "0.47/0.39/1.28", "0.03/0.05/0.10", 75.80),
        L(1, 1284, "0.42/0.20/0.81", "0.31/0.18/0.62", "0.31", "0.47/0.38/1.19", "0.03/0.05/0.11", 75.17),
        L(2, 1244, "0.44/0.21/0.82", "0.33/0.20/0.71", "0.00", "0.47/0.38/1.24", "0.03/0.05/0.14", 75.98),
    ]),
    dict(id="SiWrB", group="Silicon wrapper", model=PIXEL, bits=18, mult=3, sens=1, ch=1000,
         notes="sheet todo: investigate why some hits are found in layers 2 and 3", layers=[
        L(0, 151, "4.72/1.82/8.18", None, "2.65", "9.11/4.77/17.35", "0.85/0.66/1.32", 100.15),
        L(1, 151, "4.55/1.88/8.14", None, "1.99", "8.85/4.75/17.23", "0.83/0.66/1.32", 96.36),
    ]),
    dict(id="DCH_v2", group="Drift chamber", model=WIRE, bits=800, mult=3, sens=1, ch=1000,
         flags=["partial_layers"], notes="only the first 10 of 112 layers; sheet todo: DCH occupancy excess for MTR, probably a bug", layers=[
        L(1, 192, "0.26/1.03/1.08", "0.16/0.45/0.99", "37.50", "9.49/5.72/21.12", "0.67/0.58/1.09", 305.51),
        L(2, 192, "0.25/1.04/1.09", "0.17/0.43/1.00", "25.00", "8.97/5.72/19.85", "0.65/0.51/1.09", 297.12),
        L(3, 192, "0.29/1.20/1.40", "0.16/0.45/0.96", "30.73", "9.12/5.62/19.75", "0.66/0.50/1.09", 347.44),
        L(4, 192, "0.30/1.23/1.38", "0.18/0.45/1.01", "25.00", "9.53/5.76/20.62", "0.68/0.43/1.09", 359.42),
        L(5, 192, "0.32/1.15/1.55", "0.17/0.50/0.99", "29.69", "9.66/5.83/21.04", "0.69/0.46/1.09", 378.59),
        L(6, 192, "0.32/1.21/1.46", "0.18/0.47/1.02", "28.13", "9.97/5.81/21.02", "0.70/0.43/1.09", 377.40),
        L(7, 192, "0.33/1.25/1.50", "0.19/0.44/1.02", "27.08", "10.25/6.06/22.03", "0.73/0.45/1.09", 398.96),
        L(8, 192, "0.34/1.23/1.58", "0.19/0.48/1.04", "19.27", "10.44/6.21/21.88", "0.77/0.46/1.09", 401.36),
        L(9, 240, "0.27/1.16/1.22", "0.12/0.37/0.80", "17.08", "7.02/4.60/15.58", "0.50/0.38/0.90", 404.35),
        L(10, 240, "0.25/1.05/1.30", "0.12/0.37/0.80", "6.25", "7.22/4.59/16.22", "0.50/0.38/0.90", 368.41),
    ]),
    dict(id="ECalBarrel", group="Calorimeter", model=CELL, bits=48, mult=3, sens=1, ch=16,
         notes="sheet todo: fix cell/channel # for Calo", layers=[
        L(0, 8880, "0.23/0.04/0.29", "0.28/0.03/0.35", "0.06", "0.09/0.03/0.15", "0.00/0.00/0.01", 764.67),
        L(1, 35616, "0.03/0.01/0.03", "0.04/0.00/0.04", "0.01", "0.05/0.02/0.09", "0.00/0.00/0.01", 400.04),
        L(2, 8880, "0.03/0.01/0.04", "0.03/0.01/0.04", "0.01", "0.20/0.09/0.34", "0.01/0.01/0.02", 96.42),
        L(3, 8880, "0.02/0.00/0.02", "0.01/0.00/0.02", "0.00", "0.23/0.11/0.40", "0.01/0.01/0.03", 49.87),
        L(4, 8880, "0.01/0.00/0.01", "0.01/0.00/0.01", "0.01", "0.24/0.11/0.41", "0.01/0.01/0.04", 29.92),
        L(5, 8880, "0.01/0.00/0.01", "0.00/0.00/0.01", "0.00", "0.23/0.12/0.41", "0.01/0.01/0.04", 16.62),
        L(6, 8880, "0.00/0.00/0.01", "0.00/0.00/0.01", "0.00", "0.21/0.11/0.38", "0.01/0.01/0.04", 9.97),
        L(7, 8880, "0.00/0.00/0.01", "0.00/0.00/0.01", "0.00", "0.19/0.10/0.35", "0.01/0.01/0.04", 6.65),
        L(8, 8880, "0.00/0.00/0.01", "0.00/0.00/0.01", "0.00", "0.18/0.10/0.33", "0.01/0.01/0.03", 3.32),
        L(9, 8880, "0.00/0.00/0.01", "0.00/0.00/0.01", "0.00", "0.16/0.09/0.31", "0.01/0.01/0.02", 3.32),
        L(10, 8880, "0.00/0.00/0.01", "0.00/0.00/0.01", "0.00", "0.16/0.10/0.31", "0.01/0.01/0.02", 3.32),
    ]),
    dict(id="EMEC_turbine", group="Calorimeter", model=CELL, bits=48, mult=3, sens=1, ch=1,
         notes="sheet todo: fix cell/channel # for Calo", layers=[
        L(1, 38400, "0.49/0.13/0.73", "0.07/0.03/0.11", "0.07", "0.23/0.48/1.14", "0.01/0.06/0.03", 7015.96),
        L(-1, 38400, "0.49/0.13/0.72", "0.06/0.02/0.10", "1.37", "0.25/0.55/1.33", "0.01/0.12/0.03", 7030.33),
        L(2, 100800, "0.04/0.01/0.06", "0.02/0.01/0.03", "0.00", "0.22/0.37/1.16", "0.01/0.04/0.06", 1509.58),
        L(-2, 100800, "0.04/0.01/0.06", "0.02/0.01/0.03", "0.05", "0.22/0.36/1.11", "0.01/0.05/0.05", 1509.58),
        L(3, 462400, "0.01/0.00/0.01", "0.00/0.00/0.01", "0.00", "0.09/0.12/0.37", "0.00/0.02/0.02", 865.61),
        L(-3, 462400, "0.01/0.00/0.01", "0.00/0.00/0.01", "0.01", "0.09/0.12/0.35", "0.00/0.01/0.02", 865.61),
    ]),
    dict(id="HCalBarrel", group="Calorimeter", model=CELL, bits=48, mult=3, sens=1, ch=100,
         notes="Zbb occupancies of 25-49% vs 0.2% for IPC look wrong (cell count / normalisation?)", layers=[
        L(0, 310, "0.16/0.51/1.03", "0.11/0.37/0.73", None, "48.90/31.11/95.30", "0.61/0.88/1.97", 17.99),
        L(1, 310, "0.09/0.36/0.69", "0.07/0.30/0.65", None, "44.54/28.83/91.57", "0.53/0.75/1.81", 10.91),
        L(2, 310, "0.06/0.25/0.35", "0.05/0.24/0.34", None, "37.74/25.66/81.14", "0.51/0.75/1.71", 7.20),
        L(3, 310, "0.05/0.21/0.33", "0.04/0.22/0.31", None, "31.19/22.09/67.61", "0.47/0.65/1.25", 5.22),
        L(4, 310, "0.06/0.29/0.34", "0.05/0.31/0.30", None, "36.44/26.19/80.40", "0.63/0.85/1.92", 7.08),
        L(5, 310, "0.05/0.30/0.01", "0.04/0.27/0.01", None, "27.14/20.88/64.33", "0.59/0.73/1.93", 5.34),
        L(6, 310, "0.04/0.22/0.01", "0.03/0.18/0.01", None, "19.82/16.27/49.32", "0.56/0.67/1.84", 4.18),
        L(7, 310, "0.03/0.17/0.01", "0.02/0.13/0.01", None, "14.66/12.86/38.62", "0.54/0.64/1.78", 3.13),
        L(8, 310, "0.02/0.11/0.01", "0.01/0.13/0.01", None, "10.82/9.98/30.79", "0.53/0.61/1.81", 1.86),
        L(9, 310, "0.01/0.09/0.01", "0.01/0.10/0.01", None, "7.93/7.77/23.84", "0.52/0.61/1.75", 1.51),
        L(10, 310, "0.02/0.11/0.01", "0.01/0.12/0.01", None, "7.83/7.90/24.34", "0.70/0.88/2.15", 2.09),
        L(11, 310, "0.01/0.07/0.01", "0.01/0.09/0.01", None, "4.01/5.02/14.17", "0.66/0.85/2.13", 0.93),
        L(12, 310, "0.01/0.08/0.01", "0.00/0.05/0.01", None, "2.26/3.63/9.40", "0.64/0.90/2.09", 0.70),
    ]),
    dict(id="HCalThreePartsEndcap", group="Calorimeter", model=CELL, bits=48, mult=3, sens=1, ch=1,
         flags=["no_geometry"], notes="no geometry in the sheet", layers=[
        L(-1, None, None, None, None, None, None, None),
        L(1, None, None, None, None, None, None, None),
    ]),
    dict(id="MuonTaggerBarrel", group="Muon tagger", model=CELL, bits=32, mult=3, sens=1, ch=64,
         notes="only injection occupancies in the sheet", layers=[
        L(0, 512, None, None, "2.20", None, None, None),
        L(1, 512, None, None, "2.35", None, None, None),
    ]),
    dict(id="MuonTaggerEndcap", group="Muon tagger", model=CELL, bits=32, mult=3, sens=1, ch=64,
         notes="only injection occupancies in the sheet; sheet todo: evaluate with reco muon positions", layers=[
        L(-2, 512, None, None, "0.02", None, None, None, fpm=2),
        L(-1, 512, None, None, "0.10", None, None, None, fpm=3),
        L(1, 512, None, None, "34", None, None, None, fpm=3),
        L(2, 512, None, None, "35", None, None, None, fpm=2),
    ]),
]

IDEA = [
    dict(id="VertexBarrel", group="Vertex", model=PIXEL, bits=18, mult=3, sens=2, ch=1000,
         notes="sheet todo: values differ considerably from ALLEGRO", layers=[
        L(0, 90, None, "18.48/8.31/33.03", None, None, None, None),
        L(1, 240, None, "2.81/2.45/7.81", None, None, None, None),
        L(2, 540, None, "0.86/1.03/2.98", None, None, None, None),
        L(3, 736, None, "0.03/0.11/0.13", None, None, None, None),
        L(4, 3264, None, "0.01/0.02/0.06", None, None, None, None),
    ]),
    dict(id="VertexDisks", group="Vertex", model=PIXEL, bits=18, mult=3, sens=2, ch=1000, layers=[
        L(-1, 976, None, "0.16/0.18/0.48", None, None, None, None),
        L(-2, 992, None, "0.07/0.10/0.32", None, None, None, None),
        L(-3, 864, None, "0.09/0.11/0.31", None, None, None, None),
        L(1, 976, None, "0.18/0.18/0.49", None, None, None, None),
        L(2, 992, None, "0.08/0.11/0.32", None, None, None, None),
        L(3, 864, None, "0.09/0.12/0.32", None, None, None, None),
    ]),
]


def detector_json(det_id, label, status, subs):
    out = []
    for s in subs:
        flags = list(s.get("flags", []))
        if s["model"] != PIXEL:
            flags.append("capacity_unreliable")
            if s["mult"] == 3:
                flags.append("multiplier_inherited_from_sheet")
        out.append({
            "id": s["id"],
            "label": s["id"],
            "group": s["group"],
            "readout": {"model": s["model"], "basis": "cell", "bits_per_unit": s["bits"], "multiplier": s["mult"]},
            "sensors_per_cell": s["sens"],
            "channels_per_sensor": s["ch"],
            "cells_per_module": 9,
            "flags": flags,
            "notes": s.get("notes", ""),
            "layers": [{"id": l["id"], "n_cells": l["cells"], "fibers_per_module": l["fpm"]} for l in s["layers"]],
        })
    return {"id": det_id, "label": label, "status": status, "subdetectors": out}


def occupancy(subs, key, conv=tri):
    occ = {}
    for s in subs:
        for l in s["layers"]:
            v = conv(l[key])
            if v is not None:
                occ.setdefault(s["id"], {})[str(l["id"])] = v
    return occ


def sample(sid, source, label, occ, **kw):
    d = {"id": sid, "source": source, "working_point": "Z", "label": label, "n_events": None, "flags": [], "notes": SHEET_NOTE}
    d.update(kw)
    d["occupancy"] = occ
    return d


def main():
    # ---- global files
    dump({
        "safety_factor": 5,
        "default_technology": "fiber_5g12_8b10b",
        "default_samples": {"ipc": "LCC_V106"},
        "working_point": "Z",
    }, DATA, "settings.json")

    dump({
        "ipc": {
            "label": "Incoherent pair creation (IPC)", "kind": "background", "class": "continuous",
            "timing": "per-BC (TBC)", "rate_MHz": 52,
        },
        "sr": {
            "label": "Synchrotron radiation (SR)", "kind": "background", "class": "continuous",
            "timing": "spread (TBC)", "rate_MHz": None, "components": COMPONENTS,
            "notes": "separate samples per component (core, halo, ...; TBC), none available yet",
        },
        "beamgas": {
            "label": "Beam-gas", "kind": "background", "class": "continuous", "timing": "spread (TBC)",
            "rate_MHz": None, "notes": "no samples yet",
        },
        "touschek": {
            "label": "Touschek", "kind": "background", "class": "continuous", "timing": "spread (TBC)",
            "rate_MHz": None, "notes": "no samples yet",
        },
        "injection": {
            "label": "Top-up injection", "kind": "background", "class": "pulsed",
            "pulse": {"turns_after_injection": 700, "injection_period_s": 10},
            "rate_MHz": None, "scale": {"mean": None, "peak": None}, "components": COMPONENTS,
            "notes": "scale factors are per sample and multiply the per-event result into mean / peak rates (TBC); core/halo components TBC (sample field `component`)",
        },
        "z": {
            "label": "Z decays (physics)", "working_point": "Z", "kind": "physics", "class": "continuous", "timing": "per-BC (TBC)",
            "rate_MHz": 0.1,
            "default_selection": "mixture",
            "mixture": [
                {"sample": "Zbb", "weight": 0.7, "note": "proxy for hadronic Z"},
                {"sample": "Zmumu", "weight": 0.1, "note": "proxy for charged-lepton Z"},
            ],
            "notes": "physics signal, not background; exists only at the Z-pole working point; the weights are hardcoded in the sheet, invisible Z contributes 0",
        },
    }, DATA, "sources.json")

    dump({
        "fiber_5g12_8b10b": {
            "label": "Fiber 5.12 Gbps (8b10b)", "line_rate_Gbps": 5.12, "encoding_efficiency": 0.8,
            "power_W_per_link": None, "verified": True, "notes": "from the sheet (rounded to 4.1 there)",
        },
        "gbe_poe_1g": {
            "label": "1 GbE PoE", "line_rate_Gbps": 1.0, "encoding_efficiency": 0.94,
            "power_W_per_link": None, "verified": False, "notes": "placeholder, framing efficiency assumed",
        },
        "sipho_100g": {
            "label": "Silicon photonics 100G", "line_rate_Gbps": 100.0, "encoding_efficiency": 1.0,
            "power_W_per_link": None, "verified": False, "notes": "placeholder",
        },
    }, DATA, "technologies.json")

    # ---- ALLEGRO
    dump(detector_json("ALLEGRO_o1_v03", "ALLEGRO o1 v03", "seeded", ALLEGRO), DATA, "ALLEGRO_o1_v03", "detector.json")
    samples = [
        sample("LCC_V106", "ipc", "IPC, LCC V106", occupancy(ALLEGRO, "lcc"), lattice="LCC_V106", default=True),
        sample("GHC_V25.3-4", "ipc", "IPC, GHC V25.3-4", occupancy(ALLEGRO, "ghc"), lattice="GHC_V25.3-4"),
        sample("Zbb", "z", "Z -> bb (hadronic proxy)", occupancy(ALLEGRO, "zbb"), process="Z->bb"),
        sample("Zmumu", "z", "Z -> mumu (leptonic proxy)", occupancy(ALLEGRO, "zmm"), process="Z->mumu"),
        sample("GHC_V25.1_inj", "injection", "Injection, GHC V25.1", occupancy(ALLEGRO, "inj", mean_only),
               lattice="GHC_V25.1", n_events=1, flags=["no_stats"], scale={"mean": None, "peak": None}),
    ]
    dump({"detector": "ALLEGRO_o1_v03", "samples": samples}, DATA, "ALLEGRO_o1_v03", "samples.json")

    # ---- IDEA (partial)
    dump(detector_json("IDEA_o1_v03", "IDEA o1 v03", "partial: vertex only, version TBC", IDEA), DATA, "IDEA_o1_v03", "detector.json")
    dump({"detector": "IDEA_o1_v03", "samples": [
        sample("LCC_V106", "ipc", "IPC, LCC V106", occupancy(IDEA, "lcc"), lattice="LCC_V106", default=True),
    ]}, DATA, "IDEA_o1_v03", "samples.json")

    dump({"detectors": [
        {"id": "ALLEGRO_o1_v03", "label": "ALLEGRO o1 v03", "status": "seeded"},
        {"id": "IDEA_o1_v03", "label": "IDEA o1 v03", "status": "partial"},
    ]}, DATA, "index.json")

    # ---- reference values from the sheet, for tests/check_against_sheet
    ref = {"snapshot": "2026-01-12", "ipc_GHC_V25.3-4_Gbps_per_layer": {}, "z_Gbps_per_layer": {"VertexBarrel": {
        "0": 0.24, "1": 0.30, "2": 0.38, "3": 0.33, "4": 0.35}}}
    for s in ALLEGRO:
        for l in s["layers"]:
            if l["sheet"] is not None:
                ref["ipc_GHC_V25.3-4_Gbps_per_layer"].setdefault(s["id"], {})[str(l["id"])] = l["sheet"]
    dump(ref, TESTS, "sheet_reference_ALLEGRO.json")
    print("seeded", DATA)


if __name__ == "__main__":
    main()
