# Readout rates page

Static, read-only web page with FCC-ee readout-rate estimates per working point, detector, subdetector,
layer and sample (shown as Gbps, links, power). Not a package. Open questions: [OPEN_POINTS.md](OPEN_POINTS.md).

```
index.html, style.css, app.js      the page (all logic in app.js, no dependencies, no build step)
data/                              all numbers and assumptions, JSON (below), fetched by the page
tools/calc_ref.py                  reference rate model in Python (table per layer, or --check against the sheet)
tools/seed_from_sheet.py           one-off seed from the sheet "FCCee TDAQ numbers" (snapshot 2026-01-12)
tests/sheet_reference_ALLEGRO.json numbers read from the sheet, used by calc_ref.py --check
```

Edit a file and reload. The page fetches `data/*.json`, so it needs a web server (not `file://`); locally,
`python3 -m http.server` in this folder is enough. It is served from a clone of this repo in the EOS `www` area
(`/eos/user/a/aikoulou/www/`, URL `.../fcc-tdaq/`): edits made there are live at once, `git pull` updates it.

```bash
python3 tools/calc_ref.py --check                 # model vs sheet (passes within rounding)
python3 tools/calc_ref.py ALLEGRO_o1_v03 --sample LCC_V106
```

## Model

```
Gbps per layer = occupancy * rate[MHz] * bits * multiplier / 1000 * n_cells * SF
links          = Gbps / (line_rate * encoding_efficiency)   (the page also shows the sum rounded up per layer)
installed      = fibers_per_module * n_cells / cells_per_module   (page inputs: fibers/module per subdetector or per layer, cells/module per subdetector; defaults from detector.json)
power          = links * power_W_per_link
```

- **occupancy**: mean hits per cell per event, as a fraction (0.0265 = 2.65%), per sample, subdetector and layer.
  The page can use `p95` instead of `mean` for IPC and Z; injection always uses `mean`. `stdev` is stored but unused.
- **rate**: events per crossing-equivalent, in MHz, per source (below). IPC: one sample event per bunch crossing.
- **bits, multiplier**: readout size of one hit. Multiplier = cluster size for pixels, an inherited placeholder (3) elsewhere.
- **SF**: flat safety factor on everything.
- **Total** = sum over the enabled sources. *Peak* adds the injection window to the continuous sources (conservative,
  for sizing links); *average* weights injection by its duty (`turns_after_injection / f_rev / injection_period_s`, about 2.1%).

## Working point

Every sample belongs to a collider working point (energy). Only **Z** (Z-pole, 91.2 GeV) has samples; WW, ZH and
ttbar are listed in `working_points.json` for later. Z decays are the *physics* of the Z working point, not a background.

## data/

Occupancy numbers in the seed are rounded to 0.01% and come from the sheet; to be replaced by ROOT-derived values.

| file | content |
|---|---|
| `index.json` | `detectors`: list of `{id, label, status}`; each `id` is a folder under `data/` |
| `working_points.json` | `{id: {label, sqrt_s_GeV, beam_energy_GeV, f_bc_MHz, bunches_per_ring, f_rev_Hz, produces, status}}`; machine numbers are `null` where unknown (only Z is filled: 52 MHz, 11200, 3300 Hz) |
| `settings.json` | page defaults only: `working_point` (selected at load), `safety_factor`, `default_technology`, `default_samples {source: sample id}` |
| `sources.json` | per source (`ipc`, `sr`, `beamgas`, `touschek`, `injection`, `z`; `sr`, `beamgas` and `touschek` have no samples yet): `label`, `kind` (`background`/`physics`), `class` (`continuous`/`pulsed`), `timing` (informational), `rate_MHz` (`null` for IPC: it is `f_bc_MHz` of the working point), plus `pulse {turns_after_injection, injection_period_s}` and `scale {mean, peak}` for injection (null = not provided yet), `mixture [{sample, weight, note}]` and `default_selection` for `z` (the page can also pick one Z sample alone, counted as 100%), `components [{id, label}]` for `injection` and `sr` (core, halo; TBC), `working_point` for `z` |
| `technologies.json` | link presets `{id: {label, line_rate_Gbps, encoding_efficiency, power_W_per_link (null = unknown), verified, notes}}` |
| `penguin_facts.txt` | one fact per line (penguin currently commented out in the page) |
| `<detector>/detector.json` | `id, label, status, subdetectors[]` |
| `<detector>/samples.json` | `detector, samples[]` |

**subdetector** (`detector.json`):
`id, label, group, readout {model, basis, bits_per_unit, multiplier}, sensors_per_cell, channels_per_sensor, cells_per_module,
flags[], notes, layers [{id, n_cells, fibers_per_module}]`.
`readout.model` is `pixel_cluster`, `wire_hit` or `cell_hit`; the rate formula is the same for all, only bits and multiplier differ.
Layer `id` is a string key in the occupancy (negative = the other endcap side). `n_cells: null` = no geometry.

**sample** (`samples.json`):
`id, source (ipc|z|injection|sr), working_point, label, n_events, flags[], notes, occupancy`, optional `lattice` (beam
optics version), `process` (Z samples, e.g. Z->bb), `component` (injection/SR: core, halo, ...; null = unspecified), `scale {mean, peak}` (injection), `default` (selected at load).
The sidebar shows `process`, `lattice`, `component` and `n_events` under each selected sample.

```
"occupancy": { "<subdetector id>": { "<layer id>": { "mean": 0.2653, "stdev": 0.1393, "p95": 0.5343 } } }
```

A subdetector or layer missing from `occupancy` simply has no data for that sample.

**flags**:
`partial_layers` (not all layers present), `no_geometry`, `capacity_unreliable` (installed links undefined, hidden in the page),
`multiplier_inherited_from_sheet` (placeholder multiplier), `no_stats` (sample has too few events).

## Assumptions (all open, see OPEN_POINTS.md)

- IPC rate = `f_bc` of the working point (52 MHz at Z) and one event per crossing.
- Z physics rate 0.1 MHz, mixture 0.7 Zbb + 0.1 Zmumu (proxies; invisible decays contribute 0).
- Injection example peak: 52 MHz-equivalent, hardcoded in `app.js` (`injPeak`), not in the JSON until the scale factors exist.
- Injection: ~700 turns after each top-up, every 10 s, `f_rev` of the working point (3300 Hz at Z); mean/peak scale factors not available, the page uses an example peak.
- `cells_per_module = 9` for all subdetectors, 2 fibers per module; calo, drift chamber and muon capacity is therefore not meaningful.
- Only the fiber preset (5.12 Gbps, 8b10b) is verified; 1 GbE PoE and silicon photonics are placeholders, power unknown.
