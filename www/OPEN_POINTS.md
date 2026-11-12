# Readout rates page: open points

## Context (read first when resuming)
- **What**: static page estimating FCC-ee readout bandwidth, links and power per subdetector, layer and sample; replaces the Google Sheet "FCCee TDAQ numbers" (id `13tztOHZzISJk8RRxsn5gES2AxukquaA_7dGI2ui1mKU`, snapshot of 2026-01-12 in `tools/seed_from_sheet.py` and `tests/`). Model, JSON formats and assumptions: [README.md](README.md).
- **Where**: live at https://milin.web.cern.ch/fcc-tdaq/, served from the git clone `/eos/user/a/aikoulou/www/bib-studies` (origin `HEP-FCC/bib-studies`, folder `www/`) through the symlink `/eos/user/a/aikoulou/www/fcc-tdaq -> bib-studies/www`. Edits there are public at once; there is no build step. Commit and push only when asked (one-line messages). Local dev checkout `fcc_workdir/bib-studies` is a second clone: `git pull` there.
- **Sibling page**: ROOT-file browser at https://milin.web.cern.ch/fcc/ (repo `meleneemil/fccweb`, `/eos/user/a/aikoulou/www/fcc`): `manifest.json` from `make_manifest.py`, folder texts in `root-files/<folder>/info.json` (draft descriptions, `"draft": true`, to be corrected by the author). Both pages share a dark/light/auto toggle (`localStorage` key `fcc-theme`).
- **Check after model changes**: `python3 tools/calc_ref.py --check` (reproduces the sheet within rounding; Z within 1.5%).
- **State**: ALLEGRO complete as far as the sheet goes; IDEA has the vertex detector only. All numbers are seeded from the sheet (see Data quality), none yet derived from ROOT files.
- Retired (2026-11-12): old single-file build `fcc_workdir/tdaq-rates/` (src/, dist/, build_page.py) and `www/fcc-tdaq.bak`; the Artifact preview is no longer maintained.

Living list. Strike through or delete items as they are settled; add the date and the decision when you do.

## Decisions taken
- Z is a **physics** sample, not a background (`kind: "physics"` in `data/sources.json`); backgrounds are `kind: "background"`.
- Fibers allowed per layer = `fibers_per_module * n_cells / cells_per_module` (the sheet used `fibers_per_module * n_cells` outside the VertexBarrel, probably a copy-paste). Correct later if needed.
- Default IPC sample is LCC_V106.
- Samples belong to a collider **working point** (`working_point` in `samples.json`, `data/working_points.json`); only Z-pole (91.2 GeV) exists. The physics source is labelled "Z decays", not "Z-pole" (decided 2026-11-10).
- Gbps/fiber is a setting derived from a technology preset (`line_rate * encoding_efficiency`; fiber = 5.12 * 0.8 = 4.096, the sheet rounded to 4.1).
- Static page, read only at first, ALLEGRO first, IDEA-ready (one folder per detector under `data/`).
- Existing ROOT files of `www/fcc` were moved into `www/fcc/root-files/` on 2026-11-12 (old shared JSROOT links broke on purpose); the ROOT-files page is static (`index.php` removed).
- Hosting settled (2026-11-12): plain files, no build, git clone in the EOS www area, nothing sensitive in the files (browsable). The penguin is commented out in `index.html`/`app.js`; `data/penguin_facts.txt` is kept.

## Rate model
- [ ] Per-subdetector multiplier / formula: calo, DCH and muon inherit `multiplier = 3` (the pixel cluster size) from the sheet (flag `multiplier_inherited_from_sheet`). Each needs its own.
- [ ] Installed capacity for calo, DCH and muon is undefined (flag `capacity_unreliable`): the sheet's `fibers_per_module * cells / 9` gives absurd link counts there. The page hides capacity, the over-capacity flag and the installed-link count for them until real link counts exist.
- [ ] Fibers/module and cells/module are editable in the page (per subdetector, fibers also per layer); the JSON defaults (2 and 9) are still the sheet's placeholders.
- [ ] `cells_per_module = 9` is hardcoded in the sheet for every subdetector, including calo and DCH.
- [ ] Calo cell/channel numbers (sheet todo "fix cell/channel # for Calo").
- [ ] Occupancy is per **cell** in the formula, not per channel; the channels/sensor are not used for the rate. Confirm.
- [ ] Mean vs 95th percentile: p95 is selectable in the page for IPC and Z (done); injection always uses the mean; stdev is not used.

## Backgrounds and time structure
- [ ] Continuous sources: per-BC (IPC) vs spread (SR, beam-gas). Stored as `timing` only, not used in any calculation yet.
- [ ] Peak aggregation across pulsed sources: add as if coincident (conservative default), or `max`, or an aggregation window / period parameter?
- [ ] Injection: `scale {mean, peak}` per sample (provided by the injection people) multiplies the per-event result. Values are `null` until we have them. Window = ~700 turns after each injection, every ~10 s (duty ~2% with f_rev = 3300 Hz, to be confirmed).
- [ ] Injection samples come with "1 evt" (that is the file, not the statistics of the physics); two sets exist (GHC_V25.1 occupancy in the sheet, LCC_V106 `max #hits` / `<#hits>`, mostly empty).
- [ ] Beam-gas and Touschek: listed as continuous backgrounds in `sources.json` (and shown in the sidebar) but no samples or rates yet.
- [ ] SR: one sample per component with an on/off toggle; none available yet (`sources.json` -> `sr.components` lists placeholder core/halo; same for injection).
- [ ] Z mixture weights (0.7 Zbb + 0.1 Zmumu) are hardcoded in the sheet; invisible Z contributes 0. Make them editable (single Z samples can already be selected). Other working points (WW, ZH, ttbar) have no samples, and their `f_bc`, bunches and `f_rev` are null in `working_points.json`.
- [ ] How totals treat physics vs background: Z has an include toggle (default on); a separate background-only total is not shown.
- [ ] Injection example peak scale (52 MHz-equivalent) is hardcoded in `app.js`; move it into the JSON when the scale factors arrive.

## Data quality (seed from the sheet, 2026-01-12)
- [ ] Seed occupancies are rounded to 0.01% (transcribed from displayed cells); tiny values (calo outer layers) can be off by ~2x. Replace with ROOT-derived numbers.
- [ ] HCalBarrel Zbb occupancy is 25-49% vs ~0.2% for IPC: suspicious.
- [ ] Some rows have p95 < mean or p95 < stdev (HCal, EMEC Zmumu), a sheet artefact.
- [ ] DCH_v2: only 10 of 112 layers.
- [ ] HCalThreePartsEndcap: no geometry, no occupancy.
- [ ] MuonTaggerBarrel/Endcap: only injection occupancies, no IPC or Z; endcap values were placeholders ("evaluate with reco muon positions").
- [ ] SiWrB: no LCC_V106 column; sheet todo "why hits in layers 2 and 3".
- [ ] DCH injection mean looks like single-event granularity (25%, 37.5%...).
- [ ] IDEA: vertex only; detector version in the sheet is not stated (`IDEA_o1_v03` assumed).
- [ ] Input path later: script from the ROOT files (`hits2highLevelEstimations.py` output) or manual JSON edits. `tools/seed_from_sheet.py` is a one-off.

## Technologies
- [ ] Parameters for 1 GbE PoE and silicon photonics are placeholders (`verified: false`); power per link is `null` for all.
- [ ] How to count power: links x W/link only, or include transceivers / switches?
- [ ] More presets wanted? (lpGBT, 10/25/100 GbE...)

## UI / hosting
- [ ] Sidebar (collapsible) with inputs, sample selector, technology, detector subset (default all), unit selector (Gbps, fibers, 1 GbE ports, SiPh links, power).
- [ ] Views: overview, per-subdetector layer drill-down, per-sample comparison, table with CSV export, permalink (parameters in the URL), provenance (sample, geometry, date).
- [ ] Link back from the readout-rates page to `www/fcc` (the ROOT-files page already links here, since 2026-11-12).

## To be continued (not started)
- [ ] Custom occupancy: think how to let a user enter their own occupancy numbers (per subdetector/layer, even bypassing ours), just to play. TBC. Open: where they are typed or pasted, how the page labels it so it is never mistaken for the stored numbers, and whether it can be exported.
