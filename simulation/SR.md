# Synchrotron radiation (SR) background: inputs, weights, simulation and rates

Draft (2026-11-12). Items marked **Open** are not understood or not verified yet.
Contacts for the generation: G. Nigrelli, J. Eysermans (see also `BIB/README.md` on `/eos/project/f/fcc-ee-mdi/`).

SR is not simulated like the incoherent pair creation (IPC): the photons come from a standalone BDSIM model of the final focus
and they are only a small random slice of the beam, so a file is **not** a bunch crossing.
Read section 3 before using any SR sample.

## 1. The samples

| Production | Format | Primaries per file | Files | Location |
|---|---|---|---|---|
| LCC V106.2 `SR/core` | HEPEVT `.hepevt` | 200000 | 65184 | `/eos/project/f/fcc-ee-mdi/BIB/LCC/V106.2/SR/core` |
| LCC V106.2 `SR/halo` | HEPEVT `.hepevt` | 200000 | 21641 (21644 listed) | `/eos/project/f/fcc-ee-mdi/BIB/LCC/V106.2/SR/halo` |
| BDSIM `LCC_v2short_halo` (Jan) | edm4hep `.root` | 1000000 | 22708 | `/eos/experiment/fcc/users/j/jaeyserm/bdsim/LCC_v2short_halo` |
| BDSIM `LCC_v3short_halo_bdsim_def_G4SR` (Jan) | edm4hep `.root` | 1000000 | 34622 | `/eos/experiment/fcc/users/j/jaeyserm/bdsim/LCC_v3short_halo_bdsim_def_G4SR` |

A "primary" is a macro-particle (an electron) tracked by BDSIM. The files only keep the photons it produced that pass the cuts (section 2).
Already simulated samples (ALLEGRO_o1_v03 with the CAD beampipe, key4hep 2026-04-08, from the V106.2 `.hepevt` files):
`/eos/experiment/fcc/ee/simulation/key4hep_2026_04_08/91GeV/<detector>_mdi_v2_CAD/sim/{SR_core,SR_halo}_LCC_V106p2` (see the FCC SW MDI page).

**Open**: how the V106.2 and the `short` productions differ (lattice version, BDSIM physics, which one is the reference).

## 2. File format, numbers and assumptions

Common to all the samples (V106.2 `SR/ReadMe.md` and the metadata of the edm4hep files):
- Photons only, one beam ("one beam" in the metadata), forward going, recorded on a sampler at z = -2400 mm (`samplerName` `QD0AL`, s = 6.8 m).
- Cuts: E > 2 keV, r < 18 mm at the sampler (z = -2.4 m), and r > 9 mm at +8.4 m from the IP (`SR/ReadMe.md`; `cutRadiusMinExtrapolated_m` and `cutExtrapolationLength_m` in the metadata).
- Bunch intensity 2.02e11. Halo: 1 % of the bunch charge, weighted exponential tail for a 5 minute lifetime. Core: 99 %, matched 6-D Gaussian (eps_x = 0.7 nm, eps_y = 2.6 pm, sigma_delta = 0.1 %).
- The crossing-angle Lorentz transformation (+15 mrad) must be applied (`ddsim --crossingAngleBoost`, default 0.015 in `submit_bkg.py`).

**HEPEVT** (`output_<id>.hepevt`):

- one "event" per file
- First line = number of particles
- then one line per photon with 15 fields: `status pdg mother1 mother2 daughter1 daughter2 px py pz E m x y z t` (GeV, mm, mm/c).
- Measured: halo ~56000 photons per file (44000-73000, median E 20 keV), core ~2470 photons per file (2330-2630, median E 3 keV), all at z = -2400 mm.

**edm4hep** (`output_<id>_edm4hep.root`): 

- one `events` frame per file holding all the photons as the `MCParticles` collection (about 2e5-3e5 per file)
- a `metadata` frame with
    - `nPrimariesTotal`
    - `nEvents`
    - `nPhotons`
    - `bunchIntensity`
    - `chargeFraction`
    - `runConfig`
    - the cuts (`cutEnergyMin_GeV`, `cutRadiusMax_m`, ...)
    - `zSampler_mm` and 
    - `weightFormula = "weight = chargeFraction * bunchIntensity / sum(nPrimaries over all events used); counts per bunch crossing, one beam"`.

Photons per primary: core ~0.012 (V106.2), halo ~0.28 (V106.2), 0.27 (`v2short`), 0.21 (`v3short`) (one file each for the last two).

## 3. What a file, a sample and a set correspond to (the "event" definition)

For IPC one `.root` event (or one `.pairs` file) is one bunch crossing (BX). **For SR it is not.**

- A file (HEPEVT) or an `events` frame (edm4hep, `nEvents = 1`) is the photon output of `nPrimariesTotal` macro-particles (200000 or 1000000).
  `ddsim -N` counts these frames. It is not a BX, not a turn and not a fraction of a turn: it is a random slice of the bunch.
- A real bunch holds `chargeFraction * bunchIntensity` electrons in the component: 2.0e11 for the core, 2.02e9 for the halo.
  So the number of frames that represent one BX is

      frames per BX = chargeFraction * bunchIntensity / primaries_per_frame

  | Sample | primaries / frame | frames per BX |
  |---|---|---|
  | V106.2 core | 200000 | ~1.0e6 |
  | V106.2 halo | 200000 | ~1.0e4 |
  | `v2short` and `v3short` halo | 1000000 | ~2020 |

  Overlaying "one frame per BX", as for IPC, would be wrong by these factors. `submit_bkg.py` prints this ("BX covered", "frames per BX")
  for the edm4hep inputs (`edm4hep_helper.print_normalisation`).
- The whole set therefore covers `N_frames / frames_per_BX` BX when used with unit weights:
  V106.2 core 0.065 BX (65184 files), V106.2 halo ~2.1 BX (21641), `v2short` ~11.2 BX (22708), `v3short` ~17.1 BX (34622).
  A set with less than 1 BX (the core) needs a weight above 1, a set with more than 1 BX a weight below 1.
- Only the full set, or a subset with the weight rescaled (section 4), gives the correct counts. A few files are a statistical fraction of a BX, not "a BX".
- **Open**: SR files are slices of a BX of one beam; whether the sample must be doubled for the second beam / the +z side is not stated (the metadata says "one beam").
  The relation to orbits or turns is not documented either (the weight is per BX, not per turn).

## 4. The weights

Each simulated primary stands for `w` real electrons of that component, and every photon it produced carries that weight:

    w = chargeFraction * bunchIntensity / (total primaries of the files used)

Checked with the numbers of section 1: V106.2 core 0.99 * 2.02e11 / (65184 * 2e5) = 15.34, V106.2 halo 0.01 * 2.02e11 / (21641 * 2e5) = 0.467 (the values in `SR/ReadMe.md`),
`v2short` 0.089, `v3short` 0.0583. Multiplying any photon (or hit) count of the sample by `w` gives the count per BX per beam.

- If only n of the ntot files are used: `w = w_full * ntot / n` (as in `SR/ReadMe.md`). For the edm4hep files use `python/edm4hep_helper.py: sr_weight(path, nevents_read)`,
  it takes the constants from the metadata of any file of the sample and the number of frames that were really read.
- The weight belongs to the whole file, the photons are not tagged with the primary that produced them.
- **Open**: confirm with the authors that "counts per BX" is meant per beam and per IP, and whether the core and the halo are simply added.

## 5. How to simulate

With `submit_bkg.py` (see `simulation/README.md` for the environment, the CAD beampipe recipe and the options), one HTCondor job per file:

```sh
# lxplus, and NOT from an EOS directory
python3 $BIB_STUDIES/simulation/submit_bkg.py \
-i /eos/experiment/fcc/users/j/jaeyserm/bdsim/LCC_v2short_halo \
-t ALLEGRO_SR_halo \
-c $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml \
-s /absolute/path/to/my_BIB_steering_file.py \
--mcCollection MCParticles \
-n 10
```

- `.hepevt` inputs (V106.2) are handled the same way (`-i .../SR/core`), `ddsim -N -1` reads the single event.
- Run a few files first (`-n`, or `--nevents`), the number of hits per file is small for the core and large for the halo.
- Local test, one file: the `ddsim` command that `submit_bkg.py` writes in `<tag>/run_ddsim_<id>.sh` can be run directly (same release as the k4geo build).
- The central samples were produced with `ddsim --steeringFile steering.py --compactFile .../ALLEGRO_o1_v03_CAD.xml --inputFiles output_<id>.hepevt --numberOfEvents -1`
  (from the `runs` metadata of the output file): no `--crossingAngleBoost` appears on that command line.
  **Open**: is the boost set in their `steering.py`, or is it missing? Compare before mixing our samples with the central ones.
- Overlay: `OverlayTiming` consumes one background frame per BX. With 1e4 to 1e6 frames per BX (section 3) the SR cannot be overlaid frame by frame, **Open**: how to do it (weights, or a reduced number of frames).

## 6. From hits to a rate

`drawhits.py` normalises per event (`fill_weight = 1 / n_events`): for SR that is **per frame, not per BX**, and the SR weight is not applied anywhere in the plotting scripts yet.
With the definitions above, the hits of one component per BX and per beam are

    hits per BX = (hits per frame) * (frames per BX)        [= w * total hits of the sample]

and the rate follows with the bunch crossing rate `f_BX` (`IPC_BX_RATE_MHZ` in `python/constants.py`, 52 MHz for the IPC):

    rate = f_BX * hits per BX

Sources are added: IPC (1 event = 1 BX, already per BX), SR core, SR halo. Pixels fired = hits * cluster size, bandwidth as in `hits2highLevelEstimations.py`.
- **Open**: apply `frames per BX` in `drawhits.py` (a `--frames_per_bx` option, or read it from the metadata) so that the per-BX numbers come out of the plots and of `hits2highLevelEstimations.py`.
- **Open**: second beam / +z side (section 3), and the safety factors for SR (statistics are poor for the core: 0.065 BX).
- For the later web presentation the numbers needed per source are: frames per BX, `f_BX`, hits per frame per sub-detector, cluster size.
