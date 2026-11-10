# Simulation

Preparing and submitting simulations of background samples.

## List of available BIB generated samples

- Join egroup `fcc-ee-MDI`
- Find supported samples in the Machine-Detector Interface (MDI) EOS space:
```sh
/eos/project/f/fcc-ee-mdi/BIB/

#see readme in the folder for more info
```
also accessible through [this CERN Box link](https://cernbox.cern.ch/files/spaces/eos/project/f/fcc-ee-mdi/BIB). See readme in the folder for sample definitions and contact persons.

## Running detector simulation

### Prepare your setup

The configuration used to run the simulation for BIB studies is slightly different than the one used for physics event processing. Mainly because of the following points:
- We need a detailed modeling of the MDI elements --> we use the (slow and imperfect) CAD based beampipe
- Due to technical difficulties, there is air inside the CAD beampipe --> we use a temporary workaround setting the world volume as vacuum while waiting for a better solution
- To properly model the effect of BIB, a detailed treatment of EM processes has to be used (e.g. we enable fluorescence)

Here is a **full recipe to run the simulation** in the appropriate conditions for BIB studies:
```bash
# connect to an Alma9 machine with cvmfs mounted
source /cvmfs/sw.hsf.org/key4hep/setup.sh
#important note: this fccsetup version should be
# the same with the later one used to submit_pairs, otherwise MIGHT get ROOT or other mismatch errors
git clone https://github.com/key4hep/k4geo
cd k4geo
mkdir build install
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install -D INSTALL_BEAMPIPE_STL_FILES=ON
make install -j 8
cd ..
k4_local_repo
```

Now let's switch to the CAD beampipe, and set vacuum everywhere (ALLEGRO is taken as an example but it works the same way for other detectors):
- comment out [these lines](https://github.com/key4hep/k4geo/blob/main/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml#L34-L35)
- and un-comment [these lines](https://github.com/key4hep/k4geo/blob/main/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml#L40-L41)
- IFF running CLD:
  - also need to comment out [these lines](https://github.com/key4hep/k4geo/blob/main/FCCee/CLD/compact/CLD_o2_v08/CLD_o2_v08.xml#L401-L415) to remove the analytical compensating solenoid field which is taken from a map in the above MDI import.
  - You also need to add some material to the detector list of materials (see e.g. [here](https://github.com/key4hep/k4geo/pull/534/commits/2a2ea2591db1473d294af5c432f99aac74b8dea7#diff-f42d88422d9f50cb0863b6f08f2640a9e5cbcb9ac2ae01145642105d9fe9387d)).
- **enable detailed EM treatment in Geant4** by applying the following changes to the `ddsim` steering file
  - (if you do not already use a `ddsim` steering file, you can create the default one with `ddsim --dumpSteeringFile > mySteeringFile.py`):
  - Change the physics list to `SIM.physics.list = "FTFP_BERT_EMZ"`
  - Change the range cut: `SIM.physics.rangecut = 0.05*mm`
  - Remove the energy threshold for tracker hits: `SIM.filter.tracker = "edep0"`
  - At the bottom of the file, change the Geant4 UI configure commands to:
    ```
    SIM.ui.commandsConfigure = [
    "/cuts/setLowEdge 50 eV",
    "/process/em/lowestElectronEnergy 1 eV",
    "/process/em/auger true" ,
    "/process/em/deexcitationIgnoreCut true"]
    ```
- For some BIB (e.g. IPC), the **boost due to the crossing angle has to be applied**:
  - At the beginning of the file, use: `SIM.crossingAngleBoost = 0.015`
  - (boost depends on BIB generation => contact responsible/creator if in doubt)

Regarding steering files:
- centrally maintained in [`FCC-Config`](https://github.com/HEP-FCC/FCC-config), which is included in the `key4hep` software stack and accessible with the environment variable `$FCCCONFIG`
- Example: [IDEA_o1_v03](https://github.com/HEP-FCC/FCC-config/blob/main/FCCee/FullSim/IDEA/IDEA_o1_v03/SteeringFile_IDEA_o1_v03.py) can be used with:
```bash
ddsim --steeringFile $FCCCONFIG/share/FCC-config/FullSim/IDEA/IDEA_o1_v03/SteeringFile_IDEA_o1_v03.py ...
# or for the nightlies
ddsim --steeringFile $FCCCONFIG/FullSim/IDEA/IDEA_o1_v03/SteeringFile_IDEA_o1_v03.py ...
```
- Note: For CLD, the centrally maintained steering file lives [here](https://github.com/key4hep/CLDConfig/blob/main/CLDConfig/cld_steer.py) and can be accessed through `$CLDCONFIG`.

### Run the simulation

Example recipe for the synchrotron radiation (SR) halo background (edm4hep `.root` files, 1 event per file, 1 file ~ 1/2020 BX) in the ALLEGRO detector, in three steps:
simulate with `submit_bkg.py`, trace the tracker hits back to the MDI with `hit_traceback.py`, and compute the occupancy with `occupancy.py`.

```bash
# 1. Simulation (lxplus!!, and NOT on an EOS directory: condor submit will complain)
python3 $BIB_STUDIES/simulation/submit_bkg.py \
 -i /eos/experiment/fcc/users/j/jaeyserm/bdsim/LCC_v2short_halo \
 -t ALLEGRO_SR_halo \
 -c $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml \
 -s mySteeringFile.py \
 -k /eos/user/<u>/<user>/k4geo \
 -o <path/to/sim/output> \
 -n 10
condor_submit <tag>.cmd   # the command is printed at the end of the submission

# 2. Traceback of the VTX, DCH and SiWr hits to the origin in the MDI
# (Key4hep nightly that wrote the files; in a condor/bash script, run `set --` before `source`)
python3 $BIB_STUDIES/traceback/hit_traceback.py <path/to/sim/output> <path/to/traceback/output> --n-files 100 --merge

# 3. Occupancy of the DCH: the whole input is one readout window (31 BX x 2020 files/BX)
python3 $BIB_STUDIES/python/occupancy.py -i <path/to/sim/output> -n 62620 -s DCH --sample SR_halo \
 -d $BIB_STUDIES/detectors_dicts/ALLEGRO_o1_v03_DetectorDimensions.json --bx_per_event 1/2020
```

Notes on the recipe:
- Step 1: `mySteeringFile.py` is the BIB steering file prepared above (detailed EM treatment, `edep0` tracker filter). The crossing angle boost is applied by default (`--crossingAngleBoost 0.015`), as the SR samples are given in the beamline frame. Use `-n -1` to simulate all the files.
- Step 2: the traceback finds where in the MDI the background that hits the detector comes from. For every VTX, DCH or SiWr hit it follows the parent chain of the hit particle up to the primary SR photon, and records the earlier scatter points of the photon ancestors (the "kinks", where the photon was deflected or produced a shower). Every vertex is labelled by region (injection, MDI, the sub-detector itself or other), and a vertex-map plot of the MDI scatter points can then be drawn. The input has to be the ddsim output with all particles kept (the default of `submit_bkg.py`). Per run it writes `trees_<tag>.jsonl`, `collections_<tag>.json` and the skimmed EDM4hep file (MCParticles + traced tracker hits only), and skips truncated files (`skipped_<tag>.txt`). Process at most ~100 files per job (`--first-file` / `--n-files`), since the podio reader slows down with many files. The vertex maps are drawn from the trees with `traceback/plot_merged_maps.py`. See `traceback/NOTES.md`.
- Step 3: use `-s SiWrB`/`SiWrD` or the vertex detector names of the dictionary for the other sub-detectors, and `--per_event` instead of the window logic for overlaid/digitized files (see [Computing the occupancy](#computing-the-occupancy)). The error is the leave-one-event-out jackknife.

`submit_bkg.py` generalizes the previous `submit_pairs.py`: it submits an HTCondor job for every background file, each containing N events.
It processes the `.root` IPC and SR files in the new edm4hep format (and the `.pairs`/`.hepevt` IPC files), and leverages the embedded metadata for the weight calculations etc.

Key input arguments of `submit_bkg.py` are:
- `-i` / `--input` input folder with the background files (`.pairs`, `.hepevt` or edm4hep `.root`)
- `-o` / `--output` output folder for the simulated files
- `-t` / `--tag` tag of the dataset, used to name the submission folder and the `.cmd` / `.sh` files
- `-c` / `--compactFile` detector geometry
- `-s` / `--steering_file` path to the steering file (if not given, a default is looked up from `$FCCCONFIG`)

Other useful options are:
- `-n` / `--n_max_jobs` the maximum number of jobs to create (`-1` for all)
- `-k` / `--k4geo` path to a local k4geo build, set up in the job with `k4_local_repo` + `thisk4geo.sh`
- `-d` / `--dd4hep` path to a local DD4hep build, set up in the job before k4geo (e.g. for DD4hep changes not yet in the release)
- `--crossingAngleBoost` crossing angle boost passed to `ddsim` (default `0.015`, required for the SR samples that are given in the beamline frame)
- `--mcCollection` name of the `MCParticle` collection in the edm4hep `.root` inputs (`MCParticles` for SR, `Pairs` for IPC from guinea-pig)
- `-j` / `--job_flavor` condor job flavour (default `espresso`)
- `--seed` random seed passed to `ddsim` (default `42`, so re-running a file reproduces the same result)

Notes:
- `ddsim` writes its output in the condor scratch directory and the file is moved to the output folder at the end of the job. The output is always written in EDM4hep format and keeps all the particles (`--part.keepAllParticles`).
- Condor accepts at most 20000 jobs per submission: for larger samples one `<tag>_part<N>.cmd` file is written per slice of 20000 files, and each of them has to be submitted.
- For `.root` inputs, the number of background events per file is read and a summary of the normalisation (events, files, bunch crossings covered) is printed at submission time.
- A `<tag>.sh` script to run the same jobs locally is also created.


### Computing the occupancy

`python/occupancy.py` computes the occupancy (fired cells / total cells) of a sub-detector, in total and per layer,
with an uncertainty that takes the clustering of the hits into account (one particle fires several cells).
It can be run on the output of any of the previous steps (simulation, overlay, digitization),
and the cells per layer are taken from the detector dictionary.

Two modes, depending on what one input event represents:
- **default**: the whole input is one readout window (e.g. many raw single-BX files, which add up to a window).
  Each cell counts once, however many events fire it. The error is the leave-one-event-out jackknife, which needs at least ~10 events.
- `--per_event`: every event is one readout window (e.g. the overlaid and digitized files, where an event already holds a whole bunch train).
  The occupancy of each event is computed, and the result is the mean over the events with its standard error.

Example usage:
```sh
# inside the Key4hep environment (e.g. after source $BIB_STUDIES/setup.sh)
# raw IPC simulation: all the files are one window of 31 BX
python3 $BIB_STUDIES/python/occupancy.py -i <path/to/sim/files> -n 31 -s DCH --sample ipc \
-d $BIB_STUDIES/detectors_dicts/ALLEGRO_o1_v03_DetectorDimensions.json --bx_per_event 1

# overlaid + digitized windows: one event = one window
python3 $BIB_STUDIES/python/occupancy.py -i <path/to/overlaid_digi.root> -n 1 -e -1 -s DCH --per_event --digi --e_cut \
-d $BIB_STUDIES/detectors_dicts/ALLEGRO_o1_v03_DetectorDimensions.json \
-a $BIB_STUDIES/detectors_dicts/ALLEGRO_o1_v03_assumptions.json
```

Key input arguments are:
- `-i` / `--infilePath` input file or folder (the files are sorted, so the same ones are always taken)
- `-s` / `--subDetector` sub-detector, as named in the detector dictionary
- `-d` / `--detDictFile` detector dictionary, with the number of cells per layer (default ALLEGRO_o1_v03)
- `-n` / `--numberOfFiles` and `-e` / `--numberOfEvents` number of files and of events per file to use (`-1` for all)
- `--per_event` each event is one readout window (see above)

Other useful options are:
- `--bx_per_event` BX represented by one input event (e.g. `1` for raw IPC, `1/2020` for the SR halo), used only to warn if the window is not `--window_bx` BX long (default 31)
- `--digi`, `--e_cut` and `-a` / `--assumptions` use the digitized collection and/or apply the energy threshold, both defined in the assumptions dictionary
- `--collection` hits collection to read, e.g. `OverlayDCHCollection` for overlaid but non-digitized files
- `-p` / `--draw_hists` also draw the occupancy per layer (with the total as a dashed line)
- `--sample` name used in the plot and output file names

Notes:
- The total and per-layer occupancy are printed, and saved in a `<sample>_<N>evt[_perEvent]_<subDetector>_<input>.root` file in the current directory.
- The window length is not extrapolated: if the input covers fewer (more) BX than `--window_bx`, the occupancy is under(over)estimated, and this bias is not included in the error.
- Check the BX covered by one event before interpreting the number (printed by `submit_bkg.py` at submission: SR halo, one file is about 1/2020 BX; IPC, check the `bw_scale` metadata and the number of events per file).
- Pixel detectors: the cellID identifies a pixel but the dictionary counts sensors, so the occupancy is fired pixels / sensors (a warning is printed).
- The `drawhits.py` occupancy is instead computed per event (see the [plotting docs](../plotting/README.md)).

## Production of IPC backgrounds

The production of Incoherent Pair Creation (IPC) background simulation takes two steps.

### set_vertex_000.py

**Not needed** for files after end of 2025 (Jan fixed it ;))

<details>
<summary>Click to expand</summary>


Reset the position of particles to (0,0,0) in `.pairs` files 
created by GuineaPig. This is required as the event generator doesn't
include any B-field. Therefore, the positions are inexact, especially if
particles travel for radiuses larger than the beam pipe.
E.g. see slide 5 in
[Brieuc slides](https://indico.cern.ch/event/1559862/contributions/6608302/attachments/3107855/5508385/20250721_StatusOfBkgStudiesWrtSoftware.pdf).
This adjustment is also inexact, but more realistic.
In the future it might not be needed anymore.

Example usage:
```
set_vertex_000.py -i <regex/to/input/dirs> -o <path/to/output/dir>
```
Current default input points to A. Ciarma's IPC samples:
```
/eos/experiment/fcc/users/a/aciarma/pairs/4IP_2024may29/Z/data*
```
Note that not all the folders contain a `.pair` file, 
but only a `.dat` version of it. In that case, the `--do_dat` flag might be needed.


</details>


### submit_pairs.py

How to set up the condor (or local) submission of simulation jobs of 
IPC background files:

At the moment, `.pairs` files contain a single event
and `ddsim` can process only one of them at the time.
The  `submit_pairs.py` script generates at list of bash scripts
to automatize the submission of many single event jobs.
After the preparation is done, the command to launch the jobs
(condor or locally) is printed in the terminal.

Example usage command:
```sh
#important note: this fccsetup version should be the same with the earlier one used to compile k4geo, otherwise MIGHT get ROOT or other mismatch errors
# lxplus!! , and NOT on EOS directory (condor submit will complain)
# make sure the setup cmd matches the k4_local_repo you use (eg fccsetupnightly -r 2026-03-23)
submit_pairs.py --tag IDEA_my_test --compactFile $K4GEO/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml -n 10
```

All the available geometries are stored in the
[`k4geo`](https://github.com/key4hep/k4geo/tree/main)
repository.

If a custom variation of the standard geometry 
tha requires recompiling k4geo is needed,
specify the path to the local build with the `--k4geo` flag.

Example running on Jan's recent files 
```sh
# lxplus!! , and NOT on EOS directory (condor submit will complain)
# make sure the setup cmd matches the k4_local_repo you use (eg fccsetupnightly -r 2026-03-23)
submit_pairs.py \
-i /eos/experiment/fcc/users/j/jaeyserm/guineapig/guineapig_samples_CERN_oct25/FCCee_Z_4IP_FSR_FCCee_Z256_2T_grids8 \
-t ALLEGRO_FSR_FCCee_Z256_2T_grids8 \
-n 10 \
-c $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml \
-o /eos/home-s/sfranche/FCC/samples/bib/ipc/jaeyserm_Z_4IP_FSR_FCCee_Z256_2T_grids8 \
-s /eos/home-a/aikoulou/fcc_workdir/bib-studies/plotting/mySteeringFile.py \
--k4geo /eos/user/a/aikoulou/fcc_workdir/k4geo/
```

All the available options can be seen using the `-h` flag.
Input path `--input <your/path>` is expected to contain files with the naming format: `your/path/*_XYZ.pairs`,
where `XYZ` is an event number.

----

[What is condor??](https://htcondor.readthedocs.io/en/25.0/users-manual/managing-a-job.html) and [how to submit my own jobs](https://batchdocs.web.cern.ch/local/quick.html)

```sh
condor_q #shows idle/running jobs
condor_q -nobatch # shows jobs, not in batch groups
condor_rm 10113537.1 # stops that job..
```


### All steps together

Full example, without explanations, just to see the full recipe.

```
source .../bib-studies/setup.sh -r 2026-03-23 #or the appropriate date when you compiled your local k4 repo
cd my/k4geo/at/eos
k4_local_repo
cd back/to/afs

submit_pairs.py \
-i /eos/project/f/fcc-ee-mdi/BIB/GHC/V25.3-4/IPC/Z/ \
-t IDEA_SIM.enableDetailedShowerMode_GHC_V25.3-4_Z \
-n 4000 \
-c $K4GEO/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml \
-o /eos/user/a/aikoulou/fcc_workdir/samples/bib/ipc/jaeyserm_Z_4IP_GHC_V25.3-4_Z \
-s $FCCCONFIG/FullSim/IDEA/IDEA_o1_v03/SteeringFile_IDEA_o1_v03.py \
--k4geo /eos/user/a/aikoulou/fcc_workdir/k4geo/

condor_submit IDEA_SIM.enableDetailedShowerMode_GHC_V25.3-4_Z.cmd
```

### Merge sim output files

For IPC, as an example, multiple files (1 per event) are produced, so it's good to merge them in one, for easier access and loading later.

```sh
#example
podio-merge-files --output-file test.root /eos/home-s/sfranche/FCC/samples/bib/ipc/jaeyserm_Z_4IP_29may24_FCCee_Z256_2
T_grids8/ALLEGRO_29may24_FCCee_Z256_2T_grids8/ALLEGRO_o1_v03_99*

#note1:
- if mergine many files, might not get output for a few minutes in the very beginning

#note2:
podio-merge might run extremely slow sometimes
- check first if the running PC is not overloaded, and that EOS is not generally slow
- try to move the files off EOS (eg to local disk somewhere), and retry there
- 
```

#### merge_edm4hep_files.py

For the edm4hep samples (e.g. the SR halo, 1 event per file), `python/merge_edm4hep_files.py` merges many single-event files into chunked output files,
with one output event per input event.
Collections can be dropped at the ROOT level (without reading them), which reduces the size of the merged files.

Key input arguments are:
- `--input-dir` directory with the input files
- `--output-dir` output directory (default: `<input-dir>/../merged`)
- `--drop <regex>` drop the event collections matching the regex (can be repeated)
- `--no-mcparticles` also drop `MCParticles` (about 98% of the file size), the hit to `MCParticle` links are then not available (see the [plotting docs](../plotting/README.md))

Other useful options are:
- `--pattern` glob pattern for the input files (default `*_edm4hep_*.root`)
- `--files-per-output` number of input files per output file (`0` for all in one file)
- `--prefix` output file name prefix
- `--max-files` use only the first N input files, for testing

Example usage:
```sh
# inside the Key4hep nightly environment
python3 $BIB_STUDIES/python/merge_edm4hep_files.py --input-dir <path/to/files> --files-per-output 100 --no-mcparticles
```

## Calorimeter calibration

TODO: add docu

## Tracking studies/performance

To be exercised, instructions not checked yet: 

- clone the k4DetectorPerformance package (add it to your path with k4_local_repo (not yet in nightlies)
- switch that to true  https://github.com/HEP-FCC/FCC-config/blob/main/FCCee/FullSim/ALLEGRO/ALLEGRO_o1_v03/ctest_sim_digi_reco.sh#L11
- good to go 

See also https://github.com/key4hep/k4DetectorPerformance/tree/main/TrackingPerformance#trackingvalidation

