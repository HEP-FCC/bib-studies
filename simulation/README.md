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
# the same with the later one used to submit_bkg.py, otherwise MIGHT get ROOT or other mismatch errors (see 'Matching the software stack' below)
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
- in the compact file of the detector (e.g. `ALLEGRO_o1_v03.xml` or `IDEA_o1_v03.xml`), swap the shape-based beampipe for the CAD one. The lines are given by the name of the include, as their position changes between k4geo versions:
  - comment out the includes of `MDI_o1_v00/Beampipe_o4_v05.xml` **and** `MDI_o1_v00/BeamInstrumentation_o1_v01.xml` (the CAD collection brings its own beam instrumentation, keeping both gives a `CompSol` "already existing object" error)
  - un-comment the include of `MDI_o1_CADBased_v01/MDI_o1_CADBased_v01_collection.xml` and the line `<world material="Vacuum"/>` below it
  ```diff
  -  <include ref="../../../MDI/compact/MDI_o1_v00/Beampipe_o4_v05.xml" />
  -  <include ref="../../../MDI/compact/MDI_o1_v00/BeamInstrumentation_o1_v01.xml" />
  +  <!-- <include ref="../../../MDI/compact/MDI_o1_v00/Beampipe_o4_v05.xml" /> -->
  +  <!-- <include ref="../../../MDI/compact/MDI_o1_v00/BeamInstrumentation_o1_v01.xml" /> -->
  -  <!-- <include ref="../../../MDI/compact/MDI_o1_CADBased_v01/MDI_o1_CADBased_v01_collection.xml" />
  -  <world material="Vacuum"/> -->
  +  <include ref="../../../MDI/compact/MDI_o1_CADBased_v01/MDI_o1_CADBased_v01_collection.xml" />
  +  <world material="Vacuum"/>
  ```
- IFF running CLD (not re-checked with recent versions, see the CLD compact file of your release):
  - also remove the analytical compensating solenoid `<field>` blocks (the ones using `CompSol_min_z`), as the field is taken from a map in the above MDI import.
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

Example for processing through the ALLEGRO detector simulation an incoherent pair creation (IPC) file `pairs.pairs` (a text file with all particles' positions) generated with GuineaPig:

```bash
ddsim -N -1 \
 --inputFile /eos/experiment/fcc/users/a/aciarma/pairs/4IP_2024may29/Z/data1/pairs.pairs \
 --steeringFile mySteeringFile.py \
 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml \
 --outputFile sim_IPC_test.root
```

`simulation/submit_bkg.py` submits an HTCondor job for every background file, each simulating the N events of the file with `ddsim`.
It processes the edm4hep `.root` files (IPC from guinea-pig, SR halo) and the older `.pairs` / `.hepevt` text files (it replaces the former `submit_pairs.py`),
and leverages the metadata embedded in the `.root` files for the weight calculations.

Example usage, SR halo:
```sh
# lxplus!! , and NOT on EOS directory (condor submit will complain)
python3 $BIB_STUDIES/simulation/submit_bkg.py \
-i /eos/experiment/fcc/users/j/jaeyserm/bdsim/LCC_v2short_halo \
-t ALLEGRO_SR_halo \
-c $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml \
-s /absolute/path/to/my_BIB_steering_file.py \
-k /eos/user/<u>/<user>/k4geo \
-n 10
```

An SR file is a slice of a bunch, not a bunch crossing: see [SR.md](./SR.md) for the format, the weights, the number of frames per BX and how to get to a rate.

Example usage, IPC from guinea-pig in edm4hep format (Jan's files, 100 events = 100 BX per file):
```sh
# lxplus!! , and NOT on EOS directory (condor submit will complain)
python3 $BIB_STUDIES/simulation/submit_bkg.py \
-i /eos/experiment/fcc/users/j/jaeyserm/guineapig/FCCee_Z_LCC_V105_Z_2T_grids8_W100 \
-t IDEA_grids8_W100 \
-o <output/folder> \
-c $K4GEO/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml \
-s /absolute/path/to/my_BIB_steering_file.py \
-k /eos/user/<u>/<user>/k4geo \
--mcCollection Pairs \
-n 10
```

Key input arguments are:
- `-i` / `--input` input folder with the background files (`.pairs`, `.hepevt` or edm4hep `.root`)
- `-o` / `--output` output folder for the simulated files
- `-t` / `--tag` tag of the dataset, used to name the submission folder and the `.cmd` / `.sh` files
- `-c` / `--compactFile` detector geometry
- `-s` / `--steering_file` path to the steering file (required, there is no default). For BIB studies pass your own steering file with the detailed EM settings (see above). A relative path is converted to an absolute one, since relative paths break in the condor scratch directory

Other useful options are:
- `-n` / `--n_max_jobs` the maximum number of jobs to create (`-1` for all)
- `--nevents` the maximum number of events to simulate per job (default: all the events of each file), useful for quick tests since a full 100-event file takes hours
- `-k` / `--k4geo` path to a local k4geo build, set up in the job with `k4_local_repo` + `thisk4geo.sh`
- `-d` / `--dd4hep` path to a local DD4hep build, set up in the job before k4geo (e.g. for DD4hep changes not yet in the release)
- `--crossingAngleBoost` crossing angle boost passed to `ddsim` (default `0.015`, required for the SR samples that are given in the beamline frame)
- `--mcCollection` name of the `MCParticle` collection in the edm4hep `.root` inputs (`MCParticles` for SR, `Pairs` for IPC from guinea-pig; `Pairs0` is not simulated)
- `-j` / `--job_flavor` condor job flavour (default `tomorrow` => 24 h; a 100-event edm4hep file takes about 3-8 h with the CAD beampipe, so `espresso`/`longlunch` are too short and `workday` (8 h) leaves little margin)
- `--request_memory` condor memory request (default `12GB`; `ddsim` with the CAD beampipe reaches 7-8 GB resident and condor sees peaks up to ~10 GB)
- `--seed` random seed passed to `ddsim` (default `42`, so re-running a file reproduces the same result)
- `--stream_output` stream the job stdout to the `.out` file while the job runs (default off: condor returns it only at the end, so the log of a job killed by the wall-time limit is lost; the cost is continuous writes to the log directory, see the space note below)
- `--skip_space_check` do not check the space for the job logs (see below)

Notes:
- `ddsim` writes its output in the condor scratch directory and the file is moved to the output folder at the end of the job. The output is always written in EDM4hep format and keeps all the particles (`--part.keepAllParticles`).
- Condor accepts at most 20000 jobs per submission: for larger samples one `<tag>_part<N>.cmd` file is written per slice of 20000 files, and each of them has to be submitted.
- For `.root` inputs, the number of background events per file is read and a summary of the normalisation (events, files, bunch crossings covered) is printed at submission time.
- A `<tag>.sh` script to run the same jobs locally is also created.
- The job logs (`<ClusterId>.<ProcId>.{log,out,err}`, ~11 MB of `ddsim` stdout per 100-event job) are returned to the submission directory. If that directory runs out of space (e.g. an AFS home near its quota), the transfer fails after the job has finished and condor re-runs the whole job (up to 3 times), overwriting its output. The script therefore checks the free space (AFS quota with `fs lq`, or the filesystem) before writing the condor files and stops if the logs of all the jobs do not fit.
- The edm4hep inputs need a `ddsim` with the `--edm4hep.mcParticleCollectionName` option: it does not exist in the stable release 2025-05-29, but it is in the stable release 2026-02-01 and in the recent nightlies.
- Known limitation for the guinea-pig edm4hep files: the DD4hep edm4hep reader builds the primary vertex from the `endpoint` of the generator particles and not from their `vertex`. In these files the `endpoint` is 0, so all the primaries start at the origin and the vertex spread of the `Pairs` collection (z about +-33 mm, r < 10 mm) is not used. Momenta and times are kept.

#### Matching the software stack

The k4geo build, the local `ddsim` tests and the condor jobs must use the **same** key4hep release, otherwise ROOT or other version mismatches can appear.
1. Pick a release and source it, e.g. `source $BIB_STUDIES/setup.sh -n 2026-10-02` (a dated release does not change, so it is reproducible).
2. Build k4geo with this stack (recipe above) and `k4_local_repo` it. After a change of release, rebuild in a fresh build directory.
3. Run `submit_bkg.py` from the same shell and pass the build with `-k`. The script writes the release in the job scripts: a dated release as `setup.sh --spack -r <date>`, while for the LCG nightlies it writes `setup.sh --lcg`, which picks the newest LCG nightly **at the start of the job**. LCG nightlies are selected by weekday (Mon-Sat) and overwritten every week, so they are not suitable for production.
4. Check the header of one `<tag>/run_ddsim_*.sh` before submitting.

### Running the overlay
Before the digitization step, it is important to overlay the ddsim simulation output for which we are running our BIB studies to have a comprehensive picture of the detected occupancy for a physics
event produced at BX 0 for different subdetectors with diffreent readout windows.

The overlay is done with the `OverlayTiming` algorithm, configured in `overlay.py` and run with `k4run`.
It overlays the background onto the tracker hits of the DCH, silicon wrapper (`SiWrD`, `SiWrB`) and vertex (barrel, endcap) collections, 
with a different time window per sub-detector (`[-400, 400]` ns for the DCH, `[0, 25]` ns for the others).
Calorimeter hits are not overlaid.
The hits are saved in new `Overlay*` collections, together with the `OverlayMCParticles`.

Key input arguments of `overlay.py` are:
- `--inputFile` the physics simulation input
- `--outputFile` the overlaid EDM4hep output
- `--backgroundFilesPath` directory containing the background `.root` files (searched recursively)
- `--numEvents` the number of physics events to overlay (`-1` for all)

The bunch train (33 BX, 25 ns apart), the BX of the physics event and the number of background events per BX are set in `overlay.py`.

Note: `NumberBackground` is the number of background *events* per BX, which is not necessarily the number of BX that one background event represents.
Check the normalisation printed by `submit_bkg.py` (e.g. one SR halo file corresponds to about 1/2020 BX, while one IPC file is `bw_scale` BX) before interpreting the overlaid occupancy.

#### submit_overlay.py

`submit_overlay.py` generates the executable and the HTCondor submit file for a single `overlay.py` run
(the nightlies are used, since the `TimeWindows` property is not available in the stable releases),
and prints the command to launch it.

Example usage:
```sh
# lxplus!! , and NOT on EOS directory (condor submit will complain)
python3 $BIB_STUDIES/simulation/submit_overlay.py \
--inputFile <path/to/physics_sim.root> \
--backgroundFilesPath <path/to/background/files> \
--outputFile <path/to/overlaid.root>
```

Other useful options are:
- `--numEvents` the number of physics events to overlay (`-1` for all)
- `--jobFlavour` condor job flavour (default `tomorrow`)
- `--requestMemory` condor memory request (default `8GB`)
- `--tag` base name of the generated `.sh` / `.cmd` files, written in the current directory (default: from `--outputFile`)


### Running the digitization

To run digitization, please refer always to the FCC-config instructions for a given detector.

Currently available options:
- [ALLEGRO_o1_v02](https://github.com/HEP-FCC/FCC-config/tree/main/FCCee/FullSim/ALLEGRO/ALLEGRO_o1_v02#running-the-digitization-and-reconstruction)
- [ALLEGRO_o1_v03](https://github.com/HEP-FCC/FCC-config/tree/main/FCCee/FullSim/ALLEGRO/ALLEGRO_o1_v03#running-the-digitization-and-reconstruction)
- [IDEA_o1_v03](https://github.com/HEP-FCC/FCC-config/tree/main/FCCee/FullSim/IDEA/IDEA_o1_v03#running-the-digitization-and-reconstruction)


Note: this process run also some reconstruction by default, which can be computationally intensive and memory demanding.
For specific studies, some algorithms that aren't needed can be turned off (e.g. `--doTopoClustering=false`).

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


### `.pairs` / `.hepevt` inputs

The `.pairs` files (GuineaPig text files, one event each) and the `.hepevt` files are deprecated in favour of the edm4hep `.root` files, but `submit_bkg.py` still handles them:
it passes `-N -1` to `ddsim` and takes the file number from the file name (`<name>_XYZ.pairs`, `output_XYZ.hepevt`) or, for the older layout, from the `dataXYZ/pairs.pairs` folders.
The former `submit_pairs.py` script is superseded by `submit_bkg.py`.

----

[What is condor??](https://htcondor.readthedocs.io/en/25.0/users-manual/managing-a-job.html) and [how to submit my own jobs](https://batchdocs.web.cern.ch/local/quick.html)

```sh
condor_q #shows idle/running jobs
condor_q -nobatch # shows jobs, not in batch groups
condor_rm 10113537.1 # stops that job..
condor_q 10113537.1 -af JobStatus MemoryUsage RequestMemory HoldReason # 1 idle, 2 running, 5 held (+ why)
condor_history <user> -limit 5 # recently finished/removed jobs (check RemoveReason, e.g. wall time exceeded)
```

The job stdout/stderr are only copied back when the job ends. To look at a running job:
```sh
condor_tail 10113537.1 # tail of the stdout (e.g. which event ddsim is at)
condor_tail -stderr 10113537.1 # tail of the stderr
condor_tail -follow 10113537.1 # keep following the stdout
condor_ssh_to_job 10113537.1 # shell on the worker node, in the job scratch dir
```
The per-job log (`<ClusterId>.<ProcId>.log`, in the submit directory) shows the job events and the memory usage, updated every few minutes.


### All steps together

Full example, without explanations, just to see the full recipe.

```
source .../bib-studies/setup.sh -n 2026-10-02 #or the appropriate date when you compiled your local k4 repo
cd my/k4geo/at/eos
k4_local_repo
cd back/to/afs

submit_bkg.py \
-i /eos/experiment/fcc/users/j/jaeyserm/guineapig/FCCee_Z_LCC_V105_Z_2T_grids8_W100 \
-t IDEA_grids8_W100 \
-n 4000 \
-c $K4GEO/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml \
-o /eos/user/a/aikoulou/fcc_workdir/samples/bib/ipc/jaeyserm_Z_4IP_grids8_W100 \
-s /absolute/path/to/my_BIB_steering_file.py \
--mcCollection Pairs \
--k4geo /eos/user/a/aikoulou/fcc_workdir/k4geo/

condor_submit IDEA_grids8_W100.cmd
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

Note: a single file cannot exceed 50 GB on `/eos/user` (`sys.forced.maxsize`), and both `merge_edm4hep_files.py` and `podio-merge-files` fail there when the merged file gets that large.
Use `--files-per-output` so that every merged file stays below it (e.g. 50 ALLEGRO files of 0.84 GB each), or merge on a local disk.

## Calorimeter calibration

TODO: add docu

## Tracking studies/performance

Tracking validation (finder/fitter efficiency, resolutions, pulls vs. MC truth and vs. perfect tracks) comes from the `TrackingValidation` algorithm of [k4DetectorPerformance](https://github.com/key4hep/k4DetectorPerformance/tree/main/TrackingPerformance#trackingvalidation), which is not in the key4hep stack yet, so build it locally.
Checked with ALLEGRO_o1_v03 on the key4hep nightly 2026-10-01:

```sh
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh --spack -r 2026-10-01
git clone https://github.com/key4hep/k4DetectorPerformance.git
cd k4DetectorPerformance
k4_local_repo
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install -j 8
```

Then, in a shell with the same stack and `k4_local_repo` run from `k4DetectorPerformance` (not from a k4geo with BIB edits: use the stack's k4geo), run the [FCC-config](https://github.com/HEP-FCC/FCC-config/tree/main/FCCee/FullSim/ALLEGRO/ALLEGRO_o1_v03) sim+digi+reco script with `--runTrkValidation` (no need to edit the script):

```sh
bash $FCCCONFIG/FullSim/ALLEGRO/ALLEGRO_o1_v03/ctest_sim_digi_reco.sh --particle mu- --energy "10*GeV" --nEvents 20 --seed 1 --outputFile trkval --runTrkValidation
```

- run it in an empty directory, the script downloads the calorimeter input files into the current directory
- output: `ALLEGRO_trkval_sim.root`, `ALLEGRO_trkval_digi.root` and `ALLEGRO_sim_digi_reco_trkval.root` (the validation output: finder/perfect/fitter TTrees, resolution and efficiency vs. p graphs, pull histograms). The name does not follow `--outputFile`
- needs a machine with `libGLU.so.1` (e.g. lxplus), otherwise `libTracking.so` fails to load and the digitization algorithms are reported as unknown
- on lxplus, processes started in the background die when the ssh session closes: keep the session open or use condor

