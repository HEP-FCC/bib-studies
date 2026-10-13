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

Example for processing through the ALLEGRO detector simulation an incoherent pair creation (IPC) file `pairs.pairs` (a text file with all particles' positions) generated with GuineaPig:

```bash
ddsim -N -1 \
 --inputFile /eos/experiment/fcc/users/a/aciarma/pairs/4IP_2024may29/Z/data1/pairs.pairs \
 --steeringFile mySteeringFile.py \
 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml \
 --outputFile sim_IPC_test.root
```

simulation/submit_bkg.py serves as a generalization of the previous submit_pairs.py, which submits an HTCondor job for every background file submitted, each containing N events. 
The generalized script allows the processing of .root IPC and SR files in the new format, and leverages embedded metadata for weight calculations etc.

Example usage:
```sh
# lxplus!! , and NOT on EOS directory (condor submit will complain)
python3 $BIB_STUDIES/simulation/submit_bkg.py \
-i /eos/experiment/fcc/users/j/jaeyserm/bdsim/LCC_v2short_halo \
-t ALLEGRO_SR_halo \
-c $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml \
-k /eos/user/<u>/<user>/k4geo \
-n 10
```

Key input arguments are:
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

