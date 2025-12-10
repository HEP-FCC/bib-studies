# Simulation

Preparing and submitting simulations of background samples.

## List of available BIB generated samples

All the supported samples are stored in the Machine-Detector Interface (MDI) EOS space:
```
/eos/project/f/fcc-ee-mdi/BIB/
```

For information on the organisation of that directory, please refer to the README.md:
```
/eos/project/f/fcc-ee-mdi/BIB/README.md
```
also accessible through [this CERN Box link](https://cernbox.cern.ch/files/spaces/eos/project/f/fcc-ee-mdi/BIB/README.md).
That file contains the list of reference samples to be used for the various beam background sources
and the respective contact persons.


  
## Running detector simulation
### Prepare your setup
The configuration used to run the simulation for BIB studies is slightly different than the one used for phyiscs event processing. Mainly because of the following points:
- We need a detailed modeling of the MDI elements --> we use the (slow and imperfect) CAD based beampipe
- Due to technical difficulties, their is air inside the CAD beampipe --> we use a temporary workaround setting the world volume as vacuum while waiting for a better solution
- To properly model the effect of BIB, a detailed treatment of EM processes has to be used (e.g. we enable fluorescence)

Here is a **full recipe to run the simulation** in the appropriate conditions for BIB studies:
```bash
# connect to an Alma9 machine with cvmfs mounted
source /cvmfs/sw.hsf.org/key4hep/setup.sh
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

NB: for CLD, you further need to comment out [these lines](https://github.com/key4hep/k4geo/blob/main/FCCee/CLD/compact/CLD_o2_v08/CLD_o2_v08.xml#L401-L415) to remove the analytical compensating solenoid field which is taken from a map in the above MDI import. You also need to add some material to the detector list of materials (see e.g. [here](https://github.com/key4hep/k4geo/pull/534/commits/2a2ea2591db1473d294af5c432f99aac74b8dea7#diff-f42d88422d9f50cb0863b6f08f2640a9e5cbcb9ac2ae01145642105d9fe9387d)).

And **enable detailed EM treatment in Geant4** by applying the following changes to the `ddsim` steering file (if you do not already use a `ddsim` steering file, you can create the default one with `ddsim --dumpSteeringFile > mySteeringFile.py`):
- Change the physics list to `SIM.physics.list = "FTFP_BERT_EMZ"`
- Change the range cut: `SIM.physics.rangecut = 0.05*mm`
- Remove the energy threshold for tracker hits: `SIM.filter.tracker = "edep0"`
- At the bottom of the file, change the Geant4 UI configure commands to: `SIM.ui.commandsConfigure = ["/cuts/setLowEdge 50 eV", "/process/em/lowestElectronEnergy 1 eV", "/process/em/auger true" , "/process/em/deexcitationIgnoreCut true"]`

Finally, for some BIB (e.g. IPC), the **boost due to the crossing angle has to be applied**:
- At the beginning of the file, use: `SIM.crossingAngleBoost = 0.015`

Whether or not to apply this boost depends on how the BIB was actually generated, reach out to the contact of the sample you plan to simulate in case of doubt.

NB: centrally maintained steering files are available for IDEA in [`FCC-Config`](https://github.com/HEP-FCC/FCC-config), which is included in the `key4hep` software stack and accessible with the environment variable `$FCCCONFIG`. For example, the [IDEA_o1_v03](https://github.com/HEP-FCC/FCC-config/blob/main/FCCee/FullSim/IDEA/IDEA_o1_v03/SteeringFile_IDEA_o1_v03.py) steering file can file can be passed to the ddsim command with:
```bash
ddsim --steeringFile $FCCCONFIG/share/FCC-config/FullSim/IDEA/IDEA_o1_v03/SteeringFile_IDEA_o1_v03.py ...
# or for the nightlies
ddsim --steeringFile $FCCCONFIG/FullSim/IDEA/IDEA_o1_v03/SteeringFile_IDEA_o1_v03.py ...
```
For CLD, the centrally maintained steering file lives [here](https://github.com/key4hep/CLDConfig/blob/main/CLDConfig/cld_steer.py) and can be accessed through `$CLDCONFIG`.

### Run the simulation

Example for processing through the ALLEGRO detector simulation an incoherent pair creation (IPC) file `pairs.pairs` (a text file with all particles' positions) generated with GuineaPig:

```bash
ddsim -N -1 \
 --inputFile /eos/experiment/fcc/users/a/aciarma/pairs/4IP_2024may29/Z/data1/pairs.pairs \
 --steeringFile mySteeringFile.py \
 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml \
 --outputFile sim_IPC_test.root
```

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


### submit_pairs.py

Set up the condor (or local) submission of simulation jobs of 
IPC background files.
At the moment, `.pairs` files contain a single event
and `ddsim` can process only one of them at the time.
The  `submit_pairs.py` script generates at list of bash scripts
to automatize the submission of many single event jobs.
After the preparation is done, the command to launch the jobs
(condor or locally) is printed in the terminal.

Example usage command:
```
submit_pairs.py --tag IDEA_my_test --compactFile $K4GEO/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml -n 10
```
which will prepare the submission for 10 events (jobs),
using the `IDEA_o1_v03` geometry description.
All the available geometries are stored in the
[`k4geo`](https://github.com/key4hep/k4geo/tree/main)
repository.


If a custom variation of the standard geometry 
tha requires recompiling k4geo is needed,
specify the path to the local build with the `--k4geo` flag.

All the available options can be seen using the `-h` flag.
For example the default input path is currently:
```
/eos/user/s/sfranche/FCC/BIB/data/aciarma_4IP_2024may29/Z/
``` 
but can be modified specifying `--input <your/path>`,
which should contain files with a naming format: `your/path/*_XYZ.pairs`,
where `XYZ` is an event number.


## Calorimeter calibration

TODO: add docu

## List of samples

TODO: add list/table


