# Simulation

Preparing and submitting simulations of background samples.

## Running detector simulation

Example for processing through the ALLEGRO detector simulation 
an incoherent pair creation (IPC) file `pairs.pairs` 
(a text file with all particles' positions) generated with GuineaPig:

```
ddsim -N -1 \
 --inputFile /eos/experiment/fcc/users/a/aciarma/pairs/4IP_2024may29/Z/data1/pairs.pairs \
 --crossingAngleBoost 0.015 \
 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml \
 --outputFile sim_IPC_test.root
```
Note that a crossing angle boost of 0.015 needs to be applied.

For other detectors, it may be required to also specify a steering file, with more elaborate
options for ddsim. Centrally maintained steering files are stored in the [`FCC-Config`](https://github.com/HEP-FCC/FCC-config) repo,
which is included in the `key4hep` software stack.
For example, the [IDEA_o1_v03](https://github.com/HEP-FCC/FCC-config/blob/main/FCCee/FullSim/IDEA/IDEA_o1_v03/SteeringFile_IDEA_o1_v03.py)
steering file can file can be passed to the ddsim command with:
```
ddsim --steeringFile $FCCCONFIG/share/FCC-config/FullSim/IDEA/IDEA_o1_v03/SteeringFile_IDEA_o1_v03.py ...
```

When running on signal sample, a smearing to the vertex need also to be included, with the option:
```
--vertexSigma 0.0098 2.54e-5 0.646 0.0063
```
(These numbers may change in the future)

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
specify the path to the local build with the `--k4geo` flag [NOT TESTED YET].

All the available options can be seen using the `-h` flag.
For example the default input path is currently:
```
/eos/user/s/sfranche/FCC/BIB/data/aciarma_4IP_2024may29/Z/
``` 
but can be modified specifying `--input <your/path>`.


## Calorimeter calibration

TODO: add docu

## List of samples

TODO: add list/table


