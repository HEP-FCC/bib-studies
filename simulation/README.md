# Simulation

Helper scripts for preparing and submitting simulations of background samples.

## set_vertex_000.py
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


## submit_pairs.py

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
submit_pairs.py --tag IDEA_my_test --geo IDEA_o1_v03 -n 10
```
which will prepare the submission for 10 events (jobs),
using the `IDEA_o1_v03` geometry description.
All the available geometries are stored in the
[`k4geo`](https://github.com/key4hep/k4geo/tree/main)
repository.
To specify a different one, use the convention `DETECTOR_oX_vYY`
(note that sometimes the `oX` is omitted e.g. for ILD).


If a custom variation of the standard geometry is needed, 
uncomment  [these lines](submit_pairs.py#L76-77)
and set the path to the customized version of the `k4geo` repo.
And keep the naming convention `DETECTOR_oX_vYY(_ZZZ).xml`.

All the available options can be seen using the `-h` flag.
For example the default input path is currently:
```
/eos/user/s/sfranche/FCC/BIB/data/aciarma_4IP_2024may29/Z/
``` 
but can be modified specifying `--input <your/path>`.

