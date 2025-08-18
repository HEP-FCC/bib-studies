# Plotting

Collections of scripts for making plots from BIB simulated files.

## Draw hit maps

There are 2 scripts in this directory for drawing the hit maps
- `drawhits.py` : it plots the hit maps and some profile histograms for a given sub-detector. Example (do `--help` to see all options):
  ```
  drawhits.py -i /eos/home-s/sfranche/FCC/BIB/data/aciarma_4IP_2024may29/Z/DDSim_output/bib_v1/ -n 10 -s VertexBarrel -d ALLEGRO_o1_v03_DetectorDimensions.json -h -m -z 1 -r 2
  ```
  Note that the `--detDictFile` (`-d`) argument is required and should point to a json file produced with the `xml2json.py` script.
  The  `-z` and `-r` arguments can be used to set different binning in the z and r axis for the hit maps in mm units.
- `plot_all_subdetectors.py`: it automatically loops over all the sub-detectors and the hot collections, matching them and reading the sub-detector dimensions directly from the geometry to define plot ranges. For each collection/sub-detector it calls the drawhits.py script (see above) to plot the hit map.


## xml2json

Extract detector dimensions from the geometry XML file into a json file.


## TODO
- [ ] In plot_all_subdetectors.py and xml2json when reading the max dimensions of each sub-detector. The results seems to be off by a factor 10 (except for vertex radius). Maybe cm vs mm issue?
-  In drawhits.py add more plots:
   -  [x] energy of hits in each layer
   -  [x] occupancy per layer
   -  [x] pt of associated MC particle
      - Cannot retrieve MC particles for calorimeter hits yet
- [x] Do simulation and plots also for signal Zmumu, Zee, Zqq
- [ ] Fix the `plot_all_subdetectors.py` with new updates to `drawhits.py` or add loop over sub-detector in `drawhits.py`


