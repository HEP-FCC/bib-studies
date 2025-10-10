# Plotting

Collections of scripts for making plots from BIB simulated files.

## Draw hit maps

There are 2 scripts in this directory for drawing the hit maps
- `drawhits.py` : it plots the hit maps and some profile histograms for a given sub-detector. Example (do `--help` to see all options):
  ```sh
  #make sure simplified geometry (json) file is up to date (default used is local ALLEGRO example) 
  drawhits.py -e 1 -s VertexBarrel
  ```
  Note that the `--detDictFile` (`-d`) argument is required and should point to a json file produced with the `xml2json.py` script.
  The  `-z` and `-r` arguments can be used to set different binning in the z and r axis for the hit maps in mm units.
- `plot_all_subdetectors.py`: it automatically loops over all the sub-detectors and the hot collections, matching them and reading the sub-detector dimensions directly from the geometry to define plot ranges. For each collection/sub-detector it calls the drawhits.py script (see above) to plot the hit map.


## Draw particles info

The `drawparticles.py` script loops over the entire `MCParticles` collection making several kinematic plots.

## xml2json

Extract detector dimensions from the geometry XML file into a json file. E.g:

```sh
xml2json.py -d $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml
#or
xml2json.py -d $K4GEO/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml
```

## TODO
- [ ] In plot_all_subdetectors.py and xml2json when reading the max dimensions of each sub-detector. The results seems to be off by a factor 10 (except for vertex radius). Maybe cm vs mm issue?
-  In drawhits.py add more plots:
   -  [x] energy of hits in each layer
   -  [x] occupancy per layer
   -  [x] pt of associated MC particle
      - Cannot retrieve MC particles for calorimeter hits yet
- [x] Do simulation and plots also for signal Zmumu, Zee, Zqq
- [ ] Fix the `plot_all_subdetectors.py` with new updates to `drawhits.py` or add loop over sub-detector in `drawhits.py`


