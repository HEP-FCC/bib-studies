# Plotting

Collections of scripts for making plots from BIB simulated files.

## Draw hit maps

There are 2 scripts in this directory for drawing the hit maps
- `drawhits.py` : it plots the hit map in the (z,r) plane for a given collection of hits. Example (do `--help` to see all options):
```
drawhits.py -i /eos/home-s/sfranche/FCC/BIB/data/aciarma_4IP_2024may29/Z/DDSim_output/bib_v1/ -n 10 -c VertexBarrelCollection -m --map_binning '200,-200., 200., 200, -50., 50.'
```
- `plot_all_subdetectors.py`: it automatically loops over all the sub-detectors and the hot collections, matching them and reading the sub-detector dimensions directly from the geometry to define plot ranges. For each collection/sub-detector it calls the drawhits.py script (see above) to plot the hit map.


TODO:
1) In plot_all_subdetectors.py when reading the max dimensions of each sub-detector. The results seems to be off by a factor 10 (except for vertex radius). Maybe cm vs mm issue?
2) In drawhits.py add more plots: e.g. energy of hits in each layer,  pt of associated MC particle
3) Do simulation and plots also for signal Zmumu, Zee, Zqq


## xml2json

Extract detector dimensions from the geometry XML file into a json file.

