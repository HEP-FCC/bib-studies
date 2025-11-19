# Plotting

Collections of scripts for making plots from BIB simulated files.

## Draw hit maps

There are 3 scripts in this directory for drawing the hit maps
- `drawhits.py`: it plots the hit maps and some profile histograms for a given sub-detector.
- `hits2highLevelEstimations.py`: converts the results from `drawhits` into bandwidth, hit rate, pixel occupancy and other high-level estimations.
- **[TO BE FIXED]** `plot_all_subdetectors.py`: it automatically loops over all the sub-detectors and the hit collections, matching them. For each collection/sub-detector it calls the drawhits.py script to plot the hit map.

### drawhits.py

`drawhits.py` reads simulated events from PODIO/ROOT files and,
for a single specified sub-detector, loops over all hits, 
filling a set of summary histograms and 2D hit maps.
The script produces per-event and per-layer quantities such as zr/xy/zphi maps, 
timing distributions (including TOF-corrected), energy spectra, 
occupancy and pile-up estimators hits, plus basic MC particle kinematics; 
results are saved to a ROOT file and can also be exported as plots.

To decode the cell ID info, a detector specific `get_layer` function is used
(see for example in [`ALLEGRO.py`](../python/ALLEGRO.py)).

Key input arguments are:
- `-i` / `--infilePath` include the input file or directory (can also include wildcards).
- `-d` / `--detDictFile` path to the detector dimensions JSON produced by `xml2json.py`
- `-s` / `--subDetector` the target sub-detector name

Other useful options are:
- `-e` / `numberOfEvents` the number of events to process 
- `-n` / `numberOfFiles` the number of files to process
- `--sample` name of the sample, used as prefix of output files
- `-m` / `--draw_maps`, `-p` / `--draw_hists` whether to save maps or 1D histograms also in pdf
- `-r <width_in_mm>` , `-z  <width_in_mm>`, `--bin_width_phi  <width_in_rad>` bin widths for r/z/phi 
- `--integration_time <it_in_BX>` integration in units of bunch-crossings/events for pile-up counting
- `-a` / `--assumptions` path to the "assumptions" JSON dictionary
- `--digi` to use digitzed hits (collection defined in the assumptions dict)
- `--e_cut` to apply the energy cut (threshold defined in the assumptions dict)
- `--skip_layers` to skip creating and filling histograms per layer, useful for sub-detectors with many layers (e.g. drift chamber)

Use the `-h` / `--help` flag to print all the available arguments.

#### Examples
```sh
# Run over vertex barrel hits, producing hit maps with 1x1 mm^2 binning
drawhits.py -e 100 -s VertexBarrel -z 1 -r 1 -m -p -d $BIB_STUDIES/detectors_dicts/ALLEGRO_o1_v03_DetectorDimensions.json
```

```sh
# Run over ECal barrel digitized hits,
# applying a cut on the hit deposited energy (defined in the assumptions file)
drawhits.py -i <path_to_digitized_sim_hits_file> -e 100 -s ECalBarrel --digi --e_cut -d $BIB_STUDIES/detectors_dicts/ALLEGRO_o1_v03_DetectorDimensions.json -a $BIB_STUDIES/detectors_dicts/ALLEGRO_o1_v03_assumptions.json
```

#### drawhits.py histograms

- h_hit_x_mm_<sub-detector>: x position of all hits that matched this sub-detector
- ...

TODO: fill the list above


### hits2highLevelEstimations.py

`hits2highLevelEstimations.py` reads the per-layer and per-cell hit histograms produced by `drawhits.py`
and converts hit counts or occupancy into an estimated data bandwidth per layer, total bandwidth, 
hit rate and occupancy per layer and per cell (cell in semiconductor layers being a sensor module).
The conversion uses a chosen strategy ("hits" or "occupancy"),
a per-hit size (or per-layer sizes) and optional multipliers 
from an "assumptions" JSON dictionary, together with the detector channel counts 
from the detector JSON.
Results are written to a ROOT file and a simple 
plot of bandwidth vs layer is produced.

Key input arguments:
- `-i` / `--inputFile` path to the input ROOT file that contains the hit/occupancy histograms (output from drawhits.py).
- `-o` / `--outputFile` base name for the output ROOT file (sample name is prepended).
- `-d` / `--detDictFile` detector dimensions JSON (produced by `xml2json.py`) to get hits collection and per-layer channel counts.
- `-a` / `--assumptions` JSON with the bandwidth estimation assumptions: strategy (`hit_counts` or `occupancy`), `hit_size` (number or per-layer dict), and `multipliers` (multiplicative factors).
- `-r` / `--rate` hit rate in MHz used to scale counts to bandwidth.
- `--hitRateOccPlots` Use to also plot hit rate and occupancies per module in each layer (for vertex detector or silicon wrapper for example)

Notes:
- The assumptions JSON must define the chosen `strategy` and `hit_size`. If strategy is `occupancy`, number of cells per layer from the detector JSON are used to convert occupancy to absolute counts before applying hit size and rate.
- Use `-h` / `--help` to print all options.

Example:
```sh
# first produce the hits file
drawhits.py -e 100 -s VertexDisks --sample myhits

# then use the its output root file to estimate bandwidths
hits2highLevelEstimations.py -i myhits_hits_100evt_VertexDisks.root -d $BIB_STUDIES/detectors_dicts/ALLEGRO_o1_v03_DetectorDimensions.json -a $BIB_STUDIES/detectors_dicts/ALLEGRO_o1_v03_assumptions.json
```

## Draw particles info

The `drawparticles.py` script reads simulated events from PODIO/ROOT files and loops over all `MCParticles` to produce kinematic distributions and vertex origin maps. It categorizes particles into different collections (all, gen, sim, calorimeter, tracker) and produces per-collection histograms of key quantities like momentum, energy, eta, phi, and vertex positions. Results are saved to a ROOT file and can optionally be exported as plots.

Key input arguments are:
- `-i` / `--infilePath` input file or directory path
- `-n` / `--numberOfFiles` number of files to process
- `-e` / `--numberOfEvents` number of events per file
- `--sample` name of the sample, used as prefix for output files
- `-o` / `--outputFile` base name for the output ROOT file

## xml2json

The `xml2json.py` script reads a detector geometry XML file and extracts key dimensions, cell mappings and hits collection names into a JSON dictionary.
For each sensitive sub-detector, it retrieves properties like bounding box dimensions, detector type flags, and the number of cells per layer.

To extract the number of layers and the cells per layers, a detector specific `get_cells_map` function is used
(see for example in [`ALLEGRO.py`](../python/ALLEGRO.py)).

Key input arguments:
- `-d` / `--detGeoFile` path to the compact XML detector geometry file 

Example usage:
```sh
# After having sourced the setup.sh script, run:
xml2json.py -d $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml
# creates local .json file, with reduced geometry information 
# extracted/assumed from the detector XML, required for plotting+occupancy calculation,
# detector id, type, hits collection name, #cells per layer, max z/r>
```

```sh
# To run over the IDEA geometry
xml2json.py -d $K4GEO/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml
```

