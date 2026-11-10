# Plotting

Collections of scripts for making plots from BIB simulated files.

## Draw hit maps

There are 3 scripts in this directory for drawing the hit maps
- `drawhits.py`: it plots the hit maps and some profile histograms for a given sub-detector.
- `hits2highLevelEstimations.py`: converts the results from `drawhits` into bandwidth, hit rate, pixel occupancy and other high-level estimations.
- `plot_all_subdetectors.py`: runs `drawhits.py` (and optionally `hits2highLevelEstimations.py`) for all the sub-detectors of an assumptions file, several in parallel, with one log per sub-detector and a summary table at the end (see below).

### drawhits.py

`drawhits.py` reads simulated events from PODIO/ROOT files and,
for a single specified sub-detector, loops over all hits, 
filling a set of summary histograms and 2D hit maps.
The script produces per-event and per-layer quantities such as zr/xy/zphi maps, 
timing distributions (including TOF-corrected), energy spectra, 
occupancy and pile-up estimators hits, plus basic MC particle kinematics; 
results are saved to a ROOT file and can also be exported as plots.
With the per-layer histograms enabled (no `--skip_layers`), the hit-rate density vs |eta| of all the layers is also overlaid on one canvas,
saved in `per_layer/rate_vs_eta/` of the ROOT file (and as `<sample>_hitRateDensity_VS_eta_allLayers_<N>evt_<sub-detector>.pdf` with `-p`), to compare the layers.
A second y axis on the right shows the density scaled by the example rate assumptions (52 MHz bunch crossing rate, 5 pixels per cluster for the pixel detectors, see [`constants.py`](../python/constants.py)).

To decode the cell ID info, a detector specific `get_layer` function is used
(see for example in [`ALLEGRO.py`](../python/ALLEGRO.py)).

Key input arguments are:
- `-i` / `--infilePath` include the input file or directory (can also include wildcards).
- `-d` / `--detDictFile` path to the detector dimensions JSON produced by `xml2json.py`
- `-s` / `--subDetector` the target sub-detector name

Other useful options are:
- `-e` / `numberOfEvents` the number of events to process per file (`-1` for all events)
- `-n` / `numberOfFiles` the number of files to process (`-1` for all files)
- `--sample` name of the sample, used as prefix of output files
- `-m` / `--draw_maps`, `-p` / `--draw_hists` whether to save maps or 1D histograms also in pdf
- `-r <width_in_mm>` , `-z  <width_in_mm>`, `--bin_width_phi  <width_in_rad>` bin widths for r/z/phi 
- `--integration_time <it_in_BX>` integration window in units of bunch-crossings/events for pile-up counting: the first hit in a cell opens the window, the later events in it with a hit in the same cell count as pile-up (1 event = 1 BX, default 4, 1 = off). The pile-up is counted per cell ID, i.e. per pixel for the silicon pixel detectors (vertex, silicon wrapper), and the per-layer histograms are in `per_layer/pileup_pixel/` or `per_layer/pileup_cell/` accordingly
- `--bc_window <W>` fill the distribution of the time (in events = BX) since the previous hit in the same module (in the same cell for sub-detectors without a module field, e.g. the drift chamber): `per_layer/h_bc_window_x_layer_<collection>` (layer vs time) and, for sub-detectors with modules and without `--skip_layers`, `per_layer/bc_window_module/h_bc_window_layer<N>_<collection>` (module vs time). The y axis is 1 to W events, W+1 = more than W events (the first hit of each module/cell in the run is not counted). Default 40, -1 = off. Per-layer module histograms are skipped for layers with more than 20000 modules (calorimeters), the layer-vs-time summary is always filled.
- `-a` / `--assumptions` path to the "assumptions" JSON dictionary
- `--digi` to use digitzed hits (collection defined in the assumptions dict)
- `--collection <name>` to read a different hits collection, overriding the detector/assumptions dict (e.g. `OverlayDCHCollection` for overlaid but non-digitized files). If no CellIDEncoding is stored for it, the one of the detector dict `hitsCollection` is used
- `--e_cut` to apply the energy cut (threshold defined in the assumptions dict)
- `--skip_layers` to skip creating and filling histograms per layer, useful for sub-detectors with many layers (e.g. drift chamber)

- `--plot_primary` to plot the primary particle and direct parent of the particle creating each hit (needs the `MCParticles` collection, see below)

Use the `-h` / `--help` flag to print all the available arguments.

**Inputs without `MCParticles`** (e.g. merged BIB samples where the collection was dropped to save space):
the script still runs and fills all hit-based histograms as usual.
The hits' MC particle links cannot be resolved (calling any getter on them segfaults),
so each link is checked with `isAvailable()` before use. Hits without a link do not fill the MC-truth
histograms (`h_hit_particle_*` and the `--plot_primary` primary/parent ones), and a warning reports how many there were.
If no hit has a link, those histograms are neither drawn nor written to the output ROOT file.

**Normalisation**: all histograms are normalised per event. The number of events used is the number requested
(`-e` times the number of input files), capped to the total number of events available in the input files,
so asking for more events than exist no longer dilutes the results. With `-e -1` all events of all the input files are used.

**Total occupancy binning**: `h_occ_anyhit_tot_<collection>` and `h_occ_allhits_tot_<collection>` have fine bins at low occupancy, in units of fired cells or hits
(0, 1, 2-3 and 4-5), followed by the usual logarithmic bins, so the sparse-occupancy regime (e.g. the drift chamber) can be resolved.

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
- `--hitRateOccStats` Print mean, standard deviation and 95th percentile of the per-layer occupancy.

Notes:
- The assumptions JSON must define the chosen `strategy` and `hit_size`. If strategy is `occupancy`, number of cells per layer from the detector JSON are used to convert occupancy to absolute counts before applying hit size and rate.
- The input file must come from the current `drawhits.py` folder layout, the histograms are read from `per_layer/`, `per_layer/module_hits/` and `per_layer/occupancy_allhits/`. Older files are not supported, regenerate them with `drawhits.py`.
- The `occupancy` strategy and the `--hitRateOccStats` printout use the all-hits occupancy (`h_avg_occ_allhits_x_layer`, `h_occ_allhits_x_layer<N>`), i.e. every hit counts, also repeated hits in the same cell and event. With that choice, the `occupancy` strategy gives the same bandwidth as `hit_counts` with the same hit size.
- Use `-h` / `--help` to print all options.

Example:
```sh
# first produce the hits file
drawhits.py -e 100 -s VertexDisks --sample myhits

# then use the its output root file to estimate bandwidths
hits2highLevelEstimations.py -i myhits_hits_100evt_VertexDisks.root -d $BIB_STUDIES/detectors_dicts/ALLEGRO_o1_v03_DetectorDimensions.json -a $BIB_STUDIES/detectors_dicts/ALLEGRO_o1_v03_assumptions.json
```

### plot_all_subdetectors.py

Loops over the sub-detectors listed in the assumptions JSON (those that are also in the detector dictionary) and runs `drawhits.py` on each of them,
`-j` sub-detectors at a time. The output of each run goes to `log_<sample>_<sub-detector>.txt`, and a summary with the status and the time of each sub-detector is printed at the end.
Options not known to this script are forwarded to `drawhits.py`.

Key input arguments:
- `-i` / `--inputFile` input file or directory (required), `-n`, `-e` number of files and of events per file
- `-d` / `--detDictFile`, `-a` / `--assumptions` the detector dimensions and assumptions JSON dictionaries (default: ALLEGRO_o1_v03)
- `-j` / `--jobs` number of sub-detectors processed in parallel
- `--only`, `--exclude` comma-separated lists of sub-detectors to run over / skip
- `--skip_layers_for` comma-separated list of sub-detectors for which `--skip_layers` is added (e.g. `DCH_v2`, which has too many layers)
- `--hits2high` also run `hits2highLevelEstimations.py` on each output, the total bandwidth is added to the summary
- `--dry_run` only print the `drawhits.py` commands

Example:
```sh
plot_all_subdetectors.py -i myfile.root -e 100 -j 4 --hits2high --skip_layers_for DCH_v2 -p -m
```

## Draw MC particles info

The `drawMCparticles.py` script reads the MC particles of the events (`MCParticles` of simulated files, or another collection with `--mcCollection`) and produces kinematic distributions and production vertex maps.
It is independent of the hits: `drawhits.py` only histograms the particles that created a hit, while this script looks at all the particles. It categorizes them into all, gen (generator-level),
sim (created in the simulation), calorimeter and tracker (decaying in the calorimeter / tracker), and produces per-category histograms of energy, p<sub>T</sub>, eta, phi and PDG ID,
plus the production vertex (r, z, phi, x-y, z-r, z-phi) of the particles created in the simulation.
Results are saved to a ROOT file (one folder per category, plus `sim_vertex`) and as plots.

It can also be run on generator-level files, to inspect the BIB particles before the simulation, e.g. for the incoherent pair creation edm4hep files:
```sh
drawMCparticles.py -i <path_to_generator_file.root> --mcCollection Pairs -e 100 --sample pairs
```

Key input arguments are:
- `-i` / `--infilePath` input file or directory path (required)
- `-n` / `--numberOfFiles` number of files to process
- `-e` / `--numberOfEvents` number of events per file, -1 for all
- `--mcCollection` name of the MC particle collection (default `MCParticles`)
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

