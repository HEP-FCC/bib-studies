# Beam induced background studies

Code and scripts developed to run beam induced background (BIB) studies.
The main folders contain dedicated `READMEs` with explanations of the scripts.

## Setup

The `setup.sh` script will setup this repository paths and 
the latest *stable* release of the `key4hep` software stack:
```sh
source setup.sh
#or, for latest nightly:
source setup.sh -n

#for a given release:
source setup.sh 2025-05-29
#or for given nightly
source setup.sh -n 2025-07-22
```

## Repo structure

...

## Getting started

### Example: Hits & occupancy plots for sub-detector

```sh
#starting in an AFS/EOS folder with cvmfs+eos access
git clone g clone git@github.com:HEP-FCC/bib-studies.git
cd bib-studies
source setup.sh
cd plotting
drawhits.py -s MuonTaggerBarrel -d ALLEGRO_o1_v03_DetectorDimensions.json -e 10 -D 0
#will produce a root file with histograms
```

#### drawhits.py histograms

- h_hit_xx_mm_<sub-detector>: x position of all hits that matched this sub-detector
- ...

### Example: Updated reduced geometry json file

```sh
#see above for setup..
cd plotting
xml2json.py -d $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml
#creates local .json file, with reduced geometry information extracted/assumed from the detector XML, required for plotting+occupancy calculation
#detector id, type, hits collection name, #cells per layer, max z/r
```

### Example: other..
