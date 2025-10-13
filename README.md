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

### Example: other..
