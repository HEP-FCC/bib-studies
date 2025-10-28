# Beam induced background studies

Code and scripts developed to run beam induced background (BIB) studies,
following a _detector agnostic_ approach.

## Quick Start

```sh
# starting in an AFS/EOS folder with cvmfs+eos access, and subscription to
git clone git@github.com:HEP-FCC/bib-studies.git
cd bib-studies
source setup.sh

# create a working directory where to run the code
mkdir run
cd  run

drawhits.py -s VertexBarrel -e 10 -D 0
#will produce a root file with histograms
```

## Setup script

The `setup.sh` script will setup this repository paths and 
the latest *stable* release of the `key4hep` software stack:
```sh
source setup.sh
```
```sh
# or, for latest nightly:
source setup.sh -n
```
```sh
# for a given release:
source setup.sh 2025-05-29
```
```sh
# or for given nightly
source setup.sh -n 2025-07-22
```

## Repo structure

This repo contains four principal folders:
-  `plotting` for executable scripts concerning plots
-  `simulation` with some utility scripts and recipes for running simulations, and a list of relevant gen files.
-  `python` where some shared code is stored together with detector specific implementations of key functions used by the scripts.
-  `detectors_dicts` storing detectors related dictionaries.

The two folders with the executable scripts contain dedicated `READMEs` with more relevant info.
- [Simulation docs](./simulation/README.md)
- [Plotting docs](./plotting/README.md)



