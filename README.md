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
