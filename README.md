# Beam induced background studies

Code and scripts developed to run beam induced background (BIB) studies.
The main folders contain dedicated `READMEs` with explanations of the scripts.


## Setup

The `setup.sh` script will setup this repository paths and 
the latest *stable* release of the `key4hep` software stack:
```
source setup.sh
```
To use the latest *nightly* add the `-n` flag:
```
source setup.sh -n
```
If you need a given release, specify it as the last argument, e.g.:
```
source setup.sh 2025-05-29
```
or for nightlies, e.g.:
```
source setup.sh -n 2025-07-22
```

