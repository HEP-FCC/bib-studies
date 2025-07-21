#!/usr/bin/env bash

# Setup latest stable or nightly release.
# New stable version will appear in:
# /cvmfs/sw.hsf.org/key4hep/releases/


# Default values
BRANCH="stable"
RELEASE="2025-05-29"
SW_PATH="/cvmfs/sw.hsf.org/key4hep/setup.sh"

_help() {
  echo "Usage: ${0} [-n] [release]"
  echo "  -h: Display this help message"
  echo "  -n: Use nightly builds (for development)"
  echo "  release: Specify release, default is $RELEASE for stable and TODAY for nightly "
}

# Print help
if [[ $1 == "-h" ]]; then
  _help; return 0;
fi

# Setup the key4hep software paths
if [[ $1 == "-n" ]]; then
  BRANCH="nightly"
  SW_PATH="/cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh"
  RELEASE=$(date +"%Y-%m-%d")
  shift
fi

# Check if release is specified
if [[ $# -eq 1 ]]; then 
  RELEASE=$1
fi

# Source the path
echo "Setting up $BRANCH build (release $RELEASE)..."
source $SW_PATH -r $RELEASE

# Install the python scripts
export BIB_DIR=$PWD
export PATH=$PATH:$PWD/simulation
export PATH=$PATH:$PWD/plotting