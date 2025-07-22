#!/usr/bin/env bash

# Default values
BRANCH="stable"
RELEASE=""
COMMAND="source /cvmfs/sw.hsf.org/key4hep/setup.sh"

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
  COMMAND="source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh"
  shift
fi

# Check if release is specified
if [[ $# -eq 1 ]]; then
  RELEASE=$1
  COMMAND="$COMMAND -r $RELEASE"
fi

# Source the path
echo "Setting up $BRANCH build (release $RELEASE)..."
$COMMAND

# Install the python scripts
export BIB_DIR=$PWD
export PATH=$PATH:$PWD/simulation
export PATH=$PATH:$PWD/plotting