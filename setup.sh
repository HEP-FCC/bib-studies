#!/usr/bin/env bash

# Default values
BRANCH="stable"
RELEASE=""
SW_PATH="/cvmfs/sw.hsf.org/key4hep/setup.sh"

_help() {
  echo "Usage: ${0} [-n] [release]"
  echo "  -h: Display this help message"
  echo "  -n: Use nightly builds (for development)"
  echo "  release: Specify release, default is latest"
}

# Print help
if [[ $1 == "-h" ]]; then
  _help; return 0;
fi

# Setup the key4hep software paths
if [[ $1 == "-n" ]]; then
  BRANCH="nightly"
  SW_PATH="/cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh"
  # Check if release is specified
  if [[ $# -eq 2 ]]; then RELEASE="${2}"; fi
else
  # Check if release is specified
  if [[ $# -eq 1 ]]; then RELEASE="${1}"; fi
fi

# Source the path
if [[ -n $RELEASE]]; then
  echo "Setting up $BRANCH build (release $RELEASE)..."
  source $SW_PATH -r $RELEASE
else
  echo "Setting up $BRANCH build..."
  source $SW_PATH
fi

# Install the python scripts
export PATH=$PATH:$PWD/scripts