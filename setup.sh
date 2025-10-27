#!/usr/bin/env bash

# Print help
_help() {
  echo "Usage: ${0} [-n] [release]"
  echo "  -h: Display this help message"
  echo "  -n: Use nightly builds (for development)"
  echo "  release: Specify release, default is $RELEASE for stable and TODAY for nightly "
}

if [[ $1 == "-h" ]]; then
  _help; return 0;
fi

# Path to this setup.sh script
export BIB_STUDIES=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd -P )

# Install the repo paths
export PATH=$BIB_STUDIES/simulation:$PATH
export PATH=$BIB_STUDIES/plotting:$PATH
export PYTHONPATH=$BIB_STUDIES/python:$PYTHONPATH

# Default setup is stable release
COMMAND="source /cvmfs/sw.hsf.org/key4hep/setup.sh"

# Setup the nightly
if [[ $1 == "-n" ]]; then
  COMMAND="source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh"
  shift
fi

# Check if release is specified
if [[ $# -eq 1 ]]; then
  COMMAND="$COMMAND -r $1"
fi

echo "Setting up FCC software stack..."
$COMMAND
