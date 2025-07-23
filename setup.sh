#!/usr/bin/env bash

# Path were the setup.sh script is located
export BIB_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd -P )

# Install this repo paths
export PATH=$PATH:$BIB_DIR/simulation
export PATH=$PATH:$BIB_DIR/plotting


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
