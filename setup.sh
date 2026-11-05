#!/usr/bin/env bash

# Source this script to set up the environment of the BIB studies:
# the paths of this repository and a key4hep software stack.

_bib_help() {
  cat << EOF
Usage: source setup.sh [-n [--spack | --lcg]] [-l] [-r <release> | <release>]
  -h, --help     display this help message
  -n, --nightly  use the nightly builds (for development). By default the nightly stack of the key4hep setup (LCG on AlmaLinux 9,
                 the release is then a weekday Mon-Sat or "latest", and it is overwritten every week). A release given as a date
                 (e.g. 2026-10-02) selects the spack nightly of that day, which does not change, so it is reproducible
  --spack        with -n: force the spack nightlies (the release is a date, the latest one if not specified).
                 NOTE: they are no longer built upstream (last one 2026-10-20) and the --spack option of key4hep will soon fail
  --lcg          with -n: force the LCG nightly view (release: weekday or "latest")
  -l, --list     list the available releases and return (stable releases, or spack nightlies with -n)
  -r, --release  release to set up, or simply give it as last argument (default is the latest)

Examples:
  source setup.sh                   # latest stable release
  source setup.sh 2025-05-29        # a given stable release
  source setup.sh -n                # latest nightly (default nightly stack)
  source setup.sh -n 2026-10-02     # the spack nightly of 2026-10-02 (the --spack option is added automatically)
  source setup.sh -n Tue            # the LCG nightly of Tuesday
  source setup.sh -n -l             # list the available spack nightlies
EOF
}

# This script modifies the current shell: it has to be sourced
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  echo "This script must be sourced: source ${0}"
  exit 1
fi

# Parse the options
_bib_nightly=0; _bib_lcg=0; _bib_spack=0; _bib_list=0; _bib_release=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)    _bib_help; unset -f _bib_help; return 0 ;;
    -n|--nightly) _bib_nightly=1 ;;
    --lcg)        _bib_lcg=1 ;;
    --spack)      _bib_spack=1 ;;
    -l|--list)    _bib_list=1 ;;
    -r|--release) _bib_release="$2"; shift ;;
    -*)           echo "Unknown option: $1"; _bib_help; unset -f _bib_help; return 1 ;;
    *)            _bib_release="$1" ;;
  esac
  shift
done

# Choose the key4hep setup script and its options
if [[ $_bib_nightly -eq 1 ]]; then
  _bib_setup=/cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
  _bib_flags=()  # default stack of the nightly setup
  # a date can only be a spack nightly (the LCG ones are selected by weekday)
  if [[ $_bib_spack -eq 1 || $_bib_release =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ || ( $_bib_list -eq 1 && $_bib_lcg -eq 0 ) ]]; then
    _bib_flags=(--spack)
  fi
  [[ $_bib_lcg -eq 1 ]] && _bib_flags=(--lcg)
  if [[ $_bib_lcg -eq 1 && $_bib_spack -eq 1 ]]; then
    echo "--lcg and --spack are exclusive"; unset -f _bib_help; return 1
  fi
else
  _bib_setup=/cvmfs/sw.hsf.org/key4hep/setup.sh
  _bib_flags=()
  if [[ $_bib_lcg -eq 1 || $_bib_spack -eq 1 ]]; then
    echo "--lcg and --spack are only available for the nightlies (-n)"; unset -f _bib_help; return 1
  fi
fi

if [[ $_bib_list -eq 1 ]]; then
  if [[ $_bib_lcg -eq 1 ]]; then
    echo "The LCG nightlies are selected by weekday (Mon-Sat) or latest, there is nothing to list"; unset -f _bib_help; return 1
  fi
  echo "Available releases:"
  source "$_bib_setup" "${_bib_flags[@]}" --list-releases
  _bib_status=$?
else
  # Path to this setup.sh script
  export BIB_STUDIES=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd -P )

  # Install the repo paths
  export PATH=$BIB_STUDIES/simulation:$PATH
  export PATH=$BIB_STUDIES/plotting:$PATH
  export PYTHONPATH=$BIB_STUDIES/python:$PYTHONPATH

  [[ -n "$_bib_release" ]] && _bib_flags+=(-r "$_bib_release")
  echo "Setting up the FCC software stack: source $_bib_setup ${_bib_flags[*]}"
  source "$_bib_setup" "${_bib_flags[@]}"
  _bib_status=$?
fi

# Do not leave helper variables and functions in the shell of the user
unset -f _bib_help
unset _bib_nightly _bib_lcg _bib_spack _bib_list _bib_release _bib_setup _bib_flags
eval "unset _bib_status; return $_bib_status"
