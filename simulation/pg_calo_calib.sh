#!/usr/bin/env bash

_usage() {
    echo "usage: ${0##*/} [-h] [options]"
}

_help() {
    _usage
    cat <<EOF

Options:
 -c compact : Path to detector's .xml compact file (Default is ALLEGRO_o1_v03)
 -n events  : Number of events to generate per sample

EOF
}

# Default options
COMPACT=$K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml
N_EVENTS=1000

while getopts ":hc:n:" opt $@;
do
    case $opt in
        h) _help; exit 1;;
        c) COMPACT=${OPTARG};;
        n) N_EVENTS=${OPTARG};;
        # handle errors
        \?) _usage; echo "Unknown option: -$OPTARG" >&2; exit 1;;
        :)  _usage; echo "Missing argument for -$OPTARG" >&2; exit 1;;
        *)  _usage; echo "Unimplemented option: -$OPTARG" >&2; exit 1;;
    esac
done
shift $((OPTIND-1))

DETECTOR=${COMPACT##*/}
DETECTOR=${DETECTOR%.xml}


#####################################################
## Barrel calib:  muon @ 30 GeV gun, theta 90 deg
#####################################################

ddsim --enableGun \
      --gun.distribution uniform \
      --gun.energy "30*GeV" \
      --gun.particle mu- \
      --gun.thetaMin "89*deg" \
      --gun.thetaMax "91*deg" \
      --gun.phiMin="-0.1*rad" \
      --gun.phiMax "0.1*rad" \
      --numberOfEvents $N_EVENTS \
      --random.enableEventSeed \
      --random.seed 42 \
      --compactFile $COMPACT \
      --outputFile ${DETECTOR}_mu30GeV_theta90deg.root


#####################################################
## Endcap calib:  muon @ 30 GeV gun, theta 15-45 deg
#####################################################

ddsim --enableGun \
      --gun.distribution uniform \
      --gun.energy "30*GeV" \
      --gun.particle mu- \
      --gun.thetaMin "5*deg" \
      --gun.thetaMax "45*deg" \
      --gun.phiMin="-0.1*rad" \
      --gun.phiMax "0.1*rad" \
      --numberOfEvents $N_EVENTS \
      --random.enableEventSeed \
      --random.seed 42 \
      --compactFile $COMPACT \
      --outputFile ${DETECTOR}_mu30GeV_theta5_45deg.root
