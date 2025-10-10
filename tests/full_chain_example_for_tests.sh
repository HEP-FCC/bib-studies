#!/bin/bash

# never managed to get the repository path inside the GitHub runners... Moving to relative paths instead
#path_to_repo=${GITHUB_WORKSPACE:-$(git rev-parse --show-toplevel 2>/dev/null)}
path_to_repo=../

# set-up the environment
if echo "$KEY4HEP_STACK" | grep -q nightlies; then
    echo "Nightlies detected"
    source $path_to_repo/setup.sh -n
else
    echo "Stable release detected"
    source $path_to_repo/setup.sh
fi

# run the SIM-DIGI-RECO step if the files do not exist yet
[ -f "ALLEGRO_sim.root" ] && echo "ALLEGRO_sim.root already exists, skipping the SIM step" || ddsim --enableGun --gun.distribution uniform --gun.energy "10*GeV" --gun.particle e- --numberOfEvents 10 --outputFile ALLEGRO_sim.root --random.enableEventSeed --random.seed 42 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml
[ -f "ALLEGRO_sim_digi_reco.root" ] && echo "ALLEGRO_sim_digi_reco.root already exists, skipping the DIGI-RECO step" || k4run $FCCCONFIG/FullSim/ALLEGRO/ALLEGRO_o1_v03/run_digi_reco.py --IOSvc.Input ALLEGRO_sim.root --IOSvc.Output ALLEGRO_sim_digi_reco.root --doTopoClustering False --CreateTruthLinks.OutputLevel 5

# run local code to make sure it does not break
[ -f "ALLEGRO_o1_v03_DetectorDimensions.json" ] && echo "ALLEGRO_o1_v03_DetectorDimensions.json already exists, skipping the JSON generation step" || $path_to_repo/plotting/xml2json.py -d $K4GEO//FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml
$path_to_repo/plotting/drawhits.py -i ALLEGRO_sim_digi_reco.root -n 1 -s VertexBarrel -d ALLEGRO_o1_v03_DetectorDimensions.json -p -m -z 1 -r 2


