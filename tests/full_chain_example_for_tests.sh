#!/bin/bash

######################################################
# run the SIM-DIGI-RECO step if the files do not exist yet
######################################################
if [ -f "ALLEGRO_sim.root" ]; then 
  echo "ALLEGRO_sim.root already exists, skipping the SIM step"
else 
  ddsim --enableGun --gun.distribution uniform --gun.energy "10*GeV" \
    --gun.particle e- --numberOfEvents 10 --outputFile ALLEGRO_sim.root \
    --random.enableEventSeed --random.seed 42 \
    --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml
fi

if [ -f "ALLEGRO_sim_digi_reco.root" ]; then
  echo "ALLEGRO_sim_digi_reco.root already exists, skipping the DIGI-RECO step" 
else 
  mkdir tmp_digi && cd tmp_digi
  cp /eos/project/f/fccsw-web/www/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/* .

  k4run $FCCCONFIG/share/FCC-config/FullSim/ALLEGRO/ALLEGRO_o1_v03/run_digi_reco.py \
    --IOSvc.Input ../ALLEGRO_sim.root \
    --IOSvc.Output ../ALLEGRO_sim_digi_reco.root
  cd ../
  rm -r tmp_digi 
fi

######################################################
# run local code to make sure it does not break
######################################################
#path_to_repo=$(git rev-parse --show-toplevel)

if [ -f "ALLEGRO_o1_v03_DetectorDimensions.json" ]; then
  echo "ALLEGRO_o1_v03_DetectorDimensions.json already exists, skipping the JSON generation step" 
else
  xml2json.py -d $K4GEO//FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml
fi 

drawhits.py -i ALLEGRO_sim_digi_reco.root -e 10 -s VertexBarrel -d ALLEGRO_o1_v03_DetectorDimensions.json


