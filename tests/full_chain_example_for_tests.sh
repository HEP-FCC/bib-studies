#!/bin/bash

# never managed to get the repository path inside the GitHub runners... Moving to relative paths instead
#path_to_repo=${GITHUB_WORKSPACE:-$(git rev-parse --show-toplevel 2>/dev/null)}
path_to_repo=../

# run the SIM step if the files do not exist yet
[ -f "ALLEGRO_sim.root" ] && echo "ALLEGRO_sim.root already exists, skipping the SIM step" || ddsim --enableGun --gun.distribution uniform --gun.energy "10*GeV" --gun.particle e- --numberOfEvents 10 --outputFile ALLEGRO_sim.root --random.enableEventSeed --random.seed 42 --compactFile $K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml

# set-up the environment
if echo "$KEY4HEP_STACK" | grep -q nightlies; then
    echo "Nightlies detected"
    source $path_to_repo/setup.sh -n
    # run the DIGI-RECO step if the files do not exist yet
    [ -f "ALLEGRO_sim_digi_reco.root" ] && echo "ALLEGRO_sim_digi_reco.root already exists, skipping the DIGI-RECO step" || k4run $FCCCONFIG/FullSim/ALLEGRO/ALLEGRO_o1_v03/run_digi_reco.py --IOSvc.Input ALLEGRO_sim.root --IOSvc.Output ALLEGRO_sim_digi_reco.root --doTopoClustering False --CreateTruthLinks.OutputLevel 5

else
    echo "Stable release detected"
    source $path_to_repo/setup.sh
    # next lines to be removed when we will have a new (post 2025-05-29) stable release
    if [[ ! "$FCCCONFIG" == *share* ]]; then
      FCCCONFIG=$FCCCONFIG/share/FCC-config/
    fi
    # run the DIGI-RECO step if the files do not exist yet
    # get the files needed for calibration, noise, neighbor finding, etc
    if ! test -f ./DataAlgForGEANT.root; then  # assumes that if the last file exists, all the others as well
      wget https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/capacitances_ecalBarrelFCCee_theta.root
      wget https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/cellNoise_map_electronicsNoiseLevel_ecalB_thetamodulemerged.root
      wget https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/cellNoise_map_electronicsNoiseLevel_ecalB_thetamodulemerged_hcalB_thetaphi.root
      wget https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/cellNoise_map_endcapTurbine_electronicsNoiseLevel.root
      wget https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/cellNoise_map_electronicsNoiseLevel_ecalB_ECalBarrelModuleThetaMerged_ecalE_ECalEndcapTurbine_hcalB_HCalBarrelReadout_hcalE_HCalEndcapReadout.root
      wget https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/elecNoise_ecalBarrelFCCee_theta.root
      wget https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/lgbm_calibration-CaloClusters.onnx
      wget https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/lgbm_calibration-CaloTopoClusters.onnx
      wget https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/neighbours_map_ecalB_thetamodulemerged.root
      wget https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/neighbours_map_ecalE_turbine.root
      wget https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/neighbours_map_ecalB_thetamodulemerged_hcalB_hcalEndcap_phitheta.root
      wget https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/ALLEGRO/ALLEGRO_o1_v03/neighbours_map_ecalB_thetamodulemerged_ecalE_turbine_hcalB_hcalEndcap_phitheta.root
      wget https://fccsw.web.cern.ch/fccsw/filesForSimDigiReco/IDEA/DataAlgFORGEANT.root
    fi
    [ -f "ALLEGRO_sim_digi_reco.root" ] && echo "ALLEGRO_sim_digi_reco.root already exists, skipping the DIGI-RECO step" || k4run $FCCCONFIG/FullSim/ALLEGRO/ALLEGRO_o1_v03/run_digi_reco.py --IOSvc.Input ALLEGRO_sim.root --IOSvc.Output ALLEGRO_sim_digi_reco.root 
fi

# run local code to make sure it does not break
[ -f "ALLEGRO_o1_v03_DetectorDimensions.json" ] && echo "ALLEGRO_o1_v03_DetectorDimensions.json already exists, skipping the JSON generation step" || $path_to_repo/plotting/xml2json.py -d $K4GEO//FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml
$path_to_repo/plotting/drawhits.py -i ALLEGRO_sim_digi_reco.root -n 1 -s VertexBarrel -d ALLEGRO_o1_v03_DetectorDimensions.json -p -m -z 1 -r 2


