#!/usr/bin/env python

"""
  Script for submitting condor jobs to process .pairs, .root and .hepevt files
  through a `ddsim` simulation step.

  Generalized from submit_pairs.py
"""

import argparse
import os
import re
import stat
import sys
from pathlib import Path

import ROOT

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.append(str(PROJECT_ROOT))
from python import edm4hep_helper
# Default paths / namings
input_def_path = "/eos/home-s/sfranche/FCC/samples/bib/gen-samples/aciarma_4IP_2024may29/Z/"
output_def_folder = input_def_path+"DDSim_output/"


# Argument parser
parser = argparse.ArgumentParser('Submit condor jobs for bkg files.')
parser.add_argument('-i', '--input', default=input_def_path,
                    help='input path. Default is: '+input_def_path)
parser.add_argument('-o', '--output', default=output_def_folder,
                    help='output folder. Default is: '+output_def_folder)
parser.add_argument('-t', '--tag',
                    default='ddsim_output',
                    help='Tag of the dataset.')
parser.add_argument('-n', '--n_max_jobs', default=-1, type=int,
                    help='Maximum number of jobs.')
parser.add_argument('-c', '--compactFile', default="$K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml", type=str,
                    help='Detector geometry.')
parser.add_argument('-s', '--steering_file', default=None, type=str,
                    help='Path to steering file (if not given, attempt to find default from $FCCConfig.')
parser.add_argument('-k', '--k4geo', default=None, type=str,
                    help='Path to custom k4geo.')
parser.add_argument('--crossingAngleBoost', default=0.015, type=str,
                    help='Crossing angle boost to be applied. Required (+15 mrad) for '
                         'the SR samples, which are given in the beamline frame.')
parser.add_argument('--mcCollection', default="MCParticles", type=str,
                    help='MCParticle collection name inside edm4hep .root inputs. '
                         '"MCParticles" for the SR samples, "Pairs" for guinea-pig IPC.')
parser.add_argument('-j', '--job_flavor', default="espresso", type=str,
                    help='Job flavor for Condor submission.')
parser.add_argument('--seed', default=42, type=int,
                    help='Random seed passed to ddsim (--random.seed). Fixed by default '
                         'so re-running a file reproduces the same result.')


# Condor command content with custom JobFlavor
condor_cmd_content = """executable     = $(filename)
# for debugging:
# redirect the log file to somewhere accessible (uncomment lines below)
#Log            = $(filename).log
#Output         = $(filename).out
#Error          = $(filename).err
Log            = $(CONDOR_JOB_ID).log
Output         = $(CONDOR_JOB_ID).out
Error          = $(CONDOR_JOB_ID).err
requirements    = ( (OpSysAndVer =?= "AlmaLinux9") && (Machine =!= LastRemoteHost) && (TARGET.has_avx2 =?= True) )
max_retries    = 3
+JobFlavour    = "{1}" 
request_memory = 8GB 
RequestCpus = 1
queue filename matching files {0}
"""

# Local command content
local_cmd_content = """#!/bin/bash
SCRIPTS=$(ls {0})

for script in $SCRIPTS; do
    echo "#########################################"
    echo "Running $script..."
    echo "#########################################"
    ./$script
done
"""

# Header of executable script
fcc_cfg = os.environ["FCCCONFIG"]
if "sft-nightlies.cern.ch/lcg" in fcc_cfg:
    # On AlmaLinux9, key4hep nightlies' setup.sh defaults to sourcing the
    # rolling LCG "devkey-head" view directly from CVMFS
    fcc_ver = fcc_cfg.split("/")[6]             # get the nightly weekday (e.g. "Tue")
    exec_header = """#!/bin/bash
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh --lcg
"""
else:
    fcc_dir = "/".join(fcc_cfg.split("/")[:4])  # get software stack directory
    fcc_ver = fcc_cfg.split("/")[5]             # get the release number
    exec_header = f"""#!/bin/bash
source {fcc_dir}/setup.sh -r {fcc_ver}
"""

k4geo_path="""
# For using a local version of K4GEO
cd GEO_PATH
k4_local_repo
cd -
"""

# Path to steering files (skipped if -s option set)
steering_dict = {
    "IDEA_o1_v03":  "$FCCCONFIG/share/FCC-config/FullSim/IDEA/IDEA_o1_v03/SteeringFile_IDEA_o1_v03.py"
}

def run(args):

    tag = args.tag
    input_file_path = args.input
    output_file_path = args.output
    n_max = args.n_max_jobs
    compact = args.compactFile
    k4geo = args.k4geo
    x_angle = args.crossingAngleBoost
    job_flavor = args.job_flavor  


    # Get the short name of geometry file
    geo = compact.split("/")[-1].strip(".xml")
    #detector = re.sub("(_o[0-9]+)?_v[0-9]{2}.*", "", geo)

    # Define output storage path
    storage_path = os.path.join(output_file_path, tag)

    print("Creating output storage path:")
    print(storage_path)
    os.makedirs(storage_path, exist_ok=True)

    print("Creating submission folder:", tag)
    os.makedirs(tag, exist_ok=True)

    # Check if custom k4geo is to be used
    header = exec_header
    if k4geo != None:
        header += k4geo_path.replace("GEO_PATH",k4geo)

    # Check if a steering file is required
    steering_opt = ""
    if args.steering_file:
        steering_opt = "--steeringFile "+args.steering_file
    else:
        if geo in steering_dict:
            steering_path = steering_dict[geo]
            print("Including steering file: ", steering_path)
            steering_opt = f"--steeringFile {steering_path}"

    # sorted(): an -n subset must be the same set on every run
    items = sorted(os.listdir(input_file_path))

    # Setup the bash executables scripts
    print("Preparing submission for:")
    n_jobs = 0
    submitted_root_files = []
    exec_template_name = os.path.join(tag, "run_ddsim_FILENAME.sh")


    # NEW: Resolve the ddsim -N value for edm4hep inputs
    root_files = [os.path.join(input_file_path, i) for i in items if i.endswith(".root")]
    if n_max > 0:
        root_files = root_files[:n_max]
    if len(root_files) > 0:
        root_events = edm4hep_helper.count_events(root_files[0])   # WARNING COMMENT


    for item in items:
        if (n_max > 0) and (n_max <= n_jobs):
            break

        #TODO: the input name path could be handled a bit more gracefully

        # check if the item is a folder, following naming
        # convention used in a. ciarma's samples
        item_path = os.path.join(input_file_path, item)
        bx_id = None
        input_filename = None
        # -1 means "all events in the file", understood by the .pairs and
        # .hepevt readers but not by the edm4hep one, which needs the count.
        n_events = -1
        reader_opt = ""
        
        # 1. Check for older folder-based IPC samples
        if os.path.isdir(item_path):
            if "data" not in item:
                continue
            bx_id = item.replace("data", "")
            input_filename = os.path.join(input_file_path, item, "pairs.pairs")
            print("- "+input_filename)
            
        # 2. Check for file-based IPC samples (.pairs)
        elif item.endswith(".pairs"):
            input_filename = os.path.join(input_file_path, item)
            bx_id = re.search(r"_[0-9]+\.",item).group(0).strip("_.")
            print("- "+input_filename)
            
        # 3. Check for Synchrotron Radiation samples (.hepevt)
        elif item.endswith(".hepevt"):
            input_filename = os.path.join(input_file_path, item)
            # Extracts just the number (e.g., '100004') from 'output_100004.hepevt'
            bx_id = item.replace("output_", "").replace(".hepevt", "")
            print("- "+input_filename)

        # 4. NEW: Check for EDM4hep samples (.root), SR and IPC
        elif item.endswith(".root"):
            input_filename = os.path.join(input_file_path, item)
            # e.g. '100093_edm4hep'
            bx_id = item.replace("output_", "").replace(".root", "")

            # The edm4hep reader needs the exact event count, it does not
            # understand -N -1. Resolved for the whole sample before the loop.
            n_events = root_events
            submitted_root_files.append(input_filename)

            # SR halo calls it MCParticles, guinea-pig IPC calls it Pairs.
            reader_opt = f"--edm4hep.mcParticleCollectionName {args.mcCollection} "
            print("- "+input_filename)

        else:
            print("Skipping item:", item)
            continue

        command = header
        executable_path = exec_template_name.replace("FILENAME", bx_id)
        output_filename = os.path.join(storage_path, f"{geo}_{bx_id}_r{fcc_ver}.root")

        # for performance, write the output locally first and copy at the end
        tmp_output_filename = os.path.basename(output_filename)

        command += f"""ddsim \
            --compactFile  {compact} \
            -I {input_filename} \
            -O {tmp_output_filename} \
            -N {n_events} --crossingAngleBoost {x_angle} --random.seed {args.seed} \
            --part.keepAllParticles True  {reader_opt}{steering_opt}\n"""

        command += f"mv {tmp_output_filename} {output_filename}"

        with open(executable_path, "w") as f:
            f.write(command)

        st = os.stat(executable_path)
        os.chmod(executable_path, st.st_mode | stat.S_IEXEC)
        n_jobs += 1

    if submitted_root_files:
        edm4hep_helper.print_normalisation(submitted_root_files, root_events)

    # Setup the condor script
    condor_submit_path = f"{tag}.cmd"
    exec_pattern = exec_template_name.replace("FILENAME", "*")
    cmd_file_content = condor_cmd_content.format(exec_pattern, job_flavor)

    with open(condor_submit_path, "w") as f:
        f.write(cmd_file_content)
    submit_cmd = f"condor_submit {condor_submit_path}"

    print("To submit the condor job: ", submit_cmd)

    # Local run script
    local_submit_path = f"{tag}.sh"
    cmd_file_content = local_cmd_content.format(exec_pattern)
    with open(local_submit_path, "w") as f:
        f.write(cmd_file_content)
    print("To submit the local job: sh ", local_submit_path)


if __name__ == "__main__":
    args = parser.parse_args()
    run(args)
