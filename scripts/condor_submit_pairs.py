#!/usr/bin/env python

"""
  Script for submitting condor jobs to process .pairs files
  through a `ddsim` simulation step. These files are currently
  stored one per folder, named .../dataXXX/pairs.pairs.

  Original script from Brieuc:
  https://github.com/BrieucF/FCC_scripts/blob/main/background_studies/launch_bkg_simulation_v23.py
"""

import argparse
import os
import stat


# Condor command content
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
+JobFlavour    = "espresso"
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

# Default header of executable script
exec_header = """#!/bin/bash
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh

# For using custom geometries, a local version of K4GEO can be set
# following the pattern below
#cd /afs/cern.ch/user/your/local/k4geo
#k4_local_repo
#cd -
"""

# Default paths / namings
input_def_path = "/eos/user/s/sfranche/FCC/BIB/data/aciarma_4IP_2024may29/Z/"
output_def_folder = "DDSim_output"


# Argument parser
parser = argparse.ArgumentParser('Submit condor jobs for IPC files.')
parser.add_argument('-i', '--input', default=input_def_path,
                    help='input path. Default is: '+input_def_path)
parser.add_argument('-o', '--output', default=output_def_folder,
                    help='output folder. Default is: '+output_def_folder)
parser.add_argument('-t', '--tag',
                    help='Tag of the dataset.')
parser.add_argument('-n', '--n_max_jobs', default=-1, type=int,
                    help='Maximum number of jobs.')
parser.add_argument('-g', '--geo', default="ALLEGRO_o1_v03", type=str,
                    help='Detector geometry.')


def run(args):

    tag = args.tag
    input_file_path = args.input
    output_file_path = args.output
    n_max = args.n_max_jobs
    geo = args.geo

    detector = geo.split("_")[0]

    storage_path_parent = os.path.join(input_file_path, output_file_path)
    storage_path = os.path.join(storage_path_parent, tag)

    print("Creating output storage path:")
    print(storage_path)
    os.makedirs(storage_path, exist_ok=True)

    print("Creating submission folder:", tag)
    if not os.path.isdir(tag):
        os.mkdir(tag)

    # Setup the bash executables scripts
    print("Preparing submission for:")
    n_jobs = 0
    exec_template_name = os.path.join(tag, "run_ddsim_FILENAME.sh")

    for folder in os.listdir(input_file_path):
        if n_max > 0 and n_max <= n_jobs:
            break
        if "data" not in folder:
            continue

        bx_id = folder.replace("data", "")
        input_filename = os.path.join(input_file_path, folder, "pairs.pairs")
        print(input_filename)
        executable_path = exec_template_name.replace("FILENAME", bx_id)
        output_filename = os.path.join(storage_path, f"{geo}_{bx_id}.root")

        # for performance, write the output locally first and copy at the end
        tmp_output_filename = os.path.basename(output_filename)
        command = f"""ddsim \
            --compactFile $K4GEO/FCCee/{detector}/compact/{geo}/{geo}.xml \
            -I {input_filename} \
            -O {tmp_output_filename} \
            -N -1 --crossingAngleBoost 0.015 \
            --part.keepAllParticles True"""

        command += f"\nmv {tmp_output_filename} {output_filename}"
        with open(executable_path, "w") as f:
            f.write(exec_header)
            f.write(command)
        st = os.stat(executable_path)
        os.chmod(executable_path, st.st_mode | stat.S_IEXEC)
        n_jobs += 1

    # Setup the condor script
    condor_submit_path = f"{tag}.cmd"
    exec_pattern = exec_template_name.replace("FILENAME", "*")
    cmd_file_content = condor_cmd_content.format(exec_pattern)

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
