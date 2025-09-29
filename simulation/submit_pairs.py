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
import re
import stat

# Default paths / namings
input_def_path = "/eos/home-s/sfranche/FCC/samples/bib/gen-samples/aciarma_4IP_2024may29/Z/"
output_def_folder = input_def_path+"/DDSim_output/"


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
parser.add_argument('-c', '--compactFile', default="$K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml", type=str,
                    help='Detector geometry.')
parser.add_argument('-k', '--k4geo', default=None, type=str,
                    help='Path to custom k4geo.')
parser.add_argument('--crossingAngleBoost', default=0.015, type=str,
                    help='Crossing angle boost to be applied.')


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

# Header of executable script
fcc_cfg = os.environ["FCCCONFIG"]
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

# Path to steering files
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
    if geo in steering_dict:
        steering_path = steering_dict[geo]
        print("Including steering file: ", steering_path)
        steering_opt = f"--steeringFile {steering_path}"

    # Setup the bash executables scripts
    print("Preparing submission for:")
    n_jobs = 0
    exec_template_name = os.path.join(tag, "run_ddsim_FILENAME.sh")

    for folder in os.listdir(input_file_path):
        if (n_max > 0) and (n_max <= n_jobs):
            break
        if "data" not in folder:
            continue

        command = header
        bx_id = folder.replace("data", "")
        input_filename = os.path.join(input_file_path, folder, "pairs.pairs")
        print(input_filename)
        executable_path = exec_template_name.replace("FILENAME", bx_id)
        output_filename = os.path.join(storage_path, f"{geo}_{bx_id}_r{fcc_ver}.root")

        # for performance, write the output locally first and copy at the end
        tmp_output_filename = os.path.basename(output_filename)

        command += f"""ddsim \
            --compactFile  {compact} \
            -I {input_filename} \
            -O {tmp_output_filename} \
            -N -1 --crossingAngleBoost {x_angle} \
            --part.keepAllParticles True  {steering_opt}\n"""

        command += f"mv {tmp_output_filename} {output_filename}"

        with open(executable_path, "w") as f:
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
