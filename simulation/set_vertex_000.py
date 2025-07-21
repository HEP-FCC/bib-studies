#!/usr/bin/env python

"""
  Reset the position of particles to (0,0,0) in `.pairs` files
  created by GuineaPig. This is required as the generator doesn't
  include any B-field. Therefore, the positions are inexact.
  This adjustment is also inexact, but more realistic.

  Original script from Brieuc:
  https://github.com/BrieucF/FCC_scripts/blob/main/background_studies/set_vertex_000.py

  Content parsed by Geant4EventReaderGuineaPig:
  https://github.com/AIDASoft/DD4hep/blob/bb0dc9417100af0af2d25465182ca9f175f42531/DDG4/plugins/Geant4EventReaderGuineaPig.cpp#L211-L213C41
"""

import os
import glob
import argparse

# Path to Z folder from Andrea Ciarma's IPC files
def_path = "/eos/experiment/fcc/users/a/aciarma/pairs/4IP_2024may29/Z/data*"

# Argument parser
parser = argparse.ArgumentParser('Set particle vertex to origin.')
parser.add_argument('-i', '--input', default=def_path,
                    help='input path expression. Default is: '+def_path)
parser.add_argument('-o', '--output', required=True,
                    help='output folder path.')
parser.add_argument('-d', '--do_dat', action='store_true',
                    help='Check also for ".dat" files.')


def run(args):
    input_path = args.input
    output_path = args.output
    do_dat = args.do_dat

    for data_folder_path in glob.glob(input_path):
        bx_id = data_folder_path.split("/")[-1].replace("data", "")
        input_file = os.path.join(data_folder_path, "pairs.pairs")

        if not os.path.exists(input_file):
            if do_dat:
                input_file = input_file.replace(".pairs", ".dat")
                if not os.path.exists(input_file):
                    print(f"Couldn't find .pairs or .dat file for {bx_id}...")
                    continue
                else:
                    print(f"Processing event {bx_id} (from .dat)...")
            else:
                print(f"Skipping {bx_id}, file not produced previously")
                continue
        else:
            print(f"Processing event {bx_id}...")

        new_file_content = ""
        with open(input_file) as i_file:
            for line in i_file.readlines():
                entries = line.split(" ")
                # set the vtx to 000
                entries[4] = 0
                entries[5] = 0
                entries[6] = 0
                for entry in entries:
                    new_file_content += str(entry) + " "
            new_file_content += "\n"

        output_folder = os.path.join(output_path, f"data{bx_id}")
        os.makedirs(output_folder, exist_ok=True)
        with open(os.path.join(output_folder, "pairs.pairs"), "w") as o_file:
            o_file.write(new_file_content)


if __name__ == "__main__":
    args = parser.parse_args()
    run(args)
