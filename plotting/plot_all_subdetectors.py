#!/usr/bin/env python

"""
Run drawhits.py (and optionally hits2highLevelEstimations.py) for all the sub-detectors of an assumptions file,
several sub-detectors in parallel. One log file per sub-detector, and a summary table at the end.

Options not known to this script are forwarded to drawhits.py, e.g. -p -m --integration_time 10.
"""

import argparse
import glob
import json
import os
import re
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

from helpers import print_header, print_warning


######################################
# option parser

parser = argparse.ArgumentParser(description='plot_all_subdetectors.py',
        epilog='Example:\nplot_all_subdetectors.py -i <sim_file.root> -e 100 -j 4 --hits2high --skip_layers_for DCH_v2 -p -m',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--inputFile', required=True,
                    help='input file or directory, simulated with the **same** geometry as the detector dictionary')
parser.add_argument('-d', '--detDictFile', default='$BIB_STUDIES/detectors_dicts/ALLEGRO_o1_v03_DetectorDimensions.json',
                    help='JSON dictionary with the key detector parameters (see xml2json.py)')
parser.add_argument('-a', '--assumptions', default='$BIB_STUDIES/detectors_dicts/ALLEGRO_o1_v03_assumptions.json',
                    help='JSON dictionary with the assumptions, its keys define the sub-detectors to run over')
parser.add_argument('-n', '--numberOfFiles', type=int, default=None,
                    help='number of input files, passed to drawhits.py (its default if not set)')
parser.add_argument('-e', '--numberOfEvents', type=int, default=None,
                    help='number of events per file, passed to drawhits.py (its default if not set)')
parser.add_argument('--sample', default='ipc',
                    help='sample name, used as prefix of the output files')
parser.add_argument('-D', '--debugLevel', type=int, default=0,
                    help='debug level of drawhits.py (0:quiet, 1:default, 2:verbose)')
parser.add_argument('-j', '--jobs', type=int, default=4,
                    help='number of sub-detectors processed in parallel')
parser.add_argument('--only', default='',
                    help='comma-separated list of sub-detectors, to run only on these')
parser.add_argument('--exclude', default='',
                    help='comma-separated list of sub-detectors to skip')
parser.add_argument('--skip_layers_for', default='',
                    help='comma-separated list of sub-detectors for which to add --skip_layers (e.g. DCH_v2, too many layers)')
parser.add_argument('--hits2high', action='store_true',
                    help='also run hits2highLevelEstimations.py on the output of each sub-detector')
parser.add_argument('--logdir', default='.',
                    help='directory for the log files')
parser.add_argument('--dry_run', action='store_true',
                    help='only print the commands')

options, drawhits_extra_args = parser.parse_known_args()

split = lambda s: [x.strip() for x in s.split(',') if x.strip()]
det_dict_path = os.path.expandvars(options.detDictFile)
assumptions_path = os.path.expandvars(options.assumptions)
sample_name = options.sample
os.makedirs(options.logdir, exist_ok=True)


######################################
# sub-detectors to run over: the ones of the assumptions file, that are also in the detector dictionary

with open(assumptions_path) as f:
    assumptions_subdets = list(json.load(f).keys())
with open(det_dict_path) as f:
    dict_subdets = set(json.load(f).keys())

subdetectors = []
for sd in assumptions_subdets:
    if sd not in dict_subdets:
        print_warning(f"{sd} is in the assumptions but not in {os.path.basename(det_dict_path)}, skipped")
        continue
    if options.only and sd not in split(options.only):
        continue
    if sd in split(options.exclude):
        continue
    subdetectors.append(sd)

if not subdetectors:
    raise SystemExit("ERROR: no sub-detector to run over")


######################################
# commands

def drawhits_command(sd):
    cmd = ["drawhits.py",
           "--infilePath", options.inputFile,
           "--subDetector", sd,
           "--detDictFile", det_dict_path,
           "--assumptions", assumptions_path,
           "--sample", sample_name,
           "--debugLevel", str(options.debugLevel)]
    if options.numberOfFiles is not None:
        cmd += ["--numberOfFiles", str(options.numberOfFiles)]
    if options.numberOfEvents is not None:
        cmd += ["--numberOfEvents", str(options.numberOfEvents)]
    if sd in split(options.skip_layers_for):
        cmd += ["--skip_layers"]
    return cmd + drawhits_extra_args


def find_drawhits_output(sd, since):
    """Newest drawhits output file of this sub-detector, written after 'since' (hits2highLevelEstimations outputs excluded)"""
    files = [f for f in glob.glob(f"{sample_name}_*evt_{sd}_*.root")
             if "high-level" not in f and os.path.getmtime(f) >= since - 1]
    return max(files, key=os.path.getmtime) if files else None


def run_subdetector(sd):
    """Run drawhits.py (and hits2highLevelEstimations.py) for one sub-detector, return a summary dictionary"""
    result = {"sub_detector": sd, "ok": False, "time": 0., "bandwidth": None, "log": os.path.join(options.logdir, f"log_{sample_name}_{sd}.txt")}
    start = time.time()
    with open(result["log"], "w") as log:
        rc = subprocess.run(drawhits_command(sd), stdout=log, stderr=subprocess.STDOUT).returncode
        result["ok"] = (rc == 0)

        if rc == 0 and options.hits2high:
            root_file = find_drawhits_output(sd, start)
            if root_file is None:
                log.write("\nERROR: drawhits output file not found, hits2highLevelEstimations.py not run\n")
                result["ok"] = False
            else:
                log.write(f"\n########## hits2highLevelEstimations.py on {root_file}\n")
                log.flush()
                h2h = subprocess.run(["hits2highLevelEstimations.py", "-i", root_file, "-d", det_dict_path, "-a", assumptions_path],
                                     stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
                log.write(h2h.stdout)
                result["ok"] = (h2h.returncode == 0)
                match = re.search(r"Total bandwidth = ([0-9.eE+-]+) GB/s", h2h.stdout)
                if match:
                    result["bandwidth"] = float(match.group(1))
    result["time"] = time.time() - start
    return result


######################################
# run

if options.dry_run:
    for sd in subdetectors:
        print(" ".join(drawhits_command(sd)))
    raise SystemExit(0)

print_header(f"Running over {len(subdetectors)} sub-detectors, {options.jobs} in parallel", oneline=True)
print("Sub-detectors:", ", ".join(subdetectors))

results = []
with ThreadPoolExecutor(max_workers=max(1, options.jobs)) as pool:
    futures = {pool.submit(run_subdetector, sd): sd for sd in subdetectors}
    for future in as_completed(futures):
        r = future.result()
        results.append(r)
        print(f"  {'done  ' if r['ok'] else 'FAILED'} {r['sub_detector']:<24s} {r['time']:7.1f} s   ({len(results)}/{len(subdetectors)})")


######################################
# summary

print_header("Summary", oneline=True)
bandwidth_col = options.hits2high
print(f"{'sub-detector':<24s} {'status':<7s} {'time [s]':>9s}" + (f" {'bandwidth [GB/s]':>17s}" if bandwidth_col else "") + "  log")
for r in sorted(results, key=lambda r: subdetectors.index(r["sub_detector"])):
    bw = ""
    if bandwidth_col:
        bw = f" {r['bandwidth']:>17.2f}" if r["bandwidth"] is not None else f" {'-':>17s}"
    print(f"{r['sub_detector']:<24s} {'OK' if r['ok'] else 'FAILED':<7s} {r['time']:>9.1f}{bw}  {r['log']}")

failed = [r for r in results if not r["ok"]]
for r in failed:
    print_warning(f"{r['sub_detector']} failed, last lines of {r['log']}:")
    with open(r["log"]) as f:
        print("".join("    " + line for line in f.readlines()[-5:]))
raise SystemExit(1 if failed else 0)
