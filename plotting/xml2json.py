#!/usr/bin/env python

from collections import defaultdict
import json
import math
from optparse import OptionParser
import re

import dd4hep as dd4hepModule
from ROOT import dd4hep


######################################
# option parser
parser = OptionParser()
parser.add_option('-d', '--detGeoFile',
                  type=str, default='$K4GEO/FCCee/ALLEGRO/compact/ALLEGRO_o1_v03/ALLEGRO_o1_v03.xml',
                  help='compact XML detector geometry file')

(options, args) = parser.parse_args()

geo_file_name = options.detGeoFile


######################################
# Functions

skip_pattern = r"(supportTube)|(cryo)"
re_skip = re.compile(skip_pattern)

def get_cells(detector, n_cells = 0):
    sub_detectors = detector.children()
    if sub_detectors.size() == 0:
        # print("Counting", detector.GetName())
        # print("Counting", detector.id())
        n_cells += 1
    else:
        for d in sub_detectors:
            d_name = str(d[0])
            if re_skip.match(d_name):
                print("Skipping sub detector:", d_name)
                continue
            n_cells = get_cells(d[1], n_cells)
    return n_cells


def get_cells_map(detector, sub_det, name):

    cells_map = defaultdict(int)

    # N.B. this could be affected from naming scheme changes,
    # will need to keep track also of the geo versions.
    match name:
        case "EMEC_turbine":
            # Method 1, parse the constants from the geo file
            n_wheels = detector.constantAsLong("nWheels")

            for nw in range(n_wheels):
                # Uses "Calib" instead of ReadOut layers (need to check with a detector expert)
                r_layers = detector.constantAsLong(f"ECalEndcapNumCalibRhoLayersWheel{nw+1}")
                z_layers = detector.constantAsLong(f"ECalEndcapNumCalibZLayersWheel{nw+1}")
                n_units = detector.constantAsLong(f"nUnitCells{nw+1}")

                cells_map[f"wheel{nw+1}"] = r_layers * z_layers * n_units
                cells_map[f"wheel-{nw+1}"] = r_layers * z_layers * n_units

        case "HCalThreePartsEndcap":
            # Split only by positive and negative layers
            for de_name, de in sub_det.children():
                layer_name = str(de_name)
                if "layer-" in layer_name:
                    cells_map["layer_-1"] += get_cells(de)
                elif "layer" in layer_name:
                    cells_map["layer_1"] += get_cells(de)

        case "ECalBarrel":
            # The ECal barrel is defined as a ECalBarrel_NobleLiquid_InclinedTrapezoids_o1_v03 object
            # and has a readout implemented as a geometrical grid of the type FCCSWGridModuleThetaMerged_k4geo.
            # This (probably) means that hits happen in active material and then are associated to a cellID
            # that is defined from the grid object, without any sub-sub-...-sub-detector element associated.
            # Hence, the number of cells must be calculated from the detector and the grid parameters.

            n_layers = detector.constantAsLong("ECalBarrelNumLayers")  # Number of layers
            n_planes = detector.constantAsLong("ECalBarrelNumPlanes")  # Planes used to segment the phi angle
            r_min = detector.constantAsDouble("Bath_rmin")             # Inner radius of the active material
            r_max = detector.constantAsDouble("Bath_rmax")             # Outer radius of the active material
            dz = detector.constantAsDouble("EMBarrel_dz")              # Half length of the barrel

            # Consider only avg radius as approximation
            theta = math.atan( (r_min + r_max) * 0.5 / dz)

            # Retrieve the segmentation object
            segmentation = detector.sensitiveDetector(name).readout().segmentation().segmentation()
            d_theta = float(segmentation.parameter("grid_size_theta").value())         # Grid cell size along theta
            merged_theta_cells = segmentation.parameter("mergedCells_Theta").value()   # List of number of merged cells per layer
            phi_merged_planes = segmentation.parameter("mergedModules").value()        # List of number of merged segments (planes) along phi

            n_avg_theta_cells = (((math.pi/ 2) - theta) * 2) / d_theta  # number of cells along theta per layer, on average

            merged_theta_cells = merged_theta_cells.split(' ')  # numbers are encoded in as single string, separated with whites spaces
            phi_merged_planes = phi_merged_planes.split(' ')    # same as above

            tot_cells = 0
            # Calculate manually the number of cells per layer
            for l in range(n_layers):
                m_theta = int(merged_theta_cells[l])
                m_phi = int(phi_merged_planes[l])

                n_theta = int(n_avg_theta_cells / m_theta)
                n_phi = int(n_planes / m_phi)

                cells_map[f"layer{l}"] = n_theta * n_phi
                tot_cells += n_theta * n_phi
            print("         - Total cells calculated:",tot_cells)

        case "MuonTaggerBarrel":
            #from XML file in 2025:
            # theta = 336 bins
            # phi   = 704 bins
            # => full muon system has 336*704 cells
            #todo: FIX THIS APPROXIMATION:
            # theta binning not yet clear from XML, so very roughly assigning:
            total_cells = 2/3 * 336 * 704


            n_layers = detector.constantAsLong(f"MuonTaggerBarrelLayers")
            cells_map =  {f"layer{i}": total_cells for i in range(n_layers)}

        case "MuonTaggerEndcap":
            #todo: need to clarify with MuonTagger experts if this is correct (see also barrel above)
            #todo: FIX THIS APPROXIMATION:
            # theta binning not yet clear from XML, so very roughly assigning:
            total_cells = 1/3 * 336 * 704
            n_layers = detector.constantAsLong(f"MuonTaggerEndcapLayers")
            cells_map =  {}
            #todo: fix this fix
            #loop from -2 to 2, skipping 0, to set the layers indices for the endcaps
            for i in range(-2,3):
                if i==0: continue
                cells_map[f"layer{i}"] = total_cells

        case "VertexDisks":
            # Rename layers shifting their value by 1
            # to remove degeneracy of layer 0
            # N.B. this needs to be accounted when reading the layer number

            for de_name, de in sub_det.children():
                cells_map[str(de_name)] = get_cells(de)

            old_keys = list(cells_map.keys())
            for k in old_keys:
                ln = re.search("layer[0-9]", k).group(0) 
                ln = int(ln.strip("layer")) + 1
                if "side-1" in k:
                    ln *= -1
                new_key = f"layer{ln}"
                cells_map[new_key] = cells_map.pop(k)

        case _:
            # Loop over detector elements
            for de_name, de in sub_det.children():
                # layer_id = de.id()
                # print(f"         layer_name (id): {de_name} ({layer_id})")
                if re_skip.match(str(de_name)):
                        print("Skipping sub detector:", de_name)
                        continue
                cells_map[str(de_name)] = get_cells(de)
    return cells_map


##################################################
# Geometry part

# useful doxy:
# https://dd4hep.web.cern.ch/dd4hep/reference/classdd4hep_1_1Detector.html
# https://dd4hep.web.cern.ch/dd4hep/reference/classdd4hep_1_1SensitiveDetector.html
# https://dd4hep.web.cern.ch/dd4hep/reference/classdd4hep_1_1DetElement.html
# https://dd4hep.web.cern.ch/dd4hep/reference/classdd4hep_1_1Volume.html
# https://dd4hep.web.cern.ch/dd4hep/reference/classdd4hep_1_1Box.html
# https://dd4hep.web.cern.ch/dd4hep/reference/DD4hep_2Detector_8h_source.html

# Get detector name
detector = geo_file_name.split("/")[-1].strip(".xml")
# Load geometry
det = dd4hep.Detector.getInstance()
det.fromXML(geo_file_name)
# Get subdetectors
dict_sub = {}

for subdet_name, sens_det in det.sensitiveDetectors():
    print("==== name = ", subdet_name) 

    # Get the DetElement object
    subdet = det.detector(subdet_name)
    print("        ID = ", subdet.id())
    print("        type = ", subdet.type())
    print("        typeFlag = ", subdet.typeFlag())

    hitsCollection = str(sens_det.hitsCollection)
    print("        hits = ", hitsCollection)

    det_element_cells = get_cells_map(det, subdet, subdet_name)
    print("        det_element_cells = ", det_element_cells)

    # Get max dimensions for the subdetector in x,y,z (half-length form 0)
    volume = subdet.volume()
    box = volume.boundingBox()

    # dimensions seem off (ask Brieuc) - for the moment multiply the max by a factor 10
    max_z = 10*math.ceil(box.z()) #max z rounded up
    max_r = 10*math.ceil(math.sqrt(math.pow(box.x(),2)+math.pow(box.y(),2))) #max r rounded up
    print("        max_z = ", max_z)
    print("        max_r = ", max_r)
    if(box.x()!=box.y()): print ('WARNING: different X/Y in the detector bounding box X:',box.x(),'Y:',box.y())

    if subdet.id()>0.:
        dict_sub[str(subdet_name)]={
            'id': int(subdet.id()),
            'typeFlag': subdet.typeFlag(),
            'hitsCollection': hitsCollection,
            'det_element_cells': det_element_cells,
            'max_z': max_z, 
            'max_r': max_r,
        }

with open(geo_file_name.split('/')[-1].strip('.xml')+'_DetectorDimensions.json', 'w') as fp:
    json.dump(dict_sub, fp, indent=2)
