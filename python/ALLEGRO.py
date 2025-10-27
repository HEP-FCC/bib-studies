"""
    ALLEGRO specialized functions.
"""

from collections import defaultdict
import math
import re

from helpers import get_cells, is_endcap

def get_cells_map(detector, sub_det, name, skip_pattern = r"(supportTube)|(cryo)"):
    """
    Get the mapping of cells per layer in all sub-detectors.
     Args:
        detector: detector object (from the XML)
        sub_det:  sub-detector object
        name:     name of the sub-detector
     Parameters:
      skip_pattern: regex pattern to skip sub-detectors
    """

    re_skip = re.compile(skip_pattern)

    cells_map = defaultdict(int)

    # N.B. this could be affected from naming scheme changes,
    # will need to keep track also of the geo versions.
    match name:

        case "EMEC_turbine":
            # TODO: Fix this part, this should use detector segmentation instead!
            # Method 1, parse the constants from the geo file
            n_wheels = detector.constantAsLong("nWheels")

            for nw in range(n_wheels):
                # Uses "Calib" instead of ReadOut layers (need to check with a detector expert)
                r_layers = detector.constantAsLong(f"ECalEndcapNumReadoutRhoLayersWheel{nw+1}")
                z_layers = detector.constantAsLong(f"ECalEndcapNumReadoutZLayersWheel{nw+1}")
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
            merged_phi_planes = segmentation.parameter("mergedModules").value()        # List of number of merged segments (planes) along phi

            n_avg_theta_cells = (((math.pi/ 2) - theta) * 2) / d_theta  # number of cells along theta per layer, on average

            merged_theta_cells = merged_theta_cells.split(' ')  # numbers are encoded in as single string, separated with whites spaces
            merged_phi_planes = merged_phi_planes.split(' ')    # same as above

            tot_cells = 0
            # Calculate manually the number of cells per layer
            for l in range(n_layers):
                m_theta = int(merged_theta_cells[l])
                m_phi = int(merged_phi_planes[l])

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

        case "SiWrD":
            #loop from -2 to 2, skipping 0, to set the layers indices for the endcaps
            for i in range(-2,3):
                if i==0: continue
                cells_map[f"layer{i}"] = 7500  #hardcoding channels per layer for now

        case "SiWrB":
            # loop from 0 to 1 to set the layers indices for the barrel
            for i in range(2):
                cells_map[f"layer{i}"] = 20_000  #hardcoding channels per layer for now, https://arxiv.org/abs/2502.21223v4

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


def get_layer(cell_id, decoder, detector, dtype):
    """
    Run the decoder differently for each sub-detector
    """

    match detector:
        case "HCalThreePartsEndcap":
            # Get only the side of the hit
            # type = 0,1,2 for positive side
            # type = 3,4,5 for negative side
            tp =  decoder.get(cell_id,"type") #FIXME: Appears to be always 0?! double check using nightly sim. events
            #print(f"layer={layer}, side={side}, type={tp}")
            if tp > 2:
                return -1
            return 1

        case "EMEC_turbine":
            side = decoder.get(cell_id, "side")
            wheel = decoder.get(cell_id, "wheel") + 1
            return  wheel * side

        case "DCH_v2":
            # Number of layers per super layer could be read from geo file
            nl_x_sl = 8
            layer = decoder.get(cell_id, "layer")
            super_layer = decoder.get(cell_id, "superlayer")
            return (super_layer * nl_x_sl) + layer + 1

        case "VertexDisks" | "SiWrD":
            # shift layer number by 1 to remove degeneracy of layer 0
            layer = decoder.get(cell_id, "layer") + 1
            side = decoder.get(cell_id, "side")

            return layer * side

        case "MuonTaggerEndcap":
            # Default way: side * layer, where side should be +/- 1
            layer = decoder.get(cell_id, "layer") + 1
            #probably no side available in decoder (see xml readout part)
            theta = decoder.get(cell_id, "theta")
            #print(f"layer={layer}, theta={theta}")

            side=-1
            if theta<168:
              side = 1

            return layer * side

        case _:
            layer = decoder.get(cell_id, "layer")

            # Get the side, if available
            side = 0
            if is_endcap(dtype):
                side = decoder.get(cell_id, "side")

            if side != 0:
                layer *= side

            return layer
    return
