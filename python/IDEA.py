"""
    IDEA specialized functions.
"""

from collections import defaultdict
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
    sensor_size_map = defaultdict(float)
    sensors_per_module_map = defaultdict(int)


    # N.B. this could be affected from naming scheme changes,
    # will need to keep track also of the geo versions.
    match name:

        case "VertexBarrel" | "SiWrB":
            # Rename layers shifting their value by 1
            # to remove degeneracy of layer 0
            # N.B. this needs to be accounted when reading the layer number

            for de_name, de in sub_det.children():
                modules = 0
                sensors = 0
                for de_name2, de2 in de.children():
                    modules += 1
                    for de_name3, de3 in de2.children():
                        sensors += 1
                        area = de3.volume().solid().GetDY()*2.*10.*de3.volume().solid().GetDZ()*2.*10. # Make it to mm and get sensor area in mm2. Assuming each sensor has the same area!
                        # print("   sensor area (mm2): ", area)
                cells_map[str(de_name)] = sensors
                sensor_size_map[str(de_name)] = area
                sensors_per_module_map[str(de_name)] = int(sensors / modules)

        case "VertexDisks" | "SiWrD":
            # Rename layers shifting their value by 1
            # to remove degeneracy of layer 0
            # N.B. this needs to be accounted when reading the layer number

            for de_name, de in sub_det.children():
                modules = 0
                sensors = 0
                for de_name2, de2 in de.children():
                    modules += 1
                    for de_name3, de3 in de2.children():
                        sensors += 1
                        area = de3.volume().solid().GetDX()*2.*10.*de3.volume().solid().GetDY()*2.*10. # Make it to mm and get sensor area in mm2. Assuming each sensor has the same area!
                        # print("   sensor area (mm2): ", area)
                cells_map[str(de_name)] = sensors
                sensor_size_map[str(de_name)] = area
                sensors_per_module_map[str(de_name)] = int(sensors / modules)

            old_keys = list(cells_map.keys())
            for k in old_keys:
                ln = re.search("layer[0-9]", k).group(0) 
                ln = int(ln.strip("layer")) + 1
                if "side-1" in k:
                    ln *= -1
                new_key = f"layer{ln}"
                cells_map[new_key] = cells_map.pop(k)
                sensor_size_map[new_key] = sensor_size_map.pop(k)
                sensors_per_module_map[new_key] = sensors_per_module_map.pop(k)

        case _:
            # Loop over detector elements
            for de_name, de in sub_det.children():
                # layer_id = de.id()
                # print(f"         layer_name (id): {de_name} ({layer_id})")
                if re_skip.match(str(de_name)):
                        print("Skipping sub detector:", de_name)
                        continue
                cells_map[str(de_name)] = get_cells(de)

    return cells_map, sensor_size_map, sensors_per_module_map


def get_layer(cell_id, decoder, detector, dtype):
    """
    Run the decoder differently for each sub-detector
    """

    match detector:

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

def get_module(cell_id, decoder, detector, dtype):
    try:
        module = decoder.get(cell_id, "module")
    except:
        module = 0 # For systems that don't have modules
    return module

def get_sensor(cell_id, decoder, detector, dtype):
    sensor = decoder.get(cell_id, "sensor")
    return sensor    