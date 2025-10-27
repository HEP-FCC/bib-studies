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

    # N.B. this could be affected from naming scheme changes,
    # will need to keep track also of the geo versions.
    match name:

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
            # loop from -2 to 2, skipping 0, to set the layers indices for the endcaps
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
