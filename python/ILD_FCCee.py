"""
    ILD at FCCee specialized functions.
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

        case _:
            # Loop over detector elements
            for de_name, de in sub_det.children():

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

