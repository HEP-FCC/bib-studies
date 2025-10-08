import glob
import json
import os
import re


def path_to_list(path, ext=".root"):
    """
    Parse the input path and return a sorted list of ROOT files only. 

    We can read:
      1 - path to a ROOT file
      2 - path to a directory containing ROOT files
      3 - path with wildcard to ROOT files. eg: "path/to/dir*/*.root", 
          N.B. the quotation marks ("") are required when passing such path in the CLI
    """

    input_path = glob.glob(path)

    list_files = []
    if os.path.isdir(input_path[0]):
        for p in input_path:
            list_files += glob.glob(p+"/*"+ext)
    elif isinstance(input_path, str):
        list_files = [input_path]
    else:
        list_files = input_path

    return list_files


def sorted_n_files(list_files, n_files=-1, ext=".root"):
  """
  list n_files in the path directory ending in .root
  use sorted list to make sure to always take the same files
  """
  
  list_files = [f  for f in sorted(list_files) if f.endswith(ext)]

  if n_files > 0: #if n_files=-1 run all files
      if n_files <= len(list_files): 
          list_files = list_files[:n_files]
      else:
          raise ValueError(f"Asked to process {n_files} input files,"+
                          f"but only {len(list_files)} where found. Check your settings!")
  return list_files


def load_json(path, sub_dict=None):
    """
    Load JSON dictionary. If a "sub_dict" specified,
    only that one will be loaded, otherwise the entire JSON is passed.
    """
    out_dict = {}
    
    with open(path,"r") as f:
        json_dict = json.load(f)

        if sub_dict==None:
            return json_dict

        try:
            out_dict = json_dict[sub_dict]
        except KeyError:
            raise KeyError(f"'{sub_dict}' is not available, valid sub-dictionary entries are: "
                        +" | ".join(json_dict.keys()))
    
    return out_dict


def layer_number_from_string(layer_string):
    """
    Extract layer number form the layer string.
    """

    # check if the it is endcap / wheels / disk with any side defined
    # - should be +1 or -1
    side = 0
    try:
        side_name = re.search("(side)(_)?(-)?[0-9]+", layer_string).group(0)
        side = int(re.sub(r"[^0-9-]","",side_name))
        if abs(side) > 1:
            raise NotImplementedError(f"Found side = {side} for '{layer_string}', not supported")
    except AttributeError:
        pass

    # get the layer name
    try:
        layer_name = re.search("(layer|wheel)(_)?(-)?[0-9]+", layer_string).group(0) 
    except AttributeError:
        raise AttributeError(f"Warning: Couldn't read layer number from {layer_string}")

    # Get the layer number, sign based on the side
    ln = int(re.sub(r"[^0-9-]","",layer_name))
    if side != 0:
        ln *= side
    
    return ln
