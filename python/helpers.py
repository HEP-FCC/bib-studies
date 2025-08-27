import glob
import os


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

