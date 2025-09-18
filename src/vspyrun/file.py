import sys
import numpy as np
import os
import shutil

def get_file_name(string, directory):
    """
    Returns a list of names of all regular files in `directory` that end with ".vasp".
    By default, `directory="."` means the current folder.
    """
    files = []
    for item in os.listdir(directory):
        path = os.path.join(directory, item)
        if os.path.isfile(path) and item.endswith(f".{string}"):
            files.append(item)
    return files

def write_file(path: str, text: str) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write(text if text.endswith("\n") else text + "\n")

def prep_slurm(template, jname, outfile):
    shutil.copy(template, outfile)
    with open(outfile, 'r') as file:
        filedata = file.read()
    filedata = filedata.replace('sample_r1', f'{jname}')
    with open(outfile, 'w') as file:
        file.write(filedata)
