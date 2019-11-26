"""
some functions that will convert binary images to pngs

Jacob Long
2019-11-14
"""

import fnmatch
import multiprocessing
import os

import PIL
from PIL import Image

def bin_to_png(filepath):
    """converts a single .bin image to png"""
    # written by David Moller
    try:
        raw_data = open(filepath, 'rb').read()
        img_size = (1936, 1216)
        img = PIL.Image.frombytes('L', img_size, raw_data, 'raw')
        savename = filepath[:-4] + '.png'
        img.save(savename)
        print(savename, "converted")
    except Exception as e:
        print(f"{filepath}: {e}")

def convert_dir(directory):
    """converts all the .bin images in a directory to .png"""

    for f in os.listdir(directory):
        if f.endswith(".bin"):
            bin_to_png(os.path.join(directory, f))

def convert_dirs(directory, processes=-1):
    """calls convert dir in a new process on each
    of the subdirectories in a given directory"""

    # finding all directories in the given dir
    dirs_to_process = []
    for root, dirs, files, in os.walk(directory):
        if files:
            dirs_to_process.append(root)

    # using the max number of processes if none are specified
    if processes < 1:
        processes = multiprocessing.cpu_count()

    # opening a new process for each directory found
    with multiprocessing.Pool(multiprocessing.cpu_count()) as p:
        p.map(convert_dir, dirs_to_process)

    print("Finished converting to pngs in", directory)

if __name__ == "__main__":

    filepath = r"W:\GeneratedData\F013B2\PS2\Test"

    convert_dirs(filepath)
    