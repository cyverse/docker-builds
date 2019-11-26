"""
(formerly known as python_macro.py)
multi threshold image segmentation

A replacement for python_conv_WOPR_08-17-2019.ijm so 
that the whole process is in python.

1. find all folders in a PS2 collection
2. find all images in that folder
3. for each image, look at a 1/3 center horizontal slice
4. apply thresholds and count the amount of pixels that fall in them
5. find the average value inside each threshold
6. output a csv with these values

Jacob Long 
2019-11-1
"""

import multiprocessing
import os
import time
from collections import Counter, defaultdict
from pprint import pprint

import numpy as np
import pandas as pd
from PIL import Image


def apply_threshold(filepath):
    """
    applying a 5 thresholds to a single image.
    returns a list of 5 dictionaries.
    [
        {
            'Label': image name,
            'Area':  how many pixels met this threshold,
            'Mean':  average pixel value of pixels that met the threshold,
            'Min':   min of threshold,
            'Max':   max of threshold,
        },
        {
            ...
        }
    ]

    filepath: a filepath to a png
    """

    # reading image
    im = Image.open(filepath)
    img_array = np.asarray(im)

    # getting middle 1/3 horizontal slice of image
    width, height = im.size
    crop_line = height // 3
    img_array_cropped = img_array[crop_line:crop_line*2]

    # img = Image.fromarray(img_array_cropped)
    # img.show()

    all_pixels = []
    for row in img_array_cropped:
        all_pixels.extend(row)
    all_pixels = np.asarray(all_pixels)

    # defining thresholds
    thresholds = {
        "t1": (0,  7  ),
        "t2": (8,  10 ),
        "t3": (11, 14 ),
        "t4": (15, 19 ),
        "t5": (20, 255),
    }

    # finding how many times each pixel brightness occured in the image
    unique, counts = np.unique(all_pixels, return_counts=True)
    pixel_counts = dict(zip(unique, counts))

    # for each threshold, add the counts to an accumulator dict
    threshold_counts = Counter()
    averages_dict = {}

    for threshold_name, threshold in thresholds.items():

        average_dict = {}

        for pixel_value, count in pixel_counts.items():
            if threshold[0] <= pixel_value <= threshold[1]:

                # adding to the amount of counted pixel values that are counted in this threshold
                threshold_counts[threshold_name] += count

                average_dict[pixel_value] = count
        
        # finding the weighted average of each pixel value
        num = 0
        total = 0
        for pixel_value, count in average_dict.items():

            num += pixel_value * count
            total += count
            # print(pixel_value, count, total)

        # print("total:", total)
        try:
            averages_dict[threshold_name] = num/total
        except:
            averages_dict[threshold_name] = 0
        

    # formatting each threshold into a dictionary 
    # row to be formatted into a csv later
    output_list = []

    for key, threshold in thresholds.items():
        output_list.append({
            'Label': os.path.basename(filepath),
            'Area':  threshold_counts[key],
            'Mean':  averages_dict[key],
            'Min':   threshold[0],
            'Max':   threshold[1],
        })

    return output_list

def process_dir(filepath):
    """
    finds all pngs in a directory and outputs a 
    csv of thresholded images
    """

    # finding all png files in the subdirectory
    image_paths = [
        os.path.join(filepath, f) 
        for f in os.listdir(filepath) 
        if f.endswith(".png")
    ]

    image_dicts = []
    for ip in image_paths:
        image_dicts.extend(apply_threshold(ip))

    df = pd.DataFrame(image_dicts)


    output_dir = os.path.join(filepath, os.path.basename(filepath) + ".csv")
    df.to_csv(output_dir)
    print(output_dir, "created")

def process_collection(filepath):
    """
    finds all directories in a ps2 collection 
    and runs process_dir() on each of them
    using multiprocessing
    """

    # finding all directories to process
    for root, dirs, files in os.walk(filepath):

        # getting each directory in the main folder
        directories = [os.path.join(root, d) for d in dirs]

        with multiprocessing.Pool(multiprocessing.cpu_count()) as p:
            p.map(process_dir, directories)

        break # make sure we only go one level deep
        
if __name__ == "__main__":
    # print(apply_threshold('test_img1.png'))

    filepath = r"W:\GeneratedData\F013B2\PS2\2019-08-27"
    process_collection(filepath)
    # process_dir(r"W:\GeneratedData\F013B2\PS2\2019-08-27\2019-08-27__03-13-25-016")
    # pprint(apply_threshold(r"W:\GeneratedData\F013B2\PS2\2019-08-27\2019-08-27__03-13-25-016\bdf712f3-12fb-4738-87a0-d757abc7f356_rawData0001.png"))

    # print(apply_threshold(r"W:\GeneratedData\F013B2\PS2\2019-08-27\2019-08-27__03-13-25-016\bdf712f3-12fb-4738-87a0-d757abc7f356_rawData0052.png"))