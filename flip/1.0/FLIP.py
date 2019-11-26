#!/usr/bin/python3

"""
Fluorescence Imaging Pipeline main file

Jacob Long USDA jacob.long@usda.gov
2019-11-7
"""

import json
import os

import bin_conversion
import fluorescence_aggregation
import image_segmentation

# getting the directory that this file is running from
this_dir = os.path.dirname(os.path.abspath(__file__))

def headless(directory, output, processes=-1):
    """
    run the full pipeline on a directory without a gui
    
    1. convert bins to pngs
    2. do image segmentation (this is where the imagej macro 
        would go, but it has been rewritten)
    3. get ps2 offsets
    4. generate aggregate files
    5. generate fluorescence file 
    """

    print("Converting bins to pngs")
    bin_conversion.convert_dirs(directory, processes=processes) # converting all images in a dir from bin to png

    print("segmenting image")
    # imagej_macro(args.macro, args.directory) # old fiji macro
    image_segmentation.process_collection(directory) # new python macro

    # getting offsets
    with open(os.path.join(this_dir, "static", "ps2", "sensor_fixed_metadata.json"), 'r') as f:
        ps2_dict = json.load(f)
        ps2_x = float(ps2_dict[0]['location_in_camera_box_m']['x'])
        ps2_y = float(ps2_dict[0]['location_in_camera_box_m']['y'])

    print("generating aggregates")
    fluorescence_aggregation.generate_aggregate(
        directory,
        output,
        os.path.join(this_dir, "static", 'Plot boundaries.xlsx'),
        os.path.join(this_dir, "static", 'multithresh.json'),
        ps2_x, 
        ps2_y
    ) # generating aggregates 
    fluorescence_aggregation.generate_fluorescence(directory, output, True) # generate fluorescence files

    print("Finished processing", directory)

def run_imagej_macro():
    """
    runs a command that will call imagej with a specified macro
    """

    # trying to find a local imagej installation, but if one is not found, ask the user where it is
    imagej_path = os.path.join(this_dir, 'Fiji.app/ImageJ-win64.exe')
    if not os.path.exists(imagej_path):
        imagej_path = fd.askopenfilename(title="Could not find imagej, please select your imagej.exe")

    # asking the user where their imagej macro is.
    macro_path = fd.askopenfilename(title="Please select an imagej macro")

    # if the user cancels either of the selections, return
    if not macro_path or not imagej_path:
        return

    subprocess.call([
        imagej_path,
        "-macro",
        macro_path
    ])

def gui():
    import subprocess
    import tkinter
    import tkinter.filedialog as fd
    import tkinter.ttk

    # setting up tkinter window
    root = tkinter.Tk()
    root.title("Generate Fluorescence Aggregates")
    # root.geometry("300x125")
    # root.configure(background='white')
    root.resizable(False, False)

    # configuring grid for resizable buttons
    tkinter.Grid.columnconfigure(root, 0, weight=1)
    for y in range(3):
        tkinter.Grid.rowconfigure(root, y, weight=1)

    # creating a button for each main function
    bin_to_binary_btn = tkinter.ttk.Button(
        root, 
        text="Binary To PNG", 
        command=lambda: bin_conversion.convert_dirs(fd.askdirectory(title="Please select a ps2 collection"))
    )
    bin_to_binary_btn.grid(column=0, row=0, sticky="nsew", padx=2, pady=2, ipadx=10, ipady=7, columnspan=2)

    imagej_macro_btn = tkinter.ttk.Button(
        root, 
        text="ImageJ Macro", 
        command=run_imagej_macro
    )
    imagej_macro_btn.grid(column=0, row=1, sticky="nsew", padx=2, pady=2, ipadx=10, ipady=7)

    image_segmentation_btn = tkinter.ttk.Button(
        root, 
        text="Python Macro", 
        command=lambda: image_segmentation.process_collection(fd.askdirectory(title="Please select a ps2 collection"))
    )
    image_segmentation_btn.grid(column=1, row=1, sticky="nsew", padx=2, pady=2, ipadx=10, ipady=7)

    # getting offsets
    with open(os.path.join(this_dir, "static", "ps2", "sensor_fixed_metadata.json"), 'r') as f:
        ps2_dict = json.load(f)
        ps2_x = float(ps2_dict[0]['location_in_camera_box_m']['x'])
        ps2_y = float(ps2_dict[0]['location_in_camera_box_m']['y'])

    # calling single process
    # using a file dialog, ask the user for a collection path
    # use the Plot boundaries.xlsx and multithresh json next to this file
    # use offsets from sensor_fixed_metadata.json
    generate_aggregate_btn = tkinter.ttk.Button(
        root, 
        text="Generate Aggregate and Fluorescence", 
        command=lambda: fluorescence_aggregation.single_process(
            filepath=fd.askdirectory(title="Please select a ps2 collection"),
            plot_boundaries=os.path.join(this_dir, 'static', 'Plot boundaries.xlsx'),
            multithresh_json=os.path.join(this_dir, 'static', 'multithresh.json'),
            offset_x=ps2_x,
            offset_y=ps2_y
        ),
    )
    generate_aggregate_btn.grid(column=0, row=2, sticky="nsew", padx=2, pady=2, ipadx=10, ipady=7, columnspan=2)

    # starting gui
    root.mainloop()


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--directory', help="directory to ps2 collection")
    parser.add_argument('-o', '--output', help="directory for output files")
    parser.add_argument('-p', '--processes', help="max spawnable processes used by multiprocessing")
    # parser.add_argument('-m', '--macro', help="filepath to an imagej macro", required=True)

    args = parser.parse_args()

    # if the user gave command line arguments, then use the cli, else open a gui
    if args.directory:
    
        output = args.output

        # if processes are given, use that, else, use max possible
        processes = args.processes if args.processes else -1

        headless(args.directory, output, processes)

    else:
        gui()
