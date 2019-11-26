#!/usr/bin/python3

"""
generates a csv or json with all image names and 
the plots in which they were captured.

Jacob Long 
2019-9-13
"""

import json
import os

import pandas as pd

# getting the directory that this file is running from
this_dir = os.path.dirname(os.path.abspath(__file__))

def associate_plots(filepath, output_type="", x_offset=0, y_offset=0):
    """generates a csv or json with all image names and 
    the plots in which they were captured.
    
    filepath: a path to a directory containing subdirectories 
        with asd json files and images.

    output_type: 
        options: '', 'json', 'csv'
        file containing the folders, plots, xyz, and a list of images.
    
    """

    if not filepath:
        return

    # creating a list of dictionaries filled with foldername, x, y, and a list of image names
    associated_plot_dicts = []

    for root, dirs, files in os.walk(filepath):

        # if a json file does not exist, the images in the folder cannot be associated to a plot
        try:        
            json_dir = [f for f in files if f.endswith("_metadata.json")][0]

            image_names = [f for f in files if f.endswith(".bin") or f.endswith("raw")]

        except:
            print(f"Skipping {root}. _metadata.json not found")
            continue

        # loading in json data to a dictionary
        with open(os.path.join(root, json_dir)) as f:
            json_dict = json.load(f)

        # extracting x and y coordinates and applying offsets
        try:
            x = float(json_dict['lemnatec_measurement_metadata']['gantry_system_variable_metadata']['position x [m]']) + x_offset
            y = float(json_dict['lemnatec_measurement_metadata']['gantry_system_variable_metadata']['position y [m]']) + y_offset
            z = float(json_dict['lemnatec_measurement_metadata']['gantry_system_variable_metadata']['position z [m]'])

        except Exception as e:
            print(e)
            print(f"incorrect json file found in {filepath}")
            continue

        # finding the plot that the folder was taken in
        # reading a file defining plot boundaries
        plot_boundaries_df = pd.read_excel(os.path.join(this_dir, "static", 'Plot boundaries.xlsx'))
        plot_boundaries_df = plot_boundaries_df[['Plot', 'X Start', 'X End', 'Y Start', 'Y End']]

        plot = plot_boundaries_df.loc[
            (plot_boundaries_df['X Start'] <= x) & (plot_boundaries_df['X End'] > x) &
            (plot_boundaries_df['Y Start'] <= y) & (plot_boundaries_df['Y End'] > y), 
            'Plot'
        ]

        # assigning a plot number. if the list is empty, the images 
        # were not in the given plots, so make plot num -1
        try:
            plot_num = list(plot)[0]
        except:
            plot_num = -1

        print("Plot", plot_num, end='\r')

        # creating a dictionary to describe a plot and the images within it
        plot_dict = {
            "folder_name": os.path.basename(root),
            "x":           x,
            "y":           y,
            "z":           z,
            "plot":        plot_num,
            "image_names": image_names
        }

        associated_plot_dicts.append(plot_dict)

    
    # if the user wants a file, then create one, else just return the dictionary
    if output_type.endswith('json'):
        # with open(f'{filepath}_plot_xyz.json', 'w') as f:
        with open(os.path.join(filepath, f'{os.path.basename(filepath)}_plot_xyz.json'), 'w') as f:
            f.write(json.dumps(associated_plot_dicts, indent=4))
    
    elif output_type.endswith('csv'):
        # parsing each dictionary to get one record per image, 
        # rather than one entry in the dictionary per plot
        rows = []
        for plot_dict in associated_plot_dicts:
            for img in plot_dict['image_names']:
                rows.append(
                    {
                        "folder_name": plot_dict['folder_name'],
                        "x":           plot_dict['x'],
                        "y":           plot_dict['y'],
                        "z":           plot_dict['z'],
                        "plot":        plot_dict['plot'],
                        "image_name":  img
                    }
                )

        csv_df = pd.DataFrame(rows)
        csv_df.to_csv(os.path.join(filepath, f'{os.path.basename(filepath)}_plot_xyz.csv'), index=False)

    return associated_plot_dicts

def gui():
    import tkinter
    import tkinter.ttk
    import tkinter.filedialog as fd

    # setting up tkinter window
    root = tkinter.Tk()
    root.title("Associate Plots")
    # root.geometry("300x100")
    # root.configure(background='white')
    root.resizable(False, False)

    # configuring grid for resizable buttons
    for x in range(3):
        tkinter.Grid.columnconfigure(root, x, weight=1)
    for y in range(3):
        tkinter.Grid.rowconfigure(root, y, weight=1)


    tkinter.Label(root, text='Output Filetype:').grid(column=0, row=0, columnspan=2)
    radio_frame = tkinter.Frame(root)
    radio_frame.grid(column=0, row=1, padx=10, pady=10)

    filetype = tkinter.StringVar()
    filetype.set('csv')

    tkinter.ttk.Radiobutton(
        radio_frame,
        text='csv',
        variable=filetype,
        value='csv',
    ).grid(column=0, row=0)

    tkinter.ttk.Radiobutton(
        radio_frame,
        text='json',
        variable=filetype,
        value='json',
    ).grid(column=1, row=0)

    # getting offsets
    with open(os.path.join(this_dir, "static", "ps2", "sensor_fixed_metadata.json"), 'r') as f:
        ps2_dict = json.load(f)
        ps2_x = float(ps2_dict[0]['location_in_camera_box_m']['x'])
        ps2_y = float(ps2_dict[0]['location_in_camera_box_m']['y'])

    with open(os.path.join(this_dir, "static", "stereo", "sensor_fixed_metadata.json"), 'r') as f:
        stereo_dict = json.load(f)
        stereo_x = float(stereo_dict[0]['location_in_camera_box_m']['x'])
        stereo_y = float(stereo_dict[0]['location_in_camera_box_m']['y'])

    with open(os.path.join(this_dir, "static", "SWIR", "sensor_fixed_metadata.json"), 'r') as f:
        swir_dict = json.load(f)
        swir_x = float(swir_dict[0]['location_in_camera_box_m']['x'])
        swir_y = float(swir_dict[0]['location_in_camera_box_m']['y'])

    with open(os.path.join(this_dir, "static", "VNIR", "sensor_fixed_metadata.json"), 'r') as f:
        vnir_dict = json.load(f)
        vnir_x = float(vnir_dict[0]['location_in_camera_box_m']['x'])
        vnir_y = float(vnir_dict[0]['location_in_camera_box_m']['y'])

        # print(ps2_x, ps2_y)
        # print(stereo_x, stereo_y)
        # print(swir_x, swir_y)
        # print(vnir_x, vnir_y)

    # creating a button for each main function
    # PS2 BUTTONS
    associate_plots_btn_ps2 = tkinter.ttk.Button(
        root, 
        text="Associate Plots PS2", 
        command=lambda: associate_plots(fd.askdirectory(title="Choose a PS2 or VNIR filepath"), filetype.get(), x_offset=ps2_x, y_offset=ps2_y)
    )
    associate_plots_btn_ps2.grid(column=0, row=2, sticky="nsew", padx=2, pady=2, ipadx=10, ipady=7)

    # SWIR BUTTONS
    associate_plots_btn_stereo = tkinter.ttk.Button(
        root, 
        text="Associate Plots Stereo", 
        command=lambda: associate_plots(fd.askdirectory(title="Choose a PS2 or VNIR filepath"), filetype.get(), x_offset=stereo_x, y_offset=stereo_y)
    )
    associate_plots_btn_stereo.grid(column=0, row=3, sticky="nsew", padx=2, pady=2, ipadx=10, ipady=7)

    # SWIR BUTTONS
    associate_plots_btn_swir = tkinter.ttk.Button(
        root, 
        text="Associate Plots SWIR", 
        command=lambda: associate_plots(fd.askdirectory(title="Choose a PS2 or VNIR filepath"), filetype.get(), x_offset=swir_x, y_offset=swir_x)
    )
    associate_plots_btn_swir.grid(column=0, row=4, sticky="nsew", padx=2, pady=2, ipadx=10, ipady=7)

    # VNIR BUTTONS
    associate_plots_btn_vnir = tkinter.ttk.Button(
        root, 
        text="Associate Plots VNIR", 
        command=lambda: associate_plots(fd.askdirectory(title="Choose a PS2 or VNIR filepath"), filetype.get(), x_offset=vnir_x, y_offset=vnir_x)
    )
    associate_plots_btn_vnir.grid(column=0, row=5, sticky="nsew", padx=2, pady=2, ipadx=10, ipady=7)

    root.mainloop()

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('-f', '--filepath', help="filepath: a path to a directory containing subdirectories with asd json files and images.")
    parser.add_argument('-t', '--type', help="json or csv output filetype")
    parser.add_argument('-xo', '--x_offset')
    parser.add_argument('-yo', '--y_offset')

    args = parser.parse_args()

    # if there are arguments passed in, just use a cli
    if args.filepath and args.type:
        associate_plots(
            filepath=args.filepath, 
            output_type=args.type, 
            x_offset=args.x_offset if args.x_offset else 0, 
            y_offset=args.y_offset if args.y_offset else 0
        )

    # if no arguments are given, make a tkinter gui and have the user select some stuff
    else:
        gui()