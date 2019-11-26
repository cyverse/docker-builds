## FLIP (Fluorescence Imaging Pipeline)

### Getting Started 
---
To set up an evironment open a terminal in the project folder and type

```
$ python -m venv venv
```
to create a virtual environment 

#### For windows

To activate the vm, type 

```
$ venv\Scripts\activate 
```
You should see (venv) at the front of your command line

#### For linux

To activate the vm, type

```
$ source venv/bin/activate
```

To install dependancies into the vm, type

`pip install -r requirements.txt` 

<br/><br/>

#### CLI (command line interface)

To run the cli, run 

```
$ python FLIP.py -d <ps2 collection path>
```

This will run through the binary to png conversion, run multithreshold image segmentation, and then generate aggregate and fluorescence files all in one go.

<br/><br/>

#### GUI (graphical user interface)

To open the gui, run `python FLIP.py` with no arguments

`Binary To PNG` button converts all the .bin images in a folder to .png

`ImageJ Macro` button looks at the images generated from the binary to png button and creates a .csv with some calculations. It will try to find an imagej installation in the project folder, but if it can't, it will ask where one is. After that it will ask for the location of an imagej macro to be run

`Python Macro` button looks at the images generated from the binary to png button and applies multi threshold image segmentation. It will ask for the location of the ps2 collection

`Generate Aggregate and Fluorescence` button generates `{foldername}_aggregated.csv` and `{foldername}_fluorescence.csv` for a collection of ps2 images. Each subfolder in a collection must have the ..._metadata.json file and the {foldername}.csv to be processed.
<br/><br/>

```
2019-08-27/
    2019-08-27__00-00-09-654/
        ..._metadata.json
        2019-08-27__00-00-09-654.csv

    2019-08-27__00-00-52-305/
        .._metadata.json
        2019-08-27__00-00-52-305.csv

    2019-08-27__00-01-34-971/
        .._metadata.json
        2019-08-27__00-01-34-971.csv

    2019-08-27_aggregated.csv
    2019-08-27_fluorescence.csv
```
---
## Associate Plots
the `associate_plots.py` file will create either a .json or a .csv associating images to their plots based on `Plot boundaries.xlsx`.
<br/><br/>

#### CLI (command line interface)

```
python associate_plots.py
    -f <filepath>
    -t <csv or json>
    -xo <x_offset> (optional. 0 by default)
    -yo <y_offset> (optional. 0 by default)
```
<br/>

#### GUI (graphical user interface)

To generate {foldername}_plot_xyz.{json/csv}, run `python associate_plots.py` with no arguments. This will open a gui that will let you pick between a json and a csv file. After choosing the output filetype, press one of the buttons to open a file dialog. You can then choose a directory that contains images you want to be associated. This will look in each of the folders in a directory find the _metadata.json, and create either a csv or json file with a list of each .bat or raw image for that plot. If the direcory cannot be associated to a plot, the plot number will be -1.
