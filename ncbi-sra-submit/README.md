### NCBI SRA Submission

#### Description and Quick Start

Publish data to the NCBI Sequence Read Archive.
See https://pods.iplantcollaborative.org/wiki/pages/viewpage.action?pageId=20351132 for a description of the individual Apps and submission workflow.

#### Input File(s)

A BioProject folder is required for data submission to the NCBI Sequence Read Archive (SRA).
A BioProject folder contains 1 or more BioSample folders,
each of which contain 1 or more Library folders,
and each of those in turn contain 1 or more compressed data files for submission to the NCBI SRA.
For example, a BioProject folder may have this structure:

* BioProject/bio_sample_1/library1/data1.gz
* BioProject/bio_sample_2/library1/data2.gz
* BioProject/bio_sample_3/library1/data3.gz
* BioProject/bio_sample_3/library2/data4.gz

A BioProject folder metadata export file is required as the second input.
This file can be generated from the Data window's `Edit` menu by selecting the BioProject folder then selecting the `Save Metadata` menu item.
The appropriate metadata should be saved on each folder of the BioProject.
There are Metadata Templates available for each folder of a BioProject to facilitate adding all required and optional metadata attributes to each folder;
including a Bio Project template, Bio Sample templates, and Bio Sample Library templates.

* **Note:** All data files under a BioProject must be compressed before submitting to the SRA with this App,
            and each filename under a BioProject must be unique.
* **Note:** If a BioProject folder is not selected,
            then this App will only validate the given metadata file,
            and no data will be submitted to the NCBI SRA.

#### Parameters Used in App

By default the `Create a new Bio Project Submission` checkbox is selected.
Uncheck this checkbox to submit an update to an existing BioProject.
Note that an existing Bio Project ID must be set as a metadata attribute on the BioProject folder in the metadata export file for `Update` submissions.

#### Output File(s)

This App will generate a folder, named with the username of the submitter and the BioProject folder ID,
which contains a `submission.xml` file generated from the given metadata input file.
This XML file is also submitted to the NCBI SRA along with the compressed data files when a BioProject folder is selected as the first input.

Aspera Connect manifests and transfer log files will be returned with other job log files in the `logs` output folder.
