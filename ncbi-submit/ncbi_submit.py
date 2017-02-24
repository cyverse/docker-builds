#! /usr/bin/env python

__author__ = 'Paul Sarando'

import bioproject
import config.ncbi_submit_properties

from metadata_client import MetadataClient
from genshi.template import TemplateLoader
from lxml import etree
from argparse import ArgumentParser
from subprocess import call

import os
import shutil

usage = 'ncbi_submit.py [options]'

desc = """
Prepares files and metadata downloaded from the Discovery Environment Data Store for submission to
the NCBI Sequence Read Archive (SRA).
"""

# Parse the command-line options.
parser = ArgumentParser(usage = usage, description = desc, add_help = False)
parser.add_argument(
    '-s', '--submit-mode', dest = 'submit_mode', default = 'create',
    help = 'Specifies whether to create a new bio project or update an existing one.',
    choices = ['create', 'update']
)
parser.add_argument(
    '-i', '--private-key', dest = 'private_key_path',
    default = config.ncbi_submit_properties.private_key_path,
    help = '(optional) specify an alternative path to the id_rsa private-key file.'
)
parser.add_argument(
    '-f', '--input-dir', dest = 'input_dir',
    help = 'specify the path to the input BioProject folder to submit to SRA'
)
parser.add_argument(
    '-m', '--input-metadata', dest = 'metadata_path',
    required = True,
    help = 'specify the path to the BioProject folder metadata file.'
)
parser.add_argument(
    '-v', '--validate-metadata-only', dest = 'validate_only', action='store_true',
    help = 'when included, no data will be submitted and only the BioProject folder'
           ' metadata file will be validated.'
)
parser.add_argument(
    '-d', '--submit-dir', dest = 'submit_dir',
    help = 'specify the path to the destination BioProject SRA submission folder'
)
parser.add_argument('-?', '--help', action = 'help')
args = parser.parse_args()

# Define the objects we need.
metadata_client = MetadataClient()
loader = TemplateLoader(config.ncbi_submit_properties.templates_dir)

# Parse iPlant Data Store metadata into format usable by the submission templates
metadata = metadata_client.get_metadata(args.metadata_path)

# Create the destination submission directory
submit_dir = args.submit_dir
if not submit_dir:
    submit_dir = '{0}.{1}'.format(os.environ['IPLANT_USER'], metadata['object_id'])
if not os.path.exists(submit_dir):
    os.makedirs(submit_dir)

# The SRA project ID is required if the submission mode is not 'create'
if args.submit_mode != 'create' and 'project_id' not in metadata:
    raise Exception("Could not find SRA Project ID in Bio Project metadata for project update.")

# Generate submission.xml in the submission dir
metadata_template = loader.load('submission.xml')
submission_path = os.path.join(submit_dir, 'submission.xml')
stream = metadata_template.generate(metadata=metadata, submit_mode=args.submit_mode)
with open(submission_path, 'w') as f:
    stream.render(method='xml', out=f)

# Validate generated XML
xml_validator = bioproject.get_xml_validator()
xml_validator.validate_bioproject_xml(submission_path)

if args.validate_only or not args.input_dir:
    print "Only validated metadata, no data was submitted to the NCBI SRA."
else:
    uploader = bioproject.get_uploader(private_key_path=args.private_key_path)

    # Build the list of submission input file paths for the uploader
    input_paths = metadata_client.get_bio_project_file_paths(metadata, args.input_dir)

    uploader.upload_project(submit_dir, input_paths)
