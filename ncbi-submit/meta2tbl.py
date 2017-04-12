#! /usr/bin/env python

__author__ = 'Dennis Roberts'

import config.ncbi_submit_properties

from argparse import ArgumentParser
from genshi.template import TemplateLoader
from genshi.template.text import NewTextTemplate
from metadata_client import MetadataClient

import os

usage = 'meta2tbl.py [options]'

desc = """
Generates an SBT template file from the exported metadata for a BioProject folder.
"""

# Parse the command-line options.
parser = ArgumentParser(usage=usage, description=desc, add_help=False)
parser.add_argument(
    '-m', '--input-metadata', dest='metadata_path',
    required = True,
    help = 'specify the path to the BioProject folder metadata file.'
)
parser.add_argument(
    '-o', '--out', default='template.sbt',
    help = '(optional) specify the name of the output file (default: template.sbt)'
)
parser.add_argument('-?', '--help', action = 'help')
args = parser.parse_args()

# Create the metadata client and template loader.
metadata_client = MetadataClient(require_files=False)
loader = TemplateLoader(config.ncbi_submit_properties.templates_dir)

# Load the metadata and template.
metadata = metadata_client.get_metadata(args.metadata_path)
metadata_template = loader.load('template.sbt', cls=NewTextTemplate)

# Generate the SBT file.
stream = metadata_template.generate(metadata=metadata)
with open(args.out, 'w') as f:
    stream.render(method='text', out=f)
