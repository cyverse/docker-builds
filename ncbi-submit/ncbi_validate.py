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

usage = 'ncbi_validate.py [options]'

desc = """
Validates a submission XML file that has already been generated.
"""

# Parse the command-line options.
parser = ArgumentParser(usage = usage, description = desc, add_help = False)
parser.add_argument('-f', '--submit-file', dest = 'submit_file',
                    help = 'specify the path to the submission file')
parser.add_argument('-?', '--help', action = 'help')
args = parser.parse_args()

# Validate generated XML
xml_validator = bioproject.get_xml_validator()
xml_validator.validate_bioproject_xml(args.submit_file)

# If we get here then the validation succeeded.
print "The submission file is valid."
