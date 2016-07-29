#! /usr/bin/env python

__author__ = 'Paul Sarando'

import config.ncbi_sra_submit_properties

from metadata_client import MetadataClient
from genshi.template import TemplateLoader
from lxml import etree
from argparse import ArgumentParser
from subprocess import call

import os
import shutil

class BioProjectUploader:
    def __init__(self, ascp_cmd, private_key_path, ncbi_user, ncbi_host, ncbi_sumbit_path):
        self.ascp_cmd = ascp_cmd
        self.private_key_path = private_key_path
        self.upload_dest = '{0}@{1}:{2}'.format(ncbi_user, ncbi_host, ncbi_sumbit_path)

    def upload_project(self, submit_dir, input_paths):
        # Collect the submission files from the input paths into the submission dir
        src_files = {}
        for path in input_paths:
            filename = os.path.basename(path)
            if filename in src_files:
                raise Exception("Duplicate filenames found in input directory:\n{0}\n{1}".format(src_files[filename], path))
            src_files[filename] = path

            shutil.move(path, os.path.join(submit_dir, filename))

        ascp_cmd = self.ascp_cmd + [
            "-i", self.private_key_path,
            submit_dir,
            self.upload_dest
        ]

        try:
            retcode = call(ascp_cmd)
            if retcode != 0:
                raise Exception("Upload error: {0}".format(-retcode))

            # The file uploads were successful, so upload a 'submit.ready' file to complete the submission.
            submit_ready = "submit.ready"
            open(os.path.join(submit_dir, submit_ready), 'a').close()

            # Calling the same upload command with the same submit directory will skip all files already
            # successfully uploaded, and only upload the new 'submit.ready' file.
            retcode = call(ascp_cmd)
            if retcode != 0:
                raise Exception("Error uploading '{0}' file: {1}".format(submit_ready, -retcode))
        except OSError as e:
            raise Exception("Aspera Connect upload failed", e)

        # Clean up: Move input files back into their original directories,
        # so they are not transferred as outputs, but can be preserved as inputs
        for filename in src_files:
            shutil.move(os.path.join(submit_dir, filename), src_files[filename])

class BioProjectXmlValidator:
    def __init__(self, schemas_dir, submission_schema_path, bioproject_schema_path, biosample_schema_path):
        with open(submission_schema_path, 'r') as f:
            schema_root = etree.XML(f.read(), base_url=schemas_dir)
        schema = etree.XMLSchema(schema_root)
        self.xmlparser = etree.XMLParser(schema=schema)

        with open(bioproject_schema_path, 'r') as f:
            schema_root = etree.XML(f.read(), base_url=schemas_dir)
        self.bioproject_schema = etree.XMLSchema(schema_root)

        with open(biosample_schema_path, 'r') as f:
            schema_root = etree.XML(f.read(), base_url=schemas_dir)
        self.biosample_schema = etree.XMLSchema(schema_root)

    def validate_bioproject_xml(self, submission_path):
        with open(submission_path, 'r') as f:
            submission = etree.fromstring(f.read(), self.xmlparser)

        bioproject = submission.xpath('/Submission/Action/AddData/Data/XmlContent/Project')
        if len(bioproject) > 0:
            self.bioproject_schema.assertValid(bioproject[0])

        bio_samples = submission.xpath('/Submission/Action/AddData/Data/XmlContent/BioSample')
        for biosample in bio_samples:
            self.biosample_schema.assertValid(biosample)


usage = 'ncbi_sra_submit.py [options]'

desc = """
Prepares files and metadata downloaded from the Discovery Environment Data Store for submission to
the NCBI Sequence Read Archive (SRA).
"""

# Parse the command-line options.
parser = ArgumentParser(usage = usage, description = desc, add_help = False)
parser.add_argument('-s', '--submit-mode', dest = 'submit_mode',
                    default = 'create',
                    help = 'specify if the SRA submission is a BioProject "create" (default) or an "update" request.')
parser.add_argument('-i', '--private-key', dest = 'private_key_path',
                    default = config.ncbi_sra_submit_properties.private_key_path,
                    help = '(optional) specify an alternative path to the id_rsa private-key file.')
parser.add_argument('-f', '--input-dir', dest = 'input_dir',
                    help = 'specify the path to the input BioProject folder to submit to SRA')
parser.add_argument('-m', '--input-metadata', dest = 'metadata_path',
                    required = True,
                    help = 'specify the path to the BioProject folder metadata file.')
parser.add_argument('-v', '--validate-metadata-only', dest = 'validate_only', action='store_true',
                    help = 'when included, no data will be submitted and only the BioProject folder'
                           ' metadata file will be validated.')
parser.add_argument('-d', '--submit-dir', dest = 'submit_dir',
                    help = 'specify the path to the destination BioProject SRA submission folder')
parser.add_argument('-?', '--help', action = 'help')
args = parser.parse_args()

# Define the objects we need.
metadata_client = MetadataClient()
loader = TemplateLoader(config.ncbi_sra_submit_properties.templates_dir)

# Parse iPlant Data Store metadata into format usable by the submission templates
metadata = metadata_client.get_metadata(args.metadata_path)

# Create the destination submission directory
submit_dir = args.submit_dir
if not submit_dir:
    submit_dir = '{0}.{1}'.format(os.environ['IPLANT_USER'], metadata['sra_object_id'])
if not os.path.exists(submit_dir):
    os.makedirs(submit_dir)

# The SRA project ID is required if the submission mode is not 'create'
if args.submit_mode != 'create' and 'sra_project_id' not in metadata:
    raise Exception("Could not find SRA Project ID in Bio Project metadata for project update.")

# Generate submission.xml in the submission dir
metadata_template = loader.load('submission.xml')
submission_path = os.path.join(submit_dir, 'submission.xml')
stream = metadata_template.generate(metadata=metadata, submit_mode=args.submit_mode)
with open(submission_path, 'w') as f:
    stream.render(method='xml', out=f)

# Validate generated XML
xml_validator = BioProjectXmlValidator(config.ncbi_sra_submit_properties.schemas_dir,
                                       config.ncbi_sra_submit_properties.submission_schema_path,
                                       config.ncbi_sra_submit_properties.bioproject_schema_path,
                                       config.ncbi_sra_submit_properties.biosample_schema_path)
xml_validator.validate_bioproject_xml(submission_path)

if args.validate_only or not args.input_dir:
    print "Only validated metadata, no data was submitted to the NCBI SRA."
else:
    uploader = BioProjectUploader(config.ncbi_sra_submit_properties.ascp_cmd,
                                  args.private_key_path,
                                  config.ncbi_sra_submit_properties.ncbi_user,
                                  config.ncbi_sra_submit_properties.ncbi_host,
                                  config.ncbi_sra_submit_properties.ncbi_sumbit_path)

    # Build the list of submission input file paths for the uploader
    input_paths = metadata_client.get_bio_project_file_paths(metadata, args.input_dir)

    uploader.upload_project(submit_dir, input_paths)
