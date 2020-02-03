__author__ = 'Sarah Roberts'

import config.ncbi_submit_properties
import os
import re
import requests
import shutil

from lxml import etree
from subprocess import call

def load_schema(schemas_dir, path):
    with open(path, 'r') as f:
        schema_root = etree.XML(f.read(), base_url=schemas_dir)
    return etree.XMLSchema(schema_root);

def load_parser(schemas_dir, path):
    schema = load_schema(schemas_dir, path)
    return etree.XMLParser(schema=schema)

class BioProjectXmlValidator:
    def __init__(self, schemas_dir, schema_paths):
        self.xmlparser = load_parser(schemas_dir, schema_paths['submission'])
        self.bioproject_schema = load_schema(schemas_dir, schema_paths['bioproject'])
        self.biosample_schema = load_schema(schemas_dir, schema_paths['biosample'])
        self.genome_schema = load_schema(schemas_dir, schema_paths['genome'])

    def validate_bioproject_xml(self, submission_path):
        with open(submission_path, 'r') as f:
            submission = etree.fromstring(f.read(), self.xmlparser)

        bioproject = submission.xpath('/Submission/Action/AddData/Data/XmlContent/Project')
        if len(bioproject) > 0:
            self.bioproject_schema.assertValid(bioproject[0])

        bio_samples = submission.xpath('/Submission/Action/AddData/Data/XmlContent/BioSample')
        for biosample in bio_samples:
            self.biosample_schema.assertValid(biosample)

class SubmissionValidator:
    def __init__(self, validation_url):
        self.validation_url = validation_url

    def validate_submission(self, submission_path, report_path):
        with open(submission_path, 'r') as f:
            submission = f.read()

        # Call the submission validation endpoint.
        headers = {'Content-Type': 'text/xml'}
        r = requests.post(self.validation_url, data=submission, headers=headers)
        r.raise_for_status()

        # Write a copy of the report to disk for reference.
        with open(report_path, 'w') as out:
            out.write(r.text)

        # Parse the report.
        report = etree.fromstring(re.sub(r'^<\?xml[^>]+>\s+', '', r.text))

        # Extract and print the messages in the report for convenience.
        messages = report.xpath('//Message') or []
        if len(messages) != 0:
            print 'Messages found in Report:'
            for message in messages:
                attrs = message.attrib
                severity = attrs['severity'] if 'severity' in attrs else 'unknown-severity'
                print '\t', severity, ':', message.text
            print
            print 'Please see report.xml for details.'
            print

        # Check the report for errors.
        statuses = report.xpath('/BioSampleValidate/Action/@status') or []
        if len(statuses) == 0:
            raise Exception('status not found in downloaded report')
        if any(status == 'processed-error' for status in statuses):
            raise Exception('NCBI validation failed: please see report.xml for details')

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
                fmt = "Duplicate filenames found in input directory:\n{0}\n{1}"
                raise Exception(fmt.format(src_files[filename], path))
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

def get_xml_validator():
    return BioProjectXmlValidator(
        config.ncbi_submit_properties.schemas_dir,
        config.ncbi_submit_properties.schema_paths
    )

def get_submission_validator():
    return SubmissionValidator(config.ncbi_submit_properties.validation_url)

def get_uploader(private_key_path=None):
    key_path = config.ncbi_submit_properties.private_key_path if private_key_path is None else private_key_path
    return BioProjectUploader(
        config.ncbi_submit_properties.ascp_cmd,
        key_path,
        config.ncbi_submit_properties.ncbi_user,
        config.ncbi_submit_properties.ncbi_host,
        config.ncbi_submit_properties.ncbi_sumbit_path
    )
