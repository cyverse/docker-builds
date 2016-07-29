#! /usr/bin/env python

__author__ = 'Paul Sarando'

import config.ncbi_sra_submit_properties

from argparse import ArgumentParser
from subprocess import call

class BioProjectReportDownloader:
    def __init__(self, ascp_cmd, private_key_path, ncbi_user, ncbi_host, ncbi_sumbit_path):
        self.ascp_cmd = ascp_cmd
        self.private_key_path = private_key_path
        self.report_src = '{0}@{1}:{2}'.format(ncbi_user, ncbi_host, ncbi_sumbit_path)

    def download_report(self, sumbit_dir, download_dest):
        report_xml = '{0}/{1}/report.xml'.format(self.report_src, sumbit_dir)

        ascp_cmd = self.ascp_cmd + [
            "-i", self.private_key_path,
            report_xml,
            download_dest
        ]

        try:
            retcode = call(ascp_cmd)
            if retcode != 0:
                raise Exception("Download error: {0}".format(-retcode))
        except OSError as e:
            raise Exception("Aspera Connect download failed", e)

usage = """
ncbi_sra_report_download.py [-i <PRIVATE_KEY_PATH>] -s <SRC_SUBMIT_DIR> -d <DOWNLOAD_OUTPUT_DIR>
"""

desc = """
Downloads a report.xml file from the NCBI Sequence Read Archive (SRA).
"""

# Parse the command-line options.
parser = ArgumentParser(usage = usage, description = desc, add_help = False)
parser.add_argument('-i', '--private-key', dest = 'private_key_path',
                    default = config.ncbi_sra_submit_properties.private_key_path,
                    help = '(optional) specify an alternative path to the id_rsa'
                           ' private-key file.')
parser.add_argument('-s', '--submit-dir', dest = 'submit_dir',
                    required = True,
                    help = 'specify the name of the BioProject SRA submission'
                           ' folder where the report.xml file is located.')
parser.add_argument('-d', '--output-dir', dest = 'output_dir',
                    default = ".",
                    help = 'specify the path to the output directory where the SRA'
                           ' report.xml file will be downloaded.')
parser.add_argument('-?', '--help', action = 'help')
args = parser.parse_args()

# Define the objects we need.
downloader = BioProjectReportDownloader(config.ncbi_sra_submit_properties.ascp_cmd,
                                        args.private_key_path,
                                        config.ncbi_sra_submit_properties.ncbi_user,
                                        config.ncbi_sra_submit_properties.ncbi_host,
                                        config.ncbi_sra_submit_properties.ncbi_sumbit_path)

# Download the SRA report
downloader.download_report(args.submit_dir, args.output_dir)
