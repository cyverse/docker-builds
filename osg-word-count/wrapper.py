#!/usr/bin/env

import argparse
import csv
import itertools
import json
import subprocess

from os.path import expanduser, basename

# The maximum number of attempts to upload or download a file.
max_attempts = 3

# This is a simple context manager class designed to make it easier to read and write iRODS ticket list files.
class TicketListReader:
    def __init__(self, path):
        self.path = path

    def __enter__(self):
        self.fd = open(self.path, "rb")
        self.r = csv.reader(itertools.ifilter(lambda l: l[0] != '#', self.fd))
        return self.r

    def __exit__(self, type, value, traceback):
        self.fd.close()

# Initialize iRODS.
def init_irods(host, port):
    irods_config = {
        "irods_user_name": "anonymous",
        "irods_host": host,
        "irods_port": port,
        "irods_zone_name": ""
    }
    with open(expanduser("~/.irods/irods_environment.json"), "w") as f:
        json.dump(irods_config, f, indent=4, separators=(",", ": "))

# Download a file or directory from iRODS.
def download_file(ticket, src):
    for i in range(1, max_attempts):
        rc = subprocess.call(["iget", "-rt", ticket, src])
        if rc == 0:
            return
    raise "could not download {0} after {1} attempts".format(src, max_attempts)

# Download a set of files referenced in a ticket list file from iRODS, returning a list of downloaded files.
def download_files(ticket_list_path):
    result = []
    with TicketListReader(ticket_list_path) as tlr:
        for ticket, src in tlr:
            download_file(ticket, src)
            result.append(basename(src))
    return result

# Run the word count job.
def run_job(input_files, output_filename, error_filename):
    with open(output_filename, "w") as out, open(error_filename, "w") as err:
        rc = subprocess.call(["wc"] + input_files, stdout=out, stderr=err)
        if rc != 0:
            raise "wc returned a result code of {0}".format(rc)

# Upload a file or directory to iRODS.
def upload_file(ticket, src, dest):
    for i in range(1, max_attempts):
        rc = subprocess.call(["iput", "-rt", ticket, src, dest])
        if rc == 0:
            return
    raise "could not upload {0} to {1} after {2} attempts".format(src, dest, max_attempts)

# Upload a set of files to the directories referenced in a ticket list file to iRODS.
def upload_files(ticket_list_path, files):
    with TicketListReader(ticket_list_path) as tlr:
        for ticket, dest in tlr:
            for src in files:
                upload_file(ticket, src, dest)

if __name__ == "__main__":
    # Parse the command line.
    parser = argparse.ArgumentParser("Wrapper for word count utility.")
    parser.add_argument("--irods-host", required=True, help="The iRODS host name or IP address.")
    parser.add_argument("--irods-port", type=int, default=1247, help="The iRODS port number.")
    parser.add_argument("--input-ticket-list", required=True, help="The path to the input ticket list file.")
    parser.add_argument("--output-ticket-list", required=True, help="The path to the output ticket list file.")
    parser.add_argument("--stdout", default="out.txt", help="The name of the output file.")
    parser.add_argument("--stderr", default="err.txt", help="The name of the error output file.")
    args = parser.parse_args()

    # Initialize iRODS.
    init_irods(args.irods_host, args.irods_port)

    # Download the files from iRODS, assuming that all tickets refer to regular files for now.
    input_files = download_files(args.input_ticket_list)

    # Run the job.
    run_job(input_files, args.stdout, args.stderr)

    # Upload the output file.
    upload_files(args.output_ticket_list, [args.stdout, args.stderr])
