__author__ = 'Paul Sarando'

from argparse import ArgumentParser

import config.properties

usage = '%(prog)s [options]'

desc = """
Parses a large HT Path List file from the Discovery Environment into smaller HT Path List files,
in batches of {0} paths per output file.
""".format(config.properties.path_list_max_paths)

# Parse the command-line options.
parser = ArgumentParser(usage = usage, description = desc, add_help = False)
parser.add_argument(
    '-f', '--path-list-input', dest = 'path_list_input',
    required = True,
    help = 'Specify the path to an HT Path List file.'
)
parser.add_argument(
    '-o', '--output-prefix', dest = 'output_prefix', default=config.properties.output_prefix,
    help = 'Specify the filename prefix for each output HT Path List file,'
           ' which will be appended with a number (defaults to "{0}";'
           ' e.g. {0}1, {0}2, etc.)'.format(config.properties.output_prefix)
)
parser.add_argument('-?', '--help', action = 'help')
args = parser.parse_args()

# Parse paths from the source HT Path List
with open(args.path_list_input) as path_list_input:
    header = path_list_input.readline().strip()

    if header != config.properties.path_list_file_identifier:
        raise Exception(
            "Input file does not have a recognizable HT Path List header.\nExpected: {0}\nReceived: {1}".format(
                config.properties.path_list_file_identifier, header))

    # Write paths in batches to smaller HT Path List files
    batch = 0
    path = path_list_input.readline().strip()
    while path:
        batch += 1
        ht_list_filename = "{0}{1}".format(args.output_prefix, str(batch).zfill(4))

        with open(ht_list_filename, 'w') as ht_list_file:
            ht_list_file.write(config.properties.path_list_file_identifier)
            ht_list_file.write('\n')

            for n in range(config.properties.path_list_max_paths):
                if path:
                    ht_list_file.write(path)
                    ht_list_file.write('\n')
                else:
                    break

                path = path_list_input.readline().strip()
