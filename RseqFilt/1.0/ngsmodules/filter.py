
import argparse
import subprocess
import sys
from termcolor import colored

parser = argparse.ArgumentParser(description="Quality filtering analysis of single and paired-end sequence data")

parser.add_argument('-a', '--p1', action='store', type=str, dest='input_files_1', help='Single end input files or left '
                                                                                       'files for paired-end data '
                                                                                       '(.fastq, .fq). Multiple sample '
                                                                                       'files must be separated by comma',
                    default=None)
parser.add_argument('-b', '--p2', action='store', type=str, dest='input_files_2', help='Right files for paired-end data '
                                                                                       '(.fastq, .fq). Multiple files '
                                                                                       'must be separated by comma ',
                    default=None)
parser.add_argument('-c', '--qfmt', action='store', type=str, dest='qual_fmt', help='Quality value format '
                                                                                    '[1= Illumina 1.8, 2= Illumina 1.3,'
                                                                                    '3= Sanger]. If quality format not '
                                                                                    'provided, it will automatically '
                                                                                    'detect based on sequence data',
                    default=0)
parser.add_argument('-e', '--nb', action='store', type=str, dest='n_cont', help='Filter the reads containing given %% of '
                                                                                'uncalled bases (N)', default=101)
parser.add_argument('-f', '--adp', action='store', type=str, dest='adpt_seqs', help='Trim the adapter and '
                                                                                    'truncate the read '
                                                                                    'sequence (multiple adapter sequences'
                                                                                    ' must be separated by comma)',
                    default='NULL')
parser.add_argument('-d', '--msz', action='store', type=int, dest='min_size', help='Filter the reads which are lesser '
                                                                                      'than minimum size', default=0)
parser.add_argument('-g', '--per', action='store', type=float, dest='adpt_match', help='Truncate the read sequence if '
                                                                                       'it matches to '
                                                                           'adapter sequence equal or more than given '
                                                                           'percent (0.0-1.0) [default=0.9]', default=1)
parser.add_argument('-i', '--qthr', action='store', type=int, dest='qual_thresh', help='Filter the read sequence if '
                                                                                       'average quality of bases in '
                                                                                       'reads is lower than threshold '
                                                                                       '(1-40) [default:20]', default=20)
parser.add_argument('-n', '--trim', action='store', type=str, dest='trim_opt', help='If trim option set to True, the '
                                                                                     'reads with low quality (as '
                                                                                     'defined by option --qthr) will be '
                                                                                     'trimmed instead of discarding [True|False] '
                                                                                     '[default: False]', default='False')
parser.add_argument('-p', '--wsz', action='store', type=int, dest='wind_size', help='The window size for trimming '
                                                                                      '(5->3) the reads. This option '
                                                                                      'should always set when -trim '
                                                                                      'option is defined [default: 5]', default=5)
parser.add_argument('-r', '--mlk', action='store', type=int, dest='min_len_filt', help='Minimum length of the reads '
                                                                                       'to retain after trimming',
                    default=0)
parser.add_argument('-q', '--cpu', action='store', type=int, dest='cpu', help='Number of CPU [default:2]', default=2)
parser.add_argument('-m', '--ofmt', action='store', type=str, dest='out_fmt', help='Output file format (fastq/fasta) '
                                                                               '[default:fastq]', default='fastq')
parser.add_argument('-v', '--no-vis', action='store', type=str, dest='vis_opt', help='No figures will be produced '
                                                                                        '[True|False] [default:False]',
                    default='False')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

# print help message if no arguments provided
if len(sys.argv) < 2:
    parser.print_help()
    sys.exit(1)

results = parser.parse_args()

# run single-end filtering if --p2 is not provided
if results.input_files_2 is None:
    if results.input_files_1 is None:
        print(colored("Error: input file is missing \n", "red"))
        sys.exit(1)
    fastq_files = results.input_files_1.split(',')
    for file in fastq_files:
        # print("Filtering reads:", file)
        p1 = subprocess.Popen(["Filter_Single.py", "--p1", str(file), "--qfmt", str(results.qual_fmt),
                               "--nb", str(results.n_cont), "--adp", str(results.adpt_seqs),
                               "--msz", str(results.min_size), "--per", str(results.adpt_match),
                               "--qthr", str(results.qual_thresh), "--trim", results.trim_opt,
                               "--wsz", str(results.wind_size), "--mlk", str(results.min_len_filt),
                               "--cpu", str(results.cpu), "--ofmt", str(results.out_fmt),
                               "--no-vis", str(results.vis_opt)])

        p1.wait()
        if p1.returncode != 0:
            print(colored("Error: filtering exited with error status\n", "red"))
            sys.exit(1)
else:
    if results.input_files_2 is None:
        print(colored("Error: right read input file is missing \n", "red"))
        sys.exit(1)
    fastq_files_1 = results.input_files_1.split(',')
    fastq_files_2 = results.input_files_2.split(',')
    if len(fastq_files_1) != len(fastq_files_2):
        print(colored("Error: filtering exited with error status\nunequal number of files\n", "red"))
        sys.exit(1)
    for file1, file2 in zip(fastq_files_1, fastq_files_2):
        # print("Filtering reads:", file1, file2)
        p1 = subprocess.Popen(["Filter_Pair.py", "--p1", str(file1), "--p2", str(file2), "--qfmt",
                               str(results.qual_fmt), "--nb", str(results.n_cont), "--adp", str(results.adpt_seqs),
                               "--msz", str(results.min_size), "--per", str(results.adpt_match),
                               "--qthr", str(results.qual_thresh), "--trim", results.trim_opt,
                               "--wsz", str(results.wind_size), "--mlk", str(results.min_len_filt),
                               "--cpu", str(results.cpu), "--ofmt", str(results.out_fmt),
                               "--no-vis", str(results.vis_opt)])
        p1.wait()
        if p1.returncode != 0:
            print(colored("Error: filtering exited with error status\n", "red"))
            sys.exit(1)




