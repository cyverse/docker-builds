#!/usr/bin/python

'''
RseqFilt is an automated sequence filtering analysis tool for a single and paired-end high throughput RNA-seq data generated from Illumina sequencing platforms.
'''

import sys
from setuptools import Extension, setup


print("\nChecking for prerequisite...\n")
if sys.version_info[0] < 2 or sys.version_info[0] > 3:
    sys.stderr.write("Error: Your Python version is not compatible\nThe SRAP works best with Python 2.7 or "
                     "higher\n\n")
    sys.exit(1)
else:
    print("Found Python %s" % sys.version)

try:
    import numpy
    print("Module numpy found")
except ImportError:
    sys.stderr.write("\nError: module numpy not found\nInstall using 'pip3 install numpy'\n\n")
    sys.exit(1)

try:
    import termcolor
    print("Module termcolor found")
except ImportError:
    sys.stderr.write("\nError: module termcolor not found\nInstall using 'pip3 install termcolor'\n\n")
    sys.exit(1)

try:
    import subprocess
    print("Module subprocess found")
except ImportError:
    sys.stderr.write("\nError: module subprocess not found\nInstall using 'pip3 install subprocess'\n\n")
    sys.exit(1)

try:
    import gzip
    print("Module gzip found")
except ImportError:
    sys.stderr.write("\nError: module gzip not found\nInstall using 'pip install gzip'\n\n")
    sys.exit(1)

#   reserved for future release
'''
try:
    import MySQLdb
    print "Module MySQLdb found"
except ImportError:
    sys.stderr.write("\nError: module MySQLdb not found\nInstall using 'pip install MySQL-python'\n\n")
    sys.exit(1)
'''

try:
    import pysam
    print("Module pysam found")
except ImportError:
    sys.stderr.write("\nError: module pysam not found\nInstall using 'pip3 install pysam'\n\n")
    sys.exit(1)

try:
    import shutil
    print("Module shutil found")
except ImportError:
    sys.stderr.write("\nError: module shutil not found\nInstall using 'pip3 install shutil'\n\n")
    sys.exit(1)

try:
    import glob
    print("Module glob found")
except ImportError:
    sys.stderr.write("\nError: module glob not found\nInstall using 'pip3 install glob'\n\n")
    sys.exit(1)

try:
    import collections
    print("Module collections found")
except ImportError:
    sys.stderr.write("\nError: module collection not found\nInstall using 'pip3 install Counter'\n\n")
    sys.exit(1)

try:
    import math
    print("Module math found")
except ImportError:
    sys.stderr.write("\nError: module math not found\nInstall using 'pip3 install math'\n\n")
    sys.exit(1)

try:
    import multiprocessing
    print("Module multiprocessing found")
except ImportError:
    sys.stderr.write("\nError: module multiprocessing not found\nInstall using 'pip3 install multiprocessing'\n\n")
    sys.exit(1)

try:
    import datetime
    print("Module datetime found")
except ImportError:
    sys.stderr.write("\nError: module datetime not found\nInstall using 'pip3 install datetime'\n\n")
    sys.exit(1)

try:
    import matplotlib
    print("Module matplotlib found")
except ImportError:
    sys.stderr.write("\nError: module matplotlib not found\nInstall using 'pip3 install matplotlib'\n\n")
    sys.exit(1)

try:
    import csv
    print("Module csv found")
except ImportError:
    sys.stderr.write("\nError: module csv not found\nInstall using 'pip3 install csv'\n\n")
    sys.exit(1)

try:
    import itertools
    print("Module itertools found")
except ImportError:
    sys.stderr.write("\nError: module itertools not found\nInstall using 'pip3 install itertools'\n\n")
    sys.exit(1)



setup(
    name='filter',
    url='',
    license='MIT',
    author='Renesh Bedre',
    author_email='reneshbe@gmail.com',
    description='RseqFilt: Quality filtering analysis for RNA-seq data',
    classifiers=[
        'Intended Audience :: Developers'
        'Intended Audience :: Education',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 3',
        'Operating System :: POSIX :: Linux',
        'Environment :: Console',
        'Topic :: Scientific/Engineering :: Bioinformatics',
    ],
    scripts=['ngsmodules/filter.py',
             'ngsmodules/Filter_Single.py',
             'ngsmodules/Filter_Pair.py',
             'ngsmodules/StatisticSingle.py',
             'ngsmodules/StatisticPair.py',
             'ngsmodules/common_functions.py',
             ],
    requires=['numpy', 'python (>=3.0)', 'itertools', 'csv', 'matplotlib', 'datetime', 'multiprocessing', 'math',
              'collections', 'glob', 'shutil', 'pysam', 'subprocess', 'termcolor'],


)
