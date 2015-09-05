#!/usr/bin/python
"""Parse environment variables in postactivate file and write them
to supervisor.conf file"""
from __future__ import absolute_import, unicode_literals
import sys
import getopt
import re
import subprocess
import os

def main(argv):
    input_file = ''
    output_file = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print 'cp_p2s.py -i <postactivate> -o <supervisor.conf>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'cp_p2s.py -i <postactivate> -o <supervisor.conf>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            input_file = arg
        elif opt in ("-o", "--ofile"):
            output_file = arg

    if input_file and output_file:
        # Parse postactivate and concatenate
        vars = ''
        with open(input_file, 'r') as f:
            for line in f:
                if line.strip().startswith('export'):
                    line = re.sub('export ', '', line.strip()).strip()
                    # escape % character
                    line = re.sub('%', '%%', line)
                    if vars:  # add a comma
                        vars += ','
                    vars += line

        # Replace in supervisord
        with open(output_file, "r") as sources:
            lines = sources.readlines()
        with open(output_file, "w") as sources:
            for line in lines:
                sources.write(re.sub(r'{{ ENV_VARS }}', vars, line))


if __name__ == "__main__":
    main(sys.argv[1:])
