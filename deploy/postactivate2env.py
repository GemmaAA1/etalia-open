#!/usr/bin/python
"""To grab env var from postactivate"""

import sys
import getopt
import re
import subprocess

def main(argv):
    inputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
        print 'test.py -i <inputfile> -o <outputfile>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'postactivate2env.py -i <inputfile>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg

    if inputfile:
        print 'Open'
        vars = ''
        with open(inputfile, 'r') as f:
            for line in f:
                if line.strip().startswith('export'):
                    line = re.sub('export ', '', line.strip()).strip()
                    # escape % character
                    line = re.sub('%', '%%', line)
                    print line
                    if vars:  # add a comma
                        vars += ','
                    vars += line
        # store in environment variable
        print vars
        subprocess.call(["export ENV_FOR_SUPERVISOR=", '"'+vars+'"'], shell=True)


if __name__ == "__main__":
    main(sys.argv[1:])
