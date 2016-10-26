#!/usr/bin/env python

import argparse
import codegen
import sys
import os
import directories
import jobs


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description='Converts XML files into executable')

    parser.add_argument('--r', help = 'Remove executable (for a packaged directory  the --d option must be provided. The entire directory will then be removed from the build tree.)  generated by the XML file from the build tree. If the XML file has never been used to created an executable, this will be a no-operation', action = 'store_true')
    parser.add_argument('--d', help = 'Provide a packaging directory.',nargs = '?')
    parser.add_argument('xml file', metavar='XML File', nargs = '*', help = 'Will create an entry in the build tree for each XML file, provided the XML file is valid.')
    args = parser.parse_args()

    filename = vars(args)['xml file']
    dirname = vars(args)['d']
    if vars(args)['r'] == False:
        if dirname == None:
            fn = filename[0]
            directories.add_executable(fn)    
        else:
            directories.add_executable(dirname, filename)
    else:    
        if dirname == None:
            fn = directories.check_and_strip_name(filename[0])
            directories.detach_executable(fn)
        else:
            directories.detach_executable(dirname)

